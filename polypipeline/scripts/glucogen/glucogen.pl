#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/packages/";
use lib "$Bin/../../../GenBo/lib/GenBoDB/writeDB/";
use lib "$Bin/packages";
use Data::Dumper;
use Getopt::Long;
use Carp;
use IO::Prompt;
use Term::Menus;
use colored;
use Cwd 'abs_path';
use Net::SFTP;
use Archive::Tar;
use File::Util;
use Parallel::ForkManager;

use GBuffer;
my $buffer = new GBuffer;

my $project_name;
my $patient_names;
my $site;
my $set;
my $diabetome_project_name;
my $force;
my $checksum = 1;
my $no_exec;
my $fork = 1;
my $help;

my @sites = qw/glucogenpitie glucogenlyon glucogentoul/;

GetOptions(
	'project|genome_project=s'	=> \$project_name,
	'diabetome_project=s'		=> \$diabetome_project_name,
	'patients=s'				=> \$patient_names,
	'site|glucogen|labo=s'		=> \$site,
	'set=i'						=> \$set,
	'force'						=> \$force,
	'checksum!'					=> \$checksum,
	'no_exec'					=> \$no_exec,
	'fork=i'					=> \$fork,
	'help'						=> \$help,
) or (warn("\nError in command line arguments\n") && usage());
usage() if ($help);
die('Enter a project name') unless ($project_name);

my $project = $buffer->newProject( -name => $project_name );
warn $project_name;
my $project_desc = $project->description;
confess ('Project $project_name is not glucogen project. Check capture and description.') unless ($project_desc =~ /glucogen/i and $project->isGenome);
my $patients = $project->get_only_list_patients($patient_names);
die("No patient in project ".$project_name."\n") unless ($patients);
@$patients = sort {$a->name cmp $b->name} @$patients;
# site
$project_desc =~ /(pitie|lyon|toul)/i;
$site = lc($1) unless ($site);
die('Enter a site: pitie, lyon or toul') unless ($site);

# Keep only patients corresponding to the diabetome project
if ($diabetome_project_name) {
	my $diabetome_project = $buffer->newProject( -name => $diabetome_project_name );
	$project_desc = $diabetome_project->description;
	confess ('Diabetome project $diabetome_project_name is not diabetome. Check capture and description.') 
		unless ($project_desc =~ /glucogen/i and $project_desc =~ /diabetome/i and grep {$_->name eq 'DIABETomeV1_hg38'} @{$diabetome_project->getCaptures});
	$project_desc =~ /(pitie|lyon|toul)/i;
	$site = lc($1) unless ($site);
	confess("Sites not matching between project $project_name ($site) and diabetome project $diabetome_project_name ($1)") unless ($1 and $site eq lc($1));
	my $patients_diabetome = $diabetome_project->get_only_list_patients($patient_names);
	die("No patient in diabetome project ".$project_name."\n") unless ($patients);
	my $patients_both;
	foreach my $pat (@$patients) {
		my @pat_grep = grep {$pat->name eq $_->name} @$patients_diabetome;
		push (@$patients_both, $pat_grep[0]) if (scalar @pat_grep);
	}
	$patients = $patients_both;
	warn "Keeping only patients present in both project $project_name and diabetome project $diabetome_project_name: ". scalar @$patients;
}

# site
$site = 'toul' if ($site =~ /toulouse/i);
my @site = grep{/$site/i} @sites;
$site = uc(@site[0]);
die("Choose a site in '".join(', ',@sites)."'.") unless (scalar @site == 1);
warn $site;
# set
$project_desc =~ /set-?(\d+)/i;
$set = $1 unless ($set);
die('Enter a set number > 0') unless ($set > 0);
warn 'Set '.$set;

# Vérifie qu'il y a bien les 3 méthodes de calling SV
my @SVmethods = map{values @{$_->callingSVMethods}} @$patients;
if (grep {/canvas|manta|wisecondor/} @SVmethods == 3*scalar @SVmethods and $project->isGenome) {
	warn("Adding SV methods for genomes");
	system("add_calling_method.sh -project=$project_name -methods=canvas,manta,wisecondor");
}

# Setting validation_db to 'glucogen' (for hotspot)
unless ($project->validation_db eq 'glucogen') {
	my $pid = $project->id;
	my $dbh = $buffer->dbh();
	my $sql = "UPDATE `PolyprojectNGS`.`projects` SET `validation_db`='glucogen' WHERE `project_id`=?;";
	my $sth = $dbh->do($sql, undef, $pid);
	confess ("ERROR: Can't update validation_db to 'glucogen' for projet $project_name (id: ".$pid."):\n"
		. 'Statement: '.$sql =~ s/`project_id`=\?/`project_id`=$pid/r ."\nDB Error: ".$dbh->errstr) if ($dbh->errstr or not $sth);
}


my $dir_download = "/data-pure/workspace/download/glucogen/" . uc($site) . '/set' . $set . '/';
warn $dir_download;
confess ("No directory '$dir_download'") unless (-d $dir_download);


# Vérifie qu'il n'y a pas eu d'erreur lors du téléchargement
my $f = File::Util->new();
my @report = $f->list_dir($dir_download => { files_only => 1, files_match => qr/cnrgh_dl_report_\d{8}_\d{6}.json$/});
my $report_download = $dir_download.$report[-1];
if (scalar @report) {
	open (my $fh, '<', $report_download) or confess("Can't open '$report_download': $!");
	my $error_download = 0;
	while (<$fh>) {
	    $error_download ++ if /ERROR/;
	}
	close $fh;
	die("$error_download errors found in the lastest download report '$report_download'. Please rerun the download.") if ($error_download);
}


my $pm   = new Parallel::ForkManager($fork);
my $erreur_fork = 0;
my @patients_in_error;
$pm->run_on_finish(sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $block) = @_;
    if ($exit_signal || $exit_code != 0) {
        $erreur_fork = 1;
		my ($name, $barcode) = split(/\t/, $ident);
		push @patients_in_error, { name => $name, barcode => $barcode };
		warn "ERROR: Patient '$name' ($barcode) failed (PID: $pid, Exit: $exit_code)\n";
       	} else {
		warn "OK: Patient '$ident' completed successfully\n";
	}
});

$buffer->disconnect;
foreach my $pat (@$patients) {
	my $pid = $pm->start($pat->name."\t".$pat->barcode) and next; #unless($fork == 1);
	my $bc = $pat->barcode;
	warn $pat->name.' -> '.$bc;
	my $f = File::Util->new();
	
	# CRAM
	my @downloaded_cram_files = $f->list_dir($dir_download => { files_only => 1, files_match => qr/$bc.*\.cram(\.crai)?$/});
	my $cram = $pat->getCramFileName;
	my $idxstats = $cram =~ s/cram$/idxstats/r;
	if (scalar @downloaded_cram_files == 2) {
		my ($downloaded_cram, $downloaded_crai) = @downloaded_cram_files;
#		check_md5sum($downloaded_cram.'.md5', $dir_download) unless ($no_exec || not $checksum);
		unless (-e $cram and not $force) {
			my $cmd_mv = "mv $dir_download$downloaded_cram $cram && mv $dir_download$downloaded_crai $cram.crai";
	#		$cmd_mv =~ s/mv /cp /g;
			warn $cmd_mv."\n";
			my $exit = system($cmd_mv) unless ($no_exec);
			confess ("Error $cmd_mv") if ($exit);
		}
	}
	elsif (scalar @downloaded_cram_files != 2 and not -s $cram) {
		confess ("2 files expected (cram, crai), got ".scalar @downloaded_cram_files.': '.Dumper \@downloaded_cram_files);
	}
	unless (-s $cram.'.crai') {
		my $cmd_idxstats = $buffer->software("samtools")." index $cram";
		warn $cmd_idxstats."\n";
		system($cmd_idxstats) unless ($no_exec);
		confess("Error '$cram.crai'") unless (-s $cram.'.crai');
	}
	unless (-s $idxstats and not $force) {
		my $cmd_idxstats = $buffer->software("samtools")." idxstats $cram > $idxstats";
		warn $cmd_idxstats."\n";
		system($cmd_idxstats) unless ($no_exec);
		confess("Error '$idxstats'") unless (-s $idxstats);
	}
	
	
	# Variants/VCF
	# Extraction de l'archive
	my @archive = $f->list_dir($dir_download => { files_only => 1, files_match => qr/deliverable_VARIANT_.*$bc.*\.tar\.gz$/});
	confess ("1 file expected (.tar.gz), got ".scalar @archive.': '.Dumper \@archive) unless (scalar @archive == 1);
	my $dir_downloaded_vcf_pat = $dir_download.$archive[0] =~ s/\.tar\.gz$/\//r;
	unless (-d $dir_downloaded_vcf_pat && not $force) {
		my $cmd_untar = "tar -C $dir_download -xf ".$dir_download.$archive[0];
		warn $cmd_untar;
		warn("Extracting archive...\n");
		system ($cmd_untar);
		confess("Can't extract archive '".$dir_download.$archive[0]."'") unless (-d $dir_downloaded_vcf_pat);
		my @md5 = $f->list_dir($dir_downloaded_vcf_pat => { files_only => 1, files_match => qr/$bc.*\.md5$/});
		check_md5sum(@md5, $dir_downloaded_vcf_pat) unless ($no_exec || not $checksum);
	}
	
	my $cmd_mv;
	# Haplotypecaller4
	my @gatk4 = $f->list_dir($dir_downloaded_vcf_pat => { files_only => 1, files_match => qr/$bc.*_gatk4(_annotated|\.g)\.vcf\.gz(\.tbi|\.md5)?$/});
	my $vcf_gatk4 = $pat->getVariationsFileName('haplotypecaller4');
	unless (-e $vcf_gatk4 and not $force) {
		my @vcf = grep(/_gatk4_annotated\.vcf\.gz$/, @gatk4);
		my @tbi = grep(/_gatk4_annotated\.vcf\.gz\.tbi$/, @gatk4);
		confess("Expected 1 gatk4 vcf and vcf.tbi file in $dir_downloaded_vcf_pat") unless (scalar @vcf ==1  && scalar @tbi ==1 );
		$cmd_mv = "mv ".$dir_downloaded_vcf_pat.$vcf[0]." $vcf_gatk4 && mv ".$dir_downloaded_vcf_pat.$tbi[0]." $vcf_gatk4.tbi\n";
	}
	my $gvcf_gatk4 = $pat->gvcfFileName('haplotypecaller4');
	unless (-e $gvcf_gatk4 and not $force) {
		my @gvcf = grep(/_gatk4\.g\.vcf\.gz$/, @gatk4);
		my @tbi = grep(/_gatk4\.g\.vcf\.gz\.tbi$/, @gatk4);
		confess("Expected 1 gatk4 g.vcf and g.vcf.tbi file in $dir_downloaded_vcf_pat") unless (scalar @gvcf ==1  && scalar @tbi ==1 );
		$cmd_mv .= "mv ".$dir_downloaded_vcf_pat.$gvcf[0]." $gvcf_gatk4 && mv ".$dir_downloaded_vcf_pat.$tbi[0]." $gvcf_gatk4.tbi\n";
	}
	
	
	# Jax-CNV
	my @jaxcnv = $f->list_dir($dir_downloaded_vcf_pat => { files_only => 1, files_match => qr/$bc.*_jaxcnv\.bed\.gz(\.tbi|\.md5)?$/});
	my $vcf_jaxcnv = $pat->getVariationsFileName('jax-cnv') =~ s/\.vcf\.gz$/\.bed\.gz/r;
	unless (-e $vcf_jaxcnv and not $force) {
		my @vcf = grep(/_jaxcnv\.bed\.gz$/, @jaxcnv);
		my @tbi = grep(/_jaxcnv\.bed\.gz\.tbi$/, @jaxcnv);
		confess("Expected 1 jaxcnv vcf and vcf.tbi file in $dir_downloaded_vcf_pat") unless (scalar @vcf ==1  && scalar @tbi ==1 );
		$cmd_mv .= "mv ".$dir_downloaded_vcf_pat.$vcf[0]." $vcf_jaxcnv && mv ".$dir_downloaded_vcf_pat.$tbi[0]." $vcf_jaxcnv.tbi\n";
	}
	
	# Manta
	my @manta = $f->list_dir($dir_downloaded_vcf_pat => { files_only => 1, files_match => qr/$bc.*_manta_diploidSV\.vcf\.gz(\.tbi|\.md5)?$/});
	my $vcf_manta = $pat->getVariationsFileName('manta');
	unless (-e $vcf_manta and not $force) {
		my @vcf = grep(/_manta_diploidSV\.vcf\.gz$/, @manta);
		my @tbi = grep(/_manta_diploidSV\.vcf\.gz\.tbi$/, @manta);
		confess("Expected 1 manta vcf and vcf.tbi file in $dir_downloaded_vcf_pat") unless (scalar @vcf ==1  && scalar @tbi ==1 );
		$cmd_mv .= "mv ".$dir_downloaded_vcf_pat.$vcf[0]." $vcf_manta && mv ".$dir_downloaded_vcf_pat.$tbi[0]." $vcf_manta.tbi\n";
	}
	
	
	# Octopus
	my @octopus = $f->list_dir($dir_downloaded_vcf_pat => { files_only => 1, files_match => qr/$bc.*_octopus_annotated\.vcf\.gz(\.tbi|\.md5)?$/});
	my $vcf_octopus = $pat->getVariationsFileName('octopus');
	unless (-e $vcf_octopus and not $force) {
		my @vcf = grep(/_octopus_annotated\.vcf\.gz$/, @octopus);
		my @tbi = grep(/_octopus_annotated\.vcf\.gz\.tbi$/, @octopus);
		confess("Expected 1 octopus vcf and vcf.tbi file in $dir_downloaded_vcf_pat") unless (scalar @vcf ==1  && scalar @tbi ==1 );
		$cmd_mv .= "mv ".$dir_downloaded_vcf_pat.$vcf[0]." $vcf_octopus && mv ".$dir_downloaded_vcf_pat.$tbi[0]." $vcf_octopus.tbi\n";
	}
	
#	$cmd_mv =~ s/mv /cp /g;
	warn $cmd_mv;
	system($cmd_mv) unless ($no_exec);
	$pm->finish;
}
$pm->wait_all_children();
if ($erreur_fork) {
    my $err_mess = "\n========================================\n";
    $err_mess .= "ERROR: Processing failed for " . scalar(@patients_in_error) . " patient(s):\n";
    $err_mess .= "  - " . join("\n  - ", map{$_->{name} . ' (' . $_->{barcode} . ')'} @patients_in_error ) . "\n";
    $err_mess .= "========================================\n";
    confess($err_mess);
}


unless ($no_exec) {
	print("----------DONE----------\n");
	my $cmd_pipeline = "$Bin/../../bds_pipeline.pl -project $project_name -steps coverage,binary_depth,canvas,wisecondor,calling_wisecondor ";
	$cmd_pipeline .= "-patients=$patient_names " if ($patient_names);
	$cmd_pipeline .= "-force 1 " if ($force);
	my $cmd_cache = "$Bin/../../bds_cache.pl -project $project_name";
	$cmd_cache .= "-force 1 " if ($force);
	my $cmd_diabetome;
	if ($diabetome_project_name) {
		print("Now, run glucogen_diabetome.pl:\n");
		print($cmd_diabetome."\n");
		$cmd_diabetome = "$Bin/glucogen_diabetome.pl -diabetome_project $diabetome_project_name -genome_project $project_name -fork $fork ";
		$cmd_diabetome .= "-patients $patient_names " if ($patient_names);
		$cmd_diabetome .= "-force 1 " if ($force);
		$cmd_pipeline = "$Bin/../../bds_pipeline.pl -project $project_name -steps coverage,binary_depth ";
		$cmd_pipeline .= "-force 1 " if ($force);
		$cmd_cache =~ s/$project_name/$diabetome_project_name/;
	}
	print("Then, run coverage (and cnv for genome), then cache:\n");
	print($cmd_pipeline."\n");
	print($cmd_cache."\n");
}
print "\n";



sub check_md5sum {
	my @md5 = @_;
	my $dir = pop @md5;
	confess("Directory missing: '$dir'") unless (-d $dir);
	print("Checking md5sum...\n");
	foreach my $md5 (@md5) {
		confess("Error md5 file missing: '$dir$md5'") unless (-e $dir.$md5);
		my $exit = system("cd $dir && md5sum -c $md5");
		confess("Error while checking md5sum : $md5") if ($exit);
	}
}

sub usage {
	my $no_exit = shift @_;
	print "
$0
Move glucogen (cram and vcf) files downloaded from cnrgh in /data-pure/workspace/download/glucogen/GLUCOGEN<site>/<set>/ to the <project> directory
-----------------	
Mandatory arguments
	-project|genome_project <s>  project name
	
Optional arguments
	-diabetome_project <s>       corresponding diabetome project name
	-site <s>                    glucogen project, in ".join(',',@sites)."
	-set <i>                     set number
	-patients <s>                patient names separated with a comma
	-checksum/nochecksum         enables or disables the md5sum check [enabled]
	-force                       overwrite files if existing
	-fork <i>                    number of forks to use in parallele
	-no_exec                     do not execute the commands
	-help                        display this help message and exit

";
	exit(1) unless ($no_exit);
}
