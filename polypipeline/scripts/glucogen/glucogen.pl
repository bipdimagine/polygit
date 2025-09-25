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
my $force;
my $checksum = 1;
my $no_exec;
my $no_die;
my $fork = 1;
my $help;

my @sites = qw/glucogenpitie glucogenlyon glucogentoul/;

GetOptions(
	'project=s'				=> \$project_name,
	'patients=s'			=> \$patient_names,
	'site|glucogen|labo=s'	=> \$site,
	'set=i'					=> \$set,
	'force'					=> \$force,
	'checksum!'				=> \$checksum,
	'no_exec'				=> \$no_exec,
	'no_die'				=> \$no_die,
	'fork=i'				=> \$fork,
	'help'					=> \$help,
) or die ("Error in command line arguments");
usage() if ($help);
die('Enter a project name') unless ($project_name);
die('Enter a set number > 0') unless ($set > 0);
die('Enter a site') unless ($set > 0);


my $project = $buffer->newProject( -name => $project_name );
warn $project_name;
my $patients = $project->get_only_list_patients($patient_names);
die("No patient in project ".$project_name."\n") unless ($patients);

my $dir_download = "/data-isilon/download/glucogen/";
$site = 'toul' if ($site =~ /toulouse/i);
my @site = grep{/$site/i} @sites;
$site = @site[0];
die("Choose a site in '".join(', ',@sites)."'.") unless (scalar @site == 1);
$dir_download .= uc($site).'/set'.$set.'/';
warn $dir_download;
die ("No directory '$dir_download'") unless (-d $dir_download);


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
	die("$error_download errors found in the lastestt download report '$report_download'. Please rerun the download.") if ($error_download);
}


my $pm   = new Parallel::ForkManager($fork);
foreach my $pat (@$patients) {
	my $pid = $pm->start and next unless($fork == 1);
	my $bc = $pat->barcode;
	warn $pat->name.' -> '.$bc;
	my $f = File::Util->new();
	
	# CRAM
	my @downloaded_cram_files = $f->list_dir($dir_download => { files_only => 1, files_match => qr/$bc.*\.cram(\.crai)?$/});
	if (scalar @downloaded_cram_files == 2) {
		my ($downloaded_cram, $downloaded_crai) = @downloaded_cram_files;
	#	check_md5sum($downloaded_md5, $dir_download) unless ($no_exec || not $checksum);
		my $cram = $pat->getCramFileName;
		unless (-e $cram and not $force) {
			my $cmd_mv = "mv $dir_download$downloaded_cram $cram && mv $dir_download$downloaded_crai $cram.crai";
	#		$cmd_mv =~ s/mv /cp /g;
			warn $cmd_mv."\n";
			system($cmd_mv) unless ($no_exec);
		}
		my $idxstats = $cram =~ s/cram$/idxstats/r;
		unless (-e $idxstats and not $force) {
			my $cmd_idxstats = $buffer->software("samtools")." idxstats $cram > $idxstats";
			warn $cmd_idxstats."\n";
			system($cmd_idxstats) unless ($no_exec);
		}
	}
	else {
		warn ("2 files expected (cram, crai), got ".scalar @downloaded_cram_files.': '.Dumper \@downloaded_cram_files) if ($no_die);
		die ("2 files expected (cram, crai), got ".scalar @downloaded_cram_files.': '.Dumper \@downloaded_cram_files) unless ($no_die);
	}
	
	
	# Variants/VCF
	# Extraction de l'archive
	my @archive = $f->list_dir($dir_download => { files_only => 1, files_match => qr/deliverable_VARIANT_.*$bc.*\.tar\.gz$/});
	die ("1 file expected (.tar.gz), got ".scalar @archive.': '.Dumper \@archive) unless (scalar @archive == 1 or $no_die);
	warn ("1 file expected (.tar.gz), got ".scalar @archive.': '.Dumper \@archive) and warn"\n" and next unless (scalar @archive == 1 or not $no_die);
	my $dir_downloaded_vcf_pat = $dir_download.$archive[0] =~ s/\.tar\.gz$/\//r;
	unless (-d $dir_downloaded_vcf_pat && not $force) {
		my $cmd_untar = "tar -C $dir_download -xf ".$dir_download.$archive[0];
		warn $cmd_untar;
		warn("Extracting archive...\n");
		system ($cmd_untar);
		die("Can't extract archive '".$dir_download.$archive[0]."'") unless (-d $dir_downloaded_vcf_pat);
#		die ("'$dir_downloaded_vcf_pat' does not exist. Check that the archive was extracted.") unless (-d $dir_downloaded_vcf_pat && not $no_exec);
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
		die("Expected 1 gatk4 vcf and vcf.tbi file in $dir_downloaded_vcf_pat") unless (scalar @vcf ==1  && scalar @tbi ==1 );
		$cmd_mv = "mv ".$dir_downloaded_vcf_pat.$vcf[0]." $vcf_gatk4 && mv ".$dir_downloaded_vcf_pat.$tbi[0]." $vcf_gatk4.tbi\n";
	}
	
	my $gvcf_gatk4 = $pat->gvcfFileName('haplotypecaller4');
	unless (-e $gvcf_gatk4 and not $force) {
		my @gvcf = grep(/_gatk4\.g\.vcf\.gz$/, @gatk4);
		my @tbi = grep(/_gatk4\.g\.vcf\.gz\.tbi$/, @gatk4);
		die("Expected 1 gatk4 g.vcf and g.vcf.tbi file in $dir_downloaded_vcf_pat") unless (scalar @gvcf ==1  && scalar @tbi ==1 );
		$cmd_mv .= "mv ".$dir_downloaded_vcf_pat.$gvcf[0]." $gvcf_gatk4 && mv ".$dir_downloaded_vcf_pat.$tbi[0]." $gvcf_gatk4.tbi\n";
	}
	
	
	# Jax-CNV
	my @jaxcnv = $f->list_dir($dir_downloaded_vcf_pat => { files_only => 1, files_match => qr/$bc.*_jaxcnv\.bed\.gz(\.tbi|\.md5)?$/});
	my $vcf_jaxcnv = $pat->getVariationsFileName('jax-cnv') =~ s/\.vcf\.gz$/\.bed\.gz/r;
	unless (-e $vcf_jaxcnv and not $force) {
		my @vcf = grep(/_jaxcnv\.bed\.gz$/, @jaxcnv);
		my @tbi = grep(/_jaxcnv\.bed\.gz\.tbi$/, @jaxcnv);
		die("Expected 1 jaxcnv vcf and vcf.tbi file in $dir_downloaded_vcf_pat") unless (scalar @vcf ==1  && scalar @tbi ==1 );
		$cmd_mv .= "mv ".$dir_downloaded_vcf_pat.$vcf[0]." $vcf_jaxcnv && mv ".$dir_downloaded_vcf_pat.$tbi[0]." $vcf_jaxcnv.tbi\n";
	}
	
	# Manta
	my @manta = $f->list_dir($dir_downloaded_vcf_pat => { files_only => 1, files_match => qr/$bc.*_manta_diploidSV\.vcf\.gz(\.tbi|\.md5)?$/});
	my $vcf_manta = $pat->getVariationsFileName('manta');
	unless (-e $vcf_manta and not $force) {
		my @vcf = grep(/_manta_diploidSV\.vcf\.gz$/, @manta);
		my @tbi = grep(/_manta_diploidSV\.vcf\.gz\.tbi$/, @manta);
		die("Expected 1 manta vcf and vcf.tbi file in $dir_downloaded_vcf_pat") unless (scalar @vcf ==1  && scalar @tbi ==1 );
		$cmd_mv .= "mv ".$dir_downloaded_vcf_pat.$vcf[0]." $vcf_manta && mv ".$dir_downloaded_vcf_pat.$tbi[0]." $vcf_manta.tbi\n";
	}
	
	
	# Octopus
	my @octopus = $f->list_dir($dir_downloaded_vcf_pat => { files_only => 1, files_match => qr/$bc.*_octopus_annotated\.vcf\.gz(\.tbi|\.md5)?$/});
	my $vcf_octopus = $pat->getVariationsFileName('octopus');
	unless (-e $vcf_octopus and not $force) {
		my @vcf = grep(/_octopus_annotated\.vcf\.gz$/, @octopus);
		my @tbi = grep(/_octopus_annotated\.vcf\.gz\.tbi$/, @octopus);
		die("Expected 1 octopus vcf and vcf.tbi file in $dir_downloaded_vcf_pat") unless (scalar @vcf ==1  && scalar @tbi ==1 );
		$cmd_mv .= "mv ".$dir_downloaded_vcf_pat.$vcf[0]." $vcf_octopus && mv ".$dir_downloaded_vcf_pat.$tbi[0]." $vcf_octopus.tbi\n";
	}
	
#	$cmd_mv =~ s/mv /cp /g;
	warn $cmd_mv."\n";
	system($cmd_mv) unless ($no_exec);
	$pm->finish;
}
$pm->wait_all_children();

unless ($no_exec) {
	my $cmd_pipeline = "bds_pipeline_rocks.sh -project=$project_name -steps=coverage,binary_depth";
	my $cmd_pipeline_cnv = "bds_pipeline_rocks.sh -project=$project_name -steps=canvas,wisecondor,calling_wisecondor ";
	$cmd_pipeline_cnv .= "-patients=$patient_names " if ($patient_names);
	print("--------DONE--------\n");
	print("For diabetome project, run\n");
	print("/home/mperin/git/polygit/polypipeline/scripts/glucogen/glucogen_diabetome.pl -project=<diabetome_project> -genome_project=$project_name\n");
	print("Now, run coverage and cnv on the project:\n");
	print($cmd_pipeline."\n");
	print($cmd_pipeline_cnv."\n");
	print("/data-isilon/bipd-src/HG38/polypipeline/bds_cache.pl -project=$project_name\n");
}
print "\n";

sub check_md5sum {
	my @md5 = @_;
	my $dir = pop @md5;
	confess("Directory missing: '$dir'") unless (-d $dir);
	print("Checking md5sum...\n");
	foreach my $md5 (@md5) {
		confess("md5 file missing: '$dir$md5'") unless (-e $dir.$md5);
		my $exit = system("cd $dir && md5sum -c $md5");
		die("Error while checking md5sum : $md5") if ($exit);
	}
}

sub usage {
	my $no_exit = shift @_;
	print "
$0
-----------------	
Mandatory arguments
	-project <s>              project name
	-set <i>                  set number
	
Optional arguments
	-site <s>                 glucogen project, in ".join(',',@sites)."
	-patients <s>             patient names separated with a comma
	-checksum/nochecksum      enables or disables the md5sum check [enabled]
	-force                    overwrite files if existing
	-fork <i>                 number of forks to use in parallele
	-no_exec                  do not execute the commands
	-no_die                   do not die if the files are not found
	-help                     display this help message and exit

";
	exit(1) unless ($no_exit);
}
