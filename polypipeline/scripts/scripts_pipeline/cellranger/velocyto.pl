#!/usr/bin/env perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/";
use lib "$Bin/../../../../GenBo/lib/GenBoDB";
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages";
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moose;
use MooseX::Method::Signatures;
use file_util;
use GBuffer;
use GenBoProject;
use Term::Menus;
use Carp;
use autodie qw(system open);

  
my $bin_cecile=qq{$Bin/scripts/scripts_db_polypipeline};
my $bin_script_pipeline = qq{$Bin/scripts/scripts_pipeline};


my $projectName;
my $patients_name;
my $multi;
my $no_exec;
my $transcriptome;
my $choose_transcriptome;
my $cpu = 20;
$cpu = 32 if (`hostname` =~ /^master/);
my $force;
my $version;
my $help;

GetOptions(
	'project=s'										=> \$projectName,
	'patients=s'									=> \$patients_name,
	'multi'											=> \$multi,
	'cpu=i'											=> \$cpu,
	'transcriptome|reference=s'						=> \$transcriptome,
	'choose_transcriptome|choose_reference'			=> \$choose_transcriptome,
	'no_exec'										=> \$no_exec,
	'force'											=> \$force,
	'version=s'										=> \$version,
	'help'											=> \$help,
) || die("Error in command line arguments\n");

usage() if $help;
conferss("-project option is mandatory") unless ($projectName);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName , -version => $version);
my $patients = $project->get_only_list_patients($patients_name);
die("No patient in project $projectName") unless ($patients);
@$patients = sort {$a->name cmp $b->name} @$patients;
$multi = 1 if (grep {$_->getSampleProfile =~ /flex$/} @$patients);

my $dir_proj = $project->getProjectRootPath;
my $dir_cellranger = $project->getCountingDir('cellranger');
my $tmp_cellranger = $project->getAlignmentPipelineDir("cellranger");
$tmp_cellranger = $project->getAlignmentPipelineDir("cellranger_multi") if ($multi);
my $dir_velocyto = $project->getCountingDir('velocyto');

my @groups = map {$_->somatic_group} @$patients;
confess("No velocyto to run: no GEX library in project $projectName") unless (grep {/exp|nuclei/i} @groups);

# Vérifie que le pipeline n'ait pas déjà tourné
my @patients = @$patients;
foreach my $pat (@patients) {
	my $pname = $pat->name;
	if (-e "$dir_velocyto$pname.loom" and not $force) {
		warn "NEXT: '$pname.loom' already exists";
		@$patients = grep{$_->name ne $pname} @$patients;
		next;
	}
}
undef @patients;

# $exec version
qx{exec_singularity.sh velocyto.sif velocyto --version} =~ /^velocyto, version (\d+\.\d+\.\d+)$/ || confess("Error determining the version");
my $velocyto_version_nb = $1;
warn 'Velocyto version: '.$velocyto_version_nb;

# choose transcriptoome
my $transcriptome_dir = '/data-bipd/data-pure/public-data/10X/'.$project->genome_version.'/';
if ($choose_transcriptome or ($transcriptome and not -f $transcriptome.'/genes/genes.gtf')) {
	opendir(my $dh, $transcriptome_dir) || die "Can't opendir '$transcriptome_dir': $!";
	my @transcriptomes = sort { -M $transcriptome_dir.$a <=> -M $transcriptome_dir.$b } grep {-d $transcriptome_dir.$_ && /^refdata-(cellranger-|gex-)GRC/} readdir($dh);
	closedir ($dh);
	if (grep {/exp|adt/i} @groups) {
		confess("No transcriptome reference found in '$transcriptome_dir'") unless (scalar @transcriptomes);
		$transcriptome = $transcriptome_dir.$transcriptomes[0] if (scalar @transcriptomes == 1);
		if (scalar @transcriptomes > 1) {
			my $selected = prompt("Choose a transcritome reference for project ".$projectName.':', -menu=>\@transcriptomes);
			die unless ($selected and -d $transcriptome_dir.$selected);
			$transcriptome = $transcriptome_dir.$selected;
		}
	}
}
#my $transcriptome_gex ;
#opendir(my $dh, $transcriptome_dir) || die "Can't opendir '$transcriptome_dir': $!";
#my @transcriptomes_gex = sort { -M $transcriptome_dir.$a <=> -M $transcriptome_dir.$b } grep {-d $transcriptome_dir.$_ && /refdata-(cellranger-|gex-)GRC/} readdir($dh);
#closedir ($dh);
#confess("No transcriptome reference found in '$transcriptome_dir'") unless (scalar @transcriptomes_gex);
#$transcriptome_gex = $transcriptome_dir.$transcriptomes_gex[0] if (scalar @transcriptomes_gex == 1);
#$transcriptome_gex = $transcriptome_dir.prompt("Choose a transcritome reference for project ".$projectName.':', -menu=>\@transcriptomes_gex);


my $file_jobs = $dir_velocyto.'jobs_velocyto.txt';
open(my $jobs, ">$file_jobs") || die ("Can't open file '$file_jobs': $!");
foreach my $pat (@$patients) {
	my $pname = $pat->name;
	next if $pname =~ /_FCB$/;
	$pname =~ s/_FCA$//;
	unless ($pat->somatic_group =~ /exp|nuclei/i) {
		warn "NEXT: can not run velocyto for $pname: not a GEX library: ".$pat->somatic_group;
		next;
	}
	my $index = $project->getGenomeIndex($pat->alignmentMethod);
	$transcriptome = (qx/realpath $transcriptome/) if (-l $transcriptome);
	chomp $transcriptome;
	my $gtf = $transcriptome.'/genes/genes.gtf';
	confess("'$gtf' does not exist") unless (-f $gtf);
	
	system("add_calling_method.sh -project=$projectName -patient=$pname -method=velocyto") unless (grep (/velocyto/, @{$pat->getCallingMethods}) or $buffer->biocluster);
	
	my $dir_pat;
	my $cmd;
	my $cmd_velocyto;
	
	
	if ($multi) {
		# $tmp_cellranger/HD174_fixed/outs/per_sample_outs/HD174_fixed/sample_alignments.bam
		my $pool = $pat->family;
		# copy files in tmp unless already exist
		unless (-d "$tmp_cellranger$pool/outs/per_sample_outs/$pname") {
			$dir_pat = "$dir_cellranger/$pool/$pname";
			$dir_pat = "$dir_cellranger/$pool/outs/per_sample_outs/$pname" unless (-d $dir_pat);
			die ("No directory found for patient '$pname' in project $projectName :'$dir_pat'") unless ($dir_pat);
			$cmd = "cp -r $dir_pat $tmp_cellranger && " if (-d "$dir_pat/$pool/outs/per_sample_outs");
			$cmd = "mkdir -p $tmp_cellranger/$pool/outs/per_sample_outs/ && cp -r $dir_pat $tmp_cellranger/$pool/outs/per_sample_outs/ && " unless (-d "$dir_pat/$pool/outs/per_sample_outs");
		}
		# check gtf used by cellranger
		my $config_csv = "$tmp_cellranger$pool/outs/config.csv";
		$config_csv = "$dir_cellranger/$pool/outs/config.csv" unless (-f $config_csv);
		if (-f $config_csv) {
			open (my $config, '<', $config_csv);
			while (my $line = readline($config)) {
				chomp $line;
				$line =~ /^reference,(.*),?/ and last;
				$gtf = $1.'/genes/genes.gtf' if ($1);
				confess ("'$gtf' doesn't exist") unless (-f $gtf);
			}
			confess("GTF used by cellranger is different than the one selected for '$pname' : $gtf vs $transcriptome/genes/genes.gtf") if ($gtf ne $transcriptome.'/genes/genes.gtf');
		}
		
		$cmd_velocyto = $buffer->software('velocyto')." run -v -@ $cpu $tmp_cellranger$pool/outs/per_sample_outs/$pname/sample_alignments.bam $gtf"; # --outputfolder
		$cmd .= $cmd_velocyto;
		$cmd .= " && cp $tmp_cellranger$pool/outs/per_sample_outs/$pname/velocyto/sample_alignments_*.loom $dir_velocyto$pname.loom";
	}
	else {
		# copy files in tmp unless already exist
		unless (-d "$tmp_cellranger$pname/outs") {
			$dir_pat = "$dir_proj/$pname" if (-d "$dir_proj/$pname");
			$dir_pat = "$dir_cellranger/$pname" if (-d "$dir_cellranger/$pname");
			die ("No directory found for patient '$pname' in project $projectName") unless ($dir_pat);
			$cmd = "cp -r $dir_pat $tmp_cellranger && " if (-d "$dir_pat/outs/");
			$cmd = "mkdir -p $tmp_cellranger/$pname/outs && cp -r $dir_pat/* $tmp_cellranger/$pname/outs/ && " unless (-d "$dir_pat/outs/");
		}
		# check gtf used by cellranger
		my $cmd_file = "$tmp_cellranger$pname/_cmdline";
		$cmd_file = "$dir_cellranger/$pname/_cmdline" unless (-f $cmd_file);
		if (-f $cmd_file) {
			open (my $cmd_cellranger, '<', $cmd_file);
			my $line = readline($cmd_cellranger);
			chomp $line;
			$line =~ /--transcriptome[ =]([^\s]*) ?/;
			$gtf = $1.'/genes/genes.gtf' if ($1);
			confess ("'$gtf' doesn't exist") unless (-f $gtf);
			confess("GTF used by cellranger is different than the one selected for '$pname' : $gtf vs $transcriptome/genes/genes.gtf") if ($gtf ne $transcriptome.'/genes/genes.gtf');
			close($cmd_cellranger);
		}
		
		$cmd_velocyto = $buffer->software('velocyto')." run10x -v -@ $cpu $tmp_cellranger$pname $gtf";
		$cmd .= $cmd_velocyto;
		$cmd .= " && cp $tmp_cellranger/$pname/velocyto/$pname.loom $dir_velocyto";
	}
	print {$jobs} $cmd."\n";
	warn $cmd;
	$pat->update_software_version('velocyto', $cmd_velocyto, $velocyto_version_nb) unless ($no_exec);
	
}
	# Execution time is ~3h for a typical sample but might vary significantly by sequencing depth and cpu power.
	# ~6h avec ce script cpu 20 sur cluster bioinfo

#warn ("cat $file_jobs | $Bin/../../../../polyscripts/system_utility/run_cluster.pl -name velocyto -cpu=$cpu");
warn ("cat $file_jobs | run_cluster.pl -cpu=$cpu");
sleep(5) unless ($no_exec);
#my $exit = system ("cat $file_jobs | $Bin/../../../../polyscripts/system_utility/run_cluster.pl -name velocyto -cpu=$cpu") unless ($no_exec);
my $exit = system ("cat $file_jobs | run_cluster.pl -cpu=$cpu") unless ($no_exec);
exit($exit) if ($exit);

unless ($no_exec) {
	print "\t------------------------------------------\n";
	print("\tVelocyto:\n");
	print("\t$dir_velocyto\n");
	print "\t------------------------------------------\n\n";
}






sub usage {
	print "
$0
-------------
Obligatoires:
	project <s>                nom du projet
Optionels:
	patients <s>               noms de patients/échantillons, séparés par des virgules
	multi                      si cellranger multi a été utilisé pour les comptages
	cpu <i>                    nombre de cpu à utiliser, défaut: 20
	no_exec                    ne pas exécuter les commandes
	force                      relance le pipeline même s'il a déjà tourné
	help                       affiche ce message

";
	exit(1);
}


