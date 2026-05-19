#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/";
use lib "$Bin/../../../../GenBo/lib/GenBoDB";
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages";
use Logfile::Rotate;
use Cwd;
use PBS::Client;
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moose;
use MooseX::Method::Signatures;
use file_util;
use Class::Inspector;
use Digest::MD5::File ;
use GBuffer;
use GenBoProject;
use colored; 
use Config::Std;
use Text::Table;
use Time::Local 'timelocal';
use File::Temp qw/ tempfile tempdir /;
use Term::Menus;
use Proc::Simple;
use Storable;
use JSON::XS;
use XML::Simple qw(:strict);
use Carp;

  
my $bin_cecile=qq{$Bin/scripts/scripts_db_polypipeline};
my $bin_script_pipeline = qq{$Bin/scripts/scripts_pipeline};


my $projectName;
my $patients_name;
my $no_exec;
my $transcriptome;
my $choose_transcriptome;
my $cpu = 20;
my $force;
my $version;
my $help;

GetOptions(
	'project=s'										=> \$projectName,
	'patients=s'									=> \$patients_name,
	'cpu=i'											=> \$cpu,
	'transcriptome|reference'						=> \$transcriptome,
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

my $dir_proj = $project->getProjectRootPath;
my $dir_cellranger = $project->getCountingDir('cellranger');
my $tmp_cellranger = $project->getAlignmentPipelineDir("cellranger");
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


my $file_jobs = $dir_velocyto.'jobs_velocyto.txt';
open(my $jobs, ">$file_jobs") || die ("Can't open file '$file_jobs': $!");
foreach my $pat (@$patients) {
	my $pname = $pat->name;
	unless ($pat->somatic_group =~ /exp|nuclei/i) {
		warn "NEXT: can not run velocyto for $pname: not a GEX library: ".$pat->somatic_group;
		next;
	}
	my $transcriptome_dir = '/data-isilon/public-data/10X/'.$project->genome_version.'/';
	if ($choose_transcriptome or ($transcriptome and not -f $transcriptome.'/genes/genes.gtf')) {
		opendir(my $dh, $transcriptome_dir) || die "Can't opendir '$transcriptome_dir': $!";
		my @transcriptomes = sort { -M $transcriptome_dir.$a <=> -M $transcriptome_dir.$b } grep {-d $transcriptome_dir.$_ && /^refdata-(cellranger-|gex-)GRC/} readdir($dh);
		closedir ($dh);
		if (grep {/exp|adt/i} @groups) {
			confess("No transcriptome reference found in '$transcriptome_dir'") unless (scalar @transcriptomes);
			warn Dumper \@transcriptomes;
			$transcriptome = $transcriptome_dir.$transcriptomes[0] if (scalar @transcriptomes == 1);
			if (scalar @transcriptomes > 1) {
				my $selected = prompt("Choose a transcritome reference for project ".$projectName.':', -menu=>\@transcriptomes);
				die unless ($selected and -d $transcriptome_dir.$selected);
				$transcriptome = $transcriptome_dir.$selected;
			}
		}
	}
#	my $transcriptome_gex ;
#	opendir(my $dh, $transcriptome_dir) || die "Can't opendir '$transcriptome_dir': $!";
#	my @transcriptomes_gex = sort { -M $transcriptome_dir.$a <=> -M $transcriptome_dir.$b } grep {-d $transcriptome_dir.$_ && /refdata-(cellranger-|gex-)GRC/} readdir($dh);
#	closedir ($dh);
#	confess("No transcriptome reference found in '$transcriptome_dir'") unless (scalar @transcriptomes_gex);
#	$transcriptome_gex = $transcriptome_dir.$transcriptomes_gex[0] if (scalar @transcriptomes_gex == 1);
#	$transcriptome_gex = $transcriptome_dir.prompt("Choose a transcritome reference for project ".$projectName.':', -menu=>\@transcriptomes_gex);
	my $index = $project->getGenomeIndex($pat->alignmentMethod);
	$transcriptome = readlink $index unless ($transcriptome);
	my $gtf = $transcriptome.'/genes/genes.gtf';
	confess("'$gtf' does not exist") unless (-f $gtf);
	
	system("add_calling_method.sh -project=$projectName -patient=$pname -method=velocyto") unless (grep (/velocyto/, $pat->getCallingMethods));
	
	my $dir_pat;
	my $cmd;
	unless (-d "$tmp_cellranger$pname/outs/") {
		$dir_pat = "$dir_proj/$pname" if (-d "$dir_proj/$pname");
		$dir_pat = "$dir_cellranger/$pname" if (-d "$dir_cellranger/$pname");
		die ("No directory found for patient '$pname' in project ".$project->name) unless ($dir_pat);
		$cmd = "cp -r $dir_pat $tmp_cellranger && " if (-d "$dir_pat/outs/");
		$cmd = "mkdir $tmp_cellranger/$pname/outs && cp -r $dir_pat/* $tmp_cellranger/$pname/outs/ && " unless (-d "$dir_pat/outs/");
	}
	$cmd .= "singularity run --cleanenv -B $tmp_cellranger$pname -B $transcriptome_dir /data-beegfs/software/sif/velocyto.sif velocyto run10x -v -@ $cpu $tmp_cellranger$pname $gtf";
	$cmd .= " && cp $tmp_cellranger/$pname/velocyto/$pname.loom $dir_velocyto";
	print {$jobs} $cmd."\n";
	warn $cmd if ($no_exec);
	# Execution time is ~3h for a typical sample but might vary significantly by sequencing depth and cpu power.
	# ~6h avec ce script cpu 20
}

warn ("cat $file_jobs | run_cluster.pl -cpu=$cpu");
system ("cat $file_jobs | run_cluster.pl -cpu=$cpu") unless ($no_exec);







sub usage {
	print "
$0
-------------
Obligatoires:
	project <s>                nom du projet
Optionels:
	patients <s>               noms de patients/échantillons, séparés par des virgules
	cpu <i>                    nombre de cpu à utiliser, défaut: 20
	no_exec                    ne pas exécuter les commandes
	force                      relance le pipeline même s'il a déjà tourné
	help                       affiche ce message

";
	exit(1);
}


