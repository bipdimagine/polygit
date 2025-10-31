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
#use bds_steps;   
use file_util;
use Class::Inspector;
use Digest::MD5::File ;
use GBuffer;
use GenBoProject;
use colored; 
use Config::Std;
use Text::Table;
use File::Temp qw/ tempfile tempdir /;
use Term::Menus;
use Proc::Simple;
use Storable;
use JSON::XS;
use XML::Simple qw(:strict);
use Cwd 'abs_path';

  

my $projectName;
my $patients_name;
my @steps;
my $lane;
my $mismatch = 0;
my $feature_ref;
my $cmo_ref;
my $no_exec;
my $aggr_name;
my $chemistry;
my $create_bam;
my $cpu = 20;
my $help;

GetOptions(
	'project=s'					=> \$projectName,
	'patients=s'				=> \$patients_name,
	'steps=s{1,}'				=> \@steps,
	'lane|nb_lane=i'			=> \$lane,
	'mismatches=i'				=> \$mismatch,
	'create_bam!'				=> \$create_bam,
	'feature_ref|feature_csv=s'	=> \$feature_ref,
#	'cmo_ref|cmo_csv=s'			=> \$cmo_ref,
	'aggr_name=s'				=> \$aggr_name,
	'chemistry=s'				=> \$chemistry,
	'cpu=i'						=> \$cpu,
	'no_exec'					=> \$no_exec,
	'help'						=> \$help,
) || die("Error in command line arguments\n");

usage() if $help;
die("-project argument is mandatory") unless ($projectName);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
my $patients = $project->get_only_list_patients($patients_name);
die("No patient in project $projectName") unless ($patients);

my $run = $project->getRun();
my $run_name = $run->plateform_run_name;
my $type = $run->infosRun->{method};
my $machine = $run->infosRun->{machine};

my $exec = "cellranger";
$exec .= '-atac' if ($type eq 'atac');
$exec .= '-arc' if ($type eq 'arc');
$exec = 'spaceranger' if  ($type eq 'spatial');
warn $exec;
	

#my $dir = $project->getProjectRootPath();
my $dir = $project->getCountingDir('cellranger');
$dir = $project->getCountingDir('spaceranger') if ($type eq 'spatial');
warn $dir;

@steps = split(/,/, join(',',@steps));
my $list_steps = ['demultiplex', 'teleport', 'count', 'aggr', 'aggr_vdj', 'tar', 'cp', 'cp_web_summary', 'infos', 'all'];
#my @correct_steps = grep{/@steps/} @$list_steps;
#undef @steps if (@steps and scalar @correct_steps == 0);
unless (@steps) {
	my %Menu_1 = (
		Item_1 => {
			Text   => "]Convey[",
			Convey => $list_steps,
		},
		Select => 'Many',
		Banner => "   Select steps:"
	);
	@steps = &Menu( \%Menu_1 );
	die if ( @steps eq ']quit[' );
}
warn 'steps='.join(',',@steps);

die("cpu must be in [1;40], given $cpu") unless ($cpu > 0 and $cpu <= 40);


my @patient_names = map {$_->name} @$patients;
my @invalid_names = grep { $_ !~ /^[A-Za-z0-9_-]+$/ } @patient_names;
die ("Patient names can only contain letters, numbers, hyphens and underscores. Sample names ".join(@invalid_names, ', ')." are invalids.") if (@invalid_names);

my $hfamily;
my $prog;
foreach my $patient (@{$patients}) {
	my $sbc = $patient->barcode; 
	my $pbc = $patient->barcode2;
	my $hpatient;
	$hpatient->{name}= $patient->name();
	$hpatient->{group}= $patient->somatic_group();
	$hpatient->{bc} = $pbc;
	$hpatient->{objet} = $patient;
	push(@{$hfamily->{$sbc}},$hpatient);
#	$prog = $patient->alignmentMethod();
}
warn Dumper $hfamily;
die;



#------------------------------
# DEMULIPLEXAGE
#------------------------------
if (grep(/demultiplex|^all$/i, @steps)){
	
	my $bcl_dir = $run->bcl_dir;
	warn $bcl_dir;
	my $config = XMLin("$bcl_dir/RunInfo.xml", KeyAttr => { reads => 'Reads' }, ForceArray => [ 'reads', 'read' ]);
	my $lane_config = $config->{Run}->{FlowcellLayout}->{LaneCount};
	my $prompt = prompt("$lane_config lanes found in the RunInfo.xml. You entered -lane=$lane. Continue ? (y/n)  ", -yes) if ($lane and $lane != $lane_config);
	die if ($prompt);
	unless ($lane) {
		$lane = $lane_config;
		warn 'LaneCount=',$lane;
	}
	
	my $tmp = $project->getAlignmentPipelineDir("cellranger_demultiplex");
	warn $tmp;

	my $sampleSheet = $bcl_dir."/sampleSheet.csv";
	open (SAMPLESHEET,">$sampleSheet");
	print SAMPLESHEET "Lane,Sample,Index\n";
	
	for (my $i=1; $i<=$lane; $i++){
		foreach my $k (sort keys(%$hfamily)){
			my @pat = @{$hfamily->{$k}};
			print SAMPLESHEET $i.",".$pat[0]->{group}.",".$k."\n";
		}
	}
	close(SAMPLESHEET);
	
	my $cmd = "cd $tmp; $Bin/../../demultiplex/demultiplex.pl -dir=$bcl_dir -run=$run_name -hiseq=10X -sample_sheet=$sampleSheet -cellranger_type=$exec -mismatch=$mismatch";
	warn $cmd;
	my $exit = system ($cmd) unless ($no_exec);
	exit($exit) if ($exit);
	unless ($@ or $no_exec) {
		print "\t------------------------------------------\n";
		print("Check the demultiplex stats\n");
		print("https://www.polyweb.fr/NGS/demultiplex/$run_name/$run_name\_laneBarcode.html\n\n");
		system ("firefox https://www.polyweb.fr/NGS/demultiplex/$run_name/$run_name\_laneBarcode.html &");
		print "\t------------------------------------------\n";
	}
}




#------------------------------
# COUNT
#------------------------------
if (grep(/count|^all$/i, @steps)){
	
	$create_bam = 'false' unless ($create_bam);
	$create_bam = 'true' if ($create_bam ne 'false');
	warn "create-bam=$create_bam";
	
	 if ($chemistry) {
 		warn "chemistry=$chemistry";
		my @chemistries = qw/auto threeprime SC3Pv1 SC3Pv2 SC3Pv3 SC3Pv3HT SC3Pv4 fiveprime SC5P-PE SC5P-R2 SC5P-R2-v3 SC5PHT SC-FB SFRP MFRP MFRP-R1/;
		die ("Chemistry option '$chemistry' not valid, should be one of: ". join(', ', @chemistries)) unless (grep { $_ eq $chemistry } @chemistries);
	 }

	my $index = $project->getGenomeIndex($prog);
	my $set = $index."/probe_sets/";
	$set .= "Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv";
	$set .= "Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv" if ($project->getVersion() =~ /^MM/);
	my $tmp = $project->getAlignmentPipelineDir("cellranger_count");
	warn $tmp;
	
	
	open (JOBS, ">$dir/jobs_multi.txt");
	foreach my $k (keys(%$hfamily)){
		my @pat = @{$hfamily->{$k}};
		my $poolName = $pat[0]->{group};
		system("mkdir $dir$poolName") unless (-d "$dir$poolName");
		system("mkdir $dir$poolName") unless (-d "$dir$poolName");
		my $pobj=$pat[0]->{objet};
		my $seq_dir = $pobj->getSequencesDirectory();
		my $tmp_fastq = $tmp.'fastq/';
		my $cmd_multi = "cp $seq_dir$poolName*.fastq.gz $tmp_fastq && ";
		
		my $config_csv = "$dir$poolName.csv";
		open (CSV,">$config_csv") or confess ("Can't open '$config_csv': $@");
		print CSV "[gene-expression]\n";
		print CSV "reference,".$index."\n";
		print CSV "probe-set,".$set."\n";
		print CSV "no-bam,true\n" unless ($create_bam);
		print CSV "chemistry,$chemistry\n" unless ($create_bam);
		print CSV "\n";
		print CSV "[libraries]\n";
		print CSV "fastq_id,fastqs,feature_types\n";
		print CSV $poolName.",".$tmp_fastq.",Gene Expression\n";
		print CSV "\n";
		print CSV "[samples]\n";
		print CSV "sample_id,probe_barcode_ids,description\n";
		foreach my $p (sort { $a->name cmp $b->name } @pat){
			print CSV $p->{name}.",".$p->{bc}."\n";
#			system("mkdir $dir$pname") unless (-d "$dir$pname");
		}
		close CSV;
		$cmd_multi .= "cd $tmp && cellranger multi --id=$poolName --csv=$config_csv ";
		$cmd_multi .= "&& cp -r $tmp$poolName/outs/* $tmp$poolName/_versions $tmp$poolName/_cmdline $dir$poolName/ ";
		print JOBS $cmd_multi."\n";
	}
	close JOBS;
	
	
	my $cmd2 = "cat $dir/jobs_multi.txt | run_cluster.pl -cpu=$cpu";
	if (grep(/count|^all$/i, @steps)){
		warn $cmd2;
		system $cmd2  unless $no_exec==1;
	}
}


#------------------------------
# AGGREGATION
#------------------------------
if (grep(/aggr/, @steps)){
	my $aggr_file = $dir."jobs_aggr.txt";
	my $id = $projectName.'_aggregation';
	$id = $aggr_name if $aggr_name;
	open (JOBS_AGGR, ">$aggr_file");
	my $type = $project->getRun->infosRun->{method};
	print JOBS_AGGR "sample_id,molecule_h5\n";
	foreach my $patient (sort {$a->name cmp $b->name} @$patients) {
		print JOBS_AGGR $patient->name().",".$dir."per_sample_outs/".$patient->name()."/molecule_info.h5\n";
	}
	close JOBS_AGGR;
	my $aggr_cmd = "cd $dir && $exec aggr --id=$id --csv=$aggr_file";
	warn $aggr_cmd;
	system ($aggr_cmd)  unless $no_exec;
}



#------------------------------
# ARCHIVE / TAR
#------------------------------
if (grep(/tar|archive|^all$/i, @steps)){
	my $tar_cmd = "tar -cvzf $dir$projectName.tar.gz $dir*/outs/per_sample_outs/*/web_summary.html "
		."$dir*/outs/per_sample_outs/*/count/sample_cloupe.cloupe "
		."$dir*/outs/per_sample_outs/*/count/*_bc_matrix/* ";
	$tar_cmd .= "$dir*/outs/per_sample_outs/*/count/*.bam* " if ($create_bam and $create_bam ne 'false');
	warn $tar_cmd;
	if (-e "$dir$projectName.tar.gz" and !$no_exec) {
		my $overwrite = prompt("'archive $dir$projectName.tar.gz' already exists. Overwrite ?  (y/n) ", -yes);
		die("archive '$dir$projectName.tar.gz' already exists") unless ($overwrite);
	}
	unless ($no_exec){
		system ($tar_cmd);
		print "\t------------------------------------------\n";
	#	print "\tlink to send to the users : \n";
	#	print "\twww.polyweb.fr/NGS/$projectName/$projectName.tar.gz \n";
		print "\tArchive to send to the users : \n";
		print "\t$dir$projectName.tar.gz\n" if (-e "$dir$projectName.tar.gz");
		print "\t------------------------------------------\n\n";
	}
}



#------------------------------
# COPY to SingleCell shared directory
#------------------------------
if (grep(/^cp(_web_summar(y|ies))?$|^all$/i, @steps)){
	my $dirout = "/data-pure/SingleCell/$projectName/";
	unless (-d $dirout) {
		my $cp_cmd = "mkdir $dirout";
		warn $cp_cmd;
		system ($cp_cmd) unless ($no_exec);
	}

#	my $hgroups ;
#	foreach my $patient (@{$patients}) {
#		$hgroups->{$patient->somatic_group}++;
#		my $cmd_cp = "mkdir $dir$pool_name"
#	}
#	if (grep(/^cp$/i, @steps)) {
#		foreach my $pool_name (sort keys %$hgroups){
#			warn  "$dir$pool_name";
#			my $cmd= qq{rsync -rav  $dir/$group_name $dirout && cp $dir/$group_name.csv $dirout };
#			warn $cmd;
#			system($cmd);
#			
#		}
#	}
	
	foreach my $patient (@{$patients}) {
		my $pool_name = $patient->somatic_group;
		my $pat_name = $patient->name;
		my $web_summary = "$dir$pool_name/per_sample_outs/$pat_name/web_summary.html";
		
		if (-e $web_summary or $no_exec) {
			my $cmd_cp;
			
			# Copy ALL outs
			if (grep(/^cp$/i, @steps)) { # cp all outs
				$cmd_cp = "cp -ru $dir/$pool_name $dirout";
			}
			
			# Copy ONLY web_summary.html
			else {
				$cmd_cp = "mkdir -p $dirout$pool_name/per_sample_outs/$pat_name " unless (-d "$dirout$pool_name/per_sample_outs/$pat_name");
				$cmd_cp .= " && cp $web_summary $dirout$pool_name/per_sample_outs/$pat_name/";
			}
			
			warn $cmd_cp;
			system ($cmd_cp) unless $no_exec;
		}
		else {
			warn ("$pool_name $pat_name web summary not found: '$web_summary'\nAre you sure the counting finished successfully ?");
		}
	}
}



exit;



sub usage {
	print "
$0
-------------
Obligatoires:
	project <s>                nom du projet
	feature_ref	<s>            tableau des ADT, obligatoire seulement si step=count et qu'il y a des ADT
Optionels:
	steps <s>                  étape à réaliser: demultiplex, count, tar, aggr, aggr_vdj, cp, cp_web_summary ou all (= demultiplex, count, cp_web_summary, tar)
	patients <s>               noms de patients/échantillons, séparés par des virgules
	cpu <i>                    nombre de cpu à utiliser, défaut: 20
	lane <i>                   nombre de lanes sur la flowcell, défaut: lit le RunInfo.xml
	mismatches <i>             nombre de mismatches autorisés lors du démultiplexage, défaut: 0
	create-bam/nocreate-bam    générer ou non les bams lors du count, défaut: nocreate-bam
	aggr_name <s>              nom de l'aggrégation, lors de step=aggr ou aggr_vdj
	chemistry                  chemistry , défaut: auto (pour librairies exp et adt)
	no_exec                    ne pas exécuter les commandes
	help                       affiche ce message

";
	exit(1);
}


