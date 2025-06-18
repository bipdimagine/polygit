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

  
my $bin_cecile=qq{$Bin/scripts/scripts_db_polypipeline};
my $bin_script_pipeline = qq{$Bin/scripts/scripts_pipeline};


my $projectName;
my $patients_name;
my $step;
my $create_bam;
my $feature_ref;
my $no_exec;
my $aggr_name;
my $cmo_ref;
my $lane;
my $cpu = 20;


my $limit;
GetOptions(
	'project=s'		=> \$projectName,
	'patients=s'	=> \$patients_name,
	'step=s'		=> \$step,
	'create_bam!'	=> \$create_bam,
	'feature_ref=s'	=> \$feature_ref,
	'cmo_ref=s'		=> \$cmo_ref,
	'no_exec'		=> \$no_exec,
	'aggr_name=s'	=> \$aggr_name,
	'nb_lane=s'		=> \$lane,
	'cpu=i'			=> \$cpu,
	#'low_calling=s'	=> \$low_calling,
) || die("Error in command line arguments\n");


if ($step eq "all" || $step eq "demultiplex"){
#	die ("-bcl, -run and -nb_lane options are required") unless $bcl_dir || $run || $lane;
}
die("cpu must be in [1;40], given $cpu") unless ($cpu > 0 and $cpu <= 40);




my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
die ("Project '$projectName' is not a 10X fixed project. Use cellranger.pl instead.") unless ($project->getRun->infosRun->{method} eq 'fixed');
$patients_name = "all" unless $patients_name;
$patients_name = $project->get_only_list_patients($patients_name);
die("No patient in project $projectName") unless ($patients_name);

my $run = $project->getRun();
my $run_name = $run->plateform_run_name;
my $bcl_dir = $run->bcl_dir;
unless ($lane) {
	my $config = XMLin("$bcl_dir/RunInfo.xml", KeyAttr => { reads => 'Reads' }, ForceArray => [ 'reads', 'read' ]);
	$lane = $config->{Run}->{FlowcellLayout}->{LaneCount};
	warn 'LaneCount=',$lane;
}

my $sampleSheet = $bcl_dir."/sampleSheet.csv";
my $exec = "cellranger";

open (SAMPLESHEET,">$sampleSheet");
#my $dir = $project->getProjectRootPath();
my $dir = $project->getCountingDir('cellranger');
print SAMPLESHEET "Lane,Sample,Index\n";



my $fastq;
my %hSamples;
my %test;
my $type;
my @group;
my $hsamplesheet; 
my $prog;



my $hfamily;
foreach my $patient (@{$patients_name}) {
#	warn $patient->name;
	#my ($sbc,$pbc) = split(":",$patient->barcode);
	my $sbc = $patient->barcode; 
	my $pbc = $patient->barcode2;
	my $hpatient;
	$hpatient->{name}= $patient->name();
	$hpatient->{group}= $patient->somatic_group();
	$hpatient->{bc} = $pbc;
	$hpatient->{objet} = $patient;
	push(@{$hfamily->{$sbc}},$hpatient);
	$prog = $patient->alignmentMethod();
}
for (my $i=1;  $i<=$lane;$i++){
	foreach my $k (keys(%$hfamily)){
		my @pat = @{$hfamily->{$k}};
		print SAMPLESHEET $i.",".$pat[0]->{group}.",".$k."\n";
	}
}
close(SAMPLESHEET);


my $index = $project->getGenomeIndex($prog);
my $set = $index."/probe_sets/";
$set .= "Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv";
$set .= "Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv" if ($project->getVersion() =~ /^MM/);
my $tmp = $project->getAlignmentPipelineDir("cellranger_demultiplex");

my $cmd = "cd $tmp; /software/bin/demultiplex.pl -dir=$bcl_dir -run=$run_name -hiseq=10X -sample_sheet=$sampleSheet -cellranger_type=$exec";
if ($step eq "demultiplex" or $step eq "all"){
	my $exit = system ($cmd) unless ($no_exec);
	exit($exit) if ($exit);
	unless ($@ or $no_exec) {
		print("--------------------\n");
		print("Check the demultiplex stats\n");
		system ("firefox https://www.polyweb.fr/NGS/demultiplex/$run_name/$run_name\_laneBarcode.html &");
		print("https://www.polyweb.fr/NGS/demultiplex/$run_name/$run_name\_laneBarcode.html\n\n");
	}
}


my $tmp = $project->getAlignmentPipelineDir("cellranger_multi");
warn $tmp;

sub tmpSequencesDirectory {
	my $pat_name = shift @_;
	return "/tmp/pipeline/$projectName/".$pat_name.'/';
}
sub full_cmd {
	my ($p, $cmd) = @_;
	chomp $cmd;
	my $poolName = $p->{group};
	my $pat = $p->{objet};
	my $tmp_node = tmpSequencesDirectory($poolName);
	my $seq_dir = $pat->getSequencesDirectory;
	my $cmd1 = "mkdir -p $tmp_node && cp $seq_dir*$poolName*.fastq.gz $tmp_node && ";
	my $cmd2 = " --localcores=$cpu ";
	$cmd2 .= "&& mv $tmp$poolName/outs/* $dir$poolName/ ";
	system("mkdir $dir$poolName") unless (-d "$dir$poolName");
	$cmd2 .= "&& mv $dir$poolName/possorted_bam.bam $dir$poolName/possorted_bam.bam.bai $tmp$poolName/outs/ " if ($exec eq 'cellranger-atac' and $create_bam eq 'false');
	$cmd2 .= "&& if [[ \"$tmp_node\*.fastq.gz\" ]]; then rm $tmp_node\*.fastq.gz ; fi ";
	$cmd =~ s/$seq_dir/$tmp_node/;
	return $cmd1.$cmd.$cmd2."\n";
}



$create_bam = 'false' unless ($create_bam);
$create_bam = 'true' if ($create_bam ne 'false');
warn "create-bam=$create_bam";

open (JOBS, ">$dir/jobs_count.txt");

foreach my $k (keys(%$hfamily)){
	my @pat = @{$hfamily->{$k}};
	my $poolName = $pat[0]->{group};
	my $csv = $dir."/".$poolName.".csv";
	open (CSV,">$csv") or die ("putain");
	print CSV "[gene-expression]\n";
	print CSV "reference,".$index."\n";
	print CSV "probe-set,".$set."\n";
	print CSV "create-bam,".$create_bam."\n";
	print CSV "\n";
	print CSV "[libraries]\n";
	my $pobj=$pat[0]->{objet};
	my $fastq = tmpSequencesDirectory($poolName);
	print CSV "fastq_id,fastqs,feature_types\n";
	print CSV $poolName.",".$fastq.",Gene Expression\n";
	print CSV "\n";
	print CSV "[samples]\n";
	print CSV "sample_id,probe_barcode_ids\n";
	foreach my $p (@pat){
		print CSV $p->{name}.",";
		print CSV $p->{bc}."\n";
	}
	my $cmd = "cd $tmp && cellranger multi --id=$poolName --csv=$csv \n";
	$cmd = full_cmd($pat[0], $cmd);
	print JOBS $cmd;
	
}

close JOBS;
	
close CSV;


my $cmd2 = "cat $dir/jobs_count.txt | run_cluster.pl -cpu=$cpu";
if ($step eq "count" or $step eq "all"){
	warn $cmd2;
	my $exit = system ($cmd2) unless ($no_exec);
	exit($exit) if ($exit);
	
	# Open web summaries
	my @error;
	my $web_summaries;
	foreach my $patient (@{$patients_name}) {
		my $file = $dir."/".$patient->somatic_group."/per_sample_outs/".$patient->name."/web_summary.html";
		$web_summaries .= $file.' ' if (-e $file);
		push(@error, $file) unless (-e $file or $no_exec);
	}
	my $cmd3 = "firefox ".$web_summaries.' &';
	warn $cmd3 unless ($web_summaries);
	system($cmd3) unless ($web_summaries or $no_exec);
	die("Web summaries not found: ".join(', ', @error)) if (@error and not $no_exec);

	unless ($no_exec) {
		print("--------------------\n");
		print("Check the web summaries:\n");
		print("$dir/*/web_summary.html\n\n");
	}
}

if ($step eq "aggr" or $step eq "all"){
	my $aggr_file = $dir."/jobs_aggr.txt";
	my $id = $projectName;
	$id = $aggr_name if $aggr_name;
	open (JOBS_AGGR, ">$aggr_file");
	my $type = $project->getRun->infosRun->{method};
	if ($type =~ /vdj/) {
		print JOBS_AGGR "sample_id,vdj_contig_info,donor,origin\n" ;
	}
	else{
		print JOBS_AGGR "sample_id,molecule_h5\n";
	}	
	foreach my $patient (@{$patients_name}) {
		my $group_type = uc($patient->somatic_group());
		if ($group_type eq "VDJ") {
			print JOBS_AGGR $patient->name().",".$dir."/".$patient->name()."/"."outs/vdj_contig_info.pb\n";
		}
		else{	
			print JOBS_AGGR $patient->name().",".$dir."/".$patient->name()."/"."outs/molecule_info.h5\n";
		}
	}
	close JOBS_AGGR;
	my $aggr_cmd = "cd $dir ; $exec aggr --id=$id --csv=$aggr_file";
	warn $aggr_cmd;
	die();
	my $exit = system ($aggr_cmd) unless ($no_exec);
	exit($exit) if ($exit);
}


#foreach my $patient (@{$patients_name}) {
#	my $poolName = $patient->name();
#	my $file = $dir."/".$poolName."/outs/web_summary.html";
#	print $file."\n" if -e $file ;
#}

##commande pour copier sur /data-isilon/singleCell
#
if ($step eq "tar" or $step eq "all"){
	my $tar_cmd = "tar -cvzf $dir$projectName.tar.gz $dir*/per_sample_outs/*/web_summary.html $dir*/multi/count/raw_feature_bc_matrix "
				 ."$dir*/multi/count/raw_feature_bc_matrix/* $dir*/per_sample_outs/*/count/sample_cloupe.cloupe  "
				 ."$dir*/per_sample_outs/*/count/sample_*_feature_bc_matrix/* ";
#	my $tar_cmd = "cd $dir && tar -cvzf $dir$projectName.tar.gz */per_sample_outs/*/web_summary.html */multi/count/raw_feature_bc_matrix "
#				 ."*/multi/count/raw_feature_bc_matrix/* */per_sample_outs/*/count/sample_cloupe.cloupe  */per_sample_outs/*/count/sample_*_feature_bc_matrix/* ";
	warn $tar_cmd;
#	die ("archive '$dir$projectName.tar.gz' already exists") if -e "$dir$projectName.tar.gz";
	my $exit = system ($tar_cmd) unless ($no_exec);
	exit($exit) if ($exit);
	print "\t########################################\n";
	print "\t $dir$projectName.tar.gz\n";
	print "\t########################################\n";
}
#
if ($step eq "cp" or $step eq "all"){
	my $dirout = "/data-isilon/SingleCell/$projectName";
	
	system("mkdir $dirout") unless -e $dirout;
	my $hgroups ;
	foreach my $patient (@{$patients_name}) {
		$hgroups->{$patient->somatic_group}++ 
	}
	foreach my $group_name (keys %$hgroups){
		warn  "$dir/$group_name";
		my $cmd= "cp $dir/$group_name.csv $dirout";
		warn $cmd;
		my $exit = system ($cmd) unless ($no_exec);
		exit($exit) if ($exit);
	}
	print "\t########################################\n";
	print "\t cp to $dirout\n";
	print "\t########################################\n";
}

exit;
