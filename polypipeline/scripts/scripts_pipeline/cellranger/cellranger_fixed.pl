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
my $feature_ref;
my $create_bam = 1;
my $no_exec;
my $aggr_name;
my $cmo_ref;
my $lane;


my $limit;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'step=s'=> \$step,
	'create_bam!'	=> \$create_bam,
	'feature_ref=s' =>  \$feature_ref,
	'cmo_ref=s' =>  \$cmo_ref,
	'no_exec=s' => \$no_exec,
	'aggr_name=s' => \$aggr_name,
	#'low_calling=s' => \$low_calling,
);

die ('-project option mandatory') unless ($projectName);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
#$patients_name = "all" unless $patients_name;
$patients_name = $project->get_only_list_patients($patients_name);
die("No patient in project $projectName") unless $patients_name;
die ("Project $projectName is not in somatic mode. Activate somatic mode and check that the groups have been filled in,\
 so that they can be taken into account in the analysis.") unless ($project->isSomatic or $step eq 'demultiplex' or $step eq 'teleport');

my $dir = $project->getProjectRootPath();
warn $dir;

unless ($step) {
	my $steps = ['demultiplex', 'teleport', 'count', 'aggr', 'aggr_vdj', 'tar', 'cp', 'all'];
	$step = prompt("Select a step: ", -menu=>$steps);
}
warn 'step='.$step;

my $hfamily;
my $prog;
foreach my $patient (@{$patients_name}) {
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

my $index = $project->getGenomeIndex($prog);
my $set = $index."/probe_sets/";
$set .= "Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv" if ($project->getVersion() =~ /^HG/);
$set .= "Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv" if ($project->getVersion() =~ /^MM/);


if ($step eq "all" || $step eq "demultiplex"){
	my $run = $project->getRun();
	my $run_name = $run->plateform_run_name;
	my $bcl_dir = $run->bcl_dir;
	my $tmp = $project->getAlignmentPipelineDir("cellranger_count");
	
	my $config = XMLin("$bcl_dir/RunInfo.xml", KeyAttr => { reads => 'Reads' }, ForceArray => [ 'reads', 'read' ]);
	my $lane = $config->{Run}->{FlowcellLayout}->{LaneCount};
	warn "LaneCount=$lane";
	
	my $sampleSheet = $bcl_dir."/sampleSheet.csv";
	open (SAMPLESHEET,">$sampleSheet");
	print SAMPLESHEET "Lane,Sample,Index\n";
	for (my $i=1;  $i<=$lane;$i++){
		foreach my $k (keys(%$hfamily)){
			my @pat = @{$hfamily->{$k}};
			print SAMPLESHEET $i.",".$pat[0]->{group}.",".$k."\n";
		}
	}
	close(SAMPLESHEET);
	
	
	my $cmd = "cd $tmp; /software/bin/demultiplex.pl -dir=$bcl_dir -run=$run_name -hiseq=10X -sample_sheet=$sampleSheet -cellranger_type=cellranger";
	warn $cmd;
	system $cmd or die "impossible $cmd" unless $no_exec==1;
}



if ($step eq "count" ){
	my $cmo_set;
	if ($cmo_ref) {
		die("'$cmo_ref' does not exists") unless (-e $cmo_ref);
		open(CMO, '<', $cmo_ref) or die("Can't open '$cmo_ref'");
		while (my $line = <CMO>) {
			chomp $line;
			my $cmo_id = (split /,/ , $line)[0];
			push (@$cmo_set, $cmo_id) unless ($cmo_id eq 'id');
		}
		close(CMO);
	}
	
	open (JOBS, ">$dir/jobs_count.txt");
	
	foreach my $k (keys(%$hfamily)){
		my @pat = @{$hfamily->{$k}};
		my $poolName = $pat[0]->{group};
		next if ($poolName =~ /^(ADT|CMO)_/);
		my $csv = $dir."/".$poolName.".csv";
		open (CSV,">$csv") or die ("Can't open '$csv");
		print CSV "[gene-expression]\n";
		print CSV "reference,".$index."\n";
		print CSV "cmo-set,".abs_path($cmo_ref)."\n" if ($cmo_ref);
		print CSV "probe-set,".$set."\n";
		print CSV "create-bam,true\n" if ($create_bam);
		print CSV "create-bam,false\n" unless ($create_bam);
		print CSV "\n";
		print CSV "[libraries]\n";
		my $pobj=$pat[0]->{objet};
		my $fastq = $pobj->getSequencesDirectory();
		print CSV "fastq_id,fastqs,feature_types\n";
		print CSV $poolName.",".$fastq.",Gene Expression\n";
		my @CMO = grep {$_->name eq 'CMO_'.$pobj->family} @$patients_name;
		@CMO = grep {$_->name eq 'ADT_'.$pobj->family} @$patients_name unless (scalar @CMO);
		die("No cmo/adt 'CMO_".$pobj->family."' or 'ADT_".$pobj->family."' found") unless (scalar @CMO);
		my $cmo = $CMO[0];
		print CSV $cmo->somatic_group.",".$cmo->getSequencesDirectory.",Multiplexing Capture\n";
		print CSV "\n";
		print CSV "[samples]\n";
		print CSV "sample_id,probe_barcode_ids\n";
		foreach my $p (@pat){
			die("Can't find bc2 '".$p->{bc}."' in CMO reference\n", Dumper $cmo_set) if (not grep(/^$p->{bc}$/, @$cmo_set) and $cmo_ref);
			print CSV $p->{name}.",".$p->{bc}."\n";
		}
		close CSV;
		print JOBS "cd $dir; cellranger multi --id=$poolName --csv=$csv\n";
	}
	close JOBS;
	
	my $cmd2 = "cat $dir/jobs_count.txt | run_cluster.pl -cpu=20";
	warn $cmd2;
	system $cmd2  unless $no_exec==1;
	
	
	# Open web summaries
	my @error;
	my $cmd3 = "firefox ";
	foreach my $patient (@{$patients_name}) {
		my $file = $dir."/".$patient->name."/outs/web_summary.html";
		push(@error, $file) unless (-e $file);
		$cmd3 .= "$file " if (-e $file);
	}
	warn $cmd3;
	system $cmd3 or die("Can't open web summaries") unless ($no_exec);
	die("File(s) not found: ".join(', ', @error)) if (@error);
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
	my $aggr_cmd = "cd $dir ; cellranger aggr --id=$id --csv=$aggr_file";
	warn $aggr_cmd;
	die();
	system ($aggr_cmd) unless $no_exec==1;
}




if ($step eq "tar" or $step eq "all"){
	my $tar_cmd = "tar -cvzf $dir/$projectName.tar.gz $dir/*/outs/per_sample_outs/*/web_summary.html $dir/*/outs/per_sample_outs/*/count/sample_cloupe.cloupe $dir/*/outs/per_sample_outs/*/count/*vloupe.vloupe $dir/*/outs/per_sample_outs/*/count/*_bc_matrix/* ";
	die ("archive $dir/$projectName.tar.gz already exists") if -e $dir."/".$projectName.".tar.gz";
	system ($tar_cmd) unless $no_exec==1;
	#or die "impossible $tar_cmd";
	print "\t##########################################\n";
	print "\tlink to send to the users : \n";
	print "\twww.polyweb.fr/NGS/$projectName/$projectName.tar.gz \n";
	print "\t##########################################\n";
}



# Commande pour copier sur /data-isilon/singleCell/
if ($step eq "cp" or $step eq "all"){
	my $dirout = "/data-isilon/SingleCell/$projectName";
	
	system("mkdir $dirout") unless -e $dirout;
	my $hgroups ;
	foreach my $patient (@{$patients_name}) {
		$hgroups->{$patient->somatic_group}++ 
	}
	foreach my $group_name (keys %$hgroups){
		my $cmd= qq{rsync -ra $dir/$group_name $dirout && cp $dir/$group_name.csv $dirout };
		warn $cmd;
		system($cmd);
	}
	print "\t##########################################\n";
	print "\tcp to $dirout \n";
	print "\t##########################################\n";
}

exit;
