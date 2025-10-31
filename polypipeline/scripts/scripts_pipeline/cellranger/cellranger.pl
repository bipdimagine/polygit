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
use Time::Local 'timelocal';

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
my $lane;
<<<<<<< HEAD
=======
my $mismatch = 0;
>>>>>>> branch 'dev-mperin' of https://github.com/bipdimagine/polygit.git
my $feature_ref;
my $no_exec;
my $aggr_name;
my $chemistry;
my $create_bam;
my $limit;
my $help;

GetOptions(
	'project=s'		=> \$projectName,
	'patients=s'	=> \$patients_name,
	'step=s'		=> \$step,
	'lane=i'		=> \$lane,
<<<<<<< HEAD
=======
	'mismatches=i'	=> \$mismatch,
>>>>>>> branch 'dev-mperin' of https://github.com/bipdimagine/polygit.git
	'create_bam!'	=> \$create_bam,
	'feature_ref=s'	=> \$feature_ref,
	'aggr_name=s'	=> \$aggr_name,
	'chemistry=s'	=> \$chemistry,
	'no_exec'		=> \$no_exec,
	'limit=i'		=> \$limit,
#	'low_calling=s'	=> \$low_calling,
	'help'			=> \$help,
);

usage() if $help;
<<<<<<< HEAD
#die ('-project option required') unless $projectName;
=======
>>>>>>> branch 'dev-mperin' of https://github.com/bipdimagine/polygit.git
usage() unless ($projectName);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
$patients_name = $project->get_only_list_patients($patients_name);
die("No patient in project $projectName") unless $patients_name;

my $run = $project->getRun();
my $run_name = $run->plateform_run_name;
my $type = $run->infosRun->{method};
my $machine = $run->infosRun->{machine};
<<<<<<< HEAD
#die ("Error in sequencing machine: '$machine', expected '10X'") unless ($machine eq '10X');
=======
>>>>>>> branch 'dev-mperin' of https://github.com/bipdimagine/polygit.git
unless ($machine eq '10X') {
	my $continue = prompt( "Error in sequencing machine: '$machine', expected '10X'. Continue anyway ?  (y/n)  ", -yes_no );
	die  unless ($continue);
}

my $exec = "cellranger";
$exec .= '-atac' if ($type eq 'atac');
$exec .= '-arc' if ($type eq 'arc');
$exec = 'spaceranger' if  ($type eq 'spatial');
warn $exec;
	

my $dir = $project->getProjectRootPath();
warn $dir;
#my $tmp = $project->getAlignmentPipelineDir("cellranger_count");
#warn $tmp;

unless ($step) {
	my $steps = ['demultiplex', 'teleport', 'count', 'aggr', 'aggr_vdj', 'tar', 'cp', 'cp_web_summaries', 'all'];
	$step = prompt("Select a step: ", -menu=>$steps);
}
warn 'step='.$step;
<<<<<<< HEAD

if ($limit) {
	$limit = scalar @$patients_name + $limit if ($limit <= 0);
	$limit = 1 if ($limit <= 0);
	warn 'limit='.$limit;
}
=======
>>>>>>> branch 'dev-mperin' of https://github.com/bipdimagine/polygit.git

if ($limit) {
	$limit = scalar @$patients_name + $limit if ($limit <= 0);
#	$limit = 1 if ($limit <= 0);
	warn 'limit='.$limit;
}

<<<<<<< HEAD
if ($step eq "demultiplex" || $step eq "all"){
	
	warn $exec;
	my $bcl_dir = $run->bcl_dir;
	warn $bcl_dir;
	unless ($lane) {
		my $config = XMLin("$bcl_dir/RunInfo.xml", KeyAttr => { reads => 'Reads' }, ForceArray => [ 'reads', 'read' ]);
		$lane = $config->{Run}->{FlowcellLayout}->{LaneCount};
		warn 'LaneCount=',$lane;
	}
	
	my $tmp = $project->getAlignmentPipelineDir("cellranger_demultiplex");
	warn $tmp;

	my $sampleSheet = $bcl_dir."/sampleSheet.csv";
	open (SAMPLESHEET,">$sampleSheet");
	print SAMPLESHEET "Lane,Sample,Index\n";
	
	foreach my $patient (@{$patients_name}) {
		my $name = $patient->name();
		my $bc = $patient->barcode();
		my $bc2 = $patient->barcode2();
		for (my $i=1;  $i<=$lane;$i++){
			print SAMPLESHEET $i.",".$name.",".$bc."\n";
		}
	}
	close(SAMPLESHEET);
	
	my $cmd = "cd $tmp; $Bin/../../demultiplex/demultiplex.pl -dir=$bcl_dir -run=$run_name -hiseq=10X -sample_sheet=$sampleSheet -cellranger_type=$exec -mismatch=1";
	warn $cmd;
	system $cmd unless $no_exec;
#	system $cmd or die "Impossible $cmd" unless $no_exec;
	unless ($@ and $no_exec) {
=======
###############
# DEMULIPLEXAGE
###############
if ($step eq "demultiplex" || $step eq "all"){
	
	warn $exec;
	my $bcl_dir = $run->bcl_dir;
	warn $bcl_dir;
	unless ($lane) {
		my $config = XMLin("$bcl_dir/RunInfo.xml", KeyAttr => { reads => 'Reads' }, ForceArray => [ 'reads', 'read' ]);
		$lane = $config->{Run}->{FlowcellLayout}->{LaneCount};
		warn 'LaneCount=',$lane;
	}
	
	my $tmp = $project->getAlignmentPipelineDir("cellranger_demultiplex");
	warn $tmp;

	my $sampleSheet = $bcl_dir."/sampleSheet.csv";
	open (SAMPLESHEET,">$sampleSheet");
	print SAMPLESHEET "Lane,Sample,Index\n";
	
	foreach my $patient (@{$patients_name}) {
		my $name = $patient->name();
		my $bc = $patient->barcode();
		my $bc2 = $patient->barcode2();
		for (my $i=1;  $i<=$lane;$i++){
			print SAMPLESHEET $i.",".$name.",".$bc."\n";
		}
	}
	close(SAMPLESHEET);
	
	my $cmd = "cd $tmp; $Bin/../../demultiplex/demultiplex.pl -dir=$bcl_dir -run=$run_name -hiseq=10X -sample_sheet=$sampleSheet -cellranger_type=$exec -mismatch=$mismatch";
	warn $cmd;
	system $cmd unless $no_exec;
	unless ($@ or $no_exec) {
>>>>>>> branch 'dev-mperin' of https://github.com/bipdimagine/polygit.git
		colored::stabilo('white', "Done");
		colored::stabilo('white', "Check the demultiplex stats");
		colored::stabilo('white', "https://www.polyweb.fr/NGS/demultiplex/$run_name/$run_name\_laneBarcode.html");
	}
}


<<<<<<< HEAD
=======
###############
# TELEPORT
###############
>>>>>>> branch 'dev-mperin' of https://github.com/bipdimagine/polygit.git
if ($step eq 'teleport') {
	my $cmd = "teleport_patient.sh -project=$projectName -force=1";
	system($cmd) unless $no_exec;
}

unless ($project->isSomatic) {
	my $warn = "/!\\ Project $projectName is not in somatic mode. Activate somatic mode and check that the groups have been filled in, "
		."so that they can be taken into account in the analysis.";
	die ($warn) unless (grep{/$step/} ('demultiplex', 'teleport', 'tar', 'cp', 'cp_web_summaries'));
	warn ($warn);
	
}

<<<<<<< HEAD
die ("Project $projectName is not in somatic mode. Activate somatic mode and check that the groups have been filled in, so that they can be taken into account in the analysis.") unless ($project->isSomatic or $step eq 'demultiplex' or $step eq 'teleport');
=======
>>>>>>> branch 'dev-mperin' of https://github.com/bipdimagine/polygit.git
my @group;
foreach my $patient (@{$patients_name}) {
	my $group = $patient->somatic_group();
	push(@group,$group);
}

<<<<<<< HEAD
my @type = split /_/, $type;
for my $i ($#type) {
	warn ("No patient with group corresponding to sequencing method '$type'") unless (grep {/$type[$i]/i} @group);
}

if ($step eq "count" || $step eq "all"){
	
	if (grep {$_ =~ /adt/i} @group) {
		die("feature_ref csv required\n") unless ($feature_ref);
		die("'$feature_ref' not found") unless (-e $feature_ref);
	}
	
	my $fastq;
	my %hSamples;
	my %test;
	
	$create_bam = 'false' unless ($create_bam);
	$create_bam = 'true' if ($create_bam ne 'false');
	warn "create-bam=$create_bam";
	
	 if ($chemistry) {
 		warn "chemistry=$chemistry";
		my @chemistries = ('auto','threeprime','fiveprime','SC3Pv2','SC3Pv3','SC3Pv4','SC3Pv3HT','SC5P-PE','SC5P-R2','SC3Pv1','ARC-v1');
		die ("Chemistry option '$chemistry' not valid, should be one of: ". join(', ', @chemistries)) unless (grep { $_ eq $chemistry } @chemistries);
	 }
	
	
	
	# EXP
	my $type_exp = 1 if map {uc($_) =~ /EXP|NUCLEI/ } @group;
	if($type_exp){
		open (JOBS, ">$dir/jobs_count.txt");
		my @exp = grep { uc($_->somatic_group()) eq "EXP"} @{$patients_name};
		foreach my $e (@exp){
			my $name = $e->name;
			warn $name;
=======
unless ($type eq 'spatial') {
	my @type = split /_/, $type;
	warn Dumper \@group;
	for my $i (0..$#type) {
		warn ("No patient with group corresponding to sequencing method '$type'") unless (grep {/$type[$i]/i} @group);
	}
}

###############
# COUNT
###############
if ($step eq "count" || $step eq "all"){
	
	if (grep {$_ =~ /adt/i} @group) {
		die("feature_ref csv required\n") unless ($feature_ref);
		die("'$feature_ref' not found") unless (-e $feature_ref);
	}
	
	my $fastq;
	my %hSamples;
	my %test;
	
	$create_bam = 'false' unless ($create_bam);
	$create_bam = 'true' if ($create_bam ne 'false');
	warn "create-bam=$create_bam";
	
	 if ($chemistry) {
 		warn "chemistry=$chemistry";
		my @chemistries = ('auto','threeprime','fiveprime','SC3Pv2','SC3Pv3','SC3Pv4','SC3Pv3HT','SC5P-PE','SC5P-R2','SC3Pv1','ARC-v1');
		die ("Chemistry option '$chemistry' not valid, should be one of: ". join(', ', @chemistries)) unless (grep { $_ eq $chemistry } @chemistries);
	 }
	
	
	
	# EXP
	my $type_exp = 1 if map {uc($_) =~ /EXP|NUCLEI/ } @group;
	if($type_exp){
		open (JOBS, ">$dir/jobs_count.txt");
		my @exp = grep { uc($_->somatic_group()) eq "EXP"} @{$patients_name};
		foreach my $e (@exp){
			my $name = $e->name;
#			warn $name;
>>>>>>> branch 'dev-mperin' of https://github.com/bipdimagine/polygit.git
			my $group = $e->somatic_group();
			$fastq = $e->getSequencesDirectory();
			my $prog =  $e->alignmentMethod();
			my $index = $project->getGenomeIndex($prog);
			my $cmd = "cd $dir; $exec count --id=$name --sample=$name --fastqs=$fastq --create-bam=$create_bam --transcriptome=$index ";
			$cmd .= " --include-introns " if ($type eq "nuclei" or lc($group) eq "nuclei");
			$cmd .= " --chemistry $chemistry" if ($chemistry);
			$cmd .= "\n";
	#		warn $cmd;
			print JOBS $cmd;
		}
	}
	
	
	# VDJ
	my $type_vdj = 1 if map {uc($_) =~ /VDJ/ } @group;
	if($type_vdj){
		open (JOBS_VDJ, ">$dir/jobs_vdj.txt");
		my @vdj = grep { uc($_->somatic_group()) eq "VDJ"} @{$patients_name};
		foreach my $v (@vdj){
			my $name = $v->name(); 
#			my $vfam = $v->family();
#			my $vgroup = uc($v->somatic_group());
			my $fastq = $v->getSequencesDirectory();
			my $prog =  $v->alignmentMethod();
			my $index = $project->getGenomeIndex($prog);
			my $index_vdj = $index."_vdj";
			print JOBS_VDJ "cd $dir; cellranger vdj --sample=$name --id=$name --fastqs=$fastq --reference=$index_vdj  \n"
		}
	}


	# ADT
	my $type_adt = 1 if map {uc($_) =~ /ADT/ } @group;
	if ($type_adt){
		open (JOBS_ADT, ">$dir/jobs_count.txt");
		my @exp = grep { uc($_->somatic_group()) eq "EXP"} @{$patients_name};
		foreach my $e(@exp){
			my $ename = $e->name(); 
			my $efam = $e->family();
#			my $egroup = uc($e->somatic_group());
			my $lib_file = $dir."/".$ename."_library.csv";
			open(LIB,">$lib_file") or die "impossible $lib_file";
			my @adt = grep {$_->family() eq $efam && $_->name() =~ $ename && uc($_->somatic_group()) eq "ADT"} @{$patients_name};
			warn ("no associated ADT library") if scalar(@adt)==0 ;
			next() if scalar(@adt)==0 ; 
			my $adt_name = $adt[0]->name();
			my $lib = "fastqs,sample,library_type\n".$e->getSequencesDirectory().",".$ename.",Gene Expression\n";
			$lib .= $adt[0]->getSequencesDirectory().",".$adt_name.",Antibody Capture\n";
			print LIB $lib;
			close(LIB);
			my $prog =  $e->alignmentMethod();
			my $index = $project->getGenomeIndex($prog);
			print JOBS_ADT "cd $dir ; cellranger count --id=$ename --feature-ref=$feature_ref --transcriptome=$index  --libraries=$lib_file --create-bam=$create_bam \n"
		}
	}

	
	# SPATIAL
	my $type_spatial = 1 if map {uc($_) =~ /SPATIAL/ } @group;
 	if($type_spatial  && $step eq "count" ){
		open (JOBS_SPATIAL, ">$dir/jobs_spatial.txt");
		my @spatial = grep { uc($_->somatic_group()) eq "SPATIAL"} @{$patients_name};
		foreach my $s (@spatial){
			my $sname = $s->name(); 
			my $bc2 = $s->barcode2();
			my ($slide,$slide2,$area) = split("-",$bc2);
#			warn $sname;
			my $sfam = $s->family();
			#my $des_file = $dir."/".$sname."_spatial_descript.csv";
	
			#die("file with area, slide and path to image file is mandatory \(patientName_spatial_descript.csv\)") unless -e $des_file ;
			#open (DES, $des_file);
			#my $json;
			#my $area;
			#my $slide;
	
			#while(<DES>){
			#	warn $_;
			#	chomp($_);
			#	($json,$slide,$area)= split(",",$_);
			#}
			#close(DES);
			my $slide_final = $slide."-".$slide2;
			my $json = $dir."/".$slide_final."-".$area."-".$sname.".json";
			my $image = $dir."/".$area."-".$sname.".tif";
			my $sgroup = uc($s->somatic_group());
			my $fastq = $s->getSequencesDirectory();
			my $prog =  $s->alignmentMethod();
			my $index = $project->getGenomeIndex($prog);
			my $index_spatial = $index;
			print JOBS_SPATIAL "cd $dir ; spaceranger count --id=$sname --sample=$sname --image=$image --fastqs=$fastq --transcriptome=$index_spatial --area=$area --slide=$slide_final --loupe-alignment=$json --create-bam=$create_bam \n";
		}
	}

	
	# ATAC
	my $type_atac = 1 if map {uc($_) =~ /ATAC/ } @group;
#	warn $type_atac;
	if($type_atac ){
		warn 'ATAC';
		open (JOBS_ATAC, ">$dir/jobs_atac.txt");
		my @atac = grep { uc($_->somatic_group()) eq "ATAC"} @{$patients_name};
		foreach my $v (@atac){
			my $vname = $v->name(); 
			my $vfam = $v->family();
			my $vgroup = uc($v->somatic_group());
			my $fastq = $v->getSequencesDirectory();
			$exec = "cellranger-atac" if $type eq "atac";
			my $prog =  $v->alignmentMethod();
			my $index = $project->getGenomeIndex($prog);
			my $index_atac = $index."_atac" if $type eq "atac";
			print JOBS_ATAC "cd $dir; $exec count --sample=$vname --id=$vname --fastqs=$fastq --reference=$index_atac  \n"
		}
	}

	
<<<<<<< HEAD
	# CMO ?
=======
	# CMO
>>>>>>> branch 'dev-mperin' of https://github.com/bipdimagine/polygit.git
	my $type_cmo = 1 if map {uc($_) =~ /CMO/ } @group;
	if ($type =~ /cmo/ ){
		open (JOBS_CMO, ">$dir/jobs_cmo.txt");
#		warn $patient->somatic_group();
		my @exp = grep { uc($_->somatic_group()) eq "EXP"} @{$patients_name};
		foreach my $e(@exp){
			my $ename = $e->name(); 
			my $efam = $e->family();
			my $egroup = uc($e->somatic_group());
			my $prog =  $e->alignmentMethod();
			my $index = $project->getGenomeIndex($prog);
			my $lib_file = $dir."/".$ename."_multi.csv";
			open(LIB,">$lib_file") or die "impossible $lib_file";
			print LIB "[gene-expression]\nreference,".$index."\ncmo-set,/data-isilon/sequencing/ngs/NGS2022_6140/cmo_ref.csv\n";
			
			print LIB "[libraries]\nfastq_id,fastqs,feature_types\n";
			print LIB "$ename,".$e->getSequencesDirectory().$ename.",Gene Expression\n";
			my @cmo = grep {$_->family() eq $efam && uc($_->somatic_group()) eq "CMO"} @{$patients_name};
			my $cmo_name = $cmo[0]->name();
			print LIB "$cmo_name,".$cmo[0]->getSequencesDirectory().$ename.",Antibody Capture\n";
			print LIB "[samples]\nsample_id,cmo_ids\n";
			print LIB $ename."_B251,B251\n";
			print LIB $ename."_B252,B252\n";
			print LIB $ename."_B253,B253\n";
			close(LIB);
			print JOBS_CMO "cd $dir ; cellranger multi --id=$ename --csv=$lib_file\n";
		}
	}

	
	# ARC
	my $type_arc = 1 if map {uc($_) =~ /ARC/ } @group;
	if ($type_arc ){
		open (JOBS_ARC, ">$dir/jobs_arc.txt");
		$exec = "cellranger-arc";
		my @exp = grep { uc($_->somatic_group()) eq "EXP"} @{$patients_name};
		foreach my $e(@exp){
			my $ename = $e->name(); 
			my $efam = $e->family();
			my $egroup = uc($e->somatic_group());
			my $lib_file = $dir."/".$ename."_library.csv";
			open(LIB,">$lib_file") or die "impossible $lib_file";
			my @atac = grep {$_->family() eq $efam && uc($_->somatic_group()) eq "ARC"} @{$patients_name};
			next() if scalar(@atac == 0);
			my $atac_name = $atac[0]->name() ;
			my $prog =  $e->alignmentMethod();
			my $index = $project->getGenomeIndex($prog);
			my $index_arc = $index."_arc";
			my $lib = "fastqs,sample,library_type\n".$e->getSequencesDirectory().",".$ename.",Gene Expression\n";
			$lib .= $atac[0]->getSequencesDirectory().",".$atac_name.",Chromatin Accessibility\n";
			print LIB $lib;
			close(LIB);
			print JOBS_ARC "cd $dir ; $exec count --id=$ename --transcriptome=$index_arc  --libraries=$lib_file\n"
		}
	}
		
	close(JOBS);
	close(JOBS_VDJ);
	close(JOBS_ADT);
	close(JOBS_SPATIAL);
	close(JOBS_ATAC);
	close(JOBS_CMO);
	close(JOBS_ARC);


	my $cmd2 = "cat $dir/jobs*.txt | run_cluster.pl -cpu=20";
	$cmd2 .= " -limit=$limit" if ($limit);
	warn $cmd2;
	system ($cmd2) unless $no_exec; # || die("Can't run count") unless $no_exec;
	
<<<<<<< HEAD
	
#	colored::stabilo('white', "Done") unless $no_exec;
#	colored::stabilo('white', "Check the web summaries") unless $no_exec;
#	colored::stabilo('white', "$dir/*/outs/web_summary.html") unless $no_exec;
=======
#	unless ($no_exec) {
#		colored::stabilo('white', "Done");
#		colored::stabilo('white', "Check the web summaries");
#		colored::stabilo('white', "$dir/*/outs/web_summary.html");
#	}
>>>>>>> branch 'dev-mperin' of https://github.com/bipdimagine/polygit.git
	
	my @error;
	my $cmd3 = "firefox ";
	foreach my $patient (@{$patients_name}) {
		next if ($patient->somatic_group =~ /^ADT$/i);
		my $file = $dir."/".$patient->name."/outs/web_summary.html";
		$cmd3 .= $file.' ' if (-e $file);
		push(@error, $file) unless (-e $file or $no_exec);
<<<<<<< HEAD
=======
	}
	warn $cmd3 unless ($cmd3 eq "firefox ");
	system($cmd3) unless ($cmd3 eq "firefox " or $no_exec);
	die("Web summaries not found: ".join(', ', @error)) if (@error and not $no_exec);
}



###############
# AGGREGATION
###############
if ($step eq "aggr"){
	my @groups = map {$_->somatic_group} @$patients_name;
	warn ("Can't aggregate gene expression and vdj librairies together") if (grep {'exp'} @groups and grep {'vdj'} @groups);
#	die("No 'exp' sample to aggregate") unless (grep {'exp'} @groups);
	my $id = $aggr_name if $aggr_name;
	my $aggr_file = $dir."/jobs_aggr.txt";
	open (JOBS_AGGR, ">$aggr_file");
	
	if (grep {'exp'} @groups) {
		$id = 'aggregation_exp' unless ($id);
		my $aggr_csv = $dir."/aggr.csv";
		open (AGGR_CSV, ">$aggr_csv");
		print AGGR_CSV "sample_id,molecule_h5\n";
		foreach my $patient (@{$patients_name}) {
			warn 
			my $group_type = lc($patient->somatic_group());
			if ($group_type eq "exp") {	# ne "adt" or $_ ne "vdj"
				print AGGR_CSV $patient->name().",".$dir."/".$patient->name()."/outs/molecule_info.h5\n";
			}
		}
		close AGGR_CSV;
		print JOBS_AGGR "cd $dir ; $exec aggr --id=$id --csv=$aggr_csv";
		warn ("cat $aggr_file | run_cluster.pl -cpu=20");
		system ("cat $aggr_file | run_cluster.pl -cpu=20") or die("Can't run aggragation") unless $no_exec;
		die("Error while running cellranger aggr") unless (-d $id);
		colored::stabilo('white', "Done") if (-d "$dir/$id/" and not $no_exec);
		colored::stabilo('white', "$dir/$id/") if (-d "$dir/$id/" and not $no_exec);
	}
	
	if (grep {'vdj'} @groups) {
		$id = 'VDJ_'.$id if ($id);
		$id = 'aggregation_vdj' unless ($id);
		my $aggr_csv_vdj = $dir."/aggr_vdj.csv";
		if (-e $aggr_csv_vdj) {
			my $overwrite = prompt("'aggr_vdj.csv' already exists. Overwrite the file ?  ",-yes);
			print ("If 'aggr_vdj.csv' is completed, you can run 'echo \"cd $dir ; $exec aggr --id=$id\_VDJ --csv=$aggr_csv_vdj\" | run_cluster.pl -cpu=20'") unless ($overwrite);
			print JOBS_AGGR "cd $dir ; $exec aggr --id=$id --csv=$aggr_csv_vdj" unless ($overwrite);
			die("'$aggr_csv_vdj' already exists") unless ($overwrite);
		}
		open (AGGR_CSV_VDJ, ">$aggr_csv_vdj");
		print AGGR_CSV_VDJ "sample_id,vdj_contig_info,donor,origin\n" ;
		foreach my $patient (@{$patients_name}) {
			my $group_type = lc($patient->somatic_group());
			print AGGR_CSV_VDJ $patient->name().",$dir".$patient->name()."/outs/vdj_contig_info.pb,,\n" if ($group_type eq "vdj");
		}
		close (AGGR_CSV_VDJ);
		colored::stabilo('white', "Fill the 'donor' and 'origin' columns for each vdj sample in '$aggr_csv_vdj'.");
		colored::stabilo('white', "Then run 'echo \"cd $dir ; $exec aggr --id=$id --csv=$aggr_csv_vdj\" | run_cluster.pl -cpu=20'");
		colored::stabilo('white', "Then make an archive: 'tar -czf $dir/$projectName\_aggr.tar.gz $dir/aggregation_*/outs/web_summary.html $dir/aggregation_*/outs/count/cloupe.cloupe $dir/aggregation_*/outs/count/*_bc_matrix/* $dir/aggregation_*/outs/*/vloupe.vloupe'");
	}
	close (JOBS_AGGR);
	my $cmd_tar = "tar -czf $dir/$projectName\_aggr.tar.gz $dir/aggregation_*/outs/web_summary.html $dir/aggregation_*/outs/count/cloupe.cloupe $dir/aggregation_*/outs/count/*_bc_matrix/* $dir/aggregation_*/outs/*/vloupe.vloupe";
	colored::stabilo('white', "Then, you can make an archive: '$cmd_tar'");
#	system($cmd_tar) if (-d $dir.'aggregation_exp/' and not $no_exec);
#	colored::stabilo('white', "Archive of aggragation : $dir/$projectName\_aggr.tar.gz") if (-d "$dir/$projectName\_aggr.tar.gz" and not $no_exec);

}




###############
# AGGREGATION VDJ
###############
if ($step eq "aggr_vdj") {
	my $type = $project->getRun->infosRun->{method};
	die ("No vdj in project $projectName") if ($type !~ /vdj/);
	my @groups = map {$_->somatic_group} @$patients_name;
	die("No 'vdj' sample to aggregate") unless (grep {'vdj'} @groups);
	my $id = $projectName.'_VDJ_aggregation';
	$id = $aggr_name if $aggr_name;
#	my $aggr_file_vdj = $dir."/jobs_aggr_vdj.txt";
#	open (JOBS_AGGR_VDJ, ">$aggr_file_vdj");
	my $aggr_csv_vdj = $dir."/aggr_vdj.csv";
	if (-e $aggr_csv_vdj) {
		my $overwrite = prompt("'aggr_vdj.csv' already exists. Overwrite the file ?  ",-yes);
		print ("If 'aggr_vdj.csv' is completed, you can run 'echo \"cd $dir ; $exec aggr --id=$id --csv=$aggr_csv_vdj\" | run_cluster.pl -cpu=20'") unless ($overwrite);
		die("'$aggr_csv_vdj' already exists") unless ($overwrite);
	}
	open (AGGR_CSV_VDJ, ">$aggr_csv_vdj");
	print AGGR_CSV_VDJ "sample_id,vdj_contig_info,donor,origin\n" ;
	foreach my $patient (@{$patients_name}) {
		my $group_type = lc($patient->somatic_group());
		if ($group_type eq "vdj") {
			print AGGR_CSV_VDJ $patient->name().",".$dir."/".$patient->name()."/outs/vdj_contig_info.pb,,\n";
		}
>>>>>>> branch 'dev-mperin' of https://github.com/bipdimagine/polygit.git
	}
<<<<<<< HEAD
	warn $cmd3 unless ($cmd3 eq "firefox ");
	system($cmd3) unless ($cmd3 eq "firefox " or $no_exec);
	die("Web summaries not found: ".join(', ', @error)) if (@error and not $no_exec);
=======
	close (AGGR_CSV_VDJ);
#	print JOBS_AGGR_VDJ "cd $dir ; $exec aggr --id=$id\_VDJ --csv=$aggr_csv_vdj";
#	close (JOBS_AGGR_VDJ);
	colored::stabilo('white', "Fill the 'donor' and 'origin' columns for each vdj sample in '$aggr_csv_vdj'.");
	colored::stabilo('white', "Then run 'echo \"cd $dir ; $exec aggr --id=$id\_VDJ --csv=$aggr_csv_vdj\" | run_cluster.pl -cpu=20'");
	colored::stabilo('white', "Then make an archive: 'tar -czf $dir/$projectName\_aggr.tar.gz $dir/aggregation_*/outs/web_summary.html $dir/aggregation_*/outs/count/cloupe.cloupe $dir/aggregation_*/outs/count/*_bc_matrix/* $dir/aggregation_*/outs/*/vloupe.vloupe'");
>>>>>>> branch 'dev-mperin' of https://github.com/bipdimagine/polygit.git
}



<<<<<<< HEAD
if ($step eq "aggr"){
	my @groups = map {$_->somatic_group} @$patients_name;
	warn ("Can't aggregate gene expression and vdj librairies together") if (grep {'exp'} @groups and grep {'vdj'} @groups);
#	die("No 'exp' sample to aggregate") unless (grep {'exp'} @groups);
	my $id = $aggr_name if $aggr_name;
	my $aggr_file = $dir."/jobs_aggr.txt";
	open (JOBS_AGGR, ">$aggr_file");
	
	if (grep {'exp'} @groups) {
		$id = 'aggregation_exp' unless ($id);
		my $aggr_csv = $dir."/aggr.csv";
		open (AGGR_CSV, ">$aggr_csv");
		print AGGR_CSV "sample_id,molecule_h5\n";
		foreach my $patient (@{$patients_name}) {
			warn 
			my $group_type = lc($patient->somatic_group());
			if ($group_type eq "exp") {	# ne "adt" or $_ ne "vdj"
				print AGGR_CSV $patient->name().",".$dir."/".$patient->name()."/outs/molecule_info.h5\n";
			}
		}
		close AGGR_CSV;
		print JOBS_AGGR "cd $dir ; $exec aggr --id=$id --csv=$aggr_csv";
		warn ("cat $aggr_file | run_cluster.pl -cpu=20");
		system ("cat $aggr_file | run_cluster.pl -cpu=20") or die("Can't run aggragation") unless $no_exec;
		die("Error while running cellranger aggr") unless (-d $id);
		colored::stabilo('white', "Done") if (-d "$dir/$id/" and not $no_exec);
		colored::stabilo('white', "$dir/$id/") if (-d "$dir/$id/" and not $no_exec);
	}
	
	if (grep {'vdj'} @groups) {
		$id = 'VDJ_'.$id if ($id);
		$id = 'aggregation_vdj' unless ($id);
		my $aggr_csv_vdj = $dir."/aggr_vdj.csv";
		if (-e $aggr_csv_vdj) {
			my $overwrite = prompt("'aggr_vdj.csv' already exists. Overwrite the file ?  ",-yes);
			print ("If 'aggr_vdj.csv' is completed, you can run 'echo \"cd $dir ; $exec aggr --id=$id\_VDJ --csv=$aggr_csv_vdj\" | run_cluster.pl -cpu=20'") unless ($overwrite);
			print JOBS_AGGR "cd $dir ; $exec aggr --id=$id --csv=$aggr_csv_vdj" unless ($overwrite);
			die("'$aggr_csv_vdj' already exists") unless ($overwrite);
		}
		open (AGGR_CSV_VDJ, ">$aggr_csv_vdj");
		print AGGR_CSV_VDJ "sample_id,vdj_contig_info,donor,origin\n" ;
		foreach my $patient (@{$patients_name}) {
			my $group_type = lc($patient->somatic_group());
			print AGGR_CSV_VDJ $patient->name().",$dir".$patient->name()."/outs/vdj_contig_info.pb,,\n" if ($group_type eq "vdj");
		}
		close (AGGR_CSV_VDJ);
		colored::stabilo('white', "Fill the 'donor' and 'origin' columns for each vdj sample in '$aggr_csv_vdj'.");
		colored::stabilo('white', "Then run 'echo \"cd $dir ; $exec aggr --id=$id --csv=$aggr_csv_vdj\" | run_cluster.pl -cpu=20'");
		colored::stabilo('white', "Then make an archive: 'tar -czf $dir/$projectName\_aggr.tar.gz $dir/aggregation_*/outs/web_summary.html $dir/aggregation_*/outs/count/cloupe.cloupe $dir/aggregation_*/outs/count/*_bc_matrix/* $dir/aggregation_*/outs/*/vloupe.vloupe'");
	}
	close (JOBS_AGGR);
	my $cmd_tar = "tar -czf $dir/$projectName\_aggr.tar.gz $dir/aggregation_*/outs/web_summary.html $dir/aggregation_*/outs/count/cloupe.cloupe $dir/aggregation_*/outs/count/*_bc_matrix/* $dir/aggregation_*/outs/*/vloupe.vloupe";
	colored::stabilo('white', "Then, you can make an archive: '$cmd_tar'");
#	system($cmd_tar) if (-d $dir.'aggregation_exp/' and not $no_exec);
#	colored::stabilo('white', "Archive of aggragation : $dir/$projectName\_aggr.tar.gz") if (-d "$dir/$projectName\_aggr.tar.gz" and not $no_exec);
=======
>>>>>>> branch 'dev-mperin' of https://github.com/bipdimagine/polygit.git

<<<<<<< HEAD
}




if ($step eq "aggr_vdj") {
	my $type = $project->getRun->infosRun->{method};
	die ("No vdj in project $projectName") if ($type !~ /vdj/);
	my @groups = map {$_->somatic_group} @$patients_name;
	die("No 'vdj' sample to aggregate") unless (grep {'vdj'} @groups);
	my $id = $projectName.'_VDJ_aggregation';
	$id = $aggr_name if $aggr_name;
#	my $aggr_file_vdj = $dir."/jobs_aggr_vdj.txt";
#	open (JOBS_AGGR_VDJ, ">$aggr_file_vdj");
	my $aggr_csv_vdj = $dir."/aggr_vdj.csv";
	if (-e $aggr_csv_vdj) {
		my $overwrite = prompt("'aggr_vdj.csv' already exists. Overwrite the file ?  ",-yes);
		print ("If 'aggr_vdj.csv' is completed, you can run 'echo \"cd $dir ; $exec aggr --id=$id --csv=$aggr_csv_vdj\" | run_cluster.pl -cpu=20'") unless ($overwrite);
		die("'$aggr_csv_vdj' already exists") unless ($overwrite);
	}
	open (AGGR_CSV_VDJ, ">$aggr_csv_vdj");
	print AGGR_CSV_VDJ "sample_id,vdj_contig_info,donor,origin\n" ;
	foreach my $patient (@{$patients_name}) {
		my $group_type = lc($patient->somatic_group());
		if ($group_type eq "vdj") {
			print AGGR_CSV_VDJ $patient->name().",".$dir."/".$patient->name()."/outs/vdj_contig_info.pb,,\n";
		}
	}
	close (AGGR_CSV_VDJ);
#	print JOBS_AGGR_VDJ "cd $dir ; $exec aggr --id=$id\_VDJ --csv=$aggr_csv_vdj";
#	close (JOBS_AGGR_VDJ);
	colored::stabilo('white', "Fill the 'donor' and 'origin' columns for each vdj sample in '$aggr_csv_vdj'.");
	colored::stabilo('white', "Then run 'echo \"cd $dir ; $exec aggr --id=$id\_VDJ --csv=$aggr_csv_vdj\" | run_cluster.pl -cpu=20'");
	colored::stabilo('white', "Then make an archive: 'tar -czf $dir/$projectName\_aggr.tar.gz $dir/aggregation_*/outs/web_summary.html $dir/aggregation_*/outs/count/cloupe.cloupe $dir/aggregation_*/outs/count/*_bc_matrix/* $dir/aggregation_*/outs/*/vloupe.vloupe'");
}




=======
###############
# ARCHIVE
###############
>>>>>>> branch 'dev-mperin' of https://github.com/bipdimagine/polygit.git
if ($step eq "tar" or $step eq "all"){
	my $tar_cmd = "tar -cvzf $dir/$projectName.tar.gz $dir/*/outs/web_summary.html $dir/*/outs/cloupe.cloupe $dir/*/outs/vloupe.vloupe $dir/*/outs/*_bc_matrix/* ";
	warn $tar_cmd;
	die ("archive $dir/$projectName.tar.gz already exists") if -e "$dir/$projectName.tar.gz";
	system ($tar_cmd) or die("Can't tar files") unless $no_exec;
#	print "\t##########################################\n";
#	print "\tlink to send to the users : \n";
#	print "\twww.polyweb.fr/NGS/$projectName/$projectName.tar.gz \n";
#	print "\t##########################################\n";
	colored::stabilo('white', "Done") if (-e "$dir/$projectName.tar.gz");
	colored::stabilo('white', "$dir/$projectName.tar.gz") if (-e "$dir/$projectName.tar.gz");
}



<<<<<<< HEAD
=======
###############
# COPY TO /data-isilon/SingleCell/
###############
>>>>>>> branch 'dev-mperin' of https://github.com/bipdimagine/polygit.git
if ($step eq "cp" or $step eq "all"){
	my $dirout = "/data-isilon/SingleCell/$projectName/";
	my $cp_cmd = "mkdir $dirout" unless (-d $dirout);
	warn $cp_cmd unless (-d $dirout);
	system ($cp_cmd) unless ($no_exec or -d $dirout);
	foreach my $patient (@{$patients_name}) {
		my $name = $patient->name();
		next if ($name =~ /^ADT_/);
		my $cp_cmd = "cp -ru $dir/$name $dirout";
#		my $cp_cmd = "rsync -ra  $dir/$name $dirout";
		warn $cp_cmd;
		system ($cp_cmd) unless $no_exec;
#		system ($cp_cmd) or die("Can't copy files to '$dirout'") unless $no_exec;
	}
#	print "\t##########################################\n";
#	print "\tcp to $dirout \n";
#	print "\t##########################################\n";
	colored::stabilo('white', "Done") unless $no_exec;
	colored::stabilo('white', "cp to $dirout") unless $no_exec;

}



<<<<<<< HEAD
if ($step eq "cp_web_summaries" or $step eq "cp_web_summary"){
	my $dirout = "/data-isilon/SingleCell/$projectName/";
	unless (-d $dirout) {
		my $cp_cmd = "mkdir $dirout";
		warn $cp_cmd;
		system ($cp_cmd) unless ($no_exec);
	}
	foreach my $patient (@{$patients_name}) {
		my $name = $patient->name();
		next if ($name =~ /^ADT_/);
		confess ("$name web summary not found: '$dirout$name/outs/web_summary.html'") unless (-e "$dirout$name/outs/web_summary.html");
		my $cp_cmd = "mkdir $name ; " unless (-d $dirout.$name);
		$cp_cmd .= "mkdir $name/outs ; " unless (-d $dirout.$name.'/outs/');
		$cp_cmd .= "cp $dir$name/outs/web_summary.html $dirout$name/outs/web_summary.html";
		warn $cp_cmd;
		system ($cp_cmd) unless ($no_exec);
	}
#	print "\t##########################################\n";
#	print "\tcp web summaries to $dirout \n";
#	print "\t##########################################\n";
	colored::stabilo('white', "Done") unless $no_exec;
	colored::stabilo('white', "Web summaries copied to $dirout") unless $no_exec;

}



sub usage {
	print "
cellranger.pl
-------------
Obligatoires:
	project <s>			nom du projet
	step <s>			étape à réaliser: demultiplex, teleport, count, tar, aggr, aggr_vdj, cp, cp_web_summaries ou all (= demultiplex, count, aggr, tar, cp)
	feature_ref			tableau des ADT, seulement si step=count et qu'il y a des ADT
Optionels:
	patients <s>			noms de patients/échantillons, séparés par des virgules
	lane <i>			flowcell lane, défaut: lit le RunInfo.xml
	create-bam/nocreate-bam		générer ou non les bams lors du count, défaut: nocreate-bam
	aggr_name			noms de l'aggrégation, lors de step=aggr ou aggr_vdj
=======
###############
# COPY ONLY web_summary.html
###############
if ($step eq "cp_web_summaries" or $step eq "cp_web_summary"){
	my $dirout = "/data-isilon/SingleCell/$projectName/";
	unless (-d $dirout) {
		my $cp_cmd = "mkdir $dirout";
		warn $cp_cmd;
		system ($cp_cmd) unless ($no_exec);
	}
	foreach my $patient (@{$patients_name}) {
		my $name = $patient->name();
		next if ($name =~ /^ADT_/);
		confess ("$name web summary not found: '$dirout$name/outs/web_summary.html'") unless (-e "$dirout$name/outs/web_summary.html");
		my $cp_cmd = "mkdir $name ; " unless (-d $dirout.$name);
		$cp_cmd .= "mkdir $name/outs ; " unless (-d $dirout.$name.'/outs/');
		$cp_cmd .= "cp $dir$name/outs/web_summary.html $dirout$name/outs/web_summary.html";
		warn $cp_cmd;
		system ($cp_cmd) unless ($no_exec);
	}
#	print "\t##########################################\n";
#	print "\tcp web summaries to $dirout \n";
#	print "\t##########################################\n";
	colored::stabilo('white', "Done") unless $no_exec;
	colored::stabilo('white', "Web summaries copied to $dirout") unless $no_exec;

}



sub usage {
	print "
cellranger.pl
-------------
Obligatoires:
	project <s>			nom du projet
	step <s>			étape à réaliser: demultiplex, teleport, count, tar, aggr, aggr_vdj, cp, cp_web_summaries ou all (= demultiplex, count, aggr, tar, cp)
	feature_ref	<s>		tableau des ADT, seulement si step=count et qu'il y a des ADT
Optionels:
	patients <s>			noms de patients/échantillons, séparés par des virgules
	lane <i>			nombre de lanes sur la flowcell, défaut: lit le RunInfo.xml
	mismatches <i>		nombre de mismatches autorisés lors du démultiplexage, défaut: 0
	create-bam/nocreate-bam		générer ou non les bams lors du count, défaut: nocreate-bam
	aggr_name <s>		noms de l'aggrégation, lors de step=aggr ou aggr_vdj
>>>>>>> branch 'dev-mperin' of https://github.com/bipdimagine/polygit.git
	limit <i>			limite de jobs en simultanés, fonctionne pour step=count uniquement
	no_exec				ne pas exécuter les commandes
	help				affiche ce message

";
	exit(1);
}

