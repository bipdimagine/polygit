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

  
my $bin_cecile=qq{$Bin/scripts/scripts_db_polypipeline};
my $bin_script_pipeline = qq{$Bin/scripts/scripts_pipeline};


my $projectName;
my $patients_name;
my $step;
my $feature_ref;
my $no_exec;
my $aggr_name;
my $cmo_ref;
my $lane;


my $limit;
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'step=s'=> \$step,
#	'bcl=s' => \$bcl_dir,
#	'run=s' => \$run,
	'feature_ref=s' =>  \$feature_ref,
	'no_exec=s' => \$no_exec,
	'aggr_name=s' => \$aggr_name,
	'nb_lane=s' => \$lane,
	#'low_calling=s' => \$low_calling,
);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
#$patients_name = "all" unless $patients_name;
$patients_name = $project->get_only_list_patients($patients_name);


my $run = $project->getRun();
my $run_name = $run->plateform_run_name;
my $bcl_dir = $run->bcl_dir;
warn $bcl_dir;

if ($step eq "all" || $step eq "demultiplex"){
	#die ("-bcl, -run and -nb_lane options are required") unless $bcl_dir || $run || $lane;

unless ( $lane ){
	die(" regarde le run info vieille nouille !!!!! $bcl_dir");
}
}
#my $cmd = "cat  /data-isilon/raw-data//ILLUMINA/10X/IMAGINE/230808_A00680_0406_BHLC3KDMXY/RunInfo.xml  \| grep LaneCount \| awk -F \'\"\' \'\{print \$2\}";

#$lane = `$cmd`;
#warn $step;



my $sampleSheet = $bcl_dir."/sampleSheet.csv";
my $exec= "cellranger";

open (SAMPLESHEET,">$sampleSheet");
my $dir = $project->getProjectRootPath();
print SAMPLESHEET "Lane,Sample,Index\n";

my $fastq;
my %hSamples;
my %test;
my $type;
my @group;


foreach my $patient (@{$patients_name}) {
	my $name=$patient->name();
	warn $name;
	my $bc = $patient->barcode();
	my $bc2 = $patient->barcode2();
	my $group = $patient->somatic_group();
	push(@group,$group);
	for (my $i=1;  $i<=$lane;$i++){
		print SAMPLESHEET $i.",".$name.",".$bc."\n";
	}
}
close(SAMPLESHEET);


#my $prog =  $patient->alignmentMethod();
my $tmp = $project->getAlignmentPipelineDir("cellranger_count");
warn $tmp;



open (JOBS, ">$dir/jobs_count.txt");
foreach my $patient (@{$patients_name}) {
	my $name=$patient->name();
	my $bc = $patient->barcode();
	my $run = $patient->getRun();
	my $group = "EXP";
	$group = $patient->somatic_group() if ($patient->somatic_group());
	next() unless uc($group) eq "EXP" ;
	#to do : à regarder dans profile
	$type = $run->infosRun->{method};
	next() if $type eq "arc";
	my $prog =  $patient->alignmentMethod();
	my $index = $project->getGenomeIndex($prog);
	$fastq = $patient->getSequencesDirectory();
	$exec = "cellranger-atac" if $type eq "atac";
	$exec = "cellranger-arc" if $type eq "arc";
	$index = $project->getGenomeIndex($prog)."_atac" if $type eq "atac";
	$index = $project->getGenomeIndex($prog)."_arc" if $type eq "arc";
	my $cmd = "cd $dir; $exec count --id=$name --sample=$name --fastqs=$fastq --create-bam false --transcriptome=$index ";
	$cmd .= " --include-introns " if $type eq "nuclei";
	$cmd .= "\n";
	warn $cmd;
	print JOBS $cmd;
}	

close(JOBS);
foreach my $patient (@{$patients_name}) {
	my $run = $patient->getRun();
	my $type_vdj = 1 if map {uc($_) =~ /VDJ/ } @group;
	  
	#$run->infosRun->{method};
	if($type_vdj){
		open (JOBS_VDJ, ">$dir/jobs_vdj.txt");
		my @vdj = grep { uc($_->somatic_group()) eq "VDJ"} @{$patients_name};
		foreach my $v (@vdj){
			my $vname = $v->name(); 
			my $vfam = $v->family();
			my $vgroup = uc($v->somatic_group());
			my $fastq = $v->getSequencesDirectory();
			my $prog =  $v->alignmentMethod();
			my $index = $project->getGenomeIndex($prog);
			my $index_vdj = $index."_vdj";
		#my $vdj_file = $dir."/".$ename."_librarytest.csv";
		#my $prog =  $e->alignmentMethod();
			print JOBS_VDJ "cd $dir; cellranger vdj --sample=$vname --id=$vname --fastqs=$fastq --reference=$index_vdj  \n"
		}
	}
	close(JOBS_VDJ);
}


foreach my $patient (@{$patients_name}) {
	my $run = $patient->getRun();
	my $type = $run->infosRun->{method};
	
	my $type_spatial = 1 if map {uc($_) =~ /SPATIAL/ } @group;
 	if($type_spatial  && $step eq "count" ){
		open (JOBS_SPATIAL, ">$dir/jobs_spatial.txt");
		my @spatial = grep { uc($_->somatic_group()) eq "SPATIAL"} @{$patients_name};
		foreach my $s (@spatial){
			my $sname = $s->name(); 
			my $bc2 = $s->barcode2();
			my ($slide,$slide2,$area) = split("-",$bc2);
			warn $sname;
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
			print JOBS_SPATIAL "cd $tmp;spaceranger count --sample=$sname --image=$image --id=$sname --fastqs=$fastq --transcriptome=$index_spatial --area=$area --slide=$slide_final --loupe-alignment=$json  \n";
		}
	}
	close(JOBS_SPATIAL);
}



foreach my $patient (@{$patients_name}) {
	my $run = $patient->getRun();
	my $type = $run->infosRun->{method};
	my $type_atac = 1 if map {uc($_) =~ /ATAC/ } @group;
	if($type_atac ){
		open (JOBS_ATAC, ">$dir/jobs_atac.txt");
		my @atac = grep { uc($_->somatic_group()) eq "ATAC"} @{$patients_name};
		foreach my $v (@atac){
			my $vname = $v->name(); 
			my $vfam = $v->family();
			my $vgroup = uc($v->somatic_group());
			my $fastq = $v->getSequencesDirectory();
			my $prog =  $v->alignmentMethod();
			#my $index = $project->getGenomeIndex($prog);
			$exec = "cellranger-atac" if $type eq "atac";
			my $index = $project->getGenomeIndex($prog)."_atac" if $type eq "atac";
			print JOBS_ATAC "cd $tmp; $exec count --sample=$vname --id=$vname --fastqs=$fastq --reference=$index  \n"
		}
	}
	close(JOBS_ATAC);
}



foreach my $patient (@{$patients_name}) {
	my $run = $patient->getRun();
	my $type = $run->infosRun->{method};
	my $type_adt = 1 if map {uc($_) =~ /ADT/ } @group;
	if ($type_adt){
		open (JOBS_ADT, ">$dir/jobs_count.txt");
		warn $patient->somatic_group();
		my @exp = grep { uc($_->somatic_group()) eq "EXP"} @{$patients_name};
		foreach my $e(@exp){
			my $ename = $e->name(); 
			warn $ename;
			my $efam = $e->family();
			my $egroup = uc($e->somatic_group());
			my $lib_file = $dir."/".$ename."_library.csv";
			open(LIB,">$lib_file") or die "impossible $lib_file";
			my @adt = grep {$_->family() eq $efam && uc($_->somatic_group()) eq "ADT"} @{$patients_name};
			warn ("no associated ADT library") if scalar(@adt)==0 ;
			next() if scalar(@adt)==0 ; 
			my $adt_name = $adt[0]->name();
			my $prog =  $e->alignmentMethod();
			my $index = $project->getGenomeIndex($prog);
			my $index_vdj = $index."_vdj";
			my $lib = "fastqs,sample,library_type\n".$e->getSequencesDirectory().",".$ename.",Gene Expression\n";
			$lib .= $adt[0]->getSequencesDirectory().",".$adt_name.",Antibody Capture\n";
			print LIB $lib;
			close(LIB);
			print JOBS_ADT "cd $tmp ; cellranger count --id=$ename --feature-ref=$feature_ref --transcriptome=$index  --libraries=$lib_file\n"
		}
	}
	close(JOBS_ADT);
}


foreach my $patient (@{$patients_name}) {
	my $run = $patient->getRun();
	my $type = $run->infosRun->{method};
	my $type_cmo = 1 if map {uc($_) =~ /CMO/ } @group;
	if ($type =~ /cmo/ ){
		open (JOBS_CMO, ">$dir/jobs_cmo.txt");
		warn $patient->somatic_group();
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
			print JOBS_CMO "cd $tmp ; cellranger multi --id=$ename --csv=$lib_file\n";
		}
	}
	close(JOBS_CMO);
}




foreach my $patient (@{$patients_name}) {
	my $run = $patient->getRun();
	my $type = $run->infosRun->{method};
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
			my $index = $project->getGenomeIndex($prog)."_arc";
			my $lib = "fastqs,sample,library_type\n".$e->getSequencesDirectory().",".$ename.",Gene Expression\n";
			$lib .= $atac[0]->getSequencesDirectory().",".$atac_name.",Chromatin Accessibility\n";
			print LIB $lib;
			close(LIB);
			print JOBS_ARC "cd $tmp ; $exec count --id=$ename --transcriptome=$index  --libraries=$lib_file\n"
		}
	}
}
close(JOBS_ARC);


#my $cmd = "cd $tmp; /software/bin/demultiplex.pl -dir=$bcl_dir -run=$run -hiseq=10X -sample_sheet=$sampleSheet -cellranger_type=$exec";
warn $exec;
my $cmd = "cd $tmp; $Bin/../../demultiplex/demultiplex.pl -dir=$bcl_dir -run=$run_name -hiseq=10X -sample_sheet=$sampleSheet -cellranger_type=$exec -mismatch=1";
if ($step eq "demultiplex" or $step eq "all"){
	warn $cmd;
	system $cmd or die "impossible $cmd" unless $no_exec==1;
}
#warn "cat $dir/jobs*.txt";
#die();

my $cmd2 = "cat $dir/jobs*.txt | run_cluster.pl -cpu=20";
if ($step eq "count" or $step eq "all"){
	warn $cmd2;
	system $cmd2  unless $no_exec==1;
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
	system ($aggr_cmd)  unless $no_exec==1;
}


#foreach my $patient (@{$patients_name}) {
#	my $name = $patient->name();
#	my $file = $dir."/".$name."/outs/web_summary.html";
#	print $file."\n" if -e $file ;
#}

##commande pour copier sur /data-isilon/singleCell

#my $RefDateActuelle = date();
#my $date = $RefDateActuelle->{date};

if ($step eq "tar" or $step eq "all"){
	my $tar_cmd = "tar -cvzf $dir/$projectName.tar.gz $dir/*/outs/web_summary.html $dir/*/outs/cloupe.cloupe $dir/*/outs/vloupe.vloupe $dir/*/outs/*_bc_matrix/* ";
	warn $tar_cmd;
	die ("archive $dir/$run.tar.gz already exists") if -e $dir."/".$run.".tar.gz";
	system ($tar_cmd)  unless $no_exec==1;
	#or die "impossible $tar_cmd";
	print "\t#########################################\n";
	print "\t  link to send to the users : \n";
	print "\t www.polyweb.fr/NGS/$projectName/$projectName.tar.gz \n";
	print "\t#########################################\n";
}

if ($step eq "cp" or $step eq "all"){
	foreach my $patient (@{$patients_name}) {
		my $name = $patient->name();
		my $dirout = "/data-isilon/SingleCell/$projectName";
		my $cp_cmd = "mkdir $dirout ; cp -r $dir/$name $dirout/. ; cp -r $dir/$name $dirout/. ";
		#die ("archive $dir/$run.tar.gz already exists") if -e $dir."/".$run.".tar.gz";
		system ($cp_cmd) ;# unless $no_exec==1;
		print "\t#########################################\n";
		print "\t  cp to /data-isilon/SingleCell/$projectName  \n";
		print "\t#########################################\n";
	}
}
