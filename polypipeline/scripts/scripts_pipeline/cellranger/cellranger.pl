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
	'lane=i'					=> \$lane,
	'mismatches=i'				=> \$mismatch,
	'create_bam!'				=> \$create_bam,
	'feature_ref|feature_csv=s'	=> \$feature_ref,
	'cmo_ref|cmo_csv=s'			=> \$cmo_ref,
	'aggr_name=s'				=> \$aggr_name,
	'chemistry=s'				=> \$chemistry,
	'cpu=i'						=> \$cpu,
	'no_exec'					=> \$no_exec,
	'help'						=> \$help,
) || die("Error in command line arguments\n");

usage() if $help;
usage() unless ($projectName);

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
	
	foreach my $patient (sort @{$patients}) {
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
# TELEPORT
#------------------------------
if (grep(/teleport/, @steps)) {
#	my $cmd = "teleport_patient.sh -project=$projectName -force=1";
	my $cmd = "$Bin/../../../teleport_patients.pl -project=$projectName -force=1";
	system($cmd) unless $no_exec;
}




unless ($project->isSomatic) {
	my $warn = "/!\\ Project $projectName is not in somatic mode. Activate somatic mode and check that the groups have been filled in, "
		."so that they can be taken into account in the analysis.";
	die ($warn) unless (grep{/@steps/} ('demultiplex', 'teleport', 'tar', 'cp', 'cp_web_summaries'));
	warn ($warn);
	
}


my @group;
foreach my $pat (@{$patients}) {
	my $group = $pat->somatic_group();
	push(@group,$group);
}



#------------------------------
# COUNT
#------------------------------
if (grep(/count|^all$/i, @steps)){
	
	if (grep {$_ =~ /adt/i} @group) {
		die("feature_ref csv required\n") unless ($feature_ref);
		die("'$feature_ref' not found") unless (-e $feature_ref);
		$feature_ref = abs_path($feature_ref);
	}
	if (grep {$_ =~ /cmo/i} @group) {
		die("cmo_ref csv required\n") unless ($cmo_ref);
		die("'$cmo_ref' not found") unless (-e $cmo_ref);
		$feature_ref = abs_path($feature_ref);
	}
	
	$create_bam = 'false' unless ($create_bam);
	$create_bam = 'true' if ($create_bam ne 'false');
	warn "create-bam=$create_bam";
	
	 if ($chemistry) {
 		warn "chemistry=$chemistry";
		my @chemistries = qw/auto threeprime fiveprime SC3Pv2 SC3Pv3 SC3Pv3-polyA SC3Pv4 SC3Pv4-polyA SC3Pv3HT SC3Pv3HT SC5P-PE SC5P-PE-v3 SC5P-R2 SC5P-PE-v3 SC3Pv1 ARC-v1/;
		die ("Chemistry option '$chemistry' not valid, should be one of: ". join(', ', @chemistries)) unless (grep { $_ eq $chemistry } @chemistries);
	 }
	
	
	my $fastq;
	my $tmp = $project->getAlignmentPipelineDir("cellranger_count");
	warn $tmp;
	
	# Vérifie que les fastq sont bien nommés
	foreach my $patient (@{$patients}) {
		my $pname = $patient->name;
		my $fastq_files = $patient->fastqFiles();
		my @fastq_files = map {values %$_} @$fastq_files;
		die ("Fastq file names (sample $pname) must follow the following naming convention: [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz")
			unless (scalar (grep {/$pname\_S\d*_L\d{3}_[IR][123]_001\.fastq\.gz$/} @fastq_files) == scalar @fastq_files);
	}
	
	sub full_cmd {
		my ($pat, $cmd) = @_;
		chomp $cmd;
		my $name = $pat->name;
		my $tmp_fastq = $tmp.'fastq/';
		system("mkdir $tmp_fastq") unless (-d "$tmp_fastq");
		my $seq_dir = $pat->getSequencesDirectory;
		my $cmd1 = "cp $seq_dir$name*.fastq.gz $tmp_fastq && ";
		my $adt = 'ADT_'.$name;
		$cmd1 = "cp $seq_dir$name*.fastq.gz $seq_dir/ADT_$name*.fastq.gz $tmp_fastq && " if (grep (/$adt/, @patient_names));
#		my $cmd2 = " --localcores=$cpu ";
		my $cmd2 = "&& cp -r $tmp$name/outs/* $tmp$name/_versions $tmp$name/_cmdline $dir$name/ ";
		system("mkdir $dir$name") unless (-d "$dir$name");
#		$cmd2 .= "&& rm $dir$name/possorted_bam.bam $dir$name/possorted_bam.bam.bai" if ($exec eq 'cellranger-atac' and $create_bam eq 'false');
#		$cmd2 .= "&& if [[ -f \"$tmp_fastq\*.fastq.gz\" ]]; then rm $tmp_fastq\*.fastq.gz ; fi ";
#		$cmd2 .= "&& rm $tmp_fastq\*$name*.fastq.gz";
		$cmd =~ s/$seq_dir/$tmp_fastq/;
		return $cmd1.$cmd.$cmd2."\n";
	}
	
	# EXP
	my $type_exp = 1 if map {uc($_) =~ /(EXP|NUCLEI)/ } @group;
	system("rm $dir/jobs_count.txt") if (-f "$dir/jobs_count.txt");
	if($type_exp){
		open (JOBS, ">$dir/jobs_count.txt");
		my @exp = grep { uc($_->somatic_group()) eq "EXP" || uc($_->somatic_group()) eq "NUCLEI" } @{$patients};
		warn "EXP/NUCLEI: ".join(', ',map($_->name,@exp));
		foreach my $e (@exp){
			my $name = $e->name;
#			warn $name;
			my $group = $e->somatic_group();
			$fastq = $e->getSequencesDirectory;
			my $prog =  $e->alignmentMethod();
			my $index = $project->getGenomeIndex($prog);
			my $cmd = "cd $tmp && $exec count --id=$name --sample=$name --fastqs=$fastq --create-bam=$create_bam --transcriptome=$index ";
			$cmd .= " --include-introns true" if ($type eq "nuclei" or lc($group) eq "nuclei"); # true par défaut
			$cmd .= " --chemistry $chemistry" if ($chemistry);
			$cmd .= "\n";
			$cmd = full_cmd($e, $cmd);
	#		warn $cmd;
			print JOBS $cmd;
		}
	}
	
	
	# VDJ
	my $type_vdj = 1 if map {uc($_) =~ /VDJ/ } @group;
	system("rm $dir/jobs_vdj.txt") if (-f "$dir/jobs_vdj.txt");
	if($type_vdj){
		open (JOBS_VDJ, ">$dir/jobs_vdj.txt");
		my @vdj = grep { uc($_->somatic_group()) eq "VDJ"} @{$patients};
		warn "VDJ: ".join(', ',map($_->name,@vdj));
		foreach my $v (@vdj){
			my $name = $v->name(); 
#			my $vfam = $v->family();
#			my $vgroup = uc($v->somatic_group());
			my $fastq = $v->getSequencesDirectory;
			my $prog =  $v->alignmentMethod();
			my $index = $project->getGenomeIndex($prog);
			my $index_vdj = $index."_vdj";
			my $cmd = "cd $tmp && cellranger vdj --sample=$name --id=$name --fastqs=$fastq --reference=$index_vdj  \n";
			$cmd = full_cmd($v, $cmd);
			print JOBS_VDJ $cmd;
		}
	}


	# ADT
	my $type_adt = 1 if map {uc($_) =~ /ADT/ } @group;
	if ($type_adt){
		open (JOBS_ADT, ">$dir/jobs_count.txt");
		my @exp = grep { uc($_->somatic_group()) eq "EXP"} @{$patients};
		warn "EXP+ADT: ".join(', ',map($_->name,@exp));
		foreach my $e(@exp){
			my $ename = $e->name(); 
			my $efam = $e->family();
#			my $egroup = uc($e->somatic_group());
			my $lib_file = $dir."/".$ename."_library.csv";
			open(LIB,">$lib_file") or die "impossible $lib_file";
			my @adt = grep {$_->family() eq $efam && $_->name() =~ $ename && uc($_->somatic_group()) eq "ADT"} @{$patients};
			warn ("no associated ADT library") if scalar(@adt)==0 ;
			next() if scalar(@adt)==0 ; 
			my $adt_name = $adt[0]->name();
			my $lib = "fastqs,sample,library_type\n".$tmp.'fastq/,'.",".$ename.",Gene Expression\n";
			$lib .= $tmp.'fastq/'.",".$adt_name.",Antibody Capture\n";
			print LIB $lib;
			close(LIB);
			my $prog =  $e->alignmentMethod();
			my $index = $project->getGenomeIndex($prog);
			my $cmd = "cd $tmp && cellranger count --id=$ename --feature-ref=$feature_ref --transcriptome=$index  --libraries=$lib_file --create-bam=$create_bam";
			$cmd .= " --chemistry $chemistry" if ($chemistry);
			$cmd .= "\n";
			$cmd = full_cmd($e, $cmd);
			print JOBS_ADT $cmd;
		}
	}

	
	# SPATIAL
	my $type_spatial = 1 if map {uc($_) =~ /SPATIAL/ } @group;
	system("rm $dir/jobs_spatial.txt") if (-f "$dir/jobs_spatial.txt");
 	if($type_spatial){
		open (JOBS_SPATIAL, ">$dir/jobs_spatial.txt");
		my @spatial = grep { uc($_->somatic_group()) eq "SPATIAL"} @{$patients};
		warn "SPATIAL: ".join(', ',map($_->name,@spatial));
 		my $slide_id = prompt("Visium Slide ID: ");
 		die("Visuim slide ID should start with V1, V4, V5 oh H1") unless ($slide_id =~ /^(V[145]|H1)/);
 		my %image_type = {
 			'Brightfield image generated by the CytAssist instrument (the CytAssist image)' => 'cytaimage',
 			'Brightfield microscope image' => 'image',
 			'Dark background fluorescence microscope image' => 'darkimage',
 			'Composite colored fluorescence microscope image' => 'colorizedimage',
 		};
 		my $imagetype = prompt('Choose the type of image you have: ', -m=>%image_type);
 		die ($imagetype);
		foreach my $s (@spatial){
			my $sname = $s->name(); 
			my $bc2 = $s->barcode2();
			my $area = $bc2;
#			my ($slide,$slide2,$area) = split("-",$bc2);

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
#			my $slide_final = $slide."-".$slide2;
#			my $json = $dir."/".$slide_final."-".$area."-".$sname.".json";
			# todo: rechercher les images tif contenant le nom du sample et/ou son area
			my $image = $dir."/".$area."-".$sname.".tif";
			die("Image '$image' not found") if (! -e $image && !$no_exec);
			my $sgroup = uc($s->somatic_group());
			my $fastq = $tmp.'fastq/';
			my $prog =  $s->alignmentMethod();
			my $index_spatial = $project->getGenomeIndex($prog);
			my $set = $index_spatial."/probe_sets/";
#			my $set = "/software/distrib/spaceranger/spaceranger-3.1.3/probe_sets/";
			$set .= "Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv" if ($project->getVersion() =~ /^HG/ && $slide_id =~ /^(V1)/);
			$set .= "Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv" if ($project->getVersion() =~ /^HG/ && $slide_id =~ /^(V[45]|H1)/);
			$set .= "Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv" if ($project->getVersion() =~ /^MM/ && $slide_id =~ /^(V[145])/);
			$set .= "Visium_Mouse_Transcriptome_Probe_Set_v2.0_mm10-2020-A.csv" if ($project->getVersion() =~ /^MM/ && $slide_id =~ /^H1/);
			
			my $cmd = "cd $tmp && spaceranger count --id=$sname --sample=$sname --fastqs=$fastq --transcriptome=$index_spatial --create-bam=$create_bam ";
			$cmd .= "--$imagetype=$image  ";
			$cmd .= "--area=$area --slide=$slide_id --probe-set=$set ";
#			$cmd .= " --loupe-alignment=$json ";
			$cmd = full_cmd($s, $cmd);
			print JOBS_SPATIAL $cmd;
		}
	}

	
	# ATAC
	my $type_atac = 1 if map {uc($_) =~ /ATAC/ } @group;
	system("rm $dir/jobs_atac.txt") if (-f "$dir/jobs_atac.txt");
	if($type_atac ){
		warn 'ATAC';
		open (JOBS_ATAC, ">$dir/jobs_atac.txt");
		my @atac = grep { uc($_->somatic_group()) eq "ATAC"} @{$patients};
		warn "ATAC: ".join(', ',map($_->name,@atac));
		foreach my $a (@atac){
			my $vname = $a->name(); 
			my $vfam = $a->family();
			my $vgroup = uc($a->somatic_group());
			my $fastq = $a->getSequencesDirectory;
			$exec = "cellranger-atac" if $type eq "atac";
			my $prog =  $a->alignmentMethod();
			my $index = $project->getGenomeIndex($prog);
			my $index_atac = $index."_atac" if $type eq "atac";
			my $cmd = "cd $tmp && $exec count --sample=$vname --id=$vname --fastqs=$fastq --reference=$index_atac  \n";
			$cmd = full_cmd($a, $cmd);
			print JOBS_ATAC $cmd;
		}
	}

	
	# CMO
	my $type_cmo = 1 if map {uc($_) =~ /CMO/ } @group;
	system("rm $dir/jobs_cmo.txt") if (-f "$dir/jobs_cmo.txt");
	if ($type =~ /cmo/ ){
		open (JOBS_CMO, ">$dir/jobs_cmo.txt");
#		warn $patient->somatic_group();
		my @cmo = grep { uc($_->somatic_group()) eq "CMO"} @{$patients};
		warn "CMO: ".join(', ',map($_->name,@cmo));
		foreach my $e(@cmo){
			my $ename = $e->name(); 
			my $efam = $e->family();
			my $egroup = uc($e->somatic_group());
			my $prog =  $e->alignmentMethod();
			my $index = $project->getGenomeIndex($prog);
			my $lib_file = $dir."/".$ename."_multi.csv";
			open(LIB,">$lib_file") or die "impossible $lib_file";
			print LIB "[gene-expression]\nreference,".$index."\n$cmo_ref.csv\n";
			
			print LIB "[libraries]\nfastq_id,fastqs,feature_types\n";
			print LIB "$ename,".$tmp.'fastq/'.$ename.",Gene Expression\n";
			my @cmo = grep {$_->family() eq $efam && uc($_->somatic_group()) eq "CMO"} @{$patients};
			my $cmo_name = $cmo[0]->name();
			print LIB "$cmo_name,".$tmp.'fastq/'.$ename.",Antibody Capture\n";
			print LIB "[samples]\nsample_id,cmo_ids\n";
			print LIB $ename."_B251,B251\n";
			print LIB $ename."_B252,B252\n";
			print LIB $ename."_B253,B253\n";
			close(LIB);
			my $cmd = "cd $tmp && cellranger multi --id=$ename --csv=$lib_file\n";
			$cmd = full_cmd($e, $cmd);
			print JOBS_CMO $cmd;
		}
	}

	
	# ARC
	my $type_arc = 1 if map {uc($_) =~ /ARC/ } @group;
	system("rm $dir/jobs_arc.txt") if (-f "$dir/jobs_arc.txt");
	if ($type_arc ){
		open (JOBS_ARC, ">$dir/jobs_arc.txt");
		$exec = "cellranger-arc";
		my @arc = grep { uc($_->somatic_group()) eq "ARC"} @{$patients};
		warn "ARC: ".join(', ',map($_->name,@arc));
		foreach my $e(@arc){
			my $ename = $e->name(); 
			my $efam = $e->family();
			my $egroup = uc($e->somatic_group());
			my $lib_file = $dir."/".$ename."_library.csv";
			open(LIB,">$lib_file") or die "impossible $lib_file";
			my @atac = grep {$_->family() eq $efam && uc($_->somatic_group()) eq "ARC"} @{$patients};
			next() if scalar(@atac == 0);
			my $atac_name = $atac[0]->name() ;
			my $prog =  $e->alignmentMethod();
			my $index = $project->getGenomeIndex($prog);
			my $index_arc = $index."_arc";
			my $lib = "fastqs,sample,library_type\n".$tmp.'fastq/'.",".$ename.",Gene Expression\n";
			$lib .= $tmp.'fastq/'.",".$atac_name.",Chromatin Accessibility\n";
			print LIB $lib;
			close(LIB);
			my $cmd = "cd $tmp && $exec count --id=$ename --transcriptome=$index_arc  --libraries=$lib_file\n";
			$cmd = full_cmd($e, $cmd);
			print JOBS_ARC $cmd;
		}
	}
		
	close(JOBS);
	close(JOBS_VDJ);
	close(JOBS_ADT);
	close(JOBS_SPATIAL);
	close(JOBS_ATAC);
	close(JOBS_CMO);
	close(JOBS_ARC);


	my $cmd2 = "cat $dir/jobs*.txt | run_cluster.pl -cpu=$cpu";
	warn $cmd2;
	sleep(5) unless ($no_exec);
	my $exit = system ($cmd2) unless ($no_exec);
	exit($exit) if ($exit);
	
	# Open web summaries
	my @error;
	my $web_summaries = "";
	foreach my $patient (@{$patients}) {
		next if ($patient->somatic_group =~ /^ADT$/i);
		my $file = $dir.$patient->name."/web_summary.html";
		$web_summaries .= $file.' ' if (-e $file);
		push(@error, $file) unless (-e $file or $no_exec);
	}
	my $cmd3 = "firefox ".$web_summaries;
	$cmd3 = "google-chrome ".$web_summaries if (getpwuid($<) eq 'shanein');
	warn $cmd3 if ($web_summaries);
	system($cmd3.' &') if ($web_summaries and not $no_exec);
	die("Web summaries not found: ".join(', ', @error)) if (@error and not $no_exec);

	unless ($no_exec) {
		print "\t------------------------------------------\n";
		print("Check the web summaries:\n");
		print("$dir/*/web_summary.html\n\n");
		print "\t------------------------------------------\n\n";
	}
	
}




#------------------------------
# AGGREGATION
#------------------------------
if (grep(/aggr/, @steps)){
	my @groups = map {$_->somatic_group} @$patients;
	warn Dumper \@group;
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
		foreach my $patient (@{$patients}) {
			warn 
			my $group_type = lc($patient->somatic_group());
			if ($group_type eq "exp") {	# ne "adt" or $_ ne "vdj"
				print AGGR_CSV $patient->name().",".$dir."/".$patient->name()."/molecule_info.h5\n";
			}
		}
		close AGGR_CSV;
		print JOBS_AGGR "cd $dir ; $exec aggr --id=$id --csv=$aggr_csv";
		warn ("cat $aggr_file | run_cluster.pl -cpu=$cpu");
		system ("cat $aggr_file | run_cluster.pl -cpu=$cpu") unless $no_exec;
		die("Error while running cellranger aggr") unless (-d $id);
		if (-d "$dir/$id/" and not $no_exec) {
			print("--------------------\n");
			print("$dir/$id/\n\n");
		}
	}
	
	if (grep {'vdj'} @groups) {
		$id = 'VDJ_'.$id if ($id);
		$id = 'aggregation_vdj' unless ($id);
		my $aggr_csv_vdj = $dir."/aggr_vdj.csv";
		if (-e $aggr_csv_vdj) {
			my $overwrite = prompt("'aggr_vdj.csv' already exists. Overwrite the file ?  ",-yes);
			print ("If 'aggr_vdj.csv' is completed, you can run 'echo \"cd $dir ; $exec aggr --id=$id\_VDJ --csv=$aggr_csv_vdj\" | run_cluster.pl -cpu=$cpu'") unless ($overwrite);
			print JOBS_AGGR "cd $dir ; $exec aggr --id=$id --csv=$aggr_csv_vdj" unless ($overwrite);
			die("'$aggr_csv_vdj' already exists") unless ($overwrite);
		}
		open (AGGR_CSV_VDJ, ">$aggr_csv_vdj");
		print AGGR_CSV_VDJ "sample_id,vdj_contig_info,donor,origin\n" ;
		foreach my $patient (@{$patients}) {
			my $group_type = lc($patient->somatic_group());
			print AGGR_CSV_VDJ $patient->name().",$dir".$patient->name()."/vdj_contig_info.pb,,\n" if ($group_type eq "vdj");
		}
		close (AGGR_CSV_VDJ);
		print("--------------------\n");
		print("Fill the 'donor' and 'origin' columns for each vdj sample in '$aggr_csv_vdj'.\n");
		print("Then run 'echo \"cd $dir ; $exec aggr --id=$id --csv=$aggr_csv_vdj\" | run_cluster.pl -cpu=$cpu'\n");
		print("Then make an archive: 'tar -czf $dir/$projectName\_aggr.tar.gz $dir/aggregation_*/web_summary.html $dir/aggregation_*/count/cloupe.cloupe $dir/aggregation_*/count/*_bc_matrix/* $dir/aggregation_*/*/vloupe.vloupe'\n\n");
	}
	close (JOBS_AGGR);
#	my $cmd_tar = "tar -czf $dir/$projectName\_aggr.tar.gz $dir/aggregation_*/web_summary.html $dir/aggregation_*/count/cloupe.cloupe $dir/aggregation_*/count/*_bc_matrix/* $dir/aggregation_*/*/vloupe.vloupe\n";
#	print("Then, you can make an archive: '$cmd_tar'\n");
#	system($cmd_tar) if (-d $dir.'aggregation_exp/' and not $no_exec);
#	print("Archive of aggragation : $dir/$projectName\_aggr.tar.gz\n") if (-d "$dir/$projectName\_aggr.tar.gz" and not $no_exec);

}




#------------------------------
# AGGREGATION VDJ
#------------------------------
if (grep(/aggr_vdj/, @steps)) {
#	my $type = $project->getRun->infosRun->{method};
#	die ("No vdj in project $projectName") if ($type !~ /vdj/);
	my @groups = map {$_->somatic_group} @$patients;
	die("No 'vdj' sample to aggregate") unless (grep {'vdj'} @groups);
	my $id = $projectName.'_VDJ_aggregation';
	$id = $aggr_name if $aggr_name;
#	my $aggr_file_vdj = $dir."/jobs_aggr_vdj.txt";
#	open (JOBS_AGGR_VDJ, ">$aggr_file_vdj");
	my $aggr_csv_vdj = $dir."/aggr_vdj.csv";
	if (-e $aggr_csv_vdj) {
		my $overwrite = prompt("'aggr_vdj.csv' already exists. Overwrite the file ?  ",-yes);
		print ("If 'aggr_vdj.csv' is completed, you can run 'echo \"cd $dir ; $exec aggr --id=$id --csv=$aggr_csv_vdj\" | run_cluster.pl -cpu=$cpu'") unless ($overwrite);
		die("'$aggr_csv_vdj' already exists") unless ($overwrite);
	}
	open (AGGR_CSV_VDJ, ">$aggr_csv_vdj");
	print AGGR_CSV_VDJ "sample_id,vdj_contig_info,donor,origin\n" ;
	foreach my $patient (@{$patients}) {
		my $group_type = lc($patient->somatic_group());
		if ($group_type eq "vdj") {
			print AGGR_CSV_VDJ $patient->name().",".$dir."/".$patient->name()."/vdj_contig_info.pb,,\n";
		}
	}
	close (AGGR_CSV_VDJ);
	print "--------------------\n";
	print "Fill the 'donor' and 'origin' columns for each vdj sample in '$dir$aggr_csv_vdj'.\n";
	print "Then run 'echo \"cd $dir ; $exec aggr --id=$id\_VDJ --csv=$aggr_csv_vdj\" | run_cluster.pl -cpu=$cpu'";
	print "Then make an archive: 'tar -czf $dir/$projectName\_aggr.tar.gz $dir/aggregation_*/web_summary.html $dir/aggregation_*/count/cloupe.cloupe $dir/aggregation_*/count/*_bc_matrix/* $dir/aggregation_*/*/vloupe.vloupe'";
}



#------------------------------
# CSV with infos/paths for Francesco
#------------------------------
if (grep(/info/i, @steps)){
	system("$Bin/cellranger_infos.pl -project=$projectName -out_dir=$dir");
}



#------------------------------
# COPY to SingleCell shared directory
#------------------------------
if (grep(/^cp(_web_summar(y|ies))?$|^all$/i, @steps)){
	my $cmd_cp = "$Bin/cellranger_copy.pl -project=$projectName ";
	$cmd_cp .= "-patients=$patients_name " if ($patients_name);
	$cmd_cp .= "-all_outs " if (grep(/^cp$/i, @steps));
	$cmd_cp .= "-no_exec " if ($no_exec);
	system($cmd_cp);
}



#------------------------------
# ARCHIVE / TAR
#------------------------------
if (grep(/tar|^all$/i, @steps)){
	my $cmd_tar = "$Bin/cellranger_tar.pl -project=$projectName ";
	$cmd_tar .= "-patients=$patients_name " if ($patients_name);
	$cmd_tar .= "-create_bam " if ($create_bam and $create_bam ne 'false');
	$cmd_tar .= "-no_exec " if ($no_exec);
	system($cmd_tar);
}




#------------------------------
# Upload Imagine CloudBOX ??
#------------------------------
if (grep(/cloud/, @steps)) {
	my $archive = "$dir/$projectName.tar.gz";
#	my $archive = "test.txt";
	die ("No archive found: $archive") unless (-e $archive);
	my $username = 'melodie.perin' if (getpwuid($<) eq 'mperin');
	$username = prompt('Username: ') unless $username;
	die unless $username;
	my $password = 'cZN5k*Rtt4Qc+B6M' if ($username eq 'melodie.perin');
	$password  = prompt('Password: ', -e=>'') unless ($password);
	my $remote_dir = "Archives single cell";
	$remote_dir =~ s/ /%20/g;  # encodage des espaces
	my $upload_cmd = "curl -u $username:$password -T $archive --progress-bar https://cloudbox.institutimagine.org/remote.php/dav/files/$username/$remote_dir/";
	warn $upload_cmd;
	
}



sub usage {
	print "
$0
-------------
Obligatoires:
	project <s>                nom du projet
	feature_ref	<s>            tableau des ADT, obligatoire seulement si step=count et qu'il y a des ADT
	cmo_ref	<s>                tableau des CMO, obligatoire seulement si step=count et qu'il y a des CMO
Optionels:
	steps <s>                  étape à réaliser: demultiplex, teleport, count, tar, aggr, aggr_vdj, cp, cp_web_summary ou all (= demultiplex, count, cp_web_summary, tar)
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

