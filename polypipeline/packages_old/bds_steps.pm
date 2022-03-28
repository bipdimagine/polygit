package bds_steps;
use Term::ReadKey;
use IO::Prompt;
use Moose;
use MooseX::Method::Signatures;
use Data::Printer;
use FindBin qw($Bin);
use Time::Local;
use POSIX qw(strftime);
use Term::ANSIColor;
use Data::Dumper;
use job_bds;
use sample;

extends (qw(bds_root));

###########################################
# Variables globales pour les paths, id et versions des programmes #
###########################################






has 'fastq_files'  =>(
	is        => 'rw',
	isa       => 'HashRef',
	default   => sub{ {} },
);

#########
# Méthodes  #
#########





method job_super_types (){
	 confess();
	#my (@types) = sort{$self->{jobs_super_type}->{$a} <=> $self->{jobs_super_type}->{$b}} keys %{$self->{jobs_super_type}};
	#return \@types;
}




###########################
# Méthodes liées aux étapes  du pipeline #
###########################





method start_again (Str :$filein! ,Object :$previous! ){ 
	confess();
#	my $name = $self->patient()->name();
#	my $stepname = "restart";
#	my $project = $self->patient()->getProject();
#	my $dir_out = $project->getAlignmentPipelineDir($self->align_method)."/again";
#	mkdir $dir_out;
##	mkdir ( $project->getPipelineCallingDir("mpileup")) unless -e  $project->getPipelineCallingDir("mpileup");
#	my $fileout = $dir_out . "/" . $name . ".start.bam";
#	  
#	$filein = 	$self->patient()->getBamFile();
#	die("not found :".$filein) unless $filein;
#	
#	my $cmd = " ln  -s  $filein $fileout ";
#	warn $cmd;	
#	my $bai = $filein.".bai";
#	if (-e  $bai){
#		$cmd .= "ln  -s  $filein.bai $fileout.bai ";
#	}
#	else {
#		$cmd .= $self->samtools()." index $filein ";
#	}
#	
#	#trace
#	
#	#my $cmd = " mv $filein $fileout;samtools index $fileout";
#	$self->hash_steps->{$stepname}->{id}=13 ;
#	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>[$cmd],ppn=>1,fileout=>$fileout);
#	if ($self->unforce() && -e $fileout){
#	  	$self->add_skip_steps($stepname);
#		return ($previous,$fileout);
#		}
#	$self->add_running_steps($stepname);
#	if (defined $previous) {
#		$self->prev(job=>$job,prev=>$previous);
#	}	
#	return ($job,$fileout);
}

###méthode d'aln mem utilisée pour le pipeline miseq, miseq_primer (avec trimming)



###méthode utilisée dans le pipeline illumina (sans trimming)



method run_alignment_frag(Str :$filein){

	my $name = $self->patient()->name();
	#confess("check the code man");
	my ($dirin) = $self->patient()->getSequencesDirectory();
		my $method = $self->patient()->alignmentMethod();
	my $dirout = $self->project->getAlignmentPipelineDir("$method");
	 $filein = file_util::find_fragment_file($self->project,$name,$dirin,"gz");
	#my $ppn = 16;
	my $ppn = $self->nproc ;#if $self->nocluster;
	

	my $project_name =  $self->project->name();
	my $fileout = $dirout.$name.".align.bam";
	 my $type = "align-frag";
	 my $stepname = $self->patient->name."@".$type;
	  my $bin_dev = $self->script_dir();
	my $cmd = qq{perl $bin_dev/align.pl -file1=$filein -method=$method   -mode=frag -project=$project_name -name=$name -bam=$fileout -fork=$ppn};
	
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
	
#	
#	$self->hash_steps->{$stepname}->{id}= 2;	
#	my ($job_hisat,$job_final) = $self->construct_jobs(stepname=>$stepname,cmd=>[$cmd],ppn=>$proc);
#	return ([$job_hisat],[$bam]);
}

method alignment (Str :$filein){
	my $name = $self->patient()->name();
	my $method = $self->patient()->alignmentMethod();
	my $run = $self->patient->getRun();	
		if ($run->infosRun->{method} eq "fragment"){
			return $self->run_alignment_frag(filein=>$filein);
		}
		return $self->run_alignment_pe(filein=>$filein);
}

method alignment_fsgs (Str :$filein){
	my $method = $self->patient()->alignmentMethod();
		my $name = $self->patient()->name();
		my ($dirin) = $self->patient()->getSequencesDirectory();
		#my ($dirin_complement) = "/data-isilon/data/sequences/ILLUMINA/MISEQ/LAVOISIER/Complements_FSGS/FSGS_ComplementaryPanel_Antignac_26" ;
		#warn $dirin ; 
		#warn $dirin_complement ;
		
		my $project_name =  $self->project->name();
		my $dirout = $self->project->getAlignmentPipelineDir($method);
		
	#my $ppn = 16;
	my $ppn = $self->nproc if $self->nocluster;
	my $split = "_";
	my $ext1 = $split."R1".$split;
	my $ext2 = $ext1;
	$ext2 =~ s/1/2/;

	#my $files_pe1 = file_util::my $files_pe1 = file_util::find_file_pe($self->patient,$ext1);($self->patient->name,$dirin,$ext1);
	#p $files_pe1;
	my $files_pe_complement = file_util::find_file_pe_complement_directory($self->patient,$ext1);
	warn Dumper $files_pe_complement ; 
	die() ; 
	my $files_pe1 = file_util::find_file_pe_old($self->patient,$ext1);
	warn Dumper $files_pe1 ; 
	my $count_lane = scalar(@$files_pe1);
	my $nb_bam =1;
	my $already =0;
	my @cmds;
	my @bams;
	my $files;
	my @jobs;
	foreach my $file1 (@$files_pe1){
		my $file2 = $file1;
		$file2 =~ s/$ext1/$ext2/;
		if (exists $self->fastq_files->{$file1} or exists $self->fastq_files->{$file2}){
			die("same fastq file present in two different patient : $name $file1 $file2");
		}
		 $self->fastq_files->{$file1} ++;
		 $self->fastq_files->{$file2} ++;
		die("problem $file1 $file2 $dirin :".$name) unless -e $dirin.$file2;
		die("problem $file1 $file2 $dirin :".$name) unless -e $dirin.$file1;
		 my $bam = $dirout.$name.".F$nb_bam."."bwa.bam";
		 my $f1 = $dirin.$file1;
		 my $f2 = $dirin.$file2;
		
		 my $bin_dev = $self->script_dir();
		my $cmd = qq{perl $bin_dev/align.pl -file1=$f1 -file2=$f2 -method=$method -lane=$nb_bam  -project=$project_name -name=$name -bam=$bam -fork=$ppn};
		 my $type = "align#".$nb_bam;
		 my $stepname = $self->patient->name."@".$type;
		 my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$f1,$f2],fileout=>$bam,type=>$type,dir_bds=>$self->dir_bds);
		$self->current_sample->add_job(job=>$job_bds);
		push(@jobs,$job_bds);
			$nb_bam ++;
		if ($self->unforce() && -e $bam){
		 		$job_bds->skip();
		}
		push(@bams,$bam);
	} 
	
	my ($fileout) = $self->merge_bam(filein=>\@bams);
	return ($fileout);
}


method run_alignment_pe (Str :$filein!){
		
		my $method = $self->patient()->alignmentMethod();
		my $name = $self->patient()->name();
		my ($dirin) = $self->patient()->getSequencesDirectory();
		my $project_name =  $self->project->name();
		my $dirout = $self->project->getAlignmentPipelineDir($method);
		
	#my $ppn = 10;
	my $ppn = $self->nproc;# if $self->nocluster;

#	my $split = "_";
#	my $ext1 = $split."1".$split."C";
#	my $ext2 = $ext1;
#	$ext2 =~ s/1/2/;
	#my $files_pe1 = file_util::test($self->patient);
	my $files_pe1 = file_util::find_file_pe($self->patient,"");
	#my $files_pe1 = file_util::my $files_pe1 = file_util::find_file_pe($self->patient,$ext1);($self->patient->name,$dirin,$ext1);
	#p $files_pe1;
	my $count_lane = scalar(@$files_pe1);
	my $nb_bam =1;
	my $already =0;
	my @cmds;
	my @bams;
	my $files;
	my @jobs;
	print $self->patient->name().":\n";
	foreach my $cp (@$files_pe1){
		my $file1 = $cp->{R1};
		my $file2 =  $cp->{R2};
		print "\t $nb_bam : ".$cp->{R1}." ".$cp->{R2}."\n";
		
		if (exists $self->fastq_files->{$file1} or exists $self->fastq_files->{$file2}){
			die("same fastq file present in two different patient : $name $file1 $file2");
		}
		 $self->fastq_files->{$file1} ++;
		 $self->fastq_files->{$file2} ++;
		die("problem $file1 $file2 $dirin :".$name) unless -e $dirin.$file2;
		die("problem $file1 $file2 $dirin :".$name) unless -e $dirin.$file1;
		 my $bam = $dirout.$name.".F$nb_bam."."bwa.bam";
		 my $f1 = $dirin.$file1;
		 my $f2 = $dirin.$file2;
		
		 my $bin_dev = $self->script_dir();
		my $cmd = qq{perl $bin_dev/align.pl -file1=$f1 -file2=$f2 -method=$method -lane=$nb_bam  -project=$project_name -name=$name -bam=$bam -fork=$ppn};
		 my $type = "align#".$nb_bam;
		 my $stepname = $self->patient->name."@".$type;
		 my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$f1,$f2],fileout=>$bam,type=>$type,dir_bds=>$self->dir_bds);
		$self->current_sample->add_job(job=>$job_bds);
		push(@jobs,$job_bds);
			$nb_bam ++;
		if ($self->unforce() && -e $bam){
		 		$job_bds->skip();
		}
		push(@bams,$bam);
	} 
	
	my ($fileout) = $self->merge_bam(filein=>\@bams);
	return ($fileout);
		
}




method merge_bam (ArrayRef :$filein) {
	return $self->merge_bamba(filein=>$filein);
}

method merge_bamba (ArrayRef :$filein) {
	my $name = $self->patient->name();
	my $project = $self->patient->getProject;
	my ($f)  = File::Util->new(); 
	my$m = $self->patient->alignmentMethod();
	my $fileout = $project->getAlignmentPipelineDir($m) . "/" . $name . ".align.bam";
	my $merge_files = join( " ", @$filein );
	my $bamba =$self->project->getSoftware('sambamba');

	
	my $ppn =$self->nproc;
	$ppn = int($self->nproc/2) if $self->nocluster;
	

	my $cmd = " $bamba  merge -t $ppn $fileout $merge_files  " ;

	if (scalar(@$filein) == 1) {
		if (-e $fileout){
					$cmd = "ln -s  $merge_files $fileout  " ;#if (scalar(@files) == 1);
					
		}
		else {$cmd = "ln -s  $merge_files $fileout "}
		

	}
	 my $type = "merge-bamba";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>$filein,fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
}



method read_group_illumina (Str :$filein){
	
	my $name = $self->patient()->name();
	my $fileout;
	unless ($filein){
		 $filein = $self->patient()->getBamFileName() ;
		 $filein =~s/bwa/bwa\/moitie/;
		 die() unless -e $filein;
		 my $outputdir = $self->project->getAlignmentPipelineDir("bwa");
		 $fileout=  $outputdir."/".$self->patient()->name.".rg.bam";
	}
	else {
	$fileout = $filein;
	$fileout =~ s/bam/rg\.bam/;
	}
	my $filebai = $fileout;
	$filebai =~ s/bam/bai/;

	my $run = $self->patient()->getRun();
	my $machine = $run->machine;
	my $run_name = $run->plateform_run_name();
	my $run_date = $run->date;
	my $constructor =  $run->machine_constructor();
	my $constructormachine= $constructor."-".$machine ;
	my $plateform =  $run->plateform();
	my $bar_code= $self->patient()->barcode();
	my $project = $self->patient->getProject;
	my $project_name= $project->name();
	my $patient_name = $self->patient()->name();
	
	my $picard_path = $self->project->getSoftware('picard_path');
	
	my $ppn = 2 ;
	$ppn =1 if $self->nocluster;
	my $java = $self->project->getSoftware('java');
	$java ="java" unless -e $java;
	 
	 my $picard =  $java." -jar ".$self->project->getSoftware('picard_path');
	my $cmd  = $picard." AddOrReplaceReadGroups CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT  INPUT=".$filein." OUTPUT=".$fileout." RGDS=".$project_name." RGLB=".$run_name." RGSM=".$name." RGCN=".$plateform." RGID=".$run_name." RGPL=".$constructormachine." RGPU=1 " ;
	$cmd = $cmd."  && mv $filebai $fileout.bai ";
	
	 my $type = "read-group";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	
	return ($fileout);
}

method move_bam (Str :$filein){
	my $name = $self->patient()->name();
	my$m = $self->patient->alignmentMethod();
	my $project_name = $self->patient()->getProject->name();

	
	my $method = $self->patient()->alignmentMethod();
	my $dirout = $self->patient()->getProject->getAlignmentDir($method);
	my $fileout  = $self->patient()->getBamFileName() ;#$dirout."/".$self->patient->name.".bam";
	my $ppn = 4 ;
	$ppn = 1 if $self->nocluster;

	die() if 	$fileout eq $filein;
	my $bin_dev = $self->script_dir;
	
	my $cmd =  "perl $bin_dev/move_bam.pl -bam=$filein  -project=$project_name -patient=$name -fork=$ppn && ln -s $fileout $filein  ";
	 my $type = "move-bam";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	
	return ($fileout);
}

method rnaseq_metrics (Str :$filein){
	my $name = $self->patient()->name();
	my $project=$self->patient()->project;
	my$m = $self->patient->alignmentMethod();
	my $project_name = $self->patient()->getProject->name();

	my $dir_out= $project->getCountingDir("featureCounts")."/metrics";
	system("mkdir $dir_out && chmod a+rwx $dir_out") unless -e $dir_out;
		
	my $method = $self->patient()->alignmentMethod();
	my $fileout  =$dir_out."/$name.metrics";

	my $ppn = 2 ;
	$ppn = 1 if $self->nocluster;

	die() if 	$fileout eq $filein;
	my $refFlat = $project->refFlat_file();
	my $rRNA_file = $project->rRNA_file();
	my $opt = "";
	$opt = "RIBOSOMAL_INTERVALS=$rRNA_file" if -e $rRNA_file;
	unless (-e $refFlat){
			die("can't find $refFlat");
			
	}
	
	my $java = $project->buffer->software("java");
	my $picard =  $project->buffer->software("picard");
	
	my $cmd =  "$java -jar $picard  CollectRnaSeqMetrics I=$filein O=$fileout REF_FLAT=$refFlat STRAND=FIRST_READ_TRANSCRIPTION_STRAND $opt";
	 my $type = "metrics";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	
	return ($fileout);
}

method bam_sort  (Str :$filein){
	return $self->bam_sort_bamba(filein=>$filein);
}


method bam_sort_bamba  (Str :$filein){
	my $name = $self->patient()->name();
	my $fileout = $filein;
	$fileout =~ s/bam/sort\.bam/;
	
	my $bamba=$self->project->getSoftware('sambamba');
	
	my $patient_name = $self->patient()->name();
	my $project = $self->project->name();
#	my $tmpdir = "--tmpdir=/scratch";
	

	my $ppn = $self->nproc ;
	$ppn = int($self->nproc/2) if $self->nocluster;
		my $cmd =  $bamba." sort   -o=$fileout  --nthreads=$ppn  --out=$fileout ".$filein." ";
	
	 my $type = "sort-bam";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	
	return ($fileout);
	
}


method reorder_picard  (Str :$filein,Object :$previous){
	confess();
#	my $name = $self->patient()->name();
#	my $fileout = $filein;
#	$fileout =~ s/bam/reord\.bam/;
#	my $reference = $self->version();
#	
#	#Ajouts dans le header du bam
#	my $picard_path = $self->project->getSoftware('picard');
#	my @picard= split( "/", $picard_path);
#	my $version_picard= $picard[-1];
#
#	my $cmd =  $self->picard_exec_path()." ReorderSam   CREATE_INDEX=true INPUT=".$filein." OUTPUT=".$fileout."  REFERENCE=".$reference." ";
#
#	#trace
#	$Header_lines->{'Picard_ReorderSam'}->{'type'}="PG";
#	$Header_lines->{'Picard_ReorderSam'}->{'PN'}="ReorderSam";
#	$Header_lines->{'Picard_ReorderSam'}->{'ID'}="Picard_ReorderSam";
#	$Header_lines->{'Picard_ReorderSam'}->{'VN'}="$version_picard" ;
#	$Header_lines->{'Picard_ReorderSam'}->{'CL'}="\"".$cmd."\"" ;
#	
#	my $proc ;
#		if ($self->max_cpu()){
#			$proc = $self->max_cpu ;
#		}
#		else {$proc = 4};
#	
#	my $stepname = "reord";
#	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>[$cmd],ppn=>$proc,fileout=>$fileout,filein=>[$filein]); 
#	if ($self->restart_step (filein=>$filein, fileout=>$fileout)){
#		$self->add_skip_steps($stepname);
#		return ($previous,$fileout);
#	}
#	$self->add_running_steps($stepname);
#	if (defined $previous) {
#		$self->prev(job=>$job,prev=>$previous);
#	}	
#	return ($job,$fileout);
}


method rmdup  (Str :$filein){
	return $self->rmdup_bamba(filein=>$filein);
}


method rmdup_picard  (Str :$filein){
	confess();
#	my $name = $self->patient()->name();
#	
#	my $fileout = $filein;
#	$fileout =~ s/bam/rmdup\.bam/;
#	my $metric_file = $filein.".toto.txt";
#	
#	my $picard_path = $self->project->getSoftware('picard_path');
#	my @picard= split( "/", $picard_path);
#	my $version_picard= $picard[-1];
#	
#	my $cmd = $self->picard_exec_path()." MarkDuplicates  METRICS_FILE=$metric_file CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT INPUT=".$filein." OUTPUT=".$fileout." ";
#		
#	#trace
#	$Header_lines->{'Picard_MarkDuplicates'}->{'type'}="PG";
#	$Header_lines->{'Picard_MarkDuplicates'}->{'PN'}="MarkDuplicates";
#	$Header_lines->{'Picard_MarkDuplicates'}->{'ID'}="Picard_MarkDuplicates";
#	$Header_lines->{'Picard_MarkDuplicates'}->{'VN'}="$version_picard" ;
#	$Header_lines->{'Picard_MarkDuplicates'}->{'CL'}="\"".$cmd."\"" ;
#
#	my $stepname = "rmdup";
#	$self->hash_steps->{$stepname}->{id}=14 ;
#	
#	my $proc ;
#		if ($self->max_cpu()){
#			$proc = $self->max_cpu ;
#		}
#		else {$proc = 10};
#		
#	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>[$cmd],ppn=>$proc,fileout=>$fileout,filein=>[$filein]); 
#	if ($self->restart_step (filein=>$filein, fileout=>$fileout)){
#		$self->add_skip_steps($stepname);
#		
#		return ($previous,$fileout);
#	}	
#	$self->add_running_steps($stepname);
#	if (defined $previous) {
#			$self->prev(job=>$job,prev=>$previous);
#	}	
#	return ($job,$fileout);
}

method rmdup_nudup  (Str :$filein!){
	my $name = $self->patient()->name();
	my $fileout = $filein;
	my$m = $self->patient->alignmentMethod();
	$fileout =~ s/\.bam//;
	
	my $nudup = "/software/distrib/nudup-master/nudup.py";
	my $python ="python";
	my $bamba = $self->project->getSoftware('sambamba');
	my $tmpdir ="";
	 $tmpdir = "--tmpdir=/tmp" if $self->host eq "morgan";
	my $ppn = $self->nproc;# if $self->nocluster;
	#$tmpdir = "" if $self->nocluster;
		 	my $cmd =  "$python $nudup --paired-end $filein -o $fileout ";
		 	$fileout = $fileout.".sorted.markdup.bam";
		 my $type = "rmdup_nudup";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
		return ($fileout);
}


method rmdup_bamba  (Str :$filein!){
	my $name = $self->patient()->name();
	my $fileout = $filein;
	my$m = $self->patient->alignmentMethod();
	$fileout =~ s/bam/rmdup\.bam/;
	my $bamba = $self->project->getSoftware('sambamba');
	my $tmpdir ="";
	 $tmpdir = "--tmpdir=/tmp" if $self->host eq "morgan";
	my $ppn = $self->nproc;# if $self->nocluster;
	#$tmpdir = "" if $self->nocluster;
		 	my $cmd =  $bamba." markdup   $tmpdir  --nthreads=$ppn --overflow-list-size=2000000  ".$filein." $fileout";
		 my $type = "rmdup";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
}



method SplitNCigarReads (Str :$filein!){ 
	my $name = $self->patient()->name();
	unless ($filein){
		die();
	}
	my $fileout = $filein;
	 $fileout =~ s/bam/splincigars\.bam/;
	 
	
	 my $capture = $self->patient->getCapture();
	 my $multiplex = $capture->multiplexFile();
	 die("can't find multiplex file") unless -e $multiplex;
	my $reference = $self->reference();
	my $cmd = $self->gatk()." -I  $filein -R $reference -T SplitNCigarReads -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -o $fileout ";
		my $ppn =4;
	$ppn = 1 if $self->nocluster;
		
		 my $type = "splitncigar";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
}


method mask_primer_start_end (Str :$filein!){ 
	my $name = $self->patient()->name();
	unless ($filein){
		die();
	}
	my $fileout = $filein;
	 $fileout =~ s/bam/mask\.bam/;
	 my $capture = $self->patient->getCapture();
	 my $multiplex = $capture->multiplexFile();
	 die("can't find multiplex file") unless -e $multiplex;
	 my $sam = $filein ;
	$sam=~ s/bam/mask\.sam/;
	my $bin_dev = $self->script_dir;
	
	my $ppn =4;
	$ppn = $self->nproc if $self->nocluster;
		
	my $cmd = "$bin_dev/mask_primer_start_end.pl -file=$filein -primer=$multiplex -fork=$ppn | samtools view -Sb - > $fileout  " ;

	my $type = "mask-primer";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
}

method covariate_illumina  (Str :$filein){
	my $name = $self->patient()->name();
	my $pname = $self->project->name;
	my $fileout = $filein;
	 $fileout =~ s/bam/recal\.bam/;
	my $csv = $filein.".tmp1";
	# $fileout =~ s/bam/recal1\.bam/;
	 unless ($filein){
	 	confess();
	 }
	 my $filein_bai ,
	my $ppn = $self->nproc ;#if $self->nocluster;
	
	my $real_ppn  = 8;
	$real_ppn = int($self->nproc/2) if $self->nocluster;

	
	my $csv2 =$self->patient->getRecalFile();
	my $bed = $self->patient()->getCaptureBedFile() ;
	die("not found bed file : $bed") unless -e $bed;
	my $bin_dev = $self->script_dir;
	my $cmd = "$bin_dev/recalibration.pl -filein=$filein -fileout=$fileout -fork=$real_ppn -project=$pname -patient=$name ";
	
	my $type = "recalibration";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
}

method recalibration_table  (Str :$filein){
	my $name = $self->patient()->name();
	my $pname = $self->project->name;
	my $fileout = $filein;
	my $csv = $filein.".tmp1";
	 $fileout =~ s/bam/recal1\.bam/;
	 unless ($filein){
	 	confess();
	 }
	 my $filein_bai ,
	my $ppn = $self->nproc ;#if $self->nocluster;
	
	my $real_ppn  = 8;
	$real_ppn = int($self->nproc/2) if $self->nocluster;

	
	my $csv2 =$self->patient->getRecalFile();
	my $bed = $self->patient()->getCaptureBedFile() ;
	die("not found bed file : $bed") unless -e $bed;
	my $bin_dev = $self->script_dir;
	my $cmd = "$bin_dev/recalibration_table.pl -filein=$filein -fileout=$fileout -fork=$real_ppn -project=$pname -patient=$name && test -e $csv2 && mv $filein $fileout ";
	
	my $type = "recalibration";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($filein);
}

method realign_recal  (Str :$filein,Object :$previous){
	my $name = $self->patient()->name();
	my $pname = $self->project->name;
	 $filein = $self->current_sample->current_bam;
	unless ($filein) {
		confess() unless -e $filein;
	}
	
	my $outputdir = $self->project->getAlignmentPipelineDir("bwa");
 	my $csv = $outputdir."/".$name.".interval_list";
  	my $table = $outputdir."/".$name.".table";
	my $fileout = $filein;
	$fileout =~ s/bam/realign\.bam/;
	
 	my $ppn = $self->nproc;# if $self->nocluster;
	#info_step
	my $bin_dev = $self->script_dir;
	#my $cmd1 = "$bin_dev/realign_csv.pl -filein=$filein -fileout=$fileout -fork=$ppn -project=$pname -patient=$name ";
	#my $cmd2 = "$bin_dev/recal_table.pl -filein=$filein -fileout=$fileout -fork=$ppn -project=$pname -patient=$name ";
	my $cmd = "$bin_dev/realign_recal.pl -filein=$filein -fileout=$fileout -fork=$ppn -project=$pname -patient=$name ";
	my $type = "realign-recal";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip(); 
	}
	return ($fileout);
}

method realign_gatk  (Str :$filein,Object :$previous){
	my $name = $self->patient()->name();
	my $pname = $self->project->name;
	unless ($filein) {
		confess() unless -e $filein;
	}

	my $fileout = $filein;
	$fileout =~ s/bam/realign\.bam/;

	
	#my $ppn =16;
	my $ppn = $self->nproc ;#if $self->nocluster;
	#info_step
	my $bin_dev = $self->script_dir;
	my $cmd = "$bin_dev/realign.pl -filein=$filein -fileout=$fileout -fork=$ppn -project=$pname -patient=$name ";
	my $type = "realign";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
	
}

method callable_region_panel  (Str :$filein){ 
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();	
	my $project_name =$project->name();
	my $fileout = $filein;
	my $low_calling = "" ;
	$filein = $self->patient()->getBamFileName();# unless $filein;

	my $dirout= $project->getCallingPipelineDir("callable");
	$fileout = $dirout ."/".$name.".freeze";
	my $ppn =$self->nproc;
	my $reference = $project->genomeFasta();

	$ppn = int($self->nproc/2) if $self->nocluster;
	#$ppn =5;	
	die("-".$filein) unless $filein;
	

	my $bin_dev = $self->script_dir;
	
	my $cmd = "perl $bin_dev/callable_region.pl -project=$project_name  -patient=$name -fork=$ppn  -fileout=$fileout  -filein=$filein";
	my $type = "callable-region";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
}

method calling_panel (Str :$filein, Any :$low_calling){ 
    my $methods    =  $self->patient()->getCallingMethods();

    foreach my $m (@$methods){
        next if $m eq "seqnext";
        next if $m eq "SmCounter";
        warn $m;
        #next unless $m eq "duplicate_region_calling";
        $self->calling_generic(filein=>$filein,method=>$m,low_calling=>$low_calling);
    }
    return $filein;
}


method calling_generic  (Str :$filein, Str :$method, Any :$low_calling){ 
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();	
	my $project_name =$project->name();
	 $filein = $self->patient()->getBamFileName() unless $filein =~/bam/;
	
	my $dirout= $project->getVariationsDir($method);
	my $fileout = $dirout."/$name.vcf.gz";
	my $dirout1= $project->getCallingPipelineDir("callable");
	#$filein = $dirout1 ."/".$name.".freeze";
	my $ppn =$self->nproc;
	$ppn = int($self->nproc/2) if $self->nocluster;
	
	die("-".$filein) unless $filein;
	#my $low_calling_string = $low_calling;
#	 $low_calling_string = "-low_calling=1"  if $low_calling;
	 
	my $m = $self->patient()->alignmentMethod();
	my $dir_bam = $project->getAlignmentDir($m);
	my $bin_dev = $self->script_dir;
	
	#my $cmd = "perl $bin_dev/calling_panel.pl -project=$project_name  -patient=$name -fork=$ppn  -fileout=$fileout -method=$method -filein=$filein $low_calling_string";
	my $cmd = "perl $bin_dev/calling_panel.pl -project=$project_name  -patient=$name -fork=$ppn  -fileout=$fileout -method=$method -filein=$filein ";
	my $type = "calling-".$method;
#	$type = "lc-".$type  if $low_calling;
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
}



#method calling_merge  (Str :$filein){ 
#	my $name = $self->patient()->name();
#	my $project = $self->patient()->getProject();	
#	my $project_name =$project->name();
#	my $fileout = $filein;
#	my $low_calling = "" ;
#	$filein = $self->patient()->getBamFile();# unless $filein;
#
#	mkdir ( $project->getCallingPipelineDir("unifiedgenotyper")) unless -e  $project->getCallingPipelineDir("mpileup");
#	my $dirout= $project->getCallingPipelineDir("unifiedgenotyper");
#	$fileout = $dirout . "/$name/" .$name.".final.vcf";
#	my $ppn ;
#	
#		if ($self->max_cpu()){
#			$ppn = $self->max_cpu ;
#		}
#		else {$ppn = 8 ; 
#		};
#	die("-".$filein) unless $filein;
#	
#	my $m = $self->patient()->alignmentMethod();
#	my $dir_bam = $project->getAlignmentDir($m);
#	my $bin_dev = $self->script_dir;
#	
#	my $cmd = "perl $bin_dev/calling_illumina.pl -project=$project_name  -patient=$name -fork=$ppn  -out=$fileout -calling=unifiedgenotyper,freebayes,samtools ";
#	my $type = "calling-merge";
#	 my $stepname = $self->patient->name."@".$type;
#	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
#	$self->current_sample->add_job(job=>$job_bds);
#	if ($self->unforce() && -e $fileout){
#		 		$job_bds->skip();
#	}
#	return ($fileout);
#}



method calling_gvcf  (Str :$filein! ){ 
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();	
	my $project_name =$project->name();
	my $low_calling = "" ;
	$filein = $self->patient()->getBamFileName() ;#unless $filein !~/bam/;
	my $dir_gvcf_out  = $project->getGvcfDir("haplotypecaller");
 	my $fileout = $dir_gvcf_out."/".$name.".g.vcf.gz";
 	print  "$name \n " unless -e $dir_gvcf_out."/".$name.".g.vcf.gz";
	#probleme du fichier de log non défini dans le pipeline lancé sans pbs
#	my $cmd = "" ;
	my $ppn = $self->nproc ;# if $self->nocluster;
	my $real_ppn = int($self->nproc / 2);
	 $real_ppn = $ppn  if $self->host eq "morgan";


	die("-".$filein) unless $filein;
#	die($filein. " is empty") if (-z $filein);
	my $bin_dev = $self->script_dir;
	
	my $cmd = "perl $bin_dev/calling_individual_gvcf.pl -project=$project_name  -patient=$name -fork=$real_ppn -out=$fileout -window=1_000_000  ";
	my $type = "gvcf";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
}

method callable_regions  (Str :$filein! ){ 
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();	
	my $project_name =$project->name();
	unless ($filein){
			my $dir_gvcf_out  = $project->getGvcfDir("haplotypecaller");
 		 $filein = $dir_gvcf_out."/".$name.".g.vcf.gz";
	}
	die() unless $filein;

	my $ppn = 1;
	die("-".$filein) unless $filein;
#	die($filein. " is empty") if (-z $filein);
	my $bin_dev = $self->script_dir;
	my $dir_tmp_callable = $project->getCallingPipelineDir("callable");
	my $lmdb_file =  $name.".ok.callable";
	my $fileout =  $project->getCoverageCallable()."/".$lmdb_file;
	my $cmd = "perl $bin_dev/callable_regions_patients.pl -project=$project_name  -patient=$name   ";
	my $type = "callable";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
}


method calling_merge_low  (Str :$filein!  ){ 
	confess();
#	my $name = $self->patient()->name();
#	my $project = $self->patient()->getProject();	
#	my $project_name =$project->name();
#	my $fileout = $filein;
#	my $low_calling = "" ;
#	$filein = $self->patient()->getBamFile() unless $filein;
#	
#	mkdir ( $project->getCallingPipelineDir("unifiedgenotyper")) unless -e  $project->getCallingPipelineDir("mpileup");
#	my $dirout= $project->getCallingPipelineDir("unifiedgenotyper");
#	$fileout = $dirout . "/" .$name.".final.vcf";
#	my $proc ;
#	my $log_file = $self->log_file;
#	#probleme du fichier de log non défini dans le pipeline lancé sans pbs
##	my $cmd = "" ;
#		if ($self->max_cpu()){
#			$proc = $self->max_cpu ;
##			$cmd .= "perl $bin_dev/calling_illumina.pl -project=$project_name  -patient=$name -fork=$proc  -out=$fileout ;";
#		}
#		else {$proc = 16 ; 
##				my $log_file = $self->log_file;
##				$cmd .= "perl $bin_dev/calling_illumina.pl -project=$project_name  -patient=$name -fork=$proc -log=$log_file -out=$fileout ;";
#		};
#	die("-".$filein) unless $filein;
##	die($filein. " is empty") if (-z $filein);
#	
#	my $m = $self->patient()->alignmentMethod();
#	my $dir_bam = $project->getAlignmentDir($m);
#	my $recal_file = $dir_bam."/recal/".$name.".recal.table" ;
##	die ($recal_file. " is empty") if (-z $recal_file);
##die ($recal_file. " is absent ") unless (-e $recal_file);
#	
#	my $cmd = "perl $bin_dev/calling_illumina.pl -project=$project_name  -patient=$name -fork=$proc -log=$log_file -out=$fileout -low_calling -calling=unifiedgenotyper,freebayes,samtools ";
#	my $stepname = "calling-merge-low";
#	$self->hash_steps->{$stepname}->{id}=16 ;
#	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>[$cmd],ppn=>$proc,fileout=>$fileout,filein=>[$filein]);
#	if ($self->unforce() && -e $fileout){
#		die();
#	  	$self->add_skip_steps($stepname);
#			return ($previous,$filein);
#		}
#	$self->add_running_steps($stepname);
#	if (defined $previous) {
#		$self->prev(job=>$job,prev=>$previous);
#	}	
#	return ($job,$fileout);
}


method move_vcf  (Str :$filein! ){ 
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();	
	my $project_name =$project->name();	
	#my $fileout = $filein;
	my $ppn = 1;	
	my $method = "unifiedgenotyper";
	my $dirout= $project->getCallingPipelineDir($method);
	$filein = $dirout . "/$name/" .$name.".final.vcf" unless $filein;
	my $dirin= $project->getCallingPipelineDir($method);
	my $dir_snp= $project->getVariationsDir($method);
	my $fileout = $dir_snp."/$name.vcf.gz";
	my $bin_dev = $self->script_dir;
	my $cmd = "perl $bin_dev/move_individual_vcf.pl -project=$project_name  -patient=$name -method=unifiedgenotyper ";
	my $type = "move-vcf";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
}







method picard_stats (Str :$filein!){
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();
	my $project_name = $project->name();
	my $fileout = $self->patient->getMetricsFile();

	my $bin_dev = $self->script_dir;
	my $cmd = "perl $bin_dev/picard_stats.pl -project=$project_name  -patient=$name ";
	my $type = "stats";
	 my $stepname = $self->patient->name."@".$type;

	my $ppn =$self->nproc;
	$ppn = int($self->nproc/2) if $self->nocluster;
	 my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
		return ($fileout);
}


method lmdb_depth (Str :$filein){
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();
	my $project_name = $project->name();
 	$filein = $self->patient()->getBamFileName() ;
	
	
	
	my $fileout = $self->patient()->fileNoSqlDepth;
	
	my $ppn =$self->nproc;
	$ppn = int($self->nproc/2) if $self->nocluster;
	
	
	### Récupération de la version de samtools
	
	my $bin_dev = $self->script_dir;
	
	my $cmd = qq{perl $bin_dev/coverage_genome.pl -patient=$name  -fork=$ppn  -project=$project_name };
	
	my $type = "lmdb_depth";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);

}

method reorder_picard  (Str :$filein,Object :$previous){
    my $ppn = $self->nproc;
    my $name = $self->patient()->name();
    warn $name;
    my $fileout;
    my $project=$self->patient()->project;
    my $m = $self->patient->alignmentMethod();
    my $dir_prod = $project->getAlignmentDir($m);
    unless ($filein){
#redemarrage sans fichier bam dans dossier pipeline
         my $outputdir = $self->project->getAlignmentPipelineDir("bwa");
         my $outputdirstart = $outputdir."/start_again";

          my $fileprod = $self->patient()->getBamFileName() ;
            $filein =$outputdirstart."\/".$name.".bam";
      if (-e $fileprod){
#              $filein = $fileprod;
#              $filein =~ s/bwa/bwa\/start_again/;
            $filein =$outputdirstart."\/".$name.".bam";
              system("mkdir $outputdirstart && chmod a+rwx $outputdirstart ") unless -e $outputdirstart;
              system(" chmod a+w $fileprod ");
              system("mv $fileprod $filein");
          }
          die() unless -e $filein;
         

         $filein =  $outputdirstart."/".$self->patient()->name.".bam";
         $fileout = $self->project->getAlignmentPipelineDir("bwa")."/".$self->patient()->name().".reorder.bam";
        
    }
    else {
    $fileout = $filein;
    $fileout =~ s/bam/reord\.bam/;
    
    }
    my $filebai = $fileout;
    $filebai =~ s/bam/bai/;
    

    my $reference = $project->genomeFasta();
    warn $reference;

my $picard_path = $self->project->getSoftware('picard_path');
my $picard=$self->project->getSoftware('java')." -jar ".$self->project->getSoftware('picard_path');
my $cmd =  $picard." ReorderSam   CREATE_INDEX=true INPUT=".$filein." OUTPUT=".$fileout."  REFERENCE=".$reference." ALLOW_INCOMPLETE_DICT_CONCORDANCE=true ALLOW_CONTIG_LENGTH_DISCORDANCE=true" ;
warn $cmd;


my $type = "reorder_picard";
     my $stepname = $self->patient->name."@".$type;
    my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
    $self->current_sample->add_job(job=>$job_bds);

    if ($self->unforce() && -e $fileout){
                 $job_bds->skip();
    }
    return ($fileout);

}

##methode remplacée par le script coverage.pl


method coverage_samtools (Str :$filein){
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();
	my $project_name = $project->name();
 	$filein = $self->patient()->getBamFileName() ;

	
	my $coverage_dir = $project->getRootDir() . "/align/coverage/";
	mkdir $coverage_dir unless -e $coverage_dir;
	my $fileout = $coverage_dir . "/" . $name.".cov.gz";
	my $bed = $self->patient()->getCaptureFile();
	my $ppn =$self->nproc;
	$ppn = int($self->nproc/2) if $self->nocluster;
	
	
	### Récupération de la version de samtools
	
	my $bin_dev = $self->script_dir;
	
	my $cmd = qq{perl $bin_dev/coverage.pl -patient=$name -filein=$filein -dir=$coverage_dir -bed=$bed -fork=$ppn -name=$name -project=$project_name };
	
	my $type = "coverage";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);

}




method depthofcoverage (Str :$filein){
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();
	my $project_name = $project->name();
 	$filein = $self->patient()->getBamFileName() ;

	
	my $coverage_dir = $project->getRootDir() . "/align/coverage/depth/";
	mkdir $coverage_dir unless -e $coverage_dir;
	my $fileout = $coverage_dir . "/" . $name.".depth";
	my $bed = $self->patient()->getCaptureFile();
	my $ppn =$self->nproc;
	$ppn = int($self->nproc/2) if $self->nocluster;
	#$ppn =;
	
	### Récupération de la version de samtools
	
	my $bin_dev = $self->script_dir;
	
	my $cmd = qq{perl $bin_dev/depthofcoverage.pl -patient=$name -filein=$filein -fileout=$fileout -fork=$ppn  -project=$project_name };
	
	my $type = "depthofcoverage";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);

}



####Méthodes générales de manipulation des bam

method bam_to_fastq  (Str :$filein! ){ 
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();	
	my $project_name =$project->name();	
	
	
	my $ppn = 8;	
	my $m = $self->patient()->alignmentMethod();
	my $dir_bam = $project->getAlignmentDir($m);						# où trouver les bams
	my $dir_seq= $self->patient()->getSequencesDirectory(); 		# où ecrire les fastqs
	

	my $picard=$self->project->getSoftware('java')." -jar ".$self->project->getSoftware('picard_path');
	
	$filein = $self->patient()->getBamFileName();
	my $fastq1 = $dir_seq."/".$name."_R1_L001.fastq";
	my $fastq2=  $dir_seq."/".$name."_R2_L001.fastq";
	
	my $fileout=$fastq1;
	
	my $cmd = $picard." SamToFastq I=".$filein." FASTQ=".$fastq1." SECOND_END_FASTQ=".$fastq2;

	my $type = "bam-to-fastq";
	my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	
#	my ($fileout) = $self->zip_fastq(filein=>[$fastq1,$fastq2],R=>0);
	my ($fileout1) = $self->zip_fastq(filein=>[$fastq1,$fastq2],R=>1);
	return ($fileout1);
}

method zip_fastq (ArrayRef :$filein!, Str :$R! ){ 
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();	
	my $project_name =$project->name();	
	
	
	my $ppn = 8;	 
	my $fileout=$filein->[$R].".gz";

		 
	my $cmd="bgzip ".$filein->[$R]; 
	
	my $type = "zip-fastq#".$R;
	my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>$filein,fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
}
	 
1;
