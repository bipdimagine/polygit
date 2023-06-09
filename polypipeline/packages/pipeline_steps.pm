package pipeline_steps;

use Moo;

use Data::Printer;
use FindBin qw($Bin);
use Time::Local;
use POSIX qw(strftime);
use Term::ANSIColor;
use Data::Dumper;


###########################################
# Variables globales pour les paths, id et versions des programmes #
###########################################


my $lifescope_dir = "/data-ibrix/lifescope/projects/liferoot/";
my $recal_ext= "recal.table";

my $bin_dev = qq{$Bin/scripts/scripts_pipeline};
my $bin_cecile=qq{$Bin/scripts/scripts_db_polypipeline};

my $Header_lines ;


has 'project' => (
	is =>'rw',
	isa =>"Object",
);

has 'dbh' => (
	is =>'rw',
	isa =>"Object",
);

has 'log_file' => (
	is	=> 'ro',
	isa	=> 'Str',
); 

has 'patient' => (
	is =>'rw',
	isa =>"Object",
);

has 'max_cpu' => (
	is	=> 'rw',
	isa	=> 'Int',
); 

has 'release' => (
	is	=> 'rw',
	default => sub {
		return "HG19";
	}
); 

has 'patient_argument' => (
	is	=> 'rw',
); 


sub  goal {

		my $self = shift;
	my $type = "bwa";
	my $dirout = $self->patient()->getProject->getAlignmentDir($type);
	my $fileout1 = $dirout."/".$self->patient->name.".bam";
	return $fileout1;
}

has 'javac' => (
	is => 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		if (defined $self->max_cpu) {return $self->project->getSoftware('java').' -Xmx20g -XX:ParallelGCThreads='.$self->max_cpu;}
		else {return $self->project->getSoftware('java').' -Xmx20g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=/scratch';}
		
	}
);



has 'gatk' => (
	is => 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		return $self->javac.' -jar '.$self->project->getSoftware('gatk');
	}
);

has 'samtools' => (
	is => 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		return $self->project->getSoftware('samtools');
	}
);


has 'bwa' => (
	is => 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		return $self->project->getSoftware('bwa');
	}
);
 
has 'tmap' => (
	is => 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		return $self->project->getSoftware('tmap');
	}
);

has 'seqtk' => (
	is => 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		return $self->project->getSoftware('seqtk');
	}
);

has 'cutadapt' => (
	is => 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		return $self->project->getSoftware('cutadapt');
	}
);

has 'metricsHeader' => (
	is => 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		return $self->project->getMetricsHeader();
	}
);

has 'picard_exec_path' => (
	is => 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		return $self->javac()." -jar ".$self->project->getSoftware('picard_path');
	}
);

has 'build' => (
	is => 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		return $self->project->getVersion();
	}
);

has 'public_data' => (
	is => 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;

		return $self->project->buffer->config->{public_data}->{HG38}.'/' if $self->release eq "HG38"; 
		return $self->project->buffer->config->{public_data}->{$self->build()}.'/';
	}
);

has 'reference'  => (
	is        => 'ro',
	isa       => 'Str',
	lazy =>1,
	default   => sub {
		my $self = shift;
		return $self->public_data()."/genome/fasta/all.fa" if $self->release eq "HG38";
		return $self->public_data()."/genome/fasta/all.fa";
	},
); 

has 'hisat_reference'  => (
	is        => 'ro',
	isa       => 'Str',
	lazy =>1,
	default   => sub {
		my $self = shift;
		return $self->public_data()."genome/hisat2/genome";
	},
);

has 'cellranger_reference'  => (
	is        => 'ro',
	isa       => 'Str',
	lazy =>1,
	default   => sub {
		my $self = shift;
		return $self->public_data()."genome/10X/latest";
	},
);

has 'bwa_reference'  => (
	is        => 'ro',
	isa       => 'Str',
	lazy =>1,
	default   => sub {
		my $self = shift;
		return $self->public_data()."genome/bwa-latest/all.fa";
	},
);
has 'unforce' => (
	is        => 'rw',
	isa       => 'Int',
	default   => 1,
); 

has 'running_steps' =>(
	is        => 'rw',
	isa       => 'HashRef',
	default   => sub{ {} },
);

has 'skip_steps' =>(
	is        => 'rw',
	isa       => 'HashRef',
	default   => sub{ {} },
);

has 'fastq_files'  =>(
	is        => 'rw',
	isa       => 'HashRef',
	default   => sub{ {} },
);
my $default_file = {
	rmdup => ".1.bam"
};

### table de hash pour les steps id
has 'hash_steps' => (
	is => 'rw',
	isa       => 'HashRef',
	default   => sub{ {} },
	);
	
### table de hash pour les cmd
has 'hash_cmds' => (
	is => 'rw',
	isa       => 'HashRef',
	default   => sub{ {} },
	);
	

#########
# Méthodes  #
#########



method align_method {
	my $methods = $self->patient()->alignmentMethods();
 	warn ("more than one alignment method use default mapread") if  	scalar(@$methods) >1;
  	return "mapreads"  if  	scalar(@$methods) >1;
	return $methods->[0];
}

method start_job {
	my $date = `date`;
	chomp($date);
	my $log_file = $self->log_file;
	system ("echo ===================== $date > $log_file");
	my $job = $self->LoggingJobs(stepname=>$self->patient->name,text=>"START :".$date."\n");	
	return $job;
}

method end_job {
	my $date = `date`;
	chomp($date);
	my $job = $self->LoggingJobs(stepname=>"END_".$self->patient->name,text=>"\nEND PATIENT   :  ".$self->patient->name." ".$date."\n");
	return $job;
}
 
method final_job {
	my $date = `date`;
	chomp($date);
	my $name = $self->project->name();
	my $patient_arg = $self->patient_argument();
	my $log_file = $self->log_file;
	my $job    = PBS::Client::Job->new(
		# Job declaration options
		# Resources options
		ppn    => 1,                           # process per node
		nodes  => 1,
		script => "final_job",
		cmd    => ["perl $bin_dev/final_jobs.pl -project=$name -patients=$patient_arg -log=$log_file"],                # command to be submitted
	);
#	
#	#my $job = $self->LoggingJobs(stepname=>"COMPLETE " ,text=>"\n\x1B[32m*_*_*_*_*-*-*-*-*_*_*_*_\n ALL JOBS COMPLETED \n*_*_*_*_*-*-*-*-*_*_*_*_\n".$date."\n*_*_*_*_*-*-*-*-*_*_*_*_\n\x1B[0m");
#	return $job;
}
 
method init_log_file {
	my $dir = $self->patient->getProject()->getPipelineDir();
	my $file = $dir."/".$self->patient->name().".log";
	unlink($file) if -e $file;
	return $file;
}

method get_step_id(){
	$self->{step_id}->{$self->patient->name} = 1 unless exists  $self->{step_id}->{$self->patient->name};
	return $self->{step_id}->{$self->patient->name};
}

method LoggingJobs (Str :$stepname!,Str :$text!,Str :$step_id){ 
	my $log_file = $self->log_file;

	my $patient_name=$self->patient->name() ;
	my $analysis_id =  $self->{$patient_name}->{analysis_id};
	#confess() unless $self->hash_steps->{$stepname}->{id};
	####insérer dans cmd de soumission du jobOK la commande pour renseigner le changer le statut de l'étape (complete ou error)	
	 my $cmd =[]; 
	
	my $cmd_db="";
	my $color = "\033[0;34m";
	if ($text eq "ERROR" ){
		$color  = "\033[1;31m";
	}
	elsif ($text eq "OK" ){
			$color ="\033[1;32m";
	}
	elsif ($text =~/END/ ){
			$color ="\033[1;35m";
	}
	
	if ($step_id){
		if ($text eq "ERROR" ){
		
			push(@$cmd,qq{perl $bin_cecile/change_step_status_end_in_db.pl -analysis_id=$analysis_id -step_id=$step_id -status=4 -patient=$patient_name;},qq{perl $bin_cecile/change_analysis_status_end_in_db.pl -analysis_id=$analysis_id -status=3;});
		}
		if ($stepname eq "COMPLETE "){
			push(@$cmd, qq{perl $bin_cecile/change_analysis_status_end_in_db.pl -analysis_id=$analysis_id -status=1;});
		}
	}
	
	
	my $reset_color = "\033[0m";
	push(@$cmd,qq{echo -e "$color" `date` "$stepname\t $text $reset_color">> $log_file ;});
	$stepname.= "_OK" if $text eq "OK";
	$stepname.= "_ERROR" if $text eq "ERROR";
	my $jobsOK    = PBS::Client::Job->new(
		# Job declaration options
		# Resources options
		ppn    => 1,                           # process per node
		nodes  => 1,
		script => $stepname,
		cmd    =>[$cmd],                # command to be submitted
	);
	return $jobsOK;
}

method add_running_steps($name){
	push(@{$self->running_steps->{$self->patient->name}},$name);
}

method add_skip_steps($name){
	push(@{$self->skip_steps->{$self->patient->name}},$name);
}

method print_steps($name) {
	$name = $self->patient->name() unless defined $name;
	print colored ['black ON_BRIGHT_YELLOW'], $name.": \t";
	print colored ["green"], join(",",@{$self->running_steps->{$name}});
	print colored ["blue"], " skip=> ";
	print colored ["red"], join(",",@{$self->skip_steps->{$name}}) if exists $self->skip_steps->{$name} ;
	print color 'reset';
	print "\n";
}

method print_all_steps() {
	foreach my $name (keys %{$self->running_steps}){
		$self->print_steps($name);
	} 
}

method construct_jobs_tmp (Str :$cmd!,Str :$stepname!, Int :$ppn!) {
	my $log_file = $self->log_file;
	my $start = qq{echo  "\t\t** start =>  $stepname :" `hostname` ": $ppn  : " `date` >> $log_file ; };
	my $patient_name=$self->patient->name() ;
	
	my $toto = $cmd;
	#$toto =~ s/&&/ PUIS/g;
	#my $db_start =qq{perl $bin_cecile/change_step_status_start_in_db.pl -analysis_id=$analysis_id -step_id=$step_id -status=2 -patient=$patient_name -cmd="$toto" &&};
	#my $db_end = qq{perl $bin_cecile/change_step_status_end_in_db.pl -analysis_id=$analysis_id -step_id=$step_id -status=1 -patient=$patient_name;};
	$cmd=$cmd;
	$self->hash_cmds->{$patient_name}->{$stepname} = $cmd ;

	my $jobs = PBS::Client::Job->new(
         # Job declaration options
         # Resources options
       	 script => $stepname,
         ppn    => $ppn,  # process per node
         nodes  => 1,
		 cmd    => [$cmd],# command to be submitted 
 	);

 	my $jobsOK = $self->LoggingJobs(stepname=>$stepname,text=>"OK");
 	my $jobsERROR = $self->LoggingJobs(stepname=>$stepname,text=>"ERROR");
 	$jobs->next( { ok => [$jobsOK], fail => [$jobsERROR] } );
 	return ($jobs,$jobsOK); 	
}

method construct_jobs (Str :$cmd!,Str :$stepname!, Int :$ppn!) {
	my $log_file = $self->log_file;
	my $start = qq{echo  "\t\t** start =>  $stepname :" `hostname` ": $ppn  : " `date` >> $log_file ; };
	my $patient_name=$self->patient->name() ;
	my $analysis_id =  $self->{$patient_name}->{analysis_id};
	my $step_id= $self->get_step_id();
	my $toto = $cmd;
	$toto =~ s/&&/ PUIS/g;
	
	
	my $db_start =qq{perl $bin_cecile/change_step_status_start_in_db.pl -analysis_id=$analysis_id -step_id=$step_id -status=2 -patient=$patient_name -cmd="$toto" };
	
	#my $db_end = qq{ perl $bin_cecile/change_step_status_end_in_db.pl -analysis_id=$analysis_id -step_id=$step_id -status=1 -patient=$patient_name;};
	$cmd=$cmd ;
	$self->hash_cmds->{$patient_name}->{$stepname} = $cmd ;

	my $jobs = PBS::Client::Job->new(
         # Job declaration options
         # Resources options
       	 script => $stepname,
         ppn    => $ppn,  # process per node
         nodes  => 1,
		 cmd    => [[$start,$db_start,$cmd]],# command to be submitted 
 	);

 	my $jobsOK = $self->LoggingJobs(stepname=>$stepname,text=>"OK",step_id=>$step_id);
 	my $jobsERROR = $self->LoggingJobs(stepname=>$stepname,text=>"ERROR",step_id=>$step_id);
 	$jobs->next( { ok => [$jobsOK], fail => [$jobsERROR] } );
 	return ($jobs,$jobsOK); 	
}


method prev (Object :$job! , Object :$prev!){
	 push((@{$job->{prev}->{ok}},$prev));
}
method prev_end (Object :$job! , Object :$prev!){
	 push((@{$job->{prev}->{end}},$prev));
}
method complete (ArrayRef :$previous!,Bool :$not_run) {
	my $job = $self->final_job();
	if (defined $previous) {
		foreach my $pre (@$previous){
			$self->prev_end(job=>$job,prev=>$pre);
		}
	}
	return $job;
}

method create_filename(Str :$file!,Str :$ext!){
	$file =~ s/bam/$ext\.bam/;
	return ($file,-e $file);
}

method restart_step (Str :$filein, Str:$fileout){
	return ($self->unforce() && -e $fileout && -s $fileout >=( (-s $filein) * 0.3) );
}



###########################
# Méthodes liées aux étapes  du pipeline #
###########################





method start_again (Str :$filein! ,Object :$previous! ){ 
	my $name = $self->patient()->name();
	my $stepname = $name."_restart";
	my $project = $self->patient()->getProject();
	my $dir_out = $project->getAlignmentPipelineDir($self->align_method)."/again";
	mkdir $dir_out;
#	mkdir ( $project->getPipelineCallingDir("mpileup")) unless -e  $project->getPipelineCallingDir("mpileup");
	my $fileout = $dir_out . "/" . $name . ".start.bam";
	  
	$filein = 	$self->patient()->getBamFile();
	die("not found :".$filein) unless $filein;
	
	my $cmd = " ln  -s  $filein $fileout && ";
	warn $cmd;	
	my $bai = $filein.".bai";
	if (-e  $bai){
		$cmd .= "ln  -s  $filein.bai $fileout.bai ";
	}
	else {
		$cmd .= $self->samtools()." index $filein ";
	}
	
	#trace
	$Header_lines->{'Start_again'}->{'type'}="PG";
	$Header_lines->{'Start_again'}->{'PN'}="Start_again";
	$Header_lines->{'Start_again'}->{'ID'}="Start_again";
	$Header_lines->{'Start_again'}->{'VN'}="Start_again" ;
	$Header_lines->{'Start_again'}->{'CL'}="\"".$cmd."\"" ;
	
	#my $cmd = " mv $filein $fileout;samtools index $fileout";
	$self->hash_steps->{$stepname}->{id}=13 ;
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>1);
	if ($self->unforce() && -e $fileout){
	  	$self->add_skip_steps($stepname);
			return ($previous,$fileout);
		}
	$self->add_running_steps($stepname);
	if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	return ($job,$fileout);
}

###méthode d'aln mem utilisée pour le pipeline miseq, miseq_primer (avec trimming)

method bwa_pe_clean(Str :$filein,Object :$previous!){
	my $name = $self->patient()->name();
	my $split = "_";
	my $ext1 = $split."R1".$split;
	my $ext2 = $ext1;
	$ext2 =~ s/1/2/;
	
	my $project = $self->patient()->getProject();	
	my $dirout = $project->getAlignmentPipelineDir("bwa");
	my ($dirin) = $self->patient()->getSequencesDirectory();

	my $files_pe1 = file_util::find_patient_casava_pe_1($self->patient->name,$dirin,$ext1);
	my $count_lane = scalar(@$files_pe1);
	my $nb_bam = 1;
	my @bams;
	my @sampe_jobs;
	my $already="";
	my @cmds;
	my $fullname = $dirin.$name;
	
	my $proc ;
	if ($self->max_cpu()){
		$proc = $self->max_cpu ;
	}
	else {$proc = 16};

	my $sammem = $self->bwa()." mem -t $proc ";
	
	##Récupération de la version de bwa
	my $bwa_path = $self->project->getSoftware('bwa') ;
	
	my $cmd ="$bwa_path 2>&1 |grep Version  " ;
	my $version_bwa = `$cmd`;
	$version_bwa =~ s/Version:\ //g;
	$version_bwa="bwa-".$version_bwa;
	$version_bwa =~ s/\n//g;
	
	### Récupération de la version de samtools
	my $samtools_path = $project->getSoftware('samtools') ;
	my $cmd2 ="$samtools_path 2>&1 |grep Version  " ;
	my $version_samtools = `$cmd2`;
	$version_samtools =~ s/Version:\ //g;
	$version_samtools="samtools-".$version_samtools ;
	$version_samtools =~ s/\n//g;
	
	my $reference_bwa = $self->project->gatk_reference_bwa_file();

	foreach my $file1 (@$files_pe1){
		my $file2 = $file1;
		$file2 =~ s/$ext1/$ext2/;
		die("problem $file1 $file2 $dirin :".$name." ; R2 file doesn't exist") unless -e $dirin.$file2;
		my ($num_lane) = grep{/L0/} split($split,$file1);
		unless ($num_lane){
			$num_lane = "L1";#.$nb_bam;
		}
		warn $file1."-".$file2;
		

		my $project = $self->patient()->getProject();	
		my $project_name = $project->name;	
		
		##Récupération de la version de la date du run
		my $run = $self->patient()->getRun();
		my $run_date = $run->date;
		
		my $reference_bwa = $self->project->gatk_reference_bwa_file();
		$file1 = $dirin.$file1;
		##Récupération de l'heure et date de début d'analyse = date d'un fichier fastq
		my $datestring = strftime("%Y:%m:%d-%H:%M:%S", localtime( (stat($file1))[9] ));
		
		$file2 = $dirin.$file2;
		my $file3 = $file1.".fastq";
		my $file4 = $file2.".fastq";
		my $file5 = $file2.".se.fastq";
		my $bam = $dirout.$name.".F$nb_bam."."bwa.bam";
		my $sam= $dirout.$name.".F$nb_bam.bwa.sam" ;
#		my $cmd_seqtk = "";
		my $cmd_sampe_pourRG = "$sammem $reference_bwa $file1 $file2 -R readgroup -M >$sam && ".$self->samtools()." view -b -S $sam >$bam && rm $sam  ";
		my $readgroup = qq{\@RG\\tID:$name\_$num_lane\\tSM:$name\\tPL:ILLUMINA\\tDS:$name\_$num_lane\\tPG:bwa mem};

		#avec trimming
		my $cmd_seqtk = $self->seqtk()." trimfq $file1 -q 0.01  > $file3 && ".$self->seqtk()." trimfq -q 0.01  $file2 > $file4 ";
		my $cmd_sampe = $cmd_seqtk."$sammem $reference_bwa $file3 $file4 -R '$readgroup' -M >$sam &&".$self->samtools()." view -b -S $sam >$bam &&rm $sam && rm $file3 && rm $file4 ";
		#sans trimming
#		my $cmd_sampe = "$sammem $reference_bwa $file1 $file2 -R '$readgroup' -M | ".$self->samtools()." view -b -S - >$bam ; ";
		
		$Header_lines->{'Start_analysis'}->{'type'}="PG";
		$Header_lines->{'Start_analysis'}->{'ID'}="Start_analysis";
		$Header_lines->{'Start_analysis'}->{'VN'}="Fastq_date";
		$Header_lines->{'Start_analysis'}->{'PN'}=$datestring ;
		
		$Header_lines->{'Run_date'}->{'type'}="PG";
		$Header_lines->{'Run_date'}->{'ID'}="Run_date";
		$Header_lines->{'Run_date'}->{'PN'}=$run_date ;
		
		$Header_lines->{'Samtools'}->{'type'}="PG";
		$Header_lines->{'Samtools'}->{'ID'}="Samtools";
		$Header_lines->{'Samtools'}->{'VN'}=$version_samtools ;
		$Header_lines->{'Samtools'}->{'CL'}="samtools index" ;
		
		$Header_lines->{'Bwa'}->{'type'}="PG";
		$Header_lines->{'Bwa'}->{'PN'}="bwa mem";
		$Header_lines->{'Bwa'}->{'ID'}="Bwa";
		$Header_lines->{'Bwa'}->{'VN'}="$version_bwa" ;
		$Header_lines->{'Bwa'}->{'CL'}="\"".$cmd_sampe_pourRG."\"" ;
		
		$already ++ if -e $bam;
		push(@cmds,$cmd_sampe);
		$nb_bam++;
		push(@bams,$bam);
	} 
	
	if ($already eq $count_lane && $self->unforce() ){
		$self->add_skip_steps("align mem");
		return ([$previous],\@bams);
		die();
	}
	else {
		my $nb_bam =1;
		foreach my $cmd_sampe (@cmds){
			$self->hash_steps->{$name."_".$nb_bam."_trim_align_bwa_mem"}->{id}=1 ;
			my ($job_sampe,$job_final) = $self->construct_jobs(stepname=>$name."_".$nb_bam."_trim_align_bwa_mem",cmd=>$cmd_sampe,ppn=>$proc);
			$self->add_running_steps($name."_".$nb_bam."_trim_align_bwa_mem");
			push(@sampe_jobs,$job_sampe);
			$nb_bam ++;
		} 
	}
	return (\@sampe_jobs,\@bams);
}

###méthode utilisée dans le pipeline illumina (sans trimming)

sub construct_view_and_sort {
	my ($self,$file_out,$proc) = @_;
	my $samtools = $self->samtools;
	my $ftmp = $file_out.".tmp.bam";
	my $sambamba =$self->project->getSoftware('sambamba');
	my $cmd = "$samtools view -Sb - > $ftmp && $sambamba sort $ftmp -o $file_out -t $proc && rm  $ftmp ";
}




method run_alignment_frag(Str :$filein,Object :$previous!,Str :$method){

	my $name = $self->patient()->name();
	my ($dirin) = $self->patient()->getSequencesDirectory();
	my $dirout = $self->project->getAlignmentPipelineDir("$method");
	my $file = file_util::find_fragment_file($name,$dirin,"gz");
	my $proc ;
	if ($self->max_cpu()){
		$proc = $self->max_cpu ;
	}
	else {$proc = 16};

	my $project_name =  $self->project->name();
	my $bam = $dirout.$name.".align.bam";
	my $cmd = qq{perl $bin_dev/align.pl -file1=$file -method=$method   -project=$project_name -name=$name -bam=$bam -fork=$proc;test -f $bam };
	my $stepname = $name."_"."frg_align";
	$self->hash_steps->{$stepname}->{id}= 2;	
	my ($job_hisat,$job_final) = $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$proc);
	return ([$job_hisat],[$bam]);
}

method alignment (Str :$filein,Object :$previous!){
	my $name = $self->patient()->name();
		my $method = $self->patient()->alignmentMethod();
			my $run = $self->patient()->getRun();
			my $filebam_out = $self->project->getAlignmentPipelineDir($method) . "/" . $name . ".align.bam";
			
			if ($self->unforce() && -e $filebam_out){
				$self->add_skip_steps("alignment");
				return ($previous,$filebam_out);
			}	
		if ($run->infosRun->{method} eq "fragment"){
			return $self->run_alignment_frag(filein=>$filein,previous=>$previous,method=>$method);
		}
		return $self->run_alignment_pe(filein=>$filein,previous=>$previous,method=>$method);
}

method run_alignment_pe (Str :$filein,Object :$previous!,Str :$method){
		my $name = $self->patient()->name();
		my ($dirin) = $self->patient()->getSequencesDirectory();
		my $project_name =  $self->project->name();
		my $dirout = $self->project->getAlignmentPipelineDir($method);
		
		my $proc ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 16};
		
		
		
			my $split = "_";
	my $ext1 = $split."R1".$split;
	my $ext2 = $ext1;
	$ext2 =~ s/1/2/;
	
	my $files_pe1 = file_util::find_patient_casava_pe_1($self->patient->name,$dirin,$ext1);
	p $files_pe1;
	my $count_lane = scalar(@$files_pe1);
	my $nb_bam =1;
	my $already =0;
	my @cmds;
	my @bams;
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
		my $cmd = qq{perl $bin_dev/align.pl -file1=$f1 -file2=$f2 -method=$method -lane=$nb_bam  -project=$project_name -name=$name -bam=$bam -fork=$proc;test -f $bam };
		$already ++ if -e $bam;
		push(@cmds,$cmd);
		$nb_bam++;
		push(@bams,$bam);
	} 
my @sampe_jobs;
	if ($already eq $count_lane && $self->unforce() ){
		$self->add_skip_steps("align mem");
		return ([$previous],\@bams);
		die();
	}
	else {
		my $nb_bam =1;
		foreach my $cmd_sampe (@cmds){
			my $step_name = $name."_".$nb_bam."_align";
			
			$self->hash_steps->{$step_name}->{id}= 2;
			my ($job_sampe,$job_final) = $self->construct_jobs(stepname=>$step_name,cmd=>$cmd_sampe,ppn=>$proc);
			$self->add_running_steps($name."_".$nb_bam."_align");
			push(@sampe_jobs,$job_sampe);
			$nb_bam ++;
		} 
		
	}
	my ($job,$fileout) = $self->merge_bam(filein=>\@bams,previous=>\@sampe_jobs);
	return ($job,$fileout);
		
}




method bwa_pe_hg38(Str :$filein,Object :$previous!){
	my $name = $self->patient()->name();
	my $split = "_";
	my $ext1 = $split."R1".$split;
	my $ext2 = $ext1;
	$ext2 =~ s/1/2/;
	
	my $project = $self->patient()->getProject();	
	my $project_name = $project->name;
	my $dirout = $project->getAlignmentPipelineDir("bwa");
	my ($dirin) = $self->patient()->getSequencesDirectory();

	my $files_pe1 = file_util::find_patient_casava_pe_1($self->patient->name,$dirin,$ext1);
	
	my $fullname = $dirin.$name;
	
	my $proc ;
	if ($self->max_cpu()){
		$proc = $self->max_cpu ;
	}
	else {$proc = 24};
	
	my $sammem = $self->bwa()." mem -t $proc ";
	
	my $count_lane = scalar(@$files_pe1);
	my $nb_bam = 1;
	my @bams;
	my @sampe_jobs;
	my $already="";
	my @cmds;
	
	foreach my $file1 (@$files_pe1){
		my $file2 = $file1;
		$file2 =~ s/$ext1/$ext2/;
		if (exists $self->fastq_files->{$file1} or exists $self->fastq_files->{$file2}){
			die("same fastq file present in two different patient : $name $file1 $file2");
		}
		 $self->fastq_files->{$file1} ++;
		 $self->fastq_files->{$file2} ++;
		die("problem $file1 $file2 $dirin :".$name) unless -e $dirin.$file2;
		my ($num_lane) = grep{/L0/} split($split,$file1);
		unless ($num_lane){
			$num_lane = "L1";#.$nb_bam;
		}
		
		##Récupération de la version de bwa
		my $bwa_path = $self->project->getSoftware('bwa') ;
		my $cmd ="$bwa_path 2>&1 |grep Version  " ;
		my $version_bwa = `$cmd`;
		$version_bwa =~ s/Version:\ //g;
		$version_bwa="bwa-".$version_bwa;
		$version_bwa =~ s/\n//g;
		
		### Récupération de la version de samtools
		my $samtools_path = $project->getSoftware('samtools') ;
		my $cmd2 ="$samtools_path 2>&1 |grep Version  " ;
		my $version_samtools = `$cmd2`;
		$version_samtools =~ s/Version:\ //g;
		$version_samtools="samtools-".$version_samtools ;
		$version_samtools =~ s/\n//g;
		
		my $run = $self->patient()->getRun();
		my $run_date = $run->date;
		
		my $reference_bwa = $self->project->gatk_reference_bwa_file();
		$reference_bwa = $self->bwa_reference_hg38;
		$file1 = $dirin.$file1;
		##Récupération de l'heure et date de début d'analyse = date d'un fichier fastq
		my $datestring = strftime("%Y:%m:%d-%H:%M:%S", localtime( (stat($file1))[9] ));
		
		$file2 = $dirin.$file2;
		my $file3 = $file1.".fastq";
		my $file4 = $file2.".fastq";
		my $file5 = $file2.".se.fastq";
		my $bam = $dirout.$name.".F$nb_bam.bwa.bam";
		my $bamtmp = $dirout.$name.".F$nb_bam.bwa.tmp.bam";
		my $sam= $dirout.$name.".F$nb_bam.bwa.sam" ;
		my $cmd_sort =  $self->picard_exec_path()." SortSam   SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT   INPUT=".$bamtmp." OUTPUT=".$bam." ";
		
		
		my $cmd_sampe_pourRG = "$sammem $reference_bwa $file1 $file2 -R readgroup -M > $sam && ".$cmd_sort." &&  rm $sam ";
		my $bamsort = $dirout.$name.".F$nb_bam.bwa";
		my $readgroup = qq{\@RG\\tID:$name\_$num_lane\\tSM:$name\\tPL:ILLUMINA\\tDS:$name\_$num_lane\\tPG:bwa mem};

		my $cmd_sampe = "$sammem $reference_bwa $file1 $file2 -R '$readgroup' -M | ".$self->samtools()." view -b -S - >$bamtmp && ".$cmd_sort."&& rm $bamtmp";

		$Header_lines->{'Start_analysis'}->{'type'}="PG";
		$Header_lines->{'Start_analysis'}->{'ID'}="Start_analysis";
		$Header_lines->{'Start_analysis'}->{'VN'}="Fastq_date";
		$Header_lines->{'Start_analysis'}->{'PN'}=$datestring ;
		
		$Header_lines->{'Run_date'}->{'type'}="PG";
		$Header_lines->{'Run_date'}->{'ID'}="Run_date";
		$Header_lines->{'Run_date'}->{'PN'}=$run_date ;
		
		$Header_lines->{'Samtools'}->{'type'}="PG";
		$Header_lines->{'Samtools'}->{'ID'}="Samtools";
		$Header_lines->{'Samtools'}->{'VN'}=$version_samtools ;
		$Header_lines->{'Samtools'}->{'CL'}="samtools index" ;
		
		$Header_lines->{'Bwa'}->{'type'}="PG";
		$Header_lines->{'Bwa'}->{'PN'}="bwa mem";
		$Header_lines->{'Bwa'}->{'ID'}="Bwa";
		$Header_lines->{'Bwa'}->{'VN'}="$version_bwa" ;
		$Header_lines->{'Bwa'}->{'CL'}="\"".$cmd_sampe_pourRG."\"" ;
		$already ++ if -e $bam;
		push(@cmds,$cmd_sampe);
		$nb_bam++;
		push(@bams,$bam);
	} 

	if ($already eq $count_lane && $self->unforce() ){
		$self->add_skip_steps("align mem");
		return ([$previous],\@bams);
		die();
	}
	else {
		my $nb_bam =1;
		foreach my $cmd_sampe (@cmds){
			$self->hash_steps->{$name."_".$nb_bam."_align_bwa_mem"}->{id}= 2;
			my ($job_sampe,$job_final) = $self->construct_jobs(stepname=>$name."_".$nb_bam."_align_bwa_mem",cmd=>$cmd_sampe,ppn=>$proc);
			$self->add_running_steps($name."_".$nb_bam."_align_bwa_mem");
			push(@sampe_jobs,$job_sampe);
			$nb_bam ++;
		} 
		
	}
	return (\@sampe_jobs,\@bams);
}


method test  (Str :$filein,Object :$previous!) {
	my $name = $self->patient->name();
	

	
	my $proc ;
	if ($self->max_cpu()){
		$proc = $self->max_cpu ;
	}
	else {$proc = 16};

	my $cmd = "  " ;
	my $stepname ="test";
	$cmd = "echo $name 1  " ;#if (scalar(@files) == 1);
	my $cmd2 = "echo $name 2  YYYY" ;	
	my $cmd3 = "echo $name 3  XXXXXXX" ;	
	$self->hash_steps->{$stepname}->{id}=3 ;
	my $job = PBS::Client::Job->new(
         # Job declaration options
         # Resources options
       	 script => $stepname,
         ppn    => 1,  # process per node
         nodes  => 1,
		 cmd    => [[$cmd,$cmd2,$cmd3]],# command to be submitted 
 	);
 	
 	
	if (defined $previous) {
			$self->prev(job=>$job,prev=>$previous);
	}

	$self->add_running_steps($stepname);
	return ($job,$filein);
}

method merge_bam (ArrayRef :$filein,ArrayRef :$previous!,Bool :$not_run) {
	return $self->merge_bamba(filein=>$filein,previous=>$previous);
}

method merge_bamba (ArrayRef :$filein,ArrayRef :$previous!,Bool :$not_run) {
	my $name = $self->patient->name();
	my $project = $self->patient->getProject;
	my ($f)  = File::Util->new();
	my$m = $self->patient->alignmentMethod();
	my $fileout = $project->getAlignmentPipelineDir($m) . "/" . $name . ".align.bam";
	my $stepname = $name."_merge_bamba";

	# my $ttt equivalent a restart_step
	my $ttt = ($self->unforce() && -e $fileout && -s $fileout  > 100 );

	if ($ttt){
		warn("SKIPPING MERGE! $fileout");
		$self->add_skip_steps($stepname);
		return ($previous->[0],$fileout);
	}
	


	my $merge_files = join( " ", @$filein );
	
	
#	my $final_bam = $name . ".1.bam"; ##ne sert à rien

	my $outputdir = $project->getAlignmentPipelineDir("$m");
	my $picard_path = $self->project->getSoftware('picard_path');
	my @picard= split( "/", $picard_path);
	my $version_picard= $picard[-1];
	my $bamba =$self->project->getSoftware('sambamba');

	
	my $proc ;
	if ($self->max_cpu()){
		$proc = $self->max_cpu ;
	}
	else {$proc = 16};

	my $cmd = " $bamba  merge -t $proc $fileout $merge_files  " ;
	$Header_lines->{'Picard_Multilane_merge picard'}->{'type'}="PG";
	$Header_lines->{'Picard_Multilane_merge picard'}->{'PN'}="bamba_merge picard";
	$Header_lines->{'Picard_Multilane_merge picard'}->{'ID'}="BAMBA merge ";
	$Header_lines->{'Picard_Multilane_merge picard'}->{'VN'}="$version_picard" ;
	$Header_lines->{'Picard_Multilane_merge picard'}->{'CL'}="\"".$cmd."\"" ;
	
	if (scalar(@$filein) == 1) {
		if (-e $fileout){
					$cmd = "rm $fileout && ln -s  $merge_files $fileout  " ;#if (scalar(@files) == 1);
		}
		else {$cmd = "ln -s  $merge_files $fileout "}
		

	}
	$self->hash_steps->{$stepname}->{id}=3 ;
	my ($job,$job_next) = $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$proc);
	if (defined $previous) {
		foreach my $pre (@$previous){
			$self->prev(job=>$job,prev=>$pre);
		}
	}

	$self->add_running_steps($stepname);
	return ($job,$fileout);
}



method read_group_illumina (Str :$filein,Object :$previous){
	my $name = $self->patient()->name();
	my $fileout = $filein;
	$fileout =~ s/bam/rg\.bam/;
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
	my @picard= split( "/", $picard_path);
	my $version_picard= $picard[-1];
	my $stepname = $name."_rg"; 
	
	if ($self->restart_step (filein=>$filein, fileout=>$fileout)){
	
		$self->add_skip_steps($stepname);
		return ($previous,$fileout);
	}
	 
	
	my $cmd  = $self->picard_exec_path()." AddOrReplaceReadGroups CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT INPUT=".$filein." OUTPUT=".$fileout." RGDS=".$project_name." RGLB=".$run_name." RGSM=".$name." RGCN=".$plateform." RGID=".$run_name." RGPL=".$constructormachine." RGPU=1 " ;
	
	$cmd = $cmd." && mv $filebai $fileout.bai ";
	
	#trace
	$Header_lines->{'Picard_AddOrReplaceReadGroups'}->{'type'}="PG";
	$Header_lines->{'Picard_AddOrReplaceReadGroups'}->{'PN'}="AddOrReplaceReadGroups";
	$Header_lines->{'Picard_AddOrReplaceReadGroups'}->{'ID'}="Picard_AddOrReplaceReadGroups";
	$Header_lines->{'Picard_AddOrReplaceReadGroups'}->{'VN'}="$version_picard" ;
	$Header_lines->{'Picard_AddOrReplaceReadGroups'}->{'CL'}="\"".$cmd."\"" ;
	
	
	$self->hash_steps->{$stepname}->{id}=6 ;	
	
	my $proc ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 4};
	my ($job,$job_next) = $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$proc);
	if ($self->restart_step (filein=>$filein, fileout=>$fileout)){
		$self->add_skip_steps($stepname);
		return ($previous,$fileout);
	}
	
	$self->add_running_steps($stepname);
	if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}
		
	return ($job,$fileout);
}

method move_bam (Str :$filein,Object :$previous){
	my $name = $self->patient()->name();
	my$m = $self->patient->alignmentMethod();
	my $type = $m;
	my $project_name = $self->patient()->getProject->name();
	#$type = "bwa" if $filein =~ /bwa/;
	#$type = "tmap" if $filein =~ /tmap/;
	my $dirout = $self->patient()->getProject->getAlignmentDir($type);
	my $fileout1 = $dirout."/".$name.".bam";
	my $stepname = $name."_move_bam"; 
		if ($self->restart_step (filein=>$filein, fileout=>$self->goal)){
	
		$self->add_skip_steps($stepname);
		return ($previous,$self->goal);
	}
	my $proc ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 4};
		
	my $cmd =  "perl $bin_dev/move_bam.pl -bam=$filein  -project=$project_name -patient=$name -fork=$proc && ln -s $fileout1 $filein && chmod a-w $filein ";
	
	$self->hash_steps->{$stepname}->{id}= 10;	
	my ($job,$job_next) = $self->construct_jobs_tmp(stepname=>$stepname,cmd=>$cmd,ppn=>1);
	
	$self->add_running_steps($stepname);
	if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}
	return ($job,$fileout1);
	
	##Récupération de l'heure et date de fin d'analyse
	my $datestring = strftime("%Y:%m:%d-%H:%M:%S", localtime( (stat($fileout1))[9] ));
	#trace
	$Header_lines->{'move_bam'}->{'type'}="PG";
	$Header_lines->{'move_bam'}->{'PN'}="move_bam";
	$Header_lines->{'move_bam'}->{'ID'}="move_bam";
	$Header_lines->{'move_bam'}->{'VN'}="move_bam.pl" ;
	$Header_lines->{'move_bam'}->{'CL'}="\"".$cmd."\"" ;
	
	$Header_lines->{'End_analysis'}->{'type'}="PG";
	$Header_lines->{'End_analysis'}->{'PN'}="$datestring";
	$Header_lines->{'End_analysis'}->{'ID'}="End_analysis";
}

method bam_sort  (Str :$filein,Object :$previous){
	return $self->bam_sort_bamba(filein=>$filein,previous=>$previous);
}


method bam_sort_bamba  (Str :$filein,Object :$previous){
	my $name = $self->patient()->name();
	my $fileout = $filein;
	$fileout =~ s/bam/sort\.bam/;
	
	my $bamba=$self->project->getSoftware('sambamba');
	
	my $patient_name = $self->patient()->name();
	my $project = $self->project->name();
	my $tmpdir = "--tmpdir=/scratch";

	#trace
	$Header_lines->{'SamBamba_Sort'}->{'type'}="PG";
	$Header_lines->{'SamBamba_Sort'}->{'PN'}="sambamba";
	$Header_lines->{'SamBamba_Sort'}->{'ID'}="sambamba";
	

	my $proc ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 16};
		my $cmd =  $bamba." sort $tmpdir   -o=$fileout  --nthreads=$proc  --out=$fileout ".$filein." ";
		$Header_lines->{'SamBamba_Sort'}->{'CL'}="\"".$cmd."\"" ;
	my $stepname = $name."_bamba_sort";
	$self->hash_steps->{$stepname}->{id}=21 ;
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$proc); 
	if ($self->restart_step (filein=>$filein, fileout=>$fileout)){
		$self->add_skip_steps($stepname);
		return ($previous,$fileout);
	}
	
	$self->add_running_steps($stepname);
	if (defined $previous) {
			$self->prev(job=>$job,prev=>$previous);
	}	
	return ($job,$fileout);
}


method reorder_picard  (Str :$filein,Object :$previous){
	my $name = $self->patient()->name();
	my $fileout = $filein;
	$fileout =~ s/bam/reord\.bam/;
	my $reference = $self->version();
	
	#Ajouts dans le header du bam
	my $picard_path = $self->project->getSoftware('picard');
	my @picard= split( "/", $picard_path);
	my $version_picard= $picard[-1];

	my $cmd =  $self->picard_exec_path()." ReorderSam   CREATE_INDEX=true INPUT=".$filein." OUTPUT=".$fileout."  REFERENCE=".$reference." ";

	#trace
	$Header_lines->{'Picard_ReorderSam'}->{'type'}="PG";
	$Header_lines->{'Picard_ReorderSam'}->{'PN'}="ReorderSam";
	$Header_lines->{'Picard_ReorderSam'}->{'ID'}="Picard_ReorderSam";
	$Header_lines->{'Picard_ReorderSam'}->{'VN'}="$version_picard" ;
	$Header_lines->{'Picard_ReorderSam'}->{'CL'}="\"".$cmd."\"" ;
	
	my $proc ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 4};
	
	my $stepname = $name."_reord";
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$proc); 
	if ($self->restart_step (filein=>$filein, fileout=>$fileout)){
		$self->add_skip_steps($stepname);
		return ($previous,$fileout);
	}
	$self->add_running_steps($stepname);
	if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	return ($job,$fileout);
}

method rmdup  (Str :$filein,Object :$previous){
	return $self->rmdup_bamba(filein=>$filein,previous=>$previous);
}


method rmdup_picard  (Str :$filein,Object :$previous){
	my $name = $self->patient()->name();
	
	my $fileout = $filein;
	$fileout =~ s/bam/rmdup\.bam/;
	my $metric_file = $filein.".toto.txt";
	
	my $picard_path = $self->project->getSoftware('picard_path');
	my @picard= split( "/", $picard_path);
	my $version_picard= $picard[-1];
	
	my $cmd = $self->picard_exec_path()." MarkDuplicates  METRICS_FILE=$metric_file CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT INPUT=".$filein." OUTPUT=".$fileout." ";
		
	#trace
	$Header_lines->{'Picard_MarkDuplicates'}->{'type'}="PG";
	$Header_lines->{'Picard_MarkDuplicates'}->{'PN'}="MarkDuplicates";
	$Header_lines->{'Picard_MarkDuplicates'}->{'ID'}="Picard_MarkDuplicates";
	$Header_lines->{'Picard_MarkDuplicates'}->{'VN'}="$version_picard" ;
	$Header_lines->{'Picard_MarkDuplicates'}->{'CL'}="\"".$cmd."\"" ;

	my $stepname = $name."_rmdup";
	$self->hash_steps->{$stepname}->{id}=14 ;
	
	my $proc ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 10};
		
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$proc); 
	if ($self->restart_step (filein=>$filein, fileout=>$fileout)){
		$self->add_skip_steps($stepname);
		
		return ($previous,$fileout);
	}	
	$self->add_running_steps($stepname);
	if (defined $previous) {
			$self->prev(job=>$job,prev=>$previous);
	}	
	return ($job,$fileout);
}

method rmdup_bamba  (Str :$filein,Object :$previous){
	my $name = $self->patient()->name();
	my $fileout = $filein;
	my$m = $self->patient->alignmentMethod();
	$fileout =~ s/bam/rmdup\.bam/;
	my $bamba = $self->project->getSoftware('sambamba');
	
	my $cmd_res = $self->samtools()." view -H $filein && grep \"^\@CO\" $filein";

	my $stepname = $name."_rmdup_bamba";
	
		#trace
	$Header_lines->{'SamBamba_Rmdup'}->{'type'}="PG";
	$Header_lines->{'SamBamba_Rmdup'}->{'PN'}="sambamba";
	$Header_lines->{'SamBamba_Rmdup'}->{'ID'}="sambamba";
	
	
	$self->hash_steps->{$stepname}->{id}=22 ;
	my $proc ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 16};
	my $tmpdir = "--tmpdir=/scratch";
	my $cmd =  $bamba." markdup $tmpdir    --nthreads=$proc   ".$filein." $fileout";

	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$proc); 
	if ($self->restart_step (filein=>$filein, fileout=>$fileout)){
		$self->add_skip_steps($stepname);
		
		return ($previous,$fileout);
	}	
	$self->add_running_steps($stepname);
	if (defined $previous) {
			$self->prev(job=>$job,prev=>$previous);
	}	
	return ($job,$fileout);
}




method mask_primer (Str :$filein!,Any :$previous){ 
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
	my $cmd = "$bin_dev/mask_primer.pl -file=$filein -primer=$multiplex > $sam && samtools view -Sb $sam > $fileout && rm $sam " ;

	#trace
	$Header_lines->{'Mask_primer'}->{'type'}="PG";
	$Header_lines->{'Mask_primer'}->{'PN'}="mask_primer";
	$Header_lines->{'Mask_primer'}->{'ID'}="mask_primer";
	$Header_lines->{'Mask_primer'}->{'VN'}="mask_primer" ;
	$Header_lines->{'Mask_primer'}->{'CL'}="\"".$cmd."\"" ;
	
	my $stepname = $name."_mask_primer";
	$self->hash_steps->{$stepname}->{id}= 15;
	my $proc ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 4};
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$proc);
	 
	 $self->add_running_steps($stepname);
	  if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	 return ($job,$fileout);	 
}

method SplitNCigarReads (Str :$filein!,Any :$previous){ 
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
	my $proc ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 4};
		my $stepname = $name."_splitncigar";
	my ($job,$job_next) =  $self->construct_jobs_tmp(stepname=>$stepname,cmd=>$cmd,ppn=>$proc);
	  if ($self->restart_step (filein=>$filein, fileout=>$fileout)){
		$self->add_skip_steps($stepname);
		
		return ($previous,$fileout);
	}	
	 $self->add_running_steps($stepname);
	  if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	 return ($job,$fileout);	 
}


method mask_primer_start_end (Str :$filein!,Any :$previous){ 
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
	my $cmd = "$bin_dev/mask_primer_start_end.pl -file=$filein -primer=$multiplex > $sam && samtools view -Sb $sam > $fileout && rm $sam " ;

	#trace
	$Header_lines->{'Mask_primer'}->{'type'}="PG";
	$Header_lines->{'Mask_primer'}->{'PN'}="mask_primer";
	$Header_lines->{'Mask_primer'}->{'ID'}="mask_primer_start_end";
	$Header_lines->{'Mask_primer'}->{'VN'}="mask_primer" ;
	$Header_lines->{'Mask_primer'}->{'CL'}="\"".$cmd."\"" ;
	
	my $stepname = $name."_mask_primer_se";
	$self->hash_steps->{$stepname}->{id}= 15;
	my $proc ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 4};
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$proc);
	 
	 $self->add_running_steps($stepname);
	  if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	 return ($job,$fileout);	 
}

method covariate_illumina  (Str :$filein,Object :$previous){
	my $name = $self->patient()->name();
	my $indels_gold = $self->project->gatk_indels_gold_file();
	my $vcf = $self->project->gatk_dbsnp_file();
	my $fileout = $filein;
	my $csv = $filein.".tmp1";
	 $fileout =~ s/bam/recal1\.bam/;
	my $proc ;
	if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 8};

	my $csv2 =$self->patient->getRecalFile();
	my $bed = $self->patient()->getCaptureBedFile() ;
	die("not found bed file : $bed") unless -e $bed;
	my $reference = $self->reference();
	my $cmd1 =  $self->gatk()." -I  $filein -R $reference -knownSites $vcf --interval_padding 100 -T BaseRecalibrator -o $csv2   -nct $proc --read_filter BadCigar ";
	$cmd1 .= " -knownSites $indels_gold" if ($indels_gold);
	$cmd1 .= " -L $bed  ";

	
	my $gatk_path = $self->project->getSoftware('gatk');
	my $cmd_gatk_v = "java -jar ".$self->project->getSoftware('gatk')." --version 2>&1";
	my $version_gatk= `$cmd_gatk_v`;
	$version_gatk =~ s/\n//g ;
	
	#trace
	$Header_lines->{'GATK_BaseRecalibrator'}->{'type'}="PG";
	$Header_lines->{'GATK_BaseRecalibrator'}->{'PN'}="BaseRecalibrator";
	$Header_lines->{'GATK_BaseRecalibrator'}->{'ID'}="GATK_BaseRecalibrator";
	$Header_lines->{'GATK_BaseRecalibrator'}->{'VN'}="$version_gatk" ;
	$Header_lines->{'GATK_BaseRecalibrator'}->{'CL'}="\"".$cmd1."\"" ;
	

	my $stepname = $name."_recalibration";
	$self->hash_steps->{$stepname}->{id}= 7 ;
	

	my ($job1,$job_next1) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd1,ppn=>$proc);
	if ($self->unforce() && -e $csv2){
		$self->add_skip_steps($stepname);
		return ($previous,$filein);
	}	
	if (defined $previous) {
		$self->prev(job=>$job1,prev=>$previous);
	}	
	$self->add_running_steps($stepname);
	return ($job1,$filein);
}

method realign_gatk  (Str :$filein,Object :$previous){
	my $name = $self->patient()->name();
	my $indels_gold = $self->project->gatk_indels_gold_file();
	unless ($filein) {
		my $outputdir = $self->patient()->getProject->getAlignmentPipelineDir($self->align_method);
		my $cmd = qq{};
		my @toto = `ls -r $outputdir/$name*recal.bam`;
		chomp(@toto);
		
		$filein = $toto[0];
		warn "work on : $filein";
		confess() unless -e $filein;
	}

	my $bed = $self->patient()->getCaptureBedFile();
	die("not found bed file : $bed") unless -e $bed;
	my $fileout = $filein;
	my $csv = $filein.".interval_list";
	$fileout =~ s/bam/realign\.bam/;
	
	my $proc1 ;
	if ($self->max_cpu()){
			$proc1 = $self->max_cpu ;
		}
		else {$proc1 = 16};
	 
	my $reference = $self->reference();
	my $cmd1 = $self->gatk()." -nt $proc1  -I $filein  -R $reference  -T RealignerTargetCreator -o $csv --interval_padding 100 ";
	$cmd1 .= " -known $indels_gold" if ($indels_gold);
	$cmd1 .= " -L $bed ";
	

	my $gatk_path = $self->project->getSoftware('gatk');
	#trace
	$Header_lines->{'GATK_RealignerTargetCreator'}->{'type'}="PG";
	$Header_lines->{'GATK_RealignerTargetCreator'}->{'PN'}="RealignerTargetCreator";
	$Header_lines->{'GATK_RealignerTargetCreator'}->{'ID'}="GATK_RealignerTargetCreator";
	$Header_lines->{'GATK_RealignerTargetCreator'}->{'VN'}="$gatk_path" ;
	$Header_lines->{'GATK_RealignerTargetCreator'}->{'CL'}="\"".$cmd1."\"" ;

	
	my $cmd2 = $self->gatk()." -I $filein -R $reference -T IndelRealigner  -targetIntervals $csv -o $fileout -rf NotPrimaryAlignment --interval_padding 100" ;
	$cmd2 .= " -known $indels_gold" if ($indels_gold);
	$cmd2 .= "";

	#trace
	$Header_lines->{'GATK_IndelRealigner'}->{'type'}="PG";
	$Header_lines->{'GATK_IndelRealigner'}->{'PN'}="IndelRealigner";
	$Header_lines->{'GATK_IndelRealigner'}->{'ID'}="GATK_IndelRealigner";
	$Header_lines->{'GATK_IndelRealigner'}->{'VN'}="$gatk_path" ;
	$Header_lines->{'GATK_IndelRealigner'}->{'CL'}="\"".$cmd2."\"" ;
	
	my $stepname1 = $name."_realign1";
	
	#info_step
	$self->hash_steps->{$stepname1}->{id}= 8;
	$self->hash_steps->{$stepname1}->{stepname}= "realign1";
		if ($self->unforce() && -e $fileout){
		$self->add_skip_steps($stepname1);
		return ($previous,$fileout);
	}	

	my ($job1,$job_next1) =  $self->construct_jobs(stepname=>$stepname1,cmd=>$cmd1,ppn=>$proc1);
	$self->add_running_steps($stepname1);
	#die();
	#if ($self->restart_step (filein=>$filein, fileout=>$fileout)){
		#$self->add_skip_steps($stepname1);
		#return ($previous,$fileout);
	#}
	my $stepname2 = $name."_realign2";
	$self->hash_steps->{$stepname2}->{id}= 9;
	my $proc ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 4};
	my ($job2,$job_next2) =  $self->construct_jobs(stepname=>$stepname2,cmd=>$cmd2,ppn=>$proc);	
	$self->prev(job=>$job2,prev=>$job1);
	my $stepname =$name."_realign2";
	if ($self->restart_step (filein=>$filein, fileout=>$fileout)){
		$self->add_skip_steps($stepname);
		return ($previous,$fileout);
	}
	
	if (defined $previous) {
		$self->prev(job=>$job1,prev=>$previous);
	}
	
	$self->add_running_steps($stepname);
	
	return ($job2,$fileout);
	
}



method haplotypeCaller (Str :$filein! ,Object :$previous! ){ 
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();	
	my $pname = $project->name();
	my $proc ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 4};
	my $dir_gvcf = $project->getRawCallingDir("unifiedgenotyper"); 
	my $fileout = $dir_gvcf."/".$name.".gvcf"; 

	  my $stepname = $name."_haplo";
	  warn $fileout;
	  if ($self->unforce() && -e $fileout){
	  	$self->add_skip_steps($stepname);
			return ($previous,$filein);
		}
	my $cmd = "perl $bin_dev/haplotypeCaller.pl -project=$pname -patient=$name ";
	
	  my ($job,$job_next) =  $self->construct_jobs_tmp(stepname=>$stepname,cmd=>$cmd,ppn=>$proc);
	 
		$self->add_running_steps($stepname);
	  if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	  return ($job,$filein);
}


method calling_merge  (Str :$filein! ,Object :$previous! ){ 
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();	
	my $project_name =$project->name();
	my $fileout = $filein;
	my $low_calling = "" ;
	$filein = $self->patient()->getBamFile() unless $filein;

	mkdir ( $project->getCallingPipelineDir("unifiedgenotyper")) unless -e  $project->getCallingPipelineDir("mpileup");
	my $dirout= $project->getCallingPipelineDir("unifiedgenotyper");
	$fileout = $dirout . "/" .$name.".final.vcf";
	my $proc ;
	my $log_file = $self->log_file;
	#probleme du fichier de log non défini dans le pipeline lancé sans pbs
#	my $cmd = "" ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
#			$cmd .= "perl $bin_dev/calling_illumina.pl -project=$project_name  -patient=$name -fork=$proc  -out=$fileout ;";
		}
		else {$proc = 8 ; 
#				my $log_file = $self->log_file;
#				$cmd .= "perl $bin_dev/calling_illumina.pl -project=$project_name  -patient=$name -fork=$proc -log=$log_file -out=$fileout ;";
		};
	die("-".$filein) unless $filein;
#	die($filein. " is empty") if (-z $filein);
	
	my $m = $self->patient()->alignmentMethod();
	my $dir_bam = $project->getAlignmentDir($m);
	my $recal_file = $dir_bam."/recal/".$name.".recal.table" ;
#	die ($recal_file. " is empty") if (-z $recal_file);
#die ($recal_file. " is absent ") unless (-e $recal_file);
	
	my $cmd = "perl $bin_dev/calling_illumina.pl -project=$project_name  -patient=$name -fork=$proc -log=$log_file -out=$fileout -calling=unifiedgenotyper,freebayes,samtools ";
	my $stepname = $name."_calling_merge";
	$self->hash_steps->{$stepname}->{id}=16 ;
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$proc);
	if ($self->unforce() && -e $fileout){
	  	$self->add_skip_steps($stepname);
			return ($previous,$filein);
		}
	$self->add_running_steps($stepname);
	if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	return ($job,$fileout);
}



method calling_gvcf  (Str :$filein! ,Object :$previous! ){ 
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();	
	my $project_name =$project->name();
	my $fileout = $filein;
	my $low_calling = "" ;
	$filein = $self->patient()->getBamFile() unless $filein;
	

my $dir_gvcf_out  = $project->getGvcfDir("haplotypecaller");
 	$fileout = $dir_gvcf_out."/".$name.".g.vcf.gz";
	my $proc ;
	my $log_file = $self->log_file;
	#probleme du fichier de log non défini dans le pipeline lancé sans pbs
#	my $cmd = "" ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 16 ; 
		};
	die("-".$filein) unless $filein;
#	die($filein. " is empty") if (-z $filein);
	
	
	my $cmd = "perl $bin_dev/calling_individual_gvcf.pl -project=$project_name  -patient=$name -fork=16 -log=$log_file -out=$fileout -window=1_000_000  ";
	my $stepname = $name."_gvcf";
	$self->hash_steps->{$stepname}->{id}=21 ;
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$proc);
	if ($self->unforce() && -e $fileout){
	  	$self->add_skip_steps($stepname);
			return ($previous,$filein);
		}
	$self->add_running_steps($stepname);
	if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	return ($job,$fileout);
}

method calling_merge_low  (Str :$filein! ,Object :$previous! ){ 
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();	
	my $project_name =$project->name();
	my $fileout = $filein;
	my $low_calling = "" ;
	$filein = $self->patient()->getBamFile() unless $filein;
	
	mkdir ( $project->getCallingPipelineDir("unifiedgenotyper")) unless -e  $project->getCallingPipelineDir("mpileup");
	my $dirout= $project->getCallingPipelineDir("unifiedgenotyper");
	$fileout = $dirout . "/" .$name.".final.vcf";
	my $proc ;
	my $log_file = $self->log_file;
	#probleme du fichier de log non défini dans le pipeline lancé sans pbs
#	my $cmd = "" ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
#			$cmd .= "perl $bin_dev/calling_illumina.pl -project=$project_name  -patient=$name -fork=$proc  -out=$fileout ;";
		}
		else {$proc = 16 ; 
#				my $log_file = $self->log_file;
#				$cmd .= "perl $bin_dev/calling_illumina.pl -project=$project_name  -patient=$name -fork=$proc -log=$log_file -out=$fileout ;";
		};
	die("-".$filein) unless $filein;
#	die($filein. " is empty") if (-z $filein);
	
	my $m = $self->patient()->alignmentMethod();
	my $dir_bam = $project->getAlignmentDir($m);
	my $recal_file = $dir_bam."/recal/".$name.".recal.table" ;
#	die ($recal_file. " is empty") if (-z $recal_file);
#die ($recal_file. " is absent ") unless (-e $recal_file);
	
	my $cmd = "perl $bin_dev/calling_illumina.pl -project=$project_name  -patient=$name -fork=$proc -log=$log_file -out=$fileout -low_calling -calling=unifiedgenotyper,freebayes,samtools ";
	my $stepname = $name."_calling_merge_low";
	$self->hash_steps->{$stepname}->{id}=16 ;
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$proc);
	if ($self->unforce() && -e $fileout){
	  	$self->add_skip_steps($stepname);
			return ($previous,$filein);
		}
	$self->add_running_steps($stepname);
	if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	return ($job,$fileout);
}

method calling_merge_varscan  (Str :$filein! ,Object :$previous! ){ 
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();	
	my $project_name =$project->name();
	my $fileout = $filein;
	$filein = 	$self->patient()->getBamFile() unless $filein;
	
	mkdir ( $project->getCallingPipelineDir("unifiedgenotyper")) unless -e  $project->getCallingPipelineDir("mpileup");
	my $dirout= $project->getCallingPipelineDir("unifiedgenotyper");
	$fileout = $dirout . "/" .$name.".final.vcf";
	my $proc ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 16};
	die("-".$filein) unless $filein;
	my $log_file = $self->log_file;
	my $cmd = "perl $bin_dev/calling_mp_ug_varscan.pl -project=$project_name  -patient=$name -fork=$proc -log=$log_file -out=$fileout ";
	my $stepname = $name."_calling_merge_varscan";
	$self->hash_steps->{$stepname}->{id}=17 ;
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$proc);
	if ($self->unforce() && -e $fileout){
	  	$self->add_skip_steps($stepname);
		return ($previous,$filein);
		}
	$self->add_running_steps($stepname);
	if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	return ($job,$fileout);
}



method move_vcf  (Str :$filein! ,Object :$previous! ){ 
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();	
	my $project_name =$project->name();	
	my $fileout = $filein;
	my $ppn = 1;	
	my $cmd = "perl $bin_dev/move_individual_vcf.pl -project=$project_name  -patient=$name -method=unifiedgenotyper ";
	my $stepname = $name."_move_vcf";
	$self->hash_steps->{$stepname}->{id}=18 ;
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$ppn);

	$self->add_running_steps($stepname);
	
	if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	return ($job,$filein);	
}


method move_snp_vcf  (Str :$filein! ,Object :$previous! ){ 
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();	
	my $project_name =$project->name();	
	my $fileout = $filein;
	my $ppn = 1;
	my $cmd = "$bin_dev/move_individual_vcf.pl -project=$project_name  -patient=$name -type=snp ";
		
	my $stepname = $name."_move_only_snp";
	
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$ppn);
	$self->add_running_steps($stepname);
	if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	return ($job,$filein);	
}




#non intégrée dans le pipeline
method picard_stats (Str :$filein!,Any :$previous){
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();
	my $project_name = $project->name();
	$filein = $self->patient()->getBamFile() unless $filein;
#	warn $filein ;
	my $stat_dir = $project->getMetricsDir();
	#my $stat_dir = $project->getRootDir() . "/align/stats/";
	mkdir $stat_dir unless -e $stat_dir;
	my $fileout = $self->patient->getMetricsFile();
#	warn $fileout ;
	#my $fileout = $stat_dir . "/" . $name.".stat";
	my $metricsHeader = $project->getMetricsHeader();
#	warn $metricsHeader;
	my $metricsFile = $self->patient->getCaptureMetricsFile($metricsHeader);
	my $bed = $self->patient()->getCaptureFile();
	my $picard = $self->picard_exec_path;
	my $javac = $self->javac;
	my $reference_bwa = $self->reference;
	my $cmd = qq{$picard CalculateHsMetrics BAIT_INTERVALS=$metricsFile TARGET_INTERVALS=$metricsFile INPUT=$filein OUTPUT=$fileout REFERENCE_SEQUENCE=$reference_bwa &&};

	my $stepname = $name."_stat";
	$self->hash_steps->{$stepname}->{id}=19 ;
	my $proc ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 16};
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$proc);

	  if ($self->unforce() && -e $fileout){
			return ($previous,$filein);
		}
		$self->add_running_steps($stepname);
		
	  if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}
	
	 return ($job,$filein);
	 
}

##methode remplacée par le script coverage.pl


method coverage_samtools (Str :$filein!,Any :$previous){
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();
	my $project_name = $project->name();
	$filein = $self->patient()->getBamFile() unless $filein;
	
	my $coverage_dir = $project->getRootDir() . "/align/coverage/";
	mkdir $coverage_dir unless -e $coverage_dir;
	my $fileout = $coverage_dir . "/" . $name.".cov.gz";
	my $bed = $self->patient()->getCaptureFile();
	
	my $stepname = $name."_coverage";
	$self->hash_steps->{$stepname}->{id}= 12;
	my $proc ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 16};
	
	### Récupération de la version de samtools
	my $samtools_path = $project->getSoftware('samtools') ;
	my $cmd2 ="$samtools_path 2>&1 |grep Version  " ;
	my $version_samtools = `$cmd2`;
	$version_samtools =~ s/Version:\ //g;
	$version_samtools="samtools-".$version_samtools ;
	$version_samtools =~ s/\n//g;
	
	my $cmd = qq{/usr/bin/perl $bin_dev/coverage.pl -patient=$name -filein=$filein -dir=$coverage_dir -bed=$bed -fork=$proc -name=$name -project=$project_name };
	
	
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$proc);
	if ($self->unforce() && -e $fileout){
		return ($previous,$filein);
		}
	$self->add_running_steps($stepname);
		
	 if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	 return ($job,$filein);
}






method rehead_bam (Str :$filein!,Object :$previous){
	$filein = $self->patient()->getBamFile() unless $filein;
	my $fileout = $filein ;	
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();
	my $project_name = $project->name();
	warn $filein;
	my $cmd = "";
	
	#rehead bam
	
	$fileout =~ s/bam/header\.bam/;
	my $fileout_header =$fileout.".head";
	

	
	my $fileout_bai = $fileout;
	$fileout_bai =~ s/bam/bam.bai/;
	my $stepname = $name."_rehead_bam";
	
		if ($self->restart_step (filein=>$filein, fileout=>$self->goal)){
	
		$self->add_skip_steps($stepname);
		return ($previous,$self->goal);
	}
	my @echo_arg=("ID","PN","VN", "CL","PU","SM","DS","PL","PG");
	my $echo_line = "" ;
	foreach my $key (keys %$Header_lines){
		unless ($Header_lines->{$key} eq 'Picard_MarkDuplicates'){
			$echo_line .= "\n\@".qq{$Header_lines->{$key}->{"type"}};
			foreach my $arg (@echo_arg){
				if (exists $Header_lines->{$key}->{$arg}){
					$echo_line .= "\t".qq{$arg:$Header_lines->{$key}->{$arg} } ;
				}
			}			
		}
	}
	my @echo_line_list = split("",$echo_line);
	shift(@echo_line_list); #supprime le premier \n
	$echo_line = join("",@echo_line_list);
	$echo_line =~ s/\"/\\\"/g;
	$echo_line =~ s/&&/ PUIS /g;
	my $samtools = $project->buffer->software("samtools");
	$cmd = "echo -e '$echo_line' | perl $bin_dev/./uniqflag_bam_header.pl -file=$filein -file_out=$fileout -project=$project_name";
	warn "------";
	warn $cmd."\n";
	system("$cmd");
	#die();
	#$cmd .="/usr/bin/perl $bin_dev/./uniqflag_bam_header.pl -file=$filein  $samtools view -H  $fileout > $fileout_header && echo -E '$echo_line' >> $fileout_header && $samtools reheader $fileout_header $fileout > $fileout_reheaded && mv $fileout_reheaded $fileout && rm $fileout_header && $samtools index $fileout &&" ;
	

	#$cmd .="$samtools view -H  $fileout > $fileout_header && echo -E '$echo_line' >> $fileout_header && $samtools reheader $fileout_header $fileout > $fileout_reheaded && mv $fileout_reheaded $fileout && rm $fileout_header && $samtools index $fileout &&" ;

	
	
	$self->hash_steps->{$stepname}->{id}= 11;
	my $proc ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 4};
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$proc);
	
	$self->add_running_steps($stepname);
	if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	 return ($job,$fileout);
}

####Méthodes générales de manipulation des bam

method bam_fastq  (Str :$filein,Object :$previous){
	my $name = $self->patient()->name();
	my $fileout = $filein;
	$fileout =~ s/bam/sort\.bam/;
	my $dir1 ="/data-xfs/sequencing/ngs/NGS2014_0528/HG19/align/bwa/fastq/";
	my $fastq1 =$dir1.$name."_1.fastq";
	my $fastq2 =$dir1.$name."_2.fastq"; 
	 
	my $picard_path = $self->project->getSoftware('picard');
	my @picard= split( "/", $picard_path);
	my $version_picard= $picard[-1];
	my $comment = "multilane_merge_picard\\tpicard version : $version_picard";
	
	my $proc ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 8};
	 
	my $cmd =  $self->picard_exec_path()."/SamToFastq.jar    VALIDATION_STRINGENCY=SILENT INPUT=".$filein." FASTQ=".$fastq1." SECOND_END_FASTQ=".$fastq2."&& gzip $fastq1 && gzip $fastq2 " ;
	my $stepname = $name."_fastq";
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$proc); #illumina_util::construct_jobs($log_file,$stepname,$cmd,4);
	if ($self->restart_step (filein=>$filein, fileout=>$fileout)){
		$self->add_skip_steps($stepname);
		return ($previous,$fileout);
	}
	
	$self->add_running_steps($stepname);
	if (defined $previous) {
			$self->prev(job=>$job,prev=>$previous);
	}	
	return ($job,$filein);
}

method delete_chrM (Str :$filein,Object :$previous){
	my $name = $self->patient()->name();
	my $fileout = $filein;
	 $fileout =~ s/bam/nochrM\.bam/;
	 my $filetmp = $filein."tmp" ;
	my $cmd = $self->samtools()." view -h $filein && grep -v chrM $filein > $filetmp && samtools view -Sb $filetmp > $fileout &&".$self->samtools()." index $fileout && rm $filetmp 	";
	my $stepname = $name."_delete_chrM";
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>4); #illumina_util::construct_jobs($log_file,$stepname,$cmd,4);
	if ($self->restart_step (filein=>$filein, fileout=>$fileout)){
		$self->add_skip_steps($stepname);
		return ($previous,$fileout);
	}
	
	$self->add_running_steps($stepname);
	if (defined $previous) {
			$self->prev(job=>$job,prev=>$previous);
	}	
	return ($job,$fileout);
}

method chrM_chrMT  (Str :$filein,Object :$previous){
	my $name = $self->patient()->name();
	my $fileout = $filein;
	$fileout =~ s/bam/chrM\.bam/;
	my $filetmp = $filein."tmp"  ;
	my $cmd =  $self->samtools()." view -h  $filein >$filetmp && sed s/chrMT/chrM/ $filetmp && ".$self->samtools()." view -bS $filetmp > $fileout && rm $filetmp ";
	my $stepname = $name."_chrM";
	
	my $proc ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 8};
		
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$proc); #illumina_util::construct_jobs($log_file,$stepname,$cmd,4);
	if ($self->restart_step (filein=>$filein, fileout=>$fileout)){
		$self->add_skip_steps($stepname);
		return ($previous,$fileout);
	}
	$self->add_running_steps($stepname);
	if (defined $previous) {
			$self->prev(job=>$job,prev=>$previous);
	}	
	return ($job,$fileout);
}

method add_chr (Str :$filein!,Any :$previous){ 
	my $name = $self->patient()->name();
	unless ($filein){
		$filein = $self->patient()->getBamFile();
	}
	my $cmd = "$bin_dev/add_chr_to_bam.pl -file $filein ";

	my $stepname = $name."_add_chr";
	my $proc ;
		if ($self->max_cpu()){
			$proc = $self->max_cpu ;
		}
		else {$proc = 4};
	 my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$proc);
	 $self->add_running_steps($stepname);
	  if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	 return ($job,$filein);	 
}


1;
