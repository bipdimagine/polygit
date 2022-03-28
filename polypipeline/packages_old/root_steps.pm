package root_steps;

use Moose;
use MooseX::Method::Signatures;
use Data::Printer;
use FindBin qw($Bin);
use Time::Local;
use Term::ANSIColor;


my $picard = qq{/bip-d/soft/distrib/picard/picard-tools-1.56/};
my $javac = qq{/opt/java/latest/bin/java  -Xmx20g -XX:ParallelGCThreads=8};
my $bwa = "/bip-d/soft/bin/bwa";
my $samtools = qq{/bip-d/soft/bin/samtools};
my $reference = qq{/data-xfs/public-data/HG19/genome/fasta/all.fa};
my $reference_bwa = qq{/data-xfs/public-data/HG19/genome/bwa-latest/all.fa};
my $gatk = qq{/bip-d/soft/distrib/GATK/gatk-latest};
my $rod = qq{-B:dbsnp,VCF /data-xfs/public-data/HG19/snp/dbsnp/132/dbsnp132_chr_nopar.vcf};
my $lifescope_dir = "/data-xfs/lifescope/projects/liferoot/";
my $vcf = qq{/data-xfs/public-data/HG19/snp/dbsnp/gatk/snp.vcf.gz};
 
has 'patient' => (
	is =>'rw',
	isa =>"Object",
);


has 'argument_patient' => (
	is =>'rw',
	isa =>"Str",
);
has 'project' => (
	is =>'rw',
	isa =>"Object",
);

has 'jobs_bds' => (
	is =>'rw',
 	isa       => 'ArrayRef',
     default   => sub{[]},
);

has 'chromosome' => (
	is =>'rw',
	isa =>"Object",
);

has 'log_file' => (
      is        => 'ro',
      isa       => 'Str',
      #required => 1,
     # lazy =>1,
      #default   => \&init_log_file,
); 

has 'log_option' => (
      is        => 'ro',
      isa       => 'Str',
      
      #required => 1,
     # lazy =>1,
      default   => sub {
      	my $self =shift;
      	"-log=".$self->log_file;
      },
); 
has 'unforce' => (
      is        => 'rw',
      isa       => 'Int',
      #required => 1,
     # lazy =>1,
      default   => 1,
); 
has 'running_steps' =>(
	  is        => 'rw',
      isa       => 'HashRef',
       default   => sub{{}},
#  
);
has 'skip_steps' =>(
	  is        => 'rw',
      isa       => 'HashRef',
       default   => sub{{}},
#  
);

has 'align_method'  =>(
	  is        => 'rw',
      isa       => 'Str',
      lazy => 1,
       default   => \&init_method_align,
#  
);
has 'align_method'  =>(
	  is        => 'rw',
      isa       => 'Str',
      lazy => 1,
       default   => \&init_method_align,
#  
);

my $default_file = {
	rmdup => ".1.bam"
};

#my $running_steps;

method init_method_align {
	my $methods = $self->patient()->getProject->getAlignmentMethods();
	
 warn ("more than one alignmeent method use default mapread") if  	scalar(@$methods) >1;
  return "mapreads"  if  	scalar(@$methods) >1;
  return $methods->[0];
}

method start_job {
	my $date = `date`;
	chomp($date);
	my $log_file = $self->log_file;
	system ("echo ===================== $date > $log_file");
	my $job = $self->LoggingJobs(stepname=>$self->project->name,text=>"START :".$date);
	return $job;
}

method end_job {
	my $date = `date`;
	chomp($date);
	my $job = $self->LoggingJobs(stepname=>$self->project->name,text=>"END :".$date);
	return $job;
}
 
method init_log_file {
	my $dir = $self->patient->getProject()->getPipelineDir();
	my $file = $dir."/".$self->patient->name().".log";
	unlink($file) if -e $file;
	return $file;
}

method bds_header() {
	return qq{#!/usr/bin/env bds};
}

method bds_system(Any :$job!) {
	my $fork = $job->{ppn};
	my $name = $job->{name};
	my $cmd = $job->{cmd};
	my $tid = "tid_".$job->{id};
	my $ph = "-log=".$self->log_file;
	$cmd =~s/$ph//;
	my $str = "";
	 $str = qq{$tid := task (cpus := $fork, taskName := "$name" ) \{\n};
		#sys echo "start $i => $n"
      $str.= qq{ sys $cmd;\n};
       $str.= qq{\}\n};
       return $str;
		#sys echo "end $i => $n"	  
}
method bds_wait(Any :$job!) {
	my $name = $job->{name};
	my $cmd = $job->{cmd};
	 return qq{\n wait; print "\\n ============================\\n =====> $name  FINISH \\n ============================\\n" \n};
	
}
method bds_wait_and_exit(Any :$job!) {
	my $name = $job->{name};
	my $cmd = $job->{cmd};
	my $tid = "tid_".$job->{id};
	 return qq{\n wait; 
	 	if ($tid.isDoneOk()){
	 	print "\\n ============================\\n =====> $name  FINISH \\n ============================\\n" ;
	 	}
	 	else {
	 		print "\\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\\n =====> $name  ERROR \\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\\n" ;
	 		exit(1)
	 	}
	 	};
	
}
method LoggingJobs (Str :$stepname!,Str :$text!){ 
 my $log_file = $self->log_file;

my $commandOK = qq{echo `date` "$stepname\t $text  " >> $log_file};
$stepname.= "_OK" if $text eq "OK";
$stepname.= "_ERROR" if $text eq "ERROR";
			my $jobsOK    = PBS::Client::Job->new(

				# Job declaration options
				# Resources options
				ppn    => 1,                           # process per node
				nodes  => 1,
				script => $stepname,
				cmd    => [$commandOK],                # command to be submitted

			);
	return $jobsOK;
}
method add_running_steps($name){
	
	push(@{$self->running_steps->{$self->project->name}},$name);
}
method add_skip_steps($name){
	push(@{$self->skip_steps->{$self->project->name}},$name);
}
method print_steps($name?) {
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


method prev (Object :$job! , Object :$prev!){
	 push((@{$job->{prev}->{ok}},$prev));
	
}

method prevArray (Object :$job! , Any :$prev!){
	$job->{prev}->{ok} = $prev;
	# push((@{$job->{prev}->{ok}},$prev));
	
}

method add_job(Str :$cmd!,Str :$stepname!, Int :$ppn!,Str :$fileout! ,Str :$filein!){
	my $job ={ cmd => $cmd,
				name => $stepname,
				ppn=>$ppn,
				fileout => $fileout,
	};
	push(@{$self->jobs},$job);
}

method add_wait(Str :$stepname!){
	my $job ={ 
				cmd =>"wait",
				name => $stepname,
	};
	push(@{$self->jobs},$job);
}
method construct_jobs (Str :$cmd!,Str :$stepname!, Int :$ppn!) {
my $log_file = $self->log_file;

my $start = qq{echo  "\t\t** start =>  $stepname :" `hostname` ": $ppn  : " `date` >> $log_file ; };

	my $jobs = PBS::Client::Job->new(
         # Job declaration options
         # Resources options
       	 script=>$stepname,
         ppn       => $ppn,             # process per node
         nodes => 1,
       
         cmd       => [$start.$cmd],# command to be submitted 
 	);
 
 	my $jobsOK = $self->LoggingJobs(stepname=>$stepname,text=>"OK");
 	
 	my $jobsERROR = $self->LoggingJobs(stepname=>$stepname,text=>"ERROR");
 	
 	$jobs->next( { ok => [$jobsOK], fail => $jobsERROR } );
 	return ($jobs,$jobsOK);
}

sub return_id{
	my $self = shift;
	$self->{id} =0 unless exists $self->{id};
	$self->{id} ++;
	return $self->{id};
}
my $cid = 0;
method construct_jobs_bds (Str :$cmd!,Str :$stepname!, Int :$ppn!,Str :$type!) {
my $log_file = $self->log_file;

my $start = qq{echo  "\t\t** start =>  $stepname :" `hostname` ": $ppn  : " `date` >> $log_file ; };
my $error = qq{ || echo  "\t\t** ERROR =>  $stepname :" `hostname` ": $ppn  : " `date` >> $log_file;exit 1 };
	my $jobs = PBS::Client::Job->new(
         # Job declaration options
         # Resources options
       	 script=>$stepname,
         ppn       => $ppn,             # process per node
         nodes => 1,
       
         cmd       => [$start.$cmd." ".$self->log_option.$error],# command to be submitted 
 	);
 
 	#my $jobsOK = $self->LoggingJobs(stepname=>$stepname,text=>"OK");
 	
 	my $jobsERROR = $self->LoggingJobs(stepname=>$stepname,text=>"ERROR");
 	my $cid  =$self->return_id() ;
 		my $job ={ cmd => $cmd,
				name => $stepname,
				ppn=>$ppn,
				id=>$cid,
				
			#	fileout => $fileout,
	};
	
 	if ($type eq "wait"){
	my $jobw ={ 
				cmd =>"wait",
				name => $stepname."_wait",
				id=>$cid,
				wait=>1,
	};
	push(@{$self->jobs_bds},$jobw);
 	}
 	elsif ($type eq "wait_exit"){
 		
	my $jobw ={ 
				cmd =>"wait_exit",
				name => $stepname."_wait",
				id=>$cid,
				wait_exit=>1,
	};
	push(@{$self->jobs_bds},$job);
	push(@{$self->jobs_bds},$jobw);
 	}
 	else {
 		push(@{$self->jobs_bds},$job);
 	}
	
	
 	
 	
 	#$jobs->next( { fail => $jobsERROR } );
 	return ($jobs);
}

method run_step (Str :$cmd!, Str :$stepname, Int :$ppn, Any :$previous, Str :$fileout, Str :$filein) {
	
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$ppn);
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
 
 
 1;