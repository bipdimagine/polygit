package sample;
use Moose;
use MooseX::Method::Signatures;
use Data::Printer;
use FindBin qw($Bin);
use Time::Local;
use POSIX qw(strftime);#
use List::Util qw( uniq );
use Time::HiRes;
has 'name' => (
	is =>'ro',
);

has 'priority' => (
	is =>'ro',
	default=> sub {
		return 1;
	},
);

has 'wait_after' => (
	is =>'rw',
	default=> sub {
		return undef;
	},
);

has 'index' => (
	is =>'rw',
	default=> sub {
		return -1;
	},
);

has 'task_name' => (
	is =>'ro',
	lazy=>1,
	default=> sub {
		my $self =shift;
		my $fname = $self->name();
		$fname =~s/\-/_/g;
		$fname =~s/\./_/g;
		$fname="F_".$fname."_".int(rand(1000000));
		return $fname;
	},
);


has 'jobs' => (
	is =>'ro',
	default=> sub {
		return [];
	},
);
has 'categories' => (
	is =>'ro',
	default=> sub {
		return {};
	},
);
has 'jobs_type' => (
	is =>'ro',
	default=> sub {
		return {};
	},
);

has 'current_bam' => (
	is =>'rw',
	default=> sub {
		return undef;
	},
);

method list_categories (){
	 my @cat;
	 $self->reconstruct_dependancy();
	foreach my $j (@{$self->jobs}){
		#warn $j->name." ==> ".$j->{ppn};
		push(@cat,$j->category);
	}
	return uniq(@cat);
}

method list_hashes_categories (){
	 my @cat;
	 $self->reconstruct_dependancy();
	 my $hh;
	foreach my $j (@{$self->jobs}){
		#warn $j->name." ==> ".$j->{ppn};
		my $n = $j->category;
		push(@{$hh->{$n}->{jobs}},$j);
		$hh->{$n}->{ppn}  =0 unless exists $hh->{$n}->{ppn} ;
		$hh->{$n}->{ppn} = $j->{ppn} if  $j->{ppn} > $hh->{$n}->{ppn};
	}
	return $hh;
}
method get_jobs_category(Any :$category!){
	return $self->categories->{$category};
} 


method  delete_tmp_files() {
	my @tmp;
	foreach  my $j (@{$self->jobs}){
		next if $j->is_prod();
		next if $j->is_skip;
		push(@tmp,$j->fileout);
	}
	
	return \@tmp;
}

method  delete_prod_files() {
	my @tmp;
	foreach  my $j (@{$self->jobs}){
		next unless $j->is_prod();
		
		next if $j->is_skip;

		push(@tmp,$j->fileout);
	}
	return \@tmp;
}

method add_job(Any :$job!){
	push(@{$self->jobs},$job);
	push(@{$self->categories->{$job->category}},$job);
	$self->jobs_type->{$job->type} = $job;
	if ($job->fileout =~/\.bam/){
		$self->current_bam($job->fileout);
	}
	 delete $self->{reconstruct};
} 


method is_pending_jobs(){
	my @jobs = grep{$_->is_run} @{$self->jobs};
	return 1 if (@jobs);
	return undef;
}

method bds_sub (){
	$self->reconstruct_dependancy();
	my $jobs = $self->jobs;
	my $fname = $self->task_name;
	my $string = " void $fname () {\n";
	foreach my $j (@$jobs){
		$string.=$j->print_task;
		#$string.="\nsleep(2);\n";
	}
	$string.="}\n";
}


method reconstruct_dependancy (){
	return if exists $self->{reconstruct};
	my @un = grep{$_->is_leaf}  @{$self->jobs};
	my %fileout;
	my %filein;
	foreach my $j (@{$self->jobs}) {
		 $fileout{$j->fileout} = $j;
#	confess();
		 map{$filein{$_} = $j} @{$j->filein};
		 $j->next([]);
		  $j->prev([]);
	}

	foreach my $j (@{$self->jobs}) {
		my $fout = $j->fileout;
		my $jnext = $filein{$fout};
		$j->next([$jnext]) if $jnext;
		foreach my $fi (@{$j->filein}){
		
			next unless exists $fileout{$fi};
			$j->add_prev(previous=>$fileout{$fi});
		}
	}
}

method print_jobs() {
	$self->reconstruct_dependancy();
	my $jobs = $self->jobs;
	foreach my $j (@$jobs){
		warn $j->name();
	}
	
}

method log_for_error_jobs (){
	my $jobs = $self->jobs;
	my @error =   grep{$_->is_error()} @$jobs;
	my @files;
	foreach my $e (@error){
		push(@files,$e->run_log);
	}
	return \@files;
}

method status_jobs(){
	my $jobs = $self->jobs;
	my $nb = scalar( @$jobs);
	my @running =   grep{$_->is_running()} @$jobs;
	my @skip =   grep{$_->is_skip()} @$jobs;
	my $nbs = scalar(@skip);
	my @ok =   grep{$_->is_ok()} @$jobs;
	my @finished =   grep{$_->is_finished()} @$jobs;
	my $nb_finished = scalar(@finished);
	my $nbok = scalar(@ok);
	my @error =   grep{$_->is_error()} @$jobs;
	my $nberror = scalar(@error);
	my $color = "white";
	my $tcurrent = colored::stabilo("white","Hold",1);
	my $status = "Hold";
	my $remaining = $nb;
	my $trunning = colored::stabilo("white",0,1); 
	if (@running){
		my %t;
		map{$t{$_->category}++} @running;
		my $text ="";
		foreach my $k (keys %t){
			$text .= $k."-";
		}
		chop($text);
		$tcurrent = colored::stabilo("yellow", $text,1);
		$trunning = colored::stabilo("yellow",scalar(@running),1);
	}
	my $tskip ="-";
	if (@skip){
		$tskip = $nbs;
		
	}
	my $terror = colored::stabilo("green","0",1);
	if (@error){
		my %t;
		map{$t{$_->category}++} @error;
		$terror  =  colored::stabilo("red", join(",",keys %t),1);
		$status = colored::stabilo("red","Error",1);
	}
	my $nremaining = $nb - $nb_finished; 
	if ($nremaining == 0 ){
			$remaining = colored::stabilo("green","Done",1);
			$status = colored::stabilo("green","Done",1);
	}
	$nb_finished += scalar(@skip);
	if ($nb_finished < $nb ){
			$remaining =  colored::stabilo("magenta",$nb_finished."/".$nb,1);
	}
	return [$tcurrent,$trunning,colored::stabilo("green",$nbok,1),$terror,$remaining];
}


1;