package job_bds;

use Moo;

use Data::Printer;
use FindBin qw($Bin);
use Time::Local;
use POSIX qw(strftime);
use Term::ANSIColor;
use Data::Dumper;
use Carp;
use colored;

has 'cmd' => (
	is =>'ro',
);
has   'isLogging' =>(
	is =>'rw',
	default => sub {
		return 0;
	}
);
has 'ppn' => (
	is =>'ro',
);
has 'filein' => (
	is =>'ro',
);

has 'fileout' => (
	is =>'ro',
);

has 'type' => (
	is =>'ro',
);


has 'name' => (
	is =>'ro',
);


has 'prev' => (
	is =>'rw',
	default=> sub {
		return [];
	},
);

has 'next' => (
	is =>'rw',
	default=> sub {
		return [];
	},
);

has 'dir_bds' =>(
		is =>'ro',
);

#has 'is_skip' => (
#	is =>'rw',
#	default=> sub {
#		0;
#	},
#);


has 'category' => (
	is =>'rw',
	lazy=>1,
	default=> sub {
		my $self = shift;
		 my ($t,$nb) = split("\#",$self->type);
		 return $t;
	},
);

has 'bds_log' => (
	is =>'rw',
	lazy=>1,
	default=> sub {
		my $self = shift;
		die() unless -e $self->dir_bds;
		return $self->dir_bds."/".$self->name.".log";
	},
);
has 'bds_ok' => (
	is =>'rw',
	lazy=>1,
	default=> sub {
		my $self = shift;
		die() unless -e $self->dir_bds;
		return $self->dir_bds."/".$self->name.".ok";
	},
);

has 'bds_start' => (
	is =>'rw',
	lazy=>1,
	default=> sub {
		my $self = shift;
		return $self->dir_bds."/".$self->{name}.".start";
	},
);
has 'bds_error' => (
	is =>'rw',
	lazy=>1,
	default=> sub {
		my $self = shift;
		return $self->dir_bds."/".$self->{name}.".error";
	},
);
has 'run_log' => (
	is =>'rw',
	lazy=>1,
	default=> sub {
		my $self = shift;
		return $self->dir_bds."/".$self->{name}.".log.error";
	},
);

sub  is_running {
	my ($self) = @_;
	return if $self->is_finished();
	return -e $self->bds_start();
	return;
}

sub  is_finished {
		my ($self) = @_;
	if (-e $self->bds_error or -e $self->bds_ok){
		return 1;
	}
	return undef;
	#return  -e $self->bds_error or -e $self->bds_ok;
	
}

sub  is_pending {
		my ($self) = @_;
	if (!{$self->is_finished} && !($self->is_running)){
		return 1;
	}
	return;
}

sub  is_ok {
		my ($self) = @_;
	return  -e $self->bds_ok;
}

sub  is_error {
		my ($self) = @_;
	return  -e $self->bds_error;
}

sub is_root {
		my ($self) = @_;
		return 1 if scalar(@{$self->prev}) ==0;
}
sub is_leaf{
		my ($self) = @_;
		return 1 if  scalar(@{$self->next}) == 0;
}

sub is_first{
		my ($self) = @_;
		return undef  if $self->is_skip;
		return 1 if $self->is_root;
		my @prevs = grep{$_->is_run} @{$self->prev};
		return 1 unless @prevs;
		return;
}




sub command_bds{
		my ($self) = @_;
	my @c;
	push(@c,"date >  ".$self->bds_log);
	push(@c," touch ".$self->bds_start);
	foreach my $c (@{$self->cmd}){
		my $f = $self->run_log();
		my @t = split("&&",$c);
		 $c = join("2>>$f && ",@t);
		 $c = "$c >> $f";
		#$c = "sleep 1 ;echo ".$self->name;
		push(@c," ($c 2>>$f)  || ( touch ".$self->bds_error." ; echo 'log file here :".$f."')");
	}
	push(@c," test -e ".$self->bds_error ." &&  rm ".$self->fileout);
	push(@c," test -e ".$self->fileout ." || ( rm ".$self->bds_log."; touch ".$self->bds_error.") && touch ". $self->bds_ok." ;  sleep 5 ; ls ".$self->bds_log ." 2>/dev/null ");
	return @c;
}

sub print_task {
		my ($self) = @_;
	my $string;
	
	if ($self->is_skip){
		return "sys date > ".$self->bds_log."\n";
	}
	
	my @filein;
	if ($self->is_first){
		@filein = @{$self->filein};
		#	 $dep  =join("\",\"",@{$self->filein})."\"]" ;
	}
	else {
			@filein = map {$_->bds_log} @{$self->prev};
		
	}
	 my $dep  = "\"".$self->bds_log."\" <- [\"".join("\",\"",@filein)."\"]" ;
	 my $taskname = $self->name;
	 my $cpu = $self->ppn;
	 
  	$string = "task ( $dep , taskName := \"$taskname\" , cpus := $cpu ) {\n";	
  	foreach my $c ($self->command_bds){
  		$string .= "sys $c \n";
  	}
  	$string .= "\n}\n\n\n";
}


sub is_run{
		my ($self) = @_;
	return 1  if  $self->is_skip == 0;
	return 0;
}
sub is_type{
	my ($self,$hash) = @_;
	confess();
	my $type = $hash->{type};
	return  $self->type =~/$type/;
} 

sub print_start_status {
		my ($self) = @_;
	if ($self->is_skip){
	return colored::stabilo("white","SKIP",1);
	}
	else {
		return colored::stabilo("green","TODO",1);
	}
	
}
sub add_prev {
	my ($self,$hash) = @_;
	my $previous = $hash->{previous};
	return unless defined $previous;
	push(@{$self->prev},$previous);
} 

sub print_prev {
		my ($self) = @_;
	foreach my $p (@{$self->prev}){
			warn "\t".$p->name."\n";
	}
}
sub add_next{
		my ($self,$hash) = @_;
		confess();
	my $next = $hash->{next};
	return unless defined $next;
	push(@{$self->next},$next);
} 
sub is_skip {
	my ($self,$skip) = @_;
	
	if ($skip ){
		$self->{skip} = $skip;
	}
	return $self->{skip} if exists $self->{skip};
	 $skip =1;
	foreach my $j (@{$self->next}){
		if (!($j->is_skip())){
			$skip =0;
		}
	}
$skip =0 if scalar@{$self->next}== 0;
$self->{skip} = $skip;
return $self->{skip};
	
}

#	is =>'rw',
#	default=> sub {
#		0;
#	},
#);
sub skip {
		my ($self) = @_;
	#return undef  unless $self->prev();
	confess() if $self->name =~/replace/;
	$self->is_skip(1);
	return if $self->is_root();
	foreach my $p (@{$self->prev}){
		$p->skip();
	}
	#die();
}
sub is_prod (){
	my ($self) = @_;
	return 1 if $self->fileout !~/pipeline/;
	return undef;
}




1;
#{cmd=>$cmds,name=>$stepname,ppn=>$ppn,filein=>$filein,fileout=>$fileout,type=>$type}