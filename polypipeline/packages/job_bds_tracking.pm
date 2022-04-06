package job_bds_tracking;
use Moose;
use MooseX::Method::Signatures;
use strict;
use FindBin qw($Bin);

use lib "$Bin";
extends (qw(job_bds));


has 'software' => (
	is =>'ro',

);

has 'sample_name' => (
	is =>'ro',

);
has 'project_name' => (
	is =>'ro',

);

has 'pipeline_id' => (
	is =>'ro',

);

has   'script_tacking_dir' =>(
	is =>'ro',
	default => sub {
		return qq{$Bin/scripts/scripts_tracking};
	}
);
has   'uuid' =>(
	is =>'ro',
	required=>1,
);
sub return_start_tracking {
	my($self) = @_;
	my $step_name = $self->name;
	my $prog = $self->software;
	my $cmd_line = $self->cmd->[0];
	my $project_name = $self->project_name();
	my $name = $self->sample_name();
	my $uuid = $self->uuid;
	die() unless $uuid;
	my $cmd = "perl ".$self->script_tacking_dir."/start_step_sqlite.pl -version= -status=running -run_id=$uuid -project=$project_name -patient=$name -step=\"$step_name\" -prog=$prog -cmd=\"$cmd_line\"";
	return $cmd;
}


sub return_end_tracking {
	my($self) = @_;
	my $step_name = $self->name;
	my $prog = $self->software;
	my $name = $self->sample_name();
	my $project_name = $self->project_name();
	my $cmd_line = $self->cmd->[0];
	my $uuid = $self->uuid;
	my $cmd = "perl ".$self->script_tacking_dir."/end_step_sqlite.pl -status=finished -run_id=$uuid -project=$project_name -patient=$name -step=\"$step_name\"";
	return $cmd;
}

method command_bds (){
	my @c;
	push(@c,"date >  ".$self->bds_log);
	push(@c," touch ".$self->bds_start);
	foreach my $c (@{$self->cmd}){
		my $f = $self->run_log();
		my $start =  $self->return_start_tracking().">$f";
		my $end  =  $self->return_end_tracking().">>$f";
		my @t = split("&&",$c);
		 $c = join("2>>$f && ",@t);
		 
		 
		push(@c," ($start 2>>$f && $c 2>>$f && $end 2>> $f)  || ( touch ".$self->bds_error." ; echo 'log file here :".$f."')");
	}
	push(@c," test -e ".$self->bds_error ." &&  rm ".$self->fileout);
	push(@c," test -e ".$self->fileout ." || ( rm ".$self->bds_log."; touch ".$self->bds_error.") && touch ". $self->bds_ok." ;  sleep 5 ; ls ".$self->bds_log ." 2>/dev/null ");
	return @c;
}



1;