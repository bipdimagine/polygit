#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use Net::SSH::Perl; 
use Getopt::Long;
use Proc::Simple;
my $cmd;
my $pwd;
GetOptions(
	'cmd=s' => \$cmd,
	'pwd=s' => \$pwd,
);
my $file_lock = "/data-dragen/lock";
my $last_cmd;
my $ip_dragen = "10.1.2.10";
my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);
 my $ssh = Net::SSH::Perl->new($ip_dragen);
my $myproc = Proc::Simple->new(); 
$SIG{INT} = \&tsktsk;
$SIG{KILL} = \&tsktsk;
sub tsktsk {
	warn "kill jobs !!!!!!! wait 10 seconds";
  	$myproc->kill() ;#if $myproc->poll();
  	undef $myproc;
  	if($last_cmd){
  		sleep(10);
  		$last_cmd =~ s/  / /g;
  		$last_cmd =~ s/^\s+|\s+$//g;
  		my ($o,$e,$ee) = $ssh->cmd(qq{pkill -f \"$last_cmd\"});
  		warn qq{pkill -f \"$last_cmd\"};
  		unlink $file_lock if -e $file_lock;
  		warn " $ee $e +++++"
  	}
  	
  
}


$ssh->login("$username",$pwd);
my @cmds;

unless ($cmd){
  @cmds = <STDIN>;
	chomp(@cmds);
}
else {
  push(@cmds,split(",",$cmd));	
}

foreach my $dragen_cmd (@cmds) {
	$last_cmd = $dragen_cmd;
	my $exit = run_cmd($dragen_cmd);
	unlink $file_lock if -e $file_lock;
	die($dragen_cmd) if $exit != 0;
	$last_cmd = undef;
}

warn "end";

sub run_cmd {
	my ($cmd) = @_;
	die() if $cmd !~ /^dragen/;
	my $f =0;
	my $sleep = 5;
	my $unlock;
	while (-e $file_lock){
		my ($line) =  `cat $file_lock`;
	
		warn "lock by : ".$line if $f == 0;
		$f ++ ;
		sleep($sleep);
		if ($f%2 == 0){
			chomp($line);
			my ($user1,$cm,$t) = split(":",$line);
			my($stdout, $stderr, $exit) = $ssh->cmd(qq{ps -edf | grep $user1 | grep dragen | grep -v "grep"});
			 my @cmd = split("\n",$stdout);
			if($unlock){
				warn "delete lock";
				unlink $file_lock;
			}
			unless (@cmd){
				$unlock =1;
			}
		}
	}
	die($file_lock) if -e $file_lock;
	open(LOCK,">$file_lock");
	print LOCK "$username:$cmd:".time."\n";
	close(LOCK);
	system("chmod a+rwx $file_lock");
	#warn "run";
	$myproc = Proc::Simple->new(); 
	$myproc->kill_on_destroy(1);            # Set kill on destroy
	$myproc->signal_on_destroy("KILL");
	$myproc->start(qq{ssh $username\@}.qq{$ip_dragen $cmd});
	while ($myproc->poll()){
		sleep 10;
		
	}
	my $exit = $myproc->exit_status();
	
	$ssh->cmd("rm $file_lock");
	if ($exit == 0){
		
		return 0;
	}
	else {
		return $exit;
	}	
	
	
}
