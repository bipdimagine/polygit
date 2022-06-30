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

my $myproc = Proc::Simple->new(); 
$SIG{INT} = \&tsktsk;
$SIG{KILL} = \&tsktsk;
sub tsktsk {
	warn "kill jobs";
  	$myproc->kill() if $myproc->poll();
}

my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);
 my $ssh = Net::SSH::Perl->new("10.1.2.9");
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
	my $exit = run_cmd($dragen_cmd);
	unlink $file_lock if -e $file_lock;
	die() if $exit != 0;
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
			warn "lock";
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
	$myproc->start(qq{ssh $username\@10.1.2.9 $cmd;rm $file_lock});
	while ($myproc->poll()){
		sleep 2;
	}
	
	if ($myproc->exit_status() == 0){
		return 0;
	}
	else {
		#print $stdout."\n".$stderr."\n";
		return $myproc->exit_status;
	}	
	
	
}
