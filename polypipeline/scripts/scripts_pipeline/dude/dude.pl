#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
#use Set::IntSpan;
use GenBoNoSqlLmdb;

use Carp;
use strict;
use Set::IntSpan::Fast::XS ;
use Data::Dumper;
use GBuffer;
use Proc::Simple;
use Getopt::Long;
my $myproc = Proc::Simple->new(); 
 $myproc->redirect_output ("/dev/null", "/dev/null");
my $fork;
my $project_name;
my $debug;
GetOptions(
	'project=s' => \$project_name,
	'fork=s' => \$fork,
	'debug=s' => \$debug,
);

die("man !!!  no fork ") unless $fork; 
die("dude !!!  no project ") unless $project_name; 

my $buffer = new GBuffer;
my @list_cmd = (11,"12_controls","12cai",13,"14");
my $project = $buffer->newProject( -name 			=> $project_name );
#exit(0) if $project->isDude();
foreach my $patient (@{$project->getPatients}){
		my $dir_out= $project->getVariationsDir("dude");
		my $filegz = $dir_out."/".$patient->name.".dude.lid.gz";
		if (-e $filegz){
			unlink $filegz;
			unlink $filegz.".tbi" if -e $filegz.".tbi";
		}
}

foreach my $l (@list_cmd){
	my $time = time;
	warn " START => $l  ************* ";
	if ($debug){
		system ("$Bin/$l.pl -project=$project_name -fork=$fork") == 0  or die("");
	}
	else {
		$myproc->start("$Bin/$l.pl",             # Launch an external program
	                  "-project=$project_name","-fork=$fork"); 
		my $exit_status = $myproc->wait(); 
		if ($exit_status == 0){
			warn "end $l ==> $exit_status ".abs(time-$time);
	                  			
		}
		else {
			warn "ERROR on  $Bin/$l.pl  ==> $exit_status ";
			die();
		}
	}

}


foreach my $patient (@{$project->getPatients}){
		my $dir_out= $project->getVariationsDir("dude");
		my $filegz = $dir_out."/".$patient->name.".dude.lid.gz";
		unless (-e $filegz){
			die("you have a problem somewhere: $filegz ".$patient->name);
		}
}

print "OK \n";
exit(0);