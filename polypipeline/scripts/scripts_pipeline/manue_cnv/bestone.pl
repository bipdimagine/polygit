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
my $fork = 10;
my $project_name;
my $debug;
my $patient_name;
GetOptions(
	'project=s' => \$project_name,
	'debug=s' => \$debug,
	'fork=s' => \$fork,
	#'patient=s' => \$patient_name
);


my $buffer = new GBuffer;
my $project = $buffer->newProject( -name 			=> $project_name );
my $file_done = $project->getCNVDir()."/".$project->name.".done";
unlink $file_done if -e $file_done;
my $dir;
$dir->{manta}->{dir}= $project->getVariationsDir("manta");
$dir->{canvas}->{dir}= $project->getVariationsDir("canvas");
$dir->{wisecondor}->{dir}= $project->getVariationsDir("wisecondor");
my $load = qq{module load tcltk/8.6.9};
my $annotsv = qq{/software/distrib/AnnotSV_2.0/bin/AnnotSV};
my $prg2 = qq{/software/polyweb/poly-disk/www//cgi-bin/polymorphism-cgi/manta/Retrieve_allSV_Patient.pl };
my $pm = new Parallel::ForkManager($fork);
my $hjobs;
my $job_id = time;
$pm->run_on_finish(
    	sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    		my $j = $data->{job};
    		delete $hjobs->{$j};
    		
    }
    );	


foreach my $p (@{$project->getPatients}){
	$job_id ++;
	$hjobs->{$job_id} ++;
	 my $pid = $pm->start and next;
	my $fileout = $p->getBestOneFileName();
	unlink $fileout if -e $fileout;
	my $pname = $p->name();
	warn $pname." ".$fileout;
	my $cmd = $prg2." projectname=$project_name filename=$pname minlength=10000 maxlength=2000000000 transmission=all maxfreq=nomax dejavu=10 genotype=both select_best=1 chrom=all cytoband=all genes=all print=1 >".$fileout;

#	warn $cmd;
#	die();
	warn $cmd;
	system($cmd);
	$pm->finish(0,{job=>$job_id});
}
$pm->wait_all_children();
confess() if scalar(keys %{$hjobs});
exit(0);