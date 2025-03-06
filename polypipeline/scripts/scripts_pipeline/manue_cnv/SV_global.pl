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
my $nodejavu;
GetOptions(
	'project=s' => \$project_name,
	'debug=s' => \$debug,
	'fork=s' => \$fork,
	'nodejavu=s' => \$nodejavu,
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
my $log_dir = $project->getCNVDir();

my $load = qq{module load tcltk/8.6.9};
my $annotsv = qq{/software/distrib/AnnotSV_2.0/bin/AnnotSV};

#warn "SV ANNOT";
#launch annotsv  step 1
my $l = time;
my $lfile = "$log_dir/$l"."_sv.log";
my $cmd = qq{$Bin/annotsv.pl -project=$project_name -fork=$fork};

system("$cmd && date > $lfile");
#die() unless -e $lfile;


#my $prg = qq{/software/polyweb/poly-disk/poly-src/polymorphism-cgi/manta/Set_allSV_byPatient.pl };

#my $prg = qq{/data-xfs/dev/eollivie/polymorphism-cgi/manta/Set_allSV_byPatient.pl};

my $prg = qq{$Bin/Set_allSV_byPatient.pl };
$lfile = "$log_dir/$l"."_allsv.log";

$cmd = qq{$prg -fork=$fork -project=}.$project->name();
system($cmd." &&  date > $lfile");
#die() unless -e $lfile;

 #$prg = qq{$Bin/bestone.pl -project=$project_name -fork=$fork};
 #$lfile = "$log_dir/$l"."_svPatients.log";
#$cmd = qq{$prg projectname=}.$project->name();
#system($cmd." &&  date > $lfile");
#die() unless -e $lfile;


my $prg_gather_cnv = qq{$Bin/gather_cnv.pl };
$cmd = qq{$prg_gather_cnv -fork=$fork -project=}.$project->name();
system($cmd." &&  date > $lfile");

my $prg_gather_cnv_project = qq{$Bin/gather_cnv_project.pl };
$cmd = qq{$prg_gather_cnv_project -fork=$fork -project=}.$project->name();
system($cmd." &&  date > $lfile");

#TODO: besoin ici ??

#my $prg4 = qq{$Bin/parseBndManta.pl project=$project_name};
# $lfile = "$log_dir/$l"."_eq.log";
#system($prg4." &&  date > $lfile" );
#die() unless -e $lfile;
#
#unless ($nodejavu){
#my $prg3 = qq{$Bin/CNV_Set_DejaVu.pl };
# $lfile = "$log_dir/$l"."_dv.log";
#system($prg3." &&  date > $lfile" );
#die() unless -e $lfile;
#
#
#my $prg5 = qq{$Bin/SVeq_Set_DejaVu.pl };
# $lfile = "$log_dir/$l"."_dveq.log";
#system($prg5." &&  date > $lfile" );
#die() unless -e $lfile;
#}

#system("date > $file_done");
exit(0);



