#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/";
use lib "$Bin/../../../GenBo/lib/GenBoDB";
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../../../GenBo/lib/kyoto/";
use lib "$Bin/../../../GenBo/lib/GenBoDB/writeDB";
use lib "$Bin/../../packages";
use Logfile::Rotate;
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moose;
use GBuffer;
use GenBoProject;
use colored; 
use Config::Std;
use Text::Table;
use file_util;
use File::Temp qw/ tempfile tempdir /;; 

use Term::Menus;
 use Proc::Simple;
 use Storable;
use JSON::XS;
use Net::SSH::Perl; 

 
my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);
 my $ssh = Net::SSH::Perl->new("10.1.2.9");
$ssh->login("$username");
my $project_name;
my $patient_name;
GetOptions(
	'project=s' => \$project_name,
	'patients=s' => \$patient_name,
);
my $buffer = GBuffer->new();
my $project = $buffer->newProjectCache( -name 			=> $project_name );
my $patient = $project->getPatient($patient_name);
my $run = $patient->getRun();
my $hps =  $run->getAllPatientsInfos();
my %contr_projects;

map {$contr_projects{$_->{project}} ++} @$hps;
my $dir_out = $patient->getDragenDir("PON");
my @files;
foreach my $pr (keys %contr_projects){
	my $buffer1 = new GBuffer;
		my $project1 = $buffer1->newProjectCache( -name 			=> $pr );
		#next() if $project1->name() eq "NGS2018_2286";
		foreach my $p (grep{$_->{project} eq $pr} @$hps){
			next if $p eq $patient_name;
			my $p1 = $project1->getPatient($p->{patient});	
			push(@files,$p1->targetGCFile) if -e $p1->targetGCFile;
		}
}

die(@files) if scalar(@files) < 5;
my $ref_dragen = $project->getGenomeIndex("dragen");
my $patient_gc = $patient->targetGCFile;
my $pon_list = "$dir_out/pon.txt";
my $f1= $dir_out."/".$patient_name.".cnv.vcf.gz";
open(PON,">$pon_list");
print PON join("\n",@files);
close(PON);
my $cmd = qq{dragen -f -r $ref_dragen --output-directory $dir_out --output-file-prefix $patient_name --enable-map-align false --enable-cnv true --cnv-input $patient_gc --cnv-normals-list $pon_list --cnv-enable-gcbias-correction false};
my $exit = 0;
 $exit = system(qq{$Bin/../run_dragen.pl -cmd="$cmd"}) unless -e $f1;
#warn qq{$Bin/../run_dragen.pl -cmd="$cmd"};
#my ($out,$err,$exit) = $ssh->cmd($cmd);
die($cmd) unless $exit == 0;


my $dir_prod = $project->getVariationsDir("dragen-pon");
system("rsync -rav $f1* $dir_prod 2>/dev/null");
warn "rsync -rav $f1 $dir_prod 2>/dev/null";
exit(0);

