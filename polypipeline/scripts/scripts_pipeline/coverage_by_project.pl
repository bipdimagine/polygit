#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
 use autodie;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moose;
use MooseX::Method::Signatures;
use Try::Tiny;
use File::Temp qw/ tempfile tempdir /;
use colored;
use String::ProgressBar;
use Text::Table;

my $projectName;
my $filename;

my $name;
my $patients_name;
my $steps_name;
my $fork;
my $force;

GetOptions(
	'project=s' => \$projectName,

	'patients=s' => \$patients_name,
	'fork=s' =>\$fork,
	'force=s'=>\$force,
);
$fork =2 unless $fork;
my $buffer = GBuffer->new();

my $hsex = {
				2=>qq{F},
				1=>qq{M} ,
				'-1'=>qq{U} 
};
die("patients=(all|(toto,titi))") unless $patients_name;
my $project = $buffer->newProject( -name => $projectName );

my $patients = $project->get_list_patients($patients_name );
confess("no patients") unless scalar(@$patients);
my $errors;


	my $pr = String::ProgressBar->new( max => scalar(@$patients) );
colored::stabilo("yellow","START COVERAGE");
my $cc= 0;
my $msg;
my $tmp = "/poly-disk/tmp/";
$tmp ="/tmp" unless -e $tmp;
foreach my $p (@$patients){
 	$pr->write();
	$cc++;
     	$pr->update($cc);
	
	my $tmp = File::Temp->new( TEMPLATE => 'tempXXXXX',
                        DIR => $tmp,
                        SUFFIX => '.log');
  	my $ftmp = $tmp->filename;

	
	my $filein = 	$p->getBamFile(undef,1);
	unless ($filein){
		$errors->{$p->name()} ="No Bam File";
		next;
	}
	my $fcoverage = $p->getCoverageFile();
	my $tbi = $fcoverage.".tbi"; 
	if (-e $fcoverage && -e $tbi && !($force)){
		
		$errors->{$p->name()} ="Previously Done";
		next;
	}
	unlink $tbi if -e $tbi;
	unlink $fcoverage if -e $fcoverage;
	my $coverage_dir = $project->getAlignmentRootDir."/coverage/";

	mkdir $coverage_dir unless -e $coverage_dir;
	my $name= $p->name();

	my $fileout = $coverage_dir . "/" . $name.".cov.gz";
	my $bed     = $p->getCaptureFile();
	die("not found bed file : $bed") unless -e $bed;
	my $cmd = qq{perl $Bin/coverage.pl -filein=$filein -dir=$coverage_dir -bed=$bed -fork=$fork -patient=$name -project=$projectName 2>$ftmp 1>&2 };
	warn $cmd;
	eval {

	my $return = system("$cmd ");
	my $cov = $p->coverage();
	print "------------".colored::stabilo("blue", $p->name()." ".$cov->{mean}." ".$cov->{'30x'},1)." ". $hsex->{$p->compute_sex()} ."---------------\n";
	#warn $p->name()." ".$cov->{mean}." ".$cov->{'30x'};
	if ($return ne 0){
		$errors->{$p->name()} = "ERROR";
		
		 my @t = `cat $ftmp`;
		push(@t,"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n"."Problem with patient ".$p->name."\n"."&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
		warn @t;
		warn $cmd;
		
	}
	
	
	};
	
	
}

my @lines;
foreach my $p (@$patients){
	my @cols;
	push(@cols,colored::print_color("cyan",$p->name(),1));
	my $sex_eval = $p->compute_sex(); 
	my $sex = $p->sex(); 
	my $class = "green";
	
	if ($sex_eval ne $sex && $sex_eval ne -1){
	  push(@cols,colored::stabilo("red",$hsex->{$sex}."/".$hsex->{$sex_eval},1));
	} 
	else {
	push(@cols,colored::stabilo("green",$hsex->{$sex}."/".$hsex->{$sex_eval},1));
	}
	my $f1 = $p->getVariationsFile();
	my $nbv = "NC";
	$nbv = `zgrep -vc "#" $f1` if -e $f1;
	chomp($nbv);
	my $f2 = $p->getIndelsFile();
	my $nbi ="NC";
	$nbi = `zgrep -vc "#" $f2` if -e $f2;
	chomp($nbi);
	push(@cols,colored::print_color("blue",$nbv,1));
	push(@cols,colored::print_color("blue",$nbi,1));
	 if (exists $errors->{$p->name()}){
		push(@cols,colored::stabilo("red",$errors->{$p->name()},1));
push(@lines,\@cols);
next;

	}

	push(@cols,colored::stabilo("green","OK",1));
	my $cov = $p->coverage();
	push(@cols,colored::print_color("green",$cov->{mean},1));
	push(@cols,colored::print_color("green",$cov->{'30x'},1));

	push(@lines,\@cols);
}
my @header =("Patient","Sex/SRY","Variations","indels","Run_coverage","Mean","30x");
	my @line_null =("------","------","------","------","------","------","------","------","------","------","------");
	
	my $tb = Text::Table->new(
      @header
    );
    

	$tb->load(@lines);
	print $tb;
	print "\n";
	print "\t\t\t\t-*-*".colored::print_color("green"," Project " .$projectName,1)."-*-*-\n";


