#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use File::Basename;
use  File::Temp;
use YAML::Syck;
use colored; 
use Text::Table;
use POSIX qw{strftime};
use IO::Prompt;


 #my $yaml = YAML::Tiny->read();
my $file = $ARGV[0];

die("problem parsing results ") unless -e $file;
open(O,$file);

my @part_yaml;
my $stringy="";
my $part = 0;
my %desc;
while (<O>){
	my $s = $_;
	if ($s =~ /# Program overview/){
		$desc{overview} = $part;
	}
	if ($s =~ /# Task details/){
		$desc{tasks} = $part;
	}
	$stringy.=$s;
	if ($s =~ /^---/){
		push(@part_yaml,$stringy);
		$stringy ="";
		$part ++;
	
	}
}

close (O);
 my  $yaml = Load($part_yaml[$desc{overview}]);
print colored ['black ON_BRIGHT_YELLOW'], "General Report" ;
print "\n";
print colored ['black ON_BRIGHT_YELLOW'], "Task executed : ". $yaml->{tasksExecuted};
print "\n";
 print colored ['black ON_BRIGHT_GREEN'], "Task Failed : ".$yaml->{tasksFailed} if  $yaml->{tasksFailed} ==0;
print colored ['black ON_BRIGHT_RED'], "Task Failed : ".$yaml->{tasksFailed} if  $yaml->{tasksFailed} > 0 ;
print "\n";
print colored ['black ON_BRIGHT_YELLOW'], "Running Time : ". $yaml->{runtime};
print "\n";

$yaml = undef;
 my  $yaml_task = Load($part_yaml[$desc{tasks}]);

my $task_type;
my $task_patient;
my %error_patient;
foreach my $task (@{$yaml_task->{tasks}}){
	my @name = split(".",$task->{name});
	my ($patient,$type) = split("@",$task->{name});
	my $status = $task->{ok};
	$status = 1 if $task->{ok} eq "true";
	my @id = split(/\./,$task->{id});
	
	my ($t,$nb) = split("_",$id[-1]);
	$task_type->{$type} += $nb unless exists $task_type->{$type};
	$task_patient->{$patient}->{$type}->{status} =$status;
	$task_patient->{$patient}->{$type}->{task} =$task;
	$task_patient->{$patient}->{$type}->{elapsed} =$task->{elapsed};
}
$yaml_task = undef;
my @atasks = sort{$task_type->{$a} <=> $task_type->{$b}} keys %$task_type; 
my @header =("patients",@atasks,"total");
my @lines;
my %sum_task;

	my %nb_running = 0;
	my $failed;
foreach my $patient (sort {$a cmp $b}   keys %$task_patient){
		my @line;
		push(@line,$patient);
		my $sum =0;
	
		foreach my $type (@atasks){
				my $color = "green";
				my $text = $task_patient->{$patient}->{$type}->{elapsed};
				my @t = split(":",$text);
				$sum += $t[0]*3600+$t[1]*60+$t[2];
				$nb_running{$type} ++;
				unless (exists $task_patient->{$patient}->{$type}){
							$color = "green";	
							$text ="-";
				}
				elsif (! (exists $task_patient->{$patient}->{$type}->{status})){
						$error_patient{$patient}++;
							$color = "white";	
							$text ="-";
							$nb_running{$type} --;
				}
				elsif ($task_patient->{$patient}->{$type}->{status} ne 1){
						$error_patient{$patient}++;
							$color = "red";	
							$text ="FAILED";
							$failed =1;
				}
				
					push(@line, colored::stabilo("$color",$text,1));
					$sum_task{$type} += $t[0]*3600+$t[1]*60+$t[2] ;
		}
		$sum_task{total} += $sum;
		push(@line, colored::stabilo("cyan",strftime("\%H:\%M:\%S", gmtime($sum)),1));
		
		push(@lines,\@line);
}
my @line = ("total");
my $nbp = scalar(keys %$task_patient);
foreach my $type (@atasks){
	 $nb_running{$type} =1 if $nb_running{$type} ==0;
	
	my $mean = int($sum_task{$type}/$nb_running{$type});
	push(@line, colored::stabilo("magenta",strftime("\%H:\%M:\%S", gmtime($mean)),1));
}
push(@line, colored::stabilo("magenta",strftime("\%H:\%M:\%S", gmtime($sum_task{total})),1));
push(@lines,\@line);

	my $tb = Text::Table->new(
	@header,
    );
    	$tb->load(@lines);
	print $tb;
	print "\n";
	#exit(0) unless (keys %error_patient);
	close(STDIN);
#my $choice = prompt("run this/these step(s)   (y/n) ? ");
#	exit(0) if ($choice ne "y"); 
	exit(0) unless $failed;
	print colored::stabilo("red","ERROR LOG ",1);
	#print "\n";
	#exit(0);
	foreach my $patient (keys %$task_patient){
		 next unless exists $error_patient{$patient};
			print colored::stabilo("cyan","$patient ",1);
			print "\n";
		foreach my $type (@atasks){
				if ($task_patient->{$patient}->{$type}->{status} ne 1){
					if ($task_patient->{$patient}->{$type}->{task}->{stderr}){
						print colored::stabilo("magenta",$type.' '. $task_patient->{$patient}->{$type}->{fileout}." \n".$task_patient->{$patient}->{$type}->{filein},1)  ;
						print "\n";
						print  $task_patient->{$patient}->{$type}->{outFiles};
						print "\n";
						print colored::stabilo("magenta","----------------------STDERR-------------------" ,1)  ;
						print "\n";
						
						print  $task_patient->{$patient}->{$type}->{fileout}." \n".$task_patient->{$patient}->{$type}->{filein}."\n";
						print $task_patient->{$patient}->{$type}->{task}->{stderr}."\n";
						print colored::stabilo("magenta","----------------------STDOUT-------------------" ,1)  ;
						print "\n";
						print $task_patient->{$patient}->{$type}->{task}->{stdout}."\n";
						print "\n";
							print $task_patient->{$patient}->{$type}->{task}->{outFiles}."\n";
							print $task_patient->{$patient}->{$type}->{task}->{filein}."\n";
						print colored::stabilo("magenta","----------------------------------------------" ,1)  ;
							
					}
				}
			
		}
	}
