#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use Date::Tiny;
use POSIX;
 
use PBS::Client;
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moose;
use MooseX::Method::Signatures;
use Class::Inspector;
use Digest::MD5::File ;
use GBuffer;
use GenBoProject;
use colored; 
use Text::Table;

my $projectName;
my $list_patients;
my $log_file;

GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$list_patients,
	'log=s' => \$log_file,
	
);

if ($log_file){
	open (STDOUT,">>".$log_file);
}

my $buffer = GBuffer->new();

unless ($projectName) {
	usage();
}
my $project = $buffer->newProject( -name => $projectName );

my $patients = $project->get_list_patients($list_patients);

my $dir_prod  = $project->getGvcfDir("haplotypecaller");

my $bams;
my @lines;
my @header =("name","bam","gvcf");
my $t1 = localtime();
		
          
	

foreach my $p (@$patients){
	my @line;
	my $bam_file = $p->getBamFileName();
	push(@line,$p->name);
	if (-e $bam_file){
		my $date = POSIX::strftime( 
             "%d/%m/%y", 
             localtime( 
                 ( stat $bam_file )[8]
                 )
             );
             	my $t = (stat $bam_file )[8];
             	my $days_difference = int((time - $t) / 86400);
             	my $color = "green";
             	if ($days_difference > 0){
             		$color = "cyan";
             	}
             	if ($days_difference > 1){
             		$color = "orange";
             	}
             	push(@line, colored::stabilo("$color","+$days_difference day $date" ,1));
		
	}
	else {
			push(@line,colored::stabilo("red","XXX",1));
	}
	
	my $gvcf = $dir_prod."/".$p->name.".g.vcf.gz";
		my $gvcftbi = $dir_prod."/".$p->name.".g.vcf.gz.tbi";
	if (-e $gvcf && -e $gvcftbi){
	my $date2 = POSIX::strftime( 
             "%d/%m/%y", 
             localtime( 
                 ( stat $gvcf )[9]
                 )
             );
             my $t = (stat $gvcf )[9];
             	my $days_difference = int((time - $t) / 86400);
             	my $color = "green";
             		if ($days_difference > 0){
             		$color = "cyan";
             	}
             	if ($days_difference > 1){
             		$color = "orange";
             	}
             	push(@line, colored::stabilo("$color","+$days_difference $date2 day" ,1));
				#push(@line, colored::stabilo("green","$date2",1));
	
	}
	else {
		push(@line,colored::stabilo("red","XXX",1));
	}
	push(@lines,\@line);
}

	
	my $tb = Text::Table->new(
	@header,
    );
    	$tb->load(@lines);
	print $tb;
	print "\n";
	my $now = Date::Tiny->now;
	print $now->day."/".$now->month."/".$now->year."\n";
