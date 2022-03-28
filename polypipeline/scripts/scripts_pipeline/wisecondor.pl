#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
#use Set::IntSpan;
use GBuffer; 
use Data::Dumper;
use Getopt::Long;
use Carp;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use String::ProgressBar;
use List::Util qw(sum);



 my $project_name;
 my $fork;
 my $callable_intspan_file;
 my $patient_name;
 #my $low_calling;
 my $method;
 
 my $chr_syno ={
		1=> "chr1",
		2=>"chr2",
		3=>"chr3",
		4=>"chr4",
		5=>"chr5",
		6=>"chr6",
		7=>"chr7",
		8=>"chr8",
		9=>"chr9",
		10=>"chr10",
		11=>"chr11",
		12=>"chr12",
		13=>"chr13",
		14=>"chr14",
		15=>"chr15",
		16=>"chr16",
		17=>"chr17",
		18=>"chr18",
		19=>"chr19",
		20=>"chr20",
		21=>"chr21",
		22=>"chr22",
		X=>"chrX",
		Y=>"chrY",
		MT=>"chrM",
};

 
 my $fork = 5;
GetOptions(
	'project=s'   => \$project_name,
	"patient=s" => \$patient_name,
	"fork=s" => \$fork,
);
die("miss fork") unless $fork;


my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $wise  = $project->getSoftware('wisecondor');
my $patient = $project->getPatient($patient_name);
 my $bam = $patient->getBamFile() ;
 
 my $chr = $project->getChromosome("1");
  my $fileout = $patient->fileWiseCondor();
  #unless (-e $fileout) {
 if ($chr->fasta_name eq "1"){
 	#$bam = change_bam_file($bam);
 }
 
 my $cmd = "$wise convert $bam $fileout";
 system("$cmd");
 
 die() unless -e $fileout;

 
 
#
# 
# sub change_bam_file {
# 	my ($input_file) = @_;
# 	warn "change bam";
# 	my $dir = $project->getAlignmentPipelineDir("wise")."/$patient_name";
# 	system("mkdir -p $dir;chmod a+rwx $dir") unless -e $dir;
# 	my $samtools = $buffer->software("samtools");
# 	 	my $sambamba = $buffer->software("sambamba");
# 	
#	 my @chromosomes = `$samtools idxstats $input_file | cut -f 1`;
#	chomp(@chromosomes);
#	
#	
#	 my $pm2 = new Parallel::ForkManager($fork);
# my @tmp_file;
# 
# foreach my $chr (@chromosomes){
# 	warn $chr;
# 	next if $chr eq "*";
# 	next unless exists $chr_syno->{$chr};
#	my $bamout = $dir."/$patient_name.$chr.bam";
#
#	push(@tmp_file,$bamout);
#	next if -e $bamout;
#my $pid      = $pm2->start() and next;	
#open (BAM,"$samtools  view -h $input_file $chr|");
#my $pid = open(WRITEME, "| $samtools view -bS -   > $bamout") or die "Couldn't find file $input_file \n $!\n";
#while (<BAM>){
#	chomp();
#	my $line = $_;
#	my @t = split(" ",$line);
#	if ($line =~ /^\@/){
#		if ($line =~ /^\@SQ/){
#		
#			my ($sn,$ch) = split(":",$t[1]);
#			if (exists $chr_syno->{$ch}){
#				$t[1] = "SN:".$chr_syno->{$ch};
#			}
#			$line =  join("\t",@t);
#		
#		}
#	}
#	else {
#		
#		my $c1 = $t[2];
#		if (exists $chr_syno->{$c1}){
#			$t[2] = $chr_syno->{$c1};
#		}
#		else {
#			next;
#		}
#		my $c2 = $t[6];
#			if (exists $chr_syno->{$c2}){
#			$t[6] = $chr_syno->{$c2};
#		}
#		else {next;}
#	
#	}
#	#warn $line;
#	$line =  join("\t",@t)."\n";
#	print WRITEME $line;
#}
#close(BAM);
#close(WRITEME);
#
# 	$pm2->finish(0);
#	}
#$pm2->wait_all_children;
# my $split = join(" ",@tmp_file);
# my $final = "$dir/$patient_name.bam";
# warn "$sambamba  merge -t $fork $final $split";
# system("$sambamba  merge -t $fork $final $split");
#return $final;
#exit(0);
# }
