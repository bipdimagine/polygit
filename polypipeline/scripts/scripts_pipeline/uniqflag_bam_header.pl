#!/usr/bin/perl
use FindBin qw($Bin);

use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
#use Set::IntSpan;
use GBuffer;
use strict;
use Data::Dumper;
use threads;
use Thread::Queue;
use Getopt::Long;
use File::Temp;
use Set::IntSpan::Fast::XS;
my $bam_file;
my $file_primer;
 use Bio::DB::Sam;
 my $project_name;
 my $bam_file_out;
 
GetOptions(
	'file=s'  	=> \$bam_file,
	'file_out=s'  	=> \$bam_file_out,
	'project=s'  	=> \$project_name,
);
my $buffer = GBuffer->new();

my $samtools = $buffer->software("samtools"); 
my $sambamba = $buffer->software("sambamba"); 

my $project = $buffer->newProject( -name => $project_name );


 my $sam = Bio::DB::Sam->new(-bam  =>$bam_file);
 
 my $header = $sam->header();
 my @text =  split("\n",$header->text);
 my @new;
 while (<STDIN>){
 	chomp();
 	push(@text,$_);
 }
 my %dejavupg;
 my @headers ;
 my $fork =16;

  open (SAM,"| $sambamba view -S --format=bam  /dev/stdin -t $fork >$bam_file_out");
  
 foreach my $line (@text){
 	next unless $line =~/^@/;
 	unless  ($line =~/^\@PG/){
 		push(@headers,$line);
 		#print $line."\n";
 		next;
 	}
 	
	$line =~/ID\:(\w+)/;
	my $pg = $1;
	my $lcpg = lc($pg);
	
	 if (exists $dejavupg{$lcpg}){
	 	$dejavupg{lc($pg)} ++;
	 	my $newpg = $pg.$dejavupg{$lcpg};
	 	 $line =~ s/ID:$pg/ID:$newpg/;
	 }
	$dejavupg{$lcpg} ++;
	push(@headers,$line);
	#print $line."\n";
 }
  my $header_file = $bam_file.".header.sam";
 open (HEADER,">$header_file") || die("problem");
 
 
print HEADER join("\n",@headers);
close(HEADER);
system("$samtools reheader $header_file $bam_file > $bam_file_out && $samtools index $bam_file_out;rm $header_file");


 

