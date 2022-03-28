#!/usr/bin/perl
use CGI qw/:standard :html3/;
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";

use JSON;
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
my $filein;
my $dir;
my $file_bed;
my $project_name;
my $chr_name;
my $fork = 1;
use List::MoreUtils qw(part);
use check_utils;


my $log_file;
my $end_ext = "uni";
my $vcf_file;
GetOptions(
	'project=s'   => \$project_name,
	"log=s" =>\$log_file,
);
my $date = `date`;
chomp($date);
#$log_file = "toto.log";
if ($log_file){
open (STDOUT,">>".$log_file);
}

print  colored ['black ON_BRIGHT_MAGENTA'],"======= Check Sample Sex ==== "; 
print  color 'reset';
print  "\n";

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $ped_file = $project->getPedigreeFile();
my $build = $project->buffer->build();
my $javac 		  = "/opt/java/latest/bin/java";
my $dbsnp 		  = "/data-xfs/public-data/$build/snp/gatk/latest/dbsnp.vcf.gz";
my $genomes1k	  = "/data-xfs/public-data/$build/snp/gatk/latest/1000genomes.vcf.gz";
my $hapmap		  = "/data-xfs/public-data/$build/snp/gatk/latest/hapmap.vcf.gz";
my $reference 	  = "/data-xfs/public-data/$build/genome/fasta/all.fa";
my $gatk 		  = "/bip-d/soft/distrib/GATK/gatk-latest/GenomeAnalysisTK.jar";
my $cmd_select 	  = "/opt/java/latest/bin/java  -jar /bip-d/soft/distrib/GATK/gatk-latest/GenomeAnalysisTK.jar -T SelectVariants   -R /data-xfs/public-data/$build/genome/fasta/all.fa -env -ef";
my $cmd_mendelian = "/opt/java/latest/bin/java  -jar /bip-d/soft/distrib/GATK/gatk-latest/GenomeAnalysisTK.jar -T SelectVariants -R /data-xfs/public-data/$build/genome/fasta/all.fa -mvq 50 --mendelianViolation";

my $hfam =  $project->families();
my $sex_error = 0;
my @lHash;
foreach my $p (@{$project->getPatients}){
	my $bam = $p->getBamFile();
	my ($sex, $depth) = check_sex($bam);
	$depth = 0 unless ($depth);
	print  colored ['black ON_BRIGHT_GREEN'],$p->name()."\tpedigree: ".$p->sex."\tcompute: $sex OK\tdp: $depth" if $sex eq $p->sex ;
	print  colored ['black ON_BRIGHT_RED'],$p->name()."\t	pedigree: ".$p->sex."\tcompute: $sex  FALSE\tdp: $depth" if $sex ne $p->sex ;
	$sex_error = 1 if $sex ne $p->sex ;
	print  color 'reset';
	print  "\n";
	my $hash;
	$hash->{'patient'} = $p->name();
	$hash->{'ped_sex'} = $p->sex();
	$hash->{'bipd_status_sex'} = 'OK' if $sex eq $p->sex;
	$hash->{'bipd_status_sex'} = 'PROBLEM' if $sex ne $p->sex;
	$hash->{'bipd_sex_depth'} = $depth;
	push(@lHash, $hash);
}
if ($log_file) { close(STDOUT); }
my $jsonFile = $project->getProjectPath().'../'.$project_name.'.resume.json';
my $jsonRes = check_utils::createJson(\@lHash, $jsonFile);
open(JSON_FILE, ">$jsonFile");
print JSON_FILE $jsonRes;
close(JSON_FILE);
exit(1) if $sex_error > 0;
exit(0);


##### METHODS #####


sub check_sex {
	my ($bam) = @_;
	my @apos = ("chrY:6736371-6736371","chrY:6736443-6736443","chrY:6737845-6737845");
	my $cmd = qq{samtools depth $bam };
	my $maxDp = 0;
	foreach my $pos (@apos){
		my $cmd2 = $cmd." -r $pos | cut -f 3 ";
		my ($depth) = `$cmd2`;
		chomp($depth);
		return (1, $depth) if $depth > 10;
		$maxDp = $depth if ($depth > $maxDp);
	}
	return (2, $maxDp);
}
