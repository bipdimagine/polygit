#!/usr/bin/perl
use CGI qw/:standard :html3/;
use FindBin qw($Bin);
use strict;
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
my $filein;
my $dir;
my $file_bed;
my $project_name;
my $chr_name;
my $fork = 1;
use List::MoreUtils qw(part);
use JSON;
use check_utils;


 
my $log_file;
my $end_ext = "uni";
my $patient_name;
my $vcf_file;
GetOptions(
	'project=s'	=> \$project_name,
	"log=s" 	=> \$log_file,
	"patient=s"	=> \$patient_name,
	"vcf=s"		=> \$vcf_file,
);
my $date = `date`;
chomp($date);
if ($log_file){ open (STDOUT,">>".$log_file); }
print  colored ['black ON_BRIGHT_MAGENTA'],"======= Check DBSNP ==== "; 
print  color 'reset';
print  "\n";

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $build = $project->buffer->build();

my $javac 		  = "/opt/java/latest/bin/java";
my $dbsnp 		  = "/data-xfs/public-data/$build/snp/gatk/latest/dbsnp.vcf.gz";
my $genomes1k	  = "/data-xfs/public-data/$build/snp/gatk/latest/1000genomes.vcf.gz";
my $hapmap		  = "/data-xfs/public-data/$build/snp/gatk/latest/hapmap.vcf.gz";
my $reference 	  = "/data-xfs/public-data/$build/genome/fasta/all.fa";
my $gatk 		  = "/bip-d/soft/distrib/GATK/gatk-latest/GenomeAnalysisTK.jar";
my $cmd_select 	  = "/opt/java/latest/bin/java  -jar /bip-d/soft/distrib/GATK/gatk-latest/GenomeAnalysisTK.jar -T SelectVariants   -R /data-xfs/public-data/$build/genome/fasta/all.fa -env -ef";
my $cmd_mendelian = "/opt/java/latest/bin/java  -jar /bip-d/soft/distrib/GATK/gatk-latest/GenomeAnalysisTK.jar -T SelectVariants -R /data-xfs/public-data/$build/genome/fasta/all.fa -mvq 50 --mendelianViolation";

my $error = 0;
my @lHash;
my $jsonFile = $project->getProjectPath().'../'.$project_name.'.resume.json';
warn $jsonFile;
foreach my $p (@{$project->get_list_patients($patient_name)}){
	warn $p->name();
	my ($pourcent,$nb) = check_dbsnp($p->name);
	my $text = $p->name."\tSnp: $nb\tdbSNP: $pourcent%";
	$error = 1 if  $pourcent <= 95;
	my $color = 'black ON_BRIGHT_GREEN';
	if ($pourcent >= 90 && $pourcent < 95) { $color = 'black ON_BRIGHT_YELLOW'; }
	elsif  ($pourcent < 90){ $color = 'black ON_BRIGHT_RED'; }
	print  colored [$color],$text  ;
	print  color 'reset';
	print  "\n";
	my $hash;
	$hash->{'patient'} = $p->name;
	$hash->{'nb_snp'} = $nb;
	$hash->{'perc_dbSNP'} = $pourcent.'%';
	push(@lHash, $hash);
}

my $jsonRes = check_utils::createJson(\@lHash, $jsonFile);
open(JSON_FILE, ">$jsonFile");
print JSON_FILE $jsonRes;
close(JSON_FILE);
#exit(1) if $sex_error > 0;
exit(0);


##### METHODS #####


sub check_dbsnp {
	my ($name) = @_;
	my $line_header = `grep "#CHR" $vcf_file`;
	chomp($line_header);
	my @theader =  split(" ",$line_header);
	my $index;
	for (my $i=8;$i< @theader;$i++){
		$index = $i if $theader[$i] eq $name;
	}
	unless ($index){
		my $string;
		for (my $i=8;$i< @theader;$i++){
			$string.=$theader[$i].";"; 
		}
		print " ERROR NIT FOUND $name in vcf\n";
		print $string."\n";
		warn $string;
	}
	confess() unless $index;
	$index ++;
	
	
	my $cmd1 = qq{ grep -v "#" $vcf_file | cut -f 3,$index |grep -v "0/0" | grep -v "\\./\\." | };
	
	#my $cmd1 = qq{/bip-d/soft/distrib/vcftools/latest/bin/vcf-subset $vcf_file -c $name -e | grep -v "#" |}; 
	open ("TOTO", $cmd1 );
	my $nb = 0;
	my $rs = 0;
	my $last_chr;
	while (my $line = <TOTO>){
		$nb++;
		$rs ++ if $line =~ /rs/;
		#last if $nb > 8000;
		#last if $line =~ /chr2/;
	}
#	warn $cmd1;
#	my $cmd2 = qq{/bip-d/soft/distrib/vcftools/latest/bin/vcf-subset $vcf_file -c $name -e | grep -v "#" | grep rs | wc -l }; 
#	my ($nb) = `$cmd1`;
#	my ($rs) = `$cmd2`;
#	chomp($rs);
#	chomp($nb);
close (TOTO);
my $pourcent =0;
	 $pourcent = int( ($rs/$nb) *10000) if $nb >0; 
	return ($pourcent/100,$nb);
	return 2;
}
