#!/usr/bin/perl
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



my $log_file;
my $end_ext = "uni";
my $vcf_file;
GetOptions(
	'project=s'   => \$project_name,
	"vcf=s"	=> \$vcf_file,
	"log=s" =>\$log_file,
);

exit(1) unless -e $vcf_file;
my $date = `date`;
chomp($date);
#$log_file = "toto.log";
if ($log_file){
open (STDOUT,">>".$log_file);
}

print  colored ['black ON_BRIGHT_MAGENTA'],"======= test pedigree ==== "; 
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
my $cmd_select 	  = "/opt/java/latest/bin/java  -jar /bip-d/soft/distrib/GATK/gatk-latest/GenomeAnalysisTK.jar -T SelectVariants -R /data-xfs/public-data/$build/genome/fasta/all.fa -env -ef";
my $cmd_mendelian = "/opt/java/latest/bin/java  -jar /bip-d/soft/distrib/GATK/gatk-latest/GenomeAnalysisTK.jar -T SelectVariants -R /data-xfs/public-data/$build/genome/fasta/all.fa -mvq 50 --mendelianViolation";

my $hfam =  $project->families();
my $mendel = 0;
foreach my $fam (keys %{$hfam}){
	next unless exists $hfam->{$fam}->{parents};
	next unless exists $hfam->{$fam}->{children};
	print  colored ['black ON_BRIGHT_CYAN'],"======= Family $fam ==== "; 
	print  color 'reset';
	print  "\n";
	foreach my $p (@{$hfam->{$fam}->{children}}){
		my $mother = $hfam->{$fam}->{mother};
		my $father = $hfam->{$fam}->{father};
		$mendel += test_mendelian($p,$mother,$father);
	}
}
exit(1) if $mendel > 0;
exit(0);

#--variant res_error.vcf -o res_error.mendel.vcf 
sub test_mendelian {
	my ($p,$m,$f) =@_;
	my $dir_out = $project->getPipelineCallingDir("unifiedgenotyper");
	my $out_vcf1 = $dir_out."/".$p->name.".trio.vcf";
	my $mendel_vcf = $dir_out."/".$p->name.".mendel.vcf";
	my $cmd1 = $cmd_select." --variant $vcf_file -o $out_vcf1 ";
	$cmd1 .= "-sn ".$p->name." -sn ".$m->name()." -sn ".$f->name;

	system($cmd1." >/dev/null");
	my $ped_file = $dir_out."/".$p->name.".trio.ped";
	open (PED,">".$ped_file);
	print PED 	$m->{pedigree_line}."\n".$f->{pedigree_line}."\n".$p->{pedigree_line};
	close(PED);	
	
	my $cmd2 = $cmd_mendelian." --variant $out_vcf1 -o $mendel_vcf -ped $ped_file >/dev/null";
	
	system($cmd2);
	my ($nb_mendelian) = `cat $mendel_vcf | grep -v "#" | wc -l`;
	my ($nb_total) = `cat $out_vcf1 | grep -v "#" | wc -l`;
	
	chomp($nb_total);
	chomp($nb_mendelian);
	
	my $pourcent = int($nb_mendelian/$nb_total *10000);
	$pourcent /= 100;
	my $all_name = $p->name(); 
	warn Dumper ($m);
	$all_name .= "::".$m->name if exists $m->{name};
	$all_name .= "::-" unless exists $m->{name};
	$all_name .= "::".$p->name if exists $p->{name};
	$all_name .= "::-" unless exists $p->{name};
	print  colored ['black ON_BRIGHT_GREEN'],$all_name." ==> ".$pourcent if $pourcent < 1;
	print  colored ['black ON_BRIGHT_YELLOW'],$all_name." ==> ".$pourcent."" if $pourcent <=3 && $pourcent >=1;
	print  colored ['black ON_BRIGHT_RED'],$all_name." ==> ".$pourcent."" if $pourcent >=3;
	print  color 'reset';
	print  "\n";
	return 1 if $pourcent >3;
	return 0;
}


