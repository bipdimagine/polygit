#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use GBuffer;
use Data::Dumper;
use Getopt::Long;
use Carp;
use Bio::DB::Sam;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use List::Util qw(sum);
my $filein;
my $dir;
my $file_bed;
 
my $dp_limit = 5;
my $al = 3;
my $project_name;
my $fork;
my $patient_name;
$| =1;
my $log_file;
my $vcf_final;
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"log=s" =>\$log_file,
	"patient=s" => \$patient_name,
	"out=s"=>\$vcf_final
);

$fork =8 unless $fork;

my $date = `date`;
chomp($date);
if ($log_file){
	system ("echo  ---- start running haplotypecaller $patient_name $date >> $log_file");
	open (STDOUT,">>".$log_file);
}
print "echo  ---- start running UnifiedGenotyper $patient_name $date\n";



#my $file_bed = "/data-xfs/public-data/HG19/capture/agilent/agilent.v50.bed";
my $samtools = "/bip-d/soft/bin/samtools";
my $bcftools = "/bip-d/soft/bin/bcftools";
my $vcfutil ="vcfutils.pl";
my $vcftools = qq{/bip-d/soft/distrib/vcftools/latest/bin/vcftools };

my $javac = qq{/opt/java/latest/bin/java  };
my $varscan =" $javac -jar /bip-d/soft/distrib/VarScan/VarScan.v2.3.6.jar ";
my $dbsnp = "/data-xfs/public-data/HG19/snp/gatk/latest/dbsnp.vcf.gz";
my $genomes1k	= "/data-xfs/public-data/HG19/snp/gatk/latest/1000genomes.vcf.gz";
my $hapmap		= "/data-xfs/public-data/HG19/snp/gatk/latest/hapmap.vcf.gz";

my $gatk = qq{/bip-d/soft/distrib/GATK/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar};
my $gatk = qq{/bip-d/soft/distrib/GATK/gatk-latest/GenomeAnalysisTK.jar};
my $gatk = qq{/bip-d/soft/distrib/GATK/gatk-latest.2.7/GenomeAnalysisTK.jar};
my $vcfannotate =qq{/bip-d/soft/distrib/vcftools/latest/bin/vcf-annotate}; 

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $chrs = $project->getChromosomes();
my $reference = $project->getGenomeFasta;
my $patient = $project->getPatient($patient_name);
my $coverage_dir = $project->getAlignmentRootDir."/coverage/";
my $coverage_file = $coverage_dir."/$patient_name.cov.gz";
my $dp;
my $param_mp;
my $param_gatk;
my $param_varscan = "";
my $chr_tr;
 if ($project->isDiagnostic){
 	my $trs = $project->getCapture->transcripts_name();
	foreach my $tr (@$trs){
		my $t = $project->newTranscript($tr);
		$chr_tr->{$t->getChromosome()->name} =  Set::IntSpan::Fast::XS->new() unless exists $chr_tr->{$t->getChromosome()->name};
		$chr_tr->{$t->getChromosome()->name}->add_range($t->start()-1000,$t->end+1000);

	}

 }

if (-e $coverage_file){
	my $dp = `tabix $coverage_file mean_all:99  | cut -f 3`;
	print  "tabix $coverage_file mean_all:99 \n";
	chomp($dp);
	print "coverage :" .$dp."\n";
	
	if ($dp < 50){
		
		$dp_limit = 3;
		$al =3;
	}
	else {
		
		$param_mp = "-L 1000000 -d 100000";
		$param_gatk = "-dcov 1000000";
		$dp_limit = 10;
		$al = 10;
		 $param_varscan = " --min-coverage 10 --min-reads2 10";
		
	}
}

mkdir ( $project->getCallingPipelineDir("unifiedgenotyper")) unless -e  $project->getCallingPipelineDir("unifiedgenotyper");
my $dir_out= $project->getCallingPipelineDir("unifiedgenotyper");
$vcf_final = $dir_out."/".$patient->name.".final.vcf" unless $vcf_final;
my $bamfile = $patient->getBamFile();


my $delete_bam;
my $csv2 =$patient->getRecalFile();
#if (-e $csv2){ 
#	my $fileori = $bamfile;
#	$bamfile =  $dir_out."/".$patient->name.".bam";
#	$delete_bam = $bamfile;
#	my $cmd2 =  "$javac -jar $gatk -I $fileori --out $bamfile -R $reference  -T PrintReads  -BQSR $csv2 -nct 8 -l off >/dev/null";
#	system($cmd2);
#}	







foreach my $thr (threads->list) {
        # Ne pas rejoindre le thread principal ni nous-mêmes
        if ($thr->tid && !threads::equal($thr, threads->self)) {
       	 $thr->join;
        }
    }
unlink $delete_bam if -e $delete_bam;

my $tabfilein;
my $dir_out= $project->getCallingPipelineDir("unifiedgenotyper");



foreach my $chr_name (@chr_names)    {
	my $namef = $patient_name.".".$chr_name;
	my $vcf3 = $dir_out."/".$namef.".merge.vcf.gz";
	next unless -e  $vcf3;
	delete $chr_tr->{$chr_name} if exists $chr_tr->{$chr_name} ;
	push (@$tabfilein,$vcf3); 	  
}
die("problem no vcf file ") if scalar(@$tabfilein) == 0;
die("problem  vcf file ".Dumper $chr_tr) if scalar(keys %$chr_tr) >0;
my $vcf1 = $tabfilein->[0];
die() unless -e $vcf1;
my $vcf_final = $dir_out."/".$patient_name.".final.vcf";
my $cmd = qq{zgrep "#" $vcf1 > $vcf_final};

system ($cmd);
foreach my $file (@$tabfilein){
	my $cmd = qq{zcat $file | $vcfannotate  -f d=$dp_limit | grep PASS | grep -v "#" >> $vcf_final};

	
	system($cmd);
	unlink $file;
}


open (VCF,$vcf_final);
my $nb =0;
my $nb_pub =0;
while(my $line = <VCF>){
	next if $line =~/#/;
	warn $line;
	chomp($line);
	$line =~ /DP=(\d*);/;
	 my $dp = $1;
	#  next if $dp < 10;
	my ($chr,$pos,$rs,$ref,$alt,@data) = split(" ",$line);
	next if  length ($ref) > 1 || length($alt)>1;

	my $chromosome = $project->getChromosome($chr);
	my $id = join("_", ($chromosome->name, $pos, $alt));
	#next unless $line=~/Intersection/;
		$nb++;
	my $res =  $project->public_data->dbsnp->get_variation(chr=>$chromosome->name ,id=>$id);
	$nb_pub ++ if exists  $res->{dbname};
}
close VCF;
my $p = int ($nb_pub/$nb *1000) / 10;



############comptages_quality check
my $text = $patient_name." snp : $nb dbnsp : ".$p." %";
my $color = 'black ON_BRIGHT_GREEN';
if ($p >= 90 && $p < 95) {
		$color = 'black ON_BRIGHT_YELLOW';
	}
	elsif  ($p < 90){
		$color = 'black ON_BRIGHT_RED';
	}
	print  colored [$color],$text  ;
	
	print  color 'reset';
	print  "\n";
	my $vcf_snp ="";
	 my $vcf_snp = $dir_out."/$patient_name.snp";
	 my $cmd1 = qq{$vcftools --vcf $vcf_final --remove-indels --out $vcf_snp --recode-INFO-all --recode >/dev/null};
	 system($cmd1);
	 rename $vcf_snp.'recode.vcf' , $vcf_snp.'vcf';
	  my $vcf_indel = $dir_out."/$patient_name.indel";
	  my $cmd2 = qq{$vcftools --vcf $vcf_final --keep-only-indels --out $vcf_indel --recode-INFO-all --recode >/dev/null};

	rename $vcf_indel.'recode.vcf' , $vcf_indel.'vcf';

exit(0);


sub calling {
my ($project_name,$patient_name,$bamfile,$chr_tr) = @_;

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
 my $patient = $project->getPatient($patient_name);

 
my $nb =0;
my $filein = $bamfile;

my $dir_out= $project->getCallingPipelineDir("haplotypecaller");
  
#while ( $File->pending() ) {
#	my @delete_file;
#	my $chr_num =$File->dequeue;
#	next if $chr_num eq "";
#	my $chr = $project->getChromosome($chr_num);
#	 my $bed_file = write_bed($patient->name,$chr,$dir_out,$chr_tr);
#	 next if $bed_file eq "";
#	 push(@delete_file,$bed_file);
#	my $chr_ucsc = $chr->ucsc_name;
#	my $namef = $patient->name.".".$chr->name;
#	my $vcf1 = $dir_out."/".$namef.".mpileup.vcf";
#	push(@delete_file,$vcf1);
#	my $vcf1a = $dir_out."/".$namef.".mpileup.1.vcf";
#	push(@delete_file,$vcf1a);
#	my $vcf2 = $dir_out."/".$namef.".uni.vcf";
#	push(@delete_file,$vcf2);
#	my $vcf2a = $dir_out."/".$namef.".uni.1.vcf";
#	push(@delete_file,$vcf2a);
#	my $vcf3 = $dir_out."/".$namef.".merge.1.vcf";
#	push(@delete_file,$vcf3);
##	my $vcf4 = $dir_out."/".$namef.".merge.vcf";
##	push(@delete_file,$vcf4);
#	my $t1 = time;
#	my $param_varscan2 = $param_varscan;
#	if ($chr_num eq "MT"){
#		#$param_varscan2 .= " --min-var-freq 0.001 --p-value 0.001";
#	}


	##Haplotype Caller
	my $gvcf ;
	my $cmd_haplotypecaller = qq{$javac -jar $gatk -T HaplotypeCaller $param_gatk  
		--emitRefConfidence GVCF --variant_index_type LINEAR 
		--variant_index_parameter 128000 -L  $bed_file -rf BadCigar -R $reference 
		--dbsnp $dbsnp  -I  $filein  -o $gvcf   -ntc  };
	system ("$cmd_haplotypecaller  -l off");

				}	
	

 }

}


