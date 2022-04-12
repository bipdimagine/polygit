#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb";

use GBuffer;
use GenBoProject;
use Getopt::Long;
use Carp;
use strict;
use Storable qw(store retrieve freeze);

my $project_name;

my %hres;
my $fork= 5;
my $patient_name;
my $version;
# récupère les options  
GetOptions(
	'project=s'   => \$project_name,
	'patient=s'   => \$patient_name,
	'version=s'   => \$version,
	'fork=s' => \$fork,
);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name,-version=>"HG19");
my $patient = $project->getPatient($patient_name);
my $bamo = $patient->getBamFile();
my $chr = $project->getChromosome("chrM");
my $intspan = $chr->getIntSpanCaptureForCalling(300);
 my $dir_out = $project->getCallingPipelineDir("mt_bam_".$patient_name.".".time);
 my $bed = "$dir_out/mt.bed";
open (BED,">$bed");
print BED join("\n",$buffer->intspanToBed($chr,$intspan));
close (BED);


my $buffer_v = GBuffer->new();
my $project_v = $buffer_v->newProject( -name => $project_name,-version=>$version );
my $patient_v = $project_v->getPatient($patient_name);
my $bam = $patient_v->getBamFileName("bwa");

my $reference = $project_v->genomeFasta();
my $bazam = $buffer->software("bazam");
my $bwa = $buffer->software("bwa");
my $samtools  = $buffer->software("samtools");
my $bcftools  = $buffer->software("bcftools");
my $sambamba  = $buffer->software("sambamba");
 my @chrs_text = `$samtools idxstats $bamo | cut -f 1 | grep -v chrM | grep -v "*" `;
 chomp(@chrs_text);
 my $l = join(" ",@chrs_text);
create_bam();
create_gvcf();
create_haplotypecaller();
create_unifiedgenotyper();
create_freebayes();
create_samtools();
sub copy_dir {
	
}

sub create_bam {
return if -e $bam;
unlink $bam if -e $bam;
# my $bamo_tmp = $dir_out."/$patient_name.tmp.bam";
#my $header =$dir_out."/$patient_name.header.sam";
#system("$samtools view -H $bamo >$header");
#system("sed -i 's/16571/16569/g' $header");
warn "reheader";
my $t = time ;
#warn "samtools reheader $header $bamo > $bamo_tmp && samtools index  $bamo_tmp -@ $fork;";
#system("samtools reheader $header $bamo > $bamo_tmp && samtools index  $bamo_tmp -@ $fork;");
#warn "samtools reheader $header $bamo > $bamo_tmp && samtools index  $bamo_tmp -@ $fork;";
#warn "end header ".abs(time-$t);
my $b1 = $dir_out."/$patient_name.no_mt.bam";
my $b11 = $dir_out."/$patient_name.no.bam";
my $b1t = $dir_out."/$patient_name.tmp.mt.bam";
my $b2 = $dir_out."/$patient_name.mt.bam";
system("$sambamba slice $bamo $l > $b1");
system("$sambamba slice $bamo '*' > $b11");
system("$sambamba slice $bamo chrM > $b1t && $samtools index $b1t");
 my $ref_bwa =  $project_v->dirGenome().$project->buffer->index("bwa");
 die($ref_bwa) unless -e $ref_bwa;
 
 my $readgroup = qq{\@RG\\tID:$patient_name\\tSM:$patient_name\\tPL:ILLUMINA\\tDS:$patient_name\\tPG:bwa};
warn "$bazam -bam $b1t -L chrM:1-16572 | $bwa mem -t $fork -p -R '$readgroup' $ref_bwa  -  | samtools view -bSu - |  samtools sort -o $b2 -@ $fork";

system("$bazam -bam $b1t  | $bwa mem -t $fork -p -R '$readgroup' $ref_bwa  -  | samtools view -bSu - |  samtools sort -o $b2 -@ $fork");
warn "$samtools merge -f  $bam -h $b2 $b2 $b1  $b11 -@ $fork";
system ("$samtools merge -f  $bam -h $b2 $b2 $b1  $b11 -@ $fork");
warn "merge";
system ("$samtools index $bam -@ fork");
warn "index";
warn $bam;
}

sub create_gvcf {
	my $gvcf1 = $patient->getGvcfFile("haplotypecaller4",1);
	return unless -e $gvcf1;
	my $gvcf_tmp = "$dir_out/$patient_name.g.vcf";
	system("tabix -H $gvcf1 $l > $gvcf_tmp");
	system("sed -i 's/16571/16569/g' $gvcf_tmp");
	system("tabix  $gvcf1 $l >> $gvcf_tmp");
	
	my $gvcf_mt = "$dir_out/$patient_name.g.mt.vcf";
	system("gatk  HaplotypeCaller -I $bam -L $bed -O $gvcf_mt  -R $reference  -ERC GVCF");
	my $gvcf_out = $patient_v->gvcfFileName();
	unlink $gvcf_out if -e $gvcf_out;
	warn "*******************";
	warn "bcftools concat $gvcf_tmp $gvcf_mt -O z > $gvcf_out  && tabix -p vcf $gvcf_out";
	system("bcftools concat $gvcf_tmp $gvcf_mt -O z > $gvcf_out  && tabix -p vcf $gvcf_out");
}
sub create_haplotypecaller {
	my $vcf1 = $patient->getVariationsFile("haplotypecaller4",1);
	return unless $vcf1;
	
	my $vcf_tmp = "$dir_out/$patient_name.haplo.vcf";
	system("tabix -h $vcf1 $l > $vcf_tmp");
	my $vcf_mt = "$vcf_tmp.haplo";
	my $gvcf_out = $patient_v->gvcfFileName();
	warn $gvcf_out;
	system("gatk  GenotypeGVCFs -V $gvcf_out -L $bed -O $vcf_mt  -R $reference ");
	warn "gatk  GenotypeGVCFs -V $gvcf_out -L $bed -O $vcf_mt  -R $reference";
	#system("cat $vcf_mt | grep -v /^#/ >> $vcf_tmp ; bgzip $vcf_tmp ");
	my $dirout= $project_v->getVariationsDir("haplotypecaller4");
	
	my $fileout = $dirout."/$patient_name.vcf.gz";
	system("bcftools concat $vcf_tmp $vcf_mt -O z > $fileout  && tabix -p vcf $fileout");	
	
}
sub create_freebayes {
	my $vcf1 = $patient->getVariationsFile("freebayes",1);
	return unless -e $vcf1;
	my $vcf_tmp = "$dir_out/$patient_name.free.vcf";
	system("tabix -h $vcf1 $l > $vcf_tmp");
	my $vcf_mt = "$vcf_tmp.free";
	my $gvcf_out = $patient_v->gvcfFileName();
	my $gatk   = $buffer->software("gatk");
	my $freebayes   = $buffer->software("freebayes");
	my $freebayes_min_alternate = 0.08;
	my $correct_calling = qq{|  $Bin/../correct_calling_freebayes.pl -bam=$bam -samtools=$samtools };
	my $cmd_free = qq{  $freebayes -b $bam -f   $reference --min-coverage 20 -0  -F $freebayes_min_alternate  -t $bed $correct_calling >$vcf_mt 2>/dev/null};
	
	system($cmd_free);
	#system("cat $vcf_mt | grep -v /^#/ >> $vcf_tmp ; bgzip $vcf_tmp ");
	my $dirout= $project_v->getVariationsDir("freebayes");
	my $fileout = $dirout."/$patient_name.vcf.gz";
	system("bcftools concat $vcf_tmp $vcf_mt -O z > $fileout  && tabix -p vcf $fileout");	
	warn $fileout;
	
	#system("cat $vcf_mt | grep -v /^#/ >> $vcf_tmp ; bgzip $vcf_tmp ");
	#
	#system("mv $vcf_tmp.gz $fileout && tabix -p vcf $fileout");
	#warn $fileout;
}
sub create_samtools {
	my $vcf1 = $patient->getVariationsFile("samtools",1);
	my $vcf_tmp = "$dir_out/$patient_name.sam.vcf";
	return unless -e $vcf1;
	system("tabix -h $vcf1 $l > $vcf_tmp");
	my $vcf_mt = "$vcf_tmp.sam";
	my $gvcf_out = $patient_v->gvcfFileName();
	my $gatk   = $buffer->software("gatk");
		my $vcfutil     = $buffer->software("vcfutils");
	my $filter_samtools = "$Bin/../filter_samtools.pl";
	my $cmd_mpileup = qq{$bcftools mpileup -Ou -f $reference $bam -R $bed | bcftools call -vmO v | grep -v "ID=GL" |  $vcfutil  varFilter -a 5 -d 15   | $bcftools norm -f $reference  - | $filter_samtools > $vcf_mt 2>/dev/null};
	
	system($cmd_mpileup);
	#system("cat $vcf_mt | grep -v /^#/ >> $vcf_tmp ; bgzip $vcf_tmp ");
	my $dirout= $project_v->getVariationsDir("samtools");
	my $fileout = $dirout."/$patient_name.vcf.gz";
	system("bcftools concat $vcf_tmp $vcf_mt -O z > $fileout  && tabix -p vcf $fileout");	
	warn $fileout;
	
	#system("cat $vcf_mt | grep -v /^#/ >> $vcf_tmp ; bgzip $vcf_tmp ");
	#
	#system("mv $vcf_tmp.gz $fileout && tabix -p vcf $fileout");
	#warn $fileout;
}
sub create_unifiedgenotyper {
	my $vcf1 = $patient->getVariationsFile("unifiedgenotyper",1);
	return unless -e $vcf1;
	my $vcf_tmp = "$dir_out/$patient_name.uni.vcf";
	system("tabix -h $vcf1 $l > $vcf_tmp");
	my $vcf_mt = "$vcf_tmp.uni";
	my $gvcf_out = $patient_v->gvcfFileName();
	my $gatk   = $buffer->software("gatk");
	my $javac   = $buffer->software("java");
	my $cmd_uni = qq{$javac  -jar $gatk  -T UnifiedGenotyper   -rf BadCigar -R $reference -L $bed   -I $bam  --genotype_likelihoods_model BOTH   -o $vcf_mt };
	system($cmd_uni);
	#system("cat $vcf_mt | grep -v /^#/ >> $vcf_tmp ; bgzip $vcf_tmp ");
	my $dirout= $project_v->getVariationsDir("unifiedgenotyper");
	my $fileout = $dirout."/$patient_name.vcf.gz";
	system("bcftools concat $vcf_tmp $vcf_mt -O z > $fileout  && tabix -p vcf $fileout");	
	
	#system("cat $vcf_mt | grep -v /^#/ >> $vcf_tmp ; bgzip $vcf_tmp ");
	#
	#system("mv $vcf_tmp.gz $fileout && tabix -p vcf $fileout");
	#warn $fileout;
	
}

#  instanciation d'un objet project

#bcftools view HYP_6122GM000520.vcf.gz -t ^chrM

