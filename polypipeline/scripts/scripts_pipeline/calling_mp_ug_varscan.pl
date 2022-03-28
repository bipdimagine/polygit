#!/usr/bin/perl
use FindBin qw($Bin);
use strict;


use lib "$Bin/../../../GenBo/lib/obj-nodb/";


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
my $debug;
#$debug = "-l off"; 
my $date = `date`;
chomp($date);
if ($log_file){
	system ("echo  ---- start running haplotypecaller $patient_name $date >> $log_file");
	open (STDOUT,">>".$log_file);
}
print "echo  ---- start running UnifiedGenotyper $patient_name $date\n";



my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $build = $project->buffer->build();

#my $file_bed = "/data-xfs/public-data/$build/capture/agilent/agilent.v50.bed";
my $samtools = "/bip-d/soft/bin/samtools";
my $bcftools = "/bip-d/soft/bin/bcftools";
my $vcfutil ="/bip-d/soft/bin/vcfutils.pl";
my $vcftools = qq{/bip-d/soft/distrib/vcftools/latest/bin/vcftools };

my $javac 		= qq{/opt/java/latest/bin/java  };
my $varscan 	=" $javac -jar /bip-d/soft/distrib/VarScan/VarScan.v2.3.6.jar ";
my $dbsnp 		= "/data-xfs/public-data/$build/snp/gatk/latest/dbsnp.vcf.gz";
my $genomes1k	= "/data-xfs/public-data/$build/snp/gatk/latest/1000genomes.vcf.gz";
my $hapmap		= "/data-xfs/public-data/$build/snp/gatk/latest/hapmap.vcf.gz";

#my $gatk = qq{/bip-d/soft/distrib/GATK/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar};
my $gatk = qq{/bip-d/soft/distrib/GATK/gatk-latest/GenomeAnalysisTK.jar};
#my $gatk = qq{/bip-d/soft/distrib/GATK/gatk-latest.2.7/GenomeAnalysisTK.jar};
my $vcfannotate =qq{/bip-d/soft/distrib/vcftools/latest/bin/vcf-annotate}; 
my $chrs = $project->getChromosomes();
my $reference = $project->getGenomeFasta;
my $patient = $project->getPatient($patient_name);
my $coverage_dir = $project->getRootDir()."/align/coverage/";
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
	my $dp = `/bip-d/soft/bin/tabix $coverage_file mean_all:99  | cut -f 3`;
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




my @chr_names;
foreach my $chr (@$chrs)    {

	my $chromosome = $project->getChromosome($chr->name);
	next unless $chromosome;
	my $id = $chr->name();
	$id = 23 if ($id eq "X" );
	$id = 24 if ($id eq "Y" );
	$id = 25 if ($id eq "MT" );
	$id -= 1;
	die($chr_names[$id]) if $chr_names[$id];
	$chr_names[$id] = $chr->name();
	#push(@chr_names,$chr->name());
}
@chr_names = grep {$_ ne ""} @chr_names;
$fork = scalar(@chr_names) if $fork > scalar(@chr_names);
my $File = Thread::Queue->new;

@chr_names = grep {$_ ne ""} @chr_names;
$File->enqueue(@chr_names);
my $thr;
for (my $i=0;$i<$fork;$i++){
 	 $thr->[$i] = threads->new(\&calling,$project_name,$patient_name,$bamfile,$chr_tr);
}


foreach $thr (threads->list) {
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

my $dir_out= $project->getCallingPipelineDir("unifiedgenotyper");
  
while ( $File->pending() ) {
	my @delete_file;
	my $chr_num =$File->dequeue;
	next if $chr_num eq "";
	my $chr = $project->getChromosome($chr_num);
	 my $bed_file = write_bed($patient->name,$chr,$dir_out,$chr_tr);
	 next if $bed_file eq "";
	 push(@delete_file,$bed_file);
	my $chr_ucsc = $chr->ucsc_name;
	my $namef = $patient->name.".".$chr->name;
	my $vcf1 = $dir_out."/".$namef.".mpileup.vcf";
	push(@delete_file,$vcf1);
	my $vcf1a = $dir_out."/".$namef.".mpileup.1.vcf";
	push(@delete_file,$vcf1a);
	my $vcf2 = $dir_out."/".$namef.".uni.vcf";
	push(@delete_file,$vcf2);
	my $vcf2a = $dir_out."/".$namef.".uni.1.vcf";
	push(@delete_file,$vcf2a);
	my $vcf3 = $dir_out."/".$namef.".merge.1.vcf";
	push(@delete_file,$vcf3);
#	my $vcf4 = $dir_out."/".$namef.".merge.vcf";
#	push(@delete_file,$vcf4);
	my $t1 = time;
	my $param_varscan2 = $param_varscan;
	if ($chr_num eq "MT"){
		$param_varscan2 .= " --min-var-freq 0.02 --p-value 0.001";
	}
	my $priorityze = {
		uni => 1,
		mp => 2,
	};
	my $param_merge;
		
	## Varscan

	warn "chromsome  ... ".$chr->name;

	my $vcf_v1 = $dir_out."/".$namef.".varscan.vcf";
	push(@delete_file,$vcf_v1);

		
#	my $cmd_varscan1 = qq{$samtools mpileup -Bf $reference  $filein -l $bed_file $param_mp | $varscan mpileup2cns  --output-vcf 1 --variants 1 $param_varscan2 > $vcf_v1};
	my $cmd_varscan1 = qq{$samtools mpileup -f $reference  $filein -l $bed_file -I -Q0 -d 50000| $varscan mpileup2cns  --output-vcf 1 --variants 1 $param_varscan2 > $vcf_v1};
	warn $cmd_varscan1 ;
#	system ("$cmd_varscan1  2>/dev/null");
	warn $param_varscan2 ;
	system ("$cmd_varscan1  ");
	if (-s  $vcf_v1){
		my $vcf_varscan = change_sample_id($vcf_v1,$patient);
		push(@delete_file,$vcf_varscan);
		$param_merge ="--variant:vs $vcf_varscan ";
		$priorityze->{vs} = 3;
	}	

	## Samtools		
	my $vcf_mpileup;
		my $cmd_mpileup = qq{$samtools mpileup -ugf $reference  $filein -l $bed_file $param_mp | $bcftools view -vcg - |  $vcfutil  varFilter -a $al -d $dp_limit  | grep -v "ID=GL" > $vcf1a};
		warn $cmd_mpileup;
		system ("$cmd_mpileup  ");
	 	$vcf_mpileup = change_sample_id($vcf1a,$patient);
		push(@delete_file,$vcf_mpileup);
		$param_merge ="--variant:mp $vcf_mpileup " .$param_merge ;
		
	## Unified Genotyper
	my $vcf_uni ;
	my $cmd_unified = qq{$javac -jar $gatk -T UnifiedGenotyper $param_gatk   -L  $bed_file -rf BadCigar -R $reference --dbsnp $dbsnp  -I  $filein  -o $vcf2a   --genotype_likelihoods_model BOTH -nt 2 };
	warn $cmd_unified;
	system ("$cmd_unified  $debug");
			my @t = `cat  $vcf2a`;
		warn Dumper @t; 
	$vcf_uni = change_sample_id($vcf2a,$patient);
	push(@delete_file,$vcf_uni);
	$param_merge ="--variant:uni $vcf_uni   " .$param_merge ;
			
	## mergegrep 
	my $pr = join(",", sort  {$priorityze->{$a} <=> $priorityze->{$b} }  keys %$priorityze);
	warn Dumper $priorityze ;
	my $cmd_merge =qq{$javac -jar $gatk -T CombineVariants -R $reference  $param_merge  -o $vcf3 -genotypeMergeOptions PRIORITIZE -priority $pr };
	system($cmd_merge." $debug  " );
	warn $cmd_merge ;

	my $vcf_final = $dir_out."/".$namef.".merge.vcf.gz";
	my $vcf4 = filter_uni_vcf($vcf3,$patient);
	rename $vcf4,$vcf_final;
	rename $vcf4.".tbi",$vcf_final.".tbi";
	foreach my $f (@delete_file){
		unlink $f if -e $f;
		unlink $f.".gz" if -e $f.".gz";
		unlink $f.".idx" if -e $f.".idx";
	}	
	

 }

}

sub write_bed {
	my ($name,$chr,$dir_out,$chr_tr) = @_;
	my $span_final;
	my $project = $chr->getProject();
	my $final_span;
	my $span_chr = Set::IntSpan::Fast::XS->new("1-".$chr->length);
	if ($chr->getProject()->isDiagnostic){
		 $final_span = $chr->getIntSpanCapture(100);
		 $final_span = $final_span->union($chr_tr->{$chr->name}) if exists $chr_tr->{$chr->name};
		
	}
	else {
		my $span = $chr->getIntSpanCapture(100);
		my $file2 = "/data-xfs/public-data/$build//ccds/".$chr->name().".bed";
		if (-e $file2){
				my $span_ccds = retrieve($file2);
				$span_final = $span->union($span_ccds);
		}
		else {
			die("no ccds file");
		}
	}
	

	return "" if $final_span->is_empty();
	
	my $span4 = $span_chr->intersection($final_span);
	
	my $nb4 = scalar($span4->as_array());
	my @bed = map{$_ = $chr->ucsc_name."\t".$_} split(";",$span4->as_string({ sep => ";", range => "\t" }));


my $file_bed = $dir_out."/".$chr->name.".$name.bed";
unlink $file_bed if -e $file_bed;

open(BED,">$file_bed") or die("not open $file_bed");
print BED join("\n",@bed);
print BED "\n";
close BED;
return $file_bed;
}



sub filter_uni_vcf {
		my ($vcf_in,$patient,$filter) = @_;
		

		my $vcf_out = $vcf_in;
		$vcf_out =~s/vcf/filter\.vcf/;
		
		my $bam = $patient->getBamFile();
	     my $sam = Bio::DB::Sam->new(-bam  =>$patient->getBamFile,
                             -fasta=>$reference,
                             );
		my $project = $patient->getProject();
		die($vcf_in) unless -e $vcf_in;
		my $bam = $patient->getBamFile();
		
		my $dir_out= $patient->getProject->getCallingPipelineDir("unifiedgenotyper");
		my $span;
		foreach my $chr (@{$project->getChromosomes()}) {
			$span->{$chr->ucsc_name} = $chr->getIntSpanCapture(1);
			
		}
		
		
		if ($vcf_in =~ "gz"){
		open (VCF_IN, "zcat $vcf_in |");
	}
	else {
	open (VCF_IN, $vcf_in);
	}
	open (VCF_OUT,">".$vcf_out) || die(); 
	while (my $line = <VCF_IN>){
		if ($line =~ /#/){
			print VCF_OUT $line;
			next;
		}
			my @data = split(" ",$line);
			my $pos = $data[1];
			my $chr = $data[0];
			my $a2 = $data[4];
			if ($line =~/DP4/ && $line=~ /mp/){
 				my ($DP) = grep {$_=~ /DP\=/} split(";",$data[7]);
 				my ($DP4) = grep {$_=~ /DP4\=/} split(";",$data[7]);
 				my ($n,$z) = split("=",$DP4);
 				my ($a1,$a2,$a3,$a4) = split(",",$z);
 				my ($n1,$z1) = split("=",$DP);
 				next if ($a3+$a4) < $z1 *0.15 && $z1 > 500; 
 			}			
				
				my $nb_start =0;
				my $nb_end =0;
				my $nb_read =0;
				
				my $callback = sub {
         			my ($seqid,$pos1,$pileups) = @_;
         			return  if $pos1 ne $pos ;
         			foreach my $pileup (@$pileups){
         				 my $b     = $pileup->alignment;
         				 my $qbase  = substr($b->qseq,$pileup->qpos,1);
         				 next if  $a2 !~/$qbase/;
         				 $nb_read++;
         				 $nb_start ++ if $pileup->pos < 10;
         				 $nb_end ++ if $pileup->pos > length($b->qseq) -10;
         				 
         			}
 			};
 			
 			my $end = $pos+1;
 			$sam->fast_pileup("$chr:$pos-$end",$callback);
 			#die() if $debug;
 			next if $nb_read == 0;
 			next if ($nb_start+$nb_end)/$nb_read > 0.8;
 			
				
			print VCF_OUT $line;
	}
	
	close VCF_IN;
	close VCF_OUT;
	
	system("/bip-d/soft/bin/bgzip -f $vcf_out;/bip-d/soft/bin/tabix -p vcf $vcf_out.gz");
	return "$vcf_out.gz";

		
}


sub change_sample_id {
	my ($vcf_in,$patient,$filter) = @_;
	my $vcf_out = $vcf_in;
	 $vcf_out =~ s/vcf/sample\.vcf/;
	my $name = $patient->name;
	die($vcf_in . "=> not found") unless -e $vcf_in;
	my $bam = $patient->getBamFile();
	open (VCF_OUT,">".$vcf_out);
	if ($vcf_in =~ "gz"){
		open (VCF_IN, "zcat $vcf_in |");
	}
	else {
		open (VCF_IN, $vcf_in);
	}
	my %format;
	while (my $line = <VCF_IN>){
		chomp($line);
		if ($line =~ /#CHROM/){
			
		my @toto = split(" ",$line);
		die(scalar(@toto)." ".join(";",@toto)) if scalar(@toto) != 10;
		$toto[-1] = $name;
		$line = join("\t",@toto);
		}
		
		print VCF_OUT $line."\n"; 
	}
	close VCF_IN;
	close VCF_OUT;
	unlink $vcf_in;
	system("/bip-d/soft/bin/bgzip -f $vcf_out;/bip-d/soft/bin/tabix -p vcf $vcf_out.gz");
	return "$vcf_out.gz";
}


