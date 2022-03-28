#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use colored;
#use Set::IntSpan;
use GBuffer;
use Data::Dumper;
use Getopt::Long;
use Carp;
use calling_target;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use String::ProgressBar;
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
my $low_calling = undef;
my $chr_name;
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"log=s" =>\$log_file,
	"patient=s" => \$patient_name,
	"out=s"=>\$vcf_final,
	"low_calling"=>\$low_calling,
	"chr=s"=>\$chr_name,
);
my $verbose;
$fork =8 unless $fork;
my $debug;
#$debug = "-l off"; 
my $date = `date`;
chomp($date);
if ($log_file){
	system ("echo  ---- start running haplotypecaller $patient_name $date >> $log_file");
	open (STDOUT,">>".$log_file);
}
else {
	$verbose =1;
}
print "echo  ---- start running UnifiedGenotyper $patient_name $date\n";



my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatient($patient_name);
my $dir_out= $project->getCallingPipelineDir("unifiedgenotyper");
 $dir_out= $dir_out."/".$patient_name;
 #system ("find $dir_out -name *TMP.* -exec rm {} \\;");
my @all_vcfs;
my $pm = new Parallel::ForkManager($fork);
my @h_all_vcf;
my $pr = String::ProgressBar->new( max => scalar(@{$project->getChromosomes()}) );
my $c =0;
my $vcf_final = $dir_out."/".$patient_name.".final.$chr_name.vcf";



$pm->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    	    $pr->update($c++);
	
 		$pr->write() if $verbose;
    		push (@h_all_vcf,$data) if  $data->{vcf_final};
		    		
    		
    }
  );
my $nb = 0;

 #foreach my $chr (@{$project->getChromosomes()}){
 	my $chr = $project->getChromosome($chr_name);
 	my $intspan = $chr->getIntSpanCaptureForCalling(250);
 	my $java = $project->getSoftware('java');
 	
  #$java = "java";
	my $gatk  = $project->getSoftware('gatk');
	my $ref =  $project->genomeFasta();
	my $bam = $patient->getBamFile();
	my $summary =  calling_target::getTmpFile($dir_out,$chr->name,".sum.txt");
	my $bed =  calling_target::getTmpFile($dir_out,$chr->name,"bed");
	system ("$java -jar $gatk -T CallableLoci -R $ref -I $bam -o $bed -l off -summary $summary  --minDepth 10 --minDepthForLowMAPQ 15 -L ".$chr->ucsc_name);
 	open (GATK, " cat $bed  | grep CALLABLE |");
 	my $callable_intspan = Set::IntSpan::Fast::XS->new();
 	while (my $line  = <GATK>){
 		chomp($line);
 		my ($chr,$start,$end,$type) = split(" ",$line);
 		$callable_intspan->add_range($start, $end);
 	} 
 	 	close GATK;
 	unlink $bed;

  my $intspan_final = $intspan->intersection($callable_intspan);
 	
 	my @regions = split(",",$intspan_final->as_string());
 	$nb++;
 	if (scalar($intspan_final->as_array) == 0 ){
 		 	warn "$vcf_final";
 		 	warn "coucou ".@regions." ".join(";",@regions);
 		system ("touch $vcf_final");
 		system ("echo '##fileformat=VCFv4.2' >$vcf_final");
 		exit(0) ; 
 	}
 	
 	
 	foreach my $region (@regions){
	my $pid = $pm->start and next;
	my ($start,$end) = split("-",$region);
	my $vcf_final1 = calling_target::calling_merge( $project_name,$chr->name,$start,$end,$patient->name,$fork,{samtools=>3,freebayes=>2,unifiedgenotyper=>1},$dir_out,$Bin,$low_calling);
	my %h;
	$h{vcf_final} = $vcf_final1; 
	$h{num} = $nb; 
	$pm->finish(0,\%h);
	
	
 	}
#}
$pm->wait_all_children();
#warn "END CALLING CHROMOSOMES : $patient_name";

foreach my $h (sort{$a->{num} <=> $b->{num}}@h_all_vcf){
	push(@all_vcfs,$h->{vcf_final});
	
}
#warn Dumper @all_vcfs;
#die();
calling_target::concat_vcf(\@all_vcfs,$vcf_final,$project) ;

my $bcftools =  $buffer->getSoftware('bcftools');

open (VCF,$vcf_final);
my $nb =0;
my $nb_pub =0;
while(my $line = <VCF>){
	next if $line =~/#/;
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
my $p =0;
 if ($nb >0){
  $p = int ($nb_pub/$nb *1000) / 10;
 }
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
	 my $vcf_snp = $dir_out."/$patient_name.snp.$chr_name.vcf";
	 my $cmd1 = qq{$bcftools view -v snps,mnps $vcf_final  -o $vcf_snp};
	# my $cmd1 = qq{$vcftools --vcf $vcf_final --remove-indels --out $vcf_snp --recode-INFO-all --recode >/dev/null};
	 system($cmd1);
	  my $vcf_indel = $dir_out."/$patient_name.indel.$chr_name.vcf";
	  my $cmd2 =  qq{$bcftools view -V snps,mnps $vcf_final  -o $vcf_indel}; #qq{$vcftools --vcf $vcf_final --keep-only-indels --out $vcf_indel --recode-INFO-all --recode >/dev/null};
 	system($cmd2);
#cdwarn "end";
my $dir_out_tmp = $dir_out."/".$chr_name;

system ("find $dir_out_tmp -name *TMP.* -exec rm {} \\;");
system("rmdir $dir_out_tmp");
exit(0);


