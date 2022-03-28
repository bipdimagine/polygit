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
my $calling;
my $no_cluster;
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"log=s" =>\$log_file,
	"patient=s" => \$patient_name,
	"out=s"=>\$vcf_final,
	"low_calling"=>\$low_calling,
	"calling=s"=>\$calling,
	"no_cluster=s"=>\$no_cluster,
);

my $verbose;
$fork =8 unless $fork;
my $debug;
#$debug = "-l off"; 
my $date = `date`;

chomp($date);
my $methods = {
			"haplotypecaller"				=> { "method" => sub{calling_target::haplotypecaller(@_)}, priority =>1}, 
			"unifiedgenotyper"				=> { "method" => sub{calling_target::unifiedGenotyper(@_)}, priority =>2}, 
			"freebayes"				=> { "method" => sub{calling_target::freebayes(@_)}, priority =>3}, 
			"samtools"		=> { "method" => sub{calling_target::samtools(@_)},priority=>4},
};

foreach my $c (keys %$methods){
 	next if ($calling =~ /$c/);
 	delete $methods->{$c};
}
my $methods_name = join(",",keys %$methods);

if ($log_file && !$no_cluster){
	system ("  echo ---- start running calling with method : $methods_name => $patient_name $date >> $log_file");
	open (STDOUT,">>".$log_file);
}
else {
	$verbose =1;
}

print "  ---- start running calling : $methods_name => $patient_name $date\n";



my $buffer = GBuffer->new();

my $project = $buffer->newProject( -name => $project_name );
my $gatk  = $project->getSoftware('gatk');
my $java = $project->getSoftware('java');

	my $ref =  $project->genomeFasta();

my $patient = $project->getPatient($patient_name);
my $bam = $patient->getBamFile();
#warn $bam;
my $dir_out= $project->getCallingPipelineDir("unifiedgenotyper");


 $dir_out= $dir_out."/".$patient_name;
 system ("find $dir_out -name *TMP.* -exec rm {} \\;");
my @all_vcfs;

my @h_all_vcf;
#my $pr = String::ProgressBar->new( max => scalar(@{$project->getChromosomes()}) );
my $c =0;






#
 my $hregions = [];
 my $pm1 = new Parallel::ForkManager($fork);
 
 my $is_callable;
 
 my %all_intspan;
$pm1->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    	 # $pr->update($c++);
	
 		#$pr->write() if $verbose;
 		my $chr = $data->{chr};
 		my $ucsc = $data->{ucsc};
 			$all_intspan{$ucsc} = Set::IntSpan::Fast->new() unless exists 	$all_intspan{$ucsc};
 			foreach my $region_callable (@{$data->{array}}){
 				my $pos;
 				$pos->{chr} = $chr;
 		 	($pos->{start},$pos->{end}) = split("-",$region_callable);
 		 	$all_intspan{$ucsc}->add_range($pos->{start},$pos->{end});
 		 	$is_callable = 1;
    			push (@$hregions,$pos);
 			}
    		
    }
  );
  

 foreach my $chr (@{$project->getChromosomes()}){
 my $set = Set::IntSpan::Fast->new();
 my $intspan = $chr->getIntSpanCaptureForCalling(250);
 my @beds = $buffer->intspanToBed($chr,$intspan);
 next unless scalar @beds;

 my $bed1 =  calling_target::getTmpFile($dir_out,$chr->name,"bed");
 open(BED,">$bed1");
 print BED join("\n",@beds);
 close BED;
 my $summary =  calling_target::getTmpFile($dir_out,$chr->name,".sum.txt");
 my $bed =  calling_target::getTmpFile($dir_out,$chr->name,"bed");
 my $pid = $pm1->start and next;
 my $cmd = qq{$java -jar $gatk -T CallableLoci -R $ref -I $bam -o $bed -l off -summary $summary  --minDepth 10 --minDepthForLowMAPQ 15 -L $bed1};
 system($cmd);
  my $callable_intspan = Set::IntSpan::Fast::XS->new();
 	open (GATK, " cat $bed  | ");
 	 	while (my $line  = <GATK>){
 		chomp($line);
 	#	warn $line;
 		next unless $line =~/CALLABLE/;
 		my ($chr,$start,$end,$type) = split(" ",$line);
 		$callable_intspan->add_range($start, $end);
 	}
 	close GATK; 
 	
 	
 	my $intspan2 = $chr->getIntSpanCaptureForCalling(0);
 	my @regions_callable = split(",",$callable_intspan->as_string());
 	my $diff = $intspan2->diff($callable_intspan);
 	my $p = int(scalar($diff->as_array)/scalar($intspan2->as_array)*100);

 	my	$color = 'black ON_BRIGHT_GREEN';
 	my $text = "seems to be correct ";
 	
 	if ($p>10){
 		$color = 'black on_magenta';
 		$text= "COVERAGE WARNING PROBLEM :!!!! ";
 	}
 	if (scalar($intspan2->as_array) == 0 ){
 		$color = 'black ON_BRIGHT_RED';
 		$text= "COVERAGE PROBLEM : NO COVERAGE  => chromosome ".$chr->name;
 	}
	
	print  colored [$color],$patient->name."::".$chr->name." coverage ".$p."% ".$text ;
	
	print  color 'reset';
	print  "\n";
 	my $data;
 	$data->{chr} =  $chr->name;
 	$data->{ucsc} =  $chr->ucsc_name;
 	$data->{array} = \@regions_callable;
 	$pm1->finish(0,$data);
 	
 unlink $bed;
 unlink $summary;
 unlink $bed1;
 
 
 }
 $pm1->wait_all_children();



unless ($is_callable){
	my	$color = 'black ON_BRIGHT_RED';
 		my $text= "NO CALLABLE REGION FOR SAMPLE  ".$patient->name;
 		print  $text ;
 		exit(0);
}
	
### calling methods 
my $priority ="-priority ";
my $param_merge ="";
foreach my $type (sort {$methods->{$a}->{order} <=>$methods->{$b}->{order} }keys %$methods){
	
	my $vcf = $methods->{$type}->{method}($project,$patient,\%all_intspan,$low_calling,$fork,$verbose);
	$param_merge .=" --variant:".$type." ".$vcf;
	$priority .= $type.",";	
}


## merge results 

my $vcf_final = $dir_out."/".$patient->name.".final.vcf";
my $reference = $project->genomeFasta();
	
	my $javac = $buffer->getSoftware('java');
#$javac ="java";
my $gatk  = $buffer->getSoftware('gatk');

#my $output1  = getTmpFile($dir_out,"test","$sub.vcf");#$dir_out."/".$chr_name.".$end_ext.vcf";
#$output1 = "$dir_out/$chr_name/"

my $cmd_merge =qq{$javac -jar $gatk -T CombineVariants -R $reference  $param_merge  -o $vcf_final -genotypeMergeOptions PRIORITIZE $priority -l off};
system($cmd_merge);





my $bcftools =  $buffer->getSoftware('bcftools');
my $db_public = $project->lite_public_snps();
open (VCF,$vcf_final);
my $nb =1;
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
		my $hash =  $db_public->get($chromosome->name(),$pos);
 	#return $hash->{$self->sequence()} if $hash && exists $hash->{$self->sequence()};
	#my $res =  $project->public_data->dbsnp->get_variation(chr=>$chromosome->name ,id=>$id);
	$nb_pub ++ if $hash;
}
close VCF;
my $p =0;
 $p = int ($nb_pub/$nb *1000) / 10 if $nb ne 0;

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
	 my $vcf_snp = $dir_out."/$patient_name.snp.vcf";
	 my $cmd1 = qq{$bcftools view -v snps,mnps $vcf_final  -o $vcf_snp};
	# my $cmd1 = qq{$vcftools --vcf $vcf_final --remove-indels --out $vcf_snp --recode-INFO-all --recode >/dev/null};
	 system($cmd1);
	  my $vcf_indel = $dir_out."/$patient_name.indel.vcf";
	  my $cmd2 =  qq{$bcftools view -V snps,mnps $vcf_final  -o $vcf_indel}; #qq{$vcftools --vcf $vcf_final --keep-only-indels --out $vcf_indel --recode-INFO-all --recode >/dev/null};
 	system($cmd2);
	warn "end";
 foreach my $chr (@{$project->getChromosomes()}){
 	 my $dir_tmp = $dir_out."/".$chr->name;
 	 next unless -e $dir_tmp;
 	 system ("find $dir_tmp -name *TMP.* -exec rm {} \\;");
 	 system("rmdir $dir_tmp");
 }
#system ("find $dir_out -name *TMP.* -exec rm {} \\;");
exit(0);

