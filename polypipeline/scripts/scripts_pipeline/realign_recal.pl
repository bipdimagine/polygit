#!/usr/bin/perl
use FindBin qw($Bin $RealBin);
use strict;

use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/"; 
use GBuffer ;
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
use File::Temp;
 use Time::Elapsed qw( elapsed );
 use Time::ETA;
use Storable qw(store retrieve freeze);

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
my $region;
my $start;
my $end;
my $chr;
my $region_name;
my $fork =1;
my $window_length;
my $bamin;
my $bamout;
my $fork = 1;
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"log=s" =>\$log_file,
	"patients=s" => \$patient_name,
	"filein=s"=> \$bamin,
	"fileout=s"=> \$bamout,
	"fork=s" =>\$fork,
);
my $date = `date`;
chomp($date);

if ($log_file){
	open (STDOUT,">>".$log_file);
}

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $java = $project->getSoftware('java');
my $samtools = $project->getSoftware('samtools');
my $javac = $project->getSoftware('java');
$javac = "java" unless -e $javac;
my $gatk  = $project->getSoftware('gatk');
my $chrs = $project->getChromosomes();
my $reference = $project->getGenomeFasta;
my $bgzip = $buffer->software("bgzip");
my $tabix = $buffer->software("tabix");
my $sambamba = $buffer->software("sambamba");
my $reference  =  $project->genomeFasta();
my $patient =  $project->getPatient($patient_name);

my $reference = $project->genomeFasta();
my $outputdir = $project->getAlignmentPipelineDir("bwa");

my $indels_gold = $project->gatk_indels_gold_file();
my $vcf = $project->gatk_dbsnp_file();

my $bed = $patient->getCaptureBedFile() ;

  my $recal_file = File::Temp->new( TEMPLATE =>"$patient_name.recal.XXXXXXX",
                        DIR => $outputdir,
                        SUFFIX => ".interval_list",
                        CLEANUP => 1,
                        UNLINK =>1);
my   $csv = $recal_file->filename();     
 
 my $cmd1 = "$javac -jar $gatk  -I  $bamin  -R $reference  -T RealignerTargetCreator  --read_filter BadCigar  -nt $fork -o $csv -L $bed --interval_padding 200";
 #my $cmd1 = "$javac -jar $gatk  -I  $bamin  -R $reference  -T RealignerTargetCreator  --read_filter BadCigar  ";
$cmd1 .= " -known $indels_gold" if ($indels_gold);
system($cmd1);

 warn "end realign step 1";
 
 #recal step 1
   my $recal_file = File::Temp->new( TEMPLATE =>"$patient_name.recal.XXXXXXX",
                        DIR => $outputdir,
                        SUFFIX => ".table",
                        CLEANUP => 1,
                        UNLINK =>1);
my   $table = $recal_file->filename();     
  
my $cmd11 =  "$javac -jar $gatk  -T BaseRecalibrator  -I  $bamin -R $reference -knownSites $vcf   --read_filter BadCigar  -L $bed -nct $fork -o $table --interval_padding 200";
$cmd11 .= " -knownSites $indels_gold" if ($indels_gold);
system($cmd11);
 warn "end recal step 1";

 
 
 
	#$cmd1 .= " -L $bed ";
	my $out =  $outputdir."/".$patient_name.".temp.bam";

 my $hwindow;
 foreach my $chr (@{$project->getChromosomes()}){
my $window = 100_000_000;
 	 my $from = 1;
 	 my $to = $chr->length;;
 	 my $start;
 	 my $end;
 	 	while ($from < $to){
        $start = $from;
        $end = $from + $window;
        if ($end > $to){
            $end = $to;
        }
           $from = $end;
           	my $chr_name = $chr->name();
           	my $htype;
           	$htype->{r} = $chr->ucsc_name.":".$start."-".$end;

           	 $htype->{file} = File::Temp->new( TEMPLATE =>$chr->name().".$patient_name.recal.XXXXXXX",
                        DIR => $outputdir,
                        SUFFIX => ".bam",
                        CLEANUP => 1,
                        UNLINK =>1);
                        
                   push(@$hwindow,$htype);     
 	 	}
 }
 
 	my $cmd2 =  " $javac -jar $gatk  -I $bamin -R $reference -T IndelRealigner  -targetIntervals $csv  --BQSR $table  -rf NotPrimaryAlignment  " ;

my $pm = new Parallel::ForkManager(int($fork/2));	
my $nbf = 0;
my $files;
foreach my $h (@$hwindow){
	my $bam =  $h->{file}->filename();
	push(@$files,$bam);

	 my $pid = $pm->start and next;
	 warn $h->{r};
	  my $recal =  $cmd2." -o $bam -l off  -L ".$h->{r};
	  system($recal);
	  
	 $pm->finish();
}

$pm->wait_all_children();	

my $bams = join(" ",@$files);
my $cmd_merge;

if (scalar(@$files) == 1){
	 $cmd_merge = "chmod a+rw $bams && mv $bams $bamout && sambamba index $bamout";
}
else {
	$cmd_merge = "$sambamba  merge -t $fork $bamout $bams;rm $bams;rm $outputdir/$patient_name.recal.*" ;
}
eval{
warn $cmd_merge;
system ($cmd_merge);
};

	
