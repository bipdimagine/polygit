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
	#"chr=s"=>\$chr,
	#"start=s"=>\$start,
	#"end=s"=>\$end,
	#"region=s"=>\$region_name,
	"fork=s" =>\$fork,
	"window=s" => \$window_length,
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
	die("not found bed file : $bed") unless -e $bed;

	
	my $cmd2 = " $javac -jar $gatk -I  $bamin -R $reference --disable_indel_quals  -T PrintReads   "	;#-BQSR $csv2   -o $fileout



my $bed = $patient->getCaptureBedFile() ;
die("not found bed file : $bed") unless -e $bed;
#my $table =  "$outputdir/$patient_name.table";
#unless (-e $table){
   my $recal_file = File::Temp->new( TEMPLATE =>"$patient_name.recal.XXXXXXX",
                        DIR => $outputdir,
                        SUFFIX => ".table",
                        CLEANUP => 1,
                        UNLINK =>1);
my  $table = $recal_file->filename();     
  
my $cmd1 =  "$javac -jar $gatk  -T BaseRecalibrator  -I  $bamin -R $reference -knownSites $vcf   --read_filter BadCigar  -L $bed -nct $fork -o $table --interval_padding 50";
$cmd1 .= " -knownSites $indels_gold" if ($indels_gold);
system($cmd1);

#system("cp $table $outputdir/$patient_name.table");
#}
#warn "end phase1 $outputdir/$patient_name.table";

#my $cmd2 = " $javac -jar $gatk -I  $bamin -R $reference --disable_indel_quals  -T PrintReads  -BQSR $table  -nct $fork -o $bamout ";
#system($cmd2);
#exit(0);
my $hwindow;
 foreach my $chr (@{$project->getChromosomes()}){
my $window = 50_000_000;
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
 
 	my $cmd2 = " $javac -jar $gatk -I  $bamin -R $reference --disable_indel_quals  -T PrintReads  -BQSR $table   "	;#-BQSR $csv2   -o $fileout
 	#my $cmd2 =  " $javac -jar $gatk  -I $bamin -R $reference -T IndelRealigner  -targetIntervals $csv  --BQSR $table  -rf NotPrimaryAlignment  " ;

my $fork2 = int($fork/4);
my $pm = new Parallel::ForkManager($fork2);	
my $nbf = 0;
my $files;
foreach my $h (@$hwindow){
	my $bam =  $h->{file}->filename();
	push(@$files,$bam);

	 my $pid = $pm->start and next;
	 warn $h->{r};
	  my $recal =  $cmd2." -o $bam -l off  -nct 4 -L ".$h->{r};
	  system($recal);
	  
	 $pm->finish();
}

$pm->wait_all_children();	

my $bams = join(" ",@$files);
#eval{
my $cmd_merge = "$sambamba  merge -t $fork $bamout $bams;rm $bams;rm $outputdir/$patient_name.recal*" ;
warn $cmd_merge;
system ($cmd_merge);
#};





#warn "$Bin/join_gvcf.pl -project=$project_name -patient=$patient_name -window=$window_length  -fork=$fork";
warn "END !!!!!";
exit(0);




sub getTmpFile {
	my ($dir,$chr_name,$ext) = @_;
	die() unless $ext;
	die() unless $chr_name;
	#confess($chr_name) unless $chr_name !~/[0,1,2,3,4,5,6,7,8,9,10,11,12,13,15,15,16,17,18,19,20,21,22,X,Y,M]/;
	$dir .="/$chr_name";
	system ("mkdir -p $dir") unless -d $dir;
	 my $file_cat_tmp =  File::Temp->new( TEMPLATE =>"TMP.XXXXXXX",
                        DIR => $dir,
                        SUFFIX => ".$ext");
  return $file_cat_tmp->filename();
}