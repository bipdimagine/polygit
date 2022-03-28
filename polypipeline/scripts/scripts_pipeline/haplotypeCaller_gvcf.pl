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
use File::Temp;
 use Time::Elapsed qw( elapsed );
 use Time::ETA;
 
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
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"log=s" =>\$log_file,
	"patients=s" => \$patient_name,
	"chr=s"=>\$chr,
	"start=s"=>\$start,
	"end=s"=>\$end,
	"region=s"=>\$region_name,
	"fork=s" =>\$fork,
);
my $date = `date`;
chomp($date);

if ($log_file){
	open (STDOUT,">>".$log_file);
}

if ($region_name){
	my $reste;
	 ($chr,$reste) = split(":",$region_name);
	($start,$end) =  split("-",$reste);
}

my $other_project = [];
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );

my $samtools = $project->getSoftware('samtools');
my $javac = $project->getSoftware('java');
$javac = "java" unless -e $javac;
my $gatk  = $project->getSoftware('gatk');
my $project = $buffer->newProject( -name => $project_name );
my $chrs = $project->getChromosomes();
my $reference = $project->getGenomeFasta;
my $patients =  $project->get_list_patients($patient_name);

my $reference = $project->genomeFasta();
#my $fork = 4;
my $dir_out_gvcf= $project->getCallingPipelineDir("gvcf");
#my $dir_out_gvcf = $dir_out."/gvcf";
system ("mkdir -p $dir_out_gvcf") unless -e $dir_out_gvcf;
my $region = "$chr:$start-$end";
my $output_region = "$chr.$start.$end";
my $window =  $project->getChromosome($chr)->getWindow($start,$end,250);
#my $windows = $project->getChromosome($chr)->getWindowCaptureForCalling(250,1_000_000);
my $start_time = time;
#my ($window) = grep {$_->{start} eq $start && $_->{end} eq $end} @$windows;
die() unless $window;


my $nb = 0;
my $real_patient =[];
foreach my $patient (@$patients){
	my $outfile = $patient->getWindowGvcf($window);
	push(@$real_patient,$patient) ; 
	
}

my $bed_file = File::Temp->new( TEMPLATE =>"TMP.XXXXXXX",
                        DIR => $dir_out_gvcf,
                        SUFFIX => ".bed");
my @bed = map{$_ = $project->getChromosome($chr)->ucsc_name."\t".$_."\n"} split(";",$window->{intspan}->as_string({ sep => ";", range => "\t" }));
my $ll=0;
open(BED,">$bed_file");
foreach my $b (@bed){
	my ($c,$s,$e) = split(" ",$b);
	$ll+=abs($s-$e);
	print BED $b;
}
close BED;
my $out_gvcf;
my $total = scalar(@$real_patient);

my $now = time;

my $eta = Time::ETA->new(
    milestones =>scalar(@$real_patient),
);



my @cmds;


foreach my $patient (@$real_patient){
	my $outfile = $patient->getWindowGvcf($window);
	my $tmpoutfile = $outfile.".tmp.g.vcf";
	#my $outfile = $dir_out_gvcf."/".$patient->name().".".$window->{ext_gvcf};
	

	
	push(@$out_gvcf,$outfile);
	
	#next if -e $outfile;

	   my $bam =  $patient->getBamFile();
#	   my $line = `$samtools view -H $bam | grep "SN:chrM"`;
#	   my $version = $project->getVersion();
#	   
#	   if ($version eq "HG19b" && $line =~/16571/){
#	   		$reference =~ s/HG19b/HG19/;
#	   }
#	    if ($version eq "HG19" && $line =~/16569/){
#	    	
#	   		$reference =~ s/HG19/HG19b/;
#	   }
#	   

		#warn $reference;
  # next;
	my $recal = $patient->getRecalFile();
	my $recal_string="";
	$recal_string = "--BQSR $recal" if -e $recal;
	my $nct ="-nct 2";
	$nct ="-nct $fork" if $fork >1;
	my $cmd = qq{$javac -jar $gatk -T HaplotypeCaller -R $reference -I $bam -L $bed_file -o $tmpoutfile --emitRefConfidence GVCF    $nct -variant_index_type LINEAR -variant_index_parameter 128000 -l off};
	my $hcmd;
	$hcmd->{$outfile} = "$cmd && mv $tmpoutfile $outfile";
	push(@cmds,"$cmd && mv $tmpoutfile $outfile") unless -e $outfile; 

	    
	   
}
my $nb =0;
my $pm = new Parallel::ForkManager($fork);
$pm->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    		
    		 $eta->pass_milestone() ;
		# if ($nb % 2 == 0){
	    		my $p = int($eta->get_completed_percent());
			my $time1 = abs(time-$now);
			my $remaining = $eta->get_remaining_time();
			my $elapsed = $eta->get_elapsed_time();
			print " ** elapsed  $nb/$total :".$elapsed." remaining : ".$remaining." $region  :: $p % Done \n";
	    #}
    }
  );

foreach my $cmd (@cmds){
		$nb ++;
	my $pid = $pm->start and next;
	
	system($cmd);
	$pm->finish();

}
$pm->wait_all_children();
warn abs(time - $start_time);
if (abs(time - $start_time) < 50 ){
	warn "wait ";
	print "waiting for 180 seconds\n";
	#sleep(180);
}




unlink $bed_file;
print " END $region\n";
exit(0);
#	   $reference = $project->getGenomeFasta;
#my $string_variant = join (" --variant ",@$out_gvcf);
#my $vcf_file = $dir_out_vcf."/".$output_region.".vcf";
#my $cmd2 = qq{$javac -jar $gatk -T GenotypeGVCFs -R $reference --out $vcf_file --variant $string_variant -nt 4 };
#system($cmd2);
#print "end $chr $start $end";


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