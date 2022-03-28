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
	die("not found bed file : $bed") unless -e $bed;
	my $cmd1 = "$javac -jar $gatk  -I  $bamin  -R $reference  -T RealignerTargetCreator  --read_filter BadCigar  ";
	$cmd1 .= " -known $indels_gold" if ($indels_gold);
	#$cmd1 .= " -L $bed ";
	
	
	my $cmd2 =  " $javac -jar $gatk  -I $bamin -R $reference -T IndelRealigner   -rf NotPrimaryAlignment --interval_padding 100" ;

#my $fork = 4;
#my $dir_out_gvcf = $dir_out."/gvcf";

	
my $pm = new Parallel::ForkManager($fork);

my $total = 0;
my $bam_temp;
my $bams="";
my @bams_tmp;
foreach my $chr (@{$project->getChromosomes()}){
		 $total ++;
		     $bam_temp->{$chr->name}->{bam} = File::Temp->new( TEMPLATE =>$chr->name().".$patient_name.recal.XXXXXXX",
                        DIR => $outputdir,
                        SUFFIX => ".bam",
                        CLEANUP => 1,
                        UNLINK =>1);
                          $bam_temp->{$chr->name}->{realign} = File::Temp->new( TEMPLATE =>$chr->name().".$patient_name.recal.XXXXXXX",
                        DIR => $outputdir,
                        SUFFIX => ".interval_list",
                        CLEANUP => 1,
                        UNLINK =>1);
            $bams .=   $bam_temp->{$chr->name} ->{bam}->filename()." ";          
            push (@bams_tmp,$bam_temp->{$chr->name} ->{bam}->filename());
            
            
}

                   
#warn $cmd1;
#die();                       
my $nb;
$total += 1;
my $eta = Time::ETA->new(
    				milestones =>$total,
				);

my $lastp = 0;	
my $fork1 = int($fork/2);			
my $pm = new Parallel::ForkManager($fork1);
$pm->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    		
    		 $eta->pass_milestone() ;
		# if ($nb % 2 == 0){
	    		my $p = int($eta->get_completed_percent());
	    		if ($p ne $lastp && $p % 5 == 0){
	    			
			my $remaining = $eta->get_remaining_time();
			my $elapsed = $eta->get_elapsed_time();
			print " ** elapsed   $patient_name :".$elapsed." remaining : ".$remaining." $region  :: $p % Done \n";
			$lastp = $p;
	    		}
	    #}
    }
  );
  
warn "START computing $patient_name";
my @chromosomes = (1,16,2,22,3,21,4,20,5,19,6,18,7,'17',8,16,9,10,11,12,13,14,15,'X','Y','MT');

foreach my $chr_name (@chromosomes){
	my $chr = $project->getChromosome($chr_name);
			#next if $chr->name eq "MT";
		#	next unless $chr->name eq "18";

              my $pid = $pm->start and next;
               warn "start $chr_name realign 1";
                 my $out2 =  $bam_temp->{$chr->name} ->{bam}->filename();
               	my $csv = $bam_temp->{$chr->name}->{realign}->filename();
                my $cmd_recal1 = $cmd1." -nt 2 -o $csv -l off -L ".$chr->ucsc_name();
                system ($cmd_recal1);
                  warn "start $chr_name realign 2";
                my $recal =  $cmd2." -targetIntervals $csv -o $out2 -l off  -L ".$chr->ucsc_name();
                #$cmd2 -BQSR $recal   -o $out2 -L ".$chr->ucsc_name();
                system($recal);
                  warn "end $chr_name ";
		  		$pm->finish();
}


$pm->wait_all_children();
my $cmd_merge = "$sambamba  merge -t $fork $bamout $bams" ;
warn $cmd_merge;
system ($cmd_merge);
warn "end";
foreach my $b (@bams_tmp){
	unlink $b if -e $b;
	unlink $b.".bai" if -e $b.".bai";
	$b =~s/\.bam/\.bai/;
	unlink $b if -e $b;
 }
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