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
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"log=s" =>\$log_file,
	"patients=s" => \$patient_name,
	#"chr=s"=>\$chr,
	#"start=s"=>\$start,
	#"end=s"=>\$end,
	#"region=s"=>\$region_name,
	#"fork=s" =>\$fork,
	"window=s" => \$window_length,
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
my $java = $project->getSoftware('java');
my $samtools = $project->getSoftware('samtools');
my $javac = $project->getSoftware('java');
$javac = "java" unless -e $javac;
my $gatk  = $project->getSoftware('gatk');
my $project = $buffer->newProject( -name => $project_name );
my $chrs = $project->getChromosomes();
my $reference = $project->getGenomeFasta;
my $bgzip = $buffer->software("bgzip");
my $tabix = $buffer->software("tabix");
	my $ref =  $project->genomeFasta();
my $patient = $project->getPatientOrControl($patient_name);

my $reference = $project->genomeFasta();
#my $fork = 4;
my $dir_out_gvcf= $project->getCallingPipelineDir("gvcf");
#my $dir_out_gvcf = $dir_out."/gvcf";
$dir_out_gvcf .= "/$patient_name/";
system ("mkdir -p $dir_out_gvcf") unless -e $dir_out_gvcf;
my $region = "$chr:$start-$end";
my $output_region = "$chr.$start.$end";
my $bam =  $patient->getBamFile();

	my $recal = $patient->getRecalFile();
	my $recal_string="";
	$recal_string = "--BQSR $recal" if -e $recal;
	my $nct ="-nct 2";
	$nct="";
#	$nct ="-nct $fork" if $fork >1;
	
	my @cmds;
	my $real_fork = $fork;# int($fork /2);
	


my $dir_gvcf_out= $project->getCallingPipelineDir("gvcf");

my $dir_prod  = $patient->project->getGvcfDir("haplotypecaller");
my $final_gvcf = $dir_prod."/".$patient->name.".g.vcf.gz";
my $final_gvcf_tbi =  "$final_gvcf.tbi";

	
	
my $pm = new Parallel::ForkManager($real_fork);

warn $window_length;
$window_length =~s/_//g;# $window_length + 0.0 ;

my $file_freeze = $dir_out_gvcf."/$patient_name.regions.$window_length.freeze";

warn "$Bin/construct_regions_freeze.pl -project=$project_name -patient=$patient_name -window=$window_length -file=$file_freeze ";
unlink $file_freeze if -e $file_freeze;

system("$Bin/construct_regions_freeze.pl -project=$project_name -patient=$patient_name -window=$window_length -file=$file_freeze ") unless -e $file_freeze;
my $all_windows;
die() unless -e $file_freeze;
my $all_windows = retrieve($file_freeze);




my $total = 0;
foreach my $chr (@{$project->getChromosomes()}){
		next unless exists $all_windows->{$chr->name};
		 my $windows = $all_windows->{$chr->name};
		next unless $windows;		 
		 $total += scalar(@$windows);
}
my $nb;
$total += 1;
my $eta = Time::ETA->new(
    				milestones =>$total,
				);

my $lastp = 0;				
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

foreach my $chr (@{$project->getChromosomes()}){
			#next if $chr->name eq "MT";
		#	next unless $chr->name eq "18";
		next unless exists $all_windows->{$chr->name};
		  my $windows =$all_windows->{$chr->name} ;
		  next unless $windows;
		  foreach my $window (@$windows){
		  		$nb ++;
		  		my $outfile = $window->{outfile};#$patient->getWindowGvcf($window);
		  		#next if -e $outfile;
		  	      my $pid = $pm->start and next;
                         	if (-e $outfile){
                         		$pm->finish();
                         	}
                         	my $tmp_nam = $chr->name."_".$window->{start}."_".$window->{end}."_$nb";
                      
                        my $bed_file = File::Temp->new( TEMPLATE =>"$tmp_nam.XXXXXXX",
                        DIR => $dir_out_gvcf,
                        SUFFIX => ".bed");
                           
                   my @bed = map{$_ = $window->{chr_ucsc}."\t".$_."\n"} split(";",$window->{intspan}->as_string({ sep => ";", range => "\t" }));   
                 #  warn Dumper @bed;
                   my $ll=0;
				open(BED,">$bed_file");
				foreach my $b (@bed){
					my ($c,$s,$e) = split(" ",$b);
					$ll+=abs($s-$e);
			#		warn $b;
					print BED $b;
				}
				close BED;
				
					my $tmpoutfile = $outfile.".tmp.g.vcf";
					my $cmd = qq{$javac -jar $gatk -T HaplotypeCaller -R $reference -I $bam -L $bed_file  -o $tmpoutfile --emitRefConfidence GVCF    -variant_index_type LINEAR -variant_index_parameter 128000 -l off};
				#	warn $cmd;
					
					my $cmd2 = qq{$cmd && mv $tmpoutfile $outfile;rm $bed_file};
					#warn "start region";
					#warn $cmd;
					#sleep(5);
					system($cmd2);
					
					
					
				#	push(@cmds,"$cmd && mv $tmpoutfile $outfile && rm $bed_file") unless -e $outfile; 
				#	warn "end region $nb /$total";
		  		$pm->finish();
		  }#end window
}







$pm->wait_all_children();

warn "$Bin/join_gvcf.pl -project=$project_name -patient=$patient_name -window=$window_length  -fork=$fork";
warn("$Bin/join_gvcf.pl -project=$project_name -patient=$patient_name -window=$window_length  -fork=$fork  -file=$file_freeze && rm  $dir_out_gvcf"."/* && rmdir $dir_out_gvcf");
unlink $final_gvcf if -e $final_gvcf;
unlink $final_gvcf_tbi if -e $final_gvcf_tbi;
my $log_p = "";
 $log_p = "-log=$log_file" if $log_file;
warn "$Bin/join_gvcf.pl -project=$project_name -patient=$patient_name -window=$window_length  -fork=$fork  -file=$file_freeze && rm  $dir_out_gvcf"."/* && rmdir $dir_out_gvcf";
system("$Bin/join_gvcf.pl -project=$project_name -patient=$patient_name -window=$window_length  -fork=$fork  -file=$file_freeze && rm  $dir_out_gvcf"."/* && rmdir $dir_out_gvcf");




print "PROBLEM WITH $patient_name \n " unless -e $final_gvcf;;

die() unless -e $final_gvcf;;



unlink $file_freeze if -e $file_freeze;
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