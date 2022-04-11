#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use lib "$FindBin::Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$FindBin::Bin/../../../packages/";
use  File::Temp;
use Data::Dumper;
use Getopt::Long;
use Carp;  
use GBuffer; 
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
 use IPC::Open2;
 use Time::Elapsed qw( elapsed );
  use Time::ETA;
 use calling_target; 
my $project_name;
my $final_vcf;
my $log_file;
my $window_length;
my $list_patients;
my $fork;
my $now = time;
my $file;
my $version;
GetOptions(
	'project=s'   => \$project_name,
	"log=s" =>\$log_file,
	"vcf=s" => \$final_vcf,
	"window=s" => \$window_length,
	"patient=s"=>\$list_patients,
	"fork=s" =>\$fork,
	"file=s"=>\$file,
	"version=s" => \$version
);

my $date = `date`;
chomp($date);
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name, -version=>$version );

my $gatk  = $project->getSoftware('gatk');
my $java = $project->getSoftware('java');
my $vcffirstheader = $buffer->software("vcffirstheader");
my $vcfstreamsort = $buffer->software("vcfstreamsort");
my $vcfuniq  = $buffer->software("vcfuniq");
my $bcftools = $buffer->software("bcftools");
my $bgzip = $buffer->software("bgzip");
my $tabix = $buffer->software("tabix");
#die();
if ($log_file){
	open (STDOUT,">>".$log_file);
}
#$SIG{INT} = \&interrupt;


	my $ref =  $project->genomeFasta();



my $patient = $project->getPatient($list_patients);


die("problem with $file") unless -e $file;
my $all_windows = retrieve($file);
#delete $all_windows->{"MT"};
my $total2 =0;
foreach my $chr (@{$project->getChromosomes()}){
#	next if $chr->name eq "MT";
		next unless exists $all_windows->{$chr->name};
	my $windows = $all_windows->{$chr->name};
	next unless $windows;
	$total2+= scalar(@$windows);
}

	my $eta2 = Time::ETA->new(
    				milestones =>$total2,
				);
$| =1;
my $files = [];
	my $delete_file =[];
my $pm = new Parallel::ForkManager($fork);
my $nb =0;
my $lastp =0;
my $nb_error;
$pm->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    	$nb ++;
    		push (@$delete_file,$data->{error}) if exists $data->{error};
    		warn "error $data->{error}" if exists $data->{error};
    	 $eta2->pass_milestone() ;
	    		my $p = int($eta2->get_completed_percent());
	    		if ($p ne $lastp && $p % 5 == 0){
	    			
			my $remaining = $eta2->get_remaining_time();
			my $elapsed = $eta2->get_elapsed_time();
			my $nb_error = scalar(@$delete_file);
			print " ** checking gvcf file elapsed   $list_patients :".$elapsed." remaining : ".$remaining." :: $p % Done  :: Error detected :: $nb_error \n";
			$lastp = $p;
	    		}
				    
    	
    }
  );
  
  print "check GVCF \n";

		foreach my $chr (@{$project->getChromosomes()}){
		#	next if $chr->name eq "MT";
			#next unless $chr->name eq "21";
			next unless exists $all_windows->{$chr->name};
				my $windows =  $all_windows->{$chr->name};
				next unless $windows;
				my $intspan0 = $chr->getIntSpanCaptureForCalling(0);
				#warn $chr->name;
			#	my $chr_gvcf = $final_gvcf.".".$chr->name;
			
				foreach my $window (@{$windows}){
					
				next if $window->{intspan}->is_empty();
				
				 my $bam = $patient->getBamFile();
				 	my $outfile = $window->{outfile};
			
				push (@$files,$outfile);    		
				 my $pid = $pm->start and next;

 					my $callable_intspan = $intspan0->intersection($window->{intspan});
				 my $size  = scalar($callable_intspan->as_array);
				
				 my %h;
				my $ok_file = $outfile.".ok";
				
				if (-e  $ok_file ){
					$h{file} = $outfile;
					$pm->finish(0,\%h);
				}
				
				
		
			
			
				
				
				unless (-e $outfile){
						$h{error} = $outfile;
					print "MISSING FILE  : $outfile \n";
					$pm->finish(0,\%h);
				
					
				}
				open (IN,"cat $outfile | grep -v '#' | cut -f 1,2,4,8 |");
			
  				my $intspan = Set::IntSpan::Fast::XS->new();
  				my $pos =0;
  				my $ee;
  				my $real_start;
  				my $real_end;
  				while(my $line = <IN>){
  					chomp($line);
  					my ($chr_name,$start,$all,$end) = split(" ",$line);
  					if ($chr_name ne $chr->fasta_name){
  							$h{error} = $outfile;
  						
						print "ERROR CHROMOSOME ON FILE  : $outfile \n";
						$pm->finish(0,\%h);
  					}
  					unless ($real_start){
  					$real_start =$start;
  					
  					}
  					if ($start > $pos){
  						$ee =1;
  					}
  					if ($end =~/END/){
  						$end =~ s/END=//;
  						$intspan->add_range($start-1,$end);
  						$pos = $end;
  						 $real_end = $end; 
  					}
  					else {
  						my $end2 = $start + (length($all)) -1;
  						
  							$intspan->add_range($start,$end2);
  							$real_end = $end2;
  						#$intspan->add($start);
  						$pos++;
  					}
 					}
 					
 					my $size2 =  scalar($intspan->as_array);
 					my $diff = $callable_intspan->diff($intspan);
 					my $size3 =  scalar($diff->as_array);
 					my $intspan_global = Set::IntSpan::Fast::XS->new($window->{start}."-".$window->{end});
 					my $diff2 = $intspan_global->diff($intspan);
 					
 					if (abs($real_start - $window->{start}) > 2){
 						$h{error} = $outfile;
  						
						print "ERROR START  : $real_start  : $outfile";
						$pm->finish(0,\%h);
 					}
 					
 					
 					if ( $window->{end} >$real_end ) {
 						$h{error} = $outfile;
						print "ERROR END  : $real_end ".$window->{end}." : $outfile ".abs($real_end - $window->{end})."\n";
						$pm->finish(0,\%h);
 					}
 					
  				if ($size3 > 1000  && $outfile !~/chrM/ && $outfile !~/chrY/ ){
  				my $iter = $window->{intspan}->iterate_runs();
   					 while (my ( $from, $to ) = $iter->()) {
   					 	my $ii = Set::IntSpan::Fast::XS->new($from."-".$to);
   					 	my $iz = $ii->diff($intspan);
   					 	warn $from."-$to ".$iz->as_string()." ".$outfile;
   					 }
  				
  					if (check_callable($chr,$window,$bam,$intspan0) ){
  						$h{error} = $outfile;
  						
						print "ERROR SIZE ON FILE : $size3 : $outfile";
						warn "ERROR SIZE  ON FILE : $size3: $outfile";
						$pm->finish(0,\%h);
  					}
						
  				}	
			
				system("$bgzip -f $outfile && $tabix -p vcf $outfile.gz  ");
				
				unless (-e "$outfile.gz.tbi" ){
						$h{error} = $outfile;
						print "ERROR ON FILE : $outfile";
						warn "ERROR ON FILE : $outfile";
						system("gunzip $outfile.gz");
						push(@$delete_file,$outfile);
				}
				else {
					$h{file} = $outfile;
					system("touch $ok_file");
					unlink $outfile.".gz.tbi";
					system("gunzip $outfile.gz");
				}
				
				
		
			
				$pm->finish(0,\%h);
				#push(@$files,$outfile);
				}
		}
		$pm->wait_all_children();
			print "check file =>".elapsed(abs(time - $now))."\n";

			if (@$delete_file){
				print "error detected !!!!! ".$patient->name."\n";
				
			foreach my $f (@$delete_file){
				print " \t problem with file : $f \n";
				next unless -e $f;
				unlink $f;
				print " problem with file : $f \n";
				
			}
		}
			$now = time;
			
			my $dir_gvcf_out= $project->pipeline_tmp_dir();
		my $final_gvcf = $dir_gvcf_out."/".$patient->name.".g.vcf";
		my $gz = $final_gvcf.".gz";
		my $tbi = $final_gvcf.".gz.tbi";
				
			my $total = scalar(@$files);
			my $step = int($total/10)+1;
			my $eta = Time::ETA->new(
    				milestones =>$total,
				);
			 $nb =0;	
				open(GLOBAL, " | $vcffirstheader | $vcfstreamsort | $vcfuniq > $final_gvcf");
			
	 			my @all_vcfs;
				foreach my $outfile (@$files){
					$nb ++;
					#print "$nb / $total \r";
					
					#next unless $outfile eq "/data-isilon/sequencing/pipeline/NGS2015_0622/HG19/calling/gvcf//LAU_FAG/LAU_FAG.chr11.1099307.3380847.g.vcf";
					open (OUT,"$outfile");
					while(<OUT>){
					#	print $_;
				
						print GLOBAL $_;
					}
					close OUT;
					if ($nb%$step == 0){
						#sleep(5);
						my $p = ($nb/$total);
						my $time1 = abs(time-$now);
						if ($time1 > 5){
						my $remaining = $eta->get_remaining_time();
						my $elapsed = $eta->get_elapsed_time();
						my $p = int($eta->get_completed_percent());
						print $patient->name." complete $p % => elapsed :".$elapsed." remaining : ".$remaining."\n";
						}
					}
					 $eta->pass_milestone();
				unlink $outfile.".idx" if -e  $outfile.".idx";
				}
		close GLOBAL;
	
	#	die();
	warn $final_gvcf;
system("$bgzip -f  $final_gvcf && $tabix -p vcf $gz");	
	confess() unless -e $gz;
 	confess() unless -e $gz.".tbi";
my $dir_prod  = $patient->project->getGvcfDir("haplotypecaller4");
#system("mv $gz $dir_prod ");
system("mv $gz* $dir_prod ");
warn "mv $gz* $dir_prod ";
my @t = `$bcftools query -l $gz`;
chomp(@t);
confess('BIG BIG AND REBIG PROBLEM ') if $t[0] ne $patient->name();
print " END ".$patient->name." join \n";
$project->destroy_tmp_dir();
exit(0);
 
sub check_callable {
	 			my ($chr,$window,$bam,$intspan0) = @_;
				  my @beds = $buffer->intspanToBed($chr,$window->{intspan});
				  	 my $id  = $window->{ext_gvcf};
				  my $bed1 =  calling_target::getTmpFile($dir_gvcf_out,$chr->name,"$id.bed");
 					open(BED,">$bed1");
 					print BED join("\n",@beds);
				 close BED;
				  my $bed =  calling_target::getTmpFile($dir_gvcf_out,$chr->name,"$id.callable.bed");
				  my $summary = $bed.".summary";
				  
				  my $cmd = qq{$java -jar $gatk -T CallableLoci -R $ref -I $bam -l off -o $bed  -summary $summary  --minDepth 10 --minDepthForLowMAPQ 15 -L $bed1};
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
 					unlink $bed;
 					unlink $bed1;
 					unlink $summary;
 				my $final_intspan = $intspan0->intersection($callable_intspan);
				 my $size  = scalar($final_intspan->as_array);
				 warn $size;
				 return $size > 100;
 					
}	
DESTROY {
	$project->destroy_tmp_dir();
}
