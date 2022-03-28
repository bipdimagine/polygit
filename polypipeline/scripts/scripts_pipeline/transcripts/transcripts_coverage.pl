#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use Parallel::ForkManager;
use Storable qw(store retrieve freeze thaw);
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use Cwd 'abs_path';
use Digest::MD5 qw(md5_hex);

use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
#use lib "/software/polyweb/poly-disk/poly-src/GenBo/lib/obj-nodb/packages";
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/packages";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use Set::IntSpan::Fast::XS;
use image_coverage;
 use List::Util qw( min sum);
 use File::Temp;
#use Cache_Commons;
use polyweb_coverage;


my $fork = 1;
my $force;
my ($project_name, $patient_name);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'patient=s'    => \$patient_name,
	'force=s'  => \$force,
	
);
my $globaltime = time;
unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
#unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }
my $t = `hostname`;
my $nbErrors = 0;
my $buffer = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProject( -name => $project_name);


my $lists = $project->getListTranscripts();

my $tmp1 = "/tmp/pipeline/";
system("mkdir /tmp/pipeline;chmod a+rwx  /tmp/pipeline") unless -e $tmp1;


my $tmp   =  File::Temp::tempdir ( TEMPLATE =>"$project_name.XXXXXXXX",DIR => "$tmp1" );

my $patients = $project->get_only_list_patients($patient_name);
foreach my $chr ( @{$project->getChromosomes()}){
	$chr->getIntSpanCapture();
	$chr->length;
	
}

foreach my $p (@{$project->get_only_list_patients($patient_name)}){	
	my $array ={};
	my $no2 = $p->getTranscriptsCoverageDepth("w");
	my $f2 = $no2->filename();
	if (-e $f2){
		warn $f2;
		my $t = $no2->get("exon_no_utr_0_30");
		next if $t;
	}
	
	$no2->close();
	uri_image($lists,$p->name,$tmp);

}
warn abs(time-$globaltime);
exit(0);


#uri_image($lists,$patient_name,$tmp);
#exit(0);

	sub uri_image {
	my ($transcripts,$patient_name,$tmp) = @_;
	my $patients = $project->getPatients();
	my $patient = $project->getPatient( $patient_name);
	#my $fork = 10;
	my $nb = int(scalar(@$transcripts)/($fork*1))+1;
	
	#$nb = 10;
	#$transcripts = [values %$genes];
	
	
	my $pm = new Parallel::ForkManager($fork);
	my $iter = natatime $nb, @$transcripts;
	my @t_final;
	my $images;
	my $first;
	
	my $no = GenBoNoSqlLmdb->new(name=>$patient->name.".transcripts",dir=>$tmp,mode=>"c",is_compress=>1);
	my $f1 = $no->filename;
	warn $f1;
	$no-> put("date",time);
	$no->close();
	

	$patient->getTranscriptsCoverageDepth("d");
		my $harray = {};
	$pm->run_on_finish(
    sub { 
    	my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$h)=@_;
  
    	unless (defined($h) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
			my $type = "w";
			#$type ="c" unless $first;
			#warn $type;
			#$first = 1; 
			my $no = GenBoNoSqlLmdb->new(name=>$patient->name.".transcripts",dir=>$tmp,mode=>"w",is_compress=>1);;#$patient->getTranscriptsCoverageDepth("w");
			my $dir = $no->dir;
			
			foreach my $id (keys %{$h}){
				$no->put($id,$h->{$id});
				#warn Dumper $h->{all_min_no_utr};
				my $z = $h->{$id};
				my $debug;
				$debug =1 if $id =~ /ENST00000368474/;
#				warn Dumper $z if $debug;
					my @scores = (1,10,20,30,50,100,200);
					foreach my $type (keys %{$z->{minimum}}) {
					my $v1 = $z->{minimum}->{$type};
							foreach my $sc  (@scores) {
								$harray->{$type."_$sc"} =[] unless exists $harray->{$type."_$sc"} ;
								push(@{$harray->{$type."_$sc"}},$id) if $v1 <= $sc;
					
						}
					}
					$no->put($id."_minimum", $z->{minimum});
					delete $z->{minimum};
					$no->put($id, $z);
					
			}
		#	my $no = $chr->lmdb_image_transcripts_uri("w");
			

		$no->close();
    }
    );
  
	$project->getRuns();
	$project->getCaptures();
  	$project->buffer->dbh_deconnect();
  	$|=1;
  	my $t =time;
  	
  	die("please run NosqlDepth ") unless $project->isNoSqlDepth;
	warn "start";
	
 	 while( my @tmp = $iter->() ){
 	 		my $pid = $pm->start and next;
 	 	
			#$project->buffer->dbh_reconnect();
			my $himages ={};
			my $znb =0;
			my $dj;
			my $res;
			my $cpt =0;
			my $t = time;
 	 	foreach my $tr1  ( @tmp){ 
 	 		$cpt++;
 	 		warn $cpt."/".scalar(@tmp)." ". abs (time -$t) if $cpt%100 ==0;
 	 		my $transcript = $project->newTranscript($tr1);
 	 		#$project->buffer->dbh_deconnect();
 	 		#next;
 	 		$res->{$transcript->id} =  {};
 	 		$res->{$transcript->id} = matrix_data($patient,$transcript);
 	 	}
 	 	
 	 	$pm->finish(0,$res);
	}
	$pm->wait_all_children();
	$project->buffer->dbh_reconnect();
	

	my $no = GenBoNoSqlLmdb->new(name=>$patient->name.".transcripts",dir=>$tmp,mode=>"w",is_compress=>1);
	foreach my $k (keys %$harray){
		
#			warn $k." ".scalar(@{$harray->{$k}});
			$no->put( $k ,$harray->{$k});
	}
	$no->put( "toto","titi");
	$no->close;
	my $f2  = $project->transcriptsCoverageDir()."/".$patient->name.".transcripts";
	warn $f2;
	unlink $f2 if -e  $f2;
	warn "$f1 $f2";
	system("mv $f1 $f2;rmdir $tmp");
	sleep(5);
	
warn "read";
	
	my $no  = $patient->getTranscriptsCoverageDepth("r");
	warn Dumper $no->get("exon_no_utr_0_30");
	warn Dumper $no->get("toto");
	}






sub matrix_data {
	my ($patient,$transcript) = @_;
	
	my $maskExons = {
	  exon => 1,
      intron => 2,
 	  non_coding => 4,
 	 capture =>8,
 	 duplicate => 16,
 	  "5_utr" => 32,
 	  "3_utr" => 64,
 	   "non_coding_transcript" => 128,
};
	my $debug;
	$debug =1 if $transcript->name eq "ENST00000368474";
	my $padding = 50;
	my $exons ;
	my $intronic=1;
	my $vstart = $transcript->start -30;
	my $vend = $transcript->end +30;
	my $array =  $patient->depth($transcript->getChromosome->name,$vstart,$vend);
	#my $string = pack("C".@$array,@$array);
#if ($intronic == 0 ){
	$exons  = $transcript->getAllGenomicsParts();
	my $strand =1;
	my $capture_intspan = $transcript->getChromosome->getIntSpanCapture()->intersection($transcript->getGenomicSpan());
	 	 my $no = $transcript->project->noSqlCoverage();
	 	my $spanTr = $no->get($patient->name,$transcript->id."_spandup");	

		my $data;
		my $data2;
		


		my $data_patient;
		my $data_patient_mean;
		my $non_coding_transcript == 0;
		$non_coding_transcript	 = 1  if $transcript->getSpanCoding()->is_empty;

			my $etypes;
				my $hexons;
				my $strand = $transcript->strand;
				my $nbe = -1;
			my $minimum_global;	
foreach my $exon (sort{$a->end*$strand <=> $b->end*$strand} @$exons) {
			$nbe ++;
			
	 		my $etype = 0;
			if ($exon->isExon){
				$etype = $etype | $maskExons->{"exon"} ;
				$etype = $etype | $maskExons->{"non_coding"}  if $exon->is_noncoding ;
				$etype = $etype | $maskExons->{"5_utr"}   if $exon->is_start_utr;
				$etype = $etype | $maskExons->{"3_utr"}   if $exon->is_end_utr;
			
			
			}
			else {
			$etype = $etype | $maskExons->{"intron"} ;
			}
			$etype = $etype | $maskExons->{"non_coding_transcript"}   if  $non_coding_transcript == 1;
			my $s1 = $exon->getGenomicSpan()->intersection($capture_intspan);
			$etype = $etype | $maskExons->{"capture"}  unless $s1->is_empty;
			if ($spanTr){	
			my $span =$exon->getGenomicSpan->intersection($spanTr);
			unless ($span->is_empty){
				$etype = $etype | $maskExons->{"duplicate"} ;
			}
			}
			push(@$etypes,$etype);
			my $hexon;
	
			#my $array =  $patient->depth($transcript->getChromosome->name,$exon->start-50,$exon->end+50);
			
			
			my $is_start_utr = $exon->is_start_utr;
			my $is_no_utr = $exon->utr->is_empty;
			
			my $is_end_utr = $exon->is_end_utr;
			my $span_no_utr = $exon->intspan_no_utr;
			my $sstart = $exon->start_utr;
			my $send = $exon->end_utr;
			#my $is_no_utr = $exon->intspan_no_utr->is_empty;
				my $i =40;
				my $a1 =  ($exon->start - $vstart)-$i;
	  			my $b1 = ($exon->end-$vstart) +$i;
				my @xb =@$array[$a1 .. $b1];
				
				my @xb2 =[];
				#array no utr 
			
				unless ($exon->utr->is_empty ){
					unless ($exon->is_noncoding){
				
						$a1 =  $exon->start_utr -$vstart;
						$b1 = $exon->end_utr - $vstart;
						$a1 -= $i unless $is_start_utr ;
						$b1 += $i unless $is_end_utr ;
				
					 	@xb2 =@$array[$a1 .. $b1];
					 	
					 	 unless (scalar @xb2){
					 	 	warn $a1." ".$b1;
					 	 	warn  $exon->start_utr."-".$exon->end_utr;
					 	 	
					 	 	warn $exon->start."-".$exon->end;
					 	 	warn $exon->intspan_no_utr->as_string();
					 	 	warn $exon->utr->as_string();
					 	 	confess();
					 	 }
					 	
					}
				}
				
					
				
				
		#		warn "d ".scalar(@xb);
			for (my $i=0;$i<=30;$i+=10){	
				#my $pos = $exon->return_start_end_no_utr(padding=>50);
				#my $pos1 = $exon->return_start_end(padding=>$padding);	
				#init tab 
				splice(@xb,0,10);
			#	warn scalar(@xb);
				splice(@xb,-10);
				#splice(@$xa,-10,);
				
				
			
	  		
	  			#my @xb =@$array[$a1 .. $b1];
	  			
	  			#warn scalar(@xb);
	   			
	   			my $minb = 0;
	   			my $minb =  min(@xb);
	   			my $value_min_global;
	   			
	   			#warn $transcript->name." ".$transcript->start." ->".$exon->start." ".$vstart;
	   			my $meanb = 0;
				my $meanb = int(sum(@xb)/scalar(@xb));
				push(@{$data_patient->{min_utr}->{$i}},$minb);
				push(@{$data_patient->{mean_utr}->{$i}},$meanb);
				if ( $etype & $maskExons->{"capture"}){
					
					$minimum_global->{"exon_utr_$i"} =swap_min($minimum_global->{"exon_utr_$i"},$minb);
				}
				
				if ( $etype & $maskExons->{"capture"}){
					$minimum_global->{"capture_utr_$i"} = swap_min($minimum_global->{"capture_utr_$i"},$minb);
				}
				$minimum_global->{"all_utr_$i"} = swap_min($minimum_global->{"all_utr_$i"},$minb);;
				 
				if ($exon->is_noncoding){
					
					push(@{$data_patient->{min_no_utr}->{$i}},0);
					push(@{$data_patient->{mean_no_utr}->{$i}},0);
				}
				elsif ($exon->utr->is_empty){
						push(@{$data_patient->{min_no_utr}->{$i}},$minb);
						push(@{$data_patient->{mean_no_utr}->{$i}},$meanb);
				}
				else{
					die() unless scalar @xb2;
					unless ($is_start_utr){
						splice(@xb2,0,10);
					}
					unless  ($is_end_utr){
						splice(@xb2,-10);
					}

					
					my $minb =  min(@xb2);
					my $meanb = int(sum(@xb2)/scalar(@xb2));
					
						if ( $etype & $maskExons->{"exon"}){
					$minimum_global->{"exon_no_utr_$i"} = swap_min( $minimum_global->{"exon_no_utr_$i"},$minb);;
				}
				
				if ( $etype & $maskExons->{"capture"}){
					$minimum_global->{"capture_no_utr_$i"} = swap_min( $minimum_global->{"exon_no_utr_$i"},$minb);
				}
				
					$minimum_global->{"all_no_utr_$i"} = swap_min($minimum_global->{"exon_no_utr_$i"},$minb);
					push(@{$data_patient->{min_no_utr}->{$i}},$minb);
					
					push(@{$data_patient->{mean_no_utr}->{$i}},$meanb);
					
				}
				
			}

			
			push(@{$hexons->{name}},$exon->name);
			
	
	}
	$hexons->{nb} = scalar(@$exons);
	#return {};
	#warn $transcript->name;
	#die();
	my $hmin = -2;
	
	
	for (my $i=0;$i<=30;$i +=10){
		
		#	warn $i;	
	
		$hexons->{min_no_utr}->{$i}   = pack("w".scalar(@$exons),@{$data_patient->{min_no_utr}->{$i}});
		
		$hexons->{all_min_no_utr}->{$i} = min(@{$data_patient->{min_no_utr}->{$i}});
		
		$hexons->{min_capture_no_utr}->{$i} = min(@{$data_patient->{min_capture_no_utr}->{$i}});
		
		$hexons->{min_capture_utr}->{$i} = min(@{$data_patient->{min_capture_utr}->{$i}});
		warn $i." ".$hexons->{min_capture_no_utr}->{$i} if $debug;
		$hexons->{mean_no_utr}->{$i}   = pack("w".scalar(@$exons),@{$data_patient->{mean_no_utr}->{$i}});
		$hexons->{min_utr}->{$i}   = pack("w".scalar(@$exons),@{$data_patient->{min_utr}->{$i}});
		#$hexons->{all_min_utr}->{$i} = min(@{$data_patient->{min_utr}->{$i}});
		warn $i." ".$hexons->{all_min_utr}->{$i} if $debug;
		$hexons->{mean_utr}->{$i}   = pack("w".scalar(@$exons),@{$data_patient->{mean_utr}->{$i}});
	}
	
	$hexons->{types} = pack("w".scalar(@$exons),@{$etypes});
	$hexons->{minimum}= $minimum_global;
	return $hexons;
#return $data;	
}

sub swap_min {
	my ($v1,$v2) = @_;
	return $v2 unless $v1;
	return $v1 if $v1 < $v2;
	return $v2;
	
}
	
	