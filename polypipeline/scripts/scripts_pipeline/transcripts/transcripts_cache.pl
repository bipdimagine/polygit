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
use Statistics::Descriptive::Smoother;
use File::Temp;

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
 use List::Util qw( min sum max);
#use Cache_Commons;

my $fork = 1;
my $force;
my ($project_name, $patient_name);

GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'patient=s'    => \$patient_name,
	'force=s'  => \$force,
);



unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
#unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }
my $t = `hostname`;
warn "****  ".$t." ******";
my $nbErrors = 0;
my $buffer = new GBuffer;
$buffer->vmtouch(1);
warn "$t";
my $project = $buffer->newProjectCache( -name => $project_name);
warn "end";
#	 my $no = $project->noSqlCoverage();
#	 warn Dumper $no;
#	 die ();
my $lists = $project->getListTranscripts();
warn Dumper $lists;
my $tmp1 = "/tmp/pipeline/";
system("mkdir /tmp/pipeline;chmod a+rwx  /tmp/pipeline") unless -e $tmp1;

my $tmp   =  File::Temp::tempdir ( TEMPLATE =>"$project_name.XXXXXXXX",DIR => "$tmp1" );
#$tmp = $tmp1."/$project_name.".time*rand(time);
my $patients = $project->getPatient($patient_name);
$project->getCaptures();
$project->getRuns();
#my $dir1 = "/tmp/".$project->name();
#system("mkdir  -p $dir1;chmod a+rwx $dir1") unless -e $dir1;
my $dir2 = $project->noSqlCnvsDir();
system("rsync -av  $dir2/*  /$tmp/");
#For coverage prepare 
foreach my $chr ( @{$project->getChromosomes()}){
	$chr->getIntSpanCapture();
	$chr->length;
	
}
$project->noSqlCnvsDir($tmp);
warn "end";	
my $array;
my $hgene;


foreach my $p (@{$project->get_only_list_patients($patient_name)}){	
	my $array ={};
	my $dir3 = $p->NoSqlDepthDir;
	my $no2 = $p->getTranscriptsCoverageDepth("w");
	my $f2 = $no2->filename();
	system("touch $tmp/* ");
	if (-e $f2){
		warn $f2;
		my $t = $no2->get("exon_no_utr_0_30");
#		next if $t;
	}
	
	$no2->close();
	
	system("rsync -av ".$dir3."/".$p->name.".depth.lmdb $tmp/");
	$p->NoSqlDepthDir($tmp);
	my $hid = $p->id;
 		$array->{$hid."_high"} = [];
 		$array->{$hid."_medium"} = [];
 		$array->{$hid."_low"} = [];
 		$hgene->{high} = {};
 		$hgene->{medium} = {};
 		$hgene->{low} = {};
 		uri_image($lists,$p->name,$tmp);

}
system("rm -r $tmp") if $tmp=~/tmp/;

#warn abs(time-$globaltime);
exit(0);



sub run_on_finish_coverage {
	my ($patient,$h,$harray) =@_;
	
			#$type ="c" unless $first;
			#warn $type;
			#$first = 1; 
			warn $tmp;
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

sub run_on_finish_dude {
	
	my ($patient,$h) = @_;
	my $hid = $patient->id;
	my $no = GenBoNoSqlLmdb->new(name=>$patient->name.".dude.transcripts",dir=>$tmp,mode=>"w",is_compress=>1);#$patient->getTranscriptsDude($type);
			my $dir = $no->dir;
			
			foreach my $id (keys %{$h}){
				warn $id;
				my $v = $h->{$id}->{event};

				$no->put($id,$h->{$id});
				if ($v >= 3){
					warn "******************* CPOUOUOU HIGH ".$v." ".$id if $id =~ /ENST00000314358/;
					push (@{$array->{$hid."_high"}},$id) ;
					$hgene->{high}->{$hid."_".$h->{gene}}++;
				}
				if ($v >= 2){
					push (@{$array->{$hid."_med"}},$id) ;
					$hgene->{med}->{$hid."_".$h->{gene}}++;
				}
				if ($v >= 2){
					push (@{$array->{$hid."_low"}},$id) ;
					$hgene->{low}->{$hid."_".$h->{gene}}++;
				}
 				push (@{$array->{$hid."_med"}},$id) if ($v >= 2);
 				push (@{$array->{$hid."_low"}},$id) if ($v  >= 1);
			
			}
			#warn Dumper $array;
		#	my $no = $chr->lmdb_image_transcripts_uri("w");
			$no->close();
	return 1;
}



sub uri_image {
	my ($transcripts,$patient_name,$tmp) = @_;
	my $patients = $project->getPatients();
	my $patient = $project->getPatient( $patient_name);
	#my $fork = 10;
	my $nb = int(scalar(@$transcripts)/(100*1))+1;
	#cache declaration 
	my $no_coverage = GenBoNoSqlLmdb->new(name=>$patient->name.".transcripts",dir=>$tmp,mode=>"c",is_compress=>1);
	my $f1 = $no_coverage->filename;
	warn $f1;
	$no_coverage-> put("date",time);
	$no_coverage->close();
	$patient->getTranscriptsCoverageDepth("d");		
	
	#dude declaration 
	my $no_dude = GenBoNoSqlLmdb->new(name=>$patient->name.".dude.transcripts",dir=>$tmp,mode=>"c",is_compress=>1);
	$no_dude-> put("date",time);
	$no_dude->close;
	
	
	my $process;
	
	my $pm = new Parallel::ForkManager($fork);
	my $iter = natatime $nb, @$transcripts;
	
	my $harray = {};
	$pm->run_on_finish(
    sub { 
    	my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$h)=@_;
  
    	unless (defined($h) or $exit_code > 0) {
    			warn Dumper $h;
				print qq|No message received from child process $exit_code $pid!\n|;
				warn Dumper $process;
				return;
			}

			run_on_finish_dude($patient,$h->{dude});
			run_on_finish_coverage($patient,$h->{coverage},$harray);
			warn $h->{end_process};
		delete $process->{$h->{end_process}};
    }
    );
	
	
	$project->getRuns();
	$project->getCaptures();
  	$project->buffer->dbh_deconnect();
  	$|=1;
  	my $t =time;
  	my $id = 0;
  	die("please run NosqlDepth ") unless $project->isNoSqlDepth;
	 my $no = $project->noSqlCoverage("w");
	 $no->put($patient->name,"date",time);
	 $no->close;
	
 	 while( my @tmp = $iter->() ){
 	 	 		$id ++;
 	 		 $process->{$id} = 1;
 	 		# my ($find) = grep {$_ =~ /ENST00000382051/} @tmp;
 	 		# next unless $find;
 	 		my $pid = $pm->start and next;
 	 		
			#$project->buffer->dbh_reconnect();
			my $himages ={};
			my $znb =0;
			my $dj;
			my $res;
			my $cpt =0;
			my $t = time;
			
			my $ts = $project->newTranscripts(\@tmp);
			
			my $z;
			
#			warn "get Transcripts :".abs($t -time)." ".scalar(@$ts); 
			$project->lmdbGenBo->close();
			delete $project->{lmdbGenBo};
			my $h;
			 $t = time;
			 my $xx = 0;
			 my $t0 = time; 
			warn "start  pp$id $pid ";
		foreach my $transcript (@$ts){
			warn $transcript->id;
			#warn $xx."/".scalar(@$ts)." ".abs($t0 -time) if $xx % 1000 ==0;
			$t0 = time if $xx % 1000 ==0;
			$xx++;
 	 		$res->{coverage}->{$transcript->id} =  {};
 	 		$res->{coverage}->{$transcript->id} = matrix_data_coverage($patient,$transcript,$id);
 	 		$res->{dude}->{$transcript->id} = matrix_data_dude($patient,$transcript,$id);
 	 		
 	 		
 	 	}
 	 	warn "finish ".abs($t -time)." pp$id $pid";
 	 	$res->{end_process} = $id;
 	 	$pm->finish(0,$res);
	}
	$pm->wait_all_children();
	
	$project->buffer->dbh_reconnect();
	die(Dumper $process) if scalar %$process;
	warn "end fork !!!!";
	end_coverage($patient,$harray);
	warn "end coverage !!!!";
	end_dude($patient);
	
	
	
}


sub end_coverage {
	my ($patient,$harray) = @_;
	my $no = GenBoNoSqlLmdb->new(name=>$patient->name.".transcripts",dir=>$tmp,mode=>"w",is_compress=>1);
	my $f1 = $no->filename;
	foreach my $k (keys %$harray){
		

			$no->put( $k ,$harray->{$k});
	}
	$no->put( "toto","titi");
	$no->close;
	my $f2  = $project->transcriptsCoverageDir()."/".$patient->name.".transcripts";
	unlink $f2 if -e  $f2;
	warn "$f1 $f2";
	system("mv $f1 $f2;");
	#sleep(5);
#	my $no1  = $patient->getTranscriptsCoverageDepth("r");
	#warn Dumper $no1->get("exon_no_utr_0_30");
	#warn Dumper $no1->get("toto");
#	$no1->close;
}

sub end_dude {
	my ($patient) = @_;
	 my $no = GenBoNoSqlLmdb->new(name=>$patient->name.".dude.transcripts",dir=>$tmp,mode=>"w",is_compress=>1);
	my $hid = $patient->id;
	warn Dumper $array;
	$no->put("high",$array->{$hid."_high"});		
	$no->put("medium",$array->{$hid."_med"});	
	$no->put("low",$array->{$hid."_low"});
	$no->close();
	sleep(10);
	my $f1 = $no->filename();
	my $no2 = $patient->getTranscriptsDude("c");
	my $f2 = $no2->filename();
	$no2->close();
	warn "rsync -rav $f1 $f2; rm $f1;";
	system("rsync -rav $f1 $f2; rm $f1;rmdir $tmp");
}

sub matrix_data_coverage {
	my ($patient,$transcript,$id) = @_;
	
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
	#$debug =1 if $transcript->name eq "ENST00000368474";
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
		$hexons->{mean_no_utr}->{$i}   = pack("w".scalar(@$exons),@{$data_patient->{mean_no_utr}->{$i}});
		$hexons->{min_utr}->{$i}   = pack("w".scalar(@$exons),@{$data_patient->{min_utr}->{$i}});
		#$hexons->{all_min_utr}->{$i} = min(@{$data_patient->{min_utr}->{$i}});
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





sub between {
  my($test ,$fom, $tom)=@_;
  no warnings;
  $fom<$tom ? $test>=$fom && $test<=$tom
            : $test>=$tom && $test<=$fom;
}




sub matrix_data_dude {
	my ($patient,$transcript,$id) = @_;
	my $debug;
	#next unless $debug;
#	warn $transcript.' '.$transcript->id;
	$debug =1;
	my $chr = $transcript->getChromosome();
	my $no = $chr->get_lmdb_cnvs("r");
	my $res = $project->tabix_primers->query($transcript->getChromosome->ucsc_name,$transcript->start,$transcript->end);
	#warn Dumper $res;
	#return {};
	my $primers = $project->getPrimersByPosition($transcript->getChromosome,$transcript->start,$transcript->end);#$transcript->getPrimers();
	#warn scalar(@$primers);
	#return {};
	#warn scalar(@$primers);
	my $levels = [];
	my $scores = [];
	my $zscores;
	my $strand = $transcript->strand;
	my $positions;
	my $alert =0;
	my $nb =0;
	my $names;
	foreach my $primer (sort{$a->end*$strand <=> $b->end *$strand} @$primers) {
	#	next unless $transcript->getGenomicSpan->contains_all_range( $primer->start, $primer->end);
		push(@$names,$primer->name);
		$alert ++ if $primer->cnv_score($patient) < 0.7 &&  $primer->cnv_score($patient) > -1;
		
		my $score = int(($primer->cnv_score($patient)+0.01)*100);
#		warn $score if $score < 0;
		$score = 0 if $score < 0 ;
		push( @$scores,int($score) );
		my $level = $primer->level($patient) +1;
		
		push( @$levels,$level);
		#push( @$zscores,$primer->zscore($patient));
#		warn $patient->name." ".$primer->id." level ".$primer->level($patient)." score ".$primer->cnv_score($patient)." save ".$score if $debug;
		push( @$positions,$primer->start);
		
		
			$primer = undef;
		
	}
	#warn Dumper $scores if $debug;
	
	my $data_smoothed1 = smooth2($scores,$positions,$alert)	;
	my $data_smoothed2 = smooth_data($scores);
	my $error;
	my $max = 0;
	foreach (my $i=0;$i<@$levels;$i++){
		next if $levels->[$i] == 0 ;
	#	warn $i." ".$levels->[$i];
			if ($levels->[$i]<= 1) {
				
				if  ($levels->[$i] == 0 or between($data_smoothed1->[$i],70,140) ) {
					$max = $error if $error > $max;
					$error =0;
					next;
				} 
			}
			$error ++;
			$max = $error if $error > $max;
	}
	#$max = 99 if ($max > 0.75 * @$levels);
	my $h;
	$h->{scores} = pack("w".scalar(@$scores),@$scores);
	#warn Dumper unpack("w".scalar(@$scores),$h->{scores}) if $debug;
	$h->{event} =$max;
	
	$h->{level} = pack("w".scalar(@$levels),@$levels);
	#$h->{zscore} = pack("w".scalar(@$zscores),@$zscores);
	$h->{smooth_expo} = pack("w".scalar(@$data_smoothed1),@$data_smoothed1);
	$h->{smooth} = pack("w".scalar(@$data_smoothed2),@$data_smoothed2);
	$h->{nb} = scalar(@$scores);
	$h->{names} =$names;
	$h->{gene} = $transcript->gene_id;
	#warn Dumper $h if $debug;
	return $h;
}


sub smooth2 {
		my ($data,$positions,$alert) = @_;
	
		
	
		
		
		return $data if scalar(@$data) <5;
		my @data2;
		foreach my $v (@$data){
			$v = 100 if $v == 0;
			push(@data2,$v);
		}
		
		my $smoother = Statistics::Descriptive::Smoother->instantiate({
         method   => 'weightedexponential',
         coeff    => 0.5,
         data     => \@data2,
         samples  => $positions,
           });
my @smoothed_data = map{int($_)} $smoother->get_smoothed_data();
	return \@smoothed_data;

}
sub find_previous {
	my ($pos,$data) = @_;
	return ($pos,undef) if $pos-1 == 0;
	#warn $pos unless $data->[$pos-1] ;
	return ($pos-1,$data->[$pos-1]) if $data->[$pos-1] ne -1;
	return find_previous($pos-1,$data);
	
}
sub find_next {
	my ($pos,$data) = @_;
	return ($pos,undef) if $pos+1 >= scalar(@$data);
	return ($pos+1,$data->[$pos+1]) if $data->[$pos+1] ne -1;
	return find_previous($pos+1,$data);
	
}
sub smooth_data {
	my ($data) = shift;
#	return $data;
	my @data2;
		for (my $i = 0;$i<@$data;$i++){
		my @trio;
		my ($pos,$val) = find_previous($i,$data);
	
		if ($val){
				 push (@trio,$val);
				 	my ($pos2,$val2) = find_previous($pos,$data);
				 	 push (@trio,$val2) if $val2;
		}
		my ($posn,$valn) = find_next($i,$data);
		
		if ($valn){
			 push (@trio,$valn);
			 my ($pos2,$val2) = find_previous($posn,$data);
			 push (@trio,$val2) if $val2;
		}
			 push (@trio,$data->[$i]);
			push (@trio,$data->[$i]);
		
		  my $z = sum @trio;
		  push(@data2, int($z/(scalar(@trio) ) ) ) ;
		 
	}
	return \@data2;
}
	
	
	