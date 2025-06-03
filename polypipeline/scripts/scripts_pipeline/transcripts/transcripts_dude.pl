#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use Parallel::ForkManager;
use Storable qw(store retrieve freeze );
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
my $nbErrors = 0;
my $buffer = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProject( -name => $project_name);



my $lists = $project->getListTranscripts();

my $transcript_event;
my $array;
my $hgene;
my $tmp1 = "/tmp/pipeline/";
system("mkdir /tmp/pipeline;chmod a+rwx  /tmp/pipeline") unless -e $tmp1;


my $tmp   =  File::Temp::tempdir ( TEMPLATE =>"$project_name.XXXXXXXX",DIR => "$tmp1" );
$project->get_only_list_patients($patient_name);
$project->getCaptures();
$project->getRuns();

#if ($patient_name eq "all" ){
 foreach my $patient (sort {$a->name cmp $b->name} @{$project->get_only_list_patients($patient_name)} ){
 		
 		my $hid = $patient->id;
 		$array->{$hid."_high"} = [];
 		$array->{$hid."_medium"} = [];
 		$array->{$hid."_low"} = [];
 		$hgene->{high} = {};
 		$hgene->{medium} = {};
 		$hgene->{low} = {};
 		uri_image($lists,$patient->name,$tmp);
 		
 }
 
exit(0);


	
	sub uri_image {
	my ($transcripts,$patient_name,$tmp_dir) = @_;
	my $patients = $project->getPatients();
	my $patient = $project->getPatient( $patient_name);
	#my $fork = 10;
	my $nb = int(scalar(@$transcripts)/($fork*1))+1;
	
	#$nb = 100;
	#$transcripts = [values %$genes];
	
	
	my $pm = new Parallel::ForkManager($fork);
	my $iter = natatime $nb, @$transcripts;
	my @t_final;

	my $images;
	my $first;
	my $hid = $patient->id;
	

	my $no = GenBoNoSqlLmdb->new(name=>$patient->name.".dude.transcripts",dir=>$tmp_dir,mode=>"c",is_compress=>1);
	$no-> put("date",time);
	$no->close;
	#$patient->getTranscriptsDude("d");
	$pm->run_on_finish(
    sub { 
    	my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$h)=@_;
  
    	unless (defined($h) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
			 my $no = GenBoNoSqlLmdb->new(name=>$patient->name.".dude.transcripts",dir=>$tmp_dir,mode=>"w",is_compress=>1);#$patient->getTranscriptsDude($type);
			$first =1;
			my $dir = $no->dir;
			
			foreach my $id (keys %{$h}){
				my $v = $h->{$id}->{event};
				$no->put($id,$h->{$id});
				if ($v >= 3){
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
		#	my $no = $chr->lmdb_image_transcripts_uri("w");
			$no->close();
		
		
		
    }
    );
  
	
	
  	$project->buffer->dbh_deconnect();
  	$|=1;
  	my $t =time;
  	#my $tid = 'ENST00000358544_8';
  	die() unless $project->isNoSqlDepth;
	my $first;
 	 while( my @tmp = $iter->() ){
 	 		my $pid = $pm->start and next;
 	 	
			#$project->buffer->dbh_reconnect();
			my $himages ={};
			my $znb =0;
			my $dj;
			my $res;
 	 	foreach my $tr1  ( @tmp){ 
 	 		my $transcript = $project->newTranscript($tr1);
 	 		$res->{$transcript->id} =  {};
 	 		
 	 		$res->{$transcript->id} = matrix_data($patient,$transcript);
 	 		#warn Dumper $res->{$transcript->id} if $transcript->id eq $tid;
 	 	}
 	 	
 	 	$pm->finish(0,$res);
	}
	$pm->wait_all_children();
	$project->buffer->dbh_reconnect();
	 my $no = GenBoNoSqlLmdb->new(name=>$patient->name.".dude.transcripts",dir=>$tmp_dir,mode=>"w",is_compress=>1);
	
	$no->put("high",$array->{$hid."_high"});		
	$no->put("medium",$array->{$hid."_med"});	
	$no->put("low",$array->{$hid."_low"});
	$no->close();
	sleep(10);
	my $f1 = $no->filename();
	my $no2 = $patient->getTranscriptsDude("c");
	my $f2 = $no2->filename();
	$no2->close();
	system("rsync -rav $f1 $f2; rm $f1;rmdir $tmp_dir");
	
}


sub between {
  my($test ,$fom, $tom)=@_;
  no warnings;
  $fom<$tom ? $test>=$fom && $test<=$tom
            : $test>=$tom && $test<=$fom;
}



sub matrix_data {
	my ($patient,$transcript) = @_;
	my $debug;
	$debug =1 if $transcript->id eq "ENST00000644074_X";
	#next unless $debug;
#	warn $transcript.' '.$transcript->id;
	my $chr = $transcript->getChromosome();
	my $no = $chr->get_lmdb_cnvs("r");
	my $primers = $project->getPrimersByPosition($transcript->getChromosome,$transcript->start,$transcript->end);#$transcript->getPrimers();
	my $levels = [];
	my $scores = [];
	my $zscores;
	my $strand = $transcript->strand;
	my $positions;
	my $alert =0;
	my $nb =0;
	my $names;
	foreach my $primer (sort{$a->end*$strand <=> $b->end *$strand} @$primers) {
		push(@$names,$primer->name);
		$alert ++ if $primer->cnv_score($patient) < 0.7 &&  $primer->cnv_score($patient) > -1;
		my $score = int(($primer->cnv_score($patient)+0.01)*100);
#		warn $score if $score < 0;
		$score = 0 if $score < 0 ;
		push( @$scores,int($score) );
		my $level = $primer->level($patient) +1;
		push( @$levels,$level);
		push( @$zscores,$primer->zscore($patient));
		warn $patient->name." ".$primer->id." level ".$primer->level($patient)." score ".$primer->cnv_score($patient)." save ".$score if $debug;
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
	warn Dumper unpack("w".scalar(@$scores),$h->{scores}) if $debug;
	$h->{event} =$max;
	$h->{level} = pack("w".scalar(@$levels),@$levels);
	#$h->{zscore} = pack("w".scalar(@$zscores),@$zscores);
	$h->{smooth_expo} = pack("w".scalar(@$data_smoothed1),@$data_smoothed1);
	$h->{smooth} = pack("w".scalar(@$data_smoothed2),@$data_smoothed2);
	$h->{nb} = scalar(@$scores);
	$h->{names} =$names;
	$h->{gene} = $transcript->gene_id;
	warn Dumper $h if $debug;
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
	
	