#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use Set::IntSpan::Island;
use Set::IntSpan::Fast::XS;
use Array::IntSpan;
use lib $Bin;
use GBuffer;

use IPC::Open2;
use Data::Dumper;
use Getopt::Long;
use Carp;
use Scalar::Util qw(looks_like_number);
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use colored;
use Bio::DB::Sam;
use Parallel::ForkManager;
use String::ProgressBar;
use Set::IntSpan::Fast::XS;
#use Bio::DB::HTS;
 use JSON::XS;
use List::Util  qw(sum);
use IO::Handle;

my $filein;
my $dir;
my $file_bed;
my $name;
my $fork = 1;
my $project_name;
my $patient_name;
my $verbose;
my $use_samtools;
my $log_file;
GetOptions(
	'filein=s' => \$filein,
	'dir=s'    => \$dir,
	'name=s'   => \$name,
	'bed=s'    => \$file_bed,
	"fork=s"   => \$fork,
	"project=s" =>\$project_name,
	"patient=s" =>\$patient_name,
	"verbose=i" =>\$verbose,
	"use_samtools=i" =>\$use_samtools,
	"log=s" =>\$log_file,
);
unless($patient_name){
	$patient_name = $name;
}
unless ($name){
	$name = $patient_name;
}

if ($log_file){
	open (STDOUT,">>".$log_file);
}

my @chrs;# = (1..22, 'X', 'Y','MT' );
#my $file_bed = "/data-xfs/public-data/HG19/capture/agilent/agilent.v50.bed";
#open( BED, "zcat $file_bed |" );

colored::stabilo('cyan',"  ---- start $patient_name ----");

$SIG{INT} = \&interrupt;
$use_samtools = undef;


#die();
#
my $span;
my $span_extended;

my $buffer = GBuffer->new();
my $project = $buffer->newProject(-name=>$project_name);

my $dir_temp =  $buffer->config_path("tmp");
#my $patient = $project->get_only_list_patients($name);

foreach my $chr (@{$project->getChromosomes()}) {
	#my $chr = $project->getChromosome($chr_name);
	my $spant = $chr->getCapturesGenomicSpan();
	next if $spant->is_empty();
	push(@chrs,$chr->name);
	$span->{$chr->ucsc_name} =$spant;
	my $span_ext =  $chr->getIntSpanCaptureForCalling(250);#$chr->getExtentedGenomicSpan(250);

	$span_extended->{$chr->ucsc_name} = $span_ext;
}

my $diag;

my $samtools = $buffer->config->{software}->{samtools};
my $bgzip = $buffer->config->{software}->{bgzip};
my $tabix = $buffer->config->{software}->{tabix};
#$project->getPatientsAndControl()
my $patient = $project->getPatientOrControl($patient_name);
my $capture = $patient->getCapture();
my $transcripts_cgi = $project->bundle_transcripts() ;
#$diag =1 if $transcripts_cgi;
  $diag=1 if $project->isDiagnostic();

my $pm = new Parallel::ForkManager($fork);
my $pr = String::ProgressBar->new( max => scalar(@chrs) );
my $c =0;
$| = 1;
my $res;
my $nb_5;
my $nb_15;
my $nb_20;
my $nb_30;
my $nb_50;
my $nb_100;
my $nb_all;
my $tot;
$pm->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$h)=@_;
  	$nb_all += $h->{nb};
	$nb_5 += $h->{5}; 
	$nb_15 +=$h->{15};
	$nb_20 +=$h->{20};
	$nb_30 +=$h->{30};
	$nb_50 +=$h->{50};
	$nb_100 +=$h->{100};
	$tot += $h->{99};
      
    }
    );
  
  my %temp_file;
foreach my $chr (@chrs){
	 $temp_file{$chr} =  File::Temp->new( TEMPLATE =>"$name.$chr.".'XXXXXX',
                        DIR => $buffer->config_path("tmp"),
                        SUFFIX => '.cov');
                        
	
}   
   
$filein = $patient->getBamFile() unless $filein;
  
  $project->disconnect();
 foreach my $chr (@chrs){
 		my $pid = $pm->start and next;
 		
	my $h;
	if ($use_samtools){
		by_chr_samtools($chr,$span, $span_extended, $name, $filein,  $temp_file{$chr} );
	}
	else {
		 $h = by_chr_samtoolsfast($project,$chr,$span, $span_extended, $name, $filein, $temp_file{$chr} );
		#by_chr2($chr,$span, $span_extended, $name, $filein, $temp_file{$chr} );
	} 
	$pm->finish(0,$h);
}
$pm->wait_all_children();
colored::stabilo('magenta',"  ---- end phase 1  $patient_name ----");

warn " --------- END RAW COVERAGE------------";
$project->buffer->dbh_reconnect();

my $nb;
$dir = $project->getCoverageDir() unless $dir;
my $output   = $dir. "/$name.cov";
my $outputgz = $output . ".gz";
open( TOTO, ">$output" ) or die("can't open $output") ;
system("chmod a+w $output");
my @limit = ( 1, 5, 10, 15, 20, 50 );

#

foreach my $chr_name (@chrs) {
	my $chr = $project->getChromosome($chr_name);
	my $file =$temp_file{$chr->name}->filename;
	my $ucsc =$chr->ucsc_name;
	unlink $file unless exists $span_extended->{$ucsc}; 
	next unless exists $span_extended->{$ucsc};


	open( TITI, $file );
	while ( my $line = <TITI> ) {
		print TOTO $line;
	}
	close TITI;
	unlink $file;
}
close TOTO;

run("$bgzip -f $output;");
system("chmod a+w $outputgz");
run(" $tabix -f -s 1 -b 2 -e 2  $outputgz");

if ($@){
	die("tabix or bgzip ");
}
my $t = time;
#colored::stabilo('yellow',"  ---- start primers  $patient_name ----");
#primers($project_name,$patient_name);
#my ($project_name,$patient_name) = @_;

#colored::stabilo('yellow',"  *** end primers  $patient_name ***".int((time -$t)/60) );
if ($diag){
	warn "COMPUTE COVERAGE FOR POLYDIAG ";
	warn "compute transcripts;\n-*-*-*-*-*- you have ".scalar(@$transcripts_cgi)." in your design -*-*-*-*";
	$buffer = undef;
	$project = undef;
	$nb_all =0;
	$nb_5 =0;
	$nb_15 =0;
	$nb_20 =0;
	$nb_50=0;
	$nb_30 =0;
	$nb_100 =0;
	$tot =0;
 	my $buffer2 = GBuffer->new();
 	my $project2 = $buffer2->newProject(-name=>$project_name);
	 my $patient = $project2->getPatientOrControl($patient_name);
	  
	  my $hrequest;
	  	  foreach my $tr (@$transcripts_cgi){
	  	  	my $tr1 = $project2->newTranscript($tr);
    			my $gene = $tr1->getGene();
    			push(@{$hrequest->{$gene->id}->{transcripts}},$tr);
    			$hrequest->{$gene->id}->{start} = $gene->start;
    			$hrequest->{$gene->id}->{end} = $gene->end;
    			$hrequest->{$gene->id}->{chr} = $gene->getChromosome()->name;
    			
	  	  }

my $pm = new Parallel::ForkManager($fork);
 my $db2 = $patient->transcriptsCoverageLite('c');
$pm->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$h)=@_;
  		
      $db2->set($patient->name,$h->{gene_id}."_coverage_json",$h->{json});
      foreach my $t (keys %{$h->{transcripts}}){
      	$db2->set($patient->name,$t."_mean_exonic", $h->{transcripts}->{$t}->{mean_exonic});
    		$db2->set($patient->name,$t."_mean_coding", $h->{transcripts}->{$t}->{mean_coding});
      }
      		$nb_all +=$h->{nb_all} ;
	 		$tot += $h->{tot} ;
	 		$nb_5+=$h->{nb_5};
	 		$nb_15+=$h->{nb_15};
	 		$nb_20+=$h->{nb_20};
	 		$nb_30+=$h->{nb_30};
	 		$nb_50+=$h->{nb_50};
	 		$nb_100+=$h->{nb_100};
    }
    );

$patient->project->disconnect;
warn "disconnect -----------";
foreach my $g (keys %$hrequest){
	my $gene_id = $g;
	my $pid = $pm->start and next;
		my $buffer = GBuffer->new();
 		my $project = $buffer->newProject(-name=>$project_name);
	 	my $patient = $project->getPatientOrControl($patient_name);
	 	my $chromosome = $project->getChromosome($hrequest->{$g}->{chr});
	 	my $start = $hrequest->{$g}->{start};
	 	my $end = $hrequest->{$g}->{end};
		my $gc = GenBoCoverageTabix->new(chromosome=>$chromosome, patient=>$patient, start=>$start-100, end=>$end+100);
		my $trs = $hrequest->{$g}->{transcripts};
		my $h;
		
		  foreach my $tr (@$trs){
	  	  	my $tr1 = $project->newTranscript($tr);
	  	  	my $gc_exonic = $gc->coverage_intspan($tr1->getGenomicSpan);
    			my $gc_coding = $gc->coverage_intspan($tr1->getSpanCoding);
    			$h->{transcripts}->{$tr1->id} ->{mean_exonic} = $gc_exonic->{mean};
    			$h->{transcripts}->{$tr1->id} ->{mean_coding} = $gc_coding->{mean};
    			foreach my $z (@{$gc_coding->{array}}){
	 		$h->{nb_all} +=1;
	 		$h->{tot} += $z;
	 		$h->{nb_5}++ if $z>5;
	 		$h->{nb_15}++ if $z>15;
	 		$h->{nb_20}++ if $z>20;
	 		$h->{nb_30} ++ if $z>30;
	 		$h->{nb_50} ++ if $z>50;
	 		$h->{nb_100} ++ if $z>100;
	 	}
    			
		  }
		  $h->{json} = encode_json {start=>$gc->start,end=>$gc->end,array=>$gc->array};
		  $h->{gene_id} = $g;
	
	$pm->finish(0,$h);
} 
$pm->wait_all_children();	
}
run("gunzip  $outputgz");
open( TOTO, ">>$output" );
$nb_all =1 if $nb_all == 0;
my $mean    = $tot / $nb_all;
my $mean_5  = $nb_5 / $nb_all;
my $mean_15 = $nb_15 / $nb_all;
my $mean_20 = $nb_20 / $nb_all;
my $mean_30 = $nb_30 / $nb_all;
my $mean_50 = $nb_50 / $nb_all;
my $mean_100 = $nb_100 / $nb_all;
print TOTO "mean_all\t1\t$nb_all\n";
print TOTO "mean_all\t5\t$mean_5\n";
print TOTO "mean_all\t15\t$mean_15\n";
print TOTO "mean_all\t20\t$mean_20\n";
print TOTO "mean_all\t30\t$mean_30\n";
print TOTO "mean_all\t50\t$mean_50\n";
print TOTO "mean_all\t99\t$mean\n";
print TOTO "mean_all\t100\t$mean_100\n";
close TOTO;

run("$bgzip -f $output;");
system("chmod a+w $outputgz");

run(" $tabix -f  -s 1 -b 2 -e 2  $outputgz");
my $buffer = GBuffer->new();
 		my $project = $buffer->newProject(-name=>$project_name);
	 	my $patient = $project->getPatientOrControl($patient_name);
	warn Dumper $patient->coverage;
exit(0);

sub primers {
	my ($project_name,$patient_name) = @_;
	my $buffer2 = GBuffer->new();
 	my $project2 = $buffer2->newProject(-name=>$project_name);
	my $patient2 = $project2->getPatient($patient_name);
	my $dir_out = $project2->project_pipeline_path()."/coverage/";
	system("mkdir -p $dir_out") unless -e $dir_out;
	my $all_primers;
	my $c =0;
	
	my $chr_primer;
	
	my $samtools = $buffer2->software("samtools");
	
	
my $pm = Parallel::ForkManager->new($fork);
my %total;
$pm->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    		
    			foreach my $k (keys %$data){
    				$total{$k} = $data->{$k};
    			}
    		}
		    		
    		
  );

my $bam =  $patient2->getBamFile();
my %temp_bed;
my %regions;
foreach my $chr (@{$project2->getChromosomes}){
		warn "bed file ".$chr->name();
		
		my $bed_file =$buffer2->config_path("tmp")."/test.".$chr->name.".$patient_name.bed";
		 my $file_cat_tmp =  File::Temp->new( TEMPLATE =>"TMP.XXXXXXX",
                        DIR => $dir_out,
                        SUFFIX => ".bed");
          my $bed_file =  $file_cat_tmp->filename();
		$temp_bed{$chr->fasta_name} =  $file_cat_tmp;

	open(BED,">$bed_file");
	my $find =0;
	foreach my $p (sort {$a->start <=> $b->start} @{$chr->getPrimers}){
		print BED $p->getChromosome->fasta_name."\t".$p->start."\t".$p->end."\t".$p->id."\n";
		my $len = $p->start - $p->end;
		$regions{$p->id} = $p->getChromosome->fasta_name.";".$p->start.";".$p->end;
		$find ++;
	}
	
	close BED;
	delete $temp_bed{$chr->fasta_name} if $find ==0;
}
#my @toto = keys %regions;
#
#my $nb = scalar(@toto);
#my $zz = 0;
#foreach my $id (@toto) {
#	print "$zz/$nb\r";
#	$zz++;
#	my $pid = $pm->start and next;
#	 my ($chr1,$start1,$end1) = split(";",$regions{$id});
#	 
#	 my %res;
#	 my ($gc) = $buffer2->coverage_samtools($bam,$chr1,$start1,$end1);
#	 my $sum = $gc->sum;
#	 $res{$id."_depth"} = $sum/abs($end1-$start1);
#	$pm->finish(0,\%res);
#} 
#$pm->wait_all_children();
#die();
my @ucsc_names = map {$_->fasta_name} @{$project2->getChromosomes};
	foreach my $chr_name (@ucsc_names){
		warn "start " .$chr_name;
		next unless exists $temp_bed{$chr_name};
			warn $temp_bed{$chr_name};
			next unless -e $temp_bed{$chr_name}->filename();
			
			my $pid = $pm->start and next;
			warn "start";
		my $bed_file =	$temp_bed{$chr_name}->filename();
				warn "$samtools bedcov $bed_file $bam";
	open (OUT , "$samtools bedcov $bed_file $bam | ");

	my %res;
	while(my $line = <OUT>){
		chomp($line);
		my($chr,$start,$end,$id,$sum) = split(" ",$line);
		#$res{$id."_count"} = $sum;
		$res{$id."_depth"} = $sum/abs($end-$start);
	}
	close OUT;
	unlink $bed_file;
	colored::stabilo('white',"  ---- end   $chr_name -  $patient_name ----");
	$pm->finish(0,\%res);
	}
	warn "WAIT";
	$pm->wait_all_children();
	
	
	my $db1 = $patient2->transcriptsCoverage('c');
		foreach my $k (keys %total){
    					$db1->set($k,$total{$k});
    			}
	$db1->close();
}





sub by_chr2 {
my ($chr_num,$span,$span_extended,$name,$filein,$filetemp) = @_;
my $bam_file = $filein;
my $t = time;
my $sam = Bio::DB::Sam->new(-bam  =>$bam_file,
                            # -fasta=>"data/ex1.fa",
                             );
       my @seq_ids  = $sam->seq_ids;
      my $is_ucsc; 
    $is_ucsc =1    if ($seq_ids[0] =~ /chr/);

	my $array_intspan =  Array::IntSpan->new();
	my $out = $filetemp->filename;
	my $chr = $chr_num;
	 $chr = "chr".$chr_num  if ischrornot($filein);
	$chr="chrM" if $chr eq "chrMT";
		my $chr_ucsc = "chr" . $chr_num;
		$chr_ucsc ="chrM" if $chr_ucsc eq "chrMT";
	open (TITI,">$out",) or die("can't open $out");
	my $nb;
	my $tot;
	my $nb_5;
	my $nb_15;
	my $nb_20;
	my $nb_30;
	my $nb_100;
	return unless exists $span_extended->{$chr_ucsc};
	my $max;
	my $t = time;
	$sam->max_pileup_cnt([500000]) if $diag ==1;
foreach my $pos (split(",",$span_extended->{$chr_ucsc}->as_string)){	 
	my ($start,$end) = split("-",$pos);		
	$max = $end;
	my $coverage;
	if ($is_ucsc){
		($coverage) = $sam->features(-type=>'coverage',-seq_id=>$chr_ucsc,-start=>$start,-end=>$end-1  );
	}
	else {
		($coverage) = $sam->features(-type=>'coverage',-seq_id=>$chr,-start=>$start,-end=>$end-1  );
	}
	my $xstart ;
	my $val = -1;
	my $debug;
	$debug=1 if $chr eq "chr17";
	for (my $i=0;$i< @{$coverage->coverage};$i++){
		
		my $pos = $i+$start;	
		my $cov = 	$coverage->coverage->[$i];
		unless ($xstart){
			$xstart = $pos;
			$val = $cov;
		}
		if ($cov ne $val) {
	#		$array_intspan->set_range($xstart,$pos-1,$cov);
			$val = $cov;
			$xstart = $pos;
		}
		print TITI $chr_ucsc."\t".$pos."\t".$coverage->coverage->[$i]."\n";
	#	$array_intspan->set_range($pos,$pos,$coverage->coverage->[$i]);
			next unless $span->{$chr_ucsc}->contains($pos);
   		$nb++;
   		my $score = $coverage->coverage->[$i];
   		$tot += $score;
   		  	if ($score >= 100){
  		
   			$nb_100++;
   		}   
   		if ($score >= 30){
  		
   			$nb_30++;
   		}   
   		if ($score >= 50){
  		
   			$nb_50++;
   		}
   		if ($score >= 15){
   			$nb_15++;
   		}  
   		if ($score >= 20){
   			$nb_20++;
   		}  	 
   		if ($score >= 5){
   			$nb_5++;
   		}  	 
		
	}
			
	}  
	print TITI "mean_$chr_ucsc\t1\t$nb\n";
	print TITI "mean_$chr_ucsc\t5\t$nb_5\n";
	print TITI "mean_$chr_ucsc\t15\t$nb_15\n";
	print TITI "mean_$chr_ucsc\t20\t$nb_20\n";
	print TITI "mean_$chr_ucsc\t30\t$nb_30\n";
	print TITI "mean_$chr_ucsc\t50\t$nb_50\n";
	print TITI "mean_$chr_ucsc\t99\t$tot\n";
			print TITI "mean_$chr_ucsc\t100\t$nb_100\n";
	print "end $chr\n" unless $verbose; 
	close TITI; 
	
	#$array_intspan->consolidate(1,$max);
	my $z = time - $t;
	
		colored::stabilo('white',"  ---- end   $chr_num $patient_name  ---- ".abs(time-$t));

}

sub get_primers_array_intspan{
	my ($chr_name,$project_name) = @_;
	my $buffer3 = GBuffer->new();
 	my $project3 = $buffer3->newProject(-name=>$project_name);
	my $chr = $project3->getChromosome($chr_name);
	
	my $island;
	warn "p 0";
	foreach my $p (@{$chr->getPrimers}){
			my $set1 = Set::IntSpan::Island->new($p->start."-".$p->end);
		     $island->{$p->id."_count"} = $set1;
	}
	warn "p 1";
my $covers = Set::IntSpan::Island->extract_covers($island);
	my $foo = Array::IntSpan->new();
	foreach my $cover (@$covers) {
		next  if scalar(@{$cover->[1]}) == 0;	
		my $set = $cover->[0];	
		my $min = min $set;
		my $max = max $set;
		my @spans = spans $set;
		confess("problem ") unless scalar(@spans) ==1;
		$foo->set_range($spans[0][0],$spans[0][1],$cover->[1]);	
	}
 return $foo;	
}

sub by_chr_samtools {
my ($chr_num,$span,$span_extended,$name,$filein,$filetemp) = @_;
die();
my $bam_file = $filein;
	return unless $span_extended;
      my $is_ucsc; 
    $is_ucsc =1 ;#   if ($seq_ids[0] =~ /chr/);

	my $array_intspan =  Array::IntSpan->new();
	my $out = $filetemp->filename;
	my $chr = $chr_num;
	 $chr = "chr".$chr_num  if ischrornot($filein);
	$chr="chrM" if $chr eq "chrMT";
	my $chr_ucsc = "chr" . $chr_num;
	$chr_ucsc ="chrM" if $chr_ucsc eq "chrMT";
	#warn map{$_ = $chr_ucsc."\t".$_} split(";",$span_extended->{$chr_ucsc}->as_string({ sep => ";", range => "\t" }));
	my @bed;
	return unless $span_extended->{$chr_ucsc};
	 map{push(@bed,$chr_ucsc."\t".$_)} split(";",$span_extended->{$chr_ucsc}->as_string({ sep => ";", range => "\t" }));	

	my $nb=0;
	my $tot=0;
	my $nb_5=0;
	my $nb_15=0;
	my $nb_20 =0;
	my $nb_30=0;
	my $nb_100=0;
	my $window = 100;
	open (TITI,">$out",) or die("can't open $out");
	foreach my $b (@bed){
		my ($chr,$start,$end) = split (" ",$b);
		my $sstart =0;
		my $pos;
		for ( $sstart = $start; $sstart <$end -$window;$sstart += $window){
				push(@$pos,{start=>$sstart,end=>$sstart+$window -1});
		}
		push(@$pos,{start=>$sstart,end=>$end});
		foreach my $p (@$pos){
			my $tstart = $p->{start};
			my $tend = $p->{end};
			
		my ($gc) = $buffer->coverage_samtools($bam_file,$chr,$tstart,$tend);
		my $z =0;
		for (my $i=$tstart;$i<=$tend;$i++){
			my $depth = $gc->depth($i);
			print TITI $chr."\t".$i."\t".$depth."\n";
			next unless $span->{$chr_ucsc}->contains($i);
		
			$tot += $depth;
			$nb++;
				if ($depth >= 100){
   				$nb_100++;
   			}   
   			if ($depth >= 50){
   				$nb_50++;
   			}   
   			if ($depth >= 30){
   				$nb_30++;
   			}   
   			if ($depth >= 15){
   				$nb_15++;
   			}  	 
   			if ($depth >= 20){
   				$nb_20++;
   			}  
   			if ($depth >= 5){
   				$nb_5++;
   			}  	 	
		}
	}
	}


	print TITI "mean_$chr_ucsc\t1\t$nb\n";
	print TITI "mean_$chr_ucsc\t5\t$nb_5\n";
	print TITI "mean_$chr_ucsc\t15\t$nb_15\n";
		print TITI "mean_$chr_ucsc\t20\t$nb_20\n";
	print TITI "mean_$chr_ucsc\t30\t$nb_30\n";
	print TITI "mean_$chr_ucsc\t50\t$nb_50\n";	
	print TITI "mean_$chr_ucsc\t99\t$tot\n";
	print TITI "mean_$chr_ucsc\t100\t$nb_100\n";
	print "end $chr\n" unless $verbose; 
	close TITI; 
	colored::stabilo('white',"  ---- end   $chr_num $patient_name  ----");
}

sub by_chr_samtoolsfast {
my ($project,$chr_num,$span,$span_extended,$name,$filein,$filetemp) = @_;
$project->buffer->dbh_reconnect();
my $chr = $project->getChromosome($chr_num);

my $t = time;
my $bam_file = $filein;
	return unless $span_extended;
	
      my $is_ucsc; 
    $is_ucsc =1 ;#   if ($seq_ids[0] =~ /chr/);
	my $array_intspan =  Array::IntSpan->new();
	my $out = $filetemp->filename;

	my $chr_ucsc =  $chr->ucsc_name;
	
	$chr_ucsc ="chrM" if $chr_ucsc eq "chrMT";
	#warn map{$_ = $chr_ucsc."\t".$_} split(";",$span_extended->{$chr_ucsc}->as_string({ sep => ";", range => "\t" }));
	my @bed;
	return unless $span_extended->{$chr_ucsc};
	 map{push(@bed,$chr->fasta_name."\t".$_)} split(";",$span_extended->{$chr->ucsc_name}->as_string({ sep => ";", range => "\t" }));	
	my $nb=0;
	my $tot=0;
	my $nb_5=0;
	my $nb_15=0;
	my $nb_20=0;
	my $nb_30=0;
	my $nb_50=0;
	my $nb_100=0;
	my $window = 100;
	open (TITI,">$out",) or die("can't open $out");
	open (SAM,"$samtools depth -d 50000 -Q 1 $filein -r ".$chr->fasta_name ." | ");
	warn "$samtools depth -d 50000 -Q 1 $filein -r ".$chr->fasta_name;
	print TITI "$chr_ucsc\t0\t1\t0\n";
	my $last =1;
		my $intspan0 = Set::IntSpan::Fast->new();
	while (<SAM>){
		print TITI   $_;
		chomp();
		my ($chr,$pos,$depth) = split(" ");
		my $find;
		my $size = 1;

		
		
		 if ($span->{$chr_ucsc}->contains($pos)){
		
			$tot += $depth;
			
		#	$nb+=scalar($i1->as_array);
		
   			if ($depth >= 100){
   				$nb_100++;
   			}  
   			if ($depth >= 50){
   				$nb_50++;
   			}    
   			if ($depth >= 30){
   				$nb_30++;
   			}   
   			if ($depth >= 15){
   				$nb_15++;
   			}  	
   			if ($depth >= 20){
   				$nb_20++;
   			}  
   			if ($depth >= 5){
   				$nb_5++;
   			} 
		 }
		 $last = $pos;
		 
   		
   			 	 	
		
	}
	
	$nb = scalar($span->{$chr_ucsc}->as_array);
	
	my $h;
	$h->{nb} = $nb;
	$h->{5} =  $nb_5;
	$h->{15} =  $nb_15;
	$h->{20} =  $nb_20;
	$h->{30} =  $nb_30;
	$h->{50} =  $nb_50;
	$h->{100} =  $nb_100;
	$h->{99} =  $tot;
	my $fasta_name = $chr->fasta_name;
	print TITI "mean_$fasta_name\t1\t$nb\n";
	print TITI "mean_$fasta_name\t5\t$nb_5\n";
	print TITI "mean_$fasta_name\t15\t$nb_15\n";
	print TITI "mean_$fasta_name\t20\t$nb_20\n";
	print TITI "mean_$fasta_name\t30\t$nb_30\n";
	print TITI "mean_$fasta_name\t50\t$nb_50\n";
	print TITI "mean_$fasta_name\t99\t$tot\n";
		print TITI "mean_$fasta_name\t100\t$nb_100\n";
	print "end $chr\n" unless $verbose; 
	close TITI; 
	colored::stabilo('white',"  ---- end   $chr_num $patient_name  ---- ".abs(time-$t));
	return($h);
}


my %dejavu_file;
sub ischrornot {
	my ($file) = @_;
		return $dejavu_file{$file} if exists $dejavu_file{$file};
	my @res = `$samtools view -H $file`;

	my ($find) = grep { $_ =~ /SN:chr/ } @res;
	$dejavu_file{$file} = $find;
	return $find;
}


sub run {
	my ($cmd) = @_;
	my $return = system($cmd);
	if ($return ne 0){
		confess("error : $cmd");
	}
}

sub interrupt {
    print STDERR "Caught a control c!\n";
    map {delete $temp_file{$_}} keys %temp_file; 
    exit;  # or just about anything else you'd want to do
}

sub cache_coverage {
	my ($project_name,$fork,$patient_name) = @_;
	my $buffer1 = new GBuffer;
	my $projectP = $buffer1->newProject( -name => $project_name, -verbose =>1 );
	my $pm2 = new Parallel::ForkManager($fork);  
  	my $nbb =1;
  	warn "###################################################\n";
	warn "Prepare Cache for  COVERAGE !!!\n";
	warn "###################################################\n";

	

  my $diag;
  $diag=1 if $projectP->isDiagnostic();
  
# unless ($diag){
	
my $output   =$projectP->getCacheDir() . "/coverage_lite";
compute_coverage_diagnostic($project_name,$patient_name);
  			
	
#	preload_coverage::load_cnv_score($projectP,$project->getPatients,\@transcripts,1);
}

sub compute_coverage_diagnostic{
	my ($project_name,$patient_name) = @_;
		my $buffer = new GBuffer;
	my $project = $buffer->newProject( -name => $project_name, -verbose =>1 );
	my $patient = $project->getPatient($patient_name);
	my $no =  $project->noSqlCoverage();
	$no->clear($patient->name);
	$no->clear($patient->name."_cnv");
	
	my @paddings= (0,5,10,15,20,30);
	my @utrs = (0,1);
	my @trs = map{$project->newTranscript($_)} @{$project->bundle_transcripts()} ;
	warn Dumper @{$project->bundle_transcripts()};
	warn "1";
	preload_coverage::load_coverage_transcripts($project,[$patient],\@trs,1);
	warn "2";
	preload_coverage::load_coverage_primers($project,[$patient],\@trs,1) ;
	warn "3";
		foreach my $utr (@utrs){
					foreach my $padding (@paddings){
						warn $utr." ".$padding;
						preload_coverage::load_coverage($project,[$patient],\@trs ,$padding,$utr);
						
					
					}
		}
		$no->close();
		return;
}

