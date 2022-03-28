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
use polyweb_dude;
use Set::IntSpan::Fast::XS;
use Set::IntervalTree;

my $limit = 4;
my $fork = 1;
my $force;
my ($project_name, $patient_name);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'patient=s'    => \$patient_name,
	'force=s'  => \$force,
	
);


compute_patients($patient_name);


exit(0);

sub compute_patients {
my ($patient_name) =@_;
my $buffer = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProjectCache( -name => $project_name);
 my $lists = $project->getListTranscripts();

my $tree;

my $patient = $project->getPatient("$patient_name");
my $chr = $project->getChromosome(1);
unless ($project->isGenome){
	$limit = 2;
	#exit(0);
}

my $quality = "high";
 my %primers;
 #my $lists = $project->getListTranscripts();
 my $nbt =  scalar(@$lists);
 my $nb =0;
 my $hres;
 my $nb_tr;

 
foreach my $tra (@{$lists}) {
	warn $nb."/".$nbt if $nb%5000 == 0;
	#warn $tra;
	
	$nb++;
	#last if $nb > 200;
 	my $transcript = $project->newTranscript($tra);
 	my $gene = $transcript->getGene();
 	
 	my $debug ;
 	#$debug =1 if $gene->external_name eq "WASH7P";
 	#$debug =1 if $gene->external_name eq "MIR1302-2HG";
 	#$debug =1 if $gene->external_name eq "PLCXD1";
 	#next unless $debug;
 	#warn "coucou ";
 	my $coverage = polyweb_dude->new(patients=>[$patient],transcript=>$transcript,utr=>0,padding=>10,limit=>20,intronic=>0);

 	my $levels = $coverage->levels_matrix;
 	my $data_smoothed1 = $coverage->score_smooth_matrix;
 	my $data = $coverage->scores_matrix;
 	my $max =0;
 	my $error = 0;
 	my $span_del = Set::IntSpan::Fast::XS->new;
 	my $span_homo_del = Set::IntSpan::Fast::XS->new;
 	my $span_dup = Set::IntSpan::Fast::XS->new;
	my $nb_value = scalar(@$levels);
	my $homo_del;
 	foreach (my $i=0;$i<@$levels;$i++){
# 		warn $gene->external_name;
 #		warn Dumper $data_smoothed1 if $debug;
 #		warn Dumper $data if $debug;
 #		warn Dumper $coverage->score_smooth_expo_matrix if $debug; 
 		warn $data->[$i]->[0]." ".$data_smoothed1->[$i]->[0] if $debug;
 	#	warn $data->[$i]->[0].' '.$data_smoothed1->[$i]->[0];
 		if ($data->[$i]->[0] < 0.1 or $data_smoothed1->[$i]->[0] < 0.15){
 			$span_homo_del->add($i);
 		}
 		my $l =  $levels->[$i]->[0];
 		if ($l ==1) {
			$span_del->add($i);
		}
		
		if ($l == 2){
#				warn "coucou";
			$span_dup->add($i);
		}
		
		if ($l == 0){
			if ($data->[$i]->[0] < 0.65 or $data_smoothed1->[$i]->[0] < 0.7 ){
				$span_del->add($i);
			}
			if ($data->[$i]->[0] > 1.30 or  $data_smoothed1->[$i]->[0] > 1.40 ){
			
				$span_dup->add($i);
			}
		}
		#next if $l == 0 ;
		#next if( between($data->[$i],70,140) && $l <= 0);
	}
 #for deletion 
 warn $span_del->as_string() if $debug;
 my ($max_del,$nb_del) = Intspan_length($span_del);		
 my ($max_dup,$nb_dup) = Intspan_length($span_dup);
# warn "duplicate :".$max_dup." ".$nb_dup;
 my ($max_homo,$nb_homo) = Intspan_length($span_homo_del);
 my $chr = $transcript->getChromosome();
 warn $max_del." ".$nb_del." value: $nb_value ".$nb_dup if $debug;

 if($max_homo >0 ){
 	$nb_tr ++;
 	foreach my $i ($span_homo_del->as_array){
 		my $n = $coverage->names->[$i];
 		my $primer = $transcript->getPrimer($n);
 		next if $primer->cnv_score($patient) == -1;
 		my $p;
 		#warn  $patient->cnv_value_dude($primer->getChromosome->name,$primer->start,$primer->end) if $debug;
 		my $sd = $patient->sd_value_dude($primer->getChromosome->name,$primer->start,$primer->end) if $debug;
		#warn $sd if $debug;
		#warn $primer->getChromosome->name.":".$primer->start."-".$primer->end;
		#warn "++".$patient->getCapture->depth_controls_dude->getMean( $primer->getChromosome->name,$primer->start,$primer->end );
 		#warn "--".$patient->getNoSqlDepth->getMean( $primer->getChromosome->name,$primer->start,$primer->end );
 		
 		$p->{start} = $primer->start;
 		$p->{end} = $primer->end;
 		$p->{gt} = "1/1" if between($data->[$i]->[0],0,0.15);
 		$p->{type} = "del";# if ($data_smoothed1->[$i],0,0.15) ;
 		$p->{score} = $data->[$i]->[0];
 		$p->{score2} = $data_smoothed1->[$i]->[0]."\tA\tB\tC";
 		$p->{genes} = join(";",map{$_->external_name} @{$primer->getGenes});
 		
 		$p->{id} = $primer->id;
 		  $tree->{$chr->name} = Set::IntervalTree->new unless exists  $tree->{$chr->name};
 		$tree->{$chr->name}->insert($p, $p->{start}, $p->{end}+1);
 		push(@{$hres->{$chr->name}},$p);
 		#warn $n;
 		
 	}
 	#die() if $debug;
 }
 
 #warn $max_dup." ".$limit." ".$nb_dup." ".($nb_value *0.75)." ".$nb_del;
 #26 4 45 41.25 0
 #die();
 if (($max_del>=$limit or $nb_del > $nb_value *0.5 ) && ($nb_dup <2)) {
 	$nb_tr ++;
 	foreach my $i ($span_del->as_array){
 		my $n = $coverage->names->[$i];
 		my $primer = $transcript->getPrimer($n);
 		my $p;
 		$p->{start} = $primer->start;
 		$p->{end} = $primer->end;
 		$p->{gt} = "0/1";
 		$p->{gt} = "1/1" if between($data->[$i]->[0],0,0.15);
 		$p->{type} = "del";# if ($data_smoothed1->[$i],0,0.15) ;
 		$p->{score} = ($data->[$i]->[0]);
 		$p->{score2} = ($data_smoothed1->[$i]->[0]/100)."\tA\tB\tC";
 		$p->{genes} = join(";",map{$_->external_name} @{$primer->getGenes});
 		$p->{id} = $primer->id;
 		 $tree->{$chr->name} = Set::IntervalTree->new unless exists  $tree->{$chr->name};
 		$tree->{$chr->name}->insert($p, $p->{start}, $p->{end}+1);
 		push(@{$hres->{$chr->name}},$p);
 		#warn $n;
 		
 	}
 }
 
 elsif (($max_dup>=$limit or $nb_dup > $nb_value *0.75 ) && ($nb_del <4 && $nb_del < $nb_value *0.2 )) {
 	$nb_tr ++;
 	foreach my $i ($span_dup->as_array){
 		my $n = $coverage->names->[$i];
 		my $primer = $transcript->getPrimer($n);
 		my $p;
 		$p->{start} = $primer->start;
 		$p->{end} = $primer->end;
 		$p->{gt} = "0/1";
 		$p->{gt} = "1/1" if $data->[$i]->[0] > 1.9 ;
 		$p->{type} = "dup";# if ($data_smoothed1->[$i],0,0.15) ;
 		$p->{score} = $data->[$i]->[0];
 		$p->{genes} = join(";",map{$_->external_name} @{$primer->getGenes});
 		$p->{id} = $primer->id;
 		$p->{score2} = ($data_smoothed1->[$i]->[0]/100)."\tA\tB\tC";
 	 	$tree->{$chr->name} = Set::IntervalTree->new unless exists  $tree->{$chr->name};
 		$tree->{$chr->name}->insert($p, $p->{start}, $p->{end}+1);
 		push(@{$hres->{$chr->name}},$p);
 		#warn $n;
 		
 	}
 }
 		
#die();
}
#warn Dumper $hres;
my $dir_out= $project->getVariationsDir("dude");
my $bgzip = $buffer->software("bgzip");
my $tabix = $buffer->software("tabix");
my $file = $dir_out."/".$patient->name.".dude.lid";
my $filegz = $dir_out."/".$patient->name.".dude.lid.gz";
if (-e $filegz){
		unlink $filegz;
		unlink $filegz.".tbi";
}
#	warn $file;
my $nb_genes;
my @beds;
open (BED,">$file");

foreach my $chr (@{$project->getChromosomes}){
	my @vs = sort{$a->{start} <=> $b->{start}} @{$hres->{$chr->name}};
	my $set = Set::IntSpan::Fast::XS->new;
	
		foreach my $v (@vs) {
			$set->add_range(($v->{start}-2000),($v->{end}+2000));
		}
		my $iter = $set->iterate_runs();
	my @tt;
	my $l =0;
	my $max = -10;
    while (my ( $from, $to ) = $iter->()) {
    	my $tt = $tree->{$chr->name}->fetch($from,$to+1);
    	my $sc ;
    	my $sgt;
    	my $stype;
    	my $hgenes;
    	foreach my $t ( @$tt){
    		push(@$sc,$t->{score});
    		$sgt->{$t->{gt}} ++;
    		map{$hgenes->{$_}++} split(";",$t->{genes});
    		$stype->{$t->{type}} ++;
    	}
    	
    	my @type = sort {$stype->{$a} <=> $stype->{$b}} keys %$stype;
    	my @gt = sort {$sgt->{$a} <=> $sgt->{$b}} keys %$sgt;
    	my $mean = sum(@$sc)/scalar(@$sc);
    	my $genes = join(";",keys %$hgenes);
    	my @scores_sorted = sort { $a <=> $b } @$sc;
    	
    #	warn "sd:".$patient->sd_value_dude( $chr->name, $from, $to )." control mean ".$patient->getCapture->depth_controls_dude->getMean( $chr->name."_2", $from, $to ) . " pmean" . $patient->getNoSqlDepth->getMean( $chr->name, $from, $to )." ";
    #	warn $patient->cnv_value_dude( $chr->name, $from, $to );
    	print BED $chr->name."\t".$from."\t".$to."\t".$gt[-1]."\t".$type[-1]."\t".$genes."\t".$mean."\t".join(";", @scores_sorted)."\t0\n";
    	
    }
#	my $dj;
#	foreach my $v (@vs) {
#		next if exists $dj->{$v->{start}};
#		$dj->{$v->{start}} ++;
#		#my ($chr_name,$start,$end,$status,$type,$gene,$score1,$a,$score2,$b)
#		print BED $chr->name."\t".$v->{start}."\t".$v->{end}."\t".$v->{gt}."\t".$v->{type}."\t".$v->{genes}."\t".$v->{score2}."\n";
#		
#	}	
}
close (BED);
system("$bgzip $file && $tabix -p bed $file.gz ");
warn $filegz;
exit(0) if -e $filegz;
die();
}
#warn Dumper $hres;
sub Intspan_length{
	my ($intspan) = @_;
	my $iter = $intspan->iterate_runs();
	my @tt;
	my $l =0;
	my $max = -10;
    while (my ( $from, $to ) = $iter->()) {
    	my $d = ($to-$from) +1;
    		$l += $d;
    		$max = $d if $d > $max;
    	
    }
	#	my @tt = map{$_ = $chr->ucsc_name."\t".$_} split(";",$intspan->as_string({ sep => ";", range => "\t" }) );
		return ($max,$l);
}		

		

sub between {
  my($test ,$fom, $tom)=@_;
  no warnings;
  $fom<$tom ? $test>=$fom && $test<=$tom
            : $test>=$tom && $test<=$fom;
}

		