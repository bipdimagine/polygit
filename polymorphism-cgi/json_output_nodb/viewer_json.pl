#!/usr/bin/perl
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";

 
use GBuffer;
use GenBoStorable;
use Data::Dumper;
use Tabix;
use export_data;
use Sys::Hostname;
use Bio::DB::Sam; 
use get_variations;
use Set::IntSpan::Fast::XS;
use Storable qw/freeze thaw nfreeze nstore_fd nstore retrieve/;
#use Tabix;

my %types = ( variations => "variations" );

my $cgi    = new CGI();
my $buffer = new GBuffer;

my @cover_values = (0,5,10,15,20);

my $host = hostname;


my $project_name = $cgi->param('project');

my $type_name = "variations";
my $select_ids   = [split(",",$cgi->param('ids'))];
#creation du projet  partir du buffer
my $project = $buffer->newProject( -name => $project_name );
die( "unknown project" . $project_name ) unless $project;
my $chr;
my $no_seq = $cgi->param("no_seq");
if ($cgi->param('chromosome')){
	$chr = $project->getChromosome( $cgi->param('chromosome'));
}
else {
	$chr = $project->getChromosomes()->[0];
}
my $typeToDraw = $cgi->param('type');

#$project->getMethods($buffer->getType("variation"));


my $start =1_000_000_000;
my $end = -1;

$start = $cgi->param('start');
 $end = $cgi->param('end');
my $real_data;
###############################""
#
#  bygene
#
#################################""
my $gene_name = $cgi->param('gene') || $cgi->param('selectGene');

if ($gene_name){
	my $chr_name;
	($chr_name,$start,$end) = define_start_end($gene_name);

	 $chr = $project->getChromosome($chr_name);
	 #my $id = GenBoStorable::getStoreId( $buffer->dbh, $project->id, $chr->id,$type_name );
	 
	#$start = 288052;
	#$end = 304072;
}
elsif(scalar(@$select_ids)){
	
$real_data = get_variations::getIds($buffer,$project,$chr,$type_name,$select_ids);	
$start = $real_data->[0]->{start};
$end = $real_data->[-1]->{end};

my %gg;

foreach my $v (@$real_data){
	foreach my $g (@{$v->{genes}}){
		#warn $g;
		$gg{$g}++;
	}
}
my @ogenes = sort {$gg{$a} <=> $gg{$b}} keys %gg; 



if (scalar(@ogenes) > 0) {
	my ($chr_name,$start_genes,$end_genes) = define_start_end($ogenes[-1]);
	$start = $start_genes if $start>$start_genes;
	$end = $end_genes if $end<$end_genes;
}

#die();
}

#if  ($cgi->param('start')) {
#	$start = $cgi->param('start') if $cgi->param('start') < $start;
#	
#	$end = $cgi->param('end') if $cgi->param('end') > $end;	
#	#$real_data = get_variations::getIdsByPosition($buffer,$project,$chr,$start,$end,$type_name,$select_ids);
#	
#}






$start -= 500;
$end += 500;

my $data_out;
#($start,$end) = getStartEnd($project->buffer,$real_data->[0]->{chromosome},$start,$end,$data_out);
my $filter_patient = $cgi->param('patients_and');
my @and_patients;
if ($filter_patient) {
		  @and_patients = map {$project->getPatient($_)} split(" ",$filter_patient);
}
if ($typeToDraw ne 'transcripts') {
	($data_out->{patient},$data_out->{contig},$data_out->{contig2}) = constructVariations($project,$chr,$real_data,$start,$end);	
}
#die();
if ($typeToDraw ne 'patients') {
	$data_out->{transcript} = constructTranscripts($project->buffer,$chr->name,$start,$end,$data_out);
}
$data_out->{length} = $end-$start -1;
$data_out->{chromosome} = $real_data->[0]->{chromosome};
$data_out->{start} = 1;
$data_out->{end} = $data_out->{length};
$data_out->{gstart} = $start;
export_data::print($project,$cgi,$data_out);
exit(0);

sub define_start_end {
	my ($gene_name) = @_;
	my $chrs = [1..22,'X','Y','MT'];
	my $gene = $project->newGene($gene_name);
	
	return ($gene->getChromosome()->name,$gene->start-100,$gene->end+100);
	
	
	#return($hgene->{chr},$gene->start()-100,$gene->end()+100);
	
	
	
	
}
sub constructTranscripts {
	my ($buffer,$chr,$start,$end) = @_;
	my @out;
	my $limit_end = abs($end-$start)+50;
	#warn $start." ".$end;
	my $real_start = $start;
	#warn $chr;
	my $ochr = $project->getChromosome($chr);
	my $genes= $ochr->getGenesByPosition($start,$end);
	
	foreach my $gene ( @$genes){
		my $too_much;
		
		foreach my $tr ( @{$gene->getTranscripts}) {
			#warn $tr->name;
			#next if ($too_much && !$tr->translation); 
			my %transcript;
			$transcript{name} = $tr->name;
			$transcript{name} =  $tr->external_name if  $tr->external_name;
			$transcript{gene_name} = $gene->external_name;
			$transcript{x} = $tr->start-$real_start;
			
			next if $tr->start < $start;
			next if $tr->end > $end-50;
		
			$transcript{length} = ($tr->end-$tr->start +1) ;#$tr->length();
			my $pos;
			my $nb= 0;
			foreach my $exon (@{$tr->getExons}){
				my %hexon;
				$hexon{name} = $exon->name;
				$hexon{x} = $exon->start-$real_start;
				$hexon{length} = $exon->end-$exon->start +1;
				$hexon{strand} = $tr->strand;
				$hexon{name} = $tr->{id}."_".$nb;
				my $utr_string = $exon->utr->as_string();
				if ($utr_string){
				 ($hexon{utrstart},$hexon{utrend}) = split($utr_string); 
				 $hexon{utrstart} -= $real_start;
				 $hexon{utrend} -= $real_start;
				}
			#on ajoute les infos sur les exons  la table %transcrit du transcrit possdant ces exons
				push(@{$transcript{exon}},\%hexon);
			}
			push(@out,\%transcript);
		}#end transcript
	
	}#end gene
	
	return \@out;
}


sub constructVariations {
	my ($project,$chr,$data_variation,$start,$end) = @_;
	my %patients;
	my $traces=[];

	
	foreach my $var (@$data_variation) {
		my $var2;
		$var2->{start} = $var->{start} - ($start + 1);
		$var2->{x} = $var2->{start};
		 $var2->{length} = 1;#$var->{start};
		 $var2->{id} = $var->{id};
		  $var2->{name} = $var->{name};
		 $var2->{consequence} = $var->{"consequence!all"};
		 $var2->{"consequence!all"} = $var->{"consequence!all"};
		my @patients_names = split(";",$var->{allpatients});
		foreach my $p_name (sort(@patients_names)) {
			push(@{$patients{$p_name}},$var2); 
		}		
	}
	
	my @out;
	my $contig_span  = Set::IntSpan::Fast::XS->new();
	my $hOut;
	foreach my $pat (@{$project->getPatients()}) {
		my $p_name = $pat->name();
		my %p;
		$p{name} = $p_name;
		$p{variation} = $patients{$p_name};
		my @gtraces = grep {$_->getPatient()->name eq $p_name} @$traces;
		my @ptraces;
		($p{trace},$contig_span) = constructCover2($project,$chr,$start,$end,$pat,$contig_span);
		$hOut->{$p{name}} = \%p;
		#push(@out,\%p);
	}
	
	if (@and_patients){
		my %p;
		$p{name} = join("&",map{$_->name} @and_patients);
		$p{trace} = construct_cover_and(\@and_patients,$start,$end);
		$p{type} = "and";
		$hOut->{$p{name}} = \%p;
		#push(@out,\%p);
	}
	foreach my $patName (sort(keys(%$hOut))) { push(@out, $hOut->{$patName}); }
	
	#die();
	my @contigs;
	if ($contig_span){
		foreach my $pos (split(",",$contig_span->as_string)){	
			my %contig;
			my ($tstart,$tend) = split("-",$pos);		
			$contig{name} = "XX".$tstart;
			$contig{x} = $tstart;
			$contig{length} = ($tend-$tstart)+1;
			push(@contigs,\%contig);
		
		}
	}
	
		my $contig_agilent = constructSpanAgilent($project,$chr,$start,$end);
		return (\@out,$contig_agilent);
	
	
}

sub constructCover {
	my ($project,$chr,$start,$end,$pat,$contig_span) = @_;
	
	my $buffer = $project->buffer();
	my $bam_file = $pat->getBamFile();
	
	return([],$contig_span) unless -e $bam_file;
	my $chr_ucsc = $chr->ucsc_name();
	
	#my $command = qq{samtools  depth -r "$chr_ucsc:$start-$end" $bam_file | cut -f 3};
	my $command1 = qq{/bip-d/soft/bin/samtools view -H $bam_file | grep "SN:chr" };

	my ($tt) = `$command1`;

	unless ($tt){
		$chr_ucsc = $chr->name();
		
	}
	my $samtools = $project->buffer->{config}->{software}->{samtools};
	my $command = qq{$samtools  depth -r "$chr_ucsc:$start-$end" $bam_file | cut -f 2,3 | awk '\$2<10'};
	
	my @score;
	my %pos;
	
	 (%pos) = split(" ",`$command`);
	
	my $span_total = Set::IntSpan::Fast::XS->new("$start-$end");
	my $span10 = $span_total->diff(Set::IntSpan::Fast::XS->new(keys %pos));
	map {delete $pos{$_} if $pos{$_} <= 5  } keys %pos;
	my $span5 = $span_total->diff(Set::IntSpan::Fast::XS->new(keys %pos));
	map {delete $pos{$_} if $pos{$_} <= 1  } keys %pos;
	my $span = $span_total->diff(Set::IntSpan::Fast::XS->new(keys %pos));
    $contig_span = $contig_span->union($span5);
   my $iter2 = $span5->iterate_runs();
   my @traces;
	while (my ( $from, $to ) = $iter2->()) {
    		my $tr;
    		$tr->{x} = $from ;
    		$tr->{y} = $to;	
    				
    	$tr->{strand} = 1;	
    	push (@traces,$tr);
    			}
    			
      my $iter3 = $span10->iterate_runs();
	while (my ( $from, $to ) = $iter3->()) {
    		my $tr;
    		$tr->{x} = $from ;
    		$tr->{y} = $to;	
    				
    	$tr->{strand} = -1;	
    	push (@traces,$tr);
    			}			
   	 my $iter4 = $span->iterate_runs();
	while (my ( $from, $to ) = $iter4->()) {
    		my $tr;
    		$tr->{x} = $from ;
    		$tr->{y} = $to;	
    				
    	$tr->{strand} = 0;	
    	push (@traces,$tr);
    			}	
	
    return (\@traces,$contig_span);
}
my %patients_cover;
sub construct_cover_and {
	my ($patients,$start,$end) = @_;
	my $span_inter;
	foreach my $v (@cover_values){
		
		foreach my $p (@$patients){
			$span_inter->{$v} = $patients_cover{$p->name}->{$v}->copy unless exists $span_inter->{$v};
			$span_inter->{$v} = $span_inter->{$v}->intersection($patients_cover{$p->name}->{$v});
		
			
		}
	}
	 my @traces;
  
  foreach my $v (@cover_values){
  	 my $iter2 = $span_inter->{$v}->iterate_runs();
  	
	while (my ( $from, $to ) = $iter2->()) {
    		my $tr;
    		$tr->{x} = $from ;
    		$tr->{y} = $to;	
    				
    	$tr->{strand} = "cov".$v;	
    	push (@traces,$tr);
    		}
  }


return \@traces;
  
	
}

sub constructCover2 {
	my ($project,$chr,$start,$end,$pat,$contig_span) = @_;
	
	my $buffer = $project->buffer();
	my $coverage_dir = $project->getRootDir()."/align/coverage/";
	my $coverage_file = $coverage_dir.$pat->name().".cov.gz";
	my $tabix = $project->buffer->{config}->{software}->{tabix};
	return([],$contig_span) unless -e $coverage_file;
	my $chr_name = $chr->ucsc_name;
	my $tabix = new Tabix(-data =>$coverage_file);
	my $res = $tabix->query("$chr_name",$start,$end);
	my @data;
	eval {
	 while(my $line = $tabix->read($res)){
	 	push(@data,$line);
	 }
	};
#	my @data = `$tabix  $coverage_file $chr_name:$start-$end`;
#	chomp(@data);
	my $span;
	# = Set::IntSpan::Fast::XS->new();
  # my $span5 = Set::IntSpan::Fast::XS->new();
  # my $span10 = Set::IntSpan::Fast::XS->new();
  foreach my $v (@cover_values){
  	$span->{$v} = Set::IntSpan::Fast::XS->new();
  }
  my $size = 1;
  $size = int(scalar(@data) / 7000);
  $size = 1 if $size == 0;
  for (my $i = 0 ; $i < scalar(@data) - $size ; $i+=$size){
  	my $d = $data[$i];
  # foreach my $d (@data) {
   
   	my ($a,$b,$c) = split(" ",$d);
   	my ($a1,$b1,$c1) = split(" ",$data[$i+$size]);
  	$b-=$start;
  	$c = ($c+$c1) /2;
  	foreach my $v (@cover_values){
  	
  		if ($v == 0) {
  			$span->{0}->add_range($b,$b+$size) if $c == $v;
  		}
  		else {
  			$span->{$v}->add_range($b,$b+$size) if $c >= $v;
  		}
  		last if $c < $v;
  	}
   }
  
 $patients_cover{$pat->name} = $span;
    # $contig_span = $contig_span->union($span->{5});
     #  die();
  my @traces;
  
  foreach my $v (@cover_values){
  	 my $iter2 = $span->{$v}->iterate_runs();
  	
	while (my ( $from, $to ) = $iter2->()) {
    		my $tr;
    		$tr->{x} = $from ;
    		$tr->{y} = $to;	
    				
    	$tr->{strand} = "cov".$v;	
    	push (@traces,$tr);
    		}
  }
  
  $contig_span = $contig_span->union($span->{5});
    return (\@traces,$contig_span);
	
	
}

sub constructSpanAgilent {
	my ($project,$chr,$start,$end) = @_;
	my $span = $chr->getCapturesGenomicSpan();#getCaptureIntspan($chr,$start,$end);
	
	my $intSpan = Set::IntSpan::Fast::XS->new("$start-$end");
	my $final = $intSpan->intersection($span);
	my @contigs;
	foreach my $pos (split(",",$final->as_string)){	
			my %contig;
			my ($tstart,$tend) = split("-",$pos);		
			$contig{name} = "XX".$tstart;
			$contig{x} = $tstart-$start;
			$contig{length} = ($tend-$tstart)+1;
			push(@contigs,\%contig);
		
		}

	return \@contigs;
	
}
