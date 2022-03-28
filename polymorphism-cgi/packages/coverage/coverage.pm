package coverage;
use Storable qw/thaw/;
use strict;
use Set::IntSpan::Fast::XS;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Statistics::Descriptive; 
use Data::Dumper;
use Carp;
use List::MoreUtils qw (first_index);
my $buffer;

sub get_capture_span {
	my ($project) = @_;
	my $capture_file = $project->getCaptureFile();
	
	open (CAP,"zcat $capture_file|");

	my $gene;
	my $trs;
	my $span;
	while (<CAP>){
		chomp();
		my($chr_name,$start,$end) = split(" ",$_);
		$chr_name =~ s/chr//;
		my $chr = $project->getChromosome($chr_name);
		unless (exists $span->{$chr_name}){ 
			$span->{$chr_name} = new Set::IntSpan::Fast::XS();
		}
		$span->{$chr_name}->add_range($start,$end);
		map{$trs->{$chr_name}->{$_}++} @{$chr->getFastTranscriptsByPosition($start,$end)};
	}
	close (CAP);
	
	return($span,$trs);
}


our $span_cover;
sub get_capture {
	my ($project) = @_;
	my $capture_file = $project->getCaptureFile();
	open (CAP,"zcat $capture_file|");

	my $gene;
	my $trs;
	my $span;
	while (<CAP>){
		chomp();
		my($chr_name,$start,$end) = split(" ",$_);
		$chr_name =~ s/chr//;
		unless (exists $span->{$chr_name}){ 
			$span->{$chr_name} = new Set::IntSpan::Fast::XS();
		}
		$span->{$chr_name}->add_range($start,$end);
	}
	close (CAP);
	
	return($span);
}

sub get_exons{
	my ($chr,$tr_name,$padding) = @_;
	 $span_cover = get_capture($chr->getProject()) unless $span_cover;
	
	my $tr = thaw  $chr->getKyotoHash("ensembl","transcript")->get($tr_name);
	my $span_genomic = $tr->{genomic_span};
	my $span_coding = $tr->{span_coding};
	$span_coding = new Set::IntSpan::Fast::XS() unless $span_coding;
	my $span_utr = $span_genomic->diff($span_coding);
	
	my $exons;
	my $iter     = $span_genomic->iterate_runs();
	
	
	$exons = [];
	my $pos= [];
	while ( my ( $from, $to ) = $iter->() ) {
		my $hpos;
		$hpos->{start} = $from;
		$hpos->{end} = $to;
		push(@$pos,$hpos);
	}
	
	my $num_exon = 1;
	my $start_cds =1;
	$num_exon = scalar(@$pos) if  $tr->{strand} == -1;
	
	foreach my $hp (@$pos){
		my $from = $hp->{start};
		my $to = $hp->{end};
		
		my $hpos;
		my $ps = new Set::IntSpan::Fast::XS($from."-".$to);
		$hpos->{chromosome} = $chr->name;
		
	
		$hpos->{gstart} = $from;
		$hpos->{gend} = $to;
		$hpos->{id}                     = $num_exon;
		$num_exon+= $tr->{strand};
		$hpos->{ext} ="ex";
		if ($span_utr->contains_all_range($from,$to)){
			#$hpos->{ext} ="utr";
			$hpos->{utr} =1;
			$hpos->{ext2} = "NC";
			#next;
		}
		my $ex_span = $span_coding->intersection($ps);
		my @tt = split("-",$ex_span->as_string());
		die() if scalar(@tt)>2;

		$hpos->{start} = $tt[0];
		$hpos->{end}   = $tt[1];
		$hpos->{name}= $hpos->{ext}.$hpos->{id}.$hpos->{ext2};
		$hpos->{length}   = ($to-$from)+1;
#		$hpos->{start_cds} = $start_cds;
#		$hpos->{end_cds} = $start_cds + $hpos->{length};
		$hpos->{strand}   = $tr->{strand};
		my $ps = new Set::IntSpan::Fast::XS($from."-".$to);
		
		my $st = $ps->diff($span_cover->{$chr->{name}});
		#my @points = $st->as_array();
		$hpos->{iscover} = scalar($st->as_array());
		$hpos->{intron_5} = new Set::IntSpan::Fast::XS("1-$padding");
		my $len = $hpos->{start}-$hpos->{end}+1;
		$hpos->{intron_3} = new Set::IntSpan::Fast::XS("$len-".($len+$padding));
		$start_cds += $hpos->{length}+1;
		push( @$exons, $hpos );
	}
	my @temp = sort {$a->{start} <=> $b->{start}} @$exons;
	return (\@temp,$tr);
}



sub get_all_coverage_data {
	my($tr_name,$exons,$patient,$padding) = @_;
	
	my $file = $patient->getCoverageFile();
	my $tabix = new Tabix(-data =>$file);
  	my $sall = Statistics::Descriptive::Full->new();
	foreach my $exon (@$exons){
		 $exon->{_stat} = Statistics::Descriptive::Full->new() unless exists $exon->{stat};
		
		my $data = data($tabix,$exon->{chromosome},$exon->{start}-$padding,$exon->{end}+$padding);;
	    my $s = Statistics::Descriptive::Full->new();
		$s->add_data(@$data);
		$patient->{_data}->{$tr_name}->{$exon->{name}} = $data;
		my $m = $s->mean();
		$patient->{_mean}->{$tr_name}->{$exon->{name}} = $m;
		$patient->{_min}->{$tr_name}->{$exon->{name}} = $s->min();
		$patient->{_max}->{$tr_name}->{$exon->{name}} = $s->max();
		$exon->{_stat}->add_data(@$data);
		for ( my $i=0;$i<@$data;$i++){
			$exon->{_data}->[$i] += $data->[$i];
		}
	
		$sall->add_data(@$data);
	}
	
	$patient->{_mean_all}->{$tr_name} = $sall->mean();
	$patient->{_min_all}->{$tr_name} = $sall->min();
	$patient->{_max_all}->{$tr_name} = $sall->max();
	foreach my $exon (@$exons){
		if ($patient->{_mean_all}->{$tr_name} == 0){
			next;
		}
		
		$patient->{_meanc}->{$tr_name}->{$exon->{name}} = $patient->{_mean}->{$tr_name}->{$exon->{name}}/$patient->{_mean_all}->{$tr_name};
	}
	
}
my %dejavu_statistics;

sub get_statistics {
	my ($tr_name,$exons,$patients,$padding,$limit) = @_;
	return if exists $dejavu_statistics{$tr_name};
	 $dejavu_statistics{$tr_name}++;
	foreach my $patient (@$patients){
		get_all_coverage_data($tr_name,$exons,$patient,$padding);
	}
	my $nb = scalar(@$patients);
	my $sall = Statistics::Descriptive::Full->new();
	
	foreach my $exon (@$exons){
		$sall->add_data($exon->{_stat}->get_data());
		$exon->{_mean} = $exon->{_stat}->mean();
		$exon->{_min} = int( $exon->{_stat}->min());
		$exon->{_max} = int( $exon->{_stat}->max());
	}
	my $mean = $sall->mean();
	foreach my $exon (@$exons){
		next if $mean ==0;
		$exon->{_meanc} = $exon->{_mean}/$mean;
		$exon->{_meanall} = $mean;
	}
	
}

sub data {
	my ($tabix,$chr,$start,$end) = @_;
	my $res = $tabix->query("chr".$chr,$start,$end);
	my $previous = $start;
	my $data;
    while(my $line = $tabix->read($res)){
    		   
				my($a,$b,$c) = split(" ",$line);
				
				if (abs($b-$previous)>1){
				 	for (my $i=$previous+1;$i<$b;$i++){
					 	
				 		 push(@$data,0);
				 	}
				}
				$previous = $b;
				 push(@$data,$c);
			
				
			
    		}
    		foreach (my $i=$previous+1;$i<$end;$i++){
    			 push(@$data,0);
    		}
    return $data;		
}

sub coverage_data {
	my ($exon,$patient) = @_;
	my $file = $patient->getCoverageFile();
		my $tabix = new Tabix(-data =>$file);
		my $data;
    	my $res = $tabix->query("chr".$exon->{chromosome},$exon->{start},$exon->{end});
    	my $previous = $exon->{start};
    while(my $line = $tabix->read($res)){
    		   
				my($a,$b,$c) = split(" ",$line);
				
				if (abs($b-$previous)>1){
				 	for (my $i=$previous+1;$i<$b;$i++){
					 	
				 		 push(@$data,0);
				 	}
				}
				$previous = $b;
				 push(@$data,$c);
			
				
			
    		}
    		foreach (my $i=$previous+1;$i<$exon->{end};$i++){
    			 push(@$data,0);
    		}
    return $data;		
}
my $span_all;
my %dejavu_span;
sub get_all_span {
	my ($tr_name,$exons,$patients,$padding,$limit,$capture) = @_;
	return if exists $dejavu_span{$tr_name};
	
	 $dejavu_span{$tr_name} ++;
foreach my $patient (@$patients){
	
	foreach my $exon (@$exons){
		my $item_span;
		my $id = $exon->{name};
		next if exists $patient->{_span}->{$tr_name}->{$id};
		my $chr = $exon->{chromosome};
		($patient->{_span}->{$tr_name}->{$id}) = get_span($patient->{_data}->{$tr_name}->{$id},$limit,$padding);
		
	}#end exon
	}
	foreach my $exon (@$exons){
		my $id = $exon->{name};
			
		$exon->{_span} = get_span($exon->{_data},$limit,$exon->{start}-$padding);
			
	foreach my $patient (@$patients){
		
			$patient->{_coveragespan}->{$tr_name}->{$id} = $patient->{_span}->{$tr_name}->{$id}->diff( $exon->{_span});
			
		}
	}	
	
}

sub get_span {
	my ($data,$limit,$start) = @_;
	my $span = new Set::IntSpan::Fast::XS();
	
	my $nb =0;

	for (my $i =0;$i<@$data;$i++){
		
		$span->add($i+1) if $data->[$i] <=  $limit;

	
		
	}
	return ($span);
}

sub as_string {
	my ($span,$exon) = @_;
	my $iter     = $span->iterate_runs();
	my @string;
	my $debug;
	my $reverse;
	$reverse =1 if $exon->{strand} == -1;
	#$debug=1 if $exon->{name} eq "utr_1";
	while ( my ( $from, $to ) = $iter->() ) {
		my $pos_start;
		my $pos_end;
		my $len = $exon->{length}+20;
		warn $from." ".$to." ".$len if $debug;
		if ($reverse){
			my $temp = $to;
			$to = $len - $from;
			$from = $len - $temp;
		}
		warn $from." ".$to if $debug;
		
		my @values =($from,$to);
#		#@values= reverse (@values) if ($exon->{strand} == -1);
		for (my $i=0;$i<@values;$i++){
			
		if ($exon->{intron_5}->contains($values[$i] ) ) {
			my @array = reverse($exon->{intron_5}->as_array);
			my $index = first_index{$_ == $from } @array	;	
			$values[$i] = ($index+1);
			$values[$i] = "-";
			#$values[$i] = "+" if $reverse;	
			$values[$i] .= ($index+1);
		}
		if ($exon->{intron_3}->contains($values[$i])) {
			my $index = first_index{$_ == $values[$i] } $exon->{intron_3}->as_array;	
			$values[$i] = "+";
			#$values[$i] = "-" if $reverse;	
			$values[$i] .= ($index+1);	
		}
		}
		
#		if ($reverse){
#			push(@string,"[".$values[1]." -> ".$values[0]."]") ;
#		}
#		else {
			push(@string,"[".$values[0]." -> ".$values[1]."]") ;
#		}
	}
	warn "--".scalar($span->as_string) if $debug;
	die(join(";",@string)) if $debug;
	return(join(";",@string));
}
1;