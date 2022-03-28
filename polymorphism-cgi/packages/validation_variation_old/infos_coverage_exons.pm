package infos_coverage_exons;
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;


sub return_hash_exons_from_hash {
	my ($exon,$patient,$exons_todo,$capture_intspan,$tr1,$show_utr,$limit,$splice_5,$ph) = @_;
	
	
}

sub return_hash_exons{
	my ($exon,$patient,$exons_todo,$capture_intspan,$tr1,$show_utr,$limit,$splice_5,$ph) = @_;
	my $transcript = $tr1;
	
	my $hexons;
	my 	($mean,$intspan,$min) ;

	if ($show_utr == 1 || !($exon->isExon)){
			($mean,$intspan,$min) = $exon->mean_intspan_coverage(patient=>$patient,padding=>$splice_5,limit=>$limit);
		}
		else {
			($mean,$intspan,$min) =$exon->mean_intspan_coverage_coding(patient=>$patient,padding=>$splice_5,limit=>$limit);
		}
			
		my $s1 = $exon->getGenomicSpan()->intersection($capture_intspan);
		my $hexons;
		my $text = $exon->name();
		my $chr = $exon->getChromosome->name();
		my $start = $exon->start;
		my $end = $exon->end;
		$hexons->{name} = qq{<a href="javascript:displayInIGV('$chr',$start,$end);">$text</a>};
		$hexons->{vname} = $exon->name();
		$hexons->{mean} = int($mean);
		#$hexons->{intspan} = $intspan->as_string();
		$hexons->{min} = $min;
		$hexons->{todo} = 0; 
		$hexons->{id} = $exon->id();
		$hexons->{label}->{$patient->name()} = $patient->name().":".$tr1->getGene->external_name().":".$tr1->name.":".$exon->name.":".$tr1->getChromosome()->name.":".$exon->start.":".$exon->end;
		
		#my @nb = $intspan->as_array();
	
		$hexons->{color} = "green";
#			 
#		my $color = "tdgreen";
#		$hexons->{color} = "grey"  if $s1->is_empty();
#			if ($min <= $limit && !($s1->is_empty) ){
#				$color = "tdred";
#				$hexons->{color} = "red";
#			}
#			$hexons->{color} = "grey"  if $exon->is_noncoding()   && $show_utr ne 1;
#			
#			$hexons->{color} = "green";	
#			
#		
		my $covered = 1;
		if ($exon->isExon()){
			$covered =undef if $s1->is_empty;
		}
		else {
			$covered =undef if scalar($s1->as_array)<20;
		}
		$covered=undef  if $exon->is_noncoding()   && $show_utr ne 1 && $exon->isExon() ;

		if (exists $exons_todo->{$exon->coverage_id}) {
			
				$hexons->{type_string} = "Sanger ";
				$hexons->{color} = "blue";
				
				if ($exons_todo->{$exon->coverage_id}->{done} eq 1) {
					
				#	$htranscript->{table} = 0;
					$hexons->{type} = 0;
					$hexons->{type_string} = " OK " ;
				} 
				elsif ($exons_todo->{$exon->coverage_id}->{done} eq -1) {
					$hexons->{type} = -1;
					$hexons->{type_string} = "He Deletion" ;
				} 
				elsif ($exons_todo->{$exon->coverage_id}->{done} eq -2) {
					$hexons->{type} = -2;
					$hexons->{type_string} = "Ho Deletion" ;
				} 
			
				else {
					$hexons->{type} = 1;
					$hexons->{type_string} = "to do" ;
				}
		}
		elsif (!$covered){
				$hexons->{color} = "grey";
				$hexons->{type} = -3;
				$hexons->{type_string} = "NS" ;
		}		
		elsif ($min < $limit){
					$hexons->{color} = "red";
					$hexons->{type} = 2;
					$hexons->{type_string} = "to do" ;
				}
		
		else {
			$hexons->{type} = 99;
			$hexons->{type_string} .= " - " ;
		}
			
			
			return $hexons;
	
}



sub return_hash_exons2{
	my ($exon,$patient,$tr1,$exons_todo,$pdata,$limit) = @_;
	my $hexons;
	
		my $hexons;
		my $text = $exon->name();
		my $chr = $exon->getChromosome->name();
		my $start = $exon->start;
		my $end = $exon->end;
		$hexons->{name} = qq{<a href="javascript:displayInIGV('$chr',$start,$end);">$text</a>};
		$hexons->{vname} = $exon->name();
		$hexons->{mean} = int($pdata->{mean});
		if (exists $pdata->{raw_mean}){
		$hexons->{raw_mean} = int($pdata->{raw_mean}) ;
		$hexons->{raw_min} = int($pdata->{raw_min}) ;
		}
		#$hexons->{intspan} = $intspan->as_string();
		$hexons->{min} = $pdata->{min};
		$hexons->{todo} = 0; 
		$hexons->{id} = $exon->id();
		$hexons->{label}->{$patient->name()} = $patient->name().":".$tr1->getGene->external_name().":".$tr1->name.":".$exon->name.":".$tr1->getChromosome()->name.":".$exon->start.":".$exon->end;
		
		#my @nb = $intspan->as_array();
	
		$hexons->{color} = $pdata->{color};
		$hexons->{red} = $pdata->{red};
		my $covered = $pdata->{covered};
		$hexons->{color2} = $pdata->{color2};
		if (exists $exons_todo->{$exon->coverage_id}) {
				$hexons->{red} = 0;
				$hexons->{green} =0;
				$hexons->{type_string} = "Sanger ";
				
				if ($exons_todo->{$exon->coverage_id}->{done} eq 1) {
					$hexons->{type} = 0;
					$hexons->{type_string} = " OK " ;
					$hexons->{color} = [150,50,150];
					$hexons->{color2} = [0,0,255];
				} 
				elsif ($exons_todo->{$exon->coverage_id}->{done} eq -1) {
					$hexons->{type} = -1;
					$hexons->{type_string} = "He Deletion" ;
					$hexons->{color} = [0,250,250];
					$hexons->{color2} = [112,29,0];
				} 
				elsif ($exons_todo->{$exon->coverage_id}->{done} eq -2) {
					$hexons->{type} = -2;
					$hexons->{color} = [0,250,250];
					$hexons->{color2} = [255,0,0];
					$hexons->{type_string} = "Ho Deletion" ;
				} 
			
				else {
					$hexons->{type} = 1;
					$hexons->{color} = [50,50,200];
					$hexons->{color2} = [79,255,255];
					$hexons->{type_string} = "to do" ;
					
				}
		}
		
		
			
			return $hexons;
	
}
1;

