package polyweb_coverage;
use strict;
use Moose;
use MooseX::Method::Signatures;

use GD;
use Data::Printer;
use Data::Dumper;
 use List::Util qw( min sum);
 use Statistics::Descriptive;
 use List::MoreUtils qw(indexes);
 use URI;
  use List::Util qw( min sum );
 
has transcript => (
	is		=> 'ro',
	required=> 1,
);

 has patients => (
	is		=> 'ro',
	required=> 1,
);

  has utr => (
	is		=> 'ro',
	required=> 1,
);
 
 has limit => (
	is		=> 'ro',
	required=> 1,
);

  has padding => (
	is		=> 'ro',
	required=> 1,
);

  has intronic => (
	is		=> 'ro',
	required=> 1,
);


  has min_matrix => (
	is		=> 'rw',
	lazy =>1,
	default => sub {
		my $self = shift;
	   $self->init_matrices();
	   return $self->{min_matrix};
}
);

  has mean_matrix => (
	is		=> 'rw',
	lazy =>1,
	default => sub {
		my $self = shift;
	   $self->init_matrices();
	   return $self->{mean_matrix};
}
);

  has type_matrix => (
	is		=> 'rw',
	lazy =>1,
	default => sub {
		my $self = shift;
	   $self->init_matrices();
	   return $self->{type_matrix};
}
);

 has error => (
	is		=> 'rw',
	lazy =>1,
	default => sub {
		my $self = shift;
	   $self->init_matrices();
	   return $self->{error};
}
);

 has cgi => (
	is		=> 'rw',
	lazy =>1,
	default => sub {
		return new CGI;
}
);

has grey_hex =>   (
	is		=> 'rw',
	lazy =>1,
	default => sub {
		return  ["#cccccc","#98B4D4","#FF0000"] 
}
);



has grey_rgb =>   (
	is		=> 'rw',
	lazy =>1,
	default => sub {
		return  [[200,200,200],[152, 180, 212],[127, 255, 0]];
}
);

has red_hex =>   (
	is		=> 'rw',
	lazy =>1,
	default => sub {
		return  ["#ff9900","#ff6666","#FF0000"] 
}
);

has red_rgb =>   (
	is		=> 'rw',
	lazy =>1,
	default => sub {
		return  [[255,165,0],[255,127,80],[255,0,0]];
}
);


has green_hex =>   (
	is		=> 'rw',
	lazy =>1,
	default => sub {
		return  ["#99cc33","#99cc00","#27AE60"] ;
}
);

has green_rgb  => (
	is		=> 'rw',
	lazy =>1,
	default => sub {
		return  [[175,206,30],[152,228,16],[46, 204, 113]];
}
);

has purple_rgb  => (
	is		=> 'rw',
	lazy =>1,
	default => sub {
		return  [[107, 91, 149],[181, 1, 254]];
}
);

has purple_hex  => (
	is		=> 'rw',
	lazy =>1,
	default => sub {
		return  ["#666699","#cc00ff"];
}
);

has mask_exons  => (
	is		=> 'rw',
	lazy =>1,
	default => sub {
	return {
	  exon => 1,
      intron => 2,
 	  non_coding => 4,
 	 capture =>8,
 	 duplicate => 16,
 	  "5_utr" => 32,
 	  "3_utr" => 64,
	};
}
);

sub between {
  my($test ,$fom, $tom)=@_;
  no warnings;
  $fom<$tom ? $test>=$fom && $test<=$tom
            : $test>=$tom && $test<=$fom;
}

has names  => (
	is		=> 'rw',
	lazy =>1,
	default => sub {
				my $self = shift;
				$self->types;
				return $self->{names};
	}
);

has types  => (
	is		=> 'rw',
	lazy =>1,
	default => sub {
			my $self = shift;
			my $patient = $self->patients->[0];
			my $no = $patient->getTranscriptsCoverageDepth("r");
			my $matrix = $no->get($self->transcript->id);
			$self->{names} = $matrix->{name};
			return [unpack("w".$matrix->{nb},$matrix->{types})];
	}
);

has types_index  => (
	is		=> 'rw',
	lazy =>1,
	default => sub {
			my $self = shift;
			$self->exon_indexes;
			return $self->{types_index};
	}
);

has exons_index => (
	is		=> 'rw',
	lazy =>1,
	default => sub {
			my $self = shift;
			my $mask = $self->mask_exons->{exon} ;

			my @index = indexes {$_ & $mask }  @{$self->types};
			foreach my $ind (@index) {
				push(@{$self->{types_index}},$self->types->[$ind]);
			}
			return \@index;
	}
);

has names_index => (
	is		=> 'rw',
	lazy =>1,
	default => sub {
			my $self = shift;
			my @name;
			foreach my $ind (@{$self->exons_index}) {
				push(@name,$self->names->[$ind]);
			}
			return \@name;
	}
);

sub init_matrices {
	my ($self) = @_;
	
	my $patients = $self->patients;
	my $transcript = $self->transcript;
	my $utr = $self->utr;
	my $limit = $self->limit;
	my $padding = $self->padding;
		
my $data;
my $nbp =0;
my $error;
my $data_mean;
my $data_names;
my $type_matrix;

my $names ;
my $first  = 1;
my $patients_mean;
	foreach my $patient (sort{$a->name cmp $b->name} @$patients){
	my $no = $patient->getTranscriptsCoverageDepth("r");
	my $matrix = $no->get($transcript->id);
	#warn $padding." ".$matrix;
	#warn $transcript->id." ".$patient->name;
	#warn Dumper $matrix;
	#warn Dumper $matrix->{min_no_utr};
	#warn Dumper $matrix->{mean_no_utr};
	my $m; 
	  my $mean;
	if ($self->utr == 0){
	 $m  = [ unpack("w".$matrix->{nb},$matrix->{min_no_utr}->{$padding})];
	 $mean  = [ unpack("w".$matrix->{nb},$matrix->{mean_no_utr}->{$padding})];
	}
	else {
		$m  = [ unpack("w".$matrix->{nb},$matrix->{min_utr}->{$padding})];
		$mean  = [ unpack("w".$matrix->{nb},$matrix->{mean_utr}->{$padding})];
		
	}
	 my $ii =0;
	 foreach my $i (@{$self->exons_index}){
	 
	# for (my $i=0;$i<@$m;$i++){
	 	$data->[$ii]->[$nbp] = $m->[$i];
	 	push(@{$self->{stats}->{$patient->id}},$mean->[$i]);
	 	$data_mean->[$ii]->[$nbp] = $mean->[$i];
	 	$error ++ if $m->[$i] < $self->limit  && $self->types_index->[$ii] & $self->mask_exons->{capture} && !($self->types_index->[$ii] & $self->mask_exons->{non_coding}); 
	 	$ii++;
	 }
	$nbp++;
	}
	$self->{min_matrix} = $data;
	$self->{mean_matrix} = $data_mean;
	$self->{error} = $error;


}



sub  html_table {
 my ($self) = @_;

 my $cgi =$self->cgi;
my $patients = $self->patients;
 my $data = $self->min_matrix;
 my $data_mean = $self->mean_matrix;
 my $tid = $self->transcript->id;
 	my $out;
 	$out.=$cgi->start_div({ style=>"text-align: center;vertical-align:middle;font-size: 9px;font-family:  Verdana;margin-bottom: 0px;height:100%"});
	$out.= $cgi->start_table({class=>"table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 9px;font-family:  Verdana;margin-bottom: 0px;"});
$out.= $cgi->start_Tr();
my $gene_name = $self->transcript->getGene->external_name;
my $tr_name = $self->transcript->name;
my $chr_name = $self->transcript->getChromosome()->ucsc_name;
my $strand = "fwd";
$strand = "rev" if $self->transcript->strand == -1 ;
my $pos = "[".$self->transcript->start."-".$self->transcript->end."]";
my @pnames;
my $click_gene = qq{load_graph_gene('$tr_name');};

	my $text = qq{<div class="btn   btn-xs btn-primary " style="font-size:9px;min-width:100px!important;position:relative;" '  style="border : 1px" onClick= "$click_gene">$gene_name <br> $tr_name  <br>$chr_name:$pos ($strand)</div>};
		$out.= $cgi->th({style=>"min-width:60px;padding:4px;"},$text);
	foreach my $patient (sort{$a->name cmp $b->name} @$patients){
			my $name = $patient->name(); 
			my $temp = $self->{stats}->{$patient->id};
			my $patient_mean = int sum(@{$temp}) / scalar (@{$temp});
			my $scolor="";
			if ($patient->sex == 2){
				$scolor="background-color:rgb(255, 110, 199);color:black";
			}
			my $click_patient = qq{load_graph_transcript('$name','$tr_name');};
			my $icon = $patient->return_icon;
			push(@pnames,$name);
			my $text = qq{<div class="btn   btn-xs btn-info " style="font-size:10px;min-width:100 px !important;position:relative;;border : 1px;$scolor"  onClick= "$click_patient"> $icon<br>$name<br>$patient_mean </div>};
			$out.= $cgi->th({style=>"min-width:75px !important;padding:4px;"},$text);
	}
	
$out.= $cgi->end_Tr();






my $capture  = $self->transcript->getChromosome()->getIntSpanCapture();

for (my $i = 0 ; $i< @{$self->min_matrix};$i++){
		my $type = $self->types_index->[$i];
		my $min_line = $self->min_matrix->[$i];	
		my $mean_line = $self->mean_matrix->[$i];	
		$out.= $cgi->start_Tr();
		
		my $exon = $self->transcript->getExon( $self->names_index->[$i]);
		my $ename = $exon->name;
		my $start = $exon->start;
		my $end = $exon->end;
			my $s1 = $exon->getGenomicSpan()->intersection($capture);
		my $mean_exon = int sum(@$mean_line)/scalar(@$mean_line);
		my $mean_min_exon = int sum(@$min_line)/scalar(@$min_line);;
		my $click_exon_patients = qq{load_graph_one_exon_all_patients('$tr_name','$ename');};
		my $text = qq{<div class="btn   btn-xs btn-info " style="font-size:9px;border:1px;background-color:#5B5EA6"  onClick= "$click_exon_patients">$ename ($mean_exon) [$start-$end] </div>};
		my $l1 = qq{	<span class="label label-xs label-success" style ="font-size:08px;min-width: 50px !important;display: inline-block !important;font-size:09px;color:black;;background-color:green;color:white;" >$mean_exon<span class="label label-xs label-primary" style ="font-size:07px;float:right;padding-left:3px;position: relative;top:-2px;right:-2px">$mean_min_exon</span></span>};
		
		$out.= $cgi->td({style=>"max-height: 5px !important;min-width:60px;padding:4px;"},$text);

		foreach  (my $j=0;$j < @$min_line;$j++){
			my $name = $pnames[$j];
			my $bgcolor;
			if ( !($type & $self->mask_exons->{capture})  ){
					$bgcolor = $self->grey_hex->[1];
			}
				elsif ( $type & $self->mask_exons->{non_coding} && $self->utr == 0){
				$bgcolor = $self->grey_hex->[0];
			
			}
				elsif ($type & $self->mask_exons->{duplicate} ){
					$bgcolor = $self->purple_hex->[0];
					if ($min_line->[$j] < $self->limit){
						$bgcolor = $self->purple_hex->[1];;
					}
					
				}
			elsif ($min_line->[$j] < $self->limit){
				my $rc = $self->return_level_red($min_line->[$j] ); 
				$bgcolor = $self->red_hex->[$rc];
				if ($self->transcript->getChromosome->name eq "Y"){
					$bgcolor = "#ff6ec7";
				}
			}
			else {
				my $gc = $self->return_level_green($min_line->[$j] );
				$bgcolor = $self->green_hex->[$gc];
					#$bgcolor = "transparent";
			}
			
			my $v1 = $mean_line->[$j];
			my $v2 = $min_line->[$j];
			
				my $click = qq{load_graph_exon('$name','$tid','$ename','$start','$end');};
			my $l2 = qq{	<sup><span class="badge label-success" style ="font-size:09px;">$v2</span></sup>};
			my $l1 = qq{	<span class="label label-xs label-success" style ="font-size:08px;min-width: 50px !important;display: inline-block !important;font-size:09px;color:black;;background-color:$bgcolor;color:white;" onClick= "$click">$v1<span class="label label-xs label-primary" style ="font-size:07px;float:right;padding-left:3px;position: relative;top:-2px;right:-2px">$v2</span></span>};
			
			$out.= $cgi->td({style=>"max-height: 5px !important;min-width:60px;padding:2px;"},$l1);
		}
		$out.= $cgi->end_Tr();
 	}
	$out .= $cgi->end_table();
		$out .= $cgi->end_div();
	return $out;
}

sub return_level_red {
	my ($self,$min) = @_;
	
	
	my $max3 = $self->limit - (0.1 * $self->limit);
	if (between($min,$max3,$self->limit) ){
			return 0;		
		}	
		my $max2 = $self->limit - (0.25 * $self->limit);
		if (between($min,$max2,$self->limit) ){
			return 1;
		}		
			return 2;
}

sub return_level_green {
	my ($self,$min) = @_;
	my $limit = $self->limit;
	my $max4 = $limit + (0.1 * $limit);
		if (between($min,$limit,$max4) ){
			return 0;
		}	
			my $max5 = $limit + (0.25 * $limit);
		if (between($min,$limit,$max5) ){
			return 1;
		}	
		return 2;
}

sub image {
	 my ($self) = @_;
	 
 my $limit = $self->limit;
my $patients = $self->patients;
my $data = $self->min_matrix;
	 my $nb_line =15;
my $size = 2;

my $x = 0;
my $y = 0;

my $xr = 255;
my $xg = 10;
my $xb = 0;

my $yr = 0;
my $yg = 255;
my $yb = 0;

my $max = 20;
my $hcolors;
my $max1 = $limit * 2 ;
my $max2 = $limit - (0.25 * $limit);
my $max3 = $limit - (0.1 * $limit);


my $max4 = $limit + (0.1 * $limit);
my $max5 = $limit + (0.25 * $limit);


my $nb_col = scalar(@$patients);
my $nb_row = scalar(@$data);
my $w =$nb_col*($size+2)+1;
my $h = ($size+2)*$nb_row+1;

my $image= new GD::Image($w+1,$h+1);
		
my       $black = $image->colorAllocate(0,0,0);       


my @red_img;
if ($self->transcript->getChromosome->name eq "Y") {
	foreach my $a (@{$self->red_rgb}){
	push(@red_img,$image->colorAllocate(255, 110, 199));
	}
}
else {
foreach my $a (@{$self->red_rgb}){
	push(@red_img,$image->colorAllocate(@$a));
}
}
my @green_img;

foreach my $a (@{$self->green_rgb}){
	push(@green_img,$image->colorAllocate(@$a));
}

my @purple_img;

foreach my $a (@{$self->purple_rgb}){
	push(@purple_img,$image->colorAllocate(@$a));
}

my @grey_img;
foreach my $a (@{$self->grey_rgb}){
	push(@grey_img,$image->colorAllocate(@$a));
}

$y =2;
	my $nbe = 0;
	#foreach my $exon (sort{$a->end*$a->strand <=> $b->end*$b->strand } @$exons) {
	foreach (my $i=0;$i<@$data;$i++)	{
		my $matrix = $data->[$i];
		$x=2;
		#my $min_data = $hexon->{min};
		my @min_data = @$matrix;
		my $type = $self->types_index->[$i];
		foreach (my $j =0;$j<@$matrix;$j++) {
			my $min = $matrix->[$j];
		
			my $bgcolor;

	if ( !($type & $self->mask_exons->{capture})  ){
				#945251
				$bgcolor = $grey_img[1];
			
			}		
		elsif ( $type & $self->mask_exons->{non_coding}  && $self->utr == 0 ){
				
				$bgcolor =  $grey_img[0];
			
			}
				elsif ($type & $self->mask_exons->{duplicate} ){
					if ($min < $self->limit){
							$bgcolor =  $purple_img[1];
					}
					else {
							$bgcolor =  $purple_img[0];
					}				
			
			}
			elsif ($min < $self->limit){
				my $rc = $self->return_level_red($min ); 
				$bgcolor = $red_img[$rc];
			}
			else {
				my $gc = $self->return_level_green($min );
				$bgcolor = $green_img[$gc];
			}
			
	
	
		
		$image->filledRectangle($x,$y,$x+$size,$y+$size,$bgcolor);
		$x+=$size+2
		} #end patient
		$y+=$size+2;
		$image->rectangle(0,0,$w,$h,$black);
		$nbe ++;
	}#end exon

			
	return $image;

}


1;
 