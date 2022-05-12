package image_coverage;
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

#method image (ArrayRef :$patients, GenBoTranscript :$transcript, Int :$intronic, Int :$utr,  Int :$padding,  Int :$limit ){
	
sub gradient {
    my ( $min, $max,$num ) = @_;

    my $middle = ( $min + $max ) / 2;
    my $scale = 254 / ( $middle - $min );

    return sub {
        my $num = shift;
        return [254,0,0] if $num <= $min;    # lower boundry
        return [0,254,0] if $num >= $max;    # upper boundary

        if ( $num < $middle ) {
        	return [254,int( ( $num - $min ) * $scale ),0];
           # return sprintf "FF%02X00" => int( ( $num - $min ) * $scale );
        }
        else {
            return [254 - int( ( $num - $middle ) * $scale) ,254,0];
            #  sprintf "%02XFF00" => 255 - int( ( $num - $middle ) * $scale );
        }
    };
}
 
sub image_cnv {
	 my ($patients,$transcript,$cgi) = @_;
	 
	 return image_cnv_PROD($patients,$transcript,$cgi) ;
	 
	  if ($transcript->project->isCache()){
	  	image_cnv_no_cache($patients,$transcript,$cgi) ;
	  }
	  else {
	  	image_cnv_no_cache($patients,$transcript,$cgi);
	  }
	# image_cache_cnv($patients,$transcript,$cgi) if $transcript->project->isCache();
	 #image_cnv_no_cache($patients,$transcript,$cgi);
}


sub image_cnv_PROD {
	my ($patients,$transcript,$cgi) = @_;

my $nb_line =15;
my $size = 3;

my $x = 0;
my $y = 0;

my $xr = 255;
my $xg = 10;
my $xb = 0;

my $yr = 0;
my $yg = 255;
my $yb = 0;


my $no =  $transcript->project->noSqlCoverage();
$transcript->getChromosome->setPrimersForPatient($patients);
my @primers = sort{$a->end*$a->strand <=> $b->end*$b->strand } @{$transcript->getPrimers()};

my $nb_col = scalar(@$patients);
my $w =$nb_col*($size+2)+2;
my $h = ($size+2)*scalar(@primers)+2;

my $gradient = gradient( 1, 50);

my $image= new GD::Image($w+1,$h+1);
		


my        $grey1 = $image->colorAllocate(90,90,90);	
my        $grey2 = $image->colorAllocate(140,140,140);
my        $grey3 = $image->colorAllocate(0, 255, 0);
my        $blue = $image->colorAllocate(0,0,255);
my        $red = $image->colorAllocate(255,10,50);   
my       $white = $image->colorAllocate(255,255,255);
my       $black = $image->colorAllocate(0,0,0);
my $hcolors;
my $zz = ($size+2)*scalar(@primers)+1;
#$image->rectangle(1,1,$nb_col*($size+2)+1,($size+2)*scalar(@$primers)+1,$black);
my $data={};
$y =2;
my @color_background  = (223, 255, 216);
#foreach my $primer (sort{$a->end*$a->strand <=> $b->end*$b->strand } @$primers) {
#		foreach my $patient (sort{$a->name cmp $b->name} @$patients){
#			next if $patient->name() ne "HYP_2684";
#			
#			my $level = $primer->level($patient);
#			my $score = $primer->cnv_score($patient);
#			
#			warn $score." $level ".$patient->name." ".$level." ".$patient->count($primer->id)."\n";
#		}
#		warn "---";
#}


my $hscore;
#foreach my $patient (sort{$a->name cmp $b->name} @$patients){
#	my @ids = map{$_->id} @primers;
#	 $hscore->{$patient->name} = $no->get_bulk($patient->name,\@ids);
#	#warn Dumper $h;
#}
foreach my $patient (sort{$a->name cmp $b->name} @$patients){
	
	 $hscore->{$patient->name} = $no->get_all($patient->name()."_cnv");
	
}
foreach my $patient (sort{$a->name cmp $b->name} @$patients){
	my $pscore =0;
	my $maxscore = 0;
	my $previous = -5;
	for (my $i=0;$i<@primers;$i++){
	my $primer = $primers[$i];
	my $hash = $hscore->{$patient->name}->{$primer->id};
	die($primer->id) unless $hash;
	#my $h = $primer->level($patient,1);
	
	#warn Dumper $h; 
	#die();
	#foreach my $primer (sort{$a->end*$a->strand <=> $b->end*$b->strand } @primers) {
	#foreach my $primer (sort{$a->end*$a->strand <=> $b->end*$b->strand } @primers) {
		
		$maxscore = $pscore if $pscore > $maxscore;
		#my $hash = $hscore->{$patient->name}->{$primer->id};
		my ($score,$level);
		if ($hash){
			$score = $hash->{cnv};
			$level = $hash->{level_cnv}; 
			$hscore->{$patient->name}->{$primer->id}->{cnv} = $score;
			$hscore->{$patient->name}->{$primer->id}->{level_cnv} = $level;
			
		}
		else {
		 	$score = $primer->cnv_score($patient);
			$level = $primer->level($patient);
			$hscore->{$patient->name}->{$primer->id}->{cnv} = $score;
			$hscore->{$patient->name}->{$primer->id}->{level_cnv} = $level;
		}			
			$level = 1 if ($level == 1 || $level == 4);
			$level = 2 if ($level == 2 || $level == 3);
			
		if ($previous eq $level && $level>0){
			$pscore ++;
			next;
		}
		elsif  ($level eq 0 && abs(1-$score) > 0.2 ){
			my $start =-1;
			
			$start =0 if $i ==0;
			my $end = 2;
			my $mean =0;
			my $nb =0;
			my $nb_level_del;
			my $nb_level_dup;
			for   ( my $z = $start;$z<$end;$z++){
				my $id = $i+$z;
				next unless $primers[$id];
				next if $hscore->{$patient->name}->{$primers[$id]->id}->{level_cnv} < 0;
				 $nb ++;
				 $mean += $hscore->{$patient->name}->{$primers[$id]->id}->{cnv};
				 $nb_level_del++ if $hscore->{$patient->name}->{$primers[$id]->id}->{level_cnv} ==4;
				  $nb_level_dup++ if $hscore->{$patient->name}->{$primers[$id]->id}->{level_cnv} ==3;
				  
			}
			
			$mean /=$nb;
			$hscore->{$patient->name}->{$primer->id}->{cnv} = $score;
			$hscore->{$patient->name}->{$primer->id}->{level_cnv} = $level;
			if ($mean > 1.3){
			$hscore->{$patient->name}->{$primer->id}->{level_cnv}  = 3;
			}
			elsif  ($mean < 0.6){
				$hscore->{$patient->name}->{$primer->id}->{level_cnv}  = 4;
			}
			if (($score < 0.6) && ($previous == 2 )) {
				#$pscore +=0.5;
				next;
			}
			if (($score > 1.4) && ($previous == 1)) {
				#$pscore +=0.5;
				next;
			}
		}
		$maxscore = $pscore if $pscore > $maxscore;
		$pscore =0;
		$previous = $level;
	}
	@color_background =(127, 255, 27) if ($maxscore > 3);
	
}
#die();
	foreach my $primer ( @primers) {
			$x =2;
			foreach my $patient (sort{$a->name cmp $b->name} @$patients){
				my @colors =();
		#	my $zscore = $primer->zscore($patient);
			
			#my $level = $primer->level($patient,1);#$hscore->{$patient->name}->{$primer->id}->{level_cnv} ;
				my $level = $hscore->{$patient->name}->{$primer->id}->{level_cnv} ;
			#warn $level;
		#	warn $primer->level($patient,1)."-".$level;
			my $score = $hscore->{$patient->name}->{$primer->id}->{cnv} ;
			my $name = $patient->name;
			@colors = (223, 255, 216);
			#@colors =(172, 223, 172);
			@colors =@color_background;
			
			my $colors;
			if ($level == -2 ){
				@colors=(90,90,90);	
				#$color = $grey1;
			}
			elsif ($level == -1) {
				@colors=(140,140,140);
				}
				elsif ($level == 2){
						@colors =(52,152,255);
						@colors =(2, 71, 253);
				}
				elsif ($level==3){
					my $colorr =int (150*($score-1))+105;
					@colors =(52,152,255);
					@colors =(2, 71, 253);
				}
				elsif ($level == 1){
						@colors =(229, 40, 76);
						@colors =(255, 64, 64);
				}
				elsif ($level==4){
					my $colorr =int (50*(0.6-$score))+205;
					@colors =(229, 40, 76);
					@colors =(255, 64, 64);
					
				}
				elsif ($level ==0){
					#@colors =(205, 223, 245);
					#FFC9CF
					if ($score <=0.55){
						@colors = (254, 138, 53)
						#@colors =(222, 165, 164);
					}
						elsif ($score <=0.6){
						@colors =(246, 179, 134);
#						#@colors =(254, 209, 219);
					}
					elsif ($score <=0.7){
						@colors =(253, 218, 178);
						#rgb(253, 218, 178)
#						#@colors =(254, 209, 219);
					}
#					elsif ($score <=0.7){
#						@colors =(135, 255, 117);
#						#@colors =(254, 209, 219);
#					}
					elsif ($score >= 1.5){
						@colors = (115, 178, 217)
						#@colors =(254, 209, 219);
					}
					elsif ($score >= 1.4){
					@colors =(205, 223, 245);
				}
					
					
				}
				my $id_color = join(",",@colors);
				unless (exists $hcolors->{$id_color}){
					$hcolors->{$id_color} = $image->colorAllocate(@colors);
				}
			
				$image->rectangle($x-1,$y-1,$x+$nb_col*($size+2)+1,$y+$size+1,$black);
				$image->filledRectangle($x,$y,$x+$nb_col*($size+2),$y+$size,$hcolors->{$id_color} );
				$data->{$primer->id}->{$name}->{colors} = \@colors;
				$data->{$primer->id}->{$name}->{text} = $level;
				$data->{$primer->id}->{$name}->{cnv_score} = $score;
				$x+=$size+2
			}
		$y+=$size+2;
	}
	#next unless $exon->{exon} ==1;
$image->rectangle(0,0,$w,$h,$black);

my $re;
$re->{image} = $image;
$re->{data} = $data;
return $re;


}

sub image_cnv_no_cache {
	 my ($patients,$transcript,$cgi) = @_;

my $nb_line =15;
my $size = 3;

my $x = 0;
my $y = 0;

my $xr = 255;
my $xg = 10;
my $xb = 0;

my $yr = 0;
my $yg = 255;
my $yb = 0;

my $no =  $transcript->project->noSqlCoverage();
my $no1 =   $transcript->project->noSqlCnvs("r");
#$transcript->getChromosome()->getPrimers();


#$transcript->getChromosome->setPrimersForPatient($patients);

my @primers = sort{$a->end*$a->strand <=> $b->end*$b->strand } @{$transcript->getPrimers()};
my $nb_col = scalar(@$patients);
my $w =$nb_col*($size+2)+2;
my $h = ($size+2)*scalar(@primers)+2;

my $gradient = gradient( 1, 50);

my $image= new GD::Image($w+1,$h+1);


my        $grey1 = $image->colorAllocate(90,90,90);	
my        $grey2 = $image->colorAllocate(140,140,140);
my        $grey3 = $image->colorAllocate(0, 255, 0);
my        $blue = $image->colorAllocate(0,0,255);
my        $red = $image->colorAllocate(255,10,50);   
my       $white = $image->colorAllocate(255,255,255);
my       $black = $image->colorAllocate(0,0,0);
my $hcolors;
my $zz = ($size+2)*scalar(@primers)+1;
#$image->rectangle(1,1,$nb_col*($size+2)+1,($size+2)*scalar(@$primers)+1,$black);
my $data={};
$y =2;
my @color_background  = (223, 255, 216);



my $hscore;

for (my $i=0;$i<@primers;$i++){
	my $primer = $primers[$i];
	foreach my $patient (sort{$a->name cmp $b->name} @$patients){
		$hscore->{$patient->name}->{$primer->id}->{cnv} = $primer->cnv_score($patient);;
			$hscore->{$patient->name}->{$primer->id}->{level_cnv} = $primer->level($patient);;
	}
}
foreach my $patient (sort{$a->name cmp $b->name} @$patients){
	my $pscore =0;
	my $maxscore = 0;
	my $previous = -5;
	for (my $i=0;$i<@primers;$i++){
	my $primer = $primers[$i];
	#$primer->getPatients();
	my $hash ;#= $hscore->{$patient->name}->{$primer->id};
#	die($primer->id) unless $hash;
	#my $h = $primer->level($patient,1);
	
	#warn Dumper $h; 
	#die();
	#foreach my $primer (sort{$a->end*$a->strand <=> $b->end*$b->strand } @primers) {
	#foreach my $primer (sort{$a->end*$a->strand <=> $b->end*$b->strand } @primers) {
		
		$maxscore = $pscore if $pscore > $maxscore;
		#my $hash = $hscore->{$patient->name}->{$primer->id};
		my ($score,$level);
		$score = $primer->cnv_score($patient);
		$level = $primer->level($patient);;

			$level = 1 if ($level == 1 || $level == 4);
			$level = 2 if ($level == 2 || $level == 3);
		if ($previous eq $level && $level>0){
			$pscore ++;
			next;
		}
		elsif  ($level eq 0 && abs(1-$score) > 0.2 ){
			my $start =-1;
			
			$start =0 if $i ==0;
			my $end = 2;
			my $mean =0;
			my $nb =0;
			my $nb_level_del;
			my $nb_level_dup;
			for   ( my $z = $start;$z<$end;$z++){
				my $id = $i+$z;
				next unless $primers[$id];
				
				next if $hscore->{$patient->name}->{$primers[$id]->id}->{level_cnv} < 0;
				 $nb ++;
				 $mean += $hscore->{$patient->name}->{$primers[$id]->id}->{cnv};
				 $nb_level_del++ if $hscore->{$patient->name}->{$primers[$id]->id}->{level_cnv} ==4;
				  $nb_level_dup++ if $hscore->{$patient->name}->{$primers[$id]->id}->{level_cnv} ==3;
				  
			}
			
			$mean /=$nb;
			$hscore->{$patient->name}->{$primer->id}->{cnv} = $score;
			$hscore->{$patient->name}->{$primer->id}->{level_cnv} = $level;
			if ($mean > 1.3){
			$hscore->{$patient->name}->{$primer->id}->{level_cnv}  = 3;
			}
			elsif  ($mean < 0.6){
				$hscore->{$patient->name}->{$primer->id}->{level_cnv}  = 4;
			}
			if (($score < 0.6) && ($previous == 2 )) {
				#$pscore +=0.5;
				next;
			}
			if (($score > 1.4) && ($previous == 1)) {
				#$pscore +=0.5;
				next;
			}
		}
		$maxscore = $pscore if $pscore > $maxscore;
		$pscore =0;
		$previous = $level;
	}
	@color_background =(127, 255, 27) if ($maxscore > 3);
	
}
#die();
	foreach my $primer ( @primers) {
			$x =2;
			foreach my $patient (sort{$a->name cmp $b->name} @$patients){
				my @colors =();
		#	my $zscore = $primer->zscore($patient);
			
			#my $level = $primer->level($patient,1);#$hscore->{$patient->name}->{$primer->id}->{level_cnv} ;
				my $level = $hscore->{$patient->name}->{$primer->id}->{level_cnv} ;
				
			#warn $level;
		#	warn $primer->level($patient,1)."-".$level;
			my $score = $hscore->{$patient->name}->{$primer->id}->{cnv} ;
			my $name = $patient->name;
			@colors = (223, 255, 216);
			#@colors =(172, 223, 172);
			@colors =@color_background;
#			warn $level." ".$score." $name";
			my $colors;
			if ($level == -2 ){
				@colors=(90,90,90);	
				#$color = $grey1;
			}
			elsif ($level == -1) {
				@colors=(140,140,140);
				}
				elsif ($level == 2){
						@colors =(52,152,255);
						@colors =(2, 71, 253);
				}
				elsif ($level==3){
					my $colorr =int (150*($score-1))+105;
					@colors =(52,152,255);
					@colors =(2, 71, 253);
				}
				elsif ($level == 1){
						@colors =(229, 40, 76);
						@colors =(255, 64, 64);
				}
				elsif ($level==4){
					my $colorr =int (50*(0.6-$score))+205;
					@colors =(229, 40, 76);
					#@colors =(255, 64, 64);
					
				}
				elsif ($level ==0){
					#@colors =(205, 223, 245);
					#FFC9CF
					if ($score <=0.1){
						@colors = (0, 0,0)
						#@colors =(222, 165, 164);
					}
					elsif ($score <=0.55){
						@colors = (254, 138, 53)
						#@colors =(222, 165, 164);
					}
						elsif ($score <=0.6){
						@colors =(246, 179, 134);
#						#@colors =(254, 209, 219);
					}
					elsif ($score <=0.7){
						@colors =(253, 218, 178);
						#rgb(253, 218, 178)
#						#@colors =(254, 209, 219);
					}
#					elsif ($score <=0.7){
#						@colors =(135, 255, 117);
#						#@colors =(254, 209, 219);
#					}
					elsif ($score >= 1.5){
						@colors = (115, 178, 217)
						#@colors =(254, 209, 219);
					}
					elsif ($score >= 1.4){
					@colors =(205, 223, 245);
				}
					
					
				}
				my $id_color = join(",",@colors);
				unless (exists $hcolors->{$id_color}){
					$hcolors->{$id_color} = $image->colorAllocate(@colors);
				}
			
				$image->rectangle($x-1,$y-1,$x+$nb_col*($size+2)+1,$y+$size+1,$black);
				$image->filledRectangle($x,$y,$x+$nb_col*($size+2),$y+$size,$hcolors->{$id_color} );
				$data->{$primer->id}->{$name}->{colors} = \@colors;
				$data->{$primer->id}->{$name}->{text} = $score;
				my $cai_count= $no1->get("raw_data",$primer->id."_cai");
				$data->{$primer->id}->{$name}->{cnv_score} = $score."++".$patient->meanDepth($primer->getChromosome->name,$primer->start,$primer->end).$level."<BR>".join(":",values %$cai_count);
				if ($level == -1){
#					warn Dumper $primer->{zscore_data};
#					foreach my $p (@{$primer->getPatients}){
#						warn $p->name;
#						warn $p->name." ".$p->meanDepth($primer->getChromosome->name,$primer->start,$primer->end);
#					}
			#		warn $primer->id;
				#	die();
				}
				$x+=$size+2
			}
		$y+=$size+2;
	}
	#next unless $exon->{exon} ==1;
$image->rectangle(0,0,$w,$h,$black);

my $re;
$re->{image} = $image;
$re->{data} = $data;
return $re;
}

sub produce_matrix_data2 {
	my ($patients,$transcript,$utr,$padding) = @_;
	$padding = 20;
	my $exons ;
	my $intronic=1;
if ($intronic == 0 ){
	$exons  = $transcript->getAllGenomicsParts();
}
else {
	 $exons = $transcript->getExons();
}
my $data;
my $data2;

foreach my $exon (sort{$a->end*$a->strand <=> $b->end*$b->strand } @$exons) {
	my $hexon;
	my $data_patient_min;
	my $data_patient_mean;
	foreach my $patient (sort{$a->name cmp $b->name} @$patients){
		my $hpatient;
		
		$hpatient->{id} = $patient->id;
		
			my $pos = $exon->return_start_end_no_utr(padding=>50);
			my $xa = [];
			my $min =0;
			my $mean =0;
			
			if ($pos){
			 $xa =  $patient->depth($exon->getChromosome->name,$pos->{start},$pos->{end}) if $pos;
			 $min = min(@$xa);
			 $mean = sum(@$xa)/scalar(@$xa);

			}	
			
			push(@$data_patient_mean,$mean);
			push(@$data_patient_min,$min);
	
	}
	my $hexons;
		my $text   = pack("w".scalar(@$patients),@$data_patient_min );
		$hexons->{min_pack} = $text ;
		
		#push (@{$hexons->{min}},@$data_patient_min );
		#push (@{$hexons->{mean}},@$data_patient_mean );
		push(@{$data2->{matrix}},$hexons);
		
}

foreach my $exon (sort{$a->end*$a->strand <=> $b->end*$b->strand } @$exons) {
	push (@{$data2->{exons}},$exon->id);
}
	foreach my $patient (sort{$a->name cmp $b->name} @$patients){
			push (@{$data2->{patients}},$patient->name);
	}

return $data2;
#return $data;	
}




sub image_depth_lmdb_from_matrix2 {
	 my ($data,$limit) = @_;
	 my $nb_line =15;
my $size = 3;

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
my $max2 = $limit + 0.5 * $limit;
my $max3 = $limit + 0.2 * $limit;

my $nb_col = scalar(@{$data->{patients}});

my $w =$nb_col*($size+2)+2;
my $h = ($size+2)*scalar(@{$data->{exons}})+2;

my $gradient = gradient( 1, $limit);

my $image= new GD::Image($w+1,$h+1);
		
my       $white = $image->colorAllocate(255,255,255);
my       $black = $image->colorAllocate(0,0,0);       
my        $red = $image->colorAllocate(255,10,50);    
my        $yellow = $image->colorAllocate(255,255,0);    
my        $blue = $image->colorAllocate(0,0,255);
my        $green = $image->colorAllocate(0,255,0);
my        $green2 = $image->colorAllocate(136,235,136);
my        $green3 = $image->colorAllocate(162,255,162);
my        $green4 = $image->colorAllocate(239,255,239);
my        $red2 = $image->colorAllocate(255,239,239);
my        $grey = $image->colorAllocate(200,200,200);	
my        $grey2 = $image->colorAllocate(250,250,250);
my        $purple = $image->colorAllocate(214,119,253);

$y =2;
	foreach my $hexon (@{$data->{matrix}}){
		$x=2;
		#my $min_data = $hexon->{min};
		my @min_data = unpack("w".$nb_col,$hexon->{min_pack});
		my $i=0;
		foreach my $min (@min_data) {
			$i++;
				my $color2 = $green;
			my $color = $green;
			my $colors = [];	
			$min =0 unless $min;
		if ($min <$limit){
			$color = $red;
			$color2= $yellow;
		}		
		elsif ($min <$limit){
			 
			
				
			my $ar = $gradient->($min);
			my $id_color = join(";",@$ar);
			$hcolors->{$id_color} = $image->colorAllocate(@$ar) unless exists $hcolors->{$id_color}  ;
			#warn $patient->id." ar  ".join(";",@$ar);
			$color = $hcolors->{$id_color} ;
		}
		
		elsif ($min>=$limit){
			$color = $green;
		}
		else {
			
			my $ar = $gradient->($min);
			my $id_color = join(";",@$ar);
			$hcolors->{$id_color} = $image->colorAllocate(@$ar) unless exists $hcolors->{$id_color}  ;
			$color = $hcolors->{$id_color} ;
		}
		
		$image->filledRectangle($x,$y,$x+$size,$y+$size,$color);
		$x+=$size+2
		} #end patient
		$y+=$size+2;
		$image->rectangle(0,0,$w,$h,$black);
	}#end exon

			
	return $image;

}



sub image_depth_lmdb {
 my ($patients,$transcript,$intronic,$utr,$padding,$limit,$all,$minimum) = @_;
my $nb_line =15;
my $size = 3;

my $x = 0;
my $y = 0;

my $xr = 255;
my $xg = 10;
my $xb = 0;

my $yr = 0;
my $yg = 255;
my $yb = 0;

my $max = 20;

my $max1 = $limit * 2 ;
my $max2 = $limit + 0.5 * $limit;
my $max3 = $limit + 0.2 * $limit;

my $exons ;
if ($intronic == 1 ){
	$exons  = $transcript->getAllGenomicsParts();
}
else {
	 $exons = $transcript->getExons();
}


my $nb_col = scalar(@$patients);

my $w =$nb_col*($size+2)+2;
my $h = ($size+2)*scalar(@$exons)+2;

my $gradient = gradient( 1, $limit);

my $image= new GD::Image($w+1,$h+1);
		
my       $white = $image->colorAllocate(255,255,255);
my       $black = $image->colorAllocate(0,0,0);       
my        $red = $image->colorAllocate(255,10,50);    
my        $yellow = $image->colorAllocate(255,255,0);    
my        $blue = $image->colorAllocate(0,0,255);
my        $green = $image->colorAllocate(0,255,0);
my        $green2 = $image->colorAllocate(136,235,136);
my        $green3 = $image->colorAllocate(162,255,162);
my        $green4 = $image->colorAllocate(239,255,239);
my        $red2 = $image->colorAllocate(255,239,239);
my        $grey = $image->colorAllocate(200,200,200);	
my        $grey2 = $image->colorAllocate(250,250,250);
my        $purple = $image->colorAllocate(214,119,253);

my $hcolors;
#$image->rectangle(1,1,$nb_col*($size+2)+1,($size+2)*scalar(@$exons)+1,$black);
	my $data;
	my $alert;
	$y =2;
	my $capture_intspan = $transcript->getChromosome->getIntSpanCapture();
		my $no =  $transcript->project->noSqlCoverage();
		my $raw_data;
		
	foreach my $exon (sort{$a->end*$a->strand <=> $b->end*$b->strand } @$exons) {
		
	$x =2;
	my $s1 = $exon->getGenomicSpan()->intersection($capture_intspan);
	my $covered =1;
	if ($exon->isExon()){
			$covered =0 if $s1->is_empty;
		}
		else {
			$covered =0 if scalar($s1->as_array)<20;
		}
	unless ($exon->isExon){
		$image->filledRectangle($x,$y,$x+$nb_col*($size+2),$y+$size,$grey2);
	}
 	
  
	foreach my $patient (sort{$a->name cmp $b->name} @$patients){
	 	my $spanTr = $no->get($patients->[0]->name,$transcript->id."_spandup");	
	 	$spanTr = Set::IntSpan::Fast::XS->new() unless $spanTr; 
		my $dir = $patient->project->getVariationsDir("duplicate_region_calling")."/regions";
		my $filebed =  $patient->project->getVariationsDir("duplicate_region_calling")."/regions/".$patient->name().".dup.bed";
		my $span =$exon->getGenomicSpan->intersection($spanTr);#->intersection($exon->getGenomicSpan);
		
		if ($exon->isExon() && $exon->is_noncoding && $utr ==1 ){
			next;
		}
		my $pid =$patient->id;
		my ($mean,$intspan,$min);

			
			my $pos = $exon->return_start_end_no_utr(padding=>$padding);
			my $xa = [];
			
			
			if ($pos){
			 $xa =  $patient->depth($exon->getChromosome->name,$pos->{start},$pos->{end}) if $pos;

			$min = min(@$xa);
			$mean = sum(@$xa)/scalar(@$xa);

		}

			$data->{$exon->id}->{$patient->id}->{mean} = $mean;
			$data->{$exon->id}->{$patient->id}->{min} = $min;
			$data->{$exon->id}->{$patient->id}->{covered} = $covered;
			$data->{$exon->id}->{$patient->id}->{red} = 0;
			$data->{$exon->id}->{$patient->id}->{green} = 1;

			my $color2 = $green;
			my $color = $green;
			my $colors = [];
		unless ($span->is_empty){
			my $gc =  $no->get($patient->name,$transcript->id."_raw");
			#warn Dumper ($gc);
			 my $cover = $gc->coverage($exon->start,$exon->end);
			  $data->{$exon->id}->{$patient->id}->{raw_min} = $cover->{min};
			 $data->{$exon->id}->{$patient->id}->{raw_mean} = $cover->{mean};
			if ($exon->isExon() && $exon->is_noncoding && $utr ne 1 && $exon->intspan_no_utr->is_empty) {
				$color = $grey;
				$color2=$grey;
				 $data->{$exon->id}->{$patient->id}->{green}=0;
			}
			elsif ($cover->{min} < $limit) {
			$color = $red;
			$color2= $yellow;
			$data->{$exon->id}->{$patient->id}->{red} = 1;
			 $data->{$exon->id}->{$patient->id}->{green}=0;
			}
			else {
			$color = $purple;
			$color2=$purple;
			 $data->{$exon->id}->{$patient->id}->{green}=0;
			
			}
		}	
			
		elsif ($exon->isExon() && $exon->is_noncoding && $utr ne 1 && $exon->intspan_no_utr->is_empty) {
			$color = $grey;
			$color2=$grey;
			 $data->{$exon->id}->{$patient->id}->{green}=0;
		}
		elsif ($covered == 0){
			$color = $grey;
			$color2=$grey;
			 $data->{$exon->id}->{$patient->id}->{green}=0;
		}
		elsif ($min <$limit){
			$alert =1;
			$color = $red;
			$color2= $yellow;
			$data->{$exon->id}->{$patient->id}->{red} = 1;
			 $data->{$exon->id}->{$patient->id}->{green}=0;
			
		}
		elsif ($min>50){
			$color = $green;
		}
		else {
			my $ar = $gradient->($min);
			my $id_color = join(";",@$ar);
			$hcolors->{$id_color} = $image->colorAllocate(@$ar) unless exists $hcolors->{$id_color}  ;
			#warn $patient->id." ar  ".join(";",@$ar);
			$color = $hcolors->{$id_color} ;
			#	warn $patient->id." color ".join(";",$image->rgb($color));
		}
		

		$data->{$exon->id}->{$patient->id}->{color} = [$image->rgb($color)];
	
		$data->{$exon->id}->{$patient->id}->{color2} = [$image->rgb($color2)];
		$image->filledRectangle($x,$y,$x+$size,$y+$size,$color);
		$x+=$size+2
		} #end patient
		$y+=$size+2;
	
	}#end exon
$image->rectangle(0,0,$w,$h,$black);
binmode STDOUT;

my $re;

$re->{alert} =1 if $alert;
$re->{image} = $image;

$re->{data} = $data;
return $re;


}
	 
	
sub image {
 my ($patients,$transcript,$intronic,$utr,$padding,$limit,$all) = @_;
# die();
 $limit = 15;
my $nb_line =15;
my $size = 3;
my $alert;
my $x = 0;
my $y = 0;

my $xr = 255;
my $xg = 10;
my $xb = 0;

my $yr = 0;
my $yg = 255;
my $yb = 0;

my $max = 50;

my $max1 = $limit * 2 ;
my $max2 = $limit + 0.5 * $limit;
my $max3 = $limit + 0.2 * $limit;

my $exons ;
if ($intronic == 1 ){
	$exons  = $transcript->getAllGenomicsParts();
}
else {
	 $exons = $transcript->getExons();
}


my $nb_col = scalar(@$patients);

my $w =$nb_col*($size+2)+2;
my $h = ($size+2)*scalar(@$exons)+2;

my $gradient = gradient( 1, 50);

my $image= new GD::Image($w+1,$h+1);
		
my       $white = $image->colorAllocate(255,255,255);
my       $black = $image->colorAllocate(0,0,0);       
my        $red = $image->colorAllocate(255,10,50);    
my        $yellow = $image->colorAllocate(255,255,0);    
my        $blue = $image->colorAllocate(0,0,255);
my        $green = $image->colorAllocate(0,255,0);
my        $green2 = $image->colorAllocate(136,235,136);
my        $green3 = $image->colorAllocate(162,255,162);
my        $green4 = $image->colorAllocate(239,255,239);
my        $red2 = $image->colorAllocate(255,239,239);
my        $grey = $image->colorAllocate(200,200,200);	
my        $grey2 = $image->colorAllocate(250,250,250);
my        $purple = $image->colorAllocate(214,119,253);

my $hcolors;
#$image->rectangle(1,1,$nb_col*($size+2)+1,($size+2)*scalar(@$exons)+1,$black);
	my $data;
	$y =2;
	my $capture_intspan = $transcript->getChromosome->getIntSpanCapture();
		my $no =  $transcript->project->noSqlCoverage();
		
	foreach my $exon (sort{$a->end*$a->strand <=> $b->end*$b->strand } @$exons) {
	$x =2;
	my $s1 = $exon->getGenomicSpan()->intersection($capture_intspan);
	my $covered =1;
	if ($exon->isExon()){
			$covered =0 if $s1->is_empty;
		}
		else {
			$covered =0 if scalar($s1->as_array)<20;
		}
#	warn $exon->name();
	unless ($exon->isExon){
		$image->filledRectangle($x,$y,$x+$nb_col*($size+2),$y+$size,$grey2);
	}
	#next unless $exon->{exon} ==1;
  my $spanTr = $no->get($patients->[0]->name,$transcript->id."_spandup");	
 
  
	foreach my $patient (sort{$a->name cmp $b->name} @$patients){
	
	#	warn $spanTr;
		$spanTr = Set::IntSpan::Fast::XS->new() unless $spanTr; 
		#warn $spanTr->as_string()." ".$patient->name();
#		warn $patient->name();
		my $dir = $patient->project->getVariationsDir("duplicate_region_calling")."/regions";
		my $filebed =  $patient->project->getVariationsDir("duplicate_region_calling")."/regions/".$patient->name().".dup.bed";
#system("mkdir $dir && chmod a+rwx $dir" ) unless -e $dir;
	

		my $pid =$patient->id;
		my ($mean,$intspan,$min);
		#my $t  =$exon->cached_statistic_coverage_coding(patient=>$patient,padding=>$padding,limit=>$limit,utr=>$utr);  
		my $t  = $no->get($patient->name,$exon->id."_".$padding."_".$utr);
		if ($t){
			$mean = $t->{mean};
			$min = $t->{min};
			#die($mean);
		}
		else {
		 	confess();
		}
		
	
		 #	$patient->getCoverage($exon->getChromosome->ucsc_name,$exon->start,$exon->end);
			my $val = $min>$max ? $max : $min;
			$val -= $limit;
			#$val = 100 - $val;
			$data->{$exon->id}->{$patient->id}->{mean} = $mean;
			$data->{$exon->id}->{$patient->id}->{min} = $min;
			$data->{$exon->id}->{$patient->id}->{covered} = $covered;
			$data->{$exon->id}->{$patient->id}->{red} = 0;
			$data->{$exon->id}->{$patient->id}->{green} = 1;
			my $color2 = $green;
			my $color = $green;
			my $colors = [];
			
			my $span =$exon->getGenomicSpan->intersection($spanTr);#->intersection($exon->getGenomicSpan);
		unless ($span->is_empty){
			my $gc =  $no->get($patient->name,$transcript->id."_raw");
			#warn Dumper ($gc);
			$color = $purple;
			$color2=$purple;
			 $data->{$exon->id}->{$patient->id}->{green}=0;
			 my $cover = $gc->coverage($exon->start,$exon->end);
			 $data->{$exon->id}->{$patient->id}->{raw_min} = $cover->{min};
			 $data->{$exon->id}->{$patient->id}->{raw_mean} = $cover->{mean};
		}
		elsif ($exon->isExon() && $exon->is_noncoding && $utr ne 1){
			$color = $grey;
			$color2=$grey;
			 $data->{$exon->id}->{$patient->id}->{green}=0;
		}
		elsif ($covered == 0){
			$color = $grey;
			$color2=$grey;
			 $data->{$exon->id}->{$patient->id}->{green}=0;
		}
		elsif ($min <$limit){
			$alert =1;
			$color = $red;
			$color2= $yellow;
			$data->{$exon->id}->{$patient->id}->{red} = 1;
			 $data->{$exon->id}->{$patient->id}->{green}=0;
			
		}
		elsif ($min>15){
			$color = $green;
		}
		else {
			my $ar = $gradient->($min);
			my $id_color = join(";",@$ar);
			$hcolors->{$id_color} = $image->colorAllocate(@$ar) unless exists $hcolors->{$id_color}  ;
			#warn $patient->id." ar  ".join(";",@$ar);
			$color = $hcolors->{$id_color} ;
			#	warn $patient->id." color ".join(";",$image->rgb($color));
		}
		
#		elsif ($min>$max1){
#			$color = $green2;
#		}
#		elsif ($min>$max2){
#			$color=$green3;
#		}
#		elsif ($min>$max3){
#			$color=$green3;
#		}
#		
#		else {
#			my $red1 = int(($xr + (($val * ($yr-$xr)) / $max) ));
#			my $green1 = int(($xg + (($val * ($yg-$xg)) / $max) ));
#			my $blue1 = 50;#int(($xb + (( $val* ($yb-$xb)) / $max) ));
#			$color =  $image->colorAllocate($red1,$green1,$blue1);	
#			
#			$color = $red2;
#		}

	#	warn $patient->id." ".join(";",$image->rgb($color));
		$data->{$exon->id}->{$patient->id}->{color} = [$image->rgb($color)];
	
		$data->{$exon->id}->{$patient->id}->{color2} = [$image->rgb($color2)];
		#	die "popo ". $exon->id." ".Dumper $data->{$exon->id}->{$patient->id}->{color} if $data->{$exon->id}->{$patient->id}->{color}->[0] == 256 ;
		$image->filledRectangle($x,$y,$x+$size,$y+$size,$color);
		#$image->filledEllipse($x-1,$y-1,$size+1,$size+1,$color);
	#	$image->rectangle($x+1,$y+1,$x+$size+1,$y+$size+1,$black);
		$x+=$size+2
		} #end patient
		$y+=$size+2;
	
	}#end exon
$image->rectangle(0,0,$w,$h,$black);
binmode STDOUT;

#my $z = $image->png;
my $re;
$re->{alert} = 1 if $alert;
$re->{image} = $image;
$re->{data} = $data;
return $re;


}

sub horizontal_image_cnv {
	 my ($patient,$transcript,$cgi) = @_;
	 my $primers = $transcript->project->getPrimersByObjects($transcript);
	 my $size = 3;
	 my $x = 0;
	my $y = 0;
	my $xr = 255;
	my $xg = 10;
	my $xb = 0;
	my $yr = 0;
	my $yg = 255;
	my $yb = 0;
	my $nb_col = scalar(@$primers);
	 my $nb_line =1;
	 my $w =$nb_col*($size+2)+2;
	my $h = ($size+2)*$nb_line+2;
	
	
}

sub image_cache_cnv {
	 my ($patients,$transcript,$primers) = @_;
 $primers = $transcript->project->getPrimersByObjects($transcript) unless $primers;

my $hrun_id;
foreach my $p (@{$patients}) {
	$hrun_id->{$p->getRun->id} ++;
} 
my $nb_line =15;
my $size = 3;

my $x = 0;
my $y = 0;

my $xr = 255;
my $xg = 10;
my $xb = 0;

my $yr = 0;
my $yg = 255;
my $yb = 0;

my $nb_col = scalar(@$patients);
my $w =$nb_col*($size+2)+2;
my $h = ($size+2)*scalar(@$primers)+2;

my $gradient = gradient( 1, 50);

my $image= new GD::Image($w+1,$h+1);
		

my        $grey1 = $image->colorAllocate(90,90,90);	
my        $grey2 = $image->colorAllocate(140,140,140);
my        $grey3 = $image->colorAllocate(0, 255, 0);
my        $blue = $image->colorAllocate(0,0,255);
my        $red = $image->colorAllocate(255,10,50);   
my       $white = $image->colorAllocate(255,255,255);
my       $black = $image->colorAllocate(0,0,0);
my $hcolors;
my $zz = ($size+2)*scalar(@$primers)+1;
my $data={};
$y =2;

my $hdata;
	foreach my $primer ( @$primers) {
		 my $runs = $primer->getRuns();
		 
			foreach my $patient (@{$patients}){
				push(@{$hdata->{$patient->id}->{score}},$primer->cnv_score($patient)) if $primer->level($patient) > -1;
				push(@{$hdata->{$patient->id}->{score}},-1) if $primer->level($patient) ==  -1;
				$hdata->{$patient->id}->{level_string} .= $primer->level($patient);
				$hdata->{$patient->id}->{level} ++ if $primer->level($patient) == 1;
				$hdata->{$patient->id}->{level} ++ if $primer->level($patient) == 2;
				
			}
	}
	
		my	$max_event = 0;
		
	foreach my $patient (@$patients){
		my $string = $hdata->{$patient->id}->{level_string};
			my $count = ($string =~ tr/1//);
			$max_event = $count if $count > $max_event;
			$count = ($string =~ tr/2//);
			$max_event = $count if $count > $max_event;
		$data->{$patient->id}->{score_smooth} =  smooth_data($hdata->{$patient->id}->{score});
	}
	
	my @color_background  = (223, 255, 216);
 @color_background =(127, 255, 27) if ($max_event >=2);
 my @array_ho;
for  (my $i = 0 ; $i<@$primers;$i++ ){
	my $primer = $primers->[$i];
	my (@find) = grep{exists $hrun_id->{$_->id}} @{$primer->getRuns};
	next unless (@find);
	#foreach my $primer ( @$primers) {
			$x =2;
			foreach my $patient (sort{$a->name cmp $b->name} @$patients){
				my @colors =();
			my $level = $primer->level($patient) ;
			my $score = $primer->cnv_score($patient);
			my $vscore = $score; 
				my $rscore = $score;
				  $rscore =  $data->{$patient->id}->{score_smooth}->[$i] ;# if exists $hdata->{$patient->id}->{$primer->id};
				 $score = $rscore;
			if ($level == 0 && exists $hdata->{$patient->id}->{level}) {
				my $sd = $primer->sd($patient);
				
			
			
				 if ($rscore - $sd/2 >= 1.3){
					$level = 3;
				}	
				if ( $rscore + $sd/2   <= 0.6 ){
					$level = 4;
				}	
			}
			my $name = $patient->name;
			@colors = (223, 255, 216);
			@colors =@color_background;
			my $colors;
			if ($level == -2 ){
				@colors=(90,90,90);	
				#$color = $grey1;
			}
			elsif ($level == -1) {
				@colors=(140,140,140);
				}
				elsif ($level == 2){
						@colors =(52,152,255);
						@colors =(2, 71, 253);
				}
				elsif ($level==3){
					my $colorr =int (150*($score-1))+105;
					@colors =(52,152,255);
					@colors =(2, 71, 253);
				}
				elsif ($level == 1){
						@colors =(229, 40, 76);
						@colors =(255, 64, 64);
						@colors =(255, 0, 0) if $vscore < 0.15;
				}
				elsif ($level==4){
					my $colorr =int (50*(0.6-$score))+205;
					@colors =(229, 40, 76);
					
				}
				elsif ($level ==0){
					#FFC9CF
					if ($score <=0.1){
						@colors = (0, 0,0)
					}
					elsif ($score <=0.55){
						@colors = (254, 138, 53)
					}
						elsif ($score <=0.6){
						@colors =(246, 179, 134);
					}
					elsif ($score <=0.7){
						@colors =(253, 218, 178);
					}
					elsif ($score >= 1.5){
						@colors = (115, 178, 217)
					}
					elsif ($score >= 1.4){
					@colors =(205, 223, 245);
				}
					
					
				}
				my $id_color = join(",",@colors);
				unless (exists $hcolors->{$id_color}){
					$hcolors->{$id_color} = $image->colorAllocate(@colors);
				}
				
				

				if ( $vscore < 0.2 ) {
					my $x1 = $x;
					my $x2 = $x + ( $size + 1 );
					my $y1 = $y;
					my $y2 = $y + $size;
					push( @array_ho, [ $x1 - 1, $y1 - 1, $x2, $y2 ] );

				}
				
			
				$image->rectangle($x-1,$y-1,$x+$nb_col*($size+2)+1,$y+$size+1,$black);
				$image->filledRectangle($x,$y,$x+$nb_col*($size+2),$y+$size,$hcolors->{$id_color} );
				$data->{$primer->id}->{$name}->{colors} = \@colors;
				$data->{$primer->id}->{$name}->{text} = $level;
				$data->{$primer->id}->{$name}->{cnv_score} = $primer->cnv_score($patient)."(".$rscore."):> ".$level." z:". $primer->zscore($patient)." sd:".$primer->sd($patient);;
				$data->{$primer->id}->{$name}->{cnv_score} = $primer->cnv_score($patient)  if $level == -1;
				$x+=$size+2
			}
		$y+=$size+2;
	}
	my $point    = $image->colorAllocate( 255, 255, 251 );
	foreach my $ho (@array_ho) {
		$image->rectangle( @$ho, $point );

		#last;
	}
	#next unless $exon->{exon} ==1;
$image->rectangle(0,0,$w,$h,$black);

my $re;
$re->{image} = $image;
$re->{data} = $data;
return $re;


}
sub find_previous {
	my ($pos,$data) = @_;
	return ($pos,undef) if $pos-1 <= 0;
	warn $pos unless $data->[$pos-1] ;
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
	return $data unless $data;
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
		  push(@data2,$z/(scalar(@trio)));
		 
	}
	return \@data2;
}

1;

