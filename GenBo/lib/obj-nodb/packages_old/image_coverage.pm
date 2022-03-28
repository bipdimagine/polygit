package image_coverage;
use strict;
use Moose;
use MooseX::Method::Signatures;
use GD;
use Data::Printer;
use Data::Dumper;
 use List::Util qw( min sum);
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
	 
	
sub image {
 my ($patients,$transcript,$intronic,$utr,$padding,$limit,$all) = @_;
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
my        $purple = $image->colorAllocate(214,119,253);
my @warning ;
		  $warning[0] = $image->colorAllocate(0, 117, 243);	
		 # $warning[1] = $image->colorAllocate(255, 0, 0 );	
		  $warning[1] = $image->colorAllocate(0, 117, 243);	
my        $grey2 = $image->colorAllocate(250,250,250);


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
			$covered = -1 if $s1->is_empty;
		}
		else {
			$covered = 0 if scalar($s1->as_array)<20;
		}
	unless ($exon->isExon){
		$image->filledRectangle($x,$y,$x+$nb_col*($size+2),$y+$size,$grey2);
	}
	#next unless $exon->{exon} ==1;
  	my $spanTr = $no->get($patients->[0]->name,$transcript->id."_spandup");	
	foreach my $patient (sort{$a->name cmp $b->name} @$patients){
	
	#	warn $spanTr;
		$spanTr = Set::IntSpan::Fast::XS->new() unless $spanTr; 
		my $dir = $patient->project->getVariationsDir("duplicate_region_calling")."/regions";
		my $filebed =  $patient->project->getVariationsDir("duplicate_region_calling")."/regions/".$patient->name().".dup.bed";
		#system("mkdir $dir && chmod a+rwx $dir" ) unless -e $dir;
	

		my $pid =$patient->id;
		my ($mean,$intspan,$min);
		#my $t  =$exon->cached_statistic_coverage_coding(patient=>$patient,padding=>$padding,limit=>$limit,utr=>$utr);  
		my $t = $no->get($patient->name,$exon->id."_".$padding."_".$utr);
		if ($t){
			$mean = $t->{mean};
			$min = $t->{min};
			#die($mean);
		}
#		elsif ($patient->isNoSqlDepth()) {
#			my $xa =  $patient->depth($exon->getChromosome->name,$exon->start-$padding,$exon->end+$padding);
#			$min = min(@$xa);
#			$mean = sum(@$xa)/scalar(@$xa);
#		}
		else {
		 	confess($no->dir);
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
		elsif ($covered == -1){
			my $c = $warning[rand(2)];
			$color = $c;
			$color2=$c;
			 $data->{$exon->id}->{$patient->id}->{green}=0;
		}
		elsif ($covered == 0){
			$color = $grey;
			$color2=$grey;
			 $data->{$exon->id}->{$patient->id}->{green}=0;
		}
		elsif ($min <$limit){
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
$re->{image} = $image;
$re->{data} = $data;
return $re;


}

1;

