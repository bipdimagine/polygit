package image_primer_coverage;
use strict;
use Moose;
use MooseX::Method::Signatures;
use GD;
use Data::Printer;
use Data::Dumper;
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
sub image {
 my ($patients,$capture,$multi,$limit) = @_;
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



my @primers =   sort{$b->getChromosome->length <=> $a->getChromosome->length || $a->start <=> $b->start } @{$capture->getPrimersByMultiplex($multi)};




my $nb_col = scalar(@$patients);
my $w =$nb_col*($size+2)+2;
my $h = ($size+2)*scalar(@primers)+2;

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
my $hcolors;
#$image->rectangle(1,1,$nb_col*($size+2)+1,($size+2)*scalar(@exons)+1,$black);
	my $data;
	$y =2;
	my $dir = $capture->project->getCacheDir();
	my $no =  $capture->project->noSqlCoverage();
	

	foreach my $primer (@primers) {
	$x =2;
	my $covered =1;
	#next unless $exon->{exon} ==1;
	foreach my $patient (sort{$a->name cmp $b->name} @$patients){
		my $pid =$patient->id;
		 my ($mean,$intspan,$min) ;
		 
		  my $hash = $primer->cached_statistic_coverage($patient);
		  if ($hash){
		 	$mean = $hash->{mean};
		 	$min = $hash->{min};
		  }
		  else {
		  	confess();
		  }
			my $val = $min>$max ? $max : $min;
			$val -= $limit;
			#$val = 100 - $val;
			$data->{$primer->id}->{$patient->id}->{mean} = $mean;
			$data->{$primer->id}->{$patient->id}->{min} = $min;
			$data->{$primer->id}->{$patient->id}->{covered} = $covered;
			$data->{$primer->id}->{$patient->id}->{red} = 0;
			$data->{$primer->id}->{$patient->id}->{green} = 1;
			my $color2 = $green;
			my $color = $green;
			my $colors = [];

			$color = $red;
			$color2= $yellow;
			$data->{$primer->id}->{$patient->id}->{red} = 1;
			 $data->{$primer->id}->{$patient->id}->{green}=0;
		if ($min< $limit){
			$color = $red;
		}
		elsif ($min>50 && $min>$limit){
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
		

		$data->{$primer->id}->{$patient->id}->{color} = [$image->rgb($color)];
	
		$data->{$primer->id}->{$patient->id}->{color2} = [$image->rgb($color2)];
		$image->filledRectangle($x,$y,$x+$size,$y+$size,$color);
		#$image->filledEllipse($x-1,$y-1,$size+1,$size+1,$color);
	#	$image->rectangle($x+1,$y+1,$x+$size+1,$y+$size+1,$black);
		$x+=$size+2
		} #end patient
		$y+=$size+2;
	
	}#end exon
	#$image->rectangle(0,0,$w-1,$h-1,$black);
$image->rectangle(0,0,$w,$h,$black);
#$image->filledRectangle($w+1,$h+1,$w-2,$h+2,$black);
binmode STDOUT;

#my $z = $image->png;
my $re;
$re->{image} = $image;
$re->{data} = $data;
return $re;


}

1;

