package decode_prediction_matrix;
use strict;
use List::Util qw(max);
use POSIX qw(ceil);
use Data::Dumper;
use Compress::Zlib;
use Carp;

my $PREDICTION_TO_VAL = {
    polyphen => {
        'probably damaging' => 0,
        'possibly damaging' => 1,
        'benign'            => 2,
        'unknown'           => 3,
    },

    sift => {
        'tolerated'     => 0,
        'deleterious'   => 1,
    },
};

our $HEADER     = 'VEP';

# the format we use when 'pack'ing our predictions

my $PACK_FORMAT = 'v';

our @ALL_AAS = qw(A C D E F G H I K L M N P Q R S T V W Y);

our $NO_PREDICTION = pack($PACK_FORMAT, 0xFFFF);
my $BYTES_PER_PREDICTION = 2;

# the maximum number of distinct qualitative predictions used by any tool

my $MAX_NUM_PREDS = max( map { scalar keys %$_ } values %$PREDICTION_TO_VAL ); 

# the number of bits used to encode the qualitative prediction

my $NUM_PRED_BITS = ceil( log($MAX_NUM_PREDS) / log(2) );

throw("Cannot represent more than ".(2**6-1)." predictions") if $NUM_PRED_BITS > 6;

# a hash mapping back from a numerical value to a qualitative prediction

my $VAL_TO_PREDICTION = {
    map {
        my $tool = $_; 
        $tool => {
            map {
                $PREDICTION_TO_VAL->{$tool}->{$_} => $_ 
            } keys %{ $PREDICTION_TO_VAL->{$tool} }
        }
    } keys %$PREDICTION_TO_VAL
};

# a hash from amino acid single letter code to a numerical value

our $AA_LOOKUP = { map {$ALL_AAS[$_] => $_} 0 .. $#ALL_AAS };

# the number of valid amino acids

our $NUM_AAS = scalar(@ALL_AAS);


sub prediction_from_matrix {
	 my ($matrix, $pos, $aa) = @_;
	 $aa = uc($aa) if defined $aa;
	  confess("Invalid position: $pos") unless (defined $pos && $pos > 0);
    
   return(0,-1)  unless (defined $aa && defined $AA_LOOKUP->{$aa});
    my $offset = compute_offset($matrix,$pos, $aa);
	
    
    if ($offset + 1 > length($matrix)) {
       # warn "Offset outside of prediction matrix for position $pos and amino acid $aa?";
       return(0,-1);
        return undef;
    }
    
    my $pred = substr($matrix, $offset, $BYTES_PER_PREDICTION);
    return prediction_from_short($pred);
}

sub compute_offset{
	my($matrix,$pos,$aa) = @_;
	
    my $offset = length($HEADER) + ( ( (($pos-1) * $NUM_AAS) + $AA_LOOKUP->{$aa} ) * $BYTES_PER_PREDICTION );

    return $offset;
}

sub prediction_from_short {
    my ($val) = @_;

    # check it isn't our special null prediction

    if ($val eq $NO_PREDICTION) {
        return (0,-1);
    }

    # unpack the value as a short

    $val = unpack $PACK_FORMAT, $val;

    # shift the prediction bits down and look up the prediction string
	my $prediction = $val >> (16 - $NUM_PRED_BITS);
	
   # my $pred = $VAL_TO_PREDICTION->{$self->{analysis}}->{$val >> (16 - $NUM_PRED_BITS)};

    # mask off the top 6 bits reserved for the prediction and convert back to a 3 d.p. float

    my $prob = ($val & (2**10 - 1)) / 1000;


    return ($prediction, $prob);
}

1;