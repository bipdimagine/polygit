package NeedlemanWunsch;

# perl/NeedlemanWunsch.pm
# url: https://github.com/noporpoise/seq-align
# maintainer: Isaac Turner <turner.isaac@gmail.com>
# license: Public Domain, no warranty
# date: Oct 2012

use strict;
use warnings;

use File::Basename; #dirname
use Carp; # for reporting warnings and errors
use FileHandle; # provides autoflush
use IPC::Open2; # for opening handle to a process
use List::Util qw(max);

sub new
{
  my ($class, $buffer,@args) = @_;

  my %options = (@args);
  my $exe = $buffer->software("needleman_wunsch");
  my $cmdline = $exe;
  #my $cmdline = dirname(__FILE__)."/../bin/needleman_wunsch";
  my ($gapopen, $gapextend) = (-2, -1);
  my $timeout = 5;

  for my $key (keys %options)
  {
    $options{lc($key)} = $options{$key};
  }

  if(defined($options{'match'}) != defined($options{'mismatch'}))
  {
    carp("Cannot set only one of match/mismatch");
  }

  if(defined($options{'cmd'}))
  {
    $cmdline = $options{'cmd'};
  }

  $cmdline .= " --stdin --pretty --printscores";

  while(my ($key,$value) = each(%options))
  {
    if(grep(/^$key$/i, qw(case_sensitive nogaps nomismatches)))
    {
      # 'Flag' options -- they have no args
      if($value)
      {
        $cmdline .= " --$key";
      }
    }
    elsif($key =~ /^scoring$/i)
    {
      if(!grep($value, qw(PAM30 PAM70 BLOSUM80 BLOSUM62)))
      {
        carp("Invalid NeedlemanWunsch.pm scoring '$value'");
      }
      else
      {
        $cmdline .= " --scoring $value";
      }
    }
    elsif($key eq "timeout") {
      $timeout = $value;
    }
    elsif(grep(/^$key$/i, qw(substitution_matrix substitution_pairs
                             match mismatch gapopen gapextend wildcard)))
    {
      if(lc($key) eq "gapopen") { $gapopen = $value; }
      if(lc($key) eq "gapextend") { $gapextend = $value; }

      $cmdline .= " --$key $value";
    }
    else
    {
      carp("Unknown option '$key' => '$value' (ignored)");
    }
  }

  $cmdline .= " 2>&1";

  my ($in, $out, $err);

 # print "running '$cmdline'\n";

  my $pid = open2($in, $out, $cmdline)
    or die("Cannot run cmd: '$cmdline'");

  $in->autoflush();
  $out->autoflush();

  my $self = {
    _in => $in,
    _out => $out,
    _err => $err,
    _pid => $pid,
    _gapopen => $gapopen,
    _gapextend => $gapextend,
    _align_number => -1,
    _timeout => $timeout,
    _seq1 => undef,
    _seq2 => undef,
    _cmd => $cmdline
  };
  
  bless($self, $class);

  return $self;
}

sub destructor
{
  my ($self) = @_;
 close($self->{_in});
  close($self->{_out});

  waitpid($self->{_pid}, 1);
  
}
 sub DESTROY {
      my $self = shift;
     	close($self->{_in});
  		close($self->{_out});
  		waitpid($self->{_pid}, 1);
  }
sub read_line
{
  my ($self) = @_;

  # print "NW Waiting...\n";

  my $in = $self->{_in};
  my $next_line;

  # Reading with time out
  # http://www.perlmonks.org/?node_id=43304
  eval {
    local $SIG{ALRM} = sub { die "timeout\n" };
    alarm($self->{_timeout});
    $next_line = <$in>;
  };

  if($@ eq "timeout\n") { die("Error: timeout reading NW output"); }
  elsif($@) { die("Error: couldn't read output"); }
  alarm(0);
  
  if(defined($next_line))
  {
    chomp($next_line);
    # print "IN: '$next_line'\n";

    if($next_line =~ /^Error:/i)
    {
      croak($next_line);
    }
  }

  return $next_line;
}

sub do_alignment
{
  my ($self, $seq1, $seq2) = @_;

  my %result = ('seq1' => $seq1, 'seq2' => $seq2,
                'number' => $self->{_align_number}++);

  if($seq1 =~ /[\n\r]/ || $seq2 =~ /[\n\r]/)
  {
    croak("New lines not allowed in sequences");
  }
  elsif($seq1 eq '' || $seq2 eq '')
  {
    $result{'align1'} = $seq1;
    $result{'align2'} = $seq2;
    my $len = max(length($seq1), length($seq2));
    $result{'sep'} = '-' x $len;
    $result{'score'} = ($len>0 ? $self->{_gapopen}+$len*$self->{_gapextend} : 0);
    return \%result;
  }

  my $out = $self->{_out};

  print $out "$seq1\n$seq2\n";

  $result{'align1'} = $self->read_line();
  $result{'sep'} = $self->read_line();
  $result{'align2'} = $self->read_line();
  my $score_line = $self->read_line();

  if(!defined($result{'align1'}) || !defined($result{'sep'}) ||
     !defined($result{'align2'}) || !defined($score_line))
  {
    die("Missing lines when reading in - have you compiled?");
  }

  if($score_line =~ /(-?\d+)$/i)
  {
    $result{'score'} = $1;
  }
  else
  {
    croak("Cannot locate score in string '".$score_line."'");
  }

  # Skip remaining line
  my $skip = $self->read_line();
  chomp($skip);

  if($skip ne '')
  {
    croak("Not expecting line: '$skip'");
  }
	$self->{results} = \%result;
  return \%result;
}
sub get_next_hit
{
  my ($self) = @_;
  
  return $self->{results};
}

sub print_alignment
{
  my ($self, $hit, $out) = @_;

  if(!defined($out))
  {
    open($out, ">-");
  }

  print $out $hit->{'align1'}."\n";
  print $out $hit->{'sep'}."\n";
  print $out $hit->{'align2'}."\n";
  print $out "score: $hit->{'score'}\n\n";
}

sub test {
	my ($self,$seq1,$seq2) = @_;
	my $MATCH    =  1; # +1 for letters that match
my $MISMATCH = -1; # -1 for letters that mismatch
my $GAP      = -1; # -1 for any gap

# initialization
my @matrix;
$matrix[0][0]{score}   = 0;
$matrix[0][0]{pointer} = "none";
for (my $j = 1; $j <= length($seq1); $j++) {
    $matrix[0][$j]{score}   = $GAP * $j;
    $matrix[0][$j]{pointer} = "left";
}
for (my $i = 1; $i <= length($seq2); $i++) {
    $matrix[$i][0]{score}   = $GAP * $i;
    $matrix[$i][0]{pointer} = "up";
}

# fill
for (my $i = 1; $i <= length($seq2); $i++) {
    for (my $j = 1; $j <= length($seq1); $j++) {
        my ($diagonal_score, $left_score, $up_score);

        # calculate match score
        my $letter1 = substr($seq1, $j-1, 1);
        my $letter2 = substr($seq2, $i-1, 1);                            
        if ($letter1 eq $letter2) {
            $diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
        }
        else {
            $diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
        }

        # calculate gap scores
        $up_score   = $matrix[$i-1][$j]{score} + $GAP;
        $left_score = $matrix[$i][$j-1]{score} + $GAP;

        # choose best score
        if ($diagonal_score >= $up_score) {
            if ($diagonal_score >= $left_score) {
                $matrix[$i][$j]{score}   = $diagonal_score;
                $matrix[$i][$j]{pointer} = "diagonal";
            }
	    else {
                $matrix[$i][$j]{score}   = $left_score;
                $matrix[$i][$j]{pointer} = "left";
            }
        } 
        else {
            if ($up_score >= $left_score) {
                $matrix[$i][$j]{score}   = $up_score;
                $matrix[$i][$j]{pointer} = "up";
            }
            else {
                $matrix[$i][$j]{score}   = $left_score;
                $matrix[$i][$j]{pointer} = "left";
            }
        }
    }
}

# trace-back

my $align1 = "";
my $align2 = "";

# start at last cell of matrix
my $j = length($seq1);
my $i = length($seq2);

while (1) {
    last if $matrix[$i][$j]{pointer} eq "none"; # ends at first cell of matrix

    if ($matrix[$i][$j]{pointer} eq "diagonal") {
        $align1 .= substr($seq1, $j-1, 1);
        $align2 .= substr($seq2, $i-1, 1);
        $i--;
        $j--;
    }
    elsif ($matrix[$i][$j]{pointer} eq "left") {
        $align1 .= substr($seq1, $j-1, 1);
        $align2 .= "-";
        $j--;
    }
    elsif ($matrix[$i][$j]{pointer} eq "up") {
        $align1 .= "-";
        $align2 .= substr($seq2, $i-1, 1);
        $i--;
    }    
}

$align1 = reverse $align1;
$align2 = reverse $align2;
return ($align1,$align2);
	
}


1;
