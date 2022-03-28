#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use Data::Dumper;
use threads;
use Thread::Queue;
use Getopt::Long;
use File::Temp;
use Set::IntSpan::Fast::XS;
my $bam_file;
my $file_primer;

GetOptions(
	'file=s'  	=> \$bam_file,
	'primer=s' => \$file_primer,
);
die() unless -e $file_primer;
my $chr_intspan;
my $chr_allintspan;
$| =1;
open (PRIMER,"$file_primer");
while(my $line = <PRIMER>){
	chomp($line);
	my @data = split(" ",$line);
	#warn join("-",@data);
	unless (exists $chr_intspan->{$data[0] } ) {
		$chr_intspan->{$data[0]}->{start}= Set::IntSpan::Fast::XS->new();
		$chr_intspan->{$data[0]}->{end}= Set::IntSpan::Fast::XS->new();
	}
	$chr_intspan->{$data[0]}->{end}->add_range($data[3],$data[4]+10);
	$chr_intspan->{$data[0]}->{start}->add_range($data[1]-10,$data[2]);
	
}

my $estart;
my $eend;
my $tot;
open (TOTO , "samtools view -h $bam_file|");
my $nb;
my $problem;
while (my $bam = <TOTO> ){
	chomp($bam);
	
	if ($bam =~ /^@/){
		print $bam."\n";
		next;
	}
	
	
	my @array = split(" ",$bam);
	my $chr = $array[2];
	my $flag = $array[1];
	my $start = $array[3];
	my $seq = $array[9];
	my $ll1 = length($seq);
	my $ll = length($seq);
	my $newcigar =$array[5];
	my $cigar = parse_cigar($array[5]);
	 if ($newcigar eq "*" || $newcigar =~ /H/){
	 	print join("\t",@array)."\n";
	 	next;
	 	 
	 }
	 
	next if $newcigar eq "*";
	my $s = 0;
	$s = $cigar->[0]->{len} if $cigar->[0]->{code} eq 'S';
	my $s2 = 0;
	$s2 = $cigar->[-1]->{len} if $cigar->[-1]->{code} eq 'S';
	my $first_position = $start;
	my $align_len = $ll -$s -$s2 -1;
	my $last_position = $first_position + $align_len ;
	
	#my $intspan = = Set::IntSpan::Fast::XS->new();
	#$intspan->add_range($first_position,$last_position);
	 unless (exists $chr_intspan->{$chr}){
	 	print join("\t",@array)."\n";
	 	next;
	 }
	$tot++;
	#warn $chr;
	my @newseq = split("",$seq); 
	unless ($flag & 16){
		#forward
		if ($chr_intspan->{$chr}->{start}->contains($first_position)){
			my $z = $s;
		my $N = 0;
			for (my $i=$first_position;$i<$last_position;$i++){
				last unless ($chr_intspan->{$chr}->{start}->contains($i));
				
				
				$newseq[$z]= "N";
				$N++;
				$z++;
			}
			
			
			my $line = cigar_to_line($cigar);
			my @ct = split("",$line);
			my $pos = 0;
			my $ii =0;
			my @ct2;
			my $NN = $N;
			foreach my $c (@ct){
				next if $pos>= $N;
				if ($c eq "S"){
					$NN --;
				}
				if ($c eq "I"){
					$NN --;
				}
				if ($c eq "D"){
					$ct[$pos] = "" ;
					$N++;
					$NN++;
				}
				else {
					$ct[$pos] = "S";
				}
			
				$pos++;
			
			}
			$array[3] += $NN; 
			#my $current_code=$ct[0];
		
	
			$newcigar =rewrite_cigar(\@ct);
			#die() if (scalar(@newseq)== 247);
			for (my $i=0;$i<$s;$i++){
				$newseq[$i]= "X";
			}
			$estart++;
		}
	}
	else {
	
	if ($chr_intspan->{$chr}->{end}->contains($last_position)){
		my $z = $ll - $s2 -1;
		my $N;
		for (my $i=$last_position;$i>$start;$i--){
			last unless $chr_intspan->{$chr}->{end}->contains($i);
			$newseq[$z]= "N";
			$N++;
			$z --;
		}
		for (my $i=1;$i<=$s2;$i++){
			my $x = $i*-1;
			$newseq[$x]= "N";
			$N++;
		} 
	

	my $line = cigar_to_line($cigar);
	my @ct = split("",$line);
	my $pos = scalar(@ct);
	my $debug;
	$debug =1 if (scalar(@newseq) == 82);
	warn $newcigar if $debug ==1;
		for (my $i = 0;$i<$N;$i++){
			$pos --;
			if ($ct[$pos] eq "D"){
				$ct[$pos] = "";
				$N++; 
			}
			else {
			$ct[$pos] = "S"; 
			}
		}
		$newcigar =rewrite_cigar(\@ct);
		warn $newcigar if $debug ==1;
		#die() if $debug =1;
		$eend++;
	}
	}
	
	
	$array[9] = join("",@newseq);
	$array[5] = $newcigar;
	
	print join("\t",@array)."\n";
}
close TOTO;
warn $estart." ".$eend." ".$tot;
exit(0);

sub threePrime{
		my ($array,$pos1,$seq2,$seq1) = @_;
	 	my $cigar = parse_cigar($array->[5]);
	 	 my $start = length($seq1)+$pos1;
	 	 
	  if ($cigar->[-1]->{code} eq "S" ){
	  	
	  		$cigar->[-1]->{len} +=( (length($seq2) - $start) -1);
	  		
	  		 $array->[5] = cigar_string($cigar);
			return;	 	
	 }
		
		 
		 my $current = 1;
		 my $new_cigar;
		 my $debug;
		 foreach my $c (@$cigar){
		 	my $p = $c->{len};
		 	if ($c->{code} eq "D") {
		 		push(@$new_cigar,$c);
		 		next;
		 	}
		 	if ($current + $p < $start){
		 		$current += $p;
		 		warn "suite ".$c->{code} if $debug;
		 		push(@$new_cigar,$c);
		 		next;
		 	}
		 	
		 	unless ($c->{code} eq "S"){
		 	warn $c->{len}." ".$start." ".$current if $debug;
		  	$c->{len} = $start - $current ;
		  	push(@$new_cigar,$c);
		  	$current += $c->{len};
		 	}
		  	my $cc;
		  	$cc->{len} = ((length($seq2)) - $current)+1;
		  	$cc->{code} = "S";
		  	push(@$new_cigar,$cc);
			 
			last;
		 }
		
		# warn "3P ". cigar_string($new_cigar)." ".$array->[5]." ".length($seq2);	
		
		 $array->[5] = cigar_string($new_cigar);
	
}

sub fivePrime {
	my ($array,$pos1,$seq2) = @_;
	my $current =1;
	 my $cigar = parse_cigar($array->[5]);
	 if ($cigar->[0]->{code} eq "S" && $cigar->[0]->{len} > 10 ){
	 	return;
	 }
	my $current_pos = $pos1;
	my $new_cigar;
	
	push(@$new_cigar,{code=>"S",len=>$pos1});
	my $find;
	my $gpos = $array->[3];
#	 foreach my $c (@$new_cigar){
#		 	last if $c->{code} ne 'S' ||  $c->{code} ne 'D';
#		 	$gpos -= $c->{len};
#		}
	if ($cigar->[0]->{code} eq "S" ){
	 	$gpos -= $cigar->[0]->{len};
	 	
	 }
	 $gpos += $pos1;
	
	foreach my $c (@$cigar){
		if ($find){
			push(@$new_cigar,$c);
			next;
		}
		
		if ($c->{code} eq "D") {
			$gpos += $c->{len};
		 		#push(@$new_cigar,$c);
		 		next;
		 }
		if ($current_pos - $c->{len}  >0){			
			$current_pos -= $c->{len};
			next;
		};
		
		unless ($c->{code} eq "S"){
			
			$c->{len} = $c->{len} - $current_pos;
			push(@$new_cigar,$c) if $c->{len} > 0;
		}
		else {
			$new_cigar->[0]->{len} = $c->{len};
		}

		$find =1;

	}
	#warn "5P ". cigar_string($new_cigar)." ".$array->[5]." ".length($seq2);	
	#$array->[3] = $gpos;	 
	$array->[5] = cigar_string($new_cigar);
	
}

sub rewrite_cigar {
	my ($ct) = @_;
				my $current_code=$ct->[0];
			my $newcigar ="";
			my $len_c =0;
			foreach my $c (@$ct){
				next if $c eq "";
				if ($current_code ne $c){
						$newcigar .=$len_c.$current_code;
						$current_code = $c;
						$len_c = 0;
					}
					$len_c  ++;
			}
				$newcigar .=$len_c.$current_code;
				return $newcigar;
	
}
sub cigar_string {
	my ($cigar) = @_;
	my $st;
	foreach my $c (@$cigar){
		$st .= $c->{len}.$c->{code};
	}
	return $st;
}
sub cigar_to_line {
	my ($cigar) = @_;
	my $line="";
		foreach my $c  (@$cigar){
				for (my $i=0;$i< $c->{len};$i++){
					$line .=$c->{code};
				}
			}
	return $line;
}
sub parse_cigar {
	my ($cigar) = @_;
	my $res;
	while ( $cigar =~ m/(\d+)([A-Z])/g){
		
		my $c;
		 $c->{len} = $1;
		$c->{code} = $2;
		push(@$res,$c);
		
	}
	return $res;
	
}
