#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use Data::Dumper;
use threads;
use Thread::Queue;
use Getopt::Long;
use File::Temp;
use Set::IntSpan::Fast::XS;
use Parallel::ForkManager;

my $bam_file;
my $file_primer;
my $fork;
 use Bio::DB::Sam;
 
GetOptions(
	'file=s'  	=> \$bam_file,
	'primer=s' => \$file_primer,
	'fork=s' => \$fork,
);
$| =1;
my $chr_intspan;
my $chr_allintspan;
my $couple;
open (PRIMER,"$file_primer");
while(my $line = <PRIMER>){
	chomp($line);
	my @data = split(" ",$line);
	#warn join("-",@data);
	my $primers;
	unless (exists $chr_intspan->{$data[0] } ) {
		$chr_intspan->{$data[0]}->{start}= Set::IntSpan::Fast::XS->new();
		$chr_intspan->{$data[0]}->{end}= Set::IntSpan::Fast::XS->new();
		
	}
	  $primers->{start} = Set::IntSpan::Fast::XS->new();
	$primers->{start}->add_range($data[1]-10,$data[2]);
	$chr_intspan->{$data[0]}->{start} = $chr_intspan->{$data[0]}->{start}->union($primers->{start});
    

	 $primers->{end} = Set::IntSpan::Fast::XS->new();
	 $primers->{end}->add_range($data[3],$data[4]+10);
	 $chr_intspan->{$data[0]}->{end} = $chr_intspan->{$data[0]}->{end}->union( $primers->{end});

	$primers->{all} =   $primers->{end}->union($primers->{start});

	push(@{$couple->{$data[0]}},$primers);
	
}





my $nb;
sub get_cigar_expand {
	my $cigar_array = shift;
	my $cigar_expand;
	foreach my $c (@$cigar_array){
		for  (my $i=0;$i<$c->[1];$i++){
			$cigar_expand .=$c->[0];
		}
	}
	return $cigar_expand;
}


sub transform_cigar {
	my ($cigar) = @_;
	my $res;
	my $cigar_expand;
	my $first;
	my $len =0;
	my $end =0;
	
	while ( $cigar =~ m/(\d+)([A-Z])/g){
		
		my $c;
		 $c->{len} = $1;
		$c->{code} = $2;
		for  (my $i=0;$i<$c->{len};$i++){
			$cigar_expand .=$c->{code};
		}
		if ($c->{code} ne "I" && $c->{code} ne "D" ){
			$end += $c->{len};
		}
		
		if ($c->{code} ne "D"){
			$len += $c->{len};
		}
	#$len += $c->{len}	if 	$c->{code} ne "D";
	#$len -= $c->{len}	if 	$c->{code} eq "I";
		push(@$res,$c);
		
	}
	my $dec =0;
	if ($res->[0]->{code} eq "S"){
		$dec = $res->[0]->{len};
	}
	return ($cigar_expand,$res,$dec,$len,$end);
	
}
sub restore_cigar {
	my ($cigar_line) = @_;
	
	my $type;
	my $len;
	my $cigar;
	foreach my $c (split('',$cigar_line)){
		if ($c ne $type){
			if ($type ){
				$cigar .= $len.$type;
			
			}
				$len =1;
				$type=$c;
		}
		else {
			$len ++;
		}
		
	}
		$cigar .= $len.$type;
		return $cigar;

	
}
my $test = 0;
my $test2 = 0;
my $debug_intspan = Set::IntSpan::Fast::XS->new("110678804-110678900");
my $seq_dejavu;
my $current_chr;
sub line_forward {
	my $array = shift;
	
	my $chr =  $array->[2];
		my $flag = $array->[1];
	return  unless exists $chr_intspan->{$chr};
	if ($current_chr ne $chr){
		$seq_dejavu = {};
		$current_chr = $chr;
	}
	my $hids =  $array->[3]."_".$array->[5]."_".$array->[9];
	
	if (exists $seq_dejavu->{$hids}){
		 	$array->[3] = $seq_dejavu->{$hids}->{position};
		 	print "coucou" unless $seq_dejavu->{$hids}->{position};
    			$array->[5] =  $seq_dejavu->{$hids}->{cigar};
    			print "coucou" unless $seq_dejavu->{$hids}->{cigar};
    			$array->[9] =  $seq_dejavu->{$hids}->{mask_seq};
    			print "coucou" unless $seq_dejavu->{$hids}->{mask_seq};
    			return;
    			
	}
	return if $array->[5] eq "*";
	my $debug;
#	$debug =1 if $array->[0] =~ /M01068:122:000000000-AWTF6:1:1102:16873:28910/;
	warn Dumper $array if $debug;

	
	
		return if $array->[5] =~/H/;

	my ($cigar_expand,$cigar_array1,$dec,$len,$cend) =transform_cigar($array->[5]);
	my $count = ($cigar_expand =~ tr/M//);
	return if $count < 50;
	warn $dec." :: $count ::".$cigar_expand.$cend  if $debug;

	return if $len <30;
	my $start =  $array->[3];
	my $new_start =$start -$dec;
	my $new_end = $new_start + $cend;
	my $intspan_align = Set::IntSpan::Fast::XS->new($new_start."-".$new_end);
		my $cps = $couple->{$chr};
		my $ref_pos = $new_end;
		my $pos_key = "end";
		unless ($flag & 16){
			$ref_pos = $new_start;
			$pos_key = "start";
		}
		
	
		if ($chr_intspan->{$chr}->{$pos_key}->contains($ref_pos)){
		
			my $c;
			foreach  my $cp (@$cps){
				 if ($cp->{$pos_key}->contains($ref_pos)){
				 	$c = $cp;
				 	last;
				 }
			}
			
					warn "coucou" if $debug;
			warn "problem no primer "  unless $c;
			return -1 unless $c;
			warn $intspan_align->as_string() if $debug;
			my $intersect_start = $c->{start}->intersection($intspan_align);
			my $intersect_end = $c->{end}->intersection($intspan_align);
			my $cigar_line =$cigar_expand;
			#next unless $cigar_line =~/^S/;

			my @cigar_array = split('',$cigar_line);
			my @query_dna = split('',$array->[9]);
			
			
	
    		#	my @qpos = split("-",$ppos[0]);
    		#	my @end_pos = split("-",$ppos[1]);	
    		#	next unless $end_pos[1];
    			
    			my $len = scalar($intersect_start->as_array) ;
    			
    	
    			my $pos_dna =0;
    			my $i =0;
    			my $new_cigar;
    			
			my $decIndel =0;
    			while($pos_dna < $len){
    				if 	($cigar_array[$i] eq 'D'){
    					$decIndel ++;
    					$i++;
    					next;
    				} 
    				if 	($cigar_array[$i] eq 'I'){
    					$decIndel --;
    				}
    				$i++;
    				$new_cigar.= "S";
    				$query_dna[$pos_dna] = "N";
    				$pos_dna ++;
    				
    			}
    			#middle 
    			my $len_middle = scalar(@query_dna) - scalar($intersect_end->as_array);
    			while($pos_dna < $len_middle){
    				$new_cigar.=$cigar_array[$i] ;
    				if 	($cigar_array[$i] eq 'D'){
    					$i++;
    					next;
    				}
    				 $i++;
    				$pos_dna ++;
    				
    			}
    			
    			#end 
			 	while($pos_dna < scalar(@query_dna)){ 
			 		if 	($cigar_array[$i] eq 'D'){
    					$i++;
    					next;
    				}
    				$i++;
    				$new_cigar.= "S";
    				$query_dna[$pos_dna] = "N";
    				$pos_dna ++;
			 		
			 	}
    		
    			my $cigar = restore_cigar($new_cigar);
    			my ($cigar_expand2,$cigar_array11,$dec1,$len1) =transform_cigar($cigar);
    			my $real_start = $new_start + $dec1;
    			 if ($len1 ne length($array->[9])) {
    			 		warn $array->[9];
    			 	warn $array->[5];
    			 	warn join("\t",@$array);
    			 	warn $intersect_start->as_string;
    			 	warn $intersect_end->as_string;
    			 	warn $new_cigar;
    			 	warn $array->[9];
    			 	warn $array->[5];
    			 	warn $cigar;
    			 	return -2;
    			 	
    			 }
    			# 	warn $new_cigar;
    			 #		warn $array->[9];
    			 	#	warn $cigar;
    			 my $newseq = 	join('',@query_dna);
    		
    			 return -3 if $len1 ne length($newseq);
    			 return -4  if  length($newseq) ne length($array->[10]) ;		
    			 my $c = $newseq =~ tr/N//;
    			  return -5  if abs($c - length($newseq)) < 2 ;	
    			  
    		
    			 
    		  	$array->[3] = $new_start + $dec1+$decIndel;
    			$array->[5] = $cigar;
    			$array->[9] =  join('',@query_dna);
    				 $seq_dejavu->{$hids}->{position} = $array->[3] ;
    			 $seq_dejavu->{$hids}->{cigar} = $cigar;
    			 
    			 $seq_dejavu->{$hids}->{mask_seq} = $array->[9];
    	
    		
    			$test2 ++;# unless $intersect_end->as_string() && $pos_key eq "start";
    			warn Dumper $array if $debug;
    			die() if $debug;
    			return 1;
    		# middle alignement
	}
}


$| =1;


open(BAM,"samtools view  -h $bam_file   |");

my $pm1 = new Parallel::ForkManager($fork);

 $pm1->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
		unless (defined($data) or $exit_code>0){
				print qq|No message received from child process $exit_code $pid!\n|; 
				#die();
				return;
			}
			my $array = $data->{array};
		
			foreach my $bam (@$array){
				print join("\t",@$bam)."\n";
			}
    }
  );




my $last =1; 
while($last){
my @lines;
	while (my $bam = <BAM>){
		if ($bam =~ /^@/){
			print $bam;
			next;
			}
		chomp($bam);
		push(@lines,$bam);
		last if scalar(@lines) == 5000;
	}
	$last = undef if scalar(@lines) < 5000;
		 my $pid = $pm1->start and next;
		 warn "run 1 ".scalar(@lines);
		 my $res;	
		foreach my $bam (@lines){
			my @array = split(" ",$bam);
			my $flag = $array[1];
		 	if ($array[2] eq "*"){
		 		push(@$res,\@array);
		 			#print join("\t",@array)."\n";
					next;		 			
		 	}
			my $ok = line_forward(\@array);
			push(@$res,\@array);
		}
		
		$pm1->finish(0,{array=>$res});
}	
 $pm1->wait_all_children();
 close(BAM);

exit(0);
my $sam = Bio::DB::Sam->new(-fasta=>"/data-xfs/public-data/HG19/genome/fasta/all.fa",
                             -bam  =>$bam_file
                             );
                             
                             

  my $bam          = Bio::DB::Bam->open($bam_file);
   

  
 foreach my $chr (@{$sam->header->target_name}){	

  $sam->fetch("chr1:1-179520300",\&parse_align);
  die($chr);
 }
                     

die();

