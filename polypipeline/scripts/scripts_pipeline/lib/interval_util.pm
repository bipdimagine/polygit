package interval_util;
use FindBin qw($Bin);
use strict;
use colored;
#use Set::IntSpan;
use Carp;
use Storable qw(store retrieve freeze thaw);
use Set::IntSpan::Fast::XS;
use List::Util qw(sum uniqnum);
 use Digest::MD5 qw(md5 md5_hex md5_base64);
use Path::Tiny;
use Digest::MD5 qw(md5_hex);
use Set::IntervalTree;
use Data::Dumper;
use List::MoreUtils qw(natatime);
 use Digest::MD5::File qw (file_md5_hex);
 
sub size_intspan {
	my ($intspan,$name) =@_;
		my $iter = $intspan->iterate_runs();
	my $len =0;
    my $hpos;
	while (my ( $from, $to ) = $iter->()) {
		#warn $from." ".$to." ".abs($from-$to);
		$len += abs($from-$to)+1;
		push(@$hpos,{start=>$from,end=>$to,name=>$name=>1});
	}
	return ($hpos,$len);
}
 
 sub construct_regions_from_intspan {
 	my ($intspan,$name) =@_;
 	my $iter = $intspan->iterate_runs();
 	my $hpos =[];
 	while (my ( $from, $to ) = $iter->()) {
	#	$end = $end +1 if $end == $start;
		push(@$hpos,[$name,$from,$to]); #$from,end=>$to,name=>$name=>1});
	}
	return ($hpos);
 }

sub coverages_regions {
	my ($regions) = @_;
	
	my $big_tree = Set::IntervalTree->new;
	my @stack;
	my $hpos_bgt;
	my @values ;
#	warn ref($regions->[0]);
	my $rf = -1;# 1592940;
	@$regions = sort{$a->[2] <=> $b->[2]} @$regions;
	
	  my $it = natatime 1000000 , @$regions;
	  
	 while (my @vals = $it->()){  
	 my $tree = Set::IntervalTree->new;
	 if (ref($regions->[0]) eq "HASH"){
	 		foreach my $r (@vals){
	 			
				$tree->insert($r->{name},$r->{start},$r->{end}+1) ;
				push(@values,$r->{start},$r->{end});
	 		}
	}
	else {
		foreach my $r (@vals){
			warn Dumper $r if $r->[1] eq $rf or  $r->[2] eq $rf ;
		
			$tree->insert($r->[0],$r->[1],$r->[2]+1);
			#$r->[2]++ if $r->[2] == $r->[1];
			push(@values,($r->[1],$r->[2]+1));
		}
	}
	push(@stack,{tree=>$tree,start=>$vals[0]->[1],end=>$vals[-1]->[2]});
	 }
	warn "cov reg 1";
	my @ov = uniqnum sort{$a <=> $b} @values;
	@values = ();
	$regions = undef;
	my $last_end = -1;
	my @fr;
	my $gg;
	my $t = time;
	#warn Dumper @stack;
	my $htree = shift(@stack);
	my $limit = scalar(@ov)-1;
	my $previous = 0;
	my $ii =0;
	
	for (my $i=0;$i<$limit;$i++){
		$ii++;
	
		my $start = $ov[$i];
		my $end = $ov[$i+1];
		my $debug;
		
		my $search_pos =  $start+int(abs($end-$start)/2);
		if ($search_pos>$htree->{end}){
			$htree = shift(@stack);
		}
			die($search_pos) unless $htree;
		my $z = $htree->{tree}->fetch($search_pos,$search_pos+1);
		next unless @$z;
		push(@fr, [join(";" , @$z),$start,$end-1]);
	#	push(@fr, ["1",$start,$end-1]);
		#warn "$i/$limit ".abs(time -$t) if $ii%1000000 == 0;
	} 
	
	return \@fr;
	
	
	foreach my $r (@fr){
		print "chr1\t".$r->{start}."\t".$r->{end}."\t".$r->{name}."\n";
	}
	die();
	
}

sub coverages_regions2 {
	my ($regions,$key,$chr) = @_;
	warn $key;
	my $big_tree = Set::IntervalTree->new;
	my @stack;
	my $hpos_bgt;
	my @values ;
	#	warn ref($regions->[0]);
	my $rf = -1;# 1592940;

	@$regions = sort{$a->[2] <=> $b->[2]} @$regions;

	  my $it = natatime 1000000 , @$regions;
	  
	 while (my @vals = $it->()){  
	 my $tree = Set::IntervalTree->new;
	 if (ref($regions->[0]) eq "HASH"){
	 		foreach my $r (@vals){
	 			
				$tree->insert($r->{name},$r->{start},$r->{end}+1) ;
				push(@values,$r->{start},$r->{end});
	 		}
	}
	else {
		foreach my $r (@vals){
			warn Dumper $r if $r->[1] eq $rf or  $r->[2] eq $rf ;
		
			$tree->insert($r->[0],$r->[1],$r->[2]+1);
			#$r->[2]++ if $r->[2] == $r->[1];
			push(@values,($r->[1],$r->[2]+1));
		}
	}
	push(@stack,{tree=>$tree,start=>$vals[0]->[1],end=>$vals[-1]->[2]});
	 }
	warn "cov reg 1";
	my @ov = uniqnum sort{$a <=> $b} @values;
	@values = ();
	$regions = undef;
	my $last_end = -1;
	my @fr;
	my $gg;
	my $t = time;
	#warn Dumper @stack;
	my $htree = shift(@stack);
	my $limit = scalar(@ov)-1;
	my $previous = 0;
	my $ii =0;
	 my $f = "/tmp/".$key.".bed";
	 my $xx =0;
  open (BED,">".$f);
	for (my $i=0;$i<$limit;$i++){
		$ii++;
		$xx ++;
		my $start = $ov[$i];
		my $end = $ov[$i+1];
		
		my $debug;
		
		my $search_pos =  $start+int(abs($end-$start)/2);
		if ($search_pos>$htree->{end}){
			$htree = shift(@stack);
		}
			die($search_pos) unless $htree;
		my $z = $htree->{tree}->fetch($search_pos,$search_pos+1);
		next unless @$z;
		$end -=1  unless ($start eq $end);
		print BED "$chr\t$start\t".$end."\t".join(";" , @$z)."\n";
		#push(@fr, [join(";" , @$z),$start,$end-1]);
	#	push(@fr, ["1",$start,$end-1]);
		#warn "$i/$limit ".abs(time -$t) if $ii%1000000 == 0;
	} 
	close BED;
	
	return \@fr;
	
	
	foreach my $r (@fr){
		print "chr1\t".$r->{start}."\t".$r->{end}."\t".$r->{name}."\n";
	}
	die();
	
}
sub coverages_regions3 {
	my ($regions,$key,$chr) = @_;
	warn "coverage region ".$key;
	my $big_tree = Set::IntervalTree->new;
	my @stack;
	my $hpos_bgt;
	my @values ;
	#	warn ref($regions->[0]);
	my $rf = -1;# 1592940;

	@$regions = sort{$a->[2] <=> $b->[2]} @$regions;

	  my $it = natatime 1000000 , @$regions;
	  
	 while (my @vals = $it->()){  
	 my $tree = Set::IntervalTree->new;
	 if (ref($regions->[0]) eq "HASH"){
	 		foreach my $r (@vals){
	 			
				$tree->insert($r->{name},$r->{start},$r->{end}+1) ;
				push(@values,$r->{start},$r->{end});
	 		}
	}
	else {
		foreach my $r (@vals){
			warn Dumper $r if $r->[1] eq $rf or  $r->[2] eq $rf ;
		
			$tree->insert($r->[0],$r->[1],$r->[2]+1);
			#$r->[2]++ if $r->[2] == $r->[1];
			push(@values,($r->[1],$r->[2]+1));
		}
	}
	push(@stack,{tree=>$tree,start=>$vals[0]->[1],end=>$vals[-1]->[2]});
	 }
	warn "cov reg 1";
	my $intspan_nb = {};
	my @ov = uniqnum sort{$a <=> $b} @values;
	my $max =scalar(@ov);
	my $nn =0;
	@values = ();
	$regions = undef;
	my $last_end = -1;
	my @fr;
	my $gg;
	my $t = time;
	#warn Dumper @stack;
	my $htree = shift(@stack);
	my $limit = scalar(@ov)-1;
	my $previous = 0;
	my $ii =0;
	
	 my $xx =0;
  
	for (my $i=0;$i<$limit;$i++){
		$ii++;
		$xx ++;
		my $start = $ov[$i];
		my $end = $ov[$i+1];
		warn $i."/".$limit if $i%200000 ==0;
		my $debug;
		
		my $search_pos =  $start+int(abs($end-$start)/2);
		if ($search_pos>$htree->{end}){
			$htree = shift(@stack);
		}
			die($search_pos) unless $htree;
		my $z = $htree->{tree}->fetch($search_pos,$search_pos+1);
		next unless @$z;
		my $nb = scalar(@$z);
		unless (exists $intspan_nb->{$nb}){
			$intspan_nb->{$nb} = Set::IntSpan::Fast::XS->new();
		}
		
		$end -=1  unless ($start eq $end);
		$intspan_nb->{$nb}->add_range($start,$end);
		#print BED "$chr\t$start\t".$end."\t".join(";" , @$z)."\n";
		#push(@fr, [join(";" , @$z),$start,$end-1]);
	#	push(@fr, ["1",$start,$end-1]);
		#warn "$i/$limit ".abs(time -$t) if $ii%1000000 == 0;
	} 
	#close BED;
	
	 my $f = "/tmp/".$key.".bed";
	 open (BED,">".$f);
	warn Dumper keys %$intspan_nb;
	foreach my $k (keys %$intspan_nb){
		my $set = $intspan_nb->{$k};
		my $iter = $set->iterate_runs();
		while (my ( $from, $to ) = $iter->()) {
    		print BED "$chr\t$from\t$to\t$k\n";
		}
	}
	close BED;
	return;
	die();
	return \@fr;
	
	
	foreach my $r (@fr){
		print "chr1\t".$r->{start}."\t".$r->{end}."\t".$r->{name}."\n";
	}
	die();
	
}


sub light_coverages_regions {
	my ($tree,$values) = @_;
	

	my $regions;
	my $last_end = -1;
	
	my $limit = scalar(@$values)-1;
	my $previous = 0;
	my $ii =0;
	my $t =time;
	for (my $i=0;$i<$limit;$i++){
	
		my $start = $values->[$i];
		my $end = $values->[$i+1];
		
		my $search_pos =  $start+int(abs($end-$start)/2);
		my $z = $tree->fetch($search_pos,$search_pos+1);
		next unless @$z;
		push(@$regions, [join(";" , @$z),$start,$end-1]);
		#push(@$regions, ["1",$start,$end-1]);
		$ii ++;
		warn "$i/$limit ".abs(time -$t) if $ii%1000000 == 0;
	} 
	
	return $regions;
	
	
	
}
sub divide_by_chunks {
	my ($number,$chunksize) =@_;
	if ($chunksize >= $number){
		my $r = [1,$number];
		return [$r];
	}
	my $interval;
	#my $chunksize = int($number/$parts)+1;             # size of each interval
    my $chunkstart = 1;                        # start of interval
    my $chunkend = $chunkstart + $chunksize -1;  # end of that interval
    while ($chunkstart < $number){            # don't go beyond the range
       push(@$interval,[$chunkstart,$chunkend]);
        $chunkstart += $chunksize;           # take me to beginning of next interval
        $chunkend += $chunksize;             # also tell me where to end that
	 if ($chunkend >= $number)  {          # interval end is beyond the range
             push(@$interval,[$chunkstart,$number]);
            last;                         # we are beyond the range now
        }
}
	return $interval;
}
sub divide {
	my ($number,$parts) =@_;
	my $interval;
	my $chunksize = int($number/$parts)+1;             # size of each interval
    my $chunkstart = 1;                        # start of interval
    my $chunkend = $chunkstart + $chunksize -1;  # end of that interval
    while ($chunkstart < $number){            # don't go beyond the range
       push(@$interval,[$chunkstart,$chunkend]);
         print $chunkstart." ".$chunkend."\n";
        $chunkstart += $chunksize;           # take me to beginning of next interval
        $chunkend += $chunksize;             # also tell me where to end that
	 if ($chunkend >= $number)  {          # interval end is beyond the range
             push(@$interval,[$chunkstart,$number]);
            last;                         # we are beyond the range now
        }
}
	return $interval;
}

sub callable_patient_bam {
	my ($project_name,$patient_name,$force,$chr_name) = @_;
	warn "BAM";
		my $buffer = GBuffer->new();
 	my $sambamba = $buffer->software("sambamba");
	my @values = ("1","2","3","5","10","15","20","30") ;
	
  	my @values2 = sort{$b <=> $a} @values;
  	
	my $project = $buffer->newProject( -name => $project_name );
	my $patients = $project->get_only_list_patients($patient_name);
	my $patient=  $patients->[0];	
	my $bam = $patient->getBamFile();
	my $lmdb_file =  $patient->name.".callable";
	my $chr = $project->getChromosome($chr_name);
	
		my $intspan;
		foreach my $v (@values2){
  				$intspan->{$v} = Set::IntSpan::Fast::XS->new();# unless exists $intspan_type->{$v};
		}
		my $dir_tmp_callable = $project->getCallingPipelineDir("callable");

        my $window_size = 50000000;
		my $windows =  $chr->getWindowCaptureForCalling(250,$window_size);
		foreach my $w (@$windows){
 	  		my @bed = map{$_ = $w->{chromosome}."\t".$_."\n"} split(";",$w->{intspan}->as_string({ sep => ";", range => "\t" }));   
 	  		  my $bed_file = File::Temp->new( TEMPLATE =>"XXXXXXX",
                        DIR => $dir_tmp_callable,
                        SUFFIX => ".bed");
                        
 	  		path($bed_file)->spew(@bed);
 	 		my $qq = qq{$sambamba depth  base -F "mapping_quality >=20  " -L $bed_file   $bam -t 5 2>/dev/null  | cut -f 2,3 | };
 	 		open (BAMBA,$qq);
 	 		 while (<BAMBA>){
 	 			next if $_ =~/POS/;
 	 			chomp();
 	 			my($a,$b) = split(" ",$_);
 	 			foreach my $v (@values2){
 					if ($b >= $v){
 						$intspan->{$v}->add($a);
 						last;
 					}
  				}
 	 	
			}
 	 	
		} #end window
		  for (my $i=0;$i<@values;$i++ ){
  			for (my $j=$i+1;$j<@values;$j++ ){
  				$intspan->{$values[$i]} = $intspan->{$values[$i]}->union($intspan->{ $values[$j] } );
  			}
		  }
		  return $intspan;
		  
	
	
}


sub callable_patient_gvcf {
	my ($project_name,$patient_name,$force,$chr_name) = @_;
		my $buffer = GBuffer->new();
 

	my $project = $buffer->newProject( -name => $project_name );
	my $patients = $project->get_only_list_patients($patient_name);
	my $patient=  $patients->[0];	
	my $bam = $patient->getBamFile();
	my $dir_tmp_callable = $project->getCallingPipelineDir("callable");

	
my $lmdb_file =  $patient->name.".callable";

 my $gvcf = $patient->getGvcfFile();
 #die() unless $gvcf;
 return callable_patient_bam($project_name,$patient_name,$force) unless $gvcf;

#$lmdb->close();
my $tabix = $buffer->software("tabix");
my @values = ("1","2","3","5","10","15","20","30") ;
  	my @values2 = sort{$b <=> $a} @values;
  	my $intspan;

  	open(GVCF, "$tabix  $gvcf $chr_name | grep -v \"\#\" | cut -f 1,2,4,8,9,10 |");
  	
  	foreach my $v (@values2){
  				$intspan->{$v} = Set::IntSpan::Fast::XS->new();
  	}
  				
  	my $last;
  	while (<GVCF>){
  		chomp();
  		my @line = split(" ",$_);
  		my $chr = $line[0];
  		#warn $chr if $chr ne $last;
  		$last = $chr;
  		my $start = $line[1];
  		my $end = -1;
  		my $dp;
  	
  		if ($line[3] =~ /END/){
  			$end = $line[3];
  			$end  =~ s/END\=//;
  		}
  		my @t = split(":",$line[5]);
  		$dp = $t[1];
  		 unless( $end > 0){
  		 	$dp = $t[2];
  		 	$end = $start+(length($line[2])-1);
  		 }
		$dp = $t[2] unless $end > 0;
  			
  		foreach my $v (@values2){
 				if ($dp >= $v){
 					$intspan->{$v}->add_range($start,$end);# if $end > 0;
 					#$intspan->{$v}->add($start) unless $end > 0;
 					last;
 				}
  			}
  		
  	}
  	close (GVCF);

  	for (my $i=0;$i<@values;$i++ ){
  			for (my $j=$i+1;$j<@values;$j++ ){
  				$intspan->{$values[$i]} = $intspan->{$values[$i]}->union($intspan->{ $values[$j] } );
  			}
  	}
  return $intspan;

}



1;