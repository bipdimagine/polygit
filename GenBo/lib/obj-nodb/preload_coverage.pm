package preload_coverage;
use strict;
use Data::Dumper;
#preload_coverage::load_coverage($project,$patients,\@transcripts,$padding,$utr);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

sub computeLowTranscripts{
	my ($project,$transcripts,$no,$intronic,$utr,$padding,$print) = @_;
	print qq{<div  style="visibility: hidden">} if $print;
 	my $hp; 
 	
	foreach my $t (@$transcripts){
		my $tutr = $utr;
		warn $t->id." ".$t->isncRNA();
		my $capture_intspan = $t->getChromosome->getIntSpanCapture();
		$tutr = 1 if ($t->getChromosome->name eq "MT");  
		my $debug;
	#	$debug = 1 if $t->name eq "ENST00000339618";
		my $exons;
		if ($intronic == 1 ){
			$exons  = $t->getAllGenomicsParts();
		}
		else {
	 		my $exons1 = $t->getExons();
	 		foreach my $e (@$exons1){
	 			 if ($e->is_noncoding && $tutr ne 1){
	 			 	next;
	 			 }
	 			my $s1 = $e->getGenomicSpan()->intersection($capture_intspan);
			 	push(@$exons,$e) unless $s1->is_empty;
			 	
	 		}
	 		unless ($exons) {
	 			$exons = $exons1;
	 		}
		}
			
	my @ids = map{$_->id} sort{$a->{start}*$a->strand <=> $b->{start}*$b->strand }	@$exons;
	
	foreach my $p (@{$project->getPatients}){
		my $found;
		my @array_min;
		foreach my $eid (@ids){
		 	my $h = $no->get($p->name,$eid."_".$padding."_".$utr);
		 	push(@array_min,$h->{min});
		}
		#my $min = min(@array_min);
		$hp->{$t->id}->{$p->id} = min(@array_min);
		#$toto->{$t->id}->{$p->id} ++;#  if $found ==1;
	}
	
}


print "</div>" if $print;
return $hp;
return $hp;
	
}

sub load_coverage_transcripts {
	my ($project,$patients,$transcripts,$print) = @_;
	my $no =  $project->noSqlCoverage();
	#warn scalar(@$transcripts);
	my $cpt=0;
foreach my $patient (@$patients){
	foreach my $transcript (@$transcripts){
		my $debug ;
	#	$debug =1 if $transcript->name eq "ENST00000376265";
		$cpt++;
		warn $no->exists($patient->name,$transcript->id) if $debug;
		print"$cpt " if $print && $cpt%10 ==0;
		my $gc_gene;
		unless ($no->exists($patient->name,$transcript->id)){
			my $gc_gene;# = $transcript->getGene->return_raw_coverage_obj($patient);
			unless ($no->exists($patient->name,$transcript->getGene->id)){
				$gc_gene = $transcript->getGene->return_raw_coverage_obj($patient);
				$gc_gene->array();
				$gc_gene->patient(undef);
				$gc_gene->chromosome(undef);
				$no->put($patient->name,$transcript->getGene->id,$gc_gene);
			}
			 $gc_gene = $no->get($patient->name,$transcript->getGene->id);
			 my $cov = $gc_gene->sub_array($transcript->start,$transcript->coverage_end);
			
			my $tr_array = $gc_gene->sub_array($transcript->coverage_start,$transcript->coverage_end);
			my $gc =  GenBoCoverageTabix->new(chromosome=>$transcript->getChromosome, patient=>$patient, start=>$transcript->coverage_start, end=>$transcript->coverage_end,array=>$tr_array);
			
			#warn Dumper $gc->array;
			#$gc->array($tr_array);
			$gc->patient(undef);
			$gc->chromosome(undef);
			$no->put($patient->name,$transcript->id,$gc);
		}

	}
	
	
}
}

sub load_coverage_regions_dup{
	my ($project,$patients,$transcripts,$print) = @_;
	my $no =  $project->noSqlCoverage();
	
	my $hash_dup={};
	#$hash_dup->{$transcript->getChromosome->ucsc_name} = Set::IntSpan::Fast::XS->new();
	my $buffer1 = GBuffer->new();
	my $project1 = $buffer1->newProject(-name=>$project->name);
	my $patients2 = $project1->getPatients();
	foreach my $patient (@$patients2){
		my $filebed =  $patient->project->getVariationsDir("duplicate_region_calling")."/regions/".$patient->name().".dup.bed";
		if (-e $filebed){
	 	 	open (BED,$filebed);
	  		while(<BED>){
	  			chomp();
	  			my ($chr,$start,$end) = split(" ");
	  			unless (exists $hash_dup->{$chr} ){
	  				$hash_dup->{$chr} = Set::IntSpan::Fast::XS->new();
	  			}
	  	
	  			$hash_dup->{$chr}->add_range($start,$end);
	  		}
		}
	}
	$buffer1 = undef;
	$project1 = undef;
	foreach my $transcript (@$transcripts){
			my $chr = $transcript->getChromosome();
			next unless exists $hash_dup->{$chr->ucsc_name};
			my $spanTr  =$transcript->getGenomicSpan->intersection($hash_dup->{$transcript->getChromosome->ucsc_name}); 
			next if $spanTr->is_empty;
	foreach my $patient (@$patients){
	
					next if defined $no->exists($patient->name,$transcript->id."_raw");			
					my $gc =  GenBoCoverageSamtools->new(raw=>1,chromosome=>$transcript->getChromosome, patient=>$patient, start=>$transcript->start, end=>$transcript->end);
					$gc->array();
					$gc->patient(undef);
					$gc->chromosome(undef);
					$no->put($patient->name,$transcript->id."_raw",$gc);
					$no->put($patient->name,$transcript->id."_spandup",$spanTr);
		}
		
	}	
	#$no->close();
}

sub load_coverage_primers {
	my ($project,$patients,$transcripts,$print) = @_;
	my $primers = $project->getPrimers();


		my %hash_primers;
	 	foreach my $primer (@$primers) {
	 		$hash_primers{$primer->id} = $primer
	 	}
foreach my $patient (@$patients){
	 print"."  if $print;
	 	
	 
	 save_coverage_primers($project,$patient,\%hash_primers,$print);
		}
}

sub save_coverage_primers {
	my ($project,$patient,$hash_primers,$print) = @_;
	my $no =  $project->noSqlCoverage();
	my $nb1 =  scalar(keys %$hash_primers);
	 my $nb2 = $no->count_bulk($patient->name,[keys %$hash_primers]);
	 $| = 1;
	 print"." if $print;
		my $cpt =0;;
		
	  unless ($nb1 eq $nb2) {
	foreach my $primer (values %$hash_primers){
		 print"." if $print;
		#next if exists  $hash->{$primer->id}->{mean};
		$primer->cached_statistic_coverage($patient);
		} 
		}
	}
	 

sub load_coverage_list_primers {
	my ($project,$patients,$transcripts,$print) = @_;
	
	my $no =  $project->noSqlCoverage();
	my %hash_primers;
	#step 1
	foreach my $transcript (@$transcripts){
		my $primers = $transcript->getPrimers();
		foreach my $p (@$primers){
			#die if $p->id eq "primerchrX_41379679";
			#warn $p->id;
			$hash_primers{$p->id} = $p;
		}
	
	 foreach my $patient (@$patients){
	 	print"."  if $print;
	 	
		 save_coverage_primers($project,$patient,\%hash_primers,$print);
	 }
	
	}
	
#	die();
	print "*"  if $print;
	my $nb =0;
	foreach my $capture (@{$project->getCaptures}){
		my $list = $capture->getListPrimers();
		 foreach my $patient (@$patients){
		   foreach my $primer_id  (keys%$list) {
		   my $hash = $no->get($patient->name,$primer_id);
		   next if $hash;
		   	$nb ++;
		   	if ($print){
		   		print "." if  $nb%50 == 0;
		   	}
		   	 my ($chr_name,$start,$end) = split(" ",$list->{$primer_id});
		   	 my $chr= $patient->project->getChromosome($chr_name);
		   	 my $hash = $patient->set_statistics_coverage($chr,$primer_id,$start,$end);
		   	 $no->set($patient->name,$primer_id,$hash);
		   		#warn Dumper ($patient->set_statistics_coverage($chr,$primer_id,$start,$end));
		   }
		 }
	
	}

}

sub load_coverage_for_cache1 {
	my ($project,$patients,$transcripts,$padding,$utr,$print) = @_;
	
	my $hres;
	my $patient= $patients->[0];
	my $no =  $project->noSqlCoverage();
my %genes; 
foreach my $transcript (@$transcripts){

	my $gene = $transcript->getGene();

	push(@{$genes{$gene->id}},$transcript);
}
my @paddings= (0,5,10,15,20,30);

#warn scalar(keys %genes);
my $ng =1;
my $dj;
my $cc =0 ;
my $cn = 0;
my $test =0;
foreach my $htranscripts (values %genes){
	my $gene = $htranscripts->[0]->getGene();
	my $gene_coverage = $gene->get_coverage($patient);
	foreach my $transcript (@$htranscripts){
	
		#warn "\t".scalar( @{ $transcript->getAllGenomicsParts()});
		my $parts = $transcript->getAllGenomicsParts();
		foreach my $exon ( @{ $parts}) {
			my $start = $exon->start-50 ;
			$start =1 if $start <= 0;
			my $cover_tmp =  $gene_coverage->coverage($start,$exon->end+50);
			my $cover_exon = GenBoCoverageSamtools->new(chromosome=>$exon->getChromosome, patient=>$patient, start=>$start, end=>$exon->end+50);
			$cover_exon->{array} = $cover_tmp->{array};
			my $hpos;
			foreach my $padding (@paddings){
				
				my $pos = $exon->return_start_end(padding=>$padding);
				my $sstart = $pos->{start};
				my $send = $pos->{end} ;
				push(@{$hpos->{"$sstart;$send"}},$exon->id."_".$padding."_1");
				if ( $exon->intspan_no_utr->is_empty){
					push(@{$hpos->{"-1;-1"}},$exon->id."_".$padding."_0");
					next;
				}
				my $pos = $exon->return_start_end_no_utr(padding=>$padding);
			#	warn $exon->id if $exon->intspan_no_utr->is_empty;
				
			#	next  if $exon->intspan_no_utr->is_empty;
				 $sstart = $pos->{start};
				 $send = $pos->{end} ;
				push(@{$hpos->{"$sstart;$send"}},$exon->id."_".$padding."_0");
			}#end for padding
			foreach my $string (keys  %$hpos){
				my ($start,$end) = split(";",$string);

				my $res2;
				
				if ($start eq -1 && $end eq -1){
					$res2 = {mean=>0,intspan=>$exon->intspan_no_utr,min=>-1};
				}
				else {
				 	$res2 =  $cover_exon->coverage($start,$end);
					delete $res2->{array};
				}
				foreach my $id (@{$hpos->{$string}}){
				$hres->{$id} = $res2;
				$no->set($patient->name,$id,$res2);
				next if $id =~ /intron/;
		
					}
				
			}
			
		}#end for exon
		
		
		
		my @toto;
		 push(@toto,map{$_->id."_20"."_0"} @{ $transcript->getAllGenomicsParts()});
		 push(@toto,map{$_->id."_20"."_1"} @{ $transcript->getAllGenomicsParts()});
		 my $nb_db2 = $no->count_bulk_lite($patient->name,\@toto);
		 foreach my $ui (@toto){
		 	my $tt = $no->get($patient->name,$ui);
		 	#warn Dumper $tt ;
		 	warn $ui unless $tt;
		 	die() unless $tt;
		 }
		 
		# warn  $patient->name()."!-> ".$transcript->id."  -> ".$nb_db2." vs ".scalar(@{ $transcript->getAllGenomicsParts()})." ".scalar(@{ $transcript->getExons()});
		die() if $nb_db2 ne 2* scalar(@{ $transcript->getAllGenomicsParts()});
		
	}#end for transcript
	
	delete $gene->{coverage_obj}->{$patient->id} if exists $gene->{coverage_obj}->{$patient->id};
	
}#end for gene
return $hres;
}


sub load_coverage {
	my ($project,$patients,$transcripts,$padding,$utr,$print) = @_;
	my $no =  $project->noSqlCoverage();

my @ids; 
 my $hids;
   		my $n =0;
  foreach my $transcript (@$transcripts){

	 	foreach my $patient (@$patients) {
	 		$n++ if $no->get($patient->name,$transcript->id."_".$padding."_".$utr."_ok") ne 1;
	 	}
  }
 return if $n ==0; 

#je peux encore amelioré le systeme en ne prenant que les transcripts non ok en compte 
# mais je susi rop faineant pour le faire 
	 foreach my $transcript (@$transcripts){
	 	
	 #next if $no->get($patient->name,$transcript->id."_ok") eq 1;
	 push(@ids,map{$_->id."_".$padding."_".$utr} @{ $transcript->getAllGenomicsParts()});
	 push(@{$hids->{$transcript->id}},map{$_->id."_".$padding."_".$utr} @{ $transcript->getAllGenomicsParts()});
	}
my $nb_exons =  scalar(@ids);
foreach my $patient (@$patients) {
	#
	 print "."  if $print;
	 my $nb_db = $no->count_bulk_lite($patient->name,\@ids);
	 if ($nb_db eq $nb_exons){
	 	# tout le projet a le bon nombre d'exon donc tous les transcripts sont set a ok pour la prochaine fois;
	 	foreach my $transcript (@$transcripts){
	 		$no->set($patient->name,$transcript->id."_".$padding."_".$utr."_ok",1);
	 	}
	 }
	 #my  $nb_db =0;
	 #die();
	  unless ($nb_db eq $nb_exons){
	my $hash2;
	my $nbt = 0;
	my @paddings= (0,5,10,15,20,30);
	my $res;
	my $toto = 0;
	my $time =time;
	foreach my $transcript (@$transcripts){
		 my $nb_db2 = $no->count_bulk_lite($patient->name,$hids->{$transcript->id});
			 print "."  if $print;
		next if $nb_db2 eq scalar(@{ $transcript->getAllGenomicsParts()});
			 print "."  if $print;
	
		$toto ++;
	#	warn $transcript->name();
		#print " $nbt";
		$nbt ++;
	foreach my $exon ( @{ $transcript->getAllGenomicsParts()}) {
			 print "."  if $print;
			my $id = $exon->id."_".$padding."_".$utr;
			my ($mean,$intspan,$min)  =$exon->cached_statistic_coverage_coding(patient=>$patient,padding=>$padding,limit=>1,utr=>$utr);
	}
	
	   delete $transcript->getGene->{coverage_obj}->{$patient->id} if exists $transcript->getGene->{coverage_obj}->{$patient->id};
		
	
	}
	
	}

#	print qq{OK"\}};
#exit(0);
	}
	
}



sub load_coverage_for_cache {
	my ($project,$patients,$transcripts,$padding,$utr,$print) = @_;
	my $no =  $project->noSqlCoverage();
	
foreach my $patient (@$patients){
	 print"."  if $print;;
	 
	my @ids; 
	
	my $hash2;
	my $nbt = 0;
	my @paddings= (0,5,10,15,20,30);
	
	my $res;
	
	
	foreach my $transcript (@$transcripts){
		#print " $nbt";
		$nbt ++;
	foreach my $exon ( @{ $transcript->getAllGenomicsParts()}) {
		
		my $id = $exon->id."_".$padding."_".$utr;
			#foreach my $padding (@paddings){
		 	#my ($mean,$intspan,$min)  =  $exon->cached_statistic_coverage_coding(patient=>$patient,padding=>$padding,limit=>1,utr=>$utr);
			my ($h)  =$exon->computed_statistic_coverage_coding(patient=>$patient,padding=>$padding,limit=>1,utr=>$utr);
			$res->{$id} = $h;
			#computed_statistic_coverage_coding
			#}		
	}
	
	}
	foreach my $transcript (@$transcripts){
	 delete $transcript->getGene->{coverage_obj}->{$patient->id} if exists $transcript->getGene->{coverage_obj}->{$patient->id};
	}
	$no->put_bulk($patient->name,$res);
	}
#	print qq{OK"\}};
#exit(0);
	
	}
	
	


sub load_cnv_score {
	my ($project,$patients,$transcripts,$print) = @_;
	return;
	my $no =  $project->noSqlCoverage();
   my $nb= 0;
	#my $chr = $transcripts->[0]->getChromosome();

	#while (my $primer = $chr->next_primer){
	#	warn $primer->id;
	#}
	#die();
foreach my $tr (@$transcripts){
my $no =  $project->noSqlCoverage();
	my $kyoto_id = "cnv_".$tr->name;	
	warn "cocoi";
	$tr->getChromosome()->getPrimers();
	my @primers = sort{$a->end*$a->strand <=> $b->end*$b->strand } @{$tr->getPrimers()};
	foreach my $primer (@primers){
		warn $primer->level($patients->[0]);
	}
	die();
	my $no = $tr->getChromosome->get_lmdb_cnvs("r");
	
	
	return unless @primers;
	#	warn $tr->name()." ".scalar(@primers);
	foreach my $patient (@{$patients}){
	next if $no->get($patient->name."_cnv",$tr->id."_cnv_ok") eq 1 ;
	my @ids = map{$_->id} @primers;

	my $hash= $no->get_bulk($patient->name."_cnv",\@ids);
		
   if (scalar(keys %$hash) eq scalar(@ids)){
   		 $no->set($patient->name."_cnv",$tr->id."_cnv_ok",1);
   }
	next if scalar(keys %$hash) eq scalar(@ids);
	#	 die($patient->name);
	 print "." if $print;
	foreach my $primer (@primers){
			print "." if $nb%50 ==0 && $print;
			$nb ++;
			#warn  $primer->compute_cnv_score($patient)
			$primer->cached_cnv($patient);
		}
		
	}
	
}
$no->close();
}

1;

	
