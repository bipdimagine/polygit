package family_filter;
use strict;
use Set::Intersection;
use Moo;
 use Bio::DB::Sam;
use Data::Dumper;
use base ("patient_filter");
use Storable 'dclone';
use Set::Scalar;


sub _getFamsFromSamples {
	my ($self,$names) = @_; 
	my %fams;
	foreach my $name (@$names) {
			my $fam_name = $self->individual_pedigree->{$name}->{fam};
			push(@{$fams{$fam_name}} , $name);
	}
	return \%fams;
}

sub and_variations {
	my ($self,$names) = @_; 

	my $nb =0;
	my $fams = $self->_getFamsFromSamples($names);

	$self->and_intra_fam_variations($fams);
	$self->refresh_data();

}


sub and_intra_fam_variations {
		my ($self,$fams) = @_;
		foreach my $ped (@{$self->pedigree}) {
			my $fname = $ped->{fam};
			next unless exists $fams->{$fname};
		
			#next if scalar(@{$fams->{$fname}}) == 1;
			my $not;
			my @global ;
			my %hvars;
			my @list;
			foreach my $name (@{$fams->{$fname}}){
				#warn $name;
				my $vars = $self->get($name);
				push(@list,$vars);
			}
			my @vars = get_intersection(@list);
			my %dejavu;
			foreach my $name (@{$fams->{$fname}}){
				$dejavu{$name}++;
				$self->put($name,\@vars);
			}
		
		
			my @all = (@{$self->get_parents($ped)},@{$self->get_children($ped)});
			foreach my $name (@all){
				next if exists $dejavu{$name};
				my @vs = get_intersection(\@vars,$self->get($name));
				$self->put($name,\@vs);
			}
		
		}
		$self->update_gene();
}




sub and_inter_fam_variations{ 
		my ($self,$fams) = @_; 
					
		my $in_fam;
		my @list;
	
		foreach my $ped (@{$self->pedigree}) {
			my $fname = $ped->{fam};
			next unless exists $fams->{$fname};
			my @fams_list=();
			foreach my $name (@{$self->get_patients_pedigree($ped)}){
			
				push(@fams_list,@{$self->get($name)});
			}
			push(@list,\@fams_list);
		}
		
		my @vars = get_intersection(@list);
		#warn Dumper (@vars);
		foreach my $ped (@{$self->pedigree}) {
			my $fname = $ped->{fam};
			#next if  exists $fams->{$fname};
			foreach my $name (@{$self->get_patients_pedigree($ped)}){
				my @vs = get_intersection(\@vars,$self->get($name));
				$self->put($name,\@vs);
			}
			
		}
	$self->update_gene();
}




sub and_genes {
	my ($self,$names) = @_; 
	

	my $fams = $self->_getFamsFromSamples($names);
	
	$self->and_intra_fam_genes($fams);
	$self->refresh_data();
}

sub exclude_variations {
	my ($self,$names) = @_; 
	my $fams = $self->_getFamsFromSamples($names);
	$self->not_intra_fam_variations($fams,"all");
	$self->refresh_data();
}
sub exclude_he {
	my ($self,$names) = @_; 
	

	my $fams = $self->_getFamsFromSamples($names);
	
	$self->not_intra_fam_variations($fams,"heterozygote");
	$self->refresh_data();
}

sub exclude_ho {
	my ($self,$names) = @_; 
	my $fams = $self->_getFamsFromSamples($names);
	$self->not_intra_fam_variations($fams,"homozygote");

	$self->refresh_data();
}


sub not_intra_fam_variations {
		my ($self,$fams,$type) = @_;
		foreach my $ped (@{$self->pedigree}) {
			my $fname = $ped->{fam};
			next unless exists $fams->{$fname};
			my $not;
			my @global ;
			my %hvars;
			my @list;
			my $s = Set::Scalar->new();
			my %dejavu;
			foreach my $name (@{$fams->{$fname}}){
					my $vars;
					my $var_out;
				if ($type eq "heterozygote"){
					
					 $vars = $self->heterozygote($name);
					 $var_out  =  $self->homozygote($name);
				}
				elsif ($type eq "homozygote"){
					 $vars = $self->homozygote($name);
					 $var_out = $self->heterozygote($name);
				}
				else {
					$vars = $self->get($name) ;
					$var_out = [];
				}
			
				$s->insert (@$vars);
				$dejavu{$name}++;
				$self->put($name,$var_out);
			}
			my @all = (@{$self->get_parents($ped)},@{$self->get_children($ped)});
			foreach my $name (@all){
				next if exists $dejavu{$name};
				my $s2 = Set::Scalar->new();
				$s2->insert(@{$self->get($name)});
				my $s3 = $s2->difference($s);
				my @t = $s3->elements;
				$self->put($name,\@t);
			}
		}
		$self->update_gene();
}
sub not_inter_fam_variations{ 
		my ($self,$fams) = @_; 
		my @names;
		foreach my $ped (@{$self->pedigree}) {
			my $fname = $ped->{fam};
			next unless exists $fams->{$fname};
			push (@names,@{$self->get_parents($ped)},@{$self->get_children($ped)});
		}
		$self->SUPER::exclude_variations(\@names);
}
sub not_inter_fam_genes{
	my ($self,$fams) = @_; 
		my @names;
		foreach my $ped (@{$self->pedigree}) {
			my $fname = $ped->{fam};
			next unless exists $fams->{$fname};
			push (@names,@{$self->get_parents($ped)},@{$self->get_children($ped)});
		}
		$self->SUPER::exclude_genes(\@names);
}

sub exclude_genes {
	my ($self,$names) = @_; 
	my $fams = $self->_getFamsFromSamples($names);
	$self->not_intra_fam_genes($fams);
	$self->refresh_data();
}

sub not_intra_fam_genes {
		my ($self,$fams) = @_;
		
		foreach my $ped (@{$self->pedigree}) {
			my $fname = $ped->{fam};
			next unless exists $fams->{$fname};
			my $not;
			foreach my $name (@{$fams->{$fname}}){
				my $vars = $self->get($name);
				if (scalar(@$vars)){
					$not =1;
					last;
				}
			}
			next unless $not;
			
			my @all = (@{$self->get_parents($ped)},@{$self->get_children($ped)});
			foreach my $name (@all){
				$self->put($name,[]);
			}
		}
		$self->update_gene();
}


sub get_patients_pedigree{
	my ($self,$ped) = @_;
	my $fs;
	push(@$fs,@{$self->get_children($ped)});
	push(@$fs,@{$self->get_parents($ped)});
	return $fs;
}
sub and_inter_fam_genes{
		my ($self,$fams) = @_;
		
		my $in_fam =0;
		foreach my $ped (@{$self->pedigree}) {
			my $fname = $ped->{fam};
	
			next unless exists $fams->{$fname};
			
			my $not;
			my @global ;
			
			foreach my $name (@{$self->get_patients_pedigree($ped)}){
				my $vars = $self->get($name);
				if (scalar(@$vars)){
					$in_fam ++;
					last;
				}			
			}
			
		}
	
	#warn $in_fam ." ". scalar(keys %$fams);
		if ($in_fam != scalar(keys %$fams) ){
				$self->change_hashvar({});
				$self->delete();
				return;
		}
		
#		foreach my $ped (@{$self->pedigree}) {
#			my $fname = $ped->{fam};
#			#next if  exists $fams->{$fname};
#			foreach my $name (@{$self->get_patients_pedigree($ped)}){
#				
#					$self->put($name,[]);
#				}
#			}
#			$self->update_gene();
#		}
#	
		
}


sub and_intra_fam_genes{
		my ($self,$fams) = @_;

		foreach my $ped (@{$self->pedigree}) {
			my $fname = $ped->{fam};
			next unless exists $fams->{$fname};
			next if scalar(@{$fams->{$fname}}) == 1;
			my $not;
			my @global ;
			foreach my $name (@{$fams->{$fname}}){
				my $vars = $self->get($name);
				unless (scalar(@$vars)){
					$not++;
					last;
				}			
			}
			if ($not){
				foreach my $name (@{$ped->{samples}}){
				#foreach my $name (@{$fams->{$fname}}){
						$self->put($name,[]);
				}
			}
		}
		$self->update_gene();
}

sub update_gene {
	my ($self) = @_;
	 my $keep_var = {};
	foreach my $name (@{$self->patients}) {
		foreach my $v (@{$self->get($name)}){
			$keep_var->{$v} =undef;
		}
	}
	$self->change_hashvar($keep_var);
}
sub identical_genes {
	my ($self,$limit) = @_;
	my $nb_patient = 0;
	my %fam;
		foreach my $name (@{$self->patients}) {
				my $fam_name = $self->individual_pedigree->{$name}->{fam};
				next if exists $fam{$fam_name} ;
				 $fam{$fam_name} ++ if  scalar(@{$self->get($name)});
				
		}
		if (scalar(keys %fam) < $limit){
		
				$self->change_hashvar({});
		}	

}


sub identical_variations {
	my ($self,$limit) = @_;
	my $by_fam;
	my %var_ids;
	
	foreach my $name (@{$self->{patients}}) {
		 my $fam = $self->individual_pedigree->{$name}->{fam};
			foreach my $vid  (@{$self->get($name)}){
						$var_ids{$vid} ++ unless exists $by_fam->{$fam}->{$vid};
						 $by_fam->{$fam}->{$vid} ++;
			}	
				
		}
			my @toto =  keys %var_ids;
			foreach my $k (@toto) {
				next if $var_ids{$k}>=$limit;
				delete $var_ids{$k};
			}
			
			
		$self->change_hashvar(\%var_ids);
		return;
	
}

sub nb_parents {
	my ($self,$ped) = @ _;
	die() unless $ped;
	my @parents = ("mother","father");
	my $nb =0;
	foreach my $p (@parents) {
		next unless exists $ped->{$p};
		$nb ++;
	}
	return $nb;
}

sub nb_children {
	my ($self,$ped) = @_;
	die() unless $ped;
	my $nb;
		foreach my $name (@{$ped->{child_d}}){
		next if exists $self->attic->{$name};
		$nb++;
		}
		return $nb;
}

sub get_children {
	my ($self,$ped) = @_;
	die() unless $ped;
	my @list ;
	foreach my $name (@{$ped->{child_d}}){
		next if exists $self->attic->{$name};
		push(@list,$name);
	}
	foreach my $name (@{$ped->{child_s}}){
		next if exists $self->attic->{$name};
		push(@list,$name);
	}
	return \@list;
}

sub get_parents {
	my ($self,$ped) = @_;
	die() unless $ped;
	my @parents = ("mother","father");
	my @list ;
	foreach my $p (@parents) {
				next unless exists $ped->{$p};
				my $name = $ped->{$p};
				next if exists $self->attic->{$name};
				push(@list,$name);
				#push(@list,$self->get("heterozygote_".$name));
		}
		return \@list;
}

sub get_children_malade {
		my ($self,$ped) = @_;
		die() unless $ped;
		my @array;
		foreach my $name (@{$self->get_children($ped)}){
			if  ($self->individual_pedigree->{$name}->{status} == 2) {
				push(@array,$name);
			}
		}
		return \@array;
}

sub isMale  {
		my ($self,$name) = @_;
		return   ($self->individual_pedigree->{$name}->{sex} == 1);
}

sub isFemale  {
		my ($self,$name) = @_;
		return   ($self->individual_pedigree->{$name}->{sex} == 2);
}

sub get_parents_malade {
		my ($self,$ped) = @_;
		die() unless $ped;
		my @array;
		foreach my $name (@{$self->get_parents($ped)}){
			if  ($self->individual_pedigree->{$name}->{status} == 2) {
				push(@array,$name);
			}
		}
		return \@array;
}

sub get_children_sain {
		my ($self,$ped) = @_;
		die() unless $ped;
		my @array;
		foreach my $name (@{$self->get_children($ped)}){
			if  ($self->individual_pedigree->{$name}->{status} ne 2) {
				push(@array,$name);
			}
		}
		return \@array;
}

sub get_parents_sain {
	my ($self,$ped) = @_;
	my @array;
	die() unless $ped;
		foreach my $name (@{$self->get_parents($ped)}){
			if  ($self->individual_pedigree->{$name}->{status} ne 2) {
				push(@array,$name);
			}
		}
		return \@array;
	
}
sub delete_data_pedigree {
	my ($self,$names) = @_;
	
	foreach my $name (@$names) {
					$self->delete_patient($name);
		}
}

sub delete_variations_pedigree{
	my ($self,$ped) = @_;
	$self->delete_data_pedigree($self->get_parents($ped));
	$self->delete_data_pedigree($self->get_children($ped));
}

sub dominant {
		my $self = shift;
				
	
		my %keep_vars;
		foreach my $ped (@{$self->pedigree}) {
			my @vars_child;
			my @list;
			my $children_malade = $self->get_children_malade($ped);
			my $parent_malade = $self->get_parents_malade($ped);
			my @patient_sain = (@$children_malade,@$parent_malade);
			my $children_sain = $self->get_children_sain($ped);
			my $parent_sain = $self->get_parents_sain($ped);
			my @patients_sain = (@$children_sain,@$parent_sain);
			my @patients_malade = (@$children_malade,@$parent_malade);
			my %hsnpsain;
			#delete all variation in sain
			foreach my $name (@patients_sain){
					foreach my $v (@{$self->get($name)}){
						$hsnpsain{$v} ++;
						
					}
					$self->delete_patient($name);
			}
			
			#keep only variation in malade and not in sain
			my %hsnps;
			foreach my $name (@$parent_malade){
				 my @vars;
					foreach my $v (@{$self->get($name)}){
						next if exists ($hsnpsain{$v}) ;
						$hsnps{$v} ++;
					
					
						
					}
			}
			my $nb_parent_malade = scalar(@$parent_malade);
			#ben je crois que je vais avoir plusieur cas de figure ici 
			# first one one parent sain + one malade in this case keep only he in children
			foreach my $name (@$children_malade){
				 my @vars;
				 my $lvars ;
				 if ($nb_parent_malade == 1){
				 	$lvars = $self->get("heterozygote_".$name) ;
				 }
				 else {
				 	 	$lvars = $self->get($name);
				 }
				
				 
					foreach my $v (@$lvars){
						next if exists ($hsnpsain{$v}) ;
						$hsnps{$v} ++;
						}
				}
			
			
			
			
			
			# now keep only snp present in all malade
			my $nb = scalar(@patients_malade);
			 my @vars;
			foreach my $v (keys %hsnps){
				next if $hsnps{$v} < $nb;
				$keep_vars{$v} = undef;
				 push (@vars,$v);
			}
			foreach my $name (@patients_malade){
				$self->put($name,\@vars);
			}
			
			
		}
		$self->change_hashvar(\%keep_vars);
	$self->refresh_data();
	return;
}

sub moredenovo1 {
	my $self = shift;
	my @parents = ("mother","father");
	foreach my $ped (@{$self->pedigree}) {
		my $children_malade = $self->get_children_malade($ped);
		foreach my $name (@$children_malade){
			next if exists $self->attic->{$name};
			my @vv;
			foreach my $v (@{$self->get($name)}){
				foreach my $p (@parents) {
					next unless exists $ped->{$p};
					my $name = $ped->{$p};
				 	my $file = $self->project()->getCoverageDir()."$name.cov.gz";
				 	my $tabix = new Tabix(-data =>$file);
				 	my ($chr,$pos,$a1,$a2) = split("_",$v);
				 	my $pos2 = $pos+1;
					 my $res = $tabix->query("chr".$chr,$pos,$pos+1);
					 
					 my @data;
		 			while(my $line = $tabix->read($res)){
				 		warn $line;
		 			}
				}
			}
				
		}
	}
	die();
}


sub moredenovo {
	my $self = shift;
	my @parents = ("mother","father");
	my $aa;
		my %keep_vars;
 my $debug =1;
	foreach my $ped (@{$self->pedigree}) {
		my $children_malade = $self->get_children_malade($ped);
		foreach my $name (@$children_malade){
			next if exists $self->attic->{$name};
			my @vv;
			
			
			
			foreach my $v (@{$self->get($name)}){
				my ($chr,$pos,$a1,$a2) = split("_",$v); 
				if (length($a1) > length($a2)){
					$a2 = "-";
				}
				elsif (length($a1) < length($a2)){
						$a2 = "+";
				}
				$debug =1 if $pos == 32636947;
				
				my  $count = 0;
				foreach my $p (@parents) {
					next unless exists $ped->{$p};
					my $name2 = $ped->{$p};
					my $patient = $self->project->getPatient($name2);
					my $bam = $patient->getBamFile();
					my $depth      = 0;
 					my $positions  = 0;
 					my $sam = $self->sam($name2);
# 					        my $sam = Bio::DB::Sam->new(-bam  => $bam,
#                        
#                             );
                             
          my $limit = 2;
          my $nb_mut = 0;                   
 			my $callback = sub {
         			my ($seqid,$pos1,$pileups) = @_;
         			return  if $pos1 ne $pos ;
         			if (scalar(@$pileups) < 3 ) {
         				$nb_mut =99;
         				return;
         				
         				
         			}
         			foreach my $pileup (@$pileups){
         				 my $b     = $pileup->alignment;
         				  my $qbase  = substr($b->qseq,$pileup->qpos,1);
         				  $nb_mut ++ if $qbase eq $a2;
         				  if ($a2 eq "-"){
         				  	 $nb_mut ++ if $pileup->indel <0;
         				  }
         				  elsif ($a2 eq "+"){
         				  		 $nb_mut ++ if $pileup->indel >0;
         				  }
         				 last if $nb_mut >= $limit;
         				 
         			}
 			};
 			my $end = $pos+1;
 			$sam->fast_pileup("chr$chr:$pos-$end",$callback);
 			$count = 1 if $nb_mut >= $limit; 
 			last if $count > 0;
				}
			
				if ($count == 0) {
					$keep_vars{$v} = undef;
					push(@vv,$v);
				}
				else {
				#	warn $v;
				}
			}
			$self->put($name,\@vv);
				
		}
	}
	#warn Dumper %keep_vars if keys %keep_vars;
	
	$self->change_hashvar(\%keep_vars);
	$self->refresh_data();
	#die() if $debug;
}

sub denovo {
	my $self = shift;
	my @parents = ("mother","father");
	my %keep_vars;
	foreach my $ped (@{$self->pedigree}) {
		my @list;
		foreach my $p (@parents) {
			next unless exists $ped->{$p};
			my $name = $ped->{$p};
			next if exists $self->attic->{$name};
			push(@list,@{$self->get($name)});
			$self->delete_patient($name);
		}
		my $children_malade = $self->get_children_malade($ped);
		my $children_sain = $self->get_children_sain($ped);
		foreach my $name (@$children_sain){
			next if exists $self->attic->{$name};
			#next if  $self->individual_pedigree->{$name}->{status} == 2;
			push(@list,@{$self->get($name)});
			$self->delete_patient($name);
		}
		
		my %hash_parent;
		@hash_parent{@list} =();
		foreach my $name (@$children_malade){
			next if exists $self->attic->{$name};
			my @vv;
			foreach my $v (@{$self->get($name)}){
				next if exists $hash_parent{$v};
				push(@vv,$v);
				$keep_vars{$v} = undef;
			}
			$self->put($name,\@vv);
		}
		
				
	}	
	$self->change_hashvar(\%keep_vars);
	$self->refresh_data();

	return;
	
}
sub recessif_X {
	my $self = shift;
	#die();
	my %keep_var;
	foreach my $ped (@{$self->pedigree}) {
			my @vars_child;
			my @list;
			my $children_malade = $self->get_children_malade($ped);
			my $children_sain = $self->get_children_sain($ped);
			my %v_sain;
			foreach my $name (@$children_sain){
				if ($self->isMale($name)){
				
					foreach my $v (@{$self->get($name)}){
						$v_sain{$v} ++;
					}
				}
				else{
					foreach my $v (@{$self->get("homozygote_".$name)}){
						$v_sain{$v} ++;
					}
				}
				$self->delete_patient($name);
			}
			
			foreach my $name (@{$self->get_parents($ped)}) {
				if ($self->isMale($name)){
					foreach my $v (@{$self->get($name)}){
						$v_sain{$v} ++;
					}
					$self->delete_patient($name);
				}
			}	
			
			my %var_pats;
			foreach my $name (@$children_malade){
				
				my @ll;
					if ($self->isMale($name)){
						foreach my $v (@{$self->get($name)}){
							next if exists $v_sain{$v};
							$var_pats{$v}++;
							push(@ll,$v);
						}
					}
					else {
						foreach my $v (@{$self->get("homozygote_".$name)}){
							next if exists $v_sain{$v};
							$var_pats{$v}++;
							push(@ll,$v);
						}
					}
					
				push(@list,\@ll);
			}
			
				
			#my @toto = keys %var_pats;
			#push(@list,\@toto);	
			
			my $parent_vars;
			
			foreach my $name (@{$self->get_parents($ped)}) {
				unless ($self->isMale($name)){
					push(@list,$self->get("heterozygote_".$name));
				}
			}			
			
			my @vars = get_intersection(@list);
			foreach my $name (@{$self->get_parents($ped)}) {
				unless ($self->isMale($name)){
					$self->put($name, \@vars);
				}
				
			}

			foreach my $name (@$children_malade){
				next if exists $self->attic->{$name};
				$self->put($name,\@vars);
			}
			foreach my $v (@vars){
				$keep_var{$v} ++;
			}
		}
	$self->change_hashvar(\%keep_var);
	
	
}


sub recessif {
		my $self = shift;
		if ($self->chr eq "X") {
			$self->recessif_X();
			return;
		}
		my %keep_var;
		foreach my $ped (@{$self->pedigree}) {
			my $debug;
			$debug =1 if  $ped->{fam} eq "fam_2";
			my @vars_child;
			my @list;
			my $children_malade = $self->get_children_malade($ped);
			my $children_sain = $self->get_children_sain($ped);
			my %v_sain;
			foreach my $name (@$children_sain){
				foreach my $v (@{$self->get("homozygote_".$name)}){
					$v_sain{$v} ++;
				}
				$self->delete_patient($name);
				
			}
			my %var_pats;
			foreach my $name (@$children_malade){
				my @ll;
				foreach my $v (@{$self->get("homozygote_".$name)}){
					
					next if exists $v_sain{$v};
					$var_pats{$v}++;
					push(@ll,$v);
				}
				push(@list,\@ll);
			}
		
				
			#my @toto = keys %var_pats;
			#push(@list,\@toto);	
			
			my $parent_vars;
			
			foreach my $name (@{$self->get_parents($ped)}) {
			
				push(@list,$self->get("heterozygote_".$name));
			}			
						my @vars = get_intersection(@list);
			warn join("\n",@vars)."\n" if $debug && @vars ;
			foreach my $name (@{$self->get_parents($ped)}) {
				$self->put($name, \@vars);
			}

			foreach my $name (@$children_malade){
				next if exists $self->attic->{$name};
				$self->put($name,\@vars);
			}
			foreach my $v (@vars){
				$keep_var{$v} ++;
			}
		}
	$self->change_hashvar(\%keep_var);
	
	
}
 

sub filter_composite {
	my $self = shift;
				
		my %keep_var;
		foreach my $ped (@{$self->pedigree}) {
			my $children_malade = $self->get_children_malade($ped);
			my $children_sain = $self->get_children_sain($ped);
			my $parents = $self->get_parents($ped);
			my $nb_parents = scalar(@$parents);
			my %var_malade;
			my $found;
			foreach my $name (@$children_malade) {
			#	warn Dumper $self->heterozygote($name);
			#	$self->make_he_ho([$name],"heterozygote");
				my %temp;
				foreach my $v (@{$self->heterozygote($name)}){
					$temp{$v}  ++;
				}
				if (scalar keys(%temp) > 1 ) {
					foreach my $v (keys %temp){
							$var_malade{$v}  ++;
					}
				}
				
			}
 			
			if (scalar keys(%var_malade) <= 1) {
				$self->delete_variations_pedigree($ped);	
				# die() if $debug;
				next;
			}
			foreach my $name (@$children_malade) {
				foreach my $v (@{$self->homozygote($name)}){
					delete $var_malade{$v}  if exists $var_malade{$v};
				}
				
			}
			
			my $var_parents;
			my %all_variations_parents;
			
			foreach my $name (@$parents) {
				foreach my $v (@{$self->heterozygote($name)}){
					next unless exists $var_malade{$v};
					$all_variations_parents{$v}++;
					$var_parents->{$name} ->{$v} = undef;
				}
			}
			
			foreach my $name (@$parents) {
				foreach my $v (@{$self->homozygote($name)}){
					#delete 	$all_variations_parents{$v} if  exists $var_malade{$v};
					if (exists  $all_variations_parents{$v}){
						delete 	$all_variations_parents{$v} ;
						foreach my $name (@$parents) {
							delete $var_parents->{$name} ->{$v};
						}
					}
				
				}
			}
			# delete variation in father and mother ;
			foreach my $v (keys %all_variations_parents){
				next if $all_variations_parents{$v} ne 2;
				delete  $all_variations_parents{$v};
				foreach my $name (@$parents) {
					delete $var_parents->{$name} ->{$v};
				}
			}
			
		
			# exists at least one variation in father and one in mother 
			# distribute children variation by parents
			my $v_by_parents;
			foreach my $v (keys %var_malade){
				foreach my $name (@$parents) {
					push(@{$v_by_parents->{$name}},$v) if exists $var_parents->{$name} ->{$v};
				}
			}
			
			# only one parents so create a new one "unkwon"
			if (scalar(@$parents) < 2  ){
				foreach my $v (keys %var_malade){
					next if exists $all_variations_parents{$v};
					push(@{$v_by_parents->{"unknown1"}},$v); 
					push(@{$v_by_parents->{"unknown2"}},$v);
				} 
			
				push(@$parents,"unknown1") ;
				
		
				push(@$parents,"unknown2") if scalar(@$parents) ==1 ;
				
				
			}
			
			my $combinaisons ={};
			my $nbc = 0;
			for my $v1 (@{$v_by_parents->{$parents->[0]}}){
				for my $v2 (@{$v_by_parents->{$parents->[1]}}){
					$combinaisons->{$nbc} = [$v1,$v2];
					$nbc ++;
				}
			}
			
			
			foreach my $name (@$children_sain){
				my $vars = $self->get($name);
				foreach my $nbc (keys %$combinaisons){
					my @res = get_intersection($combinaisons->{$nbc},$vars);
					if (scalar(@res) == 2){
						delete $combinaisons->{$nbc};
					}
	
				}
			}
			 unless (%$combinaisons) {
				$self->delete_variations_pedigree($ped);	
				next;
			 }
			
			my %good_variations;
			foreach my $nbc (keys %$combinaisons){
				foreach my $v (@{$combinaisons->{$nbc}}){
					$good_variations{$v} ++;
					$keep_var{$v} ++;
				}
			}
			foreach my $name (@$children_malade) {
				my $vars = $self->get($name);
				my $res = $self->new_intersection(\%good_variations,$vars);
				if (scalar (@$res) > 1){
					$self->put($name,$res);
				}
				else {
					$self->put($name,[]);
				}
				#$self->put($name,[keys %good_variations]);
			}
			
			foreach my $name (@$children_sain) {
				my $vars = $self->get($name);
				my $res = $self->new_intersection(\%good_variations,$vars);
				$self->put($name,$res);
			}
			foreach my $name (@$parents) {
				next if $name =~ /unkwo/;
				my $vars = $self->get($name);
				my $res = $self->new_intersection(\%good_variations,$vars);
				
				$self->put($name,$res);
			}
			
		}#end pedigree
			$self->change_hashvar(\%keep_var);
		
}


sub mendelian_ho {
	
}

sub loh {

	my $self = shift;
	my @parents = ("mother","father");
	my %keep_vars;
	foreach my $ped (@{$self->pedigree}) {
		my @list;
		
		my $children_malade = $self->get_children_malade($ped);
		my $children_sain = $self->get_children_sain($ped);
		
		foreach my $name (@$children_sain){
			next if exists $self->attic->{$name};
			#next if  $self->individual_pedigree->{$name}->{status} == 2;
			push(@list,$self->get("heterozygote_".$name));
		}
		
		my %hash_parent;
		my @list2;
		foreach my $name (@$children_malade){
			next if exists $self->attic->{$name};
			push(@list,$self->get("homozygote_".$name));
		}
		
		my @vars = get_intersection(@list);
		foreach my $name (@$children_malade,@$children_sain){
			next if exists $self->attic->{$name};
				$self->put($name,\@vars);
		}
		
		foreach my $v (@vars){
				$keep_vars{$v} ++;
		}				
	}	
	$self->change_hashvar(\%keep_vars);
	$self->refresh_data();
	return;
}
#
#method compound (HashRef :$data,
#				ArrayRef :$pedigree, 
#				Str :$attic_patients ,
#				HashRef :$heho_patients,
#				ArrayRef :$patients) {
#					
#					
#					
#					}

1;