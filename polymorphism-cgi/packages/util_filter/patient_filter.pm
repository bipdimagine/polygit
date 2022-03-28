package patient_filter;
use strict;
use Set::Intersection;
#use Moose;
use base ("somatic_filter");



sub identical_variations {
	my ($self,$limit) = @_;
	my %var_ids;
	foreach my $name (@{$self->{patients}}) {
			foreach my $vid  (@{$self->get($name)}){
						$var_ids{$vid} ++;
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

sub identical_genes {
	my ($self,$limit) = @_;
	my $nb_patient = 0;
		foreach my $name (@{$self->patients}) {
			
				 $nb_patient ++ if  scalar(@{$self->get($name)});
		}
			
		if ($nb_patient < $limit){
				$self->change_hashvar({});
		}	
		
}

sub filter_composite {
	my($self,$limit_var,$limit_gene) = @_;
	
	my $debug;
		my $find;
		my $hashv;
		foreach my $name (@{$self->{patients}}) {
			my $vars = $self->get("$name");
			my $nb = scalar(@$vars);
			if ($nb >1) {
				$find ++;
				foreach my $v (@$vars){
					$hashv->{$v} = undef;
				}
			}
			else {
			
				$self->delete_key($name);
				$self->delete_key("heterozygote_".$name);
				$self->delete_key("homozygote_".$name);
				
			}
		}
		unless ($find) {
			$self->change_hashvar({});
			return;
		}
		
		if ($limit_var>0 && $find < $limit_var){
			$self->change_hashvar({});
			return;
		}
		if ($limit_gene>0 && $find < $limit_gene){
			$self->change_hashvar({});
			return;
		}
		$self->change_hashvar($hashv);
				
	$self->refresh_data();
}


#
#
#method and_variations(HashRef :$data, ArrayRef :$search_names) {
#	
#	#my ($data,$search_names) = @_;
#	my %var_ids;
#		foreach my $name (@$search_names) {
#				next unless exists $data->{$name};
#				foreach my $vid (@{$data->{$name}}){
#						$var_ids{$vid} ++;
#				}
#			}
#			my $nb = @$search_names;
#			my @ok_vars = grep {$var_ids{$_}== $nb} keys %var_ids;
#			my %hash;
#			@hash{@ok_vars} = ();
#			return \%hash;
#}
#
sub recessif {
	my $self = shift;
#	my ($data,$pedigree,$attic_patients,$hash_heho_patients,$patients) = @_;
	
		#ENSG00000162591
		my @vars;
		my $find = scalar (@{$self->patients});
		foreach my $name (@{$self->patients}) {
			my $nb_he = scalar(@{$self->get($name)});
			if ($nb_he == 1){
				#(ArrayRef :$lists,ArrayRef :$type,HashRef :$data, HashRef :$heho_patients )
				#$self->make_he_ho([$name],"homozygote");
				my $nb_ho = scalar(@{$self->get("homozygote_".$name)});
				unless ($nb_ho > 0 ) {
					$find --;
					$self->delete_key($name);
					$self->delete_key("heterozygote_".$name);
					$self->delete_key("homozygote_".$name);
					
					next;
				} 
			}
			push(@vars,@{$self->get($name)});	
		
		}
		
		my %hash;
		@hash{@vars} = ();
		$self->change_hashvar(\%hash);
}


sub exclude_he {
	my ($self,$names) = @_; 
	my @delete_keys;
	push(@delete_keys,map {"heterozygote_".$_} @$names);
	 $self->delete_variations(\@delete_keys);
	
	
}

sub exclude_ho {
	my ($self,$names) = @_; 
	
	my @delete_keys;
	push(@delete_keys,map {"homozygote_".$_} @$names);
	 $self->delete_variations(\@delete_keys);
	
}

1;