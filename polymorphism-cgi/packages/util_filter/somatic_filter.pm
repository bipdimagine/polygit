package somatic_filter;
use strict;
use Set::Intersection;
use Moo;
 use Bio::DB::Sam;
use Data::Dumper;
use base ("root_filter");
use Storable 'dclone';
use fisher;



sub dbl_evt{
	my $self = shift;
		my %keep_vars;
	foreach my $g (values %{$self->groups}) {
			my $tissue_germinal = $g->{germinal};
			my $tissue_somatic =$g->{somatic};
			my %hsnpsain;
			foreach my $name (@$tissue_germinal){
					foreach my $v (@{$self->get($name)}){
						$hsnpsain{$v} ++;
						
					}
			}
			my %hsnpcommon;
			my %hsnponly;
			foreach my $name (@$tissue_somatic){
					foreach my $v (@{$self->get("heterozygote_".$name)}){
						if (exists $hsnpsain{$v}){
							$hsnpcommon{$v} ++
						}
						else {
							$hsnponly{$v} ++;
						}
						
					}
			}
			
			my @all = (@$tissue_germinal,@$tissue_somatic);
			if (keys %hsnpcommon ==0 || keys %hsnponly == 0){
				foreach my $name (@all){
						$self->delete_patient($name);
				} 
			}
			else {
				my @snp1 = keys %hsnpcommon;
				my @snp2 = (@snp1,keys %hsnponly);
					foreach my $name (@$tissue_germinal){
						$self->put($name,\@snp1);
					}
					foreach my $name (@$tissue_somatic){
						$self->put($name,\@snp2);
					}
					foreach my $s (@snp2){
						$keep_vars{$s} = undef;
					}
				
			}
			
	}
	$self->change_hashvar(\%keep_vars);
	$self->refresh_data();
	return;
}

sub loh {
	my $self = shift;
	my %keep_vars;
	foreach my $g (values %{$self->groups}) {
			my $tissue_germinal = $g->{germinal};
			my $tissue_somatic =$g->{somatic};
			my $hsnps_he;
			my $hsnps_ho;
			my $hsnps_all;
			
			foreach my $name (@$tissue_somatic){
					foreach my $v (@{$self->get("heterozygote_".$name)}){
						$hsnps_he->{$v} = $name;
						$hsnps_all->{$v} = undef;
					}
					foreach my $v (@{$self->get("homozygote_".$name)}){
						$hsnps_ho->{$v} = $name;
						$hsnps_all->{$v} = undef;
					}
			}
			
			
         			
			my @keep_snps;
			foreach my $name (@$tissue_germinal){
				
						foreach my $v (@{$self->get("heterozygote_".$name)}){
							
							
							
							if (exists $hsnps_ho->{$v} ){
								#ok j'en ai une
								push(@keep_snps,$v);
									$keep_vars{$v} = undef;
								
							}
							elsif (exists $hsnps_he->{$v} ){
								my $debug;
								$debug =1 if $v eq "21_9683195_G_A";
								my ($t,$a,$c)  = split(":",$self->get_sample_variations($name)->{$v});
								my ($t1,$b,$d) = split(":",$self->get_sample_variations($hsnps_he->{$v})->{$v});
								my $p = fisher::fishers_exact( $a, $b, $c,$d,1);
								warn "$a/$c $b/$d" if $debug;
								warn $p if $debug;
							#	die() if $debug;
								next if $p >0.01;
 								push(@keep_snps,$v);
 								$keep_vars{$v} = undef;
							}
							
						}
							$self->put($name,\@keep_snps);
						
			}
				foreach my $name (@$tissue_somatic){
					$self->put($name,\@keep_snps);
			}
	}
	
	$self->change_hashvar(\%keep_vars);
	$self->refresh_data();
	return;
	
}

1;

