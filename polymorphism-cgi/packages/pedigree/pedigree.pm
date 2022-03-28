package pedigree;
use strict;

sub parse_ped {
	my ($file,$patients) = @_;
	
	open(FILE,$file);
	my %ped_bypatient;
	my %ped_byfam;
	while (my $line = <FILE>){
		chomp($line);
		next if $line eq "";
		my $p;
		($p->{fam},$p->{name},$p->{father},$p->{mother},$p->{sex},$p->{status}) = split(" ",$line);
		if ($p->{father} eq "0" && $p->{mother} eq "0") {
			if ($p->{sex} == 1){
				$p->{type} = "F";
			}
			else {
				$p->{type} = "M";
			}
		
		}
		else {
		
			$p->{type} = "C2";
			$p->{type} = "C1" 	if ($p->{sex} == 1);
			
		}
		$ped_bypatient{$p->{name}} = $p;
		push(@{$ped_byfam{$p->{fam}}},$p);
	}
	
	my $peds;
	
	foreach my $fam (keys %ped_byfam){
		my @child_d;
		my @child_s;
		my $pedigree;
		$pedigree->{fam} = $fam;
		foreach my $ind (@{$ped_byfam{$fam}}){
			if ($ind->{mother} ne "0" && $ind->{father} ne "0"){
				push(@child_d,$ind->{name}) if $ind->{status} == 2;
				push(@child_s,$ind->{name}) if $ind->{status} == 1;
			}
			else {
				my $namep = $ind->{name};
				my ($find) = grep {$namep eq $_} @$patients;
				next unless $find;
				push(@{$pedigree->{parents}},$ind->{name});
				$pedigree->{mother} = $ind->{name} if $ind->{sex} == 2;
				$pedigree->{father} = $ind->{name} if $ind->{sex} == 1;
			}
		}
		$pedigree->{child_d} = \@child_d;
		$pedigree->{child_s} = \@child_s;
		push(@$peds,$pedigree);
	}
	
	return (\%ped_bypatient,$peds);
}



1;