package GenBoFamily;

use strict;
use Moo;
use Data::Dumper;
extends "GenBoSampleGroup";



has pedigree_details => (
	is		 => 'rw',
	required => 1,
);

has isTrio => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return 1 if $self->mother() or  $self->father() or scalar(@{$self->getChildren})>1;
		return 1 if  scalar(@{$self->getChildren})>1;
		return undef unless $self->mother();
		return undef unless $self->father();
		return undef unless $self->children();
		return 1;
	}
);
has isSolo => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		
		return 1 if scalar @{$self->getMembers} == 1;
		return undef;
		
	
	}
);
# hash les infos des enfants pour chaque patient

has children => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my %hash;
		foreach my $childName (keys %{$self->pedigree_details->{children}}) {
			my $sex = $self->pedigree_details->{all}->{$childName}->{sex};
			my $status = $self->pedigree_details->{all}->{$childName}->{status};
			$hash{$childName} = $sex;
		}
		return \%hash;
	},
);

has mother => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return unless (exists $self->pedigree_details->{mother});
		foreach my $patient (@{$self->getPatients()}) {
			return $patient->name() if ($patient->name() eq $self->pedigree_details->{mother});	
		}
		return;
	},
);

has father => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return unless (exists $self->pedigree_details->{father});
		foreach my $patient (@{$self->getPatients()}) {
			return $patient->name() if ($patient->name() eq $self->pedigree_details->{father});	
		}
		return;
	},
);

# hash les infos des parents pour chaque patient
has parents => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my %hash;
		if ($self->mother()) { $hash{$self->mother()} = undef; }
		if ($self->father()) { $hash{$self->father()} = undef; }
		return \%hash;
	},
);

# hash avec les noms des patients sains
has healthy => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my %hash;
		foreach my $sample (keys %{$self->parents_healthy()}) { $hash{$sample} = undef; }
		foreach my $sample (keys %{$self->children_healthy()}) { $hash{$sample} = undef; }
		return \%hash;
	},
);

# hash avec les noms des patients malades
has ill => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my %hash;
		foreach my $sample (keys %{$self->parents_ill()}) { $hash{$sample} = undef; }
		foreach my $sample (keys %{$self->children_ill()}) { $hash{$sample} = undef; }
		return \%hash;
	},
);

# hash avec les noms des parents sains
has parents_healthy => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my %hash;
		foreach my $parName (keys %{$self->parents()}) {
			my $sex = $self->pedigree_details->{all}->{$parName}->{sex};
			my $status = $self->pedigree_details->{all}->{$parName}->{status};
			if ($status eq 1) { $hash{$parName} = $sex; }
		}
		return \%hash;
	},
);

# hash avec les noms des parents malades
has parents_ill => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my %hash;
		foreach my $parName (keys %{$self->parents()}) {
			my $sex = $self->pedigree_details->{all}->{$parName}->{sex};
			my $status = $self->pedigree_details->{all}->{$parName}->{status};
			if ($status eq 2) { $hash{$parName} = $sex; }
		}
		return \%hash;
	},
);

# hash avec les noms des enfants sains
has children_healthy => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my %hash;
		foreach my $childName (keys %{$self->pedigree_details->{children}}) {
			my $sex = $self->pedigree_details->{all}->{$childName}->{sex};
			my $status = $self->pedigree_details->{all}->{$childName}->{status};
			if ($status eq 1) { $hash{$childName} = $sex; }
		}
		return \%hash;
	},
);

# hash avec les noms des enfants malades
has children_ill => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my %hash;
		foreach my $childName (keys %{$self->pedigree_details->{children}}) {
			my $sex = $self->pedigree_details->{all}->{$childName}->{sex};
			my $status = $self->pedigree_details->{all}->{$childName}->{status};
			if ($status eq 2) { $hash{$childName} = $sex; }
		}
		return \%hash;
	},
);

has GeneticModel => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $z =0;
		foreach my $p (@{$self->getParents} ){
			$z ++ if $p->isIll;
		}
		return "dominant" if $z >0;
		return "recessif" ;
	
	},
);

has isDominant => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->GeneticModel eq "dominant";
	},
);

has isRecessif => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->GeneticModel ne "dominant";
	},
);


##### METHODS #####
sub getMembers {
	my $self = shift;
	my @members;
	if ($self->mother()){
		push(@members, $self->getMother);
	}
	if ($self->father()){
		
		push(@members, $self->getFather);
	}
		
		push(@members, @{$self->getChildren});
	return \@members;
}


has getMother => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return unless ($self->mother());
		return $self->project->getPatient( $self->mother);
	},
);

has getFather => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return unless ($self->father());
		return $self->project->getPatient( $self->father);
	},
);

has getParents => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my @lPat;
		if ($self->mother()) { push(@lPat, $self->getMother()); }
		if ($self->father()) { push(@lPat, $self->getFather()); }
		return \@lPat;
	},
);

has getChildren => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my @lPat;
		foreach my $patient (@{$self->getPatients()}) {
			push(@lPat, $patient) if (exists $self->children->{$patient->{name}});
		}
		return \@lPat;
	},
);

has getChild => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my @lPat;
		foreach my $patient (@{$self->getPatients()}) {
			push(@lPat, $patient) if (exists $self->children->{$patient->{name}});
		}
		confess("more than one children in this family ".$self->name()) if scalar(@lPat)>1;
		return $lPat[0];
	},
);

has getParentsHealthy => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my @lPat;
		foreach my $patient (@{$self->getParents()}) {
			push(@lPat, $patient) if ($patient->isHealthy());
		}
		return \@lPat;
	},
);

has getParentsIll => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my @lPat;
		foreach my $patient (@{$self->getParents()}) {
			push(@lPat, $patient) if ($patient->isIll());
		}
		return \@lPat;
	},
);

has getChildrenHealthy => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my @lPat;
		foreach my $patient (@{$self->getChildren()}) {
			push(@lPat, $patient) if ($patient->isHealthy());
		}
		return \@lPat;
	},
);

has getChildrenIll => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my @lPat;
		foreach my $patient (@{$self->getChildren()}) {
			push(@lPat, $patient) if ($patient->isIll());
		}
		return \@lPat;
	},
);

has getHealthy => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my @lPat;
		foreach my $patient (@{$self->getPatients()}) {
			push(@lPat, $patient) if ($patient->isHealthy());
		}
		return \@lPat;
	},
);

has getIll => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my @lPat;
		foreach my $patient (@{$self->getPatients()}) {
			push(@lPat, $patient) if ($patient->isIll());
		}
		return \@lPat;
	},
);
sub getSiblings {
		my ($self,$child) = shift;
		my $children = $self->getChildren();
		return[] if scalar(@$children) < 2;
		my @t ;
		foreach my $c (@$children){
			push(@t,$c) if $c->id ne $child->id;
		} 
		return \@t;
}

sub getHealthySiblings {
		my ($self,$child) = shift;
		my $children = $self->getSiblings();
		return[] if scalar(@$children) < 2;
		my @t ;
		foreach my $c (@$children){
			push(@t,$c) if $c->isHealthy;
		} 
		return \@t;
}

sub getIllSiblings {
		my ($self,$child) = shift;
		my $children = $self->getSiblings();
		return[] if scalar(@$children) < 2;
		my @t ;
		foreach my $c (@$children){
			push(@t,$c) if $c->isHealthy;
		} 
		return \@t;
}

has getSomatics => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my @lPat;
		foreach my $patient (@{$self->getPatients()}) {
			push(@lPat, $patient) if ($patient->isSomatic());
		}
		return \@lPat;
	},
);

has getGerminals => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my @lPat;
		foreach my $patient (@{$self->getPatients()}) {
			push(@lPat, $patient) if ($patient->isGerminal());
		}
		return \@lPat;
	},
);

has orderVntyper => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		#bvbvb
		my $sc =0;
		foreach my $p (@{$self->getMembers}){
			$sc += $p->isKestrel + $p->isadVNTR;
		}
		return $sc;
	}
);


1;