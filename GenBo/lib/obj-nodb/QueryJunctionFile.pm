package QueryJunctionFile; 

use strict;
use Moose;
use Data::Dumper;

has patient => (
	is		=> 'ro',
	isa		=> 'GenBoPatient',
	reader	=> 'getPatient',
	weak_ref=> 1,
);

has file => (
	is		=> 'rw',
	isa		=> 'Str',
	required=> 1,
	
);

has isRI => (
	is		=> 'rw',
	lazy => 1,
	default => 0,
);

has isSE => (
	is		=> 'rw',
	lazy => 1,
	default => 0,
);

has buffer => (
	is	 => 'ro',
	lazy => 1,
	default => sub{ 
		my $self = shift;
		return $self->getPatient->buffer();
	}
);

has project => (
	is	 => 'ro',
	lazy => 1,
	default => sub{ 
		my $self = shift;
		my $project;
		if ($self->noPatient()) { $project = $self->buffer->newProject( -name => 'NGS2015_0794'); }
		else { return $self->getPatient->project(); }
	}
);


sub parse_file {
	my ($self) = @_;
	return $self->parse_file_RI() if $self->isRI();
	return $self->parse_file_SE() if $self->isSE();
	confess();
}

sub parse_file_RI {
	my ($self) = @_;
	my @lJunctions;
	my ($h_header_ri, $list_res_ri) = parse_results_file($self->file());
	foreach my $h_junction (@$list_res_ri) {
		next if $h_junction->{'chr'} eq 'na';
		next if $h_junction->{'junc_ri_start'} eq 'na';
		next if $h_junction->{'junc_ri_end'} eq 'na';
		my ($h, $hchr);
		$h->{id} = $h_junction->{'chr'}.'_'.$h_junction->{'junc_ri_start'}.'_'.$h_junction->{'junc_ri_end'}.'_RI';
		my $h_chr_id;
		$h_chr_id->{$h_junction->{'chr'}} = undef;
		$h->{chromosomes_object} = $h_chr_id;
		$h->{start} = $h_junction->{'junc_ri_start'};
		$h->{end} = $h_junction->{'junc_ri_end'};
		delete $h_junction->{'junc_ri_start'};
		delete $h_junction->{'junc_ri_end'};
		$h->{annex} = $h_junction;
		push(@lJunctions, $h);
	}
	return \@lJunctions;
}

sub parse_file_SE {
	my ($self) = @_;
	my @lJunctions;
	my ($h_header_se, $list_res_se) = parse_results_file($self->file());
	foreach my $h_junction (@$list_res_se) {
		next if $h_junction->{'chr'} eq 'na';
		next if $h_junction->{'junc_se_start'} eq 'na';
		next if $h_junction->{'junc_se_end'} eq 'na';
		my $h;
		$h->{id} = $h_junction->{'chr'}.'_'.$h_junction->{'junc_se_start'}.'_'.$h_junction->{'junc_se_end'}.'_SE';
		$h->{chromosomes_object} = { $h_junction->{'chr'} => undef };
		$h->{start} = $h_junction->{'junc_se_start'};
		$h->{end} = $h_junction->{'junc_se_end'};
		#delete $h_junction->{'chr'};
		delete $h_junction->{'junc_se_start'};
		delete $h_junction->{'junc_se_end'};
		$h->{annex} = $h_junction;
		push(@lJunctions, $h);
	}
	return \@lJunctions;
}

sub parse_results_file {
	my ($file_name) = @_;
	open (FILE, $file_name);
	my ($h_header, @l_res);
	my $i = 0;
	while (<FILE>) {
		my $line = $_;
		chomp($_);
		if ($i == 0) {
			$h_header = parse_header($line);
		}
		else {
			my $h_res;
			my $nb_col = 0;
			my @l_col = @{parse_line($line)};
			if (not scalar(@l_col) == scalar keys %$h_header) {
				warn "\n$line\n\n";
				confess("\nERROR parsing file... no same nb columns...\nFile: $file_name\n\n");
			}
			if ($l_col[0] ne 'NA') {
				foreach my $res (@l_col) {
					my $cat = $h_header->{$nb_col};
					confess("\nERROR parsing file... no header...\nFile: $file_name\n\n") unless ($cat);
					$h_res->{lc($cat)} = $res;
					$nb_col++;
				}
				push(@l_res, $h_res);
			}
		}
		$i++;
	}
	close (FILE);
	return ($h_header, \@l_res);
}

sub parse_line {
	my ($line) = @_;
	chomp($line);
	my @lCol = split(' ', $line);
	return \@lCol;
}

sub parse_header {
	my ($line) = @_;
	my $h_header;
	my $nb_col = 0;
	foreach my $cat (@{parse_line($line)}) {
		$h_header->{$nb_col} = $cat;
		$nb_col++;
	}
	return $h_header;
}
1;