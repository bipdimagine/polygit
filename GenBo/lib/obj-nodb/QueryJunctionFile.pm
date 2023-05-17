package QueryJunctionFile; 

use strict;
use Moo;
use Data::Dumper;


has file => (
	is		=> 'rw',
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
		return $self->getProject->buffer();
	}
);

has project => (
	is		=> 'ro',
	reader	=> 'getProject',
	weak_ref=> 1,
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
	my ($h_header_ri, $list_res_ri) = $self->parse_results_global_file($self->file());
	return $list_res_ri;
}

sub parse_file_SE {
	my ($self) = @_;
	my @lJunctions;
	my ($h_header_se, $list_res_se) = $self->parse_results_global_file($self->file());
	return $list_res_se;
#	my ($h_header_se, $list_res_se) = $self->parse_results_file($self->file());
#	foreach my $h_junction (@$list_res_se) {
#		next if $h_junction->{'chr'} eq 'na';
#		next if $h_junction->{'junc_se_start'} eq 'na';
#		next if $h_junction->{'junc_se_end'} eq 'na';
#		my $h;
#		$h->{id} = $h_junction->{'chr'}.'_'.$h_junction->{'junc_se_start'}.'_'.$h_junction->{'junc_se_end'}.'_SE';
#		$h->{chromosomes_object} = { $h_junction->{'chr'} => undef };
#		$h->{start} = $h_junction->{'junc_se_start'};
#		$h->{end} = $h_junction->{'junc_se_end'};
#		#delete $h_junction->{'chr'};
#		delete $h_junction->{'junc_se_start'};
#		delete $h_junction->{'junc_se_end'};
#		$h->{annex} = $h_junction;
#		push(@lJunctions, $h);
#	}
#	return \@lJunctions;
}

sub parse_results_global_file {
	my ($self, $file_name) = @_;
	open (FILE, $file_name);
	my ($h_header, $h_global, @l_res);
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
					if (lc($cat) eq 'sample') {
						my ($h_tmp, $h_annex);
						my ($id, $start, $end);
						if (exists $h_res->{'junc_se_start'} and exists $h_res->{'junc_se_end'}) {
							$id = $h_res->{'chr'}.'_'.$h_res->{'junc_se_start'}.'_'.$h_res->{'junc_se_end'}.'_SE';
							$start = $h_res->{'junc_se_start'};
							$end = $h_res->{'junc_se_end'};
							delete $h_res->{'junc_se_start'};
							delete $h_res->{'junc_se_end'};
							$h_res->{type_origin_file} = 'SE';
						}
						if (exists $h_res->{'junc_ri_start'} and exists $h_res->{'junc_ri_end'}) {
							$id = $h_res->{'chr'}.'_'.$h_res->{'junc_ri_start'}.'_'.$h_res->{'junc_ri_end'}.'_RI';
							$start = $h_res->{'junc_ri_start'};
							$end = $h_res->{'junc_ri_end'};
							delete $h_res->{'junc_ri_start'};
							delete $h_res->{'junc_ri_end'};
							$h_res->{type_origin_file} = 'RI';
						}
						my $chr_id = $h_res->{'chr'};
						my $ensid = $h_res->{'ensid'};
						$h_tmp->{$chr_id} = undef;
						$h_global->{$id}->{id} = $id;
						$h_global->{$id}->{chromosomes_object} = $h_tmp;
						$h_global->{$id}->{ensid} = $h_res->{'ensid'};
						$h_global->{$id}->{gene} = $h_res->{'gene'};
						$h_global->{$id}->{chr} = $chr_id;
						$h_global->{$id}->{start} = $start;
						$h_global->{$id}->{end} = $end;
						delete $h_res->{'chr'};
						my $sample = $res;
						my $proj_name = $self->getProject->name();
						if ($ensid) { $sample =~ s/$ensid\_$chr_id\_//; }
						$sample =~ s/\_$proj_name//;
						$h_global->{$id}->{annex}->{$sample} = $h_res;
						$h_res = undef;
					}
					
					else {
						$h_res->{lc($cat)} = $res;
					}
					$nb_col++;
				}
			}
		}
		$i++;
	}
	close (FILE);

	foreach my $id (keys %{$h_global}) {
		my $is_ok;
		foreach my $patient (@{$self->getProject->getPatients()}) {
			next if not exists $h_global->{$id}->{annex}->{$patient->name()};
			$is_ok = 1;
		}
		push (@l_res, $h_global->{$id}) if $is_ok;
	}

	return ($h_header, \@l_res);
}

sub parse_results_file {
	my ($self, $file_name) = @_;
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