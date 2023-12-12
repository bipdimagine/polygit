package QueryJunctionFile; 

use strict;
use Moo;
use Carp;
use Data::Dumper;
use Bio::DB::HTS::Tabix;


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

#1. contig name
#2. first base of the splice junction (1-based)
#3. last base of the splice junction (1-based)strand (0: undefined, 1: +, 2: -)
#4. strand (0: undefined, 1: +, 2: -)
#5. intron motif: 0: noncanonical, 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
#6. 0: unannotated, 1: annotated, only if an input gene annotations file was used
#7. number of uniquely mapping reads spanning the splice junction
#8. number of multimapping reads spanning the splice junction
#9. maximum spliced alignment overhang

#Ex: 1	3277541	3283661	2	2	0	3	0	49
sub parse_dragen_file {
	my ($self, $patient, $chr) = @_;
	my $tabix = Bio::DB::HTS::Tabix->new( filename => $self->file() );
	my $iter = $tabix->query($chr->name);
	my ($h_header, $h_global, @l_res);
	my $i = 0;
	while ( my $line = $iter->next ) {
		my ($chr_id, $start, $end, $strand, $intron_motif, $is_annot, $new_j, $multiple_j, $max_j) = split(' ', $line);
		my $id = $chr_id.'_'.$start.'_'.$end.'_ALL';
		my ($h, $h_tmp);
		$h_tmp->{$chr_id} = undef;
		$h->{id} = $id;
		$h->{chromosomes_object} = $h_tmp;
#		$h->{ensid} = '';
#		$h->{gene} = '';
		$h->{chr} = $chr_id;
		$h->{start} = $start;
		$h->{end} = $end;
		$h->{annex}->{$patient->name()}->{type} = 'dragen';
		$h->{annex}->{$patient->name()}->{junc_new_count} = $new_j;
		$h->{annex}->{$patient->name()}->{junc_multiple_count} = $multiple_j;
		$h->{annex}->{$patient->name()}->{junc_normale_count} = $max_j;
		push (@l_res, $h);
	}
	close (FILE);
	return \@l_res;
}

sub parse_file {
	my ($self, $chr) = @_;
	return $self->parse_file_RI($chr) if $self->isRI();
	return $self->parse_file_SE($chr) if $self->isSE();
	confess();
}

sub parse_file_RI {
	my ($self, $chr) = @_;
	my @lJunctions;
	my ($h_header_ri, $list_res_ri) = $self->parse_results_global_file($self->file(), $chr);
	return $list_res_ri;
}

sub parse_file_SE {
	my ($self, $chr) = @_;
	my @lJunctions;
	my ($h_header_se, $list_res_se) = $self->parse_results_global_file($self->file(), $chr);
	return $list_res_se;
}

sub parse_results_global_file {
	my ($self, $file_name, $chr) = @_;
	my ($h_global, @l_res);
	my $tabix = Bio::DB::HTS::Tabix->new( filename => $self->file() );
	my $h_header = parse_header($tabix->header);
	my $found;
	foreach my $chr_id (@{$tabix->seqnames()}) {
		$found = 1 if $chr_id eq $chr->name();
	}
	return ($h_header, \@l_res) if not $found;
	my $iter = $tabix->query($chr->name);
	while ( my $line = $iter->next ) {
		my @l_col = @{parse_line($line)};
		if (not scalar(@l_col) == scalar keys %$h_header) {
			warn "\n$line\n\n";
			confess("\nERROR parsing file... no same nb columns...\nFile: $file_name\n\n");
		}
		if ($l_col[0] ne 'NA') {
			my $h_res;
			my $nb_col = 0;
			foreach my $res (@l_col) {
				my $cat = $h_header->{$nb_col};
				$h_res->{type_origin_file} = 'RI' if (lc($cat) eq 'junc_ri_start' or lc($cat) eq 'junc_ri_end');
				$h_res->{type_origin_file} = 'SE' if (lc($cat) eq 'junc_se_start' or lc($cat) eq 'junc_se_end');
				$cat = 'start' if (lc($cat) eq 'junc_se_start' or lc($cat) eq 'junc_ri_start');
				$cat = 'end' if (lc($cat) eq 'junc_se_end' or lc($cat) eq 'junc_ri_end');
				$h_res->{isCanonique} = 1 if (lc($cat) eq 'type' and lc($res) =~ /canonique/);
				$h_res->{lc($cat)} = $res;
				$nb_col++;
			}
			if ($h_res->{'end'} <= $h_res->{'start'}) {
				my $toto = $h_res->{'end'};
				$h_res->{'end'} = $h_res->{'start'};
				$h_res->{'start'} = $toto;
				$h_res->{'strand'} = -1;
			}
			
			my $id = $h_res->{'chr'}.'_'.$h_res->{'start'}.'_'.$h_res->{'end'}.'_'.$h_res->{type_origin_file};
			confess("\n\nERROR: construct junction $id for column chr. Die\n\n") if not exists $h_res->{'chr'};
			confess("\n\nERROR: construct junction $id for column start. Die\n\n") if not exists $h_res->{'start'};
			confess("\n\nERROR: construct junction $id for column end. Die\n\n") if not exists $h_res->{'end'};
			confess("\n\nERROR: construct junction $id for column type_origin_file. Die\n\n") if not exists $h_res->{'type_origin_file'};
			if (not exists $h_global->{$id}) {
				my $h_tmp;
				$h_tmp->{$h_res->{'chr'}} = undef;
				$h_global->{$id}->{id} = $id;
				$h_global->{$id}->{chromosomes_object} = $h_tmp;
				$h_global->{$id}->{ensid} = $h_res->{'ensid'};
				$h_global->{$id}->{gene} = $h_res->{'gene'};
				$h_global->{$id}->{chr} = $h_res->{'chr'};
				$h_global->{$id}->{start} = $h_res->{'start'};
				$h_global->{$id}->{end} = $h_res->{'end'};
			}
			$h_global->{$id}->{annex}->{$h_res->{'sample'}} = $h_res;
			$h_global->{$id}->{isCanonique} = 1 if exists $h_res->{isCanonique} and $h_res->{isCanonique};
		}
	}
	close (FILE);

	foreach my $id (keys %{$h_global}) {
		my $is_ok;
		foreach my $patient (@{$self->getProject->getPatients()}) {
			next if not exists $h_global->{$id}->{annex}->{$patient->name()};
			$is_ok = 1;
		}
		confess("\n\nERROR: Junction $id without sample annex. Die.\n\n") if not $is_ok;
		push (@l_res, $h_global->{$id});
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
	my @lCol = split("\t", $line);
	return \@lCol;
}

sub parse_header {
	my ($line) = @_;
	my $h_header;
	my $nb_col = 0;
	foreach my $cat (@{parse_line($line)}) {
		$cat =~ s/#//;
		$h_header->{$nb_col} = $cat;
		$nb_col++;
	}
	return $h_header;
}
1;