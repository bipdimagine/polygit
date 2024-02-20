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
	confess();
	return $self->_parse_tab_file($patient,$chr);
	
}

sub _parse_tab_file {
	my ($self, $patient, $chr) = @_;
	return [] unless -e $self->file();
	my $tabix = Bio::DB::HTS::Tabix->new( filename => $self->file() );
	my $iter = $tabix->query($chr->fasta_name);
	my ($h_header, $h_global, @l_res);
	my $i = 0;
	return [] unless $iter;
	while ( my $line = $iter->next ) {
		my ($chr_id, $start, $end, $strand, $intron_motif, $is_annot, $new_j, $multiple_j, $max_j) = split(' ', $line);
		my $id = $chr->name.'_'.$start.'_'.$end;
		my ($h, $h_tmp);
		$h_tmp->{$chr_id} = undef;
		$h->{id} = $id;
		#$h->{annex}->{$patient->name()}->{type} = 'dragen';
		$h->{annex}->{$patient->name()}->{junc_new_count_sj} = $new_j;
		$h->{annex}->{$patient->name()}->{junc_multiple_count_sj} = $multiple_j;
		$h->{annex}->{$patient->name()}->{overhang} = $max_j;
		$h->{annex}->{$patient->name()}->{intron_motif} = $intron_motif;
		$h->{annex}->{$patient->name()}->{is_sj} = 1;
		my $genes = $chr->getGenesByPosition($start,$end);
		foreach my $g (@$genes){
			my $gene_junctions = $g->get_canonic_junctions();
#			$gene_junctions->{$id}->{is_sj} = 1;
			if (exists $gene_junctions->{$id}){
				$gene_junctions->{$id}->{count}->{$patient->name} = $new_j;
				$gene_junctions->{$id}->{count_multi}->{$patient->name} = $multiple_j;
			}
		
		}
		push (@l_res, $h);
	}
	close (FILE);
	return \@l_res;
}

sub parse_SJ_file {
	my ($self, $patient, $chr,$hh) = @_;

	my $array =  $self->_parse_tab_file($patient,$chr);
	foreach my $j (@$array){
		my $id = $j->{id}; 
		next unless exists $hh->{$id};
		
		$hh->{$id}->{is_sj} = 1;
		my $pid = $patient->name;
		$hh->{$id}->{genes} = $j->{genes} if exists $j->{genes} ;
		$hh->{$id}->{type_canonic} = $j->{annex}->{$pid}->{intron_motif} unless exists $hh->{$id}->{type_canonic};
		next unless exists $hh->{$id}->{annex}->{$pid};
		
			foreach my $k (keys %{$j->{annex}->{$pid}}){
				$hh->{$id}->{annex}->{$pid}->{$k} = $j->{annex}->{$pid}->{$k};
			}
			
			 $hh->{$id}->{annex}->{$pid}->{is_sj} = 1;
			delete $hh->{$id}->{annex}->{$pid}->{intron_motif};
			
	}
	return $hh;
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
sub return_hash_from_rnaseqsea_line {
	my ($linbe,$header) = @_;
	
}
sub parse_results_global_file {
	my ($self, $file_name, $chr) = @_;
	my ($h_global, @l_res);
	my $tabix = Bio::DB::HTS::Tabix->new( filename => $self->file() );
	my $h_header = parse_header($tabix->header);
	my @header = map{lc($_)} split("\t",$tabix->header);
	my $proj_name = $self->getProject->name();
	my $found;
	my $project = $chr->project;
	foreach my $chr_id (@{$tabix->seqnames()}) {
		$found = 1 if $chr_id eq $chr->name();
	}
	return ($h_header, \@l_res) if not $found;
	
	my $iter = $tabix->query($chr->id);
	my $type = "RI";
	$type = "SE" if $tabix->header =~ /SE/; 
	
	
	while ( my $line = $iter->next ) {
		chomp($line);
		my @l_col = split("\t", $line);
		
		if (not scalar(@l_col) == scalar(@header)) {
			warn "\n$line\n\n";
			confess("\nERROR parsing file... no same nb columns...\nFile: $file_name\n\n".scalar(@l_col)."==".scalar(@header));
		}
		next if $l_col[0] eq 'NA';
		
		my $h_res;
			my $nb_col = 0;
			# parsing ligne
			for  (my $i = 0;$i<@l_col;$i++){
				$h_res->{type_origin_file} = $type;
				
				my $cat = $header[$i];
				$cat =~s/\#//;
				my $value = $l_col[$i];
				$value =~s/ //g;
				$h_res->{$cat} = $value;
			
				$cat = 'start' if (lc($cat) eq 'junc_se_start' or lc($cat) eq 'junc_ri_start');
				$cat = 'end' if (lc($cat) eq 'junc_se_end' or lc($cat) eq 'junc_ri_end');
			}
			$h_res->{isCanonique} = 1 if ($h_res->{type} =~ /canonique/);
			
			if ($type eq "SE"){
					
					$h_res->{end} = delete $h_res->{junc_se_end};
					$h_res->{start} = delete $h_res->{junc_se_start};
					$h_res->{alt_count} = $h_res->{junc_se_count};
					
				}
				else {
					$h_res->{alt_count} = $h_res->{junc_ri_count};
					$h_res->{end} = delete $h_res->{junc_ri_end};
					$h_res->{start} = delete $h_res->{junc_ri_start};
				}
			$h_res->{len} = abs(abs($h_res->{end} - $h_res->{start} +1));	
			$h_res->{canonic_count} = $h_res->{junc_normale_count};
			# junctions en reverse
			if ($h_res->{end} < $h_res->{start}) {
				my $toto = $h_res->{end};
				$h_res->{end} = $h_res->{start};
				$h_res->{start} = $toto;
				$h_res->{strand} = -1;
			}
			
			$h_res->{canonic_count} = $h_res->{junc_normale_count};
			
			my $res = $chr->genesIntervalTree->fetch($h_res->{start},$h_res->{end}+1);
			next if scalar(@$res) >= 2 &&  $h_res->{alt_count} < 5;
			next if scalar(@$res) >= 3 &&  $h_res->{alt_count} < 7;
			next if scalar(@$res) >= 4 &&  $h_res->{alt_count} < 9;	
			
			next if $h_res->{alt_count} <= 2 && $h_res->{len} >= 50000 ;
			next if $h_res->{alt_count} <= 3 && $h_res->{len} >= 100000;
			next if $h_res->{alt_count} <= 4 && $h_res->{len} >= 200000;
			
			#next if $h_res->{alt_count} <5 if h_res->{canonic_count} == 0;
			
			next if (($h_res->{alt_count}+0.001)/($h_res->{canonic_count}+0.001)) <0.01;
		
			my $ensid = $h_res->{'ensid'};
			my $sample = $h_res->{'sample'};
			my $chr_id = $h_res->{'chr'};
			if ($ensid) { $sample =~ s/$ensid\_$chr_id\_//; }
			$sample =~ s/\_$proj_name//;
			# check et regroupe same jonctions RI SE
			
			
			
			my $id = $chr_id.'_'.$h_res->{'start'}.'_'.$h_res->{'end'}.'_'.$h_res->{type_origin_file};
				warn $id if $h_res->{start} == 39550420 && $h_res->{end} == 39564909;
			confess("\n\nERROR: construct junction $id for column chr. Die\n\n") if not exists $h_res->{'chr'};
			confess("\n\nERROR: construct junction $id for column start. Die\n\n") if not exists $h_res->{'start'};
			confess("\n\nERROR: construct junction $id for column end. Die\n\n") if not exists $h_res->{'end'};
			confess("\n\nERROR: construct junction $id for column type_origin_file. Die\n\n") if not exists $h_res->{'type_origin_file'};
			if (not exists $h_global->{$id}) {
				my $h_tmp;
				$h_tmp->{$h_res->{'chr'}} = undef;
				$h_global->{$id}->{id} = $id;
				$h_global->{$id}->{chromosomes_object} = $h_tmp;
				$h_global->{$id}->{ensid} = $ensid;
				$h_global->{$id}->{gene} = $h_res->{'gene'};
				$h_global->{$id}->{chr} = $chr_id;
				$h_global->{$id}->{start} = $h_res->{'start'};
				$h_global->{$id}->{end} = $h_res->{'end'};
				$h_global->{$id}->{alt_count} = $h_res->{alt_count};
				$h_global->{$id}->{canonic_count} = $h_res->{canonic_count};
				$h_global->{$id}->{sj_id} = $chr->name.'_'.($h_res->{start}+1).'_'.($h_res->{end}-1);
				$h_global->{$id}->{rnaseqsea_type} = $type;
				
			}
			if ($ensid) { $sample =~ s/$ensid\_$chr_id\_//; }
			my $patientid = $project->getPatient($sample)->name;
			$h_global->{$id}->{annex}->{$patientid} = $h_res;
			$h_global->{$id}->{isCanonique} = 1 if exists $h_res->{isCanonique} and $h_res->{isCanonique};
	}
	close (FILE);

	foreach my $id (keys %{$h_global}) {
		my $is_ok;
		foreach my $patient (@{$self->getProject->getPatients()}) {
			next if not exists $h_global->{$id}->{annex}->{$patient->name};
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