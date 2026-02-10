package xls_export;
use strict;
use FindBin qw($Bin);
use lib "$Bin/";
use lib "$Bin/../";
use lib "$Bin/../../";
use List::Util qw( max min sum);
use Carp qw(confess croak);
use Data::Dumper;
use Moo;

use Spreadsheet::WriteExcel;
use Compress::Snappy;
use Storable qw(store retrieve freeze dclone thaw);
use session_export;


has only_main_transcripts => (
	is      => 'rw',
	default => undef,
);

has output_dir => (
	is      => 'rw',
	default => undef,
);

has workbook => (
	is      => 'rw',
	default => undef,
);

has title_page => (
	is      => 'rw',
	default => 'export.xls',
);

has header_columns => (
	is      => 'rw',
	default => undef,
);

has nb_pages => (
	is      => 'rw',
	default => 0,
);

has pages => (
	is      => 'rw',
	default => undef,
);

has can_use_hgmd => (
	is      => 'rw',
	lazy    => 1,
	default => 1,
);

has list_generic_header_junctions => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my @lLinesHeader = ('Junction', 'Type', 'Dejavu', 'Dejavu_Ratio_10', 'Dejavu_Ratio_20', 'Chr', 'Start', 'End', 'Gene', 'Description', 'Phenotypes', 'Transcript', 'Transcript_xref', 'Appris_type', 'Tr_start', 'Tr_end', 'Family', 'Patient', 'Sex', 'Status', 'Dp', 'Nb_new', 'Nb_canonique', 'Ratio', 'Score');
		return \@lLinesHeader;
	}
);

has list_generic_header_cnvs => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my @lLinesHeader = ('Variation', 'Rsname', 'Type', 'Dejavu', 'Chr', 'Start', 'End', 'Gene', 'Description', 'Phenotypes', 'Transcript', 'Transcript_xref', 'Exon_intron', 'Locus', 'Overlap');
		return \@lLinesHeader;
	}
);

has list_generic_header => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my @lLinesHeader = ('Variation', 'Rsname', 'Type', 'Dejavu', 'Chr', 'Position', 'Allele', 'Sequence', 'HGMD_Class', 'Cosmic', 'Cadd', 'Ncboost', 'ClinVar', 'Freq (%)', 'gnomad AC', 'gnomad HO', 'gnomad AN', 'Min_Pop_Freq', 'Max_Pop_Freq', 'Gene', 'Description', 'Phenotypes', 'Consequence', 'Transcript', 'Transcript_Xref', 'Appris', 'Polyphen_Score', 'Sift_Score', 'Max_Splice_ai', 'Promoter_ai', 'Alphamissense', 'Exon', 'Cdna_Pos', 'Cds_Pos', 'Protein', 'AA', 'Nomenclature', 'Prot_Nomenclature');
		return \@lLinesHeader;
	}
);

has hash_variants_global => (
	is      => 'rw',
	lazy    => 1,
	default => sub { return {}; }
);

has hash_cnvs_global => (
	is      => 'rw',
	lazy    => 1,
	default => sub { return {}; }
);

has hash_junctions_global => (
	is      => 'rw',
	lazy    => 1,
	default => sub { return {}; }
);

has hash_specific_infos => (
	is      => 'rw',
	lazy    => 1,
	default => sub { return {}; }
);

has format_types => (
	is      => 'rw',
	lazy    => 1,
	default => sub { return {}; }
);

has hash_except_category_rowspan => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $h;
		$h->{'consequence'}		= undef;
		$h->{'transcript'}		= undef;
		$h->{'transcript_xref'}	= undef;
		$h->{'appris'}			= undef;
		$h->{'polyphen_score'}	= undef;
		$h->{'sift_score'}		= undef;
		$h->{'exon'}			= undef;
		$h->{'cdna_pos'}		= undef;
		$h->{'cds_pos'}			= undef;
		$h->{'protein'}			= undef;
		$h->{'nomenclature'}	= undef;
		$h->{'prot_nomenclature'}= undef;
		$h->{'aa'}				= undef;
		$h->{'family'}			= undef;
		$h->{'patient'}			= undef;
		$h->{'sex'}				= undef;
		$h->{'status'}			= undef;
		$h->{'parent_child'}	= undef;
		$h->{'sex_status_icon'}	= undef;
		$h->{'perc'}			= undef;
		$h->{'model'}			= undef;
		$h->{'promoter_ai'}		= undef;
		$h->{'alphamissense'}	= undef;
		$h->{'max_splice_ai'}	= undef;
		return $h;
	}
);

has hash_format_categories => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $h;
		$h->{'consequence'}->{'equal'}->{'mature mirna'}   = 'get_format_red';
		$h->{'consequence'}->{'equal'}->{'splice acc/don'} = 'get_format_red';
		$h->{'consequence'}->{'equal'}->{'frameshift'}     = 'get_format_red';
		$h->{'consequence'}->{'equal'}->{'stop-gained'}    = 'get_format_red';
		$h->{'consequence'}->{'equal'}->{'(start/stop)-lost'} = 'get_format_red';
		$h->{'consequence'}->{'equal'}->{'ncrna'}         = 'get_format_orange';
		$h->{'consequence'}->{'equal'}->{'splice region'} = 'get_format_orange';
		$h->{'consequence'}->{'equal'}->{'missense'}      = 'get_format_orange';
		$h->{'consequence'}->{'equal'}->{'no-frameshift'} = 'get_format_orange';
		$h->{'consequence'}->{'equal'}->{'pseudogene'}    = 'get_format_green';
		$h->{'consequence'}->{'equal'}->{'utr'}           = 'get_format_green';
		$h->{'consequence'}->{'equal'}->{'synonymous'}    = 'get_format_green';
		$h->{'consequence'}->{'equal'}->{'upstream'}      = 'get_format_green';
		$h->{'consequence'}->{'equal'}->{'downstream'}    = 'get_format_green';
		$h->{'consequence'}->{'equal'}->{'intronic'}      = 'get_format_green';

		$h->{'clinvar'}->{'equal'}->{'drug response'}     = 'get_format_green';
		$h->{'clinvar'}->{'equal'}->{'benign'}            = 'get_format_green';
		$h->{'clinvar'}->{'equal'}->{'likely benign'}     = 'get_format_green';
		$h->{'clinvar'}->{'equal'}->{'likely pathogenic'} = 'get_format_orange';
		$h->{'clinvar'}->{'equal'}->{'pathogenic'}        = 'get_format_red';
		$h->{'clinvar'}->{'equal'}->{' pathogenic'}       = 'get_format_red';

		$h->{'polyphen'}->{'equal'}->{'benign damaging'} = 'get_format_green';
		$h->{'polyphen'}->{'equal'}->{'possibly damaging'} = 'get_format_orange';
		$h->{'polyphen'}->{'equal'}->{'probably damaging'} = 'get_format_red';

		$h->{'sift'}->{'equal'}->{'tolerated'}   = 'get_format_green';
		$h->{'sift'}->{'equal'}->{'deleterious'} = 'get_format_red';

		$h->{'cadd'}->{'>='}->{'10'} = 'get_format_green';
		$h->{'cadd'}->{'>='}->{'20'} = 'get_format_orange';
		$h->{'cadd'}->{'>='}->{'30'} = 'get_format_red';

		$h->{'model'}->{'equal'}->{'mother'}   = "get_format_pink";
		$h->{'model'}->{'equal'}->{'mother_c'} = "get_format_pink";
		$h->{'model'}->{'equal'}->{'father'}   = "get_format_blue";
		$h->{'model'}->{'equal'}->{'father_c'} = "get_format_blue";
		$h->{'model'}->{'regexp'}->{'rece'}    = "get_format_magenta";
		$h->{'model'}->{'regexp'}->{'mosa'}    = "get_format_orange";
		$h->{'model'}->{'regexp'}->{'uni'}     = "get_format_model_cyan";
		$h->{'model'}->{'regexp'}->{'denovo'}  = "get_format_red";
		$h->{'model'}->{'equal'}->{'both'}     = "get_format_default";
		$h->{'model'}->{'equal'}->{'?'}        = "get_format_default";

		$h->{'sex_status_icon'}->{'insert_image'} = "insert_image";

		$h->{'hgmd_class'}->{'equal'}->{'dm'}  = "get_format_red";
		$h->{'hgmd_class'}->{'equal'}->{'dm?'} = "get_format_orange";

		$h->{'cnv_confidence'}->{'equal'}->{'low'} = 'get_format_green';
		$h->{'cnv_confidence'}->{'equal'}->{'medium'} = 'get_format_orange';
		$h->{'cnv_confidence'}->{'equal'}->{'high'} = 'get_format_red';
		$h->{'cnv_confidence'}->{'equal'}->{'-'} = 'get_format_default';

		$h->{'score'}->{'>='}->{'0'} = 'get_format_green';
		$h->{'score'}->{'>='}->{'5'} = 'get_format_orange';
		$h->{'score'}->{'>='}->{'8'} = 'get_format_red';

		$h->{'ratio'}->{'>='}->{'10'} = 'get_format_green';
		$h->{'ratio'}->{'>='}->{'20'} = 'get_format_orange';
		$h->{'ratio'}->{'>='}->{'30'} = 'get_format_red';
		return $h;
	}
);

sub get_format_header {
	my ($self) = @_;
	return $self->format_types->{header} if ( $self->format_types() and exists $self->format_types->{header} );
	$self->{format_types}->{header} = $self->workbook->add_format(
		border    => 1,
		underline => 1,
		align     => 'center'
	);
	$self->{format_types}->{header}->set_bold();
	$self->{format_types}->{header}->set_align('center');
	$self->{format_types}->{header}->set_fg_color('silver');
	return $self->format_types->{header};
}

sub get_format_model_cyan {
	my ($self) = @_;
	return $self->format_types->{normal_cyan} if ( $self->format_types() and exists $self->format_types->{normal_cyan} );
	$self->{format_types}->{normal_cyan} = $self->workbook->add_format(  );
	$self->{format_types}->{normal_cyan}->set_color('cyan');
	return $self->format_types->{normal_cyan};
}

sub get_format_magenta {
	my ($self) = @_;
	return $self->format_types->{normal_magenta} if ( $self->format_types() and exists $self->format_types->{normal_magenta} );
	$self->{format_types}->{normal_magenta} = $self->workbook->add_format(  );
	$self->{format_types}->{normal_magenta}->set_color('magenta');
	return $self->format_types->{normal_magenta};
}

sub get_format_pink {
	my ($self) = @_;
	return $self->format_types->{normal_pink} if ( $self->format_types() and exists $self->format_types->{normal_pink} );
	$self->{format_types}->{normal_pink} = $self->workbook->add_format(  );
	$self->{format_types}->{normal_pink}->set_color('pink');
	return $self->format_types->{normal_pink};
}

sub get_format_blue {
	my ($self) = @_;
	return $self->format_types->{normal_blue} if ( $self->format_types() and exists $self->format_types->{normal_blue} );
	$self->{format_types}->{normal_blue} = $self->workbook->add_format(  );
	$self->{format_types}->{normal_blue}->set_color('blue');
	return $self->format_types->{normal_blue};
}

sub get_format_green {
	my ($self) = @_;
	return $self->format_types->{normal_green} if ( $self->format_types() and exists $self->format_types->{normal_green} );
	$self->{format_types}->{normal_green} = $self->workbook->add_format(  );
	$self->{format_types}->{normal_green}->set_color('green');
	return $self->format_types->{normal_green};
}

sub get_format_orange {
	my ($self) = @_;
	return $self->format_types->{normal_orange} if ( $self->format_types() and exists $self->format_types->{normal_orange} );
	$self->{format_types}->{normal_orange} = $self->workbook->add_format(  );
	$self->{format_types}->{normal_orange}->set_color('orange');
	return $self->format_types->{normal_orange};
}

sub get_format_red {
	my ($self) = @_;
	return $self->format_types->{normal_red} if ( $self->format_types() and exists $self->format_types->{normal_red} );
	$self->{format_types}->{normal_red} = $self->workbook->add_format(  );
	$self->{format_types}->{normal_red}->set_color('red');
	return $self->format_types->{normal_red};
}

sub get_format_default {
	my ($self) = @_;
	return $self->format_types->{normal} if ( $self->format_types() and exists $self->format_types->{normal} );
	$self->{format_types}->{normal} = $self->workbook->add_format(  );
	$self->{format_types}->{normal}->set_color('black');
	return $self->format_types->{normal};
}

sub get_format_model_cyan_merge {
	my ($self) = @_;
	return $self->format_types->{normal_cyan_merge} if ( $self->format_types() and exists $self->format_types->{normal_cyan_merge} );
	$self->{format_types}->{normal_cyan_merge} = $self->workbook->add_format(  );
	$self->{format_types}->{normal_cyan_merge}->set_color('cyan');
	return $self->format_types->{normal_cyan_merge};
}

sub get_format_magenta_merge {
	my ($self) = @_;
	return $self->format_types->{normal_magenta_merge} if ( $self->format_types() and exists $self->format_types->{normal_magenta_merge} );
	$self->{format_types}->{normal_magenta_merge} = $self->workbook->add_format(  );
	$self->{format_types}->{normal_magenta_merge}->set_color('magenta');
	return $self->format_types->{normal_magenta_merge};
}

sub get_format_pink_merge {
	my ($self) = @_;
	return $self->format_types->{normal_pink_merge} if ( $self->format_types() and exists $self->format_types->{normal_pink_merge} );
	$self->{format_types}->{normal_pink_merge} = $self->workbook->add_format(  );
	$self->{format_types}->{normal_pink_merge}->set_color('pink');
	return $self->format_types->{normal_pink_merge};
}

sub get_format_blue_merge {
	my ($self) = @_;
	return $self->format_types->{normal_blue_merge} if ( $self->format_types() and exists $self->format_types->{normal_blue} );
	$self->{format_types}->{normal_blue_merge} = $self->workbook->add_format(  );
	$self->{format_types}->{normal_blue_merge}->set_color('blue');
	return $self->format_types->{normal_blue_merge};
}

sub get_format_green_merge {
	my ($self) = @_;
	return $self->format_types->{normal_green_merge} if ( $self->format_types() and exists $self->format_types->{normal_green_merge} );
	$self->{format_types}->{normal_green_merge} = $self->workbook->add_format(  );
	$self->{format_types}->{normal_green_merge}->set_color('green');
	return $self->format_types->{normal_green_merge};
}

sub get_format_orange_merge {
	my ($self) = @_;
	return $self->format_types->{normal_orange_merge} if ( $self->format_types() and exists $self->format_types->{normal_orange_merge} );
	$self->{format_types}->{normal_orange_merge} = $self->workbook->add_format(  );
	$self->{format_types}->{normal_orange_merge}->set_color('orange');
	return $self->format_types->{normal_orange_merge};
}

sub get_format_red_merge {
	my ($self) = @_;
	return $self->format_types->{normal_red_merge} if ( $self->format_types() and exists $self->format_types->{normal_red_merge} );
	$self->{format_types}->{normal_red_merge} = $self->workbook->add_format(  );
	$self->{format_types}->{normal_red_merge}->set_color('red');
	return $self->format_types->{normal_red_merge};
}

sub get_format_default_merge {
	my ($self) = @_;
	return $self->format_types->{normal_merge} if ( $self->format_types() and exists $self->format_types->{normal_merge} );
	$self->{format_types}->{normal_merge} = $self->workbook->add_format(  );
	$self->{format_types}->{normal_merge}->set_color('black');
	return $self->format_types->{normal_merge};
}

sub open_xls_file {
	my ($self) = shift;
	my $workbook;
	my $title = $self->title_page();
	if ( $self->output_dir() ) {
		$workbook = Spreadsheet::WriteExcel->new( $self->output_dir() . '/' . $title );
	}
	else {
		print "Content-type: application/msexcel\n";
		print "Content-Disposition: attachment;filename=$title\n\n";
		$workbook = Spreadsheet::WriteExcel->new( \*STDOUT );
	}
	$self->workbook($workbook);
}

sub construct_hash_row_span {
	my ( $self, $list_lines_header, $list_lines ) = @_;
	my @listHeader = @$list_lines_header;
	my @listLines  = @$list_lines;
	my $first_col  = lc( $listHeader[0] );
	my $h_row_span;
	foreach my $h_line (@listLines) {
		foreach my $cat (@listHeader) {
			next if ( exists $self->hash_except_category_rowspan->{ lc($cat) } );
			my $id_row = $h_line->{$first_col};
			my $value  = $h_line->{ lc($cat) };
			$value = '.' unless ($value);
			$h_row_span->{$id_row}->{ lc($cat) }->{$value}++;
		}
	}
	return $h_row_span;
}

sub add_page_merged {
	my ( $self, $title, $list_lines_header, $list_lines ) = @_;
	my $h_row_span = $self->construct_hash_row_span( $list_lines_header, $list_lines );
	$self->add_page( $title, $list_lines_header, $list_lines, $h_row_span );
}

sub add_page {
	my ( $self, $title, $list_lines_header, $list_lines, $h_row_span ) = @_;
	$self->{nb_pages}++;
	my $number_page = $self->nb_pages();
	$self->{pages}->{$number_page}->{title} = $title;
	$self->{header_columns}->{$number_page} = $list_lines_header;
	$self->{pages}->{$number_page}->{lines} = $list_lines;
	if ($h_row_span) {
		$self->{pages}->{$number_page}->{use_merge}     = 1;
		$self->{pages}->{$number_page}->{hash_row_span} = $h_row_span;
	}
	else {
		$self->{pages}->{$number_page}->{use_merge} = 0;
	}
}

sub open_page {
	my ( $self, $number ) = @_;
	my $title     = $self->pages->{$number}->{title};
	if (length($title) > 30) {
		$title = substr($title, 0, 29);
	}
	my $worksheet = $self->workbook->add_worksheet($title);
	return $worksheet;
}

sub write_header_on_page {
	my ( $self, $worksheet, $number_page ) = @_;
	my $i = 0;
	my $j = 0;
	foreach my $cat ( @{ $self->header_columns->{$number_page} } ) {
		$worksheet->write( $i, $j, $cat, $self->get_format_header() );
		$j++;
	}
}

sub write_on_page {
	my ( $self, $worksheet, $number_page ) = @_;
	return $self->write_on_page_use_merge( $worksheet, $number_page ) if ( exists $self->pages->{$number_page}->{use_merge} and $self->pages->{$number_page}->{use_merge} == 1 );
	my $i = 1;
	foreach my $h_line ( @{ $self->pages->{$number_page}->{lines} } ) {
		my $j = 0;
		foreach my $cat ( @{ $self->header_columns->{$number_page} } ) {
			my $value = '.';
			$value = $h_line->{ lc($cat) } if ( exists $h_line->{ lc($cat) } );
			if ( exists $self->hash_format_categories->{ lc($cat) } ) {
				my $this_format;
				$this_format = $self->get_specific_format( lc($cat), lc($value) ) if ($value);
				if ($this_format) {
					if ( $this_format eq 'insert_image' ) {
						$worksheet->insert_image( $i, $j, $value, 20, 0, 0.4, 0.4 );
					}
					else {
						$worksheet->write( $i, $j, $value, $self->$this_format );
					}
				}
				else { $worksheet->write( $i, $j, $value ); }
			}
			else {
				$worksheet->write( $i, $j, $value );
			}
			$j++;
		}
		$i++;
	}
}

sub write_on_page_use_merge {
	my ( $self, $worksheet, $number_page ) = @_;
	my $h_intspan;
	my $i = 1;
	my $h_row_span = $self->pages->{$number_page}->{hash_row_span};
	foreach my $h_line ( @{ $self->pages->{$number_page}->{lines} } ) {
		my $var_id = $h_line->{'variation'};
		if ( exists $h_row_span->{$var_id} ) {
			foreach my $cat ( keys %{ $h_row_span->{$var_id} } ) {
				foreach my $value ( keys %{ $h_row_span->{$var_id}->{$cat} } ) {
					$h_intspan->{$var_id}->{ lc($cat) }->{$value}->{'i'} = Set::IntSpan::Fast::XS->new();
				}
			}
		}

		my $j = 0;
		my ( $last_cat_rowspan, $h_last_value_rowspan );
		foreach my $cat ( @{ $self->header_columns->{$number_page} } ) {
			my ( $can_use_i, $max_i );
			if ( exists $h_intspan->{$var_id} ) {
				if ( exists $h_intspan->{$var_id}->{ lc($cat) } ) {
					my $value = '.';
					$value = $h_line->{ lc($cat) } if ( exists $h_line->{ lc($cat) } );
					$value = '.' unless ($value);
					if (exists $h_intspan->{$var_id}->{ lc($cat) }->{$value}->{'max_i'}) {
						$max_i = $h_intspan->{$var_id}->{ lc($cat) }->{$value}->{'max_i'};
						if ( $h_intspan->{$var_id}->{ lc($cat) }->{$value}->{'i'}->contains($i) ) {
							$can_use_i = undef;
						}
						else {
							$can_use_i = 1;
						}
					}
					else {
						my $use_rowspan = $h_row_span->{$var_id}->{ lc($cat) }->{$value};
						$use_rowspan = 1 unless ($use_rowspan);
						my $min = $i + 1;
						$max_i = $i + $use_rowspan - 1;
						if ( $use_rowspan > 1 ) {
							$h_intspan->{$var_id}->{ lc($cat) }->{$value}->{'i'}->add_from_string("$min-$max_i");
						}
						$h_intspan->{$var_id}->{ lc($cat) }->{$value}->{'max_i'} = $max_i;
						$last_cat_rowspan                          = lc($cat);
						$h_last_value_rowspan->{$last_cat_rowspan} = $value;
						$can_use_i                                 = 1;
					}
				}
				elsif ($last_cat_rowspan) {
					if ( $h_intspan->{$var_id}->{ lc($last_cat_rowspan) }->{ $h_last_value_rowspan->{$last_cat_rowspan} }->{'i'}->contains($i) ) {
						$max_i = $h_intspan->{$var_id}->{ lc($last_cat_rowspan) }->{ $h_last_value_rowspan->{$last_cat_rowspan} }->{'max_i'};
						$can_use_i = undef;
					}
					else {
						$max_i     = $i;
						$can_use_i = 1;
					}
				}
				else {
					$max_i     = $i;
					$can_use_i = 1;
				}
			}
			if ($can_use_i) {
				my $value = '.';
				$value = $h_line->{ lc($cat) } if ( exists $h_line->{ lc($cat) } );
				if ( exists $self->hash_format_categories->{ lc($cat) } ) {
					my $this_format;
					$this_format = $self->get_specific_format( lc($cat), lc($value) ) if ($value);
					if ($this_format) {
						if ( $max_i > $i ) {
							if ( $this_format eq 'insert_image' ) {
								$worksheet->insert_image( $i, $j, $max_i, $j, $value, 20, 0, 0.4, 0.4 );
							}
							else {
								$this_format .= '_merge';
								$worksheet->merge_range( $i, $j, $max_i, $j, $value, $self->$this_format() );
							}
						}
						else {
							if ( $this_format eq 'insert_image' ) {
								$worksheet->insert_image( $i, $j, $value, 20, 0, 0.4, 0.4 );
							}
							else {
								$worksheet->write( $i, $j, $value, $self->$this_format() );
							}
						}
					}
					else {
						if ( $max_i > $i ) {
							$worksheet->merge_range( $i, $j, $max_i, $j, $value, $self->get_format_default_merge() );
						}
						else { $worksheet->write( $i, $j, $value ); }
					}
				}
				else {
					if ( $max_i > $i ) {
						$worksheet->merge_range( $i, $j, $max_i, $j, $value, $self->get_format_default_merge() );
					}
					else { $worksheet->write( $i, $j, $value ); }
				}
			}
			$j++;
		}
		$i++;
	}
}

sub get_specific_format {
	my ( $self, $cat, $value ) = @_;
	return unless ($value);
	return if ( $value eq '.' );
	return if ( $value eq '-' );
	return unless ( exists $self->hash_format_categories->{ lc($cat) } );
	foreach my $type_compare ( keys %{ $self->hash_format_categories->{ lc($cat) } } ) {
		if ( $type_compare eq 'equal' ) {
			next unless ( exists $self->hash_format_categories->{ lc($cat) }->{$type_compare}->{ lc($value) } );
			my $format = $self->hash_format_categories->{ lc($cat) }->{$type_compare}->{ lc($value) };
			return $format;
		}
		if ( $type_compare eq '>=' ) {
			foreach my $cat_f ( sort { $b <=> $a } keys %{ $self->hash_format_categories->{ lc($cat) }->{$type_compare} } ) {
				return $self->hash_format_categories->{ lc($cat) }->{$type_compare}->{ lc($cat_f) } if ( int($value) >= $cat_f );
			}
		}
		if ( $type_compare eq 'regexp' ) {
			foreach my $cat_f ( keys %{$self->hash_format_categories->{ lc($cat) }->{$type_compare}}) {
				return $self->hash_format_categories->{ lc($cat) }->{$type_compare}->{ lc($cat_f) } if ( $value =~ /$cat_f/ );
			}
		}
		if ( $type_compare eq 'insert_image' ) {
			return 'insert_image';
		}
	}
	return;
}

sub export {
	my $self = shift;
	$self->open_xls_file();
	foreach my $number_page ( sort { $a <=> $b } keys %{ $self->pages() } ) {
		my $worksheet = $self->open_page($number_page);
		$self->write_header_on_page( $worksheet, $number_page );
		$self->write_on_page( $worksheet, $number_page );
	}
}

sub store_cnvs_infos {
	my ( $self, $list_var_objects, $project, $list_patients ) = @_;
	$project->buffer->dbh_deconnect();
	$project->buffer->dbh_reconnect();
	print '|';
	foreach my $var (@$list_var_objects) {
		my $hash;
		$project->print_dot(5);
		my $chr      = $var->getChromosome();
		my $var_id   = $var->id();
		my $chr_h_id = $chr->id();
		$chr_h_id = '23' if ( $chr->id() eq 'X' );
		$chr_h_id = '24' if ( $chr->id() eq 'Y' );
		$chr_h_id = '25' if ( $chr->id() eq 'MT' );

		if ( exists $self->hash_cnvs_global->{$chr_h_id}->{$var_id} ) {
			$hash->{$chr_h_id}->{$var_id} = $self->hash_cnvs_global->{$chr_h_id}->{$var_id};
		}
		else {
			$hash->{$chr_h_id}->{$var_id}->{'var_id'} = $var->name();
			$hash->{$chr_h_id}->{$var_id}->{'rsname'} = '';
			$hash->{$chr_h_id}->{$var_id}->{'type'} = 'cnv' if ( $var->isCnv() );
			my $h_dejavu   = $var->deja_vu();
			my $nb_project = 0;
			my $nb_patient = 0;
			my $he         = 0;
			my $ho         = 0;
			my $st_project;

			foreach my $projName ( keys %$h_dejavu ) {
				$nb_project++;
				$st_project = $projName . ":";
				$st_project .= $h_dejavu->{$projName}->{patients};
				$he += $h_dejavu->{$projName}->{he};
				$ho += $h_dejavu->{$projName}->{ho};
			}
			$hash->{$chr_h_id}->{$var_id}->{'dejavu'} = "";
			$nb_patient = $he + $ho;
			if ( $nb_project > 0 ) {
				$hash->{$chr_h_id}->{$var_id}->{'dejavu'} = $nb_project . "/" . $nb_patient . " ho:$ho,he:$he";
			}
			$hash->{$chr_h_id}->{$var_id}->{'chr'}      = $chr->ucsc_name();
			$hash->{$chr_h_id}->{$var_id}->{'start'} = $var->start();
			$hash->{$chr_h_id}->{$var_id}->{'end'} = $var->end();
			
			
			$hash->{$chr_h_id}->{$var_id}->{'allele'} = '-';
			$hash->{$chr_h_id}->{$var_id}->{'hgmd_class'} = 'N.A.';
			if ( $self->can_use_hgmd() ) {
				$hash->{$chr_h_id}->{$var_id}->{'hgmd_class'} = '-';
				$hash->{$chr_h_id}->{$var_id}->{'hgmd_class'} = $var->hgmd_class() if $var->hgmd_class();
			}
			$hash->{$chr_h_id}->{$var_id}->{'cadd_score'} = '-';
			$hash->{$chr_h_id}->{$var_id}->{'cadd_score'} = $var->cadd_score() if defined $var->cadd_score();
			$hash->{$chr_h_id}->{$var_id}->{'cadd_score'} = '-' if ( $hash->{$chr_h_id}->{$var_id}->{'cadd_score'} eq '-1' );
			$hash->{$chr_h_id}->{$var_id}->{'ncboost_score'} = '-';
			$hash->{$chr_h_id}->{$var_id}->{'ncboost_score'} = $var->ncboost_score() if defined $var->ncboost_score();
			
			
			if ( $var->cosmic() ) {
				my @lTmpCosmic = split( ':', $var->cosmic() );
				$hash->{$chr_h_id}->{$var_id}->{'cosmic'} = $lTmpCosmic[0];
			}
			else { $hash->{$chr_h_id}->{$var_id}->{'cosmic'} = '-'; }
			$hash->{$chr_h_id}->{$var_id}->{'clinvar'} = '-';
			$hash->{$chr_h_id}->{$var_id}->{'clinvar'} = $var->text_clinvar() if ( $var->text_clinvar() ne '-5' );
			my $min_pop = '-';
			if ( $var->min_pop_name() ) {
				$min_pop = $var->min_pop_name();
				$min_pop .= ':' . sprintf("%.3f", $var->min_pop_freq()) if ( defined $var->min_pop_freq() );
			}
			my $max_pop = '-';
			if ( $var->max_pop_name() ) {
				$max_pop = $var->max_pop_name();
				$max_pop .= ':' . sprintf("%.3f", $var->max_pop_freq()) if ( defined $var->max_pop_freq() );
			}
			if (defined $var->frequency()) {
				if ($var->frequency()) { $hash->{$chr_h_id}->{$var_id}->{'freq'} = $var->percent(); }
				else { $hash->{$chr_h_id}->{$var_id}->{'freq'} = 0; }
			}
			else {
				$hash->{$chr_h_id}->{$var_id}->{'freq'} = '-';
			}
			$hash->{$chr_h_id}->{$var_id}->{'gnomad ac'} = '-';
			$hash->{$chr_h_id}->{$var_id}->{'gnomad ho'} = '-';
			$hash->{$chr_h_id}->{$var_id}->{'gnomad an'} = '-';
			$hash->{$chr_h_id}->{$var_id}->{'gnomad ac'} = $var->getGnomadAC() if ( $var->getGnomadAC() );
			$hash->{$chr_h_id}->{$var_id}->{'gnomad ho'} = $var->getGnomadHO() if ( $var->getGnomadHO() );
			$hash->{$chr_h_id}->{$var_id}->{'gnomad an'} = $var->getGnomadAN() if ( $var->getGnomadAN() );
			$hash->{$chr_h_id}->{$var_id}->{'min_pop_freq'} = $min_pop;
			$hash->{$chr_h_id}->{$var_id}->{'max_pop_freq'} = $max_pop;

			my $seq = $chr->sequence( ( $var->start() - 21 ), ( $var->end() - 1 ) );
			$seq .= '[' . $var->ref_allele() . '/' . $var->var_allele() . ']';
			$seq .= $chr->sequence( ( $var->start() + 1 ), ( $var->end() + 21 ) );
			$hash->{$chr_h_id}->{$var_id}->{'sequence'} = $seq;
			
			my $h_genes = $var->get_genes_transcripts_details_dup_del();
			foreach my $gene_id (keys %{$h_genes}) {
				my $gene = $project->newGene($gene_id);
				$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'external_name'} = $gene->external_name();
				$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'description'} = $gene->description();
				$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'phenotypes'} = $gene->phenotypes();
				foreach my $tr_id (keys %{$h_genes->{$gene_id}}) {
					my $t = $project->newTranscript($tr_id);
					$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'external_name'} = $t->external_name();
					$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'appris'} = $t->appris_type();
					$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'max_splice_ai'} = '-';
					
					foreach my $exon_intron_id (keys %{$h_genes->{$gene_id}->{$tr_id}->{'exons_introns'}}) {
						my $h_exon_intron = $h_genes->{$gene_id}->{$tr_id}->{'exons_introns'}->{$exon_intron_id};
						$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'exons_introns'}->{$exon_intron_id}->{'locus'} = $h_exon_intron->{'locus'};
						$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'exons_introns'}->{$exon_intron_id}->{'overlap'} = $h_exon_intron->{'overlap'};
					}
				}
			}
			$self->{hash_cnvs_global}->{$chr_h_id}->{$var_id} = $hash->{$chr_h_id}->{$var_id};
		}
		if ($list_patients) {
			foreach my $patient (@$list_patients) {
				$self->{hash_cnvs_global}->{$chr_h_id}->{$var_id}->{'patients'}->{$patient->name()}->{'cnv_confidence'} = '.';
				next unless exists ($var->sequencing_infos->{$patient->id()});
				foreach my $caller (keys %{$var->sequencing_infos->{$patient->id()}->{'method_calling'}}) {
					my $h_caller = $var->sequencing_infos->{$patient->id()}->{'method_calling'}->{$caller};
					$self->{hash_cnvs_global}->{$chr_h_id}->{$var_id}->{'patients'}->{$patient->name()}->{method_calling}->{$caller}->{'nb_all_ref'} = $h_caller->{'nb_all_ref'};
					$self->{hash_cnvs_global}->{$chr_h_id}->{$var_id}->{'patients'}->{$patient->name()}->{method_calling}->{$caller}->{'nb_all_mut'} = $h_caller->{'nb_all_mut'};
					$self->{hash_cnvs_global}->{$chr_h_id}->{$var_id}->{'patients'}->{$patient->name()}->{method_calling}->{$caller}->{'score'} = '-';
					$self->{hash_cnvs_global}->{$chr_h_id}->{$var_id}->{'patients'}->{$patient->name()}->{method_calling}->{$caller}->{'score'} = $h_caller->{'score'} if (exists $h_caller->{'score'});
				}
				$self->{hash_cnvs_global}->{$chr_h_id}->{$var_id}->{'patients'}->{$patient->name()}->{'he_ho'} = 'he' if $var->sequencing_infos->{$patient->id()}->{he} eq '1';
				$self->{hash_cnvs_global}->{$chr_h_id}->{$var_id}->{'patients'}->{$patient->name()}->{'he_ho'} = 'ho' if $var->sequencing_infos->{$patient->id()}->{ho} eq '1';
				#$self->{hash_cnvs_global}->{$chr_h_id}->{$var_id}->{'patients'}->{$patient->name()}->{'cnv_confidence'} = $var->cnv_confidence($patient);
			}
		}
	}
}

sub store_variants_infos {
	my ( $self, $list_var_objects, $project, $list_patients ) = @_;
	$project->buffer->dbh_deconnect();
	$project->buffer->dbh_reconnect();
	my @lCnv;
	foreach my $var (@$list_var_objects) {
		my $is_cnv;
		eval { $is_cnv = $var->isCnv(); };
		if($@) { $is_cnv = undef; }
		if ($is_cnv) {
			push(@lCnv, $var);
			next;
		}
		my $hash;
		$project->print_dot(25);
		
		
		my $chr      = $var->getChromosome();
		my $var_id   = $var->id();
		my $chr_h_id = $chr->id();
		$chr_h_id = '23' if ( $chr->id() eq 'X' );
		$chr_h_id = '24' if ( $chr->id() eq 'Y' );
		$chr_h_id = '25' if ( $chr->id() eq 'MT' );

		if ( exists $self->hash_variants_global->{$chr_h_id}->{$var_id} ) {
			$hash->{$chr_h_id}->{$var_id} = $self->hash_variants_global->{$chr_h_id}->{$var_id};
		}
		else {
			$hash->{$chr_h_id}->{$var_id}->{'var_id'} = $var_id;
			$hash->{$chr_h_id}->{$var_id}->{'rsname'} = $var->rs_name();
			$hash->{$chr_h_id}->{$var_id}->{'type'} = 'snp' if ( $var->isVariation() );
			$hash->{$chr_h_id}->{$var_id}->{'type'} = 'ins' if ( $var->isInsertion() );
			$hash->{$chr_h_id}->{$var_id}->{'type'} = 'del' if ( $var->isDeletion() );
			
			my $nb_project = $var->other_projects;
			my $nb_patient = $var->other_patients ;
			my $ho         = $var->other_patients_ho ;
			my $he         = $nb_patient - $ho;
			
			$hash->{$chr_h_id}->{$var_id}->{'dejavu'} = "";
			$nb_patient = $he + $ho;
			if ( $nb_project > 0 ) {
				$hash->{$chr_h_id}->{$var_id}->{'dejavu'} = $nb_project . "/" . $nb_patient . " ho:$ho,he:$he";
			}
			$hash->{$chr_h_id}->{$var_id}->{'chr'}      = $chr->ucsc_name();
			$hash->{$chr_h_id}->{$var_id}->{'position'} = $var->start();
			$hash->{$chr_h_id}->{$var_id}->{'allele'} = $var->ref_allele() . '/' . $var->var_allele();
			$hash->{$chr_h_id}->{$var_id}->{'hgmd_class'} = 'N.A.';
			if ( $self->can_use_hgmd() ) {
				$hash->{$chr_h_id}->{$var_id}->{'hgmd_class'} = '-';
				$hash->{$chr_h_id}->{$var_id}->{'hgmd_class'} = $var->hgmd_class() if $var->hgmd_class();
			}
			$hash->{$chr_h_id}->{$var_id}->{'cadd_score'} = '-';
			$hash->{$chr_h_id}->{$var_id}->{'cadd_score'} = $var->cadd_score() if defined $var->cadd_score();
			$hash->{$chr_h_id}->{$var_id}->{'cadd_score'} = '-' if ( $hash->{$chr_h_id}->{$var_id}->{'cadd_score'} eq '-1' );
			$hash->{$chr_h_id}->{$var_id}->{'ncboost_score'} = '-';
			$hash->{$chr_h_id}->{$var_id}->{'ncboost_score'} = $var->ncboost_score() if defined $var->ncboost_score();
			if ( $var->cosmic() ) {
				my @lTmpCosmic = split( ':', $var->cosmic() );
				$hash->{$chr_h_id}->{$var_id}->{'cosmic'} = $lTmpCosmic[0];
			}
			else { $hash->{$chr_h_id}->{$var_id}->{'cosmic'} = '-'; }
			$hash->{$chr_h_id}->{$var_id}->{'clinvar'} = '-';
			$hash->{$chr_h_id}->{$var_id}->{'clinvar'} = $var->text_clinvar() if ($var->text_clinvar() and  $var->text_clinvar() ne '-5' );
			my $min_pop = '-';
			if ( $var->min_pop_name() ) {
				$min_pop = $var->min_pop_name();
				$min_pop .= ':' . sprintf("%.3f", $var->min_pop_freq()) if ( defined $var->min_pop_freq() );
			}
			my $max_pop = '-';
			if ( $var->max_pop_name() ) {
				$max_pop = $var->max_pop_name();
				$max_pop .= ':' . sprintf("%.3f", $var->max_pop_freq()) if ( defined $var->max_pop_freq() );
			}
			if (defined $var->frequency()) {
				if ($var->frequency()) { $hash->{$chr_h_id}->{$var_id}->{'freq'} = $var->percent(); }
				else { $hash->{$chr_h_id}->{$var_id}->{'freq'} = 0; }
			}
			else {
				$hash->{$chr_h_id}->{$var_id}->{'freq'} = '-';
			}
			$hash->{$chr_h_id}->{$var_id}->{'gnomad ac'} = '-';
			$hash->{$chr_h_id}->{$var_id}->{'gnomad ho'} = '-';
			$hash->{$chr_h_id}->{$var_id}->{'gnomad an'} = '-';
			$hash->{$chr_h_id}->{$var_id}->{'gnomad ac'} = $var->getGnomadAC() if ( $var->getGnomadAC() );
			$hash->{$chr_h_id}->{$var_id}->{'gnomad ho'} = $var->getGnomadHO() if ( $var->getGnomadHO() );
			$hash->{$chr_h_id}->{$var_id}->{'gnomad an'} = $var->getGnomadAN() if ( $var->getGnomadAN() );
			$hash->{$chr_h_id}->{$var_id}->{'min_pop_freq'} = $min_pop;
			$hash->{$chr_h_id}->{$var_id}->{'max_pop_freq'} = $max_pop;

			my $seq = $chr->sequence( ( $var->start() - 21 ), ( $var->end() - 1 ) );
			$seq .= '[' . $var->ref_allele() . '/' . $var->var_allele() . ']';
			$seq .= $chr->sequence( ( $var->start() + 1 ), ( $var->end() + 21 ) );
			$hash->{$chr_h_id}->{$var_id}->{'sequence'} = $seq;
			my @lGenes = @{ $var->getGenes() };
			unless (@lGenes) {
				$hash->{$chr_h_id}->{$var_id}->{'consequence'} = 'intergenic';
			}
			foreach my $gene (@lGenes) {
				my $g_id = $gene->id();
				$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'external_name'} = $gene->external_name();
				$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'description'} = $gene->description();
				$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'phenotypes'} = $gene->phenotypes();
				eval { $var->annotation(); };
				if ($@) { next; }
				my ($max_value, $max_cat);
				$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'spliceAI'} = '-';
				my $h_score_spliceAI = $var->spliceAI_score($gene);
				my $splice_ai_txt = $var->text_max_spliceAI($gene);
				
				my @list_transcripts;
				if ($self->only_main_transcripts()) {
					foreach my $tr (@{$gene->getTranscripts()}) { push(@list_transcripts, $tr) if $tr->isMane(); }
					if (scalar(@list_transcripts) == 0) {
						@list_transcripts = @{$gene->getTranscripts()};
					} 
				}
				else {
					@list_transcripts = @{$gene->getTranscripts()};
				}
				
				foreach my $t (@list_transcripts) {
					
					next unless ( exists $var->annotation->{ $t->id() } );
					my @ok;
					my $annot_trans;
					eval { $annot_trans = $var->variationType($t) };
					if ($@) { $annot_trans = 'N.A'; }
					foreach my $cons ( split( ',', $annot_trans ) ) {
						my @lTmp = split( ';', $project->buffer->config->{ensembl_annotations}->{$cons} );
						push( @ok, $lTmp[1] );
					}
					next if ( scalar(@ok) == 0 );
					$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'max_splice_ai'} = $splice_ai_txt;
					if ($var->promoterAI_score($t)) {
						$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'promoter_ai'} = $var->promoterAI_score($t);
					}
					else {
						$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'promoter_ai'} = '-';
					}
					$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'alphamissense'} = '-';
					$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'alphamissense'} = $var->_alphamissenseTranscripts($t) if $var->isVariation();
					$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'external_name'} = $t->external_name();
					$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'appris'} = $t->appris_type();
					$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'consequence'} = join( ',', @ok );
					foreach my $cons (@ok) {
						$hash->{$chr_h_id}->{$var_id}->{'consequences'}->{$cons}= undef;
					}
					$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'cdna_position'} = $t->translate_position( $var->start() );
					$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'exon'} = $t->findExonNumber( $var->start() );
					if ( $hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'exon'} == -1 ) {
						$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'exon'} = $t->findNearestExon( $var->start(), $var->end() );
					}
					$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'nomenclature'} = $var->getNomenclature($t);
					my $prot = $t->getProtein();
					if ($prot) {
						$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'protein'} = $prot->id();
						$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'prot_nomenclature'} = $var->protein_nomenclature($prot);
						my $cds_pos = $var->getOrfPosition($prot);
						$cds_pos = '-' if (not $cds_pos or $cds_pos eq '.');
						$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'cds_position'} = $cds_pos;
						$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'protein_position'} = $var->getProteinPosition($prot);
						if ($var->isCnv() or $var->isLarge()) {
							$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'polyphen_score'} = '-';
							$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'sift_score'} = '-';
						}
						else {
							$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'polyphen_score'} = $var->polyphenScore($t);
							$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'sift_score'} = $var->siftScore($t);
						}
						my $protAA = $var->getProteinAA($prot);
						my $chanAA = $var->changeAA($prot);
						if ( $protAA and $chanAA ) {
							$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'aa'} = $protAA . '/' . $chanAA;
						}
					}
					else {
						$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'protein'} = '-';
						$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'cds_position'} = '-';
						$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'prot_nomenclature'} = '-';
						$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'protein_position'} = '-';
						$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'polyphen_score'} = '-';
						$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'sift_score'} = '-';
						$hash->{$chr_h_id}->{$var_id}->{'genes'}->{ $gene->id() }->{'transcripts'}->{ $t->id() }->{'aa'} = '-';
					}
				}
			}
			$self->{hash_variants_global}->{$chr_h_id}->{$var_id} = $hash->{$chr_h_id}->{$var_id};
		}
		if ($list_patients) {
			foreach my $patient (@$list_patients) {
				if (exists ($var->sequencing_infos->{$patient->id()})) {
					$self->{hash_variants_global}->{$chr_h_id}->{$var_id}->{'patients'}->{$patient->name()}->{method_calling}->{'max'}->{'nb_all_ref'} = $var->getNbAlleleRef($patient);
					$self->{hash_variants_global}->{$chr_h_id}->{$var_id}->{'patients'}->{$patient->name()}->{method_calling}->{'max'}->{'nb_all_mut'} = $var->getNbAlleleAlt($patient);
					$self->{hash_variants_global}->{$chr_h_id}->{$var_id}->{'patients'}->{$patient->name()}->{'he_ho'} = 'he' if $var->isHeterozygote($patient);
					$self->{hash_variants_global}->{$chr_h_id}->{$var_id}->{'patients'}->{$patient->name()}->{'he_ho'} = 'ho' if $var->isHomozygote($patient);
					$self->{hash_variants_global}->{$chr_h_id}->{$var_id}->{'patients'}->{$patient->name()}->{'model'} ='';
					if (not $var->isMei()) {
						$self->{hash_variants_global}->{$chr_h_id}->{$var_id}->{'patients'}->{$patient->name()}->{'model'} = $var->getTransmissionModelType($patient->getFamily(), $patient);
					}
				}
				else {
					$self->{hash_variants_global}->{$chr_h_id}->{$var_id}->{'patients'}->{$patient->name()}->{'he_ho'} = '';
					$self->{hash_variants_global}->{$chr_h_id}->{$var_id}->{'patients'}->{$patient->name()}->{'he_ho'} = '';
					$self->{hash_variants_global}->{$chr_h_id}->{$var_id}->{'patients'}->{$patient->name()}->{'model'} = '';
				}
			}
		}
	}
	if (@lCnv and scalar(@lCnv) > 0) {
		$self->store_cnvs_infos(\@lCnv, $project, $list_patients );
	}
}

sub store_specific_infos {
	my ( $self, $cat_name, $obj_to_store ) = @_;
	$self->{hash_specific_infos}->{datas}->{$cat_name} = $obj_to_store;
}

sub get_specific_infos_stored {
	my ( $self, $cat_name ) = @_;
	confess("\n\nERROR: $cat_name not stored with xls_export. Die.\n\n") unless ( exists $self->{hash_specific_infos}->{datas}->{$cat_name} );
	return $self->{hash_specific_infos}->{datas}->{$cat_name};
}

sub prepare_generic_datas_junctions {
	my $self = shift;
	my @lDatas;
	my ($h_patients_found);
	foreach my $chr_id ( sort { $a <=> $b } keys %{ $self->hash_junctions_global() } ) {
		foreach my $junc_id (sort keys %{$self->hash_junctions_global->{$chr_id}}) {
			my $h_junc = $self->hash_junctions_global->{$chr_id}->{$junc_id};
			my $h;
			$h->{'junction'}		= $junc_id;
			$h->{'type'}     	 	= $h_junc->{global}->{'type'};
			$h->{'dejavu'}			= $h_junc->{global}->{'dejavu'};
			$h->{'dejavu_ratio_10'}	= $h_junc->{global}->{'dejavu_ratio_10'};
			$h->{'dejavu_ratio_20'} = $h_junc->{global}->{'dejavu_ratio_20'};
			$h->{'chr'}				= $h_junc->{global}->{'chr'};
			$h->{'start'}			= $h_junc->{global}->{'start'};
			$h->{'end'}				= $h_junc->{global}->{'end'};
			foreach my $gene_name (sort keys %{$h_junc->{annotation}}) {
				my $h_junc_gene = $h_junc->{annotation}->{$gene_name};
				if (exists $h_junc_gene->{'gene_name'} and $h_junc_gene->{'gene_name'}) {
					$h->{'gene'} = $h_junc_gene->{'gene_name'}.' ['.$h_junc_gene->{'ensg'}.']';
				}
				else { $h->{'gene'} = $h_junc_gene->{'ensg'}; }
				$h->{'description'} = $h_junc_gene->{'description'};
				$h->{'phenotypes'}	= $h_junc_gene->{'phenotypes'};
				if ($h_junc_gene->{'score'} < -100) { $h->{'score'} = -100; }
				else { $h->{'score'} = $h_junc_gene->{'score'}; }
				
				foreach my $tid (sort keys %{$h_junc_gene->{transcripts}}) {
					my $hjunc_tr			 = $h_junc_gene->{transcripts}->{$tid};
					$h->{'transcript'}		 = $tid;
					$h->{'transcript_xref'}	 = $hjunc_tr->{transcript_name};
					$h->{'appris_type'}		 = $hjunc_tr->{appris_type};
					$h->{'tr_start'}		 = $hjunc_tr->{start};
					$h->{'tr_end'}			 = $hjunc_tr->{end};
					
					foreach my $patient_name (sort keys %{$h_junc->{patients}}) {
						my $j_junc_pat = $h_junc->{patients}->{$patient_name};
						$h->{'family'}		 = $j_junc_pat->{fam_name};
						$h->{'patient'}		 = $patient_name;
						$h->{'sex'}			 = $j_junc_pat->{sex};
						$h->{'status'}		 = $j_junc_pat->{status};
						$h->{'dp'}			 = $j_junc_pat->{dp};
						$h->{'nb_new'}		 = $j_junc_pat->{nb_new};
						$h->{'nb_canonique'} = $j_junc_pat->{nb_canonique};
						$h->{'ratio'}		 = $j_junc_pat->{ratio};
						push( @lDatas, $h );
					}
				}
			}
		}
	}
	return ( \@lDatas );
}

sub prepare_generic_datas_cnvs {
	my $self = shift;
	my @lDatas;
	my ($h_patients_found);
	foreach my $chr_id ( sort { $a <=> $b } keys %{ $self->hash_cnvs_global() } ) {
		foreach my $var_id ( sort keys %{ $self->hash_cnvs_global->{$chr_id} } ) {
			my $h_var = $self->hash_cnvs_global->{$chr_id}->{$var_id};
			my $h;
			$h->{'variation'}    = $h_var->{'var_id'};
			$h->{'rsname'}   	 = $h_var->{'rsname'};
			$h->{'type'}         = $h_var->{'type'};
			$h->{'dejavu'}       = $h_var->{'dejavu'};
			$h->{'chr'}          = $h_var->{'chr'};
			$h->{'start'}     	 = $h_var->{'start'};
			$h->{'end'}			 = $h_var->{'end'};
			$h->{'allele'}       = $h_var->{'allele'};
			$h->{'sequence'}     = $h_var->{'sequence'};
			$h->{'freq (%)'}     = $h_var->{'freq'};
			$h->{'gnomad ac'}    = $h_var->{'gnomad ac'};
			$h->{'gnomad ho'}    = $h_var->{'gnomad ho'};
			$h->{'gnomad an'}    = $h_var->{'gnomad an'};
			$h->{'cadd'}         = $h_var->{'cadd_score'};
			$h->{'ncboost'}      = $h_var->{'ncboost_score'};
			$h->{'cosmic'}       = $h_var->{'cosmic'};
			$h->{'clinvar'}      = $h_var->{'clinvar'};
			$h->{'min_pop_freq'} = $h_var->{'min_pop_freq'};
			$h->{'max_pop_freq'} = $h_var->{'max_pop_freq'};
			$h->{'hgmd_class'}   = $h_var->{'hgmd_class'};
			my $h_conf;
			my ($nb_pat_he, $nb_pat_ho);
			foreach my $patient_name (keys %{$h_var->{'patients'}}) {
				$h_patients_found->{$patient_name} = undef;
				if (exists $h_var->{'patients'}->{$patient_name}) {
					my $h_pat = $h_var->{'patients'}->{$patient_name};
					my $he_ho = $h_pat->{'he_ho'};
					my $cnv_confidence = $h_pat->{'cnv_confidence'};
					$h_conf->{$cnv_confidence} = undef;
					$h->{lc('patient_'.$patient_name)} = $cnv_confidence;
					$nb_pat_he++ if ($he_ho and $he_ho eq 'he');
					$nb_pat_ho++ if ($he_ho and $he_ho eq 'ho');
					$h->{lc('patient_'.$patient_name)} .= ' '.$he_ho if ($he_ho);
					foreach my $caller (keys %{$h_pat->{'method_calling'}}) {
						my $h_caller = $h_pat->{'method_calling'}->{$caller};
						my $nb_all_ref = $h_caller->{'nb_all_ref'};
						my $nb_all_mut = $h_caller->{'nb_all_mut'};
						$h->{lc('patient_'.$patient_name)} .= ' '.substr $caller, 0, 3;
						$h->{lc('patient_'.$patient_name)} .= ':'.$nb_all_ref.'/'.$nb_all_mut;
						if (exists $h_caller->{score}) {
							$h->{lc('patient_'.$patient_name)} .= ';score:'.$h_caller->{score};
						}
					}
				}
				else { $h->{lc('patient_'.$patient_name)} .= '-'; }
			}
			if (exists $h_conf->{'high'}) { $h->{'cnv_confidence'} = 'high'; }
			elsif (exists $h_conf->{'medium'}) { $h->{'cnv_confidence'} = 'medium'; }
			elsif (exists $h_conf->{'low'}) { $h->{'cnv_confidence'} = 'low'; }
			else { $h->{'cnv_confidence'} = '-'; }
			if ($nb_pat_he or $nb_pat_ho) {
				$h->{'he'} = 0;
				$h->{'ho'} = 0;
				$h->{'he'} = int($nb_pat_he) if ($nb_pat_he);
				$h->{'ho'} = int($nb_pat_ho) if ($nb_pat_ho);
			}
			if (exists $h_var->{genes}) {
				foreach my $gene_id ( sort %{ $h_var->{genes} } ) {
					my $h_gene      = $h_var->{genes}->{$gene_id};
					my $gene_name   = $h_gene->{'external_name'};
					foreach my $tr_id ( sort keys %{ $h_gene->{transcripts} } ) {
						my $h_tr = $h_gene->{'transcripts'}->{$tr_id};
						foreach my $exon_intron_id ( sort keys %{ $h_tr->{exons_introns} } ) {
							my $h_exon_intron = $h_tr->{'exons_introns'}->{$exon_intron_id};
							my $h2 = dclone(\%$h);
							$h2->{'gene'}         = $gene_id;
							$h2->{'gene'} .= '[' . $gene_name . ']' if ($gene_name);
							$h2->{'transcript'}      = $tr_id;
							$h2->{'transcript_xref'} = $h_tr->{'external_name'};
							$h2->{'description'}     = $h_gene->{'description'};
							$h2->{'phenotypes'}      = $h_gene->{'phenotypes'};
							$h2->{'exon_intron'}     = $exon_intron_id;
							$h2->{'locus'}           = $h_exon_intron->{'locus'};
							$h2->{'overlap'}         = $h_exon_intron->{'overlap'};
							push( @lDatas, $h2 );
						}
					}
				}
			}
			else {
				my $h2 = dclone(\%$h);
				$h2->{'gene'} = '-';
				push( @lDatas, $h2 );
			}
		}
	}
	if ($h_patients_found) {
		$self->list_generic_header_cnvs();
		push (@{$self->{list_generic_header_cnvs}}, 'Cnv_confidence');
		push (@{$self->{list_generic_header_cnvs}}, 'He');
		push (@{$self->{list_generic_header_cnvs}}, 'Ho');
		foreach my $patient_name (sort keys %$h_patients_found) {
			push (@{$self->{list_generic_header_cnvs}}, 'Patient_'.$patient_name);
		}
	}
	return ( \@lDatas );
}

sub prepare_generic_datas_variants {
	my $self = shift;
	my @lDatas;
	my ($h_patients_found);
	foreach my $chr_id ( sort { $a <=> $b } keys %{ $self->hash_variants_global() } ) {
		foreach my $var_id ( sort keys %{ $self->hash_variants_global->{$chr_id} } ) {
			my $h_var = $self->hash_variants_global->{$chr_id}->{$var_id};
			my $h;
			$h->{'variation'}    = $h_var->{'var_id'};
			if ($h->{'variation'} =~ /M/) {
				my @lCol = split('_', $h->{'variation'});
				$h->{'variation'} = $lCol[1].$lCol[3];
			}
			$h->{'rsname'}  	 = $h_var->{'rsname'};
			$h->{'type'}         = $h_var->{'type'};
			$h->{'dejavu'}       = $h_var->{'dejavu'};
			$h->{'chr'}          = $h_var->{'chr'};
			$h->{'position'}     = $h_var->{'position'};
			$h->{'allele'}       = $h_var->{'allele'};
			$h->{'sequence'}     = $h_var->{'sequence'};
			$h->{'freq (%)'}     = $h_var->{'freq'};
			$h->{'gnomad ac'}    = $h_var->{'gnomad ac'};
			$h->{'gnomad ho'}    = $h_var->{'gnomad ho'};
			$h->{'gnomad an'}    = $h_var->{'gnomad an'};
			$h->{'cadd'}         = $h_var->{'cadd_score'};
			$h->{'ncboost'}      = $h_var->{'ncboost_score'};
			$h->{'cosmic'}       = $h_var->{'cosmic'};
			$h->{'clinvar'}      = $h_var->{'clinvar'};
			$h->{'min_pop_freq'} = $h_var->{'min_pop_freq'};
			$h->{'max_pop_freq'} = $h_var->{'max_pop_freq'};
			$h->{'hgmd_class'}   = $h_var->{'hgmd_class'};
			my ($nb_pat_he, $nb_pat_ho);
			foreach my $patient_name (keys %{$h_var->{'patients'}}) {
				$h_patients_found->{$patient_name} = undef;
				my $h_pat = $h_var->{'patients'}->{$patient_name};
#				warn Dumper $h_pat; die;
				my $he_ho = $h_pat->{'he_ho'};
				$nb_pat_he++ if ($he_ho eq 'he');
				$nb_pat_ho++ if ($he_ho eq 'ho');
				my $model;
				$model = $h_pat->{'model'} if (exists $h_pat->{'model'});
				$model = '' if ($model eq '?');
				$h->{lc('patient_'.$patient_name)} = '';
				$h->{lc('patient_'.$patient_name)} .= $he_ho if ($he_ho);
				foreach my $caller (keys %{$h_pat->{'method_calling'}}) {
					my $h_caller = $h_pat->{'method_calling'}->{$caller};
					my $nb_all_ref = $h_caller->{'nb_all_ref'};
					my $nb_all_mut = $h_caller->{'nb_all_mut'};
					$h->{lc('patient_'.$patient_name)} .= ' '.substr $caller, 0, 3;
					$h->{lc('patient_'.$patient_name)} .= ':'.$nb_all_ref.'/'.$nb_all_mut;
				}
				$h->{lc('patient_'.$patient_name)} .= ' '.$model;
			}
			if ($nb_pat_he or $nb_pat_ho) {
				$h->{'he'} = 0;
				$h->{'ho'} = 0;
				$h->{'he'} = int($nb_pat_he) if ($nb_pat_he);
				$h->{'ho'} = int($nb_pat_ho) if ($nb_pat_ho);
			}
			
#			if (exists $h_var->{genes}) {
				foreach my $gene_id ( sort %{ $h_var->{genes} } ) {
					my $h_gene      = $h_var->{genes}->{$gene_id};
					my $gene_name   = $h_gene->{'external_name'};
					my $description = $h_gene->{'description'};
					foreach my $tr_id ( sort keys %{ $h_gene->{transcripts} } ) {
						my $h_tr = $h_gene->{'transcripts'}->{$tr_id};
						my $h2 = dclone(\%$h);
						$h2->{'gene'}         = $gene_id;
						$h2->{'gene'} .= '[' . $gene_name . ']' if ($gene_name);
						$h2->{'consequence'}     = $h_tr->{'consequence'};
						$h2->{'transcript'}      = $tr_id;
						$h2->{'transcript_xref'} = $h_tr->{'external_name'};
						$h2->{'description'}     = $h_gene->{'description'};
						$h2->{'phenotypes'}      = $h_gene->{'phenotypes'};
						$h2->{'exon'}            = $h_tr->{'exon'};
						$h2->{'cdna_pos'}        = $h_tr->{'cdna_position'};
						$h2->{'cds_Pos'}         = $h_tr->{'cds_position'};
						$h2->{'protein'}         = $h_tr->{'protein'};
						$h2->{'nomenclature'}    = $h_tr->{'nomenclature'};
						$h2->{'prot_nomenclature'}= $h_tr->{'prot_nomenclature'};
						$h2->{'polyphen_score'}  = $h_tr->{'polyphen_score'};
						$h2->{'sift_score'}      = $h_tr->{'sift_score'};
						$h2->{'max_splice_ai'}   = $h_tr->{'max_splice_ai'};
						$h2->{'alphamissense'}   = $h_tr->{'alphamissense'};
						$h2->{'promoter_ai'}     = $h_tr->{'promoter_ai'};
						$h2->{'aa'}              = $h_tr->{'aa'};
						$h2->{'appris'}          = $h_tr->{'appris'};
						push( @lDatas, $h2 );
					}
				}
#			}
#			else {
#				my $h2 = dclone(\%$h);
#				$h2->{'consequence'}   = $h_var->{'consequence'};
#				push( @lDatas, $h2 );
#			}
		}
	}
	if ($h_patients_found) {
		$self->list_generic_header();
		push (@{$self->{list_generic_header}}, 'He');
		push (@{$self->{list_generic_header}}, 'Ho');
		foreach my $patient_name (sort keys %$h_patients_found) {
			push (@{$self->{list_generic_header}}, 'Patient_'.$patient_name);
		}
	}
	return ( \@lDatas );
}

sub add_infos_patients_variants {
	my ($self, $h_var) = @_;
	my ($h, $h_patients_found);
	my ($nb_pat_he, $nb_pat_ho);
	foreach my $patient_name (keys %{$h_var->{'patients'}}) {
		$h_patients_found->{$patient_name} = undef;
		my $h_pat = $h_var->{'patients'}->{$patient_name};
		my $he_ho = $h_pat->{'he_ho'};
		$nb_pat_he++ if ($he_ho eq 'he');
		$nb_pat_ho++ if ($he_ho eq 'ho');
		my $model;
		$model = $h_pat->{'model'} if (exists $h_pat->{'model'});
		$model = '' if ($model eq '?');
		$h->{lc('patient_'.$patient_name)} = '';
		$h->{lc('patient_'.$patient_name)} .= $he_ho if ($he_ho);
		foreach my $caller (keys %{$h_pat->{'method_calling'}}) {
			my $h_caller = $h_pat->{'method_calling'}->{$caller};
			my $nb_all_ref = $h_caller->{'nb_all_ref'};
			my $nb_all_mut = $h_caller->{'nb_all_mut'};
			$h->{lc('patient_'.$patient_name)} .= ' '.substr $caller, 0, 3;
			$h->{lc('patient_'.$patient_name)} .= ':'.$nb_all_ref.'/'.$nb_all_mut;
		}
		$h->{lc('patient_'.$patient_name)} .= ' '.$model;
	}
	if ($nb_pat_he or $nb_pat_ho) {
		$h->{'he'} = 0;
		$h->{'ho'} = 0;
		$h->{'he'} = int($nb_pat_he) if ($nb_pat_he);
		$h->{'ho'} = int($nb_pat_ho) if ($nb_pat_ho);
	}
	return ($h_patients_found, $h);
}

sub save {
	my ($self) = shift;
	my $session = new session_export();
	$session->save('xls_title', $self->title_page());
	$session->save_compress('hash_variants_global', $self->hash_variants_global());
	$session->save_compress('hash_cnvs_global', $self->hash_cnvs_global());
	$session->save_compress('hash_junctions_global', $self->hash_junctions_global());
	$session->save_compress('hash_specific_infos', $self->hash_specific_infos());
	return $session->session_id();
}

sub load {
	my ( $self, $sid ) = @_;
	my $session = new session_export();
	$session->load_session( $sid );
	$self->{title_page} = $session->load('xls_title');
	$self->{hash_variants_global} = $session->load_compress('hash_variants_global');
	$self->{hash_cnvs_global} = $session->load_compress('hash_cnvs_global');
	$self->{hash_junctions_global} = $session->load_compress('hash_junctions_global');
	$self->{hash_specific_infos} = $session->load_compress('hash_specific_infos');
}

1;
