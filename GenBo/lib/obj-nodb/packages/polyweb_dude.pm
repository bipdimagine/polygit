package polyweb_dude;
use strict;
use Moose;

use GD;
use Data::Printer;
use Data::Dumper;
use List::Util qw( min sum);
use Statistics::Descriptive;
use List::MoreUtils qw(indexes);
use URI;

has transcript => (
	is       => 'ro',
	required => 1,
);

has patients => (
	is       => 'ro',
	required => 1,
);

has selected_patients => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		return [];
	}
);

has ordered_patients => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my @ps   = sort {
			(       $a->getFamily->name
				  . $a->name cmp $b->getFamily->name
				  . $b->name )
		} @{ $self->patients };

		return \@ps;
	}
);

has hash_index_patients => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $hash;
		my $index = 0;
		foreach my $p ( @{ $self->ordered_patients } ) {
			$hash->{ $p->id } = $index;
			$index++;
		}
		return $hash;
	}
);

has index_selected_patients => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my %res;
		foreach my $p ( @{ $self->selected_patients } ) {
			confess() unless exists $self->hash_index_patients->{ $p->id };
			my $index = $self->hash_index_patients->{ $p->id };
			$res{$index}++;

		}
		return \%res;
	}
);

has index_selected_ill_patients => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my %res;
		foreach my $index ( keys %{ $self->index_selected_patients } ) {
			my $p = $self->ordered_patients->[$index];
			$res{$index}++ if $p->isIll();
		}
		return \%res;
	}
);

sub return_index_selected_patients {
	my ($self) = @_;
}

has limit => (
	is       => 'ro',
	required => 1,
);

has levels_matrix => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		$self->init_matrices();
		return $self->{levels_matrix};
	}
);

has scores_matrix => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		$self->init_matrices();
		return $self->{scores_matrix};
	}
);

has score_smooth_matrix => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		$self->init_matrices();
		return $self->{score_smooth_matrix};
	}
);
has score_smooth_expo_matrix => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		$self->init_matrices();
		return $self->{score_smooth_expo_matrix};
	}
);
has error => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		$self->init_matrices();
		return $self->{error};
	}
);

has cgi => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		return new CGI;
	}
);

has grey_hex => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		return [ "#cccccc", "#98B4D4", "#FF0000" ];
	}
);

has grey_rgb => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		return [ [ 200, 200, 200 ], [ 152, 180, 212 ], [ 127, 255, 0 ] ];
	}
);

has red_hex => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		return [ "#ff9900", "#ff6666", "#FF0000" ];
	}
);

has red_rgb => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		return [
			[ 155, 15, 10 ],

			#C0392B
			#rgb(231, 76, 60)
			#rgb(229, 19, 1)
			[ 229, 19,  1 ],
			[ 195, 85,  10 ],
			[ 255, 120, 20 ]
		];
	}
);

has red_rgb_selected => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		return [
			[ 255, 165, 0 ],
			[ 255, 127, 80 ],
			[ 255, 0,   0 ],
			[ 255, 0,   131 ]
		];
	}
);

has green_hex => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		return [ "#99cc33", "#99cc00", "#27AE60" ];
	}
);

has green_rgb => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		return [ [ 60, 128, 69 ], [ 160, 218, 169 ], [ 46, 204, 113 ] ];
	}
);

has blue_rgb => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		return [
			[ 50, 50, 145 ],

			#rgb(52, 86, 139)
			#rgb(62, 101, 255)
			[ 62, 101, 255 ],

			#[ 20, 109, 255 ],
			[ 30, 30, 115 ],
			[ 10, 89, 200 ],

		];
	}
);
has blue_rgb_selected => (
	is      => 'rw',
	lazy    => 1,
	default => sub {

		#return  [[255, 0, 0],[255, 0, 0],[255, 1, 0]];
		return [ [ 127, 205, 205 ], [ 52, 152, 219 ], [ 121, 206, 255 ] ];
	}
);

has purple_rgb => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		return [ [ 107, 91, 149 ], [ 181, 1, 254 ] ];
	}
);

has purple_hex => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		return [ "#666699", "#cc00ff" ];
	}
);

has mask_exons => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		return {
			exon       => 1,
			intron     => 2,
			non_coding => 4,
			capture    => 8,
			duplicate  => 16,
			"5_utr"    => 32,
			"3_utr"    => 64,
		};
	}
);

sub between {
	my ( $test, $fom, $tom ) = @_;
	no warnings;
	$fom < $tom
	  ? $test >= $fom && $test <= $tom
	  : $test >= $tom && $test <= $fom;
}

has names => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		$self->types;
		return $self->{names};
	}
);

has types => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self    = shift;
		my $patient = $self->patients->[0];
		my $no      = $patient->getTranscriptsDude("r");
		my $matrix  = $no->get( $self->transcript->id );

		$self->{names} = $matrix->{names};
		return [];
	}
);

has types_index => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		$self->exon_indexes;
		return $self->{types_index};
	}
);

has exons_index => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $mask = $self->mask_exons->{exon};

		my @index = indexes { $_ & $mask } @{ $self->types };
		foreach my $ind (@index) {
			push( @{ $self->{types_index} }, $self->types->[$ind] );
		}
		return \@index;
	}
);

has names_index => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my @name;
		foreach my $ind ( @{ $self->exons_index } ) {
			push( @name, $self->names->[$ind] );
		}
		return \@name;
	}
);

has quality => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		return {};
	}
);

sub return_cnv_scores {
	my ( $self, $patient, $transcript ) = @_;
	my $no     = $patient->getTranscriptsDude("r");
	my $matrix = $no->get( $transcript->id );
	my $nb     = $matrix->{nb};
	return [ map { $_ / 100 } unpack( "w" . $nb, $matrix->{scores} ) ];
}

sub init_matrices {
	my ($self)     = @_;
	my $patients   = $self->patients;
	my $transcript = $self->transcript;
	my $limit      = $self->limit;

	my $data;
	my $nbp = 0;
	my $error;
	my $data_mean;
	my $data_names;
	my $type_matrix;
	$self->{scores_matrix}            = [];
	$self->{levels_matrix}            = [];
	$self->{score_smooth_expo_matrix} = [];
	$self->{score_smooth_matrix}      = [];
	my $names;
	my $first = 1;

	foreach my $patient ( @{ $self->ordered_patients } ) {
		my $no     = $patient->getTranscriptsDude("r");
		my $matrix = $no->get( $transcript->id );
		next unless $matrix;
		my $nb = $matrix->{nb};

		my $scores_matrix =
		  [ map { $_ / 100 } unpack( "w" . $nb, $matrix->{scores} ) ];
		my $levels_matrix =
		  [ map { $_ - 1 } unpack( "w" . $nb, $matrix->{level} ) ];
		my $score_smooth_expo_matrix =
		  [ map { $_ / 100 } unpack( "w" . $nb, $matrix->{smooth_expo} ) ];
		my $score_smooth_matrix =
		  [ map { $_ / 100 } unpack( "w" . $nb, $matrix->{smooth} ) ];

		$self->quality->{ $patient->name() }->{all} +=
		  ( scalar @$scores_matrix );
		my @ll;
		my $lstring = "";
		for ( my $ii = 0 ; $ii < @$scores_matrix ; $ii++ ) {

			# for (my $i=0;$i<@$m;$i++){
			$self->{scores_matrix}->[$ii]->[$nbp] = $scores_matrix->[$ii];
			$self->{levels_matrix}->[$ii]->[$nbp] = $levels_matrix->[$ii];
			$self->{score_smooth_expo_matrix}->[$ii]->[$nbp] =
			  $score_smooth_expo_matrix->[$ii];

			#if ($score_smooth_expo_matrix->[$ii] == -1) {
			my $l = "A";
			if ( $levels_matrix->[$ii] == -1 ) {
				$self->quality->{grey}++;
				$l = "G";
				$self->quality->{ $patient->name() }->{grey}++;
			}
			elsif ( $scores_matrix->[$ii] < 0.1 ) {
				$self->quality->{del_ho}++;
				$self->quality->{ $patient->name() }->{del_ho}++;
				$self->quality->{ $patient->name() }->{all}++;
				$l = "R";
			}

			#elsif ($score_smooth_expo_matrix->[$ii] <= 0.7) {
			#elsif ($scores_matrix->[$ii] <= 0.7) {
			elsif ( $levels_matrix->[$ii] == 1 or $scores_matrix->[$ii] <= 0.7 )
			{
				$self->quality->{del}++;
				$self->quality->{ $patient->name() }->{del}++;
				$l = "R";
			}

			#elsif ($score_smooth_expo_matrix->[$ii] >= 1.4) {
			#elsif ($scores_matrix->[$ii] >= 1.4) {
			if ( $levels_matrix->[$ii] == 2 or $scores_matrix->[$ii] >= 1.4 ) {
				$self->quality->{dup}++;
				$self->quality->{ $patient->name() }->{dup}++;
				$l = "B";
			}
			$lstring .= $l;
			push( @ll, $l );
			$self->{score_smooth_matrix}->[$ii]->[$nbp] =
			  $score_smooth_matrix->[$ii];

		}
		my $a = $lstring;
		$a =~ s/((.)\2+)/$2 . length($1)/ge;
		my $b = $lstring;
		$b =~ s/(.)\1+/$1/g;
		$self->quality->{noise_string_cigar}->{ $patient->name() } = $a;
		$self->quality->{noise_string_cigar_all} .= $a;
		my @ops;
		$self->quality->{noise_string_compact_all} .= $b;
		$self->quality->{noise_string_compact}->{ $patient->name() } = $b;
		$nbp++;
		$self->{nb} = $nb;
	}

}

sub cigar_quality {
	my ( $self, $patient ) = @_;
	my @pos;
	for my $op (
		grep defined,
		self->quality->{noise_string_cigar}->{ $patient->name() } =~
		/([AGRB]\d+)/g
	  )
	{
		my ( $type, $len ) = $op =~ /(\D*)(\d*)/a;
		push @pos, [ $len || 1, uc $type ];
	}
	return \@pos;
}

sub translate_rgb {
	my ( $self, $color ) = @_;
	return ( "rgb(" . join( ",", @$color ) . ")" );

}

sub html_table {
	my ($self) = @_;

	my $cgi        = $self->cgi;
	my $patients   = $self->patients;
	my $data       = $self->levels_matrix;
	my $data_score = $self->scores_matrix;
	warn Dumper $data_score;
	warn Dumper $data;
	my $tid = $self->transcript->id;
	my $out;
	my $blue = $self->blue_rgb;
	my $red  = $self->red_rgb;

	my $green = $self->green_rgb;

	$out .= $cgi->start_div(
		{
			style =>
"text-align: center;vertical-align:middle;font-size: 9px;font-family:Verdana;margin-bottom: 0px;height:600px;"
		}
	);
	$out .= $cgi->start_table(
		{
			class => "table table-striped table-bordered table-hover",
			style =>
"text-align: center;vertical-align:middle;font-size: 9px;font-family:  Verdana;margin-bottom: 0px;padding-top: 4px;padding-bottom: 4px;padding-left: 4px;padding-right: 4px;overflow-y:scroll; max-height:90%;"
		}
	);
	$out .= $cgi->start_Tr();
	my $gene_name = $self->transcript->getGene->external_name;
	my $tr_name   = $self->transcript->name;
	my $chr_name  = $self->transcript->getChromosome()->ucsc_name;
	my $strand    = "fwd";
	$strand = "rev" if $self->transcript->strand == -1;
	my $pos =
	  "[" . $self->transcript->start . "-" . $self->transcript->end . "]";
	my @pnames;
	my $click_gene = qq{load_graph_gene('$tr_name');};

	my $text =
qq{<center><div class="btn   btn-xs btn-primary " style="font-size:7px;min-width:30px,position:relative;"  onClick= "$click_gene">$gene_name <br> $tr_name  <br>$chr_name:$pos ($strand)</div></center>};
	$out .= $cgi->th(
		{
			style =>
"padding: 3px;position:sticky;z-index:10;top:0;vertical-align:middle;text-align:center;"
		},
		$text
	);

	foreach my $patient ( @{ $self->ordered_patients } ) {

		#my $index_color = 0;
		#$index_color = 1 if  ( exists $self->index_selected_patients->{$j} );
		my $name          = $patient->name();
		my $click_patient = qq{load_graph_transcript('$name','$tr_name');};
		my $scolor        = "";
		if ( $patient->sex == 2 ) {
			$scolor = "background-color:rgb(255, 110, 199);color:black";
		}
		my $icon = $patient->return_icon;
		push( @pnames, $name );
		my $fname = $patient->getFamily()->name();
		my $text =
qq{<div class="btn   btn-xs btn-info " style="font-size:7px;min-width:30px,position:relative;$scolor"  onClick= "$click_patient"> $fname <br> $icon $name </div>};
		$out .= $cgi->th(
			{
				style =>
"padding: 3px;position:sticky;z-index:10;top:0;vertical-align:middle;text-align:center;"
			},
			$text
		);
	}

	$out .= $cgi->end_Tr();

	my $capture = $self->transcript->getChromosome()->getIntSpanCapture();

	#warn $self->levels_matrix;
	#$self->transcript->getChromosome()->getPrimers();

	for ( my $i = 0 ; $i < @{ $self->levels_matrix } ; $i++ ) {
		my $min_line    = $self->levels_matrix->[$i];
		my $mean_line   = $self->scores_matrix->[$i];
		my $smooth_line = $self->score_smooth_expo_matrix->[$i];

		$out .= $cgi->start_Tr();
		my $primer = $self->transcript->getPrimer( $self->names->[$i] );
		warn $i;
		my $exons = $primer->getExons;

		my $ename = "-";
		if ($exons) {
			$ename = join( ";",
				map    { $_->name }
				  grep { $_->getTranscript->name eq $self->transcript->name }
				  @$exons );
		}
		my $start = $primer->start;
		my $end   = $primer->end;

		#my $s1 = $exon->getGenomicSpan()->intersection($capture);
		my $click_exon_patients =
		  qq{load_graph_one_exon_all_patients('$tr_name','$ename');};
		my $pid = $primer->id;
		my $text =
qq{<div class="btn   btn-xs btn-info " style="font-size:8px;min-width:30px,position:relative;" '  style="border : 1px" onClick= "$click_exon_patients">$ename [$start-$end] </div>	$pid};
		$out .= $cgi->td(
			{ style => "padding: 3px;position:sticky;z-index:9;left:0;" },
			$text );

		#my $line_levels = $data_levels->[$i];
		my $nbr = grep { $_ == 2 } @{$min_line};
		my $nbb = grep { $_ == 1 or $_ == 3 or $_ == 4 } @{$min_line};
		my $nbc   = scalar(@$min_line);
		my $nbp   = ( $nbr + $nbb ) / $nbc;
		my $nbpr  = ($nbr) / $nbc;
		my $nbpb  = ($nbb) / $nbc;
		my $nbmin = $nbpr;
		if ( $nbpr > 0 ) {
			$nbpr = $nbpr / ( $nbpr + $nbpb );
		}

		foreach ( my $j = 0 ; $j < @$min_line ; $j++ ) {
			my $name    = $pnames[$j];
			my $bgcolor = $green->[1];    # "rgb(22, 160, 133)";

			my $level  = $min_line->[$j];
			my $score  = $mean_line->[$j];
			my $smooth = $smooth_line->[$j];

			#@colors =(172, 223, 172);

			my $colors;
			if ( $level == -2 ) {
				$bgcolor = "#5A5A5A";

				#$color = $grey1;
			}
			elsif ( $level == -1 ) {
				$bgcolor = "#8C8C8C";
			}
			elsif ( $nbc > 10 && $nbp >= 0.5 && $nbpr >= 0.2 && $nbpr <= 0.8 ) {
				$bgcolor = "#DADADA";
			}
			elsif ( $level == 2 ) {

				$bgcolor = $self->translate_rgb( $blue->[1] );    #"#0247FD";
			}
			elsif ( $level == 3 ) {
				$bgcolor = $self->translate_rgb( $red->[1] );     #"#0247FD";

			}
			elsif ( $level == 1 ) {
				$bgcolor = $self->translate_rgb( $red->[1] );     #"#FF4040";
				if ( $score < 0.2 ) {
					$bgcolor = $self->translate_rgb( [ 0, 0, 0 ] );  #"#FF4040";
				}
			}
			elsif ( $level == 4 ) {
				$bgcolor = $self->translate_rgb( $red->[1] );
				;                                                    #"#FF4040";

			}
			elsif ( $level == 0 ) {

				#FFC9CF
				if ( $score <= 0.55 ) {
					$bgcolor = $self->translate_rgb( $red->[3] );

					#$bgcolor = $self->red_hex->[2];    #"#FE8A35";
				}
				elsif ( $score <= 0.6 ) {
					$bgcolor = "#F6B386";
					$bgcolor = $self->translate_rgb( $red->[3] );
				}
				elsif ( $score <= 0.7 ) {
					$bgcolor = "#FDDAB2";
					$bgcolor = $self->translate_rgb( $red->[3] );
				}

				elsif ( $score >= 1.5 ) {
					$bgcolor = $self->translate_rgb( $blue->[3] );
				}
				elsif ( $score >= 1.4 ) {
					$bgcolor = $self->translate_rgb( $blue->[3] );
				}

			}

			my $v1 = $mean_line->[$j];
			my $v2 = $min_line->[$j];

			#	my $v4 = $primer->sd($opatients[$i]);

			my $click =
			  qq{load_graph_exon('$name','$tid','$ename','$start','$end');};
			my $l2 =
qq{	<sup><span class="badge label-success" style ="font-size:09px;">$v2</span></sup>};
			my $l1 =
qq{	<span class="label label-xs label-success" style ="min-width: 55px !important;display: inline-block !important;font-size:09px;color:black;;background-color:$bgcolor;color:white" onClick= "$click">$v1  <span class="label label-xs label-primary" style ="font-size:08px;float:right">$v2</span></span>};

#my $l1 = qq{	<span class="label label-xs label-success" style ="min-width: 55px !important;display: inline-block !important;font-size:09px;color:black;;background-color:$bgcolor;color:white" onClick= "$click">$v1 <span class="label label-xs label-primary" style ="font-size:08px;float:right">$v2</span></span>};
			my $text =
qq{<div class="btn   btn-xs  " style="font-size:9px;min-width:50px,position:relative;" '  style="border : 1px">  $l1 $l2 </div>};

			$out .= $cgi->td( { style => "min-width:50px; padding:4px" }, $l1 );
		}
		$out .= $cgi->end_Tr();
	}
	$out .= $cgi->end_table();
	$out .= $cgi->end_div();
	return $out;
}

sub return_level_red {
	my ( $self, $min ) = @_;

	my $max3 = $self->limit - ( 0.1 * $self->limit );
	if ( between( $min, $max3, $self->limit ) ) {
		return 0;
	}
	my $max2 = $self->limit - ( 0.25 * $self->limit );
	if ( between( $min, $max2, $self->limit ) ) {
		return 1;
	}
	return 2;
}

sub return_level_green {
	my ( $self, $min ) = @_;
	my $limit = $self->limit;
	my $max4  = $limit + ( 0.1 * $limit );
	if ( between( $min, $limit, $max4 ) ) {
		return 0;
	}
	my $max5 = $limit + ( 0.25 * $limit );
	if ( between( $min, $limit, $max5 ) ) {
		return 1;
	}
	return 2;
}

sub control_ho {
	my ( $self, $patient, $debug ) = @_;
	my $limit       = $self->limit;
	my $patients    = $self->patients;
	my $nb_patient  = scalar(@$patients);
	my $data_levels = $self->levels_matrix;
	my @t           = keys %{ $self->index_selected_ill_patients };
	my $max         = 0;
	my $nb_samples  = scalar(@$patients);
	my @data;
	my $pindex = $self->hash_index_patients->{ $patient->id };
	confess() unless exists $self->hash_index_patients->{ $patient->id };

	for ( my $i = 0 ; $i < @$data_levels ; $i++ ) {
		my $ok;
		next if $self->scores_matrix()->[$i]->[$pindex] >= 0.15;
		if ( $nb_patient < 5 ) {
			$max += $self->control_line2( $patient, $i, $debug );
		}
		else {
			$max += $self->control_line( $self->scores_matrix()->[$i],
				$patient, $debug );
		}

		#my $primer = $self->transcript->getPrimer( $self->names->[$i] );
		#warn Dumper $primer->{cnv} if $debug;
		#die() if $debug;
	}
	return $max;

	#return $max;
}

sub control_del {
	my ( $self, $patient, $debug ) = @_;
	my $limit       = $self->limit;
	my $patients    = $self->ordered_patients;
	my $nb_patient  = scalar(@$patients);
	my $data_levels = $self->levels_matrix;
	my @t           = keys %{ $self->index_selected_patients };
	my $max         = 0;
	my $nb_samples  = scalar(@$patients);
	my @data;
	my $pindex = $self->hash_index_patients->{ $patient->id };

	confess() unless exists $self->hash_index_patients->{ $patient->id };
	for ( my $i = 0 ; $i < @$data_levels ; $i++ ) {
		my $ok;

		#foreach my $tt (keys %{$hi}) {
		next if $data_levels->[$i]->[$pindex] < 1;
		if ( $nb_patient < 5 ) {
			$max += $self->control_line2( $patient, $i, $debug );
		}
		else {
			warn "coucou " if $debug;
			$max += $self->control_line( $self->scores_matrix()->[$i],
				$patient, $debug );

		}

	}

	#warn $max;
	return $max;

	#return $max;
}

sub control_dup {
	my ( $self, $patient, $debug ) = @_;
	my $limit       = $self->limit;
	my $patients    = $self->ordered_patients;
	my $data_levels = $self->levels_matrix;
	my $nb_patient  = scalar(@$patients);
	my @t           = keys %{ $self->index_selected_patients };
	my $max         = 0;
	my $nb_samples  = scalar(@$patients);
	my @data;
	my $pindex = $self->hash_index_patients->{ $patient->id };

	for ( my $i = 0 ; $i < @$data_levels ; $i++ ) {
		next unless $data_levels->[$i]->[$pindex] == 2;
		if ( $nb_patient < 5 ) {
			$max += $self->control_line2( $patient, $i, $debug );
		}
		else {
			$max += $self->control_line( $self->scores_matrix()->[$i],
				$patient, $debug );
		}
	}

	#warn $max;
	return $max;

	#return $max;
}

sub control_line2 {
	my ( $self, $patient, $i, $debug ) = @_;
	my $primer = $self->transcript->getPrimer( $self->names->[$i] );
	my $ps     = $patient->getFamily->getMembers();
	my @v;
	my %fid;
	my $z = $primer->{cnv};
	foreach my $p (@$ps) {
		delete $z->{ $p->id };
	}
	my @t = values %{$z};
	warn Dumper @t if $debug;
	my $stat = new Statistics::Descriptive::Full;
	$stat->add_data(@t);
	my $sd   = $stat->standard_deviation();
	my $mean = $stat->mean();
	warn $stat->quantile(3) - $stat->quantile(1) if $debug;
	warn $sd . " " . $mean if $debug;

	# warn "---------------";
	return 1 if ( $mean > 0.75 && $sd < 0.3 && $mean < 1.2 );
	return 0;
}

sub control_line {
	my ( $self, $array, $patient, $debug ) = @_;
	confess() unless $patient;

	#my $child_index = $self->hash_index_patients->{$child->id};
	my $hi;

	foreach my $p ( @{ $patient->getFamily->getMembers } ) {
		my $index = $self->hash_index_patients->{ $p->id };
		$hi->{$index}++;
	}
	my @array2;
	my $i = 0;
	foreach my $a (@$array) {
		push( @array2, $a ) unless exists $hi->{$i};
		$i++;
	}

	#my @aarray = @$array;
	#splice @aarray,$index,1;
	#delete()
	my $stat = new Statistics::Descriptive::Full;
	$stat->add_data(@array2);
	my $sd   = $stat->standard_deviation();
	my $mean = $stat->mean();
	warn $sd . " " . $mean if $debug;

	#	 die() if $debug;
	return 1 if ( $mean > 0.75 && $sd < 0.25 );
	return 0;
}

sub control_transmission {
	my ( $self, $child, $mother, $father, $debug ) = @_;
	my $limit       = $self->limit;
	my $patients    = $self->ordered_patients;
	my $data_levels = $self->levels_matrix;
	my $child_index = $self->hash_index_patients->{ $child->id };
	my $father_index;
	if ($father) {
		$father_index = $self->hash_index_patients->{ $father->id };
	}
	my $mother_index;
	if ($mother) {
		$mother_index = $self->hash_index_patients->{ $mother->id };
	}
	warn $child_index . " " . $father_index . " " . $mother_index if $debug;
	my $nb_f      = 0;
	my $nb_m      = 0;
	my $nb_c      = 0;
	my $nb_r      = 0;
	my $recessive = 0;
	my $denovo    = 0;
	my $both      = 0;
	my $other     = 0;
	my $uni       = 0;

	for ( my $i = 0 ; $i < @$data_levels ; $i++ ) {
		my $father_t    = "";
		my $mother_t    = "";
		my $child_t     = "";
		my $level       = $data_levels->[$i]->[$child_index];
		my $score_child = $self->scores_matrix()->[$i]->[$child_index];

		next if ( $score_child > 0.7 and $score_child < 1.2 and $level == 0 );

		warn "\t\t--------------------->"
		  . $score_child . " "
		  . $self->scores_matrix()->[$i]->[$father_index] . " "
		  . $self->scores_matrix()->[$i]->[$mother_index]
		  if $debug;
		$nb_c++;
		$child_t = "he";
		$child_t = "ho" if $score_child < 0.2;
		$child_t = "ho" if $score_child > 1.9;
		if ($father) {
			my $score_father = $self->scores_matrix()->[$i]->[$father_index];
			if (   ( $score_father < 0.7 && $score_child < 0.7 )
				or ( $score_father > 1.3 && $score_child > 1.3 ) )
			{
				$father_t = "he";
				$father_t = "ho" if $score_father < 0.2;
				$father_t = "ho" if $score_father > 1.9;
			}

		}
		if ($mother) {
			my $score_mother = $self->scores_matrix()->[$i]->[$mother_index];
			if ($score_mother) {
				if (   ( $score_mother < 0.7 && $score_child < 0.7 )
					or ( $score_mother > 1.3 && $score_child > 1.3 ) )
				{
					$mother_t = "he";
					$mother_t = "ho" if $score_mother < 0.2;
					$mother_t = "ho" if $score_mother > 1.9;
				}
			}
		}
		warn "\t\t\t\t-+++--------------------> father : "
		  . $father_t
		  . " => mother "
		  . $mother_t
		  . " => child "
		  . $child_t
		  if $debug;
		if ( $child_t eq "ho" ) {
			if ( $father_t eq "he" && $mother_t eq "he" && $child_t eq "ho" ) {
				$recessive++;
			}
			elsif ( !($father_t) && !($mother_t) ) {
				$denovo++;
			}
			else {
				$uni++;
			}

		}
		elsif ( !($father_t) && !($mother_t) ) {
			$denovo++;
		}
		elsif ( $father_t && !($mother_t) ) {
			$nb_f++;
		}
		elsif ( $mother_t && !($father_t) ) {
			$nb_m++;
		}
		elsif ( $mother_t && $father_t ) {
			$both++;
		}
		else {
			$other++;
		}

	}
	warn "---- $nb_f ++ $nb_m ++ $recessive $denovo" if $debug;
	if ( $uni > 0 ) {
		return "unisomy";
	}
	elsif ( $recessive > 0 ) {
		return "recessive";
	}
	elsif ( $denovo > 0 && $nb_m == 0 && $nb_f == 0 ) {
		return "strict-denovo";
	}
	elsif ( $denovo > 0 ) {
		return "denovo";
	}
	elsif ( $nb_f > 0 && $nb_m == 0 ) {
		return "father";
	}
	elsif ( $nb_f == 0 && $nb_m > 0 ) {
		return "mother";
	}
	elsif ( $nb_f > 0 && $nb_m > 0 && $recessive == 0 ) {
		return "both";

	}
	else {
		return "unknown";
	}

}

sub image {
	my ($self)      = @_;
	my $limit       = $self->limit;
	my $patients    = $self->patients;
	my $data_levels = $self->levels_matrix;

	my $size = 3;

	my $x = 0;
	my $y = 0;

	my $xr = 255;
	my $xg = 10;
	my $xb = 0;

	#@color_background  =(0, 165, 145);# (141, 148, 64);
	my $yr = 0;
	my $yg = 255;
	my $yb = 0;

	my $max = 20;

	my $hcolors;
	my $type;

	my $capture = $self->selected_patients->[0]->getCapture;
	my $true_patient;
	my $true_index;
	my $index = 0;
	$max = @{ $self->ordered_patients };

	foreach my $p ( @{ $self->ordered_patients } ) {

		if ( $p->getCapture->id eq $capture->id ) {
			push( @$true_patient, $p ) if $p->getCapture->id eq $capture->id;
			$true_index->{$index}++;

			#push($true_index,$index) ;
		}
		$index++;
	}

	#	die();
	my @primer_ids;
	my $cname = $capture->name;
	for ( my $i = 0 ; $i < @$data_levels ; $i++ ) {

		#for ( my $i = @$data_levels - 1 ; $i > -1 ; $i-- ) {
		my $primer = $self->transcript->getPrimer( $self->names->[$i] );

		#		 warn $primer->name;
		my $vt = join( ":", map { $_->name } @{ $primer->getCaptures } );
		push( @primer_ids, $i ) if $vt =~ /$cname/;

		#push(@primer_ids,$i);
	}
	$max = scalar(@primer_ids) if scalar(@primer_ids) > $max;
	my $nb_col = scalar(@$true_patient);
	my $nb_row = scalar(@primer_ids) + 1;
	my $w      = $nb_col * ( $size + 2 );
	my $h      = ( ( $size + 2 ) * $nb_row );

	my $image = new GD::Image( $w, ( $h + $size ) );
	my $black = $image->colorAllocate( 0, 0, 0 );

	#my $pp    = $image->colorAllocate( 100, 100, 100 );
	my $mimosa   = $image->colorAllocate( 52,  49,  72 );
	my $point    = $image->colorAllocate( 255, 255, 251 );
	my $selected = $image->colorAllocate( 0,   185, 166 );
	$selected = $image->colorAllocate( 244, 235, 251 );
	$selected = $image->colorAllocate( 204, 195, 211 );
	my $selected2 = $image->colorAllocate( 0, 171, 169 );
	for ( my $i = 1 ; $i < 2 ; $i++ ) {
		$x = 2;
		for ( my $j = 0 ; $j < $nb_col ; $j++ ) {

		 #$image->rectangle($x-1,$y-1,$x+$nb_col*($size+2)+1,$y+$size+1,$black);
			my $color = $mimosa;

			$color = $selected
			  if ( exists $self->index_selected_patients->{$j} );
			$image->filledRectangle( $x, $y, $x + $size, $y + $size * 2,
				$color );
			if ( exists $self->index_selected_ill_patients->{$j} ) {
				$image->filledRectangle( $x, $y, $x + $size, $y + $size * 2,
					$selected2 );
				$image->rectangle(
					$x - 1, $y - 1,
					( $x + $size + 1 ),
					$y + $size * 2, $black
				);

			}
			$x += $size + 2;

		}
		$y += $size * 2 + 1;
	}
	my @array_ho;
	for ( my $ii = 0 ; $ii < @primer_ids ; $ii++ ) {

		#for ( my $i = @$data_levels - 1 ; $i > -1 ; $i-- ) {
		# next if  $primer->getCaptures->[0]->name() ne $capture->name;
		#warn  $primer->getCaptures->[0]->name();
		my $i = $primer_ids[$ii];
		$x = 2;
		my $line_levels = $data_levels->[$i];

		#warn  @$line_levels;
		#die();
		my $nbr = grep { $_ == 2 } @{$line_levels};
		my $nbb = grep { $_ == 1 or $_ == 3 or $_ == 4 } @{$line_levels};
		my $nbc   = scalar(@$line_levels);
		my $nbp   = ( $nbr + $nbb ) / $nbc;
		my $nbpr  = ($nbr) / $nbc;
		my $nbpb  = ($nbb) / $nbc;
		my $nbmin = $nbpr;
		if ( $nbpr > 0 ) {
			$nbpr = $nbpr / ( $nbpr + $nbpb );
		}

		#else ($nbpr)

		#$nbmin =   $nbpb if $nbpb < $nbpr;
		#$nbpr =1 unless $nbpr;
		for ( my $j = 0 ; $j < @$line_levels ; $j++ ) {
			next unless exists $true_index->{$j};
			my $blue = $self->blue_rgb;
			my $red  = $self->red_rgb;

			my $green       = $self->green_rgb;
			my $index_color = 0;
			$index_color = 1 if ( exists $self->index_selected_patients->{$j} );

			my @color_background = @{ $green->[$index_color] };

			my @colors = ();
			my $level  = $line_levels->[$j];

			my $score = $self->scores_matrix()->[$i]->[$j];

			@colors = @color_background;

			#			warn $level;
			my $colors;
			if ( $level == -2 ) {
				@colors = ( 90, 90, 90 );

				#$color = $grey1;
			}

			#elsif ( $nbc > 10 && $nbp >= 0.5 && $nbpr>=0.2 && $nbpr<=0.8) {
			#	$type->{blue}++;
			#	@colors = ( 140, 140, 140 );
			#}
			elsif ( $level == -1 ) {
				@colors = ( 140, 140, 140 );
			}
			elsif ( $level == 2 ) {
				@colors = ( 52, 152, 255 );
				@colors = ( 2,  71,  253 );

				#die();
				$type->{blue}++;

				#@colors = @{ $blue->[2] };
				@colors = @{ $blue->[$index_color] }
				  ;    # if ( exists $self->index_selected_patients->{$j} );

			}
			elsif ( $level == 3 ) {
				my $colorr = int( 150 * ( $score - 1 ) ) + 105;

				#@colors = ( 52, 152, 255 );
				#@colors = ( 2,  71,  253 );
				@colors = ( 245, 223, 77 );
				@colors = ( 245, 223, 77 );
				$type->{red}++;
				@colors = @{ $red->[$index_color] };
			}
			elsif ( $level == 1 ) {

				#@colors = (210, 56, 108);
				#@colors = (210, 56, 108);

				#	die();
				#F5DF4D
				$type->{red}++;
				@colors = @{ $red->[$index_color] };
				if ( $score < 0.2 ) {

					#rgb(0, 0, )
					$type->{red_ho}++;

					#@colors = (255, 227, 177);
					@colors = ( 0, 0, 0 );

					#@colors = @colors;
				}

				#@colors = (189, 61, 58);#(233, 75, 60);# @{ $red->[2] };

			}
			elsif ( $level == 4 ) {
				my $colorr = int( 50 * ( 0.6 - $score ) ) + 205;
				@colors = ( 229, 40, 76 );
				@colors = ( 255, 64, 64 );
				@colors = @{ $red->[$index_color] };
			}
			elsif ( $level == 0 ) {

				if ( $score <= 0.55 ) {
					@colors = @{ $red->[ $index_color + 2 ] };

					#@colors = (254, 138, 53)
					#@colors =(222, 165, 164);
				}
				elsif ( $score <= 0.6 ) {
					@colors = @{ $red->[ $index_color + 2 ] };

					#	@colors =(246, 179, 134);
					#						#@colors =(254, 209, 219);
				}
				elsif ( $score <= 0.7 ) {
					@colors = @{ $red->[ $index_color + 2 ] };

					#@colors =(253, 218, 178);
					#rgb(253, 218, 178)
					#						#@colors =(254, 209, 219);
				}

				#					elsif ($score <=0.7){
				#						@colors =(135, 255, 117);
				#						#@colors =(254, 209, 219);
				#					}
				elsif ( $score >= 1.5 ) {
					@colors = @{ $blue->[ $index_color + 2 ] };

					#@colors = (115, 178, 217)
					#@colors =(254, 209, 219);
				}
				elsif ( $score >= 1.4 ) {
					@colors = @{ $blue->[$index_color] };

					#@colors =(205, 223, 245);
				}

			}
			my $id_color = join( ",", @colors );
			unless ( exists $hcolors->{$id_color} ) {
				$hcolors->{$id_color} = $image->colorAllocate(@colors);
			}
			my $x1 = $x;
			my $x2 = $x + ( $size + 1 );
			my $y1 = $y;
			my $y2 = $y + $size;

			if ( $score < 0.2 ) {
				push( @array_ho, [ $x1 - 1, $y1 - 1, $x2, $y2 ] );

			}

			$image->rectangle( $x - 1, $y - 1, $x + 1, $y + $size + 1, $black );

			$image->filledRectangle( $x, $y, $x + $size, $y + $size,
				$hcolors->{$id_color} );

			#$image->filledRectangle( $x, $y, $x + $nb_col * ( $size + 2 ),
			#	$y + $size, $hcolors->{$id_color} );

			$x += $size + 2;
		}
		$y += $size + 2;
	}
	foreach my $ho (@array_ho) {
		$image->rectangle( @$ho, $point );

		#last;
	}

	#next unless $exon->{exon} ==1;
	#$image->rectangle( 0, 0, $w, $h, $black );
	return ( $image, $type );

}

1;

