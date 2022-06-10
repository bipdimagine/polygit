package table_dude;
use strict;
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../";
use lib "$Bin/../../";
use List::Util qw( max min sum);
use Carp qw(confess croak);
use Data::Dumper;

#use Bio::DB::HTS::VCF;
#use Bio::DB::HTS::Tabix;
use Tabix;
use polyweb_dude;
use List::MoreUtils qw{ natatime };
use JSON::XS;

my $VERSION = "23.09.2021";
my $level_value = { "high" => 3, "medium" => 2, "low" => 1 };

my $version_uri= "1";
my $uri_key = "_uri.".$version_uri;

sub getListGenes {
	my ( $patient, $levels, $fork,$print ) = @_;
#	die() unless $print;
	my $no3   = $patient->getGenesDude("w");
	my $lists = [];
	
	my $notok;
	my $hash_level;
	foreach my $level (@$levels) {
		$hash_level->{$level} ++;
		my $test = $no3->get( $level . "updatedate-$VERSION" );
		if ($test){
		my $list2 = $no3->get( $level . "update-$VERSION" );
		
		push( @$lists, @$list2 ) if $list2;
		}
		else {
			$notok ++;
			last;
		}
	}
	warn "first :".scalar(@$lists);
	if ( $notok ) {
		print qq{<div style="display: none">} if $print;
		get_transcripts( $patient, $fork,$print );
		foreach my $level (@$levels) {
		print "!" if $print;
		my $test = $no3->get( $level . "updatedate-$VERSION" );
		die() unless $test;
		my $list2 = $no3->get( $level . "update-$VERSION" );

		#unless($list2)
		push( @$lists, @$list2 );# if exists $hash_level->{$level};
		}
	}
	
	#warn "PRINT ".$print;
	
	print qq{</div>} if $print;
	
	return $lists;

}

sub print_by_position {
	my ($patient,$lists,$cgi,$text) = @_;
	$text .="";
	my $project = $patient->getProject;
	foreach my $chr (@{$project->getChromosomes}) {
		my @slists = grep {$_->{chromosome} eq $chr->name} @$lists;
 		next unless @slists;
 	
	print $cgi->div(
		{ class => "panel-heading" ,
			style =>"background-color:#363945;color:white"
		},
 		$chr->name." ".scalar(@slists)
	);
	
	my $nb =0;
			my @colors = ( "#F9F6FF", "#F7F7F7", "#A9A9D9" );
			my @colors = ( "#FCFCFC", "#F6F6F6", "#A9A9D9" );
	print header_dude_genes2($chr->name,scalar(@slists));	
	my $rids = [];	
	foreach my $hg ( sort { $a->{start} <=> $b->{start} } @slists ) {
		warn "1";
		$nb++;
		my $hide;
		$hide = "";
		if ( $nb > 500 && $hg->{score} < 1 && $hg->{nb_del_ho} == 0 ) {
			$hide = "display:none;";
			push( @$rids, $hg->{rid} );
		}
		my $c = $colors[ $nb % 2 ];
		$hg->{line} =~ s/XDXD/$hide/;
		$hg->{line} =~ s/COLCOL/$c/;
		print $hg->{line};
	}
		print footer_dude_genes($rids);
	}
}
sub print_dude_genes {
	my ( $patient, $lists,$print,$text) = @_;
	unless ($lists){
		print qq{<img src="data:image/png;base64,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">}."NONE ";
		return;
	}
	my $project = $patient->getProject;

	#my $lists = getListGenes( $patient, $levels,$print );

	#warn header_dude_genes();
	my $cgi = new CGI();
	my $cpt = 0;
	my $fork =1;
	my $nb = int( scalar(@$lists) / ($fork) ) + 1;
	my $pm = new Parallel::ForkManager($fork);
	my $iter = natatime $nb, @$lists;
	my $final;
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $h ) = @_;
			unless ( defined($h) or $exit_code > 0 ) {
				print
				  qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
			confess() unless  $h->{array} ;
			push( @$final, @{ $h->{array} } );
		}
	);

	
	while ( my @ltr_tmp = $iter->() ) {
		my $pid       = $pm->start and next;
		my $cpt =0;
		my %lines;
		foreach my $hg (@ltr_tmp) {
		my $gene = $project->newGene( $hg->{id} );
		$hg->{chromosome} = $gene->getChromosome()->name;
		$hg->{start} = $gene->start;
		$hg->{line} = line_dude_gene( $patient, $gene, $hg, $cgi );
	#	my $line = line_dude_gene( $patient, $gene, $hg, $cgi );
	#	$lines->{$hg->{id}} = $line;
		$cpt++;
		print "=" if ( $cpt % 50 == 0 );
		
		}
		$pm->finish( 0, { array => \@ltr_tmp } );
	}
	$pm->wait_all_children();
	warn "END +++++ ";
	
	
	
	
	print qq{<div style="display: none">} if $print;
	
#	foreach my $hg (@$final) {
#		my $gene = $project->newGene( $hg->{id} );
#		$hg->{line} = line_dude_gene( $patient, $gene, $hg, $cgi );
#		$cpt++;
#		print "=" if ( $cpt % 100 == 0 );
#
#		#$hg->{score} = $gene->score;
#	}
	my $rids = [];
	print qq{</div>} if $print ;
	#print_by_position($patient,$final,$cgi);
	#	return;
	print header_dude_genes($text);
	my $nb =0;
	my @colors = ( "#F9F6FF", "#F7F7F7", "#A9A9D9" );
	my @colors = ( "#FCFCFC", "#F6F6F6", "#A9A9D9" );
	my @finale2;
	if ($text){
		 @finale2 = sort { $b->{chr} <=> $a->{chr}  or $a->{start} <=> $b->{start}} @$final ;
	}
	else {
		@finale2 = sort { $b->{score} <=> $a->{score} } @$final;
		
	}
	foreach my $hg (  sort { $b->{score} <=> $a->{score} } @finale2 ) {
#	foreach my $hg ( sort { $b->{score} <=> $a->{score} } @$final ) {
		$nb++;
		my $hide;
		$hide = "";
	#	warn $hg->{score};
		if ($hg->{id} =~/ENSG00000120733/){
	#	warn Dumper $hg;
		#die();
		}
		if ( $nb > 500 && $hg->{score} < 1 && $hg->{nb_del_ho} == 0 ) {
			$hide = "display:none;";
			push( @$rids, $hg->{rid} );
		}
		my $c = $colors[ $nb % 2 ];
		$hg->{line} =~ s/XDXD/$hide/;
		$hg->{line} =~ s/COLCOL/$c/;
		print $hg->{line};
	}

	#warn footer_dude_genes();
	print footer_dude_genes($rids);
}
my $header_cnv_genes = [
	"<center><u>Name</u></center>",
	"<center><u>Type</u></center>",
	"<center><u>Transmission</u></center>",
	"<center><u>Nb Exons</u></center>",
	"<center><u>Locus</u></center>",
	"<center><u>Phenotypes</u></center>",
	"<center><u>Description</u></center>",
	"<center><u>Omim</u></center>",
	"<center><u>pLI</u></center>",
	"<center><u>View CNV</u></center>",
	"<center><u>View VAR</u></center>"
];

sub header_dude_genes {
	my $text = shift;
	my $level = 2;
	my $cgi   = new CGI();
	my $color = "#FF8800";

	my $html;
	$html = $cgi->start_div(
		{
			class => "panel ",
			style =>
"border-color:white;-webkit-border-radius: 3px;-moz-border-radius: 3px;border-radius: 3px;border: 1px solid black;overflow-y:auto;height:100%;"
		}
	);
	$html .= $cgi->div(
		{ class => "panel-heading" ,
			style =>"background-color:#363945;color:white;font-size:12px"
		},
qq{Dup Del - CNV Event : $text }
	);

	$html .=
	  qq{<div style="overflow-y:auto;border:solid 1px grey;">};
	$html .= $cgi->start_table(
		{
			class =>
"table table-sm table-striped table-condensed table-bordered table-primary ",
			style =>
"box-shadow: 1px 1px 6px $color;font-size: 9px;font-family:  Verdana;margin-bottom:0px;"
		}
	);

	#$header_cnv_genes ="XXXXXXXXX";
	my $nb = 0;

	$html .= $cgi->start_Tr(
		{
			style =>
"position:sticky;z-index:9;top:0;background-color:grey;color:white;"
		}
	  )
	  . $cgi->th( { style => "text-align: center;" }, $header_cnv_genes ) . ""
	  . $cgi->end_Tr();
	return $html;
}


sub footer_dude_genes {
	my ($rids) = @_;
	my $cgi    = new CGI();
	my $html   = "";
	if ( scalar(@$rids) >= 0 ) {

		my $js      = encode_json $rids;
		my $za      = "hide_tr_" . time . "_" . int( rand(50000) );
		my $nb_skip = scalar(@$rids);
		$html .= $cgi->start_Tr( { id => $za } );
		$html .= $cgi->td(
			{
				style =>
				  "box-shadow: 1px 1px 2px #555;background-color:#CECFCE;",
				colspan => scalar(@$header_cnv_genes),
				onClick => qq{showTranscripts($js,"$za");}
			},
			qq{<span class="glyphicon glyphicon-plus"></span> }
			  . "view $nb_skip Genes"
		);
		$html .= $cgi->end_Tr();

	}
	$html .= $cgi->end_table();
	$html .= qq{</div>};
	$html .= qq{</div>};
	$html .= qq{</div>};
	return $html;
}

my $nb_row = 0;

sub line_dude_gene {
	my ( $patient, $gene, $hgene_dude, $cgi ) = @_;
	my $html;
	my $h_chr_cnv;
	my $project = $patient->getProject();


	my $rid = "row_gene_" . $gene->id . "_" . time . "_" . $nb_row;
	$gene->score;
	my $nb_dup = $hgene_dude->{nb_dup};
	my $nb_del = $hgene_dude->{nb_del};
	my $nb_ho  = $hgene_dude->{nb_del_ho};
	my $nb_all = $hgene_dude->{nb_all};
	my $run_id = $patient->getCapture()->id;
	 my $nb_patients = scalar(grep {$_->getCapture->id eq $run_id} @{$patient->getProject->getPatients()});
	# $nb_patients =  scalar(@{$project->getPatients()});
	#my $nb_patients = scalar(@{$project->getPatients()});
	my $lmax = $nb_all;
	$lmax = $nb_patients if $lmax < $nb_patients;
	my $type_scale = "height";
	 $type_scale = "width" if $lmax < $nb_patients;
	$hgene_dude->{score} = $gene->score;
	

	my $width = ($nb_patients+2) *3;
	$width = ($nb_patients+2) *4 if ($width) < 50 ; 
	my $height = ($nb_all+1) *3 + 8;
	$height  = int( ($nb_all+1) *1.5 + 8)  if $height > 400;
	if ($width < 50 or $height > 400 ) {
		$width = undef;
		$height = undef;
	}
	if ( $nb_dup == 0 && $nb_del > 2 && $hgene_dude->{level} eq "high" ) {
		$hgene_dude->{score} += $nb_del * 0.5;
		$hgene_dude->{score} += $nb_ho ;
	}



	$html = $cgi->start_Tr( { id => $rid, style => "border: 1px solid;background-color:COLCOL;XDXD" } );
	$html = $cgi->start_Tr( { id => $rid, style => "border: 1px solid;background-color:;#B565A7;XDXD" } ) if $nb_ho > 0;

	if ( $nb_ho > 0 ) {

		#$hgene_dude->{score} += 0.5 * $nb_ho;
		$hgene_dude->{score} += 1 * $nb_ho if $hgene_dude->{level} eq "high";
		#$c    = $colors[2];
		
		$html = $cgi->start_Tr(
			{
				id => $rid,
				style =>
"border: 1px black solid ;background-color:#FDAC53;XDXD;box-shadow: 1px 1px 2px #555;"
			}
		) if $nb_ho > 0;

		#$gene->{score} += $nb_del * 0.3;
	}

	#$gene->score

	#$c = $colors[2] if $nb_ho > 0;

	$hgene_dude->{rid} = $rid;

	#GENE NAME
	my $in = $gene->omim_inheritance;
	$in = ""          if $in eq "-";
	$in = "X-linked " if $in =~ /X-linked/;
	$in ="<sup>".$in."</sup>"; 
	my $uc = qq{https://gnomad.broadinstitute.org/gene/} . $gene->name;
	my $oc = qq{onClick='window.open("$uc")'};

	#TODO: here
	my $gene_id = $gene->id;

	#PHENOTYPE
	my ($pheno,$nb_other_terms) =$gene->polyviewer_phentotypes();
		
	   	my $b_span_pheno = '';
   			if ($pheno) {
		   		if (length($pheno) > 70) { $pheno = substr($pheno, 0, 70).'...'; }
   				if ($nb_other_terms > 0) { $pheno .= " <span style='color:blue'>+ $nb_other_terms terms</span>"; }
   				$b_span_pheno .= qq{<a class="btn btn-xs" role="button" style="font-size:8px;background-color:#EEE;color:black;border:solid 1px black;" onclick="update_grid_gene_phenotypes(\'$gene_id\')";"><span>$pheno</span></a>};
   			}
#	my $to;
#	( $pheno, $to ) = split( /\[/, $gene->description ) unless $pheno;
	
		#$b_span_pheno .=
#qq{<a class="btn btn-xs" role="button" style="font-size:9px;background-color:#EEE;color:black;border:solid 1px black;" onclick="update_grid_gene_phenotypes(\'$gene_id\')";"><span>$pheno</span></a>} ;
	

	#OMIM
	my $b_omim;
	my $omim = $gene->omim_id();
	if ( $omim ne "" ) {
	$b_omim .=
		qq{<a class="btn btn-xs" style="font-size:9px;background-color:#EEE;color:black;border:solid 1px black;" href="http://www.omim.org/entry/$omim" role="button" target="_blank">Omim</a>};
	}
	else {
		$b_omim .=
qq{<a class="btn btn-xs" style="font-size:9px;background-color:#EEE;color:black;border:solid 1px black;" role="button">-</a>}
		  if $omim eq "";
	}

	#GTEX
	my ( $gid, $t ) = split( "_", $gene_id );

#my $b_gtex = qq{<a class="btn btn-xs" href="https://gtexportal.org/home/gene/$gid" role="button" target="_blank" style="font-size:9px;background-color:#EEE;color:black;border:solid 1px black;">Gtex</a>};

	#PLI
	my $pli  = $gene->pLI * 1.0;
	my $type = "green";
	$type = "orange" if $pli >= 0.75;
	$type = "red"    if $pli >= 0.9;
	my $b_pli =
qq{<a class="btn btn-xs" role="button" href='https://gnomad.broadinstitute.org/gene/$gid' target='_blank' style="background-color:#EEE;color:black;border:solid 1px black;"><span style="font-size:9px;background-color:#EEE;">$pli</span></a>};

	#VARIANTS
	my ( $nb_dude, $nb_var );
	$nb_var = 'View';

	my $patient_name = $patient->name();
	my $project_name = $project->name();
	my $cmd_var =
qq{dijit.byId('dialog_hgmd').show();view_var_from_proj_gene_pat('$project_name', '$gene_id', '$patient_name', '', 'all', 'nocnv');};
	my $b_var =
qq{<a class="btn btn-xs" role="button" onclick="$cmd_var" style="background-color:#3AB795;color:white;border:solid 1px black;" disabled><span style="font-size:9px;" class="glyphicon glyphicon-ban-circle"></span></a>};

	my $level_dude_text;
	my $level_dude = 'Not';
	if ( $hgene_dude->{level} eq "high" ) {
		$level_dude      = 'High';
		$level_dude_text = "high";
	}
	elsif ( $hgene_dude->{level} eq "medium" ) {
		$level_dude      = 'Med';
		$level_dude_text = "medium";
	}
	else {
		$level_dude      = 'Low';
		$level_dude_text = "low";
	}
	
	my $bname = printButton(
		$gene->score,
		[ 1, 4 ],
		$gene->external_name.' <b><u>' . $in . '</b></u>', $oc
	);

	#$gene->name = $bname;
	$html .= $cgi->td( "<center>" . $bname . "</center>" );

	my $noise  = $hgene_dude->{noise};
	my $del_ho = $hgene_dude->{nb_del_ho};

	my $level_dude_noise = '';
	if ( $noise == 0 ) {
		$level_dude_noise = qq{<button type="button" class="btn btn-xs" style="background-color:#BC243C;color:white;border:solid 1px black;font-size:10px;">0%</button>};
	}
	else {
		$level_dude_noise = qq{<button type="button" class="btn btn-xs" style="background-color:white;color:black;border:solid 1px black;font-size:9px;">$noise%</button>};
	}

	my $color_dude = 'grey';
	$color_dude = '#BC243C' if ( lc($level_dude) eq 'high' );
	$color_dude = '#EFC050' if ( lc($level_dude) eq 'med' );
	$color_dude = '#98B4D4' if ( lc($level_dude) eq 'low' );

	my $b_cnv =
qq{<button type="button" class="btn btn-xs" onclick="zoomDude(\'$project_name\','$gene_id','$patient_name', 'force_visualisation')" style="background-color:$color_dude;color:white;border:solid 1px black;font-size:9px;">$level_dude</button>};
	if ( $nb_var eq 'View' or $nb_var > 0 ) {
		$b_var =
qq{<a class="btn btn-xs" role="button" onclick="$cmd_var" style="background-color:#3AB795;color:white;border:solid 1px black;"><span style="font-size:9px;">$nb_var</span></a>};
	}

	my $no3 = $patient->getGenesDude("r");
		warn $hgene_dude->{best_transcript};
		warn $uri_key;
	my $uri_text = $no3->get( $hgene_dude->{best_transcript} . $uri_key."-$VERSION" );
	if ($uri_text) {
		my $url2 = $uri_text;
		my $size = "";
		if ($width && $height){
				#my $hg = $heigth."px";
		#my $wg = $width."px";
			$size = "width:$width"."px height:$height"."px";
		}
	
		$b_cnv =
	#qq{<div style="padding:0px;min-width:50px;min-height:50px" > $level_dude <br><img class="$class" src="$url2"  align="top" style="box-shadow: 2px 2px 9px $color_dude;vertical-align: top; "  onclick="zoomDude(\'$project_name\','$gene_id','$patient_name', 'force_visualisation')"></img></div>};

	qq{<img class="zoom1" src="$url2"  align="top" style="box-shadow: 2px 2px 9px $color_dude;vertical-align: top; "  onclick="zoomDude(\'$project_name\','$gene_id','$patient_name', 'force_visualisation')"></img>};

	}

	# LOCUS (EXOMES)
	my $chr_id     = $gene->getChromosome->id();
	my $gene_start = $gene->start();
	my $gene_end   = $gene->end();
	my $locus =
		$chr_id . ":"
	  . $gene_start . "-"
	  . $gene_end 
	  ;
	  
	  my $color = $patient->buffer->color_model($hgene_dude->{transmission});
	  my $b_transmission = qq{<button class= "btn btn-xs btn-primary " style="background-color: $color;font-size: 1em;font-family:  Verdana;color:black">}.$hgene_dude->{transmission}.qq{<button};
	  
	  
	my $ucsc =
	  qq{https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=$locus};
	my $b_locus =
qq{<button class= "btn btn-xs btn-primary " style="background-color: #EEE;font-size: 1em;font-family:  Verdana;color:black"><a href=$ucsc target='_blank' style="color:black;">}
	  . $locus
	  . qq{</a> ++ $del_ho</button>};

	#TYPE
	my ( $b_type, $b_exons );
	if ( $del_ho > 0 ) {
		$b_type =
qq{<button type="button" class="btn btn-xs" style="background-color:purple;color:white;border:solid 1px black;font-size:9px;">del HO</button>};
		$b_exons =
qq{<button type="button" class= "btn btn-xs btn-primary " style="background-color:purple;font-size: 1em;font-family:Verdana;color:white">$nb_del/$nb_all</button>};

	}
	elsif ( $nb_dup > $nb_del ) {
		$b_type =
qq{<button type="button" class= "btn btn-xs btn-primary " style="background-color: #217DBB;font-size: 1em;font-family:Verdana;color:white">gain</button>};
		$b_exons =
qq{<button type="button" class= "btn btn-xs btn-primary " style="background-color: #217DBB;font-size: 1em;font-family:Verdana;color:white">$nb_dup/$nb_all</button>};
	}
	elsif ( $nb_dup < $nb_del ) {
		$b_type =
qq{<button type="button" class= "btn btn-xs btn-primary " style="background-color: #E74C3C;font-size: 1em;font-family:Verdana;color:white">del</button>};

		$b_exons =
qq{<button type="button" class= "btn btn-xs btn-primary " style="background-color:#E74C3C;font-size: 1em;font-family:Verdana;color:white">$nb_del/$nb_all</button>};
	}
	else {
		my $nb_dup_del = $nb_dup + $nb_del;
		$b_type =
qq{<button type="button" class= "btn btn-xs btn-primary " style="background-color: grey;font-size: 1em;font-family:Verdana;color:white">gain/del ?</button>};
		$b_exons =
qq{<button type="button" class= "btn btn-xs btn-primary " style="background-color: grey;font-size: 1em;font-family:Verdana;color:white">$nb_dup_del/$nb_all</button>};
	}

	my @lbam_alamut;
	foreach my $p ( @{ $patient->getFamily->getPatients() } ) {
		push( @lbam_alamut, 'https://www.polyweb.fr/' . $p->bamUrl() );
	}
	my $string_url_bam = join( ',', @lbam_alamut );
	my $b_igv =
qq{<button type="button" class="btn btn-info btn-xs" style="background-color:white;font-size:9px;" onClick="displayInIGV('$chr_id', '$gene_start', '$gene_end');"><img style="width:16px;height:16px;" src="images/polyicons/igv_logo_32.png"></img></button>};
	my $b_alamut =
qq{<button type="button" class="btn btn-danger btn-xs" style="background-color:white;font-size:9px;" onClick="displayLocusInAlamut('$chr_id', '$gene_start', '$gene_end');"><img style="width:16px;height:16px;" src="images/polyicons/alamut_visual.png"></img></button>};

	$html .= $cgi->td( "<center>" . $b_type . "</center>" );
	$html .= $cgi->td( "<center>" . $b_transmission . "</center>" );
	$html .= $cgi->td( "<center>" . $b_exons . "</center>" );

	$html .= $cgi->td( "<center>" . $b_locus . "</center>" );

	#			$html.= $cgi->td("<center>".$tab_graph."</center>");
	$html .= $cgi->td( "<center>" . $b_span_pheno . "</center>" );
	$html .= $cgi->td( "<center>" . $gene->description() . "</center>" );

	#$html.= $cgi->td("<center>".$b_panels."</center>");
	$html .= $cgi->td( "<center>" . $b_omim . "</center>" );

	#$html.= $cgi->td("<center>".$b_gtex."</center>");
	$html .= $cgi->td( "<center>" . $b_pli . "</center>" );
	$html .= $cgi->td( "<center>" . $b_cnv . "</center>" );
	$html .= $cgi->td( "<center>" . $b_var . "</center>" );

	$html .= $cgi->end_Tr();

	#$atr->{html} = $html;
	return $html;

}

sub printButton {
	my ( $value, $types, $text, $othercss, $text_alert ) = @_;

	my $btn_class =
qq{class= "btn btn-xs btn-primary " style="background-color: #D0D0D0;font-size: 9px;font-family:  Verdana;color:black"};
	$btn_class =
qq{class= "btn btn-xs btn-primary " style="background-color: #FF8800;font-size: 9px;font-family:  Verdana;;color:white"}
	  if $value >= $types->[0];
	$btn_class =
qq{class= "btn btn-xs  btn-primary" style="background-color: #e74c3c;font-size: 9px;font-family:  Verdana;color:white"}
	  if $value >= $types->[1];
	$value = "-" unless defined $value;

	$btn_class =
qq{class= "btn btn-xs  btn-primary" style="background-color: red;font-size: 9px;font-family:  Verdana;"}
	  if $value eq "-";
	$text = $value unless $text;
	return
qq{<button onClick="alert('$text_alert')" "type="button" $btn_class $othercss>$text</button>}
	  if ($text_alert);
	return qq{<button type="button" $btn_class $othercss>$text</button>};
}

sub get_transcripts {
	my ( $patient, $fork, $print ) = @_;
	my $no2 = $patient->getTranscriptsDude("r");
	system( "vmtouch -q -t " . $no2->filename );

	#$no2->close;
	$patient->getProject->getPatients();
	my @levels = ( "high", "medium", "low" );
	my $aa = [];
	my %dj;
	foreach my $l (@levels) {
		my $array = $no2->get("$l");
#		warn $l." ".$patient->name;
		foreach my $t (@$array) {
			my $z = {};
			next if exists $dj{$t};
			$dj{$t}++;
			$z->{transcript} = $t;
			$z->{old_level}  = $l;
			push( @$aa, $z );
		}

		#my $ts = $project->newTranscripts($array);
	}
	$no2->close();
	#$print = 1;
	my $nb = int( scalar(@$aa) / ($fork) ) + 1;
	my $pm = new Parallel::ForkManager($fork);
	my $final;
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $h ) = @_;
			unless ( defined($h) or $exit_code > 0 ) {
				print
				  qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
			confess() unless  $h->{array} ;
			push( @$final, @{ $h->{array} } );
		}
	);

	my $iter = natatime $nb, @$aa;
	while ( my @ltr_tmp = $iter->() ) {
		my $pid       = $pm->start and next;
		my $new_lists = update_level( $patient, \@ltr_tmp, $print );
		$pm->finish( 0, { array => $new_lists } );
	}
	$pm->wait_all_children();
	warn "END +++++ ";

	#die();
	my $zz;

	my $id_level;
	foreach my $z (@$final) {
		my $gene_id = $z->{gene_id};

		unless ( exists $zz->{$gene_id} ) {
			$zz->{$gene_id}->{id}          = $gene_id;
			$zz->{$gene_id}->{level}       = $z->{new_level};
			$zz->{$gene_id}->{level_value} = $level_value->{ $z->{new_level} };
			$zz->{$gene_id}->{nb_dup}      = $z->{nb_dup};
			$zz->{$gene_id}->{nb_del}      = $z->{nb_del};
			$zz->{$gene_id}->{nb_del_ho}   = $z->{nb_del_ho};
			$zz->{$gene_id}->{nb_all}      = $z->{nb_all};
			$zz->{$gene_id}->{best_transcript} = $z->{transcript};
			$zz->{$gene_id}->{transcripts}->{ $z->{transcript} } =
			  $z->{new_level};
			$zz->{$gene_id}->{transmission} = $z->{transmission};
			$zz->{$gene_id}->{description}  = $z->{description};

   #push(@{$zz->{$gene_id}->{transcripts}},{$z->{transcript}=>$z->{new_level}});
			$id_level->{ $z->{new_level} }->{$gene_id}++;
			foreach my $t ( keys %{ $z->{uri} } ) {
				$zz->{$gene_id}->{uri}->{$t} = $z->{uri};
			}
		}
		else {
			$zz->{$gene_id}->{transcripts}->{ $z->{transcript} } =
			  $z->{new_level};
			next
			  if $zz->{$gene_id}->{level_value} >=
			  $level_value->{ $z->{new_level} }
			  and $zz->{$gene_id}->{nb_all} > $z->{nb_all};
			delete $id_level->{ $zz->{$gene_id}->{level} }->{$gene_id};
			$zz->{$gene_id}->{id}          = $gene_id;
			$zz->{$gene_id}->{level}       = $z->{new_level};
			$zz->{$gene_id}->{level_value} = $level_value->{ $z->{new_level} };
			$zz->{$gene_id}->{nb_dup}      = $z->{nb_dup};
			$zz->{$gene_id}->{nb_del}      = $z->{nb_del};
			$zz->{$gene_id}->{nb_del_ho}   = $z->{nb_del_ho};
			$zz->{$gene_id}->{nb_all}      = $z->{nb_all};
			$id_level->{ $z->{new_level} }->{$gene_id}++;
			$zz->{$gene_id}->{transmission}    = $z->{transmission};
			$zz->{$gene_id}->{description}     = $z->{description};
			$zz->{$gene_id}->{best_transcript} = $z->{transcript};
		}

		#$zz->{$z->{gene}}++;
	}
	warn "END COMPUTE";
	eval {
		my $no3 = $patient->getGenesDude("w");
		foreach my $l (@levels) {
			my $list1 =[];
			my $t = time;
			my $test = $no3->put( $l."updatedate-$VERSION",{date=>$t} );
			foreach my $gene_id ( keys %{ $id_level->{$l} } ) {
				$no3->put( $gene_id . $uri_key."$VERSION", $zz->{$gene_id}->{uri} );
				delete $zz->{$gene_id}->{uri};
				foreach my $t ( keys %{ $zz->{$gene_id}->{uri} } ) {

				}
				push( @$list1, $zz->{$gene_id} );
			}
			$no3->put( $l . "update-$VERSION", $list1 );
		}
		$no3->close();
	}

}

sub update_level {
	my ( $patient, $lists, $print ) = @_;
	$patient->getProject->buffer->dbh_deconnect;
	my @lPatients = @{ $patient->getProject->getPatients() };
	my @selected_patients;
	push( @selected_patients, $patient );
	my $level_high;
	my $project = $patient->getProject;

	my $hGenes_dude_tmp;
	my $xx = 0;
	
	foreach my $h1 (@$lists) {
		my $current_level = $h1->{old_level};
		$xx++;
		print "&" if $xx % 200 == 0 && $print;
		my $t = $project->newTranscript( $h1->{transcript} );
		$h1->{gene}                         = $t->getGene->external_name;
		$h1->{gene_id}                      = $t->getGene->id;
		$h1->{description}->{external_name} = $t->getGene->external_name;
		my $gene       = $t->getGene;
		my $chr_id     = $gene->getChromosome->id();
		my $gene_start = $gene->start();
		my $gene_end   = $gene->end();
		my $locus      = $chr_id . ":" . $gene_start . "-" . $gene_end;

		#$h1->{description}->{locus} = $chr_id.":".$gene_start."-".$gene_end;
		$h1->{description}->{chromosome} = $gene->getChromosome->name;
		$h1->{description}->{start}      = $gene->start;
		$h1->{description}->{end}        = $gene->end;
		$h1->{description}->{enst}       = $t->name;
		$h1->{description}->{nm}         = $t->external_name;
		my $pheno = $gene->phenotypes;
		if ( $gene->omim->{phenotype}->{omim} ) {
			my (@phenos) =
			  $gene->omim->{phenotype}->{omim} =~ /{([^}]+)}/g; # /{((.+?))}/g ;
			if (@phenos) {
				$h1->{description}->{text} = join( ";", @phenos );
			}
			else {
				my @t = split( ";", $gene->omim->{phenotype}->{omim} );
				foreach my $st (@t) {
					my @stt = split( ",", $st );
					$h1->{description}->{text} .= $stt[0] . ";";
				}
			}

		}
		else {
			$h1->{description}->{text} = $gene->description;
		}

		# FILTRE TRANSCRIPTS HIGH DUDE
		my $coverage = polyweb_dude->new(
			patients          => \@lPatients,
			transcript        => $t,
			limit             => undef,
			selected_patients => [$patient]
		);
		$coverage->init_matrices();

		my $hcov   = $coverage->quality;
		my $nb_dup = $hcov->{ $patient->name() }->{dup};
		my $nb_del = $hcov->{ $patient->name() }->{del};
		my $nb_del_ho += $hcov->{ $patient->name() }->{del_ho};
		my $debug;
		#$debug = 1 if $gene->external_name() eq "KITLG";
		warn $current_level if $debug;
		if ( $nb_del_ho > 0 ) {
			warn "HO =====" if $debug;
			my $znb = $coverage->control_ho($patient,$debug);
			$h1->{raw_nb_del_ho} = $nb_del_ho;
			$nb_del_ho = $znb;
		}

		warn "HO:".$nb_del_ho if $debug;
#		die() if $debug;
	
		$nb_del += $nb_del_ho;
		$h1->{nb_dup}    = $nb_dup + 0;
		$h1->{nb_del}    = $nb_del + 0;
		$h1->{nb_del_ho} = $nb_del_ho + 0;
		my $nb_all = $hcov->{ $patient->name() }->{all};
		$h1->{nb_all} = $nb_all + 0;
		my $nb_grey_others;
		my $nb_noise;
		my $nb_all_others;
		my $nb_ho_others;

		foreach my $other_patient ( @{ $project->getPatients() } ) {

			next if ( $other_patient->getFamily->name() eq $patient->getFamily->name() );

			#					warn $other_patient->name;
			$nb_grey_others += $hcov->{ $other_patient->name() }->{grey};
			$nb_noise +=
			  $hcov->{ $other_patient->name() }->{del} +
			  $hcov->{ $other_patient->name() }->{dup} +
			  $hcov->{ $other_patient->name() }->{del_ho};
			$nb_all_others += $hcov->{ $other_patient->name() }->{all};
			$nb_ho_others  += $hcov->{ $other_patient->name() }->{del_ho};

		}
			
		my $perc          = 0;
		my $perc_grey     = 0;
		my $perc_noise    = 0;
		my $perc_gene_dup = 0;
		my $perc_gene_del = 0;

		if ($nb_all_others) {
			$perc_grey  = ( $nb_grey_others / $nb_all_others ) * 100;
			$perc       = ( $nb_noise / $nb_all_others ) * 100;
			$perc_grey  = ( $nb_grey_others / $nb_all_others ) * 100;
			$perc_noise = ( $nb_noise / $nb_all_others ) * 100;

			#on est high et on baisse :
			$perc_gene_dup = ( $nb_dup / $nb_all ) * 100;
			$perc_gene_del = ( $nb_del / $nb_all ) * 100;
		}
		if ( $nb_dup == 0 && $nb_del == 0 ) {
			$h1->{new_level} = "low";
		}

		if ( $patient->isChild ) {
			warn "TRANS -------------------- " if $debug;
			$h1->{transmission} = $coverage->control_transmission($patient,$patient->getFamily->getMother,$patient->getFamily->getFather, $debug);
			warn "trans ".$h1->{transmission} if $debug;
			warn "---------------- -------------------- " if $debug;
		}
		warn $h1->{new_level}."-----------------" if $debug;
		if ( $nb_del_ho > 0 ) {
			if ( $perc_noise > 40 ) {
				$h1->{new_level} = "low";
			}
			elsif ( $perc_noise >= 10 ) {
				$h1->{new_level} = "medium";
			}
			elsif ( $nb_ho_others > 5 ) {
				$h1->{new_level} = "medium";
			}
			else { $h1->{new_level} = "high"; }
			next;
		}
		my $level = $current_level;
		warn " $level-----------------" if $debug;
		$h1->{new_level} = adjust_high_level( $nb_dup, $nb_del, $perc, $perc_grey, $perc_noise, $perc_gene_dup, $perc_gene_del ) if $current_level eq "high";
			warn $h1->{new_level}." $level-----------------" if $debug;
		$h1->{new_level} =   adjust_medium_level( $nb_dup, $nb_del, $perc, $perc_grey, $perc_noise, $perc_gene_dup, $perc_gene_del ) if $current_level eq "medium";
		$h1->{new_level} =   adjust_low_level( $nb_dup, $nb_del, $perc, $perc_grey, $perc_noise,$perc_gene_dup, $perc_gene_del ) if $current_level eq "low";
		warn $h1->{new_level}."-----------------" if $debug;
		#$h1->{new_level} = "low";
		
		if ( $nb_del > 0 && $nb_del < 50 && $nb_dup == 0 ) {
			my $znb = $coverage->control_del($patient,$debug);
			warn $znb if $debug;
			die() if $debug;
	#		die() if $debug;
			if ( $znb < 1 ) {
				$h1->{new_level} = "low";
			}
		}
	warn $h1->{new_level}."-----------------" if $debug;
		#	warn $nb_dup." ".$nb_del if $debug;
		if ( $nb_dup > 0 && $nb_dup < 3 && $nb_del == 0 ) {
			my $znb = $coverage->control_dup($patient,$debug);
			if ( $znb < 1 ) {
				$h1->{new_level} = "low";
			}

			#warn $nb_dup." ==> ".$znb if $debug;
		}
		warn $h1->{new_level}."-----------------" if $debug;
		die() if $debug;
			warn "OCUCU" if $debug;
				

		#					my ( $image, $type ) = $coverage->image();
		#					my $uri = "";
		#					my $uri = URI->new("data:");
		#					$uri->media_type("image/png");
		#					$uri->data( $image->png );
		#					$h1->{uri}->{$t->id}  = $uri;
		$h1->{ $t->id }->{external_name} = $t->external_name;

	}
	return $lists;
}

sub adjust_low_level {
	my ( $nb_dup, $nb_del, $perc, $perc_grey, $perc_noise, $perc_gene_dup,
		$perc_gene_del )
	  = @_;
	my $level = "low";

	#if ($perc_noise < 1) {
	#		return "medium";
	#}
	if ( $perc_noise < 2 ) {
		return "medium";
	}
	return $level;
}

sub adjust_medium_level {
	my ( $nb_dup, $nb_del, $perc, $perc_grey, $perc_noise, $perc_gene_dup,
		$perc_gene_del )
	  = @_;
	my $level = "medium";
	if ( $perc_noise < 3 ) {
		return "high";
	}
	return $level;
}

sub adjust_high_level {
	my ( $nb_dup, $nb_del, $perc, $perc_grey, $perc_noise, $perc_gene_dup,
		$perc_gene_del, $debug )
	  = @_;
	my $level = "high";

	if ( $nb_dup == $nb_del ) {
		warn "coucou " . $debug if $debug;
		if ( $perc_noise >= 20 ) {
			return "low";
		}
		else {
			return "medium";
		}
	}

	if ( $perc_gene_dup >= 60 or $perc_gene_del >= 60 ) {
		if ( $perc_noise >= 20 ) {
			return "medium";
		}
		return "high";

	}

	#CAS baisse si BCP TROP bruit
	if ( $perc_noise >= 40 ) {

		#warn "coucou ".$debug if $debug;
		return "low";
	}
	if ( $perc_noise >= 5 ) {
		warn "coucou " . $debug if $debug;
		return "medium";
	}

	return $level;

}

sub get_images {
	my ( $patient, $hg, $fork, $print ) = @_;
	print qq{<div style="display: none">} if $print;
	#$fork = 1;
	
	my $nb = int( scalar(@$hg) / ($fork) ) + 1;
	my $pm = new Parallel::ForkManager($fork);
	my $final;

	warn "start get images";
	my $iter = natatime $nb, @$hg;
	my $save = {};
	
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $h ) = @_;
			unless ( defined($h) or $exit_code > 0 ) {
				print
				  qq|No message received from child process $exit_code $pid!\n|;
				die();
				return;
			}
			return unless $h->{array};
			push( @$final, @{ $h->{array} } );
			foreach my $h ( @{ $h->{array} } ) {
				foreach my $t ( keys %{ $h->{save} } ) {
					$save->{$t} = $h->{uri}->{$t};
				}

			}
		}
	);
	$patient->project->preload_patients();
	$patient->project->buffer->dbh_deconnect();
	while ( my @ltr_tmp = $iter->() ) {
		my $pid       = $pm->start and next;
		my $new_lists = compute_uri( $patient, \@ltr_tmp, $print );
		$pm->finish( 0, { array => $new_lists } );
	}
	$pm->wait_all_children();
	warn "end save iamge";
	save_images( $patient, $save );
	print "</div>";
	warn "end";
	return $final;

}

sub save_images {
	my ( $patient, $uri ) = @_;
	eval {
		my $no3 = $patient->getGenesDude("w");
		warn "SAVE :" . scalar( keys %$uri );
		my $nb =0;
		foreach my $l ( keys %$uri ) {
			warn "s:$nb\n" if $nb %100 ==0;
			$nb ++;
			$no3->put( $l . $uri_key."-$VERSION", $uri->{$l} );
		}
		warn "close" ;
		$no3->close();
	};
}

sub compute_uri {
	my ( $patient, $hg, $print ) = @_;

	my $project = $patient->getProject();
	$project->getPatients();
	my $no3 = $patient->getGenesDude("r");
	my $cpt = 0;
	my $cpx = 0;
	foreach my $h (@$hg) {
		$cpt++;
		print "!" if $cpt % 100 && $print;

		#foreach my $ht  (keys %{$h->{transcripts}}) {

		my $uri_text;    #
		if ($print) {
			$uri_text = $no3->get( $h->{best_transcript} . $uri_key."-$VERSION" );
		
		}
		unless ($uri_text) {
			my $t        = $project->newTranscript( $h->{best_transcript} );
			my $coverage = polyweb_dude->new(
				patients          => $project->getPatients,
				transcript        => $t,
				limit             => undef,
				selected_patients => $patient->getFamily->getMembers
			);
			my ( $image, $type ) = $coverage->image();
			my $uri = URI->new("data:");
			$uri->media_type("image/png");
			$uri->data( $image->png );
			$uri_text = $uri->as_string;
			$h->{save}->{ $h->{best_transcript} }++;

			
		}
		$h->{uri}->{ $h->{best_transcript} } = $uri_text;

	
	}
	return $hg;
}

1;
