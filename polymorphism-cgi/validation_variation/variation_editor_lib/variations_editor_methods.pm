package variations_editor_methods;
use strict;
use Data::Dumper;
use Time::HiRes qw( time);

sub refine_heterozygote_composite {
	my ($project,$print_html,$list, $id, ) = @_;
		if ($project->isRocks){
		return new_refine_heterozygote_composite_score_rocks($project, $print_html,$list, $id);
		}
		else {
			return new_refine_heterozygote_composite_score_lmdb($project, $print_html,$list, $id);
		}

	
}

sub return_vp {
	 my ($project,$chr,$patient,$vid,$g) =@_;
	 my $debug;
	  if ($project->isRocks){
	  		my $no  =  $chr->rocks_polyviewer_variants("r");
			my $vp = $no->test($vid,$patient->id,$g->{id});
			die() if $debug;
			unless ($vp){
				confess();
				die($vid) unless keys %{$vp};
			}
			$vp->gene($g);
			return $vp;
	  }
	   my $no  = $chr->lmdb_polyviewer_variants( $patient, "r" );
	   my $vh = $no->get($vid);
			
		my $vp =  PolyviewerVariant->new();
			
		$vp->setOldVariant($vh,$project,$patient,$g);
		return $vp;
	  }
	my @headers = (
				"varsome", "igv",    "alamut", "var_name",
				"trio",    "gnomad", "deja_vu"
			);

sub new_refine_heterozygote_composite_score_rocks {
	my ( $project,$print_html,$list, $id ) = @_;
	
	my $cgi = $print_html->cgi;
	my $patient= $print_html->patient;
	my $final_polyviewer_all = GenBoNoSqlRocks->new(dir=>$project->rocks_directory("polyviewer_global"),mode=>"r",name=>$project->name);
	my $noP =  GenBoNoSqlRocksPolyviewerVariant->new(dir=>"/data-beegfs/test-cache/polyviewer/",mode=>"r",name=>$project->name.".polyviewer_variant");
	my $out_header;
	$out_header= $print_html->print_header("background-color:#E9DEFF");
	my $mother = $patient->getFamily->getMother();
	my $father = $patient->getFamily->getFather();
	my $hno;
	my $tsum = 0;
	my $t    = time;
	my $xp   = 0;
	my $time_decode =0;
	my $time_gene =0;
	#$t = time;
	my $current;
	my $rtime = time;
	my $total_time = time;
	my $gids ;
	my $print_time = 0;
	warn "start";
		foreach my $g (@$list) {
		foreach my $vid ( keys %{ $g->{all_variants} } ) {
			push(@$gids,$g->{chr_name}."!".$g->{all_vector_ids}->{$vid});
	
		}
		}

		$final_polyviewer_all->prepare($gids);
		warn "end prpare ".abs(time-$rtime);
	foreach my $g (@$list) {
		#last if $xp > 100;
					
		$xp++;
#		warn "$xp /".scalar(@$list);
		#warn $xp;
		print "*" if $xp % 10 == 0 && $id == 1;
		my $t1 = time;
		my ( $n, $cname ) = split( "_", $g->{id} );
		my $chr = $project->getChromosome($cname);
		if ( $current ne $cname && $current ) {

			#$project->buffer->close_lmdb();
		}

		$cname = $current;    #unless $cname;
		next unless scalar( keys %{ $g->{all_variants} } );

		my $out;

		$out .= $cgi->start_div(
			{
				class => "panel panel-primary ",
				style =>
"border-color:white;-webkit-border-radius: 3px;-moz-border-radius: 3px;border-radius: 3px;border: 1px solid black;"
			}
		);
		$out .= $cgi->start_div(
			{
				class => "panel-heading panel-face panel-grey",
				style => $print_html->bgcolor.";min-height:13px;max-height:13px;padding:10px;border:0px"
			}
		);
		my $panel_id = "panel_" . $g->{uid};
		$out .= update_variant_editor::panel_gene( $g, $panel_id, $project->name,$patient );
		$out .= $cgi->end_div();

		$out .= "\n";
		$out .= $cgi->start_div(
			{
				class => "panel-body panel-collapse collapse ",
				style => "font-size: 09px;font-family:  Verdana;",
				id    => "$panel_id"
			}
		);
		$out .= "\n";

		#$out.="<br>\n";
		$out .= $cgi->start_table(
			{
				class =>
"table table-striped table-condensed table-bordered table-hover table-mybordered",
				style =>
"vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;line-height: 25px;min-height: 25px;height: 25px;box-shadow: 3px 3px 5px #555;"
			}
		);
		$out .= "\n";
		$out .= $out_header;
		$time_gene += abs(time -$t1);
		my $debug;
		$debug=1 if $g->{external_name} eq "COL7A1";
			my $ttt =  time;
		foreach my $vid ( keys %{ $g->{all_variants} } ) {
			my $vector_id = $g->{all_vector_ids}->{$vid};
			my $stime = time;
			my $t =  time;
			my $vp = $final_polyviewer_all->get($g->{chr_name}."!".$g->{all_vector_ids}->{$vid},1);
			$vp->{transcripts} = $vp->{hgenes}->{$g->{id}}->{tr};
			$vp->{text_caller} =  $vp->{patients_calling}->{$patient->id}->{array_text_calling};
			bless $vp , 'PolyviewerVariant';
			$time_decode += abs(time -$t) ;
			#next;
			$vp->gene($g);
			
		
			
			unless ($vp) {
				confess();
#				$vh = hvariant::hash_variant_2( $patient, $vid );
			}

			if ( exists $hno->{$vid} ) {
				$vp->{composite} = 1;
			}
			
			my $style = {};
			$style = { style => "background-color: #DAEEED;opacity:0.5" } if exists $g->{all_variants}->{$vid}->{added};



		
			my $hpatients;
		
		
			#my $is_gnomad = exists $v->{value}->{ac};
			

		
			
			$print_html->variant($vp);
			##############
			# VARSOME CELL
			###############
			my $t1 = shift(@headers);
			$out .= $cgi->td( $style, $print_html->varsome() );
			
			$out .= "\n";
			##############
			# IGV CELL
			###############

			$t = shift(@headers);

			#write locus
			$out .= $cgi->td( $style, $print_html->igv);

			##############
			# ALAMUT CELL
			###############

			$t = shift(@headers);

			$out .= $cgi->td( $style,$print_html->alamut);
			$out .= "\n";

			##############
			# NAME CELL
			#
			###############
			$t = shift(@headers);

			#$name =  $v->{var_name} if exists $v->{var_name};

			$out .= $cgi->td($style,$print_html->var_name());

			$out .= "\n";
			##############
			# CELL CALLING INFOS
			###############

	
			$out .= $cgi->td( $style, $print_html->calling()) ;

			$out .= "\n";
			
			
		
			
			$t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->gnomad() );
			$out .= "\n";


			$t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->dejavu() );
			$out .= "\n";
			$t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->validations );
			
			$t = shift(@headers);
		
			$out .= "\n";
			$out .= $cgi->td( $style, $print_html->transcripts() );
			$out .= "\n";

			$out .= $cgi->td( $style, $print_html->validation_select() );
			$out .= $cgi->end_Tr();
			$out .= "\n";


		}

		$out .= $cgi->end_table();
		$out .= "\n";
		$out .= $cgi->end_div();
		$out .= $cgi->end_div();
		$out .= "<br>\n";
		$g->{out} = $out;
		$print_time += abs(time -$ttt);

	}
	$project->buffer->close_lmdb();
	
	warn "\t\t decode : ".$time_decode." gene : ".$time_gene." print_time : $print_time  total : ".abs(time -$total_time);
	return ( $list, abs(time -$total_time) );

}





sub new_refine_heterozygote_composite_score_lmdb {
	my ( $project,$print_html,$list, $id ) = @_;
	my $version_db = undef;
	my $bgcolor = qq{background-color:#607D8B};
	
	my $patient= $print_html->patient;
	my $cgi = $print_html->cgi;
	
	my $out_header;
	$out_header .= $cgi->start_Tr( { style => "background-color:#E9DEFF" } );
	foreach my $h (@headers) {
		$out_header .= $cgi->th($h);
	}
	$out_header .= $cgi->th("validations");
	$out_header .= $cgi->th("transcripts");
	$out_header .= $cgi->end_Tr();
	my $mother = $patient->getFamily->getMother();
	my $father = $patient->getFamily->getFather();
	my $hno;
	my $tsum = 0;
	my $t    = time;
	my $xp   = 0;

	my $total_time = 0;
	my $ztotal     = 0;
	#$t = time;
	my $current;
	my $rtime = 0;
	foreach my $g (@$list) {
		#last if $xp > 100;
		$xp++;

		#warn $xp;
		print "*" if $xp % 10 == 0 && $id == 1;
		my ( $n, $cname ) = split( "_", $g->{id} );
		my $chr = $project->getChromosome($cname);
		if ( $current ne $cname && $current ) {

			#$project->buffer->close_lmdb();
		}

		$cname = $current;    #unless $cname;
		next unless scalar( keys %{ $g->{all_variants} } );

		#my $no       = $chr->lmdb_polyviewer_variants( $patient, "r" );

		my $noV = $chr->get_lmdb_variations("r");
		my $no  = $chr->lmdb_polyviewer_variants( $patient, "r" );
		my $out;

		$out .= $cgi->start_div(
			{
				class => "panel panel-primary ",
				style =>
"border-color:white;-webkit-border-radius: 3px;-moz-border-radius: 3px;border-radius: 3px;border: 1px solid black;"
			}
		);
		$out .= $cgi->start_div(
			{
				class => "panel-heading panel-face panel-grey",
				style => "$bgcolor;min-height:13px;max-height:13px;padding:10px;border:0px"
			}
		);
		my $panel_id = "panel_" . $g->{uid};
		$out .=
		  update_variant_editor::panel_gene( $g, $panel_id, $project->name,
			$patient );
		$out .= $cgi->end_div();

		$out .= "\n";
		$out .= $cgi->start_div(
			{
				class => "panel-body panel-collapse collapse ",
				style => "font-size: 09px;font-family:  Verdana;",
				id    => "$panel_id"
			}
		);
		$out .= "\n";

		#$out.="<br>\n";
		$out .= $cgi->start_table(
			{
				class =>
"table table-striped table-condensed table-bordered table-hover table-mybordered",
				style =>
"vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;line-height: 25px;min-height: 25px;height: 25px;box-shadow: 3px 3px 5px #555;"
			}
		);
		$out .= "\n";
		$out .= $out_header;
		my $debug;
		$debug=1 if $g->{external_name} eq "COL7A1";
		
		foreach my $vid ( keys %{ $g->{all_variants} } ) {
		#	my $v  = $noV->get($vid);
			my $vh = $no->get($vid);
			if($debug){
			my $vv = $noV->get($vid);
			$vv->{project} = $project;
			$vv->{buffer} = $project->buffer;
				warn "\t".$vv->id;
				warn "\t ==> ".Dumper($vv->dejaVuInfosForDiag2);
				
			}			
			my $vp =  PolyviewerVariant->new();
			
			$vp->setOldVariant($vh,$project,$patient,$g,$version_db,$debug);
			#warn Dumper $vp if $debug;
			#$vp->setLmdbVariant($vh,$project,$g,$patient);
			$print_html->variant($vp);
			


			my $ttime = time;
			$rtime += abs( time - $ttime );
			unless ($vh) {
				confess();
				$vh = hvariant::hash_variant_2( $patient, $vid );
			}

			if ( exists $hno->{$vid} ) {
				$vh->{composite} = 1;

			}
			my $style = {};
			$style = { style => "background-color: #DAEEED;opacity:0.5" } if exists $g->{all_variants}->{$vid}->{added};

			my $t = time;

			$total_time += abs( time - $t );

		
			my $hpatients;
		
		
			#my $is_gnomad = exists $v->{value}->{ac};
			

			my @headers = (
				"varsome", "igv",    "alamut", "var_name",
				"trio",    "gnomad", "deja_vu"
			);
			
			
			##############
			# VARSOME CELL
			###############
			my $t1 = shift(@headers);
			$out .= $cgi->td( $style, $print_html->varsome() );
			$out .= "\n";

			##############
			# IGV CELL
			###############

			$t = shift(@headers);

			#write locus
			$out .= $cgi->td( $style, $print_html->igv);

			##############
			# ALAMUT CELL
			###############

			$t = shift(@headers);

			$out .= $cgi->td( $style,$print_html->alamut);
			$out .= "\n";

			##############
			# NAME CELL
			#
			###############
			$t = shift(@headers);

			#$name =  $v->{var_name} if exists $v->{var_name};

			$out .= $cgi->td($style,$print_html->var_name());

			$out .= "\n";
			##############
			# CELL CALLING INFOS
			###############

			
				$out .= $cgi->td( $style, $print_html->calling()) ;

			$out .= "\n";
			
			
			
			
			$t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->gnomad() );
			$out .= "\n";


			$t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->dejavu() );
			$out .= "\n";
			$t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->validations );
			
			$t = shift(@headers);
			$out .= "\n";
			$out .= $cgi->td( $style, $print_html->transcripts() );
			$out .= "\n";

			$out .= $cgi->td( $style, $print_html->validation_select() );
			$out .= $cgi->end_Tr();
			$out .= "\n";

		}

		$out .= $cgi->end_table();
		$out .= "\n";
		$out .= $cgi->end_div();
		$out .= $cgi->end_div();
		$out .= "<br>\n";
		$g->{out} = $out;

	}
	$project->buffer->close_lmdb();
	return ( $list, $total_time );

}


#sub new_refine_heterozygote_composite_score_rocks1 {
#	my ( $project,$print_html,$list, $id ) = @_;
#	
#	my $cgi = $print_html->cgi;
#	my $patient= $print_html->patient;
#	
#	my $noP =  GenBoNoSqlRocksPolyviewerVariant->new(dir=>"/data-beegfs/test-cache/polyviewer/",mode=>"r",name=>$project->name.".polyviewer_variant");
#	my $out_header;
#	$out_header= $print_html->print_header("background-color:#E9DEFF");
#	my $mother = $patient->getFamily->getMother();
#	my $father = $patient->getFamily->getFather();
#	my $hno;
#	my $tsum = 0;
#	my $t    = time;
#	my $xp   = 0;
#
#	my $total_time = 0;
#	my $ztotal     = 0;
#	#$t = time;
#	my $current;
#	my $rtime = time;
#	
#
#	
#	foreach my $g (@$list) {
#		#last if $xp > 100;
#		$xp++;
##		warn "$xp /".scalar(@$list);
#		#warn $xp;
#		print "*" if $xp % 10 == 0 && $id == 1;
#		my ( $n, $cname ) = split( "_", $g->{id} );
#		my $chr = $project->getChromosome($cname);
#		if ( $current ne $cname && $current ) {
#
#			#$project->buffer->close_lmdb();
#		}
#
#		$cname = $current;    #unless $cname;
#		next unless scalar( keys %{ $g->{all_variants} } );
#
#		my $no  =  $chr->rocks_polyviewer_variants("r");
#		my $out;
#
#		$out .= $cgi->start_div(
#			{
#				class => "panel panel-primary ",
#				style =>
#"border-color:white;-webkit-border-radius: 3px;-moz-border-radius: 3px;border-radius: 3px;border: 1px solid black;"
#			}
#		);
#		$out .= $cgi->start_div(
#			{
#				class => "panel-heading panel-face panel-grey",
#				style => $print_html->bgcolor.";min-height:13px;max-height:13px;padding:10px;border:0px"
#			}
#		);
#		my $panel_id = "panel_" . $g->{uid};
#		$out .= update_variant_editor::panel_gene( $g, $panel_id, $project->name,$patient );
#		$out .= $cgi->end_div();
#
#		$out .= "\n";
#		$out .= $cgi->start_div(
#			{
#				class => "panel-body panel-collapse collapse ",
#				style => "font-size: 09px;font-family:  Verdana;",
#				id    => "$panel_id"
#			}
#		);
#		$out .= "\n";
#
#		#$out.="<br>\n";
#		$out .= $cgi->start_table(
#			{
#				class =>
#"table table-striped table-condensed table-bordered table-hover table-mybordered",
#				style =>
#"vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;line-height: 25px;min-height: 25px;height: 25px;box-shadow: 3px 3px 5px #555;"
#			}
#		);
#		$out .= "\n";
#		$out .= $out_header;
#		my $debug;
#		$debug=1 if $g->{external_name} eq "COL7A1";
#		foreach my $vid ( keys %{ $g->{all_variants} } ) {
#			my $vector_id = $g->{all_vector_ids}->{$vid};
#			my $stime = time;
#			my $vp =  return_vp($project,$chr,$patient,$vid,$g);
#			#my $vp = $no->test($vid,$patient->id,$g->{id});
#			
#			#next;
#			$vp->gene($g);
#			
#
#			
#			unless ($vp) {
#				confess();
##				$vh = hvariant::hash_variant_2( $patient, $vid );
#			}
#
#			if ( exists $hno->{$vid} ) {
#				$vp->{composite} = 1;
#
#			}
#			my $style = {};
#			$style = { style => "background-color: #DAEEED;opacity:0.5" } if exists $g->{all_variants}->{$vid}->{added};
#
#			my $t = time;
#
#			$total_time += abs( time - $t );
#
#		
#			my $hpatients;
#		
#		
#			#my $is_gnomad = exists $v->{value}->{ac};
#			
#
#			my @headers = (
#				"varsome", "igv",    "alamut", "var_name",
#				"trio",    "gnomad", "deja_vu"
#			);
#			
#			$print_html->variant($vp);
#			##############
#			# VARSOME CELL
#			###############
#			my $t1 = shift(@headers);
#			warn $vid unless $vp->id;
#			$out .= $cgi->td( $style, $print_html->varsome() );
#			
#			$out .= "\n";
#
#			##############
#			# IGV CELL
#			###############
#
#			$t = shift(@headers);
#
#			#write locus
#			$out .= $cgi->td( $style, $print_html->igv);
#
#			##############
#			# ALAMUT CELL
#			###############
#
#			$t = shift(@headers);
#
#			$out .= $cgi->td( $style,$print_html->alamut);
#			$out .= "\n";
#
#			##############
#			# NAME CELL
#			#
#			###############
#			$t = shift(@headers);
#
#			#$name =  $v->{var_name} if exists $v->{var_name};
#
#			$out .= $cgi->td($style,$print_html->var_name());
#
#			$out .= "\n";
#			##############
#			# CELL CALLING INFOS
#			###############
#
#			
#				$out .= $cgi->td( $style, $print_html->calling()) ;
#
#			$out .= "\n";
#			
#			
#			
#			
#			$t = shift(@headers);
#			$out .= $cgi->td( $style, $print_html->gnomad() );
#			$out .= "\n";
#
#
#			$t = shift(@headers);
#			$out .= $cgi->td( $style, $print_html->dejavu() );
#			$out .= "\n";
#			$t = shift(@headers);
#			$out .= $cgi->td( $style, $print_html->validations );
#			
#			$t = shift(@headers);
#			$out .= "\n";
#			$out .= $cgi->td( $style, $print_html->transcripts() );
#			$out .= "\n";
#
#			$out .= $cgi->td( $style, $print_html->validation_select() );
#			$out .= $cgi->end_Tr();
#			$out .= "\n";
#
#		}
#
#		$out .= $cgi->end_table();
#		$out .= "\n";
#		$out .= $cgi->end_div();
#		$out .= $cgi->end_div();
#		$out .= "<br>\n";
#		$g->{out} = $out;
#
#	}
#	$project->buffer->close_lmdb();
#	return ( $list, $total_time );
#
#}
1;
