package variations_editor_methods;
use strict;
use Data::Dumper;
use Time::HiRes qw( time);

#sub refine_heterozygote_composite {
#	my ($project,$print_html,$list, $id, ) = @_;
#		if ($project->isRocks){
#		return new_refine_heterozygote_composite_score_rocks($project, $print_html,$list, $id);
#		}
#		else {
#			return new_refine_heterozygote_composite_score_lmdb($project, $print_html,$list, $id);
#		}
#
#}

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
				"mobidetails","varsome", "igv",    "alamut", "var_name",
				"trio",    "gnomad", "deja_vu"
			);

sub refine_heterozygote_composite {
	my ( $project,$print_html,$list, $id,$final_polyviewer_all) = @_;
	my $cgi = $print_html->cgi;
	my $patient= $print_html->patient;

	my $diro = $project->rocks_directory();
	my $out_header;
	$out_header= $print_html->print_header("background-color:aliceblue;color:black");
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
#	warn "start";
		foreach my $g (@$list) {
		foreach my $vid ( keys %{ $g->{all_variants} } ) {
			push(@$gids,$g->{chr_name}."!".$g->{all_vector_ids}->{$vid});
	
		}
		}
	warn "prepare ".$final_polyviewer_all->dir;
	;
	$final_polyviewer_all->prepare($gids);
	warn "end prepare";
	foreach my $g (@$list) {
					
		$xp++;
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
			my $vp = $final_polyviewer_all->get($g->{chr_name}."!".$g->{all_vector_ids}->{$vid},1);
			$vp->{transcripts} = $vp->{hgenes}->{$g->{id}}->{tr};
			$vp->{spliceAI} = $vp->{hgenes}->{$g->{id}}->{sc}->{spliceAI};
			
			#$vp->{hgenes}->{$g->{id}}->{spliceAI};
			$vp->{spliceAI_cat} = $vp->{hgenes}->{$g->{id}}->{sc}->{spliceAI_cat};
			$vp->{text_caller} =  $vp->{patients_calling}->{$patient->id}->{array_text_calling};
			bless $vp , 'PolyviewerVariant';
						$vp->gene($g);
			my $opacity;
			$opacity = 1 if exists $g->{all_variants}->{$vid}->{added};
			
			unless ($vp) {
				confess();
			}

			if ( exists $hno->{$vid} ) {
				$vp->{composite} = 1;
			}
			
			$out .= print_line_variant($vp,$print_html,$opacity)

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
	delete $final_polyviewer_all->{rocks};
	warn "\t\t decode : ".$time_decode." gene : ".$time_gene." print_time : $print_time  total : ".abs(time -$total_time);

	return ( $list, abs(time -$total_time) );

}

sub print_line_variant {
	my ($vp,$print_html,$opacity) = @_;
			
		 my $cgi = $print_html->cgi;
			
			my $style = {};
			$style = { style => "background-color: #DAEEED;opacity:0.5" }  if $opacity;

	my $out;

		
			my $hpatients;
		
		
			#my $is_gnomad = exists $v->{value}->{ac};
			

			
			$print_html->variant($vp);

			my $icon = qq{<img width="32" height="32" src="https://img.icons8.com/external-gliphyline-royyan-wijaya/32/external-laptop-laptop-collection-glyphyline-gliphyline-royyan-wijaya-15.png" alt="external-laptop-laptop-collection-glyphyline-gliphyline-royyan-wijaya-15"/>};
			 $icon   = qq{<img width="24" height="24" src="https://img.icons8.com/external-tal-revivo-filled-tal-revivo/24/external-live-preview-of-a-smart-class-education-school-filled-tal-revivo.png" alt="external-live-preview-of-a-smart-class-education-school-filled-tal-revivo"/>};
			#$icon = qq{<img width="28" height="28" src="https://img.icons8.com/external-bearicons-blue-bearicons/28/external-Clipboard-clipboards-and-notepads-bearicons-blue-bearicons-30.png" alt="external-Clipboard-clipboards-and-notepads-bearicons-blue-bearicons-30"/>};
			my $div =qq{
			<button class="external-viewer-button">
			$icon
			</button>
			<div class="external-viewer-items">
			}.
			$print_html->mobidetails().
			$print_html->gnomadurl().
			$print_html->alamuturl().
			$print_html->varsome().
        	qq{
     	 		</div>
			};
				my $t1 = shift(@headers);
			$out .= $cgi->td( $style, $div );
		
			
			
			$out .= "\n";
			##############
			# IGV CELL
			###############

			my $t = shift(@headers);

			#write locus
			$out .= $cgi->td( $style, $print_html->igv);

			##############
			
			
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




1;
