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
	my ( $project,$print_html,$list, $id,$rocksdb_pv) = @_;
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
	#$rocksdb_pv->print_html($print_html);
	eval {
		$rocksdb_pv->load_polyviewer_variant();
		};
	if ($@){
		return undef;
	};
	foreach my $g (@$list) {
		$xp++;
		print "@"  if $xp % 3 == 0 ;
		print "*" if $xp % 10 == 0 && $id == 1;
		last if $xp > 500;
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
		my $color_validation = "grey";
	
		foreach my $vid ( sort {$a cmp $b} keys %{ $g->{all_variants} } ) {
			my $id= $g->{chr_name}."!".$g->{all_vector_ids}->{$vid};
			my $vp = $rocksdb_pv->get_polyviewer_variant($id,1);
			
			#next unless $vp;
			$vp->{transcripts} = $vp->{hgenes}->{$g->{id}}->{tr};
			$vp->gene($g);
			my $opacity;
			#warn $id if exists $g->{all_variants}->{$vid}->{added};
			$opacity = 1 if exists $g->{all_variants}->{$vid}->{added};
			
			unless ($vp) {
				confess();
			}

			if ( exists $hno->{$vid} ) {
				$vp->{composite} = 1;
			}
			$out .= print_line_variant($vp,$print_html,$opacity);
			if  ($vp->clinvar_value >0 or $vp->hgmd_value> 0 or $vp->{local_value}){
				$color_validation= "warning";
			}
			elsif  ( exists $g->{pathogenic}){
				$color_validation= "danger";
			}
			
		

		}
		$out =~ s/@@@/$color_validation/;
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
			
			my $dropdown = qq{
			<div class="dropdown">
  <button class="btn btn-primary btn-xs dropdown-toggle " type="button" id="dropdownMenuButton" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false" style="font-size:10px;background-color:#C67FAE">
  $icon
  </button>
  <div class="dropdown-menu" aria-labelledby="dropdownMenuButton" style="font-size:12px;background-color:beige;color:black">
  }.
  "<li>".$print_html->mobidetails()."</li>".
  	"<li>".$print_html->gnomadurl()."</li>".
			"<li>".$print_html->alamuturl()."</li>".
			"<li>".$print_html->varsome()."</li>"
  
  .qq{</div>
</div>
};


			my $t1 = shift(@headers);
			$out .= $cgi->td( $style, $dropdown );
		
			
			
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


sub print_line_variant_all_patients {
	my ($list_vp,$list_print_html,$opacity) = @_;
	my $cgi = $list_print_html->[0]->cgi;
	my $style = {};
	$style = { style => "background-color: #DAEEED;opacity:0.5" }  if $opacity;
	my $out;
	my $hpatients;
	my $i = 0;
	foreach my $print_html (@$list_print_html) {
		$print_html->variant($list_vp->[$i]);
		$i++;
	}
	my $icon = qq{<img width="32" height="32" src="https://img.icons8.com/external-gliphyline-royyan-wijaya/32/external-laptop-laptop-collection-glyphyline-gliphyline-royyan-wijaya-15.png" alt="external-laptop-laptop-collection-glyphyline-gliphyline-royyan-wijaya-15"/>};
	$icon   = qq{<img width="24" height="24" src="https://img.icons8.com/external-tal-revivo-filled-tal-revivo/24/external-live-preview-of-a-smart-class-education-school-filled-tal-revivo.png" alt="external-live-preview-of-a-smart-class-education-school-filled-tal-revivo"/>};
	my $dropdown = qq{
		<div class="dropdown">
		<button class="btn btn-primary btn-xs dropdown-toggle " type="button" id="dropdownMenuButton" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false" style="font-size:10px;background-color:#C67FAE">
		$icon
		</button>
		<div class="dropdown-menu" aria-labelledby="dropdownMenuButton" style="font-size:12px;background-color:beige;color:black">
	};
	$dropdown .= "<li>".$list_print_html->[0]->mobidetails()."</li>";
  	$dropdown .= "<li>".$list_print_html->[0]->gnomadurl()."</li>";
	$dropdown .= "<li>".$list_print_html->[0]->alamuturl()."</li>";
	$dropdown .= "<li>".$list_print_html->[0]->varsome()."</li>";
	$dropdown .= qq{</div></div>};
	my $t1 = shift(@headers);
	$out .= $cgi->td( $style, $dropdown );
	$out .= "\n";
	my $t = shift(@headers);
	$out .= $cgi->td( $style, $list_print_html->[0]->igv);
	$t = shift(@headers);
	$out .= $cgi->td($style,$list_print_html->[0]->var_name());
	$out .= "\n";
	
	my @l_html_calling;
	foreach my $print_html (@$list_print_html) {
		#push(@l_html_calling, "<tr style='padding:3px;'><td colspan='6'><center><span style='font-size:13px;'>fam ".$print_html->patient->getFamily->name()."</span></center></td></tr>");
		push(@l_html_calling, "<tr style='padding:6px;'><td style='padding-right:10px;'><center><b>".$print_html->patient->getFamily->name()."</b></center></td><td>".$print_html->calling()."</td></tr>");
		push(@l_html_calling, "<tr style='padding:6px;'><td><br></td><td><br></td></tr>");
	}
	
	my $out_calling = qq{<center><table style='width:98%;'>};
	$out_calling .= join('', @l_html_calling);
	$out_calling .= qq{</table></center>};
	
	$out .= $cgi->td( $style, $out_calling) ;
	$out .= "\n";
	$t = shift(@headers);
	$out .= $cgi->td( $style, $list_print_html->[0]->gnomad() );
	$out .= "\n";
	$t = shift(@headers);
	$out .= $cgi->td( $style, $list_print_html->[0]->dejavu() );
	$out .= "\n";
	$t = shift(@headers);
	$out .= $cgi->td( $style, $list_print_html->[0]->validations );
	$t = shift(@headers);
	$out .= "\n";
	$out .= $cgi->td( $style, $list_print_html->[0]->transcripts() );
	$out .= $cgi->end_Tr();
	$out .= "\n";
}



sub print_results_by_genes_for_patients {
	my ($project, $list, $h_rocksdb_pv_patients,$h_genes_search) = @_;
	
	my ($print_html, $patient);
	foreach my $project_name (keys %$h_rocksdb_pv_patients) {
		foreach my $patient_name (keys %{$h_rocksdb_pv_patients->{$project_name}}) {
			$h_rocksdb_pv_patients->{$project_name}->{$patient_name}->load_polyviewer_variant();
			$h_rocksdb_pv_patients->{$project_name}->{$patient_name}->print_html();
			$print_html = $h_rocksdb_pv_patients->{$project_name}->{$patient_name}->print_html() if not $print_html;
			$patient = $h_rocksdb_pv_patients->{$project_name}->{$patient_name}->patient() if not $patient;
		}
	}
	
	my $cgi = $print_html->cgi;
	my $diro = $project->rocks_directory();
	my $out_header;
	$out_header = $print_html->print_header("background-color:aliceblue;color:black");
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
	
	foreach my $g (sort {$a->{external_name} <=> $b->{external_name}} @$list) {
		my $t1 = time;
		my ( $n, $cname ) = split( "_", $g->{id} );
		my $chr = $project->getChromosome($cname);
		$cname = $current;
		next unless scalar( keys %{ $g->{all_variants} } );

		my $out;
		$out .= $cgi->start_div(
			{
				class => "panel panel-primary",
				style => "border-color:white;width:100%;margin-right:5px;"
			}
		);
		
		my ($opacity_html, $opacity);
		if (not exists $h_genes_search->{$g->{id}} and not exists $h_genes_search->{'ALL'}) {
			$opacity_html = 'opacity:0.5;';
			$opacity = 1;
		}
		$out .= $cgi->start_div(
			{
				class => "panel-heading panel-face panel-grey",
				style => "background-color:#607D8B;height:45px;padding:10px;border:0px;$opacity_html"
			}
		);
		my $panel_id = "panel_" . $g->{uid};
		
		my $html_panel_gene = update_variant_editor::panel_gene($g, $panel_id, $project->name, $patient);
		$out .= $html_panel_gene;
		$out .= $cgi->end_div();

		$out .= "\n";
		$out .= $cgi->start_div(
			{
				class => "panel-body panel-collapse collapse",
				style => "font-size: 09px;font-family: Verdana;",
				id    => "$panel_id"
			}
		);
		$out .= "\n";
		$out .= $cgi->start_table(
			{
				class => "table table-striped table-condensed table-bordered table-hover table-mybordered",
				style => "vertical-align:middle;text-align: center;font-size: 8px;font-family:Verdana;line-height: 25px;min-height: 25px;height: 25px;box-shadow: 3px 3px 5px #555;"
			}
		);
		$out .= "\n";
		$out .= $out_header;
		$time_gene += abs(time -$t1);
		my $ttt =  time;
		my $color_validation = "grey";
		foreach my $vid (sort keys %{ $g->{all_variants} } ) {
			my $id= $g->{chr_name}."!".$g->{all_vector_ids}->{$vid};
			my ($h_fam_found, $rocksdb_pv);
			foreach my $project_name (keys %$h_rocksdb_pv_patients) {
				foreach my $patient_name (keys %{$h_rocksdb_pv_patients->{$project_name}}) {
					next if not ($h_rocksdb_pv_patients->{$project_name}->{$patient_name}->has_index($g->{chr_name}, $vid));
					$h_rocksdb_pv_patients->{$project_name}->{$patient_name}->get_polyviewer_variant($id,1);
					$rocksdb_pv = $h_rocksdb_pv_patients->{$project_name}->{$patient_name};
					my $patient = $rocksdb_pv->print_html->patient();
					if ($patient->isChild()) {
						$h_fam_found->{$project_name}->{$patient->getFamily->name}->{children}->{$patient->name} = undef;
					}
					else {
						$h_fam_found->{$project_name}->{$patient->getFamily->name}->{parents}->{$patient->name} = undef;
					}
				}
			}
			my (@list_print_html, @list_vp);
			foreach my $project_name (sort keys %$h_fam_found) {
				foreach my $fam_name (sort keys %{$h_fam_found->{$project_name}}) {
					if (exists $h_fam_found->{$project_name}->{$fam_name}->{children}) {
						my @lChilds = sort keys %{$h_fam_found->{$project_name}->{$fam_name}->{children}};
						push(@list_vp, $h_rocksdb_pv_patients->{$project_name}->{$lChilds[0]}->get_polyviewer_variant($id,1));
						push(@list_print_html, $h_rocksdb_pv_patients->{$project_name}->{$lChilds[0]}->print_html());
					}
					else {
						foreach my $parent_name (sort keys %{$h_fam_found->{$project_name}->{$fam_name}->{parents}}) {
							push(@list_vp, $h_rocksdb_pv_patients->{$project_name}->{$parent_name}->get_polyviewer_variant($id,1));
							push(@list_print_html, $h_rocksdb_pv_patients->{$project_name}->{$parent_name}->print_html());
						}
					}
				}
			}
			
			foreach my $vp (@list_vp) {
				$vp->{transcripts} = $vp->{hgenes}->{$g->{id}}->{tr};
				$vp->{gene} = $g;
				$opacity = 1 if exists $g->{all_variants}->{$vid}->{added};
				unless ($vp) {
					confess();
				}
				if ( exists $hno->{$vid} ) {
					$vp->{composite} = 1;
				}
			}
			$out .= print_line_variant_all_patients(\@list_vp, \@list_print_html, $opacity);
			if ($list_vp[0]->clinvar_value >0 or $list_vp[0]->hgmd_value> 0 or $list_vp[0]->{local_value}){
				$color_validation= "warning";
			}
			elsif (exists $g->{pathogenic}){
				$color_validation= "danger";
			}
		}
		$out =~ s/@@@/$color_validation/;
		$out .= $cgi->end_table();
		$out .= $cgi->end_div();
		$out .= $cgi->end_div();
		$g->{out} = $out;
		$print_time += abs(time -$ttt);
	}
	$project->buffer->close_lmdb();
	warn "\t\t decode : ".$time_decode." gene : ".$time_gene." print_time : $print_time  total : ".abs(time -$total_time);

	return ( $list, abs(time -$total_time) );
}

1;
