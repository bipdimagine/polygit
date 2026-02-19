package dejavu_variants;

use strict;
use Moo;
use Bit::Vector;
use Bit::Vector::Overload;
use Set::IntSpan::Fast;
use Data::Dumper;
use List::MoreUtils qw(natatime);
use Compress::Snappy;
use MIME::Base64;
use Storable qw(store retrieve freeze dclone thaw);
use JSON;


#use polyviewer_html;

#use html_polygenescout;

my @headers = ("Viewer", "var_name", "trio", "gnomad", "deja_vu");


has fork => (
	is		=> 'rw',
	lazy    => 1,
	default => 1
);

has buffer => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $buffer = new GBuffer;
		return $buffer;
	}
);

has project => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->buffer->newProject( -name => $self->project_name() );
	}
);

has project_name => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $buffer = new GBuffer;
		return $buffer->getRandomProjectName();
	}
);

has hash_projects_ids_names => (
	is		=> 'rw',
	lazy    => 1,
);


has user_name => (
	is		=> 'rw',
	lazy    => 1,
);

has pwd => (
	is		=> 'rw',
	lazy    => 1,
);

has min_promoter_ai => (
	is		=> 'rw',
	lazy    => 1,
);

has max_gnomad_ac => (
	is		=> 'rw',
	lazy    => 1,
);

has max_gnomad_ac_ho => (
	is		=> 'rw',
	lazy    => 1,
);

has hash_users_projects => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my ($h_projects, @list_hash);
		if ($self->buffer->getQuery->isUserMagic($self->user_name(), $self->pwd())) {
			@list_hash = @{$self->buffer->getQuery()->getAllProjects()};
		}
		else {
			@list_hash = @{$self->buffer->getQuery()->getProjectListForUser($self->user_name(), $self->pwd())};
		}
		foreach my $hash (@list_hash) {
			my $proj_name = $hash->{name};
			next unless ($proj_name =~ /NGS20/);
			$h_projects->{$proj_name}->{description} = $hash->{description};
			$h_projects->{$proj_name}->{name} = $proj_name;
			$h_projects->{$proj_name}->{id} = $hash->{id};
			$h_projects->{$hash->{id}} = $h_projects->{$proj_name};
		}
		return $h_projects;
	}
);

has hash_lift_variants => (
	is		=> 'rw',
	lazy    => 1
);



sub check_variants_from_gene {
	my ($self, $h_dv_rocks_ids) = @_;
	print '!';
	my $ii = 0;
	my ($h_dv_var_ids, @lVarIds, @lVar, $h_dv_var_ids_erros_lift);
	my ($total, $total_pass);
	my @list_var;
	foreach my $chr_id (sort keys %{$h_dv_rocks_ids}) {
		print '.';
		my $chr = $self->project->getChromosome($chr_id);
		foreach my $rocks_id (sort keys %{$h_dv_rocks_ids->{$chr_id}}) {
			$ii++;
			if ($ii == 50) {
				print '.';
				$ii = 0;
			}
			my $var_id = $chr->transform_rocksid_to_varid($rocks_id);
			my $var = $self->project->_newVariant($var_id);
			if ($self->max_gnomad_ac()) {
				next if defined $var->getGnomadAC() and $var->getGnomadAC() > $self->max_gnomad_ac();
			}
			if ($self->max_gnomad_ac_ho()) {
				next if defined $var->getGnomadHO() and $var->getGnomadHO() > $self->max_gnomad_ac_ho();
			}
			push(@list_var, $var);
		}
	}
	
	my ($hGenes, $hVariants);
	my $pm = new Parallel::ForkManager($self->fork());
	my $nb_errors=0;
	$pm->run_on_finish(
		sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data) = @_;
			if ($self->fork == 1) {
				$hGenes = $data->{genes};
				$hVariants = $data->{variants};
			}
			else {
				foreach my $gid (keys %{$data->{lift}}) {
					$self->{hash_lift_variants}->{$gid} = $data->{lift}->{$gid};
				}
				foreach my $gene_id (keys %{$data->{genes}}) {
					$hGenes->{$gene_id} = $data->{genes}->{$gene_id};
				}
				foreach my $var_id (keys %{$data->{variants}}) {
					$hVariants->{$var_id} = $data->{variants}->{$var_id};
				}
			}
		}
	);
	
	my $nb = int((scalar(@list_var)+1)/($self->fork));
	$nb = 20 if $nb == 0;
	my $iter = natatime($nb, @list_var);
	$self->project->disconnect();
	while( my @tmp = $iter->() ){
 	 	my $pid = $pm->start and next;
		my $hres;
		foreach my $var (@tmp) {
			$ii++;
			if ($ii == 50) {
				print '.';
				$ii = 0;
			}
			my $var_id = $var->id();
			my $rocks_id = $var->rocksdb_id();
			my $chr_id = $var->getChromosome->id();
			
			my @lGenes = @{$var->getGenes()};
			foreach my $gene (@lGenes) {
				if ($self->min_promoter_ai()) {
					my $ok;
					foreach my $tr (@{$gene->getTranscripts()}) {
						my $score = $var->promoterAI_score($tr);
						$ok = 1 if $score and abs($score) >= $self->min_promoter_ai();
						last if $ok;
					}
					next if not $ok;
				}
				$hres->{genes}->{$gene->id()}->{$var_id} = undef;
			}
			if (scalar(@lGenes) == 0) {
				$hres->{genes}->{intronic}->{$var_id} = undef;
			}
			
			# STEP 2 - PolyviewerVariant
#			my $vp = PolyviewerVariant->new();
#			$vp->setLmdbVariant($var);
			
			my $vp = PolyviewerVariant->new();
			$vp->setLmdbVariant($var);
			$vp->{hgenes} = {};
			$vp->{genes_id} = [];
			my $code = 0;
			foreach my $g (@{$var->getGenes}){
				my $h = $vp->set_gene($var,$g);
				$h->{code} = $code;
				$vp->{hgenes}->{$g->id} = $h;
				push(@{$vp->{genes_id}},$g->id);
				$code ++;
			}
			##############
			#	next;
			##############0
			$vp->{hpatients} ={};
			$vp->{patients_id} = [];
			
			$hres->{variants}->{$var_id}->{polyviewer_variant} = $vp;
			
			# STEP 3 - printHtml patient
			my @list_parquet;
			foreach my $project_id (sort keys %{$h_dv_rocks_ids->{$chr_id}->{$rocks_id}}) {
				my $project_name = $self->buffer->getQuery->getProjectNameFromId($project_id);
				next if not exists $self->hash_users_projects->{$project_name};
				my $project_type_cache = $self->{hash_projects_ids_names}->{$project_id}->{type};
				$self->{hash_projects_ids_names}->{$project_id}->{name} = $project_name;
				$self->{hash_projects_ids_names}->{$project_id}->{type} = $project_type_cache;
				my $parquet = $self->buffer->dejavu_parquet_dir().'/'.$project_name.'.'.$project_id.'.parquet';
				push(@list_parquet, "'".$parquet."'") if -e $parquet;
			}
			
			my ($h_projects_patients, $h_gnomadid) = $self->get_from_duckdb_project_patients_infos($var, \@list_parquet);
			foreach my $project_name (keys %{$h_projects_patients}) {
#				my $b = new GBuffer;
#				my $pr = $b->newProject( -name => $project_name);
				foreach my $patient_name (keys %{$h_projects_patients->{$project_name}}) {
#					warn ref($h_projects_patients->{$project_name}->{$patient_name}->{print_html});
					push(@{$hres->{variants}->{$var_id}->{polyviewer_html}}, $h_projects_patients->{$project_name}->{$patient_name}->{print_html});
					my $hh;
					$hh->{project_name} = $project_name;
					$hh->{patient_name} = $patient_name;
					$hh->{ratio} = $h_projects_patients->{$project_name}->{$patient_name}->{ratio};
					$hh->{dp} = $h_projects_patients->{$project_name}->{$patient_name}->{dp};
					$hh->{model} = $h_projects_patients->{$project_name}->{$patient_name}->{model};
					push(@{$hres->{variants}->{$var_id}->{polyviewer_html_details_proj_pat}}, $hh);
#					my $pat = $pr->getPatient($patient_name);
#					my $print_html = polyviewer_html->new( project=>$pr, patient=>$pat,header=>\@headers, bgcolor=>"background-color:#607D8B" );
#					push(@{$hVariants->{$var_id}->{polyviewer_html}}, $print_html);
				}
#				$pr = undef;
#				$b = undef;
			}
			foreach my $gid (keys %$h_gnomadid) {
				$hres->{lift}->{$gid} = $h_gnomadid->{$gid};
			}
		}
	 	$pm->finish(0, $hres);
	} 
	sleep(3); 
	$pm->wait_all_children();
	$self->project->disconnect();
	return ($hGenes, $hVariants);
	
}


sub get_table_project_patients_infos {
	my ($self, $project_name, $hash) = @_;
	my $nb_he = $hash->{he};
	my ($h_infos_patients, $h_tmp_pat);
	my $nb_pat = 0;
	foreach my $pat_id (unpack("w*",decode_base64($hash->{patients}))) {
		$nb_pat++;
		$h_infos_patients->{$nb_pat}->{id} = $pat_id;
		$h_tmp_pat->{$pat_id} = $nb_pat;
	}
	my $i = 0;
	$nb_pat = 1;
	foreach my $info (unpack("w*",decode_base64($hash->{dp_ratios}))) {
		$i++;
		if ($i == 1) { $h_infos_patients->{$nb_pat}->{dp} = $info; }
		elsif ($i == 2) {
			my $ratio = '?';
			$ratio = ($info / $h_infos_patients->{$nb_pat}->{dp}) * 100 if ($h_infos_patients->{$nb_pat}->{dp});
			my $text = 'AC:'.$info;
			if ($ratio and $ratio eq '?') { $text .= ', Ratio: ?'; }
			elsif ($ratio) { $text .= ', Ratio:'.int($ratio); }
			$h_infos_patients->{$nb_pat}->{ratio} = $text;
		}
		elsif ($i == 3) {
    		my $model;
    		if ($info == 1) { $model = 'solo'; }
    		elsif ($info == 2) { $model = 'father'; }
    		elsif ($info == 4) { $model = 'mother'; }
    		elsif ($info == 8) { $model = 'both'; }
    		elsif ($info == 16) { $model = 'is_parent'; }
    		elsif ($info == 32) { $model = 'recessif'; }
    		elsif ($info == 64) { $model = 'dominant'; }
    		elsif ($info == 128) { $model = 'denovo'; }
    		elsif ($info == 256) { $model = 'strict_denovo'; }
    		elsif ($info == 512) { $model = 'error'; }
    		else { $model = 'error2'; }
			$h_infos_patients->{$nb_pat}->{model} = $model;
			
			$i = 0;
			$nb_pat++;
		}
	}
	my $hres;
	my $b = new GBuffer;
	my $p = $b->newProject( -name => $project_name);
	my @lPat = @{$p->getPatients()};
	return undef if scalar(@lPat) == 0;
	my $found_healthy_patient;
	foreach my $pat (@lPat) {
		next if not exists $h_tmp_pat->{$pat->id};
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{name} = $pat->name;
		my $icon = $pat->small_icon();
		$icon =~ s/"/'/g;
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{status} = $icon;
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{sex} = '-';
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{sex} = 'male' if $pat->sex() eq '1';
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{sex} = 'female' if $pat->sex() eq '2';
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{status_txt} = '-';
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{status_txt} = 'healthy' if $pat->status() eq '1';
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{status_txt} = 'ill' if $pat->status() eq '2';
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{family} = $pat->getFamily->name();
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{description} = $p->description();
	#TODO: here
#		my $print_html = polyviewer_html->new( project=>$p, patient=>$pat,header=>\@headers, bgcolor=>"background-color:#607D8B" );
#		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{print_html} = $print_html;
#		if (int($h_tmp_pat->{$pat->id}) <= $nb_he) { $h_infos_patients->{$h_tmp_pat->{$pat->id}}->{heho} = 'He'; }
#		else { $h_infos_patients->{$h_tmp_pat->{$pat->id}}->{heho} = 'Ho'; }
	}
	foreach my $id (keys %{$h_infos_patients}) {
		my $pat_name = $h_infos_patients->{$id}->{name};
		next if not $pat_name;
		next if $h_infos_patients->{$id}->{status_txt} ne 'ill';
		$hres->{$pat_name} = $h_infos_patients->{$id};
	}
	$p = undef;
	$b =undef;
	
#	warn "\n";
#	warn "\n";
#	warn Dumper $hres;
	
	return $hres;
}

sub get_from_duckdb_project_patients_infos {
	my ($self, $var, $list_files) = @_;
	return if scalar(@$list_files) == 0;
	my $sql = "PRAGMA threads=3; SELECT project,chr38,chr19,pos38,pos19,he,allele,patients,dp_ratios FROM read_parquet([".join(', ', @$list_files)."])";
	
	my ($pos38, $alt38) = split('!', $var->rocksdb_id());
	my $find_pos_s = $var->start() - 5;
	my $find_pos_e = $var->start() + 5;
	my ($h_gnomadid, $h_projects_patients);
	
	if ($var->getProject->current_genome_version() eq 'HG38') {
		$sql .= " WHERE chr38='".$var->getChromosome->id()."' and pos38 BETWEEN '".$find_pos_s."' and '".$find_pos_e."';" ;
		
		my $duckdb = $self->buffer->software('duckdb');
		my $cmd = qq{set +H | $duckdb -json -c "$sql"};
		my $json_duckdb = `$cmd`;
		
		if ($json_duckdb) {
			my $decode = decode_json $json_duckdb;
			my $h_by_proj;
			foreach my $h (@$decode) {
				#next if $h->{'chr38'} ne $var->getChromosome->id;
				next if $h->{'pos38'} ne int($pos38);
				my $var_all = $h->{'allele'};
				$var_all =~ s/\+//;
				if (not $var->isDeletion) {
					next if $var_all ne $var->var_allele();
				}
				$h_gnomadid->{$var->gnomad_id()}->{chr19} = $h->{chr19};
				$h_gnomadid->{$var->gnomad_id()}->{pos19} = $h->{pos19};
				$h_gnomadid->{$var->gnomad_id()}->{ref_all} = $var->ref_allele();
				$h_gnomadid->{$var->gnomad_id()}->{var_all} = $var->var_allele();
				my $project_id = $h->{project};
				my $project_name = $self->{hash_projects_ids_names}->{$project_id}->{name};
				$h_by_proj->{$project_name} = $h;
				$h_projects_patients->{$project_name} = $self->get_table_project_patients_infos($project_name, $h);
			}
		}	
	}
	else {
		warn "HG19";
		warn $sql if $self->debug;
		confesss("HG19!!");
	}
	return ($h_projects_patients, $h_gnomadid);
}



sub print_html_gene {
	my ($self, $gene_id, $list_variants, $hVariantsDetails) = @_;
	print '|';
	my $found;
	foreach my $var_id (sort @$list_variants) {
		next if not exists $hVariantsDetails->{$var_id};
		next if not $hVariantsDetails->{$var_id}->{polyviewer_html_details_proj_pat};
		$found = 1;
	}
	return if not $found;
	
	my $gene = $self->project->newGene($gene_id);
	my @l_var_ids = keys %{$hVariantsDetails};
	
	my $pr = $self->project;
	my $pat = $self->project->getPatients->[0];
	my $print_html = polyviewer_html->new( project=>$pr, patient=>$pat,header=>\@headers, bgcolor=>"background-color:#607D8B" );

	my $g;
	$g->{id} = $gene_id;
	$g->{uid} = $gene_id;
	$g->{name} = $gene->external_name();
	$g->{external_name} = $gene->external_name();
	$g->{chr_name} = $gene->getChromosome->id();
	$g->{nb} = scalar(@$list_variants);
	$g->{omim_id} = $gene->omim_id();
	$g->{pLI} = $gene->pLI();
	$g->{phenotypes} = $gene->description();
	$g->{score} = $gene->raw_score();
	$g->{max_score} = $g->{score};
	$g->{variants} = $list_variants;
	
	my $panel_id = "panel_".$gene_id;
	my $out;
	my $bg_color = $print_html->bgcolor;
	$out .= qq{<div class="panel-heading panel-face panel-grey"	style="$bg_color;height:43px;padding:10px;border:0px;width:100%;">};
	$out .= update_variant_editor::panel_gene($g, $panel_id, $self->project->name);
	$out .= qq{</div>};
	$out .= qq{<br>};
	$out .= qq{<div class="panel-body panel-collapse collapse" style="font-size: 09px;font-family:Verdana;;" id="$panel_id">};
	$out .= qq{<table class="table table-striped table-condensed table-bordered table-hover table-mybordered" style="vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;line-height: 25px;min-height: 25px;height: 25px;box-shadow: 3px 3px 5px #555;">};
	$out .= $print_html->print_header("background-color:aliceblue;color:black");
	my $color_validation = "grey";
	my $i = 0;
	foreach my $var_id (sort @$list_variants) {
		next if not exists $hVariantsDetails->{$var_id};
		next if not $hVariantsDetails->{$var_id}->{polyviewer_html_details_proj_pat};
		$i++;
		if ($i == 20) {
			$i = 0;
			print '.';
		}
		my $opacity = undef;
		my $var = $self->project->_newVariant($var_id);
		my $polyviewer_variant = $hVariantsDetails->{$var_id}->{polyviewer_variant};
		$polyviewer_variant->{gene} = $g;
		$polyviewer_variant->{transcripts} = $polyviewer_variant->{hgenes}->{$g->{id}}->{tr};
		my @list_polyviewer_h_details = @{$hVariantsDetails->{$var_id}->{polyviewer_html_details_proj_pat}};
		$print_html->variant($polyviewer_variant);
		$out .= $self->print_line_variant_all_patients($polyviewer_variant, $print_html, \@list_polyviewer_h_details, $opacity);
	}
	$out .= qq{</table>};
	$out .= qq{</div>};
	
	$gene = undef;
	return $out;
}


sub print_line_variant_all_patients {
	my ($self,$polyviewer_variant,$print_html,$list_h_details,$opacity) = @_;
	my $style = { style => "max-height:180px;" };
	$style = { style => "background-color: #DAEEED;opacity:0.5;max-height:180px;" }  if $opacity;
	my $out;
	my $hpatients;
	my $i = 0;
	my $list_print_html;
	foreach my $h_proj_pat (@{$list_h_details}) {
		my $b = new GBuffer;
		my $pr = $b->newProject(-name => $h_proj_pat->{project_name});
		$pr->getPatients();
		$pr->getFamilies();
		$pr->project_root_path();
		my $pat = $pr->getPatient($h_proj_pat->{patient_name});
		$pat->getFamily->name();
		foreach my $this_pat (@{$pat->getFamily->getPatients()}) { $this_pat->alignmentMethods(); }
		my $this_print_html = polyviewer_html->new( project=>$pr, patient=>$pat,header=>\@headers, bgcolor=>"background-color:#607D8B" );
		$this_print_html->variant($polyviewer_variant);
		#$this_print_html->variant->{patients_calling}->{$pat->id()}->{gt} = '';
		$this_print_html->variant->{patients_calling}->{$pat->id()}->{pc} = $h_proj_pat->{ratio};
		$this_print_html->variant->{patients_calling}->{$pat->id()}->{dp} = 'DP:'.$h_proj_pat->{dp};
		$this_print_html->variant->{patients_calling}->{$pat->id()}->{model} = $h_proj_pat->{model};
		push(@$list_print_html, $this_print_html);
	}
	my $cgi = $print_html->cgi;
	my $icon = qq{<img width="32" height="32" src="https://img.icons8.com/external-gliphyline-royyan-wijaya/32/external-laptop-laptop-collection-glyphyline-gliphyline-royyan-wijaya-15.png" alt="external-laptop-laptop-collection-glyphyline-gliphyline-royyan-wijaya-15"/>};
	$icon   = qq{<img width="24" height="24" src="https://img.icons8.com/external-tal-revivo-filled-tal-revivo/24/external-live-preview-of-a-smart-class-education-school-filled-tal-revivo.png" alt="external-live-preview-of-a-smart-class-education-school-filled-tal-revivo"/>};
	my $dropdown = qq{
		<div class="dropdown">
		<button class="btn btn-primary btn-xs dropdown-toggle " type="button" id="dropdownMenuButton" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false" style="font-size:10px;background-color:#C67FAE">
		$icon
		</button>
		<div class="dropdown-menu" aria-labelledby="dropdownMenuButton" style="font-size:12px;background-color:beige;color:black">
	};
	$dropdown .= "<li>".$print_html->mobidetails()."</li>";
  	$dropdown .= "<li>".$print_html->gnomadurl()."</li>";
	$dropdown .= "<li>".$print_html->alamuturl()."</li>";
	$dropdown .= "<li>".$print_html->varsome()."</li>";
	$dropdown .= qq{</div></div>};
	$out .= $cgi->td( $style, $dropdown );
	$out .= "\n";
	
	my $var_text = $print_html->var_name();
	$var_text .= "<br><i><b>HG19: ".$self->{hash_lift_variants}->{$polyviewer_variant->gnomad_id()}->{chr19}.'-'.$self->{hash_lift_variants}->{$polyviewer_variant->gnomad_id()}->{pos19}.'</b></i>';
	$out .= $cgi->td($style, $var_text);
	$out .= "\n";
	
	my (@l_html_calling, $h_fam_done);
	foreach my $this_print_html (@$list_print_html) {
		my $html_pat = $this_print_html->calling();
		my $project_name = $this_print_html->patient->getProject->name();
		next if exists $h_fam_done->{$project_name.'_'.$this_print_html->patient->getFamily->name()};
		$h_fam_done->{$project_name.'_'.$this_print_html->patient->getFamily->name()} = undef;
		my (@l_names, @l_bam);
		foreach my $pat (@{$this_print_html->patient->getFamily->getPatients()}) {
			push(@l_names, $pat->name());
			push(@l_bam, $pat->bamUrl());
		}
		my $pnames = join(';', @l_names);
		my $f = join(';', @l_bam);
		my $gn = $this_print_html->patient->getProject->getVersion();
		my $locus;
		if ($gn =~ /HG19/) { $locus = $self->{hash_lift_variants}->{$polyviewer_variant->gnomad_id()}->{chr19}.':'.$self->{hash_lift_variants}->{$polyviewer_variant->gnomad_id()}->{pos19}.'-'.$self->{hash_lift_variants}->{$polyviewer_variant->gnomad_id()}->{pos19}; }
		else { $locus = $polyviewer_variant->locus(); }
		my $igv_b = qq{<button class='igvIcon2' onclick='launch_web_igv_js("$project_name","$pnames","$f","$locus","/","$gn")' style="color:black"></button>};
		my $project_description = $this_print_html->patient->getProject->name().', Description: '.$this_print_html->patient->getProject->description();
		push(@l_html_calling, "<tr style='padding:6px;'><td style='padding-right:10px;'><center><button onClick='alert(\"$project_description\")'>".$this_print_html->patient->getProject->name().'</button><br><b>'.$this_print_html->patient->getFamily->name()."</b><td style='padding-right:5px;'>".$igv_b."</td></center></td><td>".$html_pat."</td></tr>");
		push(@l_html_calling, "<tr style='padding:6px;'><td><br></td><td><br></td></tr>");
	}
	
	my $out_calling = qq{<center><table style='width:98%;max-height:200px;'>};
	$out_calling .= join('', @l_html_calling);
	$out_calling .= qq{</table></center>};
	
	$out .= $cgi->td( $style, $out_calling) ;
	$out .= "\n";
	$out .= $cgi->td( $style, $print_html->gnomad() );
	$out .= "\n";
	$out .= $cgi->td( $style, $print_html->dejavu() );
	$out .= "\n";
	$out .= $cgi->td( $style, $print_html->validations );
#	$out .= $cgi->td( $style, '' );
	$out .= "\n";
	$out .= $cgi->td( $style, $print_html->transcripts() );
	$out .= $cgi->end_Tr();
	$out .= "\n";
	
	$list_print_html = undef;
	return $out;
}



1;