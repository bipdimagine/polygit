#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$| = 1;
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use Data::Dumper;
use JSON;

use lib "$Bin/../GenBo/lib/obj-nodb/";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/html_polygenescout.pm";
use lib "$Bin/../GenBo";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/validation_variation";
use lib "$Bin/../cache_nodb/scripts/";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use GBuffer;
use export_data;
use JSON;
use polyweb_dude;
use GBuffer;
use GenBoProject;
use session_export;

my $cgi                = new CGI();
my $user               = $cgi->param('user');
my $annotation         = $cgi->param('annotation');
my $transcripts		   = $cgi->param('transcripts');
my $genes              = $cgi->param('genes');
my $locus              = $cgi->param('locus');
my $capture            = $cgi->param('capture');
my $phenotype          = $cgi->param('phenotype');
my $hgmd_concept       = $cgi->param('hgmd_concept');
my $use_project        = $cgi->param('use_project');
my $sort			   = $cgi->param('sort');
my $trio_project	   = $cgi->param('project');
my $trio_patient	   = $cgi->param('patient');
my $use_release	   		= $cgi->param('release');
my $use_phenotype_score	   = $cgi->param('use_phenotype');
$use_phenotype_score = lc($phenotype) if ($phenotype and not $use_phenotype_score);
$use_phenotype_score = 'intellectual disability' unless ($use_phenotype_score);

#$use_release = 'HG38';

my $buffer  = GBuffer->new();
purge_cgi_session_directory($buffer);
my $project;
if ($trio_project) {
	$project = $buffer->newProjectCache( -name => $trio_project );
}
else {
	
	#TODO: prendre projet HG38
	my $project_name = 'NGS2016_1231';
	if ($use_release eq 'HG38') { $project_name = 'NGS2024_7792'; }
	$project = $buffer->newProject( -name => $project_name );
	
	#$project = $buffer->newProject( -name => $buffer->get_random_project_name_with_this_annotations_and_genecode() );
	my ( $genecode_version, $annot_version );
	unless ($annotation) {
		$genecode_version = $buffer->getQuery()->getMaxGencodeVersion();
		$annot_version    = $buffer->getQuery()->getMaxPublicDatabaseVersion();
		$annotation       = $genecode_version . '.' . $annot_version;
	}
	$project->changeAnnotationVersion( $annotation, 1 );
}

if ($use_phenotype_score) {
	my $h_pheno_infos = $buffer->queryPhenotype->getPhenotypeInfosFromName($use_phenotype_score);
	confess("\n\nERROR: $use_phenotype_score not found in our DB. die.\n\n") unless (exists $h_pheno_infos->{$use_phenotype_score});
	my @l_pheno;
	push (@l_pheno, $use_phenotype_score);
	my $phenotype_id = $h_pheno_infos->{$use_phenotype_score}->{phenotype_id};
	$project->{phenotypes} = \@l_pheno;
	$project->{phenotypes_object}->{$h_pheno_infos->{$use_phenotype_score}->{phenotype_id}} = undef;
}

my $hRes;
my ( @lTranscripts );
my @lSearched;
my ($h_genes, $h_genes_tmp);

if ($phenotype) {
	push(@lSearched, 'Phenotype ' . $phenotype);
	my $buffer_tmp = GBuffer->new();
	foreach my $gene_name (@{ $buffer_tmp->queryPhenotype->getGenesFormPhenotypeName($phenotype) }) {
		$h_genes->{$gene_name} = $gene_name;
	}
	if ( $buffer_tmp->hasHgmdAccess($user) and $hgmd_concept ) {
		push(@lSearched, 'Phenotype ' . $phenotype . ', HGMD Concept ' . $hgmd_concept);
		foreach my $gene_name (@{ $buffer_tmp->queryHgmd->get_genes_from_concept_name($hgmd_concept)}){
			$h_genes->{$gene_name} = $gene_name;
		}
	}
	$buffer_tmp = undef;
}

if ($capture) {
	push(@lSearched, 'Capture ' . $capture);
	my ( $hCapture, $hTr );
	foreach my $this_capture ( split( ';', $capture ) ) {
		$hCapture->{ lc($this_capture) } = undef;
	}
	my $buffer_tmp  = GBuffer->new();
	my $project_tmp = $buffer_tmp->newProject( -name => $use_project );
	foreach my $c ( @{ $project_tmp->getCaptures() } ) {
		next if ( not exists $hCapture->{ lc( $c->name() ) } );
		foreach my $tr_name ( @{ $c->transcripts_name() } ) {
			$hTr->{$tr_name} = undef;
		}
	}
	$transcripts = join( ',', keys %$hTr );
	$project_tmp        = undef;
	$buffer_tmp         = undef;
}

if ($genes) {
	push(@lSearched,  'Gene(s)');
	foreach my $g (split( ',', $genes )) {
		next if (exists $h_genes->{$g});
		$h_genes->{$g} = $g;
	}
}

my ($h_found_by_locus, $h_found_by_transcript);
if ($locus) {
	push(@lSearched, 'Locus ' . $locus);
	foreach my $this_locus ( split( ',', $locus ) ) {
		my $init_this_locus = $this_locus;
		$this_locus =~ s/:/_/g;
		$this_locus =~ s/-/_/g;
		my ( $chr_name, $start, $end ) = split( '_', $this_locus );
		$end = $start + 1 unless ($end);
		if ($start>$end){
			my $t = $start;
			$start = $end;
			$end = $t;
		}
		my $chr = $project->getChromosome($chr_name);
		foreach my $gene ( @{ $chr->getGenesByPosition( $start, $end ) } ) {
			$h_found_by_locus->{$gene->external_name()}->{$init_this_locus} = undef;
			next if (exists $h_genes->{$gene->external_name()});
			$h_genes->{$gene->external_name()} = $gene->id();
			
		}
	}
}

if ($transcripts) {
	foreach my $tr_name ( split( ',', $transcripts ) ) {
		next unless ( $tr_name =~ /ENST/ );
		my ( $transcript, $gene );
		eval {
			$transcript = $project->newTranscript($tr_name);
			my $init_this_locus = $transcript->getChromosome->id().'_'.$transcript->start().'_'.$transcript->end();
			foreach my $gene ( @{ $transcript->getGenes() } ) {
				$h_found_by_transcript->{$gene->external_name()}->{$tr_name} = $init_this_locus;
				next if (exists $h_genes->{$gene->external_name()});
				$h_genes->{$gene->external_name()} = $gene->id();
			}
		};
		if ($@) { next; }
	}
}

my ($buffer_trio, $project_trio, $patient_trio);
my @headers;
if ($use_phenotype_score) {
	my $out_phenos;
	$out_phenos .= qq{<div class="input-group" style="width:100%">};
	$out_phenos .= qq{<select class="form-control" id="form_polyphenotypes_init" style="font-size:9px;height:auto;width:100%;">};
	my $list_all_phenotypes = $buffer->queryPhenotype->getAllPhenotypesName();
	foreach my $pheno (sort @$list_all_phenotypes) {
		my $value = lc($pheno);
		my $name = ucfirst($pheno);
		$name = ucfirst($name);
		my $cmd_score;
		my $selected;
		if ($value eq lc($use_phenotype_score)) {
			$selected = 'selected'
		}
		else {
			my @lArgs;
			push(@lArgs, 'user='.$user) if ($user);
			push(@lArgs, 'annotation='.$annotation) if ($annotation);
			push(@lArgs, 'transcripts='.$transcripts) if ($transcripts);
			push(@lArgs, 'genes='.$genes) if ($genes);
			push(@lArgs, 'locus='.$locus) if ($locus);
			push(@lArgs, 'capture='.$capture) if ($capture);
			push(@lArgs, 'phenotype='.$phenotype) if ($phenotype);
			push(@lArgs, 'hgmd_concept='.$hgmd_concept) if ($hgmd_concept);
			push(@lArgs, 'use_project='.$use_project) if ($use_project);
			push(@lArgs, 'sort='.$sort) if ($sort);
			push(@lArgs, 'use_phenotype='.$value);
			my $args = join('&', @lArgs);
			$cmd_score = qq{onClick="reload_score_gene('$args');"};
			#$cmd_score = qq{onClick="alert('$args');"};
		}
		$out_phenos .= qq{<option value='$value' $cmd_score $selected><span>$name</span></option>};
	}
	$out_phenos .= qq{</select>};
	$out_phenos .= qq{</div>};
	@headers = ( 'Gene Score<br>'.$out_phenos, qq{<button style="border-radius:10px" onClick="reload_and_sort('locus');"><span style="font-size:9px;" class="glyphicon glyphicon-sort-by-attributes"></span></button>  Locus}, qq{<button style="border-radius:10px" onClick="reload_and_sort('gene');"><span style="font-size:9px;" class="glyphicon glyphicon-sort-by-attributes"></span></button>  Gene Name}, 'Phenotype(s)', 'Omim', 'Gtex', 'pLi', 'Panel(s)' );
}
else {
	@headers = ( qq{<button style="border-radius:10px" onClick="reload_and_sort('locus');"><span style="font-size:9px;" class="glyphicon glyphicon-sort-by-attributes"></span></button>  Locus}, qq{<button style="border-radius:10px" onClick="reload_and_sort('gene');"><span style="font-size:9px;" class="glyphicon glyphicon-sort-by-attributes"></span></button>  Gene Name}, 'Phenotype(s)', 'Omim', 'Gtex', 'pLi', 'Panel(s)' );
}
if ( $buffer->hasHgmdAccess($user) ) {
	push(@headers, 'Concept(s)');
}
if ($trio_project and $trio_patient) {
	push(@headers, 'Gene Score');
}
push(@headers, 'DejaVu');
if ($trio_project and $trio_patient) {
	push(@headers, 'Variants '.$trio_patient);
}
if ( $buffer->hasHgmdAccess($user) ) {
	push(@headers, 'Total Mut');
	push(@headers, 'Total New Mut');
}

my $last_genecode = $buffer->getQuery->getMaxGencodeVersion();
my $proj_name = $buffer->get_random_project_name_with_this_annotations_and_genecode();
my $out;
$out .= "<div style='text-align:center;font-size:12px;width:100%;overflow-x:scroll;'>";
$out .= "<table id='table_hgmd_genes' class='table table-striped table-bordered' style='text-align:center;font-size:12px;width:100%;overflow-x:scroll;'>";
$out .= "<thead>";
$out .= $cgi->start_Tr( { style => "font-size:11px;" } );

my $hgmd_version_used = $buffer->queryHgmd->database();
my $searched = join(' - ', @lSearched);


my $gencode_version = $project->gencode_version();
my ($trio_phenotype, $used_hgmd);
if ($trio_project and $trio_patient) {
	$trio_phenotype = get_project_phenotype($trio_project);
	my $h_pheno_infos = $buffer->queryPhenotype->getPhenotypeInfosFromName($trio_phenotype);
	$used_hgmd = $h_pheno_infos->{$use_phenotype_score}->{concept};
	$out .= "<center><b>Gencode:</b> <font color='green'>$gencode_version</font> | <b>Phenotype:</b> <font color='green'>$trio_phenotype</font> | <b>HGMD Concept:</b> <font color='green'>$used_hgmd</font> | <b>Searched:</b> <font color='green'>$searched</font><b> | HGMD:</b> <font color='green'>$hgmd_version_used</font></center>";
	$out .= $cgi->start_Tr({style=>"font-size:12px;"});
	$out .= $cgi->th({colspan=>8, style=>"text-align:center;background-color:white;font-size:11px;color:#607D8B;"}, "<b>GENE Informations</b>");
	$out .= $cgi->th({colspan=>1, style=>"text-align:center;background-color:white;font-size:11px;color:#607D8B;"}, "<b>HGMD</b>") if ( $buffer->hasHgmdAccess($user) );
	#$out .= $cgi->th({colspan=>2, style=>"text-align:center;background-color:white;font-size:11px;color:#607D8B;"}, "<b>POLYWEB</b>");
	$out .= $cgi->th({colspan=>1, style=>"text-align:center;background-color:white;font-size:11px;color:#607D8B;"}, "<b>POLYWEB</b>");
	$out .= $cgi->th({colspan=>2, style=>"text-align:center;background-color:white;font-size:11px;color:#607D8B;"}, "<b>HGMD</b>") if ( $buffer->hasHgmdAccess($user) );
	$out .= $cgi->end_Tr();
	$out .= $cgi->start_Tr({style=>"font-size:11px;"});
	#@headers = ( 'Gene Score', qq{<button style="border-radius:10px" onClick="reload_and_sort('locus');"><span style="font-size:9px;" class="glyphicon glyphicon-sort-by-attributes"></span></button>  Locus},  qq{<button style="border-radius:10px" onClick="reload_and_sort('gene');"><span style="font-size:9px;" class="glyphicon glyphicon-sort-by-attributes"></span></button>  Gene Name}, 'Phenotype(s)', 'Omim', 'Gtex', 'pLi', 'Panel(s)','Concept(s)', 'DejaVu', 'Variants', 'Total Mut', 'Total New Mut' );
	@headers = ( 'Gene Score', qq{<button style="border-radius:10px" onClick="reload_and_sort('locus');"><span style="font-size:9px;" class="glyphicon glyphicon-sort-by-attributes"></span></button>  Locus},  qq{<button style="border-radius:10px" onClick="reload_and_sort('gene');"><span style="font-size:9px;" class="glyphicon glyphicon-sort-by-attributes"></span></button>  Gene Name}, 'Phenotype(s)', 'Omim', 'Gtex', 'pLi', 'Panel(s)','Concept(s)', 'Variants', 'Total Mut', 'Total New Mut' );
}
else {
	if ( $buffer->hasHgmdAccess($user) ) {
		my $h_pheno_infos = $buffer->queryPhenotype->getPhenotypeInfosFromName($use_phenotype_score);
		$used_hgmd = $h_pheno_infos->{$use_phenotype_score}->{concept};
		$out .= "<center><b>Gencode:</b> <font color='green'>$gencode_version</font> | <b>Phenotype:</b> <font color='green'>$use_phenotype_score</font> | <b>HGMD Concept:</b> <font color='green'>$used_hgmd</font> | <b>Searched:</b> <font color='green'>$searched</font><b> | HGMD:</b> <font color='green'>$hgmd_version_used</font></center>";
		#$out .= "<center>".$out_phenos."</center>";
		$out .= $cgi->start_Tr({style=>"font-size:11px;"});
		if ($trio_project and $trio_patient) {
			my $nb_col = 8;
			$nb_col++ if ($use_phenotype_score);
			$out .= $cgi->th({colspan=>$nb_col, style=>"text-align:left;background-color:white;font-size:11px;color:#607D8B;"}, "<b>GENE Informations</b>");
		}
		else {
			my $nb_col = 7;
			$nb_col++ if ($use_phenotype_score);
			$out .= $cgi->th({colspan=>$nb_col, style=>"text-align:left;background-color:white;font-size:11px;color:#607D8B;"}, "<b>GENE Informations</b>");
		}
		$out .= $cgi->th({colspan=>1, style=>"text-align:center;background-color:white;font-size:11px;color:#607D8B;"}, "<b>HGMD</b>");
		if ($trio_project and $trio_patient) {
			$out .= $cgi->th({colspan=>1, style=>"text-align:center;background-color:white;font-size:11px;color:#607D8B;"}, "<center><b>$trio_project</b></center>");
		}
		$out .= $cgi->th({colspan=>1, style=>"text-align:center;background-color:white;font-size:11px;color:#607D8B;"}, "<b>PolyWeb</b>");
		$out .= $cgi->th({colspan=>2, style=>"text-align:center;background-color:white;font-size:11px;color:#607D8B;"}, "<center><b>HGMD Informations</b></center>");
		$out .= $cgi->end_Tr();
		$out .= $cgi->start_Tr({style=>"font-size:11px;"});
	}
	else {
		$out .= "<center><b>Phenotype:</b> <font color='green'>$use_phenotype_score</font> | <b>Searched:</b> <font color='green'>$searched</font></center>";
		#$out .= "<center>".$out_phenos."</center>";
		my $out_th0;
		$out_th0 .= qq{<div class="row" style="padding:4px;">};
		$out_th0 .= qq{<div class="col-sm-6">};
		$out_th0 .= qq{<span style="font-size:12px;"><br><b>GENE Informations</b></span>};
		$out_th0 .= qq{</div>};
		$out_th0 .= qq{</div>};
		$out .= $out_th0;
	}
}

foreach my $h (@headers) {
	$out .= $cgi->th({ style => "background-color:white;font-size:12px;color:#607D8B;" }, "<center><b>" . ucfirst($h) . "</b></center>");
}

$out .= $cgi->end_Tr();
$out .= "</thead>";
$out .= "<tbody>";

if ($trio_project and $trio_patient) {
	$sort = 'gene_score' unless $sort;
}

my $h_lines_sorted;
my $nb_genes_error = 0;
my $last_gene_name;
if ($h_genes) {
	foreach my $gene_key_name ( sort keys %$h_genes ) {
		my $gene;
		my $gene_id = $h_genes->{$gene_key_name};
		eval { $gene = $project->newGene($gene_id); };
		if ( $@ or not defined($gene) ) {
			$hRes->{$gene_key_name} = 'NOT;' . $annotation;
			
			my $out_line = $cgi->start_Tr({ style => "background-color:#f5a260;font-size:11px;max-height:60px;overflow-y: auto;" });
			$out_line .= $cgi->td({ rowspan => 1, colspan => 1, style =>"max-height:60px;overflow-y:auto;padding-left:20px;font-size:12px;"}, qq{<span style="color:red;font-size:12px;">Not Found Genecode $last_genecode</span>});
			if ($gene_key_name ne $gene_id) {
				$out_line .= $cgi->td({ rowspan => 1, colspan => 1, style =>"max-height:60px;overflow-y:auto;padding-left:20px;font-size:12px;"}, qq{<span style="color:white;"><strike>$gene_key_name / $gene_id</strike></span>});
			}
			else {
				$out_line .= $cgi->td({ rowspan => 1, colspan => 1, style =>"max-height:60px;overflow-y:auto;padding-left:20px;font-size:12px;"}, qq{<span style="color:white;"><strike>$gene_key_name</strike></span>});
			}
			if ($trio_project and $trio_patient) {
				$out_line .= $cgi->td({ rowspan => 1, colspan => 1, style =>"max-height:60px;overflow-y:auto;padding-left:20px;font-size:12px;"}, qq{<span style="color:white;"><strike>-</strike></span>});
			}
			if ( $buffer->hasHgmdAccess($user) ) {
				$out_line .= $cgi->td({ rowspan => 1, colspan => 9, style =>"max-height:60px;overflow-y:auto;padding-left:20px;font-size:12px;"}, '');
			}
			else {
				$out_line .= $cgi->td({ rowspan => 1, colspan => 6, style =>"max-height:60px;overflow-y:auto;padding-left:20px;font-size:12px;"}, '');
			}
			$out_line .= $cgi->end_Tr();
			
			if ($sort and $sort eq 'locus') {
				$h_lines_sorted->{999999}->{$nb_genes_error}->{$nb_genes_error} = $out_line;
			}
			elsif ($sort and $sort eq 'gene_score') {
				$h_lines_sorted->{-999}->{$nb_genes_error} = $out_line;
			}
			else {
				$h_lines_sorted->{$gene_key_name} = $out_line;
			}
			$nb_genes_error++;
			next;
		}

		my ( $hResGene, @lAllVar );
		my $external_name = $gene->external_name();
		
		my ( $nb_concepts, $concept_list, $total_mut, $total_new );
		if ( $gene->hgmd() ) {
			$total_mut = $buffer->queryHgmd->get_gene_mut_total($external_name);
			$total_new = $buffer->queryHgmd->get_gene_new_mut_total($external_name);
			my $h_concepts = $buffer->queryHgmd->get_gene_concept($external_name);
			$concept_list = "<center><b><u>Concept(s) for gene $external_name</b></u><br><br>";
			foreach my $concept_name ( sort keys %{$h_concepts} ) {
				next if ( $h_concepts->{$concept_name}->{num_matching} eq '0' );
				$concept_list .= $concept_name . "<br>";
				$nb_concepts++;
			}
			$concept_list .= "</center>";
		}
		my $gene_name = $external_name;
		$last_gene_name                             = $external_name;
		$hResGene->{$gene_name}->{id}               = $gene->id;
		$hResGene->{$gene_name}->{external_name}    = $external_name;
		$hResGene->{$gene_name}->{pLI}              = $gene->pLI();
		$hResGene->{$gene_name}->{omim_id}          = $gene->omim_id();
		$hResGene->{$gene_name}->{omim_inheritance} = $gene->omim_inheritance();
		$hResGene->{$gene_name}->{variants}         = \@lAllVar;
		my ( $pheno, $nb_other_terms ) = $gene->polyviewer_phentotypes();
		$hResGene->{$gene_name}->{phenotypes}->{pheno} = $pheno;
		$hResGene->{$gene_name}->{phenotypes}->{nb_other_terms} = $nb_other_terms;
		my $cmd_link = qq{launch_for_gene('$gene_id');};
		$hResGene->{$gene_name}->{specific_cmd} = '';
		$hResGene->{$gene_name}->{max_score}    = $gene->score() * 1.0;
		my $description_gene = $gene->description();
		eval {
			foreach my $panel ( @{ $gene->getPanels() } ) {
				$hResGene->{$gene_name}->{panels}->{ $panel->name() }->{phenotype} = $panel->getPhenotypes()->[0]->name();
			}
		};
		if ($@) { $hResGene->{$gene_name}->{panels} = undef; }
		

		my $out_line = $cgi->start_Tr( {style =>"font-size:11px;max-height:60px;overflow-y: auto;vertical-align:center !important;"} );

		my $gene_score;
		if ($trio_project and $trio_patient) {
			$sort = 'gene_score' unless $sort;
			$gene_score = get_gene_score($gene);
			$out_line .= html_polygenescout::print_gene_score($gene_score);
			$out_line .= html_polygenescout::print_locus($gene);
			$out_line .= html_polygenescout::print_gene_basic_tables_infos($hResGene->{$gene_name}, $trio_phenotype);
			if ( $buffer->hasHgmdAccess($user) ) {
				$out_line .= html_polygenescout::print_hgmd_concepts($nb_concepts, $concept_list, $used_hgmd);
			}
#			$out_line .= html_polygenescout::print_link_dejavu($cmd_link);
			$out_line .= html_polygenescout::print_tab_variants_project($gene->getChromosome->id(), $external_name, $trio_project, $trio_patient, $h_genes,$gene->project);
			if ( $buffer->hasHgmdAccess($user) ) {
				$out_line .= html_polygenescout::print_hgmd_total_mut($total_mut, $external_name);
				$out_line .= html_polygenescout::print_hgmd_total_new($total_new, $external_name);
			}
			
		}
		else {
			if ($use_phenotype_score) {
				$sort = 'gene_score' unless $sort;
				$gene_score = get_gene_score($gene);
				#$gene_score = get_gene_score_from_phenotype($external_name, $use_phenotype_score);
				$out_line .= html_polygenescout::print_gene_score($gene_score);
			}
			$out_line .= html_polygenescout::print_locus($gene);
			if ($use_phenotype_score) {
				$out_line .= html_polygenescout::print_gene_basic_tables_infos($hResGene->{$gene_name}, $use_phenotype_score);
			}
			else {
				$out_line .= html_polygenescout::print_gene_basic_tables_infos($hResGene->{$gene_name});
			}
			if ( $buffer->hasHgmdAccess($user) ) {
				
				if ($use_phenotype_score) {
					$out_line .= html_polygenescout::print_hgmd_concepts($nb_concepts, $concept_list, $used_hgmd);
				}
				else {
					$out_line .= html_polygenescout::print_hgmd_concepts($nb_concepts, $concept_list);
				}
			}
			my $hash_locus_intervals;
			my $gene_chr = $gene->getChromosome->id();
			my $gene_start = $gene->start();
			my $gene_end = $gene->end();
			if ($gene_start > $gene_end){
				my $t = $gene_start;
				$gene_start = $gene_end;
				$gene_end = $t;
			}
			if ($transcripts) {
				my @list_transcripts  = keys %{$h_found_by_transcript->{$external_name}};
				foreach my $tr_name (sort @list_transcripts) {
					my $this_locus = $h_found_by_transcript->{$external_name}->{$tr_name};
#					my ( $locus_chr, $locus_start, $locus_end ) = split( '_', $this_locus );
#					$locus_end = $locus_start + 1 unless ($locus_end);
#					next if ($gene_chr ne $locus_chr);
#					next if ($gene_start > $locus_end);
#					next if ($gene_end < $locus_start);
#					next if ($gene_start > $locus_start and $gene_start < $locus_end and $gene_end > $locus_start and $gene_end < $locus_end);
#					my $intspan_locus = Set::IntSpan::Fast::XS->new();
#					$intspan_locus->add_range($locus_start, $locus_end);
#					my $intspan_gene = Set::IntSpan::Fast::XS->new();
#					$intspan_gene->add_range($gene_start, $gene_end);
#					my $in = $intspan_locus->intersection( $intspan_gene );
#					my $new_locus = $gene_chr.':'.$in->as_string();
#					$new_locus =~ s/_/-/;
					my $cmd_link_locus = qq{launch_for_gene('$gene_id','$this_locus', '$tr_name');};
					$hash_locus_intervals->{$tr_name} = $cmd_link_locus;
				}
			}
			if ($locus) {
				my @list_locus = keys %{$h_found_by_locus->{$external_name}};
				foreach my $this_locus (sort @list_locus) {
					$this_locus =~ s/chr//g;
					$this_locus =~ s/:/_/g;
					$this_locus =~ s/-/_/g;
					my ( $locus_chr, $locus_start, $locus_end ) = split( '_', $this_locus );
					$locus_end = $locus_start + 1 unless ($locus_end);
					if ($locus_start > $locus_end){
						my $t = $locus_start;
						$locus_start = $locus_end;
						$locus_end = $t;
					}
					next if ($gene_chr ne $locus_chr);
					next if ($gene_start > $locus_end);
					next if ($gene_end < $locus_start);
					next if ($gene_start > $locus_start and $gene_start < $locus_end and $gene_end > $locus_start and $gene_end < $locus_end);
					my $intspan_locus = Set::IntSpan::Fast::XS->new();
					$intspan_locus->add_range($locus_start, $locus_end);
					my $intspan_gene = Set::IntSpan::Fast::XS->new();
					$intspan_gene->add_range($gene_start, $gene_end);
					my $in = $intspan_locus->intersection( $intspan_gene );
					my $new_locus = $gene_chr.':'.$in->as_string();
					$new_locus =~ s/_/-/;
					my $cmd_link_locus = qq{launch_for_gene('$gene_id','$new_locus');};
					$hash_locus_intervals->{$new_locus} = $cmd_link_locus;
				}
			}
			if ($hash_locus_intervals) {
				$out_line .= html_polygenescout::print_link_dejavu_with_locus($cmd_link, $hash_locus_intervals);
			}
			else {
				$out_line .= html_polygenescout::print_link_dejavu($cmd_link);
			}
			if ( $buffer->hasHgmdAccess($user) ) {
				$out_line .= html_polygenescout::print_hgmd_total_mut($total_mut, $external_name);
				$out_line .= html_polygenescout::print_hgmd_total_new($total_new, $external_name);
			}
		}
		$out_line .= $cgi->end_Tr();
		$gene_score = 0 unless ($gene_score);
		if ($sort and $sort eq 'locus') {
			my $chr_id = $gene->getChromosome->id();
			$chr_id = '23' if ($chr_id eq 'X');
			$chr_id = '24' if ($chr_id eq 'Y');
			$chr_id = '25' if ($chr_id eq 'M');
			$chr_id = '25' if ($chr_id eq 'MT');
			$h_lines_sorted->{$chr_id}->{$gene->start()}->{$gene->end()} = $out_line;
		}
		elsif ($sort and $sort eq 'gene_score') {
			$h_lines_sorted->{$gene_score}->{$external_name} = $out_line;
		}
		else {
			$h_lines_sorted->{$external_name} = $out_line;
		}
	}
}

if ($sort and $sort eq 'locus') {
	foreach my $chr_id (sort {$a <=> $b} keys %$h_lines_sorted) {
		foreach my $start (sort {$a <=> $b} keys %{$h_lines_sorted->{$chr_id}}) {
			foreach my $end (sort {$a <=> $b} keys %{$h_lines_sorted->{$chr_id}->{$start}}) {
				$out .= $h_lines_sorted->{$chr_id}->{$start}->{$end};
			}
		}
	}
}
elsif ($sort and $sort eq 'gene_score') {
	foreach my $gene_score (sort {$b <=> $a} keys %$h_lines_sorted) {
		foreach my $external_name (sort keys %{$h_lines_sorted->{$gene_score}}) {
			$out .= $h_lines_sorted->{$gene_score}->{$external_name};
		}
	}
}
else {
	foreach my $external_name (sort keys %$h_lines_sorted) {
		$out .= $h_lines_sorted->{$external_name}; 
	}
}
$out .= "</tbody></table></div>";
$hRes = undef;
$hRes->{html} = $out;
$hRes->{gencode} = $gencode_version;
$hRes->{nb_genes} = scalar(keys %$h_genes);
$hRes->{nb_genes_error} = $nb_genes_error if ( $nb_genes_error > 0 );

print $cgi->header('text/json-comment-filtered');
my $json_encode = encode_json $hRes;
print $json_encode;
exit(0);


sub purge_cgi_session_directory {
	my $buffer = shift;
	my $dir_sessions = $buffer->config->{project_data}->{global_search};
	return unless (-d $dir_sessions);
	opendir my $dir, $dir_sessions or die "Cannot open directory: $!";
	my @files = readdir $dir;
	closedir $dir;
	foreach my $f (@files) {
		next unless ($f =~ /cgisess_/);
		my ($base, $tmp_session_id) = split('_', $f);
		#my $session = new session_export();
		#$session->tmp_dir($buffer->config->{project_data}->{global_search});
		#warn Dumper $session->load_session( $tmp_session_id );
		#$session->check_if_expired();
	}
}

sub get_project_phenotype {
	my $buffer_trio = GBuffer->new();
	my $project_trio = $buffer->newProject( -name => $trio_project );
	return $project_trio->phenotypes->[0];
}

sub get_gene_score {
	my ($gene) = @_;
	return $gene->raw_score();
}

sub get_gene_score_from_phenotype {
	my ($external_name, $use_phenotype_score) = @_;
	my $gene_trio = $project->newGene($h_genes->{$external_name});
	return $gene_trio->raw_score();
}
