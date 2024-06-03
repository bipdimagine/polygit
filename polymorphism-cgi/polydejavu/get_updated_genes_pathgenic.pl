#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;

use strict;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/validation_variation"; 
use lib "$Bin/../cache_nodb/scripts/";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use GBuffer;
use export_data;
use JSON;
use polyweb_dude;
use VcfMerge;
use GenBoNoSql;
use Set::IntervalTree;
use Spreadsheet::WriteExcel;
use Bit::Vector;
use Bit::Vector::Overload;
use Compress::Snappy;
use Storable qw(store retrieve freeze dclone thaw);
use POSIX qw(strftime);
use List::MoreUtils qw(natatime);
use CGI::Session;
use html; 
use Carp;
use Cache_Commons;
use QueryVectorFilter;
use IO::Handle;
use Time::Local;
use POSIX qw(strftime);
use session_export;

require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/html_polygenescout.pm";

my $io = IO::Handle->new();
$io->autoflush(1);


my $cgi = new CGI();
my $user_name = $cgi->param('user');
my $pwd = $cgi->param('pwd');
my $filter_only_new = $cgi->param('filter_only_new');
my $filter_concept = $cgi->param('filter_concept');
my $filter_inheritance = $cgi->param('filter_inheritance');
my $filter_search = $cgi->param('filter_search');
my $sort = $cgi->param('sort');
my $get_list_concepts = $cgi->param('get_list_concepts');
my $get_all_genes = $cgi->param('get_all_genes');
my $previous_hgmd_version = $cgi->param('previous_hgmd_db');

$sort = 'gene_score' unless ($sort);

if ($filter_concept) { $filter_concept =~ s/\+/ /g; }
$filter_concept = undef if ($filter_concept eq 'all');

$user_name = lc($user_name);
my $buffer = new GBuffer;

$buffer->queryHgmd->{previous_database} = $previous_hgmd_version if ($previous_hgmd_version);

my $nb_genes_all = 0;
my $nb_genes_selected = 0;

my ($h_genes_all, $h_genes_with_dm, $h_from_d, $h_from_go_terms, $h_from_concepts, $h_from_exp_inheritance);

if ($get_list_concepts) {
	my $h_concepts = $buffer->queryHgmd->get_hash_concepts();
	my $out_phenos;
	$out_phenos .= qq{<div class="input-group" style="width:100%">};
	$out_phenos .= qq{<select class="form-control" id="form_concepts_init" style="font-size:9px;height:auto;width:100%;">};
	$out_phenos .= qq{<option value=''><span></span></option>};
	foreach my $concept (sort keys %{$h_concepts}) {
		my $value = lc($concept);
		my $name = ucfirst($concept);
		$name = ucfirst($name);
		$out_phenos .= qq{<option value='$value'><span>$name</span></option>};
	}
	$out_phenos .= qq{</select>};
	$out_phenos .= qq{</div>};
	my $hRes;
	$hRes->{html} = $out_phenos;
	my $json_encode = encode_json $hRes;
	print $cgi->header('text/json-comment-filtered');
	print $json_encode;
	exit(0);
}

my (@l_genes_name, $h_genes_names_with_DM);
my $h = $buffer->queryHgmd->get_hash_last_released_DM();
my $h_genes_new = $buffer->queryHgmd->get_hash_new_genes_with_concept();
foreach my $hgmd_id (keys %$h) {
	my $gene_name = $h->{$hgmd_id}->{gene};
	$h_genes_names_with_DM->{$gene_name}->{$hgmd_id} = undef;
}
if ($filter_concept and $get_all_genes) {
	@l_genes_name = keys %$h_genes_new;
}
else {
	foreach my $gene_name (keys %$h_genes_names_with_DM) {
		delete $h_genes_names_with_DM->{$gene_name} unless (exists $h_genes_new->{$gene_name});
	}
	@l_genes_name = keys %$h_genes_names_with_DM;
}

foreach my $gene_name (@l_genes_name) {
	$h_genes_all->{$gene_name} = undef;
	
	my (@lConcepts);
	my $has_new_concept;
	
	if ($filter_concept and $filter_only_new) {
		next unless exists $h_genes_new->{$gene_name}->{lc($filter_concept)};
		next if ($h_genes_new->{$gene_name}->{lc($filter_concept)}->{num_matching} eq '0');
		next unless (exists $h_genes_new->{$gene_name}->{lc($filter_concept)}->{new});
		next unless ($h_genes_new->{$gene_name}->{lc($filter_concept)}->{new} == 1);
	}
	elsif ($filter_concept) {
		next if ($h_genes_new->{$gene_name}->{lc($filter_concept)}->{num_matching} eq '0');
		next unless exists $h_genes_new->{$gene_name}->{lc($filter_concept)};
	}
	elsif ($filter_only_new) {
		my $concept_ok = undef;
		foreach my $concept_name (keys %{$h_genes_new->{$gene_name}}) {
			next if ($h_genes_new->{$gene_name}->{$concept_name}->{num_matching} eq '0');
			$concept_ok = 1 if ($h_genes_new->{$gene_name}->{$concept_name}->{new} == 1);
		}
		next unless $concept_ok;
	}
	
	foreach my $concept_name (keys %{$h_genes_new->{$gene_name}}) {
		my $num_matching = $h_genes_new->{$gene_name}->{$concept_name}->{num_matching};
		next if ($num_matching eq '0');
		my $percent = $h_genes_new->{$gene_name}->{$concept_name}->{percentage};
		my $num = $h_genes_new->{$gene_name}->{$concept_name}->{num_matching};
		my $total = $h_genes_new->{$gene_name}->{$concept_name}->{total};
		$h_genes_with_dm->{$gene_name}->{concepts}->{$concept_name}->{percentage} = $percent;
		$h_genes_with_dm->{$gene_name}->{concepts}->{$concept_name}->{ratio} = $num.'/'.$total;
		$h_genes_with_dm->{$gene_name}->{concepts}->{$concept_name}->{new} = $h_genes_new->{$gene_name}->{$concept_name}->{new};
		$has_new_concept = 1 if ($h_genes_new->{$gene_name}->{$concept_name}->{new} == 1);
		push(@lConcepts, $concept_name);
	}
	
	foreach my $concept_name (@lConcepts) {
		$h_from_concepts->{$concept_name}->{$gene_name} = undef;
	}
	
	if (exists $h_genes_names_with_DM->{$gene_name}) {
		foreach my $hgmd_id (keys %{$h_genes_names_with_DM->{$gene_name}}) {
			$h_genes_with_dm->{$gene_name}->{variants}->{$hgmd_id} = $h->{$hgmd_id};
		}
	}
	my $disease = $buffer->queryHgmd->get_gene_disease($gene_name);
	$h_genes_with_dm->{$gene_name}->{disease} = $disease;
	
	my @l_go_terms_acc = split('\|', $buffer->queryHgmd->get_gene_go_terms_acc($gene_name));
	my @l_go_terms_names = split('\|', $buffer->queryHgmd->get_gene_go_terms_name($gene_name));
	
	my $i = 0;
	foreach my $go_acc (@l_go_terms_acc) {
		my $go_name = $l_go_terms_names[$i];
		$h_from_go_terms->{$go_name}->{genes}->{$gene_name} = undef;
		$i++;
	} 
	
	foreach my $this_disease (split('\|', $disease)) {
		$h_from_d->{$this_disease}->{$gene_name} = undef;
	}
	
	my $concept_exp_inheritance = 1;
	$concept_exp_inheritance = undef if ($filter_inheritance);
	my $exp_inheritance = $buffer->queryHgmd->get_gene_expected_inheritance($gene_name);
	$concept_exp_inheritance = 1 if ($filter_inheritance and lc($exp_inheritance) eq lc($filter_inheritance));
	$h_from_exp_inheritance->{$exp_inheritance}->{$gene_name} = undef;
	unless ($concept_exp_inheritance) {
		delete $h_genes_with_dm->{$gene_name};
		next;
	}
	$h_genes_with_dm->{$gene_name}->{expected_inheritance} = $exp_inheritance;
}
$nb_genes_all = keys %{$h_genes_all};

my @headers0 = ('', 'HGMD DataBase');
my @headers = ( 'Gene Score', 'Locus', 'Gene Name', 'Phenotype(s)', 'Omim', 'Gtex', 'pLi', 'Panel(s)','HGMD: Concept(s)', 'DejaVu', 'HGMD: Total Mut', 'HGMD: Total New Mut' );
my @lFields = ( 'score', 'locus', 'name', 'phenotypes', 'omim', 'gtex', 'pli', 'panels', 'hgmd', 'variants', 'hgmd_mut', 'hgmd_new' );
my $last_genecode = $buffer->getQuery->getMaxGencodeVersion();
my $proj_name = $buffer->get_random_project_name_with_this_annotations_and_genecode();
my $project = $buffer->newProject( -name => $proj_name );

my $pheno_name;
if ($filter_concept) {
	my $h_pheno_infos = $buffer->queryPhenotype->getPhenotypeInfosFromHgmdConcept($filter_concept);
	confess("\n\nERROR: $filter_concept not found in our DB. die.\n\n") unless (exists $h_pheno_infos->{$filter_concept});
	$pheno_name = $h_pheno_infos->{$filter_concept}->{name};
	my @l_pheno;
	push (@l_pheno, $pheno_name);
	$project->{phenotypes} = \@l_pheno;
	$project->{phenotypes_object}->{$h_pheno_infos->{$filter_concept}->{phenotype_id}} = undef;
}


my ($h_lines_sorted, $h_panels_found, $h_concepts_genes_new);
my $nb_errors = 0;
my $out2;
foreach my $gname (sort keys %{$h_genes_with_dm}) {
	my $has_new_concept;
	my (@lConcepts, @lConceptsNew);
	my $concept_list = "<center><b><u>Concept(s) for gene $gname</b></u><br><br>";
	foreach my $concept_name (sort keys %{$h_genes_with_dm->{$gname}->{concepts}}) {
		if ($h_genes_with_dm->{$gname}->{concepts}->{$concept_name}->{new} == 1) {
			push(@lConceptsNew, $concept_name);
			$concept_list .= $concept_name." <b><i>New!</b></i><br>";
		}
		else {
			push(@lConcepts, $concept_name);
			$concept_list .= $concept_name."<br>";
		}
	}
	$concept_list .= "</center>";
	
	$nb_genes_selected++;
	my $disease = $h_genes_with_dm->{$gname}->{disease};
	my @lDiseases = split('\|', $disease);
	my $total_mut = $buffer->queryHgmd->get_gene_mut_total($gname);
	my $total_new = $buffer->queryHgmd->get_gene_new_mut_total($gname);
	my $total_new_DM = 0;
	$total_new_DM = scalar(keys %{$h_genes_with_dm->{$gname}->{variants}}) if (exists $h_genes_with_dm->{$gname}->{variants});
	my $exp_inheritance = $h_genes_with_dm->{$gname}->{expected_inheritance};
	my $url_dv = "./get_variants_from_all_my_projects.html?gene=$gname";
	my $cmd_link = qq{document.getElementById('input_search').setAttribute('value', '$gname');launch_for_gene('$gname');};
	my $html_gene;
	
	my $gene;
	eval { $gene = $project->newGene($gname); };
	if ($@ or not $gene) {
		my $refseq_id = $buffer->queryHgmd->get_refseq_from_gene_name($gname);
		my ($refseq, $tmp_id) = split('\.', $refseq_id);
		eval {
			my $tr = $project->newTranscript($refseq);
			my @lGenes = @{$tr->getGenes()};
			$gene = $lGenes[0] if (scalar @lGenes == 1);
		};
		if ($@) {}
	}
	if (not $gene) {	
		my $class = 'hgmd_gene hgmd_gene_new';
		foreach my $concept_name (@lConcepts) { $class .= ' '.$concept_name; }
		foreach my $concept_name (@lConceptsNew) {
			$class .= ' '.$concept_name;
			$class .= ' '.$concept_name.'_new';
			$h_concepts_genes_new->{$concept_name}->{$gname} = undef;
		}
		
		my $disease_hgmd = $buffer->queryHgmd->get_disease_from_gene_name($gname);
		my $out_line = $cgi->start_Tr({style=>"font-size:11px;max-height:60pxoverflow-y: auto;", class=>$class});
		$out_line .= $cgi->td({rowspan=>1, style=>"max-height:60px;overflow-y: auto;"}, '');
		$out_line .= $cgi->td({rowspan=>1, style=>"max-height:60px;overflow-y: auto;"}, qq{<i>Not Found Genecode <span style="color:red;font-size:12px;">$last_genecode</span></i>});
		$out_line .= $cgi->td({rowspan=>1, style=>"max-height:60px;overflow-y:auto;padding-left:20px;font-size:12px;"}, qq{<span style="color:black;"><i>$gname</i></span>});
		$out_line .= $cgi->td({rowspan=>1, style=>"max-height:60px;overflow-y: auto;"}, "<i>$disease_hgmd</i>");
		$out_line .= $cgi->td({rowspan=>1, style=>"max-height:60px;overflow-y: auto;"}, '');
		$out_line .= $cgi->td({rowspan=>1, style=>"max-height:60px;overflow-y: auto;"}, '');
		$out_line .= $cgi->td({rowspan=>1, style=>"max-height:60px;overflow-y: auto;"}, '');
		$out_line .= $cgi->td({rowspan=>1, style=>"max-height:60px;overflow-y: auto;"}, '');
		my $nb_concepts = scalar(@lConcepts);
		my $nb_concepts_new = scalar(@lConceptsNew);
		if ($nb_concepts_new > 0) {
			$nb_concepts += $nb_concepts_new;
			if (lc($filter_concept) eq 'all') {
				$out_line .= qq{<td><a class="btn btn-xs" onClick="document.getElementById('span_list_panels').innerHTML='$concept_list';dijit.byId('dialog_list_panels').show();" style="min-width:30px"><span class="badge" style="opacity:0.4;background-color:red;color:white;font-size:11px;"><strike>$nb_concepts</strike></span></a></td>};
			}
			elsif (lc($filter_concept)) {
				my $is_new;
				foreach my $this_concept (@lConceptsNew) {
					$is_new = 1 if (lc($this_concept eq lc($filter_concept)));
				}
				if ($is_new == 1) {
					$out_line .= qq{<td><center><a class="btn btn-xs" onClick="document.getElementById('span_list_panels').innerHTML='$concept_list';dijit.byId('dialog_list_panels').show();" style="min-width:30px"><span class="badge" style="opacity:0.4;background-color:red;color:white;font-size:11px;"><strike>$nb_concepts</strike></span></a></center></td>};
				}
				else {
					$out_line .= qq{<td><center><a class="btn btn-xs" onClick="document.getElementById('span_list_panels').innerHTML='$concept_list';dijit.byId('dialog_list_panels').show();" style="min-width:30px"><span class="badge" style="opacity:0.4;background-color:green;color:white;font-size:11px;"><strike>$nb_concepts</strike></span></a></center></td>};
				} 
			}
			else {
				$out_line .= qq{<td><center><a class="btn btn-xs" onClick="document.getElementById('span_list_panels').innerHTML='$concept_list';dijit.byId('dialog_list_panels').show();" style="min-width:30px"><span class="badge" style="opacity:0.4;background-color:red;color:white;font-size:11px;"><strike>$nb_concepts</strike></span></a></center></td>};
			}
		}
		else { 
			$out_line .= qq{<td><center><a class="btn btn-xs" onClick="document.getElementById('span_list_panels').innerHTML='$concept_list';dijit.byId('dialog_list_panels').show();" style="min-width:30px"><span class="badge" style="opacity:0.4;background-color:green;color:white;font-size:11px;"><strike>$nb_concepts</strike></span></a></center></td>};
		}
		$out_line .= qq{<td><center><a class="btn btn-xs disabled" style="min-width:30px"><span class="badge" style="opacity:0.4;background-color:green;color:white;font-size:11px;"><strike>$total_mut</strike></span></a></center></td>};
		$out_line .= qq{<td><center><a class="btn btn-xs disabled" style="min-width:30px"><span class="badge" style="opacity:0.4;background-color:green;color:white;font-size:11px;"><strike>$total_new</strike></span></a></center></td>};
		if ($total_new_DM and $total_new_DM > 0) {
			$out_line .= qq{<td><center><a class="btn btn-xs disabled" style="min-width:30px"><span class="badge" style="opacity:0.4;background-color:red;color:white;font-size:11px;"><strike>$total_new_DM</strike></span></a></center></td>};
		}
		else {
			$out_line .= qq{<td><center><i class="fa fa-minus"></i></center></td>};
		}
		$out_line .= $cgi->end_Tr();
		
		if ($sort and $sort eq 'locus') {
			$h_lines_sorted->{999999}->{$nb_errors}->{$nb_errors} = $out_line;
		}
		elsif ($sort and $sort eq 'gene_score') {
			$h_lines_sorted->{-999}->{$nb_errors} = $out_line;
		}
		else {
			$h_lines_sorted->{$gname} = $out_line;
		}
		$nb_errors++;
	}
	else {
		my ($hResGene, @lAllVar);
		my $external_name = $gene->external_name();
		$hResGene->{$gname}->{id} = $gname;
		$hResGene->{$gname}->{external_name} = $external_name;
		$hResGene->{$gname}->{pLI} = $gene->pLI();
		$hResGene->{$gname}->{omim_id} = $gene->omim_id();
		$hResGene->{$gname}->{omim_inheritance} = $gene->omim_inheritance();
		$hResGene->{$gname}->{variants} = \@lAllVar;
		my ($pheno,$nb_other_terms) = $gene->polyviewer_phentotypes();
		$hResGene->{$gname}->{phenotypes}->{pheno} = $pheno;
		$hResGene->{$gname}->{phenotypes}->{nb_other_terms} = $nb_other_terms;
		$hResGene->{$gname}->{specific_cmd} = '';
		$hResGene->{$gname}->{max_score} = $gene->score();
		my $description_gene = $gene->description();
		my $class = 'hgmd_gene hgmd_gene_new';
		eval {
			foreach my $panel (@{$gene->getPanels()}) {
				my $name = lc($panel->name());
				$name =~ s/ /_/g;
				unless (exists $hResGene->{$gname}->{panels}->{$name}) {
					$hResGene->{$gname}->{panels}->{$name}->{phenotype} = $panel->getPhenotypes()->[0]->name();
					$h_panels_found->{$name}++;
					my $value = 'panel_'.$name;
					$class .= ' '.$value;
				}
			}
		};
		if ($@) { $hResGene->{$gname}->{panels} = undef; }
		my $panel_id = 'p_'.$gname.'_hgmd';
		$hResGene->{$gname}->{uid} = $panel_id;
		
		my $html_gene = html_polygenescout::print_gene_basic_tables_infos($hResGene->{$gname});
		foreach my $concept_name (@lConcepts) { $class .= ' '.$concept_name; }
		foreach my $concept_name (@lConceptsNew) {
			$class .= ' '.$concept_name;
			$class .= ' '.$concept_name.'_new';
			$h_concepts_genes_new->{$concept_name}->{$gname} = undef;
		}
		my $out_line = $cgi->start_Tr({style=>"font-size:11px;max-height:60px;overflow-y: auto;", class=>$class});
		
		
		
		my $gene_score;
		if ($filter_concept) {
			$gene_score = get_gene_score($gene);
			$out_line .= html_polygenescout::print_gene_score($gene_score);
		}
		$out_line .= html_polygenescout::print_locus($gene);
		$out_line .= html_polygenescout::print_gene_basic_tables_infos($hResGene->{$gname}, $pheno_name);
		
		my $nb_concepts = scalar(@lConcepts);
		my $nb_concepts_new = scalar(@lConceptsNew);
		if ($nb_concepts_new > 0) {
			$nb_concepts += $nb_concepts_new;
			if (lc($filter_concept) eq 'all') {
				$out_line .= qq{<td><a class="btn btn-xs" onClick="document.getElementById('span_list_panels').innerHTML='$concept_list';dijit.byId('dialog_list_panels').show();" style="min-width:30px"><span class="badge" style="background-color:red;color:white;font-size:11px;">$nb_concepts</span><br><span style="color:red;"><i><b>New</i></b></span></a></td>};
			}
			elsif (lc($filter_concept)) {
				my $is_new;
				foreach my $this_concept (@lConceptsNew) {
					$is_new = 1 if (lc($this_concept eq lc($filter_concept)));
				}
				if ($is_new == 1) {
					$out_line .= qq{<td><center><a class="btn btn-xs" onClick="document.getElementById('span_list_panels').innerHTML='$concept_list';dijit.byId('dialog_list_panels').show();" style="min-width:30px"><span class="badge" style="background-color:red;color:white;font-size:11px;">$nb_concepts</span><br><span style="color:red;"><i><b>New</i></b></span></a></center></td>};
				}
				else {
					$out_line .= qq{<td><center><a class="btn btn-xs" onClick="document.getElementById('span_list_panels').innerHTML='$concept_list';dijit.byId('dialog_list_panels').show();" style="min-width:30px"><span class="badge" style="background-color:green;color:white;font-size:11px;">$nb_concepts</span></a></center></td>};
				} 
			}
			else {
				$out_line .= qq{<td><center><a class="btn btn-xs" onClick="document.getElementById('span_list_panels').innerHTML='$concept_list';dijit.byId('dialog_list_panels').show();" style="min-width:30px"><span class="badge" style="background-color:red;color:white;font-size:11px;">$nb_concepts</span><br><span style="color:red;"><i><b>New</i></b></span></a></center></td>};
			}
		}
		else { 
			$out_line .= qq{<td><center><a class="btn btn-xs" onClick="document.getElementById('span_list_panels').innerHTML='$concept_list';dijit.byId('dialog_list_panels').show();" style="min-width:30px"><span class="badge" style="background-color:green;color:white;font-size:11px;">$nb_concepts</span></a></center></td>};
		}
		$out_line .= html_polygenescout::print_link_dejavu($cmd_link);
		
		$out_line .= html_polygenescout::print_hgmd_total_mut($total_mut, $external_name);
		$out_line .= html_polygenescout::print_hgmd_total_new($total_new, $external_name);
		
		
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
				$out2 .= $h_lines_sorted->{$chr_id}->{$start}->{$end};
			}
		}
	}
}
elsif ($sort and $sort eq 'gene_score') {
	foreach my $gene_score (sort {$b <=> $a} keys %$h_lines_sorted) {
		foreach my $external_name (sort keys %{$h_lines_sorted->{$gene_score}}) {
			$out2 .= $h_lines_sorted->{$gene_score}->{$external_name};
		}
	}
}
else {
	foreach my $external_name (sort keys %$h_lines_sorted) {
		$out2 .= $h_lines_sorted->{$external_name}; 
	}
}

$out2 .= "</tbody>";
$out2 .= "</table>";
$out2 .= "</div>";
$out2 .= "<br>";

my $hgmd_release_now = $buffer->queryHgmd->database();
my $hgmd_release_old = $buffer->queryHgmd->previous_database();

if ($nb_genes_selected == 0) {
	$out2 .= qq{<span style="font-size:14px;"><b><u>No new DM variants found for concept $filter_concept from $hgmd_release_old to $hgmd_release_now</b></u></span>};
}

my $out_th0;
$out_th0 .= qq{<div class="row" style="padding:4px;">};

my $text_hgmd_compare = "<center><b>Genes selected with new DM variants (From <font color='red'>$hgmd_release_old</font> To <font color='green'>$hgmd_release_now</font>";
$text_hgmd_compare .= " + ONLY IF IS A NEW CONCEPT" if ($filter_only_new);
$text_hgmd_compare .= ")</b></center>";

$out_th0 .= qq{<div class="col-sm-4">};
$out_th0 .= qq{<table>};
$out_th0 .= qq{<td>};
$out_th0 .= qq{<div class="input-group">};
$out_th0 .= qq{<label for="inputHGMDVersions">HGMD Release Database</label>};
$out_th0 .= qq{<select class="form-control" id="form_hgmd_versions">};
foreach my $old_version (reverse sort keys %{$buffer->queryHgmd->getHashOldDatabases()}) {
	my $selected;
	$selected = 'selected' if ($old_version eq $hgmd_release_old);
	my $old_version_text = $old_version;
	$old_version_text .= '-2019.1' if ($old_version_text eq 'hgmd_pro');
	$out_th0 .= qq{<option value='$old_version' $selected>From release $old_version_text </option>};
}
$out_th0 .= qq{</select>};
$out_th0 .= qq{</div>};
$out_th0 .= qq{</td>};
$out_th0 .= qq{<td style="vertical-align:bottom;">};
$out_th0 .= qq{<div>};
$out_th0 .= qq{<button onclick="change_hgmd_version();" type="submit" style="border:solid 1px black;" class="btn btn-secondary">Apply</button>};
$out_th0 .= qq{</div>};
$out_th0 .= qq{</td>};
$out_th0 .= qq{</table>};
$out_th0 .= qq{</div>};

$out_th0 .= qq{<div class="col-sm-4">};
$out_th0 .= qq{<div class="input-group" style="width:100%">};
$out_th0 .= qq{<label for="inputGroupHgmdPanels">PolyWeb Panels filter</label>};
$out_th0 .= qq{<select class="form-control" id="form_panels">};
$out_th0 .= qq{<option style="color:green;font-weight:bold;" value="hgmd_gene"><b>ALL Panels</b></option>};
foreach my $panel (sort keys %{$h_panels_found}) {
	my $name = lc($panel);
	$name =~ s/ /_/g;
	my $value = 'panel_'.$name;
	$name = ucfirst($name);
	my $nb_genes = $h_panels_found->{$panel};
	$out_th0 .= qq{<option value='$value'><span>$name [$nb_genes genes]</span></option>};
}
$out_th0 .= qq{</select>};
$out_th0 .= qq{</div>};
$out_th0 .= qq{</div>};

$out_th0 .= qq{<div class="col-sm-4">};
$out_th0 .= qq{<form>};
$out_th0 .= qq{<div class="form-group">};
$out_th0 .= qq{<label for="form_search">Quick Search</label>};
$out_th0 .= qq{<input class="form-control" type="search" name="search" value="" id="id_search" />};
$out_th0 .= qq{</div>};
$out_th0 .= qq{</form>};
$out_th0 .= qq{</div>};

$out_th0 .= qq{</div>};

my $out1;
$out1 .= $out_th0;
$out1 .= $text_hgmd_compare;
$out1 .= "<br><div style='text-align:center;font-size:12px;width:100%;overflow-x:scroll;'>";
$out1 .= "<table data-filter-control='true' data-searchable='true'  data-toggle='table' data-show-extended-pagination='true' data-cache='false' data-pagination-loop='false' data-virtual-scroll='true' data-pagination-v-align='bottom' data-pagination-pre-text='Previous' data-pagination-next-text='Next' data-pagination='true' data-page-size='20' data-page-list='[20, 50, 100, 200, 300]' data-resizable='true' id='table_hgmd_genes' class='table table-striped table-bordered' >";
$out1 .= "<thead>";
#$out1 .= $cgi->start_Tr({style=>"font-size:11px;"});
#$out1 .= $cgi->th({colspan=>8, style=>"background-color:white;font-size:11px;color:#607D8B;"}, "<b>GENE Informations</b>");
#$out1 .= $cgi->th({colspan=>1, style=>"background-color:white;font-size:11px;color:#607D8B;"}, "<center><b>HGMD</b></center>");
#$out1 .= $cgi->th({colspan=>1, style=>"background-color:white;font-size:11px;color:#607D8B;"}, "<center><b>PolyWeb</b></center>");
#$out1 .= $cgi->th({colspan=>2, style=>"background-color:white;font-size:11px;color:#607D8B;"}, "<center><b>HGMD Informations</b></center>");
#$out1 .= $cgi->end_Tr();
$out1 .= $cgi->start_Tr({style=>"font-size:11px;"});
my $z = 0;
foreach my $h (@headers){
	next unless $h;
	my $type =  $lFields[$z];
	if ($type eq 'name') {
		$out1 .= qq{<th data-searchable="true" data-filter-control-placeholder='Ex: BRCA1' data-filter-control='input' data-field='name'><center><b>$h</b></center></th>};
	}
	elsif ($type eq 'phenotypes') {
		$out1 .= qq{<th data-searchable="true" data-filter-control-placeholder='Ex: heart, autism' data-filter-control='input' data-field='phenotypes'><center><b>$h</b></center></th>};
	}
	elsif ($type eq 'locus') {
		$out1 .= qq{<th data-searchable="true" data-filter-control-placeholder='Ex: 12:48366750-48398337' data-filter-control='input' data-field='locus'><center><b>$h</b></center></th>};
	}
	else {
		$out1 .= qq{<th data-field="$type"><center><b>}.ucfirst($h).qq{</b></center></th>};
	}
	$z++;
#	$out1 .=  $cgi->th({style=>"background-color:white;font-size:11px;color:#607D8B;"}, "<center><b>".ucfirst($h)."</b></center>");
#	$z++;
}
$out1 .= $cgi->end_Tr();
$out1 .= "</thead>";
$out1 .= "<tbody>";

my $hRes;
$hRes->{html} = $out1.$out2;
$hRes->{nb_genes_errors} = $nb_errors if ($nb_errors > 0);
if ($nb_genes_selected == $nb_genes_all) { $hRes->{nb_genes_selected} = ''; }
else { $hRes->{nb_genes_selected} = $nb_genes_selected.' selected'; }
$hRes->{nb_genes_total} = $nb_genes_all.' genes';
$hRes->{filter_concept} = '';
$hRes->{filter_concept} = 'only '.ucfirst($filter_concept) if ($filter_concept);
$hRes->{filter_expected_inheritance} = '';
$hRes->{filter_expected_inheritance} = 'only '.uc($filter_inheritance) if ($filter_inheritance);
$hRes->{hgmd_version} = $hgmd_release_now;

my $json_encode = encode_json $hRes;
print $cgi->header('text/json-comment-filtered');
print $json_encode;
exit(0);



sub get_gene_score {
	my ($gene) = @_;
	return $gene->raw_score();
}
