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
use xls_export;
use session_export;

require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/html_polygenescout.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";

my $io = IO::Handle->new();
$io->autoflush(1);

my $max_dejavu = 999999999999;
my $max_dejavu_ho = 999999999999;
my $max_gnomad = 999999999999;
my $max_gnomad_ho = 999999999999;

my $cgi = new CGI();
my $user_name = $cgi->param('user');
my $pwd = $cgi->param('pwd');
my $gene_id = $cgi->param('gene');
my $gene_id_alt = $cgi->param('gene_alt');
my $only_transcript = $cgi->param('only_transcript');
$max_dejavu = $cgi->param('dejavu');
$max_dejavu_ho = $cgi->param('dejavu_ho');
$max_gnomad = $cgi->param('gnomad');
$max_gnomad_ho = $cgi->param('gnomad_ho');
my $filters_cons = $cgi->param('filters_cons');
my $only_ill = $cgi->param('only_ill');
my $only_my_projects = $cgi->param('only_my_projects');
my $models = $cgi->param('models');
my $use_session_id = $cgi->param('session_id');
my $projects_excluded = $cgi->param('proj_excluded');
my $filter_perc_allelic_min = $cgi->param('min_perc_all');
my $filter_perc_allelic_max = $cgi->param('max_perc_all');
my $only_interval = $cgi->param('only_interval');
my $export_xls = $cgi->param('export_xls');
my $only_pat_with_var = $cgi->param('only_pat_with_var');
my $only_pat_with_var_he = $cgi->param('only_pat_with_var_he');
my $only_pat_with_var_ho = $cgi->param('only_pat_with_var_ho');
#my $only_project = $cgi->param('project');
#my $only_patient = $cgi->param('patient');
#if ($only_project and $only_patient) {
#	$only_ill = undef;
#	$only_my_projects = 1;
#}

my ($only_project, $only_patient);
my ($project_dejavu, $gene_with_partial_transcrit);

$only_pat_with_var =~ s/-/_/g if ($only_pat_with_var);

$filter_perc_allelic_max = undef if ($filter_perc_allelic_max == 0);

my $h_models;
foreach my $model (split(',', $models)) {
	$h_models->{$model} = undef;
}
exit(0) unless $h_models;


supressCoreFilesFound();

my $fork = 5;
$user_name = lc($user_name);
my $buffer_init = new GBuffer;
my $can_use_hgmd = $buffer_init->hasHgmdAccess($user_name);

my $h_filters_cons;
if ($filters_cons) {
	foreach my $cons (split(',', $filters_cons)) {
		$h_filters_cons->{lc($cons)} = undef;
	}
}
my $h_annot_categories = get_hash_annot_categories();

my $hRes;
$hRes->{ok} = 1;
my @lItemsProjects;
my $hProjects = get_hash_users_projects($user_name, $pwd);
foreach my $project_name_excluded (split(',', $projects_excluded)) {
	delete $hProjects->{$project_name_excluded};
}

my @headers_validations = ("#", "varsome","alamut variant","var_name","projects / patients","gnomad","deja_vu","table_validation","table_transcript");
my @header_transcripts = ("consequence","enst","nm","ccds","appris","exon","nomenclature","codons","codons_AA", "polyphen","sift","ncboost","cadd","revel","dbscsnv","spliceAI");
my ($hVariantsIdsDejavu, $hVariantsDetails);
my ($hResGene, $hResVariants, $hResVariantsIds, $hResVariantsListPatients, $hResProjectProblem, $hResVariants_byScores, $hResVariantsRatioAll, $hResVariantsModels);
my $project_init_name;

#project preload_patients test connection DB
my $dbh_init = $buffer_init->dbh();
my $query_init = $buffer_init->getQuery();
my $project_init;

my ($h_projects_name_with_capture,$h_proj_pat_all,$h_proj_pat_ill);
if ($only_my_projects and $only_my_projects ne '1') {
	if ($only_my_projects eq 'only_genomes') {
		foreach my $p_name (@{$query_init->getSimilarProjectsByAnalyse('genome')}) {
			$h_projects_name_with_capture->{$p_name} = undef;
		}
	}
	elsif ($only_my_projects eq 'only_exomes') {
		foreach my $p_name (@{$query_init->getSimilarProjectsByAnalyse('exome')}) {
			$h_projects_name_with_capture->{$p_name} = undef;
		}
	}
	else {
		foreach my $p_name (@{$query_init->listAllProjectsNameByCapture($only_my_projects)}) {
			$h_projects_name_with_capture->{$p_name} = undef;
		}
	}
}
my @lProjectNames = reverse sort keys %{$hProjects};

foreach my $project_name (@lProjectNames) {
	my $h;
	$h->{name} = $project_name;
	$h->{description} = $hProjects->{$project_name}->{description};
	$h->{id} = $hProjects->{$project_name}->{id};
	if ($only_ill) {
		$h_proj_pat_ill->{$project_name} = get_hash_patients_ill_from_project_name($dbh_init, $query_init, $project_name);
	}
	$h_proj_pat_all->{$project_name} = get_hash_patients_all_from_project_name($dbh_init, $query_init, $project_name);
	push(@lItemsProjects, $h);
}

$project_init_name = $buffer_init->get_random_project_name_with_this_annotations_and_genecode();
#$project_init = $buffer_init->newProjectCache( -name => $project_init_name);

#TODO: here HG38
$project_init = $buffer_init->newProjectCache( -name => 'NGS2024_7792');

my $genomeFai_init = $project_init->getGenomeFai();
my $gene_init;
eval { $gene_init = $project_init->newGene($gene_id); };
if ($@) { $gene_init = $project_init->newGene($gene_id_alt); }
my $gene_init_id_for_newgene = $gene_id;
my ($gene_ensg, $gene_chr_id);
eval { ($gene_ensg, $gene_chr_id) = split('_', $gene_init->id()); };
if ($@) {
	$gene_init_id_for_newgene = $gene_id_alt;
	$gene_init = $project_init->newGene($gene_id_alt);
	($gene_ensg, $gene_chr_id) = split('_', $gene_init->id());
}

my ($use_start, $use_end);
if ($only_interval) {
	$only_interval =~ s/chr//g;
	$only_interval =~ s/:/_/g;
	$only_interval =~ s/-/_/g;
	my ( $locus_chr, $locus_start, $locus_end ) = split( '_', $only_interval );
	$locus_end = $locus_start + 1 unless ($locus_end);
	if ($locus_start > $locus_end){
		my $t = $locus_start;
		$locus_start = $locus_end;
		$locus_end = $t;
	}
	confess("\n\nERROR chromosome problem... Die...\n\n") if ($gene_init->getChromosome->id() ne $locus_chr);
	$use_start = $locus_start;
	$use_end = $locus_end;
}
else {
	$use_start = $gene_init->start();
	$use_end = $gene_init->end();
}
my $use_locus = $gene_chr_id.'_'.$use_start.'_'.$use_end;

my $chr_init = $project_init->getChromosome($gene_chr_id);
my $chr_init_id = $chr_init->id();
$gene_init->getTranscripts();
my $gene_init_id = $gene_init->id();

$project_init->getChromosomes();
my $hchr_init = $project_init->{chromosomes_object};
$project_init->{genes_object}->{$gene_init->id}++;
my $hgene_init = $project_init->{genes_object};
my $objects_init = $project_init->{objects};
foreach my $type (keys %$objects_init) {
	foreach my $id (keys %{$objects_init->{$type}}) {
		delete $objects_init->{$type}->{$id}->{buffer};
		delete $objects_init->{$type}->{$id}->{project};
	}
}

my ($filters_saved, $hResVariants_loaded, $h_count_loaded, $gene_id_loaded);
print $cgi->header('text/json-comment-filtered');
print "{\"progress\":\".";
my $h_count;
($h_count, $hResVariants, $hVariantsDetails, $hResVariantsModels) = get_variants_infos_from_projects($hResVariants_loaded, $hVariantsDetails, $hResVariantsModels, $use_locus, $only_transcript);
my $nb_var_after = scalar(keys %$hVariantsDetails);


#warn "\n\n\n";
#warn "EXPORT";



my $fsize = "font-size:10px";
my $value_red = 14;
my $value_coral = 12;
my $value_orange = 8;
my $value_yellow = 5;
my $color_red = '#CE0000';
my $color_coral = 'coral';
my $color_orange = '#EFB73E';
my $color_yellow = 'yellow';

my $session_id = save_export_xls();
print "|";
export_html($hResVariants, $session_id);
supressCoreFilesFound();
exit(0);



####################
### METHODS
####################



sub supressCoreFilesFound {
	my $cmd = "rm $Bin/core.*";
	warn $cmd;
	`$cmd`;
}

sub save_export_xls {
	print '.saveExport.';
	$project_dejavu->getProject->buffer->dbh_deconnect();
	$project_dejavu->getProject->buffer->dbh_reconnect();
#	delete $project_dejavu->{rocksPartialTranscripts};
	my $h_xls_args;
	print '_save_xls_';
	$project_dejavu->cgi_object(1);
	my (@lVarObj, $h_pubmed);
	foreach my $pos (sort keys %{$hResVariants}) {
		foreach my $var_id (sort keys %{$hResVariants->{$pos}}) {
			$project_dejavu->print_dot(50);
			my $v = $project_dejavu->_newVariant($var_id);
			if ($gene_with_partial_transcrit) {
				$v->variationTypeInterface($gene_init);
			}
			push(@lVarObj, $v);
			if ($v->hgmd_details() and not exists $h_pubmed->{$v->id}) {
				$h_pubmed->{$v->id()}->{$v->hgmd_details->{pmid}}->{url} = "https://www.ncbi.nlm.nih.gov/pubmed/".$v->hgmd_details->{pmid};
				$h_pubmed->{$v->id()}->{$v->hgmd_details->{pmid}}->{title} = $v->hgmd_details->{title};
				foreach my $pubmed_id (keys %{$v->hgmd_details->{pubmed}}) {
					if (exists $v->hgmd_details->{pubmed}->{$pubmed_id}->{title}) {
						$h_pubmed->{$v->id()}->{$pubmed_id}->{url} = "https://www.ncbi.nlm.nih.gov/pubmed/".$pubmed_id;
						$h_pubmed->{$v->id()}->{$pubmed_id}->{title} = $v->hgmd_details->{pubmed}->{$pubmed_id}->{title};
					}
				}
			}
		}
	}
	my $xls_export = new xls_export();
	$xls_export->title_page('GeneScout_'.$gene_init->external_name().'.xls');
	$project_dejavu->buffer->dbh_deconnect();
	$project_dejavu->buffer->dbh_reconnect();
	
	$xls_export->store_variants_infos(\@lVarObj, $project_dejavu);
	my ($h_patients, $h_row_span);
	foreach my $chr_id (keys %{$xls_export->{hash_variants_global}}) {
		foreach my $var_id (keys %{$xls_export->{hash_variants_global}->{$chr_id}}) {
			foreach my $project_name (keys %{$hResVariantsListPatients->{$var_id}}) {
				foreach my $pat_name (keys %{$hResVariantsListPatients->{$var_id}->{$project_name}}) {
					my $phenotypes = $hResVariantsListPatients->{$var_id}->{$project_name}->{$pat_name}->{'values'}->{'phenotypes'};
					my $description = $hResVariantsListPatients->{$var_id}->{$project_name}->{$pat_name}->{'values'}->{'description'};
					my $fam = $hResVariantsListPatients->{$var_id}->{$project_name}->{$pat_name}->{'values'}->{'fam'};
					my $name = $hResVariantsListPatients->{$var_id}->{$project_name}->{$pat_name}->{'values'}->{'name'};
					my $sex = $hResVariantsListPatients->{$var_id}->{$project_name}->{$pat_name}->{'values'}->{'sex'};
					my $status = $hResVariantsListPatients->{$var_id}->{$project_name}->{$pat_name}->{'values'}->{'status'};
					my $parent_child = $hResVariantsListPatients->{$var_id}->{$project_name}->{$pat_name}->{'values'}->{'parent_child'};
					my $sex_status_icon = $hResVariantsListPatients->{$var_id}->{$project_name}->{$pat_name}->{'values'}->{'sex_status_icon'};
					my $perc = $hResVariantsListPatients->{$var_id}->{$project_name}->{$pat_name}->{'values'}->{'percent'};
					my $model = $hResVariantsListPatients->{$var_id}->{$project_name}->{$pat_name}->{'values'}->{'model'};
					my $he_ho = $hResVariantsListPatients->{$var_id}->{$project_name}->{$pat_name}->{'values'}->{'he_ho'};
					$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'variation'} = $var_id;
					$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'project'} = $project_name;
					$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'description'} = $description;
					$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'phenotypes'} = $phenotypes;
					$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'fam'} = $fam;
					$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'name'} = $name;
					$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'sex'} = $sex;
					$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'status'} = $status;
					$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'parent_child'} = $parent_child;
					$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'sex_status_icon'} = $sex_status_icon;
					$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'perc'} = $perc.'%';
					$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'model'} = $model;
					$h_patients->{$var_id}->{$project_name}->{$pat_name}->{'he_ho'} = $he_ho;
				}	
			}
		}
	}
	
	print "|";
	$project_dejavu->buffer->dbh_deconnect();
	$project_dejavu->buffer->dbh_reconnect();
#	delete $project_dejavu->{rocksPartialTranscripts};
	print "|";
	$xls_export->store_specific_infos('projects_patients_infos', $h_patients);
	print "|";
	$xls_export->store_specific_infos('variants_pubmed', $h_pubmed);
	print "|";
	my $session_id = $xls_export->save();
	print "|";
	$project_dejavu->buffer->dbh_deconnect();
#	delete $project_dejavu->{rocksPartialTranscripts};
#	warn $session_id;
#	die;
	return $session_id;
}

sub export_html {
	my ($hResVariants, $session_id) = @_;
	
	my $class;
	$class->{rowspan} -= 1;
	$class->{rowspan} = 1 if $class->{rowspan} <=0;
	$class->{style} = "min-width:10%;padding:1px";
	
	my $out2 .= $cgi->start_div();
	
	$out2 .= qq{<table data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-v-align="both" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="100" data-page-list="[25, 50, 100, 200, 300]" data-resizable='true' id='table_variants' class='table table-striped' style='font-size:13px;'>};
	$out2 .= "<thead>";
	$out2 .= $cgi->start_Tr({style=>"background-color:#E9DEFF;$fsize"});
	foreach my $h (@headers_validations){
		my $cat = ucfirst($h);
		if ($h eq "projects / patients") {
#			if ($only_project) {
#				$out2 .= qq{<th data-field="$h" data-filter-control="input" data-filter-control-placeholder="NGS2022_4000 / Denovo">$only_project / $only_patient FAM</th>};
#			}
#			else {
				$out2 .= qq{<th data-field="$h" data-filter-control="input" data-filter-control-placeholder="NGS2022_4000 / Denovo">$cat / IGV / Alamut BAM</th>};
#			}
		}
		else {
			if (lc($cat) eq '#' or lc($cat) eq 'varsome' or lc($cat) eq 'alamut variant' or lc($cat) eq 'gnomad' or lc($cat) eq 'deja_vu') {
				$out2 .= qq{<th data-field="$h" data-sortable="true">$cat</th>};
			}
			elsif (lc($cat) eq 'table_validation') {
				$out2 .= qq{<th data-field="$h" data-filter-control="input" data-filter-control-placeholder="DM / Pathogenic">HGMD / Clinvar / Local</th>};
			}
			elsif (lc($cat) eq 'table_transcript') {
				$out2 .= qq{<th data-field="$h" data-filter-control="input" data-filter-control-placeholder="Stop / c.73A>C / exon2 / ENST00000145855">Annotations</th>};
			}
			elsif (lc($cat) eq 'var_name') {
				$out2 .= qq{<th data-field="$h" data-filter-control="input" data-sortable="true" data-filter-control-placeholder=7-15456-A-G">Annotations</th>};
			}
			else {
				$out2 .= qq{<th data-field="$h" data-filter-control="input" data-sortable="true">$cat</th>};
			}
		}
	}
	$out2 .= $cgi->end_Tr();
	$out2 .= "</thead>";
	$out2 .= "<tbody>";
	
	#warn "\n\n";
	#warn 'Nb Var Init: '.scalar keys %$hResVariants;
	
	my @lTrLines;
	my $nb_var = 0;
	my $nb_var_filtred = 0;
	
	
#	warn Dumper $hResVariants;
	
	foreach my $pos (sort keys %{$hResVariants}) {
		foreach my $var_id (sort keys %{$hResVariants->{$pos}}) {
			$nb_var++;
			next if not exists $hResVariantsListPatients->{$var_id};
			next if scalar keys %{$hResVariantsListPatients->{$var_id}} == 0;
						
			my $hvariation->{html} = $hResVariants->{$pos}->{$var_id};
			unless ($can_use_hgmd) {
				$hvariation->{html}->{hgmd} = qq{<span class="glyphicon glyphicon-ban-circle" aria-hidden="true" style='font-size:12px;color:black;'></span>};
			}
			my $class_default;
			$class_default->{style} = "max-width:350px;overflow-x:auto;vertical-align:middle;padding:5px;";
			
			my @l_pat;
			foreach my $project_name (keys %{$hResVariantsListPatients->{$var_id}}) {
				foreach my $pat_name (keys %{$hResVariantsListPatients->{$var_id}->{$project_name}}) {
					my $html_table = $hResVariantsListPatients->{$var_id}->{$project_name}->{$pat_name}->{html};
					my $percent_allele = $hResVariantsListPatients->{$var_id}->{$project_name}->{$pat_name}->{percent_allele};
					if ($filter_perc_allelic_max) {
						push(@l_pat, $html_table) if ($percent_allele < $filter_perc_allelic_max);
					}
					else {
						push(@l_pat, $html_table);
					}
				}	
			}
			
			next unless (@l_pat);
			my $nb_pat = scalar(@l_pat);
			my $table_trio;
			if (not $only_project and not $only_patient) {
				$table_trio .= qq{<div style="text-align:center;font-size:9px;padding-bottom:5px;"><div style="border:solid 0.5px black;"><b><i>Nb Patients: $nb_pat</b></i></div></div>} if ($nb_pat > 1);
			}
			$table_trio .= qq{<div style="max-height:170px;max-width:400px;overflow-y:auto;">};
			
			if ($only_project and $only_patient) { $table_trio .= join("", @l_pat); }
			else { $table_trio .= join("<br>", @l_pat); }
			$table_trio .= qq{</div>};
			
			my $out = $cgi->start_Tr();
			foreach my $h (@headers_validations){
				if ($h eq "trio" or "table_transcript"){
					$class->{style} = "min-width:200px;max-width:750px;max-height:200px;overflow-x:auto;vertical-align:middle;padding:5px;white-space: nowrap;";
				}
				elsif ($h eq 'var_name') {
					$class->{style} = "max-width:130px;overflow-x:auto;vertical-align:left;padding:5px;";
				}
				elsif ($h eq "projects / patients") {
					$class->{style} = "width:400px;overflow-x:auto;vertical-align:middle;padding:5px;";
				}
				else {
					$class->{style} = "vertical-align:middle;padding:5px;";
				}
				
				if ($h eq "projects / patients") {
					$out .= $cgi->td($class, $table_trio);
				}
				elsif ($h eq "gnomad") {
					$out .= $cgi->td($class_default,$hVariantsDetails->{$var_id}->{table_gnomad});
				}
				elsif ($h eq "deja_vu") {
					$out .= $cgi->td($class_default,$hVariantsDetails->{$var_id}->{table_dejavu});
				}
				elsif ($h eq 'table_validation') {
					$out .= $cgi->td($class_default,$hVariantsDetails->{$var_id}->{table_validation});
				}
				elsif ($h eq 'table_transcript') {
					$out .= qq{<div style="max-height:200px;overflow-y:auto;max-width:220px;overflow-x:auto;">};
					$out .= $cgi->td($class, $hVariantsDetails->{$var_id}->{table_transcript});
					$out .= qq{</div>};
				}
				elsif ($h eq '#') {
					$out.= $cgi->td($class_default,$nb_var);
				}
				elsif ($h eq 'alamut variant') {
					$out .= $cgi->td($class_default, "<center>".$hVariantsDetails->{$var_id}->{alamut_link_variant}."</center>");
				}
				elsif ($h eq 'var_name') {
					$out .= $cgi->td($class_default, "<center>".$hVariantsDetails->{$var_id}->{table_vname}."</center>");
				}
				elsif ($h eq 'varsome') {
					$out .= $cgi->td($class_default, "<center>".$hVariantsDetails->{$var_id}->{table_varsome}."</center>");
				}
				else {
					$out.= $cgi->td($class_default,$hvariation->{html}->{$h});
				}
			}
			$out .= $cgi->end_Tr();
			push(@lTrLines, $out);
		}
	}
	
	if ($gene_init->strand() == 1) { foreach my $line (@lTrLines) { $out2 .= $line; } }
	else { foreach my $line (reverse @lTrLines) { $out2 .= $line; }  }
	$h_count->{total_pass} = scalar(@lTrLines);
	$gene_init = undef;
	$gene_init = $project_dejavu->newGene($gene_init_id) unless ($gene_init);
	$gene_init->hgmd();
	$gene_init->omim();
	$gene_init->is_omim_morbid();

	$hRes->{html_gene} = get_html_gene($gene_init_id, $gene_init, $h_count);
	if ($nb_var_filtred > 0) {
		$out2 .= $cgi->start_Tr();
		$out2 .= "<td colspan='4' style='color:red;font-size:12px;'><b><u>and $nb_var_filtred variants fitred...</b></u></td>";
		$out2 .= $cgi->end_Tr();
	}
	if ($hResProjectProblem) {
		$out2 .= $cgi->start_Tr();
		$out2 .= "<td colspan='4' style='color:red;font-size:12px;'><center>---------------------------------</center></td>";
		$out2 .= $cgi->end_Tr();
		foreach my $proj_name (keys %$hResProjectProblem) {
			my $type_pb = $hResProjectProblem->{$proj_name};
			$out2 .= $cgi->start_Tr();
			$out2 .= "<td colspan='4' style='color:red;font-size:12px;'><b><u>Problem with $proj_name ($type_pb)...</b></u></td>";
			$out2 .= $cgi->end_Tr();
		}
	}
	$out2 .= "</tbody>";
	$out2 .= "</table>";
	$out2 .= "</div>";
	$out2 .= "<br>";
	
	
	$hRes->{html_variants} = $out2;
	$hRes->{project_init_name} = $project_init_name;
	$hRes->{gencode_version} = $project_dejavu->gencode_version;
	$hRes->{session_id} = $session_id;
	$hRes->{html_page_title} = 'DV '.$gene_init->external_name();
	$hRes->{html_source} = 'PolyDejaVu';
	$hRes->{html_title} = 'Gene '.$gene_init->external_name();
	
	
	my $h_annot_categories = get_hash_annot_categories();
	foreach my $this_annot (keys %{$h_annot_categories}) {
		$this_annot =~ s/ /_/g;
		if (exists $h_filters_cons->{lc($this_annot)}) {
			$h_annot_categories->{'var_'.lc($this_annot)} = 1;
		}
		else {
			$h_annot_categories->{'var_'.lc($this_annot)} = 0;
		}
	}
	$h_annot_categories->{'dejavu_'.$max_dejavu} = 1 if ($max_dejavu);
	$h_annot_categories->{'dejavu_ho_'.$max_dejavu_ho} = 1 if ($max_dejavu_ho);
	$h_annot_categories->{'gnomad_'.$max_gnomad} = 1 if ($max_gnomad);
	$h_annot_categories->{'gnomad_ho_'.$max_gnomad_ho} = 1 if ($max_gnomad_ho);
	$h_annot_categories->{'only_ill_patients'} = 1 if ($only_ill);
	if ($only_my_projects) {
		if ($only_my_projects eq '1') { $h_annot_categories->{'only_my_projects'} = 1; }
		else { $h_annot_categories->{"only_my_projects <span style='color:red;'>$only_my_projects</span>"} = 1; }
	}
	if ($models) {
		$models =~ s/,/, /g;
		$h_annot_categories->{'models '.$models} = 1;
	}
	$h_annot_categories->{'min_perc_all_'.$filter_perc_allelic_min} = 1 if ($filter_perc_allelic_min);
	$h_annot_categories->{'max_perc_all_'.$filter_perc_allelic_max} = 1 if ($filter_perc_allelic_max);
	
	my @lHeHo;
	push(@lHeHo, "<span style='color:green;'>He</span>") if ($only_pat_with_var_he);
	push(@lHeHo, "<span style='color:green;'>Ho</span>") if ($only_pat_with_var_ho);
	$h_annot_categories->{"patient_has_variant <b><span style='color:red'>".$only_pat_with_var."</span></b> (".join('+', @lHeHo).")"} = 1 if ($only_pat_with_var);
	
	$hRes->{hash_filters} = $h_annot_categories;
	
	save_html($session_id, $hRes);
	my $hRes2;
	$hRes2->{session_id} = $session_id;
	printJson($hRes2);
}

sub save_html {
	my ($session_id, $hRes) = @_;
	my $session = new session_export();
	$session->load_session( $session_id );
	$session->save( 'html_page_title', $hRes->{'html_page_title'} );
	$session->save( 'html_source', $hRes->{'html_source'} );
	$session->save( 'html_title', $hRes->{'html_title'} );
	$session->save( 'html_gene', $hRes->{'html_gene'} );
	$session->save( 'html_variants', $hRes->{'html_variants'} );
	$session->save( 'hash_filters', $hRes->{'hash_filters'} );
	$session->save( 'gencode_version', $hRes->{'gencode_version'} );
}


sub get_variants_infos_from_projects {
	my ($hResVariants_loaded, $hVariantsDetails, $hVariantsIds, $hResVariantsModels, $use_locus, $only_transcript) = @_;
	my ($h_count, $hVariantsProj);
	print '_update_dv_';
	($h_count, $hVariantsDetails, $hVariantsIds) = update_list_variants_from_dejavu($project_init_name, $gene_init_id_for_newgene, $h_proj_pat_all, $h_proj_pat_ill, $hResVariants_loaded, $hVariantsDetails, $hResVariantsRatioAll, $hResVariantsModels, $use_locus, $only_transcript);
	
	
	my $hgene;
	$hgene->{id} = $gene_init_id;
	my @l_p = keys %$hVariantsProj;
	my $hResVariants;
	
	
	my $fork = 5;
	my $pm = new Parallel::ForkManager($fork);
	
	my $nb_errors=0;
	$pm->run_on_finish(
		sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data) = @_;
			if (not $data or not $data->{start_job} or not exists $data->{start_job}) { $nb_errors++; }
			else {
				# hResVariants   $hRes->{variants}
				# $hVariantsDetails  $hRes->{variants_details}
				# $hResVariantsListPatients $hRes->{variants_list_patients}
				if (exists $data->{variants} and $data->{variants_details} and $data->{variants_list_patients}) {
					foreach my $start (keys %{$data->{variants}}) {
						foreach my $var_id (keys %{$data->{variants}->{$start}}) {
							$hResVariants->{$start}->{$var_id} = $data->{variants}->{$start}->{$var_id};
							$hVariantsDetails->{$var_id} = $data->{variants_details}->{$var_id};
							foreach my $proj_name (keys %{$data->{variants_list_patients}->{$var_id}}){
								foreach my $pat_name (keys %{$data->{variants_list_patients}->{$var_id}->{$proj_name}}) {
									$hResVariantsListPatients->{$var_id}->{$proj_name}->{$pat_name} = $data->{variants_list_patients}->{$var_id}->{$proj_name}->{$pat_name};
								}
							}
						}
					}
					foreach my $gene_id (keys %{$data->{gene}}) {
						$hResGene->{$gene_id} = $data->{gene}->{$gene_id};
					}
				}
			}
		}
	);
	
	
	print 'nbVarFiltred:'.scalar keys %{$hVariantsIds};
	
	
	foreach my $var_id (keys %{$hVariantsIds}) {
		next if (not exists $hVariantsDetails->{$var_id}->{dejavu_details});
 	 	my $pid = $pm->start and next;
 	 	my $hres;
 	 	$hres->{start_job} = 1;
		foreach my $project_name (keys %{$hVariantsIds->{$var_id}}) {
#			warn "\n\n\n";
#			warn $var_id.' -> '.$project_name;
			print ".";
			
			my $h_project_res;
			my $buffer_fork = new GBuffer;
			$buffer_fork->{gencode} = $buffer_init->gencode();
			my $project_fork = $buffer_fork->newProject( -name => $project_name );
#			my $is_cache_rocks = 1;
#			if (not $project_fork or not $project_fork->isRocks) {
#				$is_cache_rocks = undef;
#				$project_fork = undef;
#				$project_fork = $buffer_fork->newProject( -name => $project_name );
#			}
			$project_fork->{version} = 'HG38';
			$project_fork->{genome_version} = 'HG38';
			delete $project_fork->{annotation_genome_version};
			delete $project_fork->{rocks};
			delete $project_fork->{deja_vu_public_dir};
			$project_fork->changeAnnotationVersion('43.20', 1);
			$buffer_fork->{annotation_genome_version} = 'HG38';
			$buffer_fork->{public_data_version} = '20';
			my $gene_fork = $project_fork->newGene($gene_init_id);
			
			my $buffer_fork_hg19 = new GBuffer;
			my $project_fork_hg19 = $buffer_fork_hg19->newProjectCache( -name => $project_name );
			my $gene_fork_hg19;
			$gene_fork_hg19 = $project_fork_hg19->newGene($gene_init_id);
			
			my ($var_hg19, $var_id_hg19, $var_id_hg38);
			eval {
				my $var = $project_fork->_newVariant($var_id);
				my $var_id = $var->id();
				
				$var_id_hg19 = $hVariantsDetails->{$var_id}->{id_hg19};
				$var_id_hg38 = $var_id;
				my @ltmp = split('_', $var_id_hg19);

				
				my $v_pos = $project_fork_hg19->getChromosome($ltmp[0])->getVectorByPosition(int($ltmp[1]) -10000, int($ltmp[1]) +10000);
				foreach my $v_tmp (@{$project_fork_hg19->getChromosome($var->getChromosome->id())->getListVarObjects($v_pos)}) {
					if ($v_tmp->id eq $var_id_hg19) {
						$var_hg19 = $v_tmp;
					}
				}
				if (not $var_hg19) { $var_hg19 = $project_fork_hg19->_newVariant($var_id_hg19); }
				
				my $hDejaVuGlobal = $project_fork_hg19->getDejaVuInfos($var_id_hg19);
				my @l_pat_dv = split(';', $hDejaVuGlobal->{$project_name}->{patients});
			
				$hres->{variants}->{$var->start()}->{$var_id}->{value}->{id} =  $var->id;
				$hres->{variants}->{$var->start()}->{$var_id}->{html}->{id} =  $var->id;
				$hres->{variants}->{$var->start()}->{$var_id}->{value}->{type} = $var->type;
				$hres->{variants}->{$var->start()}->{$var_id}->{html}->{type} = $var->type;
				
				foreach my $p_id (keys %{$hVariantsDetails->{$var_id}->{dejavu_details}->{$project_name}}) {
					my $p_name =  $l_pat_dv[$p_id];
					$hres->{variants}->{$var->start()}->{$var_id}->{html}->{patients} .= $p_name." ";
				}
				
				$hres->{variants}->{$var->start()}->{$var_id}->{html}->{infos} .= $var->{infos};
				my $vn = $var->vcf_id;
				$vn =~ s/_/-/g;
				$vn =~ s/chr//;
				$hres->{variants}->{$var->start()}->{$var_id}->{value}->{gnomad_id} = $vn;
				$hres->{variants}->{$var->start()}->{$var_id}->{html}->{gnomad_id} = $vn;	
				$hres->{variants}->{$var->start()}->{$var_id}->{value}->{is_cnv} = 0;
				
				my $chr = $var->getChromosome();
				
				my $gene = $gene_init;
				my @lPhen;
				foreach my $pheno (@{$project_fork->getPhenotypes()}) {
					push(@lPhen, $pheno->name());
				}
				my $pheno_name = join(", ",@lPhen);
				
				update_variant_editor::vgnomad($var,$hres->{variants}->{$var->start()}->{$var_id});
				update_variant_editor::vname($var,$hres->{variants}->{$var->start()}->{$var_id});
				update_variant_editor::vspliceAI($var,$hres->{variants}->{$var->start()}->{$var_id});
			
				my $found_one;
				foreach my $p_name (@l_pat_dv) {
					next unless $p_name;
					my ($p, $p_hg19);
					eval {
						$p = $project_fork->getPatient($p_name);
						$p_hg19 = $project_fork_hg19->getPatient($p_name);
					};
					if ($@) {
						warn "\n\n";
						warn 'PolyGeneScout - Patient '.$p_name.' not found in '.$project_fork->name();
						warn "\n\n";
						next;
					}
					
					next if $only_ill and not $p->isIll();
					
					$found_one++;
					
					if ($var_hg19->vector_id) {
						update_variant_editor::vsequencing($var_hg19,$hres->{variants}->{$var->start()}->{$var_id},$p_hg19);
					}
					
					update_variant_editor::vdivers($var_hg19,$hres->{variants}->{$var->start()}->{$var_id});
					update_variant_editor::valamut_igv($var,$hres->{variants}->{$var->start()}->{$var_id},$p);
					update_variant_editor::vvarsome($hres->{variants}->{$var->start()}->{$var_id},$p);
					update_variant_editor::vhgmd($var,$hres->{variants}->{$var->start()}->{$var_id});
					
					my $html_var_name_hg38 = $hres->{variants}->{$var->start()}->{$var_id}->{html}->{var_name};
					my $html_var_name_hg19 = $hres->{variants}->{$var->start()}->{$var_id}->{html}->{var_name};
					$html_var_name_hg38 =~ s/gnomad_r2_1/gnomad_r4/;
					my $pos_hg38 = $var->start;
					my $pos_hg19 = $var_hg19->start;
					$html_var_name_hg19 =~ s/$pos_hg38/$pos_hg19/g;
					
					$hres->{variants_details}->{$var_id}->{table_vname} = "<span><b>HG38: </b></span>".$html_var_name_hg38;
					$hres->{variants_details}->{$var_id}->{table_vname} .= "<br><br><span><b>HG19: </b></span>".$html_var_name_hg19;
					
					unless (exists $hVariantsDetails->{$var_id}->{table_validation}) {
						$hres->{variants_details}->{$var_id}->{table_validation} = update_variant_editor::table_validation($p, $hres->{variants}->{$var->start()}->{$var_id}, $gene_init);
					}
					
					unless (exists $hVariantsDetails->{$var_id}->{table_dejavu}) {
						update_variant_editor::vdejavu($var,$hres->{variants}->{$var->start()}->{$var_id});
						$hres->{variants_details}->{$var_id}->{table_dejavu} = $hres->{variants}->{$var->start()}->{$var_id}->{html}->{deja_vu};
					}
					
					unless (exists $hVariantsDetails->{$var_id}->{table_gnomad}) {
						my $table_gnomad = update_variant_editor::table_gnomad($var);
						$hres->{variants_details}->{$var_id}->{table_gnomad} = $table_gnomad;
					}
									
					unless (exists $hVariantsDetails->{$var_id}->{table_transcript}) {
						$hres->{variants}->{$var->start()}->{$var_id}->{genes}->{$gene_init_id_for_newgene} = update_variant_editor::construct_hash_transcript($var, $cgi, \@header_transcripts, 2, $gene_fork);
						$hres->{variants_details}->{$var_id}->{table_transcript} = update_variant_editor::table_transcripts($hres->{variants}->{$var->start()}->{$var_id}->{genes}->{$gene_init_id_for_newgene}, \@header_transcripts, 1);
					}
					
					$hres->{variants_details}->{$var_id}->{table_varsome} = $hVariantsDetails->{$var_id}->{table_varsome};
					$hres->{variants_details}->{$var_id}->{alamut_link_variant} = $hVariantsDetails->{$var_id}->{alamut_link_variant};
					
					if (-d $p_hg19->NoSqlDepthDir() and -e $p_hg19->NoSqlDepthDir().$p_hg19->name.".depth.lmdb") {
						update_variant_editor::trio($var_hg19,$hres->{variants}->{$var->start()}->{$var_id},$p_hg19);
					}
							
					my $b_igv = $hres->{variants}->{$var->start()}->{$var_id}->{'html'}->{'igv'};
					my $chr_id = $chr->id();
					my $start = $var->start();
					
					my (@bams,@names);
					foreach my $p (@{$p->getFamily->getPatients()}){
						next unless -e $p->getBamFileName;
						push(@bams,$p->bamUrl);
						push(@names,$p->name());
					}
					my $f =  join(";",@bams);#$patient->{obj}->bamUrl;;
					my $f2 =  join(",",@bams);#$patient->{obj}->bamUrl;;
					my $v1 = $hres->{variants}->{$var->start()}->{$var_id}->{ref_allele}."/".$hres->{variants}->{$var->start()}->{$var_id}->{allele};	
					my $gn = $p->project->getVersion();
					my $project_name = $p->project->name;
					my $pnames = join(";",@names);
					my $locus = $chr_id.':'.$start.'-'.$start;
					my $gene_name = $gene->external_name();
					$hres->{variants}->{$var->start()}->{$var_id}->{'html'}->{'igv'} =qq{<button class='igvIcon2' style="width:22px;height:22px;" onclick='launch_web_igv_js("$project_name","$pnames","$f","$locus",)' style="color:black"></button>};
					$hres->{variants}->{$var->start()}->{$var_id}->{'html'}->{'alamut'} = qq{<button class="alamutView3" style="width:22px;height:22px;" onClick="httpGetLoadOnlyListBamInGene('$gene_name','$f2');"></button>};
					$hres->{variants}->{$var->start()}->{$var_id}->{html}->{pheno_name} = $pheno_name;
					
					
					if (not exists $hResVariantsListPatients->{$var_id}->{$project_fork_hg19->name()}->{$p_hg19->name()}) {
						my $perc_allele = $var->getPourcentAllele($p_hg19);
						$hResVariantsRatioAll->{$var_id}->{$project_fork_hg19->id()}->{$p_hg19->name()} = $perc_allele;
						$hres->{variants}->{$var->start()}->{$var_id}->{scaled_score_gene_used} = $gene->id;
						my ($table_trio, $model_found) = get_table_trio_from_object($var_hg19, $p_hg19, $gene_fork_hg19, $hres->{variants}->{$var->start()}->{$var_id});
						
						my $ok_model;
						$ok_model = 1 if (lc($model_found) eq 'solo' and $h_models and exists $h_models->{solo});
						$ok_model = 1 if (lc($model_found) eq "father" and exists $h_models->{any_parent});
						$ok_model = 1 if (lc($model_found) eq "father_c" and exists $h_models->{any_parent});
						$ok_model = 1 if (lc($model_found) eq "mother" and exists $h_models->{any_parent});
						$ok_model = 1 if (lc($model_found) eq "mother_c" and exists $h_models->{any_parent});
						$ok_model = 1 if (lc($model_found) eq "strict_denovo" and $h_models and exists $h_models->{denovo});
						$ok_model = 1 if (lc($model_found) eq "denovo" and $h_models and exists $h_models->{denovo});
						$ok_model = 1 if (lc($model_found) eq "recessive" and $h_models and exists $h_models->{recessive});
						$ok_model = 1 if (lc($model_found) eq "mosaic" and $h_models and exists $h_models->{any_parent});
						$ok_model = 1 if (lc($model_found) eq "uniparental" and $h_models and exists $h_models->{any_parent});
						$ok_model = 1 if (lc($model_found) eq "?");
						
						#warn "\n$ok_model";
						
						next if not $ok_model;
						
						$h_project_res->{hResVariantsModels}->{$var_id}->{$p->getProject->name()}->{$p->name} = $model_found;
						$hres->{variants_list_patients}->{$var_id}->{$project_fork_hg19->name()}->{$p_hg19->name()}->{'values'}->{'phenotypes'} = '';
						if ($project_fork_hg19->getPhenotypes()) {
							$hres->{variants_list_patients}->{$var_id}->{$project_fork_hg19->name()}->{$p_hg19->name()}->{'values'}->{'phenotypes'} = join(', ', sort @{$project_fork_hg19->phenotypes()});
						}
						$hres->{variants_list_patients}->{$var_id}->{$project_fork_hg19->name()}->{$p_hg19->name()}->{'values'}->{'description'} = $project_fork->description();
						$hres->{variants_list_patients}->{$var_id}->{$project_fork_hg19->name()}->{$p_hg19->name()}->{'values'}->{'fam'} = $p_hg19->getFamily->name();
						$hres->{variants_list_patients}->{$var_id}->{$project_fork_hg19->name()}->{$p_hg19->name()}->{'values'}->{'name'} = $p_hg19->name();
						$hres->{variants_list_patients}->{$var_id}->{$project_fork_hg19->name()}->{$p_hg19->name()}->{'values'}->{'sex'} = $p_hg19->sex();
						$hres->{variants_list_patients}->{$var_id}->{$project_fork_hg19->name()}->{$p_hg19->name()}->{'values'}->{'status'} = $p_hg19->status();
						
						my $is_heho = '-';
						eval { $is_heho = 'ho' if $var->isHomozygote($p_hg19); };
						if ($@) { $is_heho = 'pb'; }
						eval { $is_heho = 'he' if $var->isHeterozygote($p_hg19); };
						if ($@) {
							if ($is_heho eq 'ho') { $is_heho = 'ho'; }
							else { $is_heho = 'pb'; } 
						}
						
						$hres->{variants_list_patients}->{$var_id}->{$project_fork_hg19->name()}->{$p_hg19->name()}->{'values'}->{'he_ho'} = $is_heho;
						
						my $parent_child = 'solo';
						if ($p->getFamily->isTrio()) {
							$parent_child = 'mother' if ($p_hg19->isMother());
							$parent_child = 'father' if ($p_hg19->isFather());
							$parent_child = 'child' if ($p_hg19->isChild());
						}
						$hres->{variants_list_patients}->{$var_id}->{$project_fork_hg19->name()}->{$p_hg19->name()}->{'values'}->{'parent_child'} = $parent_child;
						
						my $path_polyweb;
						if (-d $Bin.'/../../polyweb/') { $path_polyweb = $Bin.'/../../polyweb/'; }
						elsif (-d $Bin.'/../../PolyWeb/') { $path_polyweb = $Bin.'/../../PolyWeb/'; }
						if (-d $Bin.'/../../../polyweb/') { $path_polyweb = $Bin.'/../../../polyweb/'; }
						elsif (-d $Bin.'/../../../PolyWeb/') { $path_polyweb = $Bin.'/../../../PolyWeb/'; }
						else {
							warn $Bin;
							warn("\n\nPath polyweb not found. die.\n\n");
						}
						my $html_icon = $p_hg19->return_icon();
						$html_icon =~ s/<img src='//;
						$html_icon =~ s/' style='padding-right:1px'>//;
						my $icon = "$path_polyweb/$html_icon";
						$hres->{variants_list_patients}->{$var_id}->{$project_fork_hg19->name()}->{$p_hg19->name()}->{'values'}->{'sex_status_icon'} = $icon;
						$hres->{variants_list_patients}->{$var_id}->{$project_fork_hg19->name()}->{$p_hg19->name()}->{'values'}->{'percent'} = $perc_allele;
						
#						$hres->{variants_list_patients}->{$var_id}->{$project_fork_hg19->name()}->{$p_hg19->name()}->{'values'}->{'model'} = '';
#						$hres->{variants_list_patients}->{$var_id}->{$project_fork_hg19->name()}->{$p_hg19->name()}->{html} = '';
						
						$hres->{variants_list_patients}->{$var_id}->{$project_fork_hg19->name()}->{$p_hg19->name()}->{'values'}->{'model'} = $model_found;
						$hres->{variants_list_patients}->{$var_id}->{$project_fork_hg19->name()}->{$p_hg19->name()}->{html} = $table_trio;
						$hres->{variants_list_patients}->{$var_id}->{$project_fork_hg19->name()}->{$p_hg19->name()}->{percent_allele} = $perc_allele;
					}
	
					#TODO: here				
	#				if (not $hResVariantsTableLocal->{$var_id} or not $hResVariantsTableLocal->{$var_id}) {
	#					$hResVariantsTableLocal->{$var_id} = update_variant_editor::table_validation($p, $h_var, $hgene);
	#				}
					$hres->{gene}->{$gene_init_id} = $hres->{variants}->{$var->start()}->{$var_id}->{'genes'}->{$gene_init_id}->[0]->{'html'};
					
				}
			
			};
			if($@) {
				warn "\n";
				warn "\n";
				warn $@;
				warn 'ERROR '.$var_id.'->'.$project_fork->name();
#				warn 'HG19 id in rocks: '.$var_id_hg19;
#				warn ' HG19 VAR: '.$var_hg19->id;
				$project_fork_hg19 = undef;
				$buffer_fork_hg19 = undef;
				$project_fork->close_rocks();
				$project_fork = undef;
				$buffer_fork = undef;
		 	 	$hres->{end_job} = 1;
		 	 	$pm->finish(0, $hres);
				next;
			}
			$project_fork_hg19 = undef;
			$buffer_fork_hg19 = undef;
			$project_fork->close_rocks();
			$project_fork = undef;
			$buffer_fork = undef;
		}
 	 	$hres->{end_job} = 1;
# 	 	warn  Dumper sort  keys %$hres;
 	 	$pm->finish(0, $hres);
	}
	sleep(3); 
	$pm->wait_all_children();
#	warn  Dumper sort  keys %$hResVariants;
	return ($h_count, $hResVariants, $hVariantsDetails, $hResVariantsModels);
}
	
	
	


sub get_html_gene {
	my ($gene_init_id, $gene_init, $h_count) = @_;
	my @lAllVar = keys %$hResVariantsIds;
	my $nb_selected = $h_count->{total_pass};
	my $total_v = $h_count->{total};
	$gene_init->{buffer} = $buffer_init;
	$gene_init->{project} = $project_init;
	$hResGene->{$gene_init_id}->{id} = $gene_init_id;
	$hResGene->{$gene_init_id}->{external_name} = "<span style='color:white;'>".$gene_init->external_name()."<span>";
	$hResGene->{$gene_init_id}->{variants} = \@lAllVar;
	my ($pheno,$nb_other_terms) = $gene_init->polyviewer_phentotypes();
	$hResGene->{$gene_init_id}->{phenotypes}->{pheno} = $pheno;
	$hResGene->{$gene_init_id}->{phenotypes}->{nb_other_terms} = $nb_other_terms;
	my $description_gene = $gene_init->description();
	$gene_init->getProject->buffer->dbh_deconnect();
	$gene_init->getProject->buffer->dbh_reconnect();
#	delete $gene_init->getProject->{rocksPartialTranscripts};
#	eval {
		foreach my $panel (@{$gene_init->getPanels()}) {
			$hResGene->{$gene_init_id}->{panels}->{$panel->name()}->{phenotype} = $panel->getPhenotypes()->[0]->name();
		}
#	};
#	if ($@) { $hResGene->{$gene_init_id}->{panels} = undef; }
	my $html_gene = update_variant_editor::panel_gene($hResGene->{$gene_init_id});
	$html_gene =~ s/CNV//;
	my $regexp1 = qq{<span class=\'badge badge-infos badge-xs \' style="color:#00C851"  >[0-9]+ </span>};
	my $regexp2 = qq{<span class=\'badge badge-infos badge-xs \' style="color:#00C851"  > $nb_selected/$total_v </span>};
	$html_gene =~ s/$regexp1/$regexp2/;
	if ($description_gene) {
		$html_gene =~ s/<\/table>//;
	 	$html_gene .= qq{<td><span style='color:white;'><i><center>$description_gene</center></i></span></td></table>};
	}
	return $html_gene;
}


sub update_list_variants_from_dejavu {
	my ($proj_name, $gene_init_id_for_newgene, $h_proj_pat_all, $h_proj_pat_ill, $hResVariants_loaded, $hVariantsDetails, $hResVariantsRatioAll, $hResVariantsModels, $use_locus, $only_transcript) = @_;
	my $time = time;
	my $buffer_dejavu = new GBuffer;
	
	$project_dejavu = $buffer_dejavu->newProject( -name => 'NGS2024_7792');
	#$project_dejavu = $buffer_dejavu->newProject( -name => $buffer_dejavu->get_random_project_name_with_this_annotations_and_genecode());
	
	
	my $gene_dejavu = $project_dejavu->newGene($gene_init_id_for_newgene);
	
	
#	warn "\n\n";
#	warn ref($gene_dejavu);
#	warn $gene_dejavu->id();
#	warn $gene_dejavu->external_name();
#	warn $project_dejavu->name();
	
	
	my $transcript_dejavu;
	foreach my $tr (@{$gene_dejavu->getTranscripts()}) {
		$transcript_dejavu = $tr if ($only_transcript and $tr->id() =~ /$only_transcript/);
		$gene_with_partial_transcrit = 1 if ($tr->is_partial_transcript());
	}
	
	# SELECT paquet de $max_variants variants
	my @lProjectNames;
	my ($total, $total_pass);
	my $hVarErrors;
	
#	warn "\n";
#	warn $gene_dejavu->id();
	my $h_dv_rocks_ids = $gene_dejavu->getChromosome->rocks_dejavu->dejavu_interval($gene_dejavu->start(), $gene_dejavu->end());
	
	my (@lVarIds);
	foreach my $rocks_id (keys %{$h_dv_rocks_ids}) {
		my $var_id = $h_dv_rocks_ids->{$rocks_id}->{hg38};
		push(@lVarIds, $var_id);
		$h_dv_rocks_ids->{$var_id}->{hg19} = $h_dv_rocks_ids->{$rocks_id}->{hg19};
		$h_dv_rocks_ids->{$var_id}->{rocks_id} = $rocks_id;
	}
	
	print '.nbRocksIds:'.scalar(keys %{$h_dv_rocks_ids});
	my $h_projects_filters_he_comp;
	if ($only_pat_with_var) {
		my $var_dv = $project_dejavu->_newVariant($only_pat_with_var);
		
		if (not $only_pat_with_var_he or not $only_pat_with_var_ho) {
			foreach my $proj_name (keys %{$h_projects_filters_he_comp}) {
				delete $h_projects_filters_he_comp->{$proj_name} if ($only_pat_with_var_ho and $h_projects_filters_he_comp->{$proj_name}->{ho} == 0);
				delete $h_projects_filters_he_comp->{$proj_name} if ($only_pat_with_var_he and $h_projects_filters_he_comp->{$proj_name}->{he} == 0);
			}
		}
	}
	
	my @lVar;
	foreach my $var_id (@lVarIds) {
		my $var_id_hg19 = $h_dv_rocks_ids->{$var_id}->{hg19};
		print '.';
		$total++;
		my $is_ok_perc = 1;
		if ($filter_perc_allelic_max and $hResVariantsRatioAll and exists $hResVariantsRatioAll->{$var_id}) {
			$is_ok_perc = 0;
			foreach my $project_id (keys %{$hResVariantsRatioAll->{$var_id}}) {
				next if ($is_ok_perc);
				foreach my $patient_id (keys %{$hResVariantsRatioAll->{$var_id}->{$project_id}}) {
					next if ($is_ok_perc);
					$is_ok_perc = 1 if ($hResVariantsRatioAll->{$var_id}->{$project_id}->{$patient_id} < $filter_perc_allelic_max);
				}
			}
		}
		next unless ($is_ok_perc);
		
		my ($var_gnomad, $var_gnomad_ho, $var_annot, $var_dejavu, $var_dejavu_ho, $var_model);
		
		next if $var_id =~ /ALU/;
		next if not $var_id =~ /[XYMT0-9]+_[0-9]+_[ATGC]+_[ATGC]+/;
		my $var = $project_dejavu->_newVariant($var_id);
		if ($hVariantsDetails and exists $hVariantsDetails->{$var_id}) {
			$var_gnomad = $hVariantsDetails->{$var_id}->{var_gnomad} if ($hVariantsDetails->{$var_id}->{var_gnomad});
			$var_gnomad_ho = $hVariantsDetails->{$var_id}->{var_gnomad_ho} if ($hVariantsDetails->{$var_id}->{var_gnomad_ho});
			$var_dejavu = $hVariantsDetails->{$var_id}->{var_dejavu} if ($hVariantsDetails->{$var_id}->{var_dejavu});
			$var_dejavu_ho = $hVariantsDetails->{$var_id}->{var_dejavu_ho} if ($hVariantsDetails->{$var_id}->{var_dejavu_ho});
			if ($only_transcript) {
				$var_annot = $var->variationTypeInterface($transcript_dejavu);
			}
			else {
				$var_annot = $hVariantsDetails->{$var_id}->{annotation} if ($hVariantsDetails->{$var_id}->{annotation});
			}
			$var_model = $hVariantsDetails->{$var_id}->{model} if ($hVariantsDetails->{$var_id}->{model});
		}
		next if ($var->isCnv() or $var->isLarge());
		
#		$var->{rocksdb_id} = $h_dv_rocks_ids->{$var_id}->{rocks_id};
#		warn Dumper $gene_dejavu->getChromosome->getDejaVuInfosForDiagforVariant($var);
#		
#		$hVariantsDetails->{$var_id}->{dejavu}->{other_projects} = $var->other_projects;
#		$hVariantsDetails->{$var_id}->{dejavu}->{other_patients} = $var->other_patients;
#		$hVariantsDetails->{$var_id}->{dejavu}->{other_patients_ho} = $var->other_patients_ho;
		
		my $not_ok;
		
		unless ($var_gnomad) {
			$var_gnomad = $var->getGnomadAC();
		}
		
		$not_ok++ if ($var_gnomad and $max_gnomad and $var_gnomad > $max_gnomad);
		next if $not_ok;
		
		unless ($var_gnomad_ho) {
			$var_gnomad_ho = $var->getGnomadHO();
		}
		$not_ok++ if ($max_gnomad_ho and $var_gnomad_ho > $max_gnomad_ho);
		next if $not_ok;
		
		
		my $is_ok_annot;
		unless ($var_annot) {
			$var->annotation();
			eval {
				if ($transcript_dejavu) { $var_annot = $var->variationTypeInterface($transcript_dejavu); }
				else { $var_annot = $var->variationTypeInterface($gene_dejavu); }
			};
			if ($@) { $var_annot = 'error'; }
		}
		
		foreach my $this_annot (split(',', $var_annot)) {
			$this_annot =~ s/ /_/g;
			$is_ok_annot ++ if (exists $h_filters_cons->{lc($this_annot)});
		}
		next unless ($is_ok_annot);

		unless ($var_dejavu) {
			$var_dejavu = $var->other_patients();
			$hVariantsDetails->{$var_id}->{var_dejavu} = $var_dejavu;
		}
		$not_ok++ if ($max_dejavu and $var_dejavu > $max_dejavu);
		next if $not_ok;
		
		unless ($var_dejavu_ho) {
			$var_dejavu_ho = $var->other_patients_ho();
			$hVariantsDetails->{$var_id}->{var_dejavu_ho} = $var_dejavu_ho;
		}
		$not_ok++ if ($max_dejavu_ho and $var_dejavu_ho > $max_dejavu_ho);
		next if $not_ok;
		
		if ($hResVariants_loaded and exists $hResVariants_loaded->{$var->start()}->{$var_id}->{table_validation}) {
			$hResVariants->{$var->start()}->{$var_id} = $hResVariants_loaded->{$var->start()}->{$var_id};
			$total_pass++;
			next;
		}
		
		
		next unless ($has_proj);
#		next unless ($ok_model);
		
		$hVariantsDetails->{$var_id}->{dejavu_details} = $hashVarId;
		
		
#		my $table_dejavu = update_variant_editor::table_dejavu($var, 'no_phenotype');
		my $table_dejavu;
		my $h_var;
		$h_var->{html}->{done_here} = 1;
		$h_var->{html}->{no_css_polydiag} = 1;
		$h_var->{value}->{id} =  $var->id;
		$h_var->{html}->{id} =  $var->id;
		$h_var->{value}->{type} = $var->type;
		$h_var->{html}->{type} = $var->type;
		
		my $vn = $var->id();
		$vn =~ s/_/-/g;
		$vn =~ s/chr//;
#		$h_var->{value}->{gnomad_id} = $vn;
#		$h_var->{html}->{gnomad_id} = $vn;
		update_variant_editor::vname2($var,$h_var);
		update_variant_editor::vhgmd($var,$h_var);
#		update_variant_editor::vclinvar($var,$h_var);
		update_variant_editor::vvarsome($h_var);
		
		$h_var->{'html'}->{'no_css_polydiag'} = 1;
		my $val1 = 'onClick="zoomHgmd';
		my $Val2 = 'onClick="zoomHgmdWithoutCss';
		$h_var->{'html'}->{'hgmd'} =~ s/$val1/$Val2/;
		#$h_var->{'html'}->{'var_name'} .= "<br><br><span style='font-size:8px;'><b>".$var->getChromosome->id().':'.$var->start().'-'.$var->end()."</span></b>";
		my $table_vname = $h_var->{'html'}->{'var_name'};
		
#		my $cmd_search_he_comp = qq{search_he_comp_with_variant('$var_id', '$user_name');};
#		my $b_search_he_comp = qq{<button onClick="$cmd_search_he_comp" style="margin-top:7px;font-size:9px;">Search He Composite</button>};
#		$table_vname .= qq{<center>$b_search_he_comp</center>};
		
		my $table_varsome = $h_var->{html}->{varsome};
		
		if ($var->hgmd) {
			$h_var->{value}->{hgmd_id} = $var->hgmd_id;
			my $n1 = $project_init_name;
		 	my $n2 = $var->hgmd_id;
		 	my $n3 = $var->id;
			my $cmd_hgmd = qq{zoomHgmdWithoutCss(\'$n1\',\'$n2\',\'$n3\')}; 
			$h_var->{html}->{hgmd} = update_variant_editor::printButton(4,[3,4], $var->hgmd->{class},qq{onClick="$cmd_hgmd"});
			if ($var->isDM()) {
				$h_var->{value}->{dm} = 1;
				#TODO: here
#				if ($var->getChromosome->is_hgmd_DM_for_gene($var->hgmd_id, $gene_dejavu)) {
#					$h_var->{value}->{dm_for_this_gene} = 1;
#				}
#				else {
					$h_var->{value}->{dm_for_this_gene} = undef;
					$h_var->{value}->{dm} = undef;
					$h_var->{value}->{hgmd} = '';
					$h_var->{html}->{hgmd} = '';
#				}
			}
		}
		else { $h_var->{value}->{dm} = ''; }
		
		if ($h_var->{value}->{clinvar_pathogenic}) {
			my $clinvar_id = $h_var->{value}->{clinvar_id};
			if ($var->getChromosome->is_clinvar_pathogenic_for_gene($clinvar_id, $gene_dejavu)) {
				$h_var->{value}->{clinvar_pathogenic_for_this_gene} = 1;
			}
			else {
				$h_var->{value}->{clinvar_pathogenic_for_this_gene} = undef;
				$h_var->{value}->{clinvar_pathogenic} = undef;
				$h_var->{value}->{clinvar} = '';
				$h_var->{html}->{clinvar} = '';
			}
		}
		
#		my $color = '#555';
#		if ($h_var->{value}->{dm} or $h_var->{value}->{clinvar_pathogenic}) { $color = "red"; }
#		elsif  ($h_var->{value}->{hgmd_id} or $h_var->{value}->{clinvar_id}) { $color = "orange"; }
		

		if ($only_transcript) {
			my @new_list;
			foreach my $htr (@{$h_var->{genes}->{$gene_init_id_for_newgene}}) {
				if ($htr->{value}->{trid} eq $only_transcript) {
					push(@new_list, $htr);
				}
			}
			$h_var->{genes}->{$gene_init_id_for_newgene} = \@new_list;
		}
		$hVariantsDetails->{$var_id}->{alamut_link_variant} = html_polygenescout::print_alamut_variant_button($var->alamut_id());
		$hVariantsDetails->{$var_id}->{table_varsome} = $table_varsome;
#		$hVariantsDetails->{$var_id}->{var_gnomad} = $var_gnomad;
#		$hVariantsDetails->{$var_id}->{var_gnomad_ho} = $var_gnomad_ho;
#		$hVariantsDetails->{$var_id}->{var_dejavu} = $var_dejavu;
#		$hVariantsDetails->{$var_id}->{var_dejavu_ho} = $var_dejavu_ho;
#		$hVariantsDetails->{$var_id}->{annotation} = $var_annot;
		$hVariantsDetails->{$var_id}->{id_hg19} = $h_dv_rocks_ids->{$var_id}->{hg19};
	}
	
	$total_pass = scalar keys %$hVariantsDetails;
	
	$buffer_dejavu = undef;
	print 'nbVar:'.$total;
	print 'nbVarPass:'.$total_pass;
	print 'nbProj:'.scalar(@lProjectNames);
	print '@time:'.abs(time) - $time;
	my $h_count;
	$h_count->{total} = $total;
	$h_count->{total_pass} = $total_pass;
	
#	delete $project_dejavu->{rocksPartialTranscripts};
	$project_dejavu->buffer->dbh_deconnect();
	
	
	
	return ($h_count, $hVariantsDetails, $hVariantsIdsDejavu);
}

sub get_table_trio_from_object {
	my ($var, $patient, $gene_used, $h_var, $no_header_project_pat, $no_header_project) = @_;
	my $gene_name =$gene_used->external_name();
	my $fam = $patient->getFamily();
	my $is_solo_trio = 'SOLO';
	$is_solo_trio = 'TRIO' if $fam->isTrio();
	my $project_name = $patient->getProject->name();
	
	my $description = $patient->getProject->description();
	my @l_users = @{$patient->getProject->get_list_emails()};
	my $patient_name = $patient->name();
	my $pheno = undef;
	$pheno = $h_var->{html}->{pheno_name} if ($h_var and exists $h_var->{html}->{pheno_name});
	
	my $color_local = 'black';
	my $color = "#c1c1c1";
	my $model;
	my $model2 = "?";
	my $patient_heho = "-";
	
	if ($var->vector_id()) {
		$patient_heho = "ho" if $var->isHomozygote($patient);
		$patient_heho = "he" if $var->isHeterozygote($patient);
		my $isMotherTransmission = $var->isMotherTransmission($fam,$patient);
		my $isFatherTransmission = $var->isFatherTransmission($fam,$patient);
		my $mode;
		eval { $model = $var->getTransmissionModel($fam,$patient,$gene_used); };
		if ($@) { $model = '?'; }
		$model2 = qq{<i class="fa fa-male  fa-2x" style="color:lightgrey"></i><i class="fa fa-female  fa-2x" style="color:lightgrey"></i>};
		$h_var->{model} = lc($model);
		if ($is_solo_trio eq 'SOLO') {
			return ('model_filter', lc($model)) if ($h_models and not exists $h_models->{solo});
			$model2 = "<center><span style='font-color:black'>SOLO</span></center>";
			$color_local = "black";
		}
		elsif (lc($model) eq "father" or (lc($model) eq 'father_c' and $isFatherTransmission)) {
			return ('model_filter', lc($model)) if (not exists $h_models->{any_parent});
			my $text_compound = ' ';
			if ($patient->getFamily->getFather && $patient->getFamily->getFather->isIll() ) {
				$model2 = "<center><img src='/icons/Polyicons/male-d.png'>$text_compound</center>";;
			}
			else {
				$model2 = "<center><img src='/icons/Polyicons/male-s.png'>$text_compound</center>";
			}
		}
		elsif (lc($model) eq "mother" or (lc($model) eq 'mother_c')) {
			return ('model_filter', lc($model)) if ($h_models and not exists $h_models->{any_parent});
			my $text_compound = ' ';
			#$text_compound = ' mother_c' if(lc($model) eq 'mother_c');
			if ($patient->getFamily->getMother && $patient->getFamily->getMother->isIll()) {
				$model2 = "<center><img src='/icons/Polyicons/female-d.png'>$text_compound</center>";
			}
			else { $model2 = "<center><img src='/icons/Polyicons/female-s.png'>$text_compound</center>"; }
		}
		elsif (lc($model) =~ "strict") {
			return ('model_filter', lc($model)) if ($h_models and not exists $h_models->{strict_denovo});
			$model2 = qq{Strict Denovo};
			$model = 'strict_denovo';
			$color_local = "white";
		}
		elsif (lc($model) =~ "denovo") {
			return ('model_filter', lc($model)) if ($h_models and not exists $h_models->{denovo});
			$model2 = qq{Denovo};
			$model = 'denovo';
			$color_local = "white";
		}
		elsif (lc($model) eq "#E74C3C") {
			return ('model_filter', lc($model)) if ($h_models and not exists $h_models->{recessive});
			$model2 = qq{Recessive};
			$color_local = "white";
		}
		elsif (lc($model) eq "mosaic") {
			return ('model_filter', lc($model)) if ($h_models and not exists $h_models->{any_parent});
			$model2 = qq{mosaic }.$var->isMosaicTransmission($fam,$patient);
			$color_local = "white";
		}
		elsif (lc($model) =~ /uniparental/) {
			return ('model_filter', lc($model)) if ($h_models and not exists $h_models->{any_parent});
			$model2 = qq{Uniparental}  if (lc($model) =~ /uniparental/);
		}
		else {
			return 'model_filter' if ($h_models and not exists $h_models->{both});
			$color_local = "white";
		}
		
		if ($is_solo_trio eq 'SOLO') { $color = "white"; }
		else { $color = update_variant_editor::color_model(lc($model)); }
	}
	
	my $patient_status;
	if ($patient->isIll()) {
		if ($patient->isMother()) { $patient_status = "<center><img style='width:20px;height:20px;' src='/icons/Polyicons/female-d.png'></center>"; }
		elsif ($patient->isFather()) { $patient_status = "<center><img style='width:20px;height:20px;' src='/icons/Polyicons/male-d.png'></center>"; }
		else {
			if ($patient->sex() eq '1') { $patient_status = "<center><img style='width:20px;height:20px;' src='/icons/Polyicons/baby-boy-d.png'></center>"; }
			else { $patient_status = "<center><img style='width:20px;height:20px;' src='/icons/Polyicons/baby-girl-d.png'></center>"; }
		}
	}
	else {
		if ($patient->isMother()) { $patient_status = "<center><img style='width:20px;height:20px;' src='/icons/Polyicons/female-s.png'></center>"; }
		elsif ($patient->isFather()) { $patient_status = "<center><img style='width:20px;height:20px;' src='/icons/Polyicons/male-s.png'></center>"; }
		else {
			if ($patient->sex() eq '1') { $patient_status = "<center><img style='width:20px;height:20px;' src='/icons/Polyicons/baby-boy-s.png'></center>"; }
			else { $patient_status = "<center><img style='width:20px;height:20px;' src='/icons/Polyicons/baby-girl-s.png'></center>"; }
		}
	}
		
	my $nb_col_span = 6;
	$nb_col_span = 7 if ($var->getProject->isGenome() && $var->isCnv);
	my $hstatus = $patient->validation_status();
	my $hval = $patient->validations();
	
	my $table_validation_id = 'tr_val_'.$var->id().'_'.$patient->getProject->name().'_'.$patient->name();
	my ($local_text, $local_text_tab, $var_id_validated, $gene_name_validated);
	if (exists $hstatus->{$patient->id()}) {
		my @lRes = @{$hstatus->{$patient->id()}}; 
		$hstatus = $lRes[0];
		my $term = $hstatus->{term};
		if ($term eq 'uncertain significance') { $term = 'uncertain sign.'; }
		if ($term eq 'likely pathogenic') { $term = 'likely path.'; }
		my $user_name = $hstatus->{user_name};
		my ($date, $hour) = split(' ', $hstatus->{date});
		($local_text_tab,$var_id_validated,$gene_name_validated) = validation_table_new($patient, $hval);
		my $style_flou;
		if ($term eq 'negative') {
			$local_text = qq{<center><div style='text-align:center;width:90%;padding-bottom:2px;padding-top:2px;color:$color_local;$style_flou'><div style='background-color:$color;border:solid 0.5px $color_local;width:95%;'><b>$term</b></div></div></center>};
		}
		elsif ($var_id_validated eq $var->id()) {
			if ($gene_name_validated eq $gene_name) {
				$style_flou = "opacity: 0.35;";
				$local_text = qq{<center><div style='text-align:center;width:90%;padding-bottom:2px;padding-top:2px;color:$color_local;$style_flou'><button onclick="document.getElementById('$table_validation_id').style.display='block';" style='background-color:$color;border:solid 0.5px $color_local;width:95%;'><b><u>$term</u></b><br>SAME Gene<br>SAME variant</button></div></center>};
			}
			else {
				$local_text = qq{<center><div style='text-align:center;width:90%;padding-bottom:2px;padding-top:2px;color:$color_local;$style_flou'><button onclick="document.getElementById('$table_validation_id').style.display='block';" style='background-color:$color;border:solid 0.5px $color_local;width:95%;'><b><u>$term</u></b><br>OTHER Gene<br>SAME variant</button></div></center>};
			}
		}
		else { 
			$local_text = qq{<center><div style='text-align:center;width:90%;padding-bottom:2px;padding-top:2px;color:$color_local;$style_flou'><button onclick="document.getElementById('$table_validation_id').style.display='block';" style='background-color:$color;border:solid 0.5px $color_local;width:95%;'><b><u>$term</u></b><br>OTHER gene<br>OTHER variant</button></div></center>};
		}
	}
	my $table_trio = qq{ <div> };
	$table_trio .= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:3px"});
	unless ($no_header_project) {
		$table_trio .= $cgi->start_Tr();
		my $users = join("<br>", @l_users);
		my $proj_text = qq{<button onclick="get_popup_users('$users');">Users</button> - <b>$project_name</b>};
		$proj_text .= qq{<sup>defidiag</sup>} if $patient->project->isDefidiag; 
		$proj_text .= " - <span style='color:red;'>$pheno</span>" if ($pheno);
		$proj_text .= "<br>$description";
		my $igv_alamut_text = $h_var->{html}->{igv}."&nbsp;".$h_var->{html}->{alamut};
		my $b_others_var;
		my $pat_name = $patient->name();
		my $var_id = $var->id();
		my $father_trans = undef;
		my $mother_trans = undef;
		my $other_trans = undef;
		
		if ($patient->isChild() and $patient->getFamily->isTrio()) {
			my ($father_name, $mother_name);
			eval { $father_name = $patient->getFamily->getFather->name(); };
			if ($@) {}
			else {
				my $cmd_f = qq{view_var_from_proj_gene_pat('$project_name','$gene_init_id','$pat_name','$var_id', 'father');};
				$father_trans = qq{<button style="text-align:middle;vertical-align:bottom;width:20px;border:solid 0.5px black;background-color:white;" onClick="$cmd_f"><img style='width:8px;height:8px;' src='/icons/Polyicons/male.png'></button>};
			}
			eval { $mother_name = $patient->getFamily->getMother->name(); };
			if ($@) {}
			else {
				my $cmd_m = qq{view_var_from_proj_gene_pat('$project_name','$gene_init_id','$pat_name','$var_id', 'mother');};
				$mother_trans = qq{<button style="text-align:middle;vertical-align:top;width:20px;border:solid 0.5px black;background-color:white;" onClick="$cmd_m"><img style='width:8px;height:8px;' src='/icons/Polyicons/female.png'></button>};
			}
		}
		if ($patient->isChild() and $patient->getFamily->isTrio()) {
			my ($father_name, $mother_name);
			if ($patient->getFamily->father()) {
				$father_name = $patient->getFamily->getFather->name();
				my $cmd_f = qq{view_var_from_proj_gene_pat('$project_name','$gene_init_id','$pat_name','$var_id', 'father');};
				$father_trans = qq{<button style="text-align:middle;vertical-align:bottom;width:20px;border:solid 0.5px black;background-color:white;" onClick="$cmd_f"><img style='width:8px;height:8px;' src='/icons/Polyicons/male.png'></button>};
			}
			if ($patient->getFamily->mother()) {
				$mother_name = $patient->getFamily->getMother->name();
				my $cmd_m = qq{view_var_from_proj_gene_pat('$project_name','$gene_init_id','$pat_name','$var_id', 'mother');};
				$mother_trans = qq{<button style="text-align:middle;vertical-align:top;width:20px;border:solid 0.5px black;background-color:white;" onClick="$cmd_m"><img style='width:8px;height:8px;' src='/icons/Polyicons/female.png'></button>};
			}
		}
		else {
			my $img_child = qq{<img style='width:14px;height:14px;' src='/icons/Polyicons/baby-boy.png'>};
			$img_child = qq{<img style='width:14px;height:14px;' src='/icons/Polyicons/baby-girl.png'>} if ($patient->sex() == 2);
			my $cmd_others = qq{view_var_from_proj_gene_pat('$project_name','$gene_init_id','$pat_name','$var_id');};
			$other_trans = qq{<button style="text-align:middle;vertical-align:top;width:28px;border:solid 0.5px black;background-color:white;" onClick="$cmd_others">$img_child</button>};
		}
		$table_trio .= $cgi->td({colspan=>($nb_col_span-1)}, $proj_text);
		
		my $child_trans_strict_denovo = undef;
		
		#if ($is_solo_trio eq 'SOLO') {
		if ($model2 eq 'Strict Denovo' or $model2 eq 'Denovo' or $is_solo_trio eq 'SOLO') {
			my $img_child = qq{<img style='width:14px;height:14px;' src='/icons/Polyicons/baby-boy.png'>};
			$img_child = qq{<img style='width:14px;height:14px;' src='/icons/Polyicons/baby-girl.png'>} if ($patient->sex() == 2);
			my $cmd_others = qq{view_var_from_proj_gene_pat('$project_name','$gene_init_id','$pat_name','$var_id');};
			$child_trans_strict_denovo = qq{<td rowspan=2><button style="text-align:middle;vertical-align:top;width:28px;border:solid 0.5px black;background-color:white;" onClick="$cmd_others">$img_child</button></td>};
			$table_trio .= $cgi->td({colspan=>1, style=>"text-align:center;vertical-align:middle;width:12px;"}, qq{<table style='width:100%;'><tr><td>$igv_alamut_text &nbsp;</td>$child_trans_strict_denovo</tr></table>});
		}
		elsif ($father_trans and $mother_trans and not $no_header_project_pat) {
			my $local_table = "<table style='width:100%;'>";
			$local_table .= qq{<tr><td rowspan=2>$igv_alamut_text &nbsp;</td><td>$father_trans</td>$child_trans_strict_denovo</tr>};
			$local_table .= qq{<tr><td>$mother_trans</td></tr>};
			$local_table .= qq{</table>};
			$table_trio .= $cgi->td({colspan=>1, style=>"text-align:center;vertical-align:middle;width:12px;"}, $local_table);
		}
		elsif ($father_trans and not $no_header_project_pat) {
			my $local_table = "<table style='width:100%;'>";
			$local_table .= qq{<tr><td rowspan=2>$igv_alamut_text &nbsp;</td><td>$father_trans</td>$child_trans_strict_denovo</tr>};
			$local_table .= qq{<tr><td><button style="text-align:middle;vertical-align:top;width:20px;border:solid 0.5px black;background-color:white;opacity:0.3;" disabled><img style='width:8px;height:8px;' src='/icons/Polyicons/baby-girl.png'></button></td></tr>};
			$local_table .= qq{</table>};
			$table_trio .= $cgi->td({colspan=>1, style=>"text-align:center;vertical-align:middle;width:12px;"}, $local_table);
		}
		elsif ($mother_trans and not $no_header_project_pat) {
			my $local_table = "<table style='width:100%;'>";
			$local_table .= qq{<tr><td rowspan=2>$igv_alamut_text &nbsp;</td><td>$mother_trans</td>$child_trans_strict_denovo</tr>};
			$local_table .= qq{<tr><td><button style="text-align:middle;vertical-align:top;width:20px;border:solid 0.5px black;background-color:white;opacity:0.3;" disabled><img style='width:8px;height:8px;' src='/icons/Polyicons/baby-boy.png'></button></td></tr>};
			$local_table .= qq{</table>};
			$table_trio .= $cgi->td({colspan=>1, style=>"text-align:center;vertical-align:middle;width:12px;"}, $local_table);
		}
		elsif ($other_trans and not $no_header_project_pat) {
			my $local_table = "<table style='width:100%;'>";
			$local_table .= qq{<tr><td>$igv_alamut_text &nbsp;</td><td>$other_trans</td>$child_trans_strict_denovo</tr>};
			$local_table .= qq{</table>};
			$table_trio .= $cgi->td({colspan=>1, style=>"text-align:center;vertical-align:middle;width:12px;"}, $local_table);
		}
#		else {
#			$table_trio .= $cgi->td({colspan=>1, style=>"text-align:center;vertical-align:middle;width:12px;"}, qq{<table style='width:100%;'><tr><td>$igv_alamut_text &nbsp;</td>$child_trans_strict_denovo</tr></table>});
#		}
		$table_trio .= $cgi->end_Tr();
		if ($no_header_project_pat) {
			$table_trio .= '</table>';
			$table_trio .= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:3px"});
		}
	}
	$table_trio .= $cgi->start_Tr({style=>"background-color:$color;"});
	$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;background-color:$color;"}, "<span style='text-align:left;'>$patient_name $local_text</span>");
	$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;background-color:$color;"}, "$patient_status");
	$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;background-color:$color;"}, "$patient_heho");
	if ($var->getProject->isGenome() && $var->isCnv) {
		$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;background-color:$color;"}, 'pr:'.$var->pr($patient));
		$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;background-color:$color;"}, 'sr:'.$var->sr($patient));
		my $cnv_score = sprintf("%.2f", log2($patient->cnv_value_dude($var->getChromosome->name,$var->start,$var->start+$var->length)));
		$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;background-color:$color;"}, 'cnv_score:'.$cnv_score);
	}
	else {
		my $perc_allele = $var->getPourcentAllele($patient);
		return 'perc_all_filter' if ($filter_perc_allelic_min and $perc_allele < $filter_perc_allelic_min);
		return 'perc_all_filter' if ($filter_perc_allelic_max and $perc_allele > $filter_perc_allelic_max);
		$perc_allele .= "%" if ($perc_allele ne '-');
		$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;background-color:$color;"}, $perc_allele);
		my $dp;
		eval { $dp = 'DP: '.$var->getDepth($patient); };
		if ($@) { $dp = '?'; }
		$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;background-color:$color;"}, $dp);
	}
	
	$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;background-color:$color;"}, $model2);
	$table_trio .= $cgi->end_Tr();
	$table_trio .= "</table></div>";
	if ($local_text_tab) {
		$table_trio .=qq{ <div onClick="document.getElementById('$table_validation_id').style.display='none';" id='$table_validation_id' style='display:none;'>$local_text_tab</div> };
	}
	return ($table_trio, lc($model));
}

sub get_hash_annot_categories {
	my $h_annot_categories;
	foreach my $cat_name (keys %{$buffer_init->config->{ensembl_annotations}}) {
		$h_annot_categories->{lc($cat_name)} = $cat_name;
		my @lOthersNames = split(';', $buffer_init->config->{ensembl_annotations}->{$cat_name});
		foreach my $other_name (@lOthersNames) {
			$other_name =~ s/ /_/g;
			$h_annot_categories->{lc($other_name)} = $cat_name;
		}
	}
	return $h_annot_categories;
}

sub get_hash_users_projects {
	my ($user_name, $pwd) = @_;
	my $h_projects;
	my @list_hash = @{$buffer_init->getQuery()->getProjectListForUser($user_name, $pwd)};
	foreach my $hash (@list_hash) {
		my $proj_name = $hash->{name};
		next unless ($proj_name =~ /NGS20/);
		$h_projects->{$proj_name}->{description} = $hash->{description};
		$h_projects->{$proj_name}->{id} = $hash->{id};
	}
	return $h_projects;
}

sub printJson {
	my ($hashRes, $test) = @_;
	my $json_encode = encode_json $hashRes;
	print ".\",";
	$json_encode =~ s/{//;
	print $json_encode;
}

sub validation_table_new {
	my ($p,$hval) = @_;
	my @headers_validations = ("var_name","trio","gene","table_transcript",);
	my @header_transcripts = ("consequence","enst","nm","nomenclature");
	my $rowspan = scalar(keys %$hval);
	my $class;
 	$class->{rowspan} -= $rowspan;
 	$class->{rowspan} = 1 if $class->{rowspan} <=0;
 	$class->{style} = "min-width:10%;padding:1px";
 	my $fam = $p->getFamily();
 	my $fin = scalar (keys %{$hval}) ;
 	my $pos =0;
 	my $patient_name = $p->name();
 	my $project_name = $p->getProject->name();
	my $out;
 	$out.= $cgi->start_table({class=>"table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 9px;font-family:  Verdana;margin-bottom: 0px;padding:2px", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
 	my ($var_valid, $gene_valid);
 	foreach my $k  (keys %{$hval}){
 		my $val = $hval->{$k}->[0];
		my $color = "#9BC8A5";
 		$color = "#E74C3C" if $val->{validation}== 5 ;
 		$color = "coral" if $val->{validation}== 4 ;
 		$color = "#217DBB" if $val->{validation}== -3 ;
 		$color = "orange" if $val->{validation}== 3 ;
	 	my $v  = $p->getProject->_newVariant($val->{polyid});
		my $btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: $color;$fsize"} ;
		my $gene_id = $val->{gene_name};
		my $gene_name = $val->{gene_name};
		my $st_date = join("-",return_date($val->{modification_date}));
 		my $term = $val->{term};
 		my $pn = $val->{user_name};
 		$out .= "<tr>";
 		$out .= $cgi->td({style=>"vertical-align:middle;padding:1px;text-align: center;vertical-align:middle", rowspan=>2},qq{ <span  class="stamp2" style="border-color:$color;color:$color">+$term+<br><span style='font-size=08px'>$pn:$st_date</span></span>});
 		$out .= $cgi->td({style=>"vertical-align:middle;padding:5px;"}, $gene_name);
 		$out .= $cgi->td({style=>"vertical-align:middle;padding:5px;"}, $val->{polyid});
	 	$out .= "</tr>";
	 	$out .= "<tr>";
	 	my @lPh = split(';', $val->{phenotypes_name});
 		$out .= $cgi->td({style=>"vertical-align:middle;padding:1px;color:red;", colspan=>2}, join("<br>", @lPh));
		$out .= "</tr>";
 		$var_valid = $val->{polyid} unless ($var_valid);
 		$gene_valid = $gene_name unless ($gene_valid);
 	}
 	$out.= $cgi->end_table();
 	$out =~ s/"/'/g;
 	return ($out, $var_valid, $gene_valid);
}

sub get_hash_patients_ill_from_project_name {
	my ($dbh, $query_init, $project_name) = @_;
	my $sql_prepare = $dbh->prepare($query_init->sql_cmd_get_ill_patients_from_project());
	$sql_prepare->execute($project_name);
	my $h = $sql_prepare->fetchall_hashref("pat_name");
	return $h;
}

sub get_hash_patients_all_from_project_name {
	my ($dbh, $query_init, $project_name) = @_;
	my $sql_prepare = $dbh->prepare($query_init->sql_cmd_get_all_patients_from_project());
	$sql_prepare->execute($project_name);
	my $h = $sql_prepare->fetchall_hashref("pat_name");
	my $hres;
	my $i = 0;
	foreach my $pat_name (keys %{$h}) {
		$i++;
		$hres->{id}->{$i} = $pat_name;
		$hres->{name}->{$pat_name} = $i;
	}
	return $hres;
}

sub return_date {
	my ($dd) = @_;
	my @amonths = ('Jan', 'Feb', 'Mar', 'Apr','May',"Jun","Jul","Aug","Sep","Oct","Nov","Dec");
	my ($date,$time) = split(" ",$dd);
    my ($year,$month,$day) =  split("-",$date);
	return ($year,$amonths[$month-1],$day);
}


