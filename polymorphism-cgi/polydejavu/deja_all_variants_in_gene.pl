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
use MIME::Base64;

require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/html_polygenescout.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";

my $io = IO::Handle->new();
$io->autoflush(1);

my $max_dejavu = 999999999999;
my $max_dejavu_ho = 999999999999;
my $max_gnomad = 999999999999;
my $max_gnomad_ho = 999999999999;
my $fork;

my $cgi = new CGI();
my $user_name = $cgi->param('user');
my $pwd = $cgi->param('pwd');
my $variant_id = $cgi->param('variant');
my $release = $cgi->param('release');
my $gene_id = $cgi->param('gene');
my $gene_id_alt = $cgi->param('gene_alt');
my $only_transcript = $cgi->param('only_transcript');
$max_dejavu = $cgi->param('dejavu');
$max_dejavu_ho = $cgi->param('dejavu_ho');
$max_gnomad = $cgi->param('gnomad');
$max_gnomad_ho = $cgi->param('gnomad_ho');
my $filters_cons = $cgi->param('filters_cons');
my $only_ill = $cgi->param('only_ill');
my $only_strict_ill = $cgi->param('only_strict_ill');
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
my $debug = $cgi->param('debug');
my $is_polybtf = $cgi->param('is_polybtf');

$fork = $cgi->param('fork');

$fork = 6 if not $fork;

$only_ill = 1 if ($only_strict_ill);

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

$user_name = lc($user_name);
my $buffer_init = new GBuffer;
my $can_use_hgmd = $buffer_init->hasHgmdAccess($user_name);

my $dir_parquet = $buffer_init->dejavu_parquet_dir();

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

my $hProjects_not_dejavu;

my $nb_variants_before_filters;
my @headers_validations = ("#", "varsome","alamut variant","var_name","projects / patients","gnomad","deja_vu","table_validation","table_transcript");
#my @header_transcripts = ("consequence","enst","nm","ccds","appris","exon","nomenclature","codons","codons_AA", "polyphen","sift","ncboost","cadd","revel","dbscsnv","spliceAI");
my @header_transcripts = ("consequence","enst","nm","ccds","appris","exon","nomenclature","codons","codons_AA", "polyphen","sift",'alphamissense',"cadd","revel","dbscsnv",'spliceAI');
my ($hVariantsIdsDejavu, $hVariantsDetails);
my ($hResGene, $hResVariants, $hResVariantsIds, $hResVariantsListPatients, $hResProjectProblem, $hResVariants_byScores, $hResVariantsRatioAll, $hResVariantsModels);
my $project_init_name;

#project preload_patients test connection DB
my $dbh_init = $buffer_init->dbh();
my $query_init = $buffer_init->getQuery();
my $project_init;

my @lProjectNames;
my ($h_projects_name_with_capture,$h_proj_pat_all);
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
	foreach my $p_name (reverse sort keys %{$hProjects}) {
		next if not exists $h_projects_name_with_capture->{$p_name};
		push(@lProjectNames, $p_name);
	}
}
else {
	@lProjectNames = reverse sort keys %{$hProjects};
}

my $hProjectsIds;
foreach my $project_name (@lProjectNames) {
	$hProjectsIds->{$hProjects->{$project_name}->{id}} = $project_name;
	my $h;
	$h->{name} = $project_name;
	$h->{description} = $hProjects->{$project_name}->{description};
	$h->{id} = $hProjects->{$project_name}->{id};
	push(@lItemsProjects, $h);
}

$project_init_name = $buffer_init->get_random_project_name_with_this_annotations_and_genecode();

my $project_name_hg38 = $buffer_init->getRandomProjectName('HG38_DRAGEN');
my $project_name_hg19 = $buffer_init->getRandomProjectName('HG19_CNG');
$project_init = $buffer_init->newProject( -name => $project_name_hg38 );

my $genomeFai_init = $project_init->getGenomeFai();
my $only_variant;
my ($gene_init, $gene_ensg, $gene_chr_id, $gene_init_id, $gene_init_id_for_newgene);
if ($gene_id) {
	eval { $gene_init = $project_init->newGene($gene_id); };
	if ($@) {
		$gene_init = $project_init->newGene($gene_id_alt);
		supressCoreFilesFound();
	}
	$gene_init_id_for_newgene = $gene_id;
	eval { ($gene_ensg, $gene_chr_id) = split('_', $gene_init->id()); };
	if ($@) {
		$gene_init_id_for_newgene = $gene_id_alt;
		$gene_init = $project_init->newGene($gene_id_alt);
		($gene_ensg, $gene_chr_id) = split('_', $gene_init->id());
		supressCoreFilesFound();
	}
}
#TODO: here
elsif ($variant_id) {
	if ($release eq 'HG19') {
		my $buffer_hg19 = new GBuffer;
		my $project_hg19 = $buffer_hg19->newProject( -name => $buffer_hg19->getRandomProjectName('HG19_CNG') );
		$variant_id =~ s/-/_/g;
		my $v_hg19 = $project_hg19->_newVariant($variant_id);
		$variant_id = $v_hg19->lift_over('HG38')->{id};
		$project_hg19 = undef;
		$buffer_hg19 = undef;
	}
	$variant_id =~ s/-/_/g;
	$only_variant = $project_init->_newVariant($variant_id);
	my @lGenes = @{$only_variant->getGenes()};
	$gene_init = $lGenes[0];
	$gene_init_id_for_newgene = $gene_init->id();
	($gene_ensg, $gene_chr_id) = split('_', $gene_init->id());
	$max_dejavu = 999999999999;
	$max_dejavu_ho = 999999999999;
	$max_gnomad = 999999999999;
	$max_gnomad_ho = 999999999999;
	$fork = 1; 
}

my $h_count;
if ($is_polybtf) {
	($hResVariants) = get_variants_infos_from_polybtf();
}
else {
	my $chr_init = $project_init->getChromosome($gene_chr_id);
	my $chr_init_id = $chr_init->id();
	$gene_init->getTranscripts();
	$gene_init_id = $gene_init->id();
	
	my ($use_locus, $only_rocks_id);
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
		next if ($gene_init->getChromosome->id() ne $locus_chr);
		$use_locus = $gene_chr_id.'_'.$locus_start.'_'.$locus_end;
	}
	elsif ($only_variant) {
		$only_rocks_id = $only_variant->rocksdb_id();
		$use_locus = $gene_chr_id.'_'.($only_variant->start() - 1).'_'.($only_variant->end() + 1);
	}
	
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
	#($h_count, $hResVariants, $hVariantsDetails, $hResVariantsModels) = get_variants_infos_from_projects($hResVariants_loaded, $hVariantsDetails, $hResVariantsModels, $use_locus, $only_transcript);
	#my $nb_var_after = scalar(keys %$hVariantsDetails);

	($hResVariants) = get_variants_infos_from_projects($use_locus, $only_transcript, $only_rocks_id);
}


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
	my $found_core;
	opendir my $dir, "$Bin/" or die "Cannot open directory: $!";
	my @files = readdir $dir;
	closedir $dir;
	foreach my $file (@files) {
		$found_core = 1 if ($file =~ /core\./);
	}
	return if not $found_core;
	my $cmd = "rm $Bin/core.*";
	`$cmd`;
}

sub save_export_xls {
	print '.saveExport.';
	$project_init->buffer->dbh_deconnect();
	$project_init->buffer->dbh_reconnect();
#	delete $project_dejavu->{rocksPartialTranscripts};
	my $h_xls_args;
	print '_save_xls_';
	$project_init->cgi_object(1);
	my (@lVarObj, $h_pubmed);
#	foreach my $pos (sort keys %{$hResVariants}) {
	foreach my $var_id (sort keys %{$hResVariants}) {
		$project_init->print_dot(50);
		my $v = $project_init->_newVariant($var_id);
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
#	}
	my $xls_export = new xls_export();
	if ($is_polybtf) {
		$xls_export->title_page('Polybtf.xls');
	}
	else {
		$xls_export->title_page('GeneScout_'.$gene_init->external_name().'.xls');
	}
	$project_init->buffer->dbh_deconnect();
	$project_init->buffer->dbh_reconnect();
	
	$xls_export->store_variants_infos(\@lVarObj, $project_init);
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
	$project_init->buffer->dbh_deconnect();
	$project_init->buffer->dbh_reconnect();
#	delete $project_dejavu->{rocksPartialTranscripts};
	print "|";
	$xls_export->store_specific_infos('projects_patients_infos', $h_patients);
	print "|";
	$xls_export->store_specific_infos('variants_pubmed', $h_pubmed);
	print "|";
	my $session_id = $xls_export->save();
	print "|";
	$project_init->buffer->dbh_deconnect();
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
	
	my @lTrLines;
	my $nb_var = 0;
	my $nb_var_filtred = 0;
	
#	foreach my $pos (sort keys %{$hResVariants}) {
	foreach my $var_id (sort keys %{$hResVariants}) {
#		next if not exists $hResVariantsListPatients->{$var_id};
#		next if scalar keys %{$hResVariantsListPatients->{$var_id}} == 0;
		$nb_var++;
					
		my $hvariation->{html} = $hResVariants->{$var_id};
		unless ($can_use_hgmd) {
			$hvariation->{html}->{hgmd} = qq{<span class="glyphicon glyphicon-ban-circle" aria-hidden="true" style='font-size:12px;color:black;'></span>};
		}
		my $class_default;
		$class_default->{style} = "max-width:350px;overflow-x:auto;vertical-align:middle;padding:5px;";
		
#		my @l_pat;
#		foreach my $project_name (keys %{$hResVariantsListPatients->{$var_id}}) {
#			foreach my $pat_name (keys %{$hResVariantsListPatients->{$var_id}->{$project_name}}) {
#				my $html_table = $hResVariantsListPatients->{$var_id}->{$project_name}->{$pat_name}->{html};
#				my $percent_allele = $hResVariantsListPatients->{$var_id}->{$project_name}->{$pat_name}->{percent_allele};
#				if ($filter_perc_allelic_max) {
#					push(@l_pat, $html_table) if ($percent_allele < $filter_perc_allelic_max);
#				}
#				else {
#					push(@l_pat, $html_table);
#				}
#			}	
#		}
#		
#		next unless (@l_pat);
#		my $nb_pat = scalar(@l_pat);
		my $table_trio;
#		if (not $only_project and not $only_patient) {
#			$table_trio .= qq{<div style="text-align:center;font-size:9px;padding-bottom:5px;"><div style="border:solid 0.5px black;"><b><i>Nb Patients: $nb_pat</b></i></div></div>} if ($nb_pat > 1);
#		}
		$table_trio .= qq{<div style="max-height:170px;max-width:400px;overflow-y:auto;">};
		
		if ($hResVariants->{$var_id}->{table_projects_patients}) {
			$table_trio .= $hResVariants->{$var_id}->{table_projects_patients};
		}
		
#		if ($only_project and $only_patient) { $table_trio .= join("", @l_pat); }
#		else { $table_trio .= join("<br>", @l_pat); }
		$table_trio .= qq{</div>};
		
#		$VAR1 = 'alamut_link_variant';
#		$VAR2 = 'annotation';
#		$VAR3 = 'dejavu_details';
#		$VAR4 = 'id_hg19';
#		$VAR5 = 'spliceAI';
#		$VAR6 = 'table_dejavu';
#		$VAR7 = 'table_gnomad';
#		$VAR8 = 'table_projects_patients';
#		$VAR9 = 'table_transcript';
#		$VAR10 = 'table_validation';
#		$VAR11 = 'table_varsome';
#		$VAR12 = 'table_vname';
#		$VAR13 = 'var_dejavu';
#		$VAR14 = 'var_dejavu_ho';
#		$VAR15 = 'var_gnomad';
#		$VAR16 = 'var_gnomad_ho';
		
				
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
			
			if ($h eq "projects / patients" or $h eq "table_projects_patients") {
				$out .= $cgi->td($class, $table_trio);
			}
			elsif ($h eq "gnomad") {
				$out .= $cgi->td($class_default,$hResVariants->{$var_id}->{table_gnomad});
			}
			elsif ($h eq "deja_vu") {
				$out .= $cgi->td($class_default,$hResVariants->{$var_id}->{table_dejavu});
			}
			elsif ($h eq 'table_validation') {
				$out .= $cgi->td($class_default,$hResVariants->{$var_id}->{table_validation});
			}
			elsif ($h eq 'table_transcript') {
				$out .= qq{<div style="max-height:200px;overflow-y:auto;max-width:220px;overflow-x:auto;">};
				$out .= $cgi->td($class, $hResVariants->{$var_id}->{table_transcript});
				$out .= qq{</div>};
			}
			elsif ($h eq '#') {
				$out.= $cgi->td($class_default,$nb_var);
			}
			elsif ($h eq 'alamut variant') {
				$out .= $cgi->td($class_default, "<center>".$hResVariants->{$var_id}->{alamut_link_variant}."</center>");
			}
			elsif ($h eq 'var_name') {
				$out .= $cgi->td($class_default, "<center>".$hResVariants->{$var_id}->{table_vname}."</center>");
			}
			elsif ($h eq 'varsome') {
				$out .= $cgi->td($class_default, "<center>".$hResVariants->{$var_id}->{table_varsome}."</center>");
			}
			else {
				$out.= $cgi->td($class_default,$hvariation->{html}->{$h});
			}
		}
		
		$out .= $cgi->end_Tr();
		push(@lTrLines, $out);
	}
#	}
	
#	warn Dumper @lTrLines;
#	die;
	
	if ($is_polybtf) {
		$hRes->{html_gene_polybtf} = "polybtf";
		foreach my $line (@lTrLines) { $out2 .= $line; }
	}
	else {
		if ($gene_init->strand() == 1) { foreach my $line (@lTrLines) { $out2 .= $line; } }
		else { foreach my $line (reverse @lTrLines) { $out2 .= $line; }  }
		$h_count->{total_pass} = scalar(@lTrLines);
		$gene_init = undef;
		$gene_init = $project_dejavu->newGene($gene_init_id) unless ($gene_init);
		$gene_init->hgmd();
		$gene_init->omim();
		$gene_init->is_omim_morbid();
	
		$hRes->{html_gene} = get_html_gene($gene_init_id, $gene_init, $h_count);
	}
	
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
	$hRes->{gencode_version} = $project_init->gencode_version;
	$hRes->{session_id} = $session_id;
	
	if ($is_polybtf) {
		$hRes->{html_page_title} = 'PolyBTF';
		$hRes->{html_source} = 'PolyBTF';
		$hRes->{html_title} = 'PolyBTF';
	}
	else {
		$hRes->{html_page_title} = 'DV '.$gene_init->external_name();
		$hRes->{html_source} = 'PolyDejaVu';
		$hRes->{html_title} = 'Gene '.$gene_init->external_name();
	}
	
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
	$h_annot_categories->{'only_strict_ill_patients'} = 1 if ($only_strict_ill);
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
	if (exists $hRes->{'html_gene_polybtf'}) { $session->save( 'html_gene', $hRes->{'html_gene_polybtf'} ); }
	elsif (exists $hRes->{'html_gene'}) { $session->save( 'html_gene', $hRes->{'html_gene'} ); }
	$session->save( 'html_variants', $hRes->{'html_variants'} );
	$session->save( 'hash_filters', $hRes->{'hash_filters'} );
	$session->save( 'gencode_version', $hRes->{'gencode_version'} );
}

sub get_variants_infos_from_polybtf {
	my ($h_count, $hVariantsProj);
	print '_polybtf_dv_';
	my ($hVariantsDetails) = update_list_variants_from_polybtf();
	return ($hVariantsDetails);
}

sub get_variants_infos_from_projects {
	my ($use_locus, $only_transcript, $only_rocks_id) = @_;
	my ($h_count, $hVariantsProj);
	print '_update_dv_';
	my ($hVariantsDetails) = update_list_variants_from_dejavu($use_locus, $only_transcript, $only_rocks_id);
	return ($hVariantsDetails);
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
	$hResGene->{$gene_init_id}->{phenotypes} = $pheno.';'.$nb_other_terms;
	my $description_gene = $gene_init->description();
	$gene_init->getProject->buffer->dbh_deconnect();
	$gene_init->getProject->buffer->dbh_reconnect();
	foreach my $panel (@{$gene_init->getPanels()}) {
		$hResGene->{$gene_init_id}->{panels}->{$panel->name()}->{phenotype} = $panel->getPhenotypes()->[0]->name();
	}
	my $html_gene = update_variant_editor::panel_gene($hResGene->{$gene_init_id});
	$html_gene =~ s/CNV//;
	my $regexp1 = qq{<span class=\'badge badge-infos badge-xs \' style="color:#00C851"  >[0-9]+ </span>};
	my $regexp2 = qq{<span class=\'badge badge-infos badge-xs \' style="color:#00C851"  > $nb_selected/$nb_variants_before_filters </span>};
	$html_gene =~ s/$regexp1/$regexp2/;
	if ($description_gene) {
		$html_gene =~ s/<\/table>//;
	 	$html_gene .= qq{<td><span style='color:white;'><i><center>$description_gene</center></i></span></td></table>};
	}
	return $html_gene;
}

sub update_list_variants_from_polybtf {
	my $time = time;
	my $h_dv_rocks_ids;
	my $sql = "SELECT * FROM '/data-isilon/public-data/repository/HG38/polybtf/20/new_hgmd_clinvar.parquet';";
	my $duckdb = $buffer_init->software('duckdb');
	my $cmd = qq{set +H | $duckdb -json -c "$sql"};
	my $json_duckdb = `$cmd`;
	if ($json_duckdb) {
		my $decode = decode_json $json_duckdb;
		
		my $h_found;
		foreach my $h_json (@$decode) {
			my $chr_id = $h_json->{'chr'};
			my $rocksid = $h_json->{'rocksid'};
			$h_found->{$chr_id}->{$rocksid} = $h_json;
		}
		
		foreach my $chr_id (sort keys %$h_found) {
			
			$buffer_init = new GBuffer;
			$project_init = $buffer_init->newProject( -name => $project_name_hg38 );
			my $chr = $project_init->getChromosome($chr_id);
			my $no = $chr->rocks_dejavu();
			foreach my $rocksid (sort keys %{$h_found->{$chr_id}}) {
#				warn $chr_id.' -> '.$rocksid;
				my $h_dv = $no->dejavu($rocksid);
				if ($h_dv) {
					$h_dv_rocks_ids->{$chr_id}->{$rocksid} = $h_dv;
					$h_dv_rocks_ids->{$chr_id}->{$rocksid}->{'polybtf'} = $h_found->{$chr_id}->{$rocksid};
				}
			}
			$no->close();
		}
	}
	
	($hVariantsDetails) = check_variants($h_dv_rocks_ids);
	print '@time:'.abs(time) - $time;
	
	warn Dumper keys %{$hVariantsDetails} if $debug;
	
	return ($hVariantsDetails);
}

sub update_list_variants_from_dejavu {
	my ($use_locus, $only_transcript, $only_rocks_id) = @_;
	my $time = time;
	my $buffer_dejavu = new GBuffer;
	$project_dejavu = $buffer_dejavu->newProject( -name => $buffer_init->getRandomProjectName());
	my $gene_dejavu = $project_dejavu->newGene($gene_init_id_for_newgene);
	
	my $h_dv_rocks_ids;
	if ($use_locus) {
		my (@lTmp) = split('_', $use_locus);
		$h_dv_rocks_ids->{$gene_dejavu->getChromosome->id()} = $gene_dejavu->getChromosome->rocks_dejavu->dejavu_interval($lTmp[1], $lTmp[2]);
	}
	else {
		$h_dv_rocks_ids->{$gene_dejavu->getChromosome->id()} = $gene_dejavu->getChromosome->rocks_dejavu->dejavu_interval($gene_dejavu->start(), $gene_dejavu->end());
	}
	($hVariantsDetails) = check_variants($h_dv_rocks_ids, $gene_dejavu);
	$buffer_dejavu = undef;
	$project_dejavu->buffer->dbh_deconnect();
	print '@time:'.abs(time) - $time;
	return ($hVariantsDetails);
}

sub check_variants {
	my ($h_dv_rocks_ids, $gene_dejavu) = @_;
	$project_init = undef;
	$buffer_init = undef;
	
	$buffer_init = new GBuffer;
	$project_init = $buffer_init->newProject( -name => $project_name_hg38 );
	
	my ($h_dv_var_ids, @lVarIds, @lVar);
	my ($total, $total_pass);
	foreach my $chr_id (sort keys %{$h_dv_rocks_ids}) {
		
		my $chr = $project_init->getChromosome($chr_id);
		foreach my $rocks_id (sort keys %{$h_dv_rocks_ids->{$chr_id}}) {
			my $nb_pat_he = 0;
			my $nb_pat_ho = 0;
			foreach my $proj_id (keys %{$h_dv_rocks_ids->{$chr_id}->{$rocks_id}}) {
				next if ($proj_id eq 'polybtf');
				$nb_pat_he += $h_dv_rocks_ids->{$chr_id}->{$rocks_id}->{$proj_id}->{he};
				$nb_pat_ho += $h_dv_rocks_ids->{$chr_id}->{$rocks_id}->{$proj_id}->{ho};
			}
			my $nb_pat_total = $nb_pat_he + $nb_pat_ho;
			my $var_id = $chr->transform_rocksid_to_varid($rocks_id);
			push(@lVarIds, $var_id);
			push(@lVar, $project_init->_newVariant($var_id));
			$h_dv_var_ids->{$var_id}->{rocks_id} = $rocks_id;
		}
	
		my $lift = liftOver->new(project=>$project_init, version=>$project_init->lift_genome_version);
		$lift->lift_over_variants(\@lVar);
	}
		
	$nb_variants_before_filters = scalar(keys %{$h_dv_var_ids});
	print '.!nbRocksIds:'.$nb_variants_before_filters.'!.';
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
	print '.use_fork'.$fork.'.';
	
	$project_init->disconnect();
	
	my $transcript_dejavu;
	if ($gene_dejavu) {
		foreach my $tr (@{$gene_dejavu->getTranscripts()}) {
			$transcript_dejavu = $tr if ($only_transcript and $tr->id() =~ /$only_transcript/);
			$gene_with_partial_transcrit = 1 if ($tr->is_partial_transcript());
		}
	}
	
	my @lProjectNames;
	my $hVarErrors;
	
	my $pm = new Parallel::ForkManager($fork);
	my $nb_errors=0;
	$pm->run_on_finish(
		sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data) = @_;
			delete $data->{start_job};
			if ($fork == 1) {
				$hVariantsDetails = $data;
				$total_pass = scalar(keys %{$hVariantsDetails}) - 2;
			}
			else {
				foreach my $var_id (keys %{$data}) {
					$hVariantsDetails->{$var_id} = $data->{$var_id};
					$total_pass++;
				}
			}
		}
	);
	
	my $nb = int((scalar(@lVar)+1)/($fork));
	$nb = 20 if $nb == 0;
	my $iter = natatime($nb, @lVar);
	while( my @tmp = $iter->() ){
 	 	my $pid = $pm->start and next;
 	 	my $hres;
 	 	$hres->{start_job} = 1;
		my $nb_i;
		foreach my $var (@tmp) {
			my $var_id = $var->id;
			
			
#			warn "\n\n";
#			warn $var_id;
			
			my $var_id_hg19 = $var->lift_over('HG19')->{id};
			
			my $gene_variant;
			
			my $is_ok_gene;
			foreach my $gene (@{$var->getGenes()}) {
				$is_ok_gene = 1 if $gene->id eq $gene_dejavu->id;
			} 
			if ($is_ok_gene) { $gene_variant = $gene_dejavu; }
			else { $gene_variant = $var->getGenes->[0]; }
			$total++;
			$nb_i++;
			if ($nb_i == 100) {
				print '.';
				$nb_i = 0;
			}
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
			if (not $is_ok_perc) {
				delete $hres->{$var_id};
				next;
			}
			warn '1 - ok perc' if $debug;
			
			my ($var_gnomad, $var_gnomad_ho, $var_annot, $var_dejavu, $var_dejavu_ho, $var_model);
			$var_id = uc($var_id);
			if ($var_id =~ /ALU/) {
				delete $hres->{$var_id};
				next;
			}
			if (not $var_id =~ /[XYMT0-9]+_[0-9]+_[ATGC]+_[ATGC]+/) {
				delete $hres->{$var_id};
				next;
			}
			warn '2 - ok id' if $debug;
			
			#my $var = $project_dejavu->_newVariant($var_id);
#			warn ref($var).' -> using: '.$var->id if $debug;
			if ($hres and exists $hres->{$var_id}) {
				$var_gnomad = $hres->{$var_id}->{var_gnomad} if ($hres->{$var_id}->{var_gnomad});
				$var_gnomad_ho = $hres->{$var_id}->{var_gnomad_ho} if ($hres->{$var_id}->{var_gnomad_ho});
				$var_dejavu = $hres->{$var_id}->{var_dejavu} if ($hres->{$var_id}->{var_dejavu});
				$var_dejavu_ho = $hres->{$var_id}->{var_dejavu_ho} if ($hres->{$var_id}->{var_dejavu_ho});
				if ($only_transcript) {
					$var_annot = $var->variationTypeInterface($transcript_dejavu);
				}
				else {
					$var_annot = $hres->{$var_id}->{annotation} if ($hres->{$var_id}->{annotation});
				}
				$var_model = $hres->{$var_id}->{model} if ($hres->{$var_id}->{model});
			}
			if ($var->isCnv() or $var->isLarge()) {
				delete $hres->{$var_id};
				next;
			}
			warn '3 - ok cnv large' if $debug;
			
	#		$var->{rocksdb_id} = $h_dv_var_ids->{$var_id}->{rocks_id};
	#		warn Dumper $gene_dejavu->getChromosome->getDejaVuInfosForDiagforVariant($var);
	#		
	#		$hres->{$var_id}->{dejavu}->{other_projects} = $var->other_projects;
	#		$hres->{$var_id}->{dejavu}->{other_patients} = $var->other_patients;
	#		$hres->{$var_id}->{dejavu}->{other_patients_ho} = $var->other_patients_ho;
			
			my $not_ok;
			
			unless ($var_gnomad) {
				$var_gnomad = $var->getGnomadAC();
			}
			
			$not_ok++ if ($var_gnomad and $max_gnomad and $var_gnomad > $max_gnomad);
			if ($not_ok) {
				delete $hres->{$var_id};
				next;
			}
			warn '4 - ok gnomad' if $debug;
			
			unless ($var_gnomad_ho) {
				$var_gnomad_ho = $var->getGnomadHO();
			}
			$not_ok++ if ($max_gnomad_ho and $var_gnomad_ho > $max_gnomad_ho);
			if ($not_ok) {
				delete $hres->{$var_id};
				next;
			}
			warn '5 - ok gnomad ho' if $debug;
			
			my $is_ok_annot;
			unless ($var_annot) {
				#$var->annotation();
				eval {
					$var->annotation();
					if ($var->{annotation}->{all}->{$var->getProject->getMaskCoding("intergenic")}) {
						$var_annot = 'intergenic';
					}
					elsif ($transcript_dejavu) { $var_annot = $var->variationTypeInterface($transcript_dejavu); }
					elsif (exists $var->{annotation}->{$gene_variant->id()}) { $var_annot = $var->variationTypeInterface($gene_variant); }
					else { $var_annot = $var->variationTypeInterface(); }
				};
				if ($@) {
					$var_annot = 'error';
					supressCoreFilesFound();
				}
			}
			warn $var_annot if $debug;
			foreach my $this_annot (split(',', $var_annot)) {
				$this_annot =~ s/ /_/g;
				$is_ok_annot ++ if (exists $h_filters_cons->{lc($this_annot)});
			}
			
			$is_ok_annot++ if $is_polybtf;
			if (not $is_ok_annot) {
				delete $hres->{$var_id};
				next;
			}
			warn '6 - ok annot' if $debug;
	
			my ($has_proj, $ok_model);
			
			warn 'rocks: '.$var->rocksdb_id if $debug;
			warn $h_dv_var_ids->{$var_id}->{rocks_id} if $debug;
			next if not exists $h_dv_rocks_ids->{$var->getChromosome->id()}->{$var->rocksdb_id};
	
			my @list_parquets;
			foreach my $proj_id (keys %{$h_dv_rocks_ids->{$var->getChromosome->id()}->{$var->rocksdb_id}}) {
				next if not exists $hProjectsIds->{$proj_id};
				my $proj_name = $hProjectsIds->{$proj_id};
				$hVariantsIdsDejavu->{$var_id}->{$proj_name} = $proj_id;
				foreach my $pat_id (@{$h_dv_rocks_ids->{$var->getChromosome->id()}->{$var->rocksdb_id}->{$proj_id}->{patients}}) {
					$hres->{$var_id}->{dejavu_details}->{$proj_name}->{$pat_id} = undef;
				}
				my $parquet = $dir_parquet.'/'.$proj_name.'.'.$proj_id.'.parquet';
				push(@list_parquets, "'".$parquet."'") if (-e $parquet);
			}
			if (not @list_parquets) {
				delete $hres->{$var_id};
				next;
			}
			
			#warn Dumper @list_parquets;

			warn '8 - ok dejavu ho' if $debug;
			
			
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
			update_variant_editor::vhgmd($var,$h_var);
			
			$h_var->{'html'}->{'no_css_polydiag'} = 1;
			my $val1 = 'onClick="zoomHgmd';
			my $Val2 = 'onClick="zoomHgmdWithoutCss';
			$h_var->{'html'}->{'hgmd'} =~ s/$val1/$Val2/;
			
			if ($var->hgmd) {
				$h_var->{value}->{hgmd_id} = $var->hgmd_id;
				my $n1 = $project_init_name;
			 	my $n2 = $var->hgmd_id;
			 	my $n3 = $var->id;
				my $cmd_hgmd = qq{zoomHgmdWithoutCss(\'$n1\',\'$n2\',\'$n3\')}; 
				$h_var->{html}->{hgmd} = update_variant_editor::printButton(4,[3,4], $var->hgmd->{class},qq{onClick="$cmd_hgmd"});
				if ($var->isDM()) {
					$h_var->{value}->{dm} = 1;
					$h_var->{value}->{dm_for_this_gene} = undef;
					$h_var->{value}->{dm} = undef;
					$h_var->{value}->{hgmd} = '';
					$h_var->{html}->{hgmd} = '';
				}
			}
			else { $h_var->{value}->{dm} = ''; }
			
			if ($h_var->{value}->{clinvar_pathogenic}) {
				my $clinvar_id = $h_var->{value}->{clinvar_id};
				if ($var->getChromosome->is_clinvar_pathogenic_for_gene($clinvar_id, $gene_variant)) {
					$h_var->{value}->{clinvar_pathogenic_for_this_gene} = 1;
				}
				else {
					$h_var->{value}->{clinvar_pathogenic_for_this_gene} = undef;
					$h_var->{value}->{clinvar_pathogenic} = undef;
					$h_var->{value}->{clinvar} = '';
					$h_var->{html}->{clinvar} = '';
				}
			}
			
			update_variant_editor::vspliceAI($var, $h_var);
			$hres->{$var_id}->{spliceAI} = $h_var->{html}->{spliceAI}->{$gene_variant->id};
							
			update_variant_editor::vhgmd($var, $h_var);
			$hres->{$var_id}->{table_validation} = update_variant_editor::table_validation_without_local($var->getProject, $h_var, $gene_variant);
			
			
			$hres->{$var_id}->{table_gnomad} = update_variant_editor::table_gnomad($var);
			$hres->{$var_id}->{table_gnomad} =~ s/gnomad_r2_1/gnomad_r4/;
			
			$hres->{$var_id}->{table_varsome} = update_variant_editor::vvarsome($h_var);
			
			$hres->{$var_id}->{table_dejavu} = update_variant_editor::vdejavu($var, $h_var);
			
			$h_var->{genes}->{$gene_variant->id} = update_variant_editor::construct_hash_transcript($var, $cgi, \@header_transcripts, 2, $gene_variant);
			
			
			if ($is_polybtf or not $is_ok_gene) {
				my $gene_name = $gene_variant->external_name();
				my $gene_id = $gene_variant->id();
				my ($pheno,$nb_other_terms) = $gene_variant->polyviewer_phentotypes();
				my $html_g = qq{<div><table>};
				$html_g .= qq{<tr><td>Gene <b>$gene_name</b> <i>[$gene_id]</i></td></tr> };
				$html_g .= qq{<tr><td>$pheno</td></tr> };
				$html_g .= qq{</table></div>};
				$hres->{$var_id}->{table_transcript} = $html_g."<br>".update_variant_editor::table_transcripts($h_var->{genes}->{$gene_variant->id}, \@header_transcripts, 1);
			}
			else {
				$hres->{$var_id}->{table_transcript} = update_variant_editor::table_transcripts($h_var->{genes}->{$gene_variant->id}, \@header_transcripts, 1);
			}
			
			
			my $gnomad_id_hg38 = $var->gnomad_id;
			my $html_vname_hg38 = update_variant_editor::vname2($var, $h_var);
			my $gnomad_id_hg19 = $var->lift_over('HG19')->{name};
			my $html_vname_hg19 = $html_vname_hg38;
			$html_vname_hg19 =~ s/gnomad_r4/gnomad_r2_1/g;
			$html_vname_hg19 =~ s/$gnomad_id_hg38/$gnomad_id_hg19/g;
			my $html_vname = qq{
				<table>
					<center>
						<tr><td><b>HG38:</b></td><td style="padding-left:10px;">$html_vname_hg38</td></tr>
						<tr><td><b>HG19:</b></td><td style="padding-left:10px;">$html_vname_hg19</td></tr>
					</center>
				</table>
			};
			$hres->{$var_id}->{table_vname} = $html_vname;
			
			if ($only_transcript) {
				my @new_list;
				foreach my $htr (@{$h_var->{genes}->{$gene_init_id_for_newgene}}) {
					if ($htr->{value}->{trid} eq $only_transcript) {
						push(@new_list, $htr);
					}
				}
				$h_var->{genes}->{$gene_init_id_for_newgene} = \@new_list;
			}
			$hres->{$var_id}->{alamut_link_variant} = html_polygenescout::print_alamut_variant_button($var->alamut_id());
			$hres->{$var_id}->{var_gnomad} = $var_gnomad;
			$hres->{$var_id}->{var_gnomad_ho} = $var_gnomad_ho;
			$hres->{$var_id}->{var_dejavu} = $var_dejavu;
			$hres->{$var_id}->{var_dejavu_ho} = $var_dejavu_ho;
			$hres->{$var_id}->{annotation} = $var_annot;
			$hres->{$var_id}->{id_hg19} = $var_id_hg19;
			
			my $table_projects_patients = get_from_duckdb_project_patients_infos($var, \@list_parquets);
			if ($table_projects_patients) {
				$hres->{$var_id}->{table_projects_patients} = $table_projects_patients;
			}
			else {
				delete $hres->{$var_id};
			}
		}
		
	 	$pm->finish(0, $hres);
	}
	sleep(3); 
	$pm->wait_all_children();
	
#	warn Dumper $hVariantsDetails;	
#	die
	
	print 'nbVarPass:'.$total_pass;
	print 'nbProj:'.scalar(@lProjectNames);
	return ($hVariantsDetails);
}

sub get_table_project_patients_infos {
	my ($project_name, $hash, $hVar_infos) = @_;
	my $locus_hg19 = $hVar_infos->{locus_hg19};
	my $locus_hg38 = $hVar_infos->{locus_hg38};
	
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
			my $ratio = ($info / $h_infos_patients->{$nb_pat}->{dp}) * 100;
			my $text = 'AC:'.$info.' ('.int($ratio).'%)';
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
	
	my $b = new GBuffer;
	my $p = $b->newProject( -name => $project_name);
		
#	if ($p->name eq 'NGS2025_08773') {
#		warn "\n";
#		warn Dumper $h_tmp_pat;
#	}
	
	my $found_healthy_patient;
	foreach my $pat (@{$p->getPatients()}) {
		if (not $pat->isIll() and $only_ill) {
			$found_healthy_patient = 1;
			delete $h_infos_patients->{$h_tmp_pat->{$pat->id}};
			next;
		}
		next if not exists $h_tmp_pat->{$pat->id};
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{name} = $pat->name;
#		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{status} = 'healthy';
#		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{status} = 'ill' if $pat->isIll();
		
		my $icon = $pat->small_icon();
		$icon =~ s/"/'/g;
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{status} = $icon;
		
		if (int($h_tmp_pat->{$pat->id}) <= $nb_he) { $h_infos_patients->{$h_tmp_pat->{$pat->id}}->{heho} = 'He'; }
		else { $h_infos_patients->{$h_tmp_pat->{$pat->id}}->{heho} = 'Ho'; }
	}
	
	return undef if not $h_infos_patients or scalar keys %$h_infos_patients == 0;
	return undef if ($only_strict_ill and $found_healthy_patient);
	
#	my $gene_name =$gene_used->external_name();
#	my $fam = $patient->getFamily();
	my $is_solo_trio = 'SOLO';
#	$is_solo_trio = 'TRIO' if $fam->isTrio();
#	my $project_name = $patient->getProject->name();
	
	my $description = $p->description();
	my @l_users = @{$p->get_list_emails()};
#	my $patient_name = $patient->name();
#	my $pheno = undef;
#	$pheno = $h_var->{html}->{pheno_name} if ($h_var and exists $h_var->{html}->{pheno_name});
	
	my $color = "#c1c1c1";
	my $model;
	my $patient_heho = "-";
	
	my $nb_col_span = 6;
#	$nb_col_span = 7 if ($var->getProject->isGenome() && $var->isCnv);
#	my $hstatus = $patient->validation_status();
#	my $hval = $patient->validations();
	
	
	my $table_trio = qq{ <div> };
	$table_trio .= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:3px"});
	
	
	my @lPhenotypes = @{$p->phenotypes()};
	my $pheno = join(', ', sort @lPhenotypes);
	
	my $version = $p->annotation_genome_version();
	my $color_version = '#85f283';
	$color_version = '#f2c37c' if $version eq 'HG38';
	
	my $users = join("<br>", @l_users);
	my $proj_text = qq{<button style="color:black;" onclick="get_popup_users('$users');">Users</button> - <b>$project_name</b>};
	$proj_text .= qq{<sup>defidiag</sup>} if $p->isDefidiag; 
	$proj_text .= " - <span style='color:#82d0f5;'>$pheno</span>" if ($pheno);
	$proj_text .= "<br>";
	$proj_text .= "<b><span style='color:$color_version;'>$version</span></b>  ";
	$proj_text .= "$description";
		
	my (@l_pat_names, @l_pat_bam);
	foreach my $id (sort keys %$h_infos_patients) {
		next if not (exists $h_infos_patients->{$id}->{name});
		push(@l_pat_names, $h_infos_patients->{$id}->{name});
		push(@l_pat_bam, $p->getPatient($h_infos_patients->{$id}->{name})->bamUrl());
	}
	my $pnames = join(';', @l_pat_names);
	my $f = join(';', @l_pat_bam);
	my $gn = $p->getVersion();
	my $chr_name = $hVar_infos->{chr_id};
	my ($locus, $start);
	if ($p->annotation_genome_version() eq 'HG19') {
		$locus = $locus_hg19;
		$start = $hVar_infos->{start_hg19};
	}
	elsif ($p->annotation_genome_version() eq 'HG38') {
		$locus = $locus_hg38;
		$start = $hVar_infos->{start_hg38};
	}
	my $a0 = $hVar_infos->{ref_all};
	my $a1 = $hVar_infos->{alt_all};
	my $igb_b = qq{<button class='igvIcon2' onclick='launch_web_igv_js("$project_name","$pnames","$f","$locus","/","$gn")' style="color:black"></button>};	
	my $alamut_b = qq{<button class="alamutView3" onClick ="displayInAlamut('$chr_name',$start,['$a0','$a1']);"></button>};

	my $no_header_project ;
	my $no_header_project_pat;
	
	$table_trio .= $cgi->start_Tr({style=>"background-color:#949292; color:white;"});
	my $b_others_var;
	my $pat_name = 'ALL';
	my $var_id = 'VARID';
	my $father_trans = undef;
	my $mother_trans = undef;
	my $other_trans = undef;
	
	my $img_child = qq{<img style='width:14px;height:14px;' src='/icons/Polyicons/baby-boy.png'>};
	#$img_child = qq{<img style='width:14px;height:14px;' src='/icons/Polyicons/baby-girl.png'>} if ($patient->sex() == 2);
	my $cmd_others = qq{view_var_from_proj_gene_pat('$project_name','$gene_init_id','$pat_name','$var_id');};
	$other_trans = qq{<button style="text-align:middle;vertical-align:top;width:28px;border:solid 0.5px black;background-color:white;" onClick="$cmd_others">$img_child</button>};

	$table_trio .= $cgi->td({colspan=>($nb_col_span)-1}, $proj_text);
	$table_trio .= $cgi->td({colspan=>1, style=>"padding:0px;width:60px;"}, '<center>'.$igb_b.' '.$alamut_b.'</center>');
	
	$no_header_project_pat = 1;
	$table_trio .= $cgi->end_Tr();
	if ($no_header_project_pat) {
		$table_trio .= '</table>';
		$table_trio .= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:3px"});
	}
	
	foreach my $nb (keys %{$h_infos_patients}) {
		next if not exists $h_infos_patients->{$nb}->{name};
		my $patient_name = $h_infos_patients->{$nb}->{name};
		my $patient = $p->getPatient($patient_name);
		my $patient_status = $h_infos_patients->{$nb}->{status};
		my $patient_heho = $h_infos_patients->{$nb}->{heho};
		my $dp = 'DP:'.$h_infos_patients->{$nb}->{dp};
		my $perc_allele = $h_infos_patients->{$nb}->{ratio};
		my $model = $h_infos_patients->{$nb}->{model};
		
		if ($model eq 'denovo' or $model eq 'strict_denovo' or $model eq 'dominant') {
			$color = 'background-color:#f54e4e';
			$model = ucfirst($model);
		}
		elsif ($model eq 'recessif') {
			$color = 'background-color:#e99ff7';
			$model = ucfirst($model);
		}
		elsif ($model eq 'father') {
			$color = 'background: linear-gradient(to right, white, #ddfbff);';
			if ($p->getPatient($patient_name)->getFamily->getFather->isIll()) { $model = qq{<img src="/icons/Polyicons/male-d.png">}; }
			else { $model = qq{<img src="/icons/Polyicons/male-s.png">}; }
		}
		elsif ($model eq 'mother') {
			$color = 'background: linear-gradient(to right, white, #ffddfd);';
			if ($p->getPatient($patient_name)->getFamily->getMother->isIll()) { $model = qq{<img src="/icons/Polyicons/female-d.png">}; }
			else { $model = qq{<img src="/icons/Polyicons/female-s.png">}; }
		}
		elsif ($model eq 'both') {
			$color = 'background: linear-gradient(to right, #ddfbff, #ffddfd);';
			if ($p->getPatient($patient_name)->getFamily->getFather->isIll()) { $model = qq{<img src="/icons/Polyicons/male-d.png">}; }
			else { $model = qq{<img src="/icons/Polyicons/male-s.png">}; }
			if ($p->getPatient($patient_name)->getFamily->getMother->isIll()) { $model .= qq{<img style="padding-left:5px;" src="/icons/Polyicons/female-d.png">}; }
			else { $model .= qq{<img style="padding-left:5px;" src="/icons/Polyicons/female-s.png">}; }
			 
		}
		elsif ($model eq 'solo') {
			$color = 'background: linear-gradient(to right, white, #c7c6c5);';
			$model = ucfirst($model);
		}
		else { $color = 'background-color:white'; }
	
		my $local_text;
		$table_trio .= $cgi->start_Tr({style=>"$color;"});
		$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, "<span style='text-align:left;'>$patient_name $local_text</span>");
		$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, "<div>".$patient_status."</div>");
		$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, "$patient_heho");
	#	if ($var->getProject->isGenome() && $var->isCnv) {
	#		$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, 'pr:'.$var->pr($patient));
	#		$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, 'sr:'.$var->sr($patient));
	#		eval {
	#			my $cnv_score = sprintf("%.2f", log2($patient->cnv_value_dude($var->getChromosome->name,$var->start,$var->start+$var->length)));
	#			$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, 'cnv_score:'.$cnv_score);
	#		};
	#		if ($@) {
	#			$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, 'cnv_score:pb');
	#		}
	#	}
	#	else {
			#my $perc_allele = $var->getPourcentAllele($patient);
			
#			return 'perc_all_filter' if ($filter_perc_allelic_min and $perc_allele < $filter_perc_allelic_min);
#			return 'perc_all_filter' if ($filter_perc_allelic_max and $perc_allele > $filter_perc_allelic_max);
#			$perc_allele .= "%" if ($perc_allele ne '-');
			$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, $perc_allele);
			$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, $dp);
	#	}
		$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, $model);
		$table_trio .= $cgi->end_Tr();
	#	if ($local_text_tab) {
	#		$table_trio .=qq{ <div onClick="document.getElementById('$table_validation_id').style.display='none';" id='$table_validation_id' style='display:none;'>$local_text_tab</div> };
	#	}
	}
	$table_trio .= "</table></div>";
	
	return ($table_trio, lc($model));
}

sub get_from_duckdb_project_patients_infos {
	my ($var, $list_files) = @_;
	return if scalar(@$list_files) == 0;
	my @list_table_trio;
	my $sql = "SELECT * FROM read_parquet([".join(', ', @$list_files)."])";
	my $find_pos_s = $var->start() - 20;
	my $find_pos_e = $var->start() + 20;
	if ($var->getProject->current_genome_version() eq 'HG38') {
		$sql .= " WHERE chr38='".$var->getChromosome->id()."' and pos38 BETWEEN '".$find_pos_s."' and '".$find_pos_e."';" ;
		my $duckdb = $buffer_init->software('duckdb');
		my $cmd = qq{set +H | $duckdb -json -c "$sql"};
		my $json_duckdb = `$cmd`;
		if ($json_duckdb) {
			my $decode = decode_json $json_duckdb;
			my $h_by_proj;
			foreach my $h (@$decode) {
				next if $h->{'chr38'} ne $var->getChromosome->id;
				my $var_start = $var->start();
				$var_start-- if $var->isInsertion() or $var->isDeletion();
				next if $h->{'pos38'} ne $var_start;
				my $var_all = $h->{'allele'};
				$var_all =~ s/\+//;
				next if $var_all ne $var->var_allele();
				my $project_id = $h->{project};
				my $project_name = $hProjectsIds->{$project_id};
				$h_by_proj->{$project_name} = $h;
			}
			my $hVar_infos;
			$hVar_infos->{locus_hg19} = $var->getChromosome->id().":".$var->lift_over('HG19')->{position}."-".$var->lift_over('HG19')->{position};;
			$hVar_infos->{locus_hg38} = $var->getChromosome->id().":".$var->start."-".$var->end;
			$hVar_infos->{start_hg19} = $var->lift_over('HG19')->{position};
			$hVar_infos->{start_hg38} = $var->start;
			$hVar_infos->{chr_id} = $var->getChromosome->id();
			$hVar_infos->{ref_all} = $var->ref_allele();
			$hVar_infos->{alt_all} = $var->var_allele();
			$hVar_infos->{alt_all} = '*' unless $var->var_allele();
			foreach my $project_name (reverse sort keys %$h_by_proj) {
				my $h = $h_by_proj->{$project_name};
				my ($table_trio, $model) = get_table_project_patients_infos($project_name, $h, $hVar_infos);
				push(@list_table_trio, $table_trio) if $table_trio;
			}
		}	
	}
	else {
		warn "HG19";
		confesss("HG19!!");
	}
	return join("<br>",@list_table_trio) if scalar(@list_table_trio) >= 1;
	return undef;
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
	my ($h_projects, @list_hash);
	
	
	
	
	if ($buffer_init->getQuery->isUserMagic($user_name, $pwd)) {
		@list_hash = @{$buffer_init->getQuery()->getAllProjects()};
		
#		my $h_project_dejavu;
#		foreach my $proj_name (@{$buffer_init->getQuery()->listProjectsForDejaVu()}) {
#			$h_project_dejavu->{$proj_name} = undef;
#		}
#		
#		
#		foreach my $h (@list_hash) {
#			next if not exists $h_project_dejavu->{$h->{name}};
#			next if not $h->{name} =~ /NGS20/;
#			my $ok;
#			my $file1 = $dir_parquet.'/'.$h->{name}.'.'.$h->{id}.'.parquet';
#			$ok = 1 if -e $file1;
#			my $file2 = $dir_parquet.'/'.$h->{name}.'.'.$h->{id}.'.parquet.no_dejavu';
#			$ok = 1 if -e $file2;
#			
#			if (not $ok) {
#				warn $h->{name};
#				warn $file1;
#			}			
#		}
#		die;
	}
	else {
		@list_hash = @{$buffer_init->getQuery()->getProjectListForUser($user_name, $pwd)};
	}
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


