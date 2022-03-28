#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;

use strict;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/validation_variation"; 
use lib "$Bin/../cache_nodb/scripts/";
use GBuffer;
use GBufferTest;
use export_data;
use JSON;
use VcfMerge;
use GenBoNoSql;
use Set::IntervalTree;
use Spreadsheet::WriteExcel;
use Compress::Snappy;
use Storable qw(store retrieve freeze dclone thaw);
use POSIX qw(strftime);
use CGI::Session;
use html; 
use Carp;
use Cache_Commons;
use QueryVectorFilter;


my $io = IO::Handle->new();
$io->autoflush(1);

my $cgi = new CGI();
my $infos					= $cgi->param('infos');
my $check_genes				= $cgi->param('check_genes');
my $only_genes				= $cgi->param('only_genes');
my $bundle_id				= $cgi->param('bundle_id');
my $get_bundles				= $cgi->param('get_bundles');
my $genes_list				= $cgi->param('genes');
my $filter_text				= $cgi->param('filter_text');
my $filter_chromosome		= $cgi->param('filter_chromosome');
my $filter_type_variation 	= $cgi->param('filter_type_variation');
my $filter_type_variation_2	= $cgi->param('filter_type_variation_2');
my $filter_attic 			= $cgi->param('filter_attic');
my $filter_attic_genes		= $cgi->param('filter_attic_genes');
my $filter_patient 			= $cgi->param('filter_patient');
my $filter_not_patient 		= $cgi->param('filter_not_patient');
my $fam_not					= $cgi->param('fam_not');
my $fam_and 				= $cgi->param('fam_and');
my $filter_he 				= $cgi->param('filter_he');
my $filter_ho 				= $cgi->param('filter_ho');
my $filter_region			= $cgi->param('filter_region');
my $filter_nbvar_regionho	= $cgi->param('filter_nbvar_regionho');
my $filter_regionho_sub_only= $cgi->param('filter_region_ho_rec_only_sub');
my $filter_diseases			= $cgi->param('filter_diseases');
my $filter_genes_intersect	= $cgi->param('filter_genes_intersect');
my $level_fam 				= $cgi->param('level_fam');
my $level_ind 				= $cgi->param('level_ind');
my $atLeast 				= $cgi->param('nb');
my $typeFilters 			= $cgi->param('mode');
my $model 					= $cgi->param('model');
my $model_2					= $cgi->param('model_2');
my $projectName 			= $cgi->param('project');
my $polyscore				= $cgi->param('polyscore');
my $dejavu					= $cgi->param('dejavu');
my $dejavu_2				= $cgi->param('dejavu_2');
my $dejavu_ho				= $cgi->param('dejavu_ho');
my $getPatientsCgi			= $cgi->param('get_patients');
my $xls_by_genes			= $cgi->param('xls_by_genes');
my $xls_by_variants			= $cgi->param('xls_by_variants');
my $xls_by_regions_ho		= $cgi->param('xls_by_regions_ho');
my $xls_outfile				= $cgi->param('xls_outfile');
my $resume_vector_results   = $cgi->param('resume_vector');
my $export_vcf_for			= $cgi->param('export_vcf_for');
my $view_all_genes			= $cgi->param('view_all_genes');
my $fork					= $cgi->param('fork');
my $no_verbose				= $cgi->param('no_verbose');
my $xls_save_session		= $cgi->param('xls_save_session');
my $xls_load_session		= $cgi->param('xls_load_session');
my $gene_atlas_view			= $cgi->param('gene_atlas_view');
my $filter_confidence		= $cgi->param('filter_confidence');
my $not_construct_genes		= $cgi->param('not_construct_genes');
my $test_with_db			= $cgi->param('test_with_db');
my $test					= $cgi->param('test');
my $path_tmp				= $cgi->param('tmp_dir');
my $export_list_var_ids		= $cgi->param('export_list_var_ids');
my $detail_project			= $cgi->param('detail_project');
my $debug					= $cgi->param('debug');
my $user_name				= $cgi->param('user_name');
my $without_stats			= $cgi->param('without_stats');
my $delete_models			= $cgi->param('delete_models');
my $filter_gnomad			= $cgi->param('gnomad');
my $filter_ratio			= $cgi->param('ratio');
my $filter_type_ratio		= $cgi->param('type_ratio');
my $filter_ncboost			= $cgi->param('ncboost');
my $panel_name				= $cgi->param('panel');
my $annot_version			= $cgi->param('annot_version');
my $keep_indels_cadd		= $cgi->param('keep_indels_cadd');
#my $filter_gnomad_test		= $cgi->param('gnomad_test');



my $hDeleteModels;
if ($delete_models) {
	foreach my $model (split(' ', $delete_models)) {
		$hDeleteModels->{$model} = undef;
	}
}

if ($xls_outfile eq 'none') {
		warn "\n\nERROR: please replace value 'none' for last -xls_outfile option.\nExample: -xls_outfile=export.xls\n\n";
	die;
}

if ($getPatientsCgi) {
	my $buffer = new GBuffer;
	my $project = $buffer->newProject( -name => $projectName );
	my @lPatNames;
	foreach my $pat (@{$project->getPatients()}) { push(@lPatNames, $pat->name()) }
	my $hash;
	$hash->{'patients'} = join(',', sort(@lPatNames));
	print $cgi->header('text/json-comment-filtered');
	print encode_json $hash;
	exit(0);
}

$typeFilters = 'individual' if ($typeFilters eq 'ind');
$typeFilters = 'familial' if ($typeFilters eq 'fam');
$typeFilters = 'somatic' if ($typeFilters eq 'somatic');
$filter_type_variation =~ s/ /,/g;	

my $buffer;
if ($test) {
	$buffer = new GBufferTest;
	$buffer->dir_project($test);
}
else {$buffer = new GBuffer; }

my $project = $buffer->newProjectCache( -name 			=> $projectName,
									    -cache 		=> '1',
									    -typeFilters 	=> $typeFilters, );
									    

if ($user_name) {
	#$buffer->getQuery->updateLastConnectionUserProject(lc($user_name), $project->id()) unless ($test or $get_bundles or $check_genes);
}

if ($annot_version) {
	$project->changeAnnotationVersion($annot_version);
}

my $hChr;					    
if ($panel_name) {
	my @lGenes;
	foreach my $gene (@{$project->getPanel($panel_name)->getGenes()}) {
		push(@lGenes, $gene->external_name());
		$hChr->{$gene->getChromosome->id()}->{$gene->id()} = undef;
	}
	$only_genes = join(',', @lGenes);
	$filter_chromosome = join(',', keys %$hChr);
}
if ($export_vcf_for) {
	$filter_chromosome = 'all';
}

my $can_use_hgmd = $project->hasHgmdAccess($user_name);

if ($check_genes) { launch_check_genes_request(); }

if ($model) { $project->model($model); }

if ($not_construct_genes) { $project->not_construct_genes(1); }

if ($filter_attic_genes) {
	my $h;
	foreach my $pat_name (split(' ', $filter_attic_genes)) {
		$h->{$pat_name} = 0;
	}
	if ($filter_attic) {
		foreach my $pat_name (split(' ', $filter_attic)) {
			$h->{$pat_name} = 1;
		}
	}
	# empeche les patients / familles d enlever le statut alors qu on la specifiquement demande, en depit de l intersect genes
	$filter_attic = join(' ', keys %$h);
	foreach my $pat_name (keys %$h) {
		if ($h->{$pat_name} == 1) { delete $h->{$pat_name}; }
	}
	my @lFilters = ($filter_patient, $filter_not_patient, $fam_not, $fam_and, $filter_he, $filter_ho);
	foreach my $type_filter (@lFilters) {
		next unless ($type_filter);
		foreach my $f (split(' ', $type_filter)) {
			if (exists $h->{$f}) { delete $h->{$f}; }
		}
	}
	$project->filter_attic_genes($h);
}

if ($filter_attic) {
	foreach my $pat_name (split(' ', $filter_attic)) {
		$project->getPatient($pat_name)->in_the_attic(1);
	}
}
 
if ($filter_patient) {
	foreach my $pat_name (split(' ', $filter_patient)) {
		$project->getPatient($pat_name)->intersected(1);
	}
}
 
if ($filter_not_patient) {
	foreach my $pat_name (split(' ', $filter_not_patient)) {
		$project->getPatient($pat_name)->excluded('all');
	}
}
 
if ($filter_he) {
	foreach my $pat_name (split(' ', $filter_he)) {
		$project->getPatient($pat_name)->excluded('he');
	}
}
 
if ($filter_ho) {
	foreach my $pat_name (split(' ', $filter_ho)) {
		$project->getPatient($pat_name)->excluded('ho');
	}
}

if ($fam_and) {
	foreach my $fam_name (split(' ', $fam_and)) {
		$project->getFamily($fam_name)->intersected(1);
	}
}

if ($xls_by_variants) 	{ $project->isXlsOutput('variants'); }
if ($xls_by_genes)    	{ $project->isXlsOutput('genes'); }
if ($xls_by_regions_ho) { $project->isXlsOutput('regions_ho'); }

if ($filter_genes_intersect) {
	my $h;
	foreach my $g_id (split(',', $filter_genes_intersect)) {
		$h->{$g_id} = undef;
	}
	$project->filter_genes_intersect($h);
}

if ($get_bundles) { return launch_bundles_request(); }

if ($bundle_id) { return launch_uniq_bundle_request(); }

if ($only_genes) {
	foreach my $name (split(',', $only_genes)) {
		if ($name =~ /capture/) {
			my @lTmp = split(':', $name);
			my $capture_name = $lTmp[-1];
			foreach my $parent_capture_name (keys %{$buffer->getAllGenesNamesInAllBundle()}) {
				next unless (exists $buffer->getAllGenesNamesInAllBundle->{$parent_capture_name}->{$capture_name});
				foreach my $gene_name (keys %{$buffer->getAllGenesNamesInAllBundle->{$parent_capture_name}->{$capture_name}->{genes}}) {
					$project->{only_genes}->{$gene_name} = undef;
				}
			}
		}
		else { $project->{only_genes}->{$name} = undef; }
	}
}

if ($test_with_db) { $test = 1; } 
my ($hFiltersChr, $hFiltersChr_var2);
foreach my $filter_name (split(',', $filter_type_variation)) {
	if ($filter_name eq 'cnv') {
		$hFiltersChr->{'large_deletion'} = undef;
		$hFiltersChr->{'large_duplication'} = undef;
	}
	elsif ($filter_name eq 'upstream_downstream') {
		$hFiltersChr->{'upstream'} = undef;
		$hFiltersChr->{'downstream'} = undef;
	}
	else { $hFiltersChr->{$filter_name} = undef; }
}

if ($filter_type_variation_2) {
	$filter_type_variation_2 =~ s/ /,/g;
	foreach my $filter_name (split(',', $filter_type_variation_2)) {
		if ($filter_name eq 'cnv') {
			$hFiltersChr_var2->{'large_deletion'} = undef;
			$hFiltersChr_var2->{'large_duplication'} = undef;
		}
		elsif ($filter_name eq 'upstream_downstream') {
			$hFiltersChr_var2->{'upstream'} = undef;
			$hFiltersChr_var2->{'downstream'} = undef;
		}
		else { $hFiltersChr_var2->{$filter_name} = undef; }
	}
}

if ($model eq 'recessif' or $model eq 'compound' or $model eq 'recessif_compound') {
	$hFiltersChr->{intergenic} = undef;
	if ($filter_type_variation_2) {
		$hFiltersChr_var2->{intergenic} = undef;
	}
}

if ($infos) {
	print "\n\n######### INFOS CHR$filter_chromosome #########\n\n";
	
	my $chr = $project->getChromosome($filter_chromosome);
	warn ref($chr);
	print "\n\n";
	print "### CHR".$chr->id()."\n\n";
	foreach my $cat_name (sort keys %{$chr->global_categories()}) {
		print "   global_categories - $cat_name: ".$chr->countThisVariants( $chr->global_categories->{$cat_name} )." variants (Size: ".$chr->global_categories->{$cat_name}->Size().")\n";
	}
	warn "\n";
	foreach my $patient (@{$chr->getPatients()}) {
		print "# Patient ".$patient->name().' ('.ref($patient)." ):\n";
		print "   all: ".$chr->countThisVariants( $patient->getVariantsVector($chr) )." variants (Size: ".$patient->getVariantsVector($chr)->Size().")\n";
		print "   he : ".$chr->countThisVariants( $patient->getHe($chr) )." variants (Size: ".$patient->getHe($chr)->Size().")\n";
		print "   ho : ".$chr->countThisVariants( $patient->getHo($chr) )." variants (Size: ".$patient->getHo($chr)->Size().")\n";
		print "\n";
	}
	print "\n\n";
	
	my @lVar = @{$chr->getListVarObjects( $chr->global_categories->{'large_duplication'} )};
	foreach my $v (@lVar) {
		warn ref($v).' -> '.$v->id();
		warn Dumper $v->annotation();
		foreach my $g (@{$v->getGenes()}) {
			warn $g->external_name().': '.$v->variationType($g);
		}
		exit(0);
	}
	
	
	exit(0);
	foreach my $gene (@{$chr->getGenes()}) {
		print "# Gene ".$gene->external_name().' ('.ref($gene)." ):\n";
		foreach my $cat_name (sort keys %{$gene->categories()}) {
			print "   categories - $cat_name: ".$chr->countThisVariants( $gene->categories->{$cat_name} )." variants (Size: ".$gene->categories->{$cat_name}->Size().")\n";
			die if ($cat_name eq 'cnv');
		}
		print "\n";
	}
	print "\n\n";
	exit(0);
}

unless ($filter_chromosome) {
	unless ($export_vcf_for or $xls_by_variants or $xls_by_genes or $xls_load_session) { 
		warn "\n\nERROR: -filter_chromosome option mandatory ! Die...\n\n";
		die;
	}
}
if ($filter_chromosome eq 'all') {
	$filter_chromosome = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT';
}
$project->set_filter_chromosome($filter_chromosome, ',') if ($filter_chromosome);
$project->level_ind('gene') if ($level_ind eq 'gene');			
$project->level_fam('gene') if ($level_fam eq 'gene' and $typeFilters eq 'familial');

my $hResumeFilters;
if ($xls_by_variants or $xls_by_genes or $xls_load_session) {
	$hResumeFilters->{typeFilters} = $typeFilters;
	$hResumeFilters->{model} = $model;
	@{$hResumeFilters->{filter_type_variation}} = split(',', $filter_type_variation);
	$hResumeFilters->{dejavu} = $dejavu;
	$hResumeFilters->{dejavu} .= ' only HO' if ($dejavu_ho);
	$hResumeFilters->{atLeast} = $atLeast;
	@{$hResumeFilters->{filter_region}} = split(',', $filter_region);
	@{$hResumeFilters->{filter_attic}} = split(',', $filter_attic);
	@{$hResumeFilters->{filter_patient}} = split(',', $filter_patient);
	@{$hResumeFilters->{filter_not_patient}} = split(',', $filter_not_patient);
	@{$hResumeFilters->{filter_he}} = split(',', $filter_he);
	@{$hResumeFilters->{filter_ho}} = split(',', $filter_ho);
	@{$hResumeFilters->{fam_not}} = split(',', $fam_not);
	@{$hResumeFilters->{fam_and}} = split(',', $fam_and);
}

if ($xls_load_session) {
	loadSessionsXLS($project, $xls_load_session);
}

if ($xls_save_session) {
	print $cgi->header('text/json-comment-filtered');
	print "{\"progress\":\".";
}
else {
	if ($export_vcf_for) {}
	elsif ($detail_project) {}
	elsif ($xls_by_regions_ho) {}
	else {
		print $cgi->header('text/json-comment-filtered');
		print "{\"progress\":\".";
	}
}


########## [BEGIN] FILTRES DE POLYQUERY PAR CHROMOSOME ########## 


my $hAllIds;
my (@lVarObj, @lHashStatsGenes, $h_all_genes_name);
my $nb_chr = 0;
foreach my $chr_id (split(',', $filter_chromosome)) {
	my $chr = $project->getChromosome($chr_id);
	$nb_chr++;
	if ($chr->not_used()) {
		print "@" unless ($export_vcf_for or $detail_project or $xls_by_regions_ho);
		next;
	}
	
	if ($debug) { warn "\n\nCHR ".$chr->id()." -> INIT - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	
	checkEssentialsCategoriesInChr($chr);
	
	#TODO: a enlever
#	my $vector_global_fam_denovo = $chr->getNewVector();
#	$chr->save_model_variants_all_patients('init');
#	$vector_global_fam_denovo += $chr->getModelVector_fam_denovo();
#	$chr->save_model_variants_all_patients('denovo');
#	$chr->load_init_variants_all_patients('init');
	
	
	if ($export_vcf_for) {}
	elsif ($detail_project) {}
	elsif ($xls_by_regions_ho) {}
	else { $project->cgi_object(1); }
	unless (check_region_filter($chr, $filter_region)) {
		print "@" unless ($export_vcf_for or $detail_project or $xls_by_regions_ho);
		next;
	}
	if ($xls_save_session) {
		print "@" unless ($export_vcf_for or $detail_project or $xls_by_regions_ho);
	}
	
	# Pour XLS, besoin de savoir si un patient possede ou non le variant, meme s il ne passe pas un filtre ou un modele pour rester coherent a la 2e page par gene
	if ($xls_by_variants) { $chr->save_model_variants_all_patients('for_xls'); }
	
	QueryVectorFilter::setInTheAttic($chr, $project->getPatientsFromListNames([split(' ', $filter_attic)]));
	print "@" unless ($export_vcf_for or $detail_project or $xls_by_regions_ho);
	
	QueryVectorFilter::filter_vector_ratio($chr, $filter_ratio, $filter_type_ratio);
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_ratio - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	
	QueryVectorFilter::filter_vector_ncboost($chr, $filter_ncboost);
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_ncboost - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	
#	QueryVectorFilter::filter_vector_global_gnomad_freq($chr, $filter_gnomad_test);
#	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER getVectorGnomadCategory - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	
	my $h_args;	
	$chr->save_model_variants_all_patients('init');
	doPolyQueryFilters_global_cat($chr, $hFiltersChr, $dejavu, $polyscore);
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER doPolyQueryFilters_global_cat - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	
	my $vector_filtered = $chr->getVariantsVector->Clone();
	my $vector_filtered_2;
	if ($hFiltersChr and $hFiltersChr_var2) {
		$chr->load_init_variants_all_patients('init');
		doPolyQueryFilters_global_cat($chr, $hFiltersChr_var2, $dejavu_2) if ($hFiltersChr_var2 or $dejavu_2 or $dejavu_ho or $filter_nbvar_regionho);
		$vector_filtered_2 = $chr->getVariantsVector->Clone();
		$h_args->{'filters_1'} = $hFiltersChr;
		$h_args->{'filters_2'} = $hFiltersChr_var2;
		$h_args->{'vector_filters_1'} = $vector_filtered;
		$h_args->{'vector_filters_2'} = $vector_filtered_2;
		if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER doPolyQueryFilters_global_cat 2 - nb Var: ".$chr->countThisVariants($vector_filtered_2); }
	}
	
	QueryVectorFilter::filter_vector_gnomad_ac($chr, $filter_gnomad) if ($filter_gnomad);
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_gnomad_ac - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	
#	QueryVectorFilter::filter_vector_ratio_allele($chr, $filter_ratio) if ($filter_ratio);
#	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_ratio_allele - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	
	if ($panel_name) {
		QueryVectorFilter::filter_genes_from_ids($chr, $hChr->{$chr->id()}, $can_use_hgmd);
		if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_genes_from_ids - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	}
	else {
		QueryVectorFilter::filter_usefull_genes_ids($chr, $hFiltersChr, $can_use_hgmd);
		if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_usefull_genes_ids - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	}
	
	# FILTRE des genes (annotations, search, selection)
	QueryVectorFilter::filter_genes_text_search($chr, $filter_text);
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_genes_text_search - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	
	QueryVectorFilter::filter_genes_only_genes_names($chr, $only_genes);
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_genes_only_genes_names - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	
	if ($hFiltersChr and $hFiltersChr_var2) {
		QueryVectorFilter::filter_genes_annotations($chr, $hFiltersChr_var2);
	}
	else {
		QueryVectorFilter::filter_genes_annotations($chr, $hFiltersChr);
	}
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_genes_annotations - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	
	# FILTRE des modeles genetiques
	if ($model and $project->typeFilters() eq 'individual') {
		QueryVectorFilter::filter_model_individual_recessif($chr) if ($model eq 'recessif');
		QueryVectorFilter::filter_model_individual_compound($chr) if ($model eq 'compound');
	}
	elsif ($model and $project->typeFilters() eq 'familial') {
		my @lModels;
		if ($model eq 'recessif_compound') {
			push(@lModels, 'recessif');
			push(@lModels, 'compound');
		}
		elsif ($model eq 'uniparental_recessif_compound') {
			push(@lModels, 'recessif');
			push(@lModels, 'compound');
			push(@lModels, 'uniparental_disomy');
		}
		else { push(@lModels, $model); }
		QueryVectorFilter::filter_models_familial_union($chr, \@lModels, $h_args);
	}
	elsif ($model and $project->typeFilters() eq 'somatic') {
		QueryVectorFilter::filter_model_somatic_loh($chr) if ($model eq 'loh');
		QueryVectorFilter::filter_model_somatic_dbl_evt($chr) if ($model eq 'dbl_evt');
		QueryVectorFilter::filter_model_somatic_only_tissues_somatic($chr) if ($model eq 'only_tissues_somatic');
	}
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER models - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	
	if ($project->filter_text()) {
		foreach my $patient (@{$chr->getPatients()}) {
			$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $chr->getVariantsVector());
		}
		foreach my $family (@{$chr->getFamilies()}) {
			$family->getVariantsVector($chr)->Intersection($family->getVariantsVector($chr), $chr->getVariantsVector());
		}
	}
	
	print "@" unless ($export_vcf_for or $detail_project or $xls_by_regions_ho);
	launch_filters_region($chr, $filter_region, 1);
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER launch_filters_region - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	
	if ($chr->getVariantsVector->is_empty()) {
		print "@" unless ($export_vcf_for or $detail_project or $xls_by_regions_ho);
		next;
	}
	print "@" unless ($export_vcf_for or $detail_project or $xls_by_regions_ho);
	
	unless ($filter_nbvar_regionho > 0) {
		QueryVectorFilter::setExcludePatient($chr, $project->getPatientsFromListNames([split(' ', $filter_he)]) , 'he');
		QueryVectorFilter::setExcludePatient($chr, $project->getPatientsFromListNames([split(' ', $filter_ho)]), 'ho');
		QueryVectorFilter::setExcludePatient($chr, $project->getPatientsFromListNames([split(' ', $filter_not_patient)]), 'all');
		QueryVectorFilter::setExcludeFamily($chr, $project->getFamiliesFromListNames([split(' ', $fam_not)]));
	}
	
	unless ($filter_nbvar_regionho == 0) {
		QueryVectorFilter::setIntersectPatient_HO_REGIONS($chr, $project->getPatientsFromListNames([split(' ', $filter_patient)]), $filter_nbvar_regionho);
		QueryVectorFilter::setIntersectFamily_REC_REGIONS($chr, $project->getFamiliesFromListNames([split(' ', $fam_and)]), $filter_nbvar_regionho);
		QueryVectorFilter::setExcludePatient_HO_REGIONS($chr, $project->getPatientsFromListNames([split(' ', $filter_not_patient)]), $filter_nbvar_regionho);
		QueryVectorFilter::setExcludeFamily_HO_REGIONS($chr, $project->getFamiliesFromListNames([split(' ', $fam_not)]), $filter_nbvar_regionho);
	}
	
	QueryVectorFilter::setIntersectPatient($chr, $project->getPatientsFromListNames([split(' ', $filter_patient)]));
	QueryVectorFilter::setIntersectPatient_GENES($chr, $project->getPatientsFromListNames([split(' ', $filter_patient)]));
	QueryVectorFilter::setIntersectFamily($chr, $project->getFamiliesFromListNames([split(' ', $fam_and)]));
	QueryVectorFilter::setIntersectFamily_GENES($chr, $project->getFamiliesFromListNames([split(' ', $fam_and)]));
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER exclude / intersect - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	
	if ($project->typeFilters() eq 'individual') {
		QueryVectorFilter::filter_atLeast($chr, $atLeast, $project->typeFilters(), $level_ind);
	}
	else {
		QueryVectorFilter::filter_atLeast($chr, $atLeast, $project->typeFilters(), $level_fam);
	}
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_atLeast - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	
	# DELETE FILTRE genetic models independant des genes (For PolyDiag)
#	QueryVectorFilter::delete_filter_model_familial_denovo($chr) 			if (exists $hDeleteModels->{'fam_denovo'});
#	QueryVectorFilter::delete_filter_model_familial_strict_denovo($chr) 	if (exists $hDeleteModels->{'fam_strict_denovo'});
#	QueryVectorFilter::delete_filter_model_familial_dominant($chr) 			if (exists $hDeleteModels->{'fam_dominant'});
#	QueryVectorFilter::delete_filter_model_familial_mosaic($chr) 			if (exists $hDeleteModels->{'fam_mosaic'});
#	QueryVectorFilter::delete_filter_model_familial_recessif($chr) 			if (exists $hDeleteModels->{'fam_recessif'});
#	QueryVectorFilter::delete_filter_model_familial_uniparental_disomy($chr)if (exists $hDeleteModels->{'fam_uniparental_disomy'});
#	QueryVectorFilter::delete_filter_model_familial_both_parents($chr)		if (exists $hDeleteModels->{'fam_both_parents'});
	
	launch_filters_region($chr, $filter_region);
	if ($debug) { warn "\nAfter launch_filters_region"; }
	
	check_variants_regions_exclude($chr);
	if ($debug) { warn "\nAfter check_variants_regions_exclude"; }
	
	if ($export_vcf_for) {
		foreach my $pat_name (split(' ', $export_vcf_for)) {
			my $patient = $chr->getPatient( $pat_name );
			foreach my $var (@{$chr->getListVarObjects($patient->getVariantsVector($chr))}) {
				my ($chr_name, $pos_vcf, $ref_vcf, $var_vcf) = split('_', $var->{check_id});
				$hAllIds->{$pat_name}->{$chr->id()}->{$pos_vcf}->{id} = $var->id();
				$hAllIds->{$pat_name}->{$chr->id()}->{$pos_vcf}->{ref_vcf} = $ref_vcf;
				$hAllIds->{$pat_name}->{$chr->id()}->{$pos_vcf}->{var_vcf} = $var_vcf;
			}
		}
	}
	elsif ($export_list_var_ids){
		push( @lVarObj, @{$chr->getStructuralVariations()} );
		if ($debug) { warn "\nAfter chr->getStructuralVariations()"; }
	}
	else {
		print "@" unless ($export_vcf_for or $detail_project or $xls_by_regions_ho);
		launchStatsGlobalChr($chr) unless ($without_stats);
		if ($debug) { warn "\nAfter launchStatsGlobalChr()"; }
	}
	print "@" unless ($export_vcf_for or $detail_project or $xls_by_regions_ho);
	if ($debug) {
		warn "\nCHR ".$chr->id()." -> nb Var: ".$chr->countThisVariants($chr->getVariantsVector());
	}
}


if ($export_list_var_ids) {
	print "\n";
	foreach my $var (@lVarObj) {
		my $start = $var->start();
		my $end = $var->end();
		my $ref = $var->ref_allele();
		my $alt = $var->var_allele();
		my $varid = $var->id();
		my @lTmp = split('_', $varid);
		my $chr_name = $lTmp[0];
		my $rsname = $var->rs_name();
		my $has_gene = 0;
		my $string;
		foreach my $gene ( @{ $var->getGenes() }) {
			$has_gene = 1;
			my $g_name = $gene->external_name();
			my $annot_gene = $var->variationType($gene);
			my $annot_trans;
			my @lTr;
			foreach my $tr (@{$gene->getTranscripts()}) {
				my $t_id = $tr->id();
				my $annot_t = $var->variationType($tr);
				push(@lTr, $t_id.','.$annot_t);
			}
			$string = "$chr_name\t$start\t$end\t$ref\t$alt\t$varid\t$rsname\tg:$g_name,$annot_gene\tt:".join(';', @lTr);
			print "$string\n";
		}
		if ($has_gene == 0) {
			$string = "$chr_name\t$start\t$end\t$ref\t$alt\t$varid\t$rsname\tg:intergenic,intergenic";
			print "$string\n";
		}
	}
	exit;
}

if ($export_vcf_for) { export_vcf($project, $export_vcf_for, $hAllIds); }
my $hashRes;
$hashRes->{'label'} = 'name';
$hashRes->{'items'} = launchStatsProject();

if ($gene_atlas_view) {
	geneAtlasView($project);
}

if    ($xls_save_session) { saveSessionXLS($project, $hashRes, $hResumeFilters); }
elsif ($xls_by_variants)  { getXls_byVar($project, undef, $hResumeFilters, $xls_outfile); }
elsif ($xls_by_genes)     {
	my @lStats = @lHashStatsGenes;
	my $h;
	$h->{by_genes} = launchStatsProjectAll_genes();
	getXls_byGenes($project, $h, $hResumeFilters, $xls_outfile);
}
else {
	$project->{lmdbOmim}->close() if (exists $project->{lmdbOmim} and $project->{lmdbOmim});
	if ($debug) {
		warn "\n\nEND of debug\n\n";
		exit;
	}
	printJson($hashRes, $test);
} 


########## [END] FILTRES DE POLYQUERY PAR CHROMOSOME ##########


sub launch_filters_region {
	my ($chr, $filter_region, $first_launch) = @_;
	return unless ($filter_region);
	my @lFilters;
	my $isIntersect;
	my $hFilters;
	foreach my $this_filter (split(" ", $filter_region)) {
		next if (exists $hFilters->{$this_filter});
		$hFilters->{$this_filter} = undef;
		my ($chrId, $start, $end, $include) = split(":", $this_filter);
		next unless ($chrId eq $chr->id());
		$chr->check_each_var_filter_region($this_filter, $first_launch);
		$isIntersect = 1 if ($include eq '0');
	}
	$chr->variants_regions_add();
	if ($chr->variants_regions_add->is_empty()) {
		if ($isIntersect) {
			foreach my $patient (@{$chr->getPatients()}) {
				$patient->getVariantsVector($chr)->Empty();
			}
			$chr->update_from_patients();
			return;
		}
	}
	if ($first_launch) {
		print "." unless ($export_vcf_for or $detail_project);
		$chr->{variants} = $chr->variants_regions_add() if ($chr->variants_regions_add());
		$chr->{variants} -= $chr->variants_regions_exclude() if ($chr->variants_regions_exclude());
		foreach my $patient (@{$chr->getPatients()}) {
			$patient->getVariantsVector($chr)->Intersection( $patient->getVariantsVector($chr), $chr->getVariantsVector() );
		}
	}
}

sub check_variants_regions_exclude() {
	my $chr = shift;
	if ($chr->variants_regions_exclude()) {
		$chr->{variants} -= $chr->variants_regions_exclude();
		foreach my $patient (@{$chr->getPatients()}) {
			$patient->getVariantsVector($chr)->Intersection( $patient->getVariantsVector($chr), $chr->getVariantsVector() );
		}
		print "." unless ($export_vcf_for or $detail_project);
	}
}

###### STATS METHODS #####



sub launchStatsProject {
	return launchStatsProjectXlsGenes() if ($project->get_xls() eq 'genes');
	return launchStatsProjectXlsVariants() if ($project->get_xls() eq 'variants');
	return launchStatsChromosomes() if ($without_stats);
	return launchStatsProjectAll();
}

sub launchStatsProjectXlsGenes {
	my $hash_stats;
	$hash_stats->{genes} = \@lHashStatsGenes;
	return $hash_stats;
}

sub launchStatsProjectXlsVariants {
	my $hash_stats;
	$hash_stats->{patients}	= launchStatsProjectAll_patients();
	return $hash_stats;
}

sub launchStatsChromosomes {
	my $hash_stats;
	warn "\n\n" if ($debug);
	warn "\n# launchStatsProjectAll_chromosomes" if ($debug);
	$hash_stats->{chromosomes} 	= launchStatsProjectAll_chromosomes();
	foreach my $chr_id (split(',', $filter_chromosome)) {
		my $chr = $project->getChromosome($chr_id);
		purge($chr);
		next if ($chr->not_used());
	}
	warn "\n# END" if ($debug);
	return $hash_stats;
}

sub launchStatsProjectAll {
	my $hash_stats;
	warn "\n\n" if ($debug);
	warn "\n# launchStatsProjectAll_chromosomes" if ($debug);
	$hash_stats->{chromosomes} 	= launchStatsProjectAll_chromosomes();
	warn "\n# launchStatsProjectAll_genes" if ($debug);
	$hash_stats->{genes} 		= launchStatsProjectAll_genes();
	warn "\n# launchStatsProjectAll_families" if ($debug);
	$hash_stats->{families} 	= launchStatsProjectAll_families();
	warn "\n# launchStatsProjectAll_patients" if ($debug);
	$hash_stats->{patients} 	= launchStatsProjectAll_patients();
	warn "\n# stats_region" if ($debug);
	$hash_stats->{regions} 	    = $project->stats_region();
	warn "\n# launchStatsProjectAll_groups" if ($debug);
	$hash_stats->{groups}		= launchStatsProjectAll_groups() if ($project->isSomaticStudy());
	warn "\n# json_all_genes_name" if ($debug);
	$hash_stats->{genes_name}	= json_all_genes_name();
	warn "\n# stats_regions_ho_rec" if ($debug);
	$hash_stats->{regions_ho_rec} = launchStatsPojectRegionsHoRec() if ($filter_nbvar_regionho);
	return export_regions_ho_xls() if ($filter_nbvar_regionho and $xls_by_regions_ho);
	warn "\n# purge" if ($debug);
	foreach my $chr_id (split(',', $filter_chromosome)) {
		my $chr = $project->getChromosome($chr_id);
		purge($chr);
		next if ($chr->not_used());
	}
	warn "\n# END" if ($debug);
	return $hash_stats;
}

sub export_regions_ho_xls {
	my $hres;
	foreach my $chr_id (split(',', $filter_chromosome)) {
		my $chr = $project->getChromosome($chr_id);
		next if ($chr->not_used());
		foreach my $h (@{launchStastChr_regionsHoRec($chr)}) {
			warn Dumper $h;
			my $start = $h->{'start'};
			my $patient_name = $h->{'name'};
			foreach my $key (keys %$h) {
				$hres->{$chr->id()}->{$start}->{$patient_name}->{$key} = $h->{$key};
			}
		}
	}
	my $projectName = $project->name();
	my $workbook;
	if ($xls_outfile) {
		$workbook = Spreadsheet::WriteExcel->new( $xls_outfile );
	}
	else {
		print "Content-type: application/msexcel\n";
		print "Content-Disposition: attachment;filename=$projectName\_regions_ho.xls\n\n";
		$workbook = Spreadsheet::WriteExcel->new( \*STDOUT );
	}
	my $hFormat = getHashFormat($workbook);
	my $xls_page = $workbook->add_worksheet($projectName);
	my $listColNames = writeHeader_byRegionsHo($workbook, $xls_page);
	my $hLen = initHashSizeColumn_byVar($listColNames);
	my @lLines;
	foreach my $chr_id (sort keys %{$hres}) {
		foreach my $start (sort {$a <=> $b} keys %{$hres->{$chr_id}}) {
			foreach my $patient_name (sort keys %{$hres->{$chr_id}->{$start}}) {
				my $h = $hres->{$chr_id}->{$start}->{$patient_name};
				my @lCol;
				push(@lCol, $h->{'chr'});
				push(@lCol, $h->{'start'});
				push(@lCol, $h->{'end'});
				push(@lCol, $h->{'length'});
				push(@lCol, $h->{'fam'});
				push(@lCol, $h->{'name'});
				push(@lCol, $h->{'nb_var_ho_before'});
				push(@lCol, $h->{'nb_genes'});
				push(@lCol, $h->{'genes'});
				push(@lLines, \@lCol);
			}
		}
	}
	my $i = 1;
	foreach my $listCol (@lLines) {
		my $j = 0;
		foreach my $value (@$listCol) {
			$xls_page->write($i, $j, $value, $hLen->{$j});
			$j++;
		}
		$i++;
	}
	exit(0);
}

sub launchStatsPojectRegionsHoRec {
	my @list;
	foreach my $chr_id (split(',', $filter_chromosome)) {
		my $chr = $project->getChromosome($chr_id);
		next if ($chr->not_used());
		foreach my $h (@{launchStastChr_regionsHoRec($chr)}) {
			push(@list, $h);
		}
	}
	return \@list;
}

sub launchStatsProjectAll_genes {
	my @lStats;
	foreach my $chr_id (split(',', $filter_chromosome)) {
		my $chr = $project->getChromosome($chr_id);
		next if ($chr->not_used());
		foreach my $gene (@{$chr->getGenes()}) {
			warn ref($gene) if ($debug);
			$project->print_dot(250);
			$gene->getVariantsVector->Intersection($gene->getVariantsVector(), $chr->getVariantsVector());
			next if ($gene->getVariantsVector->is_empty());
			#$gene->getVariantsVector()
			next unless ($gene->getVariantsVector->subset($chr->getVariantsVector()));
			warn '  -> test $gene->external_name() - si bloque: sans doute pb /data-xfs/public-data/HG19/lmdb/annotations/ (GenBoProject::liteAnnotations -> GenBoNoSqlAnnotation) - table annotations / synonyms locked ?'  if ($debug);
			warn $gene->external_name() if ($debug);
			#next if (scalar($gene->getPatients()) == 0);
			warn '  -> begin stats' if ($debug);
			my $hStats = launchStatsGene($gene);
			warn '  -> end stats' if ($debug);
			push(@lStats, $hStats );
			$h_all_genes_name->{$gene->id()}->{external_name} = uc($gene->external_name());
			$h_all_genes_name->{$gene->id()}->{patients} = $hStats->{'patients_name'};
			my @lFam = split(',', $hStats->{'families'});
			if (scalar(@lFam) > 0) {
				foreach my $famName (@lFam) {
					$h_all_genes_name->{$gene->id()}->{families}->{$famName} = undef;
				}
			}
			else { $h_all_genes_name->{$gene->id()}->{families} = undef; }
		}
		warn '-> all genes DONE' if ($debug);
	}
	return \@lStats;
}
		
sub launchStatsProjectAll_chromosomes {
	my @lStats;
	foreach my $chr_id (split(',', $filter_chromosome)) {
		my $chr = $project->getChromosome($chr_id);
		next if ($chr->not_used());
		my $hashChr = launchStatsChr($chr);
		if ($hashChr->{'variations'} == 0) { $hashChr = launchStatsChr_null($chr); }
		else { push(@lStats, $hashChr); }
	}
	if (scalar(@lStats) == 0 and $without_stats) {
		my $chr = $project->getChromosome('1');
		my $hashChr = launchStatsChr_null($chr);
		my $vector_bin = $chr->getVariantsVector->to_Bin();
		$hashChr->{vector_HEX} = $chr->compress_vector_bin($vector_bin);
		push(@lStats, $hashChr);
	}
	elsif (scalar(@lStats) == 0) { push(@lStats, launchStatsProjectAll_chromosomes_null()); }
	return \@lStats;
}

sub launchStatsProjectAll_chromosomes_null {
	my $hash;
	$hash->{id}            = 'No Result';
	$hash->{name}          = 'No Result';
	$hash->{genes}         = 0;
	$hash->{variations}    = 0;
	$hash->{substitutions} = 0;
	$hash->{deletions}     = 0;
	$hash->{insertions}    = 0;
	$hash->{cnvs}	       = 0;
	$hash->{heterozygote}  = 0;
	$hash->{homozygote}    = 0;
	$hash->{stop}          = 0;
	$hash->{coding}        = 0;
	return $hash;
}

sub launchStatsProjectAll_families {
	my @lStats;
	my $hash_global_stats_fam;
	foreach my $chr_id (split(',', $filter_chromosome)) {
		my $chr = $project->getChromosome($chr_id);
		if ($chr->not_used()) { $chr->getVariantsVector->Empty(); }
		foreach my $family (@{$project->getFamilies()}) {
			$project->print_dot(100);
			my $famName = $family->name();
			my $hStatsFam = launchStatsFamily($family, $chr);
			my @lCat1 = ('model', 'include', 'name', 'nb', 'id');
			foreach my $category (@lCat1) {
				$hash_global_stats_fam->{$famName}->{$category} = $hStatsFam->{$category};
			}
			my @lCat2 = ('homozygote', 'heterozygote', 'substitution', 'insertion', 'deletion', 'genes', 'cnv');
			foreach my $category (@lCat2) {
				unless (exists $hash_global_stats_fam->{$famName}->{$category}) { $hash_global_stats_fam->{$famName}->{$category} = 0; }
				unless ($chr->not_used()) { $hash_global_stats_fam->{$famName}->{$category} += $hStatsFam->{$category}; }
			}
		}
	}
	foreach my $famName (sort keys %{$hash_global_stats_fam}) {
		push(@lStats, $hash_global_stats_fam->{$famName});
	}
	return \@lStats;
}

sub launchStatsProjectAll_patients {
	my @lStats;
	my $hash_global_stats_pat;
	foreach my $chr_id (split(',', $filter_chromosome)) {
		my $chr = $project->getChromosome($chr_id);
		foreach my $patient (@{$project->getPatients}) {
			my $patName = $patient->name();
			next if (not exists $chr->patients_categories->{$patient->name()} and int($project->annotation_version()) < int($project->annotation_version_current()));
			$project->print_dot(1);
			my $hStatsPat = launchStatsPatient($patient, $chr);
			my @lCat1 = ('fam', 'child', 'name', 'id', 'group', 'coverage', 'status', '1x', '5x', '15x', '30x', 'include', 'sex', 'filter_heho', 'tissue', 'bam');
			foreach my $category (@lCat1) {
				$hash_global_stats_pat->{$patName}->{$category} = $hStatsPat->{$category};
			}
			next if ($chr->not_used()); 
			my @lCat2 = ('homozygote', 'heterozygote', 'substitutions', 'insertions', 'deletions', 'cnvs', 'genes', 'variations', 'composite');
			foreach my $category (@lCat2) {
				$hash_global_stats_pat->{$patName}->{$category} += $hStatsPat->{$category};
			}
		}
		foreach my $patient (@{$project->getPatients}) {
			my $patName = $patient->name();
			my @lCat = ('homozygote', 'heterozygote', 'variations', 'substitutions', 'insertions', 'deletions', 'genes', 'composite', 'cnvs');
			foreach my $category (@lCat) {
				unless (exists $hash_global_stats_pat->{$patName}->{$category}) {
					$hash_global_stats_pat->{$patName}->{$category} = 0;
				}
			}
			$hash_global_stats_pat->{$patName}->{'composite'} = 0;
		}
	}
	foreach my $patName (sort keys %{$hash_global_stats_pat}) {
		push(@lStats, $hash_global_stats_pat->{$patName});
	}
	return \@lStats;
	
}

sub launchStatsProjectAll_groups {
	my @lStats;
	my $hash_global_stats_groups;
	foreach my $chr (@{$project->getChromosomes()}) {
		next if ($chr->not_used());
		my $hashGroups = $chr->stats_groups();
		foreach my $famName (keys %{$hashGroups}) {
			my @lCat1 = ('model', 'include', 'name', 'nb', 'id');
			foreach my $category (@lCat1) {
				$hash_global_stats_groups->{$famName}->{$category} = $hashGroups->{$famName}->{$category};
			}
			my @lCat2 = ('homozygote', 'heterozygote', 'substitution', 'insertion', 'deletion');#, 'genes');
			foreach my $category (@lCat2) {
				$hash_global_stats_groups->{$famName}->{$category} += $hashGroups->{$famName}->{$category};
			}
		}
	}
	foreach my $groupName (sort keys %{$hash_global_stats_groups}) {
		push(@lStats, $hash_global_stats_groups->{$groupName});
	}
	return \@lStats;
}

sub launchStatsGlobalChr {
	my $chr = shift;
	if ($project->get_xls() eq 'variants') {
		my $hash_global;
		foreach my $patient (@{$chr->getPatients()}) {
			$project->print_dot(1);
			$hash_global = store_var_ids($chr, $patient, $hash_global);
		}
	}
	else {
		my @lIds = @{$chr->getListVarVectorIds($chr->getVariantsVector())};
		$chr->variations_type_tree();
		$chr->variations_patients_tree();
		$chr->variations_genes_tree();
		$project->print_dot(1);
		$chr->getFamilies();
		$project->print_dot(1);
	}
}

sub launchStastChr_regionsHoRec {
	my $chr = shift; 
	return unless ($filter_nbvar_regionho);
	my (@lRegions, $hRegionsDone, $hRegionsByPos);
	if ($project->typeFilters() eq 'individual') {
		my (@lPatIntersect, @lPatExcluded, $vector_intersect, $vector_excluded);
		foreach my $patient (@{$chr->getPatients()}) {
			if ($patient->intersected()) { push(@lPatIntersect, $patient); }
			if ($patient->excluded() and $patient->excluded() eq 'ho_reg') { push(@lPatExcluded, $patient); }
		}
		if (scalar @lPatIntersect > 0) {
			$vector_intersect = $chr->getVector_ho_regions_after_intersect(\@lPatIntersect, $filter_nbvar_regionho);
			return if ($vector_intersect->is_empty());
		}
		if (scalar @lPatExcluded > 0)  { $vector_excluded = $chr->getVector_ho_regions_after_exclude(\@lPatExcluded, $filter_nbvar_regionho); }
		foreach my $patient (@{$chr->getPatients()}) {
			my $vector_region_ho = $patient->getRegionHo($chr, $filter_nbvar_regionho, $filter_regionho_sub_only)->Clone();
			if ($vector_intersect) { $vector_region_ho->Intersection($vector_region_ho, $vector_intersect); }
			if ($vector_excluded)  { $vector_region_ho -= $vector_excluded; }
			next if ($vector_region_ho->is_empty());
			my @lReg = split(',', $vector_region_ho->to_Enum());
			foreach my $region (@lReg) {
				my $vector_region = Bit::Vector->new_Enum($chr->getVariantsVector->Size(), $region);
				next if ($vector_region->is_empty());
				my ($v_id_start, $v_id_end) = split('-', $region);
				my ($start, $end);
				if ($v_id_start == 0) { $start = $chr->getVarObject(($v_id_start))->start(); }
				else { $start = $chr->getVarObject($v_id_start - 1)->start(); }
				if ($v_id_end) {
					if ($v_id_end == ($chr->size_vector() - 1)) { $end = $chr->getVarObject($v_id_end)->end(); }
					else { $end = $chr->getVarObject(($v_id_end + 1))->end(); }
				}
				else {
					if ($v_id_start == ($chr->size_vector() - 1)) { $end = $chr->getVarObject($v_id_start)->end(); }
					else { $end = $chr->getVarObject(($v_id_start + 1))->end(); }
				}
				my $region_id = $patient->name().';chr'.$chr->id().':'.$start.'-'.$end;
				next if (exists $hRegionsDone->{$region_id});
				my $length = $end - $start + 1;
				my $hThisRegion;
				$hThisRegion->{id} = $region_id;
				$hThisRegion->{fam} = $patient->getFamily->name();
				$hThisRegion->{name} = $patient->name();
				$hThisRegion->{status} = 1;
				$hThisRegion->{status} = 2 if ($patient->isIll());
				$hThisRegion->{chr} = 'chr'.$chr->id();
				$hThisRegion->{start} = $start;
				$hThisRegion->{end} = $end;
				$hThisRegion->{length} = $length;
				$hThisRegion->{nb_genes} = 0;
				$hThisRegion->{has_gene_diag} = 0;
				$hThisRegion->{is_omim} = 0;
				my @lGenes;
				foreach my $gene (@{$chr->getGenes()}) {
					my $vector_gene = $gene->getVariantsVector->Clone();
					$vector_gene->Intersection($vector_gene, $vector_region);
					unless ($vector_gene->is_empty()) {
						foreach my $t_id (@{$gene->transcripts()}) {
							if (exists $project->buffer->getHashTransIdWithCaptureDiag->{$t_id}) {
								$hThisRegion->{has_gene_diag}++;
								last;
							}
						}
						if ($gene->omim_id()) { $hThisRegion->{is_omim}++; }
						push(@lGenes, $gene->external_name());
						$hThisRegion->{nb_genes}++;
					}
				}
				$hThisRegion->{genes} = join(',', @lGenes);
				$vector_region->Intersection($vector_region, $chr->saved_model->{init}->{$patient->name().'_ho'});
				$hThisRegion->{nb_var_ho_before} = $chr->countThisVariants($vector_region);
				$vector_region->Intersection($vector_region, $patient->getHo($chr));
				$hThisRegion->{nb_var_ho_after} = $chr->countThisVariants($vector_region);
				$hRegionsDone->{$region_id} = Bit::Vector->new_Enum($chr->getVariantsVector->Size(), $region);
				$hRegionsByPos->{$start}->{$patient->name()} = $hThisRegion;
				#push (@lRegions, $hThisRegion);
			}
		}
	}
	elsif ($project->typeFilters() eq 'familial') {
		my (@lFamIntersect, @lFamExcluded, $vector_intersect, $vector_excluded);
		foreach my $family (@{$chr->getFamilies()}) {
			if ($family->intersected()) { push(@lFamIntersect, $family); }
			if ($family->excluded() and $family->excluded() eq 'rec_reg') { push(@lFamExcluded, $family); }
		}
		if (scalar @lFamIntersect > 0) {
			$vector_intersect = $chr->getVector_rec_regions_after_intersect(\@lFamIntersect, $filter_nbvar_regionho);
			return if ($vector_intersect->is_empty());
		}
		if (scalar @lFamExcluded > 0)  { $vector_excluded = $chr->getVector_rec_regions_after_exclude(\@lFamExcluded, $filter_nbvar_regionho); }
		foreach my $family (@{$chr->getFamilies()}) {
			my $vector_region_rec = $family->getVectorRegionRec($chr, $filter_nbvar_regionho, $filter_regionho_sub_only);
			if ($vector_intersect) { $vector_region_rec->Intersection($vector_region_rec, $vector_intersect); }
			if ($vector_excluded)  { $vector_region_rec -= $vector_excluded; }
			next if ($vector_region_rec->is_empty());
			my @lReg = split(',', $vector_region_rec->to_Enum());
			foreach my $region (@lReg) {
				my $vector_region = Bit::Vector->new_Enum($chr->getVariantsVector->Size(), $region);
				next if ($vector_region->is_empty());
				my ($v_id_start, $v_id_end) = split('-', $region);
				my ($start, $end);
				if ($v_id_start == 0) { $start = $chr->getVarObject(($v_id_start))->start(); }
				else { $start = $chr->getVarObject($v_id_start -1 )->start(); }
				if ($v_id_end) {
					if ($v_id_end == ($chr->size_vector() - 1)) { $end = $chr->getVarObject($v_id_end)->end(); }
					else { $end = $chr->getVarObject(($v_id_end + 1))->end(); }
				}
				else {
					if ($v_id_start == ($chr->size_vector() - 1)) { $end = $chr->getVarObject($v_id_start)->end(); }
					else { $end = $chr->getVarObject(($v_id_start + 1))->end(); }
				}
				my $region_id = $family->name().';chr'.$chr->id().':'.$start.'-'.$end;
				my $length = $end - $start + 1;
				my $hThisRegion;
				$hThisRegion->{id} = $region_id;
				$hThisRegion->{fam} = $family->name();
				$hThisRegion->{chr} = 'chr'.$chr->id();
				$hThisRegion->{start} = $start;
				$hThisRegion->{end} = $end;
				$hThisRegion->{length} = $length;
				$hThisRegion->{nb_genes} = 0;
				$hThisRegion->{has_gene_diag} = 0;
				$hThisRegion->{is_omim} = 0;
				my @lGenes;
				foreach my $gene (@{$chr->getGenes()}) {
					my $vector_gene = $gene->getVariantsVector->Clone();
					$vector_gene->Intersection($vector_gene, $vector_region);
					unless ($vector_gene->is_empty()) {
						foreach my $t_id (@{$gene->transcripts()}) {
							if (exists $project->buffer->getHashTransIdWithCaptureDiag->{$t_id}) {
								$hThisRegion->{has_gene_diag}++;
								last;
							}
						}
						foreach my $t_id (@{$gene->transcripts()}) {
							if (exists $project->buffer->getOmimTranscriptsNames->{$t_id}) {
								$hThisRegion->{is_omim}++;
								last;
							}
						}
						push(@lGenes, $gene->external_name());
						$hThisRegion->{nb_genes}++;
					}
				}
				$hThisRegion->{genes} = join(',', @lGenes);
				my $vector_fam_ho_init = $chr->getNewVector();
				foreach my $patient (@{$family->getPatients()}) {
					$vector_fam_ho_init += $chr->saved_model->{init}->{$patient->name().'_ho'}
				}
				$vector_region->Intersection($vector_region, $vector_fam_ho_init);
				$hThisRegion->{nb_var_ho_before} = $chr->countThisVariants($vector_region);
				$vector_region->Intersection($vector_region, $family->getHo($chr));
				$hThisRegion->{nb_var_ho_after} = $chr->countThisVariants($vector_region);
				$hRegionsByPos->{$start}->{$family->name()} = $hThisRegion;
				$hRegionsDone->{$region_id} = Bit::Vector->new_Enum($chr->getVariantsVector->Size(), $region);
				#push (@lRegions, $hThisRegion);
			}
		}
	}
	
	foreach my $start (sort {$a <=> $b} keys %$hRegionsByPos) {
		foreach my $pat_fam_name (sort keys %{$hRegionsByPos->{$start}}) {
			my @lCommonRegions;
			my $hThisRegion = $hRegionsByPos->{$start}->{$pat_fam_name};
			foreach my $region_id (keys %{$hRegionsDone}) {
				my $v = $hRegionsDone->{$region_id}->Clone();
				$v->Intersection($v, $hRegionsDone->{$hThisRegion->{id}});
				unless ($v->is_empty()) { push(@lCommonRegions, $region_id); }
			}
			$hThisRegion->{common_regions} = join(',', @lCommonRegions);
			push (@lRegions, $hThisRegion);
		}
	}
	return \@lRegions;
}

#TODO: TEST
sub launchStatsChr_TEST {
	my $chr = shift;
	my $hash;
	my $name = $chr->id();
	if ($name eq 'X')     { $name = 23; }
	elsif ($name eq 'Y')  { $name = 24; }
	elsif ($name eq 'MT') { $name = 25; }
	$hash->{id}         = $chr->id();
	$hash->{name}       = int($name);
	$hash->{genes}      = $chr->getNbGenes();
	$hash->{variations}	= 0;
	if (($hash->{genes} + $chr->getNbGenesIntergenic()) == 0) {
		foreach my $cat_name (values %{$chr->stats_categories()}) { $hash->{$cat_name} = 0; }
		return $hash;
	}
	foreach my $type (@{$chr->variations_type_tree->fetch_window(0, $chr->size_vector()+1)}) {
		$hash->{variations}++;
		$hash->{$chr->stats_categories->{$type}}++;
	}
	return $hash;
}

sub launchStatsChr {
	my $chr = shift;
	my $hash;
	my $name = $chr->id();
	if ($name eq 'X')     { $name = 23; }
	elsif ($name eq 'Y')  { $name = 24; }
	elsif ($name eq 'MT') { $name = 25; }
	$hash->{id}         = $chr->id();
	$hash->{name}       = int($name);
	$hash->{genes}      = $chr->getNbGenes();
	$hash->{variations}	= 0;
	if (($hash->{genes} + $chr->getNbGenesIntergenic()) == 0) {
		foreach my $cat_name (values %{$chr->stats_categories()}) { $hash->{$cat_name} = 0; }
		return $hash;
	}
	foreach my $cat (keys %{$chr->stats_categories()}) {
		my $nb = $chr->countThisVariants( $chr->getCategoryVariantsVector($cat) );
		warn "CHR".$chr->id()." -> $cat: ".$nb if ($debug);
		$hash->{$chr->stats_categories->{$cat}} += $nb;
		$hash->{variations}	+= $nb;
	}
	my $vector_bin = $chr->getVariantsVector->to_Bin();
	$hash->{vector_HEX} = $chr->compress_vector_bin($vector_bin);
	my $bin = $chr->decompress_vector_bin($hash->{vector_HEX});
	unless($vector_bin eq $bin) {
		warn "\n\n";
		warn '[BEFORE decompress] Nb var: '.$chr->countThisVariants($chr->getVariantsVector());
		my $vector = $chr->getNewVector();
		$vector->from_Bin($bin);
		warn '[AFTER decompress] Nb var: '.$chr->countThisVariants($vector);
		confess("\n\nERROR: [chr$name] vector_bin and GenBoCache::decompress_vector_bin different... Die !\n\n");
	}
	return $hash;
}

sub launchStatsChr_null {
	my $chr = shift;
	$chr->not_used('1');
	my $hash;
	my $name = $chr->id();
	if ($name eq 'X')     { $name = 23; }
	elsif ($name eq 'Y')  { $name = 24; }
	elsif ($name eq 'MT') { $name = 25; }
	$hash->{id}            = $chr->id();
	$hash->{name}          = int($name);
	$hash->{genes}         = 0;
	$hash->{variations}    = 0;
	$hash->{substitutions} = 0;
	$hash->{deletions}     = 0;
	$hash->{insertions}    = 0;
	$hash->{cnvs}	       = 0;
	return $hash;
}

sub launchStatsGene {
	my $gene = shift;
	my $lNbPatForEachVar;
	my $hashStats;
	$hashStats->{'query_score'} = 0;
	my $gene_id = $gene->id();
	my $hashKyoto;
	if ($gene->is_intergenic()) {
		$hashStats->{'id'} = $gene->id();
	}
	else {
		$hashStats->{'id'} = $gene->id();
		warn '  -> stats $gene->transcripts()' if ($debug);
		my @lTrIds;
		foreach my $t_id (@{$gene->transcripts()}) {
			warn '    -> '.$t_id if ($debug);
			push(@lTrIds, $t_id) if (exists $project->buffer->getHashTransIdWithCaptureDiag->{$t_id}) ;
		}
		if (@lTrIds) {
			$hashStats->{'dejavu_capture_diag'} = join(',', @lTrIds);
			$hashStats->{'query_score'} += $project->buffer->config->{score_query_gene}->{dejavu_capture_diag};
		}
		warn '  -> stats $gene->omim_id()' if ($debug);
		if ($gene->omim_id()) {
			$hashStats->{'is_omim'} = $gene->omim_id();
			$hashStats->{'is_omim'} .= ';new' if ($gene->is_omim_new());
			$hashStats->{'is_omim_morbid'} = $gene->is_omim_morbid();
#			$hashStats->{'is_omim_morbid'} .= ';new' if ($gene->is_omim_morbid_new());
			$hashStats->{'query_score'} += $project->buffer->config->{score_query_gene}->{is_omim};
		}
	}
	warn '  -> stats $gene->external_name()' if ($debug);
	$hashStats->{'xref'}		= $gene->external_name();
	warn '  -> stats xref: '.$hashStats->{'xref'} if ($debug);
	$hashStats->{'name'} 		= $gene->ensg();
	warn '  -> stats xref: '.$hashStats->{'name'} if ($debug);
	$hashStats->{'start'} 		= $gene->start();
	warn '  -> stats xref: '.$hashStats->{'start'} if ($debug);
	$hashStats->{'end'} 		= $gene->end();
	warn '  -> stats xref: '.$hashStats->{'end'} if ($debug);
	warn '  -> stats chromosome: '.ref($gene->chromosome) if ($debug);
	$hashStats->{'chromosome'} 	= $gene->chromosome->id();
	warn '  -> stats chromosome_id: '.$hashStats->{'chromosome'} if ($debug);
	#$hashStats->{'description'} = $gene->phenotypes();
	$hashStats->{'phenotype'} = $gene->phenotypes();
	$hashStats->{'description'} = $gene->description();
	my @lIds = @{$gene->getIdsBitOn($gene->getVariantsVector())};
	$hashStats->{'vector_ids'}	= join(',', @lIds);
	my ($hashCount, $hGetFam, $hAllPat, @lGetPat);
	
	if ($can_use_hgmd and $gene->hgmd()) {
		$hashStats->{'has_hgmd'} = 2;
		
#		my $v_g = $gene->getVariantsVector->Clone();
#		$v_g->Intersection($v_g, $gene->getChromosome->getVariantsVector());
#		$v_g->Intersection($v_g, $gene->getChromosome->getVectorLmdbDm_newVersion());
#		unless ($v_g->is_empty()) { $hashStats->{'has_hgmd'} = 1; }

		$hashStats->{'has_hgmd'} .= ';'.$gene->external_name();
		
#		#$hashStats->{'hgmd_inheritance'} = $gene->hgmd_inheritance();
#		#$hashStats->{'hgmd_description'} = $gene->hgmd_description();
#		#$hashStats->{'hgmd_other_description'} = $gene->hgmd_other_description();
#		#$hashStats->{'hgmd_disease'} = $gene->hgmd_disease();
#		#$hashStats->{'hgmd_refseq'} = $gene->hgmd_refseq();
#		$hashStats->{'query_score'} += $project->buffer->config->{score_query_gene}->{has_hgmd_dm};
	}
	my ($hStats, $h_var_ids, $h_models);
	foreach my $id (@lIds){
		my $var_id = $gene->chromosome->getVarId($id);
		$h_var_ids->{$var_id} = undef;
		
		my $results_type    = $gene->chromosome->variations_type_tree->fetch(int($id), (int($id) + 1));
		my $results_patient = $gene->chromosome->variations_patients_tree->fetch(int($id), (int($id) + 1));
		
#		if ($gene->chromosome->project->typeFilters() eq 'familial') {
#			my $results_fam_model = $gene->chromosome->fam_genetic_model_tree->fetch(int($id), (int($id) + 1));
#			foreach my $model_name (@$results_fam_model) {
#				$h_models->{$model_name}++;
#			}
#		}
	
		my $var_type = $gene->chromosome->project->buffer->config->{'stats_genes'}->{$results_type->[0]};
		$hStats->{variants}->{all}++;
		$hStats->{variants}->{$var_type}++;
		map{$hStats->{patients}->{all}->{$_}++} @$results_patient;
		map{$hStats->{patients}->{$var_type}->{$_}++} @$results_patient;
		
		my $hPat;
		foreach my $pat_name (@{$results_patient}) {
			$hPat->{$pat_name} = undef;
			$hAllPat->{$pat_name} = undef;
		}
		
		my $results_gene = $gene->chromosome->variations_genes_tree->fetch(int($id), (int($id) + 1));
		foreach my $gene_info (@{$results_gene}) {
			my ($gene_ensg, $annot) = split(';', $gene_info);
			next unless ($gene_ensg eq $gene->id());
			next unless (exists $gene->chromosome->hash_filters_keeped->{$annot});
			next unless (exists $project->ensembl_annotations->{$annot});
			if (exists $project->ensembl_annotations->{$annot}) {
				my $impact = $gene->chromosome->project->buffer->config->{'stats_genes'}->{$annot};
				$hStats->{variants}->{$impact}++ if ($impact);
				map{$hStats->{patients}->{$impact}->{$_}++} @$results_patient;
			}
		}
		# TODO: a decocher - KEEP IF HAS HGMD DM VAR
		#if ($can_use_hgmd and $project->isUpdate() and $gene->is_HGMD_DM()) {
		if ($can_use_hgmd and $gene->is_HGMD_DM()) {
			my $results_hgmd = $gene->chromosome->variations_hgmd_dm_tree->fetch(int($id), (int($id) + 1));
			if (scalar(@$results_hgmd) > 0) {
				$hStats->{variants}->{hgmd_dm}++;
				map{$hStats->{patients}->{hgmd_dm}->{$_}++} @$results_patient;
			}
		}
	}
	
#	warn $gene->external_name().': '.join(', ', sort keys %$h_models)."\n";
	
	$hashStats->{'ids'} = join(',', keys %$h_var_ids);

	my @lGetPat;
	foreach my $pat_name (keys %$hAllPat) {
		my $text = $project->hash_patients_name->{$pat_name}->{id};
		push(@lGetPat, $text);
		my $fam_name = $project->hash_patients_name->{$pat_name}->{fam};
		$hGetFam->{$fam_name} = undef;
		$h_all_genes_name->{$gene->id()}->{patients_name}->{$pat_name} = undef;
	}

	$gene->{patients_found} = $hAllPat;
	$gene->{families_found} = $hGetFam;
	
	$hashStats->{'patients'} = join(',', sort @lGetPat);
	$hashStats->{'families'} = join(',', keys %$hGetFam);
	$hashStats->{'nb_fam'} = scalar(keys %$hGetFam);
	
	foreach my $cat ('all', 'substitution', 'indel', 'cnv', 'high', 'medium', 'low', 'hgmd_dm') {
		if (exists $hStats->{variants}->{$cat}) {
			$hashStats->{'v_'.$cat} = $hStats->{variants}->{$cat};
			$hashStats->{'p_'.$cat} = scalar keys %{$hStats->{patients}->{$cat}};
			if (exists $project->buffer->config->{score_query_gene}->{$cat}) {
				if ($project->buffer->config->{score_query_gene}->{method} eq 'cumul') {
					$hashStats->{'query_score'} += ($project->buffer->config->{score_query_gene}->{$cat} * $hashStats->{'v_'.$cat});
				}
				if ($project->buffer->config->{score_query_gene}->{method} eq 'uniq') {
					$hashStats->{'query_score'} += $project->buffer->config->{score_query_gene}->{$cat};
				}
			}
		}
		else {
			$hashStats->{'v_'.$cat} = 0;
			$hashStats->{'p_'.$cat} = 0;
		}
	}
	$hashStats->{'include'} = 1;
	$hashStats->{'include'} = 0 if (exists $project->filter_genes_intersect->{$gene->external_name()});
	
	# STATS des Regions Ho ou Rec
	if ($filter_nbvar_regionho > 1){
		$hashStats->{'region_ho'} = undef;
		$hashStats->{'region_ho_all'} = undef;
		$hashStats->{'region_rec'} = undef;
		$hashStats->{'region_rec_all'} = undef;
		my (@lObj, $region_method, $stat_method);
		if ($project->typeFilters() eq 'individual'){
			foreach my $patient (@{$gene->chromosome->getPatients()}) {
				my $vec_tmp = $patient->getRegionHo($gene->chromosome(), $filter_nbvar_regionho, $filter_regionho_sub_only)->Clone();
				$vec_tmp->Intersection($gene->getVariantsVector(), $patient->getRegionHo($gene->chromosome(), $filter_nbvar_regionho, $filter_regionho_sub_only));
				unless ($vec_tmp->is_empty()) {
					push(@{$hashStats->{'region_ho_all'}}, $project->hash_patients_name->{$patient->name()}->{id});
					$vec_tmp->Intersection($vec_tmp, $patient->getVariantsVector($gene->chromosome()));
					unless ($vec_tmp->is_empty()) {
						push(@{$hashStats->{'region_ho'}}, $project->hash_patients_name->{$patient->name()}->{id});
					}
				}
			}
		}
		elsif($project->typeFilters() eq 'familial'){
			foreach my $family (@{$gene->chromosome->getFamilies()}) {
				my $vec_tmp = $family->getVectorRegionRec($gene->chromosome(), $filter_nbvar_regionho, $filter_regionho_sub_only)->Clone();
				$vec_tmp->Intersection($gene->getVariantsVector(), $family->getVectorRegionRec($gene->chromosome(), $filter_nbvar_regionho, $filter_regionho_sub_only));
				unless ($vec_tmp->is_empty()) {
					push(@{$hashStats->{'region_rec_all'}}, $family->name());
					$vec_tmp->Intersection($vec_tmp, $family->getHo($gene->chromosome()));
					unless ($vec_tmp->is_empty()) {
						push(@{$hashStats->{'region_rec'}}, $family->name());
					}
				}
			}
			
		}
		if ($hashStats->{'region_ho'}) {
			$hashStats->{'region_ho'} = scalar @{$hashStats->{'region_ho'}}.';'.join('|' , @{$hashStats->{'region_ho'}});
		}
		if ($hashStats->{'region_ho_all'}) {
			$hashStats->{'region_ho_all'} = scalar @{$hashStats->{'region_ho_all'}}.';'.join('|' , @{$hashStats->{'region_ho_all'}});
		}
		if ($hashStats->{'region_rec'}) {
			$hashStats->{'region_rec'} = scalar @{$hashStats->{'region_rec'}}.';'.join('|' , @{$hashStats->{'region_rec'}});
		}
		if ($hashStats->{'region_rec_all'}) {
			$hashStats->{'region_rec_all'} = scalar @{$hashStats->{'region_rec_all'}}.';'.join('|' , @{$hashStats->{'region_rec_all'}});
		}
	}
	
	#TODO: Polyscore / MaxScore / ScaledScore
#	$hashStats->{'max_scaled_score'} = '-';
#	$hashStats->{'max_scaled_score'} = $gene->getMaxScaledScore() if ($gene->getProject->isUpdate());
	
#	if ($hashStats->{xref} eq 'MUC5B') {
#		die;
#	}
	
	return $hashStats;
}

sub launchStatsPatient {
	my ($patient, $chr) = @_;
	confess("\n\nERROR: GenBoPatientCache->stats() method need a GenBoChromosomeCache and a category object in argument. Die.\n\n") unless(ref($chr) eq 'GenBoChromosomeCache');
	return launchStatsPatient_null($patient) if ($chr->getVariantsVector->is_empty());
	#2018_06_28: Choix de patrick, on met artificiellement toutes les valeures a 0 si on exclu des regions ho dun patient (comme in the attic)
	return launchStatsPatient_null($patient) if ($patient->excluded() eq 'ho_reg');
	my $hStats;
	my $patName = $patient->name();
	$hStats->{name} 		= $patName;
	$hStats->{id} 			= int($project->hash_patients_name->{$patient->name()}->{'id'});
	$hStats->{bam} 			= $patient->bamUrl();
	$hStats->{fam} 			= $patient->family();
	$hStats->{sex} 			= $patient->sex();
	$hStats->{status} 		= $patient->status();
	$hStats->{group} 		= undef;
	$hStats->{filter_heho} 	= -1;
	if ($patient->excluded() eq 'he') {
		$hStats->{filter_heho} 	= 0;
	}
	elsif ($patient->excluded() eq 'ho') {
		$hStats->{filter_heho} 	= 1;
	}
	my $hPed = $project->pedigree_details();
	if ($hPed and scalar(keys %$hPed) > 0) {
		my $family = $patient->getFamily();
		if ($family) {
			if    ($family->getFather() and $family->getFather->name() eq $patient->name()) { $hStats->{child}	= 'F'; }
			elsif ($family->getMother() and $family->getMother->name() eq $patient->name()) { $hStats->{child}	= 'M'; }
			else { $hStats->{child}	= 'C'.$patient->sex(); }
		}
	}
	$hStats->{include} 	= 1;
	$hStats->{include} 	= 0 if ($patient->intersected());
	$hStats->{include} 	= 2 if ($patient->in_the_attic() == 1);
	$hStats->{include} 	= -1 if ($patient->excluded() eq 'all');
	if (exists $project->filter_attic_genes->{$patName}) {
		$hStats->{include} 	= 1;
	}
	my $hash_cov;
	if ($project->test()) {
		my $freeze_file = $project->getCoverageDir().$patient->name().'.freeze';
		$hash_cov = retrieve $freeze_file;
	}
	else { $hash_cov = $patient->coverage(); }
	$hStats->{'1x'} 		= $hash_cov->{'1x'};
	$hStats->{'5x'} 		= $hash_cov->{'5x'};
	$hStats->{'15x'} 		= $hash_cov->{'15x'};
	$hStats->{'30x'} 		= $hash_cov->{'30x'};
	$hStats->{'coverage'} 	= $hash_cov->{'mean'};
	$hStats->{all_variations} 	= 0;
	$hStats->{variations} 	= 0;
	foreach my $cat (keys %{$patient->stats_categories()}) {
		my $nb = $chr->countThisVariants( $patient->getCategoryVariantsVector($chr, $cat) );
		$hStats->{$patient->stats_categories->{$cat}} += $nb;
		$hStats->{variations} += $nb;
		$hStats->{all_variations} += $nb;
	}
	
	$hStats->{found} = '';
	if ($hStats->{all_variations} > 0) {
		$hStats->{found} = 'yes';
	}
	
	#$hStats->{genes} = $patient->countGenes($chr);
	if ($xls_by_variants or $xls_by_variants) { $hStats->{genes} = $patient->countGenes($chr); }
	else {
		my $nb = 0;
		foreach my $gene (@{$chr->getGenes()}) {
			next if ($gene->is_intergenic());
			next unless (exists $h_all_genes_name->{$gene->id()});
			$nb++ if (exists $gene->{patients_found}->{$patient->name()});
		}
		$hStats->{genes} = $nb;
	}
	$hStats->{composite} 	= 0;
	$hStats->{tissue} 		= 'C';
	if ($project->isSomaticStudy()) {
		$hStats->{group}  = $patient->somatic_group();
		$hStats->{tissue} = $patient->tissue();
	}
	return $hStats;
}

sub launchStatsPatient_null {
	my $patient = shift;
	my $hStats;
	my $patName = $patient->name();
	$hStats->{name} 		= $patName;
	$hStats->{id} 			= $patient->project->hash_patients_name->{$patient->name()}->{'id'};
	$hStats->{bam} 			= $patient->bamUrl();
	$hStats->{fam} 			= $patient->family();
	$hStats->{sex} 			= $patient->sex();
	$hStats->{status} 		= $patient->status();
	$hStats->{group} 		= undef;
	$hStats->{filter_heho} 	= -1;
	if ($patient->excluded() eq 'he') {
		$hStats->{filter_heho} 	= 0;
	}
	elsif ($patient->excluded() eq 'ho') {
		$hStats->{filter_heho} 	= 1;
	}
	my $hPed = $patient->project->pedigree_details();
	if ($hPed and scalar(keys %$hPed) > 0) {
		my $family = $patient->getFamily();
		if ($family) {
			if    ($family->getFather() and $family->getFather->name() eq $patient->name()) { $hStats->{child}	= 'F'; }
			elsif ($family->getMother() and $family->getMother->name() eq $patient->name()) { $hStats->{child}	= 'M'; }
			else { $hStats->{child}	= 'C'.$patient->sex(); }
		}
	}
	$hStats->{include} 	= 1;
	$hStats->{include} 	= 0 if ($patient->intersected());
	$hStats->{include} 	= 2 if ($patient->in_the_attic() == 1);
	$hStats->{include} 	= -1 if ($patient->excluded() eq 'all');
	$hStats->{include} 	= -1 if ($patient->excluded() eq 'ho_reg');
	my $hash_cov;
	if ($patient->project->test()) {
		my $freeze_file = $patient->project->getCoverageDir().$patient->name().'.freeze';
		$hash_cov = retrieve $freeze_file;
	}
	else { $hash_cov = $patient->coverage(); }
	$hStats->{'1x'} 		= $hash_cov->{'1x'};
	$hStats->{'5x'} 		= $hash_cov->{'5x'};
	$hStats->{'15x'} 		= $hash_cov->{'15x'};
	$hStats->{'30x'} 		= $hash_cov->{'30x'};
	$hStats->{'coverage'} 	= $hash_cov->{'mean'};
	$hStats->{variations} 	= 0;
	$hStats->{homozygote} 	= 0;
	$hStats->{heterozygote} = 0;
	$hStats->{substitutions}= 0;
	$hStats->{insertions} 	= 0;
	$hStats->{deletions} 	= 0;
	$hStats->{cnvs}		 	= 0;
	$hStats->{genes} 		= 0;
	$hStats->{composite} 	= 0;
	$hStats->{tissue} 		= 'C';
	if ($patient->project->isSomaticStudy()) {
		$hStats->{group}  = $patient->somatic_group();
		$hStats->{tissue} = $patient->tissue();
	}
	return $hStats;
}

sub launchStatsFamily {
	my ($family, $chr) = @_;
	my $typeFilters = $project->typeFilters();
	my $hash;
	my $nb_pat = 0;
	$nb_pat = scalar(@{$family->getPatients()});
	$hash->{name} 			= $family->name();
	$hash->{id} 			= $family->name();
	$hash->{nb} 			= $nb_pat;
	$hash->{model} 			= $family->used_model();
	$hash->{model} 			= undef if ($typeFilters eq 'individual' or $typeFilters eq 'somatic');
	$hash->{include} 		= 1;
	$hash->{include} 		= 0 if ($family->intersected());
	$hash->{include} 		= 2 if ($family->in_the_attic());
	$hash->{include} 		= -1 if ($family->excluded() eq 'all');
	$hash->{include} 		= -2 if ($family->excluded() eq 'rec_reg');
	if ($family->in_the_attic()) {
		my $need_initInclude = 1;
		foreach my $pat (@{$family->getPatients()}) { 
			unless (exists $project->filter_attic_genes->{$pat->name()}) {
				$need_initInclude = undef;
				last;
			}
		}
		$hash->{include} = 1 if ($need_initInclude);
	}
	foreach my $cat (keys %{$family->stats_categories()}) {
		$hash->{$family->stats_categories->{$cat}} += $chr->countThisVariants( $family->getCategoryVariantsVector($chr, $cat) );
	}
	my $nb = 0;
	foreach my $gene (@{$chr->getGenes()}) {
		next if ($gene->is_intergenic());
		next unless (exists $h_all_genes_name->{$gene->id()});
		$nb++ if (exists $h_all_genes_name->{$gene->id()}->{families}->{$family->name()});
	}
	$hash->{genes} = $nb;
	return $hash;
}



##### GLOBAL METHODS ######



sub purge {
	my $chr = shift;
	$project->{objects}->{proteins}    = {};
	$project->{objects}->{genes}       = {};
	$project->{objects}->{deletions}   = {};
	$project->{objects}->{insertions}  = {};
	$project->{objects}->{variations}  = {};
	$project->{objects}->{exons}       = {};
	$project->{objects}->{transcripts} = {};
	$chr->available_genes_ids(1);
	$chr->hash_filters_deleted(1);
	$chr->hash_freeze_file_genes(1);
	$chr->hash_freeze_file_all_genes(1);
	$chr->hash_freeze_file_patients(1);
	$chr->hash_filters_keeped(1);
	$chr->genes_object(1);
	delete $project->{objects}->{chromosomes}->{$chr->id()}->{fastGenes};
	delete $project->{objects}->{chromosomes}->{$chr->id()}->{fastTranscripts};
	$chr->cache_lmdb_variations->close();
	$chr->get_lmdb_patients->close();
	$chr->get_lmdb_categories->close();
	$chr->get_lmdb_genes->close() if ($chr->get_lmdb_genes());
	$chr->get_lmdb_variations->close();
	#$chr->get_lmdb_ncboost_chromosomes->close() if ($filter_ncboost);
	warn "\n\n" if ($debug);
	warn 'PURGE CHR'.$chr->id() if ($debug);
}

sub launch_check_genes_request {
	print $cgi->header('text/json-comment-filtered');
	print "{\"progress\":\".";
	my $hash_ok;
	foreach my $name (split(' ', $check_genes)) {
		$project->{only_genes}->{$name} = undef;
	}
	my $hashRes;
	$hashRes->{'label'} = 'name';
	# print only genes found in annotation
	my @lGenesOk;
	if ($project->only_genes_found()) { @lGenesOk = keys %{$project->only_genes_found()}; }
	$hashRes->{'items'} = \@lGenesOk;
	printJson($hashRes);
}
	
sub launch_uniq_bundle_request {
	my $hNames;
	foreach my $id (split(',', $bundle_id)) {
		foreach my $name (keys %{$buffer->getGenesNamesInBundle($id)}) {
			$hNames->{$name} = undef;
		}
	}
	my @lRes;
	foreach my $name (keys %$hNames) {
		my $h;
		$h->{name} = $name;
		$h->{in_analyse} = '';
		push(@lRes, $h);
	}
	my $hashRes;
	$hashRes->{'label'} = 'name';
	$hashRes->{'items'} = \@lRes;
	print $cgi->header('text/json-comment-filtered') unless ($test);
	printJson($hashRes);
}

sub launch_bundles_request {
	my $hGenesName;
	foreach my $gene_name (split(',', $genes_list)) {
		$hGenesName->{$gene_name} = undef;
	}
	my $hBundles;
	my $h_all_genes_in_bundles = $buffer->getAllGenesNamesInAllBundle();
	my @lRes;
	my @lPar_capture;
	my $nb_c = 0;
	foreach my $capture (sort keys %$h_all_genes_in_bundles) {
		my @lPar_bundles;
		my $hGenes_capture;
		my $nb_b = 0;
		foreach my $bundle (sort keys %{$h_all_genes_in_bundles->{$capture}}) {
			my @lPar_genes;
			my $h_bundle;
			my $nb_genes = scalar(keys %{$h_all_genes_in_bundles->{$capture}->{$bundle}->{genes}});
			my $text_nb_genes = 'genes';
			if ($nb_genes == 1) { $text_nb_genes = 'gene'; }
			$h_bundle->{type} = 'bundle';
			$h_bundle->{checked} = 'false';
			$h_bundle->{genes} = join(',', sort keys %{$h_all_genes_in_bundles->{$capture}->{$bundle}->{genes}});
			$h_bundle->{transcripts} = join(',', sort @{$h_all_genes_in_bundles->{$capture}->{$bundle}->{transcripts}});
			$h_bundle->{name} = $bundle.' <span id="span_bundle_'.$bundle.'"><font color="green">[' . $nb_genes . ' ' . $text_nb_genes . ']  </font></span><button style="padding-left=5px;" type="button" onclick="show_genes_in_capture(\'' . $bundle . '\', \'' . $h_bundle->{genes} . '\')">View Genes</button>';
			push(@lPar_bundles, $h_bundle);
			foreach my $gene (keys %{$h_all_genes_in_bundles->{$capture}->{$bundle}->{genes}}) {
				$hGenes_capture->{$capture}->{$gene} = undef;
			}
			$nb_b++;
		}
		my $h_capture;
		my $nb_captures = scalar(@lPar_bundles);
		my $text_nb_capt = 'captures';
		if ($nb_captures == 1) { $text_nb_capt = 'capture'; }
#		if ($project->isUpdate()) {
#			$h_capture->{name} = $capture.' <span id="span_capture_'.$capture.'"/><font color="green">[' . $nb_captures . ' ' . $text_nb_capt . ' - '.scalar(keys %{$hGenes_capture->{$capture}}).' genes]  </font></span><button style="padding-left=5px;" type="button" onclick="go_to_polydiag(\'' . $capture . '\')">DefiDiag View</button>';
#		}
#		else {
			$h_capture->{name} = $capture.' <span id="span_capture_'.$capture.'"/><font color="green">[' . $nb_captures . ' ' . $text_nb_capt . ' - '.scalar(keys %{$hGenes_capture->{$capture}}).' genes]  </font>';
#		}
		$h_capture->{type} = 'capture';
		$h_capture->{children} = \@lPar_bundles;
		$h_capture->{checked} = 'false';
		push(@lPar_capture, $h_capture);
		$nb_c++;
	}
	my $hashRes;
	$hashRes->{'label'} = 'name';
	$hashRes->{'items'} = \@lPar_capture;
	print $cgi->header('text/json-comment-filtered') unless ($test);
	printJson($hashRes); 
}

sub doPolyQueryFilters_global_cat {
	my ($chr, $hToFilter, $nb_dejavu, $polyscore_value) = @_;
	QueryVectorFilter::filter_vector_region_ho($chr, $filter_nbvar_regionho, $filter_regionho_sub_only, $project->typeFilters());
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_region_ho - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	my $vector_hgmd_dm = $chr->getNewVector();
	if ($can_use_hgmd and $project->isUpdate()) {
		$vector_hgmd_dm = $chr->getVectorLmdbDm->Clone();
		$vector_hgmd_dm->Intersection($vector_hgmd_dm, $chr->getVariantsVector());
	}
	QueryVectorFilter::filter_vector_type_variants($chr, $hToFilter);
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_type_variants - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	QueryVectorFilter::filter_vector_predicted_splice_region_global($chr, $hToFilter);
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_predicted_splice_region_global - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	QueryVectorFilter::filter_vector_cadd_variants($chr, $hToFilter, $keep_indels_cadd);
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_cadd_variants - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	QueryVectorFilter::filter_vector_freq_variants($chr, $hToFilter);
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_freq_variants - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	QueryVectorFilter::filter_vector_gnomad_ho_ac_variants($chr, $hToFilter);
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_gnomad_ho_ac_variants - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	QueryVectorFilter::filter_vector_confidence_variants($chr, $hToFilter);
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_confidence_variants - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
#	if ($can_use_hgmd and $project->isUpdate()) {}
	QueryVectorFilter::filter_vector_dejavu($chr, $nb_dejavu, $dejavu_ho, $test) if ($nb_dejavu);
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_dejavu - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	#QueryVectorFilter::filter_vector_polyscore($chr, $polyscore_value, $test) if ($polyscore_value);
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_vector_polyscore - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
}

sub geneAtlasView {
	my ($project) = @_;
	$project->geneAtlasView();
	exit(0);
}

sub check_region_filter {
	my ($chr, $filter_region) = @_;
	return 1 unless ($filter_region);
	my $inutil = 1;
	if ($filter_region) {
		foreach my $this_filter (split(" ", $filter_region)) {
			my ($chrId, $start, $end, $to_do) = split(":", $this_filter);
			if ($to_do eq '-1') {
				$inutil = 0;
			}
			elsif ($to_do eq '0') {
				$inutil = 0 if ($chrId eq $chr->id());
			}
		}
	}
	if ($inutil) {
		$chr->not_used(1);
		return 0;
	}
	return 1;
}



##### EXPORT #####



sub printJson {
	my ($hashRes, $test) = @_;
	my $json_encode = encode_json $hashRes;
	unless ($get_bundles) {
		print ".\",";
		$json_encode =~ s/{//;
	}
	print $json_encode;
	exit(0);
}

sub json_all_genes_name {
	my (@list, @lHash, $h);
	foreach my $gene_id (keys %{$h_all_genes_name}) {
		my $gene_name = $h_all_genes_name->{$gene_id}->{external_name};
		my $nb_pat = $h_all_genes_name->{$gene_id}->{patients};
		my $nb_fam = scalar keys %{$h_all_genes_name->{$gene_id}->{families}};
		push(@list, $gene_name.';'.$nb_pat.';'.$nb_fam);
	}
	if ($project->only_genes() and $project->only_genes_found()) {
		foreach my $gene_name (keys %{$project->only_genes()}) {
			unless (exists $project->only_genes_found->{uc($gene_name)}) {
				push(@list, $gene_name.';-1;-1');
			}
		}
	}
	$h->{name} = join(',', sort @list);
	push(@lHash, $h);
	return \@lHash;
}

sub loadSessionsXLS {
	my ($project, $sid) = @_;
	my $tmp_dir = $project->getTmpDir();
    my $session = new CGI::Session(undef, $sid, {Directory=>$tmp_dir});
    my $h = thaw(decompress $session->param('hash_xls'));
    $session->delete();
    if ($xls_by_variants) { getXls_byVar($project, $h, $hResumeFilters); }
    elsif ($xls_by_genes) { getXls_byGenes($project, $h, $hResumeFilters); }
	exit(0);
}

sub saveSessionXLS {
	my ($project, $hashRes, $hResumeFilters) = @_;
	my ($h, $ok);
	print ".";
	if ($xls_by_variants) {
		foreach my $patient (@{$project->getPatients()}) {
			print ".";
			$h->{by_pat}->{$patient->name()} = $patient->xls_var_id();
			foreach my $chr_id (keys %{$h->{by_pat}->{$patient->name()}}) {
				print ".";
				foreach my $id (keys %{$h->{by_pat}->{$patient->name()}->{$chr_id}}) {
					print ".";
					$h->{by_var}->{$chr_id}->{$id}->{$patient->name()} = undef;
					$ok = 1;
				}
			}
		}
	}
	elsif ($xls_by_genes) {
		my $lStats = launchStatsProjectAll_genes();
		if (scalar(@$lStats) > 0) {
			print ".";
			$h->{by_genes} = $lStats;
			$ok = 1;
		}
	}
    if ($ok) {
    	my $tmp_dir = $project->getTmpDir();
	    my $session = new CGI::Session(undef, $cgi, {Directory=>$tmp_dir});
	    $session->param('hash_xls', compress(freeze $h));
	    print '@@@';
		print "\",\"session_id\":\"";
	    print $session->id;
	    print "\"}";
    }
    exit(0);
}



##### EXPORT VCF #####




sub export_vcf {
	my ($project, $pat_name, $hVarids) = @_;
	my ($hTmpFiles, $hToSort, $hHeaders, @lLines);
	my $dir_tmp = $project->getTmpDir();
	#my $dir_tmp = '/data-xfs/dev/mbras/tmp/';
	$dir_tmp .= '/'.$project->name().'_exportVCF/';
	if (not -d $dir_tmp) {
		mkdir($dir_tmp);
		`chmod 777 $dir_tmp`;
	}
	my @lLines;
	my $patient = $project->getPatient($pat_name);
	my @lVcfFiles;
	push(@lVcfFiles, @{$patient->getVariationsFiles()}); 
	push(@lVcfFiles, @{$patient->getIndelsFiles()}); 
	push(@lVcfFiles, @{$patient->getCnvsFiles()});
	my $nb = 0; 
	foreach my $file (@lVcfFiles) {
		next unless ($file =~ /vcf/);
		if ($file =~ /.gz/) { open(FILE, "zcat $file |"); }
		else { open(FILE, "$file"); }
		while(<FILE>) {
			chomp ($_);
			my $line = $_;
			if ($line =~ /#/) {
				push(@lLines, $line."\n");
			}
			else {
				my @lTmp = split(' ', $line);
				my $chr_id = $lTmp[0];
				$chr_id =~ s/chr//;  
				if (exists $hVarids->{$pat_name}->{$chr_id}->{$lTmp[1]}) {
					my @lTmp1 = split(',', $lTmp[3]);
					my $ref_vcf = join(',', @lTmp1);
					my @lTmp2 = split(',', $lTmp[4]);
					my $var_vcf = join(',', @lTmp2);
					my $vc = $hVarids->{$pat_name}->{$chr_id}->{$lTmp[1]}->{var_vcf};
					if ($ref_vcf eq $hVarids->{$pat_name}->{$chr_id}->{$lTmp[1]}->{ref_vcf} and $var_vcf =~/$vc/) {
						$lTmp[0] = 'chr'.$chr_id;
						join("\t", @lTmp);
						delete $hVarids->{$pat_name}->{$chr_id}->{$lTmp[1]};
						push(@lLines, $line."\n");
					}
				}
			}
		}
		close (FILE);
	}
	export_data::printVcfFile ('PolyQuery_'.$project->name(), $pat_name, \@lLines);
	exit(0);
}

sub export_vcf_multi_patients {
	my ($project, $patNames, $hVarids) = @_;
	my ($hTmpFiles, $hToSort, $hHeaders, @lLines);
	my $dir_tmp = $project->getTmpDir();
	#my $dir_tmp = '/data-xfs/dev/mbras/tmp/';
	$dir_tmp .= '/'.$project->name().'_exportVCF/';
	if (not -d $dir_tmp) {
		mkdir($dir_tmp);
		`chmod 777 $dir_tmp`;
	}
	my @lPatNames = split(' ', $patNames);
	foreach my $pat_name (@lPatNames) {
		my $patient = $project->getPatient($pat_name);
		my (@lVcfFiles, @lLinesHead);
		push(@lVcfFiles, @{$patient->getVariationsFiles()}); 
		push(@lVcfFiles, @{$patient->getIndelsFiles()}); 
		push(@lVcfFiles, @{$patient->getCnvsFiles()});
		my $nb = 0; 
		foreach my $file (@lVcfFiles) {
			next unless ($file =~ /vcf/);
			my $pass_header = 0;
			$hTmpFiles->{filtered}->{$pat_name}->{$nb} = File::Temp->new( TEMPLATE => $pat_name.'.'.$nb.'.XXXXX', DIR => $dir_tmp, SUFFIX => '.vcf' );
			open (OUT, '>'.$hTmpFiles->{filtered}->{$pat_name}->{$nb}->filename());
			if ($file =~ /.gz/) { open(FILE, "zcat $file |"); }
			else { open(FILE, "$file"); }
			while(<FILE>) {
				chomp ($_);
				my $line = $_;
				if ($line =~ /#/) {
					print OUT $line."\n";
					$pass_header = 1;
				}
				else {
					my @lTmp = split(' ', $line);
					my $chr_id = $lTmp[0];
					$chr_id =~ s/chr//;  
					if (exists $hVarids->{$pat_name}->{$chr_id}->{$lTmp[1]}) {
						my @lTmp1 = split(',', $lTmp[3]);
						my $ref_vcf = join(',', @lTmp1);
						my @lTmp2 = split(',', $lTmp[4]);
						my $var_vcf = join(',', @lTmp2);
						my $vc = $hVarids->{$pat_name}->{$chr_id}->{$lTmp[1]}->{var_vcf};
						if ($ref_vcf eq $hVarids->{$pat_name}->{$chr_id}->{$lTmp[1]}->{ref_vcf} and $var_vcf =~/$vc/) {
							$lTmp[0] = 'chr'.$chr_id;
							join("\t", @lTmp);
							delete $hVarids->{$pat_name}->{$chr_id}->{$lTmp[1]};
							print OUT $line."\n";
						}
					}
				}
				$hHeaders->{$pat_name} = \@lLinesHead;
			}
			close (FILE);
			close (OUT);
			$nb++ if ($pass_header);
		}
		$hTmpFiles->{unsort}->{$pat_name} = File::Temp->new( TEMPLATE => $pat_name.'.XXXXX', DIR => $dir_tmp, SUFFIX => '.vcf.gz' );
		my $cmd = 'cat ';
		foreach my $i (keys %{$hTmpFiles->{filtered}->{$pat_name}}) {
			$cmd .= ' '.$hTmpFiles->{filtered}->{$pat_name}->{$i}->filename();
		}
		$cmd .= ' | '.$project->buffer->software('bgzip').' -c >'.$hTmpFiles->{unsort}->{$pat_name}->filename();
		`$cmd`;
	}
	my @lPatVcf;
	foreach my $pat_name (@lPatNames) {
		$hTmpFiles->{sort}->{$pat_name} = File::Temp->new( TEMPLATE => $pat_name.'.XXXXX', DIR => $dir_tmp, SUFFIX => '.sort.vcf.gz' );
		my $cmd = $project->buffer->software('vcfsort');
		$cmd   .= ' '.$hTmpFiles->{unsort}->{$pat_name}->filename();
		$cmd   .= ' | '.$project->buffer->software('bgzip').' -c >'.$hTmpFiles->{sort}->{$pat_name}->filename();
		$cmd   .= ' ; '.$project->buffer->software('tabix').' -p vcf '.$hTmpFiles->{sort}->{$pat_name}->filename();
		system($cmd);
		push (@lPatVcf, $hTmpFiles->{sort}->{$pat_name}->filename());
	}
	$hTmpFiles->{unsort}->{global} = File::Temp->new( TEMPLATE => 'GLOBAL.XXXXX', DIR => $dir_tmp, SUFFIX => '.vcf' );
	if (scalar(@lPatVcf) == 1) {
		open (OUT, 'zcat '.$lPatVcf[0].' |');
		while(<OUT>) { push(@lLines, $_); warn $_; }
		close (OUT);
		die;
	}
	else {
		my $vcf_merge = VcfMerge->new( tabix=>$project->buffer->software('tabix'), list_VCF=>\@lPatVcf, VCF_out=>$hTmpFiles->{unsort}->{global}->filename(), silent=>1);
		$vcf_merge->merge_vcf_files();
		$hTmpFiles->{sort}->{global} = File::Temp->new( TEMPLATE => 'GLOBAL.XXXXX', DIR => $dir_tmp, SUFFIX => '.sort.vcf' );
		my $cmd1 = 'more '.$hTmpFiles->{unsort}->{global}->filename()." | grep '#' >" .$hTmpFiles->{sort}->{global}->filename();
		`$cmd1`;
		foreach my $chr_id (1..22, 'X', 'Y', 'MT') {
			my $cmd2 = 'sort -k2,2n '.$hTmpFiles->{unsort}->{global}->filename()." | grep '^$chr_id\t\\|^chr$chr_id\t' | grep -v '#' >>".$hTmpFiles->{sort}->{global}->filename();
			`$cmd2`;
		}
		open (OUT, $hTmpFiles->{sort}->{global}->filename());
		while(<OUT>) { push(@lLines, $_); }
		close (OUT);
	}
	export_data::printVcfFile ('PolyQuery_'.$project->name(), join('_', @lPatNames), \@lLines);
	exit(0);
}



##### XLS METHODS #####



# hash avec les infos necessaires sur ce patient pour un export en XLS
sub store_var_ids {
	my ($chr, $patient, $hash_global) = @_;
	if ($patient->project->get_xls() eq 'variants' or $patient->project->get_xls() eq 'genes') {
		my $hash = $patient->xls_var_id();
		my $hTrans;
		my $var_patient = $chr->saved_model->{'for_xls'}->{$patient->name()};
		my $var_patient_ho = $chr->saved_model->{'for_xls'}->{$patient->name().'_ho'};
		my @lVar = @{$chr->getListVarObjects($chr->getVariantsVector())};
		foreach my $var (@lVar) {
			$patient->project->print_dot(25);
			my $id = $var->vector_id();
			next unless ($var_patient->contains($id));
			my $var_id = $var->id();
			my $chr_h_id = $chr->id();
			$chr_h_id = '23' if ($chr->id() eq 'X');
			$chr_h_id = '24' if ($chr->id() eq 'Y');
			$chr_h_id = '25' if ($chr->id() eq 'MT');
			
			if (exists $hash_global->{$chr_h_id}->{$id}) {
				$hash->{$chr_h_id}->{$id} = $hash_global->{$chr_h_id}->{$id};
			}
			else {
				$hash->{$chr_h_id}->{$id}->{'var_id'} = $var_id;
				$hash->{$chr_h_id}->{$id}->{'var_id'} .= ' ('.$var->rs_name().')' if ($var->rs_name());
				$hash->{$chr_h_id}->{$id}->{'type'} = 'snp' if ($var->isVariation());
				$hash->{$chr_h_id}->{$id}->{'type'} = 'ins' if ($var->isInsertion());
				$hash->{$chr_h_id}->{$id}->{'type'} = 'del' if ($var->isDeletion());
				my $h_dejavu = $var->deja_vu();
				my $nb_project = 0;
				my $nb_patient = 0;
				my $he = 0;
				my $ho = 0;
				my $st_project;
				foreach my $projName (keys %$h_dejavu){
					$nb_project++;
					$st_project = $projName.":";
					$st_project .= $h_dejavu->{$projName}->{patients};
					$he += $h_dejavu->{$projName}->{he};
					$ho += $h_dejavu->{$projName}->{ho};
				}
				$hash->{$chr_h_id}->{$id}->{'dejavu'} = "";
				$nb_patient = $he + $ho;
				if ($nb_project > 0) {
					$hash->{$chr_h_id}->{$id}->{'dejavu'} = $nb_project."/".$nb_patient." ho:$ho,he:$he";
				}
				$hash->{$chr_h_id}->{$id}->{'dejavu'} .= ' (only HO)' if ($dejavu_ho);
				$hash->{$chr_h_id}->{$id}->{'chr'} = $chr->ucsc_name();
				$hash->{$chr_h_id}->{$id}->{'position'} = $var->start();
				$hash->{$chr_h_id}->{$id}->{'allele'} = $var->ref_allele().'/'.$var->var_allele();
				$hash->{$chr_h_id}->{$id}->{'hgmd_class'} = 'N.A.';
				if ($can_use_hgmd) {
					$hash->{$chr_h_id}->{$id}->{'hgmd_class'} = $var->hgmd_class();
				}
				$hash->{$chr_h_id}->{$id}->{'cadd_score'} = $var->cadd_score();
				if ($var->cosmic()) {
					my @lTmpCosmic = split(':', $var->cosmic());
					$hash->{$chr_h_id}->{$id}->{'cosmic'} = $lTmpCosmic[0];
				}
				$hash->{$chr_h_id}->{$id}->{'clinvar'} = $var->text_clinvar();
				$hash->{$chr_h_id}->{$id}->{'min_pop_freq'} = $var->min_pop_name().': '.$var->min_pop_freq();
				$hash->{$chr_h_id}->{$id}->{'max_pop_freq'} = $var->max_pop_name().': '.$var->max_pop_freq();
				
				my $seq = $chr->sequence(($var->start()-21), ($var->end()-1));
				$seq .= '['.$var->ref_allele().'/'.$var->var_allele().']';
				$seq .= $chr->sequence(($var->start()+1), ($var->end()+21));
				$hash->{$chr_h_id}->{$id}->{'sequence'} = $seq;
				foreach my $gene (@{$var->getGenes()}) {
					my $g_id = $gene->id();
					$gene->chromosome($chr);
					next unless ($gene->getVariantsVector->contains($id));
					$var->annotation();
					$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'external_name'} = $gene->external_name();
					$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'description'} = $gene->description();
					foreach my $t (@{$gene->getTranscripts()}) {
						next unless (exists $var->annotation->{$t->id()});
						my @ok;
						my $annot_trans;
						eval { $annot_trans = $var->variationType($t) };
						if ($@) { $annot_trans = 'N.A'; }
						foreach my $cons (split(',', $annot_trans)) {
							my @lTmp = split(';', $patient->project->buffer->config->{ensembl_annotations}->{$cons});
							push(@ok, $lTmp[1]);
						}
						next if (scalar(@ok) == 0);
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'external_name'} = $t->external_name();
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'consequence'} = join(',', @ok);
						foreach my $cons (@ok) {
							$hash->{$chr_h_id}->{$id}->{'consequences'}->{$cons} = undef;
						}
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'cdna_position'} = $t->translate_position($var->start());
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'exon'} = $t->findExonNumber($var->start());
						if ($hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'exon'} == -1) {
							$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'exon'} = $t->findNearestExon($var->start(), $var->end());
						}
						if ($var->isCoding($t)) {
							my $prot = $t->getProtein();
							if ($prot) {
								$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'protein'} = $prot->id();
								$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'protein_xref'} = $t->{'external_protein_name'};
								$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'nomenclature'} = $var->getNomenclature($t);
								$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'cds_position'} = $var->getOrfPosition($prot);
								$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'protein_position'} = $var->getProteinPosition($prot);
								$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'polyphen_status'} = $var->polyphenStatus($prot);
								$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'sift_status'} = $var->siftStatus($prot);
								$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'polyphen_score'} = $var->polyphenScore($prot);
								$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'sift_score'} = $var->siftScore($prot);
								my $protAA = $var->getProteinAA($prot);
								my $chanAA = $var->changeAA($prot);
								if ($protAA and $chanAA) {
									$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'aa'} = $protAA.'/'.$chanAA;
								}
							}
							else {
								$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'protein'} = '';
								$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'protein_xref'} = '';
								$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'nomenclature'} = '';
								$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'cds_position'} = '';
								$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'protein_position'} = '';
								$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'aa'} = '';
							}
						}
						else {
							$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'protein'} = '';
							$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'protein_xref'} = '';
							$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'nomenclature'} = '';
							$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'cds_position'} = '';
							$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'protein_position'} = '';
							$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'aa'} = '';
						}
					}
					$hash_global->{$chr_h_id}->{$id} = $hash->{$chr_h_id}->{$id};
				}
			}
			$hash->{$chr_h_id}->{$id}->{'nb_all_ref'}->{$patient->name()} = $var->{annex}->{$patient->id()}->{nb_all_ref};
			$hash->{$chr_h_id}->{$id}->{'nb_all_mut'}->{$patient->name()} = $var->{annex}->{$patient->id()}->{nb_all_mut};
			if ($var_patient_ho->contains($id)) { $hash->{$chr_h_id}->{$id}->{'he_ho'}->{$patient->name()} = 'ho'; }
			else { $hash->{$chr_h_id}->{$id}->{'he_ho'}->{$patient->name()} = 'he' }
		}
		$hTrans = undef;
		$patient->xls_var_id($hash);
	}
	return $hash_global;
}

sub getXls_byGenes {
	my ($project, $hashRes, $hResumeFilters, $xls_outfile) = @_;
	my $projectName = $project->name();
	my $workbook;
	if ($xls_outfile) {
		$workbook = Spreadsheet::WriteExcel->new( $xls_outfile );
	}
	else {
		print "Content-type: application/msexcel\n";
		print "Content-Disposition: attachment;filename=$projectName\_byGenes.xls\n\n";
		$workbook = Spreadsheet::WriteExcel->new( \*STDOUT );
	}
	my $hFormat = getHashFormat($workbook);
	my $xls_gene = $workbook->add_worksheet($projectName);
	my $listColNames = writeHeader_byGenes($workbook, $xls_gene);
	my @lGenes = @{$hashRes->{by_genes}};
	my $hLen = initHashSizeColumn_byVar($listColNames);
	my @lLines;
	foreach my $h (@lGenes) {
		if ($project->only_genes()) {
			next unless (exists $project->{only_genes}->{$h->{xref}});
		}
		my @lCol;
		push(@lCol, $h->{name});
		push(@lCol, $h->{xref});
		push(@lCol, $h->{chromosome});
		push(@lCol, $h->{start});
		push(@lCol, $h->{end});
		push(@lCol, $h->{description});
		push(@lCol, $h->{v_all});
		push(@lCol, $h->{p_all});
		push(@lCol, $h->{v_substitution});
		push(@lCol, $h->{p_substitution});
		push(@lCol, $h->{v_insertion});
		push(@lCol, $h->{p_insertion});
		push(@lCol, $h->{v_deletion});
		push(@lCol, $h->{p_deletion});
		push(@lCol, $h->{v_coding});
		push(@lCol, $h->{p_coding});
		push(@lCol, $h->{v_silent});
		push(@lCol, $h->{p_silent});
		push(@lCol, $h->{v_utr});
		push(@lCol, $h->{p_utr});
		push(@lCol, $h->{v_splicing});
		push(@lCol, $h->{p_splicing});
		push(@lCol, $h->{v_stop});
		push(@lCol, $h->{p_stop});
		push(@lCol, $h->{v_phase});
		push(@lCol, $h->{p_phase});
		push(@lCol, $h->{v_silent});
		push(@lCol, $h->{p_silent});
		push(@lCol, $h->{v_frameshift});
		push(@lCol, $h->{p_frameshift});
		push(@lLines, \@lCol);
	}
	my $i = 1;
	foreach my $listCol (@lLines) {
		my $j = 0;
		foreach my $value (@$listCol) {
			$xls_gene->write($i, $j, $value, $hLen->{$j});
			$j++;
		}
		$i++;
	}
	my $xls_resume = $workbook->add_worksheet('FILTERS USED');
	writeResumeFiltersXls($xls_resume, $hFormat, $hResumeFilters);
	exit(0);
}

sub getXls_byVar {
	my ($project, $h, $hResumeFilters, $xls_outfile) = @_;
	my $projectName = $project->name();
	my $workbook;
	if ($xls_outfile) {
		$workbook = Spreadsheet::WriteExcel->new( $xls_outfile );
	}
	else {
		print "Content-type: application/msexcel\n";
		print "Content-Disposition: attachment;filename=$projectName\_byVar.xls\n\n";
		$workbook = Spreadsheet::WriteExcel->new( \*STDOUT );
	}
	unless ($h) {
		foreach my $patient (@{$project->getPatients()}) {
			$h->{by_pat}->{$patient->name()} = $patient->xls_var_id();
			foreach my $chr_id (keys %{$h->{by_pat}->{$patient->name()}}) {
				foreach my $id (keys %{$h->{by_pat}->{$patient->name()}->{$chr_id}}) {
					$h->{by_var}->{$chr_id}->{$id}->{$patient->name()} = undef; 
				}
			}
		}
	}
	my $hFormat = getHashFormat($workbook);
	my $nb_page = 0;
	my @lAllPat = @{$project->getPatients()};
	my $nb_all_pat_project = scalar(@lAllPat);
	while (scalar(@lAllPat) > 0) {
		my (@lThisPatients, $hPatPass);
		my $nb_pat_pass = 0;
		while ($nb_pat_pass < 200) {
			my $pat = shift(@lAllPat);
			push(@lThisPatients, $pat);
			$hPatPass->{$pat->name()} = undef;
			$nb_pat_pass++;
			unless (@lAllPat) { $nb_pat_pass = 100000; }
		}
		my (@lPatNamesHeader, @lPatNames, $hPolyphenSift);
		foreach my $pat (@lThisPatients) {
			my $aff_or_not = 'healthy';
			$aff_or_not = 'affected' if ($pat->status() eq '2');
			push(@lPatNamesHeader, $pat->name().':'.$aff_or_not);
			push(@lPatNames, $pat->name());
		}
		$nb_page++;
		my $xls_var_name = 'RESULTS part '.$nb_page;
		my $xls_var_not_merge_name = 'RESULTS NOT MERGED part '.$nb_page;
		
		my $xls_var = $workbook->add_worksheet($xls_var_name);
		my $xls_var_not_merge = $workbook->add_worksheet($xls_var_not_merge_name);
		my $listColNames = writeHeader_byVar($workbook, $xls_var, \@lPatNamesHeader);
		my $listColNames_not_merge = writeHeader_byVar($workbook, $xls_var_not_merge, \@lPatNamesHeader);
		my $i = 1;
		
		my $col_cons = scalar(@lPatNames) + 17;
		my $col_polyphen = $col_cons + 12;
		my $col_sift = $col_cons + 13;
		
		my @lLines;
		foreach my $chr_id (sort {$a <=> $b} keys %{$h->{by_var}}) {
			foreach my $id (sort {$a <=> $b} keys %{$h->{by_var}->{$chr_id}}) {
				my @lThisPatNames = keys %{$h->{by_var}->{$chr_id}->{$id}};
				my $firstPat = $lThisPatNames[0];
				my @lCol;
				push(@lCol, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{var_id});
				push(@lCol, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{type});
				push(@lCol, join(', ', sort keys %{$h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{consequences}}));
				push(@lCol, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{dejavu});
				push(@lCol, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{chr});
				push(@lCol, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{position});
				push(@lCol, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{allele});
				push(@lCol, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{sequence});
				my $nb_he = 0;
				my $nb_ho = 0;
				
				foreach my $pat_name (sort @lPatNames) {
					if (exists $h->{by_var}->{$chr_id}->{$id}->{$pat_name}) {
						my $he_ho = $h->{by_pat}->{$pat_name}->{$chr_id}->{$id}->{he_ho}->{$pat_name};
						my $nb_ref = $h->{by_pat}->{$pat_name}->{$chr_id}->{$id}->{nb_all_ref}->{$pat_name};
						my $nb_mut = $h->{by_pat}->{$pat_name}->{$chr_id}->{$id}->{nb_all_mut}->{$pat_name};
						my $ratio;
						unless (int($nb_ref) == 0 and int($nb_mut) == 0) {
							$ratio = int((int($nb_mut) / (int($nb_ref) + int($nb_mut))) * 100);
						}
						my $text = uc($he_ho).' ('.$nb_ref.'/'.$nb_mut.')';
						$text .= ' ('.$ratio.'%)' if ($ratio);
						push(@lCol, $text);
						if ($he_ho eq 'he')    { $nb_he++; }
						elsif ($he_ho eq 'ho') { $nb_ho++; }
					}
					else { push(@lCol, ''); }
				}
				push(@lCol, $nb_he);
				push(@lCol, $nb_ho);
				if (exists $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{cadd_score}) {
					if ($h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{cadd_score} eq '-1') { push(@lCol, '.'); }
					else { push(@lCol, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{cadd_score}); }
				}
				else { push(@lCol, '.'); }
				if (exists $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{cosmic}) {
					if ($h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{cosmic} eq '-1') { push(@lCol, '.'); }
					else { push(@lCol, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{cosmic}); }
				}
				else { push(@lCol, '.'); }
				my $clinvar = $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{clinvar};
				if ($clinvar eq '-5') { $clinvar = '.'; }
				push(@lCol, $clinvar);
				push(@lCol, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{min_pop_freq});
				push(@lCol, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{max_pop_freq});
				push(@lCol, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{hgmd_class});
				
				my @lGenes;
				unless (exists $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}) {
					$lCol[2] = 'Intergenic';
					my @lThisGene;
					push(@lThisGene, 'Intergenic');
					push(@lThisGene, 'Intergenic');
					push(@lThisGene, '');
					push(@lThisGene, '');
					push(@lThisGene, 'Intergenic');
					push(@lThisGene, '');
					push(@lThisGene, '-1');
					push(@lGenes, join("\t", @lThisGene));
				}
				foreach my $gene_name (keys %{$h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}}) {
					unless (exists $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{transcripts}) {
						#$lCol[2] = 'Intergenic';
						my @lThisGene;
						push(@lThisGene, $gene_name);
						push(@lThisGene, 'Intergenic');
						push(@lThisGene, '');
						push(@lThisGene, '');
						push(@lThisGene, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{description});
						push(@lThisGene, '');
						push(@lThisGene, '-1');
						push(@lGenes, join("\t", @lThisGene));
						next;
					}
					my @lTransName = keys %{$h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{transcripts}};
					foreach my $trans_name (@lTransName) {
						my @lThisGene;
						my $thisGeneName = $gene_name;
						if (exists $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{external_name}) {
							$thisGeneName .= ' ['.$h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{external_name}.']';
						}
						push(@lThisGene, $thisGeneName);
						push(@lThisGene, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{transcripts}->{$trans_name}->{consequence});
						push(@lThisGene, $trans_name);
						push(@lThisGene, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{transcripts}->{$trans_name}->{external_name});
						push(@lThisGene, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{description});
						push(@lThisGene, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{transcripts}->{$trans_name}->{exon});
						push(@lThisGene, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{transcripts}->{$trans_name}->{cdna_position});
						push(@lThisGene, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{transcripts}->{$trans_name}->{cds_position});
						push(@lThisGene, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{transcripts}->{$trans_name}->{protein});
						push(@lThisGene, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{transcripts}->{$trans_name}->{protein_xref});
						push(@lThisGene, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{transcripts}->{$trans_name}->{aa});
						push(@lThisGene, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{transcripts}->{$trans_name}->{protein_position});
						push(@lThisGene, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{transcripts}->{$trans_name}->{nomenclature});
						
						if (exists $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{transcripts}->{$trans_name}->{polyphen_score}) {
							my $status = $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{transcripts}->{$trans_name}->{polyphen_status};
							push(@lThisGene, $status.';'.$h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{transcripts}->{$trans_name}->{polyphen_score});
						}
						else { push(@lThisGene, ''); }
						if (exists $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{transcripts}->{$trans_name}->{sift_score}) {
							my $status = $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{transcripts}->{$trans_name}->{sift_status};
							push(@lThisGene, $status.';'.$h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{transcripts}->{$trans_name}->{sift_score});
						}
						else { push(@lThisGene, ''); }
						push(@lGenes, join("\t", @lThisGene));
					}
				}
				#foreach my $pat_name (@lPatNames) {
				#	delete $h->{by_pat}->{$pat_name}->{$chr_id}->{$id} if (exists $h->{by_pat}->{$pat_name}->{$chr_id}->{$id});
				#}
				#delete $h->{by_var}->{$chr_id}->{$id};
				
				my $line = join("\t", @lCol);
				my $hToMerge;
				my $done = 0;
				my ($start, $end);
				if (scalar(@lGenes) > 0) {
					foreach my $end_line (@lGenes) {
						my $this_line;
						if($done == 0) {
							$this_line = $line;
							$start = $i;
							$end = $i;
							$done++;
						}
						else {
							my $z = 0;
							while ($z < (15 + scalar(@lPatNames))) {
								$this_line .= "\t";
								$z++;
							}
							$end = $i;
						}
						$this_line .= "\t".$end_line;
						my @lCol = split("\t", $this_line);
						my $j = 0;
						
						foreach my $value (@lCol) {
							my $format = $hFormat->{line};
							if ($j == $col_polyphen) {
								$format = getHashFormatForPolyphenSift($workbook, $value, 'polyphen');
								if ($value) {
									my @lTmp = split(';', $value);
									$value = '';
									if ($lTmp[0] == 1) { $value = 'Benign (score: '.$lTmp[-1].')'; }
									elsif ($lTmp[0] == 2) { $value = 'Possibly Damaging (score: '.$lTmp[-1].')'; }
									elsif ($lTmp[0] == 3) { $value = 'Probably Damaging (score: '.$lTmp[-1].')'; }
								}
							}
							elsif ($j == $col_sift) {
								$format = getHashFormatForPolyphenSift($workbook, $value, 'sift');
								if ($value) {
									my @lTmp = split(';', $value);
									$value = '';
									if ($lTmp[0] == 1) { $value = 'Benign (score: '.$lTmp[-1].')'; }
									elsif ($lTmp[0] == 2) { $value = 'Deleterious (score: '.$lTmp[-1].')'; }
								}
							}
							elsif ($j == $col_cons) {
								$format = getHashFormatForConsequence($workbook, $value);
							}
							if ($j < (16 + scalar(@lPatNames))) {
								unless (exists $hToMerge->{$j}) {
									$hToMerge->{$j}->{value}  = $value;
									$hToMerge->{$j}->{format} = $format;
								}
							}
							else {
								$xls_var->write($i, $j, $value, $format);
								$xls_var_not_merge->write($i, $j, $value, $format);
							}
							$j++;
						}
						$i++;
						$done++;
					}
				}
				unless ($start == $end) {
					foreach my $j (sort {$a <=> $b} keys %$hToMerge) {
						$xls_var->merge_range($start, $j, $end, $j, $hToMerge->{$j}->{value}, $hFormat->{merge});
					}
					for ( my $i = $start; $i<= $end; $i++ ) {
						foreach my $j (sort {$a <=> $b} keys %$hToMerge) {
							$xls_var_not_merge->write($i, $j, $hToMerge->{$j}->{value}, $hToMerge->{$j}->{format});
						}
					}
				}
				else {
					foreach my $j (sort {$a <=> $b} keys %$hToMerge) {
						$xls_var->write($start, $j, $hToMerge->{$j}->{value}, $hToMerge->{$j}->{format});
						$xls_var_not_merge->write($start, $j, $hToMerge->{$j}->{value}, $hToMerge->{$j}->{format});
					}
				}
			}
		}
	}
		
	my $xls_resume = $workbook->add_worksheet('FILTERS USED');
	writeResumeFiltersXls($xls_resume, $hFormat, $hResumeFilters);
	exit(0);
}

sub writeResumeFiltersXls {
	my ($xls_resume, $hFormat, $hResumeFilters) = @_;
	my $x = 0;
	$xls_resume->write($x, 0, 'Ind / Fam', $hFormat->{header});
	$xls_resume->write($x, 1, $hResumeFilters->{typeFilters});
	$x++;
	$xls_resume->write($x, 0, 'Model', $hFormat->{header});
	if ($hResumeFilters->{model}) { $xls_resume->write($x, 1, $hResumeFilters->{model}); }
	else { $xls_resume->write($x, 1, 'None'); }
	$x++;
	$xls_resume->write($x, 0, 'by Var / Genes', $hFormat->{header});
	$xls_resume->write($x, 1, $hResumeFilters->{typeVarOrGenes});
	$x++;
	$xls_resume->write($x, 0, 'DejaVu', $hFormat->{header});
	if ($hResumeFilters->{dejavu}) { $xls_resume->write($x, 1, $hResumeFilters->{dejavu}); }
	else { $xls_resume->write($x, 1, 'None'); }
	$x++;
	$xls_resume->write($x, 0, 'At Least', $hFormat->{header});
	if ($hResumeFilters->{atLeast}) { $xls_resume->write($x, 1, $hResumeFilters->{atLeast}); }
	else { $xls_resume->write($x, 1, 'None'); }
	
	my @lNames;
	if ($hResumeFilters->{typeFilters} eq 'Individual') { @lNames = ('filter_region', 'filter_patient', 'filter_attic', 'filter_not_patient', 'filter_ho', 'filter_he', 'filter_type_variation'); }
	else { @lNames = ('filter_region', 'fam_not', 'fam_and', 'filter_patient', 'filter_attic', 'filter_not_patient', 'filter_ho', 'filter_he', 'filter_type_variation'); }
	foreach my $name (@lNames) {
		if (scalar @{$hResumeFilters->{$name}} > 0) {
			$x += 2;
			my @lTmp = split('_', $name);
			shift(@lTmp);
			my $this_name = join(' ', @lTmp);
			$xls_resume->write($x, 0, uc('Filter(s) '.$this_name), $hFormat->{header});
			foreach my $filter (@{$hResumeFilters->{$name}}) {
				if ($name eq 'filter_type_variation') {
					my $var_cons_new = $project->buffer->config->{ensembl_annotations}->{$filter};
					if ($var_cons_new) {
						my ($name, $name2) = split(';', $var_cons_new);
						$filter = $name;
					}
				}
				$xls_resume->write($x, 1, $filter);
				$x++;
			}
		}
	}
}

sub getHashFormat {
	my $xls = shift;
	my $hash;
	my $format_header = $xls->add_format(border => 1, underline => 1);
	$format_header->set_bold();
	$format_header->set_align('center');
	$format_header->set_fg_color('silver');
	my $format_line = $xls->add_format();
	$format_line->set_align('center');
	$format_line->set_color('black');
	my $format_merge = $xls->add_format(valign => 'vcenter');
	$format_merge->set_align('center');
	$format_merge->set_color('black');
	$hash->{header} = $format_header;
	$hash->{line} = $format_line;
	$hash->{merge}= $format_merge;
	return $hash;
}

sub getHashFormatForConsequence {
	my ($xls, $cons) = @_;
	my $format = $xls->add_format();
	$format->set_align('center');
	$format->set_color('black');
	my $h;
	foreach my $this_cons (split(',', $cons)) {
		if ($this_cons eq 'Intergenic' or $this_cons eq 'Intronic') {
	    	$h->{green} = 1;
		}
		elsif ($this_cons eq 'Synonymous' or $this_cons eq 'Utr' or $this_cons eq 'Pseudogene') {
	    	$h->{yellow} = 1;
		}
		elsif ($this_cons eq 'No-frameshift' or $this_cons eq 'Missense' or $this_cons eq 'Splice Region' or $this_cons eq 'ncRNA') {
	    	$h->{orange} = 1;
		}
		elsif ($this_cons eq '(Start/Stop)-lost' or $this_cons eq 'Stop-gained' or $this_cons eq 'Frameshift' or $this_cons eq 'Splice Acc/Don' or $this_cons eq 'Mature miRNA') {
	    	$h->{red} = 1;
		}
	}
	if    (exists $h->{red})    { $format->set_fg_color('red'); }
	elsif (exists $h->{orange}) { $format->set_fg_color('orange'); }
	elsif (exists $h->{yellow}) { $format->set_fg_color('yellow'); }
	elsif (exists $h->{green})  { $format->set_fg_color('green'); }
	return $format;
}

sub getHashFormatForPolyphenSift {
	my ($xls, $value, $type) = @_;
	my $format = $xls->add_format();
	$format->set_align('center');
	$format->set_color('black');
	return unless ($value);
	my ($status, $freq) = split(';', $value);
	if    (int($status) == 3) { $format->set_fg_color('red'); }
	elsif (int($status) == 2) {
		$format->set_fg_color('red') if ($type eq 'sift');
		$format->set_fg_color('orange') if ($type eq 'polyphen');
	}
	elsif (int($status) == 1) { $format->set_fg_color('green'); }
	return $format;
}

sub writeHeader_byVar {
	my ($xls, $xls_pat, $listPatNamesStatus) = @_;
	my $hStatus;
	foreach my $field (sort @$listPatNamesStatus) {
		my ($name, $status) = (split(':', $field));
		$hStatus->{$name} = $status;
	}
	my @listPatNames = sort keys %$hStatus;
	my $header = "Variation	Type	Consequence	Dejavu	Chr	Position	Allele	Sequence	".join("\t", @listPatNames)."	He	Ho	Cadd	Cosmic	ClinVar	Min_Pop_Freq	Max_Pop_Freq	HGMD_Class	Gene	Consequence	Transcript	Transcript_Xref	Description	Exon	Cdna_Pos	Cds_Pos	Protein	Protein_xref	AA	Protein_Pos	Nomenclature	Polyphen	Sift";
	my @lNames = split(' ', $header);
	my $hFormat = getHashFormat($xls);
	my $j = 8;
	foreach my $name (@listPatNames) {
		my $status = $hStatus->{$name};
		my $path_polyweb;
		if (-d $Bin.'/../../polyweb/') { $path_polyweb = $Bin.'/../../polyweb/'; }
		elsif (-d $Bin.'/../../PolyWeb/') { $path_polyweb = $Bin.'/../../PolyWeb/'; }
		else { warn $Bin; die; }
		my $icon = "$path_polyweb/images/polyicons/bullet_green.png";
		$icon = "$path_polyweb/images/polyicons/pill2.png" if ($status eq 'affected');
		$xls_pat->insert_image(0, $j, $icon);
		$j++;
	}
	$j = 0;
	foreach my $col_name (@lNames) {
		$xls_pat->write(0, $j, $col_name, $hFormat->{header});
		$j++;
	}
	return \@lNames;
}

sub writeHeader_byRegionsHo {
	my ($xls, $xls_pat) = @_;
	my $header = "Chr	Start	End	Length	Fam	Patient		Nb_var	Nb_genes	Genes";
	my @lNames = split(' ', $header);
	my $hFormat = getHashFormat($xls);
	my $j = 0;
	foreach my $col_name (@lNames) {
		$xls_pat->write(0, $j, $col_name, $hFormat->{header});
		$j++;
	}
	return \@lNames;
	
}

sub writeHeader_byGenes {
	my ($xls, $xls_pat) = @_;
	my $header = "Name	xref	chr	Start	End	Description	All	Pat.	sub	Pat.	ins	Pat.	del	Pat.	Coding	Pat.	Silent	Pat.	UTR	Pat.	Splicing	Pat.	Stop	Pat.	start_stop_lost	Pat.	Synonymous	Pat.	Frameshift	Pat.";
	my @lNames = split(' ', $header);
	my $hFormat = getHashFormat($xls);
	my $j = 0;
	foreach my $col_name (@lNames) {
		$xls_pat->write(0, $j, $col_name, $hFormat->{header});
		$j++;
	}
	return \@lNames;
}

sub initHashSizeColumn_byVar {
	my $listNames = shift;
	my $hash;
	my $i = 0;
	foreach my $name (@$listNames) {
		$hash->{$i} = length($name);
		$i++;
	}
	return $hash;
}

sub resizeXlsColumn {
	my ($worksheet, $hLen) = @_;
	foreach my $col (sort {$a <=> $b} keys %$hLen) {
		my $len = $hLen->{$col} + 4;
		$worksheet->set_column($col, $col,  $len);
	}
}

sub checkEssentialsCategoriesInChr {
	my $chr = shift;
	return 1 if ($chr->id() eq 'MT');
	my $hCommons = Cache_Commons::categories->{global};
	foreach my $type (keys %{$hCommons}){
		foreach my $cat (keys %{$hCommons->{$type}}) {
			next if ($cat eq 'large_deletion');
			next if ($cat eq 'large_duplication');
			next if ($cat eq 'spliceAI_high');
			next if ($cat eq 'spliceAI_medium');
			next if ($cat eq 'predicted_splice_site');
			unless (exists $chr->global_categories->{$cat}) {
				warn "\n\n";
				warn Dumper $chr->global_categories();
				confess("\n\nERROR: [CHR".$chr->id()."] global_categories->{$cat} not found. Pb cache ?? Die.\n\n");
			}
		}
	}
	return 1;
}
