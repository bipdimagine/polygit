#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;

use strict;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";

use lib "$Bin/../GenBo/lib/obj-nodb/packages/";

use lib "$Bin/../packages/export/";
use lib "$Bin/../packages/layout/";
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
use QueryVectorFilter;
use Storable qw(store retrieve freeze dclone thaw);
use POSIX qw(strftime);
use CGI::Session;
use html; 
use Carp;
use Cache_Commons;
use List::MoreUtils qw{ natatime };
use xls_export;
use session_export;
use threads;

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
my $filter_bed				= $cgi->param('filter_bed');
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
my $filter_ratio_min		= $cgi->param('min_ratio');
my $filter_ratio_max		= $cgi->param('max_ratio');
my $filter_ncboost			= $cgi->param('ncboost');
my $panel_name				= $cgi->param('panel');
my $annot_version			= $cgi->param('annot_version');
$annot_version = undef;
my $keep_indels_cadd		= $cgi->param('keep_indels_cadd');
#my $filter_gnomad_test		= $cgi->param('gnomad_test');


my $queryFilter =  new QueryVectorFilter;
$queryFilter->verbose_debug(1) if $debug;

my $hDeleteModels;
if ($delete_models) {
	foreach my $model (split(' ', $delete_models)) {
		$hDeleteModels->{$model} = undef;
	}
}

if ($xls_outfile eq 'none') {
	#modif temporaire pour export xls
		$xls_outfile = $projectName.".xls";
		
		#warn "\n\nERROR: please replace value 'none' for last -xls_outfile option.\nExample: -xls_outfile=export.xls\n\n";
	#die;
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
else {
	$buffer = new GBuffer; 
}

my $project = $buffer->newProjectCache( -name 			=> $projectName,
									    -cache 		=> '1',
									    -typeFilters 	=> $typeFilters, );
									    


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

my $h_chr_from_only_genes;
if ($only_genes) {
	foreach my $name (split(',', $only_genes)) {
#		eval {
#			my $buffer_tmp = new GBuffer;
#			my $project_tmp = $buffer_tmp->newProject( -name => $projectName );
#			my $this_gene = $project_tmp->newGene($name);
#			$h_chr_from_only_genes->{$this_gene->getChromosome->id()}->{$this_gene->id()} = undef;
#		};
#		if ($@) {}
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
my ($hFiltersChr, $hFiltersChr_var2) = get_hashes_filters($filter_type_variation, $filter_type_variation_2);

if ($model eq 'recessif' or $model eq 'compound' or $model eq 'recessif_compound') {
	$hFiltersChr->{intergenic} = undef;
	if ($filter_type_variation_2) {
		$hFiltersChr_var2->{intergenic} = undef;
	}
}

if ($infos) {
	get_infos();
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
	@{$hResumeFilters->{filter_region}} = split(' ', $filter_region);
	@{$hResumeFilters->{filter_attic}} = split(' ', $filter_attic);
	@{$hResumeFilters->{filter_patient}} = split(' ', $filter_patient);
	@{$hResumeFilters->{filter_not_patient}} = split(' ', $filter_not_patient);
	@{$hResumeFilters->{filter_he}} = split(' ', $filter_he);
	@{$hResumeFilters->{filter_ho}} = split(' ', $filter_ho);
	@{$hResumeFilters->{fam_not}} = split(' ', $fam_not);
	@{$hResumeFilters->{fam_and}} = split(' ', $fam_and);
}

my $h_filters_intersect_exclude;
$h_filters_intersect_exclude->{filter_he} = $filter_he;
$h_filters_intersect_exclude->{filter_ho} = $filter_ho;
$h_filters_intersect_exclude->{filter_not_patient} = $filter_not_patient;
$h_filters_intersect_exclude->{fam_not} = $fam_not;
$h_filters_intersect_exclude->{filter_nbvar_regionho} = $filter_nbvar_regionho;
$h_filters_intersect_exclude->{filter_patient} = $filter_patient;
$h_filters_intersect_exclude->{fam_and} = $fam_and;
$h_filters_intersect_exclude->{filter_patient} = $filter_patient;

if ($xls_load_session) {
	loadSessionsXLS($project, $xls_load_session, $hResumeFilters);
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

my @l_genome_fai = @{$project->getGenomeFai};

foreach my $chr (@{$project->getChromosomes()}) {
	$chr->not_used(1);
}

	$queryFilter->setPatients($project);
	$queryFilter->setIntheAttic($project->getPatientsFromListNames([split(' ', $filter_attic)]));
	
	$queryFilter->setIntersectedPatients( $project->getPatientsFromListNames([split(' ', $h_filters_intersect_exclude->{filter_patient})]));
	$queryFilter->setIntersectedFamilies($project->getFamiliesFromListNames([split(' ', $h_filters_intersect_exclude->{fam_and})]));
	
	$queryFilter->setFilteredHoPatients($project->getFamiliesFromListNames([split(' ', $h_filters_intersect_exclude->{filter_ho})]));
	$queryFilter->setFilteredHePatients($project->getFamiliesFromListNames([split(' ', $h_filters_intersect_exclude->{filter_he})]));
	
	$queryFilter->setExcludedPatients( $project->getPatientsFromListNames([split(' ', $h_filters_intersect_exclude->{filter_not_patient})]));
	$queryFilter->setExcludedFamilies($project->getFamiliesFromListNames([split(' ', $h_filters_intersect_exclude->{fam_not})]));

foreach my $chr_id (sort split(',', $filter_chromosome)) {
	# Cas genome MT uniquement par exemple
	
	my $chr_found;
	foreach my $h_fai (@l_genome_fai) {
		$chr_found = 1 if (exists $h_fai->{chromosomes_object}->{$chr_id});
		last if $chr_found;
	}
	if (not $chr_found)  {
		next if ($export_vcf_for or $detail_project or $xls_by_regions_ho);
		my $hashRes;
		$hashRes->{'label'} = 'name';
		$hashRes->{'items'} = launchStatsProjectAll_chromosomes_null();
		printJson($hashRes);
		exit(0);
	}
	my $chr = $project->getChromosome($chr_id);
	$chr->not_used(0);
	$nb_chr++;
	
	#delete variant with in the attic
	$queryFilter->FilterInTheAtticPatients($chr);
	
	if ($debug) { warn "\n\nCHR ".$chr->id()." -> INIT - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	
	if (not $export_vcf_for and not $detail_project and not $xls_by_regions_ho) { $project->cgi_object(1); }
	
	unless (check_region_filter($chr, $filter_region)) {
		print "@" unless ($export_vcf_for or $detail_project or $xls_by_regions_ho);
		next;
	}
	if ($xls_save_session) {
		print "@" unless ($export_vcf_for or $detail_project or $xls_by_regions_ho);
	}
	
	# Pour XLS, besoin de savoir si un patient possede ou non le variant, meme s il ne passe pas un filtre ou un modele pour rester coherent a la 2e page par gene
	if ($xls_by_variants) { $chr->save_model_variants_all_patients('for_xls'); }
	
	foreach my $p (@{$project->getPatients()}) { $p->setOrigin($chr); }
	
#	$queryFilter->setInTheAtticPatients($chr, $project->getPatientsFromListNames([split(' ', $filter_attic)]));
	print "@" unless ($export_vcf_for or $detail_project or $xls_by_regions_ho);
	
	#FILTERING ON VARIANTS 
	
	
#	launch_intersect_excude($chr, $h_filters_intersect_exclude) if not $model =~ /compound/ and not $model_2 =~ /compound/;
	
	$queryFilter->filter_vector_ratio($chr, $filter_ratio_min, 'min');
	$queryFilter->filter_vector_ratio($chr, $filter_ratio_max, 'max');
	$queryFilter->filter_vector_ncboost($chr, $filter_ncboost);
#	$queryFilter->filter_vector_global_gnomad_freq($chr, $filter_gnomad_test);
#	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER getVectorGnomadCategory - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
		
	my $h_args;	
	if ($hFiltersChr and $hFiltersChr_var2) { $chr->save_model_variants_all_patients('init'); }
	$queryFilter->filter_vector_region_ho($chr, $filter_nbvar_regionho, $filter_regionho_sub_only, $project->typeFilters());
	$queryFilter->filter_vector_type_variants($chr, $hFiltersChr);
	$queryFilter->filter_vector_cadd_variants($chr, $hFiltersChr, $keep_indels_cadd);
	if ($filter_gnomad) { $queryFilter->filter_vector_gnomad_ac($chr, $filter_gnomad) }
	else { $queryFilter->filter_vector_freq_variants($chr, $hFiltersChr); }
	$queryFilter->filter_vector_gnomad_ho_ac_variants($chr, $hFiltersChr);
	$queryFilter->filter_vector_confidence_variants($chr, $hFiltersChr);
	$queryFilter->filter_vector_dejavu($chr, $dejavu, $dejavu_ho, $test) if ($dejavu);

	# Recessif compound multi annot
	my $vector_filtered = $chr->getVariantsVector->Clone();
	my $vector_filtered_2;
	my $is_diff_hash_filters = is_differents_hash_filters($hFiltersChr, $hFiltersChr_var2, $dejavu, $dejavu_2);
	if ($is_diff_hash_filters) {
		die();
		$chr->load_init_variants_all_patients('init');
		$queryFilter->filter_vector_type_variants($chr, $hFiltersChr_var2);
		$queryFilter->filter_vector_cadd_variants($chr, $hFiltersChr_var2, $keep_indels_cadd);
		if ($filter_gnomad) { $queryFilter->filter_vector_gnomad_ac($chr, $hFiltersChr_var2) }
		else { $queryFilter->filter_vector_freq_variants($chr, $hFiltersChr_var2); }
		$queryFilter->filter_vector_gnomad_ho_ac_variants($chr, $hFiltersChr_var2);
		$queryFilter->filter_vector_confidence_variants($chr, $hFiltersChr_var2);
		$queryFilter->filter_vector_dejavu($chr, $dejavu_2, $dejavu_ho, $test) if ($dejavu_2);
		$vector_filtered_2 = $chr->getVariantsVector->Clone();
		$h_args->{'filters_1'} = $hFiltersChr;
		$h_args->{'filters_2'} = $hFiltersChr_var2;
		$h_args->{'vector_filters_1'} = $vector_filtered;
		$h_args->{'vector_filters_2'} = $vector_filtered_2;
	}
	############## FLITER VARIANTS 
	#intersection Variants level
	$queryFilter->intersectVariantsPatients($chr);
	$queryFilter->intersectVariantsFamilies($chr);
	#exclude Variants level
	$queryFilter->excludeVariantsPatients($chr);
	$queryFilter->excludeVariantsFamilies($chr);
	
	# he ho 
	$queryFilter->filterHoVariantsPatients($chr);
	$queryFilter->filterHeVariantsPatients($chr);
	
	#GENES LEVEL ?
	
	$queryFilter->filter_genes_from_ids($chr, $hChr->{$chr->id()}, $can_use_hgmd) if ($panel_name);
	$queryFilter->filter_genes_text_search($chr, $filter_text);
	$queryFilter->filter_genes_only_genes_names($chr, $only_genes);
	
	
	
	if ($is_diff_hash_filters) { $queryFilter->filter_genes_annotations($chr, $hFiltersChr_var2); }
	else { $queryFilter->filter_genes_annotations($chr, $hFiltersChr); }
	
	$queryFilter->filter_genetics_models($chr, $model, $h_args);
	
	#intersection Genes level
	$queryFilter->intersectGenesPatients($chr);
	$queryFilter->intersectGenesFamilies($chr);
	
	#exclude Genes level
	
	$queryFilter->excludeGenesPatients($chr);
	$queryFilter->excludeGenesFamilies($chr);
	
	
	
	
	if ($project->filter_text()) {
		die();
		foreach my $patient (@{$chr->getPatients()}) {
			$patient->getVariantsVector($chr)->Intersection($patient->getVariantsVector($chr), $chr->getVariantsVector());
		}
		foreach my $family (@{$chr->getFamilies()}) {
			$family->getVariantsVector($chr)->Intersection($family->getVariantsVector($chr), $chr->getVariantsVector());
		}
	}
	print "@" unless ($export_vcf_for or $detail_project or $xls_by_regions_ho);
	my $add_filter_region = launch_filters_bed($chr, $filter_bed);
	if ($add_filter_region ne '') {
		$filter_region .= $add_filter_region;
	}
	launch_filters_region($chr, $filter_region, 1);
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER launch_filters_region - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	
	if ($chr->getVariantsVector->is_empty()) {
		print "@" unless ($export_vcf_for or $detail_project or $xls_by_regions_ho);
		next;
	}
	print "@" unless ($export_vcf_for or $detail_project or $xls_by_regions_ho);
	
#	launch_intersect_excude($chr, $h_filters_intersect_exclude) if $model =~ /compound/ or $model_2 =~ /compound/;
	launch_intersect_excude($chr, $h_filters_intersect_exclude);

	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER exclude / intersect - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	
	if ($project->typeFilters() eq 'individual') {
		$queryFilter->filter_atLeast($chr, $atLeast, $project->typeFilters(), $level_ind);
	}
	else {
		$queryFilter->filter_atLeast($chr, $atLeast, $project->typeFilters(), $level_fam);
	}
	if ($debug) { warn "\nCHR ".$chr->id()." -> AFTER filter_atLeast - nb Var: ".$chr->countThisVariants($chr->getVariantsVector()); }
	
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
	elsif ($xls_save_session or $xls_by_variants or $export_list_var_ids) {
		foreach my $v_id (@{$chr->getListVarVectorIds($chr->getVariantsVector())}) {
#			my $var_id = $chr->getVarId($v_id);
#			my $var = $chr->get_lmdb_variations("r")->get($var_id);
			my $var = $chr->getProject->returnVariants($chr->name."!".$v_id);
			$var->{project} = $project;
			$var->{buffer} = $buffer;
			$var->{vector_id} = $v_id;

			push( @lVarObj, $var );
		}
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

if ($xls_save_session or $xls_by_variants) { 
	if ($xls_by_genes) {
		saveSessionXLS($project, $hashRes, $hResumeFilters);
	}
	else {
		export_xls($project, \@lVarObj, $hResumeFilters);
	}
}
#elsif ($xls_by_variants)  { getXls_byVar($project, undef, $hResumeFilters, $xls_outfile); }
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


sub launch_filters_bed {
	my ($chr, $filter_bed) = @_;
	return unless ($filter_bed);
	my $add_filter_regions = '';
	my $bed_dir = $chr->getProject->getBedPolyQueryDir();
	foreach my $bed_file (split(',', $filter_bed)) {
		print $chr->getProject->print_dot(1);
		my $bed_path = $bed_dir.'/'.$bed_file;
		confess ("\n\nERROR: BED file doesn't exist. Die...\n\n") if not -e $bed_path;
		open (FILE, $bed_path);
		while (<FILE>) {
			my $line = $_;
			chomp ($line);
			my @lTmp = split(' ', $line);
			my $chr_id = $lTmp[0];
			$chr_id =~ s/chr//;
			next if $chr->id ne $chr_id;
			$add_filter_regions .= ' '.$chr_id.':'.$lTmp[1].':'.$lTmp[2].':99';
			#TODO: ajout add ou del bed file
		}
		close (FILE);
		print $chr->getProject->print_dot(1);
	}
	$add_filter_regions = $chr->id.':1:2:99'if $add_filter_regions eq '';
	return $add_filter_regions;
}

sub is_differents_hash_filters {
	my ($h1, $h2, $dejavu, $dejavu_2) = @_;
	return 1 if $dejavu_2;
	my $h3 = $h2;
	foreach my $cat (keys %$h1) {
		delete $h3->{$cat};
	}
	my @lDifCats = keys %$h3;
	return 1 if scalar (@lDifCats > 0);
	return;
}

sub get_infos {
	
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
sub get_hashes_filters {
	my ($filter_type_variation, $filter_type_variation_2) = @_;
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
	return ($hFiltersChr, $hFiltersChr_var2);
}



sub launch_intersect_excude {
	my ($chr, $h_filters) = @_;
#	if (not $h_filters->{filter_nbvar_regionho} or $h_filters->{filter_nbvar_regionho} == 0) {
#		$queryFilter->setExcludePatients($chr, $project->getPatientsFromListNames([split(' ', $h_filters->{filter_he})]) , 'he');
#		$queryFilter->setExcludePatients($chr, $project->getPatientsFromListNames([split(' ', $h_filters->{filter_ho})]), 'ho');
#		$queryFilter->setExcludePatients($chr, $project->getPatientsFromListNames([split(' ', $h_filters->{filter_not_patient})]), 'all');
#		$queryFilter->setExcludeFamilies($chr, $project->getFamiliesFromListNames([split(' ', $h_filters->{fam_not})]));
#	}
#	if ($h_filters->{filter_nbvar_regionho} and $h_filters->{filter_nbvar_regionho} > 0) {
#		$queryFilter->setIntersectPatient_HO_REGIONS($chr, $project->getPatientsFromListNames([split(' ', $h_filters->{filter_patient})]), $h_filters->{filter_nbvar_regionho});
#		$queryFilter->setIntersectFamily_REC_REGIONS($chr, $project->getFamiliesFromListNames([split(' ', $h_filters->{fam_and})]), $h_filters->{filter_nbvar_regionho});
#		$queryFilter->setExcludePatient_HO_REGIONS($chr, $project->getPatientsFromListNames([split(' ', $h_filters->{filter_not_patient})]), $h_filters->{filter_nbvar_regionho});
#		$queryFilter->setExcludeFamily_HO_REGIONS($chr, $project->getFamiliesFromListNames([split(' ', $h_filters->{fam_not})]), $h_filters->{filter_nbvar_regionho});
#	}
#	$queryFilter->setIntersectPatients($chr, $project->getPatientsFromListNames([split(' ', $h_filters->{filter_patient})]));
#	$queryFilter->setIntersectFamilies($chr, $project->getFamiliesFromListNames([split(' ', $h_filters->{fam_and})]));
}

sub launch_filters_region {
	my ($chr, $filter_region, $first_launch) = @_;
	return unless ($filter_region);
	my @lFilters;
	my $isIntersect;
	my $hFilters;
	foreach my $this_filter (split(" ", $filter_region)) {
		print $chr->getProject->print_dot(1000);
		next if (exists $hFilters->{$this_filter});
		$hFilters->{$this_filter} = undef;
		my ($chrId, $start, $end, $include) = split(":", $this_filter);
		next unless ($chrId eq $chr->id());
		$chr->check_each_var_filter_region($this_filter, $first_launch);
		$isIntersect = 1 if ($include eq '0' or $include eq '99');
	}
	$chr->variants_regions_add();
	if ($chr->variants_regions_add->is_empty()) {
		if ($isIntersect) {
			foreach my $patient (@{$chr->getPatients()}) {
				print $chr->getProject->print_dot(5);
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




my $hpatients_genes = {};
my $hfamillies_genes = {};

sub launchStatsProjectAll {
	my $hash_stats;
	warn "\n\n" if ($debug);
		my $t =time;
	$hash_stats->{genes} 		= launchStatsProjectAll_genes();
	warn "-".abs(time -$t) if ($debug);	my $t =time;
	$hash_stats->{genes} 		= launchStatsProjectAll_genes();
	warn "-".abs(time -$t) if ($debug);
	warn "\n# launchStatsProjectAll_chromosomes" if ($debug);
	$hash_stats->{chromosomes} 	= launchStatsProjectAll_chromosomes($hash_stats->{genes});
	warn "\n# launchStatsProjectAll_genes" if ($debug);

	$t =time;
	warn "\n# launchStatsProjectAll_families" if ($debug);
	$hash_stats->{families} 	= launchStatsProjectAll_families();
	warn "\n# launchStatsProjectAll_patients" if ($debug);
	warn "- fam".abs(time -$t) if ($debug);
	$t =time;
	$hash_stats->{patients} 	= launchStatsProjectAll_patients();
	warn "-".abs(time -$t) if ($debug);
	$t =time;
	warn "\n# stats_region" if ($debug);
	$hash_stats->{regions} 	    = $project->stats_region();

	if ($filter_bed and $filter_chromosome =~ /X/) {
		my @lBed;
		foreach my $bed_file (split(',', $filter_bed)) {
			my $h_bed;
			$h_bed = {
	            'deletion' => 0,
	            'chromosome' => 'file:'.$bed_file,
	            'start' => 1,
	            'id' => 'file:'.$bed_file.':1:2:0',
	            'include' => '0',
	            'insertion' => 0,
	            'substitution' => 0,
	            'genes' => 0,
	            'end' => 1
			};
			unshift(@{$hash_stats->{regions}}, $h_bed);
		}
	}
	warn "\n# launchStatsProjectAll_groups" if ($debug);
	$hash_stats->{groups}		= launchStatsProjectAll_groups() if ($project->isSomaticStudy());
	warn "-".abs(time -$t) if ($debug);
	$t =time;
	warn "\n# json_all_genes_name" if ($debug);
	$hash_stats->{genes_name}	= json_all_genes_name();
	warn "-".abs(time -$t) if ($debug);
	$t =time;
	warn "\n# stats_regions_ho_rec" if ($debug);
#	$hash_stats->{regions_ho_rec} = launchStatsPojectRegionsHoRec() if ($filter_nbvar_regionho);
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
my %vector_buffer;





sub launchStatsProjectAll_genes_fork {
	my @lStats;
	my $fork = 3;
	my $t = time;
	my $pm = new Parallel::ForkManager($fork);
	my $i = 0;
		my $nb =0;
		my $nbg = 0;
	$pm->run_on_finish(
	sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $res) = @_;
		unless (defined($res) or $exit_code > 0) {
			print qq|No message received from child process $exit_code $pid!\n|;
			die();
			return;
		}
		
		push(@lStats,@{$res->{data}});
		warn "end";
		foreach my $k (keys %{$res->{genes}}){
			$h_all_genes_name->{$k} = $res->{genes}->{$k};
		}
		warn "finish;"
	}
);	
		my $process;
		
	
		
		$buffer->getHashTransIdWithCaptureDiag();
		
	foreach my $chr_id (split(',', $filter_chromosome)) {
		my $chr = $project->getChromosome($chr_id);
		$vector_buffer{$chr_id}->{substitution} = $chr->getVectorSubstitutions;
		$vector_buffer{$chr_id}->{insertion} =  $chr->getVectorInsertions ;
		$vector_buffer{$chr_id}->{deletion} =  $chr->getVectorDeletions;
		$vector_buffer{$chr_id}->{indel} =  $chr->getVectorInsertions + $chr->getVectorDeletions unless exists $vector_buffer{$chr_id}{indel} ;
		$vector_buffer{$chr_id}->{cnv} =  $chr->getVectorLargeDuplications() + $chr->getVectorLargeDeletions() unless exists $vector_buffer{$chr_id}{cnv} ;	
		
		next if ($chr->not_used());
		#my @genes = grep {not ($_->getCurrentVector->is_empty())} @{$chr->getGenesFromVector($chr->getVariantsVector())};
		warn "genes ";
	
		my @genes =   @{$chr->getGenes()};
		$project->disconnect;
		map{$_->enum} @genes;
		delete $project->{rocks};
		 #$chr->flush_rocks_vector();
		# $chr->flush_rocks_vector();
		 warn "purge";
		my $nb        = int( scalar(@genes) / ($fork) +1 );
		my $iter      = natatime( $nb,  @genes );
		$t = time;
		while ( my @tmp = $iter->() ) {
			$process ++;
		#	my $pid = $pm->start and next;
			
			warn "start $process";
			my $h_all_genes_name = {};
		foreach my $gene (@tmp) {	
		#foreach my $gene (@{$chr->getGenesFromVector($chr->getVariantsVector())}) {
			$nb ++;
			next if $gene->getCurrentVector->is_empty();
			warn ref($gene) if ($debug);
			$project->print_dot(50);
			next if not $gene->getCurrentVector();
			warn $gene->external_name() if ($debug);
			my $hStats = launchStatsGene($gene);
			
			if ($hStats) {
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
		}
			warn "end ".$process;
	#	$project->disconnect();
		
				warn "end ".$process;
	#	$pm->finish( 0, {data=>\@lStats,genes=>$h_all_genes_name} );
		}
		$pm->wait_all_children();
		warn '-> all genes DONE';
		warn '-> all genes DONE' if ($debug);
	}
	warn abs(time - $t);
	return \@lStats;
}


sub launchStatsProjectAll_genes {
	my @lStats;
	
	my $i = 0;
		my $nbg = 0;
		$buffer->getHashTransIdWithCaptureDiag();
	foreach my $chr_id (split(',', $filter_chromosome)) {
		my $chr = $project->getChromosome($chr_id);
		$vector_buffer{$chr_id}{substitution} = $chr->getVectorSubstitutions;
		$vector_buffer{$chr_id}{insertion} =  $chr->getVectorInsertions ;
		$vector_buffer{$chr_id}{deletion} =  $chr->getVectorDeletions;
		$vector_buffer{$chr_id}{indel} =  $chr->getVectorInsertions + $chr->getVectorDeletions ;
		$vector_buffer{$chr_id}{cnv} =  $chr->getVectorLargeDuplications() + $chr->getVectorLargeDeletions();
		
		next if ($chr->not_used());
		foreach my $gene (@{$chr->getGenesFromVector($chr->getVariantsVector())}) {
		#foreach my $gene (@{$chr->getGenes}) {
#			my $vchr = $gene->return_compact_vector( $vvv);
#			warn $vchr->Norm()." ** ".$vvv->Norm  if $debug;
#			$vsmall &= $vchr;
			next if $gene->getCurrentCompactVector->is_empty();
			$project->print_dot(50);
#			my $v_to_enum = $gene->_getVectorOrigin->to_Enum();
			#TODO: a enlever une fois vector ok
#			if ($v_to_enum =~ /0-/) {
#				next;
#			}
			$nbg ++;
			warn $gene->external_name() if ($debug);
			my $hStats = launchStatsGene($gene);
			if ($hStats) {
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
		}
		warn '-> all genes DONE' if ($debug);
	}
	return \@lStats;
}
		
sub launchStatsProjectAll_chromosomes {
	my ($stat_genes) = @_;
	my $nb_genes = scalar(@$stat_genes);
	my @lStats;
	foreach my $chr_id (split(',', $filter_chromosome)) {
		my $chr = $project->getChromosome($chr_id);
		
		next if ($chr->not_used());
		my (@ag) = grep {$_->{chromosome} eq $chr->name} @$stat_genes;
		my $nb_genes = scalar(@ag);
		my $hashChr = launchStatsChr($chr,$nb_genes);
		if ($hashChr->{'variations'} == 0) { $hashChr = launchStatsChr_null($chr); }
		else { push(@lStats, $hashChr); }
	}
	if (scalar(@lStats) == 0 and $without_stats) {
		my $chr = $project->getChromosome('1');
		my $hashChr = launchStatsChr_null($chr);
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

sub launchStatsProjectAll_groups {
	my @lStats;
	my $hash_global_stats_fam;
	foreach my $chr_id (split(',', $filter_chromosome)) {
		my $chr = $project->getChromosome($chr_id);
		if ($chr->not_used()) { $chr->getVariantsVector->Empty(); }
		foreach my $group (@{$project->getSomaticGroups()}) {
			$project->print_dot(100);
			my $famName = $group->name();
			my $hStatsFam = launchStatsFamily($group, $chr);
			my @lCat1 = ('model', 'include', 'name', 'nb', 'id');
			foreach my $category (@lCat1) {
				$hash_global_stats_fam->{$famName}->{$category} = $hStatsFam->{$category};
			}
			my @lCat2 = ('homozygote', 'heterozygote', 'substitution', 'insertion', 'deletion', 'genes', 'cnv');
			foreach my $category (@lCat2) {
				unless (exists $hash_global_stats_fam->{$famName}->{$category}) { $hash_global_stats_fam->{$famName}->{$category} = 0; }
				$hash_global_stats_fam->{$famName}->{$category} += $hStatsFam->{$category};
			}
		}
	}
	foreach my $famName (sort keys %{$hash_global_stats_fam}) {
		push(@lStats, $hash_global_stats_fam->{$famName});
	}
	return \@lStats;
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
				$hash_global_stats_fam->{$famName}->{$category} += $hStatsFam->{$category};
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
#			next if (not exists $chr->patients_categories->{$patient->name()} and int($project->annotation_version()) < int($project->annotation_version_current()));
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

#sub launchStatsProjectAll_groups {
#	my @lStats;
#	my $hash_global_stats_groups;
#	foreach my $chr_id (split(',', $filter_chromosome)) {
#		my $chr = $project->getChromosome($chr_id);
#		warn 'CHR '.$chr->id().' -> '.$chr->not_used() if ($debug);
#		my $hashGroups = $chr->stats_groups();
#		foreach my $famName (keys %{$hashGroups}) {
#			my @lCat1 = ('model', 'include', 'name', 'nb', 'id');
#			foreach my $category (@lCat1) {
#				$hash_global_stats_groups->{$famName}->{$category} = $hashGroups->{$famName}->{$category};
#			}
#			my @lCat2 = ('homozygote', 'heterozygote', 'substitution', 'insertion', 'deletion');#, 'genes');
#			foreach my $category (@lCat2) {
#				$hash_global_stats_groups->{$famName}->{$category} += $hashGroups->{$famName}->{$category};
#			}
#		}
#	}
#	foreach my $groupName (sort keys %{$hash_global_stats_groups}) {
#		push(@lStats, $hash_global_stats_groups->{$groupName});
#	}
#	return \@lStats;
#}

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
	my ($chr,$nb_gene) = @_; ;
	my $hash;
	my $name = $chr->id();
	if ($name eq 'X')     { $name = 23; }
	elsif ($name eq 'Y')  { $name = 24; }
	elsif ($name eq 'MT') { $name = 25; }
	$hash->{id}         = $chr->id();
	$hash->{name}       = int($name);
	$hash->{genes}      = $nb_gene;
	$hash->{variations}	= 0;
	
	my $v_sub = $chr->vector_global_categories('substitution')->Clone();
	$v_sub->Intersection($v_sub, $chr->getVariantsVector());
	my $v_ins = $chr->vector_global_categories('insertion')->Clone();
	$v_ins->Intersection($v_ins, $chr->getVariantsVector());
	my $v_del = $chr->vector_global_categories('deletion')->Clone();
	$v_del->Intersection($v_del, $chr->getVariantsVector());
	my $v_l_dup = $chr->vector_global_categories('large_duplication')->Clone();
	$v_l_dup->Intersection($v_l_dup, $chr->getVariantsVector());
	my $v_l_del = $chr->vector_global_categories('large_deletion')->Clone();
	$v_l_del->Intersection($v_l_del, $chr->getVariantsVector());

	$hash->{substitutions} = $chr->countThisVariants($v_sub);
	$hash->{insertions} = $chr->countThisVariants($v_ins);
	$hash->{deletions} = $chr->countThisVariants($v_del);
	$hash->{cnvs} = $chr->countThisVariants($v_l_del);
	$hash->{cnvs} = $chr->countThisVariants($v_l_dup);
	$hash->{variations}	= $chr->countThisVariants( $chr->getVariantsVector() );
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

my %global_hcat;
sub launchStatsGene {
	my ($gene) = @_;
	
	my $vv = $gene->getCurrentVector() & $gene->getChromosome->getVariantsVector();
	$gene->setCurrentVector($vv);
	my $t1 = $gene->compact_vector;
	my $va = $gene->getCurrentCompactVector();
	return if $va->is_empty;
	
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
		my $chr = $gene->getChromosome();
	$hashStats->{'name'}  = $gene->name();
	$hashStats->{'xref'}  = $gene->external_name();
	$hashStats->{'start'}  = $gene->start();
	$hashStats->{'end'}  = $gene->end();
	$hashStats->{chromosome} = $chr->name;
	#Interval_Copy
	#$hashStats->{'description'} = $gene->phenotypes();
	
	$hashStats->{'phenotype'} = $gene->polyquery_phenotypes();
	
	$hashStats->{'description'} = $gene->description() if ($gene->description());
	#my ($pheno,$nb_other_terms) = $gene->polyviewer_phentotypes();
	#$hashStats->{'description'} = "$pheno + $nb_other_terms terms";

	#
	my ($h_vector_type, $hStats, $h_var_ids, $h_models, $hAllPat);
	
	my $vv = $gene->getCurrentVector() & $gene->getChromosome->getVariantsVector();
	$gene->setCurrentVector($vv);
	
	my $t1 = $gene->compact_vector;
	my $va = $gene->getCurrentCompactVector();
	return if $va->is_empty;
	my $djvf;
	foreach my $patient (@{$chr->getProject->getPatients()}) {
		
		my $vv = $va & $gene->getCompactVectorPatient($patient);
		next if $vv->is_empty;
		$hAllPat->{$patient->name()} = undef;
		$hStats->{variants}->{all} = $chr->countThisVariants($vv);
		$hStats->{patients}->{all}->{$patient->id()} ++;
		$hStats->{fam}->{all}->{$patient->getFamily->id()} ++;
		$hpatients_genes->{$patient->id()}++;
		
	}
	
	
	foreach my $f   (keys %{$hStats->{fam}->{all}}){
			$hfamillies_genes->{$f} ++;
		
	}
	
	$h_vector_type->{low}  = $gene->getCompactVector("low") &  $va  ;
	
	$hStats->{variants}->{low} = $chr->countThisVariants($h_vector_type->{low}) ;
	
	$h_vector_type->{medium} =  $gene->getCompactVector("medium") &  $va  ;
	$hStats->{variants}->{medium} = $chr->countThisVariants($h_vector_type->{medium}) ;
	
	$h_vector_type->{high} =  $gene->getCompactVector("high") &  $va  ;
	$hStats->{variants}->{high} = $chr->countThisVariants($h_vector_type->{high}) ;
	$h_vector_type->{all} = $va;
	$hStats->{variants}->{all} = $chr->countThisVariants($h_vector_type->{all}) ;
	
	$h_vector_type->{substitution} = $gene->getCompactVector("substitution") & $va;
	$hStats->{variants}->{substitution} = $chr->countThisVariants($h_vector_type->{substitution}) ;
	
	$h_vector_type->{indel} = $gene->getCompactVector("indel") & $va;
	$hStats->{variants}->{indel} = $chr->countThisVariants($h_vector_type->{indel}) ;
	
	$h_vector_type->{cnv} =$gene->getCompactVector("cnv") & $va;
	$hStats->{variants}->{cnv} = $chr->countThisVariants($h_vector_type->{cnv}) ;

	$hashStats->{'vector_ids'}	= join(',', @{$gene->getIdsBitOn($gene->getCurrentVector())});
	if ($can_use_hgmd and $gene->hgmd()) {
		$hashStats->{'has_hgmd'} = 2;
		$hashStats->{'has_hgmd'} .= ';'.$gene->external_name();
		$hashStats->{'query_score'} += $project->buffer->config->{score_query_gene}->{has_hgmd_dm};
	}
	my $hGetFam;
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
#	if ($filter_nbvar_regionho > 1){
#		confess();
#		$hashStats->{'region_ho'} = undef;
#		$hashStats->{'region_ho_all'} = undef;
#		$hashStats->{'region_rec'} = undef;
#		$hashStats->{'region_rec_all'} = undef;
#		my (@lObj, $region_method, $stat_method);
#		if ($project->typeFilters() eq 'individual'){
#			foreach my $patient (@{$gene->chromosome->getPatients()}) {
#				my $vec_tmp = $patient->getRegionHo($gene->chromosome(), $filter_nbvar_regionho, $filter_regionho_sub_only)->Clone();
#				$vec_tmp->Intersection($gene->getCurrentVector(), $patient->getRegionHo($gene->chromosome(), $filter_nbvar_regionho, $filter_regionho_sub_only));
#				unless ($vec_tmp->is_empty()) {
#					push(@{$hashStats->{'region_ho_all'}}, $project->hash_patients_name->{$patient->name()}->{id});
#					$vec_tmp->Intersection($vec_tmp, $patient->getCurrentVector($gene->chromosome()));
#					unless ($vec_tmp->is_empty()) {
#						push(@{$hashStats->{'region_ho'}}, $project->hash_patients_name->{$patient->name()}->{id});
#					}
#				}
#			}
#		}
#		elsif($project->typeFilters() eq 'familial'){
#			foreach my $family (@{$gene->chromosome->getFamilies()}) {
#				my $vec_tmp = $family->getVectorRegionRec($gene->chromosome(), $filter_nbvar_regionho, $filter_regionho_sub_only)->Clone();
#				$vec_tmp->Intersection($gene->getCurrentVector(), $family->getVectorRegionRec($gene->chromosome(), $filter_nbvar_regionho, $filter_regionho_sub_only));
#				unless ($vec_tmp->is_empty()) {
#					push(@{$hashStats->{'region_rec_all'}}, $family->name());
#					$vec_tmp->Intersection($vec_tmp, $family->getHo($gene->chromosome()));
#					unless ($vec_tmp->is_empty()) {
#						push(@{$hashStats->{'region_rec'}}, $family->name());
#					}
#				}
#			}
#			
#		}
#		if ($hashStats->{'region_ho'}) {
#			$hashStats->{'region_ho'} = scalar @{$hashStats->{'region_ho'}}.';'.join('|' , @{$hashStats->{'region_ho'}});
#		}
#		if ($hashStats->{'region_ho_all'}) {
#			$hashStats->{'region_ho_all'} = scalar @{$hashStats->{'region_ho_all'}}.';'.join('|' , @{$hashStats->{'region_ho_all'}});
#		}
#		if ($hashStats->{'region_rec'}) {
#			$hashStats->{'region_rec'} = scalar @{$hashStats->{'region_rec'}}.';'.join('|' , @{$hashStats->{'region_rec'}});
#		}
#		if ($hashStats->{'region_rec_all'}) {
#			$hashStats->{'region_rec_all'} = scalar @{$hashStats->{'region_rec_all'}}.';'.join('|' , @{$hashStats->{'region_rec_all'}});
#		}
#	}

	return $hashStats;

}

sub launchStatsPatient {
	my ($patient, $chr) = @_;
	confess("\n\nERROR: GenBoPatientCache->stats() method need a GenBoChromosomeCache and a category object in argument. Die.\n\n") unless(ref($chr) eq 'GenBoChromosomeCache');
	return launchStatsPatient_null($patient) if ($chr->getVariantsVector->is_empty());
	#2018_06_28: Choix de patrick, on met artificiellement toutes les valeures a 0 si on exclu des regions ho dun patient (comme in the attic)
	return launchStatsPatient_null($patient) if ($patient->excluded() eq 'ho_reg');
	return launchStatsPatient_null($patient) if ($patient->in_the_attic);
	my $hStats;
	my $patName = $patient->name();
	$hStats->{name} 		= $patName;
	$hStats->{id} 			= int($project->hash_patients_name->{$patient->name()}->{'id'});
	$hStats->{bam} 			= "";
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
	
	
	my $v_all = $patient->getVectorOrigin($chr) & $chr->getVariantsVector();

	my $v_he = $patient->getVectorOriginHe($chr) & $v_all;
	my $v_ho = $patient->getVectorOriginHo($chr) & $v_all;
	
	my $v_substitution = $chr->getNewVector();
	my $v_insertion = $chr->getNewVector();
	my $v_deletion = $chr->getNewVector();
	my $v_cnvs = $chr->getNewVector();
	
	$v_substitution = $vector_buffer{$chr->id}{substitution} & $v_all if $vector_buffer{$chr->id}{substitution};
	$v_insertion = $vector_buffer{$chr->id}{insertion} & $v_all if $vector_buffer{$chr->id}{insertion};
	$v_deletion = $vector_buffer{$chr->id}{deletion} & $v_all if $vector_buffer{$chr->id}{deletion};
	$v_cnvs = $vector_buffer{$chr->id}{cnv} & $v_all if $vector_buffer{$chr->id}{cnv};
	
	my $nb_all = $chr->countThisVariants($v_all);
	$hStats->{variations} = $nb_all;
	$hStats->{all_variations} = $nb_all;
	
	$hStats->{homozygote} =  $chr->countThisVariants($v_ho);
	$hStats->{heterozygote} = $chr->countThisVariants($v_he);
	$hStats->{substitutions} = $chr->countThisVariants($v_substitution);
	$hStats->{insertions} = $chr->countThisVariants($v_insertion);
	$hStats->{deletions} = $chr->countThisVariants($v_deletion);
	$hStats->{cnvs} = $chr->countThisVariants($v_cnvs);
	$hStats->{found} = '';
	$hStats->{genes} = $hpatients_genes->{$patient->id};
	if ($hStats->{all_variations} > 0) {
		$hStats->{found} = 'yes';
	}

	#$hStats->{genes} = $patient->countGenes($chr);
	if ($xls_by_variants or $xls_by_variants) { $hStats->{genes} = $patient->countGenes($chr); }
	
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

	my $v_all = $chr->getNewVector();
	my $v_ho = $chr->getNewVector();
	my $v_he = $chr->getNewVector();
	my $v_substitution = $chr->vector_global_categories('substitution')->Clone();
	my $v_insertion = $chr->vector_global_categories('insertion')->Clone();
	my $v_deletion = $chr->vector_global_categories('deletion')->Clone();
	my $v_large_deletion = $chr->vector_global_categories('large_deletion')->Clone();
	my $v_large_duplication = $chr->vector_global_categories('large_duplication')->Clone();
	foreach my $patient (@{$family->getPatients()}) {
		$v_he += $patient->getVectorOriginHe($chr);
		$v_ho += $patient->getVectorOriginHo($chr);
		$v_all += $v_he;
		$v_all += $v_ho;
	}
	$v_all->Intersection($v_all, $chr->getVariantsVector());
	$v_he->Intersection($v_he, $v_all);
	$v_ho->Intersection($v_ho, $v_all);
	$v_substitution->Intersection($v_substitution, $v_all);
	$v_insertion->Intersection($v_insertion, $v_all);
	$v_deletion->Intersection($v_deletion, $v_all);
	$v_large_deletion->Intersection($v_large_deletion, $v_all);
	$v_large_duplication->Intersection($v_large_duplication, $v_all);

	$hash->{substitution} = $chr->countThisVariants($v_substitution);
	$hash->{insertion} = $chr->countThisVariants($v_insertion);
	$hash->{deletion} = $chr->countThisVariants($v_deletion);
	$hash->{cnv} = $chr->countThisVariants($v_large_deletion);
	$hash->{cnv} += $chr->countThisVariants($v_large_duplication);
	$hash->{heterozygote} = $chr->countThisVariants($v_he);
	$hash->{homozygote} = $chr->countThisVariants($v_ho);
	$hash->{genes} = $hfamillies_genes->{$family->id};
#	my $nb = 0;
#	foreach my $gene (@{$chr->getGenes()}) {
#		next unless $gene->getVectorOrigin();
#		next if ($gene->is_intergenic());
#		my $v_gene = $gene->getCurrentVector->Clone();
#		next if ($v_gene->is_empty());
#		$v_gene->Intersection($v_gene, $v_all);
#		$nb++ if (not $v_gene->is_empty()) ;
#	}
#	$hash->{genes} = $nb;
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
#	$chr->cache_lmdb_variations->close();
#	$chr->get_lmdb_patients->close();
#	$chr->get_lmdb_categories->close();
#	$chr->get_lmdb_genes->close() if ($chr->get_lmdb_genes());
#	$chr->get_lmdb_variations->close();
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
		elsif (-d $Bin.'/../../../polyweb/') { $path_polyweb = $Bin.'/../../../polyweb/'; }
		elsif (-d $Bin.'/../../../PolyWeb/') { $path_polyweb = $Bin.'/../../../PolyWeb/'; }
		else { warn "\n\n$Bin\n\n"; die; }
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
			next if ($cat eq 'large_insertion');
			next if ($cat eq 'spliceAI_high');
			next if ($cat eq 'spliceAI_medium');
			next if ($cat eq 'predicted_splice_site');
			next if ($cat eq 'junction');
			unless (exists $chr->global_categories->{$cat}) {
				warn "\n\n";
				warn Dumper $chr->global_categories();
				confess("\n\nERROR: [CHR".$chr->id()."] global_categories->{$cat} not found. Pb cache ?? Die.\n\n");
			}
		}
	}
	return 1;
}

#sub loadSessionsXLS {
#	my ($project, $sid) = @_;
#	my $tmp_dir = $project->getTmpDir();
#    my $session = new CGI::Session(undef, $sid, {Directory=>$tmp_dir});
#    my $h = thaw(decompress $session->param('hash_xls'));
#    $session->delete();
#    if ($xls_by_variants) { getXls_byVar($project, $h, $hResumeFilters); }
#    elsif ($xls_by_genes) { getXls_byGenes($project, $h, $hResumeFilters); }
#	exit(0);
#}

sub saveSessionXLS {
	my ($project, $hashRes, $hResumeFilters) = @_;
	
	die;
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

sub export_xls {
	my ($project, $lVar, $hResumeFilters) = @_;
	my $xls_export = new xls_export();
	$xls_export->title_page('PolyQuery_'.$project->name().'.xls');
	$xls_export->store_variants_infos($lVar, $project, $project->getPatients());
	if ($xls_save_session) {
		my $session_id = $xls_export->save();
	    print '@@@';
		print "\",\"session_id\":\"";
	    print $session_id;
	    print "\"}";
		exit(0);
	}
	elsif ($xls_outfile) {
		$xls_export->output_dir($xls_outfile);
	}
	create_xls_variants($xls_export, $hResumeFilters);
}

sub loadSessionsXLS {
	my ($project, $sid, $hResumeFilters) = @_;
	if ($xls_by_genes) {
		return loadSessionsXLS_byGenes($project, $sid, $hResumeFilters);
	}
	my $xls_export = new xls_export();
	$xls_export->load($sid);
	
	create_xls_variants($xls_export, $hResumeFilters);
}

sub create_xls_variants {
	my ($xls_export, $hResumeFilters) = @_;
	my ($list_datas_annotations) = $xls_export->prepare_generic_datas_variants();
	my ($list_datas_annotations_cnvs) = $xls_export->prepare_generic_datas_cnvs();
	if ($detail_project) {
		$xls_export->add_page_merged('Variants Merged', $xls_export->list_generic_header(), $list_datas_annotations);
	}
	$xls_export->add_page('Variants Not Merged', $xls_export->list_generic_header(), $list_datas_annotations);
	if (scalar @$list_datas_annotations_cnvs > 0) {
		$xls_export->add_page('Cnvs', $xls_export->list_generic_header_cnvs(), $list_datas_annotations_cnvs);
	}
	$xls_export->export();
	my $xls_resume = $xls_export->workbook->add_worksheet('FILTERS USED');
	writeResumeFiltersXls($xls_resume, undef, $hResumeFilters);
	exit(0);
}


sub loadSessionsXLS_byGenes {
	my ($project, $sid) = @_;
	my $tmp_dir = $project->getTmpDir();
    my $session = new CGI::Session(undef, $sid, {Directory=>$tmp_dir});
    my $h = thaw(decompress $session->param('hash_xls'));
    $session->delete();
    if ($xls_by_variants) { getXls_byVar($project, $h, $hResumeFilters); }
    elsif ($xls_by_genes) { getXls_byGenes($project, $h, $hResumeFilters); }
	exit(0);
}
