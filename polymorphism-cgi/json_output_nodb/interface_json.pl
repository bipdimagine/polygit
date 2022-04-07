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
use lib "$Bin/../packages/validation_variation"; 
use GBuffer;
use export_data;
use JSON;
use VcfMerge;
use GenBoNoSql;
use Spreadsheet::WriteExcel;
use Compress::Snappy;
use Storable qw(store retrieve freeze dclone thaw);
use POSIX qw(strftime);
use CGI::Session;
use html; 


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
my $filter_diseases			= $cgi->param('filter_diseases');
my $filter_genes_intersect	= $cgi->param('filter_genes_intersect');
my $level_fam 				= $cgi->param('level_fam');
my $level_ind 				= $cgi->param('level_ind');
my $atLeast 				= $cgi->param('nb');
my $typeFilters 			= $cgi->param('mode');
my $model 					= $cgi->param('model');
my $projectName 			= $cgi->param('project');
my $dejavu					= $cgi->param('dejavu');
my $dejavu_2				= $cgi->param('dejavu_2');
my $dejavu_ho				= $cgi->param('dejavu_ho');
my $getPatientsCgi			= $cgi->param('get_patients');
my $xls_by_genes			= $cgi->param('xls_by_genes');
my $xls_by_variants			= $cgi->param('xls_by_variants');
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
my $mosaique_min_p			= $cgi->param('mosaique_min_p');


#	my $buffer = new GBuffer;
#	my $project = $buffer->newProjectCache( -name => $projectName );
#	my $chr = $project->getChromosome($filter_chromosome);
#	foreach my $pat (@{$project->getPatients()}) {
#		next unless ($pat->name() eq '16TL00381');
#		#next unless ($pat->name() eq 'GUI_REF');
#		warn "\n";
#		warn $pat->name();
#		warn "\n\n";
#		warn '--------------------';
#		warn "\n\n";
#		warn $pat->is_cache_ok();
#		warn "\n\n";
#		warn '--------------------';
#		warn "\n\n";
#		warn $pat->is_cache_ok_2();
#		die;
#	}

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

my $buffer = new GBuffer;
my $project;
if ($test) {
	$project = $buffer->newProjectCache( -name 			=> $projectName,
										 -test	 		=> '1',
										 -cache 		=> '1',
										 -typeFilters 	=> $typeFilters, );
}
else {
	$project = $buffer->newProjectCache( -name 			=> $projectName,
										 -cache 		=> '1',
										 -typeFilters 	=> $typeFilters, );
}

if ($check_genes) {
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

if ($xls_by_variants) { $project->isXlsOutput('variants'); }
if ($xls_by_genes)    { $project->isXlsOutput('genes'); }

if ($filter_genes_intersect) {
	my $h;
	foreach my $g_id (split(',', $filter_genes_intersect)) {
		$h->{$g_id} = undef;
	}
	$project->filter_genes_intersect($h);
}

if ($filter_nbvar_regionho) {
	$project->filter_nbvar_regionho($filter_nbvar_regionho);
}

if ($get_bundles) {
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
		$h_capture->{name} = $capture.' <span id="span_capture_'.$capture.'"/><font color="green">[' . $nb_captures . ' ' . $text_nb_capt . ' - '.scalar(keys %{$hGenes_capture->{$capture}}).' genes]</font></span>';
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

if ($bundle_id) {
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
		$hFiltersChr->{'large_insertion'} = undef;
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
			$hFiltersChr_var2->{'large_insertion'} = undef;
		}
		elsif ($filter_name eq 'upstream_downstream') {
			$hFiltersChr->{'upstream'} = undef;
			$hFiltersChr->{'downstream'} = undef;
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
$project->filter_text($filter_text) if ($filter_text);
$project->nb_max_dejavu($dejavu) if ($dejavu);
$project->nb_max_dejavu_2($dejavu_2) if ($dejavu_2);
$project->dejavu_ho(1) if ($dejavu_ho);

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

if ($xls_save_session) { print $cgi->header('text/html'); }
else {
	if ($export_vcf_for) {}
	elsif ($detail_project) {}
	elsif (not $test) {
		print $cgi->header('text/json-comment-filtered');
		print "{\"progress\":\".";
	}
}

my $hAllIds;
my @lVarObj;
my $nb_chr = 0;
foreach my $chr_id (split(',', $filter_chromosome)) {
	my $chr = $project->getChromosome($chr_id);
	$nb_chr++;
	$chr->project->nb_at_least($atLeast) if ($atLeast and ($level_fam eq 'gene' or $level_ind eq 'gene'));
	if ($chr->not_used()) {
		print "@" unless ($export_vcf_for or $detail_project);
		next;
	}
	if ($mosaique_min_p) { $chr->mosaique_min_p($mosaique_min_p); }
	if ($export_vcf_for) {}
	elsif ($detail_project) {}
	else { $project->cgi_object(1); }
	unless (check_region_filter($chr, $filter_region)) {
		print "@" unless ($export_vcf_for or $detail_project);
		next;
	}
	if ($xls_save_session) {
		print "@" unless ($export_vcf_for or $detail_project);
	}
	# raccourci double clic famille, intersect sur lui meme
	if ($fam_and) {
		my @lFamAnd = split(' ', $fam_and);
		if (scalar(@lFamAnd) == 1) {
			my $h_attic;
			foreach my $fam (@{$chr->getFamilies()}) {
				next if ($fam->name() eq $lFamAnd[0]);
				foreach my $pat (@{$fam->getPatients()}) {
					$h_attic->{$pat->name()} = undef;
				}
			}
			if ($filter_attic) {
				foreach my $pat_name (split(' ', $filter_attic)) { $h_attic->{$pat_name} = undef; }
			}
			$filter_attic = join(' ', keys %$h_attic);
			$fam_and = undef;
		}
	}
	
	# Pour XLS, besoin de savoir si un patient possede ou non le variant, meme s il ne passe pas un filtre ou un modele pour rester coherent a la 2e page par gene
	if ($xls_by_variants) { $chr->save_model_variants_all_patients('for_xls'); }
	
	$chr->setInTheAttic($project->getPatientsFromListNames([split(' ', $filter_attic)]));
	print "@" unless ($export_vcf_for or $detail_project);
	$chr->doPolyQueryFilters($hFiltersChr, $hFiltersChr_var2);
	print "@" unless ($export_vcf_for or $detail_project);
	launch_filters_region($chr, $filter_region, 1);
	if ($chr->getVariantsVector->is_empty()) {
		print "@" unless ($export_vcf_for or $detail_project);
		$chr->purge();
		next;
	}
	print "@" unless ($export_vcf_for or $detail_project);
	$chr->setExcludePatient($project->getPatientsFromListNames([split(' ', $filter_he)]) , 'he');
	$chr->setExcludePatient($project->getPatientsFromListNames([split(' ', $filter_ho)]), 'ho');
	$chr->setExcludePatient($project->getPatientsFromListNames([split(' ', $filter_not_patient)]), 'all');
	$chr->setExcludeFamily($project->getFamiliesFromListNames([split(' ', $fam_not)]));
	$chr->setIntersectPatient($project->getPatientsFromListNames([split(' ', $filter_patient)]));
	$chr->setIntersectPatient_GENES($project->getPatientsFromListNames([split(' ', $filter_patient)]));
	$chr->setIntersectFamily($project->getFamiliesFromListNames([split(' ', $fam_and)]));
	$chr->setIntersectFamily_GENES($project->getFamiliesFromListNames([split(' ', $fam_and)]));
	launch_filters_atLeast($chr);
	launch_filters_region($chr, $filter_region);
	check_variants_regions_exclude($chr);
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
	}
	else {
		print "@" unless ($export_vcf_for or $detail_project);
		$chr->launchStats();
	}
	$chr->purge() unless ($gene_atlas_view);
	print "@" unless ($export_vcf_for or $detail_project);
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
$hashRes->{'items'} = $project->stats();

if ($gene_atlas_view) {
	geneAtlasView($project);
}

if    ($xls_save_session) { saveSessionXLS($project, $hashRes, $hResumeFilters); }
elsif ($xls_by_variants)  { getXls_byVar($project, undef, $hResumeFilters, $xls_outfile); }
elsif ($xls_by_genes)     {
	my @lStats = @{ $project->stats_genes() };
	my $h;
	$h->{by_genes} = \@lStats;
	getXls_byGenes($project, $h, $hResumeFilters, $xls_outfile);
}
else                      { printJson($hashRes, $test); } 



###### GLOBAL METHODS #####




sub geneAtlasView {
	my ($project) = @_;
	
	
	$project->geneAtlasView();
	
	exit(0);
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
		print ".";
		my @lStats = @{ $project->stats_genes() };
		if (scalar(@lStats) > 0) {
			print ".";
			$h->{by_genes} = \@lStats;
			$ok = 1;
		}
	}
    if ($ok) {
    	my $tmp_dir = $project->getTmpDir();
	    my $session = new CGI::Session(undef, $cgi, {Directory=>$tmp_dir});
	    $session->param('hash_xls', compress(freeze $h));
	    print '@@@';
	    print $session->id;
    }
    exit(0);
}

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

sub launch_filters_dejavu {
	my $chr = shift;
	return unless ($dejavu);
	$chr->project->dejavu_ho(1) if ($dejavu_ho);
	$chr->dejaVuFilter($dejavu);
}

#sub launch_filters_atLeast {
#	my $chr = shift;
#	if ($atLeast and $atLeast >= 2) {
#		$chr->project->nb_at_least($atLeast);
#		if ($level_ind eq 'gene') {
#			if ($typeFilters eq 'familial') { $chr->atLeastFilter_genes_ind_fam(); }
#			else { $chr->atLeastFilter_genes_ind(); }
#		}
#		else {
#			if ($typeFilters eq 'familial') { $chr->atLeastFilter_var_ind_fam(); }
#			else { $chr->atLeastFilter_var_ind(); }
#		}
#	}
#	if ($atLeastFam and $atLeastFam >= 2) {
#		$chr->project->nb_at_least_fam($atLeastFam);
#		if ($level_fam eq 'gene') { $chr->atLeastFilter_genes_fam(); }
#		else { $chr->atLeastFilter_var_fam(); }
#	}
#	return;
#}

sub launch_filters_atLeast {
	my $chr = shift;
	if ($atLeast and $atLeast >= 2) {
		$chr->project->nb_at_least($atLeast);
		if ($typeFilters eq 'familial') {
			if ($level_fam eq 'gene') { $chr->atLeastFilter_genes_fam(); }
			else { $chr->atLeastFilter_var_fam(); }
		}
		else {
			if ($level_fam eq 'gene') { $chr->atLeastFilter_genes_ind(); }
			else { $chr->atLeastFilter_var_ind(); }
		}
	}
	return;
}

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



##### XLS METHODS #####



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
	my $hFormat = getHashFormat($workbook);
	my (@lPatNamesHeader, @lPatNames, $hPolyphenSift);
	foreach my $pat (@{$project->getPatients()}) {
		my $aff_or_not = 'healthy';
		$aff_or_not = 'affected' if ($pat->status() eq '2');
		push(@lPatNamesHeader, $pat->name().':'.$aff_or_not);
		push(@lPatNames, $pat->name());
	}
#	my $need_decompress;
	unless ($h) {
		foreach my $patient (@{$project->getPatients()}) {
			$h->{by_pat}->{$patient->name()} = $patient->xls_var_id();
			foreach my $chr_id (keys %{$h->{by_pat}->{$patient->name()}}) {
				foreach my $id (keys %{$h->{by_pat}->{$patient->name()}->{$chr_id}}) {
					$h->{by_var}->{$chr_id}->{$id}->{$patient->name()} = undef; 
				}
			}
		}
#		foreach my $chr_id (sort {$a <=> $b} keys %{$h->{by_var}}) {
##			my $hRes;
##			$hRes->{by_pat} = $h->{by_pat};
##			$hRes->{by_var}->{$chr_id} = $h->{by_var}->{$chr_id};
##			push(@$listHash, $hRes);
#			$h->{by_var}->{$chr_id} = undef;
##			$hRes = undef;
#		}
	}
	
	my $xls_var = $workbook->add_worksheet('RESULTS');
	my $xls_var_not_merge = $workbook->add_worksheet('RESULTS NOT MERGED');
	my $listColNames = writeHeader_byVar($workbook, $xls_var, \@lPatNamesHeader);
	my $listColNames_not_merge = writeHeader_byVar($workbook, $xls_var_not_merge, \@lPatNamesHeader);
	my $i = 1;
	
	my $col_cons = scalar(@lPatNames) + 16;
	my $col_polyphen = $col_cons + 12;
	my $col_sift = $col_cons + 13;
	
#	foreach my $h_tmp (@$listHash) {
#		my $h;
#		if ($need_decompress) { $h = thaw(decompress $h_tmp); }
#		else { $h = $h_tmp; }
#		$h_tmp = undef;
		my @lLines;
		foreach my $chr_id (sort {$a <=> $b} keys %{$h->{by_var}}) {
#			my $i = 1;
#			my $worksheet_name = 'CHR'.$chr_id;
#			$worksheet_name = 'CHRX'  if ($chr_id == 23);
#			$worksheet_name = 'CHRY'  if ($chr_id == 24);
#			$worksheet_name = 'CHRMT' if ($chr_id == 25);
#			my $xls_var = $workbook->add_worksheet($worksheet_name);
#			my $listColNames = writeHeader_byVar($workbook, $xls_var, \@lPatNames);
			#my $hLen = initHashSizeColumn_byVar($listColNames);
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
				
				#warn Dumper (sort @lPatNames);
				#warn Dumper $h->{by_var}->{$chr_id}->{$id};
				
				#die;
				
				foreach my $pat_name (sort @lPatNames) {
					if (exists $h->{by_var}->{$chr_id}->{$id}->{$pat_name}) {
						my $he_ho = $h->{by_pat}->{$pat_name}->{$chr_id}->{$id}->{he_ho};
						my $nb_ref = $h->{by_pat}->{$pat_name}->{$chr_id}->{$id}->{nb_all_ref};
						my $nb_mut = $h->{by_pat}->{$pat_name}->{$chr_id}->{$id}->{nb_all_mut};
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
				push(@lCol, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{clinvar});
				push(@lCol, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{min_pop_freq});
				push(@lCol, $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{max_pop_freq});
				
				my @lGenes;
				foreach my $gene_name (keys %{$h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}}) {
					unless (exists $h->{by_pat}->{$firstPat}->{$chr_id}->{$id}->{genes}->{$gene_name}->{transcripts}) {
						$lCol[2] = 'Intergenic';
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
				foreach my $pat_name (@lPatNames) {
					delete $h->{by_pat}->{$pat_name}->{$chr_id}->{$id} if (exists $h->{by_pat}->{$pat_name}->{$chr_id}->{$id});
				}
				delete $h->{by_var}->{$chr_id}->{$id};
				
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
							while ($z < (14 + scalar(@lPatNames))) {
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
							if ($j < (15 + scalar(@lPatNames))) {
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
#		if ($need_decompress) { $h = undef; }
#	}
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
	my $header = "Variation	Type	Consequence	Dejavu	Chr	Position	Allele	Sequence	".join("\t", @listPatNames)."	He	Ho	Cadd	Cosmic	ClinVar	Min_Pop_Freq	Max_Pop_Freq	Gene	Consequence	Transcript	Transcript_Xref	Description	Exon	Cdna_Pos	Cds_Pos	Protein	Protein_xref	AA	Protein_Pos	Nomenclature	Polyphen	Sift";
	my @lNames = split(' ', $header);
	my $hFormat = getHashFormat($xls);
	my $j = 8;
	foreach my $name (@listPatNames) {
		my $status = $hStatus->{$name};
		my $path_polyweb;
		if (-d $Bin.'/../../polyweb/') { $path_polyweb = $Bin.'/../../polyweb/'; }
		elsif (-d $Bin.'/../../PolyWeb/') { $path_polyweb = $Bin.'/../../PolyWeb/'; }
		if (-d $Bin.'/../../../polyweb/') { $path_polyweb = $Bin.'/../../../polyweb/'; }
		elsif (-d $Bin.'/../../../PolyWeb/') { $path_polyweb = $Bin.'/../../../PolyWeb/'; }
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


