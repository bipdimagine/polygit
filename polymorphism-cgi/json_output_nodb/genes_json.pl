#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/kyoto";
#use lib "/bip-d/soft/distrib/tabix/latest/perl";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/util_filter";
use filter_cache qw/get_combination_param_for_cache get_data delete_variations  new_intersection  /;  
use connect;
use GBuffer;
use GenBoStorable;
use Data::Dumper; 
use Bio::DB::Sam;
use get_variations;
#use Cache::File;
use family_filter;  
use patient_filter;
use somatic_filter;
use util_file;
use Storable qw/freeze thaw nfreeze nstore_fd nstore retrieve/;
use Digest::MD5::File qw( file_md5_hex); 
use layout;
use export_excel; 
use export_data;
#use Set::Object; 
use Set::Intersection;
use GenBoFilter;
use GenBoProjectQuery;
use Set::IntSpan::Fast::XS ;
use Tabix;
use Set::Intersection;
use JSON;

my $nb_genes_all;
my $debug_genes;
my %formatter = (
	"formatSnp"         => \&formatter::formatSnp,
	formatMethods       => \&formatter::formatMethods,
	formatReferenceGene => \&formatter::formatReferenceGene,
	formatFilter        => \&formatter::formatFilter,
	formatValid    => \&formatter::formatValid,
	formatPolyphen => \&formatter::formatPolyphen,
	formatHapmap  => \&formatter::formatHapmap
);

my @filter_param_non_statistic = (
   "polyphen3",
   "polyphen2",
   'polyphen1',
   "polyphen0",
   'sift2',
   'sift1',
   'sift0',
   'dbsnp_none',
   'dbsnp_1p',
   'dbsnp_clinical',
   'evs',
   'evs_1p',
   '1000genomes_1p',
   '1000genomes',
    'intronic',
    'score0',
    'score1',
    'score2',
    'bipd',
    'pseudogene',
    'intergenic',
    'non-frameshift',
    'ncrna',
    'maturemirna'
    );
our $apph;

my %hfilter_param_non_statistic;
@hfilter_param_non_statistic{@filter_param_non_statistic} = undef;

our @STATISTICS_PARAMS = ("homozygote","heterozygote");

my @special_filters = (   "essential_splicing",
						  "splicing",
						  "phase",	
						  "coding",
						  "silent",
						  "phase",
						  "non-frameshift",
						  "frameshift",
						  "stop");
my %hspecial_filters;
@hspecial_filters{@special_filters} = undef;	


				  

my %types = ( variations => "variations" );

my $cgi    = new CGI();
my $buffer = new GBuffer;
my $delete_genes; 
my $project_name = $cgi->param('project');
my $project = $buffer->newProject( -name => $project_name );
die( "unknown project" . $project_name ) unless $project;
my $ped_file = $project->getRootDir()."/../".$project->name().".ped";

my $pedigree ={};
my $ped_fam =[];

if ($cgi->param('get_patients')) {
	my @lPatNames;
	foreach my $pat (@{$project->getPatients()}) { push(@lPatNames, $pat->name()) }
	my $hash;
	$hash->{'patients'} = join(',', sort(@lPatNames));
	print $cgi->header('text/json-comment-filtered');
	print encode_json $hash;
	exit(0);
}
	


#
############
## FILTER 
#####
my $filter_name= $cgi->param('filter_name');
my $user_name= $cgi->param('user');
my $run_cache = $cgi->param('run_cache');
my $dejavu_ho = $cgi->param('dejavu_ho');
if ($filter_name) {
my $query = $buffer->queryFilter();
my $user_id =$query->get_user_id(user_name=>$user_name);	
die("no user") unless $user_id;
my $filter_id =$query->get_filter_id(filter_name=>$filter_name,project_id=>$project->id, user_id=> $user_id);	
my $params = $query->getParam(filter_id=>$filter_id);

foreach my $name (keys %$params) {
	#next if ($name =~ /region/);
	$params->{$name}->{value} =~ s/\+/ /g;
	$cgi->param(-name=>$name,-value=>$params->{$name}->{value});
	
}
}




# Load gene_atlas file
my $file = $buffer->{config}->{gene_atlas_diseases}->{dir}.$buffer->{config}->{gene_atlas_diseases}->{genes};
	my $genes_atlas;
	if (-e $file) {
 		$genes_atlas = retrieve($file);
	}





my $trio_recessif_model =1;



#creation du projet a partir du buffer






my $filter_he = $cgi->param('filter_he');
my $filter_ho = $cgi->param('filter_ho');

#http://10.200.27.102/cgi-bin/polymorphism-cgi//json_output_nodb/genes_json.pl?project=NGS2016_1027&stat=all&filter_type_variation=dbsnp_none+1000genomes+1000genomes_1p+evs+evs_1p+silent+intergenic+intronic+pseudogene+loh+dbl_evt+polyphen_sift&mode=ind&level_ind=variation&level_fam=variation
my @filter_type_variations  = grep {$_ ne "dbsnp_none"} split(" ",$cgi->param('filter_type_variation'));


my ($polyphen_sift) = grep {$_ eq "polyphen_sift"} @filter_type_variations;
my @list_polyphen ;
my @list_sift;
if ($polyphen_sift){
	my @new_list;
	
		foreach my $f (@filter_type_variations) {
			if ($f =~ /polyphen[0-3]/) {
				push(@list_polyphen,$f);
			}
			elsif ($f =~ /sift[0-2]/) {
				push(@list_sift,$f);
			}
			else {
				push(@new_list,$f);
			}
			}
		 @filter_type_variations = @new_list;
}


#args for variation filter
my %hfilter_type_variations;

my $data = load_data();
my $delete_bipd;
my $dejavu =  $cgi->param('dejavu');
if ($dejavu > 0){
	$delete_bipd = 1;
#	$hfilter_type_variations{bipd} = 1;
	push(@filter_type_variations,"bipd");
}



#filter pour les patients par nb 


my $filter_nb = $cgi->param('nb');

my $filter_diseases = $cgi->param('filter_diseases');


my $filter_patient = $cgi->param('filter_patient');
my $exclude_patient = $cgi->param('filter_not_patient');
my $attic_patients = $cgi->param('filter_attic')."";

my $patient_and = {};
 map{$patient_and->{$_} = undef} split(" ",$cgi->param('filter_patient'));

my $patient_not = {};
 map{$patient_not->{$_} = undef} split(" ",$cgi->param('filter_not_patient'));

my $patient_attic = {};
 map{$patient_attic->{$_} = undef} split(" ",$cgi->param('filter_attic'));
	
my $fam_and = {};
 map{$fam_and->{$_} = undef} split(" ",$cgi->param('fam_and'));
 
 my $fam_not = {};
 map{$fam_not->{$_} = undef} split(" ",$cgi->param('fam_not'));


my $stat = $cgi->param('stat');
my $filter_gene_level_ind = undef;
$filter_gene_level_ind = 1 if ($cgi->param('level_ind') eq "gene") ;

my $filter_gene_level_fam = undef;
$filter_gene_level_fam = 1 if ($cgi->param('level_fam') eq "gene") ;

my $filter_text = $cgi->param('filter_text') ;
my $model = $cgi->param('model') ;
my $filter_loh = $cgi->param('filter_loh');
my $filter_dblevt = $cgi->param('filter_dblevt');

my %selected_chromosome ;



my @filter_region =  split(" ",$cgi->param('filter_region')) ;
push(@filter_region, split(" ",$cgi->param('filter_regions'))) ;
my $gene_ontology;
my @search_text = split(" ",$filter_text);
if ($filter_text) {
	$gene_ontology = find_gene_ontology(\@search_text);
}

my $regions_set;
my %exclude_regions;
if (@filter_region) {
	foreach my $r (@filter_region){
	
		my ($chr,$start,$end,$include) = split(":",$r);
		
		next if $chr eq ""; 
		$include = 0 unless $include;
		
		my $set = Set::IntSpan::Fast::XS->new($start."-".$end);
		my $item;
		
		$item->{chromosome} = $chr;
		$item->{start} = $start;
		$item->{end} = $end;
		$item->{set} = $set;
		$item->{include} = $include;
		if ($include==0) {
		if (exists $selected_chromosome{$chr}) {
			
			$selected_chromosome{$chr} = $selected_chromosome{$chr}->union($set);
		}
		else {
			
			$selected_chromosome{$chr} = Set::IntSpan::Fast::XS->new($start."-".$end); 
		}
		}
		elsif ($include == -1) {
			
			if (exists $exclude_regions{$chr}){
				
				$exclude_regions{$chr} = $exclude_regions{$chr}->union($set);
			}
			else {
				
				$exclude_regions{$chr} = Set::IntSpan::Fast::XS->new($start."-".$end); 
			} 
			
		}
		push(@$regions_set,$item);
	
	}
}

my @filter_chromosome;
my $param_filter_chromosome = $cgi->param('filter_chromosome');
if ($param_filter_chromosome eq 'getAll') {
	@filter_chromosome = (1..22, 'X', 'Y', 'MT');
}
else { @filter_chromosome = split(" ", $param_filter_chromosome); }

if (@filter_chromosome) { 

	foreach my $chr1 (keys %selected_chromosome) {
		my $find = grep {$_ eq $chr1 }@filter_chromosome;
		unless ($find){
			delete $selected_chromosome{$chr1};
		}
		
	}
	foreach my $chr (@filter_chromosome){
		unless (exists $selected_chromosome{$chr}){
			my $set = Set::IntSpan::Fast::XS->new();
			$set->add_range(0,$set->POSITIVE_INFINITY);
			$selected_chromosome{$chr} = $set;
		}		
	}
		

}


%selected_chromosome = undef if $filter_chromosome[0] eq "all";

$filter_patient = undef if $filter_patient eq "all";
$exclude_patient = undef if $exclude_patient eq "all";

my $all_data;

my $patients =  $data->{patients};


push(@STATISTICS_PARAMS,@$patients);

$all_data = $data->{chromosomes};
my $toto;
###################
# store data for statistics
#######################""
my $stats_by_patients;
my $stats_by_chromosome;
my $stats_by_families;
#my $p = new patient_filter(data=>$all_data,patients=>$patients);
my $fm;
my $type =0;
$type =  $data->{type_cache} if exists $data->{type_cache};


my $mode = $cgi->param('mode');

if ($mode eq "fam"){
	
  $fm =  new family_filter(patients=>$patients,data=>$all_data,project=>$project, type=>$type );
 
}
else {
	$fm =  new patient_filter(patients=>$patients,data=>$all_data,project=>$project,type=>$type);
}
my $ps = $fm->patients();

my $pair_file = $project->getRootDir()."/../".$project->name().".ped";
#if (-e $pair_file){
#	parse_pair($pair_file);
#}
if (-e $ped_file){
	($pedigree,$ped_fam) = parse_ped($ped_file);
	$fm->pedigree($ped_fam);
	$fm->individual_pedigree($pedigree);
}
else {
	#/* modificvation ICI *******/
	my $pedigree;
		$pedigree->{fam} = "none";
		$pedigree->{parents} = [];
		$pedigree->{samples} = $patients;
		$pedigree->{child_d} = [];
		$pedigree->{child_m} = $patients;
		$fm->pedigree([$pedigree]);
}

if ($project->isSomaticStudy()) {
	$fm->somatic_file($project->getSomaticFile);
	#($pedigree) = parse_somatic($project->getSomaticFile);
	
}

$fm->set_attic($attic_patients);
my $compt =0;
foreach my $chr (keys %$all_data){
	
	  my $filter_start = -1 ;
	  my $filter_end = -1;
	  my $global_patients;
	  my $global_families;
	  my $global_chr;
	  ### 
	  # initialise les regions a filtrer 
	  ###
	  my $filter_set;
	
	if (scalar(keys %selected_chromosome) > 0){
	 if (exists $selected_chromosome{$chr}){
	 	$filter_set = 	 $selected_chromosome{$chr};
	 	}
		else {# && !(exists $selected_chromosome{$chr})) { 
			
			
			delete $all_data->{$chr};
			next;
		}
	}
	next unless exists $all_data->{$chr}->{genes};
	$fm->chr($chr);
	my $file =  $project->getCacheGenesKyotoFile($chr);
	
	 my $nb_bipd;
	 my $kyoto_cache_filter_genes = $project->getCacheFilterGenesKyotoFile($chr);
	 my $gene_no_data = {};
	if (-e $kyoto_cache_filter_genes ){
		my $fgenes;

	
		my $db1 = new KyotoCabinet::DB;
		if (!$db1->open($kyoto_cache_filter_genes, $db1->ONOLOCK | $db1->OREADER)){
					warn  " [$kyoto_cache_filter_genes] open error: %s\n". $db1->error;
					confess();
			}
			
	
		my $toto =  $db1->get($data->{prefiltering}) ;
		$all_data->{$chr}->{genes} = thaw $toto  if  $toto ;
		$db1->close;

	}
	
### boucle des genes;	
#	warn "1 ".scalar keys %{$all_data->{$chr}->{genes}};
#	warn "2 ".scalar keys %$gene_no_data;#."-".scalar keys %{$all_data->{$chr}->{genes}};
	foreach my $gid (keys  %{$all_data->{$chr}->{genes}}) { 
	$compt++;
	$fm->gid($gid);
	if (exists $hfilter_type_variations{intergenic} && $fm->gene_data()->{intergenic} ){
			delete $all_data->{$chr}->{genes}->{$gid};
			next;
	}	
	#####
	# create hash all variations if not exists
	#####
	my @filter =  keys %{$fm->gene_data()};

	my $debug;
	if (scalar (@filter) ==0){
		delete $all_data->{$chr}->{genes}->{$gid};
		next;
	}
	
	#$fm->gid($gid);
	next unless keys %{$fm ->gene_data};
	#$fm->set_gid();
	#die();
	if (exists $delete_genes->{$gid}){
		delete $all_data->{$chr}->{genes}->{$gid} ;
		next;
	}
		#####
		# filtre gene atlas diseases
		####
		if ($filter_diseases) {
				my $gname =  $fm->gene_data()->{name};
				unless (exists $genes_atlas->{$gname}->{ids}->{$filter_diseases}) {
					$fm->delete();
			 		next;
				}
					
		}
		
		####
		# filtre by regions
		###
		
		if (exists $selected_chromosome{$chr}){
			return_genes_annot($fm->gene);
			my $start = $fm->gene->{start};
			my $end = $fm->gene->{end};
			my $tempset =  Set::IntSpan::Fast::XS->new($start."-".$end);	
			my $in = $tempset->intersection($selected_chromosome{$chr});
			if ($in->is_empty()){
				$fm->delete();
			 	next;
			}		
		}
	
	 ######
	  ## here optimization delete all key not useful , 1 no statistic on this key and not part of a filter 
	  #####
	  	if ($polyphen_sift && @list_sift && @list_polyphen){
			$fm->prepare_sift_polyphen(\@list_polyphen,\@list_sift);
					 
		}
	  
#	  foreach my $filter (keys %hfilter_param_non_statistic){
#	  	delete $all_data->{$chr}->{genes}->{$gid}->{data}->{$filter};
#	  } 
#	
	
		########################
		# Patient in the attic  
		########################
		if($attic_patients){
			$fm->in_the_attic();
			if ($fm->isEmpty()){
				$fm->delete();
				next;
			}
		}
		my @delete_keys;	
		if (scalar(@filter_type_variations)){
			
			if ($delete_bipd){
				###########################
				# filter bipd on the fly
				###########################
				#$fm->delete_key("bipd");
			
				#warn "start";
				my $ids;
				if ($dejavu_ho){
				
					$ids = $project->getDejaVuHo($chr,$fm->gene_data->{homozygote},$dejavu);
					my %diff;
					if ($fm->gene_data->{heterozygote}){
						
						foreach my $t (@{$fm->gene_data->{heterozygote}}){
								$diff{$t} ++;						
							
						}
						if ($fm->gene_data->{homozygote}){
							foreach my $t (@{$fm->gene_data->{homozygote}}){
								delete $diff{$t}  if exists $diff{$t};						
							
						}
						}
					}
				
					push(@$ids,@{$project->getDejaVuHo($chr,[keys %diff],$dejavu-1)});
				}
				else {
					$ids = $project->getDejaVu($chr,[keys  %{$fm->hashvar}],$dejavu);
				}
			
				$fm->put("bipd",$ids);
				#$all_data->{$chr}->{genes}->{$gid}->{data}->{bipd} = \@ids;
				

				if (scalar(@$ids) == scalar( keys %{$fm->hashvar}) ) {
				
					$fm->delete();
					next;
				}
			
			}
			
		
		
			my @filters_step1 = grep {! exists $hspecial_filters{$_}} @filter_type_variations;
			
			my @filters_step2 = grep { exists $hspecial_filters{$_}} @filter_type_variations;
			#next;
		#	push(@delete_keys,@filters_step1);
#			my $result = delete_variations(\@filters_step1,$all_data->{$chr}->{genes}->{$gid},1);
#			next;
			 $fm->delete_variations(\@filters_step1);
			if ($fm->isEmpty()) {
				
				$fm->delete();
				next;
			}
			
			$fm->delete_special_filters(\@filters_step2);
		
			if ($fm->isEmpty()) { 
			#	die();
				$fm->delete();
				next;
			}
 
			$fm->refresh_data();
			
		}
		#die();
	# prepare heterozygote et homozygote by patient 
		#next;
	
	#my $gene = $all_data->{$chr}->{genes}->{$gid};
 
	
	
	########################
	# exclude heterozygote for one patient  
	########################
	
	if ($filter_he)  {
		#
		## for heterozygote 
		# 	
		
		my @t_patients =  split(" ",$filter_he);
		$fm->exclude_he(\@t_patients);
		
	}
	########################
	# exclude homozygote for one patient  
	########################
	
	if ($filter_ho) {
		my @t_patients =  split(" ",$filter_ho);
		$fm->exclude_ho(\@t_patients);
		push(@delete_keys,map {"homozygote_".$_} @t_patients);
		
	}
	
	
	################################
	# refresh data for he ho   
	################################
	if ($filter_ho || $filter_he) {
			if ($fm->isEmpty()) {
				$fm->delete();
				next;
			}
	
		
		$fm->refresh_data();
	}
	
 	########################
	# exclude patient simple  
	########################
	
	if ($exclude_patient) {		
		if ($filter_gene_level_ind) {
			$fm->exclude_genes([keys %$patient_not]);
		}
		else {
		
			$fm->exclude_variations([keys %$patient_not]);
			
		}
		
		if ($fm->isEmpty()) {
		
			$fm->delete();
			next;
		}
	
		$fm->refresh_data();
	}

	
		##################"
	# NOT  famille    
	#############
	
	if (keys %$fam_not){
		if ($filter_gene_level_fam) {
		
			$fm->not_inter_fam_genes($fam_not);
		}
		else {
			$fm->not_inter_fam_variations($fam_not);
		}
		
		if ($fm->isEmpty()) {
			$fm->delete();
			next;
		}
		$fm->refresh_data();
	}

	##################"
	# AND  famille    
	#############
	
	if (keys %$fam_and){
		if ($filter_gene_level_fam) {
		
			$fm->and_inter_fam_genes($fam_and);
		}
		else {
			$fm->and_inter_fam_variations($fam_and);
		}
		
		if ($fm->isEmpty()) {
			$fm->delete();
			next;
		}
		$fm->refresh_data();
	}
	  ##################"
	# AND  patient simple  
	#############
	if ($filter_patient) {
	
		
	
		if ($filter_gene_level_ind) {
		
			$fm->and_genes([keys %$patient_and]);
		}
		else {
			$fm->and_variations([keys %$patient_and]);
		}
		
		if ($fm->isEmpty()) {
			$fm->delete();
			next;
		}
		$fm->refresh_data();
		
	}
	###################
	# recessif model heterozygote composite or at least one homozygote;  
	###################
	if ($model){
		my $keep_var;
		if ($model eq "recessif"){
			 $keep_var = $fm->recessif();
		}
		elsif  ($model eq "denovo"){
			 $keep_var = $fm->denovo();
		}
		elsif ($model eq "strict-denovo"){
			 $keep_var = $fm->denovo();
				$fm->moredenovo();
		}
		elsif ($model eq "dominant"){
			$keep_var = $fm->dominant();
		}
		elsif ($model eq "compound"){
				$fm->filter_composite();
			#$fm->filter_composite($nb,$filter_identical_gene);
		}
		
		if ($fm->isEmpty()) {
			
			$fm->delete();
			next;
		}
		
		$fm->refresh_data();
	}

	##################"
	# AND  famille    
	#############
	
	if (keys %$fam_and){
		if ($filter_gene_level_fam) {
		
			$fm->and_inter_fam_genes($fam_and);
		}
		else {
			$fm->and_inter_fam_variations($fam_and);
		}
		
		if ($fm->isEmpty()) {
			$fm->delete();
			next;
		}
		$fm->refresh_data();
	}
	
			
 
	
	
  ##################"
	# AND  patient simple  
	#############
	if ($filter_patient) {
		
		
	
		if ($filter_gene_level_ind) {
		
			$fm->and_genes([keys %$patient_and]);
		}
		else {
			$fm->and_variations([keys %$patient_and]);
		}
		
		if ($fm->isEmpty()) {
			$fm->delete();
			next;
		}
		$fm->refresh_data();
		
	}
	###################
	# SOMATIQUE FILTER 
	###################
	if ($filter_loh) {
	$fm->loh();
		if ($fm->isEmpty()) {
				$fm->delete();
				die() if exists  $fm->{data}->{$chr}->{genes}->{$gid};
				next;
		}
		$fm->refresh_data();
	}
	
	if ($filter_dblevt) {
	$fm->dbl_evt();
		if ($fm->isEmpty()) {
				$fm->delete();
				die() if exists  $fm->{data}->{$chr}->{genes}->{$gid};
				next;
		}
		$fm->refresh_data();
	}
	

	 

	
	if ($filter_nb > 0){
		my $filter_gene_level = $filter_gene_level_ind;
		
		if ($mode eq "fam"){
			 $filter_gene_level = $filter_gene_level_fam;
		}
		if ($filter_gene_level) {
			$fm->identical_genes($filter_nb);
		}
		else {
			$fm->identical_variations($filter_nb);
		}
		if ($fm->isEmpty()) {
				$fm->delete();
				next;
		}					
		$fm->refresh_data();
	}
	
if ($filter_text) {
			my $gname = $fm->gene->{name};
			my $gid = $fm->gene->{id};
			return_genes_annot($fm->gene);
			
			
			my $find;
			foreach my $search (@search_text) {
			my $text = lc($genes_atlas->{$gname}->{text}).lc($fm->gene->{description})." ".lc($fm->gene->{name})." ".lc($fm->gene->{xref});
			my $r = index($text, lc($search) );
			
			if ($r ne -1) {
				$find =1;
				last;
			};
			}
			if ($gene_ontology){
				
					my $xref = $fm->gene->{xref};
				
					if (exists $gene_ontology->{$xref}) {
						$find =1 ;
						$fm->gene->{description} .= ";".join(";",@{$gene_ontology->{$xref}});
					}
			}
			
			
			unless ($find) {
				$fm->delete();
				next;
			}
			
			
		}

	
		
	###########################################################
	# exclude regions a la fin pour pouvoir faire les stats
	###########################################################
	
	
	if (exists $exclude_regions{$chr}){
			my $start = $$fm->gene_data()->{start};
			my $end = $fm->gene_data()->{end};
			my $tempset =  Set::IntSpan::Fast::XS->new($start."-".$end);
			my $in = $tempset->intersection($exclude_regions{$chr});
			
			if (!($in->is_empty())){
				my $g = $fm->gene_data();
				foreach my $region (@$regions_set){	
					next if $chr ne $region->{chromosome};
				
					
					my $in2 = $tempset->intersection($region->{set});
					next if $in2->is_empty;
					$region->{genes}++;
					
					map{$region->{hvar}->{$_}++} keys %{$fm->hashvar };
					map{$region->{hcoding}->{$_}++} @{$fm->get("coding")};
					map{$region->{hcoding}->{$_}++} @{$fm->get("stop")}  ;
					map {$region->{insertion}->{$_}++}	@{$fm->get("insertion")};
					map {$region->{subs}->{$_}++}	@{$fm->get("substitution")};
					map {$region->{deletion}->{$_}++} @{$fm->get("deletion")};;
			}
				$fm->delete();
			 	next;
			}
		}
			
		my $gene = $fm->gene_data();
		###########################
		###### PREPARE STATISTIC 
		###########################
		
		my %dejavu_fam;
		my $ped = $fm->individual_pedigree();
		
		
		my @stats = ("substitution","insertion","deletion","homozygote","large_deletion");
	
	 	   	foreach my $stat (@stats){
	    			foreach my $v (@{$fm->get($stat)}){
	    				my $vstat = $stat ;
	    				$vstat = "deletion" if $stat eq "large_deletion";
	    				$global_patients->{$vstat}->{$v}=undef ;
	    				$global_families->{$vstat}->{$v}=undef ;
	    			}
	    		}
		
		 my @stats = ("coding","stop","frameshift","phase");
	     foreach my $stat (@stats){
	     	map {$global_patients->{coding}->{$_}++} @{$fm->get($stat)};
	     }
	     

		foreach my $patient (@$patients) {
			next unless exists $fm->{gene_data}->{$patient};
			next unless scalar(@{$fm->{gene_data}->{$patient}});
		
			my $toto = $fm->get($patient);
			my $fam = $ped->{$patient}->{fam};
		
		
			my $nb =  scalar(@$toto);

			next if $nb == 0;#exists $gene->{data}->{$patient};
			$stats_by_families->{$fam}->{genes} ++ unless exists $dejavu_fam{$fam} ;
			$dejavu_fam{$fam} ++ ;
			$stats_by_patients->{$patient}->{gene} ++;
			$stats_by_patients->{$patient}->{composite} ++ if $nb> 1;
			foreach my $v (@$toto){
				$global_patients->{$patient}->{$v}=undef;
				$global_families->{$fam}->{$v} = undef;
			
			}
			
			
	
		 } #end statistic for patients
		 
		 
		
	
	
	   $nb_genes_all ++;
		$stats_by_chromosome->{$chr}->{gene} ++; 
		my @status = ("utr","splicing","silent");
		foreach my $stat (@status){
			$fm->get($stat);
		}
		
	}# end genes
my @t = keys %{$global_patients->{452}};
	my @stats = ("substitution","insertion","deletion","homozygote","coding");

	foreach my $stat (@stats){
		my $nb = scalar(keys %{$global_patients->{$stat}});
		
		#$stats_by_chromosome->{$chr}->{variation} += $nb unless $stat eq "homozygote";
		
		$stats_by_chromosome->{$chr}->{$stat} = $nb;
	} 
	
	$stats_by_chromosome->{$chr}->{variation} = $stats_by_chromosome->{$chr}->{insertion} +$stats_by_chromosome->{$chr}->{deletion} + $stats_by_chromosome->{$chr}->{substitution} ;
	if (scalar( keys %{$all_data->{$chr}->{genes}}) == 0 ){
		delete $all_data->{$chr};
	}
my $ped = $fm->individual_pedigree();
my $all_vars_fam;

#
#compute stats by patients and by family
#



foreach my $fam (@{$fm->pedigree()}){
			my $h= {};
		my @stats = ("substitution","insertion","deletion");
		
			foreach my $patient (@{$fam->{samples}} ){
					my @toto = keys %{$global_patients->{$patient}};
					foreach my $t (@toto){
						foreach my $s (@stats){
								$stats_by_patients->{$patient}->{$s} ++ if exists $global_patients->{$s}->{$t};
								$h->{$s}->{$t} ++ if exists $global_patients->{$s}->{$t};
						}
							$stats_by_patients->{$patient}->{homozygote} ++ if exists $fm->heho_patients->{$patient}->{homozygote}->{$t};
				}
			}
		
		my @stats = ("substitution","insertion","deletion","homozygote");
		my $fname = $fam->{fam};
			foreach my $s (@stats){
				$stats_by_families->{$fname}->{$s}  += scalar(keys %{$h->{$s}});
			}
}
}#end chromosome
#warn $compt;

if ($stat eq "all"){
	my %hdata;
	$hdata{chromosomes} = by_chromosomes2($all_data,$patients);

	$hdata{genes} =  by_genes($all_data,$patients);
	$hdata{patients} =  by_patients($all_data,$patients,$stats_by_patients);
	$hdata{families} =  by_families($all_data,$patients,$stats_by_families);
	
	$hdata{diseases} = by_diseases($all_data,$patients);
	export_data::print($project,$cgi,\%hdata);
	exit(0);
}

if ($stat eq "variations"){
	my $var = all_variations($all_data);
	my @var_datas;
	my $ert =0;
	foreach my $chr_name (keys %$var){
		my $select_ids = [keys %{$var->{$chr_name}}];
		$ert += scalar(@$select_ids);;
		my $chr = $project->getChromosome($chr_name);
		die("ici") unless $chr;
		my $temp_type_name;
		my $tdata = get_variations::getIds($buffer,$project,$chr,$temp_type_name,$select_ids);
		my $user = my $filter_name= $cgi->param('user');
		export_data::update_deja_vu($project,$tdata,$user);
		#update_validations( $project, $data);
		push(@var_datas,@$tdata);
	}	
	export_data::variation_report_xls($project,\@var_datas);
	exit(0);
}

if ($stat eq 'vcf') {
	my $patientName = $cgi->param('only_patient');
	my $hVarIds = byVarIdsForVcf($all_data, $patientName);
	my $lVcfLines = getVcflines($hVarIds, $patientName);
	export_data::printVcfFile($project_name, $patientName, $lVcfLines);
}

my $out;

if ($stat){
	
	if ($stat eq "patient"){
	
		$out = by_patients($all_data,$patients,,$stats_by_patients);
	}
	elsif ($stat eq "region"){
		by_regions($all_data,$regions_set);
		exit(0);
	}
	elsif ($stat eq "button"){
		#by_buttons(\@filter_type_variations,$filter_composite);
		exit(0);
	}
	elsif ($stat eq "diseases"){
		$out = by_diseases($all_data,$patients);
		
	}

	
	else {
		$out = by_chromosomes2($all_data);
	}
	
}
else {
	
	$out = by_genes($all_data,$patients);
}


if ($run_cache == 1){
	
	exit(0);
}

	warn "end";
export_data::print($project,$cgi,$out);

exit(0);



my %hash_patients;

sub nb_patients {
	my ($variations,$data,$patients,$list) = @_;
#	my %hash;
#	@hash{@$variations} = ();
	my $nb =0;
	my %res;
	foreach my $name (@$patients) {
		next unless exists $data->{data}->{$name};
		my $t = new_intersection($variations,$data->{data}->{$name});
	
		#my @t = get_intersection($variations,$data->{data}->{$name});
		$res{$name} ++ if @$t;
	} 
	return scalar (keys %res) unless $list;
	return keys %res;
}

sub return_go {
	my ($xref) = @_;
	my $apph = $buffer->get_go_db();
	return "" unless $apph;
	my $terms = $apph->get_terms({product=>$xref});
	my $text;
	foreach my $term (@$terms){
		$text .=$term->name.";"
	}
	return $text;
}
sub return_genes_annot{
	my ($gene) = @_;
	return if exists $gene->{description};
	my $zz = $project->liteAnnotations->get("annotations",$gene->{id});;
	#my $zz = thaw $project->kyotoGenes->{$gene->{id}};
	#return unless $zz;
	unless ($zz){
		$gene->{start} = $gene->{end}-1;
		$gene->{end} = $gene->{end}+1;
		return;
	}
	$zz = {} unless ($zz);
	$gene->{xref} = $zz->{external_name};
	$gene->{description} = $zz->{description};
	$gene->{start} = $zz->{start};
	$gene->{end} = $zz->{end};
	
}

sub by_genes {
	my ($data,$patient) = @_;
	
	my $out;
	my $stop = 30;

	unless ($cgi->param("filter_chromosome") || $cgi->param('xls') ) {
		my $nb;
		foreach  my $c (keys (%$data)) {
			$nb += scalar(keys %{$data->{$c}->{genes}});
			last if $nb >200;
		}
		$stop = 1 if $nb> 200;  
	}
	
	foreach my $c (sort {$a cmp $b} keys (%$data)) {
		last if $stop == 0;
		$stop --;
		my $genes = $data->{$c}->{genes};
		my $variations = $data->{$c}->{variations};
		foreach my $gene (sort {$a->{start} <=> $b->{start}} values %$genes){
			my %item;
			$item{id} = $gene->{id};
			return_genes_annot($gene);
			#my $zz = thaw $project->kyotoGenes->{ $gene->{id}};
		#	warn Dumper($zz);
			#die();
			$item{name} = $gene->{name};
			my $name = $item{name};
			$item{xref} = $gene->{xref};
			$item{chromosome} = $c;
			$item{coding_cover} = $gene->{cover};
			$item{start} = $gene->{start};
			$item{end} = $gene->{end};
			$item{description} = $gene->{description};
			$item{gene_atlas_diseases} = "";
			if (exists $genes_atlas->{$item{name}}) {
				$item{gene_atlas_diseases} = 1;
						
			}
#			if ($filter_text){
#				my $text = return_go($gene->{xref});
#				$item{description} .= $text;
#			}
			$item{v_all} = scalar(keys %{$gene->{hash_all_variations}});
			$item{p_all} = nb_patients($gene->{hash_all_variations},$gene,$patient);
			foreach my $p (@$patient){
				next unless $gene->{data}->{$p};			
				$item{$p} = scalar(@{$gene->{data}->{$p}});
				}
			my @status = ("substitution","insertion","deletion","coding","utr","splicing","stop","silent","phase","frameshift");
			foreach my $st (@status){
				my $kv = "v_".$st;
				my $kp = "p_".$st;
				$item{$kv} = 0;
				$item{$kp} = 0;
				if (exists $gene->{data}->{$st}){
					$item{$kv} = scalar(@{$gene->{data}->{$st}});
					my %hash;
					@hash{@{$gene->{data}->{$st}}} = ();
					$item{$kp} = nb_patients(\%hash,$gene,$patients);
				}
				else {
					$gene->{data}->{$st} =[];
				}
			}
			$item{"coding_consequence"} = $item{v_coding}+$item{v_phase}+$item{v_stop}+$item{v_frameshift};
				my %hash;
				@hash{(@{$gene->{data}->{coding}},@{$gene->{data}->{frameshift}},@{$gene->{data}->{stop}},@{$gene->{data}->{phase}})} = ();
			$item{p_coding_consequence} = nb_patients(\%hash,$gene,$patients);
			
			$item{v_focus} = $item{v_phase}+$item{v_stop};
			my %hash2;
			@hash{@{$gene->{data}->{stop}},@{$gene->{data}->{phase}}} = ();
			$item{p_focus} = nb_patients(\%hash2,$gene,$patients);
			$item{ids} = join(",",keys %{$gene->{hash_all_variations}});
			push(@$out,\%item);
		}
	}
	
	unless ($out){
		my %item;
		$item{id} =500;
		$item{name} = "Nada";
			$item{xref} = "Nil";
			$item{chromosome} = "Zip";;
		
			$item{description} = "!!! --- Desperate empty --- !!!";
		
			$item{coding_cover} = -99;
			$item{start} = "zero";
			$item{end} = "zilch";
		
			$item{v_all} = -1;
			
			$item{p_all} = -1;
			my @status = ("substitution","insertion","deletion","coding","utr","splicing","stop","silent","phase");
			foreach my $st (@status){
				my $kv = "v_".$st;
				my $kp = "p_".$st;
				$item{$kv} = 0;
				$item{$kp} = 0;
			}
			
			
		push(@$out,\%item);

	}
	return $out;
	export_data::print($project,$cgi,$out);
	exit(0);
}


sub by_diseases {
	my ($data,$patient) = @_;
	
	
	
	my $out;
	my %cat;
	my %sub_cat;
	
	my $genes_cat;
	my $patients_cat;
	foreach my $c ( keys (%$data)) {
	
		my $genes = $data->{$c}->{genes};
	
		my $variations = $data->{$c}->{variations};
		foreach my $gene (values %$genes){
		my $name =  $gene->{name};
		
		if (exists $genes_atlas->{$name}) {
			
			my @pp = nb_patients($gene->{hash_all_variations},$gene,$patient,1);
			my %hgenes;
				foreach my $id (keys %{$genes_atlas->{$name}->{ids}}){
				
					foreach my $p (@pp){
						$patients_cat->{$id}->{$p} ++;
					}
					$genes_cat->{$id}++;
				}				
			}
		}
	}

##################
# construct tree
##################
return construct_tree($genes_cat,$patients_cat);
exit(0);	

}
sub construct_tree {
	my ($gene_cat,$tree_cat) = @_;
	my $file2 = $buffer->{config}->{gene_atlas_diseases}->{dir}. $buffer->{config}->{gene_atlas_diseases}->{tree};
	my $tree = retrieve($file2);
	my %root;
	$root{name} = "GeneAtlas Diseases";
	$root{type} = "root";
	$root{id} = "root";
	foreach my $t (sort {$a->{name} cmp $b->{name}} values %$tree){
		
		my $item = add_children($t,$gene_cat,$tree_cat);
		next unless $item;
		
		$item->{type} = "category";
		push(@{$root{children}},$item);
	}

return 	[\%root];
export_data::print(undef,$cgi,[\%root]);
	
	exit(0);
}

my $uuid =0;
sub add_children {
	my ($cat,$genes_cat,$patients_cat) = @_;
	my $id = $cat->{id};
	return  unless exists $genes_cat->{$id}; 
	my %item;
	change_category_name($cat,$genes_cat,$patients_cat);
	$item{name} = $cat->{name};	
	$item{sid} = $cat->{id};
	$item{id} = $uuid++;
	$item{type} = $cat->{type};
	return \%item unless exists $cat->{children};
	foreach my $tt (values %{$cat->{children}}){
		my $its = add_children($tt,$genes_cat,$patients_cat);
		push(@{$item{children}},$its) if $its;
	}
	
	return \%item;
}	

sub change_category_name {
	my ($cat,$genes_cat,$patients_cat) = @_;
	my $id = $cat->{id};
	$cat->{name} .= "(g:".$genes_cat->{$id}." p:".scalar(keys %{$patients_cat->{$id}}).")"	;
}

sub intersection_for_stats {
	my ($hash1,$hash2,$hashres) = @_;
	foreach my $k (keys %$hash2){
		next unless exists $hash1->{$k};
		$hashres->{$k} =undef;
		delete $hash2->{$k};		
	}
}

sub by_families {
		my ($data,$patients,$stats) = @_;
		my @out;
		foreach my $fam (@{$fm->pedigree()}){
			my $n = $fam->{fam}; 
			my %item;
			$item{id} = $n;
			$item{name} = $n;
			my @astats = ("substitution","insertion","deletion","homozygote");
			foreach my $s (@astats){
				$item{$s} += $stats->{$n}->{$s};
			}
			#$item{variations} += $stats->{$n}->{variations};
			$item{genes} += $stats->{$n}->{genes};
			$item{include} = 1;
				$item{nb} = scalar(@{$fam->{samples}});
			$item{include} = 1;
			$item{include} = 0	if (exists $fam_and->{$n});
			$item{include} = -1	if (exists $fam_not->{$n});
			my $not_all =1;
			foreach my $s (@{$fam->{samples}} ){
				unless (exists $patient_attic->{$s} ) {
					$not_all =undef;
					last;
				}
				
			}
			$item{include} =2 if $not_all;
			$item{model} = $model;
				$item{nb} = scalar(@{$fam->{samples}});
				push(@out,\%item);
		}	
		return \@out;
}

sub by_patients {
	my ($data,$patients,$stats) = @_;
	
	my $coverage_dir = $project->getRootDir()."/align/coverage/";
	
	
	
	
	
	my @out;
	my @and_status = split(" ",$filter_patient);

	my @not_status = split(" ",$exclude_patient);
	my @attic = split(" ",$attic_patients);
	my @list_patient = join(" ",@$patients);
	my $fams;
	 
	 
	 
	foreach my $patient (@$patients){
		my $coverage_file = $coverage_dir.$patient.".cov.gz";
	
		my %item;
		$item{coverage} = -1;
		
		get_coverage($coverage_file,\%item);
		$item{id} = $patient;
		$item{name} = $patient;
		$item{variations} = 0;
		$item{genes} =0;
		$item{homozygote} = 0;
		$item{heterozygote} = 0;
		
		$item{filter_heho} = -1;
	
		$item{fam} = $pedigree->{$patient}->{fam};
		$item{status} = $pedigree->{$patient}->{status};
		$item{child} = 0;
		if ($fm->patients_groups){
			$item{group} = $fm->patients_groups->{$patient}->{group};
			$item{tissue} = "C";
			$item{tissue} = "T" if $fm->patients_groups->{$patient}->{status} ==2;
		}
		
		if ($pedigree->{$patient}->{father} eq "0" && $pedigree->{$patient}->{mother} eq "0") {
			if ($pedigree->{$patient}->{sex} == 1){
				$item{child} = "F";
			}
			else {
				$item{child} = "M";
			}
		
		}
		else {
		
			$item{child} = "C2";
			$item{child} = "C1" 	if ($pedigree->{$patient}->{sex} == 1);
			
		}
		$item{sex} = $pedigree->{$patient}->{sex};
		$item{composite} += $stats->{$patient}->{composite};
		$item{deletions} += $stats->{$patient}->{deletion} +  $stats->{$patient}->{large_deletion};
		$item{genes} +=$stats->{$patient}->{gene};
		
		$item{insertions}+= $stats->{$patient}->{insertion};
		$item{substitutions} += $stats->{$patient}->{substitution};
		$item{homozygote} += $stats->{$patient}->{homozygote};
		
		$item{variations} += $item{deletions}+$item{insertions}+$item{substitutions} ;
		$item{heterozygote} += $item{variations} - $item{homozygote};
		
		$item{include} = 1;
		if (exists $patient_and->{$patient}) {
		$item{include}=0;
		}
		elsif (exists $patient_not->{$patient}) {
			$item{include}=-1 ;
		}
		elsif (exists $patient_attic->{$patient}) {
			$item{include}=2 ;
		}
		if ($filter_he =~ /$patient/ && $filter_he) {
		
			$item{filter_heho} = 0;
		}
		
		elsif ($filter_ho =~ /$patient/  ) {
			
			$item{filter_heho} = 1;
		}
	
		push(@out,\%item);
		
	
		
	}
	@out = sort {$a->{fam} cmp $b->{fam} || $a->{child} cmp $b->{child}} @out;
	return \@out;
	#export_data::print($project,$cgi,\@out);
	exit(0);
}

sub get_coverage {
	my ( $coverage_file,$item) = @_;
	return unless -e $coverage_file;
	eval {
		my $tabix = new Tabix(-data =>$coverage_file);
		 my $res = $tabix->query("mean_all");
		 my @data;
		 while(my $line = $tabix->read($res)){
    		
				my($a,$b,$c) = split(" ",$line);
				 if ($b == 99){
				 	$b= "coverage";
				 
				 }
				 else {$b.="x";$c *=100;};
				 $item->{"$b"} = int($c);
				 
				 
			}
	 }

}
sub by_regions {
	my ($data,$region_set) = @_;
	my $data2; 
	foreach my $region (@$region_set){
		my $chr = $region->{chromosome};
		my $set = $region->{set};
		my %item;
		$item{id} = $chr.$region->{start}.$region->{end};
		$item{name} = $chr.$region->{start}.$region->{end};
		$item{chromosome} = $chr;
		$item{start} = $region->{start};
		$item{end} = $region->{end};
		$item{include} = $region->{include};
		if ($region->{include} == -1){
		$item{variations} = -1*scalar(keys %{$region->{hvar}});
		$item{coding} = -1*scalar(keys %{$region->{hcoding}});
		$item{stop} =  -1*scalar(keys %{$region->{hstop}});
		$item{insertion} = -1*scalar(keys %{$region->{insertion}});
		$item{deletion} = -1*scalar(keys %{$region->{deletion}});
		$item{substitution} = -1*scalar(keys %{$region->{subs}});
		$item{genes} = -1*$region->{genes};
	
		}
		else {
		my $nb_genes = 0; 
		my %hvar;
		my %hcoding;
		my %hstop;
		my %insert_dejavu;
		my %subs_dejavu; 
		my %deletion_dejavu;
		foreach my $g (values %{$data->{$chr}->{genes}}) {
			my $tempset =  Set::IntSpan::Fast::XS->new($g->{start}."-".$g->{end});	
			my $in = $tempset->intersection($set);
			next if ($in->is_empty);
			$nb_genes++;
			map{$hvar{$_}++} keys %{$g->{hash_all_variations}};
			map{$hcoding{$_}++} @{$g->{data}->{coding}}  if  exists  $g->{data}->{coding};
			map{$hcoding{$_}++} @{$g->{data}->{stop}}  if  exists  $g->{data}->{stop};
			map{$hstop{$_}++} @{$g->{data}->{stop}}  if  exists  $g->{data}->{stop};
			if ($g->{data}->{insertion}){ 
					map {$insert_dejavu{$_}++}	@{$g->{data}->{insertion}};
				}
			if ($g->{data}->{substitution}){
					map {$subs_dejavu{$_}++}	@{$g->{data}->{substitution}};	
			}
			if ($g->{data}->{deletion}){
					map {$deletion_dejavu{$_}++} @{$g->{data}->{deletion}};# if $gene->{data}->{deletion}
				}
			
			#$item{coding} += scalar(get_intersection($g->{all_variations},$g->{coding})) if exists $g->{coding};
		}
		$item{variations} = scalar(keys %hvar);
		$item{coding} = scalar(keys %hcoding);
		$item{stop} =  scalar(keys %hstop);
		$item{insertion} =  scalar(keys %insert_dejavu);
		$item{deletion} =  scalar(keys %deletion_dejavu);
		$item{substitution} =  scalar(keys %subs_dejavu);
		$item{genes} = $nb_genes;
		}
		push(@$data2,\%item);
	}
	#my @data3 = sort {$a->{id} <=> $b->{id}} @$data2; 
	
	export_data::print($project,$cgi,$data2);
	exit(0);
	
	
}



sub all_variations{
	my ($data) = @_;
	#confess();
		my $hvar;
	foreach my $c (sort {$a <=> $b} keys (%$data)) {
	
	
		foreach my $g (values %{$data->{$c}->{genes}}){
			map{$hvar->{$c}->{$_}++} keys %{$g->{hash_all_variations}};
		}
	}
return $hvar;
}

sub by_chromosomes2 {
	my ($data,$patients) = @_;
	my %chr;
	my $data2; 
	if (scalar (keys %$data) == 0) {
		my %item;
		$item{id} = "1";
		$item{name} = "Game Over";
		$item{genes} ="...";
		$item{value} = "none";
		push(@$data2,\%item);
	}
		my @coverage_value = (0,5,10,20,50);
	
	foreach my $c (sort {$a <=> $b} keys (%$data)) {
	
		my %item;
		
		$item{id} = $c;
		$item{id} = 23 if ($c eq 'X');
		$item{id} = 24 if ($c eq 'Y');
		$item{id} = 25 if ($c eq 'MT');
		$item{name} = $c;
		$item{genes} = $stats_by_chromosome->{$c}->{gene};
		$item{variations} = 0;
		$item{coding} = 0;
		$item{stop} =0;
	
		$item{variations} = $stats_by_chromosome->{$c}->{variation};
		$item{insertions} = $stats_by_chromosome->{$c}->{insertion};
 		$item{deletions} = $stats_by_chromosome->{$c}->{deletion};
 		$item{substitutions} = $stats_by_chromosome->{$c}->{substitution};
 		$item{homozygote} =  $stats_by_chromosome->{$c}->{homozygote};
 		$item{coding} =   $stats_by_chromosome->{$c}->{coding};
 		$item{heterozygote} =$item{variations} - $item{homozygote};

		
		push(@$data2,\%item);
	}

	my @data3 = sort {$a->{id} <=> $b->{id}} @$data2; 
	return \@data3;
	export_data::print($project,$cgi,\@data3);
	exit(0);
}


sub by_buttons {
	my ($var,$comp) = @_;
	my @data;
	foreach my $param (@$var){
		my %item;
		$item{id} = $param;
		$item{name} = $param;
		
		push(@data,\%item);
	}
	unless ($comp) {
		my %item;
	
		$item{id} = "filter_composite";
		$item{name} = "filter_composite";
		push(@data,\%item);
	}
	
	export_data::print($project,$cgi,\@data);
	exit(0);
}

sub cover_chromosome {
	my($chr,$project,$all_coding_span,$item) = @_;
	my %nbs;
	my $len =0;
foreach my $patient(@{$project->getPatients}) {
	my $cover_file = util_file::get_cover_file({project=>$project,patient_name=>$patient->name,coverage=>"all"});
	next unless -e $cover_file;
	my $span_cover = retrieve($cover_file);
	 $len += scalar($all_coding_span->{$chr}->as_array());
	 foreach my $cv (keys %{$span_cover->{$chr}}) {
	 	$nbs{$cv} += scalar($all_coding_span->{$chr}->intersection($span_cover->{$chr}->{$cv})->as_array);
	 }
}


foreach my $cv (keys %nbs) {
	
	$item->{"cover_".$cv} = int(($nbs{$cv}/$len)*10000)/100;
}

return 1;
}


sub get_genes_thread {
	my ($project,$chr_id) = @_;
	my $buffer1 = new GBuffer;
	my $id = GenBoStorable::getStoreId( $buffer1->dbh, $project->id, $chr_id,"test" );
	return unless $id;
	
	return (GenBoStorable::getStore( $buffer1->dbh, $id ));
}


sub get_genes {
	my ($project,$chr_id) = @_;
	
	my $id = GenBoStorable::getStoreId( $project->buffer->dbh, $project->id, $chr_id,"test" );
	return unless $id;
	
	return (GenBoStorable::getStore( $project->buffer->dbh, $id ));
} 
sub get_variations {
	my ($project) = @_;
	my $id = GenBoStorable::getStoreId( $project->buffer->dbh, $project->id, $project->id,"variations");
	return (GenBoStorable::getStore( $project->buffer->dbh, $id ));
} 
 
 

sub get_cnvs {
	my ($project) = @_;
	my $id = GenBoStorable::getStoreId( $project->buffer->dbh, $project->id, $project->id,"cnvs");
	return (GenBoStorable::getStore( $project->buffer->dbh, $id ));
} 
 
 

exit(0);





sub variations_statistic_gene_test {
	my ( $g, $items, $variations ) = @_;
	my @types = ( "Coding", "UTR", "silent", "Splicing", "Intergenic" );
	
	$items->{p_all} = returnNbPatientstest(@$variations);
	$items->{v_all} = scalar(@$variations);

	foreach my $type (@types) {
		getInfostest( $variations, $items, $type, $g, "" );
	}
	

}    

sub returnNbPatientstest {
	my (@vs) = @_;
	my %uniqPatients;
	foreach my $vv (@vs) {
		map { $uniqPatients{ $_ }++ } split(";",$vv->{allpatients});
	}
	return scalar( keys %uniqPatients );
}



sub getInfostest {
	my ( $array, $items, $type, $gene, $suf ) = @_;	
	my (@tab) = grep { lc( $_->{$gene->{name}."_consequence"}) eq lc($type) } @$array;
	$items->{ "p_" . lc($type) . $suf } = returnNbPatientstest(@tab);
	$items->{ "v_" . lc($type) . $suf } = scalar(@tab);
}


sub parse_somatic {
	my ($file) = @_;
	open(FILE,$file);
	my %ped_bypatient;
	my %ped_byfam;
	while (my $line = <FILE>){
		chomp($line);
		next if $line eq "";
		my $p;
		($p->{fam},$p->{name},$p->{status}) = split(" ",$line);
		$ped_bypatient{$p->{name}} = $p;
	}
	

	
	return (\%ped_bypatient);
}

sub parse_ped {
	my ($file) = @_;
	open(FILE,$file);
	my %ped_bypatient;
	my %ped_byfam;
	while (my $line = <FILE>){
		chomp($line);
		next if $line eq "";
		my $p;
		($p->{fam},$p->{name},$p->{father},$p->{mother},$p->{sex},$p->{status}) = split(" ",$line);
		$stats_by_families->{$p->{fam}}={};
		$ped_bypatient{$p->{name}} = $p;
		push(@{$ped_byfam{$p->{fam}}},$p);
	}
	my $peds;
	
	foreach my $fam (keys %ped_byfam){
		my @child_d;
		my @child_s;
		my $pedigree;
		$pedigree->{fam} = $fam;
		$pedigree->{parents} = [];
		foreach my $ind (@{$ped_byfam{$fam}}){
			if ($ind->{mother} ne "0" || $ind->{father} ne "0"){
				push(@child_d,$ind->{name}) if $ind->{status} == 0;
				push(@child_d,$ind->{name}) if $ind->{status} == 2;
				push(@child_s,$ind->{name}) if $ind->{status} == 1;
			}
			else {
				my $namep = $ind->{name};
				my ($find) = grep {$namep eq $_} @$patients;
				next unless $find;
				push(@{$pedigree->{parents}},$ind->{name});
				$pedigree->{mother} = $ind->{name} if $ind->{sex} == 2;
				$pedigree->{father} = $ind->{name} if $ind->{sex} == 1;
			}
		}
		$pedigree->{child_d} = \@child_d;
		$pedigree->{child_s} = \@child_s;
		
		$pedigree->{samples} = [@child_s,@child_d,@{$pedigree->{parents}}];
		
		push(@$peds,$pedigree);
	}
	return (\%ped_bypatient,$peds);
}


sub load_data {

 my $types_for_pre_filtering = {"evs"=>1,"1000genomes"=>1,"dbsnp"=>1,"intronic"=>1,"1000genomes_1p"=>1,"evs_1p"=>1,"pseudo"=>1,"intergenic"=>1,"dbsnp_none"=>1};
 
 my @array_prefiltering;
@hfilter_type_variations{@filter_type_variations} = ();
foreach my $filter_name (@filter_type_variations) {
	delete $hfilter_param_non_statistic{$filter_name};
	push(@array_prefiltering,$filter_name) if exists $types_for_pre_filtering->{$filter_name};
}

my $string_prefiltering = join(";",sort{$a cmp $b} @array_prefiltering);

my @param_set;


my ($data,$kyoto) = get_data($project,join(";",@param_set));
$data->{prefiltering} = $string_prefiltering;
@filter_type_variations = keys %hfilter_type_variations if $kyoto;

return $data;
}


sub find_gene_ontology { 
	my ($text) = @_;
	
 my $apph = $buffer->get_go_db();
 return unless $apph;
my $genes;
foreach my $search (@$text){ 

 my $terms = $apph->get_terms({search=>"$search*",search_fields=>"name,synonym,definition"});

 foreach my $term (@$terms){
 
 	my $products = $apph->get_deep_products({term=>$term->name});
 	foreach my $g (@$products){
 		push(@{$genes->{$g->symbol}},$term->name);
	
 	}
 	
 }
}

 return $genes;
}

sub byVarIdsForVcf {
	my ($all_data, $patientName) = @_;
	my $hdata;
	foreach my $chr (sort(keys(%$all_data))) {
		foreach my $geneId (keys(%{$all_data->{$chr}->{'genes'}})) {
			foreach my $varId (@{$all_data->{$chr}->{'genes'}->{$geneId}->{data}->{$patientName}}) {
				$hdata->{$varId} = undef;
			}
		}
	}
	return $hdata;
}

sub getVcflines {
	my ($hVarIds, $patientName) = @_;
	my @lVcfLines_out;
	my $nbFound = 0;
	foreach my $file (@{$project->getPatient($patientName)->getVariationsFiles()}) {
		warn $file;
		if ($file =~ /.gz/) { open(FILE, "zcat $file |"); }
		else { open(FILE, "$file"); }
		while(<FILE>) {
			my $line = $_;
			if ($line =~ /#/) { push(@lVcfLines_out, $line); }
			my @lFields = split("\t", $line);
			my @lVarAll = split(',', $lFields[4]);
			my $alreadyWritten;
			foreach my $var_all (@lVarAll) {
				my $id = $lFields[0].'_'.$lFields[1].'_'.$lFields[3].'_'.$var_all;
				$id =~ s/chr//;
				if (exists($hVarIds->{$id})) {
					$nbFound++;
					unless ($alreadyWritten){
						push(@lVcfLines_out, $line);
						$alreadyWritten = 1;
					}
					delete($hVarIds->{$id});
				}
			}
		}
		close(FILE);
	}
	my $patId = $project->getPatient($patientName)->id();
	foreach my $ind (@{$project->getPatient($patientName)->getIndels()}) {
		if (exists($hVarIds->{$ind->id()})) {
				my $hAnnex = $ind->annex();
				if (exists ($hAnnex->{$patId}->{'info_vcf'})) {
				my $id = $ind->id();
				my @lTmp = split('_', $id);
				my $line = 'chr'.$lTmp[0];
				$line .= "\t".$ind->start();
				$line .= "\t.\t".$lTmp[2];
				$line .= "\t".$lTmp[3];
				$line .= "\t.\t.\t.\tGT:AD:DP";
				$line .= "\t0/1" if ($hAnnex->{$patId}->{'he'} eq '1');
				$line .= "\t1/1" if ($hAnnex->{$patId}->{'ho'} eq '1');
				$line .= ":".$hAnnex->{$patId}->{'nb_all_ref'};
				$line .= ",".$hAnnex->{$patId}->{'nb_all_mut'};
				$line .= ":".$hAnnex->{$patId}->{'dp'};
				$line .= "\n";
				push(@lVcfLines_out, $line);
				delete($hVarIds->{$id});
				$nbFound++;
			}
		}
	}
	warn 'Nb Found: '.$nbFound;
	warn 'Nb Left: '.scalar(keys(%$hVarIds));
	@lVcfLines_out = sort(@lVcfLines_out);
	return \@lVcfLines_out;
}
