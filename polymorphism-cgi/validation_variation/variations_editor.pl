#!/usr/bin/perl 
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);

#use lib "/software/polyweb/poly-src/GenBo/lib/obj-nodb/";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag";

use GBuffer;

#use lib "/bip-d/soft/distrib/tabix/latest/perl";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/coverage";
use lib "$Bin/../packages/validation_variation";
use lib "$Bin/../packages/cache";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../GenBo/lib/obj-nodb/polyviewer";
use draw_cnv;
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/utility.pm";
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone retrieve);
use html;
use polyviewer_css;
use Carp;
use strict;
use Data::Dumper;
use GenBo;
use File::Temp qw/ tempfile tempdir /;
#require "$Bin/../GenBo/lib/obj-nodb/GBuffer.pm";
use List::MoreUtils qw{ natatime };
use Getopt::Long;
use Carp;
use Set::Intersection;
use Tabix;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Spreadsheet::WriteExcel;
use POSIX;
use validationQuery;
use QueryValidationAcmg;
use Date::Tiny;
use JSON::XS;
use List::MoreUtils qw{part firstval lastval};
use constant mm => 25.4 / 72;
use constant in => 1 / 72;
use constant pt => 1;
use Time::HiRes qw ( time alarm sleep );
use Storable qw/thaw freeze/;
use LMDB_File qw(:flags :cursor_op :error);
use lib "$Bin/variation_editor_lib";
use hvariant;
use variations_editor_methods;
use annotations_methods;
use max_score;
use MCE::Loop;
use Sys::Hostname;
use CGI::Cache;
use table_dude;
use Digest::MD5::File qw(md5 md5_hex file_md5_hex url_md5_hex file_md5);
use Capture::Tiny ':all';
use PolyviewerVariant;
#use HTML::TableExtract ;
use polyviewer_html;
use HTML::TableExtract ;
#use HTML::Minifier;
#use CHI;
use Storable qw( freeze );
use xls_export;
use session_export;
use session_export_test;
use HTML::Packer;
use Proc::Simple;
use GenBoNoSqlRocksTinyPolyviewerVariant;
##########################################
#VAriaable definition 

######################################
# VERSION
##############################
my $VERSION = "21.09.2021";

my $VERSION_UPDATE_VARIANT_EDITOR = "9dd1dc46f32240dad92ea588d85d9f9d1"; #23/09/2021;
my $dev;
$dev = 1 if $ENV{SERVER_NAME} eq  "10.200.27.103";
my $cgi = new CGI();

my $buffer_polyviewer ={};


my $hgnomad_ac_ho = {
	"1" => "0",
	"2" => "5",
	"3" => "50",
	"4" => "1000",
	"5" => "5000",
	"6" => "-1",
};

my $hgnomad_ac = {
	"1" => "5",
	"2" => "50",
	"3" => "100",
	"4" => "1000",
	"5" => "10000",
	"6" => "-1",
};

my $hsample_dv = {
	"1" => "5",
	"2" => "10",
	"3" => "50",
	"4" => "100",
	"5" => "1000",
	"6" => "-1",
};

my $hsample_dv_ho = {
	"1" => "5",
	"2" => "10",
	"3" => "50",
	"4" => "100",
	"5" => "500",
	"6" => "-1",
};

my $hin_this_run = {
	"1" => 0 ,
	"2" =>10,
	"3" => 25,
	"4" => 50,
	"5" => 75,
	"6" =>100,
};
my @headers =
  ( "Viewer", "igv", "var_name", "trio", "gnomad", "deja_vu" );

my @header_transcripts = (
	"consequence", "enst",     "nm",           "ccds",
	"appris",      "exon",     "nomenclature", "codons",
	"codons_AA",   "polyphen", "sift",         "cadd",
	"revel",       "dbscsnv",  "spliceAI"
);

my @header_transcripts_cnv = ("consequence", "enst", "nm", "ccds", "appris", "start", "end");

  ; #qq{background: #CFD8DC url('https://everything.typepad.com/photos/uncategorized/2008/03/28/crosshatch.png') 0 0 repeat};# #607d8b";#"#607D8B";
my $bgcolor2 = "background-color:#E17F63";
my $dir_tmp2 = tempdir( CLEANUP => 1, DIR => "/tmp/" );

my $host  = hostname; 
my $datestring = localtime();
my $string =  "--> $datestring :: ";
my @params = $cgi->param();
my $tparam;
foreach my $op (@params ){
	$tparam .= $op."=".$cgi->param($op).";";
}

my $print      = $cgi->param('print');
my $purge_gene = $cgi->param('purge_gene');

my $gene_id_filtering;
my $gene_name_filtering = $cgi->param('gene_name');

my $keep_pathogenic     = $cgi->param('keep_pathogenic');
my $filter_quality      = $cgi->param('allele_quality');




################################
open(LOG,">>/data-isilon/tmp/log/".$host.".log");
print LOG $string.$tparam."\n";
close LOG;
################################
################################
$| = 1;

$SIG{INT} = sub {
	warn "I'll die\n";
	exit 1;
};



my $project_name = $cgi->param('project');

my $user         = $cgi->param('user_name');
my $panel_name  = $cgi->param('panel');
my $patient_name = $cgi->param('patients');
my $TEST = $cgi->param('TEST');
my $buffer       = GBuffer->new();

my $project = $buffer->newProjectCache(
	-name        => $project_name,
	-typeFilters => 'individual',
	-cgi_object  => 1
);
#my $clinvar_lmdb =  $project->getChromosome->get_lmdb_database("clinvar",$project->type_public_db);

my $version_db =  $project->public_database_version;
if($version_db>=13){
	$VERSION = $VERSION."-".$version_db."1";
	
}
if($version_db == 15){
	$VERSION = $VERSION."-".$version_db."clinvar";
	$buffer->public_data_version(16);
	$project->public_database_version(16);
	
}

#if($version_db > 17){
#	$VERSION = $VERSION."-".$version_db."clinvar.inh";
#}
else {
	$version_db = undef;
}
my $patient = $project->getPatient($patient_name);
my $fam = $patient->getFamily();
my $print_html = polyviewer_html->new( project => $project, patient => $patient,header=> \@headers,bgcolor=>"background-color:#607D8B" );
$print_html->init();
if ($patient->isMale){
	$VERSION .= "11" ;#if $patient->isMale;
}

unless ($cgi->param('phenotype')){
	if ($project->phenotypes){
		$cgi->param(-name=>'phenotype',-value=>$project->phenotypes->[0]);
	}
}
#
if ($gene_name_filtering){
	my $gene;
	eval {
	 $gene = $project->newGene($gene_name_filtering);
	};
	if ($gene) {
		$gene_name_filtering = $gene->external_name;
	}
	else {
		html::print_cgi_header($cgi);
		print polyviewer_css::css1();
		print polyviewer_css::css2();
		error(" ERROR GENE NAME <HR> WHAT DO YOU MEAN BY : \"".$gene_name_filtering."\"");
		exit(0);
	}
}
my $force;
my $no_cache;
$force=1 unless $project->isGenome();
if ($cgi->param('force') == 1) {
 $dev = 2;	
}

$no_cache = $patient->get_lmdb_cache("w");
my $keys = return_uniq_keys($patient,$cgi); 
my $level_dude = 'high,medium';


my $cache_dude_id = "$level_dude::" . $project_name . "-" . $patient->name . "-" . $VERSION;
my $cache_id = join( ";", @$keys );


if ($project->isDiagnostic){
	$cache_id.="121222";
}
if ($cgi->param('only_DM')){
	$cache_id.="o";
}


my $text = $no_cache->get_cache($cache_id."====");
unless ($text) {
	$cache_id = md5_hex($cache_id);
	$text = $no_cache->get_cache($cache_id);
}
$text = "" if $dev;
$text= "";
my $html_dude = "<!--DUDE-->";
my $cache_icon;

 $cache_icon = qq{<span class="glyphicon glyphicon-floppy-remove" aria-hidden="true" style="text-align:right;font-size:10px;color:red"></span>};
if($text){

		$|= 1;
		my @toto = split("\n<!--SPLIT-->\n",$text);
		my $t = time;
		my $html;
	#	for (my $i = 0;$i<110;$i++){
		#	$html.= $toto[$i];
	#	}
		$html = $text;
		my @toto = split("\n<!--SPLIT-->\n",$text);
		my $t = time;
		my $html;
		#
		for (my $i = 0;$i<@toto;$i++){
			if($project->genome_version_generic =~/HG38/){
				$toto[$i] =~ s/\?dataset\=gnomad_r2_1//g;
			}
			$html.= $toto[$i];
		}
		
		print $html;
		print"<br>".abs(time-$t)."<br>";
		print"<br><div style='float:right;'><small>$cache_icon</small></div><br>"; 
	  	$no_cache->close();
	  	exit(0);
}

$project->cgi_user($user);
$cgi->param(-name=>'user_name',-value=>'');
#my @keys =  $cgi->Vars;


##############################
# PARAMETRES  GLOBALL 
##############################
my $annot = $cgi->param('annot');
my $only_DM = $cgi->param('only_DM');

my $statistics = {};

#                           #
###### Transmission  ########
#                           #



my $vfreq = $cgi->param('frequence');

#                           #
##### #Gnomad and DejaVu ###
#                           #
my $limit_ac = $hgnomad_ac->{ $cgi->param('ac') };
my $limit_ac_ho = $hgnomad_ac_ho->{ $cgi->param('ach') };
my $limit_sample_dv = $hsample_dv->{ $cgi->param('dv') }; ;
my $limit_sample_dv_ho =  $hsample_dv_ho->{ $cgi->param('dv_ho') };


#                           #
######### In this run #######
#                           #

my $in_this_run = $cgi->param('in_this_run');
 unless ($in_this_run){
 	$in_this_run = 100 ;
 }
else {
	$in_this_run = $hin_this_run->{$in_this_run};
	$project->in_this_run_patients();
}


#                           #
######### Annotations #######
#                           #

my ( @annots, @tconsequences );

if ( $annot =~ /\+/ ) {
	@annots = split( '\+', $annot );
}
else {
	@annots = split( " ", $annot );
}
my $maskcoding = 0;
confess() unless $annot;

my $annotation_filer =[];
foreach my $a (@annots) {
	foreach my $this_a ( split( ',', $buffer->get_genbo_annotation_name($a) ) )
	{
#		warn $this_a;
		$maskcoding = $maskcoding | $project->getMaskCoding($this_a);
		push( @tconsequences, $this_a );
	}
}







die() unless $patient;


my $vquery = QueryValidationAcmg->new(
	dbh      => $buffer->dbh,
	database => $project->validation_db()
);
warn "no cache ******************** ".$project->name." ".$patient->name unless $cgi->param('export_xls');           
my $version;
$version->{gencode}->{version} = $project->gencode_version();
$version->{gnomad}->{version} = $project->getChromosome('1')->rocksdb_version('gnomad');
$version->{hgmd}->{version} = $project->getChromosome('1')->rocksdb_version('hgmd');
$version->{cadd}->{version} = $project->getChromosome('1')->rocksdb_version('cadd');
$version->{clinvar}->{version} = $project->getChromosome('1')->rocksdb_version('clinvar');

my $t    = time;
$patient->getLatestValidationStatus($user);
$project->getChromosomes();
$project->validations;
my $ztime =  "";
my $fam  = $patient->getFamily;
my $hgmd = $buffer->queryHgmd->getHGMD($user);

#WARNING

my $stdout = tee_stdout {
	html::print_cgi_header($cgi) unless $cgi->param('export_xls');
	print polyviewer_css::css1() unless $cgi->param('export_xls');
	print polyviewer_css::css2() unless $cgi->param('export_xls');
	my $t = $cgi->param('patients');
	print polyviewer_css::css3($t) unless $cgi->param('export_xls');
};
 
 #print  qq{<div> DEV MODE NO CACHE </div>} if $dev;;
#compute date cache and other 
my $date = date_cache_bam($project);


if ( $project->isDefidiag ) {
	$panel_name = "ACMG-Actionable" if ( $patient->isParent() );
}

#$panel_name="ACMG-Actionable";
my $panel;
$panel = $project->getPanel($panel_name) if $panel_name && lc($panel_name) ne "all" && lc($panel_name) ne "all panels genes";


my @list_genes;
my $hash_genes_panel;
my $start_vector ;
foreach my $c (@{$project->getChromosomes}){
		$start_vector->{$c->name} = $c->getVectorOrigin();
}


	
if ($project->isDiagnostic ){
	my @transcripts_cgi = @{$project->bundle_transcripts() } ;
	foreach my $c (@{$project->getChromosomes}){
		$start_vector->{$c->name} = $c->getNewVector();
	}
#	foreach my $ts (@transcripts_cgi){
#		my $t = $project->newTranscript($ts);
#		my $v = $t->getGene->getVectorPatient($patient);
#		
#		$start_vector->{$t->getChromosome->name} += $v;
#		$hash_genes_panel->{$t->getGene->id} = undef;
#	}
}



if ($panel) {
	$panel->getGenes();
	$hash_genes_panel = $panel->{genes_object};
	my $hPhenoIds;
	foreach my $phenotype ( @{ $panel->getPhenotypes() } ) {
		$hPhenoIds->{ $phenotype->id() } = undef;
	}
	$project->{phenotypes_object} = $hPhenoIds;
}

my $phenotype_name = $cgi->param('phenotype');
if ($phenotype_name) {
	my $hPhenoIds;
	foreach my $panel ( @{ $project->getPanels() } ) {
		foreach my $phenotype ( @{ $panel->getPhenotypes() } ) {
			if ( $phenotype->name() eq $phenotype_name ) {
				$hPhenoIds->{ $phenotype->id() } = undef;
			}
		}
		if ($panel_name && lc($panel_name) eq "all panels genes") {
			$panel->getGenes();
			foreach my $ensg (keys %{$panel->{genes_object}}) {
				$hash_genes_panel->{$ensg} = undef;
			}
		}
	}
	$project->{phenotypes_object} = $hPhenoIds;
}

$t = time;
$buffer->disconnect();

#my $list_saved;
	#warn "here";

##################################
################## GET VECTORS 
##################################
#constructChromosomeVectorsPolyDiagTest($project, $patient,$statistics );
#


my $rocksdb_pv =  GenBoNoSqlRocksTinyPolyviewerVariant->new(mode=>"r",patient=>$patient,project=>$project,print_html=>$print_html);
$t = time;
my (  $list, $id_by_genes_id) =  getListVariantsFromDuckDB($project,$patient,$statistics);
	$ztime .= ' vectors:' . ( abs( time - $t ) );
	warn "duckdb :: ".$ztime if $print;
	
	$t = time;

##################################
################## GET GENES 
##################################

my ($genes) = run_annnotations( $list, $id_by_genes_id);

export_xls($patient, $genes) if $cgi->param('export_xls');

$ztime .= ' ' . scalar(@$genes) . '_genes:' . ( abs( time - $t ) );
warn "annot ++ == ++ " . abs( time - $t ) if $print;

$statistics->{genes} = scalar(@$genes);




#$statistics->{variations} = $statistics->{variants};
$project->buffer->dbh_reconnect();
unless ( $project->buffer->getQuery->isUserMagic($user) ) {
	$ztime = undef;
}


$project->buffer->dbh_reconnect();


#warn $genes->[0]->{js_id};
my $sct = time;

warn "end max_score :".abs(time -$sct);
	my $vgenes =[];
	my $nb2 = 0;
	my $limit =3000; 
	my $hchr; 

foreach my $g (sort{$b->{max_score} <=> $a->{max_score}} @$genes)	{
		if ($gene_name_filtering) {
			if ($gene_name_filtering eq $g->{name} or $gene_id_filtering eq $g->{id}) { 
				$hchr->{$g->{chr_name}} ++;
				push(@$vgenes,$g);
			}
		}
		else {
			$hchr->{$g->{chr_name}} ++;
			push(@$vgenes,$g);
		}
		$nb2 ++;
		last if $nb2 > $limit;
	}
$genes = $vgenes;	
$genes = [] unless $genes;

$t = time;

$ztime .= ' polycyto:' . ( abs( time - $t ) );
$t     = time;
my $stdout_nav_bar = tee_stdout {
	print"<br><div style='float:right;'>$cache_icon</div><br>"; 
    update_variant_editor::printNavBar( $patient, $genes, $statistics, $version,$date, $user, $ztime,update_variant_editor::compose_string_filtering($cgi) );
	#update_variant_editor::print_hotspot( $patient, $panel );
};
my $stdoutcnv;

if (not $patient->isGenome() ) {
	$stdoutcnv = tee_stdout {
		update_variant_editor::print_cnv_exome( $patient, $level_dude, $panel );
	};
}

#warn "end";
	
	
$t     = time;
my $stdout_end = tee_stdout {
	
	warn "genes:".scalar(@$genes);
	if (@$genes){
	$genes = refine_heterozygote_composite_score_fork( $project, $genes,$hchr ,$buffer_polyviewer) ;
	}
	else {
		if ($gene_name_filtering ) {
		print qq{<div class="knockout">No Variation Found for gene $gene_name_filtering</div>};
	}
	else {
		print qq{<div class="knockout">No Variation Found</div>};
	}
	}
	warn "hetero " . abs( time - $t ) if $print;
	$ztime .= ' hetero:' . ( abs( time - $t ) );
	print $ztime;

};
#warn "--> $cache_id <-- put";
$no_cache->put_cache_text($cache_id,$stdout.$stdout_nav_bar.$stdoutcnv.$stdout_end,2400) ;#unless $dev;
#$no_cache->get_cache($cache_id);
$no_cache->close();
exit(0);


sub error {
	my ($text) = @_;
	print qq{</div>};
	print qq{<div class="knockout" style="background-color:red" >$text</div>}
	  ;    # if scalar(@$genes) == 0;
	print
	  qq{<div id ="logo" > contact : bioinformatique\@users-imagine.fr</div>}
	  ;    # if scalar(@$genes) == 0;
	warn "error !!!!!!!!! $text";
	warn $TEST;
	die() if $TEST;
	exit(0);
}



sub calculate_max_score {
	my ( $project, $list,$no ) = @_;
	my $fork = 1;
	my $nb        = int( scalar(@$list) / ($fork) +1 );
	my $iter      = natatime( $nb,  @$list );
	my $vid        = time;
		my $genes  = max_score::calculate($project,$patient,$list,$no,$rocksdb_pv);

	return $genes;
	
}




sub refine_heterozygote_composite_score_fork {
		my ( $project, $genes,$hchr ) = @_;
	$t = time;
	$| =1;
	print qq{<div style="display: none">};
	print "refine ";
	my $fork      = scalar (keys %$hchr);
	$fork = 1;
	my $res_genes = [];
	my $vid        = 0;
	my $hrun;
	
	#$t = time;
	my $current = 1;
	my $vres;
	my @res_final;
	my @toto;
	my $wtime = 0;
	my $maxgene =100;
	my $ngene =0;
	
	my $final_polyviewer_all ;
	#if ($project->isRocks){
		my $diro = $project->rocks_directory();
	
	
	
	
	#TODO: essayer d'integrer XLS Store variant ici pour le forker
	
	$project->buffer->dbh_deconnect();
	
	#delete $project->{validations_query1};
		my $t   = time;
		my $res;
		$project->buffer->dbh_reconnect();
		( $res->{genes}, $res->{total_time} ) = variations_editor_methods::refine_heterozygote_composite( $project,$print_html, \@$genes, $vid,$rocksdb_pv);
		$res->{run_id} = $vid;
		warn "===============================";
	
	
	$project->buffer->dbh_reconnect();
	
	warn "end hetero " if $print;
	warn "....";
	print qq{</div>};
	warn "....";
	my $nb_genes = scalar(@{$res->{genes}});
	
	foreach my $g (@{$res->{genes}}){
		print $g->{out} . "\n";
		delete  $g->{out};
		#last if $g->{max_score} < 5 && $nb_genes > 300;
	}
	
	#$final_polyviewer_all->deactivate_cache();
	warn "end";
	error("Hey it looks like you're having an error  !!! ") if scalar keys %$vres;	
	return ;
	exit(0);
	print qq{</div>};
	return $res_genes;
}


sub run_annnotations {
	my ( $list,$id_by_genes_id ) = @_;
	unless ($cgi->param('export_xls')) {
		print qq{<div style="display: none">};
		print "annotations";
	}
	my $hgenes;
	$project->getChromosomes();
	$project->buffer->dbh_deconnect();
	$project->disconnect;	
	my $final_polyviewer_all;
	 $final_polyviewer_all = GenBoNoSqlRocks->new(cache=>1,dir=>$project->rocks_directory."/patients/",mode=>"r",name=>$patient->name,vmtouch=>1);
	 my $t   = time;
		$final_polyviewer_all->prepare($list);
		my $nb = 0;
		my $res,
		my $agenes;
		foreach my $id (@$list) {
			my $debug;
			$nb++;
			unless ($cgi->param('export_xls')) {
				print "." if $nb % 100 == 0;
			}
			my $hg = $final_polyviewer_all->get($id);
				foreach my $gene ( @{$hg->{array}}) {
					my $gid = $gene->{id};					
					next unless exists $id_by_genes_id->{$gid}->{$id};
					unless ( exists $hgenes->{$gid} ) {
						$hgenes->{$gid} = $gene;
						$hgenes->{$gid}->{score_father_mother} = $hgenes->{$gid}->{score_mother} + $hgenes->{$gid}->{score_father};
						$hgenes->{$gid}->{max_score} = $hgenes->{$gid}->{score} + ( $hgenes->{$gid}->{score_father_mother} >  $hgenes->{$gid}->{score_biallelic} ? $hgenes->{$gid}->{score_father_mother} : $hgenes->{$gid}->{score_biallelic} );
						$hgenes->{$gid}->{max_score} -=( $hgenes->{$gid}->{denovo_rare} * 0.3 ) if $hgenes->{$gid}->{denovo_rare} > 2;
						$hgenes->{$gid}->{max_score} += 2 if $hgenes->{$gid}->{cnv_dup} > 1 or $hgenes->{$gid}->{cnv_del} > 1;
					}
					else {
						push(@{ $hgenes->{$gid}->{variants} }, @{ $gene->{variants} } );
					}
					foreach my $k ( keys %{ $gene->{variant} } ) {
						$hgenes->{$gid}->{all_variants}->{$k} = $gene->{variant}->{$k};
						$hgenes->{$gid}->{all_vector_ids}->{$k} = $gene->{vector_ids}->{$k};
					}
				}
			}
			

	print qq{</div>} unless ($cgi->param('export_xls'));
	$project->buffer->dbh_reconnect();
	return calculate_max_score($project,[ values %$hgenes ],$final_polyviewer_all);
}






sub construct_panel_vector {
	my ( $panel, $hashVector ) = @_;
	my $hashintspan;
	foreach my $gene ( @{ $panel->getGenes } ) {
		my $chr = $gene->getChromosome();
		next if $chr->size_vector() == 0;

		#warn $chr->name()." ".$chr->size_vector() if $chr->name eq "MT";
		$hashVector->{ $chr->name } = $gene->getVectorOrigin
		  unless exists $hashVector->{ $chr->name };
		$hashVector->{ $chr->name } += $gene->getVectorOrigin;

	}
}






sub getListVariantsFromDuckDB {
	my ( $project, $patient,$statistics ) = @_;
	
my $h_transmissions = {
    solo          => 1 << 0,  # 2^0 = 1
    father        => 1 << 1,  # 2^1 = 2
    mother        => 1 << 2,  # 2^2 = 4
    both          => 1 << 3,  # 2^3 = 8
    is_parent     => 1 << 4,  # 2^4 = 16
    recessif      => 1 << 5,  # 2^5 = 32
    dominant      => 1 << 6,  # 2^6 = 64
    denovo        => 1 << 7,  # 2^7 = 128
    strict_denovo => 1 << 8,  # 2^8 = 256
    error         => 1 << 9,  # 2^9 = 512
    mosaic         => 1 << 10,  # 2^10 = 1024
    uniparental     => 1 << 11,  # 2^11 = 2048
    
};
	my $filter_transmission;
	my $xtime =time;
	my $list_genes;
	my $mask_transmission  |= $h_transmissions->{solo};
	$mask_transmission   |= $h_transmissions->{is_parent};
	$mask_transmission   |= $h_transmissions->{denovo} if $cgi->param('denovo');	
	$mask_transmission   |= $h_transmissions->{strict_denovo} if $cgi->param('strict_denovo');
	$mask_transmission   |= $h_transmissions->{recessive} if $cgi->param('recessive');
	$mask_transmission   |= $h_transmissions->{uniparental} if $cgi->param('recessive');
	$mask_transmission   |= $h_transmissions->{mosaic} if $cgi->param('mosaic');
	$mask_transmission   |= $h_transmissions->{both} if $cgi->param('both');
	$mask_transmission   |= $h_transmissions->{mother} if $cgi->param('xor_mother');
	$mask_transmission   |= $h_transmissions->{father} if $cgi->param('xor_father');
	$mask_transmission   |= $h_transmissions->{mother} if $cgi->param('xor');
	$mask_transmission   |= $h_transmissions->{father} if $cgi->param('xor');
	
	
	my $hashVector_panel = {};
	my $hashVector       = {};

	#construct_panel_vector( $panel, $hashVector_panel ) if $panel;
	my $list_transcript;
	my $trio = $patient->getFamily->isTrio;
	my $list_variants =[];
	my $hash_variants_DM = {};
	



	delete $project->{rocks};

	my $finalVector = {};
	
	#######################
	# FILERING patient trans and ratio
	#######################
	my $suffix = "patient_".$patient->id;
	#CAST(patient_55661_transmission AS INTEGER)
	my $asql_patient = [];
	my $alt = 0;
	if($cgi->param('alt') ){
		$alt = $cgi->param('alt');
	}
	
	push(@$asql_patient ,$suffix."_alt >".$alt);#." and ".$suffix."_transmission  & $mask_transmission <> 0" ;
	if($cgi->param('ratio') ){
		push(@$asql_patient ,$suffix."_ratio > ".$cgi->param('ratio'));
	}
	
	if($cgi->param('ref') ){
		push(@$asql_patient ,$suffix."_ref > ".$cgi->param('ref'));
	}
	push(@$asql_patient ,$suffix."_transmission  & $mask_transmission <> 0");
	
	my $sql_patient = "(".join(" and ",@$asql_patient).")";
	#push(@column_patient,"patient_".$c."_ref");
	#push(@column_patient,"patient_".$c."_alt");
	#push(@column_patient,"patient_".$c."_ratio");
	#push(@column_patient,"patient_".$c."_type");
	#push(@column_patient,"patient_".$c."_transmission");
	
	#######################
	# FILERING Genes
	#######################
	#my @column_gene = ("gene_name","gene_mask");
	
	my $sql_gene = "gene_name != '-' "."and gene_mask & $maskcoding <> 0";
	my $gene;
	if ($gene_name_filtering) {
		$gene = $project->newGene($gene_name_filtering);
		$gene_id_filtering = $gene->id();
		$sql_gene .= "gene_name = '".$gene->id."' and gene_mask & maskcoding <> 0";
	}
	$sql_gene = "(".$sql_gene.")";
	
	#######################
	# FILERING frequence 
	#######################
	#my @column_frequences_gnomad = ("gnomad_ac","gnomad_an","gnomad_min_pop_name","gnomad_min_pop_freq","gnomad_max_pop_name","gnomad_max_pop_freq","gnomad_ho","getGnomadAC_Male");
	#my @column_frequences_dejavu = ("other_projects","other_patients","other_patients_ho","similar_projects","similar_patients","similar_patients_ho","in_this_run_patients");
	my $asql_frequence;


	push(@{$asql_frequence}, "variant_gnomad_ac < $limit_ac ") if $limit_ac > 0 ;
	push(@{$asql_frequence}, "variant_gnomad_ho < $limit_ac_ho ")  if $limit_ac_ho >0 ;;
	push(@{$asql_frequence}, "variant_other_patients < $limit_sample_dv ") if $limit_sample_dv > 0;
	push(@{$asql_frequence}, "variant_other_patients_ho < $limit_sample_dv_ho ") if $limit_sample_dv_ho > 0; 
	my $sql_frequence = join(" and ",@$asql_frequence);
	$sql_frequence = "(".$sql_frequence.")";
	my $sql_only_dm;
	if ($cgi->param('only_DM')){
		$sql_frequence = "";
		$sql_gene = "";
		$sql_only_dm = qq{(variant_isDM  = 1 or variant_isClinvarPathogenic =1)};
	}
	
	my $sql_where_and;
	
	push(@$sql_where_and,$sql_frequence) if $sql_frequence;
	push(@$sql_where_and,$sql_gene) if $sql_gene;
	push(@$sql_where_and,$sql_patient) if $sql_patient;
	push(@$sql_where_and,$sql_only_dm) if $sql_only_dm;
	my $where = join (" and ",@$sql_where_and);
	
	if ($cgi->param('keep_pathogenic') == 1 ){
		#$where .= "or (variant_keepPathogenic = 1) ";
		
	}
	#!!!!!!!!!!!!!!!!
	
	#!!!!!!!!!!!!!!!!
	#add in this run 
	#add only DM
	#add only DM
	#add vlocal 
	#!!!!!!!!!!!!!!!!
	#!!!!!!!!!!!!!!!!
	#get_join_parquet($project,$sql_frequence,$sql_patient,$sql_gene,$suffix);
	return get_rocksdb_mce_polyviewer_variant($project,$where,$suffix);
#
 	
	return ( $finalVector, $list_variants, $hash_variants_DM,$list_genes);
}


sub get_rocksdb_mce_polyviewer_variant {
	my ($project,$where,$suffix) = @_;
	warn "rocksdb  ";
	my $parquet = $project->parquet_cache_variants();
	#$parquet = "/data-beegfs/tmp/new/NGS2025_09289.variants.parquet";
	my $dir_parquet = $project->parquet_cache_dir;
	my $diro = $project->rocks_directory();
	
	my $sql =qq{select variant_index,variant_rocksdb_id,gene_name from '$parquet' where  $where ; };

	my $cmd = qq{duckdb -json -c "$sql"};
	warn $cmd;
 	my $t = time;
 	my $res =`$cmd`;
	my $array_ref  = decode_json $res;
	my	$list_variants = [];
 	my	$hash_variants_DM = {};
 	
 	my $nbv = 0;
 	 warn "sql ".abs(time -$t);
 	 $t =time;
 	 my $id_by_genes_id;
 	 my %ids ;
 	foreach my $a (@$array_ref){
 		my ($c,$b) = split("!",$a->{variant_rocksdb_id});
 		my $index = $a->{variant_index};
 		$rocksdb_pv->add_index($c,$index);
 		$id_by_genes_id->{$a->{gene_name}}->{$a->{variant_index}} ++;
 		#$id_by_patients_id->{$c}->{$index.":".$patient->id} ++;
 		$nbv ++;
 	
 	}
	return ($rocksdb_pv->indexes,$id_by_genes_id);
}



sub  date_cache_bam {
	my($project) = @_;
	#my $cno = $project->getChromosomes()->[0]->lmdb_hash_variants("r");



my $date;
#( $date->{cache} ) = utility::return_date_from_file( $cno->filename );
( $date->{bam} )   = utility::return_date_from_file( $patient->getBamFileName );
return $date;
}

sub export_xls {
	my ($patient, $genes) = @_;
	my @lVarObj;
	my $project = $patient->getProject();
	foreach my $h_gene (@$genes) {
		foreach my $var_id (keys %{$h_gene->{all_variants}}) {
			my $var = $project->_newVariant($var_id);
			push(@lVarObj, $var);
		}
	}
	my $xls_export = new xls_export();
	$xls_export->title_page('PolyViewer_'.$patient->name().'.xls');
	$xls_export->store_variants_infos(\@lVarObj, $project, $project->getPatients());
	my ($list_datas_annotations) = $xls_export->prepare_generic_datas_variants();
	my ($list_datas_annotations_cnvs) = $xls_export->prepare_generic_datas_cnvs();
	$xls_export->add_page_merged('Variants Merged', $xls_export->list_generic_header(), $list_datas_annotations);
	$xls_export->add_page('Variants Not Merged', $xls_export->list_generic_header(), $list_datas_annotations);
	if (scalar @$list_datas_annotations_cnvs > 0) {
		$xls_export->add_page('Cnvs', $xls_export->list_generic_header_cnvs(), $list_datas_annotations_cnvs);
	}
	$xls_export->export();
	exit(0);
}

sub save_session_for_test {
	my ($patient, $genes) = @_;
	my $session_id = $cgi->param('session_test_id') if $cgi->param('session_test_id');
	confess("\n\nNo session_test_id param. Die.\n\n") unless $session_id;
	
	my $session_base_name = $cgi->param('session_base_name') if $cgi->param('session_base_name');
	confess("\n\nNo session_base_name param. Die.\n\n") unless $session_base_name;
	
	my @lVarObj;
	my $project = $patient->getProject();
	foreach my $h_gene (@$genes) {
		foreach my $var_id (keys %{$h_gene->{all_variants}}) {
			my $var = $project->_newVariant($var_id);
			push(@lVarObj, $var);
		}
	}
	my $xls_export = new xls_export();
	$xls_export->store_variants_infos(\@lVarObj, $project, $project->getPatients());
	my ($list_datas_annotations) = $xls_export->prepare_generic_datas_variants();
	my ($list_datas_annotations_cnvs) = $xls_export->prepare_generic_datas_cnvs();
	
	my $session_test;
	$session_test = new session_export() if ($session_base_name eq 'dev');
	$session_test = new session_export_test() if ($session_base_name eq 'prod');
	$session_test->load_session($session_id);
	$session_test->save( $session_base_name.'_annot_var', $list_datas_annotations );
	$session_test->save( $session_base_name.'_annot_cnvs', $list_datas_annotations_cnvs );
	
	exit(0);
}

sub update_score_clinvar {
	my ($id) = @_;
	my $vh = $project->_newVariant($id);
	return 1 if $vh->isClinvarPathogenic();
	return 0;
	
}



sub return_uniq_keys {
my ($patient,$cgi) = @_;
	
my %hkeys = $cgi->Vars;
my @keys;
my $string;
foreach my $k  (sort {$a cmp $b} keys %hkeys){
	next if $k =~ /force/;
	next if $k =~ /user/;
	push(@keys,"$k");
	my $c = $hkeys{$k};
	$c =~ s/\"//g;
	$c =~ s/\+/ /g;
	push(@keys,$c);
}

$project->validations_query(1);


my @key2;
foreach my $chr  (@{$project->getChromosomes}){
		my $no = $chr->lmdb_polyviewer_variants( $patient, "r" );
		my @st = (stat($no->filename));
		 push(@key2, ($st[9].$st[11].$st[12]));
		my $no2 = $chr->lmdb_polyviewer_variants_genes( $patient, "r" );
		@st = (stat($no2->filename));
		 push(@key2, ($st[9].$st[11].$st[12]));
}
#warn join("-",@key2);
push(@keys,@key2);

#push(@keys,file_md5_hex($Bin."/variations_editor.pl") );
my @key3;
push(@key3,$VERSION);
push(@key3,$VERSION_UPDATE_VARIANT_EDITOR );
push(@keys,@key3);

warn join("-",@key3);
my $stv = $patient->get_string_validations();
unless ($patient->isGenome ) {
	$stv .= ':::'.$patient->get_string_identification();
}
#warn $stv;
#keep compatibilty 
if  ($stv ){
push(@keys,"validation".":".md5_hex($stv));
}
else{
push(@keys,encode_json ({}));
push(@keys,encode_json ({}));
push(@keys,encode_json ({}));
#push(@keys,encode_json ($h2));	
}




return \@keys;
}

sub  date_cache_bam {
	my($project) = @_;
	
	#my $cno = $project->getChromosomes()->[0]->lmdb_hash_variants("r");

my $align_file = $patient->getAlignFileName;
my $type_align = 'bam';
$type_align = 'cram' if $patient->isCram();

my $date;
#( $date->{'cache'} ) = "toto";#utility::return_date_from_file( $cno->filename );
( $date->{'align-'.$type_align} )   = utility::return_date_from_file( $align_file );
return $date;
}

sub get_vector_from_duckdb {
	my ($patient,$chr,$limit) =@_;
	my $vector = $chr->getNewVector();
	my $parquet = $chr->project->parquet_cache_variants();
	my $col = "patient_".$patient->id."_ratio";
	my $v2 = $chr->ucsc_name;
	my $sql =qq{select variant_vector_id from '$parquet' where  $col > $limit and variant_chromosome='$v2'};
 	my $cmd = qq{duckdb   -json -c "$sql"};
 	my $res =`$cmd`;
 	return $vector unless $res;
 	my $array_ref  = decode_json $res;
 	foreach my $a (@$array_ref){
 		my $z = $a->{variant_vector_id};
 		 $vector->Bit_On($z);
 	}
 	return $vector;
}
