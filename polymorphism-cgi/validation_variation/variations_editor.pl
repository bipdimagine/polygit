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

use Sys::Hostname;
use CGI::Cache;
use table_dude;
use Digest::MD5::File qw(md5 md5_hex file_md5_hex url_md5_hex file_md5);
use Capture::Tiny ':all';
use PolyviewerVariant;
#use HTML::TableExtract ;
use polyviewer_html;
use HTML::TableExtract ;

#use CHI;
use Storable qw( freeze );
use xls_export;
use session_export;
use session_export_test;
#use Proc::Simple;

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


#my $myproc = Proc::Simple->new();        # Create a new process object
 
#$myproc->start("$Bin/vmtouch.pl project=".$cgi->param('project')." patient=".$cgi->param('patients'));
#warn "$Bin/vmtouch.pl project=".$cgi->param('project')." patient=".$cgi->param('patients');
#die();
#allele count
my $hgnomad_ac_ho = {
	"1" => "gnomad_ho_ac_5",
	"2" => "gnomad_ho_ac_10",
	"3" => "gnomad_ho_ac_50",
	"4" => "gnomad_ho_ac_100",
	"5" => "gnomad_ho_ac_1000",
	"6" => "gnomad_ho_ac_all",
};

my $hgnomad_ac = {
	"1" => "gnomad_ac_5",
	"2" => "gnomad_ac_10",
	"3" => "gnomad_ac_50",
	"4" => "gnomad_ac_100",
	"5" => "gnomad_ac_1000",
	"6" => "gnomad_ac_all",
};

my $hsample_dv = {
	"1" => "sdv_5",
	"2" => "sdv_10",
	"3" => "sdv_50",
	"4" => "sdv_100",
	"5" => "sdv_1000",
	"6" => "sdv_all",
};

my $hsample_dv_ho = {
	"1" => "shodv_5",
	"2" => "shodv_10",
	"3" => "shodv_50",
	"4" => "shodv_100",
	"5" => "shodv_500",
	"6" => "shodv_all",
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
  ( "varsome", "igv", "alamut", "var_name", "trio", "gnomad", "deja_vu" );

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


my $limit_ratio         = $filter_quality;
my @array_ratio = (-1,10,20,40,80,100);
my $limit_ratio1 = firstval { $_ >= $filter_quality } @array_ratio;
my $limit_ratio2 = lastval { $_ <= $filter_quality } @array_ratio;

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
#my $no_cache = GenBoNoSqlLmdb->new(
#			dir         => '/data-beegfs/tmp/lmdb',
#			mode        => "w",
#			name        => "$project_name",
#			is_compress => 1,
#			#vmtouch     => $self->buffer->vmtouch
#);
#system("chmod a+w ".$no_cache->filename);

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

my $text = $no_cache->get_cache($cache_id);
#$dev=1;
$text = "" if $dev;
#$text = "";
my $html_dude = "<!--DUDE-->";
my $cache_icon;

 $cache_icon = qq{<span class="glyphicon glyphicon-floppy-remove" aria-hidden="true" style="text-align:right;font-size:10px;color:red"></span>};
if($text){
		my $dataset = "gnomad_r2_1";
		$dataset = "gnomad_r4" if $project->getVersion() =~ /HG38/;
	
		$text =~ s/onClick='editor\(1,2\);'/onClick='load_polyviewer_export_xls\(1\);'/g;
		
		my $regex = qq{href\='https:\/\/gnomad\.broadinstitute\.org\/variant\/(.+)' target};
		$text =~ s/\?dataset\=gnomad_r2_1//g;
		$text =~ s/$regex/href\='https:\/\/gnomad\.broadinstitute\.org\/variant\/$1\?dataset\=$dataset' target/g;
		
		my $regex2 = qq{https:\/\/gnomad\.broadinstitute\.org\/region\/([0-9]+-[0-9]+)};
		$text =~ s/$regex2/https:\/\/gnomad\.broadinstitute\.org\/region\/$1\?dataset\=$dataset/g;
		
		my $regex3 = qq{https:\/\/gnomad\.broadinstitute\.org\/gene\/([0-9A-Za-z]+)};
		$text =~ s/$regex3/https:\/\/gnomad\.broadinstitute\.org\/gene\/$1\?dataset\=$dataset/g;
		
	 	$cache_icon = qq{<span class="glyphicon glyphicon-floppy-saved" aria-hidden="true" style="text-align:right;font-size:10px;color:green"></span>};
	 	
		my $regexp_varsome = qq{varsome.com\/variant\/hg19\/};
		$text =~ s/$regexp_varsome/varsome.com\/variant\/hg38\//g if $project->getVersion() =~ /HG38/;
	 	$|= 1;
		my @toto = split("\n<!--SPLIT-->\n",$text);
		my $t = time;
		my $html;
		for (my $i = 0;$i<110;$i++){
			$html.= $toto[$i];
		}
		print $html;
		print"<br>".abs(time-$t)."<br>";
		print $text;
  		
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
my $limit_ac;
my $limit_ac_ho;
my $limit_sample_dv;
my $limit_sample_dv_ho;

unless ( $hgnomad_ac_ho->{ $cgi->param('ach') } ) {
	$limit_ac_ho        = $hgnomad_ac_ho->{ $cgi->param('ac') };
	$limit_ac           = $limit_ac_ho;
	$limit_sample_dv    = $hsample_dv->{ $cgi->param('dv') };
	$limit_sample_dv_ho = $hsample_dv_ho->{ $cgi->param('dv_ho') };
	

}
else {
	$limit_ac = $hgnomad_ac->{ $cgi->param('ac') };
	warn " default value for ac !!!!!!!!!!!!!!!!!!" unless $limit_ac;
	$limit_ac           = "gnomad_ho_ac_100" unless $limit_ac;
	$limit_ac_ho        = $hgnomad_ac_ho->{ $cgi->param('ach') };
	$limit_sample_dv    = $hsample_dv->{ $cgi->param('dv') };
	$limit_sample_dv_ho = $hsample_dv_ho->{ $cgi->param('dv_ho') };
}

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
my $v = $project->public_database_version;
$version->{gencode}->{version} = "45";#$project->get_public_data_version("gencode");
$version->{gnomad}->{version} = "-";#-$project->get_public_data_version("gnomad");
$version->{hgmd}->{version} = "-";#$project->get_public_data_version("hgmd");
$version->{cadd}->{version} = "-";#$project->get_public_data_version("cadd");
$version->{clinvar}->{version} = "-";#$project->get_public_data_version("clinvar");
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



##################################
################## GET VECTORS 
##################################
#constructChromosomeVectorsPolyDiagTest($project, $patient,$statistics );
	warn "start";
	
	my ( $vectors, $list, $hash_variants_DM,$list_genes) = constructChromosomeVectorsPolyDiagFork( $project, $patient,$start_vector,$statistics);
	$ztime .= ' vectors:' . ( abs( time - $t ) );
	warn $ztime if $print;
$t = time;

##################################
################## GET GENES 
##################################

my ($genes) = fork_annnotations( $list, [], $maskcoding,$vectors);
save_session_for_test($patient, $genes) if $cgi->param('test_mode');
export_xls($patient, $genes) if $cgi->param('export_xls');

warn "------------";
$ztime .= ' ' . scalar(@$genes) . '_genes:' . ( abs( time - $t ) );
warn "annot " . abs( time - $t ) if $print;

$statistics->{genes} = scalar(@$genes);






#$statistics->{variations} = $statistics->{variants};
$project->buffer->dbh_reconnect();
unless ( $project->buffer->getQuery->isUserMagic($user) ) {
	$ztime = undef;
}


$project->buffer->dbh_reconnect();


warn "max_score";
#warn $genes->[0]->{js_id};
my $sct = time;
$genes = calculate_max_score($project,$genes);
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
	#warn $genes->[0]->{js_id};

$t = time;




print print_hotspot($patient);

$ztime .= ' polycyto:' . ( abs( time - $t ) );
$t     = time;
my $stdout_nav_bar = tee_stdout {
	print"<br><div style='float:right;'>$cache_icon</div><br>"; 
    update_variant_editor::printNavBar( $patient, $genes, $statistics, $version,$date, $user, $ztime,update_variant_editor::compose_string_filtering($cgi) );
};
my $stdoutcnv;
unless ( $patient->isGenome() ) {
		 $stdoutcnv = $no_cache->get_cache($cache_dude_id);
		 
		if ($stdoutcnv){
			print $stdoutcnv;
		}
		else {
			$stdoutcnv = tee_stdout {
		update_variant_editor::print_cnv_exome( $patient, $level_dude, $panel );
		};
			 $no_cache->put_cache_text($cache_dude_id,$stdoutcnv,2400) ;#unless $dev; 
		}
}
	
	
#$myproc->kill();
$t     = time;
#my $exit_status = $myproc->wait(); 
#warn "ok *** $exit_status ";
#exit(0);
my $stdout_end = tee_stdout {
#warn "genes:".scalar(@$genes);
if (@$genes){
	warn "refine";
$genes = refine_heterozygote_composite_score_fork( $project, $genes,$hchr ) ;
#warn "genes fin:".scalar(@$genes);

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

$no_cache->put_cache_text($cache_id,$stdout.$stdout_nav_bar.$stdoutcnv.$stdout_end,2400) ;#unless $dev;
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
	my ( $project, $list ) = @_;
	my $fork = 1;
	my $nb        = int( scalar(@$list) / ($fork) +1 );
	my $iter      = natatime( $nb,  @$list );
	my $vid        = time;
	my $hrun;
	my $pm        = new Parallel::ForkManager($fork);
	#$t = time;
	
	my $final =[];
	
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $h ) = @_;

			unless ( defined($h) or $exit_code > 0 ) {
				print
				  qq|No message received from child process $exit_code $pid!\n|;
				die();
				return;
			}
			my $id = $h->{run_id};
			delete $hrun->{ $h->{run_id} };
			push(@$final,@{$h->{genes}});
		}
	);
	
	#TODO: essayer d'integrer XLS Store variant ici pour le forker
	$project->disconnect();
	
	
	#$project->buffer->dbh_deconnect();
	
	#delete $project->{validations_query1};
	while ( my @tmp = $iter->() ) {
		$vid++;
		$hrun->{$vid}++;
		my $pid = $pm->start and next;
			my $res;
			$project->disconnect;
			my $genes  = max_score::calculate($project,$patient,\@tmp);
			$project->disconnect;
			$res->{run_id} = $vid;
			$res->{genes} = $genes;
		$pm->finish( 0, $res );
	}
	$pm->wait_all_children();
	return $final;
	
}


sub calculate_max_score_toto {
	my ( $project, $list ) = @_;
	$project->buffer->close_lmdb();
	
	my $mother = $patient->getFamily->getMother();
	my $father = $patient->getFamily->getFather();


	my $hno;
	my $tsum = 0;
	my $t    = time;
	my $xp   = 0;

	my $total_time = 0;
	my $ztotal     = 0;

	my $h_all_variants_validations = $patient->getProject->validations();
	my $no_dude                    = $patient->getGenesDude();
	$no_dude = undef if $project->isGenome();
	my $final_polyviewer_all = GenBoNoSqlRocks->new(dir=>$project->rocks_directory."/patients/",mode=>"r",name=>$patient->name);
	my @ids = map{$_->{id}} @{$list};
	$final_polyviewer_all->prepare(\@ids);
	foreach my $hgene (@$list) {
		my $gid = $hgene->{id};
		my $debug ;
		my $gene = $project->newGenes( [$gid] )->[0];
		my ( $n, $chr_name ) = split( "_", $gid );
		my $chr = $gene->getChromosome();
		$hgene->{chr_name} = $chr_name;
		
		if ( $no_dude && -e $no_dude->filename ) {
			$hgene->{level_dude} = $no_dude->get($gid);
		}
		else {
			$hgene->{level_dude} = -1;
		}
		my $global_gene = $final_polyviewer_all->get($gid); #$chr->get_polyviewer_genes($patient,$gid);
		foreach my $k ( keys %{$global_gene} ) {
			next if $k eq "penality";
			next if $k eq "denovo_rare";
			$hgene->{$k} = $global_gene->{$k};
		}
		$hgene->{score} = $gene->score;
		my $class;
		$class->{biallelic} = [];
		$class->{mother}    = [];
		$class->{father}    = [];
	
		foreach my $k ( keys %{ $hgene->{all_variants} } ) {
			if ($version_db){
				$hgene->{score} += update_score_clinvar($k);
			}
			#my $pub = $db->get_with_sequence($self->start,$self->alternate_allele);
			if ( exists $h_all_variants_validations->{ $gid . '!' . $k } ) {
				my $score_validation = $h_all_variants_validations->{ $gid . '!' . $k }->[0]->{validation};
				$hgene->{score} += 0.5 if ( $score_validation == 3 );
				$hgene->{score} += 2   if ( $score_validation == 4 );
				$hgene->{score} += 3   if ( $score_validation == 5 );
			}
			
			push(@{ $class->{biallelic} },$hgene->{all_variants}->{$k}->{score}) if exists $hgene->{all_variants}->{$k}->{biallelic};
			push( @{ $class->{mother} }, $hgene->{all_variants}->{$k}->{score} ) if exists $hgene->{all_variants}->{$k}->{mother};
			push( @{ $class->{father} }, $hgene->{all_variants}->{$k}->{score} ) if exists $hgene->{all_variants}->{$k}->{father};
		}
		if (scalar( @{ $class->{mother} } ) > 0 && scalar( @{ $class->{father} } ) == 0 && exists $hgene->{father}->{id} ) {
			my $nid = $hgene->{father}->{id};
			$hgene->{all_variants}->{$nid}->{father} = 1;
			$hgene->{all_variants}->{$nid}->{score} = $hgene->{father}->{score};
			$hgene->{all_variants}->{$nid}->{added}++;
			push( @{ $class->{father} }, ( $hgene->{father}->{score} - 2 ) );

		}
		elsif (scalar( @{ $class->{father} } ) > 0 && scalar( @{ $class->{mother} } ) == 0 && exists $hgene->{mother}->{id} )
		{
			my $nid = $hgene->{mother}->{id};
			$hgene->{all_variants}->{$nid}->{mother} = 1;
			$hgene->{all_variants}->{$nid}->{score} = $hgene->{mother}->{score};
			push( @{ $class->{mother} }, ( $hgene->{mother}->{score} - 2 ) );
			$hgene->{all_variants}->{$nid}->{added}++;
		}

		my $score_father    = max( @{ $class->{father} } );
		my $score_mother    = max( @{ $class->{mother} } );
		my $score_biallelic = max( @{ $class->{biallelic} } );
		if ( $score_father + $score_mother > $score_biallelic && $patient->getFamily->isTrio()) {
			$hgene->{max_score} = $hgene->{score} + $score_father + $score_mother;
		}
		else {
			
			$hgene->{max_score} = $hgene->{score} + $score_biallelic + $hgene->{penality} ;
		}
	}

#	$project->buffer->close_lmdb();
	
}


sub refine_heterozygote_composite_score_fork {
	my ( $project, $genes,$hchr ) = @_;
	$t = time;
	$| =1;
	print qq{<div style="display: none">};
	print "refine";
	my $fork      = scalar (keys %$hchr);
 	$fork =8;	
	$fork = 16 if $project->isGenome();
	
	$fork=10;
	my $pm        = new Parallel::ForkManager($fork);
	$pm        = new Parallel::ForkManager($fork);
	#@$vgenes = @$vgenes[0..500];
	my $nb        = int( scalar(@$vgenes) / ($fork) +1 );
	my $iter      = natatime( $nb,  @$vgenes );
	my $res_genes = [];
	my $vid        = 0;
	my $hrun;
	
	#$t = time;
	my $current = 1;
	my $vres;
	my @res_final;
	my @toto;
	my $wtime = 0;
	my $maxgene =1000;
	my $ngene =0;
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $h ) = @_;

			unless ( defined($h) or $exit_code > 0 ) {
				print
				  qq|No message received from child process $exit_code $pid!\n|;
				die();
				return;
			}
			my $id = $h->{run_id};
			delete $hrun->{ $h->{run_id} };
			warn "==>" . abs( time - $h->{time} );# if $print;
			$h->{genes} = [] unless $h->{genes};
			$vres->{$id} = $h->{genes};
			while (exists $vres->{$current}){
				print qq{</div>} if $current  == 1 ;
				 my $xtime =time;
					foreach my $g (@{ $vres->{$current}}){
						last if $ngene > $maxgene;
						push(@toto,$g->{name});
						#warn $g->{out};
						print $g->{out};
						print  "\n<!--SPLIT-->\n";
						$ngene ++;
						
					#	push(@res_final,$g->{out})
					}
					delete $vres->{$current};
					$current ++;
					$wtime += abs($xtime - time);
				
			}
			push( @$res_genes, @{ $h->{genes} } );
		}
	);
	
	#TODO: essayer d'integrer XLS Store variant ici pour le forker
	
	$project->buffer->dbh_deconnect();
	
	#delete $project->{validations_query1};
	while ( my @tmp = $iter->() ) {
		$vid++;
		$hrun->{$vid}++;
		my $pid = $pm->start and next;
		my $t   = time;
		my $res;
		$res->{tmp}  = \@tmp;
		$res->{time} = time;
		$project->buffer->dbh_reconnect();
		( $res->{genes}, $res->{total_time} ) = variations_editor_methods::refine_heterozygote_composite( $project,$print_html, \@tmp, $vid);
		$res->{run_id} = $vid;
		$pm->finish( 0, $res );
	}
	$pm->wait_all_children();
	$project->buffer->dbh_reconnect();
	error("Hey it looks like you're having an error  !!! ")
	  if scalar keys %$hrun;
	warn "end hetero " if $print;
	warn "....";
	print qq{</div>};
	warn "....";
	while (exists $vres->{$current}){
				print qq{</div>} if $current  == 1 ;
				 my $xtime =time;
					foreach my $g (@{ $vres->{$current}}){
						last if $ngene > $maxgene;
						print $g->{out} . "\n";
							$ngene ++;
					}
					delete $vres->{$current};
					$current ++;
				$wtime += abs($xtime - time);	
				
			}
		print "<br>". $wtime;
		warn "***** ".$wtime;
	error("Hey it looks like you're having an error  !!! ") if scalar keys %$vres;	
	return ;
	exit(0);
	print qq{</div>};
	return $res_genes;

}




sub fork_annnotations {
	my ( $list, $list_saved, $maskcoding,$vector ) = @_;
	unless ($cgi->param('export_xls')) {
		print qq{<div style="display: none">};
		print "annotations";
	}
	my $fork = 6;
	$fork = 15 if $project->isGenome();
	$fork=10;
	#ici $fork= 20;
	my $nb   = int( (scalar(@$list) +1) / ($fork-1)  );
	$nb = 1 if scalar(@$list) < $fork;
	my $pm   = new Parallel::ForkManager($fork);
	my $hrun = {};
	my $hvariations;
	my $hgenes;

	my $tsum;
	my $tsum_finih = 0;
	my $ttsum;
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $h ) = @_;
			my $t = time;
			unless ( defined($h) or $exit_code > 0 ) {
				print
				  qq|No message received from child process $exit_code $pid!\n|;
				error("Hey it looks like you're having an error  !!! ");
				return;
			}

			my $id = $h->{run_id};
			warn  abs( time - $h->{ttime});
			$ttsum += abs( time - $h->{ttime} );
			delete $h->{run_id};
			delete $hrun->{$id};
			foreach my $gene ( @{ $h->{genes} } ) {
				my $gid = $gene->{id};
				#unless ($keep_pathogenic){
				#die($gene->{name}) if  exists $gene->{DM};
				next
				  unless ( $gene->{mask} & $maskcoding or exists $gene->{DM} );

				#}

				unless ( exists $hgenes->{$gid} ) {

					#my $chr =
					$hgenes->{$gid} = $gene;
					$hgenes->{$gid}->{score_father_mother} = $hgenes->{$gid}->{score_mother} + $hgenes->{$gid}->{score_father};

					$hgenes->{$gid}->{max_score} = $hgenes->{$gid}->{score} + ( $hgenes->{$gid}->{score_father_mother} >  $hgenes->{$gid}->{score_biallelic} ? $hgenes->{$gid}->{score_father_mother} : $hgenes->{$gid}->{score_biallelic} );
					$hgenes->{$gid}->{max_score} -=( $hgenes->{$gid}->{denovo_rare} * 0.3 ) if $hgenes->{$gid}->{denovo_rare} > 2;
					$hgenes->{$gid}->{max_score} += 2 if $hgenes->{$gid}->{cnv_dup} > 1 or $hgenes->{$gid}->{cnv_del} > 1;
					foreach my $k ( keys %{ $gene->{variant} } ) {
						$hgenes->{$gid}->{all_variants}->{$k} = $gene->{variant}->{$k};
						$hgenes->{$gid}->{all_vector_ids}->{$k} = $gene->{vector_ids}->{$k};
					}
					next;
				}
				push( @{ $hgenes->{$gid}->{variants} }, @{ $gene->{variants} } );

				foreach my $k ( keys %{ $gene->{variant} } ) {
					$hgenes->{$gid}->{all_variants}->{$k} = $gene->{variant}->{$k};
					$hgenes->{$gid}->{all_vector_ids}->{$k} = $gene->{vector_ids}->{$k};
				}

			}
			
			$tsum += ( abs( time - $t ) );

		}
	);
	my $iter = natatime( $nb, @$list );
	my $id   = time;
	$project->getChromosomes();
	$project->buffer->dbh_deconnect();
	#$project->buffer->close_lmdb();
	my $tt = time;
	my $final_polyviewer_all ;
	if ($project->isRocks){
	   $final_polyviewer_all = GenBoNoSqlRocks->new(dir=>$project->rocks_directory."/patients/",mode=>"r",name=>$patient->name);
	}
	while ( my @tmp = $iter->() ) {
		$id++;
		$hrun->{$id}++;
		my $pid = $pm->start and next;
		my $t   = time;
		my $res = annotations_methods::annotations( $project,$patient, \@tmp, $list_saved, $maskcoding,$final_polyviewer_all,$hash_genes_panel,$hash_variants_DM );
		$res->{run_id} = $id;
		$res->{ttime}  = time;
		$pm->finish( 0, $res );

	}
	$pm->wait_all_children();
	print qq{</div>} unless ($cgi->param('export_xls'));
	warn '----- after Annotation: '
	  . $tsum . ' final :'
	  . abs( $tt - time ) . ' :: '
	  . $ttsum
	  if $print;
	warn '----- nb genes ' . scalar( values %$hgenes ) if $print;
	error("Hey it looks like you're having an error  !!! ")
	  if scalar keys %$hrun;
	die("-------- PROBLEM ----------- ") if scalar keys %$hrun;
	#refine_heterozygote_composite_score($hgenes);
	#my @genes = sort {$b->{max_score} <=> $a->{max_score}} values %$hgenes;
	$project->buffer->dbh_reconnect();
	#warn 'end annotations';
	return ( [ values %$hgenes ], $hvariations );
}

sub return_litedb {
	my ($project) = @_;
	return GenBoNoSql->new( dir => $project->lmdb_cache_dir(), mode => 'r' );
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


sub constructChromosomeVectorsPolyDiagFork {
	my ( $project, $patient,$startVector,$statistics ) = @_;
	unless ($cgi->param('export_xls')) {
		print qq{<div style="display: none">};
		print "vector";
	}
	my $filter_transmission;
	my $xtime =time;
	my $list_genes;
	$filter_transmission->{denovo}          = 1 if $cgi->param('denovo');	
	#$filter_transmission->{"strict_denovo"} = 1 if $cgi->param('denovo');
	$filter_transmission->{"strict_denovo"} = 1 if $cgi->param('strict_denovo');
	$filter_transmission->{'denovo/?'}      = 1 if $cgi->param('denovo');
	
	$filter_transmission->{recessive}       = 1 if $cgi->param('recessive');
	$filter_transmission->{both}            = 1 if $cgi->param('both');
	$filter_transmission->{xor_mother}      = 1 if $cgi->param('xor_mother') or $cgi->param('xor');
	$filter_transmission->{xor_father}      = 1 if $cgi->param('xor_father') or $cgi->param('xor');
	#$filter_transmission->{xor}      = 1 if  $cgi->param('xor');
	if (exists $filter_transmission->{xor_mother} && exists $filter_transmission->{xor_father} ){
		$filter_transmission->{xor} =1;
		delete $filter_transmission->{xor_mother};
		delete $filter_transmission->{xor_father};
	}
	if ( exists $filter_transmission->{denovo} && exists $filter_transmission->{strict_denovo}){
		delete $filter_transmission->{strict_denovo};
	}
	
	$filter_transmission->{xor}             = 1 if  $cgi->param('xor');
	warn "\n";
	
	my $hashVector_panel = {};
	my $hashVector       = {};
	construct_panel_vector( $panel, $hashVector_panel ) if $panel;
	my $list_transcript;
	my $trio = $patient->getFamily->isTrio;
	my $list_variants =[];
	my $hash_variants_DM = {};
	my $gene;

	if ($gene_name_filtering) {
		$gene = $project->newGene($gene_name_filtering);
		#error("GENE NAME $gene_name_filtering NOT FOUND !!!!!!") unless $gene;
		$gene_id_filtering = $gene->id();
		
	}
#	my $gene_exception = $project->newGene("PCDH19");
	my $fork = 24;
	#my $pm   = new Parallel::ForkManager($fork);
	my $hrun;
#	$pm->run_on_finish(
#		sub {
#			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $h ) = @_;
#			unless ( defined($h) or $exit_code > 0 ) {
#				warn
#				  qq|No message received from child process $exit_code $pid!\n|;
#				die();
#				return;
#			}
#		
#			
#			my $chr = $h->{chromosome};
#			foreach my $k  (keys %{$h->{statistics}}){
#				$statistics->{$k} += $h->{statistics}->{$k};
#				
#			}
#			unless ( exists $hashVector->{$chr} ) {
#				$hashVector->{$chr} = $h->{vector};
#			}
#			else {
#				$hashVector->{$chr} &= $h->{vector};
#			}
#			my $id = $h->{run_id};
#			delete $hrun->{ $h->{run_id} };
#			error("PROBLEM !!!!") unless exists  $h->{run_id};
#			push( @$list_variants, @{ $h->{list_variants} } ) if  $h->{list_variants};
#			map { $hash_variants_DM->{$_}++ } @{ $h->{list_variants_DM} };
#
#		}
#		
#	);
	my $id = time;
	$project->disconnect();
	delete $project->{rocks};
	my $hno;
	foreach my $chr ( @{ $project->getChromosomes } ) {
		#$hno->{$chr->name} = GenBoNoSqlRocksVector->new(chromosome=>$chr->name,dir=>$project->rocks_directory("vector"),mode=>"r",name=>$chr->name);
	}
	my $finalVector = {};
#	warn $project->rocks_directory("vector");
#	warn "------";
	foreach my $chr ( @{ $project->getChromosomes } ) {

		if ($gene) {
			next if ( $gene->getChromosome()->name ne $chr->name );
		}
		if ($panel) {
			next unless exists $hashVector_panel->{ $chr->name };
		}
		$id++;
		$hrun->{$id} ++;
		#my $pid = $pm->start and next;
		my $xxt = time;
		#$project->buffer->dbh_reconnect();
		#$chr->rocks_vector->rocks;
		
		#my $no = GenBoNoSqlRocksVector->new(chromosome=>$chr->name,dir=>"/tmp/vector",mode=>"r",name=>$chr->name); #$chr->flush_rocks_vector();
		my $no = $chr->flush_rocks_vector();
		$no->prepare_vector([$limit_ac,$limit_ac_ho,$limit_sample_dv,$limit_sample_dv_ho,"intergenic","dm",$patient->name]);
		my $res = {};
		my $debug;
		my $statistics = {};
		
		my $hashVector = {};
		
		if ($gene) {

			$hashVector->{ $chr->name } = $gene->getVectorOrigin();

		}
		elsif ($panel) {
			$hashVector->{ $chr->name } =
			  $hashVector_panel->{ $chr->name }->Clone;

		}
		else {
			$hashVector->{ $chr->name } = $startVector->{$chr->name};
		}
		
		
		print "=" unless $cgi->param('export_xls');
			my $testid = 4201;
		my $debug;
		
		#	$debug =1 if $chr->name eq  "12";
		if ($panel) {
			#next unless  (exists $hashVector->{$chr->name});
			$hashVector->{ $chr->name } &= $no->get_vector_chromosome($limit_ac);    #if $limit_ac ne "all";
			
			#$hashVector->{ $chr->name } &=
			 # $chr->getVectorScore($limit_ac_ho);    # if $limit_ac_ho ne "all";
			$hashVector->{ $chr->name } &= $no->get_vector_chromosome($limit_sample_dv) ;    # if $limit_sample_dv ne "all";
			$hashVector->{$chr->name} &= $no->get_vector_chromosome($limit_sample_dv_ho);
			$hashVector->{ $chr->name } -= $no->get_vector_chromosome("intergenic");

		}

		else {
			#$hashVector->{ $chr->name } &= $chr->getVariantsVector();
			$hashVector->{ $chr->name } = $no->get_vector_chromosome($limit_ac);    # if $limit_ac ne "gnomad_ac_all";;
			$hashVector->{ $chr->name } &= $no->get_vector_chromosome($limit_ac_ho);    # if $limit_ac_ho ne "all";
			$hashVector->{ $chr->name } &= $no->get_vector_chromosome($limit_sample_dv);    #  if $limit_sample_dv ne "all";
			 $hashVector->{$chr->name} &= $no->get_vector_chromosome($limit_sample_dv_ho);
			#$hashVector->{$chr->name} &= $chr->getVectorLargeDeletions;
			#$hashVector->{ $chr->name } &= $vquality if $vquality;
			$hashVector->{ $chr->name } -= $no->get_vector_chromosome("intergenic");
		
			
		}
		my $vannotations =  $chr->getNewVector ;
		$no->prepare(\@tconsequences);
		foreach my $annot  (@tconsequences){
			 $vannotations += $no->get_vector_chromosome("$annot");
		}
		$hashVector->{ $chr->name } &=$vannotations;
		
		
		#$hashVector->{ $chr->name } &= $no->get_vector_chromosome($chr);
		$statistics->{variations} += $patient->countThisVariants( $hashVector->{ $chr->name } );

		my $vDM = $chr->vectorDM();
		$vDM += $chr->vectorClinvarPathogenic();
		#TODO: ajout susceptibilité ici;
		$vDM += $chr->vectorVariantsForcedViewing();
		$vDM &= $patient->getVectorOrigin($chr);
		#$vDM &= $hashVector->{ $chr->name };
		$statistics->{DM} += $patient->countThisVariants($vDM);
		my $vtr               = $chr->getNewVector();
		#$hashVector->{$chr->name}  |= $v1;
		if ($trio) {
			my @list_transmission = keys %$filter_transmission;
			
			
			if ( $patient->getFamily()->isDominant() ) {
				my $none = undef;
				if (exists $filter_transmission->{xor_mother}){			
							$vtr |= $patient->getFamily()->getVector_individual_mother( $chr, $patient );
							$none =1;
				}
				if (exists $filter_transmission->{xor_father}){		
							$vtr |= $patient->getFamily()->getVector_individual_father( $chr, $patient );
							$none =1;
				}
				if (exists $filter_transmission->{both}){
				  $vtr = $patient->getVectorOrigin($chr);
				  $none =1;
				}
				unless ($none){
					$vtr |= $patient->getFamily()->getVector_family_dominant($chr);
				}
				
	
				

			}
			else {
				foreach my $tr ( keys %$filter_transmission ) {
						print "__";
					if ( $tr eq "both" ) {
						 $vtr = $patient->getVectorOrigin($chr);
						last;
					}
					if ( $tr eq "denovo") {
							$vtr |= $patient->getFamily()->getVectorDenovoTransmission( $chr, $patient );
					}
					
					elsif ( $tr eq "strict_denovo" ) {
							$vtr |= $patient->getFamily()->getVector_individual_strict_denovo( $chr, $patient );

					}
					elsif ( $tr eq "recessive" ) {
						$vtr |= $patient->getFamily()->getVectorRecessiveTransmission( $chr, $patient );
					}
					elsif ( $tr eq "xor" ) {
						#my $vtr2 = $vtr->Clone();
						my $vtr2 = $chr->getNewVector();
						$vtr2 |= $patient->getFamily()->getVectorMotherTransmission( $chr, $patient ) if $patient->getFamily->getMother();
						$vtr2 |= $patient->getFamily()->getVectorFatherTransmission( $chr, $patient )if $patient->getFamily->getFather();
						$vtr2 -= $patient->getFamily->getFather->getVectorHo( $chr,$patient ) if $patient->getFamily->getFather();
						$vtr2 -= $patient->getFamily->getMother->getVectorHo( $chr, $patient ) if $patient->getFamily->getMother();
						if ($chr->name eq "X" && $patient->isMale ) {
							$vtr2 -=  $patient->getFamily()->getVectorRecessiveTransmission( $chr, $patient ) if $patient->getFamily->getMother;
						}
						$vtr |= $vtr2;
					}
					elsif ( $tr eq "xor_mother" ) {

						$vtr |= $patient->getFamily()->getVector_individual_mother( $chr, $patient );
					}
					elsif ( $tr eq "xor_father" ) {

						$vtr |= $patient->getFamily()->getVector_individual_father( $chr, $patient );
					}
				}
				
#				if ($chr->name eq $gene_exception->getChromosome->name and $patient->getVectorOrigin($chr)){
#					my $vector = $gene_exception->getVectorOrigin();
#				 	$vector &= $patient->getVectorOrigin($chr);
#					my $father = $patient->getFamily->getFather();
#					my $mother = $patient->getFamily->getMother();
#					if ($father){
#						my $v1 = $father->getVectorOrigin($chr);
#						$vector &= $v1;
#						if ($mother){
#							my $v2  = $mother->getVectorOrigin($chr);
#							$vector -= $v2;
#						}
#						 $vtr |= $vector;
#					}
#					
#				}
			
				
			}
			$hashVector->{ $chr->name } &= $vtr;
			$hashVector->{ $chr->name } &= $patient->getVectorOrigin($chr);

		}
		
		#  die($chr->name) if $debug && $hashVector->{ $chr->name }->contains(7710) && $debug;
		$hashVector->{ $chr->name } &= $hashVector_panel->{ $chr->name } if exists $hashVector_panel->{ $chr->name };
		  $statistics->{variations} += $patient->countThisVariants( $hashVector->{ $chr->name } );
		if($only_DM){
				#keep only pathogenic variation if cgi with option $only_DM
				$hashVector->{ $chr->name } = $vDM;
		}
		if ($keep_pathogenic) {
			$hashVector->{ $chr->name } += $vDM;
		}
		$res->{vector}        = $hashVector->{ $chr->name };

		###############
		# AFFINE RATIO 
		##############
		if($limit_ratio>0){
			print "XX";
			my $vquality  =  $chr->getNewVector();
			if ($limit_ratio1 < 100){
				my $vector_ratio_name = $patient->name . "_ratio_" . $limit_ratio1;
				 $vquality = $no->get($vector_ratio_name);
			}
			
			if ($limit_ratio2 != $limit_ratio1) {
				die();
				my $no = $chr->lmdb_polyviewer_variants( $patient, "r" );
			#	my $no2 = $chr->lmdb_polyviewer_variants_genes( $patient, "r" );
				my $vector_ratio_name = $patient->name . "_ratio_" . $limit_ratio2;
				$vector_ratio_name = $patient->name . "_ratio_all" if $limit_ratio2 == -1;
				 my $vquality2 = $no->get($vector_ratio_name);
				 $vquality2 -= $vquality;
				
				$vquality2 &=  $res->{vector};
				$res->{vector} &= $vquality;	
				foreach my $id ( @{ to_array( $vquality2, $chr->name ) } ) {
					print "!";
					my $av = $no->get($id);
					my ($a1,$a2)= split(":",$no->get($av->{id})->{value}->{ratio}->[0]);
					$a2 =~ s/%//;
					next if $a2 <  $limit_ratio;
					my ($i1,$i2) = split("!",$id);
					$res->{vector}->Bit_On($i2);	
					
		 			}
		 		
			}
			else {
				$res->{vector} &= $vquality;	
			}
		}
	
		###############
		# in this run 
		##############
		if ($in_this_run < 100 ){
			print " in_this_run $in_this_run ";
			my $no = $chr->lmdb_polyviewer_variants( $patient, "r" );
			foreach my $id ( @{ to_array($res->{vector}, $chr->name ) } ) {
				my $av =  $no->get($id);
				
				my $vobj   = $project->returnVariants($id);#$chr->cache_lmdb_variations->get_varid($id);
				
				if ($vobj->in_this_run_ratio * 100 > $in_this_run) {
					my ($i1,$i2) = split("!",$id);
					$res->{vector}->Bit_Off($i2);	
				}
				
			}
		}
		

#		my $vLocal = $chr->getVectorLocalValidations();
#		
#
#		$vLocal &= $patient->getVectorOrigin($chr);
#
#		#Limitation gnomad ho ac a 1000 si projet exome ou genome
#		if ( $project->isExome() || $project->isGenome() ) {
#			$vLocal &= $chr->getVectorScore("gnomad_ho_ac_1000");
#		}
#
#		$res->{vector} |= $vLocal;
		
		###############
		#¥¥¥ END Chromosome
		##############
		my $xxx = 0;
		#push(@$list_variants,join($chr->name."!",split(",",$res->{vector}->to_Enum)));
	
		to_array($res->{vector}, $chr->name,$list_variants);
		if ($keep_pathogenic){
			to_hash($vDM, $chr->name,$hash_variants_DM);
		}
		push(@$list_genes,@{$chr->getGenesIdFromVector($res->{vector})});

		
				
			$finalVector->{$chr->name} = $res->{vector} ;
#			warn "end chr".$chr->name;
			delete $project->{rocks};
		}
	#$pm->wait_all_children();
	warn "END !!!!!!!!! ".abs(time -$xtime);
	#error("ARGGG Problem") if keys %$hrun;
	print qq{</div>} unless $cgi->param('export_xls');
	return ( $finalVector, $list_variants, $hash_variants_DM,$list_genes);
}



sub to_hash {
	my ( $v, $name,$list) = @_;
	my @pos1 = split(",", $v->to_Enum);
	$list ={} unless $list;
	foreach my $pos (@pos1){
		my ($a,$b) = split("-",$pos);
		unless ($b) {
			 $list->{$name."!".$a} ++
		}
		else {
			for my $i ($a..$b) {
				 $list->{$name."!".$i} ++
			}	
		}
	}
	return $list;

}


sub to_array {
	my ( $v, $name,$list) = @_;
	my @pos1 = split(",", $v->to_Enum);
	$list =[] unless $list;
	foreach my $pos (@pos1){
		my ($a,$b) = split("-",$pos);
		unless ($b) {
			push(@$list,$name."!".$a);
		}
		else {
			for my $i ($a..$b) {
				push(@$list,$name."!".$i);
			}	
		}
	}
	return $list;
#	my $set  = Set::IntSpan::Fast::XS->new( $v->to_Enum );
#	my $iter = $set->iterate_runs();
#	my @t;
#	while ( my ( $from, $to ) = $iter->() ) {
#		for my $member ( $from .. $to ) {
#			if ($name) {
#				push( @t, $name . "!" . $member );
#			}
#			else {
#				push( @t, $member );
#			}
#		}
#	}
#	return \@t;
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



foreach my $chr  (@{$project->getChromosomes}){
		my $no = $chr->lmdb_polyviewer_variants( $patient, "r" );
		my @st = (stat($no->filename));
		 push(@keys, ($st[9].$st[11].$st[12]));
		my $no2 = $chr->lmdb_polyviewer_variants_genes( $patient, "r" );
		@st = (stat($no2->filename));
		 push(@keys, ($st[9].$st[11].$st[12]));
}
#push(@keys,file_md5_hex($Bin."/variations_editor.pl") );
push(@keys,$VERSION);
push(@keys,$VERSION_UPDATE_VARIANT_EDITOR );
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



my $date;
( $date->{cache} ) = "toto";#utility::return_date_from_file( $cno->filename );
( $date->{bam} )   = "tutu";#utility::return_date_from_file( $patient->getBamFileName );
return $date;
}


sub print_hotspot {
	my ($patient) = @_;
	my $hotspots = $patient->hotspot;
	return "" unless $hotspots;

	my $out ="";
	
	#$out .=  $cgi->start_div({class=>"panel-heading panel-alert alert ",style=>" min-height:13px;max-height:13px;padding:1px;border:1px"});
	my $label_id = "hs_polyviewer_".$patient->id;
	my $panel_id = "pa_polyviewer_".$patient->id;
	$out .= qq{<div class="btn  btn-warning btn-xs " style="position:relative;bottom:1px;min-width:150px;" onClick='collapse("$panel_id","$label_id")'>  <span id= "$label_id" class="glyphicon glyphicon-triangle-right  "   style="float:left;"></span> HOTSPOT &nbsp</div>};
	#$out .=$cgi->span({class=>"label label-success"},qq{<span class='badge badge-primary badge-xs'  >-</span>});
	$out .=  $cgi->start_div({class=>"panel-body panel-collapse  collapse",style=>"width:50%;font-size: 09px;font-family:  Verdana;",id=>$panel_id});
	my $div_alert;	
	my $s_id = $patient->{name};

my $t = time;


#$out .=  $cgi->start_div({class=>"panel-heading panel-warning warning ",style=>" min-height:13px;max-height:13px;padding:1px;border:1px"});
#	$out .= qq{<div class="btn  btn-success btn-xs " style="position:relative;bottom:1px;min-width:150px;" onClick='collapse("$panel_id","$label_id")'>  <span id= "$label_id" class="glyphicon glyphicon-triangle-right  "   style="float:left;"></span> $text &nbsp</div>};
	   		#	$out .=$cgi->span({class=>"label label-success"},qq{<span class='badge badge-primary badge-xs'  >$nbv</span>});
#		my $nbv = scalar (keys %{$hotspot->{results}->{$s_id}});
#		$out .=$cgi->span({class=>"label label-success"},qq{<span class='badge badge-primary badge-xs'  >$nbv</span>});	
	my @header = ("ID","NAME","PROT","A","C","G","T","DEL","COV");	 
	$out .= $cgi->start_table({class=>"table table-striped table-condensed table-bordered table-hover table-mybordered",style=>"font-size: 9px;font-family:  Verdana;"});
foreach my $g (keys %$hotspots){
	#$out .=  $cgi->start_div({class=>"panel panel-info" });
	 #panel heading
	 
	 # $out.= $cgi->end_div();
		#REF	POS	COV	A	C	G	T	DEL	REFSKIP	SAMPLE
	#my $var_obj = $self->cache_lmdb_variations->get($vid);
	#  panel table
	$out.= $cgi->start_Tr();
	$out.=$cgi->th({colspan=>(scalar(@header)+1),style=>"background-color:#217DBB;color:white;font-size:12px"},$g);
	$out.= $cgi->end_Tr();
	$out.= $cgi->start_Tr();
	$out.=$cgi->th({style=>"background:#E0E0FF"},["igv",@header]);
	$out.= $cgi->end_Tr();
	
	my @bams;
	my @names;
	foreach my $p (@{$patient->getFamily->getPatients()}){
		push(@bams,$patient->bamUrl);
		push(@names,$patient->name());
	}
					
	my $f =  join(";",@bams);#$patient->{obj}->bamUrl;;
	 my $pnames = join(";",@names);
	foreach my $hotspot (@{$hotspots->{$g}}){
		my @td;
		my $chr = $project->getChromosome($hotspot->{REF});
		#my $var_obj = $chr->cache_lmdb_variations->get($hotspot->{GENBO_ID});
		my $var_obj = $project->_newVariant($hotspot->{GENBO_ID});
		#warn $chr->cache_lmdb_variations->get(0);
		my $style ={};
		 $style = {style=>"background-color:#D2386C;color:white"} if $var_obj && defined $var_obj->vector_id() && $var_obj->existsPatient($patient);
		 $out.= $cgi->start_Tr($style);
		 my $nba;
		 my $chrn = $chr->name;
		 my $start = $hotspot->{POS};
		 my $l = $chr->name.":".$start;
		 my $gn = $project->getVersion();
		 my $project_name = $project->name; 
		 my $v1 = "?/?";#.$hvariation->{allele};	
		# launch_web_igv_js
		my $text =qq{<button class='igvIcon2' onclick='launch_web_igv_js("$project_name","$pnames","$f","$l","$v1","$gn")' style="color:black"></button>};
		#my $text =qq{<button dojoType="dijit.form.Button"   iconClass='igvIcon' onclick='view_web_igv_bam("dialog_igv", "div_igv", "$l", "$f", "$pnames")' style="color:black"></button>};
		 
		 
		$out.=$cgi->td($text);
		foreach my $h (@header){
			if ($h eq  $hotspot->{A_ALT}){
				my $pc = int(($hotspot->{$h}/$hotspot->{COV})*1000)/10;
				my $color = "#f2dedc";
				$color = "#F7BFB9" if $pc>2;
				$color = "#E9897E" if $pc>5;
				$out.=$cgi->td({style=>"background-color:$color;color:black"},"$pc% (".$hotspot->{$h}.")");
			}
			elsif ($h eq  $hotspot->{A_REF}){
				my $pc = int(($hotspot->{$h}/$hotspot->{COV})*1000)/10;
				$out.=$cgi->td({style=>"background-color:#c7eadd;color:black"},"$pc% (".$hotspot->{$h}.")");
			}
			elsif ($h eq  "DEL" && $hotspot->{A_ALT} eq "-"){
				my $pc = int(($hotspot->{$h}/$hotspot->{COV})*1000)/10;
				my $color = "#f2dedc";
				$color = "#F7BFB9" if $pc>2;
				$color = "#E9897E" if $pc>5;
				$out.=$cgi->td({style=>"background-color:#F7BFB9;color:black"},"$pc% (".$hotspot->{$h}.")");
			}
			else {
				my $pc = int(($hotspot->{$h}/$hotspot->{COV})*1000)/10;
				$out .= $cgi->td($hotspot->{$h});
			}
			#push(@td, $hotspot->{$h});
		}
	
		#$out.=$cgi->td(\@td);	
		$out.= $cgi->end_Tr();
	}
	
	}
	$out.= $cgi->end_table();	
	$out.= $cgi->end_div();	#$out.="<!-- 3 -->";	
	return $out;
}
