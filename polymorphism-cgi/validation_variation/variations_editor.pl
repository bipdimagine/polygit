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
use Compress::Snappy;
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

#use PDF::API2;
#use PDF::Table;
use constant mm => 25.4 / 72;
use constant in => 1 / 72;
use constant pt => 1;
use Time::HiRes qw ( time alarm sleep );
use Storable qw/thaw freeze/;
use Digest::MD5 qw(md5 md5_hex md5_base64);
use LMDB_File qw(:flags :cursor_op :error);
use Proc::Simple;
use lib "$Bin/variation_editor_lib";
use hvariant;
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
my $myproc = Proc::Simple->new();        # Create a new process object
 
$myproc->start("$Bin/vmtouch.pl project=".$cgi->param('project')." patient=".$cgi->param('patients'));
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
	"5" => "sdv_500",
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

my $bgcolor = qq{background-color:#607D8B}
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

my $version_db =  $project->public_database_version;
if($version_db>=13){
	$VERSION = $VERSION."-".$version_db."1";
}
my $patient = $project->getPatient($patient_name);
my $fam = $patient->getFamily();
my $print_html = polyviewer_html->new( project => $project, patient => $patient );
$print_html->init();

if ($patient->isMale){
	$VERSION .= "11" ;#if $patient->isMale;
}

unless ($cgi->param('phenotype')){
	if ($project->phenotypes){
		$cgi->param(-name=>'phenotype',-value=>$project->phenotypes->[0]);
	}
}
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
if ($force == 1) {
 $dev = 2;	
}
$no_cache = $patient->get_lmdb_cache("w");

my $keys = return_uniq_keys($patient,$cgi);
my $level_dude = 'high,medium';


my $cache_dude_id =
  "$level_dude::" . $project_name . "-" . $patient->name . "-" . $VERSION;
my $cache_id = join( ";", @$keys );
my $text = $no_cache->get_cache($cache_id);
$text = "" if $dev;

my $html_dude = "<!--DUDE-->";
my $cache_icon;

 $cache_icon = qq{<span class="glyphicon glyphicon-floppy-remove" aria-hidden="true" style="text-align:right;font-size:10px;color:red"></span>};
if($text){
		$text =~ s/onClick='editor\(1,2\);'/onClick='load_polyviewer_export_xls\(1\);'/g;
	 	$cache_icon = qq{<span class="glyphicon glyphicon-floppy-saved" aria-hidden="true" style="text-align:right;font-size:10px;color:green"></span>};
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
my $javascript_id = int( time + rand(400000) );
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

foreach my $a (@annots) {
	foreach my $this_a ( split( ',', $buffer->get_genbo_annotation_name($a) ) )
	{
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
$version->{gencode} = $project->get_gencode_description();
$version->{gnomad}  = $project->get_public_data_description("gnomad-exome");
$version->{hgmd}    = $project->get_public_data_description("hgmd");
$version->{cadd}    = $project->get_public_data_description("cadd");
$version->{clinvar} = $project->get_public_data_description("clinvar");
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
$panel = $project->getPanel($panel_name)
  if $panel_name && lc($panel_name) ne "all";
my $hash_genes_panel;
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
	}
	$project->{phenotypes_object} = $hPhenoIds;
}


$t = time;
$buffer->disconnect();

#my $list_saved;




##################################
################## GET VECTORS 
##################################

	my ( $vectors, $list, $hash_variants_DM ) = constructChromosomeVectorsPolyDiagFork( $project, $patient,$statistics );
	$ztime .= ' vectors:' . ( abs( time - $t ) );
	warn $ztime if $print;



$t = time;
##################################
################## GET GENES 
##################################
my ($genes) = fork_annnotations( $list, [], $maskcoding );

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

calculate_max_score($project,$genes);
warn "end";
	my $vgenes =[];
	my $nb2 = 0;
	my $limit =10000; 
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
	
	
$myproc->kill();

#my $exit_status = $myproc->wait(); 
#warn "ok *** $exit_status ";
#exit(0);
my $stdout_end = tee_stdout {
warn "genes:".scalar(@$genes);
if (@$genes){
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

warn "write cache";

$no_cache->put_cache_text($cache_id,$stdout.$stdout_nav_bar.$stdoutcnv.$stdout_end,2400) ;#unless $dev;

$no_cache->close();
 warn "END";
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
	die() if $TEST;
	exit(0);
}

sub calculate_max_score {
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
	foreach my $hgene (@$list) {
		my $gid = $hgene->{id};
		my $debug ;
		my $gene = $project->newGenes( [$gid] )->[0];
		my ( $n, $chr_name ) = split( "_", $gid );
		$hgene->{chr_name} = 
		my $chr = $gene->getChromosome();
		$hgene->{chr_name} = $chr_name;
		if ( $no_dude && -e $no_dude->filename ) {
			$hgene->{level_dude} = $no_dude->get($gid);
		}
		else {
			$hgene->{level_dude} = -1;
		}

		my $global_gene = $chr->lmdb_polyviewer_genes($patient)->get($gid);
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
			warn "mother" if $debug;
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
	#$fork=20;
	my $pm        = new Parallel::ForkManager($fork);
#	my $t = time;
#	
#	foreach my $chr_name (keys %$hchr){
#		my $pid = $pm->start and next;
#		print "!";
#		my $no       = $project->getChromosome($chr_name)->lmdb_polyviewer_variants( $patient, "r" );
#		$no->test(1);
#		$no->lmdb($no->name);
#		print "!";
#		$pm->finish();
#		#$no->lmdb($chr->name);
#	}
#	$pm->wait_all_children();
	#warn "cp ".abs(time-$t);
	$fork = 10;
	 $pm        = new Parallel::ForkManager($fork);
	my $nb        = int( scalar(@$vgenes) / ($fork*2) + 1 );
	
	my $iter      = natatime( $nb,  @$vgenes );
	my $res_genes = [];
	my $vid        = 0;
	my $hrun;
	
	#$t = time;
	my $current = 1;
	my $vres;
	my @res_final;
	my @toto;
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
			warn "==>" . abs( time - $h->{time} ) if $print;
			$h->{genes} = [] unless $h->{genes};
			$vres->{$id} = $h->{genes};
			while (exists $vres->{$current}){
				print qq{</div>} if $current  == 1 ;
					foreach my $g (@{ $vres->{$current}}){
						push(@toto,$g->{name});
						print $g->{out} . "\n<!--SPLIT-->\n";
					#	push(@res_final,$g->{out})
					}
					delete $vres->{$current};
					$current ++;
					
				
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

		( $res->{genes}, $res->{total_time} ) =
		  new_refine_heterozygote_composite_score_old( $project, \@tmp, $vid );

		$res->{run_id} = $vid;
		$pm->finish( 0, $res );
	}
	$pm->wait_all_children();
	$project->buffer->dbh_reconnect();
	error("Hey it looks like you're having an error  !!! ")
	  if scalar keys %$hrun;
	warn "end hetero " if $print;
	print qq{</div>};
	while (exists $vres->{$current}){
				print qq{</div>} if $current  == 1 ;
					foreach my $g (@{ $vres->{$current}}){
						print $g->{out} . "\n";
					}
					delete $vres->{$current};
					$current ++;
					
				
			}
	error("Hey it looks like you're having an error  !!! ") if scalar keys %$vres;	
	#print scalar (keys %{$vres});
	#warn  scalar (keys %{$vres});
	#return;
	return ;
	exit(0);
	print qq{</div>};
	return $res_genes;

}




sub fork_annnotations {
	my ( $list, $list_saved, $maskcoding ) = @_;
	unless ($cgi->param('export_xls')) {
		print qq{<div style="display: none">};
		print "annotations";
	}
	my $fork = 6;
	$fork = 10 if $project->isGenome();
	#ici $fork= 20;
	my $nb   = int( scalar(@$list) / ($fork) + 0.5 );
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
			$ttsum += abs( time - $h->{ttime} );
			delete $h->{run_id};
			delete $hrun->{$id};
			foreach my $gene ( @{ $h->{genes} } ) {

				my $gid = $gene->{id};
				my $debug;
				$debug = 1 if $gene->{name} eq "SDHA";

				#unless ($keep_pathogenic){
				#die($gene->{name}) if  exists $gene->{DM};
				next
				  unless ( $gene->{mask} & $maskcoding or exists $gene->{DM} );

				#}
				#warn Dumper  $gene if $gene->{name} eq "SDHA";

				unless ( exists $hgenes->{$gid} ) {

					#my $chr =
					$hgenes->{$gid} = $gene;
					$hgenes->{$gid}->{score_father_mother} = $hgenes->{$gid}->{score_mother} + $hgenes->{$gid}->{score_father};

					$hgenes->{$gid}->{max_score} = $hgenes->{$gid}->{score} + ( $hgenes->{$gid}->{score_father_mother} >  $hgenes->{$gid}->{score_biallelic} ? $hgenes->{$gid}->{score_father_mother} : $hgenes->{$gid}->{score_biallelic} );
					$hgenes->{$gid}->{max_score} -=( $hgenes->{$gid}->{denovo_rare} * 0.3 ) if $hgenes->{$gid}->{denovo_rare} > 2;
					$hgenes->{$gid}->{max_score} += 2 if $hgenes->{$gid}->{cnv_dup} > 1 or $hgenes->{$gid}->{cnv_del} > 1;
					foreach my $k ( keys %{ $gene->{variant} } ) {
						$hgenes->{$gid}->{all_variants}->{$k} =
						  $gene->{variant}->{$k};
					}
					next;
				}
				push( @{ $hgenes->{$gid}->{variants} },
					@{ $gene->{variants} } );

				foreach my $k ( keys %{ $gene->{variant} } ) {
					$hgenes->{$gid}->{all_variants}->{$k} =
					  $gene->{variant}->{$k};
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

	while ( my @tmp = $iter->() ) {
		$id++;
		$hrun->{$id}++;

		my $pid = $pm->start and next;
		my $t   = time;
		my $res = annotations2( $project, \@tmp, $list_saved, $maskcoding );

		$res->{run_id} = $id;
		$res->{ttime}  = time;
		$pm->finish( 0, $res );

	}
	$pm->wait_all_children();
	print qq{</div>} unless ($cgi->param('export_xls'));
	warn '----- after Annotation: '
	  . $tsum . ' '
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

sub annotations2 {
	my ( $project, $list, $list_saved, $maskcoding ) = @_;
	my $cgi = new CGI();

	#$project->buffer->dbh_reconnect();
	print "." unless ($cgi->param('export_xls'));
	my $tsum = 0;

	#my $vs =  $project->myflushobjects($list,"variants");;
	my $tglobal = time;
	my $res     = {};
	my $agenes  = [];
	my $e;
	my $nb = 1;
	my $mtime_flush;
	my $mtime_gene;
	my $mtime;

	#my $lite = return_litedb($project);
	#foreach my $v (@$vs) {
	$javascript_id += int( rand(10000) );

	#	$project->setListVariants($list);
	#	my $t = time;
	my $list2 = [];
	my $dchr;
	foreach my $id (@$list) {
		my $debug;
		my ( $cname, $xid ) = split( "!", $id );
		$nb++;
		unless ($cgi->param('export_xls')) {
			print "." if $nb % 100 == 0;
		}
		my $t = time;

#my $chr = $project->getChromosome($cname);
#CHANG LMDB POLYVIEWER : $chr->{$cname} = $project->getChromosome($cname)->lmdb_hash_variants("r") unless exists $dchr->{$cname};

		$dchr->{$cname} =
		  $project->getChromosome($cname)->lmdb_polyviewer_variants_genes( $patient, "r" )
		  unless exists $dchr->{$cname};

		$tsum += abs( $t - time() );

		my $hg = $dchr->{$cname}->get($id);
		unless ($hg){
			$hg = hvariant::hash_variant($patient,$id);
		}
		
		#TODO: faire EXPORT XLS ici pour simplifier MAJ

		if ($hg) {
			if ($hash_genes_panel) {
				foreach my $ag ( @{ $hg->{array} } ) {
					next unless exists $hash_genes_panel->{ $ag->{id} };

					#die();
					push( @$agenes, $ag );    # if $ag->{mask} & $maskcoding;
				}
			}
			elsif ( exists $hash_variants_DM->{$id} ) {
				foreach my $a ( @{ $hg->{array} } ) {
					$a->{DM}++;
					push( @$agenes, $a );
				}

				#push(@$agenes, map{$_->{DM}++} @{$hg->{array}});
				#warn Dumper @$agenes;
			}
			else {
				#push(@$agenes, @{$hg->{array}});
				#warn Dumper $hg->{array};
				push( @$agenes,grep { $_->{mask} & $maskcoding } @{ $hg->{array} } );
			}

			next;
		}
		push( @$list2, $id );
	}
	$res->{genes} = $agenes;
#	warn "\t end -->  "
#	  . $tsum
#	  . ' global '
#	  . abs( time - $tglobal ) . ' '
#	  . scalar( keys %$dchr )
#	  if $print;
	return $res unless @$list2;

	###RETURN ####
	warn Dumper @$list2;
	foreach my $z (@$list2){
		my $vobj = $project->returnVariants($z);
		#warn $vobj->name;
		warn Dumper $vobj->gnomad;
	}
	confess();

}

sub refine_heterozygote_composite_score_one {
	my ( $project, $list,$id,$out_header ) = @_;
	
	

	my $out_header;
		$out_header .= $cgi->start_Tr( { style => "background-color:#E9DEFF" } );
		foreach my $h (@headers) {
			$out_header .= $cgi->th($h);
		}
		$out_header .= $cgi->th("validations");
		$out_header .= $cgi->th("transcripts");
		$out_header .= $cgi->end_Tr();
		


	#	warn "genes ".scalar(@$list);

#	$project->buffer->close_lmdb();

	#return;
	#$project->buffer->dbh_reconnect();
	my $mother = $patient->getFamily->getMother();
	my $father = $patient->getFamily->getFather();

#return $list unless $mother or $father;
#	my $no  = GenBoNoSqlLmdb->new(dir=>$dir_tmp,mode=>"w",name=>"tutu",is_compress=>1);
#if ($patient->isChild && $patient->getFamily->isTrio()){
	my $hno;
	my $tsum = 0;
	my $t    = time;
	my $xp   = 0;

	my $total_time = 0;
	my $ztotal     = 0;

	
	
	#$t = time;
	my $current;
	my $rtime =  0;
	foreach my $g ( @$list) {
		#last if $xp > 100;
		$xp++;
		#warn $xp;
		print "*" if $xp % 10 == 0 && $id ==1;
#		warn "*" if $xp % 10 == 0;
		my ( $n, $cname ) = split( "_", $g->{id} );
		my $chr = $project->getChromosome($cname);
		if ( $current ne $cname && $current ) {
			#$project->buffer->close_lmdb();
		}

		$cname = $current;    #unless $cname;
		next unless scalar(keys %{$g->{all_variants}});
		my $no       = $chr->lmdb_polyviewer_variants( $patient, "r" );
		my $out;

		$out .= $cgi->start_div(
			{
				class => "panel panel-primary ",
				style =>
"border-color:white;-webkit-border-radius: 3px;-moz-border-radius: 3px;border-radius: 3px;border: 1px solid black;"
			}
		);
		$out .= $cgi->start_div(
			{
				class => "panel-heading panel-face panel-grey",
				style =>
"$bgcolor;min-height:13px;max-height:13px;padding:10px;border:0px"
			}
		);
		my $panel_id = "panel_" . $g->{uid};
		$out .= update_variant_editor::panel_gene( $g, $panel_id, $project->name, $patient );
		$out .= $cgi->end_div();

		$out .= "\n";
		$out .= $cgi->start_div(
			{
				class => "panel-body panel-collapse collapse ",
				style => "font-size: 09px;font-family:  Verdana;",
				id    => "$panel_id"
			}
		);
		$out .= "\n";

		#$out.="<br>\n";
		$out .= $cgi->start_table(
			{
				class =>
"table table-striped table-condensed table-bordered table-hover table-mybordered",
				style =>
"vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;line-height: 25px;min-height: 25px;height: 25px;box-shadow: 3px 3px 5px #555;"
			}
		);
		$out .= "\n";
		$out .= $out_header;
		
		foreach my $vid ( keys %{ $g->{all_variants} } ) {
			my $v;
			$v = $no->get($vid);
			if ($v->{value}->{is_cnv} == 1) {
				$out .= "<tr style='background-image: linear-gradient(to right, white, #f9e1d8);'>" if ($v->{value}->{type} eq 'large_deletion');
				$out .= "<tr style='background-image: linear-gradient(to right, white, #d8eef9);'>" if ($v->{value}->{type} eq 'large_duplication');
			}
			else { $out .= $cgi->start_Tr(); }
			
#			$out .= $cgi->start_Tr();
		
			#$v = $lite->get($chr->name(),$vid."@".$patient->name);
			my $ttime = time;
			 $rtime += abs(time -$ttime);
			unless ($v){
				$v = hvariant::hash_variant_2($patient,$vid);
			}
			#warn keys %$v;
			#die();

			warn "$vid -- " unless $v;
			confess()       unless $v;
			next            unless $v;

			die($vid) unless $v;

			if ( exists $hno->{$vid} ) {
				$v->{composite} = 1;

				#$v = $hno->{$vid};#if exists $hno->{$vid};

			}
			my $style = {};
			$style = { style => "background-color: #DAEEED;opacity:0.5" } if exists $g->{all_variants}->{$vid}->{added};
			update_variant_editor::alamut_link_js( $v, $patient );
			update_variant_editor::table_validation( $patient, $v, $g );
			$v->{id} = $vid;
			update_variant_editor::table_dejavu_live( $v, $patient->project,$patient,$g );
			my $t = time;

			$total_time += abs( time - $t );

			$v->{id} = $vid;


			update_variant_editor::trio( undef, $v, $patient, $project ) unless $patient->isChild;
			
	

			foreach my $t (@headers) {
				if ($t =~ 'igv' and $v->{value}->{type} ne 'substitution') {
					my $chrom = $v->{value}->{chromosome};
					my $start = $v->{value}->{start};
					my $end = $start;
					my $locus1 = $chrom.':'.$start;
					my $locus2 = $chrom.':'.$start;
					if ($v->{value}->{type} eq 'deletion' or $v->{value}->{type} eq 'large_deletion' or $v->{value}->{type} eq 'large_duplication') {
						my $start2 = $start - 100;
						my $end2 = $v->{value}->{end} + 100;
						$locus2 = $chrom.':'.$start2.'-'.$end2;
					}

					$v->{html}->{$t} =~ s/$locus1/$locus2/ unless ($v->{html}->{$t} =~ /$locus2/);
				}
				if ($t eq 'var_name') {
					if (exists $v->{value}->{manta}->{is_imprecise} and $v->{value}->{manta}->{is_imprecise} == 1) {
						$v->{html}->{$t} .=  qq{<br><span style='font-size:7px'><b><u>IS IMPRECISE</b></u></span>};
					}
					if (exists $v->{value}->{cnv_details_genes} and scalar keys %{$v->{value}->{cnv_details_genes}} > 1) {
						my @lGenesNames;
						push(@lGenesNames, "<table class='table table-striped table-condensed table-bordered table-hover table-mybordered'>");
						push(@lGenesNames, "<tr>");
						push(@lGenesNames, "<td style='text-align:center;'><b>Gene Name</b></td>");
						push(@lGenesNames, "<td style='text-align:center;'><b>Nb Panel(s)</b></td>");
						push(@lGenesNames, "<td style='text-align:center;'><b>Description</b></td>");
						push(@lGenesNames, "</tr>");
						foreach my $g_id (keys %{$v->{value}->{cnv_details_genes}}) {
							push(@lGenesNames, "<tr>");
							my $cmd = qq{"view_var_from_proj_gene_pat('$project_name', '$g_id', '$patient_name', '', '', '');"};
							my $disabled = '';
							$disabled = 'disabled' if ($g_id eq $g->{id});
							my $this_g = $project->newGene($g_id);
							
							my $nb_panels = scalar(@{$this_g->getPanels()});
							my $panels_text = '';
							foreach my $p (@{$this_g->getPanels()}) {
								$panels_text += $p->name()."\n";
							}
							my $cmd_alert_panel = qq{alert("$panels_text");};
							push(@lGenesNames, "<td style='text-align:center;'><button onclick=$cmd $disabled>".$this_g->external_name()."</button></td>");
							push(@lGenesNames, "<td style='text-align:center;'>".$nb_panels."</td>");
							push(@lGenesNames, "<td style='text-align:center;'>".$this_g->description()."</td>");
							push(@lGenesNames, "</tr>");
						}
						push(@lGenesNames, "</table>");
						#TODO: here
						my $genes_text = join("", @lGenesNames);
						my $id_info = 'b_multi_genes'.$v->{value}->{id}.'_'.$g->{id}.'_'.$patient->name();
						
						$v->{html}->{$t} .= qq{<br><br><button id="$id_info" type="button" class="btn btn-xs  btn-primary" style="background-color: #9796C4;font-size: 7px;font-family:  Verdana;color:white">Multi Genes !</button>};
						$v->{html}->{$t} .= qq{<div data-dojo-type="dijit/Tooltip" data-dojo-props="connectId:'$id_info',position:['after']"><div style="width:450px;height:auto;text-align:center;">$genes_text</div></div>};
					}
				}
				if ($t eq 'trio') {
					my $cnv_score;
					my $max_dp = $v->{value}->{max_dp};
					my $max_dp_text = $max_dp;
					$max_dp_text = '-' if ($max_dp_text == -1);
					my $text_caller = join("<br>", @{$v->{value}->{ngs}});
					my $id_info = 'table_infos_'.$v->{value}->{id}.'_'.$g->{id}.'_'.$patient->name();
					my $vtype = '';
					$vtype = '<span style="text-align:center;font-size:7px"><b><i>Large Deletion</b></i></span><br><br>' if ($v->{value}->{type} eq 'large_deletion');
					$vtype = '<span style="text-align:center;font-size:7px"><b><i>Large Duplication</b></i></span><br><br>' if ($v->{value}->{type} eq 'large_duplication');
					$v->{html}->{$t} =~ s/<table/$vtype<table id="$id_info"/;
					$v->{html}->{$t} .= qq{<div data-dojo-type="dijit/Tooltip" data-dojo-props="connectId:'$id_info',position:['after']"><div style="min-width:300px;width:auto;height:auto;"><center><b><u>Patient $patient_name</b></u><br><br>};

					if (exists $v->{value}->{manta}->{is_imprecise} and $v->{value}->{manta}->{is_imprecise} == 1) {
						$text_caller .=  qq{<br><br><b><u>IS IMPRECISE (MANTA)</b></u>};
					}

					$v->{html}->{$t} .= qq{<b><u>Calling Method(s):</b></u><br> $text_caller </center></div></div>};
					$v->{html}->{$t} .= qq{$cnv_score} if ($cnv_score);
					
					
					
				}
				$out .= $cgi->td( $style, $v->{html}->{$t} );
				$out .= "\n";
			}
			
			update_variant_editor::check_is_hgmd_dm_for_gene($patient, $v, $g); #TODO: if DM a faire
			update_variant_editor::check_is_clinvar_pathogenic_for_gene($patient, $v, $g); #TODO: if CLINVAR a faire
			$out .= $cgi->td( $style,
			update_variant_editor::table_validation( $patient, $v, $g ) );
			if ($v->{value}->{type} eq 'large_deletion' or $v->{value}->{type} eq 'large_duplication') {
				$out .= $cgi->td($style, update_variant_editor::table_transcripts_cnv($v, $v->{genes}->{ $g->{id} }, $g->{id}, \@header_transcripts_cnv));
			}
			else {
				$out .= $cgi->td($style, update_variant_editor::table_transcripts($v->{genes}->{ $g->{id} }, \@header_transcripts));
			}
			$out .= $cgi->td($style, update_variant_editor::validation_select( $patient, $v, $g ) );
			$out .= $cgi->end_Tr();
			$out .= "\n";

		}
		
		#$chr->lmdb_polyviewer_variants( $patient, "close" );
		
		$out .= $cgi->end_table();
		$out .= "\n";
		$out .= $cgi->end_div();
		$out .= $cgi->end_div();
		$out .= "<br>\n";
		$g->{out} = $out;

		#print $out;
	}
	$project->buffer->close_lmdb();
#	warn "refine get: $rtime => step2 :".abs(time-$t);
	#$no->close();
	return ( $list, $total_time );

	die();
}

sub new_refine_heterozygote_composite_score_old {
	my ( $project, $list, $id, $out_header ) = @_;

	my $out_header;
	$out_header .= $cgi->start_Tr( { style => "background-color:#E9DEFF" } );
	foreach my $h (@headers) {
		$out_header .= $cgi->th($h);
	}
	$out_header .= $cgi->th("validations");
	$out_header .= $cgi->th("transcripts");
	$out_header .= $cgi->end_Tr();
	my $mother = $patient->getFamily->getMother();
	my $father = $patient->getFamily->getFather();

	my $hno;
	my $tsum = 0;
	my $t    = time;
	my $xp   = 0;

	my $total_time = 0;
	my $ztotal     = 0;

	#$t = time;
	my $current;
	my $rtime = 0;
	foreach my $g (@$list) {

		#last if $xp > 100;
		$xp++;

		#warn $xp;
		print "*" if $xp % 10 == 0 && $id == 1;
		my ( $n, $cname ) = split( "_", $g->{id} );
		my $chr = $project->getChromosome($cname);
		if ( $current ne $cname && $current ) {

			#$project->buffer->close_lmdb();
		}

		$cname = $current;    #unless $cname;
		next unless scalar( keys %{ $g->{all_variants} } );

		#my $no       = $chr->lmdb_polyviewer_variants( $patient, "r" );

		my $noV = $chr->get_lmdb_variations("r");
		my $no  = $chr->lmdb_polyviewer_variants( $patient, "r" );
		my $out;

		$out .= $cgi->start_div(
			{
				class => "panel panel-primary ",
				style =>
"border-color:white;-webkit-border-radius: 3px;-moz-border-radius: 3px;border-radius: 3px;border: 1px solid black;"
			}
		);
		$out .= $cgi->start_div(
			{
				class => "panel-heading panel-face panel-grey",
				style =>
"$bgcolor;min-height:13px;max-height:13px;padding:10px;border:0px"
			}
		);
		my $panel_id = "panel_" . $g->{uid};
		$out .=
		  update_variant_editor::panel_gene( $g, $panel_id, $project->name,
			$patient );
		$out .= $cgi->end_div();

		$out .= "\n";
		$out .= $cgi->start_div(
			{
				class => "panel-body panel-collapse collapse ",
				style => "font-size: 09px;font-family:  Verdana;",
				id    => "$panel_id"
			}
		);
		$out .= "\n";

		#$out.="<br>\n";
		$out .= $cgi->start_table(
			{
				class =>
"table table-striped table-condensed table-bordered table-hover table-mybordered",
				style =>
"vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;line-height: 25px;min-height: 25px;height: 25px;box-shadow: 3px 3px 5px #555;"
			}
		);
		$out .= "\n";
		$out .= $out_header;

		foreach my $vid ( keys %{ $g->{all_variants} } ) {
		#	my $v  = $noV->get($vid);
			my $vh = $no->get($vid);
			my $vp =  PolyviewerVariant->new();
			$vp->setOldVariant($vh,$project,$patient,$g);
			#$vp->setLmdbVariant($vh,$project,$g,$patient);
			$print_html->variant($vp);
			


			my $ttime = time;
			$rtime += abs( time - $ttime );
			unless ($vh) {
				confess();
				$vh = hvariant::hash_variant_2( $patient, $vid );
			}

			if ( exists $hno->{$vid} ) {
				$vh->{composite} = 1;

			}
			my $style = {};
			$style = { style => "background-color: #DAEEED;opacity:0.5" } if exists $g->{all_variants}->{$vid}->{added};

			my $t = time;

			$total_time += abs( time - $t );

		
			my $hpatients;
		
		
			#my $is_gnomad = exists $v->{value}->{ac};
			

			my @headers = (
				"varsome", "igv",    "alamut", "var_name",
				"trio",    "gnomad", "deja_vu"
			);
			
			
			##############
			# VARSOME CELL
			###############
			my $t1 = shift(@headers);
			$out .= $cgi->td( $style, $print_html->varsome() );
			$out .= "\n";

			##############
			# IGV CELL
			###############

			$t = shift(@headers);

			#write locus
			$out .= $cgi->td( $style, $print_html->igv);

			##############
			# ALAMUT CELL
			###############

			$t = shift(@headers);

			$out .= $cgi->td( $style,$print_html->alamut);
			$out .= "\n";

			##############
			# NAME CELL
			#
			###############
			$t = shift(@headers);

			#$name =  $v->{var_name} if exists $v->{var_name};

			$out .= $cgi->td($style,$print_html->var_name());

			$out .= "\n";
			##############
			# CELL CALLING INFOS
			###############

			
				$out .= $cgi->td( $style, $print_html->calling()) ;

			$out .= "\n";
			
			
			
			
			$t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->gnomad() );
			$out .= "\n";


			$t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->dejavu() );
			$out .= "\n";
			$t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->validations );
			
			$t = shift(@headers);
			$out .= "\n";
			$out .= $cgi->td( $style, $print_html->transcripts() );
			$out .= "\n";

			$out .= $cgi->td( $style, $print_html->validation_select() );
			$out .= $cgi->end_Tr();
			$out .= "\n";

		}

		$out .= $cgi->end_table();
		$out .= "\n";
		$out .= $cgi->end_div();
		$out .= $cgi->end_div();
		$out .= "<br>\n";
		$g->{out} = $out;

	}
	$project->buffer->close_lmdb();
	return ( $list, $total_time );

	die();
}

sub new_refine_heterozygote_composite_score_one {
	my ( $project, $list, $id, $out_header ) = @_;

	my $out_header;
	$out_header .= $cgi->start_Tr( { style => "background-color:#E9DEFF" } );
	foreach my $h (@headers) {
		$out_header .= $cgi->th($h);
	}
	$out_header .= $cgi->th("validations");
	$out_header .= $cgi->th("transcripts");
	$out_header .= $cgi->end_Tr();
	my $mother = $patient->getFamily->getMother();
	my $father = $patient->getFamily->getFather();

	my $hno;
	my $tsum = 0;
	my $t    = time;
	my $xp   = 0;

	my $total_time = 0;
	my $ztotal     = 0;

	#$t = time;
	my $current;
	my $rtime = 0;
	foreach my $g (@$list) {

		#last if $xp > 100;
		$xp++;

		#warn $xp;
		print "*" if $xp % 10 == 0 && $id == 1;
		my ( $n, $cname ) = split( "_", $g->{id} );
		my $chr = $project->getChromosome($cname);
		if ( $current ne $cname && $current ) {

			#$project->buffer->close_lmdb();
		}

		$cname = $current;    #unless $cname;
		next unless scalar( keys %{ $g->{all_variants} } );
		my $no  = $chr->lmdb_polyviewer_variants( $patient, "r" );
		my $noV = $chr->get_lmdb_variations("r");
		my $out;

		$out .= $cgi->start_div(
			{
				class => "panel panel-primary ",
				style =>
"border-color:white;-webkit-border-radius: 3px;-moz-border-radius: 3px;border-radius: 3px;border: 1px solid black;"
			}
		);
		$out .= $cgi->start_div(
			{
				class => "panel-heading panel-face panel-grey",
				style =>
"$bgcolor;min-height:13px;max-height:13px;padding:10px;border:0px"
			}
		);
		my $panel_id = "panel_" . $g->{uid};
		$out .=
		  update_variant_editor::panel_gene( $g, $panel_id, $project->name,
			$patient );
		$out .= $cgi->end_div();

		$out .= "\n";
		$out .= $cgi->start_div(
			{
				class => "panel-body panel-collapse collapse ",
				style => "font-size: 09px;font-family:  Verdana;",
				id    => "$panel_id"
			}
		);
		$out .= "\n";

		#$out.="<br>\n";
		$out .= $cgi->start_table(
			{
				class =>
"table table-striped table-condensed table-bordered table-hover table-mybordered",
				style =>
"vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;line-height: 25px;min-height: 25px;height: 25px;box-shadow: 3px 3px 5px #555;"
			}
		);
		$out .= "\n";
		$out .= $out_header;

		foreach my $vid ( keys %{ $g->{all_variants} } ) {
			my $v = $noV->getHash( $vid, $patient );
			next;
			my $vh = $no->get($vid);
			foreach my $sid ( keys %{ $v->{patients} } ) {
				$v->{patients}->{$sid}->{model} =
				  $vh->{patients}->{$sid}->{model};
			}

			if ( $v->{value}->{is_cnv} == 1 ) {
				$out .=
"<tr style='background-image: linear-gradient(to right, white, #f9e1d8);'>"
				  if ( $v->{value}->{type} eq 'large_deletion' );
				$out .=
"<tr style='background-image: linear-gradient(to right, white, #d8eef9);'>"
				  if ( $v->{value}->{type} eq 'large_duplication' );
			}
			else { $out .= $cgi->start_Tr(); }
			my $ttime = time;
			$rtime += abs( time - $ttime );
			unless ($v) {
				confess();
				$v = hvariant::hash_variant_2( $patient, $vid );
			}

			warn "$vid -- " unless $v;
			confess()       unless $v;
			next            unless $v;

			die($vid) unless $v;

			if ( exists $hno->{$vid} ) {
				$v->{composite} = 1;

				#$v = $hno->{$vid};#if exists $hno->{$vid};

			}
			my $style = {};
			$style = { style => "background-color: #DAEEED;opacity:0.5" }
			  if exists $g->{all_variants}->{$vid}->{added};
			update_variant_editor::alamut_link_js( $v, $patient );

			#update_variant_editor::table_validation( $patient, $v, $g );
			$v->{id} = $vid;

 #update_variant_editor::table_dejavu_live( $v, $patient->project,$patient,$g );
			my $t = time;

			$total_time += abs( time - $t );

			$v->{id} = $vid;

#update_variant_editor::trio( undef, $v, $patient, $project ) unless $patient->isChild;
#update_variant_editor::check_is_hgmd_dm_for_gene($patient, $v, $g); #TODO: if DM a faire
#update_variant_editor::check_is_clinvar_pathogenic_for_gene($patient, $v, $g); #TODO: if CLINVAR a faire
			my @headers = (
				"varsome", "igv",    "alamut", "var_name",
				"trio",    "gnomad", "deja_vu"
			);
			my $t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->varsome($v) );
			$out .= "\n";
			$t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->igv($v) );
			$out .= "\n";
			$t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->alamut($v) );
			$out .= "\n";
			$t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->var_name($v) );
			$out .= "\n";
			$out .= "\n";
			$t = shift(@headers);

			#$out .= $cgi->td( $style,$v->{html}->{$t});
			#	if ($v->{name} eq "2-179312317-ins-100"){
			#		warn Dumper $v;
			#		die();
			#	}
			if ( exists $v->{value}->{large_evt} ) {
				$out .= $cgi->td( $style, $print_html->calling_cnv( $v, $g ) );
			}
			else {
				$out .= $cgi->td( $style, $print_html->calling( $v, $g ) );
			}

			$out .= "\n";
			$t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->gnomad($v) );
			$out .= "\n";
			$t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->dejavu($v) );
			$out .= "\n";

			$t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->validations( $v, $g ) );
			$t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->transcripts( $v, $g ) );
			$out .= "\n";

			$out .= $cgi->td( $style,
				update_variant_editor::validation_select( $patient, $v, $g ) );
			$out .= $cgi->end_Tr();
			$out .= "\n";

		}

		$out .= $cgi->end_table();
		$out .= "\n";
		$out .= $cgi->end_div();
		$out .= $cgi->end_div();
		$out .= "<br>\n";
		$g->{out} = $out;

	}
	$project->buffer->close_lmdb();
	return ( $list, $total_time );

	die();
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
	my ( $project, $patient,$statistics ) = @_;
	unless ($cgi->param('export_xls')) {
		print qq{<div style="display: none">};
		print "vector";
	}
	my $filter_transmission;
	$filter_transmission->{denovo}          = 1 if $cgi->param('denovo');	
	$filter_transmission->{"strict_denovo"} = 1 if $cgi->param('denovo');
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
	$filter_transmission->{xor}             = 1 if  $cgi->param('xor');
	warn "\n";
	
	my $hashVector_panel = {};
	my $hashVector       = {};
	construct_panel_vector( $panel, $hashVector_panel ) if $panel;
	my $list_transcript;
	my $trio = $patient->getFamily->isTrio;
	my $list_variants;
	my $hash_variants_DM = {};
	my $gene;

	if ($gene_name_filtering) {
		$gene = $project->newGene($gene_name_filtering);
		#error("GENE NAME $gene_name_filtering NOT FOUND !!!!!!") unless $gene;
		$gene_id_filtering = $gene->id();
		
	}
		#my @exception = ("PCDH19") ;
	my $gene_exception = $project->newGene("PCDH19");
	my $fork = 12;
	my $pm   = new Parallel::ForkManager($fork);
	my $hrun;
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $h ) = @_;

			unless ( defined($h) or $exit_code > 0 ) {
				print
				  qq|No message received from child process $exit_code $pid!\n|;
				die();
				return;
			}
			my $chr = $h->{chromosome};
			foreach my $k  (keys %{$h->{statistics}}){
				$statistics->{$k} += $h->{statistics}->{$k};
				
			}
			unless ( exists $hashVector->{$chr} ) {
				$hashVector->{$chr} = $h->{vector};
			}
			else {
				$hashVector->{$chr} &= $h->{vector};
			}
			my $id = $h->{run_id};
			delete $hrun->{ $h->{run_id} };
			error("PROBLEM !!!!") unless exists  $h->{run_id};
			push( @$list_variants, @{ $h->{list_variants} } ) if  $h->{list_variants};
			map { $hash_variants_DM->{$_}++ } @{ $h->{list_variants_DM} };
			
			
		}
		
	);
	my $id = time;
	$project->buffer->dbh_deconnect();

	
	

	foreach my $chr ( @{ $project->getChromosomes } ) {
		#warn $chr->name;
		if ($gene) {
			next if ( $gene->getChromosome()->name ne $chr->name );
		}
		if ($panel) {
			next unless exists $hashVector_panel->{ $chr->name };
		}
		$id++;
		$id ++;
		$hrun->{$id} ++;
		my $pid = $pm->start and next;
		$project->buffer->dbh_reconnect();
		my $debug;
		$debug =1 if $chr->name eq "12";
		my $statistics = {};
		my $res;
		$res->{run_id} = $id . "_" . $chr->name;
		my $hashVector = {};
		$res->{list_variants_DM} = [];
		
		if ($gene) {

			$hashVector->{ $chr->name } = $gene->getVectorOrigin();

		}
		elsif ($panel) {
			$hashVector->{ $chr->name } =
			  $hashVector_panel->{ $chr->name }->Clone;

		}
		else {
			$hashVector->{ $chr->name } = $chr->getVariantsVector();
		}
		
		#my $vquality = $chr->getVectorScore($vector_ratio_name);
		print "=" unless $cgi->param('export_xls');
			my $testid = 4201;
		my $debug;

		#	$debug =1 if $chr->name eq  "12";
		if ($panel) {

			#next unless  (exists $hashVector->{$chr->name});
			$hashVector->{ $chr->name } &= $chr->getVectorScore($limit_ac);    #if $limit_ac ne "all";
			
			#$hashVector->{ $chr->name } &=
			 # $chr->getVectorScore($limit_ac_ho);    # if $limit_ac_ho ne "all";
			$hashVector->{ $chr->name } &=
			  $chr->getVectorScore($limit_sample_dv)
			  ;    # if $limit_sample_dv ne "all";
			 $hashVector->{$chr->name} &= $chr->getVectorScore($limit_sample_dv_ho);
			#$hashVector->{ $chr->name } &= $vquality if $vquality;
			#$hashVector->{ $chr->name } &= $vquality if $vquality;

			#$hashVector->{$chr->name} &= $chr->getVectorLargeDeletions;
			$hashVector->{ $chr->name } -= $chr->getVectorScore("intergenic");

		}

		else {
			$hashVector->{ $chr->name } &= $chr->getVariantsVector();
			warn "\n 1 ==> " . $hashVector->{$chr->name}->contains($testid)  if $debug;
			$hashVector->{ $chr->name } = $chr->getVectorScore($limit_ac);    # if $limit_ac ne "gnomad_ac_all";;
			warn "\n 1 ==> " . $hashVector->{$chr->name}->contains($testid)  if $debug;
			# warn $limit_sample_dv_ho;
			$hashVector->{ $chr->name } &= $chr->getVectorScore($limit_ac_ho);    # if $limit_ac_ho ne "all";
			warn "\n 1 ==> " . $hashVector->{$chr->name}->contains($testid)  if $debug;
			$hashVector->{ $chr->name } &= $chr->getVectorScore($limit_sample_dv);    #  if $limit_sample_dv ne "all";
			warn "\n 1 ==> " . $hashVector->{$chr->name}->contains($testid)  if $debug;
			 $hashVector->{$chr->name} &= $chr->getVectorScore($limit_sample_dv_ho);
			 #$hashVector->{$chr->name} &= $chr->getVectorLargeDeletions;
			#$hashVector->{ $chr->name } &= $vquality if $vquality;
			$hashVector->{ $chr->name } -= $chr->getVectorScore("intergenic");
			warn "\n 1 ==> " . $hashVector->{$chr->name}->contains($testid)  if $debug;
			
		}
	#$hashVector->{ $chr->name } -=  $chr->getVectorLargeDeletions();
	#$hashVector->{ $chr->name } -=  $chr->getVectorLargeDuplications();

  #warn $patient->countThisVariants( $patient->getVectorOrigin($chr)) if $debug;
  
  warn "\n 1 ==> " . $hashVector->{$chr->name}->contains($testid)  if $debug;
  warn "\n\n+".$patient->getVectorOrigin($chr)->contains($testid) if $debug;
  #warn $chr->name;
		$hashVector->{ $chr->name } &= $patient->getVectorOrigin($chr);
		
		$patient->countThisVariants( $patient->getVectorOrigin($chr) )
		  if $debug;

   warn "2" if 	$debug && $hashVector->{$chr->name}->contains($testid) && $debug;

		$statistics->{variations} += $patient->countThisVariants( $hashVector->{ $chr->name } );

		my $vDM = $chr->vectorDM();
		$vDM += $chr->vectorClinvarPathogenic();
		$vDM &= $patient->getVectorOrigin($chr);
		#$vDM &= $hashVector->{ $chr->name };
		$statistics->{DM} += $patient->countThisVariants($vDM);
		#$hashVector->{$chr->name}  |= $v1;

		if ($trio) {
			my @list_transmission = keys %$filter_transmission;
			my $vtr               = $chr->getNewVector();
			
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
				#my @type = ["both","strict_denovo","denovo","recessive","xor","xor_mother","xor_father"];
				
				foreach my $tr ( keys %$filter_transmission ) {
						print "__";
					if ( $tr eq "both" ) {
						 $vtr = $patient->getVectorOrigin($chr);
						last;
					}
				#$vtr |= $patient->getFamily()->getVector_family_dominant($chr);
					if ( $tr eq "strict_denovo" ) {
						$vtr |= $patient->getFamily()->getVectorStrictDenovoTransmission( $chr, $patient );

					}
					elsif ( $tr eq "denovo" ) {
						$vtr |= $patient->getFamily()->getVectorDenovoTransmission( $chr, $patient );

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
						#	 warn "trio $tr =>".$hashVector->{$chr->name}->contains($testid)  if 	$debug && $debug;	
							$vtr2 -=  $patient->getFamily()->getVectorRecessiveTransmission( $chr, $patient ) if $patient->getFamily->getMother;
						}
						$vtr |= $vtr2;
					}
					elsif ( $tr eq "xor_mother" ) {
#						#my $v1 = $patient->getVectorOrigin($chr);
#						$v1 &= $patient->getFamily()->getVectorMotherTransmission( $chr, $patient ) if $patient->getFamily->getMother;
#						$v1 -= $patient->getFamily->getFather->getVectorOrigin( $chr,$patient ) if $patient->getFamily->getFather();
#						$v1 -= $patient->getFamily()->getVectorRecessiveTransmission( $chr, $patient );
					#	warn "mother";
						$vtr |= $patient->getFamily()->getVector_individual_mother( $chr, $patient );
					}
					elsif ( $tr eq "xor_father" ) {
					#	warn "father";
					#	my $v1 = $patient->getVectorOrigin($chr);
					#	$v1 &= $patient->getFamily()->getVectorFatherTransmission( $chr, $patient ) if $patient->getFamily->getFather;
					#	$v1 -= $patient->getFamily->getMother->getVectorOrigin( $chr ) if $patient->getFamily->getMother();
					#	$v1 -= $patient->getFamily()->getVectorRecessiveTransmission( $chr, $patient );
						$vtr |= $patient->getFamily()->getVector_individual_father( $chr, $patient );
					}
				}
				
				if ($chr->name eq $gene_exception->getChromosome->name and $patient->getVectorOrigin($chr)){
					my $vector = $gene_exception->getVectorOrigin();
				 	$vector &= $patient->getVectorOrigin($chr);
					my $father = $patient->getFamily->getFather();
					my $mother = $patient->getFamily->getMother();
					if ($father){
						my $v1 = $father->getVectorOrigin($chr);
						$vector &= $v1;
						if ($mother){
							my $v2  = $mother->getVectorOrigin($chr);
							$vector -= $v2;
						}
						 $vtr |= $vector;
					}
					
				}
			
				
			}
			
			$hashVector->{ $chr->name } &= $vtr;

		}
		#  die($chr->name) if $debug && $hashVector->{ $chr->name }->contains(7710) && $debug;
		$hashVector->{ $chr->name } &= $hashVector_panel->{ $chr->name } if exists $hashVector_panel->{ $chr->name };
		  $statistics->{variations} += $patient->countThisVariants( $hashVector->{ $chr->name } );
		$res->{chromosome} = $chr->name;
		if($only_DM){
				#keep only pathogenic variation if cgi with option $only_DM
				$hashVector->{ $chr->name } = $vDM;
		}
		if ($keep_pathogenic) {
			$vDM &= $hashVector->{ $chr->name };
			$hashVector->{ $chr->name } |= $vDM;
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
				 $vquality = $chr->getVectorScore($vector_ratio_name);
			}
			
			if ($limit_ratio2 != $limit_ratio1) {
				
				my $no = $chr->lmdb_polyviewer_variants( $patient, "r" );
				my $no2 = $chr->lmdb_polyviewer_variants_genes( $patient, "r" );
				my $vector_ratio_name = $patient->name . "_ratio_" . $limit_ratio2;
				$vector_ratio_name = $patient->name . "_ratio_all" if $limit_ratio2 == -1;
				 my $vquality2 = $chr->getVectorScore($vector_ratio_name);
				 $vquality2 -= $vquality;
				
				$vquality2 &=  $res->{vector};
				$res->{vector} &= $vquality;	
				#warn $vquality;
				#warn $vquality;
				foreach my $id ( @{ to_array( $vquality2, $chr->name ) } ) {
					print "!";
					my $av = $no->get($id);
					#warn Dumper $av;
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
		
		
		$res->{statistics}    = $statistics;
		$res->{list_variants} = [];

		my $vLocal = $chr->getVectorLocalValidations();
		

		$vLocal &= $patient->getVectorOrigin($chr);

		#Limitation gnomad ho ac a 1000 si projet exome ou genome
		if ( $project->isExome() || $project->isGenome() ) {
			$vLocal &= $chr->getVectorScore("gnomad_ho_ac_1000");
		}

		$res->{vector} |= $vLocal;
		my $xxx;
		#
		foreach my $id ( @{ to_array( $res->{vector}, $chr->name ) } ) {
			$xxx++;
			print "++" if $xxx % 100 ==0;
			push( @{ $res->{list_variants} }, $id );
		}
		if ($keep_pathogenic){
		$vDM |= $vLocal if $keep_pathogenic;
		foreach my $id ( @{ to_array( $vDM, $chr->name ) } ) {
		#	warn $id;
		#	my $vobj   = $project->returnVariants($id);
		#	warn $vobj->name;
			push( @{ $res->{list_variants_DM} }, $id );
		}
		}
		$res->{run_id} = $id;
		
		$pm->finish( 0, $res );
	}
	$pm->wait_all_children();
	error("ARGGG Problem") if keys %$hrun;
	$project->buffer->dbh_reconnect();
	print qq{</div>} unless $cgi->param('export_xls');
	return ( $hashVector, $list_variants, $hash_variants_DM );
}

sub listVariants {
	my ($vectors) = @_;
	my @list_variants;

	my $already;
	my $not_saved = 0;
	foreach my $k ( keys %$vectors ) {
		print "@";

		foreach my $id ( @{ to_array( $vectors->{$k}, $k ) } ) {
			my ( $a, $b ) = split( "!", $id );
		
			push( @list_variants, $id );
		}
	}
	return ( \@list_variants, $already, $not_saved );
}

sub to_array {
	my ( $v, $name ) = @_;
	my $set  = Set::IntSpan::Fast::XS->new( $v->to_Enum );
	my $iter = $set->iterate_runs();
	my @t;
	while ( my ( $from, $to ) = $iter->() ) {
		for my $member ( $from .. $to ) {
			if ($name) {
				push( @t, $name . "!" . $member );
			}
			else {
				push( @t, $member );
			}
		}
	}
	return \@t;
}

sub return_max_variant_score {
	my ( $vector, $gene, $patient ) = @_;
	my $no  = $gene->getChromosome->lmdb_polyviewer_genes( $patient, "r" );
	my $chr = $gene->getChromosome();
	my @a;
	my $avs = to_array($vector);
	my $hmax;
	foreach my $vid ( @{$avs} ) {

		my $hg =
		  $chr->getPolyviewer_score( $patient, "s:" . $vid . ":" . $gene->id );
		unless ($hg) {

			my $gid = $chr->cache_lmdb_variations->get_varid($vid);
			my $xid = $chr->name() . "!" . $vid;

			#warn "s:".$vid.":".$gene->id."-".$patient->name;

			#warn $xid;
			my $v = $gene->project->returnVariants($xid);
			foreach my $p ( @{ $v->getPatients() } ) {

				#warn  "\t".$p->name();
			}
			my $vv = $no->get($xid);

			die();
		}
		die( $vid . " " . $gene->getChromosome->name . " " . $patient->name )
		  unless $hg;
		die() unless $hg->{id};
		$hmax = $hg unless $hmax;
		$hmax = $hg if $hg->{score} > $hmax->{score};
	}

	#warn "end";

	return $hmax;
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
	my $cno = $project->getChromosomes()->[0]->lmdb_hash_variants("r");



my $date;
( $date->{cache} ) = utility::return_date_from_file( $cno->filename );
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
