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

require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";

my $io = IO::Handle->new();
$io->autoflush(1);


my $cgi = new CGI();
my $user_name = $cgi->param('user');
my $pwd = $cgi->param('pwd');
my $sort_by_captures = $cgi->param('sort_by_captures');

$user_name = lc($user_name);
my $buffer_init = new GBuffer;
my $can_use_hgmd = $buffer_init->hasHgmdAccess($user_name);



my $hRes;
$hRes->{ok} = 1;
my @lItemsProjects;
my $hProjects = get_hash_users_projects($user_name, $pwd);

my @header_transcripts = ("consequence","enst","nm","ccds","appris","exon","nomenclature","codons","codons_AA", "polyphen","sift","ncboost","cadd","revel","dbscsnv");
my ($hVariantsIdsDejavu, $hVariantsDetails);
my ($hResGene, $hResVariants, $hResVariantsIds, $hResVariantsListPatients, $hResProjectProblem, $hResVariants_byScores);
my $project_init_name;

#project preload_patients test connection DB
my $dbh_init = $buffer_init->dbh();
my $query_init = $buffer_init->getQuery();
my $sql_prepare;
my $project_init;


my $path_dv_date_test =  $buffer_init->deja_vu_public_dir('HG38',"variations").'/1.g.rocks/chunk.json';
my $date_res = `ls -ltrh $path_dv_date_test`;
chomp($date_res);
my @ltmp = split(' ', $date_res);
my $date_dejavu = 'Last DV update: '.$ltmp[-3].' '.$ltmp[-4];

my $out =  $cgi->start_div({style=>"font-size:10px;"});
$out .= qq{<center>};
$out .= qq{<br>};
$out .=  $cgi->end_div();
$out .=  $cgi->start_div();
$out .= "<table class='table table-striped' style='font-size:13px;'>";
$out .= "<thead>";
$out .= $cgi->start_Tr({style=>"background-color:#E9DEFF;font-size:10px;"});
#$out .=  $cgi->th("<center>#</center>");
$out .=  $cgi->th("<center>Found DV Parquet</center>");
$out .=  $cgi->th("<center>Status</center>");
$out .=  $cgi->th("<center>Name</center>");
$out .=  $cgi->th("<center>Id</center>");
$out .=  $cgi->th("<center>Description</center>");
$out .= $cgi->end_Tr();
$out .= "</thead>";
$out .= "<tbody>";
my ($h_p, $h_c);
my $hCapturesNames;
my $hCapturesNamesProject;
my @lProjectNames = reverse sort keys %{$hProjects};
my $i_p = 0;

my $dir_pr = $buffer_init->config_path("root","dejavu")."/HG38/variations/rocks/";
my $ro_projects = GenBoNoSqlRocks->new(mode=>"r",dir=>$dir_pr, name=>"projects");

foreach my $project_name (@lProjectNames) {
	my $id = $hProjects->{$project_name}->{id};
	my $description = $hProjects->{$project_name}->{description};
	$i_p++;
	
	$h_p->{$project_name} = '';
	$h_p->{$project_name} .=  $cgi->start_Tr();
	my ($html_glyph, $status, $style);
	if ($ro_projects->get_raw($project_name)) {
		$html_glyph = qq{<span style='color:green;' class='glyphicon glyphicon-ok'></span>};
		$status = 'OK';
	}
	else {
		my $h_this_proj = $buffer_init->getQuery()->getProjectByName($project_name);
		if ($h_this_proj->{is_somatic} == 1 or $h_this_proj->{dejavu} != 1) {
			$html_glyph = qq{<span style='color:orange;' class='glyphicon glyphicon-ok'></span>};
			$status = 'Excluded';
			$status .= ' Project' if $h_this_proj->{dejavu} != 1;
			$status .= ' (somatic)' if $h_this_proj->{is_somatic} == 1;
			$style = qq{style="color:orange;"};
		}
		else {
			$html_glyph = qq{<span style='color:red;' class='glyphicon glyphicon-exclamation-sign'></span>};
			$status = 'Not in DV';
			$style = qq{style="color:red;"};
		}
	}
	$h_p->{$project_name} .=  $cgi->td("<center>$html_glyph</center>");
	$h_p->{$project_name} .=  $cgi->td("<center><span $style>$status</span></center>");
	$h_p->{$project_name} .=  $cgi->td("<center><span $style>$project_name</span></center>");
	$h_p->{$project_name} .=  $cgi->td("<center><span $style>$id</span></center>");
	$h_p->{$project_name} .=  $cgi->td("<center><span $style>$description</span></center>");
	$h_p->{$project_name} .= $cgi->end_Tr();
}

if ($sort_by_captures) {
	foreach my $c (sort keys %$h_c) {
		foreach my $project_name (sort keys %{$h_c->{$c}}) {
			$out .= $h_p->{$project_name};
		}
	}
}
else {
	foreach my $project_name (@lProjectNames) {
		$out .= $h_p->{$project_name};
	}
}
$out .= "</tbody>";
$out .= "</table></center>";

$ro_projects->close();

my $hCapturesNames_for_sort;
foreach my $c_name (keys %$hCapturesNames) { $hCapturesNames_for_sort->{lc($c_name)} = $c_name; }
my @lCapturesNames = sort keys %$hCapturesNames_for_sort;

my $out_captures;
$out_captures .= qq{<div class="input-group" style="width:100%">};
$out_captures .= qq{<select class="form-control" id="form_my_captures" style="font-size:9px;height:auto;width:100%;">};
$out_captures .= qq{<option value=''><span></span></option>};
foreach my $capture_for_sort (sort @lCapturesNames) {
	my $capture = $hCapturesNames_for_sort->{$capture_for_sort};
	my $name = lc($capture);
	$name =~ s/ /_/g;
	my @lProjects = sort keys %{$hCapturesNamesProject->{$capture}};
	my $value = $name;
	$value .= ';'.$lProjects[-1];
	$name = ucfirst($name);
	my $nb_proj = $hCapturesNames->{$capture};
	$out_captures .= qq{<option value='$value'><span>$name [$nb_proj projects]</span></option>};
}
$out_captures .= qq{</select>};
$out_captures .= qq{</div>};

my $h_pheno;
foreach my $pheno_id (@{$buffer_init->queryPhenotype->getAllPhenotypes()}) {
	my $h = $buffer_init->queryPhenotype->getPhenotypeInfos($pheno_id);
	my $name = $h->{$pheno_id}->{name};
	my $concept = $h->{$pheno_id}->{concept};
	$h_pheno->{$name} = $name;
	if ($concept and $concept ne '') { $h_pheno->{$name} .= ';'.$concept; }
}

my $out_phenos;
$out_phenos .= qq{<div class="input-group" style="width:100%">};
$out_phenos .= qq{<select class="form-control" id="form_phenotypes" style="font-size:9px;height:auto;width:100%;">};
$out_phenos .= qq{<option value=''><span></span></option>};
foreach my $phenotype (sort keys %$h_pheno) {
	my $name = lc($phenotype);
	$name = ucfirst($name);
	my $value = $h_pheno->{$phenotype};
	$out_phenos .= qq{<option value='$value'><span>$name</span></option>};
}
$out_phenos .= qq{</select>};
$out_phenos .= qq{</div>};


$hRes->{date_dejavu} = $date_dejavu;
$hRes->{html_projects} = $out;
$hRes->{list_projects}= \@lProjectNames;
$hRes->{html_list_captures}=$out_captures;
$hRes->{html_list_phenotypes}=$out_phenos;
print $cgi->header('text/json-comment-filtered');
my $json_encode = encode_json $hRes;
print $json_encode;
exit(0);



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
