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

my $out =  $cgi->start_div({style=>"font-size:10px"});
$out .= qq{<center>};
$out .= qq{<button type="button" onClick="get_lists_projects('sort_by_projects')" class="btn">SORT by names</button> };
$out .= qq{<button type="button" onClick="get_lists_projects('sort_by_captures')" class="btn">SORT by captures</button> };
$out .= qq{<button type="button" onClick="ckeck_list_all_projects()" class="btn btn-success">Select ALL</button> };
$out .= qq{<button type="button" onClick="unckeck_list_all_projects()" class="btn btn-danger">Unselect ALL</button> };
$out .= qq{</center><br>};
$out .=  $cgi->end_div();
$out .=  $cgi->start_div();
$out .= "<table class='table table-striped' style='font-size:13px;'>";
$out .= "<thead>";
$out .= $cgi->start_Tr({style=>"background-color:#E9DEFF;font-size:10px"});
#$out .=  $cgi->th("<center>#</center>");
$out .=  $cgi->th("<center>Select</center>");
$out .=  $cgi->th("<center>Name</center>");
$out .=  $cgi->th("<center>Description</center>");
$out .=  $cgi->th("<center>Capture(s)</center>");
$out .=  $cgi->th("<center>Nb Patients</center>");
$out .= $cgi->end_Tr();
$out .= "</thead>";
$out .= "<tbody>";
my ($h_p, $h_c);
my $hCapturesNames;
my $hCapturesNamesProject;
my @lProjectNames = reverse sort keys %{$hProjects};
my $i_p = 0;
foreach my $project_name (@lProjectNames) {
	$i_p++;
	my $description = $hProjects->{$project_name}->{description};
	my $id = $hProjects->{$project_name}->{id};

	my $nb_pat = 0;
	my $patients = 	$query_init->getPatients($id);
	my %captures_id;
	foreach my $p (@$patients){
		$captures_id{$p->{capture_id}} ++;
		$nb_pat++;
	}
	my $hCaptures;
	foreach my $cid (keys %captures_id){
		my $capt =  $query_init->getCaptureInfos($cid);
		my $validation_db = $capt->{'validation_db'};
		$hCapturesNamesProject->{$capt->{name}}->{$project_name} = undef;
		$hCapturesNames->{$capt->{name}}++ if ($validation_db and not $validation_db eq '');
		$hCaptures->{lc($capt->{name})}++;
	}
	
	my $captures = join("<br>", sort keys %$hCaptures);
	my $button_id = "b_proj_$project_name";
	my $button_html = qq{<input type="checkbox" dojoType="dijit.form.CheckBox" id="$button_id" checked/>};
	
	my $captures_id = join(",", sort keys %$hCaptures);
	$h_c->{$captures_id}->{$project_name} = undef;
	
	$h_p->{$project_name} = '';
	$h_p->{$project_name} .=  $cgi->start_Tr();
	#$h_p->{$project_name} .=  $cgi->td("<center>$i_p</center>");
	$h_p->{$project_name} .=  $cgi->td("<center>$button_html</center>");
	$h_p->{$project_name} .=  $cgi->td("<center>$project_name</center>");
	$h_p->{$project_name} .=  $cgi->td("<center>$description</center>");
	$h_p->{$project_name} .=  $cgi->td("<center>$captures</center>");
	$h_p->{$project_name} .=  $cgi->td("<center>$nb_pat</center>");
	$h_p->{$project_name} .=  $cgi->end_Tr();
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
$out .= "</table>";


my $hCapturesNames_for_sort;
foreach my $c_name (keys %$hCapturesNames) {
	$hCapturesNames_for_sort->{lc($c_name)} = $c_name;
}
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



$hRes->{html_projects} = $out;
$hRes->{list_projects}= \@lProjectNames;
$hRes->{list_captures}= \@lCapturesNames;
#$hRes->{list_panels}= \@lCapturesNames;
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
