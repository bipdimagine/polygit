#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/kyoto";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/util_filter";
use GBuffer;
use Data::Dumper; 
use JSON;


my $hashRes;
my $cgi = new CGI();
my $user = $cgi->param('user');
my $pwd    = $cgi->param('pwd');
my $origin_project = $cgi->param('origin_project');
my $force_annot_version = $cgi->param('force_annot_version');
my $defidiag = $cgi->param('defidiag');
my $buffer = new GBuffer;
my $query = $buffer->getQuery();

my $lp2 = $query->listAllProjects();
my $lp = $lp2;
if ($defidiag) {
	@$lp = grep{$_->[3] eq "defidiag"} @$lp2;
}
my %projects;
my $nb;

foreach my $p (@$lp){
	my $capture_name = $p->[1];
	$hashRes->{global}->{patients}->{$capture_name} ++;
	$hashRes->{global}->{patients}->{all} ++;
	$hashRes->{global}->{patients}->{target}++ if $capture_name ne "exome" && $capture_name ne "ciliome" && $capture_name ne "genome" && $capture_name ne "rnaseq";
	$projects{$p->[0]}++;
} 
my @lAllProj = keys %projects;
$hashRes->{global}->{projects}->{all} = scalar(@lAllProj);

 $lp2 = $query->listAllProjects($user,$pwd);
 $lp = $lp2;
if ($defidiag) {
	@$lp = grep{$_->[3] eq "defidiag"} @$lp2;
}
 %projects = ();
foreach my $p (@$lp){
	my $capture_name = $p->[1];
	$hashRes->{user}->{patients}->{$capture_name} ++;
	$hashRes->{user}->{patients}->{target}++ if $capture_name ne "exome" && $capture_name ne "ciliome" && $capture_name ne "genome" && $capture_name ne "rnaseq";
	$hashRes->{user}->{patients}->{all} ++;
	$projects{$p->[0]}++;
}


my $last_genecode = $buffer->getQuery->getMaxGencodeVersion();
my $last_annot = $buffer->getQuery()->getMaxPublicDatabaseVersion();
my $last_version_annot = $last_genecode.'.'.$last_annot;
$hashRes->{global}->{projects}->{last_annot_version} = $last_version_annot;
if ($origin_project) {
	my $projectTmp = $buffer->newProject(-name => $origin_project);
	$hashRes->{global}->{projects}->{annot_version} = $projectTmp->annotation_version();
}
elsif ($force_annot_version) { $hashRes->{global}->{projects}->{annot_version} = $force_annot_version; }
else { $hashRes->{global}->{projects}->{annot_version} = $last_version_annot; }

foreach my $cp  (keys %{$hashRes->{global}->{patients}} ){
	$hashRes->{user}->{patients}->{$cp} = 0 unless exists $hashRes->{user}->{patients}->{$cp};
}
$hashRes->{user}->{projects}->{all} = scalar(keys %projects);

my @lItems;
my $project = $buffer->newProject(-name => $lAllProj[0]);
my $hGenecode = $buffer->getQuery->getHashAllGenesReleasesAnntotations();

foreach my $genecode_id (sort {$a <=> $b} keys %{$hGenecode}) {
	next if $genecode_id > $last_genecode;
	my $genecode_version = $hGenecode->{$genecode_id}->{name};
	next if ($genecode_version < 19);
	my $h = undef;
	my $version_hgmd = $buffer->get_public_data_version("hgmd",$last_annot);
	my $version_clinvar = $buffer->get_public_data_version("clinvar",$last_annot);
	my $version_cosmic = $buffer->get_public_data_version("cosmic",$last_annot);;
	
	$h->{id} = $genecode_version.'.'.$last_annot;
	$h->{value}  = '<b>Hgmd </b>'.$version_hgmd if ($version_hgmd);
	$h->{value} .= '<b> - Clinvar </b>'.$version_clinvar if ($version_clinvar);
	$h->{value} .= '<b> - Cosmic </b>'.$version_cosmic if ($version_cosmic);
	push(@lItems, $h);
}
$hashRes->{global}->{projects}->{all_annot_version} = \@lItems;

#my @lCaptures = sort @{$query->listAllCaptureAnalyse()};
#foreach my $captName (sort @lCaptures) {
#	my $nb = scalar(@{$query->listAllPatientsIdByCapture($captName)});
#	$hashRes->{patients}->{$captName} = $nb;
#	$hashRes->{patients}->{all} += $nb;
#}
#my @lProjects;
#foreach my $captName (sort @lCaptures) {
#	my @lThisProjects = @{$query->listAllProjectsNameByCapture($captName)};
#	my $nb = scalar(@lThisProjects);
#	push (@lProjects, @lThisProjects);
#	$hashRes->{projects}->{$captName} = $nb;
#	$hashRes->{projects}->{all} += $nb;
#}

printJson($hashRes);


sub printJson {
	my ($hashRes) = @_;
	my $cgi = new CGI;
	print $cgi->header('text/json-comment-filtered');
	print encode_json $hashRes;
	exit(0);
}