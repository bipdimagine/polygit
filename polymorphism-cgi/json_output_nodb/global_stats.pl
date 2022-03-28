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

my $buffer = new GBuffer;
my $query = $buffer->getQuery();

my @lCaptures = sort @{$query->listAllCaptureAnalyse()};
foreach my $captName (sort @lCaptures) {
	my $nb = scalar(@{$query->listAllPatientsIdByCapture($captName)});
	$hashRes->{'nb_patients_capture|'.$captName} = $nb;
}
my @lProjects;
foreach my $captName (sort @lCaptures) {
	my @lThisProjects = @{$query->listAllProjectsNameByCapture($captName)};
	my $nb = scalar(@lThisProjects);
	push (@lProjects, @lThisProjects);
	$hashRes->{'nb_projects_capture|'.$captName} = $nb;
	$hashRes->{'nb_projects'} += $nb;
}
my @lYears;
my $existsHiseq;
my $lastYear = 2000;
foreach my $projName (@lProjects) {
	next if ($projName =~ /TEST/);
	$existsHiseq = 1 if ($projName eq 'HISEQ');
	my @lTmp = split('_', $projName);
	$lTmp[0] =~ s/[A-Z]+//;
	push(@lYears, $lTmp[0]);
	$hashRes->{'nb_projects_year|'.$lTmp[0]}++;
	$lastYear = $lTmp[0] if (int($lTmp[0]) > $lastYear);
}
$hashRes->{'nb_projects_year|'.$lastYear}++ if ($existsHiseq);

$hashRes->{'nb_patients'} = scalar(@{$buffer->listAllPatientsId()});
#foreach my $year (sort(@lYears)) {
#	$hashRes->{'nb_patients_year|'.$year} = scalar(@{$query->listAllPatientsIdByYear($year)});
#}
my $nbWithfDate = 0; 
my $nbUndefDate = 0; 
my $hashDate = $query->listAllPatientsIdByYear();
foreach my $patId (keys %$hashDate) {
	my @lTmp = split('-', $hashDate->{$patId}->{'creation_date'});
	if ($lTmp[0] eq '0000') {
		$nbUndefDate++;
		my @lTmp2 = split('_', $hashDate->{$patId}->{'project_name'});
		$lTmp2[0] =~ s///;
		$lTmp2[0] =~ s/[A-Z]+//;
		$hashRes->{'nb_patients_year|'.$lTmp2[0]}++;
	}
	else {
		$nbWithfDate++;
		$hashRes->{'nb_patients_year|'.$lTmp[0]}++;
	}
}
warn 'Date defined: '.$nbWithfDate;
warn 'No Date defined: '.$nbUndefDate;
printJson($hashRes);


sub printJson {
	my ($hashRes) = @_;
	my $cgi = new CGI;
	print $cgi->header('text/json-comment-filtered');
	print encode_json $hashRes;
	exit(0);
}