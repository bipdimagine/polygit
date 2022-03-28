#!/usr/bin/perl

use strict;
use CGI qw/:standard :html3/;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/";
use lib "$Bin/../../GenBo/lib/obj-nodb/";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use GBuffer; 
use Data::Dumper;
use Getopt::Long;
use Parallel::ForkManager;
use Tabix;
use JSON;
use export_data;

my $limitRes = 1000000000;
my @listSnp;
my @listHashRes;
my $hashRes;
my $buffer = GBuffer->new();
my $query = $buffer->getQuery();
my $cgi = new CGI;
my $varId = $cgi->param('input');
my $user = $cgi->param('login');
my $pass = $cgi->param('pwd');
my $exceptProjName = $cgi->param('thisproject');

unless ($varId) {
	my $hash;
	$hash->{'varId'} = "No result...";
	push(@listHashRes, $hash);
	printJson(\@listHashRes);
}

my ($hProjAuthorized);
if ($user and $pass) {
	foreach my $hash (@{$query->getProjectListForUser($user, $pass)}) { $hProjAuthorized->{$hash->{'name'}} = undef; }
}
else {
	foreach my $project_name (@{$buffer->getListAllProjectName()}) { $hProjAuthorized->{$project_name} = undef; }
}
my @lAuthProj = keys(%$hProjAuthorized);
my $projectTmp = $buffer->newProject(-name => $lAuthProj[0]);

if ($varId =~ /1kg_/) { $varId =~ s/1kg_//; }
if ($varId =~ /evs_/) { $varId =~ s/evs_//; }
if (($varId =~ /rs/) or (($varId =~ /TMP_ESP/))) {
	my $rsName = $varId;
	$varId = $projectTmp->convertRsNameToVarId($rsName);
}
if ($varId =~ /\//) {
	my @lTmpVarFields = split('_', $varId);
	my @lVarAll = split('/', $lTmpVarFields[-1]);
	foreach my $varAll (@lVarAll) {
		my $thisVarId = $lTmpVarFields[0].'_'.$lTmpVarFields[1].'_'.$lTmpVarFields[2].'_'.$varAll;
		push(@listSnp, $thisVarId);
	}
}	
else { push(@listSnp, $varId); }

if (scalar(@listSnp) > 0) {
	foreach my $varId (@listSnp) {
		last if (scalar(keys(%$hashRes)) >= $limitRes);
		my $hash = $projectTmp->getDejaVuInfos($varId);
		if ($hash) { 
			my $hPatiens;
			my $var = $projectTmp->_newVariant($varId);
			foreach my $projName (sort keys %$hash) {
				#next if ($projName eq 'NGS2014_0387');
				if ($exceptProjName) {
					next if ($projName eq $exceptProjName);
				} 
				my $thisProjectTmp = $buffer->newProject(-name => $projName);
				my $h_heho;
				foreach my $string (split(';', $hash->{$projName}->{string})) {
					my ($patName, $heho) = split(':', $string);
					$h_heho->{$patName} = $heho;
				}
				foreach my $patName (split(';', $hash->{$projName}->{patients})) {
					my @lTmp = split('_', $varId);
					my @traces_id;
					my $chr_name = $lTmp[0];
					my $var_all = $lTmp[-1];
					push (@traces_id, $chr_name);
					my $thisProjName = $projName;
					my $thisPatName  = $patName;
					my $other_project;
					if (exists($hProjAuthorized->{$projName})) {
						$other_project = 'n';
					}
					else {
						$other_project = 'y';
					}
					my $id = $thisProjName.'_'.$thisPatName;
					$hashRes->{$id}->{'other_project'} = $other_project;
					$hashRes->{$id}->{'traces_id'} = \@traces_id;
					$hashRes->{$id}->{'traces_name'} = \@traces_id;
					$hashRes->{$id}->{'poly_id'} = $varId;
					$hashRes->{$id}->{'variation_id'} = $varId;
					$hashRes->{$id}->{'project'} = $thisProjName;
					$hashRes->{$id}->{'name'} = $thisPatName;
					$hashRes->{$id}->{'id'} = $thisProjName.'_'.$thisPatName;
					$hashRes->{$id}->{'var_allele'} = $var_all;
					$hashRes->{$id}->{'hohe'} = 'ho' if ($h_heho->{$patName} == 1);
					$hashRes->{$id}->{'hohe'} = 'he' if ($h_heho->{$patName} == 2);
				}
				$thisProjectTmp = undef;
			}
		}
	}
}

if (scalar(keys(%$hashRes)) == 0) { 
	my $hash;
	$hash->{'id'} = "No result...";
	push(@listHashRes, $hash);
}
else {
	foreach my $id (sort(keys(%$hashRes))) { push(@listHashRes, $hashRes->{$id}); }
}
printJson(\@listHashRes);



sub printJson {
	my ($listHash) = @_;
	my $hash;
	$hash->{'identifier'} = 'id';
	$hash->{'label'} = 'id';
	$hash->{'items'} = $listHash;
	print $cgi->header('text/json-comment-filtered');
	print encode_json $hash;
	exit(0);
}