#!/usr/bin/perl
$|=1;
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

my @listHashRes;
my $hashRes;
my $buffer = GBuffer->new();
my $query = $buffer->getQuery();
my $cgi = new CGI;
my $varId = $cgi->param('input');
my $user = $cgi->param('login');
my $pass = $cgi->param('pwd');
my $build_use = $cgi->param('build');
my $patient_name = $cgi->param('patient');
my $project_name = $cgi->param('project');
my $export_html_bootstrap = $cgi->param('export_html_bootstrap');
my $in_this_run = $cgi->param('in_this_run');

$build_use = 'HG19' unless ($build_use);

$build_use = 'HG38';

print $cgi->header('text/json-comment-filtered');
print "{\"progress\":\".";
my $fam;
my $p;
my $h_this_run_patients;
my $solo;
my $solo_fam = {};
my $var;

my $projectTmp;
unless ($varId) {
	my $hash;
	$hash->{'varId'} = "No result...";
	push(@listHashRes, $hash);
	printJson(\@listHashRes);
}
elsif ($varId =~ /;/) {
	$varId = split(';', $varId)->[0];
}
if ($project_name and $project_name =~ /NGS/){
	$projectTmp = $buffer->newProjectCache(-name => $project_name);
	$var = $projectTmp->_newVariant($varId);
	if($projectTmp->isDefidiagSolo){
		$solo =1 ;
		foreach my $patient (@{$var->getPatients}){
			$solo_fam->{$patient->getFamily->name()} ++;
		}
	}
	if ($in_this_run) {
		my $in_this_run_patients =  $projectTmp->in_this_run_patients();
		foreach my $k (keys %{$in_this_run_patients}) {
			next unless ($k =~ /NGS20/);
			my $run_proj_name = $k;
			foreach my $run_pat_id (keys %{$in_this_run_patients->{$run_proj_name}}) {
				my $run_pat_name = $in_this_run_patients->{$run_proj_name}->{$run_pat_id};
				$h_this_run_patients->{$run_proj_name}->{$run_pat_name} = $run_pat_id;
			}
		}
	}
}
else {
	my $project_init_name = $buffer->get_random_project_name_with_this_annotations_and_genecode();
	$projectTmp = $buffer->newProject(-name => $project_init_name);
	$var = $projectTmp->_newVariant($varId);
}

unless ($build_use) {
	die();
	my $hash;
	$hash->{'varId'} = "Please give build option...";
	push(@listHashRes, $hash);
	printJson(\@listHashRes);
}

my $hProjAuthorized;
foreach my $hash (@{$query->getProjectListForUser($user, $pass)}) { $hProjAuthorized->{$hash->{'name'}} = undef; }

my $id = 0;
my $hDejaVuGlobal = $var->dejavu_hash_projects_patients();

foreach my $projName (sort keys %$hDejaVuGlobal) {
	next if ($in_this_run and not exists $h_this_run_patients->{$projName});
	
	print '.';
	my $thisProject = $buffer->newProject(-name => $projName);
	delete $thisProject->{version};
	#delete $thisProject->{genomeFai};
	#delete $thisProject->{genomeFasta};
	#delete $thisProject->{dirGenome};
	#delete $thisProject->{genome_version};
	#warn $build_use;
	$thisProject->version($build_use);
	my @lPheno;
	foreach my $pheno_obj (@{$thisProject->getPhenotypes()}) {
		push(@lPheno, $pheno_obj->name());
	}
	my @lMails;
	foreach my $hashOwners (@{$query->getOwnerProject_byName($projName)}) { push(@lMails, $hashOwners->{'email'});  }
	
	if (not exists $hDejaVuGlobal->{$projName}->{patients}) {
		my $h = $var->infos_dejavu_parquet($thisProject);
		$hDejaVuGlobal->{$projName}->{patients} = $h->{patients_infos} if $h and exists $h->{patients_infos};
	}
	
	foreach my $patName (sort keys %{$hDejaVuGlobal->{$projName}->{patients}}) {
		next if ($in_this_run and not exists $h_this_run_patients->{$projName}->{$patName});
		my $pp = $thisProject->getPatient($patName);
		next if exists $solo_fam->{$pp->getFamily()->name} ;
		
		#die();
		my $hash;
		$hash->{id} = $id;
		$hash->{project} = $projName;
		if (@lPheno and not $export_html_bootstrap) {
			$hash->{project} .= ';'.join(' ', sort @lPheno);
		}
		$hash->{name} = $patName;
		$hash->{contacts} = join(' - ', sort @lMails);
		if (exists $hProjAuthorized->{$projName}) {
			$hash->{implicated} = 'X';
		}
		else {
			$hash->{implicated} = '';
		}
		$hash->{ho_he} = 'todo';
#		$hash->{ho_he} = 'he' if ($ho_he eq '2');
#		$hash->{ho_he} = 'ho' if ($ho_he eq '1');

		$hash->{variation_id} = $varId;
		$hash->{poly_id} = $varId;
		my ($chr, $pos, $ref, $mut) = split('_', $varId);
		$hash->{traces_id} = ($chr);
		$hash->{traces_name} = ($chr);
		my $locus = 'chr'.$var->getChromosome->id().':'.$var->start().'-'.$var->end();
		
		my $patient;
		eval { $patient = $thisProject->getPatient($patName); };
		if ($@) {
			$hash->{igv} = 'Not Found';
		}
		else {
			my @lTmpBam = split('ngs', $thisProject->getPatient($patName)->getBamFiles->[0]); 
			my $bam_file = '/NGS'.$lTmpBam[1];
			my $cmd_igv = "IGV;launch_web_igv_js('$projName','$patName','$bam_file','$locus')";
			$hash->{igv} = $cmd_igv;
			my $run = $patient->getRun();
			my @lDate = split(' ',$run->date());
			$hash->{run_id} = $run->id(). ' '.$lDate[0];
			$hash->{machine} = $run->machine();
			my @lCaptures;
			foreach my $capture (@{$patient->getCaptures()}) { push(@lCaptures, $capture->name()); }
			$hash->{captures} = join(', ', @lCaptures);
			$hash->{family} = $patient->family();
			$hash->{sex} = $patient->sex();
			$hash->{status} = $patient->status();
			$hash->{ho_he} = $hDejaVuGlobal->{$projName}->{patients}->{$patName}->{heho} if exists $hDejaVuGlobal->{$projName}->{patients}->{$patName}->{heho};
			$hash->{model} = $hDejaVuGlobal->{$projName}->{patients}->{$patName}->{model} if exists $hDejaVuGlobal->{$projName}->{patients}->{$patName}->{model};
			$hash->{dp} = $hDejaVuGlobal->{$projName}->{patients}->{$patName}->{dp} if exists $hDejaVuGlobal->{$projName}->{patients}->{$patName}->{dp};
			$hash->{ratio} = $hDejaVuGlobal->{$projName}->{patients}->{$patName}->{ratio} if exists $hDejaVuGlobal->{$projName}->{patients}->{$patName}->{ratio};
			my $sex = $patient->sex();
			my $status = $patient->status();
			if ($patient->isFather()) { $hash->{sex_status} = 'F'.$status; }
			elsif ($patient->isMother()) { $hash->{sex_status} = 'M'.$status; }
			else { $hash->{sex_status} = 'C'.$sex.$status; }
		}
		push(@listHashRes, $hash);
		$id++;
	}
	$thisProject = undef;
}


if (scalar(@listHashRes) == 0) { 
	my $hash;
	$hash->{'varId'} = "No result...";
	push(@listHashRes, $hash);
}
printHtml_bootstrap(\@listHashRes) if ($export_html_bootstrap);
printJson(\@listHashRes);



sub printHtml_bootstrap {
	my ($listHash) = @_;
	my $cgi = new CGI();
	my $hash;
	my @lHeaders = ('Acc', 'Project', 'Family', 'Patient', 'Sex/Status', 'He/Ho', "Ratio", "Model", 'Run Id', 'Machine', 'Capture', 'Contacts', 'IGV');
	my @lCat = ('implicated', 'project', 'family', 'name', 'sex_status', 'ho_he', 'ratio', 'model', 'run_id', 'machine', 'captures', 'contacts', 'igv');
	my $out = "<table class='table table-bordered' style='font-size:11px;'>";
	$out .= "<thead>";
	$out .= $cgi->start_Tr();
	foreach my $header (@lHeaders){
		$out .=  $cgi->th("<center><b><u>".ucfirst($header)."</center></b></u>");
	}
	$out .= $cgi->end_Tr();
	$out .= "</thead>";
	$out .= "<tbody>";
	my $class_td;
	$class_td->{style} = "max-height:50px;overflow-y:auto;";
	foreach my $h (@$listHash) {
		$out .= $cgi->start_Tr();
		foreach my $cat (@lCat) {
			my $value = $h->{$cat};
			next unless $value;
			if ($cat eq 'implicated') { $value = format_implicated_icon($value); }
			if ($cat eq 'project') { $value = format_project_name_phenotype($value); }
			if ($cat eq 'sex_status') { $value = format_patient_icon($value); }
			if ($cat eq 'contacts') { $value = format_mails_list($value); }
			if ($cat eq 'igv') { $value = format_igv_icon($value); }
			$out .= $cgi->td($class_td, "<center>".$value."</center>");
		}
		$out .= $cgi->end_Tr();
	}
	$out .= "</tbody>";
	$out .= "</table>";
	$hash->{'html_table'} = $out;
	print ".\",";
	my $json_encode = encode_json $hash;
	$json_encode =~ s/{//;
	print $json_encode;
	exit(0);
}

sub format_project_name_phenotype {
	my $value = shift;
	my $pheno_name;
	my $buffer_pheno = GBuffer->new();
	my $project_pheno = $buffer->newProject(-name => $value);
	print ".";
	my @lPheno;
	foreach my $pheno_obj (@{$project_pheno->getPhenotypes()}) {
		push(@lPheno, $pheno_obj->name());
	}
	$project_pheno = undef;
	$buffer_pheno = undef;
	if (@lPheno) {
		my $pheno_name = join('<br>', sort @lPheno);
		$value .= qq{<br><center><span style='color:red;'>$pheno_name</span></center>};
	}
	return $value;
}

sub format_mails_list {
	my $value = shift;
	$value =~ s/ - /<br>/g;
	return $value;
}

sub format_implicated_icon {
	my $value = shift;
	if ($value eq 'X') { return qq{<span class="glyphicon glyphicon-ok" aria-hidden="true" style='font-size:12px;color:green;'></span>}; }
	return qq{<span class="glyphicon glyphicon-remove" aria-hidden="true" style='font-size:12px;color:red;'></span>};
}

sub format_igv_icon {
	my $value = shift;
	my ($name, $cmd) = split(';', $value);
	my $b = qq{<button onClick="$cmd">IGV</button>};
	return $b;
}

sub format_patient_icon {
	my $value = shift;
	if ($value eq "C21") { return "<center><img src='/icons/Polyicons/baby-girl-s.png'></center>"; }
	if ($value eq "C22") { return "<center><img src='/icons/Polyicons/baby-girl-d.png'></center>"; }
	if ($value eq "C11") { return "<center><img src='/icons/Polyicons/baby-boy-s.png'></center>"; }
	if ($value eq "C12") { return "<center><img src='/icons/Polyicons/baby-boy-d.png'></center>"; }
	if ($value eq "F1") { return "<center><img src='/icons/Polyicons/male-s.png'></center>"; }
	if ($value eq "F2") { return "<center><img src='/icons/Polyicons/male-d.png'></center>"; }
	if ($value eq "M1") { return "<center><img src='/icons/Polyicons/female-s.png'></center>"; }
	if ($value eq "M2") { return "<center><img src='/icons/Polyicons/female-d.png'></center>"; }
	if ($value eq "F") { return "<center><img src='/icons/Polyicons/male.png'></center>"; }
	if ($value eq "M") { return "<center><img src='/icons/Polyicons/female.png'></center>"; }
	if ($value eq "1") { return "<center><img src='/icons/Polyicons/male.png'></center>"; }
	if ($value eq "2") { return "<center><img src='/icons/Polyicons/female.png'></center>"; }
	return "<center><img src='/icons/Polyicons/12-em-check.png'></center>";
}

sub printJson {
	
	my ($listHash) = @_;
	my $hash;
	$hash->{'identifier'} = 'id';
	$hash->{'label'} = 'id';
	$hash->{'items'} = $listHash;
	my $json_encode = encode_json $hash;
	print ".\",";
	$json_encode =~ s/{//;
	print $json_encode;
	exit(0);
}
