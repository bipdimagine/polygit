#!/usr/bin/perl
$|=1;

use strict;
use CGI qw/:standard :html3/;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/";
use lib "$Bin/../../GenBo/lib/obj-nodb/";
require "$Bin/../../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use GBuffer; 
use Data::Dumper;
use Getopt::Long;
use Parallel::ForkManager;
use Tabix;
use JSON;
use export_data;
use Math::Combinatorics;
use Spreadsheet::WriteExcel;
use POSIX qw(strftime);
  use Set::Intersection;
 use  Algorithm::Combinatorics qw(combinations);
my @listSnp;
my @listHashRes;
my $hashRes;
my $buffer = GBuffer->new();
my $query = $buffer->getQuery();
my $cgi = new CGI;
my $input = $cgi->param('input');
my $user = $cgi->param('login');
my $pass = $cgi->param('pwd');
my $thisproject = $cgi->param('thisproject');
my $onlyho = $cgi->param('onlyho');
my $show_large_deletion = $cgi->param('show_large_del');
my $build_use = $cgi->param('build');
my $export_xls = $cgi->param('xls');
my $outfile = $cgi->param('outfile');
my $filters_consequence =  $cgi->param('f_cons');
my $filters_consequence_1 =  $cgi->param('f_cons1');
my $filters_consequence_2 =  $cgi->param('f_cons2');
my $filters_var_id =  $cgi->param('f_id');
my $origin_project = $cgi->param('origin_project');
my $annot_version = $cgi->param('force_annot_version');


if (not $input or $onlyho) {
	my $hash;
	$hash->{'id'} = "No result...";
	push(@listHashRes, $hash);
	printJson(\@listHashRes);
}
unless ($build_use) {
	my $hash = getHashInError();
	push(@listHashRes, $hash);
	printJson(\@listHashRes);
}

my $hProjAuthorized;
foreach my $hash (@{$query->getProjectListForUser($user, $pass)}) { $hProjAuthorized->{$hash->{'name'}} = undef; }
my @lAuthProj = keys(%$hProjAuthorized);
my $projectTmp;
if ($origin_project) {
	$projectTmp = $buffer->newProject(-name => $origin_project);
}
else {
	$projectTmp = $buffer->newProject(-name => $lAuthProj[0]);
	if ($annot_version) { $projectTmp->changeAnnotationVersion( $annot_version, 1 ); }
	else { $projectTmp->changeAnnotationVersion( $projectTmp->buffer->getQuery->getCurrentGenomeProjectReleasesAnntotations(), 1 ); }
}

$projectTmp->version($build_use);
my $tabix = $projectTmp->getSoftware('tabix');
my $type;


if ( $input =~ 'id:' ) {
	$input =~ s/id://;
	$input =~ s/ //;
	if ( $input =~ /:/ ) {
		$input =~ s/:/_/;
		$input =~ s/-/_/;
		$input = 'region:'.uc($input);
	}
	elsif ( $input =~ /-/ ) {
#		my @lTmp = split('-', $input);
#		$input = 'geneId:'.uc($lTmp[-1]);
		$input = 'geneId:'.uc($input);
	}
	elsif ( $input =~ /[Ee][Nn][Ss][Tt]/ ) { $input = 'transId:'.uc($input); }
	elsif ( $input =~ /[Ee][Nn][Ss][Gg]/ ) { $input = 'geneId:'.uc($input); }
	elsif ( $input =~ /,/ ) { $input = 'varId:'.$input; }
	elsif ( $input =~ /_/ ) { $input = 'varId:'.$input; }
	elsif ( $input =~ /[Rr][Ss][0-9]+/ ) { $input = 'varId:'.lc($input); }
	else { $input = 'geneId:'.uc($input); }
}

my $originInput = $input;

if ( $input =~ 'region:' ) {
	my $hash;
	$hash->{'id'} = "Analyse only available for Transcript or Gene input...";
	push(@listHashRes, $hash);
	printJson(\@listHashRes);
}

if ( $input =~ 'varId:' ) {
	$type = 'var';
	my $varId = $input;
	$varId =~ s/varId://;
	my @lIds = split(',', $varId);
	foreach my $varId (@lIds) { 
		if ($varId =~ /1kg_/) { $varId =~ s/1kg_//; }
		if ($varId =~ /evs_/) { $varId =~ s/evs_//; }
		if (($varId =~ /rs/) or (($varId =~ /TMP_ESP/))) {
			my $rsName = $varId;
			$varId = $projectTmp->convertRsNameToVarId($rsName);
			unless ($varId) {
				my $hash;
				$hash->{'id'} = "...";
				$hash->{'implicated'} = "...";
				$hash->{'var_1'} = "ERROR: $rsName variant doesn't exist...";
				$hash->{'var_2'} = "...";
				$hash->{'consequence_1'} = "...";
				$hash->{'consequence_2'} = "...";
				$hash->{'nbProj'} = "...";
				$hash->{'nbPat'} = "...";
				push(@listHashRes, $hash);
				printJson(\@listHashRes);
				exit(0);
			}
			if ($varId =~ /,/) {
				my $hash;
				$hash->{'id'} = "...";
				$hash->{'implicated'} = "...";
				$hash->{'var_1'} = "ERROR: multiple varIds found for $rsName ($varId). Please, choose one of them !";
				$hash->{'var_2'} = "...";
				$hash->{'consequence_1'} = "...";
				$hash->{'consequence_2'} = "...";
				$hash->{'nbProj'} = "...";
				$hash->{'nbPat'} = "...";
				push(@listHashRes, $hash);
				printJson(\@listHashRes);
				exit(0);
			}
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
	}
}

if ( $input =~ 'transId:' ) {
	my $transId = $input;
	$transId =~ s/transId://;
	if ($transId =~ /ENST/) {
		my $trans;
		eval { $trans = $projectTmp->newTranscript($transId); };
		if($@) {
			my $hash;
			$hash->{'id'} = "ERROR: $transId transcript doesn't exist...";
			push(@listHashRes, $hash);
			printJson(\@listHashRes);
		}
		my $chr = $trans->getChromosomes()->[0]->name();
		my $start = int($trans->start());
		my $end = int($trans->end());
		$input = 'region:'.$chr.'_'.$start.'_'.$end;
	}
}

if ( $input =~ 'geneId:' ) {
	my @lOutput;
	my $geneId = $input;
	$geneId =~ s/geneId://;
	my @lGenesId;
	@lGenesId = split(',', $geneId);
	my @lGenesId2;
	foreach my $geneId (@lGenesId) {
		my $gene;
		eval { $gene = $projectTmp->newGene($geneId);  };
		if($@) {
			if (scalar(@lGenesId) == 1) {
				my $hash;
				$hash->{'id'} = "ERROR: $geneId gene doesn't exist...";
				push(@listHashRes, $hash);
				printJson(\@listHashRes);
			}
		}
		else {
			my $chr = $gene->getChromosomes()->[0]->name();
			my $start = int($gene->start());
			my $end = int($gene->end());
			push(@lOutput, 'region:'.$chr.'_'.$start.'_'.$end);
			push(@lGenesId2, $gene->id());
		}
	}
	$input = join(',', @lOutput);
	$originInput = join(',', @lGenesId2);
}

my $hVar;
my @lTmp = split(',', $input);
$input = undef;
unless ($type eq 'var)') {
	foreach my $input (sort(@lTmp)) {
		if ($input =~ 'region:') {
			my $posId = $input;
			$posId =~ s/region://;
			@listSnp = @{ $projectTmp->getDejaVuIdsFromInterval($posId) };
		}
	}
}

unless ($export_xls) {
	$projectTmp->cgi_object(1);
	print $cgi->header('text/json-comment-filtered');
	print "{\"progress\":\".";
}

my $nbSnpsSearching = scalar(@listSnp);
my $hRes;
if (scalar(@listSnp) > 1) {
	foreach my $varId (@listSnp) {
		$projectTmp->print_dot(1);
		my $hDejaVuId = $projectTmp->getDejaVuInfos($varId);
		foreach my $projName (keys %$hDejaVuId) {
			next if ($hDejaVuId->{$projName}->{he} == 0);
			_getInfosForVar($originInput, $varId);
			foreach my $patString (split(';', $hDejaVuId->{$projName}->{string})) {
				my ($patName, $heho) = split(':', $patString);
				if ($heho eq '2') {
					my $projPatName = $projName.'|'.$patName;
					$hashRes->{$projPatName}->{'project_patient'} = $projPatName;
					push(@{$hashRes->{$projPatName}->{'var_he'}}, $varId);
				}
			}
		}
	}
}
my $lHashHeComp = getHeComposite1($hashRes);
if (scalar(@$lHashHeComp) == 0) { 
	my $hash;
	$hash->{'id'} = "...";
	$hash->{'implicated'} = "...";
	$hash->{'var_1'} = "No result...";
	$hash->{'var_2'} = "...";
	$hash->{'consequence_1'} = "...";
	$hash->{'consequence_2'} = "...";
	$hash->{'nbProj'} = "...";
	$hash->{'nbPat'} = "...";
	push(@listHashRes, $hash);
}
else { @listHashRes = @$lHashHeComp; }
if ($export_xls) { exportXls($originInput, \@listHashRes); }
else { printJson(\@listHashRes); }



sub _isASnp {
	my ($varId) = shift;
	my $var = $projectTmp->_newVariant($varId);
	return 1 if ($var->isVariation());
	return undef;
}

sub _getInfosForVar {
	my ($input, $varId) = @_;
	my $var = $projectTmp->_newVariant($varId);
	my @lGenesObj = @{$var->getGenes()};
	my @lTransObj = @{$var->getTranscripts()};
	my $objAnalyse = undef;
	if ($input =~ /geneId:/) {
		my $tmp = $input;
		$tmp =~ s/geneId://;
		foreach my $gene (@lGenesObj) {
			my $name = $gene->name();
			if ($tmp =~ /$name/) { $objAnalyse = $gene; }
			else {
				my $name2 = $gene->external_name();
				if ($tmp =~ /$name2/) { $objAnalyse = $gene; }
			}
		}
		unless($objAnalyse) { return undef; }
	}
	elsif ($input =~ /transId:/) {
		my $tmp = $input;
		$tmp =~ s/transId://;
		foreach my $trans (@lTransObj) { if ($trans->name() eq $tmp) { $objAnalyse = $trans; } }
	}
	eval {
		my @lConsOk;
		foreach my $cons (split(',', lc($var->variationType($objAnalyse)))) {
			my @lCons = split(';', $projectTmp->buffer->config->{ensembl_annotations}->{$cons});
			push(@lConsOk, $lCons[1]);
		}
		$hVar->{$varId}->{'consequence'} = join(', ', @lConsOk);
	};
	if ($@) { $hVar->{$varId}->{'consequence'} = 'Error annotation...'; }
	my $polyphenScore = "-";
	if (($objAnalyse) and ($objAnalyse->isTranscript())) { $polyphenScore = $var->polyphenScore($objAnalyse); }
	eval { $hVar->{$varId}->{'polyphen'} = $var->polyphenStatus($objAnalyse) . "+" . $polyphenScore; };
	if ($@) { $hVar->{$varId}->{'polyphen'} = '4+-'; }
	my $siftScore = "-";
	if (($objAnalyse) and ($objAnalyse->isTranscript())) { eval { $siftScore = $var->siftScore($objAnalyse); } }
	eval { $hVar->{$varId}->{'sift'} = $var->siftStatus($objAnalyse) . "+" . $siftScore; };
	if ($@) { $hVar->{$varId}->{'sift'} = '4+-'; }
	$hVar->{$varId}->{'rsName'} = $var->rs_name() if ($var->rs_name());
	foreach my $trans (@{$var->getTranscripts()}) { $hVar->{$varId}->{'transcript'}->{$trans->id()} = undef; }
	if ($var->getProject->hasHgmdAccess($user)) {
		$hVar->{$varId}->{'hgmd_class'} = '';
		if ($var->hgmd()) {
			my $h;
			$h->{obj} = $var;
			my ($class, $cmd) = update::hgmd($var->getProject,$h,'return_cmd');
			$hVar->{$varId}->{'hgmd_class'} = $class.';'.$cmd;
			if ($var->isNewHgmd()) {
				$hVar->{$varId}->{'hgmd_class'} .= ";New!";
			}
		}
	}
	else {
		$hVar->{$varId}->{'hgmd_class'} = "N.A.";
	}
	$var = undef;
}

sub _getPatProjHeFromHash {
	my ($hash) = shift;
	my @lRes;
	foreach my $proj (keys %$hash) {
		if ($thisproject eq $proj) { next(); }
		foreach my $patName (keys(%{$hash->{$proj}})) { if (int($hash->{$proj}->{$patName}->{'annex'}->{'he'}) > 0) { push(@lRes, $proj.'|'.$patName); }; }
	}
	return (\@lRes);
}

sub getHeComposite {
	my ($hash) = @_;
	my %combine;
	my %vars;
	
	my $apatients = [];
	my %he;
	foreach my $project_patient (keys(%$hash)) {
		next if scalar($hash->{$project_patient}->{'var_he'}) < 2;
	#	my @all =  combine(2, @{$hash->{$project_patient}->{'var_he'}});
		my $iter = combinations($hash->{$project_patient}->{'var_he'},2);
		my $nb1=0;
		while (my $p = $iter->next) {
			$nb1++;
 		}
 	#	die();
		my @vars = sort @{$hash->{$project_patient}->{'var_he'}};
		my $nb=0;
		for (my $i =0; $i<@vars;$i++){
			for (my $j =$i+1; $j<@vars;$j++){
				$he{$vars[$i]."-".$vars[$j]} ++;
				$nb ++;
			}
		}
		warn $nb1." ".$nb;
	}
warn scalar keys 	%he;
die();
	my @vars2;
	my %transcripts;
	
	#
	foreach my $v (keys %vars){
	#	next if scalar(@{$vars{$v}}) <2;
		push(@vars2,$v);
	}
		
		 	my $c = Math::Combinatorics->new( count => 2,       #treated as int
                                             data =>  \@vars2,#arrayref or anonymous array
                                           );
           my $nb;                                
           while (my @t = $c->next_combination()){
           	my @intersection = get_intersection($vars{$t[0]}, $vars{$t[1]});
           	next unless @intersection;
           	my $string = join("-",sort{$a <=> $b} @t);
           	$combine{$string} = \@intersection;
           	$nb ++;
           } 
           die();
}

sub getHeComposite1 {
	my ($hash) = @_;
	#my @lFilters;
	my $hashDetails;
	my $hashHeComp;
	my $hashId;
	foreach my $project_patient (keys(%$hash)) {
		$projectTmp->print_dot(50);
		if (scalar($hash->{$project_patient}->{'var_he'}) >= 2) {
			my @lFields = split('\|', $project_patient);
			my $projName = $lFields[0];
			my $patName = $lFields[1];
			my $myProj;
			if (exists($hProjAuthorized->{$projName})) { $myProj = 'X'; }
			#my @all_combination;
 			#push(@all_combination, combine(2, @{$hash->{$project_patient}->{'var_he'}}));
 			my @vars = sort @{$hash->{$project_patient}->{'var_he'}};
 			for (my $i =0; $i<@vars;$i++){
				for (my $j =$i+1; $j<@vars;$j++){
					my $combHeComp = [$vars[$i],$vars[$j]];
					
 			#foreach my $combHeComp (sort(@all_combination)){
 					my $hasCommonTrans;
 					foreach my $transName (keys(%{$hVar->{$combHeComp->[0]}->{'transcript'}})) {
 						if (exists($hVar->{$combHeComp->[1]}->{'transcript'}->{$transName})) { $hasCommonTrans = '1'; }
 					}
 					next unless $hasCommonTrans;
 					my $couple = join('|', sort(@$combHeComp));
					my @coupleId;
					$hashHeComp->{$couple}->{'nbProj'}->{$projName} = undef;
					$hashHeComp->{$couple}->{'nbPat'}->{$projName}->{$patName} = undef;
					$hashHeComp->{$couple}->{'var_1'} = $combHeComp->[0];
					my $id = $combHeComp->[0];
					if (exists ($hVar->{$combHeComp->[0]}->{'rsName'})) {
						$hashHeComp->{$couple}->{'var_1'} = $hVar->{$combHeComp->[0]}->{'rsName'};
						push(@coupleId, $id.';'.$hVar->{$combHeComp->[0]}->{'rsName'});
					}
					else { push(@coupleId, $id); }
					$hashHeComp->{$couple}->{'consequence_1'} = $hVar->{$combHeComp->[0]}->{'consequence'};
					$hashHeComp->{$couple}->{'hgmd_class_1'} = $hVar->{$combHeComp->[0]}->{'hgmd_class'};
					
					$id = $combHeComp->[1];
					$hashHeComp->{$couple}->{'var_2'} = $combHeComp->[1];
					if (exists ($hVar->{$combHeComp->[1]}->{'rsName'})) {
						$hashHeComp->{$couple}->{'var_2'} = $hVar->{$combHeComp->[1]}->{'rsName'};
						push(@coupleId, $id.';'.$hVar->{$combHeComp->[1]}->{'rsName'});
					}
					else { push(@coupleId, $id); }
					$hashHeComp->{$couple}->{'consequence_2'} = $hVar->{$combHeComp->[1]}->{'consequence'};
					$hashHeComp->{$couple}->{'hgmd_class_2'} = $hVar->{$combHeComp->[1]}->{'hgmd_class'};
					@coupleId = sort(@coupleId);
					$id = join('|', @coupleId);
					$hashHeComp->{$couple}->{'id'} = $id;
					$hashHeComp->{$couple}->{'implicated'} = '';
					if ($myProj) { $hashHeComp->{$couple}->{'implicated'} = 'X'; }
					else { $hashHeComp->{$couple}->{'implicated'} = 'no'; }
					
					my @lCons;
					foreach my $cons (split(', ', $hashHeComp->{$couple}->{'consequence_1'})) {
						push(@lCons, $cons.'.1');
					}
					foreach my $cons (split(', ', $hashHeComp->{$couple}->{'consequence_2'})) {
						push(@lCons, $cons.'.2');
					}
					push(@lCons, $id);
					$hashHeComp->{$couple}->{'allConsequences'} = join('_', @lCons);
					$hashHeComp->{$couple}->{'pass_filter'} = '';
					#my $coupleCons = $hashHeComp->{$couple}->{'consequence_1'} . '_' . $hashHeComp->{$couple}->{'consequence_2'};
					#$hashHeComp->{$couple}->{'allConsequences'} = $coupleCons.'_'.$id;
					
#					if ($hashHeComp->{$couple}->{'id'} =~ /rs929387/ and $hashHeComp->{$couple}->{'id'} =~ /rs846266/) {
#	 					warn $couple;
#	 					warn Dumper $hashHeComp->{$couple};
#	 					die;
#					}
 			}#end for j
 			}#end for j
		}
	}
	my @lHashHeComp;
	foreach my $coupleId (sort(keys(%$hashHeComp))) {
		my @lProjects;
		foreach my $projName (sort(keys(%{$hashHeComp->{$coupleId}->{'nbProj'}}))) {
#			if (exists($hProjAuthorized->{$projName})) {
				push(@lProjects, "$projName");
#			}
#			else {
#				warn "$projName not authorized..."; 
#				push(@lProjects, "X");
#			}
		}
		$hashHeComp->{$coupleId}->{'nbProj'} = join('; ', @lProjects);
		$hashHeComp->{$coupleId}->{'implicated'} = '.';
		my $nbPat = 0;
		my @lPatients;
		foreach my $projName (sort(keys(%{$hashHeComp->{$coupleId}->{'nbPat'}}))) {
			$nbPat += scalar(keys(%{$hashHeComp->{$coupleId}->{'nbPat'}->{$projName}}));
#			if (exists($hProjAuthorized->{$projName})) {
				my $res = $projName.' ['.join(', ', sort(keys(%{$hashHeComp->{$coupleId}->{'nbPat'}->{$projName}}))).']';
				push(@lPatients, $res);
				$hashHeComp->{$coupleId}->{'implicated'} = 'X' if (exists($hProjAuthorized->{$projName}));
#			}
#			else {
#				my @lPat;
#				foreach my $pat (sort(keys(%{$hashHeComp->{$coupleId}->{'nbPat'}->{$projName}}))) { push(@lPat, 'X'); }
#				my $res = 'X ['.join(', ', @lPat).']';
#				push(@lPatients, $res);
#			}
		}
		$hashHeComp->{$coupleId}->{'nbPat'} = $nbPat . '|' . join('; ', sort(@lPatients));
		push(@lHashHeComp, $hashHeComp->{$coupleId});
	}
	return \@lHashHeComp;
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

sub getHashFormat {
	my $xls = shift;
	my $hash;
	my $format_header = $xls->add_format(border => 1, underline => 1);
	$format_header->set_bold();
	$format_header->set_align('center');
	$format_header->set_fg_color('silver');
	my $format_line = $xls->add_format();
	$format_line->set_align('center');
	$format_line->set_color('black');
	my $format_line_green = $xls->add_format();
	$format_line_green->set_align('center');
	$format_line_green->set_color('green');
	my $format_line_orange = $xls->add_format();
	$format_line_orange->set_align('center');
	$format_line_orange->set_color('orange');
	my $format_line_red = $xls->add_format();
	$format_line_red->set_align('center');
	$format_line_red->set_color('red');
	my $format_merge = $xls->add_format(valign => 'vcenter');
	$format_merge->set_align('center');
	$format_merge->set_color('black');
	$hash->{header} = $format_header;
	$hash->{line} = $format_line;
	$hash->{line_green} = $format_line_green;
	$hash->{line_orange} = $format_line_orange;
	$hash->{line_red} = $format_line_red;
	$hash->{merge}= $format_merge;
	return $hash;
}

sub writeHeader_byVar {
	my ($xls, $xls_pat) = @_;
	my $header = "Access VarId/RsName_1 Consequence_1 VarId/RsName_2 Consequence_2 Nb_Project(s) Nb_Patient(s)";
	my @lNames = split(' ', $header);
	my $hFormat = getHashFormat($xls);
	my $j = 0;
	foreach my $col_name (@lNames) {
		$xls_pat->write(2, $j, $col_name, $hFormat->{header});
		$j++;
	}
}

sub getHashFilters {
	my $hFilters;
	$hFilters->{consequence} = undef;
	if ($filters_consequence) {
		$filters_consequence   =~ s/_/ /g;
		foreach my $name (split(',', $filters_consequence)) {
			$hFilters->{consequence}->{$name} = undef;
		}
	}
	elsif ($filters_consequence_1 and $filters_consequence_2) {
		$filters_consequence_1 =~ s/_/ /g;
		$filters_consequence_2 =~ s/_/ /g;
		foreach my $name (split(',', $filters_consequence_1)) {
			$hFilters->{consequence_1}->{$name} = undef;
		}
		foreach my $name (split(',', $filters_consequence_2)) {
			$hFilters->{consequence_2}->{$name} = undef;
		
		}
	}
	return $hFilters;	
}

sub exportXls {
	my ($originInput, $listHash) = @_;
	my $date = strftime '%d/%m/%Y', localtime;
	my $hFilters = getHashFilters();
	$originInput =~ s/geneId://;
	$originInput =~ s/transId://;
	$originInput =~ s/varId://;
	$originInput =~ s/region://;
	my $workbook;
	if ($outfile) {
		$workbook = Spreadsheet::WriteExcel->new( $outfile );
	}
	else {
		print "Content-type: application/msexcel\n";
		print "Content-Disposition: attachment;filename=DejaVu_HeComp_$originInput.xls\n\n";
		$workbook = Spreadsheet::WriteExcel->new( \*STDOUT );
	}
	my $hFormat = getHashFormat($workbook);
	my $xls_var = $workbook->add_worksheet($originInput);
	$xls_var->write(0, 0, 'DATE', $hFormat->{header});
	$xls_var->write(0, 1, $date, $hFormat->{line});
	writeHeader_byVar($workbook, $xls_var);
	my $i = 3;
	my $format = $hFormat->{line};
	my $format_green = $hFormat->{line_green};
	my $format_orange = $hFormat->{line_orange};
	my $format_red = $hFormat->{line_red};
	$xls_var->write(0, 3, 'FILTERS', $hFormat->{header});
	my $z = 4;
	foreach my $filter_name (keys %{$hFilters->{consequence}}) {
		$xls_var->write(0, $z, $filter_name, $hFormat->{line});
		$z++;
	}
	foreach my $h (@$listHash) {
		my $pass_id = 0;
		if ($filters_var_id) {
			next unless ($h->{id} =~ /$filters_var_id/);
			if ($filters_consequence_2) {
				$pass_id = 1 if ($h->{var_1} eq $filters_var_id);
				$pass_id = 1 if ($h->{var_2} eq $filters_var_id);
			}
		}
		next if ($filters_var_id and $pass_id == 0);
		if ($filters_consequence) {
			my $ok;
			foreach my $cons (split(', ', $h->{consequence_1})) {
				$ok = 1 if (exists $hFilters->{consequence}->{$cons});
			}
			foreach my $cons (split(', ', $h->{consequence_2})) {
				$ok = 1 if (exists $hFilters->{consequence}->{$cons});
			}
			next unless ($ok);
		}
		elsif ($filters_consequence_1 and  $filters_consequence_2) {
			my $ok;
			my @lCons1 = (split(', ', $h->{consequence_1}));
			my @lCons2 = (split(', ', $h->{consequence_2}));
			foreach my $cons1 (@lCons1) {
				if (exists $hFilters->{consequence_1}->{$cons1}) {
					foreach my $cons2 (@lCons2) {
						$ok = 1 if (exists $hFilters->{consequence_2}->{$cons2});
					}
				}
			}
			foreach my $cons1 (@lCons1) {
				if (exists $hFilters->{consequence_2}->{$cons1}) {
					foreach my $cons2 (@lCons2) {
						$ok = 1 if (exists $hFilters->{consequence_1}->{$cons2});
					}
				}
			}
			next unless ($ok);
		}
		my $format_acces = $format_red;
		my $access = 'NO';
		if ($h->{implicated} eq 'X') {
			$access = 'YES';
			$format_acces = $format_green;
		}
		$xls_var->write($i, 0, $access, $format_acces);
		$xls_var->write($i, 1, $h->{var_1}, $format);
		$xls_var->write($i, 2, $h->{consequence_1}, $format);
		$xls_var->write($i, 3, $h->{var_2}, $format);
		$xls_var->write($i, 4, $h->{consequence_2}, $format);
		my @lProj = split('; ', $h->{nbProj});
		$xls_var->write($i, 5, scalar(@lProj), $format);
		my @lPat = split('|', $h->{nbPat});
		$xls_var->write($i, 6, $lPat[0], $format);
		$i++;
	}
	exit(0);
}

