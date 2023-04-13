#!/usr/bin/perl
$|=1;

use strict;
use CGI qw/:standard :html3/;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/";
use lib "$Bin/../../GenBo/lib/obj-nodb/";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
require "$Bin/../../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";
use GBuffer; 
use Data::Dumper;
use Getopt::Long;
use Parallel::ForkManager;
use Tabix;
use JSON;
use export_data;
use CGI::Session;
use Storable qw(store retrieve freeze dclone thaw);
use Compress::Snappy;
use Spreadsheet::WriteExcel;
use POSIX qw(strftime);

my @listSnp;
my @listHashRes;
my $hashRes;
my $buffer = GBuffer->new();
my $query = $buffer->getQuery();
my $cgi = new CGI;
my $input = $cgi->param('input');
my $details = $cgi->param('details');
my $user = $cgi->param('login');
my $pass = $cgi->param('pwd');
my $onlyStrictHo = $cgi->param('onlyStrictHo');
my $onlyHo = $cgi->param('onlyHo');
my $show_large_deletion = $cgi->param('show_large_del');
my $build_use = $cgi->param('build');
my $export_xls = $cgi->param('xls');
my $export_xls_details = $cgi->param('xls_details');
my $filters_sift = $cgi->param('f_sift');
my $filters_polyphen = $cgi->param('f_polyphen');
my $filters_consequence = $cgi->param('f_cons');
my $filters_ids =  $cgi->param('f_ids');
my $filters_polyquery_annot = $cgi->param('filters_annot');
my $filters_freq = $cgi->param('filters_freq');
my $max_dejavu = $cgi->param('max_dejavu');
my $origin_project = $cgi->param('origin_project');
my $annot_version = $cgi->param('force_annot_version');
my $xls_outfile = $cgi->param('xls_outfile');
my $xls_load_session = $cgi->param('xls_load_session');

unless ($input) {
	my $hash = getHashInError();
	push(@listHashRes, $hash);
	printJson(\@listHashRes);
}
unless ($build_use) {
	my $hash = getHashInError();
	push(@listHashRes, $hash);
	printJson(\@listHashRes);
}


my $h_projects_phenotypes;
my $hProjAuthorized;
foreach my $hash (@{$query->getProjectListForUser($user, $pass)}) { $hProjAuthorized->{$hash->{'name'}} = undef; }
my @lAuthProj = keys(%$hProjAuthorized);

my $projectTmp;
if ($origin_project) {
	$projectTmp = $buffer->newProject(-name => $origin_project);
}
else {
	foreach my $project_name (@lAuthProj) {
		my $buffer_test = GBuffer->new();
		my $project_test = $buffer_test->newProject(-name => $project_name);
		if ($project_test->annotation_genome_version() eq $build_use and ($project_test->isDiagnostic() or $project_test->isGenome() or $project_test->isExome())) {
			
			$project_test = undef;
			$buffer_test = undef;
			$projectTmp = $buffer->newProject(-name => $project_name);
			last;
		}
	}
	die unless ($projectTmp);
	if ($annot_version) { $projectTmp->changeAnnotationVersion( $annot_version, 1 ); }
	else {
		my $max_gencode = $projectTmp->buffer->getQuery->getMaxGencodeVersion();
		my $max_annot = $projectTmp->buffer->getQuery->getMaxPublicDatabaseVersion();
		$projectTmp->changeAnnotationVersion( $max_gencode.'.'.$max_annot, 1 );
	}
}

$projectTmp->cgi_object(1);
$projectTmp->version($build_use);
my $tabix = $projectTmp->getSoftware('tabix');

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
if ($xls_load_session) {
	loadSessionsXLS($xls_load_session);
	exit(0);
}

if ( $input =~ 'varId:' ) {
	$show_large_deletion = 1;
	my $varId = $input;
	$varId =~ s/varId://;
	my @lIds = split(',', $varId);
	foreach my $varId (@lIds) {
#		if ($varId =~ /1kg_/) { $varId =~ s/1kg_//; }
#		elsif ($varId =~ /evs_/) { $varId =~ s/evs_//; }
#		elsif (($varId =~ /rs/) or (($varId =~ /TMP_ESP/))) {
#			my $rsName = $varId;
#			my $res = $projectTmp->convertRsNameToIntervalPositions($rsName);
#			if ($res) {
#				$res =~ s/:/_/;
#				$res =~ s/-/_/;
#				$input = 'region:'.uc($res);
#				next;
#			}
#			my $res2 = $projectTmp->buffer->find_variant_id_from_rsname($varId);
#			unless ($res2) {
#				my $hash = getHashInError();
#				$hash->{'varId'} = "ERROR: $varId variant doesn't exist...";
#				push(@listHashRes, $hash);
#				printJson(\@listHashRes);
#			}
#			my @lCol = split('_', $res2);
#			if (length($lCol[2]) == length($lCol[3])) { push(@listSnp, $res2); }
#			else { $input = 'region:'.$lCol[0].'_'.($lCol[1] - 10).'_'.($lCol[1] + 10); }
#		}
#		elsif ($varId =~ /,/) {
#			foreach my $this_id (split(',', $varId)) {
#				push(@listSnp, $this_id);
#			}
#		}
		
		if ($varId =~ /1kg_/) { $varId =~ s/1kg_//; }
		if ($varId =~ /evs_/) { $varId =~ s/evs_//; }
		if (($varId =~ /rs/) or (($varId =~ /TMP_ESP/))) {
			my $rsName = $varId;
			
			my $res = $projectTmp->convertRsNameToIntervalPositions($rsName);
			unless ($res) {
				my $hash = getHashInError();
				$hash->{'varId'} = "ERROR: $rsName variant doesn't exist...";
				push(@listHashRes, $hash);
				printJson(\@listHashRes);
			}
			$res =~ s/:/_/;
			$res =~ s/-/_/;
			$input = 'region:'.uc($res);
		}
		if ($varId =~ /,/) {
			foreach my $this_id (split(',', $varId)) {
				push(@listSnp, $this_id);
			}
		}
		else {
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
}

if ( $input =~ 'transId:' ) {
	my $transId = $input;
	$transId =~ s/transId://;
	if ($transId =~ /ENST/) {
		my $trans;
		eval { $trans = $projectTmp->newTranscript($transId); };
		if($@) {
			my $hash = getHashInError();
			$hash->{'varId'} = "ERROR: $transId transcript doesn't exist...";
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
#	warn $geneId;
	my @lGenesId;
	@lGenesId = split(',', $geneId);
	my @lGenesId2;
	foreach my $geneId (@lGenesId) {
		my $gene;
		eval {
			$gene = $projectTmp->newGene($geneId);
			unless (ref($gene)) {
				my $hash = getHashInError();
				$hash->{'varId'} = "ERROR: $geneId gene doesn't exist... (Annotation Version used: ".$projectTmp->annotation_version().")";
				push(@listHashRes, $hash);
				printJson(\@listHashRes);
			}
		};
		if($@) {
			if (scalar(@lGenesId) == 1) {
				my $hash = getHashInError();
				$hash->{'varId'} = "ERROR: $geneId gene doesn't exist... (Annotation Version used: ".$projectTmp->annotation_version().")";
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

unless ($export_xls) {
	print $cgi->header('text/json-comment-filtered');
	print "{\"progress\":\".";
}
if ($export_xls and not $xls_load_session and not $xls_outfile) {
	print $cgi->header('text/json-comment-filtered');
	print "{\"progress\":\".";
}

my $hPolyqueryFilters;
foreach my $annot_name (split(' ', $filters_polyquery_annot)) {
	if (exists $buffer->config->{ensembl_annotations}->{$annot_name}) {
		my @lTmp = split(';', $buffer->config->{ensembl_annotations}->{$annot_name});
		$hPolyqueryFilters->{$lTmp[1]} = undef;
	}
	else { $hPolyqueryFilters->{$annot_name} = undef; }
}

my $hFreqFilters;
foreach my $freq (split(' ', $filters_freq)) {
	$hFreqFilters->{$freq} = undef;
}

my @lTmp = split(',', $input);
$input = undef;
foreach my $input (sort(@lTmp)) {
	if ($input =~ 'region:') {
		my $posId = $input;
		$posId =~ s/region://;
		$posId =~ s/:/_/;
		$posId =~ s/-/_/;
		@listSnp = @{$projectTmp->getDejaVuIdsFromInterval($posId)};
	}
	my $nbSnpsSearching = scalar(@listSnp);
	my $hRes;
	if (scalar(@listSnp) > 0) {
		foreach my $varId (@listSnp) {
			$projectTmp->print_dot(50);
			my $hashVarId = $projectTmp->getDejaVuInfos($varId);
			if ($details eq 'true') {
				foreach my $hDetails (_getVarTranscriptsDetails($varId)) { $hashRes->{$hDetails->{'transcript'}} = $hDetails; }
			}
			elsif ($hashVarId) { 
				my $hashVar = _getVarGeneralDetails($originInput, $varId, $hashVarId);
				if ($hashVar) {
					my $nb_ok;
					if ($hPolyqueryFilters) {
						foreach my $annot_name (split(', ', $hashVar->{'consequence'})) {
							$nb_ok++ unless (exists $hPolyqueryFilters->{$annot_name});
						}
					}
					else { $nb_ok++; }
					if ($max_dejavu and $max_dejavu < int($hashVar->{'nbProj_value'})) { $nb_ok = undef; }
					if ($filters_freq and exists $hFreqFilters->{$hashVar->{'cat_freq'}}) { $nb_ok = undef; }
					$hashRes->{$varId} = $hashVar if ($nb_ok);
				}
			}
		}
	}
}

if (scalar(keys(%$hashRes)) == 0) { 
	my $hash = getHashInError();
	push(@listHashRes, $hash);
}
else { @listHashRes = values(%$hashRes); }
if ($export_xls) {
	if ($xls_outfile) {
		exportXls($originInput, \@listHashRes);
	}
	else {
		saveSessionXLS(\@listHashRes);
	}
}
else { printJson(\@listHashRes); }




sub loadSessionsXLS {
	my ($sid) = @_;
	my $tmp_dir = $projectTmp->buffer()->getDataDirectory("cache").'/tmp/';
    my $session = new CGI::Session(undef, $sid, {Directory=>$tmp_dir});
    my $h = thaw(decompress $session->param('list'));
    $session->delete();
    exportXls($originInput, $h->{list});
	exit(0);
}

sub saveSessionXLS {
	my ($listHashRes) = @_;
	my ($h, $ok);
	print ".";
	$h->{list} = $listHashRes;
	my $tmp_dir = $projectTmp->buffer()->getDataDirectory("cache").'/tmp/';
	unless (-d $tmp_dir) {
		`mkdir $tmp_dir`;
		`chmod 777 $tmp_dir`;
	}
	my $session = new CGI::Session(undef, $cgi, {Directory=>$tmp_dir});
	$session->param('list', compress(freeze $h));
	my $sid = $session->id;
	print '@@@';
	print $sid;
    exit(0);
}

sub getHashInError {
	my $hash;
	$hash->{'varId'} = "No Result !";
	$hash->{'implicated'} = ".";
	$hash->{'type'} = ".";
	$hash->{'chr'} = ".";
	$hash->{'pos'} = ".";
	$hash->{'locus'} = ".";
	$hash->{'gene'} = ".";
	$hash->{'transcript'} = ".";
	$hash->{'protein'} = "0";
	$hash->{'nbProj'} = "0";
	$hash->{'nbPat'} = "0";
	$hash->{'nbHo'} = "0";
	$hash->{'nbHe'} = "0";
	$hash->{'consequence'} = ".";
	$hash->{'polyphen'} = ".";
	$hash->{'sift'} = ".";
	$hash->{'freqDbSnp'} = ".";
	$hash->{'ncboost'} = ".";
	$hash->{'cadd_score'} = ".";
	$hash->{'found_phenotypes'} = '.';
	$hash->{'gnomad_ac'} = ".";
	$hash->{'gnomad_ho'} = ".";
	$hash->{'gnomad_max'} = ".";
	$hash->{'gnomad_min'} = ".";
	$hash->{'gnomad_an'} = ".";
	return $hash;
}

sub _getVarGeneralDetails {
	my ($input, $varId, $hash_input) = @_;
	my ($proj, $pat, $nbHo, $nbHe, $myProj, $nbProj, $pheno) = _getPatProjHoHeNbFromHash($hash_input);
	my @lTmp = split('_', $varId);
	my $hash;
	my $var = $projectTmp->_newVariant($varId);
	if ($var->isVariation()) { $hash->{'type'} = 'snp'; }
	if ($var->isInsertion()) { $hash->{'type'} = 'ins'; }
	if ($var->isDeletion()) {
		if ($var->checkLargeDeletion_newVariant()) {
			$hash->{'type'} = 'large_del';
			return unless ($show_large_deletion);
		}
		else { $hash->{'type'} = 'del'; }
	}
	$hash->{'chr'} = $var->getChromosomes()->[0]->name();
	$hash->{'pos'} = $var->start();
	$hash->{'locus'} = 'chr'.$var->getChromosomes()->[0]->name().':'.$var->start();
	my $cadd_score = $var->cadd_score();
	if ($cadd_score and $cadd_score ne '-1' and $cadd_score ne '-') {
		if (int($cadd_score) < 10) { $hash->{'cadd_score'} = '0'.$cadd_score; }
		else { $hash->{'cadd_score'} = $cadd_score; }
	}
	$hash->{'cadd_score'} = '.' unless ($hash->{'cadd_score'});
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
	#warn $objAnalyse;
	my @lFilters;
	$hash->{'varId'} = $varId;
	if ($var->rs_name()) { $hash->{'varId'} .= ';'.$var->rs_name(); }
	elsif ($var->isInsertion()) { $hash->{'varId'} .= ';'.$var->name(); }
	elsif ($var->isDeletion()) { $hash->{'varId'} .= ';'.$var->name(); }
	$hash->{'id'} .= $varId;
	$hash->{'rsName'} = '';
	if ($var->rs_name()) {
		$hash->{'rsName'} .= $var->rs_name();
		push(@lFilters, 'rsname');
	}
	else { push(@lFilters, 'notRsname'); }
	$hash->{'nbProj'} = $proj;
	$hash->{'nbProj_value'} = $nbProj;
	$hash->{'implicated'} = '';
	if ($myProj) {
		$hash->{'implicated'} = 'X';
		push(@lFilters, 'myVar');
	}
	else { push(@lFilters, 'notMyVar'); }
	$hash->{'nbPat'} = $pat;
	$hash->{'nbHo'} = $nbHo;
	$hash->{'nbHe'} = $nbHe;
	
	if ($var->getProject->hasHgmdAccess($user)) {
		$hash->{'hgmd_class'} = '.';
		if ($var->hgmd()) {
			my $h;
			$h->{obj} = $var;
			my ($class, $cmd) = update::hgmd($var->getProject,$h,'return_cmd');
			$hash->{'hgmd_class'} = $class.';'.$cmd;
			if ($var->isNewHgmd()) {
				$hash->{'hgmd_class'} .= ";New!";
			}
		}
	}
	else {
		$hash->{'hgmd_class'} = "N.A.";
	}
	
	if ($nbHo > 0) { push(@lFilters, 'ho'); }
	if ($nbHe > 0) {
		push(@lFilters, 'he');
		if ($onlyStrictHo) { return; }
	}
	if (($nbHo == 0) and ($nbHe > 0)) {
		push(@lFilters, 'strictHe');
		if ($onlyHo) { return; }
	}
	if (($nbHe == 0) and ($nbHo > 0)) { push(@lFilters, 'strictHo'); }
	if ($input =~ /geneId:/) {
		if ($objAnalyse->external_name()) { $hash->{'gene'} = $objAnalyse->external_name(); }
		else { $hash->{'gene'} = $objAnalyse->name(); }
	}
	else {
		my @lGenes;
		foreach my $gene (@lGenesObj) { push(@lGenes, $gene->external_name()); }
		if (scalar(@lGenes) == 0) { $hash->{'gene'} = 0; }
		else { $hash->{'gene'} = join('; ', sort(@lGenes)); }
		
	}
	if ($input =~ /transId:/) {
		if ($objAnalyse->external_name()) { $hash->{'transcript'} = $objAnalyse->external_name().' (='.$objAnalyse->name().')'; } 
		else { $hash->{'transcript'} = $objAnalyse->name(); }
		if ($export_xls) {
			my $tr_nomenclature = $var->getNomenclature($objAnalyse);
			if (length($tr_nomenclature) > 30) {
				$tr_nomenclature = substr $tr_nomenclature, 0, 30;
				$tr_nomenclature .= '...';
			}
			$hash->{'transcript'} .= ' ['.$tr_nomenclature.']';
		}
	}
#	else {
#		my @lTrans;
#		foreach my $trans (@lTransObj) {
#			my $text = $trans->name();
#			$text .= ' ['.$var->getNomenclature($trans).']';
#			#warn $text;
#			push(@lTrans, $text);
#		}
#		if (scalar(@lTransObj) == 0) { $hash->{'transcript'} = 0; }
#		else { $hash->{'transcript'} = join('; ', sort(@lTrans)); }
#	}
#	my $hProt;
#	foreach my $trans (@lTransObj) {
#		foreach my $prot (@{$trans->getProteins()}) {
#			unless (exists $hProt->{$prot->external_name()}) {
#				$hProt->{$prot->external_name()} = undef;
#				if ($export_xls) {
#					$hProt->{$prot->external_name()} = $var->protein_nomenclature($prot);
#				}
#			}
#		}
#	}
#	if (scalar(keys(%$hProt)) == 0) { $hash->{'protein'} = 0; }
#	else {
#		if ($export_xls) {
#			$hash->{'protein'} = undef;
#			foreach my $prot_name (sort keys(%$hProt)) {
#				if ($hash->{'protein'}) { $hash->{'protein'} .= ';'; }
#				$hash->{'protein'} .= $prot_name.' ['.$hProt->{$prot_name}.']';
#			}
#		}
#		else {
#			$hash->{'protein'} = join(';', sort keys(%$hProt));
#		}
#	}
	$hash->{'polyphen'} = '';
	$hash->{'freqDbSnp'} = '.';
	if ($var->frequency()) {
		my $this_freq = $var->frequency();
		if ($this_freq eq '-1') { $hash->{'freqDbSnp'} = ' '; }
		else { $hash->{'freqDbSnp'} = sprintf("%.5f", ($this_freq * 100)).' %'; }
		$hash->{'cat_freq'} = $var->categorie_frequency();
	}
	eval {
		my @lConsOk;
		foreach my $cons (split(',', lc($var->variationType($objAnalyse)))) {
			my @lCons = split(';', $projectTmp->buffer->config->{ensembl_annotations}->{$cons});
			push(@lConsOk, $lCons[1]);
		}
		$hash->{'consequence'} = join(', ', @lConsOk);
	};
	if ($@) { $hash->{'consequence'} = 'Error annotation...'; }
	my $polyphenScore = "-";
	if (($objAnalyse) and ($objAnalyse->isTranscript())) { $polyphenScore = $var->polyphenScore($objAnalyse); }
	eval {
		my $polyStatus = $var->polyphenStatus($objAnalyse);
		push(@lFilters, 'polyphen'.$polyStatus);
		$hash->{'polyphen'} = $polyStatus . "+" . $polyphenScore;
	};
	if ($@) { $hash->{'polyphen'} = '4+-'; }
	my $siftScore = "-";
	if (($objAnalyse) and ($objAnalyse->isTranscript())) { eval { $siftScore = $var->siftScore($objAnalyse); } }
	eval {
		my $siftStatus = $var->siftStatus($objAnalyse);
		push(@lFilters, 'sift'.$siftStatus);
		$hash->{'sift'} = $siftStatus . "+" . $siftScore;
	};
	if ($@) { $hash->{'sift'} = '4+-'; }
	my $ncboost_score = $var->ncboost_score();
	if ($ncboost_score eq '-') { $hash->{'ncboost'} = '-'; }
	if ($ncboost_score eq ';') { $hash->{'ncboost'} = '.'; }
	else { $hash->{'ncboost'} =  $var->ncboost_category().';'.$ncboost_score;}
	$hash->{'alamut'} = $var->alamut_id();
	$hash->{'filters'} = join(',', sort(@lFilters));
	$hash->{'found_phenotypes'} = '.';
	$hash->{'found_phenotypes'} = $pheno if ($pheno);
	$hash->{'gnomad_ac'} = '.';
	$hash->{'gnomad_ac'} = $var->getGnomadAC() if ($var->getGnomadAC());
	$hash->{'gnomad_ho'} = '.';
	$hash->{'gnomad_ho'} = $var->getGnomadHO() if ($var->getGnomadHO());
	$hash->{'gnomad_max'} = '.';
	$hash->{'gnomad_max'} = $var->max_pop_name().' '.sprintf("%.4f", $var->max_pop_freq) if ($var->max_pop_name());
	$hash->{'gnomad_min'} = '.';
	$hash->{'gnomad_min'} = $var->min_pop_name().' '.sprintf("%.4f", $var->min_pop_freq) if ($var->min_pop_name());
	$hash->{'gnomad_an'} = '.';
	$hash->{'gnomad_an'} = $var->getGnomadAN() if ($var->getGnomadAN());
	$var = undef;
	if ($export_xls) {
		foreach my $proj_name (sort keys %{$hash_input}) {
			my $thisProject;
			$thisProject = $buffer->newProject(-name => $proj_name) if $export_xls_details;
			foreach my $pat_info (split(';', $hash_input->{$proj_name}->{string})) {
				my ($pat_name, $he_ho) = split(':', $pat_info);
				$hash->{my_patients}->{$proj_name}->{$pat_name}->{status} = '';
				if ($export_xls_details) {
					my $patient;
					eval { $patient = $thisProject->getPatient($pat_name) };
					if ($@) { next; }
					$hash->{my_patients}->{$proj_name}->{$pat_name}->{status} = $patient->status();
				}
				if    ($he_ho eq '1') { $he_ho = 'Ho'; }
				elsif ($he_ho eq '2') { $he_ho = 'He'; }
				$hash->{my_patients}->{$proj_name}->{$pat_name}->{heho} = $he_ho;
				#TODO: ici
			}
			unless (exists $hash->{projects}->{$proj_name}) {
				my $hmails;
				foreach my $h (@{$buffer->getQuery()->getOwnerProject_byName($proj_name)}) {
					my @lContact = split('@', $h->{email});
					$hmails->{$lContact[0]} = undef;
				} 
				$hProjAuthorized->{projects}->{$proj_name}->{contacts} = join(';', sort keys %$hmails);
			}
			$thisProject = undef;
		}
	}
	return $hash;
}

sub _getVarTranscriptsDetails {
	my ($varId) = shift;
	my $var = $projectTmp->_newVariant($varId);
	my $hashRes =  _getVarObjectDetails($var, 'transcript', $var->getTranscripts());
	return values(%$hashRes);
}

sub _getVarObjectDetails {
	my ($var, $type, $listObjects) = @_;
	my $hash;
	foreach my $obj (@$listObjects) {
		my $name = $obj->name();
		eval { $hash->{$name}->{'protein'} = $obj->getProtein()->name().'+'.$obj->getProtein()->external_name().')'; };
		if ($@) { $hash->{$name}->{'protein'} = ''; }
		$hash->{$name}->{'transcript'} = $obj->name().'+'.$obj->external_name();
		my $gene = $obj->getGenes()->[0];
		$hash->{$name}->{'gene'} = $gene->name().';'.$gene->external_name();
		$hash->{$name}->{'description'} = $gene->description();
		eval {
			my @lConsOk;
			foreach my $cons (split(',', lc($var->variationType($obj)))) {
				my @lCons = split(';', $projectTmp->buffer->config->{ensembl_annotations}->{$cons});
				push(@lConsOk, $lCons[1]);
			}
			$hash->{$name}->{'consequence'} = join(', ', @lConsOk);
		};
		if ($@) { $hash->{$name}->{'consequence'} = 'Error annotation...'; }
		my $polyphenScore = "-";
		if ($obj) { eval { $polyphenScore = $var->polyphenScore($obj); } }
		eval { $hash->{$name}->{'polyphen'} = $var->polyphenStatus($obj) . "+" . $polyphenScore; };
		if ($@) { $hash->{$name}->{'polyphen'} = '4+-'; }
		my $siftScore = "-";
		if ($obj) { eval { $siftScore = $var->siftScore($obj); } }
		eval { $hash->{$name}->{'sift'} = $var->siftStatus($obj) . "+" . $siftScore; };
		if ($@) { $hash->{$name}->{'sift'} = '4+-'; }
		
		if ($var->isCoding()) {
			eval { $hash->{$name}->{'nomenclature'} = $var->getNomenclature($obj); };
			eval { $hash->{$name}->{'cdna_pos'} = $obj->translate_position($var->start()) };
			eval { $hash->{$name}->{'cds_pos'} = $var->getOrfPosition($obj->getProtein()); };
			eval { $hash->{$name}->{'prot_pos'} = $var->getProteinPosition($obj->getProtein()); };
			eval { $hash->{$name}->{'exon'} = $obj->findExonNumber($var->start()) };
		}
	}
	return $hash;
}

sub _getPatProjHoHeNbFromHash {
	my ($hash) = shift;
	#warn Dumper $hash; die;
	my ($myProj, $hPatiens);
	my (@lProjects, @lPatients);
	my $nbPat = 0;
	my $nbHo  = 0;
	my $nbHe  = 0;
	my $h_phenotypes;
	foreach my $projName (keys %$hash) {
		my ($authorized, $res);
		if (exists ($hProjAuthorized->{$projName})) {
			$authorized = 1;
			$myProj = 1;
			push(@lProjects, $projName);
		}
		else { push(@lProjects, 'X'); }
		my @lPatTmp;
		foreach my $patString (split(';', $hash->{$projName}->{string})) {
			$nbPat++;
			my ($patName, $heho) = split(':', $patString);
			#warn $patString;
			$nbHo++ if ($heho eq '1');
			$nbHe++ if ($heho eq '2');
			$hPatiens->{$projName}->{$patName} = undef;
#			if ($authorized) {
				if ($heho eq '1') { push(@lPatTmp, '(Ho)'.$patName); }
				else { push(@lPatTmp, $patName); }
#			}
#			else { push(@lPatTmp, 'X') }
		}
		$res = $projName.' [';
#		if ($authorized) { $res = $projName.' ['; }
#		else { $res = 'X ['; }
		$res .= join(', ', @lPatTmp).']';
		push(@lPatients, $res);
		
		my $pheno_name;
		if (exists $h_projects_phenotypes->{$projName}) {
			$pheno_name = $h_projects_phenotypes->{$projName};
		}
		else {
			my $buffer_pheno = GBuffer->new();
			my $project_pheno = $buffer_pheno->newProject(-name => $projName);
			my @lPheno;
			foreach my $pheno_obj (@{$project_pheno->getPhenotypes()}) {
				push(@lPheno, $pheno_obj->name());
			}
			$pheno_name = join(', ', sort @lPheno) if (@lPheno);
			$project_pheno = undef;
			$buffer_pheno = undef;
			$h_projects_phenotypes->{$projName} = $pheno_name;
		}
		$h_phenotypes->{$pheno_name}++ if ($pheno_name);
	}
	my $nb_proj = scalar(@lProjects);
	my $pheno;
	foreach my $pheno_name (sort keys %$h_phenotypes) {
		$pheno .= ', ' if ($pheno);
		$pheno .= ucfirst($pheno_name).' ('.$h_phenotypes->{$pheno_name}."/$nb_proj)";
	}
	my $proj = join('; ', sort(@lProjects));
	my $pat = $nbPat . '|' . join('; ', sort(@lPatients));
	return ($proj, $pat, $nbHo, $nbHe, $myProj, $nb_proj, $pheno);
}

sub printJson {
	my ($listHash) = @_;
	my $hash;
	if ($details eq 'true') {
		$hash->{'identifier'} = 'transcript';
		$hash->{'label'} = 'transcript';
	}
	else {
		$hash->{'identifier'} = 'varId';
		$hash->{'label'} = 'varId';
	}
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
	my ($xls, $xls_pat, $i) = @_;
	my $header = "Access VarId RsName Type Chr Position Gene(s) Transcript(s) Protein(s) Nb_Project(s) Nb_Patient(s) Nb_Ho Nb_He Consequence Polyphen_Status Sift_Status %_Freq_dbSnp Ncboost Cadd";
	my @lNames = split(' ', $header);
	my $hFormat = getHashFormat($xls);
	my $j = 0;
	foreach my $col_name (@lNames) {
		$xls_pat->write($i, $j, $col_name, $hFormat->{header});
		$j++;
	}
}

sub getHashFilters {
	my $hFilters;
	$hFilters->{ids} = undef;
	$hFilters->{sift} = undef;
	$hFilters->{polyphen} = undef;
	$hFilters->{consequence} = undef;
	foreach my $name (split(',', $filters_ids)) {
		$hFilters->{ids}->{$name} = undef;
	}
	foreach my $name (split(',', $filters_sift)) {
		$hFilters->{sift}->{$name} = undef;
	}
	foreach my $name (split(',', $filters_polyphen)) {
		$hFilters->{polyphen}->{$name} = undef;
	}
	foreach my $name (split(',', $filters_consequence)) {
		$hFilters->{consequence}->{$name} = undef;
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
	if ($xls_outfile) {
		$workbook = Spreadsheet::WriteExcel->new( $xls_outfile );
	}
	else {
		print "Content-type: application/msexcel\n";
		print "Content-Disposition: attachment;filename=DejaVu_$originInput.xls\n\n";
		$workbook = Spreadsheet::WriteExcel->new( \*STDOUT );
	}
	my $hFormat = getHashFormat($workbook);
	my $xls_var = $workbook->add_worksheet($originInput);
	my $xls_dejavu = $workbook->add_worksheet('PATIENTS');
	$xls_var->write(0, 0, 'DATE', $hFormat->{header});
	$xls_var->write(0, 1, $date, $hFormat->{line});
	my $i = 0;
	my @lFilters_cons = keys %{$hFilters->{consequence}};
	my @lFilters_sift = keys %{$hFilters->{sift}};
	my @lFilters_polyphen = keys %{$hFilters->{polyphen}};
	my $z = 4;
	if ($filters_ids) {
		$xls_var->write($i, 3, 'VAR IDS', $hFormat->{header});
		foreach my $filter_name (split(',', $filters_ids)) {
			$xls_var->write($i, $z, $filter_name, $hFormat->{line});
			$z++;
		}
		$z = 4;
		$i++;
		
	}
	if ($onlyHo) {
		$xls_var->write($i, 3, 'ONLY HO', $hFormat->{header});
		$xls_var->write($i, 4, 'Yes', $hFormat->{line});
		$i++;
	}
	if ($onlyStrictHo) {
		$xls_var->write($i, 3, 'HO EXCLUSIVE', $hFormat->{header});
		$xls_var->write($i, 4, 'Yes', $hFormat->{line});
		$i++;
	}
	if (scalar(@lFilters_cons) > 0) {
		$xls_var->write($i, 3, 'FILTERS', $hFormat->{header});
		foreach my $filter_name (sort keys %{$hFilters->{consequence}}) {
			$xls_var->write($i, $z, $filter_name, $hFormat->{line});
			$z++;
		}
		$z = 4;
		$i++;
	}
	if (scalar(@lFilters_sift) > 0) {
		$xls_var->write($i, 3, 'SIFT', $hFormat->{header});
		foreach my $filter_name (sort keys %{$hFilters->{sift}}) {
			$filter_name = 'None' if ($filter_name eq '0');
			$filter_name = 'Benign' if ($filter_name eq '1');
			$filter_name = 'Deleterious' if ($filter_name eq '2');
			$filter_name = 'Error' if ($filter_name eq '4');
			$xls_var->write($i, $z, $filter_name, $hFormat->{line});
			$z++;
		}
		$z = 4;
		$i++;
	}
	if (scalar(@lFilters_polyphen) > 0) {
		$xls_var->write($i, 3, 'POLYPHEN', $hFormat->{header});
		foreach my $filter_name (sort keys %{$hFilters->{polyphen}}) {
			$filter_name = 'None' if ($filter_name eq '0');
			$filter_name = 'Benign' if ($filter_name eq '1');
			$filter_name = 'Possibly Damaging' if ($filter_name eq '2');
			$filter_name = 'Probably Damaging' if ($filter_name eq '3');
			$filter_name = 'Error' if ($filter_name eq '4');
			$xls_var->write($i, $z, $filter_name, $hFormat->{line});
			$z++;
		}
		$z = 4;
		$i++;
	}
	$i++;
	writeHeader_byVar($workbook, $xls_var, $i);
	$i++;
	my $format = $hFormat->{line};
	my $format_green = $hFormat->{line_green};
	my $format_orange = $hFormat->{line_orange};
	my $format_red = $hFormat->{line_red};
	foreach my $h (@$listHash) {
		
#		warn Dumper $h; die;
		
		my @lTmp_s = split('\+', $h->{sift});
		my @lTmp_p = split('\+', $h->{polyphen});
		if ($filters_ids) {
			my $ok;
			foreach my $name (split(';', $h->{varId})) {
				$ok = 1 if (exists $hFilters->{ids}->{$name});
			}
			next unless ($ok);
		}
		if ($filters_sift) {
			next unless (exists $hFilters->{sift}->{$lTmp_s[0]});
		}
		if ($filters_polyphen) {
			next unless (exists $hFilters->{polyphen}->{$lTmp_p[0]});
		}
		if ($filters_consequence) {
			my $ok;
			foreach my $name (split(', ', $h->{consequence})) {
				$ok++ if (exists $hFilters->{consequence}->{$name});
			}
			next unless ($ok == scalar(keys %{$hFilters->{consequence}}));
		}
		my $format_acces = $format_red;
		my $access = 'NO';
		if ($h->{implicated} eq 'X') {
			$access = 'YES';
			$format_acces = $format_green;
		}
		$xls_var->write($i, 0, $access, $format_acces);
		my ($varId, $rsName) = split(';', $h->{varId});
		$xls_var->write($i, 1, $varId, $format);
		$xls_var->write($i, 2, $rsName, $format);
		$xls_var->write($i, 3, $h->{type}, $format);
		$xls_var->write($i, 4, $h->{chr}, $format);
		$xls_var->write($i, 5, $h->{pos}, $format);
		$xls_var->write($i, 6, $h->{gene}, $format);
		$xls_var->write($i, 7, $h->{transcript}, $format);
		$xls_var->write($i, 8, $h->{protein}, $format);
		my @lProj = split('; ', $h->{nbProj});
		$xls_var->write($i, 9, scalar(@lProj), $format);
		$xls_var->write($i, 10, (int($h->{nbHo}) + int($h->{nbHe})), $format);
		$xls_var->write($i, 11, $h->{nbHo}, $format);
		$xls_var->write($i, 12, $h->{nbHe}, $format);
		$xls_var->write($i, 13, $h->{consequence}, $format);
		my ($polyphen, $format_polyphen);
		if ($lTmp_p[0] eq '1') {
			$polyphen = 'Benign';
			$format_polyphen = $format_green;
		}
		elsif ($lTmp_p[0] eq '2') {
			$polyphen = 'Possibly Damaging';
			$format_polyphen = $format_orange;
		}
		elsif ($lTmp_p[0] eq '3') {
			$polyphen = 'Probably Damaging';
			$format_polyphen = $format_red;
		}
		elsif ($lTmp_p[0] eq '4') {
			$polyphen = 'Error...';
			$format_polyphen = $format;
		}
		else { $format_polyphen = $format; }
		$xls_var->write($i, 14, $polyphen, $format_polyphen);
		my ($sift, $format_sift);
		if ($lTmp_s[0] eq '1') {
			$sift = 'Benign';
			$format_sift = $format_green;
		}
		elsif ($lTmp_s[0] eq '2') {
			$sift = 'Deleterious';
			$format_sift = $format_red;
		}
		elsif ($lTmp_s[0] eq '4') {
			$sift = 'Error...';
			$format_sift = $format;
		}
		else { $format_sift = $format; }
		$xls_var->write($i, 15, $sift, $format_sift);
		$xls_var->write($i, 16, $h->{freqDbSnp}, $format);
		$xls_var->write($i, 17, $h->{ncboost}, $format);
		$xls_var->write($i, 18, $h->{cadd_score}, $format);
		$i++;
	}
	$xls_dejavu->write(0, 0, 'DATE', $hFormat->{header});
	$xls_dejavu->write(0, 1, $date, $hFormat->{line});
	$xls_dejavu->write(1, 0, 'Access', $hFormat->{header});
	$xls_dejavu->write(1, 1, 'VarId', $hFormat->{header});
	$xls_dejavu->write(1, 2, 'RsName', $hFormat->{header});
	$xls_dejavu->write(1, 3, 'Projects', $hFormat->{header});
	$xls_dejavu->write(1, 4, 'Contacts', $hFormat->{header});
	$xls_dejavu->write(1, 5, 'Patients', $hFormat->{header});
	my $i = 2;
	foreach my $h (@$listHash) {
		next unless (exists $h->{my_patients});
		my $format_acces = $format_red;
		my $access = 'NO';
		if ($h->{implicated} eq 'X') {
			$access = 'YES';
			$format_acces = $format_green;
		}
		my ($varId, $rsName) = split(';', $h->{varId});
		foreach my $proj_name (keys %{$h->{my_patients}}) {
			my $contacts = $hProjAuthorized->{projects}->{$proj_name}->{contacts};
			$xls_dejavu->write($i, 0, $access, $format_acces);
			$xls_dejavu->write($i, 1, $varId, $format);
			$xls_dejavu->write($i, 2, $rsName, $format);
			$xls_dejavu->write($i, 3, $proj_name, $format);
			$xls_dejavu->write($i, 4, $contacts, $format);
			my $j = 5;
			foreach my $pat_name (keys %{$h->{my_patients}->{$proj_name}}) {
				my $he_ho = $h->{my_patients}->{$proj_name}->{$pat_name}->{heho};
				my $status = $h->{my_patients}->{$proj_name}->{$pat_name}->{status};
				if ($status == 1) {
					$status = "healthy";
				}
				if ($status == 2) {
					$status = "ill";
				}
				my $text = "$pat_name:$status:$he_ho";
				$xls_dejavu->write($i, $j, $text, $format);
				$j++;
			}
			$i++;
		}
	}
	exit(0);
}

