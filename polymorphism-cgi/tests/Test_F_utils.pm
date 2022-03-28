package Test_F_utils;

use strict;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use JSON;
use File::Compare;
use Data::Dumper;
use Spreadsheet::WriteExcel;
use Term::ANSIColor;
use Getopt::Long;
use Test::More qw(no_plan);
my $projectName;
my $hNewKeys;
my $hDelKeys;
my $hDiffRef;
my $hToTest;
my $hResumeErrors;
my $debug;

sub launchTestF_compareTwoJsonText {
	my ($jsonObs, $jsonExp, $testName, $this_debug) = @_;
	$debug = 1 if ($this_debug);
	compareTwoJsonText($jsonObs, $jsonExp);
	testAll($testName);
	resumeErrors($testName);
	printErrors($testName);
}

sub printErrors {
	my ($testName) = @_;
	if (scalar(keys(%$hResumeErrors)) > 0) {
		print "\n\n\n";
		print colored ['black ON_BRIGHT_YELLOW'], "[$testName] Resume ERRORS:\n\n";
		foreach my $testName (sort(keys(%$hResumeErrors))) {
			foreach my $error (sort(@{$hResumeErrors->{$testName}})) {
				warn colored ['black ON_BRIGHT_YELLOW'], $error;
			}
		}
		print color 'reset';
		print "\n\n";
		if ($debug) { warn "\n\n##### DEBUG MODE - DEBUG MODE - DEBUG MODE - DEBUG MODE - DEBUG MODE #####\n\n"; }
		die;
	}
	else {
		print "\n\n\n";
		print colored ['black ON_BRIGHT_GREEN'], "[$testName] All TESTS OK !!";
		print color 'reset';
		print "\n\n";
	}
}

sub resumeErrors {
	my ($testName) = @_;
	warn "\n".'Nb ok: '.scalar(keys %$hToTest);
	warn 'Nb diff prod/dev: '.scalar(keys %$hDiffRef);
	warn 'Nb new key: '.scalar(keys %$hNewKeys);
	warn 'Nb del key: '.scalar(keys %$hDelKeys);
	foreach my $key (sort(keys(%$hToTest))) {
		unless ($hToTest->{$key}->{'exp'} eq $hToTest->{$key}->{'obs'}) {
			push(@{$hResumeErrors->{$testName}}, ' * ['.$testName.'] '.$key.' - PROD: '.$hToTest->{$key}->{'exp'}.'   DEV: '.$hToTest->{$key}->{'obs'});
		}
	}
	foreach my $key (sort(keys(%$hDiffRef))) {
		push(@{$hResumeErrors->{$testName}}, ' * ['.$testName.'] '.$key.' - PROD: '.$hDiffRef->{$key}->{'exp'}.'   DEV: '.$hDiffRef->{$key}->{'obs'});
	}
	foreach my $key (sort(keys(%$hNewKeys))) {
		push(@{$hResumeErrors->{$testName}}, ' * ['.$testName.'] '.$key.' - NEW KEY in dev...');
	}
	foreach my $key (sort(keys(%$hDelKeys))) {
		push(@{$hResumeErrors->{$testName}}, ' * ['.$testName.'] '.$key.' - DELETED KEY in dev...');
	}
}

sub testAll {
	my ($testName) = @_;
	foreach my $key (sort(keys(%$hToTest))) {
		ok($hToTest->{$key}->{'exp'} eq $hToTest->{$key}->{'obs'}, '['.$testName.'] '.$key);
	}
	foreach my $key (sort(keys(%$hDiffRef))) {
		ok($hDiffRef->{$key}->{'exp'} eq $hDiffRef->{$key}->{'obs'}, '['.$testName.'] Pb Type Ref - '.$key);
	}
	foreach my $key (sort(keys(%$hNewKeys))) {
		ok($hNewKeys->{$key}->{'exp'} eq $hNewKeys->{$key}->{'obs'}, '['.$testName.'] New Key - '.$key);
	}
	foreach my $key (sort(keys(%$hDelKeys))) {
		ok($hDelKeys->{$key}->{'exp'} eq $hDelKeys->{$key}->{'obs'}, '['.$testName.'] Deleted Key - '.$key);
	}
}

sub compareTwoJsonText {
	my ($jsonObs, $jsonExp) = @_;
	my $hashObs = decodeFirstJsonText($jsonObs);
	my $hashExp = decodeFirstJsonText($jsonExp);
	if ($debug) {
		my @lTmp = [1..20];
		$hashObs->{'label'} = \@lTmp;
		$hashObs->{'NEW TATA'} = 'je suis nouveau !';
		$hashExp->{'OlD TOTO'} = 'je suis oubliŽ !';
	}
	my $lKeysOk = compareHashes($hashObs, $hashExp, 'HASH->');
}

sub compareHashes {
	my ($hashObs, $hashExp, $oldKeys) = @_;
	my $hKeysObs;
	foreach my $name (keys(%$hashObs)) { $hKeysObs->{$name} = undef; }
	my $hKeysExp;
	foreach my $name (keys(%$hashExp)) { $hKeysExp->{$name} = undef; }
	my @lKeysOk;
	foreach my $name (keys(%$hKeysObs)) {
		if (exists($hKeysExp->{$name})) {
			my $refObs = getVarRefType($hashObs->{$name});
			my $refExp = getVarRefType($hashExp->{$name});
			if ($refObs eq $refExp) {
				delete($hKeysObs->{$name});
				delete($hKeysExp->{$name});
				push(@lKeysOk, $name);	
			}
			else {
				$hDiffRef->{$oldKeys.'{'.$name.'}'}->{'exp'} = $refExp;
				$hDiffRef->{$oldKeys.'{'.$name.'}'}->{'obs'} = $refObs;
				delete($hKeysObs->{$name});
				delete($hKeysExp->{$name});
			}
		}
	}
	checkNewOrDelKeys($oldKeys, $hKeysObs, $hKeysExp);
	foreach my $key (sort @lKeysOk) {
		if (getVarRefType($hashObs->{$key}) eq 'HASH') {
			compareHashes($hashObs->{$key}, $hashExp->{$key}, $oldKeys.'{'.$key.'}->');
		}
		elsif (getVarRefType($hashObs->{$key}) eq 'ARRAY') {
			compareLists($oldKeys.'{'.$key.'}->', $hashExp->{$key}, $hashObs->{$key});
		}
		elsif (getVarRefType($hashObs->{$key}) eq 'STRING') {
			$hToTest->{$oldKeys.'{'.$key.'}'}->{'exp'} = $hashExp->{$key};
			$hToTest->{$oldKeys.'{'.$key.'}'}->{'obs'} = $hashObs->{$key};
		}
	}
}

sub compareLists {
	my ($oldKeys, $lExp, $lObs) = @_;
	my $hashExp = listToHach($lExp);
	my $hashObs = listToHach($lObs);
	compareHashes($hashObs, $hashExp, $oldKeys.'LIST->');
}

sub listToHach {
	my ($list) = @_;
	my $hash;
	foreach my $elmt (@$list) {
		if (getVarRefType($elmt) eq 'HASH') {
			if (exists($elmt->{'id'})) { $hash->{$elmt->{'id'}} = $elmt; }
		}
		elsif (getVarRefType($elmt) eq 'STRING') {
			$hash->{$elmt} = $elmt;
		}
		elsif (getVarRefType($elmt) eq 'ARRAY') {
			warn 'Probleme: list de liste... QUE FAIRE POUR COMPARER ???';
			die;
		}
		else {
			warn 'Probleme: cas non prevu ('.getVarRefType($elmt).')... QUE FAIRE POUR COMPARER ???';
			die;
		}
	}
	return $hash;
}

sub checkNewOrDelKeys {
	my ($oldKeys, $hKeysObs, $hKeysExp) = @_;
	foreach my $name (sort(keys(%$hKeysObs))) {
		$hNewKeys->{$oldKeys.'{'.$name.'}'}->{'obs'} = 'NEW';
		$hNewKeys->{$oldKeys.'{'.$name.'}'}->{'exp'} = undef;
	}
	foreach my $name (sort(keys(%$hKeysExp))) {
		$hDelKeys->{$oldKeys.'{'.$name.'}'}->{'obs'} = undef;
		$hDelKeys->{$oldKeys.'{'.$name.'}'}->{'exp'} = 'DELETED';
	}
}

sub getVarRefType {
	my ($variable) = @_;
	my $ref = ref($variable);
	unless ($ref) { $ref = 'STRING'; }
	return $ref;
}

sub decodeFirstJsonText {
	my ($jsonText) = @_;
	my $hash;
	if ($jsonText =~ /Content-Type/) {
		my @lJson = split("\n", $jsonText);
		$hash = decodeJsonText($lJson[-1]);
	}
	else { $hash = decodeJsonText($jsonText); }
	return $hash;
}

sub decodeJsonText {
	my ($jsonText) = @_;
	my $hash = decode_json($jsonText);
	return $hash;
}

1;