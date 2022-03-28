#!/usr/bin/perl
package check_utils;

use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib $Bin;
use JSON;
use Data::Dumper;
use Term::ANSIColor;
use Text::Table;
use colored;
my $errorValue = 0;
my $warnValue = 0;

sub createJson {
	my ($listHash, $jsonFile) = @_;
	my $printjson;
	my $hash;
	$hash->{'identifier'} = 'patient';
	$hash->{'label'} = 'patient';
	unless (-e $jsonFile) { $hash->{'items'} = $listHash; }
	else {
		my $listHashJson = getListHashFromExistedJsonFile($jsonFile);
		my $listHashMerged = mergeTwoListHashJson($listHash, $listHashJson);
		$hash->{'items'} = $listHashMerged;
	}
	my $cgi = new CGI;
	$printjson = $cgi->header('text/json-comment-filtered');
	$printjson .= encode_json $hash;
	return $printjson;
}


sub getListHashFromExistedJsonFile {
	my ($jsonFile) = shift;
	open (JSON_FILE, "<$jsonFile");
	while(<JSON_FILE>) {
		chomp($_);
		my $line = $_;
		if ($line =~ /{/) {
			my $globalHash = decode_json $line;
			my @listHashJson = @{$globalHash->{'items'}};
			close (JSON_FILE);
			return \@listHashJson;
		}
	}
	return;
}

sub mergeTwoListHashJson {
	my ($listHash, $listHashJson) = @_;
	my $globalHash;
	foreach my $hash (@$listHash) { $globalHash->{$hash->{'patient'}} = $hash; }
	foreach my $hash (@$listHashJson) {
		my $patName = $hash->{'patient'};
		foreach my $key (keys(%$hash)) {
			unless (exists($globalHash->{$patName}->{$key})) { $globalHash->{$patName}->{$key} = $hash->{$key}; }
		}
	}
	my @lHash;
	foreach my $key (keys(%$globalHash)) { push(@lHash, $globalHash->{$key}); }
	return \@lHash;
}

sub getHashFromJsonFile {
	my ($jsonFile) = shift;
	my @listHashJson;
	open (JSON_FILE, "<$jsonFile");
	while(<JSON_FILE>) {
		chomp($_);
		my $line = $_;
		if ($line =~ /{/) {
			my $globalHash = decode_json $line;
			@listHashJson = @{$globalHash->{'items'}};
		}
	}
	close (JSON_FILE);
	my $globalHash;
	foreach my $hash (@listHashJson) { $globalHash->{$hash->{'patient'}} = $hash; }

	return $globalHash;
}

sub getResumeFromJsonFile {
	my ($jsonFile,$pname) = shift;
	print "\n\n### RESUME PROJECT RESULTS BY PATIENT:\n\n";
	print "patient [ped_sex] | bip_sex_check | plink_sex_check | plink_mendelian | freq_dbSNP | cov_15x | cov_mean\n\n";
	my $hash = getHashFromJsonFile($jsonFile);
	foreach my $patName (sort(keys(%$hash))) {
		#if ($hash->{$patName}->{'relation'}) {
			printField_patient($hash, $patName);
			print "\t";
			printField_bip_sex_check($hash, $patName);
			print "\t";
			printField_plink_sex_check($hash, $patName);
			print "\t";
			printField_plink_mendelian($hash, $patName);
			print "\t";
			printField_freq_dbSNP($hash, $patName);
			print "\t";
			printField_cov_15x($hash, $patName);
			print "\t";
			printField_cov_mean($hash, $patName);
			print "\n";
		#}
	}

	print "\nRemark: PLINK_sex says --\n  - M if value >= 0.8\n  - F if value < 0.2\n  - No answer (= ??? in tab) if other value\n\n";
	return ($warnValue, $errorValue);
}
sub printErrorPedigree{
	my ($project) = @_;
	 my $errors = $project->testPedigree();
	my $ped_file = $project->getPedigreeFile();
	if (keys %$errors){
			colored::print_color("red","------------------------------------------------------------------ ");
			colored::print_color("red","problem in pedigree  : $ped_file correct it ");
			colored::print_color("red","------------------------------------------------------------------ ");
			print "\n";
			foreach my $e (keys %$errors){	
				colored::print_color("red","--".$e);
				foreach my $t (@{$errors->{$e}}){
					print "\t";
					colored::print_color("red","*".$t);
				}
				
			}
			return undef;

	} 
	else {
			colored::stabilo("green","Your pedigree seems to be OK ... let's go ");
			
	}
	return 1;
	
}
sub printTable {
	my ($hash,$pname) = @_;
	print "\n\n### RESUME PROJECT RESULTS BY PATIENT: $pname\n\n";
	my @header =("fam","patient","bip_sex_check","plink_sex_check","plink_mendelian","freq_dbSNP","nb SNP","indel dbsnp","nb Indel","cov_15x","cov_mean","cache");
	my @line_null =("------","------","------","------","------","------","------","------","------","------","------");
	
	my $tb = Text::Table->new(
      @header
    );
    
    my @lines;
	my $current_fam;
	foreach my $patName (sort {$hash->{$a}->{family} cmp $hash->{$b}->{family}  || $hash->{$a}->{relation_code} <=>$hash->{$b}->{relation_code} }(keys(%$hash))) {
		#print $patName;
		#if ($hash->{$patName}->{'relation'}) {
		#	print " OK !\n";
		
			next unless exists $hash->{$patName}->{indb};
				if ($current_fam ne $hash->{$patName}->{family}){
				$current_fam =  $hash->{$patName}->{family};
				
					push(@lines,[$hash->{$patName}->{family}."[".$hash->{$patName}->{mendelian_errors}->{nb_snp}."]",@line_null]);
			}
				my @res;
			push(@res," ") ;
			push(@res,printField_patient($hash->{$patName},$patName));
		
			push(@res,printField_sex_check($hash->{$patName}->{check_sex}->{bipd}, $hash->{$patName}->{sex}));
			push(@res,printField_sex_check($hash->{$patName}->{check_sex}->{plink}, $hash->{$patName}->{sex}));
			push(@res,printField_plink_mendelian($hash->{$patName},$hash->{$patName}->{relation_code}));
			push(@res,printField_freq_dbSNP($hash->{$patName}->{snp},"snp"));
			push(@res,printField_nbSNP($hash->{$patName}->{snp}));
			push(@res,printField_freq_dbSNP($hash->{$patName}->{indel},"indel"));
			push(@res,printField_nbSNP($hash->{$patName}->{indel}));
			push(@res,printField_cov_15x($hash->{$patName}->{coverage}));
			push(@res,printField_cov_mean($hash->{$patName}->{coverage}));
			push(@res,printField_cache($hash->{$patName}->{cache}));
			push(@lines,\@res);
		
		}
	$tb->load(@lines);
	print $tb;
	print "\n";
		print "\t\t\t\t-*-*".colored::print_color("green"," Project " .$pname,1)."-*-*-\n";
	return ($warnValue, $errorValue);
}


sub printField_patient {
	my ($hash,$name) = @_;
		my $tsex ='M';
		$tsex ='F' if $hash->{'sex'} eq "2";
	return "$name [$tsex]";
}
sub printField_cache {
	my ($hash) = @_;
	
	my $text = "OK";
	my $color = "green";
	unless ($hash){
		$hash = {};
		$hash->{date} = 999;
	}
	if ( $hash->{date} == 999){
		$text = "CACHE NOT RUN !!!!]";
		 $color = "red";
	}
	elsif (abs($hash->{variations}) > 5 || abs($hash->{indels}) > 50  || $hash->{date} ne 1){
		$text = "ERROR [".$hash->{snps}.",".$hash->{indels}.",".$hash->{date}."]";
		 $color = "red";
	}
	return colored::stabilo($color,$text,1);
}

sub printField_sex_check {
	my ($hash, $sex) = @_;
	my $text ='N.A.';
	my $color = "yellow";
	if ($hash->{sex} == $sex){
		$text = 'OK ['.$hash->{score}.']';
		$color = "green";
		}
		elsif ($hash->{sex} == 0){
			$text = '?? ['.$hash->{score}.']';
		}
		else {
			$color = "red";
			$text = 'NOK ['.$hash->{score}.']';
		}
		return colored::stabilo($color,$text,1);
}


sub printField_plink_mendelian {
	my ($hash,$relation) = @_;
	
	my $text ='N.A.' ;
	my $color = "yellow";
	return colored::stabilo($color,$text,1) unless defined $hash->{mendelian_errors}->{percent}; 
	my $text = '';
	$text .= '[FATHER]' if $relation == 0;
	$text .= '[MOTHER] ' if $relation == 1;
	
	my $pc = $hash->{mendelian_errors}->{percent};
	$text = "[$pc] ".$text;
	 my $text1 = "ERROR";
	 $color = "red";
	
	if ($pc< 1){
	 $text1 = "OK";
		 $color = "green";
	}
	elsif ($pc<2) {
		 $text1 = "WARN";
		 $color = "yellow";
	
	}
	
		return colored::stabilo($color,$text1." ".$text,1);
}

sub printField_freq_dbSNP {
	my ($hash,$type) = @_;
	my $text ='N.A.';
	my $color = "yellow";
	my $limit1 = 96;
	my $limit2 = 92;
	if ($type eq "indel"){
		$limit1 = 70;
	 $limit2 = 50;
	}
	return colored::stabilo($color,$text,1) unless exists $hash->{percent}; 
	
	my $pc = $hash->{percent};
	 $text = "[$pc%]";
	 my $text1 = "";
	 $color = "red";
	if ($pc>$limit1){
		 $text1 = "";
		 $color = "green";
		 
	}
	elsif ($pc>$limit2) {
		 $text1 = "";
		 $color = "yellow";
	}
	return colored::stabilo($color,$text1." ".$text,1);
}
sub printField_nbSNP {
	my ($hash) = @_;
	my $text ='N.A.' ;
	my $color = "yellow";
	return colored::stabilo($color,$text,1)  unless $hash->{nb};
	
	my $zscore = abs($hash->{zscore});
	my $pc = $hash->{nb};
	 $text = "[$pc]";
	 my $text1 = "";
	 $color = "red";
	 
		if ($zscore< 2 ){
		 $color = "green";
	}
	elsif ($zscore<2.5) {
		 $color = "yellow";
	}
	return colored::stabilo($color,$text1." ".$text,1);

	
}

sub printField_cov_15x {
	my ($hash) = @_;
	my $text ='N.A.' ;
	my $color = "yellow";
	return colored::stabilo($color,$text,1)  unless $hash->{15};
	
	my $pc = $hash->{15};
	 $text = "[$pc%]";
	 my $text1 = "";
	 $color = "red";
	 
	if ($pc>96){
		 $text1 = "";
		 $color = "green";
		 
	}
	elsif ($pc>92) {
		 $text1 = "";
		 $color = "yellow";
	}
	return colored::stabilo($color,$text1." ".$text,1);
	
}

sub printField_cov_mean {
	my ($hash) = @_;
	my $text ='N.A.' ;
	my $color = "red";
	return colored::stabilo($color,$text,1) unless $hash->{mean};
	
	my $pc = $hash->{mean};
	
		 $text = "[$pc]";
	 my $text1 = "";
	 $color = "red";
	 
	if ($pc>96){
		 $text1 = "OK";
		 $color = "green";
		 
	}
	elsif ($pc>80) {
		 $text1 = "WARN";
		 $color = "yellow";
	}
	return colored::stabilo($color,$text1." ".$text,1);
	
}



1;