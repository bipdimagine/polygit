#!/usr/bin/perl

use strict;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use JSON;
use File::Compare;
use Data::Dumper;
use Test::More;
use Spreadsheet::WriteExcel;
use Term::ANSIColor;
use Getopt::Long;
use Test_F_utils;

my $projectName = 'NGS2013_0318';
my $hNewKeys;
my $hDelKeys;
my $hDiffRef;
my $hToTest;
my $hResumeErrors;
my $debug;

GetOptions(
	'project=s' => \$projectName,
	'debug=s' 	=> \$debug,
);

if ($debug) { warn "\n\n##### DEBUG MODE - DEBUG MODE - DEBUG MODE - DEBUG MODE - DEBUG MODE #####\n\n"; }
Test_F_genes_json_indiv_default();
Test_F_genes_json_indiv_default_coding();
Test_F_genes_json_indiv_all();
Test_F_genes_json_fam_default();
Test_F_genes_json_fam_default_coding();
Test_F_genes_json_fam_all();





########## METHODS ##########





sub Test_F_genes_json_indiv_default {
	my $args = ' stat=all';
	$args .= ' filter_type_variation=dbsnp+dbsnp_1p+dbsnp_none+1000genomes+1000genomes_1p+evs+evs_1p+silent+intergenic+intronic+pseudogene+loh+dbl_evt+polyphen_sift';
	$args .= ' mode=ind';
	$args .= ' level_ind=variation';
	$args .= ' level_fam=variation';
	Test_F_genes_json($args, 'indiv_default');
}

sub Test_F_genes_json_indiv_default_coding {
	my $args = ' stat=all';
	$args .= ' filter_type_variation=dbsnp+dbsnp_1p+dbsnp_none+1000genomes+1000genomes_1p+evs+evs_1p+utr+splicing+essential_splicing+intergenic+intronic+ncrna+maturemirna+pseudogene+loh+dbl_evt+polyphen_sift';
	$args .= ' mode=ind';
	$args .= ' level_ind=variation';
	$args .= ' level_fam=variation';
	Test_F_genes_json($args, 'indiv_default_coding');
}

sub Test_F_genes_json_indiv_all {
	my $args = ' stat=all';
	$args .= ' filter_type_variation=filter_type_variation=loh+dbl_evt+polyphen_sift';
	$args .= ' mode=ind';
	$args .= ' level_ind=variation';
	$args .= ' level_fam=variation';
	Test_F_genes_json($args, 'indiv_all');
}

sub Test_F_genes_json_fam_default {
	my $args = ' stat=all';
	$args .= ' filter_type_variation=dbsnp+dbsnp_1p+dbsnp_none+1000genomes+1000genomes_1p+evs+evs_1p+silent+intergenic+intronic+pseudogene+loh+dbl_evt+polyphen_sift';
	$args .= ' mode=fam';
	$args .= ' level_ind=variation';
	$args .= ' level_fam=variation';
	Test_F_genes_json($args, 'fam_default');
}

sub Test_F_genes_json_fam_default_coding {
	my $args = ' stat=all';
	$args .= ' filter_type_variation=dbsnp+dbsnp_1p+dbsnp_none+1000genomes+1000genomes_1p+evs+evs_1p+utr+splicing+essential_splicing+intergenic+intronic+ncrna+maturemirna+pseudogene+loh+dbl_evt+polyphen_sift';
	$args .= ' mode=fam';
	$args .= ' level_ind=variation';
	$args .= ' level_fam=variation';
	Test_F_genes_json($args, 'fam_default_coding');
}

sub Test_F_genes_json_fam_all {
	my $args = ' stat=all';
	$args .= ' filter_type_variation=filter_type_variation=loh+dbl_evt+polyphen_sift';
	$args .= ' mode=fam';
	$args .= ' level_ind=variation';
	$args .= ' level_fam=variation';
	Test_F_genes_json($args, 'fam_all');
}


sub Test_F_genes_json {
	my ($args, $testName) = @_;
	my $cmd1 = "/bip-d/perl/polymorphism-cgi/json_output_nodb/genes_json.pl project=$projectName" . $args;
	print "\n### TESTS $testName\n";
	print " -> Check PROD environment\n";
	my $jsonObs = `$cmd1`;
	my $cmd2 = "perl $Bin/../json_output_nodb/genes_json.pl project=$projectName" . $args;
	print " -> Check DEV environment\n";
	my $jsonExp = `$cmd2`;
	Test_F_utils::launchTestF_compareTwoJsonText($jsonObs, $jsonExp, $testName, $debug);
}

