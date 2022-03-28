#!/usr/bin/perl
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
use Getopt::Long;
use GBuffer;
use GenBoProject;

my $fork = 1;
my ($project_name, $chr_name);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'chr=s'        => \$chr_name,
);

unless ($project_name) { die("\n\nERROR: -project option missing... Die...\n\n"); }
unless ($chr_name) { die("\n\nERROR: -chr option missing... Die...\n\n"); }


my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $project_name );
my $dir_log = $project->getCacheBitVectorDir()."/log/";


my @lChr;
if ($chr_name eq 'all') { @lChr = (1..22, 'X', 'Y', 'MT'); }
else { push(@lChr, $chr_name); }

foreach my $chr_name (@lChr) {
	warn "\n##### START CHR$chr_name #####\n";
	
	my $cmd1 = "$RealBin/scripts/cache_store_ids.pl -project=$project_name -chr=$chr_name -fork=$fork";
	eval { system("$cmd1") and die "\n\nERROR: $project_name -> cache_store_ids -> chr$chr_name. Die.\n\n"};
	if ($@) {
		die;
	};
	if (not -e $project->getCacheBitVectorDir()."/lmdb_cache/$chr_name.dv.freeze") {
		warn "\n\nNOT FOUND: ".$project->getCacheBitVectorDir().'/lmdb_cache/$chr_name.dv.freeze'."\n";
		warn "\nERROR: $project_name -> cache_store_ids -> chr$chr_name. Die.\n\n";
		die;
	}
	else{ warn "### cache_store_ids OK"; }
	
	my $cmd3 = "$RealBin/scripts/cache_store_annotations.pl -project=$project_name -chr=$chr_name -fork=$fork";
	eval { system("$cmd3") and die "\n\nERROR: $project_name -> cache_store_annotations -> chr$chr_name. Die.\n\n"};
	if ($@) {
		die;
	};
	if (not -e $project->getCacheBitVectorDir()."/lmdb_cache/$chr_name/genes_index") {
		warn "\n\nNOT FOUND: ".$project->getCacheBitVectorDir()."/lmdb_cache/$chr_name/genes_index"."\n";
		warn "\nERROR: $project_name -> cache_store_annotations -> chr$chr_name. Die.\n\n";
		die;
	}
	else{ warn "### cache_store_annotations OK"; }
	
	my $cmd4 = "$RealBin/scripts/cache_check_step.pl -project=$project_name -chr=$chr_name -fork=$fork -step=store_annotations";
	eval { system("$cmd4") and die "\n\nERROR: $project_name -> cache_check_step store_annotations -> chr$chr_name. Die.\n\n"};
	if ($@) {
		die;
	};
	if (not -e "$dir_log/check_store_annotations.$chr_name.ok") {
		warn "\n\n$dir_log/check_store_annotations.$chr_name.ok\n";
		warn "\nERROR: $project_name -> cache_check_step store_annotations -> chr$chr_name. Die.\n\n";
		die;
	}
	else{ warn "### check cache_store_annotations OK"; }
	
	my $cmd5 = "$RealBin/scripts/cache_strict_denovo.pl -project=$project_name -chr=$chr_name -fork=$fork";
	eval { system("$cmd5") and die "\n\nERROR: $project_name -> cache_strict_denovo -> chr$chr_name. Die.\n\n"};
	if ($@) {
		die;
	};
	if (not -e $project->getCacheBitVectorDir()."/strict-denovo/$chr_name.lite") {
		warn "\n\nERROR: $project_name -> cache_strict_denovo -> chr$chr_name. Die.\n\n";
		die;
	}
	else{ warn "### cache_strict_denovo OK"; }
	
	my $cmd6 = "$RealBin/scripts/cache_check_step.pl -project=$project_name -chr=$chr_name -fork=$fork -step=strict_denovo";
	eval { system("$cmd6") and die "\n\nERROR: $project_name -> cache_check_step strict_denovo -> chr$chr_name. Die.\n\n"};
	if ($@) {
		die;
	};
	if (not -e "$dir_log/check_strict_denovo.$chr_name.ok") {
		warn "\n\nERROR: $project_name -> cache_check_step strict_denovo -> chr$chr_name. Die.\n\n";
		die;
	}
	else{ warn "### check cache_strict_denovo OK"; }
	
	if ($project->isSomaticStudy()) {
		my $cmd7 = "$RealBin/scripts/cache_loh.pl -project=$project_name -chr=$chr_name -fork=$fork";
		eval { system("$cmd7") and die "\n\nERROR: $project_name -> cache_somatic_loh -> chr$chr_name. Die.\n\n"};
		if ($@) {
			die;
		};
		if (not -e $project->getCacheBitVectorDir()."/somatic_loh/$chr_name.lite") {
			warn "\n\nERROR: $project_name -> cache_somatic_loh -> chr$chr_name. Die.\n\n";
			die;
		}
		else{ warn "### loh OK"; }
		
		my $cmd8 = "$RealBin/scripts/cache_check_step.pl -project=$project_name -chr=$chr_name -fork=$fork -step=loh";
		eval { system("$cmd8") and die "\n\nERROR: $project_name -> cache_check_step loh -> chr$chr_name. Die.\n\n"};
		if ($@) {
			die;
		};
		if (not -e "$dir_log/check_loh.$chr_name.ok") {
			warn "\n\nERROR: $project_name -> cache_check_step loh -> chr$chr_name. Die.\n\n";
			die;
		}
		else{ warn "### check loh OK"; }
	}
	
	warn "\n##### FINISHED CHR$chr_name #####\n";
}


