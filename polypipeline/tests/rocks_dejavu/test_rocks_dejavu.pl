#!/usr/bin/perl 
use strict;
use FindBin qw($Bin);
use Data::Dumper;
use RocksDB;
use Getopt::Long;
use Carp;

use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use GBuffer;
use GenBoNoSqlRocksGenome;
use rocks_dejavu_tests;

my ($verbose, $global_check, $old_path, $new_path, $chr_names);
GetOptions(
	'old=s' => \$old_path,
	'new=s' => \$new_path,
	'chr=s' => \$chr_names,
	'verbose=s' => \$verbose,
	'global_check=s' => \$global_check,
);

my $test_obj = new rocks_dejavu_tests;
$test_obj->verbose(1) if $verbose;
$test_obj->global_check() if $global_check;
$test_obj->set_directory_path_origin($old_path) if ($old_path);
$test_obj->set_directory_path_new($new_path);
if ($chr_names) {
	my @lChr = split(',', $chr_names);
	$test_obj->{use_chromosomes} = \@lChr;
}
my $is_all_ok = $test_obj->launch_test_all_chromosomes();
confess() if not $is_all_ok;
exit(0);


