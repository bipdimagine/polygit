#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin/../../lib/";
use Logfile::Rotate; 
use DBI; 
use DBD::mysql;
use Data::Dumper;
use Getopt::Long;

my $out;
my $in;
GetOptions(
		'in=s'       => \$in,
         'out=s'   => \$out,
);
my $file = "$out/poly-src.tar.gz";

system ("tar -zcv --exclude='*test*' --exclude='core*' --exclude='*svn*' -f $file $in");

 my $log = new Logfile::Rotate( File => $file,
	 								Count => 20000,
	 								Gzip  => 'lib',
	 								 Flock  => 'no',
	 								);

$log->rotate();								
unlink($file);




