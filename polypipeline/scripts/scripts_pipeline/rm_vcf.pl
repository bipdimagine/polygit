#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use GBuffer;
use Data::Dumper;
use Getopt::Long;
use Carp;
use Logfile::Rotate; 
use File::Basename;


my $vcf_file = $ARGV[0];


die("$vcf_file not found ") unless -e $vcf_file;

backup($vcf_file);
exit(0);



sub backup {
	my ($final_gz) = @_;
	my $dir =  dirname($final_gz);
	my $dirb = $dir."/backup";
	unless (-e $dirb){
		mkdir $dirb unless -e $dirb;
		system("chmod a+w $dirb");
	}
	if (-e $final_gz){
	my $log = new Logfile::Rotate( File => $final_gz,	
	 								Count => 5,
	 								Gzip  => 'no',
	 								 Flock  => 'no',
	 								 Dir => $dirb
	 								);
		$log->rotate();	 								
		}
		unlink $final_gz;
		unlink $final_gz.".tbi" if -e $final_gz.".tbi";
}