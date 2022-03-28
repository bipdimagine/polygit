#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use Getopt::Long;



my $name;
#my $password;
my $dir_in;
#my $type;
my $projects;

GetOptions(
	'name=s' => \$name,
	'projects=s' => \$projects,
);

foreach my $project (split(",",$projects)){
	system("add_calling_method.sh -project=$project -method=canvas,manta,wisecondor");
#	system("../../../bds_cache.pl -project=$project -yes=1 -force=1 > $project.log");

}
