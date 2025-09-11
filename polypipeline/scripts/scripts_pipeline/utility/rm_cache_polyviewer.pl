#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
#use lib "/software/polyweb/poly-disk/poly-src/GenBo/lib/obj-nodb/packages";
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/packages";
use GBuffer;
use Set::IntSpan::Fast::XS;
use Getopt::Long;
use JSON::XS;
use File::Copy;
use POSIX qw(strftime);

my $fork = 1;
my ($project_name, $patient_name,$vid);
my $file;

GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
);

my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => $project_name);

my $dir = $project->getCacheDir();


# Récupère tous les fichiers *.cache
my @files = glob("$dir/*.cache");

# Supprime-les
for my $file (@files) {
	
    if (-f $file) {   # vérifie que c'est bien un fichier
     	my $date = strftime("%Y%m%d", localtime);
    	my $rotated = "$file.$date";
    	system ("mv $file $rotated ");
    	system("gzip", $rotated) == 0 or warn "Impossible de compresser $rotated: $!";
    }
}


exit(0);