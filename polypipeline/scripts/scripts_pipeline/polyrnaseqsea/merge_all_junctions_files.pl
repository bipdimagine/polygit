#!/usr/bin/perl 
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use Parallel::ForkManager;
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
use Getopt::Long;
use Carp;
use GBuffer;



my ($project_name);
GetOptions(
	'project=s'  => \$project_name,
);

confess("\n\nERROR: -project option mandatory. Die...\n\n") unless $project_name;


my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $project_name );

my ($dirout_method, @l_RI, @l_SE);
foreach my $patient (@{$project->getPatients()}) {
	my $align_method = $patient->alignmentMethods->[0];
	$dirout_method = $project->project_path."junctions/".$align_method."/rnaseqsea";
	my $fileRI = $dirout_method.'/'.$patient->name().'_RI.txt.gz';
	my $fileSE = $dirout_method.'/'.$patient->name().'_SE.txt.gz';
	confess("\n\nERROR: file $fileRI not found. Die...\n\n") unless -e $fileRI;
	confess("\n\nERROR: file $fileSE not found. Die...\n\n") unless -e $fileSE;
	push(@l_RI, $fileRI);
	push(@l_SE, $fileSE);
}


my $dirout_analyse = $project->project_path . "/analysis/AllRes/";
unless (-d $dirout_analyse) {
	`mkdir $dirout_analyse`;
	`chmod 775 $dirout_analyse`;
}
my $fileout_analyse_RI = $dirout_analyse."/allResRI.txt";
my $fileout_analyse_SE = $dirout_analyse."/allResSE.txt";
my $fileout_analyse_RI_gz = $fileout_analyse_RI.".gz";
my $fileout_analyse_SE_gz = $fileout_analyse_SE.".gz";

`rm $fileout_analyse_RI` if -e $fileout_analyse_RI;
`rm $fileout_analyse_SE` if -e $fileout_analyse_SE;
`rm $fileout_analyse_RI_gz` if -e $fileout_analyse_RI_gz;
`rm $fileout_analyse_SE_gz` if -e $fileout_analyse_SE_gz;
`rm $fileout_analyse_RI_gz.tbi` if -e $fileout_analyse_RI_gz.'.tbi';
`rm $fileout_analyse_SE_gz.tbi` if -e $fileout_analyse_SE_gz.'.tbi';

my $first_RI = $l_RI[0];
my $cmd1 = "zmore $first_RI | grep '#' >$fileout_analyse_RI";
`$cmd1`;

foreach my $f (@l_RI) {
	my $cmd_m1 = "zgrep -v '#' $f >>$fileout_analyse_RI";
	`$cmd_m1`;	
}

my @l_cmd_gz1;
push(@l_cmd_gz1, "(grep '^#' $fileout_analyse_RI; grep -v '^#' $fileout_analyse_RI");
push(@l_cmd_gz1, 'sort -t"`printf \'\t\'`" -k4,4 -k5,5n)');
push(@l_cmd_gz1, "bgzip > $fileout_analyse_RI_gz;");
my $cmd_gz1 = join(' | ', @l_cmd_gz1);
`$cmd_gz1`;

my $cmd_t1 = "tabix -s 4 -b 5 $fileout_analyse_RI_gz";
`$cmd_t1`;

my $first_SE = $l_SE[0];
my $cmd2 = "zmore $first_SE | grep '#' >$fileout_analyse_SE";
`$cmd2`;

foreach my $f (@l_SE) {
	my $cmd_m2 = "zgrep -v '#' $f >>$fileout_analyse_SE";
	`$cmd_m2`;	
}

my @l_cmd_gz2;
push(@l_cmd_gz2, "(grep '^#' $fileout_analyse_SE; grep -v '^#' $fileout_analyse_SE");
push(@l_cmd_gz2, 'sort -t"`printf \'\t\'`" -k5,5 -k6,6n)');
push(@l_cmd_gz2, "bgzip > $fileout_analyse_SE_gz;");
my $cmd_gz2 = join(' | ', @l_cmd_gz2);
`$cmd_gz2`;

my $cmd_t2 = "tabix -s 5 -b 6 $fileout_analyse_SE_gz";
`$cmd_t2`;

`rm $fileout_analyse_RI` if -e $fileout_analyse_RI;
`rm $fileout_analyse_SE` if -e $fileout_analyse_SE;

my $cmd_cp = "cp $fileout_analyse_RI_gz* $fileout_analyse_SE* $dirout_method/.";
`$cmd_cp`;

