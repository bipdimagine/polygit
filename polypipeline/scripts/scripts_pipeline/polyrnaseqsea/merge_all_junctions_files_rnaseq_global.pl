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


my $dirout_analyse = $project->project_path . "/analysis/AllRes/";
my $dirin_analyse = $project->project_path . "/analysis/AllRes/resPat/";
unless (-d $dirin_analyse) {
	die("\n\nERROR $dirout_analyse not found. Die\n\n");
}


`gunzip $dirin_analyse/allResRI.txt.gz`;
`gunzip $dirin_analyse/allResSE.txt.gz`;

open (OUT, ">$dirout_analyse/allResRI.txt");
print OUT "#Junc_RI\tENSID\tGene\tChr\tJunc_RI_Start\tJunc_RI_End\tJunc_RI_Count\tJunc_Normale\tJunc_Normale_Start\tJunc_Normale_End\tJunc_Normale_Count\tType\tscore\tminMoyNcountParJunc\tmaxMoyNcountParJunc\tSample\n";
open (RES, "$dirin_analyse/allResRI.txt");
while (<RES>) {
	my $line = $_;
	my @lcol = split("\t", $line);
	print OUT $line if ($lcol[0] =~ /_/ and $lcol[1] =~ /ENSG/);
}
close(RES);
close (OUT);

my @l_cmd_gz1;
push(@l_cmd_gz1, "(grep '^#' $dirout_analyse/allResRI.txt; grep -v '^#' $dirout_analyse/allResRI.txt");
push(@l_cmd_gz1, 'sort -t"`printf \'\t\'`" -k4,4 -k5,5n)');
push(@l_cmd_gz1, "bgzip > $dirout_analyse/allResRI.txt.gz;");
my $cmd_gz1 = join(' | ', @l_cmd_gz1);
`$cmd_gz1`;
`tabix -s 4 -b 5 $dirout_analyse/allResRI.txt.gz`;



#### ------------------------------------ ####



open (OUT2, ">$dirout_analyse/allResSE.txt");
print OUT2 "#Junc_SE\tENSID\tGene\tChr\tJunc_SE_Start\tJunc_SE_End\tJunc_SE_Count\tNbreSkippedExons\tSkippedExonsID\tTranscritsID\tJunc_Normale\tJunc_Normale_Count\tscore\tType\tminMoyNcountParJunc\tmaxMoyNcountParJunc\tSample\n";
open (RES2, "$dirin_analyse/allResSE.txt");
while (<RES2>) {
	my $line = $_;
	my @lcol = split("\t", $line);
	print OUT2 $line if ($lcol[0] =~ /_/ and $lcol[1] =~ /ENSG/);
}
close(RES2);
close (OUT2);

my @l_cmd_gz2;
push(@l_cmd_gz2, "(grep '^#' $dirout_analyse/allResSE.txt; grep -v '^#' $dirout_analyse/allResSE.txt");
push(@l_cmd_gz2, 'sort -t"`printf \'\t\'`" -k4,4 -k5,5n)');
push(@l_cmd_gz2, "bgzip > $dirout_analyse/allResSE.txt.gz;");
my $cmd_gz2 = join(' | ', @l_cmd_gz2);
`$cmd_gz2`;
`tabix -s 4 -b 5 $dirout_analyse/allResSE.txt.gz`;

`rm $dirout_analyse/allResRI.txt`;
`rm $dirout_analyse/allResSE.txt`;
`rm $dirout_analyse/allResRI.rds`;
`rm $dirout_analyse/allResSE.rds`;
`rm -r $dirin_analyse`;


foreach my $dir_name ('BamJuncs', 'AllresRI', 'AllresSE', 'AllRes') {
	my $dir = $project->project_path . "/".$dir_name."/";
	`rm -r $dir` if -d $dir;
}




