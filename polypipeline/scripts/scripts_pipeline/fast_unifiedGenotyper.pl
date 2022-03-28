#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";

use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use colored;
use calling;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
my $filein;
my $dir;
my $file_bed;
my $project_name;
my $chr_name;
my $fork = 1;
use List::MoreUtils qw(part);


$| = 1;

my $recal_ext = qq{};
my $log_file;
my $end_ext = "uni";
my $patient_name;
my $dir_out;
my $remove_par_regions;
my $keep;
my @temp_file;
my $somatic = undef;

GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"chr=s"	=> \$chr_name,
	"log=s" =>\$log_file,
	"ext=s" => \$end_ext,
	"patient=s" => \$patient_name,
	"dir_out=s" => \$dir_out,
	"keep=s" => \$keep,
	"somatic=s" => \$somatic,
	"remove_par_regions=s" => \$remove_par_regions
);
my $date = `date`;
chomp($date);

if ($log_file){
	open (STDOUT,">>".$log_file);
}
$SIG{INT} = \&interrupt;



	
my $other_project = [];
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $chr = $project->getChromosome($chr_name);
my $patients =  $project->get_list_patients($patient_name);
#"patients" => $patients, "fork" => $fork ,
my $final =  $project->getCallingPipelineDir("unifiedgenotyper")."/".$chr_name.".$end_ext.vcf";
#my $vcf = calling::calling_merge  ( $project_name ,$chr_name,$patient_name,$fork,{freebayes=>1,unifiedgenotyper=>2,samtools=>3});
$dir_out= $project->getCallingPipelineDir("unifiedgenotyper");
$dir_out .="/all/";
unless (-e $dir_out){
	system("mkdir -p $dir_out;chmod a+rwx $dir_out");
}
else {
system ("find $dir_out/$chr_name -name TMP.* -exec rm {} \\;");
}
my $vcf = calling::calling_merge  ( $project_name ,$chr_name,$patient_name,$fork,{unifiedgenotyper=>1,freebayes=>2},$dir_out,$somatic);

rename $vcf ,$final;
system ("find $dir_out/$chr_name -name TMP.* -exec rm {} \\;");
exit(0);
