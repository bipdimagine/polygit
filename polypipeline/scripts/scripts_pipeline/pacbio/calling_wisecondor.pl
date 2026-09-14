#!/usr/bin/env perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use GBuffer; 
use Data::Dumper;
use Getopt::Long;
use Carp;
use Storable qw(store retrieve freeze);
use String::ProgressBar;
use List::Util qw(sum);


my $project_name;
my $fork;
my $patient_name;
my $force;
 
my $fork = 5;
GetOptions(
	'project=s'   => \$project_name,
	"patient=s" => \$patient_name,
	"fork=s" => \$fork,
	"force=s" =>\$force,
);
die("miss fork") unless $fork;


my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );

my $patient = $project->getPatient($patient_name);
my $bam = $patient->getBamFile() ;
my $npz =  $patient->fileWiseCondor();
my $run = $patient->getRun();
my $ref = $project->get_wisecondor_reference;
my $wise_sif  = "wisecondor.sif";
my $singularity = "run_singularity.sh";
my $wise  = "$singularity $wise_sif WisecondorX";
#warn $ref;
#die();
#5kb.nova.npz
#$ref = "/data-isilon/public-data/repository/HG38/wisecondor/novaseqx/5kb/5kb.npz";
my $blfile = $project->public_data_root . "/". $project->annotation_genome_version . "/wisecondor/blacklist.bed";
#my $blfile = "/data-pure/public-data/blacklist/blacklist.bed";
#$blfile = "/data-pure/public-data/blacklist/blacklist.spectre.bed";

$blfile = $project->buffer()->config_path("root","public_data")."/repository/HG38//blacklist/encode.blacklist.1.bed";
my $blacklist = "";
$blacklist = "--blacklist ".$blfile if -e $blfile;
# warn $blacklist;
#'"/data-beegfs/npz/reference/ref5Kb.npz";


my $dir = $project->buffer()->config_path("root","public_data").'/repository/HG38/wisecondor-ref/new_version/';
if ($run->machine eq "REVIO" or $run->machine eq "SEQUEL"){
	$dir .= "revio/";
}
else {
	$dir .= "novaseq/"
}
$ref = $dir."reference.5k.npz";
warn $ref;
my $out = $project->getCallingPipelineDir("wiseCondor")."/".$patient->name;
my $cnd_wise ="exec_singularity.sh wisecondor.gemini.sif wisecondorx ";
my $outbed1 = $out."_bins.bed";
my $prod_file1 = $project->getVariationsDir("wisecondor")."/".$patient->name."_bins.bed.gz";
my $outbed2 = $out."_aberrations.bed";
my $prod_file2 = $project->getVariationsDir("wisecondor")."/".$patient->name."_aberrations.bed.gz";
if (-e $prod_file2){
	if ($force){
		unlink $prod_file2;
	}
	else {
		warn "already done: $prod_file2";
		exit(0);
	}
}

my $cmd = "$cnd_wise predict  $npz $ref $out  --beta 1 --bed  ".$blacklist;
warn $cmd;
system("$cmd");
my $bgzip = "bgzip";#$buffer->software("bgzip");
my $tabix = "tabix";#$buffer->software("tabix");
system("$bgzip $outbed1 && mv $outbed1.gz $prod_file1 ; $tabix -f -p bed -S 1 $prod_file1 ");
system("$bgzip $outbed2 && mv $outbed2.gz $prod_file2 ; $tabix -f -p bed -S 1 $prod_file2 ");
warn $ref;
exit(0);
 