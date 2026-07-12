#!/usr/bin/env perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
#use Set::IntSpan;
use GBuffer; 
use Data::Dumper;
use Getopt::Long;
use Carp;
use Storable qw(store retrieve freeze);
use String::ProgressBar;
use List::Util qw(sum);


 my $project_name;
 my $fork;
 my $callable_intspan_file;
 my $patient_name;
 #my $low_calling;
 my $method;
 
 my $chr_syno ={
		1=> "chr1",
		2=>"chr2",
		3=>"chr3",
		4=>"chr4",
		5=>"chr5",
		6=>"chr6",
		7=>"chr7",
		8=>"chr8",
		9=>"chr9",
		10=>"chr10",
		11=>"chr11",
		12=>"chr12",
		13=>"chr13",
		14=>"chr14",
		15=>"chr15",
		16=>"chr16",
		17=>"chr17",
		18=>"chr18",
		19=>"chr19",
		20=>"chr20",
		21=>"chr21",
		22=>"chr22",
		X=>"chrX",
		Y=>"chrY",
		MT=>"chrM",
};

 
my $fork = 5;
GetOptions(
	'project=s'   => \$project_name,
	"patient=s" => \$patient_name,
	"fork=s" => \$fork,
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
 
 
 my $repository = $project->buffer()->config_path("root","public_data").'/repository/HG38/wisecondor-ref/new_version/novaseq';
my $dir = $project->buffer()->config_path("root","public_data").'/repository/HG38/wisecondor-ref/new_version/';
#$ref = $repository ."/reference.10000.default.npz";
if ($run->machine eq "REVIO" or $run->machine eq "SEQUEL"){
	$dir .= "revio/";
}
else {
	$dir .= "novaseq/"
}
 $ref = $dir."reference.10k.npz";
warn $ref;
 my $out = $project->getCallingPipelineDir("wiseCondor")."/".$patient->name;
 my $cnd_wise ="exec_singularity.sh wisecondor.gemini.sif wisecondorx ";
 
 my $cmd = "$cnd_wise predict  $npz $ref $out  --beta 1 --bed  ".$blacklist;
 warn $cmd;
 system("$cmd");
 my $bgzip = "bgzip";#$buffer->software("bgzip");
 my $tabix = "tabix";#$buffer->software("tabix");
my $outbed1 = $out."_bins.bed";

my $prod_file = $project->getVariationsDir("wisecondor")."/".$patient->name."_bins.bed.gz";
system("$bgzip $outbed1 && mv $outbed1.gz $prod_file ; $tabix -f -p bed -S 1 $prod_file ");

$outbed1 = $out."_aberrations.bed";
$prod_file = $project->getVariationsDir("wisecondor")."/".$patient->name."_aberrations.bed.gz";
system("$bgzip $outbed1 && mv $outbed1.gz $prod_file ; $tabix -f -p bed -S 1 $prod_file ");
warn $ref;
 exit(0);
 