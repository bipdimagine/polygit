#!/usr/bin/perl
use FindBin qw($Bin $RealBin);
use strict;

use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/"; 
use GBuffer ;
use Data::Dumper;
use Getopt::Long;
use Carp;
use Bio::DB::Sam;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use List::Util qw(sum);
use File::Temp;
 use Time::Elapsed qw( elapsed );
 use Time::ETA;
use Storable qw(store retrieve freeze);

my $filein;
my $dir;
my $file_bed;
 
my $dp_limit = 5;
my $al = 3;
my $project_name;
my $fork;
my $patient_name;
$| =1;
my $log_file;
my $vcf_final;
my $region;
my $start;
my $end;
my $chr;
my $padding;
my $fileout;
my $bamin;
my $gvcf;
my $bamout;
GetOptions(
	'project=s'   => \$project_name,
	"patients=s" => \$patient_name,
	"padding=s" => \$padding,
	"bamin=s" =>\$bamin,
	"bamout=s" => \$bamout,
	"gvcf=s" => \$gvcf,
	"fork=s" => \$fork,
	);
die() unless $fork;
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatient($patient_name);
my $elprep = $buffer->software("elprep5");#$project->getSoftware('java');	
$bamin = $patient->getBamFile() unless $bamin;
my $type_elprep = "filter";
my $tmp_dir = "";
my $target_bed;

my $dir_prod  = $patient->project->getGvcfDir("haplotypecaller");
my $final_gvcf = $dir_prod."/".$patient->name.".g.vcf.gz";
unlink $final_gvcf if -e $final_gvcf;
my $tmpdir = $project->getCallingPipelineDir("elprep5_gvcf");
$tmpdir.= "/".$patient->name.".".time."/";
unless (-e $tmpdir){
		system("mkdir -p $tmpdir; chmod a+rwx $tmpdir");
}
my $bam2 = $tmpdir."/".$patient->name.".out.bam";
my $dirout= $project->getCallingPipelineDir("elprep5_gvcf");
my $window_length = 5_000_000;
unless ($project->isGenome){
	my @abed ;
	
	foreach my $chr (@{$project->getChromosomes}){
		#next if $chr->name ne 21; 
		my $windows = $chr->getWindowCaptureForCalling(250,$window_length);
		my $intspan = Set::IntSpan::Fast::XS->new();
		 foreach my $window (@$windows){
		 	$intspan = $intspan->union($window->{intspan});
		 }
		#my $intspan = $chr->getIntSpanCaptureForCalling($padding);
		push (@abed, $buffer->intspanToBed($chr,$intspan));
	}
	
	my $bedfile= $tmpdir."/".$patient->name.".bed";
	open(BED,">$bedfile");
	print BED join("\n",@abed);
	close(BED);
	system("bgzip $bedfile");
	$bedfile = $bedfile.".gz";
	system("tabix -p bed  $bedfile");
	
	$target_bed = " --target-regions $bedfile  ";
}
 unlink $final_gvcf if -e $final_gvcf;
my $ref =  $project->dirGenome().$project->buffer->index("elprep");
my $cmd = qq { $elprep filter $bamin $bamout --nr-of-threads $fork -mark-duplicates --sorting-order coordinate  --reference $ref   --haplotypecaller  $final_gvcf  $target_bed};
;
my $t = time; 
system($cmd);
my $tabix = $buffer->software("tabix");
system("$tabix -p vcf $final_gvcf");
warn "TIME : ".abs(time - $t);
 