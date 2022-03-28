#!/usr/bin/perl
use FindBin qw($Bin $RealBin);
use strict;

use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
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
use List::MoreUtils qw(natatime);

my $project_name;
my $patient_name;
my $fileout;
my $filein;

GetOptions(
	'project=s'   => \$project_name,
);



my $other_project = [];
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
#my $patients = $project->getPatient($patient_name);

my $xhmm =$buffer->software("xhmm");

my $xhmm_parameter = $buffer->software("xhmm_parameter");

my $capture = $project->getCapture();
my $dir_out= $project->getCallingPipelineDir("xhmm");
my $bed_file = $dir_out."/capture.bed";


	my @all;
foreach my $chr (@{$project->getChromosomes}){
		warn "start ".$chr->name;
		my $intspan = $capture->getIntSpanForChromosome($chr,50);
		push(@all,$buffer->intspanToBed($chr,$intspan));
			
}
	open(BED,">".$bed_file);
print BED join("\n",@all);
close(BED);
my $reference = $project->genomeFasta();
my $javac = $project->getSoftware('java');
my $gatk  = $project->getSoftware('gatk');
my $fttemp = $dir_out."/gatk_tmp.txt";;
my $extreme_gc = $dir_out."/extreme_gc_targets.txt";
my $cmd0 = qq{$javac  -Xmx48g -jar $gatk -T GCContentByInterval -L $bed_file -R $reference -o $fttemp && cat $fttemp | awk '{if (\$2 < 0.1 || \$2 > 0.9) print $1}' >$extreme_gc && rm $fttemp };
system($cmd0);
my $DATA_RD = $dir_out."/DATA.RD.txt";

my $cmd1 = qq{$xhmm  --mergeGATKdepths  -o $DATA_RD };
my $coverage_dir = $project->getRootDir() . "/align/coverage/depth/";

foreach my $patient (@{$project->getPatients}){
	my $file = $coverage_dir."/".$patient->name.".depth";
	$cmd1.=" --GATKdepths $file";
}
#warn $cmd1;
system($cmd1);

my $filtered_targets = $dir_out."/DATA.filtered_centered.RD.txt.filtered_targets.txt";
my $filtered_sample =  $dir_out."/DATA.filtered_centered.RD.txt.filtered_samples.txt";
my $filtered_centered =  $dir_out."/DATA.filtered_centered.RD.txt";


my $cmd2 = qq{$xhmm  --matrix -r $DATA_RD  --centerData --centerType target  -o $filtered_centered --excludeTargets $extreme_gc  --outputExcludedTargets  $filtered_targets --outputExcludedSamples $filtered_sample  --minTargetSize 10 --maxTargetSize 10000 --minMeanTargetRD 10 --maxMeanTargetRD 500 --minMeanSampleRD 25 --maxMeanSampleRD 200 --maxSdSampleRD 150};
#--excludeTargets ./extreme_gc_targets.txt --excludeTargets ./low_complexity_targets.txt \
system($cmd2);
my $pca = $dir_out."/DATA.RD_PCA";
my $cmd3 = qq{$xhmm --PCA -r $filtered_centered  --PCAfiles $pca };
warn $cmd3;
system($cmd3);

my $pca_nom = $dir_out."/DATA.PCA_normalized.txt";
my $cmd4 = qq{$xhmm --normalize -r $filtered_centered  --PCAfiles $pca  --normalizeOutput $pca_nom --PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7};
warn $cmd4;system ($cmd4);

my $sample_zscores = $dir_out."/DATA.PCA_normalized.filtered.sample_zscores.RD.txt";
my $exclude_sample_zscores = $dir_out."/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt";
my $exclude_target_zscores =$dir_out.."/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt";
my $cmd5 = qq{$xhmm --matrix -r $pca_nom --centerData --centerType sample --zScoreData -o $sample_zscores --outputExcludedTargets $exclude_target_zscores   --outputExcludedSamples $exclude_sample_zscores --maxSdTargetRD 30};
warn $cmd5;system($cmd5);


my $DATA_RD_SAME_FILTERED =  $dir_out."/DATA.same_filtered.RD.txt";
my$cmd6 =qq{$xhmm  --matrix -r $DATA_RD --excludeTargets $filtered_targets --excludeTargets $exclude_target_zscores --excludeSamples  $filtered_sample --excludeSamples $exclude_sample_zscores -o  $DATA_RD_SAME_FILTERED};
warn $cmd6;system ($cmd6);
my $final_dir = $project->getVariationsDir("xhmm");
my $xcnv = $final_dir."/".$project->name.".xcnv";
my $aux = $dir_out.'/DATA.aux_xcnv';
my $DATA = $dir_out.'/DATA';
my $cmd7 =qq{$xhmm --discover -p $xhmm_parameter -r $sample_zscores -R  $DATA_RD_SAME_FILTERED -c $xcnv -a $aux -s $DATA };
warn $cmd7;
system($cmd7);
warn $xcnv;


my $vcf = $final_dir."/".$project->name.".vcf";

my $cmd8 = qq{$xhmm --genotype -p $xhmm_parameter -r $sample_zscores -R  $DATA_RD_SAME_FILTERED  -g $xcnv -F $reference -v $vcf};
warn $cmd8;
system($cmd8);
die() unless -e $vcf;
system "rm $dir_out/* && rmdir $dir_out";
warn $vcf;