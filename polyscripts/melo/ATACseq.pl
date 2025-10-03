#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../polypipeline/packages/";
use lib "$Bin/../../polypipeline/dragen/scripts/";
use dragen_util; 
use Data::Dumper;
use Getopt::Long;
use Carp;
use IO::Prompt;
use Net::SFTP;
use File::Find::Rule ;
use autodie;
use file_util;
use Term::Menus;


use GBuffer;
my $buffer = new GBuffer;

my $project_name = 'NGS2025_09521'; # Test ATACseq : NGS2025_09521
my $patient_name;
my @steps;
my $cpu = 20;
my $force;
my $no_exec;
GetOptions(
	'project=s'		=> \$project_name,
	'patients=s'	=> \$patient_name,
	'steps=s{1,}'	=> \@steps,
	'cpu=i'			=> \$cpu,
	'force'			=> \$force,
	'no_exec'		=> \$no_exec,
) || die("Error in command line");


my $project = $buffer->newProjectCache( -name => $project_name );
my $patients = $project->getPatients();

@steps = split(/,/, join(',',@steps));
unless (@steps) {
	my $liste_steps = ['dragen-align','flexbar','preprocessing'];
	my %Menu_1 = (
		Item_1 => {
			Text   => "]Convey[",
			Convey => $liste_steps,
		},
		Select => 'Many',
		Banner => "   Select steps:"
	);
	@steps = &Menu( \%Menu_1 );
	die if ( @steps eq ']quit[' );
}
warn 'steps='.join(',',@steps);

if (grep(/^flexbar$/i, @steps)){
	my $flexbar = $buffer->getSoftware('flexbar');
	my $output_dir = $project->getAlignmentPipelineDir('flexbar');
	warn $output_dir;
	
	open(my $fh, ">$output_dir/jobs_flexbar.txt");
	foreach my $pat (@$patients) {
		my $pname = $pat->name;
		
		my $fastq_dir = $output_dir.'input_fastq/';
		system("mkdir -p $fastq_dir;chmod a+rwx $fastq_dir") unless -e $fastq_dir;
		my ($fastq1,$fastq2) = dragen_util::get_fastq_file($pat,$fastq_dir);
	#	warn $fastq1;
	#	warn $fastq2;
		
		my $prefix = "$output_dir$pname";
		my $output_log = "$output_dir$pname.flexbar.log";
		
		my $cmd = "$flexbar \
		--reads $fastq1 \
		--reads2 $fastq2 \
		--target $prefix \
		--adapter-preset Nextera \
		--adapter-trim-end RIGHT \
		--adapter-revcomp ON \
		--adapter-revcomp-end LEFT \
		--zip-output GZ \
		--output-log $output_log \
		--threads $cpu";
		warn $cmd;
		print {$fh} $cmd."\n"; 
		
	}
	close($fh);
	system("cat $output_dir/jobs_flexbar.txt | run_cluster.pl -cpu=$cpu\n") unless ($no_exec);
}


if (grep(/dragen_align/i, @steps)){
	my $ref_dragen = $project->getGenomeIndex("dragen");
	foreach my $pat (@$patients) {
		my $pname = $pat->name;
		my $dragen_dir = $pat->getDragenDirName('pipeline')."$pname/";
		my $dragen_dir_log = $pat->getDragenDirName('pipeline')."$pname/";
		my $flexbar_dir = $project->getAlignmentPipelineDir('flexbar');
		my $runid = $pat->getRun()->id;
		my $log_pipeline = $dragen_dir_log."/".$pname.".pipeline.".time.".log"; 
		my $ok_pipeline = $dragen_dir_log."/".$pname.".ok.pipeline.".time; 
	
		my $dragen_cmd = "dragen \
		-f -r $ref_dragen \
		--output-directory $dragen_dir \
		--intermediate-results-dir /staging/tmp \
		--output-file-prefix $pname \
		-1 $flexbar_dir$pname\_1.fastq.gz \
		-2 $flexbar_dir$pname\_2.fastq.gz \
		--RGID $runid \
		--RGSM $pname \
		--enable-duplicate-marking true \
		--remove-duplicates true\
		> $log_pipeline 2> $log_pipeline.err \
		&& touch $ok_pipeline";
#		--enable-duplicate-marking true \
#		--remove-duplicates true\
		my $cmd = "$Bin/../run_dragen.pl -cmd=\"$dragen_cmd\"";
		system($cmd) unless ($no_exec);
		die("Error dragen, see $log_pipeline.err") unless (-e $ok_pipeline);
#		my $cmd_cp = "cp $dragen_dir.$pname.bam $dragen_dir.$pname.bam.bai "
		my $cmd_mv = ("perl $Bin/scripts/dragen_move.pl -project=$project_name -patient=$pname -command=align");
		warn $cmd_mv;
		system($cmd_mv) unless ($no_exec);
	}	
}


if (grep(/preprocessing/i, @steps)){
	my $cmd_cluster;
	my $output_dir = $project->getAlignmentDir('dragen-align/preprocessed');
	warn $output_dir;
	open(my $fh, ">$output_dir/jobs_preprocessing.txt");
	foreach my $pat (@$patients) {
		my $pname = $pat->name;
		my $inbam = $pat->getBamFile('dragen-align');
		my $tmp_picard = $project->getAlignmentPipelineDir('Picard-rmdup');
		my $tmp_bam = $tmp_picard.$pname.'.rmdup.bam';
		my $metrics_file = $tmp_bam =~ s/bam$/txt/r;
		my $outbam = $pat->getBamFileName('dragen-align/preprocessed');
		next if (-e $outbam and not $force);
		
		# Remove non-uniquely/low quality mapped reads, not properly mapped reads, mitochondiral reads
		my $cmd_samtools = "samtools view $inbam -q 10 -f 0x2 -e 'rname != \"chrM\"' -o $tmp_bam -@ $cpu ";
		# Remove PCR duplicates (already done by dragen)
		my $cmd_rmdup = $buffer->getSoftware('java').' '.$buffer->getSoftware('picard')." MarkDuplicates I=$tmp_bam O=$outbam M=$metrics_file REMOVE_DUPLICATES=true ";
		my $cmd_bai = "samtools index $outbam -@ $cpu ";
		print {$fh} "$cmd_samtools && $cmd_rmdup && $cmd_bai\n";
		
		# Remove alignements within problematic regions ?
	}
	system("cat $output_dir/jobs_preprocessing.txt | run_cluster.pl -cpu=$cpu") unless ($no_exec);
}


if (grep(/^qc$/i, @steps)){
	my $bam_files = join( ' ', map {$_->getBamFileName('dragen-align/preprocessed')} @$patients );
#	$bam_files =~ s/bam$/cram/;
	my $dir_bam = $project->getAlignmentDir('dragen-align/preprocessed');
	my $singularity_deeptools = "singularity exec -B $dir_bam /software/distrib/deeptools/deeptools.sif ";
	my $cmd1 = $singularity_deeptools." plotFingerprint --bamfiles $bam_files --binSize=1000 "
		."--plotFile $dir_bam/fingerprint.svg --outQualityMetrics $dir_bam/metrics.tsv --smartLabels -p $cpu";
	my $cmd2 = $singularity_deeptools."multiBamSummary bins --bamfiles $bam_files --smartLabels "
		."--outFileName $dir_bam/coverage_matrix.npz --binSize 5000 -p $cpu";
	my $cmd3 = $singularity_deeptools."plotCorrelation --corData $dir_bam/coverage_matrix.npz --plotFile $dir_bam/correlation_heatmap.svg "
		."--whatToPlot heatmap --corMethod pearson --plotNumbers "; # --outFileCorMatrix $dir_bam/correlation_bin.txt 
	my $cmd3_min = $cmd3 =~ s/correlation_heatmap/correlation_heatmap_min095/r ."-min 0.95 ";
	warn $cmd1;
	warn $cmd2;
	warn $cmd3;
	warn $cmd3_min;
	warn "\n";
	my $cmd = "$cmd1\n$cmd2 && $cmd3 && $cmd3_min";
	warn $cmd;
#	system("echo \"$cmd\" | run_cluster.pl -cpu=$cpu") unless ($no_exec);
		
}










