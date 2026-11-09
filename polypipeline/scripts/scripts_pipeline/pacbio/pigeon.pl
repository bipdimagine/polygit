#! /usr/bin/env perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use GBuffer; 
use Getopt::Long;
use Carp;
use autodie qw(system);


my $project_name;
my $patient_name;
my $force;
my $fork = 0;
GetOptions(
	'project=s'	=> \$project_name,
	'patient=s'	=> \$patient_name,
	'fork=i'	=> \$fork,
	'force'		=> \$force,
) || confess("Error in command line arguments\n");

die("Error: -project is mandatory") unless ($project_name);
die("Error: -patient is mandatory") unless ($patient_name);

my $buffer			= GBuffer->new();
my $project			= $buffer->newProject( -name => $project_name );
my $patient			= $project->getPatient($patient_name);
my $profile 		= $patient->getSampleProfile;
my $capture_name	= $patient->getCapture->name;
#my $kit = $patient->getChemistry;

my $dedup_fasta 	= $project->getAlignmentPipelineDir('isoseq_groupdedup') .  $patient_name .'.dedup.fasta';
confess ("Error: isoseq groupdedup fasta '$dedup_fasta' not found") unless (-f $dedup_fasta);
my $dir_collapse	= $project->getCallingPipelineDir('isoseq_collapse');
my $collpase_gff	= $dir_collapse . $patient_name .'.collapsed.gff';
confess ("Error: Input gff '$collpase_gff' not found") unless (-f $collpase_gff);
my $abundance		= $dir_collapse . $patient_name .'.collapsed.abundance.txt';
confess ("Error: isoseq collapse abundance '$abundance' not found") unless (-f $abundance);
my $flnc_count		= $dir_collapse . $patient_name .'.collapsed.flnc_count.txt';
confess ("Error: isoseq collapse fasta '$flnc_count' not found") unless (-f $flnc_count);
my $collapse_group = $dir_collapse .  $patient_name .'.collapsed.group.txt';
confess ("Error: isoseq collapse group text '$collapse_group' not found") unless (-f $collapse_group);

my $dir_tmp = $project->getCallingPipelineDir('pigeon');
my $dir_out = $project->getCountingDir('pigeon');
my $out_make_seurat = $dir_out . 'make-seurat/' . $patient_name . '/filtered/isoforms_seurat/matrix.mtx';
if (-e  $out_make_seurat){
	if ($force){
		system("rm -r $dir_out$patient_name* ${dir_out}make-seurat/$patient_name");
	}
	else {
		warn "Already done: $out_make_seurat";
		exit(0);
	}
}

my $ref_fasta 		= $project->genomeFasta;
my $gtf_annotation 	= $project ->dirGenome . 'genes/genes.gtf';
$gtf_annotation 	= $project->gtf_file() unless (-f $gtf_annotation);
confess("Error: GTF file not found: $gtf_annotation") unless (-e $gtf_annotation);




#---------------
# pigeon prepare
#---------------

warn ("\n".'-'x 20 ,"\npigeon prepare\n".'-'x 20,"\n\n");

my $collpase_gff_sorted = $collpase_gff =~ s/\.(\w+)$/\.sorted.$1/r || confess();
warn $collpase_gff_sorted;
my $gtf_annotation_sorted = $gtf_annotation =~ s/\.(\w+)$/\.sorted.$1/r || confess();
warn $gtf_annotation_sorted;
unless (-f $collpase_gff_sorted and -f $gtf_annotation_sorted) {
	# pigeon prepare <gencode.annotation.gtf> <reference.fa> <cage.bed> <intropolis.tsv> <collapsed.gff>
	my $cmd_prepare = qq{exec_singularity.sh smrttools.sif pigeon prepare };
	$cmd_prepare .= "$collpase_gff " unless (-e $collpase_gff_sorted);
	$cmd_prepare .= "$gtf_annotation " unless (-e $gtf_annotation_sorted);
	warn $cmd_prepare."\n";
	system($cmd_prepare);
}
die("Error running pigeon prepare: no outputs '$collpase_gff_sorted' and/or '$gtf_annotation_sorted'") unless (-f $collpase_gff_sorted and -f $gtf_annotation_sorted);
$collpase_gff = $collpase_gff_sorted;
$gtf_annotation = $gtf_annotation_sorted;


#---------------
# pigeon classify
#---------------

warn ("\n".'-'x 20 ,"\npigeon classify\n".'-'x 20,"\n\n");

#my $dir_classify = $project->getCallingPipelineDir('pigeon_classify');

my $cmd_classify;
if ($profile =~ /^singlecell/ or $capture_name eq 'transcriptome_10X') {	# SC
	# pigeon classify <collapsed.sorted.gff> <annotations.gtf> <reference.fa> --fl abundance.txt 	# sc
	$cmd_classify = qq{cd $dir_tmp && exec_singularity.sh smrttools.sif pigeon classify $collpase_gff $gtf_annotation $ref_fasta --fl $abundance --out-dir $dir_out -j $fork };
}
else {	# Bulk
	# pigeon classify <collapsed.sorted.gff> <annotations.gtf> <reference.fa> --fl flnc_count.txt 	# bulk
	$cmd_classify = qq{cd $dir_tmp && exec_singularity.sh smrttools.sif pigeon classify $collpase_gff $gtf_annotation $ref_fasta --fl $flnc_count --out-dir $dir_out -j $fork };
}
warn $cmd_classify."\n";
system($cmd_classify);

my $classification = $dir_out . $patient_name .'_classification.txt';
die("Error running pigeon classify: no output '$classification'") unless (-e $classification);


#---------------
# pigeon filter
#---------------

warn ("\n".'-'x 20 ,"\npigeon filter\n".'-'x 20,"\n\n");

# pigeon filter <classification.txt> (--isoforms <collapsed.sorted.gff>)
my $cmd_filter = qq{exec_singularity.sh smrttools.sif pigeon filter $classification --isoforms $collpase_gff -j $fork };
warn $cmd_filter."\n";
system($cmd_filter);

#my $classification_filtered = $dir_out . $patient_name . '_classification.filtered_lite_classification.txt';
my $classification_filtered = $classification =~ s/_classification\.txt$/_classification.filtered_lite_classification.txt/r || confess ("Error getting filtered classification file name");
die("Error running pigeon filter: no output '$classification_filtered'") unless (-e $classification_filtered);


#---------------
# pigeon report
#---------------

warn ("\n".'-'x 20 ,"\npigeon report\n".'-'x 20,"\n\n");

my $dir_report = $project->getCallingPipelineDir('pigeon_report');
#my $saturation = $dir_report . $patient_name .'.saturation.txt';
my $saturation = $dir_out . $patient_name .'.saturation.txt';

my $cmd_report;
if ($profile =~ /^singlecell/ or $capture_name eq 'transcriptome_10X') {	# SC
	# pigeon report <classification.filtered_lite_classification.txt> <saturation.txt> 	# sc
	$cmd_report = qq{exec_singularity.sh smrttools.sif pigeon report $classification_filtered $saturation -j $fork };
}
else {	# Bulk
	# pigeon report --exclude-singletons <classification.filtered_lite_classification.txt> <saturation.txt> 	# bulk
	$cmd_report = qq{exec_singularity.sh smrttools.sif pigeon report --exclude-singletons $classification_filtered $saturation -j $fork };
}
warn $cmd_report."\n";
system($cmd_report);

die("Error running pigeon report: no output '$saturation'") unless (-e $saturation);


#---------------
# pigeon make-seurat
#---------------

warn ("\n".'-'x 20 ,"\npigeon make-seurat\n".'-'x 20,"\n\n");

my $dir_make_seurat_raw = $dir_out . 'make-seurat/' . $patient_name . '/raw/';
system("mkdir -p $dir_make_seurat_raw --mode 775") unless (-d $dir_make_seurat_raw);
my $dir_make_seurat_filtered = $dir_out . 'make-seurat/' . $patient_name . '/filtered/';
system("mkdir $dir_make_seurat_filtered --mode 775") unless (-d $dir_make_seurat_filtered);

# pigeon make-seurat --dedup <dedup.fasta> --group <collapse.group.txt> -d <output_dir> <classification.filtered_lite_classification.txt>
my $cmd_make_seurat = qq{exec_singularity.sh smrttools.sif pigeon make-seurat --dedup $dedup_fasta --group $collapse_group -d $dir_make_seurat_raw $classification --out-prefix $patient_name -j $fork };
warn $cmd_make_seurat."\n";
system($cmd_make_seurat);

warn "\nfiltered\n";
my $cmd_make_seurat_filtered = qq{exec_singularity.sh smrttools.sif pigeon make-seurat --dedup $dedup_fasta --group $collapse_group -d $dir_make_seurat_filtered $classification_filtered --out-prefix $patient_name -j $fork };
warn $cmd_make_seurat_filtered."\n";
system($cmd_make_seurat_filtered);

#my $out_make_seurat = $dir_make_seurat . 'isoforms_seurat/matrix.mtx';
die("Error running pigeon make-seurat: no output '$out_make_seurat'") unless (-e $out_make_seurat);


exit(0);


