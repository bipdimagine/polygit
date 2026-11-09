#! /usr/bin/env perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use GBuffer; 
use Getopt::Long;
use Carp;
use autodie qw(system open);
use File::Basename;
use JSON::PP;

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
my $bam_in			= $project->getAlignmentPipelineDir('isoseq_refine') . $patient->name.'.flnc.bam';
confess ("Error: Input bam '$bam_in' not found") unless (-f $bam_in);

my $dir_out			= $project->getAlignmentPipelineDir('isoseq_correct');
my $bam_out			= $dir_out . $patient->name.'.flnc-bccorr.bam';
my $report_file		= $dir_out . $patient->name.'.flnc-bccorr.report.json';

if (-e  $bam_out){
	if ($force){
		unlink $bam_out;
		unlink $bam_out.'pbi';
	}
	else {
		warn "Already done: $bam_out";
		check_report($report_file);
		exit(0);
	}
}

my $path_cellranger = `realpath \$(which cellranger)`;
chomp $path_cellranger;
$path_cellranger	=~ s/bin\/cellranger//;
warn ($path_cellranger);

# 3M-3pgex-may-2023_TRU.txt.gz -> 3'GEM-X -> SC3Pv4
# 3M-febrary-2018.txt.gz -> 3' NextGEM -> SC3Pv3
# For Poly-A capture, the v9 pipeline will use the inclusion list, 3M-3pgex-may-2023_TRU.txt.gz or 3M-febrary-2018_TRU.txt.gz.
# For Capture Sequence 1, the v9 pipeline will use the barcode inclusion list file, 3M-3pgex-may-2023_NXT.txt.gz and 3M-february-2018_NXT.txt.gz,

# https://kb.10xgenomics.com/s/article/115004506263-What-is-a-barcode-inclusion-list-formerly-barcode-whitelist
my $kit = $patient->getChemistry;
my $bc10x			= $path_cellranger."/lib/python/cellranger/barcodes/3M-5pgex-jan-2023.txt.gz" if ($kit =~ /^SC5P.*v3$/);
$bc10x				= $path_cellranger."/lib/python/cellranger/barcodes/3M-3pgex-may-2023_TRU.txt.gz" if ($kit eq 'SC3Pv4');
$bc10x				= $path_cellranger."/lib/python/cellranger/barcodes/3M-february-2018_TRU.txt.gz" if ($kit =~ /^SC3Pv3/);

# Rev comp des barcodes
if ($kit =~ /^SC3P/) {
	my $bc10x_fwd = $bc10x;
	open(my $in, "gunzip -c $bc10x_fwd |");
	my $t = time;
	$bc10x = $dir_out . basename($bc10x_fwd)=~s/\.txt\.gz$/_RC.$t.txt.gz/r; #'3M-3pgex-may-2023_TRU_RC.txt.gz';
	open(my $out, "| gzip > $bc10x");
	while (my $line = <$in>) {
		chomp($line);
		my $rc = reverse($line);
		$rc =~ tr/ACGTNacgtn/TGCANtgcan/;
		print $out "$rc\n";
	}
	close($in);
	close($out);
}
warn $bc10x;
die("No cellranger barcode list found for chemistry '$kit'") unless (-f $bc10x);

my $cmd = qq{cd $dir_out && exec_singularity.sh lima-isoseq.sif isoseq correct --barcodes $bc10x $bam_in $bam_out -j $fork};
print $cmd."\n";
system($cmd);
die("Error running isoseq correct: no output '$bam_out'") unless (-e $bam_out);


check_report($report_file);

exit(0);



sub check_report {
	my ($json_file, $min_yield_fraction) = @_;
	$min_yield_fraction //= 0.95; # Seuil par défaut : 85% (0.85)
	confess("Error: no output JSON report file '$json_file'") unless (-e $json_file);

	open(my $fh, '<', $json_file) or die "Can't open '$json_file' : $!";
	my $json_text = do { local $/; <$fh> };
	close($fh);
	my $data = eval { decode_json($json_text) };
	confess("Error: failed to parse JSON in '$json_file': $@") if $@;
	my $yield_fraction;
	if (exists $data->{attributes} && ref($data->{attributes}) eq 'ARRAY') {
		for my $attr (@{$data->{attributes}}) {
			if ($attr->{id} eq 'yield_fraction') {
				$yield_fraction = $attr->{value};
			}
		}
	}
	confess("Error: Can't extract 'yield_fraction' in '$json_file'") unless (defined $yield_fraction);

	warn(sprintf("IsoSeq Correct summary: %.3f yield fraction", $yield_fraction));
	confess(sprintf("Error control: IsoSeq Correct yield (%.3f) is below threshold (%.3f).", $yield_fraction, $min_yield_fraction)) if ($yield_fraction < $min_yield_fraction);
}


