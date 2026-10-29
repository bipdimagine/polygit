#! /usr/bin/env perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
#use Set::IntSpan;
use GBuffer; 
use Getopt::Long;
use IO::Uncompress::Gunzip;
use Carp;
use autodie qw(system open);


my $method;
my $project_name;
my $force;
GetOptions(
	"project=s" => \$project_name,
	"force=s" 	=> \$force,
) || confess("Error in command line arguments\n");

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );

foreach my $patient (@{$project->getPatients}){
	warn $patient->name;
	my $run = $patient->getRun();
	my $run_name = $run->plateform_run_name;
	my $bcl_dir = $run->bcl_dir."/no_sample_id/*/fastq_pass/*";
	#system("mkdir $bcl_dir") unless -e bcl_tmp;
	my $fastq_dir = $run->fastq_dir();
	system("mkdir $fastq_dir") unless -e $fastq_dir;
	
	my $bc = $patient->barcode();
	warn qq{find  /$bcl_dir -name "*pass*${bc}_*gz"};
	my @files =  `find  /$bcl_dir -name "*pass*${bc}_*gz"`;
	chomp(@files);
	die("no file for bc $bc ".$patient->name) unless @files;
	
	#my $model = basecall_model($files[0]);
	my $final_fastq = $fastq_dir."/".$patient->name.".fastq.gz";
	my $ubam_out = $fastq_dir."/".$patient->name.".bam";
	if ( -e $ubam_out ) {
		if ($force) {
			unlink $ubam_out;
		}
		else {
			warn "Already done: $ubam_out";
			warn "\n";
			next;
		}
	}
	
	system("cat ".join(" ",@files)." >".$final_fastq);
	warn $final_fastq;
	system("samtools import -0 $final_fastq -@ 20 -O BAM -o ${ubam_out} && rm $final_fastq");
	warn $ubam_out;
	warn "\n";
}

exit(0);



sub basecall_model {
	my ($file) = @_;
	
	my $gz = IO::Uncompress::Gunzip->new($file)
    or die "Impossible d'ouvrir $file: \n";

my $model;
while (my $line = <$gz>) {

    next unless $line =~ /^\@/;

    if ($line =~ /basecall_model_version_id=([^\s]+)/) {
		$model = $1;
        print "Model: $1\n";

        last;

    }

}
 if ($model =~ /sup/ ){
		return "sup";
	}
	die("-".$model."- define model");
	return $model;
}