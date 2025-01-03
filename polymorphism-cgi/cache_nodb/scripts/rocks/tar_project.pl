#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use Parallel::ForkManager;
use Storable qw(store retrieve freeze thaw);
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use Cwd 'abs_path';
use Digest::MD5 qw(md5_hex);
use lib "/software/polyweb/poly-disk/poly-src/GenBo/lib/obj-nodb/";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Net::FTP;
use Carp;
use Bio::DB::HTS;
use GBuffer;
use Bit::Vector::Overload;
use Sys::Hostname;
use List::Util qw(min max sum);
use Archive::Tar;
use Digest::MD5 qw(md5 md5_hex md5_base64);
use File::Basename;

my $fork = 1;
my ($project_name, $chr_name,$annot_version);
my $ok_file;
my $patient_name;

GetOptions(
	'project=s'    => \$project_name,
	'patient=s'    => \$patient_name,
);

my $nbErrors = 0;
my $buffer = new GBuffer;

my $project = $buffer->newProject( -name =>$project_name);

my $cmd_tar = qq{tar -rf $project_name.fastq.tar -C }; 
my $md5file =  $project_name.".md5";
foreach my $p (@{$project->getPatients}){
	my $fastq = $p->fastqFiles();
	foreach my $f (@$fastq){
		
		my $basename1 = basename($f->{R1});
		my $basename2 = basename($f->{R2});
		my $dir = dirname($f->{R1});
		system($cmd_tar." $dir $basename1 $basename2 ");
		system(" md5sum ".$f->{R1}." ".$f->{R2}." >> $md5file ");
		
	}
}

system("tar -rf $project_name.fastq.tar $md5file");
warn "write";
die();
#;
#
##foreach my $variation (@{$patient->getVariations}) {
##	warn $variation->getGnomadAC;
##	foreach my $gene (@{$variation->getGenes()}) {
##		warn $gene->external_name." ".$variation->name();
##	}
##}
##die();
#my $bam = $patient->getBamFile();
#my $ref = $project->getGenomeFasta();
#my $output = $patient->gvcfFileName("haplotypecaller_andreii");
#$output =~ s/\.gz//;
#
#warn $output;
#
#my $cmd = qq{/software/distrib/GATK/gatk-4/gatk HaplotypeCaller -R $ref  --input $bam --output $output --emit-ref-confidence GVCF --java-options '-Xmx20G -Xms20g'};
#print $cmd."\n";