#!/usr/bin/env perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use GBuffer; 
use Data::Dumper;
use Getopt::Long;
use Carp;
#use Storable qw(store retrieve freeze);
use Term::ANSIColor;
#use Thread::Queue;
#use Set::IntSpan::Fast::XS;
#use String::ProgressBar;
#use List::Util qw(sum);
use autodie qw(system open);


my $project_name;
my $patient_name;
my $no_exec;
my $force;
my $threads = 20;

GetOptions(
	'project=s'				=> \$project_name,
	"patient=s"				=> \$patient_name,
	"no_exec"				=> \$no_exec,
    'force=s'				=> \$force,
	"cpu|threads|fork=i"	=> \$threads,
) || confess("\nError in command line arguments");


my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatient($patient_name);
warn $patient->name;
warn Dumper $patient->uBams_revio();

my $ubam_dir = $patient->getSequencesDirectory();
mkdir ($ubam_dir, 0775, ) unless -d $ubam_dir;
my $ubam = $patient->uBam_filename();
warn $ubam;

my $ubams_revio = $patient->uBams_revio();
my $samtools = $buffer->software("samtools");
my $ubam_size_sum = 0;
my $nb_reads_revio = 0;
foreach my $bam (@$ubams_revio){
	$ubam_size_sum += -s $bam;
	
	# Check ubam RG
	open(my $fh, "samtools view -H $bam |") or die $!;
	my $find;
	while (<$fh>) {
    	next unless /^\@RG/;
		$find = 1;
    	if (/\tSM:([^\t]+)/) {
			my $sm = $1;
			my $pname = $patient->name;
			if ($sm ne $patient->name ){
				print colored("ERROR RG/name $bam - ubam: $sm - DB: $pname", 'red'), "\n";
				print "If you want to change the RG/name in bam, use '$Bin/change_bam_rg.pl -project $project_name -patient $patient'\n";
				print "Else, change patient name in database/polyproject\n";
				die();
			}
		}
    }
	close($fh);
	die('No @RG found in header') unless $find;
	
	my $pbi = $bam.".pbi";
	warn "run_singularity.sh pbbam.sif pbindexdump $pbi"; 
	open(my $fh, "-|", "run_singularity.sh pbbam.sif pbindexdump $pbi") or confess "Impossible de lancer pbdump: $!";
	while (my $line = <$fh>) {
		if ($line =~ /"numReads"\s*:\s*(\d+)/) {
			warn $1;
			$nb_reads_revio += $1;
			last;
		}
		# sécurité : on ne lit que les 30 premières lignes
		last if $. > 30;
	}
	close($fh);
	
}
#warn "Sum ubam sizes: " . format_bits($ubam_size_sum);
warn "Nombre de reads total Revio : $nb_reads_revio \n";	

# Check reads nb if ubam already exists
if (-e $ubam) {
	unless ($force) {
#	my $ubam_size = -s $ubam;
#	warn "Merged ubam size: " . format_bits($ubam_size);
#	my $size_ratio = abs($ubam_size_sum - $ubam_size) / $ubam_size_sum;
#	warn "Merged BAM size does not match the sum of Revio BAM sizes: " . format_bits($ubam_size) . " vs " . format_bits($ubam_size_sum) . ". Force merge.";
#	$force = 1 if ($size_ratio < 0.1);
		my ($nb_reads_merged) = `samtools idxstats -@ $threads $ubam` =~ /^\*\t\d+\t\d+\t(\d+)/m;
		warn "Nombre de reads ubam mergé : $nb_reads_merged \n";	
		chomp($nb_reads_merged);
		if ($nb_reads_merged != $nb_reads_revio) {
			warn (colored("Le nombre de reads ne correspond pas: ubams Revio $nb_reads_revio vs ubam mergé $nb_reads_merged. Merge de nouveau", 'yellow'));
			$force = 1;
		}
	}
	if ($force){
		unlink $ubam;
		unlink $ubam.".bai" if -e $ubam.".bai";
	}
	else {
		warn "Already done: $ubam";
		exit();
	}
}

# Merge
my $cmd  = qq{$samtools merge --write-index -f -o $ubam -@ $threads }.join(" ",@$ubams_revio);#.qq{ && $samtools index $ubam};
#my $cmd  = qq{run_singularity.sh pbbam.sif pbmerge -o $ubam -j $threads }.join(" ",@$ubams_revio);
warn $cmd;
system($cmd) unless ((-e $ubam and not $force) or $no_exec);
 
# Check sizes
#my $ubam_size = -s $ubam;
#my $size_ratio = abs($ubam_size_sum - $ubam_size) / ($ubam_size_sum + $ubam_size);
#if ($size_ratio >= 0.1) {
#	system("rm $ubam");
#	confess("ERROR $project_name $patient_name" . "\n" . "Total ubam size: " . format_bits($ubam_size_sum) . " vs merged ubam size: " . format_bits($ubam_size));
#}

if (-e $ubam) {
	my ($nb_reads_merged) = `samtools idxstats $ubam` =~ /^\*\t\d+\t\d+\t(\d+)/m;
	warn "Nombre de reads ubam mergé : $nb_reads_merged \n";	
	chomp($nb_reads_merged);
	if ($nb_reads_merged != $nb_reads_revio) {
		#unlink $ubam;
		confess (colored("ERROR Le nombre de reads ne correspond pas: ubams Revio $nb_reads_revio vs ubam mergé $nb_reads_merged", 'red'));
	}
	else {
		warn colored("Done! Output: $ubam", 'green') . "\n";
	}
}

exit(0);


sub format_bits {
    my ($bits) = @_;
    
    return "0 B" if $bits == 0;
    
    my @units = ('', 'K', 'M', 'G', 'T', 'P', 'E');
    my $unit_index = 0;
    my $size = $bits;
    
    while ($size >= 1024 && $unit_index < $#units) {
        $size /= 1024;
        $unit_index++;
    }
    
    # Formater avec une précision adaptée
    if ($size == int($size)) {
        return sprintf("%d%s", $size, $units[$unit_index]);
    } else {
        return sprintf("%.1f%s", $size, $units[$unit_index]);
    }
}
