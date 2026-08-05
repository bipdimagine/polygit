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
use Term::ANSIColor;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use String::ProgressBar;
use List::Util qw(sum);
use autodie qw(system open);



 my $project_name;
 my $fork;
 my $callable_intspan_file;
 my $patient_name;
 #my $low_calling;
 my $method;
 my $force;
 
 
my $fork = 5;
GetOptions(
	'project=s'   => \$project_name,
	"patient=s" => \$patient_name,
	"fork=s" => \$fork,
    'force'		=> \$force,
);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatient($patient_name);

warn Dumper $patient->uBams_revio();

my $ubam_dir = $patient->getSequencesDirectory();
mkdir ($ubam_dir, 0775, ) unless -d $ubam_dir;
my $ubam = $patient->uBam_filename();
warn $ubam;

my $samtools = $buffer->software("samtools");
my $find = "SM:".$patient->name;
foreach my $bam (@{$patient->uBams_revio()}){
	
	# Check ubam RG
	open(my $fh, "samtools view -H $bam |") or die $!;
my $find;
while (<$fh>) {

    next unless /^\@RG/;
	$find =1;
    if (/\tSM:([^\t]+)/) {
    	my $sm = $1;
    	my $pname = $patient->name;
		if ($sm ne $patient->name ){
		print colored("$bam - $sm - $pname", 'red'), "\n";
		die();
		}

    }

}

close($fh);
die("No @RG found in header") unless $find;
	
}
warn $ubam;
my $size1 ;
my $nreads = 0;
foreach my $file (@{$patient->uBams_revio()}){
	die($project->name) unless -e $file;
	my $v =  -s $file;
	warn $v." ".$file;
	 $size1 += $v;

	my $pbi = $file.".pbi";
	warn "run_singularity.sh pbbam.sif pbindexdump $pbi"; 
	open(my $fh, "-|", "run_singularity.sh pbbam.sif pbindexdump $pbi") or die "Impossible de lancer pbdump: $!";
	while (my $line = <$fh>) {
		if ($line =~ /"numReads"\s*:\s*(\d+)/) {
    	    $nreads += $1;
			last;
		}
		# sécurité : on ne lit que les 20 premières lignes
		last if $. > 30;
	}
	close($fh);
	print "Nombre de reads : $nreads : \n";	
}

# Check reads nb if ubam already exists
if (-e $ubam ){
	unless ($force) {
		my ($nb_reads_merged) = `samtools idxstats -@ $threads $ubam` =~ /^\*\t\d+\t\d+\t(\d+)/m;
		warn "Nombre de reads ubam mergé : $nb_reads_merged \n";	
		chomp($nb_reads_merged);
		if ($nb_reads_merged != $nb_reads_revio) {
			warn ("Le nombre de reads ne correspond pas: ubams Revio $nb_reads_revio vs ubam mergé $nb_reads_merged. Merge de nouveau");
			$force = 1;
		}
	}
	if ($force){
		unlink $ubam;
		unlink $ubam.".bai" if -e $ubam.".bai";
	}
	else {
		warn "already done";
		exit();
	}
}

# Merge
my $cmd  = qq{$samtools merge -f -o $ubam -@ $fork --write-index  }.join(" ",@{$patient->uBams_revio()});

system($cmd) ;#unless -e $ubam;
my ($v) = `samtools idxstats $ubam` =~ /^\*\t\d+\t\d+\t(\d+)/m;
chomp($v);

unless ($nreads == $v){
	warn $nreads." ".$v;
	unlink $ubam;
	die();
}

 




