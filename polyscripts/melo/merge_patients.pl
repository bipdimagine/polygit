#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/";
use lib "$Bin/../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../GenBo/lib/obj-nodb/packages/";
use lib "$Bin/../../GenBo/lib/GenBoDB/writeDB/";
use lib "$Bin/packages";
use Data::Dumper;
use Getopt::Long;
use Carp;
use IO::Prompt;
use colored;
use LWP::UserAgent;
use GBuffer;


my $project_name;
my $dest_project;
my $patient_names;
my $rm_bam;
my $no_exec;
my $help;

GetOptions(
	'project=s'		=> \$project_name,
	'target=s'		=> \$dest_project,
	'patients=s'	=> \$patient_names,
	'rm_bam'		=> \$rm_bam,
	'no_exec'		=> \$no_exec,
	'help'			=> \$help,
);

usage() unless ($project_name and $dest_project);
usage() if ($help);

my @patient_names = split(',', $patient_names) if ($patient_names ne "all");

my $buffer = new GBuffer;
my $project = $buffer->newProject(-name=>$project_name);
#warn("Project ".$project_name);
my $patients = $project->get_only_list_patients(@patient_names);
#my $patients = $project->getPatients;
die("No patient in project ".$project_name."\n") unless $patients;

#my $merge;
#foreach my $pat (@$patients){
#	my $dest_project = $buffer->getQuery->getProjectDestination($pat->id);
#	push(@{$merge->{$dest_project}},$pat) if ($dest_project && $dest_project ne $project_name)
#}
#die("No patient to merge in project $project_name\n") unless $merge;
my $merge->{$dest_project} = $patients;

foreach my $pdes (keys %$merge){
	my $buffer2 = GBuffer->new();
	my $newProject = $buffer2->newProject(-name=>$pdes);
	warn("Project $pdes");
	
	open(my $fh, '>', $project->getAlignmentDirName .'jobs_merge.txt') or confess("Can't open ".$project->getAlignmentDirName .'jobs_merge.txt');
	warn $project->getAlignmentDirName .'jobs_merge.txt';
	
	my $cmd_mv_bam;
	foreach my $pat (@{$merge->{$pdes}}){
		my $pname = $pat->name;
		my $pat_dest = $newProject->getPatient($pat->name)
			or confess ("Can't find patient '$pname' to merge in the destination project $pdes ");
		
		# Vérifie que les infos des 2 patients à merger sont identiques (méthodes alignement et calling)
		die ("Alignment methods for the two patients to merge are not the same: $project_name".$pat->alignmentMethod." ; $pdes ".$pat_dest->alignmentMethod) if ($pat->alignmentMethod ne $pat_dest->alignmentMethod);
		die ("Calling methods for the two patients to merge are not the same: $project_name ".join(', ',@{$pat->callingMethods})." ; $pdes ".join(', ',@{$pat_dest->callingMethods})) if (join(', ',@{$pat->callingMethods}) ne join(', ',@{$pat_dest->callingMethods}));
		die ("Captures for the two patients to merge are not the same: $project_name".$pat->getCapture->name." ; $pdes ".$pat_dest->getCapture->name) if ($pat->getCapture->name ne $pat_dest->getCapture->name);
		
		# Merge les bam, puis trie le bam obtenu
		my $samtools = $buffer->software("samtools");
		my ($bam1, $bam2);
		eval {$bam1 = $pat_dest->getBamFile};
		eval {$bam2 = $pat->getBamFile} ;
		confess ("Can't find bam file:\n$@") if ($@);
		my $merged_bam = $bam1 =~ s/bam/merged\.bam/r;
		my $cmd_samtools_merge = "$samtools merge -f $merged_bam $bam1 $bam2";
		my $cmd_samtools_index = "$samtools index $merged_bam";
		warn ("'$merged_bam' already exits, overwritting") if (-e $merged_bam);
#		warn $cmd_samtools_merge;
#		warn $cmd_samtools_index;
#		system("echo \"$cmd_samtools_merge -@ 20 ; $cmd_samtools_index -@ 20 \" | run_cluster.pl -cpu=20");
		print {$fh} "$cmd_samtools_merge -@ 20 ; $cmd_samtools_index -@ 20 ; \n";
		
		# Efface les bams partiels ou les renomme
		if ($rm_bam) {
			$cmd_mv_bam->{$merged_bam} = "rm $bam1 $bam2 $bam1.bai $bam2.bai ; mv $merged_bam $bam1 ; mv $merged_bam.bai $bam1.bai";
		}
		else {
#			my $bam1part = $bam1 =~ s/bam/part\.bam/r;
#			my $bam2part = $bam2 =~ s/bam/part\.bam/r;
			$cmd_mv_bam->{$merged_bam} = "mv $bam1 $bam1.part ; mv $bam1.bai $bam1.bai.part ; mv $bam2 $bam2.part ; mv $bam2.bai $bam2.bai.part ; mv $merged_bam $bam1 ; mv $merged_bam.bai $bam1.bai";
		}
#		warn $cmd_mv_bam;
#		system($cmd_mv_bam);
#		print {$fh} "$cmd_mv_bam\n";

		#todo: déplacer fastq ?
	}
	close($fh);
	system('cat '. $project->getAlignmentDirName ."jobs_merge.txt | run_cluster.pl -cpu=20") unless ($no_exec);
	my @errors;
	foreach my $file (keys %$cmd_mv_bam){
		system($cmd_mv_bam->{$file}) if (-e $file and not $no_exec);
		push(@errors, $file=~s/\.merged\.bam//r) unless (-e $file);
	}
	die ("Error while merging: ".join(', ',@errors)) if (scalar(@errors) and not $no_exec);
	
	# Relance les pipelines pour le calling, coverage, etc.
	my @patients_merged = map {$_->name} @{$merge->{$pdes}};
#	my $cmd_dragen = "dragen_pipeline.sh -project=$pdes -patient=".join(',',@patients_merged)." -force=1";
	my $cmd_bds = "bds_pipeline.sh -project=$pdes -patient=".join(',',@patients_merged)." -steps=coverage,binary_depth -force=1";
	my $cmd_calling = "bds_calling.sh -project=$pdes -patient=all -force=1";
#	my $cmd_cache = "bds_cache.sh -project=$pdes -force=1";
#	my $cmd_cache_splice = "bds_cache_rna_junctions.sh -project=$pdes -force=1";
	unless ($no_exec) {
		colored::stabilo('white', "--------DONE--------");
		colored::stabilo('yellow', "Now, run coverage and calling on the merged patients:");
		colored::stabilo('yellow', $cmd_bds);
		colored::stabilo('yellow', $cmd_calling);
	}
}





sub usage {
	print "
merge_patients.pl
-----------------	
Mandatory arguments
	-project <s>		project name containing the complementary bam
	-target <s>		destination project to merge with
	
Optional arguments
	-patients <s>		patient names separated with a comma
	-rm_bam			delete partial bam file that have been merged, default false
	-no_exec		do not run the commands
	-help				

";
#	-no_pipeline		do not run the pipelines for coverage, calling and cache after the merge
	exit(1);
}


