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

my $buffer = new GBuffer;
my $project = $buffer->newProject(-name=>$project_name);
#warn("Project ".$project_name);
my $patients = $project->get_only_list_patients($patient_names);
die("No patient in project ".$project_name."\n") unless $patients;

#die("No patient to merge in project $project_name\n") unless $merge;
#my $merge;
#foreach my $p (@$patients) {
#	push(@{$merge->{$dest_project}}, $p);
#}
my $merge->{$dest_project} = $patients;


foreach my $pdes (keys %$merge){
	my $buffer2 = GBuffer->new();
	my $newProject = $buffer2->newProject(-name=>$pdes);
	warn("Project $pdes");
	
	die ("Genome release for the two projects to merge are different: $project_name ".$project->getVersion." ; $pdes ".$newProject->getVersion) if ($project->getVersion ne $newProject->getVersion);
	
	open(my $fh, '>', $project->getAlignmentDirName .'jobs_merge.txt') or confess("Can't open ".$project->getAlignmentDirName .'jobs_merge.txt');
	warn $project->getAlignmentDirName .'jobs_merge.txt';
	
	my $cmd_mv_bam;
	my $cmd_lien_fastq;
	foreach my $pat (@{$merge->{$pdes}}){
		my $pname = $pat->name;
		warn $pname;
		my $pat_dest = $newProject->getPatient($pat->name)
			or confess ("Can't find patient '$pname' to merge in the destination project $pdes ");
		
		# Vérifie que les infos des 2 patients à merger sont identiques
		# todo: sex, family, status, calling,  ?
		die ("Sex for the two patients to merge are different: $project_name ".$pat->sex." ; $pdes ".$pat_dest->sex) if ($pat->sex ne $pat_dest->sex);
		die ("Status for the two patients to merge are different: $project_name ".$pat->status." ; $pdes ".$pat_dest->status) if ($pat->status ne $pat_dest->status);
		die ("Captures for the two patients to merge are different: $project_name ".$pat->getCapture->name." ; $pdes ".$pat_dest->getCapture->name) if ($pat->getCapture->name ne $pat_dest->getCapture->name);
		die ("Sequencing machine for the two patients to merge are different: $project_name ".$pat->getRun->machine." ; $pdes ".$pat_dest->getRun->machine) if ($pat->getRun->machine ne $pat_dest->getRun->machine);
		die ("Sequencing plateform method for the two patients to merge are different: $project_name ".$pat->getRun->plateform." ; $pdes ".$pat_dest->getRun->plateform) if ($pat->getRun->plateform ne $pat_dest->getRun->plateform);
		die ("Alignment method for the two patients to merge are different: $project_name ".$pat->alignmentMethod." ; $pdes ".$pat_dest->alignmentMethod) if ($pat->alignmentMethod ne $pat_dest->alignmentMethod);
#		die ("Calling method(s) for the two patients to merge are different: $project_name ".join(', ',@{$pat->callingMethods})." ; $pdes ".join(', ',@{$pat_dest->callingMethods})) if (join(', ',@{$pat->callingMethods}) ne join(', ',@{$pat_dest->callingMethods}));
		
		# Merge les bam, puis trie le bam obtenu
		my $samtools = $buffer->software("samtools");
		my ($bam1, $bam2);
		eval {$bam1 = $pat_dest->getBamFile};
		confess ("Can't find bam file:\n$@") if ($@);
		eval {$bam2 = $pat->getBamFile} ;
		confess ("Can't find bam file:\n$@") if ($@);
		my $merged_bam = $bam1 =~ s/bam/merged\.bam/r;
		my $cmd_samtools_merge = "$samtools merge -f $merged_bam $bam1 $bam2";
		my $cmd_samtools_index = "$samtools index $merged_bam";
		warn ("'$merged_bam' already exits, overwritting") if (-e $merged_bam);
		print {$fh} "$cmd_samtools_merge -@ 20 ; $cmd_samtools_index -@ 20 ; \n";
		warn ("$cmd_samtools_merge -@ 20 ; $cmd_samtools_index -@ 20 ; \n");
		
		# Efface ou  renomme les bams partiels et mergés
		if ($rm_bam) {
			$cmd_mv_bam->{$bam1} = "rm $bam1 $bam2 $bam1.bai $bam2.bai ; mv $merged_bam $bam1 ; mv $merged_bam.bai $bam1.bai";
			warn ($cmd_mv_bam->{$bam1}."\n");
		}
		else {
			$cmd_mv_bam->{$bam1} = "mv $bam1 $bam1.sans_preseq ; mv $bam1.bai $bam1.bai.sans_preseq ;\n mv $bam2 $bam2.sans_preseq ; mv $bam2.bai $bam2.bai.sans_preseq ;\n mv $merged_bam $bam1 ; mv $merged_bam.bai $bam1.bai";
			warn ($cmd_mv_bam->{$bam1}."\n");
		}

		# Créé des liens des fastq complémentaires dans le run du projet d'origine
		my $fastq_files = $pat->fastqFiles;
		my @fastq_files = map {values %$_} @$fastq_files;
		my $dir_fastq_dest = $pat_dest->getSequencesDirectory;
		foreach my $fastq (@fastq_files) {
			$cmd_lien_fastq->{$fastq} = "ln -s $fastq $dir_fastq_dest" unless (-e $dir_fastq_dest.$fastq);
		}
		
	}
	close($fh);
	system('cat '. $project->getAlignmentDirName ."jobs_merge.txt | run_cluster.pl -cpu=20") unless ($no_exec);
	
	# Efface ou  renomme les bams partiels et mergés
	my @errors;
	foreach my $file (keys %$cmd_mv_bam){
		warn $cmd_mv_bam->{$file};
		system($cmd_mv_bam->{$file}) if (-e $file and not $no_exec);
		push(@errors, $file) unless (-e $file);
	}
	
	# Créé des liens des fastq complémentaires dans le run du projet d'origine
	foreach my $file (keys %$cmd_lien_fastq) {
		warn $cmd_lien_fastq->{$file};
		system ($cmd_lien_fastq->{$file}) unless (-e $file or $no_exec);
	}
	warn Dumper \@errors if (scalar(@errors));
	die ("Error while merging: ".join(', ',@errors)) if (scalar(@errors) and not $no_exec);
	
	# Relance les pipelines pour le calling, coverage, etc.
	my @patients_merged = map {$_->name} @{$merge->{$pdes}};
	my $cmd_dragen = "dragen_pipeline.sh -project=$pdes -patient=".join(',',@patients_merged)." -force=1";
	my $cmd_bds = "bds_pipeline.sh -project=$pdes -patient=".join(',',@patients_merged)." -force=1";
	my $cmd_calling = "bds_calling.sh -project=$pdes -patient=all -force=1";
#	my $cmd_cache = "bds_cache.sh -project=$pdes -force=1";
#	my $cmd_cache_splice = "bds_cache_rna_junctions.sh -project=$pdes -force=1";
	unless ($no_exec) {
		colored::stabilo('white', "--------DONE--------");
		colored::stabilo('yellow', "Now, run coverage and calling on the merged patients:");
		colored::stabilo('yellow', $cmd_dragen);
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


