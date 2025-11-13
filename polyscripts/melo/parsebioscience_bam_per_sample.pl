#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/";
use lib "$Bin/../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../GenBo/lib/obj-nodb/packages/";
use lib "$Bin/../../GenBo/lib/GenBoDB/writeDB/";
use lib "$Bin/packages";use Data::Dumper;
use Getopt::Long;
use Carp;
use IO::Prompt;
use colored;
use Cwd 'abs_path';
use Net::SFTP;
use Archive::Tar;

#use Spreadsheet::Reader::ExcelXML;

use GBuffer;
my $buffer = new GBuffer;


my $project_name;	# 'NGS2024_7793' -> TWIST1_EBV_Pool1 ; 'NGS2024_8252' -> TWIST2_EBV_Pool2
my $patient_names;
my $force;
my $no_exec;
my $help;

GetOptions(
	'project=s'		=> \$project_name,		# Projet
	'patients=s'	=> \$patient_names,		# Patient(s)/sample(s)
	'force'			=> \$force,
	'no_exec'		=> \$no_exec,
	'help'			=> \$help,		
);

usage() unless ($project_name);
usage() if ($help);

my $project = $buffer->newProject(-name=>$project_name);

# Récupère le nom des libraires
my $patients = $project->getPatients;
#my @pat_names = map {$_->name} @$patients;
#warn Dumper \@pat_names;
warn("No patient/library in project ".$project->name."\n") unless $patients;
my @libraries = grep {$_->name =~ /^TWIST2_EBV_Pool2_/} @$patients;
@libraries = map {$_->name} @libraries;
#@libraries = map {$_->name =~ s/^TWIST2_EBV_Pool2_BC/TWIST2_EBV_Pool2_unmapped_BC/r} @libraries;
print(scalar(@libraries)." libraries:\n\t". join("\n\t",@libraries). "\n\n");

# Récupère le nom des échantillons
my @samples;
@samples = split /,/, $patient_names if ($patient_names);
@samples = ('AG_HD1_NS','AG_HD1_IFNa','AG_HD2_NS','AG_HD2_IFNa','AG_HD3_NS','AG_HD3_IFNa','LH','VO','Gauth','Lantaz') unless ($patient_names);
#@samples = ('Lantaz');
die ("No patient/sample") unless (scalar @samples);
print(scalar(@samples)." samples:\n\t". join("\n\t",@samples). "\n\n");

# todo: read sample names from excel loading table
#my $sample_sheet = $project->getRun->sample_sheet;
#my $sample_sheet = $project->getRun->infosRun->{document};

my $dir = $project->getProjectRootPath;
#$dir =~ s/NGS2024_8252/NGS2024_8252.genome/;
warn $dir;
my $samtools = $buffer->software("samtools");
print ("\n");

# Erreur s'il n'existe pas $lib/all-sample/barcode_headAligned_anno.bam
my $errors;
foreach my $lib (@libraries) {
	my $bam_lib = "$dir$lib/process/barcode_headAligned_anno.bam";
	$errors->{$lib} = "Error: can't find $bam_lib" unless (-e $bam_lib);
}
confess (join("\n", values %$errors)) if ($errors);

# Extract sample bam from each library
my $jobs;
my $bam_patient_lib;
my @bam_merged;
my $cpu = 2;
foreach my $pat (@samples) {
	foreach my $lib (@libraries) {
		my $bam_extracted = "$dir$lib/$pat/$pat\_$lib\_by_bc.bam";
		$bam_patient_lib->{$pat}->{$lib} = $bam_extracted;
		my $cmd_bc = "grep -v bc_wells $dir$lib/$pat/DGE_unfiltered/cell_metadata.csv | cut -d ',' -f 1 > $dir$lib/$pat/$pat\_bc.csv";
#		my $cmd_filter = "$samtools view -@ $cpu -D CB:$dir$lib/$pat/$pat\_bc.csv $dir$lib//process/barcode_headAligned_anno.bam -h -b -o $bam_extracted";
		my $cmd_filter = qq{$samtools view -@ $cpu $dir$lib//process/barcode_headAligned_anno.bam -h | grep -E '^@|"\$(sed 's/^/^/' $dir$lib/$pat/$pat\_bc.csv | tr '\\n' '|' | sed 's/\|\$//')"' | $samtools view -@ $cpu -b -o $bam_extracted};
		$cmd_filter =~ s/barcode_headAligned_anno\.bam/barcode_headAligned\.out\.bam/ if (-e "$dir$lib/process/barcode_headAligned.out.bam");
#		my $cmd_filter = "$samtools view -@ $cpu -d pS:$pat $dir$lib/process/barcode_headAligned_anno.bam -h -b -o $bam_extracted";
		$jobs->{'extract'} .= "$cmd_bc ; " if ($cmd_bc and (not -e $bam_extracted or $force));
		$jobs->{'extract'} .= "$cmd_filter\n" unless (-e $bam_extracted and not $force);
	}
}

chomp $jobs->{'extract'};
print ("Extract sample bam from each library:\n".$jobs->{'extract'}."\n\n") if ($jobs->{'extract'});
system("echo \"".join("\n",$jobs->{'extract'})."\" | run_cluster.pl -cpu=$cpu") unless ($no_exec or not $jobs->{'extract'});
sleep(3) unless ($no_exec);

# Merge sample bam from each library
$cpu = 5;
foreach my $pat (@samples) {
	foreach my $lib (@libraries) {
		my $bam_extracted = "$dir$lib/$pat/$pat\_$lib\_by_bc.bam";
		$errors->{$pat}->{$lib} = "Error while extracting the reads of sample '$pat' from library '$lib': missing '$bam_extracted'." unless (-e $bam_extracted or $no_exec);
	}
	unless ($errors->{$pat}) {
		my $bam_merged = "$dir$pat/$pat\_by_bc.bam";
		my $cmd_mkdir = "mkdir $dir$pat/ ;" unless (-d "$dir$pat/");
		my $cmd_merge = "$samtools merge -@ $cpu -f $bam_merged ".join(' ', values %{$bam_patient_lib->{$pat}});# $dir*/$pat/$pat*.bam";
		my $cmd_sort = "$samtools sort -@ $cpu $bam_merged -o $dir$pat/$pat.sorted.bam";
		my $cmd_mv = "mv $dir$pat/$pat.sorted.bam $bam_merged";
		my $cmd_index = "$samtools index -@ $cpu $bam_merged";
		$jobs->{'merge'} .= "$cmd_mkdir $cmd_merge ; $cmd_sort ; $cmd_mv ; $cmd_index \n" unless (-e $bam_merged and not $force);
	}
}

chomp $jobs->{'merge'};
print ("Merge sample bam from each library:\n".$jobs->{'merge'}."\n\n") if ($jobs->{'merge'});
system("echo \"".join("\n",$jobs->{'merge'})."\" | run_cluster.pl -cpu=$cpu") unless ($no_exec or not $jobs->{'merge'});

#exit if ($no_exec);

# Vérifie que les fichiers de sortie existent
my @ok;
foreach my $pat (@samples) {
	my $bam_merged = "$dir$pat/$pat\_by_bc.bam";
	$errors->{$pat} = "Error while merging bam for sample '$pat': missing '$bam_merged'" unless (-e $bam_merged or $no_exec);
	push(@ok,$pat) if (-e $bam_merged and keys %$jobs);
}
# Message pour les bam bien ou non mergés
colored::stabilo('green', "Sample reads extracted from libraries and merged") if (@ok);
foreach my $pat_ok (@ok) {
	colored::stabilo('green', "\t$pat_ok");
	print "\t$dir$pat_ok/$pat_ok\.bam\n";
}
colored::stabilo('green', "Sample bams already extracted. Maybe you wanted to use -force ?") unless (keys %$jobs or $force);

# Erreur s'il n'existe pas les bam mergés
if (keys %$errors) {
	colored::stabilo('red',"Errors");
	foreach my $pat_error (keys %$errors) {
		colored::stabilo('red', "\t$pat_error");
	}
	confess(join("\n", values %$errors));
}

sub usage {
	print "
parsebioscience_bam_per_sample.pl
-----------------	
Mandatory arguments
	-project <s>		project name containing the complementary bam
	
Optional arguments
	-patients <s>		samples names separated with a comma (not used for libraries)
	-force			overwrite existing files
	-no_exec		do not run the commands
	-help				

";
	exit(1);
}