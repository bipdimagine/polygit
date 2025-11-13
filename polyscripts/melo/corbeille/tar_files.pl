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

use GBuffer;
my $buffer = new GBuffer;


my $project_names;
my $patient_names;
my $file_types;
my $archive_dir = "/data-isilon/download/";
my $no_die;
my $archive_name;
my $help;

GetOptions(
	'project=s'		=> \$project_names,		# Projet(s)
	'patients=s'	=> \$patient_names,		# Patient(s)
	'files=s'		=> \$file_types,		# fastq, bam, vcf, htlv1 or all
	'archive_dir=s'	=> \$archive_dir,
	'archive_name=s'=> \$archive_name,
	'nodie'		=> \$no_die,			
	'help'			=> \$help,		
);

usage() unless ($project_names and $file_types);
usage() if ($help);
confess ("Directory \"$archive_dir\" does not exit") unless (-d $archive_dir);


# Récupère les fichiers
$file_types = [split(',', $file_types)];
my $files;
my @patient_names = split(',', $patient_names) if ($patient_names ne "all");

$project_names = [split(',', $project_names)];
foreach my $project_name (@$project_names) {
	my $project = $buffer->newProject(-name=>$project_name);
	colored::stabilo("white", "Project ".$project->name);
	
	my $patients = $project->get_only_list_patients($patient_names);
	warn("No patient found in project ".$project->name."\n") unless $patients;
	
	foreach my $pat (@$patients) {
		my $pat_name = $pat->name;
		@patient_names = grep {$_ ne $pat_name} @patient_names;
		colored::stabilo("yellow", $pat_name);
		
		# todo: bcl ?
		
		#fastq
		if (grep {$_ =~ /fastq/i or $_ =~ /all/i} @$file_types) {
			my $fastq_files;
			print "fastq:\n";
			eval {
				foreach my $fastq (@{$pat->fastqFiles}) {
					foreach my $r (keys %$fastq) {
						push(@$fastq_files, ($fastq->{$r})) if ($fastq->{$r} =~ /.fastq.gz$/);
					}
				}
				print "\t". join("\n\t", @$fastq_files) ."\n";
				push(@$files, @$fastq_files);
			};
			if ($@) {
				print "\t";
				colored::stabilo("red","No fastq file found");
				die unless ($no_die);
				warn $@;
			}
		}
		
		# bam
		if (grep {$_ =~ /bam/i or $_ =~ /all/i} @$file_types) {
			print "bam:\n";
			my $bam_files = $pat->getBamFiles;
			if (scalar @$bam_files) {
				print "\t". join("\n\t", @$bam_files) ."\n";
				push(@$files, @$bam_files);
				# *.bam.bai
				foreach my $bam (@$bam_files) {
					$bam .= '.bai';
				}
				print "\t". join("\n\t", @$bam_files) ."\n";
				push(@$files, @$bam_files);
			}
			else {
				colored::stabilo("red","No bam file found");
				print "\n";
				die unless ($no_die);
			};
		}
		
		# vcf
		if (grep {$_ =~ /vcf/i or $_ =~ /all/i} @$file_types) {
			print "vcf:\n";
			my $vcf_files = $pat->getVariationsFiles;
			if (scalar @$vcf_files) {
				print "\t". join("\n\t", @$vcf_files) ."\n";
				push(@$files, @$vcf_files);
				
				# *.vcf.tbi
				foreach my $vcf (@$vcf_files) {
					$vcf .= '.tbi';
				}
				print "\t". join("\n\t", @$vcf_files) ."\n";
				push(@$files, @$vcf_files);
			}
			else {
				colored::stabilo("red","No vcf file found");
				print "\n";
				die unless ($no_die);
			}
		}
		
		# htlv1
		if (grep {$_ =~ /htlv1/i} @$file_types) {
			print "htlv1:\n";
			my $htlv1_dir_path = $project->getVariationsDir("htlv1_calling");
			opendir(my $htlv1_dir, $htlv1_dir_path) or confess ("Can not open dir '$htlv1_dir_path': $!");
			my $htlv1_files; # = grep { -f $_ && /^$pat_name-/ } readdir($htlv1_dir);
			while (my $file = readdir($htlv1_dir)){
				if (-f "$htlv1_dir_path$file" & $file =~ /^$pat_name/) {
					push (@$htlv1_files, "$htlv1_dir_path$file");
				}
			}
			close($htlv1_dir);
			if (scalar @$htlv1_files) {
				print "\t". join("\n\t", @$htlv1_files) ."\n";
				push(@$files, @$htlv1_files);
			}
			if (scalar @$htlv1_files < 4) {
				colored::stabilo("magenta", 4-scalar(@$htlv1_files)." file(s) missing:");
				if (grep {$_ =~ /-clonalityResults.txt$/ } @$htlv1_files) {
					print "\t\t$htlv1_dir$pat_name-clonalityResults.txt\n";
				}
				if (grep {$_ =~ /-mergedIS.txt$/ } @$htlv1_files) {
					print "\t\t$htlv1_dir$pat_name-mergedIS.txt\n";
				}
				if (grep {$_ =~ /-mergedIS.xls$/ } @$htlv1_files) {
					print "\t\t$htlv1_dir$pat_name-mergedIS.xls\n";
				}
				if (grep {$_ =~ /-SIMPLIFIED_mergedIS.txt$/ } @$htlv1_files) {
					print "\t\t$htlv1_dir$pat_name-SIMPLIFIED_mergedIS.txt\n";
				}
			}
			elsif(not scalar @$htlv1_files) {
				colored::stabilo("red","No htlv1 file found"); # (*-clonalityResults.txt, *-mergedIS.txt, *-mergedIS.xls, *-SIMPLIFIED_mergedIS.txt)
				print "\n";
				die unless ($no_die);
			}
		}
		
		print "\n";
	}
	print "\n";
}


if (@patient_names) {
	print "Patient(s) not found in project(s) join(', ',$project_names):\n";
	foreach my $name (@patient_names) {
		colored::stabilo('red', "\"$name\"");
	}
	print "\n";
	die unless ($no_die);
}
print "Total: ". scalar @$files ." files\n";
print "\n";

die ("No file to copy\n") unless $files;



# Create the archive
$archive_name = join('-',@$project_names) .'-'. join('_',@$file_types) .'.tar' unless ($archive_name);
$archive_name =~ s/.tar.gz$/.tar/ unless (lc($file_types) eq 'htlv1');
$archive_name .= '.tar' unless ($archive_name =~ /.tar$/);
$archive_name .= '.gz' if (lc($file_types) eq 'htlv1' && /.tar.gz$/);

my $choice = prompt("Make an archive '$archive_name' of these files in '". abs_path($archive_dir) ."/' ?  (y/n) ");
die() if ($choice ne "y");

my $cmd = "tar -cvf $archive_dir/$archive_name ". join(" \\\n", @$files);
$cmd =~ s/-cvf/-cvzf/ if (lc($file_types) eq 'htlv1');
print "\n";
#warn $cmd;
system $cmd;
confess('Error while making the archive. The archive was not created') unless (-e "$archive_dir/$archive_name");
if (-e "$archive_dir/$archive_name") {
	colored::stabilo('white', "--------DONE--------");
	colored::stabilo('white', "Archive created: '$archive_dir/$archive_name'");
}
exit;



sub usage {
	print "
tar_files.pl
-------------	
Mandatory arguments
	-project <s>		project names separated with a comma
	-files <s>		file types separated with a comma (fastq, bam, vcf, htlv1 or all (all=fastq,bam,vcf))
	
Optional arguments
	-patients <s>		patient names separated with a comma
	-archive_dir <s>	path to create the archive, default '/data-isilon/download/'
	-archive_name <s>	archive name
	-nodie			do not die if file(s) or patient(s) are not found
	-help				

";
	exit(1);
}



#my $choice = prompt("Copy these files to ". abs_path($outDir) ."/    (y/n) ? ");
#die() if ($choice ne "y"); 
#
#
#my $cmd = "cp -v ". join(' ', @$files) ." $outDir";
#print "\n";
##warn $cmd;
#system $cmd;
#print "\n----------\nDone !\n----------\n";


