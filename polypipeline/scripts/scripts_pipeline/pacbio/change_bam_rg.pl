#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use warnings;
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use GBuffer; 
use Data::Dumper;
use Getopt::Long;
use Carp;
use Term::ANSIColor;
use autodie qw(system);

my $project_name;
my $patient_name;
my $no_exec;
my $force;
my $threads = 20;

GetOptions(
	'project=s'		=> \$project_name,
	'patient=s'		=> \$patient_name,
	'no_exec'		=> \$no_exec,
	'force'			=> \$force,
	"cpu|threads|fork=i"	=> \$threads,
) || confess("\nError in command line arguments");

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatient($patient_name);

my $ubam_dir = $patient->getSequencesDirectory();
mkdir($ubam_dir, 0775) unless -d $ubam_dir;
my $ubam = $patient->uBam_filename();

# Répertoire temporaire pour les fichiers intermédiaires
my $dir_tmp = $project->getAlignmentPipelineDir('change_bam_rg');

my $samtools = $buffer->software("samtools");
my $pname = $patient->name;

warn "Patient: $pname" . "\n";
warn "Output BAM: $ubam" . "\n";
warn "Temp directory: $dir_tmp" . "\n";
warn "Input BAMs: " . join(", ", @{$patient->uBams_revio()}) . "\n\n";


# Vérifications préliminaires
warn "=== Preliminary Check ===" . "\n\n";
# Vérifier le nombre de reads total des ubams Revio
my $nb_reads_revio = 0;
foreach my $bam (@{$patient->uBams_revio()}){
	my $pbi = $bam.".pbi";
	#warn "run_singularity.sh pbbam.sif pbindexdump $pbi"; 
	open(my $fh, "-|", "run_singularity.sh pbbam.sif pbindexdump $pbi") or confess "Impossible de lancer pbdump: $!";
	while (my $line = <$fh>) {
		if ($line =~ /"numReads"\s*:\s*(\d+)/) {
			$nb_reads_revio += $1;
			last;
		}
		# sécurité : on ne lit que les 30 premières lignes
		last if $. > 30;
	}
	close($fh);
}
#warn "Total number of reads from Revio ubams : $nb_reads_revio \n";	


# Le BAM de sortie existe-t-il déjà avec le bon SM ?
if (-e $ubam && !$force) {
	warn colored("Output BAM already exists: $ubam", 'yellow');
	
	# Vérifier le SM du BAM existant
	my $existing_sm;
	open(my $fh_check, "$samtools view -H $ubam |") or die "Cannot read header of $ubam: $!";
	while (<$fh_check>) {
		if (/^\@RG/) {
			if (/\tSM:([^\t]+)/) {
				$existing_sm = $1;
				last;
			}
		}
	}
	close($fh_check);
	
	# Vérifier le nombre de reads dans le BAM existant
	my ($nb_reads_merged) = `samtools idxstats -@ $threads $ubam` =~ /^\*\t\d+\t\d+\t(\d+)/m;
#	warn "Nombre de reads ubam mergé : $nb_reads_merged \n";	
	
	if (defined $existing_sm && $existing_sm eq $pname) {
		warn colored("  Existing SM: '$existing_sm' (matches patient name)", 'green'), "\n";
		if ($nb_reads_merged == $nb_reads_revio) {
			warn colored("\nOutput BAM already exists with correct SM and size. ", 'green'), "Use --force to overwrite.\n\n";
			exit(0);
		}
		warn colored("  Number of reads not matching: sum BAMs Revio $nb_reads_revio vs merged BAM $nb_reads_merged", 'red'), "\n";
	} else {
		warn colored("  Existing SM: '$existing_sm' (does NOT match patient name '$pname')", 'red'), "\n";
	}
	warn colored("  Will proceed with regeneration.\n\n", 'yellow');
} elsif (-e $ubam && $force) {
	warn colored("Output BAM already exists: $ubam", 'yellow'), "\n";
	warn colored("--force specified: will overwrite existing BAM.\n\n", 'yellow');
}

# Phase 1: Analyser tous les BAMs et collecter les mismatches
my @mismatches;
my @final_bams;
my @temp_files_to_clean;

warn "=== Phase 1: Analysis ===" . "\n\n";

foreach my $bam (@{$patient->uBams_revio()}) {
	warn "Processing: $bam" . "\n";
	
	# Extraire l'en-tête
	my $bam_basename = (split(/\//, $bam))[-1];
	my $header_file = "$dir_tmp/${bam_basename}.header.tmp";
	system("$samtools view -H $bam > $header_file");
	push @temp_files_to_clean, $header_file;
	
	my $current_sm;
	open(my $fh, '<', $header_file) or die "Cannot open $header_file: $!";
	while (<$fh>) {
		if (/^\@RG/) {
			if (/\tSM:([^\t]+)/) {
				$current_sm = $1;
				last;
			}
		}
	}
	close($fh);
	
	if (defined $current_sm && $current_sm ne $pname) {
		warn colored("  SM mismatch: '$current_sm' -> '$pname'", 'red');
		push @mismatches, {
			bam => $bam,
			old_sm => $current_sm,
			new_sm => $pname,
		};
	} else {
		warn colored("  SM OK: '$current_sm'", 'green');
	}
	warn "\n";
}

# Phase 2: Demander confirmation si des mismatches ont été trouvés
if (@mismatches) {
	warn "=== Phase 2: Confirmation Required ===" . "\n\n";
	warn colored("Found " . scalar(@mismatches) . " BAM file(s) with SM mismatch:", 'red') . "\n";
	warn "-" x 80, "\n";
	
	foreach my $m (@mismatches) {
		warn sprintf("  %-50s  SM: %s -> %s\n", 
			$m->{bam}, 
			colored($m->{old_sm}, 'red'), 
			colored($m->{new_sm}, 'green')
		);
	}
	warn "-" x 80, "\n\n";
	
	unless ($force) {
		print colored("Do you want to proceed with these changes? [y/N]: ", 'yellow');
		my $response = <STDIN>;
		chomp($response);
		
		if (lc($response) ne 'y' && lc($response) ne 'yes') {
			warn "Aborted by user.";
			# Nettoyer les fichiers d'en-tête temporaires
			foreach my $tmp_file (@temp_files_to_clean) {
				unlink($tmp_file) if -e $tmp_file;
			}
			exit(1);
		}
	}
	warn "Proceeding with changes...\n\n";
} else {
	warn "No SM mismatches found. All BAMs are correct.\n\n";
}

# Phase 3: Appliquer les corrections si nécessaire
warn "=== Phase 3: Processing ===" . "\n\n";

foreach my $bam (@{$patient->uBams_revio()}) {
	my $bam_basename = (split(/\//, $bam))[-1];
	my $header_file = "$dir_tmp/${bam_basename}.header.tmp";
	my $new_header_file = "$dir_tmp/${bam_basename}.header.new";
	my $processed_bam = "$dir_tmp/${bam_basename}.rg_fixed.bam";
	
	# Vérifier si ce BAM nécessite une modification
	my $mismatch = (grep { $_->{bam} eq $bam } @mismatches)[0];
	
	if ($mismatch) {
		warn "Fixing: $bam";
		
		# Relire l'en-tête et modifier le SM
		my @new_header_lines;
		open(my $fh, '<', $header_file) or die "Cannot open $header_file: $!";
		while (<$fh>) {
			if (/^\@RG/) {
				s/\tSM:[^\t]+/\tSM:$pname/;
			}
			push @new_header_lines, $_;
		}
		close($fh);
		
		# Écrire le nouvel en-tête
		open(my $fh_out, '>', $new_header_file) or die "Cannot write $new_header_file: $!";
		print $fh_out @new_header_lines;
		close($fh_out);
		push @temp_files_to_clean, $new_header_file;
		
		# Créer un BAM temporaire avec le nouvel en-tête
		my $cmd = "$samtools reheader $new_header_file $bam > $processed_bam";
		warn "  Reheader: $cmd";
		system($cmd) unless $no_exec;
		
		# Indexer le BAM corrigé
#		system("$samtools index $processed_bam") unless $no_exec;
		
		push @final_bams, $processed_bam;
		push @temp_files_to_clean, $processed_bam, "$processed_bam.bai";
		
		warn "  Done!\n";
	} else {
		warn "Skipping (SM OK): $bam";
		push @final_bams, $bam;
		warn "\n";
	}
}

# Phase 4: Merge des BAMs
warn "=== Phase 4: Merging ===" . "\n\n";
warn "Merging " . scalar(@final_bams) . " BAM files...";
my $cmd = "$samtools merge --write-index -f -o $ubam -@ $threads " . join(" ", @final_bams);
warn "$cmd";
system($cmd) unless (-e $ubam or $no_exec);

# Indexer le BAM final
#warn "Indexing final BAM...";
#system("$samtools index $ubam") unless (-e "$ubam.bai" or $no_exec);

# Nettoyer tous les fichiers temporaires
warn "\n=== Phase 5: Cleanup ===" . "\n\n";
warn "Cleaning temporary files...";
foreach my $tmp_file (@temp_files_to_clean) {
	if (-e $tmp_file) {
#		warn "  Removing: $tmp_file";
		unlink($tmp_file) unless $no_exec;
	}
}

# Phase 4: Vérification du nombre de reads
warn "=== Phase 4: Final check ===" . "\n\n";
unless ($no_exec and ! -e $ubam) {
	my ($nb_reads_merged) = `samtools idxstats $ubam` =~ /^\*\t\d+\t\d+\t(\d+)/m;
	warn "Nombre de reads ubam mergé : $nb_reads_merged \n";	
	chomp($nb_reads_merged);
	if ($nb_reads_merged != $nb_reads_revio) {
		#unlink $ubam;
		confess (colored("ERROR Le nombre de reads ne correspond pas: ubams Revio $nb_reads_revio vs ubam mergé $nb_reads_merged", 'red'));
	}
}

warn colored("Done! Output: $ubam", 'green') . "\n" unless $no_exec;
warn colored("Dry run: no merged bam file produced", 'yellow'), "\n" if $no_exec;



