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
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use List::Util qw(shuffle);
require "$RealBin/Cache_Commons.pm";
use Sys::Hostname;

 my $host = hostname();


#warn "*_*_*_*_*_".$host."*_*_*_*_*_";
#use Cache_Commons;



my $fork = 1;
my ($project_name, $patient_name, $fileout);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'patient=s'    => \$patient_name,
	'fileout=s'    => \$fileout,
);

open(FILE, ">$fileout");
my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => $project_name );
$project->getChromosomes();
my $patient = $project->getPatient($patient_name);
#$patient->use_not_filtred_junction_files(1);


if (not $project->is_human_genome()) {
	print FILE 'PASS not Human release';
	close (FILE);
	exit(0); 
}

my $hType_patients;
$hType_patients = $project->get_hash_patients_description_rna_seq_junction_analyse() if (-d $project->get_path_rna_seq_junctions_analyse_description_root());

my $h_chr_vectors;
my $j_total = 0;
foreach my $chr (@{$project->getChromosomes()}) {
	my $vector_patient = $patient->getJunctionsVector($chr);
	$j_total += $chr->countThisVariants($vector_patient);
	$h_chr_vectors->{$chr->id()} = $vector_patient->Clone();
}

if ($j_total >= 10000) {
	my $no_cache = $patient->get_lmdb_cache("r");
	my $cache_vectors_enum_id = 'splices_linked_'.$patient->name().'_'.$j_total.'_chr_vectors_enum';
	my $h_res_v_enum = $no_cache->get_cache($cache_vectors_enum_id);
	if ($h_res_v_enum) {
		foreach my $chr_id (keys %$h_res_v_enum) {
			my $v_filters = $project->getChromosome($chr_id)->getNewVector();
			$v_filters->from_Enum($h_res_v_enum->{$chr_id}->{min8});
			$h_chr_vectors->{$chr_id}->Intersection($h_chr_vectors->{$chr_id}, $v_filters);
		}
	}
	else { die; }
}

my $nb = 0;
my @lJunctions;
foreach my $chr (@{$project->getChromosomes()}) {
	my @lScores = (5..10);
	@lScores = reverse sort {$a <=> $b} @lScores;
	foreach my $score (@lScores) {
		last if ($nb == 1500);
		foreach my $junction (@{$chr->getListVarObjects($h_chr_vectors->{$chr->id()})}) {
			next if ($junction->isCanonique($patient));
			next if ($junction->junction_score_without_dejavu_global($patient) < $score);
			push(@lJunctions, $junction);
			$nb++;
			last if ($nb == 1500);
		}
	}
}

if (scalar(@lJunctions) == 0) {
	print FILE '0 junction';
	close (FILE);
	exit(0); 
}

if (not $hType_patients or ($hType_patients and exists $hType_patients->{$patient->name()}->{pat})) {
	my $pm = new Parallel::ForkManager($fork);
	foreach my $junction (@lJunctions) {
		my $pid = $pm->start and next;
		$junction->createListSashimiPlots($patient);
		print FILE 'Ok junction '.$junction->id()."\n";
		$pm->finish();
	}
	$pm->wait_all_children();
}
else {
	print FILE "Patient $patient_name is a control\n";
}

close (FILE);
exit(0);
