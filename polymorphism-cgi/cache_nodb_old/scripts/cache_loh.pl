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
use Cache_Commons;



my $fork = 1;
my ($project_name, $chr_name);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'chr=s'        => \$chr_name,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }



warn "\n### CACHE: somatic loh model step\n";
my $nbErrors = 0;
my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => $project_name, -cache => '1', -typeFilters => 'familial' );
unless ($project->isSomaticStudy()) {
	warn "-> NOT necessary (not a somatic study)\n";
	exit(0);
}
$project->getPatients();
my $chr = $project->getChromosome($chr_name);
mkdir ($project->getCacheBitVectorDir().'/somatic_loh') unless (-d $project->getCacheBitVectorDir().'/somatic_loh');
my $freeze = $project->getCacheBitVectorDir()."/lmdb_cache/$chr_name.dv.freeze";
my $hAllVar = retrieve $freeze;
if (scalar keys %{$hAllVar} == 0) {
	my $cmd = "touch ".$project->getCacheBitVectorDir()."/somatic_loh/".$chr->id().".lite";
	`$cmd`;
	exit();
}
my $fasta_file = $project->getGenomeFasta();
my $hResults;
my $var_tmp = $chr->getNewVector();
my $var_somatic_ho  = $chr->getNewVector();
my $var_somatic_he  = $chr->getNewVector();
my $var_patient = $chr->getNewVector();
my $var_group = $chr->getNewVector();
my $var_global = $chr->getNewVector();
my $no_var_infos = $chr->get_lmdb_variations();
foreach my $group (@{$chr->getSomaticGroups()}) {
	$var_tmp->Empty();
	$var_group->Empty();
	$var_somatic_ho->Empty();
	$var_somatic_he->Empty();
	$group->used_model('loh');
	my @lPatSomatic;
	foreach my $patient (@{$group->getPatientsSomatic()}) {
		$patient->used_model('loh');
		push(@lPatSomatic, $patient);
		$var_somatic_he += $patient->getHe($chr);
		$var_somatic_ho += $patient->getHo($chr);
	}
	foreach my $patient (@{$group->getPatientsGerminal()}) {
		$patient->used_model('loh');
		$var_patient->Empty();
		$var_tmp->Intersection( $var_somatic_ho, $patient->getHe($chr) );
		$var_patient += $var_tmp;
		$var_patient -= $patient->getHo($chr);
		$var_tmp->Intersection( $var_somatic_he, $patient->getHe($chr) );
		foreach my $id (@{$chr->getListVarVectorIds($var_tmp)}) {
			my $var_id = $chr->getVarId($id);
			my $sample_var_1 = $no_var_infos->get($var_id)->{patients_details}->{$patient->name()}->{he_ho_details};
			unless ($sample_var_1) {
				warn "\n\nERROR: pb with $var_id and patient ".$patient->name()."... confess\n\n";
				confess ();
			}
			my ($t,$a,$c)  = split(":", $sample_var_1);
			my ($t1,$b,$d);
			my $already;
			foreach my $pat_somatic (@lPatSomatic) {
				next if ($already);
				next unless ($pat_somatic->getHe($chr)->contains($id));
				my $sample_var_2 = $no_var_infos->get($var_id)->{patients_details}->{$pat_somatic->name()}->{he_ho_details};
				unless ($sample_var_2) {
					warn "\n\nERROR: pb with $var_id (somatic patients: ".join(', ', @lPatSomatic).")... confess\n\n";
					confess ();
				}
				($t1,$b,$d) = split(":", $sample_var_2);
				$a = 0 unless ($a);
				$b = 0 unless ($b);
				$c = 0 unless ($c);
				$d = 0 unless ($d);
				my $p = fisher::fishers_exact( $a, $b, $c, $d, 1);
				next if $p > 0.01;
				$var_patient->Bit_On($id);
				$already = 1;
			}
		}
		$patient->getVariantsVector($chr)->Intersection( $patient->getVariantsVector($chr), $var_patient );
		$var_group  += $var_patient;
		$var_global += $var_patient;
	}
	foreach my $patient (@{$group->getPatientsSomatic()}) {
		$patient->getVariantsVector($chr)->Intersection( $patient->getVariantsVector($chr), $var_group );
		$hResults->{$patient->name()} = $patient->getVariantsVector($chr);
	}
}
$no_var_infos->close();

warn 'store 1/1: nosql loh';
my $nosql = GenBoNoSql->new(dir => $ project->getCacheBitVectorDir().'/somatic_loh', mode => 'c');
$nosql->put_bulk($chr->id(), $hResults);
$nosql->close();



