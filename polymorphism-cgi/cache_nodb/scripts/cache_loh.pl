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
my ($project_name, $chr_name,$annot_version);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'chr=s'        => \$chr_name,
	'annot_version=s'    => \$annot_version,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }



warn "\n### CACHE: somatic loh model step\n";
my $nbErrors = 0;
my $buffer = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProjectCache( -name => $project_name, -cache => '1', -typeFilters => 'familial' );
if ($annot_version) {
	$project->changeAnnotationVersion($annot_version);
}
unless ($project->isSomaticStudy()) {
	my $nosql = GenBoNoSql->new(dir => $ project->getCacheBitVectorDir().'/somatic_loh', mode => 'c');
	my $t = time;
	$nosql->put($chr_name, "time",$t);
	$nosql->close();
	
	warn "-> Create NOT necessary (not a somatic study)\n";
	exit(0);
}
$project->getPatients();
my $chr = $project->getChromosome($chr_name);
mkdir ($project->getCacheBitVectorDir().'/somatic_loh') unless (-d $project->getCacheBitVectorDir().'/somatic_loh');
my $fasta_file = $project->getGenomeFasta();
my $hResults;
my $var_tmp = $chr->getNewVector();
my $var_somatic_ho  = $chr->getNewVector();
my $var_somatic_he  = $chr->getNewVector();
my $var_patient = $chr->getNewVector();
my $var_group = $chr->getNewVector();
my $var_global = $chr->getNewVector();
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
			my $var = $chr->getVarObject($id);
			my $var_id = $var->id();
			
			my $t = 'he';
			$t = 'ho' if $var->isHomozygote($patient);
			my $a = $var->getNbAlleleRef($patient);#	unless ($nb_all_ref){$nb_all_ref = 0;}
			my $c = $var->getNbAlleleAlt($patient);#	unless ($nb_all_mut){$nb_all_mut = 0;}
			my $already;
			foreach my $pat_somatic (@lPatSomatic) {
				next if ($already);
				next unless ($pat_somatic->getHe($chr)->contains($id));
				my $t1 = 'he';
				$t1 = 'ho' if $var->isHomozygote($pat_somatic);
				my $b = $var->getNbAlleleRef($pat_somatic);#	unless ($nb_all_ref){$nb_all_ref = 0;}
				my $d = $var->getNbAlleleAlt($pat_somatic);#	unless ($nb_all_mut){$nb_all_mut = 0;}
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

my $rocks4 = $chr->rocks_vector("w");
foreach my $patient (@{$project->getPatients()}) {
	my $vloh_pat = $chr->getNewVector();
	$vloh_pat += $hResults->{$patient->name} if exists $hResults->{$patient->name};
	warn $patient->name().': '.$vloh_pat->Norm().'/'.$vloh_pat->Size();
	$rocks4->put_batch_vector_transmission($patient,"som_loh",$vloh_pat);	
}
$rocks4->close();

exit(0);




