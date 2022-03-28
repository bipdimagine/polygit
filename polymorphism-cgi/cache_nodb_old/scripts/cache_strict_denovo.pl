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



warn "\n### CACHE: strict-denovo model step\n";
my $nbErrors = 0;
my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => $project_name, -cache => '1', -typeFilters => 'familial' );
$project->getPatients();
my $chr = $project->getChromosome($chr_name);
mkdir ($project->getCacheBitVectorDir().'/strict-denovo') unless (-d $project->getCacheBitVectorDir().'/strict-denovo');
my $freeze = $project->getCacheBitVectorDir()."/lmdb_cache/$chr_name.dv.freeze";
unless (-e $freeze) {
	die ("\n\nERROR: $freeze doesn't exist... Die\n\n");
}
my $hAllVar = retrieve $freeze;
if (scalar keys %{$hAllVar} == 0) {
	my $cmd = "touch ".$project->getCacheBitVectorDir()."/strict-denovo/".$chr->id().".lite";
	`$cmd`;
	exit();
}
my $fasta_file = $project->getGenomeFasta();
my $hResults;
my $hErrors;

my $nbErrors = 0;
my $pm = new Parallel::ForkManager($fork);
$pm->run_on_finish (
	sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hResGlobal) = @_;
		if (defined($hResGlobal)) {
			foreach my $fam_name (keys %{$hResGlobal}) {
				delete $hErrors->{$fam_name};
				foreach my $pat_name (keys %{$hResGlobal->{$fam_name}}) {
					$hResults->{$pat_name} = $hResGlobal->{$fam_name}->{$pat_name};
				}
			}
		}
		else {
			$nbErrors++;
			print qq|No message received from child process $pid!\n|;
		}
	}
);
foreach my $family (@{$project->getFamilies()}) {
	$hErrors->{$family->name()} = undef;
}

foreach my $family (@{$project->getFamilies()}) {
	#next unless $family->name eq '19';
	$pm->start() and next;
	$buffer->dbh_reconnect();
	my $hTmp;
	$hTmp->{$family->name()} = undef;
	foreach my $patient (@{$family->getPatients()}) {
		$hTmp->{$family->name()}->{$patient->name()} = $chr->getNewVector();
	}
	if (scalar (@{$family->getChildrenIll()}) == 0 or scalar (@{$family->getHealthy()}) == 0) {
		$pm->finish(0, $hTmp);
		next;
	}
	else {
		my $nb_init = $chr->countThisVariants( $family->getVariantsVector($chr) );
		my $vector_denovo = $family->getModelVector_fam_denovo($chr)->Clone();
		$hResults->{$family->name()}->{global}->{denovo} = $vector_denovo->Clone();
		my $nb_denovo = $chr->countThisVariants( $hResults->{$family->name()}->{global}->{denovo} );
		my @lVarIds = @{$chr->getListVarVectorIds($vector_denovo)};
		foreach my $patient (@{$family->getHealthy()}) {
			my $hVarDeleted;
			my @lBamFiles = @{$patient->getBamFiles()};
			foreach my $bam_file (@lBamFiles) {
				#TODO: mettre HTS plus tard
				my $sam = Bio::DB::Sam->new(-bam=>$bam_file, -fasta=>$project->getGenomeFasta());
				foreach my $vector_id (@lVarIds) {
					my $var = $chr->getVarObject($vector_id);
					next if (exists $hVarDeleted->{$vector_id});
					my $to_keep;
					my $limit = 2;
					if ($var->isLarge()) {
						$hVarDeleted->{$vector_id} = undef;
						next;
					}
					elsif ($var->isVariation) {
						$to_keep = check_substitution($sam, $var, $limit);
					}
					elsif ($var->isDeletion) {
						my $chr = 'chr'.$var->getChromosome->id();
						my $start = $var->start() - 2;
						my $end = $var->end() + 2;
						my $locus = $chr.':'.$start.'-'.$end;
						$to_keep = check_del($bam_file, $var->type_public_db(), $locus, $limit);
					}
					elsif ($var->isInsertion) {
						my $chr = 'chr'.$var->getChromosome->id();
						my $length = 3;
						if ($length > 3) { $length = $var->length(); }
						my $start = $var->start() - $length - 1;
						my $end = $var->start() + $length;
						my $locus = $chr.':'.$start.'-'.$end;
						$limit = 1;
						$to_keep = check_ins($bam_file, $var->type_public_db(), $locus, $limit);
					}
					unless ($to_keep) { $hVarDeleted->{$vector_id} = undef; }
				}
			}
			my $vector_to_del = Bit::Vector->new_Enum($vector_denovo->Size(), join(',', keys %{$hVarDeleted}));
			$vector_denovo -= $vector_to_del;
		}
		my $vector_strict_denovo = $vector_denovo->Clone();
		foreach my $patient (@{$family->getPatients()}) {
			$hTmp->{$family->name()}->{$patient->name()} = $patient->getVariantsVector($chr)->Clone();
			$hTmp->{$family->name()}->{$patient->name()}->Intersection($hTmp->{$family->name()}->{$patient->name()}, $vector_strict_denovo);
		}
		$pm->finish(0, $hTmp);
	}
}
$pm->wait_all_children();
$buffer->dbh_reconnect();

if ($nbErrors > 0) {
	warn ("\n\nERROR: $nbErrors found in fork... DIE...\n\n");
	warn Dumper $hErrors;
	die();
}
else {
	warn 'store 1/1: nosql strict-denovo';
	my $nosql = GenBoNoSql->new(dir => $chr->project->getCacheBitVectorDir().'/strict-denovo', mode => 'c');
	$nosql->put_bulk($chr->id(), $hResults);
	$nosql->close();
}



sub check_substitution {
	my ($sam, $var, $limit) = @_;
	my $chr_name = $chr->id();
	my $start = $var->start();
	my $end = $var->end();
	my $all_var = $var->var_allele();
	my $count = 0;
	my $nb_mut = 0;
	my $callback = sub {
		my ($seqid, $pos1, $pileups) = @_;
		return if ($pos1 ne $start);
		if (scalar(@$pileups) < 3) {
			$nb_mut = 99;
			return;
		}
		foreach my $pileup (@$pileups){
			my $b     = $pileup->alignment;
			my $qbase  = substr($b->qseq,$pileup->qpos,1);
			$nb_mut ++ if ($qbase eq $all_var);
			last if ($nb_mut >= $limit);
		}
	};
	$sam->fast_pileup("chr$chr_name:$start-$end", $callback);
	return if ($nb_mut >= $limit);
	return 1;
}

sub check_del {
	my ($bam_file, $type, $locus, $limit) = @_;
	my $cmd = $buffer->getSoftware('sambamba')." depth base $bam_file -L $locus";
	my $res = `$cmd`;
	my @lRes = split("\n", $res);
	my $name_to_check = 'DEL';
	my $col_to_check;
	my $i = 0;
	foreach my $name (split(' ', shift(@lRes))) {
		if ($name eq $name_to_check) {
			$col_to_check = $i;
			last;
		}
		else { $i++; }
	}
	foreach my $res_pos (@lRes) {
		my @lCol = split(' ', $res_pos);
		my $value = $lCol[$col_to_check];
		return if ($value >= $limit)
	}
	return 1;
}

sub check_ins {
	my ($bam_file, $type, $locus, $limit) = @_;
	my $cmd = $buffer->getSoftware('samtools')." mpileup $bam_file -r $locus";
	my $res = `$cmd`;
	my @lRes = split("\n", $res);
	my $count = 0;
	foreach my $line (@lRes) {
		my @lCol = split(' ', $line);
		my $this_count = ($lCol[4] =~ tr/\+//);
		$count += $this_count;
		return if ($count >= $limit);
	}
	return 1;
}