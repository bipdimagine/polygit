#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use Getopt::Long;
use Data::Dumper;
use Parallel::ForkManager;
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use GBuffer;

my $fork = 1;
my ($project_name, $chr_name);
GetOptions(
	'project=s' => \$project_name,
	'chr=s'     => \$chr_name,
	'fork=s'     => \$fork,
);
die("\n\nERROR: -project option not defined. Die...\n\n") unless ($project_name);
die("\n\nERROR: -chr option not defined. Die...\n\n") unless ($chr_name);

my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => $project_name );

my $nb_ok = 0;
my $nb_error = 0;
my $chr = $project->getChromosome($chr_name);
my $chr_name = 'chr'.$chr->id();
warn "\n\n### CHECKING $chr_name\n";
my $pm = new Parallel::ForkManager($fork);
$pm->run_on_finish(
	sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hRes) = @_;
		unless (defined($hRes) or $exit_code > 0) {
			print qq|No message received from child process $exit_code $pid!\n|;
			return;
		}
		$nb_ok += $hRes->{'ok'};
		$nb_error += $hRes->{'error'};
		if ($nb_error > 0) {
			my $this_patient_name = $hRes->{patient_name};
			die("\n\nERROR: found error QueryVcf for $this_patient_name in $chr_name ! Die...\n\n");
		}
	}
);

my $hPat;
foreach my $patient (@{$chr->getPatients()}) {
	my $v = $patient->getVariantsVector($chr)->Clone();
	$v->Intersection($v, $chr->global_categories->{substitution});
	$hPat->{$patient->name()}->{nb_cache} = $chr->countThisVariants($v);
	$hPat->{$patient->name()}->{vcf_file} = $patient->getVariationsFiles->[0];
}

foreach my $patient_name (keys %{$hPat}) {
	my $pid = $pm->start and next;
	my ($hRes, $nb_sub_cache, $nb_sub_bcftools);
	$hRes->{patient_name} = $patient_name;
	$nb_sub_cache = $hPat->{$patient_name}->{nb_cache};
	my $file_vcf = $hPat->{$patient_name}->{vcf_file};
	my $cmd = "bcftools norm --multiallelics - $file_vcf | bcftools view --types snps - | bcftools view --exclude 'GT=\"0/0\"' - | bcftools view --exclude 'GT=\"./.\"' - | bcftools stats -t '$chr_name' - | grep 'number of SNPs:'";
	my $res = `$cmd`;
	my @lTmp = split(' ', $res);
	$nb_sub_bcftools = $lTmp[-1];
	my $is_ok = 'ERROR';
	if ($nb_sub_cache == $nb_sub_bcftools) {
		$is_ok = 'OK';
		$hRes->{'ok'}++;
	}
	else { $hRes->{'error'}++; }
	warn  "[$is_ok] patient $patient_name\tCache:$nb_sub_cache\tBcftools:$nb_sub_bcftools\n";
	$pm->finish( 0, $hRes );
}
$pm->wait_all_children();
warn "\n\n### RESUME\n";
warn "Nb OK:\t$nb_ok\n";
warn "Nb ERROR:\t$nb_error\n";
exit(0);
