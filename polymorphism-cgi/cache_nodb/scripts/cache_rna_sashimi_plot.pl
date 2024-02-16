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


my $score = 0;
my $NB_MAX_PLOTS = 300;

if (not $project->is_human_genome()) {
	print FILE 'PASS not Human release';
	close (FILE);
	exit(0); 
}



my $h_junctions_already_scored;
my $pm1       = new Parallel::ForkManager($fork);
my $nbErrors = 0;
$pm1->run_on_finish(
	sub {
		my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $hres ) = @_;
		unless ( defined($hres) or $exit_code > 0 ) {
			$nbErrors++;
			print qq|No message received from child process $exit_code $pid!\n|;
			return;
		}
		delete $hres->{toto};
		foreach my $score (keys %$hres) {
			foreach my $chr_id (keys %{$hres->{$score}}) {
				foreach my $vid (keys %{$hres->{$score}->{$chr_id}}) {
					$h_junctions_already_scored->{$score}->{$chr_id}->{$vid} = undef;
				}	
			}
		}
	}
);

foreach my $chr (@{$project->getChromosomes()}) {
	my $pid = $pm1->start and next;
	$buffer->dbh_deconnect();
	$buffer->dbh_reconnect();
	my $hres;
	my $nb_done;
	$hres->{toto} = undef;
	foreach my $junction (@{$chr->getListVarObjects($patient->getJunctionsVector($chr))}) {
		my $j_pos = $chr->id().'_'.$junction->start().'_'.$junction->end();
		next if ($junction->isCanonique());
		my $j_score_with_dv;
		my $j_score = $junction->junction_score_without_dejavu_global($patient);
		next if ($j_score < $score);
		my $score_f = $junction->junction_score($patient);
		next if ($score_f < $score);
		$hres->{$score_f}->{$chr->id}->{$junction->vector_id()} = undef;
		$nb_done++;
	}
	$chr = undef;
	$pm1->finish(0,$hres);
}
$pm1->wait_all_children();

my $i = 0;
my @lJunctions;
foreach my $score (sort {$b <=> $a} keys %$h_junctions_already_scored) {
	foreach my $chr_id (sort {$a <=> $b} keys %{$h_junctions_already_scored->{$score}}) {
		foreach my $vector_id (sort {$a <=> $b} keys %{$h_junctions_already_scored->{$score}->{$chr_id}}) {
			my $junction = $project->getChromosome($chr_id)->getVarObject($vector_id);
			push(@lJunctions, $junction);
			$i++;
			last if $i == $NB_MAX_PLOTS;
		}
		last if $i == $NB_MAX_PLOTS;
	}
	last if $i == $NB_MAX_PLOTS;
}

my $pm = new Parallel::ForkManager($fork);
foreach my $junction (@lJunctions) {
	my $pid = $pm->start and next;
	$buffer->dbh_deconnect();
	$buffer->dbh_reconnect();
	$junction->can_create_sashimi_plots(1);
	$junction->getListSashimiPlotsPathFiles($patient);
	print FILE 'Ok junction '.$junction->id()."\n";
	$pm->finish();
}
$pm->wait_all_children();


if (scalar(@lJunctions) == 0) {
	print FILE '0 junction';
}

close (FILE);
exit(0);
