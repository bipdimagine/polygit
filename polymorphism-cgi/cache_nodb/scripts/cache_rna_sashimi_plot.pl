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
use Set::IntSpan::Fast::XS;

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
my $bam_file = $patient->getBamFiles->[0];
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
	my $vector_junctions = $patient->getJunctionsVector($chr);
	if (exists $chr->patients_categories->{$patient->name.'_ratio_10'}) {
		$vector_junctions &= $chr->patients_categories->{$patient->name.'_ratio_10'};
	}
	if (exists $chr->global_categories->{'dejavu_20_r10'}) {
		$vector_junctions &= $chr->global_categories->{'dejavu_20_r10'};
	}
	foreach my $junction (@{$chr->getListVarObjects($vector_junctions)}) {
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


my $pipeline_dir = $project->project_pipeline_path();


my $i = 0;
my (@lJunctions, $h_intspan);
foreach my $score (sort {$b <=> $a} keys %$h_junctions_already_scored) {
	foreach my $chr_id (sort {$a <=> $b} keys %{$h_junctions_already_scored->{$score}}) {
		my $chr = $project->getChromosome($chr_id);
		foreach my $vector_id (sort {$a <=> $b} keys %{$h_junctions_already_scored->{$score}->{$chr_id}}) {
			my $junction = $chr->getVarObject($vector_id);
			push(@lJunctions, $junction);
			$h_intspan->{$chr->fasta_name()} = Set::IntSpan::Fast::XS->new() if not exists $h_intspan->{$chr->fasta_name()} ;
			my $s = ($junction->start() - 10000);
			$s = 1 if $s < 1;
			my $e = ($junction->end() + 10000);
			$e = $chr->length() -1 if $junction->end() >= $chr->length();
			$h_intspan->{$chr->fasta_name()}->add_range($s, $e);
			$i++;
			last if $i == $NB_MAX_PLOTS;
		}
		last if $i == $NB_MAX_PLOTS;
	}
	last if $i == $NB_MAX_PLOTS;
}


my $bam_tmp = $pipeline_dir."/".$patient_name.".tmp.bam";
my $bam_tmp_sort = $pipeline_dir."/".$patient_name.".sort.tmp.bam";
my $bam_rmdup = $pipeline_dir."/".$patient_name.".rmdup.bam";

if (not -e $bam_tmp_sort) { 
	chdir($pipeline_dir);
	
	my $bed_file = $pipeline_dir."/".$patient_name.".regions.bed";
	open (BED, ">$bed_file");
	foreach my $chr_id (sort keys %$h_intspan) {
		foreach my $in (split(',', $h_intspan->{$chr_id}->as_string())) {
			my ($s, $e) = split('-', $in);
			print BED "$chr_id\t$s\t$e\n";
		}
	}
	close (BED);
	
	my $bam_regions = $pipeline_dir."/".$patient_name.".regions.bam";
	my $cmd_regions = $buffer->software('samtools')." view -@ $fork -b -M -L $bed_file $bam_file >$bam_regions";
#	warn "\nR/\n".$cmd_regions;
	`$cmd_regions`;
	
	my $bam_sort = $pipeline_dir."/".$patient_name.".sort.bam";
	my $cmd_bam = $buffer->software('samtools')." sort -@ $fork -n $bam_regions >$bam_sort";
#	warn "\n0/\n".$cmd_bam;
	`$cmd_bam`;
	
	my $cmd_bam_1a = $buffer->software('samtools')." fixmate -@ $fork -m -O bam $bam_sort $bam_tmp";
#	warn "\n1A/\n".$cmd_bam_1a;
	`$cmd_bam_1a`;
	
	my $cmd_bam_1b = $buffer->software('samtools')." sort -@ $fork $bam_tmp >$bam_tmp_sort";
#	warn "\n1b/\n".$cmd_bam_1b;
	`$cmd_bam_1b`;
}

if (not -e $bam_rmdup) { 
	my $cmd_bam_2 = $buffer->software('samtools')." markdup -@ $fork -r ".$bam_tmp_sort." ".$bam_rmdup;
#	warn "\n2/\n".$cmd_bam_2;
	`$cmd_bam_2`;
	my $cmd_bam_3 = $buffer->software('samtools')." index $bam_rmdup";
#	warn "\n3/\n".$cmd_bam_3;
	`$cmd_bam_3`;
}



my $pm = new Parallel::ForkManager($fork);
foreach my $junction (@lJunctions) {
	my $pid = $pm->start and next;
	$buffer->dbh_deconnect();
	$buffer->dbh_reconnect();
#	warn $junction->id;
	$junction->can_create_sashimi_plots(1);
	$junction->getListSashimiPlotsPathFiles($patient, $bam_rmdup);
	print FILE 'Ok junction '.$junction->id()."\n";
	$pm->finish();
}
$pm->wait_all_children();

if (-e $bam_tmp or -e $bam_rmdup) {
	my $cmd_rm_bam_tmp1 = "rm $pipeline_dir/$patient_name.*.bed";
	my $cmd_rm_bam_tmp2 = "rm $pipeline_dir/$patient_name.*.bam";
	my $cmd_rm_bam_tmp3 = "rm $pipeline_dir/$patient_name.*.bai";
	`$cmd_rm_bam_tmp1`;
	`$cmd_rm_bam_tmp2`;
	`$cmd_rm_bam_tmp3`;
}


if (scalar(@lJunctions) == 0) {
	print FILE '0 junction';
}

close (FILE);
exit(0);


