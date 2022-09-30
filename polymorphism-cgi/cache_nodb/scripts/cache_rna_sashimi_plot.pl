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

my $hType_patients;
$hType_patients = $project->get_hash_patients_description_rna_seq_junction_analyse() if (-d $project->get_path_rna_seq_junctions_analyse_description_root());


my @lJunctions;
foreach my $chr (@{$project->getChromosomes()}) {
	my $vector_patient = $patient->getJunctionsVector($chr);
	foreach my $junction (@{$chr->getListVarObjects($vector_patient)}) {
		next if ($junction->isCanonique($patient));
		next if ($junction->get_ratio_new_count($patient) == 1);
		#next if ($junction->get_nb_new_count($patient) < 3);
		next if ($junction->get_dp_count($patient) < 3);
		next if ($junction->length() > 10000);
		next if ($junction->get_percent_new_count($patient) < 20);
		$junction->dejavu_percent_coordinate_similar(96);
		my $nb_dejavu_run = $junction->dejavu_nb_int_this_run_patients(20);
		next if ($nb_dejavu_run >= 4);
		#my $nb_dejavu_pat = $junction->dejavu_nb_patients() - $nb_dejavu_run;
		my $nb_dejavu_pat = $junction->dejavu_nb_patients();
		next if ($nb_dejavu_pat > 30);
		push(@lJunctions, $junction);
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
