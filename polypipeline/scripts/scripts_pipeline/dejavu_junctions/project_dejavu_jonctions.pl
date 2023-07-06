#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use JSON;
use Data::Dumper;
use CGI qw/:standard :html3/;
use Set::IntSpan::Fast::XS;
use Set::IntervalTree;
use List::Util qw[min max];
use FindBin qw($Bin);
use Storable qw(store retrieve freeze);
use Number::Format qw(:subs);
use File::Basename;
use Getopt::Long;
use JSON;

use GBuffer;
use GenBoProject;
use GenBoCache;
use GenBoNoSqlIntervalTree;


my $cgi = new CGI;
my ($name);
my $fork = 20;

GetOptions(
	'project=s'    => \$name,
	'fork=s'    => \$fork,
);


my $hash_project_phenotypes;
my $b1 = GBuffer->new;
my $project = $b1->newProjectCache( -name => $name );


if (-d $project->get_path_rna_seq_junctions_analyse_all_res()) {
	my $ok;
	my $path = $project->get_path_rna_seq_junctions_analyse_all_res();
	my $se_file = $path.'/allResSE.txt';
	my $ri_file = $path.'/allResRI.txt';
	$ok = 1 if (-e $se_file);
	$ok = 1 if (-e $ri_file);
	if ($ok) {
		if ($project->is_human_genome()) {
			foreach my $pheno_name (@{$project->phenotypes()}) {
				$pheno_name =~ s/ /_/g;
				$hash_project_phenotypes->{lc($pheno_name)}->{$name} = undef;
			}
		}
	}
}
my $ok_dragen;
foreach my $patient (@{$project->getPatients()}) { 
	my $dragen_file = $project->getVariationsDir('dragen-sj').'/'.$patient->name().'.SJ.out.tab.gz';
	$ok_dragen = 1 if (-e $dragen_file);
	my $star_file = $project->getJunctionsDir('star').'/'.$patient->name().'.SJ.tab.gz';
	$ok_dragen = 1 if (-e $star_file);
}
if ($ok_dragen) {
		if ($project->is_human_genome()) {
		foreach my $pheno_name (@{$project->phenotypes()}) {
			$pheno_name =~ s/ /_/g;
			$hash_project_phenotypes->{lc($pheno_name)}->{$name} = undef;
		}
	}
}

my $dv_dir_path = $project->DejaVuJunction_path();
my $dir_dv_proj = $dv_dir_path.'/projects/';
unless (-d $dir_dv_proj) {
	`mkdir $dir_dv_proj`;
	`chmod 777 $dir_dv_proj`;
}
$dir_dv_proj .= $name;
unless (-d $dir_dv_proj) {
	`mkdir $dir_dv_proj`;
	`chmod 777 $dir_dv_proj`;
}

print "LAUNCH $name\n";
my ($h_proj_junctions, $h_proj_junctions_canoniques);


my $hType_patients;
$hType_patients = $project->get_hash_patients_description_rna_seq_junction_analyse() if (-d $project->get_path_rna_seq_junctions_analyse_description_root());

my $h_vector;
my $pm = new Parallel::ForkManager($fork);
my $nbErrors = 0;
$pm->run_on_finish(
	sub {
		my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $hres ) = @_;
		unless ( defined($hres) or $exit_code > 0 or not exists $hres->{done} ) {
			$nbErrors++;
			print qq|No message received from child process $exit_code $pid!\n|;
			warn Dumper $hres;
			return;
		}
		delete $hres->{done};
		
		my $chr_id = $hres->{'chr'};
		foreach my $jid (keys %{$hres->{junctions}}) {
			$h_proj_junctions->{$chr_id}->{$jid} = $hres->{junctions}->{$jid};
		}
		foreach my $jid (keys %{$hres->{canoniques}}) {
			$h_proj_junctions_canoniques->{$chr_id}->{$jid} = $hres->{canoniques}->{$jid};
		}
	}
);


foreach my $chr (@{$project->getChromosomes}) {
	$pm->start and next;
	my $chr_id = $chr->id();
	my $hres;
	$hres->{'chr'} = $chr_id;
	my @lJunctions = @{$chr->getJunctions()};
	my ($h_proj_junctions_canoniques_tmp, $h_proj_junctions_tmp);
	foreach my $junction (@lJunctions) {
		my $start = $junction->start();
		$start =~ s/ //g;
		my $end = $junction->end();
		$end =~ s/ //g;
		my $junction_id = $chr_id.'_'.$start.'_'.$end.'_junction';
		$junction_id =~ s/ //g;
		foreach my $pat_j (@{$junction->getPatients()}) {
			my $type = $junction->getTypeDescription($pat_j);
			my $gene_name = $junction->annex->{$pat_j->name()}->{ensid};
			my $count_new_junction = $junction->get_nb_new_count($pat_j);
			my $count_normal_junction = $junction->get_nb_normal_count($pat_j);
			my $score = int($junction->get_percent_new_count($pat_j));
			if ($junction->isCanonique()) {
				if (not exists $h_proj_junctions_canoniques_tmp->{$junction_id}) {
					$h_proj_junctions_canoniques_tmp->{$junction_id}->{start} = $start;
					$h_proj_junctions_canoniques_tmp->{$junction_id}->{end} = $end;
				}
				$h_proj_junctions_canoniques_tmp->{$junction_id}->{dejavu}->{$project->name()}->{$pat_j->name()}->{count_junctions} = $count_new_junction;
				$h_proj_junctions_canoniques_tmp->{$junction_id}->{dejavu}->{$project->name()}->{$pat_j->name()}->{count_normal} = $count_normal_junction;
				$h_proj_junctions_canoniques_tmp->{$junction_id}->{dejavu}->{$project->name()}->{$pat_j->name()}->{score} = $score;
				$h_proj_junctions_canoniques_tmp->{$junction_id}->{dejavu}->{$project->name()}->{$pat_j->name()}->{type} = $type;
				$h_proj_junctions_canoniques_tmp->{$junction_id}->{dejavu}->{$project->name()}->{$pat_j->name()}->{gene_name} = $gene_name;
			}
			else {
				if (not exists $h_proj_junctions_tmp->{$junction_id}) {
					$h_proj_junctions_tmp->{$junction_id}->{start} = $start;
					$h_proj_junctions_tmp->{$junction_id}->{end} = $end;
				}
				$h_proj_junctions_tmp->{$junction_id}->{dejavu}->{$project->name()}->{$pat_j->name()}->{count_junctions} = $count_new_junction;
				$h_proj_junctions_tmp->{$junction_id}->{dejavu}->{$project->name()}->{$pat_j->name()}->{count_normal} = $count_normal_junction;
				$h_proj_junctions_tmp->{$junction_id}->{dejavu}->{$project->name()}->{$pat_j->name()}->{score} = $score;
				$h_proj_junctions_tmp->{$junction_id}->{dejavu}->{$project->name()}->{$pat_j->name()}->{type} = $type;
				$h_proj_junctions_tmp->{$junction_id}->{dejavu}->{$project->name()}->{$pat_j->name()}->{gene_name} = $gene_name;
			}
		}
	}
	$hres->{canoniques} = $h_proj_junctions_canoniques_tmp;
	$hres->{junctions} = $h_proj_junctions_tmp;
	$hres->{done} = 1;
	$pm->finish( 0, $hres );
}
$pm->wait_all_children();
die if $nbErrors > 0;

my $json_encode = encode_json $h_proj_junctions;
open (JSON1, '>'.$dir_dv_proj.'/'.$name.'.json');
print JSON1 $json_encode;
close (JSON1);

open (JSON2, '>'.$dir_dv_proj.'/'.$name.'.canoniques.json');
if ($h_proj_junctions_canoniques) {
	my $json_c_encode = encode_json $h_proj_junctions_canoniques;
	print JSON2 $json_c_encode;
}
close (JSON2);

print "$dir_dv_proj/ -> Done!\n";
print "$name -> Done!\n";
