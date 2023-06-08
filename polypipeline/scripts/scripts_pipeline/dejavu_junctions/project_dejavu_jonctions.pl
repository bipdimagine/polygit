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

GetOptions(
	'project=s'    => \$name,
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

foreach my $this_patient (@{$project->getPatients()}) {
	if (($hType_patients and exists $hType_patients->{$this_patient->name()}->{pat}) or not $hType_patients) {
		my @lJunctions = @{$this_patient->getJunctions()};
		foreach my $junction (@lJunctions) {
			$junction->getPatients();
			my $type = $junction->getTypeDescription($this_patient);
			my $chr_id = $junction->getChromosome->id();
			my $start = $junction->start();
			$start =~ s/ //g;
			my $end = $junction->end();
			$end =~ s/ //g;
			my $gene_name = $junction->annex->{$this_patient->name()}->{ensid};
			my $count_new_junction = $junction->get_nb_new_count($this_patient);
			my $count_normal_junction = $junction->get_nb_normal_count($this_patient);
			my $score = int($junction->get_percent_new_count($this_patient));
			my $junction_id = $chr_id.'_'.$start.'_'.$end.'_junction';
			$junction_id =~ s/ //g;
			
			if ($junction->isCanonique($this_patient)) {
				if (not exists $h_proj_junctions_canoniques->{$chr_id}->{$junction_id}) {
					$h_proj_junctions_canoniques->{$chr_id}->{$junction_id}->{start} = $start;
					$h_proj_junctions_canoniques->{$chr_id}->{$junction_id}->{end} = $end;
				}
				$h_proj_junctions_canoniques->{$chr_id}->{$junction_id}->{dejavu}->{$project->name()}->{$this_patient->name()}->{count_junctions} = $count_new_junction;
				$h_proj_junctions_canoniques->{$chr_id}->{$junction_id}->{dejavu}->{$project->name()}->{$this_patient->name()}->{count_normal} = $count_normal_junction;
				$h_proj_junctions_canoniques->{$chr_id}->{$junction_id}->{dejavu}->{$project->name()}->{$this_patient->name()}->{score} = $score;
				$h_proj_junctions_canoniques->{$chr_id}->{$junction_id}->{dejavu}->{$project->name()}->{$this_patient->name()}->{type} = $type;
				$h_proj_junctions_canoniques->{$chr_id}->{$junction_id}->{dejavu}->{$project->name()}->{$this_patient->name()}->{gene_name} = $gene_name;
			}
			else {
				if (not exists $h_proj_junctions->{$chr_id}->{$junction_id}) {
					$h_proj_junctions->{$chr_id}->{$junction_id}->{start} = $start;
					$h_proj_junctions->{$chr_id}->{$junction_id}->{end} = $end;
				}
				$h_proj_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$project->name()}->{$this_patient->name()}->{count_junctions} = $count_new_junction;
				$h_proj_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$project->name()}->{$this_patient->name()}->{count_normal} = $count_normal_junction;
				$h_proj_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$project->name()}->{$this_patient->name()}->{score} = $score;
				$h_proj_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$project->name()}->{$this_patient->name()}->{type} = $type;
				$h_proj_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$project->name()}->{$this_patient->name()}->{gene_name} = $gene_name;
			}
		}
	}
}

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
