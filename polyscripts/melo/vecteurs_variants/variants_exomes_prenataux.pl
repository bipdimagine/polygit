#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
#print "$RealBin\n";
use lib "$RealBin";
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use lib "$RealBin/../../../polypipeline/dragen/scripts/";
use lib "$RealBin/../../../polypipeline/packages/";
use dragen_util; 
use file_util; 
use Data::Dumper;
use Getopt::Long;
use Carp;
use IO::Prompt;
use Net::SFTP;
use File::Find::Rule;
use Parallel::ForkManager;
use GBuffer;

my $projects_name = 'NGS2024_8172';
my $fork = 1; # nb de job en parallele
GetOptions(
	'projects=s'		=> \$projects_name,
#	'patient=s'			=> \$patient_name,
	'fork=i'			=> \$fork,
);

$fork = 40 if ($fork > 40 or $fork <= 0);

my $nb_var;

my $file = "/home/mperin/git/polygit/polyscripts/melo/vecteurs_variants/exomes_prenat-variants-2.csv";
warn $file;
open(my $fh, ">$file") or confess("Can't open '$file'");
print {$fh} "project;family;child;HGMD DM;ClinVar pathogenic;CADD > 25;denovo;strict denovo;recessif;compound;freq < 1%;DéjàVu < 50 (sample);";
print {$fh} "freq < 1 & DéjàVu < 50;HGMD DM & freq < 1 & DéjàVu < 50;ClinVar pathogenic & freq < 1 & DéjàVu < 50;";
print {$fh} "CADD > 25 & freq < 1 & DéjàVu < 50;denovo & freq < 1 & DéjàVu < 50;strict denovo & freq < 1 & DéjàVu < 50;recessif & freq < 1 & DéjàVu < 50;compound & freq < 1 & DéjàVu < 50;";
#print {$fh} "Union = DM OR patho OR CADD > 25 OR denovo OR recessif OR compound AND freq < 1% AND DéjàVu < 50 (sample);";
#print {$fh} "Intersection = DM AND patho AND CADD > 25 AND (denovo OR recessif OR compound) AND freq < 1% AND DéjàVu < 50 (sample)\n";
print {$fh} "(DM OR patho OR CADD > 25) AND (strict denovo OR recessif OR compound) AND freq < 1% AND DéjàVu < 50;";
print {$fh} "Comment\n";
close ($fh);

my @projects_name = split /,\s*/, $projects_name;
warn scalar @projects_name.' project(s)';

my $pm = new Parallel::ForkManager($fork);
foreach my $project_name (@projects_name) {
	my $pid = $pm->start and next;
	my $buffer = new GBuffer;
	my $project = $buffer->newProjectCache(-name=>$project_name) or confess ("can\'t open project cache '$project_name'");
	warn $project->name;
	
	open(my $fh, ">>$file") or confess("Can't open '$file'");
#	open(my $fh, ">/home/mperin/git/polygit/polyscripts/melo/vecteurs_variants/$project_name-variants.csv") or confess("Can't open 'exomes_prenat-variants.csv'");

	my $patients;
	my $families = $project->getFamilies;
	warn scalar @$families.' families';
	foreach my $fam (@$families) {
#		next unless ($fam->name eq 'TOU');
		warn $fam->name;
		my @pat = keys(%{$fam->children_ill});
		warn ("no ill child found: $project_name, ".$fam->name) unless (scalar @pat);
		next unless (scalar @pat);
		foreach my $child (@pat) {
			my $pat = $project->getPatient($child);
			my $pat_name = $pat->name;
			warn $pat_name;
			$nb_var->{$project_name}->{$pat_name}->{'commentaire'} = 'Missing both parents' unless ($fam->mother || $fam->father);
			$nb_var->{$project_name}->{$pat_name}->{'commentaire'} = 'Missing mother' if ($fam->father && not $fam->mother);
			$nb_var->{$project_name}->{$pat_name}->{'commentaire'} = 'Missing father' if ($fam->mother && not $fam->father);
			warn $nb_var->{$project_name}->{$pat_name}->{'commentaire'} if (exists $nb_var->{$project_name}->{$pat_name}->{'commentaire'});
			foreach my $chr (@{$project->getChromosomes}) {
	#			next unless ($chr->name eq '1');
	#			next unless ($chr->name eq 'MT' or $chr->name eq '1');
				warn $chr->name;
				next unless ($chr->getNewVector()->Size);
	
				my $v0 = $pat->getVectorVariants($chr);
				my $v1 = $v0 & $chr->vectorDM;
				my $v2 = $v0 & $chr->vectorClinvarPathogenic;
	#			warn Dumper $project->hash_cadd_filters;
				my $v3 = $v0 & $chr->get_vector_cadd_variants({'cadd_25'=>'1'});# unless (scalar @{$chr->getGenes} == 0);
				my $v4 = $fam->getVectorDenovoTransmission($chr, $pat);
				my $v8 = $fam->getVectorStrictDenovoTransmission($chr, $pat);
				my $v5 = $fam->getModelVector_fam_recessif($chr, $pat);
				my $v6 = $chr->getNewVector;
#				if ($fam->isTrio) {
					foreach my $gene (@{$chr->getGenes}) {
						$v6 += $fam->getModelVector_fam_compound($gene);
					}
#				}
				my $var_freq = $chr->getNewVector();
#				warn Dumper $chr->vector_global_categories;
				foreach my $filter_name ('freq_01', 'freq_001', 'freq_0001', 'freq_none') {
#					warn $filter_name.' -> '.$chr->vector_global_categories->{$filter_name}->Norm();
				    $var_freq += $chr->vector_global_categories->{$filter_name};
				}
#				warn '$var_freq -> '.$var_freq->Norm."\t".'&= $v0 -> '.($v0 & $var_freq)->Norm;
				$var_freq &= $v0;
	#			my $v8 = $chr->get_vector_dejavu($v0, 50);
				my $v_dv = $v0 & $chr->lmdb_score_impact->get('sdv_50');
	#			warn '$v_dv -> '.$v_dv->Norm."\t".'&= $v0 -> '.($v0 & $v_dv)->Norm;
				
				$nb_var->{$project_name}->{$pat_name}->{'DM'} += $v1->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'ClinVar patho'} += $v2->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'CADD 25'} += $v3->Norm if ($v3);
				$nb_var->{$project_name}->{$pat_name}->{'denovo'} += $v4->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'strict denovo'} += $v8->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'recessif'} += $v5->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'compound'} += $v6->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'freq 1%'} += ($v0 & $var_freq)->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'dejavu'} += $v_dv->Norm;
				my $v_f = $var_freq & $v_dv;
				$nb_var->{$project_name}->{$pat_name}->{'freq & dejavu'} += $v_f->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'DM & filtres'} += ($v1 & $v_f)->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'Clinvar & filtres'} += ($v2 & $v_f)->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'CADD & filtres'} += ($v3 & $v_f)->Norm if ($v3);
				$nb_var->{$project_name}->{$pat_name}->{'denovo & filtres'} += ($v4 & $v_f)->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'strict denovo & filtres'} += ($v8 & $v_f)->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'recessif & filtres'} += ($v5 & $v_f)->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'compound & filtres'} += ($v6 & $v_f)->Norm;
				
				my $v_union = $v1 + $v2 + $v4 + $v8 + $v5 + $v6 ;
				$v_union += $v3 if ($v3);
				$v_union &= $v0 & $var_freq & $v_dv;
	#			$v_union &= $v_dv if ($v_dv);
	#			$v_union = $chr->get_vector_dejavu($v_union, 50);
				
				my $v_intersec = $v4 + $v8 + $v5 + $v6;
				$v_intersec &= $v0 & $v1 & $v2 & $var_freq & $v_dv ;
				$v_intersec &= $v3 if ($v3);
	#			$v_intersec &= $v_dv if ($v_dv);
	#			$v_intersec = $chr->get_vector_dejavu($v_intersec, 50);

				my $v_impact = $v1 + $v2;
				$v_impact += $v3 if ($v3);
				my $v_transmission = $v8 + $v5 + $v6;
				my $v_ensemble = $v0 & $v_impact & $v_transmission & $var_freq & $v_dv;
				
				$nb_var->{$project_name}->{$pat_name}->{'union'} += $v_union->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'intersection'} += $v_intersec->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'ensemble'} += $v_ensemble->Norm;
				
			}
			print {$fh} $project_name.';'.$fam->name.';'.$pat_name.';'
				.$nb_var->{$project_name}->{$pat_name}->{'DM'}.';'
				.$nb_var->{$project_name}->{$pat_name}->{'ClinVar patho'}.';'
				.$nb_var->{$project_name}->{$pat_name}->{'CADD 25'}.';'
				.$nb_var->{$project_name}->{$pat_name}->{'denovo'}.';'
				.$nb_var->{$project_name}->{$pat_name}->{'strict denovo'}.';'
				.$nb_var->{$project_name}->{$pat_name}->{'recessif'}.';'
				.$nb_var->{$project_name}->{$pat_name}->{'compound'}.';'
				.$nb_var->{$project_name}->{$pat_name}->{'freq 1%'}.';'
				.$nb_var->{$project_name}->{$pat_name}->{'dejavu'}.';'
#				.$nb_var->{$project_name}->{$pat_name}->{'union'}.';'
#				.$nb_var->{$project_name}->{$pat_name}->{'intersection'}."\n";
				.$nb_var->{$project_name}->{$pat_name}->{'freq & dejavu'}.';'
				.$nb_var->{$project_name}->{$pat_name}->{'DM & filtres'}.';'
				.$nb_var->{$project_name}->{$pat_name}->{'Clinvar & filtres'}.';'
				.$nb_var->{$project_name}->{$pat_name}->{'CADD & filtres'}.';'
				.$nb_var->{$project_name}->{$pat_name}->{'denovo & filtres'}.';'
				.$nb_var->{$project_name}->{$pat_name}->{'strict denovo & filtres'}.';'
				.$nb_var->{$project_name}->{$pat_name}->{'recessif & filtres'}.';'
				.$nb_var->{$project_name}->{$pat_name}->{'compound & filtres'}.';'
				.$nb_var->{$project_name}->{$pat_name}->{'ensemble'}.';'
				.$nb_var->{$project_name}->{$pat_name}->{'commentaire'};
			print {$fh} "\n";
		}
	}
	close ($fh);
	warn Dumper $nb_var->{$project_name};
	
	$project->disconnect;
	$project = undef;
	$buffer = undef;
	$pm->finish;
}
$pm->wait_all_children();

#close $fh;
warn Dumper $nb_var;









