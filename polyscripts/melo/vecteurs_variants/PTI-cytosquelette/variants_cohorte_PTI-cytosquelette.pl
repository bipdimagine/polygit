#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
#print "$RealBin\n";
use lib "$RealBin";
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
use lib "$RealBin/../../../../polypipeline/dragen/scripts/";
use lib "$RealBin/../../../../polypipeline/packages/";
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


# Faris ZECO
# cohorte de patients PTI-cytosquelette

# CADD >= 20 
# spliceAI >= 0.5
# MAF < 0.01%
# vues dans gnomAD: He < 10; Ho = 0
# tous types de variants
# uniquement les patients (pas de parents si trios)


my $splice_AI = 0.5;
my $dv_projects = 1;
my $dv_pat = 5;
my $dv_pat_ho = 0;

my $projects_file = '/home/mperin/git/polygit/polyscripts/melo/vecteurs_variants/PTI-cytosquelette/Liste_projets_Polyweb.csv';
my $projects_file = '/home/mperin/git/polygit/polyscripts/melo/vecteurs_variants/PTI-cytosquelette/Projets_polyweb_PTI_cytosquelette_round_2.csv';
my $projects_file = '/home/mperin/git/polygit/polyscripts/melo/vecteurs_variants/PTI-cytosquelette/Projets_polyweb_PTI_cytosquelette_09_2025.txt';
my $genes_file = '/home/mperin/git/polygit/polyscripts/melo/vecteurs_variants/PTI-cytosquelette/Gènes_cytosquelette.csv';
my $genes_file = '/home/mperin/git/polygit/polyscripts/melo/vecteurs_variants/PTI-cytosquelette/Liste_complète_gènes_cytosquelette.txt';
my $genes_file = '/home/mperin/git/polygit/polyscripts/melo/vecteurs_variants/PTI-cytosquelette/All_genes_final_list_01_10_25.csv';
my $fork = 1;
GetOptions(
	'projects_file=s'		=> \$projects_file,
	'genes_file=s'			=> \$genes_file,
	'fork=i'				=> \$fork,
) || die ("Error in command line arguments\n");

die ("fork should be [1;40], given $fork") if ($fork > 40 or $fork <= 0);
warn 'fork='.$fork;

# Charge la liste des noms de projets
die ("Can't find file '$projects_file'") unless (-e $projects_file);
warn $projects_file;
my @projects_list;
open(my $fh_projects, "<", "$projects_file") || confess("Can't open '$projects_file': $!");
while (my $l = readline($fh_projects)) {
	chomp $l;
	next unless ($l);
	push(@projects_list, $l);
}
close($fh_projects);
warn scalar @projects_list.' project(s) in list';
#warn Dumper \@projects_list;

# Charge la liste des correctifs de noms de gènes
my $correction_file = '/home/mperin/git/polygit/polyscripts/melo/vecteurs_variants/PTI-cytosquelette/correction_gene_names.tsv';
open(my $fh_correction, "<", "$correction_file") || confess("Can't open '$correction_file': $!");
my $correction;
while (my $l = readline($fh_correction)) {
	next if ($l =~ /^#/);
	chomp $l;
	my @l = split("\t",$l);
	$correction->{$l[0]} = $l[1];
}
close($fh_correction);
#warn Dumper $correction;

# Charge le noms de gènes, et les corrige si besoin
die ("Can't find file '$genes_file'") unless (-e $genes_file);
warn $genes_file;
my @genes_list;
open(my $fh_genes, "<", "$genes_file") || confess("Can't open '$genes_file': $!");
while (my $l = readline($fh_genes)) {
	chomp $l;
	next unless ($l);
	my @grep = grep{/^$l$/i} keys %$correction;
	$l = $correction->{$grep[0]} if (scalar @grep);
	push(@genes_list, $l);
}
close($fh_genes);
warn scalar @genes_list.' gene(s) in list';
#warn Dumper \@genes_list;


# Ecrit le header des résultats
my $file = "/home/mperin/git/polygit/polyscripts/melo/vecteurs_variants/PTI-cytosquelette/variants.csv";
warn $file;
open(my $fh, ">$file") or confess("Can't open '$file'");
print {$fh} "project;family;child;";
print {$fh} "(CADD > 25 OR spliceAI > $splice_AI) AND freq < 0.01% AND gnomAD He < 5 AND gnomAD Ho = 0 AND in gene list AND dejavu other projects < $dv_projects AND DejaVu other patients < $dv_pat AND DejaVu other patients ho < $dv_pat_ho;";
print {$fh} "id, gene(s), nb all ref, nb all alt, depth;";
print {$fh} "link(s)";
print {$fh} ";& cov >= 30";
print {$fh} "\n";
close ($fh);

my $file_genes = $file =~ s/variants(-\d)?\.csv/nb_variants_par_genes$1\.csv\.tmp/r;
warn $file_genes;
open(my $stats_genes, ">$file_genes") or confess("Can't open '$file_genes'");
print {$stats_genes} "project;gene;number of patients;patients\n";
close ($stats_genes);


# log file
my $file_err = $file =~ s/variants.csv$/log.error/r;
warn $file_err;
system ("rm $file_err") if (-e $file_err);



my $nb_var;
my $count_genes;
my $pm = new Parallel::ForkManager($fork);
foreach my $project_name (@projects_list) {
	
	my $pid = $pm->start and next;# unless ($fork <= 1);
	my $buffer = new GBuffer;
	my $project = $buffer->newProjectCache(-name=>$project_name) or confess ("can\'t open project cache '$project_name'");
	warn $project->name;
	
	open(my $fh, ">>$file") or confess("Can't open '$file'");
	open(my $stats_genes, ">>$file_genes") or confess("Can't open '$file_genes'");
	open(my $err, ">>$file_err") or confess("Can't open '$file_err'");
	
	unless ($project->isExome or $project->isGenome) {
		my @analyse_type = map {$_->analyse} @{$project->getCaptures};
		print {$err} "$project_name\t" . join(', ', @analyse_type ) ."\n\n";
		warn "$project_name is not a exome or genome project : ".join(', ', @analyse_type );
		next;
	}
	
	my $h_genes_ok;
	
#	my $genes = $project->newGenes(\@genes_list);
	my $genes;
	my $genes_not_found;
	foreach my $gene_name (@genes_list) {
#		warn $gene_name;
		my $gene;
		eval{
			$gene = $project->newGene($gene_name);
			unless ($gene) {
				push(@$genes_not_found, $gene_name);
				print {$err} "$project_name\t$gene_name not found\tgencode ".$project->buffer->getQuery()->getGencodeVersion( $project->id )."\n\n\n";
				next;
			}
		};
		if ($@) {
			print {$err} "$project_name\t$gene_name\n$@\n\n\n";
			next;
		}
#		warn $gene->external_name."\t".$gene->getChromosome->name;
		push(@$genes, $gene);
	}
	warn "Genes not found:\n". Dumper $genes_not_found if ($genes_not_found); 
#	warn Dumper map {$_->external_name} @$genes;
	
	warn "Phenotypes $project_name: ". join(', ', @{$project->phenotypes});
	print {$err} "No phenotype for $project_name\n\n" unless (scalar @{$project->phenotypes});
	
	my $patients;
	my $families = $project->getFamilies;
	warn scalar @$families.' families';
	foreach my $fam (sort @$families) {
		
		warn $fam->name;
		my @pat = keys(%{$fam->children_ill});
		warn ("no ill child found: $project_name, ".$fam->name) unless (scalar @pat);
#		warn Dumper @pat;
		next unless (scalar @pat);
		foreach my $child (sort @pat) {
			my $pat = $project->getPatient($child);
			my $pat_name = $pat->name;
			warn $pat_name;
			$nb_var->{$project_name}->{$pat_name}->{'commentaire'} = 'Missing both parents' unless ($fam->mother || $fam->father);
			$nb_var->{$project_name}->{$pat_name}->{'commentaire'} = 'Missing mother' if ($fam->father && not $fam->mother);
			$nb_var->{$project_name}->{$pat_name}->{'commentaire'} = 'Missing father' if ($fam->mother && not $fam->father);
			warn $nb_var->{$project_name}->{$pat_name}->{'commentaire'} if (exists $nb_var->{$project_name}->{$pat_name}->{'commentaire'});
			foreach my $chr (@{$project->getChromosomes}) {
				
				warn 'chr '.$chr->name;
				next unless ($chr->getNewVector()->Size);
	
				my $v0 = $pat->getVectorVariants($chr)->Clone();
				
				my $var_cadd = $v0 & $chr->get_vector_cadd_variants({'cadd_25'=>'1'});
				my $var_genomad_ac = $v0 & $chr->getVectorScore('gnomad_ac_5');
				my $var_genomad_ho = $v0 & $chr->getVectorScore('gnomad_ho_ac_0');
				
#				warn $v0->Norm();
#				warn $var_cadd->Norm();
#				warn $var_genomad_ac->Norm();
#				warn $var_genomad_ho->Norm();
				
				my $var_freq = $chr->getNewVector();
				foreach my $filter_name ('freq_0001', 'freq_none') {
				    $var_freq += $chr->vector_global_categories->{$filter_name};
				}
				$var_freq &= $v0;
#				warn $var_freq->Norm();
				
				my $v_genes = $chr->getNewVector();
				foreach my $gene (@$genes) {
					next unless ($gene->getChromosome->name eq $chr->name);
					$h_genes_ok->{$gene->external_name} = $gene->id;
					
#					warn $gene->external_name .' -> from '.$gene->start.' to '.$gene->end;
#					my $this_v = $chr->getVectorByPosition($gene->start,$gene->end);
					my $this_v = $gene->getVectorOrigin();
#					warn $this_v->to_Enum();
#					foreach my $var (@{$chr->getListVarObjects($this_v)}) {
#						warn $var->id().' -> '.$var->vector_id;	
#					}
#					die if $this_v->bit_test(8524);
					
					$v_genes += $this_v;
				}
				$v_genes &= $v0;
#				warn $v_genes->Norm();
				
				my $var_spliceAI_in_genes = $chr->getNewVector();
				foreach my $var (@{$chr->getListVarObjects($v_genes)}) {
					foreach my $gene (@{$var->getGenes}) {
						my @ensg = split(/_/, $gene->id);
						my $motif = '^'.$ensg[0].'$|^'.$gene->external_name.'$';
#						warn $motif;
#						warn grep (/$motif/i, @genes_list);
						if (grep (/$motif/i, @genes_list) and $var->max_spliceAI_score($gene) >= $splice_AI) {
							$var_spliceAI_in_genes->Bit_On($var->vector_id);
						}
					}
				}
				
				my $v_intersection = $var_cadd & $var_genomad_ac & $var_genomad_ho & $var_freq & $v_genes;
#				warn $v_intersection->Norm();
				$v_intersection += $var_spliceAI_in_genes & $var_genomad_ac & $var_genomad_ho & $var_freq & $v_genes;
#				warn $v_intersection->Norm();
				
#				my $i=0;
#				foreach my $v ($v_genes,$var_cadd,$var_spliceAI_in_genes,$var_genomad_ac,$var_genomad_ho,$var_freq) {
#					warn $i;
#					die if $v->bit_test(8524);
#					warn $v->bit_test(8524) if ($chr->name == 6); # vector id 6_122734703_C_T = 8524
#					$i++;
#				}
				
				$nb_var->{$project_name}->{$pat_name}->{'CADD'} += $var_cadd->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'spliceAI'} += $var_spliceAI_in_genes->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'freq'} += $var_freq->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'genomAD he'} += $var_genomad_ac->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'genomAD ho'} += $var_genomad_ho->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'in genes list'} += $v_genes->Norm;
				$nb_var->{$project_name}->{$pat_name}->{'intersection'} += $v_intersection->Norm;
				warn $v_intersection->Norm.' variants';
				
				
				my $infos;
				my $link;
				foreach my $var (@{$chr->getListVarObjects($v_intersection)}) {
					warn $var->id;
#					warn $infos;
#					warn $var->vector_id;
					$nb_var->{$project_name}->{$pat_name}->{'intersection dejavu other projects'} ++ if ($var->other_projects <= $dv_projects);
					$nb_var->{$project_name}->{$pat_name}->{'intersection dejavu other patients'} ++ if ($var->other_patients <= $dv_pat);
					$nb_var->{$project_name}->{$pat_name}->{'intersection dejavu other patients ho'} ++ if ($var->other_patients_ho <= $dv_pat_ho);
					next unless ($var->other_projects <= $dv_projects && $var->other_patients <= $dv_pat && $var->other_patients_ho <= $dv_pat_ho);
					$nb_var->{$project_name}->{$pat_name}->{'intersection dejavu other'} ++;
					
					$infos .= "\n" if ($infos or $nb_var->{$project_name}->{$pat_name}->{'intersection details'});
					$infos .= $var->id.',';
					my $gname;
					foreach my $g (@{$var->getGenes}) {
						my @ensg = split(/_/, $g->id);
						my $motif = '^'.$ensg[0].'$|^'.$g->external_name.'$';
#						warn grep (/$motif/i, @genes_list);
						if (grep (/$motif/i, @genes_list)) {
							$gname = $g->external_name;
							$infos .= $gname.',';
							$count_genes->{'intersection'}->{$project_name}->{$gname}->{$pat_name} = undef;
						};
					}
					unless ($gname) {
						warn "\n\n";
						warn $project->name;
						warn $pat->name;
						warn $var->id;
						warn $var->vector_id;
						foreach my $g (@{$var->getGenes}) {
							warn $g->external_name;
							warn $g->id;
							warn $g->external_name .' -> from '.$g->start.' to '.$g->end;
						}
						warn $v_genes->to_Enum();
						print {$err} $project_name.' '.$pat->name.': '.$var->id.' '.$var->vector_id.":\nGene(s) ".join( ', ', map {$_->name.' '.$_->external_name} @{$var->getGenes}) ." not in gene list\n\n";
						die("Gene(s) ".join( ', ', map {$_->name.' '.$_->external_name} @{$var->getGenes}) ." not in gene list");
						warn ("\n\n");
					}
					
					$infos .= $var->getNbAlleleRef($pat).',';
					$infos .= $var->getNbAlleleAlt($pat).',';
					$infos .= $var->getDepth($pat);
					$nb_var->{$project_name}->{$pat_name}->{'intersection dejavu cov'} ++ if ($var->getDepth($pat) >= 30);
					$link .= "\n" if ($nb_var->{$project_name}->{$pat_name}->{'intersection links'} or $link);
					$link .= 'https://www.polyweb.fr/polyweb/vector/detailProject.html?'
						.'project='.$project_name.'&chromosome='.$chr->name.'&genename='.$gname
						.'&type=ngs&vector_ids='.$var->vector_id.'&type_cache=undefined';
				}
#				warn $infos if ($infos);
				$nb_var->{$project_name}->{$pat_name}->{'intersection details'} .= $infos;
				$nb_var->{$project_name}->{$pat_name}->{'intersection links'} .= $link;
#				warn $nb_var->{$project_name}->{$pat_name}->{'intersection details'};
				$chr->disconnect;
			}
			$nb_var->{$project_name}->{$pat_name}->{'intersection dejavu other projects'} = 0 unless (exists $nb_var->{$project_name}->{$pat_name}->{'intersection dejavu other projects'});
			$nb_var->{$project_name}->{$pat_name}->{'intersection dejavu other patients'} = 0 unless (exists $nb_var->{$project_name}->{$pat_name}->{'intersection dejavu other patients'});
			$nb_var->{$project_name}->{$pat_name}->{'intersection dejavu other patients ho'} = 0 unless (exists $nb_var->{$project_name}->{$pat_name}->{'intersection dejavu other patients ho'});
			$nb_var->{$project_name}->{$pat_name}->{'intersection dejavu other'} = 0 unless (exists $nb_var->{$project_name}->{$pat_name}->{'intersection dejavu other'});
			$nb_var->{$project_name}->{$pat_name}->{'intersection dejavu cov'} = 0 unless (exists $nb_var->{$project_name}->{$pat_name}->{'intersection dejavu cov'});
			
			print {$fh} $project_name.' '.$project->description.';'.$fam->name.';'.$pat_name.';'
				.$nb_var->{$project_name}->{$pat_name}->{'intersection dejavu other'}.';'
				.'"'.$nb_var->{$project_name}->{$pat_name}->{'intersection details'}.'";'
				.'"'.$nb_var->{$project_name}->{$pat_name}->{'intersection links'}.'";'
				.'"'.$nb_var->{$project_name}->{$pat_name}->{'intersection dejavu cov'}.'";';
			print {$fh} "\n";
		}
	}
	foreach my $g (keys %{$count_genes->{'intersection'}->{$project_name}}) {
		my @p = keys %{$count_genes->{'intersection'}->{$project_name}->{$g}};
		print {$stats_genes} $project_name.';'.$g.';'.scalar @p.';'.join(',',@p)."\n";
	}
	
	close ($fh);
	close ($stats_genes);
	close ($err);
	warn Dumper $nb_var->{$project_name};
	warn Dumper $count_genes->{'intersection'}->{$project_name};
	
	$project->disconnect;
	$project = undef;
	$buffer = undef;
	$pm->finish;
}
$pm->wait_all_children();

warn Dumper $nb_var;





open($stats_genes, "<$file_genes") or confess("Can't open '$file_genes'");
foreach my $line (readline($stats_genes)) {
	chomp $line;
	next if ($line =~ /^project;gene/);	# /^project;gene;number of patients;patients$/
	my ($project,$gene,$count,$pats) = split(/;/,$line);
	$count_genes->{'total'}->{$gene} += $count;
	}
close ($stats_genes);

my $file_stats = $file_genes =~ s/\.tmp$//r;
open(my $fh_stats, ">$file_stats") or confess("Can't open '$file_stats'");
print {$fh_stats} "# (CADD > 25 OR spliceAI > 0,5) AND freq < 0,01% AND gnomAD He < 5 AND gnomAD Ho = 0 AND in gene list AND DéjàVu other project < $dv_projects AND DéjàVu other patients < $dv_pat AND DéjàVu other patients ho < $dv_pat_ho\n";
print {$fh_stats} "genes;nb de patients trouvés\n";
foreach my $g (keys %{$count_genes->{'total'}}) {
	print {$fh_stats} $g.';'.$count_genes->{'total'}->{$g}."\n";
}
close ($fh_stats);

warn Dumper $stats_genes;




