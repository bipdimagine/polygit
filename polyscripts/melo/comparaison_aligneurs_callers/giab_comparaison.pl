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
use Cwd 'abs_path';
use Parallel::ForkManager;
use List::Util qw(sum);
use GBuffer;

# Pré-process: restreindre tous les vcf au même bed
# bash /home/mperin/git/polygit/polyscripts/melo/comparaison_aligneurs_callers/restrict_vcf.sh

#perl /home/mperin/git/polygit/polyscripts/melo/giab_comparaison.pl 

my @original_argv = @ARGV;
my $cmdline = "$0 " . join(' ', @original_argv);

my $project_name;
my $test_patient_name;
my $control_patient_name = "HG002";
my $depth = 30;
my $help;
GetOptions(
	'project=s'			=> \$project_name,
	'test_patient=s'	=> \$test_patient_name,
	'control_patient=s'	=> \$control_patient_name,
	'depth=i'			=> \$depth,
	'help'				=> \$help,
);
if ($help) {
	usage();
	exit;
}
unless ($project_name) {
	my $prompt_projects = {"NGS2025_08659     RENOME_V5_Hg19         dragen-align"	=> 'NGS2025_08659',
						   "NGS2025_08660     RENOME_V5_Hg19         bwa"			=> 'NGS2025_08660',
						   "NGS2025_08679     Rapid_V4_TubeF_UMI     dragen-align"	=> 'NGS2025_08679',
						   "NGS2025_08680     Rapid_V4_TubeF_UMI     bwa"			=> 'NGS2025_08680',
						   "NGS2025_08964     Rapid_V4_TubeF         dragen-align"	=> 'NGS2025_08964',
						   };
	$project_name = prompt ("Select a test project:", -menu=>$prompt_projects);
	die unless ($project_name);
}
$depth = 0 if ($depth < 0);
warn "depth = $depth";



my $buffer = new GBuffer;
my $project = $buffer->newProject(-name=>$project_name) || confess ("Can't open project '$project_name': $!");
warn $project->name;

# test
$test_patient_name = 'GIAB_V5' if ($project_name =~ /^NGS2025_086(59|60)$/ and not $test_patient_name);
$test_patient_name = 'GIA_HG0_6123GM002954_STB' if ($project_name =~ /^NGS2025_08(679|680|964)$/ and not $test_patient_name);
die ("--test_patient_name argument required") unless ($test_patient_name);
my $test = $project->getPatient($test_patient_name) || confess ("Can't find patient '$test_patient_name' in project '$project_name': $!");
warn 'Test patient: '.$test->name."  id:".$test->id;
my $methods = $test->alignmentMethods();
confess("More than one alignment method") if scalar(@$methods) > 1;
my $capture = $test->getCapture->name;
warn "Alignment: ".$methods->[0];
warn "Capture: ".$capture;


# GIAB NIST
my $control = $project->getPatient($control_patient_name) || confess ("Can't find patient '$control_patient_name' in project '$project_name': $!");
warn 'Control patient: '.$control->name."  id:".$control->id;


# Méthodes de calling
my $callers = $test->getCallingMethods;
my @callers = sort {$a cmp $b} @$callers;
@callers = grep { $_ !~ /^(melt|duplicate_region_calling)$/ } @callers;
my $nb_callers = scalar @callers;


# Variants trouvés pour patient test par les méthodes de calling
my $all_var = $test->getStructuralVariations;
#$nb_alleles->{'at least 1 caller'}->{'total'} = scalar @$all_var;
warn ("nb var all callers (including indels and low dp) for test patient $test_patient_name: ".scalar @$all_var);
die ("No variant found for test patient '$test_patient_name'. Check the vcf files and calling methods registered") unless (scalar @$all_var);

# Variants GIAB NIST (ref)
my $all_var_giab = $control->getStructuralVariations;
#$nb_alleles->{'nist giab'}->{'total'} = scalar @$all_var_giab;
warn "nb var giab (including indels and low dp) for $control_patient_name: ".scalar @$all_var_giab;
die ("No variant found for '$control_patient_name'. Check the vcf files and calling methods registered") unless (scalar @$all_var_giab);


my @types = ('all variants','snp','indel');
my $nb_alleles;
my $low_dp;
my $lists;
my $Venn;
foreach my $var (@$all_var) {
	my $var_id = $var->id;
	my $annex = $var->annex;
	next unless (exists $annex->{$test->id});
	my $hvar = $annex->{$test->id};
	die unless (exists $hvar->{'method_calling'});
	
	# Filtre les variants peu couverts
	if ($hvar->{'dp'} < $depth) {
		my @methods = keys %{$hvar->{'method_calling'}};
		die unless (scalar @methods);
		my $max_dp = $hvar->{'dp'};
		if (scalar @methods > 1) {
			foreach my $method (@methods) {
				my $dp = $hvar->{'method_calling'}->{$method}->{'nb_all_ref'} + $hvar->{'method_calling'}->{$method}->{'nb_all_mut'} + $hvar->{'method_calling'}->{$method}->{'nb_all_other_mut'};
				$max_dp = $dp if ($dp > $max_dp);
			}
		}
		if ($max_dp < $depth) {
			$low_dp->{'at least 1 caller'} ++;
#			warn Dumper $hvar if $var->isDeletion;
			next;
		}
	}
	
	# Variants trouvés par au moins une méthode de calling
	my @callers_var = keys %{$hvar->{method_calling}};
	my $search = '^('.join('|',@callers).')$';
	my $nb_callers_var = grep ( /$search/, @callers_var );
	if ($nb_callers_var) {
		$nb_alleles->{'snp'}->{'at least 1 caller'}->{'total'} ++ if ($var->isVariation);
		$nb_alleles->{'indel'}->{'at least 1 caller'}->{'total'} ++ if ($var->isInsertion or $var->isDeletion);
		$nb_alleles->{'cnv'}->{'at least 1 caller'}->{'total'} ++ if ($var->isLargeDeletion or $var->isLargeDeletion or $var->isLargeDuplication);
		die ($var_id) unless (($var->isVariation) or ($var->isInsertion or $var->isDeletion));
		die ($var_id) if (($var->isVariation) and ($var->isInsertion or $var->isDeletion));
		$nb_alleles->{'all variants'}->{'at least 1 caller'}->{'total'} ++;
		foreach my $type (@types) {
			my $cond;
			$cond = 1 if ($type eq 'all variants');
			$cond = $var->isVariation if ($type eq 'snp');
			$cond = ($var->isInsertion + $var->isDeletion) if ($type eq 'indel');
			if ($cond) {
				$nb_alleles->{$type}->{'at least 1 caller'}->{'TP'} ++ if (exists $annex->{$control->id});
				$nb_alleles->{$type}->{'at least 1 caller'}->{'FP'} ++ unless (exists $annex->{$control->id});
			}
		}
	}
	
	# Variants trouvés par chacune des méthodes de calling
	foreach my $caller (@callers_var) {
		$nb_alleles->{'all variants'}->{$caller}->{'total'} ++;
		$nb_alleles->{'snp'}->{$caller}->{'total'} ++ if ($var->isVariation);
		$nb_alleles->{'indel'}->{$caller}->{'total'} ++ if ($var->isInsertion or $var->isDeletion);
		$nb_alleles->{'cnv'}->{$caller}->{'total'} ++ if ($var->isLargeDeletion or $var->isLargeInsertion or $var->isLargeDuplication);
		die ($var_id) unless (($var->isVariation) or ($var->isInsertion or $var->isDeletion));
		die ($var_id) if (($var->isVariation) and ($var->isInsertion or $var->isDeletion));
		foreach my $type (@types) {
			my $cond;
			$cond = 1 if ($type eq 'all variants');
			$cond = $var->isVariation if ($type eq 'snp');
			$cond = ($var->isInsertion + $var->isDeletion) if ($type eq 'indel');
			if ($cond) {
				if (exists $annex->{$control->id}) {
					$nb_alleles->{$type}->{$caller}->{'TP'} ++;
					push (@{$nb_alleles->{$type}->{$caller}->{'dp'}->{'TP'}}, $hvar->{'dp'});
					push(@{$lists->{$type}->{'TP'}->{$caller}}, $var->id);
				}
				unless (exists $annex->{$control->id}) {
					$nb_alleles->{$type}->{$caller}->{'FP'} ++;
					push(@{$lists->{$type}->{'FP'}->{$caller}}, $var->id);
					push (@{$nb_alleles->{$type}->{$caller}->{'dp'}->{'FP'}}, $hvar->{'dp'});
				}
			}
		}
	}
	
	# Variants trouvés par toutes les méthodes de calling
	if ($nb_callers_var == $nb_callers) {
		$nb_alleles->{'all variants'}->{"all $nb_callers callers"}->{'total'} ++;
		$nb_alleles->{'snp'}->{"all $nb_callers callers"}->{'total'} ++ if ($var->isVariation);
		$nb_alleles->{'indel'}->{"all $nb_callers callers"}->{'total'} ++ if ($var->isInsertion or $var->isDeletion);
		$nb_alleles->{'cnv'}->{"all $nb_callers callers"} ++ if ($var->isLargeDeletion or $var->isLargeInsertion or $var->isLargeDuplication);
		foreach my $type (@types) {
			my $cond;
			$cond = 1 if ($type eq 'all variants');
			$cond = $var->isVariation if ($type eq 'snp');
			$cond = ($var->isInsertion + $var->isDeletion) if ($type eq 'indel');
			if ($cond) {
				$nb_alleles->{$type}->{"all $nb_callers callers"}->{'TP'} ++ if (exists $annex->{$control->id});
				$nb_alleles->{$type}->{"all $nb_callers callers"}->{'FP'} ++ unless (exists $annex->{$control->id});
			}
		}
	}
	
#	foreach my $i (0 .. $#callers) {
#		$Venn->{'TP'}->{$var_id} |= (1 << $i) if ($var->isVariation && exists $hvar->{'method_calling'}->{$callers[$i]} && exists $annex->{$control->id});
#		$Venn->{'FP'}->{$var_id} |= (1 << $i) if ($var->isVariation && exists $hvar->{'method_calling'}->{$callers[$i]} && not exists $annex->{$control->id});
#	}
}

# Variants du GIAB NIST (ref) & FN
foreach my $var (@$all_var_giab) {
	my $var_id = $var->id;
	my $annex = $var->annex;
	my $hvar = $annex->{$control->id};
	
	if ($hvar->{'dp'} < $depth) {
		$low_dp->{'nist giab'} ++;
	}
	
	$nb_alleles->{'snp'}->{'nist giab'}->{'total'} ++ if ($var->isVariation);
	$nb_alleles->{'indel'}->{'nist giab'}->{'total'} ++ if ($var->isInsertion or $var->isDeletion);
	$nb_alleles->{'cnv'}->{'nist giab'}->{'total'} ++ if ($var->isLargeDeletion or $var->isLargeDeletion or $var->isLargeDuplication);
	$nb_alleles->{'all variants'}->{'nist giab'}->{'total'} ++;
	
	foreach my $type (@types) {
		my $cond;
		$cond = 1 if ($type eq 'all variants');
		$cond = $var->isVariation if ($type eq 'snp');
		$cond = ($var->isInsertion + $var->isDeletion) if ($type eq 'indel');
		next unless($cond);
#		warn $var_id;
		my @callers_var = keys %{$var->annex->{$test->id}->{'method_calling'}};
		my $search = '^('.join('|',@callers).')$';
		my $nb_callers_var = scalar grep ( /$search/, @callers_var);
		unless ($nb_callers_var) {
			$nb_alleles->{$type}->{"all $nb_callers callers"}->{'FN'} ++ ;
#			warn "FN 'all $nb_callers callers'";
#			warn  Dumper $hvar
		}
		$nb_alleles->{$type}->{"at least 1 caller"}->{'FN'} ++ unless ($nb_callers_var == $nb_callers);
#		warn "FN 'at least 1 caller'" unless ($nb_callers_var == $nb_callers);
			
		foreach my $caller (@callers) {
			unless (exists $annex->{$test->id}->{'method_calling'}->{$caller}) {
				$nb_alleles->{$type}->{$caller}->{'FN'} ++ ;
				push (@{$nb_alleles->{$type}->{$caller}->{'dp'}->{'FN'}}, $hvar->{'dp'});
				push(@{$lists->{$type}->{'FN'}->{$caller}}, $var->id);
#				warn "FN '$caller'";
			}
		}
#		warn Dumper $annex if (exists $annex->{$test->id} && exists $annex->{$test->id}->{'method_calling'});
#		warn $var->type if (exists $annex->{$test->id} && exists $annex->{$test->id}->{'method_calling'});
#		warn defined $annex->{$test->id}->{'method_calling'} if (exists $annex->{$test->id} && exists $annex->{$test->id}->{'method_calling'});
#		warn exists $annex->{$test->id}->{'method_calling'} if (exists $annex->{$test->id} && exists $annex->{$test->id}->{'method_calling'});
#		die if (exists $annex->{$test->id} && exists $annex->{$test->id}->{'method_calling'});
	
#		foreach my $i (0 .. $#callers) {
#			unless (exists $annex->{$test->id} && exists $annex->{$test->id}->{'method_calling'} && exists $annex->{$test->id}->{'method_calling'}->{$callers[$i]}) {
#				$Venn->{$type}->{'TP_giab'}->{$var_id} |= (1 << $i) if (exists $annex->{$test->id}->{'method_calling'}->{$callers[$i]} );
#				$Venn->{$type}->{'FN'}->{$var_id} |= (1 << $i) unless (exists $annex->{$test->id}->{'method_calling'}->{$callers[$i]} );
#			}
#		}
	}
}

#warn "\nnb variants < ".$depth."X: \n".Dumper($low_dp)."\n";
warn "nb variants < ".$depth."X excluded: \n\t$test_patient_name: ".$low_dp->{'at least 1 caller'};#."\n\t$control_patient_name: ".$low_dp->{'nist giab'}."\n";
warn Dumper $low_dp;



# Ecrit les résultats dans des fichiers tsv/list

# Stats par caller
my $dir = $project_name.'_'.$capture.'_'.$methods->[0].'_'.$depth.'X/';
system ("mkdir $dir") unless (-d $dir);
my $file_count = $dir.'stats_by_caller.tsv';
warn $file_count;
open (my $out, '>', $file_count) || confess ("Can't open '$file_count': $!");
print {$out} "##cmdline: $cmdline\n";
print {$out} "##$project_name\n";
print {$out} "##capture: ".$capture ."\n";
print {$out} "##alignment method: ".$methods->[0] ."\n";
print {$out} "##depth >= $depth\n";
my @columns = ('SNP','indel','total',"SNP\nTP","SNP\nFP","SNP\nFN","SNP\nAccuracy","SNP\nSensitivity","indel\nTP","indel\nFP","indel\nFN","indel\nAccuracy","indel\nSensitivity");
print {$out} '#caller'."\t\"".join("\"\t\"", @columns) ."\"\n";
my @lines = @callers;
push(@lines, 'at least 1 caller', "all $nb_callers callers", 'nist giab');
foreach my $caller (@lines) {
	print {$out} $caller;
	foreach my $type (@types) {
		foreach my $group ('TP', 'FP', 'FN') {
			if (exists $nb_alleles->{$type}->{$caller}->{'dp'} && exists $nb_alleles->{$type}->{$caller}->{'dp'}->{$group}) {
				$nb_alleles->{$type}->{$caller}->{'dp'}->{$group} = sum(@{$nb_alleles->{$type}->{$caller}->{'dp'}->{$group}}) / scalar @{$nb_alleles->{$type}->{$caller}->{'dp'}->{$group}};
			}
		}
		unless ($caller eq 'nist giab' or !$nb_alleles->{$type}->{$caller}->{'TP'} ) {
			$nb_alleles->{$type}->{$caller}->{'Accuracy'} = $nb_alleles->{$type}->{$caller}->{'TP'} /( $nb_alleles->{$type}->{$caller}->{'TP'} + $nb_alleles->{$type}->{$caller}->{'FP'} );
			$nb_alleles->{$type}->{$caller}->{'Sensitivity'} = $nb_alleles->{$type}->{$caller}->{'TP'} /( $nb_alleles->{$type}->{$caller}->{'TP'} + $nb_alleles->{$type}->{$caller}->{'FN'} );
		}
	}
	foreach my $col ('snp','indel','all variants') {
		print {$out} "\t".$nb_alleles->{$col}->{$caller}->{'total'};
	}
	foreach my $type ('snp','indel') {
		foreach my $col ('TP','FP','FN','Accuracy','Sensitivity') {
			print {$out} "\t".$nb_alleles->{$type}->{$caller}->{$col};
		}
	}
	print {$out} "\n";
}
close($out);



# Listes des variants par catégorie et par caller
foreach my $cat ('TP','FP','FN') {
	foreach my $type ('snp','indel') {
	system("mkdir $dir$type") unless (-d $dir.$type);
	system("mkdir $dir$type/$cat") unless (-d $dir.$type.'/'.$cat);
		foreach my $caller (@callers) {
			my $file_liste_var = $dir.$type.'/'.$cat.'/variants_'.$caller.'.list';
	#		warn $file_liste_var;
			open (my $fh, '>', $file_liste_var);
			print {$fh} join("\n", @{$lists->{$type}->{$cat}->{$caller}})."\n" if (exists $lists->{$type}->{$cat}->{$caller});
			close($fh);
		}
	}
}

# Table comptage pour diagramme de Venn
sub venn {
	my ($hvariants, $legend) = shift @_;
	my $venn_table;
	foreach my $variant (keys %$hvariants) {
		my $presence_mask = $hvariants->{$variant};
		$venn_table->{$presence_mask}++;
	}
	if ($legend) {	# Affichage de la légende des masques
#		print "\nLegend (Mask -> Callers):\n";
#		foreach my $mask (sort { $a <=> $b } keys %$venn_table) {
#		my @active_callers;
#			for my $i (0 .. $#callers5) {
#				push @active_callers, $callers5[$i] if ($mask & (1 << $i));
#			}
#			printf "%05b -> %s\n", $mask, join(", ", @active_callers);
#		}
#		print "\n";
		# Légende générale
		print "#Legend (Mask -> Callers):\n";
		for my $i (0 .. $#callers) {
			my $mask = 1 << $i;
			printf "#%05b -> %s\n", $mask, $callers[$i];
		}
	}
	print "#Combination\tCount\n";
	foreach my $mask (sort { $a <=> $b } keys %$venn_table) {
		printf "%05b\t%d\n", $mask, $venn_table->{$mask};
	}
	return $venn_table;
}

#foreach my $cat (keys %$Venn) {
#	system("mkdir $dir$cat") unless (-d $dir.$cat);
#	my $file_venn = $dir.$cat.'/venn_tables.tsv';
#	open(my $fh, '>', $file_venn) || die ("Can't open file '$file_venn': $!");
#	select $fh; # all print will go to $fh by default
#	venn($Venn->{$cat}, 0);
#	close($fh);
#}
#select STDOUT;






warn Dumper $nb_alleles;
print "\n\n";

#warn Dumper $lists->{'indel'}->{'FN'};










sub usage {
	print
"
giab_comparaison.pl
-------------
Optionels:
	project <s>		nom du projet, défaut: liste à choix
	test_patient <s>	nom du patient/échantillon à comparer
	control_patient <s>		nom de patient/échantillon de control/référence auquel comparer test_patient [HG002]
	depth <i>		ne garde que les variants dont la profondeur de couverture est supérieure ou égale à cette valeur, défaut: 30
	help			affiche ce message

";
}



# todo: récupérer les listes des tous les variants id du projet par caller (5 + giab) -> 6 listes -> diagramme de Venn à 6 ensembles ?








#warn "\n\n";
#warn $test->getProject->name();
#
#warn ref ($test).' -> '.$test->name();
#warn "\n\n";
#foreach my $v (@{$test->getStructuralVariations}) {
#	if ($v->name eq '16-15892492-G-A') {
#		warn "\n\n";
#		warn ref ($v).' -> '.$v->id;
#		warn $v->name;
#	}
#	if ($v->name eq '16-15892492-C-A') {
#		warn "\n\n";
#		warn ref ($v).' -> '.$v->id;
#		warn $v->name;
#	}
#}

#tabix /data-isilon/sequencing/ngs/NGS2025_08659/HG19_MT/variations/freebayes/GIAB_V5.vcf.gz chr16:15892480-15892499
#chr16	15892481	.	C	A	6872.42	.	AB=0;ABP=0;AC=2;AF=1;AN=2;AO=229;CIGAR=1X;DP=229;DPB=229;DPRA=0;EPP=320.567;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=59.8166;MQMR=0;NS=1;NUMALT=1;ODDS=322.067;PAIRED=0.982533;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=8102;QR=0;RO=0;RPL=108;RPP=4.61283;RPPR=0;RPR=121;RUN=1;SAF=131;SAP=13.3366;SAR=98;SRF=0;SRP=0;SRR=0;TYPE=snp	GT:DP:AD:RO:QR:AO:QA:GL	1/1:229:0,229:0:0:229:8102:-728.558,-68.9359,0
#chr16	15892492	.	CAGAAG	AAAAAA	95.8026	.	AB=0.230769;ABP=133.967;AC=1;AF=0.5;AN=2;AO=48;CIGAR=1X1M1X2M1X;DP=208;DPB=213;DPRA=0;EPP=98.736;EPPR=114.54;GTI=0;LEN=6;MEANALT=12;MQM=59.125;MQMR=60;NS=1;NUMALT=1;ODDS=22.0594;PAIRED=0.979167;PAIREDR=0.979167;PAO=30;PQA=292;PQR=0;PRO=0;QA=1159;QR=4447;RO=144;RPL=47;RPP=98.736;RPPR=220.158;RPR=1;RUN=1;SAF=48;SAP=107.241;SAR=0;SRF=23;SRP=147.835;SRR=121;TYPE=complexGT:DP:AD:RO:QR:AO:QA:GL	0/1:208:144,48:144:4447:48:1159:-63.6519,0,-333.341
#[mperin@sanger comparaison_aligneurs_callers]$ perl giab_comparaison.pl -project=NGS2025_08659
#depth = 30 at giab_comparaison.pl line 62.
#NGS2025_08659 at giab_comparaison.pl line 67.
#GIAB_V5 id:147310 at giab_comparaison.pl line 74.
#Control id:147889 at giab_comparaison.pl line 78.
#
#
#NGS2025_08659 at giab_comparaison.pl line 92.
#GenBoPatient -> GIAB_V5 at giab_comparaison.pl line 94.
#
#
#
#
#GenBoVariation -> 16_15892492_C_A at giab_comparaison.pl line 104.
#16-15892492-C-A at giab_comparaison.pl line 105.
#
#
#GenBoVariation -> 16_15892494_G_A at giab_comparaison.pl line 99.
#16-15892492-G-A at giab_comparaison.pl line 100.
#
#
#GenBoVariation -> 16_15892497_G_A at giab_comparaison.pl line 99.
#16-15892492-G-A at giab_comparaison.pl line 100.
#Died at giab_comparaison.pl line 120.
#[mperin@sanger comparaison_aligneurs_callers]$ 
