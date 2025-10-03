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
use List::MoreUtils qw(natatime);
use GBuffer;
use Time::HiRes qw(gettimeofday tv_interval);


# Elodie BAL

# cat medium impact -> medium:maturemirna,stop,phase,frameshift,essential_splicing,ncrna,nonsynonymous,nonframeshift,splice_site,predicted_splice_site
# (gnomad) allelecount 10, ho 5
# dejavu 10, ho 5
# cov 15
# seulement les transcripts principaux (MANE,CDDS)


# J'ai regardé des listes que j'avais sortie de Polyquiery par genes et j'obtenais entre 300 et 400 par patients avec les filtres suivants:
# FILTER(S) TYPE VARIATION
# Synonymous
# Utr
# Intergenic
# Up/Downstream
# Intronic
# ncRNA
# Pseudogene
# freq_001
# freq_01
# freq_05
# freq_1
# J'excluais en effet les ncRNA (difficile d'interprétation sur une première approche)
# et je mettais un filtre de fréquence de 10-4 sur gnomad. Je préfère en fréquence qu'en nombre de fois vu car tout dépend du dénominateur. 
# Et j'avais mis un déjà vu à <=10 aussi. Mais j'aimerai voir si le filtre gnomad seul ne suffit pas avec tous les autre critères de couverture.
# J'avais remarqué que dans cette liste il y avait des faux positifs et je souhaiterai donc rajouter un filtre de coverage : >=15 me semble bien pour commencer
# et aussi le nombre de reads avec mutation pour éviter les faux positifs. je pensais à 30%.



my $out_dir = '/home/mperin/git/polygit/polyscripts/melo/vecteurs_variants/Elodie_Bal/';
my $project_name = 'NGS2022_5039';
my $fork = 1;
my $ac = 10;
my $ac_ho = 5;
my $dv = 10;
my $dv_ho = 5;
my $cov = 15;
my $prop_alt = 0.3;	# proportion d'alleles alt / depth

my $dv_projects = 1;
my $dv_pat = 5;
my $dv_pat_ho = 0;


GetOptions(
	'project=s'				=> \$project_name,
	'out_dir=s'				=> \$out_dir,
	'fork=i'				=> \$fork,
) || die ("Error in command line arguments\n");

die ("fork should be [1;40], given $fork") if ($fork > 40 or $fork <= 0);
#$fork = 5 if ($fork > 5);
#warn ("fork decreased to 5 for memory usage") if ($fork > 5);
warn 'fork='.$fork;
die ("$out_dir doesn't exists") unless (-d $out_dir);


my $nb_var;
my $filters;


my $buffer = new GBuffer;
my $project = $buffer->newProjectCache(-name=>$project_name) or confess ("Can\'t open project cache '$project_name'");
warn $project->name;


my $file_nb_var = $out_dir."nb_variants.csv";
warn $file_nb_var;
open(my $fh, ">$file_nb_var") or confess("Can't open '$file_nb_var': $!");
print {$fh} "project;family;affected child;";
print {$fh} "(gnomad ac he <= $ac or ho <= $ac_ho) & freq <= 10^-4 & medium impact & cov >= $cov & prop all alt >= $prop_alt & main transcripts;";
print {$fh} "& (dejavu he <= $dv or ho <= $dv_ho);";
#print {$fh} "Comment";
print {$fh} "\n";
close ($fh);

my $file_var = $out_dir."variants.csv";
warn $file_var;
open(my $fh, ">$file_var") or confess("Can't open '$file_var': $!");
print {$fh} "project name;project description;family;affected child;variant id;db freq;type;GnomAD HO;GnomAD AC;dejavu;dejavu ho;nb all alt;depth;gene(s);transcript(s);link";
print {$fh} "\n";
close ($fh);



my @list_children;
my $families = $project->getFamilies;
warn scalar @$families.' families';# & '.scalar @{$project->getPatients()}.' patients';
my @families = sort {$a->name cmp $b->name} @$families;
foreach my $fam (@families) {
	my @children = keys(%{$fam->children_ill});
	warn ("no ill child found: $project_name, ".$fam->name) unless (scalar @children);
	next unless (scalar @children);
	push (@list_children, @children);
}
warn Dumper sort @list_children;
	
$project->getPatients();
my $pm = new Parallel::ForkManager($fork);

$buffer->disconnect();

foreach my $child (sort @list_children) {
	my $pid = $pm->start and next;
	$buffer->dbh_reconnect();
	
	my $pat = $project->getPatient($child);
	my $pat_name = $pat->name;
	my $fam = $pat->getFamily;
	warn $pat_name;
	$nb_var->{$project_name}->{$pat_name}->{'commentaire'} = 'Missing both parents' unless ($fam->mother || $fam->father);
	$nb_var->{$project_name}->{$pat_name}->{'commentaire'} = 'Missing mother' if ($fam->father && not $fam->mother);
	$nb_var->{$project_name}->{$pat_name}->{'commentaire'} = 'Missing father' if ($fam->mother && not $fam->father);
	warn $nb_var->{$project_name}->{$pat_name}->{'commentaire'} if (exists $nb_var->{$project_name}->{$pat_name}->{'commentaire'});
	
	open(my $fh, ">>$file_var") or confess("Can't open '$file_var': $!");
	foreach my $chr (@{$project->getChromosomes}) {
		
		warn "chr ".$chr->name;
		next unless ($chr->getNewVector()->Size);
		my $v0 = $pat->getVectorVariants($chr)->Clone();
#		warn 'total variants: '.$pat->getVariantsVector($chr)->Norm;
		my $v_intersection = $v0->Clone();
		
		## GnomAD
		# todo: en freq
		my $v_gnomad_ac_he = $v0 & $chr->getVectorScore("gnomad_ac_$ac");
		my $v_gnomad_ac_ho = $v0 & $chr->getVectorScore("gnomad_ho_ac_$ac_ho");
		$v_intersection = $v_gnomad_ac_he + $v_gnomad_ac_ho;


		## Frequence
		my $v_freq = $chr->getNewVector()->Clone();
#		warn Dumper $chr->vector_global_categories;
		foreach my $filter_name ('freq_0001', 'freq_none') {	# prendre la freq en question + tous celles en dessous
			$v_freq += $chr->vector_global_categories->{$filter_name};
		}
		$v_freq &= $v0;
		$v_intersection &= $v_freq;


		## Impact factor
		my $v_medium_impact = $chr->getNewVector()->Clone();
		foreach my $gene (@{$chr->getGenes()}) {
			$chr->getProject->print_dot(1);
#			foreach my $filter (keys %{$gene->categories_intspan()}) {
#				$filters->{$filter} += 1;
#			}
				my $var_gene_annot = $gene->getCategoriesVariantsVector({
				# medium impact (not the good one) -> maturemirna,stop,phase,frameshift,essential_splicing,ncrna,nonsynonymous,nonframeshift,splice_site,predicted_splice_site
				# high impact -> Mature miRNA;Splice Acc/Don;Frameshift;Stop-gained;(Start/Stop)-lost = maturemirna,essential_splicing,frameshift,stop,phase
				# medium impact -> ncRNA;Splice Region;Missense;No-frameshift = ncrna,splice_site,nonsynonymous,coding,nonframeshift
				# catégories à éliminer
#				'essential_splicing'    => 0,
#				'frameshift'            => 0,
#				'maturemirna'           => 0,
				'ncrna'                 => 1,
#				'non-frameshift'        => 0,
#				'nonframeshift'         => 0,
#				'nonsynonymous'         => 0,
#				'phase'                 => 0,
#				'predicted_splice_site' => 0,
#				'splice_site'           => 0,
#				'splicing'              => 0,
#				'stop'                  => 0,
#				'coding'                => 0,
				'downstream'            => 1,
				'intergenic'            => 1,
				'intronic'              => 1,
				'pseudogene'            => 1,
				'silent'                => 1,
				'upstream'              => 1,
				'upstream_downstream'   => 1,
				'utr'                   => 1,
			});
			$v_medium_impact += $var_gene_annot;
		}
		$v_medium_impact &= $v0;
		$v_intersection &= $v_medium_impact;
		
		
		## Couverture et transcripts
		my $v_cov = $chr->getNewVector()->Clone();
		my $v_main_transcript = $chr->getNewVector()->Clone();
		foreach my $var (@{$chr->getListVarObjects($v_intersection)}) {
			my $vid = $var->vector_id;
			next unless ($var->getDepth($pat) >= $cov);
			next unless ($var->getNbAlleleAlt($pat)/$var->getDepth($pat) >= $prop_alt);
			$v_cov->Bit_On($vid);
			my $transcripts = $var->getTranscripts;
			foreach my $transcript (@$transcripts) {
				next unless ($transcript->isMain);
				$v_main_transcript->Bit_On($vid);
			}
		}
		$v_intersection &= $v_cov;
		$v_intersection &= $v_main_transcript;
		
		
		## DéjàVu
		my $v_dejavu_he = $chr->get_vector_dejavu($v_intersection, $dv, 0);
		my $v_dejavu_ho = $chr->get_vector_dejavu($v_intersection, $dv_ho, 1);
		my $v_dejavu = $v_dejavu_he + $v_dejavu_ho;
		my $v_intersection_dv = $v_intersection & $v_dejavu;
		
		
		$nb_var->{$project_name}->{$pat_name}->{'intersection'} += $v_intersection->Norm;
		$nb_var->{$project_name}->{$pat_name}->{'intersection & dv'} += $v_intersection_dv->Norm;
		
		
		my $link;
		foreach my $var (@{$chr->getListVarObjects($v_intersection)}) {
			warn $var->id;
			my $infos = $project_name.';'.$project->description.';'.$fam->name.';'.$pat_name.';';
			$infos .= $var->id.';';
			$infos .= $var->frequency.';';
			$infos .= $var->variationTypeInterface.';';
			$infos .= $var->getGnomadHO().';';
			$infos .= $var->getGnomadAC().';';
			$infos .= $var->nb_dejavu.';';
			$infos .= $var->nb_dejavu_ho.';';
			$infos .= $var->getNbAlleleAlt($pat).';';
			$infos .= $var->getDepth($pat).';';
			my $transcripts = $var->getTranscripts;
			my @gene_names;
			my @gene_id;
			foreach my $transcript (@$transcripts) {
				next unless ($transcript->isMain);
				my @gnames = map {$_->external_name} @{$transcript->getGenes()};
				my @gid = map {$_->id} @{$transcript->getGenes()};
				warn Dumper \@gnames;
				push (@gene_names, @gnames);
				push (@gene_id, @gid);
			}
#			warn join(', ', @gene_names ).';';
			my @transcript_names = map {$_->name} @$transcripts;
			$infos .= join(', ', @gene_names ).';';
			$infos .= join(', ', @transcript_names ).';';
			
			$link = 'https://www.polyweb.fr/polyweb/vector/detailProject.html?'
				.'project='.$project_name.'&chromosome='.$chr->name#.'&genename='.$var->getGenes->[0]->name
				.'&type=ngs&vector_ids='.$var->vector_id.'&type_cache=undefined';
			$infos .= $link;
			print {$fh} $infos."\n";
		}
		
		$chr->disconnect;
	}
	$nb_var->{$project_name}->{$pat_name}->{'intersection'} = 0 unless (exists $nb_var->{$project_name}->{$pat_name}->{'intersection'});

	close($fh);
	
	open(my $fh, ">>$file_nb_var") or confess("Can't open '$file_nb_var'");
	print {$fh} $project_name.';'.$fam->name.';'.$pat_name.';'
		.'"'.$nb_var->{$project_name}->{$pat_name}->{'intersection'}.'";'
		.'"'.$nb_var->{$project_name}->{$pat_name}->{'intersection & dv'}.'";';
#		.$nb_var->{$project_name}->{$pat_name}->{'commentaire'};
	print {$fh} "\n";
	close ($fh);
	$pm->finish;
}
$pm->wait_all_children();



#$project->disconnect;
#$project = undef;
#$buffer = undef;

#warn Dumper $nb_var;


























