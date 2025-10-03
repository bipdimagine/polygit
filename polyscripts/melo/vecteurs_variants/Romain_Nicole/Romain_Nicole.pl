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

#Je te contacte au sujet de la cohorte de patients présentant un trouble spécifique du langage oral (TSLO) que vous aviez étudiés avec Clothilde, Valérie et Vincent (https://pubmed.ncbi.nlm.nih.gov/39948625/) après un séquençage de génome à Imagine.
# 
#Dans le cadre de ma thèse de sciences, il y a une analyse complémentaire que j'aimerais réaliser sur ces patients.
#Pour toute les familles, à l'exception de la famille DLD-12 (qui a déjà un diagnostic moléculaire), je voulais savoir s'il te serait possible d'extraire la liste de tous les gènes dans lesquels il y a des variants perte-de-fonction présents à l'état hétérozygote ?
#La nomenclature spécifique du variant ne m'intéresse pas vraiment (nonsense, frameshift, splicing etc.) tant que le variant est tronquant.
# 
#Dans l'idéal j'aimerais avoir pour chaque individu de cette cohorte:
#statut: atteint/malade
#génération: parent/enfant
#liste de tous les gènes porteurs de variants tronquants à l'état hétérozygote. La nomenclature du variant n'est pas indispensable mais si cela est extrêmement facilement récupérable sans temps de calcul supplémentaire, je ne dis pas non ^^
#La cerise sur le gâteau ça serait même d'avoir l'identifiant d'ensembl.org pour chaque gène. Dans mon dataset pour ma thèse, j'utilise le ENSG...... pour identifier chaque gène. Au pire je peux faire un mapping ENSG_id - gene_name si besoin, donc ne t'embête pas si ce n'est pas évident. 
# 
#Dans l'idée j'aimerais pouvoir proposer des paires de variants tronquants dans des paires de gènes qui seraient responsables de pathologies humaines, alors que chacun de ces variants serait individuellement bénin. J'ai calculé un score équivalent à une sorte de "pLI pour des paires de gènes" et j'aimerais l'appliquer à des données de patients. Ca serait une preuve très forte de la pertinence de mon score si j'arrivais à montrer des paires de gènes candidats dans des patients actuellement sans diagnostic.



my $out_dir = '/home/mperin/git/polygit/polyscripts/melo/vecteurs_variants/Romain_Nicole/';
my $project_name;
my $fork = 1;
# NGS2025_08802 # cache rocks
# NGS2023_7068 NGS2023_7051 NGS2023_6816 NGS2022_5624 NGS2022_5335 NGS2022_5247 NGS2021_4769 NGS2021_4030 NGS2021_3758

my @args = @ARGV;
GetOptions(
	'project=s'				=> \$project_name,
	'out_dir=s'				=> \$out_dir,
	'fork=i'				=> \$fork,
) || die ("Error in command line arguments\n");

die ("No project specified") unless ($project_name);
die ("fork should be [1;40], given $fork") if ($fork > 40 or $fork <= 0);
warn 'fork='.$fork;
die ("$out_dir doesn't exists") unless (-d $out_dir);

warn Dumper @args;
if ($project_name eq 'NGS2025_08802') {
	my $cmd_hg38 = 'perl /home/mperin/git/polygit-hg38/polyscripts/melo/vecteurs_variants/Romain_Nicole/Romain_Nicole.pl '.join(' ', @args);
	warn $cmd_hg38;
	my $exit = system($cmd_hg38);
	exit($exit)
}


my $buffer = new GBuffer;
my $project = $buffer->newProjectCache(-name=>$project_name) or confess ("Can\'t open project cache '$project_name'");
warn $project->name;



my $families = $project->getFamilies;
warn scalar @$families.' families';# & '.scalar @{$project->getPatients()}.' patients';
my @families = sort {$a->name cmp $b->name} @$families;
$buffer->disconnect();
my $pm = new Parallel::ForkManager($fork);
foreach my $fam (@families) {
	my $pid = $pm->start and next unless ($fork == 1);
	$buffer->dbh_reconnect();
	my $fname = $fam->name;
	warn $fname;
	my $file_out = $out_dir."$project_name-$fname.csv";
	warn $file_out;
	open(my $fh, ">$file_out") or confess("Can't open '$file_out': $!");
	my $children = $fam->getChildren();
	my @children = sort {$a->name cmp $b->name} @$children;
	warn scalar @children. 'children';
	warn (join(', ', map {$_->name} @children));
#	my $child = $fam->getChild;
	my $mother = $fam->getMother();
	my $father = $fam->getFather();
	print {$fh} "gene;variant;transmission;type;".join(';', map {$_->name} @children).';'.$mother->name.';'.$father->name."\n";
	
	my $fpat = $fam->getPatients;
	warn join(', ', map {$_->name} @$fpat);
	
	foreach my $chr (@{$project->getChromosomes}) {
		warn "chr ".$chr->name;
		warn $chr->getNewVector()->Size;
		next unless ($chr->getNewVector()->Size);
		
		my $v_fam = $children->[0]->getVectorVariants($chr)->Clone() 
					+ $mother->getVectorVariants($chr)->Clone() 
					+ $father->getVectorVariants($chr)->Clone();
		$v_fam += $_->getVectorVariants($chr)->Clone() for @$children;
		warn scalar @{$chr->getGenes()};
		foreach my $gene (@{$chr->getGenes()}) {
			$chr->getProject->print_dot(1);
			my $var_gene_annot = $gene->getCategoriesVariantsVector({
				# catégories à éliminer
				'essential_splicing'    => 1,
#				'frameshift'            => 0,
				'maturemirna'           => 1,
				'ncrna'                 => 1,
				'non-frameshift'        => 1,
				'nonframeshift'         => 1,
				'nonsynonymous'         => 1,
				'phase'                 => 1, #
				'predicted_splice_site' => 1,
				'splice_site'           => 1,
				'splicing'              => 1,
				'essential_splicing'	=> 1,
#				'stop'                  => 0,
				'coding'                => 1,
				'downstream'            => 1,
				'intergenic'            => 1,
				'intronic'              => 1,
				'pseudogene'            => 1,
				'silent'                => 1,
				'upstream'              => 1,
				'upstream_downstream'   => 1,
				'utr'                   => 1,
			});

##			warn Dumper keys %{$gene->compact_vector};
#			my $var_gene_annot = $gene->getVectorOriginCategory('frameshift')
#								+ $gene->getVectorOriginCategory('stop');
##								+ $gene->getVectorOriginCategory('phase');
			$var_gene_annot &= $v_fam;
			next unless ($var_gene_annot->Norm);
			foreach my $var (@{$chr->getListVarObjects($var_gene_annot)}) {
				warn $var->id;
				my $infos = $gene->name.';'.$var->id.';';
				$infos .= join(',', map {$var->getTransmissionModelType($fam,$_,$gene)} @children).';';
#				$infos .= $var->variationType($gene).';';
				$infos .= $var->variationTypeInterface($gene).';';
				$infos .= join(';', map {$var->getNbAlleleAlt($_).'/'.$var->getDepth($_)} @children).';';
#				warn '$child->depth: '.Dumper $child->depth($chr->name,$var->start,$var->end);
#				warn '$var->getDepth: '.$var->getDepth($child);	# prend la valeur DP dans vcf
#				if ($var->getDepth($child) != $child->meanDepth($chr->name,$var->start,$var->end)) {
#					warn '$child->depth +/-10: '.Dumper $child->depth($chr->name,$var->start-10,$var->end+10);
#					warn '$child->meanDepth: '.$child->meanDepth($chr->name,$var->start,$var->end);
#				}
#				warn "\n\n";
				$infos .= $var->getNbAlleleAlt($mother).'/'.$var->getDepth($mother).';';
				$infos .= $var->getNbAlleleAlt($father).'/'.$var->getDepth($father);
				print {$fh} $infos."\n";
			}
		}
	}
	close ($fh);
	$pm->finish;
}
$pm->wait_all_children();


	

		
















