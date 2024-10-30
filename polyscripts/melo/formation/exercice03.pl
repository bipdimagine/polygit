#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use Data::Dumper;
use Getopt::Long;
use Carp;

use GBuffer;


my $project_name = 'NGS2015_0794';

my $buffer = new GBuffer;
my $project = $buffer->newProject(-name => $project_name);

my $alignment_method = $project->alignmentMethods;
#print join(', ', @$alignment_method) . "\n";
my $variations_dir = $project->getVariationsDir(@$alignment_method);
#print $variations_dir . "\n";	# /data-isilon/sequencing/ngs/NGS2015_0794/HG19/variations/bwa/

my $project_path = $project->project_path;
#print $project_path . "\n";	# /data-isilon/sequencing/ngs/NGS2015_0794/HG19/
$variations_dir = $project_path."variations/";


open(my $output, '>', 'output_ex3.txt');

sub recurrent_open_dir {
	my ($dir_path) = @_;
	opendir(my $dir, $dir_path);
	while (my $element = readdir($dir)) {
		
		unless ($element eq '.' or $element eq '..') {
			my $path = $dir_path.$element;
			if (-d $path) {
				$path = $path.'/';
				print {$output} "dir\t" . $path . "\n";
				recurrent_open_dir($path);
			}
			else {print {$output} "file\t" . $path . "\n";}
		}
	}
	closedir($dir);
}

recurrent_open_dir($variations_dir);
close ($output);





# 1/ A partir des objets, trouver le dossier appelé "variations" possèdant tous
#    les dossiers des fichiers de calling du projet NGS2015_0794

# 2/ à partir de ce dossier, me répertorier l'ensemble des dossiers / fichiers
#    (y compris ceux DANS les dossiers) et me les lister
#    (recherche Google du type opendir)

# 3/ dans un fichier de sortie, me donner le nom de tous les répertoires et
#    dossiers observés (path complet) ainsi que de me dire si c'est un fichier
#    ou un dossier
#    (recherche Google test if is file / directory)

