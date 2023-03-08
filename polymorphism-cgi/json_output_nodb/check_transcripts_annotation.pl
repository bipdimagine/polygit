#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use Data::Dumper;
use JSON;

use lib "$Bin/../GenBo/lib/obj-nodb/";
use GBuffer;
use GenBoProject;


my $cgi = new CGI();
my $annotation	= $cgi->param('annotation');
my $file_name	= $cgi->param('file');
my $transcripts	= $cgi->param('tr');
my $genes	= $cgi->param('genes');

my $buffer = GBuffer->new();
my $project = $buffer->newProject(-name=>'NGS2015_0794');
my ($genecode_version, $annot_version);
unless ($annotation) {
	$genecode_version = $buffer->getQuery()->getMaxGencodeVersion();
	$annot_version = $buffer->getQuery()->getMaxPublicDatabaseVersion();
	$annotation = $genecode_version.'.'.$annot_version;
}

$project->changeAnnotationVersion($annotation, 1);
warn 'release: '.$project->annotation_version();

my (@lTranscripts, @lGenes);
if ($file_name) {
	open (FILE, $file_name);
	while (<FILE>){
		chomp($_);
		push(@lTranscripts, $_);
	}
	close (FILE);
}
elsif ($transcripts) {
	@lTranscripts = split(',', $transcripts);
}
elsif ($genes) {
	@lGenes = split(',', $genes);
}

my $hRes;
if ($transcripts) {
	foreach my $tr_name (@lTranscripts) {
		next unless ($tr_name =~ /ENST/);
		my ($transcript, $gene);
		eval {
			$transcript = $project->newTranscript($tr_name);
			$hRes->{$tr_name} = 'OK;'.$annotation;
		};
		if ($@) {
			$hRes->{$tr_name} = 'NOT;'.$annotation;
		}
		if ($transcript) {
			eval {
				$gene = $transcript->getGene;
				$hRes->{$tr_name} .= ';'.$gene->external_name();
			};
			if ($@) {
				$hRes->{$tr_name} .= ";<font color='orange'>?</font>";
			}
		}
		else {
			$hRes->{$tr_name} .= ';';
		}
	}
}
if ($genes) {
	foreach my $gene_name (@lGenes) {
		eval {
			my $gene = $project->newGene($gene_name);
			if ($gene->isGene()) {
				my $locus = 'chr'.$gene->getChromosome->id().':'.$gene->start().'-'.$gene->end();
				$hRes->{$gene_name} = 'OK;'.$genecode_version.';'.$gene->id().';'.$gene->external_name().';'.$locus;
			}
			else {
				$hRes->{$gene_name} = 'NOT; ---> '.$gene_name.' is not a GENE';
			}
		};
		if ($@) {
			$hRes->{$gene_name} = 'NOT;'.$annotation;
		}
	}
}

print $cgi->header('text/json-comment-filtered');
my $json_encode = encode_json $hRes;
print $json_encode;
exit(0);
