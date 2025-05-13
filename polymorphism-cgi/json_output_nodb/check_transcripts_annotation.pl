#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use Data::Dumper;
use JSON;
use Carp;

use lib "$Bin/../GenBo/lib/obj-nodb/";
use GBuffer;
use GenBoProject;


my $cgi = new CGI();
my $annotation	= $cgi->param('annotation');
my $file_name	= $cgi->param('file');
my $transcripts	= $cgi->param('tr');
my $genes	= $cgi->param('genes');
my $release	= $cgi->param('release');

confess() unless $release;

my $buffer = GBuffer->new();
my $project_name = $buffer->getRandomProjectName($release, $annotation);
my $project = $buffer->newProject(-name=>$project_name);

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
				my $remap = $transcript->remap_status();
				if ($transcript->is_partial_transcript()) {
					$remap = 'partial';
					my $h_patial_infos = $transcript->hash_partial_infos->{intspan};
					my $is_frameshift;
					my @lNt = sort {$a <=> $b} keys %{$h_patial_infos};
					if ($lNt[0] == 0) { shift(@lNt); }
					my @lNt_text;
					foreach my $nt (@lNt) {
						if ($nt % 3 > 0) { $is_frameshift = 1; }
						push (@lNt_text, '+'.$nt.'nt');
					}
					$remap .= " FRAMESHIFT " if ($is_frameshift);
					$remap .= " [".join(', ', @lNt_text)."]";
				}
				$hRes->{$tr_name} .= ';'.$remap;
			};
			if ($@) {
				warn "\n\n";
				warn Dumper $@;
				warn Dumper $transcript->hash_partial_infos;
				warn "\n\n";
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
				$hRes->{$gene_name} = 'OK;'.$annotation.';'.$gene->id().';'.$gene->external_name().';'.$locus;
			}
			else {
				$hRes->{$gene_name} = 'NOT; ---> '.$gene_name.' is not a GENE';
			}
		};
		if ($@) {
			#warn Dumper $@;
			$hRes->{$gene_name} = 'NOT;'.$annotation;
		}
	}
}

print $cgi->header('text/json-comment-filtered');
my $json_encode = encode_json $hRes;
print $json_encode;
exit(0);
