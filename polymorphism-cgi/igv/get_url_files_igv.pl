#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/obj-nodb";
use GBuffer;
use JSON;

my $cgi = new CGI();

my $project_name = $cgi->param('project');
my $buffer = GBuffer->new();
my $project = $buffer->newProject(-name=>$project_name);

my $url_gene_bed = $project->get_gtf_genes_annotations_igv();
$url_gene_bed =~ s/\/data-isilon//;

my $url_genome_fasta = $project->getGenomeFasta();
$url_genome_fasta =~ s/\/data-isilon//;

my $hash;
$hash->{'url_gene_bed'} = $url_gene_bed;
$hash->{'url_genome_fasta'} = $url_genome_fasta;
#$hash->{'url_cytoband'} = $project->();

print $cgi->header('text/json-comment-filtered');
print encode_json $hash;
exit(0);