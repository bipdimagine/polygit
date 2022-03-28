#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
#use lib "/software/polyweb/poly-src/GenBo/lib/obj-nodb/";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag";

#use lib "/bip-d/soft/distrib/tabix/latest/perl";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/coverage";
use lib "$Bin/../packages/validation_variation"; 
use lib "$Bin/../packages/cache"; 

use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use draw_cnv; 
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/utility.pm";
use html; 
use infos_coverage_exons;
use JSON::XS;
use image_coverage;
#use Set::;
use Carp;
use export_data;
use strict;
use Set::IntSpan::Fast::XS;
use Data::Dumper;
use GBuffer;
use Getopt::Long;
use Carp;
use Set::Intersection;
use Tabix;
use Storable qw/thaw/;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use coverage;
use Spreadsheet::WriteExcel;
use POSIX;
use validationQuery;
use Date::Tiny;
use JSON::XS;
use List::MoreUtils qw{part};
#use PDF::API2;
#use PDF::Table;
use constant mm => 25.4 / 72;
use constant in => 1 / 72;
use constant pt => 1;
use Time::HiRes qw ( time alarm sleep );
my $cgi          = new CGI();

#CD000062
my $buffer = GBuffer->new();
my $project_name = $cgi->param('project');
warn $project_name;
#$project_name = "NGS2018_2181";
my $project;
#if ($cgi->param('cnv_coverage') ne 1){
$project = $buffer->newProject( -name 			=> $project_name, -typeFilters=>'individual' ,-cgi_object=>1);

my $is_polyquery_view = $cgi->param('polyquery_view');
my $gene_id = $cgi->param('gene');
my $vid = $cgi->param('vid');
my $hid = $cgi->param('hid');



my  $CSS = qq{
	<style type="text/css"> 
.bs-callout {
    padding: 20px;
    margin: 20px 0;
    border: 1px solid #eee;
    border-left-width: 5px;
    border-radius: 3px;
}
.bs-callout h4 {
    margin-top: 0;
    margin-bottom: 5px;
}
.bs-callout p:last-child {
    margin-bottom: 0;
}
.bs-callout code {
    border-radius: 3px;
}
.bs-callout+.bs-callout {
    margin-top: -5px;
}
.bs-callout-default {
    border-left-color: #777;
}
.bs-callout-default h4 {
    color: #777;
}
.bs-callout-primary {
    border-left-color: #428bca;
}
.bs-callout-primary h4 {
    color: #428bca;
}
.bs-callout-success {
    border-left-color: #5cb85c;
}
.bs-callout-success h4 {
    color: #5cb85c;
}
.bs-callout-danger {
    border-left-color: #d9534f;
}
.bs-callout-danger h3 {
    color: #d9534f;
}
.bs-callout-danger h4 {
    color: black;
}
.bs-callout-warning {
    border-left-color: #f0ad4e;
}
.bs-callout-warning h4 {
    color: #f0ad4e;
}
.bs-callout-info {
    border-left-color: #5bc0de;
}
.bs-callout-info h4 {
    color: #5bc0de;
}
</style>
};

html::print_header_polydiag($cgi) unless ($is_polyquery_view);
print $CSS;

if ($gene_id) {
	my $query = $buffer->queryHgmd();
	my $h = $query->getDataHGMDPro_gene_infos($gene_id);
	
	my $gene_name = $h->{$gene_id}->{gene};
	my $disease = $h->{$gene_id}->{disease};
	my $inheritance = $h->{$gene_id}->{inheritance};
	
	print qq{
		<div class="bs-callout bs-callout-danger">
	  	<h3><i class="fa fa-dot-circle"></i>&nbsp<big>$gene_name - $disease : $inheritance</big>&nbsp<i class="fa fa-dot-circle"></i></h3> 
	};
	
	my $omimid = '.';
	$omimid = $h->{$gene_id}->{omimid} if (exists $h->{$gene_id}->{omimid});
	my $hgncID = '.';
	$hgncID = $h->{$gene_id}->{hgncID} if (exists $h->{$gene_id}->{hgncID});
	my $entrezID = '.';
	$entrezID = $h->{$gene_id}->{entrezID} if (exists $h->{$gene_id}->{entrezID});
	my $gdbID = '.';
	$gdbID = $h->{$gene_id}->{gdbid} if (exists $h->{$gene_id}->{gdbid});
	my $date = $h->{$gene_id}->{gene_date};
	my $chrom = $h->{$gene_id}->{chrom};
	my $refseq = $h->{$gene_id}->{refseq};
	my $alt_names = $h->{$gene_id}->{altsymbol};
	$alt_names =~ s/\|/, /g;
	my $go_terms_name = $h->{$gene_id}->{go_terms_name};
	$go_terms_name =~ s/\|/, /g;
	
	print qq{
		<h4>OMIM:$omimid HGNC:$hgncID ENTREZ:$entrezID GDB:$gdbID </h4>
		<h4>DATE:  $date </h4>
		<h4>ALT NAMES: $alt_names </h4> 
		<h4>Refseq: $refseq   </h4>
	 	GO Terms Name: $go_terms_name
	</div>
	};
}


elsif ($vid) {
	$project->changeAnnotationVersion($buffer->get_polybtf_default_release(), 1);
	my $data = $buffer->queryHgmd()->getDataHGMDPro("$hid");
	my $acc_num  =  $data->{acc_num};
	my $v = $project->_newVariant($vid);
	$v->hgmd_id($acc_num);
	my $d =  $data->{disease};;
	my $tag =  $data->{tag};
	my $date =  $data->{new_date};
	my $pmid  =  $data->{pmid};
	my $pubmed = "https://www.ncbi.nlm.nih.gov/pubmed/".$pmid;
	my $title =  "<a href=\"$pubmed\" target=\"_blank\" >".$data->{title}."</a>";
	my $gene = $data->{gene};
	my $comment =  $data->{comment};
	my $Inheritance =  $data->{inheritance};
	my $name = $v->name;
	my $pos = $v->getChromosome()->ucsc_name.":".$v->start."-".$v->start;	
	my $cyto = $data->{chrom};
	my $refseq = $data->{refseq};
	my $genedesc = $data->{genename};
	my $omim = "https://www.omim.org/entry/".$data->{omimid};
	my $omimh =  "<a href=\"$omim\" target=\"_blank\" >".$data->{omimid}."</a>";
	my $dbsnp = "https://www.ncbi.nlm.nih.gov/snp/".$data->{dbsnp};
	my $dbsnph =   "<a href=\"$dbsnp\" target=\"_blank\">".$data->{dbsnp}."</a>";
	my $hgvs = $data->{hgvs}." ".$data->{hgvsAll}." ".$data->{amino};
	my $refseq = "https://www.ncbi.nlm.nih.gov/nuccore/".$data->{refseq};
	my $refseqh =   "<a href=\"$refseq\" target=\"_blank\">".$data->{refseq}."</a>";
	my $releases = $v->hgmd_releases();
	my $gene_name;
	eval { $gene_name = $v->hgmd_gene_name(); };
	
	print qq{
	<div class="bs-callout bs-callout-danger">
	  <h3><i class="fa fa-dot-circle"></i>&nbsp<big>$tag  $d : $Inheritance&nbsp<br><br>&nbspGene: $gene_name</big>&nbsp<i class="fa fa-dot-circle"></i></h3> 
	};
	
	print qq{  <h4>  Ref : $title</h4>};
	foreach my $pmid (keys %{$data->{pubmed}}){
		my $pubmed = "https://www.ncbi.nlm.nih.gov/pubmed/".$pmid;
		my $title =  "<a href=\"$pubmed\" target=\"_blank\" >".$data->{pubmed}->{$pmid}->{title}."</a>";
		print qq{  <h4>  Ref : $title</h4>};
	}
	
	print qq{
	  <h4>OMIM: $omimh dbsnp:$dbsnph HGMD: $acc_num ($date)</h4>
	  <h4>$cyto: $pos </h4> 
	  <h4>$gene: $genedesc </h4>
	  <h4>Refseq: $refseqh   </h4>
	  <h4>Release(s): $releases   </h4>
	 Comment : $comment
	</div>
	};
}

exit(0);


