#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;

use strict;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/validation_variation"; 
use lib "$Bin/../cache_nodb/scripts/";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use GBuffer;
use export_data;
use JSON;
use polyweb_dude;
use VcfMerge;
use GenBoNoSql;
use Set::IntervalTree;
use Spreadsheet::WriteExcel;
use Bit::Vector;
use Bit::Vector::Overload;
use Compress::Snappy;
use Storable qw(store retrieve freeze dclone thaw);
use POSIX qw(strftime);
use List::MoreUtils qw(natatime);
use CGI::Session;
use html; 
use Carp;
use Cache_Commons;
use QueryVectorFilter;
use IO::Handle;
use xls_export;
use session_export;

require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/html_polygenescout.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";

my $io = IO::Handle->new();
$io->autoflush(1);


my $cgi = new CGI();
my $gene_id = $cgi->param('gene');

my $buffer_init = new GBuffer;
my $project_init_name = $buffer_init->get_random_project_name_with_this_annotations_and_genecode();
my $project_init = $buffer_init->newProjectCache( -name => $project_init_name);
my $genomeFai_init = $project_init->getGenomeFai();
my $gene = $project_init->newGene($gene_id);

my $html = qq{<div class="input-group" style="width:100%">};
$html .= qq{<select class="form-control" id="form_transcripts_exons_introns" style="font-size:9px;height:auto;width:100%;">};
$html .= qq{<option value=''><span>ALL Gene</span></option>};
my $h_genes_transcripts_exons;
foreach my $tr (@{$gene->getTranscripts()}) {
	my $tr_all_id = $tr->id();
	my @lTmp = split('_', $tr_all_id);
	my $tr_id = $lTmp[0];
	my $locus_tr = $gene->getChromosome->id().':'.$tr->start().'-'.$tr->end();
	$h_genes_transcripts_exons->{$tr_id}->{'locus'} = $locus_tr;
	$h_genes_transcripts_exons->{$locus_tr} = $tr_id;
	$html .= qq{<option value='$locus_tr'><span>$tr_id</span></option>};
	
	my ($h_e, $h_i);
	foreach my $exon (@{$tr->getExons()}) {
		my $locus = $gene->getChromosome->id().':'.$exon->start().'-'.$exon->end();
		$h_genes_transcripts_exons->{$tr_id}->{exons}->{$exon->id()} = $locus;
		$h_genes_transcripts_exons->{$locus} = $exon->id()  if (not exists $h_genes_transcripts_exons->{$locus});
		my $id = $exon->id();
		$id =~ s/$tr_all_id//;
		$id =~ s/_//;
		my $this_html = qq{<option value='$locus'><span>$tr_id $id</span></option>};
		my $nb_exon = $id;
		$nb_exon =~ s/ex//;
		$h_e->{$nb_exon} = $this_html;
	}
	foreach my $intron (@{$tr->getIntrons()}) {
		my $locus = $gene->getChromosome->id().':'.$intron->start().'-'.$intron->end();
		$h_genes_transcripts_exons->{$tr_id}->{introns}->{$intron->id()} = $locus;
		$h_genes_transcripts_exons->{$locus} = $intron->id() if (not exists $h_genes_transcripts_exons->{$locus});
		my $id = $intron->id();
		$id =~ s/$tr_all_id//;
		$id =~ s/_//;
		my $this_html = qq{<option value='$locus'><span>$tr_id $id</span></option>};
		my $nb_intron = $id;
		$nb_intron =~ s/intron//;
		$h_i->{$nb_intron} = $this_html;
	}
	foreach my $i (sort {$a <=> $b} keys %$h_e) { $html .= $h_e->{$i}; }
	foreach my $i (sort {$a <=> $b} keys %$h_i) { $html .= $h_i->{$i}; }
}
$html .= qq{</select>};
$html .= qq{</div>};


my $hRes;
$hRes->{html_list} = $html;
$hRes->{transcripts_exons} = $h_genes_transcripts_exons;
print $cgi->header('text/json-comment-filtered');
my $json_encode = encode_json $hRes;
print $json_encode;
exit(0);