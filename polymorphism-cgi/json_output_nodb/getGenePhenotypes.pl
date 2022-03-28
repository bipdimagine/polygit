#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";
use GBuffer;
use Digest::MD5 qw(md5_hex);
use JSON;


my $h_excepted_words = {
	'of' => 1,
	'to' => 1,
	'on' => 1,
	'for' => 1,
	'all' => 1,
	'the' => 1,
	'and' => 1,
	'&' => 1,
	'&&' => 1,
	'in' => 1,
	'or' => 1,
	'with' => 1,
	'without' => 1,
	'i' => 1,
	'ii' => 1,
	'iii' => 1,
	'iv' => 1,
	'v' => 1,
	'vi' => 1,
};

my $cgi = new CGI();
my $projectName = $cgi->param('project');
my $geneId = $cgi->param('gene');

my $buffer = new GBuffer;
unless ($projectName) {
	$projectName = $buffer->get_random_project_name_with_this_annotations_and_genecode();
}
my $project = $buffer->newProject(-name => $projectName);
my $gene = $project->newGene($geneId);

warn ref($gene);

my $hGenePheno = $gene->hash_phenotypes();
my @lCat = sort keys %$hGenePheno;

my ($h, @list_pheno);
$h->{id} = $geneId;
$h->{name} = $geneId;
$h->{description} = '-';
$h->{description} = $gene->description();
my $out = qq{<table class="table" style="border: solid black 1px;">};
$out .= qq{<tbody>};
$out .= qq{<tr>};
$out .= qq{<td><span style="color:green;font-size:12px;">DESCRIPTION</span></td>};
$out .= qq{<td><span id='span_gene_phenotype_description' style='font-size:12px;'>}.$h->{description}.qq{</span></td>};
$out .= qq{</tr>};
foreach my $cat (sort keys %{$hGenePheno}) {
	my $pheno = join (" ;", sort @{$hGenePheno->{$cat}});
	$h->{$cat} = $pheno;
	push(@list_pheno, $pheno);
	my $td_id = 'span_gene_phenotype_'.$cat;
	my $cat_txt = uc($cat);
	$out .= qq{<tr>};
	$out .= qq{<td><span style="color:green;font-size:12px;">$cat_txt</span></td>};
	$out .= qq{<td><span id='$td_id' style='font-size:12px;'>$pheno</span></td>};
	$out .= qq{</tr>};
}
$out .= qq{</tbody>};
$out .= qq{</table>};
$h->{html} = $out;

my $h_words;
foreach my $cat (@lCat) {
	my $txt = $h->{$cat};
	$txt =~ s/\n/ /g;
	$txt =~ s/[;\.,()]/ /g;
	$txt =~ s/\s+/ /g;
	$txt =~ s/\[.+\]//g;
	$txt = lc($txt);
	foreach my $word (split(' ', $txt)) {
		next unless ($word =~ /[a-z]/);
		next if (exists $h_excepted_words->{$word});
		$h_words->{$word}->{nb}++;
		$h_words->{$word}->{cat}->{$cat}++;
	}
}

my @l_h_words;
foreach my $word (keys %$h_words) {
	my $h;
	$h->{x} = $word;
	$h->{value} = $h_words->{$word}->{nb};
	$h->{category} = join('/', sort keys %{$h_words->{$word}->{cat}});
	push(@l_h_words, $h)
}
$h->{txt_all} = '-';
if (scalar(@l_h_words) > 0) { $h->{txt_all} = encode_json \@l_h_words; }
					
my @lItems;
push(@lItems, $h);

my $hashRes;
$hashRes->{'label'} = 'line';
$hashRes->{'items'} = \@lItems;
print $cgi->header('text/json-comment-filtered');
print encode_json $hashRes;
exit(0);