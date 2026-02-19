#!/usr/bin/perl
$|=1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb/polyviewer";
use lib "$Bin/../GenBo/lib/obj-nodb/packages/";
use lib "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/";
use lib "$Bin/../validation_variation/variation_editor_lib/";
use lib "$Bin/../packages/validation_variation/";
use lib "$Bin/../polydejavu/";

use GBuffer;
use html;
use polyviewer_css;
use Data::Dumper;
use polyviewer_html;
use GenBoNoSqlRocksTinyPolyviewerVariant;
use variations_editor_methods;
use update_variant_editor;
use PolyviewerVariant;
use dejavu_variants;
use Storable qw/freeze thaw nfreeze nstore_fd nstore retrieve/;
use JSON;


my $cgi = new CGI();

print $cgi->header('text/json-comment-filtered');
print "{\"progress\":\".";

my $login = $cgi->param('login');
my $pwd   = $cgi->param('pwd');
my $promoter_ai_value = $cgi->param('promoter_ai_value');

my $dejavu_variants = new dejavu_variants();
$dejavu_variants->user_name($login);
$dejavu_variants->pwd($pwd);
$dejavu_variants->hash_users_projects();
exit(0) if not $dejavu_variants->hash_users_projects();

print '.nb_proj.'.scalar(keys %{$dejavu_variants->hash_users_projects()});
$dejavu_variants->fork(8);

my $max_gnomad_ac = 30;
my $max_gnomad_ac_ho = 10;
my $max_dejavu = 30;
$promoter_ai_value = 0.7 if not $promoter_ai_value;

$dejavu_variants->max_gnomad_ac($max_gnomad_ac) if $max_gnomad_ac;
$dejavu_variants->max_gnomad_ac_ho($max_gnomad_ac) if $max_gnomad_ac_ho;
$dejavu_variants->min_promoter_ai($promoter_ai_value) if $promoter_ai_value;


my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $buffer->getRandomProjectName() );
my $parquet_promoter_ai = $project->get_promoterAI_filtred_parquet();

print '.duckdb.';
my @list_table_trio;
#my $sql_promoter_ai = "PRAGMA threads=6; SELECT rocksdb_id FROM read_parquet(['$parquet_promoter_ai'])";
my $sql_promoter_ai = "PRAGMA threads=6; SELECT rocksdb_id, geneid FROM read_parquet(['$parquet_promoter_ai']) WHERE ABS(promoterAI) >= $promoter_ai_value";
#warn $sql_promoter_ai;

my ($h_res_duck, $h_genes_only);
my $duckdb = $buffer->software('duckdb');
my $cmd = qq{set +H | $duckdb -json -c "$sql_promoter_ai"};
my $json_duckdb = `$cmd`;
if ($json_duckdb) {
	my $decode = decode_json $json_duckdb;
	my $i = 0;
	foreach my $h (@$decode) {
		if ($i == 1000) {
			print '.';
			$i = 0;
		}
		$h_res_duck->{$h->{rocksdb_id}}->{geneid} = $h->{geneid};
		$h_genes_only->{$h->{geneid}} = undef;
	}
}

$project->getChromosomes();
my $h_vector;
my $n = 0;
print '.vector.'.scalar(keys %$h_res_duck).'.';
foreach my $rocksid (keys %$h_res_duck) {
	$n++;
	if ($n == 100) {
		print '.';
		$n = 0;
	}
	my @l_tmp = split('!', $rocksid);
	my $chr_id = $l_tmp[0];
	my $pos = int($l_tmp[1]);
	
	if (not exists ($h_vector->{$chr_id})) {
		my $chr = $project->getChromosome($chr_id);
		$h_vector->{$chr_id} = Bit::Vector->new($chr->end() + 1);
	}
	$h_vector->{$chr_id}->Bit_On($pos);
}

print '.convert.';
my $ok = 0;
my $h_rocks_to_view;
foreach my $chr_id (sort keys %$h_vector) {
	print '.';
	#next if $chr_id ne '10';
	
	my $chr = $project->getChromosome($chr_id);
	my @list_ids = $h_vector->{$chr_id}->Index_List_Read();
	my $i = 0;
	foreach my $id (@list_ids) {
		$i++;
		if ($i == 500) {
			print '.';
			$i = 0;
		}
		my $res = $chr->rocks_dejavu->dejavu_interval($id -1 , $id +1);
		foreach my $dv_rocks_id (keys %{$res}) {
			my $nb_pat_he = 0;
			my $nb_pat_ho = 0;
			my $has_in_my_projects;
			foreach my $proj_id (keys %{$res->{$dv_rocks_id}}) {
				$has_in_my_projects = 1 if exists $dejavu_variants->hash_users_projects->{$proj_id};
				next if ($proj_id eq 'polybtf');
				$nb_pat_he += $res->{$dv_rocks_id}->{$proj_id}->{he};
				$nb_pat_ho += $res->{$dv_rocks_id}->{$proj_id}->{ho};
			}
			next if not $has_in_my_projects;
			my $nb_total = $nb_pat_he + $nb_pat_ho;
			next if $nb_total > $max_dejavu;
			if (exists $h_res_duck->{$chr_id.'!'.$dv_rocks_id}) {
				$h_res_duck->{$chr_id.'!'.$dv_rocks_id}->{dejavu} = $res->{$dv_rocks_id};
				$h_rocks_to_view->{$chr_id}->{$dv_rocks_id} = $res->{$dv_rocks_id};
				$ok++;
			}
		}
	}
}
print '.after_dv.'.$ok++.'.';


print '.checkvar.';
my ($hGenes, $hVariantsDetails) = $dejavu_variants->check_variants_from_gene($h_rocks_to_view);
print '...html...nbVar:'.scalar(keys %{$hVariantsDetails}).'.';
print '.nbGenes:'.scalar(keys %{$hGenes}).'.';

my $html;
$html .= qq{<table>};
$html .= qq{<tr>};
foreach my $gene_id (sort keys %{$hGenes}) {
	my @l_gene_id_tmp = split('_', $gene_id);
	next if ($h_genes_only and not exists $h_genes_only->{$gene_id} and not exists $h_genes_only->{$l_gene_id_tmp[0]});
	print '.';
	my @list_variants = keys %{$hGenes->{$gene_id}};
	$html .= $dejavu_variants->print_html_gene($gene_id, \@list_variants, $hVariantsDetails);
}
$html .= qq{</tr>};
$html .= qq{</table>};

my $hRes;
$hRes->{html} = $html;
my $json_encode = encode_json $hRes;
print ".\",";
$json_encode =~ s/{//;
print $json_encode;
exit(0);



