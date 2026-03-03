#!/usr/bin/perl
$|=1;
use POSIX qw(SIGPIPE);
$SIG{SIGPIPE} = 'IGNORE';

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
use List::MoreUtils qw(natatime);
use JSON;
use MCE::Loop;


my $fork = 5;

my $cgi = new CGI();

my $outfile;
my $launch_job = $cgi->param('launch_job');

if ($launch_job) {
	use POSIX qw(setsid);
	use JSON;
	use File::Path qw(make_path);

	my $job_id = $cgi->param('job_id');
	my $job_dir = "/tmp/polyweb_jobs";
	make_path($job_dir) unless -d $job_dir;
	
	# =========================
	# MODE POLLING (attente résultat)
	# =========================
	if ($job_id) {
	    my $file = "$job_dir/$job_id.json";
	    print $cgi->header('application/json');
	
	    if (-e $file) {
	        open(my $fh, "<", $file);
	        local $/;
	        my $json = <$fh>;
	        close $fh;
	        print $json;
	    }
	    else {
	        print encode_json({ status => "running" });
	    }
	    exit(0);
	}
	
	# =========================
	# MODE LANCEMENT JOB
	# =========================
	
	$job_id = time() . "_" . $$ . "_" . int(rand(10000));
	$outfile = "$job_dir/$job_id.json";
	
	print $cgi->header('application/json');
	print encode_json({ status => "started", job_id => $job_id });
	
	my $pid = fork();
	exit(0) if $pid;   # le parent sort immédiatement
	
	# ===== ENFANT =====
	setsid();
	open STDIN,  '<', '/dev/null';
	open STDOUT, '>', '/dev/null';
	open STDERR, '>', '/dev/null';
}

print $cgi->header('text/json-comment-filtered');
print "{\"progress\":\".";

my $login = $cgi->param('login');
my $pwd   = $cgi->param('pwd');
my $promoter_ai_value = $cgi->param('promoter_ai_value');
my $max_gnomad_ac = $cgi->param('gnomad');
my $max_gnomad_ac_ho = $cgi->param('gnomad_ho');
my $max_dejavu = $cgi->param('dejavu');
my $max_dejavu_ho = $cgi->param('dejavu_ho');
my $min_ratio = $cgi->param('min_perc_all');
my $only_my_projects = $cgi->param('only_my_projects');
my $only_ill = $cgi->param('only_ill');
my $only_strict_ill = $cgi->param('only_strict_ill');
my $models = $cgi->param('models');

my $dejavu_variants = new dejavu_variants();
$dejavu_variants->user_name($login);
$dejavu_variants->pwd($pwd);

if ($only_my_projects) {
	$dejavu_variants->hash_users_projects() if $only_my_projects;
	exit(0) if not $dejavu_variants->hash_users_projects();
	print '.nb_proj.'.scalar(keys %{$dejavu_variants->hash_users_projects()});
}
else {
	$dejavu_variants->{hash_users_projects}->{all} = 1;
	print '.all_projects.';
}

$dejavu_variants->fork(10);

$dejavu_variants->max_gnomad_ac($max_gnomad_ac) if $max_gnomad_ac;
$dejavu_variants->max_gnomad_ac_ho($max_gnomad_ac) if $max_gnomad_ac_ho;
$dejavu_variants->min_promoter_ai($promoter_ai_value) if $promoter_ai_value;
$dejavu_variants->min_ratio($min_ratio) if $min_ratio;
$dejavu_variants->only_ill_patients(1) if $only_ill;
$dejavu_variants->only_strict_ill_patients(1) if $only_strict_ill;


if ($models) {
	my $h_models;
	foreach my $model_name (split(',', $models)) {
		$h_models->{$model_name} = 1;
	}
	$dejavu_variants->{models} = $h_models
}


my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $buffer->getRandomProjectName() );
my $parquet_promoter_ai = $project->get_promoterAI_filtred_parquet();

print '.duckdb.';
my @list_table_trio;
my $sql_promoter_ai = "PRAGMA threads=6; SELECT rocksdb_id, geneid FROM read_parquet(['$parquet_promoter_ai']) WHERE ABS(promoterAI) >= $promoter_ai_value";
my ($h_res_duck, $h_genes_only);
my $duckdb = $buffer->software('duckdb');
my $cmd = qq{set +H | $duckdb -json -c "$sql_promoter_ai"};
my $json_duckdb = `$cmd`;
print '.|.';
if ($json_duckdb) {
	my $decode = decode_json $json_duckdb;
	my $i = 0;
	foreach my $h (@$decode) {
		$i++;
		if ($i == 100000) {
			print '.';
			$i = 0;
		}
		$h_res_duck->{$h->{rocksdb_id}}->{geneid} = $h->{geneid};
		$h_genes_only->{$h->{geneid}} = undef;
	}
}
print '.@.';

my $h_vector;
my $n = 0;
print '.vector.'.scalar(keys %$h_res_duck).'.';
foreach my $rocksid (keys %$h_res_duck) {
	$n++;
	if ($n == 15000) {
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
	my $chr = $project->getChromosome($chr_id);
	my @list_ids = $h_vector->{$chr_id}->Index_List_Read();
	my $i = 0;
	foreach my $id (@list_ids) {
		$i++;
		if ($i == 2500) {
			print '.';
			$i = 0;
		}
		my $res = $chr->rocks_dejavu->dejavu_interval($id -1 , $id +1);
		foreach my $dv_rocks_id (keys %{$res}) {
			my $nb_pat_he = 0;
			my $nb_pat_ho = 0;
			my $has_in_my_projects;
			foreach my $proj_id (keys %{$res->{$dv_rocks_id}}) {
				$has_in_my_projects = 1 if $only_my_projects and exists $dejavu_variants->hash_users_projects->{$proj_id};
				next if ($proj_id eq 'polybtf');
				$nb_pat_he += $res->{$dv_rocks_id}->{$proj_id}->{he};
				$nb_pat_ho += $res->{$dv_rocks_id}->{$proj_id}->{ho};
			}
			next if $only_my_projects and not $has_in_my_projects;
			my $nb_total = $nb_pat_he + $nb_pat_ho;
			next if $nb_total > $max_dejavu;
			next if $nb_pat_ho > $max_dejavu_ho;
			if (exists $h_res_duck->{$chr_id.'!'.$dv_rocks_id}) {
				$h_res_duck->{$chr_id.'!'.$dv_rocks_id}->{dejavu} = $res->{$dv_rocks_id};
				$h_rocks_to_view->{$chr_id}->{$dv_rocks_id} = $res->{$dv_rocks_id};
				$ok++;
			}
		}
	}
}
print '.after_dv.'.$ok++.'.';


$dejavu_variants->fork(5);
print '.checkvar.';

my ($hGenes, $hVariantsDetails) = $dejavu_variants->check_variants_from_gene($h_rocks_to_view);
print '...html...nbVar:'.scalar(keys %{$hVariantsDetails}).'.';
my $nb_genes = scalar(keys %{$hGenes});
print '.nbGenes:'.$nb_genes.'.';


MCE::Loop->init(
   max_workers => $fork, chunk_size => 'auto'
);
my @results = mce_loop {
	my ($mce, $chunk_ref, $chunk_id) = @_;
	my $list;
	foreach my $gene_id (@$chunk_ref) {
		my @l_gene_id_tmp = split('_', $gene_id);
		next if ($h_genes_only and not exists $h_genes_only->{$gene_id} and not exists $h_genes_only->{$l_gene_id_tmp[0]});
		print '.';
		my @list_variants = keys %{$hGenes->{$gene_id}};
		my $this_html;
		eval {
			$this_html = $dejavu_variants->print_html_gene($gene_id, \@list_variants, $hVariantsDetails);
		};
		if ($@) {
			$this_html = qq{ERROR with $gene_id};
		}
		push(@$list, $this_html);
    }
	MCE->gather($list);
} sort keys %{$hGenes};	

my @list_html_genes;
foreach my $t (@results) {
	push(@list_html_genes,@$t);
}

my $html;
$html .= qq{<table>};
$html .= qq{<tr>};
$html .= join('', @list_html_genes);
$html .= qq{</tr>};
$html .= qq{</table>};

if ($launch_job) {
	my $hRes;
	$hRes->{status} = "finished";
	$hRes->{html}   = $html;
	
	open(my $out, ">", $outfile);
	print $out encode_json($hRes);
	close $out;
	
	POSIX::_exit(0);
}
else {
	my $hRes;
	$hRes->{html} = $html;
	my $json_encode = encode_json $hRes;
	print ".\",";
	$json_encode =~ s/{//;
	print $json_encode;
	
	POSIX::_exit(0);
}


