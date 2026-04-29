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
use Set::IntSpan;

my $cgi = new CGI();

my $fork = 10;
my $launch_job = $cgi->param('launch_job');

my $outfile;
my $min_cadd_for_ncboost = 25;

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
my $ncboost_value = $cgi->param('ncboost_value');
my $max_gnomad_ac = $cgi->param('gnomad');
my $max_gnomad_ac_ho = $cgi->param('gnomad_ho');
my $max_dejavu = $cgi->param('dejavu');
my $max_dejavu_ho = $cgi->param('dejavu_ho');
my $min_ratio = $cgi->param('min_perc_all');
my $only_my_projects = $cgi->param('only_my_projects');
my $only_ill = $cgi->param('only_ill');
my $only_strict_ill = $cgi->param('only_strict_ill');
my $models = $cgi->param('models');
my $region = $cgi->param('region');

my $dejavu_variants = new dejavu_variants();
$dejavu_variants->user_name($login);
$dejavu_variants->pwd($pwd);

$dejavu_variants->hash_users_projects() if $only_my_projects;
exit(0) if not $dejavu_variants->hash_users_projects();
print '.nb_proj.'.scalar(keys %{$dejavu_variants->hash_users_projects()});
if (not $only_my_projects) {
	$dejavu_variants->{hash_users_projects}->{all} = 1;
	print '.all_projects.';
}

$dejavu_variants->fork($fork);

$dejavu_variants->max_dejavu($max_dejavu) if $max_dejavu;
$dejavu_variants->max_dejavu_ho($max_dejavu_ho) if $max_dejavu_ho;
$dejavu_variants->max_gnomad_ac($max_gnomad_ac) if $max_gnomad_ac;
$dejavu_variants->max_gnomad_ac_ho($max_gnomad_ac_ho) if $max_gnomad_ac_ho;
$dejavu_variants->min_ratio($min_ratio) if $min_ratio;
$dejavu_variants->only_ill_patients(1) if $only_ill;
$dejavu_variants->only_strict_ill_patients(1) if $only_strict_ill;

my $h_models;
if ($models) {
	foreach my $model_name (split(',', $models)) {
		$h_models->{$model_name} = 1;
	}
	$dejavu_variants->{models} = $h_models
}

my $buffer = new GBuffer;
my $project_name = $buffer->getRandomProjectName();
my $project = $buffer->newProject( -name => $project_name );


my $filters_cons = $cgi->param('filters_cons');
my $h_filters_cons;
if ($filters_cons) {
	foreach my $cons (split(',', $filters_cons)) {
		$dejavu_variants->{hash_filters_cons}->{lc($cons)} = undef;
	}
}

my $path_ncboost_global = $project->ncboost_parquet_path();
my $path_ncboost_dejavu = $project->ncboost_parquet_dejavu_filter_path();
my $parquet_promoter_ai = $project->get_promoterAI_parquet();
my $parquet_promoter_ai_filtred = $project->get_promoterAI_filtred_parquet();

my ($h_res_duck, $h_genes_only, $h_ncboost_values);
my ($h_res_duck_promoter_ai, $h_genes_only_promoter_ai);
my ($h_res_duck_ncboost);

### PART 1 - Select variants with promoterAI / NCBOOST


my $h_rocks_to_view;
if ($promoter_ai_value) {
	launch_promoter_ai ($dejavu_variants, $promoter_ai_value);
	print '.promoter_ai.'.scalar(keys %{$h_res_duck});

	my $h_vector;
	my $n = 0;
	print '.intspan.'.scalar(keys %$h_res_duck).'.';
	
	foreach my $rocksid (keys %$h_res_duck) {
		my @l_tmp = split('!', $rocksid);
		my $chr_id = $l_tmp[0];
		my $pos = int($l_tmp[1]);
		my $h_gad = $project->getChromosome($chr_id)->rocksdb('gnomad')->value($l_tmp[1].'!'.$l_tmp[2]);
		if ($h_gad and exists $h_gad->{ac} and $h_gad->{ac} > $max_gnomad_ac) {
			delete $h_res_duck->{$rocksid};
			next;
		}
		if ($h_gad and exists $h_gad->{ho} and $h_gad->{ho} > $max_gnomad_ac_ho) {
			delete $h_res_duck->{$rocksid};
			next;
		}
		$h_vector->{$chr_id}->{$pos} = 1;
		$n++;
		print '.' if $n % 10000 == 0;
	}
	
	
	### PART 2 - Intersect variants filters and DejaVu
	
	print '.convert.'.$n.'.';
	my $ok = 0;
	my @lChr_vec = sort keys %$h_vector;
	
	MCE::Loop->init(
	   max_workers => $fork,
	   chunk_size => 1,
	);
	my @results_convert = mce_loop {
		my ($mce, $chunk_ref, $chunk_id) = @_;
		my $b = new GBuffer;
		my $pr = $b->newProject(-name => $project_name);
		my ($local_h_res, $local_h_rocks, $local_ok);
		foreach my $chr_id (@$chunk_ref) {
			print '.';
			my $chr = $pr->getChromosome($chr_id);
			#my @list_ids = $h_vector->{$chr_id}->Index_List_Read();
			my @list_ids = sort keys %{$h_vector->{$chr_id}};
			
			my $i = 0;
			my $no = $chr->rocks_dejavu();
			
			foreach my $id (@list_ids) {
				$i++;
				if ($i == 2500) {
					print '.';
					$i = 0;
				}
				my $res = $no->dejavu_interval($id -1 , $id +1);
				foreach my $dv_rocks_id (keys %{$res}) {
					my $nb_pat_he = 0;
					my $nb_pat_ho = 0;
					my $has_in_my_projects;
					$has_in_my_projects = 1 if not $only_my_projects;
					foreach my $proj_id (keys %{$res->{$dv_rocks_id}}) {
						$has_in_my_projects = 1 if exists $dejavu_variants->{hash_users_projects}->{$proj_id};
						next if ($proj_id eq 'polybtf');
						$nb_pat_he += $res->{$dv_rocks_id}->{$proj_id}->{he};
						$nb_pat_ho += $res->{$dv_rocks_id}->{$proj_id}->{ho};
					}
					next if $only_my_projects and not $has_in_my_projects;
					my $nb_total = $nb_pat_he + $nb_pat_ho;
					next if $nb_total > $max_dejavu;
					next if $nb_pat_ho > $max_dejavu_ho;
					if (exists $h_res_duck->{$chr_id.'!'.$dv_rocks_id}) {
						$local_h_res->{$chr_id.'!'.$dv_rocks_id}->{dejavu} = $res->{$dv_rocks_id};
						$local_h_rocks->{$chr_id}->{$dv_rocks_id} = $res->{$dv_rocks_id};
						$local_ok++;
					}
				}
			}
			$no->close();
		}
	    MCE->gather({
	        h_res   => $local_h_res,
	        h_rocks => $local_h_rocks,
	        ok      => $local_ok
	    });
	} @lChr_vec;	
		
	foreach my $data (@results_convert) {
	    $ok += $data->{ok};
	    if ($data->{h_res}) {
	        foreach my $k (keys %{$data->{h_res}}) {
	            $h_res_duck->{$k} = $data->{h_res}->{$k};
	        }
	    }
	    if ($data->{h_rocks}) {
	        foreach my $chr_id (keys %{$data->{h_rocks}}) {
	            foreach my $id (keys %{$data->{h_rocks}->{$chr_id}}) {
	                $h_rocks_to_view->{$chr_id}->{$id} = $data->{h_rocks}->{$chr_id}->{$id};
	            }
	        }
	    }
	}
	print '.after_dv.'.$ok++.'.';
	MCE::Loop->finish();
}
if ($ncboost_value) {
	launch_ncboost ($dejavu_variants, $ncboost_value, $max_dejavu, $max_dejavu_ho);
	$dejavu_variants->hash_ncboost_values($h_ncboost_values);
	print '.ncboost.'.scalar(keys %{$h_res_duck_ncboost});
	my $h_proj;
	foreach my $id (keys %{$dejavu_variants->hash_users_projects()}) {
		next if $id =~ /NGS/;
		$h_proj->{$id} = undef;
	}
	my $i = 0;
	my $found = 0;
	foreach my $id (keys %{$h_res_duck_ncboost}) {
		$i++;
		if ($i == 1000) {
			print '.';
			$i = 0;
		}
		my @ltmp = split('!', $id);
#		if (exists $dejavu_variants->hash_users_projects->{all}) {
#			my $h_dv = $project->getChromosome($ltmp[0])->rocks_dejavu->dejavu($ltmp[1].'!'.$ltmp[2]);
#			my ($h_proj);
#			my $nb_he = 0;
#			my $nb_ho = 0;
#			foreach my $proj_id (keys %{$h_dv}) {
#				$h_proj->{$proj_id} = undef;
#				$nb_he += $h_dv->{$proj_id}->{he};
#				$nb_ho += $h_dv->{$proj_id}->{ho};
#			}
#			next if (($nb_he + $nb_ho) >= $max_dejavu);
#			next if ($nb_ho >= $max_dejavu_ho);
#			$h_rocks_to_view->{$ltmp[0]}->{$ltmp[1].'!'.$ltmp[2]} = $h_proj;
#			$found++;
#		}
#		else {
			$h_rocks_to_view->{$ltmp[0]}->{$ltmp[1].'!'.$ltmp[2]} = $h_proj;
			$found++;
#		}
	}
	print '.now.'.$found++.'.';
}

if ($region) {
	if (not $promoter_ai_value and not $ncboost_value) {
		print '.only_region.';
		my ($chr_filter, $start_filter, $end_filter) = split('-', $region);
		$dejavu_variants->{only_chromosome} = $chr_filter;
		my $sql_pos;
		if ($start_filter and $end_filter) {
			$sql_pos = "and pos38 >= $start_filter and pos38 <= $end_filter";
		}
		my $i = 0;
		my $sql_parquets = $dejavu_variants->sql_projects_parquet();
		
		my $sql = "
			PRAGMA threads=$fork;
			WITH base AS ( SELECT * FROM $sql_parquets WHERE chr38='$chr_filter' $sql_pos ),
			agg AS (
			    SELECT 
			        chr38, pos38, allele,
			        SUM(he) AS sum_he,
			        SUM(ho) AS sum_ho
			    FROM base
			    GROUP BY chr38, pos38, allele
			    HAVING (SUM(he) + SUM(ho)) <= $max_dejavu AND SUM(ho) <= $max_dejavu_ho
			)
			SELECT 
			    b.project, b.chr38, b.pos38, b.allele, b.he, b.ho
			FROM base b
			JOIN agg USING (chr38, pos38, allele);
		";
		
		
		my $duckdb = $dejavu_variants->buffer->software('duckdb');
		open(my $fh, "-|", "$duckdb -csv -c \"$sql\"") or die "duckdb failed";
		while (my $line = <$fh>) {
		    chomp $line;
		    my ($project_id,$this_chr38,$this_pos38,$allele,$he,$ho) = split(/,/, $line);
	    	next if $project_id eq 'project';

	    	next if not $allele =~ /[ATGC]+/;

	    	my $rocksid = sprintf("%010d", $this_pos38).'!'.$allele;
	    	$h_rocks_to_view->{$this_chr38}->{$rocksid}->{$project_id} = undef;
	    	$h_rocks_to_view->{$this_chr38}->{$rocksid}->{he} += $he;
	    	$h_rocks_to_view->{$this_chr38}->{$rocksid}->{ho} += $ho;
	    	$i++;
	    	print '.' if ($i % 100000 == 0);
		}
		close($fh);
		print '.found.'.$i.'.';
		$i = 0;
		my $chr = $dejavu_variants->project->getChromosome($chr_filter);
		my $no = $chr->rocksdb('gnomad');
		
		my @l_var_chr = keys %{$h_rocks_to_view->{$chr_filter}};
		print ".begin_prepare_rocks.";
		$no->prepare(\@l_var_chr);
		print ".end_prepare_rocks.";
		
		foreach my $rocksid (@l_var_chr) {
			my $h_gad = $no->value($rocksid);
			if ($h_gad and exists $h_gad->{ac} and $h_gad->{ac} > $max_gnomad_ac) {
				delete $h_rocks_to_view->{$chr_filter}->{$rocksid};
				next;
			}
			if ($h_gad and exists $h_gad->{ho} and $h_gad->{ho} > $max_gnomad_ac_ho) {
				delete $h_rocks_to_view->{$chr_filter}->{$rocksid};
				next;
			}
			
			delete $h_rocks_to_view->{$chr_filter}->{$rocksid}->{he};
			delete $h_rocks_to_view->{$chr_filter}->{$rocksid}->{ho};
			$i++;
		}
		$no->close();
		print '.found_filtred.'.$i.'.';
	}
}


### PART 3 - check variants from calling variables

my ($hGenes, $hVariantsDetails) = $dejavu_variants->check_variants_from_gene($h_rocks_to_view);
print '...html...nbVar:'.scalar(keys %{$hVariantsDetails}).'.';
my $nb_genes = scalar(keys %{$hGenes});
print '.nbGenes:'.$nb_genes.'.';

#warn Dumper $hVariantsDetails;
#die;

### PART 4 - print HTML 

my ($h_html_genes, $h_phenos);
MCE::Loop->init(
   max_workers => $fork,
   chunk_size => 1,
   gather => sub {
        my ($data) = @_;
        foreach my $pheno_name (keys %{$data->{phenotypes}}) {
        	$h_phenos->{$pheno_name}->{tag} = $data->{phenotypes}->{$pheno_name}->{tag};
        }
        delete $data->{phenotypes};
        foreach my $score (sort {$b <=> $a} keys %$data) {
	        foreach my $gene_id (sort keys %{$data->{$score}}) {
	        	$h_html_genes->{$score}->{$gene_id} = $data->{$score}->{$gene_id} if $data->{$score}->{$gene_id};
	        }
        }
   }
);
mce_loop {
	my ($mce, $chunk_ref, $chunk_id) = @_;
	my $hres;
	foreach my $gene_id (@$chunk_ref) {
		next if $gene_id eq 'intronic';
		my @l_gene_id_tmp = split('_', $gene_id);
		print '.';
		my @list_variants = keys %{$hGenes->{$gene_id}};
		my ($this_html, $this_h_pheno, $max_score_gene);
		eval {
			($this_html, $this_h_pheno, $max_score_gene) = $dejavu_variants->print_html_gene($gene_id, \@list_variants, $hVariantsDetails);
			$hres->{$max_score_gene}->{$gene_id} = $this_html;
			foreach my $pheno_name (keys %$this_h_pheno) {
				$hres->{phenotypes}->{$pheno_name}->{tag} = $this_h_pheno->{$pheno_name};
			}
		};
		if ($@) {
			$hres->{$gene_id} = qq{ERROR with $gene_id};
		}
    }
	MCE->gather($hres);
} sort keys %{$hGenes};	
MCE::Loop->finish();

my @lPhenos = sort keys %$h_phenos;
my $html;

if ($dejavu_variants->alert_too_much_results()) {
	$html .= "<div style='width:100%;overflow-x:auto;'><table><tr>";
	$html .= "<td><b><i><span class='glyphicon glyphicon-alert' style='color:red'></span><span style='color:red;'> Too much results... partial results !</span></b></i>&nbsp;&nbsp;</td>";
	$html .= "</tr></table></div><br>"
}

if ($dejavu_variants->alert_ncboost_min_cadd_25()) {
	$html .= "<div style='width:100%;overflow-x:auto;'><table><tr>";
	$html .= "<td><b><i><span class='glyphicon glyphicon-alert' style='color:red'></span><span style='color:red;'> Only variants with cadd score >= 25 for ncboost filter !</span></b></i>&nbsp;&nbsp;</td>";
	$html .= "</tr></table></div><br>"
}

$html .= "<div style='width:100%;overflow-x:auto;'><table><tr>";
$html .= "<td><b>View phenotype</b>&nbsp;&nbsp;</td>";
my $cmd_all = qq{show_phenotype('');};
$html .= "<td><button type='button' class='btn btn-outline-primary' onClick=\"$cmd_all\" style='margin-right:5px;border: solid 0.5 black;font-size:12px;'><b><span style='color:green;'>All</span></b></button></td>";
if ($h_phenos) {
	foreach my $pheno (@lPhenos) {
		my $pheno_tag = $h_phenos->{$pheno}->{tag};
		my $cmd = qq{show_phenotype('$pheno_tag');};
		$html .= "<td><button type='button' class='btn btn-outline-primary' onClick=\"$cmd\" style='margin-right:5px;border: solid 0.5 black;font-size:12px;'>$pheno <i>(<b><pan id='span_nb_".$pheno."'>?</span></b>)</i></button></td>";
	}
}
$html .= "</tr></table></div><br>";
$html .= qq{<table id='table_genes' data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-virtual-scroll="true" data-pagination-v-align="both" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="50" data-page-list="[25, 50, 100, 200, 300]" data-resizable='true' class='table table-striped' style='font-size:13px;'>};
$html .= "<thead>";
$html .= $cgi->start_Tr({style=>"background-color:#E9DEFF;"});
$html .= qq{<th data-field="gene" data-filter-control="input" data-filter-control-placeholder="Gene name, description, ..."</th>};
$html .= $cgi->end_Tr();
$html .= "</thead>";
$html .= "<tbody>";
foreach my $score (sort {$b <=> $a} keys %$h_html_genes) {
	foreach my $gene_id (sort keys %{$h_html_genes->{$score}}) {
		my $html_gene = $h_html_genes->{$score}->{$gene_id};
		$html .= "<tr><td>".$html_gene."</td></tr>";
	}
}
$html .= "</tbody>";
$html .= "</table>";

if ($launch_job) {
	my $hRes;
	$hRes->{status} = "finished";
	$hRes->{html}   = $html;
	$hRes->{phenotypes} = join(',', @lPhenos);
	
	open(my $out, ">", $outfile);
	print $out encode_json($hRes);
	close $out;
}
else {
	my $hRes;
	$hRes->{html} = $html;
	my $json_encode = encode_json $hRes;
	print ".\",";
	$json_encode =~ s/{//;
	print $json_encode;
	exit(0);
}




sub launch_ncboost {
	my ($dejavu_variants, $ncboost_value, $max_dejavu, $max_dejavu_ho) = @_;
	$dejavu_variants->min_ncboost($ncboost_value);
	my $h_projects_parquet;
	if (exists $dejavu_variants->hash_users_projects->{all}) {
		$dejavu_variants->{alert_ncboost_min_cadd_25} = 1;
		$h_projects_parquet->{all} = $project->deja_vu_public_projects_parquet()."/NGS*.parquet";
	}
	else {
		foreach my $id (keys %{$dejavu_variants->hash_users_projects()}) {
			my $proj_name = $dejavu_variants->hash_users_projects->{$id}->{name};
			my $proj_id = $dejavu_variants->hash_users_projects->{$id}->{id};
			my $parquet = $project->deja_vu_public_projects_parquet()."/".$proj_name.".".$proj_id.".parquet";
			$h_projects_parquet->{$proj_id} = $parquet if -e $parquet;
		}
	}
	
	my @l_files;
	foreach my $file (values %{$h_projects_parquet}) {
		push(@l_files, $file);
	}
	my @lChr = (1..22, 'X', 'Y');
	MCE::Loop->init(
		max_workers => 1,
		chunk_size => '1',
		gather => sub {
	        my ($data) = @_;
	        print '|';
	        foreach my $id (keys %{$data->{h_res_duck}}) {
	        	$h_res_duck_ncboost->{$id} = $data->{h_res_duck}->{$id};
	        }
	        foreach my $id (keys %{$data->{h_ncboost_values}}) {
	        	$h_ncboost_values->{$id} = $data->{h_ncboost_values}->{$id};
	        }
	    }
	);
	mce_loop {
		my ($mce, $chunk_ref, $chunk_id) = @_;
		my $hres;
		
		my @list_sql_annot;
		foreach my $annot (sort keys %{$dejavu_variants->{hash_filters_cons}}) {
			push(@list_sql_annot, "(annotation='".lc($annot)."')");
		}
		my $sql_annot = "(".join(' OR ', @list_sql_annot).")";
		
		my ($chr_filter, $start_filter, $end_filter) = split('-', $region) if $region;
		foreach my $chr_id (@$chunk_ref) {
			next if $chr_filter and $chr_filter ne $chr_id;
			my $sql_region_end;
			if ($start_filter and $end_filter) {
				$sql_region_end = "AND a.pos >= $start_filter AND a.pos <= $end_filter";
			}
			my $sql;
			if ($dejavu_variants->is_magic_user()) {
				$sql = "
					PRAGMA threads=$fork;
					SELECT 
					    a.pos, a.rocksdb_id, a.score AS ncboost
					FROM read_parquet('/data-isilon/public-data/repository/HG38/ncboost/20260301/parquet/dejavu_filter/chr=$chr_id/*.parquet') a
					WHERE a.score >= $ncboost_value 
					  AND dejavu <= $max_dejavu 
					  AND dejavu_ho <= $max_dejavu_ho
					  AND a.gnomad_ac <= $max_gnomad_ac
					  AND a.gnomad_ho <= $max_gnomad_ac_ho
					  AND $sql_annot
					  $sql_region_end
					  ;
				 ";
			}
			else {
				my $sql_parquets = $dejavu_variants->sql_projects_parquet();
				$sql = "
					PRAGMA threads=$fork;
	
					WITH b_filtered AS (
					    SELECT pos38
					    FROM $sql_parquets
					    WHERE chr38 = '$chr_id'
					      AND max_ratio > 0
					      AND allele IN ('A','T','C','G')
					    GROUP BY pos38
					)
					
					SELECT 
					    a.pos, a.rocksdb_id, a.score AS ncboost
					FROM read_parquet('/data-isilon/public-data/repository/HG38/ncboost/20260301/parquet/dejavu_filter/chr=$chr_id/*.parquet') a
					JOIN b_filtered b
					    ON a.pos = b.pos38
					WHERE a.score >= $ncboost_value 
					  AND dejavu <= $max_dejavu 
					  AND dejavu_ho <= $max_dejavu_ho
					  AND a.gnomad_ac <= $max_gnomad_ac
					  AND a.gnomad_ho <= $max_gnomad_ac_ho
					  AND $sql_annot
					  $sql_region_end
					  ;
				";
			}
			
			
##			if (exists $dejavu_variants->hash_users_projects->{all}) {
##				$sql .= "SELECT pos, rocksdb_id, score AS ncboost 
##						FROM read_parquet('$path_ncboost_dejavu/chr=$chr_id/data_0.parquet')
##						WHERE
##							score >= $ncboost_value 
##							AND dejavu <= $max_dejavu 
##							AND dejavu_ho <= $max_dejavu_ho
##							AND gnomad_ac <= $max_gnomad_ac
##							AND gnomad_ho <= $max_gnomad_ac_ho
##							AND $sql_annot
##							AND (cadd >= 25);";
##			}
##			else {
#				$sql .= "WITH b_filtered AS (";
#				my @lproj_sql;
#				my $zz = 0;
#				foreach my $file (@l_files) {
#					
#					my $sql_part = "SELECT project,chr38,chr19,pos38,pos19,he,allele,patients,dp_ratios FROM read_parquet('$file')
#						WHERE
#							chr38='$chr_id'
#							and max_ratio > 0
#							and (allele='A' or allele='T' or allele='C' or allele='G')";
#					push(@lproj_sql, $sql_part);
#				}
#				
#				$sql .= join(' UNION ALL ', @lproj_sql);
#				$sql .= ")";
#				$sql .= "SELECT 
#							a.pos, a.rocksdb_id, a.score AS ncboost 
#							FROM read_parquet('$path_ncboost_dejavu/chr=$chr_id/data_0.parquet') a
#							WHERE
#								a.score >= $ncboost_value 
#								AND dejavu <= $max_dejavu 
#								AND dejavu_ho <= $max_dejavu_ho
#								AND a.gnomad_ac <= $max_gnomad_ac
#								AND a.gnomad_ho <= $max_gnomad_ac_ho
#								AND $sql_annot
#								AND a.pos IN (SELECT pos38 FROM b_filtered)
#								$sql_region_end;";
##			}
			
			my $i = 0;
			my $duckdb = $buffer->software('duckdb');
			open(my $fh, "-|", "$duckdb -csv -c \"$sql\"") or die "duckdb failed";
			while (my $line = <$fh>) {
			    chomp $line;
			    my ($pos, $id, $ncboost) = split(',', $line);
			    next if $id eq 'rocksdb_id';
				$i++;
				if ($i == 20000) {
					print '.';
					$i = 0;
				}
				$hres->{h_res_duck}->{$id}->{geneid} = $line;
				$hres->{h_ncboost_values}->{$id} = $ncboost;
			}
			close($fh);
		}
		MCE->gather($hres);
	} @lChr;
}

sub launch_promoter_ai {
	my ($dejavu_variants, $promoter_ai_value) = @_;
	my ($chr_filter, $start_filter, $end_filter) = split('-', $region) if $region;
	$dejavu_variants->min_promoter_ai($promoter_ai_value);
	my $sql_promoter_ai = "PRAGMA threads=$fork; SELECT rocksdb_id, geneid FROM read_parquet(['$parquet_promoter_ai_filtred']) WHERE ABS(promoterAI) >= $promoter_ai_value";
	my $i = 0;
	my $duckdb = $buffer->software('duckdb');
	open(my $fh, "-|", "$duckdb -csv -c \"$sql_promoter_ai\"") or die "duckdb failed";
	while (my $line = <$fh>) {
	    chomp $line;
	    my ($id, $geneid) = split(/,/, $line);
	    next if $id eq 'rocksdb_id';
	    if ($region) {
	    	my @ltmp = split('!', $id);
	    	next if $chr_filter ne $ltmp[0];
	    	if ($start_filter and $end_filter) {
	    		next if int($ltmp[1]) < $start_filter;
	    		next if int($ltmp[1]) > $end_filter;
	    	}
	    }
		$h_res_duck->{$id} = $geneid;
		$h_genes_only->{$geneid} = undef;
		$i++;
		if ($i == 100000) {
			print '.';
			$i = 0;
		}
	}
	close($fh);
}

sub get_hash_annot_categories {
	my $h_annot_categories;
	foreach my $cat_name (keys %{$buffer->config->{ensembl_annotations}}) {
		$h_annot_categories->{lc($cat_name)} = $cat_name;
		my @lOthersNames = split(';', $buffer->config->{ensembl_annotations}->{$cat_name});
		foreach my $other_name (@lOthersNames) {
			$other_name =~ s/ /_/g;
			$h_annot_categories->{lc($other_name)} = $cat_name;
		}
	}
	return $h_annot_categories;
}

