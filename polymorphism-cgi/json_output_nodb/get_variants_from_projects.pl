#!/usr/bin/env perl
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
use xls_export;


## LIMIT RAM utilisee - 20go
#use BSD::Resource; 
#my $limit_ram = 20000 * 1024 * 1024;
#setrlimit(RLIMIT_AS, $limit_ram, $limit_ram) or die "\n\nTmpossible fixer limite RAM: $!";

my $cgi = new CGI();

my $fork = 8;
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


my $login = $cgi->param('login');
my $pwd   = $cgi->param('pwd');
my $load_html = $cgi->param('load_html');
my $promoter_ai_value = $cgi->param('promoter_ai_value');
my $ncboost_value = $cgi->param('ncboost_value');
my $max_gnomad_ac = $cgi->param('gnomad');
my $max_gnomad_ac_ho = $cgi->param('gnomad_ho');
my $max_dejavu = $cgi->param('dejavu');
my $max_dejavu_ho = $cgi->param('dejavu_ho');
my $min_ratio = $cgi->param('min_perc_all');
my $min_reads = $cgi->param('min_reads');
my $only_ill = $cgi->param('only_ill');
my $only_strict_ill = $cgi->param('only_strict_ill');
my $models = $cgi->param('models');
my $region = $cgi->param('region');
my $keep_pathogenic = $cgi->param('keep_pathogenic');
my $only_genes = $cgi->param('only_genes');
my $use_phenotype = $cgi->param('use_phenotype');
my $exclude_projects = $cgi->param('exclude_projects');
my $only_project = $cgi->param('only_project');
my $only_patient = $cgi->param('only_patient');
my $view_all_projects = $cgi->param('view_all_projects');
my $export_xls = $cgi->param('export_xls');
my $session_id = $cgi->param('session_id');

my $dejavu_variants = new dejavu_variants();
$dejavu_variants->view_all_projects(1) if $view_all_projects;
$dejavu_variants->user_name($login);
$dejavu_variants->pwd($pwd);


if ($load_html) {
	my $buffer = new GBuffer;
	print $cgi->header(-type => 'text/html; charset=utf-8');
	my $file = $buffer->config_path("root","global_search").'/'.$load_html.'.html';
	open my $fh, '<', $file or die $!;
	print while (<$fh>);
	close $fh;
	unlink($file);
	exit(0);
}

if ($export_xls) {
	export_xls($dejavu_variants, $session_id);
	exit(0);
}

print $cgi->header('text/json-comment-filtered');
print "{\"progress\":\".";

if ($exclude_projects) {
	my @lExcl = split(',', $exclude_projects);
	foreach my $p_name (@lExcl) {
		my $project_name = 'NGS20'.$p_name;
		$dejavu_variants->{hash_exclude_projects}->{$project_name} = undef;
	}
}
if ($only_project) {
	my $patient_name = undef;
	$patient_name = $only_patient if $only_patient;
	$dejavu_variants->{hash_only_project}->{$only_project} = $patient_name;
}

exit(0) if not $dejavu_variants->hash_users_projects();
print '.nb_proj.'.scalar(keys %{$dejavu_variants->hash_users_projects()});

$dejavu_variants->fork($fork);

$dejavu_variants->use_phenotype($use_phenotype) if $use_phenotype;
$dejavu_variants->max_dejavu($max_dejavu) if $max_dejavu;
$dejavu_variants->max_dejavu_ho($max_dejavu_ho) if $max_dejavu_ho;
$dejavu_variants->max_gnomad_ac($max_gnomad_ac) if $max_gnomad_ac;
$dejavu_variants->max_gnomad_ac_ho($max_gnomad_ac_ho) if $max_gnomad_ac_ho;
$dejavu_variants->min_ratio($min_ratio) if $min_ratio;
$dejavu_variants->only_ill_patients(1) if $only_ill;
$dejavu_variants->only_strict_ill_patients(1) if $only_strict_ill;
$dejavu_variants->min_reads($min_reads) if $min_reads;
$dejavu_variants->only_genes(1) if $only_genes;

my $h_models;
if ($models) {
	foreach my $model_name (split(',', $models)) {
		$h_models->{$model_name} = 1;
	}
	$dejavu_variants->{models} = $h_models
}

my $buffer = new GBuffer;
my $project_name = $buffer->getRandomProjectName('HG38_DRAGEN', '46', '21');
my $project = $buffer->newProject( -name => $project_name );

my $h_errors_found;
my $h_only_genes;
if ($only_genes) {
	$only_genes =~ s/ //g;
	$only_genes =~ s/[,]+/,/g;
	foreach my $gene_id (split(',', $only_genes)) {
		next if not $gene_id;
		my @tests;
		push(@tests, $gene_id);
		push(@tests, uc($gene_id));
		push(@tests, lc($gene_id));
		my $found;
		foreach my $gene_id_test (@tests) {
			next if $found;
			my $gene;
			eval {
				$gene = $project->newGene($gene_id_test);
				$h_only_genes->{$gene->id()} = $gene->getChromosome->id().'-'.($gene->start - 500).'-'.($gene->end + 500);
				$h_only_genes->{$gene->external_name()} = $gene->getChromosome->id().'-'.($gene->start - 500).'-'.($gene->end + 500);
				$found = 1;
			};
			if($@) {
				$found = undef;
			}
		}
		if (not $found) {
			$h_errors_found->{$gene_id} = "gene not found in gencode ".$project->gencode_version();
		}
	}
	$dejavu_variants->{hash_only_genes} = $h_only_genes;
}

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

my ($h_res_duck, $h_res_duck_only_pos, $h_genes_only, $h_ncboost_values);
my ($h_res_duck_promoter_ai, $h_genes_only_promoter_ai);
my ($h_res_duck_ncboost);

### PART 1 - Select variants with promoterAI / NCBOOST / REGION


my $h_rocks_to_view;

if ($promoter_ai_value) {
	launch_promoter_ai ($dejavu_variants, $promoter_ai_value);
	print '.promoter_ai.'.scalar(keys %{$h_res_duck});
	my $h_vector;
	my $n = 0;
	print '.intspan.'.scalar(keys %$h_res_duck).'.';
	my $ok = 0;
	foreach my $rocksid (keys %$h_res_duck) {
		$n++;
		print '.' if $n % 10000 == 0;
		my @l_tmp = split('!', $rocksid);
		my $chr_id = $l_tmp[0];
		my $pos = int($l_tmp[1]);
		my $dv_rocks_id = $l_tmp[1].'!'.$l_tmp[2];
		my $h_dv = pass_frequence_or_clinvar_pathogenic($chr_id, $dv_rocks_id);
		next if not $h_dv;
		$h_rocks_to_view->{$chr_id}->{$dv_rocks_id} = $h_dv->{details};
		$ok++;
	}
	print '.after_dv.'.$ok++.'.';
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
		$h_rocks_to_view->{$ltmp[0]}->{$ltmp[1].'!'.$ltmp[2]} = $h_proj;
		$found++;
	}
	print '.now.'.$found++.'.';
}

if ($region or $only_genes) {
	if (not $promoter_ai_value and not $ncboost_value) {
		
		my @l_regions_tmp;
		if ($region) { push(@l_regions_tmp, $region); }
		if ($only_genes) {
			foreach my $gid (keys %{$h_only_genes}) {
				push(@l_regions_tmp, $h_only_genes->{$gid});
			}
		}
		
		my $h_regions;
		foreach my $this_region (@l_regions_tmp) {
			my ($chr_id, $start, $end) = split('-', $this_region);
			if (not exists $h_regions->{$chr_id}) { $h_regions->{$chr_id} = Set::IntSpan::Fast::XS->new(); }
			$h_regions->{$chr_id}->add_range($start, $end) if $start and $end;
		}
		
		my @l_regions;
		foreach my $chr_id (keys %$h_regions) {
			my @l_this_region;
			push(@l_this_region, $chr_id);
			my $string = $h_regions->{$chr_id}->as_string();
			foreach my $interval (split(',', $string)) {
				my ($start, $end) = split('-', $interval);
				push(@l_this_region, $start.'-'.$end);
			}
			push(@l_regions, \@l_this_region);
		}
		
		my ($fork2, $fork_sql);
		if ($dejavu_variants->is_magic_user() or scalar(@l_regions) == 1) {
			$fork2 = 1;
			$fork_sql = $dejavu_variants->fork();
		}
		else {
			$fork2 = scalar(@l_regions);
			$fork2 = 10 if $fork2 >= 10;
			$fork_sql = 1;
		}
		MCE::Loop->init(
		    max_workers => $fork2,
		    chunk_size  => 'auto',
		    gather      => sub {
		        my ($data) = @_;
		        return unless $data;
		        foreach my $this_chr38 (keys %{$data->{h_rocks_to_view}}) {
		        	#$h_rocks_to_view
		        	foreach my $rocksid (keys %{$data->{h_rocks_to_view}->{$this_chr38}}) {
		        		foreach my $cat (keys %{$data->{h_rocks_to_view}->{$this_chr38}->{$rocksid}}) {
		        			if ($cat eq 'he') { $h_rocks_to_view->{$this_chr38}->{$rocksid}->{he} += $data->{h_rocks_to_view}->{$this_chr38}->{$rocksid}->{he}; }
		        			elsif ($cat eq 'ho') { $h_rocks_to_view->{$this_chr38}->{$rocksid}->{ho} += $data->{h_rocks_to_view}->{$this_chr38}->{$rocksid}->{ho}; } 
		        			else {
		        				my $project_id = $cat;
		        				$h_rocks_to_view->{$this_chr38}->{$rocksid}->{$project_id} = undef;
		        			}
		        		}
		        	}
		        }
		        if (exists $data->{pathogenic}) {
			        foreach my $rocksid (keys %{$data->{pathogenic}}) {
			        	$dejavu_variants->{hash_variant_pathogenic}->{$rocksid} = undef;
			        }
		        }
		    }
		);
		
		my $worker = sub {
			my ($mce, $chunk_ref, $chunk_id) = @_;
			my $hres;
			
			foreach my $list_this_region (@$chunk_ref) {
				my $chr_filter = shift(@$list_this_region);
				print '.only_region.chr'.$chr_filter.'.';
				$dejavu_variants->{only_chromosome} = $chr_filter;
				my ($sql_pos, @list_sql_pos, $found_pos);
				foreach my $interval (@$list_this_region) {
					my ($start_filter, $end_filter) = split('-', $interval);
					$end_filter = $start_filter +1 if not $end_filter;
					push(@list_sql_pos, "(pos38 >= $start_filter and pos38 <= $end_filter)");
					$found_pos++;
				}
				if ($found_pos) {
					$sql_pos = 'and ('.join(' OR ', @list_sql_pos).')';
				}
				
				my $i = 0;
				my $sql_parquets = $dejavu_variants->sql_projects_parquet();
				
				my $this_fork = $dejavu_variants->fork();
				my $sql = "
					PRAGMA threads=$fork_sql;
					WITH base AS ( SELECT project, chr38, pos38, allele, he, ho FROM $sql_parquets WHERE chr38='$chr_filter' $sql_pos ),
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
			    	my $rocksid = sprintf("%010d", $this_pos38).'!'.$allele;
			    	$hres->{h_rocks_to_view}->{$this_chr38}->{$rocksid}->{$project_id} = undef;
			    	$hres->{h_rocks_to_view}->{$this_chr38}->{$rocksid}->{he} += $he;
			    	$hres->{h_rocks_to_view}->{$this_chr38}->{$rocksid}->{ho} += $ho;
			    	$i++;
			    	print '.' if ($i % 100000 == 0);
				}
				close($fh);
				
				$i = 0;
				my $chr = $dejavu_variants->project->getChromosome($chr_filter);
				my $no = $chr->rocksdb('gnomad');
				
				my @l_var_chr = keys %{$hres->{h_rocks_to_view}->{$chr_filter}};
				foreach my $rocksid (@l_var_chr) {
					eval {
						my $h_gad = $no->value($rocksid);
						if ($h_gad and exists $h_gad->{ac} and $h_gad->{ac} > $max_gnomad_ac) {
							delete $hres->{h_rocks_to_view}->{$chr_filter}->{$rocksid};
							next;
						}
						if ($h_gad and exists $h_gad->{ho} and $h_gad->{ho} > $max_gnomad_ac_ho) {
							delete $hres->{h_rocks_to_view}->{$chr_filter}->{$rocksid};
							next;
						}				
						delete $hres->{h_rocks_to_view}->{$chr_filter}->{$rocksid}->{he};
						delete $hres->{h_rocks_to_view}->{$chr_filter}->{$rocksid}->{ho};
						$i++;
					};
					if ($@) {
						$i++;
						next;
					}
				}
				$no->close();
				
				#ADD CLINVAR PATHOGENIC
				if ($keep_pathogenic) {
					my $sql_clinvar = "
						PRAGMA threads=$fork_sql;
						WITH dejavu AS (
							SELECT *, chr38 || '!' || LPAD(CAST(pos38 AS VARCHAR), 10, '0') || '!' || allele AS 'index', LPAD(CAST(pos38 AS VARCHAR), 10, '0') || '!' || allele AS 'rocksid' FROM $sql_parquets WHERE chr38='$chr_filter' $sql_pos
						)
						SELECT d.project, c.index, d.rocksid
							FROM dejavu d
							JOIN '/data-pure/public-data/repository/HG38/clinvar/20250824/parquet/clinvar.csv' c ON d.index = c.index
							WHERE c.clinvar_class = 'pathogenic'
							GROUP BY d.project, c.index, d.rocksid;
					";
					
					open(my $fh2, "-|", "$duckdb -csv -c \"$sql_clinvar\"") or die "duckdb failed";
					while (my $line = <$fh2>) {
					    chomp $line;
					    my ($project_id, $genomic_rocksid, $rocksid) = split(/,/, $line);
				    	next if $project_id eq 'project';
				    	my $res = $project->getChromosome($chr_filter)->rocks_dejavu->dejavu($rocksid);
						my ($nb_pat_he, $nb_pat_ho) = get_dv_he_ho_from_request($res);
				    	$hres->{h_rocks_to_view}->{$chr_filter}->{$rocksid}->{$project_id} = undef;
				    	$hres->{h_rocks_to_view}->{$chr_filter}->{$rocksid}->{he} = $nb_pat_he;
				    	$hres->{h_rocks_to_view}->{$chr_filter}->{$rocksid}->{ho} = $nb_pat_ho;
				    	
				    	$hres->{pathogenic}->{$rocksid} = undef;
				    	$i++;
				    	print '.' if ($i % 100000 == 0);
					}
					close($fh2);
#					print '.with_clinvar.'.$i.'.';
				}
			}
		    MCE->gather($hres);
		};
		MCE::Loop->run($worker, [ @l_regions ]);
		MCE::Loop->finish();
	}
}


### PART 3 - check variants from calling variables
my ($hGenes, $hVariantsDetails) = $dejavu_variants->check_variants_from_gene($h_rocks_to_view);

my $nb_var = scalar(keys %{$hVariantsDetails});
print '...html...nbVar:'.$nb_var.'.';
my $nb_genes = scalar(keys %{$hGenes});
if ($only_genes) {
	foreach my $gene_id (keys %$hGenes)	{
		delete $hGenes->{$gene_id} if not exists $dejavu_variants->{hash_only_genes}->{$gene_id};
	}
}
$nb_genes = scalar(keys %{$hGenes});
print '.nbGenes:'.$nb_genes.'.';

### PART 4 - print HTML  GENES
my ($h_genes_out_html, $h_phenos);
MCE::Loop->init(
   max_workers => $fork,
   chunk_size => 'auto',
   gather => sub {
        my ($data) = @_;
        foreach my $gene_id (keys %{$data}) {
        	$h_genes_out_html->{$gene_id} = $data->{$gene_id};
        }
   }
);
mce_loop {
	my ($mce, $chunk_ref, $chunk_id) = @_;
	my $hres;
	foreach my $gene_id (@$chunk_ref) {
		my @list_variants = keys %{$hGenes->{$gene_id}};
		$hres->{$gene_id} = $dejavu_variants->print_html_gene_only($gene_id, \@list_variants, $hVariantsDetails);
	}
	MCE->gather($hres);
} sort keys %{$hGenes};	
MCE::Loop->finish();

#warn Dumper $h_genes_out_html;
#die;

#Le HTML est dans $hres->{$gene_id}->{html} !!

### PART 4 - print HTML VARIANTS
my $nb_var_todo = 0;
my @l_var_todo;
foreach my $gene_id (sort keys %{$hGenes}) {
	foreach my $var_id (keys %{$hGenes->{$gene_id}}) {
		my $h;
		next if not exists $hVariantsDetails->{$var_id};
		next if not $hVariantsDetails->{$var_id}->{polyviewer_html_details_proj_pat};
		$h->{var_id} = $var_id;
		$h->{h_gene} = $h_genes_out_html->{$gene_id};
		push(@l_var_todo, $h);
		$nb_var_todo++;
	}
}
print '.nbVar:'.$nb_var_todo.'.';

#list_html_variants
MCE::Loop->init(
   max_workers => $fork,
   chunk_size => 'auto',
   gather => sub {
        my ($data) = @_;
        print '|';
        foreach my $gene_id (keys %{$data->{genes}}) {
        	foreach my $html_var (@{$data->{genes}->{$gene_id}}) {
        		push(@{$h_genes_out_html->{$gene_id}->{list_html_variants}}, $html_var);
        	}
        }
        foreach my $pheno_name (keys %{$data->{phenotypes}}) {
        	my $name = $data->{$pheno_name}->{name};
        	$h_phenos->{$pheno_name}->{name} = $name;
        	foreach my $gene_id (keys %{$data->{$pheno_name}->{genes}}) {
        		$h_phenos->{$pheno_name}->{genes}->{$gene_id} = undef;
        	}
        }
   }
);
mce_loop {
	my ($mce, $chunk_ref, $chunk_id) = @_;
	print '.';
	my $hres;
	foreach my $h (@$chunk_ref) {
		my $gene_id = $h->{h_gene}->{id};
		my ($this_out, $this_h_pheno) = $dejavu_variants->print_html_variant_only($h->{var_id}, $h->{h_gene}, $hVariantsDetails);
		push(@{$hres->{genes}->{$gene_id}}, $this_out);
		foreach my $pheno (keys %$this_h_pheno) {
			my $name = $this_h_pheno->{$pheno};
			$hres->{phenotypes}->{$pheno}->{name} = $name;
			$hres->{phenotypes}->{$pheno}->{genes}->{$gene_id} = undef;
		}
	}
	MCE->gather($hres);
} @l_var_todo;	
MCE::Loop->finish();



#warn Dumper $h_genes_out_html;
my $h_html_genes;

foreach my $gene_id (keys %$h_genes_out_html) {
	next if not exists $h_genes_out_html->{$gene_id}->{list_html_variants};
	
	#TODO: recuperer le bg_color
	my $bg_color = "background-color:#607D8B";
	my $nb_ok = scalar(@{$h_genes_out_html->{$gene_id}->{list_html_variants}});
	my $panel_id = "panel_".$gene_id;
	$panel_id .= "_he_comp" if $dejavu_variants->only_genes() and $dejavu_variants->hash_only_project() and $dejavu_variants->hash_only_patients();
	my $out = qq{<div class="panel-heading panel-face panel-grey"	style="$bg_color;height:43px;padding:10px;border:0px;width:100%;">};
	$out .= $h_genes_out_html->{$gene_id}->{html};
	my $table_id = 'table_'.$panel_id;
	$out .= qq{</div>};
	$out .= qq{<div style="height:3px;" loading="lazy"></div>};
	$out .= qq{<div loading="lazy" class="panel-body panel-collapse collapse" style="font-size: 09px;font-family:Verdana;" id="$panel_id" loading="lazy">};
	if ($dejavu_variants->only_genes() and $dejavu_variants->hash_only_project() and $dejavu_variants->hash_only_patients() and $nb_ok >= 5) {
		$out .= qq{<table loading="lazy" id='$table_id' data-filter-control='true' data-virtual-scroll="true" data-toggle="table" data-show-extended-pagination="false" data-cache="false" data-pagination-loop="true" data-pagination-v-align="both" data-pagination-h-align="left" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="5" data-page-list="[5, 10, 20]" data-resizable='true' class='table table-striped table-borderless table-bootstraptable' style='font-size:9px;'>};
	}
	elsif ($nb_ok >= 20) {
		$out .= qq{<table loading="lazy" id='$table_id' data-filter-control='true' data-virtual-scroll="true" data-toggle="table" data-show-extended-pagination="false" data-cache="false" data-pagination-loop="true" data-pagination-v-align="top" data-pagination-h-align="left" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="10" data-page-list="[5, 10, 20]" data-resizable='true' class='table table-striped table-borderless table-bootstraptable' style='font-size:9px;'>};
	}
	else {
		$out .= qq{<table class="table table-striped table-borderless" style="vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;line-height: 25px;min-height: 25px;height: 25px;box-shadow: 3px 3px 5px #555;">};
	}
	$out .= "<thead><tr style='background-color:aliceblue;color:black'>";
	$out .= qq{<th data-field="1"></th>};
	$out .= qq{<th data-field="2" data-sortable="true" data-filter-control="input" data-filter-control-placeholder="">Position</th>};
	$out .= qq{<th data-field="3" data-filter-control="input" data-filter-control-placeholder="">Project / Patient</th>};
	$out .= qq{<th data-field="4">gnomAD</th>};
	$out .= qq{<th data-field="5">DejaVu</th>};
	$out .= qq{<th data-field="6" data-filter-control="input" data-filter-control-placeholder="">Clinvar</th>};
	$out .= qq{<th data-field="7" data-filter-control="input" data-filter-control-placeholder="">Annotation</th>};
	$out .= "</tr></thead>";
	$out .= "<tbody>";
	
	$out .= join(' ', @{$h_genes_out_html->{$gene_id}->{list_html_variants}});
	
	$out .= "</tbody>";
	$out .= qq{</table>};
	$out .= qq{</div>};
	
	my $score = $h_genes_out_html->{$gene_id}->{max_score};
	$h_html_genes->{$score}->{$gene_id} = $out;
}




#my ($h_html_genes, $h_phenos);
#MCE::Loop->init(
#   max_workers => $fork,
#   chunk_size => 'auto',
#   gather => sub {
#        my ($data) = @_;
#        foreach my $pheno_name (keys %{$data->{phenotypes}}) {
#        	$h_phenos->{$pheno_name}->{tag} = $data->{phenotypes}->{$pheno_name}->{tag};
#        }
#        delete $data->{phenotypes};
#        foreach my $score (sort {$b <=> $a} keys %$data) {
#	        foreach my $gene_id (sort keys %{$data->{$score}}) {
#	        	$h_html_genes->{$score}->{$gene_id} = $data->{$score}->{$gene_id} if $data->{$score}->{$gene_id};
#	        }
#        }
#   }
#);
#mce_loop {
#	my ($mce, $chunk_ref, $chunk_id) = @_;
#	my $hres;
#	foreach my $gene_id (@$chunk_ref) {
#		next if $gene_id eq 'intronic';
#		if ($only_genes) {
#			next if not exists $h_only_genes->{$gene_id};
#		}
#		
#		my @l_gene_id_tmp = split('_', $gene_id);
#		my @list_variants = keys %{$hGenes->{$gene_id}};
#		print '.'.scalar(@list_variants).'.';
#		my ($this_html, $this_h_pheno, $max_score_gene);
#		eval {
#			($this_html, $this_h_pheno, $max_score_gene) = $dejavu_variants->print_html_gene($gene_id, \@list_variants, $hVariantsDetails);
#			$hres->{$max_score_gene}->{$gene_id} = $this_html;
#			foreach my $pheno_name (keys %$this_h_pheno) {
#				$hres->{phenotypes}->{$pheno_name}->{tag} = $this_h_pheno->{$pheno_name};
#			}
#		};
#		if ($@) {
#			$hres->{'-99'}->{$gene_id} = qq{ERROR with $gene_id};
#		}
#    }
#	MCE->gather($hres);
#} sort keys %{$hGenes};	
#MCE::Loop->finish();


my $export_session_id;
if ($nb_genes > 0 and $nb_var > 0) {
	my $hash_session_global; 
	print '.xls.prepare.';
	($hash_session_global->{hash_variants_global}, $hash_session_global->{hash_specific_infos}->{datas}->{projects_patients_infos}) = return_hashes_variants_in_session_export($dejavu_variants, $hGenes, $hVariantsDetails);
	print '.session.'; 
	$export_session_id = $dejavu_variants->xls_export_session->save_with_storable($hash_session_global);
	undef $hash_session_global;
	undef $hVariantsDetails;
	undef $hGenes;
	print '.saved';
}
else {
	$export_session_id = 'no_results_'.time;
}
my $gencode_version = $project->gencode_version();
my ($gencode, $annot_version) = split('\.', $project->annotation_version());

my @lPhenos = sort keys %$h_phenos;

my $html_file = $buffer->config_path("root","global_search").'/'.$export_session_id.'.html';
my $html;

#if ($dejavu_variants->alert_too_much_results()) {
#	$html .= "<div style='width:100%;overflow-x:auto;background-color:white !important;'><table><tr>";
#	$html .= "<td><b><i><span class='glyphicon glyphicon-alert' style='color:red'></span><span style='color:red;'> Too much results... partial results !</span></b></i>&nbsp;&nbsp;</td>";
#	$html .= "</tr></table></div><br>"
#}
#
#if ($dejavu_variants->alert_ncboost_min_cadd_25()) {
#	$html .= "<div style='width:100%;overflow-x:auto;'><table><tr>";
#	$html .= "<td><b><i><span class='glyphicon glyphicon-alert' style='color:red'></span><span style='color:red;'> Only variants with cadd score >= 25 for ncboost filter !</span></b></i>&nbsp;&nbsp;</td>";
#	$html .= "</tr></table></div><br>"
#}
if ($h_errors_found) {
	$html .= "<div style='width:100%;overflow-x:auto;'><table><tr>";
	$html .= "<td><b><pan style='color:red;'>ERRORS: </span></b>&nbsp;&nbsp;</td>";
	foreach my $val (keys %$h_errors_found) {
		my $msg = $val.': '.$h_errors_found->{$val};
		$html .= "<td><button type='button' class='btn btn-outline-primary' style='color:white;background-color:red;margin-right:5px;border: solid 0.5 black;font-size:12px;'>$msg</button></td>";
	}
	$html .= "</tr></table></div><br>";
}


if ($dejavu_variants->only_genes() and $dejavu_variants->hash_only_project() and $dejavu_variants->hash_only_patients()) { $html .= "<br>"; }
else {
	$html .= "<div style='width:100%;overflow-x:auto;'><table><tr>";
	$html .= "<td><b><nobr>Score Phenotype</nobr></b>&nbsp;&nbsp;</td>";
	$html .= "<td><button type='button' class='btn btn-outline-danger' style='margin-right:5px;border: solid 0.5 black;font-size:12px;'><b><span style='color:red;'><nobr>$use_phenotype</nobr></span></b></button></td>";
	
	$html .= "<td><b>Gencode</b>&nbsp;&nbsp;</td>";
	$html .= "<td><button type='button' class='btn btn-outline-success' style='margin-right:5px;border: solid 0.5 black;font-size:12px;'><b><span style='color:green;'>$gencode_version</span></b></button></td>";
	foreach my $cat_name (sort keys %{$buffer->public_data->{$annot_version}}) {
		next if lc($cat_name) eq 'polyscore';
		next if $cat_name eq 'cldb';
		next if $cat_name eq 'dbscsnv_rf';
		next if $cat_name eq 'cytoband';
		next if $cat_name eq 'hgmd';
		next if $cat_name =~ /gnomad-.*/;
		next if $cat_name =~ /-hg19/;
		my $cat_version = $buffer->public_data->{$annot_version}->{$cat_name}->{version};
		$html .= "<td>&nbsp;&nbsp;<b>$cat_name</b>&nbsp;&nbsp;</td>";
		$html .= "<td><button type='button' class='btn btn-outline-primary' style='margin-right:5px;border: solid 0.5 black;font-size:12px;'><b>$cat_version</b></button></td>";
	}
	$html .= "</tr></table></div>";
}
undef $project;
undef $buffer;
undef $dejavu_variants;

$html .= "<div style='width:100%;margin-top:10px;overflow-x:auto;'><table><tr>";
$html .= "<td><b>Export Results</b>&nbsp;&nbsp;</td>";
my $cmd_export_xls = qq{export_xls('$export_session_id');};
$html .= "<td><button type='button' class='btn btn-outline-success' onClick=\"$cmd_export_xls\" style='margin-right:5px;border: solid 0.5 black;font-size:12px;'><b><span style='color:green;'>XLS</span></b></button></td>";
$html .= "<td><b>&nbsp;&nbsp;View phenotype</b>&nbsp;&nbsp;</td>";
my $cmd_all = qq{show_phenotype('');};
$html .= "<td><button type='button' class='btn btn-outline-success' onClick=\"$cmd_all\" style='margin-right:5px;border: solid 0.5 black;font-size:12px;'><b><span style='color:green;'>All</span></b></button></td>";
if ($h_phenos) {
	foreach my $pheno (@lPhenos) {
		my $pheno_tag = $h_phenos->{$pheno}->{name};
		my $cmd = qq{show_phenotype('$pheno_tag');};
		$html .= "<td><button type='button' class='btn btn-outline-primary' onClick=\"$cmd\" style='margin-right:5px;border: solid 0.5 black;font-size:12px;'>$pheno <i>(<b><span id='span_nb_".$pheno."'>?</span></b>)</i></button></td>";
	}
}
$html .= "</tr></table></div>";
undef $h_phenos; 


my $nb_genes = scalar keys %$h_html_genes;
my $data_search = 'false';
$data_search = 'true' if $nb_genes <= 50;

$html .= qq{<table id='table_genes' data-filter-control='$data_search' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-virtual-scroll="true" data-pagination-v-align="both" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="20" data-page-list="[20, 50, 100, 200, 300]" data-resizable='true' class='table' style='font-size:13px;'>};
$html .= "<thead>";
$html .= $cgi->start_Tr({style=>"background-color:#E9DEFF;"});
if ($nb_genes > 50) { $html .= qq{<th data-field="gene"</th>}; }
else { $html .= qq{<th data-field="gene" data-filter-control="input" data-filter-control-placeholder="Gene name, description, ..."</th>}; }
$html .= $cgi->end_Tr();
$html .= "</thead>";
$html .= "<tbody>";
foreach my $score (sort {$b <=> $a} keys %$h_html_genes) {
	foreach my $gene_id (sort keys %{$h_html_genes->{$score}}) {
		my $html_gene = $h_html_genes->{$score}->{$gene_id};
		$html .= "<tr style='font-size:11px;'><td>".$html_gene."</td></tr>";
	}
}
$html .= "</tbody>";
$html .= "</table>";
$html =~ s/glyphicon-minus/fa fa-minus/g;

undef $h_html_genes;

print '.write_html.';
open (HTML, ">$html_file");
print HTML $html;
close (HTML);
undef $html;
print '.END';

if ($launch_job) {
	my $hRes;
	$hRes->{status} = "finished";
	$hRes->{session_id} = $export_session_id;
	$hRes->{phenotypes} = join(',', @lPhenos);
	if ($h_errors_found) {
		$hRes->{errors}  = join(', ', values %$h_errors_found);
	}
	open(my $out, ">", $outfile);
	print $out encode_json($hRes);
	close $out;
	exit(0);
}
else {
	my $hRes;
	$hRes->{session_id} = $export_session_id;
	if ($h_errors_found) {
		$hRes->{errors}  = join(', ', values %$h_errors_found);
	}
	my $json_encode = encode_json $hRes;
	print ".\",";
	$json_encode =~ s/{//;
	print $json_encode;
	exit(0);
}


#TODO: erreur certains VAR IDS ne passent pas les filtres (ex: 19_35739335_G_A dans KMT2B et filtres par defaut)
sub export_xls {
	my ($dejavu_variants, $session_id) = @_;
	$dejavu_variants->xls_export_session->load_with_storable($session_id);
	my $list_datas_annotations = $dejavu_variants->xls_export_session->prepare_generic_datas_variants();
	my $list_header = $dejavu_variants->xls_export_session->list_generic_header();
	$dejavu_variants->xls_export_session->add_page('Results', $list_header, $list_datas_annotations);
	my (@list_datas_patients, $h_by_patients);
	eval { $h_by_patients = $dejavu_variants->xls_export_session->get_specific_infos_stored('projects_patients_infos'); };
	if ($@) {}
	else {
		my $h_done;
		foreach my $var_id (sort keys %{$h_by_patients}) {
			foreach my $h_infos (@{$h_by_patients->{$var_id}}) {
				next if exists $h_done->{$var_id.'-'.$h_infos->{project_name}.'-'.$h_infos->{patient_name}.'-'.$h_infos->{'gene'}};
				my $h;
				$h->{'variation'} = $var_id;
				$h->{'project'} = $h_infos->{project_name};
				$h->{'patient'} = $h_infos->{patient_name};
				$h->{'model'} = $h_infos->{model};
				$h->{'he_ho'} = $h_infos->{heho};
				$h->{'gnomad ac'} = $h_infos->{'gnomad ac'};
				$h->{'gnomad ho'} = $h_infos->{'gnomad ho'};
				$h->{'gene(s)'} = $h_infos->{'gene'};
				$h->{'consequence'} = $h_infos->{'consequence'};
				$h->{'mane transcript'} = $h_infos->{'transcript'};
				$h->{'nomenclature'} = $h_infos->{'nomenclature'};
				$h->{'prot_nomenclature'} = $h_infos->{'prot_nomenclature'};
				my ($ac, $ratio) = split(', ', $h_infos->{ratio});
				my $dp = $h_infos->{dp};
				$ac =~ s/Reads://;
				$ratio =~ s/Ratio://;
				$h->{'nb reads'} = $ac;
				$h->{'dp'} = $dp;
				$h->{'ratio'} = $ratio.'%';
				$h_done->{$var_id.'-'.$h_infos->{project_name}.'-'.$h_infos->{patient_name}.'-'.$h_infos->{'gene'}} = undef;
				push(@list_datas_patients, $h);
			}
		}
		my @lLinesHeaderPatients = ('Variation', 'Project', 'Patient', 'Model', 'He_Ho', 'Dp', 'Nb Reads', 'Ratio', 'gnomAD AC', 'gnomAD HO', 'Gene(s)', 'Mane Transcript', 'Consequence', 'Nomenclature', 'Prot_Nomenclature');
		$dejavu_variants->xls_export_session->add_page('Projects Patients', \@lLinesHeaderPatients, \@list_datas_patients);
	}
	$dejavu_variants->xls_export_session->export();
	exit(0);
}

sub return_hashes_variants_in_session_export {
	my ($dejavu_variants, $hgenes, $hVariantsDetails) = @_;
	my (@lVar, $h_variants, $h_patients);
	foreach my $gene_id (keys %{$hGenes}) {
		foreach my $var_id (keys %{$hGenes->{$gene_id}}) {
			my $var = $dejavu_variants->project->_newVariant($var_id);
#			push(@lVar, $var);
			next if not exists $hVariantsDetails->{$var_id}->{can_export};
			my $pv = $hVariantsDetails->{$var_id}->{polyviewer_variant};
			my $hres = { %{$pv} };
			my $h;
			$h->{'chr'} = $hres->{chromosome};
			$h->{'var_id'} = $hres->{'id'};
			$h->{'rsname'} = $var->rs_name();
			$h->{'position'} = $hres->{'start'};
			$h->{'cadd_score'} = '-';
			$h->{'cadd_score'} = $var->cadd_score() if defined $var->cadd_score();
			$h->{'clinvar'} = '-';
			$h->{'clinvar'} = $var->text_clinvar() if ($var->text_clinvar() and  $var->text_clinvar() ne '-5' );
			$h->{'freq'} = '-';
			$h->{'freq'} = $var->percent() if defined $var->frequency() and $var->frequency();
			$h->{'gnomad ac'} = $var->getGnomadAC();
			$h->{'gnomad an'} = $var->getGnomadAN();
			$h->{'gnomad ho'} = $var->getGnomadHO();
			$h->{'ncboost_score'} = '-';
			$h->{'ncboost_score'} = $var->ncboost_score() if defined $var->ncboost_score();
			$h->{'cosmic'} = '-';
			$h->{'cosmic'} = $var->cosmic() if $var->cosmic();
			$h->{'hgmd_class'} = undef;
			$h->{'allele'} = $hres->{'allele'};
			$h->{'sequence'} = $hres->{'ref_allele'}.'/'.$hres->{'allele'};
			$h->{'dejavu'} = 'proj:'.$hres->{'dejavu_other_projects'}.', pat:'.$hres->{'dejavu_other_patients'}.' (ho:'.$hres->{'dejavu_other_patients_ho'}.')';
			$h->{'cadd'} = undef;
			$h->{'freq (%)'} = undef;
			$h->{'max_pop_freq'} = $hres->{'gnomad_max_pop_name'}.':'.$hres->{'gnomad_max_pop'};
			$h->{'min_pop_freq'} = $hres->{'gnomad_min_pop_name'}.':'.$hres->{'gnomad_min_pop'};
			my $h_genes_cons;
			foreach my $g_id (keys %{$hres->{hgenes}}) {
				my $g = $var->getProject->newGene($g_id);
				my $cons_g = $var->variationTypeInterface($g);
				$h->{genes}->{$g_id}->{external_name} = '-';
				eval {
					$h->{genes}->{$g_id}->{external_name} = $g->external_name();
				};
				$h->{genes}->{$g_id}->{description} = '-';
				eval { $h->{genes}->{$g_id}->{description} = $g->description(); };
				if ($@) {}
				$h->{genes}->{$g_id}->{phenotypes} = '-';
				eval { $h->{genes}->{$g_id}->{phenotypes} = $g->phenotypes(); };
				if ($@) {}
				foreach my $htr (@{$hres->{hgenes}->{$g_id}->{tr}}) {
					next if lc($htr->{consequence}) eq 'intergenic';
					$htr->{external_name} = $htr->{name};
					$htr->{max_splice_ai} = $htr->{spliceAI}.' '.$htr->{spliceAI_cat};
					$htr->{promoter_ai} = $htr->{promoterAI_score};
					$h->{genes}->{$g_id}->{transcripts}->{$htr->{name}} = $htr;
					my $t = $var->getProject->newTranscript($htr->{name});
					$htr->{cdna_position} = $t->translate_position( $var->start() );
					$htr->{polyphen_score} = $var->polyphenScore($t);
					$htr->{sift_score} = $var->siftScore($t);
					my $prot = $t->getProtein();
					if ($prot) {
						my $cds_pos = $var->getOrfPosition($prot);
						$cds_pos = '-' if (not $cds_pos or $cds_pos eq '.');
						$htr->{cds_position} = $cds_pos;
						my $prot_nom;
						eval { $prot_nom = $var->protein_nomenclature($prot); };
						if ($@) { $prot_nom = '-'; }
						$htr->{'prot_nomenclature'} = $prot_nom;
						$htr->{protein} = $prot->id();
						my $protAA = $var->getProteinAA($prot);
						my $chanAA = $var->changeAA($prot);
						$htr->{aa} = $protAA.'/'.$chanAA if ( $protAA and $chanAA );
						$htr->{protein_position} = $var->protein_nomenclature($prot);
					}
					if ($t->isMane()) {
						push(@{$h_genes_cons->{external_name}}, $g->external_name().' ('.$g->id.')');
						push(@{$h_genes_cons->{transcript}}, $t->name(). ' ('.$g->external_name().')');
						push(@{$h_genes_cons->{consequence}}, $htr->{consequence});
						if ($prot) {
							push(@{$h_genes_cons->{nomenclature}}, $htr->{'nomenclature'});
							push(@{$h_genes_cons->{prot_nomenclature}}, $htr->{'prot_nomenclature'});
						}
						else {
							push(@{$h_genes_cons->{nomenclature}}, '-');
							push(@{$h_genes_cons->{prot_nomenclature}}, '-');
							
						}
					}
				}
			}
			my $chr_id = $h->{'chr'};
			$chr_id = 23 if lc($chr_id) eq 'x';
			$chr_id = 24 if lc($chr_id) eq 'y';
			$chr_id = 25 if lc($chr_id) eq 'm';
			$chr_id = 25 if lc($chr_id) eq 'mt';
			$h_variants->{$chr_id}->{$var_id} = $h; 
			if (exists $hVariantsDetails->{$var_id}->{polyviewer_html_details_proj_pat}) {
				foreach my $h_tmp (@{$hVariantsDetails->{$var_id}->{polyviewer_html_details_proj_pat}}) {
					$h_tmp->{'gnomad ac'} = $var->getGnomadAC();
					$h_tmp->{'gnomad ho'} = $var->getGnomadHO();
					if (exists $h_genes_cons->{'consequence'}) {
						$h_tmp->{'consequence'} = join(' ; ', @{$h_genes_cons->{'consequence'}}) if exists $h_genes_cons->{'consequence'};
					}
					else {
						$h_tmp->{'consequence'} = $var->variationTypeInterface();
					}
					$h_tmp->{'gene'} = '-';
					$h_tmp->{'gene'} = join(' ; ', @{$h_genes_cons->{'external_name'}}) if exists $h_genes_cons->{'external_name'};
					$h_tmp->{'transcript'} = '-';
					$h_tmp->{'transcript'} = join(' ; ', @{$h_genes_cons->{'transcript'}}) if exists $h_genes_cons->{'transcript'};
					$h_tmp->{'nomenclature'} = '-';
					$h_tmp->{'nomenclature'} = join(' ; ', @{$h_genes_cons->{'nomenclature'}}) if exists $h_genes_cons->{'nomenclature'};
					$h_tmp->{'prot_nomenclature'} = '-';
					$h_tmp->{'prot_nomenclature'} = join(' ; ', @{$h_genes_cons->{'prot_nomenclature'}}) if exists $h_genes_cons->{'prot_nomenclature'};
					push(@{$h_patients->{$var_id}}, $h_tmp);
				}
			}
			undef $var;
		}
	}
	$dejavu_variants->project->buffer->dbh_deconnect();
#	$dejavu_variants->xls_export_session->store_variants_infos(\@lVar, $dejavu_variants->project());
#	$dejavu_variants->xls_export_session->{hash_variants_global} = $h_variants;
#	$dejavu_variants->xls_export_session->store_specific_infos('projects_patients_infos', $h_patients);
	return ($h_variants, $h_patients);	
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
			my $ncboost_filters_dir = $project->ncboost_parquet_dejavu_filter_path();
			if ($dejavu_variants->is_magic_user()) {
				$sql = "
					PRAGMA threads=$fork;
					SELECT 
					    a.pos, a.rocksdb_id, a.score AS ncboost
					FROM read_parquet('$ncboost_filters_dir/chr=$chr_id/*.parquet') a
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
					FROM read_parquet('$ncboost_filters_dir/chr=$chr_id/*.parquet') a
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

sub pass_frequence_or_clinvar_pathogenic {
	my ($chr_id, $rocks_id) = @_;
	if (is_rocks_id_clinvar_pathogenic($project->getChromosome($chr_id)->rocksdb('clinvar'), $rocks_id)) {
		my $h_dv = get_dejavu_he_ho($chr_id, $rocks_id);
		return if not $h_dv;
		return $h_dv;
	}
	return if not pass_gnomad_filter($chr_id, $rocks_id);
	my $h_dv = get_dejavu_he_ho($chr_id, $rocks_id);
	return if not $h_dv;
	return if $h_dv->{total} > $max_dejavu;
	return if $h_dv->{ho} > $max_dejavu_ho;
	return $h_dv;
}

sub pass_gnomad_filter {
	my ($chr_id, $rocks_id) = @_;
	my $h_gad = $project->getChromosome($chr_id)->rocksdb('gnomad')->value($rocks_id);
	return if $h_gad and exists $h_gad->{ac} and $h_gad->{ac} > $dejavu_variants->max_gnomad_ac();
	return if $h_gad and exists $h_gad->{ho} and $h_gad->{ho} > $dejavu_variants->max_gnomad_ac_ho();
	return 1;
}

sub get_dejavu_he_ho {
	my ($chr_id, $rocks_id) = @_;
	my ($id, $alt) = split('!', $rocks_id);
	my ($res, $nb_pat_he, $nb_pat_ho);
	if ($alt eq 'A' or $alt eq 'T' or $alt eq 'G' or $alt eq 'C') {
		$res = $project->getChromosome($chr_id)->rocks_dejavu->dejavu($rocks_id);
		return if not $res;
		return if not pass_dejavu_filter_only_my_projects($res);
		($nb_pat_he, $nb_pat_ho) = get_dv_he_ho_from_request($res);
		my $h;
		$h->{details} = $res;
		$h->{he} = $nb_pat_he;
		$h->{ho} = $nb_pat_ho;
		$h->{total} = $nb_pat_he + $nb_pat_ho;
		return $h;
	}
	else {
		$res = $project->getChromosome($chr_id)->rocks_dejavu->dejavu_interval($id -1 , $id +1);
		return if not $res;
		foreach my $dv_rocks_id (keys %{$res}) {
			if (exists $h_res_duck->{$chr_id.'!'.$dv_rocks_id}) {
				return if not pass_dejavu_filter_only_my_projects($res->{$dv_rocks_id});
				($nb_pat_he, $nb_pat_ho) = get_dv_he_ho_from_request($res->{$dv_rocks_id});
				my $h;
				$h->{details} = $res->{$dv_rocks_id};
				$h->{he} = $nb_pat_he;
				$h->{ho} = $nb_pat_ho;
				$h->{total} = $nb_pat_he + $nb_pat_ho;
				die if $h->{total} == 0;
				return $h;
			}
		}
	}
	return;
}

sub pass_dejavu_filter_only_my_projects {
	my ($res) = @_;
	foreach my $proj_id (keys %{$res}) {
		return 1 if exists $dejavu_variants->{hash_users_projects}->{$proj_id};
	}
	return;
}


sub get_dv_he_ho_from_request {
	my ($res) = @_;
	my $nb_pat_he = 0;
	my $nb_pat_ho = 0;
	foreach my $proj_id (keys %{$res}) {
		next if ($proj_id eq 'polybtf');
		$nb_pat_he += $res->{$proj_id}->{he};
		$nb_pat_ho += $res->{$proj_id}->{ho};
	}
	return ($nb_pat_he, $nb_pat_ho);
}

sub is_rocks_id_clinvar_pathogenic {
	my ($no_clinvar_pathogenic, $rocks_id) = @_;
	my $h_clinvar_pathogenic = $no_clinvar_pathogenic->clinvar($rocks_id);
	return if not $h_clinvar_pathogenic;
	return if not exists $h_clinvar_pathogenic->{score};
	return 1 if $h_clinvar_pathogenic->{score} >= 5;
	return;
}

