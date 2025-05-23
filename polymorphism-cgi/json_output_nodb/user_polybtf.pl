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
use lib "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/";

require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/html_polygenescout.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";

use GBuffer;
use export_data;
use JSON;
use polyweb_dude;
use VcfMerge;
use GenBoNoSql;
use Set::IntervalTree;
use Spreadsheet::WriteExcel;
use Compress::Snappy;
use Storable qw(store retrieve freeze dclone thaw);
use POSIX qw(strftime);
use CGI::Session;
use html; 
use Carp;
use Cache_Commons;
use QueryVectorFilter;
use IO::Handle;
use update_variant_editor;
use Getopt::Long;
use File::Basename;
use MIME::Base64;


my $cgi = new CGI();
my $login = $cgi->param('login');
my $pwd = $cgi->param('pwd');
my $count = $cgi->param('count');

#my @header_transcripts = ("consequence","enst","nm","nomenclature");
my @header_transcripts = ("consequence","enst","nm","ccds","appris","exon","nomenclature","codons","codons_AA", "polyphen","sift",'alphamissense',"cadd","revel","dbscsnv",'spliceAI');


my $buffer = GBuffer->new;
my $annotation = $buffer->config->{polybtf_default_releases}->{default};
my $btf_path = $buffer->get_polybtf_path($annotation);


print $cgi->header('text/json-comment-filtered');
print "{\"progress\":\".";
my $file = get_json_file_name();
load_json($file) if -e $file;


my $project_name_hg38 = $buffer->getRandomProjectName('HG38_DRAGEN', $annotation);
my $project = $buffer->newProject( -name => $project_name_hg38 );
#warn $login.' '.$pwd;
my $res = $buffer->getQuery()->getProjectListForUser($login, $pwd);

my $h_liftover;

#my $h_json;
#if ($annotation eq '43.20') {
#	my $json_btf = '/data-isilon/polycache/polybtf/43.20/new_hgmd_clinvar.json';
#	open (JSON, $json_btf);
#	my $json = <JSON>;
#	$h_json = decode_json $json;
#	close (JSON);
#	my $buffer_hg19 = GBuffer->new;
#	my $project_name_hg19 = $buffer->getRandomProjectName('HG19', $annotation);
#	my $project_hg19 = $buffer_hg19->newProject( -name => $project_name_hg19 );
#	my @lVar;
#	foreach my $var_id (keys %$h_json) {
#		my $var = $project_hg19->_newVariant($var_id);
#		my @ltmp = split ('_', $var->id);
#		next if not $ltmp[-2] =~ /[ATGC]+/;
#		next if not $ltmp[-1] =~ /[ATGC]+/;
#		#warn ref($var). ' -> '.$var->id;
#		push(@lVar, $var);
#	}
#	my @lVarId_38;
#	warn $project_hg19->lift_genome_version;
#	my $lift = liftOver->new(project=>$project_hg19, version=>$project_hg19->lift_genome_version);
#	$lift->lift_over_variants(\@lVar);
#	foreach my $var (@lVar) {
#		#warn ref($var).' -> hg19:'.$var->id.' -> hg38:'.$var->lift_over('HG38')->{id};#.' -> rocks:'.$var->rocksdb_id();
#		push(@lVarId_38, $var->lift_over('HG38')->{id});
#	}
#	$project_hg19 = undef;
#	$buffer_hg19 = undef;
#	
#	my @lH;
#	foreach my $var_id (@lVarId_38) {
#		my $h;
#		my $var = $project->_newVariant($var_id);
#		$h->{'chr'} = $var->getChromosome->id;
#		$h->{'pos'} = $var->start;
#		$h->{'id'} = $var->id;
#		$h->{'rocksid'} = $var->rocksdb_id;
#		push (@lH, $h);
#	}
#	my $json_encode = encode_json \@lH;
#	open (JSON, ">/data-isilon/public-data/repository/HG38/polybtf/20/lift_polybtf.json");
#	print JSON $json_encode;
#	close (JSON);
#	
#	die;
#}

my @list_ids;
my ($h_proj, $h_years);
foreach my $h (@$res) {
	next if not $h->{name} =~ /NGS20/;
	my $annot = $buffer->getQuery->getPublicDatabaseVersion($h->{id});
	next if $annot >= 20;	
	push(@list_ids, int($h->{id}));
	$h_proj->{$h->{id}}->{name} = $h->{name};
	$h_proj->{$h->{id}}->{description} = $h->{description};
	my @lTmp = split('_', $h->{name});
	my $year = $lTmp[0];
	$year =~ s/NGS//;
	$h_years->{$year} = undef;
}

my @lbtf;
foreach my $year (keys %{$h_years}) {
	my $file = $btf_path.'/polybtf.'.$year.'.parquet';
	$file =~ s/\/\//\//g;
	push (@lbtf, "'".$file."'") if -e $file; 
}

exit(0) if scalar(@lbtf) == 0;
exit(0) if scalar(@list_ids) == 0;

my $ids = join(', ', sort @list_ids);
my $sql = "PRAGMA threads=6;";
if ($count) {
	$sql .= "SELECT count(DISTINCT rocksid) FROM read_parquet([".join(',', @lbtf)."]) WHERE project in ($ids);";
}
else {
	$sql .= "SELECT * FROM read_parquet([".join(',', @lbtf)."]) WHERE project in ($ids);";
}

my ($h_var_already_filtered, $out2, $hVarHtml, $hProjVar);
my $duckdb = $buffer->software('duckdb');
my $cmd = qq{set +H | $duckdb -json -c "$sql"};
my $json_duckdb = `$cmd`;
if ($json_duckdb) {
	my $decode = decode_json $json_duckdb;
	if ($count) {
		my $hres;
		$hres->{count} = $decode->[0]->{'count(DISTINCT rocksid)'};
		printJson($hres);
		exit(0);
	}
	else {
		my $i = 0;
		foreach my $h (@$decode) {
			$i++;
			if ($i ==50) {
				print '.';
				$i = 0;
			}
			my $project_id = $h->{project};
			my $chr38 = $h->{chr38};
			my $pos38 = $h->{pos38};
			my $rocksid = $h->{rocksid};
			next if exists $h_var_already_filtered->{$chr38.'-'.$rocksid};
			
			my $allele = $h->{allele};
			my $hgmd = $h->{class};
#			my $hgmd_old = $h->{hgmd_old},
			my $clinvar = $h->{clnsig};
			my $clinvar_old = $h->{clnsig_OLD};
			my $transmission_id = $h->{transmissions};
			my $hgmd_gene = $h->{gene};
			my $hgmd_phen = $h->{phen};
			$hgmd_phen =~ s/_/ /g;
			$hgmd_phen =~ s/%2C//g;
			
			my $chr = $project->getChromosome($chr38);
			my $var_id = $chr->transform_rocksid_to_varid($rocksid);
			my $var = $project->_newVariant($var_id);
			
			my $ok = 1;
			$ok = undef if $var->getGnomadAC() > 50;
			$h_var_already_filtered->{$chr38.'-'.$rocksid} = undef if not $ok;
			next if not $ok;
			
			$ok = undef if $var->getGnomadHO() > 10;
			$h_var_already_filtered->{$chr38.'-'.$rocksid} = undef if not $ok;
			next if not $ok;
			
			$ok = undef if $var->other_patients() > 50;
			$h_var_already_filtered->{$chr38.'-'.$rocksid} = undef if not $ok;
			next if not $ok;
			
			$ok = undef if $var->other_patients_ho() > 10;
			$h_var_already_filtered->{$chr38.'-'.$rocksid} = undef if not $ok;
			next if not $ok;
			
			
			if (not exists $hVarHtml->{$var_id}) {
				my $h_var;
				$h_var->{html}->{done_here} = 1;
				$h_var->{html}->{no_css_polydiag} = 1;
				$h_var->{value}->{id} =  $var->id;
				$h_var->{html}->{id} =  $var->id;
				$h_var->{value}->{type} = $var->type;
				$h_var->{html}->{type} = $var->type;
				my $vn = $var->id();
				$vn =~ s/_/-/g;
				$vn =~ s/chr//;
				
				$hVarHtml->{$var_id}->{hgmd} = qq{<span style='color:red;'><b><i>$hgmd</b></i></span><br><span style='color:green;'><i>$hgmd_phen</i></span><br><span style='color:blue;'><i>$hgmd_gene</i></span>};
				
				if ($clinvar_old ne $clinvar) {
					my $color_old = 'grey';
					$color_old = 'red' if lc($clinvar_old) eq 'pathogenic';
					my $color = 'grey';
					$color = 'red' if lc($clinvar) eq 'pathogenic';
					$hVarHtml->{$var_id}->{clinvar} = qq{<s><span style='color:$color_old;'><i>$clinvar_old</i></span></s><br><i class="bi bi-caret-down-fill" style="font-size: 1.5rem; font-weight: bold;"><br><span style='color:$color;'><b><i>$clinvar</b></i></span>};
				}
				else {
					my $color = 'grey';
					$color = 'red' if lc($clinvar) eq 'pathogenic';
					$hVarHtml->{$var_id}->{clinvar} = qq{<span style='color:$color;'><b><i>$clinvar</b></i></span>};
				}
				
				$hVarHtml->{$var_id}->{table_vname} = update_variant_editor::vname2($var, $h_var);
				$hVarHtml->{$var_id}->{table_gnomad} = update_variant_editor::table_gnomad($var);
				$hVarHtml->{$var_id}->{table_gnomad} =~ s/gnomad_r2_1/gnomad_r4/;
				$hVarHtml->{$var_id}->{table_varsome} = update_variant_editor::vvarsome($h_var);
				$hVarHtml->{$var_id}->{table_dejavu} = update_variant_editor::vdejavu($var, $h_var);
				
				$hVarHtml->{$var_id}->{table_gene} = undef;
				my @lGenes = @{$var->getGenes()};
				foreach my $gene (@lGenes) {
					$h_var->{genes}->{$gene->id} = update_variant_editor::construct_hash_transcript($var, $cgi, \@header_transcripts, 2, $gene);
					
					my $is_good_gene;
					$is_good_gene = 1 if $gene->external_name() eq $hgmd_gene;
					
					my $html_gene;
					if ($is_good_gene) { 
						$html_gene = qq{<div style="opacity: 1">};
						$html_gene .= '<span>Gene: <b>'.$gene->external_name().'</b> | '.$gene->id().'</span><br>';
						$html_gene .= update_variant_editor::table_transcripts($h_var->{genes}->{$gene->id}, \@header_transcripts, 1);
						$html_gene .= '</div>';
						$html_gene .= '</div>';
						$html_gene .= '<br>' if $hVarHtml->{$var_id}->{table_gene};
						$hVarHtml->{$var_id}->{table_gene} = $html_gene.$hVarHtml->{$var_id}->{table_gene};
					}
					else {
						$html_gene = '<br>' if $hVarHtml->{$var_id}->{table_gene};
						$html_gene .= qq{<div style="opacity: 0.33">};
						$html_gene .= '<span>Gene: <b>'.$gene->external_name().'</b> | '.$gene->id().'</span><br>';
						$html_gene .= update_variant_editor::table_transcripts($h_var->{genes}->{$gene->id}, \@header_transcripts, 1);
						$html_gene .= '</div>';
						$html_gene .= '</div>';
						$hVarHtml->{$var_id}->{table_gene} .= $html_gene;
					}
				}
			}
			my $ok;
			my @list_parquets;
			my $parquet = $buffer->dejavu_parquet_dir().'/'.$h_proj->{$project_id}->{name}.'.'.$project_id.'.parquet';
			if (-e $parquet) {
				push(@list_parquets, "'".$parquet."'");
				my $html = get_from_duckdb_project_patients_infos($var, \@list_parquets);
				if ($html) {
					$hVarHtml->{$var_id}->{table_patients} .= "<br>" if exists $hVarHtml->{$var_id}->{table_patients};			
					$hVarHtml->{$var_id}->{table_patients} .= $html;
					$ok = 1;	
				}
			}
			elsif (-e $parquet.'.no_dejavu') {
				warn "\n\n";
				warn $parquet.'.no_dejavu';
				die;
			}
			else {
				warn $parquet;
				die;
			}
			if (not $ok) {
				$hVarHtml->{$var_id}->{table_patients} .= qq{<span><b>}.$h_proj->{$project_id}->{name}.qq{</b> (pb)</span>};
			}
		}
		
		my $fsize = "font-size:10px";

		$out2 .= $cgi->start_div();
		$out2 .= qq{<table data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-v-align="both" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="25" data-page-list="[25, 50, 100, 200, 300]" data-resizable='true' id='table_variants' class='table table-striped' style='font-size:13px;'>};
		$out2 .= "<thead>";
		$out2 .= $cgi->start_Tr({style=>"background-color:#E9DEFF;$fsize"});
		$out2 .= qq{<th data-field="project_name" data-filter-control="input" data-sortable="true">Project(s) / Patient(s)</th>};
		$out2 .= qq{<th data-field="var_name" data-filter-control="input" data-sortable="true" data-filter-control-placeholder=7-15456-A-G">Var Name</th>};
		$out2 .= qq{<th data-field="varsome">Varsome</th>};
		$out2 .= qq{<th data-field="gnomad">gnomAD</th>};
		$out2 .= qq{<th data-field="dejavu">DejaVu</th>};
		$out2 .= qq{<th data-field="hgmd" data-filter-control="input" data-sortable="true">HGMD</th>};
		$out2 .= qq{<th data-field="clinvar" data-filter-control="input" data-sortable="true">Clinvar</th>};
		$out2 .= qq{<th data-field="annotation" data-filter-control="input" data-filter-control-placeholder="Stop / c.73A>C / exon2 / ENST00000145855">Annotations</th>};
		$out2 .= $cgi->end_Tr();
		$out2 .= "</thead>";
		$out2 .= "<tbody>";
		foreach my $var_id (sort keys %{$hVarHtml}) {
			$out2 .= $cgi->start_Tr();
			$out2 .= qq{<td><div style="max-height:250px;max-width:400px;overflow-y:auto;">}.$hVarHtml->{$var_id}->{table_patients}.qq{</div></td>};
			$out2 .= qq{<td><center>}.$hVarHtml->{$var_id}->{table_vname}.qq{<br>[HG19: }.$h_liftover->{$var_id}.qq{]</center></td>};
			$out2 .= qq{<td><center>}.$hVarHtml->{$var_id}->{table_varsome}.qq{</center></td>};
			$out2 .= qq{<td><center>}.$hVarHtml->{$var_id}->{table_gnomad}.qq{</center></td>};
			$out2 .= qq{<td><center>}.$hVarHtml->{$var_id}->{table_dejavu}.qq{</center></td>};
			
			#HGMD
			$out2 .= qq{<td><center>}.$hVarHtml->{$var_id}->{hgmd}.qq{</center></td>};
			$out2 .= qq{<td><center>}.$hVarHtml->{$var_id}->{clinvar}.qq{</center></td>};
			#CLINVAR
			
			
			$out2 .= qq{<td><div style="max-height:250px;overflow-y:auto;">}.$hVarHtml->{$var_id}->{table_gene}.qq{</div></td>};
			$out2 .= $cgi->end_Tr();
		}
		$out2 .= "</tbody>";
		$out2 .= "</table>";
		$out2 .= "</div>";
		$out2 .= "<br>";
		
	}
}

#print $out2;
my $hRes;
$hRes->{html} = $out2;
$hRes->{var_filtred} = scalar(keys %$h_var_already_filtered);
printJson($hRes);

 


sub get_table_project_patients_infos {
	my ($project_name, $hash, $hVar_infos) = @_;
	my $locus_hg19 = $hVar_infos->{locus_hg19};
	my $locus_hg38 = $hVar_infos->{locus_hg38};
	
	my $nb_he = $hash->{he};
	my ($h_infos_patients, $h_tmp_pat);
	my $nb_pat = 0;
	foreach my $pat_id (unpack("w*",decode_base64($hash->{patients}))) {
		$nb_pat++;
		$h_infos_patients->{$nb_pat}->{id} = $pat_id;
		$h_tmp_pat->{$pat_id} = $nb_pat;
	}
	my $i = 0;
	$nb_pat = 1;
	foreach my $info (unpack("w*",decode_base64($hash->{dp_ratios}))) {
		$i++;
		if ($i == 1) { $h_infos_patients->{$nb_pat}->{dp} = $info; }
		elsif ($i == 2) {
			my $ratio = ($info / $h_infos_patients->{$nb_pat}->{dp}) * 100;
			my $text = 'AC:'.$info.' ('.int($ratio).'%)';
			$h_infos_patients->{$nb_pat}->{ratio} = $text;
		}
		elsif ($i == 3) {
    		my $model;
    		if ($info == 1) { $model = 'solo'; }
    		elsif ($info == 2) { $model = 'father'; }
    		elsif ($info == 4) { $model = 'mother'; }
    		elsif ($info == 8) { $model = 'both'; }
    		elsif ($info == 16) { $model = 'is_parent'; }
    		elsif ($info == 32) { $model = 'recessif'; }
    		elsif ($info == 64) { $model = 'dominant'; }
    		elsif ($info == 128) { $model = 'denovo'; }
    		elsif ($info == 256) { $model = 'strict_denovo'; }
    		elsif ($info == 512) { $model = 'error'; }
    		else { $model = 'error2'; }
			$h_infos_patients->{$nb_pat}->{model} = $model;
			
			$i = 0;
			$nb_pat++;
		}
	}
	
	my $b = new GBuffer;
	my $p = $b->newProject( -name => $project_name);
	
	my @lPat = @{$p->getPatients()};
	return undef if scalar(@lPat) == 0;
	
	my $found_healthy_patient;
	foreach my $pat (@lPat) {
		next if not exists $h_tmp_pat->{$pat->id};
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{name} = $pat->name;
		my $icon = $pat->small_icon();
		$icon =~ s/"/'/g;
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{status} = $icon;
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{sex} = '-';
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{sex} = 'male' if $pat->sex() eq '1';
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{sex} = 'female' if $pat->sex() eq '2';
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{status_txt} = '-';
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{status_txt} = 'healthy' if $pat->status() eq '1';
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{status_txt} = 'ill' if $pat->status() eq '2';
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{family} = $pat->getFamily->name();
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{description} = $p->description();
		if (int($h_tmp_pat->{$pat->id}) <= $nb_he) { $h_infos_patients->{$h_tmp_pat->{$pat->id}}->{heho} = 'He'; }
		else { $h_infos_patients->{$h_tmp_pat->{$pat->id}}->{heho} = 'Ho'; }
	}
	
	return undef if not $h_infos_patients or scalar keys %$h_infos_patients == 0;
	
#	my $gene_name =$gene_used->external_name();
#	my $fam = $patient->getFamily();
	my $is_solo_trio = 'SOLO';
#	$is_solo_trio = 'TRIO' if $fam->isTrio();
#	my $project_name = $patient->getProject->name();
	
	my $description = $p->description();
	my @l_users = @{$p->get_list_emails()};
#	my $patient_name = $patient->name();
#	my $pheno = undef;
#	$pheno = $h_var->{html}->{pheno_name} if ($h_var and exists $h_var->{html}->{pheno_name});
	
	my $color = "#c1c1c1";
	my $model;
	my $patient_heho = "-";
	
	my $nb_col_span = 6;
#	$nb_col_span = 7 if ($var->getProject->isGenome() && $var->isCnv);
#	my $hstatus = $patient->validation_status();
#	my $hval = $patient->validations();
	
	
	my $table_trio = qq{ <div> };
	$table_trio .= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:3px"});
	
	
	my @lPhenotypes = @{$p->phenotypes()};
	my $pheno = join(', ', sort @lPhenotypes);
	
	my $version = $p->annotation_genome_version();
	my $color_version = '#85f283';
	$color_version = '#f2c37c' if $version eq 'HG38';
	
	my $users = join("<br>", @l_users);
	my $proj_text = qq{<button style="color:black;" onclick="get_popup_users('$users');">Users</button> - <b>$project_name</b>};
	$proj_text .= qq{<sup>defidiag</sup>} if $p->isDefidiag; 
	$proj_text .= " - <span style='color:#82d0f5;'>$pheno</span>" if ($pheno);
	$proj_text .= "<br>";
	$proj_text .= "<b><span style='color:$color_version;'>$version</span></b>  ";
	$proj_text .= "$description";
		
	my (@l_pat_names, @l_pat_bam);
	foreach my $id (sort keys %$h_infos_patients) {
		next if not (exists $h_infos_patients->{$id}->{name});
		push(@l_pat_names, $h_infos_patients->{$id}->{name});
		push(@l_pat_bam, $p->getPatient($h_infos_patients->{$id}->{name})->bamUrl());
	}
	my $pnames = join(';', @l_pat_names);
	my $f = join(';', @l_pat_bam);
	my $gn = $p->getVersion();
	my $chr_name = $hVar_infos->{chr_id};
	my ($locus, $start);
	if ($p->annotation_genome_version() eq 'HG19') {
		$locus = $locus_hg19;
		$start = $hVar_infos->{start_hg19};
	}
	elsif ($p->annotation_genome_version() eq 'HG38') {
		$locus = $locus_hg38;
		$start = $hVar_infos->{start_hg38};
	}
	my $a0 = $hVar_infos->{ref_all};
	my $a1 = $hVar_infos->{alt_all};
	my $igb_b = qq{<button class='igvIcon2' onclick='launch_web_igv_js("$project_name","$pnames","$f","$locus","/","$gn")' style="color:black"></button>};	
	my $alamut_b = qq{<button class="alamutView3" onClick ="displayInAlamut('$chr_name',$start,['$a0','$a1']);"></button>};

	my $no_header_project ;
	my $no_header_project_pat;
	
	$table_trio .= $cgi->start_Tr({style=>"background-color:#949292; color:white;"});
	my $b_others_var;
	my $pat_name = 'ALL';
	my $var_id = 'VARID';
	my $father_trans = undef;
	my $mother_trans = undef;
	my $other_trans = undef;
	
	my $img_child = qq{<img style='width:14px;height:14px;' src='/icons/Polyicons/baby-boy.png'>};
	#$img_child = qq{<img style='width:14px;height:14px;' src='/icons/Polyicons/baby-girl.png'>} if ($patient->sex() == 2);
#	my $cmd_others = qq{view_var_from_proj_gene_pat('$project_name','$gene_init_id','$pat_name','$var_id');};
	my $cmd_others = '';
	$other_trans = qq{<button style="text-align:middle;vertical-align:top;width:28px;border:solid 0.5px black;background-color:white;" onClick="$cmd_others">$img_child</button>};

	$table_trio .= $cgi->td({colspan=>($nb_col_span)-1}, $proj_text);
	$table_trio .= $cgi->td({colspan=>1, style=>"padding:0px;width:60px;"}, '<center>'.$igb_b.' '.$alamut_b.'</center>');
	
	$no_header_project_pat = 1;
	$table_trio .= $cgi->end_Tr();
	if ($no_header_project_pat) {
		$table_trio .= '</table>';
		$table_trio .= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:3px"});
	}
	
	foreach my $nb (keys %{$h_infos_patients}) {
		next if not exists $h_infos_patients->{$nb}->{name};
		my $patient_name = $h_infos_patients->{$nb}->{name};
		my $patient = $p->getPatient($patient_name);
		my $patient_status = $h_infos_patients->{$nb}->{status};
		my $patient_heho = $h_infos_patients->{$nb}->{heho};
		my $dp = 'DP:'.$h_infos_patients->{$nb}->{dp};
		my $perc_allele = $h_infos_patients->{$nb}->{ratio};
		my $model = $h_infos_patients->{$nb}->{model};
		
		if ($model eq 'denovo' or $model eq 'strict_denovo' or $model eq 'dominant') {
			$color = 'background-color:#f54e4e';
			$model = ucfirst($model);
		}
		elsif ($model eq 'recessif') {
			$color = 'background-color:#e99ff7';
			$model = ucfirst($model);
		}
		elsif ($model eq 'father') {
			$color = 'background: linear-gradient(to right, white, #ddfbff);';
			if ($p->getPatient($patient_name)->getFamily->getFather->isIll()) { $model = qq{<img src="/icons/Polyicons/male-d.png">}; }
			else { $model = qq{<img src="/icons/Polyicons/male-s.png">}; }
		}
		elsif ($model eq 'mother') {
			$color = 'background: linear-gradient(to right, white, #ffddfd);';
			if ($p->getPatient($patient_name)->getFamily->getMother->isIll()) { $model = qq{<img src="/icons/Polyicons/female-d.png">}; }
			else { $model = qq{<img src="/icons/Polyicons/female-s.png">}; }
		}
		elsif ($model eq 'both') {
			$color = 'background: linear-gradient(to right, #ddfbff, #ffddfd);';
			if ($p->getPatient($patient_name)->getFamily->getFather->isIll()) { $model = qq{<img src="/icons/Polyicons/male-d.png">}; }
			else { $model = qq{<img src="/icons/Polyicons/male-s.png">}; }
			if ($p->getPatient($patient_name)->getFamily->getMother->isIll()) { $model .= qq{<img style="padding-left:5px;" src="/icons/Polyicons/female-d.png">}; }
			else { $model .= qq{<img style="padding-left:5px;" src="/icons/Polyicons/female-s.png">}; }
			 
		}
		elsif ($model eq 'solo') {
			$color = 'background: linear-gradient(to right, white, #c7c6c5);';
			$model = ucfirst($model);
		}
		else { $color = 'background-color:white'; }
	
		my $local_text;
		$table_trio .= $cgi->start_Tr({style=>"$color;"});
		$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, "<span style='text-align:left;'>$patient_name $local_text</span>");
		$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, "<div>".$patient_status."</div>");
		$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, "$patient_heho");
	#	if ($var->getProject->isGenome() && $var->isCnv) {
	#		$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, 'pr:'.$var->pr($patient));
	#		$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, 'sr:'.$var->sr($patient));
	#		eval {
	#			my $cnv_score = sprintf("%.2f", log2($patient->cnv_value_dude($var->getChromosome->name,$var->start,$var->start+$var->length)));
	#			$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, 'cnv_score:'.$cnv_score);
	#		};
	#		if ($@) {
	#			$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, 'cnv_score:pb');
	#		}
	#	}
	#	else {
			#my $perc_allele = $var->getPourcentAllele($patient);
			
#			return 'perc_all_filter' if ($filter_perc_allelic_min and $perc_allele < $filter_perc_allelic_min);
#			return 'perc_all_filter' if ($filter_perc_allelic_max and $perc_allele > $filter_perc_allelic_max);
#			$perc_allele .= "%" if ($perc_allele ne '-');
			$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, $perc_allele);
			$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, $dp);
	#	}
		$table_trio .= $cgi->td({style=>"text-align:center;vertical-align:middle;"}, $model);
		$table_trio .= $cgi->end_Tr();
	#	if ($local_text_tab) {
	#		$table_trio .=qq{ <div onClick="document.getElementById('$table_validation_id').style.display='none';" id='$table_validation_id' style='display:none;'>$local_text_tab</div> };
	#	}
	}
	$table_trio .= "</table></div>";
	
	return ($h_infos_patients, $table_trio, lc($model));
}

sub get_from_duckdb_project_patients_infos {
	my ($var, $list_files) = @_;
	return if scalar(@$list_files) == 0;
	my @list_table_trio;
	my $sql = "PRAGMA threads=6; SELECT project,chr38,chr19,pos38,pos19,he,allele,patients,dp_ratios FROM read_parquet([".join(', ', @$list_files)."])";
	
	my $find_pos_s = $var->start() - 20;
	my $find_pos_e = $var->start() + 20;
	
	my $h_projects_patients;
	
	if ($var->getProject->current_genome_version() eq 'HG38') {
		$sql .= " WHERE chr38='".$var->getChromosome->id()."' and pos38 BETWEEN '".$find_pos_s."' and '".$find_pos_e."';" ;
		
		my $duckdb = $buffer->software('duckdb');
		my $cmd = qq{set +H | $duckdb -json -c "$sql"};
		my $json_duckdb = `$cmd`;
		if ($json_duckdb) {
			my $decode = decode_json $json_duckdb;
			my $h_by_proj;
			foreach my $h (@$decode) {
				#next if $h->{'chr38'} ne $var->getChromosome->id;
				my $var_start = $var->start();
				$var_start-- if $var->isInsertion() or $var->isDeletion();
				next if $h->{'pos38'} ne $var_start;
				
#				my $var_all = $h->{'allele'};
#				$var_all =~ s/\+//;
				
#				if (not $var->isDeletion) {
#					next if $var_all ne $var->var_allele();
#				}
				
				my $project_id = $h->{project};
				my $project_name = $h_proj->{$project_id}->{name};
				$h_by_proj->{$project_name} = $h;
			}
			
			my $hVar_infos;
			$hVar_infos->{locus_hg19} = $var->getChromosome->id().":".$var->lift_over('HG19')->{position}."-".$var->lift_over('HG19')->{position};;
			$hVar_infos->{locus_hg38} = $var->getChromosome->id().":".$var->start."-".$var->end;
			$hVar_infos->{start_hg19} = $var->lift_over('HG19')->{position};
			$hVar_infos->{start_hg38} = $var->start;
			$hVar_infos->{chr_id} = $var->getChromosome->id();
			$hVar_infos->{ref_all} = $var->ref_allele();
			$hVar_infos->{alt_all} = $var->var_allele();
			$hVar_infos->{alt_all} = '*' unless $var->var_allele();
			foreach my $project_name (reverse sort keys %$h_by_proj) {
				my $h = $h_by_proj->{$project_name};
				my ($h_infos_patients, $table_trio, $model) = get_table_project_patients_infos($project_name, $h, $hVar_infos);
				push(@list_table_trio, $table_trio) if $table_trio;
				foreach my $id (sort keys %$h_infos_patients) {
					my $p_name = $h_infos_patients->{$id}->{'name'};
					$h_projects_patients->{$project_name}->{$p_name} = $h_infos_patients->{$id};
				}
			}
			$h_liftover->{$var->id} = $var->lift_over('HG19')->{id};
		}	
	}
	else {
		warn "HG19";
		warn $sql;
		confesss("HG19!!");
	}
	if (scalar(@list_table_trio) >= 1) {
		my $html = join("<br>",@list_table_trio);
		return ($h_projects_patients, $html);
	}
	return undef;
}

sub printJson {
	my ($hashRes, $test) = @_;
	my $json_encode = encode_json $hashRes;
	print ".\",";
	$json_encode =~ s/{//;
	save_json($json_encode);
	print $json_encode;
	exit(0);
}

sub get_json_file_name {
	my $file = $btf_path.'/users/';
	if ($buffer->getQuery->isLoginSTAFF($login)) { $file .= 'staff.json'; }
	else { $file .= $login.'.json'; }
	return $file;
} 

sub save_json {
	my ($json) = @_;
	return if $count;
	my $file = get_json_file_name();
	open (FILE, '>'.$file);
	print FILE print "{\"progress\":\".";
	print FILE ".\",".$json;
	close (FILE);
}

sub load_json {
	my ($file) = @_;
	return if $count;
	open (FILE, $file);
	my $json = <FILE>;
	close (FILE);
	chomp($json);
	print $json;
	exit(0);
}

