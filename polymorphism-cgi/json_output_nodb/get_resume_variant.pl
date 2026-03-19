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
use List::MoreUtils qw(natatime);
use MIME::Base64;
use Storable qw(store retrieve freeze dclone thaw);
use JSON;
use MCE::Loop;


my $cgi = new CGI();

print $cgi->header('text/json-comment-filtered');
print "{\"progress\":\".";

my $variant = $cgi->param('variant');
my $release = $cgi->param('release');

$release = 'HG38' if not $release;

my $buffer = new GBuffer;
my $project_name = $buffer->getRandomProjectName($release);
my $project = $buffer->newProject( -name => $project_name );

if ($variant =~ /-/) {
	my @lTmp = split('-', $variant);
	if (length($lTmp[-2]) == length($lTmp[-1])) {
		$variant =~ s/-/_/g;
	}
	else {
		$variant = $lTmp[0].'_'.($lTmp[1] + 1).'_'.$lTmp[2].'_'.$lTmp[3];
	}
}
my $var = $project->_newVariant($variant);

my $h;
$h->{span_release} = $release;
$h->{id} = $var->id();
$h->{span_gnomad_id} = $var->name();
$h->{span_genomic_rocks_id} = $var->genomic_rocksdb_id();

my $vp = PolyviewerVariant->new();
$vp->setLmdbVariant($var);
$vp->{hgenes} = {};
$vp->{genes_id} = [];
my $code = 0;
foreach my $g (@{$var->getGenes}){
	my $h = $vp->set_gene($var,$g);
	$h->{code} = $code;
	$vp->{hgenes}->{$g->id} = $h;
	push(@{$vp->{genes_id}},$g->id);
	$code ++;
}
$vp->{hpatients} = {};
$vp->{patients_id} = [];

my $patient = $project->getPatients->[0];
my $print_html = polyviewer_html->new( project=>$project, patient=>$patient, bgcolor=>"background-color:#607D8B" );

$print_html->variant($vp);
$h->{html_mobidetails} = $print_html->mobidetails();
$h->{html_gnomadurl} = $print_html->gnomadurl();
$h->{html_alamuturl} = $print_html->alamuturl();
$h->{html_varsome} = $print_html->varsome();
$h->{html_gnomad} = $print_html->gnomad();



### ANNOTATIONS

my @l_html_tr;
foreach my $g (@{$var->getGenes()}) {
	$vp->{gene} = $g;
	$vp->{transcripts} = $vp->{hgenes}->{$g->{id}}->{tr};
	$print_html->variant($vp);
	my $html_tr;
	$html_tr .= "<div style='border: double 1px black;background-color:#7FAB5E;color:white;text-align:left;font-size:13px;>'><span style='padding-left:15px;'><b>".$g->external_name().'</b>  <i>'.$g->id().'</i>'.'</span></div>'; 
	$html_tr .= $print_html->transcripts();
	
	#$html_tr =~ s/max-height:500px;/height:350px;/;
	$html_tr =~ s/max-height:500px;//;
	$html_tr =~ s/max-height:150px;//g;
	push(@l_html_tr, $html_tr);
}
$h->{html_transcripts} = join("<br>", @l_html_tr);

### DEJAVU

my $h_dv = $var->getChromosome->rocks_dejavu->dejavu($var->rocksdb_id);

my (@list_parquet, $h_projects_ids_names);
foreach my $project_id (sort keys %{$h_dv}) {
	print '.';
	my $project_name = $var->buffer->getQuery->getProjectNameFromId($project_id);
	$h_projects_ids_names->{$project_id} = $project_name;
	my $parquet = $var->buffer->dejavu_parquet_dir().'/'.$project_name.'.'.$project_id.'.parquet';
	push(@list_parquet, "'".$parquet."'") if -e $parquet;
}
my ($h_projects_patients, $h_gnomadid) = get_from_duckdb_project_patients_infos($var, \@list_parquet);
$var->getChromosome->rocks_dejavu->close();
#my $html_dv = qq{<table style="width:100%;border: solid 1px black;">};

my ($last_dejavu, $locus_hg38, $locus_hg19);
if ($release eq 'HG38') {
	$h->{span_locus_hg19} = 'chr'.$h_gnomadid->{$var->gnomad_id()}->{chr19}.':'.($h_gnomadid->{$var->gnomad_id()}->{pos19});
	$h->{span_locus_hg38} = 'chr'.$var->getChromosome->id().':'.$var->start();
	$locus_hg38 = $var->getChromosome->id().':'.($var->start() - 100).'-'.($var->start() + 100);
	$locus_hg19 = $h_gnomadid->{$var->gnomad_id()}->{chr19}.':'.($h_gnomadid->{$var->gnomad_id()}->{pos19} - 100).'-'.($h_gnomadid->{$var->gnomad_id()}->{pos19} + 100);
}
else {
	$h->{span_locus_hg19} = 'chr'.$var->getChromosome->id().':'.$var->start();
	$h->{span_locus_hg38} = 'chr'.$h_gnomadid->{$var->gnomad_id()}->{chr38}.':'.($h_gnomadid->{$var->gnomad_id()}->{pos38});
	$locus_hg38 = $h_gnomadid->{$var->gnomad_id()}->{chr38}.':'.($h_gnomadid->{$var->gnomad_id()}->{pos38} - 100).'-'.($h_gnomadid->{$var->gnomad_id()}->{pos38} + 100);
	$locus_hg19 = $var->getChromosome->id().':'.($var->start() - 100).'-'.($var->start() + 100);
}

my $type_control = "select";
$type_control = "input" if (scalar keys %{$h_projects_patients} > 100);



my $html_dv;
$html_dv .= "<div style='border: double 1px black;background-color:#5E8AAB;color:white;text-align:left;font-size:13px;>'><span style='padding-left:15px;'><b>DejaVu Project(s) / Sample(s)</b></span></div>";
$html_dv .= qq{<table id="table_dejavu" data-sort-name="project" data-sort-order="desc" data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-virtual-scroll="true" data-pagination-v-align="bottom" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="10" data-page-list="[10, 25, 50, 100, 200, 300]" data-resizable='true' class="table table-striped table-condensed table-bordered table-hover table-mybordered" style="height:400px;overflow-y: scroll;width:100%;vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;line-height: 25px;min-height: 25px;height: 25px;box-shadow: 3px 3px 5px #555;">};
$html_dv .= "<thead>";
$html_dv .= "<tr>";
$html_dv .= qq{<th data-field="release" data-filter-control="$type_control" data-sortable="true">Release</th>};
$html_dv .= qq{<th data-field="phenotypes" data-filter-control="$type_control" data-sortable="true">Phenotype(s)</th>};
$html_dv .= qq{<th data-field="project" data-filter-control="$type_control" data-sortable="true">Project Name</th>};
$html_dv .= qq{<th data-field="users" data-filter-control="input" data-sortable="true">Users</th>};
$html_dv .= qq{<th data-field="family" data-filter-control="input" data-sortable="true">Family Name</th>};
$html_dv .= qq{<th data-field="patient" data-filter-control="input" data-sortable="true">Patient Name</th>};
$html_dv .= qq{<th data-field="status" data-filter-control="input" data-sortable="true">Status</th>};

$html_dv .= qq{<th data-field="heho" data-filter-control="select" data-sortable="true">He / Ho</th>};
$html_dv .= qq{<th data-field="ac" data-filter-control="input" data-sortable="true">Allele Count</th>};
$html_dv .= qq{<th data-field="ratio" data-filter-control="input" data-sortable="true">Ratio (%)</th>};
$html_dv .= qq{<th data-field="dp" data-filter-control="input" data-sortable="true">DP</th>};
$html_dv .= qq{<th data-field="model" data-filter-control="$type_control" data-sortable="true">Transmission</th>};
$html_dv .= qq{<th data-field="igv">IGV</th>};
$html_dv .= "</tr>";
$html_dv .= "</thead>";
$html_dv .= "<tbody>";
foreach my $project_name (sort {$b <=> $a} keys %{$h_projects_patients}) {
	my $users = get_list_mails($buffer, $project_name);
	print '.';
	foreach my $patient_name (sort keys %{$h_projects_patients->{$project_name}}) {
		my $pat_url = $h_projects_patients->{$project_name}->{$patient_name}->{align_url};
		
		my $this_release = $h_projects_patients->{$project_name}->{$patient_name}->{release};
		my $this_locus = $locus_hg38;
		$this_locus = $locus_hg19 if ($this_release =~ /HG19/);
		
		$html_dv .= qq{<tr>};
		$html_dv .= qq{<td>}.$this_release.qq{</td>};
		$html_dv .= qq{<td>}.$h_projects_patients->{$project_name}->{$patient_name}->{phenotypes}.qq{</td>};
		$html_dv .= qq{<td>}.$project_name.qq{</td>};
		$html_dv .= qq{<td><button style="color:black;" onclick="get_popup_users('$users');">Users</button></td>};
		$html_dv .= qq{<td>}.$h_projects_patients->{$project_name}->{$patient_name}->{family}.qq{</td>};
		$html_dv .= qq{<td>}.$patient_name.qq{</td>};
		$html_dv .= qq{<td>}.$h_projects_patients->{$project_name}->{$patient_name}->{small_icon}.qq{</td>};
		$html_dv .= qq{<td>}.$h_projects_patients->{$project_name}->{$patient_name}->{heho}.qq{</td>};
		$html_dv .= qq{<td>}.$h_projects_patients->{$project_name}->{$patient_name}->{ac}.qq{</td>};
		$html_dv .= qq{<td>}.$h_projects_patients->{$project_name}->{$patient_name}->{ratio}.qq{</td>};
		$html_dv .= qq{<td>}.$h_projects_patients->{$project_name}->{$patient_name}->{dp}.qq{</td>};
		$html_dv .= qq{<td>}.$h_projects_patients->{$project_name}->{$patient_name}->{model}.qq{</td>};
		
		my $cmd = qq{view_web_igv_bam_simple('div_igv','$this_locus','$pat_url','$patient_name')};
		my $b_igv = qq{<button class='igvIcon2' onClick="$cmd"></button>};
		$html_dv .= qq{<td>}.$b_igv.qq{</td>};
		$html_dv .= qq{</tr>};
		
		my @last_igv;
		push(@last_igv, $project_name);
		push(@last_igv, $project_name.' - '.$patient_name);
		push(@last_igv, $pat_url);
		push(@last_igv, $this_locus);
		$last_dejavu = join(';', @last_igv);
	}
}

$html_dv .= "</tbody>";
$html_dv .= qq{</table>};
$h->{html_dejavu} = $html_dv;
$h->{last_dejavu} = $last_dejavu;




my $json_encode = encode_json $h;
print ".\",";
$json_encode =~ s/{//;
print $json_encode;

exit(0);


sub get_list_mails {
	my ($buffer, $project_name) = @_;
	my $h_proj = $buffer->getQuery->getProjectByName($project_name);
	my $proj_id = $h_proj->{id};
	my @lUsers;
	foreach my $h (@{ $buffer->getQuery()->getOwnerProject( $proj_id ) }) { push( @lUsers, $h->{email} ); }
	foreach my $g (@{ $buffer->getQuery()->getGroups( $proj_id ) } ){ push( @lUsers, $g ); }
	return join("<br>", sort @lUsers);
}

sub get_from_duckdb_project_patients_infos {
	my ($var, $list_files) = @_;
	return if scalar(@$list_files) == 0;
	my ($h_gnomadid, $h_projects_patients, $h_by_proj);
	my $iter = natatime(120, @$list_files);
	while( my @tmp = $iter->() ){
		print '|';
		my $sql = "PRAGMA threads=6; SELECT project,chr38,chr19,pos38,pos19,he,allele,patients,dp_ratios FROM read_parquet([".join(', ', @tmp)."])";
		my ($posVar, $altVar) = split('!', $var->rocksdb_id());
		if ($var->getProject->current_genome_version() eq 'HG38') {
			$sql .= " WHERE chr38='".$var->getChromosome->id()."' and pos38=$posVar;" ;
		}
		else {
			$sql .= " WHERE chr19='".$var->getChromosome->id()."' and pos19=$posVar;" ;
		}
		my $duckdb = $buffer->software('duckdb');
		my $cmd = qq{set +H | $duckdb -json -c "$sql"};
		my $json_duckdb = `$cmd`;
		if ($json_duckdb) {
			print '|';
			my $decode = decode_json $json_duckdb;
			MCE::Loop->init(
				max_workers => 6,
				chunk_size => 10,
				gather => sub {
					my ($data) = @_;
					foreach my $proj_name (keys %{$data->{h_by_proj}}) {
						$h_by_proj->{$proj_name} = $data->{h_by_proj}->{$proj_name};
					}
					foreach my $proj_name (keys %{$data->{h_projects_patients}}) {
						$h_projects_patients->{$proj_name} = $data->{h_projects_patients}->{$proj_name};
					}
					foreach my $vid (keys %{$data->{h_gnomadid}}) {
						$h_gnomadid->{$vid} = $data->{h_gnomadid}->{$vid};
					}
				}
			);
			mce_loop {
				my ($mce, $chunk_ref, $chunk_id) = @_;
				my $hres;
				foreach my $h (@$chunk_ref) {
					if ($var->getProject->current_genome_version() eq 'HG38') {
						next if $h->{'pos38'} ne int($posVar);
					}
					else {
						next if $h->{'pos19'} ne int($posVar);
					}
					my $var_all = $h->{'allele'};
					$var_all =~ s/\+//;
					if (not $var->isDeletion) {
						next if $var_all ne $var->var_allele();
					}
					if ($var->getProject->current_genome_version() eq 'HG38') {
						$hres->{h_gnomadid}->{$var->gnomad_id()}->{chr19} = $h->{chr19};
						$hres->{h_gnomadid}->{$var->gnomad_id()}->{pos19} = $h->{pos19};
						
					}
					else {
						$hres->{h_gnomadid}->{$var->gnomad_id()}->{chr38} = $h->{chr38};
						$hres->{h_gnomadid}->{$var->gnomad_id()}->{pos38} = $h->{pos38};
					}
					$hres->{h_gnomadid}->{$var->gnomad_id()}->{ref_all} = $var->ref_allele();
					$hres->{h_gnomadid}->{$var->gnomad_id()}->{var_all} = $var->var_allele();
					my $project_id = $h->{project};
					my $project_name = $h_projects_ids_names->{$project_id};
					$hres->{h_by_proj}->{$project_name} = $h;
					$hres->{h_projects_patients}->{$project_name} = get_table_project_patients_infos($project_name, $h);
				}
				MCE->gather($hres);
			} $decode;	
			MCE::Loop->finish();
		}
	}
	return ($h_projects_patients, $h_gnomadid);
}

sub get_table_project_patients_infos {
	my ($project_name, $hash) = @_;
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
			my $ratio = '?';
			$ratio = ($info / $h_infos_patients->{$nb_pat}->{dp}) * 100 if ($h_infos_patients->{$nb_pat}->{dp});
			$h_infos_patients->{$nb_pat}->{ac} = $info;
			my $text_ratio;
			if ($ratio and $ratio eq '?') { $text_ratio .= '?'; }
			elsif ($ratio) { $text_ratio = int($ratio); }
			$h_infos_patients->{$nb_pat}->{ratio} = $text_ratio;
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
	my $hres;
	my $b = new GBuffer;
	my $p = $b->newProject( -name => $project_name);
	my $this_release = $p->annotation_genome_version();
	$this_release = 'HG19' if ($this_release =~ /HG19/);
	$this_release = 'HG38' if ($this_release =~ /HG38/);
	
	my $pr_phenotypes = join(', ', sort @{$p->phenotypes});
	my @lPat = @{$p->getPatients()};
	return undef if scalar(@lPat) == 0;
	my $found_healthy_patient;
	my $found_he;
	
	foreach my $pat (@lPat) {
		next if not exists $h_tmp_pat->{$pat->id};
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{name} = $pat->name;
		my $icon = $pat->small_icon();
		$icon =~ s/"/'/g;
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{phenotypes} = $pr_phenotypes;
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{release} = $this_release;
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{status} = $icon;
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{sex} = '-';
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{sex} = 'male' if $pat->sex() eq '1';
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{sex} = 'female' if $pat->sex() eq '2';
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{status_txt} = '-';
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{status_txt} = 'healthy' if $pat->status() eq '1';
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{status_txt} = 'ill' if $pat->status() eq '2';
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{small_icon} = $pat->small_icon();
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{family} = $pat->getFamily->name();
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{align_url} = $pat->alignmentUrl();
		
		my $is_heho;
		if (int($h_tmp_pat->{$pat->id}) <= $nb_he) {
			$is_heho = 'He';
			$found_he = 1;
		}
		else { $is_heho = 'Ho'; }
		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{heho} = $is_heho;
	}
	foreach my $id (keys %{$h_infos_patients}) {
		my $pat_name = $h_infos_patients->{$id}->{name};
		next if not $pat_name;
		$hres->{$pat_name} = $h_infos_patients->{$id};
	}
	$p = undef;
	$b =undef;
	return $hres;
}
