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

my $io = IO::Handle->new();
$io->autoflush(1);

my $cgi = new CGI();
my $update = $cgi->param('update');
my $user_name = $cgi->param('user_name');
my $pwd = $cgi->param('pwd');
my $first_elem = $cgi->param('first_elem');
my $elem_per_page = $cgi->param('elem_per_page');
my $max_dejavu = $cgi->param('max_dejavu');
my $max_gnomad = $cgi->param('max_gnomadac');
my $max_gnomad_ho = $cgi->param('max_gnomadacho');
my $only_clinvar = $cgi->param('only_clinvar');
my $release = $cgi->param('release');
my $force_db_annotion = $cgi->param('force_db_annot');
my $multi_release = $cgi->param('multi_release');
my $print = $cgi->param('print');
my $view_others = $cgi->param('view_others');
my $debug = $cgi->param('debug');
my $use_fork = $cgi->param('fork');

$release = $force_db_annotion if ($force_db_annotion);

my $buffer = new GBuffer;
if ($update) {
	$user_name = 'all';
	$pwd = 'all';
	my $dir_out = $buffer->get_polybtf_path($release);
	unless (-d $dir_out) {
		my $cmd = "mkdir $dir_out; chmod 777 $dir_out;";
		`$cmd`;
	}
	$dir_out .= '/users/';
	unless (-d $dir_out) {
		my $cmd = "mkdir $dir_out; chmod 777 $dir_out;";
		`$cmd`;
	}
}
$user_name = lc($user_name);

my $h_polybtf_infos = $buffer->polybtf_infos();
my $date_now = $h_polybtf_infos->{date_now};

unless ($release) {
	$release = $buffer->get_polybtf_default_release();
	if ($release =~ /,/) {
		$multi_release = $release;
	}
}

my @lparams;
push(@lparams, $user_name) if ($user_name);
push(@lparams, $max_dejavu) if ($max_dejavu);
push(@lparams, $max_gnomad) if ($max_gnomad);
push(@lparams, $max_gnomad_ho) if ($max_gnomad_ho);
push(@lparams, $only_clinvar) if ($only_clinvar);
push(@lparams, $view_others) if ($view_others);
my $cache_params = join(";",@lparams).".table";

my $no_cache;
if ($print) { $no_cache = $buffer->get_lmdb_cache_btf_view($release, $user_name); }
else { $no_cache = $buffer->get_lmdb_cache_btf_resume($release, $user_name); }
my $html = $no_cache->get_cache($cache_params);
$no_cache->close();
if ($html) {
	my $h_resume = decode_json $html;
	my $h_update = $buffer->polybtf_infos();
	$h_resume->{date_last_days} = $h_update->{date_last_days};
	my $json = encode_json $h_resume;
	print $cgi->header('text/json-comment-filtered');
	print $json;
	exit(0);
}

if ($print) { $no_cache = $buffer->get_lmdb_cache_btf_view($release, $user_name, "w"); }
else { $no_cache = $buffer->get_lmdb_cache_btf_resume($release, $user_name, "w"); }

if (not $update and not $release and not $multi_release and $print) {
	if (exists $buffer->config->{polybtf_default_releases}->{default_muti}) {
		$multi_release = $buffer->config->{polybtf_default_releases}->{default_muti};
	}
}

my $hAllVAlidations;
foreach my $g_v_ids (keys %{$buffer->validations_query->getAllValidations(1)}) {
	my @lHash_val = @{$buffer->validations_query->getAllValidations(1)->{$g_v_ids}};
	my $h_last = $lHash_val[0];
	$hAllVAlidations->{$g_v_ids}->{last_status} = $h_last->{term};
	$hAllVAlidations->{$g_v_ids}->{last_gene} = $h_last->{gene_name};
	$hAllVAlidations->{$g_v_ids}->{last_date} = $h_last->{modification_date};
	$hAllVAlidations->{$g_v_ids}->{last_user} = $h_last->{user_name};
	$hAllVAlidations->{$g_v_ids}->{last_project_name} = $h_last->{project_name};
	$hAllVAlidations->{$g_v_ids}->{last_patient_name} = $h_last->{sample_name};
	$hAllVAlidations->{$g_v_ids}->{last_patient_id} = $h_last->{sample_id};
	foreach my $h_val (@lHash_val) {
		$hAllVAlidations->{$g_v_ids}->{patients_ids}->{$h_val->{sample_id}}++;
		unless (exists $hAllVAlidations->{$g_v_ids}->{projects_patients}->{$h_val->{project_name}.'_'.$h_val->{sample_name}}) {
			$hAllVAlidations->{$g_v_ids}->{projects_patients}->{$h_val->{project_name}.'_'.$h_val->{sample_name}} = $h_val->{term};
		}
	}
}

my $can_use_hgmd = $buffer->hasHgmdAccess($user_name);
$only_clinvar = 1 unless ($can_use_hgmd);
if ($user_name eq 'all' and $pwd eq 'all') {
	$only_clinvar = undef;
	$can_use_hgmd = 1;
}

$first_elem = 1 unless ($first_elem);
$elem_per_page = 500 unless ($elem_per_page);
#$elem_per_page = 200 unless ($elem_per_page);
#$elem_per_page = 30 unless ($elem_per_page);

$max_dejavu = 20 unless ($max_dejavu);
$max_gnomad = 30 unless ($max_gnomad);
$max_gnomad_ho = 10 unless ($max_gnomad_ho);

my $value_red = 14;
my $value_coral = 12;
my $value_orange = 8;
my $value_yellow = 5;


my $color_red = '#CE0000';
my $color_coral = 'coral';
my $color_orange = '#EFB73E';
my $color_yellow = 'yellow';

my $dir_polybtf_path = $buffer->get_polybtf_path($release);
unless (-d $dir_polybtf_path) {
	`mkdir $dir_polybtf_path`;
	`chmod 777 $dir_polybtf_path`;
	`mkdir $dir_polybtf_path/users/`;
	`chmod 777 $dir_polybtf_path/users/`;
}

warn "\n--- Before get_hash_by_user ---\n" if ($debug);
my $h_my_projects = get_hash_users_projects($user_name, $pwd);

my @list_other_releases;
if ($multi_release) {
	@list_other_releases = reverse sort {$a <=> $b} split(',', $multi_release);
	$release = shift @list_other_releases;
}
die unless ($release);

print $cgi->header('text/json-comment-filtered');
print "{\"progress\":\".";


my @headers_validations = ("varsome","igv","alamut","var_name","trio","gnomad","deja_vu","table_validation","table_transcript","validation_select");
#my @headers_validations = ("varsome","igv","alamut","var_name","trio","gnomad","deja_vu","table_validation","table_transcript");
my @header_transcripts = ("consequence","enst","nm","ccds","appris","exon","nomenclature","codons","codons_AA", "polyphen","sift","ncboost","cadd","revel","dbscsnv");

exit(0) unless ($h_my_projects);

my ($hFound, $hResProj) = get_hash_by_user($user_name, $pwd, $h_my_projects, $release);

my $hAllVarIds = get_hash_new_hgmd_clinvar($h_my_projects);


#warn Dumper $hAllVarIds; die;

my $hResProjGlobal;
foreach my $score (keys %$hResProj) {
	foreach my $project_name (keys %{$hResProj->{$score}}) {
		$hResProjGlobal->{$score}->{$project_name}->{$release} = undef;
	}
}

if ($multi_release) {
	foreach my $other_release (@list_other_releases) {
		my ($hFound_2, $hResProj_2) = get_hash_by_user($user_name, $pwd, $h_my_projects, $other_release);
		foreach my $score_2 (keys %$hResProj_2) {
			foreach my $project_name (keys %{$hResProj_2->{$score_2}}) {
				$hResProjGlobal->{$score_2}->{$project_name}->{$other_release} = undef;
			}
		}
		my $hAllVarIds_2 = get_hash_new_hgmd_clinvar($h_my_projects, $other_release);
		foreach my $var_id (keys %$hAllVarIds_2) {
			if (exists $hAllVarIds->{$var_id}) {
				foreach my $type (keys %{$hAllVarIds_2->{$var_id}}) {
					$hAllVarIds->{$var_id}->{$type} = undef;
				}
			}
			else { $hAllVarIds->{$var_id} = $hAllVarIds_2->{$var_id}; }
		}
	}
}

my $hVarFoundMyProjects;

my @lProJName = reverse sort keys %$h_my_projects;

my $fsize = "font-size:10px";
my $color = "#9BC8A5";
my $class;
$class->{rowspan} -= 1;
$class->{rowspan} = 1 if $class->{rowspan} <=0;
$class->{style} = "min-width:10%;padding:1px";

my $out = "";
$out .= "<div class='row' hidden>";
$out .= "<div class='col-sm-2'></div>";
$out .= "<div class='col-sm-8'>";
$out .= "<center>";
$out .= "<nav aria-label='Page navigation example'>";
$out .= "<ul class='pagination'>";


my $max_projects = scalar(@lProJName);
my $j = 1;
my $nb_page = 0;
my ($h_page, $h_page_infos, $h_part, $last_page);
while ($j < $max_projects) {
	$nb_page++;
	$h_page->{$nb_page} = 'check_new_hgmd_clinvar('.$j.')';
	$h_page_infos->{$nb_page} = 'From '.$lProJName[$j];
	$h_part->{$j} = $nb_page;
	$j += $elem_per_page;
	$h_page_infos->{$nb_page} .= ' to '.$lProJName[$j];
	$last_page = $nb_page;
}

my $this_page = $h_part->{$first_elem};
my $previous_page = 0;
if ($this_page > 0) { $previous_page = $this_page - 1; }
my $next_page = $this_page + 1;

if ($this_page > 1) {
	$out .= "<li class='page-item'><a class='page-link' href='#' onclick='".$h_page->{$this_page - 1}."'>Previous</a></li>";
}

my ($moins_3, $moins_2, $moins_1, $middle, $plus_3, $plus_2, $plus_1);
if ($this_page < 8) {
	$moins_3 = 1;
	$moins_2 = 2;
	$moins_1 = 3;
	$middle = 4;
	$plus_1 = 5;
	$plus_2 = 6;
	$plus_3 = 7;
}
else {
	$moins_3 = $this_page - 3;
	$moins_2 = $this_page - 2;
	$moins_1 = $this_page - 1;
	$middle = $this_page;
	$plus_1 = $this_page + 1;
	if ($plus_1 > $last_page) { $plus_1 = undef; }
	$plus_2 = $this_page + 2;
	if ($plus_2 > $last_page) { $plus_2 = undef; }
	$plus_3 = $this_page + 3;
	if ($plus_3 > $last_page) { $plus_3 = undef; }
}

if ($moins_3 > 1) {
	$out .= "<li class='page-item'><a class='page-link' href='#' onclick='".$h_page->{1}."'>1..</a></li>";
}

if ($moins_3 == $this_page) { $out .= "<li class='page-item active'><a class='page-link' href='#' onclick='".$h_page->{$moins_3}."'><span id='span_moins_3'>".$moins_3."</span></a></li>"; }
elsif ($moins_3 and $moins_3 < $last_page) { $out .= "<li class='page-item'><a class='page-link' href='#' onclick='".$h_page->{$moins_3}."'><span id='span_moins_3'>".$moins_3."</span></a></li>"; }

if ($moins_2 == $this_page) { $out .= "<li class='page-item active'><a class='page-link' href='#' onclick='".$h_page->{$moins_2}."'>".$moins_2."</a></li>"; }
elsif ($moins_2 and $moins_2 <= $last_page) { $out .= "<li class='page-item'><a class='page-link' href='#' onclick='".$h_page->{$moins_2}."'>".$moins_2."</a></li>"; }

if ($moins_1 == $this_page){ $out .= "<li class='page-item active'><a class='page-link' href='#' onclick='".$h_page->{$moins_1}."'>".$moins_1."</a></li>"; }
elsif ($moins_1 and $moins_1 <= $last_page) { $out .= "<li class='page-item'><a class='page-link' href='#' onclick='".$h_page->{$moins_1}."'>".$moins_1."</a></li>"; }

if ($middle == $this_page) { $out .= "<li class='page-item active'><a class='page-link' href='#' onclick='".$h_page->{$middle}."'>".$middle."</a></li>"; }
elsif ($middle and $middle <= $last_page) { $out .= "<li class='page-item'><a class='page-link' href='#' onclick='".$h_page->{$middle}."'>".$middle."</a></li>"; }

my $has_last_button = 1;
if ($plus_1 and $plus_1 == $this_page) { $out .= "<li class='page-item active'><a class='page-link' href='#' onclick='".$h_page->{$plus_1}."'>".$plus_1."</a></li>"; }
elsif ($plus_1 and $plus_1 <= $last_page) { $out .= "<li class='page-item'><a class='page-link' href='#' onclick='".$h_page->{$plus_1}."'>".$plus_1."</a></li>"; }
else { $has_last_button = undef; }

if ($plus_2 and $plus_2 == $this_page) { $out .= "<li class='page-item active'><a class='page-link' href='#' onclick='".$h_page->{$plus_2}."'>".$plus_2."</a></li>"; }
elsif ($plus_2 and $plus_2 <= $last_page) { $out .= "<li class='page-item'><a class='page-link' href='#' onclick='".$h_page->{$plus_2}."'>".$plus_2."</a></li>"; }
else { $has_last_button = undef; }


if ($plus_3 and $plus_3 == $this_page) { $out .= "<li class='page-item active'><a class='page-link' href='#' onclick='".$h_page->{$plus_3}."'>".$plus_3."</a></li>"; }
elsif ($plus_3 and $plus_3 <= $last_page) { $out .= "<li class='page-item'><a class='page-link' href='#' onclick='".$h_page->{$plus_3}."'>".$plus_3."</a></li>"; }

if ($plus_3 >= $last_page) { $has_last_button = undef; }


if ($has_last_button) { $out .= "<li class='page-item'><a class='page-link' href='#' onclick='".$h_page->{$last_page}."'>..".$last_page."</a></li>"; }

if ($next_page <= $last_page) {
	$out .= "<li class='page-item'><a class='page-link' href='#' onclick='".$h_page->{$this_page + 1}."'>Next</a></li>";
}

$out .= "</ul>";
$out .= "</nav>";
$out .= "</center>";
$out .= "</div>";
$out .= "<div class='col-sm-2'></div>";

$out .= "<div class='btn-group mr-2' role='group' aria-label='First group'>";
$out .= "<button type='button' class='btn btn-danger' onClick='select_tr(\"tr_project\", \"3\")'>High</button>";
$out .= "<button type='button' class='btn btn-warning' onClick='select_tr(\"tr_project\", \"2\")'>Medium</button>";
$out .= "<button type='button' class='btn ' onClick='select_tr(\"tr_project\", \"1\")'>Others</button>";
$out .= "</div>";
$out .= "</div>";

$out .= "<div class='row'>";
$out .= "<div class='col-sm-12'>";
$out .= "<table id='table_projects' class='table' style='padding:5px;font-size:13px;'>";
$out .= "<tbody>";

my $last_elem = $first_elem + $elem_per_page;
my $first_project_done;
my $nb_proj = 0;
my $nb_proj_no_strict_new_genes = 0;
my $out2;
my $hCollapseIds;

#warn Dumper $hResProjGlobal; die;

foreach my $score (reverse sort {$a <=> $b} keys %{$hResProjGlobal}) {
	next if ($score eq '');
	foreach my $proj_name (sort keys %{$hResProjGlobal->{$score}}) {
		foreach my $release (sort keys %{$hResProjGlobal->{$score}->{$proj_name}}) {
			print '.';
			$nb_proj++;
			
			next if ($nb_proj < $first_elem);
			last if ($nb_proj > $last_elem);
			
			my $dir = $dir_polybtf_path.'/../'.$release;
			my $dir_proj = $dir.'/'.$proj_name;
			my $json_file = $dir_proj.'/'.$proj_name.'_new_public_db.json';
			
			open (FILE, $json_file);
			my $json_encode = <FILE>;
			close (FILE);
			my $hProject = decode_json $json_encode;
			
			my $json_file_resume = $dir_proj.'/'.$proj_name.'_new_public_db.resume.json';
			open (FILE, $json_file_resume);
			my $json_encode_resume = <FILE>;
			close (FILE);
			my $hResResume = decode_json $json_encode_resume;
			
			my $h_scaled_score = $hProject->{res_by_score};
			my $hgenes = $hProject->{hgenes};
			my ($nb_red, $nb_coral, $nb_orange, $nb_yellow, $nb_grey);
			if ($can_use_hgmd) {
				$nb_red = $hResResume->{nb_red};
				$nb_coral = $hResResume->{nb_coral};
				$nb_orange = $hResResume->{nb_orange};
				$nb_yellow = $hResResume->{nb_yellow};
				$nb_grey = $hResResume->{nb_grey};
			}
			else {
				$nb_red = $hResResume->{nb_red_no_hgmd};
				$nb_coral = $hResResume->{nb_coral_no_hgmd};
				$nb_orange = $hResResume->{nb_orange_no_hgmd};
				$nb_yellow = $hResResume->{nb_yellow_no_hgmd};
				$nb_grey = $hResResume->{nb_grey_no_hgmd};
			}
			my $max_alert_project = $hResResume->{max_alert};
			
			my $collapse_id = 'collapse_'.$proj_name;
			my $tr_id = 'tr_'.$proj_name;
			$collapse_id .= '_'.$release;
			$collapse_id =~ s/\./_/;
			$tr_id .= '_'.$release;
			$tr_id =~ s/\./_/;
			$hCollapseIds->{$collapse_id} = undef;
			
			my @lScores = sort {$a <=> $b} keys %{$h_scaled_score};
			
			if ($only_clinvar) {
				my $has_new_clinvar;
				foreach my $scale_score (@lScores) {
					foreach my $gene_id (keys %{$h_scaled_score->{$scale_score}}) {
						my $has_new_clinvar_gene;
						my @lVarIds = keys %{$h_scaled_score->{$scale_score}->{$gene_id}};
						foreach my $var_id (@lVarIds) {
							if (exists $hAllVarIds->{$var_id}->{clinvar}) {
								$has_new_clinvar = 1;
								$has_new_clinvar_gene = 1;
							}
	#						else {
	#							my $hvariation = $h_scaled_score->{$scale_score}->{$gene_id}->{$var_id};
	#							if ($hvariation->{scaled_score}->{$gene_id} >= $value_red) { $nb_red--; }
	#							elsif ($hvariation->{scaled_score}->{$gene_id} >= $value_coral) { $nb_coral--; }
	#							elsif ($hvariation->{scaled_score}->{$gene_id} >= $value_orange) { $nb_orange--; }
	#							elsif ($hvariation->{scaled_score}->{$gene_id} >= $value_yellow) { $nb_yellow--; }
	#							else { $nb_grey--; }
	#						}
						}
						delete $h_scaled_score->{$scale_score}->{$gene_id} unless ($has_new_clinvar_gene);
					}
				}
				next unless ($has_new_clinvar);
			}
			
			#TODO: here
			my $h_pass_var_local;
			my $has_variants = 0;
			if ($print) {
				foreach my $scale_score (reverse @lScores) {
					foreach my $gene_id (sort keys %{$h_scaled_score->{$scale_score}}) {
						foreach my $var_id (sort keys %{$h_scaled_score->{$scale_score}->{$gene_id}}) {
							my $g_v_id = $gene_id.'!'.$var_id;
							my $hvariation = $h_scaled_score->{$scale_score}->{$gene_id}->{$var_id};
							
							foreach my $gene_name (keys %{$hvariation->{html}->{'table_validation_select'}}) {
								foreach my $patient_name (keys %{$hvariation->{html}->{'table_validation_select'}->{$gene_name}}) {
									if (exists $hAllVAlidations->{$g_v_id}->{projects_patients}->{$proj_name.'_'.$patient_name}) {
										$h_pass_var_local->{$var_id}++;
										next;
									}
									$has_variants = 1;
									last;
								}
							}
						}
						last if ($has_variants == 1);
					}
					last if ($has_variants == 1);
				}
				
				next unless ($has_variants);
				
				my $class_tr = "tr_project 1";
				$class_tr = " tr_project 2" if ($nb_orange > 0);
				$class_tr = " tr_project 3" if ($nb_red > 0);
				
				$out2 .= "<tr id='".$tr_id."' class='".$class_tr."' data-toggle='collapse' href='#$collapse_id' aria-expanded='false' aria-controls='$collapse_id'>";
				$out2 .= "<td>";
				
				$out2 .= "<div class='row'>";
				$out2 .= "<div class='col-sm-12' style='text-align:left;background-color:#194b7f;'>";
				$out2 .= "<div class='btn-toolbar style='padding-left:30px;' role='toolbar' aria-label='Toolbar with button groups'>";
				$out2 .= "<div class='btn-group' role='group' aria-label='First group'>";
				
				$out2 .= "<button type='button' style='min-width:50px;font-size:11px;background-color:#194b7f;' class='btn'></button>";
				
#				my @lB_none;
#				if ($nb_red > 0) {
#					$out2 .= "<button type='button' style='min-width:40px;border:solid 1px #5A6268;font-size:11px;' class='btn btn-danger'>".$nb_red."</button>";
#				}
#				else {
#					push(@lB_none, "<button type='button' style='min-width:40px;font-size:11px;background-color:#194b7f;' class='btn'></button>");
#				}
#				if ($nb_coral > 0) {
#					$out2 .= "<button type='button' style='min-width:40px;border:solid 1px #5A6268;background-color:coral;font-size:11px;' class='btn'>".$nb_coral."</button>";
#				}
#				else {
#					push(@lB_none, "<button type='button' style='min-width:40px;font-size:11px;background-color:#194b7f;' class='btn'></button>");
#				}
#				if ($nb_orange > 0) {
#					$out2 .= "<button type='button' style='min-width:40px;border:solid 1px #5A6268;font-size:11px;' class='btn btn-warning'>".$nb_orange."</button>";
#				}
#				else {
#					push(@lB_none, "<button type='button' style='min-width:40px;font-size:11px;background-color:#194b7f;' class='btn'></button>");
#				}
#				if ($nb_yellow > 0) {
#					$out2 .= "<button type='button' style='min-width:40px;border:solid 1px #5A6268;background-color:yellow;color:black;font-size:11px;' class='btn'>".$nb_yellow."</button>";
#				}
#				else {
#					push(@lB_none, "<button type='button' style='min-width:40px;font-size:11px;background-color:#194b7f;' class='btn'></button>");
#				}
#				if ($nb_grey > 0) {
#					$out2 .= "<button type='button' style='min-width:40px;border:solid 1px #5A6268;color:black;font-size:11px;' class='btn'>".$nb_grey."</button>";
#				}
#				else {
#					push(@lB_none, "<button type='button' style='min-width:40px;font-size:11px;background-color:#194b7f;' class='btn'></button>");
#				}
#				$out2 .= join('', @lB_none);
				
				if ($multi_release) {
					$out2 .= "<button type='button' style='background-color:#194b7f;color:yellow;font-size:11px;min-width:100px;max-width:100px;' class='btn'>release ".$release."</button>";
				}
				
				$out2 .= "<button type='button' style='min-width:50px;font-size:11px;background-color:#194b7f;' class='btn'></button>";
				
				$out2 .= "<button type='button' style='background-color:#194b7f ;color:white;font-size:13px;' class='btn'><b>".$hResResume->{project_name}."</b></button>";
				if (exists $hResResume->{project_phenotypes} and $hResResume->{project_phenotypes} ne '') {
					$out2 .= "<button type='button' style='background-color:#194b7f;color:orange;font-size:11px;min-width:250px;max-width:250px;' class='btn'>".$hResResume->{project_phenotypes}."</button>";
				}
				else {
					$out2 .= "<button type='button' style='background-color:#194b7f;color:white;font-size:11px;min-width:250px;max-width:250px;' class='btn'>no phenotype</button>";
				}
				$out2 .= "<button type='button' style='background-color:#194b7f;color:white;font-size:11px;min-width:250px;max-width:500px;' class='btn'>".$hResResume->{project_description}."</button>";
				
				$out2 .= "<div class='btn-group' role='group' aria-label='First group'>";
				my $astyle = "border:solid 1px #FFFFFF;background-color:#FFFFFF";
				my $glyph = "";
				if    ($max_alert_project == 1) { $glyph = qq{<span class="glyphicon glyphicon-star-empty text-default" aria-hidden="true"></span>} }
				elsif ($max_alert_project == 2) { $glyph = qq{<span class="glyphicon  glyphicon-alert text-alert" aria-hidden="true" style="color:red"></span>} }
				elsif ($max_alert_project == 3) {
					#$glyph = qq{<img src="https://img.icons8.com/color/24/000000/treatment-plan.png" style="float:right">}
				}
				elsif ($max_alert_project == 4) {
					#$astyle = "border:solid 1px #5A6268;background-color:#CC6ED0";
					#$glyph = qq{<img src="https://img.icons8.com/color/24/000000/treatment-plan.png" style="float:right;height:20px;width:20px;">}
				}
				if ($max_alert_project >= 4) {
					#$out2 .= "<button type='button' style='min-width:46px;$astyle;' class='btn'>$glyph</button>";
				}
				$out2 .= "</center></div></div></div>";
				
				$out2 .= "</div>";
				$out2 .= "</div>";
				$out2 .= "</td>";
				$out2 .= "</tr>";
				unless ($first_project_done) {
					$out2 .= "<tr class='collapse in' id='$collapse_id'>";
					$first_project_done = $proj_name;
				}
				else {
					$out2 .= "<tr class='collapse' id='$collapse_id'>";
				}
				#$out2_header .= "<td colspan='4'>";
				$out2 .= "<td>";
				$out2 .= "<div>";
				
				$out2 .= "<div class='col-sm-2'></div>";
			}
			my ($h_var_table_html, $h_var_table_html_printed);
			foreach my $scale_score (reverse @lScores) {
				foreach my $gene_id (sort keys %{$h_scaled_score->{$scale_score}}) {
					my @lVarIds = sort keys %{$h_scaled_score->{$scale_score}->{$gene_id}};
					if ($print) {
						foreach my $var_id (@lVarIds) {
							next if (exists $h_pass_var_local->{$var_id});
							my ($last_date_validation, $last_status_validation, $last_gene_name_validation, $color_last_validation, $last_user_validation, $last_project_validation, $last_patient_validation);
							my $g_v_id = $gene_id.'!'.$var_id;
							#TODO: here
							next if (exists $hAllVAlidations->{$g_v_id} and $hAllVAlidations->{$g_v_id}->{projects}->{$proj_name});
							if (exists $hAllVAlidations->{$g_v_id}) {
								$last_date_validation = $hAllVAlidations->{$g_v_id}->{last_date};
								$last_status_validation = $hAllVAlidations->{$g_v_id}->{last_status};
								$last_gene_name_validation = $hAllVAlidations->{$g_v_id}->{last_gene};
								$last_user_validation = $hAllVAlidations->{$g_v_id}->{last_user};
								$last_project_validation = $hAllVAlidations->{$g_v_id}->{last_project_name};
								$last_patient_validation = $hAllVAlidations->{$g_v_id}->{last_patient_name};
								$color_last_validation = 'blue';
								$color_last_validation = 'green' if (lc($last_status_validation) eq 'benign' or lc($last_status_validation) eq 'likely benign');
								$color_last_validation = 'green' if (lc($last_status_validation) eq 'uncertain significance');
								$color_last_validation = 'red' if (lc($last_status_validation) eq 'likely pathogenic' or lc($last_status_validation) eq 'pathogenic');
							}
							if ($only_clinvar) {
								next unless (exists $hAllVarIds->{$var_id}->{clinvar});
							}
							my $hvariation = $h_scaled_score->{$scale_score}->{$gene_id}->{$var_id};
							
							#next if ($hvariation->{value}->{dm} == 1 and $hvariation->{value}->{clinvar_pathogenic} == 1);
							
							unless (exists $h_var_table_html->{$var_id}) {
								unless ($can_use_hgmd) {
									$hvariation->{html}->{hgmd} = qq{<span class="glyphicon glyphicon-ban-circle" aria-hidden="true" style='font-size:12px;color:black;'></span>};
								}
								foreach my $h (@headers_validations){
#									if ($h eq "trio" or "table_transcript"){
#										$class->{style} = "min-width:10%;vertical-align:middle;padding:5px;";
#									}
#									elsif ($h eq "igv" or "alamut"){
#										$class->{style} = "max-width:50px;vertical-align:middle;padding:5px;";
#									}
#									else {
#										$class->{style} = "min-width:5%;vertical-align:middle;padding:5px;";
#									}
									if ($h eq 'scaled_score') {
										my $color = 'black';
										if ($hvariation->{$h}->{$gene_id} >= $value_red) { $color = $color_red; }
										elsif ($hvariation->{$h}->{$gene_id} >= $value_coral) { $color = $color_coral; }
										elsif ($hvariation->{$h}->{$gene_id} >= $value_orange) { $color = $color_orange; }
										elsif ($hvariation->{$h}->{$gene_id} >= $value_yellow) { $color = $color_yellow; }
										my $score = sprintf("%.1f", $hvariation->{$h}->{$gene_id});
										my $b = "<button type='button' class='btn ' style='background-color:white;font-size:10px;'>";
										$b .= "<span style='color:$color;'>".$score."<span>";
										$b .= "</button>";
#										$h_var_table_html->{$var_id}->{$h} = $cgi->td($class,$b);
										$h_var_table_html->{$var_id}->{$h} = $b;
									}
									elsif ($h eq 'table_validation') {
#										if ($can_use_hgmd) { $h_var_table_html->{$var_id}->{$h} = $cgi->td($class,$hvariation->{html}->{$h}); }
#										else { $h_var_table_html->{$var_id}->{$h} = $cgi->td($class,$hvariation->{html}->{table_validation_hgmd_no_access}); }
										if ($can_use_hgmd) { $h_var_table_html->{$var_id}->{$h} = $hvariation->{html}->{$h}; }
										else { $h_var_table_html->{$var_id}->{$h} = $hvariation->{html}->{table_validation_hgmd_no_access}; }
										if ($last_date_validation and $last_status_validation) {
											$h_var_table_html->{$var_id}->{$h} .= "<center><div style='margin-top:5px;border:solid 1px $color_last_validation;'>";
											$h_var_table_html->{$var_id}->{$h} .= "<span style='margin-top:5px;color:$color_last_validation;font-size:9px;'><b>$last_project_validation - $last_patient_validation</b></span><br>"; 
											$h_var_table_html->{$var_id}->{$h} .= "<span style='margin-top:5px;color:$color_last_validation;font-size:9px;'><b><u>Last Status:</b></u> $last_status_validation (gene: $last_gene_name_validation)</span><br>"; 
											$h_var_table_html->{$var_id}->{$h} .= "<span style='color:$color_last_validation;font-size:9px;'><b><u>Last Date:</b></u> $last_date_validation</span><br>";
											$h_var_table_html->{$var_id}->{$h} .= "<span style='color:$color_last_validation;font-size:9px;'><b><u>Last User:</b></u> $last_user_validation</span></div></center>"; 
										}
									}
									elsif ($h eq 'validation_select') {
										my $html_val_sel = "<div>";
										$html_val_sel .= "<table style='font-size:11px;border:solid 1px #194b7f;'>";
										$html_val_sel .= "<thead>";
										$html_val_sel .= $cgi->start_Tr({style=>"background-color:#194b7f;color:white;"});
										$html_val_sel .= $cgi->th('<center>Gene</center>');
										$html_val_sel .= $cgi->th('<center>Patient</center>');
										$html_val_sel .= $cgi->th('<center>Validation</center>');
										$html_val_sel .= $cgi->end_Tr();
										$html_val_sel .= "</thead>";
										$html_val_sel .= "<tbody>";
										
										
#										if ($var_id eq '20_57429596_G_C') {
#											warn "\n\n\n";
#											warn Dumper $hvariation->{html}->{'table_validation_select'};
#											warn "\n\n\n";
#											die;
#										}
										
										foreach my $gene_name (keys %{$hvariation->{html}->{'table_validation_select'}}) {
											my $nbPat = scalar keys %{$hvariation->{html}->{'table_validation_select'}->{$gene_name}};
											$html_val_sel .= $cgi->start_Tr({style=>"background-color:#E9DEFF;$fsize"});
											$html_val_sel .= "<td rowspan='$nbPat'><button style='color:black;'>$gene_name</button></td>";
											foreach my $patient_name (keys %{$hvariation->{html}->{'table_validation_select'}->{$gene_name}}) {
												$html_val_sel .= "<td><button style='color:black;'>".$patient_name."</button></td>";
												if (exists $hAllVAlidations->{$g_v_id}->{projects_patients}->{$proj_name.'_'.$patient_name}) {
													my $this_status = $hAllVAlidations->{$g_v_id}->{projects_patients}->{$proj_name.'_'.$patient_name};
													my $color = 'black';
													$color = 'green' if (lc($this_status) eq 'benign' or lc($this_status) eq 'likely benign');
													$color = 'orange' if (lc($this_status) eq 'uncertain significance');
													$color = 'red' if (lc($this_status) eq 'likely pathogenic' or lc($this_status) eq 'pathogenic');
													$html_val_sel .= "<td style='min-width:100px;'><center><button style='color:white;background-color:".$color.";'>".$this_status."</button></center></td>";
												}
												else {
													$html_val_sel .= $hvariation->{html}->{'table_validation_select'}->{$gene_name}->{$patient_name};
												}
											}
										}
										$html_val_sel .= "</tbody>";
										$html_val_sel .= "</table>";
										$html_val_sel .= "</div>";
										$h_var_table_html->{$var_id}->{$h} = $html_val_sel;
									}
									else {
										$h_var_table_html->{$var_id}->{$h} = $hvariation->{html}->{$h};
									}
								}
							}
							
							if (exists $hvariation->{genes_pathogenic}->{$gene_id}) {
								my $html_table_transcript_gene;
								$html_table_transcript_gene .= $hvariation->{html}->{$gene_id}->{table_transcript};
								foreach my $h_tr (sort @{$hvariation->{genes}->{$gene_id}}) {
									my $value_cons = $h_tr->{value}->{consequence};
									my $html_cons = $h_tr->{html}->{consequence};
									$h_var_table_html->{$var_id}->{hash_table_transcript_consequences}->{$gene_id}->{$value_cons} = $html_cons;
								}
								$h_var_table_html->{$var_id}->{hash_table_transcript}->{$scale_score}->{$gene_id} = $html_table_transcript_gene;
							}
						}
					}
				}
			}
			if ($print) {
				$out2 .=  $cgi->start_div({style=>"background-color:white;"});
				$out2 .= "<table class='table' style='font-size:13px;'>";
				$out2 .= "<thead>";
				$out2 .= $cgi->start_Tr({style=>"background-color:#E9DEFF;$fsize"});
				foreach my $h (@headers_validations){
					if ($h eq 'scaled_score') { $out2 .= $cgi->th('Old_diag_score'); }
					else { $out2 .=  $cgi->th(ucfirst($h)); }
				}
				$out2 .= $cgi->end_Tr();
				$out2 .= "</thead>";
				$out2 .= "<tbody>";
			}
			foreach my $scale_score (reverse @lScores) {
				foreach my $gene_id (reverse sort {$a <=> $b} keys %{$h_scaled_score->{$scale_score}}) {
					my @lVarIds = sort keys %{$h_scaled_score->{$scale_score}->{$gene_id}};
					foreach my $var_id (@lVarIds) {
						if ($only_clinvar) {
							next unless (exists $hAllVarIds->{$var_id}->{clinvar});
						}
						$hVarFoundMyProjects->{$var_id} = undef;
						if ($print) {
							next if (exists $h_var_table_html_printed->{$var_id});
							my $class_tr;
							$class_tr->{style} = "background:white";
							$out2 .= $cgi->start_Tr($class_tr);
							foreach my $h (@headers_validations){
								if ($h eq "trio" or "table_transcript" or "table_validation"){
									$class->{style} = "max-height:155px;overflow-y:auto;vertical-align:middle;padding:5px;";
								}
								elsif ($h eq "igv" or "alamut"){
									$class->{style} = "vertical-align:middle;padding:5px;";
								}
								else {
									$class->{style} = "vertical-align:middle;padding:5px;";
								}
								if ($h eq 'table_transcript') {
									my $html_table_transcript = qq{<div style="max-height:300px;overflow-y:auto;overflow-x:hidden;">};
									$html_table_transcript .= "<table class='table' style='font-size:13px;'>";
									$html_table_transcript .= "<tbody>";
									my $hGenes_done;
									foreach my $scale_score (keys %{$h_var_table_html->{$var_id}->{hash_table_transcript}}) {
										foreach my $gene_id (keys %{$h_var_table_html->{$var_id}->{hash_table_transcript}->{$scale_score}}) {
											next if (exists $hGenes_done->{$gene_id});
											
											my $this_collapse_id = 'tr_'.$proj_name.'_'.$var_id.'_'.$gene_id;
											my $class_tr_gene->{style} = "background-color:#4A4F53;height:25px;padding:0px;white-space: nowrap;";
											$html_table_transcript .= qq{<tr style="background-color:#4A4F53;" data-toggle='collapse' href='#$this_collapse_id' aria-expanded='false' aria-controls='$this_collapse_id'>};
											$hgenes->{$gene_id}->{max_score} = $scale_score;
											my @listConsHtml;
											foreach my $cons (sort keys %{$h_var_table_html->{$var_id}->{hash_table_transcript_consequences}->{$gene_id}}) {
												push(@listConsHtml, $h_var_table_html->{$var_id}->{hash_table_transcript_consequences}->{$gene_id}->{$cons});
											}
											$html_table_transcript .= $cgi->td($class_tr_gene, update_variant_editor::panel_gene_short_button($hgenes->{$gene_id}, $proj_name, \@listConsHtml));
											$html_table_transcript .= $cgi->end_Tr();
											$html_table_transcript .= qq{<tr class='collapse' id=$this_collapse_id style="background-color:#F3F3F3;border:solid 3px #4A4F53;border-bottom:solid 2px #4A4F53;">};
											$html_table_transcript .= $cgi->td($class, $h_var_table_html->{$var_id}->{hash_table_transcript}->{$scale_score}->{$gene_id});
											$html_table_transcript .= $cgi->end_Tr();
											$html_table_transcript .= qq{<tr style="height:3px">};
											$html_table_transcript .= $cgi->td($class, "");
											$html_table_transcript .= $cgi->end_Tr();
											$hGenes_done->{$gene_id} = undef;
											
										}
									}
									$html_table_transcript .= "</tbody>";
									$html_table_transcript .= "</table>"; 
									$html_table_transcript .= "</div>";
									$out2 .= $cgi->td($class, $html_table_transcript);
								}
								else {
									$out2 .= $cgi->td($class, $h_var_table_html->{$var_id}->{$h});
								}
							}
							$out2 .= $cgi->end_Tr();
							$h_var_table_html_printed->{$var_id} = 1;
						}
					}
				}
			}
			if ($print) {
				$out2 .= "</tbody>";
				$out2 .= "</table>";
				$out2 .= "</div>";
				$out2 .= "<br>";
			}
			
			foreach my $chr_id (keys %{$buffer->{lmdb}}) {
				foreach my $db (keys %{$buffer->{lmdb}->{$chr_id}}) {
					foreach my $type (keys %{$buffer->{lmdb}->{$chr_id}->{$db}}) {
						$buffer->{lmdb}->{$chr_id}->{$db}->{$type}->close();
						delete $buffer->{lmdb}->{$chr_id}->{$db}->{$type};
					}
					delete $buffer->{lmdb}->{$chr_id}->{$db};
				}
				delete $buffer->{lmdb}->{$chr_id};
			}
	
			if ($out2) {
				$out2 .= "</div>";
				$out2 .= "</td>";
				$out2 .= "</tr>";
			}
		}
	}
}

if ($nb_proj_no_strict_new_genes > 0 and not $view_others) {
	$out2 .= "<tr id='tr_others_projects'>";
	$out2 .= "<td>";
	$out2 .= "<div class='row'>";
	$out2 .= "<div class='col-sm-10' style='text-align:left;'>";
	$out2 .= "<div class='btn-toolbar style='padding-left:30px;' justify-content-between' role='toolbar' aria-label='Toolbar with button groups'>";
	$out2 .= "<div class='btn-group' role='group' aria-label='First group'>";
	$out2 .= "<button type='button' style='background-color:#5A6268;color:white;' class='btn' onclick='check_new_hgmd_clinvar_others_projects()'>+ $nb_proj_no_strict_new_genes projects</button>";
	$out2 .= "</div>";
	$out2 .= "</div>";
	$out2 .= "</div>";
	$out2 .= "</div>";
	$out2 .= "</tr>";
}


$out .= $out2;
$out .= "</tbody>";
$out .= "</table>";
$out .= "</div>";
$out .= "</div>";

#die;

$buffer = new GBuffer;
my $hashRes;
my (@lReleases, $hReleasesVersions);
foreach my $release (reverse sort keys %{$buffer->config->{polybtf_releases}}) {
	my @lNews = split(',', $buffer->config->{polybtf_releases}->{$release});
	my $text = $release.'|<b><u>'.$release.'</b></u>';
	foreach my $soft_name (sort @lNews) {
		next if($soft_name eq 'hgmd' and not $can_use_hgmd);
		my ($genecode, $genes_annot_version) = split('\.', $release);
		my $soft_version = $buffer->public_data->{$genes_annot_version}->{$soft_name}->{version};
		$text .= ' - '.$soft_name.' ('.$soft_version.')';
		$hReleasesVersions->{$release}->{$soft_name} = $soft_version;
	}
	push(@lReleases, $text);
}
$hashRes->{releases_available} = join(';', @lReleases);

my ($nb_new, $current_version_text);
if ($only_clinvar) {
	if ($release and $release ne $buffer->getQuery->getCurrentGenomeProjectReleasesAnntotations()) { 
		$current_version_text = 'Clinvar: '.$hReleasesVersions->{$release}->{clinvar};
	}
	else {
		$current_version_text = 'Clinvar: '.$buffer->queryClinvarPathogenic->last_release_name();
	}
}
else {
	if ($release and $release ne $buffer->getQuery->getCurrentGenomeProjectReleasesAnntotations()) { 
		$current_version_text = 'HGMD: '.$hReleasesVersions->{$release}->{hgmd}.'<br>Clinvar: '.$hReleasesVersions->{$release}->{clinvar};
	}
	else {
		my $hgmd_version = $buffer->queryHgmd->database();
		$hgmd_version =~ s/hgmd_pro-//;
		$current_version_text = 'HGMD: '.$hgmd_version.'<br>Clinvar: '.$buffer->queryClinvarPathogenic->last_release_name();
	}
}

$hashRes->{news} = join(', ', @{$h_polybtf_infos->{news}});
$hashRes->{date_release} = $h_polybtf_infos->{date_release};
$hashRes->{date_now} = $h_polybtf_infos->{date_now};
$hashRes->{date_last_days} = $h_polybtf_infos->{date_last_days};

$nb_new = scalar(keys %$hAllVarIds);
$hashRes->{nb_new} = $nb_new;
$hashRes->{nb_new_in_projects} = scalar(keys %$hVarFoundMyProjects);
$hashRes->{max_dejavu} = $max_dejavu;
$hashRes->{max_gnomadac} = $max_gnomad;
$hashRes->{max_gnomad_ho} = $max_gnomad_ho;
$hashRes->{current_version} = $current_version_text;
$hashRes->{html_table} = $out if ($print);
$hashRes->{first_project} = $first_project_done;
if ($can_use_hgmd) {
	$hashRes->{last_hgmd_release} = $buffer->queryHgmd->database();
	$hashRes->{last_hgmd_release} =~ s/hgmd_pro-//; 
}
else { $hashRes->{last_hgmd_release} = qq{<span class="glyphicon glyphicon-ban-circle" aria-hidden="true" style='font-size:12px;color:black;'></span>}; }
$hashRes->{last_clinvar_release} = $buffer->queryClinvarPathogenic->last_release_name();

if ($release and $release ne $buffer->getQuery->getCurrentGenomeProjectReleasesAnntotations()) { $hashRes->{releases_used} = $release; }
else { $hashRes->{releases_used} = $buffer->getQuery->getCurrentGenomeProjectReleasesAnntotations(); }
if ($multi_release) { $hashRes->{releases_used} .= ", ".join(', ', @list_other_releases); }

my $json_encode = encode_json $hashRes;
$no_cache->put_cache($cache_params, $json_encode, 2400); 
$no_cache->close();

print ".\",";
$json_encode =~ s/{//;
print $json_encode;




exit(0);



sub get_hash_by_user {
	my ($user_name, $pwd, $h_my_projects, $release) = @_;
	warn "--- ENTER get_hash_by_user ---\n\n" if ($debug);
	warn $release if ($debug);
	my ($hFound, $hResProj);
	my $dir;
	if ($release) { $dir = $dir_polybtf_path.'/../'.$release.'/users/'; }
	else { $dir = $dir_polybtf_path.'/users/'; }
	my $json_user_file = $dir.'/'.$user_name.'.json';
	
	warn $json_user_file if ($debug);
	
	if (-e $json_user_file) {
		warn "\n--- Exists $json_user_file ---\n" if ($debug);
		warn $release if ($debug);
		$hFound = get_hash_new_hgmd_clinvar_with_dejavu_gnomad_filters(\{}, $release);
		
		open (FILE, $json_user_file);
		my $json_encode = <FILE>;
		close (FILE);
		$hResProj = decode_json $json_encode;
	}
	else {
		warn "\n--- NOT Exists $json_user_file ---\n" if ($debug);
		warn $release if ($debug);
		($hFound, $hResProj) = get_projects_hashes_from_user($user_name, $pwd, $h_my_projects, $release);
		
		my $json_encode;
		unless ($hResProj) { $json_encode = encode_json {}; }
		else { $json_encode = encode_json $hResProj; }
		
		unless (-d $dir) {
			my $cmd1 = "mkdir $dir";
			`$cmd1`;
			my $cmd2 = "chmod 777 $dir";
			`$cmd2`;
		}
		open (FILE, ">$json_user_file");
		print FILE $json_encode;
		close (FILE);
	}
	warn "\n--- Return hFound hResProj---\n" if ($debug);
	return ($hFound, $hResProj);
}

sub get_hash_users_projects {
	my ($user_name, $pwd) = @_;
	warn "\n--- Enter get_hash_users_projects ---\n" if ($debug);
	return get_hash_all_projects() if ($user_name eq 'all' and $pwd eq 'all');
	my $h_projects;
	my @list_hash = @{$buffer->getQuery()->getProjectListForUser($user_name, $pwd)};
	foreach my $hash (@list_hash) {
		my $proj_name = $hash->{name};
		next unless ($proj_name =~ /NGS20/);
		$h_projects->{$proj_name} = $hash->{description};
	}
	return $h_projects;
}

sub get_hash_all_projects {
	warn "\n--- Enter get_hash_all_projects ---\n" if ($debug);
	my $h_projects;
	my @list_hash = @{$buffer->getQuery()->getAllProjects()};
	foreach my $hash (@list_hash) {
		my $proj_name = $hash->{name};
		next unless ($proj_name =~ /NGS20/);
		$h_projects->{$proj_name} = $hash->{description};
	}
	return $h_projects;
}

sub get_projects_hashes_from_user {
	my ($user_name, $pwd, $h_my_projects, $release) = @_;
	
	warn "\n--- Enter get_projects_hashes_from_user ---\n" if ($debug);
	warn $release if ($debug);
	my $hResProj;
#	my $h_my_projects = get_hash_users_projects($user_name, $pwd);
	my $hAllVarIds = get_hash_new_hgmd_clinvar($h_my_projects, $release);
	my $hFound = get_hash_new_hgmd_clinvar_with_dejavu_gnomad_filters($hAllVarIds, $release);
	
	my @lProJName = sort keys %$hFound;
	
	warn "\n--- Before FORK ---\n" if ($debug);
	my $fork = 1;
	$fork = $use_fork if ($use_fork);
	my $pm = new Parallel::ForkManager($fork);
	$pm->run_on_finish(
		sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
			my $project_name = $data->{project_name};
			my $max_score = $data->{max_score};
			$hResProj->{$max_score}->{$project_name} = undef;
		}
	 );
	
	my $i = 0;
	my $last = $first_elem + $elem_per_page;
	foreach my $proj_name (@lProJName) {
		next unless (exists $h_my_projects->{$proj_name});
		$i++;
		next if ($i < $first_elem);
		#last if ($i == $last);
		my $pid = $pm->start and next;
		
		my $dir = $dir_polybtf_path;
		if ($release) { $dir .= '/../'.$release.'/'; }
		my $dir_proj = $dir.'/'.$proj_name;
		my $json_file = $dir_proj.'/'.$proj_name.'_new_public_db.json';
		my $json_file_resume = $dir_proj.'/'.$proj_name.'_new_public_db.resume.json';
		if (-e $json_file and -e $json_file_resume) {
			print '.';
			open (FILE, $json_file_resume);
			my $json_encode = <FILE>;
			close (FILE);
			my $hResResume = decode_json $json_encode;
			$pm->finish(0, $hResResume);
			next;
		}
		$buffer = undef;
		$buffer = new GBuffer;
		my $project = $buffer->newProjectCache( -name => $proj_name );
		$project->cgi_object(1);
		$project->print_dot(1);
		$project->getCacheDir();
		warn 'POLYCACHE: '.$project->getCacheDir() if ($debug);
		my $last_annot_release = $buffer->get_polybtf_default_release();
		
		my ($genecode_version, $annot_version) = split('\.', $last_annot_release);
		$buffer->public_data_version($annot_version);
		
		my $max_alert_project = 0;
		my @lVar = sort keys %{$hFound->{$proj_name}};
		my $out2;
		my $out2_header;
		my ($hChrUsed, $h_scaled_score, $hgenes);
		foreach my $var_id (@lVar) {
			my ($is_dm, $is_clinvar, $is_local) = (0, 0, 0);
			$is_dm = 1 if exists $hAllVarIds->{$var_id}->{hgmd};
			$is_clinvar = 1 if exists $hAllVarIds->{$var_id}->{clinvar};
			if (exists $hAllVarIds->{$var_id}->{'local'}) {
				$is_local = 1;
				$is_local = 0 if (exists $hAllVarIds->{$var_id}->{'local'}->{$proj_name});
			}
			next if ($is_dm == 0 and $is_clinvar == 0 and $is_local == 0);
			if ($only_clinvar) {
				next unless (exists $hAllVarIds->{$var_id}->{clinvar});
			}
			my $chr_id = $hFound->{$proj_name}->{$var_id}->{chr_id};
			$hChrUsed->{$chr_id} = undef;
			my ($h_local_scaled_score_done);
			my $v;
			eval { $v = $project->getVariant($var_id); };
			if ($@) {
				last;
			}
			foreach my $g (@{$v->getGenes()}) {
				next if ($g->getVariantsVector->is_empty());
				my $h_families_done;
				foreach my $family (@{$project->getFamilies()}) {
					my @lTrios;
					foreach my $patient (@{$family->getChildrenIll()}) {
						next unless ($patient);
						next unless (exists $hFound->{$proj_name}->{$var_id}->{patients}->{$patient->name()});
						$v = $project->getVariant($var_id) unless ($v);
						delete $v->{scale_score} if (exists $v->{scale_score});
						if (exists $hAllVarIds->{$var_id}->{hgmd}) { $v->{isDM} = 1; }
						if (exists $hAllVarIds->{$var_id}->{clinvar}) { $v->{isClinvarPathogenic} = 1; }
						my $hvariation;
						eval { $hvariation = update_variant_editor::construct_hash_variant( $project, $v, undef, $patient); };
						if ($@) { next; }
						# NO CSS ADDED from others scripts, like view.pl (HGMD).
						$hvariation->{html}->{no_css_polydiag} = 1;
						$hvariation->{obj} = $v;
						my $max_score = -99;
						my $max_tr;
						foreach my $tr (@{$g->getTranscripts()}) {
							my $this_score = $v->scaledScoreVariant($tr, $patient, undef);
							if ($this_score > $max_score) {
								$max_score = $this_score;
								$max_tr = $tr;
							}
						}
						my $this_score = $max_score;
						my $tr = $max_tr;
						$hvariation->{scaled_score}->{$g->id()} = $this_score;
						$hvariation->{scaled_score_gene_used} = $g->external_name();
						$hvariation->{scaled_score_transcript_used} = $tr->id();
						unless (exists $hgenes->{$g->id}){
							$hgenes->{$g->id}->{project_name} = $proj_name;
							$hgenes->{$g->id}->{name} = $g->external_name;
							$hgenes->{$g->id}->{description} = $g->description;
							$hgenes->{$g->id}->{phenotypes} = $g->phenotypes;
							$hgenes->{$g->id}->{score_mother} = 0;
							$hgenes->{$g->id}->{score_father} = 0;
							$hgenes->{$g->id}->{score_biallelic} = -10;
							$hgenes->{$g->id}->{score} = $g->score;
							$hgenes->{$g->id}->{id} = $g->id;
							$hgenes->{$g->id}->{omim_inheritance} = $g->omim_inheritance;
							$hgenes->{$g->id}->{external_name} = $g->external_name;
							$hgenes->{$g->id}->{pLI} = $g->pLI;
							$hgenes->{$g->id}->{omim_id} = $g->omim_id;
							$hgenes->{$g->id}->{panels} = $buffer->queryPanel()->getPanelsForGeneName($g->external_name);
							$hgenes->{$g->id}->{js_id} = $proj_name."_".$g->id;
						}
						$hgenes->{$g->id}->{score_variant}->{$v->id} = $this_score;
						if ($v->isDM_for_gene($g) or $v->isClinvarPathogenic_for_gene($g) or exists $hAllVarIds->{$var_id}->{local}) {
							$hgenes->{$g->id}->{pathogenic} ++;
							$hvariation->{genes_pathogenic}->{$g->id()}++;
						}
						$hgenes->{$g->id}->{clinvar_hgmd} ++ if  $v->hgmd or $v->clinvar;
						#$hgenes->{$g->id}->{denovo_rare} ++ if $v->getGnomadAC< 10 &&  $v->isDenovoTransmission($patient->getFamily,$patient);
						my $val = $v->score_validations($g);
						if ($val){
							$hgenes->{$g->id}->{validations} = $val->{validation};
						}
						push (@{$hgenes->{$g->id}->{variants}}, $v->id);
						if ($patient->isChild && $patient->getFamily->isTrio()){
							 if  ($v->isMotherTransmission($family, $patient)){
							 	$hgenes->{$g->id}->{score_mother} = $this_score;
							 }
							 elsif  ($v->isFatherTransmission($family, $patient)){
							 	$hgenes->{$g->id}->{score_father} = $this_score;
							 }
							 else {
							 	$hgenes->{$g->id}->{score_biallelic} = $this_score;
							 }
						}
						else {
							$hgenes->{$g->id}->{score_biallelic} = $this_score;
						}
						my $this_alert = 0;
						if ($hgenes->{$g->id}->{validations} > 2) { $this_alert = 3; }
						if ($hgenes->{$g->id}->{validations} > 4) { $this_alert = 4; }
						$max_alert_project = $this_alert if ($max_alert_project < $this_alert);
						delete $v->{hgmd};
						delete $v->{isNewHgmd};
						delete $v->{hgmd_details};
						delete $v->{hgmd_id};
						delete $v->{hgmd_class};
						delete $v->{isDM};
						delete $v->{hgmd_inheritance};
						delete $v->{hgmd_hgvs};
						delete $v->{hgmd_phenotype};
						delete $v->{hgmd_disease};
						delete $v->{hgmd_releases};
						delete $v->{hgmd_pubmed};
						delete $v->{clinvar};
						delete $v->{text_clinvar};
						delete $v->{score_clinvar};
						delete $v->{isClinvarPathogenic};
						delete $v->{clinical};
						$v->{isDM} = 1 if (exists $hAllVarIds->{$var_id}->{hgmd});
						$v->{isClinvarPathogenic} = 1 if (exists $hAllVarIds->{$var_id}->{clinvar});
						update_variant_editor::vhgmd($v,$hvariation);
						update_variant_editor::vclinvar($v,$hvariation);
						if (exists $hAllVarIds->{$var_id}->{hgmd}) { $hvariation->{html}->{hgmd} .= "<img src='images/polyicons/new.jpeg'>"; }
						if (exists $hAllVarIds->{$var_id}->{clinvar}) { $hvariation->{html}->{clinvar} .= "<img src='images/polyicons/new.jpeg'>"; }
						if ($can_use_hgmd and $v->hgmd() and $v->hgmd_phenotype()) {
							$hvariation->{hgmd_phenotype} = $v->hgmd_phenotype();
						}
						$hvariation->{html}->{hgmd_no_access} = qq{<span class="glyphicon glyphicon-ban-circle" aria-hidden="true" style='font-size:12px;color:black;'></span>};
						update_variant_editor::table_validation($patient,$hvariation,$g);
						$hvariation->{html}->{table_validation} =~ s/Local/Local at $date_now/;
						$hvariation->{html}->{table_validation} =~ s/view_variation_validation\(/view_variation_validation_for_project\('$proj_name', /;
						if (exists $hAllVarIds->{$var_id}->{local}) {
							$hvariation->{html}->{table_validation} =~ s/pathogenic/pathogenic <img src='images\/polyicons\/new.jpeg'/;
						}
						
						$h_local_scaled_score_done->{$hvariation->{scaled_score}} = undef;
						unless (exists $h_families_done->{$family->name()}) {
							update_variant_editor::trio($v, $hvariation, $patient);
							my $validation = $patient->getLatestValidationStatus($user_name);
							my $term = "";
							my $display = "visibility:hidden;";
							my $date2 = "";
							my $patient_name = $patient->name;
							my $stamp;
							if ($validation){
								$term = $validation->{term};
								$date2 = join("-",return_date($validation->{modification_date}));
								$display = "";
								$stamp  = qq{ <span id="stamp_$patient_name" class="stamp is-approved_red" style="$display"><span>$term</span><br><small>$date2</small></span>};
								push(@lTrios, "<table><td><b><span style='font-size:10px;'>Fam ".$patient->getFamily->name().'</b></span></td><td>'.$stamp.'</td></table><br>'.$hvariation->{html}->{trio});
							}
							else {
								push(@lTrios, "<b><span style='font-size:10px;'>Fam ".$patient->getFamily->name().'</b><br>'.$hvariation->{html}->{trio});
							}
						}
						$hvariation->{html}->{trio} = join('<br>', @lTrios);
						if (exists $hgenes->{$g->id}->{pathogenic} and $hvariation->{genes_pathogenic}->{$g->id()}) {
							$hvariation->{html}->{$g->id()}->{'table_transcript'} = update_variant_editor::table_transcripts($hvariation->{genes}->{$g->id()}, \@header_transcripts);
						}
						$h_families_done->{$family->name()} = undef;
						delete $hvariation->{obj};
						$h_scaled_score->{$this_score}->{$g->id()}->{$var_id} = $hvariation;
						my $html_var_select = update_variant_editor::validation_select($patient, $v, $g);
						$html_var_select =~ s/validation_acmg\(/validation_acmg_for_project\('$proj_name',/;
						$hvariation->{html}->{table_validation_select}->{$g->external_name()}->{$patient->name()} = $html_var_select;
						$v = undef;
					}
				}
			}
			$v = undef;
		}
		foreach my $chr_id (keys %$hChrUsed) { purge($project->getChromosome($chr_id)); }
				
		my $nb_red = 0;
		my $nb_coral = 0;
		my $nb_orange = 0;
		my $nb_yellow = 0;
		my $nb_grey = 0;
		my $nb_red_no_hgmd = 0;
		my $nb_coral_no_hgmd = 0;
		my $nb_orange_no_hgmd = 0;
		my $nb_yellow_no_hgmd = 0;
		my $nb_grey_no_hgmd = 0;
		my ($varIds_colors, $varIds_colors_no_hgmd);
		foreach my $scale_score (keys %{$h_scaled_score}) {
			my $color;
			if ($scale_score >= $value_red) { $color = 'red'; }
			elsif ($scale_score >= $value_coral) { $color = 'coral'; }
			elsif ($scale_score >= $value_orange) { $color = 'orange'; }
			elsif ($scale_score >= $value_yellow) { $color = 'yellow'; }
			else { $color = 'grey'; }
			foreach my $g_id (keys(%{$h_scaled_score->{$scale_score}})) {
				foreach my $var_id (keys(%{$h_scaled_score->{$scale_score}->{$g_id}})) {
					$varIds_colors->{$var_id}->{$color}++;
					if (exists $hAllVarIds->{$var_id}->{clinvar}) {
						$varIds_colors_no_hgmd->{$var_id}->{$color.'_no_hgmd'}++;
					}
				}
			}
		}
		foreach my $var_id (keys(%{$varIds_colors})) {
			if (exists $varIds_colors->{$var_id}->{'red'}) { $nb_red++; }
			elsif (exists $varIds_colors->{$var_id}->{'coral'}) { $nb_coral++; }
			elsif (exists $varIds_colors->{$var_id}->{'orange'}) { $nb_orange++; }
			elsif (exists $varIds_colors->{$var_id}->{'yellow'}) { $nb_yellow++; }
			else { $nb_grey++; }
		}
		foreach my $var_id (keys(%{$varIds_colors_no_hgmd})) {
			if (exists $varIds_colors_no_hgmd->{$var_id}->{'red_no_hgmd'}) { $nb_red_no_hgmd++; }
			elsif (exists $varIds_colors_no_hgmd->{$var_id}->{'coral_no_hgmd'}) { $nb_coral_no_hgmd++; }
			elsif (exists $varIds_colors_no_hgmd->{$var_id}->{'orange_no_hgmd'}) { $nb_orange_no_hgmd++; }
			elsif (exists $varIds_colors_no_hgmd->{$var_id}->{'yellow_no_hgmd'}) { $nb_yellow_no_hgmd++; }
			elsif (exists $varIds_colors_no_hgmd->{$var_id}->{'grey_no_hgmd'}) { $nb_grey_no_hgmd++; }
		}
		
		my ($hRes, $hResResume);
		$hResResume->{max_score} = -999;
		if ($h_scaled_score) {
			my @lScores = sort {$a <=> $b} keys %{$h_scaled_score};
			$hResResume->{max_score} = $lScores[-1];
		}
		$hResResume->{max_alert} = $max_alert_project;
		$hResResume->{nb_red} = $nb_red;
		$hResResume->{nb_coral} = $nb_coral;
		$hResResume->{nb_orange} = $nb_orange;
		$hResResume->{nb_yellow} = $nb_yellow;
		$hResResume->{nb_grey} = $nb_grey;
		$hResResume->{nb_red_no_hgmd} = $nb_red_no_hgmd;
		$hResResume->{nb_coral_no_hgmd} = $nb_coral_no_hgmd;
		$hResResume->{nb_orange_no_hgmd} = $nb_orange_no_hgmd;
		$hResResume->{nb_yellow_no_hgmd} = $nb_yellow_no_hgmd;
		$hResResume->{nb_grey_no_hgmd} = $nb_grey_no_hgmd;
		$hResResume->{project_name} = $project->name();
		$hResResume->{project_description} = $project->description();
		$hResResume->{project_phenotypes} = join(', ', @{$project->phenotypes()});
		$hRes->{res_by_score} = $h_scaled_score;
		$hRes->{hgenes} = $hgenes;
		
		unless (-e $json_file) {
			my $json_res = encode_json $hRes;
			my $json_res_resume = encode_json $hResResume;
			
			unless (-d $dir_proj) {
				my $cmd1 = "mkdir $dir_proj";
				`$cmd1`;
				my $cmd2 = "chmod 777 $dir_proj";
				`$cmd2`;
			}
			open (FILE, ">$json_file");
			print FILE $json_res;
			close (FILE);
			
			open (FILE, ">$json_file_resume");
			print FILE $json_res_resume;
			close (FILE);
		}
		$pm->finish(0, $hResResume);
	}
	$pm->wait_all_children();
	warn "\n--- After FORK ---\n" if ($debug);
	
	warn "\n--- Return hFound hResProj---\n" if ($debug);
	return ($hFound, $hResProj);
}

sub get_hash_polyweb_users {
	my $dbh = $buffer->dbh();
	my $sql = qq{select u.user_id as id from bipd_users.USER u; };
	my $sth = $buffer->dbh->prepare($sql) || die();
	$sth->execute();
	my $res = $sth->fetchall_hashref("id");
	my $h_users;
	foreach my $uid (sort keys %$res){
		my $sql2 = qq{ SELECT u.nom_responsable as name  FROM bipd_users.`USER` u where u.user_id='$uid'; };
        my $sth2 = $buffer->dbh->prepare($sql2) || die();
        $sth2->execute();
        my $res2 = $sth2->fetchall_hashref('name');
        my @lUsers = keys %{$res2};
        $h_users->{lc($lUsers[0])} = undef;
	}
	return $h_users;
}

sub get_hash_new_hgmd_clinvar_strict {
	my ($hAllVarIds) = @_;
	my $hAllVarIds_new;
	my $h_polyweb_users = get_hash_polyweb_users();
	my $project = $buffer->newProject( -name => 'NGS2015_0794' );
	foreach my $var_id (keys %$hAllVarIds) {
		my $v = $project->_newVariant($var_id);
		next if ($v->isDM() and not exists $hAllVarIds->{$var_id}->{hgmd});
		next if ($v->isClinvarPathogenic() and not exists $hAllVarIds->{$var_id}->{clinvar});
		if (exists $hAllVarIds->{$var_id}->{'local'} and not exists $hAllVarIds->{$var_id}->{'hgmd'} and not exists $hAllVarIds->{$var_id}->{'clinvar'}) {
			next if ($v->isDM());
			next if ($v->isClinvarPathogenic());
		}
#		if ($v->isDM()) {
#			my $hgmd_id = $v->hgmd_id();
#			next unless ($hgmd_id);
#			my $hgmd_details = $buffer->queryHgmd()->getDataHGMDPro($hgmd_id);
#			my $author = lc($hgmd_details->{'author'}) if (exists $hgmd_details->{'author'});
#			next if ($author and exists $h_polyweb_users->{$author});
#		}
		$hAllVarIds_new->{$var_id} = $hAllVarIds->{$var_id};
	}
	return $hAllVarIds_new;
}

sub get_hash_new_hgmd_clinvar {
	my ($h_my_projects, $release) = @_;
	warn "\n--- Enter get_hash_new_hgmd_clinvar ---\n" if ($debug);
	my $dir = $dir_polybtf_path;
	my $json_global_file;
	if ($release) { $json_global_file = $dir.'/../'.$release.'/new_hgmd_clinvar.json'; }
	else { $json_global_file = $dir.'/new_hgmd_clinvar.json'; }
	warn $json_global_file if ($debug);
	
	my ($hAllVarIds, $nb_new_hgm, $nb_new_clinvar);
	
	if (-e $json_global_file) {
		open (FILE, "$json_global_file");
		my $json_encode = <FILE>;
		close (FILE);
		$hAllVarIds = decode_json $json_encode;
	}
	else {
		my $h_new = $buffer->queryHgmd->get_hash_last_released_DM();
		my $nb_acc_num = scalar (keys %{$h_new});
		
		foreach my $acc_num (keys %{$h_new}) {
			print '.';
			next unless $h_new->{$acc_num}->{tag} eq 'DM';
			my $chr_id = $h_new->{$acc_num}->{chrom};
			my $pos = $h_new->{$acc_num}->{pos};
			my $ref = $h_new->{$acc_num}->{ref};
			my $alt = $h_new->{$acc_num}->{alt};
			my $var_id = $chr_id.'_'.$pos.'_'.$ref.'_'.$alt;
			$hAllVarIds->{$var_id}->{hgmd} = 1;
			$nb_new_hgm++;
		}
		foreach my $var_id (@{$buffer->queryClinvarPathogenic->getAllVarIds_onlyLastRelease()}) {
			$hAllVarIds->{$var_id}->{clinvar} = 1;
			$nb_new_clinvar++;
		}
		foreach my $gene_var_id (keys %{$buffer->validations_query->getAllValidations()}) {
			my $last_h_validation = $buffer->validations_query->getAllValidations_last_months('5', '4')->{$gene_var_id}->[0];
			my $validation = $last_h_validation->{'validation'};
			my $validation_project_name = $last_h_validation->{'project_name'};
			next unless ($validation >= 5);
			my ($gene_id, $var_id) = split('!', $gene_var_id);
			$hAllVarIds->{$var_id}->{'local'}->{$validation_project_name} = 1;
			$nb_new_clinvar++;
		}
		
		my $hAllVarIds_new = get_hash_new_hgmd_clinvar_strict($hAllVarIds);
		$hAllVarIds = $hAllVarIds_new;
		
		my $json_encode = encode_json $hAllVarIds;
		open (FILE, ">$json_global_file");
		print FILE $json_encode;
		close (FILE);
	}
	return $hAllVarIds;
}

sub get_hash_new_hgmd_clinvar_with_dejavu_gnomad_filters {
	my ($hAllVarIds, $release) = @_;
	warn "\n--- Enter get_hash_new_hgmd_clinvar_with_dejavu_gnomad_filters ---\n" if ($debug);
	warn $hAllVarIds  if ($debug);
	warn $release  if ($debug);
	my $hFound;
	my $dir = $dir_polybtf_path;
	my $json_global_file;
	if ($release) { $json_global_file = $dir.'/../'.$release.'/found_projects.json'; }
	else { $json_global_file = $dir.'/found_projects.json'; }
	warn $json_global_file  if ($debug);
	
	if (-e $json_global_file) {
			open (FILE, $json_global_file);
			my $json_encode = <FILE>;
			close (FILE);
			$hFound = decode_json $json_encode;
	}
	else {
		die unless ($hAllVarIds);
		die if (scalar keys %$hAllVarIds == 0);
		my @lProj = sort keys %{$h_my_projects};
		my $buffer = new GBuffer;
		#my $project = $buffer->newProject( -name => $lProj[-1] );
		my $project = $buffer->newProject( -name => 'NGS2015_0794' );
		if ($release) {
			$project->changeAnnotationVersion($release, 1);
		}
		foreach my $var_id (keys %$hAllVarIds) {
			my @lTmp = split('_', $var_id);
			next unless (scalar @lTmp == 4);
			my $var = $project->_newVariant($var_id);
			
			my $type = 'snps';
			$type = 'insertions' if (length($lTmp[2]) < length($lTmp[3]));
			$type = 'deletions' if (length($lTmp[2]) > length($lTmp[3]));
			my $gnomad_ac = $buffer->get_gnomad($lTmp[0], $type, $lTmp[1], $lTmp[3])->{populations}->{all}->{AC};
			next if ($gnomad_ac > $max_gnomad);
			next if ($var->getGnomadHO() > $max_gnomad_ho);
			my $chr_id = $var->getChromosome->id();
			my $hDejaVu = $project->getDejaVuInfos($var->id());
			
			next unless ($hDejaVu);
			next if ($var->nb_deja_vu_samples() > $max_dejavu);
			my $have_one_project;
			foreach my $proj_name (keys %$hDejaVu) {
#				$hFoundByVar->{$var_id} = undef;
				my $patients = $hDejaVu->{$proj_name}->{patients};
				$have_one_project++;
				foreach my $pat_name (split(';', $patients)) {
					$hFound->{$proj_name}->{$var_id}->{patients}->{$pat_name} = undef;
					$hFound->{$proj_name}->{$var_id}->{chr_id} = $chr_id;
				}
			}
		}
		$buffer = undef;
		
		my $json_res = encode_json $hFound;
		unless (-d $dir) {
			my $cmd1 = "mkdir $dir";
			`$cmd1`;
			my $cmd2 = "chmod 777 $dir";
			`$cmd2`;
		}
		open (FILE, ">$json_global_file");
		print FILE $json_res;
		close (FILE);
	}
	return $hFound;
}

sub return_date {
	my ($dd) = @_;
	my @amonths = ('Jan', 'Feb', 'Mar', 'Apr','May',"Jun","Jul","Aug","Sep","Oct","Nov","Dec");
	my ($date,$time) = split(" ",$dd);
    my ($year,$month,$day) =  split("-",$date);
	return ($year,$amonths[$month-1],$day);
}

sub purge {
	my $chr = shift;
	$chr->project->{objects}->{proteins}    = {};
	$chr->project->{objects}->{genes}       = {};
	$chr->project->{objects}->{deletions}   = {};
	$chr->project->{objects}->{insertions}  = {};
	$chr->project->{objects}->{variations}  = {};
	$chr->project->{objects}->{exons}       = {};
	$chr->project->{objects}->{transcripts} = {};
	$chr->available_genes_ids(1);
	$chr->hash_filters_deleted(1);
	$chr->hash_freeze_file_genes(1);
	$chr->hash_freeze_file_all_genes(1);
	$chr->hash_freeze_file_patients(1);
	$chr->hash_filters_keeped(1);
	$chr->genes_object(1);
	delete $chr->project->{objects}->{chromosomes}->{$chr->id()}->{fastGenes};
	delete $chr->project->{objects}->{chromosomes}->{$chr->id()}->{fastTranscripts};
	$chr->cache_lmdb_variations->close();
	$chr->get_lmdb_patients->close();
	$chr->get_lmdb_categories->close();
	$chr->get_lmdb_genes->close() if ($chr->get_lmdb_genes());
	$chr->get_lmdb_variations->close();
}
