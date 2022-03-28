#!/usr/bin/perl
$|=1;
use CGI qw/:standard :html3/;
use JSON;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/";
use lib "$Bin/../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../GenBo/lib/obj-nodb/packages/cache/polydiag/";
use Data::Dumper;
use GBuffer;
use QueryVcf;
use Getopt::Long;
use update;
use QueryVectorFilter;




my $cgi = new CGI();
my $user		= $cgi->param('user');
my $pwd			= $cgi->param('pwd');
my $sort_by		= $cgi->param('sort_by');
my $sort_by_asc_desc = $cgi->param('sort_by_asc_desc');
my $only_project = $cgi->param('project');
my $only_patients = $cgi->param('patient');
my $only_genes = $cgi->param('genes');
my $interface_view = $cgi->param('interface_view');
my $only_new_hgmd = $cgi->param('only_new_hgmd');

print_cgi_header_html($cgi);

unless ($user) {
	die("\n\nuser option mandatory ! Die...\n\n");
}

unless ($pwd) {
	die("\n\npwd option mandatory ! Die...\n\n");
}

my $buffer = GBuffer->new();
my $query = $buffer->getQuery();


my ($h_only_patients, $h_only_genes, $h_genes_by_captures);
if ($only_patients) {
	foreach my $pat_name (split(',', $only_patients)) {
		$h_only_patients->{$pat_name} = undef;
	}
}

confess("\n\nERROR: DIE - only one patient in argument...\n\n") if (scalar keys %{$h_only_patients} > 1);

if ($only_genes =~ /capture/) {
	foreach my $label1 (keys %{$buffer->getAllGenesNamesInAllBundle()}) {
		foreach my $capture_name (keys %{$buffer->getAllGenesNamesInAllBundle->{$label1}}) {
			foreach my $gene_name (keys %{$buffer->getAllGenesNamesInAllBundle->{$label1}->{$capture_name}->{genes}}) {
				$h_genes_by_captures->{$capture_name}->{$gene_name} = undef;
			}
		}
	}
}
if ($only_genes) {
	foreach my $gene_name (split(',', $only_genes)) {
		if ($gene_name =~ /capture/) {
			my ($tmp, $capture_name) = split(':', $gene_name);
			foreach my $this_gene_name (keys %{$h_genes_by_captures->{$capture_name}}) {
				$h_only_genes->{$this_gene_name} = undef;
			}
		}
		else { $h_only_genes->{$gene_name} = undef; }
	}
}

unless ($sort_by) { $sort_by = 'pos'; }
unless ($sort_by_asc_desc) { $sort_by_asc_desc = 'asc'; }

unless ($buffer->queryHgmd->getHGMD($user) == 1) {
	my $out = qq{<tbody>};
	$out .= $cgi->start_div( { class => "panel panel-danger", style => "background-color:red;overflow:auto;border: 1px solid black;color-adjust: exact;padding:2px;" } );
	$out .= qq{<span style="color:white;">ERROR: This user $user is not authorized to view HGMD DB... </span>};
	$out .= qq{</div>};
	$out .= qq{</tbody>};
	print $out;
	exit(0);
}

unless ($interface_view) {
	$interface_view = 'by_captures_genes';
}

my $hBundles;
my $hBundlesTmp = $query->getAllGenesNamesInAllBundle();
foreach my $bundle_name (keys %{$hBundlesTmp}) {
	$hBundles->{$bundle_name} = undef;
	foreach my $capture_name (keys %{$hBundlesTmp->{$bundle_name}}) {
		$hBundles->{$capture_name} = undef;
	}
}

my $hProjectUser = $query->getProjectHashForUser($user, $pwd);
my @lProjectsForUser = keys %{$hProjectUser};

my ($hProjAuthorized, $hProjAuthorized_by_analyse, $last_project_name);
unless (@lProjectsForUser > 0) {
	my $out = qq{<tbody>};
	$out .= $cgi->start_div( { class => "panel panel-danger", style => "background-color:red;overflow:auto;border: 1px solid black;color-adjust: exact;padding:2px;" } );
	$out .= qq{<span style="color:white;">ERROR: No project found with this user ($user) and this password or error in user name / password... Reload the page ;)</span>};
	$out .= qq{</div>};
	$out .= qq{</tbody>};
	print $out;
	exit(0);
}

foreach my $project_name (@lProjectsForUser) {
	next unless ($project_name =~ /NGS20/);
	$hProjAuthorized->{$project_name} = undef;
	$last_project_name = $project_name;
}

if (not $hProjAuthorized or not exists $hProjAuthorized->{$only_project}) {
	my $out = qq{<tbody>};
	$out .= $cgi->start_div( { class => "panel panel-danger", style => "background-color:red;overflow:auto;border: 1px solid black;color-adjust: exact;padding:2px;" } );
	$out .= qq{<span style="color:white;">ERROR: No project $only_project found with this user $user and this password... Reload the page ;) </span>};
	$out .= qq{</div>};
	$out .= qq{</tbody>};
	print $out;
	exit(0);
}

foreach my $analyse (@{$query->listAllCaptureAnalyse()}) {
	foreach my $project_name (@{$query->listAllProjectsNameByAnalyse($analyse)}) {
		if (exists $hProjAuthorized->{$project_name}) {
			$hProjAuthorized->{$project_name} = $analyse;
			$hProjAuthorized_by_analyse->{$analyse}->{$project_name} = undef;
		}
	}
}

my $hVcf;
my $project;
if ($only_project) { $project = $buffer->newProjectCache( -name => $only_project); }
else {
	my $out = qq{<tbody>};
	$out .= $cgi->start_div( { class => "panel panel-danger", style => "background-color:red;overflow:auto;border: 1px solid black;color-adjust: exact;padding:2px;" } );
	$out .= qq{<span style="color:white;">ERROR: no project selected !</span>};
	$out .= qq{</div>};
	$out .= qq{</tbody>};
	print $out;
	exit(0);
}
$project->cgi_object(1);

my ($out_dialog);
my (@lVar, $hToSort, $hRes, $hGenesNewHgmd, $hGenePhen, $hPatAtticByVarPhen);
print "<div hidden>";
my $h_old_hgmd_acc_num = $project->buffer->queryHgmd->hashOldAccNum();
foreach my $chr (@{$project->getChromosomes()}) {
	#print "<div hidden>";
	my (@lVarHgmdDm, $h_var_project);
	my $chr_name = $chr->name();
	
	my $vector_chr = $chr->getVariantsVector->Clone();
	$vector_chr->Empty();
	if ($only_genes) {
		QueryVectorFilter::filter_genes_only_genes_names($chr, $only_genes);
		$vector_chr += $chr->getVariantsVector(); 
	}
	else {
		foreach my $gene (@{$chr->getGenes()}) {
			if ($gene->is_HGMD_DM()) {
				$vector_chr += $gene->getVariantsVector();
			}
		}
	}
	if ($vector_chr->is_empty()) {
		print "@";
		#print "</div>";
		next;
	}
	
	foreach my $var (@{$chr->getListVarObjects($vector_chr)}) {
		$project->print_dot(100);
		delete $var->{hgmd};
		delete $var->{isNewHgmd};
		delete $var->{hgmd_details};
		delete $var->{hgmd_id};
		delete $var->{hgmd_class};
		delete $var->{isDM};
		delete $var->{hgmd_inheritance};
		delete $var->{hgmd_hgvs};
		delete $var->{hgmd_phenotype};
		delete $var->{hgmd_disease};
		delete $var->{hgmd_releases};
		delete $var->{hgmd_pubmed};
		if ($var->isDM()) {
			if ($only_new_hgmd) {
				if (not $var->isNewHgmd()) { $vector_chr->Bit_Off($var->vector_id()); }
			}
		} 
		else { $vector_chr->Bit_Off($var->vector_id()); }
	}
	my $nb_var = $chr->countThisVariants($vector_chr);
	if ($vector_chr->is_empty()) {
		print "@";
		#print "</div>";
		#print qq{<button type="button" class="btn btn-secondary progress_chr">chr$chr_name ($nb_var)</button>};
		next;
	}
	
	foreach my $var (@{$chr->getListVarObjects($vector_chr)}) {
		#warn "\n\n";
		#warn $var->id().' -> '.ref($var);
		my $var_id = $var->id();
		my ($h_var_info, $hProjDejavu_user) = get_var_informations($var);
		
		if ($h_var_info and $hProjDejavu_user) {
			foreach my $project_name (keys %{$h_var_info->{polyquery_patients_attic}}) {
				$hPatAtticByVarPhen->{$var_id}->{$project_name} = $h_var_info->{polyquery_patients_attic}->{$project_name};
			}
			foreach my $analyse (keys %$hProjDejavu_user) {
				my $hgmd_phen = $h_var_info->{hgmd_phen};
				my $hgmd_ad_ar = $h_var_info->{hgmd_ad_ar};
				my $hgmd_new = $h_var_info->{hgmd_new};
				my @lGenes = split('\+', $h_var_info->{gene});
				foreach my $gene_name (@lGenes) {
					#$hGenePhen->{$gene_name}->{$hgmd_phen}++;
					$hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos} = $h_var_info;
					foreach my $project_name (keys %{$hProjDejavu_user->{$analyse}}) {
						foreach my $pat_name (keys %{$hProjDejavu_user->{$analyse}->{$project_name}}) {
							my $he_ho = $hProjDejavu_user->{$analyse}->{$project_name}->{$pat_name}->{he_ho};
							my $status = $hProjDejavu_user->{$analyse}->{$project_name}->{$pat_name}->{status};
							my $family = $hProjDejavu_user->{$analyse}->{$project_name}->{$pat_name}->{family};
							my $image = $hProjDejavu_user->{$analyse}->{$project_name}->{$pat_name}->{img_child_mother_father};
							$hGenesNewHgmd->{$gene_name}->{has_new_hgmd} = 1 if ($hgmd_new);
							$hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{projects}->{$project_name}->{$pat_name}->{he_ho} = $he_ho;
							$hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{projects}->{$project_name}->{$pat_name}->{status} = $status;
							$hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{projects}->{$project_name}->{$pat_name}->{family} = $family;
							$hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{projects}->{$project_name}->{$pat_name}->{img_child_mother_father} = $image;
							$hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{var_infos} = $h_var_info;
						}
					}
				}
			}
		}
	}
	print "@";
	#print "</div>";
	#print qq{<button type="button" class="btn btn-success progress_chr">chr$chr_name ($nb_var)</button>};
}
print "</div>";
#exit(0);

if ($interface_view eq 'by_captures_genes') { draw_interface_by_captures_genes($hRes, $project->name()); }
elsif ($interface_view eq 'by_variants') { draw_interface_by_variants($hRes, $project->name()); }
else { draw_interface_by_captures_genes($hRes, $project->name()); }


sub draw_interface_by_variants {
	my ($hRes, $project_name) = @_;
	
	my $hVar;
	my $hPhen_by_varid;
	
	my $hVarInfos;
	foreach my $analyse (sort keys %{$hRes}) {
		foreach my $gene_name (keys %{$hRes->{$analyse}}) {
			my $modif_gene_name = $gene_name;
			$modif_gene_name =~ s/\./_/g;
			foreach my $phen_name (keys %{$hRes->{$analyse}->{$gene_name}}) {
				foreach my $var_id (keys %{$hRes->{$analyse}->{$gene_name}->{$phen_name}}) {
					$hVar->{$var_id} = undef;
					$hPhen_by_varid->{$var_id}->{$gene_name} = $phen_name;
					$hVarInfos->{$var_id} = $hRes->{$analyse}->{$gene_name}->{$phen_name}->{$var_id}->{var_infos};
				}
			}
		}
	}
	
	my $buffer = GBuffer->new();
	my $project = $buffer->newProjectCache( -name => $only_project);
	my $vquery;
	if ($project->validation_db()) {
		$vquery = validationQuery->new(
			dbh          => $buffer->dbh,
			capture_name => $project->validation_db()
		);
	}

	my $out = qq{<table data-toggle="table" class="display table table-bordered table-hover table-responsive" style="width:100%;height:600px;">};
	$out .= qq{<thead>};
	$out .= qq{<tr>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>HGMD_Phen</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>HGMD_Class</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Gene</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Consequence</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Patient(s)</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Rs_Name</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Var_Name</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>HGMD_ID</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>OMIM_ID</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Freq</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Freq_HO</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Max_POP</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Min_POP</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Clinvar</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>DejaVu</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Ngs</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Ratio</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Caller</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Genomique</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Transcript</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Exon</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Nomenclature</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Consequence</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Codons</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Codons_AA</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Polyphen</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Sift</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Cadd</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>IGV_View</b></th>};
	$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>Alamut</b></th>};
	$out .= qq{</tr>};
	$out .= qq{</thead>};
	$out .= qq{<tbody style="width:100%">};
	
	foreach my $var_id (sort keys %{$hVar}) {
		my $var = $project->_newVariant($var_id);
		next unless ($var->isDM());
		my $patient = $project->getPatient($only_patients);
		my @lTrans = @{$var->getTranscripts()};
		my $nb_trans = scalar(@lTrans);
		my $i = 0;
		foreach my $tr (@lTrans) {
			$i++;
			
			my $h_res = update::construct_variant($project, $var, $tr, $patient, $vquery);
			my $var_type = $project->buffer->config->{'stats_genes'}->{$h_res->{consequence}};
			my $color_type_transcipt = 'green';
			if ($var_type eq 'medium') { $color_type_transcipt = 'orange'; }
			if ($var_type eq 'high') { $color_type_transcipt = 'red'; }
			
			my $cons_gene = $var->variationType($tr->getGene());
			my $cons_gene_interface = $var->variationTypeInterface($tr->getGene());
			my $var_type_gene = $project->buffer->config->{'stats_genes'}->{$cons_gene};
			my $color_type_gene = 'green';
			if ($var_type_gene eq 'medium') { $color_type_gene = 'orange'; }
			if ($var_type_gene eq 'high') { $color_type_gene = 'red'; }
			
			if ($i == 1) {
				my $h;
				$h->{obj} = $var;
				my ($class, $cmd_hgmd) = update::hgmd($project,$h,'return_cmd');
				my $class_name = $hVarInfos->{$var_id}->{hgmd_class};
				my $hgmd_class = qq{<button onClick="$cmd_hgmd">$class_name</button>};

				my $clinvar_text;
				$clinvar_text = $var->text_clinvar() unless ($var->text_clinvar() eq '-5');
				
				my $dejavu = $hVarInfos->{$var_id}->{dejavu};
				$dejavu =~ s/\//<br>/;
				my $dejavu_url = "http://www.polyweb.fr/polyweb/polydejavu/dejavu.html?input=".$var_id;
				my $dejavu_href = "<a href='".$dejavu_url."' target='_blank'>".$dejavu."</a>";
				
				my $omim_href = "<a href='https://www.omim.org/entry/".$hVarInfos->{$var_id}->{omimid}."' target='_blank'>".$hVarInfos->{$var_id}->{omimid}."</a>";
				
				my $rsname_href;
				if ($var->rs_name()) {
					$rsname_href = "<a href='https://www.ncbi.nlm.nih.gov/snp/?term=".$var->rs_name()."' target='_blank'>".$var->rs_name()."</a>";
				}
				
				my $polyphen = "<img  width=16 height=16 src='/icons/Polyicons/difficulty_3_b.png'>";
				if ( $h_res->{polyphen} == 3 ) { $polyphen = "<img  width=16 height=16 src='/icons/Polyicons/difficulty_2_1.png'>"; }
				if ( $h_res->{polyphen} == 2 ) { $polyphen = "<img  width=16 height=16 src='/icons/Polyicons/difficulty_1.png'>"; }
				if ( $h_res->{polyphen} == 1 ) { $polyphen = "<img  width=16 height=16 src='/icons/Polyicons/difficulty_0.png'>"; }
				
				my $sift = "<img  width=16 height=16 src='/icons/Polyicons/difficulty_3_b.png'>";
				if ( $h_res->{sift} == 2 ) { $sift = "<img  width=16 height=16 src='/icons/Polyicons/difficulty_2_1.png'>"; }
				if ( $h_res->{sift} == 1 ) { $sift = "<img  width=16 height=16 src='/icons/Polyicons/difficulty_0.png'>"; }
				
				$out .= '<tr>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center; color:red;'><b>".$hVarInfos->{$var_id}->{hgmd_phen}.' ['.$hVarInfos->{$var_id}->{hgmd_ad_ar}.']</b></td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$hgmd_class.'</td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center; color:$color_type_gene;'><b>".$h_res->{gene}.'</b></td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center; color:$color_type_gene;'>".$cons_gene_interface.'</td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>";
				my $p = 0;
				foreach my $this_patient (@{$var->getPatients()}) {
					$p++;
					if ($p > 1) { $out .= "<br>"; }
					if (exists $h_only_patients->{$this_patient->name()}) {
						$out .= "<span style='color:red;'>";
						$out .= $this_patient->name();
						if    ($this_patient->status eq '1') { $out .= qq{ <img src='/icons/Polyicons/bullet_green.png'>}; }
						elsif ($this_patient->status eq '2') { $out .= qq{ <img src='/icons/Polyicons/pill2.png'>}; }
						$out .= "</span>";
					}
					else {
						$out .= $this_patient->name();
						if    ($this_patient->status eq '1') { $out .= qq{ <img src='/icons/Polyicons/bullet_green.png'>}; }
						elsif ($this_patient->status eq '2') { $out .= qq{ <img src='/icons/Polyicons/pill2.png'>}; }
					}
				}
				$out .= '</b></td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$rsname_href.'</td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$h_res->{id}.'</td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$hVarInfos->{$var_id}->{hgmd_id}.'</td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$omim_href.'</td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$var->frequency().'</td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$var->frequency_homozygote().'</td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$var->max_pop_name().':'.$var->max_pop_freq().'</td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$var->min_pop_name().':'.$var->min_pop_freq().'</td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$clinvar_text.'</td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$dejavu_href.'</td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$h_res->{ngs}.'</td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$h_res->{ratio}.'</td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$h_res->{caller}.'</td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$h_res->{genomique}.'</td>';
				$out .= "<td rowspan='1'><span style='color:$color_type_transcipt;'>".$h_res->{transcript}.'</span></td>';
				$out .= "<td rowspan='1'><span style='color:$color_type_transcipt;'>".$h_res->{exon}.'</span></td>';
				$out .= "<td rowspan='1'><span style='color:$color_type_transcipt;'>".$h_res->{nomenclature}.'</span></td>';
				$out .= "<td rowspan='1'><span style='color:$color_type_transcipt;'>".$h_res->{consequence_interface}.'</span></td>';
				$out .= "<td rowspan='1'><span style='color:$color_type_transcipt;'>".$h_res->{codons}.'</span></td>';
				$out .= "<td rowspan='1'><span style='color:$color_type_transcipt;'>".$h_res->{codons_AA}.'</span></td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$polyphen.'</td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$sift.'</td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$var->cadd_score().'</td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$h_res->{igv}.'</td>';
				$out .= "<td rowspan='$nb_trans' style='vertical-align : middle;text-align:center;'>".$h_res->{alamut}.'</td>';
				$out .= '</tr>';
			}
			else {
				$out .= '<tr>';
				#$out .= "<td colspan='11'></td>";
				$out .= "<td rowspan='1'><span style='color:$color_type_transcipt;'>".$h_res->{transcript}.'</span></td>';
				$out .= "<td rowspan='1'><span style='color:$color_type_transcipt;'>".$h_res->{exon}.'</span></td>';
				$out .= "<td rowspan='1'><span style='color:$color_type_transcipt;'>".$h_res->{nomenclature}.'</span></td>';
				$out .= "<td rowspan='1'><span style='color:$color_type_transcipt;'>".$h_res->{consequence_interface}.'</span></td>';
				$out .= "<td rowspan='1'><span style='color:$color_type_transcipt;'>".$h_res->{codons}.'</span></td>';
				$out .= "<td rowspan='1'><span style='color:$color_type_transcipt;'>".$h_res->{codons_AA}.'</span></td>';
				#$out .= "<td colspan='8'></td>";
				$out .= '</tr>';
			}
		}
	}
	$out .= qq{</table>};
	print $out;
	
}




sub draw_interface_by_captures_genes {
	my ($hRes, $project_name) = @_;
	unless ($hRes) {
		my $out = qq{<tbody>};
		$out .= $cgi->start_div( { class => "panel panel-warning", style => "background-color:orange;overflow:auto;border: 1px solid black;color-adjust: exact;padding:2px;" } );
		$out .= qq{<span style="color:white;">No results found...</span>};
		$out .= qq{</div>};
		$out .= qq{</tbody>};
		print $out;
		exit(0);
	}
	foreach my $analyse (sort keys %{$hRes}) {
		my $out;
		my $h_proj_pat_seen;
		my ($hPhen, $hProj, $hVarIds);
		my ($nb_gene, $nb_phen, $nb_var, $nb_proj);
		foreach my $gene_name (keys %{$hRes->{$analyse}}) {
			$nb_gene++;
			foreach my $phen_name (keys %{$hRes->{$analyse}->{$gene_name}}) {
				my $modif_gene_name = $gene_name;
				$modif_gene_name =~ s/\./_/g;
				push(@{$hPhen->{$phen_name}}, $modif_gene_name);
				push(@{$hPhen->{'all_genes'}}, $modif_gene_name);
				foreach my $var_id (keys %{$hRes->{$analyse}->{$gene_name}->{$phen_name}}) {
					$hVarIds->{$var_id} = undef;
					foreach my $proj_name (keys %{$hRes->{$analyse}->{$gene_name}->{$phen_name}->{$var_id}->{projects}}) {
						push(@{$hProj->{$proj_name}}, $modif_gene_name);
					}
				}
			}
		}
		$nb_proj = scalar keys %{$hProj};
		$nb_phen = scalar keys %{$hPhen};
		my $text_other_proj = join("\n", sort keys %{$hProj});
		my $text_other_phen = join("\n", sort keys %{$hPhen});
		
		my $collapse_analyse_proj_id = 'collapse_table_proj_'.$analyse.'_'.$project_name;
		my $collapse_analyse_phen_id = 'collapse_table_phen_'.$analyse.'_'.$project_name;
		
		
		$nb_var = scalar keys %{$hVarIds};
		$hVarIds = undef;
		my $collapse_id = 'table_'.$analyse.'_'.$project_name;
		$out .= $cgi->start_div( { class => "panel panel-sucess", style => "background-color:#99D3DF;overflow:auto;border: 1px solid black;color-adjust: exact;padding:2px;" } );
		$out .=  $cgi->start_div( { class => "btn-group btn-xs", style => "float:left;min-width:180px;" } );
			$out .=qq{<div class="btn  btn-info btn-xs " style="background-color:#88BBD6;position:relative;bottom:1px;min-width:300px;" data-toggle="collapse" data-target="#$collapse_id"> $analyse &nbsp </div>};
		$out .= $cgi->end_div();
		$out .=  $cgi->start_div( { class => "btn-group btn-xs", style => "float:right;min-width:180px;" } );
			$out .= qq{<div class="btn  btn-info btn-xs " style="background-color:#88BBD6;position:relative;bottom:1px;" data-toggle="collapse" data-target="#$collapse_id"> $nb_gene gene(s)&nbsp </div>};
			
			$out .= qq{<div class="btn  btn-info btn-xs " style="background-color:#88BBD6;position:relative;bottom:1px;" data-toggle="collapse" data-target="#$collapse_analyse_phen_id">};
			$out .= qq{<a href="#" data-toggle="tooltip" style="padding-left:5px;" data-placement="left" title="$text_other_phen"><span style="color:white;"><i><u> $nb_phen phenotype(s)&nbsp </i></u></span></a>};
			$out .= qq{</div>};
			
			$out .= qq{<div class="btn  btn-info btn-xs " style="background-color:#88BBD6;position:relative;bottom:1px;" data-toggle="collapse" data-target="#$collapse_analyse_proj_id">};
			$out .= qq{<a href="#" data-toggle="tooltip" style="padding-left:5px;" data-placement="left" title="$text_other_proj"><span style="color:white;"><i><u> $nb_proj project(s)&nbsp </i></u></span></a>};
			$out .= qq{</div>};
			
			$out .= qq{<div class="btn  btn-info btn-xs " style="background-color:#88BBD6;position:relative;bottom:1px;" data-toggle="collapse" data-target="#$collapse_id"> $nb_var variant(s)&nbsp </div>};
		$out .= $cgi->end_div();
		$out .= $cgi->end_div();
		
		
		$out .= qq{<div id="$collapse_analyse_proj_id" class="collapse" style="width:100%">};
		$out .= qq{<table data-toggle="table" class="display table table-bordered table-hover table-responsive" style="width:100%">};
		foreach my $proj_name (sort keys %{$hProj}) {
			my $nb_genes = scalar @{$hProj->{$proj_name}};
			
			my $description = $hProjectUser->{$proj_name}->{description};
			#my $all_genes = join(',', @{$hPhen->{'all_genes'}});
			#my $cmd_click_analyse = qq{collapse_close();};
			my $cmd_click_analyse .= qq{collapse_them('$collapse_analyse_proj_id', '$collapse_id', '};
			$cmd_click_analyse .= join(',', @{$hProj->{$proj_name}});
			$cmd_click_analyse .= qq{');};
			
			$out .= qq{<tr>};
			$out .= qq{<td><div class="btn  btn-info btn-xs " style="background-color:#88BBD6;position:relative;bottom:1px;" onClick="$cmd_click_analyse"> View </div></td>};
			$out .= qq{<td>$nb_genes gene(s)</td>};
			$out .= qq{<td>$proj_name</td>};
			$out .= qq{<td>$description</td>};
			$out .= qq{</tr>};
		}
		$out .= qq{</table>};
		$out .= qq{</div>};
		
		$out .= qq{<div id="$collapse_analyse_phen_id" class="collapse" style="width:100%">};
		$out .= qq{<table data-toggle="table" class="display table table-bordered table-hover table-responsive" style="width:100%">};
		foreach my $phen_name (sort keys %{$hPhen}) {
			my $nb_genes = scalar @{$hPhen->{$phen_name}};
			
			#my $all_genes = join(',', @{$hPhen->{'all_genes'}});
			#my $cmd_click_analyse = qq{collapse_close();};
			my $cmd_click_analyse .= qq{collapse_them('$collapse_analyse_phen_id', '$collapse_id', '};
			$cmd_click_analyse .= join(',', @{$hPhen->{$phen_name}});
			$cmd_click_analyse .= qq{');};
			
			$out .= qq{<tr>};
			$out .= qq{<td><div class="btn  btn-info btn-xs " style="background-color:#88BBD6;position:relative;bottom:1px;" onClick="$cmd_click_analyse"> View </div></td>};
			$out .= qq{<td>$nb_genes gene(s)</td>};
			$out .= qq{<td>$phen_name</td>};
			$out .= qq{</tr>};
		}
		$out .= qq{</table>};
		$out .= qq{</div>};
		$text_other_phen = undef;
		$hPhen = undef;
		
		$out .= qq{<div id="$collapse_id" class="collapse" style="width:100%">};
		foreach my $gene_name (sort keys %{$hRes->{$analyse}}) {
			$nb_phen = scalar keys %{$hRes->{$analyse}->{$gene_name}};
			$nb_var = 0;
			foreach my $phen_name (keys %{$hRes->{$analyse}->{$gene_name}}) {
				$nb_var += scalar keys %{$hRes->{$analyse}->{$gene_name}->{$phen_name}};
			}
			my $collapse_id_2 = 'table_'.$analyse.'_'.$gene_name.'_'.$project_name;
			$collapse_id_2 =~ s/\./_/g;
			$out .= $cgi->start_div( { class => "panel panel-sucess", style => "background-color:#C7D8C6;overflow:auto;border: 1px solid black;color-adjust: exact;padding:2px;" } );
			$out .=  $cgi->start_div( { class => "btn-group btn-xs", style => "float:left;padding-left:40px;min-width:180px;" } );
				$out .=qq{<div class="btn  btn-info btn-xs " style="background-color:#7A9D96;position:relative;bottom:1px;min-width:200px;" data-toggle="collapse" data-target="#$collapse_id_2"> Gene $gene_name&nbsp </div>};
				#$out .=  qq{<span class='badge badge-xs ' style="font-size:10px;padding:5px;">$hbilanGene{$gname}->{"totalVar"} </span> };
			$out .= $cgi->end_div();
			
			if (exists $hGenesNewHgmd->{$gene_name}) {
				$out .=  $cgi->start_div( { class => "btn-group btn-xs", style => "float:left;padding-left:5px;min-width:180px;" } );
					$out .=qq{<div class="btn  btn-info btn-xs text-danger" style="background-color:#FF7373;position:relative;bottom:1px;" data-toggle="collapse" data-target="#$collapse_id_2"><span style="padding-left:3px;padding-right:3px;"><i>HGMD New!</i></span></div>};
				$out .= $cgi->end_div();	
			}
			
			$out .=  $cgi->start_div( { class => "btn-group btn-xs", style => "float:right;min-width:180px;padding-right:40px;" } );
			my @lOther_phen;
			my $nb_phen_write = 0;
			foreach my $phen_name (sort keys %{$hRes->{$analyse}->{$gene_name}}) {
				if ($nb_phen_write < 2) {
					$out .= qq{<div class="btn  btn-info btn-xs " style="background-color:#7A9D96;position:relative;bottom:1px;" data-toggle="collapse" data-target="#$collapse_id_2">};
					$out .= qq{<a href="#" data-toggle="tooltip" style="padding-left:5px;max-width:550px;" data-placement="left" title="$phen_name"><span style="color:white;"><i><u> $phen_name&nbsp </i></u></span></a>};
					$out .= qq{</div>};
				}
				else { push(@lOther_phen, $phen_name); }
				$nb_phen_write++;
			}
			if (@lOther_phen) {
				my $phen_name = '+ '.scalar(@lOther_phen).' phen';
				$text_other_phen = join("\n", @lOther_phen);
				$out .= qq{<div class="btn  btn-info btn-xs " style="background-color:#7A9D96;position:relative;bottom:1px;" data-toggle="collapse" data-target="#$collapse_id_2">};
				$out .= qq{<a href="#" data-toggle="tooltip" style="padding-left:5px;max-width:550px;" data-placement="left" title="$text_other_phen"><span style="color:white;"><i><u> $phen_name&nbsp </i></u></span></a>};
				$out .= qq{</div>};
			}
			$out .=qq{<div class="btn  btn-info btn-xs " style="background-color:#7A9D96;position:relative;bottom:1px;" data-toggle="collapse" data-target="#$collapse_id_2"> $nb_var variant(s)&nbsp </div>};
			$out .= $cgi->end_div();
			
			
			$out .= $cgi->end_div();
			
			$out .= qq{<div id="$collapse_id_2" class="collapse" style="width:100%">};
			foreach my $hgmd_phen (sort keys %{$hRes->{$analyse}->{$gene_name}}) {
				
				$out .= qq{<table data-toggle="table" class="display table table-bordered table-hover table-responsive" id="table_hgmd" style="width:100%">};
				$out .= qq{<thead>};
				$out .= qq{<tr>};
				$out .= qq{<th align='center' style="color:#6290D8;border-top: 1px solid black;border-left: 1px solid black;border-bottom: 1px solid black;"><b>$analyse</b></th>};
				$out .= qq{<th align='center' style="color:#449D44;border-top: 1px solid black;border-bottom: 1px solid black;"><b>Gene $gene_name</b></th>};
				$out .= qq{<th align='center' style="color:red;border-top: 1px solid black;border-right: 1px solid black;border-bottom: 1px solid black;"><b>$hgmd_phen</b></th>};
				$out .= qq{</tr>};
				$out .= qq{</thead>};
				$out .= qq{<tbody style="width:100%">};
				
				my $bg_color;
				my $nb_var = 0;
				foreach my $var_id (sort keys %{$hRes->{$analyse}->{$gene_name}->{$hgmd_phen}}) {
					$nb_var++;
					if ($nb_var % 2) { $bg_color = 'background-color:#F5F5F5'; }
					else { $bg_color = ''; }
					
					my $hgmd_ad_ar = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{hgmd_ad_ar};
					#my $hgmd_phen = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{hgmd_phen};
					my $type = uc(substr($hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{type}, 0, 3));
					my $ids = $var_id;
					
					if ($hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{rs_name}) {
						my $rsname = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{rs_name};
						my $rsname_href = "<a href='https://www.ncbi.nlm.nih.gov/snp/?term=".$rsname."' target='_blank'>".$rsname."</a>";
						$ids .= ' ; '.$rsname_href; 
					}
					
					$ids   .= ' ; HGMD: '.$hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{hgmd_id};
					if ($hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{omimid}) {
						my $omim_id = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{omimid};
						my $omim_href = "<a href='https://www.omim.org/entry/".$omim_id."' target='_blank'>".$omim_id."</a>";
						$ids   .= ' ; OMIM: '.$omim_href;
					}
					
					my $annot = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{annot};
					my $dna = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{hgmd_hgvs};
					my $prot = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{hgmd_prot};
					my $dejavu = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{dejavu};
					my $hgmd_class = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{hgmd_class};
					my $hgmd_versions = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{hgmd_versions};
					my $hgmd_disease = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{hgmd_disease};
					my $hgmd_new = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{hgmd_new};
					
					
					my $dejavu_url = "http://www.polyweb.fr/polyweb/polydejavu/dejavu.html?input=".$var_id;
					my $dejavu_href = "<a href='".$dejavu_url."' target='_blank'>".$dejavu."</a>";
			
			
			
					my $row_span_var = get_row_span_variant_td($hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{projects});
					my $nb_max_projects = scalar keys %{$hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{projects}};
					my $nb_project = 0;
					my $is_first_line = 1;
					foreach my $project_name (sort keys %{$hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{projects}}) {
						$out .= qq{<tr>};
						$nb_project++;
						my $text_project_name;
						if (lc($analyse) eq 'exome' or lc($analyse) eq 'genome' or lc($analyse) =~ /cilliome/ or lc($analyse) =~ /ciliome/ or lc($analyse) =~ /agilent/) {
							my ($chr_id, $pos, $ref, $all) = split('_', $var_id);
							#my $this_genes = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{gene};
							my $url_filters = "?project=".$project_name."&chromosome=".$chr_id."&gene=".$gene_name."&patients=".join(',', keys %{$hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{projects}->{$project_name}});
							my $interface_url = "http://www.polyweb.fr/polyweb/vector/detailProject.html".$url_filters;
							
							#my $patients_attic = $hPatAtticByVarPhen->{$var_id}->{$project_name};
							#my $url_filters = "?project=".$project_name."&type=ngs&use_args=filter_type_variation=;only_genes=".$this_genes.";mode=ind;level_ind=variation;level_fam=variation;filter_attic=".$patients_attic."+"; 
							#my $interface_url = "http://www.polyweb.fr/polyweb/vector/gene.html".$url_filters;
							$text_project_name = "<a href='".$interface_url."' target='_blank'>".$project_name."</a>";
						}
						else {
							my $interface_url = "http://www.polyweb.fr/polyweb/coverage.html?project=".$project_name.'&type=ngs';
							$text_project_name = "<a href='".$interface_url."' target='_blank'>".$project_name."</a>";
						}
						my $row_span_proj = get_row_span_project_td($hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{projects}->{$project_name});
						#$out .= qq{<td rowspan="$row_span_proj" style="min-width:150px;">};
						
						my $border_proj = "border-left: 1px solid black;";
						if ($nb_project == 1) { $border_proj .= "border-top: 1px solid black;"; }
						if ($nb_project == $nb_max_projects) { $border_proj .= "border-bottom: 1px solid black;"; }
						$out .= qq{<td rowspan="$row_span_proj" style="width:15%;$bg_color;$border_proj">};
						$out .= qq{$text_project_name};
						my $description = $hProjectUser->{$project_name}->{description};
						$out .= qq{<br>$description};
						$out .= qq{</td>};
						
						my $max_nb_patients = scalar keys %{$hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{projects}->{$project_name}};
						my $nb_patient = 0;
						foreach my $patient_name (sort keys %{$hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{projects}->{$project_name}}) {
							$nb_patient++;
							if ($nb_patient > 1) {
								$out .= qq{<tr>};
							}
							
							my $image = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{projects}->{$project_name}->{$patient_name}->{img_child_mother_father};
							
							my $family = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{projects}->{$project_name}->{$patient_name}->{family};
							if ($family and $family ne $patient_name) { $family = "[$family] "; }
							
							my $status = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{projects}->{$project_name}->{$patient_name}->{status};
							if    ($status eq '1') { $status = qq{<img src='/icons/Polyicons/bullet_green.png'>}; }
							elsif ($status eq '2') { $status = qq{<img src='/icons/Polyicons/pill2.png'>}; }
							else { $status = ''; }
							
							my $he_ho = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{projects}->{$project_name}->{$patient_name}->{he_ho};
							if ($he_ho eq 'ho') { $he_ho = qq{<span style="color:red">HO</span>}; }
							if ($he_ho eq 'he') { $he_ho = qq{<span style="color:green">HE</span>}; }
							
							my $border_pat = "";
							if ($is_first_line == 1 and $nb_patient == 1) { $border_pat .= "border-top: 1px solid black;"; }
							if ($nb_project == $nb_max_projects and $nb_patient == $max_nb_patients) { $border_pat .= "border-bottom: 1px solid black;"; }
							
							#my $b_igv = qq{<button type="button" class="btn btn-primary btn-xs" style="font-size:9px;"><i>IGV View</i></button>};
							my $b_igv;
							
							$out .= qq{<td style="width:20%;$bg_color;$border_pat">};
							$out .= qq{$b_igv [$he_ho] $image $status $family - $patient_name  };
							$out .= qq{</td>};
							
							if ($is_first_line == 1) {
								$out .= qq{<td rowspan="$row_span_var" style="width:65%;border-top: 1px solid black;border-bottom: 1px solid black;border-right: 1px solid black;$bg_color;">};
								#$out .= qq{<span style="color:red">$hgmd_phen</span>};
								my $nb_annot;
								my @lAnnot = split(';', $annot);
								foreach my $this_annot (@lAnnot) {
									$nb_annot++;
									#my ($gene_name, $gene_annot) = split(':', $this_annot);
									$out .= qq{<br>} if ($nb_annot > 1);
									$out .= qq{<b><u>Annot:</b></u> <span style="color:blue"> $this_annot  </span>};
	#								my @lText_other_phen;
	#								foreach my $phen (sort keys %{$hGenePhen->{$gene_name}}) {
	#									push(@lText_other_phen, '['.$hGenePhen->{$gene_name}->{$phen}.' var] '.$phen);
	#								}
	#								my $nb_others_hgmd_phen = scalar(@lText_other_phen);
	#								if ($nb_others_hgmd_phen > 2) {
	#									my $text_other_phen = join("\n", @lText_other_phen);
	#									$out .= qq{<a href="#" data-toggle="tooltip" style="padding-left:5px;" data-placement="top" title="$text_other_phen"><span style="color:orange;"><i>[$nb_others_hgmd_phen HGMD Phen Found]</i></span></a>};
	#								}
								}
								$out .= qq{<br><b><u>$type:</b></u> $ids};
								$out .= qq{<br><b><u>HGVS:</b></u> $dna};
								$out .= qq{<br><b><u>Prot:</b></u> $prot} if ($prot);
								$out .= qq{<br><b><u>DejaVu:</b></u> $dejavu_href};
								$out .= qq{<br><b><u>HGMD Class:</b></u> $hgmd_class};
								$out .= qq{<br><b><u>HGMD Transmission Mode:</b></u> $hgmd_ad_ar};
								$out .= qq{<br><b><u>HGMD Disease:</b></u> $hgmd_disease};
								$out .= qq{<br>};
								if ($hgmd_new) {
									$out .=qq{<b><div class="btn  btn-info btn-xs text-danger" style="background-color:#FF7373;position:relative;bottom:1px;padding-left:5px;" ><div style="padding-left:5px;padding-right:5px;>  </div><span style="padding-left:3px;padding-right:3px;"><i>New!</i></span></div><u>HGMD Versions:</b></u> $hgmd_versions</div>};
								}
								else {
									$out .= qq{<b><u>HGMD Versions:</b></u> $hgmd_versions};
								}
								my @lPubmedIds = sort keys %{$hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{hgmd_pubmed}};
								@lPubmedIds = reverse(@lPubmedIds);
								foreach my $pubmed_id (@lPubmedIds) {
									my $author = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{hgmd_pubmed}->{$pubmed_id}->{author};
									my $year = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{hgmd_pubmed}->{$pubmed_id}->{year};
									my $pubmed_title = $hRes->{$analyse}->{$gene_name}->{$hgmd_phen}->{$var_id}->{infos}->{hgmd_pubmed}->{$pubmed_id}->{title}.' ['.$year.', '.$author.' and al.]';
									my $pubmed_href = "<a href='https://www.ncbi.nlm.nih.gov/pubmed/?term=".$pubmed_id."' target='_blank'>".'<i><u>'.$pubmed_title."</i></u></a>";
									$out .= qq{<br><b><u>Pubmed:</b></u> $pubmed_href};
								}
								$out .= qq{</td>};
								$is_first_line = 0;
							}
							$out .= qq{</tr>};
							print $out;
							$out = '';
						}
					}
					
				}
				print qq{</table>};
			}
			print qq{</div>};
		}
		print qq{</div>};
	}
	print qq{</tbody>};
}

sub get_row_span_project_td {
	my $this_h = shift;
	return scalar keys %{$this_h};
}

sub get_row_span_variant_td {
	my $this_h = shift;
	my $nb = 0;
	foreach my $project_name (keys %{$this_h}) {
		foreach my $patient_name (keys %{$this_h->{$project_name}}) {
			$nb++;
		}
	}
	return $nb;
}

sub get_var_informations {
	my ($var) = @_;
	my $var_id = $var->id();
	my $is_ok;
	my ($h_var_info, $hProjDejavu_user);
	my $nb_project = 0;
	my $nb_patients = 0;
	my $h_dejavu = $project->getDejaVuInfos($var_id);
	foreach my $project_name (keys %{$h_dejavu}) {
		next unless ($project_name =~ /NGS20/);
		next if ($only_project and $project_name ne $only_project);
		$nb_project++;
		$nb_patients += $h_dejavu->{$project_name}->{nb};
		next unless exists ($hProjAuthorized->{$project_name});
		$is_ok = 1;
		my @lpat_attic;
		unless (exists $hProjectUser->{$project_name}->{ped}) {
			my ($lHPat, $hPed) = pedigree_details($hProjectUser->{$project_name}->{id});
			$hProjectUser->{$project_name}->{ped}->{list_hash} = $lHPat;
			$hProjectUser->{$project_name}->{ped}->{hash_ped} = $hPed;
		}
		my $analyse = $hProjAuthorized->{$project_name};
		my $hAllPatNames;
		foreach my $h (@{$hProjectUser->{$project_name}->{ped}->{list_hash}}) { $hAllPatNames->{$h->{name}} = undef; }
		foreach my $pat_he_ho (split(';', $h_dejavu->{$project_name}->{string})) {
			my ($pat_name, $he_ho) = split(':', $pat_he_ho);
			next if ($only_patients and not exists $h_only_patients->{$pat_name});
			delete $hAllPatNames->{$pat_name};
			foreach my $h (@{$hProjectUser->{$project_name}->{ped}->{list_hash}}) { 
				if ($h->{name} eq $pat_name) {
					$hProjDejavu_user->{$analyse}->{$project_name}->{$pat_name}->{status} = $h->{status};
					$hProjDejavu_user->{$analyse}->{$project_name}->{$pat_name}->{family} = $h->{family};
					my $image;
					if (exists $hProjectUser->{$project_name}->{ped}->{hash_ped}->{$h->{family}}->{father} and $hProjectUser->{$project_name}->{ped}->{hash_ped}->{$h->{family}}->{father} eq $pat_name) {
						$image = "<img src='http://www.polyweb.fr/icons/Polyicons/male.png'>";
					}
					elsif (exists $hProjectUser->{$project_name}->{ped}->{hash_ped}->{$h->{family}}->{mother} and $hProjectUser->{$project_name}->{ped}->{hash_ped}->{$h->{family}}->{mother} eq $pat_name) {
						$image = "<img src='http://www.polyweb.fr/icons/Polyicons/female.png'>";
					}
					else {
						if   ($hProjectUser->{$project_name}->{ped}->{hash_ped}->{$h->{family}}->{all}->{$pat_name}->{sex} == 1) { $image = "<img src='http://www.polyweb.fr/icons/Polyicons/baby-boy.png'>"; }
						elsif($hProjectUser->{$project_name}->{ped}->{hash_ped}->{$h->{family}}->{all}->{$pat_name}->{sex} == 2) { $image = "<img src='http://www.polyweb.fr/icons/Polyicons/baby-girl.png'>"; }
					}
					$hProjDejavu_user->{$analyse}->{$project_name}->{$pat_name}->{img_child_mother_father} = $image;
				}
			}
			$he_ho = 'he' if ($he_ho eq '2');
			$he_ho = 'ho' if ($he_ho eq '1');
			$hProjDejavu_user->{$analyse}->{$project_name}->{$pat_name}->{he_ho} = $he_ho;
		}
		$h_var_info->{polyquery_patients_attic}->{$project_name} = join('+', keys %$hAllPatNames);
	}
	return unless ($is_ok);
	next unless ($var->isDM());
	
	my $h_no_hgmd = $project->buffer->queryHgmd->getDataHGMDPro($var->hgmd_id());
	$h_var_info->{type} = $var->type_public_db();
	$h_var_info->{hgmd_id} = $var->hgmd_id();
	$h_var_info->{rs_name} = $var->rs_name();
	$h_var_info->{dejavu} = 'Proj: '.$nb_project.' / Pat: '.$nb_patients;
	my $annot;
	foreach my $gene (@{$var->getGenes()}) {
		unless (exists $h_var_info->{gene}) { $h_var_info->{gene} = $gene->external_name(); }
		else { $h_var_info->{gene} .= '+'.$gene->external_name(); }
		if ($annot) { $annot .= ";".$gene->external_name().': '.$var->variationTypeInterface($gene); }
		else { $annot = $gene->external_name().': '.$var->variationTypeInterface($gene); }
	}
	$h_var_info->{annot} = $annot;
	$h_var_info->{hgmd_new} = $var->isNewHgmd();
	$h_var_info->{hgmd_hgvs} = $var->hgmd_hgvs();
	$h_var_info->{hgmd_prot} = $h_no_hgmd->{INFO}->{PROT};
	$h_var_info->{hgmd_class} = $var->hgmd_class();
	if ($var->hgmd()) {
		$h_var_info->{hgmd_phen} = $var->hgmd_phenotype();
	}
	$h_var_info->{hgmd_disease} = $var->hgmd_disease();
	$h_var_info->{hgmd_ad_ar} = $var->hgmd_inheritance();
	$h_var_info->{omimid} = $h_no_hgmd->{PRO}->{omimid};
	foreach my $pubmed_id (keys %{$var->hgmd_pubmed()}) {
		$h_var_info->{hgmd_pubmed}->{$pubmed_id}->{title} = $var->hgmd_pubmed->{$pubmed_id}->{title};
		$h_var_info->{hgmd_pubmed}->{$pubmed_id}->{author} = $var->hgmd_pubmed->{$pubmed_id}->{author};
		$h_var_info->{hgmd_pubmed}->{$pubmed_id}->{year} = $var->hgmd_pubmed->{$pubmed_id}->{year};
	}
	
	my $h_hgmd_version_current = $project->buffer->config->{'hgmd_db_current_version'};
	$h_var_info->{hgmd_versions} = '(current) '.join(', ', sort values %{$project->buffer->config->{'hgmd_db_current_version'}});
	if (exists $project->buffer->queryHgmd->hashOldAccNum->{$var->hgmd_id()}) {
		$h_var_info->{hgmd_versions} .= ', '.join(', ', sort keys %{$project->buffer->queryHgmd->hashOldAccNum->{$var->hgmd_id()}})
	}
	return ($h_var_info, $hProjDejavu_user);
}

sub pedigree_details {
	my $project_id = shift;
	my $listHash = $query->getPatients($project_id);
  	my ($hTmp, $hPed);
  	foreach my $hFam (@$listHash) {
  		my $fam_name = $hFam->{family};
  		$hTmp->{$fam_name}->{father} = $hFam->{father} if ($hFam->{father});
  		$hTmp->{$fam_name}->{mother} = $hFam->{mother} if ($hFam->{mother});
  		$hTmp->{$fam_name}->{'all'}->{$hFam->{name}}->{'sex'} = $hFam->{sex};
  		$hTmp->{$fam_name}->{'all'}->{$hFam->{name}}->{'status'} = $hFam->{status};
  	}
  	foreach my $fam_name (keys %$hTmp) {
  		if (exists $hTmp->{$fam_name}->{father}) { $hPed->{$fam_name}->{father} = $hTmp->{$fam_name}->{father}; }
  		if (exists $hTmp->{$fam_name}->{mother}) { $hPed->{$fam_name}->{mother} = $hTmp->{$fam_name}->{mother}; }
  		foreach my $pat_name (keys %{$hTmp->{$fam_name}->{all}}) {
  			if (exists $hTmp->{$fam_name}->{father} and $hTmp->{$fam_name}->{father} eq $pat_name) {
  				$hPed->{$fam_name}->{'all'}->{$pat_name}->{'sex'} = $hTmp->{$fam_name}->{all}->{$pat_name}->{sex};
		  		$hPed->{$fam_name}->{'all'}->{$pat_name}->{'status'} = $hTmp->{$fam_name}->{all}->{$pat_name}->{status};
  			}
  			elsif (exists $hTmp->{$fam_name}->{mother} and $hTmp->{$fam_name}->{mother} eq $pat_name) {
  				$hPed->{$fam_name}->{'all'}->{$pat_name}->{'sex'} = $hTmp->{$fam_name}->{all}->{$pat_name}->{sex};
		  		$hPed->{$fam_name}->{'all'}->{$pat_name}->{'status'} = $hTmp->{$fam_name}->{all}->{$pat_name}->{status};
  			}
  			else {
  				$hPed->{$fam_name}->{'children'}->{$pat_name} = undef;
  				$hPed->{$fam_name}->{'all'}->{$pat_name}->{'sex'} = $hTmp->{$fam_name}->{all}->{$pat_name}->{sex};
		  		$hPed->{$fam_name}->{'all'}->{$pat_name}->{'status'} = $hTmp->{$fam_name}->{all}->{$pat_name}->{status};
  			}
  		}
  	}
#	warn Dumper $hPed;
#	die;
	return ($listHash, $hPed);
}

sub get_line_variant {
	my ($no, $var_id) = @_;
	my $out;
	my $h_dejavu = $project->getDejaVuInfos($var_id);
	next unless ($h_dejavu);
	my $is_ok;
	my $nb_dejavu_my_projects = 0;
	my $nb_dejavu_my_patients = 0;
	my $nb_dejavu_my_patients_he = 0;
	my $nb_dejavu_my_patients_ho = 0;
	my $nb_dejavu_my_exome = 0;
	my $nb_dejavu_my_genome = 0;
	my $nb_dejavu_my_diag = 0;
	my $last_project_name;
	my $resume_dejavu;
	foreach my $project_name (keys %{$h_dejavu}) {
		next unless exists ($hProjAuthorized->{$project_name});
		next unless ($project_name =~ /NGS20/);
		$is_ok = 1;
		my $analyse = 'diag';
		$nb_dejavu_my_projects++;
		$nb_dejavu_my_patients += $h_dejavu->{$project_name}->{'nb'};
		$nb_dejavu_my_patients_he += $h_dejavu->{$project_name}->{'he'};
		$nb_dejavu_my_patients_ho += $h_dejavu->{$project_name}->{'ho'};
		if    ($analyse eq 'exome')  { $nb_dejavu_my_exome++; }
		elsif ($analyse eq 'genome') { $nb_dejavu_my_genome++; }
		elsif ($analyse eq 'diag')   { $nb_dejavu_my_diag++; }
		if ($resume_dejavu) { $resume_dejavu .= '|'; }
		$resume_dejavu .= $project_name.';'.$h_dejavu->{$project_name}->{'patients'};
		$last_project_name = $project_name;
	}
	if ($is_ok) {
		my $var = $project->_newVariant($var_id);
		my $type = $var->type_public_db();
		my $h_no_hgmd = $no->get($var_id);
		my $hgmd_id = $h_no_hgmd->{ID};
		$out .= qq{<tr>};
		$out .= "<td align='center' style='font-size:10px;'>".$var_id;
		if ($var->rs_name) { $out .= "<br>".$var->rs_name(); }
		$out .= "<br>".$hgmd_id."</td>";
		
		my $dejavu_url = "../polydejavu/dejavu.html?input=".$var_id;
		my $text_dejavu = "<a href='".$dejavu_url."' target='_blank'>"."Diag:$nb_dejavu_my_diag Others:".($nb_dejavu_my_exome+$nb_dejavu_my_genome)." He:$nb_dejavu_my_patients_he Ho:$nb_dejavu_my_patients_ho</a>";
		$out .= "<td align='center' style='font-size:10px;'>".$text_dejavu."</td>";
		
		my $annot;
		foreach my $gene (@{$var->getGenes()}) {
			if ($annot) { $annot .= "<br>".$gene->external_name().': '.$var->variationTypeInterface($gene); }
			else { $annot = $gene->external_name().': '.$var->variationTypeInterface($gene); }
		}
	
		$out .= "<td align='center' style='font-size:10px;'>".$annot."</td>";
		$out .= "<td align='center' style='font-size:10px;'>".$h_no_hgmd->{INFO}->{DNA}."<br>".$h_no_hgmd->{INFO}->{PROT}."</td>";
		$out .= "<td align='center' style='font-size:10px;'>".$h_no_hgmd->{INFO}->{CLASS}."</td>";
		
		my $text_publi = "<a href='https://www.ncbi.nlm.nih.gov/pubmed/?term=".$h_no_hgmd->{PRO}->{pmid}."' target='_blank'>".'<i><u>[Pubmed] '.$h_no_hgmd->{PRO}->{title}."</i></u></a>";
		
		my $text_phen =	$h_no_hgmd->{INFO}->{PHEN};
		$text_phen   .= "<br>";
		$text_phen   .= qq{<a href="" target="_blank"><span id="span_hgmd_pub_link" onClick="alert('link publi')">$text_publi</span></a>};
		$out .= "<td align='center' style='font-size:10px;'>".$text_phen."</td>";
		
		$out .= qq{</tr>};
	}
	return $out;
}

sub print_cgi_header_html {
	my ($cgi) = @_;
	print $cgi -> header;
    print qq{
    	<script type="text/javascript" src="//maxcdn.bootstrapcdn.com/bootstrap/3.3.2/js/bootstrap.min.js""></script>
    };
	print qq{		
		<style>
			table-mybordered > thead > tr > th, .table-mybordered > thead > tr > th, table-mybordered > tbody > tr > th, .table-mybordered > tbody > tr > th, table-mybordered > tfoot > tr > th, .table-mybordered > tfoot > tr > th, table-mybordered > thead > tr > td, .table-mybordered > thead > tr > td, table-mybordered > tbody > tr > td, .table-mybordered > tbody > tr > td, table-mybordered > tfoot > tr > td, .table-mybordered > tfoot > tr > td {
				border: 1px solid #95A5A6;
		  		vertical-align:middle;
		   		text-align:center; 
		   		padding:2px;
			}
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
			.badge1 {
			  padding: 1px 9px 2px;
			  font-size: 12.025px;
			  font-weight: bold;
			  white-space: nowrap;
			  color: #ffffff;
			  background-color: #999999;
			  -webkit-border-radius: 9px;
			  -moz-border-radius: 9px;
			  border-radius: 9px;
			}
			
			.hide-bullets {
			    list-style:none;
			    margin-left: -40px;
			    margin-top:20px;
			}
			
			.thumbnail {
			    padding: 0;
			}
			
			.carousel-inner>.item>img, .carousel-inner>.item>a>img {
			    width: 100%;
			}
		</style>	
	}	
}

