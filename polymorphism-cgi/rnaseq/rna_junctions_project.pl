#!/usr/bin/perl
$|=1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo";
use lib "$Bin/../../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";

use connect;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use JSON;
use xls_export;
use session_export;

my $cgi    = new CGI;
my $project_name = $cgi->param('project');

my $buffer = GBuffer->new;
my $project = $buffer->newProjectCache( -name => $project_name );

print $cgi->header('text/json-comment-filtered');
print "{\"progress\":\".";

my $hType_patients;
$hType_patients = $project->get_hash_patients_description_rna_seq_junction_analyse() if (-d $project->get_path_rna_seq_junctions_analyse_description_root());

my $h_resume;
foreach my $patient (@{$project->getPatients}) {
	#$patient->use_not_filtred_junction_files(0);
	if (($hType_patients and exists $hType_patients->{$patient->name()}->{pat}) or not $hType_patients) {
		$h_resume->{$patient->name()}->{nb_junctions_all} = 0;
		$h_resume->{$patient->name()}->{nb_junctions_ri} = 0;
		$h_resume->{$patient->name()}->{nb_junctions_se} = 0;
		foreach my $chr (@{$project->getChromosomes()}) {
			$h_resume->{$patient->name()}->{nb_junctions_all} += $chr->countThisVariants($patient->getJunctionsVector($chr));
			$h_resume->{$patient->name()}->{nb_junctions_ri} += $chr->countThisVariants($patient->getVectorJunctionsRI($chr)) if ($patient->getVectorJunctionsRI($chr));
			$h_resume->{$patient->name()}->{nb_junctions_se} += $chr->countThisVariants($patient->getVectorJunctionsSE($chr)) if ($patient->getVectorJunctionsSE($chr));
		}
	}
	else {
		$h_resume->{$patient->name()}->{nb_junctions_all} = 'CONTROL';
		$h_resume->{$patient->name()}->{nb_junctions_ri} = 'CONTROL';
		$h_resume->{$patient->name()}->{nb_junctions_se} = 'CONTROL';
	}
	$h_resume->{$patient->name()}->{bam_file} = $patient->bamUrl();
	print '.';
}

my $html_table_patients = qq{<center><span style='font-size:17px;color:#23527C;'>PROJECT RESUME</span></center><br><table data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="10" data-page-list="[10, 20, 50, 100, 200, 300]" data-resizable='true' id='table_id_patients' class='table table-striped' style='font-size:13px;'>};
$html_table_patients .= qq{<thead>};
$html_table_patients .= qq{<th data-field="name" data-sortable="true"><center><b>Patient Name</b></center></th>};
$html_table_patients .= qq{<th data-field="nb_junctions_all" data-sortable="true"><center><b>Nb Junctions ALL</b></center></th>};
$html_table_patients .= qq{<th data-field="nb_junctions_ri" data-sortable="true"><center><b>Nb Junctions RI</b></center></th>};
$html_table_patients .= qq{<th data-field="nb_junctions_se" data-sortable="true"><center><b>Nb Junctions SE</b></center></th>};
$html_table_patients .= qq{<th data-field="igv_bam" data-sortable="true"><center><b>IGV</b></center></th>};
$html_table_patients .= qq{</thead>};
$html_table_patients .= qq{<tbody>};

my (@l_tr_patient, @l_tr_control);
foreach my $patient_name (sort keys %$h_resume) {
	my $all = $h_resume->{$patient_name}->{nb_junctions_all};
	
	my $s1000 = $h_resume->{$patient_name}->{nb_junctions_score_1000};
	my $s100 = $h_resume->{$patient_name}->{nb_junctions_score_100};
	my $s10 = $h_resume->{$patient_name}->{nb_junctions_score_10};
	my $s1 = $h_resume->{$patient_name}->{nb_junctions_score_1};
	
	my $ri = $h_resume->{$patient_name}->{nb_junctions_ri};
	my $se = $h_resume->{$patient_name}->{nb_junctions_se};
	my $ri_perc = '0%';
	$ri_perc = int(($ri / $all) * 100).'%' if $ri > 0;
	my $se_perc = '0%';
	$se_perc = int(($se / $all) * 100).'%' if $se > 0;
	my $bam_file = $h_resume->{$patient_name}->{bam_file};
	my $tr = qq{<tr style="text-align:center;font-size:12px;">};
	$tr .= qq{<td>$patient_name</td>};
	$tr .= qq{<td><b">$all</b></td>};
	$tr .= qq{<td>$ri ($ri_perc)</td>};
	$tr .= qq{<td>$se ($se_perc)</td>};
	
	my $bam_pat_file = $h_resume->{$patient_name}->{bam_file};
	my $cmd_bam = qq{onclick="add_bam_igv('$bam_pat_file')";};
	$tr .= qq{<td><button $cmd_bam>Add BAM IGV</button></td>};
	$tr .= qq{</tr>};
	if (lc($all) =~ /control/) { push(@l_tr_control, $tr); }
	else { push(@l_tr_patient, $tr); }
}
foreach my $this_tr (@l_tr_patient) { $html_table_patients .= $this_tr; }
foreach my $this_tr (@l_tr_control) { $html_table_patients .= $this_tr; }
$html_table_patients .= qq{</tbody>};
$html_table_patients .= qq{</table>};
print '.';

my $hash;
$hash->{html} = $html_table_patients;
if ($hType_patients) {
	my @lPatOk;
	foreach my $patname (sort keys %$hType_patients) {
		next if (not exists $h_resume->{$patname});
		push(@lPatOk, $patname) if exists $hType_patients->{$patname}->{pat};
	}
	$hash->{patients} = join(';', @lPatOk);
}
else {
	$hash->{patients} = join(';', sort keys %$h_resume);
}

printJson($hash);
exit(0);



sub printJson {
	my ($hashRes) = @_;
	my $json_encode = encode_json $hashRes;
	print ".\",";
	$json_encode =~ s/{//;
	print $json_encode;
	exit(0);
}
