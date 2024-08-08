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
use GBuffer;
use export_data;
use JSON;
use polyweb_dude;
use VcfMerge;
use GenBoNoSql;
use Set::IntervalTree;
use Spreadsheet::WriteExcel;
use Bit::Vector;
use Bit::Vector::Overload;
use Compress::Snappy;
use Storable qw(store retrieve freeze dclone thaw);
use POSIX qw(strftime);
use List::MoreUtils qw(natatime);
use CGI::Session;
use html; 
use Carp;
use Cache_Commons;
use QueryVectorFilter;
use IO::Handle;
use xls_export;
use session_export;


my $cgi = new CGI();
my $use_session_id = $cgi->param('session_id');
my $export_xls = $cgi->param('export_xls');
my $merged = $cgi->param('merged');
my $only_transcript = $cgi->param('only_transcript');


if ($only_transcript) { $only_transcript =~ s/_.+//; }

return load_xls($use_session_id) if ($export_xls);

my $session = new session_export();
$session->load_session( $use_session_id );

my $hRes;
$hRes->{'html_variants'} = 'Session expired... Please, Rrload your filters...';

warn $session->load('html_variants');

eval { $hRes->{'html_gene'} = $session->load('html_gene'); };
if ($@) { warn 'no html_gene'; }
eval { $hRes->{'html_variants'} = $session->load('html_variants'); };
if ($@) { warn 'no html_variants'; }
eval { $hRes->{'html_page_title'} = $session->load('html_page_title'); };
if ($@) { warn 'no html_page_title'; }
eval { $hRes->{'html_source'} = $session->load('html_source'); };
if ($@) { warn 'no html_source'; }
eval { $hRes->{'html_title'} = $session->load('html_title'); };
if ($@) { warn 'no html_title'; }
eval { $hRes->{'hash_filters'} = $session->load('hash_filters'); };
if ($@) { warn 'no hash_filters'; }
eval { $hRes->{'gencode_version'} = $session->load('gencode_version'); };
if ($@) { warn 'no gencode_version'; }

unless ($hRes->{'html_variants'}) {
	$hRes->{'html_variants'} = '<div>Session expired... Please, reload your filters...</div>';
	$hRes->{'html_page_title'} = 'Expired';
	$hRes->{'html_title'} = 'Expired';
}

print $cgi->header('text/json-comment-filtered');
my $json_encode = encode_json $hRes;
print $json_encode;
exit(0);





sub load_xls() {
	my $use_session_id = shift;
	my $xls_export = new xls_export();
	$xls_export->load($use_session_id);
	my ($list_datas_annotations) = $xls_export->prepare_generic_datas_variants();
	my ($list_datas_annotations_cnvs) = $xls_export->prepare_generic_datas_cnvs();
	my ($h_by_patients, @list_datas_patients);
	eval { $h_by_patients = $xls_export->get_specific_infos_stored('projects_patients_infos'); };
	if ($@) {}
	else {
		foreach my $var_id (keys %{$h_by_patients}) {
			foreach my $project_name (keys %{$h_by_patients->{$var_id}}) {
				foreach my $patient_name (keys %{$h_by_patients->{$var_id}->{$project_name}}) {
					my $h;
					$h->{'variation'} = $h_by_patients->{$var_id}->{$project_name}->{$patient_name}->{'variation'};
					$h->{'project'} = $h_by_patients->{$var_id}->{$project_name}->{$patient_name}->{'project'};
					$h->{'description'} = $h_by_patients->{$var_id}->{$project_name}->{$patient_name}->{'description'};
					$h->{'phenotypes'} = $h_by_patients->{$var_id}->{$project_name}->{$patient_name}->{'phenotypes'};
					$h->{'perc'} = $h_by_patients->{$var_id}->{$project_name}->{$patient_name}->{'perc'};
					$h->{'patient'} = $h_by_patients->{$var_id}->{$project_name}->{$patient_name}->{'name'};
					$h->{'family'} = $h_by_patients->{$var_id}->{$project_name}->{$patient_name}->{'fam'};
					$h->{'sex_status_icon'} = $h_by_patients->{$var_id}->{$project_name}->{$patient_name}->{'sex_status_icon'};
					$h->{'sex'} = $h_by_patients->{$var_id}->{$project_name}->{$patient_name}->{'sex'};
					$h->{'status'} = $h_by_patients->{$var_id}->{$project_name}->{$patient_name}->{'status'};
					$h->{'parent_child'} = $h_by_patients->{$var_id}->{$project_name}->{$patient_name}->{'parent_child'};
					$h->{'model'} = '-';
					$h->{'model'} = $h_by_patients->{$var_id}->{$project_name}->{$patient_name}->{'model'} if (lc($h->{'parent_child'}) eq 'child');
					$h->{'he_ho'} = $h_by_patients->{$var_id}->{$project_name}->{$patient_name}->{'he_ho'};
					#$self->{hash_variants_global}->{$chr_h_id}->{$var_id}->{'patients'}->{$patient->name()}->{'he_ho'}
					push(@list_datas_patients, $h);
				}
			}
		}
	}
	
	my ($h_pubmed, @list_datas_pubmed);
	eval { $h_pubmed = $xls_export->get_specific_infos_stored('variants_pubmed'); };
	if ($@) {}
	else {
		foreach my $var_id (keys %{$h_pubmed}) {
			foreach my $pubmed_id (keys %{$h_pubmed->{$var_id}}) {
				my $h;
				$h->{'variation'} = $var_id;
				$h->{'url'} = $h_pubmed->{$var_id}->{$pubmed_id}->{url};
				$h->{'title'} = $h_pubmed->{$var_id}->{$pubmed_id}->{title};
				push(@list_datas_pubmed, $h);
			}
		}
	}
	my (@list_datas_annotations_with_patients, $hProjFound, $h_by_tr);
	my @lHeaderWithProj = @{$xls_export->list_generic_header()};
	
	my $hdone_tr;
	my $isHGMD_export;
	my $nb_v = 0;
	foreach my $hvar (@$list_datas_annotations) {
		$nb_v++;
		my $enst = $hvar->{transcript};
		my ($hproj, $hpat);
		my ($var_id, $rs) = split(' ', $hvar->{variation});
		foreach my $proj_name (keys %{$h_by_patients->{$var_id}}) {
			foreach my $pat_name (keys %{$h_by_patients->{$var_id}->{$proj_name}}) {
				$isHGMD_export = 1 if ($proj_name eq 'HGMD' and $pat_name eq 'HGMD');
				my $hvar_new = dclone($hvar);
				$hvar_new->{project} = $proj_name;
				$hvar_new->{patient} = $pat_name;
				$hvar_new->{sex} = 'male';
				$hvar_new->{sex} = 'female' if $h_by_patients->{$var_id}->{$proj_name}->{$pat_name}->{sex} == 2;
				$hvar_new->{status} = $h_by_patients->{$var_id}->{$proj_name}->{$pat_name}->{status};
				$hvar_new->{perc} = $h_by_patients->{$var_id}->{$proj_name}->{$pat_name}->{perc};
				$hvar_new->{model} = '-';
				$hvar_new->{model} = $h_by_patients->{$var_id}->{$proj_name}->{$pat_name}->{model} if $h_by_patients->{$var_id}->{$proj_name}->{$pat_name}->{model} ne '?';
				$hvar_new->{he_ho} = '-';
				$hvar_new->{he_ho} = $h_by_patients->{$var_id}->{$proj_name}->{$pat_name}->{he_ho} if $h_by_patients->{$var_id}->{$proj_name}->{$pat_name}->{he_ho};
				push(@list_datas_annotations_with_patients,  dclone($hvar_new));
				push(@{$h_by_tr->{$enst}},  dclone($hvar_new));
			}
		}
	}
	
	if ($merged) {
		$xls_export->add_page_merged('Variants Merged', $xls_export->list_generic_header(), \@list_datas_annotations_with_patients);
	}
	else {
		my @lHeaderWithPat = @{$xls_export->list_generic_header()};
		push (@lHeaderWithPat, 'Project', 'Patient', 'Sex', 'Status', 'Perc', 'He_Ho', 'Model') if (not $isHGMD_export);
		if (not $only_transcript and $nb_v > 5000) {
			$xls_export->add_page('ALL TRANSCRIPTS', \@lHeaderWithPat, \@list_datas_annotations_with_patients);
		}
		foreach my $enst (sort keys %$h_by_tr) {
			next if ($only_transcript and not $enst =~ /$only_transcript/);
			next if $nb_v > 5000;
			$xls_export->add_page($enst, \@lHeaderWithPat, $h_by_tr->{$enst});
		}
	}
	if (scalar @$list_datas_annotations_cnvs > 0) {
		$xls_export->add_page('Cnvs', $xls_export->list_generic_header_cnvs(), $list_datas_annotations_cnvs);
	}
	
	my @lLinesHeaderPatients = ('Variation', 'Project', 'Description', 'Family', 'Patient', 'Parent_child', 'Sex', 'Status', 'Perc', 'He_Ho', 'Model', 'Phenotypes');
	$xls_export->add_page('Patients', \@lLinesHeaderPatients, \@list_datas_patients) if @list_datas_patients;
	
	my @lLinesHeaderPubmed = ('Variation', 'Url', 'Title');
	$xls_export->add_page('Pubmed', \@lLinesHeaderPubmed, \@list_datas_pubmed) if @list_datas_pubmed;
	
	$xls_export->export();
	exit(0);
}
