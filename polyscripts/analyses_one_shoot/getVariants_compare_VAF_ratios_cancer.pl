#!/usr/bin/perl

#/data-isilon/bipd-src/mbras/git_repository/polyscripts/analyses_one_shoot/getVariants_compare_VAF_ratios_cancer.pl -p=NGS2022_6029 -ph=HPO95_97_E_T -pi=HPO95_97_A_T_1 -max_vaf_healthy=60 -min_vaf_ill=65 -o=/data-isilon/bipd-src/mbras/

#/data-isilon/bipd-src/mbras/git_repository/polyscripts/analyses_one_shoot/getVariants_compare_VAF_ratios_cancer.pl -p=NGS2022_6029 -ph=HPO7679_E_T -pi=HPO76_79_A_T -max_vaf_healthy=60 -min_vaf_ill=65 -o=/data-isilon/bipd-src/mbras/

#/data-isilon/bipd-src/mbras/git_repository/polyscripts/analyses_one_shoot/getVariants_compare_VAF_ratios_cancer.pl -p=NGS2022_6029 -ph=HPO69_E_c -pi=HPO69_A_T -max_vaf_healthy=60 -min_vaf_ill=65 -o=/data-isilon/bipd-src/mbras/

use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/";
use lib "$Bin/../../GenBo/lib/obj-nodb";
use lib "$Bin/../../GenBo/lib/obj-nodb/packages";
use strict;
use QueryOnlyVcf;
use Getopt::Long;
use Data::Dumper;
use String::ProgressBar;
use Compress::Snappy;
use Parallel::ForkManager;
use session_export;
use xls_export;
use Storable qw/thaw freeze/;

my $fork = 1;
my ($use_freq_1_1000, $use_freq_1_100, $max_gnomad_ho);
my ($project_name, $patient_name_healthy, $patient_name_ill,$max_vaf_healthy, $min_vaf_ill, $outdir);
GetOptions(
	'project|p=s' => \$project_name,
	'patient_healthy|ph=s' => \$patient_name_healthy,
	'patient_ill|pi=s' => \$patient_name_ill,
	'max_vaf_healthy=s' => \$max_vaf_healthy,
	'min_vaf_ill=s' => \$min_vaf_ill,
	
	'max_gnomad_ho=s' => \$max_gnomad_ho,
	'use_freq_1_100=s' => \$use_freq_1_100,
	'use_freq_1_1000=s' => \$use_freq_1_1000,
	
	'outdir|o=s' => \$outdir,
	'fork|f=s' => \$fork,
);
die "\n\nNo -project or -p option... Die...\n\n" unless ($project_name);
die "\n\nNo -patient_healthy or -ph option... Die...\n\n" unless ($patient_name_healthy);
die "\n\nNo -patient_ill or -pi option... Die...\n\n" unless ($patient_name_ill);
die "\n\nNo -max_vaf_healthy  option... Die...\n\n" unless ($max_vaf_healthy);
die "\n\nNo -min_vaf_ill option... Die...\n\n" unless ($min_vaf_ill);
die "\n\nNo -outdir or -o option... Die...\n\n" unless ($outdir);


my $buffer = GBuffer->new;
my $project = $buffer->newProjectCache( -name => $project_name );
my $xls_export = new xls_export();

my $patient_healthy = $project->getPatient($patient_name_healthy);
my $patient_ill = $project->getPatient($patient_name_ill);

my (@lVar, @lPat);
push(@lPat, $patient_healthy);
push(@lPat, $patient_ill);

foreach my $chr (@{$project->getChromosomes()}) {
	my $res;
	my $vector_common = $chr->getNewVector();
	$vector_common->Intersection($patient_healthy->getVectorHe($chr), $patient_ill->getVectorHe($chr));
	my @lVarLocal;
	foreach my $var (@{$chr->getListVarObjects($vector_common)}) {
		if ($use_freq_1_1000) {
			next if $var->frequency() > 0.0001;
		}
		elsif ($use_freq_1_100) {
			next if $var->frequency() > 0.001;
		}
		my $vaf_p_healthy = $var->getRatio($patient_healthy);
		next if ($vaf_p_healthy > $max_vaf_healthy);
		my $vaf_p_ill = $var->getRatio($patient_ill);
		next if ($vaf_p_ill < $min_vaf_ill);
		if ($max_gnomad_ho) {
			my $gho = $var->getGnomadHO();
			if ($gho =~ /0-9/) {
				next if $gho > $max_gnomad_ho;
			}
		}
		push(@lVar, $var);
	}
}

$xls_export->output_dir($outdir);
my $title_ext;
$title_ext .= "_max_freq_1_100" if $use_freq_1_100;
$title_ext .= "_max_freq_1_1000" if $use_freq_1_1000;
$title_ext .= "_max_gnomad_ho_$max_gnomad_ho" if $max_gnomad_ho;
my $title_xls = $project->name().'_VAF_analyse_'.$patient_ill->name().'_MinVaf_'.$min_vaf_ill.'_vs_'.$patient_healthy->name().'_MaxVaf_'.$max_vaf_healthy.$title_ext.'.xls';
$xls_export->title_page($title_xls);
#$xls_export->store_variants_infos(\@lVar, $project);
$xls_export->store_variants_infos(\@lVar, $project, \@lPat);
my ($list_datas_annotations) = $xls_export->prepare_generic_datas_variants();
my ($list_datas_annotations_cnvs) = $xls_export->prepare_generic_datas_cnvs();
$xls_export->add_page_merged('Variants Merged', $xls_export->list_generic_header(), $list_datas_annotations);
$xls_export->add_page('Variants Not Merged', $xls_export->list_generic_header(), $list_datas_annotations);
if (scalar @$list_datas_annotations_cnvs > 0) {
	$xls_export->add_page('Cnvs', $xls_export->list_generic_header_cnvs(), $list_datas_annotations_cnvs);
}
$xls_export->export();
exit(0);




