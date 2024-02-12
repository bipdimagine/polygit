#!/usr/bin/perl
$|=1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo";
use lib "$Bin/../../../../GenBo/lib/obj-nodb";
use lib "$Bin/../../../../GenBo/lib/obj-nodb/packages";
require "$Bin/../../../../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";

use connect;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use JSON;
use xls_export;
use session_export;
use List::MoreUtils qw{ natatime };
use Parallel::ForkManager;
use JSON::XS;
use DateTime;


my ($project_name, $force);
GetOptions(
	'project=s' => \$project_name,
	'force=s' => \$force,
);

my $buffer = GBuffer->new;
my $project = $buffer->newProject( -name => $project_name );

my $project_path = $project->getProjectPath();
my $json_analise_file = $project_path.'/'.$project_name.'.splices.analyse.json';
if (-e $json_analise_file and not $force) {
	print $json_analise_file;
	exit(0);
}

my $hash;

# ANALYSE
my $dt = DateTime->now;
$hash->{analyse}->{date} = $dt->dmy;
$hash->{analyse}->{path} = $project_path.'/analysis/';

# PROJECT
$hash->{project}->{gencode}->{release} = $project->getVersion();
if ($hash->{project}->{gencode}->{release} =~ /HG19/) {
	$hash->{project}->{gencode}->{version} = $project->gencode_version();
}
elsif ($hash->{project}->{gencode}->{release} eq 'MM38') {
	$hash->{project}->{gencode}->{version} = 'M25';
}
elsif ($hash->{project}->{gencode}->{release} eq 'MM39') {
	$hash->{project}->{gencode}->{version} = 'M32';
}
$hash->{project}->{gencode}->{gtf} = $project->gtf_file();
$hash->{project}->{gencode}->{rds} = $project->rds_gencode_file();
$hash->{project}->{gencode}->{junctions_canoniques_rds} = $project->rds_junctions_canoniques_gencode_file();

# SOFTWARES
$hash->{softwares}->{samtools} = $project->getSoftware('samtools');
$hash->{softwares}->{sambamba} = $project->getSoftware('sambamba');
$hash->{softwares}->{bedtools} = $project->getSoftware('bedtools');
$hash->{softwares}->{picard} = $project->getSoftware('picard');

# POLYRNASEQSEA PATHS
$hash->{rnaseqseapaths}->{script} = "$Bin/junctions/RNAseqSEA_capt_js_dev.r";
$hash->{rnaseqseapaths}->{biblio} = "$Bin/junctions/biblio_RNAseqSEA_capt_js.r";

# PATIENTS
foreach my $patient (@{$project->getPatients()}) {
	$hash->{inputs}->{bam}->{$patient->name()} = $patient->getBamFile();
}

my $json_encode = encode_json $hash;
open (JSON, ">$json_analise_file");
print JSON $json_encode;
close (JSON);

`chmod 775 $json_analise_file`;

print $json_analise_file;

exit(0);

