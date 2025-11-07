#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/kyoto";

#use lib "/bip-d/soft/distrib/tabix/latest/perl";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/coverage";
use lib "$Bin/../packages/validation_variation"; 
use lib "$Bin/../packages/cache"; 
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/utility.pm";

use draw_cnv; 
use html; 
use infos_coverage_exons;
use JSON::XS;
use image_coverage;
#use Set::;
use Carp;
use export_data;
use strict;
use Set::IntSpan::Fast::XS;
use Data::Dumper;
use GBuffer;
use Getopt::Long;
use Carp;
use Set::Intersection;
use Tabix;
use Storable qw/thaw/;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use coverage;
use Spreadsheet::WriteExcel;
use POSIX;
use validationQuery;
use Date::Tiny;
use JSON::XS;
use List::MoreUtils qw{part};
use URI;

#
#
my $buffer = GBuffer->new();
my $cgi          = new CGI();
my $project_name = $cgi->param('project');
my $patient_name = $cgi->param('patient');
my $position = $cgi->param('position');
#
my $project = $buffer->newProject(-name=>$project_name);;
my $patients = $project->get_only_list_patients($patient_name);
#my $patients = $project->getPatients();
my @bams;
my $server = 'https://'.$ENV{HTTP_HOST};
if ($server eq "10.200.27.108:8010") {
	$server = "defidiag.polyweb.fr"
}
elsif ($server eq "10.200.27.108") {
	$server = "www.polyweb.fr"
}
elsif ($server eq "10.200.27.103") {
	$server = "darwin.bipd.fr"
}

foreach my $p (@$patients){
	 my $z = "{";
	 	my $name = $p->name;
	 	my $bam = $p->bamUrl();
		$z .= qq{"name":"$name",} .qq{ type: "alignment",format: "bam"};
		$z.= qq{,"url":"http://$server/$bam"};
		$z.= qq{,"indexURL":"http://$server/$bam.bai"};
          $z .= "}"; 
		push(@bams,$z);
	
}

my $text_track_bam = join(",\n",@bams);



#my $t = $ENV{'HTTP_X_FORWARDED_FOR'};#." + ".$ENV{HTTP_REFERER}." ".$cgi->self_url;

print $cgi -> header;
print qq{
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no">
    <meta name="description" content="">
    <meta name="author" content="">
    <link rel="shortcut icon" href="https://igv.org/web/img/favicon.ico">

    <title>igv.js</title>

    <!-- IGV JS-->
    <script src="https://cdn.jsdelivr.net/npm/igv@3.0.2/dist/igv.min.js"></script>s

</head>

<body>

<p>
 $server
<h1>1000 Genomes VCF file</h1>

<div id="igvDiv" style="padding-top: 10px;padding-bottom: 10px; border:1px solid lightgray"></div>

<script type="text/javascript">

    document.addEventListener("DOMContentLoaded", function () {

        var options =
            {
                // Example of fully specifying a reference genome.  We could alternatively use  genome: "hg19"
                reference:
                    {
                        id: "hg19",
                        fastaURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/1kg_v37/human_g1k_v37_decoy.fasta",
                        cytobandURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/b37/b37_cytoband.txt"
                    },
                locus: "8:128,750,969-128,751,025",
                tracks:
                    [
                        {
                            name: "Genes",
                            type: "annotation",
                            format: "bed",
                            url: "https://s3.amazonaws.com/igv.broadinstitute.org/annotations/hg19/genes/refGene.hg19.bed.gz",
                            indexURL: "https://s3.amazonaws.com/igv.broadinstitute.org/annotations/hg19/genes/refGene.hg19.bed.gz.tbi",
                            order: Number.MAX_VALUE,
                            visibilityWindow: 300000000,
                            displayMode: "EXPANDED"
                        },
                        {
                            name: "Phase 3 WGS variants",
                            type: "variant",
                            format: "vcf",
                            url: "https://s3.amazonaws.com/1000genomes/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz",
                            indexURL: "https://s3.amazonaws.com/1000genomes/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz.tbi"
                        },
							$text_track_bam
                    ]

            };

        var igvDiv = document.getElementById("igvDiv");

        igv.createBrowser(igvDiv, options)
            .then(function (browser) {
                console.log("Created IGV browser");
            })

    })

</script>

</body>

</html>
};



