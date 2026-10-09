#!/usr/bin/env perl
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use  File::Temp;
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use CGI qw/:standard :html3/;
use Storable qw(store retrieve freeze);
use List::MoreUtils qw(natatime);
use JSON::XS;
use GBuffer;  
use MCE::Loop;

 my $project_name;
my $final_vcf;
my $log_file;
my $patient_name;
my $fork;
my $cgi = new CGI();

my $project_name = $cgi->param("project");
my $position = $cgi->param("position");
my @patients_name = split(",",$cgi->param("patients"));
warn Dumper @patients_name;
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
print $cgi->header(

    -type    => 'text/html',

    -charset => 'utf-8'

);

my $header = qq{
	<head>
    <meta charset="utf-8">
    <title>pileup</title>

</head>

<body>

<h1>VCF variants and CRAM alignments</h1>

<div id="igvDiv" style="padding-top: 10px;padding-bottom: 10px; border:1px solid lightgray"></div>

<script type="module">

   import igv from "https://cdn.jsdelivr.net/npm/igv\@3.0.2/dist/igv.esm.min.js"

    const options =
        \{
            
            locus: "$position",
             showCytobandNames: true,
 	    reference: \{
        	id: "hg38",
        	name: "HG38",
        	fastaURL: "https://www.polyweb.fr/NGS/genome/HG38_DRAGEN/fasta/all.fa",
        	indexURL: "https://www.polyweb.fr/NGS/genome/HG38_DRAGEN/fasta/all.fa.fai"

   		 \},
            tracks:
                [
			
};
print $header."\n";
foreach my $name (@patients_name) {
	my $patient = $project->getPatient($name);
	my $bam = $patient->alignmentUrl();
	my $type = "bam";
	$type = "cram" if $bam =~ /\.cram/;
	my $type_index = "bai";
	$type_index = "crai" if $bam =~ /\.cram/;
	print qq{
		 {
                        type: 'alignment',
                        format: '${type}',
                        url: '${bam}',
                         indexURL: '${bam}.${type_index}',
                        showSoftClips: true,    
                        name: '${name}',
                        height: 200,
						showAlignments: false
          },
	};
}

print qq {
	          ]

        \}

    var igvDiv = document.getElementById("igvDiv")

    igv.createBrowser(igvDiv, options)
        .then(function (browser) {
            console.log("Created IGV browser")
        })

</script>

</body>

</html>
};
