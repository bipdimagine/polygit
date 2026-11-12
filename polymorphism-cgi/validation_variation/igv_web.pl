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

my $final_vcf;
my $log_file;
my $fork;
my $cgi = new CGI();

my $project_name = $cgi->param("project");
my $position = $cgi->param("position");
my @patients_name = split(",",$cgi->param("patients"));
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );

my ($genome, $genome_id, $genome_name, $fasta_url, $fasta_index_url, $gencode_url, $gencode_name);

my $locus;
$locus = qq{locus: "$position",};

if ($project->getVersion() =~  /HG19/) {
	$genome = 'hg19';
	$genome_id = "HG19_MT";
	$genome_name = "HG19_MT";
	$fasta_url = "/NGS/genome/HG19_MT/fasta/all.fa";
	$fasta_index_url = "/NGS/genome/HG19_MT/fasta/all.fa.fai";
	if ($project->gencode_version() ne '-1') {
		$gencode_name = 'Gencode-v'.$project->gencode_version();
		$gencode_url = "/NGS/repository/HG19/annotations/gencode.v".$project->gencode_version().'/gtf/';
		$gencode_url .= 'gencode.v'.$project->gencode_version().'.annotation.gtf.gz';
	}
}
elsif ($project->getVersion() =~  /HG38/) {
	$genome = 'hg38';
	$genome_id = "HG38";
	$genome_name = "HG38";
	$fasta_url = "/NGS/genome/HG38_DRAGEN/fasta/all.fa";
	$fasta_index_url = "/NGS/genome/HG38_DRAGEN/fasta/all.fa.fai";
	if ($project->gencode_version() ne '-1') {
		$gencode_name = 'Gencode-v'.$project->gencode_version();
		$gencode_url = "/NGS/repository/HG38/annotations/gencode.v".$project->gencode_version().'/gtf/';
		$gencode_url .= 'gencode.v'.$project->gencode_version().'.annotation.gtf.gz';
		
	}
}
else {
	$genome = lc($project->getVersion());
	$genome_id = uc($project->getVersion());
	$genome_name = uc($project->getVersion());
	my $fasta_file = $project->getGenomeFasta();
	$fasta_file =~ s/.+\/genome//;
	$fasta_url = "/NGS/genome/".$fasta_file;
	$fasta_index_url = "/NGS/genome/".$fasta_file.".fai";
}




my $height_patient_track = '500';
$height_patient_track = '250' if scalar(@patients_name) > 1;

my $html;
$html .= qq{
	<head>
    <meta charset="utf-8">
    <title>PolyWeb IGV</title>

</head>

<body>

<div id="igvDiv" style="padding-top: 10px;padding-bottom: 10px;"></div>

<script type="module">
	import igv from "https://cdn.jsdelivr.net/npm/igv\@3.8.10/dist/igv.esm.min.js"
    const options = \{
		$locus
		showCytobandNames: true,
};

$html .= qq{
	reference: {
       	id: "$genome_id",
       	name: "$genome_name",
       	fastaURL: "$fasta_url",
      	indexURL: "$fasta_index_url"
    },
};

$html .= qq{
		tracks: [        
};

if ($genome eq 'hg19' or $genome eq 'hg38') {
	#in error, pas le bon path
#	my $gencode_version = $project->gencode_version;
#	$gtf_file = $project->gtf_file();
#	$gtf_file =~ s/.+\/repository\///;
#	$gtf_file = '/NGS/genome/'.$gtf_file;
#	
#	$html .= qq{
#		{
#			name: 'Gencode $gencode_version',
#			format: 'refgene',
#			url: '$gtf_file',
#			order: '1000000',
#			indexed: 'false',
#		},
#	};
	$html .= qq{
		{
			name: 'Refseq All',
			format: 'refgene',
			url: 'https://hgdownload.soe.ucsc.edu/goldenPath/$genome/database/ncbiRefSeq.txt.gz',
			infoURL: 'https://www.ncbi.nlm.nih.gov/gene/?term=$$',
			order: '1000000',
			indexed: 'false',
		},
	};
	
}

my $h_pat_done;
foreach my $name (@patients_name) {
	my $patient = $project->getPatient($name);
	my $bam = $patient->alignmentUrl();
	my $type = "bam";
	$type = "cram" if $bam =~ /\.cram/;
	my $type_index = "bai";
	$type_index = "crai" if $bam =~ /\.cram/;
	
	$html .= qq{
		 {
				type: 'alignment',
				format: '${type}',
				url: '${bam}',
				indexURL: '${bam}.${type_index}',
				showSoftClips: true,    
				name: '${name}',
				height: $height_patient_track,
				showAlignments: true,
          },
	};
	
#	foreach my $vcf_file (@{$patient->getVariationsFiles()}) {
#		$vcf_file =~ s/.+ngs\///;
#		$vcf_file = "/NGS/".$vcf_file;
#		my $text_name = $patient->name().' Variants';
#		$html .= qq{
#			 {
#					format: 'vcf',
#					url: '$vcf_file',
#					name: '$text_name',
#					height: 60,
#	          },
#		};
#	}
}
	
if ($gencode_url) {
	$html .= qq{
		 {
				format: 'gtf',
				url: '${gencode_url}',
				name: '${gencode_name}',
          },
	};
}

$html .= qq{
	          ]

        \}

    var igvDiv = document.getElementById("igvDiv");

    igv.createBrowser(igvDiv, options)
        .then(function (browser) {
            console.log("Created IGV browser")
        })

</script>

</body>

</html>
};


# Autoriser l'accès depuis n'importe quel domaine (ou remplacez * par votre domaine exact)
print "Access-Control-Allow-Origin: *\n";
print "Access-Control-Allow-Methods: GET, POST, OPTIONS\n";
print "Access-Control-Allow-Headers: Content-Type\n";
print "Content-type: text/html; charset=utf-8\n\n"; # Votre entête JSON habituel
print $html;
exit(0);
