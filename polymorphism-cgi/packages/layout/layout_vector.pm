package layout_vector;
use strict;
use Data::Dumper;


my $layout;
$layout->{variations}->{ngs} = [
#	{
#		name   => "Row",
#		get    => "getRow",
#		width  => "3%",
#		styles => 'text-align:center;'
#	},
	{
		name      => "Name",
		field     => "name",
		formatter => "formatSnp",
		width     => "10%",
		styles    => 'text-align:center;'
	},
	{
		name   => "Position",
		field  => "position",
		width  => "10%",
		styles => 'text-align:center;'
	},
	{
		name   => "Allele",
		field  => "text",
		width  => "7%",
		styles => 'text-align:center;'
	},
	{
		name      => "Type",
		field     => "structural_type",
		width     => "5%",
		styles    => 'text-align:center;',
		hidden => 'true',
	},
	{
		name      => "DejaVu",
		field     => "bipd_db",
		formatter => "formatDejavu",
		width     => "10%",
		styles    => 'text-align:center;'
	},
	{
		name      => "Ho/He",
		field     => "homo_hetero",
		formatter => "formatHomozygote",
		width     => "5%",
		styles    => 'text-align:center;'
	},
	{
		name   => "Chr",
		field  => "chromosome",
		width  => "5%",
		styles => 'text-align:center;',
		hidden => 'true',
	},
	{
		name   => "Position",
		field  => "start",
		width  => "10%",
		styles => 'text-align:center;',
		hidden => 'true',
	},
	{
		name   => "Patients",
		field  => "nb_patients",
		width  => "5%",
		styles => 'text-align:center;'
	},

	{
		name   => "Consequence",
		field  => "consequence!all",
		formatter => "formatConsequence",
		width  => "10%",
		styles => 'text-align:center;'
	},
	{
		name      => "Gene ",
		field     => "external_names",
		formatter => "formatReferenceGene",
		width  	  => "10%",
		styles    => 'text-align:center;'
	},
#	{
#		name      => "Poly<br>Score",
#		field     => "polyweb_score",
#		formatter => "formatPolywebScoreDetailProject",
#		width     => "5%",
#		styles    => 'text-align:center;'
#	},
	{
		name      => "Polyphen<br>status",
		field     => "polyphen_status",
		formatter => "formatPolyphen",
		width     => "5%",
		styles    => 'text-align:center;'
	},
	{
		name      => "Sift<br>status",
		field     => "sift_status",
		formatter => "formatSift",
		width     => "5%",
		styles    => 'text-align:center;'
	},
#	{
#		name      => "PolyWeb<br>status",
#		field     => "polyweb_score",
#		formatter => "formatPolywebScoreDetailProject",
#		width     => "5%",
#		styles    => 'text-align:center;'
#	},
	{
		name      => "SpliceAI",
		field     => "spliceAI",
		formatter => "formatSpliceAI",
		width     => "5%",
		styles    => 'text-align:center;'
	},
	{
		name      => "Cadd<br>score",
		field     => "cadd_score",
		width     => "5%",
		styles    => 'text-align:center;'
	},
	{
		name      => "Ncboost",
		field     => "ncboost",
		formatter => "formatNcboost",
		width     => "5%",
		styles    => 'text-align:center;'
	},
	{
		name      => "Varsome",
		field     => "varsome",
		formatter => "formatVarsome",
		width     => "5%",
		styles    => 'text-align:center;'
	},
	{
		name      => "gAD AC",
		field     => "gnomad_ac",
		formatter => "formatGnomadLink",
		width     => "5%",
		styles    => 'text-align:center;'
	},
	{
		name      => "gAD HO",
		field     => "gnomad_ho",
		formatter => "formatGnomadLink",
		width     => "5%",
		styles    => 'text-align:center;'
	},
	{
		name      => "gAD AN",
		field     => "gnomad_an",
		formatter => "formatGnomadLink",
		width     => "5%",
		styles    => 'text-align:center;'
	},
	{
		name      => "DB Freq",
		field     => "freq",
		width     => "5%",
		styles    => 'text-align:center;'
	},
#	{
#		name      => "Is<br>Clinical",
#		field     => "isClinical",
#		width     => "5%",
#		styles    => 'text-align:center;'
#	},
	{
		name		=> "Gene<br>Omim",
		field		=> "gene_in_omim",
		width		=> "4%",
		styles		=> 'text-align:center;font-size : 09px;',
		formatter	=> "format_is_omim",
	}, 
	{
		name      => "HGMD<br>Class",
		field     => "hgmd_class",
		formatter => "formatHgmd",
		width     => "5%",
		styles    => 'text-align:center;'
	},
];
$layout->{variations}->{cnv} = [
	{
		name   => "Row",
		get    => "getRow",
		width  => 2,
		styles => 'text-align:center;'
	},
	{
		name      => "Name",
		field     => "name",
		formatter => "formatCnv",
		width     => 12,
		styles    => 'text-align:center;'
	},
	{
		name      => "DGV",
		field     => "linkdgv",
		formatter => "formatDgv",
		width     => 2,
		styles    => 'text-align:center;'
	},
	{
		name   => "Chr",
		field  => "chromosome",
		width  => 4,
		styles => 'text-align:center;'
	},
	{
		name   => "start",
		field  => "start",
		width  => 8,
		styles => 'text-align:center;'
	},
	{
		name   => "end",
		field  => "end",
		width  => 8,
		styles => 'text-align:center;'
	},
	{
		name   => "Patients",
		field  => "nb_patients",
		width  => 3,
		styles => 'text-align:center;'
	},
	{
		name      => "method",
		field     => "method",
		formatter => "formatMethods",
		width     => 3,
		styles    => 'text-align:center;'
	},
	{
		name   => "Consequence",
		field  => "consequence!all",
		styles => 'text-align:center;'
	},
	{
		name      => "Gene",
		field     => "genes",
		formatter => "formatReferenceGene",
		width     => 10,
		styles    => 'text-align:center;'
	},
	{
		name      => "Valid",
		field     => "valid",
		formatter => "formatValid",
		width     => 3,
		styles    => 'text-align:center;'
	}

];
$layout->{variations}->{array} = [
	{
		name   => "Row",
		get    => "getRow",
		width  => 2,
		styles => 'text-align:center;'
	},
	{
		name   => "Name",
		field  => "name",
		width  => 12,
		styles => 'text-align:center;'
	},
	{
		name      => "type",
		field     => "structural_type",
		width     => 3,
		styles    => 'text-align:center;'
	},
	{
		name   => "Nomenclature",
		field  => "nomenclature",
		width  => 10,
		styles => 'text-align:center;'
	},
	{
		name   => "Allele",
		field  => "text",
		width  => 3,
		styles => 'text-align:center;'
	},
	{
		name   => "Position",
		field  => "start",
		width  => 8,
		styles => 'text-align:center;'
	},
	{
		name   => "Nb traces",
		field  => "lectures",
		width  => 3,
		styles => 'text-align:center;'
	},
	{
		name      => "method",
		field     => "method",
		formatter => "formatMethods",
		width     => 3,
		styles    => 'text-align:center;'
	},
	{
		name   => "Consequence",
		field  => "consequence",
		styles => 'text-align:center;'
	},
	{
		name      => "Gene",
		field     => "genes",
		formatter => "formatReferenceGene",
		width     => 10,
		styles    => 'text-align:center;'
	},
	{
		name      => "Filter",
		field     => "filter",
		formatter => "formatFilter",
		width     => 5,
		styles    => 'text-align:center;'
	},
	{
		name      => "Valid",
		field     => "valid",
		formatter => "formatValid",
		width     => 3,
		styles    => 'text-align:center;'
	},
	{
		name      => "polyphen<br>status",
		field     => "polyphen_status",
		formatter => "formatPolyphen",
		width     => 5,
		styles    => 'text-align:center;'
	}

];

$layout->{variations}->{classic} = [
	{
		name   => "Row",
		get    => "getRow",
		width  => 2,
		styles => "text-align:center;"
	},
	{
		name      => "Name",
		field     => "name",
		formatter => "formatSnp",
		width     => 12,
		styles    => 'text-align:center;'
	},
	{
		name      => "type",
		field     => "structural_type",
		width     => 3,
		styles    => 'text-align:center;'
	},
	{
		name   => "Nomenclature",
		field  => "nomenclature",
		width  => 10,
		styles => 'text-align:center;'
	},
	{
		name   => "Allele",
		field  => "text",
		width  => 3,
		styles => 'text-align:center;'
	},
	{
		name   => "Position",
		field  => "start",
		width  => 8,
		styles => 'text-align:center;'
	},
	{
		name   => "contig",
		field  => "group",
		width  => 5,
		styles => 'text-align:center;'
	},
	{
		name   => "Nb traces",
		field  => "lectures",
		width  => 3,
		styles => 'text-align:center;'
	},
	{
		name      => "method",
		field     => "method",
		formatter => "formatMethods",
		width     => 3,
		styles    => 'text-align:center;'
	},
	{
		name   => "Consequence",
		field  => "consequence!all",
		styles => 'text-align:center;'
	},
	{
		name      => "Gene",
		field     => "genes",
		formatter => "formatReferenceGene",
		width     => 10,
		styles    => 'text-align:center;'
	},
	{
		name      => "Filter",
		field     => "filter",
		formatter => "formatFilter",
		width     => 5,
		styles    => 'text-align:center;'
	},
	{
		name      => "Valid",
		field     => "valid",
		formatter => "formatValid",
		width     => 3,
		styles    => 'text-align:center;'
	},
	{
		name      => "polyphen<br>status",
		field     => "polyphen_status",
		formatter => "formatPolyphen",
		width     => 5,
		styles    => "text-align:center;"
	}

];
########"
#### indels
#############
$layout->{indels}->{classic} = [
	{
		name   => "Row",
		get    => "getRow",
		width  => 2,
		styles => 'text-align:center;'
	},
	{
		name   => "Name",
		field  => "name",
		width  => 15,
		styles => 'text-align:center;'
	},
	{
		name   => "Nomenclature",
		field  => "nomenclature",
		width  => 15,
		styles => 'text-align:center;'
	},
	{
		name   => "sequence",
		field  => "text",
		width  => 3,
		styles => 'text-align:center;'
	},
	{
		name   => "Position",
		field  => "start",
		width  => 6,
		styles => 'text-align:center;'
	},
	{
		name   => "Nb traces",
		field  => "lectures",
		width  => 3,
		styles => 'text-align:center;'
	},
	{
		name   => "Bipd",
		field  => "indel_bipd",
		formatter => "formatMethods",
		width  => 5,
		styles => 'text-align:center;'
	},
	{
		name   => "Polyphred",
		field  => "indel_polyphred",
		formatter => "formatMethods",
		width  => 5,
		styles => 'text-align:center;'
	},
	{
		name   => "Consequence",
		field  => "consequence!all",
		width  => 7,
		styles => 'text-align:center;'
	},
	{
		name      => "Valid",
		field     => "valid",
		formatter => "formatValid",
		width     => 3,
		styles    => 'text-align:center;'
	}
];

$layout->{indels}->{ngs} = [
	{
		name   => "Row",
		get    => "getRow",
		width  => 2,
		styles => 'text-align:center;'
	},
	{
		name   => "Name",
		field  => "name",
		width  => 15,
		styles => 'text-align:center;'
	},
	{
		name   => "Nomenclature",
		field  => "nomenclature",
		width  => 15,
		styles => 'text-align:center;'
	},
	{
		name   => "sequence",
		field  => "text",
		width  => 3,
		styles => 'text-align:center;'
	},
	{
		name      => "Ho/He",
		field     => "homo_hetero",
		formatter => "formatHomozygote",
		width     => 4,
		styles    => 'text-align:center;'
	},
	{
		name   => "Position",
		field  => "start",
		width  => 6,
		styles => 'text-align:center;'
	},
	{
		name   => "Patients",
		field  => "lectures",
		width  => 3,
		styles => 'text-align:center;'
	},
	{
		name   => "maq",
		field  => "maq",
		formatter => "formatMethods",
		width  => 5,
		styles => 'text-align:center;'
	},
	{
		name      => "Gene",
		field     => "genes",
		formatter => "formatReferenceGene",
		width     => 10,
		styles    => 'text-align:center;'
	},
	{
		name   => "Consequence",
		field  => "consequence!all",
		width  => 7,
		styles => 'text-align:center;'
	},
	{
		name      => "Valid",
		field     => "valid",
		formatter => "formatValid",
		width     => 3,
		styles    => 'text-align:center;'
	}
];

$layout->{indels}->{array}   = $layout->{indels}->{classic};

########"
#### cover
#############
$layout->{cover}->{ngs} = [
	{
		name   => "Row",
		get    => "getRow",
		width  => 2,
		styles => 'text-align:center;'
	},
	{
		name  => "Name",
		field => "name",

		width  => 5,
		styles => 'text-align:center;'
	},
	{
		name  => "% cover",
		field => "maq",

		width  => 5,
		styles => 'text-align:center;'
	},
	{
		name  => "total",
		field => "score2",

		width  => 5,
		styles => 'text-align:center;'
	},
	{
		name   => "forward",
		field  => "score3",
		width  => 5,
		styles => 'text-align:center;'
	},
	{
		name   => "reverse",
		field  => "score4",
		width  => 5,
		styles => 'text-align:center;'
	}
];

$layout->{cover}->{classic} = [
	{
		name   => "Row",
		get    => "getRow",
		width  => 2,
		styles => 'text-align:center;'
	},
	{
		name  => "Name",
		field => "name",

		width  => 5,
		styles => 'text-align:center;'
	},
	{
		name  => "base",
		field => "base",

		width  => 5,
		styles => 'text-align:center;'
	},
	{
		name   => "Bipd",
		field  => "indel_bipd",
		formatter => "formatMethods",
		width  => 5,
		styles => 'text-align:center;'
	},
	{
		name   => "Polyphred",
		field  => "indel_polyphred",
		formatter => "formatMethods",
		width  => 5,
		styles => 'text-align:center;'
	},
	{
		name      => "valid",
		field     => "valid",
		formatter => "formatValid",
		width     => 2,
		styles    => 'text-align:center;'
	}

];

$layout->{cover}->{array} = [
	{
		name   => "Row",
		get    => "getRow",
		width  => 2,
		styles => 'text-align:center;'
	},
	{
		name  => "Name",
		field => "name",

		width  => 5,
		styles => 'text-align:center;'
	},
	{
		name  => "base",
		field => "base",

		width  => 5,
		styles => 'text-align:center;'
	},
	{
		name  => "score",
		field => "score",

		width  => 5,
		styles => 'text-align:center;'
	},
	{
		name      => "valid",
		field     => "valid",
		formatter => "formatValid",
		width     => 2,
		styles    => 'text-align:center;'
	}

];
########"
#### cover
#############
$layout->{annex}->{ngs_pedigree} = [
#	{
#		name   => "Row",
#		get    => "getRow",
#		width  => 2,
#		styles => 'text-align:center;'
#	},
	{
		name  => "Priority",
		field => "priority",
		width  => 4,
		styles => 'font-size:9px;text-align:center;',
		hidden => 'true'
	},
	{
		name  => "Fam",
		field => "fam",
		width  => 4,
		formatter => "formatColorCell",
		styles => 'font-size:9px;text-align:center;'
	},
	{
		name  => "Name",
		field => "name",
		width  => 7,
		formatter => "formatColorCell",
		styles    => 'font-size:8px;text-align:center;'
	},
	{
		name  => "Ped",
		field => "type",
		width  => 3,
		formatter => "formatChild",
		styles => 'font-size:9px;text-align:center;'
	},
	{
		name  => "Base",
		field => "base",
		width  => 3,
		formatter => "formatColorCell",
		styles => 'font-size:10px;text-align:center;'
	},
	
	{
		name  => "He Ho",
		field => "hohe",
		width  => 2,
		formatter => "formatColorCell",
		styles => 'font-size:8px;text-align:center;'
	},
	{
		name  => "Ref/Alt",
		field => "perc",
		width  => 4,
		formatter => "formatColorCell",
		styles => 'font-size:8px;text-align:center;'
	},
#	{
#		name  => "ref",
#		field => "ref",
#		width  => 2,
#		styles => 'font-size:8px;text-align:center;'
#	},
#	
#	{
#		name  => "alt",
#		field => "alt",
#		width  => 2,
#		styles => 'font-size:8px;text-align:center;'
#	},
		{
		name  => "met",
		field => "met",
		width  => 2,
		styles => 'font-size:8px;text-align:center;'
	},
	{
		name      => "Transmission",
		field     => "transmission",
		width     => 7,
		formatter => "formatColorCell",
		#formatter => "formatTransmission",
		styles    => 'font-size:9px;text-align:center;'
	},
#	{
#		name      => "Score",
#		field     => "scale_score",
#		width     => 7,
#		formatter => "formatColorCell",
#		styles    => 'font-size:9px;text-align:center;'
#	},
	{
		name      => "In Ho region",
		field     => "in_region_ho",
		width     => 4,
		formatter => "formatRegionHo2",
		styles    => 'font-size:9px;text-align:center;'
	}

];    
$layout->{annex}->{ngs} = [
{
		name  => "Fam",
		field => "fam",
		width  => 6,
		formatter => "formatColorCell",
		styles => 'font-size:9px;text-align:left;'
	},
	{
		name  => "Name",
		field => "name",
		width  => 8,
		formatter => "formatColorCell",
		styles    => 'font-size:8px;text-align:left;'
	},
	{
		name  => "Ped",
		field => "type",
		width  => 2,
		formatter => "formatChild",
		styles => 'font-size:9px;text-align:left;'
	},
	{
		name  => "Base",
		field => "base",
		width  => 3,
		formatter => "formatColorCell",
		styles => 'font-size:9px;text-align:center;'
	},
	{
		name  => "St",
		field => "status",
		width  => 2,
		formatter => "formatColorCell",
		styles => 'font-size:9px;text-align:center;display:none;'
	},
	{
		name  => "Ho/He",
		field => "hohe",
		width  => 2,
		formatter => "formatColorCell",
		styles => 'font-size:8px;text-align:center;'
	},
	{
		name  => "Ref/Alt",
		field => "perc",
		width  => 8,
		formatter => "formatColorCell",
		styles => 'font-size:8px;text-align:center;'
	},
#	{
#		name  => "met",
#		field => "met",
#		width  => 2,
#		styles => 'font-size:8px;text-align:center;'
#	},
	{
		name      => "Transmission",
		field     => "transmission",
		width     => 7,
		formatter => "formatColorCell",
		#formatter => "formatTransmission",
		styles    => 'font-size:9px;text-align:center;'
	},
#	{
#		name      => "Score",
#		field     => "scale_score",
#		width     => 7,
#		formatter => "formatColorCell",
#		styles    => 'font-size:9px;text-align:center;'
#	},
#	{
#		name      => "valid",
#		field     => "valid",
#		formatter => "formatValid",
#		width     => 3,
#		styles    => 'font-size:7px;text-align:left;'
#	},
	{
		name      => "In Ho region",
		field     => "in_region_ho",
		width     => 4,
		formatter => "formatRegionHo2",
		styles    => 'font-size:9px;text-align:center;'
	}

]; 
$layout->{annex}->{classic} = [
	{
		name   => "Row",
		get    => "getRow",
		width  => 2,
		styles => 'text-align:center;'
	},
	{
		name  => "Name",
		field => "name",

		width  => 5,
		styles => 'text-align:center;'
	},
	{
		name  => "base",
		field => "base",

		width  => 5,
		styles => 'text-align:center;'
	},
	{
		name      => "method",
		
	},

	{
		name      => "valid",
		field     => "valid",
		formatter => "formatValid",
		width     => 2,
		styles    => 'text-align:center;'
	}

];

$layout->{annex}->{cnv} = [
	{
		name   => "Row",
		get    => "getRow",
		width  => 2,
		styles => 'text-align:center;'
	},
	{
		name  => "Name",
		field => "name",

		width  => 5,
		styles => 'text-align:center;'
	},
	{
		name      => "method",
		
	},
	{
		name      => "score",
		field     => "score2",
		width     => 5,
		styles    => 'text-align:center;'
	},
	
	{
		name      => "valid",
		field     => "valid",
		formatter => "formatValid",
		width     => 2,
		styles    => 'text-align:center;'
	}

];
$layout->{annex}->{array} =$layout->{annex}->{ngs};

##############
# consequence
##############

$layout->{consequence}->{ngs} = [
	{
		field  => "gene",
		name   => "Gene",
		width  => 15,
		styles => 'text-align:center;'
	},
	{
		field  => "description",
		name   => "Description",
		width  => 35,
		styles => 'text-align:center;'
	},
	{
		field     => "transcript",
		name      => "Transcript",
		formatter => "formatTranscriptName",
		width     => 17,
		styles    => 'text-align:center;'
	},
	{
		field     => "protein",
		name      => "Protein",
		formatter => "formatProteinName",
		width     => 10,
		styles    => 'text-align:center;'
	},
	{
		field  => "exon",
		name   => "exon",
		width  => 4,
		styles => 'text-align:center;'
	},
	{
		field  => "cdna_position",
		name   => "cdna pos",
		width  => 3,
		styles => 'text-align:center;'
	},
	{
		field  => "cds_position",
		name   => "cds pos",
		width  => 3,
		styles => 'text-align:center;'
	},
	
	
	{
		field  => "protein_position",
		name   => "prot pos",
		width  => 3,
		styles => 'text-align:center;'
	},
	{
		field  => "nomenclature",
		name   => "nomenclature",
		width  => 7,
		styles => 'text-align:center;'
	},
	{
		field  => "consequence",
		name   => "Consequence",
		width  => 7,
		styles => 'text-align:center;'
	},
	
	{
		field     => "polyphen_status",
		name      => "polyphen<br>status",
		formatter => "formatPolyphen",
		width     => 5,
		styles    => 'text-align:center;'
	},
	{
		field     => "sift_status",
		name      => "sift<br>status",
		formatter => "formatSift",
		width     => 5,
		styles    => 'text-align:center;'
	},
	
];

$layout->{consequence}->{cnv} = $layout->{consequence}->{array} = $layout->{consequence}->{classic} =
  $layout->{consequence}->{ngs};


$layout->{genes}->{ngs} =  [
				{
					name=>"Row",
					get=>"getRow",
					width=>4,
					styles=>'text-align:center;',
					rowSpan=>2
				}, {
					name=>"Name",
					field=>"name",
					#formatter=>"formatReferenceGene",
					width=>10,
					styles=>'text-align:center;',
					rowSpan=>2
				}, {
					name=>"xref",
					field=>"xref",
					width=>7,
					styles=>'text-align:center;',
					rowSpan=>2
				},{
					name=>"chr",
					field=>"chromosome" ,
					width=>2,
					styles=>'text-align:center;',
					rowSpan=>2
				}, {
					name=>"Start",
					field=>"start",
					width=>7,
					styles=>'text-align:center;',
					rowSpan=>2
				}, {
					name=>"End",
					field=>"end",
					width=>7,
					styles=>'text-align:center;',
					rowSpan=>2
				}, {
					name=>"Description",
					field=>"description",
					width=>50,
					styles=>'text-align:center;',
					rowSpan=>2
				},
				{
					name=>"All",
					field=>"v_all",
					width=>3,
					styles=>'text-align:center;',
					rowSpan=>1
				},
				 {
					name=>"Pat.",
					field=>"p_all",
					width=>3,
					styles=>'text-align:center;',
					rowSpan=>1
					},	
				
				{
					name=>"sub",
					field=>"v_substitution",
					width=>3,
					styles=>'text-align:center;',
					rowSpan=>1
				},
				 {
					name=>"Pat.",
					field=>"p_substitution",
					width=>3,
					styles=>'text-align:center;',
					rowSpan=>1
					},	
					
				
				{
					name=>"ins",
					field=>"v_insertion",
					width=>3,
					styles=>'text-align:center;',
					rowSpan=>1
				},
				{
					name=>"Pat.",
					field=>"p_insertion",
					width=>3,
					styles=>'text-align:center;',
					rowSpan=>1
					},	
				{
					name=>"del",
					field=>"v_deletion",
					width=>3,
					styles=>'text-align:center;',
					rowSpan=>1
				},
				{
					name=>"Pat.",
					field=>"p_deletion",
					width=>3,
					styles=>'text-align:center;',
					rowSpan=>1
					},				
				 {
					name=>"Coding",
					field=>"v_coding",
					width=>3,
					styles=>'text-align:center;',
					rowSpan=>1
				}, ,
					{
					name=>"Pat.",
					field=>"p_coding",
					width=>3,
					styles=>'text-align:center;',
					colSpan=>1
				},
				
				{
					name=>"Silent",
					field=>"v_silent",
					width=>5,
					styles=>'text-align:center;',
					rowSpan=>1
				},
				{
					name=>"Pat.",
					field=>"p_silent",
					width=>3,
					styles=>'text-align:center;',
					colSpan=>1
				},
				{
					name=>"UTR",
					field=>"v_utr",
					width=>3,
					styles=>'text-align:center;',
					rowSpan=>1
				},
				{
					name=>"Pat.",
					field=>"p_utr",
					width=>3,
					styles=>'text-align:center;',
					colSpan=>1
				},
				{
					name=>"Splicing",
					field=>"v_splicing",
					width=>3,
					styles=>'text-align:center;',
					rowSpan=>1
				},
				
				
				
				{
					name=>"Pat.",
					field=>"p_splicing",
					width=>3,
					styles=>'text-align:center;',
					colSpan=>1
				},	
				,
				{
					name=>"Stop",
					field=>"v_stop",
					width=>3,
					styles=>'text-align:center;',
					rowSpan=>1
				},	
					{
					name=>"Pat.",
					field=>"p_stop",
					width=>3,
					styles=>'text-align:center;',
					colSpan=>1
				},	
				{
					name=>"start_stop_lost",
					field=>"v_phase",
					width=>3,
					styles=>'text-align:center;',
					rowSpan=>1
				},	
					{
					name=>"Pat.",
					field=>"p_phase",
					width=>3,
					styles=>'text-align:center;',
					colSpan=>1
				},	
					{
					name=>"Synonymous",
					field=>"v_silent",
					width=>3,
					styles=>'text-align:center;',
					rowSpan=>1
				},	
					{
					name=>"Pat.",
					field=>"p_silent",
					width=>3,
					styles=>'text-align:center;',
					colSpan=>1
				},	
					{
					name=>"Frameshift",
					field=>"v_frameshift",
					width=>3,
					styles=>'text-align:center;',
					rowSpan=>1
				},	
					{
					name=>"Pat.",
					field=>"p_frameshift",
					width=>3,
					styles=>'text-align:center;',
					colSpan=>1
				},	
				
			];
			
#,"utr","splicing","stop","silent","phase","frameshift");

$layout->{patient}->{ngs} = [
				{ name=> "N", get=> "getRow", width=>2, styles=> 'text-align:center;'},		
				{ name=> "name",field=> "name",  width=> 4, styles=> 'text-align:center;'},		
				{ name=> "Sub",field=> "substitutions", width=> 3, styles=> 'text-align:center;'},
				{ name=> "Del",field=> "deletions", width=> 3, styles=> 'text-align:center;'},
				{ name=> "Ins",field=> "insertions", width=> 3, styles=> 'text-align:center;'},
				{ name=> "ho",field=> "homozygote", width=> 3, styles=> 'text-align:center;'},
				{ name=> "he",field=> "heterozygote", width=> 3, styles=> 'text-align:center;'},		
				{ name=> "Genes",field=> "genes", width=> 5, styles=> 'text-align:center;'},	
				{ name=> "Comp",field=> "composite", width=> 3, styles=> 'text-align:center;'},
				{ name=> "Cov",field=> "coverage", width=> 3, styles=> 'text-align:center;'},
				{ name=> "5x",field=> "5x", width=> 3, styles=> 'text-align:center;'},
				{ name=> "15x",field=> "15x", width=> 3, styles=> 'text-align:center;'},
						
				{ field=> "filter_heho",name=> 'he',formatter=>"formathe", width=> 3, styles=> 'text-align: center;'},
				{ field=> "filter_heho",name=> 'ho',formatter=>"formatho", width=> 3, styles=> 'text-align: center;'},
				{ field=> 'include', name=> 'Sl.',formatter=>"formatInclude", width=> 3, styles=> 'text-align: center;'},	
				{ field=> "fam",name=> 'fam', width=> 3, styles=> 'text-align: center;'},
				{ field=> "child",name=> 'Ped', width=> 3, formatter=>"formatChild",styles=> 'text-align: center;'},
				{ field=> "status",name=> 'st', width=> 3,formatter=>"formatDisease", styles=> 'text-align: center;'},	
			
				
						 			
		];
			
			
$layout->{genes}->{cnv}=$layout->{genes}->{classic} = $layout->{genes}->{array} = $layout->{genes}->{ngs};			

sub returnLayout {
	my ( $type, $project ) = @_;
	#warn  $project->projectType();
	my $projectType = $project->projectType();
	$projectType = "ngs_pedigree" if ($project->isFamilial() && $type eq "annex");
	my $obj_name = "variation";
	$obj_name = "cnv" if $project->is_cnv;
	if ($type eq "variations" || $type eq "annex" ) {
	
		my $methods = ["calling"];#$project->getCallingMethods();
		
		return addMethods2Layout($methods,$layout->{$type}->{$projectType});
	}	

	
	return $layout->{$type}->{$projectType};
}
sub returnLayout_nodb {
	my ( $type, $project ) = @_;
	
	my $projectType = $project->projectType();
	
	$projectType = "ngs_pedigree" if ($project->isFamilial() && $type eq "annex");
	my $obj_name = "variation";
	$obj_name = "cnv" if $project->is_cnv;
	if ($type eq "variations" || $type eq "annex" ) {
	
		#my $methods = $project->getCallingMethods();
		my $methods = ["calling"];
		return addMethods2Layout($methods,$layout->{$type}->{$projectType});
	}	
	

	return $layout->{$type}->{$projectType};
}
sub addMethods2Layout {
	my ($methods,$layout) = @_;
	
	my @final_layout;
	foreach my $col (@$layout){
		if  ($col->{name} eq "method") {
			
			foreach my $method (@$methods){
		
				my $m = {
					name      => $method,
					field     => $method,
					formatter => "formatMethods",
					width     => 4,
					styles    => 'text-align:left;font-size:9px;'
				};
				push(@final_layout,$m);# unless $col->name ne "method";
				
			}
			next;
		}
		push(@final_layout,$col);# unless $col->name ne "method";
		
	} 
	
	return \@final_layout;
}
1;
