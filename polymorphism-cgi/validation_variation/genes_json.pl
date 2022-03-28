#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/kyoto";
use lib "$Bin/../packages/export";
use Carp;
use export_data;
use strict;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use JSON::XS;
use validationQuery;
my $buffer = GBuffer->new();
my $cgi    = new CGI();
my $project_name = $cgi->param('project');
#my $project_name ;
#GetOptions(
#	'project=s'   => \$project_name
#	);
#warn $project_name ;
#my $outfile = "/data-xfs/dev/cfourrag/polymorphism-cgi/cecile/essai_gene_2.xls" ;

my $project = $buffer->newProject(-name=>$project_name);
my $bundle = $cgi->param('bundle');
my $xls=$cgi->param('xls');
#my $xls =1;
my $capture = $project->getCaptures();
my	$transcripts_cgi = $project->bundle_transcripts() ;
my @out;
my $item;


if ($bundle){

	my $bundle = $project->bundles;
	my %bundle;
	
	foreach my $bb (values %$bundle){
		my $b = $bb->[0];
		my $item;
		$item->{label} = $b->{name};
		$item->{description} = $b->{description}; 
		
		$item->{transcript} = scalar(@$bb);#;
		$item->{id} = $b->{ID};
		$item->{active} = 1;
		
		push(@out,$item);	
	}
}

else{
	foreach my $tr (@$transcripts_cgi){
		my $debug;
		$debug =1  if $tr eq "ENST00000378016";
		warn $tr if $tr eq "ENST00000378016";
		
		my $tr1 = $project->newTranscript($tr);
		warn $tr1->genomic_span->as_string if $debug;
		my $item;
		$item->{label} = $tr1->getGene->external_name();
		$item->{bundle} = join(";",@{$project->return_bundle($tr)});
		$item->{bundles}= $project->return_bundle($tr);#;
		$item->{description} = $tr1->getGene->description(); 
		$item->{description} =~ s/\[.*\]//;
		warn "coucou" if $debug;
		$item->{nb_exons} = "".scalar(@{$tr1->getExons()})."";

		warn $item->{nb_exons} if $debug;
		$item->{transcript} = $tr;
		$item->{id} = $tr;
		$item->{ccds} = $tr1->{ccds_name};
		$item->{ccds} =~ s/\..*//;
		$item->{refseq} = $tr1->{external_name};
		$item->{refseq} =~ s/\..*//;
		$item->{protein} = $tr1->{external_protein_name};
		$item->{start} = $tr1->{start};
		$item->{end} = $tr1->{end};
		$item->{chr} = $tr1->getChromosome()->name();
		$item->{coucou}="coucou";
		push(@out,$item);		
	}
}

my @out2 = sort{$a->{label} cmp $b->{label}} @out;
unless ($xls){
		export_data::print($project,$cgi,\@out2);
		}
else {
		printTableGenesXls() ;
		}

#if ($xls) {printTableGenesXls()};
#printTableGenesXls();

exit(0);


sub printTableGenesXls {
	print "Content-type: application/msexcel\n";
	print "Content-Disposition: attachment;filename=".$project->name()."-export_genes".".xls\n\n";
	
	my $workbook  = Spreadsheet::WriteExcel->new(\*STDOUT);
#	my $workbook = Spreadsheet::WriteExcel->new($outfile) ;
	my $worksheet = $workbook->add_worksheet("ENST");

	my $format_header = $workbook->add_format(border => 1, underline => 1);
	$format_header->set_bold();
	$format_header->set_align('center');
	$format_header->set_fg_color('silver');
	my $format_line = $workbook->add_format();
	$format_line->set_align('left');
	$format_line->set_color('black');
	my $format_line_bis = $workbook->add_format();
	$format_line_bis->set_align('center');
	$format_line_bis->set_color('black');
	$format_line_bis->set_fg_color('white');
	my $format_line_2 = $workbook->add_format();
	$format_line_2->set_align('center');
	$format_line_2->set_color('black');
	$format_line_2->set_fg_color('silver');
	my $format_line_3 = $workbook->add_format();
	$format_line_3->set_align('left');
	$format_line_3->set_color('black');
	$format_line_3->set_fg_color('silver');
	
	my @infos = ("name","transcript","chr","start","end","refseq","ccds","nb_exon","group","description");
	my $col = 0; 
	my $row =0;
	
	for (my $i=0;$i<@infos;$i++){
		$worksheet->write($row,$i,$infos[$i],$format_header);
	}
		
	$row++;
	foreach my $item_gene (@out2) {
		$worksheet->write($row,0,$item_gene->{label},$format_line_3);
		$worksheet->write($row,1,$item_gene->{transcript},$format_line_bis);
		$worksheet->write($row,2,$item_gene->{chr},$format_line_2);
		$worksheet->write($row,3,$item_gene->{start},$format_line);
		$worksheet->write($row,4,$item_gene->{end},$format_line_3);
		$worksheet->write($row,5,$item_gene->{refseq},$format_line);
		$worksheet->write($row,6,$item_gene->{ccds},$format_line_3);
		$worksheet->write($row,7,$item_gene->{nb_exons},$format_line_bis);
		$worksheet->write($row,8,$item_gene->{bundle},$format_line_3);
		$worksheet->write($row,9,$item_gene->{description},$format_line);
		$row++;
		}
	

		$workbook->close();
		exit(0) ;
}

