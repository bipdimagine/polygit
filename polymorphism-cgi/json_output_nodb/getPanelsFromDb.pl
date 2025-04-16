#!/usr/bin/perl
$|=1;
use CGI qw/:standard :html3/;


use strict;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";
use GBuffer;
use Digest::MD5 qw(md5_hex);
use Spreadsheet::WriteExcel;
use JSON;

my $cgi = new CGI();
my $projectName = $cgi->param('project');
my $export_panel_name = $cgi->param('panel_name');
my $export_xls = $cgi->param('export_xls');
my $export_panels_with_gene = $cgi->param('panels_with_gene');
my $for_coverage_panel = $cgi->param('for_coverage_panel');

my $buffer = new GBuffer;
my $project = $buffer->newProject(-name => $projectName);
 

unless ($export_xls) {
	print $cgi->header('text/json-comment-filtered');
	print "{\"progress\":\".";
}

my ($hPhenotypesProject, $hPanelsPhenotypes, $hPanels, $hPhenotypes, $hPanelsNbGenes);
foreach my $phenotype (@{$project->getPhenotypes()}) {
	$hPhenotypesProject->{$phenotype->name()} = undef;
}

my $query_panels = $buffer->queryPanel();
foreach my $panel_id (@{$query_panels->getAllPanelsIds()}) {
	print "." unless ($export_xls);
	my $hInfosPanel = $query_panels->getPanelInfos($panel_id);
	next unless ($hInfosPanel->{$panel_id}->{current} == 1);
	$hPanels->{$panel_id} = $hInfosPanel->{$panel_id};
}

my $query_phenotypes = $buffer->queryPhenotype();
foreach my $phenotype_id (@{$query_phenotypes->getAllPhenotypes()}) {
	print "." unless ($export_xls);
	my $hInfosPheno = $query_phenotypes->getPhenotypeInfos($phenotype_id);
	my $phenotype_name = $hInfosPheno->{$phenotype_id}->{name};
	$hPhenotypes->{$phenotype_name}->{'All Genes'}->{id} = $phenotype_id.'_'.9999998;
	$hPhenotypes->{$phenotype_name}->{'All Genes'}->{name} = 'All Genes';
	$hPhenotypes->{$phenotype_name}->{'All Genes'}->{genes}->{'All Genes'} = 'All Genes';
	$hPhenotypes->{$phenotype_name}->{'All Genes'}->{nb_genes} = 'all';
	
	$hPhenotypes->{$phenotype_name}->{'All Panels Genes'}->{id} = $phenotype_id.'_'.9999999;
	$hPhenotypes->{$phenotype_name}->{'All Panels Genes'}->{name} = 'All Panels Genes';
	$hPhenotypes->{$phenotype_name}->{'All Panels Genes'}->{genes}->{'All Panels Genes'} = 'All Panels Genes';
	$hPhenotypes->{$phenotype_name}->{'All Panels Genes'}->{nb_genes} = 'all';
	$hPhenotypes->{$phenotype_name}->{'All Panels Genes'}->{description} = 'All Panels Genes';
	my $only_current_panels = 1;
	foreach my $panel_id (@{$query_phenotypes->getPanelsId($phenotype_id, $only_current_panels)}) {
		print "." unless ($export_panel_name);
		my $panel_name = $hPanels->{$panel_id}->{name};
		next unless ($panel_name);
		my $correct_panel_name;
		if ($export_panel_name) {
			$correct_panel_name = $panel_name;
			$correct_panel_name =~ s/-/_/g;
			$correct_panel_name =~ s/ /_/g;
			$correct_panel_name =~ s/\./_/g;
			next unless ($export_panel_name eq $correct_panel_name);
		}
		print "." unless ($export_panel_name);
		my $hGenes;
		if ($export_panel_name) {		
			$hGenes = $query_panels->getGenesForPanels($panel_id);
			$panel_name = $correct_panel_name;
		} 
		$hPhenotypes->{$phenotype_name}->{$panel_name} = $hPanels->{$panel_id};
		$hPhenotypes->{$phenotype_name}->{$panel_name}->{id} = $phenotype_id.'_'.$panel_id;
		$hPhenotypes->{$phenotype_name}->{$panel_name}->{source} = $hPanels->{$panel_id}->{source};
		$hPhenotypes->{$phenotype_name}->{$panel_name}->{description} = $hPanels->{$panel_id}->{description};
		$hPhenotypes->{$phenotype_name}->{$panel_name}->{creator} = $hPanels->{$panel_id}->{creator};
		$hPhenotypes->{$phenotype_name}->{$panel_name}->{creator_email} = $hPanels->{$panel_id}->{creator_email};
		if ($export_panel_name) {
			foreach my $ensg (keys %$hGenes) {
				my $gene_name = $hGenes->{$ensg}->{'name'};
				$hPhenotypes->{$phenotype_name}->{$panel_name}->{genes}->{$gene_name} = undef;
			}
			$hPhenotypes->{$phenotype_name}->{$panel_name}->{nb_genes} = scalar keys %{$hGenes};
		}
		else {
			$hPhenotypes->{$phenotype_name}->{$panel_name}->{nb_genes} = $query_panels->getNbGenesForPanels($panel_id);
		}
	}
}

# EXPORT XLS -> Genes in Panel
if ($export_panel_name and $export_xls) {
	foreach my $pheno_name (keys %{$hPhenotypes}) {
		foreach my $panel_name (keys %{$hPhenotypes->{$pheno_name}}) {
			next unless (lc($panel_name) eq lc($export_panel_name));
			my $panel_name_text = lc($panel_name);
			$panel_name_text =~ s/ //g;
			my $file_name = 'panel_'.$panel_name_text.'.xls';
			print "Content-type: application/msexcel\n";
			print "Content-Disposition: attachment;filename=$file_name\n\n";
			my $workbook = Spreadsheet::WriteExcel->new( \*STDOUT );
			my $xls_page = $workbook->add_worksheet('PANEL');
			$xls_page->write(0, 0, 'Panel '.$export_panel_name);
			my $i = 1;
			foreach my $gene_name (sort keys %{$hPhenotypes->{$pheno_name}->{$panel_name}->{genes}}) {
				$xls_page->write($i, 0, $gene_name);
				$i++;
			}
			exit(0);
		}
	}
	exit(0);
}

# VIEW PANELS with Gene
if ($export_panels_with_gene) {
	my (@lItems, $h);
	my $panel_list = "<font style='color:red;'><b><u>Gene ".uc($export_panels_with_gene)."</b></u></font><br><br>";
	foreach my $pheno_name (sort keys %{$hPhenotypes}) {
		my @lPanels_ok;
		foreach my $panel_name (keys %{$hPhenotypes->{$pheno_name}}) {
			if (exists $hPhenotypes->{$pheno_name}->{$panel_name}->{genes}->{uc($export_panels_with_gene)}) {
   				push(@lPanels_ok, $panel_name);
			}
		}
		if (scalar @lPanels_ok > 0) {
			$panel_list .= "<b><u>Phenotype: ".$pheno_name."</b></u><br>";
			foreach my $panel_name (sort @lPanels_ok) { $panel_list .= '  - '.$panel_name."<br>"; }
			$panel_list .= "<br>";
		}
	}
	$h->{genes_table_html} = $panel_list;
	push(@lItems, $h);
	my $hashRes;
	$hashRes->{'label'} = 'id';
	$hashRes->{'items'} = \@lItems;
	print $cgi->header('text/json-comment-filtered');
	print encode_json $hashRes;
	exit(0);
}

my $h_items;
my @lItems;
foreach my $pheno_name (sort keys %{$hPhenotypes}) {
	foreach my $panel_name (sort keys %{$hPhenotypes->{$pheno_name}}) {
		my $h;
		$h->{phenotype_project} = 0;
		if (exists $hPhenotypesProject->{$pheno_name}) { $h->{phenotype_project} = 1; }
		$h->{phenotype_name} = $pheno_name;
		$h->{name} = $panel_name;
		$h->{id} = $pheno_name.'_'.$panel_name;
		$h->{description} = $hPhenotypes->{$pheno_name}->{$panel_name}->{description};
		$h->{creator} = $hPhenotypes->{$pheno_name}->{$panel_name}->{creator};
		$h->{creator_email} = $hPhenotypes->{$pheno_name}->{$panel_name}->{creator_email};
		$h->{source} = $hPhenotypes->{$pheno_name}->{$panel_name}->{source};
		$h->{nb_genes} = $hPhenotypes->{$pheno_name}->{$panel_name}->{nb_genes};
		$h->{genes_table_html} = join(', ', sort keys %{$hPhenotypes->{$pheno_name}->{$panel_name}->{genes}});
		if ($for_coverage_panel) { $h_items->{lc($panel_name)} = $h; }
		else { push(@lItems, $h); }
	}
}
if ($for_coverage_panel) { 
	foreach my $panel_name (sort keys %$h_items) {
		push(@lItems, $h_items->{$panel_name});
	}	
}

my $hashRes;
$hashRes->{'label'} = 'id';
$hashRes->{'items'} = \@lItems;
my $json_encode = encode_json $hashRes;
print ".\",";
$json_encode =~ s/{//;
print $json_encode;
exit(0);