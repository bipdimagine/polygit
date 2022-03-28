#!/usr/bin/perl
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use GBuffer;
use GenBoStorable;
use Data::Dumper;
use export_excel;  
use export_data;
use Sys::Hostname;
use get_variations;
use Storable qw/freeze thaw nfreeze nstore_fd nstore retrieve/;
my %types = ( variations => "variations" );

my $cgi    = new CGI();
my $buffer = new GBuffer;



my $host = hostname;



my $project_name = $cgi->param('project');

my $type_name = $cgi->param('type');
my $type2    = $cgi->param('type2');
my $select_ids   = [split(",",$cgi->param('ids'))];
#creation du projet à partir du buffer
my $project = $buffer->newProject( -name => $project_name );
die( "unknown project" . $project_name ) unless $project;
#$project->getMethods($buffer->getType("variation"));
if ($type_name eq "variations" && $project->is_cnv){
	$type_name = "cnvs";
}

my $temp_type_name = $type_name;

my $chr = $project->getChromosome( $cgi->param('chromosome'));
die() unless $chr;
my $data = get_variations::getIds($buffer,$project,$chr,$temp_type_name,$select_ids);


my $genes;
foreach my $var (@$data){
	map {$genes->{$_}++} @{$var->{genes}};
	
	#warn join("\n",keys %$var);
}


my @data_out;
foreach my $gene_name (keys %$genes){
	my $gene = $project->newGene($gene_name);
	my $transcripts = $gene->getTranscripts();

	foreach my $tr (@$transcripts) {
		my %items;
		$items{name} = $tr->id();
		
		my $xref;
#		my (@ext_reftr) = grep {$_->dbname() =~/refseq/i} @{$tr->get_all_DBEntries()};
		$items{label} =$tr->id();
		$items{label} = $xref . " (". $tr->external_name . ") " if $tr->external_name ;
		
		push(@data_out,\%items);
		 
	}

}

export_data::printWithLabel($project,$cgi,\@data_out);
exit(0);
