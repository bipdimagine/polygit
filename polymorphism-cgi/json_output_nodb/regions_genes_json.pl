#!/usr/bin/perl
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/util_filter";
#use filter_cache qw/get_genes_all/;  
use connect;
use GBuffer;
use GenBoStorable;
use Data::Dumper;
use Bio::DB::Sam;
use util_file;
use Storable qw/freeze thaw nfreeze nstore_fd nstore retrieve/;
use export_data;
use threads;
use threads::shared;
#use Set::Object; 
use Set::Intersection;



my $cgi    = new CGI();
my $buffer = new GBuffer;

my $project_name = $cgi->param('project');
my $project = $buffer->newProject( -name => $project_name );
die( "unknown project" . $project_name ) unless $project;

my $select_ids   = $cgi->param('ids');

my @ids = split(",",$select_ids);

#filter pour les patients par nb 

my $chr_name = $cgi->param('chr');

my $start = $cgi->param('start');
my $end = $cgi->param('end');


my $chr = $project->getChromosome($chr_name);

my $gg = $chr->getGenesByPosition($start,$end); 

my @data_out;
foreach my $gb (@$gg){
		my %items;
		my $g = $gb->id;
		$g =~ s/_.*//;
	
		$items{label} = $g."(".$gb->external_name.")";
		$items{name} = $gb->id;
		push(@data_out,\%items);
}
export_data::printWithLabel($project,$cgi,\@data_out);
exit(0);
		
		