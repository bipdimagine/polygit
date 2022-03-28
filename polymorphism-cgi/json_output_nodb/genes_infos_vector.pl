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
use filter_cache qw/get_data/;  
use patient_filter;
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



my $cgi = new CGI();
my $project_name = $cgi->param('project');
my $chr_name = $cgi->param('chromosome');
my $gene_id = $cgi->param('gene');
my $var_ids = $cgi->param('ids');

my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $project_name );
my $ensemblId = $gene_id.'_'.$chr_name;
my $hAnnot = $project->liteAnnotations->get("annotations", $ensemblId);

my $hash;
$hash->{complete_name} = $hAnnot->{name}."(".$hAnnot->{external_name}.")";
$hash->{name} = $hAnnot->{name};
$hash->{position} = $chr_name.":".$hAnnot->{start}."-".$hAnnot->{end};
$hash->{description} = $hAnnot->{description};
$hash->{variations} = scalar(split(',', $var_ids));

my @list;
push (@list, $hash);
export_data::print($project,$cgi,\@list);
exit(0);

