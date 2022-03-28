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



my $cgi    = new CGI();
my $buffer = new GBuffer;

my $project_name = $cgi->param('project');
my $project = $buffer->newProject( -name => $project_name );
die( "unknown project" . $project_name ) unless $project;


my $fm;

my ($data,$kyoto) = get_data($project,());

my $mode = $cgi->param('mode');

my $fm =  new patient_filter(patients=>$data->{patients},data=>$data->{chromosomes},project=>$project,type=>1);



my $select_ids   = $cgi->param('ids');

my @ids = split(",",$select_ids);


#filter pour les patients par nb 

my $chromosome = $cgi->param('chromosome');
$fm->chr($chromosome);
my @out;

foreach my $gid (keys  %{$data->{chromosomes}->{$chromosome}->{genes}}) { 
	$fm->gid($gid);
	my $find = undef;
	my $item;
	foreach my $v (@ids){
		if (exists $fm->hashvar->{$v}){
			$find ++;
			#last;
		}
	}
	if ($find == scalar(@ids)){
		my $hgene = $data->{chromosomes}->{$chromosome}->{genes}->{$gid};
			$item->{complete_name} = $hgene->{name}."(".$hgene->{xref}.")";
		$item->{name} = $hgene->{name};
		$item->{position} = $chromosome.":".$hgene->{start}."-".$hgene->{end};
		$item->{description} = $hgene->{description};
			$item->{variations} = $find;
		push(@out,$item);
		last;
	}
	

}

export_data::print($project,$cgi,\@out);
exit(0);
#my $patients =  $data->{patients};
#$all_data = $data->{chromosomes};


