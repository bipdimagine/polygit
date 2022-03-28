#!/usr/bin/perl
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";


use connect;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use export_data;
use Carp;
use JSON;

my $cgi    = new CGI;

#chargement du buffer 
my $buffer = GBuffer->new;

my $project_name   = $cgi->param('project');
my $project = $buffer->newProject( -name => $project_name );
die( "unknown project" . $project_name ) unless $project;


my $search = $cgi->param("name");
$search = uc($search);
	
 my $z = $project->liteAnnotations->get_like("annotations","$search*");
 $z ={} unless $z;
my $test = [];
foreach my $zz (keys %$z){
		next unless $zz =~/ENSG/;
		my $g = $z->{$zz};
		
		my $item;
		$item->{name} = $g->{external_name}."-".$g->{genbo_id};
	#	my $chr = $g->{chromosome};
		my ($ensg,$chr) = split("_",$g->{genbo_id});
		my $text = $g->{external_name}." - ".$g->{description}."-";
		$item->{label} = "<img src='/icons/Polyicons/8-em-check.png'></img> <b>$text</b> - chr$chr ";
		push(@$test,$item);
	}
#
print $cgi->header('text/json-comment-filtered');
		print encode_json $test;
exit(0);