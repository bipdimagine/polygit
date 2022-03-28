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

#ce script est utilisé pour renvoyer les noms des projets d'un utilisateur après sa connection à l'interface ainsi que pour telecharger des fichiers via l'interface ou encore exporter les tableaux de données au format xls

my $cgi    = new CGI;

#chargement du buffer 
my $buffer = GBuffer->new;

encode_json_table_references($buffer);


sub encode_json_table_references {
	my $res = $buffer->getQueryPolypipeline()->getReferencesList();
	warn Dumper $res ;
	export_data::print_simpleJson($cgi,$res, "reference_id");
}
