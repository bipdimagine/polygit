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
my $bundle = $cgi->param('bundle');
my $project = $buffer->newProjectCache(-name=>$project_name);

my $panel_name = $cgi->param('panel');
my @out;
unless ($panel_name) {
#$project->setPanel("$panel") if $panel;
my $captures = $project->getCaptures();
my $vquery = validationQuery->new(dbh=>$buffer->dbh,capture_name=>$captures->[0]->validation_db());

my $dejavu;
my $transcripts_cgi  = $project->bundle_transcripts();
foreach my $capture (@$captures){
my $bundles = $buffer->getQuery()->getBundleFromCapture($capture->id);
foreach my $b (keys %$bundles){
	my $item;
	$item->{label} = $bundles->{$b}->{name};
	$item->{name} = $bundles->{$b}->{name};
	$item->{type} = "group";
	my $ts = $buffer->getQuery()->getTranscriptsFromBundle($bundles->{$b}->{name});
	my %dejavu;
	foreach my $t (keys %$ts){
		my $gene_name = $ts->{$t}->{gene};
		next if exists $dejavu{$gene_name};
		$dejavu{$gene_name} ++;
		my $gene;
		
	 	$gene->{name} = $gene_name;
	  	$gene->{type} = "gene";
	  	$gene->{toto} = "gene";
	  	 push(@{$item->{children}},$gene);
	}
	push(@out,$item);
}
}

export_data::print($project,$cgi,\@out);
exit(0);
}

my $panel = $project->getPanel("$panel_name");
my $bundles = $panel->getBundles();
 
foreach my $bundle (@$bundles){
	my $item;
	$item->{label} = $bundle->name;
	$item->{name} = $bundle->name;
	$item->{type} = "group";
	my $transcripts = $bundle->getTranscripts();
	foreach my $t (@$transcripts){
		my $gene_name = $t->getGene->name;
		my $gene;
	 	$gene->{name} = $gene_name;
	  	$gene->{type} = "gene";
	  	$gene->{toto} = "gene";
	  	 push(@{$item->{children}},$gene);
	}
	push(@out,$item);
	
}
export_data::print($project,$cgi,\@out);
exit(0);