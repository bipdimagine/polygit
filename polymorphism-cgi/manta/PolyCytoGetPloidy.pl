#!/usr/bin/perl

use Carp;
use strict;
use JSON;
use Data::Dumper;
use CGI qw/:standard :html3/;
use Set::IntSpan::Fast::XS;
use Set::IntervalTree;
use List::Util qw[min max];
use FindBin qw($Bin);
use Storable qw(store retrieve freeze);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
#
use GBuffer;
use GenBoProject;
use GenBoCache;

# Instanciation d'un objet buffer pour gérer les connexions aux bases de données  
my $buffer = GBuffer->new();		

# recupere les options  
my $cgi = new CGI;
my $project = $cgi->param('project');
my $patient_name = $cgi->param('patient');

die("\n\nNo -project option... Die...\n\n") unless ($project);

#  instanciation d'un objet project_cache et des objets patient et chromosome
my $project = $buffer->newProjectCache( -name => $project, -typeFilters=>"");
my $patient = $project->getPatient($patient_name);

my $hres;
my @listHashRes;

foreach my $chr ( @{$project->getChromosomes()} )
{
	#next  if ($chr->not_used());
	
	my $chrname = $chr->name();
	next if $chrname eq "MT";
	
	my $id = $chrname;
	$id = 23 if ($chrname eq "X");
	$id = 24 if ($chrname eq "Y");
	
	
	my ($score,$z) = $patient->ploidy_value($chr);
	my $style="box-shadow: 1px 1px 6px grey;background-color:white; color:black; border:1px solid orange; font-size:16px;";
	
	# gain
	$style="box-shadow: 1px 1px 6px black;background-color:blue; color:white; border:2px solid orange; font-size:18px;" if ($score > 0.1);
	
	# perte
	$style="box-shadow: 1px 1px 6px black;background-color:#DD4132; color:white; border:2px solid gray; font-size:18px;" if ($score < -0.5);
	

    $hres->{$id}->{'chr'} = int($id);
    $hres->{$id}->{'chrmean'} = "<td style=\"border-top-color:white;padding:5px;\"> <button class=\"btn btn-default btn-xs\" id=bchr".$id." onclick=\"launch_plotChr(".$id.");\" style=\"".$style."\">".$chrname."</button> </td>";
}

# pour le json
if (scalar(keys(%{ $hres})) == 0) 
{ 
				my $hash;
				$hash->{'chr'} = "-";
				$hash->{'chrmean'} = "-";
			 	 push(@listHashRes, $hash);
}
else
{
	foreach my $id ( keys %{$hres} ) 
	{
			push( @listHashRes, { %{ $hres->{$id} } } );
	}
}

printJson( \@listHashRes );
exit(0);


###############################################################################
#
#	methodes
#
#################################################################################


sub printJson {
	my ($listHash) = @_;
	my $hash;
	my @t = sort {$a->{'chr'} <=> $b->{'chr'}}  @$listHash;
	
	$hash->{'items'}      = \@t;
	
	print $cgi->header('text/json-comment-filtered');
	print encode_json $hash;
	print "\n";
}
	


