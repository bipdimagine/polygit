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

use GBuffer;
use GenBoProject;
use GenBoCache;

# Instanciation d'un objet buffer pour gérer les connexions aux bases de données  
my $buffer = GBuffer->new();		

# recupere les options  
my $cgi = new CGI;
my $project_name = $cgi->param('project');
my $patient_name = $cgi->param('patient');

die("\n\nNo -project option... Die...\n\n") unless ($project_name);

#  instanciation d'un objet project_cache et des objets patient et chromosome
my $project = $buffer->newProjectCache( -name => $project_name, -typeFilters=>"");
my $child = $project->getPatient($patient_name);


next unless $child->isChild;	# only care of child
my $fam = $child->getFamily();
my $mother = $fam->getMother();
my $father = $fam->getFather();

next unless $father;	# child with parents
next unless $mother;

my $hres;
my @listHashRes;
my $type;

my $dir = $project->getSVDir()."/UPD";

if (-e "$dir/".$child->name.".json")
{


	open(JSON,"$dir/".$child->name.".json");

	my $tt = <JSON>;
	chomp($tt);
	my $hash = decode_json $tt;

  foreach my $chr (keys %$hash)
  { 
	next if $chr eq "X" && $child->sex == 1;
	next if $chr eq "Y" && $child->sex == 2;
	next if $chr eq "MT";
	
	my $name = $child->name;
	my $style;
	my $etiq;
	my $id = $chr;


	my $max  = -10;
	my $show = 0;


# Detection de larges deletions

	if (( $hash->{$chr}->{"TRANSMIS_BYF"} > 10 )	&& ( $hash->{$chr}->{"TRANSMIS_BYM"} > 10 )) # pas d'unidisomie
	{
		if (( $hash->{$chr}->{"ISO_UPD_M"} > 10 ) && ( $hash->{$chr}->{"HoFrom_M"} > 50 ) )
		{
			$style="box-shadow: 1px 1px 2px black;background-color:lightblue; color:white; border:2px solid white; font-size:16px;";
			$show = 1;
			$type = "Suspicion of large deletion "; # deletion des alleles paternels
		}
		if (( $hash->{$chr}->{"ISO_UPD_F"} > 10 ) && ( $hash->{$chr}->{"HoFrom_F"} > 50 ) )
		{
			$style="box-shadow: 1px 1px 2px black;background-color:lightblue; color:white; border:2px solid white; font-size:16px;";
			$show = 2;
			$type = "Suspision of large deletion "; # deletion des alleles maternels
		}
	}
	else
	{
		if ( $hash->{$chr}->{"TRANSMIS_BYF"} < 10 )		# si moins de 10% des variants du père transmis : unidisomie Maternelle
		{
			$style="box-shadow: 1px 1px 2px black;background-color:pink; color:white; border:2px solid white; font-size:16px;";
			$show = 2;
			if ( $hash->{$chr}->{"ISO_UPD_M"} > 30 )		# si plus de 30% des variants He de la mère sont Ho chez l'enfant : isodisomie 
			{
				$type = " Suspicion of UPD (Isodisomy) ";
			}
			else
			{
				$type =  " Suspicion of UPD (Heterodisomy) ";
			}
		}
		else
		{
			if ( $hash->{$chr}->{"TRANSMIS_BYM"} < 10 )		# si moins de 10% des variants de la mère transmis : unidisomie Paternelle
			{
				$style="box-shadow: 1px 1px 2px black;background-color:lightblue; color:white; border:2px solid white; font-size:16px;";
				$show = 1;
				my $val = $hash->{$chr}->{"TRANSMIS_BYM"};
				if ( $hash->{$chr}->{"ISO_UPD_F"} > 30 )		# si plus de 30% des variants He du père sont Ho chez l'enfant : isodisomie 
				{
					$type = " Suspicion of UPD (Isodisomy) ";
				}
				else
				{
					$type = " Suspicion of UPD (Heterodisomy) ";
				}
			}
		}
	}
	
	
	if($id eq "X") {$id=23;}
   	if($id eq "Y") {$id=24;}
   	
	
	if($show == 1)
	{
		$hres->{$id}->{'chr'} = int($id);
   		$hres->{$id}->{'value'} = "<td style=\"border-top-color:white;padding:5px;\"> <label style=\"font-size:12px; color:orange;\">".$type."<\label><button class=\"btn btn-default btn-xs\" id=bchr".$id." onclick=\"launch_plot_BA(".$id.");\" style=\"".$style."\">".$chr."<img src=\"../../images/polyicons/icons8-person-24.png\"></button> </td>";
	}
	if($show == 2)
	{
		$hres->{$id}->{'chr'} = int($id);
   		$hres->{$id}->{'value'} = "<td style=\"border-top-color:white;padding:5px;\"> <label style=\"font-size:12px; color:orange;\">".$type."<\label> <button class=\"btn btn-default btn-xs\" id=bchr".$id." onclick=\"launch_plot_BA(".$id.");\" style=\"".$style."\">".$chr."<img src=\"../../images/polyicons/icons8-person-female-24.png\"></button> </td>";
	}
	
 } # boucle sur les chromosome

} # fin du if
else
{
	$type = "UPD not performed for this patient. Contact Bioinformatics plateforme";
	$hres->{"0"}->{'chr'} = "-";
   	$hres->{"0"}->{'value'} = "<td style=\"border-top-color:white;padding:5px;\"> <label style=\"font-size:12px; color:orange;\">".$type."<\label></td>";
}

# pour le json
if (scalar(keys(%{$hres})) == 0) 
{ 
		my $hash;
		$hash->{'chr'} = "0";
		$hash->{'value'} = "<td style=\"border-top-color:white;padding:5px;\"> <label style=\"font-size:12px; color:orange;\"> UniParental Disomy : No results <\label>";
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
	
	$hash->{'items'} = \@t;
	
	print $cgi->header('text/json-comment-filtered');
	print encode_json $hash;
	print "\n";
}
	


