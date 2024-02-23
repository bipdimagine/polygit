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
#my $patient_name = $cgi->param('patient');

die("\n\nNo -project option... Die...\n\n") unless ($project_name);

#  instanciation d'un objet project_cache et récupération des données familles
my $project = $buffer->newProjectCache( -name => $project_name, -typeFilters=>"");

# boucle sur les patients du projets
my $listPatients = $project->getPatients();

foreach my $patient (@$listPatients)
{
			my $prg = qq{$Bin/../../polypipeline/scripts/scripts_pipeline/upd/getUPD.pl};
			system("perl $prg -project=$project_name -patient=".$patient->name." 2>/dev/null >/dev/null");

}

foreach my $child (@$listPatients)
{
	next unless($child->isChild);

	my $fam = $child->getFamily();
	my $mother = $fam->getMother();
	my $father = $fam->getFather();

	my $dir = $project->getSVDir()."/UPD";
	
	if ($child->isChild && $father && $mother)  # child with both parents
	{
			open(JSON,"$dir/".$child->name.".json");

			my $tt = <JSON>;
			chomp($tt);
			my $hash = decode_json $tt;

 			foreach my $chr (keys %$hash)
 			{ 
 		
				next if $chr eq "X" && $child->sex == 1;
				next if $chr eq "Y";
				next if $chr eq "MT";
	
				my $name = $child->name;
				
				# Detection de larges deletions

				if (( $hash->{$chr}->{"TRANSMIS_BYF"} > 10 )	&& ( $hash->{$chr}->{"TRANSMIS_BYM"} > 10 )) # pas d'unidisomie
				{
					if (( $hash->{$chr}->{"ISO_UPD_M"} > 10 ) && ( $hash->{$chr}->{"HoFrom_M"} > 50 ) )
					{
						print $project_name." ".$name."/chr".$chr." : Suspicion of Paternal Large Deletion\n"; # deletion des alleles paternels
					}
					if (( $hash->{$chr}->{"ISO_UPD_F"} > 10 ) && ( $hash->{$chr}->{"HoFrom_F"} > 50 ) )
					{
						print $project_name." ".$name."/chr".$chr." : Suspicion of Maternal Large Deletion\n"; # deletion des alleles maternels
					}
				}
				else
				{
					if ( $hash->{$chr}->{"TRANSMIS_BYF"} < 10 )		# si moins de 10% des variants du père transmis : unidisomie Maternelle
					{
						if ( $hash->{$chr}->{"ISO_UPD_M"} > 30 )		# si plus de 30% des variants He de la mère sont Ho chez l'enfant : isodisomie 
						{
							print $project_name." ".$name." chr".$chr." : Suspicion of Maternal Isodisomy\n";
						}
						else
						{
							print $project_name." ".$name."/chr".$chr." : Suspicion of  Maternal Heterodisomy\n";
						}
					}
					else
					{
						if ( $hash->{$chr}->{"TRANSMIS_BYM"} < 10 )		# si moins de 10% des variants de la mère transmis : unidisomie Paternelle
						{
							if ( $hash->{$chr}->{"ISO_UPD_F"} > 30 )		# si plus de 30% des variants He du père sont Ho chez l'enfant : isodisomie 
							{
								print $project_name." ".$name."/chr".$chr." : Suspicion of Paternal Isodisomy\n";
							}
							else
							{
								print $project_name." ".$name."/chr".$chr." : Suspicion of Paternal Heterodisomy\n";
							}
						}
					}
				}
			 } # boucle sur les chromosome
	} # fin du si enfant avec parents
} # fin de la boucle sur les patients
exit(0);



