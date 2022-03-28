#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";

use GBuffer;
use GenBoProject;
use Getopt::Long;
use Carp;
use strict;
use Set::IntSpan::Fast::XS;
use Data::Dumper;
use Data::Printer;
use Storable qw(store retrieve freeze);
use file_util;
use Term::Menus;
use IO::Prompt;

#------------------------------------------------------------------------------------------------------------------
#   Calculer le pourcentage de bases couvertes  (> a un seuil)  sur l'ensemble des exons codants
#------------------------------------------------------------------------------------------------------------------

# Instanciation d'un objet buffer pour gérer les connexions aux bases de données  
my $buffer = GBuffer->new();		
my $projectName;
my $cov_min;
my $transcriptName;

my %hres;



# récupère les options  
GetOptions(
	'project=s'   => \$projectName,
	'cov_min=s' => \$cov_min,
	'transcript=s' => \$transcriptName
);

die("\n\nNo -project option... \nUsage : -project= -cov_min= -transcript=  (all or ENSXXX or ENSXXX,ENSXXX,ENSXXX )  Die...\n\n") unless ($projectName);
die("\n\nNo -transcript option... \nUsage : -project= -cov_min= -transcript=  (all or ENSXXX or ENSXXX,ENSXXX,ENSXXX ) Die...\n\n") unless ($transcriptName);
$cov_min = 500 unless($cov_min);

#  instanciation d'un objet project
my $project = $buffer->newProject(-name=>$projectName);




#  creation de la liste des patients du projet
my $patients = $project->getPatients();
confess("No patients") unless scalar(@$patients);



# recuperer la liste des transcripts
my @list_transcript;
if ($transcriptName eq "all")
{
	@list_transcript = @{$project->bundle_transcripts()};
}
else
{
	@list_transcript=split(/,/,$transcriptName);
}

my $nbPatients=scalar(@{$patients});
my $nbTranscript = scalar(@list_transcript);

my $nbok_on_project;
my $nbtot;	# nombtre total de base sur le projet

foreach my $transName (@list_transcript)
{
	#  obtenir le transcript sur le projet
	my $transcript = $project->newTranscript($transName);
	confess("transcript doesn't exist") unless ($transcript);
	
	my $pok_on_transcript;
	
	foreach my $pat (@{$patients})
	{
		my $nbok;	
		
		my $hcov = $transcript->getGene()->get_coverage($pat)->coverage_intspan($transcript->getSpanCoding);
		

			foreach my $val (@{$hcov->{array}})
			{
				$nbok++ if ($val > $cov_min);
				$nbtot++;
			}

			my $pok = ($nbok / $hcov->{nb})*100;
			print $transName."  ".$pat->name()."   ";
			printf("%1.2f\n", $pok); 	
		
			$pok_on_transcript += $pok;
			$nbok_on_project += $nbok;
			
	}
	print $transName." : pourcentage bases couvertes au minimum en ".$cov_min."X sur l'ensemble des patients = ";	
	printf("%1.2f\n", ($pok_on_transcript/$nbPatients)); 	
	print "\n";
}

print "Pourcentage de bases couvertes au minimum en ".$cov_min."X sur l'ensemble des transcripts et sur l'ensemble des patients = ";
printf("%1.2f\n", ($nbok_on_project/$nbtot)*100); 	





 	