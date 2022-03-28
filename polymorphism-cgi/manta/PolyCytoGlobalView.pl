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
use Number::Format qw(:subs);

use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";

use GBuffer;
use GenBoProject;
use GenBoCache;

use layout;
use export_excel;
use export_data;

###########################################################################################################
#  1) recupere les infos des fichier bed pour les CNV et les BND
#  2) construit une table de hash unique
#  3) construit le Json resultant pour affichage
############################################################################################################

my $cgi = new CGI;
my $TheProjectName = $cgi->param('projectname');
my $thePatientName     = $cgi->param('filename');

# pour construire le projet
my $buffer  = GBuffer->new();
my $project = $buffer->newProjectCache( -name => $TheProjectName );
# et les paths
my $pathCNV = $project->getCNVDir();
my $pathBND = $project->getSVeqDir();

# pour les resultats
my $hres;
my @listHashRes;

# pour lire les fichiers bed 
my $fd1;
my $bedFileCNV  = $pathCNV.$thePatientName.".CNV.bed";
#open($fd1, '<' , $bedFileCNV) or die("open: $!"); 
open($fd1, '<' , $bedFileCNV) if (-f $bedFileCNV); 

my $fd2;
my $bedFileBND  = $pathBND.$thePatientName.".SVeq.bed";
open($fd2, '<' , $bedFileBND) or die("open: $!");


# lecture des fichiers et stockage dans une table de hash
my $id=1;
my $ligne;
while( defined( $ligne  = <$fd1> ) ) 
{
 	  chomp $ligne;
 	  last if ($ligne eq "no results");
   
  	 my ($chr,$pos,$event_id,$cytob,$djv) = split(/ /,$ligne);
   
 	 my ($type,$c,$d,$f) = split(/_/,$event_id);
  	 $c = "X" if ($c == 23);
   	$c = "Y" if ($c == 24);
   
  	 $hres->{$id}->{'id'} = $id;
     $hres->{$id}->{'chr'} = $chr;
   	 $hres->{$id}->{'pos'} = $pos;
  	 $hres->{$id}->{'cytob'} = $cytob;
   	 $hres->{$id}->{'djv'} = $djv;
  	 $hres->{$id}->{'event'} = $type."(".$c.$cytob.") : ".format_number($d)." - ".format_number($f);
  	 $id++;
}


while( defined( $ligne  = <$fd2> ) ) 
{
   chomp $ligne;
  
   my ($chr,$pos,$event_id,$cytob,$djv) = split(/ /,$ligne);
   my ($c1,$p1,$c2,$p2) = split(/_/,$event_id);
   
   $c1 = "X" if ($c1 == 23);
   $c1 = "Y" if ($c1 == 24);
   $c2 = "X" if ($c2 == 23);
   $c2 = "Y" if ($c2 == 24);
   
   my $type;
   $type  = "t(".$c1.";".$c2.")" if ($c1 != $c2);
 	$type  = "inv(".$c1.")" if ($c1 == $c2);
 	
   $hres->{$id}->{'id'} = $id;
   $hres->{$id}->{'chr'} = int($chr);
   $hres->{$id}->{'pos'} = int($pos); 
   $hres->{$id}->{'cytob'} = $cytob;
   $hres->{$id}->{'djv'} = $djv;
   $hres->{$id}->{'event'} = $type."(".$cytob.") : ".format_number($p1)." - ".format_number($p2);
   $id++;
}

# pour le json
if (scalar(keys(%{ $hres})) == 0) 
{ 
				my $hash;
				$hash->{'id'} = "-";
				$hash->{'pos'} = "-";
				$hash->{'chr'} = "-";
				$hash->{'event'} = "No results";
				$hash->{'type'} = "-";
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
	my @t;

	@t = sort { $a->{chr} <=> $b->{chr} or $a->{pos} <=> $b->{pos} } @$listHash;
			
	$hash->{'identifier'} = 'id';
	$hash->{'label'}      = 'id';
	$hash->{'items'}      = \@t;
	
	print $cgi->header('text/json-comment-filtered');
	print encode_json $hash;
	print "\n";
}

