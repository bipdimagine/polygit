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

use File::Basename;

use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";

use GBuffer;
use GenBoProject;
use GenBoCache;
use GenBoNoSqlIntervalTree;
#use layout;


#######################################################################
#  1) liste  l'ensemble des fichiers  project.dejavu
#  2) retrieve chaque table de hash correspondante et l'integre dans une table de hash globale 
#  3) pour tout les project_dejavu plus récents que le dejavuglobal : retrieve la hash et l'integre dans la hash globale
#  4) freeze de la table de hash gllobale
########################################################################

#my $cgi = new CGI;

# lister les fichiers allSV
#my $cmd = "ls /data-xfs/Manue/Test_SV/DejaVu/TransLoc/*.store";
#my @res = `$cmd`;

my $halldejavu;
my $nbPatient;
my $buffer = GBuffer->new();
my @releases = ("HG19");	

# pour acceder aux dejavu de chaque projet
my $dir = $buffer->config->{'deja_vu_SV'}->{root}.$releases[0]."/".$buffer->config->{'deja_vu_SV'}->{SVeq};
my $cmd = "ls $dir/projects/*.SVeqDejavu";
my @res = `$cmd`;

foreach my $TransLocFile (@res)
{
	chomp($TransLocFile);
	my $filename = basename($TransLocFile);
	my ($projectname,$rien) = split(/\./,$filename);
	
	warn $projectname;
	
	my $hdejavu = retrieve($TransLocFile) or die "Can't retrieve datas from ".$TransLocFile." !\n";

	foreach my $id (keys %{$hdejavu})
	{
			$halldejavu->{$id}->{$projectname} = $hdejavu->{$id};
	}
	
}

#my $file_alldejavu = "/data-xfs/Manue/Test_SV/DejaVu/TransLoc/TranslocDejavu.all";
#store(\ %{$halldejavu}, $file_alldejavu) or die "Can't store $file_alldejavu!\n";

#Pour stocker le dejavu global
my $file_alldejavu = $dir."/SVeqDejavu.all";
store(\ %{$halldejavu}, $file_alldejavu) or die "Can't store $file_alldejavu!\n";

 	