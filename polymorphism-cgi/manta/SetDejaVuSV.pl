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
use layout;
#use export_excel; 
#use export_data;


#######################################################################
#  1) liste  l'ensemble des fichiers  project.dejavu
#  2) retrieve chaque table de hash correspondante et l'integre dans une table de hash globale 
#  3) pour tout les project_dejavu plus récents que le dejavuglobal : retrieve la hash et l'integre dans la hash globale
#  4) freeze de la table de hash gllobale
########################################################################

#my $cgi = new CGI;

# lister les fichiers allSV
my $cmd = "ls /data-xfs/Manue/Test_SV/DejaVu/newVersion/*.dejavu";
my @res = `$cmd`;

my $nbok;
my $nbSV;

my $halldejavu;
#my $nbCNV=0;


my $htree_dejavu;
my $ids;
my $total;

foreach my $CNVFile (@res)
{
	chomp($CNVFile);
	my $filename = basename($CNVFile);
	my ($projectname,$rien) = split(/\./,$filename);
	
	warn $projectname;

	my $hdejavu = retrieve($CNVFile) or die "Can't retrieve datas from ".$CNVFile." !\n";
 	#my $no = GenBoNoSqlIntervalTree->new(dir=>"/data-xfs/Manue/Test_SV/DejaVu/newVersion/" ,mode=>"c");
  
 
	foreach my $type (keys %{$hdejavu})
	{
		foreach my $num (keys %{$hdejavu->{$type}})
		{
			
				foreach my $id (keys %{$hdejavu->{$type}->{$num}})
				{		
					$ids->{$num}->{$type}->{$id} ++;
					$total->{$id}->{$projectname} = $hdejavu->{$type}->{$num}->{$id};
					#$halldejavu->{$type}->{$num}->{$id}->{$projectname} = $hdejavu->{$type}->{$num}->{$id};
					#$nbCNV++;
				}
		}
	}	
}

 my $lmdb = GenBoNoSqlLmdb->new(dir=>"/data-xfs/Manue/Test_SV/DejaVu/newVersion/",mode=>"c",name=>"dejavu_sv",is_compress=>1);
 
 foreach my $id (keys %{$total}){
	 $lmdb->put($id,$total->{$id}); 	
 }
$lmdb->close();


my $no = GenBoNoSqlIntervalTree->new(dir=>"/data-xfs/Manue/Test_SV/DejaVu/newVersion/" ,mode=>"c");

foreach my $chr (keys %{$ids})
{
		foreach my $type  (keys %{$ids->{$chr}}) 
		{
				my $tree;
				foreach my $id  (keys %{$ids->{$chr}->{$type}})
				{
								my ( $t, $c, $d, $f ) = split( /_/,$id);
								push(@$tree,[$id,$d,$f]);
				}
				$no->put("sv_dv_interval_tree",$chr."_".$type,$tree);
		}
}
$no->close();




#warn $nbCNV;



#my $file_alldejavu = "/data-xfs/Manue/Test_SV/DejaVu/newVersion/dejavu_allCNV";
#store(\ %{$halldejavu}, $file_alldejavu) or die "Can't store $file_alldejavu!\n";



 	