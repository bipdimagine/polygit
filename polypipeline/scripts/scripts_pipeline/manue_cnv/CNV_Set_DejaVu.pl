#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
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


use GBuffer;
use GenBoProject;
use GenBoCache;
use GenBoNoSqlIntervalTree;


#######################################################################
#  1) liste  l'ensemble des fichiers  project.dejavu
#  2) retrieve chaque table de hash correspondante et l'integre dans une table de hash globale 
#  3) pour tout les project_dejavu plus r??cents que le dejavuglobal : retrieve la hash et l'integre dans la hash globale
#  4) freeze de la table de hash gllobale
########################################################################

#my $cgi = new CGI;
my $buffer = GBuffer->new();
my @releases = ("HG19");

# lister les fichiers allSV
system("$Bin/new_dejavu_CNV.pl ");
my $dir = $buffer->config->{'deja_vu_SV'}->{root}.$releases[0]."/".$buffer->config->{'deja_vu_SV'}->{CNV};#;
my $cmd = "ls $dir/projects/*.dejavu";
my @res = `$cmd`;
my %exclude;
$exclude{NGS2014_0001} = 1;
my $projects_dv;
 map{$projects_dv->{$_} ++} grep { !( exists $exclude{$_} ) } @{ $buffer->listProjectsForDejaVu() };
 #warn Dumper $projects_dv;
# die();
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
	next unless exists $projects_dv->{$projectname};
	#warn $projectname;

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

 my $lmdb = GenBoNoSqlLmdb->new(dir=>$dir,mode=>"c",name=>"dejavu_sv",is_compress=>1);
 
 foreach my $id (keys %{$total}){
	 $lmdb->put($id,$total->{$id}); 	
 }
$lmdb->close();


my $no = GenBoNoSqlIntervalTree->new(dir=>$dir ,mode=>"c");

	 
foreach my $chr (keys %{$ids})
{
		foreach my $type  (keys %{$ids->{$chr}}) 
		{
				warn $type;
				my $tree;
				foreach my $id  (keys %{$ids->{$chr}->{$type}})
				{
								my ( $t, $c, $d, $f ) = split( /_/,$id);
								my $text;
								push(@$tree,[$id,$d,$f]);
								
				}
				$no->put("sv_dv_interval_tree",$chr."_".$type,$tree);
		}


}
$no->close();





#warn $nbCNV;



#my $file_alldejavu = "/data-xfs/Manue/Test_SV/DejaVu/newVersion/dejavu_allCNV";
#store(\ %{$halldejavu}, $file_alldejavu) or die "Can't store $file_alldejavu!\n";



 	