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

my $dir = $buffer->config->{'deja_vu_SV'}->{root}.$releases[0]."/".$buffer->config->{'deja_vu_SV'}->{CNV};#;

#my $nodejavu = GenBoNoSqlDejaVu->new( dir => $dir, mode => "r" );
#5-42628373-cnv-del-2751
#5-42626542-cnv-del-6398
#1-152276681-cnv-dup-4916
#warn Dumper $nodejavu->get_cnv("5",42628373,42631373,"DEL");
#14-105415224-del-523
#2-179308111-cnv-del-1037
#5-140563576-cnv-del-9979
#1-152276681-cnv-dup-4916
#warn Dumper $nodejavu->get_cnv("X",215988,216489,"DUP");


#die();


my $cmd = "ls $dir/projects/*.dejavu";
my @res = `$cmd`;

my $nbok;
my $nbSV;

my $halldejavu;
#my $nbCNV=0;


my $htree_dejavu;
my $ids;
my $total;
my $xx =0;
my $projects_dv;
my %exclude;
$exclude{NGS2014_0001} = 1;
 map{$projects_dv->{$_} ++} grep { !( exists $exclude{$_} ) } @{ $buffer->listProjectsForDejaVu() };
 
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
	$xx ++;
	#last if $xx > 50;
}






my $nodejavu = GenBoNoSqlDejaVu->new( dir => $dir, mode => "c" );
warn $dir;

	 
foreach my $chr (keys %{$ids})
{
	$nodejavu->create_table($chr);
	my $sth = $nodejavu->dbh($chr)->prepare(
		'insert into  __DATA__(_key,_value,start,end,variation_type,ho,projects)  values(?,?,?,?,?,?,?) ;') or die $DBI::errstr;
		$sth->execute();
		foreach my $type  (keys %{$ids->{$chr}}) 
		{
				warn $type;
				my $tree;
				foreach my $id  (keys %{$ids->{$chr}->{$type}})
				{
								my ( $t, $c, $d, $f ) = split( /_/,$id);
								my $text;
								$sth->execute($id,$nodejavu->encode($total->{$id}),$d,$f,$type,0,0);
								push(@$tree,[$id,$d,$f]);
								my $id2 = $nodejavu->dbh($chr)->sqlite_last_insert_rowid();
								#warn Dumper $total->{$id};
	  							$nodejavu->sth_insert_position_cached($chr)->execute($id2,$d,$f) ;
				}
		}
$nodejavu->dbh($chr)->do(qq{CREATE UNIQUE INDEX if not exists _key_idx  on __DATA__ (_key);});
$nodejavu->dbh($chr)->do(qq{CREATE  INDEX if not exists _start_idx  on __DATA__ (start);});
$nodejavu->dbh($chr)->do(qq{CREATE  INDEX if not exists _end_idx  on __DATA__ (end);});
$nodejavu->dbh($chr)->do(qq{CREATE  INDEX if not exists _type_idx  on __DATA__ (variation_type);});

}
$nodejavu->close();





#warn $nbCNV;



#my $file_alldejavu = "/data-xfs/Manue/Test_SV/DejaVu/newVersion/dejavu_allCNV";
#store(\ %{$halldejavu}, $file_alldejavu) or die "Can't store $file_alldejavu!\n";



 	