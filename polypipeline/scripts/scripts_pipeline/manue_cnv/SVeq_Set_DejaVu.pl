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
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../../packages/export";
use lib "$Bin/../../../../packages/layout";

use GBuffer;
use GenBoProject;
use GenBoCache;
use GenBoNoSqlIntervalTree;



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
my @releases = ("HG19","HG38");	

# pour acceder aux dejavu de chaque projet
my $dir = $buffer->config->{'deja_vu_SV'}->{root}.$releases[0]."/".$buffer->config->{'deja_vu_SV'}->{SVeq};
my $cmd = "ls $dir/projects/*.SVeqDejavu";
my @res = `$cmd`;
my $nb;
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
	
	$nb ++;
	last if $nb >50;
}

#my $file_alldejavu = "/data-xfs/Manue/Test_SV/DejaVu/TransLoc/TranslocDejavu.all";
#store(\ %{$halldejavu}, $file_alldejavu) or die "Can't store $file_alldejavu!\n";

#Pour stocker le dejavu global
my $file_alldejavu = $dir."/SVeqDejavu.all";
my @pos;
foreach my $id (keys %$halldejavu) {
	my ($c1,$p1,$c2,$p2) = split("_",$id);
	push(@pos,join("\t",$c1,$p1,$p1+1,$id));
	push(@pos,join("\t",$c2,$p2,$p2+1,$id));
}


my $tmp_dir = ".";
my $name = "tmp.".time.".".rand(time);
my $f1 = "$tmp_dir/"."$name.bed";
my $f2 = "$tmp_dir/"."$name.out.bed";
 my $f3 = "$tmp_dir/"."$name.error.bed";
 
  
 write_file("$tmp_dir/"."$name.bed", @pos);
 
 system("/software/distrib/ucsc_util/liftOver $f1 $RealBin/../hg19ToHg38.over.chain $f2 $f3 ");
 
   open (BED,"$f2");
  
   my @unsave;
   my $newTotal;
   my $boundary;

   while (my $line = <BED>){
   	chomp($line);
   	my ($chr,$start,$end,$id) = split(" ",$line);
   
   	
   	$chr =~ s/chr//;
   	$chr = "MT" if ($chr eq "M");
   	push({$boundary->{$id}->{chr}},$chr);
   	push({$boundary->{$id}->{start}},$start);
   
   }
   
    unlink $f1;
 unlink $f2;
 unlink $f3;
 
 warn Dumper $boundary;
 
   
   die();
   foreach my $id (keys %$boudary){
   	next if scalar @($boundary->{$id}->{start}) < 2;
   	my $new_id = $boundary->{$id}->{chr}->[0]."_".$boundary->{$id}->{start}->[0]."_".$boundary->{$id}->{chr}->[1]."_".$boundary->{$id}->{start}->[1];
   		$newTotal->{$new_id} = $halldejavu->{$id};
   }
   
 $dir = $buffer->config->{'deja_vu_SV'}->{root}.$releases[1]."/".$buffer->config->{'deja_vu_SV'}->{SVeq};
 $cmd = "ls $dir/projects/*.SVeqDejavu";
 @res = `$cmd`;

foreach my $TransLocFile (@res)
{
	chomp($TransLocFile);
	my $filename = basename($TransLocFile);
	my ($projectname,$rien) = split(/\./,$filename);
	
	
	my $hdejavu = retrieve($TransLocFile) or die "Can't retrieve datas from ".$TransLocFile." !\n";

	foreach my $id (keys %{$hdejavu})
	{
			$newTotal->{$id}->{$projectname} = $hdejavu->{$id};
	}
	
}
   
 
 



die();
store(\ %{$halldejavu}, $file_alldejavu) or die "Can't store $file_alldejavu!\n";

 	