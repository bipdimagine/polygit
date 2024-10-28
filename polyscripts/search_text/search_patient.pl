#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
use GBuffer;
use Set::IntSpan::Fast::XS;
use Getopt::Long;
use Search::Indexer;
use Text::Levenshtein::XS qw(distance);
my $search;
GetOptions(
	'patient=s'       => \$search,
);

my $buffer = new GBuffer;

my $dbh = $buffer->dbh;

my $sql = qq{SELECT p.name patient ,pr.name,p.patient_id  pid FROM PolyprojectNGS.patient as p,PolyprojectNGS.projects as pr where pr.project_id=p.project_id;};
my $sth = $dbh->prepare($sql);
$sth->execute();
my $h = $sth->fetchall_hashref('pid');
my $res;
my $xx; 
my $hid;
foreach my $pid (keys %$h) {
 my $p = $h->{$pid}->{patient};
 push(@{$hid->{$p}},$pid);
if ($p eq $search){
		 $res->{$p} = -1;
		 next;
}
if ($p =~/$search/){
	$res->{$p} = 0;
		 next;
}
 my $distance = distance($search,$p);
 next if $distance > 5;
 $res->{$p} = $distance;
}
unless (keys %$res){
	print $search ." NOT FOUND \n";
	exit(0);
}
my @sort = sort{$res->{$a} <=> $res->{$b}} keys %$res;
my $min = $res->{$sort[0]};

my $nearest;
foreach my $p (@sort){
	#next if $res->{$p} ne -1;
	last if $res->{$p}>$min;
	push(@$nearest,$p)
	#warn $p." ".$res->{$p};
	#last;
}
my $text;


#if ($min == -1) {
#	
#	foreach my $p (@sort){
#		next if $res->{$p} ne -1;
#		foreach my $pid (@{$hid->{$p}}){
#			 push(@$text,$p."\t".$h->{$pid}->{project}."\t".$min);
#		}
#	}
#	
#		
#}
#exit(0);
#my $nearest;
#foreach my $p (@sort){
#	next if $res->{$p} ne -1;
#	#last if $res->{$p}>$min;
#	push(@$nearest,$p)
#	#warn $p." ".$res->{$p};
#	#last;
#}
#my $text;

foreach my $near (@$nearest){
	foreach my $pid (@{$hid->{$near}}){
	 push(@$text,$search."\t".$near."\t".$h->{$pid}->{name});
	}
}
if (scalar(@$text)>2){
	print join("\t",@$text)." ***\n" ;
	
}
else {
	print join("\t",@$text)."\n" ;
}
	

#my $indexer = Search::Indexer->new();
#my $z = 0;
#foreach my $p (keys %$h){
#	$indexer->add($z++, $h->{$p}->{patient});
#	
#}
#
#my $search ;
## Ajout de documents à l'index
##$indexer->build_index();
#
## Recherche de termes
#my @results = $indexer->search("");
#
#warn Dumper @results;
#
#
