#!/usr/bin/perl
use FindBin qw($Bin);
use lib $Bin;
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use GBuffer;
#use lib "/bip-d/perl/ensembl64/ensembl/modules";
use Bio::SearchIO;
use strict;
use Getopt::Long;
use Carp;
use JSON::XS;
use Getopt::Long; 
use List::MoreUtils qw(uniq);
#use ensembl_buffer;
use Storable qw/freeze thaw nfreeze nstore_fd nstore/;
use Set::IntSpan::Fast::XS ;
use Sys::Hostname;
use Storable qw/retrieve store/;
use  Set::IntSpan::Island;
use Parallel::ForkManager;
use Array::IntSpan;
use Data::Dumper;
use Array::IntSpan;
use GenBoNoSqlLmdb;
use GenBoNoSqlAnnotation;
use GenBoNoSqlDejaVu;
#use Bio::Ensembl::Registry;
use Storable qw(dclone);
use List::MoreUtils qw(uniq natatime);
use Set::IntervalTree;
use GenBoGene;
use GenBoTranscript;
use GenBoProtein;
use GenBoExon;
use GenBoIntron;
 use Storable;
my $sqliteDir ;#=  "/tmp/lmdb/annotation_75";
my $version;
GetOptions(
	'version=s' => \$version,
);


my $dir = "/data-isilon/public-data/repository/HG19/annotations/gencode.v$version/lmdb/";

my $annot2 =  GenBoNoSqlLmdb->new(name=>"genbo",dir=>$dir,mode=>"r",is_compress=>1);
my $annot_main =  GenBoNoSqlLmdb->new(name=>"main_transcripts",dir=>$dir,mode=>"c",is_compress=>1);
my $annot =  GenBoNoSqlAnnotation->new(dir=>$dir,mode=>"r");
my $type ="gene";
 my $z = $annot->get_like("synonyms","*".$type."*");

 foreach my $id (values %$z) {
 	warn $id;
 	my $obj = $annot2->get($id);
 	my $hg;
 	foreach my $tid (keys %{$obj->{transcripts_object}}){
 		my $tr = $annot2->get($tid);
 		my $tags = join(";",keys %{$tr->{tag}});
 		warn $tags;
 		 if ($tags =~ /CCDS/ or $tags =~ /appris_principal/ or $tags =~ /appris_alternative_1/ or $tags =~ /MANE/) {
 		 	$annot_main->put($tid,1);
 		 	$hg->{$tid} ++;
 		 }
 	}

	unless  ($hg) {
		foreach my $tid (keys %{$obj->{transcripts_object}}){
			$annot_main->put($tid,1);
				$hg->{$tid} ++;
		}
	}
			$annot_main->put($id,$hg);
 }
 
 $annot_main->close();