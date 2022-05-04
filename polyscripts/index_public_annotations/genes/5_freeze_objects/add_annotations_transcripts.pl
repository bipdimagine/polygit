#!/usr/bin/perl
use FindBin qw($Bin);
use lib $Bin;
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../packages/";
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
use Bio::ToolBox;
use Bio::ToolBox::parser::gff;
use Bio::ToolBox::GeneTools qw(:all);
 
 
my $sqliteDir ;#=  "/tmp/lmdb/annotation_75";
my $version;
GetOptions(
	'version=s' => \$version,
);


my $dir =  "/tmp/lmdb/$version/annotations/";
#my $notools =  GenBoNoSqlLmdb->new(name=>"genetools",dir=>$dir,mode=>"r",is_compress=>1); 
my $annot2 =  GenBoNoSqlLmdb->new(name=>"genbo",dir=>$dir,mode=>"r",is_compress=>1);
my $annot_main =  GenBoNoSqlLmdb->new(name=>"main_transcripts",dir=>$dir,mode=>"c",is_compress=>1);
my $short_annot_main =  GenBoNoSqlLmdb->new(name=>"shortlist_main_transcripts",dir=>$dir,mode=>"c",is_compress=>1);
#my $annot =  GenBoNoSqlAnnotation->new(dir=>$dir,mode=>"r");
my $annot =  GenBoNoSqlAnnotation->new(dir=>$dir,mode=>"r");
my $type ="gene";
 my $z = $annot->get_like("synonyms","*"."gene*");
my $total;
my %all;
 foreach my $id (values %$z) {
 	my $obj = $annot2->get($id);
 	my $hg;
 	my $tr_mains;
 	my $tr_synonyms;
 	my $tlevel;
 	foreach my $tid (keys %{$obj->{transcripts_object}}){
 		my $tr = $annot2->get($tid);
 		
 		my $span = $tr->getSpanCoding;
 		$span = $tr->getGenomicSpan if $span->is_empty;
 		my $to = $annot->get("annotations",$tid);
 		my $tl = $to->{'transcript_support_level'};
 		$tl = 6 if $tl eq "NA";
		$tl = 6 if $tl eq "Missing";
 		push(@{$tlevel->{$tl}},$tid);
 		my $tags = join(";",keys %{$tr->{tag}});
 		 if ($tags =~ /CCDS/ or $tags =~ /appris_principal/ or $tl == 1) {
 		 	
 		 	push(@{$tr_mains->{$span->as_string}},$tr);
 		 }
 		# if ( $tags =~ /appris_alternative_1/ && $tl >= 2){
 		 #	push(@{$tr_mains->{$span->as_string}},$tr);
 		 #}
 	}
 	unless  ($tr_mains) {
 		foreach my $tid (keys %{$obj->{transcripts_object}}){
 		my $tr = $annot2->get($tid);
 		
 		my $span = $tr->getSpanCoding;
 		$span = $tr->getGenomicSpan if $span->is_empty;
 		my $to = $annot->get("annotations",$tid);
 		my $tl = $to->{'transcript_support_level'};
 		$tl = 6 if $tl eq "NA";
		$tl = 6 if $tl eq "Missing";
 		#push(@{$tlevel->{$tl}},$tid);
 		my $tags = join(";",keys %{$tr->{tag}});
 	
 		 if ( $tags =~ /appris_alternative/ and $tl <= 3){
 		 	push(@{$tr_mains->{$span->as_string}},$tr);
 		 
		}
 		}
 	}

	unless  ($tr_mains) {
		
 		my @t  =sort {$a <=> $b} keys %$tlevel;
 		my $min = $t[0];
# 		warn scalar(@{$tlevel->{$min}}).' '.$min;
 		
		foreach my $tid ( @{$tlevel->{$min}}) {
#			warn $tid;

			my $tr = $annot2->get($tid);
			if ( $tr->{biotype} =~/processed_pseudogen/){
				next;
			}
			warn scalar(@{$tlevel->{$min}}).' '.$min.' '.$tid;
 			my $span = $tr->getSpanCoding;
 			$span = $tr->getGenomicSpan if $span->is_empty;
			#$annot_main->put($tid,1);
			#	$hg->{$tid} ++;
				push(@{$tr_mains->{$span->as_string}},$tr);
		}
	}
	unless  ($tr_mains) {
		next;
	}
	foreach my $k (keys %$tr_mains){
		my $id = $tr_mains->{$k}->[0]->id;
		$hg->{$id} ++;
		$annot_main->put($id,1);
		$all{$id}++;
		$total ++;	
	}
	$annot_main->put($id,$hg);
	
		#	$annot_main->put($id,$hg);
 }
 $annot_main->put('transcripts',[keys %all]);
 warn $short_annot_main->filename;
 warn $total;
 $short_annot_main->close();
 $annot_main->close();