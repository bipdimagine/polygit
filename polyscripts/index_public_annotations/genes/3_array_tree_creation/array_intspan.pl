#!/usr/bin/perl
use FindBin qw($Bin);
use lib $Bin;
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
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
use GenBoNoSql;
use GenBoNoSqlAnnotation;
#use Bio::Ensembl::Registry;
use Storable qw(dclone);
use List::MoreUtils qw(uniq);
use Set::IntervalTree;
use Date::Tiny;
use Digest::MD5::File qw( file_md5_hex );
#use Tree::Interval::Fast;
#use Tree::Interval::Fast::Interval;
#use Tree::R;
use GenBoNoSqlIntervalTree;
use DBI;


my $sqliteDir ;
my $out;
my $type;
my $version;
GetOptions(
	'version=s' => \$version,
	'type=s' => \$type,
);
my $sqliteDir =  "/tmp/lmdb/$version/annotations";
die(" -dir= ") unless $sqliteDir;

warn "work on $sqliteDir";

intervaltree_freeze();
array_intspan() if $type eq "all";
exit(0);
rtree_freeze();
die();
warn "end freeze";
#die();
#array_intspan();

sub intervaltree_freeze {
	my $annot =  GenBoNoSqlAnnotation->new(dir=>$sqliteDir,mode=>"r");
	
		my $type = "gene";
 	my $z = $annot->get_like("annotations","*".$type."*");
 	my @genes_id = values %{$z};
 my $hintspan;
 my $htree;
 my $array_tree;
 my $array_limit;
 my $freeze;
 my $dd =0;
 
 foreach my $gene (@genes_id){
 	$dd ++;
 	warn $dd if $dd%100 == 0 ;
 	my $isncrna;
 	my $chr = $gene->{chromosome};

 	$hintspan->{genes}->{$chr} = [] unless exists $hintspan->{genes}->{$chr};
 	$hintspan->{transcripts}->{$chr} = [] unless exists $hintspan->{transcripts}->{$chr};
 	my $isncrna;
 	
 	foreach my $tid (@{$gene->{transcripts}}){
 		my ($id,$n) = split("_",$tid);
 		my $tr = $annot->get("annotations",$id."_".$chr);
 		$isncrna = 1 if $tr->{biotype} =~ /pseudo/ || $tr->{biotype} =~ /RNA/;
 		my $limit =11;
		$limit = 0 if $isncrna;
		$limit= 0 if $chr eq 'MT';
		my $start1 = $tr->{start}-$limit;
		my $end1 = $tr->{end}+$limit;
		push(@{$hintspan->{transcripts}->{$chr}},[$tr->{genbo_id},$start1,$end1+1]);
		my $intspan = $tr->{genomic_span};
		if ( $tr->{span_coding}){
		my $iter = $tr->{span_coding}->iterate_runs();
		my $protid  = $tr->{protein};
		unless ($protid =~/_/){
			$protid.="_".$tr->{chromosome};
		}
		while (my ( $from, $to ) = $iter->()) { 
			push(@{$hintspan->{exons}->{$chr}},["c:".$tr->{genbo_id},$from,$to+1]);
			push(@{$hintspan->{proteins}->{$chr}},[$protid,$from,$to+1]);
		}
		$intspan = $tr->{genomic_span}->diff($tr->{span_coding});
		}
		my $iter2 = $intspan->iterate_runs();
		while (my ( $from, $to ) = $iter2->()) { 
			push(@{$hintspan->{exons}->{$chr}},["u:".$tr->{genbo_id},$from,$to+1]);
		}
 	}
 	my $limit = 1000;
		$limit = 15 if $isncrna;
		$limit= 0 if $chr eq 'MT';
	my $gstart = $gene->{start} - $limit;
	$gstart =1 if $gstart < 0;
	my $gend = $gene->{end} + $limit;
	push(@{$hintspan->{genes}->{$chr}},[$gene->{genbo_id},$gstart,$gend+1]);
	
 }
  my $no = GenBoNoSqlIntervalTree->new(dir=>$sqliteDir,mode=>"c");
 foreach my $type  (keys %$hintspan){
 	foreach my $chr (keys %{$hintspan->{$type}}){
 
 		 $no->put($type,$chr,$hintspan->{$type}->{$chr});
 	}
 }	
}



sub array_intspan {
	
	my $annot =  GenBoNoSqlAnnotation->new(dir=>$sqliteDir,mode=>"r");

my $type = "gene";
 my $z = $annot->get_like("annotations","*".$type."*");
 #warn Dumper $annot->get("annotations","NP001257972");
 my @genes_id = values %{$z};
 my $hintspan;
 my $htree;
 my $array_tree;
 my $array_limit;
 my $freeze;
 foreach my $gene (@genes_id){
 	my $chr = $gene->{chromosome};
 	$hintspan->{genes}->{$chr} = Set::IntSpan::Fast::XS->new() unless exists $hintspan->{genes}->{$chr};
 	$htree->{genes}->{$chr} = Set::IntervalTree->new  unless exists $htree->{genes}->{$chr};
 
 	my $isncrna;
 	
 	foreach my $tid (@{$gene->{transcripts}}){
 		my ($id,$n) = split("_",$tid);
 		my $tr = $annot->get("annotations",$id."_".$chr);
 		$isncrna = 1 if $tr->{biotype} =~ /pseudo/ || $tr->{biotype} =~ /RNA/;
 		$hintspan->{transcripts}->{$chr} = Set::IntSpan::Fast::XS->new() unless exists $hintspan->{transcripts}->{$chr};
 		$htree->{transcripts}->{$chr} = Set::IntervalTree->new  unless exists $htree->{transcripts}->{$chr};
 		my $limit =11;
		$limit = 0 if $isncrna;
		$limit= 0 if $chr eq 'MT';
		my $start1 = $tr->{start}-$limit;
		my $end1 = $tr->{end}+$limit;
		$hintspan->{transcripts}->{$chr}->add_range($start1,$end1);
		$htree->{transcripts}->{$chr}->insert($tr->{genbo_id},$start1,$end1+1);
		push(@{$array_tree->{transcripts}->{$chr}},[$tr->{genbo_id},$start1,$end1]);
		push(@{$array_limit->{transcripts}->{$chr}},($start1,$end1+1));
		my $set1 = Set::IntSpan::Island->new($start1."-".$end1);
		$freeze->{transcripts}->{$chr}->{$tr->{genbo_id}} = $set1;
		
 	}
 	my $limit = 500;
		$limit = 15 if $isncrna;
		$limit= 0 if $chr eq 'MT';
	my $gstart = $gene->{start} - $limit;
	$gstart =1 if $gstart < 0;
	my $gend = $gene->{end} + $limit;
	$hintspan->{genes}->{$chr}->add_range($gstart,$gend);
	$htree->{genes}->{$chr}->insert($gene->{genbo_id},$gstart,$gend+1);
	push(@{$array_limit->{genes}->{$chr}},($gstart,$gend));
	my $set1 = Set::IntSpan::Island->new($gstart."-".$gend);
	$freeze->{genes}->{$chr}->{$gene->{genbo_id}} = $set1;
	push(@{$array_tree->{genes}->{$chr}},[$gene->{genbo_id},$gstart,$gend+1]);
 }

warn "STEP 1";
 my $no = GenBoNoSql->new(dir=>$sqliteDir,mode=>"c");
 
 foreach my $type  (keys %$hintspan){
 	warn $type;
 	#next unless $type eq "transcripts";
 	foreach my $chr (keys %{$hintspan->{$type}}){
 			my $covers = Set::IntSpan::Island->extract_covers($freeze->{$type}->{$chr});
 			my $foo = Array::IntSpan->new();
	foreach my $cover (@$covers) {
		
		next  if scalar(@{$cover->[1]}) == 0;	
		my $set = $cover->[0];	
		my $min = min $set;
		my $max = max $set;
		my @spans = spans $set;
		confess("problem ") unless scalar(@spans) ==1;
		$foo->set_range($spans[0][0],$spans[0][1],$cover->[1]);	
	}
	 $no->put("intspan_".$type,$chr,$foo);
 			
 	}
 }
warn "END";
}






sub rtree_sqlite {
	my $annot =  GenBoNoSqlAnnotation->new(dir=>$sqliteDir,mode=>"r");

	my $db_genes = DBI->connect("dbi:SQLite:dbname=$sqliteDir/genes_rtree_db.lite","","");
	my $db_transcripts = DBI->connect("dbi:SQLite:dbname=$sqliteDir/transcripts_rtree_db.lite","","");
	 $db_genes->do(qq{
		CREATE VIRTUAL TABLE positions USING rtree_i32(
   id,              -- Integer primary key
   start, end,      -- Minimum and maximum X coordinate
   chr_start, chr_end       -- Minimum and maximum Y coordinate
);
	});
	
	my $create_table = qq{
		CREATE TABLE data (id INTEGER PRIMARY KEY,name TEXT);
	
	};
	
	$db_genes->do($create_table);
	
	
	$db_transcripts->do(qq{
		CREATE VIRTUAL TABLE positions USING rtree_i32(
   id,              -- Integer primary key
   start, end,      -- Minimum and maximum X coordinate
   chr_start, chr_end       -- Minimum and maximum Y coordinate
);
	});
	$db_transcripts->do($create_table);
	
	my $sth_transcripts_position = $db_transcripts->prepare("INSERT INTO positions VALUES (?,?,?,?,?)");
	my $sth_transcripts_name = $db_transcripts->prepare("INSERT INTO data VALUES (?,?)");
	my $sth_genes_position = $db_genes->prepare("INSERT INTO positions VALUES (?,?,?,?,?)");
	my $sth_genes_name = $db_genes->prepare("INSERT INTO data VALUES (?,?)");
	
	my $type = "gene";
 	my $z = $annot->get_like("annotations","*".$type."*");
 	my @genes_id = values %{$z};
 my $dd =0;
 my $dd_tr =0;
 foreach my $gene (@genes_id){
 	$dd ++;
 	warn $dd if $dd%100 == 0;
 	my $isncrna;
 	my $chr = $gene->{chromosome};
 	my $lat = $chr;
 	$lat = 23 if ($lat eq 'X');
 	$lat = 24 if ($lat eq 'Y');
 	$lat = 25 if ($lat eq 'MT');
	my $yy = ($lat*10)+1;
	my $y = $lat*10;
	
	foreach my $tid (@{$gene->{transcripts}}){
		$dd_tr ++;
 		my ($id,$n) = split("_",$tid);
 		my $tr = $annot->get("annotations",$id."_".$chr);
 		$isncrna = 1 if $tr->{biotype} =~ /pseudo/ || $tr->{biotype} =~ /RNA/;

 		my $limit =11;
		$limit = 0 if $isncrna;
		$limit= 0 if $chr eq 'MT';
		my $start1 = $tr->{start}-$limit;
		my $end1 = $tr->{end}+$limit;
		$sth_transcripts_position->execute($dd_tr, $start1,$end1,$y,$yy);
		$sth_transcripts_name->execute($dd_tr,$tr->{genbo_id});
		#$sth_transcripts_position->execute;
 	}
 	my $limit = 500;
		$limit = 15 if $isncrna;
		$limit= 0 if $chr eq 'MT';
	my $gstart = $gene->{start} - $limit;
	$gstart =1 if $gstart < 0;
	my $gend = $gene->{end} + $limit;
	$sth_genes_position->execute($dd, $gstart,$gend,$y,$yy);
	$sth_genes_name->execute($dd,$gene->{genbo_id});
	#$sth_genes_position->;
 }
	

}

