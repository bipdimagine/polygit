#!/usr/bin/perl
use FindBin qw($Bin);
use lib $Bin;
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "/data-isilon/bipd-src/pnitschk/git2-polyweb/polygit/GenBo/lib/obj-nodb";
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
 use GenBoNoSqlRocksAnnotation;
my ($version,$fork,$use_dir,$genome_version);
GetOptions(
	'version=s' => \$version,
	'fork=s' => \$fork,
	'use_dir=s' => \$use_dir,
	'genome=s' => \$genome_version,
);

$fork = 3 unless ($fork);
my $sqliteDir =  "/tmp/lmdb/$version.$genome_version/annotations";
$sqliteDir = $use_dir if ($use_dir);
my $rocks_dir = "/data-isilon/public-data/repository/$genome_version/annotations/gencode.v$version/rocksdb/";
confess()  unless -e $rocks_dir;

my $annot =  GenBoNoSqlAnnotation->new(dir=>$sqliteDir,mode=>"r");

 my $z = $annot->get_like("annotations","*"."gene*");
 my $nb;

die(" -version= ") unless $version;
die(" -version= $sqliteDir ") unless -e $sqliteDir;
#my $in = "/data-xfs/public-data/HG19/sqlite/75/annotations/";
		my $hashTypeObject = {
		 
		 	'transcript'		=> 'GenBoTranscript',
		 	'protein'			=> 'GenBoProtein',
		 	'gene'				=> 'GenBoGene',
		 	'exon'				=> 'GenBoExon',
		 	'intron'				=> 'GenBoIntron',
 	
		};


my $rocks = GenBoNoSqlRocksAnnotation->new(dir=>$rocks_dir."/",mode=>"w",name=>"genbo");	

my $no3 = GenBoNoSqlIntervalTree->new(dir=>$rocks_dir,mode=>"w");
my $annot4 =  GenBoNoSqlAnnotation->new(dir=>$sqliteDir,mode=>"w");
 my $puid =0;
 my $guid =0; 
   ##################### 
   ##  END TRANSCRIPTS WORK ON GENES AND PROTEIN NO FORK
   ###################
	my (@types) = ("gene","protein");
	my @main_transcripts;
	
	
	
	my $tree_global;
	my $debug;
foreach my $type (@types){
	warn "START $type";
	$debug =1 if $type eq "gene";
 my $z = $annot->get_like("annotations","*".$type."*");
 
 foreach my $id (keys %$z){
 	warn $id." ".$z->{$id}->{external_name} if $debug && $id =~ /ENSG00000264462/;
 	#warn $id;
 	#next unless $id =~ /ENSG00000291050_21/;
 	#$gene->{id} = "toto";
 	unless (exists $z->{$id}->{id}) {
 		$z->{$id}->{id} = $z->{$id}->{genbo_id};
 	}
 	my $chr;
 	my $name;
 	($name,$chr) = split("_",$z->{$id}->{genbo_id});
 		
 		die("----".Dumper $z->{$id}) unless $z->{$id}->{genbo_id};
 	$z->{$id}->{name} = $name  unless $z->{$id}->{name};;
 	$z->{$id}->{chromosome} = $chr unless $z->{$id}->{chromosome};
 	my $obj;
 	eval {
 	 $obj = create($z->{$id},$type);
 	};
 	if ($@){
 		
 		warn Dumper $z->{$id};
 		warn $id;
 		die();
 	}
 	my $rocks_uid = $obj->id;
 	$puid ++;
 	if ($type eq "gene"){
 		my $gene_array=[];
 		warn $id." ".$obj->{external_name} if $debug && $id =~ /ENSG00000264462/;
 		my $count =0;
 		foreach my $t (sort {$a cmp $b} keys %{$obj->{transcripts_object}}){
 			my $tid = $t;
 			$tid .= "_".$z->{$id}->{chromosome} unless $tid=~/_/;
 			my $tr = $rocks->get("!".$tid);
 			die($tid) unless $tr;
 			my $annotations_array = $tr->functional_annotations;#$no3->get_data("functional_annotations",$tid);
 			my $chr = $z->{$id}->{chromosome};
 			foreach my $tt (@{$annotations_array}) {
 				my @a = ([$tid,$tt->[0]],$tt->[1],$tt->[2]);
 				push(@$gene_array,\@a);
 				push(@{$tree_global->{$chr}},[[$obj->id,$tid,$tt->[0]],$tt->[1],$tt->[2]]);
 			}
 			$obj->{functional_annotations} = $gene_array;
 		
 			#push(@{$obj->{array_transcripts_tree}},[$tr->id,$tr->start,$tr->end]);
 			
 		}
		$rocks->put_raw($obj->{external_name},$rocks_uid);
 		my ($a,$b) = split("_",$rocks_uid);
     	$rocks->put_raw("g:".$b.":".$puid,$rocks_uid);
     	$rocks->put_raw($a,$rocks_uid);
     	
     	$obj->{functional_annotations} = $gene_array;
     	 
 		#$no3->put("functional_annotations",$obj->{id},$gene_array);
 		die() unless $obj->{strand} ;
 		my $limit = 10000;
 
		my $toto = return_main_transcript($rocks,$obj);
		push(@main_transcripts,@$toto);
 		
 	}
 
 	
 	
 
 		$obj->{object_type} = $type;
 	
 		if ($type  eq "protein"){
 				$obj->{genes_object}->{$z->{$id}->{gene}."_".$z->{$id}->{chromosome}} = undef;
 				die($id) unless scalar(keys %{$obj->{genes_object}});
 				my ($a,$b) = split("_",$rocks_uid);
     	 		$rocks->put_raw("p:".$b.":".$puid++,$rocks_uid);
 				$rocks->put_raw($a,$rocks_uid);
 		}
 		my $tid = $obj->{id};
 		$id.=" ".$obj->{id} unless $id=~/$tid/;
 		
 		
 		
 		#foreach my $a (split(" ",$obj->{id})){
 		#	 $rocks->put_raw($a,$rocks_uid);
 			
 		#}
 			
		$rocks->put("!".$rocks_uid,$obj);
 }#end objects	

if ($type eq "gene"){
	  my $no3 = GenBoNoSqlIntervalTree->new(dir=>$rocks_dir,mode=>"c");
	  foreach my $chr (keys %$tree_global) {
	  	$no3->put("annotations",$chr,$tree_global->{$chr});
	  }
	  $no3->close();
	}
	if ($type eq "gene"){
	  my $no3 = GenBoNoSqlIntervalTree->new(dir=>$rocks_dir,mode=>"r");
	  my $t = time;
	  foreach my $chr (keys %$tree_global){
	  	my $tree = $no3->get("annotations",$chr);
	  }
	  warn "end sqlite => ".abs(time -$t);
			    
	  }
	  
#$annot3->close();

}		

$annot4->close();

warn "put main";

foreach my $tid (@main_transcripts){
	my $t = $rocks->genbo($tid);
	die() unless $t;
	$t->{isMain} = 1;
	$rocks->put("!".$tid,$t);
	
}
$rocks->start_iter("g:");

	

my $rocks = GenBoNoSqlRocksAnnotation->new(dir=>$rocks_dir."/",mode=>"r",name=>"genbo");	

my $xx =0;
while (my $id = $rocks->next("g:")){
	my $oi = $rocks->genbo($id);
	$xx++;
	warn Dumper $oi->functional_annotations;
	last if $xx > 10;
}
$rocks->close();

sub create {
	my ($hash,$type) = @_;
	#return  unless $hash->{name};
		my $chr;
 	my $name;
 	($name,$chr) = split("_",$hash->{genbo_id});
 		#
 		die("----".Dumper $hash) unless $hash->{genbo_id};
 		$hash->{name} = $name  unless $hash->{name};;
 		$hash->{chromosome} = $chr unless $hash->{chromosome};
	
	my $ts = "_".$hash->{chromosome};
	unless ($hash->{id}=~/$ts/){
		 $hash->{id} =~s/_//;
		 $hash->{id} = $hash->{id}."_". $hash->{chromosome};
		  $hash->{genbo_id} = $hash->{id}."_". $hash->{chromosome};
	}
	$hash->{remap_status} = 0 unless exists $hash->{remap_status} ;
	my $obj = $hashTypeObject->{$type}->new($hash);
	
	unless (exists $obj->{chromosomes_object}) {
 		$obj->{chromosomes_object}->{$hash->{chromosome}} = undef;
 		die() unless $obj->{chromosomes_object};
 		die() unless $hash->{chromosome};
 	}
 	
	return $obj;
	#warn $obj->sequence;
	
	
}









sub return_main_transcript {
	my ($rocks,$gene) = @_; 
	my $htr;
	my @kh = keys %{$gene->{transcripts_object}};
	if (scalar (@kh) == 1 ){
		return [$kh[0]];
	}
	
 	foreach my $tid (sort {$a cmp $b} @kh ){
 		$htr->{$tid} = $rocks->genbo($tid);
		 		
 	}
 	
 	my $res = 	mane_ccds_principal([values %$htr]);
 	$res = alternative([values %$htr]) unless $res;
 	$res = last_chance([values %$htr]) unless $res;
 	$res = longest([values %$htr]) unless $res;
 	my $tmain = [];
 	foreach my $k (keys %$res) {
		my $id = $res->{$k}->[0];
		push(@$tmain,$id);
	}
	$gene->{main_transcripts} = $tmain 	
 	
}

sub mane_ccds_principal {
	my ($transcripts) = @_;
	my $tr_mains = {};
	foreach my $tr (@$transcripts){
		my $span = $tr->getSpanCoding;
 		$span = $tr->getGenomicSpan if $span->is_empty;
 		my $tl = $tr->transcript_support_level;
 		die() if $tl eq "NA";
		die()  if $tl eq "Missing";
 		my $tags = join(";",keys %{$tr->{tag}});
 		if ($tags =~ /MANE/){
 			push(@{$tr_mains->{$span->as_string."*"}},$tr->id);
 			next;
 		}
 		 if ( $tags =~ /CCDS/ or $tags =~ /appris_principal/ or $tl == 1) {
 		 	
 		 	push(@{$tr_mains->{$span->as_string}},$tr->id);
 		 }
		}
	return $tr_mains;
}

sub alternative {
	my ($transcripts) = @_;
	my $tr_mains = {};
	foreach my $tr (@$transcripts){
 		
 		my $span = $tr->getSpanCoding;
 		$span = $tr->getGenomicSpan if $span->is_empty;
 		my $tl = $tr->transcript_support_level;
 		$tl = 6 if $tl eq "NA";
		$tl = 6 if $tl eq "Missing";
 		my $tags = join(";",keys %{$tr->{tag}});
 		 if ( $tags =~ /appris_alternative/ and $tl <= 3){
 		 	push(@{$tr_mains->{$span->as_string}},$tr->id);
 		 
		}
 		}
	return $tr_mains;
}

sub last_chance {
	my ($transcripts) = @_;
	my $tr_mains = {};

	my @t  =sort {$a->transcript_support_level <=> $b->transcript_support_level} @{$transcripts} ;
 	my $min = $t[0]->transcript_support_level;
	
	foreach my $tr (sort {$a->transcript_support_level <=> $b->transcript_support_level} @$transcripts){
		last if $tr->transcript_support_level > $min;
 		next if ( $tr->{biotype} =~/processed_pseudogen/);
 		my $span = $tr->getSpanCoding;
 		$span = $tr->getGenomicSpan if $span->is_empty;
		push(@{$tr_mains->{$span->as_string}},$tr->id);
 		
 		}
	return $tr_mains;
}

sub longest {
	my ($transcripts) = @_;
	my $tr_mains = {};

	my @t  =sort {$a->length <=> $b->length} @{$transcripts} ;
	my $tr = $t[-1];
	my $span = $tr->getSpanCoding;
 	$span = $tr->getGenomicSpan if $span->is_empty;
	push(@{$tr_mains->{$span->as_string}},$tr->id);
	return $tr_mains;
}
#sub get_void_object {
#	my ($self,$type) = @_;
#	return dclone ($self->{void}->{$type}) if exists $self->{void}->{$type};
#	my $typeObj = $hashTypeObject->{$type};
#	
#	$self->{void}->{$type} =$typeObj->new(id=>"titi");
#		return dclone ($self->{void}->{$type} );
#}
#
#
#sub createObject {
#	my ($self, $type, $hash) = @_;
#		my $hashTypeObject = $self->hashTypeObject();
#	confess("\n\nERROR: No type defined to create object !! die...\n\n") unless exists $hashTypeObject->{$type};
#	$hash->{project} = $self;
#
#	my $typeObj = $hashTypeObject->{$type};
#
#	my $z = 0;
#	$z = time;
#	my $object;
#	if ($type eq "variations" or $type eq "deletions" or $type eq "insertions" or $type eq "large_duplication" or $type eq "large_deletions") {
#		$object = $self->get_void_object($type);
#		foreach my $k (keys %$hash){
#			$object->{$k} = $hash->{$k};
#		}
#	  }
#	else {
#		$object = $typeObj -> new ($hash);
#	}
#	  
##	#my $zz = dclone($object);
#	 if ($type eq "variations"){
#	 	my $a =abs(time -$z);
#	 	$self->time_test($self->time_test+$a);
#	 }
#	 
#	#$object->{project} = undef;
#	$self->{objects}->{$type}->{$hash->{id}} = $object;
#	#$object->{project} = $self->encoder->encode($object);
#}