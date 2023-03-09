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
my ($version,$fork,$use_dir);
GetOptions(
	'version=s' => \$version,
	'fork=s' => \$fork,
	'use_dir=s' => \$use_dir,
);
$fork = 3 unless ($fork);
my $sqliteDir =  "/tmp/lmdb/$version/annotations";
$sqliteDir = $use_dir if ($use_dir);

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

my (@types) = ("transcript","protein","gene");
 #(@types) = ("gene");


 # my $z = $annot->get_like("annotations","ENSG00000169093_Y");
    #my @transcripts_id = keys %{$z};
    #warn Dumper @transcripts_id;
    #warn Dumper $z;
    #die();
    
my %genes_annot_array;
my $first = 1;
my $pm = Parallel::ForkManager->new($fork);

my $mode = "w";
$mode ="c" if $first;
my $synonym =  GenBoNoSqlAnnotation->new(dir=>$sqliteDir,mode=>$mode);
	my $annot4 =  GenBoNoSqlAnnotation->new(dir=>$sqliteDir,mode=>$mode);
	my $annot2 =  GenBoNoSqlLmdb->new(name=>"genbo",dir=>$sqliteDir,mode=>$mode,is_compress=>1);
	$annot2->create();
	#$annot2->delete("date");
	$annot2->close;
	
	$annot4->put("synonyms","datex",time);
	$annot4->close();
	my $no3 = GenBoNoSqlIntervalTree->new(dir=>$sqliteDir,mode=>"c");
	$no3->put("functional_annotations","date",time);
	$no3->close;
	 $first =1;
$pm->run_on_finish(
    sub { 
    	my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$h)=@_;
  
    	unless (defined($h) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
		#$mode ='w';	
		#my $mode = 'c' if $first;
		#$first = undef;
	  	my $annot2 =  GenBoNoSqlLmdb->new(name=>"genbo",dir=>$sqliteDir,mode=>"w",is_compress=>1);
	  	my $annot4 =  GenBoNoSqlAnnotation->new(dir=>$sqliteDir,mode=>"w");
	  	 my $no3 = GenBoNoSqlIntervalTree->new(dir=>$sqliteDir,mode=>"w");
	  	
     	my $objs = $h->{objs};
     	 warn "start save ".$h->{nb};
     	foreach my $obj (@$objs){
     		$no3->put("functional_annotations",$obj->{id},$obj->{annotations_array});
     		
     		delete $obj->{annotations_array};
     	 	$annot2->put($obj->{id},$obj);
     	 	
     	 	if (exists $obj->{total_index}){
     	 		my $id1 =  $obj->{total_index};
     	 		delete $obj->{total_index};
     	 		$annot4->put("synonyms",$id1."  ".$obj->{id},$obj->{id});
     	 	}
     	 	else {
     	 		$annot4->put("synonyms",$obj->{id},$obj->{id});
     	 		#	$annot4->put("annotations.objects",$obj->{id},$obj);
     	 	}
     	 	
     	}
     	
     	 $annot2->close();
     	 $annot4->close();
     	 warn "end save ".$h->{nb};
    }
   
   );
   
   ##################### 
   ##  Work on transcripts
   ###################
   
  
    $z = $annot->get_like("annotations","*"."transcript"."*");
 
   $annot->close();
    my @transcripts_id = keys %{$z};
     my $iter = natatime 9000, @transcripts_id;
      $nb = 0;
     
     while( my @tmp = $iter->() ){
     		$nb ++;
     	my $pid = $pm->start and next;
     		warn "start $nb/".scalar(@tmp);
     	
     		my $objs;
     		my $annot =  GenBoNoSqlAnnotation->new(dir=>$sqliteDir,mode=>"r");
     	 foreach my $id (@tmp){
     	 	
     	 		unless (exists $z->{$id}->{id}) {
 					$z->{$id}->{id} = $z->{$id}->{genbo_id};
 				}
 			
 				warn 'cpoucou $id' if $z->{$id}->{genbo_id} eq 'ENST00000628650_2'; 
     	 	my $obj = create($z->{$id},"transcript");
     	 	warn Dumper $obj if $z->{$id}->{genbo_id} eq 'ENST00000628650_2'; 
     	 	die() if $z->{$id}->{genbo_id} eq 'ENST00000628650_2';
     	 	#next;
     	 $obj->{object_type} ="transcript";
 			

 		
 		die($id) unless exists $obj->{genes_object};
 		 
 		#$obj->{genes_object}->{$z->{$id}->{gene_kyoto_id}} =undef;
 		$obj->{total_index} = $id;
 		my $exons = create_exons($obj);
 		my $introns = create_introns($obj);
 #		warn $introns:
 		my $hgene =  $annot->get("annotations",$z->{$id}->{gene_kyoto_id});
 		 annotations_array($obj	,$hgene);
 		
 		unless (keys %{$obj->{genes_object}}){
 			warn Dumper $z->{$id};
 			confess();
 		}
 		die($id) unless $z->{$id}->{gene};
 		die($id) unless $obj->{genes_object};
 		die($id) unless scalar(keys %{$obj->{genes_object}});
 		push(@$objs,$obj);
 		push(@$objs,@$exons);
 		 push(@$objs,@$introns);
     	 }
     	 my $h;
     	 $h->{objs} = $objs;
     	 
     	 $h->{nb} = $nb;
       	 warn "end $nb";
     	 $pm->finish(0,$h); 
     }
$pm->wait_all_children;




warn "end transcripts";


my $array_chr;
my $no3 = GenBoNoSqlIntervalTree->new(dir=>$sqliteDir,mode=>"w");
my $annot4 =  GenBoNoSqlAnnotation->new(dir=>$sqliteDir,mode=>"w");
 #my $no_lmdb_tree =  GenBoNoSqlLmdb->new(name=>"intervals.tree.lmdb",dir=>$sqliteDir,mode=>"c",is_compress=>1);
   ##################### 
   ##  END TRANSCRIPTS WORK ON GENES AND PROTEIN NO FORK
   ###################
  	my $annot2 =  GenBoNoSqlLmdb->new(name=>"genbo",dir=>$sqliteDir,mode=>"w",is_compress=>1);
	my (@types) = ("gene","protein");
	 #(@types) = ("protein");
	my $tree_global;

 	my $annot_read =  GenBoNoSqlLmdb->new(name=>"genbo",dir=>$sqliteDir,mode=>"r",is_compress=>1);
	
foreach my $type (@types){
	warn "START $type";
	
	$first = undef;
 my $z = $annot->get_like("annotations","*".$type."*");
 
 foreach my $id (keys %$z){
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
 	if ($type eq "gene"){
 		my $gene_array=[];
 		foreach my $t (keys %{$obj->{transcripts_object}}){
 			my $tid = $t;
 			$tid .= "_".$z->{$id}->{chromosome} unless $tid=~/_/;
 			my $tr = $annot_read->get($tid);
 			die($tid) unless $tr;
 			my $annotations_array = $no3->get_data("functional_annotations",$tid);
 			my $chr = $z->{$id}->{chromosome};
 			#warn Dumper $annotations_array;
 			foreach my $tt (@{$annotations_array}) {
 				my @a = ([$tid,$tt->[0]],$tt->[1],$tt->[2]);
 				push(@$gene_array,\@a);
 				push(@{$tree_global->{$chr}},[[$obj->id,$tid,$tt->[0]],$tt->[1],$tt->[2]]);
 			}
 			
 			
 			
 			#delete $obj->{array_annotations_tree};
 			push(@{$obj->{array_transcripts_tree}},[$tr->id,$tr->start,$tr->end]);
 			
 		}
 	#	warn $obj->{id};
 		$no3->put("functional_annotations",$obj->{id},$gene_array);
 		die() unless $obj->{strand} ;
 		my $limit = 10000;
 		 my $tree = Set::IntervalTree->new;
	 foreach my $a (@{$obj->{array_transcripts_tree}}){
	 	#foreach my $a (@{$tr->annotations_array}){
	 		$tree->insert(@$a);
	 #	}
	 }
 		

 		
 	}
 
 	
 	
 
 	$obj->{object_type} = $type;
 	
 		if ($type  eq "protein"){
 				$obj->{genes_object}->{$z->{$id}->{gene}."_".$z->{$id}->{chromosome}} = undef;
 				die($id) unless scalar(keys %{$obj->{genes_object}});
 			#	my $tid = $z->{$id}->{transcript};
 			#		$tid .= "_".$z->{$id}->{chromosome} unless $tid=~/_/;
 			#	$obj->{transcripts_object}->{$tid} = undef;
 			
 		}
 		my $tid = $obj->{id};
 		$id.=" ".$obj->{id} unless $id=~/$tid/;
 		$annot4->put("synonyms",$id,$tid);
 		$annot2->put($obj->{id},$obj);

 
 }	
$annot2->close();
if ($type eq "gene"){
	  my $no3 = GenBoNoSqlIntervalTree->new(dir=>$sqliteDir,mode=>"c");
	  foreach my $chr (keys %$tree_global){
	  	$no3->put("annotations",$chr,$tree_global->{$chr});
	  }
	  $no3->close();
	}
	if ($type eq "gene"){
	  my $no3 = GenBoNoSqlIntervalTree->new(dir=>$sqliteDir,mode=>"r");
	  my $t = time;
	  foreach my $chr (keys %$tree_global){
	  	my $tree = $no3->get("annotations",$chr);
	  }
	  warn "end sqlite => ".abs(time -$t);
			    
	  }
	  
#$annot3->close();

}		

$annot4->close();

	
	





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
	my $obj = $hashTypeObject->{$type}->new($hash);
	
	unless (exists $obj->{chromosomes_object}) {
 		$obj->{chromosomes_object}->{$hash->{chromosome}} = undef;
 		die() unless $obj->{chromosomes_object};
 		die() unless $hash->{chromosome};
 	}
 	
	return $obj;
	#warn $obj->sequence;
	
	
}

sub create_introns {
	my ($obj) = @_;
	my $out_objs;
	my $span_genomic = $obj->genomic_span();
	my $span1 =  new Set::IntSpan::Fast::XS($obj->start."-".$obj->end);
	my $span_intronic = $span1->diff($span_genomic);
	#my $span_coding = $self->span_coding();
	
	my $exons;
	
	
	my $hpos;
	my $iter     = $span_intronic->iterate_runs();
	return [] if $span_intronic->is_empty; 
	my $hreturn;
	my $pos= [];
	while ( my ( $from, $to ) = $iter->() ) {
		my $hpos;
		$hpos->{start} = $from;
		$hpos->{end} = $to;
	
		push(@$pos,$hpos);
	}
	my $num_exon = 1;
	my $start_cds =1;
	$num_exon = scalar(@$pos)+1 if  $obj->strand == -1;
	
	foreach my $hp (@$pos){
		my $from = $hp->{start};
		my $to = $hp->{end};
		
		my $hpos;
		my $ps = new Set::IntSpan::Fast::XS($from."-".$to);
		
		$hpos->{chromosomes_object}->{$obj->{chromosome}} =undef;
		$hpos->{transcripts_object}->{$obj->id}= undef;
	
		$hpos->{gstart} = $from;
		$hpos->{gend} = $to;
		$hpos->{id}= $num_exon;
		$num_exon+= $obj->strand;
		$hpos->{ext} ="intron";
		$hpos->{start} = $from;
		$hpos->{end}   = $to;
		$hpos->{name}= $hpos->{ext}."[".$hpos->{id}."-".($num_exon)."]";#.$hpos->{ext2};
		
		$hpos->{id}= $obj->id.$hpos->{ext}.$hpos->{id};#.$hpos->{ext2};
		$hreturn->{$hpos->{id}} = undef;
		$hpos->{length}   = ($to-$from)+1;
		$hpos->{strand}   = $obj->strand();
		$hpos->{intspan} = $ps;
		$hpos->{utr} = $ps->diff($obj->getSpanCoding);
		
		#$hpos->{name} .= "NC" if $hpos->{utr}->equals($hpos->{intspan});
		
		my $len = $hpos->{start}-$hpos->{end}+1;
		push( @$exons, $hpos );
	}
	if (scalar(@$exons) == 1){
		 $exons->[0]->{utr}->empty ;
		
	}
	my @temp = sort {$a->{start} <=> $b->{start}} @$exons;
	my $objs;
	foreach my $hash (@$exons){
		my $obj_intron = $hashTypeObject->{intron}->new($hash);
		push(@{$obj->{array_tree_introns}},[$obj_intron->{id},$obj_intron->start,$obj_intron->end+1]);
		$obj_intron->{chromosomes_object}->{$obj->{chromosome}} = undef;
		$obj_intron->{transcripts_object}->{$obj->id} = undef;
		$obj->{introns_object}->{$obj_intron->{id}} = undef;
		push(@$out_objs,$obj_intron);
		
	}

	return $out_objs;
	
}

sub create_exons {
	my ($obj,$annot2,$annot3) = @_;
	my $out_objs;
	my $span_genomic = $obj->genomic_span();
	#my $span_coding = $self->span_coding();
	
	my $exons;
	
	
	my $hpos;
	my $iter     = $span_genomic->iterate_runs();
	my $hreturn;
	my $pos= [];
	while ( my ( $from, $to ) = $iter->() ) {
		my $hpos;
		$hpos->{start} = $from;
		$hpos->{end} = $to;
	
		push(@$pos,$hpos);
	}

	my $num_exon = 1;
	my $start_cds =1;
	$num_exon = scalar(@$pos) if  $obj->strand == -1;

	foreach my $hp (@$pos){
		my $from = $hp->{start};
		my $to = $hp->{end};
		
		my $hpos;
		my $ps = new Set::IntSpan::Fast::XS($from."-".$to);
		#$hpos->{chromosomes_object}->{$obj->getChromosome->id} =undef;
		$hpos->{transcripts_object}->{$obj->id}= undef;
		
		$hpos->{gstart} = $from;
		$hpos->{gend} = $to;
		$hpos->{id}= $num_exon;
		$num_exon+= $obj->strand;
		$hpos->{ext} ="ex";
		$hpos->{start} = $from;
		$hpos->{end}   = $to;
		$hpos->{name}= $hpos->{ext}.$hpos->{id};#.$hpos->{ext2};
		$hpos->{id}= $obj->id.$hpos->{ext}.$hpos->{id};#.$hpos->{ext2};
	
		$hreturn->{$hpos->{id}} = undef;
		$hpos->{length}   = ($to-$from)+1;
		$hpos->{strand}   = $obj->strand();
		$hpos->{intspan} = $ps;
		$hpos->{utr} = $ps->diff($obj->getSpanCoding);
		$hpos->{name} .= "NC" if $hpos->{utr}->equals($hpos->{intspan});
		
		my $len = $hpos->{start}-$hpos->{end}+1;
		push( @$exons, $hpos );
	}
	my $objs;
	foreach my $exon (@$exons){
		
		my $exon_obj = $hashTypeObject->{exon}->new($exon);
		push(@{$obj->{array_tree_exons}},[$exon_obj->{id},$exon_obj->start,$exon_obj->end+1]);
		$exon_obj->{chromosomes_object}->{$obj->{chromosome}} = undef;
		$exon_obj->{transcripts_object}->{$obj->id} = undef;
		$obj->{exons_object}->{$exon_obj->{id}} = undef;
		push(@$out_objs,$exon_obj);
		
	}

	return $out_objs;
	
}

sub test_all {
	
	my @chrs = (1..22,"X","Y","MT");
	my $no = GenBoNoSqlIntervalTree->new(dir=>$sqliteDir,mode=>"r");
	my $no_test = GenBoNoSqlLmdb->new(name=>"genbo",dir=>$sqliteDir,mode=>"r",is_compress=>1);
	my @types = ("transcripts","genes","proteins");
	foreach my $t (@types){
		
	foreach my $chr (@chrs){
		
		my $array = $no->get_data("$t",$chr);
		foreach my $a (@$array){
			my $id =$a->[0];
			warn $id;
			my $z = $no_test->get($id);
			die($id) unless $z; 
			
		}
		
	}
	
	}
	
	
	
}



 sub transform_instpan_to_array {
 	my ($intspan,$type) = @_;
 	my $array = [];
 	my $iter = $intspan->iterate_runs();
 	 while (my ( $from, $to ) = $iter->()) {
 	 	push(@$array,["$type",$from,$to+1]);
 	 }
 	 return @$array;
 	 
 }
 


sub annotations_array {

	my ($self,$hgene) = @_;
	my $debug;
	$debug = 1   if $self->name eq "ENST00000445884";
	 #my $intspan = $self->getGene->getGenomicSpan($self->getGenomicSpan);
	 my $array =[];
	#exonic
	#$self->transform_instpan_to_interval_tree($self->getGenomicSpan(),$tree,"exonic");
	 my $span1 =  new Set::IntSpan::Fast::XS($self->start."-".$self->end);
	 
	 #splice site 
	 my $i1 = $self->getSpanSpliceSite()->intersection($span1);
	 push(@$array,transform_instpan_to_array($i1,"splice_site"));
	 
	 #essential splicing
	  $i1 = $self->getSpanEssentialSpliceSite()->intersection($span1);
	  push(@$array,transform_instpan_to_array($i1,"essential_splicing"));
	#intronic
	my $type_intronic = "intronic";
	$type_intronic = "non_coding_transcript_intron_variant" if $self->isncRNA ;
	$type_intronic = "non_coding_transcript_intron_variant" if $self->ispseudoGene;
	push(@$array,transform_instpan_to_array($self->intronic_span(),$type_intronic));
	
	
	
	if ( $self->isncRNA() ) {
			if ($self->span_mature()->is_empty){
				
				  push(@$array,transform_instpan_to_array($self->getGenomicSpan,"ncrna"));
			}	
			else {
				 push(@$array,transform_instpan_to_array($self->span_mature(),"maturemirna"));
					my $i1 = $self->getGenomicSpan->diff($self->span_mature);
					 push(@$array,transform_instpan_to_array($i1,"non_coding_transcript_exon_variant"));
			}	
		
	}
	elsif ( $self->ispseudoGene ) {
			  		push(@$array,transform_instpan_to_array($self->getGenomicSpan,"pseudogene"));
					push(@$array,transform_instpan_to_array($self->getGenomicSpan,"non_coding_transcript_exon_variant"));
	}
	else {
		push(@$array,transform_instpan_to_array($self->getSpanCoding,"coding"));
		my $i1 = $self->getGenomicSpan->diff($self->getSpanCoding);
		push(@$array,transform_instpan_to_array($i1,"utr"));
	
	}
		if ($self->{chromosome} ne "MT"){
			if ($self->strand == 1){
				push(@$array,["upstream",$self->start-5001,$self->start]);
				push(@$array,["downstream",$self->end+1,$self->end+5001]);
			} 
			else {
				push(@$array,["downstream",$self->start-5001,$self->start]);
				push(@$array,["upstream",$self->end+1,$self->end+5001]);
			}
		}
	
	warn $self->name() unless $array; 
	#die() if $debug;
	die() unless $array;
	$self->{annotations_array} = $array;
	return;
	return $array;
	
	 
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