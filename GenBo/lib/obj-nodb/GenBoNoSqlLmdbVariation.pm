package GenBoNoSqlLmdbVariation;
use strict;
use warnings;
use DBD::SQLite;
use Moose;
use MooseX::Method::Signatures;
use Data::Dumper;
use Set::IntSpan::Fast::XS;
use Module::Load;
#use Compress::Zstd;
use Compress::Snappy;
use Storable qw/thaw freeze dclone/;
#use Gzip::Faster;
extends "GenBoNoSqlLmdb";
use POSIX;

 my $himpact_sorted = {
	"high" => "4",
	"moderate" =>"3",
	"low" =>"1",
};

sub nb_keys {
	my ($self) = @_;
	my $db_name = $self->name();
	return 0 unless -e $self->filename();
	my $cursor2 = $self->lmdb($db_name."_index")->stat();
	return $cursor2->{entries};
}

has array_for_bit => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $for_bit = $self->read_string("for_bit");
		unless ($for_bit){
				confess() if $self->mode() ne "c";
		 $for_bit = "isPanel,isLargeDeletion,isLargeDuplication,isSVDeletion,isMnp,isProtein,isSVDuplication,isRegulatoryRegion,isGene,isInsertion,isComplex,isVariant,isVariation,isTranscript,isPublic,svcaller,isCnv";
		 $self->save_string("for_bit",$for_bit);
		}
		my @toto = split(",",$for_bit);
		return \@toto;
	},
);
has array_seq_by_patient => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $array_seq_by_patient = $self->read_string("array_seq_by_patient");
		unless ($array_seq_by_patient){
			confess() if $self->mode() ne "c";
			$array_seq_by_patient = "ratio:w,depth:w,nalt:w,nref:w,genotype:w,norm_depth:w,dude:w,sr1:w,sr0:w,pr0:w,pr1:w,sr0_raw:w,:w,sr1_rawf:w,sr1_rawr:w,sr_quality_start:w,sr_quality_end:w,event_quality:w,text:s";
		#	$array_seq_by_patient = "ratio:w,depth:w,nalt:w,nref:w,genotype:w,infos_text:s";
			$self->save_string("array_seq_by_patient",$array_seq_by_patient);
		}
		my @toto =  split(",",$array_seq_by_patient);
		return \@toto;
	},
);

has array_for_string => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $for_string = $self->read_string("for_string");
		
		unless ($for_string){
			confess() if $self->mode() ne "c";
		 	$for_string =  "name,gnomad_id,id,vcf_sequence,type_object,type_public_db,structuralTypeObject,ref_allele,min_pop_name,check_id,vcf_id,max_pop_name,var_allele,structuralType,validation_method,rs_name,sequence";
		 	 $self->save_string("for_string",$for_string);
		}
		my @toto = split(",",$for_string);
		return \@toto;
	},
);

has array_for_string2 => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $for_string = $self->read_string("for_string");
		
		unless ($for_string){
			confess() if $self->mode() ne "c";
		 	$for_string =  "name,gnomad_id,id,vcf_sequence,ref_allele,check_id,vcf_id,var_allele,rs_name,sequence";
		 	 $self->save_string("for_string",$for_string);
		}
		my @toto = split(",",$for_string);
		return \@toto;
	},
);
has array_for_index_string => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $for_string = $self->read_string("for_string");
		
		unless ($for_string){
			confess() if $self->mode() ne "c";
		 	$for_string =  "type_object,type_public_db,structuralTypeObject,min_pop_name,max_pop_name,structuralType,validation_method";
		 	 $self->save_string("for_string",$for_string);
		}
		my @toto = split(",",$for_string);
		return \@toto;
	},
);

has array_for_transcript_name => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $for_string = $self->read_string("for_transcript_name");
		
		unless ($for_string){
			confess() if $self->mode() ne "c";
		 	$for_string =  "enst,nm,ccds,appris,consequence,impact_score_text";
		 	 $self->save_string("for_transcript_name",$for_string);
		}
		my @toto = split(",",$for_string);
		return \@toto;
	},
);

has array_for_transcript_keys => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $for_string = $self->read_string("for_transcript_keys");
		unless ($for_string){
			confess() if $self->mode() ne "c";
			 $for_string = "enst:ti,nm:ti,ccds:ti,appris:ti,consequence:ti,impact_score_text:ti,st_name,prot_pos,codons,codons_AA,sift:wi,polyphen:wi,exon:ti,nomenclature,impact_score,main";
		 	 $self->save_string("for_transcript_keys",$for_string);
		}
		
		my @toto = split(",",$for_string);
		return \@toto;
	},
);
has array_for_transcript_keys_by_type => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my @toto;
		foreach my $t (@{$self->array_for_transcript_keys}){
			my ($k,$v) = split(":",$t);
			push(@toto,[$k,$v]);
		}
		return \@toto;
	}
);
has array_for_integer => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $for_integer = $self->read_string("for_integer");
		unless ($for_integer){
			confess() if $self->mode() ne "c";
		 	$for_integer = "start,end,strand,vcf_position,length,gaMall,ganall,gacall,ghoall,sum_mask_coding,cadd_score,ngs_score";
		 	 $self->save_string("for_integer",$for_integer);
		}
		my @toto = split(",",$for_integer);
		return \@toto;
	},
);
has array_for_float => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $for_float = $self->read_string("for_float");
		unless ($for_float){
			confess() if $self->mode() ne "c";
			 $for_float = "frequency,max_pop_freq,min_pop_freq,ncboost_score,dbscsnv_ada,dbscsnv_rf,revel_score";
			 $self->save_string("for_float",$for_float);
		}
		my @toto =  split(",",$for_float);
		return \@toto;
	},
);
has array_for_annex => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		
		my $for_annex =  $self->read_string("for_annex");
		unless ($for_annex){
			confess() if $self->mode() ne "c";
		 	$for_annex ="nb_all_ref:w,ho:w,he:w,nb_all_mut:w,score:w,dp:w,is_cas_1_2:w,nb_all_other_mut:w,sr:a,pr:a,method:a";
		 	 $self->save_string("for_annex",$for_annex);
		}
		my @toto =  split(",",$for_annex);
		return \@toto;
	},
);


has array_for_annotation => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $for_annot =  $self->read_string("for_annot");
		unless ($for_annot){
		confess() if $self->mode() ne "c";
		 $for_annot ="prot_position:w,transcript_position:w,pos_codon:w,orf_position:w,orf_end:w,aa_mut:a,seq_mut:a,codon:a,aa:a,codon_mut:a,seq_orf:a";
		 $self->save_string("for_annot",$for_annot);
		}
		my @toto =  split(",",$for_annot);
		return \@toto;
	},
);

has array_for_dejavu => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $for_annot =  $self->read_string("for_dejavu");
		unless ($for_annot){
		confess() if $self->mode() ne "c";
		 $for_annot ="other_project,other_patients,other_patients_ho,similar_projects,similar_patients,similar_patients_ho,this_run_patients";
		 $self->save_string("for_dejavu",$for_annot);
		}
		my @toto =  split(",",$for_annot);
		return \@toto;
	},
);

#has array_for_gnomad => (
#	is      => 'rw',
#	lazy    => 1,
#	default => sub {
#		my $self = shift;
#		my $for_annot =  $self->read_string("for_gnomad");
#		unless ($for_annot){
#		confess() if $self->mode() ne "c";
#		 $for_annot =",gnomad_ac:w,gnomad_an:w,gnomad_ho:w,gnomad_ho_male:w,gnomad_min_pop:w,gnomad_max_pop:w,gnomad_min_pop_name:ti,gnomad_max_pop_name:ti";
#		 $self->save_string("for_dejavu",$for_annot);
#		}
#		my @toto =  split(",",$for_annot);
#		return \@toto;
#	},
#);

	
	


has dictionary => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
			if (@{$self->array_dictionary}){
			my $hash ={};
			for (my $i=0;$i<@{$self->array_dictionary};$i++){
				$hash->{$self->array_dictionary->[$i]} = $i;
			}
			return $hash;
		}
		return {};
	},
);
has last_index => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return 0;
	},
);
has array_dictionary => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		
		return $self->load_dictionary();
	},
);
sub get_object_index {
	my ($self,$o) = @_;
	#$o = "undef" unless defined $o;
	return $self->dictionary->{$o} if exists  $self->dictionary->{$o};
	 $self->dictionary->{$o} = $self->last_index;
	 $self->last_index($self->last_index + 1);
	 return   $self->dictionary->{$o};
	
}
my $last_index = 0;



 sub encode_annex {
 	my ($self,$vannex) = @_;

 	my $st;
	my $st_annex="";
	my $nb =0;
	my $current = undef;
 	
 	foreach my $at (@{$self->array_for_annex}){
				my ($a,$t) = split(":",$at);
				$current = $t unless $current;
				my $value = undef;
				 
				 $value =  $vannex->{$a} if exists  $vannex->{$a};
				unless ( defined $value){
						$t = "P";
						$value =  undef;
					}
				elsif ($t eq "a"){
					$t = "a".length($value);
				}
			  if ($t eq "w") {
						$value = int($value);
			}
				$st_annex .= $t.""; 
				$st .= pack("$t",$value);
				delete $vannex->{$a};
			}
 		return ($st_annex,$st);
 }
 
 sub return_code_numeric {
 	my ($self,$v) = @_;
 	if($v>=0){
 	return "C" if ( $v < 255);
 	return "I" if ( $v < 65535);
 	return"w";
 	}
 	else {
 	return "c" if ($v>-127);
 	return "i" if ( $v<-32767);
 	return "l";# if ( $v<-32767);
 	}
 }

sub encode_calling_information {
	my ($self,$variation,$debug) = @_;
	 	my @patients = keys %{$variation->{seq}};
	 	
	 	foreach my $pid (@patients){
	 		my @code;
	 		my $string_seq ="";
	 			foreach my $a (@{$self->array_seq_by_patient}){
	 				my ($k,$v) = split(":",$a);
	 				unless  (exists $variation->{seq}->{$pid}->{$k}) {
	 					$string_seq .= pack("A","!");
	 					 push(@code,"A");
	 					 next;
	 				}
	 				if ($variation->{seq}->{$pid}->{$k} eq "-") {
	 					$string_seq .= pack("A","-");
	 					 push(@code,"A");
	 					 next;
	 				}
	 				
	 				my $value = $variation->{seq}->{$pid}->{$k};
	 				if ($v eq "w"){
	 					my $cc = $self->return_code_numeric($value);
	 					$string_seq .= pack("$cc",$value);
	 					 push(@code,"$cc");
	 					 next;
	 				}
	 				elsif ($v eq "s"){
	 					my $cc = "A".length($value);
	 					$string_seq .= pack("$cc",$value);
	 					 push(@code,"$cc");
	 					 
	 				}
	 				else{confess()}
	 			}
	 			#push(@{$variation->{s1}},[$pid,$string_seq] );
	 			push(@{$variation->{s1}},[$pid,join("",@code),$string_seq] );
	 	}
	 #	die() if $debug;
}

sub decode_calling_information {
	my ($self,$variation) = @_;
		my @tab = @{$variation->{s1}};
		foreach my $line ( @{$variation->{s1}}){
			my $pid = $line->[0];
			my @values = unpack($line->[1],$line->[2]);
			die() if scalar(@values) ne scalar(@{$self->array_seq_by_patient});
			
			for (my $i=0 ;$i<  @values;$i++) {
				next if  $values[$i] eq "!"; 
				my ($k,$v) = split(":",$self->array_seq_by_patient->[$i]);
				$variation->{seq}->{$pid}->{$k} = $values[$i];
			}
			
		}
		delete $variation->{s1};
#		use Test::More tests => 1000;
#		warn "TEST 4 ";
	#	warn Dumper  $variation->{seq1};
	#	warn Dumper  $variation->{seq};
	#	die();
}



sub decode_annotation_polyviewer{
	my ($self,$age) = @_;
	my $annotation_polyviewer = {};
	return $annotation_polyviewer unless $age;
	
	my $nb = scalar(@{$self->array_for_transcript_name});
	foreach my $hgene (@{$age}){
		my $gid = unpack("w",$hgene->{id});
		my $gn = $self->array_dictionary->[$gid];
		my @v =  unpack("c W",$hgene->{sai});
		$annotation_polyviewer->{$gn}->{spliceAI_score} = $v[0]/100;
		$annotation_polyviewer->{$gn}->{spliceAI_cat} = $self->array_dictionary->[$v[1]];
		delete $hgene->{sai};
		my $string = delete $hgene->{tr};
		my @array = split("!",$string);
		my @res;
		foreach my $st2 (@array){
			my @array2 = split(";",$st2);
			die(scalar(@array2)." ".scalar(@{$self->array_for_transcript_keys})) if (scalar(@array2)) ne scalar(@{$self->array_for_transcript_keys});
			my $tr = {};
			for (my $i=0;$i<@array2;$i++){
				my $k = $self->array_for_transcript_keys_by_type->[$i]->[0];
				my $si = $self->array_for_transcript_keys_by_type->[$i]->[1];
				my $v = $array2[$i];
				
				next if ($v eq -2);
				if ($v eq -1){
					$tr->{$k}  = undef;
					next;
				}
				unless ($si){
						$tr->{$k} = $v;
				}
				elsif ($si eq "ti"){
					$tr->{$k}  = $self->array_dictionary->[$v];
					
				}
				elsif ($si eq "wi"){
					#my $id1 = $self->get_object_index($value);
					$tr->{$k} = $v;
				}
				else {
					confess();
				}
				
			}
			push(@res,$tr);
		}
		$annotation_polyviewer->{$gn}->{trs} = \@res ;

	}
	return $annotation_polyviewer;
 
}	
my $first;
sub encode {
		my ($self,$variation) = @_;
		my $debug;
#		$debug = 1 if "Y_10054455_G_GGTGATGTGCTCATTCATCTCACAGATTTGAAATGTTCTTTTCATTGACCAGTTTGGATAGAGTCCTTTTGTAGGATTTGTTTTGCGATATTTGTGAGCCCTTTGAAGCCTATGGTGAAAAAAGAAATATCTTCACATAAAAACTAGACAGAAGCTTTCTGAGAAACTTCGTTGTGATGTATGCATTCATCCCAAAGAGTTGAGCCTGTCTTTGGATTGAGCAATTTTGAAAGAGTTCTTTTGTAGAATGTTCAAAGGGATATTTGGGATCCTTTTTTGGCCTATGGTGAAAAAGGAAATGTCTTCATATAAAAACTAGACAGAAACATTCTGAGAAACTTCTT";
		my $hv;
		$first ++;
		#jete check_id
		#my $for_bit = "isPanel,isLargeDeletion,isLargeDuplication,isSVDeletion,isMnp,isProtein,isSVDuplication,isRegulatoryRegion,isGene,isInsertion,isComplex,isVariant,isVariation,isTranscript,isPublic";
		$variation->isInsertion;
		$variation->isPublic;
		$variation->length;
		my $vid = $variation->{id};
		my $null;
		my $nb = 0;
		my $st_byte="";
		### BIT
		foreach my $id (@{$self->array_for_bit}){
			unless (defined  $variation->{$id}){
				die($id);
			}
			
			#$variation->{bit} .= pack("C",$variation->{$id});
			$variation->{values}->{$id} = $variation->{$id};
			$st_byte.=$variation->{$id};
			# $st_byte = pack ("b12",'000101011100');
			$nb++;
			delete $variation->{$id};
		}
		$variation->{bit} .= pack("b".$nb,$st_byte);
		$variation->{st_bit} = "b".$nb;
		
		
		### INTEGER
		my $st_integer="";
		$nb =0;
		foreach my $id (@{$self->array_for_integer}) {
			if (!(defined $variation->{$id}) or $variation->{$id} eq '-' or $variation->{$id} eq "" or  $variation->{$id} == -1) {
				if($nb>0){
					$st_integer .="w$nb ";
				}
				$st_integer .="A ";
				my $c ="U";
			
				$c = "-" if $variation->{$id} && $variation->{$id} eq '-';
				$variation->{$id} = $c;
				$variation->{values}->{$id}  = $c unless $c eq "U";
				$variation->{int} .= pack("A",$c);
				delete $variation->{$id} ;
				$nb=0;
				$null =1;
				next;
			}
			$variation->{values}->{$id} = $variation->{$id};
			warn $id." ".$variation->{$id} if $variation->{$id} < 0 ;
			$variation->{int} .= pack("w",int($variation->{$id}));
			delete $variation->{$id};
			$nb++;
		}
		if($nb>0){
			$st_integer .="w$nb";
		}
		
		$variation->{st_int} =$st_integer;
		
		#FLOAT 
		my $st_float="";
		foreach my $id (@{$self->array_for_float}){
			if (!(defined $variation->{$id})){
				$variation->{values}->{$id} = $variation->{$id};
				$st_float .="P";
				$variation->{float} .= pack("P",undef);
				delete $variation->{$id};
				next;
			}
			
			if ( $variation->{$id} eq "-") {
				$variation->{values}->{$id} = $variation->{$id};
				$st_float .="A";
				$variation->{float} .= pack("A","-");
				delete $variation->{$id};
				next;
			}
			#$st.="f";
			$variation->{values}->{$id} = $variation->{$id};
			#die($variation->{$id})  if $variation->{values}->{$id} == 0 && $variation->{$id} > 0;
			$st_float .="F";
			$variation->{float} .= pack("F",$variation->{$id});
			delete $variation->{$id};
		}
		$variation->{st_float} =$st_float;
		
		###########
		#STRING INDEX ###
		###########
		my @int1;
		foreach my $id (@{$self->array_for_index_string}){
		#	next if $id =~ /id/;
			#warn $id." ".$variation->{$id}."\n";
			
			 $variation->{$id} = "" unless $variation->{$id};
			 #warn $variation->{$id}." ".$sid;
		# my $sid = $self->get_object_index($variation->{$id});
			#delete $variation->{$id};
			#$variation->{f} .= pack("s",$variation->{$id});
			#push(@{$variation->{f}} ,$sid);
		}
		#$variation->{fsti} = pack("w".scalar(@int1),@int1);
		###########
		#STRING ###
		###########
		my $t;
		foreach my $id (@{$self->array_for_string}){
		#next unless exists $variation->{$id};
		#	next if $id =~ /id/;
			#warn $id." ".$variation->{$id}."\n";
			 $variation->{$id} = "" unless $variation->{$id};
			push(@{$variation->{fs}},$variation->{$id});
			
			$variation->{values}->{$id} = $variation->{$id}."";
			delete $variation->{$id}
		}
		
		$variation->{fst} = join("!",@{$variation->{fs}});
		delete $variation->{fs};
		
	#foreach my $id (@{$self->array_for_string2}){
	#				delete $variation->{$id}
	#	}
		
	
		                                                                                                   
		#die();
	##############	
	#	seq deatil
	#############	
	#$variation->{values}->{seq2} = dclone $variation->{seq};
	#$self->encode_calling_information($variation,$debug);

	my $xst;

	delete $variation->{annex}  unless keys %{$variation->{annex}};	
	delete $variation->{annex};
	#object id 
	############
	# DEJAVU 
	##############
	my $string;
	foreach my $a (@{$self->array_for_dejavu}){
		$string .= pack("w",$variation->{value}->{$a});
		delete $variation->{value}->{$a};
	} 
	$variation->{st_dejavu} = $string;
	delete $variation->{dejaVuInfosForDiag2};
	my $hash_oid;
	my $array =[];
	my $toto;
foreach my $o (keys %{$variation->{annotation}}) {
		my $id = $self->get_object_index($o);
		$toto = 1 if (exists $variation->{annotation}->{$o}->{coding}->{sequences});
		delete $variation->{annotation}->{$o}->{coding} unless (exists $variation->{annotation}->{$o}->{coding}->{sequences});
		$variation->{annotation}->{$id} = delete $variation->{annotation}->{$o};
}

	$variation->{ano} = join(";",@$array);
$debug =1 if $vid eq "Y_2718882_C_CT";	
if (exists $variation->{annotation_polyviewer})	{
my @genes;
foreach my $o (keys %{$variation->{annotation_polyviewer}}) {
	my $hgenes ={};
	my $id = $self->get_object_index($o);
	$hgenes->{id} = pack("w",$id);
	$variation->{annotation_polyviewer}->{$id} = delete $variation->{annotation_polyviewer}->{$o};
	my $nb  = int($variation->{annotation_polyviewer}->{$id}->{spliceAI_score} *100);
	$hgenes->{sai} = pack("c",$nb);
	my $tid= $self->get_object_index($variation->{annotation_polyviewer}->{$id}->{spliceAI_cat});
	$hgenes->{sai} .= pack("w",$tid);
	my $array =  delete  $variation->{annotation_polyviewer}->{$id}->{trs};
	my $nx = scalar(@{$self->array_for_transcript_name});
	my $key = "st_name;prot_pos;codons;codons_AA;sift;polyphen;exon;nomenclature";
	my @akey = split(";",$key);
	my $string_annot ="";
	foreach my $tr (@$array){
			my $pack_string ="";
			my @tt;
			foreach my $ki  (@{$self->array_for_transcript_keys_by_type}){
				my $k = $ki->[0];
				my $si = $ki->[1];
				unless (exists $tr->{$k}){
					push(@tt, -2);
					next;
				}
				unless ($tr->{$k}){
					push(@tt, -1);
					next;
				}
				my $value = delete $tr->{$k};
				unless ($si){
						push(@tt, $value);
				}
				elsif ($si eq "ti"){
					my $id1 = $self->get_object_index($value);
					push(@tt, $id1);
				}
				elsif ($si eq "wi"){
					push(@tt, $value);
				}
				else {
					confess();
				}
			}
			$string_annot .= join(";",@tt)."!";
			}
		
	$hgenes->{tr} = $string_annot;
	push(@{$variation->{age}},$hgenes);  
}

$variation->{annotation_polyviewer} ={};
delete $variation->{annotation_polyviewer};
}	
delete $variation->{transcripts_object};
################
# MAXAI 
##############
	
	# foreach my $g  (keys %{$variation->{AI_max_cat}}){
	 #	my $zid = $self->get_object_index($g);
	 #	$variation->{AI}->{$zid} = $variation->{AI_max_cat}->{$g};
	# }
	
	if (exists $variation->{vannot}){
		my $i =0;
		$variation->{transcripts_object} = {};
		map {$variation->{transcripts_object}->{$_} = 1} keys %{$variation->{vannot}};
	}
	else {
		$variation->{transcripts_object} = {};
	}
	
	my $oo =[];
#	$variation->{values}->{transcripts_object} = dclone $variation->{transcripts_object};
	foreach my $t (keys %{$variation->{transcripts_object}}){
		push(@$oo, $self->get_object_index($t));
	}
	
	$variation->{tro} = join(";",@$oo);
	
	 $oo =[];
	
	 my $a =  dclone $variation->{genes_object};
	 delete $variation->{genes_object};
	 $variation->{genes_object} = {};
	foreach my $t (@{$a}){
		$variation->{genes_object}->{$t} =1;
		push(@$oo, $self->get_object_index($t));
	}
	$variation->{geo} = join(";",@$oo);	
	# $variation->{values}->{genes_object} = dclone $variation->{genes_object};
	 
		 $oo =[];
	foreach my $t (keys %{$variation->{references_object}}){
		$variation->{references_object}->{$t} = 1;
		push(@$oo, $self->get_object_index($t));
	}
	$variation->{refo} = join(";",@$oo);	
#	 $variation->{values}->{references_object} = dclone $variation->{references_object};
			 $oo =[];
	
	foreach my $t (keys %{$variation->{patients_object}}){
		$variation->{patients_object}->{$t} =1;
		push(@$oo, $self->get_object_index($t));
	}
	$variation->{values}->{patients_object} = dclone $variation->{patients_object};		 
	$variation->{pato} = join(";",@$oo);	
	delete $variation->{transcripts_object};
	delete $variation->{genes_object};
	delete $variation->{patients_object};
	delete $variation->{references_object};
	
	
	delete $variation->{vannot};
	delete $variation->{position};
	delete $variation->{intspan};
	delete $variation->{AI};
	delete $variation->{position};
	delete $variation->{patients_details};
	delete $variation->{heho_string};
	delete $variation->{gnomad};
	delete $variation->{AI_max};

	delete $variation->{seq_details};
	$variation->{values} = {};
	delete $variation->{values};
	delete $variation->{AI_max_cat};
	delete $variation->{validation_ngs};
	delete $variation->{line_infos};
	delete $variation->{chromosome_object};
	
	## 
	#test 
	delete $variation->{impact};
	delete $variation->{nomenclature};

	delete $variation->{impact};

	delete $variation->{seq};
	delete $variation->{value};
	delete $variation->{in_this_run_patients};
	delete $variation->{validation_sanger};

	#delete $variation->{fst};
	#die();
	return  $self->_compress(freeze ({data=>$variation}));
}

sub save_dictionary{
	my ($self) = @_;
	my $db_name = $self->name();
	my $array;
	my $i=0;
	foreach my $k (sort {$self->dictionary->{$a} <=> $self->dictionary->{$b}} keys %{$self->dictionary} ){
		die($k." ".$i." ".$self->dictionary->{$k}) if $i ne $self->dictionary->{$k};
		push(@{$self->array_dictionary},$k);
		$i++;
	}
	$self->lmdb($db_name)->put("array_dictionary",   $self->_compress(freeze ($self->array_dictionary) ));
	return  1;
} 

sub load_dictionary {
	my ($self) = @_;
	my $db_name = $self->name();
	my $code = $self->lmdb($db_name)->get("array_dictionary");
	return [] unless $code;
	my $array =  thaw($self->_uncompress($code));
	
	return  $array if $array;
	return [];
} 
sub save_string {
	my ($self,$id,$string) = @_;
	my $db_name = $self->name();
	$self->lmdb($db_name)->put($id,   $self->_compress($string));
	return 1;
} 

sub read_string {
	my ($self,$id) = @_;
	my $db_name = $self->name();
	my $code = $self->lmdb($db_name)->get($id);
	return unless $code;
	my $string =  $self->_uncompress($code);
	return $string;
} 

sub decode_annex {
	my ($self,$st_annex,$st,$hash,) = @_;
	my @values = unpack($st_annex,$st);
		my $i =0;
		confess(Dumper (@values) ." ".$st_annex." ".$st." ") if scalar(@values) ne scalar @{$self->array_for_annex};
		
				foreach my $av (@{$self->array_for_annex}) {
						my ($a,$n) = split(":",$av);
						unless ( defined $values[$i]){
							$i++;
							next;
						}
						$hash->{$a} = $values[$i];
					$i++;
				}
				
}



sub decode {
	my ($self,$code) = @_;
	return undef unless $code;
	my $obj = thaw ($self->_uncompress($code));
	my $variation = $obj->{data};	
	###############
	#objects 
	###############
	$variation->{transcripts_object} ={};
	foreach my $i (split(";",$variation->{tro}) ){
		my $v = $self->array_dictionary->[$i];
		$variation->{transcripts_object}->{$v} =1;
	}
	
	delete $variation->{tro};
		if (exists $variation->{values}){
		die() unless is_deeply($variation->{values}->{transcripts_object}, $variation->{transcripts_object}, 'data structures should be the same');
	}
	$variation->{genes_object} ={};
	foreach my $i (split(";",$variation->{geo})){
			my $v = $self->array_dictionary->[$i];
			$variation->{genes_object}->{$v} =1;
	}
	delete  $variation->{geo};
		$variation->{references_object} ={};
	foreach my $i (split(";",$variation->{refo})){
			my $v = $self->array_dictionary->[$i];
			$variation->{references_object}->{$v} =1;
	}
	delete  $variation->{refo};
	foreach my $i (split(";",$variation->{pato})){
		my $v = $self->array_dictionary->[$i];
		$variation->{patients_object}->{$v} =1;
	}
	delete  $variation->{pato};
	

	if (exists $variation->{values}){
		die() unless is_deeply($variation->{values}->{genes_object}, $variation->{genes_object}, 'data structures should be the same');
	}
	if (exists $variation->{values}){
		die() unless is_deeply($variation->{values}->{references_object}, $variation->{references_object}, 'data structures should be the same');
	}
	if (exists $variation->{values}){
		die() unless is_deeply($variation->{values}->{patients_object}, $variation->{patients_object}, 'data structures should be the same');
	}
	
	
	
	### BIT
	my @bits = split("",unpack($variation->{st_bit},$variation->{bit}) );
	my $i=0;
	foreach my $id (@{$self->array_for_bit}){
		$variation->{$id} = $bits[$i];
		$i++;
	}	
	delete $variation->{st_bit};
	delete  $variation->{bit};
	#FLOAT 
	my @float = unpack($variation->{st_float},$variation->{float});
	$i=0;
	
	foreach my $id (@{$self->array_for_float}){
		$variation->{$id} = $float[$i];
		$i++;
	}
	delete $variation->{st_float};
	delete  $variation->{float};
	#INTEGER
	my @integer = unpack($variation->{st_int},$variation->{int});
	
	$i=0;
	
	foreach my $id (@{$self->array_for_integer}){
		$integer[$i] = undef if $integer[$i] eq "U";
		$variation->{$id} = $integer[$i];
		$i++;
	}
	delete $variation->{st_int};
	delete  $variation->{int};
	# STRING 
	$i=0;
	my @ast = split("!",$variation->{fst});
	 die(scalar(@ast)." ".scalar @{$self->array_for_string}) if scalar(@ast) ne scalar @{$self->array_for_string};
	foreach my $id (@{$self->array_for_string}){
		$variation->{$id} = $ast[$i];
		$i++;
	}
	
	delete $variation->{fst};
	if (exists $variation->{values}){
	foreach my $k (@{$self->array_for_bit}){
		die("bit".Dumper ($variation->{values})." ".Dumper($variation)) if $variation->{values}->{$k} ne $variation->{$k};
	}
	foreach my $k (@{$self->array_for_integer}){
		next unless $variation->{values}->{$k};
		die("integer $k :".Dumper ($variation->{values})." ".Dumper($variation)) if $variation->{values}->{$k} ne $variation->{$k};
	}
		foreach my $k (@{$self->array_for_float}){
				next unless $variation->{values}->{$k};
				next if $variation->{values}->{$k} eq "-";
				
		die ("float".Dumper ($variation->{values})." ".Dumper($variation)) if $variation->{values}->{$k} ne $variation->{$k};
	}
	}
	
	###################################
	# ANNEX 		patients_object
	#####################################
	delete $variation->{annex};




###################################		
#		ANNTOTATION
##################################		
foreach my $id (keys %{$variation->{annotation}}) {
		my $o = $self->array_dictionary->[$id];
		$variation->{annotation}->{$o} = delete $variation->{annotation}->{$id};
}


#####################################
#		hash polyviewer Annotations
#####################################	
if (exists $variation->{age})	{
	$variation->{annotation_polyviewer} = $self->decode_annotation_polyviewer(delete($variation->{age}));
 $variation->{age} = {};
delete $variation->{age};
}	

########################################
# DEJAVU
######################################
	my $string;
	my $nb = scalar(@{$self->array_for_dejavu});
	my @t = unpack("w$nb",$variation->{st_dejavu});
	for (my $i = 0;$i< @t;$i++){
		$variation->{value}->{$self->array_for_dejavu->[$i]} = $t[$i];
	} 
	delete $variation->{st_dejavu};



 if (exists $variation->{values}){
 	die(Dumper($variation->{annotation})." ".Dumper($variation->{values}->{annotation})) unless is_deeply($variation->{values}->{annotation}, $variation->{annotation}, 'data structures should be the same');
 }

		
 if (exists $variation->{values}){
 
	foreach my $k (@{$self->array_for_bit}){
		die("bit".Dumper ($variation->{values})." ".Dumper($variation)) if $variation->{values}->{$k} ne $variation->{$k};
	}
	foreach my $k (@{$self->array_for_integer}){
		next unless $variation->{values}->{$k};
		#warn $k." ". $variation->{values}->{$k}." :: ".$variation->{$k};
		die("integer $k :".Dumper ($variation->{values})." ".Dumper($variation)) if $variation->{values}->{$k} ne $variation->{$k};
	}
		foreach my $k (@{$self->array_for_float}){
			
				next if !($variation->{values}->{$k}) && !($variation->{$k});
				next if $variation->{values}->{$k} eq "-";
			#warn $k." ". $variation->{values}->{$k}." :: ".$variation->{$k};
				
		die ("float".Dumper ($variation->{values})." ".Dumper($variation)) if $variation->{values}->{$k} ne $variation->{$k};
	}
#	use Test::More tests => 1000;
	warn "TEST 4 ";
	die() unless is_deeply($variation->{values}->{seq2}, $variation->{seq}, 'data structures should be the same');
	die() unless is_deeply($variation->{values}->{annex}, $variation->{annex}, 'data structures should be the same');
	}
		
	#use Test::More tests => 1;
	#die() unless is_deeply($variation->{annex2}, $variation->{annex}, 'data structures should be the same');		
	
	return $variation;
	
}


sub getPolyViewerVariant {
	my ($self,$v,$vp,$project,$gene,$patient) = @_;
	
	############################
	# General
	############################
	#$vp->gene($gene);
	$vp->id($v->id);
	$vp->start($v->start);
	$vp->end($v->end);
	$vp->ref_allele($v->ref_allele);	
	$vp->allele($v->alternate_allele);	
	$vp->gnomad_id($v->gnomad_id);
	$vp->chromosome($v->getChromosome->name);
	$vp->name($v->name);
	$vp->type($v->type);
	
	#######################
	# HGMD CLINVAR
	#######################
		
		$vp->hgmd($v->hgmd_class);
		$vp->hgmd_id($v->hgmd_id);
		if (exists $v->genes_pathogenic_DM->{$gene->{id}} && $v->genes_pathogenic_DM->{$gene->{id}}->{DM} ){
				$vp->dm($v->isDM);
				$vp->hgmd_phenotype($v->hgmd_phenotype);
			}
		
		
		$vp->clinvar_id($v->clinvar_id);
		$vp->clinvar($v->clinvar_class);
		if (exists $v->genes_pathogenic_DM->{$gene->{id}} && $v->genes_pathogenic_DM->{$gene->{id}}->{pathogenic} ){
				$vp->clinvar_pathogenic($v->isClinvarPathogenic);
		}
	
	############################
	# Calling
	############################
	
		$vp->isCnv($v->isCnv);
		
		my $hpatients;
		#if ($v->isDude){
		#	$vp->setDudeValues($v->getChromosome,$patient,$v);
		#}
		#elsif ($v->isCnv){
		#	$vp->setCnvValues($v->getChromosome,$patient,$v);
		#}
		#else {
			foreach my $p (@{$patient->getFamily()->getMembers}) {
				$hpatients->{ $p->id }->{gt} = $v->getSequencingGenotype($p);
				$hpatients->{ $p->id }->{pc} = $v->getRatio($p);
				$hpatients->{ $p->id }->{dp} = $v->getDP($p);
				$hpatients->{ $p->id }->{model} = $v->getTransmissionModelType($p->getFamily(),$p);
			}
		#}
	
	$vp->patients_calling($hpatients);
	##################
	# gnomad
	##################
		
	$vp->gnomad_ac($v->getGnomadAC);
	$vp->gnomad_an($v->getGnomadAN);
	$vp->gnomad_min_pop_name($v->min_pop_name);
	$vp->gnomad_min_pop($v->min_pop_freq);
	$vp->gnomad_max_pop_name($v->max_pop_name);
	$vp->gnomad_max_pop($v->max_pop_freq);
	$vp->gnomad_ho($v->getGnomadHO);
	$vp->gnomad_ho_male($v->getGnomadAC_Male);
	#########
	# DEJAVU
	#########
	$vp->dejavu_other_projects($v->other_projects());
	$vp->dejavu_other_patients($v->other_patients());
	$vp->dejavu_other_patients_ho($v->other_patients_ho());
	$vp->dejavu_similar_projects( $v->similar_projects());
	$vp->dejavu_similar_patients($v->similar_patients());
	$vp->dejavu_similar_patients_ho($v->similar_patients_ho());
	$vp->dejavu_this_run_patients($v->in_this_run_patients);# = '-';
	#############
	# SCORE
	############
	#spliceAI
	$vp->spliceAI($v->max_spliceAI_score($gene));
	$vp->spliceAI_cat($v->max_spliceAI_categorie($gene));
	$vp->cadd($v->cadd_score);
	$vp->revel($v->revel_score);
	$vp->rf( $v->dbscsnv_rf);
	$vp->ada( $v->dbscsnv_ada);
	################
	# transcripts
	##################
	my $tgene = delete $v->{annotation_polyviewer}->{$gene->{id}};
	
	
	#$hvariation->{genes} = delete $v->{annotation_polyviewer};
	
	my $gene_id = $gene->id;
	
	#$vp->transcripts([]);
	$vp->transcripts($tgene->{trs});	
	return $vp;
	
} 

sub transforminHash {
	my ($self,$v,$patient) = @_;
	
	my $project = $patient->project;
	my $hvariation;
	$hvariation->{value} = delete $v->{value};
	$hvariation->{value}->{id} =  $v->id;
	$hvariation->{value}->{name} =  $v->name;
	$hvariation->{value}->{type} = $v->type;
	$hvariation->{value}->{patients} = "";
	
#	foreach my $p (@{$v->getPatients}){
#		$hvariation->{value}->{patients} .= $p->name.";";
#		
#		
#	}
	$hvariation->{value}->{gnomad_id}  = $v->gnomad_id;
	
	$hvariation->{value}->{is_cnv} = 0;
	$hvariation->{value}->{is_cnv} = 1 	if ($v->isCnv());
	$hvariation->{value}->{allele} =$v->getSequence;
	$hvariation->{value}->{ref_allele} =$v->ref_allele;
	$hvariation->{value}->{start} =$v->start;
	$hvariation->{value}->{end} =$v->end;
	$hvariation->{value}->{chromosome} = $v->getChromosome()->name;
		
	####################
	#gnomad 
	#####################
	$hvariation->{value}->{ac} = $v->getGnomadAC;
	$hvariation->{value}->{max_pop} = "-";
	$hvariation->{value}->{max_pop} =$v->max_pop_name.":".sprintf("%.4f", $v->max_pop_freq ) if $v->max_pop_name;
	
	$hvariation->{value}->{min_pop} = "-";
	$hvariation->{value}->{min_pop} =$v->min_pop_name.":".sprintf("%.4f", $v->min_pop_freq ) if $v->min_pop_name;
	$hvariation->{value}->{an} = $v->getGnomadAN;
	$hvariation->{value}->{ac_ho} = $v->getGnomadHO;
	$hvariation->{value}->{ac_ho_Male} = $v->getGnomadAC_Male;
	
	####################
	#DEJAVU
	#####################

	if ($v->isCnv()) {
		$hvariation->{value}->{is_cnv} = 1;
		#$hvariation->{value}->{sd_value_controls} = sprintf("%.2f",$patient->sd_value_dude($chr_id,$start,$end));
		$hvariation->{value}->{cnv_details_genes} = $v->get_genes_transcripts_details_dup_del();
		if ($v->getProject->isGenome()) {
		#	$hvariation->{value}->{cnv_confidence}->{$patient->name()} = $v->cnv_confidence($patient);
		}
	
	}
#	$v->getTranscripts();
##############
## SCORE
#################
	 $hvariation->{value}->{cadd} = $v->cadd_score;
	 $hvariation->{value}->{revel} = $v->revel_score;
	 $hvariation->{value}->{rf} = $v->dbscsnv_rf();
	 $hvariation->{value}->{ada} = $v->dbscsnv_ada();
	 $hvariation->{value}->{large_evt} = 1 if $v->isLargeDeletion or $v->isLargeDuplication;
	 $hvariation->{genes} = delete $v->{annotation_polyviewer};
	return $hvariation;
}

sub getHash{
	my ($self,$key,$patient) = @_;
	
	my $db_name = $self->name();
	my $variation =   $self->decode($self->lmdb($db_name)->get($self->lmdb_key($key)));
	warn $key unless $variation;
	$variation->{buffer} = $patient->buffer;
	$variation->{project} = $patient->project;
	return $self->transforminHash($variation,$patient);
}

sub gethashVariant{
	
}

sub close {
	my ($self) = @_;
	$self->save_dictionary() if $self->mode() ne "r";
	 $self->SUPER::close(); 
	
}

1;