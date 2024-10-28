package save;
use FindBin qw($RealBin);
use strict;
use Moo;
use Data::Dumper;
use Storable qw(store retrieve freeze thaw fd_retrieve);

use Data::Dumper;
use lib "$RealBin";
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
use GenBoNoSqlLmdb; 
use GenBoNoSqlLmdbInteger;
use GenBoNoSql;
 use Statistics::Descriptive;
 use Set::IntSpan::Fast::XS;
 # use RocksDB;
use List::Util qw(sum max min);
use Storable qw/thaw freeze/;
 #use Compress::Snappy;
has name =>(
	is		=> 'ro',
	required=>1,
);

has version =>(
	is		=> 'ro',
);
has db_type =>(
	is		=> 'ro',
	required=>1,
);

has chromosome =>(
	is		=> 'ro',
	required=>1,
);

has mode =>(
	is		=> 'ro',
	default => sub {
		return 'w';
	}
);
has integer =>(
	is		=> 'ro',
	default => sub {
		return '0';
	}
);
has compress =>(
	is		=> 'ro',
	default => sub {
		return '1';
	}
);
has temp =>(
	is		=> 'ro',
	default => sub {
		return '0';
	}
);

has dir_tmp =>(
	is		=> 'ro',
	lazy => 1,
	default => sub {
		my ($self) = @_;
		#return "/tmp/";
		return  "/tmp/lmdb/";
			
	},
);

has dir_lmdb =>(
	is		=> 'rw',
	lazy => 1,
	default => sub {
		my ($self) = @_;
		#return "/tmp/";
		return  "/tmp/lmdb/";
			
	},
);


has final_dir_lmdb =>(
	is		=> 'ro',
	lazy => 1,
	default => sub {
		my ($self) = @_;
		my $dir = $self->dir_lmdb."".$self->name."/";
		$dir = $self->dir_tmp."".$self->name."/" if $self->temp == 1;
		system("mkdir -p $dir") unless -e $dir;
		return $dir ;
			
	},
);

sub move_lmdb {
	my ($self,$type) = @_;
	my $file = $self->dir_tmp_lmdb."/".$type."/".$self->chromosome()."*";
	my $dir = $self->final_dir_lmdb."/".$type;
	system("mkdir -p $dir") unless -e $dir;
	#warn "mv $file $dir/";
	system("mv $file $dir/");
	
}


has dir_rocks =>(
	is		=> 'ro',
	lazy => 1,
	default => sub {
		my ($self) = @_;
		my $dir = "/tmp/rocks/".$self->name."/";
		system("mkdir -p $dir") unless -e $dir;
		return $dir ;
			
	},
);
has dir_intspan =>(
	is		=> 'ro',
	lazy => 1,
	default => sub {
		my ($self) = @_;
		my $dir = $self->dir_lmdb."/intspan/".$self->name;
		system("mkdir -p $dir") unless -e $dir;
		return $dir ;
			
	},
);

sub file_intspan {
	my ($self,$type) = @_;
	my $chr = $self->chromosome();
	my $dir =$self->dir_intspan."/$type/";
	mkdir $dir unless -e $dir;
	return $self->dir_intspan."/$type/"."$chr.intspan.freeze";
}

sub get_intspan {
	my ($self,$type) = @_;
	return $self->{intspan}->{$type} if exists $self->{intspan}->{$type};
	$self->{intspan}->{$type} =  Set::IntSpan::Fast::XS->new();
	return  Set::IntSpan::Fast::XS->new();
}


sub get_saved_intspan {
	my ($self,$type) = @_;
	my $chr = $self->chromosome();
	my $file = $self->file_intspan($type);
	return Set::IntSpan::Fast::XS->new() unless -e $file;
	return retrieve($file);
	
}



sub save_intspan {
		my ($self) = @_;
		my $chr = $self->chromosome();
		 	my @types =  keys %{$self->{intspan}};#("snps","insertions","deletions");
		 	foreach my $type (@types){
		 		my $file =  $self->file_intspan($type);
#		 		warn $file;
		 		 store($self->{intspan}->{$type}, $file) if $self->{intspan}->{$type};
		 	}
}



sub close {
	my ($self) = @_;
	foreach my $t (keys %{$self->{db}}){
		 $self->{db}->{$t}->close();
		 delete $self->{$t};
	}
	
}
sub nb_keys {
	my ($self) = @_;
	return  $self->get_snp_db("snps")->nb_keys;
}

sub get_lmdb_db {
	my ($self,$type) = @_;
	my $mode = $self->mode;
	#die("no mode ") unless $mode;
	return $self->{db}->{$type} if exists $self->{db}->{$type};
	my $dir;
	if ($self->mode eq 'r'){
			$dir  = $self->final_dir_lmdb."/$type";
	}
	else {
	 	$dir  = $self->final_dir_lmdb."/$type";
	
	}
	mkdir $dir unless -e $dir;
	
	if ($self->mode eq 'w' or $self->mode eq 'c'){
		if ($type eq 'relation_variant_gene') { $self->{db}->{$type}= GenBoNoSqlLmdb->new(dir=>$dir,mode=>$self->mode,name=>$self->chromosome,is_compress=>$self->compress); }
		elsif ($self->integer == 0) { $self->{db}->{$type}= GenBoNoSqlLmdb->new(dir=>$dir,mode=>$self->mode,name=>$self->chromosome,is_compress=>$self->compress); }
		elsif ($self->integer == 1) { $self->{db}->{$type}= GenBoNoSqlLmdb->new(dir=>$dir,mode=>$self->mode,name=>$self->chromosome,is_integer=>1,is_compress=>$self->compress); }
	}
	else {
		if ($type eq 'relation_variant_gene') { $self->{db}->{$type}= GenBoNoSqlLmdb->new(dir=>$dir,mode=>$self->mode,name=>$self->chromosome,is_compress=>1,is_index=>1); }
		else { $self->{db}->{$type}= GenBoNoSqlLmdb->new(dir=>$dir,mode=>$self->mode,name=>$self->chromosome,is_compress=>1,integer=>1,is_index=>1) ; }
	}
	return $self->{db}->{$type};
}


sub get_rocks_db {
	my ($self,$type,$mode) = @_;
	return $self->{rocks}->{$type} if exists   $self->{rocks}->{$type};
	mkdir $self->dir_rocks."/$type/" unless -e $self->dir_rocks."/$type/";
	warn  $self->dir_rocks."/$type/";
	my $db;
	if ($self->mode eq "r"){
	 $db = RocksDB->new($self->dir_rocks."/$type/".$self->chromosome, { read_only => 1,max_open_files=> 256 });
	}
	else {
	 $db = RocksDB->new($self->dir_rocks."/$type/".$self->chromosome, { create_if_missing => 1,max_open_files=> 256 });
	}
	  $self->{rocks}->{$type} = $db;
	  return $db;
}

sub get_rocks {
	my ($self,$type,$pos) = @_;
	my $z = $self->get_rocks_db($type)->get($pos);
	return thaw $z if $z;
	return undef;
	
}

sub save_caad {
	my ($self,$variation,$debug) = @_;
	my $id_pos = $variation->{start};
	my $alternate_allele = $variation->{alternate_allele};
	my $db = $self->get_snp_db_lite("caad");

	my $h = $db->get($id_pos);
		unless ($h){
			$h ={};
		}
	
	return unless $alternate_allele;
	$h->{$alternate_allele}->{c} = [ $variation->{caad_r},$variation->{caad_p}];
	$db->put($id_pos,$h);
		my $hh;
	$hh->{caad} = [ $variation->{caad_r},$variation->{caad_p}];
	my $dbr = $self->get_rocks_db("caad");
	$dbr->put($id_pos."_$alternate_allele",freeze($hh));
#	my $db2 = $self->get_snp_db_lite("caad2");
#	
#	my $hh;
#	$hh->{caad} = [ $variation->{caad_r},$variation->{caad_p}];
#	#$hh->{c} =  $variation->{caad_r}.";".$variation->{caad_p};
#	#$hh->{cr} =  $variation->{caad_r};
#	#$hh->{cp} = $variation->{caad_p};
#	
}


#sub save {
#		my ($self,$variation,$debug) = @_;
#		my $db = $self->get_snp_db($variation->{type});
#		my $id_pos = $variation->{start};
#		my $h = $db->get($id_pos);
#		warn Dumper $h if $debug;
#		warn "---------2--------------" if $debug;
#		return unless $id_pos;
#		#	my $debug;
#		#	$debug =1 if $id_pos == 36042179  or $id_pos == 33355899;
#		#	warn Dumper $variation if $debug;
#		unless ($h){
#			$h ={};
#		#	warn "new";
#		}
#		else {
#			#warn "ok";
#		}
#		my $alternate_allele = $variation->{alternate_allele};
#		return unless $alternate_allele;
#			
#		$h->{$alternate_allele}->{clinical} = $variation->{clinical} if exists $variation->{clinical} ;
#		 
#
#		
#		$h->{$alternate_allele}->{db_frequence}->{$self->name} = $variation->{frequence}  if exists $variation->{frequence};
#		$h->{$alternate_allele}->{db_frequence_hom}->{$self->name} = $variation->{frequence_hom} if exists $variation->{frequence_hom};
#		$h->{$alternate_allele}->{db_rsname}->{$self->name} =  $variation->{rs} if $variation->{rs} ;
#		
#		$h->{$alternate_allele}->{db_clinical}->{$self->name} =  $variation->{clinical} if exists $variation->{clinical} ;		
#
#		unless (exists $h->{$alternate_allele}->{name}){
#			$h->{$alternate_allele}->{name} = $variation->{rs} if $variation->{rs} ;
#		}
#		$h->{$alternate_allele}->{clinical} = $variation->{clinical} if exists $variation->{clinical} ;	
#		$h->{$alternate_allele}->{db}->{$self->name} =1;
#		$db->put($id_pos,$h);
#		$self->update_caad($id_pos,$variation->{type},$h);
#		$self->save_lite($id_pos,$variation->{type},$h);
#		#my $db2 = $self->get_snp_db_lite($variation->{type});
#		if ($debug){
#			warn Dumper $db->get($id_pos);
#		}
#}

sub open_no_sql{
	my ($self,$type) = @_;
	return $self->{sqlite}->{$type} if exists  $self->{sqlite}->{$type};
		 $self->{sqlite}->{$type} = GenBoNoSql->new(dir =>"/tmp/sqlite/".$type, mode => 'c');

		  return  $self->{sqlite}->{$type};
}

sub save_relation_variant_gene {
	my ($self,$relation) = @_;
	$self->save_relation_variant_gene_lmdb($relation);
	return;
}

sub save_object {
	my ($self,$variation) = @_;
	if ($self->db_type eq "sqlite"){
		$self->save_object_lite($variation);
		return;
	}
	if ($self->db_type eq "lmdb"){
		$self->save_object_lmdb($variation);
		return;
	}
	confess("db_type must be sqlite or lmdb and not ".$self->db_type);
}

sub save_object_lite {
	my ($self,$variation) = @_;
	my $db = $self->open_no_sql($variation->{type});
	my $id = $variation->{start};#."_".$variation->{alternate_allele};
	my $fh =  $db->get($self->chromosome,$id);
	my $alt = $variation->{alternate_allele};
	my $h;
	#$h->{start} = $variation->{start};
	#$h->{alternate_allele} = $variation->{alternate_allele};
 	if (exists $variation->{clinical}){
		$h->{clinical} = $variation->{clinical} if exists $variation->{clinical};
		$h->{sig} = $variation->{sig} if exists $variation->{sig};
		$h->{clnid} = $variation->{id} if exists $variation->{id};
 	}
 	
 	$h->{AC} = $variation->{AC}  if exists $variation->{AC};
 	$h->{AN} = $variation->{AN}  if exists $variation->{AN};
 		$h->{AF} = $variation->{AF}  if exists $variation->{AF};
 	$h->{Hom} = $variation->{Hom}  if exists $variation->{Hom};
 	$h->{FILTER} = $variation->{FILTER}  if exists $variation->{FILTER};
 	#$h->{INFO}  =  $variation->{INFO}  if exists $variation->{INFO};
 	#AC_AN_Hom_FILTER
	$h->{frequence} = $variation->{frequence}  if exists $variation->{frequence};
	$h->{frequence_hom} = $variation->{frequence_hom}  if exists $variation->{frequence_hom};
	$h->{rsname} =  $variation->{rs} if $variation->{rs} ;
	#$h->{db} = $self->name;
	$h->{short_db} = $self->get_short_db($self->name);
	#my $lmdb = $self->get_lmdb_db($type);
	$fh->{$alt} = $h;
	$db->put($self->chromosome,$id,$fh);
	
	
}

sub save_relation_variant_gene_lmdb {
	my ($self,$relation) = @_;
	my $type = 'relation_variant_gene';
	my $lmdb = $self->get_lmdb_db($type);
	my $type_id = $relation->{type_id};
	delete $relation->{type_id};
	$relation->{gene} = 'none' unless ($relation->{gene});
	my $id = $relation->{$type_id}.'_'.$relation->{gene};
	$lmdb->put($id, $relation);
	my $id2 = $relation->{$type_id};
	my $fh = $lmdb->get($id2);
	unless ($fh) { $fh= {}; }
	$fh->{$relation->{gene}} = undef;
	$lmdb->put($id2, $fh);
}

sub save_object_lmdb {
	my ($self,$variation) = @_;
	my $type = $variation->{type};
	my $lmdb = $self->get_lmdb_db($type);
	my $id = $variation->{start};
	my $fh =  $lmdb->get($id);
	$self->get_intspan($type)->add($id);
	unless ($fh){
		$fh= {};
	}
	#."_".$variation->{alternate_allele};
	my $alt = $variation->{alternate_allele};
	my $h;
	#$h->{start} = $variation->{start};
	#$h->{alternate_allele} = $variation->{alternate_allele};
	my $field =[keys %$variation];
	
	foreach my $f (@$field){
		$h->{$f} = $variation->{$f} if exists $variation->{$f};
	}
	if (exists $variation->{fields}){
		foreach my $kf (keys %{$variation->{fields}}){
			$h->{$kf} = $variation->{fields}->{$kf};
		}
	
	}
	
	$h->{short_db} = $self->get_short_db($self->name);
	#my $lmdb = $self->get_lmdb_db($type);
	
	$fh->{$alt} = $h;
	$lmdb->put($id,$fh);
	#$self->get_rocks_db($type)->put($id,freeze ($fh));
}
sub hash_snp {
	my ($self,$variation) = @_;
	my $vcfRefAllele = $variation->{ref_all};
	 if ($vcfRefAllele eq 'N'){
			confess();
	 }
if (length( $variation->{ht}) > 1){
return;
}
	if (length($vcfRefAllele) == 1){
			$variation->{sequence} =  $variation->{ht};
			$variation->{alternate_allele} = $variation->{ht};
	}
	else {
		my @alls = split("",$variation->{alt_all});
		my @refs = split("",$variation->{ref_all});
		my $first_pos = 0;
		for (my $i =0;$i<@alls;$i++){
			if ($alls[$i] eq $variation->{ht} && $refs[$i] ne $alls[$i]  ){
				$first_pos = $i;
				die();
			}
		} 
		
		$variation->{alternate_allele} = $variation->{ht};
		$variation->{start} +=  $first_pos;
	}
	$variation->{type} = "snps";
	return $variation;
}

sub add_snp {
	my ($self,$variation) = @_;
	my $vcfRefAllele = $variation->{ref_all};
	 if ($vcfRefAllele eq 'N'){
			confess();
	 }
if (length( $variation->{ht}) > 1){
return;
}
	if (length($vcfRefAllele) == 1){
			$variation->{sequence} =  $variation->{ht};
			$variation->{alternate_allele} = $variation->{ht};
	}
	else {
		my @alls = split("",$variation->{alt_all});
		my @refs = split("",$variation->{ref_all});
		my $first_pos = 0;
		for (my $i =0;$i<@alls;$i++){
			if ($alls[$i] eq $variation->{ht} && $refs[$i] ne $alls[$i]  ){
				$first_pos = $i;
			}
		} 
		
		$variation->{alternate_allele} = $variation->{ht};
		$variation->{start} +=  $first_pos;
	}
	$variation->{type} = "snps";
	$self->save_object($variation);
}

sub add_relation_variant_gene {
	my ($self, $relation) = @_;
	$self->save_relation_variant_gene($relation);
}

sub hash_deletion {
	my ( $self,$variation,$debug) = @_;
	$variation->{type} = 'deletions';
	my $vcfRefAllele = $variation->{ref_all};
	my $first = substr($vcfRefAllele, 0, 1); 
	$variation->{sequence} =   "-";
	$variation->{alternate_allele} = $variation->{ht};
	my $dpos = index($vcfRefAllele,$variation->{ht});
	$variation->{sequence} =   "-";
	confess($vcfRefAllele."-->".$variation->{ht}) if $dpos == -1;
	$variation->{start} +=   $dpos;
	return $variation;
	#$self->save_object($variation);
}

sub add_deletion {
	my ( $self,$variation,$debug) = @_;
	$variation->{type} = 'deletions';
	my $vcfRefAllele = $variation->{ref_all};
	my $first = substr($vcfRefAllele, 0, 1); 
	$variation->{sequence} =   "-";
	$variation->{alternate_allele} = $variation->{ht};
	my $dpos = index($vcfRefAllele,$variation->{ht});
	$variation->{sequence} =   "-";
	confess($vcfRefAllele."-->".$variation->{ht}) if $dpos == -1;
	$variation->{start} +=   $dpos;
	$self->save_object($variation);
}
sub hash_insertion {
my 	( $self,$variation) = @_;
	$variation->{type} = 'insertions';
	my $vcfRefAllele = $variation->{ref_all};
	my $vcfVarAllele = $variation->{alt_all};
	my $first = substr($vcfRefAllele, 0, 1); 
    $variation->{alternate_allele} = $variation->{ht};
	my $dpos = index($vcfVarAllele,$variation->{ht});
	if ($dpos == -1){
							#cas de l'insertion au milieu de la sequence 
							# exemple : TATTTTA ref :TTTTTA
								confess($vcfVarAllele." ref :".$vcfRefAllele." ".$variation->{ht}." ") if $dpos == -1;
		}
		else {
					
								$variation->{start} +=  $dpos ;
			}
		return $variation;
}


sub add_insertion {
my 	( $self,$variation) = @_;
	$variation->{type} = 'insertions';
	my $vcfRefAllele = $variation->{ref_all};
	my $vcfVarAllele = $variation->{alt_all};
	my $first = substr($vcfRefAllele, 0, 1); 
    $variation->{alternate_allele} = $variation->{ht};
	my $dpos = index($vcfVarAllele,$variation->{ht});
	if ($dpos == -1){
							#cas de l'insertion au milieu de la sequence 
							# exemple : TATTTTA ref :TTTTTA
								confess($vcfVarAllele." ref :".$vcfRefAllele." ".$variation->{ht}." ") if $dpos == -1;
		}
		else {
					
								$variation->{start} +=  $dpos ;
			}
			$self->save_object($variation);
}

sub get_corrected_freq{
	my ($self,$freq) = @_;
	my $frequence;
	my @correct_freq;
	 sub outlier_filter {
 	my ($stat,$value) = @_;
 	#warn $value;
 	 #my $z = abs($stat->mean - $value)/$stat->standard_deviation();
 	 #warn $z;
 	 #return $z;
 	 };
 	 my $nbf = scalar(@$freq);
 	 if ($nbf==0){
	 	 $frequence = ".";
	 	 @correct_freq = ();
	 }
	  if ($nbf==0){
	 	 $frequence = ".";
	 	 @correct_freq = ();
	 }
	 			elsif ($nbf == 1){
	 	 					$frequence = $freq->[0];
	 	 					@correct_freq = @$freq;
	 	 				}
	 	 				else{
	 	 					my $stat = Statistics::Descriptive::Full->new();
	 	 					$stat->add_data(@$freq);
	 	 					
	 	 					my $sd = $stat->standard_deviation();
	 	 					$sd = 0 unless $sd;
	 	 					die() unless defined $sd;
	 	 					if ($sd < 0.05){
	 	 						$frequence = $stat->mean();
	 	 						@correct_freq = @$freq;
	 	 					}
	 	 					else {
	 	 						my $mean;
	 	 						if ($nbf  == 2){
	 	 							$frequence = sum(@$freq) - max(@$freq);
	 	 							@correct_freq = min(@$freq);
	 	 						}
	 	 						else {
	 	 							$Statistics::Descriptive::Min_samples_number =3;
	 	 							
	 	 							$stat->set_outlier_filter( \&outlier_filter );
	 	 			                my @data = $stat->get_data_without_outliers();
	 	 			                $frequence = sum(@data) /scalar(@data);
	 	 			                @correct_freq = @data;
	 	 							
	 	 						}
	 	 					}
	 	 				}
return ($frequence,\@correct_freq);	 	 				
}
sub get_short_db {
	my ($self,$db) = @_;
	return "ge" if $db eq "gnomad-exome";
	return "gg" if $db eq "gnomad-genome";
	return "ex" if $db eq "exac";
	return "db" if $db eq "dbsnp";
	return "ev" if $db eq "evs";
	return "co" if $db eq "cosmic";
	return "cl" if $db eq "clinvar";
	#return "lc" if $db eq "cldb";
	return $db;
}



1;