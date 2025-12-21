package GenBoNoSqlRocksTinyPolyviewerVariant;
use FindBin qw($RealBin);
use Moo; 
use strict;
use Data::Dumper;
use JSON::XS;
use Digest::MD5 qw(md5_hex);
use Carp;
use File::Basename;
use GenBoNoSqlRocks;
use MCE::Loop;
use MCE::Shared;
use Storable qw/thaw freeze retrieve store/;
use Sereal qw(sereal_encode_with_object sereal_decode_with_object write_sereal_file read_sereal_file);
use List::Util qw(sum);
use POSIX qw(log);
has columns_file =>(
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my ($self) = @_;
		my $file_freeze = $self->dir()."/".$self->project->name.".store";
		#die($file_freeze) unless -e $file_freeze;
		return $file_freeze;
	}
); 
##  


my $htr = {
 'dbscsnv' => "INTEGER",
 'impact_score' => "INTEGER",
 'nomenclature' => 'VARCHAR',
 'end' => "INTEGER",
 'alphamissense' =>"INTEGER",
 'cadd' => "INTEGER",
 'codons' => 'VARCHAR',
 'name' => 'ENST00000450305',
 'consequence' => 'Upstream',
 'exon' => 'VARCHAR',
 'mane' => 'BOOLEAN',
 'spliceAI_cat' => 'VARCHAR',
 'impact_score_text' => 'VARCHAR',
 'main' => "BOOLEAN",
 'start' => "INTEGER",
 'revel' => "INTEGER",
 'polyphen' => "REAL",
 'sift' => "REAL",
 'prot' => 'VARCHAR',
 'appris' =>'VARCHAR',
 'spliceAI' => "INTEGER",
 'promoterAI' => "VARCHAR",
 'codons_AA' => 'VARCHAR',
 'nm' => 'VARCHAR'
};

has columns_transcripts  =>(
	is		=> 'rw',
	lazy    => 1,
	default => sub {
my @column_transcripts = ();
push(@column_transcripts,sort{$a cmp $b} keys %$htr);
return \@column_transcripts;
}
);
has columns_patients  =>(
	is		=> 'rw',
	lazy    => 1,
	default => sub {
	my $hpat = {
														'pc' => 'INTEGER',
                                                       'id' => 'INTEGER',
                                                       'gt' => 'VARCHAR',
                                                       'model' => 'INTEGER',
                                                       'pr' => 'INTEGER',
                                                       'sr' => 'INTEGER',
                                                       'name' => 'VARCHAR',
                                                       'dp' => 'INTEGER',
                                                       'array_text_calling' => "VARCHAR",#[

	};
my @column_patient = ();
push(@column_patient, sort{$a cmp $b} keys %$hpat);
return \@column_patient;
}
);
has columns_patients_sv  =>(
	is		=> 'rw',
	lazy    => 1,
	default => sub {
	my $hpat = {
														'pc' => 'INTEGER',
                                                       'id' => 'INTEGER',
                                                       'gt' => 'VARCHAR',
                                                       'model' => 'INTEGER',
                                                       'pr' => 'INTEGER',
                                                       'sr' => 'INTEGER',
                                                       'name' => 'VARCHAR',
                                                       'dp' => 'INTEGER',
                                                       'array_text_calling' => "VARCHAR",#[
													   'log2_ratio' => "VARCHAR",#[
														'norm_depth' => "INTEGER",
														'norm_depth_after' => "INTEGER",
														'norm_depth_before' => "INTEGER",
	};
	
my @columns_patients_sv = ();
push(@columns_patients_sv, sort{$a cmp $b} keys %$hpat);
return \@columns_patients_sv;
}
);

has columns_global  =>(
	is		=> 'rw',
	lazy    => 1,
	default => sub {
	my $hglobal = {          
                 'cadd' => '-',
                 'allele' => 'C',
                 'clinvar' => undef,
                 'gnomad_ho_male' => 4712,
                 'hgmd_id' => undef,
                 'dejavu_other_patients_ho' => 322,
                 'rocksdb_id' => '0000010146!1',
                 'gnomad_max_pop' => '0.654008448123932',
                 'gnomad_max_pop_name' => 'asj',
                 'dejavu_other_patients' => 1179,
                 'gnomad_min_pop' => '0.333333343267441',
                 'isSrPr' => undef,
                 'gnomad_min_pop_name' => 'mid',
                 'end' => 10147,
                 'gnomad_an' => 16096,
                 'ada' => '-',
                 'isJunction' => 0,
                 'gnomad_id' => '1-10146-AC-A',
                 'dejavu_this_run_patients' => 7,
                 'dejavu_similar_patients' => 0,
                 'hgmd' => undef,
                 'isCnv' => 0,
                 'type' => 'deletion',
                 'chromosome' => '1',
                 'name' => '1-10146-AC-A',
                 'rf' => '-',
                 'ref_allele' => 'C',
                 'dejavu_other_projects' => 749,
                 'dejavu_similar_projects' => 0,
                 'reference' => 'GenBoDeletion',
                 'id' => '1_10147_AC_A',
                 'gnomad_ho' => 2613,
                 'start' => 10147,
                 'dejavu_similar_patients_ho' => 0,
                 'global_vector_id' => '1!0',
                 'gnomad_ac' => 10095
	};
	my @column_global;
	push(@column_global,sort{$a cmp $b} keys %$hglobal); 
	return \@column_global;
}
);

has print_html =>(
	is		=> 'rw',
);

has mode =>(
	is		=> 'ro',
	required => 1,
);

has project =>(
	is		=> 'ro',
	required => 1,
);
has patient =>(
	is		=> 'rw',
	#required => 1,
);
has dir =>(
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my ($self) = @_;
		
		my $diro2 = $self->project->tiny_rocks_cache_dir;
		
		system("mkdir $diro2") unless -d $diro2;
		return $diro2;
	}
); 



sub add_index {
	my ($self,$chr,$id) = @_;
	$self->{list_index}->{$chr}->{$id} ++
}

sub indexes{
	my ($self) = @_;
	my @all_variants;
	for my $chr (keys %{$self->{list_index}}) {
    push @all_variants, keys %{ $self->{list_index}->{$chr} };
	}
	return \@all_variants;
}

sub rocksdb1 {
	my ($self,$chr) = @_;
	
	confess() if $self->mode eq "temp";
	return $self->{db}->{$chr} if exists $self->{db}->{$chr};
	if ($self->mode eq "c") {
		$self->{db}->{$chr} =GenBoNoSqlRocks->new(dir=>$self->dir,mode=>"c",name=>$chr->name,compression=>"snappy",pipeline=>1);
	}
	else {
		$self->{db}->{$chr} = GenBoNoSqlRocks->new(dir=>$self->dir,mode=>$self->mode,name=>$chr,vmtouch=>1);
	}
	return $self->{db}->{$chr};
}

sub rocksdb {
	my ($self,$chr,$id) = @_;
	confess() unless $id;
	confess() if $self->mode eq "temp";
	return $self->{db}->{$chr}->{$id} if exists $self->{db}->{$chr}->{$id};
	my $dir = $self->dir."/".$id;
	if ($self->mode eq "c") {
		
		$self->{db}->{$chr}->{$id} =GenBoNoSqlRocks->new(dir=>$dir,mode=>"c",name=>$chr->name,compression=>"snappy",pipeline=>1);
	}
	else {
		$self->{db}->{$chr}->{$id} = GenBoNoSqlRocks->new(dir=>$dir,mode=>$self->mode,name=>$chr,vmtouch=>1);
	}
	return $self->{db}->{$chr}->{$id};
}



 
sub columns  {
	my ($self,$value) =@ _;
	return $self->{columns}->{$value} if exists $self->{columns}->{$value};
	my $file_freeze = $self->columns_file();
	if ($self->mode eq "c"){
		$self->{columns}->{genes} = $self->columns_transcripts;
		$self->{columns}->{variants} = $self->columns_global;
		$self->{columns}->{patients} = $self->columns_patients;
		$self->{columns}->{patients_sv} = $self->columns_patients_sv;
		store ($self->{columns}, $file_freeze) or die "Impossible de sauver: $!";;
	}
	else {
	 	die() unless -e $file_freeze;
 	 	$self->{columns}  = retrieve($file_freeze);
 		return $self->{columns}->{$value};
 	}
 	
}

sub cache_polyviewer_variant {
	my ($self,$chr,$ids) =@_;
	$self->rocksdb($chr,$self->patient->getFamily->id)->prepare($ids);
	foreach my $id (@$ids){
		$self->{rcache}->{$id} = $self->get_polyviewer_variant($id);
		#$self->{print_cache}->{$id} = $self->print_line_variant($self->{rcache}->{$id});
	} 
}

sub clean {
	my ($self) = @_;
	return delete $self->{rcache};
}

sub close {
	my ($self) = @_;
	foreach my $chr (keys %{$self->{db}}){
		$self->_close($chr);
	}
	
}
sub return_keys {
	my ($self) = @_;
	return [keys %{$self->{rcache}}];
}
sub get_cache {
	my ($self) = @_;
	return $self->{rcache}
}
sub set_cache {
	my ($self,$hash) = @_;
	$self->{rcache} = $hash;
}

sub _close {
	my ($self,$chr) =@_;
	foreach my $id (keys %{$self->{db}->{$chr}}){
		 $self->{db}->{$chr}->{$id}->close();
		 delete $self->{db}->{$chr}->{$id};
	}
		
	delete $self->{db}->{$chr};
	return 1;
}

sub get_polyviewer_variant {
	my ($self,$id,$debug) =@_;
	return $self->{rcache}->{$id} if exists $self->{rcache}->{$id};
	#warn "00"  if $debug;;
	#return if $debug;
	my $test = "1!75318";
	my ($chr,$a)  =split("!",$id);
	my $v = $self->rocksdb($chr,$self->patient->getFamily->id)->get($id);
	die($id) unless $v;
	my $vp = PolyviewerVariant->new() ;
     	foreach my $c (@{$self->columns("variants")}){
     		$vp->{$c} = shift @{$v->[0]};
     	}
     	
     
     	#1-215853720-T-C
     	my $genes = $v->[1];
     	foreach my $g (@$genes){
     		
     		my $gene_id= shift @$g;
     		my $ht ;
     		foreach my $c (@{$self->columns("genes")}){
     			$ht->{$c} = shift @{$g};
     		}
     		push(@{$vp->{hgenes}->{$gene_id}->{tr}},$ht);
     	}
     #	$debug =1 if  $vp->{name} eq "1-1392838-del-29892";
     	my $patients = $v->[2];
     	my $key_patient = "patients";
     	$key_patient = "patients_sv" if $vp->{isSrPr};
     	foreach my $p  (@$patients) {
     		my $fam = shift @$p;
     		my $patient_id = shift @$p;
     		
     	foreach my $c (@{$self->columns($key_patient)}){
			#warn $c if $vp->{isSrPr};
			#warn Dumper @{$p} if $vp->{isSrPr};
     			$vp->{patients_calling}->{$patient_id}->{$c} = shift @{$p};
     	  }
     	  }
     	  $self->update_clinvar_id($chr,$vp) if $vp->{clinvar};
     	  warn Dumper $vp unless $vp->name;
     	  warn $id  unless $vp->name;;
     	  die() unless $vp->name;
     	  
     	return $vp;
}
#
sub update_clinvar_id {
	my ($self,$chr,$vp) = @_;
	my $pub = $self->clinvar_rocks_db($chr)->clinvar($vp->rocksdb_id);
	$vp->{clinvar_id} = $pub->{clinvar_id};
}
sub clinvar_rocks_db {
	my ($self,$chr,$vp) = @_;
	return $self->{"clinvarrocks".$chr} if exists  $self->{"clinvarrocks".$chr};
	my $chromosome = $self->project->getChromosome($chr);
	$self->{"clinvarrocks".$chr} =   $chromosome->rocksdb("clinvar");
	 return $self->{"clinvarrocks".$chr};
	
}
sub load_polyviewer_variant {
	my ($self) =@_;
	my @chr = sort {$a cmp $b} keys %{$self->{list_index}};
 	$self->{rcache} = {};
 	my $error;
 	my $jobs;
 	my $shared_hash = MCE::Shared->hash();
 	MCE::Loop::init {
    chunk_size => 'auto',
    max_workers => 'auto',
    gather => sub {
        my ($mce,$data,$chrs) = @_;
        delete $shared_hash->{$mce};
       		foreach my $c (@$chrs){
				$jobs->{$c} ++;
			}
        	@{$self->{rcache}}{keys %$data} = values %$data;
    	},
    	   on_post_exit => sub {
        my ($mce, $pid, $exit_code, $ident) = @_;
        if ($exit_code != 0) {
            warn "?? Worker $pid (ident=$ident) exited with error $exit_code\n";
            $error ++;
        } else {
            print "? Worker $pid (ident=$ident) exited normally\n";
        }
    }
	};

	 mce_loop {
  	  my ($mce, $chra) = @_;
  	  my $chr = $chra->[0];
  	  my $x;
  		my $hash_vp;
  		warn "start ".MCE->chunk_id;
  		$shared_hash->{MCE->chunk_id} ++;
  	  foreach my $chr (@$chra){
    	 my $nbv =0;
 	   	$self->cache_polyviewer_variant($chr,[keys %{$self->{list_index}->{$chr}}]);
  	  }
  	   MCE->gather(MCE->chunk_id,$self->clean,$chra);
   	  $self->close();
	
		} @chr;
 		MCE::Loop->finish;
 		confess()  if %$shared_hash;
 		confess("Argh,  Something went wrong !!! ") if $error;
 		confess("Argh,  Something went wrong !!! ") if scalar (keys %$jobs) != scalar(@chr);
}

sub transform_polyviewer_variant {
	my ($self,$chr,$index,$vp,$sereal) = @_ ;
	die() unless $vp;
	 my $array;
	 my $res1;
       	foreach my $k (@{$self->columns_global}) {
				my $value = $vp->{$k};
				my $vp;
				push(@$res1,$value);
				
		}
		$array->[0] = $res1;
		#1-153043318-del-22892
        foreach my $g (keys %{$vp->{hgenes}}){
			
			foreach my $t (@{$vp->{hgenes}->{$g}->{tr}}){
				my $res ;
				$t->{rocksdb_id} = $chr->name."!".$vp->{rocksdb_id};
				$t->{gene} = $g;
				$t->{chromosome} = $chr->ucsc_name;
				$t->{vector_id} = $index;
				push(@$res,$g);
				foreach my $k (@{$self->columns_transcripts}){
					my $value = $t->{$k};
					$value ="-" unless $value;
					if ($value eq "-"){
						$value = -1 if ($htr->{$k} ne "VARCHAR" );
					}
					
					push(@$res,$value);
				}
				push(@{$array->[1]},$res);
			}
		}
		my $array_c = $self->columns_patients;
		$array_c = $self->columns_patients_sv if $vp->{isSrPr};
		my $debug;
		$debug = 1  if  $vp->name eq "11-68282943-del-128811";
	  	warn "OK" if  $vp->name eq "1-144896379-del-10338";
	  	my $fam ={};
		foreach my $p (keys %{$vp->{patients_calling}}) {
			
			my $hash = $vp->{patients_calling}->{$p};
			my $patient= $chr->project->getPatient($p);
			#my $name =  $vp->{patients_calling}->{$p}->{name};
			$hash->{name} = $patient->name;
			#next unless $name;
		#	my $patient= $project->getPatient($name);
			my $f = $patient->getFamily->id;
			my $res;
			push(@{$res},$patient->getFamily->name);
			push(@{$res},$p);
			
			if (defined $vp->{isSrPr} && !(exists $vp->{level})) {
			my $a = $patient->mean_normalize_depth($chr->name,$vp->start,$vp->end);
			$hash->{norm_depth} = $a;#int(sum(@$a)/scalar(@$a));
			$a = $patient->mean_normalize_depth($chr->name,$vp->start-5050,$vp->start-50);
			$hash->{norm_depth_before} = $a;# int(sum(@$a)/scalar(@$a));
			$a = $patient->mean_normalize_depth($chr->name,$vp->end+50,$vp->end+5050);
			$hash->{norm_depth_after} =  $a;#int(sum(@$a)/scalar(@$a));
			my $m =  int(($hash->{norm_depth_before} + $hash->{norm_depth_after})/2);
			$hash->{norm_depth} += 0.001;
			$m += 0.001;
			 my $v = log($hash->{norm_depth}/$m)/log(2);
	 		$v =-2 if $v < -2; 
    		$hash->{log2_ratio} =    sprintf("%.2f",$v);
			}
			foreach my $k (@$array_c){
				my $value = $hash->{$k};
				push(@{$res},$value);
			}
			push(@{$fam->{$f}},$res);
			#push(@{$array->[2]},$res);
		}
		#return $array unless $sereal;
		my $hf = {};
		foreach my $f (keys %{$fam}){
			$hf->{$f} =  $self->encode([$array->[0],$array->[1],$fam->{$f}]);
		} 
		return $hf;
		#return $self->encode($array);
}

sub encode {
	my ( $self, $value ) = @_;
	return sereal_encode_with_object( $self->sereal_encoder, $value );
}

has sereal_encoder => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;

#return Sereal::Encoder->new({compress=>Sereal::SRL_ZSTD,compress_threshold=>0});
		return Sereal::Encoder->new(
			{ compress => Sereal::SRL_ZSTD, compress_threshold => 0 } );
		return 0;
	},
);


sub print_line_variant {
	my ($self,$vp,$opacity) = @_;
		my $print_html = $self->print_html();
		 my $cgi = $print_html->cgi;
			my @headers = (
				"mobidetails","varsome", "igv",    "alamut", "var_name",
				"trio",    "gnomad", "deja_vu"
			);
			my $style = {};
			$style = { style => "background-color: #DAEEED;opacity:0.5" }  if $opacity;

	my $out;

		
			my $hpatients;
		
		
			#my $is_gnomad = exists $v->{value}->{ac};
			

			
			$print_html->variant($vp);

			my $icon = qq{<img width="32" height="32" src="https://img.icons8.com/external-gliphyline-royyan-wijaya/32/external-laptop-laptop-collection-glyphyline-gliphyline-royyan-wijaya-15.png" alt="external-laptop-laptop-collection-glyphyline-gliphyline-royyan-wijaya-15"/>};
			 $icon   = qq{<img width="24" height="24" src="https://img.icons8.com/external-tal-revivo-filled-tal-revivo/24/external-live-preview-of-a-smart-class-education-school-filled-tal-revivo.png" alt="external-live-preview-of-a-smart-class-education-school-filled-tal-revivo"/>};
			#$icon = qq{<img width="28" height="28" src="https://img.icons8.com/external-bearicons-blue-bearicons/28/external-Clipboard-clipboards-and-notepads-bearicons-blue-bearicons-30.png" alt="external-Clipboard-clipboards-and-notepads-bearicons-blue-bearicons-30"/>};
			
			my $dropdown = qq{
			<div class="dropdown">
  <button class="btn btn-primary btn-xs dropdown-toggle " type="button" id="dropdownMenuButton" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false" style="font-size:10px;background-color:#C67FAE">
  $icon
  </button>
  <div class="dropdown-menu" aria-labelledby="dropdownMenuButton" style="font-size:12px;background-color:beige;color:black">
  }.
  "<li>".$print_html->mobidetails()."</li>".
  	"<li>".$print_html->gnomadurl()."</li>".
			"<li>".$print_html->alamuturl()."</li>".
			"<li>".$print_html->varsome()."</li>"
  
  .qq{</div>
</div>
};

			my $t1 = shift(@headers);
			$out .= $cgi->td( $style, $dropdown );
		
			
			
			$out .= "\n";
			##############
			# IGV CELL
			###############

			my $t = shift(@headers);

			#write locus
			$out .= $cgi->td( $style, $print_html->igv);

			##############
			
			
			##############
			# NAME CELL
			#
			###############
			$t = shift(@headers);

			#$name =  $v->{var_name} if exists $v->{var_name};
		
			$out .= $cgi->td($style,$print_html->var_name());
			$out .= "\n";
			##############
			# CELL CALLING INFOS
			###############

	
			$out .= $cgi->td( $style, $print_html->calling()) ;

			$out .= "\n";
			
			return;
		
			
			$t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->gnomad() );
			$out .= "\n";


			$t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->dejavu() );
			$out .= "\n";
			$t = shift(@headers);
			$out .= $cgi->td( $style, $print_html->validations );
			
			$t = shift(@headers);
		
			$out .= "\n";
			$out .= $cgi->td( $style, $print_html->transcripts() );
			$out .= "\n";

			$out .= $cgi->td( $style, $print_html->validation_select() );
			$out .= $cgi->end_Tr();
			$out .= "\n";
}


1;
