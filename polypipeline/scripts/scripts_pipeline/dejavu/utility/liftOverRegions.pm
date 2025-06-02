package liftOverRegions;
use strict;
use warnings;
use IPC::Run3;
use Exporter 'import';
use Data::Dumper;
use Moo;
use JSON::XS;
# Déclare les fonctions exportées
our @EXPORT_OK = qw(lift_over_variants);
#####


has project => (
	is		=> 'ro',
	required=> 1,
);

has version => (
	is		=> 'ro',
	required=> 1,
);
has regions => (
	is		=> 'ro',
	default => sub {
		return ;
	},
);

has file_regions => (
	is		=> 'rw',
	#isa		=> 'Str',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $fht2 = File::Temp->new(
    	DIR    =>  $self->tmp_dir,     # Optionnel : répertoire pour le fichier
    	SUFFIX => '.bed'      # Optionnel : suffixe du fichier
	);
	return $fht2;
	},
);
has tmp_dir => (
	is		=> 'rw',
	#isa		=> 'Str',
	lazy	=> 1,
	default	=> sub {
		return "/tmp";
	},
);

has fh_write_regions  => (
	is		=> 'rw',
	#isa		=> 'Str',
	lazy	=> 1,
	default	=> sub {
	my $self = shift;	
	 my $fh1;
	open  $fh1, '>', $self->file_regions->filename or die "Impossible d'ouvrir le fichier temporaire: $!\n";
	return $fh1;
	},
);
has chain_file  => (
	is		=> 'ro',
	#isa		=> 'Str',
	lazy	=> 1,
	default	=> sub {
	my $self = shift(@_);
	 return $self->project->liftover_chain_file($self->version);
	},
);

has hash_regions  => (
	is		=> 'ro',
	#isa		=> 'Str',
	lazy	=> 1,
	default	=> sub {
	return {};
	},
);
	
sub add_region {
	  my ($self,$region) = @_;
	  my $fh = $self->fh_write_regions;
	  my $id = $region->{id};
	  $id =  $region->{rocksid} unless $id;
	  confess() unless $id;
	  $self->hash_regions->{$id} =  $region;
	  print  $fh  join("\t",$region->{chromosome},$region->{start},$region->{end},encode_json($region))."\n";
}

sub add_region_id  {
	  my ($self,$region) = @_;
	  my $fh = $self->fh_write_regions;
	  my $id = $region->{id};
	  confess() unless $id;
	  die() if exists $self->hash_regions->{$id};
	  $self->hash_regions->{$id} =  $region;
	  print  $fh  join("\t",$region->{chromosome},$region->{start},$region->{end},$id)."\n";
}

sub _liftOver_regions{
	 my ($self,$name,$opt) = @_;
	 $opt = "" unless $opt;
	   close $self->fh_write_regions;
	   delete $self->{fh_write_regions};
	   
		my $fht2 = File::Temp->new(
    	DIR    => $self->tmp_dir,     # Optionnel : répertoire pour le fichierliftOver_regions
    	SUFFIX => '.bed'      # Optionnel : suffixe du fichier
		);
		
		my $fileout = $fht2->filename;
		if ($name){
			$fileout =$self->tmp_dir."/".$name.time.".tmp";
		}
		my $fht21 = File::Temp->new(
    	DIR    => $self->tmp_dir,     # Optionnel : répertoire pour le fichier
    	SUFFIX => '.bed'      # Optionnel : suffixe du fichier
		);
		
		my $fileoutsort = $fht21->filename;
		if ($name){
			$fileoutsort =$self->tmp_dir."/".$name.time.".sort.tmp";
		}
	   my $cmd = $self->project->buffer->software("liftOver")." ".$self->file_regions->filename." ". $self->chain_file." $opt ".$fileout." /dev/stderr >/dev/null 2>/dev/null";
	   system($cmd." && sort -k1,1V -k2,2n $fileout > $fileoutsort && rm $fileout"  );
	   return $fileoutsort;
}

sub liftOver_regions {
	   my ($self,$name) = @_;
		my $fileoutsort = $self->_liftOver_regions($name);
		
	   return $self->parse_bed_region($fileoutsort);
}

sub parse_bed_region {
	my ($self,$fileout) = @_;
	open(my $fh, '<', $fileout) or die "Impossible d'ouvrir le fichier '$fileout' : $!";
	# Lis le fichier ligne par ligne
	my $tab = {};
	my $hv;
	while (my $line = <$fh>) {
		
		chomp $line;  # Supprime le caractère de fin de ligne (\n)
   	 	next if $line =~/^#/;
   	 	my @t = split("\t",$line);
   	 	my $h = decode_json($t[3]);
   	 	$h->{LIFT}->{start} =  $t[1];
   	 	$h->{LIFT}->{end} =  $t[2];
   	 	$h->{LIFT}->{chromosome} =  $t[0];
   	 	push(@{$tab->{$t[0]}},$h);
	}
	close($fh);
	unlink $fileout;
	return $tab;
}



sub liftOver_regions_cnv {
	   my ($self,$name) = @_;
		my $fileoutsort = $self->_liftOver_regions($name,"-multiple -minMatch=0.95 ");
	   return $self->parse_bed_cnv($fileoutsort);
}

sub parse_bed_cnv {
	my ($self,$fileout) = @_;
	open(my $fh, '<', $fileout) or die "Impossible d'ouvrir le fichier '$fileout' : $!";
	# Lis le fichier ligne par ligne
	my $tab = {};
	my $hv;
	while (my $line = <$fh>) {
		chomp $line;  # Supprime le caractère de fin de ligne (\n)
   	 	next if $line =~/^#/;
   	 	my @t = split("\t",$line);
   	 	my $id = $t[3];
   	 	die() unless exists $self->hash_regions->{$id};
   	 	my $sv = $self->hash_regions->{$id};
   	 	$self->hash_regions->{$id}->{LIFT}->{MULTI} ++;
   	 	$self->hash_regions->{$id}->{LIFT}->{start} =  $t[1];
   	 	$self->hash_regions->{$id}->{LIFT}->{end} =  $t[2];
   	 	$self->hash_regions->{$id}->{LIFT}->{chromosome} =  $t[0];
	}
	close($fh);
	unlink $fileout;
	return $self->hash_regions
}

sub add_region_bnd  {
	  my ($self,$region) = @_;
	  my $fh = $self->fh_write_regions;
	  my $id = $region->{id};
	  confess() unless $id;
	  die() if exists $self->hash_regions->{$id};
	  $self->hash_regions->{$id} =  $region;
	  print  $fh  join("\t",$region->{chromosome1},$region->{position1},$region->{position1},$id)."\n";
	  print  $fh  join("\t",$region->{chromosome2},$region->{position2},$region->{position2},$id)."\n";
}

sub liftOver_bnd {
	   my ($self,$name) = @_;
	   my $fileoutsort = $self->_liftOver_regions($name);
	   return $self->parse_bed_bnd($fileoutsort);
}


sub parse_bed_bnd {
	my ($self,$fileout) = @_;
	open(my $fh, '<', $fileout) or die "Impossible d'ouvrir le fichier '$fileout' : $!";
	# Lis le fichier ligne par ligne
	my $tab = {};
	my $vh;
	while (my $line = <$fh>) {
		chomp $line;  # Supprime le caractère de fin de ligne (\n)
   	 	next if $line =~/^#/;
   	 	my @t = split("\t",$line);
   	 	my $id = $t[3];
   	 	die() unless exists $self->hash_regions->{$id};
   	 	my $sv = $self->hash_regions->{$id};
   	 	push(@{$self->hash_regions->{$id}->{LIFT}},{position=>$t[1],chromosome=>$t[0]});
   	 	$self->hash_regions->{$id}->{NB_LIFT} ++;
	}
	close($fh);

	unlink $fileout;
return $self->hash_regions;
}

sub add_variant {
	 my ($self,$v) = @_;
	  my $fh = $self->fh_write_regions;
	  my $vcf_id = $v->theoric_vcf_id();
	  print  $fh, join("\t",$v->getChromosome->ucsc_name,$vcf_id->{position},($vcf_id->{position}+1),$v->id);
}

sub lift_over_variant {
	   my ($self,$variation,$version,$key) = @_;
	   my $id = $variation->id;
	   my $cmd = 
		$self->lift_over_variants([$variation],$key);
}


# Fonction exportable pour effectuer le lift-over d'un variant
sub lift_over_variants {
    my ($self,$variations,$key,$db) = @_;
    my $version = $self->version;
    my $project= $self->project;
	 $key = "lift_over_".$self->version unless $key;
    # Crée une entrée au format VCF
    my $vcf_input = "##fileformat=VCFv4.2\n";
    $vcf_input   .= "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    my $bed_input="";
    my @string;
    my $vpos;
    my $res = {};
    foreach my $v (@$variations){
		my ($c,$pos,$allele_ref,$allele_var) = split("_",$v->vcf_id);
		
		my $vcf_id = $v->theoric_vcf_id();
		if ($v->getChromosome->name eq "MT") {
  			$res->{$v->id}->{chromosome} = $v->getChromosome->ucsc_name;
   	 		$res->{$v->id}->{position_vcf} = $v->start;
   	 		next;
  		}
  		push(@string,join("\t", $v->getChromosome->ucsc_name,$vcf_id->{position},($vcf_id->{position}+1),$v->id ));
    	#push(@string, join("\t", $v->getChromosome->ucsc_name,$vcf_id->{position}, '.', $vcf_id->{ref},$vcf_id->{alt}, '.', '.',$v->id ));
    	$vpos->{$v->id} = length($allele_ref) - length($allele_var);
    	
	}
	$vcf_input .= join("\n",@string);
	$bed_input .= join("\n",@string);
	
  # 	my $res =  run_crossmap ($project,$vcf_input,$version);
  

   #run_crossmap($project,$vcf_input,$version,$res);
   run_liftOver($project,$bed_input,$res);
   foreach my $v (@$variations){
   	my $id = $v->id;
   	die() unless exists $res->{$id};
   	$v->{$key} = delete $res->{$id};
   	my $pvcf = $v->{$key}->{position_vcf} ;
	$v->{$key}->{position} = $v->{$key}->{position_vcf} ;
	$v->{$key}->{position} +=  1 unless $v->isVariation;
	my ($c,$s,$a,$b) = split("_",$v->vcf_id);
	$v->{$key}->{vcf_id} = join("_",$v->{$key}->{chromosome}, $v->{$key}->{position_vcf},$a,$b);
	 ($c,$s,$a,$b) = split("_",$v->id);
	my $chr = $v->project->getChromosome($v->{$key}->{chromosome});
	$v->{$key}->{id} = join("_",$chr->name, $v->{$key}->{position},$a,$b);
	($c,$s,$a,$b) = split("-",$v->name);
	$v->{$key}->{name} = join("-",$chr->name,$v->{$key}->{position},$a,$b);
   }
}



sub run_liftOver {
	my ($self,$bed,$res) = @_;
	my $version = $self->version;
	my $project = $self->project;
	my $chain_file = $self->project->liftover_chain_file($version);
	die "Fichier de chaîne manquant !" unless -e $chain_file;
	
	my $fht = File::Temp->new(
    	DIR    => "/tmp",     # Optionnel : répertoire pour le fichier
    	SUFFIX => '.bed'      # Optionnel : suffixe du fichier
	);

	# Récupérer le nom du fichier
	my $fileout = $fht->filename;
		my $fht2 = File::Temp->new(
    	DIR    => "/tmp",     # Optionnel : répertoire pour le fichier
    	SUFFIX => '.bed'      # Optionnel : suffixe du fichier
	);

	# Récupérer le nom du fichier
	my $fileout2 = $fht->filename;
    # Vérifie si le fichier de chaîne existe
    # Prépare la commande CrossMap
    my $vg = "HG38_DRAGEN";
    $vg = "HG19_MT" if $version eq "HG19";
    my $fasta  = $project->buffer()->config_path("public_data") . "/genome/"
		  . $vg . "/fasta/all.fa";
		  
	my @cmd = (
    	$project->buffer->software("liftOver"), "/dev/stdin", $chain_file,$fileout,"/dev/stderr"
	);
	
	my $stdout;
	my $stderr;
	run3 \@cmd, \$bed, \$stdout, \$stderr;
 	parse_bed($fileout,$res);
 return 1;
}



sub parse_bed {
	my ($fileout,$res) = @_;
	open(my $fh, '<', $fileout) or die "Impossible d'ouvrir le fichier '$fileout' : $!";

	# Lis le fichier ligne par ligne
	while (my $line = <$fh>) {
		chomp $line;  # Supprime le caractère de fin de ligne (\n)
   	 	next if $line =~/^#/;
   	 	my @t = split("\t",$line);
   	 	
   	 	my $id = $t[3];
   	 	$res->{$id}->{chromosome} = $t[0];
   	 	$res->{$id}->{chromosome} =~ s/chrM/chrMT/;
   	 	$res->{$id}->{chromosome} =~ s/chr//;
   	 	$res->{$id}->{position_vcf} = $t[1];
	}

# Ferme le fichier
close($fh);
}

1;
