package liftOver;
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

	
sub add_region {
	  my ($self,$region) = @_;
	  my $fh = $self->fh_write_regions;
	  print  $fh  join("\t",$region->{chromosome},$region->{start},$region->{end},encode_json($region))."\n";
}


sub liftOver_regions {
	   my ($self) = @_;
	   close $self->fh_write_regions;
	   delete $self->{fh_write_regions};
		my $fht2 = File::Temp->new(
    	DIR    => $self->tmp_dir,     # Optionnel : répertoire pour le fichier
    	SUFFIX => '.bed'      # Optionnel : suffixe du fichier
		);
		
		my $fileout = $fht2->filename;
		my $fht21 = File::Temp->new(
    	DIR    => $self->tmp_dir,     # Optionnel : répertoire pour le fichier
    	SUFFIX => '.bed'      # Optionnel : suffixe du fichier
		);
		
		my $fileoutsort = $fht21->filename;
		
	   my $cmd = $self->project->buffer->software("liftOver")." ".$self->file_regions->filename." ". $self->chain_file." ".$fileout." /dev/stderr >/dev/null 2>/dev/null";
	   warn $cmd." && sort -k1,1V -k2,2n $fileout > $fileoutsort && rm $fileout";
	   system($cmd." && sort -k1,1V -k2,2n $fileout > $fileoutsort && rm $fileout"  );
	   warn "end lift";
	   return $self->parse_bed_region($fileoutsort);
}


sub parse_bed_region {
	my ($self,$fileout) = @_;
	open(my $fh, '<', $fileout) or die "Impossible d'ouvrir le fichier '$fileout' : $!";

	# Lis le fichier ligne par ligne
	my $tab = {};
	while (my $line = <$fh>) {
		chomp $line;  # Supprime le caractère de fin de ligne (\n)
   	 	next if $line =~/^#/;
   	 	my @t = split("\t",$line);
   	 	my $h = decode_json($t[3]);
   	 	$h->{$self->{version}}->{start} =  $t[1];
   	 	$h->{$self->{version}}->{end} =  $t[2];
   	 	$h->{$self->{version}}->{chromosome} =  $t[0];
   	 	push(@{$tab->{$t[0]}},$h);
	}
	close($fh);
return $tab;
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
    my $fasta  = $project->buffer()->config->{'public_data'}->{root} . "/genome/"
		  . $vg . "/fasta/all.fa";
		  
	my @cmd = (
    	$project->buffer->software("liftOver"), "/dev/stdin", $chain_file,$fileout,"/dev/stderr"
	);
	
	my $stdout;
	my $stderr;
	warn join(" ",@cmd);
	run3 \@cmd, \$bed, \$stdout, \$stderr;
 	parse_bed($fileout,$res);
 return 1;
}



sub parse_bed {
	my ($fileout,$res) = @_;
	open(my $fh, '<', $fileout) or die "Impossible d'ouvrir le fichier '$fileout' : $!";

	# Lis le fichier ligne par ligne
	while (my $line = <$fh>) {
		warn $line;
		chomp $line;  # Supprime le caractère de fin de ligne (\n)
   	 	next if $line =~/^#/;
   	 	my @t = split("\t",$line);
   	 	
   	 	my $id = $t[3];
   	 	$res->{$id}->{chromosome} = $t[0];
   	 	$res->{$id}->{position_vcf} = $t[1];
	}

# Ferme le fichier
close($fh);
}

1;
