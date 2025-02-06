package liftOver;
use strict;
use warnings;
use IPC::Run3;
use Exporter 'import';
use Data::Dumper;
use Bio::DB::Fasta;
# Déclare les fonctions exportées
our @EXPORT_OK = qw(lift_over_variants);
#####

# Fonction exportable pour effectuer le lift-over d'un variant
sub lift_over_variants {
    my ($project,$variations,$version,$key) = @_;
	 $key = "lift_over_$version" unless $key;
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
   run_liftOver($project,$bed_input,$version,$res);
   foreach my $v (@$variations){
	   	my $id = $v->id;
	   	next unless exists $res->{$id};
	   	die() unless exists $res->{$id};
	   	$v->{$key} = delete $res->{$id};
	   	my $pvcf = $v->{$key}->{position_vcf} ;
		$v->{$key}->{position} = $v->{$key}->{position_vcf} ;
		$v->{$key}->{position} +=  1 unless $v->isVariation;
		my ($c,$s,$a,$b) = split("_",$v->vcf_id);
		$v->{$key}->{vcf_id} = join("_",$v->{$key}->{chromosome},$v->{$key}->{position_vcf},$a,$b);
		($c,$s,$a,$b) = split("_",$v->id);
		if (not $v->{$key}->{chromosome} =~ /^[0-9XYMT]+$/) {
			$v->{$key}->{id} = join("_",$v->{$key}->{chromosome},$v->{$key}->{position},$a,$b);
			($c,$s,$a,$b) = split("-",$v->name);
			$v->{$key}->{name} = join("-",$v->{$key}->{chromosome},$v->{$key}->{position},$a,$b);
		} 
		else {
			my $chr = $v->project->getChromosome($v->{$key}->{chromosome});
			$v->{$key}->{id} = join("_",$chr->name, $v->{$key}->{position},$a,$b);
			($c,$s,$a,$b) = split("-",$v->name);
			$v->{$key}->{name} = join("-",$chr->name,$v->{$key}->{position},$a,$b);
		}
   }
}

sub lift_over_variant {
	   my ($variation,$version,$key) = @_;
	   my $id = $variation->id;
		lift_over_variants($variation->project,[$variation],$version,$key);
}
sub run_liftOver {
	my ($project,$bed,$version,$res) = @_;
	my $chain_file = $project->liftover_chain_file($version);
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
	#warn join(" ",@cmd);
	run3 \@cmd, \$bed, \$stdout, \$stderr;
	my $db = Bio::DB::Fasta->new($fasta);
 parse_bed($fileout,$res,$db);
 return 1;
}
sub run_crossmap {
	my ($project,$vcf,$version,$res) = @_;
	my $chain_file = $project->liftover_chain_file($version);
	die "Fichier de chaîne manquant !" unless -e $chain_file;
	
	my $fht = File::Temp->new(
    	DIR    => "/tmp",     # Optionnel : répertoire pour le fichier
    	SUFFIX => '.vcf'      # Optionnel : suffixe du fichier
	);

	# Récupérer le nom du fichier
	my $fileout = $fht->filename;
    # Vérifie si le fichier de chaîne existe
    # Prépare la commande CrossMap
    my $vg = "HG38_DRAGEN";
    $vg = "HG19_MT" if $version eq "HG19";
    my $fasta  = $project->buffer()->config->{'public_data'}->{root} . "/genome/"
		  . $vg . "/fasta/all.fa";
		  
	my @cmd = (
    $project->buffer()->software("singularity"), 'run',
    '-B', $project->buffer()->config->{'singularity'}->{mount_data},
    $project->buffer()->config->{'singularity'}->{dir}.$project->buffer()->config->{'singularity'}->{crossmap},
    'CrossMap.py', 'vcf',
    $chain_file,
    '-',   # Lecture de VCF depuis stdin
    $fasta,
    $fileout    # Écriture de VCF vers stdout
	);
	my $stdout;
	my $stderr;
	run3 \@cmd, \$vcf, \$stdout, \$stderr;
	warn $stdout;
	warn $stderr;

 parse_vcf($fileout,$res);
}
#sub run_picard {
#	my ($project,$vcf,$version) = @_;
#	my $chain_file = $project->liftover_chain_file($version);
#	die "Fichier de chaîne manquant !" unless -e $chain_file;
#	
#	my $fht = File::Temp->new(
#    	DIR    => "/tmp",     # Optionnel : répertoire pour le fichier
#    	SUFFIX => '.vcf'      # Optionnel : suffixe du fichier
#	);
#
#	# Récupérer le nom du fichier
#	my $fileout = $fht->filename;
#    # Vérifie si le fichier de chaîne existe
#    # Prépare la commande CrossMap
#    my $vg = "HG38_DRAGEN";
#    my $java      = $project->getSoftware('java');
#    $vg = "HG19_MT" if $version eq "HG19";
#    my $fasta  = $project->buffer()->config->{'public_data'}->{root} . "/genome/"
#		  . $vg . "/fasta/all.fa";
#	my $picard = $java . " -jar " . $project->getSoftware('picard_path');
#	my @cmd = (
#	$java , "-jar" , $project->getSoftware('picard_path'),
#    "LiftoverVcf" , 'I=/dev/stdin',
#    "O=$fileout", 'REJECT=/dev/null',
#    "R=$fasta", "CHAIN=$chain_file",
#	);
#	my $stdout;
#	my $stderr;
#	run3 \@cmd, \$vcf, \$stdout, \$stderr;
#	
#	 parse_vcf($fileout);
#	
#}

sub parse_bed {
	my ($fileout,$res,$db) = @_;
	open(my $fh, '<', $fileout) or die "Impossible d'ouvrir le fichier '$fileout' : $!";

	# Lis le fichier ligne par ligne
	while (my $line = <$fh>) {
		#warn $line;
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
sub parse_vcf {
	my ($fileout,$res) = @_;
	open(my $fh, '<', $fileout) or die "Impossible d'ouvrir le fichier '$fileout' : $!";

	# Lis le fichier ligne par ligne
	while (my $line = <$fh>) {
    	chomp $line;  # Supprime le caractère de fin de ligne (\n)
   	 	next if $line =~/^#/;
   	 	my @t = split("\t",$line);
   	 	my $id = $t[7];
   	 	$res->{$id}->{chromosome} = $t[0];
   	 	$res->{$id}->{position_vcf} = $t[1];
   	 	
	}

# Ferme le fichier
close($fh);
}
1;
