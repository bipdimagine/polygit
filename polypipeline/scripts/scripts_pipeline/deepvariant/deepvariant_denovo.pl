#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use FindBin qw($Bin);
#use lib "/software/polyweb/poly-disk/poly-src/polygit/GenBo/lib/obj-nodb/";
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages/";
use  File::Temp;
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use colored;
use Storable qw(store retrieve freeze);


my $project_name;
my $final_vcf;
my $log_file;
my $patient_name;
my $fork;

GetOptions(
	'project=s'   => \$project_name,
	"log=s" =>\$log_file,
	"vcf=s" => \$final_vcf,
	"patient=s"=>\$patient_name,
	"fork=s" =>\$fork,
);


my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );

foreach my $f (@{$project->getFamilies}){
	my $parents = $f->getParents;
	next unless @$parents;
	my $variants_child = get_variants($f->getChildren);
	change_vcf_parent($f->getMother,$variants_child) if $f->getMother();
	change_vcf_parent($f->getFather,$variants_child) if $f->getFather();
	
}
sub change_vcf_parent {
	my ($p,$variants) = @_;
	my $vcf = $p->getVariationsFile("deepvariant");
	my $dir_vcf_out= $project->getCallingPipelineDir("deepvariant");
	my $vcf_out = $dir_vcf_out."/".$p->name.".vcf";
	
	
    warn "-> in : $vcf \n";
    warn  "    $vcf_out...\n";
my $nb_change =0;
    my $fh_in;
    if ($vcf =~ /\.gz$/) {
        open($fh_in, "gzip -dc $vcf |") or die "Impossible d'ouvrir $vcf : $!\n";
    } else {
        open($fh_in, '<', $$vcf) or die "Impossible d'ouvrir $vcf : $!\n";
    }
    
    open(my $fh_out, '>', $vcf_out) or die "Impossible de créer $vcf_out : $!\n";

    while (my $ligne = <$fh_in>) {
        chomp $ligne;
        
        # Si c'est une ligne d'en-tête, on l'écrit telle quelle
        if ($ligne =~ /^#/) {
            print $fh_out "$ligne\n";
            next;
        }
		if ($ligne !~ /RefCall/) {
            print $fh_out "$ligne\n";
            next;
        }
        my @colonnes = split(/\t/, $ligne);
        my $chrom = $colonnes[0];
        my $pos   = $colonnes[1];
        my $ref   = $colonnes[3];
        my $alt   = $colonnes[4];

        my $cle = "$chrom:$pos:$ref:$alt";
        unless (exists $variants->{$cle}){
        	  print $fh_out "$ligne\n";
           	next;
        }
        # Si ce variant existe chez l'enfant
            
           my $format = $colonnes[8];
            
            if (defined $format) {
                my @champs_format = split(/:/, $format);
                my $index_gt = -1;
                my $index_ad = -1; # Ajout de la recherche de l'index AD
                for my $i (0 .. $#champs_format) {
                    $index_gt = $i if $champs_format[$i] eq 'GT';
                    $index_ad = $i if $champs_format[$i] eq 'AD';
                }

                if ($index_gt != -1) {
                    for my $s (9 .. $#colonnes) {
                        my @champs_echantillon = split(/:/, $colonnes[$s]);
                        my $gt = $champs_echantillon[$index_gt];
                        
                        # Vérifier que c'est bien 0/0 ou 0|0
                        if (defined $gt && ($gt eq '0/0' || $gt eq '0|0' || $gt eq '1/1' || $gt eq '1|1')) {
                            
                            my $modification_autorisee = 0;

                            # Si le champ AD est présent
                            if ($index_ad != -1 && defined $champs_echantillon[$index_ad] && $champs_echantillon[$index_ad] ne '.' ) {
                                my $ad = $champs_echantillon[$index_ad];
                                my @profondeurs = split(/,/, $ad);
                                
                                # profondeurs[0] = reads REF, profondeurs[1] = reads ALT
                                my $reads_alt = $profondeurs[1];
                                my $reads_ref = $profondeurs[0];
                                # VERIFICATION : Nombre de reads ALT > 1
                                if (defined $reads_alt && $reads_alt ne '.' && $reads_alt > 1) {
                                    $modification_autorisee = 1;
                                    $modification_autorisee = 2 if (($reads_alt/($reads_ref+$reads_alt)) > 0.8);
                                }
                            }

                            # Si les conditions sont remplies, on transforme en 0/1
                            if ($modification_autorisee) {
                                $champs_echantillon[$index_gt] = '0/1';
                                $champs_echantillon[$index_gt] = '1/1' if $modification_autorisee == 2 ;
                                $colonnes[$s] = join(":", @champs_echantillon);
                                $nb_change++;
                                  $colonnes[6] ="RefCall_dnv";
                                 $ligne = join("\t", @colonnes);
                                
                            }
                        }
                   }
			}
        }
        # Écrire la ligne (modifiée ou non)
      
     print $fh_out "$ligne\n";
    }
    close($fh_in);
    close($fh_out);
    system("bgzip -f $vcf_out ; cp ${vcf_out}.gz $vcf;tabix -f -p vcf $vcf");
    warn "  end \n\n";
    warn $vcf_out;
	
	
}




sub get_variants {
	my ($children) = @_;
	my $hash;
	foreach my $p (@$children){
		my $fh_enfant;
		my $vcf_enfant = $p->getVariationsFile("deepvariant");
		if ($vcf_enfant =~ /\.gz$/) {
    		open($fh_enfant, "gzip -dc $vcf_enfant |") or die "Impossible d'ouvrir $vcf_enfant : $!\n";
		} else {
    		open($fh_enfant, '<', $vcf_enfant) or die "Impossible d'ouvrir $vcf_enfant : $!\n";	
		}
		while (my $ligne = <$fh_enfant>) {
    chomp $ligne;
    next if $ligne =~ /^#/; # Ignorer l'en-tête

    my @colonnes = split(/\t/, $ligne);
    my $chrom = $colonnes[0];
    my $pos   = $colonnes[1];
    my $ref   = $colonnes[3];
    my $alt   = $colonnes[4];

    # Création d'une clé unique pour chaque variant
    my $cle = "$chrom:$pos:$ref:$alt";
    $hash->{$cle} = 1;
}
close($fh_enfant);
}
return $hash;
}
