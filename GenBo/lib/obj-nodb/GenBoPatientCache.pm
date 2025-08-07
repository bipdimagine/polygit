package GenBoPatientCache;
use strict;
use Storable qw(retrieve thaw);
use Moo;
use Carp;
use Bit::Vector;
use Bit::Vector::Overload;
use Data::Dumper;
use FindBin qw($Bin);
use Storable qw(store retrieve freeze dclone thaw);
extends 'GenBoPatient', 'GenBoCache';

 

# lien avec l'objet GenBoChromosomeCache
has chromosome => (
	is		=> 'rw',
);

# tag si patient in the attic
has in_the_attic => (
	is		=> 'rw',
	default => undef,
);

# tag si patient exclue entierement ou he ou ho
has excluded => (
	is		=> 'rw',
	default => undef,
);

# tag si patient intersect
has intersected => (
	is		=> 'rw',
	default => undef,
);

# ids de tous les genes/intergenic_regions que le patient contient
has hash_list_genes_ids => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub { {} },
);

# nom du model genetique utilise sur ce patient
sub used_model {
	my $self = shift;
	return if (scalar keys %{$self->hash_models_genetics_used()} == 0);
	my @lModels = sort keys (%{$self->hash_models_genetics_used()});
	return join(' + ', @lModels);
}

# table de hash pour l'export en vcf
has xls_var_id => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub { {} },
);

has hash_models_genetics_used => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub { {} },
);

# vecteur temporaire util pour certains modeles genetiques
has model_vector_var_ok => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub { {} },
);

# vector des variants exclus par les differents GenBoPatientCache
has variants_excluded => (
	is      => 'rw',
	lazy    => 1,
	default => sub { {} },
);

has fam_recessif_variants => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->getNewVector();
	}
);

has fam_recessif_categories => (
	is      => 'rw',
	lazy    => 1,
	default => undef,
);

has fam_recessif_global_categories => (
	is      => 'rw',
	lazy    => 1,
	default => undef,
);



has calling_methods_categories => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $h;
		foreach my $chr (@{$self->project->getChromosomes()}) {
			my $no = $chr->get_lmdb_calling_methods("r");
			$h->{$chr->id()} = $no->get($self->name());
			$no->close();
			delete $h->{$chr->id()}->{'index_lmdb'} if (exists $h->{$chr->id()}->{'index_lmdb'});
		}
		return $h;
	},
);

has stats_categories => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->project->buffer->config->{stats_patients};
	},
);



##### METHODS #####



# recupere les anenx d'un variant appartenant a ce patient
sub getVariantAnnex {
	my ($self, $varId) = @_;
	return $self->chromosome->getVariantVcfInfos($varId)->{$self->name()}->{'he_ho_details'};
}

# hash avec les infos necessaires sur ce patient pour un export en XLS
sub store_var_ids {
	my ($self, $chr) = @_;
	if ($self->project->get_xls() eq 'variants' or $self->project->get_xls() eq 'genes') {
		my $hash = $self->xls_var_id();
		my $hTrans;
		my $var_patient = $chr->saved_model->{'for_xls'}->{$self->name()};
		my $var_patient_ho = $chr->saved_model->{'for_xls'}->{$self->name().'_ho'};
		foreach my $var (@{$chr->getListVarObjects($chr->getVariantsVector())}) {
			$self->project->print_dot(25);
			my $id = $var->vector_id();
			next unless ($var_patient->contains($id));
			my $var_id = $var->id();
			my $chr_h_id = $chr->id();
			$chr_h_id = '23' if ($chr->id() eq 'X');
			$chr_h_id = '24' if ($chr->id() eq 'Y');
			$chr_h_id = '25' if ($chr->id() eq 'MT');
			$hash->{$chr_h_id}->{$id}->{'var_id'} = $var_id;
			$hash->{$chr_h_id}->{$id}->{'var_id'} .= ' ('.$var->rs_name().')' if ($var->rs_name());
			$hash->{$chr_h_id}->{$id}->{'type'} = 'snp' if ($var->isVariation());
			$hash->{$chr_h_id}->{$id}->{'type'} = 'ins' if ($var->isInsertion());
			$hash->{$chr_h_id}->{$id}->{'type'} = 'del' if ($var->isDeletion());
			my $h_dejavu = $var->deja_vu();
			my $nb_project = 0;
			my $nb_patient = 0;
			my $he = 0;
			my $ho = 0;
			my $st_project;
			foreach my $projName (keys %$h_dejavu){
				$nb_project++;
				$st_project = $projName.":";
				$st_project .= $h_dejavu->{$projName}->{patients};
				$he += $h_dejavu->{$projName}->{he};
				$ho += $h_dejavu->{$projName}->{ho};
			}
			$hash->{$chr_h_id}->{$id}->{'dejavu'} = "";
			$nb_patient = $he + $ho;
			if ($nb_project > 0) {
				$hash->{$chr_h_id}->{$id}->{'dejavu'} = $nb_project."/".$nb_patient." ho:$ho,he:$he";
			}
			$hash->{$chr_h_id}->{$id}->{'dejavu'} .= ' (only HO)' if ($chr->project->dejavu_ho());
			$hash->{$chr_h_id}->{$id}->{'chr'} = $chr->ucsc_name();
			$hash->{$chr_h_id}->{$id}->{'position'} = $var->start();
			$hash->{$chr_h_id}->{$id}->{'allele'} = $var->ref_allele().'/'.$var->var_allele();
			$hash->{$chr_h_id}->{$id}->{'nb_all_ref'} = $var->{annex}->{$self->id()}->{nb_all_ref};
			$hash->{$chr_h_id}->{$id}->{'nb_all_mut'} = $var->{annex}->{$self->id()}->{nb_all_mut};
			$hash->{$chr_h_id}->{$id}->{'cadd_score'} = $var->cadd_score();
			if ($var->cosmic()) {
				my @lTmpCosmic = split(':', $var->cosmic());
				$hash->{$chr_h_id}->{$id}->{'cosmic'} = $lTmpCosmic[0];
			}
			$hash->{$chr_h_id}->{$id}->{'clinvar'} = $var->text_clinvar();
			$hash->{$chr_h_id}->{$id}->{'min_pop_freq'} = $var->min_pop_name().': '.$var->min_pop_freq();
			$hash->{$chr_h_id}->{$id}->{'max_pop_freq'} = $var->max_pop_name().': '.$var->max_pop_freq();
			
			my $seq = $chr->sequence(($var->start()-21), ($var->end()-1));
			$seq .= '['.$var->ref_allele().'/'.$var->var_allele().']';
			$seq .= $chr->sequence(($var->start()+1), ($var->end()+21));
			$hash->{$chr_h_id}->{$id}->{'sequence'} = $seq;
			if ($var_patient_ho->contains($id)) {
				$hash->{$chr_h_id}->{$id}->{'he_ho'} = 'ho';
			}
			else { $hash->{$chr_h_id}->{$id}->{'he_ho'} = 'he' }
			foreach my $gene (@{$var->getGenes()}) {
				my $g_id = $gene->id();
				next unless ($gene->getVariantsVector->contains($id));
				$var->annotation();
				$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'external_name'} = $gene->external_name();
				$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'description'} = $gene->description();
				foreach my $t (@{$gene->getTranscripts()}) {
					next if ($var->end() < $t->start() or $t->end() < $var->start());
					next unless (exists $var->annotation->{$t->id()});
					my @ok;
					my $annot_trans;
					eval { $annot_trans = $var->variationType($t) };
					if ($@) { $annot_trans = 'N.A'; }
					foreach my $cons (split(',', $annot_trans)) {
						my @lTmp = split(';', $self->project->buffer->config->{ensembl_annotations}->{$cons});
						push(@ok, $lTmp[1]);
					}
					next if (scalar(@ok) == 0);
					$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'external_name'} = $t->external_name();
					$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'consequence'} = join(',', @ok);
					foreach my $cons (@ok) {
						$hash->{$chr_h_id}->{$id}->{'consequences'}->{$cons} = undef;
					}
					$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'cdna_position'} = $t->translate_position($var->start());
					$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'exon'} = $t->findExonNumber($var->start());
					if ($hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'exon'} == -1) {
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'exon'} = $t->findNearestExon($var->start(), $var->end());
					}
					if ($var->isCoding($t)) {
						my $prot = $t->getProtein();
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'protein'} = $prot->id();
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'protein_xref'} = $t->{'external_protein_name'};
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'nomenclature'} = $var->getNomenclature($t);
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'cds_position'} = $var->getOrfPosition($prot);
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'protein_position'} = $var->getProteinPosition($prot);
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'polyphen_status'} = $var->polyphenStatus($prot);
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'sift_status'} = $var->siftStatus($prot);
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'polyphen_score'} = $var->polyphenScore($prot);
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'sift_score'} = $var->siftScore($prot);
						my $protAA = $var->getProteinAA($prot);
						my $chanAA = $var->changeAA($prot);
						if ($protAA and $chanAA) {
							$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'aa'} = $protAA.'/'.$chanAA;
						}
					}
					else {
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'protein'} = '';
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'protein_xref'} = '';
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'nomenclature'} = '';
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'cds_position'} = '';
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'protein_position'} = '';
						$hash->{$chr_h_id}->{$id}->{'genes'}->{$gene->id()}->{'transcripts'}->{$t->id()}->{'aa'} = '';
					}
				}
			}
		}
		$hTrans = undef;
		$self->xls_var_id($hash);
	}
}
sub setOrigin {
	my ($self,$chr) = @_;
	unless ($chr->size_vector()) {
		$self->{origin}->{all}->{$chr->name} = Bit::Vector->new(0);
		$self->{origin}->{ho}->{$chr->name} = Bit::Vector->new(0);
		$self->{origin}->{he}->{$chr->name} = Bit::Vector->new(0);
	}
	else {
		#my $no_patients = $chr->get_vector_patient($self);
		my $hh = $chr->get_vector_patient($self,['all','ho','he']);
		#$hh =  $no_patients->get($self->id()) unless $hh;;
		confess() unless  $hh;
		$self->{origin}->{all}->{$chr->name} = $hh->{all};
		$self->{origin}->{ho}->{$chr->name} = $hh->{ho};
		$self->{origin}->{he}->{$chr->name} = $hh->{he};
	}
	return;
	
}

sub _getRocksVector {
	my ($self,$chr,$type) = @_;
	return $self->{origin}->{$type}->{$chr->name}->Clone if exists $self->{origin}->{$type}->{$chr->name};
	$self->{origin}->{$type}->{$chr->name} = $chr->rocks_vector("r")->get_vector_patient($self,"$type");
	return $self->{origin}->{$type}->{$chr->name}->Clone;
}

sub getVectorOrigin {
	my ($self,$chr) = @_;
	confess("no chr") unless $chr;
	return $self->{origin}->{all}->{$chr->name}->Clone if exists $self->{origin}->{all}->{$chr->name};
	if ($self->project->isRocks){
		return $self->_getRocksVector($chr,"all");
	}
	$self->setOrigin($chr);
	#warn   $self->{origin}->{all}->{$chr->name};
	return $self->{origin}->{all}->{$chr->name}->Clone;
}

sub getVectorOriginCategory  {
	my ($self,$chr,$type) = @_;
	confess() unless $chr;
	confess() unless $type;
	return $self->{origin}->{$type}->{$chr->name}->Clone if exists $self->{origin}->{$type}->{$chr->name};
	if ($self->project->isRocks){
		my $vector_all = $self->getVectorOrigin($chr)->Clone();
		if ($vector_all->Norm() == 0) {
			$self->{origin}->{$type}->{$chr->name} = $vector_all;
		}
		else {
			$self->{origin}->{$type}->{$chr->name} = $self->_getRocksVector($chr,$type);
		}
		return $self->{origin}->{$type}->{$chr->name};
	}
	$self->setOrigin($chr);
	return $self->{origin}->{$type}->{$chr->name}->Clone;
}

sub getVectorOriginHe {
	my ($self,$chr) = @_;
	return $self->getVectorOriginCategory($chr, 'he');
}


sub getVectorOriginHo {
	my ($self,$chr) = @_;
	return $self->getVectorOriginCategory($chr, 'ho');
}

sub getVectorOriginJunctionsRI {
	my ($self,$chr) = @_;
	return $self->getVectorOriginCategory($chr, 'ri');
}

sub getVectorOriginJunctionsSE {
	my ($self,$chr) = @_;
	return $self->getVectorOriginCategory($chr, 'se');
}

sub getVectorOriginJunctionsN {
	my ($self,$chr) = @_;
	return $self->getVectorOriginCategory($chr, 'n');
}

sub getVectorOriginJunctionsD {
	my ($self,$chr) = @_;
	return $self->getVectorOriginCategory($chr, 'd');
}

sub getVectorOriginJunctionsA {
	my ($self,$chr) = @_;
	return $self->getVectorOriginCategory($chr, 'a');
}

sub getVectorOriginJunctionsDA {
	my ($self,$chr) = @_;
	return $self->getVectorOriginCategory($chr, 'da');
}

sub getVectorOriginJunctionsNDA {
	my ($self,$chr) = @_;
	return $self->getVectorOriginCategory($chr, 'nda');
}

sub getVectorOriginJunctionsRatio {
	my ($self,$chr,$ratio) = @_;
	confess() if not $ratio =~ /^[1-9]0$/; 
	return $self->getVectorOriginCategory($chr, 'junc_ratio_'.$ratio);
}



# recupere les infos de DejaVu pour un variant appartenant a ce patient
sub getDejaVuInfos {
	my ($self, $var_id) = @_;
	confess();
	my $nb_project = 0;
	my $nb_patient = 0;
	my $he = 0;
	my $ho = 0;
	my $list = $self->project->getDejaVuInfos($var_id);
	delete $list->{ $self->project->name() };
	my $st_project;
	foreach my $oproject (keys %$list){
		$nb_project ++;
		$st_project=$oproject.":";
		$st_project .= $list->{$oproject}->{patients};
		$he += $list->{$oproject}->{he};
		$ho += $list->{$oproject}->{ho};
	}
	my $dejavu = "-";
	$nb_patient = $he + $ho;
	if ($nb_project > 0) { $dejavu = "Proj:$nb_project | Pat:$nb_patient | ho:$ho, he:$he"; }
	return $dejavu;
}

# creer un objet GenBoVariation / GenBoInsertion / GenBoDeletion appartenant a ce patient, a partir de son ID.
sub getGenBoVariant {
	my ($self, $var_id) = @_;
	return $self->project->_newVariant_cache($var_id, $self);
}

# methode pour mettre le patient in the attic
sub setInTheAttic {
	my ($self) = @_;
	$self->in_the_attic(1);
}

# renvoie un hash avec le comptage de chaque categories (nb snps, ins, del, intronic, silent, utr, etc...)
sub getHashCountTypeVariants {
	my ($self, $chr) = @_;
	confess("\n\nERROR: GenBoPatientCache->getHashCountTypeVariants() method need a GenBoChromosomeCache and a category object in argument. Die.\n\n") unless($chr);
	my $hash;
	my @categories = ('all_variations', 'substitution', 'insertion', 'deletion', 'large_deletion', 'ho', 'he','large_duplication');
	my $hashTypeVariants = $self->getHashTypeVariants($chr);
	foreach my $type (@categories) {
		if (exists $hashTypeVariants->{$type}) {
			$hash->{$type} = $self->countThisVariants( $hashTypeVariants->{$type} );
		}
		else {
			$hash->{$type} = 0;
		}
	}
	$hash->{'cnv'} = 0;
	if (exists $hash->{'large_deletion'}) {
		$hash->{'cnv'} += $hash->{'large_deletion'};
	}
	if (exists $hash->{'large_duplication'}) {
		$hash->{'cnv'} += $hash->{'large_duplication'};
	}
	$hash->{'he'} = $hash->{'all_variations'} -$hash->{'ho'}; 
	#warn "\n\n";
	#warn "CHR".$chr->id();
	#warn Dumper $hash;;
	#die if ($chr->id() eq '1');
	return $hash;
}

# renvoie un hash avec le vector de chaque annotation
sub getHashTypeVariants {
	my ($self, $chr) = @_;
	confess("\n\nERROR: GenBoPatientCache->getHashTypeVariants() method need a GenBoChromosomeCache and a category object in argument. Die.\n\n") unless($chr);
	my $hash;
	my @categories = ('substitution', 'insertion', 'deletion', 'large_deletion', 'ho', 'he','large_duplication');
	foreach my $category (@categories) {
		my $var_tmp = $self->getCategoryVariantsVector($chr, $category);
		if (exists $hash->{$category}) {
			$hash->{$category} += $var_tmp;
		}
		else { $hash->{$category} = $var_tmp; }
		if (($category eq 'substitution') or ($category eq 'deletion') or ($category eq 'large_deletion') or ($category eq 'insertion')or ($category eq 'large_duplication')) {
			if (exists $hash->{'all_variations'}) {
				$hash->{'all_variations'} += $var_tmp;
			}
			else { $hash->{'all_variations'} = $var_tmp; }
		}
	}
	return $hash;
}

# renvoie le vector des HO
sub getHo {
	my ($self, $chr_obj) = @_;
	confess("\n\nERROR: GenBoPatientCache->getHo() method need a GenBoChromosomeCache object in argument. Die.\n\n") unless($chr_obj);
	return $self->getVectorOriginHo($chr_obj);
	return $self->getCategoryVariantsVector( $chr_obj, 'ho' );
}

# renvoie le vector des HE
sub getHe {
	my ($self, $chr_obj) = @_;
	confess("\n\nERROR: GenBoPatientCache->getHe() method need a GenBoChromosomeCache object in argument. Die.\n\n") unless($chr_obj);
	return $self->getVectorOriginHe($chr_obj)->Clone;
	return $self->getCategoryVariantsVector( $chr_obj, 'he' );
}
sub getVectorHe {
	my ($self, $chr_obj) = @_;
	return $self->getHe($chr_obj);
}
sub getVectorHo {
	my ($self, $chr_obj) = @_;
	return $self->getHo($chr_obj);
}
# renvoie le vector avec TOUS les variants
sub getAll {
	my ($self, $chr_obj) = @_;
	confess("\n\nERROR: GenBoPatientCache->getAll() method need a GenBoChromosomeCache object in argument. Die.\n\n") unless($chr_obj);
	return $self->getVariantsVector( $chr_obj );
}

sub getVariantsFromCallingMethod {
	my ($self, $chr_obj, $calling_method) = @_;
	confess("\n\nERROR: GenBoPatientCache->getVariantsFromCallingMethod() method need a GenBoChromosomeCache object in argument. Die.\n\n") unless($chr_obj);
	confess("\n\nERROR: GenBoPatientCache->getVariantsFromCallingMethod() method need a calling_method name object in argument. Die.\n\n") unless($calling_method);
	if (exists $self->calling_methods_categories->{$chr_obj->id()}->{$calling_method}) { return $self->calling_methods_categories->{$chr_obj->id()}->{$calling_method}; }
	return $chr_obj->getNewVector();
} 

sub delete_variants {
	my ($self, $chr, $variants) = @_;
	confess();
	confess("\n\nchr missing...\n\n") unless ($chr);
	confess("\n\nvector missing...\n\n") unless ($variants);
	$chr->patients_categories->{$self->name()} -= $variants;
}

sub delete_variants_from_gene {
	my ($self, $chr, $gene) = @_;
	confess("\n\nchr missing...\n\n") unless ($chr);
	confess("\n\ngene missing...\n\n") unless ($gene);
	$chr->patients_categories->{$self->name()} -= $gene->getVariantsVector();
}

sub getIntspan {
	my $self = shift;
	return Set::IntSpan::Fast::XS->new( $self->getVariantsVector( $self->chromosome() )->to_Enum() );
}

# renvoies le nombre de genes de ce patient
sub countGenes {
	my ($self, $chr) = @_;
	return 0 if ($self->excluded() and $self->excluded() eq '1');
	my $nb = 0;
	my $vector_tmp = $chr->getNewVector();
	foreach my $gene (@{$chr->getGenes()}) {
		next if ($gene->is_intergenic());
		if ($gene->patients_found()) {
			$nb++ if (exists $gene->patients_found->{$self->name()});
		}
		else {
			next if (exists $gene->patients_never_found->{$self->name()});
			$vector_tmp->Empty();
			$vector_tmp->Intersection($self->getVariantsVector($chr), $gene->getVariantsVector());
			unless ($vector_tmp->is_empty()) {
				$nb++;
			}
		}
	}
	return $nb;
}

sub getVariantsVector {
	my ($self, $chr_obj) = @_;
	return $self->getVectorOrigin($chr_obj);
	confess();
	confess("\n\nERROR: GenBoPatientCache->getVariantsVector() method need a GenBoChromosomeCache object in argument. Die.\n\n") unless($chr_obj);
#	unless (exists $chr_obj->patients_categories->{$self->name()}) {
#		confess() #deleteX 
#		$chr_obj->patients_categories->{$self->name()} = $chr_obj->getNewVector();
#		$chr_obj->patients_categories->{$self->name().'_he'} = $chr_obj->getNewVector();
#		$chr_obj->patients_categories->{$self->name().'_ho'} = $chr_obj->getNewVector();
#		return $chr_obj->patients_categories->{$self->name()};
#	}
#	if ($chr_obj->variants_intersected()) {	
#			ceonfess() #deleteX 
#		$chr_obj->patients_categories->{$self->name()}->Intersection( $chr_obj->patients_categories->{$self->name()}, $chr_obj->variants_intersected() );
#		$chr_obj->patients_categories->{$self->name().'_he'}->Intersection( $chr_obj->patients_categories->{$self->name().'_he'}, $chr_obj->variants_intersected() );
#		$chr_obj->patients_categories->{$self->name().'_ho'}->Intersection( $chr_obj->patients_categories->{$self->name().'_ho'}, $chr_obj->variants_intersected() );
#	}
#	if ($chr_obj->variants_excluded() and $chr_obj->variants_excluded() eq 'all') {
#		confess() #deleteX 
#		$chr_obj->patients_categories->{$self->name()} -= $chr_obj->variants_excluded();
#		$chr_obj->patients_categories->{$self->name().'_he'} -= $chr_obj->variants_excluded();
#		$chr_obj->patients_categories->{$self->name().'_ho'} -= $chr_obj->variants_excluded();
#	}
	#warn $chr_obj->getVariantsVector();

	my $v1 = $self->getVectorOrigin();
	$v1 &= $chr_obj->getVariantsVector();
	return $v1;
}

sub getCategoryVariantsVector {
	my ($self, $chr_obj, $cat) = @_;
	confess("\n\nERROR: GenBoPatientCache->getCategoryVariantsVector() method need a GenBoChromosomeCache and a category object in argument. Die.\n\n") unless($chr_obj);
	confess("\n\nERROR: GenBoPatientCache->getCategoryVariantsVector() method need a GenBoChromosomeCache and a category object in argument. Die.\n\n") unless($cat);
	my $var = $self->getVariantsVector($chr_obj)->Clone();
	$var->Intersection($var, $chr_obj->getVariantsVector());
	$var += $self->getVariantsVector( $chr_obj );
	if ($cat eq 'he' or $cat eq 'ho') {
		unless (exists $chr_obj->patients_categories->{$self->name().'_'.$cat}) {
			warn $var->Size();
			warn $chr_obj->getVariantsVector()->Size();
			$var->Empty();
			return;
		}
		$chr_obj->patients_categories->{$self->name().'_'.$cat}->Intersection($chr_obj->patients_categories->{$self->name().'_'.$cat}, $chr_obj->getVariantsVector());
		$var->Intersection( $var, $chr_obj->patients_categories->{$self->name().'_'.$cat} );
	}
	else { $var->Intersection( $var, $chr_obj->getCategoryVariantsVector($cat) ); }
	return $var;
}

sub getStatsCategoryVariantsVector {
	my ($self, $chr_obj, $cat) = @_;
	confess("\n\nERROR: GenBoPatientCache->getStatsCategoryVariantsVector() method need a GenBoChromosomeCache and a category object in argument. Die.\n\n") unless($chr_obj);
	confess("\n\nERROR: GenBoPatientCache->getStatsCategoryVariantsVector() method need a GenBoChromosomeCache and a category object in argument. Die.\n\n") unless($cat);
	if ($self->in_the_attic()) {
		unless (exists $self->{stats_categories}->{$chr_obj->id()}->{$cat}) {
			$self->{stats_categories}->{$chr_obj->id()}->{$cat} = $chr_obj->getNewVector();
			$self->{stats_categories}->{$chr_obj->id()}->{$cat} += $chr_obj->patients_categories->{$self->name().'_he'};
			$self->{stats_categories}->{$chr_obj->id()}->{$cat} += $chr_obj->patients_categories->{$self->name().'_ho'};
			$self->{stats_categories}->{$chr_obj->id()}->{$cat}->Intersection( $self->{stats_categories}->{$chr_obj->id()}->{$cat}, $chr_obj->getCategoryVariantsVector($cat) );
		}
	}
	else {
		unless (exists $self->{stats_categories}->{$chr_obj->id()}->{$cat}) {
			$self->{stats_categories}->{$chr_obj->id()}->{$cat} = $self->getCategoryVariantsVector($chr_obj, $cat);
		}
	}
	return $self->{stats_categories}->{$chr_obj->id()}->{$cat};
}

sub getVariantsCategory {
	my $self = shift;
	confess();
}

sub getVectorVariants {
	my ($self,$chr) = @_;
	return $self->{vector}->{variants}->{$chr->id} if exists  $self->{vector}->{variants}->{$chr->id};
	my $hash = $chr->get_vector_patient($self,['all']);
	$self->{vector}->{variants}->{$chr->id} = $hash->{all};
	return $self->{vector}->{variants}->{$chr->id};
}

sub getVectorVariations {
	my($self,$chr) = @_;
	confess() unless $chr;
	return  $self->getVectorVariants($chr);
}

sub getVectorPublicVariations {
	my($self,$chr) = @_;
	my $vector = $self->getVectorVariations($chr)->Clone();
	$vector &= $chr->getVectorPublicVariations();
}

sub getVectorSubstitutions {
	my($self,$chr) = @_;
	#return  $patient->{global_categories}->{$self->id()}->{substitution}
	return  $self->{global_categories}->{substitution}->{$chr->id} if exists $self->{vector}->{substitution}->{$chr->id};
	$self->{vector}->{substitution}->{$chr->id} =  $self->getVectorVariants($chr)->Clone;
	$self->{vector}->{substitution}->{$chr->id} &= $chr->getVectorSubstitutions();
	return 	$self->{vector}->{substitution}->{$chr->id};
}

sub getVectorInsertions {
	my($self,$chr) = @_;
	return  $self->{vector}->{insertion}->{$chr->id} if exists $self->{vector}->{insertion}->{$chr->id};
	$self->{vector}->{insertion}->{$chr->id} =  $self->getVectorVariants($chr)->Clone;
	$self->{vector}->{insertion}->{$chr->id}  &= $chr->getVectorInsertions();
	return 	$self->{vector}->{insertion}->{$chr->id};
}

sub getVectorDeletions {
	my($self,$chr) = @_;
	return  $self->{vector}->{deletion}->{$chr->id} if exists $self->{vector}->{deletion}->{$chr->id};
	$self->{vector}->{deletion}->{$chr->id} =  $self->getVectorVariants($chr)->Clone;
	$self->{vector}->{deletion}->{$chr->id} &= $chr->getVectorDeletions();
	return 	$self->{vector}->{deletion}->{$chr->id};
	confess();
}

sub getVectorIndels {
	my($self,$chr) = @_;
	return  $self->{vector}->{indel}->{$chr->id} if exists $self->{vector}->{indel}->{$chr->id};
	$self->{vector}->{indel}->{$chr->id} =  $self->getVectorDeletions($chr)->Clone;
	$self->{vector}->{indel}->{$chr->id} +=  $self->getVectorInsertions($chr);
	return 	$self->{vector}->{indel}->{$chr->id};
}

sub countPublicVariations {
	my($self,$chr) = @_;
	if ($chr){
		my $v = $self->getVectorPublicVariations($chr);
		return $self->countThisVariants($v);
	}
	else {
		my $nb = 0;
		foreach my $chr (@{$self->project->getChromosomes}){
			my $v = $self->getVectorPublicVariations($chr);
			$nb+= $self->countThisVariants($v);
		}
		return $nb;
	}
}

sub countSubstitutions {
	my($self,$chr) = @_;
	if ($chr){
		my $v = $self->getVectorSubstitutions($chr);
		return $self->countThisVariants($v);
	}
	else {
		my $nb = 0;
		foreach my $chr (@{$self->project->getChromosomes}){
			my $v = $self->getVectorSubstitutions($chr);
			$nb+= $self->countThisVariants($v);
		}
		return $nb;
	}
}

sub countVariations{
	my($self,$chr) = @_;
	if ($chr){
		my $v = $self->getVectorVariants($chr);
		return $self->countThisVariants($v);
	}
	else {
		my $nb = 0;
		foreach my $chr (@{$self->project->getChromosomes}){
			my $v = $self->getVectorVariants($chr);
			$nb+= $self->countThisVariants($v);
		}
		return $nb;
	}
}

sub countIndels {
	my($self,$chr) = @_;
	if ($chr){
		my $v = $self->getVectorIndels($chr);
		return $self->countThisVariants($v);
	}
	else {
		my $nb = 0;
		foreach my $chr (@{$self->project->getChromosomes}){
			my $v = $self->getVectorIndels($chr);
			$nb+= $self->countThisVariants($v);
		}
		return $nb;
	}
}
 
sub countHomozygote {
	my($self,$chr) = @_;
	if ($chr){
		my $v = $self->getVectorOriginHo($chr);
		return $v->Norm()
	}
	else {
		my $nb = 0;
		foreach my $chr (@{$self->project->getChromosomes}){
			my $v = $self->getVectorOriginHo($chr);
			$nb+= $v->Norm();
		}
		return $nb;
	}
}

sub countHeterozygote {
	my($self,$chr) = @_;
	if ($chr){
		my $v = $self->getVectorOriginHe($chr);
		return $v->Norm();
	}
	else {
		my $nb = 0;
		foreach my $chr (@{$self->project->getChromosomes}){
			my $v = $self->getVectorOriginHe($chr);
			$nb+= $v->Norm();
		}
		return $nb;
	}
}

sub getVectorLargeDeletions {
	my ($self,$chr) = shift;
	die() unless $chr;
	my $v = $chr->getVectorLargeDeletions->Clone();
	$v &= $self->getOriginVector();
}
sub getVectorLargeDuplications {
	my ($self,$chr) = shift;
	die() unless $chr;
	my $v = $chr->getVectorLargeDuplications->Clone();
	$v &= $self->getOriginVector();
	return $v;
}


sub setVariants {
	my ($self, $type) = @_;
	foreach my $chr (@{$self->project->getChromosomes()}) {
		my $vector = $chr->getNewVector();
		if ($type eq 'variations') {
			$vector->Intersection( $self->getVariantsVector($chr), $chr->getVectorSubstitutions() );
		}
		elsif ($type eq 'insertions') {
			$vector->Intersection( $self->getVariantsVector($chr), $chr->getVectorInsertions() );
		}
		elsif ($type eq 'deletions') {
			$vector->Intersection( $self->getVariantsVector($chr), $chr->getVectorDeletions() );
		}
		elsif ($type eq 'large_deletions') {
			$vector->Intersection( $self->getVariantsVector($chr), $chr->getVectorLargeDeletions() );
		}
		elsif ($type eq 'large_duplications') {
			$vector->Intersection( $self->getVariantsVector($chr), $chr->getVectorLargeDuplications() );
		}
		foreach my $var (@{$chr->getListVarObjects($vector)}) {
			$self->{$var->type_object()}->{$var->id()} = undef;
			unless (exists $self->project->{objects}->{$type}->{$var->id()}) {
				$self->project->{objects}->{$type}->{$var->id()} = $var;
			}
		}
	}
}

# call by getGenes return all genes for this patient it will be smarter if it's not set but getGenes just if the vector change 
# we have to discuss that
sub setGenes {
	my ($self) = @_;
	my $hGenesId;
	foreach my $chr (@{$self->project->getChromosomes()}) {
	my $vector = $self->getVariantsVector($chr);
	my @enum =  split(",",$vector->to_Enum);
		foreach my $en (@enum){
			my ($start,$end) = split("-",$en);
			$end = $start +1 unless $end;
			my $results = $chr->genes_tree->fetch($start,$end+1);
			foreach my $r (@$results){
				$hGenesId->{$r} = undef;
			}
		}
	}
	return $hGenesId;
}

sub getRegionHo_intspan {
	my ($self, $chr, $nbVar, $filter_regionho_sub_only) = @_;
	confess() unless ($nbVar);
	my $vecReg = Bit::Vector->new($chr->size_vector());
	my $chrname = $chr->name();
	my $spanHo = Set::IntSpan::Fast::XS->new();	# creation du span pour stocker la liste des regions Ho du patient sur le chr
	return $vecReg if ($spanHo eq "Y");
	    		
	# creation de nouveaux vecteurs
	my $vecSubHe = $self->getHe($chr);
	my $vecSubHo = $self->getHo($chr);

	# si option ok - seulement les SUB pour le filtre
	if ($filter_regionho_sub_only) {	
		$vecSubHe->Intersection($vecSubHe, $chr->getVectorSubstitutions());
		$vecSubHo->Intersection($vecSubHo, $chr->getVectorSubstitutions());
	}

	# recuperation des regions Ho  etendues aux zones sans variants	
	$vecReg->Complement($vecSubHe);
	
    # pour ne conserver que les regions dont le nombre de variants Ho est > a filtre 
    # et redelimiter la region Ho (au premier / dernier variant Ho rencontre dans la region)
    my @H_index = split(/,/,$vecReg->to_Enum());         
			
	foreach my $reg (@H_index) {
		my @index=split(/-/,$reg);
		if (scalar(@index) == 2) {
		    my $start = $index[0];
		    my $end =   $index[-1];
		    my $len = ($end - $start + 1);
		    next if ($len < $nbVar);
			my $sub_vector = $chr->getNewVector($len);
			$sub_vector->Interval_Copy($vecReg, 0, $start, $len);
			my $sub_ho = $chr->getNewVector($len);
			$sub_ho->Interval_Copy($vecSubHo, 0, $start, $len);
			$sub_vector->Intersection($sub_vector, $sub_ho);
			if ($self->countThisVariants($sub_vector) >= $nbVar) {
				$spanHo->add_range($start,$end);
			}
		}
	}
	return $spanHo;
}

has hash_regionHo_vector => (
	is 		=> 'rw',
	lazy	=> 1,
	default => sub { {} },
);

has hash_stats_regionHo_vector => (
	is 		=> 'rw',
	lazy	=> 1,
	default => sub { {} },
);

sub getRegionHo {
	my ($self, $chr, $nbVar, $filter_regionho_sub_only, $compute) = @_;
	confess("\n\nERROR: nedd nbVar ! Die;\n\n") unless ($nbVar);
	return $self->hash_regionHo_vector->{$chr->id()}->{$nbVar} if (exists $self->hash_regionHo_vector->{$chr->id()}->{$nbVar});
	my $key = "region_ho_25_".$self->name;
	if ($nbVar == 25 && $self->project->isRocks && (! defined $compute)){
		$self->{hash_regionHo_vector}->{$chr->id()}->{$nbVar} = $chr->rocks_vector->get_vector_transmission($self,"region_ho_25");
		return $self->{hash_regionHo_vector}->{$chr->id()}->{$nbVar};
	}
	my $intspan_regionHo = $self->getRegionHo_intspan($chr, $nbVar, $filter_regionho_sub_only);
	if ($intspan_regionHo->is_empty()) {
		$self->{hash_regionHo_vector}->{$chr->id()}->{$nbVar} = $chr->getNewVector();
		return $chr->getNewVector();
	}
	my $v_ho = Bit::Vector->new_Enum($chr->size_vector(), join(',', $intspan_regionHo->as_array()));
	$self->{hash_regionHo_vector}->{$chr->id()}->{$nbVar} = $v_ho;
	return $v_ho;
}

sub getHashStatsRegionsHo{
	my ($self, $chr) = @_;
	confess();
	return $self->hash_stats_regionHo_vector->{$chr->id()} if (exists $self->hash_stats_regionHo_vector->{$chr->id()});
	my $v_ho = $self->getRegionHo($chr);
	my $hStatsRegionsHo;
	my $name = $self->name();
	my @lRegionHo_vec = split(',',$v_ho->to_Enum());
	foreach my $coord_vec (@lRegionHo_vec){
		my @lPositions_vec = split('-', $coord_vec);
		my $vec_start = $lPositions_vec[0];
		my $var_start = $chr->getVarObject($vec_start);
		my $vec_end = $lPositions_vec[1];
		my $var_end = $chr->getVarObject($vec_end);
		my $region = $vec_start.'-'.$vec_end;
		$hStatsRegionsHo->{$region}->{'coord'} = 'chr'.$chr->id().':'.$var_start->start().'-'.$var_end->end();
		$hStatsRegionsHo->{$region}->{'length'} = $var_end->end() - $var_start->start();
		my $thisRegionVecHo = Bit::Vector->new_Enum($chr->size_vector(), $coord_vec);
		my $varInThisRegionHo = $chr->getNewVector();
		$varInThisRegionHo->Intersection($chr->getPatient($name)->getHo($chr), $thisRegionVecHo);
		$hStatsRegionsHo->{$region}->{'nb_var'} = $chr->countThisVariants($varInThisRegionHo);
	}
	$self->{hash_stats_regionHo_vector}->{$chr->id()} = $hStatsRegionsHo;
	return $hStatsRegionsHo;
}
sub getVectorClinvar{
	my ($self,$chr) = @_;
	return Bit::Vector->new(0) unless ($chr->size_vector());
	my $vector_hgmd =  $self->getVectorOrigin($chr); 
	my $v2 =  $chr->vectorClinvarPathogenic;
	return  $vector_hgmd->Shadow() unless $v2;
	$vector_hgmd &= $v2;
	return $vector_hgmd;
}


sub getVectorDM{
	my ($self,$chr) = @_;
	return Bit::Vector->new(0) unless ($chr->size_vector());
	die() unless $chr;
	my $vector_hgmd =  $self->getVectorOrigin($chr); 
	#warn $self->getChromosome->lmdb_score_impact->get("dm")->to_Enum;
	my $v2 =  $chr->vectorDM;
	return  $vector_hgmd->Shadow() unless $v2;
	$vector_hgmd &= $v2;
	
	return $vector_hgmd;
	
}

sub getVectorRatio {
	my ($self, $chr, $filter_name) = @_;
	return Bit::Vector->new(0) unless ($chr->size_vector());
	my $vector_ok = $self->getVariantsVector($chr)->Clone();
	$vector_ok->Intersection($vector_ok, $chr->getVectorScore($self->name().'_'.$filter_name));
	return $vector_ok;
}

sub getJunctionsVector {
	my ($self, $chr) = @_;
	my $vector = $self->getVectorOrigin($chr)->Clone();
	return $vector if ($vector->is_empty);
	$vector->Intersection($vector, $chr->getJunctionsVector());
	return $vector;
}

sub getVectorJunctionsRI {
	my ($self, $chr) = @_;
	my $vector_ri = $self->getVectorOriginJunctionsRI($chr)->Clone();
	$vector_ri->Intersection($vector_ri, $chr->getJunctionsVector());
	return $vector_ri;
}

sub getVectorJunctionsSE {
	my ($self, $chr) = @_;

	my $vector_se = $self->getVectorOriginJunctionsSE($chr)->Clone();
	$vector_se->Intersection($vector_se, $chr->getJunctionsVector());
	return $vector_se;
}

sub getVectorJunctionsN {
	my ($self, $chr) = @_;
	my $vector = $self->getVectorOriginJunctionsN($chr)->Clone();
	$vector->Intersection($vector, $chr->getJunctionsVector());
	return $vector;
}


sub getVectorJunctionsD {
	my ($self, $chr) = @_;
	my $vector = $self->getVectorOriginJunctionsD($chr)->Clone();
	$vector->Intersection($vector, $chr->getJunctionsVector());
	return $vector;
}

sub getVectorJunctionsA {
	my ($self, $chr) = @_;
	my $vector = $self->getVectorOriginJunctionsA($chr)->Clone();
	$vector->Intersection($vector, $chr->getJunctionsVector());
	return $vector;
}

sub getVectorJunctionsDA {
	my ($self, $chr) = @_;
	my $vector = $self->getVectorOriginJunctionsDA($chr)->Clone();
	$vector->Intersection($vector, $chr->getJunctionsVector());
	return $vector;
}

sub getVectorJunctionsNDA {
	my ($self, $chr) = @_;
	my $vector = $self->getVectorOriginJunctionsNDA($chr)->Clone();
	$vector->Intersection($vector, $chr->getJunctionsVector());
	return $vector;
}

sub getVectorJunctionsRatio {
	my ($self, $chr, $ratio) = @_;
	my $vector = $self->getVectorOriginJunctionsRatio($chr, $ratio)->Clone();
	$vector->Intersection($vector, $chr->getJunctionsVector());
	return $vector;
}


sub setJunctions {
	my $self = shift;
	my $h;
	foreach my $chr (@{$self->getProject->getChromosomes()}) {
		next if $chr->not_used();
		my $vector = $self->getJunctionsVector($chr);
		foreach my $junction (@{$chr->getListVarObjects($vector)}) {
			$h->{$junction->id()} = undef;
			unless (exists $self->project->{objects}->{junctions}->{$junction->id()}) {
				$self->project->{objects}->{junctions}->{$junction->id()} = $junction;
			}
		}
	}
	return $h;
}

sub hasJunctions {
	my $self = shift;
	my $has_junctions;
	foreach my $chr (@{$self->getProject->getChromosomes()}) {
		my $v = $self->getJunctionsVector($chr);;
		next if $v->is_empty();
		$has_junctions = 1;
		last;
	}
	return $has_junctions;
}

1;