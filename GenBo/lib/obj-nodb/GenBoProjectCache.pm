package GenBoProjectCache;
use strict;
use Storable qw(retrieve);
use Moose;
use Bit::Vector;
use Bit::Vector::Overload;
use Data::Dumper;
use Config::Std;
use JSON;
use GenBoChromosomeCache;
use GenBoFamilyCache;
use GenBoVariationCache;
use GenBoDeletionCache;
use GenBoInsertionCache;
use GenBoLargeDeletionCache;
use GenBoLargeDuplicationCache;
use GenBoSomaticGroupCache;
use GenBoPanelCache;
use GenBoBundleCache;
use Storable qw(store retrieve freeze dclone thaw);
extends 'GenBoProject', 'GenBoCache';






# dossier avec les fichiers necessaires aux Tests_F
has dir_files_test => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $dir = $self->buffer->config->{'public_data'}->{'HG19'}.$self->buffer->config->{'public_data'}->{'test_f'};
		return $dir;
	},
);

# tag true - pour redireger le cache_vector et le test dejavu en figer
has test => (
	is		=> 'rw',
	lazy	=> 1,
	default => undef,
);

# type de models (ind, fam, somatic)
has typeFilters => (
	is		=> 'rw',
	lazy	=> 1,
	default => 'individual',
);

# path du fichier freeze global_infos (ex: avec les infos MD5 de chaque VCF)
has global_infos => (
	is => 'ro',
	lazy => 1,
	default => sub {
		my $self = shift;
		my $file = $self->getCacheBitVectorDir().'/global_infos.freeze';
		my ($hInfos, $hRes);
		confess() unless -e $file;
		return retrieve $file;
	}
);

# mode par variation ou par gene en IND
has level_ind => (
	is		=> 'rw',
	lazy	=> 1,
	default => 'variation',
);

# mode par variation ou par gene en FAM
has level_fam => (
	is		=> 'rw',
	lazy	=> 1,
	default => 'variation',
);

# valeur du quick search
has filter_text => (
	is		=> 'rw',
	lazy	=> 1,
	default => undef,
);

# hash avec la liste des genes qui passent le Quick Search
has gene_search_ok => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $text = $self->filter_text();
		$text = uc($text);
		$text =~ s/\+/ /;
		my @lText = split(' ', $text);
		my @lOkText;
		my $i = 0;
		while ($i < scalar(@lText)) {
			my $this_text = $lText[$i];
			if ($this_text eq 'NOT') { $i++; }
			else { push(@lOkText, $lText[$i]); }
			$i++;
		}
		my $hash_ok = $self->filter_gene_ontology( \@lOkText );
		foreach my $ok (@lOkText) {
			$hash_ok->{$ok} = undef unless (exists $hash_ok->{$ok});
		}
		return $hash_ok;
	},
);

# hash avec la liste des genes qui NE passent PAS le Quick Search
has gene_search_not => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $text = $self->filter_text();
		$text = uc($text);
		my @lText = split(' ', $text);
		my @lNotText;
		my $i = 0;
		while ($i < scalar(@lText)) {
			my $this_text = $lText[$i];
			if ($this_text eq 'NOT') {
				$i++;
				push(@lNotText, $lText[$i]) if ($i < scalar(@lText));
			}
			$i++;
		}
		return unless (@lNotText);
		my $hash_not = $self->filter_gene_ontology( \@lNotText );
		foreach my $not (@lNotText) {
			$hash_not->{$not} = undef unless (exists $hash_not->{$not});
		}
		return $hash_not;
	},
);

# filtre sur les annotations
has hash_ensembl_annotations => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->buffer->config->{'ensembl_annotations'};
	},
);

# filtre sur le variant (je supprime les variants compris dans leurs vecteurs - les autres filtres sont justes ignores)
has hash_var_type_filters => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->buffer->config->{'variants_type_filters'};
	},
);

# comme ceux du dessus mais sur les frequences
has hash_frequence_filters => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->buffer->config->{'frequence_filters'};
	},
);

# comme ceux du dessus mais sur les confidence
has hash_confidence_filters => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->buffer->config->{'confidence_filters'};
	},
);

# comme ceux du dessus mais sur les predictions
has hash_prediction_filters => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->buffer->config->{'predictions_type_filters'};
	},
);

# comme ceux du dessus mais sur les cadd score
has hash_cadd_filters => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->buffer->config->{'cadd_filters'};
	},
);

# annotations que l on garde meme apres filtre sur Polyphen et Sift (sauf si l utilisateur les decoche de l interface)
has hash_prediction_except_filters => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->buffer->config->{'predictions_except_filters'};
	},
);

# si autre que 'no' (ex: variants or genes), alors l'interface renvoie un fichier XLS
has get_xls  => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> 'no',
);

# sous partie des stats (ici par region)
has stats_region => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub { [] },
);

# sous partie des stats (ici par genes)
has stats_genes => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub { [] },
);

# retourne un hash avec les relations gene / maladies reference dans gene atlas
has gene_atlas => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $file = $self->buffer->{config}->{gene_atlas_diseases}->{dir};
		$file .= $self->buffer->{config}->{gene_atlas_diseases}->{genes};
		my $genes_atlas;
		if (-e $file) {
	 		return $genes_atlas = retrieve($file);
		}
		return;
	},
);

# hash avec les chr_id que l'on veut regarder avec cette requete (desormais chr par chr le lancement de l'interface)
has filter_chromosome => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> undef,
);

# hash avec les gene id que l on a demander a intersecter
has filter_genes_intersect => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub { {} },
);

# hash avec les patient id que l on a demander a mettre en in the attic uniquement pour faire le gene intersect et qui n auront pas le logo intheattic en resultat
has filter_attic_genes => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub { {} },
);

# utile dans la methode add_children pour la construction de l arbre de gene atlas
has uuid => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> 0,
);

# hash avec l'ensemble des filtres par regions sur chaque chromosome
has filter_region => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub { {} },
);

# renvoie HASH avec infos du projet
#has infosProject => (
#	is		=> 'ro',
#	lazy	=> 1,
#	default => sub {
#		my $self = shift;
#		#return retrieve $self->dir_files_test().'/infosProject.freeze' if ($self->test());
#		warn "jesuis la";
#		warn $self->buffer->getQuery;
#		die();
#		return $self->buffer->getQuery->getProjectByName( $self->name() );
#	},
#);

# genes a regarder uniquement si utilise une capture de diag
has only_genes => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> undef,
);

# si option onl_genes utilise, renvoie hash avec les genes 
has only_genes_found => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return undef unless ($self->only_genes());
		my $h;
		foreach my $gene_name (keys %{$self->only_genes()}) {
			if ($self->if_exists_liteObject($gene_name, 'gene')) {
				$h->{uc($gene_name)} = undef;
			}
		}
		return $h;
	},
);

# hash qui caracterisent chaque type d'objet (lors de sa creation)
sub hashTypeObject {
	my $self = shift;
	my $hashTypeObject = {
		'phenotypes'	=> 'GenBoPhenotype',
		'panels'		=> 'GenBoPanelCache',
		'bundles'		=> 'GenBoBundle',
		'captures'		=> 'GenBoCapture',
	 	'variations'	=> 'GenBoVariationCache',
	 	'deletions'		=> 'GenBoDeletionCache',	
	 	'insertions'	=> 'GenBoInsertionCache',
	 	'large_deletions'	=> 'GenBoLargeDeletionCache',	
	 	'large_duplications'=> 'GenBoLargeDuplicationCache',
	 	'references'	=> 'GenBoReference',
	 	'transcripts'	=> 'GenBoTranscriptCache',
	 	'proteins'		=> 'GenBoProtein',
	 	'positions'		=> 'Position',
	 	'exons'			=> 'GenBoExon',
	 	'primers'		=> 'GenBoPrimer',
	 	'genes'			=> 'GenBoGeneCache',
		'chromosomes'	=> 'GenBoChromosomeCache',
		'introns'		=> 'GenBoIntron',
		'mnps'			=> 'GenBoMnp',
		'multiplexes'	=> 'GenBoPcrMultiplex',
	 	'part_chromosomes' => 'GenBoPartChromosome',	
		'families'		   => 'GenBoFamilyCache',
		'patients'		   => 'GenBoPatientCache',
		'somatic_groups'   => 'GenBoSomaticGroupCache',
	};
	return $hashTypeObject;
}

# tag true if using interface filters -> may be some genes was deleted. So getGenes can be in errors if using as usually
has used_interface_filters => (
	is => 'rw',
	lazy => 1,
	default => undef,	
);

has hash_patients_name => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		my $h;
		my $id = 1;
		foreach my $patient (@{$self->getPatients()}) {
			$h->{$patient->name()}->{id} = $id;
			$h->{$patient->name()}->{fam} = $patient->family();
			$h->{by_id}->{$id} = $patient->name();
			$id++;
		}
		return $h;
	},
);

# nom du model utilise
has model => (
	is      => 'rw',
	lazy    => 1,
	default => undef,
);

has is_research_by_genes => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return 1 if ($self->level_ind() eq 'gene');
		return 1 if ($self->level_fam() eq 'gene');
		return;
	},
);

has hashTypeStats => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $hash = {
			'substitution' => undef,
			'insertion' => undef,
			'deletion' => undef,
			'large_deletion' => undef,
			'silent' => undef,
			'utr' => undef,
			'splicing' => undef,
			'coding' => undef,
			'stop' => undef,
#			'phase' => undef,
			'frameshift' => undef,
			'he' => undef,
			'ho' => undef,
		};
		return $hash;
	},
);

has dejavu_ho => (
	is      => 'rw',
	lazy    => 1,
	default => undef,
);

# filtre mosaique min pvalue
has mosaique_min_p => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		$self->buffer->config->{mosaic}->{min_pvalue};
	},
);



##### METHODS #####




sub isXlsOutput {
	my ($self, $type) = @_;
	$self->get_xls($type);
	$self->{hashTypeObject}->{'transcripts'}     = 'GenBoTranscript';
	$self->{hashTypeObject}->{'variations'}      = 'GenBoVariation';
	$self->{hashTypeObject}->{'insertions'}      = 'GenBoInsertion';
	$self->{hashTypeObject}->{'deletions'}       = 'GenBoDeletion';
	$self->{hashTypeObject}->{'large_deletions'} = 'GenBoDeletion';
}


# dossier TMP, notamment util pour exporter en XLS
sub getTmpDir {
	my $self = shift;
	my $dir = $self->buffer()->getDataDirectory("cache").'/tmp/';
	unless (-d $dir) {
		`mkdir $dir`;
		`chmod 777 $dir`;
	}
	return $dir;
}

# methode pour recupere un gene suivant une ontology
sub filter_gene_ontology {
	my ($self, $listText) = @_;
	my $apph = $self->buffer->get_go_db();
	return unless $apph;
	my $hGenes;
	foreach my $search (@$listText) { 
		my $terms = $apph->get_terms({search=>"$search*",search_fields=>"name,synonym,definition"});
		foreach my $term (@$terms){
			my $products = $apph->get_deep_products({term=>$term->name});
			foreach my $g (@$products){
				push(@{$hGenes->{$g->symbol}},$term->name);
			}
		}
	}
	return $hGenes;
}

# parse l'option filter_chromosome de l'interface et rempli le hash filter_chromosome de GenBoProjectCache
sub set_filter_chromosome {
	my ($self, $filters, $separator) = @_;
	unless ($separator) { $separator = ','; }
	my $hash;
	foreach my $chrId (split( $separator, $filters )) {
		$hash->{$chrId} = undef;
	}
	$self->filter_chromosome( $hash );
}

# parse l'option filter_region de l'interface et rempli le hash filter_region et filter_chromosome de GenBoProjectCache
sub set_filter_region {
	my ($self, $filters, $separator) = @_;
	unless ($separator) { $separator = ','; }
	my $hashChromosomes = $self->filter_chromosome();
	my $hashFilters;
	my @lFilters;
	foreach my $filter (split( $separator, $filters )) {
		my ($chrId, $start, $end, $type) = split( ':', $filter );
		$filter = "$chrId:$start:$end:exclude" if ($type eq '-1');
		push(@{$hashFilters->{$chrId}}, $filter);
		$hashChromosomes->{$chrId} = 0;
	}
	$self->filter_region( $hashFilters );
	$self->filter_chromosome( $hashChromosomes );
}

# methode pour construire l'arbe de gene atlas
sub construct_tree {
	my $self = shift;
	my $file2 = $self->buffer->{config}->{gene_atlas_diseases}->{dir};
	$file2 .= $self->buffer->{config}->{gene_atlas_diseases}->{tree};
	my $tree = retrieve($file2);
	my @List;
	my %root;
	$root{name} = "GeneAtlas Diseases";
	$root{type} = "root";
	$root{id}   = "root";
	$root{children} = \@List;
	$root{item} = ();
	foreach my $cat (sort {$a->{name} cmp $b->{name}} values %{$tree}){
		my $item = $self->add_children($cat);
		next unless $item;
		$item->{type} = "category";
		push(@{$root{children}}, $item);
	}
	return [\%root];
}

# methode utilise pour la creation d'un embranchement de l'arbre de gene atlas
sub add_children {
	my ($self, $cat) = @_;
	my $id = $cat->{id};
	my $genes_cat    = $self->atlas_genes_cat();
	my $patients_cat = $self->atlas_patients_cat();
	return  unless exists $genes_cat->{$id}; 
	my %item;
	$self->change_category_name($cat, $genes_cat, $patients_cat);
	$item{name} = $cat->{name};	
	$item{sid}  = $cat->{id};
	$item{id}   = $self->{uuid}++;
	$item{type} = $cat->{type};
	return \%item unless exists $cat->{children};
	foreach my $tt (values %{$cat->{children}}){
		my $its = $self->add_children($tt);
		push(@{$item{children}},$its) if $its;
	}
	return \%item;
}

# methode pour donner le nom de la categorie d'un sous ensemble de l'arbre de gene atlas
sub change_category_name {
	my ($self, $cat, $genes_cat, $patients_cat) = @_;
	my $id = $cat->{id};
	$cat->{name} .= "(g:".$genes_cat->{$id}." p:".scalar(keys %{$patients_cat->{$id}}).")"	;
}

sub _newVariant_cache{
	my ($self, $id, $v_id) = @_;
	my $hash;
	my $type;
	my $strucType;
	my ($chr_name,$start,$ref,$alt) = split('_', $id);
	$hash->{id} = $id;
	$hash->{annex} = undef;
	$hash->{line_infos} = "";
	$hash->{start} = $start;
	#if ($db_idpos->get($id)) { $hash->{start} = $db_idpos->get($id); }
	#else { warn "Problem... ".$id; }
	#$db_idpos->close();
	$hash->{end} = $hash->{start};
	$hash->{strand} = 1;
	$hash->{vcf_id} = $id;
	
	if (length($ref) == length($alt)) {
		$type = "variations";
		$strucType = "snp";
		$hash->{ref_allele} = $ref;
		$hash->{var_allele} = $alt;
		$hash->{isVariation} = 1;
	}
	elsif (length($ref) > length($alt)) {
		$type = "deletions";
		$strucType = "del";
		my $t = 0;
		my @lRef = split("", $ref);
		while ($t < length($alt)) {
			shift(@lRef);
			$t++;
		}
		$hash->{ref_allele} = join("", @lRef);
		$hash->{var_allele} = "-";
		$hash->{end} = $hash->{start}+length($hash->{ref_allele})-1;
		$hash->{isDeletion} = 1;
	}
	elsif (length($ref) < length($alt)) {
		$type = "insertions";
		$strucType = "ins";
		my $t = 0;
		my @lAlt = split("", $alt);
		while ($t < length($ref)) {
			shift(@lAlt);
			$t++;
		}
		$hash->{ref_allele} = $lAlt[0];
		$hash->{var_allele} = join("", @lAlt);
		$hash->{end} = $start;
		$hash->{isInsertion} = 1;
	}
	$hash->{structuralTypeObject} = $type;
	$hash->{structuralType} = $strucType;
	my $chr = $self->getChromosome($chr_name);
	$hash->{chromosomes_object} = { $chr->id() => undef };
	$hash->{chromosome} = $chr;
	$hash->{vector_id} = $v_id;
	my $obj = $self->flushObject($type, $hash);
}

sub getVariantFromId {
	my ($self, $id) = @_;
	my @lFieldsId = split('_', $id);
	my $chrName = $lFieldsId[0];
	my $refAll = $lFieldsId[2];
	my $varAll = $lFieldsId[3];
	if (length($refAll) == length($varAll))   {
		return $self->getChromosome($chrName)->{variations_object}->{$id};
	}
	elsif (length($refAll) < length($varAll)) {
		return $self->getChromosome($chrName)->{insertions_object}->{$id};
	}
	elsif (length($refAll) > length($varAll)) {
		return $self->getChromosome($chrName)->{deletions_object}->{$id};
	}
	die("\n\nERROR: no variant found with ID $id. Exit...\n\n");
}

sub setVariants {
	my ($self, $type) = @_;
	
	my $method;
	if ($type eq 'variations')         { $method = 'getVariations'; }
	elsif ($type eq 'insertions')      { $method = 'getInsertions'; }
	elsif ($type eq 'deletions')       { $method = 'getDeletions'; }
	elsif ($type eq 'large_deletions') { $method = 'getLargeDeletions'; }
	elsif ($type eq 'large_duplications') { $method = 'getLargeDuplications'; }
	my $h;
	foreach my $chr (@{$self->getChromosomes()}) {
		next if ($chr->not_used());
		foreach my $var (@{$chr->$method()}) {
			$h->{$var->id()} = undef;
		}
	}

	return $h;
}




sub setVariations {
	my $self = shift;
	return $self->setVariants('variations');
}

sub setInsertions {
	my $self = shift;
	return $self->setVariants('insertions');
}

sub setDeletions {
	my $self = shift;
	return $self->setVariants('deletions');
}

sub setLargeDeletions {
	my $self = shift;
	return $self->setVariants('large_deletions');
}


sub setGenes {
	my $self = shift;
	my $hGenes;
	foreach my $chr (@{$self->getChromosomes()}) {
		next if ($chr->not_used());
		foreach my $g (@{$chr->getGenes()}) {
			$hGenes->{$g->id()} = undef;
		} 
	}
	return $hGenes;
}

#setBundleMethods : use this method to limit the getTranscripts 
#trick is to increment the  genes_object on the chromosome object so You will never call setgenes
# and same thing for transcripts_object 
#may be I have to test 

 sub setBundle {
 	my $self = shift;
 	my @transcripts_cgi = @{$self->bundle_transcripts() } ;
 	
 	 $self->setListTranscripts(\@transcripts_cgi);

 }
 
 sub setListTranscripts {
 	my ($self,$list) = @_;
 	confess() unless $list;
 	delete   $self->{genes_object};
 	delete  $self->{transcripts_object};
 		foreach my $tr (@$list) {
 			my  $tr1 = $self->newTranscript($tr);
		$self->{transcripts_object}->{$tr1->id} = undef;
		 my $gene = $tr1->getGene();

		  $self->{genes_object}->{$gene->id} ++;
		 $gene->getChromosome()->{genes_object}->{$gene->id} = undef;
		  $gene->getChromosome()->{transcripts_object}->{$tr1->id} = undef;
		 $gene->{transcripts_object}->{$tr1->id} =undef;
 		}
 }
 
sub _newVariant{
	my ($self,$id,$patient) = @_;
	my $can_create_variant = 1;
	my ($chr_name,$start,$ref,$alt) = split('_', $id);
	my $chr = $self->getChromosome($chr_name);
	my $obj =  $chr->flushObjectVariantCache($id, $can_create_variant);
	unless ($obj) {
#		confess();
		$obj = $self->SUPER::_newVariant($id);
	}
	#die();
	return $obj;
	confess();
	
}

sub getVariant{
	my ($self,$id,$patient) = @_;
	my ($chr_name,$start,$ref,$alt) = split('_', $id);
	my $chr = $self->getChromosome($chr_name);
	return $chr->flushObjectVariantCache($id);
	confess();
}



sub getPrimersByPosition {
	my ($self,$chr,$start,$end) = @_;
	my $tabix = $self->tabix_primers;
	return unless $tabix;
	my $no = $chr->get_lmdb_cnvs("r");

	my $res = $tabix->query_full($chr->ucsc_name,$start,$end);
	return [] unless $res;
		# return {mean=>0,} unless defined $res->{_};
		 my @data;
		my $objs =[];
		 while(my $line = $res->next){
		 	
				my($a,$b,$c,$pid) = split(" ",$line);
				my $o;
				if (exists  $self->{objects}->{primers}->{$pid}){
					 $o = $self->{objects}->{primers}->{$pid};
				}
				else 
				{
				$o = $no->get($pid);
				confess() unless $o;
				bless $o , 'GenBoPrimerCache';
				$o->{project} =  $self;
				$o->{buffer} = $self->buffer;
				}
				
				$self->{objects}->{primers}->{$o->{id}} = $o;
#				$self->{objects}->{primers}->{$o->{id}} = $o;
				push(@$objs,$o);
				
				
		 }
		 return $objs;
}

sub getPrimersByObjects {
	my ($self,$obj) = @_;
	my $chr = $obj->getChromosome();
	my $no = $chr->get_lmdb_cnvs("r");
	my $tabix = $self->tabix_primers;
	return unless $tabix;
	my $res = $tabix->query_full($chr->ucsc_name,$obj->start,$obj->end);
	
	#return [] unless $res->get();
		# return {mean=>0,} unless defined $res->{_};
		 my @data;
		my $objs =[];
		 while(my $line = $res->next){
				my($a,$b,$c,$pid) = split(" ",$line);
				my $o;
				if (exists  $self->{objects}->{primers}->{$pid}){
					 $o = $self->{objects}->{primers}->{$pid};
				}
				else 
				{
				$o = $no->get($pid);
				confess() unless $o;
				bless $o , 'GenBoPrimerCache';
				$o->{project} =  $self;
				$o->{buffer} = $self->buffer;
				}
				
			
				$obj->{primers_object}->{$o->id} = 0;
				$self->{objects}->{primers}->{$o->{id}} = $o;
				push(@$objs,$o);
				
				
		 }
		 return $objs;
	
	
}

#Patrick add 

sub get_only_list_patients {
	my ($self,$patients_name,$separator) = @_;
	$separator = "," unless $separator;
	my $patients;
	if (not $patients_name or $patients_name eq 'all'){
		$patients = $self->getPatients();
	}
	else {
		my %names;
		map{$names{$_}++} split($separator,$patients_name);
		
		foreach my $patient (@{$self->getPatients()}){
			if (exists $names{$patient->name}){
				push(@$patients,$patient);
				delete $names{$patient->name};
			}
			else {
				delete $self->{objects}->{patients}->{$patient->id};
				delete $self->patients_object->{$patient->id};
			}
			
		}
			
		}
		my $v ;
		foreach my $chr (@{$self->getChromosomes}){
			#my $v1 = $chr->variants();
			my $v1;
			foreach my $patient (@$patients){
				unless ($v1){
					$v1 =   $patient->getVectorOrigin($chr)->Clone;
				}
				else {
					$v1 += $patient->getVectorOrigin($chr);
				}
			}
			$chr->variants($v1);
		}
	return $patients;
}


sub returnVariants {
		my ($self, $id, $type) = @_;
		my ($chr_name,$vid) = split("!",$id);
		my $chr = $self->getChromosome($chr_name);
		 my $gid = $chr->cache_lmdb_variations->get_varid($vid);
		 #warn $gid." ".$vid." ".$id;
		 my $obj = $chr->cache_lmdb_variations->get($gid);
		  $obj->{global_vector_id} = $id;
		 $obj->{vector_id} = $vid;
		 
		 $obj->{project} =  $self;
		$obj->{buffer} = $self->buffer;
		return $obj;
		
}

sub setListVariants {
	my ($self,$list) = @_;
	unless ($list){
		foreach my $chr (@{$self->getChromosomes}){
			my $vector = $chr->getVariantsVector();
			my $set = Set::IntSpan::Fast::XS->new($vector->to_Enum);
			my $iter = $set->iterate_runs();
			while (my ( $from, $to ) = $iter->()) {
   				for my $member ($from .. $to) {
   					push(@$list,$chr->name."!".$member);
   				}
    		}
		}
	}
		$self->{list_variant} = $list;
	
	
}
sub getVariantsList {
	my ($self,$chr_query) = @_;
	my $list;

	
	foreach my $chr (@{$self->getChromosomes}){
			if ($chr_query) {
				next if $chr->name ne $chr_query->name;
			}
			my $vector = $chr->getVariantsVector();
			my $set = Set::IntSpan::Fast::XS->new($vector->to_Enum);
			my $iter = $set->iterate_runs();
			while (my ( $from, $to ) = $iter->()) {
   				for my $member ($from .. $to) {
   					push(@$list,$chr->name."!".$member);
   				}
    		}
		}
		return $list;
}

sub nextVariant {
	my ($self) = @_;
	$self->setListVariants() unless exists $self->{list_variant};
	return unless scalar(@{$self->{list_variant}});
	my $id = shift(@{$self->{list_variant}} );
	my $var_obj = $self->returnVariants($id);
	
	my $ref = ref($var_obj);
	if ($ref eq 'GenBoVariation'){
					bless $var_obj , 'GenBoVariationCache';
		}
		elsif  ($ref eq 'GenBoLargeDeletion'){
					bless $var_obj , 'GenBoLargeDeletionCache';
					
		}
		elsif  ($ref eq 'GenBoLargeInsertion'){
					bless $var_obj , 'GenBoLargeDuplicationCache';
		}
		elsif  ($ref eq 'GenBoDeletion'){
					bless $var_obj , 'GenBoDeletionCache'; 
		}
		elsif  ($ref eq 'GenBoInsertion'){
					bless $var_obj , 'GenBoInsertionCache';
		}
		elsif  ($ref eq 'GenBoLargeDuplication'){
					bless $var_obj , 'GenBoLargeDuplicationCache';
		}
		elsif  ($ref eq 'GenBoLargeDeletion'){
					bless $var_obj , 'GenBoLargeDeletionCache';
		}
				
		return $var_obj;
	
}




sub myflushobjects {
	my ($self, $ids, $type) = @_;

	my $array_ids;
	if (ref($ids) eq 'HASH' ){
		if (exists $ids->{none}) {
			return [];
		}
		$array_ids =[keys %$ids];
	} 
	elsif  (ref($ids) eq 'ARRAY' ){
			
		$array_ids = $ids;
	}
	else {
		return[];
		confess($array_ids);
	}
	my @objs = (); 
		foreach  my $id (@$array_ids) {
			
			confess("id".$id) unless $id;
			if (exists $self->{lmdb_id}->{$id}) {
				$id =  $self->{lmdb_id}->{$id};
			}
			unless (exists  $self->{objects}->{$type}->{$id}){
				if ($type =~ /gene/ or $type =~ /transcript/ or  $type =~/protein/ or $type =~ /exon/ or $type  =~ /intron/) {
					
					my $obj = $self->lmdbGenBo->get( $id );
					unless ($obj) {
						my $syno =  $self->liteAnnotations->get("synonyms", $id);
						$obj = $self->lmdbGenBo->get( $syno );
					}
					#cas avec des transcript id sans _ impossible a retrouver
					confess($id." ".$type) unless $obj;
					bless $obj , 'GenBoGeneCache' if ($type eq "genes");
					bless $obj , 'GenBoTranscriptCache' if ($type eq "transcripts");
					$obj->{project} =  $self;
					$obj->{buffer} = $self->buffer;
					$self->{objects}->{$type}->{$id} = $obj;
					
				}
				elsif ($type =~ /variant/ or $type eq 'variations' or $type eq 'deletions' or $type eq 'insertions' or $type eq 'large_duplications' or $type eq 'large_duplications' ){
					
					#confess() if 
					my $vector_id;
					my $chr;
					my $real_id;
					if ($id =~/!/){
						my ($chr_name,$vid) = split("!",$id);
						$real_id = $id;
						$vector_id = $vid;
						 $chr = $self->getChromosome($chr_name);
						unless (exists $self->{lmdb_id}->{$id}){
						$self->{lmdb_id}->{$id} =  $chr->cache_lmdb_variations->get_varid($vid);
					#$id = $chr->cache_lmdb_variations->get_varid($vid);
						}
						$id = $self->{lmdb_id}->{$id};
					}
					else {
					confess("bizarre comme id $id");
					my ($chr_name,$vid) = split("_",$id);
				
					 $chr = $self->getChromosome($chr_name);
					}
					my $var_obj = $chr->cache_lmdb_variations->get($id);
					
					my $ref = ref($var_obj);
					$var_obj->{vector_id}= $vector_id;
				if ($ref eq 'GenBoVariation'){
					bless $var_obj , 'GenBoVariationCache';
					$self->{objects}->{variations}->{$id}= $var_obj;
				}
				elsif  ($ref eq 'GenBoLargeDeletion'){
					bless $var_obj , 'GenBoLargeDeletionCache';
					$self->{objects}->{large_deletions}->{$id}= $var_obj;
					
				}
				elsif  ($ref eq 'GenBoLargeInsertion'){
					bless $var_obj , 'GenBoLargeDuplicationCache';
					$self->{objects}->{large_duplications}->{$id}= $var_obj;
				}
				elsif  ($ref eq 'GenBoDeletion'){
					bless $var_obj , 'GenBoDeletionCache'; 
					$self->{objects}->{deletions}->{$id}= $var_obj;
				}
				elsif  ($ref eq 'GenBoInsertion'){
					bless $var_obj , 'GenBoInsertionCache';
					$self->{objects}->{insertions}->{$id}= $var_obj;
				}
					elsif  ($ref eq 'GenBoLargeDuplication'){
					bless $var_obj , 'GenBoLargeDuplicationCache';
					$self->{objects}->{insertions}->{$id}= $var_obj;
				}
					elsif  ($ref eq 'GenBoLargeDeletion'){
					bless $var_obj , 'GenBoLargeDeletionCache';
					$self->{objects}->{insertions}->{$id}= $var_obj;
				}
				elsif  ($ref ne 'GenBoVariationCache' &&  $ref ne 'GenBoInsertionCache' && $ref ne 'GenBoDeletionCache' && $ref ne 'GenBoLargeDuplicationCache' && $ref ne 'GenBoLargeDeletionCache') {
					confess("$ref =+>".$var_obj);
				}
				$var_obj->{project} =  $self;
				$var_obj->{buffer} = $self->buffer;
				
				$self->{objects}->{$type}->{$id}= $var_obj;
				
				}
				elsif ($type eq 'runs') {$self->getRunFromId($id); }
				elsif ($type eq 'patients') {
					$self->setPatients(); 
					confess($id) unless exists $self->{objects}->{$type}->{$id};
				}
				elsif ($type eq 'captures') {
					#$self->setCaptures(); 
					$self->createObject($type,{id=>$id});
					confess() unless exists $self->{objects}->{$type}->{$id};
				}
				elsif ($type eq 'panels') {
					$self->createObject($type,{id=>$id});
					confess() unless exists $self->{objects}->{$type}->{$id};
				}
				elsif ($type eq 'bundles') {
					$self->createObject($type,{id=>$id});
					confess() unless exists $self->{objects}->{$type}->{$id};
				}
				elsif ($type eq 'primers') {
					$self->createObject($type,{id=>$id});
					confess() unless exists $self->{objects}->{$type}->{$id};
				}
				elsif ($type eq 'phenotypes') {
					$self->createObject($type,{id=>$id});
					confess() unless exists $self->{objects}->{$type}->{$id};
				}
				else {
					confess("je fais quoi ici $type");
				}
				
			}
			else {
				#warn "coucou" if $type eq "variants";
			}
			#confess($id) unless (exists  $self->{objects}->{$type}->{$id});
			push(@objs, $self->{objects}->{$type}->{$id});

		}
	#}
	
	return \@objs;
}

has interval_tree_vector_transcripts => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		warn $self->dir_lite_cache_tree_objects();
		my $no = GenBoNoSqlIntervalTree->new(dir=>$self->dir_lite_cache_tree_objects(),mode=>"r");
		return $no;
	},
);




1;
