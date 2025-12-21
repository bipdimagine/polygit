package GenBoVariant;

use strict;
use Moo;

use FindBin qw($Bin);
use lib "$Bin/packages/";

use Parallel::ForkManager;
use Data::Dumper;
use Config::Std;
use GenBoCapture;
use GenBoGenomic;
use Position;
use List::Util;
use JSON::XS;
use Carp;
use Compress::Snappy;
 use List::Util qw( max min sum);
use Storable qw/thaw freeze/;
use Carp;
use liftOver;
use MIME::Base64;
#use bytes;
extends "GenBoGenomic";



has isVariant => (
	is		=> 'ro',
	default	=> 1,
);

has isLarge => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> 0
);

has isCnv => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> 0
);

# id in vector
has vector_id => (
	is	=> 'ro',
	lazy => 1,
	default => undef,
);

has isVariation => (
	is		=> 'ro',
	default	=> 0,
);

has isLargeDeletion => (
	is		=> 'ro',
	
	default	=> 0,	
);
#
#has isLargeDuplication => (
#	is		=> 'rw',
#	default	=> 0,
#);


has vcf_position => (
	is		=> 'rw',
);
has vcf_first_base => (
	is		=> 'rw',
	
);
#has vcf_sequence => (
#	is		=> 'rw',
#);

has ref_allele => (
	is		=> 'rw',
	#required=> 1,
);

has var_allele => (
	is		=> 'rw',
	#required=> 1,
);

has structuralType => (
	is		=> 'rw',
	reader	=> 'getStructuralType',
	#required=> 1,
);

has validation_sanger => (
	is		=> 'rw',
	default	=> "0",
);

has validation_method => (                
	is		=> 'rw',
	default	=> "",
);
has validation_ngs => (
	is		=> 'rw',
	default	=> "0",
);

has is_in_pseudoautosomal_region=> (
 is		=> 'rw',
lazy =>1,
default => sub {
	my $self = shift;
	my $chr = $self->getChromosome();
	return $chr->isPseudoAutosomal($self->start,$self->end);
}
);

# lift_over

sub lift_over {
	my ($self,$version) = @_;
	confess() unless $version;
	my $origin_version = $self->project->genome_version_generic();
	if ($origin_version eq $version) {
		die();
	}
	my $key = "lift_over_$version";
	return $self->{$key} if exists $self->{$key};
	my $lift = liftOver->new(project=>$self->project,version=>$version);
	$lift->lift_over_variant($self, $version);
	return $self->{$key};
	return $self->{$key} unless $self->{$key};
	return $self->{$key};
}

#################
#
# public data
#
##################

#############
#SPLICE AI
#############




sub spliceAI {
	my ($self) = @_;
	return $self->{spliceAI_rocks} if exists $self->{spliceAI_rocks};
	$self->{spliceAI_rocks} = {};
	my $spliceAI =  $self->getChromosome->rocksdb("spliceAI")->spliceAI($self->rocksdb_id); 
	die() if $spliceAI && $spliceAI eq "-";
	if ($spliceAI){
		$self->{spliceAI_rocks} =  $spliceAI;
	}
	return $self->{spliceAI_rocks};
	
}

sub spliceAI_score {
	my ($self,$gene,$debug) = @_;
	delete $self->{spliceAI_rocks};
	return {} unless exists $self->spliceAI->{uc($gene->external_name)};
	my $sp = $self->spliceAI->{uc($gene->external_name)};
	my $hash = {};
	my $score = ["AG","AL","DG","DL"];
	foreach my $s (@$score){
			$hash->{$s} = $sp->{$s};
	}

	return $hash;
	
	
}


sub alphamissense {
	return "-";
}
sub max_spliceAI_score {
	my ($self,$gene) = @_;
#	warn $self->id." ".$self->project->name." ".$self->spliceAI;
#	confess() if $self->spliceAI eq 0;
	return 0 unless exists $self->spliceAI->{uc($gene->external_name)};
	return $self->spliceAI->{uc($gene->external_name)}->{max};
	#return $self->getChromosome()->score_gene_spliceAI($v,$gene->external_name);
}

sub max_spliceAI_categorie {
	my ($self,$gene) = @_;
	return "-" unless exists $self->spliceAI->{uc($gene->external_name)};
	return $self->spliceAI->{uc($gene->external_name)}->{max_cat};
}

sub text_max_spliceAI {
	my ($self,$gene) = @_;
	return "-" unless exists $self->spliceAI->{uc($gene->external_name)};
	return $self->max_spliceAI_categorie($gene).":".$self->max_spliceAI_score($gene)
}




sub setFamilies {
	my $self = shift;
	my $hash ={};
	foreach my $patient (@{$self->getPatients}){
		my $fam = $patient->getFamily;
		$hash->{$fam->id} = undef if $fam;
		
	}
	return $hash;
}


sub return_extension_categorie_frequency {
	my ($self,$value) = @_;
	my $text = "none";
	return $text unless $value ;
	$text = "none" unless $value;
	if    ($value eq "-")     { $text = "none"; }
	elsif    ($value == -1)     { $text = "none"; }
		elsif ($value <= 0.0001) { $text = "0001"; }
		elsif ($value <= 0.001)  { $text = "001"; }
		elsif ($value <= 0.01)   { $text = "01"; }
		elsif ($value <= 0.05)   { $text = "05"; }
		else                    { $text = "1"; }		
		return $text;
		
}




sub alternate_allele {
	my ($self )= @_;
	return $self->sequence();
} 




#########
## HGMD
#######
has hgmd => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return unless $self->project->isHGMD;
		my $db =  $self->getChromosome->rocksdb("hgmd");
		return $db->hgmd($self->rocksdb_id);
		my $pub = $db->get_with_sequence($self->start,$self->alternate_allele);
		return $pub;
	}
);

has isNewHgmd => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return unless ($self->hgmd_id());
		return if ($self->project->buffer->queryHgmd->hashOldAccNum and exists $self->project->buffer->queryHgmd->hashOldAccNum->{$self->hgmd_id()});
		return 1;
	}
);

has hgmd_details => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->project->buffer->queryHgmd->getDataHGMDPro($self->hgmd_id()) if $self->hgmd_id;
		return;
	}
);

has hgmd_id => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $pub = $self->hgmd();
		if ($pub){
			return $pub->{hgmd_id};
		}
		return undef;
	},
);

has hgmd_class => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $pub = $self->hgmd();
		if ($pub){
			return $pub->{class};
		}
		return;
	},
);


has isDM => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return 	unless $self->project->isHGMD;
		return unless $self->hgmd_id();
		return 1 if ($self->hgmd_class() eq 'DM');
		return;
	},
);

has genes_pathogenic_DM => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hash = {};
		my $chr = $self->getChromosome;
		foreach my $gene (@{$self->getGenes}){
			 $hash->{$gene->id}->{DM} = undef;
			 $hash->{$gene->id}->{pathogenic} = undef;
			 if ($self->isDM && !($self->project->isHGMD)){
				eval { $hash->{$gene->id}->{DM} = $chr->is_hgmd_DM_for_gene($self->hgmd_id(), $gene); };
				if ($@) { $hash->{$gene->id}->{DM} = 1; }
			 }
			 if ($self->isClinvarPathogenic && $self->clinvar_id) {
				eval {  $hash->{$gene->id}->{pathogenic} = $self->is_clinvar_pathogenic_for_gene($self->clinvar_id(), $gene); };
				if ($@) { $hash->{$gene->id}->{pathogenic} = 1; }
			}
		}
		
		return $hash;
	},
);

sub isDM_for_gene {
	my ($self, $gene) = @_;
	return $self->genes_pathogenic_DM->{$gene->id}->{DM};
#	return unless ($self->isDM);
#	confess unless ($gene);
#	my $res;
#	eval { $res = $self->getChromosome->is_hgmd_DM_for_gene($self->hgmd_id(), $gene); };
#	if ($@) { $res = 1; }
#	return $res;
}



has hgmd_inheritance => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->hgmd_details->{inheritance} if (exists $self->hgmd_details->{inheritance});
		return $self->hgmd_details->{expected_inheritance} if (exists $self->hgmd_details->{expected_inheritance});
		confess("\n\nERROR: HGMD key inheritance not found. Die.\n\n");
	},
);

has hgmd_hgvs => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->hgmd_details->{hgvs} if (exists $self->hgmd_details->{hgvs});
		confess("\n\nERROR: HGMD key hgvs not found. Die.\n\n");
	},
);


		
has hgmd_phenotype => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->hgmd->{phen} if ($self->hgmd && exists $self->hgmd->{phen});
		confess("\n\nERROR: HGMD key phen not found. Die.\n\n");
	},
);

has hgmd_disease => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->hgmd_details->{disease} if (exists $self->hgmd_details->{disease});
		confess("\n\nERROR: HGMD key disease not found. Die.\n\n");
	},
);


has hgmd_releases => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my @lReleases;
		push(@lReleases, $self->project->buffer->queryHgmd->database());
		foreach my $hgmd_version (keys %{$self->project->buffer->queryHgmd->getHashOldDatabases()}) {
			my $h_details = $self->project->buffer->queryHgmd->getDataHGMDPro($self->hgmd_id(), $hgmd_version);
			push(@lReleases, $hgmd_version) if ($h_details);
		}
		my $releases = join(', ', reverse sort (@lReleases));
		return $releases;
	},
);

has hgmd_pubmed => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $hPubmed;
		$hPubmed->{$self->hgmd_details->{pmid}}->{title} = $self->hgmd_details->{title};
		$hPubmed->{$self->hgmd_details->{pmid}}->{author} = $self->hgmd_details->{author};
		$hPubmed->{$self->hgmd_details->{pmid}}->{year} = $self->hgmd_details->{year};
		$hPubmed->{$self->hgmd_details->{pmid}}->{fullname} = $self->hgmd_details->{fullname};
		if (exists $self->hgmd_details->{pubmed}) {
			foreach my $pubmed_id (keys %{$self->hgmd_details->{pubmed}}) {
				foreach my $cat (keys %{$self->hgmd_details->{pubmed}->{$pubmed_id}}) {
					$hPubmed->{$pubmed_id}->{$cat} = $self->hgmd_details->{pubmed}->{$pubmed_id}->{$cat};
				}
			}
		}
		return $hPubmed;
	},
);

has hgmd_gene_name => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $db =  $self->getChromosome->get_lmdb_database("hgmd",'relation_variant_gene');
		my $pub = $db->get($self->hgmd_id());
		my $gene_name = join(',', sort keys %$pub);
		return $gene_name;
	},
);

#####¥###########
## clinvar
############### 

has clinvar => (
		is		=> 'rw',
		lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $db =   $self->getChromosome->rocksdb("clinvar");
		my $pub = $db->clinvar($self->rocksdb_id);
		$pub ={} unless $pub;
		return $pub;
	}
);

sub isClinvarPathogenic {
	my $self = shift;
	return 1 if $self->clinvar and exists $self->clinvar->{score} and $self->clinvar->{score} >= 5;
	return;
}


sub clinvar_id {
	my $self = shift;
	return $self->clinvar->{clinvar_id} if exists $self->clinvar->{clinvar_id};
	return;
}

sub score_clinvar {
	my $self = shift;
	return $self->clinvar->{score} if exists $self->clinvar->{score};
	return;
}

sub text_clinvar {
	my ($self) = @_;
	return $self->clinvar->{clnsig} if (exists $self->clinvar->{clnsig});
	return;
}
sub is_clinvar_pathogenic_for_gene {
	my ($self, $gene) = @_;
	my $name = $gene->external_name;
	return unless $self->isClinvarPathogenic;
	return $self->clinvar->{genes}->{$name} if exists $self->clinvar->{genes}->{$name};
	return;
}

sub isClinvarPathogenic_for_gene {
	my ($self, $gene) = @_;
	confess();
}
has clinvar_class => (
		is		=> 'rw',
		lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->clinvar->{clnsig};
		confess();
		my $v = $self->score_clinvar();
		return "" unless $v;
		if ($v == 0){
	  		 return "Uncertain"
	  	}
	  	elsif ($v ==1 ){return "drug response"}
	  	elsif ($v ==2 ){return  "Benign"}
	  	elsif ($v ==3 ){return  "Likely benign"}
	  	elsif ($v ==4 ){return "Likely pathogenic";}
	  	elsif ($v ==5 ){return  "pathogenic"	}
	  	elsif ($v ==-2 ){return  "uncertain_significance";}
	  	elsif ($v ==-1 ){return  "not_provided";}
	  	return  "other";
	},
);


#clinical local database

has clinical_local => (
		is		=> 'rw',
		lazy	=> 1,
	default	=> sub {
	my $self = shift;
	return undef;
	my $db =  $self->buffer->get_lmdb_database("cldb",$self->getChromosome->name,$self->type_public_db);
	my $pub =  $db->get_with_sequence($self->start,$self->alternate_allele);
	#my $db1 =  $self->buffer->get_lmdb_clinical_local_db("cldb");
	#my $pub = $db1->get($self->id);
	#$pub = $db->get_with_sequence($self->start,$self->alternate_allele) if $pub; 
	
	return $pub;
	}
);


has score_clinical_local => (
		is		=> 'rw',
		lazy	=> 1,
	default	=> sub {
	my $self = shift;
	my $pub = $self->clinical_local;
	return 5 if $pub;
	my $v = 0;
	if ($pub){
		$pub->{sig} =~s/,/;/g;
	   	my $v = max( split(";",$pub->{sig}));
	   	return $v;
	}
	return undef;
	},
);

has pmid_local => (
		is		=> 'rw',
		lazy	=> 1,
	default	=> sub {
	my $self = shift;
	my $pub =  $self->clinical_local;
	return ""  unless $pub;
	return ""  unless $pub->{PMID};
	my $id = $pub->{PMID};
	$id =~s/PMID://;
	return $id
	
	},
);


has comment_clinical_local => (
		is		=> 'rw',
		lazy	=> 1,
	default	=> sub {
	my $self = shift;
	my $pub =  $self->clinical_local;
	return $pub->{comment} if $pub;
	return "";
	
	},
);

has text_clinical_local => (
		is		=> 'rw',
		lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $v = $self->score_clinical_local();
		return "" unless $v;
		if ($v == 0){
	  		 return "Uncertain"
	  	}
	  	elsif ($v ==1 ){return "drug response"}
	  	elsif ($v ==2 ){return  "Benign"}
	  	elsif ($v ==3 ){return  "Likely benign"}
	  	elsif ($v ==4 ){return "Likely pathogenic";}
	  	elsif ($v ==5 ){return  " pathogenic"	}
	  	elsif ($v ==-2 ){return  " uncertain_significance";}
	  	elsif ($v ==-1 ){return  " not_provided";}
	  	return  "other";
	},
	);


###################
# GNOMAD METHODS
###################



sub gnomad {
	my ($self,$type) = @_;
	
	unless (exists $self->{gnomad_rocks}){
		if (exists $self->{gnomad} && %{$self->{gnomad}}){
		 my $a = delete $self->{gnomad};
					 
		 	$self->{gnomad_rocks}->{min_pop} = $a->{minpop};
			 $self->{gnomad_rocks}->{min} = $a->{populations}->{$a->{minpop}}->{F};
		 	$self->{gnomad_rocks}->{max_pop} = $a->{maxpop}; ;
		 	$self->{gnomad_rocks}->{public} = $a->{populations}->{$a->{maxpop}}->{F};
			$self->{gnomad_rocks}->{xy}= $a->{populations}->{all}->{AC_male};
		 	$self->{gnomad_rocks}->{an} = $a->{populations}->{all}->{AN};
		 	$self->{gnomad_rocks}->{ac} = $a->{populations}->{all}->{AC};
		 	$self->{gnomad_rocks}->{ho} = $a->{populations}->{all}->{HO};
		 	$self->{gnomad_rocks}->{rs} = $a->{rsname};
		 	$self->{gnomad_rocks}->{public} = 1;
		 	$a =undef;
			return $self->{gnomad_rocks};
		}
	else {
		my $gnomad =  $self->getChromosome->rocksdb("gnomad")->value($self->rocksdb_id); 
			$self->{gnomad_rocks} = {} ;
		if ($gnomad){
			$self->{gnomad_rocks} =  $gnomad;
			$self->{gnomad_rocks}->{public} = 1;
		
		}
		}
		
	}
	return if $type eq "cache";
	confess("no type for gnomad od call with cache") unless $type;
	return $self->{gnomad_rocks}->{$type};
	
}
sub min_pop_name {
	my ($self) = @_;
	return undef unless  $self->gnomad("public");
	return $self->gnomad("min_pop");
}

sub min_pop_freq {
	my ($self) = @_;
	return "-" unless  $self->gnomad("public");
	return $self->gnomad("min");
}


sub max_pop_freq {
	my ($self) = @_;
	return "-" unless  $self->gnomad("public");
	return $self->gnomad("max");
}
sub max_pop_name {
	my ($self) = @_;
	return undef unless  $self->gnomad("public");
	return $self->gnomad("max_pop");
}


sub isPublic {
	my ($self) = @_;
	return undef unless  $self->gnomad("public");
	return $self->gnomad("public");
}

sub getGnomadHO {
	my ($self) = @_;
	return undef unless  $self->gnomad("public");
	#return undef if  $self->gnomad("ho") eq "-";
	unless  ($self->getChromosome->isAutosomal($self->start,$self->end)){
		return $self->gnomad("xy");
	}
	return $self->gnomad("ho");
}
sub getGnomadAC {
	my ($self) = @_;
	return undef unless  $self->gnomad("public");
	return $self->gnomad("ac");
}
sub getGnomadAN {
	my ($self) = @_;
	return $self->gnomad("an");
}
sub getGnomadAC_Male {
	my ($self) = @_;	
	return undef unless  $self->gnomad("public");
	return $self->gnomad("xy");
}
sub getGnomadAC_XY {
	my ($self) = @_;	
	return undef unless  $self->gnomad("public");
	return $self->gnomad("xy");
}
sub frequency {
	my ($self) = @_;
	return undef unless  $self->gnomad("public");
	return undef if $self->getGnomadAN == 0;
	return ($self->getGnomadAC/$self->getGnomadAN);
}

sub frequency_homozygote {
	my ($self) = @_;	
	return undef unless  $self->gnomad("public");
	return undef  if $self->getGnomadAN == 0;
		return 0 if $self->getGnomadAN == 0;
	return ($self->getGnomadHO()*2 / $self->getGnomadAN());
}
sub rs_name {
	my ($self) = @_;	
	return "-" unless  $self->gnomad("public");
	return "rs".$self->gnomad("rs");
}



###################
## END GNOMAD
##################

has isClinical =>(
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $sc = 0;
		$sc = $self->score_clinvar() if $self->score_clinvar();
		
		return 1 if $sc > 3 or $self->clinical_local;
	#	return 1 if $self->public_data()->[2] &   $self->buffer()->mask_database("clinvar");
	}
);


has hotspot => (
	is		=> 'ro',
	
	default	=> "0",
		lazy	=> 1,
	default	=> sub {
		my $self = shift;
			my $hash = $self->getChromosome()->get_lmdb_hotspot($self->type_public_db)->get_with_sequence($self->start,$self->alternate_allele);
			return undef unless $hash ;
			return 1;
	},
);


#has cadd_score => (
#	is		=> 'ro',
#	
#	default => sub {
#		return "-" ;
#	},
#	
#);

has dbscsnv => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my ($self) = @_;
		my $score = $self->getChromosome->rocksdb("dbscSNV")->value($self->rocksdb_id);
	 	return $score if $score;
	 	return {};
	},

);
sub dbscsnv_ada {
	my ($self) = @_;
	return "-" unless exists $self->dbscsnv->{ada};
	return $self->dbscsnv->{ada};
}
sub dbscsnv_rf {
	my ($self) = @_;
	return "-" unless exists $self->dbscsnv->{rf};
	return $self->dbscsnv->{rf};
}


has cadd_score => (
	is		=> 'ro',
	lazy =>1,
	default => sub {
		my $self = shift;
	 	my $score = $self->getChromosome->rocksdb("cadd")->value($self->rocksdb_id);
	 	return $score->{cadd_score} if $score;
	 	return "-";
	 	},
);

has promoterAI => (
	is		=> 'ro',
	lazy =>1,
	default => sub {
		my $self = shift;
		return if $self->project->getVersion() =~ /HG19/;
		return if $self->getChromosome->name() eq "MT";
		my $rocksid = $self->getChromosome->id().'!'.$self->rocksdb_id;
	 	my $res = $self->getChromosome->rocksdb("promoterAI")->get_raw($rocksid);
	 	return if not $res;
	 	my $h;
	 	foreach my $r (split(',', $res)) {
		 	my @l = split(';', $r);
		 	$h->{$l[0]}->{transcript_id} = $l[0];
		 	$h->{$l[0]}->{strand} = $l[1];
		 	$h->{$l[0]}->{tss_score} = $l[2];
		 	$h->{$l[0]}->{score} = $l[3];
	 	}
	 	return $h;
	 },
);

has rocksdb_id => (
	is		=> 'ro',
	lazy => 1,
	default => sub {
		my $self = shift;
		confess();
		return "1-A-A";
		confess($self->gnomad_id);
	},
	
);

sub genomic_rocksdb_id {
	my $self = shift;
	return $self->getChromosome->name."!".$self->rocksdb_id;
	
}

has rocksdb_index => (
   is              => 'rw',
        lazy    => 1,
        default => sub { return "-" ; },
);
has ncboost_score => (
        is              => 'rw',
        lazy    => 1,
        default => sub { return "-" ; },
);

has ncboost_category => (
        is              => 'rw',
        lazy    => 1,
        default => sub { return "-" ; },
);

sub promoterAI_tss_score {
	my ($self, $transcript) = @_;
	confess("\n\nERROR: transcript object mandatory. Die \n\n") if not $transcript;
	return if not $self->promoterAI();
	my ($tr_id, $chr_id) = split('_', $transcript->id());
	return if not exists $self->promoterAI->{$tr_id};
	return $self->promoterAI->{$tr_id}->{tss_score};
}

sub promoterAI_score {
	my ($self, $transcript) = @_;
	confess("\n\nERROR: transcript object mandatory. Die \n\n") if not $transcript;
	return if not $self->promoterAI();
	my ($tr_id, $chr_id) = split('_', $transcript->id());
	return if not exists $self->promoterAI->{$tr_id};
	return $self->promoterAI->{$tr_id}->{score};
}

sub revel_score {
	return "-" ;
}

sub get_infos_database {
	my ($self,$database) = @_;
	die() unless $database;
	my $hash = $self->getChromosome()->get_lmdb_database($database,$self->type_public_db)->get($self->start);
	return undef unless $hash;
	return undef unless exists $hash->{$self->alternate_allele};
	return $hash->{$self->alternate_allele};
}

has population_frequencies=> (
		is		=> 'rw',
	lazy	=> 1,
	default => sub {
			my $self = shift;
			
			return  $self->buffer->get_populations_frequencies($self->getChromosome()->name,$self->type_public_db,$self->start,$self->alternate_allele);
	}
	
	
	
);





has cosmic =>(
	is		=> 'ro',
	lazy=> 1,
	default=> sub {
		my $self = shift;
		my $cosmic =  $self->getChromosome->rocksdb("cosmic")->cosmic($self->rocksdb_id);
		return $cosmic;
		my $hash = $self->getChromosome()->get_lmdb_database("cosmic",$self->type_public_db)->get_with_sequence($self->start,$self->alternate_allele);
		return undef unless $hash;
		$hash->{frequence} ="-" unless $hash->{frequence};
		$hash->{rs} = "" unless exists $hash->{rs};
		$hash->{rsname} = $hash->{rs} unless exists $hash->{rsname};
		return $hash->{rsname}.":".$hash->{frequence};
	},
);

has categorie_frequency => (
	is		=> 'ro',	
	#default	=> "0",
		lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $text = $self->return_extension_categorie_frequency($self->frequency);
		return "freq_".$text;
		
	},
);
#has categorie_frequency_ho=> (
#		is		=> 'rw',
#	
#		lazy	=> 1,
#	default	=> sub {
#		my $self = shift;
#	 	return 0;
#	},
#);


has categorie_frequency_ho=> (
		is		=> 'rw',
	
		lazy	=> 1,
	default	=> sub {
	 	my $self = shift;
	 	my $freq = $self->frequency_homozygote();

		my $text = $self->return_extension_categorie_frequency($freq);
		return "freq_ho_".$text;
	},
);


sub categorie_frequency_prediction { 
			my ($self,$g) = @_;
			die() unless $g;
			my $key_polyphen = "polyphen".$self->polyphenStatus($g);
			my $key_sift = "sift".$self->siftStatus($g);
			my $key_polyphen_sift;
					if (($key_polyphen eq 'polyphen3') and ($key_sift eq 'sift2')) {
						$key_polyphen_sift = 'prediction_3';
					}
					elsif (($key_polyphen eq 'polyphen3') or ($key_sift eq 'sift2')) {
						$key_polyphen_sift = 'prediction_2';
					}
					elsif ($key_polyphen eq 'polyphen2') {
						$key_polyphen_sift = 'prediction_1';
					}
					else {
						$key_polyphen_sift = 'prediction_0';
					}
					return $key_polyphen_sift;
	}
	
has name => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->gnomad_id;
	},
);





has patients_object => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {{}},
);



has kyoto_id => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		confess("rien a faire ici");
	},
);
has gnomad_id => (
	is		=> 'ro',
	lazy=> 1,
	default=> sub {
	my $self = shift;
	my $vn=$self->vcf_id;
	confess unless $vn;
	$vn =~ s/_/-/g;
	$vn=~ s/chr//;
	return $vn;
	}
	);




has vcf_id => (
	is		=> 'rw',
	#required	=> 1,
);

sub theoric_vcf_id {
	my ($self) = @_;
	confess($self);
}
has check_id => (
	is		=> 'rw',
	#required	=> 1,
);

has annotation => (
	is		=> 'ro',
	lazy=> 1,
	default=> sub {
		my $self = shift;
		my $annot = $self->init_annotation();
		my $sum = 0;
		foreach my $cat (keys %{$self->project->maskImpact()}) {
			$sum += $self->project->maskImpact->{$cat};
		}
		$self->sum_mask_coding($sum);
		return $annot
	}
);

has sum_mask_coding => (
	is	    => 'rw',
	lazy    => 1,
	default => 0,
);

has sequence =>(
	is		=> 'ro',
	lazy=> 1,
	default=> sub {my $self = shift; return $self->{var_allele}},
	
);

has alleles => (
	is		=> 'ro',
	lazy	=> 1,
	default=> sub {
		my $self = shift;
		return $self->getChromosome->sequence($self->start(), $self->end())."/".$self->sequence();
	},
);

has alamut_id => (
	is		=> 'ro',
	lazy	=> 1,
	default=> sub {
		my $self = shift;
		my @lTmp = split('/', $self->alleles());
		my $build = 'GRCh37';
		$build = 'GRCh38' if ($self->project->getVersion() =~ /HG38/);
		my $alamut_id = 'chr'.$self->getChromosome->id().':'.$build.':'.$self->start().$lTmp[0].">".$lTmp[1];
		return $alamut_id;
	},
);





#has nb_dejavu =>(
#	is	    => 'ro',
#	lazy    => 1,
#	default => sub {
#		my $self = shift;
#		confess();
#		my $no_dejavu = GenBoNoSqlDejaVu->new( dir => $self->project->deja_vu_lite_dir , mode => "r");
#		my $res = $no_dejavu->get_hohe_nbprojects($self->getChromosome->id(), $self->id());
#		return 0 unless ($res);
#		return $res->[1];
#	},
#);
#
#has nb_dejavu_ho =>(
#	is	    => 'ro',
#	lazy    => 1,
#	default => sub {
#		my $self = shift;
#		confess();
#		my $no_dejavu = GenBoNoSqlDejaVu->new( dir => $self->project->deja_vu_lite_dir , mode => "r");
#		my $res = $no_dejavu->get_hohe_nbprojects($self->getChromosome->id(), $self->id());
#		return 0 unless ($res);
#		return $res->[0];
#	},
#);


####
# method position
###

sub getCodons {
	my ($self,$obj) = @_;
	my $id = $obj->id;
	
	return "" unless exists $self->annotation->{ $id }->{coding}->{sequences};
	return $self->annotation->{ $id }->{coding}->{sequences}->{codon}."/".$self->annotation->{ $id }->{coding}->{sequences}->{codon_mut};
}

sub position {
	my ($self,$obj) = @_;
	confess() unless $obj;
	return $self->{position}->{$obj->id} if exists $self->{position}->{$obj->id};
	if ($obj->isChromosome){
		my $start = $self->start();
		my $end = $self->end();
		my $strand = $self->strand();
		$self->{position}->{$obj->id} = new Position({start=>$start,end=>$end,strand=>1});
		return $self->{position}->{$obj->id};
		
	}
	if ($obj->isTranscript()) {
		my $pos = $self->position($self->getChromosome());
		my $start = $obj->translate_position($pos->start);
		my $end = $obj->translate_position($pos->end);
		my $strand = 1;
		if ($start > $end){
			my $t = $start;
			$start = $end;
			$end = $t;
			$strand = -1;
		}
		$self->{position}->{$obj->id} = new Position({start=>$start,end=>$end,strand=>$strand});
		return $self->{position}->{$obj->id};	
	}
	elsif ($obj->isProtein()) {
		my $start = $self->getProteinPosition($obj);
		return new Position({start=>$start,end=>$start,strand=>1});	
	}
	else {
		confess();
	}
}

sub getOrfPosition {
	my ( $self, $prot ) = @_;
	confess ($self->id) unless $prot;
	confess($prot) if $prot->isTranscript();
	return -1 unless exists $self->annotation()->{ $prot->id }->{coding}->{sequences};
	return  $self->annotation()->{ $prot->id }->{coding}->{sequences}->{orf_position};
}

sub getProteinAA {
	my ( $self, $prot ) = @_;
	
	confess() if $prot->isTranscript();
	return unless exists $self->annotation()->{ $prot->id }->{coding}->{sequences};
	return  $self->annotation()->{ $prot->id }->{coding}->{sequences}->{aa};
}

###
# method isSomething
##

sub isCoding {
	my ( $self, $obj ) = @_;
	return $self->return_mask_testing($obj,"coding");	
}

sub isUpstream {
	my ( $self, $obj ) = @_;
	return $self->return_mask_testing($obj,"upstream");
}

sub isDownstream {
	my ( $self, $obj ) = @_;
	return $self->return_mask_testing($obj,"downstream");
}

sub isIntergenic {
	my ( $self, $obj ) = @_;
	return $self->return_mask_testing($obj,"intergenic");
}

sub isSilent {
	my ( $self, $obj ) = @_;
	return $self->return_mask_testing($obj,"silent");
	
}

sub isStop {
	my ( $self, $obj ) = @_;
	return $self->return_mask_testing($obj,"stop");
}

sub isPhase {
	my ( $self, $obj ) = @_;
	return $self->return_mask_testing($obj,"phase");
}

sub isEssentialSplicing {
	my ( $self, $obj ) = @_;
	return $self->return_mask_testing($obj,"essential_splicing");
}

sub isMatureMiRNA {
	my ( $self, $obj ) = @_;
	return $self->return_mask_testing($obj,"maturemirna");
}

sub isFrameshift {
	my ( $self, $obj ) = @_;
	return $self->return_mask_testing($obj,"frameshift");
}

sub isNonSynonymous {
	my ( $self, $obj ) = @_;
	return $self->return_mask_testing($obj,"nonsynonymous");
}

sub isUtr {
	my ( $self, $obj ) = @_;
	return $self->return_mask_testing($obj,"utr");
}

sub isConsequence {
	my ( $self, $obj,$cat) = @_;
	
	my $id = "all";
	if ($obj ne "all"){
		$id = $obj->id();
	}
	my $mask = $self->{annotation}->{$id}->{mask};
	confess() unless $mask;
	return ($mask & $self->project->getMaskCoding($cat));
}

############# 
# ANNEX
#############

sub annex {
	confess();
}

#has annex => (
#	is =>'rw',
#	#required => 1,
#);



sub isHighImpact {
	my ( $self, $obj ) = @_;
	return 1 if ($self->isStop($obj));
	return 1 if ($self->isPhase($obj));
	return 1 if ($self->isEssentialSplicing($obj));
	return 1 if ($self->isMatureMiRNA($obj));
	return 1 if ($self->isFrameshift($obj));
	return;
}

sub return_mask_testing {
	my ($self,$obj,$type) = @_;
	my $id = "all";
	if ($obj){
		#confess("call only on transcript") unless $obj->isTranscript();
		$id = $obj->id ;
	}
	
	unless  ($self->annotation()->{$id}->{mask}){
		
		warn $id." ".$self->getChromosome->name." ".$self->id;
		warn $self->annotation()->{$id}->{mask}." $type ".$self->getProject()->getMaskCoding($type);
		confess() ;
	}
#	warn $type." $id ".$self->id unless exists $self->annotation()->{$id}->{mask};
	return $self->annotation()->{$id}->{mask} & $self->getProject()->getMaskCoding($type);
}


=head2 changeAA
	Title   : changeAA
 	Usage   : $variation->changeAA();
 	Function: Get the changed AA cause by the variation
 	Returns : The changed AA (string)
 	Args	: A GenBoProtein to applicate the variation
	Note    : 
=cut

sub changeAA {
	my ( $self, $prot ) = @_;
	confess() if $prot->isTranscript();
	return "" unless exists $self->annotation->{ $prot->id }->{coding}->{sequences};
	return  $self->annotation->{ $prot->id }->{coding}->{sequences}->{aa_mut};
}


sub protein_nomenclature {
	my ( $self, $prot ) = @_;
		confess() unless $prot->isProtein();
	
		
		my $pos = $self->getProteinPosition($prot);
		return "-" unless $self->changeAA($prot);
			return "-" unless $self->getProteinAA($prot);
		#	warn $self->id."-".$self->getProteinAA($prot)."-".$pos."-".$self->changeAA($prot) unless $self->getProteinAA($prot);
		return $self->getProteinAA($prot).$pos.$self->changeAA($prot);
		
		#$hvariation->{codons_AA} =   $v->getProteinAA($prot).$hvariation->{prot}.$v->changeAA($prot);
		my $dec =  $-[0];
		$pos += $dec;
		
		return "Ins".substr($self->changeAA($prot),$dec).$pos;
		 
}

sub init_annotation {
	my $self = shift;
	my $debug  ;
#	$debug = 1 if $self->id eq "X_8095173_T_A";
#	warn "init annotation ".$self->id();
	#$self->set_origin_database() unless exists $self->{origindb}; ## a verifié 
	my $project = $self->getProject();
	my $annot;
	unless (@{ $self->getGenes()} ) {
#		warn "intergenic ".$self->getChromosome->name." :: ".$self->start;
		$self->{intergenic}++ unless @{ $self->getGenes() };
		$annot->{all}->{mask}  = $project->getMaskCoding("intergenic");
		#warn $self->id;
		$annot->{$self->id}->{mask}  = $project->getMaskCoding("intergenic");
		return $annot;
		
	}
	my $span = $self->getGenomicSpan();
	$annot->{all}->{mask} = 0;

	foreach my $tr ( @{ $self->getTranscripts() } ) {
		my $gid = $tr->getGene()->id();
		$annot->{$gid}->{mask} =  0 unless exists $annot->{$gid}->{mask};
		$annot->{$tr->id}->{mask} =  0;
		my $trid = $tr->id;
		
		###
		# test exonic
		###
		my $score_spliceAI = $self->max_spliceAI_score($tr->getGene());
		if ($score_spliceAI && $score_spliceAI ne "-" && $score_spliceAI >= $self->project->buffer->config->{spliceAI}->{medium}) {
							$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("predicted_splice_site");
							#last;
		}
		
		
		my $span2 = $span->intersection( $tr->getGenomicSpan );
		if ( $span2->is_empty ) {
			### intronic
			### test splice
			my $span_splice = $span->intersection( $tr->getSpanSpliceSite() );
			my $span_essential =   $span->intersection( $tr->getSpanEssentialSpliceSite() );
			
			$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("intronic");
			if (($tr->strand == 1 && $tr->start>$self->start) or ($tr->strand == -1 && $tr->end<$self->start)){
				$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("upstream");
			}
			if (($tr->strand == 1 && $tr->end<$self->start) or ($tr->strand == -1 && $tr->start>$self->start)){
				$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("downstream");
			}
			unless ($span_essential->is_empty){
				$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("essential_splicing");
			}
			else {
				if (not $span_splice->is_empty){
						$annot->{$trid}->{mask} = $annot->{$trid}->{mask} |$project->getMaskCoding("splice_site");
				}
				else {
					
				}
			}
			


		}
		else {
				
			### exonic
	
			
			if ($tr->protein() && $tr->protein !~ /ENST/ ) {
				#$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $maskCoding->{exonic};
				my $span_splice = $span->intersection( $tr->getSpanSpliceSite() );
				
				unless ( $span_splice->is_empty ) {
					$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("splice_site");
				}
			
				my $span_coding = $span->intersection( $tr->getSpanCoding );
				if ( $span_coding->is_empty ) {
					$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("utr");
				
				}
				else {

				
					#$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $maskCoding->{nonsynonymous};
					$self->annotation_coding($tr,$annot,$span_coding);

				}

			}
			else {
				#confess() if $self->id() eq "X_140008327_T_C";
				#non coding rna
				if ( $tr->isncRNA() ) {
					my $span_mature = $span->intersection($tr->span_mature());
					if ($span_mature->is_empty) {
						$annot->{$trid}->{mask} = $annot->{$trid}->{mask} |$project->getMaskCoding("ncrna");
					}
					else {
						$annot->{$trid}->{mask} = $annot->{$trid}->{mask} |$project->getMaskCoding("maturemirna");
					}
					
				}
				elsif ( $tr->ispseudoGene ) {
					$annot->{$trid}->{mask} = $annot->{$trid}->{mask} |$project->getMaskCoding("pseudogene");
					
				}
				else {
					$annot->{$trid}->{mask} = $annot->{$trid}->{mask} |$project->getMaskCoding("ncrna");	
				}
				
				#else {
				#	warn $tr->name();
				#	confess("problem");
				#}
			}
		}
		$annot->{$gid}->{mask} = $annot->{$gid}->{mask} | $annot->{$trid}->{mask};
		$annot->{all}->{mask} = $annot->{all}->{mask} | $annot->{$trid}->{mask};
	}

	
	foreach my $g  ( @{ $self->getGenes() } ) {
		next if exists  $annot->{$g->id}->{mask};
		$annot->{$g->id}->{mask} =  $project->getMaskCoding("intronic");
	}
	foreach my $g  ( @{ $self->getGenes() } ) {
		my $trgs = $g->getTranscripts();
		my %trs;
		map{$trs{$_->id}++} @{$self->getTranscripts()};
		foreach my $trgene (@$trgs){
			next if exists $annot->{$trgene->id}->{mask};
			$annot->{$trgene->id}->{mask} =   $project->getMaskCoding("intronic");
			$annot->{all}->{mask} = $annot->{all}->{mask} | $annot->{$trgene->id}->{mask};
		}
	}
	return $annot;

}







##### METHODS #####


has isCosmic => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->cosmic();
	},
);
  

has public_data_id => (
	is		=> 'ro',
	lazy=> 1,
	default=> sub {
		my $self = shift;
		my $pos = $self->start();
		my $seq = $self->sequence();
		my $id = join("_",($self->getChromosome()->name, $pos, "$seq"));
		return $id;
	},
);




 has public_db  =>(
	is		=> 'ro',
	lazy=> 1,
	default=> sub {
		my $self = shift;
		confess();
	},

);




has isDude => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
	return undef;
	},
);


has percent =>(
	is		=> 'ro',
	lazy=> 1,
	default=> sub {
		my $self = shift;
	
		return $self->frequency *100;
	},
);

sub isFoundBySVCaller {
	my ($self, $patient) = @_;
	confess() unless $patient;
	return $self->sequencing_infos->{$patient->id}->{max}->[4] == 1;
}

sub isDudeCaller {
	my ($self, $patient,$debug) = @_;
	return $self->sequencing_infos->{$patient->id}->{dud} if exists $self->sequencing_infos->{$patient->id}->{dud};
	return;
}


sub getNGSScore{
	my ($self) = @_;
	return $self->{ngs_score} if exists $self->{ngs_score};
	my $score =0;
	 my $max = -1;
	 my $filter = "None";

	foreach my $p (@{$self->getPatients()}){
		my $pid = $p->id;
		my $dp = $self->getDP($p);
		my $ratio = $self->getPourcentAllele($p);
		if ($ratio < 10){
			next;
		}
		if ($dp >= 30 && $ratio > 20 ){
			$score = 2;
			last;
		}
		elsif ($dp >= 10 ){
			$score = 1;
		}
	}
	$self->{ngs_score}  = $score;
	return $self->{ngs_score};
}

sub getSequence{
	my $self = shift;
	return $self->{var_allele};
}


sub effectImpact {
	my ( $self, $obj ) = @_;
	return $self->{impact}->{$obj->id} if exists $self->{impact}->{$obj->id};
	$self->annotation() unless exists $self->{annotation};
	my $id = "all";
	my $text = "unknown";
	if ($obj){
	if ($obj->isProtein){
			$id = $obj->getTranscript()->id 
		}
		else {
			$id= $obj->id;
		}
	}
	my $mask = $self->{annotation}->{$id}->{mask};
	my @levels = ("high","moderate","low");
	unless ($mask) {
	 foreach my $t (@{$self->getTranscripts}){
	 	warn $t->name;
	 }
	 warn $self->id()." id obj ".$id." ".$obj->id." ".$obj->isTranscript."-" unless $mask;
	 warn $obj->start." ".$obj->end;
	 confess();
	}
	
	foreach my $level (@levels){
		if ($mask &  $self->getProject()->getMaskCoding($level)) {
			$text = $level;
			last;
		}
	}
	 $self->{impact}->{$obj->id} = $text;
	return $text;
	
}



sub isEffectImpact {
	my ( $self, $obj,$level,$debug ) = @_;
	
	$self->annotation() unless exists $self->{annotation};
	my $id = "all";
	my $text = "unknown";
	if ($obj){
		if ($obj->isProtein){
			$id = $obj->getTranscript()->id 
		}
		else {
			$id= $obj->id;
		}
	}
	my $mask = $self->{annotation}->{$id}->{mask};
	return $mask &  $self->getProject()->getMaskCoding($level);#$self->getProject()->getMaskCoding($level);

	
}

sub variationTypeInterface {
	my ($self, $obj) = @_;
	my @lRes;
	foreach my $text (split(',', $self->variationType($obj))) {
		if (exists $self->project->buffer->config->{ensembl_annotations}->{$text}) {
			my @lTmp = split(';', $self->project->buffer->config->{ensembl_annotations}->{$text});
			push(@lRes, $lTmp[0]);
		}
		else { push(@lRes, $text); }
	}
	return join(',', @lRes);
}


sub get_codon_text {
	my ($self, $strand) = @_;
	return $self->{codon}->{$strand} if exists $self->{codon}->{$strand};
	if ($self->isDeletion){
		$self->{codon}->{$strand}  =  $self->delete_sequence."/".$self->sequence();
	}
	else {
		$self->{codon}->{$strand}  =  $self->getChromosome()->sequence($self->start,$self->end)."/".$self->sequence();
	}
	if ($strand == -1 ){
		$self->{codon}->{$strand}  =  BioTools::complement_sequence($self->getChromosome()->sequence($self->start,$self->end))."/".BioTools::complement_sequence($self->sequence());
	}
	return $self->{codon}->{$strand};
}

sub variationType {
	my ( $self, $obj ) = @_;
	$self->annotation() unless exists $self->{annotation};
	my $id = "all";
	my $text = "intergenic";
	if ($obj){
	if ($obj->isProtein){
			$id = $obj->getTranscript()->id 
		}
		else {
			$id= $obj->id;
		}
	}
	my $mask = $self->{annotation}->{$id}->{mask};
	unless($mask){
		warn 'ID wanted: '.$id;
		warn ref($obj);
		warn $obj->id;
		warn $self->id;
		warn $self->name();
		warn Dumper $self->init_annotation();
		warn Dumper $self->{annotation};
		confess();
	}
	$text = "";
	
	if ($mask &  $self->getProject()->getMaskCoding("exonic")) {
		if ( $mask & $self->getProject()->getMaskCoding("stop")) {
			$text = "stop";
			
		}
		elsif ( $mask & $self->getProject()->getMaskCoding("phase")){
			$text = "phase";
		}
		elsif ( $mask & $self->getProject()->getMaskCoding("essential_splicing")){
			$text = "essential_splicing";
		}
		elsif ( $mask & $self->getProject()->getMaskCoding("nonsynonymous")){
			$text = "coding";
		}
		elsif ( $mask &$self->getProject()->getMaskCoding("silent")){
		#	if ( $mask & $self->getProject()->getMaskCoding("splice_site") ){
		#		$text = "splicing";
		#	}
		#	else {
				$text = "silent";
		#	}
		}
		elsif ( $mask & $self->getProject()->getMaskCoding("frameshift")){
			$text = "frameshift";
		}
		elsif ( $mask & $self->getProject()->getMaskCoding("nonframeshift")){
			$text = "non-frameshift";
		}
		
		elsif (  $mask & $self->getProject()->getMaskCoding("utr") ) {
			$text = "utr";
		}
		else {
			confess();
		}
		if ( $mask & $self->getProject()->getMaskCoding("splice_site")){
			$text = "splicing,".$text;
		}
	} # end exonic
	elsif ( $mask & $self->getProject()->getMaskCoding("splice_site") ){
		$text = "splicing";
	}
	elsif ( $mask & $self->getProject()->getMaskCoding("essential_splicing")){
		$text = "essential_splicing";
	}
	
	elsif (  $mask &$self->getProject()->getMaskCoding("maturemirna") ) {

		$text = "maturemirna";
	}
	elsif (  $mask &$self->getProject()->getMaskCoding("ncrna") ) {

		$text = "ncrna";
	}
	elsif (  $mask & $self->getProject()->getMaskCoding("pseudogene") ) {

		$text = "pseudogene";
	}
	elsif ( $mask & $self->getProject()->getMaskCoding("essential_splicing")){
		$text = "essential_splicing";
	}
	elsif (  $mask &$self->getProject()->getMaskCoding("upstream") ) {
		$text = "upstream";
	}
	elsif (  $mask &$self->getProject()->getMaskCoding("downstream") ) {
		$text = "downstream";
	}
	elsif (  $mask &$self->getProject()->getMaskCoding("intronic") ) {
		$text = "intronic";
	}
	else {
		$text = "intergenic";
	}
	if ( $mask & $self->getProject()->getMaskCoding("predicted_splice_site") ){
		$text = $text.",predicted_splice_site";
	}
	#$text = "essential_splicing,".$text if $mask &  $maskCoding->{exonic} && $mask & $maskCoding->{essential_splicing};
	return $text;
}

=head2 polyphenStatus
	Title   : polyphenStatus
 	Usage   : $variation->polyphenStatus($protein);
 	Function: Get the polyphen status calculate by the polyphen programm
 	Returns : The polyphen status (int)
 	Args	: A GenBoProtein
	Note    : 
=cut

sub polyphenStatus {
	my ( $self, $obj ) = @_;
	my $mask = $self->getPredictionMask($obj);
	$mask = 0 unless $mask;
	return 3 if $mask & $self->getProject()->getMaskPrediction('polyphen_probably damaging');
	return 2 if $mask & $self->getProject()->getMaskPrediction('polyphen_possibly damaging');
	return 1 if $mask & $self->getProject()->getMaskPrediction( 'polyphen_benign');
	return 0 ;
	#
	# Non coding
	# Stop 
	# Silent
	# Possibly damaging
	# Probably damaging
	# Begnin 
	# Unknown
	#
	
	
}

sub getProteinPosition {
	my ( $self, $obj ) = @_;
	my $transcript = $obj;
	 if ($obj->isProtein()){
	 	# TODO: rustine pour correction de la prot NP001257972_19 et de son transcript NM001271043 (parfois ecris NM001271043_19_19 ...)
	 	if ($obj->{transcript} =~ /_/) {
	 		my ($tr_id, $chr_id) = split('_', $obj->{transcript});
	 		$obj->{transcript} = $tr_id;
	 		my $h;
	 		$h->{$tr_id.'_'.$chr_id} = undef;
	 		$obj->{transcripts_object} = $h;
	 	}
	 	$transcript = $obj->getTranscript($obj->{transcript});
	 	return -1 unless $transcript;
	 }
	  if ($obj->isTranscript()){
	 	#$transcript = $obj->getTranscript;
	 	return -1 unless $obj->getProtein;
	 }
	#confess() if $prot->isTranscript();
	$self->annotation() unless exists $self->{annotation};
	return unless exists $self->annotation->{ $transcript->id }->{coding}->{sequences};
	return  $self->annotation->{ $transcript->id }->{coding}->{sequences}->{prot_position};
}





sub rfPredScore { return undef; }

sub isTodo{
		my ($self,$patient) = @_;
		return unless exists $patient->{validation_ngs}->{$self->id} ;
		return $patient->{validation_ngs}->{$self->id} == -3;
}

sub is_ngs_rejected{
		my ($self,$patient) = @_;
		return unless exists $patient->{validation_ngs}->{$self->id} ;
		return $patient->{validation_ngs}->{$self->id} == -1;
}
sub is_ngs_validated{
		my ($self,$patient) = @_;
		return unless exists $patient->{validation_ngs}->{$self->id} ;
		return $patient->{validation_ngs}->{$self->id} > 0;
}
sub is_ngs_he{
		my ($self,$patient) = @_;
		return unless exists $patient->{validation_ngs}->{$self->id} ;
		return $patient->{validation_ngs}->{$self->id} == 2;
}
sub is_ngs_ho{
		my ($self,$patient) = @_;
		return unless exists $patient->{validation_ngs}->{$self->id} ;
		return $patient->{validation_ngs}->{$self->id} == 1;
}

sub is_sanger_rejected{
		my ($self,$patient) = @_;
		return unless exists $patient->{validation_sanger}->{$self->id} ;
		return $patient->{validation_sanger} ->{$self->id}== -5;
}
sub is_sanger_validated{
		my ($self,$patient) = @_;
			return unless exists $patient->{validation_sanger}->{$self->id} ;
		return $patient->{validation_sanger}->{$self->id} > 0;
}
sub is_sanger_he{
		my ($self,$patient) = @_;
			return unless exists $patient->{validation_sanger}->{$self->id} ;
		return $patient->{validation_sanger}->{$self->id} == 2;
}
sub is_sanger_ho{
		my ($self,$patient) = @_;
		return unless exists $patient->{validation_sanger}->{$self->id} ;
		return $patient->{validation_sanger}->{$self->id} == 3;
}
sub is_edited{
	my ($self,$patient) = @_;
	
	return undef unless (exists  $patient->{validation_sanger}->{$self->id} or exists  $patient->{validation_ngs}->{$self->id});
	return 1 if   $patient->{validation_sanger}->{$self->id} ne 0 or   $patient->{validation_ngs}->{$self->id} ne 0;
	return undef;
}



###########
# Various method for scoring
###########


sub score_variants_solo {
	my ($self,$child,$score,$tr,$debug) = @_;	
	my $pc = $self->getPourcentAllele($child);
	$pc = 0 if $pc eq "-";
	my $dp = $self->getDP($child);
	$dp = 0  if $dp eq "-";
	$score += 0.4 if $pc > 70;
	$score += 0.3 if  ($pc > 90);
	$score += 0.3 if $dp > 10;
		my $n = $self->getGnomadHO();
		$n = 0 unless $n;
		$n = 0 if $n eq "-";
		$score += 0.5 if $n == 0;
		$score -= 0.4 if $n > 3;
		return $score*2;
}

sub score_variants_trio_dominant {
	my ($self,$child,$score,$tr,$debug) = @_;
	my $fam = $child->getFamily();
	return ($score -3) unless $self->isDominantTransmission($fam,$child);
	return $score*2;
		
}


sub score_variants_trio {
	my ($self,$child,$score,$tr,$debug) = @_;	
	
	#return $self->{score}->{trio}->{$child->name} if exists $self->{score}->{trio}->{$child->name};
	my $fam = $child->getFamily();
	return $self->score_variants_trio_dominant($child,$score) if $fam->isDominant();
	#my $model = $self->getTransmissionModel($fam,$child);
	
	return ($score -2)  if $self->isBothTransmission($fam,$child,$debug);
	
	my $pc = $self->getPourcentAllele($child);
	$pc = 0 if $pc eq "-";
	my $dp = $self->getDP($child);
	$dp = 0  if $dp eq "-";
	
	if ($self->isMosaicTransmission($fam,$child)){
			#my $pcp = $self->getPourcentAllele($parent);
			warn "mosaic  " if $debug;
			$score += 0.3 if $pc > 40;
			return $score *2 ;
	}
	
	#TODO: ATTENTEION il peut avoir frere ou soeur -> le faire sur le vecteur FAMILLE et pas du PATIENT.
	if ($self->isRecessiveTransmission($fam)){
		warn "recessive  $score" if $debug;
		$score += 0.4 if $pc > 70;
		$score += 0.3 if  ($pc > 90);
		$score += 0.3 if $dp > 10;
		warn "recessive  $score" if $debug;
		my $n = $self->getGnomadHO();
		$n = 0 unless $n;
		$score += 0.5 if $n == 0;
		$score -= 0.4 if $n > 3;
		warn "===>Trio $score " if $debug;
		
		return $score*2;
	}
	
	
	my $isMotherTransmission = $self->isMotherTransmission($fam);
	my $isFatherTransmission = $self->isFatherTransmission($fam);

	
	#Uniparental disomy
	
	if ($self->isUniparentalDisomyTransmission($fam,$child)){
			warn "isUniparentalDisomyTransmission  " if $debug;
		#	warn $self->name();
		#$score +=20;
			$score += 0.3 if $pc > 80 ;
			return $score *=1.6;
	}
	
	if  ($isMotherTransmission or $isFatherTransmission) {
		
		warn "**** mother  " if $debug;
		#IDeal $self->getFamily->isFatherTransmission($fam,$child)){
	 	$score += 0.3 if $pc > 15;
	 	$score += 0.4 if $pc > 20;
		$score += 0.3 if $dp > 10;
		
		my $n = $self->getGnomadHO();
		$n = 0 unless $n;
		$n =0 if $n eq "-";
		$score += 0.25 if $n == 0;
		$score -= 0.4 if $n > 3;
	
		
		return $score;

	}
	
	if ($self->isMotherTransmission($fam,$child) == 1 ){
			$score -=2 unless $isMotherTransmission;
			return $score;
		}
		if ($self->isFatherTransmission($fam,$child) == 1){
			$score -=2 unless $isFatherTransmission;
			return $score;
		}
	#TODO: est ce que tous les freres et soeurs malades  malades ont aussi le meme modele recessif ?
	#TODO: oui -> petit plus au score
	#TODO: non -> le penaliser
	#if ($self->getModelVector_fam_denovo_in_all_children($fam) or $self->isStrictDenovoTransmission($fam)) {
	if ($self->isDenovoTransmission($fam)){# or $self->isStrictDenovoTransmission($fam)) {
			warn "denovo ".$score if $debug;
			
			my $n = $self->getGnomadAC();
			$n = 0 unless $n;
			$score += 0.5 if $n == 0;
			$score -= 0.4 if $n > 30;
			$score -= 0.75 if ($self->project->isGenome && ($self->isLargeDeletion or $self->isLargeDuplication)); 
			if ($self->isStrictDenovoTransmission($fam) && $fam->getFather && $fam->getMother){
					$score += 0.3 if   $pc > 20  ;
					$score += 0.4 if   $pc > 30  ;
					$score += 0.3 if $dp > 10;
				}
			elsif (!($self->isStrictDenovoTransmission($fam)) && $fam->getFather && $fam->getMother){
				$score -= 1;
			}
				else {
					
					$score += 0.1 if   $pc > 20  ;
					$score += 0.2 if   $pc > 30  ;
					#$score += 0.1 if $dp > 10;
					$score -= 0.5;# if   $pc > 20  ;
				}
				
				$score += 0.4 if ($pc > 20 && ($pc < 75 && $self->getChromosome->name ne "X"));
				$score += 0.3 if $dp > 10;
				
			return  $score*2 ;#if $fam->getFather && $fam->getMother;
	}
	
	
	return $score - 1;
			
}


sub score_prediction_refined {
	my ($self,$patient,$gene,$debug) = @_;
		my $score =0;
		my $cadd_score = $self->cadd_score();
#		
		$cadd_score = -1 unless $cadd_score;
		$cadd_score = -1 if $cadd_score eq "-";
		$cadd_score = -1 if $cadd_score == -1;
		my $alphamissense = $self->alphamissense($gene);
		
		if ($cadd_score > 50){
			$score =  0.5;
		}
		elsif ($cadd_score > 30 ){
			$score = 0.3
		}
		elsif ($cadd_score > 25){
			
			$score  = 0.2;
		}
		elsif ($cadd_score < 10 && $cadd_score > 0 ){
			$score  = -0.2;
		}
		
		my $scoreAI = $self->max_spliceAI_score($gene);
		$scoreAI = -1  if $scoreAI eq "-";
		if ($scoreAI>0.9 && !($self->isEssentialSplicing($gene)) ){
			$score += 0.3;
		}
		elsif ($scoreAI>0.7  && !($self->isEssentialSplicing($gene)) ){
			$score += 0.2;
		}
		elsif ($scoreAI>0.5  && !($self->isEssentialSplicing($gene))){
			$score += 0.1;
		}
		return $score;
}

sub score_refined {
	my ($self,$patient,$debug) = @_;
		$debug= undef;
		my $refined_score = 0;
		
		my $ac = $self->getGnomadHO;
		$ac = 0 unless $ac;
		$ac = 0 if $ac eq "-";
		my $p = $self->getPourcentAllele($patient);
	
		$p = 0 if $p eq "-";
		$p = 50 if $self->project->isSomatic;
		my$bad;
		if ($p < 10 ){
			$refined_score -= 1 ;
			#return $refined_score;
		}
#		warn "\t\t a $p ".$refined_score if $debug;
		if ($p < 15 ){
			$refined_score -= 1 ;
			#return $refined_score;
		}
#		warn "\t\t b ".$refined_score if $debug;
		if ($p < 20 ){
			$refined_score -= 0.5 ;
			#return $refined_score;
		}
		
		warn "\t\t c ".$refined_score if $debug;
		
		my $dd = $self->getDP($patient);
		$dd =0 if $dd eq "-";
			$dd =0 unless $dd ;
		if ($dd <5 ){
			$refined_score -=  1;
		}
		if ($dd <10 ){
			$refined_score -=  1;
		}
		if ($dd <15  ){
			$refined_score -= 0.5;
		}
		if ($dd <20 && $p< 30  ){
			$refined_score -= 1;
		}
		warn "\t\t".$refined_score if $debug;
		 if ($ac > 10 ){
							$refined_score -= 0.5 if $p > 20 ;
			}
		if ($ac > 20 ){
							$refined_score -= 0.2;# if $v->{max_pc} > 15;
	}	
		if ($ac > 50 ){
							$refined_score -= 0.5;# if $v->{max_pc} > 15;
	}
		if ($ac > 100 ){
							$refined_score -= 1;# if $v->{max_pc} > 15;
	}
	warn "\t\t depth  ". $refined_score." ==> $p ".$dd if $debug;
	my $pr = $self->other_projects;
	my $op = $self->other_patients_ho;
	warn $op if $debug;
	$pr = 0 unless $pr;
	$op = 0 unless $op;
		if ($pr >= 5 or $op >= 10){
						$refined_score  -= 0.3;# if $v->{max_pc} > 15;
		}
		if ($pr > 10 or $op > 20){
						$refined_score  -= 0.3;# if $v->{max_pc} > 15;
		}
		if ($pr > 70 or $op >50 ){
						$refined_score  -= 0.5;# if $v->{max_pc} > 15;
		}
		if ($pr > 200 or $op >= 100  ){
				$refined_score  -= 1 ;# if $v->{max_pc} > 15;
		}
		warn "\t\t".$refined_score if $debug;
		
		return $refined_score;
}
sub score_validations {
	my ($self,$gene) = @_;
	my $all_validations = $self->project->validations;
	my $zid = $gene->id."!".$self->id;
	return $all_validations->{$zid}->[0] if exists  $all_validations->{$zid};
	return undef;
}

sub scaledScoreVariant{
	my ($self, $tr,$patient,$vquery,$debug) = @_;
	$vquery = undef;
	#$debug = 1;
	# $debug = 1 if $self->id eq "7_158708074_CG_C";# && $patient->name eq "GEF2000057";
	#	if ($self->isClinvarPathogenic or $self->isDM && $self->scaled_score_frequence >= 2 ){
	#	return 5;
	#}
	#return $self->{scale_score}->{$tr->id}->{$patient->id} if exists  $self->{scale_score}->{$tr->id}->{$patient->id};
	my $gene;
	if ($tr->isGene){
		$gene = $tr;
	}
	elsif ($tr->isTranscript) {
		$gene = $tr->getGene;
	}
	
	my $score = $self->scoreVariant( $tr,$patient,$debug);
	my $scaled_score = 1;
	$scaled_score = 1;
	$scaled_score =  2 if $score > 50;
	$scaled_score = 3 if $score >= 150;
	$scaled_score = 4 if $score >= 200;
	warn "\t 1- ".$score." ".$scaled_score if $debug;
	
	$scaled_score += $self->score_refined($patient);
	warn "\t 2- ".$scaled_score." ".$self->score_refined($patient,$debug) if $debug;

	
	$scaled_score += $self->score_prediction_refined($patient,$gene);
	warn "\t 3- ".$scaled_score." ".$self->score_refined($patient,$debug) if $debug;
	
	if ($self->isSrPr){
	my $dp_before = $self->getNormDPBefore($patient);
	my $dp_after = $self->getNormDPAfter($patient);
	my $mean = ($dp_before+$dp_after)/2;
	$mean +=0.01 if $mean == 0;
	my $dp = $self->getNormDP($patient);
	$dp = 0,000001 if $dp == 0;
	my $ratio = $dp/$mean;
	my $logratio = $self->buffer->log2($ratio);
	$scaled_score -= 0.5 if ($mean < 10 );
	$dp += 0.0001;
	#warn $dp_after." ".$dp_before." ".$dp." ".$logratio;
	#$scaled_score -= 0.5 if $self->event_align_quality($patient) < 20;
	if ($self->isDeletion) {
		$scaled_score -= 0.5 if $logratio < 0.2 ;
		$scaled_score -= 0.5 if $logratio < 0.6 ;
	}
	elsif ($self->isLargeDuplication) {
		$scaled_score -= 0.5 if $logratio > -0.1 ;
		$scaled_score -= 0.5 if $logratio > -0.25 ;
	}
	else {
		 	my $sr1 = $self->sr1($patient);
			 $sr1 = 0 if $sr1 eq "-";
			my $pr1 = $self->pr1($patient);
			$pr1 = 0 if $pr1 eq "-";
			$scaled_score -= 0.5 if $sr1 < 5;
			$scaled_score -= 0.5 if $pr1 < 5;
	}
	
		
		
		
	}
	my $val = $self->score_validations($tr);
	
	if ($val){
			my $score = $val->{validation};
			$scaled_score = $val->{validation};
	}
	
	
	if ($patient->isChild && $patient->getFamily->isTrio()) {
		$scaled_score = $self->score_variants_trio($patient,$scaled_score,$tr,$debug);
		
	}
	else {
		$scaled_score = $self->score_variants_solo($patient,$scaled_score,$tr,$debug);
	}
	warn "\t trio- ".$scaled_score." :: ".$patient->name if $debug;
	if ($self->is_clinvar_pathogenic_for_gene($gene) or $self->isDM_for_gene($gene)  ){
		my $gac  = $self->getGnomadAC ;
		$gac = 0 unless $gac;
		$scaled_score ++;
		if ($self->project->isHGMD){
		$scaled_score += 0.5 	if ($self->is_clinvar_pathogenic_for_gene($gene) && $self->isDM_for_gene($gene));
		$scaled_score += 0.2 	if ($self->is_clinvar_pathogenic_for_gene($gene) && $self->isDM_for_gene($gene) && $gac < 1000);
		}
		else {
				$scaled_score += 0.5 	if ($self->is_clinvar_pathogenic_for_gene($gene) );
				$scaled_score += 0.2 	if ($self->is_clinvar_pathogenic_for_gene($gene) && $gac < 1000);
		}
		$scaled_score +=0.5 	if  $gac < 100 ;
		$scaled_score ++ 	if  $gac < 30;
	}
	
	warn "\t 4- ".$scaled_score if $debug;
#	die() if $debug;
	$self->{scale_score}->{$tr->id}->{$patient->id} = $scaled_score;
	return $self->{scale_score}->{$tr->id}->{$patient->id};
}

sub scaledScoreVariantPolydiag{
	my ($self, $tr,$patient,$query,$debug) = @_;
	$debug = undef;
	 #$debug = 1 if $self->id eq "1_209961988_A_C";# && $patient->name eq "GEF2000057";
	#	if ($self->isClinvarPathogenic or $self->isDM && $self->scaled_score_frequence >= 2 ){
	#	return 5;
	#}
	#return $self->{scale_score}->{$tr->id}->{$patient->id} if exists  $self->{scale_score}->{$tr->id}->{$patient->id};
	

	my $score = $self->scoreVariant ( $tr,$patient,$debug);
	my $scaled_score = 1;
	$scaled_score = 1;
	$scaled_score =  2 if $score > 50;
	$scaled_score = 3 if $score >= 150;
	$scaled_score = 4 if $score >= 200;
	warn "\t 1- ".$score." ".$scaled_score if $debug;
	$scaled_score += $self->score_refined($patient);
	warn "\t 2- ".$scaled_score." ".$self->score_refined($patient,$debug) if $debug;
	my $val = $self->score_validations($tr);

	if ($val){
			my $score = $val->{validation};
			$scaled_score = $val->{validation};
	}
	
	if ($patient->isChild && $patient->getFamily->isTrio()){
		my $sc = $scaled_score;
		$sc = $self->score_variants_trio($patient,$scaled_score,$tr,$debug);
	
		if ($sc<=$scaled_score ){
			$scaled_score -=1;
		}
	}
	
	warn "\t trio- ".$scaled_score." :: ".$patient->name if $debug;
	if ($self->is_clinvar_pathogenic_for_gene($tr->getGene()) or $self->isDM_for_gene($tr->getGene())  ){
		my $gac  = $self->getGnomadAC ;
		$gac = 0 unless $gac;
		$scaled_score ++;
		my $gene = $tr->getGene();
		if ($self->project->isHGMD){
		$scaled_score += 0.5 	if ($self->is_clinvar_pathogenic_for_gene($gene) && $self->isDM_for_gene($gene));
		$scaled_score += 0.2 	if ($self->is_clinvar_pathogenic_for_gene($gene) && $self->isDM_for_gene($gene) && $gac < 1000);
		}
		else {
				$scaled_score += 0.5 	if ($self->is_clinvar_pathogenic_for_gene($gene) );
				$scaled_score += 0.2 	if ($self->is_clinvar_pathogenic_for_gene($gene) && $gac < 1000);
		}
		$scaled_score +=0.5 	if  $gac < 500 ;
		$scaled_score ++ 	if  $gac < 80;
		
	}
		
	warn "\t 3- ".$scaled_score if $debug;
#	die() if $debug;
	$self->{scale_score}->{$tr->id}->{$patient->id} = $scaled_score;
	
	return $self->{scale_score}->{$tr->id}->{$patient->id};
}


has score_frequence_public => (
	is		=> 'rw',
	#required=> 1,
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $score =0;
		unless ($self->frequency()) {
				
					$score  += 100;
					return $score;
		}
		if ($self->frequency() < 0.0001){
					$score  += 100;
		}	
		elsif($self->frequency() < 0.001){
					$score  += 75;
		}
		elsif($self->frequency() < 0.01){
					$score  += 50;
		}
		elsif($self->frequency() < 0.05){
					$score  -=150;
		}
		else{
					$score  -=1000;
		}
		return $score;
	}
);

has scaled_score_cadd => (
	is		=> 'rw',
	#required=> 1,
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $score =0;
		my $cadd_score = $self->cadd_score();
		return 0 unless $cadd_score;
		return 0 if $cadd_score eq "-";
		return 0 if $cadd_score  == -1;
		if ($cadd_score > 50){
			$score  += 100;
		}
		if ($cadd_score > 30 ){
			$score  += 60;
		}
			if ($cadd_score > 25){
			$score  += 20;
		}
		if ($cadd_score < 10 ){
			$score -=20;
		}
		return $score;
	}
);

has scaled_score_frequence_public =>(
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $freq = $self->frequency();
		my $freq_score = 1;
		my $freq_level = 4;
		unless ($freq) {
			$freq =0;
		}
		if ($freq < 0 ){
							$freq_score = 4; 
							$freq_level = 1; 
		}
		unless ($freq) {
				$freq_score = 4; 
				$freq_level = 1; 
		}
		elsif ($freq <= 0.001){
							$freq_score = 4; 
							$freq_level = 1; 
		}
		
		elsif ($freq <= 0.01){
				$freq_score = 2; 
				$freq_level= 3; 
		}
		
		return {freq_score=>$freq_score,freq_level=>$freq_level};
	}
);


has scaled_score_frequence_dejavu => (
	is		=> 'rw',
	#required=> 1,
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
			my $hvariation;
				my $freq_score = 1;
				my $freq_level = 4;
				my $nb = $self->other_projects;
				$nb = 0 unless $nb;
				if ($nb == 0){
					$freq_score = 4; 
					$freq_level =1;
					
				}
				elsif ($nb <= $self->project->buffer->{config}->{dejavu}->{rare}  ){
							$freq_score = 3; 
							$freq_level =2;
					
				}
				elsif ($nb <=  $self->project->buffer->{config}->{dejavu}->{occasional}){
						$freq_score = 2; 
						$freq_level =3;
					
				}
				else {
						$freq_score = 1;
						$freq_level =4; 
				}
					$hvariation->{freq_score} = $freq_score;
					$hvariation->{freq_level} = $freq_level;
					return {freq_score=>$freq_score,freq_level=>$freq_level};
		}	
);

has scaled_score_frequence =>(
	is		=> 'rw',
	#required=> 1,
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
			my $hvariation;
			my $scaled_score =  $self->scaled_score_frequence_public();
			my $score_dj = $self->scaled_score_frequence_dejavu;
			if ($scaled_score->{freq_score} == 4 or $scaled_score->{freq_score} == 2) {
						$scaled_score = $self->scaled_score_frequence_dejavu;
					
				}
			return $scaled_score;
		}	
);


has score_frequence_dejavu => (
	is		=> 'rw',
	#required=> 1,
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
	my $score =0;
	my $nbPat_sim = ($self->similar_patients/$self->total_similar_patients)*100;
	#warn $nbPat_sim;
	#warn $nbPat_sim;
				if ($nbPat_sim > 60){
					$score -= 100;
				}elsif ( $nbPat_sim> 50){
					$score -= 70;
				}elsif ($nbPat_sim > 40){
					$score -= 40;
				}
			
			
		return $score;
	}
);

has score_validation_public => (
	is		=> 'rw',
	#required=> 1,
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $score =0;
		return 		$score  unless $self->score_clinvar();
		if ($self->score_clinvar()  >= 5)  {
							$score += 50;
		}
		elsif ($self->score_clinvar()  == 2)  {
				$score -= 100;
		}
		return $score;
	}
);






sub score_patient_quality {
		my ($self,$patient) = @_;
		my $score;
		if ($self->isHomozygote($patient)){
					$score += 50;
		}
		else
		{
					$score += 30;
		}
		my $nb_all_ref = $self->getNbAlleleRef($patient);#	unless ($nb_all_ref){$nb_all_ref = 0;}
		my $nb_all_mut = $self->getNbAlleleAlt($patient);#	unless ($nb_all_mut){$nb_all_mut = 0;}
		$nb_all_ref = 0 if $nb_all_ref eq "-";
		$nb_all_mut = 0 if $nb_all_mut eq "-";
		my $nb_read = $nb_all_ref + $nb_all_mut;

## fonction du nbre de read
		if($nb_read < 10){
					$score -=100;
				}
				elsif($nb_read < 30){
					$score -=50;
			}
			return $score;
}

sub score_transcript_consequence {
	my ($self,$tr,$debug) = @_;
	my $score = 0;
	$score  -= 100  unless $tr->isMain();
#	warn $tr->name() if $debug;
		
	if($self->isEffectImpact($tr,'high')) {
					$score  += 150;
					
		}
	elsif($self->isEffectImpact($tr,'moderate')){
				$score  += 50;
	}
	warn "# ==> ".$score if $debug;
	my $geneCsq = $self->variationType($tr);
	warn $geneCsq if $debug;
	if ($geneCsq =~/splic/ ){
		my $ename = lc($self->nearest_exons($tr));
		$score -= 150 if $ename =~ /nc/
	}
	if($geneCsq eq 'non-frameshift' or $geneCsq eq 'silent'){
					$score  += 10;
	}
		warn "# ==> ".$score if $debug;
	if ($geneCsq =~/splic/ && $geneCsq !~ /cod/){
			if ($self->dbscsnv_ada ne "-" or  $self->dbscsnv_rf ne "-"){
				warn "# ".$score if $debug;
			my $a = 0;
			warn "# ".$score if $debug;
			$a = $self->dbscsnv_ada if 	$self->dbscsnv_ada ne "-";
			$a = ($a > $self->dbscsnv_rf?$a: $self->dbscsnv_rf) if 	$self->dbscsnv_rf ne "-";
			warn "# ".$a if $debug;
			$score += 50 if $a>0.95;
			$score += 25 if $a>0.75;
			$score -= 100 if $a< 0.5;
			}
			else {
				$score -=20;
			}
			my $h_score_spliceAI = $self->spliceAI_score($tr->getGene());
			
			
	}
	$score += 50 if  $self->max_spliceAI_score($tr->getGene()) > 0.9;
	$score += 50 if  $self->max_spliceAI_score($tr->getGene()) > 0.7;
	#$score += 100 if  $self->max_spliceAI_score($tr->getGene()) > 0.9;
#	die($geneCsq." ".$self->id) if $self->max_spliceAI_score($tr->getGene()) > 0.7;
	
	
	#else {
	#	$score  -= 100;
	#}
	
	warn "# return ".$score if $debug;
	return $score;	
} 





sub scoreVariant{
	my ($self, $obj,$patient,$debug) = @_;
	
	$debug= undef;
	#return  $self->{score_variant}->{$tr->id}->{$patient->id}  if exists  $self->{score_variant}->{$tr->id}->{$patient->id} ;
	my $score =-150;
	if ($obj->isGene){
		foreach my $tr (@{$self->getTranscripts}){
			next if $tr->getGene->id ne $obj->id;
			my $tscore = $self->score_transcript_consequence($tr,$debug);
			
			#warn "transcript  score : ".$tscore if $debug;
			$score = ($score>$tscore?$score:$tscore);
			last if $score == 150;
		}
		
	}
	elsif ($obj->isTranscript){
			$score = $self->score_transcript_consequence($obj,$debug);
	}
	else {
		confess();
	}
	warn "\t\t\t score transcript ". $score if $debug;
#	die() if $debug;
	#my $score += $self->score_transcript_consequence($tr);
	$self->isEffectImpact($obj,'high',$debug);
	$score += $self->score_patient_quality($patient);
	#warn "\t\t\t uality ". $score if $debug;
	$score +=  $self->score_frequence_public();
	#warn "\t\t\t frequence ". $score if $debug;
	#warn $score if $debug;
	$score += $self->scaled_score_cadd();
	#warn "\t\t\t cadd ". $score if $debug;
	$self->{score_variant}->{$obj->id}->{$patient->id} = $score;
	#warn $score if $debug;
	return $score;
}
sub getNomenclatureForUtr {
	my ( $self, $transcript) = @_;
	my $pos_transcript = $transcript->translate_position($self->start);
		 my $seqv = $self->sequence();
		 my $seqr = $self->getChromosome()->sequence($self->start,$self->start);
		 die() unless  $seqr;
		if ($transcript->strand() == -1 ){
				 $seqv = BioTools::complement_sequence($seqv);
				 $seqr = BioTools::complement_sequence($seqr);
			}
			if ($pos_transcript   < $transcript->orf_start()){
			return "c.".($pos_transcript-$transcript->orf_start())."".$seqr.">".$seqv if ($self->isVariation);
			return "c.".($pos_transcript-$transcript->orf_start())."ins".$seqv  if ($self->isInsertion);
			return "c.".($pos_transcript-$transcript->orf_start())."+".length($seqr)."del"  if ($self->isDeletion);
		}
		else {
			return "c.*".($pos_transcript-$transcript->orf_end())."".$seqr.">".$seqv if ($self->isVariation);
			return "c.*".($pos_transcript-$transcript->orf_end())."ins".$seqv if ($self->isInsertion);
			return "c.*".($pos_transcript-$transcript->orf_end())."+".length($seqr)."del"  if ($self->isDeletion);
		}
}

sub getPositionForNomenclatureIntronic {
	my ( $self, $transcript, $debug ) = @_;
# >>>
	if ($self->start  < $transcript->start){
		return (-1,"") unless 		$transcript->genomic_orf_start();
		my $p1 = ($self->start - $transcript->genomic_orf_start());
		$p1 = $p1 * $transcript->strand();
		return($p1,"")
		
	}
	if ($self->end  > $transcript->end){
		return (-1,"") unless 		$transcript->genomic_orf_end();
		my $p1 = ($self->end - $transcript->genomic_orf_end());
		$p1 = $p1 * $transcript->strand();
		return($p1,"")
		
	}
	my ($dist,$nearest,$ref_pos) = $transcript->computeNearestExon($self->start,$self->end);
	
	my $pos = $transcript->translate_coding_position($ref_pos);
	if ($pos == -1 && $transcript->genomic_orf_start() ){
		
		if ($ref_pos   < $transcript->genomic_orf_start()){
			$pos = $ref_pos-$transcript->genomic_orf_start();
		}
		else {
			$pos = $ref_pos-$transcript->genomic_orf_end()
		} 
	}

#	warn $transcript->strand;
		 return ($dist,$pos);
}

sub getNomenclature {
		my ( $self, $transcript, $debug ) = @_;
	return $self->{nomenclature}->{$transcript->id} if exists $self->{nomenclature}->{$transcript->id};
	$self->{nomenclature}->{$transcript->id} = $self->constructNomenclature($transcript, $debug );
	return $self->{nomenclature}->{$transcript->id};
}


sub getNeededCoord {
        my $self = shift;
        my $start = $self->start();
        my $end = $self->end();
        if ($self->isInsertion() or $self->isDeletion()) {
                my $max = length($self->ref_allele());
                $max = length($self->var_allele()) if (length($self->var_allele()) > $max);
                $start = $start - $max - 20;
                $end = $end + $max + 20;
        }
        return ($start, $end);
}

sub getReference {
	my $self = shift;
	return $self->getChromosome->getReferences($self->start()-1, $self->end()+1)->[0];
}

sub nearest_exons{
	my ($self,$tr1) = @_;
#	confess();
	confess() unless $tr1;
	return $self->{exons}->{$tr1->id} if exists $self->{exons}->{$tr1->id};
		 $self->{exons}->{$tr1->id} = $tr1->findExonNumber($self->start,$self->end);
		 $self->{exons}->{$tr1->id} = $tr1->findNearestExon($self->start,$self->end) if  $self->{exons}->{$tr1->id} == -1;
		return  $self->{exons}->{$tr1->id};
		
}
sub isRecessive {
	my ($self,$fam) = @_;
		confess("code it by yourself");
		
}
sub isFatherTransmission {
	my ($self,$fam) = @_;
	confess("code it by yourself")
}

sub isMotherTransmission {
	my ($self,$fam) = @_;
	confess("code it by yourself")
}

sub isBothTransmission{
	my ($self,$fam) = @_;
		confess("code it by yourself ".$self);
}
sub isDenovo {
	my ($self,$fam) = @_;
		confess("code it by yourself")
}
sub isStrictDenovo{
	my ($self,$fam) = @_;
		confess("code it by yourself")
}

sub getGeneticModel {
		my ($self,$fam) = @_;
		my $vector = 
		confess("code it by yourself")
}

sub returnTextGenotype {
	 my ($self,$value) = @_;
	 return "he" if $value == 1;
	 return "ho" if $value == 2;
	 return "ref" if $value == 0;
	  confess();
}


sub isHomozygote {
	my ($self, $patObj) = @_;
	
	return $self->getSequencingGenotype($patObj) eq "ho";
}
sub isHeterozygote {
	my ($self, $patObj) = @_;
	return $self->getSequencingGenotype($patObj) eq "he";
}

sub text_heho {
 my ($self,$patient) = @_;
 return unless $self->existsPatient($patient);
 return $self->getSequencingGenotype($patient).":".$self->getNbAlleleRef($patient).":".$self->getNbAlleleAlt($patient);
}

	#$seq_infos->{$method_name}->{nb_ref} = $nb_ref;
	#	$seq_infos->{$method_name}->{nb_alt} = $nb_alt;



###LMDB_DEPTH_INFOS
has dp_infos =>(
	is =>'ro',	
	lazy=>1,
	default=> sub {
		my $self = shift;
		 my $hash = {}; 
		 return {} unless $self->isSrPr();
		foreach my $pat (@{$self->getPatients}){
			foreach my $patient (@{$pat->getFamily()->getMembers}){
					my $pid = $patient->id;
				unless ($patient->hasBamFile) {
					$hash->{$pid} = ["-","-","-"];
					next;
				}

			#next unless $patient->isGenome;
				my $mean_dp =  int($patient->meanDepth($self->getChromosome->name, $self->start, $self->end));
				my $norm_depth =  int($patient->mean_normalize_depth($self->getChromosome->name, $self->start, $self->end+1));
				my $norm_depth_before = int($patient->mean_normalize_depth($self->getChromosome->name, $self->start-510, $self->start-10));
				my $norm_depth_after = int($patient->mean_normalize_depth($self->getChromosome->name, $self->end+10, $self->end+510));
				
				my $dude = "-";
				if($self->project->isGenome){
					$dude = "-";
					#$dude = int($patient->cnv_value_dude($self->getChromosome->name,$self->start,$self->start+$self->length)*100);
				}
				$hash->{$pid} = [$mean_dp,$norm_depth,$dude,$norm_depth_before,$norm_depth_after];
			}
						
		}
			return $hash;
	},
);



sub getMeanDP {
		my ($self,$patient,$method) = @_;
		my $pid = $patient->id;
		if ($self->isSrPr){
			warn "cuicui";
			return $self->split_read_infos->{$pid}->[6];
		}
		else {
			$self->dp_infos->{$pid}->[0];
		}
}

sub getNormDP {
		my ($self,$patient,$method) = @_;
		my $pid = $patient->id;
		return $self->dp_infos->{$pid}->[1];
}
sub getNormDPBefore {
		my ($self,$patient,$method) = @_;
		my $pid = $patient->id;
		return $self->dp_infos->{$pid}->[3];
}
sub getNormDPAfter {
		my ($self,$patient,$method) = @_;
		my $pid = $patient->id;
		return $self->dp_infos->{$pid}->[4];
}

sub getLog2Ratio {
		my ($self,$patient,$method) = @_;
		my $m = int($self->getNormDPBefore($patient,$method)+ $self->getNormDPAfter($patient,$method))/2;
		return "-1" if $m ==0;
		return  $self->buffer->log2($self->getNormDP($patient,$method)/$m);
}


sub getCNVDude {
	my ($self,$patient) = @_;
	my $pid = $patient->id;
	return "??" unless exists  $self->dp_infos->{$pid};
	return "??" if scalar @{$self->dp_infos->{$pid}} == 0;
	return "??" if $self->dp_infos->{$pid}->[2] eq "-";
	warn  $self->id;
	confess();
	return  $self->buffer->log2($self->dp_infos->{$pid}->[2]/100);
}







sub existsPatient {
	my ($self,$patient) = @_;
	$self->getPatients() unless exists  $self->{patients_object};
	return exists $self->{patients_object}->{$patient->id};
}

has test_patients =>(
	is =>'ro',	
	lazy=>1,
	default=> sub {my $self = shift; my %toto; 
			foreach my $p (@{$self->getPatients}){$toto{$p->name}++}
			return \%toto;
	},
);

sub getGenotype {
	my ($self,$patient) = @_;
	unless (exists $self->test_patients->{$patient->name}){
		return $self->{'ref_allele'}." ".$self->{'ref_allele'};
	}
	elsif ($self->isHomozygote($patient)){
		return $self->sequence." ".$self->sequence;
	}
	elsif ($self->isHeterozygote($patient)){
		return $self->{'ref_allele'}." ".$self->sequence;
	}
	confess();
}

########################
#SEQUENCING INFOS
########################

#TODO: cas 1/2
sub check_if_1_2_cas {
	my ($self, $patient) = @_;
	return  $self->sequencing_infos->{$patient->id}->{max}->[3];
}




sub getDepth {
	my ($self,$patient,$method) = @_;
	return $self->getDP($patient,$method);
}

sub getDP {
		my ($self,$patient,$method) = @_;
		
		my $pid = $patient->id;
		unless ($self->existsPatient($patient)){
			return int($patient->meanDepth($self->getChromosome->name, $self->start, $self->end+1) ) ;
		}
		confess() unless exists $self->sequencing_infos->{$pid};
		my $res;
		unless ($method){
			$method = "max";
		}
		elsif ($self->project->return_calling_methods_short_name($method)) {
			$method =$self->project->return_calling_methods_short_name($method);
		}
		$res = int($self->sequencing_infos->{$pid}->{$method}->[0] + $self->sequencing_infos->{$pid}->{$method}->[1]);
		return $res;
}


sub getNbAlleleAlt {
		my ($self,$patient,$method) = @_;
		my $pid = $patient->id;
			unless ($self->existsPatient($patient)){
			return "-";
		}
		unless ($method) {
			$method = "max";
		}
		elsif ($self->project->return_calling_methods_short_name($method)) {
			$method =$self->project->return_calling_methods_short_name($method);
		}
		my $res = $self->sequencing_infos->{$pid}->{$method}->[1];
		return $res;
}

sub getNbAlleleRef {
		my ($self,$patient,$method) = @_;
	my $pid = $patient->id;
			unless ($self->existsPatient($patient)){
			return "-";
		}
		unless ($method) {
			$method = "max";
		}
		elsif ($self->project->return_calling_methods_short_name($method)) {
			$method = $self->project->return_calling_methods_short_name($method);
		}
#		warn "\n";
#		warn Dumper  $self->sequencing_infos;
#		warn 'PID: '.$pid;
#		warn $method;
		my $res = $self->sequencing_infos->{$pid}->{$method}->[0];
		return $res;
}


sub getSequencingGenotype {
	 my ($self,$patient,$method) = @_;
        my $pid = $patient->id;
       return "" unless exists $self->{patients_object}->{$patient->id};
       	unless ($method) {
			$method = "max";
		}
		elsif ($self->project->return_calling_methods_short_name($method)) {
			$method =$self->project->return_calling_methods_short_name($method);
		}
        my $res = $self->sequencing_infos->{$pid}->{$method}->[2];
			return $res;
}


#sub getSequencingInfos {
#sub getTextSequencingInfos	{
#        my ($self,$patient) = @_;
#        my $pid = $patient->id;
#        my @string;
#        foreach my $m (keys %{$self->sequencing_infos->{$pid}}) {
#        	next if $m eq "max";
#        	my $v = $self->sequencing_infos->{$pid}->{$m};
#        	push(@string,$v->[0].":".$v->[1]."(".$v->[2]."/".$v->[3].")");
#        }
#        return \@string;
#}
#
#sub getTextSequencingMethods {
#	 my ($self,$patient) = @_;
#      my $pid = $patient->id;
#      my @string;
#      foreach my $v (@{$self->sequencing_infos->{$pid}->{values}}){
#        	push(@string,$v->[0]);
#        }
#        return \@string;
#        
#}
#sub getTextSequencingRatio {
#        my ($self,$patient) = @_;
#        my $pid = $patient->id;
#
#        my @string;
#        foreach my $v (@{$self->sequencing_infos->{$pid}->{values}}){
#        	my $pc = "-";
#        
#        	 $pc = sprintf("%.0f", ($v->[3]/($v->[2]+$v->[3]))*100 ) if ($v->[2]+$v->[3]) > 0;
#        	push(@string,$v->[0].":$pc%");
#        }
#        return join(";",@string);
#}



sub getPourcentAllele {
        my ($self,$patient,$method) = @_;
        my $pid = $patient->id;
        
		unless (exists $self->sequencing_infos->{$pid}){
			return 0;
		}
	 	unless ($method) {
			$method = "max";
		}
		elsif ($self->project->return_calling_methods_short_name($method)) {
			$method =$self->project->return_calling_methods_short_name($method);
		}
		confess() unless exists $self->sequencing_infos->{$pid};
		my $alt =$self->sequencing_infos->{$pid}->{$method}->[1];
		$alt = 0 unless $alt;
		my $sum = (int($self->sequencing_infos->{$pid}->{$method}->[0])+int($alt));
		return 100 if $sum == 0;
		my $v = $self->sequencing_infos->{$pid}->{$method}->[1];
		$v = 0.00001 unless $v; 
		return (($self->sequencing_infos->{$pid}->{$method}->[1]/$sum)*100);
}

sub getRatio {
        my ($self,$patient,$method) = @_;
		return sprintf("%.0f", ($self->getPourcentAllele($patient,$method)));
}

sub array_sequencing_text {
	  my ($self,$patient) = @_;
	   my $pid = $patient->id;
	   return [] unless exists $self->sequencing_infos->{$pid};
	   my @res;
	   foreach my $m (keys %{$self->sequencing_infos->{$pid}}){
	   	next if $m eq "max";
	   	my $a = $self->sequencing_infos->{$pid}->{$m};
	   		push(@res,$m.":".$a->[2]."(".$a->[0]."/".$a->[1].")");
	   }
	 return \@res;
}


sub getMethods {
	my ($self,$patient) = @_;
	my @ms;
	delete $self->{sequencing_infos} if exists $self->sequencing_infos->{$patient->id}->{values};
	   foreach my $m (keys %{$self->sequencing_infos->{$patient->id}}){
	   	 	next if $m eq "max";
	   		push(@ms,$m);
	   }	   	
	   return \@ms;
	   
}


has sequencing_infos =>(
	is =>'ro',	
	lazy=>1,
	default=> sub {
			my $self = shift;
			 my $hash; 
			 $hash = {};
			foreach my $patient (@{$self->getPatients}){
				my $pid = $patient->id;
				# $hash->{$pid}->{ok} =1;
				 my @methods = sort {$a cmp $b} keys %{ $self->{annex}->{$patient->id}->{method_calling} };
				 my $mr = -1;
				 my $ma = -1;
				 my $genotype;
				 $genotype->{ho} = 0;
				 $genotype->{he} = 0;
				 my $sv = 0;
				foreach my $method (@methods){
					#next unless exists $self->sequencing_details()->{$patient->id}->{method_calling}->{$method};
					next unless exists $self->{annex}->{$patient->id}->{method_calling}->{$method}->{nb_all_ref};
					my $all_annex = $self->{annex}->{$patient->id}->{method_calling}->{$method};
					
					my $nb_ref =$all_annex->{nb_all_ref};
					$nb_ref = 0 if $all_annex->{nb_all_ref} and $all_annex->{nb_all_ref} eq "?";
					my $nb_alt =  $all_annex->{nb_all_mut};
					$nb_alt = 0 if $all_annex->{nb_all_mut} and $all_annex->{nb_all_mut} eq "?";
					my $nb_all_other_mut = 0;
					$nb_ref += $all_annex->{nb_all_other_mut} if (exists $all_annex->{nb_all_other_mut});
					my $method_name = $self->project->return_calling_methods_short_name($method);
					if ($self->project->is_sv_calling_methods($method)){
						$sv = 1;
					}
					my $type = "he";
					$type = "ho" if $all_annex->{ho};
					$genotype->{$type} ++;
					$hash->{$pid}->{$method_name} = [$nb_ref,$nb_alt,$type];
					$mr = $nb_ref if $mr < $nb_ref;
					$ma = $nb_alt if $ma < $nb_alt;
				}
				my $g = "he";#
				$g = "ho" if ($genotype->{ho} > $genotype->{he});
				my $case = 0;
				$case = 1  if exists $self->{annex}->{$patient->id}->{is_cas_1_2};
				$hash->{$pid}->{max} = [$mr,$ma,$g,$case,$sv];
				
			}
			return $hash;
	},
);

sub sequencing_details {
	my ($self,$patient) = @_;
	#confess();
	my $pid = $patient->id;
	return $self->{seq_details}->{$pid} if exists $self->{seq_details}->{$pid};
	#return $self->{seq_infos}->{$pid} if exists $self->{seq}->{$pid};
	my @methods;
	my @p_ratio;
	my $seq_infos;
	foreach my $method (@{$patient->callingMethods}){
			next unless exists $self->annex()->{$patient->id}->{method_calling}->{$method};
			next unless exists $self->annex()->{$patient->id}->{method_calling}->{$method}->{nb_all_ref};
			$self->annex()->{$patient->id}->{method_calling}->{$method}->{nb_all_mut} = "?" unless defined $self->annex()->{$patient->id}->{method_calling}->{$method}->{nb_all_mut};
			my $all_annex = $self->annex()->{$patient->id}->{method_calling}->{$method};
			my $nb_all_other_mut = 0;
			$nb_all_other_mut = $all_annex->{nb_all_other_mut} if (exists $all_annex->{nb_all_other_mut});
			
			my $nb_ref =$all_annex->{nb_all_ref};
			$nb_ref = 0 if $all_annex->{nb_all_ref} eq "?";
		
			my $nb_alt =  $all_annex->{nb_all_mut};
			$nb_alt = 0 if $all_annex->{nb_all_mut} eq "?";
			my $type;
			my $method_name = substr $method,0,3;
			if ($method eq "dragen-sv"){
				$method_name = "dsv";
			}
			elsif ($method eq "dragen-cnv"){
				$method_name = "dcnv";
			}
			
			push(@methods,$method_name);
			my $sequence_info = "he("; 
			$type = "he";
			my $pc ="-";		
			if ($self->annex()->{$patient->id}->{nb_all_ref} eq "?" ){
				$sequence_info = "??";
			}
			elsif ($nb_ref eq "?" ){
				
				warn Dumper $all_annex;
				die();
			}
			elsif ($nb_alt eq "?" ){
					warn Dumper $self->annex()->{$patient->id};
					warn $method;
					warn $self->id;
				die();
				
			}
			else {
				$type = "ho" if $all_annex->{ho};
				$sequence_info = "ho(" if $all_annex->{ho};
			#	warn $nb_ref." ".$nb_alt;
				$nb_ref += $nb_all_other_mut;
				my $sum = $nb_ref + $nb_alt;
				if ($sum >0){
		 			$pc = int ($nb_alt *100/($sum));
		 			push(@p_ratio,$pc);
				}
			$sequence_info .= $nb_ref."/".$nb_alt.")";
	
		}
		$sequence_info = $method_name.":".$sequence_info;
		#$pc = $method_name.":".$pc."%";
		$seq_infos->{$method_name}->{pourcent} = $pc;
		$seq_infos->{$method_name}->{status} = $type;
		$seq_infos->{$method_name}->{nb_ref} = $nb_ref;
		$seq_infos->{$method_name}->{nb_alt} = $nb_alt;
		$seq_infos->{$method_name}->{text} = $sequence_info;
		push(@{$self->{seq_infos}->{$pid}->{m}},$method_name);
		push(@{$self->{seq_infos}->{$pid}->{g}},$type);
		push(@{$self->{seq_infos}->{$pid}->{nr}},$nb_ref);
		push(@{$self->{seq_infos}->{$pid}->{na}},$nb_alt);
		#$seq_infos->{$method_name}->{text2} = "$method_name:";
	#	die($sequence_info)  if $self->isCnv;
	}
	
	$self->{seq_details}->{$pid} = $seq_infos;
   return $seq_infos;
	
}



sub getTranscripts {
	my ($self,$gene) = @_;
	if ($gene){
		my $ids = $self->transcripts_object();
		my $id_genes = $gene->transcripts_object();
		my $results;
		foreach my $tid (keys %$id_genes){
			$results->{$tid}= undef if exists $ids->{$tid};
		}
		return $self->getProject()->myflushobjects($results, "transcripts");
	}
	return $self->SUPER::getTranscripts();
	die();
}

sub check_no_var_allele_in_bam {
	my ($self, $patient, $limit) = @_;
	confess("\n\nERROR: need GenBoPatient object for GenBoVariant::check_no_var_allele_in_bam() method. Die.\n\n") unless ($patient);
	$limit = 2 unless ($limit);
	foreach my $bam_file (@{$patient->getBamFiles()}) {
		my $sam = $patient->bio_db_sam();
		my $chr_name = $self->getChromosome()->fasta_name();
		my $start = $self->start();
		my $end = $self->end();
		my $all_var = $self->var_allele();
		$all_var = '+' if ($self->isInsertion());
		$all_var = '-' if ($self->isDeletion());
		my $count = 0;
		my $nb_mut = 0;
		my $callback = sub {
			my ($seqid, $pos1, $pileups) = @_;
			return if ($pos1 ne $start);
			if (scalar(@$pileups) < 3) {
				$nb_mut = 99;
				return;
			}
			foreach my $pileup (@$pileups){
				my $b = $pileup->alignment;
				my $qbase = substr($b->qseq,$pileup->qpos,1);
				$nb_mut++ if ($qbase eq $all_var);
				if ($all_var eq "-"){ $nb_mut ++ if $pileup->indel < 0; }
				elsif ($all_var eq "+") { $nb_mut ++ if $pileup->indel > 0; }
				last if ($nb_mut >= $limit);
			}
		};
		$sam->fast_pileup("chr$chr_name:$start-$end", $callback);
		return if ($nb_mut >= $limit);
		return 1;
	}
	return;
}

#########
# DejaVu
###########
has in_this_run_patients => (
        is              => 'rw',
        #required=> 1,
        lazy    => 1,
        default => sub {
                my $self = shift;
                
                #if (exists  $self->dejaVuInfosForDiag2("in_this_run_patients");}){
                #	return  $self->dejaVuInfosForDiag2->{in_this_run_patients};
                #}
                return scalar(@{$self->getPatients});
        }
);

has in_this_run_ratio => (
        is              => 'rw',
        #required=> 1,
        lazy    => 1,
        default => sub {
                my $self = shift;
                return 0 if $self->in_this_run_patients  <= 1;
                my $nb;
                # if (exists  $self->dejaVuInfosForDiag2->{total_in_this_run_patients}){
                  
                #	$nb =   $self->dejaVuInfosForDiag2->{total_in_this_run_patients};
                #}
                #else {
               		$nb = scalar(@{$self->project->getPatients})-1;
                #}
                confess() if $nb == 0;
                return  ($self->in_this_run_patients / $nb);
        }
);

sub similar_projects {

                my $self = shift;
              return  $self->dejaVuInfosForDiag2("similar_projects");
        }


sub  similar_patients{
                my $self = shift;
                return $self->dejaVuInfosForDiag2("similar_patients")
        }
sub similar_patients_ho {
                my $self = shift;
  				return $self->dejaVuInfosForDiag2("similar_patients_ho");
 }
 
 
sub other_projects  {
        my ($self) = @_;
     	return $self->dejaVuInfosForDiag2("other_projects");
}
sub other_projects_ho  {
        my ($self) = @_;
     	return $self->dejaVuInfosForDiag2("other_projects_ho");
}
sub other_patients{
          my $self = shift;
           return $self->dejaVuInfosForDiag2("other_patients");
     }

sub other_patients_ho {
                my $self = shift;
             return $self->dejaVuInfosForDiag2("other_patients_ho");
        }

sub exome_projects {
                my $self = shift;
               return $self->dejaVuInfosForDiag2("exome_projects");
}

sub exome_patients {
      
                my $self = shift;
               return  $self->dejaVuInfosForDiag2("exome_patients")
}


sub  exome_patients_ho {
                my $self = shift;
                return   $self->dejaVuInfosForDiag2("exome_patients_ho");
 }


sub total_exome_projects  {
                my $self = shift;
                  return   $self->dejaVuInfosForDiag2("total_exome_projects");
        }


sub total_exome_patients {
                my $self = shift;
                   return   $self->dejaVuInfosForDiag2("total_exome_patients");
        }


sub total_similar_projects {
              my $self = shift;
               return   $self->dejaVuInfosForDiag2("total_similar_projects");
}
        
sub total_similar_patients  {
                my $self = shift;
               return   $self->dejaVuInfosForDiag2("total_similar_patients");
        }
        
sub dejavu_hash_projects_patients {
	my ($self) = @_;
	my $hres;
	my $h_dv = $self->getChromosome->rocks_dejavu->dejavu($self->rocksdb_id);
	warn $self->getChromosome->rocks_dejavu->dir();
	warn Dumper $h_dv;
	foreach my $proj_id (keys %{$h_dv}) {
		my $proj_name = $self->buffer->getProjectNameFromId($proj_id);
		
		my $h_pat_proj;
		foreach my $h (@{$self->buffer->getQuery->getPatients($proj_id)}) {
			$h_pat_proj->{$h->{patient_id}} = $h->{name};
		}
		foreach my $pat_id (@{$h_dv->{$proj_id}->{patients}}) {
			$hres->{$proj_name}->{patients}->{$h_pat_proj->{$pat_id}} = undef;
		}
		$hres->{$proj_name}->{nb_he} = $h_dv->{$proj_id}->{he};
		$hres->{$proj_name}->{nb_ho} = $h_dv->{$proj_id}->{ho};
	}
	
	
	return $hres;
}

sub dejaVuInfosForDiag2 {
	my ($self,$key) = @_;
	unless (exists $self->{array_dejavu}) {
		#my $hash = $self->getChromosome->getDejaVuPhenotypeInfosForDiagforVariant($self);
		my $hash = $self->getChromosome->getDejaVuInfosForDiagforVariant($self);
		
		$self->{array_dejavu} = $self->buffer->hash_to_array_dejavu($hash);
	}

	return $self->{array_dejavu} unless $key;
	my $index = $self->buffer->index_dejavu($key);
	return  $self->{array_dejavu}->[$index];
}

	
#	unless  (exists $self->{dejaVuInfosForDiag2}){
#		$self->{dejaVuInfosForDiag2} = $self->getProject->getDejaVuInfosForDiag($self->id());
#	}
#	if (exists $self->{dejaVuInfosForDiag2}){
#		$self->{array_dejavu} 
#	}
#	if (exists $self->{array_dejavu}){
#		my $index = $self->buffer->index_dejavu($key);
#		return  $self->{array_dejavu}->[$index]
#	}
#	
#	if (exisst $self->{dejaVuInfosForDiag2})
#	if (exists array)
#}
     
#has dejaVuInfosForDiag2 => ( 
#	is              => 'rw',
#	lazy    => 1,
#	default => sub {
#		my $self = shift;
#		return $self->getProject->getDejaVuInfosForDiag($self->id());
#	}
#);
     


sub string_nomenclature {
	my($self,$st,$seqv) =@_;
	if(length ($seqv) > 50 ){
		return "$st".length($seqv);
	} 
	else 
	{
		return "$st".$seqv;
	}
}

#################
#
# PURGE METHODS
#
##################
sub purge_public_data {
	my $self = shift;
	delete $self->{scaled_score_frequence_public};
	delete $self->{score_frequence_public};
	delete $self->{gnomad};
	delete $self->{isPublic};
	delete $self->{name};
}


sub purge_deja_vu {
	 my $self = shift;
	 delete $self->{dejaVuInfosForDiag};
}

#################
#
#
# SPLIT READ INFORMATION
#
##################

has split_read_infos =>(
	is =>'ro',	
	lazy=>1,
	default=> sub {
		my $self = shift;
		 my $hash; 
		return {} unless $self->isSrPr();
		foreach my $pat (@{$self->getPatients}){
			foreach my $patient (@{$pat->getFamily()->getMembers}){
				my $pid = $patient->id;
				my $pr0 ="-";
				my $pr1 ="-";
				my $sr0;
				my $sr1;
				my $srq_end = -1 ;
				my $srq_start = -1 ;
				my $equality = -1;
				#$srq_end = $patient->mean_align_quality($self->getChromosome, $self->start+$self->length,$self->start+$self->length);
				#$srq_start = $patient->mean_align_quality($self->getChromosome, $self->start, $self->start);
				#$equality = $patient->mean_align_quality($self->getChromosome, $self->start, $self->end);
				my $mean_dp =  int($patient->meanDepth($self->getChromosome->name, $self->start, $self->end));
				#$self->{seq}->{$pid}->{sr_quality_end} = $patient->mean_align_quality($self->getChromosome, $self->start+$self->length,$self->start+$self->length);
				unless ($self->existsPatient($patient)){
					
					my ($a,$b,$c) = (0,0,0);#$patient->sr_raw($self->getChromosome, $self->start);
					$sr0  =  $a;
					$sr1 = int(($b+$c)/2); 
		 			($a,$b,$c) = (0,0,0);#$patient->sr_raw($self->getChromosome, $self->end);
					$sr0  +=  $a;
					$sr1 += int(($b+$c)/2); 
					$sr0  = int($sr0/2);
				}
				else {
					unless ( exists $self->{annex}->{$patient->id()}->{sr}) {
					#	warn "ATTENTION PAS  DE SR "." ".$self->start." ".$self->getChromosome->name." ".Dumper( $self->annex()->{$patient->id()});
					#	die();
					
						$hash->{$pid} = ["-1","-1","-1","-1",$srq_end,$srq_start,$equality];
						
					}
					else {
			#	confess(Dumper ($self->annex)." ".$self->start." ".$self->getChromosome->name) unless  exists $self->annex()->{$patient->id()}->{sr};
					($sr0,$sr1) = split(",",$self->{annex}->{$patient->id()}->{sr});
					if (exists  $self->{annex}->{$patient->id()}->{pr}){
						($pr0,$pr1) = split(",",$self->{annex}->{$patient->id()}->{pr});
					}
					}
				}
				$hash->{$pid} = [$sr0,$sr1,$pr0,$pr1,$srq_end,$srq_start,$equality];
			}
						
		}
			return $hash;
	},
);

sub sr0 {
	my ($self, $patient) = @_;
	my $pid = $patient->id;
		my $v = $self->split_read_infos->{$pid}->[0];
	$v = 0 unless $v ;
	$v = "-1" if $v eq "-";
	return $v;
}
sub sr1 {
	my ($self, $patient,$debug) = @_;
	my $pid = $patient->id;
	my $v = $self->split_read_infos->{$pid}->[1];
	$v = 0 unless $v ;
	$v = "-1" if $v eq "-";
	return $v;
}
sub pr0 {
	my ($self, $patient) = @_;
	my $pid = $patient->id;
	my $v = $self->split_read_infos->{$pid}->[2];
	$v = 0 unless $v ;
	$v = "-1" if $v eq "-";
	return $v;
}
sub pr1 {
	my ($self, $patient) = @_;
	confess($patient) unless $patient;
	my $pid = $patient->id;
		my $v = $self->split_read_infos->{$pid}->[3];
	$v = 0 unless $v ;
	$v = "-1" if $v eq "-";
	return $v;
}

sub sr_align_quality {
	my ($self, $patient) = @_;
	confess();
	my $pid = $patient->id;
	return int(($self->split_read_infos->{$pid}->[4]+$self->split_read_infos->[5])/2);
}

sub event_align_quality {
	my ($self, $patient) = @_;
		confess();
	my $pid = $patient->id;
	return $self->split_read_infos->{$pid}->[6];
}

sub sr {
	my ($self, $patient) = @_;
	my $pid = $patient->id;
	return "-" unless $self->isSrPr;
	my $sr0 = $self->sr0($patient);
	return "-" if not $sr0 or $sr0 eq "-";
	my $sr1 = $self->sr1($patient);
	$sr0 = '-' unless ($sr0);
	$sr1 = '-' unless ($sr1);
	return  $sr0.",".$sr1;

}

sub pr {
	my ($self, $patient) = @_;
	my $pid = $patient->id ;
	return "-"  unless $self->isSrPr;
	return "-" if $self->pr0($patient) eq "-";
	return  $self->pr0($patient).",". $self->pr1($patient);
}

sub is_manta {
	my ($self) = @_;
	#my $pid = $p->id;
	return $self->isSrPr();
	
}

has is_forced_viewing => (
	is              => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		foreach my $gene (@{$self->getGenes()}) {
			next if not $gene->hash_variants_forced_viewing();
			next if not exists $gene->hash_variants_forced_viewing->{$self->id()};
			return $gene->hash_variants_forced_viewing->{$self->id()};
		}
		return;
	}
);

my $h_parquet_models_id;
$h_parquet_models_id->{1} = 'solo';
$h_parquet_models_id->{2} = 'father';
$h_parquet_models_id->{4} = 'mother';
$h_parquet_models_id->{8} = 'both';
$h_parquet_models_id->{16} = 'is_parent';
$h_parquet_models_id->{32} = 'recessif';
$h_parquet_models_id->{64} = 'dominant';
$h_parquet_models_id->{128} = 'denovo';
$h_parquet_models_id->{256} = 'strict_denovo';
$h_parquet_models_id->{512} = 'error';

sub infos_dejavu_parquet {
	my ($self, $project) = @_;
	$project = $self->getProject() if not $project;
	my $dir_parquet = $project->buffer->dejavu_parquet_dir();	
	my $parquet = $dir_parquet.'/'.$project->name().'.'.$project->id().'.parquet';
	return undef if not -e $parquet;
#	confess("\n\nERROR $parquet file doesn't exist. die \n\n") if not -e $parquet;
	
	my $sql = "SELECT chr38 as chr, pos38 as pos, * FROM read_parquet(['".$parquet."'])";

	my $find_pos_s = $self->start() - (20 + $self->length());
	my $find_pos_e = $self->start() + (20 + $self->length());
	if ($self->getProject->current_genome_version() eq 'HG38') {
		$sql .= " WHERE chr38='".$self->getChromosome->id()."' and pos38 BETWEEN '".$find_pos_s."' and '".$find_pos_e."';" ;
	}
	elsif ($self->getProject->current_genome_version() eq 'HG19') {
		$sql = "SELECT chr19 as chr ,pos19 as pos,  * FROM read_parquet(['".$parquet."'])" ;
		$sql .= " WHERE chr19='".$self->getChromosome->id()."' and pos19 BETWEEN '".$find_pos_s."' and '".$find_pos_e."';" ;
	}
	else { confess("\n\nERROR VAR PARQUET DEJAVU !\n\n"); }
	my $duckdb = $project->buffer->software('duckdb');
	my $cmd = qq{set +H | $duckdb -json -c "$sql"};
	my $json_duckdb = `$cmd`;
	if ($json_duckdb) {
		my $decode = decode_json $json_duckdb;
		my $h_by_proj;
		my $rocks = $self->rocksdb_id();
		my ($pos, $alt) = split('!', $rocks);
		$pos = int($pos);
		foreach my $h (@$decode) {
			
			next if $h->{'chr'} ne $self->getChromosome->id;
			next if $h->{'pos'} ne $pos;
			next if $h->{'allele'} ne $alt;
			my $nb_he = $h->{he};
			my ($h_infos_patients, $h_tmp_pat);
			my $nb_pat = 0;
			foreach my $pat_id (unpack("w*",decode_base64($h->{patients}))) {
				$nb_pat++;
				$h_infos_patients->{$nb_pat}->{id} = $pat_id;
				$h_tmp_pat->{$pat_id} = $nb_pat;
			}
			my $i = 0;
			$nb_pat = 1;
			foreach my $info (unpack("w*",decode_base64($h->{dp_ratios}))) {
				$i++;
				if ($i == 1) { $h_infos_patients->{$nb_pat}->{dp} = $info; }
				elsif ($i == 2) { $h_infos_patients->{$nb_pat}->{ratio} = $info; }
				elsif ($i == 3) {
		    		my $model = 'error';
		    		$model = $h_parquet_models_id->{$info} if exists $h_parquet_models_id->{$info};
					$h_infos_patients->{$nb_pat}->{model} = $model;
					$i = 0;
					$nb_pat++;
				}
			}
			
			foreach my $pat (@{$project->getPatients()}) {
				next if not exists $h_tmp_pat->{$pat->id};
				$h->{patients_infos}->{$pat->name}->{patient_id} = $pat->id();
				$h->{patients_infos}->{$pat->name}->{patient_parquet_id} = $h_tmp_pat->{$pat->id};
				if (int($h_tmp_pat->{$pat->id}) <= $nb_he) { $h->{patients_infos}->{$pat->name}->{heho} = 'he'; }
				else { $h->{patients_infos}->{$pat->name}->{heho} = 'ho'; }
				$h->{patients_infos}->{$pat->name}->{model} = $h_infos_patients->{$h_tmp_pat->{$pat->id}}->{model};
				$h->{patients_infos}->{$pat->name}->{ratio} = $h_infos_patients->{$h_tmp_pat->{$pat->id}}->{ratio};
				$h->{patients_infos}->{$pat->name}->{dp} = $h_infos_patients->{$h_tmp_pat->{$pat->id}}->{dp};
			}
			return $h;
		}
	}	
	return undef;
}

1;
