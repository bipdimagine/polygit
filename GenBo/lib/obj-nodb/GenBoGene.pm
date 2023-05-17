package GenBoGene;

use strict;
use Vcf;
use Moo;
use Data::Dumper;
use Config::Std;
 use List::Util qw( max );
extends "GenBoGenomic";

 
has isGene => (
	is		=> 'ro',
	default	=> 1,
);

has kyotoId => (
	is		=> 'ro',
#	required=> 1,
);

has ensg => (
	is		=> 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		if ($self->is_intergenic) { return $self->name(); }
		my ($ensg, $chr_id) = split('_', $self->name());
		return $ensg;
	},
);
	
has external_name => (
	is		=> 'ro',
	#required=> 1,
);


has strand => (
	is		=> 'ro',
	required=> 1,
);

has transcripts => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { [] },
);

has description => (
	is		=> 'ro',
	lazy	=> 1,
	default => 'NA',
);

has pLI => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		my $z = $self->project->lmdbPLI->get($self->name());
		return "-" unless $z;
		 return $z;
		 },
);

has hpo => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		my $z = $self->project->lmdbHPO_genes_to_phenotypes->get($self->external_name());
		return undef unless $z;
		return $z;
	},
);

has hpo_phenotypes_ids => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		my @list_ids; 
		@list_ids = keys %{$self->hpo()} if ($self->hpo());
		return \@list_ids
	},
);

has hpo_phenotypes_names => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		my @list_names; 
		@list_names = values %{$self->hpo()} if ($self->hpo());
		return \@list_names
	},
);

has is_HGMD_DM => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return 1 if (exists $self->chromosome->hash_hgmd_class_DM_genes_name->{$self->external_name()});
		return;
	},
);


	
 
has annotations_tree => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->project->liteIntervalTree->get("functional_annotations",$self->id);
		
		#return  $self->project->liteIntspan->get("intspan_genes",$self->name);
	}
);


sub getMainTranscripts{
	my ($self) = @_;
#	if ($self->project->isDiagnostic){
#		my (@t) = grep {$_->getGene->id eq $self->id}  @{$self->project->getListGenBoTranscripts()};
#		return \@t if @t;
#	}
	my $trs = $self->project->lmdbMainTranscripts->get($self->id);
	return $self->project->myflushobjects($trs,"transcripts");
	
	
	
}

###### SET OBJECTS #####
sub getTanscriptsAnnotations {
	my ($self,$start,$end) = @_;
	my $h ={};
	my $res = $self->annotations_tree->fetch($start,$end+1);
	unless (@$res){
#		# TODO: cas des genes avec aucun transcript (exemple: NGS2019_2365 ENSG00000215791_1 -> AL645728.2 0 transcript et 1 sur Ensembl mais SANS Protein)
#		if (scalar(@{$self->getTranscripts()}) == 0) {
#			return $h;
#		}
#		# TODO: cas gene avec des transcripts mais aucun avec une protein
#		# TODO: exemple NGS2015_0794, gene RP11-33B1.1, ensembl -> http://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000245958;r=4:120375946-120473180
#		my $nbProt = 0;
#		foreach my $tr (@{$self->getTranscripts()}) {
#			$nbProt++ if ($tr->getProtein());
#		}
#		if ($nbProt == 0) {
#			return $h;
#		}
		my $nearest = $self->annotations_tree->fetch_nearest_up($end);
		unless ($nearest){
			$nearest = $self->annotations_tree->fetch_nearest_down($start);
		} 
		$res=[$nearest];
	}
	unless (@$res) {
		die($self->id." ".$self->name." ".$start." ".$end);
	}
	my $bug;
	foreach my $r  (@$res){
	
		$bug=1 unless $r->[0];
		push(@{$h->{$r->[0]}},$r->[1]);
	}
	confess() if $bug;
	return $h;
}



## pourquoi cette fois c'est transcripts et pas transcriptIds tout court ????

sub xref {
	my $self = shift;
	die();
	return $self->external_name();
}

sub setTranscripts {
	my $self = shift;
	my $hTranscriptsId = {};
	my $lTranscriptsId = $self->transcripts();
	foreach my $id (@$lTranscriptsId) {
			unless ($id =~/_/){
				$id.="_".$self->getChromosome->name();
			}
			$hTranscriptsId->{$id} = undef;
	}
	return $hTranscriptsId;
}

has omim_id => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		return $self->omim->{omim_id} if exists $self->omim->{omim_id};
		return undef;
		 },
);

has phenotypes  => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;

		my $pheno = join(' ;', @{$self->array_phenotypes});
		return $pheno;
	},
);

has hash_phenotypes => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		my $hPhen;
		if (exists $self->omim->{version}){
			$hPhen->{omim} = $self->omim->{phenotypes}->{omim} if (exists $self->omim->{phenotypes}->{omim});
			$hPhen->{hgmd} = $self->hgmd_disease() if ($self->hgmd and $self->hgmd_disease);
			$hPhen->{ddg2p} = $self->omim->{phenotypes}->{ddg2p} if (exists $self->omim->{phenotypes}->{ddg2p});
			$hPhen->{orphanet} = $self->omim->{phenotypes}->{orphanet} if (exists $self->omim->{phenotypes}->{orphanet});
			$hPhen->{hpo} = $self->hpo_phenotypes_names() if ($self->hpo());
		}
		else {
			my @lOmim = split(";",$self->omim->{phenotype}->{omim}) if (exists $self->omim->{phenotype}->{omim} and $self->omim->{phenotype}->{omim});
			$hPhen->{omim} = \@lOmim if (scalar @lOmim > 0);
			$hPhen->{hgmd} = $self->hgmd_disease() if ($self->hgmd and $self->hgmd_disease);
			my @lEnsembl = split(";", $self->omim->{phenotype}->{ensembl}) if (exists $self->omim->{phenotype}->{ensembl});
			$hPhen->{ensembl} = \@lEnsembl if (scalar @lEnsembl > 0);
			$hPhen->{hpo} = $self->hpo_phenotypes_names() if ($self->hpo());
		}
		return $hPhen;
	},
);

has array_phenotypes  => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		my @lPhen;
#		foreach my $cat (keys %{$self->hash_phenotypes()}) {
#			push(@lPhen, @{$self->hash_phenotypes->{$cat}});
#		}
		if (exists $self->omim->{version}){
			push(@lPhen, @{$self->omim->{phenotypes}->{omim}} ) if (exists $self->omim->{phenotypes}->{omim});
			push(@lPhen, @{$self->hgmd_disease} ) if ($self->hgmd and $self->hgmd_disease);
			push(@lPhen, @{$self->omim->{phenotypes}->{ddg2p}} ) if (exists $self->omim->{phenotypes}->{ddg2p});
			push(@lPhen, @{$self->omim->{phenotypes}->{orphanet}}) if (exists $self->omim->{phenotypes}->{orphanet});
			push(@lPhen,  sort @{$self->hpo_phenotypes_names()}) if ($self->hpo());
			return \@lPhen;
		}
		else {
			push(@lPhen, split(";",$self->omim->{phenotype}->{omim}) )   if (exists $self->omim->{phenotype}->{omim} and $self->omim->{phenotype}->{omim});
			push(@lPhen,  @{$self->hgmd_disease}) if ($self->hgmd and $self->hgmd_disease);
			push(@lPhen,   split(";",$self->omim->{phenotype}->{ensembl})) if (exists   $self->omim->{phenotype}->{ensembl});
			push(@lPhen,  sort @{$self->hpo_phenotypes_names()}) if ($self->hpo());
		}
		#
		return \@lPhen;
	},
);

sub polyquery_phenotypes {
	my ($self) = @_;
	my @lPhen;
	push(@lPhen, 'Omim: '.join(', ', @{$self->omim->{phenotypes}->{omim}}) ) if (exists $self->omim->{phenotypes}->{omim});
	push(@lPhen, 'Hgmd: '.join(', ', @{$self->hgmd_disease}) ) if ($self->hgmd and $self->hgmd_disease);
	push(@lPhen, 'Ddg2p: '.join(', ', @{$self->omim->{phenotypes}->{ddg2p}}) ) if (exists $self->omim->{phenotypes}->{ddg2p});
	push(@lPhen, 'Hpo: '.join(', ', @{$self->hpo_phenotypes_names()})) if ($self->hpo());
	return join(' ; ',@lPhen);
}

sub polyviewer_phentotypes {
	my ($self) = @_;
	my $array_phenotypes = $self->array_phenotypes;
	 my $pheno;
	my $nb_other_terms = 0;
	if (@$array_phenotypes){
	   #	$pheno = $self->array_phenotypes->[0] ;
	    foreach my $p (@{$self->array_phenotypes}){
	    	$pheno .= $p;
	    	last if length($pheno)>100;
	    	$pheno .= " | ";
	    }
	   	$nb_other_terms = scalar(@$array_phenotypes) - 1;
	 }
	 else {
	   	$pheno = $self->omim->{title} if $self->omim->{title} ;
	   	my $to;
	 	if ( $self->description){
	   	($pheno,$to) = split(/\[/,$self->description) unless $pheno;
	 	} #&& $self->description;
	 	$pheno = "-" unless $pheno;
	   }
	return ($pheno,$nb_other_terms);
}

has short_phenotypes  => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		#return "-" unless exists $self->omim->{phenotype};
		my $pheno = "";
		if (exists $self->omim->{version}){
			if ($self->omim->{phenotypes}->{omim}){
				$pheno = join (";",@{$self->omim->{phenotypes}->{omim}});
				return $pheno;
			}
			if ($self->omim->{phenotypes}->{ddg2p}){
				$pheno = "ddg2p:".join (";",@{$self->omim->{phenotypes}->{ddg2p}});
				return $pheno;
			}
			if ($self->omim->{phenotypes}->{orphanet}){
				$pheno = "oprh:".join (";",@{$self->omim->{phenotypes}->{orphanet}});
				return $pheno;
			}
			if (@{$self->hgmd_disease}){
				$pheno = join (";",@{$self->hgmd_disease});
				return $pheno;
			}
			return $self->omim->{title} if $self->omim->{title};
			return ""; 
		
		}
		else {
		if ($self->omim->{phenotype}->{omim}){
			$pheno = join (";",$self->omim->{phenotype}->{omim});
		}
		if (@{$self->hgmd_disease}){
			#$pheno = "HGMD : ";
			$pheno = join (";",@{$self->hgmd_disease});
			return $pheno;
		}
			
		 return " Ens: ".$self->omim->{phenotype}->{ensembl} if exists   $self->omim->{phenotype}->{ensembl};
		 return " Omim: ".$self->omim->{phenotype}->{omim} if exists   $self->omim->{phenotype}->{omim};
		return $pheno;
		}
		},
);


has omim_inheritance => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		return "-" unless exists $self->omim->{inheritance}->{omim};
		my $in="";
		$in = "X-linked " if exists $self->omim->{inheritance}->{omim}->{"X-linked"};
		$in = "Y-linked " if exists $self->omim->{inheritance}->{omim}->{"Y-linked"};
		$in .=  "recessive " if exists  $self->omim->{inheritance}->{omim}->{recessive};
		$in .=  "dominant " if exists  $self->omim->{inheritance}->{omim}->{dominant};
		return $in;
		 },
);

has inheritance =>(
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		my $hi = $self->hgmd_inheritance;
		my $ho = $self->omim_inheritance;
		return $hi if ($hi eq $ho);
		return $hi."/".$ho ;#if ($hi eq $ho);
		
		 },
);
has omim => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		my $hash = $self->project->lmdbOmim->get($self->name);
		$hash = $self->project->lmdbOmim->get($self->id) unless $hash;
		delete    $hash->{phenotype}->{ensembl} if  exists $hash->{phenotype}->{ensembl} && $hash->{phenotype}->{ensembl} eq ".";
		unless (exists $hash->{is_morbid}) {
			eval {
				$hash->{is_morbid} = 1 if (exists $self->buffer->hash_genes_omim_morbid->{$self->external_name()});
			};
			if ($@) { $hash->{is_morbid} = undef; }
		}
		if (exists $hash->{phenotypes}->{omim_root}){
			$hash->{phenotypes}->{omim} = $hash->{phenotypes}->{omim_root};
			#die();
			delete $hash->{phenotypes}->{omim_root};
		}
		return {} unless $hash;
		return $hash;
	 },
);

has is_omim_morbid => (
	is		=> 'ro',
	lazy	=> 1,
	default =>  sub { 
		my $self = shift;
		return 1 if (exists $self->omim->{is_morbid} and $self->omim->{is_morbid});
		return;
	 },
);

has is_omim_new => (
	is		=> 'ro',
	lazy	=> 1,
	default =>  sub { 
		my $self = shift;
		return 1 if (exists $self->omim->{'new'} and $self->omim->{'new'});
		return;
	 },
);

has is_omim_morbid_new => (
	is		=> 'ro',
	lazy	=> 1,
	default =>  sub { 
		my $self = shift;
	
		return 1 if (exists $self->omim->{'new_morbid'} and $self->omim->{'new_morbid'});
		return;
	 },
);

has hgmd => (
	is		=> 'ro',
	lazy => 1,
	default => sub { 
		my $self = shift;
		my $query = $self->buffer->queryHgmd();
		my $h;
		 $h = $query->getDataHGMDPro_gene_infos($self->external_name());
		return unless ($h);
		return $h->{$self->external_name()};
	 },
);

has hgmd_inheritance => (
	is		=> 'ro',
	lazy => 1,
	default => sub { 
		my $self = shift;
		return "" unless $self->hgmd();
		
		my $in="";
		
		$in = "X-linked " if $self->hgmd->{expected_inheritance}  && $self->getChromosome->name() eq "X";
		
		$in = "Y-linked " if $self->hgmd->{expected_inheritance} &&  $self->getChromosome->name() eq "YL";
		$in .=  "recessive " if $self->hgmd->{expected_inheritance} =~ /R/;
		$in .=  "dominant "if $self->hgmd->{expected_inheritance} =~ /D/;
		$in = "-" unless $in;
		return $in;
		
	 },
);

has hgmd_description => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		return unless $self->hgmd();
		return $self->hgmd->{genename};
	 },
);


has hgmd_other_description => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		return unless $self->hgmd();
		my @list;
		if ($self->hgmd->{altname}) {
			@list = split('\|', $self->hgmd->{altname});
		}
		return \@list;
	 },
);

has hgmd_disease => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		return [] unless $self->hgmd();
		my @list;
		if ($self->hgmd->{disease}) {
			@list = split('\|', $self->hgmd->{disease});
		}
		return \@list;
	 },
);
has summary_disease => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		return "" unless $self->hgmd();
		my %word;
		my $text;
		foreach my $d (@{$self->hgmd_disease}){
			$d =~s/\,/ /;
			$d =~s/\;/ /;
			$d =~s/&/ /;
			$d =~s/with/ /;
			$d =~s/of/ /;
			$d =~s/and/ /;
			map{$word{lc($_)}++} split(" ",$d);
		}
		my $limit = 4;
		
		foreach my $k (sort{$word{$b} <=> $word{$a}}keys %word){
			next if length($k) <3;
			 next if lc($k) eq "syndrome"; 
			 next if lc($k) eq "association"; 
			  next if lc($k) eq "recessive"; 
			  $text.= $k." ";
			$limit --;
			last if $limit ==0;
		}
		warn "text : ". $text;
		warn join(";",@{$self->hgmd_disease})."\n";
		warn "-----------";
		
	 },
);
has hgmd_refseq => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		return unless $self->hgmd();
		return $self->hgmd->{refseq};
	 },
);
has hgmd_go_terms_name => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		return unless $self->hgmd();
		my @list;
		if ($self->hgmd->{go_terms_name}) {
			@list = split('\|', $self->hgmd->{go_terms_name});
		}
		return \@list;
	 },
);


has hgmd_go_terms_acc => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		return unless $self->hgmd();
		my @list;
		if ($self->hgmd->{go_terms_acc}) {
			@list = split('\|', $self->hgmd->{go_terms_acc});
		}
		return \@list;
	 },
);

has hash_variants_dejavu => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		return $self->getProject->get_deja_vu_from_position2($self->getChromosome->id(), $self->start(), $self->end);
	 },
);

sub setPhenotypes {
	my $self = shift;
	my %hids;
	foreach my $panel (@{$self->getPanels}){
		map {$hids{$_->id}++;} @{$panel->getPhenotypes}; 
	}
	return \%hids;
}

sub setPanels {
	my $self = shift;
	my %hids;
	my $query = $self->buffer->queryPanel();
	my $hs = $query->getPanelsForGeneName($self->external_name);
	my %hIds;
	foreach my $k (keys %$hs){
		$hIds{$hs->{$k}->{panel_id}} ++;
	}
	
	return \%hIds;
}

#$ogene->buffer->queryPanel()->getPanelsForGeneName($ogene->external_name);


sub getExons {
	my ($self) = @_;
	my $exons;
	foreach my $tr (@{$self->getMainTranscripts}){
		push(@$exons, @{$tr->getExons});
	}
	return $exons;
}

 

has panels_name => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		my $panels = $self->buffer->queryPanel()->getPanelsForGeneName($self->external_name);
		return $panels;
		
	 },
); 

has  hpanels => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub { 
		my $self = shift;
		my $panels = $self->buffer->queryPanel()->getPanelsForGeneName($self->external_name);
		
		return $panels;
		
	 },
); 

sub score {
	my ($self) = @_;
	return $self->{score} if exists $self->{score};
	my $score =0;
	$score += 0.5 if $self->pLI ne "-" && $self->pLI > 0.95;
	$score += 1 if $self->is_omim_morbid();

	my $polyScore = $self->project->lmdbpolyScore->get($self->id);
	my $phenotypes = $self->getProject->getPhenotypes();
	my @pscore;
	if ($phenotypes){
		foreach my $ph (@$phenotypes){
			my $add = 0;
			# $add = 0.5 if $self->pLI ne "-" && $self->pLI ==1;
			push(@pscore,($polyScore->{$ph->id}->{score}+$add)) if exists $polyScore->{$ph->id}->{score};
		}
	}
	else {
		map {push(@pscore,$polyScore->{$_}->{score})} keys %$polyScore;
	}
	#warn Dumper(@pscore);

	$score += max(@pscore) if @pscore;
	$self->{score} = $score;
	return $score;
}
 
sub raw_score {
	my ($self) = @_;
	my $debug;
	$debug = 1 if  $self->external_name eq "MED12";
		my $score = 0;
		$score += 0.5 if $self->pLI ne "-" && $self->pLI > 0.95;
		my $gpheno = $self->phenotypes();
		
		my $pscore =0;
		$score += 1 if $self->is_omim_morbid();
		warn $score if $debug;
		my $apheno = $self->getPhenotypes();
	#	warn $score if $debug;
		foreach my $pheno (@{$self->getProject->getPhenotypes}){
			warn $pheno->name() if $debug;
			foreach my $k  (@{$pheno->keywords()}) {
				my $lk = lc($k);
				if ($gpheno =~ /$lk/){
					$score += 1 ;
					last;
				}
			}
			warn $score if $debug;
			if (exists $pheno->statistic_genes->{$self->id}){
				$score += 1;
				my $p = ($pheno->statistic_genes->{$self->id} / $pheno->nb_panels) *100;
				$score += 0.5 if $p >= 25;
				$score += 0.5 if $p >= 50;
				$score += 0.5 if $p >= 75;
				$score += 0.5 if $p >= 90;
				$score += 0.5 if $p >= 100;
			}
			else {
				my @o = grep{$pheno->id ne $_->id;} @$apheno; 
				$score += 0.2 * (scalar(@o));
			}
			warn $score if $debug;
		} 
		warn $score if $debug;
		#die() if $debug;
		
		
		return $score;
}

1;