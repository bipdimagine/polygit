package PolyviewerVariant;
use strict;
use FindBin qw($Bin);
use Moo;
use Data::Dumper;



has impact_sorted => (
	is		=> 'ro',
	default => sub {
		return {
	"high" => "4",
	"moderate" =>"3",
	"low" =>"1",
};
	}
		
);

sub parseGnomadTable {
	my ($self,$html,$h) = @_;
	my $hpatients;
		#
		# HASH PATIENTS POUR LE CALLING 
		#
		my $te = HTML::TableExtract->new();
			$te->parse( $html);
			my $table      = $te->first_table_found;
			my @rows = $table->rows;
			my $zmodel;
			shift (@rows);
		
		foreach my $l (@rows) {
		#next if lc($l->[0]) eq "ac";
			my $z =0;	
			my $hemi;
			$hemi =1  if $self->chromosome eq "X" or $self->chromosome eq "Y";
			
			my $an = pop(@$l);
			$self->gnomad_an($an);
			next unless $an;
			$self->gnomad_ac($l->[$z]);
			$z++;
			$self->gnomad_ho($l->[$z]);
			
			$z++;
			if($hemi){
				$self->gnomad_ho_male($l->[$z]);
				$z++;
			}
			my $v1 ;
			my $v2;
			 ($v1,$v2) = split(" ",$l->[$z]) if $l->[$z];
			$self->gnomad_min_pop_name($v1);
			$self->gnomad_min_pop($v2);
			$z++;
			
		
			 ($v1,$v2) = split(" ",$l->[$z]);
			$z++;
			$self->gnomad_max_pop_name($v1);
			$self->gnomad_max_pop($v2);
			
					}
			
			return;		
} 

sub parseDejaVuTable {
	my ($self,$html,$h) = @_;
	my $hpatients;
		#
		# HASH PATIENTS POUR LE CALLING 
		#
		my $te = HTML::TableExtract->new();
			$te->parse( $html);
			my $table      = $te->first_table_found;
			my @rows = $table->rows;
			shift(@rows);
			my $l = shift(@rows);
			shift(@$l);
			$self->dejavu_other_projects(shift(@$l));# = $v->other_projects();
			$self->dejavu_other_patients(shift(@$l));# = $v->other_patients();
		$self->dejavu_other_patients_ho(shift(@$l));# = $v->other_patients_ho();
		$l = shift(@rows);
		shift(@$l);
			
		$self->dejavu_similar_projects(shift(@$l) );#= $v->similar_projects();
		$self->dejavu_similar_patients(shift(@$l));# = $v->similar_patients();
		$self->dejavu_similar_patients_ho(shift(@$l));# = $v->similar_patients_ho();
	
		if (@rows){
			die();
		$l = shift(@rows);
		shift(@$l);
		$self->dejavu_this_run_patients(shift(@$l));# = '-';
		}	
		return;	
			
	
} 

sub parseTrioTableVariant {
	my ($self,$html,$project) = @_;
	my $hpatients;
	#
	# HASH PATIENTS POUR LE CALLING 
	#
	
	my $te = HTML::TableExtract->new();
	$te->parse( $html);
	my $table      = $te->first_table_found;
	my @rows = $table->rows;
	my $zmodel ="?";
		foreach my $l (@rows) {
		next if $l->[0] eq "pat";
		my $p = $project->getPatient( $l->[0] );
		$hpatients->{ $p->id }->{gt} = $l->[2];
		$l->[3] =~ s/%// if $l->[3];
		$hpatients->{ $p->id }->{pc} = $l->[3];
		$l->[4] = '-' unless ($l->[4]);
		$l->[4] =~ s/\-//;
		$hpatients->{ $p->id }->{dp} = $l->[4];
		if ($p->isMother ){
			if  ($hpatients->{ $p->id }->{gt} eq "-"){
				$hpatients->{ $p->id }->{model} = "-m" ;
			}
			else {
				$zmodel = "mother";
				$hpatients->{ $p->id }->{model} = "+m" ;#"+m"; 
			}
		}
		elsif ($p->isFather){
			if  ($hpatients->{ $p->id }->{gt} eq "-"){
				$hpatients->{ $p->id }->{model} = "-f" ;
			}
			else {
				$zmodel = "father";
				$hpatients->{ $p->id }->{model} = "+f" ;#"+m"; 
			}
		}
		else {
			if ($hpatients->{ $p->id }->{gt} eq '-') { $hpatients->{ $p->id }->{model} = '-'; }
			else {					
				$hpatients->{ $p->id }->{model} = $l->[5];
				$zmodel = "?" unless $zmodel;
				unless ( $hpatients->{ $p->id }->{model}) {$hpatients->{ $p->id }->{model} = $zmodel;}
				if ( $hpatients->{ $p->id }->{model} eq " ") {$hpatients->{ $p->id }->{model} = $zmodel;}
			}
		}
	}
	return $hpatients;		
} 

sub check_is_hgmd_dm_for_gene {
	my ($self,$hvariation,$project,$gene) = @_;
	return $hvariation->{value}->{dm_for_this_gene} if (exists $hvariation->{value}->{dm_for_this_gene});
	
	if ($hvariation->{value}->{dm}) {
		my $chr = $project->getChromosome($hvariation->{value}->{chromosome});
		my $hgmd_id = $hvariation->{value}->{hgmd_id};
		
		my $g = $project->newGene($gene->{id});
		if ($chr->is_hgmd_DM_for_gene($hgmd_id, $g)) {
			$hvariation->{value}->{dm_for_this_gene} = 1;
			return 1;
		}
		else {
			$hvariation->{value}->{dm_for_this_gene} = undef;
			$hvariation->{value}->{dm} = undef;
			$hvariation->{value}->{hgmd} = '';
			$hvariation->{html}->{hgmd} = '';
			return undef;
		}
	}
	return;
}

sub check_is_clinvar_pathogenic_for_gene  {
	my ($self,$hvariation,$project,$gene) = @_;
	return $hvariation->{value}->{clinvar_pathogenic_for_this_gene} if (exists $hvariation->{value}->{clinvar_pathogenic_for_this_gene});
	if ($hvariation->{value}->{clinvar_pathogenic}) {
		#my $chr = $project->getChromosome($hvariation->{value}->{chromosome});
		#my $clinvar_id = $hvariation->{value}->{clinvar_id};
		#my $g = $project->newGene($gene->{id});
		#if ($chr->is_clinvar_pathogenic_for_gene($clinvar_id, $g)) {
			$hvariation->{value}->{clinvar_pathogenic_for_this_gene} = 1;
			return 1;
		#}
		#else {
		#	$hvariation->{value}->{clinvar_pathogenic_for_this_gene} = undef;
		#	$hvariation->{value}->{clinvar_pathogenic} = undef;
		#	$hvariation->{value}->{clinvar} = '';
		#	$hvariation->{html}->{clinvar} = '';
		#}
	}
	return;
}


has id => (
	is		=> 'rw',
	
);

has name => (
	is		=> 'rw',
	
);

has start => (
	is		=> 'rw',
	
);

has length => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my ($self) = @_;
		return abs($self->start- $self->end) +1;
	}
		
);
has end => (
	is		=> 'rw',
);

has chromosome => (
	is		=> 'rw',

);
has allele => (
	is		=> 'rw',

);
has ref_allele => (
	is		=> 'rw',

);


has locus => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my ($self) = @_;
		return $self->chromosome.":".$self->start."-".$self->end;
	}
		
);

has isSrPr => (
	is		=> 'rw',
	default => sub {
		undef;
	}
);
has isCnv => (
	is		=> 'rw',
	default => sub {
		undef;
	}
);

has isJunction => (
	is		=> 'rw',
);

has patients_calling => (
	is		=> 'rw',
	default => sub {
		return {}
	}
);
has type  => (
	is		=> 'rw',
);

has other_genes  => (
	is		=> 'rw',
);

has rocksdb_id => (
	is		=> 'rw',
);
#GNOMAD 


has gnomad_id => (
	is		=> 'rw',
);

has gnomad_ac => (
	is		=> 'rw',
);

has gnomad_an => (
	is		=> 'rw',
);
has gnomad_ho => (
	is		=> 'rw',
);
has gnomad_ho_male => (
	is		=> 'rw',
);

has gnomad_min_pop => (
	is		=> 'rw',
);
has gnomad_min_pop_name => (
	is		=> 'rw',
);
has gnomad_max_pop_name => (
	is		=> 'rw',
);
has gnomad_max_pop => (
	is		=> 'rw',
);
#################################
#dejavu 
#################################
has dejavu_other_projects => (
	is		=> 'rw',
);
has dejavu_other_patients => (
	is		=> 'rw',
);
has dejavu_other_patients_ho => (
	is		=> 'rw',
);
has dejavu_similar_projects => (
	is		=> 'rw',
);
has dejavu_similar_patients => (
	is		=> 'rw',
);
has dejavu_similar_patients_ho => (
	is		=> 'rw',
);

has dejavu_this_run_patients => (
	is		=> 'rw',
);


has transcripts => (
	is		=> 'rw',
);


has gene => (
	is		=> 'rw',
	
);
#####################################
# SCORE
####################################
has revel => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		return "-";
	}
	
);

has ada => (
	is		=> 'rw',
	default => sub {
		return "-";
	}
	
);
has rf => (
	is		=> 'rw',
	default => sub {
		return "-";
	}
	
);

has cadd => (
	is		=> 'rw',
	default => sub {
		return "-";
	}
	
);
has spliceAI => (
	is		=> 'rw',
	
);
has spliceAI_cat => (
	is		=> 'rw',
	
);
###########################
#Validation
###########################

has hgmd_dm_for_gene => (
	is		=> 'rw',
);
has clinvar_pathogenic_for_gene => (
	is		=> 'rw',
);
has dm => (
	is		=> 'rw',
);
has clinvar_pathogenic => (
	is		=> 'rw',
);
has hgmd_id => (
	is		=> 'rw',
);
has hgmd => (
	is		=> 'rw',
);
has clinvar => (
	is		=> 'rw',
);
has clinvar_id => (
	is		=> 'rw',
);
has hgmd_phenotype => (
	is		=> 'rw',
);
has text_caller => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		return [];
	}
);

has other_gene => (
	is		=> 'rw',
);
has cnv_details_genes => (
	is		=> 'rw',
);


has hgmd_value => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my ($self) = @_;
		return -1 unless $self->hgmd();
		my $nb = 1;
		return 4 if $self->hgmd eq 'DM';
		return 3 if $self->hgmd eq 'DM?';	
		return 2;
	}
);

has clinvar_value => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
	my ($self) = @_;
	return -1 unless $self->clinvar();
	if ($self->clinvar =~ /pathogenic/i){
		return 4 if $self->clinvar =~ /likely/i;
		return 5;
	};
	return 1;
	}		
);
has isDude => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
	
	return undef;
	}		
);

has is_variant_forced_viewing  => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		return undef;
	}		
);

has isMantaImprecise  => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		return undef;
	}		
);

has variants_same_position  => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		return undef;
	}		
);
has log2_ratio  => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		return 0;
	}		
);
sub setCnvValues {
	my ($self,$chr,$patient,$variant)  =@_;
	
	foreach my $p (@{$patient->getFamily()->getMembers}){
		$self->patients_calling->{ $p->id }->{gt} = $variant->getSequencingGenotype($p);
		$self->patients_calling->{ $p->id }->{pr} = $variant->pr($p);
		$self->patients_calling->{ $p->id }->{sr} = $variant->sr($p);
		$self->patients_calling->{ $p->id }->{norm_depth} = $variant->getNormDP($p);
		$self->patients_calling->{ $p->id }->{dude_score} = $variant->getCNVDude($p);
		$self->patients_calling->{ $p->id }->{model}  = $variant->getTransmissionModelType($p->getFamily(),$p);
	}
}

sub setDudeValues {
	my ($self,$chr,$patient,$variant)  =@_;
	my $primers = $patient->getProject->getPrimersByPosition($chr,$self->start,$self->end);
	foreach my $primer (@$primers){
		foreach my $p (@{$patient->getFamily()->getMembers}){
			$self->patients_calling->{$p->id}->{norm_depth} = $p->cnv_region_ratio_norm($chr->name,$primer->start,$primer->end);
			$self->patients_calling->{$p->id}->{dude_score} = $primer->cnv_score_log2($p);
			$self->patients_calling->{$p->id}->{model} = $variant->getTransmissionModelType($p->getFamily(),$p);
			$self->patients_calling->{$p->id}->{gt} = $variant->getSequencingGenotype($p);
		}
	}
}
sub init {
	my ($self) =@_;
	my $hh;
	$hh->{id} = undef;
	$hh->{gt} = undef;
	$hh->{pc} = undef;
	$hh->{dp} = undef;
	$hh->{model} =undef;
	$hh->{pr} = undef;
	$hh->{sr} =undef;
	$hh->{norm_depth} = undef;
	$hh->{dude_score} = undef;
	$hh->{norm_depth_before} = undef;			
	$hh->{norm_after_after} = undef;			
}

sub set_patient_cache {
		my ($self,$vh,$p) =@_;
		
		my $hh = $self->init();
		$hh->{id} = $p->id;
		$hh->{pr} = undef;
		$hh->{sr} = undef;
		unless ($vh->existsPatient($p)){
			$hh->{norm_depth} = $vh->getNormDP($p);
			$hh->{dp} = $vh->getDP($p);
			if ($vh->isSrPr or $vh->isDude){
				$hh->{norm_depth_before} =  $vh->getNormDPBefore($p);
				$hh->{norm_depth_after} = $vh->getNormDPAfter($p);
			}
			return $hh;
		}
		if ($vh->isJunction()) {
			$hh->{name} = $p->name;
			$hh->{isRI} = $vh->isRI($p);
			$hh->{isSE} = $vh->isSE($p);
			$hh->{is_sj} = $vh->is_sj($p);
			$hh->{is_ri_aval} = $vh->is_ri_aval($p);
			$hh->{is_ri_amont} = $vh->is_ri_amont($p);
			$hh->{isRegtools} = $vh->isRegtools($p);
			$hh->{get_nb_new_count} = $vh->get_nb_new_count($p);
			$hh->{get_canonic_count} = $vh->get_canonic_count($p);
			$hh->{get_dp_count} = $vh->get_dp_count($p);
			$hh->{ratio} = $vh->ratio($p);
			$hh->{get_percent_new_count} = $vh->get_percent_new_count($p);
			$hh->{multiple_align_count} = $vh->multiple_align_count($p);
			$hh->{is_sj} = $vh->is_sj($p);
		}
		elsif ($vh->isSrPr){
			$hh->{gt} = $vh->getSequencingGenotype($p);
			$hh->{norm_depth} = $vh->getNormDP($p);
			$hh->{dude_score} = $vh->getCNVDude($p);
			$hh->{pr} = $vh->pr($p);
			$hh->{sr} = $vh->sr($p);
			$hh->{norm_depth_before} =  $vh->getNormDPBefore($p);
			$hh->{norm_depth_after} = $vh->getNormDPAfter($p);
			#$hh->{log2_ratio} =  $vh->getLog2Ratio($p);
			
		}
		elsif ($vh->isDude){
			my $primers = $p->getProject->getPrimersByPosition($vh->getChromosome,$vh->start,$vh->end);
			my $primer = $primers->[0]; 
			$hh->{norm_depth} = $p->cnv_region_ratio_norm($vh->getChromosome,$primer->start,$primer->end);
			$hh->{dude_score} = $primer->cnv_score_log2($p);
			$hh->{gt} = $vh->getSequencingGenotype($p);
			$hh->{norm_depth_before} =  $vh->getNormDPBefore($p);
			$hh->{norm_after_after} = $vh->getNormDPAfter($p);
		
		}
		else {
			$hh->{name} = $p->name;
			$hh->{id} = $p->id;
			$hh->{gt} = $vh->getSequencingGenotype($p);
			$hh->{pc} = $vh->getRatio($p);
			$hh->{dp} = $vh->getDP($p);
			#$hh->{model} = "-";#$vh->getTransmissionModelType($p->getFamily(),$p);
			
		}
		$hh->{array_text_calling} = $vh->array_sequencing_text($p) if not $vh->isJunction();
		return $hh;
}

sub set_patient_cache_1 {
	my ($self,$vh,$patient) =@_;
	my $apatients;
	my $hh = $self->init();
	
	
	#init
		
		foreach my $p (@{$patient->getFamily()->getMembers}){
			my $hh = $self->init();
			$hh->{id} = $p->id;
			$hh->{pr} = undef;
			$hh->{sr} = undef;
			
			if ($vh->isSrPr){
			$hh->{gt} = $vh->getSequencingGenotype($p);
			$hh->{norm_depth} = $vh->getNormDP($p);
			$hh->{dude_score} = $vh->getCNVDude($p);
			$hh->{pr} = $vh->pr($p);
			$hh->{sr} = $vh->sr($p);
			
		}
		elsif ($vh->isDude){
			my $primers = $patient->getProject->getPrimersByPosition($vh->getChromosome,$vh->start,$vh->end);
			my $primer = $primers->[0]; 
			$hh->{norm_depth} = $p->cnv_region_ratio_norm($vh->getChromosome,$primer->start,$primer->end);
			$hh->{dude_score} = $primer->cnv_score_log2($p);
			$hh->{gt} = $vh->getSequencingGenotype($p);
		}
		else {
				$hh->{name} = $p->name;
				$hh->{id} = $p->id;
				$hh->{gt} = $vh->getSequencingGenotype($p);
				$hh->{pc} = $vh->getRatio($p);
				$hh->{dp} = $vh->getDP($p);
				#$hh->{model} = "-";#$vh->getTransmissionModelType($p->getFamily(),$p);
				
			}
			$hh->{array_text_calling} = $vh->array_sequencing_text($p);
			push(@$apatients,$hh);
		}
		return $apatients;
}
sub set_patient {
	my ($self,$vh,$patient) =@_;
	my $hpatients;
		if ($vh->isJunction()) {
			
		}
		elsif ($vh->isDude){
			$self->setCnvValues($vh->getChromosome,$patient,$vh);
		}
		elsif ($vh->isCnv){
			$self->setCnvValues($vh->getChromosome,$patient,$vh);
		}
		else {
			foreach my $p (@{$patient->getFamily()->getMembers}) {
				$hpatients->{ $p->id }->{gt} = $vh->getSequencingGenotype($p);
				$hpatients->{ $p->id }->{pc} = $vh->getRatio($p);
				$hpatients->{ $p->id }->{dp} = $vh->getDP($p);
				$hpatients->{ $p->id }->{model} = "-";#$vh->getTransmissionModelType($p->getFamily(),$p);
			}
		}
		return $hpatients;
}


sub return_specific_value {
	my ($self,$value) =@_;
	unless (defined $value){
		return -99;
	}
	if ($value eq "-"){
		return -99;
	}
	if ($value == -1 ){
		return -99;
	}
	return $value;
}


sub set_gene {
		my ($self,$v,$gene) =@_;
		my  $transcripts = $v->getTranscripts();
		my $all_transcripts = [];
		foreach my $tr1 (sort { ( ($b->isMane <=> $a->isMane) or $self->impact_sorted->{$v->effectImpact($b)} <=>  $self->impact_sorted->{$v->effectImpact($a)})  or ($a->appris_level <=> $b->appris_level)} @$transcripts) {
				next if $tr1->getGene->id ne $gene->id;
				my $htr = {};
	
				my $enst = $tr1->name;
				$htr->{name} = $tr1->name();
				$htr->{nm} = $tr1->external_name;
				$htr->{ccds} = $tr1->ccds_name;
				$htr->{ccds} = "-" if ($tr1->ccds_name && length($tr1->ccds_name) <3);
				$htr->{appris} = $tr1->appris_type;
				$htr->{main}=  0 ;
		 		$htr->{main}  = 1 if $tr1->isMain();
		 		$htr->{mane}  = 0 ;
		 		$htr->{mane}  = 1 if $tr1->isMane();
				my @coding_infos = ("sift","polyphen","prot","codons","codons_AA","impact_score_text","impact_score","consequence");
			
				#initialize
				foreach my $c (@coding_infos){
					$htr->{$c} = "-";
				}
				my $sc = $v->effectImpact($tr1);
				$htr->{impact_score_text} = $sc;
				$htr->{impact_score} = $self->impact_sorted->{$sc};	
				$htr->{consequence} = $v->variationTypeInterface($tr1);;	
				$htr->{nomenclature} =  $v->getNomenclature($tr1);
				$htr->{nomenclature} = "-" unless $htr->{nomenclature} ;
				$htr->{exon} = $tr1->findExonNumber($v->start, $v->end);
				$htr->{exon} = $tr1->findNearestExon($v->start, $v->end) if $htr->{exon} == -1;
				$htr->{start} = "";
				$htr->{end} = "";
				$htr->{cadd} = $self->return_specific_value($v->cadd_score);
				$htr->{dbscsnv} = $self->return_specific_value($v->dbscsnv_rf);
				$htr->{dbscsnv} = $self->return_specific_value($v->dbscsnv_ada) if $self->return_specific_value($v->dbscsnv_ada) > $htr->{dbscsnv} ;
				$htr->{revel} = $self->return_specific_value($v->revel_score($tr1));
				
				$htr->{spliceAI} = $self->return_specific_value($v->max_spliceAI_score($gene));
				 
				$htr->{spliceAI_cat} = $v->max_spliceAI_categorie($gene);
				
				$htr->{sift} = -99;
				$htr->{polyphen} = -99;
				$htr->{alphamissense} = -99;
				
				my $r1 = $tr1->exons_introns_tree()->fetch($v->start,$v->start+1);
				if (@$r1){
					if ($r1->[0]->{type} eq "intron"){
						$htr->{start} = "intron".$r1->[0]->{nb};
					}
					else {
						$htr->{start} = $r1->[0]->{name};
					}
				}
				 $r1 = $tr1->exons_introns_tree()->fetch($v->end,$v->end+1);
				 if (@$r1){
					
					if ($r1->[0]->{type} eq "intron"){
						$htr->{end} = "intron".$r1->[0]->{nb};
					}
					else {
						$htr->{end} = $r1->[0]->{name};
					}
				}
				 
				if ($v->isCoding($tr1)) {
					my $prot = $tr1->getProtein();
					if ($v->isLargeDeletion or $v->isLargeDuplication){
						foreach my $c (@coding_infos){
							$htr->{$c} = "-";
						}
						$htr->{impact_score_text} ="high";
						$htr->{impact_score} = 4;
						$htr->{consequence} = $v->variationTypeInterface($tr1);;	
						
					}
					else {
						$htr->{prot} = $v->getProteinPosition($prot);
						$htr->{codons} = $v->getCodons($tr1);
						$htr->{codons_AA} = $v->protein_nomenclature($prot);
						$htr->{sift} =  $self->return_specific_value($v->siftScore($tr1));
						$htr->{polyphen} =  $self->return_specific_value($v->polyphenScore($tr1));
						
						$htr->{alphamissense} = $self->return_specific_value($v->alphamissense($tr1));
						
						
					}
				
				}
					push(@$all_transcripts,$htr);
		}
		my $score;
		if (exists $v->genes_pathogenic_DM->{$gene->id} && $v->genes_pathogenic_DM->{$gene->id}->{DM} ){
				$score->{dm} =$v->isDM;
				$score->{hgmd_phenotype} = $v->hgmd_phenotype;
			}
		
		
		
		if (exists $v->genes_pathogenic_DM->{$gene->id} && $v->genes_pathogenic_DM->{$gene->id}->{pathogenic} ){
				$score->{clinvar_pathogenic} = $v->isClinvarPathogenic;
		}
		$score->{spliceAI} = $v->max_spliceAI_score($gene);
		
		$score->{spliceAI_cat} = $v->max_spliceAI_categorie($gene);
		
		
		return {sc=>$score,tr=>$all_transcripts};
}


sub set_gene_junction {
		my ($self,$j,$gene) =@_;
		my  $transcripts = $gene->getTranscripts();
		my $all_transcripts = [];
		foreach my $tr1 (sort { ( ($b->isMane <=> $a->isMane))  or ($a->appris_level <=> $b->appris_level)} @$transcripts) {
				next if $tr1->getGene->id ne $gene->id;
				my $htr = {};
	
				my $enst = $tr1->name;
				$htr->{name} = $tr1->name();
				$htr->{nm} = $tr1->external_name;
				$htr->{ccds} = $tr1->ccds_name;
				$htr->{ccds} = "-" if ($tr1->ccds_name && length($tr1->ccds_name) <3);
				$htr->{appris} = $tr1->appris_type;
				$htr->{main}=  0 ;
		 		$htr->{main}  = 1 if $tr1->isMain();
		 		$htr->{mane}  = 0 ;
		 		$htr->{mane}  = 1 if $tr1->isMane();
				my @coding_infos = ("sift","polyphen","prot","codons","codons_AA","impact_score_text","impact_score","consequence");
			
				#initialize
				foreach my $c (@coding_infos){
					$htr->{$c} = "-";
				}
				$htr->{exon} = $tr1->findExonNumber($j->start, $j->end);
				$htr->{exon} = $tr1->findNearestExon($j->start, $j->end) if $htr->{exon} == -1;
				$htr->{start} = "";
				$htr->{end} = "";
				my $r1 = $tr1->exons_introns_tree()->fetch($j->start,$j->start+1);
				if (@$r1){
					if ($r1->[0]->{type} eq "intron"){
						$htr->{start} = "intron".$r1->[0]->{nb};
					}
					else {
						$htr->{start} = $r1->[0]->{name};
					}
				}
				 $r1 = $tr1->exons_introns_tree()->fetch($j->end,$j->end+1);
				 if (@$r1){
					
					if ($r1->[0]->{type} eq "intron"){
						$htr->{end} = "intron".$r1->[0]->{nb};
					}
					else {
						$htr->{end} = $r1->[0]->{name};
					}
				}
				 
				push(@$all_transcripts,$htr);
		}
		my $score;
		return {sc=>$score,tr=>$all_transcripts};
}


has get_genes_transcripts_details_dup_del => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		my $h;
		my $start1 = $self->start();
		my $end1 = $self->start() + $self->length();
		my $intspan_v = $self->getGenomicSpan();
		foreach my $g (@{$self->getGenes()}) {
			foreach my $t (@{$g->getTranscripts()}) {
				my $intspan_t = $t->getGenomicSpan();
				my $inter1 = $intspan_v->intersection( $intspan_t );
				next if $inter1->is_empty();
				my $t_id = $t->id();
				$h->{$g->id()}->{$t_id}->{nm} = 'NM';
				$h->{$g->id()}->{$t_id}->{ccds} = $t->ccds_name();
			
				#$h->{$g->id()}->{$t_id}->{appris} = $t->{tag};
				foreach my $e (@{$t->getExons()}) {
					my $intspan_e = $e->getGenomicSpan();
					my $inter2 = $intspan_v->intersection( $intspan_e );
					my $e_id = $e->id();
					$e_id =~ s/$t_id//;
					next if ($inter2->is_empty());
					my @lTmp = split('-', $inter2->as_string());
					my $overlap = $lTmp[-1] - $lTmp[0] + 1;
					next if ($overlap < 1);
					$h->{$g->id()}->{$t_id}->{positions}->{$e->start()} = $e_id; 
					my $perc = sprintf("%.2f", ($overlap / $e->length) * 100);
					$h->{$g->id()}->{$t_id}->{exons_introns}->{$e_id}->{locus} = $e->getChromosome->id().':'.$e->start().'-'.$e->end();
					$h->{$g->id()}->{$t_id}->{exons_introns}->{$e_id}->{overlap} = $overlap.' nt ('.$perc.'%)';
				}
				foreach my $i (@{$t->getIntrons()}) {
					my $intspan_i = $i->getGenomicSpan();
					my $inter2 = $intspan_v->intersection( $intspan_i );
					my $i_id = $i->id();
					$i_id =~ s/$t_id//;
					next if ($inter2->is_empty());
					my @lTmp = split('-', $inter2->as_string());
					my $overlap = $lTmp[-1] - $lTmp[0] + 1;
					next if ($overlap < 1);
					$h->{$g->id()}->{$t_id}->{positions}->{$i->start()} = $i_id; 
					my $perc = sprintf("%.2f", ($overlap / $i->length) * 100);
					$h->{$g->id()}->{$t_id}->{exons_introns}->{$i_id}->{locus} = $overlap.' nt ('.$perc.'%)';
					$h->{$g->id()}->{$t_id}->{exons_introns}->{$i_id}->{overlap} = $overlap.' nt ('.$perc.'%)';
				}
			}
		}
		return $h;
	},
);
sub setParquetVariant {
	my ($self,$vh,$patient) = @_; 
	my $project = $patient->getProject();
	$self->id($vh->{variant_index});
	$self->start($vh->{variant_start});
	$self->end($vh->{variant_start});
	
	#$self->ref_allele($vh->ref_allele);	
	#$self->allele($vh->alternate_allele);
	
#		if ($self->allele && length ($self->allele) > 50 ){
#			my $s = "+";
#			$s="-" if $vh->isDeletion;
#			$self->allele($s.length($self->allele));
#		}

		#######################
		# les IOs
		#######################

		$self->gnomad_id($vh->{variant_gnomad_id});
		$self->rocksdb_id($vh->{variant_rocksdb_id});
		$self->chromosome($vh->{variant_chromosome});
		$self->name($vh->{variant_gnomad_id});
		$self->type($vh->{variant_type});
		#$self->{reference} = ref($vh);
		#######################
		#hgmd et clinvar
		#######################
		
		
		

		$self->clinvar_id($vh->{variant_clinvar_id});
		$self->clinvar($vh->{variant_clinvar_class}) ;
		
		
		################
		# Calling
		##################
		
		#$self->isJunction(0);
		#$self->isSrPr($vh->isSrPr);
		#if ($vh->isCnv){
		#	$self->isCnv(1);
		#}
		#else {
			$self->isCnv(0);
		#}

		
		#################
		# gnomad
		##################
		
		$self->gnomad_ac($vh->{variant_gnomad_ac});
		$self->gnomad_an($vh->{variant_gnomad_an});
		$self->gnomad_min_pop_name($vh->{gnomad_min_pop_name});
		$self->gnomad_min_pop($vh->{gnomad_min_pop_freq});
		$self->gnomad_max_pop_name($vh->{gnomad_max_pop_name});
		$self->gnomad_max_pop($vh->{gnomad_min_pop_freq});
		$self->gnomad_ho($vh->{variant_gnomad_ho});
		$self->gnomad_ho_male($vh->{variant_getGnomadAC_Male});
		
	
	#########
	# DEJAVU
	#########
	 
	$self->dejavu_other_projects($vh->{variant_other_projects});
	$self->dejavu_other_patients($vh->{variant_other_patients});
	$self->dejavu_other_patients_ho($vh->{variant_other_patients_ho});
	$self->dejavu_similar_projects( $vh->{variant_similar_projects});
	$self->dejavu_similar_patients($vh->{variant_similar_patients});
	$self->dejavu_similar_patients_ho($vh->{variant_similar_patients_ho});
	$self->dejavu_this_run_patients($vh->{in_this_run_patients});# = '-';
	
	$self->text_caller([]);

	foreach my $p (@{$patient->getFamily()->getMembers}){
		my $suffix = "patient_".$patient->id;
		next unless exists $vh->{$suffix."_alt"};
	#	next if $vh->{$suffix."_alt"} == -1;
		
		$self->{patients_calling}->{$patient->id}->{model} = $vh->{$suffix."_transmission"} ;
		$self->{patients_calling}->{$patient->id}->{gt} = $vh->{$suffix."_type"} ;
		$self->{patients_calling}->{$patient->id}->{sr} = undef ;
		$self->{patients_calling}->{$patient->id}->{pr} = undef ;
		$self->{patients_calling}->{$patient->id}->{id} =  $patient->id;
		$self->{patients_calling}->{$patient->id}->{dp} = $vh->{$suffix."_alt"} +$vh->{$suffix."_ref"} ;
		$self->{patients_calling}->{$patient->id}->{pc} = $vh->{$suffix."_ratio"} ;
		$self->{patients_calling}->{$patient->id}->{array_text_calling} = [];	

	}	
	
	#push(@column_patient,"patient_".$c."_ref");
	#push(@column_patient,"patient_".$c."_alt");
	#push(@column_patient,"patient_".$c."_ratio");
	#push(@column_patient,"patient_".$c."_type");
	#push(@column_patient,"patient_".$c."_transmission");
	
}
sub setTranscriptFromParquet {
	my ($self,$vh,$project) =@_;
	$self->$self->transcripts($vh);
}

sub setLmdbVariant {
	my ($self,$vh,$project,$gene,$patient) = @_; 
		if ($project){
			$vh->{project} = $project;
			$vh->{buffer} = $project->buffer;
		}
		#$self->gene($gene);
		$self->id($vh->id);
		$self->start($vh->start);
		$self->end($vh->end);
		$self->ref_allele($vh->ref_allele);	
		$self->allele($vh->alternate_allele);
		if ($self->allele && length ($self->allele) > 50 ){
			my $s = "+";
			$s="-" if $vh->isDeletion;
			$self->allele($s.length($self->allele));
		}
		$self->gnomad_id($vh->gnomad_id);
		$self->rocksdb_id($vh->rocksdb_id);
		$self->chromosome($vh->getChromosome->name);
		$self->name($vh->name);
		$self->type($vh->type);
		$self->{reference} = ref($vh);
		$self->gnomad_id($vh->gnomad_id);
		$self->rocksdb_id($vh->rocksdb_id);
		#######################
		#hgmd et clinvar
		#######################
		
		
		

		$self->clinvar_id($vh->clinvar_id);
		$self->clinvar($vh->clinvar_class) ;
		
		
		################
		# Calling
		##################
		$self->isJunction(0);
		$self->isSrPr($vh->isSrPr);
		if ($vh->isCnv){
			$self->isCnv(1);
		}
		else {
			$self->isCnv(0);
		}
	
		
		if($patient) {
			my $hpatients = $self->set_patient($patient);
			$self->patients_calling($hpatients);
		};
		if($gene){
			my $h = $self->set_gene($gene);
			$self->transcripts($h->{tr});
			foreach my $k (keys %{$h->{sc}}){
				$self->{$k} = $h->{sc}->{$k};
			}
		}
		
		#################
		# gnomad
		##################
		
		$self->gnomad_ac($vh->getGnomadAC);
		$self->gnomad_an($vh->getGnomadAN);
		$self->gnomad_min_pop_name($vh->min_pop_name);
		$self->gnomad_min_pop($vh->min_pop_freq);
		$self->gnomad_max_pop_name($vh->max_pop_name);
		$self->gnomad_max_pop($vh->max_pop_freq);
		$self->gnomad_ho($vh->getGnomadHO);
		$self->gnomad_ho_male($vh->getGnomadAC_Male);
	
	#########
	# DEJAVU
	#########
	 
	$self->dejavu_other_projects($vh->other_projects());
	$self->dejavu_other_patients($vh->other_patients());
	$self->dejavu_other_patients_ho($vh->other_patients_ho());
	$self->dejavu_similar_projects( $vh->similar_projects());
	$self->dejavu_similar_patients($vh->similar_patients());
	$self->dejavu_similar_patients_ho($vh->similar_patients_ho());
	$self->dejavu_this_run_patients($vh->in_this_run_patients);# = '-';
	$self->text_caller([]);
	$self->hgmd($vh->hgmd_class);
		$self->hgmd_id($vh->hgmd_id);
	########
	# SCORE
	########
#	 $self->revel($vh->revel_score());
#	 $self->rf($vh->dbscsnv_rf());	
#	 $self->ada($vh->dbscsnv_rf());
#	 $self->cadd($vh->cadd_score);
	 my $ag;
	 my $genes = $vh->getGenes;
   	$self->{other_genes} = {}; 
	 if (scalar (@$genes) >1 ){
	 	foreach my $g (@{$genes}){
	 		$self->{other_genes}->{$g->id} = $g->external_name;
		}
	 }
	
}

sub setLmdbJunction {
	my ($self,$vh,$project,$gene,$patient) = @_; 
		if ($project){
			$vh->{project} = $project;
			$vh->{buffer} = $project->buffer;
		}
		$self->id($vh->id);
		$self->start($vh->start);
		$self->end($vh->end);
		$self->rocksdb_id($vh->rocksdb_id);
		$self->chromosome($vh->getChromosome->name);
		$self->name($vh->name);
		$self->type($vh->type);
		$self->{reference} = ref($vh);
		$self->rocksdb_id($vh->rocksdb_id);
		$self->isJunction(1);
		
		if($patient) {
			my $hpatients = $self->set_patient_cache($vh, $patient);
			$self->patients_calling($hpatients);
		};
		if($gene){
			my $h = $self->set_gene_junction($gene);
			$self->transcripts($h->{tr});
			foreach my $k (keys %{$h->{sc}}){
				$self->{$k} = $h->{sc}->{$k};
			}
		}
		
#	 my $ag;
	 my $genes = $vh->getGenes;
	$self->{other_genes} = {}; 
	 if (scalar (@$genes) >1 ){
	 	foreach my $g (@{$genes}){
	 		$self->{other_genes}->{$g->id} = $g->external_name;
		}
	 }
	
}



sub update_clinvar {
	my ($self,$project,$gene) = @_;
		my $vh = $project->_newVariant($self->id);
		 	
	
		$self->clinvar_id($vh->clinvar_id);
		return unless $vh->clinvar_id;
		$self->clinvar($vh->clinvar_class) ;
		if (exists $vh->genes_pathogenic_DM->{$gene->{id}} && $vh->genes_pathogenic_DM->{$gene->{id}}->{pathogenic} ){
				$self->clinvar_pathogenic($vh->isClinvarPathogenic);
		}
	
}

sub setOldVariant {
	my ($self,$vh,$project,$patient,$gene,$update,$debug) = @_; 
	
		$self->gene($gene);
		$self->id($vh->{value}->{id});
		$self->start($vh->{value}->{start});
		$self->end($vh->{value}->{end});
		$self->ref_allele($vh->{value}->{ref_allele});	
		$self->allele($vh->{value}->{allele});	
		$self->gnomad_id($vh->{value}->{gnomad_id});
		$self->chromosome($vh->{value}->{chromosome});
		$self->isCnv($vh->{value}->{is_cnv});
		if ($self->isCnv) {
			foreach my $caller (@{$vh->{value}->{caller}}) {
				$self->isDude(1) if ($caller =~ /dud/);
			}
			if (exists $vh->{value}->{manta} and exists $vh->{value}->{manta}->{is_imprecise}) {
				$self->isMantaImprecise(1) if ($vh->{value}->{manta}->{is_imprecise});
			}
		}
		$self->type($vh->{value}->{type});
		$self->name($vh->{value}->{var_name});
		
		$self->check_is_hgmd_dm_for_gene($vh,$project,$gene);
		$self->hgmd_dm_for_gene($vh->{value}->{dm_for_this_gene});
		$self->check_is_clinvar_pathogenic_for_gene($vh,$project,$gene);
		
		$self->clinvar_pathogenic_for_gene($vh->{value}->{clinvar_pathogenic_for_this_gene});
		
		$self->dm($vh->{value}->{dm});
		$self->clinvar_pathogenic($vh->{value}->{clinvar_pathogenic});
		$self->hgmd_id($vh->{value}->{hgmd_id});
		$self->clinvar_id($vh->{value}->{clinvar_id});
		$self->hgmd_phenotype($vh->{value}->{hgmd_phenotype});
		if ($vh->{value}->{hgmd_id}){
		my ($v1,$v2) = split(":",$vh->{value}->{hgmd});
		$self->hgmd($v1);
		}
		if ($vh->{value}->{clinvar_id}){
				$self->clinvar($vh->{value}->{clinvar});
		}
		if (exists $vh->{value}->{cnv_details_genes}) {
			$self->cnv_details_genes($vh->{value}->{cnv_details_genes});
		}
		$self->update_clinvar($project,$gene) if $update; 
		unless ($self->name)
		{
			my ($a,$b,$c,$d) = split("-",$self->gnomad_id);
			my $l1 = length($c) -1;
			my $l2 = length($d) -1;
			if ($l1>15){  
				$self->name($a."- ".$b."-del-".$l1);
			}
			elsif ($l2>15){
				$self->name($a."-".$b."-ins-".$l2);
			}
			else {
				$self->name($self->gnomad_id);
			}
		}
				
		##
		# gnomad
		######
		$self->parseGnomadTable($vh->{html}->{gnomad},$vh);


	
	#########
	# DEJAVU
	#########
	if ($project->isDiagnostic){
	my $h = $project->getDejaVuInfosForDiag($self->id);
	
	$self->dejavu_other_projects($h->{other_projects});
	$self->dejavu_other_patients($h->{other_patients});
	$self->dejavu_other_patients_ho($h->{other_patients_ho});
	$self->dejavu_similar_projects( $h->{similar_projects});
	$self->dejavu_similar_patients($h->{similar_patients});
	$self->dejavu_similar_patients_ho($h->{similar_patients_ho});
	$self->dejavu_this_run_patients($h->{in_this_run_patients});# = '-';
	}
	else {
	$self->parseDejaVuTable($vh->{html}->{deja_vu},$vh);
	$self->dejavu_this_run_patients($vh->{value}->{this_run_patients});# = '-';
	}
	#########
	# caller
	#########
	$self->text_caller($vh->{value}->{ngs});
	
	
	my @atr;
	my $tgene = delete $vh->{genes}->{$gene->{id}};
	foreach my $l (@$tgene) {
		
		delete $l->{html};
		push(@atr,$l->{value}) ;
		
	}
	
	my @k = keys %{$vh->{genes}};
	my @g;
	delete $project->{liteAnnotations};
	foreach my $gid (@k){
		my ($g,$c) = split("_",$gid);
		my $gene = $project->newGene($gid);
		push(@g,$gene->external_name);	
	}
	
   $self->other_genes(\@g); 
   
   
  # delete $vh->{genes};
   
	$self->transcripts(\@atr);
	if ($self->isDude) {
		my $vh = $project->_newVariant($self->id);
		$self->setDudeValues($vh->getChromosome,$patient,$vh);
	}
	elsif ($self->isCnv){
		my $vh = $project->_newVariant($self->id); 
		$self->setCnvValues($vh->getChromosome,$patient,$vh);
	}
	else {
		$self->patients_calling($self->parseTrioTableVariant($vh->{html}->{trio},$project));
	}
	
	 $self->revel($self->transcripts->[0]->{revel}) if $self->transcripts->[0]->{revel};
	 $self->rf($self->transcripts->[0]->{rf}) if $self->transcripts->[0]->{rf};
	 $self->ada($self->transcripts->[0]->{ada}) if $self->transcripts->[0]->{ada};
	 $self->cadd($self->transcripts->[0]->{cadd}) if $self->transcripts->[0]->{cadd}; 
	 if (exists $vh->{value}->{spliceAI}->{$gene->{id}}){
	 my	($v1,$v2) = split(":",$vh->{value}->{spliceAI}->{$gene->{id}});
	 	$self->spliceAI_cat($v1);
	 	$self->spliceAI($v2);
	 }
	 my @lVar_same_pos = @{$self->get_variants_same_position($project->getChromosome($self->chromosome()), $vh->{vector_id}, $patient)};
 	 if (scalar(@lVar_same_pos) > 0) {
 		foreach my $other_var (@lVar_same_pos) {
 			next if ($other_var->isJunction());
 			push(@{$self->{patients_calling}->{$patient->id()}->{var_same_pos}->{name}}, $other_var->name());
			my $heho_other = "ho" if $other_var->isHomozygote($patient);
			$heho_other = "he" if $other_var->isHeterozygote($patient);
			$heho_other .= '('.$other_var->ref_allele().'/'.$other_var->var_allele().')';
 			push(@{$self->{patients_calling}->{$patient->id()}->{var_same_pos}->{gt}}, $heho_other);
 		}
 	}
}

sub get_variants_same_position {
	my ($self, $chr, $vector_id, $init_patient) = @_;
	my @lVar_same_pos;
	my $check_before_done = 0;
	my $check_after_done = 0;
	my $i = $vector_id;
	my $j = $vector_id;
	my $var = $chr->getVarObject($vector_id);
	$check_before_done = 1 if ($vector_id == 0);
	$check_after_done = 1 if ($vector_id == ($chr->getVariantsVector->Size() -1) );
	while ($check_before_done == 0) {
		$i--;
		my $var_before = $chr->getVarObject($i);
		if ($var_before->start() == $var->start()) {
			foreach my $patient (@{$var_before->getPatients()}) {
				next if ($patient->name() ne $init_patient->name());
				push(@lVar_same_pos, $var_before);
			}
		}
		else { $check_before_done = 1; }
		$check_before_done = 1 if ($i == 0);
	}
	while ($check_after_done == 0) {
		$j++;
		my $var_after = $chr->getVarObject($j);
		if ($var_after->start() == $var->start()) {
			foreach my $patient (@{$var_after->getPatients()}) {
				next if ($patient->name() ne $init_patient->name());
				push(@lVar_same_pos, $var_after);
			}
		}
		else { $check_after_done = 1; }
		$check_after_done = 1 if ($vector_id == ($chr->getVariantsVector->Size() - 1));
	}
	return \@lVar_same_pos;
}


1;

