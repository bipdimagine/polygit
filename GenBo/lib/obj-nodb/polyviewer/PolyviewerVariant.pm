package PolyviewerVariant;
use strict;
use FindBin qw($Bin);
use Moose;
use Data::Dumper;



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
			
		#	warn $self->gnomad_an() if $self->name eq "2-179306434-G-GACTCTAGCCTGC";
		#	die() if $self->name eq "2-179306434-G-GACTCTAGCCTGC";
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
has gnomad_id => (
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

has isCnv => (
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

sub setLmdbVariant {
	my ($self,$vh,$project,$gene,$patient) = @_; 
		$vh->{project} = $project;
		$vh->{buffer} = $project->buffer;
		$self->gene($gene);
		$self->id($vh->id);
		$self->start($vh->start);
		$self->end($vh->end);
		$self->ref_allele($vh->ref_allele);	
		$self->allele($vh->alternate_allele);	
		$self->gnomad_id($vh->gnomad_id);
		$self->chromosome($vh->name);
		$self->name($vh->name);
		$self->type($vh->type);
		#######################
		#hgmd et clinvar
		#######################
		
		
		$self->hgmd($vh->hgmd_class);
		$self->hgmd_id($vh->hgmd_id);
		if (exists $vh->genes_pathogenic_DM->{$gene->{id}} && $vh->genes_pathogenic_DM->{$gene->{id}}->{DM} ){
				$self->dm($vh->isDM);
				$self->hgmd_phenotype($vh->hgmd_phenotype);
				#$self->hgmd($vh->hgmd_class);
				#$self->hgmd_id($vh->hgmd_id);
			}
		
		
		$self->clinvar_id($vh->clinvar_id);
		$self->clinvar($vh->clinvar_class);
		if (exists $vh->genes_pathogenic_DM->{$gene->{id}} && $vh->genes_pathogenic_DM->{$gene->{id}}->{pathogenic} ){
				
			#	$self->check_is_clinvar_pathogenic_for_gene($vh,$project,$gene);
				#$self->clinvar_pathogenic_for_gene($vh->{value}->{clinvar_pathogenic_for_this_gene});
				$self->clinvar_pathogenic($vh->isClinvarPathogenic);
		}
		
		
		################
		# Calling
		##################
		$self->isCnv($vh->isCnv);
		
		my $hpatients;
		if ($vh->isDude){
			$self->setDudeValues($vh->getChromosome,$patient,$vh);
		}
		elsif ($vh->isCnv){
			$self->setCnvValues($vh->getChromosome,$patient,$vh);
		}
		else {
			foreach my $p (@{$patient->getFamily()->getMembers}) {
				$hpatients->{ $p->id }->{gt} = $vh->getSequencingGenotype($p);
				$hpatients->{ $p->id }->{pc} = $vh->getRatio($p);
				$hpatients->{ $p->id }->{dp} = $vh->getDP($p);
				$hpatients->{ $p->id }->{model} = $vh->getTransmissionModelType($p->getFamily(),$p);
			}
		}
		

		$self->patients_calling($hpatients);
		################################
		#
		# CNV
		#
		####################@
		$self->isCnv($vh->isCnv);
		if ($self->isCnv){
		$self->isDude(1); ### test
		
		$self->patients_calling($self->parseTrioTableCnv($vh->{html}->{trio},$project));
			if ($self->isDude){
				my $chr = $project->getChromosome($self->chromosome());
				my $primers = $project->getPrimersByPosition($chr,$self->start,$self->end);
				my $hpatients;
				  map {$hpatients->{$_->id} = $_} @{$project->getPatients()};
				foreach my $primer (@$primers){
					foreach my $ph (keys %{$self->patients_calling}){
						my ($patient) = $hpatients->{$ph};
							$self->patients_calling->{ $patient->id }->{norm_depth} = $patient->cnv_region_ratio_norm($chr->name,$primer->start,$primer->end);
							$self->patients_calling->{ $patient->id }->{dude_score} = $primer->cnv_score($patient);
					}
				}
			}
			
			
			
		
		}
		
		#################
		# gnomad
		##############@
		
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
}
sub update_clinvar {
	my ($self,$project,$gene) = @_;
		my $vh = $project->_newVariant($self->id);
		 	
		 #	die();
	
		$self->clinvar_id($vh->clinvar_id);
		return unless $vh->clinvar_id;
		$self->clinvar($vh->clinvar_class) ;
		if (exists $vh->genes_pathogenic_DM->{$gene->{id}} && $vh->genes_pathogenic_DM->{$gene->{id}}->{pathogenic} ){
				$self->clinvar_pathogenic($vh->isClinvarPathogenic);
		}
	
}

sub setOldVariant {
	my ($self,$vh,$project,$patient,$gene,$update) = @_; 
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
	$self->parseDejaVuTable($vh->{html}->{deja_vu},$vh);
	$self->dejavu_this_run_patients($vh->{value}->{this_run_patients});# = '-';
	
	#########
	# caller
	#########
	$self->text_caller($vh->{value}->{ngs});
	
	#warn Dumper $vh->{transmission_model};
	
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
	#	warn Dumper $vh->{genes}->{$gid};
	#	die();
	#	warn $gid;
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

