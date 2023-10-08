package GenBoCnv;

use strict;
use Moo;
use Data::Dumper;
use Config::Std;
use GenBoCapture;
use Carp;
use List::Util qw[min max];
extends "GenBoVariant";



sub length {
	my ($self) = @_;
	return $self->allele_length;
}

has allele_length => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return abs($self->start - $self->end)+1;
	},
);



has name_cnv => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $al = $self->allele_length;
		confess();
		#$al = "*" if $al == -1;
		warn $self->start()." :: ".$self->sv_type;
		die();
		my $id = $self->getChromosome->id().'-'.$self->start().'-'.lc($self->sv_type).'-'.$self->allele_length;
		$id  = $self->getChromosome->id().'-'.$self->start().'-'.lc($self->mei_type).'-'.$self->allele_length if $self->isMei;
		return $id;
	},
);



has mei_type => (
	is		=> 'ro',
	default	=> sub {
		my $self = shift;
		return "";
	},
);

########
# isXXX Methods
########
has isSV => (
	is		=> 'ro',
	default	=> 1,
);
has isCnv => (
	is		=> 'ro',
	default	=> 1,
);

has isLarge => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> 1,
);
has isMei => (
	is		=> 'ro',
	default	=> undef,
);
has isLargeDuplication => (
	is		=> 'rw',
	default	=> undef,
);
has isLargeDeletion => (
	is		=> 'rw',
	default	=> undef,
);
has isLargeInsertion => (
	is		=> 'rw',
	default	=> undef,
);

has isDude => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		my $is_dude = 0;
		foreach my $patient (@{$self->getPatients()}) {
			$is_dude = 1 if $self->isDudeCaller($patient);
		}
		return $is_dude;
	},
);


############
# INTSPAN 
#############

#has intspan => (
#	is		=> 'ro',
#	reader	=> 'getGenomicSpan',
#	lazy	=> 1,
#	default => sub {
#		my $self = shift;
#		my $intSpan = Set::IntSpan::Fast::XS->new();
#		$intSpan->add_range($self->start()-100, $self->start() + $self->length() +100);
#    	return $intSpan;
#	},
#);
#
#has strict_intspan => (
#	is		=> 'ro',
#	reader	=> 'getStrictGenomicSpan',
#	lazy	=> 1,
#	default => sub {
#		my $self = shift;
#		my $intSpan = Set::IntSpan::Fast::XS->new();
#		$intSpan->add_range($self->start(),$self->event_end());
#    	return $intSpan;
#	},
#);




sub annotation_coding {
	my ( $self, $tr, $annot ) = @_;
	my $span =  $self->getGenomicSpan()->intersection( $tr->getSpanCoding );
	my $prot  = $tr->getProtein();
	my $gid   = $tr->getGene->id();
	my $trid = $tr->id;
	my $namep = $prot->name();
	my $pos   = $self->start();
	my $seq   = $self->sequence();
	my @array = $span->as_array();
	my $project = $self->getProject();
	my $start_tr = $array[0];
	my $end_tr = $array[-1];
	my $tres = scalar(@array) % 3;
	my $pos_transcript = $tr->translate_position($pos);
	$end_tr = $start_tr unless $end_tr; 
	my $consequence =  $tr->codonsConsequenceForVariations($self,$start_tr,$end_tr);
	$annot->{ $tr->id }->{coding}->{sequences} = $consequence;
	$annot->{ $prot->id }->{coding}->{sequences} = $consequence;
	unless ( $self->getGenomicSpan()->intersection($tr->getSpanCodonStartEnd())->is_empty){
		$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("phase");
	}
	$annot->{$trid}->{mask} = $annot->{$trid}->{mask} | $project->getMaskCoding("frameshift");
	
}

sub polyphenScore {
	return '-';
}

sub siftScore {
	return '-';
}
sub polyphenStatus {
	my ( $self, $obj ) = @_;
	return 0 ;
}

sub getRatio {
        my ($self,$patient,$method) = @_;
        my $pid = $patient->id;
        return 0;
}
sub siftStatus {
	my ( $self, $obj ) = @_;
	return 0 ;
}


sub getPourcentAllele {
		my ($self,$patient,$method) = @_;
			return $self->getRatio($patient,$method);
}




sub getTransmissionModelType {
		my ($self,$fam,$child,$gene) = @_;
		return "-" unless (exists $self->{patients_object}->{$child->id});
		if ($child->isMother ){
			return "-m" unless exists $self->{patients_object}->{$child->id};
			return "+m"; 
		}
		if ($child->isFather ){
			return "-f" unless exists $self->{patients_object}->{$child->id};
			return "+f"; 
		}
		#my $is_mosaic = $self->isMosaicTransmission($fam,$child);
		if ($self->isRecessiveTransmission($fam,$child) == 1 ){
			return "recessive";
		}
		elsif ($self->isUniparentalDisomyTransmission($fam,$child) == 1  && $self->getPourcentAllele($child) > 75 ){
			
			return "Uniparental Disomy";	
		}
		elsif  ($self->isMosaicTransmission($fam,$child)  ){
			return "Mosaic";
		}
		elsif ($self->isBothTransmission($fam,$child) == 1  ){
			return "Both"
		}
		elsif ($self->isFatherTransmission($fam,$child) ==1){
			if ($gene and $self->isCompoundTransmission($fam,$child,$gene,'mother') == 1){
				return "Father_c";
			}
			return "Father";
		}
		elsif ($self->isMotherTransmission($fam,$child) == 1){
			if ($gene and $self->isCompoundTransmission($fam,$child,$gene,'father') == 1){
				return "Mother_c";
			}
			return "Mother";
		}
		elsif  ($self->isStrictDenovoTransmission($fam,$child) == 1  ){
			return "strict_denovo";
		}
		elsif  ($self->isDenovoTransmission($fam,$child) == 1  ){
			
			if ($fam->father && $fam->mother){
				return "denovo";
			}
			else {
				return "denovo/?";
			}
			
		}
		
		else {
			return "?";
		}
		confess("code it by your self")
}


sub getNGSScore{
	return 1;
}
sub spliceAI {
	my ($self) = @_;
	return {};
}


sub getTransmissionModel_cnv {
	my ($self,$fam,$child) = @_;
	return $self->{trans_cnv}->{$child->name} if exists $self->{trans_cnv}->{$child->name};
		
	my $mother = $fam->getMother();
	my $father = $fam->getFather(); 
	unless ($mother && $father){
		$self->{trans_cnv}->{$child->name} = "-";
		return $self->{trans_cnv}->{$child->name};
	}
	my $h ={};
	my $list = $self->list_same_cnv(0.9);
	foreach my $id (keys %$list){
		$h->{mother}->{he} ++ if $mother && $mother->getVectorOriginHe($self->getChromosome)->contains($id);
		$h->{father}->{he} ++  if $father && $father->getVectorOriginHe($self->getChromosome)->contains($id);
		$h->{mother}->{ho} ++ if $mother && $mother->getVectorOriginHo($self->getChromosome)->contains($id);
		$h->{father}->{ho} ++  if $father && $father->getVectorOriginHo($self->getChromosome)->contains($id);
	}
	my ($manta) = grep {$_ =~ /manta/} split("\n",Dumper $self->{annex});
	
	if ($self->project->isGenome && !($self->is_manta)){
		if ($mother && !exists $h->{mother} ){
			my $sd = $mother->sd_value_dude($self->getChromosome->name,$self->start,$self->start+$self->length);
			my $cnv_value_dude = $mother->cnv_value_dude($self->getChromosome->name,$self->start,$self->start+$self->length);
			my $limit = 30;
			$limit= 200 if $mother->getVectorOriginHo($self->getChromosome)->contains($self->vector_id);
			if ($sd < $limit && $self->limit_cnv_value_ho($cnv_value_dude) ){
				$h->{mother}->{ho} ++;
			}
			elsif ($sd < $limit && $self->limit_cnv_value($cnv_value_dude)  ){
				$h->{mother}->{he} ++;
			}
		}
		if ($father && !exists $h->{father} ){
			my $sd = $father->sd_value_dude($self->getChromosome->name,$self->start,$self->start+$self->length);
			my $cnv_value_dude = $father->cnv_value_dude($self->getChromosome->name,$self->start,$self->start+$self->length);
			#warn $sd." ".$cnv_value_dude;
			my $limit = 30;
			$limit= 200 if $father->getVectorOriginHo($self->getChromosome)->contains($self->vector_id);
			if ($sd < $limit && $self->limit_cnv_value_ho($cnv_value_dude) ){
				$h->{father}->{ho} ++;
			}
			elsif ($sd < $limit && $self->limit_cnv_value($cnv_value_dude)  ){
				$h->{father}->{he} ++;
			}
			
		}
	}
	if (!exists $h->{mother} && !exists $h->{father}){
		my $h2 ={};
		my $lis2 = $self->list_same_cnv(0.75);
		foreach my $id (keys %$list){
			$h2->{mother}->{he} ++ if $mother && $mother->getVectorOriginHe($self->getChromosome)->contains($id);
			$h2->{father}->{he} ++  if $father && $father->getVectorOriginHe($self->getChromosome)->contains($id);
			$h2->{mother}->{ho} ++ if $mother && $mother->getVectorOriginHo($self->getChromosome)->contains($id);
			$h2->{father}->{ho} ++  if $father && $father->getVectorOriginHo($self->getChromosome)->contains($id);
		}
		if (!exists $h2->{mother} && !exists $h2->{father}){ $self->{trans_cnv}->{$child->name} = "strict_denovo"; }
		else { $self->{trans_cnv}->{$child->name} = "denovo"; }
	}
	elsif (!exists $h->{mother} && exists $h->{father}){
		$self->{trans_cnv}->{$child->name} = "father";
	}
	elsif (!exists $h->{father} && exists $h->{mother}){
		$self->{trans_cnv}->{$child->name} = "mother";
	}
	
	elsif ((exists $h->{mother}->{he} or exists $h->{father}->{he}) && $self->isHomozygote($child)){
		$self->{trans_cnv}->{$child->name} = "recessive";

	}
	else {
		$self->{trans_cnv}->{$child->name} = "both";
	}
	#warn $self->{trans_cnv}->{$child->name};
	return $self->{trans_cnv}->{$child->name};
}

sub isRecessiveTransmission {
	my ($self,$fam,$child) = @_;
	return 1 if $self->getTransmissionModel_cnv($fam,$child) eq "recessive";
	
}

sub isFatherTransmission {
	my ($self,$fam,$child) = @_;
	return 1 if $self->getTransmissionModel_cnv($fam,$child) eq "father";
}
sub isMotherTransmission {
	my ($self,$fam,$child) = @_;
	return 1 if $self->getTransmissionModel_cnv($fam,$child) eq "mother";
}
sub isBothTransmission {
	my ($self,$fam,$child) = @_;
	return 1 if $self->getTransmissionModel_cnv($fam,$child) eq "both";
}

sub isDenovoTransmission {
	my ($self,$fam,$child) = @_;
	return 1 if $self->getTransmissionModel_cnv($fam,$child) eq "strict_denovo";
	return 1 if $self->getTransmissionModel_cnv($fam,$child) eq "denovo";
}

sub isStrictDenovoTransmission {
	my ($self,$fam,$child) = @_;
	return 1 if $self->getTransmissionModel_cnv($fam,$child) eq "strict_denovo";
}

sub score_variants_trio_dominant{
	my ($self,$fam,$child,$score) = @_;

	$fam = $child->getFamily();

	my $mother = $fam->getMother();
	my $father = $fam->getFather(); 

	if ($mother && $mother->isHealthy && $self->isMotherTransmission($fam,$child)){
		$score -=2 ;
		return $score;
	}
	if ($father && $father->isHealthy && $self->isFatherTransmission($fam,$child)){
		$score -=2 ;
		return $score;
	}
	if ($self->isDenovoTransmission($fam,$child)){
		return $score;
	}
	return $score*2;
}


sub score_variants_trio {
	my ($self,$child,$score,$tr,$debug) = @_;	
	#return $self->{score}->{trio}->{$child->name} if exists $self->{score}->{trio}->{$child->name};
	my $fam = $child->getFamily();
	return $self->score_variants_trio_dominant($fam,$child,$score) if $fam->isDominant();
	#TODO: ATTENTION il peut avoir frere ou soeur -> le faire sur le vecteur FAMILLE et pas du PATIENT.
	if ($self->isRecessiveTransmission($fam,$child)){
		$score += 1;
		return $score*2;
	}
	
	
	my $isMotherTransmission = $self->isMotherTransmission($fam,$child);
	my $isFatherTransmission = $self->isFatherTransmission($fam,$child);


	
	if  ($isMotherTransmission or $isFatherTransmission) {
		#IDeal $self->getFamily->isFatherTransmission($fam,$child)){
	 	$score += 1;
		return $score;

	}
	
	if ($self->isStrictDenovoTransmission($fam,$child)){
			$score += 1;
			return  $score*2 ;
	}
	
	if ($self->isDenovoTransmission($fam,$child)){
			
			return  $score*2 ;#if $fam->getFather && $fam->getMother;
	}
	
	
	return $score - 1;
			
}

sub score_refined {
	my ($self,$patient,$debug) = @_;
		my $refined_score = 0;
	
	return $refined_score unless ($self->getProject->isGenome());
	
	my $pr = $self->other_projects + $self->similar_projects();
	if ($debug ){
		warn "++++ --------> ".$pr;
		#die();
	}
			# return {mean=>0,} unless defined $res->{_};
			
	my @data;
	if($self->project->isGenome){
	my $sd = $patient->sd_value_dude($self->getChromosome->name,$self->start,$self->start+$self->length);
	my $cnv_value_dude = $patient->sd_value_dude($self->getChromosome->name,$self->start,$self->start+$self->length);
	if ($sd > 40 ){
		$refined_score  -= 0.5;
	}
	if ($cnv_value_dude > 0.7 && $cnv_value_dude<1.3 ){
		$refined_score  -= 0.5;
	}
	if ($sd < 20 && ($cnv_value_dude < 0.6 or $cnv_value_dude>1.4 )){
		$refined_score  += 1;
	}
	}
	$pr = 0 unless $pr;
		if ($pr >= 5){
						$refined_score  -= 0.5;# if $v->{max_pc} > 15;
		}
		if ($pr > 10 ){
						$refined_score  -= 1;# if $v->{max_pc} > 15;
		}
		if ($pr > 15  ){
						$refined_score  -= 1;# if $v->{max_pc} > 15;
		}
		if ($pr > 20  ){
				$refined_score  -= 1 ;# if $v->{max_pc} > 15;
		}
			if ($pr > 50  ){
				$refined_score  -= 2 ;# if $v->{max_pc} > 15;
		}
		
		
		return $refined_score;
}

has is_imprecise => (
	is		=> 'rw',
	lazy	=> 1,
	default => undef,
);

sub list_same_cnv {
	my ($self, $seuil) = @_;
	my $tree = $self->return_interval_tree;
	my $start = $self->start();
	my $end = $self->start() + $self->length();
	$end++ if ($start == $end);
	my $results = $tree->fetch($start,$end);
	my $hreturns;
	foreach my $r (@$results){
		#next if $r->[2] == $self->vector_id;
		if ($self->buffer->areSameCNV($start,$end,$r->[0],$r->[1],$seuil)){
			$hreturns->{$r->[2]} ++;
		}
	}
	return $hreturns;
}


sub dejaVuInfosForDiag2 {
	my ($self,$key) = @_;
	unless (exists $self->{array_dejavu}) {
		my $hash = $self->getProject->getDejaVuInfosForDiag($self->id());
		$self->{array_dejavu} = $self->return_dejavu;
	}
	return $self->{array_dejavu} unless $key;
	my $index = $self->buffer->index_dejavu($key);
	return  $self->{array_dejavu}->[$index];
}
sub return_dejavu {
	my $self = shift;
	my $project = $self->project;	
	return $self->SUPER::dejaVuInfosForDiag2() if ($self->isDude());
	my $chr = $self->getChromosome()->name();
	my $in_this_run_patients =  $self->project->in_this_run_patients();
	my $no = $project->dejavuSV();

	if ($self->length == 0) {
		warn $self->end;
		warn $self->length;
		warn Dumper $self->getPatientsInfos();
		warn $self->sequence;
		foreach my $p (@{$self->getPatients}){
			warn $p->name;
		}
	}
	my $h = $no->get_cnv($self->getChromosome->name,$self->start,$self->start+$self->length,$self->sv_type,'all',75);
	my $similar = $project->similarProjects();
	my $exomes = $project->exomeProjects();
	my $pe =  $project->countExomePatients();
	my $ps =  $project->countSimilarPatients();
	my $res;
	$res->{similar_projects} = 0;
	$res->{other_projects} = 0;
	$res->{exome_projects} = 0;
	$res->{other_patients} = 0;
	$res->{exome_patients} = 0;
	$res->{similar_patients} = 0;
	$res->{other_patients_ho} = 0;
	$res->{exome_patients_ho} = 0;
	$res->{similar_patients_ho} = 0;
	$res->{total_in_this_run_patients} = $in_this_run_patients->{total} + 0;
	if ($res->{total_in_this_run_patients} == 0 ){
		$res->{total_in_this_run_patients} = scalar(@{$project->getPatients});
	}
	$res->{in_this_run_patients} = 0;
	$res->{in_this_run_patients} += scalar(@{$self->getPatients});
	return $self->buffer->hash_to_array_dejavu($res) unless ($h);
	foreach my $p (keys %$h) {
		my $nhe = scalar(keys %{$h->{$p}});
		my $nho = 0;
		next if $p eq $project->name();
		#IN THIS RUN 
		
		#IN EXOME 	
		if (exists $exomes->{$p}){
			$res->{exome_projects}  ++;
			$res->{exome_patients}   += $nhe;
			$res->{exome_patients_ho}   += $nho;
		}
		#in similar;
		if (exists $similar->{$p}){
			$res->{similar_projects}  ++;
			$res->{similar_patients} += $nhe;
			$res->{similar_patients_ho} += $nho;
		}
		else {
			$res->{other_projects} ++;
			$res->{other_patients}+= $nhe;
			$res->{other_patients_ho}+= $nho;
		}
	}	
	$res->{total_exome_projects} =  scalar(keys %{$project->exomeProjects()});
	$res->{total_exome_patients} =  $project->countExomePatients();
	$res->{total_similar_projects} =  scalar(keys %{$project->similarProjects()});
	$res->{total_similar_patients} =  $project->countSimilarPatients();
	return  $self->buffer->hash_to_array_dejavu($res);
}


sub get_dude_scores {
	my ($self, $patient) = @_;
	my @lRes;
	my @lDudeFiles = @{$patient->getDudeFiles()};
	foreach my $file (@lDudeFiles) {
		
		my $tabix = Bio::DB::HTS::Tabix->new( filename => $file );#new Tabix( -data => $file );
		my $res = $tabix->query_full($self->getChromosome()->name(), ($self->start()-1), ($self->start + $self->length() + 1));
		while ( my $line = $res->next ) {
			my ($chr_name,$start,$end,$status,$type,$gene,$score1,$score2,$score3,$score4) = split(" ",$line);
			my $h_pr;
			$h_pr->{id} = $chr_name.':'.$start.'-'.$end;
			$h_pr->{status} = $status;
			$h_pr->{type} = $type;
			$h_pr->{score1} = $score1;
			push(@lRes, $h_pr);
		}
	}
	return \@lRes;
}

has get_genes_dude_all_patients => (
	is		=> 'rw',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		my @lGenes = @{$self->getGenes()};
		my $h_fam_genes_dude;
		return $h_fam_genes_dude unless (@lGenes);
		foreach my $pat (@{$self->getProject()->getPatients()}) {
			my $no = $pat->getGenesDude();
			foreach my $g (@lGenes) {
				$h_fam_genes_dude->{$pat->name()}->{$g->id()} = $no->get($g->id());
			}
			$no->close();
		}
		return $h_fam_genes_dude;
	},
);

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

sub cnv_confidence {
	my ($self, $patient) = @_;
	
	my $sr = $self->sr($patient);
	my $pr = $self->pr($patient);
	my @lPr = split(',', $pr);
	my @lSr = split(',', $sr);
	return "low" unless @lPr;
	return "low" unless @lSr;
	return "low" unless $lPr[1];
	return "low" unless $lSr[1];
	my $sum_alt_sr_pr = $lPr[1] + $lSr[1];
	my $cnv_len = $self->length();
	my $cnv_confidence = 1;
	if ($sr eq '-' or $pr eq '-') { $cnv_confidence = 1; }
	elsif ($cnv_len >= 5000 and $sum_alt_sr_pr > 8 and $lPr[1] > 0 and $lSr[1] > 0) { $cnv_confidence = 3; }
	elsif ($cnv_len >= 1000 and $sum_alt_sr_pr > 13 and $lPr[1] > 0 and $lSr[1] > 0) { $cnv_confidence = 3; }
	elsif ($cnv_len >= 500 and $sum_alt_sr_pr > 8) { $cnv_confidence = 2; }
	my $cnv_dude_score = $patient->cnv_value_dude($self->getChromosome->name, $self->start, $self->start+$self->length);
	if ($self->type() eq 'large_deletion') {
		if ($cnv_dude_score < 0.7) { $cnv_confidence += 1; }
		elsif ($cnv_dude_score > 1.3) { $cnv_confidence -= 2; }
		elsif ($cnv_dude_score > 1.1) { $cnv_confidence -= 1; }
	}
	elsif ($self->type() eq 'large_duplication') {
		if ($cnv_dude_score > 1.3) { $cnv_confidence += 1; }
		elsif ($cnv_dude_score < 0.7) { $cnv_confidence -= 2; }
		elsif ($cnv_dude_score < 0.9) { $cnv_confidence -= 1; }
	}
	return 'high' if ($cnv_confidence >= 3);
	return 'medium' if ($cnv_confidence == 2);
	return 'low' if ($cnv_confidence == 1);
	return '-';
}


has structural_type => (
	is		=> 'ro',
	lazy 	=> 1,
	default	=> sub {
		my $self = shift;
		return  lc($self->sv_type);
	},
);
1;