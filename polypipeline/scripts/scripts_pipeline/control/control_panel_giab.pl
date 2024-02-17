#!/usr/bin/perl
use FindBin qw($Bin $RealBin);
use strict;

use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages/"; 
use GBuffer ;
use Data::Dumper;
use Getopt::Long;
use Carp;
use Bio::DB::Sam;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use List::Util qw(sum);
use File::Temp;
 use Time::Elapsed qw( elapsed );
 use Time::ETA;
use Storable qw(store retrieve freeze);
use JSON::XS;

my $buffer = GBuffer->new();

my $project_name;
my $p1name;
my $p2name;
GetOptions(
	'project=s' => \$project_name,
	'control=s' => \$p1name,
	'giab=s' => \$p2name,
	#'low_calling=s' => \$low_calling,
);
#NGS2021_3851


my $project = $buffer->newProjectCache( -name => "$project_name" );

my $dir_out = $project->getVariationsDir("control_giab");
my $p;
unless ($p1name){
	die() if scalar(@{$project->getPatients})ne 2;
	($p) = grep{$_->name ne $p2name} @{$project->getPatients};
}
#my $p = $project->getPatient("$p1name");
my $giab =  $project->getPatient("$p2name");
my $GIAB_DIR = "/data-isilon/public-data/repository/HG19/GIAB";
if ($p2name =~ /001/){
	$GIAB_DIR = "/data-isilon/public-data/repository/HG19/GIAB/HG001";
}
elsif ($p2name =~ /002/){
	$GIAB_DIR = "/data-isilon/public-data/repository/HG19/GIAB/HG002";
}
else {
	die($p2name." not foud HG001 or HG002");
}
die() unless -e $GIAB_DIR;
my $hintspan_giab = bed_to_intspan("$GIAB_DIR/bed.gz");
my $hintspan_global;
my $hintspan_project;
my $nb_bases;
foreach my $chr (@{$project->getChromosomes}){
	next if $chr->name eq "Y";
	next unless exists $hintspan_giab->{$chr->name};
	$hintspan_project->{$chr->name} = $chr->getIntSpanCapture(20);
	$hintspan_global->{$chr->name} = $hintspan_project->{$chr->name}->intersection($hintspan_giab->{$chr->name});
	my @arary = $hintspan_global->{$chr->name}->as_array();
	$nb_bases += scalar(@arary);
}



my @errors;

my $vs = $p->getIndels();
my $nb =0;
my $nbf = 0;
my %ids;

my $vector_fale_positive;
my $nb_total_indels;
foreach my $v (@$vs){
	next unless exists $hintspan_giab->{$v->getChromosome->name};
	 next unless $hintspan_global->{$v->getChromosome->name}->contains($v->start);
	# next if $p->depth($v->getChromosome->name,$v->start,$v->end)->[0] < 10;
	$nb_total_indels ++;
	next if $p->meanDepth($v->getChromosome->name,$v->start+5,$v->end+5) < 15;
	next if $v->sequencing_details($p)->{hap}->{pourcent} < 20;
	
		 $nb ++;
	warn $v->name;
	#$ids{$v->id}= $v->vector_id;
	if (scalar(@{$v->getPatients}) == 1){
		#next if $v->sequencing_details($p)->{hap}->{nb_ref} + $v->sequencing_details($p)->{hap}->{nb_alt} < 15;
		
		$nbf++;
		print_variant($v,\@errors);	
	
	#	delete $ids{$v->id};
	}
	
	
}

my $fpi = $nbf; 
print "false postive Indel $nbf / $nb \n";
warn "nb_total control indels : ".$nb_total_indels;

 $vs = $p->getVariations();
 $nb =0;
 $nbf = 0;
my $nb_total_vars;
warn scalar(@$vs);
foreach my $v (@$vs){
		
	next unless exists $hintspan_giab->{$v->getChromosome->name};
	next unless $hintspan_global->{$v->getChromosome->name}->contains($v->start);
	 $nb_total_vars ++;
	 #	 next if $p->depth($v->getChromosome->name,$v->start,$v->end)->[0] < 50;
	
	 	next if $p->meanDepth($v->getChromosome->name,$v->start-2,$v->end+2) < 15;
	#	warn Dumper ($v->sequencing_details($p))." ".$v->id if $v->sequencing_details($p)->{hap}->{pourcent} == undef;
		next if $v->sequencing_details($p)->{hap}->{pourcent} < 1;
		#if ($v->id eq  "9_215057_T_C"){
		#my $patients =  $v->getPatients;
		#warn scalar(@$patients);
		#warn $patients->[0]->name();
		#die();
	#}
	 $nb ++;
	
	#warn $v->name;
	

	if (scalar(@{$v->getPatients}) == 1){
		#next if $v->sequencing_details($p)->{hap}->{nb_ref} + $v->sequencing_details($p)->{hap}->{nb_alt} < 15;

		$nbf++;
		print_variant($v,\@errors);	
	#	delete $ids{$v->id};
	}
	else {
		$ids{$v->id}= $v->vector_id;
	}
	
	
}

my $fps = $nbf; 
print "false positive Substitutions $nbf / $nb \n";
warn "nb_total control variations : ".$nb_total_vars;

#die();

warn "end";
my $hi;
my $nb2v =0;
my $nb2vindel =0;
my $nb2fv =0;
my $nb2findel =0;
foreach my $chr (@{$project->getChromosomes}) {
	next if $chr->name eq "Y";
	next unless exists $hintspan_giab->{$chr->name};
	warn $chr->name();
	my $v1 = $chr->getVectorByIntspan( $hintspan_global->{$chr->name});
	
	my $vgiab = $giab->getVectorOrigin($chr);
	my $vp = $p->getVectorOrigin($chr);
	#$v1 &= $vgiab;
	$v1 &= $vgiab;
	
	my $z = to_array($v1,$chr->name);
	
	foreach my $vid (@$z){
		my $v = $project->returnVariants($vid);
		next unless $hintspan_global->{$v->getChromosome->name}->contains($v->start);
		
		#my $f = grep {$_->name eq $giab->name} @{$v->getPatients}; 
		#next unless $f;
		if (scalar(@{$v->getPatients}) == 2){
			next if $v->sequencing_details($p)->{hap}->{pourcent} < 1;
			#next if $p->depth($v->getChromosome->name,$v->start,$v->end)->[0] < 15;
		}
		if($v->isVariation) {
					$nb2v++;
		}
		else {
			$nb2vindel ++;
		}
		
		next if $p->depth($chr->name,$v->start,$v->end)->[0] < 10;
		
		
		#next if $v->sequencing_details($p)->{hap}->{pourcent} < 20;
		
		if (scalar(@{$v->getPatients}) == 1){
				print_variant_false($v,\@errors);
					if($v->isVariation){
						$nb2fv++;
					}
					else {
						$nb2findel ++;
					}
				
				#delete $ids{$v->id};
		}
	} 
	
}
my $results;
$results->{global}->{indels}->{fp} = $fpi;
$results->{global}->{indels}->{fn} = $nb2findel;
$results->{global}->{indels}->{nb} = $nb_total_indels;  
$results->{global}->{substitution}->{fp} = $fps;
$results->{global}->{substitution}->{fn} = $nb2fv;
$results->{global}->{substitution}->{nb} = $nb_total_vars;  

warn $nb_bases;
warn  $nb_total_vars;
$results->{patient}= $p->name();
$results->{date}= time;
$results->{stats}= "variations ".($nb_total_vars+$nb_total_indels)." bases $nb_bases <br> ";
$results->{errors}= \@errors;
my $file = $dir_out."/".$project->name.".giab.json";
warn "write $file"; 
open (JSON,">$file");
print JSON encode_json $results;
close(JSON);
print "false negative Substitutions $nb2fv /  $nb2v \n";
print "false negative Indels $nb2findel /  $nb2vindel \n";
print "nb de bases : $nb_bases\n";




 sub bed_to_intspan {
 	my ($bed) = @_;
 	my $tabix = $buffer->software("tabix");
 	my $hintspan;
 	my @list = `tabix -l $bed`;
 	chomp(@list);
 	foreach my $chr (@list){
 		my $c = $project->getChromosome($chr);
 		my $cname = $c->name();
 		open (BED,"$tabix $bed $chr | ");
 		$hintspan->{$cname} = Set::IntSpan::Fast::XS->new();
 		while(<BED>){
 			chomp();
 			my ($c,$start,$end) = split(" ",$_);
 			$hintspan->{$cname}->add_range($start,$end);
 		}
 		
 	}
 	
 	return $hintspan;
 	
 }
sub to_array {
	my ( $v, $name ) = @_;
	my $set  = Set::IntSpan::Fast::XS->new( $v->to_Enum );
	my $iter = $set->iterate_runs();
	my @t;
	while ( my ( $from, $to ) = $iter->() ) {
		for my $member ( $from .. $to ) {
			if ($name) {
				push( @t, $name . "!" . $member );
			}
			else {
				push( @t, $member );
			}
		}
	}
	return \@t;
}

sub print_variant_false {
	my ($v,$tab) = @_;
	my @vs;
	my $vh;
	$vh->{chromosome} = $v->getChromosome->fasta_name;
	$vh->{event} = "false negative";
	$vh->{id} = $v->id;
	$vh->{type} = "Indel";
	$vh->{type} = "SNP" if $v->isVariation;
	$vh->{vtype} = $v->type;
	$vh->{locus} = $v->getChromosome->fasta_name.":".$v->start."-".$v->end;
	push(@vs,$v->getChromosome->name);
	#push(@vs,$v->start);
	push(@vs,$v->id);
	my @t= split("_",$v->id);
	#unless ($v->type){
		my $chr = $v->getChromosome();
		my $vid = $v->vector_id;
#		$vid++;
#		 my $gid = $chr->cache_lmdb_variations->get_varid($vid);
#		 #warn $gid." ".$vid." ".$id;
#		 my $obj = $chr->cache_lmdb_variations->get($gid);
#		 $vid-=2;
#		 warn $vid;
#		 my $gid = $chr->cache_lmdb_variations->get_varid($vid);
#		 #warn $gid." ".$vid." ".$id;
#		 my $obj = $chr->cache_lmdb_variations->get($gid);
#		 warn $obj->sequence;
#		 warn $v->sequence;
	#	 die();
		
		
	#} 
	my $seqref = $v->getChromosome()->getSequence($v->start,$v->start);
	my $seqalt = $v->sequence();
	$seqalt = $v->alternate_allele if $seqalt eq "-";
	$vh->{text_sequence} = $v->getChromosome()->getSequence($v->start-5,$v->start-1)."_".$seqref."[".$seqalt."]_".$v->getChromosome()->getSequence($v->start+1,$v->start+21);
	
	push(@vs,$v->getChromosome()->getSequence($v->start-5,$v->start-1)."_".$seqref."[".$seqalt."]_".$v->getChromosome()->getSequence($v->start+1,$v->start+21));
	$vh->{giab} = $v->sequencing_details($giab)->{hap}->{nb_ref}."/".$v->sequencing_details($giab)->{hap}->{nb_alt};
	push(@vs,$v->sequencing_details($giab)->{hap}->{nb_ref}."/".$v->sequencing_details($giab)->{hap}->{nb_alt});
	$vh->{control} = "./.";
	push(@vs,"./.");
	$vh->{control_depth} = $p->depth($v->getChromosome->name,$v->start,$v->end)->[0];
	
	push(@vs,$p->depth($v->getChromosome->name,$v->start,$v->end)->[0]);
	push(@$tab,$vh);
	print join("\t",@vs)."\n";
	
}

sub print_variant {
	my ($v,$tab) = @_;
	my @vs;
	my $vh;
	
	$vh->{chromosome} = $v->getChromosome->name;
	$vh->{event} = "false positive";
	$vh->{id} = $v->id;
	$vh->{type} = "Indel";
	$vh->{type} = "SNP" if $v->isVariant;
	$vh->{locus} = $v->getChromosome->fasta_name.":".$v->start."-".$v->end;
	push(@vs,$v->getChromosome->name);
	#push(@vs,$v->start);
	push(@vs,$v->id);
	my @t= split("_",$v->id);
	
	my $seqref = $v->getChromosome()->getSequence($v->start,$v->start);
	my $seqalt = $v->sequence();
	$vh->{vtype} = $v->type;
	$seqalt = $v->alternate_allele if $seqalt eq "-";
	$vh->{text_sequence} = $v->getChromosome()->getSequence($v->start-5,$v->start-1)."_".$seqref."[".$seqalt."]_".$v->getChromosome()->getSequence($v->start+1,$v->start+21);
	push(@vs,$v->getChromosome()->getSequence($v->start-5,$v->start-1)."_".$seqref."[".$seqalt."]_".$v->getChromosome()->getSequence($v->start+1,$v->start+21));
	$vh->{giab} = "./.";
	$vh->{control} = $v->sequencing_details($p)->{hap}->{nb_ref}."/".$v->sequencing_details($p)->{hap}->{nb_alt};
	$vh->{control_depth} = $p->depth($v->getChromosome->name,$v->start,$v->end)->[0];
	warn "coucou ".$v->type;
	#unless ($v->type){
	warn "cuicuiu ->".$v->vector_id;
		my $chr = $v->getChromosome();
	#	my $vid = $v->vector_id;
	#	$vid++;
	#	 my $gid = $chr->cache_lmdb_variations->get_varid($vid);
	#	 warn $gid." ==> ".$vid." =+> ".$v->id;
	#	 unless ($gid)
		 
	#	 my $obj = $chr->cache_lmdb_variations->get($gid);
	#	  warn "--++++--";
	#	 warn $obj->sequence;
	#	 warn $v->sequence;
	#	 $vid-=2;
	#	 my $gid = $chr->cache_lmdb_variations->get_varid($vid);
	#	 #warn $gid." ".$vid." ".$id;
	#	 my $obj = $chr->cache_lmdb_variations->get($gid);
	#	 warn $obj->sequence;
	#	 warn $v->sequence;
	
	push(@$tab,$vh);
 	push(@vs,"./.");
 	
	push(@vs,$v->sequencing_details($p)->{hap}->{nb_ref}."/".$v->sequencing_details($p)->{hap}->{nb_alt});
	push(@vs,$p->depth($v->getChromosome->name,$v->start,$v->end)->[0]);
	print join("\t",@vs)."\n";
	
}
