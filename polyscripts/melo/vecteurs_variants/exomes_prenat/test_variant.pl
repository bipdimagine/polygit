#!/usr/bin/perl

use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
use GBuffer;
use Set::IntSpan::Fast::XS;
use Getopt::Long;

my $fork = 1;
my $cmd;
my ($project_name, $patient_name, $vid);
my $file;

GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'cmd=s'  => \$cmd,
	'variant=s'  => \$vid,
	'patient=s'  => \$patient_name,
);

my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => "NGS2023_7054");

my $child1= $project->getPatient("23A2_PATIENT_Het");
my $child2= $project->getPatient("23A3_BROTH_Het");
my $father= $project->getPatient("23NA1_FATH_WT");
my $mother= $project->getPatient("23A1_MOTH_Het");

my $fork = 1;
my $pm   = new Parallel::ForkManager($fork);
	my $hrun;
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $h ) = @_;

			unless ( defined($h) or $exit_code > 0 ) {
				print
				  qq|No message received from child process $exit_code $pid!\n|;
				die();
				return;
			}
			
		print join("\n",@{$h->{lines}}) if $h->{lines};
		}
		
	);
	my $id = time;
	$project->buffer->dbh_deconnect();

my $limit_cad = 20;

foreach my $chr (@{$project->getChromosomes}) {
	my $pid = $pm->start and next;
	warn "start ".$chr->name();
	my $v1 = $child1->getVectorHe($chr);
	$v1 &= $child2->getVectorHe($chr);
	 
	
	
	my $vf = $father->getVectorHe($chr) - $mother->getVectorOrigin($chr);
	my $vh1 =  $v1 & $vf;
	
	my $vm = $mother->getVectorHe($chr) - $father->getVectorOrigin($chr);
	my $vh2 =  $v1 & $vm;
	
	
	
	
	
	my $array1 = to_array($vh1,$chr->name);
	my $array2 = to_array($vh2,$chr->name);
	
	my $nb = 0;
	my $hgenes;
	foreach my $id (@$array1){
		my $v = $project->returnVariants($id);
		#warn $v->name();
		my $keep = 0;
		foreach my $g (@{ $v->getGenes()}){
			 $keep = 1 if $v->max_spliceAI_score($g) > 0.8;
		}
		next if  ($v->cadd_score < $limit_cad && $keep == 0);
		foreach my $g (@{ $v->getGenes()}){
			$hgenes->{$g->id}->{nbf} ++;
			$hgenes->{$g->id}->{external_name} = $g->external_name;
			push(@{$hgenes->{$g->id}->{variants_father}},$id);
		}
		
		 
		$nb ++;
		warn $nb."F" if $nb % 100 ==0; 
	} 
	$nb = 0;
	warn "mother";
	foreach my $id (@$array2){
		my $v = $project->returnVariants($id);
		#warn $v->name();
		my $keep = 0;
		foreach my $g (@{ $v->getGenes()}){
			next unless exists $hgenes->{$g->id};
			$keep = 1 if $v->max_spliceAI_score($g) > 0.8;
		}
		next if  ($v->cadd_score < $limit_cad && $keep ==0);
		foreach my $g (@{ $v->getGenes()}){
			next unless exists $hgenes->{$g->id};
			$hgenes->{$g->id}->{nbm} ++;
			push(@{$hgenes->{$g->id}->{variants_mother}},$id);
		}
		
		$nb ++;
		warn $nb."M" if $nb % 200 ==0; 
	} 
	warn "genes";
	
	foreach my $g (keys %$hgenes){
		delete $hgenes->{$g} unless $hgenes->{$g}->{nbm};
		
	}
	my $line;
	foreach my $g (keys %$hgenes){
		push(@$line,$g."\t".$hgenes->{$g}->{external_name}."\t".join(";",@{$hgenes->{$g}->{variants_mother}})."\t".join(";",@{$hgenes->{$g}->{variants_father} }));
	}
	warn $line;
	$pm->finish( 0, {lines=>$line} );
}
$pm->wait_all_children();

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
