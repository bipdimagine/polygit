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

my $buffer = GBuffer->new();
my $project = $buffer->newProjectCache( -name => "NGS2019_2761" );


my $patient = $project->getPatient("AJ_SON");
my $patient_control =$project->getPatient("AJ_SON_HC");

my $hintspan ;

 my $chr = $project->getChromosome(10);
 
my $vector_control = $patient_control->getVectorOrigin($chr);
 $vector_control -= $patient->getVectorOrigin($chr);
 my $intspan = 	$chr->getCapturesGenomicSpan();
  my $vector_control_intersect = $patient_control->getVectorOrigin($chr);
  $vector_control_intersect &= $patient->getVectorOrigin($chr);
  my $nb =0;
   foreach my $v_id (@{$chr->getListVarVectorIds($vector_control_intersect)}) {
   	$nb ++;
   }
   warn $nb;
 
 my $nb_30 = 0; 
 foreach my $v_id (@{$chr->getListVarVectorIds($vector_control)}) {
 	my $v = $chr->getVarObject($v_id);
 	
	next unless $intspan->contains($v->start); 
	$nb_30 ++;
 	next if  $patient->meanDepth($chr->name,$v->start,$v->end) < 10;
 	
 	warn $v_id." ".$chr->name." ".$v->start." ".$v->end." ".$patient->meanDepth($chr->name,$v->start,$v->end).' '.$v->sequence;
 }
 warn $nb_30;

 my $vector_patient = $patient->getVectorOrigin($chr);
 
$vector_patient -= $patient_control->getVectorOrigin($chr);



foreach my $v_id (@{$chr->getListVarVectorIds($vector_patient)}) {
	my $v = $chr->getVarObject($v_id);
	#warn $patient->meanDepth($chr->name,$v->start,$v->end);
	#warn $v->start if $patient->meanDepth($chr->name,$v->start,$v->end) > 30;
}