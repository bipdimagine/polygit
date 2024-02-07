#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
#use lib "/software/polyweb/poly-disk/poly-src/GenBo/lib/obj-nodb/packages";
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/packages";
use GBuffer;
use Set::IntSpan::Fast::XS;
use Getopt::Long;
use JSON::XS;

my $fork = 1;
my ($project_name, $patient_name,$vid);
my $file;

GetOptions(
'fork=s'       => \$fork,
'project=s'    => \$project_name,
'patient=s'    => \$patient_name,
);

my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => $project_name);
my $patient = $project->getPatient($patient_name);

unless ($patient->isChild){
	my $all_res ={};
	save_json($patient,$all_res);
	exit(0);
} 
#warn $project_name;
my $fam = $patient->getFamily();
	my $mother = $fam->getMother();
	my $father = $fam->getFather();
	save_json($patient,{}) unless $father;
	save_json($patient,{}) unless $mother;			# il faut les deux parents

  	foreach my $child (@{$fam->getChildren})
  	{
  		next if $child->name ne $patient->name();
		my $all_res;
		foreach my $chr (@{$project->getChromosomes}) {
			next if  $chr->name eq "MT";
			# snp he chez un des parents absent chez l'autre et qui se retrouve ho chez l'enfant (isodisomie)
			my $father_he = $father->getVectorHe($chr)->Clone();
			$father_he -= $mother->getVariantsVector($chr);
			my $child_ho = $child->getVectorHo($chr)->Clone();
			$child_ho &= $father_he ;
			my $nb1 = $chr->countThisVariants($father_he);
			my $nb2 = $chr->countThisVariants($child_ho);
			
			$nb1 +=1 ;
			my $sc1 = int(($nb2/$nb1)*100);  # % de snp he du père retrouvés ho chez l'enfant 

			my $mother_he = $mother->getVectorHe($chr)->Clone();
			$mother_he -= $father->getVariantsVector($chr);
			my $child_ho2 = $child->getVectorHo($chr)->Clone();
			$child_ho2 &= $mother_he ;
			my $nb3 = $chr->countThisVariants($mother_he);
			my $nb4 = $chr->countThisVariants($child_ho2);
			$nb3++;
			my $sc2 = int(($nb4/$nb3)*100); # % de snp he de la mère retrouvés ho chez l'enfant

			next if $child->sex == 1 && $chr->name eq "X"; 

			# snp ho chez un des parents absent chez l'autre et qui se retrouve ho chez l'enfant (isodisomie + heterodisomie)
			my $father_all = $father->getVectorHo($chr)->Clone();
			$father_all -= $mother->getVariantsVector($chr);
			my $child_all1 = $child->getVectorHo($chr)->Clone();
			$child_all1 &= $father_all;
			my $nb5 = $chr->countThisVariants($father_all);
			my $nb6 = $chr->countThisVariants($child_all1);
			$nb5 ++;
			my $sc3 = int(($nb6/$nb5)*100)  ; # % de snp ho du père retrouvés ho chez l'enfant

			my $mother_all = $mother->getVectorHo($chr)->Clone();
			$mother_all -= $father->getVariantsVector($chr);
			my $child_all2 = $child->getVectorHo($chr)->Clone();
			$child_all2 &= $mother_all;
			my $nb7 = $chr->countThisVariants($mother_all);
			my $nb8 = $chr->countThisVariants($child_all2);
			$nb7 ++;
			my $sc4 = int(($nb8/$nb7)*100); # % de snp ho de la mère retrouvés ho chez l'enfant

			# snp present que chez le père transmis à l'enfant 
			my $father_only = $father->getVariantsVector($chr);
			$father_only -= $mother->getVariantsVector($chr);
			my $child_only_F  = $child->getVariantsVector($chr);
			$child_only_F &= $father_only;

			my $nb_only_F = $chr->countThisVariants($father_only); 
			my $nb_child_only_F =  $chr->countThisVariants($child_only_F);
			$nb_only_F ++;


			# snp present que chez la mère transmis à l'enfant
			my $mother_only = $mother->getVariantsVector($chr);
			$mother_only -= $father->getVariantsVector($chr);
			my $child_only_M  = $child->getVariantsVector($chr);
			$child_only_M &= $mother_only;

			my $nb_only_M = $chr->countThisVariants($mother_only); 
			my $nb_child_only_M =  $chr->countThisVariants($child_only_M);
			$nb_only_M ++;

			my $transmissionF= int(($nb_child_only_F/$nb_only_F)*100);
			my $transmissionM= int(($nb_child_only_M/$nb_only_M)*100);

			# region Ho : nb snp ho consécutifs chez enfant provenant d'un des parents
			my $lmax_HoF = _max($child_ho->to_Enum());
			my $lmax_HoM = _max($child_ho2->to_Enum());
			$all_res->{$chr->name}->{TRANSMIS_BYF} = int(($nb_child_only_F/$nb_only_F)*100);
			$all_res->{$chr->name}->{TRANSMIS_BYM} = int(($nb_child_only_M/$nb_only_M)*100);
			$all_res->{$chr->name}->{ISO_UPD_F} = $sc1;
			$all_res->{$chr->name}->{ISO_UPD_M} = $sc2;
			$all_res->{$chr->name}->{ISOHE_UPD_F} = $sc3;
			$all_res->{$chr->name}->{ISOHE_UPD_M} = $sc4;
			$all_res->{$chr->name}->{HoFrom_F} = $lmax_HoF;
			$all_res->{$chr->name}->{HoFrom_M} = $lmax_HoM;
			
		} # chr
		save_json($patient,$all_res);
		
  	} # childs

sub _max {
my ($string) = @_;
my $max =0;
my $set = Set::IntSpan::Fast::XS->new();

foreach my $s (split(",",$string)){
	my ($a,$b) = split("-",$s);
	next unless $b;
	$set->add_range($a-5, $b+5);
}
my $iter = $set->iterate_runs();
my @tt;
my $l =0;

while (my ( $from, $to ) = $iter->()) {
    	$max = abs($from - $to) if $max < abs($from - $to);
}
return $max;

}
sub save_json{
my ($patient,$hash) =@_;
my $file = $patient->upd_file();
open(JSON,">$file");
print JSON  encode_json $hash;
close JSON;
exit(0);
}
exit(0);
	