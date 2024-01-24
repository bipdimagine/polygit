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
);

my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => $project_name);
warn $project_name;
my $fams = $project->getFamilies();

foreach my $fam ($fams){
exit(0) if scalar(@$fams) > 1;
my $fam = $fams->[0];
my $mother = $fam->getMother();
my $father = $fam->getFather();
next unless $father;
next unless $mother;

foreach my $child (@{$fam->getChildren}){
my $all_res;
foreach my $chr (@{$project->getChromosomes}){
	#next unless $chr->name eq "X"; 
	
	my $father_he = $father->getVectorHe($chr)->Clone();
	$father_he -= $mother->getVariantsVector($chr);
	my $child_ho = $child->getVectorHo($chr)->Clone();
	$child_ho &= $father_he ;
	my $nb1 = $chr->countThisVariants($father_he);
	my $nb2 = $chr->countThisVariants($child_ho);
	warn _max($child_ho->to_Enum());
	#my $var_par_1_he = $lParents[0]->getVectorHe($chr)->Clone();
	$nb1 +=1 ;
	my $sc1 = int(($nb2/$nb1)*100);
	my $mother_he = $mother->getVectorHe($chr)->Clone();
	$mother_he -= $father->getVariantsVector($chr);
	my $child_ho2 = $child->getVectorHo($chr)->Clone();
	$child_ho2 &= $mother_he ;
	warn _max($child_ho2->to_Enum());
	my $nb3 = $chr->countThisVariants($mother_he);
	my $nb4 = $chr->countThisVariants($child_ho2);
	$nb3++;
	my $sc2 = int(($nb4/$nb3)*100);
	next if $child->sex == 1 && $chr->name eq "X"; 
	
	my $father_all = $father->getVectorHo($chr)->Clone();
	$father_all -= $mother->getVariantsVector($chr);
	my $child_all1 = $child->getVectorHo($chr)->Clone();
	$child_all1 &= $father_all;
	my $nb5 = $chr->countThisVariants($father_all);
	my $nb6 = $chr->countThisVariants($child_all1);
	$nb5 ++;
	my $sc3 = int(($nb6/$nb5)*100)  ;
	
	my $mother_all = $mother->getVectorHo($chr)->Clone();
	$mother_all -= $father->getVariantsVector($chr);
	my $child_all2 = $child->getVectorHo($chr)->Clone();
	$child_all2 &= $mother_all;
	my $nb7 = $chr->countThisVariants($mother_all);
	my $nb8 = $chr->countThisVariants($child_all2);
	$nb7 ++;
	my $sc4 = int(($nb8/$nb7)*100);
	
	my $father_only = $father->getVariantsVector($chr);
	 $father_only -= $mother->getVariantsVector($chr);
	 my $child_only_F  = $child->getVariantsVector($chr);
	$child_only_F &= $father_only;
	
	my $nb_only_F = $chr->countThisVariants($father_only); 
	my $nb_child_only_F =  $chr->countThisVariants($child_only_F);
	$nb_only_F +=1;
	
	
	
	
	my $mother_only = $mother->getVariantsVector($chr);
	$mother_only -= $father->getVariantsVector($chr);
	my $child_only_M  = $child->getVariantsVector($chr);
	$child_only_M &= $mother_only;
	
	my $nb_only_M = $chr->countThisVariants($mother_only); 
	my $nb_child_only_M =  $chr->countThisVariants($child_only_M);
	$nb_only_M ++;
	
	
	#FREQ ALLELIQUE
	my $vector_40 =  $child->getVectorRatio($chr,"ratio_40");
	my $father_fa = $father->getVectorHe($chr);
	$father_fa -= $mother->getVariantsVector($chr);
	$father_fa  &= $child->getVariantsVector($chr);
	my $nb_father_fa =  $chr->countThisVariants($father_fa);
	$father_fa &= $vector_40;
	my $nb_father_fa_40 =  $chr->countThisVariants($father_fa);
	
	
	my $mother_fa = $mother->getVectorHe($chr);
	$mother_fa -= $father->getVariantsVector($chr);
	$mother_fa  &= $child->getVariantsVector($chr);
	my $nb_mother_fa =  $chr->countThisVariants($mother_fa);
	$mother_fa &= $vector_40;
	my $nb_mother_fa_40 =  $chr->countThisVariants($mother_fa);
	$nb_father_fa ++;
	$nb_mother_fa ++;
	
	$all_res->{$child->name}->{$chr->name}->{RATIO_40_F} = int(($nb_father_fa_40/$nb_father_fa)*100);
	$all_res->{$child->name}->{$chr->name}->{RATIO_40_M} = int(($nb_mother_fa_40/$nb_mother_fa)*100);
	
	$all_res->{$child->name}->{$chr->name}->{ONLY_F} = int(($nb_child_only_F/$nb_only_F)*100);
	$all_res->{$child->name}->{$chr->name}->{ONLY_M} = int(($nb_child_only_M/$nb_only_M)*100);
	$all_res->{$child->name}->{$chr->name}->{UPD_F} = $sc1;
	$all_res->{$child->name}->{$chr->name}->{UPD_M} = $sc2;
	$all_res->{$child->name}->{$chr->name}->{nb_UPD_F} = $sc1;
	$all_res->{$child->name}->{$chr->name}->{nb_UPD_F} = "$nb2 / $nb1 ";
	$all_res->{$child->name}->{$chr->name}->{nb_UPD_M} = "$nb4 / $nb3 ";
	$all_res->{$child->name}->{$chr->name}->{HE_UPD_F} = $sc3;
	$all_res->{$child->name}->{$chr->name}->{nb_HE_UPD_F} = "$nb5 / $nb6 ";
	$all_res->{$child->name}->{$chr->name}->{HE_UPD_M} = $sc4;
	$all_res->{$child->name}->{$chr->name}->{nb_HE_UPD_M} = "$nb7 / $nb8 ";
	print $project->name."\t".$chr->name.":F->".$sc1. " = $sc3 ".$all_res->{$child->name}->{$chr->name}->{nb_UPD_F}."\tM->".$sc2." = $sc4 ".$all_res->{$child->name}->{$chr->name}->{nb_UPD_M}."\n" if ($sc1+$sc2) > 10 or ($sc4+$sc3) > 30;
	print "HE => ".$project->name."\t".$chr->name.":F->".$sc3. " F : ".$all_res->{$child->name}->{$chr->name}->{nb_HE_UPD_F}."\tM->".$sc4. "M : ".$all_res->{$child->name}->{$chr->name}->{nb_HE_UPD_M}."\n" if ($sc4+$sc3) > 30;
}
my $json = encode_json $all_res;
my $dir = $project->getSVDir()."/UPD";
system("mkdir -p $dir ; chmod a+rwx $dir") unless -e $dir;

foreach my $c (keys %$all_res){
	open(JSON,">$dir/".$c.".json");
	print JSON  encode_json $all_res->{$c};
	close JSON;
}
}
}

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

exit(0);
