#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;

use strict;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/validation_variation"; 
use lib "$Bin/../cache_nodb/scripts/";
use GBuffer;
use JSON;



my $cgi = new CGI();
my $project_name	= $cgi->param('project');
my $patient_1		= $cgi->param('patient_1');
my $patient_2		= $cgi->param('patient_2');
my $chr_id			= $cgi->param('chr');
my $interval = $cgi->param('interval_length');
my $only_ho = $cgi->param('only_ho');
my $only_rare = $cgi->param('only_rare');


$interval = 100 unless ($interval);

my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => $project_name );

my $chr = $project->getChromosome($chr_id);
my $v_sum = $chr->getNewVector();
my $v_intersect = $chr->getNewVector();
my $v_intersect_he = $chr->getNewVector();
my $v_intersect_ho = $chr->getNewVector();


my $p1 = $project->getPatient($patient_1);
my $p2 = $project->getPatient($patient_2);

$v_sum += $p1->getVectorOrigin($chr);
$v_sum += $p2->getVectorOrigin($chr);

if ($only_rare) {
	my $v_freq = $chr->getNewVector();
	$v_freq += $chr->global_categories->{freq_01};
	$v_freq += $chr->global_categories->{freq_001} if (exists $chr->global_categories->{freq_001});
	$v_freq += $chr->global_categories->{freq_0001} if (exists $chr->global_categories->{freq_0001});
	$v_sum->Intersection($v_sum, $v_freq);
	$v_intersect->Intersection($v_intersect, $v_freq);
	$v_intersect_he->Intersection($v_intersect_he, $v_freq);
	$v_intersect_ho->Intersection($v_intersect_ho, $v_freq);
}

#$nb_interval = 1 unless ($project->isGenome());

$v_intersect->Intersection($p1->getVectorOrigin($chr), $p2->getVectorOrigin($chr));
$v_intersect_he->Intersection($p1->getVectorOriginHe($chr), $p2->getVectorOriginHe($chr));
$v_intersect_ho->Intersection($p1->getVectorOriginHo($chr), $p2->getVectorOriginHo($chr));

print $cgi->header('text/json-comment-filtered');
print "{\"progress\":\".";
my (@l_pos, @l_values_intersect, @l_values_all, @l_values_intersect_he, @l_values_intersect_ho);

my $chr_length = $chr->length();


my $intSpanSum = Set::IntSpan::Fast::XS->new();
my $intSpanCommon = Set::IntSpan::Fast::XS->new();
my $intSpanDiff = Set::IntSpan::Fast::XS->new();


my $i = 1;
foreach my $var (@{$chr->getListVarObjects($v_sum)}) {
	print '@' if ($i == 0);
	if ($i == 10000) {
		print '.';
		$i = 1;
	}
	my $vector_id = $var->vector_id();
	my $pos = $var->start();
	$intSpanSum->add($pos);
	
	if ($only_ho) {
		if ($v_intersect_ho->contains($vector_id)) {
			$intSpanCommon->add($pos);
		}
		else {
			$intSpanDiff->add($pos);
		}
		
	}
	else {
		if ($v_intersect->contains($vector_id)) {
			$intSpanCommon->add($pos);
		}
		else {
			$intSpanDiff->add($pos);
		}
	}
	
	
	$i++;
}

my $max_value = 0;
my $min_value = 0;
my $pos = 0;
$i = 1;
my $j = 1;
my $n = 0;
my $sum = 0;
my $sum2 = 0;

my (@l_x, @l_y, @l_x2, @l_y2);
my $h_pos_inter;
foreach my $pos ($intSpanSum->as_array()) {
	my $nb_inter = ($pos/$interval) +1;
	#warn $pos.' -> '.$nb_inter;
	if ($intSpanCommon->contains($pos)) {	
		$h_pos_inter->{int($nb_inter)}->{common}++;
	}
	else {
		$h_pos_inter->{int($nb_inter)}->{diff}--;
	}
}

foreach my $pos (sort {$a <=> $b} keys %$h_pos_inter) {
	push(@l_x, $pos);
	if (exists $h_pos_inter->{$pos}->{common}) { push(@l_y, $h_pos_inter->{$pos}->{common}); }
	else { push(@l_y, undef); }
	
	
	push(@l_x2, $pos);
	if (exists $h_pos_inter->{$pos}->{diff}) { push(@l_y2, $h_pos_inter->{$pos}->{diff}); }
	else { push(@l_y2, undef); }
}



my $hashRes;
$hashRes->{'label'} = 'name';
$hashRes->{'x_common'} = \@l_x;
$hashRes->{'y_common'} = \@l_y;
$hashRes->{'x_diff'} = \@l_x2;
$hashRes->{'y_diff'} = \@l_y2;
$hashRes->{'interval'} = $interval;

my $json_encode = encode_json $hashRes;
print ".\",";
$json_encode =~ s/{//;
print $json_encode;
exit(0);


