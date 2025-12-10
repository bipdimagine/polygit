#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);

#use lib "/software/polyweb/poly-src/GenBo/lib/obj-nodb/";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";

use GBuffer;

#use lib "/bip-d/soft/distrib/tabix/latest/perl";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/coverage";
use lib "$Bin/../packages/validation_variation";
use lib "$Bin/../packages/cache";

use lib "$Bin/../GenBo/lib/obj-nodb/packages";


use draw_cnv;
use Carp;
use strict;
use Data::Dumper;
use GenBo;
use validationQuery;
use QueryValidationAcmg;
use Date::Tiny;
use JSON::XS;

my $cgi = new CGI();


my $buffer = GBuffer->new();		

# recupere les options  (ici -project echantillon_condition1 echantillon_condition2) 
my $project_name = $cgi->param('project');
my $patient_name = $cgi->param('patient');
$patient_name = $cgi->param('patients') unless $patient_name;
my $chr_name_cgi = $cgi->param('chr');

my $coef = 0.5;
my $mod = 40;
die("\n\nNo -project option... Die...\n\n") unless ($project_name);
my $chr_name =  $chr_name_cgi;
$chr_name = "X" if ($chr_name_cgi == 23);
$chr_name = "Y" if ($chr_name_cgi == 24);
#  instanciation d'un objet project_cache et des objets patient et chromosome
my $project = $buffer->newProjectCache( -name => $project_name, -typeFilters=>"");
$coef = 0.15 if $project->isDiagnostic;
$mod = 1 if $project->isDiagnostic;
my $patient = $project->getPatient("$patient_name");
my $chr = $project->getChromosome($chr_name);
my ($fx,$fy) = get_dude_tab($patient);
my $hres = {};
$hres->{items}->{FX} = $fx;
$hres->{items}->{FY} = $fy;
$hres->{items}->{CHR} = $chr_name_cgi;
$hres->{items}->{PATIENT} = $chr_name_cgi;
$hres->{items}->{NO_DATA} = 1 unless $fx; 

print $cgi->header('text/json-comment-filtered');
print encode_json $hres;#[$hres];
exit(0);

sub get_dude_tab {
	my ($patient) = @_;
	my $primers = $chr->getPrimers();
	return unless @$primers;
	my $data =[] ;
	my $data_mother=[];
	my $data_father=[];
	my $positions;
	
	my $fx ={};
	my $fy ={};
	
	my $father = $patient->getFamily->getFather();
	my $mother = $patient->getFamily->getMother();
	foreach my $p (@$primers){
		#if ($p->cnv_score($patient) == -1 or $p->level($patient) == -1){
		#	push(@$positions,$p->start);
		#	push(@$data,log2(1));
		#}
		
		next if $p->cnv_score($patient) == -1;
		next if $p->level($patient) == -1;
		next if $p->level($patient) == -3;
		my $genes = $p->getGenes();
		next unless @$genes;
		next unless @{$genes->[0]->getMainTranscripts()};
		my $v = log2($p->cnv_score($patient));
		unless($project->isDiagnostic){
		if ($v<-.75 or $v > 0.7 ){
			my $st;
			my $stats= new Statistics::Descriptive::Full;
			my $score =0;
			next if $p->sd($patient) > 0.6;
			next if $p->statistics->quantile(1)< 0.7;
			 next if $p->statistics->quantile(3)> 1.3;
		}
		}
		push(@$positions,$p->start);
		
		push(@$data,$v);
		push(@$data_mother,log2($p->cnv_score($mother))) if $mother; ;
	 	push(@$data_father,log2($p->cnv_score($father))) if $father;
	}
	return if @$data < 3;
	#return 
  my $smoothed_data = smooth_data($data,$positions);
  my $smoothed_data_father = smooth_data($data_father,$positions)  if $father;
  my $smoothed_data_mother = smooth_data($data_mother,$positions) if $mother;
   my $unsmooth = smooth_data_3($data);
 # my $unsmooth = {};
	#warn Dumper $data_smoothed1;
	my $z =0;
	my $red =0;
	my $blue = 0;

	for (my $i=0;$i<@$positions;$i++){
		if ($smoothed_data->[$i] > -0.3 && $smoothed_data->[$i] < 0.3){
			$z++;
			next unless $z%$mod == 0;
		}
		push(@{$fx->{all}},$positions->[$i]);
		
		my $data_local = $smoothed_data->[$i];
		my $data_local_father = $smoothed_data_father->[$i];
		my $data_local_mother = $smoothed_data_mother->[$i];
		if (exists $unsmooth->{$i}){
			$data_local = $data->[$i];
			$data_local_mother = $data_mother->[$i];
			$data_local_father = $data_father->[$i];
		}
		$data_local =  sprintf("%.2f", $data_local);
		
		push(@{$fy->{all}},$data_local);
		
		if ($data_local > 0.6){
			push(@{$fx->{dup}},$positions->[$i]);
			push(@{$fy->{dup}},$data_local);
			push(@{$fy->{dup_father}},$data_local_father);
			push(@{$fy->{dup_mother}},$data_local_mother);
			$blue++;
		}
		elsif ($data_local < -0.75){
			push(@{$fx->{del}},$positions->[$i]);
			push(@{$fy->{del}},$data_local);
			push(@{$fy->{del_father}},$data_local_father);
			push(@{$fy->{del_mother}},$data_local_mother);
			$red++;
		}
		else {
			push(@{$fy->{all_father}},$data_local_father);
			push(@{$fy->{all_mother}},$data_local_mother);
		}
	}
	$fy->{red} = $red;
	$fy->{blue} = $blue;
	unless ($mother){
		$fy->{all_mother} = [] ;
		$fy->{del_mother}= [] ;
		$fy->{dup_mother}= [] ;
	}
	unless ($father){
		$fy->{all_father} = [] ;
		$fy->{del_father}= [] ;
		$fy->{dup_father}= [] ;
	}
	return ($fx,$fy);
	
}



sub smooth_data {
	my ($data,$positions) = @_;
	my $smoother = Statistics::Descriptive::Smoother->instantiate({
         method   => 'weightedexponential',
         coeff    => $coef,
         data     => $data,
         samples  => $positions,
  });
  my @smoothed_data = $smoother->get_smoothed_data();
  return \@smoothed_data 
}
sub smooth_data_3 {
	my ($data) = shift;
	my $res;
	for (my $i = 5;$i<scalar(@$data)-5;$i++){
		next if $data->[$i] > 0.6; 
		next if  $data->[$i] < 0.3; 
		my $mean = ($data->[$i-5]+$data->[$i-4]+$data->[$i-3]+$data->[$i-2]+$data->[$i-1])/5;
		
		my $mean2 = ($data->[$i+5]+$data->[$i+4]+$data->[$i+3]+$data->[$i+2]+$data->[$i+1])/5;
		if (abs($mean-$mean2) < 0.15 && ($mean -1) < 0.1  && ($mean2 -1) < 0.1){
			$res->{$i}++;
		} 
	}
	return $res;
}

sub log2 {
	 my $n = shift;
	 return -1.5 if $n <0.01;
	 my $v = log($n)/log(2);
	 $v = -1.5 if $v < -1.5; 
	 
    return $v;
}

sub printJson {
	my ($listHash) = @_;
	my $hash;
	my @t;
	
	#@t = sort { $a->{POS} <=> $b->{POS}} @$listHash;
	$hash->{'identifier'} = 'PLOT';
	$hash->{'items'} = \@$listHash;
	
	print $cgi->header('text/json-comment-filtered');
	print encode_json $hash;
	print "\n";
}

