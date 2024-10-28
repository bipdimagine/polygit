#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
use Data::Dumper;

use GBuffer;

my $project_name = 'NGS2015_0794';

my $filename = '/home/mperin/test.ped';
open(my $pedigree, '<', $filename) or die "Could not open file '$filename' $!";


my $buffer = new GBuffer;
warn ref($buffer)."\n";

my $project = $buffer->newProject(-name => $project_name);
warn ref($project)."\n";


my $table1;
my $table2;
my $table3;

while(my $patient = readline($pedigree)){
	chomp $patient; # retire le "\n" à la fin
	(my $family, my $name, my $father, my $mother, my $sex, my $status) = split("\t", $patient);
#	print join(", ", $name, $family, $status, $sex)."\n";
	
	if ($status == 1) {$status = 'sain'}
	else {$status = 'malade'}

	if ($sex == 1){$sex = "garçon"}
	else {$sex = "fille"}
	
	my $relationship;
	if ($father == 0 and $mother == 0) {$relationship = "parent"}
	else {$relationship = "enfant"}

	my $patient = $project->getPatient($name);
	my $bamURL = $patient->getBamFile;
	
	my $infos;
#	@$infos = ($status, $sex, $relationship, $bamURL);
	$infos->{'status'} = $status;
	$infos->{'sex'} = $sex;
	$infos->{'srelationship'} = $relationship;
	$infos->{'bamURL'} = $patient->getBamFile;
	
	my $tab1_family->{$family} = $infos;
	$table1->{$name} = $tab1_family;
	
	my $tab2_patient->{$name}=$infos;
	$table2->{$family}->{$name} = $tab2_patient->{$name};
	
#	@$infos = ($family, $sex, $relationship, $bamURL);
	$infos->{'family'} = $family;
	delete ($infos->{'status'});
	my $tab3_patient->{$name} = $infos;
	$table3->{$status}->{$name} = $tab3_patient->{$name};
}

close($pedigree);

print "table1:\n";
print Dumper $table1;
print "\n\n";

print "table2:\n";
print Dumper $table2;
print "\n\n";

print "table3:\n";
print Dumper $table3;
print "\n\n";


