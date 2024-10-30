#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/.;/../../GenBo/lib/obj-nodb/";
use Data::Dumper;
use Getopt::Long; # librairie pour récuper les options depuis le terminal
use Carp; # librairie pour `confess`

use GBuffer;

my $filename;
GetOptions(
	'file=s'    => \$filename,
);

# comme `die` mais avec plus d'infos
confess("\n\nERROR: -file mandatory. Die.\n\n") unless ($filename);


# @list
# $list pareil que \@list

#my $project_name = 'NGS2015_0794';

open(my $vcf, '<', $filename) or die "Could not open file '$filename' $!";

my $table;

while(my $line = readline($vcf)){
	chomp $line;
	
#	if ($line =~ /^#[^#]/) { # test si $line commence par un seul # (pas par ##)
#		print $line, "\n";
#	}

# regex:
# /NGS20.+/
	
	if (not $line =~ /#/) { # si $line ne commence pas par #
		my ($chr, $pos, $ID_rs, $ref, $alt, $qual, 
			$filter, $info, $format, $patient)  = split("\t", $line);

		my @format = split(":", $format);
		my @patient = split(':', $patient);
		
		my $ID = substr($chr, 3) . '_' . $pos . '_' . $ref . '_' . $alt;
		
		# récupère la position de GT et AD
		my $GT_index, my $AD_index;
		foreach my $i (0 .. $#format){
			if ($format[$i] eq 'GT') {$GT_index = $i}
			if ($format[$i] eq 'AD') {$AD_index = $i}
		}
		# Détermine si le patient est ho ou he
		my $GT = $patient[$GT_index];
		my $he_ho;
		if ($GT eq '0/0') {$he_ho = 'ho ref'}
		elsif ($GT eq '0/1') {$he_ho = 'he'}
		elsif ($GT eq '1/1') {$he_ho = 'ho var'}
		$table->{$ID}->{'he/ho'} = $he_ho;
		
		# Détermine le nombre d'allèles de reférence et variant
		my $AD = $patient[$AD_index];
		(my $nb_ref, my $nb_alt) = split(',', $AD);
		$table->{$ID}->{'nb_ref'} = $nb_ref;
		$table->{$ID}->{'nb_alt'} = $nb_alt;
		
		# Autres infos
		foreach my $i (0 .. $#format) {
			$table->{$ID}->{'details'}->{$format[$i]} = $patient[$i];
			
#		print $ID . "\t";
#		print $GT . "\t";
#		print $he_ho . "\t\t";
#		print $nb_ref . "\t" . $nb_alt ."\n";
		}
	}
}

close($vcf);
print(Dumper $table);
print("\n\n\n");

