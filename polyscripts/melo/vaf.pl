#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
use Data::Dumper;
use Getopt::Long;
use Carp;
use List::Util qw(sum);

use GBuffer;

my $buffer = new GBuffer;
my $project = $buffer->newProject(-name=>'NGS2024_7797'); # NGS2024_7600
print $project->name() ."\n";

my $calling_methods = $project->getCallingMethods;
# print Dumper $calling_methods;

my $calling_method = $calling_methods->[0];
# print $calling_method ."\n";

my $vcf_dir = $project->getVariationsDir($calling_method);
# print $vcf_dir ."\n";

my $hash;

my $patients = $project->getPatients;
my $vaf;
my $ratio_allele;

foreach my $pat(@$patients) {
#	print $pat->name . "\n";
#	my $vcf_file = $pat->getVariationsFile($calling_method);
	my $vcf_file = '/data-isilon/sequencing/ngs/NGS2024_7797/ref_mutee/' . $pat->name .'/variants/medaka.annotated.vcf.gz';
#	print $vcf_file ."\n";
	
	open(my $vcf, "zcat $vcf_file |") or confess("Can't open file $vcf_file");
	
	while(my $line = readline($vcf)) {
		chomp $line;
		if ($line =~ /#/) { # headers
#			print $line . "\n";
		}

		if (not $line =~ /#/) { # not an header
#			print $line . "\n";
			my ($chr, $pos, $ID_rs, $ref, $alt, $qual, 
				$filter, $info, $format, $patient)  = split("\t", $line);
	
			my $var_ID = $chr . '_' . $pos . '_' . $ref . '_' . $alt;
			$hash->{$var_ID}->{'chr_id'} = $chr;
			$hash->{$var_ID}->{'start'} = $pos;
			$hash->{$var_ID}->{'end'} = $pos + abs( length($alt) - length($ref) ) + 1;  # $pos - length($alt) + length($ref) +1;
			$hash->{$var_ID}->{'ref_all'} = $ref;
			$hash->{$var_ID}->{'var_all'} = $alt;
			
			my @info = split(";", $info);
			foreach my $i (@info){
				my ($key, $val) = split("=", $i);
				$hash->{$var_ID}->{'annex'}->{$pat->id}->{$key} = $val;
			}
			
			my @format = split(":", $format);
			my @patient = split(':', $patient);
			foreach my $i (0 .. scalar @format -1) {
				$hash->{$var_ID}->{'annex'}->{$pat->id}->{$format[$i]} = $patient[$i];
			}
			
#			print Dumper $hash;
			
			$vaf->{$pat->name}->{$var_ID}->{'start'} = $pos;
			$vaf->{$pat->name}->{$var_ID}->{'ref_all'} = $ref;
			$vaf->{$pat->name}->{$var_ID}->{'var_all'} = $alt;
			
			# taux revertants = SR alt1 fwd + alt1 rev / (DP-AR)
			my @SR = split(',', $hash->{$var_ID}->{'annex'}->{$pat->id}->{'SR'});
			my @AR = split(',', $hash->{$var_ID}->{'annex'}->{$pat->id}->{'AR'});
			my $DP = $hash->{$var_ID}->{'annex'}->{$pat->id}->{'DP'};
			
			$vaf->{$pat->name}->{$var_ID}->{'VAF'} = ($SR[2] + $SR[3]) / ($DP - sum @AR); # reads alt (fwd + rev) / (tot - ambiguous reads)
#			die if ($SR_alt / ($DP - $AR) != ($SR[2] + $SR[3]) / ($DP - sum @AR));
#			print Dumper $vaf;


			$ratio_allele->{$pat->name}->{$var_ID}->{'start'} = $pos;
			$ratio_allele->{$pat->name}->{$var_ID}->{'ref_all'} = $ref;
			$ratio_allele->{$pat->name}->{$var_ID}->{'var_all'} = $alt;
			
			my @SC = split(',', $hash->{$var_ID}->{'annex'}->{$pat->id}->{'SC'});
			$ratio_allele->{$pat->name}->{$var_ID}->{'ratio_ref'} = ($SC[0]+$SC[1]) / sum @SC;
			$ratio_allele->{$pat->name}->{$var_ID}->{'ratio_alt'} = ($SC[2]+$SC[3]) / sum @SC;
		}
	}
	
	close($vcf);
}

print "\n";

# VAF
print "VAF:\n";
print "#patient \tpos \tref_all \tvar_all \tvaf\n";
foreach my $pat (sort(keys(%$vaf))) {
	foreach my $var (sort(keys (%{$vaf->{$pat}}))) {
		print $pat ."\t"
			. $vaf->{$pat}->{$var}->{'start'} ."\t" 
			. $vaf->{$pat}->{$var}->{'ref_all'} ."\t"
			. $vaf->{$pat}->{$var}->{'var_all'} ."\t"
			. $vaf->{$pat}->{$var}->{'VAF'} ."\n";
	}
}

print "\n";

# Ratio allélique
print "Ratio allélique:\n";
print "#patient \tpos \tref_all \tvar_all \tratio ref \tratio alt \n";
foreach my $pat (sort (keys %$ratio_allele)) {
	foreach my $var (keys %{$ratio_allele->{$pat}}) {
		print $pat ."\t"
			. $ratio_allele->{$pat}->{$var}->{'start'} ."\t" 
			. $ratio_allele->{$pat}->{$var}->{'ref_all'} ."\t"
			. $ratio_allele->{$pat}->{$var}->{'var_all'} ."\t"
			. $ratio_allele->{$pat}->{$var}->{'ratio_ref'} ."\t"
			. $ratio_allele->{$pat}->{$var}->{'ratio_alt'} ."\n";
	}
}




# DP:	"Depth of reads at pos"
# DPS:	"Depth of reads at pos by strand (fwd, rev)"
# DPSP:	"Depth of reads spanning pos +-25"
# SR:	"Depth of spanning reads by strand which best align to each allele (ref fwd, ref rev, alt1 fwd, alt1 rev, etc.)"
# AR:	"Depth of ambiguous spanning reads by strand which align equally well to all alleles (fwd, rev)"
# SC:	"Total alignment score to each allele of spanning reads by strand (ref fwd, ref rev, alt1 fwd, alt1 rev, etc.) aligned with parasail: match 5, mismatch -4, open 5, extend 3"

# Variant Allele Frequency (VAF) = SR alt1 fwd + alt1 rev / (DP-AR)

# exemple:
# CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  F_CBE_1
# CBE     253     .       G       A       33.384  PASS    AR=4,13;DP=300;DPS=150,150;DPSP=300;SC=35878,34077,36373,34626;SR=46,39,100,98  GT:GQ   1:33
# => VAF = (100+98) / (300-(4+13)) = 0,699646643

