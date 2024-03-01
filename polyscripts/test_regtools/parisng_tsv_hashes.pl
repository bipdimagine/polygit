#!/usr/bin/perl
use strict;
use warnings;
use Text::CSV;
use Data::Dumper;


my $input_file = '/home/mhamici/parsing_regtools/results/chr1.tsv';

my $csv = Text::CSV->new({ sep_char => "\t", auto_diag => 1, binary => 1 });

open my $fh, '<', $input_file or die "Cannot open $input_file: $!";
my $header = <$fh>;
chomp $header;
$csv->parse($header);

my @column_names = $csv->fields();

my %DA_pos;
my %DA_neg;
my %NDA_pos;
my %NDA_neg;
my %D_pos;
my %D_neg;
my %A_pos;
my %A_neg;
my %N_pos;
my %N_neg;

# Process each line and create a hash for each line
while (<$fh>) {
    chomp;
    $csv->parse($_);
    my %data_hash;
    @data_hash{@column_names} = $csv->fields();

    # Use the value in the "name" column as the key for the hash of hashes
    my $key = $data_hash{'name'};
    
    # Filter the lines by anchor and strand into different hashes
    if ($data_hash{"anchor"} eq "DA" && $data_hash{"strand"} eq "+"){
    	$DA_pos{$key} = \%data_hash;
    } elsif ($data_hash{"anchor"} eq "DA" && $data_hash{"strand"} eq "-") {
    	$DA_neg{$key} = \%data_hash;
    } elsif ($data_hash{"anchor"} eq "NDA" && $data_hash{"strand"} eq "+") {
    	$NDA_pos{$key} = \%data_hash;
    } elsif ($data_hash{"anchor"} eq "NDA" && $data_hash{"strand"} eq "-") {
    	$NDA_neg{$key} = \%data_hash;
    } elsif ($data_hash{"anchor"} eq "D" && $data_hash{"strand"} eq "+") {
    	$D_pos{$key} = \%data_hash; 	
    } elsif ($data_hash{"anchor"} eq "D" && $data_hash{"strand"} eq "-") {
    	$D_neg{$key} = \%data_hash;
    } elsif ($data_hash{"anchor"} eq "A" && $data_hash{"strand"} eq "+") {
    	$A_pos{$key} = \%data_hash;
    } elsif ($data_hash{"anchor"} eq "A" && $data_hash{"strand"} eq "-") {
    	$A_neg{$key} = \%data_hash;
    } elsif ($data_hash{"anchor"} eq "N" && $data_hash{"strand"} eq "+") {
    	$N_pos{$key} = \%data_hash;
    } elsif ($data_hash{"anchor"} eq "N" && $data_hash{"strand"} eq "-") {
    	$N_neg{$key} = \%data_hash;
    }

}

close $fh;

#Find the DA to the D+ junctions
for my $D_pos_junction (keys %D_pos) {
	my $start = $D_pos{$D_pos_junction}{"start"};
	
	my @matching_junctions = grep { $DA_pos{$_}{"start"} == $start } keys %DA_pos;
	if (@matching_junctions) {
		print "The D_pos junction $D_pos_junction that starts with $start have the following matching DA junction(s)";
		print join (', ',@matching_junctions);
		print "\n";
	} else {
		print "The D_pos junction $D_pos_junction that starts with $start has no matching DA junction \n";
	}
}

#Find the DA to the D- junctions
for my $D_neg_junction (keys %D_neg) {
	my $end = $D_neg{$D_neg_junction}{"end"};
	my @matching_junctions = grep { $DA_neg{$_}{"end"} == $end } keys %DA_neg;
	if (@matching_junctions) {
		print "The D_neg junction $D_neg_junction that ends with $end have the following matching DA junction(s)";
		print join (', ',@matching_junctions);
		print "\n";
	} else {
		print "The D_neg junction $D_neg_junction that ends with $end has no matching DA junction \n";
	}
}

#Find the DA to the A+ junctions
for my $A_pos_junction (keys %A_pos) {
	my $end = $A_pos{$A_pos_junction}{"end"};
	my @matching_junctions = grep { $DA_pos{$_}{"end"} == $end } keys %DA_pos;
	if (@matching_junctions) {
		print "The A_pos junction $A_pos_junction that ends with $end have the following matching DA junction(s)";
		print join (', ',@matching_junctions);
		print "\n";
	} else {
		print "The A_pos junction $A_pos_junction that ends with $end has no matching DA junction \n";
	}
}

#Find the DA to the A- junctions
for my $A_neg_junction (keys %A_neg) {
	my $start = $A_neg{$A_neg_junction}{"start"};
	my @matching_junctions = grep { $DA_neg{$_}{"start"} == $start } keys %DA_neg;
	if (@matching_junctions) {
		print "The A_neg junction $A_neg_junction that startss with $start have the following matching DA junction(s)";
		print join (', ',@matching_junctions);
		print "\n";
	} else {
		print "The A_neg junction $A_neg_junction that starts with $start has no matching DA junction \n";
	}
}

#Find the DA to the NDA+ junctions ### voir si c'est vraiment nécessaire de faire pos et neg séparés
for my $NDA_pos_junction (keys %NDA_pos) {
	my $start = $NDA_pos{$NDA_pos_junction}{"start"};
	my $end = $NDA_pos{$NDA_pos_junction}{"end"};
	my @matching_junctions = grep { 
		$DA_pos{$_}{"start"} >= $start &&
		$DA_pos{$_}{"end"}<= $end 
	} keys %DA_pos;
	if (@matching_junctions) {
		print "The NDA_pos junction $NDA_pos_junction that starts $start and ends with $end have the following DA junction(s)";
		print join (', ',@matching_junctions);
		print "\n";
	} else {
		print "The NDA_pos junction $NDA_pos_junction that starts $start and ends with $end has no matching DA junction \n";
	}	
}

#Find the DA to the NDA- junctions ### voir si c'est vraiment nécessaire de faire pos et neg séparés
for my $NDA_neg_junction (keys %NDA_neg) {
	my $start = $NDA_neg{$NDA_neg_junction}{"start"};
	my $end = $NDA_neg{$NDA_neg_junction}{"end"};
	my @matching_junctions = grep { 
		$DA_neg{$_}{"start"} >= $start &&
		$DA_neg{$_}{"end"}<= $end 
	} keys %DA_neg;
	if (@matching_junctions) {
		print "The NDA_neg junction $NDA_neg_junction that starts $start and ends with $end have the following DA junction(s)";
		print join (', ',@matching_junctions);
		print "\n";
	} else {
		print "The NDA_neg junction $NDA_neg_junction that starts $start and ends with $end has no matching DA junction \n";
	}	
}

#Find the reverse DA to the N+ junctions  ### quoi faire avec less else?
for my $N_pos_junction (keys %N_pos) {
	my $start = $N_pos{$N_pos_junction}{"start"};
	my $end = $N_pos{$N_pos_junction}{"end"};
	my @matching_junctions = grep {
		$DA_neg{$_}{"start"} == $start &&
		$DA_neg{$_}{"end"} == $end
	} keys %DA_neg;
	if (@matching_junctions) {
		print "The N_pos junction $N_pos_junction that starts with $start and ends with $end has the following inverted DA- junction(s)";
		print join (', ',@matching_junctions);
		print "\n";
	} else {
		print "The N_pos junction $N_pos_junction that starts with $start and ends with $end has no matching inverted DA- junction \n";
	}
}

#Find the reverse DA to the N- junctions  ### quoi faire avec less else?
for my $N_neg_junction (keys %N_neg) {
	my $start = $N_neg{$N_neg_junction}{"start"};
	my $end = $N_neg{$N_neg_junction}{"end"};
	my @matching_junctions = grep {
		$DA_pos{$_}{"start"} == $start &&
		$DA_pos{$_}{"end"} == $end
	} keys %DA_pos;
	if (@matching_junctions) {
		print "The N_neg junction $N_neg_junction that starts with $start and ends with $end has the following inverted DA+ junction(s)";
		print join (', ',@matching_junctions);
		print "\n";
	} else {
		print "The N_neg junction $N_neg_junction that starts with $start and ends with $end has no matching inverted DA+ junction \n";
	}
}