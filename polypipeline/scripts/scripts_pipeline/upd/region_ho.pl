#!/usr/bin/env perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
use GBuffer;
use Getopt::Long;
use JSON::XS;
use MCE::Loop;
use MCE::Shared;

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

 	my $regions = MCE::Shared->array();
#warn $project_name;
my $fam = $patient->getFamily();
my $mother = $fam->getMother();
my $father = $fam->getFather();

my $all_res;
#unless ($fam->isTrio){

	
	MCE::Loop::init {
    chunk_size => '3',  # chaque worker recevra un seul élément
    max_workers => '10', 

};
$project->getChromosomes;
$project->disconnect;
mce_loop {
    my ($mce, $rall,$chunk_id) = @_;
   	warn $chunk_id;
   foreach my $chr (@$rall){
	next if $chr->name eq 'X' && $patient->isMale();
	next if $chr->name eq 'Y';
	next if $chr->name eq 'MT';
   my($nb,$start,$end) = max2($patient,$chr);
   if ($nb > 5000){
	 push(@$regions,{chr=>$chr->name,start=>$start,end=>$end,nb=>$nb});
	}
  
   
	}
	# MCE->gather($chunk_id,$hash);
} @{$project->getChromosomes};

my $list =[];
foreach my $a (@$regions){
	push(@$list,$a);
}
save_json($patient,$list);
#}
sub save_json{
my ($patient,$regions) =@_;
my $file = $patient->upd_file().".ho";

warn $file;
open(JSON,">$file");
print JSON  encode_json $regions;
close JSON;
exit(0);
}

sub max2 {
	my ($patient,$chr)  = @_;
	my $project = $patient->project;
	my $parquet = $project->parquet_cache_variants();
	my $chr_name = $chr->ucsc_name;
	my $threshold = 0.97;
	my $table_patient = "patient_".$patient->id."_type";
	my $cp2 = "patient_".$patient->id."_type";
	#CAST(patient_55661_transmission AS INTEGER)
	my $asql_patient = [];
	my $asql_patient_only_transmission = [];
	my $t = time;
	my $sql = qq{ SELECT  variant_start,$table_patient  FROM '$parquet' a where  ${table_patient} <> 0 and variant_type=1 and variant_chromosome='${chr_name}' order by variant_start;};
	
	my $cmd = qq{duckdb -column -noheader -c "$sql"};
	open my $fh, '-|', $cmd or die "Cannot execute $cmd: $!";
	my @snp;
while (my $line = <$fh>) {

	
		chomp($line);
	 	my ($pos, $value) = split(" ",$line);
	 	next if $chr->blacklist->contains($pos);
	    push @snp, [$pos, $value];

	}
close($fh);	
# ------------------------------------------------------------

# Recherche de la plus longue région >= 95 % homozygote

# 2 = homozygote

# 1 = hétérozygote

# ------------------------------------------------------------

my $left = 0;

my $n_hom = 0;

my $n_het = 0;

my ($best_left, $best_right) = (0, 0);

my $best_length = 0;

for my $right (0 .. $#snp) {

    my ($pos, $value) = @{$snp[$right]};

    if ($value == 2) {

        $n_hom++;

    }

    elsif ($value == 1) {

        $n_het++;

    }

    # Tant que la proportion d'homozygotes est < 95 %,

    # on raccourcit la fenêtre par la gauche.

    while ($left <= $right) {

        my $length = $right - $left + 1;

        my $fraction = $n_hom / $length;

        last if $fraction >= $threshold;

        my ($left_pos, $left_value) = @{$snp[$left]};

        if ($left_value == 2) {

            $n_hom--;

        }

        elsif ($left_value == 1) {

            $n_het--;

        }

        $left++;

    }

    my $length = $right - $left + 1;

    if ($length > $best_length) {

        $best_length = $length;

        $best_left  = $left;

        $best_right = $right;

    }

}

# ------------------------------------------------------------

# Résultat

# ------------------------------------------------------------

my $start = $snp[$best_left][0];

my $end   = $snp[$best_right][0];

my $hom = 0;

my $het = 0;

for my $i ($best_left .. $best_right) {

    if ($snp[$i][1] == 2) {

        $hom++;

    }

    elsif ($snp[$i][1] == 1) {

        $het++;

    }

}
$best_length = 0.0001 if $best_length ==0;
my $fraction = $hom / $best_length;

return($hom,$start,$end);

print "===================== $chr_name ===================\n";

print "Plus longue région >= ", $threshold * 100,

      "% homozygote\n";

print "========================================\n";

print "Start       : $start\n";

print "End         : $end\n";

print "Length      : ", $end - $start + 1, " bp\n";

print "Nb SNP      : $best_length\n";

print "Nb homo     : $hom\n";

print "Nb hetero   : $het\n";

printf "Fraction hom: %.2f%%\n", $fraction * 100;

print "========================================\n";

print "========================================\n";
	warn $cmd;
	
}
