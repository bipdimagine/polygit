#!/usr/bin/env perl

use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";
use strict;
use warnings;
use Data::Dumper;
use Carp;
use GBuffer;
use JSON::XS; 
use CGI qw/:standard :html3/;


my $project_name;
my $patient_name;
my $fork = 10;
my $force;

my $cgi = new CGI();
my $buffer = GBuffer->new();

my $project_name = $cgi->param("project");
my $force = $cgi->param("force");

my $project = $buffer->newProject( -name => $project_name );
my $out =  $project->getVariationsDir("crispressoWGS");
my $json_file = $out."/offtarget.json";
unlink $json_file if $force;
warn  $json_file;

if (-e $json_file){


open(my $fh, "<", $json_file)

    or die "Impossible d'ouvrir $json_file : $!";

local $/;

my $json = <$fh>;

close($fh);

print $cgi->header(

    -type => 'application/json',

    -charset => 'utf-8'

);
 print $json;
 exit(0);
 
}

#  --flexiguide_seq GCTGAACTTTTATGTCTCTA --flexiguide_homology 80
my $crispresso = "exec_singularity.sh crispresso2.sif  CRISPRessoWGS";
#my $bam = $patient->getAlignmentFile();
my $ref               = $project->genomeFasta();
my $bed = $project->project_path()."/crispresso.bed";



my $input = $project->project_path()."/crispor.tsv";
open(my $IN,  "<", $input)  or die "Cannot open $input: $!";

my $window = 44;
my $half   = int($window / 2);

my $header;
my $annot;
while (my $line = <$IN>) {

    chomp $line;
    next if $line =~ /^\s*$/;
    next if $line =~ /^#/;
    my @f = split(/\t/, $line, -1);

    my $chrom  = $f[7];
    my $start  = $f[8];
    my $end    = $f[9];
    my $strand = $f[10];

    next unless defined $chrom && $chrom ne "";
    next unless $start =~ /^\d+$/ && $end =~ /^\d+$/;
	$annot->{$chrom.$start.$end}->{annot} = $f[-1];
	$annot->{$chrom.$start.$end}->{grna} = $f[1];
	$annot->{$chrom.$start.$end}->{strand} = $strand;
	warn $strand;
}

close $IN;




open (BED ,"$bed");
my $dejavu;
my $controls = {1825=>1 , 1629=>1,1827=>1};

while (my $line = <BED>){
	chomp ($line);
	my ($chr,$start,$end,$pname,@a ) = split(" ",$line);
	$pname =~ s/:/_/;
	my ($c1,$p2) = split("_",$pname);
	my ($s1,$e1) = split("-",$p2);
	# ./CRISPRessoWGS_on_1623/CRISPResso_on_chr18_21612650-21612672/Alleles_frequency_table.zip
	foreach my $p (@{$project->getPatients}){
		my $patient_name = $p->name;
	my $alleles = $out."/CRISPRessoWGS_on_${patient_name}/CRISPResso_on_${pname}/Alleles_frequency_table.zip";
	my @lres= `zgrep -w \"MODIFIED\" $alleles`;
	chomp(@lres);
	if (@lres) {
		#warn Dumper @res;
		#
		foreach my $line (@lres){
		my @res = split(" ",$line);	
		
		#$dejavu->{$res[1]}->{nb} ++;
		$dejavu->{$res[0]}->{chr} =$chr;
		$dejavu->{$res[0]}->{id} =$pname;
		$dejavu->{$res[0]}->{start} =$s1;
		$dejavu->{$res[0]}->{end} =$e1;
		$dejavu->{$res[0]}->{reference} =$res[1];
		$dejavu->{$res[0]}->{sequence} = $res[0];
		my $m = int($res[8]*10);
		$m /= 10;
		
		$dejavu->{$res[0]}->{patients}->{$patient_name} = $res[7]." (".$m."%)";
		
		#warn join("\t",$chr,$start,$end,$res[0],$res[1],$res[6],$res[7],$res[8])."\n";
		}
		
	}
	}
}

close(BED);
my @patients = grep {!(exists $controls->{$_})} map { $_->name } @{$project->getPatients};

my @rows;
#KO 1637 , 1628,1826,1824 
#HE 1323 1833 1624 
#CONTROL: 1825 1629 1827 

foreach my $d (sort { $a->{chr} cmp $b->{chr} ||  $a->{start} cmp $b->{start}} values %$dejavu) {
	my $key = $d->{sequence};
   # my $d = $dejavu->{$key};

    my @patient_values;
	 my $nbp ;
	 my $co;
	 foreach my $c (keys %$controls){
			$co ++ if exists $d->{patients}->{$c};
		}
		next if $co;
    foreach my $p (@patients) {
		
        if (exists $d->{patients}->{$p}) {
			
            push @patient_values, $d->{patients}->{$p};
            $nbp++;
        }
        else {
            push @patient_values, '-';
        }
    }
     #."/".BioTools::complement_sequence($v->sequence());
    unless ($key eq "ATGCTTTGTATCTTTAGAGACATAAAATTACAACCAAGAGGAAA" ){
	#	next ;
		#warn $d->{reference};
		#die();
	}
	#warn $d->{reference};
	#warn $key; 
	my $ref = $d->{reference}."";
	if( $annot->{$d->{chr}.$d->{start}.$d->{end}}->{strand} eq "-"){
		$key = BioTools::complement_sequence($key);
		$ref =  BioTools::complement_sequence($ref);
	}
    push @rows, {
        chr         => $d->{chr},
        start       => $d->{start},
        end         => $d->{end},
        nb_sample   => $nbp,
        patients    => \@patient_values,
        reference   => $ref,
        sequence    => $key,
        annotations => $annot->{$d->{chr}.$d->{start}.$d->{end}}->{annot} // '',
        grna => $annot->{$d->{chr}.$d->{start}.$d->{end}}->{grna} // '',
        genome => lc($project->genome_version_generic),
        strand =>  $annot->{$d->{chr}.$d->{start}.$d->{end}}->{strand},
    };
}


my $data = {
    patients => \@patients,
    rows     => \@rows,
};

open(my $fh, ">", $json_file)

    or die "Impossible d'ouvrir $json_file : $!";

print $fh encode_json($data);

#
#print join("\t","chr","start","end","nb_sample",@patients,"reference","sequence","annotations")."\n";
#for my $key (keys %$dejavu){
#	print $dejavu->{$key}->{chr}."\t".$dejavu->{$key}->{start}."\t".$dejavu->{$key}->{end}."\t".$dejavu->{$key}->{nb}."\t";
#	my @t;
#	foreach my $p (@patients){
#		unless( exists $dejavu->{$key}->{patients}->{$p}){
#			push(@t,'-');
#			next;
#		}
#		push(@t,$dejavu->{$key}->{patients}->{$p});
#	}
#	print join("\t",@t,$dejavu->{$key}->{reference},$key,$annot->{$dejavu->{$key}->{chr}.$dejavu->{$key}->{start}.$dejavu->{$key}->{end}})."\n";
#	
#}