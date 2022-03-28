#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
#use Set::IntSpan;
use GBuffer;
use Data::Dumper;
use Getopt::Long;
use Carp;
use Bio::DB::HTS;
use Bio::DB::Sam;
use colored;



my $file2;
my $file1;
my $patient_name;
my $buffer1 = GBuffer->new();
my $run_id;
GetOptions(
	'run=s'   => \$run_id,
	
);

my $query = $buffer1->getQuery();
my $res = $query->getAllPatientsFromRunId($run_id);
foreach my $r (@$res){
my $buffer = GBuffer->new();	
my $error;

my $project = $buffer->newProject( -name 			=> $r->{project});	

my $patient = $project->getPatient($r->{patient});

my $file1 = $patient->getBamFile();
my $hfile        = Bio::DB::HTSfile->open($file1);                            
my $header       = $hfile->header_read;
 $patient_name = $patient->name();
warn $r->{patient};
warn $file1;
my $PG = parse_header_PG($header);
colored::stabilo("blue","NO elprep $patient_name !!!!") unless exists $PG->{elprep};
 next unless exists $PG->{elprep};
my $RG = parse_header_RG($header);

warn "test RG ".$RG->{SM}." <=> ".$patient_name ;

if ($RG->{SM} ne $patient_name){
	 colored::stabilo("red","problem RG $patient_name ".Dumper($RG));
	 $error ++ ;
}

warn "test elprep bam" ;
my $helprep = $PG->{elprep};
die() unless $helprep;
 colored::stabilo("red","probnlem read group elprep ")if $helprep->{CL} !~/SM:$patient_name/;;
 $error ++ if $helprep->{CL} !~/SM:$patient_name/;


my @bams = grep{/\.bam/} split(" ",$helprep->{CL});
foreach my $f (@bams){
	my @ff = split("/",$f);
	 colored::stabilo("red","ERROR ") if $ff[-1] !~/$patient_name/;
	 $error ++ if $helprep->{CL} !~/SM:$patient_name/;
	 
}




warn "test fastq" ;

my $hbwa = $PG->{bwa};
die() unless $hbwa;

my @fastq = grep{/fastq/} split(" ",$hbwa->{CL});

foreach my $f (@fastq){
	my @ff = split("/",$f);
	 colored::stabilo("red","ERROR ") if $ff[-1] !~/$patient_name/;
 	die("problem $file1 \n $patient_name <=> ". $ff[-1]) if $ff[-1] !~/$patient_name/;
}

colored::stabilo("green","OK $patient_name !!!!") unless $error;
colored::stabilo("red","ERROR  $patient_name !!!!") if $error;

}

exit(0);



sub parse_header_PG {
	my ($header) = @_;
	my @lines = split("\n",$header->text);
	my $hash;
	foreach my $l (@lines){
		chomp($l);
		my ($id,@b) = split("\t",$l);
		next if $id ne "\@PG";
		my $h1;
		foreach my $bb (@b){
			my ($k,@v) = split(":",$bb);
			$h1->{$k} = join(":",@v);
		}
	my $pname = $h1->{PN};
	$pname = $h1->{ID} unless $pname;
	die() unless $pname;
	$hash->{$pname} = $h1;
	}
	return $hash;
}
	
sub parse_header_RG {
	my ($header) = @_;
	my @lines = split("\n",$header->text);
	my $hash;
	foreach my $l (@lines){
		chomp($l);
		my ($id,@b) = split(" ",$l);
		next if $id ne "\@RG";
		my $h1;
		foreach my $bb (@b){
			my ($k,@v) = split(":",$bb);
			$h1->{$k} = join(":",@v);
		}
		return $h1;
	}
}