#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/kyoto";

#use lib "/bip-d/soft/distrib/tabix/latest/perl";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/coverage";
use JSON::XS;
#use Set::;
use Carp;
use export_data;
use strict;
use Set::IntSpan::Fast::XS;
use Data::Dumper;
use GBuffer;
use Getopt::Long;
use Carp;
use Set::Intersection;
use Tabix;
use Storable qw/thaw/;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use coverage;
use Spreadsheet::WriteExcel;
use POSIX;
use validationQuery;
use Date::Tiny;
use lib "$Bin/../packages/validation_variation"; 
use draw_cnv; 

	


my $nb_gene_by_patient = 3;
my $nb_exon_by_genes = 10;

my $buffer = GBuffer->new();

my $cgi          = new CGI();

my $project_name = $cgi->param('project');

my $project = $buffer->newProject(-name=>$project_name);
my $captures = $project->getCaptures();
my $capture = $captures->[0];
my $capture = $project->getCaptures()->[0];



my $vquery;
my $seq_type = $cgi->param('exome');
my $max = 50;


if ($project->isDiagnostic){
	$max = 75;
	$vquery = validationQuery->new(dbh=>$buffer->dbh,capture_name=>$capture->validation_db());
}

my $patient_name = $cgi->param('patients');
my $gene =  $cgi->param('gene');
my $only_red =1 if $cgi->param('only_red');
my $only_orange =1 if $cgi->param('only_orange');
$patient_name ="all" unless $patient_name;
my $patients = $project->get_list_patients($patient_name,",");


my $transcript_id =  $cgi->param('transcript');
my $transcript =  $project->newTranscript($transcript_id) ;



my $splice_5 = 10;
my $limit    = $cgi->param('limit') +0;


my %dejavu_transcripts;
my $res;
my @exons = sort{$a->start*$a->strand <=> $b->end*$b->strand } @{$transcript->getExons()};

my $xx =0;
my $nb=0;

foreach my $exon (@exons){
			my $line;
			$line->{id} = $exon->name()."";
				$line->{transcript} = $transcript_id."";
			my $zz =0;
		foreach my $patient (sort{$a->name <=> $b->name} @{$patients}){
			my $name = $patient->name();
			my ($mean,$intspan,$min) = $exon->mean_intspan_coverage_coding(patient=>$patient,padding=>$splice_5,limit=>$limit);
			$line->{$name} = int($mean);
			my $hjson;
			$hjson->{mean} = int($mean);
			$hjson->{min} = int($min);
			$hjson->{flag} = 0;
			if ($min< $limit) {
					$hjson->{flag} = -1;
					$nb++;
			}
			$line->{$name} =encode_json $hjson;
			
		}
		push(@$res,$line);
}
$res->[0]->{flag} = 1;
$res->[0]->{flag} = 0 if ($nb == 0);
warn $res->[0]->{flag} ;
my $nb_patient = scalar(@{$patients});

capture_problem($res) if $only_red;
coverage_problem($res) if $only_orange;
export_data::print($project,$cgi,$res);
exit(0);

sub capture_problem {
	my ($res) = @_;
	my $filtered_res= [];
foreach my $line (@$res){
	my $nb =0;
	foreach my $patient (sort{$a->name <=> $b->name} @{$patients}){
		my $name = $patient->name();
		my $z = decode_json $line->{$name};
		$nb  ++ if $z->{flag} == -1;
	}
	next unless ($nb >= $nb_patient*0.9);
	push(@$filtered_res,$line);
}
unless (scalar(@$filtered_res)){
	add_empty_line($filtered_res);
}
export_data::print($project,$cgi,$filtered_res);
exit(0);
}

sub coverage_problem {
	my ($res) = @_;
	my $filtered_res=[];
foreach my $line (@$res){
	my $nb =0;
	foreach my $patient (sort{$a->name <=> $b->name} @{$patients}){
		my $name = $patient->name();
		my $z = decode_json $line->{$name};
		$nb  ++ if $z->{flag} == -1;
	}
	next if ($nb == 0);
	next if ($nb >= $nb_patient*0.9);
	warn $nb;
	push(@$filtered_res,$line);
}
unless (scalar(@$filtered_res)){
	add_empty_line($filtered_res);
}
export_data::print($project,$cgi,$filtered_res);
exit(0);
}

sub add_empty_line{
	my ($filtered_res) = @_;
		my $line;
		$line->{id} = "none";
			$line->{transcript} = $transcript_id."";
	foreach my $patient (sort{$a->name <=> $b->name} @{$patients}){
			my $name = $patient->name();
		
			$line->{$name} ="ex1";
			my $hjson;
			$hjson->{mean} ="..";
			$hjson->{min} = "..";
			$hjson->{flag} = 0;
			$line->{$name} =encode_json $hjson;
		
		}
			push(@$filtered_res,$line);
}
