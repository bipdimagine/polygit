#!/usr/bin/perl
use CGI;
use strict;
 
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/GenBoDB/writeDB";
use lib "$Bin/../packages/electro";

#use ABI;  
use GenBoCacheQuery;  
use electro;  
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use Carp;
use JSON::XS;
use util_file;
use Storable qw/retrieve/;
use Bio::DB::Sam;
use Set::IntSpan::Fast::XS ;
my $cgi    = new CGI;

#ce script est utilise pour tracer les electrophoregrammes des sequences pour le sequencage classique ou les histogrammes pour les puces de resequencage

#chargement du buffer
my $buffer = GBuffer->new;

use Data::Dumper;


my @bases = ( "A", "T", "C", "G" );
my %tab;
my @mut;
my %traces;


#recuperation des arguments passes au cgi

my $projectName = $cgi->param('project');
my $variation_id =  $cgi->param('variation_id');
my $poly_id =  $cgi->param('poly_id');
die() unless $poly_id;
my $patient_name =  $cgi->param('patient_name');


my $med = 5;


#si projet est de type ngs alors on ne va pas plus loin
# pour l'instant du reste je ne sais pas pourquoi 


#si les donneees ne sont pas dans le cache

my $project;
if ($projectName eq $patient_name and $patient_name =~ /^[0-9]+_[0-9]+$/) {
	($projectName, $patient_name) = split('_', $projectName);
}
unless ($projectName =~ /[A-Za-z]/) {
	my $projId = $projectName;
	my $patId  = $patient_name;
	$patient_name = undef;
	$projectName = getProjectNameFromId($projId);
	$project = $buffer->newProject(-name=>$projectName);
	foreach my $patObj (@{$project->getPatients()}) { if ($patObj->id() eq $patId) { $patient_name = $patObj->name(); } }
	#warn "### PROJECT $projectName and PATIENT $patient_name";
}
else { $project = $buffer->newProject(-name=>$projectName); }
die ("unknown project ".$projectName ) unless $project;
die ("unknown patient ".$projectName ) unless $patient_name;
my $variation;
my $contig;
my $trace;
#


my $data;
		
		my $data2 = getCoverAlignment($poly_id);
		$data2->{type} = 3;
		push(@$data,$data2);

printJson($data);

exit(0);
 
 sub getCoverAlignment {
 	my ($poly_id) = @_;
 	my ($chr,$pos,$all1,$all2) = split("_",$poly_id);
 	$chr =~ s/chr//;
 	my $chr = $project->getChromosome($chr);;
 	my $patient = $project->getPatient($patient_name);
 	

	my $bam_file = $patient->getBamFile();

	return unless -e $bam_file;	
	

	my $ll = 10;
	
	my $dir_g = $project->buffer->{config}->{public_data}->{$project->{golden_path}};
	
	my $sam = Bio::DB::Sam->new(-bam  =>$bam_file,
								-expand_flags => 1,
								#-fasta => $dir_g.'/genome/fasta/all.fa.seq',
								);
							
	my @alignments1 = $sam->get_features_by_location(-seq_id => $chr->ucsc_name,
                                                 -start  => $pos,
                                                 -end    => $pos);
  
 	unless (@alignments1){
 		
 		@alignments1 = $sam->get_features_by_location(-seq_id => $chr->name,
                                                 -start  => $pos,
                                                 -end    => $pos);
 	}
   
my %aligns;
my $left = 10;
my $right = 20;
my %inserts;
my %dejavu;

my @alignments;
my (@find) = grep {$_->score > 30} @alignments1;
@alignments1 = @find if (scalar(@find) >= 50);


my @seq_align;
my $seq_ref = $chr->sequence($pos-$left,$pos+$left);
my $nb =0;

foreach my $aa (sort{$a->start <=> $b->start} @alignments1){
	my $seq = $aa->query->seq->seq;
	next if $dejavu{$seq};
	$nb ++;
	my ($ref,$matches,$query) = $aa->padded_alignment;
	my @tab_query = split("",$query);
	my $temp_start = $aa->start();
	
	my $cigars = $aa->cigar_array();
	shift(@$cigars) if ($cigars->[0]->[0] eq "H" );
 	if ($cigars->[0]->[0] eq "S" ){
 		$temp_start = $temp_start- $cigars->[0]->[1] ;
 	
 	}
 	
 	
 	if ($aa->cigar_str =~ /I/){
 		my $pos_ins = 0;
 		my $seq_ins = $query;
 		my $i_intspan =  Set::IntSpan::Fast::XS->new();
 		foreach my $cig (@$cigars){
 			if ($cig->[0] eq "I"){
 				$i_intspan->add_range($pos_ins,($pos_ins+$cig->[1])-1);
 			}
 			$pos_ins += $cig->[1];
 		}
 	my @seq_ins_before = @tab_query;
 	my $iter = $i_intspan->iterate_runs();  
    while (my ( $from, $to ) = $iter->()) {
    	my $ins = join("",@seq_ins_before[$from..$to]);
      	$seq_ins_before[$from-1] = $seq_ins_before[$from-1].lc($ins);
    }  
 
    my $n_intspan = Set::IntSpan::Fast::XS->new();
    $n_intspan->add_range(0,length $query);
	my $diff = $n_intspan->diff($i_intspan);
	
	my @tdiff =  $diff->as_array();
	@tab_query = @seq_ins_before[@tdiff];
	foreach (my $i=0;$i<@tab_query;$i++){
		$tab_query[$i] = " " if $tab_query[$i] eq undef;
	}
 	} 
 	my $len = $pos - $temp_start ;
 	$len -= $left;
 	if ($len < 0 ){
 		my $sharp = " ";
 		for (my $s=0;$s<abs($len);$s++){
 			unshift(@tab_query,$sharp);
 			$query = $sharp.$query;
 		}
 		
 		$len=0;
 	}
	for  (my $sp =0;$sp<$left*2;$sp++){
 		push(@tab_query," ");
 	}
 	my @cut_query = @tab_query[$len..scalar(@tab_query)];

	push(@seq_align,\@cut_query);
	$dejavu{$seq}++;
	
	
	push(@alignments,$aa);
	#last;
}

 my $output_data;  
 my @tref_seq = split("",$seq_ref);  
  $output_data->{reference_sequence} = \@tref_seq;
 
 
   
  
  my $all = [];
  my $nb =0;
  foreach my $align (@seq_align){ 
     
 	
    push(@{$all->[$nb]}, @$align[0..($left*2)]);
    $nb++;
  }
  
  my @final = sort {$a->[$left] cmp $b->[$left]} @$all;
    $output_data->{traces_sequence} = \@final; 
	$output_data->{base_position} = $left;   
	#die();
	return $output_data; 


 }
 
 
sub printJson {
	my ($data)= @_;	
		print $cgi->header('text/json-comment-filtered');
		print encode_json $data;
	exit(0)
}

sub getProjectNameFromId {
	my $projId = shift;
	my $query = $buffer->getQuery();
	my $patNmae = $query->getProjectNameFromId($projId);
	return $patNmae;
}

