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
use List::MoreUtils qw/pairwise/;

my $nb_gene_by_patient = 3;
my $nb_exon_by_genes = 10;

my $buffer = GBuffer->new();

my $cgi          = new CGI();
my $project_name = $cgi->param('project');

my $project = $buffer->newProjectCache(-name=>$project_name);
my $runs = $project->getRuns();



my $patient_name = $cgi->param('patient');
my $cgi_transcript =  $cgi->param('transcript');
my $cgi_exon = $cgi->param('exon');

$patient_name ="all" unless $patient_name;
my $patients = $project->get_list_patients($patient_name," ");


my $captures = $project->getCaptures();


my $patient_name = $cgi->param('patients');
my $cgi_transcript =  $cgi->param('transcripts');
my $type = $cgi->param('type');

my $vp =  $project->getPatient($patient_name);
my $run = $vp->getRun();

my $splice_5 = 50;
my $limit    = 0;

$limit = $cgi->param('limit');
$splice_5 = $cgi->param('span');

my $start_cgi = $cgi->param('start');
my $end_cgi = $cgi->param('end');

$splice_5 = 1 unless $splice_5;
$splice_5 = 100;
#$splice_5= 10;
my $exons_red;
my $common_red;
my $data; 
my $all_htranscripts;
my $patient_exons;
my $tr = $project->newTranscript($cgi_transcript);

my $cgi    = new CGI();
my @res;
my $seq;
my $pdata;
my $estart;
my $eend;
my $is_exon;
my $nbp;
my @names;
my $utr;
my $length_exon;
push(@names,$patient_name);
return_coverage($cgi->param('start'),$cgi->param('end'));


exit(0);

sub return_array_depth {
	my ($patient,$chr_name,$start,$end) = @_;
	if ($type eq "raw"){
		return $patient->depth($chr_name,$start,$end);
	}
	else {
			return $patient->normalize_depth($chr_name,$start,$end);
	}
}

sub return_coverage{
	my ($start,$end) =@_;
	warn $start." ".$end;
	my %data_coverage;
	my @res;
	my $nb_patient_in_run = 0;
	my $vexon;
	my $pdata;
		foreach my $exon (@{$tr->getAllGenomicsParts}){
		next if $exon->name() ne $cgi_exon;
		$vexon = $exon;
		last;
		}
		my $genomic_start =   $start-$splice_5;
		my $genomic_end =   $end+$splice_5;
		$pdata = return_array_depth($vp,$tr->getChromosome->name,$genomic_start,$genomic_end) ;
		foreach my $patient (@{$run->getPatients}){
			next if $patient->name eq $vp->name;
			$data = return_array_depth($patient,$tr->getChromosome->name,$genomic_start,$genomic_end);
			$nb_patient_in_run ++;
		 	@res = pairwise{$a +$b} @res,@$data;
		 	last if $nb_patient_in_run > 10;
		}

	

my @seq = split("",$tr->getChromosome()->getSequence($genomic_start,$genomic_end));
$data_coverage{$patient_name} = $pdata;
	
$data_coverage{all} = \@res;

	
my $nb_point = scalar(@res);	
my $out;
  $out->{exon_name} = $cgi_exon;
  $out->{transcript_name} = $tr->name;
   $out->{gene_name} = $tr->getGene()->external_name();
	my $data;
	my $pos= $genomic_start ;
	my $len =0;
	my $big;
	$big =1 if (scalar(@$pdata)) > 999;
	my $csv;
	my $starti = 0;
	my @sequences;
	for (my $i = 0;$i<scalar(@$pdata);$i++){
		$len ++;
		my $covm = int($data_coverage{all}->[$i] /$nb_patient_in_run);
		my $cov = $data_coverage{$patient_name}->[$i];
		my $value; 
		my $valuem; 
		push(@{$csv}, [$starti,$cov+0,$covm+0]);
		$starti++;
		push(@sequences,shift @seq);
	if ($big){
		 $value = $cov;
	 	$valuem = $covm;
	}
	else {
	 #	$value = { y=> $cov, name=>"Position:". $pos."[".$seq[$i]."]"};
	 	#$valuem = { y=> $covm, name=>"Position:". $pos."[".$seq[$i]."]"};
	}
	$pos++;
	
   #	push(@{$data->[0]->{data}},$value);
#	push(@{$data->[1]->{data}},$valuem);

   }
   
my $span =  Set::IntSpan::Fast::XS->new($genomic_start."-".$genomic_end);
my $exonic_span  = $tr->getSpanCoding->intersection($span);
my $intronic_span = $span->diff($tr->getGenomicSpan);
 my $span2 = $tr->getGenomicSpan->intersection($span);
 
my $utr_span = $span2->diff($tr->getSpanCoding);
my $legende;
push(@$legende, @{return_code_highChart($intronic_span,$genomic_start,"#E6E2EB","intronic")});
push(@$legende, @{return_code_highChart($utr_span,$genomic_start,"#FAFAFA","utr")});
 push(@$legende, @{return_code_highChart($exonic_span,$genomic_start,"#DFF0D8","exonic")});

	   $out->{legende} = $legende;

  push(@names,$project->name());
   $out->{sequences}=\@sequences;
  $out->{data}=$data;
  $out->{len} = $len;
   $out->{names} = \@names;
   $out->{labels} = ['x',@names];
  $out->{start} = $splice_5;
  $out->{genomic_start} = $genomic_start;
  $out->{end} =  - 1;
  $out->{csv} = $csv;
  my $exonic = return_start_end_intspan($exonic_span,$genomic_start);
  $out->{exonic} = $exonic;
  $out->{utr} =return_start_end_intspan($utr_span,$genomic_start);
   $out->{intronic} =return_start_end_intspan($intronic_span,$genomic_start);
   $out->{exon_start} = $exonic->[0]->{from};
   $out->{exon_end} = $exonic->[0]->{to};
  my $strand = "forward";
   $strand = "reverse" if $tr->strand == -1;
   my $tpos = $tr->getChromosome->name.":".$tr->start."-".$tr->end;
   my $epos = "[".$genomic_start."-".$genomic_end."]";
     $epos = "[".$vexon->start."-".$vexon->end."]" if $vexon;
     my $vname = "";
     $vname = $vexon->name() if $vexon;
  $out->{title} = "normalized average depth :".$tr->getGene->external_name." ".$tr->name." ".$tr->external_name." $tpos $strand"." <br> ".$vname." ".$epos ;#if $vexon;
export_data::print($project,$cgi,$out);

exit(0);
	
	
	
}

sub return_start_end_intspan{
	my ($span,$genomic_start) = @_;
	my $tab= [];
	my $iter = $span->iterate_runs();
    while (my ( $from, $to ) = $iter->()) {
    		my $pos;
    		push(@$tab,{from=>abs($from-$genomic_start),to=>abs($to-$genomic_start)});
    }
    return $tab;
}
sub return_code_highChart{
	my ($span,$genomic_start,$color,$text) = @_;
	my $tab =[];
	
	my $iter = $span->iterate_runs();
    while (my ( $from, $to ) = $iter->()) {
	my $res;
	$res->{color}->{linearGradient}->{x1} = 0;
	$res->{color}->{linearGradient}->{x2} = 0;
	$res->{color}->{linearGradient}->{y1} = 0;
	$res->{color}->{linearGradient}->{y2} = 5;
	$res->{color}->{stops} = [
                [0, $color],
                [1, $color]
                ] ;
	
		$res->{from} = $from - $genomic_start-0.5;
		$res->{to} = $res->{from}+ ($to - $from)+1;
		$res->{label}->{text} = "$text";
		push(@$tab,$res);
    }
    return $tab;
}
