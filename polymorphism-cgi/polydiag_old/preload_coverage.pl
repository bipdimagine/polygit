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
use lib "$Bin/../packages/validation_variation"; 
use html; 

#use Set::;
use Storable qw/store thaw retrieve/;
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
use infos_coverage_exons;
use preload_coverage;

use JSON::XS;
$|=1;
my $buffer = GBuffer->new();


                
my $cgi          = new CGI();
print $cgi->header('text/json-comment-filtered');   
print qq{ \{"status":"};
my $out;
my $prg =  $cgi->url(-relative=>1);
my $url_image = url(-absolute=>1);
$url_image =~ s/$prg/image_coverage.pl/;

my $project_name = $cgi->param('project');

my $utr =$cgi->param('utr')+0;
my $intronic = $cgi->param('intronic')+0;
my $limit = $cgi->param('limit');
my $padding = $cgi->param('span');

my $project = $buffer->newProject(-name=>$project_name);
my $patient_name = $cgi->param('patient_cached');
my $patients = $project->get_only_list_patients($patient_name,",");
$out.= html::print_cadre($cgi,"Exons Coverage ");
my $gene_id =  $cgi->param('gene');
my @transcripts_cgi;
if ($gene_id) {
	my @real_gene_id;
	foreach my $gid (split(",",$gene_id)){
			my $gene = $project->newGene($gid);
			
			push(@transcripts_cgi, map{$_->id} @{$gene->getTranscripts});
	}
	
	
}
else {
my $cgi_transcript =  $cgi->param('transcripts');
if ($cgi_transcript eq "all"){
	@transcripts_cgi = @{$project->bundle_transcripts() } ;
}
else {
	@transcripts_cgi = split(",",$cgi_transcript);
}
}

my @transcripts = sort{$a->getGene->external_name cmp $b->getGene->external_name} map{$project->newTranscript($_)} @transcripts_cgi ;

preload_coverage::load_coverage_transcripts($project,$patients,\@transcripts,1);
preload_coverage::load_coverage($project,$patients,\@transcripts,$padding,$utr,1);
print "@";
preload_coverage::load_coverage_primers($project,$patients,\@transcripts,1) if $project->isDiagnostic();
preload_coverage::load_coverage_regions_dup($project,$patients,\@transcripts,1);# if $project->isDiagnostic();
preload_coverage::load_coverage_list_primers($project,$patients,\@transcripts,1) unless $project->isDiagnostic();
print qq{OK"\}} ;
exit(0);
my $no =  $project->noSqlCoverage();
foreach my $patient (@$patients){
	foreach my $transcript (@transcripts){
		
		unless ($no->exists($patient->name,$transcript->id)){
			my $gc = $transcript->return_raw_coverage_obj($patient);
			$gc->patient(undef);
			$gc->chromosome(undef);
			$no->put($patient->name,$transcript->id,$gc);
		}
		unless ($no->exists($patient->name,$transcript->getGene->id)){
		my $gc = $transcript->getGene->return_raw_coverage_obj($patient);
			$gc->patient(undef);
			$gc->chromosome(undef);
			$no->put($patient->name,$transcript->getGene->id,$gc);
	}
	}
}
exit(0);
foreach my $patient (@$patients){
	print $patient->name();
	 print".";
	my @ids; 
	foreach my $transcript (@transcripts){
	 push(@ids,map{$_->id."_".$padding."_".$utr} @{ $transcript->getAllGenomicsParts()});
	}
	#eval{
	 my $hash= $no->get_bulk($patient->name,\@ids);
	  unless (scalar(keys %$hash) eq scalar(@ids)){
	my $hash2;
	foreach my $transcript (@transcripts){
	 #my $gc = $transcript->getGene->get_coverage($patient);#->coverage($transcript->start,$transcript->end);
	foreach my $exon ( @{ $transcript->getAllGenomicsParts()}) {
		my $id = $exon->id."_".$padding."_".$utr;
		#$patient->mean($exon);
		#$patient->minimum($exon);
		#next if exists $hash->{$id};
		# $exon->getTranscript()->getGene->get_coverage($patient)->coverage($exon->start,$exon->end);
		 my ($mean,$intspan,$min)  =  $exon->cached_statistic_coverage_coding(patient=>$patient,padding=>$padding,limit=>1,utr=>$utr);
		}
	
	}
	foreach my $transcript (@transcripts){
	 delete $transcript->getGene->{coverage_obj}->{$patient->id} if exists $transcript->getGene->{coverage_obj}->{$patient->id};
	}
	}
#	print qq{OK"\}};
#exit(0);
	my $primers = $project->getPrimers();
	my @pids = map {$_->id_coverage} @$primers;
	 my $hash= $no->get_bulk($patient->name,\@pids);
	 print".";
	   unless (scalar(keys %$hash) eq scalar(@pids)){
	foreach my $primer (@$primers){
		next if exists $hash->{$primer->id}->{mean};
		#$patient->mean($primer);
		 #$hash->{min} = $patient->minimum($primer);
		$primer->cached_statistic_coverage($patient);
		
	} 
}
	}
	
print qq{OK"\}};
exit(0);
	response({status=>"ok"});
sub response {
	my ($rep) = @_;
	print $cgi->header('text/json-comment-filtered');
	print encode_json $rep;
	exit(0);
}