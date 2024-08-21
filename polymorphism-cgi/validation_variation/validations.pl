#!/usr/bin/perl
use CGI;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo";
use lib "$Bin/../../GenBo/lib/GenBoDB";
use lib "$Bin/../../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/validations";
#use lib "$Bin/../GenBo/lib/obj";
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use Carp;
use JSON::XS;
use GenBoFilter;
use validationQuery;

my $cgi    = new CGI;

#chargement du buffer 




my $buffer = GBuffer->new;
$buffer->dbh->{AutoCommit} = 0;
my $project_name = $cgi->param('project');
my $project = $buffer->newProject(-name=>$project_name);
die() unless $project;
my $capture = $project->getCapture();

my $vquery = validationQuery->new(dbh=>$buffer->dbh,capture_name=>$capture->validation_db());
my $samples = $cgi->param('samples');
my $patients = $project->get_list_patients($samples);

#my $query1 = qq{ use Polyexome_HG19; };
#$buffer->dbh->do($query1);
#$buffer->{version_polyproject} = "Polyexome_HG19";

my %validations = ( 
					valid=>1,
					unvalid=>-1,
					notsequenced =>0
);
					
my $vid  = $cgi->param('variation_id');

#

#warn Dumper $last_type;
my $contig;

my $poly_id  =  $cgi->param('poly_id');
my $variation;
($variation->{chr_name},$variation->{pos},$variation->{ref},$variation->{all}) = split("_",$poly_id);
my $vcf_id;
my $chr = $project->getChromosome($variation->{chr_name});
$variation->{pos_vcf} = $variation->{pos};
$variation->{type} = "variation";
my $seq_ref = $chr->getSequence($variation->{pos_vcf},$variation->{pos_vcf});

if (length($variation->{ref}) < length ($variation->{all}) ){
	$variation->{ref} = $seq_ref;
	$variation->{type} = "deletion";
	#$variation->{pos_vcf} -= length($variation->{ref} );
}
elsif (length($variation->{ref}) > length ($variation->{all}) ) {
	$variation->{ref} = $seq_ref;
	$variation->{all} = $seq_ref.$variation->{all}; 
	$variation->{type} = "insertion";
	$variation->{pos_vcf} -= length($variation->{ref} );
}
else {
	$variation->{type} = "variation";
}



my $vcf_id = join("_",$chr->ucsc_name,$variation->{pos_vcf},$variation->{ref},$variation->{all});


my $patients = $project->get_list_patients($samples);
 
my $value = $cgi->param("value");
my $user = $cgi->param('user');
my $tabix =  $buffer->{config}->{software}->{tabix};
my $samtools = $buffer->{config}->{software}->{samtools};
my %dejavu;
foreach my $patient (@$patients){
	
	my  ($v) = grep {$_->id eq $poly_id} @{$patient->getStructuralVariations()};
	$dejavu{$v->id} ++;
	die() if scalar keys(%dejavu) >1;
	die("no variant with this id") unless $v;
	my $vid = $vquery->createVariation(polyid=>$v->id,vcfid=>$v->vcf_id,genboid=>$vid,
									version=>$project->getVersion);
									
	my $vcf_files;
	if ($v->isVariation()){
		 $vcf_files = $patient->getVariationsFiles();
	}
	my ($vcf_chr,$vcf_pos,$vcf_ref,$vcf_all) = split("_",$v->vcf_id);
	my $tabix_position = $vcf_chr.":".$vcf_pos."-".$vcf_pos; 
	my $vcf_line = "";
	foreach my $vcf_file (@$vcf_files){
	die() unless -e $vcf_file;
	
	 #$vcf_pos = $v->{line_infos}->{POS};
	
	#my $seq = $vcf_all;
	my @lines = grep {$_ =~ /$vcf_all/} `$tabix $vcf_file  $tabix_position`;
	next unless scalar(@lines);
	die ("2 lines $tabix $vcf_file  $tabix_position") if scalar(@lines) ne 1;
	 $vcf_line .= $lines[-1];
	chomp($vcf_line);
	}
	my $bam_file=$patient->getBamFile();
	my $st_pos = $chr->ucsc_name.":".$v->start."-".$v->end;
	my $bam_lines2;# = `$samtools view -h $bam_file  $st_pos`;
	$bam_lines2="";
	my $res = $vquery->createValidation (variation_id=>$vid,project=>$project->name(),
									  sample=>$patient->name(),user=>$user,
									  	vcf_line=>$vcf_line,validation_status =>$value,
									  	,bam_line=>$bam_lines2,method=>"iontorrent");
sendError("Insertion Problem") unless $res == 1;
}
$buffer->dbh->commit();	
sendOK("ok");	

exit(0);


#
#my $chr = $variation->getChromosome();
#my $pos = $variation->position($chr);
#
#my $chr_sequence = $chr->sequence($pos->start,$pos->end);
#
#my $vcf_position = $pos->start-1;
#
#
#my $tabix_position = $variation->getTabixPosition();
#
#
#
#if ($variation->isIndel()){
#	$vcf_file = $patients->[0]->getIndelsFile();
#}
#my $tabix =  $buffer->{config}->{software}->{tabix};
#my $samtools = $buffer->{config}->{software}->{samtools};
#warn qq{$tabix $vcf_file  $tabix_position};
#my @lines = `$tabix $vcf_file  $tabix_position`;
#confess() if scalar(@lines) ne 1;
#my $vcf_line = $lines[0];
#chomp($vcf_line);
#
#my $bam_file=$patient->getBamFile();
#my $st_pos = $chr->ucsc_name.":".$pos->start()."-".$pos->end;
#my $bam_lines2 = `$samtools view -h $bam_file  $st_pos`;
#
# 
#
#my (@args) = split(" ",$vcf_line);
#my $polyid = join("_",($chr->ucsc_name(),$pos->start(),$chr_sequence,$variation->sequence()));
#my $vcfid = join("_",($args[0],$args[1],$args[3],$args[4]));
#my $vid = $vquery->createVariation(polyid=>$polyid,vcfid=>$vcfid,genboid=>$variation->id,
#									version=>$project->getGoldenPath);
#sendError("Insertion Problem") unless $vid;
#									
#my $res = $vquery->createValidation (variation_id=>$vid,project=>$project->name(),
#									  sample=>$patient->name(),user=>$user,
#									  	vcf_line=>$vcf_line,validation_status =>$value,
#									  	,bam_line=>$bam_lines2,method=>"iontorrent");
#
#sendError("Insertion Problem") unless $res == 1;
#sendOK("ok");		
#exit(0);



	




sendOK("ok");
sub sendOK {
	my ($text) = @_;
	my $resp;
	$resp->{status}  = "OK";
	$resp->{message} = $text;
	response($resp);
}

sub sendError {
	my ($text) = @_;
	my $resp;
	$resp->{status}  = "error";
	$resp->{message} = $text;
	response($resp);
}

sub response {
	my ($rep) = @_;
	print $cgi->header('text/json-comment-filtered');
	print encode_json $rep;
	exit(0);
}
