#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use GBuffer;
#use lib "/bip-d/perl/ensembl64/ensembl/modules";
use Bio::SearchIO;
use strict;
use Getopt::Long;
use Carp;
use JSON::XS;
use Getopt::Long; 
use List::MoreUtils qw(uniq);
#use ensembl_buffer;
use Storable qw/freeze thaw nfreeze nstore_fd nstore/;
use Set::IntSpan::Fast::XS ;
use Sys::Hostname;
use Storable qw/retrieve store/;
use  Set::IntSpan::Island;
use Parallel::ForkManager;
use Array::IntSpan;
use Data::Dumper;
use Array::IntSpan;
use GenBoNoSqlLmdb;
use GenBoNoSqlAnnotation;
use GenBoNoSqlDejaVu;
#use Bio::Ensembl::Registry;
use Storable qw(dclone);
use List::MoreUtils qw(uniq natatime);
use Try::Tiny;

#ENSG00000234585
my $buffer = new GBuffer;
my $project_name= "NGS2021_3830";
my $gencode;
GetOptions(
	'project=s' => \$project_name,
	'gencode=s' => \$gencode,
	#'low_calling=s' => \$low_calling,
);


my $project = $buffer->newProject( -name 			=> $project_name );
$project->gencode_version($gencode) if $gencode;
warn "GENCODE :".$gencode;
#my $tr = $project->newTranscript("ENST00000377745");

my $trs =  $project->getCapture->transcripts_name;
#my @tt = `cat list_transcriupts.txt`;
#my $transcript = $project->newTranscript("ENST00000621098");
#warn $transcript->getProtein->name();
#my $gene = $project->newGene("NR2E3");
#warn $transcript->name();

#warn $tr->name;
#die();

foreach my $tr_name (@$trs){
	chomp($tr_name);
	try {
	my $transcript = $project->newTranscript($tr_name);
	print $transcript->name ."  OK\n";
	}
	catch {
		print $tr_name."    Problem \n";
	};
	
}
exit(0);
