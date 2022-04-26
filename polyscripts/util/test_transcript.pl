#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/obj-nodb/";
#use lib "//software/polyweb/poly-disk/poly-src/GenBo/lib/obj-nodb/";
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
my $project_name= "NGS2021_4598";
my $gencode;
GetOptions(
	'project=s' => \$project_name,
	'gencode=s' => \$gencode,
	#'low_calling=s' => \$low_calling,
);
      

my $project = $buffer->newProjectCache( -name 			=> $project_name );


#$gencode= 34;
warn $project->gencode_version();
$project->gencode_version($gencode) if $gencode;
my $query = $buffer->getQuery();
my $capture = $project->getCapture();
my $hquery = $query->getCaptureTranscripts( $capture->id );
foreach my $t (@{$hquery->{transcripts_name}}){
	
	eval {
		warn $t;
	my $tr = $project->newTranscript("$t");
	
	};

}
die();
