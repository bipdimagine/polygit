#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";

use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use colored;
use File::Basename;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;


my $buffer = GBuffer->new();
my $version = 33.1;

my @versions = keys %{$buffer->gencode};

foreach my $version (sort {$a <=> $b}@versions){
print "test version $version : ";
my $dir_gencod =  $buffer->public_data_root."/HG19/". $buffer->gencode->{$version}->{directory};
my $no = GenBoNoSqlAnnotation->new(dir=> $dir_gencod,mode=>"r");

my $capture_name = "BONEome";
my $capture_id =  $buffer->getQuery()->getCaptureId($capture_name);

my $hquery =  $buffer->getQuery()->getCaptureTranscripts($capture_id);


#$self->{directory}->{$version}->{$database} = $self->public_data_root."/".$self->annotation_genome_version."/".$self->buffer->gencode->{$version}->{directory};
#
my $error ;
foreach my $tr_id (@{$hquery->{transcripts_name}}) {
					unless ($no->get("synonyms", $tr_id)) {
						$error->{$tr_id} ++;
						#warn Dumper $hquery->{transcripts}->{$tr_id};
					}
}
unless (keys %$error){
	print "OK \n";
	next;
}

print "ERROR  \n";
foreach my $tr_id (keys %$error) {
	print "\t\t $tr_id ".join(" - ",@{$hquery->{transcripts}->{$tr_id}})."\n";
}

}