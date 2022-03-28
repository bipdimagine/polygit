#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use GBuffer;
use GenBoProject;
use Getopt::Long;
use Data::Dumper;
#use illumina_util;
use IO::Prompt;
use Sys::Hostname;
 use File::Find::Rule ;
use Text::Table;
use Term::Twiddle;
my $project_name;
my $project_name_origin;
my $filename;

my $patients_name;
my $user_file;
my $steps_name;
my $force;
my $type;
my $fastq_ext;
GetOptions(
	'project=s' => \$project_name,
);


my $buffer = GBuffer->new();
my $dbh = $buffer->dbh;
$dbh->{AutoCommit} = 0;
my $project = $buffer->newProject( -name => $project_name );

my $fam = $project->getPatients()->[0]->getFamily();

my $child = $fam->getChild();
my $name = $child->name;
my $url = qq{https://defidiag.polyweb.fr/cgi-bin/polymorphism-cgi//validation_variation/variations_editor.pl?project=$project_name&panel=&dv_ho=2&ac=2&edit_mode=1&never=1&this=6&impact=2&frequence=2&allele_quality=5&report_mode=1&project_summary=1&all_coverage=1&table_variations=1&sanger_variations=1&validated_variations=1&all_variations=1&denovo=1&recessive=1&xor=1&both=1&span=20&limit=30&intronic=0&utr=0&transcripts=all&user_name=pnitschk&patients=$name&annot=splicing+essential_splicing+coding+stop+phase+maturemirna+frameshift+non-frameshift&ach=2&dv=2};
print $url."\n";