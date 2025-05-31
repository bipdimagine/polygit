#!/usr/bin/perl
$|=1;

use Carp;
use strict;
use JSON;
use Data::Dumper;
use CGI qw/:standard :html3/;
use Set::IntSpan::Fast::XS;
use Set::IntervalTree;
use List::Util qw[min max];
use FindBin qw($Bin);
use Storable qw(store retrieve freeze);
use Number::Format qw(:subs);

use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use Digest::MD5::File qw(dir_md5_hex file_md5_hex url_md5_hex file_md5);
use GBuffer;
use GenBoProject;
use GenBoCache;
use Time::HiRes;
use layout;
use export_excel;
use export_data;
use Capture::Tiny ':all';


my $cgi = new CGI;
my $projectname = $cgi->param('project');
my $patientname     = $cgi->param('patient');
my $eid    = $cgi->param('id');
my $type    = $cgi->param('type');



# pour récupérer les objets project et patient
my $buffer = GBuffer->new();
my $project = $buffer->newProjectCache( -name => $projectname);
my $patient = $project->getPatient($patientname);
my $dirsv = $project->getCacheSV(). "/rocks/";
my $rockssv = GenBoNoSqlRocks->new(dir=>"$dirsv",mode=>"r",name=>"sv");
my $dir = $project->getCacheCNV(). "/rocks/";
my $rocks = GenBoNoSqlRocks->new(dir=>"$dir",mode=>"r",name=>"cnv");

	

my $dir = $project->getCacheCNV(). "/rocks/";
my $cnv = $rocks->get($eid);
$cnv = $rockssv->get($eid) unless $cnv;
print $cgi->header('text/json-comment-filtered');
my $h;
$h->{project} = $projectname;
my @g = sort{$b->{score} <=> $a->{score}} @{$cnv->{genes}};
$h->{genes} =\@g;;
print encode_json $h;
print "\n";
