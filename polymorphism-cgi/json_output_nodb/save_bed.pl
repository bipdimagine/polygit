#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;

use strict;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb/packages/";
use lib "$Bin/../packages/export/";
use lib "$Bin/../packages/layout/";
use lib "$Bin/../packages/validation_variation"; 
use lib "$Bin/../cache_nodb/scripts/";
use GBuffer;
use GBufferTest;
use export_data;
use JSON;
use VcfMerge;
use GenBoNoSql;
use Set::IntervalTree;
use Spreadsheet::WriteExcel;
use Compress::Snappy;
use Storable qw(store retrieve freeze dclone thaw);
use POSIX qw(strftime);
use CGI::Session;
use html; 
use Carp;
use Cache_Commons;
use QueryVectorFilter;
use xls_export;
use session_export;


my $io = IO::Handle->new();
$io->autoflush(1);

my $cgi = new CGI();
my $project_name = $cgi->param('project_name');
my $bed_file_name = $cgi->param('bed_file_name');
my $content = $cgi->param('content');
my $is_last_file = $cgi->param('is_last');

my @lLines = split(",",$content);

my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $project_name );
my $path = $project->getBedPolyQueryDir();
my $file = $path.'/'.$bed_file_name;

open(FILE, ">$file");
foreach my $line (@lLines) {
	$line =~ s/_/\t/g;
	print FILE $line."\n";
}
close(FILE);

my $hash;
if ($is_last_file and $is_last_file == 1) {
	my $global_file = $bed_file_name;
	$global_file =~ s/\.part[0-9]+//;
	my $type_file = $global_file.'.part*';
	
	my $cmd1 = "rm $path/$global_file";
	$hash->{cmd1} = $cmd1;
	warn $cmd1;
	system($cmd1) if -e $path.'/'.$global_file;
	my $cmd2 = "cat $path/$type_file >$path/$global_file";
	$hash->{cmd2} = $cmd2;
	warn $cmd2;
	system($cmd2);
	my $cmd3 = "rm $path/$type_file";
	$hash->{cmd3} = $cmd3;
	warn $cmd3;
	system($cmd3);
}

$hash->{'done'} = 'OK';
$hash->{'file'} = $bed_file_name;
print $cgi->header('text/json-comment-filtered');
my $json_encode = encode_json $hash;
print $json_encode;
exit(0);
	
	
