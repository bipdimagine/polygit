#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;

use strict;
use Data::Dumper;
use FindBin qw($RealBin);
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use GBuffer;
use JSON;

my $cgi = new CGI();
my $project_name = $cgi->param('project');
my $chr_name = $cgi->param('chr');

my $buffer = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProject( -name => $project_name );
my $dir_log = $project->getCacheBitVectorDir()."/log/";

my @lPb;
foreach my $step_name ('store_ids', 'store_annotations', 'loh') {
	next if ($step_name eq 'loh' and not $project->isSomaticStudy());
	my $h;
	my $cmd = "$RealBin/cache_check_step.pl -project=$project_name -chr=$chr_name -step=$step_name -no_verbose=1";
	eval {
		`$cmd`;
	};
	if (not is_exists_ok_file($dir_log, $step_name, $chr_name)) { push(@lPb, $step_name); }
}

my (@lItems, $hStatus);
if (scalar(@lPb) == 0) {
	$hStatus->{status} = "<font color='green'><b>Check OK !</b></font>";
}
else { $hStatus->{status} = "<font color='red'><b>ERROR: </b>".join(', ', @lPb)."</font>"; }
push(@lItems, $hStatus);

my $hJson;
$hJson->{'label'} = 'status';
$hJson->{'items'} = \@lItems;
print $cgi->header('text/json-comment-filtered');
print encode_json $hJson;

sub is_exists_ok_file {
	my ($dir_log, $step_name, $chr_name) = @_;
	return 1 if (-e "$dir_log/check_$step_name.$chr_name.ok");
	return;
}