#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;

use strict;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/validation_variation"; 
use GBuffer;
use export_data;
use JSON;
use VcfMerge;
use GenBoNoSql;
use Spreadsheet::WriteExcel;
use Compress::Snappy;
use Storable qw(store retrieve freeze dclone thaw);
use POSIX qw(strftime);
use CGI::Session;
use html; 


my $cgi = new CGI();
my $projectName		= $cgi->param('project');
my $type			= $cgi->param('type');
my $patients		= $cgi->param('patients');

my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $projectName );

my $hPatients;
foreach my $pat_name (split(',', $patients)) {
	$hPatients->{$pat_name} = undef;
}

my (@list, $hRes);
foreach my $capture (@{$project->getCaptures}){
	my $h;
	$h->{'fam'} = $project->name();
	$h->{'name'} = $capture->name();
	$h->{'type'} = 'bed';
	$h->{'type_desc'} = 'Capture (BED)';
	$h->{'file'} = $capture->gzFileName();
	if (-e $h->{'file'}) { $h->{'enabled'} = 'yes'; }
	else { $h->{'enabled'} = 'no' };
	$h->{'file'} =~ s/.+\/public-data\//\/NGS\/PUBLIC\//;
	#push(@list, $h); 
}

my @lType;
if ($type eq 'all' or $type eq 'bam') { push(@lType, 'bam'); }

foreach my $patient (@{$project->getPatients()}) {
	if ($patients) {
		next unless (exists $hPatients->{$patient->name()});
	}
	foreach my $this_type (@lType) {
		my $h;
		$h->{'fam'} = 'All';
		$h->{'fam'} = $patient->family() if ($patient->family());
		$h->{'name'} = $patient->name();
		if ($this_type eq 'bam') {
			$h->{'type'} = 'bam';
			$h->{'type_desc'} = 'Align (BAM)';
			$h->{'file'} = $patient->getBamFile();
			if (-e $h->{'file'}) { $h->{'enabled'} = 'yes'; }
			else { $h->{'enabled'} = 'no' };
			$h->{'file'} =~ s/.+\/ngs\//\/NGS\//;
			push(@list, $h); 
		}
	}
}
$hRes->{'label'} = 'name';
$hRes->{'items'} = \@list;
print $cgi->header('text/json-comment-filtered');
print encode_json $hRes;
exit(0);