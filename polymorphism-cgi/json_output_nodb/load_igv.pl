#!/usr/bin/perl
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use connect;
use GBuffer;
use GenBoStorable;
use Data::Dumper;
use util_file;
use Storable qw/freeze thaw nfreeze nstore_fd nstore retrieve/;
use export_data;
#use threads;
#use threads::shared;
#use Getopt::Long;
#use Set::Object; 
use Set::Intersection;



my $cgi    = new CGI();

my $buffer = new GBuffer;
my $project_name = $cgi->param('project');

#my ($project_name) ;
#
#GetOptions( 'project=s' => \$project_name,
#			);
			

my $project = $buffer->newProject( -name => $project_name );
die( "unknown project" . $project_name ) unless $project;

my $server =  'https://'.$ENV{HTTP_HOST};
#warn Dumper $ENV;

my $root_dir = $project->getRootDir();
my $ln_dir = $root_dir."/align"; 



open (OUT,">".$ln_dir."/igv.xml");

my $genome = lc($project->getVersion());


if ($genome =~ /hg19/) {$genome = 'hg19';}
elsif ($genome =~ /hg38/) {$genome = 'hg38';}
elsif ($genome =~ /mt/) {$genome = 'hg38';}
else {$genome = $project->getVersion();}




my $patients = $project->getPatients();
my $filter_patient = $cgi->param('patients');
if ($filter_patient) {
	my @pp;
	foreach my $p (@$patients){
		my $pname = $p->name();
		push(@pp,$p) if $filter_patient=~ /$pname/; 
	}
	$patients = \@pp;
}
print OUT qq{<?xml version="1.0" encoding="UTF-8"?>
<Global genome="$genome" version="1">
   <Resources>
};
#warn $server;
$server = "www.polyweb.fr" if $server eq "bipd";
$server = "www.polyweb.fr" if $server =~/10\.200\.27/;
my $base_url =  $server;
my $temp_url = $cgi->url();
warn $temp_url;
if ($temp_url =~ /^http/){
	$temp_url =~ s/http:\/\///;
	my @url = split("/",$temp_url);

	$base_url= $url[0];
}

 system ("chmod  a+w $ln_dir/igv.xml");  
my $header = qq{<Resource path="};

my $url1 = "http://$server/NGS/".$project->name()."/".$project->getVersion()."/align/";
my $url = "http://$server/NGS/".$project->name()."/".$project->getVersion()."/align/";
my $url_variations = "http://$server/NGS/".$project->name()."/".$project->getVersion()."/variations/";


 my $url_capture =  "http://$server/public-data/".$project->version."/capture/";
 
 my %beds;
 
foreach my $patient (@$patients){
	my $capture = $patient->getCapture();
 	my $file = $capture->gzFileName();
 	my $t = $url_capture. $capture->type . "/" . $capture->file_name . ".gz";
 	$beds{$t} ++;
 	
 }

my $footer= qq{"/>};

foreach my $patient (@$patients){
	 my $file = $patient->getBamFile(undef,1);

	 next unless -e $file;

	 my ($p,$a) = split ("align",$file);
	 #my $file_vcfs = $patient->getVariationsFiles();
	 	# my ($p1,$vf) = split ("variations",$file_vcf);
	
        print OUT $header.$url1.$a.$footer."\n";
        #  print OUT $header.$url_variations.$vf.$footer."\n";
}

foreach my $k (keys %beds){
	 print OUT $header.$k.$footer."\n";
}

print OUT  qq{</Resources>
</Global>};
my $data_out;
my $data;
$data->{id} = 1;

$data->{url} = $url."igv.xml";
push(@$data_out,$data);
export_data::print($project,$cgi,$data_out);

exit(0);

