#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/kyoto";

#use lib "/bip-d/soft/distrib/tabix/latest/perl";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/coverage";
use lib "$Bin/../packages/validation_variation"; 
use lib "$Bin/../packages/cache"; 
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/utility.pm";

use draw_cnv; 
use html; 
use infos_coverage_exons;
use JSON::XS;
use image_coverage;
#use Set::;
use Carp;
use export_data;
use strict;
use Set::IntSpan::Fast::XS;
use Data::Dumper;
use GBuffer;
use Getopt::Long;
use Carp;
use Set::Intersection;
use Tabix;
use Storable qw/thaw/;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use coverage;
use Spreadsheet::WriteExcel;
use POSIX;
use validationQuery;
use Date::Tiny;
use JSON::XS;
use List::MoreUtils qw{part};

my $buffer = GBuffer->new();

my $cgi          = new CGI();
my $project_name = $cgi->param('project');
my $transcript = $cgi->param('transcript');
my $vid = $cgi->param('variation_id');

my $force_cache =  $cgi->param('force_cache');
my $project = $buffer->newProject(-name=>$project_name);
 my $variation = $project->_newVariant($vid);
my $patients = $project->getPatients();
my $CSS;
my $out_global;
$| =1;
html::print_header_polydiag($cgi,"$vid");

print $out_global;

print $cgi->start_div({class=>"panel panel-warning",style=>"font-size: 11px;font-family:  Verdana;"});
	
print $cgi->div({class=>"panel-heading"},'<i class="fa fa-flag fa-lg" style="text-align: center;font-color: white;">'."</i> Current Project : $project_name $vid");
print "<br>";
print $cgi->start_table({class=>"table table-striped  table-condensed table-bordered ",style=>"font-size: 8px;font-family:  Verdana;"});
print $cgi->start_Tr({class=>"success"});
		print $cgi->th("Patient");
		print $cgi->th( "ngs");
		print $cgi->th("ratio");
		print $cgi->th("raw depth");
	#	print $cgi->th("A");
	#	print $cgi->th("T");
	#	print $cgi->th("C");
#		print $cgi->th("G");
	#	print $cgi->th("Del");
print $cgi->end_Tr();


my $vquery = validationQuery->new(dbh=>$buffer->dbh,capture_name=>$project->validation_db());


my @header = ("ngs","ratio");

my $samtools = $buffer->software('samtools');
my $region = $variation->getChromosome->ucsc_name.":".$variation->start."-".$variation->end;	
my @letters = ("A","T","C","G","\\*");
my $ps = $project->in_this_run_patients2();
delete $ps->{total};
delete $ps->{nb_patients};
unless (exists $ps->{$project_name}){
my $i=0;
foreach my $patient (@{$project->getPatients}){
	$ps->{$project_name}->{$i} = $patient->name;
	$i++;
}
}
warn Dumper $ps;
foreach my $pr2 (keys %$ps){
	
	my @listp = map{$ps->{$pr2}->{$_}} keys %{$ps->{$pr2}};
	my $buffer2 = GBuffer->new();
	my $project2 = $buffer->newProject(-name=>$pr2);
	my $patients2 = $project2->getPatients();
	foreach my $p (@$patients2){
		next unless exists $ps->{$pr2}->{$p->id};
		my $d = $p->depth($variation->getChromosome->name,$variation->start,$variation->start)->[0];
		#my $bam = $p->getBamFile;
		#my ($res) = `$samtools depth -d 50000   $bam -r $region  | cut -f 3`;
		#chomp($res);
		#my $d2 =$res;;
		my $hvariation = utility::return_hash_variant($project,$vid,$transcript,$p,$vquery);
		if ($hvariation){
		print $cgi->start_Tr({class=>"infos",style=>"vertical-align:middle;background-color:yellowgreen"});
		}
		else {
		foreach my $i (@header){
			$hvariation->{$i} = "-";
		}
		print $cgi->start_Tr({style=>"vertical-align:middle"});
	}
	print $cgi->td($p->name);
	foreach my $i (@header){
		print $cgi->td($hvariation->{$i}) ;
	}
	  
	print $cgi->td($d);
	foreach my $l (@letters){
		
	#	my @c = $pileup =~ /$l/g;
		#print $cgi->td(scalar(@c));
	}
	
	
	print $cgi->end_Tr();
	}
}
print $cgi->end_table();
print $cgi->end_div();
exit(0);
my $similar = $project->similarProjects();
my $hres = $project->getDejaVuInfos($vid);

foreach my $p (keys %$hres){
	delete $hres->{$p} unless exists $similar->{$p};
}
foreach my $p (keys %$hres){

	project_table($p);
	
}


exit(0);


sub project_table {
	my ($pname) = @_;
	my $buffer = GBuffer->new();
	my $project = $buffer->newProject(-name=>$pname);
	
print $cgi->start_div({class=>"panel panel-info",style=>"font-size: 11px;font-family:  Verdana;"});
print $cgi->div({class=>"panel-heading"},'<i class="fa fa-flag fa-lg" style="text-align: center;font-color: white;">'."</i> $pname ");
print "<br>";
print $cgi->start_table({class=>"table table-striped  table-condensed table-bordered ",style=>"font-size: 8px;font-family:  Verdana;"});
print $cgi->start_Tr({class=>"success"});
		print $cgi->th("Patient");
		print $cgi->th( "ngs");
		print $cgi->th("ratio");
print $cgi->end_Tr();
	
	my $patients = $project->getPatients();
	foreach my $p (@$patients){
	my $hvariation = utility::return_hash_variant($project,$vid,$transcript,$p,$vquery);
	if($hvariation){
	print $cgi->start_Tr({class=>"infos"},{style=>"vertical-align:middle"});
	}
	else {
		next;
		foreach my $i (@header){
			$hvariation->{$i} = "-";
		}
		print $cgi->start_Tr({style=>"vertical-align:middle"});
	}
	print $cgi->td($p->name);
	foreach my $i (@header){
		
	print $cgi->td($hvariation->{$i}) ;
	}
	
	
	print $cgi->end_Tr();
}
print $cgi->end_table();
print $cgi->end_div();
	
}
#
#my $similar = $project->similarProjects();
#	my $similar = $project->similarProjects();
#	my $hres = $project->getDejaVuInfosForDiag($vid);
#	
#my $hdj = $project->getDejaVuInfos($vid);
