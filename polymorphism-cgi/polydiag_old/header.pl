#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/kyoto";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
#use lib "/bip-d/soft/distrib/tabix/latest/perl";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/coverage";
use lib "$Bin/../packages/validation_variation"; 
use html; 
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

use image_coverage;




my $buffer = GBuffer->new();

my $cgi          = new CGI();
html::print_header_polydiag($cgi);
my $project_name = $cgi->param('project');

my $project = $buffer->newProject(-name=>$project_name);
my $captures = $project->getCaptures();
my $similar = $project->similarProjects();

my $vquery = validationQuery->new(dbh=>$buffer->dbh,capture_name=>$project->validation_db());
my $patient_name = $cgi->param('patients');
#my $edit_mode = $cgi->param('edit_mode');
$patient_name ="all" ;#unless $patient_name;
my $patients = $project->get_list_patients($patient_name,",");;
my $cgi_transcript =  $cgi->param('transcripts');
my $user =  $cgi->param('user');
$cgi_transcript = "all";

my @transcripts_cgi ;
if ($cgi_transcript eq "all"){
	my %tr;
	foreach my $capture (@$captures){
		next unless $capture->transcripts_name();
		map{$tr{$_}++} @{$capture->transcripts_name()} ;
	}
	@transcripts_cgi = keys %tr;
}
else {
	@transcripts_cgi = split(",",$cgi_transcript);
}
my $means;
my $total;
my $coverage;



my ($date,$time) = split(" ",$project->creation_date);
my ($year,$month,$day) =  return_date($project->creation_date);


	

my $h1 = "info";
my $h2 = "primary";



my $div1 = qq{<div class="btn-group" role="group">};
my $style_b ="font-size: 10px;letter-spacing:2px;font-family: Gill Sans, Verdana;text-transform:uppercase; " ;
my $style_a ="font-size: 8px;font-family:Verdana;letter-spacing:1px;	padding-left: 10px;text-transform:uppercase; ";
print $cgi->start_div({class=>"btn-group  btn-toolbar-xs btn-group-justified",role=>"group"});

	print $cgi->start_div({class=>"btn-group btn-group-xs", role=>"group"});
		print start_button({class=>"btn btn-$h1",style=>$style_b});
			print "PolyDiag &nbsp;";
				print $cgi->span({class=>"badge ",style=>$style_a},"<i class='fa fa-user''></i>$user");
		print end_button();
	print  $cgi->end_div;
	
	print $cgi->start_div({class=>"btn-group btn-group-sm", role=>"group"});
		print start_button({class=>"btn btn-$h2",style=>$style_b});
			print $project->name()."&nbsp;&nbsp;";
			print $cgi->span({class=>"badge ",style=>$style_a},"$day-$month-$year");
		print end_button();
	print  $cgi->end_div;
	
		print $cgi->start_div({class=>"btn-group btn-group-xs", role=>"group"});
		print start_button({class=>"btn  btn-$h1 btn-block",style=>$style_b});
			print "amplification &nbsp;";
			print $cgi->span({class=>"badge ",style=>$style_a},join("<br>",map{$_->name}@$captures));
		print end_button();
	print  $cgi->end_div;
		
		print $cgi->start_div({class=>"btn-group btn-group-sm", role=>"group"});
		print start_button({class=>"btn  btn-$h2",style=>$style_b});
			print "TRANSCRIPTS &nbsp;";
			print $cgi->span({class=>"badge ",style=>$style_a, id=>"span_transcripts"},scalar(@transcripts_cgi));
		print end_button();
	print  $cgi->end_div;
	
			print $cgi->start_div({class=>"btn-group btn-group-xs", role=>"group"});
		print start_button({class=>"btn  btn-$h1",style=>$style_b});
			print "run&nbsp;&nbsp;";
			print   $cgi->span({class=>"badge ",style=>$style_a},scalar(@{$project->getRuns}));
		print end_button();
	print  $cgi->end_div;
	
			print $cgi->start_div({class=>"btn-group btn-group-sm", role=>"group"});
		print start_button({class=>"btn  btn-$h2",id=>"bubble_sample",style=>$style_b});
			print "samples &nbsp";
			print  $cgi->span({class=>"badge",style=>$style_a,jsId=>"span_patients",id=>"span_patients"},scalar(@{$project->getPatients}));
		print end_button();
	print  $cgi->end_div;
	


	
print $cgi->end_div;





print qq{
	 <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.2/jquery.min.js"></script>
    <!-- Include all compiled plugins (below), or include individual files as needed -->

  
};

exit(0);

sub start_button {
	my ($arg) = @_;
	my $class = $arg->{class};
	my $id = $arg->{id};
	my $style = $arg->{style};

	return qq{<button class="$class" style="$style"  id ="$id">};
}
sub end_button{
	return "</button>";
}
	sub return_date {
		my ($dd) = @_;
		 my @amonths = ('Jan', 'Feb', 'Mar', 'Apr','May',"Jun","Jul","Aug","Sep","Oct","Nov","Dec");
		my ($date,$time) = split(" ",$dd);
	    my ($year,$month,$day) =  split("-",$date);
		return ($year,$amonths[$month-1],$day);
	}