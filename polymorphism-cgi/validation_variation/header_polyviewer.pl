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

use Capture::Tiny ':all';

#exit(0);
$| =1;
my $version = "2";
my $cgi          = new CGI();
my $buffer = GBuffer->new();
my $project_name = $cgi->param('project');
my $panel_name =  $cgi->param('panel');
my $project = $buffer->newProject(-name=>$project_name);
my $user =  $cgi->param('user');
my $no_cache = $project->get_lmdb_cache_summary("r");
my $cache_id = "$project_name;$user;$version".".title";
my $header = $no_cache->get_cache($cache_id);

if ($header){
	print $header;
	exit(0);
}
my $stdout = tee_stdout {
html::print_header_polydiag($cgi);
print qq {
<style type="text/css"> 
.colorgraph1 {
  height: 3px;
  border-top: 0;
  background: #c4e17f;
  border-radius: 5px;
  margin:0px;
 # background-image: -webkit-linear-gradient(left, #c4e17f, #c4e17f 12.5%, #f7fdca 12.5%, #f7fdca 25%, #fecf71 25%, #fecf71 37.5%, #f0776c 37.5%, #f0776c 50%, #db9dbe 50%, #db9dbe 62.5%, #c49cde 62.5%, #c49cde 75%, #669ae1 75%, #669ae1 87.5%, #62c2e4 87.5%, #62c2e4);
  #background-image: -moz-linear-gradient(left, #c4e17f, #c4e17f 12.5%, #f7fdca 12.5%, #f7fdca 25%, #fecf71 25%, #fecf71 37.5%, #f0776c 37.5%, #f0776c 50%, #db9dbe 50%, #db9dbe 62.5%, #c49cde 62.5%, #c49cde 75%, #669ae1 75%, #669ae1 87.5%, #62c2e4 87.5%, #62c2e4);
  #background-image: -o-linear-gradient(left, #c4e17f, #c4e17f 12.5%, #f7fdca 12.5%, #f7fdca 25%, #fecf71 25%, #fecf71 37.5%, #f0776c 37.5%, #f0776c 50%, #db9dbe 50%, #db9dbe 62.5%, #c49cde 62.5%, #c49cde 75%, #669ae1 75%, #669ae1 87.5%, #62c2e4 87.5%, #62c2e4);
#  background-image: linear-gradient(to right,  rgb(54, 57, 69) 14.5%, rgb(54, 57, 69) 15%, rgb(16, 121, 178) 14%, rgb(16, 121, 178) 57%, rgb(150, 45, 62) 57%, rgb(150, 45, 62));
   background-image:linear-gradient(to right, rgb(54, 57, 69) , #F0F0F0 50%,rgb(54, 57, 69))
}
.container1 {
  max-width: 50%;  
  margin: 40px auto;
}

.hr-text {
  line-height: 1em;
  position: relative;
  outline: 0;
  border: 0;
  color: black;
  text-align: center;
  height: 1.5em;
  opacity: .5;
  &:before {
    content: '';
    // use the linear-gradient for the fading effect
    // use a solid background color for a solid bar
    background: linear-gradient(to right, transparent, #818078, transparent);
    position: absolute;
    left: 0;
    top: 50%;
    width: 100%;
    height: 1px;
  }
  &:after {
    content: attr(data-content);
    position: relative;
    display: inline-block;
    color: black;

    padding: 0 .5em;
    line-height: 1.5em;
    // this is really the only tricky part, you need to specify the background color of the container element...
    color: #818078;
    background-color: #fcfcfa;
  }
}

</style>
};
};


my $captures = $project->getCaptures();
my $similar = $project->similarProjects();
my $ver1 = $project->genome_version;

 my $ver = $project->annotation_version;
my $machine = $project->getRuns->[0]->machine();

my $vquery = validationQuery->new(dbh=>$buffer->dbh,capture_name=>$project->validation_db());
my $patient_name = $cgi->param('patients');
#my $edit_mode = $cgi->param('edit_mode');
$patient_name ="all" ;#unless $patient_name;
my $patients = $project->get_list_patients($patient_name,",");
my $cgi_transcript =  $cgi->param('transcripts');

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
# my $cno = $project->getChromosomes()->[0]->lmdb_hash_variants("r");
my $date_cache =  "";#utility::return_date_from_file($cno->filename);

	

my $h1 = "info";
my $h2 = "primary";
my @color =("#2d3e50","#FECF71","#2d3e50","#F0776C","#2d3e50","#BFD641") ;


my $fcolor = "white";
my @tt = split("",$project_name);
my $r = $tt[-1];
#$r =  int(rand(7));
my $ccolor = "rgb(150, 45, 62)";
my $bcolor = "rgb(16, 121, 178)";

if ($r ==1){
	$ccolor ="#A20025" ; $bcolor= "#F0A30A";$fcolor="white";
}
elsif ($r ==2){
	$ccolor ="#50394C"; $bcolor= "#F4E1D2";$fcolor="black";
}
elsif ($r ==3){
	$ccolor ="#618685"; $bcolor= "#36486B";$fcolor="white";
}
elsif ($r ==4){
$ccolor ="#b8a9c9"; $bcolor= "#622569";$fcolor="white";
}
elsif ($r ==5){
	$ccolor ="#3e4444"; $bcolor= "#82b74b";$fcolor="white";
}
elsif ($r ==6){
	$ccolor =" #96897f"; $bcolor= "#625750";$fcolor="white";
	 }
 $stdout .= tee_stdout {
print qq{<hr class="colorgraph1" data-content="AND">};
my $div1 = qq{<div class="btn-group" role="group">};
my $style_b ="font-size: 14px;letter-spacing:2px;font-family: Gill Sans, Verdana;text-transform:uppercase; " ;
my $style_a ="font-size: 8px;font-family:Verdana;letter-spacing:1px;	padding-left: 10px;text-transform:uppercase; ";
print $cgi->start_div({class=>"btn-group  btn-toolbar-xs btn-group-justified",role=>"group",style=>"box-shadow: 0 6px 4px -4px black;"});

#

	print $cgi->start_div({class=>"btn-group btn-group-xs", role=>"group"});
		print start_button({class=>"btn btn-$h1",style=>$style_b.";background-color:#363945;border-width:0px;color:white;"});
			print "<img src='/icons/Polyicons/logo_polyweb_white2.png' width='140px' height='29px' style='margin:2px'><br>
			<span style='font-size: 7px;font-family:Verdana;position:relative;top:-1px;color:#F5DF4D;letter-spacing:normal'>- University of Paris<span style='color:#FF6F61'> - Imagine Institute -</span></span>
			";
			#print "<tr><td style='font-size: 6px;'>Imagine Institute ; Uni. of Paris</td></tr></table>";
		print end_button();
	print  $cgi->end_div;
	
	
	print $cgi->start_div({class=>"btn-group btn-group-xs", role=>"group",style=>"background-color:$ccolor"});
		print start_button({class=>"btn btn-$h2",style=>$style_b.";font-size:12px;border-width:0px;background-color:$ccolor".";color:white"});
		if ($project->isDefidiag){
			print "Polyviewer<BR><i>Defidiag</i>";
		}
		else {
			print "Polyviewer";
		}
		print end_button();
	print  $cgi->end_div;
	
	print $cgi->start_div({class=>"btn-group btn-group-xs", role=>"group",style=>"background-color:$ccolor"});
		print start_button({class=>"btn btn-$h2",style=>$style_b.";font-size:12px;border-width:0px;background-color:$ccolor".";color:white"});
			if ($project->isGenome){
					print "WGS<BR>$ver1"; 
			}
			elsif ($project->isExome){
					print "WES<BR><i>$ver1</i>"; 
			}
			else{
					print "Panel<BR><i>$ver1</i>"; 
			}
		print end_button();
	print  $cgi->end_div;
	
		print $cgi->start_div({class=>"btn-group btn-group-xs", role=>"group",style=>"background-color:$ccolor"});
			print start_button({class=>"btn btn-$h2",style=>$style_b.";font-size:12px;border-width:0px;background-color:$ccolor".";color:white"});
				print $project->name()."<BR> Samples: ".scalar(@{$project->getPatients});
			print end_button();
		print  $cgi->end_div;
		print $cgi->start_div({class=>"btn-group btn-group-xs", role=>"group",style=>"background-color:$bcolor"});
			print start_button({class=>"btn btn-$h2",style=>$style_b.";border-width:0px;background-color:$bcolor".";color:$fcolor"});
			print $project->description();
		print end_button();
	print  $cgi->end_div;
		print $cgi->start_div({class=>"btn-group btn-group-xs", role=>"group",style=>"background-color:$bcolor"});
		print start_button({class=>"btn btn-$h2",style=>$style_b.";font-size:12px;border-width:0px;background-color:$bcolor".";color:$fcolor"});
			print "$day-$month-$year";
		print end_button();
	print  $cgi->end_div;
	#USER
			print $cgi->start_div({class=>"btn-group btn-group-xs", role=>"group",style=>"background-color:$bcolor"});
		print start_button({class=>"btn btn-$h2",style=>$style_b.";font-size:12px;border-width:0px;background-color:$bcolor".";color:$fcolor"});
			print qq{<i class='fa fa-user'></i>&nbsp;$user};
		print end_button();
	print  $cgi->end_div;
	
	


	
print $cgi->end_div;
 };
  $no_cache->put_cache_text($cache_id,$stdout,2400) ;#unless $dev; 
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