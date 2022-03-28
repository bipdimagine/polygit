#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
#use lib "/software/polyweb/poly-src/GenBo/lib/obj-nodb/";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag";

#use lib "/bip-d/soft/distrib/tabix/latest/perl";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/coverage";
use lib "$Bin/../packages/validation_variation"; 
use lib "$Bin/../packages/cache"; 

use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use draw_cnv; 
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/utility.pm";
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
use List::MoreUtils qw{ natatime };
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
#use PDF::API2;
#use PDF::Table;
use constant mm => 25.4 / 72;
use constant in => 1 / 72;
use constant pt => 1;
use Time::HiRes qw ( time alarm sleep );



my $cgi          = new CGI();
my $project_name = $cgi->param('project');
html::print_header_polydiag($cgi);
my $url = url(-absolute=>1);

print qq{
	<form action="$url" method="post" id ="formid" name="formid">
    <textarea class="form-control" id="regions" cols="86" rows ="20"  name="regions" style="height:100px;width:300px">
    1:10525-129332
    2:1-43953
    3:64070-114243
    </textarea>
  
  	 <button type="submit" class="btn btn-primary">Submit</button>
	</form>
};
#1:10525-129332
#2:1-43953
#3:64070-114243
my $regions = $cgi->param('regions');
my @reg = split(" ",$regions);
	 my $tabix  = new Tabix(-data =>$Bin."/methyl_all.txt.gz");
	 
foreach my $r (@reg){
	my ($chr,$reste) = split(":",$r);
	my ($start,$end) = split("-",$reste);
	  my $res = $tabix->query($chr,$start,$end);
	 print qq{<div><br><br><span class="label label-danger">$r</span><br><div>};
	  	print $cgi->start_table({class=>"table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 10px;font-family:  Verdana;", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
	  	 while(my $line = $tabix->read($res)) {
	  	 	print qq{<tr>};
	  	 	my @l = split(" ",$line);
	  	 	foreach my $a (@l){
	  	 			print $cgi->td({style=>"text-align: center;"},$a);
	  	 	}
	  	 
	  	 	print qq{</tr>};
	  	 }
	  	 print $cgi->end_table();

}


exit(0);


