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
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
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
use File::stat;
#use Time::Piece;
use List::MoreUtils qw{ natatime };
use image_coverage;
use warnings;
use utf8;
use QueryValidationAcmg;
use CGI::Cache;
#use CGI::Carp;
use Digest::MD5::File qw(dir_md5_hex file_md5_hex url_md5_hex file_md5);
use Statistics::Descriptive;
use Number::Format;
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/utility.pm";
use feature ':5.10';
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
my $w = "500px";
my $format_number = new Number::Format;
my $VERSION_SCRIPT =1;
my $hsex = {
				2=>qq{  <i class="fa fa-venus fa-2x" > </i>},
				1=>qq{  <i class="fa fa-mars  fa-2x" > </i>} ,
				'-1'=>qq{  <i class="fa fa-minus  fa-2x" > </i>} 
};

	my $hsex1 = {
				2=>qq{  <i class="fa fa-venus fa-1x" > </i>},
				1=>qq{  <i class="fa fa-mars  fa-1x" > </i>} ,
				'-1'=>qq{  <i class="fa fa-minus  fa-1x" > </i>} 
};	
my $CSS = qq{
	<style type="text/css"> 
	\@import "compass/css3";
 

* {
  margin:0;
  padding:0;
}


html, body {
  height: 100%;
}


.table1 {
  background-color:#ffffff; 
  width:100%;

  
 
}

ul li {
  float:left;
  width:110px;
  margin-right: 5px;
  margin-left: 5px;
  margin-top: 5px;
  margin-bottom: 5px;
  text-align:center;

	
  list-style-type: none;
}




.top {
  background-color:#EAE9E4;
  height:35px;
   font-size:11px;
   border:1px solid black;
  color:black !important;
  .h11 {
    color:black !important;
    padding-top:10px;
  }
}

.circle {
  width:30px;
  height:30px;
  border-radius:30px;
  font-size:1px;
  color:#fff;
  line-height:30px;
  text-align:center;
  background:#ff;
  margin-left:73px;
  margin-top:-4px;
   box-shadow: 1px 1px 2px 1px rgba(0, 0, 0, .4);
}

.circle2M {
  width:30px;
  height:30px;
  border-radius:30px;
  font-size:1px;
  color:#fff;
  line-height:30px;
  text-align:center;
  background:#ff;
  margin-bottom:3px;
   box-shadow: 1px 1px 2px 1px rgba(0, 0, 255, .4);
}

.circle2F {
  width:30px;
  height:30px;
  border-radius:30px;
  font-size:1px;
  color:#fff;
  line-height:30px;
  text-align:center;
  background:#ff;
  margin-bottom:3px;
   box-shadow: 1px 1px 2px 1px rgba(255, 0, 0, 0.4);
}

.bottom {
margin-top:10px;
  
  p {
     padding:0px;
    span {
     //font-weight:bold; 
    }
  }
}

.p1 {
     padding:0px;
     color:black;
  
      
     margin-bottom: 1px;
     border-bottom-style: solid;
     border-bottom-width: 1px;
    span {
     font-weight:bold; 
    }
  }
  
  .p1danger {
     padding:0px;
     color:white;
      background-color:#DD4124;
      //background-color:#E15D44;
     margin-bottom: 1px;
     border-bottom-style: solid;
     border-bottom-width: 1px;
    span {
     font-weight:bold; 
    }
  }
  
.p1 warning {
     padding:0px;
     color:black;
     margin-bottom: 1px;
     border-bottom-style: solid;
     border-bottom-width: 1px;

    span {
     //font-weight:bold; 
    }
  }
.blue1 {
    box-shadow: 1px 1px 2px 1px rgba(0, 0, 255, .4);
}
.pink1 {
    box-shadow: 1px 1px 2px 1px rgba(255, 0, 0, .4);
}



.blue {
  background-color:#C6E7F6;   
    box-shadow: 1px 1px 2px 1px rgba(0, 0, 0, .4);
}

.white {
 color:#FFFFFF; 
}

.duplicate {
  background-color:#FFD662;
  box-shadow: 1px 1px 2px 1px rgba(0, 0, 0, .4);
}
.duplicate1 {
  box-shadow: 1px 1px 2px 1px rgba(255, 132, 14, .4);
}
.pink {
  background-color:#FFEBED;
  box-shadow: 1px 1px 2px 1px rgba(0, 0, 0, .4);
}


.danger {
  background-color:#DD4132;
  box-shadow: 1px 1px 3px 2px rgba(0, 0, 0, 0.4);
}
.danger1 {
    box-shadow: 1px 1px 3px 2px rgba(255, 0, 0, 0.7);
}

.success {
  background-color:#00A591;
   background-color:#577284;
   color:white;
  box-shadow: 1px 1px 3px 2px rgba(0, 0, 0, 0.4);
}
.primary {
  background-color:#00A591;
   background-color:#F0EDE5;
   color:white;
  box-shadow: 1px 1px 3px 2px rgba(0, 0, 0, 0.4);
}

.success1 {
  //  box-shadow: 1px 1px 3px 2px rgba(0, 255, 0,0.7);
}

.aler {
  background-color:#E15D44;
  background-color:#FA9A85;
  box-shadow: 1px 1px 3px 2px rgba(0, 0, 0, 0.4);
}
.aler1 {
    box-shadow: 1px 1px 3px 2px rgba(250, 154,133, 0.7);
}

.warning {
  background-color:#FFD662;
  box-shadow: 1px 1px 3px 2px rgba(0, 0, 0, 0.4);
}


.warning1 {
    box-shadow: 1px 1px 3px 2px rgba(255, 214,98, 0.7);
}

.colorgraph {
  height: 3px;
  border-top: 0;
  background: #c4e17f;
  border-radius: 5px;
  background-image: -webkit-linear-gradient(left, #c4e17f, #c4e17f 12.5%, #f7fdca 12.5%, #f7fdca 25%, #fecf71 25%, #fecf71 37.5%, #f0776c 37.5%, #f0776c 50%, #db9dbe 50%, #db9dbe 62.5%, #c49cde 62.5%, #c49cde 75%, #669ae1 75%, #669ae1 87.5%, #62c2e4 87.5%, #62c2e4);
  background-image: -moz-linear-gradient(left, #c4e17f, #c4e17f 12.5%, #f7fdca 12.5%, #f7fdca 25%, #fecf71 25%, #fecf71 37.5%, #f0776c 37.5%, #f0776c 50%, #db9dbe 50%, #db9dbe 62.5%, #c49cde 62.5%, #c49cde 75%, #669ae1 75%, #669ae1 87.5%, #62c2e4 87.5%, #62c2e4);
  background-image: -o-linear-gradient(left, #c4e17f, #c4e17f 12.5%, #f7fdca 12.5%, #f7fdca 25%, #fecf71 25%, #fecf71 37.5%, #f0776c 37.5%, #f0776c 50%, #db9dbe 50%, #db9dbe 62.5%, #c49cde 62.5%, #c49cde 75%, #669ae1 75%, #669ae1 87.5%, #62c2e4 87.5%, #62c2e4);
  background-image: linear-gradient(to right, #c4e17f, #c4e17f 12.5%, #f7fdca 12.5%, #f7fdca 25%, #fecf71 25%, #fecf71 37.5%, #f0776c 37.5%, #f0776c 50%, #db9dbe 50%, #db9dbe 62.5%, #c49cde 62.5%, #c49cde 75%, #669ae1 75%, #669ae1 87.5%, #62c2e4 87.5%, #62c2e4);
}
/* ----------------------------
 * msg
 * ----------------------------
 */
.msg {
    background: #fefefe;
    color: #666666;
    font-weight: bold;
    font-size: small;
    padding: 6px;
    margin-right: 10px;
    padding-left: 6px;
    border-top: solid 3px #CCCCCC;
    border-radius: 5px;
    margin-bottom: 1px;
    -webkit-box-shadow: 0 10px 10px -5px rgba(0,0,0,.08);
       -moz-box-shadow: 0 10px 10px -5px rgba(0,0,0,.08);
            box-shadow: 0 10px 10px -5px rgba(0,0,0,.08);
}
.msg-clear {
    border-color: #fefefe;
    -webkit-box-shadow: 0 7px 10px -5px rgba(0,0,0,.15);
       -moz-box-shadow: 0 7px 10px -5px rgba(0,0,0,.15);
            box-shadow: 0 7px 10px -5px rgba(0,0,0,.15);
}
.msg-info {
    border-color: #b8dbf2;
}
.msg-y {
    border-color: #FECF71;
}
.msg-rp {
    border-color: #BD3D3A;
}
.msg-green {
    border-color: #BFD641;
}
.msg-purple {
    border-color: #BE9EC9;
}
.msg-coral {
    border-color: #FF6F61;
}

.msg-violet {
    border-color: #6B5B95;
}
.msg-brown {
    border-color: #6C4F3D;
}


/* ----------------------------
 * Callouts
 * ----------------------------
 */




.callout-mage h1,
h2,
h3,
h4 {
	font-weight: 300;
	line-height: 1.4;
}

.callout-mage {
	padding: 15px;
	background-color: #1079B2;
	color: #fff;
	 box-shadow: 1px 1px 3px 2px rgba(0, 0, 0, 0.4);
}
.callout-mage h1 {
     //margin-left: -20px;
     margin-top: 0;
     position: relative;
     width: 90%;
    color: white;
    text-shadow: 0 2px 0 black;
    font-size:20px;
}


.mcontainer {
  display: flex; /* or inline-flex */
  flex-direction: row;
  justify-content:center;
}



/*STAMP */
.stamp1 {

	color: black;
	font-size: 2.1rem;
	font-weight: 700;
	border: 0.25rem solid #555;
	display: inline-block;
	padding: 0.25rem 1rem;
	text-transform: uppercase;
	border-radius: 1rem;
	font-family: 'Courier';
	-webkit-mask-image: url('https://s3-us-west-2.amazonaws.com/s.cdpn.io/8399/grunge.png');
  -webkit-mask-size: 944px 604px;
  mix-blend-mode: multiply;
  color: #BC243C;
	border: 0.5rem solid #BC243C;
	-webkit-mask-position: 13rem 6rem;
	transform: rotate(0deg);
  border-radius: 0;
}
.stamp2 {

	color: black;
	font-size: 1.2rem;
	font-weight: 500;
	border: 0.25rem solid #555;
	display: inline-block;
	padding: 0.25rem 1rem;
	text-transform: uppercase;
	border-radius: 1rem;
	font-family: 'Courier';
	-webkit-mask-image: url('https://s3-us-west-2.amazonaws.com/s.cdpn.io/8399/grunge.png');
  -webkit-mask-size: 944px 604px;
  mix-blend-mode: multiply;
  color: #BC243C;
	border: 0.5rem solid #BC243C;
	-webkit-mask-position: 13rem 6rem;
	transform: rotate(0deg);
  border-radius: 0;
}

.stamp3 {

	color: black;
	font-size: 1.2rem;
	font-weight: 500;
	border: 0.25rem solid #555;
	display: inline-block;
	padding: 0.25rem 1rem;
	text-transform: uppercase;
	border-radius: 1rem;
	font-family: 'Courier';
	-webkit-mask-image: url('https://s3-us-west-2.amazonaws.com/s.cdpn.io/8399/grunge.png');
  -webkit-mask-size: 944px 604px;
  mix-blend-mode: multiply;
  color:#4B0082;
	border: 0.2rem solid #4B0082;
	-webkit-mask-position: 13rem 6rem;
	transform: rotate(0deg);
  border-radius: 0;
}

	</style>
};

#CGI::Cache::setup();

my $buffer = GBuffer->new();

 my $fsize = "font-size:10px";
my $cgi          = new CGI();
my $project_name = $cgi->param('project');

my $project = $buffer->newProjectCache( -name 			=> $project_name, -typeFilters=>'individual' ,-cgi_object=>1);



$project->validations_query(1);
my $user = $cgi->param('user_name');
$project->cgi_user($user);





 if ($cgi->param('xls')  ){
	#die();
	#html::print_cgi_header($cgi,$CSS,1,$project_name." - Variations Editor");
	#print $CSS;
	print "Content-type: application/msexcel\n";
	print "Content-Disposition: attachment;filename=".$project->name().".xls\n\n";
	table_patients_xls($project);
	exit(0);
}


#my @array_control_panels = ("control_mendel","control_quality","control_sex","control_duplicate");
my $list_control_panels = "control_design,control_mendel,control_quality,control_sex,control_duplicate,control_muc1";


#html::print_cgi_header($cgi,$CSS,1,$project_name." - Variations Editor");
if ($cgi->param('print')  ){
	#die();
	html::print_cgi_header($cgi,$CSS,1,$project_name." - Variations Editor");
	print $CSS;
	table_patients_printer($project);
	exit(0);
}

html::print_header_polydiag($cgi);
print $CSS;



my $hgmd = $buffer->queryHgmd->getHGMD($user);
$hgmd = 0 unless $hgmd;
#my $project = $buffer->newProject(-name=>$project_name);




#die();

my $t = time;


 $t = time;
my $hv;
my $stat; 

my $dev;
#$dev = 1 if $ENV{SERVER_NAME} eq  "10.200.27.103";


 $t = time;
 $|=1;
 my $key_quality = args_quality($project);
 
 #warn $key_quality;
 # my $key_quality = args_quality($project);
 my $no_cache = $project->get_lmdb_cache_summary("r");
 push(@$key_quality,"muc1") if $project->getCaptures->[0]->analyse =~ /renom/i;
 my $header = $no_cache->get_cache(join(";",@$key_quality).".header");
 warn "HEADER ==> ".$header;
 warn join(";",@$key_quality);
 
 $header = undef  if $dev or $cgi->param('force') == 1;
 #$header = undef;
 $no_cache->close();
$project->getChromosomes();
$project->getPatients();
my $hmendel;
my ($gstats,$lstats,$patient_value) ;
 warn "HEADER ==> ".$header;
 my $cache_icon ="";
 unless ($header){
 	print qq{<div style="display: none">};
 	($gstats,$lstats,$patient_value) = statistics_projects($project);
	 $hmendel = get_mendelian_statistics($project);
	foreach my $run (@{$project->getRuns()}) {
		$header->{$run->id} =  table_run_header($run);
 	}
 	print qq{</div>};
 	$no_cache = $project->get_lmdb_cache_summary("w");
 	my $id = join(";",@$key_quality);
 	$no_cache->put_cache_hash($id.".header",$header,2400); 
 	 $no_cache->close();
 	  $cache_icon .= qq{<span class="glyphicon glyphicon-floppy-remove" aria-hidden="true" style="text-align:right;color:red"></span>};
 }
 else {
 	  $cache_icon .= qq{<span class="glyphicon glyphicon-floppy-saved" aria-hidden="true" style="text-align:right;font-size:10px;color:green"></span>};
 	
 }



#my $key_validation = args_validation($project);
 #push(@$key_quality,"muc1") if $project->getCaptures->[0]->analyse =~ /renom/i;
$no_cache = $project->get_lmdb_cache_summary("r");
my $htable = $no_cache->get_cache(join(";",@$key_quality).".table");
#$htable = undef;
 $no_cache->close();
$htable = undef if $dev or $cgi->param('force') == 1;;

unless ($htable){
		print qq{<div style="display: none">};
 	 ($gstats,$lstats,$patient_value) = statistics_projects($project) unless $patient_value;
	 $hmendel = get_mendelian_statistics($project) unless $hmendel;
	foreach my $run (@{$project->getRuns()}) {
		$htable->{$run->id} =  table_patients($run);
 	}
 	$no_cache = $project->get_lmdb_cache_summary("w");
 	my $id = join(";",@$key_quality);
 
 	$no_cache->put_cache_hash($id.".table",$htable,2400); 
 	 $no_cache->close();
 	print qq{</div>};
 	 $cache_icon .= qq{&nbsp; <span class="glyphicon glyphicon-floppy-remove" aria-hidden="true" style="text-align:right;font-size:10px;color:red"></span>};
}
else {
 	  $cache_icon .= qq{<span class="glyphicon glyphicon-floppy-saved" aria-hidden="true" style="text-align:right;font-size:10px;color:green"></span>};
 	
 }


 $| =1;
	
		
print qq{<center><div class="row" style="height:auto;width:98%;">};


print qq{<div class="col-sm-12">};
foreach my $run (@{$project->getRuns()}) {
	print $header->{$run->id};
	print $htable->{$run->id};
}
print"<br><div style='float:right;'><small>$cache_icon</small></div><br>";
print qq{</div>};

#print qq{<div class="col-sm-4">};
#print qq{</div>};

print qq{</div></center>};
exit(0);

 

 

	sub return_date {
		my ($dd) = @_;
		 my @amonths = ('Jan', 'Feb', 'Mar', 'Apr','May',"Jun","Jul","Aug","Sep","Oct","Nov","Dec");
		my ($date,$time) = split(" ",$dd);
	    my ($year,$month,$day) =  split("-",$date);
		return ($year,$amonths[$month-1],$day);
	}
	
 sub get_array {
 	my ($vector) = @_;
 	my $set = Set::IntSpan::Fast->new($vector->to_Enum);
 	return [$set->as_array()];
 }	
	
sub check_level {
	my ($data) = @_;
	my $level = "success";
	foreach my $v (@{$data->{data}}){
		my @zs = values %$v;
		foreach my $z (@zs){
			if ($z->{type} eq "danger"){
				$level = "danger";
				last;
			}
			if ($z->{type} eq "warning"){
				$level = "warning";
			}
		}
		last if $level eq "danger";
	}
 return $level;
}


 
 sub print_line_patient {
 	my ($p,$nb_row) = @_;
 		 #@title = ("Fam","Patient","Sex","mendelian error","Cov","30x","Sub","Indels","%He","validation");
 		 my $line;
 		 my $sex_eval = $p->compute_sex(); 
 		my $pname = "check_".$p->name();
			my $class ={};
			my $hval = $p->validations();
			my $hval_cnv =  $p->validations_cnv;
			if (keys %$hval){
		#	 $nb_row += scalar(keys %$hval);
			}
			$class = {rowspan=>$nb_row} if $nb_row >1;
			$class->{style}= "vertical-align:middle";
			my $fam = $p->getFamily();
			my $ir = $p->return_icon." ".$p->name;
			my $cov_sry = $p->coverage_SRY();
				#$line->{Fam}  = $cgi->td($class,$p->name);
				
				$line->{Sex}  =  $cgi->td($class,$p->return_icon);
			
			
			#print $cgi->td($class,$hsex->{$sex});
			my $style_btn_name= qq{style ="background-color:#FFEBED;color:black"}  ;
			$style_btn_name= qq{style ="background-color:#C6E7F6;color:black"}  if $p->sex() == 1;
#			warn $sex_eval ." ".$p->sex();
#			die();
			if ($sex_eval ne $p->sex() && $sex_eval ne -1){
				# $class->{class}= "danger";
				  $style_btn_name= qq{style ="background-color:#E74C3C;color:white"};
				 
			} 
			
			
				
			#$line->{Patient}  = $cgi->td($class,$p->return_icon."&nbsp;".$p->name);
			
			my $cov = $p->coverage();
			#global value
			
			my $sex = $p->sex(); 
			if ($sex_eval ne $sex && $sex_eval ne -1){
				# $class->{class}= "danger";
			} 
				#$line->{Fam}  = $cgi->td($class,$p->name);
				
				$line->{Sex}  =  $cgi->td($class,$p->return_icon);
			
			
			#print $cgi->td($class,$hsex->{$sex});
			 $sex = $p->sex(); 
			if ($sex_eval ne $sex && $sex_eval ne -1){
			#	 $class->{class}= "danger";
			} 
			$line->{Sex} = $cgi->td($class,$p->return_icon);
			$line->{Sex2}  = $cgi->td($class,$hsex->{$sex_eval}."<small>(".$cov_sry.")</small>");
			
		
			
			
			 my $v = $cov->{mean};
			 $v= 0 unless $v;
			 my $btn_class = qq{class= "btn btn-xs btn-success " style = "$fsize" };
			$btn_class = qq{class= "btn  btn-xs btn-warning " style = "$fsize" } if $v < 15;
			$btn_class = qq{class= "btn  btn-xs btn-alert " style = "$fsize" } if $v < 10;
			$line->{Cov}  = $cgi->td($class,qq{<button type="button" $btn_class >$v</button>});
			$v = $cov->{'30x'};
			 $v= 0 unless $v;
			$line->{"30x"}  = $cgi->td($class,qq{<button type="button" $btn_class >$v</button>}); $cgi->td($class,qq{<button type="button" $btn_class >$v %</button>});
		
			my $style ={};

			my $hstatus = $p->getLatestValidationStatus($user);
  		
			
			#print $cgi->td($class,$td);
		#	my $mtime = stat("/etc/passwd")->mtime; 
				$line->{"validation"}  = $cgi->td($class,"-");

			if (exists $hstatus->{term} && $hstatus->{term} eq "negative"){
				my $term = $hstatus->{term};
				my $date2 ="";
				my $xx;
				($date2,$xx) = split(" ",$hstatus->{modification_date});
				my $u = $hstatus->{user_name}; 
				my $text = qq{ <span  class="stamp1"><span>-$term-</span>&nbsp;-&nbsp;<small>$date2</small>&nbsp-&nbsp;$u</span>};
			
				 $line->{"validation"}  =  $cgi->td($text);
				
			}
			else {
				my @text;
				if (keys %$hval){
					push(@text,validation_table_new($p,$hval));
			 
				}
			
				if (keys %$hval_cnv){
					push(@text,validation_table_cnv($p,$hval_cnv));
				}
				if(@text) {
					 $line->{"validation"}  = $cgi->td({style=>"background-color:#E74C3C;padding:3px"},join("<BR>",@text));
				}
				
			
			}
			my $patname = $p->name();
			
			
			
			$line->{Patient} = $cgi->td($class,qq{<button type="button" class ="btn btn-xs btn-primary " onclick="select_patient_interface_in_summary('$patname');"  $style_btn_name >$ir</button>});
			
			$line->{"PolySplice"} = $cgi->td({style=>"vertical-align:middle;"}, "<a type='button' class='btn btn-xs' style='background-color:white;border: solid 1px red;' disabled><span style='color:red;' class='glyphicon glyphicon-remove'></span></button></a>");
		
			$line->{"PolyViewer"} = $cgi->td({style=>"vertical-align:middle;"}, "<a type='button' class='btn btn-xs' style='background-color:white;border: solid 1px green' onclick=\"select_patient_interface_in_summary('$patname');\"><span style='color:green;' class='glyphicon glyphicon-ok'></span></button></a>");
			if ($project->isDiagnostic) {
				$line->{"Variations Editor"} = $cgi->td({style=>"vertical-align:middle;"}, "<a type='button' class='btn btn-xs' style='background-color:white;border: solid 1px green' onclick=\"select_patient_interface_in_summary('$patname','polydiag');\"><span style='color:green;' class='glyphicon glyphicon-ok'></span></button></a>");
			}
			else {
				$line->{"Variations Editor"} = $cgi->td({style=>"vertical-align:middle;"}, "<a type='button' class='btn btn-xs' style='background-color:white;border: solid 1px red;' disabled><span style='color:red;' class='glyphicon glyphicon-remove'></span></button></a>");
			}
			if ($project->isGenome()) {
				$line->{"PolyDupDel"} = $cgi->td({style=>"vertical-align:middle;"}, "<a type='button' class='btn btn-xs' style='background-color:white;border: solid 1px red;' disabled><span style='color:red;' class='glyphicon glyphicon-remove'></span></button></a>");
				$line->{"PolyCyto"} = $cgi->td({style=>"vertical-align:middle;"}, "<a type='button' class='btn btn-xs' style='background-color:white;border: solid 1px green' onclick=\"select_patient_interface_in_summary('$patname','polycyto');\"><span style='color:green;' class='glyphicon glyphicon-ok'></span></button></a>");
			}
			else {
				$line->{"PolyDupDel"} = $cgi->td({style=>"vertical-align:middle;"}, "<a type='button' class='btn btn-xs' style='background-color:white;border: solid 1px green' onclick=\"select_patient_interface_in_summary('$patname','polydupdel');\"><span style='color:green;' class='glyphicon glyphicon-ok'></span></button></a>");
				$line->{"PolyCyto"} = $cgi->td({style=>"vertical-align:middle;"}, "<a type='button' class='btn btn-xs' style='background-color:white;border: solid 1px red;' disabled><span style='color:red;' class='glyphicon glyphicon-remove'></span></button></a>");
			}
			
			my ($hasJunctions, $hasJunctionsAndDNA);
			$hasJunctions = 1 if ($p->hasJunctions());
			if ($hasJunctions) {
				$line->{"PolySplice"} = $cgi->td({style=>"vertical-align:middle;"}, "<a type='button' class='btn btn-xs' style='background-color:white;border: solid 1px green' onclick=\"select_patient_interface_in_summary('$patname','polysplice');\"><span style='color:green;' class='glyphicon glyphicon-ok'></span></button></a>");
				my $nb_j = 0;
				foreach my $chr (@{$project->getChromosomes()}) {
					$nb_j += $chr->countThisVariants($p->getJunctionsVector($chr));
					my $v_dna = $chr->getCategoryVariantsVector('substitution')->Clone();
					$v_dna += $chr->getCategoryVariantsVector('insertion')->Clone();
					$v_dna += $chr->getCategoryVariantsVector('deletion')->Clone();
					$hasJunctionsAndDNA = 1 if ($chr->countThisVariants($v_dna) > 0);
					last if ($hasJunctionsAndDNA);
				}
				$line->{"Nb Junctions"} = $cgi->td($class,qq{<button type="button" $btn_class >$nb_j</button>});
				if (not $hasJunctionsAndDNA) {
					$line->{"PolyViewer"} = $cgi->td({style=>"vertical-align:middle;"}, "<a type='button' class='btn btn-xs' style='background-color:white;border: solid 1px red;' disabled><span style='color:red;' class='glyphicon glyphicon-remove'></span></button></a>");
					$line->{"Variations Editor"} = $cgi->td({style=>"vertical-align:middle;"}, "<a type='button' class='btn btn-xs' style='background-color:white;border: solid 1px red;' disabled><span style='color:red;' class='glyphicon glyphicon-remove'></span></button></a>");
					$line->{"PolyDupDel"} = $cgi->td({style=>"vertical-align:middle;"}, "<a type='button' class='btn btn-xs' style='background-color:white;border: solid 1px red;' disabled><span style='color:red;' class='glyphicon glyphicon-remove'></span></button></a>");
					$line->{"PolyCyto"} = $cgi->td({style=>"vertical-align:middle;"}, "<a type='button' class='btn btn-xs' style='background-color:white;border: solid 1px red;' disabled><span style='color:red;' class='glyphicon glyphicon-remove'></span></button></a>");
				}
			}
			
			##
			my $res_muc = $p->vntyperResults();
			if (@$res_muc){
				my $text = qq{ <span  class="stamp1"><span>MUC1</span></span>};
			 	$line->{"MUC1"}  = $cgi->td({style=>"vertical-align:middle"},"$text");
			}
			elsif (-e  $project->getVariationsDir("vntyper")."/muc1/" && !(-e $p->vntyperTsv)){
				my $text = qq{ <span  class="stamp1"><span>PROBLEM !!!!!</span></span>};
				$line->{"MUC1"}  = $cgi->td({style=>"vertical-align:middle"},"$text");
			}
			else {
				my $text = qq{ <span  class="stamp3"><span>NONE</span></span>};
				$line->{"MUC1"}  = $cgi->td({style=>"vertical-align:middle"},"$text");
				# $line->{"MUC1"}  = $cgi->td({style=>"vertical-align:middle"},"-");
				 #$line->{"MUC1"}  = $cgi->td({style=>"vertical-align:middle"},"-");
			}
			
			
			return $line;
 }


sub header_run {
	my ($run,$nb_run,$total_run) = @_;
	$nb_run =1;
	$total_run = 1;
	
		my $stats = {};
		 my $info ;
		push(@$info , $run->description());
		
			#	$info .=$run->getCapture()->type;
		
	#$info .= $run->version;
 	 my $cno = $project->getChromosomes()->[0]->lmdb_hash_variants("r");
 	
 	 my $version;
 	 $version->{gencode} = $project->gencode_version();
	$version->{gnomad} = $buffer->description_public_lmdb_database("gnomad-exome")->{version};
	$version->{hgmd} = $buffer->description_public_lmdb_database("hgmd")->{version};
	$version->{cadd} = $buffer->description_public_lmdb_database("cadd")->{version};
	$version->{clinvar} = $buffer->description_public_lmdb_database("clinvar")->{version};
	
	 my $date_cache =  utility::return_date_from_file($cno->filename);
	 
#
	my $captures;
	foreach my $c (@{$run->getPatients}){
		my $capt = $c->getCapture();
		$captures->{$capt->id} = $capt;
#		my $cov = $c->coverage();
#		foreach my $k (keys %$cov){
#				$stats->{$k} += $cov->{$k};
#		}
	}
	my $bcolor = qq{background-color:#1079B2};
	$bcolor = "background-color:#BD3D3A" if $project->isSomatic;
	my ($date1,$hour) = split(" ",$run->date);
	my($y,$m,$d) = split("-",$date1);
	$y =~s/20//;
	my $date = "$d/$m/$y"; 
	my $capt = values %$captures;
	foreach my $c (values %$captures){
			push(@$info , $project->description);
			push(@$info , $c->type);
	}
	
	 my $fcolor = "white";
	 my $text_info ="";
		$text_info .= qq{SOMATIC<span class="glyphicon glyphicon-minus" aria-hidden="true"></span> } if $project->isSomatic;
		 $text_info .= join(qq{<span class="glyphicon glyphicon-minus" aria-hidden="true"></span>},@$info);
		$text_info .= qq{<span class="glyphicon glyphicon-minus" aria-hidden="true"></span>SOMATIC } if $project->isSomatic;
	my $nb = scalar(@{$run->getPatients});
		 
		#$btn_class = qq{class= "label  label-xs label-warning " style="font-size: 12px;" } if int($stats->{"30x"}/$nb)  < 95;
		#$btn_class = qq{class= "label  label-xs label-alert "  style="font-size: 12px;" } if int($stats->{"30x"}/$nb)  < 90;
		#warn $nb;
		#die();
		 my $text30 = '';# = "30X :" .$format_number->round($lstats->{"30X"}->{$run->id}->mean,1)."% <i>(".$format_number->round($gstats->{"30X"}->mean,1).")</i>";;
					 
		 my $text100 = '';# = "100X :".$format_number->round($lstats->{"100X"}->{$run->id}->mean,1)."%<i>(".$format_number->round($gstats->{"100X"}->mean,1).")</i>" if  $lstats->{"100X"}->{$run->id} ;
		
		  my $text15 = '';# = "15X :" .$format_number->round($lstats->{"15X"}->{$run->id}->mean,1)."%<i>(".$format_number->round($gstats->{"15X"}->mean,1).")</i>";;
		  my $textmean = '';# = "Cov :" .$format_number->round($lstats->{"mean"}->{$run->id}->mean,1)." <i>(".$format_number->round($gstats->{mean}->mean,1)."&nbsp;&#177;&nbsp;". $format_number->round($gstats->{mean}->standard_deviation,0).")</i>";
		  
 		 my $capture = $project->getCaptures()->[0];
 		my $ver1 = $project->genome_version;
 		my $ver = $project->annotation_version;
 		my $machine = $run->machine();
# 		my $capture  = qq{<span class="label label-info " style="font-size: 12px;">}. $capt->name."</span>" ;

	 
 my $texta;
 foreach my $name (keys %$version){
			my $v = $version->{$name};
			$texta .= qq{<span class="label " style="font-size: 12px;$bcolor">$name <span class="badge badge-alert" style="font-size: 10px;">$v </span></span>	};
		}
		
	my $cmd =qq{project_printer("$project_name",$user)};
my  $out = qq{
  	<div  style="position: relative;">
       	$texta
       </div>
       <div class="callout-mage text-center fade-in-b" style="padding-bottom: 10px;padding-top: 5px;$bcolor">
       
       	
					<h1>
					$text_info
            		</h1>
				
        <div class="mcontainer" center>
  			
    		
				 <div class="col-sm msg msg-y" > 
				 <i class="fa fa-simplybuilt"></i>&nbsp;$machine
				 </div>
				 <div class="col-sm msg msg-rp" > 
				 <i class="fa fa-calendar"></i>&nbsp;$date
				 </div>
				 <div class="col-sm msg msg-green" > 
				 <i class="fa fa-laptop"></i> &nbsp;$date
				 </div>
				 <div class="col-sm msg msg-violet" > 
				 <i class="fa fa-line-chart"></i>&nbsp;$textmean
				 </div>
				 <div class="col-sm msg msg-purple" > 
				 $text15
				 </div>
   				<div class="col-sm msg msg-coral" > 
				 $text30
				 </div>
				<div class="col-sm msg msg-brown" > 
				  $ver1
				 </div>
				
				 <div class="col-sm msg msg-green" onClick="project_printer(\'$project_name\',\'$user\')" > 
				 <img src="https://img.icons8.com/officel/30/000000/send-to-printer.png">
					 Print
				 </div>
   			</div> 
 
};
;

if ($project->isDiagnostic){
	my @color = ("#C4E17F","#F7FDCA","#FECF71","#F0776C","#DB9DBE","#C49CDE");
	my @color = ("aliceblue","#DDDDDD");
	$out.= qq{<div>};
	$out .= $cgi->start_table({class=> "table table-striped table-bordered ",style=>"color:black;text-align: center;vertical-align:middle;font-size: 10px;font-family:  Verdana;box-shadow: 1px 1px 2px 1px YYY;margin-top: 5px;margin-bottom: 0px;"});
	
my $n=0;
foreach my $m (@{$project->callingMethods}) {
	my $c = $color[$n%(scalar(@color))];
	$out.= "<td style=\"background-color:$c\">".$m."</td>\n";
	$n ++;
}
$out.= qq{</table></div>};
}
$out.= qq{</div>};

$out.= qq{
   <hr class="colorgraph" style="margin-top:1px;margin-bottom:1px">
 };
 return $out;
}


sub table_mendelian {
	
	 my ($run) = @_;
	 my $out;
	
	 my $fams = $run->getFamilies;
	 
	  	$out = $cgi->start_div({class=>"panel-body panel-primary panel-collapse collapse  ",style=>"font-size: 09px;font-family:  Verdana;border: 5px coral; border-color: coral;",id=>"control_mendel_".$run->id});	
	   $out .= qq{
		<div class="container">
		<div class="row">
	   };
	   #	 	<div class="wrap">
#		<div class="table1">
#	
	   my $nb;
	   my $greport = {};
	 foreach my $fam (@$fams){
	 
	 	next unless  $fam->isTrio();
	 	$greport->{fam} ++;
	 	my $report= {};
	 	 my $error = 0;
	 	my $warning = 0;
	 		$nb ++;
	 		my $name = $fam->name;
	 		my $c1 = "pink";
	 		my $c2= $c1."1";
	 		my $father = $fam->getFather;
	 		
	 		my $mother = $fam->getMother;
	 		
	 		$out .= qq{
	 	 <div class="col-xs-6 col-md-4 col-lg-3">
			
	 		};
	 	$out .= $cgi->start_table({class=> "table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 10px;font-family:  Verdana;box-shadow: 1px 1px 2px 1px YYY ;", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
	 	
	 
		#	<table align="center" style="border: 1px  solid black;">
		$out .= $cgi->start_Tr({style=>"font-weight: bold;XXX "});
		#	$out .= $cgi->start_Tr({style=>"background-color:red;color:white; "});
		my $img = qq{<img src="https://img.icons8.com/offices/24/000000/family.png">};
				$out .= qq{
					<td colspan="2"> $img $name </td>
					</tr>
				};
				
	 		
	 		$out .= $cgi->start_Tr();
	 		
	 		my ($text,$w,$e);
	 		 $out .=  parent_box($father,"F",$report);
	 		$out .=  parent_box($mother,"M",$report);
	 		
	 		$out .= $cgi->end_Tr();
	 		foreach my $c (@{$fam->getChildren}){
	 				$out .= $cgi->start_Tr();
	 			$out .= children_box($c,$report);
	 					$out .= $cgi->end_Tr();
	 			
	 		}
	 	
	 		$out .= qq{</table></div>};
	 	#last if $nb ==2;
	 	
	 	
	 	 my $style = "background-color:#CAECD7;color:black; ";
	 	  $style = "background-color:#1079B2;color:white; ";
	 	  $style = "background-color:#00A591;color:white; ";
	 	  	my $shadow = "rgba(121,199,83, .4)";
	  if (exists $report->{warning} ){
	  		$style = "background-color:#FE840E;color:black; ";
	  	  	$greport->{warning} ++;
	  	  	$shadow = "rgba(254, 132,14, .6)";
	  }
	  if (exists $report->{error}) {
	  	$style = "background-color:#F2552C;color:black; ";
	  	$style = "background-color:#DD4132;color:black;";
	  		$greport->{error} ++;
	  		$shadow = "rgba(221, 65,36, .7)";
	 
	  }
	  $out =~ s/XXX/$style/;
	 	  $out =~ s/YYY/$shadow/;
	 	
	 	}
	$out .= qq{
	 	</div>
		</div>
		</div>
	   };
	  
	  my $out2;
	   my $disabled = "disabled";
	 	my $btn = "";
	 	my $text = qq{  <img src="https://img.icons8.com/material-sharp/20/000000/genealogy.png"> Mendelian Control };
	 	my $label;
	 	 if (exists $greport->{fam}){
	 	 	 	 $btn = " btn-success";
	 	 	  	$text = qq{ <img src="https://img.icons8.com/material-sharp/22/000000/genealogy.png">  Mendelian Control };
	 	 	  $label = qq{	<span class="badge badge-success">0</span>};
	 	 	  $disabled ="";
	 	 }
	 	 my $error = 0;
	 	 my $warning =0;
	 	  $label ="";
	 	 
	 	if (exists $greport->{warning}){
	  	$text = qq{ <img src="https://img.icons8.com/material-sharp/20/000000/genealogy.png"> Mendelian control };
	  	 $btn = " btn-warning";
	  	 $warning = $greport->{warning};
	  	 $label = qq{	<span class="badge badge-warning">$warning</span>};
	 }
	 if (exists $greport->{error}){
	  	$text = qq{ <img src="https://img.icons8.com/material-sharp/20/000000/genealogy.png"> Mendelian control };
	  	 $btn = " btn-danger";
	  	 $error = $greport->{error};
	  	 $label = qq{	<span class="badge badge-danger">$error</span>};
	 }
	 my $run_id = $run->id;
	
	  $out2 =  qq{<div class="btn   btn-xs $btn " style="position:relative;bottom:1px;min-width:200px;border-color:black;" onClick='collapse_panel("control_mendel","$list_control_panels","$run_id")'  style="border : 1px"  $disabled>  $text $label </div>};
	  
	 return ($out2,$out);
	 
	 }


sub table_design {
	
	 my ($run) = @_;
	 #return ("","");
	 return ("","") unless $project->isDiagnostic();
	 my $out;
	 my $capture = $run->getCapture();
	 my  @transcripts_cgi = @{$project->bundle_transcripts() } ;
	 my $problem;
	 my $nbp;
	 my $capture_intspan_keep;
	 foreach my $tr_name  (@transcripts_cgi){
	 	my $transcript =  $project->newTranscript($tr_name);
	 	my $chr = $transcript->getChromosome;
	 	my $capture_intspan = $transcript->getChromosome->getIntSpanCapture(100);
	 	unless (exists $capture_intspan_keep->{$chr->name}){
	 		$capture_intspan_keep->{$chr->name} = $transcript->getChromosome->getIntSpanCapture();
	 	}
	 	
	 	foreach my $exon (sort{$a->start <=> $b->start} @{$transcript->getExons}){
	 		$capture_intspan_keep->{$chr->name}->remove_range($exon->start-1000,$exon->end+1000);
	 		next if $exon->is_noncoding;
	 		my $s1 = $exon->getGenomicSpan()->intersection($capture_intspan);
	 		 if ($s1->is_empty){
	 		 push (@{$problem->{problem}->{$transcript->name}},$exon);
	 		 $nbp ++;
	 		 }
	 		
	 		 
	 		#$covered = -1 if $s1->is_empty;
	 	}
	 }
	 
	 
	 my $hgenes;
	 foreach my $chr_name (keys %{$capture_intspan_keep}){
	 	my $iter = $capture_intspan_keep->{$chr_name}->iterate_runs();
#	 	warn $capture_intspan_keep->{$chr_name}->as_string();
	 	my $chr = $project->getChromosome($chr_name);
		my @tt;
    	while (my ( $from, $to ) = $iter->()) {
    		my $genes = $chr->getGenesByPosition($from,$to);
    		foreach my $g (@$genes){
    			$hgenes->{$g->id}->{obj} = $g;
    			$hgenes->{$g->id}->{nb} ++;
    		 	#warn $g->external_name;
    		}
    	}
    }
    
	 	

	  my 	$style = "background-color:#DD4132;color:black;";
	my	  		$shadow = "rgba(221, 65,36, .7)";
	 	$out = $cgi->start_div({class=>"panel-body panel-primary  panel-collapse collapse ",style=>"font-size: 09px;font-family:  Verdana;border: 5px coral; border-color: coral;",id=>"control_design_".$run->id});	
	   $out .= qq{
		<div class="container">
		<div class="row">
	   };
	   my @gs = grep{$hgenes->{$_}->{nb}>1} keys %{$hgenes};
	#   @gs =();
	   if (@gs){
	   	$nbp ++;
	   	$out .= qq{
	 	 <div class="col-xs-6 col-md-4 col-lg-3">
			
	 		};
	 		$out .= $cgi->start_table({class=> "table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 10px;font-family:  Verdana;box-shadow: 1px 1px 2px 1px YYY ;", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
	 		$out .= $cgi->start_Tr({style=>"font-weight: bold;XXX "});
	 		$out .= qq{
					<td colspan="2"> Potential missing genes  </td>
					</tr>
				};
				
		#my @gs = grep{$hgenes->{$_}->{nb}>1} keys %{$hgenes};
		my $max =0;
	   foreach my $gid (@gs){
	   			$out .= $cgi->start_Tr();
				$out .= qq{<td colspan="2">}.$hgenes->{$gid}->{obj}->external_name." miss : ".$hgenes->{$gid}->{nb}." exons </td>"; 
				$max = $hgenes->{$gid}->{nb} if $hgenes->{$gid}->{nb} > $max;
				#$out .= qq{<td colspan="2">}.$exon->name." ".$exon->getChromosome->name.":".$exon->start."-".$exon->end."</td>"; 
				$out .= $cgi->end_Tr();
	   }
	   $out .= qq{</table></div>};
	   }
	   
	   
	#$out = $cgi->start_div({class=>"panel-body panel-primary panel-collapse collapse  ",style=>"font-size: 09px;font-family:  Verdana;border: 5px coral; border-color: coral;",id=>"control_design_".$run->id});

		foreach my $tid (keys %{$problem->{problem}}){
			
			my $transcript = $problem->{problem}->{$tid}->[0]->getTranscript();
			
				
	 		$out .= qq{
	 	 <div class="col-xs-6 col-md-4 col-lg-3">
			
	 		};
	 		$out .= $cgi->start_table({class=> "table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 10px;font-family:  Verdana;box-shadow: 1px 1px 2px 1px YYY ;", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
	 		$out .= $cgi->start_Tr({style=>"font-weight: bold;XXX "});
	 		my $gname = $transcript->getGene->external_name;
	 		my $tname = $transcript->name;
	 		$out .= qq{
					<td colspan="2"> $gname :  $tname </td>
					</tr>
				};
			foreach my $exon (@{$problem->{problem}->{$tid}}){
				$out .= $cgi->start_Tr();
				$out .= qq{<td colspan="2">}.$exon->name.": &nbsp;[".$exon->getChromosome->ucsc_name.":".$exon->start."-".$exon->end."]</td>"; 
				#$out .= qq{<td colspan="2">}.$exon->name." ".$exon->getChromosome->name.":".$exon->start."-".$exon->end."</td>"; 
					$out .= $cgi->end_Tr();
			}
				$out .= qq{</table></div>};
		}
		$out =~ s/XXX/$style/g;
	 	$out =~ s/YYY/$shadow/g;
	    $out .= qq{
	 	</div>
		</div>
		</div>
	   };	  
		

	
	
	  my $out2;
	  my $disabled = "disabled";
	  my $btn = "";
	  my $text = qq{  <img src="https://img.icons8.com/material-sharp/20/000000/genealogy.png"> Panel design };
	 my $label;
	 $nbp += 0;
	 if ( $nbp > 0){
	 	$disabled = "";
	  	 $btn = " btn-danger";
	 }
	 	$label = qq{	<span class="badge badge-danger"> $nbp</span>};
	 my $run_id = $run->id;
	
	  $out2 =  qq{<div class="btn   btn-xs $btn " style="position:relative;bottom:1px;min-width:200px;border-color:black;" onClick='collapse_panel("control_design","$list_control_panels","$run_id")'  style="border : 1px"  $disabled>  $text $label </div>};
	 return ($out2,$out);
	 
 		
	 }


sub parent_box {
	my ($p,$stype,$r) = @_;
		my $out;
		my $class = "circle2M";
		 $class = "circle2F" if $stype eq "M";
		
	unless ($p) {
		 my $icon = qq{<div class ="$class"> <img src="https://img.icons8.com/ios/32/000000/decision.png"> </div>};
		$out .= qq{
 			
 			<td align="center">  $icon   </td>
 		
 			
	 		};
	 		return $out;
	}
		return "" unless $p;
		my $v =$hmendel->{$p->name};
		my $type = "btn-success";
		$v = "?" if  $v eq "";
		if ($v eq "?"){
			$type = "btn-danger";
 			$r->{error}++;
		}
 		elsif ($v>3) {
 			$type = "btn-danger";
 			$r->{error}++;
 		}
 		elsif ($v >1){
 			$type = "btn-warning";
 			$r->{warning}++;
 		}
	
		 my $btn_class = qq{class= "btn btn-xs $type " style="font-size:10px;padding-left: 1px;padding-right: 1px;"};
		 	my $name = $p->name;
		  my $value = qq{<button type="button" $btn_class  > $name - $v <span>&#37;</span> </button>};
		 my $icon = qq{<div class ="$class">}. $p->return_icon."  </div>";
		$out .= qq{
 			
 			<td align="center">  $icon $value </td>
 		
 			
	 		};
	 		return $out;
}
	 	
sub children_box {
	my ($p,$r) = @_;
		my $class = "circle2M";
		my $v =$hmendel->{$p->name};
		 $class = "circle2F" if $p->sex ne 1;
		 	my $type = "btn-success";
			$v = "?" if  $v eq ""; 
			if ($v eq "?"){
			$type = "btn-danger";
 			$r->{error}++;
		}	
 		elsif ($v>3) {
 			$type = "btn-danger";
 			$r->{error}++;
 		}
 		elsif ($v >1){
 			$type = "btn-warning";
 			$r->{warning}++;
 		}
 		
		 my $btn_class = qq{class= "btn btn-xs $type"  style="font-size:10px" };
		 	my $name = $p->name;
		  my $value = qq{<button type="button" $btn_class > $name - $v <span>&#37;</span> </button>};
		my $out;
	
		my $icon = qq{<div class ="$class">}. $p->return_icon."  </div>";

		$out .= qq{
 				<td align="center" colspan="2"> $icon  $value</td>};
 			return $out;
 				
}

sub get_mendelian_statistics {
	my ($project) = @_;
my $no = $project->noSqlQuality("r");
my $data = $no->get($project->name,"mendelian");
my $hmendel = {};
return $hmendel unless $data;
#warn scalar @{$data->{data}};
#die();
foreach my $line (@{$data->{data}}){
	next unless  $line;
	my $name = $line->{sample}->{text};
	my $patient = $project->getPatient($name);
		my $value =  $line->{mendelian}->{text};
	$value =~ s/\[//;
		$value =~ s/\]//;
	my ($t,$p,$n) = split(" ",$value);
	$p =~ s/%//;
	if ($patient->isParent) {
	
		
		 $p = $format_number->round($p/scalar(@{$patient->getFamily->getChildren}),1) if $p;
	} 

	$hmendel->{$name} = $p;

#	$hmendel->{$name} = $line->{mendelian}->{text};
	
}
	return $hmendel;
}

sub statistics_projects {
	my ($run) = @_;
	 my $capture = $project->getCaptures()->[0];
	my $query = $project->buffer->getQuery->listAllProjectsNameByCaptureId($capture->id());
	
	my $stats;
	my $nb = 0;
	my $this_stats;
	my $gstat;
	my $lstat;	
	my $patient_value;
	push(@$query,$project->name()) unless $query;

	my $nbpp =0;
	if (scalar(@$query) > 0){
		foreach my $p (@$query){
			next  if $nbpp > 10 && $p ne $project->name;
			next if $p ne $project->name;
			my $buffer2 = GBuffer->new();;
			my $project2 =  $buffer2->newProject( -name 			=> $p);
			next unless $project2->existsnoSqlQuality;
			my $no = $project2->noSqlQuality("r");
			foreach my $v ("coverage_stats","statistics_variations"){
			my $data = $no->get($project2->name,"$v");
		     
				next unless $data;
					
					$nbpp ++;
				
				foreach my $p_value (@{$data->{data}}){
					$nb ++;
					#my $patient = $p_value->{patients}->{text};
					my $patient = $project2->getPatientOrControl($p_value->{patients}->{text});
					next if $patient->is_control;
					foreach my $k (keys %$p_value){
						
						next if $k eq "patients";
						next unless $p_value->{$k}->{text};
						$p_value->{$k}->{text} =~s/%//;
						$patient_value->{$patient->id}->{$k} = $p_value->{$k}->{text};

						push(@{$stats->{$k}},$p_value->{$k}->{text});
					#	warn  $p." ".$project->name; 
						if ($p eq $project->name){
								
								push(@{$this_stats->{$k}->{$patient->getRun->id}},$p_value->{$k}->{text});
						}
					
					}
				
				
				}
			}
		}
	}
	unless ($this_stats){
		$gstat->{"15X"} = Statistics::Descriptive::Full->new();
		$gstat->{"30X"} = Statistics::Descriptive::Full->new();
			$gstat->{"100X"} = Statistics::Descriptive::Full->new();
			return ($gstat,$lstat,$patient_value);
		
	}
	
	
	
	foreach my $k (keys %$stats){
		
		my $stat = Statistics::Descriptive::Full->new();
		$gstat->{$k} = Statistics::Descriptive::Full->new();
	
		$gstat->{$k} ->add_data(@{$stats->{$k}});
		
	}
	
	foreach my $k (keys %$this_stats){
		foreach my $r (keys %{$this_stats->{$k}} ) {
			
			$lstat->{$k}->{$r} = Statistics::Descriptive::Full->new();
			$lstat->{$k}->{$r} ->add_data(@{$this_stats->{$k}->{$r}});
		
		}
	}
	
	return ($gstat,$lstat,$patient_value);
	
	
	} 


sub table_quality {
	
	 my ($run) = @_;

	 my $out;
	 my $error = 0;
	 
	  $out = $cgi->start_div({class=>"panel-body panel-primary panel-collapse collapse ",style=>"font-size: 09px;font-family:  Verdana;border: 5px coral; border-color: coral;",id=>"control_quality_".$run->id});	
	 $out .= qq{
	 	<div class="wrap">
<div class="table1">
<ul>
	 };
	 
	 my $sum = {};
	 my $nb = scalar(@{$run->getPatients});

	 
	  foreach my $p (sort{$a->name cmp $b->name } @{$run->getPatients}){
	  	
	  	my $icon = $p->return_icon;
	  	my $name = $p->name;
	  	my $c1 = "XXX";
	  	my $cov = $p->coverage();
	  	my $mean = $cov->{mean};
			my $c2= $c1."1";

	  	$out .= qq{
  <li class="$c2">
  <div class="top $c1">
    <h11>$name</h11>
    <div class="circle">$icon</div>
  </div>
  <div class="bottom" style="text-align:left">
	  	};
	  	my $level = 0;
	  	my $red = 0;
	  	$red =1 if ($patient_value->{$p->id}->{"15X"} and ($patient_value->{$p->id}->{"15X"}*1.0 < 95));
	  		#warn " red $red" if $p->name eq "GEF1800528";
	  	foreach my $v ("mean","15X","30X","100X","snp","indel","%he","%public") {
	  		next unless $lstats->{$v}->{$run->id};
	  		my $value = $patient_value->{$p->id}->{$v};
	  		$value =0 unless $value;
	  		my $mean  = int($gstats->{$v}->mean);
	  		$mean = 0 unless $mean;
	  		my $std  = $gstats->{$v}->standard_deviation;
	  		my $p1;
	  		$p1 = "p1";
	  		my $color = "#009B77";
	  		if ($v eq "%he") {
	  			
	  			if ($value < 50) {
	  				$color = "#FFD662";
	  				 $level = 1 if $level<1;
	  			}
	  			if ($value >75) {
	  				$color = "#FFD662";
	  				 $level = 1 if $level<1;
	  			}
	  			if ($value >80) {
	  				$color = "#FA9A85";
	  				 $level = 2 if $level<2;
	  			}
	  			if ($value >85) {
	  				$color = "#DD4132";
	  				 $level = 3;# if $level<3;
	  			}
	  		}
	  		elsif ($v eq "15X" or $v eq "30X" or $v eq "100X" or $v eq "mean" or $v eq "%public" ) {
	  			#warn $v." $value $mean ::  $std". ($mean - 3*$std)." $red" if $p->name eq "GEF1800528";
	  			if ($value <  ($mean - 3*$std) && $red == 1)  {
	  			$color = "#DD4132";
	  			 $level = 3 if $level<3;
	  		} 
	  			
	  		elsif ($value <  ($mean - 2*$std) && $red ==1 )  {
	  			#$p1 = "p1danger";
	  			$color = "#FA9A85";
	  			 $level = 2 if $level<2;
	  		} 
	  		elsif ($value <  ($mean - 2*$std) )  {
	  			#$p1 = "p1danger";
	  			$color = "#FA9A85";
	  			 $level = 1 if $level<1;
	  		} 
	  		}
	  		
	  		else {
	  		if ($value <  ($mean - 5*$std) || $value >( 5*$std + $mean))  {
	  		
	  			#$p1 = "p1danger";
	  			$color = "#DD4132";
	  			 $level = 3 if $level<3;
	  		} 
	  		if ($value <  ($mean - 4*$std) || $value >( 4*$std + $mean))  {
	  			#warn $value ." ".$std." ".$mean;
	  			#$p1 = "p1danger";
	  			$color = "#DD4132";
	  			 $level = 3 if $level<3;
	  		} 
	  			 
	  		elsif ($value <  ($mean - 3*$std) || $value >( 3*$std + $mean))  {
	  			#$p1 = "p1danger";
	  			$color = "#FA9A85";
	  			 $level = 2 if $level<2;
	  		} 
	  		
	  		}
	  		
	  		my $lmean = $format_number->round($lstats->{$v}->{$run->id}->mean*1.0,1);
	  		$out .= qq{<p class="$p1"><i class="fa fa-circle" style="color:$color;margin-right: 5px;margin-left: 2px; "></i>$v <span> $value <span> <span><i>($lmean)</i> </span>};
	  		
	  	}

     $out .= qq{
  </div>
</li>
	  	};
	  	my $switch = "primary";
	  	
	  	
	  	if ($level ==1) {
	  		 $error = 1 if $error<1;
	  		$switch = "warning";
	  	}
	  	elsif ($level ==2) {
	  		 $error = 2 if $error<2;
	  		$switch = "aler";
	  	}
	  	elsif ($level ==3) {
	  		 $error = 3 if $error<3;
	  		$switch = "danger";
	  	} 
	  	$out =~ s/XXX/$switch/g;
	  }
	  	
$out .= qq{	  	
</ul>
</div>
</div>
</div>
	 };
	 my $out2;
	 	my $btn = " btn-primary";
	 	my $text = qq{ <img src="https://img.icons8.com/flat_round/20/000000/warranty-card.png"> Quality Control};
	 	my $color = "#00A591";
	if ($error == 1 ){
		 $color = "#FFD662";
	} 	
	elsif ($error ==2) {
		 $color = "#FA9A85";
	}
	elsif ($error ==3) {
		 $color = "#DD4132";
	}
	 	
	 if ($error > 0){
	  	$text = qq{  <img src="https://img.icons8.com/flat_round/20/000000/warranty-card.png"> Quality Control };
	  	 $btn = " btn-danger";
	 }
	  my $run_id = $run->id;
	  $out2 =  qq{<div class="btn  btn-info btn-xs " style="position:relative;bottom:1px;min-width:200px;;border-color:black;background-color:$color" onClick='collapse_panel("control_quality","$list_control_panels","$run_id")'>  $text <span class="badge badge-danger">$error</span> </div>};
	 return ($out2,$out);
		 
}

sub construct_identito_vigilence {
	my ($p) = @_;
	  my @iv = split("",$p->identity_vigilance());
	
	   return ("",0) unless @iv;
	   
	   unless ($p->identity_vigilance_vcf()){
	   	system("$Bin/../../polypipeline/scripts/scripts_pipeline/identito_vigilence.pl -project=".$project->name.">/dev/null 2>/dev/null");
	   }
	   my @iv_vcf ;
	   @iv_vcf = split("",$p->identity_vigilance_vcf()) if $p->identity_vigilance_vcf();
	   my $level;
	  unless (@iv_vcf){
	  		 for (my $i=0;$i<@iv;$i++){
	  		 	$iv_vcf[$i] = "";
	  		 }
	  }
	  elsif (scalar(@iv_vcf) >  scalar(@iv)) {
	  		$level = pop(@iv_vcf);
	  }
	  my $out;
	  my $error = 0;
	    
	 for (my $i=0;$i<@iv;$i++){
	 	if ($iv[$i] eq $iv_vcf[$i]){
	 		$out.=$cgi->span({style=>"background-color:#00A591;color:white;font-size:10px"},$iv[$i]);
	 	}
	 	elsif ($iv[$i] eq "4"){
	 		$out.=$cgi->span({style=>"background-color:#6C4F3D;color:white;font-size:10px"},$iv_vcf[$i]);
	 	}
	 	elsif ($iv_vcf[$i] eq ""){
	 		
	 		$out.=$cgi->span({style=>"background-color:grey;color:white;font-size:10px"},$iv[$i]);
	 	}
	 	else {
	 		$error ++;
	 		$out.=$cgi->span({style=>"background-color:#E94B3C;color:white;font-size:10px"},$iv_vcf[$i]);
	 	}
	 	
	 }	
	 if ($level){
	 	$error = 0;
		$error = 2 if $level eq "e";
	  	$error = 1 if $level eq "w";
	 }
	 else {
	 	$error = 2 if $error >=2;
	 }
	 
	
	 return ($out,$error);
	
	
}

sub table_sex {
	
	 my ($run) = @_;
	 my $out;
	 my $error = 0;
	
	  $out = $cgi->start_div({class=>"panel-body panel-primary panel-collapse collapse ",style=>"font-size: 09px;font-family:  Verdana;border: 5px coral; border-color: coral;",id=>"control_sex_".$run->id});	
	 $out .= qq{
	 	<div class="wrap">
<div class="table1">
<ul>
	 };
	 
	  foreach my $p (sort{$a->name cmp $b->name } @{$run->getPatients}){
	  	print ".";
	  	my $icon = $p->return_icon;
	  	my $name = $p->name;
	  	my $c1 = "pink";
	  	$c1 = "blue" if $p->sex ==1;
	  	my $cov = $p->coverage();
	  	my $mean = $cov->{mean};
	  	my $x30 = $cov->{"30x"};
	  	my $cov_sry = $p->coverage_SRY();
	  	 my $sex_eval = $p->compute_sex(); 
	  	# warn $p->name.' '.$p->compute_sex.' '.$p->coverage_SRY();
	  	 my $color = "#009B77";
	  	 if ($sex_eval ne $p->sex() && $sex_eval ne -1){
				# $class->{class}= "danger";
				 # $style_btn_name= qq{style ="background-color:#E74C3C;color:white"};
				  $c1 = "danger";
				  $error ++;
	  				$color = "#DD4132";
			} 
			my $c2= $c1."1";
	  	my $text2 = qq{<i class="fa fa-circle" style="color:$color;margin-right: 5px;margin-left: 2px; "></i>}.$hsex1->{$sex_eval};
	  	my $iv ="";
	  	my $iverror = 0;
	  		   ($iv,$iverror) = construct_identito_vigilence($p) if $p->identity_vigilance;
	  	if ($iverror == 2){
	  			  $c1 = "danger";
				  $error ++;
	  	} 
	  	if ($iverror == 1){
	  			  $c1 = "warning";
				  $error ++;
	  	} 
	  	
	  	my $text3 = "<small>(".$cov_sry.")</small>";
	  	$out .= qq{
  <li class="$c2">
  <div class="top $c1">
    <h11>$name</h11>
    <div class="circle">$icon</div>
  </div>
  <div class="bottom">
    <p>SRY : <span>$text2</span> $text3</p>
	<p> IV <span> $iv </p>
  </div>
</li>
	  	};
	  	
	  }
	  	
$out .= qq{	  	
</ul>
</div>
</div>
</div>
	 };
	 my $out2;
	 	my $btn = " btn-primary";
	 	my $text = qq{ <img src="https://img.icons8.com/material/14/000000/mars-symbol.png">/<img src="https://img.icons8.com/material/14/000000/venus-symbol.png"> Gender control };
	 if ($error > 0){
	  	$text = qq{  <img src="https://img.icons8.com/material/14/000000/mars-symbol.png">/<img src="https://img.icons8.com/material/14/000000/venus-symbol.png">Gender control };
	  	 $btn = " btn-danger";
	 }
	  my $run_id = $run->id;
	  $out2 =  qq{<div class="btn  btn-info btn-xs $btn" style="position:relative;bottom:1px;min-width:200px;;border-color:black" onClick='collapse_panel("control_sex","$list_control_panels","$run_id")'>  $text <span class="badge badge-danger">$error</span> </div>};
	 return ($out2,$out);
		 
}

sub table_control {
	 my ($run) = @_;
	 my $controls = $run->getPatientsControl;
	 	my $samtools = $buffer->software("samtools");
	 	unless (@$controls) {
	 		
	 			my $out1 =  qq{<div  type ="button"  style="position:relative;bottom:1px;min-width:200px;border-color:black;display: none" class="btn  btn-xs"  disabled> <img src="https://img.icons8.com/ios/20/000000/empty-test-tube.png"> Control (Blanc) &nbsp;&nbsp</div>};
	 			
	 			return ($out1,"");
	 	}
	 	
	 	
	#return "" unless $controls;
	if (@$controls){
		my $out_table ="";
		my $style = "success";
		
		foreach my $c (@$controls){
		$out_table .= $cgi->td($c->name);
		
		my $cov = $c->coverage();
	
		my $bf = $c->getBamFileName;
		my $reads;
		my $c1;
		my $nb = 0;
		foreach my $patient (@{$project->getPatients}){
		my $bam = $patient->getBamFile;	
		my $res = ` $samtools idxstats $bam | awk '{sum+=\$4;sum2+=\$3} END {print sum";"sum2}'`;
		chomp($res);
		my @z = split(";",$res);
		$reads += ($z[0]+$z[1]);
		$c1 += $patient->coverage()->{mean};
		$nb ++;
		last if $nb>10;
		}

		my $meancov = int($c1/$nb);
		my $meanread = int($reads/$nb);
		
		
		unless (-e $bf){
			$out_table .= $cgi->td('-');
		}
		else {
		my $bam = $c->getBamFile;	
		my $res = ` $samtools idxstats $bam | awk '{sum+=\$4;sum2+=\$3} END {print sum";"sum2}'`;
		chomp($res);
		my @z = split(";",$res);
		my $readscontrol = ($z[0]+$z[1]);
		my $meancontrol = $c->coverage()->{mean};
		my $c15x = $c->coverage()->{"15x"};
		my $c30x = $c->coverage()->{"30x"};
		my $v =int( ($meancontrol/$meancov)*1000)/10;
	
		my $r = int( ($readscontrol/$meanread)*1000)/10;
		if ($v > 0.25){
			$style= "warning";
		}
		if ($r > 1){
			$style= "warning";
		}
		if ($v> 0.5){
			$style= "danger";
		}
		
		if ($r > 2){
			$style= "danger";
		}
		$out_table .= $cgi->td($readscontrol);
		$out_table .= $cgi->td($meanread);
		$out_table .= $cgi->td($r."%");
		$out_table .= $cgi->td($meancontrol);
		$out_table .= $cgi->td($meancov);
		$out_table .=  $cgi->td($v."%");
		$out_table .= $cgi->td($c15x."%");
		$out_table .= $cgi->td($c30x."%");
		my $p = int($z[1]/($z[0]+$z[1])*1000)/10;
		$out_table .= $cgi->td($p."%");
		
		}
	
			}	
	  
	my $out1;	
	my $nb = scalar(@$controls);
	 my $run_id = $run->id;
	$out1 =  qq{<div class="btn  btn-info btn-xs btn-$style" style="position:relative;bottom:1px;min-width:200px;border-color:black" onClick='collapse_panel("control_panel","$list_control_panels","$run_id")'> <img src="https://img.icons8.com/ios/20/000000/empty-test-tube.png"></span>Control (Blanc)  &nbsp;&nbsp;<span class="badge badge-info">$nb</span></div>};
	my $out;
	$out .= $cgi->start_div({class=>"row"});
	 my $pstyle = "panel-primary";#.$style;
	$out .= $cgi->start_div({class=>"panel-body panel-primary panel-collapse collapse ",style=>"font-size: 09px;font-family:  Verdana;",id=>"control_panel_".$run->id});	
	$out .= $cgi->start_div({class=>"col-xs-6"});
	$out .= $cgi->start_table({class=> "table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 10px;font-family:  Verdana;", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
	$out .= $cgi->start_Tr({class=>"$style"});
	$out .= $cgi->th( {style=>"text-align: center;"} ,"name") ;
	$out .= $cgi->th( {style=>"text-align: center;"} ,"Nb Reads Control") ;
	$out .= $cgi->th( {style=>"text-align: center;"} ,"Nb Reads Project") ;
	$out .= $cgi->th( {style=>"text-align: center;"} ,"Ratio %") ;
	$out .= $cgi->th( {style=>"text-align: center;"} ,"Coverage Control") ;
	$out .= $cgi->th( {style=>"text-align: center;"} ,"Coverage Project") ;
	$out .= $cgi->th( {style=>"text-align: center;"} ,"Ratio %") ;
	$out .= $cgi->th( {style=>"text-align: center;"} ,"15x") ;
	$out .= $cgi->th( {style=>"text-align: center;"} ,"30x") ;
	$out .= $cgi->th( {style=>"text-align: center;"} ,"% aligned reads") ;
	$out .= $cgi->end_Tr();
	$out .= $out_table;
	$out .= $cgi->end_table();
	$out .= $cgi->end_div();
	$out .= $cgi->end_div();
	$out .= $cgi->end_div();
	
	return ($out1,$out);
}
}

sub table_muc1 {
	 my ($run) = @_;
	 my $controls = $run->getPatientsControl;
	 my $patients ;
	 	
	#return "" unless $controls;
		my $out_table ="";
		my $style = "success";
		
	 my $nb =0;
		foreach my $patient (@{$project->getPatients}){
			my $t = $patient->vntyperResults;
			unless (@{$t}){
				if(-e $patient->vntyperTsv()){
			 	$out_table .= $cgi->start_Tr({class=>""});
			 	$out_table .= $cgi->td($patient->name);
			 	$out_table .= $cgi->td(["-","-","-","-","-","-","-","-","-","-","-"]);
			 	#,"-","-","-","-","-","-","-","-","-","-","-"]);
				$out_table .= $cgi->end_Tr({class=>""});
				}
				else {
				$out_table .= $cgi->start_Tr({class=>""});
			 	$out_table .= $cgi->td([$patient->name,"PROBLEM","PROBLEM","PROBLEM","PROBLEM","PROBLEM","PROBLEM","PROBLEM","PROBLEM","PROBLEM","PROBLEM"]);
				$out_table .= $cgi->end_Tr({class=>""});
				}
			 }
			next unless @{$t};
			$out_table .= $cgi->start_Tr({class=>""});	
			$out_table.= $cgi->td({rowspan=>scalar(@$t)},$patient->name);
			foreach my $tt (@$t){
				$out_table .= $cgi->td($tt);
				$out_table .= $cgi->end_Tr({class=>""});
			}
			#$out_table .= $cgi->end_Tr({class=>""});	
			$nb ++;
		}
	
	my $out1;	
	#my $nb =  1;
	my $run_id = $run->id;
	$out1 =  qq{<div class="btn  btn-info btn-xs btn-$style" style="position:relative;bottom:1px;min-width:200px;border-color:black;background-color:#C49CDE;color:black" onClick='collapse_panel("control_muc1","$list_control_panels","$run_id")'> <img src="https://img.icons8.com/fluency-systems-filled/20/null/biotech.png"/></span>MUC1 &nbsp;&nbsp;<span class="badge badge-info">$nb</span></div>};
	unless (-e  $project->getVariationsDir("vntyper")."/muc1/" ){
		return ("","");
	}
	my $out;
	$out .= $cgi->start_div({class=>"row"});
	 my $pstyle = "panel-primary";#.$style;
	$out .= $cgi->start_div({class=>"panel-body panel-primary panel-collapse collapse ",style=>"font-size: 09px;font-family:  Verdana;",id=>"control_muc1_".$run->id});	
	$out .= $cgi->start_div({class=>"col-xs-6"});
	$out .= $cgi->start_table({class=> "table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 10px;font-family:  Verdana;", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
	$out .= $cgi->start_Tr({class=>"$style"});
	
	$out .= $cgi->th( {style=>"text-align: center;"} ,"name") ;
		$out .= $cgi->th( {style=>"text-align: center;"} ,"Motif") ;
	$out .= $cgi->th( {style=>"text-align: center;"} ,"Variant") ;
	$out .= $cgi->th( {style=>"text-align: center;"} ,"POS") ;
	$out .= $cgi->th( {style=>"text-align: center;"} ,"REF") ;
	$out .= $cgi->th( {style=>"text-align: center;"} ,"ALT") ;
	$out .= $cgi->th( {style=>"text-align: center;"} ,"Motif_sequence (RevCom)") ;
	$out .= $cgi->th( {style=>"text-align: center;"} ,"Estimated_Depth_AlternateVariant") ;
	$out .= $cgi->th( {style=>"text-align: center;"} ,"Estimated_Depth_Variant_ActiveRegion") ;
	$out .= $cgi->th( {style=>"text-align: center;"} ,"Depth_Score") ;
	$out .= $cgi->th( {style=>"text-align: center;"} ,"Confidence") ;
	$out .= $cgi->th( {style=>"text-align: center;"} ,"Kmer") ;
	$out .= $cgi->end_Tr();
	$out .= $out_table;
	$out .= $cgi->end_table();
	$out .= $cgi->end_div();
	$out .= $cgi->end_div();
	$out .= $cgi->end_div();
	
	return ($out1,$out);
}

sub table_duplicate {
	 my ($run) = @_;
	 my $hash_dup;
	 foreach my $patient (@{$run->getPatients}){

 		my $filebed =  $patient->project->getVariationsDir("duplicate_region_calling")."/regions/".$patient->name().".dup.bed";
		if (-e $filebed){
	  	open (BED,$filebed);
	  	while(<BED>){
	  		chomp();
	  		my ($chr,$start,$end) = split(" ");
	  		unless (exists $hash_dup->{$chr} ){
	  			$hash_dup->{$chr} = Set::IntSpan::Fast::XS->new();
	  		}
	  	
	  		$hash_dup->{$chr}->add_range($start,$end);
	  	}
		}
 
 	}
 
 my $he;
 my @warnings;
 my $ws;
  my $htr ;
   unless ($project->bundle_transcripts()){
   	my $out1 = qq{<div class="btn  btn-xs btn-dark "  style="position:relative;bottom:1px;min-width:200px;;border-color:black" disabled>  Regions Dups  &nbsp;&nbsp;<span class="badge badge-dark">0</span></div>};
	return ($out1,"");
   }
 	map{$htr->{$_} ++} @{$project->bundle_transcripts() } ;
 	
 	
 	
 	
 	
 	
 	
 foreach my $cn (keys %$hash_dup){
 	my $chr = $project->getChromosome($cn);
 	my $iter = $hash_dup->{$cn}->iterate_runs();
    while (my ( $from, $to ) = $iter->()) {
    	
 	my $ts  = $chr->getTranscriptsByPosition($from,$to);
 	
 	foreach my $t (@$ts){
 		next unless exists $htr->{$t->id};
 		foreach my $e (@{$t->getExons}){
 			next if exists $he->{$e->name."_".$t->id};
 			
 			my $span = $hash_dup->{$cn}->intersection($e->getGenomicSpan);
 			my $l1 = scalar($e->getGenomicSpan->as_array);
 			my $l2 = scalar($span->as_array);
 			my $p = int (($l2/$l1 *100));
 			unless ($span->is_empty){
 				 my $i1 = $span->iterate_runs();
 				 my @st;
 				  while (my ( $from, $to ) = $i1->()) {
 				  	push(@st,"[".abs($from-$e->start)."-".abs($to-$e->start)."]");
 				   }
 				 $he->{$e->name."_".$t->id} ++;
 				 $ws->{$t->getGene->external_name()}->{tr}->{$t->name}->{$e->name}->{pos} = join(" ",@st);
 				 $ws->{$t->getGene->external_name()}->{tr}->{$t->name}->{$e->name}->{percent} =$p."%";
 				  $ws->{$t->getGene->external_name()}->{row} ++;
 				  
 			 	# my $warning = $t->getGene->external_name()." ".$t->name." ".$e->name." ".$span->as_string." ".$e->getGenomicSpan->as_string;
 			}
 			  
 			}
 		}
    }
 }
my $out;
my $out1;
my $nb = scalar( keys %{$ws});

if ($nb ==0) {
	$out1 = qq{<div class="btn  btn-xs btn-dark "  style="position:relative;bottom:1px;min-width:200px;;border-color:black" disabled>  Regions Dups  &nbsp;&nbsp;<span class="badge badge-dark">0</span></div>};
	return ($out1,"");
}
 my $run_id = $run->id;	 
$out1 .= qq{<div class="btn  btn-xs btn-danger" style="position:relative;bottom:1px;min-width:200px;;border-color:light-grey; background-color:#FFD662;color:black" onClick='collapse_panel("control_duplicate","$list_control_panels","$run_id")'>   <i class="fa fa-bell" aria-hidden="true"></i> Regions Dups  &nbsp;&nbsp;<span class="badge badge-danger">$nb</span></div>};
$out = $cgi->start_div({class=>"panel-body panel-primary panel-collapse collapse  ",style=>"font-size: 09px;font-family:  Verdana;border: 5px coral; border-color: coral;",id=>"control_duplicate_".$run->id});		
#$out .= $cgi->start_div({class=>"panel-body panel-collapse collapse panel-alert",style=>"font-size: 09px;font-family:  Verdana;",id=>"control_duplicate"});
	
	
	
	 $out .= qq{
	 	<div class="wrap">
<div class="table1">
<ul>
	 };
	 
	foreach my $g (keys %{$ws}){
	  foreach my $t (keys %{$ws->{$g}->{tr}}){
	  	  	my $c1 = "duplicate";
	  	  	my $c2= $c1."1";
	  	$out .= qq{
  <li class="$c2">
  <div class="top $c1">
    <h11>$g</h11>
    <h11>$t</h11>
 
  </div>
  <div class="bottom">
	  	};
	  $out .= qq{ <p > <table class=" table-bordered"> };
		foreach my $e (keys %{$ws->{$g}->{tr}->{$t}}){  	
				my $p = $ws->{$g}->{tr}->{$t}->{$e}->{percent};
			my $pos = $ws->{$g}->{tr}->{$t}->{$e}->{pos};
			  $out .= qq{ <tr>};
			   $out .= $cgi->td([$e,$p,$pos]);
			     $out .= qq{ </tr>};
		
			
   			# $out .= qq{ <p  	style="border-bottom-style: solid;border-bottom-width: 1px;padding-bottom: 1px; border-color:light-grey;color:black" ><div class="btn  btn-xs btn-danger">$e  <span>$p</span> $pos </div></p>};
    	}
    
 $out .= qq{ </table></p >  };
	  	}
	  	
	  	$out .= "  </div></li>";#end transcript
	  }
	  	
$out .= qq{	  	
</ul>
</div>
</div>
</div>
	 };
	
	
	 return ($out1,$out);
	
	
$out .= $cgi->start_table({class=>"table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 10px;font-family:  Verdana;", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
$out .= qq{<thead>};
		$out .= $cgi->start_Tr({class=>"danger"});
		$out .= $cgi->th( {colspan=>"5",style=>"text-align: center;"} ,"Region Dups");
	$out .= $cgi->end_Tr();
		$out .= $cgi->start_Tr({class=>"success"});
		$out .= $cgi->th( {style=>"text-align: center;"} ,"Gene");
		$out .= $cgi->th( {style=>"text-align: center;"} ,"Transcript");
		$out .= $cgi->th( {style=>"text-align: center;"} ,"exon");
		$out .= $cgi->th( {style=>"text-align: center;"} ,"% dup");
		$out .= $cgi->th( {style=>"text-align: center;"} ,"position");
	$out .= $cgi->end_Tr();
		
foreach my $g (keys %{$ws}){
	$out .= $cgi->start_Tr();
	$out .= $cgi->td({class=>"warning",rowspan=>$ws->{$g}->{row}},$g);
	foreach my $t (keys %{$ws->{$g}->{tr}}){
		my $rs = scalar (keys %{$ws->{$g}->{tr}->{$t}});
		$out .= $cgi->td({class=>"warning",rowspan=>$rs}, $t);
		foreach my $e (keys %{$ws->{$g}->{tr}->{$t}}){
				$out .= $cgi->td({class=>"warning"}, $e);
				$out .= $cgi->td({class=>"warning"}, $ws->{$g}->{tr}->{$t}->{$e}->{percent});
				$out .= $cgi->td({class=>"warning"}, $ws->{$g}->{tr}->{$t}->{$e}->{pos});
			
					$out .= $cgi->end_Tr();
		}
		$out .= $cgi->end_Tr();
		
	}
	$out .= $cgi->end_Tr();

}

$out .= $cgi->end_table();
 $out .= $cgi->end_div();
  $out .= $cgi->end_div();
	 return ($out1,$out);
	
}
sub table_patients_printer2 {
	 my ($project) = @_;
	 my $out;
		$out.="<body onload='this.print();'>";
		$out.="<body '>";
	$out .= $cgi->start_table({class=>"table table-striped table-condensed table-bordered table-hover ",style=>"vertical-align:middle;text-align: left;font-size: 8px;font-family:  Verdana;"});
		 my $capture = $project->getCaptures()->[0];
 		my $ver1 = $project->genome_version;
 		my $ver = $project->annotation_version;
 		my $machine = $project->getRuns->[0]->machine();
 		
		$out .= $cgi->start_Tr();
		my  $div .=   $cgi->start_div();
		$div.=   $cgi->h3($project->name.":".$project->description);
		$div.=   $cgi->h4(join("-", return_date($project->creation_date)  )."&nbsp;". $machine."&nbsp;".$capture->name);
		 $div .=   $cgi->end_div();
		$out .= $cgi->th({colspan=>13},$div);
		$out .= $cgi->end_Tr();
	my @infos = ("Familly","Patient","Sex","SRY","Status","Cov","30X","Analyse","validation","nb variants","user");
		$out .= $cgi->start_Tr();
		for (my $i=0;$i<@infos;$i++){
			$out .= $cgi->th($infos[$i]);
		}   
	 $out .= $cgi->end_Tr();
	 my $hsex = {
				2=>qq{  <i class="fa fa-venus" > </i>},
				1=>qq{  <i class="fa fa-mars" > </i>} ,
				'-1'=>qq{  <i class="fa fa-minus" > </i>} 
		};
	 foreach my $patient (sort {$a->name cmp $b->name} @{$project->getPatients}){
	 	$out .= $cgi->start_Tr();
	 	$out .= $cgi->td({style=>"border: 1px solid black;"},$patient->getFamily->name);
	 	$out .= $cgi->td({style=>"border: 1px solid black;"},$patient->name);
	 	my $sex_eval = $patient->compute_sex(); 
		my $sex = $patient->sex(); 
		my $class ;
		$class->{style} ="border:1px";
		if ($sex_eval ne $sex && $sex_eval ne -1){
			$class->{class} ="danger";
		} 
		my $cov_sry = $patient->coverage_SRY();
		$out .= $cgi->td({style=>"border: 1px solid black;"},,$hsex->{$sex});
		$out .= $cgi->td({style=>"border: 1px solid black;"},,$hsex->{$sex_eval}."<small>(".$cov_sry.")</small>");
		
			if ($patient->isMother){
				$out .= $cgi->td({style=>"border: 1px solid black;"},,qq{<i class="fa fa-female fa-1x" aria-hidden="true" style="color:pink"></i>});
			}
			elsif ($patient->isFather){
				$out .= $cgi->td({style=>"border: 1px solid black;"},,qq{<i class="fa fa-male fa-1x" aria-hidden="true" style="color:blue"></i>});
			}
			else {
					$out .= $cgi->td({style=>"border: 1px solid black;"},,qq{<i class="fa fa-child fa-1x" aria-hidden="true" style="color:red"></i>});
				
			}
		my $cov = $patient->coverage();
		$out .= $cgi->td({style=>"border: 1px solid black;"},$cov->{mean});
		$out .= $cgi->td({style=>"border: 1px solid black;"},$cov->{"30x"}."%");
	  	  	
  		my $hval = $project->validations_query->getValidationPatient($patient);
  		
  		if (keys %$hval){
  			my @values = values %$hval;
  			
  			my $val = $values[0]->[0];
  			my $st_date = join("-",return_date($val->{modification_date}));
  			$out .= $cgi->td({style=>"border: 1px solid black;"},$st_date);
			 my $term = $buffer->getValidationTerm($val->{validation});
			 $out .= $cgi->td({style=>"border: 1px solid black;"},$term);
			my $pn = $val->{user_name};
			 $out .= $cgi->td({style=>"border: 1px solid black;"},scalar(@values));
			  $out .= $cgi->td({style=>"border: 1px solid black;"},$pn);
  		
  			
  		}
  		
  		else {
  			$out .= $cgi->td({style=>"border: 1px solid black;"},"");
  			$out .= $cgi->td({style=>"border: 1px solid black;"},"");
  			$out .= $cgi->td({style=>"border: 1px solid black;"},"");
  			$out .= $cgi->td({style=>"border: 1px solid black;"},"");
  			
  		}
  	$out .= $cgi->end_Tr();
	  	
	}
	$out .= $cgi->end_table();
	print $out;
	exit(0);
	 
}
sub table_patients_printer {
	 my ($project) = @_;
	 my $out;
		$out.="<body onload='this.print();'>";
		$out.="<body '>";
	$out .= $cgi->start_table({class=>"table table-striped table-condensed table-bordered table-hover ",style=>"vertical-align:middle;text-align: left;font-size: 8px;font-family:  Verdana;"});
		 my $capture = $project->getCaptures()->[0];
 		my $ver1 = $project->genome_version;
 		my $ver = $project->annotation_version;
 		my $machine = $project->getRuns->[0]->machine();
 		
		$out .= $cgi->start_Tr();
		my  $div .=   $cgi->start_div();
		$div.=   $cgi->h3($project->name.":".$project->description);
		$div.=   $cgi->h4(join("-", return_date($project->creation_date)  )."&nbsp;". $machine."&nbsp;".$capture->name);
		 $div .=   $cgi->end_div();
		$out .= $cgi->th({colspan=>13},$div);
		$out .= $cgi->end_Tr();
	my @infos = ("Familly","Patient","Sex","SRY","Status","Cov","30X","Analyse","validation","nb variants","user");
		$out .= $cgi->start_Tr();
		for (my $i=0;$i<@infos;$i++){
			$out .= $cgi->th($infos[$i]);
		}   
	 $out .= $cgi->end_Tr();
	 my $hsex = {
				2=>qq{  <i class="fa fa-venus" > </i>},
				1=>qq{  <i class="fa fa-mars" > </i>} ,
				'-1'=>qq{  <i class="fa fa-minus" > </i>} 
		};
	 foreach my $patient (sort {$a->name cmp $b->name} @{$project->getPatients}){
	 	$out .= $cgi->start_Tr();
	 	$out .= $cgi->td({style=>"border: 1px solid black;"},$patient->getFamily->name);
	 	$out .= $cgi->td({style=>"border: 1px solid black;"},$patient->name);
	 	my $sex_eval = $patient->compute_sex(); 
		my $sex = $patient->sex(); 
		my $class ;
		$class->{style} ="border:1px";
		if ($sex_eval ne $sex && $sex_eval ne -1){
			$class->{class} ="danger";
		} 
		my $cov_sry = $patient->coverage_SRY();
		$out .= $cgi->td({style=>"border: 1px solid black;"},,$hsex->{$sex});
		$out .= $cgi->td({style=>"border: 1px solid black;"},,$hsex->{$sex_eval}."<small>(".$cov_sry.")</small>");
		
			if ($patient->isMother){
				$out .= $cgi->td({style=>"border: 1px solid black;"},,qq{<i class="fa fa-female fa-1x" aria-hidden="true" style="color:pink"></i>});
			}
			elsif ($patient->isFather){
				$out .= $cgi->td({style=>"border: 1px solid black;"},,qq{<i class="fa fa-male fa-1x" aria-hidden="true" style="color:blue"></i>});
			}
			else {
					$out .= $cgi->td({style=>"border: 1px solid black;"},,qq{<i class="fa fa-child fa-1x" aria-hidden="true" style="color:red"></i>});
				
			}
		my $cov = $patient->coverage();
		$out .= $cgi->td({style=>"border: 1px solid black;"},$cov->{mean});
		$out .= $cgi->td({style=>"border: 1px solid black;"},$cov->{"30x"}."%");
	  	  	
  		my $hval = $project->validations_query->getValidationPatient($patient);
  		
  		if (keys %$hval){
  			my @values = values %$hval;
  			
  			my $val = $values[0]->[0];
  			my $st_date = join("-",return_date($val->{modification_date}));
  			$out .= $cgi->td({style=>"border: 1px solid black;"},$st_date);
			 my $term = $buffer->getValidationTerm($val->{validation});
			 $out .= $cgi->td({style=>"border: 1px solid black;"},$term);
			my $pn = $val->{user_name};
			 $out .= $cgi->td({style=>"border: 1px solid black;"},scalar(@values));
			  $out .= $cgi->td({style=>"border: 1px solid black;"},$pn);
  		
  			
  		}
  		
  		else {
  			$out .= $cgi->td({style=>"border: 1px solid black;"},"");
  			$out .= $cgi->td({style=>"border: 1px solid black;"},"");
  			$out .= $cgi->td({style=>"border: 1px solid black;"},"");
  			$out .= $cgi->td({style=>"border: 1px solid black;"},"");
  			
  		}
  	$out .= $cgi->end_Tr();
	  	
	}
	$out .= $cgi->end_table();
	print $out;
	exit(0);
	 
}


my @header_transcripts = ("consequence","enst","nm","ccds","appris","exon","nomenclature","codons","codons_AA", "polyphen","sift","cadd","revel","dbscsnv",'spliceAI');

sub merge_range{
	my ($array,$row,$col,$merge,$value) = @_;
	push(@$array,{row=>$row,col=>$col,merge=>$merge,value=>$value});
}

sub printPatientLineXls {
	my ($patient) = @_;
my $hval = $patient->validations();
return  unless %$hval;
my $line;
	my @header_transcripts = ("consequence","enst","nm","ccds","appris","exon","nomenclature","codons","codons_AA", "polyphen","sift","cadd","revel","dbscsnv",'spliceAI');

my $line_patient = [];

	my $row =1;

 	foreach my $k  (keys %{$hval}){
 		my $col =0;
			 my $val = $hval->{$k}->[0];
			  my $v  = $project->_newVariant($val->{polyid});
			  my $no       = $v->getChromosome->lmdb_polyviewer_variants( $patient, "r" );
			 my $hvariant = $no->get($val->{polyid});
		
			 
			# update_variant_editor::vgnomad($v,$hvariant);
		
			  my $gene ;
			 		if ($val->{gene_id}){
			 			 $gene = $project->newGene($val->{gene_id}); 
			 		}
			 		else {
			 			if ($project->isDiagnostic){
			 				my $l = list_genes($patient->project);
			 				($gene) = grep {exists $l->{$_->id}} @{ $v->getGenes()};
			 			}
			 			else {
			 				die();
			 			}
			 		}
			 		#my @transcripts = @{$gene->getTranscripts}; 
			 		my $htranscripts = $hvariant->{genes}->{$gene->id};
			 		my $nb_merge = scalar(@$htranscripts);
			 		#warn $nb_merge;
			 		my $st_date = join("-",return_date($val->{modification_date}));
			 		my $term = $val->{term};
			 		my $pn = $val->{user_name};
			 		
			 		merge_range($line_patient,$row,$col,$nb_merge,$patient->name);$col++;
			 		merge_range($line_patient,$row,$col,$nb_merge,$term." - $pn - $st_date");$col++;
			 		merge_range($line_patient,$row,$col,$nb_merge,$v->name);$col++;
			 		merge_range($line_patient,$row,$col,$nb_merge,$v->getChromosome()->name.":".$v->start."-".$v->end);$col++;
			 		merge_range($line_patient,$row,$col,$nb_merge,join(";",@{$hvariant->{value}->{ngs}}) );$col++;
			 		##############
					####GNOMAD 
					##############
					update_variant_editor::vgnomad($v,$hvariant);
			 		merge_range($line_patient,$row,$col,$nb_merge,$hvariant->{value}->{ac});$col++;
			 		merge_range($line_patient,$row,$col,$nb_merge,$hvariant->{value}->{ho});$col++;
			 		
			 		merge_range($line_patient,$row,$col,$nb_merge,$hvariant->{value}->{an});$col++;
			 		merge_range($line_patient,$row,$col,$nb_merge,$hvariant->{value}->{max_pop});$col++;
			 		merge_range($line_patient,$row,$col,$nb_merge,$hvariant->{value}->{min_pop});$col++;
			 		##############
					####DEJAVU 
					##############
					merge_range($line_patient,$row,$col,$nb_merge,$hvariant->{value}->{other_project}."::".$hvariant->{value}->{other_patients}."::".$hvariant->{value}->{other_patients_ho});$col++;
					merge_range($line_patient,$row,$col,$nb_merge,$hvariant->{value}->{similar_projects}."::".$hvariant->{value}->{similar_patients}."::".$hvariant->{value}->{similar_patients_ho});$col++;
					merge_range($line_patient,$row,$col,$nb_merge,,$hvariant->{value}->{this_run_patients});$col++;
					
#			 		$worksheet->merge_range($row,$col,$row+$nb_merge,$col,$v->gnomadAC,$bold_merge);$col++;
					 merge_range($line_patient,$row,$col,$nb_merge,$gene->external_name);$col++;
					 ###################
					 # VALIDATION 
					 #####################
					
					 
					 update_variant_editor::check_is_hgmd_dm_for_gene($patient->getProject(), $hvariant, $gene);
					 update_variant_editor::check_is_clinvar_pathogenic_for_gene($patient->getProject(), $hvariant, $gene);
					 $hvariant->{value}->{hgmd} = "-" if $hvariant->{value}->{hgmd} =~ /minus/;
					 $hvariant->{value}->{clinvar} = "-" if $hvariant->{value}->{clinvar} =~ /minus/;
					 merge_range($line_patient,$row,$col,$nb_merge,$hvariant->{value}->{hgmd});$col++;
					 merge_range($line_patient,$row,$col,$nb_merge,$hvariant->{value}->{clinvar});$col++;
			 		
			 		#$worksheet->merge_range($row,$col,$row+$nb_merge,$col,$gene->external_name,$bold_merge);$col++;
			 
			 		my $hvariation = update::construct_variant( $project, $v, $gene->getTranscripts->[0], $patient, undef );
		
			 		foreach my $t (@$htranscripts){
			 			my $vcol = $col;
			 			foreach my $v1 (@header_transcripts){
			 				merge_range($line_patient,$row,$vcol++,0,$t->{value}->{$v1});
			 			}
			 			$row++;
			 			
			 		}
			 		$row ++;
			 	#	$worksheet->write($row,$ind,$v->getChromosome()->name.":".$v->start."-".$v->end);$ind++;
			 	#	$worksheet->write($row,$ind,$gene->external_name);$ind++;
			 		
			 #		$row++;
 	}
	return $line_patient;
		
}

sub printWorkSheet {
		my ($patient,$workbook) = @_;
my $hval = $patient->validations();
return unless %$hval;
my $line_patient;
my $bold_merge       = $workbook->add_format(valign      => 'vcentre',
                                        align       => 'centre',
                                        );
	my $worksheet = $workbook->add_worksheet($patient->name);
	my $row =1;
	my @header_transcripts = ("consequence","enst","nm","ccds","appris","exon","nomenclature","codons","codons_AA", "polyphen","sift","cadd","revel","dbscsnv",'spliceAI');
 	foreach my $k  (keys %{$hval}){
 		my $col =0;
			 my $val = $hval->{$k}->[0];
			  my $v  = $project->_newVariant($val->{polyid});
			  my $no       = $v->getChromosome->lmdb_polyviewer_variants( $patient, "r" );
			 my $hvariant = $no->get($val->{polyid});
		
			 
			# update_variant_editor::vgnomad($v,$hvariant);
		
			  my $gene ;
			 		if ($val->{gene_id}){
			 			 $gene = $project->newGene($val->{gene_id}); 
			 		}
			 		else {
			 			if ($project->isDiagnostic){
			 				my $l = list_genes($patient->project);
			 				($gene) = grep {exists $l->{$_->id}} @{ $v->getGenes()};
			 			}
			 			else {
			 				die();
			 			}
			 		}
			 		#my @transcripts = @{$gene->getTranscripts}; 
			 		my $htranscripts = $hvariant->{genes}->{$gene->id};
			 		my $nb_merge = scalar(@$htranscripts);
			 			 warn $val->{polyid}." ".$row." ".@$htranscripts;
			 		#warn $nb_merge;
			 		$worksheet->merge_range($row,$col,$row+$nb_merge,$col,$patient->name,$bold_merge);$col++;
			 		$worksheet->merge_range($row,$col,$row+$nb_merge,$col,$v->name,$bold_merge);$col++;
			 		$worksheet->merge_range($row,$col,$row+$nb_merge,$col,$v->getChromosome()->name.":".$v->start."-".$v->end,$bold_merge);$col++;
			 		$worksheet->merge_range($row,$col,$row+$nb_merge,$col,join(";",@{$hvariant->{value}->{ngs}}),$bold_merge);$col++;
			 		##############
					####GNOMAD 
					##############
			 		 update_variant_editor::vgnomad($v,$hvariant);
			 		 $worksheet->merge_range($row,$col,$row+$nb_merge,$col,$hvariant->{value}->{ac},$bold_merge);$col++;
					 $worksheet->merge_range($row,$col,$row+$nb_merge,$col,$hvariant->{value}->{an},$bold_merge);$col++;
					 $worksheet->merge_range($row,$col,$row+$nb_merge,$col,$hvariant->{value}->{ac_ho},$bold_merge);$col++; 
					 $worksheet->merge_range($row,$col,$row+$nb_merge,$col,$hvariant->{value}->{max_pop},$bold_merge);$col++; 
					 $worksheet->merge_range($row,$col,$row+$nb_merge,$col,$hvariant->{value}->{min_pop},$bold_merge);$col++; 
					 
					##############
					####DEJAVU 
					##############
					 $worksheet->merge_range($row,$col,$row+$nb_merge,$col,$hvariant->{value}->{other_project}."::".$hvariant->{value}->{other_patients}."::".$hvariant->{value}->{other_patients_ho},$bold_merge);$col++;
					 $worksheet->merge_range($row,$col,$row+$nb_merge,$col,$hvariant->{value}->{similar_projects}."::".$hvariant->{value}->{similar_patients}."::".$hvariant->{value}->{similar_patients_ho},$bold_merge);$col++;
					 $worksheet->merge_range($row,$col,$row+$nb_merge,$col,$hvariant->{value}->{this_run_patients},$bold_merge);$col++;
					 	 
					 
#			 		$worksheet->merge_range($row,$col,$row+$nb_merge,$col,$v->gnomadAC,$bold_merge);$col++;

					 
					 ###################
					 # VALIDATION 
					 #####################
					$worksheet->merge_range($row,$col,$row+$nb_merge,$col,$hvariant->{value}->{hgmd},$bold_merge);$col++;
					$worksheet->merge_range($row,$col,$row+$nb_merge,$col,$hvariant->{value}->{clinvar},$bold_merge);$col++;	 		
			 		
			 		#$worksheet->merge_range($row,$col,$row+$nb_merge,$col,$gene->external_name,$bold_merge);$col++;
			 
			 		my $hvariation = update::construct_variant( $project, $v, $gene->getTranscripts->[0], $patient, undef );
		
			 		foreach my $t (@$htranscripts){
			 			my $vcol = $col;
			 			foreach my $v1 (@header_transcripts){
			 			$worksheet->write($row,$vcol++,$t->{value}->{$v1});
			 			}
			 			$row++;
			 			
			 		}
			 		$row ++;
			 	#	$worksheet->write($row,$ind,$v->getChromosome()->name.":".$v->start."-".$v->end);$ind++;
			 	#	$worksheet->write($row,$ind,$gene->external_name);$ind++;
			 		
			 #		$row++;
 	}
	
}
sub table_patients_xls {
	 my ($project) = @_;
	my $workbook  = Spreadsheet::WriteExcel->new(\*STDOUT);
		my $bold_merge       = $workbook->add_format(valign      => 'vcentre',
                                        align       => 'centre',
                                        );
	my $worksheet = $workbook->add_worksheet($project->name);
	my @header_transcripts = ("consequence","enst","nm","ccds","appris","exon","nomenclature","codons","codons_AA", "polyphen","sift","cadd","revel","dbscsnv",'spliceAI');
		my $row = 0;
	foreach my $patient (@{$project->getPatients}){
		my $line_patient=printPatientLineXls($patient);
		next unless $line_patient;
	
		my $col =0;
		#$worksheet->merge_range($row,1,$row+$line_patient->[-1]->{row}+1,$col,$patient->name,$bold_merge);
		foreach my $l (@$line_patient){
			if ($l->{merge}>0){
			$worksheet->merge_range($l->{row}+$row,$l->{col}+$col,$l->{row}+$l->{merge}+$row,$l->{col}+$col,$l->{value},$bold_merge);
			}
			else {
#				warn $l->{row}." ".$l->{col};
				$worksheet->write($l->{row}+$row,$l->{col}+$col,$l->{value});
			}
		}
		$row+= $line_patient->[-1]->{row}+1;
		#last;
		
	}
	$workbook->close();

		exit(0);
}


sub table_patients_xls2 {
	 my ($project) = @_;
		my $workbook  = Spreadsheet::WriteExcel->new(\*STDOUT);
	my $worksheet = $workbook->add_worksheet();
	my $bg_color;
	my @colors = ("red","orange","blue","green","cyan","gray");

	foreach my $c (@colors){
		$bg_color->{$c} = $workbook->add_format(valign      => 'vcentre',
                                        align       => 'centre',
                                        #bg_color => 'white',
                                        color => $c,
                                        bold=>1,
                     );                                 
	}
	$bg_color->{strike} = $workbook->add_format(valign      => 'vcentre',
                                        align       => 'centre',
                                        #bg_color => 'white',
                                        color => "gray",
                                        bold=>1,
                                        font_strikeout=>1,
                     );    
       my $col = 0; 
		my $row =0;
		my @infos = ("Familly","Patient","Sex","SRY","Status","Cov","30X","Analyse","user","validation");
		
		for (my $i=0;$i<@infos;$i++){
			$worksheet->write($row,$i,$infos[$i]);
		}        
	$row ++;
	
	foreach my $patient (sort {$a->name cmp $b->name} @{$project->getPatients}){
		$col =0;
		$worksheet->write($row,$col++,$patient->getFamily->name);
		$worksheet->write($row,$col++,$patient->name);
		$worksheet->write($row,$col++,$patient->sex);
		my $cov_sry = $patient->coverage_SRY();
	  	my $sex_eval = $patient->compute_sex();
	  	$worksheet->write($row,$col++,$cov_sry); 
	  	$worksheet->write($row,$col++,$patient->status); 
	  	my $cov = $patient->coverage();
	  	$worksheet->write($row,$col++,$cov->{"mean"});
	  	$worksheet->write($row,$col++,$cov->{"30x"});
  		my $hval = $project->validations_query->getValidationPatient($patient);
  		
  		if (keys %$hval){
  			my @values = values %$hval;
  			
  			my $val = $values[0]->[0];
  			my $st_date = join("-",return_date($val->{modification_date}));
			 my $term = $buffer->getValidationTerm($val->{validation});
			my $pn = $val->{user_name};
  			$worksheet->write($row,$col++,$st_date);
  			$worksheet->write($row,$col++,scalar(@values));
  			$worksheet->write($row,$col++,$term);
  			$worksheet->write($row,$col++,$pn);
  			
  		}
  	$row ++;
	  	
	}
	$workbook->close();
	exit(0);
}
sub table_patients {
		 my ($run) = @_;
		 my $out;
		 $out .= $cgi->start_table({class=>"table  table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 9px;font-family:  Verdana;padding:10px", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
		
		#### print header table
		 my $col_hgmd =3;
	 	$col_hgmd = 2 unless $hgmd ==1;

		my $html_select_all = qq{<input id="check_all" style="border:solid 2px white;margin:2px;" type="checkbox" aria-label="..."  onchange="select_all_patients(this)" checked></input>};
		
		#Nb Junctions
		my $hasJunctions;
		foreach my $p (@{$project->getPatients()}) { $hasJunctions = 1 if ($p->hasJunctions()); }
		my (@title, @lcol);
		if ($project->isDiagnostic){
			if ($hasJunctions) {
				@title = ("Fam",$html_select_all,"Patient","View","Print","Cov","30x","Nb Junctions",'PolyViewer','Variations Editor','PolyDupDel','PolyCyto','PolySplice');
				@lcol = ("Fam","Select","Patient","View","Print","Cov","30x","Nb Junctions",'PolyViewer','Variations Editor','PolyDupDel','PolyCyto','PolySplice');
			}
			else {
				@title = ("Fam",$html_select_all,"Patient","View","Print","Cov","30x",'PolyViewer','Variations Editor','PolyDupDel','PolyCyto','PolySplice');
				@lcol = ("Fam","Select","Patient","View","Print","Cov","30x",'PolyViewer','Variations Editor','PolyDupDel','PolyCyto','PolySplice');
			}
		}
		else {
			if ($hasJunctions) {
				@title = ("Fam","Patient","View","Print","Cov","30x","Nb Junctions",'PolyViewer','Variations Editor','PolyDupDel','PolyCyto','PolySplice');
				@lcol = ("Fam","Patient","View","Print","Cov","30x","Nb Junctions",'PolyViewer','Variations Editor','PolyDupDel','PolyCyto','PolySplice');
			}
			else {
				@title = ("Fam","Patient","View","Print","Cov","30x",'PolyViewer','Variations Editor','PolyDupDel','PolyCyto','PolySplice');
				@lcol = ("Fam","Patient","View","Print","Cov","30x",'PolyViewer','Variations Editor','PolyDupDel','PolyCyto','PolySplice');
			}
		}
		my $isfam;
		$isfam = 1; 
		
		push(@title,"MUC1") if (-e  $project->getVariationsDir("vntyper")."/muc1/" );
		$out .= $cgi->start_Tr({style=>"background-color:#1079B2;color:white"});
		foreach my $p (@title){
			$out .= $cgi->th({style=>"text-align: center;min-width:5%"},$p);
		}
		push(@title,"validation");
		$out .= $cgi->th({style=>"text-align: center;min-width:70%"},"validation");
		$out .= $cgi->end_Tr();
		
		#end header
		
		my @colors = ("#F4F4F4","#DEDFDE");
		@colors = ("#F9F6FF","#FFFFFF");
		#F0EAD6
		my $nbf =0;
		my $nb = 0;
		
		foreach my $fam (sort{$a->name cmp $b->name }@{$run->getFamilies}){
		my $color = $colors[$nbf%2]; 
		$nbf ++;
		$nb ++;
  		my $nb_members = scalar(@{$fam->getMembers});
  		foreach my $p1 (@{$fam->getMembers}){
  			my $hval = $project->validations_query->getValidationPatient($p1);
  			
  		}
  		
  		$out.= $cgi->start_Tr({style =>"background-color:$color"});
  		my $pname = "check_".$fam->name();
  			
  	#	 if ($nb_members>1){

  		 	
  	#	 $out .= $cgi->td({style=>"vertical-align:middle"},qq{<input id="$pname" type="checkbox" aria-label="..." onClick="selection(event,this)"></input>});



  		 $out .= $cgi->td({rowspan=>$nb_members,style=>"vertical-align:middle"},$fam->name );# if $isfam ;
			
		foreach my $p (@{$fam->getMembers}){
			my $pp = $p->name;
			my $cmd = qq{printer('$pp');};
			my $cmd2 = qq{printer2('$pp','1');};
			my $line = print_line_patient($p,0);
			$line->{Select} = $cgi->td({style=>"vertical-align:middle"},qq{<input id="$pname" type="checkbox" aria-label="..." onClick="select_patient(this)"  checked></input>});
			$line->{View} = $cgi->td({style=>"vertical-align:middle"},'<a type="button" class="btn btn-xs btn-info" onclick="'.$cmd.'"><i class="fa fa-clipboard pull-left  "></i>View</button></a>');
			$line->{Print} = $cgi->td({style=>"vertical-align:middle"},'<a type="button" class="btn btn-xs btn-success" onclick="'.$cmd2.'"><i class="fa fa-print pull-left  "></i>Print</button></a>');
			foreach my $col (@lcol){
				if ($line->{$col}) { $out .= $line->{$col}; }
				else { $out .= '-'; }
			}
			$out .= $cgi->end_Tr;
			$out .= $cgi->start_Tr({style =>"background-color:$color"}) if $nb_members >1;
		}
		$out.= $cgi->end_Tr  if $nb_members >1;
	}
	$out .= $cgi->end_table();
$out .= $cgi->end_div();
		 return $out;
		 
		 
}

sub table_run_header {
	 my ($run) = @_;
	 my $out ;
	 $out .= header_run($run);
	  $out .=  $cgi->start_table({style=>"margin-bottom: 5px;margin-right: 10px;margin-left: 100px;margin-top: 5px;"});
		my ($b,$p) =   table_sex($run);
			my ($b1,$p1) =   table_duplicate($run);
			my ($b2,$p2) =   table_control($run);
			my ($b3,$p3) =  table_mendelian($run);
			my ($b4,$p4) =   table_quality($run);
			my ($b5,$p5) =   table_design($run);
			my ($b6,$p6) = table_muc1($run);
				$out.=  $cgi->start_Tr;
		 $out .= $cgi->th({style=>"margin-bottom: 5px;margin-right: 10px;margin-left: 100px;margin-top: 5px;"},$b);
		  $out .= $cgi->th({style=>"margin-bottom: 5px;margin-right: 10px;margin-left: 100px;margin-top: 5px;"},$b4);
		  $out .= $cgi->th({style=>"margin-bottom: 5px;margin-right: 10px;margin-left: 100px;margin-top: 5px;"},$b3);
		  $out .= $cgi->th({style=>"margin-bottom: 5px;margin-right: 10px;margin-left: 100px;margin-top: 5px;"},$b2);
		   $out .= $cgi->th({style=>"margin-bottom: 5px;margin-right: 10px;margin-left: 100px;margin-top: 5px;"},$b1);
		    $out .= $cgi->th({style=>"margin-bottom: 5px;margin-right: 10px;margin-left: 100px;margin-top: 5px;"},$b5);
		       $out .= $cgi->th({style=>"margin-bottom: 5px;margin-right: 10px;margin-left: 100px;margin-top: 5px;"},$b6);
		   $out .= $cgi->end_Tr;
			 $out .= $cgi->end_table();
			 
	 $out .= $p1;
	$out .= $p;
	$out .= $p2;
	$out .= $p3;
	$out .= $p4;
	$out .= $p5;
	$out .= $p6;
	return $out;
}
	
  
 my $list_genes;
sub list_genes {
	my ($project) = @_;
	map {$list_genes->{$_->getGene->id}++}  @{$project->getListGenBoTranscripts};
	return $list_genes;
	
}
sub validation_table_cnv {
	my ($p,$hval) = @_;
	
	my @headers_validations = ("type","chromosome","start","end");
	my $rowspan = scalar(keys %$hval);
	my $class;
	$class->{rowspan} -= $rowspan;
	$class->{rowspan} = 1 if $class->{rowspan} <=0;
	$class->{style} = "min-width:10%;padding:1px";
	my $fam = $p->getFamily();
	my $fin = scalar (keys %{$hval}) ;
	my $pos =0;
	my $out;
	#$out.= $cgi->start_table({class=>"table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 9px;font-family:  Verdana;margin-bottom: 0px;padding:2px", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
	#$out.= $cgi->th({colspan=>7},"CNV");
	$out .= update_variant_editor::print_validated_cnv($p,$hval);
	return $out;
		foreach my $k  (keys %{$hval}){
			$out.= update_variant_editor::print_cnv($p, undef,$hval->{$k});
			last;
			 	$out.= $cgi->start_Tr();
			 		my $val = $hval->{$k}->[0];
			 		my $color = "#9BC8A5";
			 		$color = "#E74C3C" if $val->{validation}== 5 ;
			 		$color = "coral" if $val->{validation}== 4 ;
			 		$color = "#217DBB" if $val->{validation}== -3 ;
			 		$color = "orange" if $val->{validation}== 3 ;
			 		my $btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: $color;$fsize"} ;
			 		my $st_date = join("-",return_date($val->{modification_date}));
			 		my $term = $buffer->getValidationTerm($val->{validation});
			 		
			 		$out.= $cgi->td($class,qq{<button type="button" $btn_class >$term</button>});
					 my $pn = $val->{user_name};
			 		 $out.= $cgi->td($class,qq{<button type="button" $btn_class >$pn</button>});
			 		 $out.=  $cgi->td($class,qq{<button type="button" $btn_class >$st_date</button>});
			 		 my ($type,$chr,$start,$end) = split("_",$k);
			 		 $out.=  $cgi->td($class,qq{<button type="button" $btn_class >$type</button>});
			 		  $out.=  $cgi->td($class,qq{<button type="button"  $btn_class >$chr</button>});
			 		   $out.=  $cgi->td($class,qq{<button type="button" $btn_class >$start</button>});
			 		    $out.=  $cgi->td($class,qq{<button type="button" $btn_class >$end</button>});
			 		    
			$out.=  $cgi->end_Tr();
			 		#$pos++;
			 		#print $cgi->start_Tr() if $pos < ($fin-3); 
			 		
			 	}
			 	$out.= $cgi->end_table();
			 	warn $out;
			 	return $out;
}



sub validation_table_new {
		my ($p,$hval) = @_;
			print "*";
				my @headers_validations = ("igv_web","var_name","trio","gene","table_transcript",);
				my @header_transcripts = ("consequence","enst","nm","nomenclature");
			
				my $rowspan = scalar(keys %$hval);
				my $class;
			 	$class->{rowspan} -= $rowspan;
			 	$class->{rowspan} = 1 if $class->{rowspan} <=0;
			 	$class->{style} = "min-width:10%;padding:1px";
			 	my $fam = $p->getFamily();
			 	my $fin = scalar (keys %{$hval}) ;
			 	my $pos =0;
			 	my $out;
			 	$out.= $cgi->start_table({class=>"table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 9px;font-family:  Verdana;margin-bottom: 0px;padding:2px", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
			 	foreach my $k  (keys %{$hval}){
			 		
			 
			 		my $val = $hval->{$k}->[0];
			 			my $color = "#9BC8A5";
			 		$color = "#E74C3C" if $val->{validation}== 5 ;
			 		$color = "coral" if $val->{validation}== 4 ;
			 		$color = "#217DBB" if $val->{validation}== -3 ;
			 		$color = "orange" if $val->{validation}== 3 ;
			 		 	my $v  = $project->_newVariant($val->{polyid});
			 		my $btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: $color;$fsize"} ;
			 	
			 		my $gene ;
			 		if ($val->{gene_id}){
			 			 $gene = $project->newGene($val->{gene_id}); 
			 		}
			 		else {
			 			if ($project->isDiagnostic){
			 				my $l = list_genes($p->project);
			 				($gene) = grep {exists $l->{$_->id}} @{ $v->getGenes()};
			 			}
			 			else {
			 				die();
			 			}
			 		}
			 		my $hvariation = update::construct_variant( $project, $v, $gene->getTranscripts->[0], $p, undef );
					$hvariation->{'igv_web'} =~ s/></>Igv</g;
			 		my $gn = $gene->external_name;
			 		 $hvariation->{gene} = qq{<button type="button" $btn_class onClick="update_grid_gene_phenotypes('ENSG00000110367_11')">$gn</button>};
			 		update::trio($project,undef,$hvariation,$p,$cgi,undef);
			 		$hvariation->{table_transcript} = update::construct_table_transcript($v, $cgi,\@header_transcripts,3,$gene,1); 
		
					my $st_date = join("-",return_date($val->{modification_date}));
			 		my $term = $val->{term};
			 		my $pn = $val->{user_name};
			 		$out.=  $cgi->td({style=>"vertical-align:middle;padding:1px;text-align: center;vertical-align:middle"},qq{ <span  class="stamp2" style="border-color:$color;color:$color">+$term+<br><span style='font-size=08px'>$pn:$st_date</span></span>});
			 		 
					foreach my $h (@headers_validations){
							if ($h eq "trio" or "table_transcript"){
								$class->{style} = ";vertical-align:middle;padding:1px";
							}
							else {
									$class->{style} = ";vertical-align:middle;padding:1px";
							}
							
							$out.=  $cgi->td($class,$hvariation->{$h});
					}
					
			 		 
			 		$out.=  $cgi->end_Tr();
			 		#$pos++;
			 		#print $cgi->start_Tr() if $pos < ($fin-3); 
			 		
			 	}
			 	$out.= $cgi->end_table();
			 	return $out;
}

sub args_quality {
	my ($project) = @_;
	my @z;
	if ($project->existsnoSqlQuality) {
		my $t = stat($project->quality_dir."/".$project->name.".lite")->[9];
		push(@z,$t);
	}
#	else {confess()};
	foreach my $p1 (sort{$a->name cmp $b->name} @{$project->getPatients}){
  			my $string = $p1->get_string_identification;
  			push(@z,$string);
  		}
	push(@z,"$VERSION_SCRIPT");#file_md5_hex($Bin."/summary_panel.pl"));	
	#push(@z,"11");	
  	return \@z;	
}	

sub args_validation {
	my ($project) = @_;
	my @z;
	foreach my $p1 (sort{$a->name cmp $b->name} @{$project->getPatients}){
  			my $stv =$p1->get_string_validations();
			next unless $stv;
			#warn $stv;
  			push(@z,$stv);
  			
  		}
  	push(@z,file_md5_hex($Bin."/summary_panel.pl"));	
  	return \@z;	
}	
