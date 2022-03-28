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
use List::MoreUtils qw/pairwise/;

my $nb_gene_by_patient = 3;
my $nb_exon_by_genes = 10;

my $buffer = GBuffer->new();
;

my $cgi          = new CGI();
my $project_name = $cgi->param('project');

my $cgi          = new CGI();
#<head>
#<meta charset="utf-8">
#<meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
print $cgi->header(-meta=>{'name'=>'viewport',
				    'content'=>'width=device-width, initial-scale=1, shrink-to-fit=no'},
				    	-meta=>{charset=>"utf-8"});
print qq{
<style>

.main-section{
	width:100%;
	margin:10 auto;
	text-align: center;
	padding: 10px 5px;
}
.dashbord{
	width:20%;
	display: inline-block;
	background-color:#34495E;
	color:#fff;
	margin-top: 50px; 
}
.icon-section i{
	font-size: 30px;
	padding:10px;
	border:1px solid #fff;
	border-radius:50%;
	margin-top:-25px;
	margin-bottom: 10px;
	background-color:#34495E;
}
.icon-section p{
	margin:0px;
	font-size: 20px;
	padding-bottom: 10px;
}
.detail-section{
	background-color: #2F4254;
	padding: 5px 0px;
}
.dashbord .detail-section:hover{
	background-color: #5a5a5a;
	cursor: pointer;
}
.detail-section a{
	color:#fff;
	text-decoration: none;
}
.dashbord-success .icon-section,.dashbord-success .icon-section i{
	background-color: #16A085;
}
.dashbord-success .detail-section{
	background-color: #149077;
}
.dashbord-warning .icon-section,.dashbord-warning .icon-section i{
	background-color: #F39C12;
}
.dashbord-warning .detail-section{
	background-color: #DA8C10;
}
.dashbord-default .icon-section,.dashbord-default .icon-section i{
	background-color: #2980B9;
}
.dashbord-default .detail-section{
	background-color:#2573A6;
}
.dashbord-danger .icon-section,.dashbord-danger .icon-section i{
	background-color:#E74C3C;
}
.dashbord-danger .detail-section{
	background-color:#CF4436;
}
.dashbord-primary .icon-section,.dashbord-primary .icon-section i{
	background-color:#8E44AD;
}
.dashbord-primary .detail-section{
	background-color:#803D9B;
}
.dashbord-grey .icon-section,.dashbord-grey .icon-section i{
	background-color:#CCCCCC;
}
.dashbord-grey .detail-section{
	background-color:#AAAAAA;
}

</style>


<body>
<!-- Basic Card with Image, Title and Description -->
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">
<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.11.0/umd/popper.min.js" integrity="sha384-b/U6ypiBEHpOf/4+1nzFpr53nxSS+GLCkfwBdFNTxtclqqenISfwAzpKaMNFNmj4" crossorigin="anonymous"></script>
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0-beta/js/bootstrap.min.js" integrity="sha384-h0AbiXch4ZDo7tp9hKZ4TsHbi047NrKGLO3SEJAg45jXxnGIfYzk4Si90RDIqNm1" crossorigin="anonymous"></script>
<link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css">
 <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js"></script>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js"></script>
<link href="https://fonts.googleapis.com/icon?family=Material+Icons" rel="stylesheet">

<!-- Latest compiled and minified CSS -->

<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap-table/1.11.1/bootstrap-table.min.css">

<!-- Latest compiled and minified JavaScript -->
<script src="https://cdnjs.cloudflare.com/ajax/libs/bootstrap-table/1.11.1/bootstrap-table.min.js"></script>


<style>
.modal-body {
    max-width: 100%;
    overflow-x: auto;
     overflow-y: auto;
}
</style>

<body>
};

print $cgi->self_url;
my $project;
unless ($project_name){
	my $res = $buffer->listProjects;
	my $resume;
	push(@{$resume->{columns}},{field=>"project",title=>"project"});
#	$h->{columns}= {field=>"projects",title=>"projects"};
	
	foreach my $p (@$res){
		my $a = $cgi->a({href=>$cgi->url."?project=".$p},$p);
		push (@{$resume->{data}},{project=>{text=>$a,type=>"default"}});
	}
	print_html_table($resume);
	exit(0);
}
else {
 $project = $buffer->newProject(-name=>$project_name);
}


my @types = ("files","mendelian","identity","duplicate_regions","coverage_stats","bam_stats","coverage_transcripts","statistics_variations");
my $no = $project->noSqlQuality("r");
my $icons = {
	"statistics_variations" => qq{<i class="fa fa-signal" aria-hidden="true"></i>},
	"variants_transcripts" => qq{<i class="fa fa-table" aria-hidden="true"></i>},
	"identity" => qq{<i class="fa fa-user-plus" aria-hidden="true"></i>},
	"files" => qq{<i class="fa fa-file" aria-hidden="true"></i>},
	"mendelian" => qq{<i class="fa fa-users" aria-hidden="true"></i>},
	"duplicate_regions" => qq{<i class="fa fa-copy" aria-hidden="true"></i>},
	"coverage_transcripts" => qq{<i class="fa fa-line-chart" aria-hidden="true"></i>},
	"bam_stats" => qq{<i class="fa fa-line-chart" aria-hidden="true"></i>},
	"coverage_stats" => qq{<i class="fa fa-line-chart" aria-hidden="true"></i>},
};
my $pname = $project->name();
print qq{
		<div class="main-section">
};
foreach my $type (@types){
	my $data = $no->get($project->name,$type);
		my $h = "#".$type;
	my $icon = $icons->{$type};
	unless ($data){
		print qq{
	<div class="dashbord dashbord-grey">
			<div class="icon-section">
				$icon<br>
				<small>$type</small>
				<p>Not Compute</p>
			</div>
			<div class="detail-section" data-toggle="modal" data-target="$h" >
			<span >More Info </span>
			</div>
		</div>
	};
	next;
	}
my $level =  check_level($data);
#die(); 
	print qq{
	<div class="dashbord dashbord-$level">
			<div class="icon-section">
				$icon<br>
				<small>$pname</small>
				<p>$type</p>
			</div>
			<div class="detail-section" data-toggle="modal" data-target="$h" >
			<span >More Info </span>
			</div>
		</div>
	};

	print_html_table2($type,$data);
	
}
print qq{</div>};
exit(0);




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

sub print_table_header {
 	my ($header) = @_;
 	my $text =    $cgi->start_thead({class=>"thead-inverse" });
     $text .=$cgi->start_Tr();
 #     $text .= $cgi->start_thead({class=>"thead-inverse"});
     foreach my $col (@$header){
     	$text .= $cgi->th({class=>"col-md-2;padding:1px",},$col->{title});
     	
     }
     $text .= $cgi->end_thead();
     $text .= $cgi->end_Tr();
 }
 
 sub print_td {
 	my ($l) = @_;
 	my $b = $l->{text};
 	$l->{text} ="-" unless $l->{text};
 	if ($l->{type} !~/default/){
 		#my $b = qq{<button class="btn btn-}.$l->{type}.qq{">}.$l->{text}."</button>";
 		my $b = qq{<div style="padding: 0px;border:0 px;margin:1px" class="alert bg-}.$l->{type}.qq{">}.$l->{text}."</div>";
 		 return  $cgi->td({rowspan=> $l->{rowspan},style=>"vertical-align: middle;padding:1px"},$b) if exists $l->{rowspan};
 		 return  $cgi->td({style=>"vertical-align: middle;padding:1px"},$b) ;
 		
 	}
 	
 	return $cgi->td({class=>"bg-".$l->{type},style=>"vertical-align: middle;padding:1px",rowspan=>$l->{rowspan}},$l->{text})  if $l->{rowspan};
 	return  $cgi->td({class=>"bg-".$l->{type},style=>"vertical-align: middle;padding:px"},$l->{text});
 	
 }
 
  sub print_html_table2 {
	my ($name,$h) = @_;
	my $hn = "#".$name;
	print $cgi->start_div({class=>"modal fade",id=>"$name",role=>"dialog"});
	print $cgi->start_div({class=>"modal-dialog modal-lg",style=>" width:98%; max-width:98%;"});
	print $cgi->start_div({class=>"modal-content"});
	print qq{
		 <div class="modal-header">
        <h4 class="modal-title">$name</h4>
        <button type="button" class="close" data-dismiss="modal">&times;</button>
      </div>
	};
	
	print $cgi->start_div({class=>"modal-body"});

	print $cgi->start_div({class=>"row"});
	
	print $cgi->start_table({class=>"table table-sm table-striped table-bordered table-hover ",style=>"text-align: center;vertical-align:middle;font-size:10px;font-family:  Verdana;padding:3px"});
	
	my $columns = $h->{columns};
	print print_table_header($columns);
	my $ids ;
	map{push(@$ids,$_->{field})} @$columns;
	my $data = $h->{data};
	my $rowspan_id;
	foreach my $line (@$data){
		
		print $cgi->start_Tr({style=>"padding:1px"});
		
		foreach my $id (@$ids){
				my $l = $line->{$id};
				if (exists  $l->{rowspan}){
					my $b = $cgi->button({class=>"btn btn-".$l->{type}},$l->{text});
					print print_td($l) unless exists $rowspan_id->{$id};
				
					 $rowspan_id->{$id}->{nb} ++;
					  $rowspan_id->{$id}->{limit} = $l->{rowspan};
					next;
				}
						print print_td($l)
				#print $cgi->td({class=>"bg-".$l->{type}},$l->{text});
			
				
		}
		
	print $cgi->end_Tr() ;
		foreach my $k (keys %$rowspan_id){
					delete $rowspan_id->{$k} if $rowspan_id->{$k}->{nb} >= $rowspan_id->{$k}->{limit};
		}
	
		
	}
	print $cgi->end_table();
	
	 print $cgi->end_div();
	 print $cgi->end_div(); 
       print $cgi->end_div();
         print $cgi->end_div(); 
           print $cgi->end_div();  
  }
  sub print_html_table{
  	my ($h) =@_;
  	print qq{
  		<div class="container">
  <div class="row">
    <div class="col-sm">
    <div class="col-2">
    </div>
   <div class="col-6">
  	};
  	
  	print $cgi->start_table({class=>"table table-sm table-striped table-bordered table-hover ",style=>"text-align: center;vertical-align:middle;font-size:10px;font-family:  Verdana;padding:3px"});
	
	my $columns = $h->{columns};
	print print_table_header($columns);
	my $ids ;
	map{push(@$ids,$_->{field})} @$columns;
	my $data = $h->{data};
	my $rowspan_id;
	foreach my $line (reverse @$data){
		
		print $cgi->start_Tr({style=>"padding:1px"});
		
		foreach my $id (@$ids){
				my $l = $line->{$id};
				if (exists  $l->{rowspan}){
					my $b = $cgi->button({class=>"btn btn-".$l->{type}},$l->{text});
					print print_td($l) unless exists $rowspan_id->{$id};
				
					 $rowspan_id->{$id}->{nb} ++;
					  $rowspan_id->{$id}->{limit} = $l->{rowspan};
					next;
				}
						print print_td($l)
				#print $cgi->td({class=>"bg-".$l->{type}},$l->{text});
			
				
		}
		
	print $cgi->end_Tr() ;
		foreach my $k (keys %$rowspan_id){
					delete $rowspan_id->{$k} if $rowspan_id->{$k}->{nb} >= $rowspan_id->{$k}->{limit};
		}
	
		
	}
	print $cgi->end_table();
	print qq{
		</div>
		</div>
		</div>
		</div>
	}
  }
  