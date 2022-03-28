package html;
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;


sub print_cadre2 {
	my ($cgi,$title) =  @_;
	my $out2 = $cgi->start_fieldset({class=>"box"});
	$out2.= $cgi->legend($title);
	return $out2;
	
}

sub print_cadre {
	my ($cgi,$title) =  @_;
	my $out2= $cgi->start_div({class=>"panel panel-primary",style=>"font-size: 12px;font-family:  Verdana;"});
	$out2 .= $cgi->div({class=>"panel-heading",syle=>"padding=100px;"},$title)."";
	
	
	return $out2;
	
}

sub end_cadre {
	my $cgi = shift;
	return "</div><br>";
}

sub end_cadre2 {
	my $cgi = shift;
		return $cgi->end_fieldset({class=>"box"})."<div><br></div>";
}

sub print_header_polydiag_printer {
	my ($cgi,$title) = @_;
	print $cgi -> header;
	$title = "PolyDiag" unless $title; 
	print '<link rel="stylesheet" type="text/css" media="print" href=""/polyweb/css/bootstrap.css",">';		
	print '<link rel="stylesheet" type="text/css" media="print" href="/polyweb/css/polydiag/polydiagPrint.css">';	
		
		print $cgi -> start_html(-title=>$title,
						-style=>{
							-src=>['/polyweb/css/polydiag/polydiagPrint.css',
										"/polyweb/css/bootstrap.css",
										"//maxcdn.bootstrapcdn.com/font-awesome/4.3.0/css/font-awesome.min.css"
										],
							},
                             -script=>{
                             	#-src=>  "//code.jquery.com/jquery-1.11.2.min.js",
                             },
                             );
   # print qq{<script type="text/javascript" src="//maxcdn.bootstrapcdn.com/bootstrap/3.3.2/js/bootstrap.min.js""></script> };
                            
	print qq{		
		<style>
		
		  table {   table-layout: fixed;
  }
    
 table th, table td { overflow: hidden; width:50px;}
 
		</style>					
	<style media="print" type="text/css">
  table {   table-layout: fixed;
    word-wrap: break-word;
    overflow: hidden;
    border: 3px solid yellow;

    }
    
 table th, table td { overflow: hidden; width:1em;}
 
 .table-nowrap {
	table-layout:fixed;
}
.table-condensed > thead > tr > th, .table-condensed > tbody > tr > th, .table-condensed > tfoot > tr > th, .table-condensed > thead > tr > td, .table-condensed > tbody > tr > td, .table-condensed > tfoot > tr > td {
    padding: 2px;
 
}
/*
	.table-mybordered > thead > tr > th, .table-mybordered > thead > tr > th, .table-mybordered > tbody > tr > th, .table-mybordered > tbody > tr > th, table-mybordered > tfoot > tr > th, .table-mybordered > tfoot > tr > th, table-mybordered > thead > tr > td, .table-mybordered > thead > tr > td, table-mybordered > tbody > tr > td, .table-mybordered > tbody > tr > td, table-mybordered > tfoot > tr > td, .table-mybordered > tfoot > tr > td {

   padding:2px;
   font-family: Verdana, "Times New Roman", serif; 
   font-size :0.7em;
  
   
}
*/
.table-nowrap td {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;

}
</style>
	};

#		<style>
#	.table tbody  {
#    	font-size: 5px !important;
#    	padding:20px;
#    }
#        .table tbody tr {
#        	padding:50px;
#    }
#        .table tbody th {
#        	padding:50px;
#    }
#      .table tbody td {
#        	padding:50px;
#    }
#		.yellow  {
#    	background-color: yellow !important;
#    	color:black;
#    	   -webkit-print-color-adjust:exact;
#
#    	    	padding:20px;
#    }
#    	.green  {
#    	background-color: #00FF00 !important;
#    	color:black;
#    	   -webkit-print-color-adjust:exact;
#    }
#	.blue  {
#    	background-color: blue !important;
#    	    	color:white;
#    	
#    	color:white;
#    	   -webkit-print-color-adjust:exact;
#    }
#    .cyan  {
#    	background-color: cyan !important;
#    	
#    	color:white;
#    	   -webkit-print-color-adjust:exact;
#    }
#    </style>
#};	
#	<style media="print" type="text/css">
#		body {-webkit-print-color-adjust: exact  !important;}
#	
#   .panel-primary .panel-heading {
#      background-color: #DFF0D7 !important;
#      -webkit-print-color-adjust:exact;
#    }
#    
# 
#    
#    .table tbody  {
#    	font-size: 5px !important;
#    }
#       .default th{
#      background-color: #ECECEC !important;
#      -webkit-print-color-adjust:exact;
#    }
#    
#    .table tbody tr {
#    	color:black;
#    	 	padding:1px;
#    }
#    .table tbody tr > th {
#    	 	padding:1px;
#    }
#.table tbody tr > td.danger {
#  background-color: #E74C3C !important;
#    	   -webkit-print-color-adjust:exact;
#    	  padding:1px;
#    }
#    
#    .table tbody tr > td.yellow {
#
#  background-color: yellow !important;
#      -webkit-print-color-adjust:exact;
#       	padding:1px;
#    }
#   .table tbody tr > td.green {
#  		background-color:  #00FF00 !important;
#      -webkit-print-color-adjust:exact;
#       padding:1px;
#    }
#    
#</style>
	
}

sub print_header_polydiag {
	my ($cgi,$title) = @_;
	print $cgi -> header;
	$title = "PolyDiag " unless $title; 
		print $cgi -> start_html(-title=>$title,
						-style=>{
						
							-src=>['/polyweb/css/polydiag/polydiag.css',
										"//maxcdn.bootstrapcdn.com/bootswatch/3.3.4/flatly/bootstrap.min.css",
										"//maxcdn.bootstrapcdn.com/font-awesome/4.3.0/css/font-awesome.min.css",
										#"//cdnjs.cloudflare.com/ajax/libs/bootstrap-table/1.7.0/bootstrap-table.min.css"
										
										],
							},
							
							 -script=>{
                             	-src=>  "//code.jquery.com/jquery-1.11.2.min.js",
                             },
                             );
    print qq{<script type="text/javascript" src="//maxcdn.bootstrapcdn.com/bootstrap/3.3.2/js/bootstrap.min.js""></script> };
    print qq{
    	<link rel="stylesheet" href="//cdnjs.cloudflare.com/ajax/libs/bootstrap-table/1.7.0/bootstrap-table.min.css">
    	<script src="//cdnjs.cloudflare.com/ajax/libs/bootstrap-table/1.7.0/bootstrap-table.min.js"></script>
    	<script src="//cdnjs.cloudflare.com/ajax/libs/bootstrap-table/1.7.0/locale/bootstrap-table-zh-CN.min.js"></script>
    };               
					
print qq{		
		<style>
		table-mybordered > thead > tr > th, .table-mybordered > thead > tr > th, table-mybordered > tbody > tr > th, .table-mybordered > tbody > tr > th, table-mybordered > tfoot > tr > th, .table-mybordered > tfoot > tr > th, table-mybordered > thead > tr > td, .table-mybordered > thead > tr > td, table-mybordered > tbody > tr > td, .table-mybordered > tbody > tr > td, table-mybordered > tfoot > tr > td, .table-mybordered > tfoot > tr > td {
  border: 1px solid #95A5A6;
  vertical-align:middle;
   text-align:center; 
   padding:2px;
}


.bs-callout {
    padding: 20px;
    margin: 20px 0;
    border: 1px solid #eee;
    border-left-width: 5px;
    border-radius: 3px;
}
.bs-callout h4 {
    margin-top: 0;
    margin-bottom: 5px;
}
.bs-callout p:last-child {
    margin-bottom: 0;
}
.bs-callout code {
    border-radius: 3px;
}
.bs-callout+.bs-callout {
    margin-top: -5px;
}
.bs-callout-default {
    border-left-color: #777;
}
.bs-callout-default h4 {
    color: #777;
}
.bs-callout-primary {
    border-left-color: #428bca;
}
.bs-callout-primary h4 {
    color: #428bca;
}
.bs-callout-success {
    border-left-color: #5cb85c;
}
.bs-callout-success h4 {
    color: #5cb85c;
}
.bs-callout-danger {
    border-left-color: #d9534f;
}
.bs-callout-danger h4 {
    color: #d9534f;
}
.bs-callout-warning {
    border-left-color: #f0ad4e;
}
.bs-callout-warning h4 {
    color: #f0ad4e;
}
.bs-callout-info {
    border-left-color: #5bc0de;
}
.bs-callout-info h4 {
    color: #5bc0de;
}
.badge1 {
  padding: 1px 9px 2px;
  font-size: 12.025px;
  font-weight: bold;
  white-space: nowrap;
  color: #ffffff;
  background-color: #999999;
  -webkit-border-radius: 9px;
  -moz-border-radius: 9px;
  border-radius: 9px;
}
		</style>	
<script type="text/javascript">


	function zoomCoverage(tr){
	  	var start_url = "/cgi-bin/polymorphism-cgi/polydiag_old/coverage_exons.pl";
	 	var url = construct_url(start_url,tr);
	 	dijit.byId("popup_coverage").show();
	 	dijit.byId("div_popup_coverage").setHref(url);
	}


	function zoomCoverage(tr){
	  	var start_url = "/cgi-bin/polymorphism-cgi/polydiag_old/coverage_exons.pl";
	 	var url = construct_url(start_url,tr);
	 	dijit.byId("popup_coverage").show();
	 	dijit.byId("div_popup_coverage").setHref(url);
	}

    jQuery(function (\$) {
        
        \$('.panel-heading span.clickable').on("click", function (e) {
            if (\$(this).hasClass('panel-collapsed')) {
                // expand the panel
                \$(this).parents('.panel').find('.panel-body').slideDown();
                \$(this).removeClass('panel-collapsed');
                \$(this).find('i').removeClass('glyphicon-chevron-down').addClass('glyphicon-chevron-up');
            }
            else {
                // collapse the panel
                \$(this).parents('.panel').find('.panel-body').slideUp();
                \$(this).addClass('panel-collapsed');
                \$(this).find('i').removeClass('glyphicon-chevron-up').addClass('glyphicon-chevron-down');
            }
        });
    });
</script>		
};
}

sub print_header_genes_viewer {
	my ($cgi,$title) = @_;
	print $cgi -> header;
	$title = "PolyDiag " unless $title; 
		print $cgi -> start_html(-title=>$title,
						-style=>{
						
							-src=>['polyweb/css/polydiag/polydiag.css',
										"//maxcdn.bootstrapcdn.com/bootswatch/3.3.4/flatly/bootstrap.min.css",
										"//maxcdn.bootstrapcdn.com/font-awesome/4.3.0/css/font-awesome.min.css",
										#"//cdnjs.cloudflare.com/ajax/libs/bootstrap-table/1.7.0/bootstrap-table.min.css"
									 "//ajax.googleapis.com/ajax/libs/dojo/1.10.4/dojo/resources/dojo.css",
			 "//ajax.googleapis.com/ajax/libs/dojo/1.13.0/dijit/themes/claro/claro.css",
			  "//ajax.googleapis.com/ajax/libs/dojo/1.13.0/dojox/grid/resources/claroGrid.css",
			  "//ajax.googleapis.com/ajax/libs/dojo/1.13.0/dojox/editor/plugins/resources/css/Save.css",
			 "//ajax.googleapis.com/ajax/libs/dojo/1.13.0/dojox/widget/Dialog/Dialog.css",	
			 "//ajax.googleapis.com/ajax/libs/dojo/1.13.0//dijit/themes/claro/claro.css",
     		 "//ajax.googleapis.com/ajax/libs/dojo/1.13.0//dijit/themes/claro/document.css",
										],
							},
							
							 -script=>{
							 #	-src=>  "//code.jquery.com/jquery-1.11.2.min.js",
							 	-src=>"//ajax.googleapis.com/ajax/libs/dojo/1.13.0/dojo/dojo.js",
							 	#-src=>"http://darwin.bipd.fr/aptana/polyweb/js/polyjs/Utility/phenotypes.js",
                             	#-src=>  "//code.jquery.com/jquery-1.11.2.min.js",
                             },
                             );
                             
   # print qq{<script type="text/javascript" src="//maxcdn.bootstrapcdn.com/bootstrap/3.3.2/js/bootstrap.min.js""></script> };
     print qq{<script src="//ajax.googleapis.com/ajax/libs/dojo/1.13.0/dojo/dojo.js" djConfig="parseOnLoad:true" type="text/javascript"></script> };
    print qq {<script src="/../../aptana/polyweb/js/polyjs/Utility/phenotypes.js"> </script>};
    print qq{<script type="text/javascript" src="http://darwin.bipd.fr/aptana/polyweb/js/polyjs/include.polyviewer.js" djConfig="parseOnLoad:true"></script>};
	print qq {
		<script src="https://cdn.anychart.com/releases/v8/js/anychart-base.min.js"></script>
		 <script src="https://cdn.anychart.com/releases/v8/js/anychart-tag-cloud.min.js"></script>	
	};
	print qq{
	  <link rel="stylesheet" href="//maxcdn.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css">
	    <link rel="stylesheet" href="//maxcdn.bootstrapcdn.com/bootswatch/3.4.1/flatly/bootstrap.min.css">
 
	};
print qq {		
		<style>
		table-mybordered > thead > tr > th, .table-mybordered > thead > tr > th, table-mybordered > tbody > tr > th, .table-mybordered > tbody > tr > th, table-mybordered > tfoot > tr > th, .table-mybordered > tfoot > tr > th, table-mybordered > thead > tr > td, .table-mybordered > thead > tr > td, table-mybordered > tbody > tr > td, .table-mybordered > tbody > tr > td, table-mybordered > tfoot > tr > td, .table-mybordered > tfoot > tr > td {
  	border: 1px solid #95A5A6;
  vertical-align:middle;
   text-align:center; 
   padding:3px;
}

		</style>	
<script type="text/javascript">
	dojo.require("dojo.data.ItemFileWriteStore");
    
</script>	
 	
};
#<body class="claro" style ="font-size : 100%;">
}

sub print_cgi {
	my ($cgi,$out,$print,$title) = @_;
	print qq{
		<script>
		dojo.connect(_container_, "onLoad", function(){
   /* your script */
});
		alert("cuicui");
window.addEventListener('load', function () {
  alert("xx");
});
window.addEventListener("DOMContentLoaded", function(){
  alert("xxX");
});
document.onload = alert("cocuou");
$(document).ready(function(){ alert("jq") }) 
dojo.addOnLoad(function() {	
	dijit.byId('waiting').show();
});
dijit.byId('waiting').show();
</script>
	};
	return;
	print_header_polydiag2($cgi,$title);
	if ($print){
		print_header_polydiag_printer($cgi,$title);
	}
	else {
		print_header_polydiag($cgi,$title);
	}
	
	print $out;
	
	exit(0);
}

sub print_cgi_header {
	my ($cgi,$out,$print,$title) = @_;
	
	if ($print){
		print_header_polydiag_printer($cgi,$title);
	}
	else {
		print $cgi -> header;
	#die();
		return;
		print_header_polydiag($cgi,$title);
	}
}

sub print_header_table_version {
	my ($cgi,$title) = @_;
	print $cgi -> header;
	$title = "PolyDiag " unless $title; 
		print $cgi -> start_html(-title=>$title,
						-style=>{
						
							-src=>['/polyweb/css/polydiag/polydiag.css',
										"//maxcdn.bootstrapcdn.com/bootswatch/3.3.4/flatly/bootstrap.min.css",
										"//maxcdn.bootstrapcdn.com/font-awesome/4.3.0/css/font-awesome.min.css",
										#"//cdnjs.cloudflare.com/ajax/libs/bootstrap-table/1.7.0/bootstrap-table.min.css"
										
										],
							},
							
							 -script=>{
                             	-src=>  "//code.jquery.com/jquery-1.11.2.min.js",
                             },
                             );
    print qq{<script type="text/javascript" src="//maxcdn.bootstrapcdn.com/bootstrap/3.3.2/js/bootstrap.min.js""></script> };
    print qq{
    	<link rel="stylesheet" href="//cdnjs.cloudflare.com/ajax/libs/bootstrap-table/1.7.0/bootstrap-table.min.css">
    	<script src="//cdnjs.cloudflare.com/ajax/libs/bootstrap-table/1.7.0/bootstrap-table.min.js"></script>
    	<script src="//cdnjs.cloudflare.com/ajax/libs/bootstrap-table/1.7.0/locale/bootstrap-table-zh-CN.min.js"></script>
    };               
					
	print qq{		
			<style>
			table-mybordered > thead > tr > th, .table-mybordered > thead > tr > th, table-mybordered > tbody > tr > th, .table-mybordered > tbody > tr > th, table-mybordered > tfoot > tr > th, .table-mybordered > tfoot > tr > th, table-mybordered > thead > tr > td, .table-mybordered > thead > tr > td, table-mybordered > tbody > tr > td, .table-mybordered > tbody > tr > td, table-mybordered > tfoot > tr > td, .table-mybordered > tfoot > tr > td {
		  border: 1px solid #95A5A6;
		  vertical-align:middle;
		   text-align:center; 
		   padding:2px;
		}
	
	
		.pricing .card {
		  background: #FFFFFF;
		  border: none;
		  border-radius: 1rem;
		  transition: all 0.2s;
		  box-shadow: 0 0.5rem 1rem 0 rgba(0, 0, 0, 0.1);
		}
		
		.pricing hr {
		  margin: 1.5rem 0;
		}
		
		.pricing .card-title {
		  margin: 0.5rem 0;
		  font-size: 0.9rem;
		  letter-spacing: .1rem;
		  font-weight: bold;
		}
		
		.pricing .card-price {
		  font-size: 3rem;
		  margin: 0;
		}
		
		.pricing .card-price .period {
		  font-size: 0.8rem;
		}
		
		.pricing ul li {
		  margin-bottom: 1rem;
		}
		
		.pricing .text-muted {
		  opacity: 0.7;
		}
		
		.pricing .btn {
		  font-size: 80%;
		  border-radius: 5rem;
		  letter-spacing: .1rem;
		  font-weight: bold;
		  padding: 1rem;
		  opacity: 0.7;
		  transition: all 0.2s;
		}
		
		/* Hover Effects on Card */
		
		\@media (min-width: 992px) {
		  .pricing .card:hover {
		    margin-top: -.25rem;
		    margin-bottom: .25rem;
		    box-shadow: 0 0.5rem 1rem 0 rgba(0, 0, 0, 0.3);
		  }
		  .pricing .card:hover .btn {
		    opacity: 1;
		  }
		}	

		.bs-callout {
		    padding: 20px;
		    margin: 20px 0;
		    border: 1px solid #eee;
		    border-left-width: 5px;
		    border-radius: 3px;
		}
		.bs-callout h4 {
		    margin-top: 0;
		    margin-bottom: 5px;
		}
		.bs-callout p:last-child {
		    margin-bottom: 0;
		}
		.bs-callout code {
		    border-radius: 3px;
		}
		.bs-callout+.bs-callout {
		    margin-top: -5px;
		}
		.bs-callout-default {
		    border-left-color: #777;
		}
		.bs-callout-default h4 {
		    color: #777;
		}
		.bs-callout-primary {
		    border-left-color: #428bca;
		}
		.bs-callout-primary h4 {
		    color: #428bca;
		}
		.bs-callout-success {
		    border-left-color: #5cb85c;
		}
		.bs-callout-success h4 {
		    color: #5cb85c;
		}
		.bs-callout-danger {
		    border-left-color: #d9534f;
		}
		.bs-callout-danger h4 {
		    color: #d9534f;
		}
		.bs-callout-warning {
		    border-left-color: #f0ad4e;
		}
		.bs-callout-warning h4 {
		    color: #f0ad4e;
		}
		.bs-callout-info {
		    border-left-color: #5bc0de;
		}
		.bs-callout-info h4 {
		    color: #5bc0de;
		}
		.badge1 {
		  padding: 1px 9px 2px;
		  font-size: 12.025px;
		  font-weight: bold;
		  white-space: nowrap;
		  color: #ffffff;
		  background-color: #999999;
		  -webkit-border-radius: 9px;
		  -moz-border-radius: 9px;
		  border-radius: 9px;
		}
				</style>	
		<script type="text/javascript">
		
			function zoomCoverage(tr){
			  	var start_url = "/cgi-bin/polymorphism-cgi/polydiag_old/coverage_exons.pl";
			 	var url = construct_url(start_url,tr);
			 	dijit.byId("popup_coverage").show();
			 	dijit.byId("div_popup_coverage").setHref(url);
			}
		
		    jQuery(function (\$) {
		        
		        \$('.panel-heading span.clickable').on("click", function (e) {
		            if (\$(this).hasClass('panel-collapsed')) {
		                // expand the panel
		                \$(this).parents('.panel').find('.panel-body').slideDown();
		                \$(this).removeClass('panel-collapsed');
		                \$(this).find('i').removeClass('glyphicon-chevron-down').addClass('glyphicon-chevron-up');
		            }
		            else {
		                // collapse the panel
		                \$(this).parents('.panel').find('.panel-body').slideUp();
		                \$(this).addClass('panel-collapsed');
		                \$(this).find('i').removeClass('glyphicon-chevron-up').addClass('glyphicon-chevron-down');
		            }
		        });
		    });
		</script>		
		}
}


1;