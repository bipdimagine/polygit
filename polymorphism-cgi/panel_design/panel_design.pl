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

use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);
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
my $buffer = GBuffer->new();

 my $fsize = "font-size:10px";
my $cgi          = new CGI();
my $gencode = $cgi->param('gencode');
my $name = $cgi->param('name');

my $project = $buffer->newProjectCache( -name 			=> "NGS2021_4169");
$project->gencode_version($gencode);# if $gencode;
 	 my $capture_intspan_keep;

	 my $out;
	# my  @transcripts_cgi = @{$project->bundle_transcripts() } ;
	 my $padding = 1;
	 my $bed_file= "/data-isilon/public-data/design/$name.bed.gz";
	 my $hcapture_intspan ;#=  Set::IntSpan::Fast::XS->new();
	 open (BED,"zcat $bed_file | ");
	 	
	foreach my $l (<BED>){ 
		chomp ($l);
		my ($chr_name,$start,$end,@reste) = split(" ",$l);
		my $chr = $project->getChromosome($chr_name);
		$hcapture_intspan->{$chr->name} = Set::IntSpan::Fast::XS->new() unless exists  $hcapture_intspan->{$chr->name};
		$capture_intspan_keep->{$chr->name} = Set::IntSpan::Fast::XS->new()  unless exists  $capture_intspan_keep->{$chr->name};
		$hcapture_intspan->{$chr->name}->add_range($start-$padding,$end+$padding);
		$capture_intspan_keep->{$chr->name}->add_range($start-$padding,$end+$padding);
	 }
	 
	
	my @list_transcripts = `cat /data-isilon/public-data/design/$name.list`;
	 chomp(@list_transcripts);
	 
	 my $problem;
	 my $nbp;
	
	 #ENST00000450892
	 my $notfound;
	 foreach my $tr_name  (@list_transcripts){
	 #	next 
	  my $debug ;
	  $debug = 1 if $tr_name eq "ENST00000464845";
	#  warn $tr_name;
	  	my $transcript;
	   eval {
	 	 $transcript =  $project->newTranscript($tr_name);
	   };
	   push(@$notfound,$tr_name) unless $transcript;
	   next unless $transcript;
	# 	 $transcript->getChromosome->
	 	my $chr = $transcript->getChromosome;
	 	my $capture_intspan = $hcapture_intspan->{$chr->name};
	 	
	 #	die() if $debug ;
	 	
	 	foreach my $exon (sort{$a->start <=> $b->start} @{$transcript->getExons}){
	 		warn $exon->start."-".$exon->end if $debug; 
	 		my $span =  Set::IntSpan::Fast::XS->new($exon->start."-".$exon->end );
	 		$capture_intspan_keep->{$chr->name}->remove_range($exon->start-100,$exon->end+100);
	 		next if $exon->is_noncoding;
	 		my $s1 = $exon->getGenomicSpan()->intersection($capture_intspan);
	 		 if ($s1->is_empty){
	 		 push (@{$problem->{problem}->{$transcript->name}},$exon);
	 		 $nbp ++;
	 		 }
	 		
	 		 
	 		#$covered = -1 if $s1->is_empty;
	 	}
	 }
	# 12: 44,152,747-44,183,346
	 my $hgenes;
	 foreach my $chr_name (keys %{$capture_intspan_keep}){
	 	my $iter = $capture_intspan_keep->{$chr_name}->iterate_runs();
	 	my $chr = $project->getChromosome($chr_name);
		my @tt;
    	while (my ( $from, $to ) = $iter->()) {
    		#next if ($from-$to) < 200;
    		#warn $chr_name." ".$from."  ".$to." ".abs($from-$to);
    		my $genes = $chr->getGenesByPosition($from,$to);
    		
    		foreach my $g (@$genes){
    			my $htranscripts = {};
    			
    				my $trs = $chr->getTranscriptsByPosition($from,$to);
    				
    				foreach my $t (@$trs){
    						my $span =  Set::IntSpan::Fast::XS->new($from."-".$to);
    					 	my $xx = $t->getGenomicSpan()->intersection($span);
    					 	next if $xx->is_empty;
    					 	#warn $xx->as_string() if $t->name eq "ENST00000450892";
    							$htranscripts->{$t->id} = $xx->as_string;
    				}
    				#warn $chr_name." ".$from."  ".$to." ".abs($from-$to);
    			if (%$htranscripts){
    			$hgenes->{$g->id}->{obj} = $g;
    			$hgenes->{$g->id}->{nb} ++;
    			foreach my $t (keys %$htranscripts){
    				push(@{$hgenes->{$g->id}->{transcripts}->{$t}} , $htranscripts->{$t});
    				#$hgenes->{$g->id}->{transcripts} .= $t.":".$htranscripts->{$t}."<br>";
    			}
    			}
    		 	#warn $g->external_name;
    		}
    	}
    }

		

	  my 	$style = "background-color:#DD4132;color:black;";
	my	  		$shadow = "rgba(221, 65,36, .7)";
	#my $out="";
	# 	$out = $cgi->start_div({class=>"panel-body panel-primary  panel-collapse collapse ",style=>"font-size: 09px;font-family:  Verdana;border: 5px coral; border-color: coral;",id=>"control_design_".$run->id});	
	$out .=$out .= $cgi->start_table({class=> "table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 10px;font-family:  Verdana;box-shadow: 1px 1px 2px 1px YYY ;", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
	 		$out .= $cgi->start_Tr({style=>"font-weight: bold;XXX "});
	 		$out .= qq{
					<td colspan="2"> Transcript Problem </td>
					</tr> 
				};
		foreach my $t (@$notfound){
			$out .= $cgi->start_Tr();
				$out .= qq{<td colspan="2">$t</td></tr>};
		}
		$out.="</table>";
				
	   $out .= qq{
		<div class="container">
		<div class="row">
	   };
	   my @gs = grep{$hgenes->{$_}->{nb}>1} keys %{$hgenes};
	#   @gs =();
	   if (@gs){
	   	$nbp ++;
	   	$out .= qq{
	 	 <div>
			
	 		};
	 		
	 		$out .= $cgi->start_table({class=> "table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 10px;font-family:  Verdana;box-shadow: 1px 1px 2px 1px YYY ;", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
	 		$out .= $cgi->start_Tr({style=>"font-weight: bold;XXX "});
	 		$out .= qq{
					<td colspan="2"> Potential missing genes </td>
					</tr> 
				};
				
		#my @gs = grep{$hgenes->{$_}->{nb}>1} keys %{$hgenes};
		my $max =0;
	   foreach my $gid (@gs){
	   			$out .= $cgi->start_Tr();
				$out .= qq{<td colspan="2">}.$hgenes->{$gid}->{obj}->external_name." miss : ".$hgenes->{$gid}->{nb}." exons <br>";
				foreach my $t (keys %{$hgenes->{$gid}->{transcripts}}) {
					$out.= $t." ".join(",",@{$hgenes->{$gid}->{transcripts}->{$t}})."<br>";
				}
				$out.="</td>";
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
	 	 <div >
			
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
	html::print_cgi_header($cgi,"",1," - PolyDiag");
	print "<h1> GENCODE : $gencode</h1>";
	 # $out2 =  qq{<div class="btn   btn-xs $btn " style="position:relative;bottom:1px;min-width:200px;border-color:black;" onClick='collapse_panel("control_design","$list_control_panels","$run_id")'  style="border : 1px"  $disabled>  $text $label </div>};
	  print $out;
	  
	  