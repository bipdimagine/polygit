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


my $types = {
							he=>2,
							ho=>3,
							"-"=>0,
							nc=>-5
							
						};

my $exon_types = {
							he=>-1,
							ho=>-2,
							'-'=>0,
							ok=>1,
						};						
						my $hvalidation = {
							2=>"he",
							3=>"ho",
							#-3=>"td",
							-5=>"nc",
							
						};
my $nb_gene_by_patient = 3;
my $nb_exon_by_genes = 10;

my $buffer = GBuffer->new();

my $cgi          = new CGI();
my $project_name = $cgi->param('project');

my $project = $buffer->newProject(-name=>$project_name);
my $capture = $project->getCapture();

my $vquery = validationQuery->new(dbh=>$buffer->dbh,capture_name=>$capture->validation_db());
my $patient_name = $cgi->param('patients');
$patient_name ="all" unless $patient_name;
my $patients = $project->get_list_patients($patient_name,",");
my $gene =  $cgi->param('gene');

my $cgi_transcript =  $cgi->param('transcripts');
my @transcripts_cgi ;
if ($cgi_transcript eq "all"){
	@transcripts_cgi = @{$capture->transcripts_name()} ;
}
else {
	@transcripts_cgi = split(",",$cgi_transcript);
}

my $chromosomes = $project->getChromosomes();
my $print =  $cgi->param('print');

my $splice_5 = 10;

my $htranscripts;
my $todo_patients =[];
foreach my $patient (sort{$a->name cmp $b->name} @{$patients}){
	my $variations_todo = $vquery->get_variations_todo(project_name=>$project_name,sample_name=>$patient->{name});#project_name=>$project_name,id=>$hexons->{label}->{$patient->name()});
	my $exons_todo = $vquery->get_exons(project_name=>$project_name,sample_name=>$patient->{name});
	if ($variations_todo){
		push(@$todo_patients,$patient);
	}
	elsif ($exons_todo) {
		 	push(@$todo_patients,$patient);
		
	}
	my $gvariations =0;
	foreach my $tr_name (@transcripts_cgi){ 
		my $tr_obj = $project->newTranscript($tr_name);
		my $tr = $htranscripts->{$tr_name};
		my $exons = $tr_obj->getExons();
		my $nb_exons =0;
		$tr->{lname} = $tr_obj->getGene()->name." ".$tr_obj->name." :".$tr_obj->getGene()->external_name;
		
		my @variations_tr;
			foreach my $v (values %$variations_todo){
					my $chr_v = $project->getChromosome($v->{chromosome});
 				next if ($chr_v->name ne $tr_obj->getChromosome->name) ;
		#		next if $v->{chromosome} ne $tr_obj->getChromosome->ucsc_name;
				next if $v->{validation} ne -3;
				
				next if  $v->{validation_sanger} ne 0; 
				next unless $tr_obj->span_position_genomic->contains($v->{start});
				my ($pos,$exon_name) = split("_", $tr_obj->findNearestExon($v->{start}));
				$v->{exon} = $exon_name;
				$v->{pos_exon} = $pos;
				
				push(@variations_tr,$v);
			}
		foreach my $exon (@$exons) {
			my $id_exon = $patient->name().":".$tr_obj->getGene()->external_name.":".$tr_obj->name.":".$exon->name; 
			$id_exon .= ":".$tr_obj->getChromosome()->name.":".$exon->start.":".$exon->end;
			$exon->{label}->{$patient->name()} = $patient->name().":".$tr_obj->getGene->external_name().":".$tr_obj->name.":".$exon->name.":".$tr_obj->getChromosome()->name.":".$exon->start.":".$exon->end;
		
			my $todo_coverage;
			if (exists $exons_todo->{$exon->coverage_id}){
				if ($exons_todo->{$exon->coverage_id}->{done} == 0){
					$htranscripts->{$tr_name}->{$exon->name()}->{resequencing}->{$patient->name}->{coverage} = $id_exon ;
			 		$htranscripts->{$tr_name}->{$exon->name()}->{resequencing}->{$patient->name}->{level} = 1 ;
			 		$htranscripts->{$tr_name}->{has_todo}=1;
			 		$htranscripts->{$tr_name}->{$exon->name()}->{nb_resequencing}->{$exon->id}++;
			 	
					
				}
			}
			my $variations;
			
			foreach my $v (@variations_tr){
				next if $v->{exon} ne $exon->name();
				$htranscripts->{$tr_name}->{has_todo} = 1;
				 $v->{exon_start} = $v->{start} - $exon->{gstart};
					push(@$variations,$v);
					$htranscripts->{$tr_name}->{$exon->name()}->{variations_resequencing}->{$v->{id}}->{$patient->{name}} =$v;
					$htranscripts->{$tr_name}->{$exon->name()}->{nb_resequencing}->{$v->{id}}++;
			}
				 $htranscripts->{$tr_name}->{$exon->name()}->{resequencing}->{$patient->name}->{variations} = $variations if $variations;
			
		} # end exons

	}
}#end patient	
html() unless $cgi->param('xls');
xls();



sub html{
my $out="";
my $nb_patients = scalar(@$todo_patients);
if ($print){
$out.="<body onload='this.print();'>";
}
my $header_title;
$header_title = "Todo :  Sample : ".$nb_patients ." he : HEterozygous sanger confirmed &nbsp; ho : HOmozygous sanger cofirmed &nbsp; nc : Not Confrmed in sanger" ;
unless($print){
$header_title .= $cgi->start_div({class=>"btn-group",role=>"group",'aria-label'=>"...",style=>"float:right"});

$header_title .='<a type="button" class="btn btn-xs btn-success" onclick="todo_print()"><i class="fa fa-print  pull-left  "></i>Print</button></a>';
$header_title.='<a type="button" class="btn btn-xs btn-info" onclick="todo_xls()"><i class="fa fa-file-excel-o"></i> Export</button></a>';
$header_title.=' <a type="button" class="btn btn-xs btn-danger" onclick="show_save();"><i class="fa fa-floppy-o"></i> Save</button></a>';
$header_title.=$cgi->end_div();
}
$out .= html::print_cadre($cgi,$header_title);
my $zz =0;
my $tdid =0;
my $exjsid =0;
my $nb = scalar(@$patients);
#$out.= 	;
my @class = ("success","info");
my @colors =("#FFD827","#EAF3F3");
my @colors =("#E8EDFF","#EAF3F3");
#my @colors =("#18BC9C","#217DBB");
my $xt=0;
			$out.=  $cgi->start_table({class=>"table  table-condensed table-bordered table-mybordered",style=>"border-color:red;vertical-align: middle;font-size: 8px;font-family:  Verdana;width:1800px"});
#			$out.= $cgi->start_Tr();
#			$out.= $cgi->th({colspan=>3},"Genes Infos");
#			$out.= $cgi->th({colspan=>3},"Exons Infos");
#			$out.= $cgi->th({colspan=>2},"Validation type");
#			$out.= $cgi->th({colspan=>scalar(@$patients)},"Samples");
#			$out.= $cgi->end_Tr();
			my @ths = ("gene","transcript","chr","exon","start","end","type","Allele","Position") ;
			my $th_string = "";
			$th_string.= $cgi->start_Tr({class=>"success",style=>"border-color:red;"});
			$th_string.= $cgi->th({style=>"width:0.5em"},"gene ");
			$th_string.= $cgi->th({style=>"width:0.1em"},"chr");
			$th_string.= $cgi->th({style=>"width:1.4em"},"transcripts");
			
			$th_string.= $cgi->th({style=>"width:0.05em;"},"exon ");
			$th_string.= $cgi->th("start");
			$th_string.= $cgi->th("end");
			$th_string.= $cgi->th({style=>"width:0.4em"},"type");
			$th_string.= $cgi->th("Allele");
			$th_string.= $cgi->th("Position");
			my $colspan = 9;
	foreach my $patient (@$todo_patients){
		$colspan ++;
		push(@ths,$patient->{name});
			$th_string.= $cgi->th($patient->{name});
		}
		$th_string.= $cgi->end_Tr();
	
		my @transcripts;
			
foreach my $tr_name (@transcripts_cgi){ 

	next unless exists 	 $htranscripts->{$tr_name}->{has_todo};
	push (@transcripts,$project->newTranscript($tr_name));
}
foreach my $tr_obj (sort{$a->getGene->external_name cmp $b->getGene->external_name} @transcripts){
		$out .= $th_string;
		$xt++;
	#$out.= $cgi->start_div({class=>"panel panel-default",style=>"font-size: 10px;font-family:  Verdana;"});
	#$out.= $cgi->div({class=>"panel-heading",style=>"font-size: 11px;font-family:  Verdana;background-color:#DFF0D8",id=>},$tr_obj->getGene()->external_name."::".$tr_name." ".$tr_obj->getChromosome->ucsc_name.":".$tr_obj->start."-".$tr_obj->end);
	my $tr_name = $tr_obj->name();
$out.= $cgi->start_Tr({nowrap=>1,style=>"vertical-align: middle;"."background-color:".$colors[$xt%2].";"});
	my $htr =  $htranscripts->{$tr_name};
	my @exons = grep {exists $htr->{$_->name()}->{resequencing}} @{$tr_obj->getExons()};
	my $span = scalar(@exons);
	my $nbv ;
	map{ $nbv += scalar (keys %{$htr->{$_->name()}->{nb_resequencing}}) } @exons;
	$span = $nbv;
	#$span -=1;
	my $st = {rowspan=>$span,nowrap=>1,style=>"vertical-align: middle;"};
	$out.= $cgi->td($st,$tr_obj->getGene()->external_name);	
	$out.= $cgi->td({rowspan=>$span,nowrap=>1,style=>"vertical-align: middle;max-width:10px"},$tr_obj->getChromosome->name);
	
	$out.= $cgi->td($st,$tr_name);
	
		
			
			foreach my $exon (sort {$a->start <=> $b->end} @exons){
			
				
				next unless exists $htr->{$exon->name()}->{resequencing};
				
				my $nbv = 	scalar (keys %{$htr->{$exon->name()}->{nb_resequencing}});
				my $hexon = $htr->{$exon->name()}->{resequencing};
				
				my $line = 0;
				my $hvariations = $htranscripts->{$tr_name}->{$exon->name()}->{variations_resequencing};
				my $st = {rowspan=>$nbv,nowrap=>1,style=>"vertical-align: middle;"};
				

			
			
				$out.= $cgi->td($st,$exon->name);
				$out.= $cgi->td($st,$exon->start);
				$out.= $cgi->td($st,$exon->end);
			#	$out.= $cgi->td($st,"[".$exon->start."-".$exon->start."]");
				 if (exists $htr->{$exon->name()}->{nb_resequencing}->{$exon->id}){
				 	$out .= $cgi->td("exons").$cgi->td("&nbsp;"). $cgi->td("&nbsp");
				 	foreach my $patient (@$todo_patients){
				 		 if (exists $hexon->{$patient->name}->{coverage}){
				 		 		$tdid++;
				 		 		my $td_uniq_id = "td_todo_".$tdid;
				 				my $value = 	$exon->{label}->{$patient->name()};
								my $label =$exon->name; 
								$exjsid ++;
								my $vz = "ex_".$exjsid;
								my $text;
								unless ($print){
								 $text =qq{<div  style="vertical-align:middle;padding:1px;border-radius: 5px; -moz-border-radius: 5px; -webkit-border-radius: 5px; border: 2px solid #FFCF13;" >};
								 $text .= qq{<select onChange="put_validation_exons2(this,'$td_uniq_id','$value')" style="color:black;" >};
								foreach my $type (keys %$exon_types){
									my $a = $exon_types->{$type};
									my $st ="";
									$st="selected"  if $a == 0;
										$text .=  qq{<option  value="$a"  $st>$type</option>};
							}
				 		
						
								}
								else {
										$text = "";
								foreach my $type (keys %$exon_types){
									next if $type eq "-";
									$text .= "<i class='a fa-square-o'>$type</i> "
								}
								}
								
								$out .= $cgi->td( {id=>$td_uniq_id},$text);
				 		 }
				 		 else {
				 		 	$out .= $cgi->td("-");
				 		 }
				 	}
				 	$out.= $cgi->end_Tr();
				 	$out.= $cgi->start_Tr({nowrap=>1,style=>"vertical-align: middle;"."background-color:".$colors[$xt%2].";"});
				 }
				foreach my $vid (keys %$hvariations){
				
					my $outemp ="";
					my $td_all;
					foreach my $patient (@$todo_patients){
						$tdid++;
						my $td_uniq_id = "td_".$tdid;
						if (exists $hvariations->{$vid}->{$patient->{name}}){
							my $v = $hvariations->{$vid}->{$patient->{name}};
							my $vsid =  $hvariations->{$vid}->{$patient->{name}}->{validation_id};
						
							my $v_text;
							unless ($print){
							 $v_text = qq{<div  style="vertical-align:middle;padding:1px;border-radius: 5px; -moz-border-radius: 5px; -webkit-border-radius: 5px; border: 2px solid #18BC9C;" >};
							 $v_text .= qq{<select   name="$vsid" onChange="put_variations_status(this,'$td_uniq_id')" style="color:black;padding:1px" >};
								$v_text .=  qq{<option  value="none">-</option>};
								foreach my $type (keys %$types){
										$v_text .=  qq{<option  value="$type" ">$type</option>};
									}
				 			$v_text .= "</select></div>";
							}
							else {
								$v_text = "";
								foreach my $type (keys %$types){
									$v_text .= "<i class='a fa-square-o'>$type</i> "
								}
							}
							$outemp.= $cgi->td({id=>$td_uniq_id},$v_text);
							my $v = $hvariations->{$vid}->{$patient->{name}};
							$td_all =  $cgi->td("variant").$cgi->td($v->{ref}."/".$v->{alt}).$cgi->td($v->{start});
						}
						else {
							$outemp.= $cgi->td("-");
						}
					}
					$out .= $td_all.$outemp;
					$out.= $cgi->end_Tr();
					
					$out.= $cgi->start_Tr({nowrap=>1,style=>"vertical-align: middle;"."background-color:".$colors[$xt%2].";"});
				}
		
				
			
			}#end for exons
		$out.= $cgi->start_Tr();
			
				
		$out.= $cgi->td({colspan=>$colspan,class=>"primary"},"&nbsp;");
			$out.= $cgi->end_Tr();	

}#end for transcript
	$out.=$cgi->end_table();
	$out.= $cgi->end_div();
$out .= html::end_cadre($cgi);


html::print_cgi($cgi,$out,$print);
exit(0);
#html::print_cgi($cgi,$out);

}




sub xls{

	print "Content-type: application/msexcel\n";
	print "Content-Disposition: attachment;filename=".$project->name().".xls\n\n";
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
                     
	my $row_patient = 1;
	my $col = 0; 
	my $row =0;
	my $nb_patients = scalar(@$todo_patients);
#print $cgi->start_table({caption=>$tr->{name}, class=>"table1"});
#print header 
$worksheet->write($row,$col++,"Gene");
$worksheet->write($row,$col++,"chr");
$worksheet->write($row,$col++,"transcript");
$worksheet->write($row,$col++,"exon");
$worksheet->write($row,$col++,"start");
$worksheet->write($row,$col++,"end");
$worksheet->write($row,$col++,"type");
$worksheet->write($row,$col++,"allele");
$worksheet->write($row,$col++,"position");
$worksheet->write($row,$col++,"sequence");
	foreach my $patient (@$todo_patients){
	$worksheet->write($row,$col++,$patient->{name});
		}
$row++;		
$col=0;
foreach my $tr_name (@transcripts_cgi){ 
	next unless exists 	 $htranscripts->{$tr_name}->{has_todo};
	my $tr_obj = $project->newTranscript($tr_name);
		my $htr =  $htranscripts->{$tr_name};
		

		
			#$worksheet->write($row,$col++,);
		my $exons = $tr_obj->getExons();
	
	
		
		my @texons = grep {exists $htr->{$_->name()}->{resequencing}} @{$tr_obj->getExons()};
			
		foreach my $exon (@texons) {
		next unless exists $htr->{$exon->name()}->{resequencing};
		my $hvariations = $htranscripts->{$tr_name}->{$exon->name()}->{variations_resequencing};
		if (exists $htr->{$exon->name()}->{nb_resequencing}->{$exon->id}){
			my $hexon = $htr->{$exon->name()}->{resequencing};
			$worksheet->write($row,$col++,$tr_obj->getGene->external_name);
			$worksheet->write($row,$col++,$tr_obj->getChromosome()->name());
			$worksheet->write($row,$col++,$tr_obj->{name});
			$worksheet->write($row,$col++,$exon->name);
			$worksheet->write($row,$col++,$exon->start);
			$worksheet->write($row,$col++,$exon->end);
			$worksheet->write($row,$col++,"exon");
			$worksheet->write($row,$col++,"-");
			$worksheet->write($row,$col++,"-");
			$worksheet->write($row,$col++,"-");
			foreach my $patient (@$todo_patients){
				 		 if (exists $hexon->{$patient->name}->{coverage}){
				 		 	$worksheet->write($row,$col++,"TODO");
				 		 }
				 		 else{
				 		 	$worksheet->write($row,$col++,"OK");
				 		 }
				 		 
			}
			$row ++;
			$col=0;
			}
			
			foreach my $vid (keys %$hvariations){
				my @vs = values %{$hvariations->{$vid}};
				my $v = $vs[0] ;
				my $chr = $project->getChromosome($v->{chromosome});
				$worksheet->write($row,$col++,$tr_obj->getGene->external_name);
				$worksheet->write($row,$col++,$tr_obj->getChromosome()->name());
				$worksheet->write($row,$col++,$tr_obj->{name});
				$worksheet->write($row,$col++,$exon->name);
				$worksheet->write($row,$col++,$exon->start);
				$worksheet->write($row,$col++,$exon->end);
				$worksheet->write($row,$col++,"variation");
				$worksheet->write($row,$col++,$v->{ref}."/".$v->{alt});
				$worksheet->write($row,$col++,$v->{start});
				$worksheet->write($row,$col++,$chr->sequence($v->{start}-5,$v->{end}+5));
					foreach my $patient (@$todo_patients){
						
						if (exists $hvariations->{$vid}->{$patient->{name}}){
							$worksheet->write($row,$col++,"TODO");
							
							}
						else {
							$worksheet->write($row,$col++," ");
						}#if
					}#pat
					$row ++;
					$col =0;
				}#var
			
			
		}#end for exon		

 }#end for trnascript
	
	
	
	
                   
} 

	
	