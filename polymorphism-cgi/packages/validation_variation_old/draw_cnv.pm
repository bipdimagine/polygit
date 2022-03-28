package draw_cnv;
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;


my $id =0;
sub table_cnv_transcripts{
	my ($cgi,$patient,$tr) =@_;
	
	my $project = $patient->getProject();
my $max = 50;
my $print_out ="";
my $xr = 255;
my $xg = 10;
my $xb = 0;

my $yr = 0;
my $yg = 255;
my $yb = 0;
my $jsid=0;
my $nid=0;
my $cnv_patients;
my $tab_ids;

my $tooltip_name = "bubble_".$patient->id."_".$tr->id."_".$id;	

	$print_out.= $cgi->start_Tr();
	$print_out.= $cgi->th({scope=>"col", abbr=>"Configurations", class=>"nobg"}, "Gene ".$tr->getGene()->external_name);
	foreach my $primer (sort{$b->getChromosome->length <=> $a->getChromosome->length || $a->start <=> $b->start } @{$tr->getPrimers()} ){
		my $text_exons;
		foreach my $e (@{$primer->getExons}){
			$text_exons.= $e->name()." ";
			
		}
		my $pn = $patient->{name};
			my $exon = $primer->{$patient->name}->{exon};
			my $value = "";#$patient->{name}.":".$exon->{label};
			my $mean = 0;#int($primer->mean($patient));
			my $min = 0;#$primer->min($patient);
			my $text = "<small>".$primer->mean($patient)."</small>";
		
			#$val = 100 - $val;
			
			my $label =$exon; 
			my $text = $mean."(".$min.")";
		
			my $en = $exon;
			
			my $mean = 0;#int($primer->statistics()->mean()*100)/100;
			#my $text = $exon->{mean}."<small>(".$exon->{min}.")<br>[".$exon->{intspan}."]</small>";
			$jsid++;
			my $check;
			$text = $cgi->p($text);
				#$text = qq{<div class="openGraph">$text</div>} ;
				my $pn = $patient->name();
				my $en = $primer->name();
				my $click = qq{load_graph_primer('$pn','$en');};
				my $plex = $primer->multiplex;
				my $score = $primer->cnv_score($patient);
				
				my $sdp= 0;#int($primer->statistics()->standard_deviation()*100)/100;
				my $sdm = 0;#int($patient->statistics_multiplex($primer->multiplex)->standard_deviation()*100)/100;;
				my $meanm = 0;#int($patient->statistics_multiplex($primer->multiplex)->mean()*100)/100;
				my $zscore = $primer->zscore($patient);
				my $itext = $score." +/- ".$sdp." -  $mean ; multi: $meanm +/-  $sdm ; zscore  ".$zscore;
				$id++;
				#my $id = $tr->id."_".$patient->id()."_".$nid;
					my $level = $primer->level($patient);
				my $text_bubble = "<b>".$patient->name." - ".$tr->getGene->external_name()."-".$text_exons."</b>  Mutliplex : ".$primer->multiplex."  level:".$level;
				$text_bubble .= "; $itext";
				$itext="";
			
				if ($level == -2 ){
			#	if($sdm+$meanm>1.5 || $meanm-$sdm <0.5){
					$print_out.= $cgi->td({onClick=>$click,id=>"td_$id", onMouseOver=>"change_label('$tooltip_name','$text_bubble')", class=>"red",style=>"background-color:rgb(90,90,90);font-size:8px;"},$itext);
				}
				elsif ($level == -1) {
					
			#	elsif($sdp+$mean>1.5 || $mean-$sdp < 0.5){
					$print_out.= $cgi->td({onClick=>$click, onMouseOver=>"change_label('$tooltip_name','$text_bubble')",id=>"td_$id",class=>"red",style=>"background-color:rgb(120,140,140);font-size:8px;"},$itext);
				
				}
			#	elsif ($score - $sdp >= 1.4){
				elsif ($level == 2){
					my $color =int (200*($score-0.5))+50;
					
					$print_out.= $cgi->td({onClick=>$click, class=>"red", onMouseOver=>"change_label('$tooltip_name','$text_bubble')",id=>"td_$id",style=>"background-color:rgb(0,0,255);font-size:8px;"},$itext );
				}
				#elsif ($score + $sdp  <= 0.6){
				elsif($level == 1){
					
					my $color =int (50*(0.6-$score))+205;
					$print_out.= $cgi->td({onClick=>$click, class=>"red", onMouseOver=>"change_label('$tooltip_name','$text_bubble')",id=>"td_$id",style=>"background-color:rgb($color,0,0);font-size:8px;"},$itext);
				}
				elsif ($level==3){
					my $color =int (150*($score-1))+105;
					$print_out.= $cgi->td({onClick=>$click, onMouseOver=>"change_label('$tooltip_name','$text_bubble')",id=>"td_$id",style=>"background-color:rgb(0,50,$color);font-size:8px;"},$itext);
					#$print_out.= $cgi->td({onClick=>$click, class=>"red", onMouseOver=>"change_label('$tooltip_name','$text_bubble')",id=>"td_$id",style=>"background-color:rgb(0,0,0);font-size:8px;"},$itext);
				}
				elsif ($level==4){
					my $color =int (50*(0.6-$score))+205;
					$print_out.= $cgi->td({onClick=>$click, class=>"red", onMouseOver=>"change_label('$tooltip_name','$text_bubble')",id=>"td_$id",style=>"background-color:rgb($color,50,0);font-size:8px;"},$itext);
					#$print_out.= $cgi->td({onClick=>$click, class=>"red", onMouseOver=>"change_label('$tooltip_name','$text_bubble')",id=>"td_$id",style=>"background-color:rgb(0,0,0);font-size:8px;"},$itext);
				}
				
				else {
						$print_out.= $cgi->td({onClick=>$click, class=>"red", onMouseOver=>"change_label('$tooltip_name','$text_bubble')",id=>"td_$id",style=>"background-color:rgb(220,220,220);font-size:8px;"},$itext );
				}
			
			push(@$tab_ids,"'td_$id'");
	}#end primer
	$print_out.= $cgi->end_Tr();
#	
#
#	my $string_connect = join(",",@$tab_ids);
#	
#	$print_out.= qq{
#	<div id="$tooltip_name" data-dojo-type="dijit/Tooltip" data-dojo-props="connectId:[$string_connect],position:['above']">
#	</div>
#	};
	return $print_out;
}
1;
