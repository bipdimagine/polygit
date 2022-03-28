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
use File::stat;
#use Time::Piece;

use image_coverage;

my $w = "500px";




my $buffer = GBuffer->new();

my $cgi          = new CGI();
my $project_name = $cgi->param('project');
my $user = $cgi->param('user_name');
my $hgmd = $buffer->queryHgmd()->getHGMD($user);
#my $project = $buffer->newProject(-name=>$project_name);
my $project = $buffer->newProjectCache( -name 			=> $project_name, -typeFilters=>'individual' ,-cgi_object=>1);

my $ss = $project->getStructuralVariations();
my $nbc;
my $nbcl;
my $nbdm;
foreach my $v (@$ss){
	my $find;
	
	if ($v->text_clinvar =~/pathogenic/i){
		$find->{clinvar} ++;
		
	}
		if ($v->isDM ){
			$find->{dm} ++;
	}
		if ($v->clinical_local ){
			$find->{local} ++;
	}
	if ($find){
		foreach my $pn (@{$v->project->getPatientsAndControl}){
				$nbcl->{$pn->name} ++ if exists $find->{local};
				$nbdm->{$pn->name} ++  if exists $find->{dm};
				$nbc->{$pn->name} ++ if exists $find->{clinvar};
		}
	}
	
}

my $captures = $project->getCaptures();
my $similar = $project->similarProjects();

my $vquery = validationQuery->new(dbh=>$buffer->dbh,capture_name=>$project->validation_db());

my $patient_name = $cgi->param('patients');
#my $edit_mode = $cgi->param('edit_mode');
$patient_name ="all" ;#unless $patient_name;
my $patients = $project->get_list_patients($patient_name,",");;
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

	html::print_header_polydiag($cgi);
my $hsex = {
				2=>qq{  <i class="fa fa-venus" > </i>},
				1=>qq{  <i class="fa fa-mars" > </i>} ,
				'-1'=>qq{  <i class="fa fa-minus" > </i>} 
};
	
#	print $cgi->start_ol({class=>"circle-list"});
 my $nb_run = scalar (@{$project->getRuns()});
 my $nb_r =0;
 my $hash_dup;
 my $htr ;
 	map{$htr->{$_} ++} @{$project->bundle_transcripts() } ;
 foreach my $patient (@{$project->getPatientsAndControl}){

 my $filebed =  $patient->project->getVariationsDir("duplicate_region_calling")."/regions/".$patient->name().".dup.bed";
#system("mkdir $dir && chmod a+rwx $dir" ) unless -e $dir;
	
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
 foreach my $cn (keys %$hash_dup){
 	my $chr = $project->getChromosome($cn);
 	my $iter = $hash_dup->{$cn}->iterate_runs();
    while (my ( $from, $to ) = $iter->()) {
    	
 	my $ts  = $chr->getTranscriptsByPosition($from,$to);
 	foreach my $t (@$ts){
 		next unless exists $htr->{$t->id};
 		foreach my $e (@{$t->getExons}){
 			next if exists $he->{$e->name};
 			
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
 				 $he->{$e->name} ++;
 				 $ws->{$t->getGene->external_name()}->{tr}->{$t->name}->{$e->name}->{pos} = join(" ",@st);
 				 $ws->{$t->getGene->external_name()}->{tr}->{$t->name}->{$e->name}->{percent} =$p."%";
 				  $ws->{$t->getGene->external_name()}->{row} ++;
 				  
 			 	# my $warning = $t->getGene->external_name()." ".$t->name." ".$e->name." ".$span->as_string." ".$e->getGenomicSpan->as_string;
 			}
 			  
 			}
 		}
    }
 }
print qq{<div>};

if (keys %{$ws}){
print qq{<div class="btn  btn-info btn-xs btn-danger" style="position:relative;bottom:1px;min-width:200px;" onClick='collapse("projects_panel","projects_label")'>  <span id= "projects_label" class="glyphicon glyphicon-triangle-right  " aria-hidden="true"  style="float:left;"></span> <i class="fa fa-bell" aria-hidden="true"></i> Regions Dups</div>};
	
print $cgi->start_div({class=>"panel-body panel-collapse collapse panel-alert",style=>"font-size: 09px;font-family:  Verdana;",id=>"projects_panel"});	
	
print $cgi->start_table({class=>"table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 10px;font-family:  Verdana;", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
print qq{<thead>};
		print $cgi->start_Tr({class=>"danger"});
		print $cgi->th( {colspan=>"5",style=>"text-align: center;"} ,"Region Dups");
	print $cgi->end_Tr();
		print $cgi->start_Tr({class=>"success"});
		print $cgi->th( {style=>"text-align: center;"} ,"Gene");
		print $cgi->th( {style=>"text-align: center;"} ,"Transcript");
		print $cgi->th( {style=>"text-align: center;"} ,"exon");
		print $cgi->th( {style=>"text-align: center;"} ,"% dup");
		print $cgi->th( {style=>"text-align: center;"} ,"position");
	print $cgi->end_Tr();
		
foreach my $g (keys %{$ws}){
	print $cgi->start_Tr();
	print $cgi->td({class=>"warning",rowspan=>$ws->{$g}->{row}},$g);
	foreach my $t (keys %{$ws->{$g}->{tr}}){
		my $rs = scalar (keys %{$ws->{$g}->{tr}->{$t}});
		print $cgi->td({class=>"warning",rowspan=>$rs}, $t);
		foreach my $e (keys %{$ws->{$g}->{tr}->{$t}}){
				print $cgi->td({class=>"warning"}, $e);
				print $cgi->td({class=>"warning"}, $ws->{$g}->{tr}->{$t}->{$e}->{percent});
				print $cgi->td({class=>"warning"}, $ws->{$g}->{tr}->{$t}->{$e}->{pos});
			
					print $cgi->end_Tr();
		}
		print $cgi->end_Tr();
		
	}
	print $cgi->end_Tr();

}

print $cgi->end_table();
 print $cgi->end_div();
  print $cgi->end_div();
}
 
	foreach my $run (@{$project->getRuns()}){
		print $cgi->start_div({class=>"panel panel-primary",style=>"font-size: 11px;font-family:  Verdana;"});
		#print $cgi->start_li();
		my ($year,$month,$day) =  return_date($run->date);
		#warn return_date($run->date);
		$nb_r++;
		
		print $cgi->div({class=>"panel-heading"},'<i class="fa fa-flag fa-lg" style="text-align: center;font-color: white;"></i> &nbsp;'."run: $nb_r /$nb_run &nbsp;&nbsp;&nbsp;&nbsp;<i class='fa fa-calendar'></i> &nbsp; $day-$month-$year&nbsp;&nbsp;&nbsp;<i class='fa  fa-desktop'></i> &nbsp; ".$run->machine);
		print "<br>";
		
#		print $cgi->p("$day-$month-$year ".$run->machine);
		#print $cgi->start_table({class=>"sum"});
		print $cgi->start_table({class=>"table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 10px;font-family:  Verdana;", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
		my @title2 = (" ","Variant Validation","validation date","responsable");
		print qq{<thead>};
		print $cgi->start_Tr({class=>"success"});

#		if ($project->isFamilial()){
#				print $cgi->th( {colspan=>"4",style=>"text-align: center;"} ,"General") ;
#		}
#		else {
#			
#		}
	 my $col_hgmd =3;
	 $col_hgmd = 2 unless $hgmd ==1;

		print $cgi->th( {colspan=>"6",style=>"text-align: center;"} ,"General") ;

		print $cgi->th( {colspan=>"2",style=>"text-align: center;"} ,"Coverage");
		print $cgi->th( {colspan=>"3",style=>"text-align: center;"} ,"Clinical");
		print $cgi->th( {colspan=>"2",style=>"text-align: center;"} ,"Calling");
		print $cgi->th({colspan=>$col_hgmd,style=>"text-align: center"},"Variations Validation ");
		print $cgi->th({colspan=>2,style=>"text-align: center"},"Exons Validation ");
		print $cgi->th({colspan=>2,style=>"text-align: center"},"Report");
		print $cgi->th({colspan=>3,style=>"text-align: center"},"Latest modification");
		print $cgi->end_Tr();
		print qq{</thead>};
		print $cgi->start_Tr({class=>"warning"});
		my @title = ("Patient","Sex","SRY","Cov","30x","clinical","Sub","Indels","Sanger","Ngs","ToDo", "Sanger","ToDo","View","Print","Analyse","User","Validation");
		 @title = ("Fam","Patient","Sex","SRY","status","Cov","30x","hgmd","clinvar","local","Sub","Indels","Sanger","Ngs","ToDo", "Sanger","ToDo","View","Print","Analyse","User","Validation");# if ($project->isFamilial());
		 @title = ("Fam","Patient","Sex","SRY","status","Cov","30x","clinvar","local","Sub","Indels","Sanger","Ngs","ToDo", "Sanger","ToDo","View","Print","Analyse","User","Validation") unless $hgmd == 1;

		#	print $cgi->th();
		print $cgi->td({style=>"text-align: center;"},qq{<input id="check_all" type="checkbox" aria-label="..."  onchange="select_all(this)"></input>});
		foreach my $p (@title){
			print $cgi->td({style=>"text-align: center;"},$p);
		}
	print $cgi->end_Tr();
	
	
		#TODO: here
		my @lPatients;
		foreach my $p (@{$run->project->getPatientsAndControl()}) {
			push (@lPatients, $p) if ($p->is_control());
		}
		foreach my $p (sort{$a->getFamily->name cmp $b->getFamily->name or $a->name cmp $b->name} @{$run->getPatients()}){
			push (@lPatients, $p);
		}
			
		foreach my $p (@lPatients){
			my $pname = "check_".$p->name();
			if ($p->is_control) { print $cgi->start_Tr({style=>"background-color:#d1e5ff"}); }
			else { print $cgi->start_Tr(); }
			print $cgi->td(qq{<input id="$pname" type="checkbox" aria-label="..." onClick="selection(event,this)"></input>});
			print $cgi->td($p->getFamily->name);
			print $cgi->td($p->name);
			my $cov = $p->coverage();
			
			my $class ={};
			unless ($p->is_control()) {
				my $sex_eval = $p->compute_sex(); 
				my $sex = $p->sex(); 
				if ($sex_eval ne $sex && $sex_eval ne -1){
					 $class ={class=>"danger"};
				} 
				my $cov_sry = $p->coverage_SRY();
				print $cgi->td($class,$hsex->{$sex});
				print $cgi->td($class,$hsex->{$sex_eval}."<small>(".$cov_sry.")</small>");
			}
			else {
				print $cgi->td($class,"-");
				print $cgi->td($class,"-");
			}
		
			my $color = "black";
			unless ($p->is_control()) {
				$color = "red" if $p->isIll;
			}
			if ($p->is_control()) {
				print $cgi->td($class,qq{<span>Control</span>});
			}				
			elsif ($p->isChild){
				print $cgi->td($class,qq{<i class="fa fa-child fa-1x" aria-hidden="true" style="color:$color"></i>});
			}
			elsif ($p->isMother){
				print $cgi->td($class,qq{<i class="fa fa-female fa-1x" aria-hidden="true" style="color:$color"></i>});
			}
			elsif ($p->isFather){
				print $cgi->td($class,qq{<i class="fa fa-male fa-1x" aria-hidden="true" style="color:$color"></i>});
			}
			print $cgi->td($cov->{mean});
			
		
			#	warn (stat $filename)[9]; 
				#warn Dumper($cov);
				#$cov->{'30x'} = $coverage->{$p->id}->{30}/$coverage->{$p->id}->{nb}*100;
				
			print $cgi->td($cov->{'30x'});
			
			if ($p->is_control()) {
				print $cgi->td('-');
				print $cgi->td('-');
				print $cgi->td('-');
				print $cgi->td('-');
				print $cgi->td('-');
				print $cgi->td('-');
				print $cgi->td('-');
				print $cgi->td('-');
				print $cgi->td('-');
				print $cgi->td('-');
				print $cgi->td('-');
				print $cgi->td('-');
				print $cgi->td('-');
				print $cgi->td('-');
				print $cgi->td('-');
			}
			else {
				my $style ={};
				 $style = {style =>'font-size: 10px;font-family:  Arial;text-align:center;  vertical-align:middle;background-color: coral;' } if $nbdm->{$p->name}>= 1 ;
				print $cgi->td($style,$nbdm->{$p->name}) if $hgmd ==1;
				print $cgi->td($nbc->{$p->name});
				print $cgi->td($nbcl->{$p->name});
				my $nbi =0;
				my $f1s = $p->getIndelsFiles();
					my $nbv = "0";
				foreach my $f1 (@$f1s){
				#zcat IMB.vcf.gz IMD4.vcf.gz | grep -v "#" | cut -f 1,2 | sort -u
				
				my $nbvt += `zgrep -vc "#" $f1` if -e $f1;
				chomp($nbvt);
				$nbi += $nbvt;
				}
				#my $vs = $p->getVariations();
				#warn $p->name();
				if ($nbi == 0){
				
					my $bcftools = $buffer->software("bcftools");
					$bcftools ="bcftools" unless -e $bcftools; 
					 my $f1s = $p->getVariationsFiles();
					  $nbv = $p->countSubstitutions();
					 $nbi =$p->countIndels();
					# my $files = join(" ",@$f1s);
					 
					 #my $cmd = 	qq{$bcftools  concat $files -a -d all  2>/dev/null| $bcftools stats - | grep "^SN" | grep "SNPs\\|indels" | cut -f 4};
					 #($nbv,$nbi) = `$cmd`;
				}
				else {
				my $f1s = $p->getVariationsFiles();
					my $nbv = "0";
				foreach my $f1 (@$f1s){
				#zcat IMB.vcf.gz IMD4.vcf.gz | grep -v "#" | cut -f 1,2 | sort -u
				
				my $nbvt += `zgrep -vc "#" $f1` if -e $f1;
				chomp($nbvt);
				$nbv += $nbvt;
				
				}
				}
				print $cgi->td($nbv);
				
	
				
				
				print $cgi->td($nbi);
			
				
				#patient validated ?
				my $text_vs = "No Db";
				eval{
				my $variations_sanger = $vquery->get_variations_sanger(project_name=>$project->name,sample_name=>$p->name);
				my $text ="";
				my $color = {style =>'font-size: 10px;font-family:  Arial;text-align:center;  vertical-align:middle;' };
				 $text = "<span class='badge alert-success' >".scalar(keys %{$variations_sanger}).'</span>' ;
				print $cgi->td($color,$text);
				my $variations_validated = $vquery->get_variations_ions(project_name=>$project->name,sample_name=>$p->name);
				$text = "<span class='badge alert-info' >".scalar(keys %{$variations_validated}).'</span>' ;
				#$color = "black";
				print $cgi->td($color, $text);
				my $variations_todo = $vquery->get_variations_todo(project_name=>$project->name,sample_name=>$p->name,uniq=>1);
				my %nb_todo;
				
			
				$text = "<span class='badge alert-warning' >".scalar(keys %{$variations_todo}).'</span>' ;
				my %users;
				print $cgi->td($color,$text);
				
				$text ="-";
				
	
					my $exons_todo = $vquery->get_exons(project_name=>$project_name,sample_name=>$p->name);
					
					$text = "<span class='badge alert-success' >".scalar( grep {$_->{done} ==1} values %$exons_todo).'</span>' ;
					print $cgi->td($color,$text);
					$text = "<span class='badge alert-warning' >".scalar( grep {$_->{done} ==0} values %$exons_todo).'</span>' ;
					print $cgi->td($color,$text);
				
	#			my $test = $vquery->get_report(project=>$project->name,sample=>$p->name);
	#			
	#			
	#			my @validations = ("-","-");
	#			
	#			if ($test){
	#				@validations= ();
	#				push(@validations,$test->{creation_date});
	#				push(@validations,$test->{user_name});
	#			}
	#			
	#			foreach my $va (@validations){
	#				print $cgi->td($va);
	#			}
				};
				if ($@){
						print $cgi->td("-");
						print $cgi->td("-");
						print $cgi->td("-");
						print $cgi->td("-");
						print $cgi->td("-");
				}
				my $pp = $p->name;
				my $cmd = qq{printer('$pp');};
				my $cmd2 = qq{printer2('$pp','1');};
				my $td =[
				'<a type="button" class="btn btn-xs btn-info" onclick="'.$cmd.'"><i class="fa fa-clipboard pull-left  "></i>View</button></a>',
				'<a type="button" class="btn btn-xs btn-success" onclick="'.$cmd2.'"><i class="fa fa-print pull-left  "></i>Print</button></a>'
				
				];
				
				
				print $cgi->td($td);
			#	my $mtime = stat("/etc/passwd")->mtime;
				eval{
				my ($users,$last_user,$date) = $vquery->get_validations_users(project_name=>$project->name,sample_name=>$p->name);
				my $st_date ="-";
				my $text ="-";
					if ($date){
	
							$text=$last_user;
							my @d = split(" ",$date);
							
							 $st_date = join("-",return_date($date));
					}
					else {
						$text = "-";
					}
					print $cgi->td(join("-",($year,$month,$day)));
					print $cgi->td($text);
					print $cgi->td($st_date);
				};
	#			if ($@){
	#				 print $cgi->td("-");
	#				  print $cgi->td("-");
	#				 print $cgi->td(localtime((stat  $p->getVariationsFile())[9])->ymd);
	#		
	#			}
			}
			print $cgi->end_Tr();
		}
		
		print $cgi->end_table();
		print $cgi->div({class=>"panel-footer"});
		print $cgi->end_div();
	}

#	print $cgi->end_ol;
print qq{
	 <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.2/jquery.min.js"></script>
    <!-- Include all compiled plugins (below), or include individual files as needed -->
 
  
};

	sub return_date {
		my ($dd) = @_;
		 my @amonths = ('Jan', 'Feb', 'Mar', 'Apr','May',"Jun","Jul","Aug","Sep","Oct","Nov","Dec");
		my ($date,$time) = split(" ",$dd);
	    my ($year,$month,$day) =  split("-",$date);
		return ($year,$amonths[$month-1],$day);
	}