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
use List::MoreUtils qw{ natatime };
use image_coverage;
use warnings;
use utf8;
use QueryValidationAcmg;
use feature ':5.10';
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";
my $w = "500px";




my $buffer = GBuffer->new();


my $cgi          = new CGI();
html::print_header_polydiag($cgi);
my $project_name = $cgi->param('project');
my $user = $cgi->param('user_name');

my $hgmd = $buffer->queryHgmd->getHGMD($user);
#my $project = $buffer->newProject(-name=>$project_name);
my $project = $buffer->newProjectCache( -name 			=> $project_name, -typeFilters=>'individual' ,-cgi_object=>1);


#die();
my $panel_name = $cgi->param('panel');
#my @transcripts_cgi ;
my $panel;
if ($panel_name){
 $panel = $project->getPanel("$panel_name");

}
else {
	confess();
	# @transcripts_cgi = @{$project->bundle_transcripts() } ;
}
my $t = time;


 $t = time;
my $hv;
my $stat; 


 $t = time;
 $|=1;
print qq{<div style="display: none">};
my $no = $project->noSqlQuality("r");
my $data = $no->get($project->name,"mendelian");
my $hmendel;
foreach my $line (@{$data->{data}}){
	my $name = $line->{sample}->{text};
	my $value =  $line->{mendelian}->{text};
	$value =~ s/\[//;
		$value =~ s/\]//;
	my ($t,$p,$n) = split(" ",$value);
	$p =~ s/%//;
	$hmendel->{$name} = $p;

#	$hmendel->{$name} = $line->{mendelian}->{text};
	
}



	
my $t1 =time;

my $pm = new Parallel::ForkManager(15);





$project->getChromosomes();
$project->getPatients();

 $project->disconnect();
$pm->run_on_finish(
    sub { 
    	my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$h)=@_;
  
    	unless (defined($h) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
		
		foreach my $t (keys %$h){
			foreach my $p (keys %{$h->{$t}}){
				$stat->{project}->{$t}->{$p} += $h->{$t}->{$p} ;
			}
		}
	
		
    }
    );
    
#warn Dumper $data;
foreach my $chr (@{$project->getChromosomes}){
	next if ($chr->not_used());
my $pid = $pm->start and next;
foreach my $p (@{$project->getPatients}){


#	next;
	print ". ===>".$chr->name;
	
#	warn $chr->lmdb_score_impact->get("gnomad_ho_ac_5");	

	my $vo = $p->getVectorOrigin($chr);
		
		print "*";
#	my $v =  $vo & $panel->getVectorDM($chr); 
#	warn  $chr->countThisVariants($panel->getVectorDM($chr))." ".$chr->countThisVariants($v)." ".$chr->countThisVariants($chr->vectorDM);
	#my $v2 = $vo & $panel->getVectorClinvarPathogenic($chr); 
	#my $v3 = $vo & $panel->getVector($chr); 
	my $v4 = $p->getVectorDeletions($chr) | $p->getVectorInsertions($chr);
#
	#warn $p->name." ".$chr->countThisVariants($v);
	#$stat->{dm}->{$p->name} +=  $chr->countThisVariants($v);
	#$stat->{clinvar}->{$p->name} +=  $chr->countThisVariants($v2);
	#$stat->{total}->{$p->name} +=  $chr->countThisVariants($v3);

	$stat->{dm}->{$p->name} +=  $p->countThisVariants($p->getVectorDM($chr));
	$stat->{total}->{$p->name} += $p->countThisVariants($p->getVectorOrigin($chr));
	$stat->{total_indels}->{$p->name} += $p->countThisVariants($v4);
	$stat->{total_substitutions}->{$p->name} += $p->countThisVariants($p->getVectorSubstitutions($chr));
	$stat->{total_he}->{$p->name} += $p->countThisVariants($p->getVectorHe($chr));
	
	
	$stat->{clinvar}->{$p->name} +=  $p->countThisVariants($p->getVectorClinvar($chr));
	print ".clinvar";
}
$pm->finish(0,$stat);
print "end";
}
$pm->wait_all_children();
$project->buffer->dbh_reconnect();
#warn abs(time-$t1);
print qq{</div>};
	
#die();
my $nbc;
my $nbcl;
my $nbdm;
my $nb_variations;


	print qq{<div style="visibility: hidden">};
	
	my $capture = $project->getCaptures()->[0];
my $query = $project->buffer->getQuery();


#my $stats = compute_hgmd_clinvar(\@transcripts_cgi);
	print qq{</div>};


my $captures = $project->getCaptures();
my $vquery = QueryValidationAcmg->new(dbh=>$buffer->dbh,database=>$project->validation_db());

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

	
my $hsex = {
				2=>qq{  <i class="fa fa-venus" > </i>},
				1=>qq{  <i class="fa fa-mars" > </i>} ,
				'-1'=>qq{  <i class="fa fa-minus" > </i>} 
};
	
#	print $cgi->start_ol({class=>"circle-list"});
 my $nb_run = scalar (@{$project->getRuns()});
 my $nb_r =0;

my $fam = $project->getFamilies();
print $cgi->start_div({class=>"panel panel-primary",style=>"font-size: 11px;font-family:  Verdana;"});

		 ($year,$month,$day)  =("-","-","-");# =return_date($run->date);
		my @types = ("mendelian","identity","coverage_stats","coverage_transcripts","variations");
		 $no = $project->noSqlQuality("r");
		my $nlevel = 0;
		my $qb="";
		foreach my $type (@types) {
		my $data = $no->get($project->name,$type);
		my $level =  check_level($data);
		my $h = "#".$type;
		$qb .= qq{<div class="detail-section" data-toggle="modal" data-target="$h" ><a type="button" class="btn  btn-$level btn-xs" onclick="open_report_quality()" style="float: right;"><i class="fa fa-clipboard"></i>&nbsp$type</button></a></div>};
		}
		print "<br>";
			print $cgi->div({class=>"panel-heading"},'<i class="fa fa-flag fa-lg" style="text-align: center;font-color: white;"></i> &nbsp;'."run: $nb_r /$nb_run &nbsp;&nbsp;&nbsp;&nbsp;<i class='fa fa-calendar'></i> &nbsp; $day-$month-$year&nbsp;&nbsp;&nbsp;<i class='fa  fa-desktop'></i> &nbsp; "."&nbsp".$qb);
	
			print $cgi->start_table({class=>"table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 9px;font-family:  Verdana;", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
		
		#### print header table
		 my $col_hgmd =3;
	 $col_hgmd = 2 unless $hgmd ==1;

		print $cgi->th( {colspan=>"9",style=>"text-align: center;"} ,"General") ;

		print $cgi->th( {colspan=>"2",style=>"text-align: center;"} ,"Coverage");
		print $cgi->th( {colspan=>"2",style=>"text-align: center;"} ,"Clinical");
		print $cgi->th( {colspan=>"3",style=>"text-align: center;"} ,"Calling");
		print $cgi->th({colspan=>4,style=>"text-align: center;"},"Validation");
		print $cgi->end_Tr();
		print qq{</thead>};

		my  @title = ("-","Fam","SV","Patient","Sex","status","SRY","mendelian error","Cov","30x","hgmd","clinvar","Sub","Indels","%He","validation");# if ($project->isFamilial());
		 @title = ("-","Fam","SV","Patient","Sex","status","SRY", "mendelian error","status","Cov","30x","clinvar","Sub","Indels","%He","validation") unless $hgmd == 1;
		 
		
		print $cgi->start_Tr({class=>"warning"});
		print $cgi->td({style=>"text-align: center;"},qq{<input id="check_all" type="checkbox" aria-label="..."  onchange="select_all(this)"></input>});
		foreach my $p (@title){
			print $cgi->td({style=>"text-align: center;"},$p);
		}
		print $cgi->end_Tr();
		
		#end header
		
		my @colors = ("#F4F4F4","#DEDFDE");
		@colors = ("#F9F6FF","#FFFFFF");
		my $nb =0;
	
		
		foreach my $fam (sort{$a->name cmp $b->name }@{$project->getFamilies}){
			
		my $color = $colors[$nb%2]; 
		$nb ++;
  		my $nb_members = scalar(@{$fam->getMembers});
  		foreach my $p1 (@{$fam->getMembers}){
  			my $hval = $vquery->getValidationPatient($p1);
  			
  			 $nb_members += scalar(keys %$hval);
  		}
  		
  		print $cgi->start_Tr({style =>"background-color:$color"});
  		my $pname = "check_".$fam->name();
  		 if ($nb_members>1){
  		 		print $cgi->td({rowspan=>$nb_members,"style"=>"vertical-align:middle"},qq{<input id="$pname" type="checkbox" aria-label="..." onClick="selection(event,this)"></input>});
  		 		print $cgi->td({rowspan=>$nb_members,"style"=>"vertical-align:middle"},qq{<img src="https://img.icons8.com/office/24/000000/family.png">} );
  		 		print $cgi->td({rowspan=>$nb_members,"style"=>"vertical-align:middle"},$fam->name);
  		 }
  		 else {
  		 		print $cgi->td({rowspan=>$nb_members,"style"=>"vertical-align:middle"},qq{<input id="$pname" type="checkbox" aria-label="..." onClick="selection(event,this)"></input>});
  		 		print $cgi->td({rowspan=>$nb_members,"style"=>"vertical-align:middle"},"-"); 
  		 		print $cgi->td({rowspan=>$nb_members,"style"=>"vertical-align:middle"},"-"); 
  		 }

  		
		foreach my $p (@{$fam->getMembers}){
			print_line_patient($p,0);
			print $cgi->start_Tr({style =>"background-color:$color"});
		}
		print $cgi->end_Tr;
	}
	print $cgi->end_table();
print $cgi->end_div();
exit(0);


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
 	 my $fsize = "font-size:10px";
 		my $pname = "check_".$p->name();
			my $class ={};
			my $hval = $vquery->getValidationPatient($p);
			if (keys %$hval){
			 $nb_row += scalar(keys %$hval);
			}
			$class = {rowspan=>$nb_row} if $nb_row >1;
			$class->{style}= "vertical-align:middle";
			my $fam = $p->getFamily();
			if ($fam->isTrio && $p->isChild){
					my $url = "http://defidiag.polyweb.fr/polyweb/html/manta/Url.Trio_SV_Editor.html?project=".$project->name."&filename=".$p->name;
					print $cgi->td($class,qq{<a type="button" class= "btn btn-xs  btn-primary" href="$url" role="button" target="_blank"><img src="https://img.icons8.com/color/30/000000/biotech.png">&nbspPolyCyto</a>});
			}
			elsif ($fam->isTrio && !($p->isChild)){
					#my $url = "http://defidiag.polyweb.fr/polyweb/html/manta/Url.Trio_SV_Editor.html?project=".$project->name."&filename=".$p->name;
					print $cgi->td($class,qq{.});
			}
			else {
				my $url = "http://defidiag.polyweb.fr/polyweb/html/manta/Url.Basic_SV_Editor.html?project=".$project->name."&filename=".$p->name;
					print $cgi->td($class,qq{<a type="button" class= "btn btn-xs  btn-primary" href="$url" role="button" target="_blank"><img src="https://img.icons8.com/color/30/000000/biotech.png">&nbspPolyCyto</a>});
				
			}
			
			print $cgi->td($class,$p->name);
		
			my $cov = $p->coverage();
			warn $cov;
			#global value
			my $sex_eval = $p->compute_sex(); 
			my $sex = $p->sex(); 
			if ($sex_eval ne $sex && $sex_eval ne -1){
				 $class->{class}= "danger";
			} 
			my $cov_sry = $p->coverage_SRY();
			
				if ($p->isChild){
			#	print $cgi->td($class,qq{<i class="fa fa-child fa-1x" aria-hidden="true" style="color:$color"></i>});
				print $cgi->td($class,qq{<img src="https://img.icons8.com/color/24/000000/boy.png">}) if $p->isMale;
				print $cgi->td($class,qq{<img src="https://img.icons8.com/color/24/000000/girl.png">}) if $p->isFemale;
				
			}
			elsif ($p->isMother){
				#print $cgi->td($class,qq{<i class="fa fa-female fa-1x" aria-hidden="true" style="color:$color"></i>});
				print $cgi->td($class,qq{<img src="https://img.icons8.com/office/16/000000/businesswoman.png">});
			}
			elsif ($p->isFather){
				
				#print $cgi->td($class,qq{<i class="fa fa-male fa-1x" aria-hidden="true" style="color:$color"></i>});
				print $cgi->td($class,qq{<img src="https://img.icons8.com/office/16/000000/person-male.png">});
			}
			
			my $color = "black";
			$color = "red" if $p->isIll;
			if ($p->isIll){
				print $cgi->td($class,qq{<img src="https://img.icons8.com/office/24/000000/treatment-plan.png">});
			}
			else {
				print $cgi->td($class,qq{<img src="https://img.icons8.com/office/16/000000/checked.png">});
			}
			#print $cgi->td($class,$hsex->{$sex});
			 $sex = $p->sex(); 
			if ($sex_eval ne $sex && $sex_eval ne -1){
				 $class->{class}= "danger";
			} 
			
			print $cgi->td($class,$hsex->{$sex_eval}."<small>(".$cov_sry.")</small>");
			
			my $v =$hmendel->{$p->name};
			
			my $btn_class = qq{class= "btn btn-xs btn-success "  style = "$fsize" };
			if ($p->isParent()){
				my $c = $p->getFamily->getChildren();
				 if ($v) {
				 	$v = $v/ scalar(@$c);
				 }
			}
			unless ($v){
				 $btn_class = qq{class= "btn btn-xs btn-light " style = "$fsize"  };
				 $v= "-";
			}
			else {
				$btn_class = qq{class= "btn  btn-xs btn-warning " style = "$fsize" } if $v >1;
				$btn_class = qq{class= "btn  btn-xs btn-alert " style = "$fsize" } if $v > 3;
				
			}
			print $cgi->td($class,qq{<button type="button" $btn_class >$v <span>&#37;</span> </button>});
			
			
			 $v = $cov->{mean};
			 $v= 0 unless $v;
			 $btn_class = qq{class= "btn btn-xs btn-success " style = "$fsize" };
			$btn_class = qq{class= "btn  btn-xs btn-warning " style = "$fsize" } if $v < 15;
			$btn_class = qq{class= "btn  btn-xs btn-alert " style = "$fsize" } if $v < 10;
			print $cgi->td($class,qq{<button type="button" $btn_class >$v</button>});
			$v = $cov->{'30x'};
			 $v= 0 unless $v;
			print $cgi->td($class,qq{<button type="button" $btn_class >$v %</button>});
		
			my $style ={};
		 #	$style = {style =>'font-size: 10px;font-family:  Arial;text-align:center;  vertical-align:middle;background-color:#82CFFD;'} if $stat->{project}->{dm}->{$p->name}  >= 1 && $p->isFather;
		 #		$style = {style =>'font-size: 10px;font-family:  Arial;text-align:center;  vertical-align:middle;background-color:#F6CCDA;'} if $stat->{project}->{dm}->{$p->name}  >= 1 && $p->isMother;
			 $style = {style =>'font-size: 10px;font-family:  Arial;text-align:center;  vertical-align:middle;background-color: coral;'} if $stat->{project}->{dm}->{$p->name}  >= 1 && $p->isChild();
			
			 $v = $stat->{project}->{dm}->{$p->name};
			
			 $btn_class = qq{class= "btn btn-xs btn-primary " style="background-color: #D0D0D0;$fsize"};
			$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: coral;$fsize"}  if $v > 0  && $p->isChild() ;
			print $cgi->td($class,qq{<button type="button" $btn_class >$v</button>}) if $hgmd ==1;	
			 $v = $stat->{project}->{clinvar}->{$p->name};
			$btn_class = qq{class= "btn btn-xs btn-primary " style="background-color: #D0D0D0;$fsize"};
			$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: coral;$fsize"}  if $v > 0  && $p->isChild() ;
			print $cgi->td($class,qq{<button type="button" $btn_class >$v</button>});
			
			$btn_class = qq{class= "btn btn-xs btn-primary " style="background-color: #9BC8A5;$fsize"};
			$v = $stat->{project}->{total_substitutions}->{$p->name};
			print $cgi->td($class,qq{<button type="button" $btn_class >$v</button>});
			
			$v =$stat->{project}->{total_indels}->{$p->name};
			print $cgi->td($class,qq{<button type="button" $btn_class >$v</button>});
			$v = int(($stat->{project}->{total_he}->{$p->name}/$stat->{project}->{total}->{$p->name}) *1000)/10;
			print $cgi->td($class,qq{<button type="button" $btn_class >$v %</button>});
		
			 
			my $text_vs = "No Db";

			#		print $cgi->td("-");
			#		print $cgi->td("-");
			#		print $cgi->td("-");
			#		print $cgi->td("-");
			#		print $cgi->td("-");
			my $pp = $p->name;
			my $cmd = qq{printer('$pp');};
			my $cmd2 = qq{printer2('$pp','1');};
			my $td =[
			'<a type="button" class="btn btn-xs btn-info" onclick="'.$cmd.'"><i class="fa fa-clipboard pull-left  "></i>View</button></a>',
			'<a type="button" class="btn btn-xs btn-success" onclick="'.$cmd2.'"><i class="fa fa-print pull-left  "></i>Print</button></a>'
			
			];
			
			#print $cgi->td($class,$td);
		#	my $mtime = stat("/etc/passwd")->mtime; 
			
			if (keys %$hval){
				
				my @headers_validations = ("igv_web","var_name","trio","gene","table_transcript","nomenclature");
				my @header_transcripts = ("consequence","enst","nm");
			
				my $rowspan = scalar(keys %$hval);
			 	$class->{rowspan} -= $rowspan;
			 	$class->{rowspan} = 1 if $class->{rowspan} ==0;
			 	my $fam = $p->getFamily();
			 	my $fin = scalar (keys %{$hval}) ;
			 	my $pos =0;
			 	my $out;
			 	$out.= $cgi->start_table({class=>"table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 9px;font-family:  Verdana;margin-bottom: 0px;", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
			 	foreach my $k  (keys %{$hval}){
			 		
			 		my $val = $hval->{$k}->[0];
			 			my $color = "#9BC8A5";
			 		$color = "#E74C3C" if $val->{validation}== 5 ;
			 		$color = "coral" if $val->{validation}== 4 ;
			 		$color = "orange" if $val->{validation}== 3 ;
			 		
			 		$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: $color;$fsize"} ;
			 		my $v  = $project->_newVariant($val->{polyid});
			 		my $gene = $project->newGene($val->{gene_id});
			 		my $hvariation = update::construct_variant( $project, $v, $gene->getTranscripts->[0], $p, undef );
			 		my $gn = $val->{gene_name};
			 		 $hvariation->{gene} = qq{<button type="button" $btn_class >$gn</button>};
			 		update::trio($project,undef,$hvariation,$p,$cgi,undef);
			 		$hvariation->{table_transcript} = update::construct_table_transcript($v, $cgi,\@header_transcripts,3,$gene); 
		
			 		my $st_date = join("-",return_date($val->{modification_date}));
			 		my $term = $buffer->getValidationTerm($val->{validation});
			 		
			 		# if $v > 0  && $p->isChild() ;
					$out.= $cgi->td($class,qq{<button type="button" $btn_class >$term</button>});
					 my $pn = $val->{user_name};
			 		 $out.= $cgi->td($class,qq{<button type="button" $btn_class >$pn</button>});
			 		 $out.=  $cgi->td($class,qq{<button type="button" $btn_class >$st_date</button>});
			 		 
					foreach my $h (@headers_validations){
							$out.=  $cgi->td($class,$hvariation->{$h});
					}
					
			 		 
			 		$out.=  $cgi->end_Tr();
			 		#$pos++;
			 		#print $cgi->start_Tr() if $pos < ($fin-3); 
			 		
			 	}
			 	$out.= $cgi->end_table();
			 	
			 	print $cgi->td({style=>"background-color:#E74C3C;padding:3px"},$out);
			 
			}
			else {
				print  $cgi->td($class,"-");
				# print   $cgi->td($class,"-");
				#print  $cgi->td($class,"-");

			}
			
			#if ($project->isTrio){
				
			#}
		
			print $cgi->end_Tr();
 }

 

	