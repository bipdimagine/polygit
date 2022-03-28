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
 use List::MoreUtils qw{ natatime };
#use Set::;
use Storable qw/store thaw retrieve/;
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
use lib "$Bin/../packages/validation_variation"; 
use draw_cnv; 
use infos_coverage_exons;
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use preload_coverage;

use image_coverage;
$|=1;
my $buffer = GBuffer->new();

        
                
my $cgi          = new CGI();
my $out;
my $prg =  $cgi->url(-relative=>1);
my $url_image = url(-absolute=>1);
$url_image =~ s/$prg/image_coverage.pl/;


my $order = $cgi->param('order');

my $project_name = $cgi->param('project');
my $only_low;
 $only_low = 1 if $cgi->param('only_low') ==1;
my $project = $buffer->newProjectCache(-name=>$project_name);
my $panel_name = $cgi->param('panel_name');
$panel_name ="DI_DI44" unless  $panel_name;
#$project->setPanel("$panel") if $panel;
my $utr =$cgi->param('utr')+0;
my $intronic = $cgi->param('intronic')+0;

my $limit = $cgi->param('limit');
my $padding = $cgi->param('span');

my @transcripts_cgi;
my $gene_id =  $cgi->param('gene');
my $panel;
if ($gene_id) {
	my @real_gene_id;
	foreach my $gid (split(",",$gene_id)){
		my $ensg_id = $project->getEnsgIDs($gid);
		foreach my $geid (@$ensg_id){
			my $gene = $project->newGene($geid);
			
			push(@transcripts_cgi, map{$_->id} @{$gene->getTranscripts});
		}
	}
	
	
}
else {
	my $cgi_transcript =  $cgi->param('transcripts');
	if ($panel_name and $panel_name ne 'all'){
	 	$panel = $project->getPanel("$panel_name");
	 	@transcripts_cgi = @{$panel->main_transcripts() } ;
	}
	elsif ($cgi_transcript eq "all"){
		push(@transcripts_cgi, 'all');
	}
	else {
		@transcripts_cgi = split(",",$cgi_transcript);
	}
}

#my $capture = $project->getCaptures()->[0];
my $vquery;
if ($project->isDiagnostic){
	$vquery = validationQuery->new(dbh=>$buffer->dbh,capture_name=>$project->validation_db());
}
my $transcripts_exons_todo ={};
my $transcripts_exons_validated = {};

foreach my $patient (@{$project->getPatients}){
	my $exons_todo = {};
#	 $exons_todo = $vquery->get_exons(project_name=>$project_name,sample_name=>$patient->{name}) if $vquery;
	foreach my $v (values %$exons_todo){
	
		if  ($v->{done} == 0){
			$transcripts_exons_todo->{$v->{transcript}}  ++ ;
		}
		else {
				$transcripts_exons_validated->{$v->{transcript}}  ++ ;
		}
		
	}
}
my $patient_names = $cgi->param('patients');

my @transcripts;
if (scalar(@transcripts_cgi) == 1 and $transcripts_cgi[0] eq 'all') {
	@transcripts = @{$project->getTranscripts()};
}
else {
	@transcripts = @{$project->newTranscripts(\@transcripts_cgi)};	
}
 
if ($order eq "name") {
	 @transcripts = sort{$a->getGene->external_name cmp $b->getGene->external_name} @transcripts; 
}
else {
	 @transcripts = sort{$a->getChromosome->id <=> $b->getChromosome->id || $a->start <=> $b->start} @transcripts; 
}



my $nbgenes =0;
my $nbt = scalar(@transcripts);
my $nbtp = 0;#scalar(@transcripts);
#@transcripts = splice (@transcripts ,$cgi->param('start'),$cgi->param('step')) if $cgi->param('step') && !$only_low ; 
$out.= $cgi->start_table({class=>"table table-striped table-condensed table-bordered  table-mybordered",style=>"font-size: 8px;font-family:  Verdana;"});
my $col = 12;

$col = 10 if scalar(@{$project->getPatients()}>15);
$col = 6 if scalar(@{$project->getPatients()}>30);
$col = 5 if scalar(@{$project->getPatients()}>40);

my $nb=0;

$url_image.="?project=".$project_name."&limit=$limit&span=$padding&utr=$utr&intronic=$intronic&patients=$patient_names&transcript=";
my @ths;
my@tds;
html::print_cgi_header($cgi);
my $no =  $project->noSqlCoverage();

my $patients = $project->getPatients();
my $nb =0;
my @problem_transcripts;

 if ($only_low) {
 	#my $list_transcripts = getLowTranscripts();
 	
 	#@transcripts = grep {exists $list_transcripts->{$_->id}} @transcripts;
 }
my $images = uri_image($project,\@transcripts);
@transcripts = grep {exists $images->{$_->id}} @transcripts;


$out.= html::print_cadre($cgi," Coverage $panel_name  : Transcripts ($nbtp/$nbt) ");
print qq {
	<div class="btn-group  btn-sm clearfix pull-right" style="position:relative;top:-14px;float:left">
         	<button type="button" class="btn btn-alert btn-xs">Genes <span class="badge badge-success">$nbgenes </span></button>	
         	<button type="button" class="btn btn-alert btn-xs">Transcripts <span class="badge"> $nbt</span></button>
         	<button type="button" class="btn btn-danger btn-xs">Problem <span class="badge">$nbtp</span></button>
         </div>
};
for (my $i=0;$i<@transcripts;$i++){

if ($i%$col == 0 && $i>0){
		$out.=print_lines(\@ths,\@tds);
		@ths=();
		@tds=();
	
	}


 	 my $tname= $transcripts[$i]->name;
# 	my $class = "nobg";
# $class ="th1" if $z%2 ==0;
 my $text = $transcripts[$i]->getGene()->external_name()."<br>".$tname;
 if (exists $transcripts_exons_todo->{$tname} ){
 			$text .= "<br>";
 			$text .= " validated :".$transcripts_exons_validated->{$tname} if exists $transcripts_exons_validated->{$tname} ;
 			$text .= " todo:".$transcripts_exons_todo->{$tname} if (exists $transcripts_exons_todo->{$tname} ) ;
 
 	}
 	elsif (exists $transcripts_exons_validated->{$tname} ){
 			$text .= "<br> validated :".$transcripts_exons_validated->{$tname};
 			$text .= " todo:".$transcripts_exons_todo->{$tname} if (exists $transcripts_exons_todo->{$tname} ) ;
 
 	}
	push(@ths,$cgi->th({class=>"success"},$text));
	my $tr_id = $transcripts[$i]->id;
 	my $url2 = $url_image.$tr_id;
 	$url2 = $images->{$transcripts[$i]->id};
	 my $img = qq{<div style="padding:0px;min-width:36px;min-height:36px" ><img src="$url2"  align="top" style="box-shadow: 2px 2px 3px #aaaaaa;"  ></img></div>}; 
	 push(@tds,$cgi->td({onClick=>"zoomCoverage('$tr_id')",style=>"padding:0px;background-color:#DFF0D8"},$img));
	
	
 }
 $out.=print_lines(\@ths,\@tds);
$out.=$cgi->end_table();
$out.= html::end_cadre($cgi);
print $out;
#html::print_cgi($cgi,$out);
exit(0);

sub getLowTranscriptsOri {
	my $list_transcripts2 = $no->get($project->name,"minimum-".$padding."-".$utr."-".$intronic) ;#$no->get($project->name,"minimum-".$padding."-".$utr) ;
	my $list_transcripts2;
	unless ($list_transcripts2){
		
		$list_transcripts2 = preload_coverage::computeLowTranscripts($project,\@transcripts,$no,$intronic,$utr,$padding,1);
		$no->put($project->name,"minimum-".$padding."-".$utr."-".$intronic,$list_transcripts2) ;
		#$no->put($project->name,"minimum-".$padding."-".$utr,$list_transcripts2) ;
	}
	my $patients =  $project->get_list_patients($patient_names,",");
my $list_transcripts;
foreach my $t (@transcripts){
	my @z;
	foreach my $p (@$patients){
		warn $t->genomic_span->as_string;
		my $n = scalar($t->genomic_span->as_string);
		
		my $array =   $p->depthIntspan($t->getChromosome->name,$t->genomic_span);
		warn $t->id;
		push(@z,$list_transcripts2->{$t->id}->{$p->{id}})
	}
	my $min = min(@z);
	#warn $min;
	next if $min >= $limit;
	$list_transcripts->{$t->id} ++; 
}
return $list_transcripts;
}

sub getLowTranscripts {

my $patients =  $project->get_list_patients($patient_names,",");
my $list_transcripts;
foreach my $t (@transcripts){
		warn $t->name;
	my @z;
	my $find;
	foreach my $p (@$patients){
		my $n = scalar($t->genomic_span->as_array);
		
		my $array =   $p->depthIntspan($t->getChromosome->name,$t->getIntspan($padding,$utr,$intronic));
		my $min = min(@$array);
		if ($min<=$limit){
			$list_transcripts->{$t->id} ++; 
			last;
		}
		
	}
	
}
return $list_transcripts;
}
sub uri_image {
	my ($projects,$transcripts) = @_;
	
	my $fork = 5;
	my $nb = int(scalar(@$transcripts)/($fork*2))+1;
#	warn $nb;
	my $genes ;
	foreach my $t (@$transcripts){
		my $gene = $t->getGene;
		next if exists $genes->{$gene->id};
		$genes->{$gene->id} = $gene;
	}
	$nbgenes = scalar(keys %$genes);
	#$transcripts = [values %$genes];
	my $pm = new Parallel::ForkManager($fork);
	my $iter = natatime $nb, @$transcripts;
	my @t_final;
	print qq{<div style="visibility: hidden">};
	my $images;
	$pm->run_on_finish(
    sub { 
    	my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$h)=@_;
  
    	unless (defined($h) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
			$nbtp += $images->{nb};
			delete $images->{nb};
		foreach my $k (keys %{$h}){
			$images->{$k} = $h->{$k};
		}
    }
    );
  
	
  	$project->buffer->dbh_deconnect();
  	$|=1;
  	my $t =time;
  	
 	 while( my @tmp = $iter->() ){
 	 		my $pid = $pm->start and next;
 	 	
			$project->buffer->dbh_reconnect();
			my $himages ={};
			my $znb =0;
			my $dj;
			my $nbtp = 0;
 	 	foreach my $tr1  ( @tmp){ 
 	 			print "." if $znb %20 == 0;
 	 			$znb ++;
 	 			#my $gene = $tr1->getGene;
 	 			#next if exists $dj->{$gene->id};
 	 			 #$dj->{$gene->id} ++;
 	 			my $res;
 	 			my $no = $tr1->getChromosome->lmdb_image_transcripts_uri("r");
 	 			my $z  = $no->get($tr1->id);
 	 			$nbtp ++ if $z->{alert};
 	 			next unless  $z->{alert};
 	 			if ($z ){
 	 				$himages->{$tr1->id} = $z->{uri};
 	 				next;
 	 			}
 	 	
 	 		if ($project->isNoSqlDepth){	
				 $res  = image_coverage::image_depth_lmdb ($patients, $tr1,$intronic,$utr, $padding, $limit );
			
 	 		}
 	 		else {
 	 			 $res  = image_coverage::image ($patients, $tr1,$intronic,$utr, $padding, $limit );
 	 		}
 	 		next unless exists $res->{alert};
			my $uri = URI->new("data:");
			$uri->media_type("image/png");
			$uri->data($res->{image}->png);
			$himages->{$tr1->id} = $uri;
		
 	 	}
 	 		$himages->{nb} = $nbtp;
 	 	$pm->finish(0,$himages);
	}
	$pm->wait_all_children();
	$project->buffer->dbh_reconnect();
	print qq{</div>};
	
	return ($images);
}


sub print_lines{
	my ($ths,$tds) = @_;
	my $out = $cgi->start_Tr().join("\n",@$ths).$cgi->end_Tr();
	 $out .= $cgi->start_Tr().join("\n",@$tds).$cgi->end_Tr();
	$ths=[];
	$tds=[];
	return $out;
}