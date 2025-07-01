#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/";
use lib "$Bin/../../GenBo";
use lib "$Bin/../../GenBo/lib/obj-nodb";
use lib "$Bin/../../GenBo/lib/obj-nodb/packages/";
use lib "$Bin/../../GenBo/lib/obj-nodb/packages/cache/polydiag/";

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
use Storable qw/thaw/;
use coverage;
use validationQuery;
use lib "$Bin/../../../packages/validation_variation";
use draw_cnv;
use infos_coverage_exons;
use lib "$Bin/../../../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../../../GenBo/lib/obj-nodb/packages/cache/polydiag/";
use preload_coverage;
use update_variant_editor;
use table_dude;
use List::Util qw(
  uniq
);

use polyweb_dude;
$| = 1;
my $buffer = GBuffer->new();
my $CSS    = qq{
	<style type="text/css"> 

 
   .container_grid {
  		display: flex;
  		flex-wrap:  wrap ;
  		background-color:rgb(240, 234, 214)
	}
	.item_grid {
  align-self: flex-start;
  margin: auto;
  padding :1px;
}
.circle_tt {
  width:30px;
  height:30px;
  border-radius:30px;
  font-size:12px;
  color:#fff;
  line-height:30px;
  text-align:center;
  background:#ff;
  box-shadow: 1px 1px 2px 1px rgba(0, 0, 0, .4);
  position: relative;
  top: -8px;
  left: -5px;
}
</style>

};

my $array_quality = [ "high", "high", "medium", "low", "all" ];

my $cgi = new CGI();
my $out;
my $prg = $cgi->url( -relative => 1 );
my $url_image = url( -absolute => 1 );
$url_image =~ s/$prg/image_coverage.pl/;
if ( $cgi->param('test_cache') ) {
	html::print_cgi_header($cgi);
	exit(0);
}
my $order         = $cgi->param('order');
my $value_quality = 4;
$value_quality = $cgi->param('quality') if $cgi->param('quality');
$value_quality = 4 if $value_quality >=3;
my $patient_name = $cgi->param('patients');
my $quality      = $array_quality->[$value_quality];

#$quality ='medium';
my $project_name = $cgi->param('project');
my $only_low;
$only_low = 1 if $cgi->param('only_low') == 1;
my $project    = $buffer->newProjectCache( -name => $project_name );
my $panel_name = $cgi->param('panel_name');
$panel_name = "DI_DI44" unless $panel_name;

#$project->setPanel("$panel") if $panel;
my $utr      = $cgi->param('utr') + 0;
my $intronic = $cgi->param('intronic') + 0;

my $limit   = $cgi->param('limit');
my $padding = $cgi->param('span');

my @transcripts_cgi;
my $gene_id = $cgi->param('gene');
my $panel;
$panel_name = "ACMG-Defidiag";

my $force_visualisation = $cgi->param('force_visualisation');

my $transcripts_arg = $cgi->param('transcripts');
my $genes_arg       = $cgi->param('genes');
my $filter;
my @selected_patients;

if ($value_quality and not $genes_arg) {
	my $p   = $cgi->param('patients');
	my $pat = $project->getPatient($p);
	
	  
	my @lTypes;
	if ($value_quality == 1) { @lTypes = ('high'); }
	if ($value_quality == 2) { @lTypes = ('high','medium'); }
	if ($value_quality == 3) { @lTypes = ('high','medium', 'low'); }
	if ($value_quality == 4) { @lTypes = ('high','medium', 'low'); }
	my $this_types = join(', ', @lTypes);
	if ($project->isDiagnostic()) { @lTypes = ('high','medium', 'low'); }
	my $hRes = update_variant_editor::get_hash_genes_dude($pat, 'ids', \@lTypes,undef,1);
#	warn Dumper $hRes;
	if ($project->isDiagnostic()) { 
		foreach my $g_id (keys %{$hRes}) {
			if ($value_quality == 1 and not exists $hRes->{$g_id}->{'high'}) {
				delete $hRes->{$g_id};
			}
			if ($value_quality == 2) {
				next if (exists $hRes->{$g_id}->{'high'});
				next if (exists $hRes->{$g_id}->{'medium'});
				delete $hRes->{$g_id};
			}
		}
	}
	$genes_arg = join(';', keys %{$hRes});
	if (not $genes_arg or $genes_arg eq 'no_result') {
		html::print_cgi_header($cgi);
		print $CSS;
		$out .= $cgi->start_table({class=>"table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 10px;font-family:  Verdana;", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
		$out .=  $cgi->td({style=>"text-align: center;vertical-align:middle"},'');
		$out .=  $cgi->td({style=>"text-align: center;vertical-align:middle"},"<span><b>No result for $this_types... Sorry..</b></span>");
		$out .=  $cgi->end_table();
		$out .=  "</div></div>";
		print $out;
		exit(0);
	}
}

if ($genes_arg) {
	if ( $cgi->param('patients') ) {
		my $p   = $cgi->param('patients');
		my $pat = $project->getPatient($p);
		foreach my $m ( @{ $pat->getFamily->getMembers } ) {
			push( @selected_patients, $m );
		}

	}
	my $ps = $project->getPatients();
	my @lgenesIds = split( ";", $genes_arg );
	@lgenesIds = split( ",", $genes_arg ) if ($genes_arg =~ /,/);
	my $limit = 3;
	$limit = 2 if scalar(@lgenesIds) > 200;
	$limit = 1 if scalar(@lgenesIds) > 500;
	foreach my $ga ( @lgenesIds ) {
		my $gene = $project->newGene($ga);
		next unless $gene;
	
		
		my $tt =0;
		my $no = $ps->[0]->getTranscriptsDude("r");
		foreach my $t ( sort { $b->length <=> $a->length }
			@{ $gene->getTranscripts } )
		{

			my $matrix = $no->get( $t->id );
			next unless $matrix;
			#if ($project->isDiagnostic)
			#next unless ($matrix && ($project->isDiagnostic));
			#	next unless $t->ccds_name;
			$tt++;
			push( @transcripts_cgi, $t->id );
			#last if $tt > $limit ;
			
		}
		
	}
}

#my
else {
	if ( $quality eq "all" ) {
		@transcripts_cgi =
		  @{ $project->getListTranscripts( { transcripts => $transcripts_arg } )
		  };
	}
	else {

		$filter = 1;
		foreach my $patient ( @{ $project->getPatients } ) {
			#next if $patient->name() ne "F504_1169";
			if ( $patient_name ne 'all' ) {
				next if $patient_name && $patient->name ne $patient_name;
			}
			next if $patient->status ne "2";
			push( @selected_patients, $patient );
			my $no = $patient->getTranscriptsDude("r");
			my $t  = $no->get($quality);
			next unless $t;
			push( @transcripts_cgi, @{$t} );

		}
		@transcripts_cgi = uniq(@transcripts_cgi);
	}
}

my $transcripts_exons_todo      = {};
my $transcripts_exons_validated = {};
my $patient_names = $cgi->param('patients');
my @transcripts;
if ( scalar(@transcripts_cgi) == 1 and $transcripts_cgi[0] eq 'all' ) {
	@transcripts = @{ $project->getTranscripts() };
}
else {
	@transcripts = @{ $project->newTranscripts( \@transcripts_cgi ) };
}

@transcripts = sort {
		 $a->getChromosome->id <=> $b->getChromosome->id
	  || $a->start <=> $b->start
} @transcripts;

my $nbgenes = 0;
my $nbt     = scalar(@transcripts);

my $nbtp = 0;    #scalar(@transcripts);

#@transcripts = splice (@transcripts ,$cgi->param('start'),$cgi->param('step')) if $cgi->param('step') && !$only_low ;

my $col = 12;

$col = 10 if scalar( @{ $project->getPatients() } > 15 );
$col = 6  if scalar( @{ $project->getPatients() } > 30 );
$col = 5  if scalar( @{ $project->getPatients() } > 40 );
my $size   = "";
html::print_cgi_header($cgi);
print $CSS;
my $images = uri_image( $project, \@transcripts, $force_visualisation);
#warn scalar(@transcripts);
my $nb_images;
foreach my $t (@transcripts) {
	if (   exists $images->{ $t->id }->{type}->{red}
		or exists $images->{ $t->id }->{type}->{red} )
	{
		$nb_images++;
	}
}
if ( $nb_images < $col ) {
	$size =
		"max-width:"
	  . ( scalar( @{ $project->getPatients } ) * 5 + (130) ) * $nb_images
	  . "px";
}

$out .= $cgi->start_table(
	{
		class =>
"table table-striped table-condensed table-bordered  table-mybordered",
		style => "font-size: 9px;font-family:Verdana;$size"
	}
);
my $nb = 0;

$url_image .=
	"?project="
  . $project_name
  . "&limit=$limit&span=$padding&utr=$utr&intronic=$intronic&patients=$patient_names&transcript=";

my @ths;
my @tds;

my $no = $project->noSqlCoverage();

my $patients = $project->getPatients();
my $nb       = 0;
my @problem_transcripts;

@transcripts = grep { exists $images->{ $_->id } } @transcripts;
my $nbc = 0;
for ( my $i = 0 ; $i < @transcripts ; $i++ ) {
	if ( $nbc % $col == 0 && $nbc > 0 ) {
		$out .= print_lines( \@ths, \@tds );
		@ths = ();
		@tds = ();

	}

	my $tname = $transcripts[$i]->name;
	# 	my $class = "nobg";
	# $class ="th1" if $z%2 ==0;
	my $chr_name = $transcripts[$i]->getChromosome()->name();
	my $text =
	  qq{<div class="circle_tt">$chr_name</div>}
	  . $cgi->start_table(
		{ style => "position:relative;top:-15px;left:2px;" } );

#$text .=  $cgi->start_Tr().$cgi->td(qq{<div class="circle">21</div>}).$cgi->end_Tr();
	$text .= $cgi->start_Tr()
	  . $cgi->td(
		{
			class => "btn btn-xs btn-primary",
			style => "background-color: #EEEFEE;font-size: 9px;font-family:Verdana;color:black"
		},
		$transcripts[$i]->getGene()->external_name()
	  ) . $cgi->end_Tr();
	  
	$text .= $cgi->start_Tr() . $cgi->td($tname) . $cgi->end_Tr();
	$text .= $cgi->start_Tr() . $cgi->td($transcripts[$i]->getChromosome()->name().":".$transcripts[$i]->start."-".$transcripts[$i]->end) . $cgi->end_Tr();
	$text .= $cgi->start_Tr() . $cgi->td($transcripts[$i]->getGene()->description) . $cgi->end_Tr();

#  $text .=  $cgi->start_Tr().$cgi->td($transcripts[$i]->getChromosome()->name().":".$transcripts[$i]->start."-".$transcripts[$i]->end).$cgi->end_Tr();
	$text .= $cgi->end_table();

# <button type="button" class="btn btn-xs btn-primary " style="background-color: #FF8800;font-size: 7px;font-family:  Verdana;;color:white" onclick="window.open(&quot;https://gnomad.broadinstitute.org/gene/ENSG00000172264&quot;)">MACROD2 <b><u></u></b></button>
# my $text = qq{<div class="circle">21</div><div>}.$transcripts[$i]->getGene()->external_name()."</div><br>".$tname."<br>". $transcripts[$i]->getChromosome()->name().":".$transcripts[$i]->start."-".$transcripts[$i]->end."<br><small>".$transcripts[$i]->getGene()->description."</small>";
#	$text.= qq{<div class="circle">21</div>};
	my $box = "box-shadow: 1px 1px 2px 1px rgba(100, 100, 100, .4);";

	#my $color = "background-color:rgb(63, 105, 170)";
	my $color = "background-color:#34495D";
	$color = "background-color:#DEDEDE;color:black";
	
	if (   exists $images->{ $transcripts[$i]->id }->{type}->{red}
		&& exists $images->{ $transcripts[$i]->id }->{type}->{blue} )
	{
		$box   = "box-shadow: 1px 1px 2px 1px rgba(157, 94, 183, .7);";
		$color = "background-color:#9D5EB7";
	}
	elsif ( exists $images->{ $transcripts[$i]->id }->{type}->{red} ) {
		$box   = "box-shadow: 1px 1px 2px 1px rgba(255, 0, 0, .4);";
		$color = "background-color:#FF6F61";
	}
	elsif ( exists $images->{ $transcripts[$i]->id }->{type}->{blue} ) {
		$box   = "box-shadow: 1px 1px 2px 1px rgba(0, 0, 255, .4);";
		$color = "background-color: #3398DB";
	}
	else {
		$color = "background-color: grey;";
		#next;
	}
	$nbc++;

	push( @ths, $cgi->th( { class => "success", style => "$color" }, $text ) );
	my $tr_id = $transcripts[$i]->id;
	my $url2  = $url_image . $tr_id;
	$url2 = $images->{ $transcripts[$i]->id }->{uri};

	my $img =
qq{<div style="padding:0px;min-width:36px;min-height:36px" ><img src="$url2"  align="top" style="$box "  ></img></div>};
	push(
		@tds,
		$cgi->td(
			{
				onClick => "zoomCnv('$tr_id','50')",
				style   => "padding:0px;background-color:#ECF0F1;"
			},
			$img
		)
	);

}
#warn Dumper @tds;
$out .= print_lines( \@ths, \@tds );
$out .= $cgi->end_table();
$out .= html::end_cadre($cgi);
print $out;

#html::print_cgi($cgi,$out);
exit(0);
sub uri_image {
	my ( $projects, $transcripts, $force_visualisation) = @_;
	my $patients = $project->getPatients();
	my $fork     = 16;
	my $nb       = int( scalar(@$transcripts) / ( $fork * 2 ) ) + 1;
	print qq{<div style="visibility: hidden">};
	#	warn $nb;
	my $genes;
	my $zz =0;
	foreach my $t (@$transcripts) {
		$zz++;
		print "*" if $zz % 100 ==0;
		my $gene = $t->getGene;
		next if exists $genes->{ $gene->id };
		$genes->{ $gene->id } = $gene;
	}
	$nbgenes = scalar( keys %$genes );

	#$transcripts = [values %$genes];
	my $pm   = new Parallel::ForkManager($fork);
	my $iter = natatime $nb, @$transcripts;
	my @t_final;

	my $images;
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $h ) = @_;

			unless ( defined($h) or $exit_code > 0 ) {
				print
				  qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
			$nbtp += $images->{nb};
			delete $images->{nb};

			foreach my $k ( keys %{$h} ) {
				if ( $k eq "nb" ) {
					$images->{$k} = $h->{$k};
					next;
				}
				$images->{$k}->{type} = $h->{$k}->{type};
				$images->{$k}->{uri}  = $h->{$k}->{uri};
			}
		}
	);
	$project->disconnect();
	$| = 1;
	my $t = time;
	my $pp;
	if (not $patient_name or $patient_name eq 'all') {
		@selected_patients = @{$project->getPatients()};
		$pp = $selected_patients[0];
		$force_visualisation = 1;
	}
	else { $pp = $project->getPatient($patient_name); }
	while ( my @tmp = $iter->() ) {
		my $pid = $pm->start and next;
		my $himages = {};
		my $znb     = 0;
		my $dj;
		my $nbtp = 0;

		foreach my $tr1 (@tmp) {
#			warn $tr1;
			print "." if $znb % 20 == 0;
			
			$znb++;
			my $coverage = polyweb_dude->new(
				patients          => [$pp],
				transcript        => $tr1,
				limit             => $limit,
				selected_patients => \@selected_patients
			);
			$coverage->init_matrices;

			unless ($force_visualisation) {
				if ($patient_name && $patient_name ne "all"){
					next  unless (exists $coverage->quality()->{$patient_name}->{dup} or  exists $coverage->quality()->{$patient_name}->{del});
				}
				if ($patient_name ne "all"){
					next unless (exists $coverage->quality->{del} or exists $coverage->quality()->{dup});
				}
			}
			
			  if ($filter) {
				next
				  if ( $coverage->quality->{dup} > 5
					&& $coverage->quality->{del} > 5 );
				next if ( $coverage->quality->{grey} > 10 );    # &&  $coverage->quality->{del} >2);
			}

			$coverage = polyweb_dude->new(
				patients          => $patients,
				transcript        => $tr1,
				limit             => $limit,
				selected_patients => \@selected_patients
			);
			$coverage->init_matrices;


			#	next if $coverage->error == 0 && $only_low ==1 ;
			$nbtp++;
			my ( $image, $type ) = $coverage->image
			  ; #image_coverage::image_depth_lmdb_from_matrix($tr1,$patients,,$matrix,30);

			my $uri = URI->new("data:");

			$uri->media_type("image/png");
			$uri->data( $image->png );

			$himages->{ $tr1->id }->{uri}  = $uri;
			$himages->{ $tr1->id }->{type} = $type;
			next;

		}
		$himages->{nb} = $nbtp;
		$pm->finish( 0, $himages );
	}
	$pm->wait_all_children();
	$project->buffer->dbh_reconnect();
	print qq{</div>};

	return ($images);
}

sub debug_matrix {
	my ( $patient, $transcript ) = @_;
	my $exons = $transcript->getExons();
	foreach my $exon ( sort { $a->end * $a->strand <=> $b->end * $b->strand }
		@$exons )
	{
		my $pos  = $exon->return_start_end_no_utr( padding => 50 );
		my $xa   = [];
		my $min  = 0;
		my $mean = 0;

		if ($pos) {
			$xa = $patient->depth( $exon->getChromosome->name,
				$pos->{start}, $pos->{end} )
			  if $pos;
			$min  = min(@$xa);
			$mean = int( sum(@$xa) / scalar(@$xa) );
		}
		warn $min . " " . $mean;
	}
}

sub print_lines {
	my ( $ths, $tds ) = @_;
	my $out = $cgi->start_Tr() . join( "\n", @$ths ) . $cgi->end_Tr();
	$out .= $cgi->start_Tr() . join( "\n", @$tds ) . $cgi->end_Tr();
	$ths = [];
	$tds = [];
	return $out;
}
