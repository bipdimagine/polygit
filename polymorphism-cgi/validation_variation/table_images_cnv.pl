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
use Parallel::ForkManager;
use html;
use List::MoreUtils qw{ natatime };

#use Set::;
use Storable qw/store thaw retrieve freeze/;
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
use image_coverage;
use preload_coverage;
use Digest::MD5::File qw( md5_hex file_md5_hex file_md5);

$| = 1;
my $cgi = new CGI();
#
html::print_cgi_header($cgi);
my $cc = $cgi->param('no_cache');
exit(0) if $cc == 1;
my %images;
my $prg = $cgi->url( -relative => 1 );
my $url_image = url( -absolute => 1 );
$url_image =~ s/$prg/image_cnv.pl/;

my $project_name  = $cgi->param('project');
my $order         = $cgi->param('order');
my $patient_names = $cgi->param('patients');

#ProjectCache or not that's the questions , do you have run dude or not ?
my $buffer  = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );

##
my $no_cache = $project->get_lmdb_cache_cnv("w");
my $version = 1;
my $key = return_uniq_keys($cgi);	
my $cache_id = md5_hex("polydiag_cnv_".join(";",@$key).".$version");
my $text = $no_cache->get_cache($cache_id);
if ($text){
	$| =1;
	print $text;
	exit(0);
}
##

if ( $project->isDude ) {
	my $valid_cache = $cgi->param('test_cache');
	exit(0) if $valid_cache;
	$buffer  = GBuffer->new();
	$project = $buffer->newProjectCache( -name => $project_name );
}


#CAS MAJ annotations dans le cache. Evite de relancer le CNV_LITE qui n est pas lie a l annotation 
delete $project->{noSqlCnvsDir};
$project->noSqlCnvsDir();

my $panel = $cgi->param('panel');
$project->setPanel("$panel") if $panel;
my $patients = $project->get_only_list_patients( $patient_names, "," );
my $out;
my $utr      = $cgi->param('utr') + 0;
my $intronic = $cgi->param('intronic') + 0;

my $limit   = $cgi->param('limit');
my $padding = $cgi->param('span');

my @transcripts_cgi;

my $cgi_transcript = $cgi->param('transcripts');
if ( $cgi_transcript eq "all" ) {
	@transcripts_cgi = @{ $project->bundle_transcripts() };
}
else {
	@transcripts_cgi = split( ",", $cgi_transcript );
}

my $capture = $project->getCaptures()->[0];
my $vquery;

my $runs = $project->getRuns();
my $res;
foreach my $run (@$runs) {
	my $infos = $run->getAllPatientsInfos();
	my @p =
	  map { $_->{patient} } grep { $_->{project} eq $project_name } @$infos;
	$res->{ $run->id } = \@p;
}

if ( $project->isDude ) {
	@transcripts_cgi = select_transcripts( $project, \@transcripts_cgi );
}
else {
	$out .= html::print_cadre( $cgi, "CNV" );
	$out .= $cgi->start_table(
		{
			class => "table table-condensed table-bordered table-mybordered",
			style => "font-size: 08px;font-family:  Verdana;"
		}
	);
	$out .= $cgi->start_div( { class => "container" } );
	$out .=
	  '<br>Please contact BioInformatic Platform to add CNV Caller.<br><br>';
	$out .= $cgi->end_div();
	$out .= $cgi->end_table();
	$out .= html::end_cadre( $cgi, "CNV" );
	print $out;
	exit(0);
}
my @transcripts = sort{$a->getGene->external_name cmp $b->getGene->external_name} @{ $project->newTranscripts( \@transcripts_cgi ) };



$out .= html::print_cadre( $cgi, "CNV" );
$out .= $cgi->start_table(
	{
		class => "table table-condensed table-bordered table-mybordered",
		style => "font-size: 08px;font-family:  Verdana;"
	}
);

#$out.= $cgi->start_div({class=>"container"});
#$out.= $cgi->start_div({class=>"row"});

my $col = 12;
$col = 10 if scalar( @{ $project->getPatients() } > 20 );
$col = 7  if scalar( @{ $project->getPatients() } > 30 );
$col = 5  if scalar( @{ $project->getPatients() } > 40 );
my $nb = 0;
my @ths;
my @tds;

#for (my $i=0;$i<@transcripts_by_ordered_chrom;$i++){
for ( my $i = 0 ; $i < @transcripts ; $i++ ) {
	if ( $i % $col == 0 && $i > 0 ) {
		$out .= print_lines( \@ths, \@tds );
		@ths = ();
		@tds = ();

	}
	my $tname     = $transcripts[$i]->name;
	my $gene_name = $transcripts[$i]->getGene->external_name;
	my $chr_name  = $transcripts[$i]->getChromosome->name;
	my $tr_id     = $transcripts[$i]->id;

	#last unless $transcripts[$z];
	my $td = $cgi->start_table(
		{
			class => "table table-condensed table-bordered table-mybordered",
			style => "font-size: 08px;font-family:  Verdana;"
		}
	);

	foreach my $r ( keys %$res ) {

		#my $patient_names = join(",",@{$res->{$r}});
		my $url_image_tmp =
			$url_image
		  . "?project="
		  . $project_name
		  . "&limit=$limit&span=$padding&utr=$utr&intronic=$intronic&patients=$patient_names&run_id=$r&transcript=";
		my $url2 = $url_image_tmp . $transcripts[$i]->id();
		my $toto = $images{ $transcripts[$i]->id() };

		#$ret->{image}->png;
		$td .= $cgi->td(
			{
				style =>
'background-color:#FCF8E3;background-image:../../images/polyicons/layer_his_add.png',
				onClick => "zoomCnv('$tr_id','$r','$panel')"
			},
qq{<img class="load" src="$toto"  align="top" style="box-shadow: 2px 2px 3px #aaaaaa;" lowsrc="../../images/polyicons/layer_his_add.png"></img>}
			  . ""
		);
		########################################
		last; #attention ici c'et spou reviter de rpeter la meme image 2 fois je c alcul l'image pour tous les patients du projet et non pas uniquemùent du run 
		########################################
		
	}
	$td .= $cgi->end_table();
	push( @tds, $cgi->td($td) );

	push(
		@ths,
		$cgi->th(
			{ class => 'warning' },
			qq{ $gene_name &nbsp; chr:$chr_name<br>$tname}
		)
	);

}
$out .= print_lines( \@ths, \@tds );
$out .= $cgi->end_table();
$out .= html::end_cadre( $cgi, "CNV" );
##
$no_cache->put_cache_text($cache_id,$out,2400);
$no_cache->close();
exit(0) if $cgi->param('pipeline');

print $out;

#html::print_cgi($cgi,$out);
exit(0);

sub print_lines {
	my ( $ths, $tds ) = @_;
	my $out = $cgi->start_Tr() . join( "\n", @$ths ) . $cgi->end_Tr();
	$out .= $cgi->start_Tr() . join( "\n", @$tds ) . $cgi->end_Tr();
	$ths = [];
	$tds = [];
	return $out;
}


sub select_transcripts {
	my ( $project, $list_transcripts ) = @_;

	my $fork = 3;
	my $nb   = int( scalar(@$list_transcripts) / ( $fork + 1 ) );
	$nb =1 if $nb ==0;
	my $pm   = new Parallel::ForkManager($fork);
	my $iter = natatime $nb, @$list_transcripts;
	my @t_final;
	print qq{<div style="visibility: hidden">};
	my $proc;
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $h ) = @_;

			unless ( defined($h) or $exit_code > 0 ) {
				print
				  qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
			my $id = delete $h->{proc};
			warn "del->".$id;
			delete $proc->{$id};
			foreach my $k ( keys %{$h} ) {
				$images{$k} = $h->{$k};

				#push(@t_final,$k);
			}

			#push(@total,@{$h->{array}}) if $h->{array};
			#push(@total_vars,@{$h->{vars}}) if $h->{vars};
		}
	);

	$project->buffer->dbh_deconnect();
	$| = 1;
	my $t = time;
	my $id = time;
	
	while ( my @tmp = $iter->() ) {
		$id++;
		$proc->{$id} ++;
		my $pid         = $pm->start and next;
		my $transcripts = $project->newTranscripts( \@tmp );
		my $himages;
		$project->buffer->dbh_reconnect();
		my $final;

		foreach my $tr (@$transcripts) {
			my $tr_objs;
			print ".";
			my $debug;
		#	$debug = 1 ;#if ( $tr->name eq "ENST00000254457" );
			my $primers = $project->getPrimersByObjects($tr);
			my $find;
			my $max_length = 0;
			my @tab_string;
			warn "---> " . scalar(@$primers) if $debug;

			#	die() if $debug;
			foreach my $primer (@$primers) {
				my @tab;
				if ( $project->isDude ) {
					@tab = values( %{ $primer->{level} } );
				}
				else {
					foreach my $patient ( @{ $project->getPatients } ) {
						push( @tab, $primer->level($patient) );

					}
				}
				push( @$tr_objs, \@tab );

				my $string = join( "", @tab );

				#my $string = "ThisXlineXhasXsomeXx'sXinXit";
				my $count  = ( $string =~ tr/1// );
				my $count2 = ( $string =~ tr/2// );
				my $m      = max(@tab);

				$find= 1;
				last;
				#warn $m;
				if ( $m > 0 ) {
					$find = 1;
					last;
				}
			}

			#warn scalar(@{$tr_objs}) if $debug;
			#die() if $debug;
			if ( $find == 1 ) {

				#transpose

				#warn Dumper $tr_objs->[0] if $debug;
  
				#image_cnv_no_cache($patients,$transcript,$cgi);
				my $res;
				$res =
				  image_coverage::image_cache_cnv( $patients, $tr, $primers,
					$max_length );
				my $uri = URI->new("data:");
				$uri->media_type("image/png");
				$uri->data( $res->{image}->png );
				$himages->{ $tr->id } = $uri;
			}
		}
		$himages->{proc} = $id;
		$pm->finish( 0, $himages );
	}
	$pm->wait_all_children();
	$project->buffer->dbh_reconnect();

	#warn abs(time -$t);

	print "</div>\n";
	 if (keys %$proc){
	 	print "<h1>ERROR !!!!!!! </h1>";
	 	warn Dumper %$proc;
	 	exit(0);
	 }
	return keys %images;
}

###
#
####
sub return_uniq_keys {
my ($cgi) = @_;
	
my %hkeys = $cgi->Vars;
my @keys;
my $string;
foreach my $k  (sort {$a cmp $b} keys %hkeys){
	next if $k =~ /force/;
	next if $k =~ /user/;
	next if $k =~ /pipeline/;
	push(@keys,"$k");
	my $c = $hkeys{$k};
	$c =~ s/\"//g;
	$c =~ s/\+/ /g;
	push(@keys,$c);
}
my $dir_out   =$project->noSqlCnvsDir;
my $f2 = "$dir_out/primers.lite";
my @st = (stat($f2));
push(@keys, ($st[9].$st[11].$st[12]));
$f2 = "$dir_out/1";
my @st2 = (stat($f2));
push(@keys, ($st2[9].$st2[11].$st2[12]));
return \@keys;
}
#

