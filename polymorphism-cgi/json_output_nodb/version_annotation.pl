#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/validation_variation"; 
use lib "$Bin/../cache_nodb/scripts/";

require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
require "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update.pm";
use html; 
use connect;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use export_data;
use Carp;
use JSON;
use Time::Piece;
#use Text::CSV::Slurp; 

#ce script est utilisé pour renvoyer les noms des projets d'un utilisateur après sa connection à l'interface ainsi que pour telecharger des fichiers via l'interface ou encore exporter les tableaux de données au format xls

my $cgi    = new CGI;
#chargement du buffer 
my $buffer = GBuffer->new;

my $last_annot = $buffer->getQuery()->getMaxPublicDatabaseVersion();

my $version = $cgi->param('version');
my $genome = $cgi->param('genome');
my $color = $cgi->param('color');
my ($gencode,$annot,$z) = split(/\-/,$version);
if ($z){
	$gencode = $gencode.".".$annot;
	$annot = $z;
}
my $h_infos_last_annot = get_infos_from_annotion($last_annot);
my $h_infos_annot = get_infos_from_annotion($annot);


#write HTML CODE 
html::print_header_table_version($cgi);

my $color_title = qq{style="color:black"'};
my $color_card_project;
my $color_card_last = qq{style="background-color:rgb(39, 174, 96);color:white"'};
if ($color) {
	$color_card_project = qq{style="background-color:$color;color:white"'};
}

if ($last_annot == $annot) {
	print qq{
		<section class="pricing py-5" style="padding:30px;">
		  <div class="container">
		    <div class="row">
		      <!-- Free Tier -->
		      <div class="col-lg-4"></div>
		      <div class="col-lg-4">
		        <div class="card mb-5 mb-lg-0" $color_card_project>
		          <div class="card-body" style="padding-top:5px;">
		            <h5 class="card-title text-muted text-uppercase text-center" $color_title>Project Version</h5>
		            <h6 class="card-price text-center">$version</h6>
		            <hr>
		            <ul class="fa-ul">
		            
						<table class="card-table table" style="padding-right:15px">
							<thead>
								<tr>
									<th scope="col"></th>
									<th scope="col"><b>Database</b></th>
									<th scope="col"><b>Version</b></th>
									<th scope="col"><b>Date</b></th>
								</tr>
							</thead>
							<tbody>
	};
}
else {
	
	
	print qq{
		<section class="pricing py-5" style="padding:30px;">
		  <div class="container">
		    <div class="row">
		      <!-- Free Tier -->
		      <div class="col-lg-4">
		        <div class="card mb-5 mb-lg-0" $color_card_project>
		          <div class="card-body" style="padding-top:5px;">
		            <h5 class="card-title text-muted text-uppercase text-center" $color_title>Project Version</h5>
		            <h6 class="card-price text-center">$version</h6>
		            <hr>
		            <ul class="fa-ul">
		            
						<table class="card-table table" style="padding-right:15px">
							<thead>
								<tr>
									<th scope="col"></th>
									<th scope="col"><b>Database</b></th>
									<th scope="col"><b>Version</b></th>
									<th scope="col"><b>Date</b></th>
								</tr>
							</thead>
							<tbody>
	};
}

foreach my $cat_name (sort keys %$h_infos_last_annot){
	my ($cat_version, $cat_date, $glyph_ok_not);
	my ($li_start, $li_end);
	if (exists $h_infos_annot->{$cat_name}) {
		$cat_version = $h_infos_annot->{$cat_name}->{version};
		$cat_date = $h_infos_annot->{$cat_name}->{date};
		$glyph_ok_not = "glyphicon-ok";
		if ($h_infos_annot->{$cat_name}->{version} ne $h_infos_last_annot->{$cat_name}->{version}) {
			$glyph_ok_not = "glyphicon-remove";
			$li_start = "<span class='text-muted'>";
			$li_end = "</span>";
		}
	}
	else {
		$cat_version = '<span class="glyphicon glyphicon-ban-circle"></span>';
		$cat_date = '<span class="glyphicon glyphicon-ban-circle"></span>';
		$glyph_ok_not = "glyphicon-remove";
		$li_start = "<span class='text-muted'>";
		$li_end = "</span>";
	}
	print qq{
		<tr>
			<td>$li_start<span class="glyphicon $glyph_ok_not" float="left"></span></td>
			<td>$li_start $cat_name $li_end</td>
			<td>$li_start $cat_version $li_end</td>
			<td>$li_start $cat_date $li_end</td>
		</tr>
	};
	
	#<li class="text-muted"><span class="fa-li"><i class="fas fa-times"></i></span>Monthly Status Reports</li>
}

print qq{
    </tbody>
  </table>
};

if ($last_annot == $annot) {
	print qq{
		            </ul>
		          </div>
		        </div>
		      </div>
		    </div>
		  </div>
		</section>
	};
}
else {
	print qq{
		            </ul>
		          </div>
		        </div>
		      </div>
			  <div class="col-lg-4">
			  	<center>
			  		<span class="glyphicon glyphicon-arrow-right" style="padding-top:50%;font-size:100px;"></span>
		  		</center>
			  </div>
		      <div class="col-lg-4">
		        <div class="card mb-5 mb-lg-0" $color_card_last>
		          <div class="card-body" style="padding-top:5px;">
		            <h5 class="card-title text-muted text-uppercase text-center" $color_title>Last Version Availabled</h5>
		            <h6 class="card-price text-center">$gencode-$last_annot</h6>
		            <hr>
		            <ul class="fa-ul">
		            
						<table class="card-table table" style="padding-right:15px">
							<thead>
								<tr>
									<th scope="col"></th>
									<th scope="col"><b>Database</b></th>
									<th scope="col"><b>Version</b></th>
									<th scope="col"><b>Date</b></th>
								</tr>
							</thead>
							<tbody>
	};
	
	
	foreach my $cat_name (sort keys %$h_infos_last_annot){
		my $cat_version = $h_infos_last_annot->{$cat_name}->{version};
		my $cat_date = $h_infos_last_annot->{$cat_name}->{date};
		my $glyph_ok_not = "glyphicon-ok";
		print qq{
			<tr>
				<td><span class="glyphicon $glyph_ok_not" float="left"></span></td>
				<td>$cat_name</td>
				<td>$cat_version</td>
				<td>$cat_date</td>
			</tr>
		};
	}
	
	print qq{
		    </div>
		  </div>
		</section>
	};
}

exit(0);




sub get_infos_from_annotion {
	my ($annot) = @_;
	my $lines;
	foreach my $database (keys %{$buffer->public_data->{$annot}}){
		my $file = $buffer->public_data_root."/HG19/"."/".$buffer->public_data->{$annot}->{$database}->{config}->{directory}."/description.json";
		$file = $buffer->public_data_root."/HG19/"."/".$buffer->public_data->{$annot}->{$database}->{config}->{directory}."/version.json" unless -e $file;
		next unless -e $file;
		my $h;
		$h->{name} = $database;
		$h->{version} = $buffer->public_data->{$annot}->{$database}->{version};
		open(JSON, $file);
	 	my $desc = decode_json <JSON>;
		$h->{date} = $desc->{date};
		$lines->{$database} =$h 
	}
	my $h_infos;
	foreach my $l (sort{$a->{name} cmp $b->{name}} values %$lines){
		my $cat_name = $l->{name};
		$h_infos->{$cat_name}->{version} = $l->{version};
		$h_infos->{$cat_name}->{date} = $l->{date};
	}
	return $h_infos;
}
