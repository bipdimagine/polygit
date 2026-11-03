#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/validation_variation"; 
use lib "$Bin/../cache_nodb/scripts/";


use connect;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use export_data;
use Carp;
use JSON;
use Time::Piece;
use Compress::Zlib qw(compress);

#ce script est utilisé pour renvoyer les noms des projets d'un utilisateur après sa connection à l'interface ainsi que pour telecharger des fichiers via l'interface ou encore exporter les tableaux de données au format xls

my $cgi    = new CGI;
my $site = {
	"psl" => "pitie",
	"nck" => "necker",
	"lyon" => "lyon",
	"bourgogne" => "dijon",
	"rouen" => "rouen",
	"strasbourg" => "strasbourg",
	"igbmc" => "strasbourg",
	"lyon" => "lyon",
	"bordeaux" => "bordeaux",
	"ap-hm" => "marseille"
	
};
#chargement du buffer 
my $buffer = GBuffer->new;

my $login = $cgi->param('login');
my $pwd   = $cgi->param('pwd');
my $keep_rna = $cgi->param('keep_rna');
my $only_rna = $cgi->param('only_rna');
my $quick_request = $cgi->param('quick_request');
getProjectLists($buffer,$login,$pwd);


 
sub getProjectLists {
	my ( $buffer, $login, $pwd ) = @_;
	my $isUserMagic = $buffer->getQuery()->isUserMagic($login);
	my $res = $buffer->getQuery()->getProjectListForUser($login, $pwd );
	
	my $h_capture_id_found;
	
	my @list_datas;
	foreach my $pro (sort @$res) {	
		my $debug;

		my $find;
		my (@c, $h_type);
		$pro->{'type_analyse'} = join(', ', sort keys %$h_type);
		
		next if $pro->{'type_analyse'} eq 'rnaseq' and not $keep_rna;
		next if $pro->{'type_analyse'} eq 'singlecell';
		$pro->{'type_analyse'} = 'target' if $pro->{'type_analyse'} ne "exome" && $pro->{'type_analyse'} ne "genome";
		next if not $pro->{name} =~ /NGS20[0-9][0-9]_[0-9]+/;
		
		my @samples_names;
		my %captures_id;
		my $is_rna;
		if (exists $pro->{patient_details}) {
			foreach my $patient_info (split(';', $pro->{patient_details})) {
				my ($pat_id, $pat_name, $pat_family, $pat_status, $pat_type, $pat_capture_id) = split(',', $patient_info);
				$captures_id{$pat_capture_id} ++;
				if ($quick_request) {
					push(@samples_names, $pat_family.';'.$pat_name.';'.$pat_status);
				}
				else {
					push(@samples_names, $pat_name);
				}
				$is_rna = 1 if $pat_type eq 'rna';
			}
		}
		else {
			my @samples =  @{$buffer->getQuery->getPatients($pro->{id})};
			foreach my $p (@samples){
				$captures_id{$p->{capture_id}} ++;
				my $pat_name = $p->{name};
				if ($quick_request) {
					push(@samples_names, $p->{family}.';'.$pat_name.';'.$p->{status});
				}
				else {
					push(@samples_names, $pat_name);
				}
				$is_rna = 1 if $p->{type} eq 'rna';
			}
			
		}
		next if $is_rna and not $keep_rna;
		
		if (not $quick_request) {
			foreach my $cid (sort keys %captures_id){
				if (not exists $h_capture_id_found->{$cid}) {
					my $capt =  $buffer->getQuery()->getCaptureInfos($cid);
					$h_capture_id_found->{$cid}->{name} = $capt->{name};
					$h_capture_id_found->{$cid}->{analyse} = $capt->{analyse};
				}
				push(@c, $h_capture_id_found->{$cid}->{name});
				$h_type->{$h_capture_id_found->{$cid}->{analyse}} = undef;
			}
		}
		$pro->{capture_name} = join(", ",@c);
		if (not exists $pro->{genome}) {
			$pro->{genome} = $buffer->getQuery->getReleaseGenome($pro->{id});
		}
		
		my @list_project;
		push(@list_project, $pro->{name});
		push(@list_project, $pro->{description});
		my $nb_samples = scalar(@samples_names);
		
		my $project_name = $pro->{name};
		my $patients_text = join(",", @samples_names);
		#my $cmd_resume = qq{window.open('resume_project?project=$project_name', '_blank')};
		my $cmd_resume = qq{go_to_resume_project_page('$project_name', '$patients_text')};
		my $b_patients = qq{<center><button onclick="$cmd_resume" style='border-radius:3px;background-color:white;color:black;border: solid 1px black;'>$nb_samples</button></center>};
		push(@list_project, $b_patients);
		if ($quick_request) {
			if ($is_rna) { push(@list_project, 'RNA'); }
			else { push(@list_project, 'DNA'); }
		}
		else {
			push(@list_project, ucfirst($pro->{type_analyse}));
		}
		push(@list_project, $pro->{capture_name});
		push(@list_project, $pro->{genome});
		
		if ($quick_request) {
			push(@list_datas, \@list_project);
			next;
		}
		
		my ($url_query, $url_diag, $url_viewer, $url_resume);
		my ($b_query, $b_diag, $b_viewer, $b_resume);
		my $disabled;
		
		my $is_public = $buffer->getQuery->isProjectPublic($pro->{id});
		next if not $isUserMagic and not $is_public;

		if ($pro->{'cache'} eq 'rocks') {
			$url_query = 'https://www.polyweb.fr/polyweb/vector/gene.html?project='.$pro->{'name'};
			$url_diag = 'https://www.polyweb.fr/polyweb/coverage.html?project='.$pro->{'name'};
			$url_viewer = 'https://www.polyweb.fr/polyweb/polyviewer.html?project='.$pro->{'name'};
			
		}
		elsif ($pro->{'cache'} eq 'resume') {
			$disabled = 'disabled';
		}
		else {
			$url_query = 'https://'.$ENV{HTTP_HOST}."/".$buffer->config->{polyweb_url}->{polyweb_OLD}.'polyweb/vector/gene.html?project='.$pro->{'name'};
			$url_diag = 'https://'.$ENV{HTTP_HOST}."/".$buffer->config->{polyweb_url}->{polyweb_OLD}.'polyweb/coverage.html?project='.$pro->{'name'};
			$url_viewer = 'https://'.$ENV{HTTP_HOST}."/".$buffer->config->{polyweb_url}->{polyweb_OLD}.'polyweb/polyviewer.html?project='.$pro->{'name'};
		}
		if ($is_public) {
			#$b_resume .= qq{<button onclick="$cmd_resume" class='btn' style='background-color:white;color:orange;border:solid 1px orange;border-radius:3px;'>Resume</button>};
			if ($disabled) {
				$b_query .= qq{<b><i>--- NO INTERFACE ---</b></i>};
			}
			else {
				$b_query .= qq{<button onclick="window.open('$url_query', '_blank')"  class='btn' style="background-color:cornflowerblue;color:white;border:solid 1px cornflowerblue;border-radius:3px;" >PolyQuery</button>};
				if ($pro->{'type_analyse'} eq 'target') {
					$b_diag .= qq{<button onclick="window.open('$url_diag', '_blank')"  class='btn' style="background-color:orange;color:white;border:solid 1px orange;border-radius:3px;">PolyDiag</button>};
				}
				$b_viewer .= qq{<button onclick="window.open('$url_viewer', '_blank')"  class='btn' style="background-color:darkseagreen;color:white;border:solid 1px darkseagreen;border-radius:3px;">PolyViewer</button>};
			}
		}
		else {
			$b_query .= qq{<span style='color:red;'>--- OFFLINE ---</span>};
			$b_viewer .= qq{<button onclick="window.open('$url_viewer', '_blank')"  class='btn' style="background-color:darkseagreen;color:white;border:solid 1px darkseagreen;border-radius:3px;">PolyViewer</button>};
		}
		
		my @l_interfaces;
		#push(@l_interfaces, $b_resume);
		push(@l_interfaces, $b_query);
		push(@l_interfaces, $b_diag);
		push(@l_interfaces, $b_viewer);
		if ($pro->{cache} eq 'lmdb' and $is_public) { push(@l_interfaces, " <span style='color:red;font-size:9px;'><i><b>--OLD--</i></b></span>"); }
		push(@list_project, '<nobr><center>'.join('  ', @l_interfaces).'</center></nobr>');
		
		push(@list_project, $patients_text);
		push(@list_datas, \@list_project);
	}
	my $hRes;
	$hRes->{results} = \@list_datas;
	my $json_encode = encode_json $hRes;
	my $compressed_data = compress($json_encode);
	
	print "Access-Control-Allow-Origin: *\n";
	print "Access-Control-Allow-Methods: GET, POST, OPTIONS\n";
	print "Access-Control-Allow-Headers: Content-Type\n";
	print "Content-Encoding: deflate\n"; 
	print "Content-type: application/json; charset=utf-8\n";
	print "Content-Length: " . length($compressed_data) . "\n\n";
	print $compressed_data;
	exit(0);
}


