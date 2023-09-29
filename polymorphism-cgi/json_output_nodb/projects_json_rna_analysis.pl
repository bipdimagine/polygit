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
#use Text::CSV::Slurp; 

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

#si le paramètre type est egal à login alors un utilisateur cherhce à se connecter
#on recupère alors la valeur des paramètres login et pwd et on appelle la fonction getProjetLogin pour renvoyer la liste des projets associés à ce login et mdp
my $login = $cgi->param('login');
my $pwd   = $cgi->param('pwd');
my $project_query   = $cgi->param('project');
my $action = $cgi->param('action');
my $only_diag = $cgi->param('only_diag');
my $export_xls = $cgi->param('export_xls');
my $viewer = $cgi->param('viewer');
my $defidiag = $cgi->param('defidiag');
getProjectListsRNA($buffer,$login,$pwd,$project_query) if $action eq "list";

sub get_path_rna_seq_junctions_root {
	my ($h_proj) = @_;
	my $path = $buffer->getDataDirectory("root")."/ngs/".$h_proj->{name}."/".$h_proj->{relname}."/analysis/";
	return $path;
}

sub get_path_rna_seq_polyrna_root {
	my ($h_proj) = @_;
	my $path = $buffer->getDataDirectory("root")."/ngs/".$h_proj->{name}."/".$h_proj->{relname}."/polyRNA/";
	return $path if -d $path;
	my @lPotentialRelease = ('HG19', 'HG19_CNG', 'HG19_MT', 'HG38', 'HG38-ERCC', 'MM38', 'MM39');
	foreach my $rel2 (@lPotentialRelease) {
		my $alt_path = $buffer->getDataDirectory("root")."/ngs/".$h_proj->{name}."/".$rel2."/polyRNA/";
		return $alt_path if -d $alt_path;
	}
	return $path;
}
 
sub getProjectListsRNA {
	my ( $buffer, $login, $pwd, $project_name ) = @_;
	print $cgi->header('text/json-comment-filtered');
	print "{\"progress\":\".";
	
	my $is_BIPD_login = $buffer->getQuery->isLoginSTAFF($login);
	my $list_proj = $buffer->getQuery->getListProjectsRnaSeqFromLoginPwd($login, $pwd, $project_name);
	my ($hDone, $h_ok);
	
	
	
	
	foreach my $h (@$list_proj) {
		print '.';
		my $h_users;
		foreach my $this_user_name (split(',', $h->{username})) { $h_users->{$this_user_name} = undef; }
		unless ($is_BIPD_login) {
			next unless exists $h_users->{$login};
		}
		my $name = $h->{name};
		next if exists $hDone->{name};
		my $path = get_path_rna_seq_junctions_root($h);
		
		$h->{button} = '';
		my @lbuttons;
		my @lPatNames = split(';', $h->{patient_name});
		if (-d $path) {
			my $ok;
			foreach my $patient_name (@lPatNames) { 
				my $dragen_file = $path.'/../variations/dragen-sj/'.$patient_name.'.SJ.out.tab.gz';
				if (-e $dragen_file) {
					$ok = 1;
					last;
				}
			}
			my $se_file = $path.'/allResSE.txt';
			my $ri_file = $path.'/allResRI.txt';
			$ok = 1 if (-e $se_file);
			$ok = 1 if (-e $ri_file);
			$se_file = $path.'/AllRes/allResSE.txt';
			$ri_file = $path.'/AllRes/allResRI.txt';
			$ok = 1 if (-e $se_file);
			$ok = 1 if (-e $ri_file);
			push(@lbuttons, '1::'.$h->{name}) if ($ok);
			$hDone->{$name} = undef if $ok;
		}
		else {
			my $ok_dragen;
			foreach my $patient_name (@lPatNames) { 
				my $dragen_file = $path.'/../variations/dragen-sj/'.$patient_name.'.SJ.out.tab.gz';
				if (-e $dragen_file) {
					$ok_dragen = 1;
					last;
				}
			}
			push(@lbuttons, '1::'.$h->{name}) if ($ok_dragen);
			$hDone->{$name} = undef if $ok_dragen;
		}
		
		my $path_polyRNA = get_path_rna_seq_polyrna_root($h);
		if (-d $path_polyRNA) {
			push(@lbuttons, '2::'.$h->{name});			
			$hDone->{$name} = undef;
		}
		$h->{button} = join(';', @lbuttons);
		$hDone->{$name} = undef;
		next unless $h->{button};
		$h_ok->{$name} = $h;
	}
	my @list_res;
	foreach my $name (reverse sort keys %$h_ok) {
		push (@list_res, $h_ok->{$name});
	}
	print_simpleJson($cgi, \@list_res);
	exit(0);
}



sub print_simpleJson {
	my ($cgi,$data,$type_identifier) = @_;
	
	$type_identifier = "name" unless defined $type_identifier;
	my %all;
	#$all{identifier} = "name";
	$all{label}      = $type_identifier;
	$all{items}      = $data;	
	
	my $json_encode = encode_json \%all;
	print ".\",";
	$json_encode =~ s/{//;
	print $json_encode;
	exit(0);
}



sub checkAuthentification {
	my ( $buffer, $login, $pwd,$project ) = @_;
	my $res = GenBoProjectQueryNgs::getAuthentificationForUser($buffer->dbh,$project,$login,$pwd);
	GenBoProjectQueryNgs::writeLatestInfos($buffer->dbh,$project,$login,$pwd);
	die();
	
	my $items;
	$items->{name} = "BAD_LOGIN";
	$items->{name} = "OK" if $res>0;

	export_data::print_simpleJson($cgi,[$items]);
	
	exit(0);
	
}
	