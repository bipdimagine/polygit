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


 
sub getProjectListsRNA {
	my ( $buffer, $login, $pwd, $project_name ) = @_;
	my $is_BIPD_login = $buffer->getQuery->isLoginSTAFF($login);
	my $list_proj = $buffer->getQuery->getListProjectsRnaSeq($project_name);
	
	my @list_res;
	foreach my $h (@$list_proj) {
		my $h_users;
		foreach my $this_user_name (split(',', $h->{username})) { $h_users->{$this_user_name} = undef; }
		unless ($is_BIPD_login) {
			next unless exists $h_users->{$login};
		}
		my $name = $h->{name};
		my $b1 = GBuffer->new;
		my $p1 = $b1->newProject( -name => $name );
		$h->{button} = '';
#		if (-d $p1->get_path_rna_seq_junctions_root()) {
#			my $project_name = $h->{name};
#			#my $path_function_ok = $p1->get_path_rna_seq_junctions_root().'/AllRes/';
#			#$project_name = 'disabled' if (not -d $path_function_ok);
#			$h->{button} = '1::'.$project_name;
#		}
#		$h->{button} = '2::'.$h->{name} if (-d $p1->get_path_rna_seq_polyrna_root());
		if (-d $p1->get_path_rna_seq_junctions_root()) {
			my $ok;
			foreach my $pat (@{$p1->getPatients()}) {
				eval { $pat->getJunctionsAnalysePath() };
				if ($@) { next; }
				$ok = 1 if ($pat->junction_RI_file_filtred());
				$ok = 1 if ($pat->junction_SE_file_filtred());
				$ok = 1 if ($pat->junction_RI_file());
				$ok = 1 if ($pat->junction_SE_file());
				last if $ok;
			}
			$h->{button} = '1::'.$h->{name} if ($ok);
		}
		$h->{button} = '2::'.$h->{name} if (-d $p1->get_path_rna_seq_polyrna_root());
		push (@list_res, $h);
	}
	export_data::print_simpleJson($cgi, \@list_res);
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
	