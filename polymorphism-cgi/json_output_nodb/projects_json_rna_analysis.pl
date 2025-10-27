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
		
		my $b1 = GBuffer->new;
		my $p1 = $b1->newProject( -name => $name );
		$h->{button_splices} = '';
		$h->{button_rna} = '';
		$h->{link_analyse} = '';
#		if (-d $p1->get_path_rna_seq_junctions_root()) {
#			my $project_name = $h->{name};
#			#my $path_function_ok = $p1->get_path_rna_seq_junctions_root().'/AllRes/';
#			#$project_name = 'disabled' if (not -d $path_function_ok);
#			$h->{button} = '1::'.$project_name;
#		}
#		$h->{button} = '2::'.$h->{name} if (-d $p1->get_path_rna_seq_polyrna_root());
		my $path = $p1->get_path_rna_seq_junctions_analyse_all_res();
		my (@lbuttons_splices, @lbuttons_rna);
		
		my $type_cache = $h->{'cache'};
		
		if (-d $path or -d $p1->getJunctionsDir('regtools')) {
			my $ok;
			foreach my $patient (@{$p1->getPatients()}) { 
				my $dragen_file = $p1->getVariationsDir('dragen-sj').'/'.$patient->name().'.SJ.out.tab.gz';
				$ok = 1 if (-e $dragen_file);
#				$ok = 1 if ($patient->regtools_file() and -e $patient->regtools_file());
			}
			my $se_file = $path.'/allResSE.txt' if (-e $path.'/allResSE.txt');
			$se_file = $path.'/allResSE.txt.gz' if (-e $path.'/allResSE.txt.gz');
			my $ri_file = $path.'/allResRI.txt' if (-e $path.'/allResRI.txt');
			$ri_file = $path.'/allResRI.txt.gz' if (-e $path.'/allResRI.txt.gz');
			$ok = 1 if (-e $se_file);
			$ok = 1 if (-e $ri_file);
			$ok = 1 if (-e $path.'/regtools/'.$p1->name().'_regtools.tsv.gz');
			push(@lbuttons_splices, '1::'.$h->{name}.'::'.$type_cache) if ($ok);
			$hDone->{$name} = undef if $ok;
		}
#		else {
#			my $ok_dragen;
#			foreach my $patient (@{$p1->getPatients()}) { 
#				my $dragen_file = $p1->getVariationsDir('dragen-sj').'/'.$patient->name().'.SJ.out.tab.gz';
#				$ok_dragen = 1 if (-e $dragen_file);
#			}
#			push(@lbuttons_splices, '1::'.$h->{name}) if ($ok_dragen);
#			$hDone->{$name} = undef if $ok_dragen;
#		}
		
		my $path_polyrna = $p1->get_path_rna_seq_polyrna_root();
		$path_polyrna =~ s/HG19_MT/HG19/;
		$path_polyrna =~ s/HG19_CNG/HG19/;
		$path_polyrna =~ s/HG19/HG38/;
		my $path_polyrna_HG38 = $path_polyrna;
		$path_polyrna_HG38 =~ s/HG38/HG38_CNG/;
		my $path_polyrna_HG38_CNG = $path_polyrna_HG38;
		if (-d $p1->get_path_rna_seq_polyrna_root() or -d $path_polyrna_HG38 or -d $path_polyrna_HG38_CNG) {
			push(@lbuttons_rna, '2::'.$h->{name});			
			$hDone->{$name} = undef;
		}
		
		my @lAnalysis;
		if (-d $p1->get_path_rna_seq_junctions_root()) {
			$hDone->{$name} = undef;
			
			my $path_DB = $p1->get_path_rna_seq_junctions_root();
			my @lPotentialPath;
			push(@lPotentialPath, $path_DB);
			
			
			if ($path_DB =~ /HG19/) {
				push(@lPotentialPath, $p1->buffer()->config_path("root","project_data")."/".$p1->getProjectType()."/".$p1->name()."/HG38/analysis/");
				push(@lPotentialPath, $p1->buffer()->config_path("root","project_data")."/".$p1->getProjectType()."/".$p1->name()."/HG38_CNG/analysis/");
			}
			elsif ($path_DB =~ /HG38/) {
				push(@lPotentialPath, $p1->buffer()->config_path("root","project_data")."/".$p1->getProjectType()."/".$p1->name()."/HG19/analysis/");
				push(@lPotentialPath, $p1->buffer()->config_path("root","project_data")."/".$p1->getProjectType()."/".$p1->name()."/HG19_MT/analysis/");
				push(@lPotentialPath, $p1->buffer()->config_path("root","project_data")."/".$p1->getProjectType()."/".$p1->name()."/HG19_CNG/analysis/");
			}
			elsif ($path_DB =~ /MM38/) {
				push(@lPotentialPath, $p1->buffer()->config_path("root","project_data")."/".$p1->getProjectType()."/".$p1->name()."/MM39/analysis/");
			}
			elsif ($path_DB =~ /MM39/) {
				push(@lPotentialPath, $p1->buffer()->config_path("root","project_data")."/".$p1->getProjectType()."/".$p1->name()."/MM38/analysis/");
			}
			foreach my $path (@lPotentialPath) {
				next unless -d $path;
				opendir my $dir, $path or die "Cannot open directory: $!";
				my @files = readdir $dir;
				closedir $dir;
				foreach my $path_2 (sort @files) {
					my $new_path = $path.'/'.$path_2;
					if (-d $new_path) {
						my $link = $new_path;
						$link =~ s/.+ngs/\/NGS/; 
						my $name_analysis = $path_2;
						$name_analysis =~ s/_NGS20.+//;
						push(@lAnalysis, '3::'.$name_analysis.'::'.$link.'/index.html') if (-e $new_path.'/index.html');
					}
				}
			}
		}
		
		my @ltmp = split(' ', $h->{creation_date});
		$h->{creation_date} = $ltmp[0];
		$h->{button_splices} = join(';', @lbuttons_splices);
		$h->{button_rna} = join(';', @lbuttons_rna);
		$h->{link_analyse} = $lAnalysis[-1] if (scalar @lAnalysis > 0);
		$hDone->{$name} = undef;
		
		my $is_ok;
		$is_ok = 1 if $h->{button_splices};
		$is_ok = 1 if $h->{button_rna};
		$is_ok = 1 if $h->{link_analyse};
		next unless $is_ok;
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
	