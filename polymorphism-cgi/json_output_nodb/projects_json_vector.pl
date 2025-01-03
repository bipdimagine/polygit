#!/usr/bin/perl
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";


use connect;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use export_data;
use Carp;
use JSON;
use Digest::MD5 qw(md5_hex);

#ce script est utilisé pour renvoyer les noms des projets d'un utilisateur après sa connection à l'interface ainsi que pour telecharger des fichiers via l'interface ou encore exporter les tableaux de données au format xls

my $cgi    = new CGI;

#chargement du buffer 
my $buffer = GBuffer->new;


#si le paramètre type est egal à login alors un utilisateur cherhce à se connecter
#on recupère alors la valeur des paramètres login et pwd et on appelle la fonction getProjetLogin pour renvoyer la liste des projets associés à ce login et mdp
my $login = $cgi->param('login');
my $pwd   = $cgi->param('pwd');
my $project_query   = $cgi->param('project');
my $action = $cgi->param('action');
my $only_diag = $cgi->param('only_diag');
getProjectLists($buffer,$login,$pwd) if $action eq "list";
checkAuthentification($buffer,$login,$pwd,$project_query) if $action eq "check";
checkMD5Files($buffer, $project_query) if $action eq 'checkMD5';

sub getProjectLists {
	my ( $buffer, $login, $pwd ) = @_;
	#renvoit tous les origin associés au login et au mdp spécifié
	
	my $res = $buffer->getQuery()->getProjectListForUser($login, $pwd );
	my $db_name = $buffer->{config}->{server}->{status};
	
	my $nb_project = scalar(@{$buffer->listProjects()});
	#pour chacun des id de projet 
	if ($project_query){
	 foreach my $project (@$res) {
	 	next if ($project_query && $project_query ne $project->{name});
	 	my $buffer2 = GBuffer->new;
	 	my $project2 = $buffer2->newProject( -name => $project->{name} );
	 	my $captures = $project2->getCaptures();
	 	my $runs = $project2->getRuns();
	 	$project->{nb_run} = scalar(@$runs);
	 	$project->{capture} = join(" , ",map{$_->name}@$captures); 

	 	#my $husers = $buffer->getQuery()->getOwnerProject($project->{id});
	 	my $ph = $buffer->getQuery()->getProjectByName($project->{name});
	 	my $sm =  $buffer->getQuery()->getSequencingMachines($project->{id});
	 	
	    my $sm_string = join(";",@$sm);
		
		$project->{genome_version} = $ph->{version};#join(",",map{$_->{name}}@$release);
		$project->{db_name} = $db_name; 
		my ($t,$h) = split("_",$project->{dbname});
		$project->{nb_total} = $nb_project; 
		$project->{pedigree} = 0;
		$project->{somatic} = 0;
		$project->{pedigree} = 1 if ($project2->pedigree_details());
		$project->{somatic} = 1 if ($project2->somatic_details());
			$project->{sequencing} = $sm_string;
		if (lc($project->{dbname}) eq "polyexome" || lc($project->{dbname}) eq "polyrock"){
			$project->{dbname} = $project->{dbname}."_".$ph->{version};
		}
		
	 }
	export_data::print_simpleJson($cgi,$res) if $project_query;
	}
	my @res2;
	foreach my $project (@$res) {

		my $husers = $buffer->getQuery()->getOwnerProject($project->{id});
		my $ph = $buffer->getQuery()->getProjectByName($project->{name});
		
#	 	my $captures = $project2->getCaptures();
#	 	next unless @{$captures->{0}->transcripts_name};
		
		$project->{genome_version} = $ph->{version};#join(",",map{$_->{name}}@$release);
		$project->{db_name} = $db_name; 
		if (lc($project->{dbname}) eq "polyexome" || lc($project->{dbname}) eq "polyrock"){
			$project->{dbname} = $project->{dbname}."_".$ph->{version};
		}
		my $nb_p = 0;
		my ($name_database,$junk) = split("_",$project->{dbname});
#		unless (exists $database_nb->{$name_database}){
#			my $sql = qq{select p.name as NAME from PolyprojectNGS.polydb pb , PolyprojectNGS.databases_projects dp , PolyprojectNGS.projects p 
#	 						where pb.name="$name_database" 
#							and dp.db_id=pb.db_id and p.project_id=dp.project_id and p.type_project_id=3};
#			
#			my @pp = keys %{connect::return_hashref($buffer->dbh,$sql,"NAME")};
#			
#			$database_nb->{$name_database}->{NB} = scalar(@pp); 
#			
#			
#		}
		$project->{nb_total} = 0;#$database_nb->{$name_database}->{NB};
		$project->{users} = join("\n",map{$_->{email}} @$husers);		
		my $patients = 	$buffer->getQuery()->getPatients($project->{id});
		my %captures_id;
		foreach my $p (@$patients){
			$captures_id{$p->{capture_id}} ++;
		}
		my $find;
		my @c;
		my @t;
		my $nb =0;
		foreach my $cid (keys %captures_id){
			my $capt =  $buffer->getQuery()->getCaptureInfos($cid);
			push(@c,$capt->{name});
			push(@t,$capt->{type});
			$find=1 if $capt->{analyse} ne "exome";
			#$find=1 if $capt->{transcripts} ne "";
			$nb += split(";",$capt->{transcripts} );
			
			
		}
		if ($only_diag){
			next unless $find;
			$project->{capture_type} = join(",",@t);
			$project->{nb_transcripts} = join(",",$nb);
			my %machines;
			my $runs = $buffer->getQuery()->getRuns($project->{id});
			$project->{nb_runs} = scalar(@$runs);
			foreach my $rid (@$runs){
				my $run =  $buffer->getQuery()->getRunInfosById($rid->{id});
				$machines{$run->[0]->{machine} } ++;
			}
			$project->{machines} = join(",",keys %machines);
			
		}
		$project->{capture_name} = join(",",@c);
		
		#on compte un certains nombres d'objets pour afficher des infos générelaes relatives au projet
		map {$project->{patient_name}.=$_->{name}.";"} @$patients;
		$project->{patient}  = scalar(@$patients);#GenBoProjectQueryNgs::countObjects($buffer->dbh(),$project->{dbname},$buffer->getType("patient")->id,$project->{id});
		$project->{substitution}  = 0;#GenBoProjectQueryNgs::countObjects($buffer->dbh(),$project->{dbname},$buffer->getType("variation")->id,$project->{id});
		$project->{insertion}  = 0;#GenBoProjectQueryNgs::countObjects($buffer->dbh(),$project->{dbname},$buffer->getType("insertion")->id,$project->{id});
		$project->{deletion}   = 0;#GenBoProjectQueryNgs::countObjects($buffer->dbh(),$project->{dbname},$buffer->getType("deletion")->id,$project->{id});
		#$items{reference}     = 0;#GenBoProjectQuery::countObjects($buffer->dbh(),$buffer->getType("reference")->id,$var->{id});
		$project->{coverage} = 0;#getCoverage($project);
		push(@res2,$project);
		
		
	}
	
	export_data::print_simpleJson($cgi,\@res2);
	
	exit(0);
}

sub getCoverage {
	my ($project) = @_;
	my $projectName = $project->name();
	my $build = $project->getVersion();
	my $coverage_dir = $buffer->config_path("project_data")."/". $buffer->{config}->{project_data}->{ngs}."/$projectName/$build/align/coverage/";
	my @files = `ls $coverage_dir/*.cov.gz 2>/dev/null`;
	chomp (@files);
	
	my $tabix = $buffer->{config}->{software}->{tabix};;
	my $nb;
	my $total = 0;
	foreach my $f (@files){
		$nb++;
		my ($data) = grep {$_ =~ /mean_all\t99/} `$tabix  $f mean_all`;
		my @t = split(" ",$data);
		$total += $t[-1];	
	}
	return -1 if $nb == 0;
	
	my $rr = int $total/$nb;
	
	return ($rr);
	
}

sub checkAuthentification {
	my ( $buffer, $login, $pwd,$project ) = @_;
	my $res = GenBoProjectQueryNgs::getAuthentificationForUser($buffer->dbh,$project,$login,$pwd);
	my $items;
	$items->{name} = "BAD_LOGIN";
	$items->{name} = "OK" if $res>0;
	export_data::print_simpleJson($cgi,[$items]);
	exit(0);
}

sub checkMD5Files {
	my ($buffer, $project_name) = @_;
	my @lItems;
	my $project = $buffer->newProjectCache( -name 	=> $project_name,
										    -cache 	=> '1', );
	my $hInfos = $project->global_infos->{'check'};
	foreach my $var_ind (sort keys %{$hInfos->{vcf}}) {
		foreach my $pat_name (sort keys %{$hInfos->{vcf}->{$var_ind}}) {
			my $j = 1;
			foreach my $file (sort keys %{$hInfos->{vcf}->{$var_ind}->{$pat_name}}) {
				my $hash;
				$hash->{file} = $file;
				unless (-e $file) {
					$hash->{status} = 'ERROR';
					push(@lItems, $hash);
				}
				else {
					my $md5_kct = $hInfos->{vcf}->{$var_ind}->{$pat_name}->{$file};
					my $md5_now = md5_hex($file);
					if ($md5_kct eq $md5_now) { $hash->{status} = 'OK'; }
					else {
						$hash->{status} = 'ERROR';
						push(@lItems, $hash);
					}
				}
			}
		}
	}
	my $hash;
	$hash->{label} = 'file';
	$hash->{name} = 'file';
	$hash->{items} = \@lItems;
	export_data::print_json($cgi, $hash);
	exit(0);
}
	