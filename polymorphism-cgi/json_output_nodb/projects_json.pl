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

use connect;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use export_data;
use Carp;
use JSON;
#use Text::CSV::Slurp; 

#ce script est utilisé pour renvoyer les noms des projets d'un utilisateur après sa connection à l'interface ainsi que pour telecharger des fichiers via l'interface ou encore exporter les tableaux de données au format xls

my $cgi    = new CGI;

#chargement du buffer 
my $buffer = GBuffer->new;


#si le paramètre type est egal à login alors un utilisateur cherhce à se connecter
#on recupère alors la valeur des paramètres login et pwd et on appelle la fonction getProjetLogin pour renvoyer la liste des projets associés à ce login et mdp
my $login = $cgi->param('login');
my $pwd   = $cgi->param('pwd');
my $project_query   = $cgi->param('project');
my $project_name   = $cgi->param('project_name');
my $action = $cgi->param('action');
my $only_diag = $cgi->param('only_diag');
my $export_xls = $cgi->param('export_xls');
my $viewer = $cgi->param('viewer');
my $mycode   = $cgi->param('mycode');
checkOtpAuth($buffer, $mycode, $login) if $action eq "otp";
getProjectLists($buffer,$login,$pwd) if $action eq "list";
checkAuthentification($buffer,$login,$pwd,$cgi->param('project')) if $action eq "check";
checkHgmdAccess($buffer,$login,$pwd) if $action eq "check_hgmd_access";



sub checkOtpAuth {
	my ( $buffer, $mycode, $login ) = @_;
	my $hRes;
	if ($login =~ /[a-zA-Z]/) {
		my $need_otp = $buffer->use_otp_for_login($login);
		if ($need_otp) {
			my $secret_code = $buffer->google_auth_secret_pwd($login);
			my $auth = $buffer->google_auth();
			$auth->secret32( $secret_code );
			if ($auth->verify($mycode)) { $hRes->{otp_verif} = 'yes'; }
			else { $hRes->{otp_verif} =  'no'; }
		}
		else {
			$hRes->{otp_verif} = 'yes';
		}
	}
	else {
		$hRes->{otp_verif} = 0;
	}
	my $json_encode = encode_json $hRes;
	print $cgi->header('text/json-comment-filtered');
	print $json_encode;
	exit(0);
}

sub checkHgmdAccess {
	my ( $buffer, $login, $pwd ) = @_;
	my $hRes->{hgmd_access} = $buffer->hasHgmdAccess($login);
	my $json_encode = encode_json $hRes;
	print $cgi->header('text/json-comment-filtered');
	print $json_encode;
	exit(0);
}

sub getProjectLists {
	my ( $buffer, $login, $pwd ) = @_;
	#renvoit tous les origin associés au login et au mdp spécifié
	
	my $res = $buffer->getQuery()->getProjectListForUser($login, $pwd );
	if ($res->[0]->{otp} == 1 ) {
		#$res->[0]->{otp_qr_code} = $buffer->google_auth_qr_code();
		$res->[0]->{otp_code} = $buffer->google_auth_issuer().': '.$buffer->google_auth_key_id();
		#$res->[0]->{otp_need} == 1
		
	}
	$res->[0]->{otp_need} = $buffer->use_otp_for_login($login);
	my $db_name = $buffer->{config}->{server}->{status};
	
	my $nb_project = scalar(@{$buffer->listProjects()});
	
	my $h_proj_user;
	foreach my $project ( @$res) {
		$h_proj_user->{$project->{name}}->{high} = 0;
		$h_proj_user->{$project->{name}}->{medium} = 0;
		$h_proj_user->{$project->{name}}->{others} = 0;
		$h_proj_user->{$project->{name}}->{all} = 0;
	}
	#pour chacun des id de projet 
	if ($project_query){
	 foreach my $project ( @$res) {
	 	next if ($project_query && $project_query ne $project->{name}); 		
	 	
	 	my $buffer2 = GBuffer->new;
	 	my $project2 = $buffer2->newProject( -name => $project->{name} );
	 	 	
	 	my $captures = $project2->getCaptures();
	 	my $runs = $project2->getRuns();
	 	$project->{nb_run} = scalar(@$runs);
	 	$project->{capture} = join(" , ",map{$_->name} @$captures); 
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
		$project->{pedigree} = 1 if -e $buffer->config_path("root","project_data")."/".$project->{type}."/".$project->{name}."/".$project->{name}.".ped";
		$project->{somatic} = 1 if -e $buffer->config_path("root","project_data")."/".$project->{type}."/".$project->{name}."/".$project->{name}.".somatic";
			$project->{sequencing} = $sm_string;
		if (lc($project->{dbname}) eq "polyexome" || lc($project->{dbname}) eq "polyrock"){
			$project->{dbname} = $project->{dbname}."_".$ph->{version};
		}
		$project->{type_seq} = "diagnostic";
		#$project->{type} = ;
	
		$project->{isDude} = 0;
		$project->{isDude} = 1 if ($project2->isDude());
		#warn 
		unless ($buffer->getQuery()->isUserMagic($login)){
		$buffer->getQuery()->writeLatestInfos($project2->id,$login);
		$viewer ="polyDiag" unless  $viewer;
		$buffer->getQuery->updateLastConnectionUserProject($login, $project2->id(),$viewer);# unless ($test or $get_bundles or $check_genes);
		}
	 }
	 if ($project_name){
	 	my $buffer2 = GBuffer->new;
	 	my $project2 = $buffer2->newProject( -name => $project_name ); 	
	 	my $captures = $project2->getCaptures();
	 	foreach my $a (@$res){
	 		next unless $a->{name} ne $project_name;
	 		$a->{sequencing_type} = "panel";
	 		$a->{sequencing_type} = "genome" if $project2->isGenome;
	 		$a->{sequencing_type} = "exome" if $project2->isExome;
	 		
	 	}
	 }
	export_data::print_simpleJson($cgi,$res) if $project_query;
	}
	my (@res2, $h_proj_new_path);
	my @test;
	my @test2;
	foreach my $project (sort {$b->{name} cmp $a->{name}} @$res) {
		next unless	 $project->{name} =~ /NGS/;
		push(@test, "\"".$project->{name}."\"");
		push(@test2,$project->{id});
	}
	 my $phname = $buffer->getQuery()->fast_getProjectByName(join(",",@test));
	 my $phpatients = $buffer->getQuery()->fast_getPatients(join(",",@test2));
	# warn Dumper $phpatients;
	 my $cptinfos;
	foreach my $project (sort {$b->{name} cmp $a->{name}} @$res) {
		next unless	 $project->{name} =~ /NGS/;
		my $husers = $buffer->getQuery()->getOwnerProject($project->{id});
		my $groups = $buffer->getQuery()->getUserGroupsForProject($project->{id});
	 	$project->{genome} = $buffer->getQuery()->getBuildFromProjectName($project->{name});

		$project->{users} = "";
		if (@$groups){
			$project->{users} = join(";",@$groups);
		
		}
		else {
			$project->{users} .= join("\n",map{$_->{email}} @$husers);
		$project->{users} .= scalar(@$husers)." users " if scalar(@$husers) > 6;	
		}
		
		my $ph;
		 $ph =  $phname->{$project->{name}};
	
		my $patients =$phpatients->{$project->{id}};
	
		$project->{genome_version} = $ph->{version};#join(",",map{$_->{name}}@$release);
		
		my $nb_p = 0;
		$project->{nb_total} = 0;#$database_nb->{$name_database}->{NB};
		
		if ($project->{users} eq "") {
			$project->{users} = qq{--??!!!!!!!!??--};
		}
		my %captures_id;
		foreach my $p (@$patients){
			$captures_id{$p->{capture_id}} ++;
		}
		my $find;
		my @c;
		my @t;
		my $nb =0;
	
		foreach my $cid (keys %captures_id){
			my $capt;
			if (exists $cptinfos->{$cid}){
				$capt = $cptinfos->{$cid};
			}
			else {
			 $capt =  $buffer->getQuery()->getCaptureInfos($cid);
			}
			push(@c,$capt->{name});
			push(@t,$capt->{type});
			$find=1 if $capt->{analyse} ne "exome" && $capt->{analyse} ne "genome";
			#$find=1 if $capt->{analyse} ne "other";
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
		$project->{patient}  = scalar(@$patients);
		 $project->{new_pathogenic} ="-";
#		next;
#	 	my $list_resume = get_values_new_hgmd($buffer, $project->{name});
#		my $versions = $buffer->config->{polybtf_default_releases}->{default};
#		my @lVersions = (split(',', $versions));
#		my $i = 0;
#	 	foreach my $resume_new_pathogenic (@$list_resume) {
#	 		if (exists $project->{new_pathogenic}) { $project->{new_pathogenic} .= ';'; }
#			my ($nb_red, $nb_coral, $nb_orange, $nb_yellow, $nb_grey) = split(';', $resume_new_pathogenic);
#			if ($resume_new_pathogenic eq '?') { $project->{new_pathogenic} .= 'F'; }
#			elsif ($nb_red)    { $project->{new_pathogenic} .= 'A::'.$project->{name}.'::New+++'; }
#			elsif ($nb_coral)  { $project->{new_pathogenic} .= 'B::'.$project->{name}.'::New++'; }
#			elsif ($nb_orange) { $project->{new_pathogenic} .= 'C::'.$project->{name}.'::New+'; }
#			elsif ($nb_yellow) { $project->{new_pathogenic} .= 'D::'.$project->{name}.'::New+'; }
#			elsif ($nb_grey)   { $project->{new_pathogenic} .= 'E::'.$project->{name}.'::New'; }
#			else               { $project->{new_pathogenic} .= 'F'; }
#			$project->{new_pathogenic} .= '::'.$lVersions[$i];
#			$i++;
#	 	}
		push(@res2, $project);
		
	}
	if ($export_xls) { 
 		export_xls($res);
 		exit(0);
	}
	
	export_data::print_simpleJson($cgi,\@res2);
	
	exit(0);
}

sub get_values_new_hgmd {
	my ($buffer, $project_name) = @_;
	my @lValues;
	my $versions = $buffer->config->{polybtf_default_releases}->{default};
	foreach my $version (split(',', $versions)) {
		my $h = $buffer->get_polybtf_project_resume($project_name, $version);
		unless ($h) {
			push(@lValues, '?');
			next;
		}
		my $color_red = '#CE0000';
		my $color_coral = 'coral';
		my $color_orange = '#EFB73E';
		my $color_yellow = 'yellow';
		my $values;
		if ($buffer->hasHgmdAccess($login)) {
			$values = $h->{nb_red};
			$values .= ';'.$h->{nb_coral};
			$values .= ';'.$h->{nb_orange};
			$values .= ';'.$h->{nb_yellow};
			$values .= ';'.$h->{nb_grey};
		}
		else {
			$values = $h->{nb_red_no_hgmd};
			$values .= ';'.$h->{nb_coral_no_hgmd};
			$values .= ';'.$h->{nb_orange_no_hgmd};
			$values .= ';'.$h->{nb_yellow_no_hgmd};
			$values .= ';'.$h->{nb_grey_no_hgmd};
		}
		push(@lValues, $values);
	}
	return \@lValues;
}

sub export_xls {
	my ($res) = @_;
	$| = 1;
	print "Content-type: application/msexcel\n";
	print "Content-Disposition: attachment;filename=polyweb_list_projects.xls\n\n";
	my $workbook  = Spreadsheet::WriteExcel->new(\*STDOUT);
	my $worksheet = $workbook->add_worksheet('List Projects');
	my @header = ("Name", "Description", "Capture", "Capture Type", "Nb Runs", "Genes", "Nb Patient(s)", "Patient(s) name(s)", "Sequencers", "Users");
	my $nb_col = scalar(@header);
	my $format = $workbook->add_format( valign => 'vcentre', align => 'centre' );
	my $format_header = $workbook->add_format( valign => 'vcentre', align => 'centre' );
	$format_header->set_bold();
	$worksheet->write( 0, 0, \@header, $format_header );
	my $i = 1;
	foreach my $h_proj (@$res) {
		my @lToPrint;
		$lToPrint[0] = $h_proj->{name} if (exists $h_proj->{name});
		$lToPrint[1] = $h_proj->{description} if (exists $h_proj->{description});
		$lToPrint[2] = $h_proj->{capture_name} if (exists $h_proj->{capture_name});
		$lToPrint[3] = $h_proj->{capture_type} if (exists $h_proj->{capture_type});
		$lToPrint[4] = $h_proj->{nb_runs} if (exists $h_proj->{nb_runs});
		$lToPrint[5] = $h_proj->{genes} if (exists $h_proj->{genes});
		$lToPrint[6] = $h_proj->{patient} if (exists $h_proj->{patient});
		if (exists $h_proj->{patient_name}) {
			$h_proj->{patient_name} =~ s/;/, /g;
			$lToPrint[7] = $h_proj->{patient_name};
		}
		$lToPrint[8] = $h_proj->{sequencers} if (exists $h_proj->{sequencers});
		$lToPrint[9] = $h_proj->{users} if (exists $h_proj->{users});
		$worksheet->write( $i, 0, \@lToPrint, $format );
		$i++;
	}
 	$workbook->close();
	exit(0);
}

sub getCoverage {
	my ($project) = @_;
	my $projectName = $project->name();
	my $build = $project->getVersion();
	my $coverage_dir = $buffer->config_path("root","project_data")."/". $buffer->{config}->{project_data}->{ngs}."/$projectName/$build/align/coverage/";
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
	GenBoProjectQueryNgs::writeLatestInfos($buffer->dbh,$project,$login,$pwd);
	die();
	
	my $items;
	$items->{name} = "BAD_LOGIN";
	$items->{name} = "OK" if $res>0;

	export_data::print_simpleJson($cgi,[$items]);
	
	exit(0);
	
}
	