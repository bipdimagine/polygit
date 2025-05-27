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
getProjectListsDefidiag($buffer,$login,$pwd) if $action eq "list" && !($project_query) ;
checkAuthentification($buffer,$login,$pwd,$cgi->param('project')) if $action eq "check";


 
sub getProjectListsDefidiag {
	my ( $buffer, $login, $pwd ) = @_;
	print $cgi->header('text/json-comment-filtered');
	print "{\"progress\":\".";
	
	my $res2 = $buffer->getQuery()->getProjectListForUser($login, $pwd );
	print '.';
	my $res = $res2;
	if ($defidiag){
		@$res = grep{$_->{validation_db}  eq "defidiag"} @$res2;
	}
	else {
		@$res = grep{$_->{validation_db}  ne "defidiag"} @$res2;
	}
	my $vquery = QueryValidationAcmg->new(dbh=>$buffer->dbh,database=>"ACMG");
	
	my @ids = map{$_->{id}} @$res;
	foreach my $pro (@$res) {
		print '.';
		if($defidiag) {
			#next if $pro->{validation_db}  ne "defidiag";
		}
		my $debug;
		$debug=1 if $pro->{name} eq "NGS2020_3165";
		my $pheno = $buffer->queryPhenotype->getPhenotypeInfosForProject($pro->{id});
		$pro->{gencode} = $buffer->getQuery->getGencodeVersion($pro->{id});
		$pro->{genome} = $buffer->getQuery->getReleaseGenome($pro->{id});
		$pro->{annotation} = $buffer->getQuery->getPublicDatabaseVersion($pro->{id});
	 	$pro->{build} = $buffer->getQuery()->getBuildFromProjectName($pro->{name});
		my $husers = $buffer->getQuery()->getOwnerProject($pro->{id});
		my $hstatus = $vquery->getLatestStatusProject($pro->{id});
		my $nb_validation = scalar(keys %$hstatus);
		
	#	die() if $hstatus;
	
		if ($pheno){
			$pro->{phenotype} = join(";",keys %$pheno);
		}
		else {
			$pro->{phenotype} = "-";
		}
		if (exists $pro->{group}){
			$pro->{sites} = lc($pro->{group});
		}
		else{
		my $usites;
		foreach my $u (@$husers) {
			foreach my $k (keys %$site){
				if ($u->{email} =~ /$k/) {
					$usites->{$site->{$k}}++;
				}
			}
		}
			$pro->{sites} = join(";",keys %$usites);
		}
	
		my ($latest_view,$since) = $buffer->getQuery->getLastConnectionUserProject($login,$pro->{id});
		my $gencode = $buffer->getQuery->getGencodeVersion($pro->{id});

		my @samples =  @{$buffer->getQuery->getPatients($pro->{id})};
		
		$pro->{samples} = scalar(@samples);
		$pro->{validations} = $nb_validation;
		my @trio = grep{$_->{status} == 2 && ($_->{mother} or $_->{father})} @samples;
		my $dir = $buffer->getDataDirectory("cache")."/".$pro->{genome}.'.'.$gencode.".".$pro->{annotation}."/".$pro->{name}."/vector/lmdb_cache/1";
		my $dir_polyviewer =   $buffer->getDataDirectory("cache")."/".$pro->{genome}.'.'.$gencode.".".$pro->{annotation}."/".$pro->{name}."/vector/lmdb_cache/1";
		my $dir_polyviewer_rocks =   $buffer->getDataDirectory("cache")."/rocks/".$pro->{genome}.'.'.$gencode.".".$pro->{annotation}."/".$pro->{name}."/polyviewer_objects.rocksdb/IDENTITY";
		if ($pro->{genome} eq 'MT') {
			$dir_polyviewer =   $buffer->getDataDirectory("cache")."/".$pro->{genome}.'.'.$gencode.".".$pro->{annotation}."/".$pro->{name}."/vector/lmdb_cache/MT";
		}
		$pro->{polyviewer} = 2;
		$pro->{polyviewer} = 1 if @samples && -e "$dir_polyviewer/".$samples[0]->{name}.".variants-genes.polyviewer";
		$pro->{polyviewer} = 1 if @samples && -e $dir_polyviewer_rocks;
		
		if (-e $dir_polyviewer_rocks and $gencode >= 46) {
			$pro->{genome} .= '-rocks';
		}
		
#		warn $pro->{polyviewer} if $pro->{name} eq "NGS2020_3165";
		my $t = (stat $dir )[9];
		my ($dc,$h) = split(" ", $pro->{creation_date});
		my $date = POSIX::strftime( 
             "%y/%m/%d", 
             localtime( 
               		$t
                 )
             );
             
		my $days_difference = int((time - $t) / 86400);
		
		
		my $vtrio;
		if (@trio){
		$pro->{child} = $trio[0]->{name};
		if (scalar(@trio) > 1){
		$pro->{child} = scalar @trio;
		}
		my $hv = $vquery->getValidationPatientName($pro->{child},$pro->{name});
		
		map {$pro->{patient_name}.=$_->{name}.";".$_->{family}.";"} @samples;
		
		$pro->{child_sex} = $trio[0]->{sex};
		$pro->{trio} = "trio";
		$vtrio = scalar(@trio);
		}
		else {
			$pro->{trio} = "solo";
			my @children = grep{$_->{status} == 2} @samples;
			@children = @samples unless (@children);
			$pro->{child} = $children[0]->{name};
			$pro->{child_sex} = $children[0]->{sex};
			$vtrio = 0;
			map {$pro->{patient_name}.=$_->{name}.";".$_->{family}.";"} @samples;
		}

		my @t = sort{$a <=> $b} keys %{$buffer->public_data};
		my $latest = $t[-1]; 
		$pro->{users} = join("\n",map{$_->{email}} @$husers);
		$pro->{users} = scalar(@$husers)." users " if scalar(@$husers) > 6;	
		
		
		$pro->{version} = abs($latest-$pro->{annotation})."::".$pro->{gencode}."-".$pro->{annotation};
		$pro->{creation_date} = $dc;
		$pro->{description} =$pro->{description};
		$pro->{annotation_date} =localtime($t)->ymd;
		$pro->{annotation_since} = $days_difference;
		$pro->{name} = $pro->{name};
		$pro->{trio} = $pro->{trio};
		$pro->{b1} = $pro->{name};
		$pro->{b2} = $pro->{name};
		#if ($pro->{child}){
			my $hv = $vquery->getValidationProjectName($pro->{name});
			$pro->{acmg} = "*";
			if (scalar(@$hv)) {
			   my $text = $buffer->getValidationTerm($hv->[0]->{validation});
			   my $toto;
			   ($hv->[0]->{modification_date},$toto) = split(" ",$hv->[0]->{modification_date});
			   
			   $pro->{acmg} = "$text"."_".$hv->[0]->{modification_date}."_".$hv->[0]->{user_name}."_".$hv->[0]->{validation};
			}
			
		#}
		if ($latest_view){
		my $h;
		($pro->{latest_view},$h) =split(" ",$latest_view);
		$pro->{latest_view} = $pro->{latest_view};
		
		$pro->{since} =$since;
		}
	
		$pro->{new_pathogenic} ="";
			my $list_resume = get_values_new_hgmd($buffer, $pro->{name});
		my $versions = $buffer->config->{polybtf_default_releases}->{default};
		my @lVersions = (split(',', $versions));
		my $i = 0;
	 	foreach my $resume_new_pathogenic (@$list_resume) {
	 		if ( $pro->{new_pathogenic}) { $pro->{new_pathogenic} .= ';'; }
			my ($nb_red, $nb_coral, $nb_orange, $nb_yellow, $nb_grey) = split(';', $resume_new_pathogenic);
			if ($resume_new_pathogenic eq '?') { $pro->{new_pathogenic} .= 'F'; }
			elsif ($nb_red)    { $pro->{new_pathogenic} .= 'A::'.$pro->{name}.'::New+++'; }
			elsif ($nb_coral)  { $pro->{new_pathogenic} .= 'B::'.$pro->{name}.'::New++'; }
			elsif ($nb_orange) { $pro->{new_pathogenic} .= 'C::'.$pro->{name}.'::New+'; }
			elsif ($nb_yellow) { $pro->{new_pathogenic} .= 'D::'.$pro->{name}.'::New+'; }
			elsif ($nb_grey)   { $pro->{new_pathogenic} .= 'E::'.$pro->{name}.'::New'; }
			else               { $pro->{new_pathogenic} .= 'F'; }
			$pro->{new_pathogenic} .= '::'.$lVersions[$i];
			$i++;
	 	}
		
		my $patients = 	$buffer->getQuery()->getPatients($pro->{id});
		my %captures_id;
		foreach my $p (@$patients){
			$captures_id{$p->{capture_id}} ++;
		}
		my $find;
		my @c;
		foreach my $cid (keys %captures_id){
			my $capt =  $buffer->getQuery()->getCaptureInfos($cid);
			push(@c,$capt->{name});
		}
		$pro->{capture_name} = join(",",@c);
		$pro->{button} = $pro->{polyviewer}."::".$pro->{name}.'::'.$pro->{genome};

		
	}
	
	my @res2 = sort{$a->{annotation_since} <=> $b->{annotation_since}} @$res;
	print_simpleJson($cgi,\@res2);
	
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
	my $coverage_dir = $buffer->{config}->{project_data}->{root}."/". $buffer->{config}->{project_data}->{ngs}."/$projectName/$build/align/coverage/";
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
	