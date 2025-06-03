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
	my @lHTML;
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
		my (@c, $h_type);
		foreach my $cid (keys %captures_id){
			my $capt =  $buffer->getQuery()->getCaptureInfos($cid);
			push(@c,$capt->{name});
			$h_type->{$capt->{analyse}} = undef;
		}
		$pro->{capture_name} = join(",",@c);
		$pro->{button} = $pro->{polyviewer}."::".$pro->{name}.'::'.$pro->{genome};

		$pro->{'type_analyse'} = join(', ', sort keys %$h_type);
		next if $pro->{'type_analyse'} eq 'rnaseq';
		next if $pro->{'type_analyse'} eq 'singlecell';
		$pro->{'type_analyse'} = 'target' if $pro->{'type_analyse'} ne "exome" && $pro->{'type_analyse'} ne "genome";

		my $out = "<td>".$pro->{'b1'}."</td>";
		$out .= "<td>".$pro->{'description'}."</td>";
		my @lPat = split(';', $pro->{'patient_name'});
		my $patient_text = join(", ", sort @lPat);
		$out .= qq{<td><center><button onClick="alert('Patients: $patient_text');">}.$pro->{'samples'}."</button></center></td>"; #TODO: ajouter bouton
		$out .= "<td>".$pro->{'type_analyse'}."</td>";
		$out .= "<td>".$pro->{'capture_name'}."</td>";
		$out .= "<td>".$pro->{'phenotype'}."</td>";
		$out .= "<td><center>".$pro->{'genome'}."</center></td>";
		$out .= "<td><center>".$pro->{'gencode'}.'.'.$pro->{'annotation'}."</center></td>";
		
		my $base_url;
		$base_url = '//'.$buffer->config->{polyweb_url}->{polyweb_rocks_extension} if ($pro->{'cache'} eq 'rocks');
		
		my $url_query = $base_url.'/polyweb/vector/gene.html?project='.$pro->{'name'};
		$url_query =~ s/\/\//\//g;
		$out .= qq{<td><center><button onclick="window.open('$url_query', '_blank')" class='btn btn-success btn-sm'>PolyQuery</button></center></td>};
		
		if ($pro->{'type_analyse'} eq 'target') {
			my $url_diag = $base_url.'/polyweb/coverage.html?project='.$pro->{'name'};
			$url_diag =~ s/\/\//\//g;
			$out .= qq{<td><center><button onclick="window.open('$url_diag', '_blank')" class='btn btn-warning btn-sm'>PolyDiag</button></center></td>};
		}
		else {
			if ($pro->{'cache'} eq 'rocks') {
				my $url_viewer = $base_url.'/polyweb/polyviewer.html?project='.$pro->{'name'};
				$url_viewer =~ s/\/\//\//g;
				$out .= qq{<td><center><button onclick="window.open('$url_viewer', '_blank')" class='btn btn-primary btn-sm'>PolyViewer</button></center></td>};
			}
			elsif (@samples && -e "$dir_polyviewer/".$samples[0]->{name}.".variants-genes.polyviewer") {
				my $url_viewer = $base_url.'/polyweb/polyviewer.html?project='.$pro->{'name'};
				$url_viewer =~ s/\/\//\//g;
				$out .= qq{<td><center><button onclick="window.open('$url_viewer', '_blank')" class='btn btn-primary btn-sm'>PolyViewer</button></center></td>};
			}
			else {
				$out .= "<td></td>";
			}
		}
		push(@lHTML, $out);
	}
	
	my $out2 = $cgi->start_div();
	$out2 .= qq{<table data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-v-align="top" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="20" data-page-list="[10, 20]" data-resizable='true' id='table_projects' class='table table-striped' style='font-size:13px;'>};
	$out2 .= "<thead>";
	$out2 .= $cgi->start_Tr({style=>"background-color:#E9DEFF;font-size:10px"});
	$out2 .= qq{<th data-field="name" data-sortable="true" data-filter-control="input" data-filter-control-placeholder="NGS2022_4000">Name</th>};
	$out2 .= qq{<th data-field="description" data-sortable="true" data-filter-control="input" data-filter-control-placeholder="">Description</th>};
	$out2 .= qq{<th data-field="samples" data-sortable="true" data-filter-control="input" data-filter-control-placeholder="Patient Name">Samples</th>};
	$out2 .= qq{<th data-field="type" data-filter-control="select">Type</th>};
	$out2 .= qq{<th data-field="capture" data-sortable="true" data-filter-control="input">Capture Name</th>};
	$out2 .= qq{<th data-field="phenotype" data-filter-control="input">Phenotype(s)</th>};
	$out2 .= qq{<th data-field="genome" data-sortable="true" data-filter-control="input">Genome</th>};
	$out2 .= qq{<th data-field="annotation" data-sortable="true" data-filter-control="input">Annotation</th>};
	$out2 .= qq{<th data-field="polyquery">PolyQuery</th>};
	$out2 .= qq{<th data-field="polyviewer">PolyViewer/Diag</th>};
	$out2 .= $cgi->end_Tr();
	$out2 .= "</thead>";
	$out2 .= "<tbody>";
	
	foreach my $line (@lHTML) { $out2 .= "<tr>$line</tr>"; }
	$out2 .= "</tbody>";
	$out2 .= "</table>";
	$out2 .= "</div>";
	$out2 .= "<br><br>";
	$out2 .= get_imagine_descarte_hml_div(); 
	
	my $hRes;
	$hRes->{html} = $out2;
	
	printJson($hRes);
	exit(0);
	
}


sub printJson {
	my ($hashRes, $test) = @_;
	my $json_encode = encode_json $hashRes;
	print ".\",";
	$json_encode =~ s/{//;
	print $json_encode;
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

sub get_imagine_descarte_hml_div {
	my $html = qq{
		<br>
		<center>
		<table>
			<td style="width=50%;">
				<div style="height:80px;padding-top:3px;">
					<div style="height:80px;">
						<a href="http://www.institutimagine.org/fr/" target="_blank">
							<img src="install/images/imagine_2.png" class="img-fluid" alt="" height="80px" whidth="auto">
							<div class="mask rgba-white-light"></div>
						</a>
					</div>
				</div>
			</td>
			<td style='padding-left:150px;'>			
				<div style="height:80px;padding-top:3px;">
						<div style="height:80px;">
							<a href="https://www.parisdescartes.fr/" target="_blank">
								<img src="install/images/logo_universite_paris_descartes.png" class="img-fluid" alt="" height="80px" whidth="auto">
								<div class="mask rgba-white-light"></div>
							</a>
						</div>
				</div>
			</td>
		</table>
		</center>
		<br><br><br>
	};
	return $html;
}
	