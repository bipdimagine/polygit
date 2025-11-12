package QueryMoose;

use strict;
use Carp;
use Moo;

use Data::Dumper;


has config => (
	is		=> 'ro',
	reader	=> 'getConfig',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $config = $self->all_config;
		return $$config{polyprod};
	},
);

has buffer => (
	is       => 'rw',
	required => 1,
	weak_ref => 1,
);


has all_config => (
	is		=> 'ro',
	required=> 1,
);


sub dbh {
	my ($self) = @_;
	return $self->buffer->dbh;
}
sub getDbh {
	my ($self) = @_;
	return $self->buffer->dbh;
}




##### METHODS #####

sub insert_cache_type_rocks {
	my ($self, $project_id) = @_;
	my $sql = qq { UPDATE `PolyprojectNGS`.`projects` SET `cache`='rocks' WHERE `project_id`=$project_id; };
	my $dbh = $self->getDbh();
	$dbh->do($sql) or confess();
	return 1;	
}

sub delteLastConnectionUser {
	my ($self, $user_name) = @_;
	my $user_id = $self->getUserId($user_name);
	confess unless ($user_id);
	my $dbh = $self->getDbh();
	my $sql = $self->sql_cmd_delete_last_connection_user;
	my $sth = $dbh->prepare($sql);
	$sth->execute($user_id);
	return 1;
}

sub addLastConnectionUser {
	my ($self, $user_name) = @_;
	my $user_id = $self->getUserId($user_name);
	confess unless ($user_id);
	my $dbh = $self->getDbh();
	my $sql = $self->sql_cmd_insert_last_connection_user;
	my $sth = $dbh->prepare($sql);
	$sth->execute($user_id);
	return 1;
}


sub getLastConnectionUserProject {
	my ($self, $user_name, $project_id) = @_;
	my $user_id = $self->getUserId($user_name);
	confess unless ($user_id);
	my $dbh = $self->getDbh();
	my $sql = $self->sql_cmd_get_last_connection_user_project;
	my $sth = $dbh->prepare($sql);
	$sth->execute($user_id, $project_id);
	my $h = $sth->fetchall_hashref('user_id');
	if (exists $h->{$user_id}->{date}) { return ($h->{$user_id}->{date},$h->{$user_id}->{since}) }
	return;
}

sub deleteLastConnectionUserProject {
	my ($self, $user_name, $project_id) = @_;
	my $user_id = $self->getUserId($user_name);
	confess unless ($user_id);
	my $dbh = $self->getDbh();
	my $sql = $self->sql_cmd_delete_last_connection_user_project;
	my $sth = $dbh->prepare($sql);
	$sth->execute($user_id, $project_id);
	return 1;
}

sub addLastConnectionUserProject {
	my ($self, $user_name, $project_id) = @_;
	my $user_id = $self->getUserId($user_name);
	confess unless ($user_id);
	my $dbh = $self->getDbh();
	my $sql = $self->sql_cmd_insert_last_connection_user_project;
	my $sth = $dbh->prepare($sql);
	$sth->execute($user_id, $project_id);
	return 1;
}

sub updateLastConnectionUserProject {
	my ($self, $user_name, $project_id,$type) = @_;
	my $user_id = $self->getUserId($user_name);
	confess unless ($user_id);
	my $dbh = $self->getDbh();
	my $sql = $self->sql_cmd_insert_last_connection_user_project;
	my $sth = $dbh->prepare($sql);
	$sth->execute($user_id, $project_id,$type);
	return 1;
}

sub getUserId {
	my ($self, $user_name) = @_;
	my $dbh = $self->getDbh();
	my $sql = $self->sql_cmd_get_user_id;
	my $sth = $dbh->prepare($sql);
	$sth->execute($user_name);
	my $h = $sth->fetchall_hashref('LOGIN');
	if (exists $h->{$user_name}) { return $h->{$user_name}->{USER_ID}; }
	return;
	
}

sub getLastConnectionUser {
	my ($self, $user_name) = @_;
	my $dbh = $self->getDbh();
	my $sql = $self->sql_cmd_get_last_connection_user;
	my $sth = $dbh->prepare($sql);
	$sth->execute($user_name);
	my $h = $sth->fetchall_hashref('LOGIN');
	if (exists $h->{$user_name}) { return $h->{$user_name}->{date}; }
	return;
}

sub getProjectListReleaseAnnotation {
	my ($self, $proj_name) = @_;
	my @res;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_get_project_release_annotation);
	$sth->execute($proj_name);
	my $toto = $sth->fetchall_hashref("id");
	my @list_ids = keys %$toto;
	my @list;
	foreach my $id (sort @list_ids) {
		push(@list, $toto->{$id}->{name});
	}
	return \@list;
}

sub getProjectListReleaseGene {
	my ($self, $proj_name) = @_;
	my @res;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_get_project_release_gene);
	$sth->execute($proj_name);
	my $toto = $sth->fetchall_hashref("id");
	my @list_ids = keys %$toto;
	my @list;
	foreach my $id (sort @list_ids) {
		push(@list, $toto->{$id}->{name});
	}
	return \@list;
}

sub getTranscriptsTooLongFromBundles {
	my $self = shift;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_getTranscriptsTooLongFromBundles);
	$sth->execute() || die();
	return $sth->fetchall_hashref("ENSEMBL_ID");
}

sub getTranscriptsInfosTranscriptsTable {
	my ($self, $enst_id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_getTranscriptsInfosTranscriptsTable);
	$sth->execute($enst_id) || die();
	return $sth->fetchall_hashref("ENSEMBL_ID");
}

sub getAllBundlesIdsFromTranscriptsId {
	my ($self, $tr_id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_getAllBundlesIdsFromTranscriptsId);
	$sth->execute($tr_id) || die();
	return $sth->fetchall_hashref("bundle_id");
}

sub getPatientProjectInfo {
	my ($self, $proj_id) = @_;
	my @res;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_getPatientProjectInfo);
	$sth->execute($proj_id);
	while (my $id = $sth->fetchrow_hashref ) {
		push(@res, $id);
	}
	return \@res;
}

sub getPatientSomaticProjectInfo {
	my ($self, $proj_id) = @_;
	my @res;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_getPatientSomaticProjectInfo);
	$sth->execute($proj_id);
	while (my $id = $sth->fetchrow_hashref ) {
		push(@res, $id);
	}
	return \@res;
}

sub listCaptureDiagnosticFromTranscriptId {
	my ($self, $trans_id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_list_capture_diag_from_transcrip_id);
	$sth->execute($trans_id) || die();
	my $toto = $sth->fetchall_hashref("name");
	my @list = keys %$toto; 
	return \@list;
}

sub hashTransIdWithCaptureDiag {
	my $self = shift;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_transcrip_id_with_capture_diag);
	$sth->execute() || die();
	return $sth->fetchall_hashref("trans_id");
}

sub listProjects {
	my ($self,$project_name) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_list_projects_name);
	$sth->execute() || die();
	my @row;
	my @list ;
	while (@row = $sth->fetchrow_array) {  # retrieve one row
    		next if $row[1] !~/NGS/;
    		
    		push(@list,$row[1]);
	}
		return \@list;

}

sub listProjectsWithoutDejaVu {
	my $self = shift;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_list_projects_name_without_dejavu);
	$sth->execute() || die();
	my $toto = $sth->fetchall_hashref("name");
	my @list = grep {$_ =~ /NGS/ or $_ =~ /ROCK/} keys %$toto; 
	return \@list;
} 

sub listProjectsForDejaVu {
	my $self = shift;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_list_projects_name_for_dejaVu);
	$sth->execute() || die();
	my $toto = $sth->fetchall_hashref("name");
	my @list = grep {$_ =~ /NGS/ or $_ =~ /ROCK/} keys %$toto; 
	return \@list;
}

sub listProjectsExomeForDejaVu {
	my $self = shift;
	
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_list_projects_name_exome_for_dejaVu);
	$sth->execute() || confess();
	my $toto = $sth->fetchall_hashref("name");
	my @list = grep {$_ =~ /NGS/ or $_ =~ /ROCK/} keys %$toto; 
	return \@list;
}


sub listProjectsIdExomeForDejaVu {
	my $self = shift;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_list_projects_id_exome_for_dejaVu);
	$sth->execute() || confess();
	return [keys %{$sth->fetchall_hashref("id")} ];
}

sub listAllProjects {
	my ($self,$user,$pwd) = @_;
	my $dbh = $self->getDbh();
	my $sth;
	if ($user){
		 $sth = $dbh->prepare($self->sql_allprojects_infos_for_user);
	$sth->execute($user) || die();
	}
	else {
		 $sth = $dbh->prepare($self->sql_allprojects_infos);
		$sth->execute() || die();
	}

	my $toto = $sth->fetchall_arrayref();
	
	my @list = grep {$_->[0] =~ /NGS/ or $_->[0] =~ /ROCK/} @$toto; 

	return \@list;
}

sub hashAllPatientsInfos {
	my ($self, $by) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_hash_ngs_patients_infos);
	$sth->execute() || die();
	unless ($by) { $by = "patient_id"; }
	my $hash = $sth->fetchall_hashref($by);
	return $hash;
}

sub listAllPatientsId {
	my $self = shift;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_list_ngs_patients_id);
	$sth->execute() || die();
	my $toto = $sth->fetchall_hashref("patient_id");
	my @list = keys %$toto; 
	return \@list;
}

sub listAllPatientsIdByYear {
	my ($self) = shift;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_list_ngs_patients_id_byYear);
	$sth->execute() || die();
	return $sth->fetchall_hashref("patient_id");
}

sub listAllCaptureAnalyse {
	my ($self) = shift;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_list_all_capture_analyse);
	$sth->execute() || die();
	my $toto = $sth->fetchall_hashref("analyse");
	my @list = keys %$toto; 
	return \@list;
}

sub listAllPatientsIdByCapture {
	my ($self, $captureName) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_list_ngs_patients_id_byCapture);
	$sth->execute($captureName) || die();
	my $toto = $sth->fetchall_hashref("patient_id");
	my @list = keys %$toto; 
	return \@list;
}

sub listAllProjectsNameByCapture {
	my ($self, $captureName) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_list_ngs_projects_name_byCapture);
	$sth->execute($captureName) || die();
	my $toto = $sth->fetchall_hashref("name");
	my @list = keys %$toto; 
	return \@list;
}

sub listAllProjectsNameByCaptureId {
	my ($self, $captureId) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_list_ngs_projects_name_byCaptureId);
	$sth->execute($captureId) || die();
	my $toto = $sth->fetchall_hashref("name");
	my @list = keys %$toto; 
	return \@list;
}

sub listAllProjectsNameByAnalyse {
	my ($self, $analyse) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_list_ngs_projects_name_byAnalyse);
	$sth->execute($analyse) || die();
	my $toto = $sth->fetchall_hashref("name");
	my @list = keys %$toto; 
	return \@list;
}

sub getRuns {
	my ($self, $project_id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_select_run);
	$sth->execute($project_id);
	my $s = $sth->fetchall_hashref('id');
	return [values %$s];
}
sub getRunInfosById {
	my ($self,$run_id) = @_;
	my $dbh = $self->getDbh();
	my $sql = qq{	select r.run_id as id ,creation_date as date , r.plateform_run_name as plateform_run_name,  description as description, p.name as plateform ,m.name as machine,m.type as constructor,s.name as method,
r.document as document,r.name as run_name from PolyprojectNGS.run r,PolyprojectNGS.run_plateform rp , PolyprojectNGS.plateform p , PolyprojectNGS.run_machine rm, PolyprojectNGS.sequencing_machines m , 
PolyprojectNGS.run_method_seq rs, PolyprojectNGS.method_seq s
where r.run_id=? and r.run_id=rp.run_id and rp.plateform_id=p.plateform_id and rm.run_id=r.run_id and rm.machine_id=m.machine_id and rs.run_id=r.run_id and rs.method_seq_id=s.method_seq_id ;};
	my $sth = $dbh->prepare($sql);
	$sth->execute($run_id);
	my $s = $sth->fetchall_hashref('id');
	return [values %$s];
}
sub getSequencingMethodFromRunId {
		my ($self,$run_id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare(q{
    SELECT name
    FROM PolyprojectNGS.run_method_seq rs
    LEFT JOIN PolyprojectNGS.method_seq ms
        ON rs.method_seq_id = ms.method_seq_id
    WHERE rs.run_id = ?
});
$sth->execute($run_id);
my ($name) = $sth->fetchrow_array;  # une seule valeur
return $name;
	
}
sub getProjectNameFromId{
	my ($self, $project_id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_project_by_id());
	$sth->execute($project_id);
	my $s = $sth->fetchall_arrayref({});
	return $s->[0]->{'name'};
}

sub getListAllReleasesAnntotations {
	my ($self, $project_id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_get_list_releases_annotations());
	$sth->execute();
	return sort {$a <=> $b} keys %{$sth->fetchall_hashref("annotation_release")};
}

sub getListProjectReleasesAnntotations {
	my ($self, $project_id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_get_list_project_releases_annotations());
	$sth->execute($project_id);
	return sort {$a <=> $b} keys %{$sth->fetchall_hashref("annotation_release")};
}

sub getReleaseGenome {
	my ($self, $project_id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare("SELECT name FROM PolyprojectNGS.releases r , PolyprojectNGS.project_release p where p.release_id = r.release_id and p.project_id = ?;");
	$sth->execute($project_id);
	my @list = keys %{$sth->fetchall_hashref("name")};
	return $list[0];
}

sub getReleaseGeneFromBundle {
	my ($self, $capture_id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_release_gene_capture_id);
	$sth->execute($capture_id);
	return  [sort {$a <=> $b} keys %{$sth->fetchall_hashref("release_id")}];
	
}

sub getHashAllGenesReleasesAnntotations {
	my ($self, $project_id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_get_current_genes_releases_annotations());
	$sth->execute();
	my $hash = $sth->fetchall_hashref("rel_gene_id");
	return $hash;
}

sub getCurrentGenesReleasesAnntotations {
	my ($self, $project_id) = @_;
	my $hash = $self->getHashAllGenesReleasesAnntotations();
	my @lIds = sort {$a <=> $b} keys %$hash;
	return $hash->{$lIds[-1]}->{name};
}


sub getCurrentDiagProjectReleasesAnntotations {
	my ($self, $project_id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_get_current_diag_project_releases_annotations());
	$sth->execute();
	my @list = keys %{$sth->fetchall_hashref("name")};
	return $list[0];
}

sub getCurrentGenomeProjectReleasesAnntotations {
	my ($self, $project_id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_get_current_genome_project_releases_annotations());
	$sth->execute();
	my @list = keys %{$sth->fetchall_hashref("name")};
	return $list[0];
}

sub getReleaseAnnotationId {
	my ($self, $name) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_get_release_annotation_id());
	$sth->execute($name);
	my $res = $sth->fetchall_arrayref({});
	return $res->[0]->{rel_annot_id};
}

sub getPatientsFromRunId{
	my ($self, $project_id,$run_id) = @_;
	my $dbh = $self->getDbh();
	my $sql = qq{	select patient_id as id  from PolyprojectNGS.patient where project_id=? and run_id =?;};
	my $sth = $dbh->prepare($sql);
	$sth->execute($project_id,$run_id);
	my $s = $sth->fetchall_hashref('id');
	return [values %$s];
}

sub getAllPatientsFromRunId{
	my ($self, $run_id) = @_;
	my $dbh = $self->getDbh();
	my $sql = qq{
		select p.name as patient, p.control as control, pr.name  as project ,p.father,p.mother,p.status,p.patient_id as id,p.sex as sex,p.family as family, cs.name as capture  
			,p.type as type , pr.project_id as project_id, pr.validation_db as validation_db  from PolyprojectNGS.patient p, PolyprojectNGS.projects pr, PolyprojectNGS.capture_systems cs 
				where p.run_id=? and pr.project_id=p.project_id and p.capture_id=cs.capture_id;
	};
	
	my $sth = $dbh->prepare($sql);
	$sth->execute($run_id);
	my $res = $sth->fetchall_arrayref({});
 	return $res ;
}

sub getDuplicateVariations {
	my ($self,$pid,$ver,$db) = @_;
	my $sql = qq{select o.genbo_id as id,nb as NB from DEJAVU_STATIC s,ORIGIN_GENBO o where origin_id=? and  o.type_genbo_id in (8,15,16) and s.genbo_id=o.genbo_id and nb>1};
		my $sth = $self->dbh->prepare($sql);
	$sth->execute($pid);
	my $d = $sth->fetchall_hashref("id");
	my %dd = %$d;
	return $d;
}
sub countPatients {
	my ($self,$projects) = @_;
		return 0 unless @$projects;
	my $case = join(",",map {qq{\"$_\"} } @$projects);
	my $sql = qq{select count(p.patient_id) as nb from PolyprojectNGS.patient as p , PolyprojectNGS.projects as pr where p.project_id = pr.project_id and pr.name IN  ($case)};
	
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($sql);
	$sth->execute();
	
	return ($sth->fetchrow_array) ;
}

sub getSimilarProjects {
	my ($self,$cid) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_get_similar_projects);
	$sth->execute($cid);
	return [keys %{$sth->fetchall_hashref("name")} ];
}
sub getSimilarProjectsByValidation_db {
	my ($self,$vdb) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_get_similar_projects_by_validation_db);
	$sth->execute($vdb);# || confess();
	
	return [keys %{$sth->fetchall_hashref("name")} ];
}
sub getSimilarProjectsIdByValidation_db {
	my ($self,$vdb) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_get_similar_projects_id_by_validation_db);
	$sth->execute($vdb);# || confess();
	
	return [keys %{$sth->fetchall_hashref("id")} ];
}

sub getSimilarProjectsByAnalyse {
	my ($self,$vdb) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_get_similar_projects_by_analyse);
	$sth->execute($vdb);
	return [keys %{$sth->fetchall_hashref("name")} ];
}

sub getSimilarProjectsIdByAnalyse {
	my ($self,$vdb) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_get_similar_projects_id_by_analyse);
	$sth->execute($vdb);
	return [keys %{$sth->fetchall_hashref("id")} ];
}


sub getSimilarProjectsByPhenotype {
	my ($self,$pheno) = @_;
	confess();
	return [] unless ($pheno);
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare(qq{select pr.name as name   FROM PolyprojectNGS.projects as pr ,PolyprojectNGS.phenotype as p , PolyprojectNGS.phenotype_project as pp where pp.project_id=pr.project_id and pp.phenotype_id= p.phenotype_id and p.name=?;});
	warn qq{select pr.name as name   FROM PolyprojectNGS.projects as pr ,PolyprojectNGS.phenotype as p , PolyprojectNGS.phenotype_project as pp where pp.project_id=pr.project_id and pp.phenotype_id= p.phenotype_id and p.name=?;};
	warn $pheno;
	$sth->execute($pheno) or confess();
	return [keys %{$sth->fetchall_hashref("name")} ];
}

sub getSimilarProjectsIdByPhenotype {
	my ($self,$pheno) = @_;
	confess();
	return [] unless ($pheno);
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare(qq{select pr.project_id as id   FROM PolyprojectNGS.projects as pr ,PolyprojectNGS.phenotype as p , PolyprojectNGS.phenotype_project as pp where pp.project_id=pr.project_id and pp.phenotype_id= p.phenotype_id and p.name=?;});
	$sth->execute($pheno) or confess();
	return [keys %{$sth->fetchall_hashref("id")} ];
}
sub getSimilarProjectsIdByPhenotypeId {
	my ($self,$pheno) = @_;
	return [] unless ($pheno);
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare(qq{select p.project_id as id   FROM PolyPhenotypeDB.phenotype_project as p  where  p.phenotype_id=?;});
	$sth->execute($pheno) or confess();
	return [keys %{$sth->fetchall_hashref("id")} ];
}

sub getPhenotypesForProject {
	my ($self,$id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare(qq{select p.name as name   FROM PolyprojectNGS.phenotype as p , PolyprojectNGS.phenotype_project as pp where pp.project_id=? and pp.phenotype_id= p.phenotype_id;});
	$sth->execute($id);
	return [keys %{$sth->fetchall_hashref("name")} ];
}


sub getSequencingMachines {
	my ($self, $pid) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_cmd_sequencing_machine());
	$sth->execute($pid);
	my $s = $sth->fetchall_hashref("name");
	my @lRes = keys(%$s);
	@lRes = sort(@lRes);
	return \@lRes;
}

#query getPublicDatabaseVersion
sub getPublicDatabaseVersion {
	my ($self, $pid) = @_;
	my $dbh = $self->getDbh();
	return $dbh->selectrow_array(qq{SELECT max(rg.version_id) as gencode FROM PolyprojectNGS.project_release_public_database pr, PolyprojectNGS.release_public_database as rg where rg.version_id = pr.version_id and project_id = $pid;});
	#my $s = $sth->fetchall_hashref("$pid");
}

sub getMaxPublicDatabaseVersion {
	my ($self) = @_;
	my $dbh = $self->getDbh();
	return $dbh->selectrow_array(qq{SELECT version_id as gencode FROM PolyprojectNGS.release_public_database as rpd where rpd.default_rocks=1;});
	#my $s = $sth->fetchall_hashref("$pid");
}

sub insertPublicDatabaseVersion {
	my ($self,$pid,$dbversion) = @_;
	my $dbversion_id = $self->existsPublicDatabaseVersionId($dbversion);
	confess() unless $dbversion_id;
	my $dbh = $self->getDbh();
	$dbh->do(qq{insert into  PolyprojectNGS.project_release_public_database (project_id, version_id) VALUES ($pid, $dbversion_id)}) or confess();
	return 1;
}



sub existsPublicDatabaseVersionId {
	my ($self,$version_id) = @_;
	my $dbh = $self->getDbh();
	return $dbh->selectrow_array(qq{SELECT version_id  as version_id from  PolyprojectNGS.release_public_database as rg where version_id = $version_id;});
	#my $s = $sth->fetchall_hashref("$pid");
}


#cahe history
sub ListCacheHistoryVersion {
	my ($self,$pid) = @_;
	my $dbh = $self->getDbh();
	my $sql = qq{SELECT version_annotation  as version_annotation from  PolyprojectNGS.project_release_cache_history  where project_id = $pid order by update_date;};
	#my $s = $sth->fetchall_hashref("$pid");
	return $dbh->selectall_arrayref($sql);
}

#sub selectCacheVersion {
#	my ($self,$pid,$version) = @_;
#	die() unless $version;
#	my $dbh = $self->getDbh();
#	$dbh->do(qq{INSERT INTO `PolyprojectNGS`.`project_release_cache_history` (`project_id`, `version_annotation`,`creation_date`,`update_date`) VALUES ($pid, $version,NOW(),NOW() ) 
#				on DUPLICATE KEY update update_date = NOW();}) or confess();
#	return 1;
#}

sub insertHistoryCacheVersion {
	my ($self,$pid,$version) = @_;
	die() unless $version;
	my $dbh = $self->getDbh();
	my $sql = qq{INSERT INTO `PolyprojectNGS`.`project_release_cache_history` (`project_id`, `version_annotation`,`creation_date`,`update_date`) VALUES ($pid, "$version", NOW(), NOW() ) on DUPLICATE KEY update update_date = NOW();};
	$dbh->do($sql) or confess();
	return 1;
}

#query gencode

sub getGencodeVersion {
	my ($self, $pid) = @_;
	my $dbh = $self->getDbh();
	return $dbh->selectrow_array(qq{SELECT max(name) as gencode FROM PolyprojectNGS.project_release_gene pr, PolyprojectNGS.release_gene as rg where rg.rel_gene_id = pr.rel_gene_id and project_id = $pid;});
	#my $s = $sth->fetchall_hashref("$pid");
}

sub getMaxGencodeVersion {
	my ($self) = @_;
	my $dbh = $self->getDbh();

	return $dbh->selectrow_array(qq{SELECT name as gencode FROM PolyprojectNGS.release_gene as rg where rg.default=1;});

	#my $s = $sth->fetchall_hashref("$pid");
}

sub getGencodeId {
	my ($self,$gencode) = @_;
	my $dbh = $self->getDbh();
	return $dbh->selectrow_array(qq{SELECT rel_gene_id  as gencode_id from  PolyprojectNGS.release_gene as rg where name = $gencode;});
	#my $s = $sth->fetchall_hashref("$pid");
}

sub insertGencodeVersion {
	my ($self,$pid,$gencode) = @_;
	my $gencode_id = $self->getGencodeId($gencode);
	confess() unless $gencode_id;
	#die($gencode_id);
	my $dbh = $self->getDbh();
	$dbh->do(qq{insert into  PolyprojectNGS.project_release_gene (project_id, rel_gene_id) VALUES ($pid, $gencode_id)}) or confess();
	return 1;
}



sub prepare {
	my ($self,$query) = @_;
	return $self->{prepare}->{$query} if exists $self->{prepare}->{$query};
	my $dbh = $self->getDbh();
	$self->{prepare}->{$query} = $dbh->prepare($query);
	return $self->{prepare}->{$query};
}

sub getUserGroupsForProject {
		my ($self, $pid) = @_;
	  	my $query =  qq{select ug.name FROM PolyprojectNGS.ugroup_projects G, bipd_users.UGROUP ug  
		where G.project_id=? and ug.ugroup_id=G.UGROUP_ID ;};
		my $sth = $self->prepare($query);
		$sth->execute($pid);
		my $res = $sth->fetchall_hashref("name");
		delete $res->{STAFF};
		return [keys %$res];
}

sub fast_getProjectByName {
	my ($self, $name, $verbose) = @_;
	confess() unless $name;
	my $dbh = $self->getDbh();
	my $config = $self->getConfig();
	#$self->project_name($name);
	my $query = qq{select o.project_id as id, o.name as name,o.creation_date as creation_date,o.description as description,o.validation_db as validation_db,  pt.type_project_id as projectTypeId, pt.name as projectType,r.name as version, o.somatic as is_somatic from PolyprojectNGS.projects o, PolyprojectNGS.project_release pr , PolyprojectNGS.releases r, PolyprojectNGS.project_types pt
		 where o.name IN ($name) 
		   and pr.project_id=o.project_id and pr.release_id = r.release_id and pr.default=1 and o.type_project_id=pt.type_project_id;};
		 
		   
	my $sth = $dbh->prepare($query);
	$sth->execute();
	my $res = $sth->fetchall_hashref("name");
	
	return $res;

}
 ###
sub fast_getPatients {
	my ($self, $project_id) = @_;
	my $dbh = $self->getDbh();
	my $query = qq{select * from PolyprojectNGS.patient g where project_id IN ($project_id) };
	
	my $sth = $dbh->prepare($query);
	$sth->execute();
	my $s = $sth->fetchall_hashref('patient_id');
	my $res;
	foreach my $p (values %$s){
		my $pid = $p->{project_id};
		
		push(@{$res->{$pid}},$p);
	}
	#die($project_id);
	return $res;
}

sub getProjectByName {
	my ($self, $name, $verbose) = @_;
	confess("-".$name."-") unless $name;
	my $dbh = $self->getDbh();
	my $config = $self->getConfig();
	#$self->project_name($name);
	my $sth = $self->prepare($self->sql_cmd_project_by_name);
	$sth->execute("$name");
	my $res = $sth->fetchall_hashref("dbname");

	if (scalar(keys %$res) ==0 ){ confess("no project found $name"); }
	if (scalar(keys %$res) >1 && ! exists $config->{name}){
		while (my ($k, $v)=each(%$res)){
			print "  K: $k  -  V: $v\n";
			while (my ($k2, $v2)=each(%$v)){
				print "     K2: $k  -  V2: $v\n";
			}
		}
		confess("2 different databases for the same project\n");
	}
	my ($dbname) = keys %$res;

	my $s = $res->{$dbname};
	$config->{root}= $dbname;
	if (lc($dbname) eq "polyexome" || lc($dbname) eq "polyrock"){ $config->{name} = $dbname."_".$res->{$dbname}->{version}; }
	else { $config->{name} = $dbname; }#."_".$res->{$dbname}->{version}; }
	$dbname = $config->{name};
	my $sql2 = qq{use $dbname;};
	my $user = getlogin();
	$self->connectDB($sql2);

	return $s;
}

sub connectDB {
	my ( $self,$sql) = @_;
	#I just did something for queryMoOsePOlyprojectNGS and I really don't know why 
	return;
}

sub update_software_version {
	my ($self, $version_id,$patient_id, $cmd,$verbose) = @_;
	confess("-".$version_id."-") unless $version_id;
	confess("-".$patient_id."-") unless $patient_id;
	my $dbh = $self->getDbh();
	my $config = $self->getConfig();
	my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);
	$dbh->do(qq{insert into  PolyprojectNGS.version_patients  (version_id, patient_id,modification_date,cmd) VALUES ($version_id, $patient_id,NOW(),"$username:$cmd" )}) or confess();
	return 1;
}


sub getLatestSoftwareVersion {
	my ($self, $name, $verbose) = @_;
	confess("-".$name."-") unless $name;
	my $dbh = $self->getDbh();
	my $config = $self->getConfig();
	#$self->project_name($name);
	my $sth = $self->prepare("SELECT version_id,version  FROM PolyprojectNGS.version_software where name=? ORDER BY version_id DESC LIMIT 0, 1");
	$sth->execute("$name");
	my $res = $sth->fetchall_arrayref({});
	return ($res->[0]->{version_id},$res->[0]->{version}) if $res;
	return (undef,undef);

}


sub getLatestSoftwareVersionByPatient {
	my ($self, $name,$patient_id) = @_;
	confess("-".$name."-") unless $name;
	my $dbh = $self->getDbh();
	my $config = $self->getConfig();
	#$self->project_name($name);
	my $sth = $self->prepare("SELECT vs.version_id as version_id ,vs.version as version  FROM PolyprojectNGS.version_software as vs , PolyprojectNGS.version_patients as vp where name=? and vs.version_id=vp.version_id and vp.patient_id= ? ORDER BY vs.version_id DESC LIMIT 0, 1");
	$sth->execute("$name",$patient_id);
	my $res = $sth->fetchall_arrayref({});
	return ($res->[0]->{version_id},$res->[0]->{version}) if $res;
	return (undef,undef);

}

sub getCaptureInfos {
	my ($self, $capture_id) = @_;
	my $dbh = $self->getDbh();
	my $config = $self->getConfig();
	my $sth = $dbh->prepare($self->sql_capture_infos);
	$sth->execute($capture_id);
	my $res = $sth->fetchall_arrayref({});
	my %hashCaptureFilesName;
	for (@$res) { $hashCaptureFilesName{$_->{'filename'}} = 1; }
	if (scalar(keys(%hashCaptureFilesName)) > 1) {
		warn "Multiple capture filed declared ! Exit...\n";
		foreach my $name (keys(%hashCaptureFilesName)) { print "  -> $name\n"; }
		die();
	}
	return $$res[0];
}

sub getUmiInfos {
	my ($self, $capture_id) = @_;
	my $dbh = $self->getDbh();
	my $config = $self->getConfig();
	
	my $sth = $dbh->prepare($self->sql_umi);
	$sth->execute($capture_id);
	my $res = $sth->fetchall_arrayref({});
	my %hashCaptureFilesName;
	#for (@$res) { $hashCaptureFilesName{$_->{'filename'}} = 1; }
	if (scalar(keys(%hashCaptureFilesName)) > 1) {
		warn "Multiple capture filed declared ! Exit...\n";
		foreach my $name (keys(%hashCaptureFilesName)) { print "  -> $name\n"; }
		die();
	}
	return $$res[0];
}

sub getCaptureId {
	my ($self, $capture_name) = @_;
	my $dbh = $self->getDbh();
	my $config = $self->getConfig();
	
	my $sth = $dbh->prepare($self->sql_capture_infos_by_name);
	$sth->execute($capture_name);
	my $res = $sth->fetchall_arrayref({});
	return $res->[0]->{capture_id};
}	
	

sub getPatients {
	my ($self, $project_id) = @_;
	my $dbh = $self->getDbh();

	my $sth = $self->prepare($self->sql_get_patients);
	$sth->execute($project_id);
	my $s = $sth->fetchall_hashref('name');
	return [values %$s];
}

sub getGroups {
	my ($self, $project_id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_get_groups);
	$sth->execute($project_id);
	my $s = $sth->fetchall_hashref('pname');
	return [values %$s];
}
sub getPatient {
confess();
}

sub getGenboIdPatient {
confess();
}


sub getOriginMethods {
	my ($self, $project_id, $type) = @_;
	my $dbh = $self->getDbh();
	my $config = $self->getConfig();
	$type = uc($type);
	
	my $sth = $dbh->prepare($self->sql_get_origin_methods);
	$sth->execute($project_id,$type);
	my $s = $sth->fetchall_hashref("methodname");
	my @lRes = sort(keys(%$s));
	
	return \@lRes;
}

sub getMethods{
		my ($self,%arg) = @_;
	my $patient_name = $arg{patient_name};
	my $project_id = $arg{project_id};
	my $type = $arg{type};
	my $dbh = $self->getDbh();
	my $config = $self->getConfig();
	my $sth = $dbh->prepare($self->sql_get_methods);
	$sth->execute($patient_name,$project_id,$type);
	my $s = $sth->fetchall_hashref('name');
	my @lRes = sort(keys(%$s));
	return \@lRes;
}

sub getCallingMethods{
	my ($self,%arg) = @_;
	my $patient_name = $arg{patient_name};
	my $project_id = $arg{project_id};
	my $type = "SNP";
	$type = $arg{type} if exists $arg{type};
	return $self->getMethods(type=>$type, patient_name=>$patient_name,project_id=>$project_id);

}
sub getAlignmentMethods {
	my ($self,%arg) = @_;
	return $self->getMethods(type=>"ALIGN", patient_name=>$arg{patient_name},project_id=>$arg{project_id});
}




sub getOwnerProject {
	my ($self,$pid) = @_;
	my $dbh = $self->getDbh;
	my $sth = $dbh->prepare( $self->sql_cmd_owner_project );
	$sth->execute($pid);
	my $res = $sth->fetchall_arrayref({});
	
	$sth = $dbh->prepare( $self->sql_cmd_owner_group_project );
	$sth->execute($pid);
	my $res2 = $sth->fetchall_arrayref({});
	push(@$res,@$res2);
	return $res;
}


sub getOwnerProject_byName {
	my ($self,$projName) = @_;
	my $dbh = $self->getDbh;
	my $sth = $dbh->prepare( $self->sql_cmd_owner_project_by_name );
	$sth->execute($projName);
	my $res = $sth->fetchall_arrayref({});
	return $res;
}

sub getBundleTranscripts {
	return [];
}



sub getCaptureTranscriptsbyName {
	my ($self,$capture_id) = @_;
	my $dbh = $self->getDbh;
	my $sth = $dbh->prepare($self->sql_cmd_capture_transcripts_by_name);
	
	$sth->execute($capture_id);
	my $res = $sth->fetchall_arrayref({});
	my $res2 = {};
	$res2->{bundle} = {};
	$res2->{transcripts} = {};

	foreach my $h (@$res) {
		$res2->{transcripts}->{$h->{ENSEMBL_ID}} ++;
	#	push(@{$res2->{bundle}->{$h->{name}}},$h); 
	#	push (@{$res2->{transcripts}->{$h->{ENSEMBL_ID}}},$h->{name});
	}
	#	$res2->{transcripts_name} = [keys %{$res2->{transcripts}}];
	return $res2;
}

sub getTranscriptsFromBundle {
	my ($self,$bundle_name) = @_;
	my $dbh = $self->getDbh;
	my $query = qq{SELECT tr.ENSEMBL_ID as id , tr.GENE as gene ,tr.* FROM 
	PolyprojectNGS.transcripts as tr , PolyprojectNGS.bundle_transcripts as bt , PolyprojectNGS.bundle as b where
 	bt.transcript_id= tr.ID and bt.bundle_id = b.bundle_id and b.name= ?;
	};
	
	my $sth = $dbh->prepare($query );
	$sth->execute($bundle_name);
	my $res = $sth->fetchall_hashref("id");
	return $res;
}
sub getBundleFromCapture {
	my ($self,$capture_id) = @_;
	my $dbh = $self->getDbh;
	my $query = qq{
		SELECT b.bundle_id as id, b.name  FROM 
	PolyprojectNGS.capture_systems as c, PolyprojectNGS.capture_bundle as bt , PolyprojectNGS.bundle as b where
	bt.bundle_id= b.bundle_id and bt.capture_id = c.capture_id and c.capture_id= ?
	};
	my $sth = $dbh->prepare($query );
	$sth->execute($capture_id);
	my $res = $sth->fetchall_hashref("id");
	return $res;
}

sub getCaptureTranscripts {
	my ($self,$capture_id) = @_;
	my $dbh = $self->getDbh;
	my $sth = $dbh->prepare($self->sql_cmd_capture_transcripts );

	$sth->execute($capture_id);
	my $res = $sth->fetchall_arrayref({});
	my $res2 = {};
	$res2->{bundle} = {};
	$res2->{transcripts} = {};
	my $tr;
	
	foreach my $h (@$res){
		push(@{$res2->{bundle}->{$h->{name}}},$h); 
		push (@{$res2->{transcripts}->{$h->{ENSEMBL_ID}}},$h->{name});
		
	
	}
		$res2->{transcripts_name} = [keys %{$res2->{transcripts}}];
	return $res2;
}


sub getCaptureGenesTranscripts {
	my ($self,$capture_id) = @_;
	my $dbh = $self->getDbh;
	my $sth = $dbh->prepare($self->sql_cmd_capture_transcripts );

	$sth->execute($capture_id);
	my $res = $sth->fetchall_arrayref({});
	my $res2 = {};
	$res2->{bundle} = {};
	$res2->{transcripts} = {};
	my $tr;
	
	foreach my $h (@$res){
		push(@{$res2->{bundle}->{$h->{name}}},$h); 
		push (@{$res2->{transcripts}->{$h->{GENE}}->{$h->{ENSEMBL_ID}}},$h->{name});
		
	
	}
	#$res2->{transcripts_name} = [keys %{$res2->{transcripts}}];
	$res2->{transcripts_name} = $res2->{transcripts};
	return $res2;
}

sub getAllProjects {
	my ($self) = @_;
	
	my $dbh = $self->getDbh;
	my $type_db = $self->getConfig()->{type_db};
	my $sth = $dbh->prepare($self->sql_list_all_projects());
	$sth->execute();
	my $res = $sth->fetchall_hashref("id");
#	my @lProj = sort keys %{$res};
#	return \@lProj;
	my @toto = sort{$a->{name} cmp $b->{name}} values %$res;
	return \@toto; 
}

sub getProjectListForUser {
	my ($self,$login,$pwd)=@_;
	my $dbh = $self->getDbh;
	my $type_db = $self->getConfig()->{type_db};

	
	my $sth = $dbh->prepare($self->sql_list_project_for_user);
	
	$sth->execute($login,$pwd,$login,$pwd);
	my $res = $sth->fetchall_hashref("id");
	my $res_group = $self->getProjectHashForGroup($login,$pwd);
	if ($res_group) {
		foreach my $project_id (keys %{$res_group}) {
			$res->{$project_id} = $res_group->{$project_id} unless (exists $res->{$project_id});
		}
	}
	
	#my @toto = sort{$a->{type} cmp $b->{type} || $a->{name} cmp $b->{name}} values %$res;
	my @toto =   values %$res;
	
	return \@toto; 
}


sub getProjectHashForGroup {
	my ($self,$login,$pwd)=@_;
	my $dbh = $self->getDbh;
	my $sth = $dbh->prepare($self->sql_list_project_for_group);
	$sth->execute($login,$pwd);
	my $res = $sth->fetchall_hashref("id");
	return $res;
}

sub writeLatestInfos {
	my ($self,$project_id,$login,$pwd)=@_;
	my $dbh = $self->getDbh;
	my $type_db = $self->getConfig()->{type_db};
	my $query = qq{	  	
			select h.* from PolyprojectNGS.project_release_cache_history h   where  h.project_id=$project_id order by update_date DESC limit 1;
	};
	
	
	my $sth = $dbh->prepare($query);
	$sth->execute();
	my $res = $sth->fetchall_arrayref();
	my $id = $res->[0]->[0];
	my $query2 = qq{	  	
			update PolyprojectNGS.project_release_cache_history h  set latest_polyweb_activity = NOW() ,latest_user='$login' where id=$id;
	};
	$dbh->do($query2);
	return 1;
}

sub getProjectHashForUser {
	my ($self,$login,$pwd)=@_;
	my $dbh = $self->getDbh;
	my $type_db = $self->getConfig()->{type_db};
	my $sth = $dbh->prepare($self->sql_list_project_for_user);
	$sth->execute($login,$pwd,$login,$pwd);
	return $sth->fetchall_hashref("name");
}

sub getProjectsList {
	my ($self,$login,$pwd)=@_;
	my $dbh = $self->getDbh;
	my $type_db = $self->getConfig()->{type_db};
	my $sth = $dbh->prepare($self->sql_list_projects);
	$sth->execute();
	my $res = $sth->fetchall_hashref("name");
	my @toto = sort{ $a->{name} cmp $b->{name}} values %$res;
	return \@toto; 
}

sub getUserType {
	my ($self,$login,$pwd)=@_;
	my $dbh = $self->getDbh;
	my $type_db = $self->getConfig()->{type_db};
		my $sql = qq{	  	
			select BU.login login,BU.equipe_id
			from   bipd_users.USER BU 
			  where 
			  BU.LOGIN=? and BU.PW=password(?);
	};
	my $sth = $dbh->prepare($self->sql_list_project_for_user);
	$sth->execute($login,$pwd,$login,$pwd);
	my $res = $sth->fetchall_hashref("id");
	my @toto = sort{$a->{type} cmp $b->{type} || $a->{name} cmp $b->{name}} values %$res;
	return \@toto; 
}

sub isUserMagic {
	my ($self,$login,$pwd)=@_;
	my $dbh = $self->getDbh;
	my $type_db = $self->getConfig()->{type_db};
		my $sql = qq{	  	
			select BU.login login,BU.equipe_id
			from   bipd_users.USER BU 
			  where 
			  BU.LOGIN=? ;
	};
	my $sth = $dbh->prepare($sql);
	
	$sth->execute($login);
	my $res = $sth->fetchall_arrayref();
	return unless $res;
	
	return 1 if $res->[0]->[1] && $res->[0]->[1] eq '6'; 
}
sub getProjectDestination {
	my ($self,$pid)=@_;
	my $dbh = $self->getDbh;
	my $type_db = $self->getConfig()->{type_db};
	my $cmd = qq{SELECT pr.name as name FROM PolyprojectNGS.patient p, PolyprojectNGS.projects pr where patient_id=? and p.project_id_dest=pr.project_id;};
	my $sth = $dbh->prepare($cmd);
	$sth->execute($pid);
	my $res = $sth->fetchall_hashref("name");
	return unless $res;
	my @res1 = keys %$res;
	
	die() if scalar(@res1) > 1;
	return $res1[0];

}

sub isBipd{
	my ($self,$id)=@_;
	my $dbh =  $self->getDbh;
	my $sql = qq{select nb_project, count(distinct g) as nb_patient, sum(he) as he , sum(ho) as ho from  (select  d.nb as nb_project, d.genbo_id as did, g.genbo_id  as g , if (ra.he>=1,1,0) as he, if (ra.ho>=1,1,0) as ho
from DEJAVU_STATIC d, RELATION r,GENBO g, RELATION_ANNEX ra
 where d.NB >1 and d.GENBO_ID = $id and r.genbo2_id= d.genbo_id and g.genbo_id=r.genbo_id
  and g.type_genbo_id = 5
  and ra.relation_id=r.relation_id  group by ra.relation_id,d.genbo_id ) as tbl ;
	};
	#my $sql = qq{select  d.nb as nb_project ,count(g.genbo_id) as nb_patient, sum(ra.he) as he , sum(ra.ho) as ho  from DEJAVU_STATIC d, RELATION r,GENBO g, RELATION_ANNEX ra where d.NB >1 and d.GENBO_ID = $id and r.genbo2_id= d.genbo_id and g.genbo_id=r.genbo_id and g.type_genbo_id = 5  and ra.relation_id=r.relation_id and ra.he>=0 and ra.ho>=0 group by d.genbo_id;};
	my $sth = $dbh->prepare($sql);
	$sth->execute();
	my $col1;
	my $col2;
	my $ho;
	my $he;
	$sth->bind_columns(\$col1,\$col2,\$he,\$ho);
	$sth->fetch();
	#my @res = $sth->fetchrow_array();
 	return ($col1,$col2,$he,$ho);
}

sub isBIPDForOneProject{
	my ($self,$pid,$id)=@_;
	my $dbh =  $self->getDbh;
	confess() unless $pid;
	my $sql = qq{select nb_project, count(distinct g) as nb_patient, sum(he) as he , sum(ho) as ho from  (select  d.nb as nb_project, d.genbo_id as did, g.genbo_id  as g , if (ra.he>=1,1,0) as he, if (ra.ho>=1,1,0) as ho
from DEJAVU_STATIC d, RELATION r,GENBO g, RELATION_ANNEX ra
 where d.NB >1 and d.GENBO_ID = $id and r.genbo2_id= d.genbo_id and g.genbo_id=r.genbo_id
  and g.type_genbo_id = 5
  and g.origin_id=$pid
  and ra.relation_id=r.relation_id  group by ra.relation_id,d.genbo_id ) as tbl ;
	};
	my $sth = $dbh->prepare($sql);
	$sth->execute();
	my $col1;
	my $col2;
	my $ho;
	my $he;
	$sth->bind_columns(\$col1,\$col2,\$he,\$ho);
	$sth->fetch();
	#my @res = $sth->fetchrow_array();
 	return ($col1,$col2,$he,$ho);
}

sub getBundleGenesAllProjectsUsers {
	my ($self, $user) = @_;
	my $dbh =  $self->getDbh;
	confess() unless $user;
	my $sth = $dbh->prepare($self->sql_cmd_getBundleGenesAllProjectsUsers);
	$sth->execute($user);
	my $res = $sth->fetchall_hashref("id");
	return $res; 
}

sub getGenesNamesInBundle {
	my ($self, $bundle_id) = @_;
	my $dbh =  $self->getDbh;
	confess() unless $bundle_id;
	my $sth = $dbh->prepare($self->sql_cmd_getGenesNamesInBundle);
	$sth->execute($bundle_id);
	my $res = $sth->fetchall_hashref("gene_name");
	return $res; 
}

sub getAllGenesNamesInAllBundle {
	my $self = shift;
	my $dbh =  $self->getDbh;
	my $res;
	my $emps = $dbh->selectall_arrayref( $self->sql_cmd_getGenesNamesAllBundle, { Slice => {} } );
	foreach my $h (@$emps) {
		my $capture = $h->{capture};
		my $bundle = $h->{name};
		my $id = $h->{id};
		my $description = $h->{description};
		my $gene_name = $h->{gene_name};
		my $transcript_id = $h->{transcript_id};
		#$res->{$capture}->{$bundle}->{id} = $id;
		#$res->{$capture}->{$bundle}->{description} = $description;
		push(@{$res->{$capture}->{$bundle}->{genes}->{$gene_name}}, $transcript_id);
		push(@{$res->{$capture}->{$bundle}->{transcripts}}, $transcript_id);
	} 
	return $res;
}

sub getOmimGenesTranscriptsNames {
	my $self = shift;
	my $dbh =  $self->getDbh;
	my $res;
	my $emps = $dbh->selectall_arrayref( $self->sql_cmd_getOmimGenesTranscriptsNames, { Slice => {} } );
	foreach my $h (@$emps) {
		my $capture = $h->{capture};
		my $bundle = $h->{name};
		my $id = $h->{id};
		my $description = $h->{description};
		my $gene_name = $h->{gene_name};
		my $transcript_id = $h->{transcript_id};
		push(@{$res->{$capture}->{$bundle}->{genes}->{$gene_name}}, $transcript_id);
		push(@{$res->{$capture}->{$bundle}->{transcripts}}, $transcript_id);
	} 
	return $res;
}

sub get_hash_projects_ids_genes_databases_version {
	my $self = shift;
	my $dbh =  $self->getDbh;
	my $sth = $dbh->prepare($self->sql_cmd_get_proj_ids_genes_datatbase_version);
	$sth->execute();
	my $res = $sth->fetchall_hashref("project_id");
}

1;