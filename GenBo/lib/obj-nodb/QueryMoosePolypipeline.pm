package QueryMoosePolypipeline;

use strict;
use Vcf;
use Moo;

use Data::Dumper;
use Config::Std;
extends "QueryMoose";



##Requettes pour récupérer les informations des tables de la base Polypipeline

has sql_cmd_list_all_analysis => (
	is	 => 'ro',
	isa	 => 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{select pa.analysis_id as analysis_id, pa.project_id as project_id, pa.project_name as project_name, 
			pa.patients_names as patients_names, psa.status_name as status, pa.user as user, pa.type as type, pa.date as date, pa.patients_ids as patients_ids, pa.pipeline_steps as pipeline_steps
			 from Polypipeline.Analysis pa, Polypipeline.Status_Analysis psa 
			 where psa.status_id= pa.status;};
	return $sql
	},
);

has sql_cmd_list_all_projects => (
	is	 => 'ro',
	isa	 => 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{select distinct project_name from Polypipeline.Analysis ;};
	return $sql
	},
);

has sql_cmd_list_all_versions => (
	is	 => 'ro',
	isa	 => 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift;
#		my $sql = qq{select version_id from Polypipeline.Versions;};
		my $sql = qq{select * from Polypipeline.Versions;};
	return $sql
	},
);

has sql_cmd_list_all_programs => (
	is	 => 'ro',
	isa	 => 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{select * from Polypipeline.Programs;};
	return $sql
	},
);

has sql_cmd_list_all_references => (
	is	 => 'ro',
	isa	 => 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{select * from Polypipeline.References;};
	return $sql
	},
);

has sql_cmd_list_all_steps => (
	is	 => 'ro',
	isa	 => 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{select * from Polypipeline.Step_infos;};
	return $sql
	},
);

##Requêtes dans les tables Analysis_* sur selection d'un analysis_id

has sql_cmd_list_steps_by_analysis => (
	is	 => 'ro',
	isa	 => 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift ;
#		my $sql = qq{select * from Polypipeline.Analysis_Steps where analysis_id=? ;};
		my $sql = qq{select pas.analysis_id as analysis_id, ps.step_name as step_name, pas.step_order as step_order, pas.patient as patient, pss.status_name as status, pas.start as start, pas.end as end, pas.command as command , pss.status_name as status_name
			from Polypipeline.Steps_infos ps, Polypipeline.Analysis_Steps pas,  Polypipeline.Status_Steps pss 
			where  pas.analysis_id=? 
			and ps.step_id=pas.step_id 
			and pss.status_id=pas.status;
		};
		warn $sql ;
	return $sql
	},
);

has sql_cmd_list_versions_by_analysis => (
	is	 => 'ro',
	isa	 => 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift ;
		my $sql = qq{SELECT pav.analysis_id as analysis_id, pv.version_name as version_name, pv.version_path as version_path 
			FROM Polypipeline.Analysis_version pav, Polypipeline.Versions pv 
			where pav.analysis_id=? and pv.version_id=pav.version_id ;};
	return $sql
	},
);

has sql_cmd_list_programs_by_analysis => (
	is	 => 'ro',
	isa	 => 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift ;
		my $sql = qq{SELECT pap.analysis_id as analysis_id, pp.program_name as program_name, pp.program_id as program_id, pp.program_path as program_path, pp.program_md5 as program_md5 
			FROM Polypipeline.Analysis_programs pap, Polypipeline.Programs pp 
			where pap.analysis_id=? and pp.program_id=pap.program_id ;};
	return $sql
	},
);

has sql_cmd_list_references_by_analysis => (
	is	 => 'ro',
	isa	 => 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift ;
		my $sql = qq{SELECT par.analysis_id as analysis_id, pr.ref_name as ref_name, pr.ref_id as ref_id, pr.ref_path as ref_path, pr.ref_md5 as ref_md5 
			FROM Polypipeline.Analysis_references par, Polypipeline.References pr 
			where par.analysis_id=? and pr.ref_id=par.ref_id ;};
	return $sql
	},
);

##Méthodes 
sub getVersionsList{
	my $self= shift ;
	my $dbh = $self->getDbh ;
	my $sth = $dbh->prepare($self->sql_cmd_list_all_versions);
	$sth->execute() || die();
	my $res = $sth->fetchall_hashref("version_id");
	warn Dumper $res ;
	my @list = sort{$a->{version_name} cmp $b->{version_name} } values %$res;
	warn Dumper @list ;
	return \@list;

};


sub getReferencesList{
	my $self= shift ;
	my $dbh = $self->getDbh ;
	my $sth = $dbh->prepare($self->sql_cmd_list_all_references);
	$sth->execute() || die();
	my $res = $sth->fetchall_hashref("ref_id");
	warn Dumper $res ;
	my @list = sort{$a->{ref_name} cmp $b->{ref_name} || $a->{ref_path} cmp $b->{ref_path}} values %$res;
	warn Dumper @list ;
	return \@list;

};

sub getProgramsList{
	my $self= shift ;
	my $dbh = $self->getDbh ;
	my $sth = $dbh->prepare($self->sql_cmd_list_all_programs);
	$sth->execute() || die();
	my $res = $sth->fetchall_hashref("program_name");
	warn Dumper $res ;
	my @list = sort{$a->{program_name} cmp $b->{program_name} } values %$res;
	warn Dumper @list ;
	return \@list;

};

sub getAnalysisList{
	my $self= shift ;
	my $dbh = $self->getDbh ;
	my $sth = $dbh->prepare($self->sql_cmd_list_all_analysis);
	$sth->execute() || die();
	my $res = $sth->fetchall_hashref("analysis_id");
	warn Dumper $res ;
	my @list = sort{$a->{analysis_id} cmp $b->{analysis_id} } values %$res;
	warn Dumper @list ;
	return \@list;

};

sub getProjectsList{
	my $self= shift ;
	my $dbh = $self->getDbh ;
	my $sth = $dbh->prepare($self->sql_cmd_list_all_projects);
	$sth->execute() || die();
	my $res = $sth->fetchall_hashref("project_name");
#	warn Dumper $res ;
	my $hProject ;
	my @list_project_names = sort keys %$res;
	foreach my $project_name (@list_project_names){
		warn $project_name ;
#		my $sql = qq{SELECT * FROM Polypipeline.Analysis where project_name="$project_name";};
my $sql2 = qq{select pa.analysis_id as analysis_id,  pa.project_name as project_name, psa.status_name as status
			 from Polypipeline.Analysis pa, Polypipeline.Status_Analysis psa 
			 where psa.status_id= pa.status and pa.project_name="$project_name";};
		warn $sql2 ;
#		my $sth2 = $dbh->prepare($sql);
my $sth2 = $dbh->prepare($sql2);
		$sth2->execute() || die();
		my $res1 = $sth2->fetchall_hashref("analysis_id");
		warn Dumper $res1 ;
		my @list_analysis_id = sort keys %$res1;
	#		warn join(" ", @list_analysis_id);
#		my $sth3 = $dbh->prepare($sql);
my $sth3 = $dbh->prepare($sql2);
		$sth3->execute() || die();
		my $res2 = $sth3->fetchall_hashref("status");
#		warn Dumper $res2 ;
		my @list_status = sort keys %$res2;
		my $final_status = $list_status[0];
#		warn $final_status ;
#		warn join(" ", @list_status);
		$hProject->{$project_name}->{"status"}=$final_status;
		$hProject->{$project_name}->{"analysis_id"}=join(" ", @list_analysis_id);
		$hProject->{$project_name}->{"project_name"}=$project_name;
		
	}
	warn Dumper $hProject ;
#	warn join(" ", @list_project_ids) ;
	my @list = sort{$a->{project_name} cmp $b->{project_name} } values %{$hProject};
	warn @list ;
	return \@list;

};

sub getAnalysisSteps{
	my ($self, $analysis_id) = @_ ;
	my $dbh = $self->getDbh ;
	my $sth = $dbh->prepare($self->sql_cmd_list_steps_by_analysis);
	$sth->execute($analysis_id) ;
	my @keys=('patient','step_order');
#	my $res = $sth->fetchall_hashref("command");
	my $res = $sth->fetchall_hashref(\@keys);
	warn Dumper $res ;
	my @list ;
	foreach my $patient (sort keys %$res){
		warn $patient ;
		warn Dumper $res->{$patient} ;
		my @list_by_patient = sort{$a->{step_order} <=> $b->{step_order} } values %{$res->{$patient}};
		push(@list, @list_by_patient );
#		warn Dumper @list ;
#		warn \@list ;
	}
	
#	my @list = sort{$a->{step_order} cmp $b->{step_order} || $a->{patient} cmp $b->{patient} } values %$res->{$patient} ;
	warn Dumper @list ; #list est une liste de table de hash
	return \@list;
};

sub getAnalysisVersions{
	my ($self, $analysis_id) = @_ ;
	my $dbh = $self->getDbh ;
	my $sth = $dbh->prepare($self->sql_cmd_list_versions_by_analysis);
	$sth->execute($analysis_id);
	my $res = $sth->fetchall_hashref("version_name");
	warn Dumper $res ;
	my @list = sort{$a->{version_name} cmp $b->{version_name} } values %$res;
	warn Dumper @list ;
	return \@list;
};

sub getAnalysisPrograms{
	my ($self, $analysis_id) = @_ ;
	my $dbh = $self->getDbh ;
	my $sth = $dbh->prepare($self->sql_cmd_list_programs_by_analysis);
	$sth->execute($analysis_id) || die();
	my $res = $sth->fetchall_hashref("program_id");
	warn Dumper $res ;
	my @list = sort{$a->{program_name} cmp $b->{program_name} } values %$res;
	warn Dumper @list ;
	return \@list;
};

sub getAnalysisReferences{
	my ($self, $analysis_id) = @_ ;
	my $dbh = $self->getDbh ;
	my $sth = $dbh->prepare($self->sql_cmd_list_references_by_analysis);
	$sth->execute($analysis_id) || die();
	my $res = $sth->fetchall_hashref("ref_id");
	warn Dumper $res ;
	my @list = sort{$a->{ref_name} cmp $b->{ref_name} } values %$res;
	warn Dumper @list ;
	return \@list;
};

1;

