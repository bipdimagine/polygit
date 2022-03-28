package QueryMoosePolyProject;

use strict;
use Vcf;
use Moose;
use MooseX::Method::Signatures;
use Data::Dumper;
use Config::Std;
extends "QueryMoose";



has sql_cmd_list_all_capture_analyse => (
	is	 => 'ro',
	isa	 => 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
				select type as analyse from Polyproject.capture_systems;
		};
	return $sql
	},
);

has sql_cmd_list_ngs_patients_id_byCapture => (
	is	 => 'ro',
	isa	 => 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
				select distinct genbo.GENBO_ID as patient_id from Polyexome_HG19.GENBO genbo, Polyproject.projects proj, Polyproject.capture_systems capt 
					where genbo.TYPE_GENBO_ID=4 and proj.capture_id=capt.capture_id and capt.type=? and genbo.ORIGIN_ID=proj.project_id and proj.type_project_id=3 and not proj.name REGEXP 'TEST';
		};
	return $sql
	},
);

has sql_cmd_list_ngs_projects_name_byCapture => (
	is	 => 'ro',
	isa	 => 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
				select distinct proj.name from Polyexome_HG19.GENBO genbo, Polyproject.projects proj, Polyproject.capture_systems capt 
					where genbo.TYPE_GENBO_ID=4 and proj.capture_id=capt.capture_id and capt.type=? and genbo.ORIGIN_ID=proj.project_id and proj.type_project_id=3 and not proj.name REGEXP 'TEST';
		};
	return $sql
	},
);

has sql_cmd_list_ngs_projects_name_byCaptureId => (
	is	 => 'ro',
	isa	 => 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
				select distinct proj.name from Polyexome_HG19.GENBO genbo, Polyproject.projects proj, Polyproject.capture_systems capt 
					where genbo.TYPE_GENBO_ID=4 and proj.capture_id=capt.capture_id and capt.capture_id=? and genbo.ORIGIN_ID=proj.project_id and proj.type_project_id=3 and not proj.name REGEXP 'TEST';
		};
	return $sql
	},
);

# TODO: � faire...
#has sql_cmd_list_ngs_projects_id => (
#	is	 => 'ro',
#	isa	 => 'Str',
#	lazy =>1,
#	default	=> sub {
#		my $self = shift;
#		my $sql = qq{
#				SELECT genbo.GENBO_ID as patient_id FROM Polyexome_HG19.GENBO genbo, Polyproject.projects proj 
#					where genbo.TYPE_GENBO_ID=4 and genbo.ORIGIN_ID=proj.project_id and proj.type_project_id =3 and not proj.name REGEXP 'TEST';
#		};
#	return $sql
#	},
#);

has sql_cmd_list_ngs_patients_id => (
	is	 => 'ro',
	isa	 => 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
				SELECT genbo.GENBO_ID as patient_id FROM Polyexome_HG19.GENBO genbo, Polyproject.projects proj 
					where genbo.TYPE_GENBO_ID=4 and genbo.ORIGIN_ID=proj.project_id and proj.type_project_id =3 and not proj.name REGEXP 'TEST';
		};
	return $sql
	},
);

has sql_cmd_list_ngs_patients_id_byYear => (
	is	 => 'ro',
	isa	 => 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
				SELECT genbo.GENBO_ID as patient_id, genbo.DATE, proj.name as project_name FROM Polyexome_HG19.GENBO genbo, Polyproject.projects proj 
					where genbo.TYPE_GENBO_ID=4 and genbo.ORIGIN_ID=proj.project_id and proj.type_project_id =3 and not proj.name REGEXP 'TEST';
		};
	return $sql
	},
);

has sql_cmd_list_projects_name => (
is		=> 'ro',
	isa		=> 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
				select name from Polyproject.projects where type_project_id =3 order by project_id;
		};
	return $sql
	},
);
has sql_cmd_owner_project => (
	is		=> 'ro',
	isa		=> 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
		select BU.email as email, BU.PRENOM_U as firstname FROM Polyproject.user_projects U, Polyproject.projects O, bipd_users.USER BU
		where U.project_id=O.project_id and O.project_id=? and U.USER_ID=BU.USER_ID and BU.equipe_id != 6 ;
		};
	return $sql
		
	},
);

has sql_cmd_project_by_id => (
	is		=> 'ro',
	isa		=> 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
			select name as name from Polyproject.projects where project_id=?;
		};
		return $sql;
	},
);

has sql_cmd_project_by_name => (
	is		=> 'ro',
	isa		=> 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{	  	
			select o.project_id as id,o.name as name,o.creation_date as creation_date,o.description as description, pt.name as projectType , db.name as dbname , r.name as version
			from Polyproject.projects o , Polyproject.project_types pt, Polyproject.databases_projects dp,Polyproject.polydb db, Polyproject.project_release pr , Polyproject.releases r 
			  where o.name= ? 
			  and dp.project_id=o.project_id and db.db_id=dp.db_id
			  and pr.project_id=o.project_id and pr.release_id = r.release_id and pr.default=1 and pt.type_project_id = o.type_project_id;
		};
		return $sql;
	},
);



has sql_cmd_sequencing_machine =>(
	is		=> 'ro',
	isa		=> 'Str',
	lazy =>1,
	default	=> sub {
	my $sql = qq{select m.name as name,m.type as type 
				 from Polyproject.projects_machines pm , Polyproject.sequencing_machines m 
				where pm.project_id=? and pm.machine_id=m.machine_id };
	return $sql
	},
);

has sql_get_patients =>(
	is		=> 'ro',
	isa		=> 'Str',
	lazy =>1,
	default	=> sub {
	my $sql = qq{
		select s.name as name, p.project_id as project_id , sample_id as genbo_id, 1 as run_id, sample_id as patient_id, p.capture_id as capture_id from Polyproject.samples s , Polyproject.projects p
		 where s.project_id=? and p.project_id=s.project_id
	};
	return $sql
	},
);

has sql_list_project_for_user =>(
	is		=> 'ro',
	isa		=> 'Str',
	lazy =>1,
	default	=> sub {
		my $sql = qq{	  	
			select o.project_id as id,o.name as name, pt.name as type , db.name as dbname ,o.description as description, BU.EQUIPE_ID team
			from Polyproject.projects o , Polyproject.databases_projects dp,Polyproject.polydb db,  Polyproject.user_projects  up ,  bipd_users.USER BU, Polyproject.project_types pt
			  where 
			  	up.project_id=o.project_id and ((up.USER_ID=BU.USER_ID AND BU.LOGIN=? and BU.PW=password(?)) or 
			   and dp.project_id=o.project_id and db.db_id=dp.db_id  
			   and o.type_project_id = pt.type_project_id order by pt.name,o.name;
	};
	
	return $sql
	},
);

has sql_capture_infos =>(
	is		=> 'ro',
	isa		=> 'Str',
	lazy =>1,
	default	=> sub {
	my $sql = qq{SELECT c.* FROM Polyproject.capture_systems c where c.capture_id= ?;};

	return $sql
	},
);

has sql_cmd_owner_project_by_name => (
	is		=> 'ro',
	isa		=> 'Str',
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
		select BU.email as email, BU.PRENOM_U as firstname FROM Polyproject.user_projects U, Polyproject.projects O, bipd_users.USER BU
		where U.project_id=O.project_id and O.name=? and U.USER_ID=BU.USER_ID and BU.equipe_id != 6 ;
	};
	return $sql
	},
);

has sql_get_origin_methods =>(
	is		=> 'ro',
	isa		=> 'Str',
	lazy =>1,
	default	=> sub {
	my $sql = qq{
		SELECT distinct m.name as methodname 
    		FROM Polyproject.projects pr, Polyproject.methods m , Polyproject.project_methods pm
        		where pr.project_id=? and pr.project_id=pm.project_id and m.method_id=pm.method_id and m.type=?;
	};
	return $sql
	},
);

has sql_get_methods =>(
	is		=> 'ro',
	isa		=> 'Str',
	lazy =>1,
	default	=> sub {
	my $sql = qq{
		select m.* from Polyproject.samples p ,  Polyproject.project_methods rm ,Polyproject.methods m where p.name= ? and p.project_id= ? and p.project_id=rm.project_id and m.method_id=rm.method_id and m.type=?;
	};
	
	return $sql
	},
);
sub connectDB {
	my( $self,$sql)= @_;
	#confess();
	$self->getDbh()->do($sql);
	#I just did something for queryMoOsePOlyprojectNGS and I really don't know why 
	return;
}
1;