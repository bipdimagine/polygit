package QueryMoosePolyProjectNgs;

use strict;
use Vcf;
use Moo;

use Data::Dumper;
use Config::Std;
extends "QueryMoose";




has sql_cmd_delete_last_connection_user => (
	is	 => 'ro',
	lazy => 1,
	default	=> sub {
		my $sql = qq{
			DELETE FROM `PolyprojectNGS`.`user_last_connection` WHERE `user_id`=?;
		};
		return $sql;
	},
);

has sql_cmd_insert_last_connection_user => (
	is	 => 'ro',
	lazy => 1,
	default	=> sub {
		my $sql = qq{
			INSERT INTO `PolyprojectNGS`.`user_last_connection` (`user_id`) VALUES (?) ON DUPLICATE KEY UPDATE date=now();
		};
		return $sql;
	},
);

has sql_cmd_get_last_connection_user_project => (
	is	 => 'ro',
	lazy => 1,
	default	=> sub {
		my $sql = qq{
			SELECT *,DATEDIFF(current_date,u.date) as since FROM PolyprojectNGS.user_projects_last_connection u where u.user_id=? and u.project_id=?;
		};
		return $sql;
	},
);

has sql_cmd_insert_last_connection_user_project => (
	is	 => 'ro',
	lazy => 1,
	default	=> sub {
		my $sql = qq{
			INSERT INTO `PolyprojectNGS`.`user_projects_last_connection` (`user_id`, `project_id`,`viewer`, `date`) VALUES (?, ?,?, NOW()) on DUPLICATE KEY update date=NOW();
		};
		return $sql;
	},
);

has sql_cmd_delete_last_connection_user_project => (
	is	 => 'ro',
	lazy => 1,
	default	=> sub {
		my $sql = qq{
			DELETE FROM `PolyprojectNGS`.`user_projects_last_connection` WHERE `user_id`=? and `project_id`=?;
		};
		return $sql;
	},
);


has sql_cmd_get_ill_patients_from_project => (
	is	 => 'ro',
	lazy => 1,
	default	=> sub {
		my $sql = qq{
			SELECT pat.name as pat_name FROM PolyprojectNGS.patient as pat, PolyprojectNGS.projects as proj
				where pat.project_id = proj.project_id and proj.name=? and pat.status=2;
		};
		return $sql;
	},
);

has sql_cmd_get_user_id => (
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{
			SELECT USER_ID, LOGIN FROM bipd_users.USER where LOGIN=?;
		};
		return $sql;
	},
);

has sql_cmd_get_last_connection_user => (
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{
			SELECT b.LOGIN, p.user_id, p.date FROM bipd_users.USER as b, PolyprojectNGS.user_last_connection as p
				where b.LOGIN=? and b.USER_ID=p.user_id;
		};
		return $sql;
	},
);

has sql_cmd_get_project_release_annotation => (
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{
			SELECT ra.rel_annot_id as id, ra.name as name FROM PolyprojectNGS.projects as p, PolyprojectNGS.project_release_annotation as pra, PolyprojectNGS.release_annotation as ra
				where p.name=? and p.project_id=pra.project_id and pra.rel_annot_id=ra.rel_annot_id;
		};
		return $sql;
	},
);

has sql_cmd_get_project_release_gene => (
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{
			SELECT rg.rel_gene_id as id, rg.name as name FROM PolyprojectNGS.projects as p, PolyprojectNGS.project_release_gene as prg, PolyprojectNGS.release_gene as rg
    			where p.name=? and p.project_id=prg.project_id and prg.rel_gene_id=rg.rel_gene_id;
		};
		return $sql;
	},
);

has sql_cmd_get_release_annotation_id => (
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{ SELECT rel_annot_id  FROM PolyprojectNGS.release_annotation where name=?; };
		return $sql;
	},
);

has sql_cmd_get_current_genes_releases_annotations => (
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{ SELECT * FROM PolyprojectNGS.release_gene; };
		return $sql;
	},
);

has sql_cmd_get_current_diag_project_releases_annotations => (
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{ SELECT name FROM PolyprojectNGS.release_annotation where diag=2; };
		return $sql;
	},
);

has sql_cmd_get_current_genome_project_releases_annotations => (
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{ SELECT name FROM PolyprojectNGS.release_annotation where genome=2; };
		return $sql;
	},
);

has sql_cmd_get_list_project_releases_annotations => (
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{ SELECT ra.name as annotation_release FROM PolyprojectNGS.project_release_annotation as pra, PolyprojectNGS.projects as p, PolyprojectNGS.release_annotation as ra where p.project_id=pra.project_id and pra.rel_annot_id=ra.rel_annot_id and p.project_id=?; };
		return $sql;
	},
);

has sql_cmd_get_list_releases_annotations => (
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{ SELECT ra.name as annotation_release FROM PolyprojectNGS.release_annotation as ra; };
		return $sql;
	},
);

has sql_cmd_getTranscriptsTooLongFromBundles =>(
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{ SELECT * FROM PolyprojectNGS.transcripts where LENGTH(ENSEMBL_ID) > 15; };
		return $sql
	},
);

has sql_cmd_getTranscriptsInfosTranscriptsTable =>(
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{ SELECT * FROM PolyprojectNGS.transcripts where ENSEMBL_ID=?; };
		return $sql
	},
);

has sql_cmd_getAllBundlesIdsFromTranscriptsId =>(
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{ SELECT * FROM PolyprojectNGS.bundle_transcripts where transcript_id=?; };
		return $sql
	},
);


has sql_cmd_list_capture_diag_from_transcrip_id => (
        is       => 'ro',
        
        lazy =>1,
        default => sub {
                my $self = shift;
                my $sql = qq{
                                SELECT cs.name  
                                    FROM PolyprojectNGS.transcripts as tr, PolyprojectNGS.bundle_transcripts as b_tr, PolyprojectNGS.capture_bundle as cb, PolyprojectNGS.bundle as b, PolyprojectNGS.capture_systems as cs
                                        where tr.ensembl_id=? and b_tr.transcript_id=tr.id and b_tr.bundle_id=cb.bundle_id and cb.capture_id=cs.capture_id and b.bundle_id=cb.bundle_id and cs.analyse not like "exome" and cs.analyse not like "genome";
                };
        return $sql
        },
);

has sql_cmd_transcrip_id_with_capture_diag => (
        is       => 'ro',
        lazy =>1,
        default => sub {
                my $self = shift;
                my $sql = qq{
								SELECT tr.ensembl_id as trans_id
                                    FROM PolyprojectNGS.transcripts as tr, PolyprojectNGS.bundle_transcripts as b_tr, PolyprojectNGS.capture_bundle as cb, PolyprojectNGS.bundle as b, PolyprojectNGS.capture_systems as cs
                                        where b_tr.transcript_id=tr.id and b_tr.bundle_id=cb.bundle_id and cb.capture_id=cs.capture_id and b.bundle_id=cb.bundle_id and cs.analyse not like "exome" and cs.analyse not like "genome";
                };
        return $sql
        },
);

has sql_cmd_list_all_capture_analyse => (
	is	 => 'ro',

	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
				select analyse from PolyprojectNGS.capture_systems;
		};
	return $sql
	},
);

has sql_cmd_list_ngs_patients_id_byCapture => (
	is	 => 'ro',

	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
				select distinct pat.patient_id from PolyprojectNGS.patient pat, PolyprojectNGS.projects proj, PolyprojectNGS.capture_systems capt
					where pat.capture_id=capt.capture_id and capt.analyse=? and pat.project_id=proj.project_id and proj.type_project_id=3 and not proj.name REGEXP 'TEST';
		};
	return $sql
	},
);

has sql_cmd_list_ngs_projects_name_byCapture => (
	is	 => 'ro',

	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
				select distinct proj.name from PolyprojectNGS.patient pat, PolyprojectNGS.projects proj, PolyprojectNGS.capture_systems capt
					where pat.capture_id=capt.capture_id and capt.name=? and pat.project_id=proj.project_id and proj.type_project_id=3 and not proj.name REGEXP 'TEST';
		};
	return $sql
	},
);

has sql_cmd_list_ngs_projects_name_byCaptureId => (
	is	 => 'ro',

	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
				select distinct proj.name from PolyprojectNGS.patient pat, PolyprojectNGS.projects proj, PolyprojectNGS.capture_systems capt
					where pat.capture_id=capt.capture_id and capt.capture_id=? and pat.project_id=proj.project_id and proj.type_project_id=3 and not proj.name REGEXP 'TEST';
		};
	return $sql
	},
);

has sql_cmd_list_ngs_projects_name_byAnalyse => (
	is	 => 'ro',

	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
				select distinct proj.name from PolyprojectNGS.patient pat, PolyprojectNGS.projects proj, PolyprojectNGS.capture_systems capt
					where pat.capture_id=capt.capture_id and capt.analyse=? and pat.project_id=proj.project_id and proj.type_project_id=3 and not proj.name REGEXP 'TEST';
		};
	return $sql
	},
);

has sql_cmd_list_projects_name => (
	is	 => 'ro',

	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
				select project_id,name,description from PolyprojectNGS.projects where type_project_id=3 order by project_id;
		};
	return $sql
	},
);

has sql_cmd_list_projects_name_exome_for_dejaVu => (
	is	 => 'ro',

	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
				select distinct (pr.name) as name from PolyprojectNGS.projects as pr ,PolyprojectNGS.capture_systems as cs, PolyprojectNGS.patient as p where p.project_id = pr.project_id and pr.dejaVu=1 AND p.capture_id = cs.capture_id and (cs .analyse="exome" or cs.analyse = "genome") ;;
		};
	return $sql
	},
);

has sql_cmd_list_projects_name_for_dejaVu => (
	is	 => 'ro',

	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
				select name from PolyprojectNGS.projects where type_project_id=3 and dejaVu=1;
		};
	return $sql
	},
);

has sql_cmd_hash_ngs_patients_infos => (
	is	 => 'ro',

	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
				select pat.patient_id, proj.name as project, pat.name, pat.family, pat.father, pat.mother, pat.sex, pat.status, pat.capture_id, cs.analyse from PolyprojectNGS.patient pat, PolyprojectNGS.projects proj, PolyprojectNGS.capture_systems cs
					where pat.project_id=proj.project_id and proj.type_project_id =3 and pat.capture_id=cs.capture_id and proj.dejaVu=1;
		};
	return $sql
	},
);

has sql_cmd_list_ngs_patients_id => (
	is	 => 'ro',

	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
				select patient_id from PolyprojectNGS.patient pat, PolyprojectNGS.projects proj
					where pat.project_id=proj.project_id and proj.type_project_id =3 and not proj.name REGEXP 'TEST';
		};
	return $sql
	},
);

has sql_cmd_list_ngs_patients_id_byYear => (
	is	 => 'ro',

	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
				select pat.patient_id, pat.creation_date, proj.name as project_name from PolyprojectNGS.patient pat, PolyprojectNGS.projects proj
					where pat.project_id=proj.project_id and proj.type_project_id=3 and not proj.name REGEXP 'TEST';
		};
	return $sql
	},
);

has sql_cmd_owner_project => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
		select BU.email as email, BU.PRENOM_U as firstname,BU.LOGIN as login FROM PolyprojectNGS.user_projects U, PolyprojectNGS.projects O, bipd_users.USER BU
		where U.project_id=O.project_id and O.project_id=? and U.USER_ID=BU.USER_ID and BU.equipe_id != 6 ;
	};
	return $sql
	},
);

has sql_cmd_owner_project_by_name => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
		select BU.email as email, BU.PRENOM_U as firstname FROM PolyprojectNGS.user_projects U, PolyprojectNGS.projects O, bipd_users.USER BU
		where U.project_id=O.project_id and O.name=? and U.USER_ID=BU.USER_ID and BU.equipe_id != 6 ;
	};
	return $sql
	},
);

has sql_list_all_projects => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $sql = qq{select distinct(project_id) as id, name, description from PolyprojectNGS.projects where name regexp 'NGS';};
		return $sql
	},
);

has sql_list_project_for_user =>(
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $sql = qq{	  	
			select distinct(o.project_id) as id,o.name as name, pt.name as type , db.name as dbname ,o.description as description, BU.EQUIPE_ID team, o.validation_db,o.creation_date 
			from PolyprojectNGS.projects o , PolyprojectNGS.databases_projects dp,PolyprojectNGS.polydb db,  PolyprojectNGS.user_projects  up ,  bipd_users.USER BU, PolyprojectNGS.project_types pt
			  where 
			  	up.project_id=o.project_id and (
			  	(up.USER_ID=BU.USER_ID AND BU.LOGIN=? and BU.password_txt=password(?)) or 
			  	( BU.LOGIN=? and BU.password_txt=password(?) and BU.EQUIPE_ID=6 ) )
			   and dp.project_id=o.project_id and db.db_id=dp.db_id  
			   and o.type_project_id = pt.type_project_id order by o.creation_date asc;
	};
	return $sql
	},
);

has sql_list_project_for_group =>(
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $sql = qq{
			SELECT distinct(O.project_id) as id, O.name as name, pt.name as type, db.name as dbname, O.description as description, BU.EQUIPE_ID team, O.validation_db,O.creation_date,UGG.NAME as "group"  FROM bipd_users.UGROUP as UGG, bipd_users.UGROUP_USER as UGU, bipd_users.USER as BU, PolyprojectNGS.ugroup_projects as UGP, PolyprojectNGS.projects O, PolyprojectNGS.polydb db , PolyprojectNGS.project_types pt, PolyprojectNGS.databases_projects dp  
				where 
					BU.LOGIN=? 
		            and BU.password_txt=password(?) 
					and UGU.USER_ID = BU.USER_ID 
					and UGU.UGROUP_ID = UGP.UGROUP_ID
					and UGU.UGROUP_ID = UGG.UGROUP_ID
		            and UGP.PROJECT_ID = O.project_id 
					and dp.project_id = O.project_id 
		            and db.db_id = dp.db_id
		            and O.type_project_id = pt.type_project_id 
						order by pt.name, O.name;
		};
		return $sql
	},
);

has sql_list_projects =>(
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $sql = qq{	  	
			select o.name as name  
			from PolyprojectNGS.projects o where NAME like 'NGS%';
	};
	
	return $sql
	},
);


has sql_cmd_project_by_id => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $sql = qq{
			select name as name from PolyprojectNGS.projects where project_id=?;
		};
		return $sql;
	},
);

has sql_cmd_project_by_name => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
	my $self = shift;
	my $sqlCmd = "select o.project_id as id, o.name as name,o.creation_date as creation_date,o.description as description,o.validation_db as validation_db,  pt.type_project_id as projectTypeId, pt.name as projectType, db.name as dbname, r.name as version, o.somatic as is_somatic";
	$sqlCmd .= " from PolyprojectNGS.projects o , PolyprojectNGS.databases_projects dp, PolyprojectNGS.polydb db, PolyprojectNGS.project_release pr , PolyprojectNGS.releases r, PolyprojectNGS.project_types pt";
	$sqlCmd .= " where o.name= ? ";
	$sqlCmd .= " and dp.project_id=o.project_id and db.db_id=dp.db_id";
	$sqlCmd .= " and pr.project_id=o.project_id and pr.release_id = r.release_id and pr.default=1 and o.type_project_id=pt.type_project_id;";
	return $sqlCmd;
	},
);
has sql_cmd_sequencing_machine =>(
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
	my $sql = qq{
		SELECT distinct sm.name as name 
    		FROM PolyprojectNGS.projects pr, PolyprojectNGS.patient pa, PolyprojectNGS.run_machine rm, PolyprojectNGS.sequencing_machines sm 
        		where pr.project_id= ? and pr.project_id=pa.project_id and pa.run_id=rm.run_id and rm.machine_id=sm.machine_id;
	};
	return $sql
	},
);

has sql_get_patients =>(
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
	my $sql = qq{
		select * from PolyprojectNGS.patient g where project_id=?  
	};
	return $sql
	},
);

has sql_get_groups =>(
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
	my $sql = qq{
		select g.*,patient.patient_id,patient.name as pname from PolyprojectNGS.patient,PolyprojectNGS.patient_groups as pg ,PolyprojectNGS.group as g where project_id=? and patient.patient_id=pg.patient_id and pg.group_id=g.group_id
	};
	return $sql
	},
);

sub getPatient {
	my ($self, $project_id, $genbo_id) = @_;
	confess();
	my $dbh = $self->getDbh();
	my $config = $self->getConfig();
	my $cmd = "select * from PolyprojectNGS.patient where project_id=$project_id and genbo_id=$genbo_id sort by creation_date asc";
	my $query = qq{ $cmd };
	my $sth = $dbh->prepare($query);
	$sth->execute();
	my $s = $sth->fetchall_arrayref({});
	return $s;
}

sub getGenboIdPatient {
	my ($self, $project_id, $patient_name) = @_;
	confess();
	my $dbh = $self->getDbh();
	my $config = $self->getConfig();
	my $cmd = "select genbo_id from PolyprojectNGS.patient where project_id=$project_id and name='$patient_name'";
	my $query = qq{ $cmd };
	my $sth = $dbh->prepare($query);
	$sth->execute();
	my $s = $sth->fetchall_hashref('genbo_id');
	my @tmp = keys(%$s);
	return int($tmp[0]);
}

has sql_get_origin_methods =>(
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
	my $sql = qq{
		SELECT distinct m.name as methodname 
    		FROM PolyprojectNGS.projects pr, PolyprojectNGS.patient pa, PolyprojectNGS.run r, PolyprojectNGS.run_methods rm, PolyprojectNGS.methods m 
        		where pr.project_id=? and pr.project_id=pa.project_id and pa.run_id=rm.run_id and rm.method_id=m.method_id and m.type=?;
	};
	return $sql
	},
);
has sql_get_similar_projects =>(
is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
	my $sql = qq{SELECT pr.name as name FROM PolyprojectNGS.patient p,PolyprojectNGS.projects pr  where p.capture_id = ? and p.project_id = pr.project_id group by pr.project_id};

	return $sql
	},
);

has sql_get_similar_projects_by_validation_db =>(
is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
	my $sql = qq{SELECT pr.name as name ,p.patient_id as patient_id FROM PolyprojectNGS.patient p,PolyprojectNGS.projects pr ,PolyprojectNGS.capture_systems cs  where cs.validation_db = ? and p.capture_id=cs.capture_id  and p.project_id = pr.project_id group by pr.project_id};
	return $sql
	},
);
has sql_get_similar_projects_by_analyse =>(
is		=> 'ro',
	
	lazy =>1,
	
	default	=> sub {
		my $sql = qq{
				select distinct (pr.name) as name from PolyprojectNGS.projects as pr ,PolyprojectNGS.capture_systems as cs, PolyprojectNGS.patient as p where p.project_id = pr.project_id and pr.dejaVu=1 AND p.capture_id = cs.capture_id and cs .analyse=?  ;
		};
	return $sql;#return $sql
	},
);
has sql_get_methods =>(
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
	my $sql = qq{
		select m.* from PolyprojectNGS.patient p  ,PolyprojectNGS.methods m , PolyprojectNGS.patient_methods pm 
		where p.name= ? and p.project_id= ? and pm.patient_id = p.patient_id and m.method_id=pm.method_id  and m.type rlike ?;
	};
	return $sql
	},
);

has sql_cmd_select_run =>(
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
	my $sql = qq{
		select distinct(run_id) as id  from PolyprojectNGS.patient where project_id=?;
	};
	return $sql
	},
);

#has sql_panel_all_names => (
#	is		=> 'ro',
#	
#	lazy =>1,
#	default	=> sub {
#		my $sql = qq{ SELECT name FROM PolyPanelDB.panel; };
#		return $sql;
#	},
#);

#has sql_panel_all_ids => (
#	is		=> 'ro',
#	
#	lazy =>1,
#	default	=> sub {
#		my $sql = qq{ SELECT panel_id FROM PolyPanelDB.panel; };
#		return $sql;
#	},
#);




has sql_panel_infos => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $sql = qq{ SELECT * FROM PolyPanelDB.panel where panel_id=?; };
		return $sql;
	},
);

has sql_capture_infos =>(
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
	my $sql = qq{SELECT c.*,pr.name as gencode_version FROM  PolyprojectNGS.capture_systems c, PolyprojectNGS.release_gene pr where c.capture_id=? and c.rel_gene_id = pr.rel_gene_id ;};

	return $sql
	},
);

has sql_umi =>(
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
	my $sql = qq{SELECT u.name as name , u.mask as mask FROM  PolyprojectNGS.capture_systems c, PolyprojectNGS.umi u where c.capture_id=? and c.umi_id = u.umi_id ;};

	return $sql;
	},
);



has sql_capture_infos_by_name =>(
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
	my $sql = qq{SELECT c.capture_id FROM  PolyprojectNGS.capture_systems c where c.name=?;};

	return $sql
	},
);

has sql_allprojects_infos_for_user =>(
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
	my $sql = qq{select p.name as name,analyse,pa.name,p.validation_db as type from 
					PolyprojectNGS.projects p ,PolyprojectNGS.patient pa , PolyprojectNGS.run r ,PolyprojectNGS.capture_systems c, 
					PolyprojectNGS.user_projects pu , bipd_users.USER u
					where p.project_id = pa.project_id and pa.run_id=r.run_id and pa.capture_id = c.capture_id and u.login = ? and p.project_id=pu.PROJECT_ID and u.USER_ID=pu.USER_ID;
	};
	return $sql
	},
);
has sql_allprojects_infos =>(
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
	my $sql = qq{select p.name as name,analyse,pa.name,p.validation_db as type from PolyprojectNGS.projects p ,PolyprojectNGS.patient pa , PolyprojectNGS.run r ,PolyprojectNGS.capture_systems c where p.project_id = pa.project_id and pa.run_id=r.run_id and pa.capture_id = c.capture_id;
	};
	return $sql
	},
);

has sql_cmd_capture_transcripts =>(
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
	my $sql = qq{select b.*,tr.* from  PolyprojectNGS.capture_bundle cb, PolyprojectNGS.bundle b, PolyprojectNGS.bundle_transcripts bt,PolyprojectNGS.transcripts tr where cb.capture_id = ? and cb.bundle_id = b.bundle_id and bt.bundle_id = cb.bundle_id and bt.transcript_id = tr.ID and tr.deprecated=0 order by tr.GENE;	};

	return $sql
	},
);

has sql_release_gene_capture_id =>(
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
	my $sql = qq{ SELECT rg.name as release_id  
		FROM PolyprojectNGS.release_gene as rg, PolyprojectNGS.capture_systems as cp, 
		PolyprojectNGS.capture_bundle as cb , PolyprojectNGS.bundle_release_gene  as br where
		  rg.rel_gene_id = br.rel_gene_id and br.bundle_id = cb.bundle_id and cp.capture_id=cb.capture_id
		    and cp.capture_id=?;
	};
	return $sql;
	}
);

has sql_cmd_capture_transcripts_by_name =>(
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		
	my $sql = qq{SELECT t.* FROM 
PolyprojectNGS.capture_systems as c , 
PolyprojectNGS.capture_bundle as cb,
PolyprojectNGS.transcripts as t, 
PolyprojectNGS.bundle_transcripts as bt
 where c.name = ? and cb.capture_id = c.capture_id and cb.bundle_id = bt.bundle_id  and bt.transcript_id=t.ID;};

	return $sql
	},
);

has sql_cmd_getBundleGenesAllProjectsUsers => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $sql = qq{
			SELECT distinct capture_systems.name as capture, bundle.bundle_id as id, bundle.name as name, bundle.description as description
			    FROM PolyprojectNGS.user_projects, bipd_users.USER, PolyprojectNGS.patient, PolyprojectNGS.capture_bundle, PolyprojectNGS.bundle, PolyprojectNGS.bundle_transcripts, PolyprojectNGS.transcripts, PolyprojectNGS.capture_systems
			        where
			            USER.LOGIN=?
			            and user_projects.USER_ID=USER.USER_ID
			            and patient.project_id=user_projects.PROJECT_ID
			            and patient.capture_id=capture_bundle.capture_id
			            and capture_bundle.bundle_id=bundle.bundle_id
			            and capture_bundle.capture_id=capture_systems.capture_id
			            and bundle_transcripts.bundle_id=capture_bundle.bundle_id
			            and transcripts.ID=bundle_transcripts.transcript_id;
		};
		return $sql;
	},
);

has sql_cmd_getGenesNamesInBundle => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $sql = qq{
			SELECT distinct transcripts.GENE as gene_name
			    FROM PolyprojectNGS.bundle_transcripts, PolyprojectNGS.transcripts
			        where
			            bundle_transcripts.bundle_id=?
			            and transcripts.ID=bundle_transcripts.transcript_id;
		};
		return $sql;
	},
);

has sql_cmd_getGenesNamesAllBundle => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $sql = qq{
			SELECT distinct capture_systems.name as capture, bundle.bundle_id as id, bundle.name as name, bundle.description as description, transcripts.GENE as gene_name, transcripts.ENSEMBL_ID as transcript_id
			    FROM PolyprojectNGS.capture_bundle, PolyprojectNGS.bundle, PolyprojectNGS.bundle_transcripts, PolyprojectNGS.transcripts, PolyprojectNGS.capture_systems
			        where
			            capture_bundle.bundle_id=bundle.bundle_id
			            and capture_bundle.capture_id=capture_systems.capture_id
			            and bundle_transcripts.bundle_id=capture_bundle.bundle_id
			            and transcripts.ID=bundle_transcripts.transcript_id
                        and capture_systems.name!='agilent_50_v5';
		};
		return $sql;
	},
);

has sql_cmd_getOmimGenesTranscriptsNames => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $sql = qq{
			SELECT distinct capture_systems.name as capture, bundle.bundle_id as id, bundle.name as name, bundle.description as description, transcripts.GENE as gene_name, transcripts.ENSEMBL_ID as transcript_id
			    FROM PolyprojectNGS.capture_bundle, PolyprojectNGS.bundle, PolyprojectNGS.bundle_transcripts, PolyprojectNGS.transcripts, PolyprojectNGS.capture_systems
			        where
			            capture_bundle.bundle_id=bundle.bundle_id
			            and capture_bundle.capture_id=capture_systems.capture_id
			            and bundle_transcripts.bundle_id=capture_bundle.bundle_id
			            and transcripts.ID=bundle_transcripts.transcript_id
                        and (capture_systems.name = 'Omim' or capture_systems.name = 'Omim_hg38');
		};
		return $sql;
	},
);

sub getBundleTranscripts {
	my ($self,$project_id) = @_;
	my $dbh = $self->getDbh;
	my $sql = qq{select b.*,tr.* from  PolyprojectNGS.project_bundle cb, PolyprojectNGS.bundle b, PolyprojectNGS.bundle_transcripts bt,PolyprojectNGS.transcripts tr where cb.project_id = ? and cb.bundle_id = b.bundle_id and bt.bundle_id = cb.bundle_id and bt.transcript_id = tr.ID and tr.deprecated=0 order by tr.GENE;	};
	my $sth = $dbh->prepare($sql );
	$sth->execute($project_id);
	my $res = $sth->fetchall_arrayref({});
	my $res2 = {};
	$res2->{bundle} = {};
	$res2->{transcripts} = {};
	my $tr;
	
	
	foreach my $h (@$res){
		push(@{$res2->{bundle}->{$h->{name}}},$h); 
		push (@{$res2->{transcripts}->{$h->{ENSEMBL_ID}}},$h->{name});
		push (@{$res2->{transcripts_name}},$h->{ENSEMBL_ID});
	}
	return $res2;
}
	
sub isDejaVu{
	my ($self,$id)=@_;
	my $sql = qq{select nb_project, count(distinct g) as nb_patient, sum(he) as he , sum(ho) as ho from  (select  d.nb as nb_project, d.genbo_id as did, g.genbo_id  as g , if (ra.he>=1,1,0) as he, if (ra.ho>=1,1,0) as ho
from DEJAVU_STATIC d, RELATION r,GENBO g, RELATION_ANNEX ra
 where d.NB >1 and d.GENBO_ID = $id and r.genbo2_id= d.genbo_id and g.genbo_id=r.genbo_id
  and g.type_genbo_id = 5
  and ra.relation_id=r.relation_id  group by ra.relation_id,d.genbo_id ) as tbl ;
	};
	#my $sql = qq{select  d.nb as nb_project ,count(g.genbo_id) as nb_patient, sum(ra.hne) as he , sum(ra.ho) as ho  from DEJAVU_STATIC d, RELATION r,GENBO g, RELATION_ANNEX ra where d.NB >1 and d.GENBO_ID = $id and r.genbo2_id= d.genbo_id and g.genbo_id=r.genbo_id and g.type_genbo_id = 5  and ra.relation_id=r.relation_id and ra.he>=0 and ra.ho>=0 group by d.genbo_id;};
	my $sth =  $self->getDbh()->prepare($sql);
	#warn $sql;
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

sub DejaVuForOneProject{
	my ($self,$pid,$id)=@_;
	confess() unless $pid;
	my $sql = qq{select nb_project, count(distinct g) as nb_patient, sum(he) as he , sum(ho) as ho from  (select  d.nb as nb_project, d.genbo_id as did, g.genbo_id  as g , if (ra.he>=1,1,0) as he, if (ra.ho>=1,1,0) as ho
from DEJAVU_STATIC d, RELATION r,GENBO g, RELATION_ANNEX ra
 where d.NB >1 and d.GENBO_ID = $id and r.genbo2_id= d.genbo_id and g.genbo_id=r.genbo_id
  and g.type_genbo_id = 5
  and g.origin_id=$pid
  and ra.relation_id=r.relation_id  group by ra.relation_id,d.genbo_id ) as tbl ;
	};
	my $sth =  $self->getDbh()->prepare($sql);
	#warn $sql;
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

has sql_cmd_getPatientProjectInfo => (
        is       => 'ro',
        
        lazy =>1,
        default => sub {
                my $self = shift;
		        my $query = qq{
		                SELECT distinct a.project_id,a.run_id,a.name,a.patient_id,
			                a.project_id_dest,a.family,a.father,a.mother,a.sex,a.status,a.bar_code,a.flowcell,
			                r.description as desRun, r.document,r.name as nameRun,
			                r.file_name as FileName,r.file_type as FileType,                
			                S.name as macName,
			        		GROUP_CONCAT(DISTINCT f.name ORDER BY f.name DESC SEPARATOR ' ') as 'plateformName',
			                ms.name as methSeqName,
			                C.name as capName,
			                a.creation_date as cDate
			                FROM PolyprojectNGS.patient a,
				                PolyprojectNGS.run r,
				                PolyprojectNGS.run_machine rm,
				                PolyprojectNGS.sequencing_machines S,
				                PolyprojectNGS.run_plateform rp,
				                PolyprojectNGS.plateform f,
				                PolyprojectNGS.run_method_seq rs,
				                PolyprojectNGS.method_seq ms,
				                PolyprojectNGS.capture_systems C
				                where a.project_id=?
					                and a.run_id=r.run_id
					                and r.run_id=rm.run_id
					                and rm.machine_id=S.machine_id
					                and r.run_id=rp.run_id
					                and rp.plateform_id=f.plateform_id
					                and r.run_id=rs.run_id
					                and rs.method_seq_id=ms.method_seq_id
					                and a.capture_id=C.capture_id
					                and a.control <> 1 
					        		group by a.patient_id
					                order by a.run_id;      
		        };
				return $query;
        },
);

has sql_cmd_getPatientSomaticProjectInfo => (
        is       => 'ro',
        
        lazy =>1,
        default => sub {
                my $self = shift;
		        my $query = qq{
		                SELECT distinct a.project_id,a.run_id,a.name,a.patient_id,
			                a.project_id_dest,a.family,a.father,a.mother,a.sex,a.status,a.bar_code,a.flowcell,
			                r.description as desRun, r.document,r.name as nameRun,
			                r.file_name as FileName,r.file_type as FileType,                
			                S.name as macName,
			        		GROUP_CONCAT(DISTINCT f.name ORDER BY f.name DESC SEPARATOR ' ') as 'plateformName',
			                ms.name as methSeqName,
			                C.name as capName,
			                a.creation_date as cDate,
			                g.name as somatic_group
			                FROM PolyprojectNGS.patient a,
				                PolyprojectNGS.run r,
				                PolyprojectNGS.run_machine rm,
				                PolyprojectNGS.sequencing_machines S,
				                PolyprojectNGS.run_plateform rp,
				                PolyprojectNGS.plateform f,
				                PolyprojectNGS.run_method_seq rs,
				                PolyprojectNGS.method_seq ms,
				                PolyprojectNGS.capture_systems C,
				                PolyprojectNGS.patient_groups pg,
				                PolyprojectNGS.group g
				                where a.project_id=?
					                and a.patient_id=pg.patient_id
					                and pg.group_id=g.group_id
					                and a.run_id=r.run_id
					                and r.run_id=rm.run_id
					                and rm.machine_id=S.machine_id
					                and r.run_id=rp.run_id
					                and rp.plateform_id=f.plateform_id
					                and r.run_id=rs.run_id
					                and rs.method_seq_id=ms.method_seq_id
					                and a.capture_id=C.capture_id
					        		group by a.patient_id
					                order by a.run_id;      
		        };
				return $query;
        },
);

has sql_cmd_get_proj_ids_genes_datatbase_version => (
        is       => 'ro',
        
        lazy =>1,
        default => sub {
                my $self = shift;
		        my $query = qq{
		        	SELECT p.project_id, p.name, prg.rel_gene_id, prpd.version_id FROM PolyprojectNGS.projects as p, PolyprojectNGS.project_release_gene as prg, PolyprojectNGS.project_release_public_database as prpd
						where p.project_id=prg.project_id and p.project_id=prpd.project_id;
		        };
				return $query;
        },
);

has sql_cmd_get_alamut_api_key_from_username => (
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{
			SELECT USER.LOGIN as login, USER.EQUIPE_ID as equipe_id, ALAMUT_LICENCE.NAME as licence_alamut, ALAMUT_LICENCE.API_KEY as api_key, ALAMUT_LICENCE.INSTITUTION_KEY as institution_key  FROM bipd_users.USER, bipd_users.EQUIPE_ALAMUT_LICENCE, bipd_users.ALAMUT_LICENCE
				where USER.login=? and USER.EQUIPE_ID=EQUIPE_ALAMUT_LICENCE.EQUIPE_ID and EQUIPE_ALAMUT_LICENCE.ALAMUT_ID=ALAMUT_LICENCE.ALAMUT_ID;
		};
		return $sql;
	},
);

has sql_cmd_get_projects_control_giab => (
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{
			SELECT name FROM PolyprojectNGS.projects where validation_db='control_giab';
		};
		return $sql;
	},
);

has sql_cmd_get_quick_projects_list_RNA => (
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{
			SELECT p.project_id as id,  p.name as name, p.description, p.creation_date as cDate, po.name as dbname, r.name as relname,
				GROUP_CONCAT(DISTINCT pp.version_id ORDER BY pp.version_id  SEPARATOR ' ') as 'ppversionid',
				GROUP_CONCAT(DISTINCT rg.name ORDER BY rg.name  SEPARATOR ' ') as 'relGene',
				GROUP_CONCAT(DISTINCT U.login ORDER BY U.login DESC SEPARATOR ',') as 'username'
				FROM
					PolyprojectNGS.projects p 
					LEFT JOIN PolyprojectNGS.databases_projects dp
					ON p.project_id =dp.project_id
					LEFT JOIN PolyprojectNGS.polydb po
					ON dp.db_id = po.db_id
					LEFT JOIN PolyprojectNGS.project_release pr
					ON p.project_id=pr.project_id
					LEFT JOIN PolyprojectNGS.releases r
					ON pr.release_id=r.release_id
			        LEFT JOIN PolyprojectNGS.project_release_public_database pp
			        ON p.project_id = pp.project_id
			        LEFT JOIN PolyprojectNGS.project_release_gene pg
			        ON p.project_id = pg.project_id
			        LEFT JOIN PolyprojectNGS.release_gene rg
			        ON rg.rel_gene_id=pg.rel_gene_id
			        LEFT JOIN PolyprojectNGS.user_projects up
					ON p.project_id = up.project_id
					LEFT JOIN bipd_users.`USER` U
					ON up.user_id = U.user_id
			        LEFT JOIN PolyprojectNGS.ugroup_projects gp
					ON p.project_id = gp.project_id
			        LEFT JOIN bipd_users.UGROUP ug
			        ON gp.ugroup_id=ug.ugroup_id
					LEFT JOIN bipd_users.EQUIPE E
					ON U.equipe_id = E.equipe_id
					LEFT JOIN bipd_users.UNITE T
					ON E.unite_id = T.unite_id
					WHERE p.type_project_id=3 and dp.db_id !=2
			        GROUP BY p.project_id;
		};
		return $sql;
	},
);

has sql_cmd_check_project_RNA=> (
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{
			SELECT p.project_id as id,  p.name as name, p.description, p.creation_date as cDate, po.name as dbname, r.name as relname,
				GROUP_CONCAT(DISTINCT pp.version_id ORDER BY pp.version_id  SEPARATOR ' ') as 'ppversionid',
				GROUP_CONCAT(DISTINCT rg.name ORDER BY rg.name  SEPARATOR ' ') as 'relGene',
				GROUP_CONCAT(DISTINCT U.login ORDER BY U.login DESC SEPARATOR ',') as 'username'
				FROM
					PolyprojectNGS.projects p 
					LEFT JOIN PolyprojectNGS.databases_projects dp
					ON p.project_id =dp.project_id
					LEFT JOIN PolyprojectNGS.polydb po
					ON dp.db_id = po.db_id
					LEFT JOIN PolyprojectNGS.project_release pr
					ON p.project_id=pr.project_id
					LEFT JOIN PolyprojectNGS.releases r
					ON pr.release_id=r.release_id
			        LEFT JOIN PolyprojectNGS.project_release_public_database pp
			        ON p.project_id = pp.project_id
			        LEFT JOIN PolyprojectNGS.project_release_gene pg
			        ON p.project_id = pg.project_id
			        LEFT JOIN PolyprojectNGS.release_gene rg
			        ON rg.rel_gene_id=pg.rel_gene_id
			        LEFT JOIN PolyprojectNGS.user_projects up
					ON p.project_id = up.project_id
					LEFT JOIN bipd_users.`USER` U
					ON up.user_id = U.user_id
			        LEFT JOIN PolyprojectNGS.ugroup_projects gp
					ON p.project_id = gp.project_id
			        LEFT JOIN bipd_users.UGROUP ug
			        ON gp.ugroup_id=ug.ugroup_id
					LEFT JOIN bipd_users.EQUIPE E
					ON U.equipe_id = E.equipe_id
					LEFT JOIN bipd_users.UNITE T
					ON E.unite_id = T.unite_id
					WHERE p.type_project_id=3 and dp.db_id !=2 and p.name=?
			        GROUP BY p.project_id;
		};
		return $sql;
	},
);

has sql_cmd_get_quick_patients_list_from_project_id => (
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{
			SELECT
		        a.patient_id,
		        a.name as name,
		        a.type as type,
		        r.run_id,
				C.name as capName,
				C.analyse as capAnalyse,
		        GROUP_CONCAT(DISTINCT S.name ORDER BY S.name ASC SEPARATOR ' ') as 'macName',
		        GROUP_CONCAT(DISTINCT f.name ORDER BY f.name ASC SEPARATOR ' ') as 'plateformName',
		        GROUP_CONCAT(DISTINCT ms.name ORDER BY ms.name ASC SEPARATOR ' ') as 'methSeqName',
		        GROUP_CONCAT(DISTINCT case M.type when 'ALIGN' THEN M.name ELSE NULL END ORDER BY M.name ASC SEPARATOR ' ') as 'methAln',
		        GROUP_CONCAT(DISTINCT case M.type when 'SNP' THEN M.name ELSE NULL END ORDER BY M.name ASC SEPARATOR ' ') as 'methCall'
		        FROM PolyprojectNGS.patient a
		        LEFT JOIN PolyprojectNGS.run r
		        ON a.run_id = r.run_id
		        LEFT JOIN PolyprojectNGS.run_machine rm
		        ON r.run_id = rm.run_id
		        LEFT JOIN PolyprojectNGS.sequencing_machines S
		        ON rm.machine_id=S.machine_id
		        LEFT JOIN PolyprojectNGS.run_plateform rp
		        ON r.run_id=rp.run_id
		        LEFT JOIN PolyprojectNGS.plateform f
		        ON rp.plateform_id=f.plateform_id
		        LEFT JOIN PolyprojectNGS.run_method_seq rs
		        ON r.run_id=rs.run_id
		        LEFT JOIN PolyprojectNGS.method_seq ms
		        ON rs.method_seq_id=ms.method_seq_id
		        LEFT JOIN PolyprojectNGS.patient_methods pm
		        ON a.patient_id = pm.patient_id
		        LEFT JOIN PolyprojectNGS.methods M
		        ON pm.method_id=M.method_id
		        LEFT JOIN PolyprojectNGS.capture_systems C
		        ON a.capture_id=C.capture_id
				WHERE a.project_id = ?
				GROUP BY patient_id
		        ORDER BY name;
		};
		return $sql;
	},
);

has sql_cmd_get_projects_ids_with_patients_type_rna => (
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{
			SELECT distinct(project_id) FROM PolyprojectNGS.patient
				where type='rna'; 
		};
		return $sql;
	},
);

has sql_cmd_get_projects_ids_with_patients_type_rna_with_project_name => (
	is	 => 'ro',

	lazy => 1,
	default	=> sub {
		my $sql = qq{
			SELECT pr.project_id as id, pr.name as name, pr.creation_date as cDate FROM PolyprojectNGS.patient as pa, PolyprojectNGS.projects as pr
				where pa.type='rna' and pa.project_id=pr.project_id and pr.name=?;
		};
		return $sql;
	},
);

sub isLoginSTAFF {
	my ($self, $name) = @_;
	my $dbh = $self->getDbh();
	my $sql = qq{SELECT (EQUIPE_ID) FROM bipd_users.USER where login=?};
	my $sth = $dbh->prepare($sql);
	$sth->execute($name);
	my $res = $sth->fetchall_arrayref({});
	return 1 if ($res->[0]->{EQUIPE_ID} == 6);
	return;
}

sub getListProjectsRnaSeqFromLoginPwd {
	my ($self, $login, $pwd,  $project_name) = @_;
	
	my $is_BIPD_login = $self->isLoginSTAFF($login);
	my $h_found;
	
	my @lProj = @{$self->getListProjectsRnaSeq($project_name)};
	foreach my $hpr (@lProj) {
		$h_found->{$hpr->{name}} = undef;
	}
	if (not $is_BIPD_login) {
		my $res_group = $self->getProjectHashForGroup($login,$pwd);
		if ($res_group) {
			foreach my $project_id (keys %{$res_group}) {
				$res_group->{$project_id}->{username} = $login;
				push(@lProj, $res_group->{$project_id}) unless exists $h_found->{$project_id};
			}
		}
	}
	return \@lProj;
}

sub getListProjectsRnaSeq {
	my ($self, $project_name) = @_;
	my @l_res;
	my $dbh = $self->getDbh();
	my ($sql, $sql2);
	if ($project_name) {
		$sql = $self->sql_cmd_check_project_RNA();
		$sql2 = $self->sql_cmd_get_projects_ids_with_patients_type_rna_with_project_name();
	}
	else {
		$sql = $self->sql_cmd_get_quick_projects_list_RNA();
		$sql2 = $self->sql_cmd_get_projects_ids_with_patients_type_rna();
	}
	my $sth = $dbh->prepare($sql);
	my $sth2 = $dbh->prepare($sql2);
	if ($project_name) {
		$sth->execute($project_name);
		$sth2->execute($project_name);
	}
	else {
		$sth->execute();
		$sth2->execute();
	}
	my $h = $sth->fetchall_hashref('id');
	my $h2 = $sth2->fetchall_hashref('project_id');
	
	
	my @l_projects_ids = keys %$h;
	foreach my $pr_id (keys %$h2) {
		push(@l_projects_ids, $pr_id) unless (exists $h->{$pr_id});
	}
	
	foreach my $project_id (sort {$b <=> $a} @l_projects_ids) {
		next if ($project_id == 0);
		my $sql2 = $self->sql_cmd_get_quick_patients_list_from_project_id();
		my $sth2 = $dbh->prepare($sql2);
		$sth2->execute($project_id);
		my $h_patients = $sth2->fetchall_hashref('name');
		my ($is_rna, $h_captures);
		foreach my $patient_name (keys %$h_patients) {
			$is_rna = 1 if ($h_patients->{$patient_name}->{type} eq 'rna');
			$is_rna = 1 if ($h_patients->{$patient_name}->{capAnalyse} eq 'rnaseq');
			$h_captures->{$h_patients->{$patient_name}->{capName}} = undef;
		}
		next unless ($is_rna);
		my @l_versions = split(' ', $h->{$project_id}->{ppversionid});
		@l_versions = sort {$a <=> $b} @l_versions;
		my $max_annot = $self->getMaxPublicDatabaseVersion();
		$h->{$project_id}->{version} = abs($max_annot - $l_versions[-1]).'::'.$h->{$project_id}->{relGene}.'-'.$l_versions[-1];
		$h->{$project_id}->{samples} = scalar(keys %$h_patients);
		$h->{$project_id}->{patient_name} = join(';', keys %$h_patients);
		$h->{$project_id}->{capture_name} = join(';', keys %$h_captures);
		$h->{$project_id}->{creation_date} = $h->{$project_id}->{cDate};
		$h->{$project_id}->{genome} = $h->{$project_id}->{relname};
		$h->{$project_id}->{project_id} = $project_id;
		push(@l_res, $h->{$project_id});
	}
	return \@l_res;
}

sub getListProjectsControlGIAB {
	my ($self) = @_;
	my $dbh = $self->getDbh();
	my $config = $self->getConfig();
	my $sql = $self->sql_cmd_get_projects_control_giab();
	my $sth = $dbh->prepare($sql);
	$sth->execute();
	my @l = keys %{$sth->fetchall_hashref('name')};
	return \@l;
}

sub getAlamutApiKeyFromUserName {
	my ($self, $user_name) = @_;
	my $dbh = $self->getDbh();
	my $config = $self->getConfig();
	my $sql = $self->sql_cmd_get_alamut_api_key_from_username();
	my $sth = $dbh->prepare($sql);
	$sth->execute($user_name);
	my $h = $sth->fetchall_hashref('login');
	return $h->{$user_name};
}

has sql_cmd_get_runid_from_samplename => (
	is       => 'ro',
	
	lazy =>1,
	default => sub {
		my $self = shift;
		my $query = qq{
			SELECT * FROM PolyprojectNGS.patient where name=?;
		};
		return $query;
	},
);

sub getHashRunIdFromSampleName {
	my ($self, $sample_name) = @_;
	my $dbh = $self->getDbh();
	my $sql = $self->sql_cmd_get_runid_from_samplename();
	my $sth = $dbh->prepare($sql);
	$sth->execute($sample_name);
	my $h = $sth->fetchall_hashref('run_id');
	return $h;
}

has sql_cmd_get_runid_infos => (
	is       => 'ro',
	
	lazy =>1,
	default => sub {
		my $self = shift;
		my $query = qq{
			SELECT * FROM PolyprojectNGS.run where run_id=?;
		};
		return $query;
	},
);

sub getHashRunIdInfos {
	my ($self, $run_id) = @_;
	my $dbh = $self->getDbh();
	my $sql = $self->sql_cmd_get_runid_infos();
	my $sth = $dbh->prepare($sql);
	$sth->execute($run_id);
	my $h = $sth->fetchall_hashref('run_id');
	return $h;
}


1;