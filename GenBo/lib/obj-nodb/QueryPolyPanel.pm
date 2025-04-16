package QueryPolyPanel;
use strict;
use Vcf;
use Moo;

use Data::Dumper;
use Config::Std;
use Carp;

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
has database => (
is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return "PolyPanelDB";
	}
);
has all_config => (
	is		=> 'ro',
	required=> 1,
);

has dbh => (
	is		=> 'ro',
	required=> 1,
	reader	=> 'getDbh',
);

###########
#SQL 
###########

##########
#SQL PANEL 
#########
has sql_panel_all_ids => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $db = $self->database();
		my $sql = qq{ SELECT panel_id FROM $db.panel where current=1; };
		return $sql;
	},
);

has sql_panel_all_names => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $db = $self->database();
		my $sql = qq{ SELECT name FROM $db.panel where current=1; };
		return $sql;
	},
);

has sql_panel_infos => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $db = $self->database();
		my $sql = qq{
			SELECT p.panel_id, p.name, p.description, p.source, p.version, p.current, p.validation_db, p.source_id, s.name as source, s.link, p.creation_date, p.update_date, p.creator, p.creator_email 
				FROM PolyPanelDB.panel p, PolyPanelDB.source s  
					where p.panel_id=? and p.source_id=s.source_id;
		 };
		return $sql;
	},
);

has sql_panel_by_gene_name => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $db = $self->database();
		my $sql = qq{SELECT panel.name as pname, panel.panel_id as panel_id, phenotype.name as phenotype, phenotype.phenotype_id as phenotype_id  
		FROM $db.gene as gene,$db.gene_bundle ,$db.panel_bundle,$db.panel,PolyPhenotypeDB.phenotype_panel,PolyPhenotypeDB.phenotype
		 where gene.name=? and gene.gene_id = gene_bundle.gene_id and gene_bundle.bundle_id=panel_bundle.bundle_id and panel_bundle.panel_id = panel.panel_id and phenotype_panel.panel_id =panel.panel_id and phenotype_panel.phenotype_id = phenotype.phenotype_id group by panel.panel_id};
		return $sql;
	},
);

has sql_genes_by_panel_id => (
	is		=> 'ro',
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $db = $self->database();
		my $sql = qq{
			 SELECT ensg.ensg_vid as ensg, gene.name as name
				FROM PolyPanelDB.gene, PolyPanelDB.gene_bundle, PolyPanelDB.panel_bundle,PolyPanelDB.panel, PolyPanelDB.ensg, PolyPanelDB.gene_ensg
			 		where gene_ensg.gene_id=gene.gene_id and gene_ensg.ensg_id=ensg.ensg_id and gene.gene_id=gene_bundle.gene_id and gene_bundle.bundle_id=panel_bundle.bundle_id and panel_bundle.panel_id=panel.panel_id and panel.panel_id=? and panel.current=1;
		};
		return $sql;
	},
);

has sql_genes_by_panel_name => (
	is		=> 'ro',
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $db = $self->database();
		my $sql = qq{
			 SELECT ensg.ensg_vid as ensg, gene.name as name
				FROM PolyPanelDB.gene, PolyPanelDB.gene_bundle, PolyPanelDB.panel_bundle,PolyPanelDB.panel, PolyPanelDB.ensg, PolyPanelDB.gene_ensg
			 		where gene_ensg.gene_id=gene.gene_id and gene_ensg.ensg_id=ensg.ensg_id and gene.gene_id=gene_bundle.gene_id and gene_bundle.bundle_id=panel_bundle.bundle_id and panel_bundle.panel_id=panel.panel_id and panel.name=? and panel.current=1;
		};
		return $sql;
	},
);

has sql_count_genes_by_panel_name => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $db = $self->database();
		my $sql = qq{
			 SELECT count(DISTINCT gene.gene_id) as nb
				FROM PolyPanelDB.gene, PolyPanelDB.gene_bundle, PolyPanelDB.panel_bundle,PolyPanelDB.panel, PolyPanelDB.gene_ensg
			 		where  gene_ensg.gene_id=gene.gene_id and gene.gene_id=gene_bundle.gene_id and gene_bundle.bundle_id=panel_bundle.bundle_id and panel_bundle.panel_id=panel.panel_id and panel.panel_id=? and panel.current=1;
		};
		return $sql;
	},
);



has sql_panel_by_gene_id => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $db = $self->database();
		my $sql = qq{SELECT panel.name as pname
		FROM $db.gene as gene,$db.gene_bundle ,$db.panel_bundle,$db.panel  
		 where gene.name=? and gene.gene_id = gene_bundle.gene_id and gene_bundle.bundle_id=panel_bundle.bundle_id and panel_bundle.panel_id = panel.panel_id};
		return $sql;
	},
);
###########
#SQL BUNDLE 
###########
has sql_bundle_ids_from_panel_id => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
			my $self = shift;
		my $db = $self->database();
		my $sql = qq{ SELECT b.* FROM $db.panel as p, $db.bundle as b, $db.panel_bundle as pb where p.panel_id=pb.panel_id and pb.bundle_id=b.bundle_id and p.panel_id=?; };
		return $sql;
	},
);

has sql_bundle_infos => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
			my $self = shift;
		my $db = $self->database();
		my $sql = qq{ SELECT * FROM $db.bundle where bundle_id=?; };
		return $sql;
	},
);

has sql_bundle_genes => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $db = $self->database();
		my $sql = qq{ SELECT g.* FROM $db.bundle as b, $db.gene_bundle as gb, $db.gene as g where b.bundle_id=gb.bundle_id and gb.gene_id=g.gene_id and b.bundle_id=?; };
		return $sql;
	},
);
has sql_bundle_genes_ensg => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $db = $self->database();
		my $sql = qq{ SELECT ensg.ensg_vid as gid FROM $db.bundle as b, $db.gene_bundle as gb, $db.gene as g,$db.gene_ensg as gensg, $db.ensg  where b.bundle_id=gb.bundle_id and gb.gene_id=g.gene_id and b.bundle_id=? and gensg.gene_id=g.gene_id and ensg.ensg_id = gensg.ensg_id;; };
		return $sql;
	},
);

has sql_genes_omim_morbid => (
	is		=> 'ro',
	
	lazy =>1,
	default	=> sub {
		my $self = shift;
		my $db = $self->database();
		my $sql = qq{ SELECT g.name FROM PolyPanelDB.bundle as b, PolyPanelDB.gene_bundle as gb, PolyPanelDB.gene as g where b.name="Omim-morbid" and b.bundle_id=gb.bundle_id and gb.gene_id=g.gene_id; };
		return $sql;
	},
);


##END SQL BUNDLE##


##### METHODS FOR BUNDLE #####

sub getBundlesIdsFromPanelId {
	my ($self, $panel_id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_bundle_ids_from_panel_id);
	$sth->execute($panel_id);
	my $s = $sth->fetchall_hashref('bundle_id');
	return [keys %$s];
} 

sub getBundleInfos {
	my ($self, $id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_bundle_infos);
	$sth->execute($id);
	my $h = $sth->fetchall_hashref('bundle_id');
	return $h;
}

sub getBundleGenes {
	my ($self, $id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_bundle_genes);
	$sth->execute($id);
	my $h = $sth->fetchall_hashref('name');
	return [keys %$h];
}

sub getGenesOmimMorbid {
	my ($self) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_genes_omim_morbid);
	$sth->execute();
	my $h = $sth->fetchall_hashref('name');
	return $h;
}

sub getBundleGenesByEnsg {
	my ($self, $id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_bundle_genes_ensg);
	$sth->execute($id);
	my $h = $sth->fetchall_hashref('gid');
	return [keys %$h];
}

###########################
##### METHODS FOR PANEL #####
sub getAllPanelsIds {
	my ($self, $id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_panel_all_ids);
	$sth->execute();
	my $s = $sth->fetchall_hashref('panel_id');
	return [keys %$s];
} 

sub getAllPanelsNames {
	my ($self, $id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_panel_all_names);
	$sth->execute();
	my $s = $sth->fetchall_hashref('name');
	return [keys %$s];
} 

sub getPanelInfos {
	my ($self, $id) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_panel_infos);
	$sth->execute($id);
	my $h = $sth->fetchall_hashref('panel_id');
	return $h;
}


sub getPanelsForGeneName {
	my ($self,$name) = @_;
	my $dbh = $self->getDbh();
	my $sth = $dbh->prepare($self->sql_panel_by_gene_name);
	my $s;
	eval {
		$sth->execute($name) ;
		$s = $sth->fetchall_hashref("pname");
		
	};
	if ($@){
		confess();
	}
	
	return $s;
}

sub getGenesForPanels {
	my ($self,$id) = @_;
	my $dbh = $self->getDbh();
	unless (exists $self->{genes_prepare}){
		$self->{genes_prepare} = $dbh->prepare($self->sql_genes_by_panel_id);
	}
	#warn $self->sql_genes_by_panel_name;
	$self->{genes_prepare}->execute($self,$id);
	my $s = $self->{genes_prepare}->fetchall_hashref("ensg");
	return $s;
}

sub getNbGenesForPanels {
	my ($self,$id) = @_;
	my $dbh = $self->getDbh();
	my $db = $self->database();
	my $sql = $self->sql_count_genes_by_panel_name();
	my $sth = $dbh->prepare($sql);
	$sth->execute($id);
	my $s = $sth->fetchall_hashref("nb");
	my @l = keys %$s;
	return $l[0];
}



sub getAllGenes {
	my ($self,$name) = @_;
	my $dbh = $self->getDbh();
	my $db = $self->database();
	my $sth = $dbh->prepare("SELECT ensg_vid as id FROM $db.ensg");
	$sth->execute();
	my $s = $sth->fetchall_hashref("id");
	return [keys %$s];
	
}
###########################

#######TRASH 


#sub getAllGenes_Panel {
##	sub getAllGenes_Panel {
#	my ($self) = @_;
#	my $dbh = $self->getDbh();
#	my $sth = $dbh->prepare("SELECT gene.name,panel.name 
#		FROM PolyPanelDB.gene as gene,PolyPanelDB.gene_bundle ,PolyPanelDB.panel_bundle,PolyPanelDB.panel  
#		 where gene.gene_id = gene_bundle.gene_id and gene_bundle.bundle_id=panel_bundle.bundle_id and panel_bundle.panel_id = panel.panel_id"
#		 );
#		
#	my $s = $sth->fetchall_hashref();
#}

1;