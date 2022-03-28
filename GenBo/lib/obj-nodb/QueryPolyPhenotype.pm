package QueryPolyPhenotype;
use strict;
use Vcf;
use Moose;
use MooseX::Method::Signatures;
use Data::Dumper;
use Config::Std;

has config => (
	is		=> 'ro',
	isa		=> 'HashRef',
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
		return "PolyPhenotypeDB";
	}
);


has all_config => (
	is		=> 'ro',
	#isa		=> 'HashRef',
	required=> 1,
);

has dbh => (
	is		=> 'ro',
	isa		=> 'DBI::db',
	required=> 1,
	reader	=> 'getDbh',
);

###########
#SQL 
###########

##########
#SQL PANEL 
#########

sub getPhenotypesIdFromProjectId {
	my ($self, $project_id) = @_;
	my $dbh = $self->getDbh();
	my $db = $self->database;
	my $sql = qq{SELECT phenotype_id FROM $db.phenotype_project where project_id = ?; }; 
	my $sth = $dbh->prepare($sql);
	$sth->execute($project_id);
	my $s = $sth->fetchall_hashref('phenotype_id');
	return [keys %$s];
} 

sub getPhenotypesIdFromPanelId {
	my ($self, $panel_id) = @_;
	my $dbh = $self->getDbh();
	my $db = $self->database;
	my $sql = qq{SELECT phenotype_id FROM $db.phenotype_panel where panel_id = ?; }; 
	 my $sth = $dbh->prepare($sql);
	$sth->execute($panel_id);
	my $s = $sth->fetchall_hashref('phenotype_id');
	return [keys %$s];
} 

sub getProjectsId {
	my ($self, $id) = @_;
	my $dbh = $self->getDbh();
	my $db = $self->database;
	my $sql = qq{SELECT project_id FROM $db.phenotype_project where phenotype_id = ?; }; 
	 my $sth = $dbh->prepare($sql);
	$sth->execute($id);
	my $s = $sth->fetchall_hashref('project_id');
	return [keys %$s];
} 

sub getPanelsId {
	my ($self, $id) = @_;
	my $dbh = $self->getDbh();
	my $db = $self->database;
	my $sql = qq{SELECT panel_id FROM $db.phenotype_panel where phenotype_id = ?; }; 
	 my $sth = $dbh->prepare($sql);
	$sth->execute($id);
	my $s = $sth->fetchall_hashref('panel_id');
	return [keys %$s];
} 
sub getPhenotypeInfosForProject {
	my ($self, $project_id) = @_;
	my $dbh = $self->getDbh();
	my $db = $self->database();
	my $sql = qq{ SELECT p.* FROM $db.phenotype p, $db.phenotype_project pp where pp.project_id = ? and p.phenotype_id = pp.phenotype_id; };
	my $sth = $dbh->prepare($sql);
	$sth->execute($project_id);
	my $h = $sth->fetchall_hashref('name');
	return $h;
}

sub getPhenotypeInfos {
	my ($self, $id) = @_;
	my $dbh = $self->getDbh();
	my $db = $self->database();
	my $sql = qq{ SELECT * FROM $db.phenotype where phenotype_id=?; };
	my $sth = $dbh->prepare($sql);
	$sth->execute($id);
	my $h = $sth->fetchall_hashref('phenotype_id');
	return $h;
}

sub getPhenotypeInfosFromName {
	my ($self, $name) = @_;
	my $dbh = $self->getDbh();
	my $db = $self->database();
	my $sql = qq{ SELECT * FROM $db.phenotype where name=?; };
	my $sth = $dbh->prepare($sql);
	$sth->execute($name);
	my $h = $sth->fetchall_hashref('name');
	return $h;
}

sub getPhenotypeInfosFromHgmdConcept {
	my ($self, $concept) = @_;
	my $dbh = $self->getDbh();
	my $db = $self->database();
	my $sql = qq{ SELECT * FROM $db.phenotype where concept=?; };
	my $sth = $dbh->prepare($sql);
	$sth->execute($concept);
	my $h = $sth->fetchall_hashref('concept');
	return $h;
}

sub getKeywords {
	my ($self, $id) = @_;
	my $dbh = $self->getDbh();
	my $db = $self->database();
	my $sql = qq{ SELECT k.name as keyword FROM $db.keyword as k,$db.phenotype_keyword as pk where pk.phenotype_id=? and k.keyword_id = pk.keyword_id; };
	my $sth = $dbh->prepare($sql);
	$sth->execute($id);
	my $h = $sth->fetchall_hashref('keyword');
	return [keys %$h];
}

sub getAllPhenotypes {
	my ($self, $id) = @_;
	my $dbh = $self->getDbh();
	my $db = $self->database();
	my $sql = qq{ SELECT phenotype_id as id FROM $db.phenotype; };
	my $sth = $dbh->prepare($sql);
	$sth->execute();
	my $h = $sth->fetchall_hashref('id');
	return [keys %$h];
}

sub getAllPhenotypesName {
	my ($self, $id) = @_;
	my $dbh = $self->getDbh();
	my $db = $self->database();
	my $sql = qq{ SELECT name FROM $db.phenotype; };
	my $sth = $dbh->prepare($sql);
	$sth->execute();
	my $h = $sth->fetchall_hashref('name');
	return [keys %$h];
}

sub getGenesFormPhenotypeName {
	my ($self, $name) = @_;
	my $dbh = $self->getDbh();
	my $db = $self->database();
	my $sql = qq{
		SELECT g.name FROM PolyPhenotypeDB.phenotype as p, PolyPhenotypeDB.phenotype_panel as pp, PolyPanelDB.panel as pa, PolyPanelDB.panel_bundle as pb, PolyPanelDB.gene_bundle as gb, PolyPanelDB.gene as g
			where p.phenotype_id=pp.phenotype_id and pp.current=1
				and pp.panel_id=pa.panel_id 
					and pa.panel_id=pb.panel_id 
						and pb.bundle_id=gb.bundle_id
							and gb.gene_id=g.gene_id
								and p.name=?;
	};
	my $sth = $dbh->prepare($sql);
	$sth->execute($name);
	my $h = $sth->fetchall_hashref('name');
	return [keys %$h];
	
}

1;
