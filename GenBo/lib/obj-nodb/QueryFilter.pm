=head1 NAME

GenBoFilter : use for save and retreive filter parameter

=head1 SYNOPSIS


=head1 DESCRIPTION

GenBoCacheQuery provides a set of functions to create a cache (stored informations in the table CACHE_ELECTRO in the database)

=head1 METHODS

=cut


package QueryFilter;

use strict;
use Moo;

use Data::Dumper;



has database =>(
	is		=> 'ro',
	required=> 1,
);

has dbh => (
	is		=> 'ro',
	required=> 1,
);


sub get_filter_id {
	my ($self,%arg) = @_;
	my $user_id = $arg{user_id};
	my $filter_name = $arg{filter_name};
	my $project_id = $arg{project_id};
	my $dbh = $self->dbh();
	my $db = $self->database();	
	my $sql = qq{select fu.filter_id as id from $db.filters f,$db.filters_users fu where FILTER_NAME='$filter_name' and PROJECT_ID=$project_id and f.filter_id=fu.filter_id and fu.USER_ID=$user_id};
	return $dbh->selectall_arrayref($sql)->[0][0];
	#return connect::returnOneVal($dbh,$sql);
}

sub get_user_id  {
	my ($self,%arg) = @_;
	my $user_name = $arg{user_name};
	my $sql = qq{SELECT USER_ID as id FROM bipd_users.USER U WHERE U.LOGIN = '$user_name'};
	return $self->dbh->selectall_arrayref($sql)->[0][0];
}
sub newFilter {
	my ($self,%arg) = @_;
	my $user_id = $arg{user_id};
	my $filter_name = $arg{filter_name};
	my $project_id = $arg{project_id};
	
	my $id = $self->get_filter_id(filter_name=>$filter_name, project_id=>$project_id, user_id => $user_id);
	if ($id){
		$self->delete_param(filter_id=>$id);
		return $id;
	}
	my $db = $self->database;
	my $query = qq{
		call $db.new_filter('$filter_name',$project_id,$user_id);
	};
	$self->dbh->do($query);
	return $self->get_filter_id(filter_name=>$filter_name, project_id=>$project_id, user_id => $user_id);
}

sub addParam {
		my ($self,%arg) = @_;
	my $filter_id = $arg{filter_id};
	my $param_name = $arg{param_name};
	my $param_value = $arg{param_value};
	my $db = $self->database;

	 
	my $query = qq{
		insert into $db.filters_param (FILTER_ID,PARAM_NAME,PARAM_VALUE) values ($filter_id,'$param_name','$param_value')
		 ON DUPLICATE KEY UPDATE PARAM_VALUE='$param_value'
		;
	};

	$self->dbh->do($query) ;
}
sub delete_param {
	my ($self,%arg) = @_;
	my $filter_id = $arg{filter_id};
		my $db = $self->database;
		my $query = qq{
			delete from $db.filters_param where FILTER_ID = $filter_id;
		};
		$self->dbh->do($query);
		
}
sub delete_filter {
	my ($self,%arg) = @_;
	my $filter_id = $arg{filter_id};
	my $user_id = $arg{user_id};
	my $db = $self->database;
	my $query2 = qq{
			call $db.delete_filter($filter_id,$user_id);
		};
		
	
		$self->dbh->do($query2);
		return;
}
sub getParam{
	my ($self,%arg) = @_;
	my $filter_id = $arg{filter_id};
	my $db = $self->database;
	my $query = qq{select PARAM_NAME as name, PARAM_VALUE as value from $db.filters_param  where FILTER_ID= $filter_id}; 
	my $sth = $self->dbh->prepare($query);
	$sth->execute();
	return $sth->fetchall_hashref("name");
}
sub getAllFilterName {
	my ($self,%arg) = @_;
	my $project_id = $arg{project_id};
	my $user_id = $arg{user_id};
	my $db = $self->database;
	
	my $query = qq{select f.filter_name as name , f.filter_id as id, creation_date as date from $db.filters f ,$db.filters_users fu where f.filter_id=fu.filter_id and fu.USER_ID=$user_id and PROJECT_ID=$project_id order by creation_date};
	my $sth = $self->dbh->prepare($query);
	$sth->execute();
	return($sth->fetchall_hashref("id"));
	
}

1;