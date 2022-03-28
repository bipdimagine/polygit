package QueryClinvarPathogenic;
use strict;
use Vcf;
use Moose;
use MooseX::Method::Signatures;
use Data::Dumper;
use Config::Std;
use JSON;
use GBuffer;

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

has dbh => (
	is		=> 'ro',
	isa		=> 'DBI::db',
	required=> 1,
	reader	=> 'getDbh',
);

has build => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return 'HG19';
	}
);

has last_release_id => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
        my $dbh =  $self->getDbh;
        my $sql = qq{ SELECT * FROM Clinvar.clinvar_release; };
        my $sth = $dbh->prepare($sql);
        $sth->execute();
        my $res = $sth->fetchall_hashref("release_id");
        my @lIds = sort {$a <=> $b} keys %$res;
        $self->last_release_name($res->{$lIds[-1]}->{'release_name'});
		return $lIds[-1];
	}
);

has last_release_name=> (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
        my $dbh =  $self->getDbh;
        my $sql = qq{ SELECT * FROM Clinvar.clinvar_release; };
        my $sth = $dbh->prepare($sql);
        $sth->execute();
        my $res = $sth->fetchall_hashref("release_id");
        my @lIds = sort {$a <=> $b} keys %$res;
        $self->last_release_id($lIds[-1]);
		return $res->{$lIds[-1]}->{'release_name'};
	}
);

sub getAllVarIds {
	my ($self) = @_;
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT `var_id`, GROUP_CONCAT( `release_id` SEPARATOR ',' ) as releases
					FROM Clinvar.clinvar_pathogenic 
						GROUP BY `var_id`; };
	my $sth = $dbh->prepare($sql);
	$sth->execute();
	my $res = $sth->fetchall_hashref("var_id");
	return keys %$res;
}

sub getAllVarIds_onlyLastRelease {
	my ($self) = @_;
	my $dbh =  $self->getDbh;
	my $last_release_id = $self->last_release_id();
	my $sql = qq{ SELECT `var_id`, GROUP_CONCAT( `release_id` SEPARATOR ',' ) as releases
					FROM Clinvar.clinvar_pathogenic 
						GROUP BY `var_id`
							HAVING releases=?; };
	my $sth = $dbh->prepare($sql);
	$sth->execute($last_release_id);
	my $res = $sth->fetchall_hashref("var_id");
	my @lIds = keys %$res;
	return \@lIds;
}

1;