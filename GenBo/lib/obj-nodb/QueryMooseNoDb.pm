package QueryMooseNoDb;

use strict;
use Vcf;
use Moose;
use MooseX::Method::Signatures;
use Data::Dumper;
use Config::Std;
extends "QueryMoose";

has projects => (
	is		=> 'ro',
	isa		=> 'Config::Std::Hash',
	lazy =>1,
	default	=> sub {
		my $self = shift; 
		my $config;
		my $file_list_project =  $self->all_config()->{project_data}->{root}."/".$self->all_config()->{nodb}->{projects};
		read_config $file_list_project => $config ;
		return $config;
	},
);



has project_id =>(
	is		=> 'rw',
	isa		=> 'Str',
);

has config_patients => (
	is		=> 'ro',
	isa		=> 'Config::Std::Hash',
	lazy =>1,
	default	=> sub {
	my $self = shift; 
	confess() unless $self->project_id;
	my $file_patients = $self->all_config()->{project_data}->{root}."/ngs/".$self->project_id."/sample.list";
	my $config;
	read_config $file_patients => $config ;
	return $config;
	}
);

sub getProjectByName {
	my ($self, $name, $verbose) = @_;
	confess("no poject with $name defined") unless exists $self->projects()->{$name};

	return {
          'version' => $self->projects()->{$name}->{version},
          'projectType' => 'ngs',
          'projectTypeId' => '3',
          'name' => $name,
          'id' => $name,
          'dbname' => 'nodb'
        };
}

sub getPatients {
	my ($self, $project_id) = @_;
	$self->project_id($project_id);
	my @array;
	foreach my $c (keys %{$self->config_patients}){
		my $p = {
			'origin' => $project_id,
                  'name' => $c,
                      'genbo_id' => $c,
                      'patient_id' => $c,
                      'capture_id' => $self->config_patients()->{$c}->{capture},
                  'run_id' => '1',
                  'project_id' => $project_id,
                  
		};
		push(@array,$p);
		
	}

	return \@array;
	
}


sub getSequencingMachines {
	my ($self, $project_id) = @_;
	$self->project_id($project_id);
	my @res = split(",",$self->projects->{$project_id}->{machine});
	return \@res;
	
}
sub getOriginMethods {
	my ($self, $project_id, $type) = @_;
	$self->project_id($project_id);
	$type = lc($type);
	warn $type;
	return [split(",",$self->projects->{$project_id}->{$type})];
}
sub getProjectListForUser {
	my ($self,$login,$pwd)=@_;
	my @res;
	foreach my $p (keys %{$self->projects}){
		my $hp;
		$hp->{type} = "ngs";
		$hp->{name} = $p;
		$hp->{genome} = $self->projects->{$p}->{version};
		$hp->{description} = $self->projects->{$p}->{description};
		push (@res,$hp);
	}
	return \@res; 
}
1;