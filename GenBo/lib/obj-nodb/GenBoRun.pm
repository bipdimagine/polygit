package GenBoRun;
use FindBin qw($Bin);

use strict;
use Moose;
use MooseX::Method::Signatures;
use String::ProgressBar;
use Data::Dumper;
use Config::Std;
extends 'GenBo';

sub setPatients {
	my $self = shift;
	my $pids; 
	map{$pids->{$_->id} ++ } @{$self->getProject()->getPatientsAndControl()};
	
	my $query = $self->buffer->getQuery();
	my $res = $query->getPatientsFromRunId($self->getProject->id,$self->id);
	
	my %hids;
	foreach my $h (@$res){
		next unless exists ($pids->{$h->{id}});
		$hids{$h->{id}} = undef;
	}
	return \%hids;
}


sub getFamilies {
		my $self = shift;
		my @res;
		foreach my $fam (sort{$a->name cmp $b->name }@{$self->project->getFamilies}){
				my ($nb) = grep {$_->getRun->id eq $self->id} @{$fam->getMembers};  
				push(@res,$fam) if $nb;
			
		}
		return \@res;
	
}
sub getPatientsControl {
	my $self = shift;
	my @z = grep{$_->is_control && $_->getRun->id eq $self->id } @{$self->getProject()->getPatientsAndControl} ;
	return \@z;
}

has lmdb_cache_dir => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
	my $self = shift;
	my $dir_root = $self->project->getCoverageDir();
	my $dir_out = $dir_root."/".$self->id();
	unless (-e  $dir_out){
		system ("mkdir  $dir_out");
		system("chmod a+rwx $dir_out");
	}
	return $dir_out;
	}
);

sub get_cai_count {
	my ($self,$mode) = @_;
	return $self->_get_lmdb($mode,"cai_count");
}
sub _get_lmdb {
	my($self,$mode,$fname) = @_;
	die() unless $fname;
	$mode = "r" unless $mode;
	my $dir_out = $self->lmdb_cache_dir();
	if ($mode eq "r" && -z "$dir_out/$fname"){
		return undef;
	}
	
	my $no2 = GenBoNoSqlLmdb->new(dir=>$dir_out,mode=>$mode,name=>$fname,is_compress=>0,vmtouch=>$self->buffer->cache);
	$no2->clean_files() if ($mode eq "c");
	$no2->create() if ($mode eq "c");
	return $no2;
}


sub getAllPatientsInfos {
	my $self = shift;
	return $self->{patients_infos} if exists $self->{patients_infos};
	
	my $pids; 
	map{$pids->{$_->id} ++ } @{$self->getProject()->getPatients()};
	my $query = $self->buffer->getQuery();
	my $res = $query->getAllPatientsFromRunId($self->id);
	
	my @lPat;
	foreach my $h (@$res) {
		 next if  $h->{patient} =~ /GIAB/;
		#if (exists $pids->{$h->{id}}) {
			push(@lPat, $h);
		#}
	}
	$self->{patients_infos} = \@lPat;
	
	return $self->{patients_infos};
}

has infosRun => (
	is		=> 'ro',
	#isa		=> 'HashRef[Str]',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $query = $self->getProject->buffer->getQuery();
		my $res = $query->getRunInfosById( $self->id );
		return $res->[0];
	},
);

has description => (
	is		=> 'ro',

	lazy 	=> 1,
	default => sub {
		my $self = shift;
	
		return $self->infosRun->{description};
	},
);

has plateform_run_name => (
	is		=> 'ro',

	lazy 	=> 1,
	default => sub {
		my $self = shift;
		return $self->infosRun->{plateform_run_name};
	},
);


has date => (
	is		=> 'ro',

	lazy 	=> 1,
	default => sub {
		my $self = shift;
	
		return $self->infosRun->{date};
	},
);

has machine => (
	is		=> 'ro',

	lazy 	=> 1,
	default => sub {
		my $self = shift;
	
		return $self->infosRun->{machine};
	},
);
has machine_constructor => (
	is		=> 'ro',

	lazy 	=> 1,
	default => sub {
		my $self = shift;
	
		return $self->infosRun->{constructor};
	},
);
has plateform => (
	is		=> 'ro',

	lazy 	=> 1,
	default => sub {
		my $self = shift;
	
		return $self->infosRun->{plateform};
	},
);

sub setCaptures {
	my $self = shift;
	my %captIds;
	foreach my $patient (@{$self->getPatients()}) {
		foreach my $thisCapture (@{$patient->getCaptures()}) {
			$captIds{$thisCapture->id()} = undef;
		}
	}
	return \%captIds;
}

1;