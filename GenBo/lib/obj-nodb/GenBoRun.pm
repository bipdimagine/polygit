package GenBoRun;
use FindBin qw($Bin);

use strict;
use Moo;

use String::ProgressBar;
use Data::Dumper;
use Config::Std;
use File::Path qw(make_path);
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
	
	my $hCaptNames;
	foreach my $capture (@{$self->getProject->getCaptures()}) {
		$hCaptNames->{$capture->name()} = undef;
	}
	
	my @lPat;
	foreach my $h (@$res) {
		 next if  $h->{patient} =~ /GIAB/;
		 next if not exists $hCaptNames->{$h->{capture}};
		 
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
		warn Dumper $res;
		return $res->[0];
	},
);
has sequencing_method => (
	is		=> 'ro',
	#isa		=> 'HashRef[Str]',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $query = $self->getProject->buffer->getQuery();
		return $query->getSequencingMethodFromRunId($self->id);
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

has hash_fastq_screen_patients_files => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $h;
		foreach my $patient (@{$self->getPatients()}) {
			my $patient_name = $patient->name();
			$h->{html}->{$patient_name} = $patient->fastq_screen_file_html();
			$h->{specie}->{$patient_name} = $patient->fastq_screen_file_specie();
		}
		return $h;
	},
);

has hash_fastq_screen_patients_html_files => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $h = $self->hash_fastq_screen_patients_files();
		return if not $h;
		return if not exists $h->{html};
		return $h->{html};
	},
);

has hash_fastq_screen_patients_specie_files => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $h = $self->hash_fastq_screen_patients_files();
		return if not $h;
		return if not exists $h->{specie};
		return $h->{specie};
	},
);

has hash_fastq_screen_patients_specie => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $h = $self->hash_fastq_screen_patients_specie_files();
		return if not $h;
		my $h2;
		foreach my $patient (@{$self->getPatients()}) {
			next unless exists $h->{$patient->name()};
			next unless  $h->{$patient->name()};
			next if not exists $h->{$patient->name()} or not -e $h->{$patient->name()};

			open (F, $h->{$patient->name()});
			my $specie_found = <F>;
			chomp($specie_found);
			close (F);
			$h2->{$patient->name()} = $specie_found;
		}
		return $h2;
	},
);

has fastq_dir => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self        = shift;
		my $machine     = $self->machine;
		my $run_name    = $self->plateform_run_name();
		my $constructor = $self->machine_constructor();
		my $plateform   = $self->plateform();
		my $path        = $self->buffer()->getDataDirectory("sequences-isilon");
	#ici c'ets just epour renviyer su rl'isilon si le projet existe sur l'isilon sinon on va sur le pure '
		if ($path){
		my $seq_dir2 =
			$path . "/"
		  . $constructor . "/"
		  . $machine . "/"
		  . $plateform
		  . "/$run_name.saved/";
		return $seq_dir2 if -e $seq_dir2;
			my $seq_dir =
			$path . "/"
		  . $constructor . "/"
		  . $machine . "/"
		  . $plateform
		  . "/$run_name/";
		return $seq_dir if -e $seq_dir;
		}
		$path        = $self->buffer()->getDataDirectory("sequences");
		my $seq_dir2 =
			$path . "/"
		  . $constructor . "/"
		  . $machine . "/"
		  . $plateform
		  . "/$run_name.saved/";
		return $seq_dir2 if -e $seq_dir2;
		my $seq_dir =
			$path . "/"
		  . $constructor . "/"
		  . $machine . "/"
		  . $plateform
		  . "/$run_name/";
		return $seq_dir;
		

	},
);




has bcl_dir => (
	is      => 'ro',
	lazy    => 1,
	default => sub {
		my $self        = shift;
		my $machine     = $self->machine;
		my $run_name    = $self->run_name();
		my $constructor = $self->machine_constructor();
		my $plateform   = $self->plateform();
		my $path        = $self->buffer()->getDataDirectory("bcl");

		my $seq_dir =
			$path . "/"
		  . $constructor . "/"
		  . $machine . "/"
		  . $plateform
		  . "/$run_name/";
		return $seq_dir;

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
has run_name => (
	is		=> 'ro',

	lazy 	=> 1,
	default => sub {
		my $self = shift;
		return $self->infosRun->{run_name};
	},
);
has sample_sheet => (
	is		=> 'ro',

	lazy 	=> 1,
	default => sub {
		my $self = shift;
		return connect::uncompressData($self->infosRun->{document});
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