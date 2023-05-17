package GenBoProjectCacheTest;
use strict;
use Data::Dumper;
use Storable qw(store retrieve freeze dclone thaw);
use Moo;
use FindBin qw($Bin);
extends 'GenBoProjectCache';




has infosProject => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $h = {
          'is_somatic' => '0',
          'version' => 'HG19',
          'name' => 'TESTS_F',
          'description' => 'Project TEST',
          'dbname' => 'Polyexome',
          'creation_date' => '2018-09-25 15:00:00',
          'projectTypeId' => '3',
          'projectType' => 'ngs',
          'id' => '1'
        };
		return $h;
	},
);

has pedigree_details => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $hPed = {
          	'GUI' => {
                     'father' => 'GUI_YOS',
                     'children' => {
                                     'GUI_REF' => undef,
                                     'GUI_NAT' => undef
                                   },
                     'mother' => 'GUI_YAE',
                     'all' => {
                                'GUI_REF' => {
                                               'status' => '2',
                                               'sex' => '1'
                                             },
                                'GUI_YOS' => {
                                               'status' => '1',
                                               'sex' => '1'
                                             },
                                'GUI_YAE' => {
                                               'status' => '1',
                                               'sex' => '2'
                                             },
                                'GUI_NAT' => {
                                               'status' => '1',
                                               'sex' => '1'
                                             }
                              }
                   },
          'LEF' => {
                     'father' => 'LEF_JER',
                     'children' => {
                                     'LEF_MEL' => undef
                                   },
                     'mother' => 'BOU_CHR',
                     'all' => {
                                'LEF_MEL' => {
                                               'status' => '2',
                                               'sex' => '2'
                                             },
                                'BOU_CHR' => {
                                               'status' => '1',
                                               'sex' => '2'
                                             },
                                'LEF_JER' => {
                                               'status' => '1',
                                               'sex' => '1'
                                             }
                              }
                   }
		};
		$self->pedigree(1);
		return $hPed;
	},
);

# path du fichier freeze global_infos (ex: avec les infos MD5 de chaque VCF)
has global_infos => (
	is => 'ro',
	lazy => 1,
	default => sub {
		my $self = shift;
		my $file = $self->getCacheBitVectorDir() . '/global_infos.freeze';
		my ($hInfos, $hRes);
		confess() unless -e $file;
		return retrieve $file;
	}
);

# dossier avec le cache vector
sub getCacheBitVectorDir {
	my $self = shift;
	return $self->buffer->dir_project() . '/datas/project_test/vector/';
}

sub getDejaVuLmdbDirTest {
	my $self = shift;
	return $self->buffer->dir_project() . '/datas/project_test/dejavu/';
}

sub get_lmdb_dejavu_tests_f_only {
	my ($self, $mode) = @_;
	unless ($mode) { $mode = 'r'; }
	my $dir_new_dejavu = $self->getDejaVuLmdbDirTest();
	my $no_dejavu_TEST = GenBoNoSqlLmdb->new(dir=>$dir_new_dejavu, mode=>$mode, is_index=>1, name=>'dejavu_test', is_compress=>1); 
	return $no_dejavu_TEST;
}

sub setPatients {
	my $self = shift;
	my (@res, %names);
	push(@res, { 'status' => '2', 'flowcell' => 'A', 'genbo_id' => '0', 'patient_id' => '1', 'capture_id' => '100', 'project_id' => '1', 'name' => 'GUI_REF', 'origin' => 'GUI_REF', 'mother' => 'GUI_YAE', 'sex' => '1', 'description' => '', 'run_id' => '891', 'panel_id' => '1', 'bar_code' => '', 'creation_date' => '2018-09-25 15:03:00', 'father' => 'GUI_YOS', 'project_id_dest' => '0', 'family' => 'GUI' });
	push(@res, { 'status' => '1', 'flowcell' => 'A', 'genbo_id' => '0', 'patient_id' => '2', 'capture_id' => '100', 'project_id' => '1', 'name' => 'GUI_YOS', 'origin' => 'GUI_YOS', 'mother' => '', 'sex' => '1', 'description' => '', 'run_id' => '891', 'panel_id' => '1', 'bar_code' => '', 'creation_date' => '2018-09-25 15:03:00', 'father' => '', 'project_id_dest' => '0', 'family' => 'GUI' });
	push(@res, { 'status' => '2', 'flowcell' => 'A', 'genbo_id' => '0', 'patient_id' => '3', 'capture_id' => '100', 'project_id' => '1', 'name' => 'LEF_MEL', 'origin' => 'LEF_MEL', 'mother' => 'BOU_CHR', 'sex' => '2', 'description' => '', 'run_id' => '891', 'panel_id' => '1', 'bar_code' => '', 'creation_date' => '2018-09-25 15:03:00', 'father' => 'LEF_JER', 'project_id_dest' => '0', 'family' => 'LEF' });
	push(@res, { 'status' => '1', 'flowcell' => 'A', 'genbo_id' => '0', 'patient_id' => '4', 'capture_id' => '100', 'project_id' => '1', 'name' => 'BOU_CHR', 'origin' => 'BOU_CHR', 'mother' => '', 'sex' => '2', 'description' => '', 'run_id' => '891', 'panel_id' => '1', 'bar_code' => '', 'creation_date' => '2018-09-25 15:03:00', 'father' => '', 'project_id_dest' => '0', 'family' => 'LEF' });
	push(@res, { 'status' => '1', 'flowcell' => 'A', 'genbo_id' => '0', 'patient_id' => '5', 'capture_id' => '100', 'project_id' => '1', 'name' => 'LEF_JER', 'origin' => 'LEF_JER', 'mother' => '', 'sex' => '1', 'description' => '', 'run_id' => '891', 'panel_id' => '1', 'bar_code' => '', 'creation_date' => '2018-09-25 15:03:00', 'father' => '', 'project_id_dest' => '0', 'family' => 'LEF' });
	push(@res, { 'status' => '1', 'flowcell' => 'A', 'genbo_id' => '0', 'patient_id' => '6', 'capture_id' => '100', 'project_id' => '1', 'name' => 'GUI_YAE', 'origin' => 'GUI_YAE', 'mother' => '', 'sex' => '2', 'description' => '', 'run_id' => '891', 'panel_id' => '1', 'bar_code' => '', 'creation_date' => '2018-09-25 15:03:00', 'father' => '', 'project_id_dest' => '0', 'family' => 'GUI' });
	push(@res, { 'status' => '1', 'flowcell' => 'A', 'genbo_id' => '0', 'patient_id' => '7', 'capture_id' => '100', 'project_id' => '1', 'name' => 'GUI_NAT', 'origin' => 'GUI_NAT', 'mother' => 'GUI_YAE', 'sex' => '1', 'description' => '', 'run_id' => '891', 'panel_id' => '1', 'bar_code' => '', 'creation_date' => '2018-09-25 15:03:00', 'father' => 'GUI_YOS', 'project_id_dest' => '0', 'family' => 'GUI' });
	foreach my $h (@res){
		$h->{id} = $h->{patient_id};
		$names{$h->{id}} = undef;
		$h->{project} = $self;
		next if exists $self->{objects}->{patients}->{$h->{id}};
		$self->{objects}->{patients}->{$h->{id}} = $self->flushObject('patients', $h);
	}
	return \%names;
}

1;