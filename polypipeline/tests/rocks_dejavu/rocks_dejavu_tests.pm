package rocks_dejavu_tests;
use strict;
use Moo;
use FindBin qw($Bin);
use Data::Dumper;
use RocksDB;
use Getopt::Long;
use Carp;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use GBuffer;
use GenBoNoSqlRocksGenome;
use Term::ANSIColor;
use MCE::Loop;
use MCE::Flow;

use Test::Simple tests => 27;


has confess_as_soon_as_possible => (
	is      => 'rw',
	lazy    => 1,
	default => 1,
);

has verbose => (
	is      => 'rw',
	lazy    => 1,
	default => 0,
);

has use_chromosomes => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my @list = (1..22, 'X', 'Y', 'MT');
		return \@list;
	},
);

has path_origin => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $buffer = new GBuffer;
		confess("\n\nERROR: release version not defined. Die\n\n") if not $self->release();
		return $buffer->config_path("root","dejavu").'/'.$self->release()."/variations/rocks/";
	},
);

has path_new => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		confess("\n\nERROR: NEW PATH TO COMPARE version not defined. Die\n\n")
	},
);

has release => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return 'HG19' if lc($self->{path_new}) =~ /hg19/;
		return 'HG38' if lc($self->{path_new}) =~ /hg38/;
		confess("\n\nERROR: release version not defined. Die\n\n");
	},
);

has limit_ratio_lost => (
	is      => 'rw',
	lazy    => 1,
	default => 2,
);

has limit_ratio_new => (
	is      => 'rw',
	lazy    => 1,
	default => 2,
);

has resume_tests => (
	is      => 'rw',
	lazy    => 1,
	default => sub { {} },
);

has print_resume_text => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		print "\n";
		print "# RESUME PATHS / RATIOS\n";
		print "-> old path: ".$self->path_origin()."\n";
		print "-> new path: ".$self->path_new()."\n";
		print "-> cutoff ratio (lost variants): ".$self->limit_ratio_lost()."%\n";
		print "-> cutoff ratio (new variants): ".$self->limit_ratio_new()."%\n";
		print "\n";
	},
);

sub set_directory_path_origin {
	my ($self, $path) = @_;
	confess("\n\nERROR: $path doesn't exists. Die\n\n") if not -d $path;
	$self->{path_origin} = $path;
}

sub set_directory_path_new {
	my ($self, $path) = @_;
	confess("\n\nERROR: $path doesn't exists. Die\n\n") if not -d $path;
	$self->{path_new} = $path;
}

sub global_check {
	my ($self) = @_;
	$self->confess_as_soon_as_possible(0);
}

sub launch_test_all_chromosomes {
	my ($self) = @_;
	confess("\n\nERROR: ".$self->path_origin()." doesn't exists. Die\n\n") if not -d $self->path_origin();
	confess("\n\nERROR: ".$self->path_new()." doesn't exists. Die\n\n") if not -d $self->path_new();
	$self->print_resume_text();

	my $list_chr = $self->use_chromosomes();
	my (@files_ok, $h_part_files);
	foreach my $chr_name (@$list_chr) {
		opendir my $dir, $self->{path_new}.'/'.$chr_name.'.g.rocks/' or die "Cannot open directory: $!";
		my @files = readdir $dir;
		closedir $dir;
		foreach my $part_path (@files) {
			next if $part_path eq '.';
			next if $part_path eq '..';
			next if $part_path eq 'chunk.json';
			next if not $part_path =~ /\.[0-9]+\.[0-9]+.+/;
			push(@files_ok, $part_path);
			$h_part_files->{$part_path} = $chr_name;
		}
	}
	
	my $h_res_by_chr;
	my $nb_total_path1 = 0;
	my $nb_total_path2 = 0;
	my $nb_total_lost = 0;
	my $nb_total_new = 0;
	MCE::Loop->init(
		max_workers => 'auto', chunk_size => 1,
		gather => sub {
	        my ($mce, $data) = @_;
	        my $chr_name = $data->{chr};
	        my $part_name = $data->{part_name};
			$nb_total_path1 += $data->{nb1};
			$nb_total_path2 += $data->{nb2};
			$nb_total_lost += $data->{nb_lost};
			$nb_total_new += $data->{nb_new};
			$h_res_by_chr->{$chr_name}->{nb1} += $data->{nb1};
			$h_res_by_chr->{$chr_name}->{nb2} += $data->{nb2};
			$h_res_by_chr->{$chr_name}->{nb_lost} += $data->{nb_lost};
			$h_res_by_chr->{$chr_name}->{nb_new} += $data->{nb_new};
			if ($data->{in_error}) {
				$h_res_by_chr->{$chr_name}->{errors}->{$part_name} = $data;
				$h_res_by_chr->{$chr_name}->{errors}->{$part_name}->{ratio_lost} = $self->get_ratio_lost_variants($data->{nb1}, $data->{nb2}, $data->{nb_lost});
				$h_res_by_chr->{$chr_name}->{errors}->{$part_name}->{ratio_new_var} = $self->get_ratio_new_variants($data->{nb1}, $data->{nb2}, $data->{nb_new});
				confess ("\n\nERROR for $part_name\n\n") if $self->confess_as_soon_as_possible();
			}
	    }
	);
	
	mce_loop {
		my ($mce, $chunk_ref, $chunk_id) = @_;
		if (ref($chunk_ref) ne "ARRAY") {
			my $chr_name = $h_part_files->{$chunk_ref};
			my $hash = $self->check_this_part($chr_name, $chunk_ref);
			MCE->gather($chunk_id, $hash);
		}
		else {
			foreach my $region (@$chunk_ref){
				my $chr_name = $h_part_files->{$region};
				my $hash = $self->check_this_part($chr_name, $region);
				MCE->gather($chunk_id, $hash);
			}
		}
	} sort @files_ok;
	
	my $ratio_lost = $self->get_ratio_lost_variants($nb_total_path1, $nb_total_path2, $nb_total_lost);
	my $ratio_new_var = $self->get_ratio_new_variants($nb_total_path1, $nb_total_path2, $nb_total_new);
	if ($self->verbose()) {
		print "\n\n";
		print "# TOTAL var in ref DV: ".$nb_total_path1."\n";
		print "# TOTAL var in new DV: ".$nb_total_path2."\n";
		print "# LOST var (in ref DV and not in new DV): ".$nb_total_lost.' ('.$ratio_lost.'%)'."\n";
		print "# NEW var (in new DV and not in ref DV):: ".$nb_total_new.' ('.$ratio_new_var.'%)'."\n";
	}
	
	my $is_all_ok = 0;
	if ($ratio_lost > $self->limit_ratio_lost()) {
		print "\n-> ERROR !\n" if $self->verbose();
	}
	elsif ($ratio_new_var > $self->limit_ratio_new()) {
		print "\n-> ERROR !\n" if $self->verbose();
	}
	else {
		print "\n-> OK!\n" if $self->verbose();
		$is_all_ok = 1;
	}
	print "\n";
	
	my $nb_errors = 0;
	foreach my $chr_name (@$list_chr) {
		my $nb_1 = $h_res_by_chr->{$chr_name}->{nb1};
		my $nb_2 = $h_res_by_chr->{$chr_name}->{nb2};
		my $nb_lost = $h_res_by_chr->{$chr_name}->{nb_lost};
		my $nb_new = $h_res_by_chr->{$chr_name}->{nb_new};
		my $ratio_lost = $self->get_ratio_lost_variants($nb_1, $nb_2, $nb_lost);
		my $ratio_new_var = $self->get_ratio_new_variants($nb_1, $nb_2, $nb_new);
		my $status_errors = 'OK';
		if (exists $h_res_by_chr->{$chr_name}->{errors}) {
			$status_errors = 'ERROR';
			$nb_errors += scalar(keys %{$h_res_by_chr->{$chr_name}->{errors}});
		}
		my $text = "chr$chr_name - nb_old:$nb_1 - nb_new:$nb_2 - LOST:$nb_lost ($ratio_lost%) - NEW:$nb_new ($ratio_new_var%)";
		ok( $status_errors eq 'OK', $text) or do {
			print "\n";
			print colored(['white on_red'], 'ERROR for chromosome '.$chr_name);
			print "\n";
			foreach my $part_name (keys %{$h_res_by_chr->{$chr_name}->{errors}}) {
				my $this_log = '  -> part '.$part_name;
				$this_log .= ' - LOST '.$h_res_by_chr->{$chr_name}->{errors}->{$part_name}->{ratio_lost}.'% (- '.$h_res_by_chr->{$chr_name}->{errors}->{$part_name}->{nb_lost}.' var)' if ($h_res_by_chr->{$chr_name}->{errors}->{$part_name}->{ratio_lost} > $self->limit_ratio_lost());
				$this_log .= ' - GAIN '.$h_res_by_chr->{$chr_name}->{errors}->{$part_name}->{ratio_new_var}.'% (+ '.$h_res_by_chr->{$chr_name}->{errors}->{$part_name}->{nb_new}.' var)' if ($h_res_by_chr->{$chr_name}->{errors}->{$part_name}->{ratio_new_var} > $self->limit_ratio_new());
				print colored(['white on_red'], $this_log);
				print "\n";
			}
			print "\n";
			confess("\n\nERROR for chr$chr_name\n\n") if $self->confess_as_soon_as_possible();
		};
	}
	ok( $nb_errors == 0, "CHECK parts chromosomes with problems") or do {
		print "\n";
		print colored(['white on_red'], "ERROR GLOBAL RESULTS -> founds $nb_errors part(s) chromosome(s) with problem(s) !!");
		print "\n\n";
		confess("\n\nERROR founds $nb_errors part(s) chromosome(s) with problem(s) !! \n\n");
		die;
	};
	ok( $is_all_ok == 1, "GLOBAL DEJAVU - nb_old:$nb_total_path1 - nb_new:$nb_total_path2 - LOST:$nb_total_lost ($ratio_lost%) - NEW:$nb_total_new ($ratio_new_var%)") or do {
		confess("\n\nERROR founds $nb_errors part(s) chromosome(s) with problem(s) !! \n\n");
	};
		
	return 1;
}

sub check_this_part {
	my ($self, $chr_name, $part_path) = @_;
	my ($dir_to_check_A, $dir_to_check_B, $nb, $nb2, $nb_lost, $nb_new, $nb_pb_long_reads);
	$dir_to_check_A = $self->path_origin()."/$chr_name.g.rocks/$part_path/";
	$dir_to_check_B = $self->path_new()."/$chr_name.g.rocks/$part_path/";
	($nb, $nb2, $nb_lost, $nb_new) = check_this_directory($part_path, $dir_to_check_A, $dir_to_check_B);
	my $ratio_lost = $self->get_ratio_lost_variants($nb, $nb2, $nb_lost);
	my $ratio_new_var = $self->get_ratio_new_variants($nb, $nb2, $nb_new);
	my $in_error = 0;
	my $log = 'part:'.$part_path.' found in OLD_OK:'.$nb." - NEW_FROM_PURE_PARQUET: ".$nb2." | OLD not_in NEW: ".$nb_lost.' ['.sprintf("%.2f", $ratio_lost).'%]'.' - NEW not_in OLD: '.$nb_new.' ['.sprintf("%.2f", $ratio_new_var).'%]';
	if ($nb_lost > 10 and $ratio_lost > $self->limit_ratio_lost()) {
		$in_error = 1;
	}
	elsif ($nb_new > 10 and $ratio_new_var > $self->limit_ratio_new()) {
		$in_error = 1;
	}
	print $log."\n" if $self->verbose();
	
	ok($in_error == 0, $log) or do {
		print "\n";
		print colored(['red on_bright_yellow'], 'ERROR '.$log);
		print "\n\n";
	};
	
	my $data;
	$data->{'chr'} = $chr_name;
	$data->{'part_name'} = $part_path;
	$data->{'nb1'} = $nb;
	$data->{'nb2'} = $nb2;
	$data->{'nb_lost'} = $nb_lost;
	$data->{'nb_new'} = $nb_new;
	$data->{'in_error'} = $in_error;
	return $data;
}

sub get_ratio_lost_variants {
	my ($self, $nb1, $nb2, $nb_lost) = @_;
	return 0 if $nb1 == 0;
	return 0 if $nb_lost == 0;
	my $ratio_lost = ($nb_lost/$nb1 *100);
	return sprintf("%.2f",$ratio_lost);
}

sub get_ratio_new_variants {
	my ($self, $nb1, $nb2, $nb_new) = @_;
	return 0 if $nb2 == 0;
	return 0 if $nb_new == 0;
	my $ratio_new = ($nb_new/$nb2 *100);
	return sprintf("%.2f",$ratio_new);
}

sub check_this_directory {
	my ($part_path, $dir_to_check_A, $dir_to_check_B) = @_;
	my $hashA = {};
	my $rocks = RocksDB->new($dir_to_check_A,{ read_only => 1 });
	my $iter = $rocks->new_iterator->seek_to_first;
	my $nb =0;
	while (my ($key, $value) = $iter->each) {
	    $nb ++;
	    $hashA->{$key} ++;
	}
	my $hashB = {};
	my $rocks2 = RocksDB->new($dir_to_check_B,{ read_only => 1 });
	my $iter2 = $rocks2->new_iterator->seek_to_first;
	my $nb2 =0;
	my $notA = 0;
	my $nb3 ; 
	while (my ($key, $value) = $iter2->each) {
	    $nb2 ++;
	    unless (exists $hashA->{$key}) {
			$notA++;
	    }
	    $hashB->{$key} ++;
	}
	my $notB = 0;
	my $iter3 = $rocks->new_iterator->seek_to_first;
	while (my ($key, $value) = $iter3->each) {
	    unless (exists $hashB->{$key}) {
			$notB++;
	    }
	}
	return ($nb, $nb2, $notB, $notA);
}

1;