#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use lib "$RealBin/../../../../../../../GenBo/lib/obj-nodb/";
use Carp;
use GBuffer;



my $h_projects_no_dv;

my $buffer = new GBuffer;
foreach my $pname (@{$buffer->getQuery->listProjectsWithoutDejaVu()}) {
	$h_projects_no_dv->{$pname} = undef;
}



my $dir = $buffer->dejavu_parquet_dir();
$dir =~ s/\/\//\//;
my @files = glob( $dir.'*' );
foreach my $file (@files ) {
	next if -d $file;
	my $file_short = $file;
	$file_short =~ s/$dir//;
	my @ltmp = split('\.', $file_short);
	next if $ltmp[-1] ne 'parquet' and $ltmp[-1] ne 'no_dejavu';
	
	if ($ltmp[0] =~ /^NGS20[0-9][0-9]_[0-9]+$/) {
		if (exists $h_projects_no_dv->{$ltmp[0]}) {
			next if ($ltmp[-1] eq 'no_dejavu');
			del_dejavu($file);
		}
		else {
			next if ($ltmp[-1] eq 'parquet');
			add_dejavu($file);
		}
	}
	else {
		next if ($ltmp[-1] eq 'no_dejavu');
		del_dejavu($file);
	}
}


sub del_dejavu {
	my ($file) = @_;
	my $file_out = $file.'.nodejavu';
	my $cmd = "mv $file $file_out";
	print $cmd."\n";
	`$cmd`;
}

sub add_dejavu {
	my ($file) = @_;
	my $file_out = $file;
	$file_out =~ s/\.no_dejavu//;
	my $cmd = "mv $file $file_out";
	print $cmd."\n";
	`$cmd`;
}