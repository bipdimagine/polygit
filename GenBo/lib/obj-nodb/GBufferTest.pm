package GBufferTest;
use strict;
use Moo;
use FindBin qw($Bin);
use GenBoProjectCacheTest;
use Carp;
extends 'GBuffer';
 

has dir_project => (
	is 		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		confess("\n\nERROR: you nedd to define dir_project... Die...\n\n");
	},
);

sub newProjectCache {
	my $self = shift;
	my $args = _checkArguments(@_);
	my $name = 'undef';
	my $release = 'undef';
	my $typeFilters = 'individual';
	my $test = 0;
	if (exists $args->{-test}){ $test = 1; }
	if (exists $args->{-name}){ $name = $args->{-name}; }
	if (exists $args->{-release}){ $release = $args->{-release}; }
	if (exists $args->{-typeFilters}){ $typeFilters = $args->{-typeFilters}; }
	my $project = GenBoProjectCacheTest -> new ( 	name => $name,
												buffer => $self,
												release => $release,
												cache => 1,
												typeFilters => $typeFilters );
	return $project;
}

sub _checkArguments {
	my $index;
    for ($index=0; $index < @_; $index += 2) {
        my $key = $_[$index];
        unless ($key =~ /^\-/o) {
            confess ("Please, could you be so kind as to check your arguments for method 'construct'? I have the impression you wrote '$key' instead of '-$key' -- didn't you?\n");
            return undef;
        }
    }
    my %args = @_;
    # use lowercase in %args's keys -- usefull in case nameSpace was written instead of namespace!
    foreach my $key (keys %args) {
        $args{lc($key)} = $args{$key};
    }
   return \%args;
}

1;