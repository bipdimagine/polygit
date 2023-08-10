package GenBoNoSqlDejaVuJunctionsPhenotype;
# ABSTRACT: An embedded, NoSQL SQLite database with SQL indexing
use Moo;
use strict;
use warnings;
use Data::Dumper;
extends "GenBoNoSqlDejaVuJunctions";


has phenotype_name =>(
	is		=> 'ro'
);

has extension =>(
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return "dejavu_junctions_".$self->phenotype_name().".lite";
	}
);


1;