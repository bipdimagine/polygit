package GenBoNoSqlDejaVuJunctionsCanoniques;
# ABSTRACT: An embedded, NoSQL SQLite database with SQL indexing
use Moo;
use strict;
use warnings;
use Data::Dumper;
extends "GenBoNoSqlDejaVuJunctions";


has extension =>(
	is		=> 'rw',
		default => sub {
		return "dejavu_junctions_canoniques.lite";
	}
);


1;