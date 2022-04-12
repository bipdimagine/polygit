#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../../../../../lib/obj-nodb";
use strict; 
use Bio::SearchIO;
use strict;

use Data::Printer;
use Getopt::Long;
use Carp;
use JSON;

use Data::Dumper;
use Getopt::Long; 
use Storable qw/retrieve freeze thaw nfreeze nstore_fd nstore/;
use GBuffer;


my $release;
GetOptions(
	'release=s'   => \$release,
);
die("please add -release=") unless $release;

my $buffer = new GBuffer;


my $new_dir = $buffer->config->{'public_data_annotation'}->{root} . '/annotations/'.$release."/";
die("\n\n$new_dir doesn't exist... Die...\n\n") unless (-d $new_dir);
$buffer->lmdb_public_dir($new_dir);

my $json_file = $new_dir.'/clinvar/version.json';
die("\n\n$json_file doesn't exist... Die...\n\n") unless (-e $json_file);

open (FILE, $json_file);
my $json_encode = <FILE>;
close (FILE);

my $hVersion = decode_json $json_encode;

my $clinvar_version = $hVersion->{version};
warn "\n### USING CLINVAR Version: $clinvar_version\n";

my $release_db_id;
my $query = $buffer->getQuery();
my $sql_get_id = qq{SELECT * FROM Clinvar.clinvar_release where release_name=?;};
my $dbh = $query->getDbh();
my $sth = $dbh->prepare($sql_get_id);
$sth->execute($clinvar_version);
my $h = $sth->fetchall_hashref('release_name');
if (exists $h->{$clinvar_version}) {
	$release_db_id = $h->{$clinvar_version}->{'release_id'};
}
else {
	die("\n\nERROR: $release_db_id not found in Clinvar DB ! Please add clinvar_id for this version first. Die...\n\n") ;
}
die("\n\nERROR: release_db_id empty ! Die...\n\n") unless (-e $release_db_id);

warn "\n### USING CLINVAR DB ID Found: $release_db_id\n";




my @lSqls;
foreach my $chr_id (1..22, 'X', 'Y', 'MT') {
	warn "\n";
	foreach my $type ('snps', 'insertions', 'deletions') {
		my $db =  $buffer->get_lmdb_database("clinvar", $chr_id, $type);
		my $nb_all = 0;
		my $nb_pathogenic = 0;
		while (my ($k, $v) = $db->next_key_value()) {
			last unless ($k);
			my $pos = $k;
			foreach my $alt_ref (keys %{$v}) {
				$nb_all++;
				my @lSig = split(';', $v->{$alt_ref}->{sig});
				#sig == 5 -> pathogenic
				foreach my $sig (@lSig) {
					if ($sig == 5) {
						my $var_id = $chr_id.'_'.$pos.'_'.$v->{$alt_ref}->{ref_all}.'_'.$alt_ref;
						my $sql = qq{insert into Clinvar.clinvar_pathogenic (var_id, release_id) values ("$var_id", $release_db_id);};
						push(@lSqls, $sql);
						$nb_pathogenic++;
					}
				}
			}
		}
		#warn Dumper $h_var_pathogenic;
		warn "chr$chr_id - $type - Nb pathogenic $nb_pathogenic / $nb_all\n";
		$db->close();
	}
}
#die;
warn "\n\n";
warn '# START INSERT DB:';
warn "\n\n";
foreach my $sql (@lSqls) {
	#warn $sql;
	$dbh->do($sql);
}
	
warn "\n\n";
warn "### FINISHED :)";
warn "\n\n";

