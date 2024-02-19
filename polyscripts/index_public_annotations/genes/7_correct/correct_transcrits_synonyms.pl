#!/usr/bin/perl
use FindBin qw($Bin);
use lib $Bin;
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use strict;
use Getopt::Long; 
use Data::Dumper;
use GenBoNoSql;
use GenBoNoSqlAnnotation;
use Carp;


require "$Bin/../packages/parse_gff.pm";

my $version;
my $genome_version = "HG19";
GetOptions(
	'version=s' => \$version,
);


#my $lite_synonyms_dir = "/data-isilon/bipd-src/mbras/";
#my $lite_synonyms = "/data-isilon/bipd-src/mbras/synonyms.annot.search";
#my $file = "/data-isilon/public-data/repository/HG19/annotations/gencode.v43/tabix/gencode.v43lift37.metadata.RefSeq.gz";

my $lite_synonyms_dir = "/data-isilon/public-data/repository/HG19/annotations/gencode.v".$version."/lmdb/";
my $file = "/data-isilon/public-data/repository/HG19/annotations/gencode.v".$version."/tabix/gencode.v".$version."lift37.metadata.RefSeq.gz";


my $lite_synonyms_file = $lite_synonyms_dir.'/synonyms.annot.search';
$lite_synonyms_file =~ s/\/\//\//g;

confess("\n\nERROR: $file doesn't exists. die...\n\n") if (not -e $file);
confess("\n\nNo backup found of original lmbd. Use before: \ncp $lite_synonyms_file $lite_synonyms_file.origin \n\n") if (not -e $lite_synonyms_file.'.origin');


my $no = GenBoNoSqlAnnotation->new(
	dir  => $lite_synonyms_dir,
	mode => "w"
);

my $h_refseq;
open (READ, "zcat $file | ") || die "Cannot open $file: $!.\n";
while (my $line = <READ>){
	chomp($line);
	my ($enst_1, $nm_1, @t) = split(" ", $line);
	my ($enst, $tmp) = split('\.', $enst_1);
	my ($nm, $tmp2) = split('\.', $nm_1);
	$h_refseq->{$enst}->{$nm} = undef;
}
close (READ);

my $nb_updated = 0;
foreach my $enst_id (keys %{$h_refseq}) {
	next if scalar keys %{$h_refseq->{$enst_id}} == 1;
	my ($key, $value) = $no->get_key_value( "synonyms", $enst_id);
	my $new_key = join(' ',sort keys %{$h_refseq->{$enst_id}});
	my @lold = split(' ', $key);
	foreach my $old_v (@lold) {
		next if exists $h_refseq->{$enst_id}->{$old_v};
		$new_key .= ' '.$old_v;
	}
	
	if ($key and $new_key and $key ne $new_key) {
		print "\n";
		print "delete\n";
		print "$key -> $value\n";
		$no->delete_bulk("synonyms",$key);
		
		print "add\n";
		print "$new_key -> $value\n";
		$no->put("synonyms",$new_key,$value);
		$nb_updated++;
	}
}
print "\nNb transcripts updated: $nb_updated\n\n";

#$no->delete_bulk("synonyms",'NR_165197 ENST00000594239 ENST00000594239_X  genboid:ENST00000594239_X transcript chr:X CCDS14757  ENST00000594239_X');
#$no->put("synonyms",'NM_001099857.5 NR_165197 ENST00000594239 ENST00000594239_X  genboid:ENST00000594239_X transcript chr:X CCDS14757  ENST00000594239_X','ENST00000594239_X');

