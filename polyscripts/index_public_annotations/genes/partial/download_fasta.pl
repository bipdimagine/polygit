#!/usr/bin/perl
use strict;
use Getopt::Long;
#include the module
use Bio::DB::HTS::Faidx;
 use Data::Dumper;
#create the index object
use FindBin qw($RealBin);
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
use GBuffer;
use Bio::DB::HTS::Tabix;
use FindBin qw($RealBin);

my $version;
GetOptions(
	'version=s'  => \$version,
);
my $file;
my $url;
my $ln;
$ln->{transcript} = "transcripts.fa.gz";
$ln->{protein} = "proteins.fa.gz";
$ln->{pc} = "transcripts.pc.fa.gz";
$url->{HG19} = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_".$version."/GRCh37_mapping/";
$file->{HG19}->{transcript} = "gencode.v".$version."lift37.transcripts.fa.gz";
#$file->{HG19}->{transcript}->{ln} = "tra";
$file->{HG19}->{pc} ="gencode.v".$version."lift37.pc_transcripts.fa.gz";
$file->{HG19}->{protein}= "gencode.v".$version."lift37.pc_translations.fa.gz";
$url->{HG38} = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_".$version."/";

$file->{HG38}->{transcript} = "gencode.v".$version.".transcripts.fa.gz";
$file->{HG38}->{pc} = "gencode.v".$version.".pc_transcripts.fa.gz";
$file->{HG38}->{protein}= "gencode.v".$version.".pc_translations.fa.gz";

warn "DOWNLOAD ==> ";
my @genomes= ("HG19","HG38");
my $root_dir = "/data-isilon/public-data/repository/HG19/annotations/gencode.v".$version."/";
my $fasta_dir = $root_dir."fasta";
system("mkdir $fasta_dir && chmod g+rwx $fasta_dir") unless -e $fasta_dir;
foreach my $genome (@genomes){
	my $dd = $fasta_dir."/".$genome;
	system("mkdir -p $dd && chmod g+rwx $dd") unless -e $dd;
	
	foreach my $type (keys %{$file->{$genome}}){
		warn $type;
		my $f = $file->{$genome}->{$type};
		my $f1 = $dd."/".$file->{$genome}->{$type};
		my $tmp = $f;
		$tmp =~ s/\.gz//;
		my $lna = $ln->{$type};
		system ("cd $dd; wget ".$url->{$genome}.$f." && gunzip $f1 && bgzip $tmp -\@10  && ln -s $f $lna && samtools faidx $lna");
		my $lnf = $dd."/".$lna;
#		warn "ln -s $f1 ".$lnf;
#		system("ln -s $f1 ".$lnf);
#		warn "samtools faidx ".$lnf;
#		system("samtools faidx ".$lnf);
#		warn $RealBin."/change_fai.pl ".$lnf;
		system($RealBin."/change_fai.pl -file=".$lnf);
		die($lnf.".fai") unless -e $lnf.".fai";
		
	}
	
}
my $dd = $fasta_dir."/HG38";
my $file_gff = qq{https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$version/gencode.v$version.annotation.gff3.gz};
system ("cd $dd; wget $file_gff");