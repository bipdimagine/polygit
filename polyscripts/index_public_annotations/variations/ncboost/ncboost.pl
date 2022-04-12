#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../../../../../lib/obj-nodb";
use strict; 
use strict;
use Getopt::Long;
use Carp;
use GenBoNoSqlLmdbScore;
use Data::Dumper;
use Getopt::Long; 
use Storable qw/retrieve freeze thaw nfreeze nstore_fd nstore/;
use GBuffer;
use Sys::Hostname;
use Parallel::ForkManager;
use Vcf;
use JSON::XS;
use Digest::MD5::File qw( file_md5_hex );
use Date::Tiny;
 use POSIX;
 use Set::IntSpan::Fast::XS;
 use GenBoBitVector;
#use GenBoNoSql;
#require "$Bin/dbsnp.pm";
 require "$Bin/../packages/save.pm";
my $allsnps;


#my $chromosomes = [2];
#my $dir = "/data-xfs/public-data/HG19/snp/exac/latest/";
my $pm = new Parallel::ForkManager(8);
my $buffer = GBuffer->new();
my $description = {
	score_description=>["a","b"],
	pack_string=>"w2",
	factor=>[100,100],
};



my $version;
GetOptions(
	'version=s'   => \$version,
);
die("please add -version=") unless $version;
my $prg = "ncboost";
my $dir_vcf  = "/public-data/repository/HG19/$prg/$version/tabix/";
my @files = `ls $dir_vcf/*.gz`;
chomp(@files);
die() unless @files;
my $file = $files[0];
die() unless -e $file;
my @chr_list = `tabix -l $file`;
chomp (@chr_list);
warn Dumper @chr_list;
my $chromosomes = \@chr_list;
  
  warn "work on : $file ";
  my $md5_1 ;#= file_md5_hex( $file);
  warn $file." ".$md5_1;

die("clinar file doesn t exists " .$file ) unless -e $file;
my $dir_out  = "/public-data/repository/HG19/$prg/$version/";

if (-e $dir_out."/lmdb"){
	die("hey the output directory already exists !!! $dir_out/lmdb");
}
system("mkdir  $dir_out/lmdb && chmod a+rwx $dir_out/lmdb");


open(JSON,">$dir_out/lmdb/description.json") or die();
print JSON encode_json $description;
close JSON;

my $hrun;
my $id = time;
my $nb;

$pm->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    	warn "end".$data->{id}; 
		delete $hrun->{$data->{id}};		    
    	
    }
  );

  foreach my $chr (@$chromosomes){
  	#next unless $chr eq 21;
	$id ;
	$hrun->{$id."_".$chr}++;
	my $pid = $pm->start and next;
	warn "start $chr";
	run_chr($chr);
	$pm->finish(0,{id=>$id."_".$chr});
}#end chr
  $pm->wait_all_children;


if (keys %$hrun){
	warn Dumper $hrun;
	warn "ERROR ON YOUR LMDB DATABASE $prg";
	die();
}
my $hversion;
$hversion->{name} = "$prg";
$hversion->{version} = "$version";
$hversion->{file} = $file;
  my $d = Date::Tiny->now;

$hversion->{date} =  $d->as_string;
$hversion->{md5sum} = $md5_1;
open(JSON,">$dir_out/lmdb/version.json") or die();
print JSON encode_json $hversion;
close JSON;


sub run_chr {
	my ($chr) = @_;

	my $span = Set::IntSpan::Fast::XS->new();

	my $lmdb = GenBoNoSqlLmdbScore->new(dir=>$dir_out."/lmdb",mode=>"c",name=>$chr,is_integer=>1,is_compress=>0);
	#score_description=>$description->{score},pack_string=>$description->{pack_string},factor=>$description->{factor});
	#my $lmdb = GenBoNoSqlLmdbChromosomeBitVector->new(dir=>$dir_out."/lmdb",mode=>"c",name=>$chr,is_integer=>1,is_compress=>0);#score_description=>$description->{score},pack_string=>$description->{pack_string},factor=>$description->{factor});
	
	open (TABIX,"tabix $file $chr | ");
	my $nb;
	while (<TABIX>){
		chomp();
		$nb++;
		my @value = split(" ",$_);
		my $v1 = ceil($value[4]*100);
		my $v2 = ceil($value[5]*100);
		die() if $value[4] > 1;
		die() if $value[5] > 1;
		$lmdb->put_score($value[1],{a=>$value[4],b=>$value[5]});
		my $h =  $lmdb->get_score($value[1]);
		die() if $value[0] ne $chr;
		warn $chr." ".$nb if $nb%1000000 == 0;
		
		#print $value[1]." ".$value[4]." ".$value[5]."\n";
		
	}
	$lmdb->close;
	
	
}