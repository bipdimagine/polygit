#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
#use Set::IntSpan;
use GenBoNoSqlLmdb;

use Carp;
use strict;
use Set::IntSpan::Fast::XS ;
use Data::Dumper;
use GBuffer;
use Proc::Simple;
use Getopt::Long;
use File::Util;
my $myproc = Proc::Simple->new(); 

$myproc->redirect_output ("/dev/null", "/dev/null");
my $fork;
my $project_name;
my $debug;
my $type;
my $set;

GetOptions(
	'project=s' => \$project_name,
	'debug=s' => \$debug,
	'name=s' => \$type,
	'set=s' => \$set,
);

my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $project_name );


my $dir_in = "/data-isilon/download/defidiag/www.cnrgh.fr/data/".$type;
$dir_in = "/data-isilon/download/devodecode/www.cnrgh.fr/data/".$type if lc($type) =~/devo/;
$dir_in = "/data-isilon/download/malar/www.cnrgh.fr/data/".$type if lc($type) =~/malrar/;
$dir_in = "/data-isilon/download/alpsseq/www.cnrgh.fr/data/".$type if lc($type) =~/alps/;
my @notfound;
my %dejavu;
my $mv;
foreach my $patient (@{$project->getPatients}){
	my $bc = $patient->barcode();
	my @files = `find $dir_in -name *$bc*`;
	chomp(@files);
	my ($bam) = grep{$_=~ /.bam$/} @files;
#	my ($bai) = grep{$_=~ /.bai$/} @files;
	my $bai = $patient->getBamFileName().".bai"; 
	next if -e $bai;
	next  if !(-e $bam) && -e $bai;
	push(@notfound,$bc) unless $bam ;
	warn "$bc " unless $bam ;
}
#warn Dumper @notfound;
#die();

my $nb;
foreach my $bc (@notfound){
	warn "download $bc => ".$nb."/".scalar(@notfound)."/".scalar(@{$project->getPatients});
	$nb ++;
	system("./download_one_file.pl -name=$type -set=$set -bc=$bc");
}

warn "CHECK DOWNLOAD !!! ";
system("./check_download.pl -name=$type -set=$set");

foreach my $patient (@{$project->getPatients}){
	my $bc = $patient->barcode();
	
	#warn $bc;
	my @files = `find $dir_in -name *$bc*`;
	chomp(@files);
	my ($bam) = grep{$_=~ /.bam$/} @files;
	my ($bai) = grep{$_=~ /.bai$/} @files;
	
	#my $md5out = `md5sum $dir_in/$bamout | cut -f 1 -d " "`;
	#chomp($md5out);
	#warn $md5out;
	
	#die("$md5out- $md5in-") if $md5out ne $md5in;
	my $bai = $patient->getBamFileName().".bai"; 
	next if -e $bai;
	next  if !(-e $bam) && -e $bai;
	push(@notfound,$bc) unless $bam ;
	warn "$bc " unless $bam ;
	next  unless -e $bam;
	
	#warn "OK $bam";
	#
 	confess() if exists $dejavu{$bam};
	$dejavu{$bam} ++;
	push(@$mv,{out_bam=>$patient->getBamFileName(),out_bai=>$patient->getBamFileName().".bai",in_bam=>$bam,in_bai=>$bai,bc=>$bc,name=>$patient->name});
}

foreach my $bc (@notfound) {
	
}


my @cmds;

foreach my $f (@$mv){
	my $header = $dir_in."/".$f->{name}."_".$f->{bc}.".txt";
	my $cmd = " samtools view -H ".$f->{in_bam}." | sed  's/SM:".$f->{bc}."/SM:".$f->{name}."/' > $header && ";
	$cmd .= "samtools reheader $header ".$f->{in_bam}." > ".$f->{out_bam}." && samtools index ".$f->{out_bam}." && test -e ".$f->{out_bai} ."|| echo ".$f->{in_bam}." &&  rm ".$f->{in_bam}."*" ;
	push(@cmds,$cmd);

}
print join("\n",@cmds);
print "\n";