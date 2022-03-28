#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../lib/";
use lib "$Bin/../../../lib/GenBoDB";
use lib "$Bin/../../../lib/obj-lite/";
use lib "$Bin/../../../lib/GenBoDB/writeDB";
use lib "$Bin/packages";
use lib $Bin;
#use Set::IntSpan;

use Data::Dumper;
use Getopt::Long;
use Carp;

use Storable qw(store retrieve freeze);
use Term::ANSIColor;
#use Bio::DB::Sam; 
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
my $filein;
my $dir;
my $file_bed;
my $name;
my $build = "HG38";
my $fork = 1;
GetOptions(
	'filein=s'   => \$filein,
	'dir=s'   => \$dir,
	'name=s'=> \$name,
	'bed=s'   => \$file_bed,
	"build=s"  => \$build,
	"fork=s"  => \$fork,
);

#my $file_bed = "/data-xfs/public-data//$build//capture/agilent/agilent.v50.bed";
my $samtools = "/bip-d/soft/bin/samtools";
my $bcftools = "/bip-d/soft/bin/bcftools";
open(BED,"$file_bed");
my $span;
my $span_extended;

while (my $line = <BED>){ 
	chomp($line);
		my @data = split(" ",$line);
		$span->{$data[0]} =  Set::IntSpan::Fast::XS->new() unless exists $span->{$data[0]} ;
		$span->{$data[0]}->add_range($data[1],$data[2]);
		$span_extended->{$data[0]} =  Set::IntSpan::Fast::XS->new() unless exists $span_extended->{$data[0]} ;
		$span_extended->{$data[0]}->add_range($data[1]-150,$data[2]+150);
		
	} 
close BED;


my @chrs = (1..22,'X','Y');

my $bdext = "$dir/$name.bed.extend.bed";
open (BEDEXT,">$bdext") || die();




foreach my $num (@chrs){
	my $chr = "chr".$num;
	next unless exists $span_extended->{$chr};
	my $iter = $span_extended->{$chr}->iterate_runs;
	 while (my ( $from, $to ) = $iter->()) {
        print BEDEXT $chr."\t".$from."\t".$to."\t".abs($from-$to)."\n";
    }
  #  delete  $span_extended->{$chr};
}

#die("problem") if scalar(keys %$span_extended);


;# unless $fork;
#
##creation du projet à partir du buffer
#
#die( "unknown project" . $project_name ) unless $project;

my $File = Thread::Queue->new;
#my @variation_id = map{$_->id} @{$project->getVariations()};


$File->enqueue(@chrs);

my $thr;
for (my $i=0;$i<$fork;$i++){
 	 $thr->[$i] = threads->new(\&by_chr,$span,$span_extended,$name,$filein,$dir);
}


foreach $thr (threads->list) {
        # Ne pas rejoindre le thread principal ni nous-mêmes
        if ($thr->tid && !threads::equal($thr, threads->self)) {
       	 $thr->join;
        }
    }

my $tabfilein;
foreach my $chr (@chrs)    {
	next unless -e  "$dir/$name.$chr.bcf";
	push (@$tabfilein,"$dir/$name.$chr.bcf");
   	  
}

my $cmd = "$bcftools cat ".join(" ",@$tabfilein).">$dir/$name.bcf";
system ("$cmd");

foreach my $file (@$tabfilein){
	warn $file;
	unlink $file;
}
unlink($bdext);
exit(0);   

sub by_chr {
my ($span,$span_extended,$name,$filein,$dir) = @_;
my $bam_file = $filein;
                             
while ( $File->pending() ) {
	my $chr_num =$File->dequeue;
	my $chr = $chr_num;
	 $chr = "chr".$chr_num  if ischrornot($filein);
	my $chr_ucsc = "chr".$chr_num;
	next unless exists $span_extended->{$chr_ucsc};
	my $cmd = "$samtools mpileup -Rugf /data-xfs/public-data//$build//genome/fasta/all.fa -r $chr  $filein -l $bdext -C50 | $bcftools view -bvcg - >$dir/$name.$chr_num.bcf";
	system ("$cmd");
	warn "end $chr";
	
	}
}




sub ischrornot {
	my ($file) = @_;
	my @res = `samtools view -H $file`;

	my ($find) = grep {$_=~ /SN:chr/} @res;
	return $find;
}
