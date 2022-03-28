#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
#use Set::IntSpan;
use GBuffer;
use Data::Dumper;
use Getopt::Long;
use Carp;



my $file2;
my $file1;
my $patient_name;
my $buffer = GBuffer->new();
GetOptions(
	'file1=s'   => \$file1,
	'file2=s'   => \$file2,
	'patient=s'   => \$patient_name,
);
my $samtools = $buffer->software("samtools");
my ($a) = `samtools idxstats $file1 | awk '{ sum += \$3 } END { print sum }'`;
my ($b) = `samtools idxstats $file2 | awk '{ sum += \$3 } END { print sum }'`;
chomp($a);
chomp($b);

problem("problem with bam file not same amount of read $file1 $file2") if $a ne $b;

warn "verif 1 OK sale sequence => $a $b";

my @t1 = `$samtools view -H $file1 | grep "PN:bwa"`;

chomp(@t1);
my @cn1 = grep{/gz/} map{split(" ",$_)} @t1;

my $tt1 = join("-",@cn1);
my @t2 = `$samtools view -H $file2 | grep "PN:bwa"`;
chomp(@t2);
my @cn2 = grep{/gz/} map{split(" ",$_)} @t2;
my $tt2 = join("-",@cn2);
problem("\n".$tt1." \n$tt2 ") if $tt1 ne $tt2;
warn "verif 2 OK same sequence ";
problem($tt2) unless $tt2 =~/$patient_name/;

warn "verif 3 OK same sequence";

my ($z) = `$samtools view -H $file2 | grep "^\@RG" | grep "SM:$patient_name" | wc -l`;
chomp $z;
problem("problem read group $z") if $z ne 1;
my $nb = 3000 + int(rand(1000));
my $cmd_head = "head -".$nb;
my ($md1) = `samtools view -f 19 $file1 | cut -f1 | $cmd_head | sort -u | md5sum`;
my ($md2) = `samtools view -f 19 $file2 | cut -f1 | $cmd_head | sort -u | md5sum`;
chomp($md1);
chomp($md2);
warn $md1." $nb ".$md2;
 $nb = 4000 +  int(rand(1000));;
 $cmd_head = "head -".$nb;
if ($md1 ne $md2){
 ($md1) = `samtools view -f 19 $file1 | cut -f1 | $cmd_head | sort -u | md5sum`;
 ($md2) = `samtools view -f 19 $file2 | cut -f1 | $cmd_head | sort -u | md5sum`;
chomp($md1);
chomp($md2);
}
problem($md1." ".$md2) if $md1 ne $md2; 
warn "test read name OK \n $md1 $md2 \n $file1 $file2";
#my @aa = `samtools view -f 3  $file1 | cut -f1 | head -200 | sort -b`;
#chomp(@aa);
#my @bb = `samtools view -f 3 $file2 | cut -f1 | head -200 | sort -b`;
#chomp(@bb);
#for (my $i =0;$i<@aa;$i++){
#	print $i." ".$aa[$i]." ".$bb[$i]."\n" if $aa[$i] ne $bb[$i];
#	
#}
warn "Elprep OK !!!! ";

sub problem {
	my ($tex) = @_;
	warn $tex;
	system("mv $file2 $file2.bad");
	die($tex);
}
