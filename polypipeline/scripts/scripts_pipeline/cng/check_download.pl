#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use Getopt::Long;
use Data::Dumper;
use  Parallel::ForkManager;
use Term::ANSIColor;

my $name;
#my $password;
my $dir_in;
#my $type;
my $set;

GetOptions(
	'name=s' => \$name,
	'set=s' => \$set,
);

my @list = `cat /software/polyweb/poly-disk/poly-src/defidiag/project/download/list.txt`;
chomp(@list);
my ($l) = grep{$_ =~ /$name/} @list;
die() unless $l;

my($n,$password,$type) = split(" ",$l);

die("type = devodecode or defidiag") unless $type;
my $dir = "/data-isilon/download/$type/www.cnrgh.fr/data/$name/alignements/$set/";

my @bams = `ls $dir/*.bam`;
chomp(@bams);
warn "verif download MD5  bam";
verif_download_md5(\@bams);
 warn  colored(['bright_green on_black'],"Download MD5 OK");
compute_md5(\@bams);

my $problem = compare_md5(\@bams);
 if(@$problem){
 	warn  colored(['bright_red on_black'],"_______________________________________________")."\n" ;
 	warn  colored(['bright_red on_black'],"-- START DOWNLOAD BAM : ".scalar(@$problem)."--")."\n" ;
 	warn  colored(['bright_red on_black'],"_______________________________________________")."\n" ;
 	sleep(2);
	redownload($problem);# if @$problem;
	compute_md5(\@bams);
	$problem = compare_md5(\@bams);
	warn  colored(['bright_red on_black'],"START COMPLETE DOWNLOAD BAM : ".scalar(@$problem)) ;
 	sleep(2);
	if ($problem){
		redownload($problem);
		compute_md5(\@bams);
		$problem = compare_md5(\@bams);
		if ($problem){
			warn Dumper $problem;
			warn  colored(['bright_red on_black'],"-----------------------------------");
			warn  colored(['bright_red on_black'],"!!!!!!!! YOU HAVE PROBLEM !!!!!!!! ");
			warn  colored(['bright_red on_black'],"-----------------------------------");
			die();
		}
	}
 }
 warn  colored(['bright_green on_black'],"################################")."\n";
 warn  colored(['bright_green on_black'],"#Everything seems OK let's go  #")."\n";
 warn  colored(['bright_green on_black'],"################################")."\n";
exit(0);
sub compute_md5 {
	my ($bams) = @_;
	open (RUN , ">$Bin/$name.$set.check ");
	my $nb=0;
	foreach my $bam (@$bams){
		
		next if -e $bam.".local.md5";
		$nb ++;
		print RUN "md5sum $bam > $bam.local.md5\n";
	}
	close (RUN);
	return if $nb ==0;
	system("cat $name.$set.check | run_cluster.pl -cpu=5");
}

sub verif_download_md5 {
	my ($bams) = @_;
	foreach my $bam (@$bams){
		next if -e -e  "$bam.md5";
			my $bam2 = $bam;
			$bam2 =~ s/\/data-isilon\/download\/$type\///;
			my $cmd_bam="cd /data-isilon/download/$type;wget -c -N --recursive --no-parent --user $name --no-check-certificate  --password   $password https://$bam2.md5 -R cram >/dev/null";  
			warn "download $bam.md5";
			system($cmd_bam);
		die("$bam.md5") unless -e "$bam.md5";
	}
}

sub compare_md5 {
	my ($bams) = @_;
	my @problem;
	foreach my $bam (@$bams){
		my @m1 = `cat $bam.local.md5 $bam.md5 | cut -f 1 -d " "`;
		chomp (@m1);
		die("not md5 $bam") if scalar(@m1)!=2;
		system("touch $bam.ok") if $m1[0] eq $m1[1]; 
		next if $m1[0] eq $m1[1];
		push(@problem,$bam);
	}
	return \@problem;
}

sub redownload {
	my ($bams) = @_;
	return unless @$bams;
	open (RUN , ">$Bin/$name.$set.download ");
	foreach my $bam (@$bams){
		#unlink $bam;
		if (-e "$bam.reload"){
			warn "big probleme with $bam";
			unlink $bam;
		}
		unlink $bam.'.local.md5';
		my $bam2 = $bam;
		$bam2 =~ s/\/data-isilon\/download\/$type\///;
		my $cmd_bam="cd /data-isilon/download/$type;wget -c -N --recursive --no-parent --user $name --no-check-certificate  --password   $password https://$bam2 -R cram ";  
		print RUN $cmd_bam."\n";
		system("touch $bam.reload");
		#warn $cmd_bam;
	}
	system("cat $name.$set.download | run_cluster.pl -cpu=5 -limit=3");
}
	
