#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use Getopt::Long;
use Data::Dumper;
use  Parallel::ForkManager;

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
warn 'end';
my @reload;

my $pm = new Parallel::ForkManager(5);
$pm->run_on_finish(
    sub { 
    	my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$h)=@_;
  
    	unless (defined($h) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				die();
				return;
			}
			if (exists $h->{reload}){
				push(@reload, $h->{reload});
			}
    }
	);

foreach my $bam (@bams){
	#next unless $bam =~/C001NBU/;
	#warn $bam;
	next if -e "$bam.ok";
	my $pid = $pm->start and next;
	next unless -e "$bam.md5";
	my @t = `cat $bam.md5`;
	chomp(@t);
	my ($md5_2,$f) = split(" ",$t[0]);
	my $s = `md5sum $bam`;
	chomp($s);
	my ($md5,$f2) = split(" ",$s);
	my $res = {};
	if ($md5 ne $md5_2){
		warn $md5.' '.$md5_2;
		warn $bam;
		$bam =~ s/\/data-isilon\/download\/$type\///;
		my $cmd_bam="wget -c -N --recursive --no-parent --user $name --no-check-certificate  --password   $password https://$bam -R cram";  
		warn $cmd_bam;
		#system($cmd_bam);
		$res->{reload} = $cmd_bam;
		#push(@reload,$bam);
	}
	else {
		system("touch $bam.ok");
	}
	$pm->finish(0,$res);
	}
	$pm->wait_all_children();
	
	
my $pm2 = new Parallel::ForkManager(2);
foreach my $cmd (@reload) {
	my $pid = $pm2->start and next;
	 system($cmd);
	$pm2->finish(0,{});
	}
	$pm2->wait_all_children();

exit(0);

