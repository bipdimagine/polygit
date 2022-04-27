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
 my @list = `cat /software/polyweb/poly-disk/poly-src/defidiag/project/download/list.txt`;;
chomp(@list);
my ($line) = grep {$_ =~ /$type/} @list; 
die($type."\n".@list) unless $line;


my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $project_name );
my($n,$password,$dir) = split(" ",$line);

my $dir_in = "/data-isilon/download/$dir/www.cnrgh.fr/data/".$type."/";
die($dir_in) unless -e $dir_in;
my @notfound;
my %dejavu;
my $mv;
my $all_bams;
foreach my $patient (@{$project->getPatients}){
	my $bc = $patient->barcode();
	next if -e $patient->getBamFileName();
	my @files = `find $dir_in -name *$bc*`;
	warn "find $dir_in -name *$bc*";
	chomp(@files);
	my ($bam) = grep{$_=~ /.bam$/} @files;
	next if (-e $bam);
	push(@notfound,$bc) ;
	warn "not found ==> $bc " unless $bam ;
}
#warn Dumper @notfound;
#die();
warn "====++====== download miss bam ====++======> ".scalar(@notfound);
reload_bams_bc(\@notfound) if @notfound;
my @list_bams_to_check;
foreach my $patient (@{$project->getPatients}){
	my $bc = $patient->barcode();
	next if -e $patient->getBamFileName();
	my @files = `find $dir_in -name *$bc*`;
	chomp(@files);
	my ($bam) = grep{$_=~ /.bam$/} @files;
	warn $bam;
	next if -e "$bam.ok";
	die() unless -e $bam;
	unless (-e  "$bam.md5"){
			my $bam2 = $bam;
			$bam2 =~ s/\/data-isilon\/download\/$type\///;
			my $cmd_bam="cd /data-isilon/download/$type;wget -c -N --recursive --no-parent --user $type --no-check-certificate  --password   $password https://$bam2.md5 -R cram";  
			system($cmd_bam);
	}
	push(@list_bams_to_check,$bam) ;
	#	warn "not found ==> $bc " unless $bam ;
}

warn "====++====== check  bam ====++====== ".scalar(@list_bams_to_check);



my $reload  = check_bams(\@list_bams_to_check,1);
warn "====++====== reload md5  bam ====++====== ".scalar(@$reload);

 if (@$reload){
	reload_bams($reload);
	my $reload2  = check_bams($reload,1);
	confess(Dumper $reload2) if @$reload2;
 }



foreach my $patient (@{$project->getPatients}){
	my $bc = $patient->barcode();
	next if -e $patient->getBamFileName();
	#warn $bc;
	my @files = `find $dir_in -name *$bc*`;
	chomp(@files);
	my ($bam) = grep{$_=~ /.bam$/} @files;
	#my ($bai) = grep{$_=~ /.bai$/} @files;
	
	#my $md5out = `md5sum $dir_in/$bamout | cut -f 1 -d " "`;
	#chomp($md5out);
	#warn $md5out;
	
	#die("$md5out- $md5in-") if $md5out ne $md5in;
	my $bai = $patient->getBamFileName().".bai"; 
	next if -e $bai && -e $patient->getBamFileName();
	next  if !(-e $bam) && -e $bai && $patient->getBamFileName();
	push(@notfound,$bc) unless $bam ;
	warn "$bc " unless $bam ;
	next  unless -e $bam;
	
	#warn "OK $bam";
	#
 	confess() if exists $dejavu{$bam};
	$dejavu{$bam} ++;
	push(@$mv,{out_bam=>$patient->getBamFileName(),out_bai=>$patient->getBamFileName().".bai",in_bam=>$bam,in_bai=>$bai,bc=>$bc,name=>$patient->name});
}
warn Dumper @notfound if @notfound;
die() if @notfound;


my @cmds;

foreach my $f (@$mv){
	my $header = $dir_in."/".$f->{name}."_".$f->{bc}.".txt";
	my $cmd = " samtools view -H ".$f->{in_bam}." | sed  's/SM:".$f->{bc}."/SM:".$f->{name}."/' > $header && ";
	#$cmd .= "samtools reheader -\@5 $header ".$f->{in_bam}." > ".$f->{out_bam}."  && samtools index -\@5 ".$f->{out_bam}."  && test -e ".$f->{out_bai} ."|| echo ".$f->{in_bam}." &&  rm ".$f->{in_bam}."*" ;
	 $cmd .= "samtools reheader $header ".$f->{in_bam}." > ".$f->{out_bam}."  && samtools index -\@5 ".$f->{out_bam}."  && test -s ".$f->{out_bai} ."|| test -s ".$f->{in_bam}." &&  rm ".$f->{in_bam}."*" ;
	
	push(@cmds,$cmd);

}

print join("\n",@cmds);



sub check_bams {
	my ($bams,$force) = @_;
	my $pm = new Parallel::ForkManager(5);
	my @reload;
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

foreach my $bam (@$bams){
	#next unless $bam =~/C001NBU/;
	#warn $bam;
	next if -e "$bam.ok";
	warn "$bam.md5" unless -e "$bam.md5";
	unless (-e  "$bam.md5"){
			my $bam2 = $bam;
			$bam2 =~ s/\/data-isilon\/download\/$type\///;
			my $cmd_bam="cd /data-isilon/download/$type;wget -c -N --recursive --no-parent --user $type --no-check-certificate  --password   $password https://$bam2.md5 -R cram";  
			system($cmd_bam);
	}
	die("$bam.md5") unless -e "$bam.md5";
	my $pid = $pm->start and next;
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
		#warn $bam;
		#die();
		unlink $bam if ($force);
		#my $cmd_bam="/data-isilon/download/$type;wget -c -N --recursive --no-parent --user $name --no-check-certificate  --password   $password https://$bam -R cram";  
		$res->{reload} = $bam;
		#push(@reload,$bam);
	}
	else {
		system("touch $bam.ok");
	}
	$pm->finish(0,$res);
	}
	$pm->wait_all_children();
	
	return \@reload;
}

sub reload_bams_bc {
	my ($bcs) = @_;
	my $pm2 = new Parallel::ForkManager(3);
	foreach my $bc (@$bcs) {
	#my $pid = $pm2->start and next;
	print ("./download_one_file.pl -name=$type -set=$set -bc=$bc \n");
	 #system(qq{"/data-isilon/download/$type;wget -c -N --recursive --no-parent --user $type --no-check-certificate  --password   $password https://$bam -R cram});
	#$pm2->finish(0,{});
	}
	#$pm2->wait_all_children();
	die();
}
sub reload_bams {
	my ($bams) = @_;
	my $pm2 = new Parallel::ForkManager(3);
	foreach my $bam (@$bams) {
	my $pid = $pm2->start and next;
	 #system(qq{"/data-isilon/download/$type;wget -c -N --recursive --no-parent --user $type --no-check-certificate  --password   $password https://$bam -R cram});
	 my $bam2 = $bam;
		$bam2 =~ s/\/data-isilon\/download\/$type\///;
		my $cmd_bam="cd /data-isilon/download/$type;wget -c -N --recursive --no-parent --user $type --no-check-certificate  --password   $password https://$bam2.md5 -R cram";  
		system($cmd_bam);
	$pm2->finish(0,{});
	}
	$pm2->wait_all_children();
	
}

print "\n";