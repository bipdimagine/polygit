#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use GBuffer;
use GenBoProject;
use Getopt::Long;
use Data::Dumper;
#use illumina_util;
use IO::Prompt;
use Sys::Hostname; 
 use File::Find::Rule ;
use Text::Table;
use Term::Twiddle;

my $project_name_1;
my $chr_name;
my $ppn;
my $set;
my $name;
my $version;

GetOptions(
	'project=s' => \$project_name_1,
	'chr=s' => \$chr_name,
	'fork=s' => \$ppn,
	'set=s' => \$ppn,
	'fork=s' => \$ppn,
	'version=s' => \$version,
	#'fork=s' => \$ppn,
);
$ppn=20;
my @list = split(",",$project_name_1);
foreach my $project_name (@list){
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name,-version=> $version );
my $chr = $project->getChromosome($chr_name);
$chr_name = $chr->name;
my $hash_cmd;
my $cmd_perl = "/usr/bin/perl $Bin/../../../../polymorphism-cgi/cache_nodb/scripts/rocks";
my $tmp_dir = "/data-beegfs/tmp/log/".$project_name;
system("mkdir -p $tmp_dir") unless -e $tmp_dir;

push(@$hash_cmd, {
	cmd => "$cmd_perl/cache_store_ids.pl -project=$project_name -chr=$chr_name -fork=$ppn ",
}
);
push(@$hash_cmd, {
	cmd => "$cmd_perl/cache_store_annotations.pl -project=$project_name -chr=$chr_name -fork=$ppn ",
}
);



push(@$hash_cmd, {
	cmd => "$cmd_perl/cache_strict_denovo.pl  -project=$project_name -chr=$chr_name -fork=$ppn ",
}
);

push(@$hash_cmd, {
	cmd => "$cmd_perl/polyviewer.pl -project=$project_name -chr=$chr_name -fork=$ppn ",
}
);

my $id = int(rand(1000) + time);;
my $nb = 0;
foreach my $hc (@$hash_cmd){
	 $id += int(rand(50));
	 $nb++;
	 my $file_tmp = $tmp_dir."/$id".".".$chr_name.".".$nb;
	my $cmd = $hc->{cmd};
	#next if -e  $hc->{fileout};
	my $error = $file_tmp.".error";
	my $ok = $file_tmp.".ok";
	my $time = time;
	warn "**** start ****";
	warn $cmd;
	system($cmd."   \|\| touch $error  \&\& touch $ok");
	
	warn "end ==> ".abs(time-$time);
	
	if (-e $error){
		warn $error;
		die($cmd);
		unlink $error;
	}
	die($name) unless -e $ok;
	unlink $ok
} 
}
