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

my $project_name;
my $chr_name;
my $ppn;
my $set;
my $name;
my $version;

GetOptions(
	'project=s' => \$project_name,
	'chr=s' => \$chr_name,
	'fork=s' => \$ppn,
	'set=s' => \$ppn,
	'fork=s' => \$ppn,
	'version=s' => \$version,
	#'fork=s' => \$ppn,
);
my @list;
if($name){
 @list = `cat ../../../../../defidiag/project/$set/$name.txt`;
 chomp(@list);
}
elsif($project_name){
	push(@list,$project_name);
}
die() unless @list;


foreach my $project_name (@list){
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name,-version=> $version );
my $chr = $project->getChromosome($chr_name);
$chr_name = $chr->name;
my $hash_cmd;
my $cmd_perl = "/usr/bin/perl $Bin/../../../../polymorphism-cgi/cache_nodb/scripts";
push(@$hash_cmd, {
	cmd => "$cmd_perl/cache_store_ids.pl -project=$project_name -chr=$chr_name -fork=$ppn",
	fileout => $project->getCacheBitVectorDir()."/lmdb_cache/".$chr_name.".dv.freeze",
}
);

push(@$hash_cmd, {
	cmd => "$cmd_perl/cache_store_annotations.pl -project=$project_name -chr=$chr_name -fork=$ppn ",
	fileout => $project->getCacheBitVectorDir()."/lmdb_cache/".$chr_name."/genes_index",
}
);


push(@$hash_cmd, {
	cmd => "$cmd_perl/cache_check_step.pl -step=store_annotations  -project=$project_name -chr=$chr_name -fork=$ppn ",
	fileout => $project->getCacheBitVectorDir()."/log/check_store_annotations.$chr_name.ok",
}
);


push(@$hash_cmd, {
	cmd => "$cmd_perl/cache_strict_denovo.pl  -project=$project_name -chr=$chr_name -fork=$ppn ",
	fileout => $project->getCacheBitVectorDir()."/strict-denovo/$chr_name.lite",
}
);

push(@$hash_cmd, {
	cmd => "$cmd_perl/cache_check_step.pl -project=$project_name -step=strict_denovo -project=$project_name -chr=$chr_name -fork=$ppn && touch  ".$project->getCacheBitVectorDir()."/strict-denovo/$chr_name.lite.check",
	fileout => $project->getCacheBitVectorDir()."/strict-denovo/$chr_name.lite.check",
}
);
push(@$hash_cmd, {
	cmd => "$cmd_perl/polyviewer.pl -project=$project_name -chr=$chr_name -fork=$ppn  ",
	fileout => $chr->lmdb_cache_dir()."/lmdb.ok",
}
);
my $tmp_dir = "/data-beegfs/tmp/log/".$project_name;
system("mkdir -p $tmp_dir") unless -e $tmp_dir;
my $id = int(rand(1000) + time);;
my $nb = 0;
foreach my $hc (@$hash_cmd){
	 $id += int(rand(50));
	 $nb++;
	 my $file_tmp = $tmp_dir."/$id".".".$chr_name.".".$nb;
	my $cmd = $hc->{cmd};
	warn $hc->{fileout};
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
	
	die($cmd." ".$hc->{fileout}) unless -e  $hc->{fileout};
	die("cocucuou") unless -e $ok;
	unlink $ok
} 
}
