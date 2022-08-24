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
my $solo;
GetOptions(
	'project=s' => \$project_name,
	'fork=s' => \$ppn,
	'set=s' => \$set,
	'name=s' => \$name,
	'solo=s' => \$solo,
);
my @list;
my $dir2 = "";
if($solo){
	$dir2 ="solo";
}
$set = "set".$set unless $set=~/set/;
if($name){
 @list = `cat ../../../../defidiag/project/$name/$set.txt`;
 chomp(@list);
}
elsif($project_name){
	push(@list,split(",",$project_name));
}
die() unless @list;

foreach my $project_name (@list){
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );

foreach my $chr (@{$project->getChromosomes}){
	my $fork =20;
	
	if ($chr->karyotypeId<3){
		$fork =40;
	}
	elsif ($chr->karyotypeId<15) {
		$fork = 20;
	}
	elsif ($chr->karyotypeId<22) {
		$fork = 10;
	}
	elsif ($chr->karyotypeId==23) {
		$fork = 20;
	}
	else {
		$fork = 5;
	}
	my $cmd ="$Bin/launch_chr_cache.pl -project=$project_name -fork=$fork -chr=";
	print $cmd.$chr->name." :::$fork\n";
}
unless ($project->isGenome){
	
	my $bin_dev = "$Bin/../../scripts_pipeline";
	my $ppn = 20;
	foreach my $patient (@{$project->getPatients}){
		my $name = $patient->name;
		my $no1 = $patient->getTranscriptsDude();
		my $cmd = qq{perl $bin_dev/transcripts/transcripts_cache.pl -patient=$name  -fork=$ppn  -project=$project_name  && perl $bin_dev/transcripts/genes_level_dude.pl -patient=$name  -fork=$ppn  -project=$project_name};
		print $cmd." :::$ppn\n";
	}
	

}


}
