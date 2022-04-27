#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb";

use GBuffer;
use GenBoProject;
use Getopt::Long;
use Carp;
use strict;
use Storable qw(store retrieve freeze);

my $project_name;

my %hres;

my $patient_name;
my $version;
# récupère les options  
GetOptions(
	'project=s'   => \$project_name,
	'version=s'   => \$version,
);
die("version man ?") unless $version;
my $fork= 5;
my $bb = GBuffer->new();
my $pp = $bb->newProject( -name => $project_name);
if ($pp->genome_version() eq "HG19_MT"){
#	system("touch $project_name.ok.done");
#	exit(0);
}

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name,-version=>"HG19" );

my $buffer_v = GBuffer->new();
my $project_v = $buffer_v->newProject( -name => $project_name,-version=>$version );
 my $dir_out = $project->getCallingPipelineDir("cmd");
 open(CMD,">$dir_out/cmd1.list");
foreach my $p (@{$project->getPatients}){
	my $patient_name = $p->name;
	print CMD "$Bin/mito_remapping.pl -project=$project_name -patient=$patient_name -version=$version -fork=10 \n";
}
close (CMD);

my $dir1;
my $dir_v;
system("cat $dir_out/cmd1.list | run_cluster.pl -cpu=10");
#warn "***** END MAPPING *******";
warn scalar(@{$project->getPatients});
exit(0);
unlink "$dir_out/cmd1.list";
warn " copy dude";
 $dir1  = $project->getVariationsDir("dude");
 $dir_v  = $project_v->getVariationsDir("dude");
 warn "rsync -ra $dir1/ $dir_v/ 2>/dev/null";
system("rsync -ra $dir1/ $dir_v/ 2>/dev/null");
warn " copy coverage";
$dir1  = $project->getCoverageDir();
$dir_v  = $project_v->getCoverageDir;
system("rsync -rav $dir1/ $dir_v/");
warn " copy cache";
$dir1  = $project->getCacheDir();
$dir_v  = $project_v->getCacheDir();
system("rsync -ra $dir1/ $dir_v/;rm $dir_v/*.cache*");
warn " copy quality $dir_v";
$dir1  = $project->quality_dir();
$dir_v  = $project_v->quality_dir();
system("rsync -rav $dir1/ $dir_v/");
warn " compute cache MT";
update_version($project->id,319);
my $cmd ="echo $Bin/../../scripts_cache/step/launch_chr_cache.pl -project=$project_name -fork=5 -chr=MT -version=$version| run_cluster.pl -cpu=5";
warn $cmd;
system($cmd);
warn " compute dejavu";
my $cmd2 = "echo perl  $Bin/../../../../polymorphism-cgi/cache_nodb/scripts/cache_lite_dejavu.pl -project=$project_name -version=$version";
system($cmd2." | run_cluster.pl -cpu=5 ");
warn " compute polydiag";
my $root_cmd = "perl $Bin/../../../../polymorphism-cgi//cache_nodb_old/scripts/cache_polydiag.pl";
 open(CMD,">$dir_out/cmd2.list");
foreach my $p (@{$project->getPatients}){
	my $cmd = " $root_cmd -fork=2 -project=$project_name -version=$version  -patient=".$p->name;
	print CMD $cmd."\n";
}	
close (CMD);
warn "cat $dir_out/cmd2.list | run_cluster.pl -cpu=5";
system("cat $dir_out/cmd2.list | run_cluster.pl -cpu=5");
unlink "$dir_out/cmd2.list";


system("touch $project_name.done");
exit(0);


sub update_version{
	my ($project_id,$version_id) = @_;
	my $dbh = $buffer->dbh;
	$dbh->{AutoCommit} = 0;
	my $query = qq{
		UPDATE PolyprojectNGS.project_release SET release_id=$version_id WHERE project_id=$project_id;
	};
	$dbh->do($query) ;
	$dbh->commit;
}


#polydude,polydiag,dejavu