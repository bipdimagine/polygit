#!/usr/bin/env perl

use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use lib "$RealBin/../GenBo/lib/obj-nodb/";
warn "$RealBin/../GenBo/lib/obj-nodb/";
#use lib "/data-isilon/bipd-src/pnitschk/git-master/polygit/GenBo/lib/obj-nodb";

use GBuffer;
use Set::IntSpan::Fast::XS;
use Getopt::Long;
use Digest::MD5 qw(md5 md5_hex md5_base64);
my ($project_name, $patient_name,$vid);
my $file;
my $release;
GetOptions(
	'project=s'    => \$project_name,
	'release=s' => \$release,
);

my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $project_name);
warn $project;
my $path = $project->project_path ;

if ($path !~ /HG19/) {
	die($path." not in HG19");
}
my $dir_backup ="/data-pure/sequencing/ngs/backup_HG19/";
$dir_backup ="/data-pure/sequencing/ngs/backup_HG19/".$project->name."/". $project->genome_version ."/";
my $dir_backup_cache ="/data-pure/sequencing/ngs/backup_HG19/".$project->name."/".$project->genome_version."/cache/";
my $dir_backup_cache_tiny ="/data-pure/sequencing/ngs/backup_HG19/".$project->name."/".$project->genome_version."/cache_tiny/";
system("mkdir -p $dir_backup");
system("mkdir  $dir_backup_cache");
system("mkdir  $dir_backup_cache_tiny");
#warn "rclone move $path $dir_backup --progress --transfers 8  --checks 8";
#die();
clean($project);
system ("rclone move $path $dir_backup --progress --transfers 8  --checkers 8 --copy-links && rclone delete $path --transfers 8 --checkers 8");

sub clean {
		my ($project) =@_;
		#my $tr = $project->rocks_cache_2_root_dir()."/".$project->name;
		my $tr2 = $project->tiny_rocks_cache_dir();
		my $cache_directory_actual = $project->getCacheDir();
		
		system("rclone copy ${cache_directory_actual}/ ${dir_backup_cache_tiny} --progress --transfers 8  --checkers 8 && rclone delete ${cache_directory_actual}/ --transfers 8 --checkers 8") if -e ${cache_directory_actual};
		
		
		system("rclone copy ${tr2}/ ${dir_backup_cache_tiny} --progress --transfers 8  --checkers 8 && rclone delete ${tr2}/ --transfers 8 --checkers 8") if -e $tr2;
		
}