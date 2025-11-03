#!/usr/bin/perl
use CGI qw/:standard :html3/;
use FindBin qw($Bin); 
use strict;
use Getopt::Long;
use List::Util qw(shuffle); 
use Email::Send;
use Email::Send::Gmail;
use Email::Simple::Creator;
use Archive::Tar;
use File::Path qw(remove_tree); 
use lib "$Bin/../../../../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../../../../../GenBo/lib/test/rocks_dejavu/";
use rocks_dejavu_tests;
use Data::Dumper;
use POSIX 'strftime';
 my @chromosomes = (1..22,'X','Y','MT');
 my $filetmp =  "/data-beegfs/tmp/cmd.".time;



open(CMD,">$filetmp");
 my $filename ;
 my $test_obj;
 foreach my $version (("HG19","HG38")){
 	  $test_obj->{$version} = new rocks_dejavu_tests;
	  $test_obj->{$version}->release($version); 	 
 foreach my $chr (@chromosomes){
 		my $date = strftime("%d-%m-%Y", localtime());
 		my $time = time;
 		my $file = "rocks-".$version."-$date.$time";
 		my $dir = $test_obj->{$version}->path_root_prod();
 		$filename->{$version} = $file;
 		print CMD "$Bin/dejavu_rocks.pl -chr=$chr -fork=5 -version=$version -file=$file -dir=$dir\n"; 
 	}
 }
 
 system("cd /data-beegfs/tmp/ ; cat $filetmp | /software/bin/run_cluster.pl -cpu=20 -limit=40 -noprint=1 && touch $filetmp.ok");
 #warn  $filetmp;
 die() unless "$filetmp.ok";


 my @chrs = (1..22, 'X', 'Y','MT' );
 #die() unless -e "$filetmp.ok";
 	foreach my $version (("HG19","HG38")){
 	 	my $dir_final = $test_obj->{$version}->path_root_prod."/".$filename->{$version};
		$test_obj->{$version}->set_directory_path_new($dir_final);
		foreach my $chr (@chrs){
			warn "START ".$chr; 
			$test_obj->{$version}->{use_chromosomes} = [$chr];
			my $is_all_ok = $test_obj->{$version}->launch_test_all_chromosomes();
			die("problem ") if not $is_all_ok;
		}
	}
foreach my $version (("HG19","HG38")){
	my $dir = $test_obj->{$version}->path_root_prod()."/".$filename->{$version} ;
	system("$Bin/list_project.pl -dir=$dir");
}
foreach my $version (("HG19","HG38")){
my $prod_path = $test_obj->{$version}->path_origin;
my $dir_final = $test_obj->{$version}->path_root_prod."/".$filename->{$version};
if ( -l $prod_path ) {  # -l = is symlink
 	my $real = $test_obj->{$version}->path_real();
 	my @t = split("/",$real);
 	my $name = pop(@t);
 	my $dir = join("/",@t);
 	my $tar = $real.".tar";
 	system("tar -cvf $tar -C $dir $name > /dev/null 2>&1");
    unlink $prod_path or warn "cannot unlink $prod_path: $!";
    system("ln -s $dir_final $prod_path");
    unlink $prod_path.".phenotype" or warn "cannot unlink $prod_path: $!";
     system("ln -s $dir_final.phenotype $prod_path.phenotype");
    system("zstd -T20 --rm $real.tar && rm -r $real");
}
else {
	confess();
}
}

exit(0);
#
#foreach my $version (("HG19","HG38")){
#my $backup_path = ;
#my $file = $test_obj->path_root_prod."/variations/rocks.backup.tar.gz";
#system ("tar -zcvf $file -C $backup_path rocks.backup");
#}
 	
exit(0);
 	