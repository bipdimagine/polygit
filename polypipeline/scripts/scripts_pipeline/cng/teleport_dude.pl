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
my $project_name_trio;
my $project_name_origin;
my $filename;

my $name;
my $patients_name;
my $steps_name;
my $force;
my $type;
my $fastq_ext;
GetOptions(
	'project_trio=s' => \$project_name_trio,
	'project_origin=s' => \$project_name_origin,
);




my $buffer_solo = GBuffer->new();
my $project_solo = $buffer_solo->newProject( -name => $project_name_origin );
my $patients_solo = $project_solo->getPatients();
warn scalar(@$patients_solo);
my $buffer_trio = GBuffer->new();
my $project_trio = $buffer_trio->newProject( -name => $project_name_trio );
my $patients_trio = $project_trio->getPatients();
warn scalar(@$patients_trio);

my @common ;



my $primer_trio  = $project_trio->noSqlCnvsDir."/primers.bed.gz";
my $primer_solo  = $project_solo->noSqlCnvsDir."/primers.bed.gz";
warn $primer_solo;
unless (-e $primer_trio){
		my $cmd = "rsync -rav ".$project_solo->noSqlCnvsDir."/ ".$project_trio->noSqlCnvsDir."/";
		warn $cmd;
		system($cmd);
		#push(@allcmd,$cmd);
}

#foreach my $cmd (@allcmd){
#	system($cmd);
#}
exit(0);
#my $primers_dir = $project->getCacheDir() . "/coverage_lite/";
#my $ped_ori =  $project->getPedigreeFile();
#my @no_moove;
#my $moove;
#my $dbh = $buffer->dbh;
#$dbh->{AutoCommit} = 0;
#my @allcmd;
#foreach my $p (@{$project->getPatients}){
#	my $dest_project = $buffer->getQuery->getProjectDestination($p->id);
#	if (!$dest_project || $dest_project eq $projectName){
#		push(@no_moove,$p);
#	} 
#	else {
#		push(@{$moove->{$dest_project}},$p);
#	}
#}
#IO::Prompt::hand_print("hummm humm \n\n ");
#if (scalar@no_moove){
#	
#	IO::Prompt::hand_print("First it's seems that you don't want  move this patients ?\n");
#	foreach my $p (@no_moove){
#		colored::stabilo("white",$p->name()); 
#	}
#}
#
#IO::Prompt::hand_print(" let's start with  ?\n");
#
#foreach my $pdes (keys %$moove){
#	colored::stabilo("yellow","Project ".$pdes); 
#	my $pname = $project->name;
#	
#	foreach my $p (@{$moove->{$pdes}}){
#		print "\t";
#		colored::stabilo("green",$p->name); 
#		
#	}
#	die() unless  (prompt(" still sure for this one (y/n)",-yes));
#	my $buffer2 = GBuffer->new();
#	my $newProject = $buffer2->newProject(-name=>$pdes);
#	my $new_project_id = $newProject->id();
#	
#	my $project_id = $project->id;
#	my $sql = qq{
#		UPDATE `PolyprojectNGS`.`patient` SET `project_id`=$new_project_id  WHERE `patient_id`=? and`project_id`=$project_id;
#		};
#foreach my $p (@{$moove->{$pdes}}){
#		my $pid = $p->id;
#		colored::stabilo("blue"," move : ".$p->name); 
#		my $sth= $dbh->prepare($sql);
#		push(@allcmd,move_patient($project,$newProject,$p));
#			my $new_primers_out = $newProject->noSqlCnvsDir."/primers.bed.gz";
#			my $primers  = $project->noSqlCnvsDir."/primers.bed.gz";
#		if (-e $primers && -e $new_primers_out){
#			colored::stabilo("red","----------------CHECK DUDE ---------------");
#		}
#		else {
#			if (-e $primers) {
#				my $cmd = "rsync -rav".$newProject->noSqlCnvsDir."/ ".$newProject->noSqlCnvsDir."/";
#			}
#		}
#		$sth->execute($pid);
#		$sth->finish;
#	}
#}
#IO::Prompt::hand_print("Commit database change ?");
#unless (prompt("(yes/no)",-yes)){
#	die(" back to the future");
#}
#colored::stabilo("magenta"," beam me up scotty "); 
#	
#	my $zz =0;
#my $spinner = new Term::Twiddle;
#colored::stabilo("red"," warning don't kill or stop here I'm MOOVING the files !!!!"); 
#$spinner->start;
#  
#
#	foreach my $cmd (@allcmd){
#		system($cmd);
#	}
#$spinner->stop;
#$dbh->commit(); #if prompt("(yes/no)",-yes);
#colored::stabilo("green","------------------------------------"); 
#IO::Prompt::hand_print(" Done \n");
#colored::stabilo("green","------------------------------------"); 
#print "\n";
#colored::stabilo("yellow","you can now run the cache and the quality  :".join(",", keys %$moove));
#colored::stabilo("yellow","lists of projects :".join(",",keys %$moove));
#colored::stabilo("green","$Bin/after_all.pl -projects=".join(",",keys %$moove));
#exit(0);
#	






sub move_patient{
	my ($project,$newProject,$patient) =@_;
	my $project_in = $project->name;
	my $projectName = $newProject->name;	
	my $pn = $patient->name();
 	my $dir1 = $project->getProjectPath();
 	my $dir2 = $newProject->getProjectPath();
 	my $dir3 = $project->getSequencesRootDir();

 #my @files = `find $dir1 -name $pn.*`;
  my @files = File::Find::Rule->file()
                                  ->name( "$pn.*" )
                                  ->in( $dir1);
                                  
  my (@bam) = grep {/bam$/} @files;
   my (@vcf) = grep {/vcf.gz$/} @files;
   warn scalar @vcf;
  
 if (scalar(@bam) ne 1 || scalar(@vcf) < 2){
 	colored::stabilo("red","----------------$pn problem with files---------------");
	colored::stabilo("red","----------------".join(",",@bam)." bam file  ---------------");
	colored::stabilo("red","---------------- ".join(",",@vcf)." vcf  file  ---------------");
	colored::stabilo("red","-----------------------------------------------------"); 
	warn Dumper @files;
 	
 }
 
my @files2 = File::Find::Rule->file()
                                  ->name( "$pn.*" )
                                  ->in( $dir3);


push(@files,@files2);
my $bam;
my $vcf;
my @cmd;
foreach my $f (@files){
	my $f2 = $f;
	$bam ++ if $f2 =~ /bam$/;
	$vcf ++ if $f2 =~ /vcf.gz$/;
	$f2 =~   s/$project_in/$projectName/;
	warn ("file exists $f2") if -e $f2;
	die ("file exists $f2 $f") if -e $f2;
	my @dir = split("/",$f2);
	pop(@dir);
	
	my $dir_out = join("/",@dir);
	system("mkdir -p $dir_out") unless -e $dir_out;
	
	die("no directory  $dir_out") unless -e $dir_out;
	push(@cmd,"ln -s $f $f2");
	
} 

if ($bam ne 1 || $vcf < 2){
	colored::stabilo("red","----------------$pn problem with files---------------");
	colored::stabilo("red","----------------$bam bam file  ---------------") if $bam ne 1; 
	colored::stabilo("red","---------------- $vcf vcf  file  ---------------") if $vcf < 2; 
	colored::stabilo("red","-----------------------------------------------------"); 
	warn Dumper @cmd;
	die();
}
	
	
 return @cmd;
}
exit(0);