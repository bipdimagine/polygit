#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use POSIX;
use Term::Spinner;
use Data::Dumper;
use File::Temp qw/tempfile tempdir /;
use File::Slurp;
use String::ProgressBar;
use Term::StatusBar;
use Parallel::ForkManager;
 use Digest::MD5 qw(md5 md5_hex md5_base64);
use Term::ANSIColor;
 use warnings;
use IPC::Run3 'run3';
use Time::Piece;
use Time::Seconds qw/ ONE_DAY /;
use lib "$Bin/../../GenBo/lib/";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../../GenBo/lib/kyoto/";
use lib "$Bin/../../GenBo/lib/GenBoDB/writeDB";
use lib "$Bin/../packages";
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moose;
use GBuffer;
use GenBoProject;
use colored; 
use Text::Table;
use file_util;
use File::Temp qw/ tempfile tempdir /;; 
use Term::Menus;
use Proc::Simple;

my @running_jobs;
my $cluster_jobs;
my $runnings_jobs;
my $sentences = [];
my $file_cluster = "/tmp/tmp.".time."cmd";
my $script_perl = $Bin."/scripts/";
my $script_pipeline = $Bin."/../scripts/scripts_pipeline/";
my $num_jobs =0;
$SIG{INT} = \&tsktsk;
$SIG{KILL} = \&tsktsk;
sub tsktsk {
	warn "kill jobs";
  	foreach my $j (@running_jobs){
  		$j->kill() if $j->poll();
  	}
  	foreach my $c (keys %$runnings_jobs){
		unlink "slurm-".$c."out";
		system("scancel ".$c);
	}
}

my $project_name;
my $patients_name;
my $limit;
GetOptions(
	'project=s' => \$project_name,
	'patients=s' => \$patients_name,
	#'low_calling=s' => \$low_calling,
);
my $steps = {
				"dragen-alignment"=> \&run_align,
				"move"=>  \&run_move,
				"genotype"=>  \&run_genotype,
				"dragen-sv"=>  \&run_sv,
				"dragen-target"=>  \&run_target,
				"dragen-pon"=>  \&run_pon,
				"lmdb-depth"=>  \&run_lmdb_depth,
				"run_coverage"=>  \&run_coverage,
				"run_dragen_cnv_coverage" =>\&run_dragen_cnv_coverage,
				
};

my $buffers;
my $projects;
my @apatients_name = split(":",$patients_name);
foreach my $pname (split(",",$project_name)){
	my $buffer = GBuffer->new();
	my $project = $buffer->newProject( -name => $pname );
	$project->isGenome;
	$project->get_only_list_patients($apatients_name[0]);
	push(@$projects,$project);
	push(@$buffers,$buffer);
}
 $project_name =~ s/\,/\./g ;
 warn $projects->[0]->buffer;
my $dir_log = $projects->[0]->buffer->config->{project_pipeline}->{bds}."/".$project_name.".dragen.".time;
system("mkdir $dir_log && chmod a+rwx $dir_log");




 $|=1;
 
system("clear");# unless $noprint;

my $start_time = time;
my $jobs =[];

####### Alignement


run_align($projects);

run_move($projects);
run_gvcf($projects);
run_dragen_cnv_coverage($projects);
run_genotype($projects);
#run_sv();
#run_target();
run_pon($projects);
run_coverage($projects);
run_lmdb_depth($projects);
exit(0);



####### Genotype

run_genotype();
run_lmdb_depth();
run_sv_pon();

warn "end";
exit(0);



sub running_text {
	my ($text,$row) = @_;
	print   "\033[".$row.";0H",colored(['bright_white'],"$text :").colored(['bright_cyan on_black']," RUNNING ")."\n";
}



########################################################################################################
#########################---------------------  PIPELINE CLUSTER --------------##########################
########################################################################################################

### MOVE 
sub run_move {
	my ($projects) = @_;
	my $jobs =[];

foreach my $project (@$projects){
	my $ppn = 5;
	my $projectName = $project->name;
	my $patients = $project->getPatients($patients_name);
	foreach my $patient (@$patients){
		my $dir_pipeline = $patient->getDragenDir("pipeline");
		my $prefix = $patient->name;
		my $bam_pipeline = $dir_pipeline."/".$prefix.".bam";
		next if -e $patient->getBamFileName("dragen-align");
		die() unless  -e $bam_pipeline;
		my $cmd = "perl $script_perl/dragen_move.pl -project=$projectName -patient=".$patient->name;
		my $job;
		$job->{name} = $patient->name.".lmdb";
		$job->{cmd} =$cmd;
		$job->{cpus} = $ppn;
		push(@$jobs,$job);
	}
}
	
	my $text = "$num_jobs- MOVE BAM";
	my $limit = 5;
	steps_cluster("MOVE BAM  ",$jobs,$limit);

	
}


### LMDB 
sub run_lmdb_depth {
	my ($projects) = @_;
	my $jobs =[];

foreach my $project (@$projects){
	
	my $ppn = 10;
	my $project_name = $project->name;
	my $patients = $project->getPatients($patients_name);
	foreach my $patient (@$patients){
		my $name = $patient->name;
		my $fileout = $patient->fileNoSqlDepth;
		next if -e $fileout;
		my $cmd;
		 $cmd = qq{perl $script_pipeline/coverage_genome.pl -patient=$name  -fork=$ppn  -project=$project_name  };
		if ($project->isGenome){
			$cmd .= qq{ && perl $script_pipeline/coverage_statistics_genome.pl -patient=$name  -fork=$ppn  -project=$project_name};
		}
		else {
		#	 $cmd = qq{perl $script_pipeline/coverage.pl -patient=$name -filein=$filein -dir=$coverage_dir -bed=$bed -fork=$ppn -name=$name -project=$project_name :::$ppn};
		}
		my $job;
		$job->{name} = $patient->name.".lmdb";
		$job->{cmd} =$cmd;
		$job->{cpus} = $ppn;
		push(@$jobs,$job);
	}
	}
	steps_cluster("LMDBDepth ",$jobs);
	
}


sub run_coverage {
	my ($projects) = @_;
	my $jobs =[];

foreach my $project (@$projects){
	next if $project->isGenome;
	$jobs= [];
	my $ppn = 10;
	my $project_name = $project->name;
	my $coverage_dir = $project->getRootDir() . "/align/coverage/";
	system("mkdir $coverage_dir && chmod a+rwx $coverage_dir") unless -e $coverage_dir;
	my $patients = $project->getPatients($patients_name);
	foreach my $patient (@$patients){
		my $name = $patient->name;
		
		my $fileout = $coverage_dir . "/" . $name.".cov.gz";
		next if -e $fileout;
		my $bed = $patient->getCaptureFile();
		my $filein = $patient->getBamFileName();
		
		 my $cmd = qq{perl $script_pipeline/coverage.pl -patient=$name -filein=$filein -dir=$coverage_dir -bed=$bed -fork=$ppn -name=$name -project=$project_name };
		my $job;
		$job->{name} = $patient->name.".coverage";
		$job->{cmd} =$cmd;
		$job->{cpus} =$ppn;
		push(@$jobs,$job);
	}
}
	steps_cluster("Coverage  ",$jobs);
}



########################################################################################################
#########################---------------------  PIPELINE DRAGEN --------------##########################
########################################################################################################

##########################
### ALIGN 
##########################

sub run_align {
my ($projects) = @_;
my $jobs =[];

foreach my $project (@$projects){
my $patients = $project->getPatients($patients_name);
my $projectName = $project->name;
foreach my $patient (@$patients){
	my $dir_pipeline = $patient->getDragenDir("pipeline");
	my $prefix = $patient->name;
	my $bam_pipeline = $dir_pipeline."/".$prefix.".bam";
	#next if -e $patient->getBamFileName("dragen-align");
	next if -e $patient->getBamFileName();
	next if -e $bam_pipeline;
	my $job;
	$job->{name} = $patient->name.".aln";
	$job->{cmd} = "perl $script_perl/dragen_align_calling.pl -project=$projectName -patient=".$patient->name;
	$job->{out} =  $dir_pipeline."/".$prefix.".bam";
	push(@$jobs,$job);
	}
}
	
	my $text = "$num_jobs-  DRAGEN ALIGN";
	steps_system("Dragen Align",$jobs);
	return;
}

##########################
### GENOTYPE 
##########################
sub run_genotype {
my ($projects) = @_;
	my $jobs =[];
foreach my $project (@$projects){
	my $projectName = $project->name;
	die();
 my @ps;
 my $patients = $project->getPatients();
 foreach my $p (@$patients) {
 	my $dir_out= $project->getVariationsDir("dragen-calling");
 	my $f = $dir_out."/".$p->name.".vcf.gz";
 	next if -e $f;
 	push(@ps,$p->name);
 }
 if(@ps){
	my $cmd_genotype = "perl $script_perl/dragen_genotype.pl -project=$projectName -patient=".join(",",@ps);
	push(@$jobs,{cmd=>$cmd_genotype,name=>$project->name.".genotype"});
 }
}
 	steps_system("Dragen Genotype ",$jobs);
}
##########################
#         GVCF
##########################
sub run_gvcf {
my ($projects) = @_;
my $jobs =[];

foreach my $project (@$projects){
	my $projectName = $project->name;
	my $patients = $project->getPatients($patients_name);
	my $cmd = qq{perl $script_perl/dragen_gvf.pl -project=$projectName};
	foreach my $patient (@$patients){
		my $fileout = $patient->gvcfFileName("dragen-calling");;
		next if -e $fileout;
			my $cmd = qq{perl $script_perl/dragen_gvcf.pl -project=$projectName -patient=}.$patient->name;
			push(@$jobs,{name=>$patient->name.".gvcf" ,cmd=>$cmd,fileout=>$fileout});
	}
}
steps_system("GVCF",$jobs);
}

##########################
#         CNV et COVERAGE
##########################
sub run_dragen_cnv_coverage {
	my ($projects) = @_;
	my $jobs =[];

foreach my $project (@$projects){
	my $projectName = $project->name;
	my $cmd = qq{perl $script_perl/dragen_cnv_coverage.pl -project=$projectName};
	my $patients = $project->getPatients($patients_name);
	foreach my $patient (@$patients){
		my $dir = $patient->project->getTargetCountDir();
		my $fileout = $dir."/".$patient->name.".target.counts.gc-corrected.gz";
		next if -e $fileout;
		my $cmd = qq{perl $script_perl/dragen_cnv_coverage.pl -project=$projectName -patient=}.$patient->name;
		push(@$jobs,{name=>$patient->name.".cnv_cov" ,cmd=>$cmd,fileout=>$fileout});
	}
}
steps_system("CNV_COV",$jobs);
}
##########################
#         CNV AND SV
##########################

sub run_pon {
	my ($projects) = @_;
	my $jobs =[];

foreach my $project (@$projects){
	my $projectName = $project->name;
	my $dir_prod = $project->getVariationsDir("dragen-pon");
	my $patients = $project->getPatients($patients_name);
	my $cmd = qq{perl $script_perl/dragen_cnv.pl -project=$projectName};
	foreach my $patient (@$patients){
		my $fileout = $dir_prod."/".$patient->name.".cnv.vcf.gz";
		next if -e $fileout;
			my $cmd = qq{perl $script_perl/dragen_cnv_pon.pl -project=$projectName -patient=}.$patient->name;
				push(@$jobs,{name=>$patient->name.".pon" ,cmd=>$cmd,fileout=>$dir_prod."/".$patient->name.".cnv.vcf.gz"});
	}
}
steps_system("PON",$jobs);
}

sub run_sv {
	my ($projects) = @_;
	my $jobs =[];

foreach my $project (@$projects){
	my $projectName = $project->name;
	my $dir_prod2 = $project->getVariationsDir("dragen-sv");
	my $patients = $project->getPatients($patients_name);
	foreach my $patient (@$patients){
		#push(@$jobs,{name=>$patient->name.".sv", cmd=>"sleep 3",out=> $dir_prod2."/".$patient->name.".sv.vcf.gz"});
	 	next if -e $dir_prod2."/".$patient->name.".sv.vcf.gz";
		my $cmd = qq{perl $script_perl/dragen_sv.pl -project=$projectName -patient=}.$patient->name;
		push(@$jobs,{name=>$patient->name.".sv", cmd=>$cmd,out=> $dir_prod2."/".$patient->name.".sv.vcf.gz"});
	}
}
	steps_system("SV",$jobs);
}

sub run_target {
	my ($projects) = @_;
	my $jobs =[];

foreach my $project (@$projects){
	my $projectName = $project->name;
	my $dir_prod = $project->getVariationsDir("dragen-target");
	my $patients = $project->getPatients($patients_name);
	 foreach my $patient (@$patients){
	 	#	push(@$jobs,{name=>$patient->name.".target", cmd=>"sleep 3",out=> $patient->targetGCFile()});
	 	next if -e $patient->targetGCFile();
		my $cmd = qq{perl $script_perl/dragen_target.pl -project=$projectName -patient=}.$patient->name;
		push(@$jobs,{name=>$patient->name.".target", cmd=>$cmd,out=> $patient->targetGCFile()});
	}
}
	steps_system("Target",$jobs);
}


##########################
#STEP RUN SUB 
##########################



sub steps_system {
	my ($name , $jobs) = @_;
	$num_jobs ++;
	my $text = "$num_jobs- $name";
	running_text($text,$num_jobs);
	run_system($jobs,$num_jobs);
	text_system($text,$jobs);
}

##########################
#GENERIC SUB 
##########################

sub run_system {
	my ($jobs,$row) = @_; 
	$row ++;
	return 0 unless @$jobs;
	my $nb_jobs =  scalar(@$jobs);
	my $status = new Term::StatusBar (
                    label => 'jobs Done : ',
                   showTime=>1,
                   subTextAlign =>"center",
                    totalItems => $nb_jobs,  ## Equiv to $status->setItems(10)
                   startRow =>$row,
 );
$status->start();
	
my $exit;
my $failed =0;
my $pending = scalar(@$jobs);
my $running=1;
my $ok=0;
my $line = 4;
$line = ($row -1) + $line if $row > 0;
my $line2 = $line + 2;
foreach my $hcmd (@$jobs){
		my $t1 = time;
		my $myproc = Proc::Simple->new(); 
		print  "\033[".($line2+1).";0H",colored(['bright_magenta on_black'],"running job: ".$hcmd->{name});
		push(@running_jobs,$myproc);
		$myproc->redirect_output ($dir_log."/".$hcmd->{name}.".log", $dir_log."/".$hcmd->{name}.".err");
		push(@running_jobs,$myproc);
		$myproc->start($hcmd->{cmd});
		$running = 1;
		$pending --;
		my $FH = $status->{fh};
		print $FH "\033[$line;0H", (' 'x(5))."total: ".colored(['bright_cyan on_black'],$nb_jobs)." Pending: ".colored(['bright_blue on_black'],$pending)." Running: ".colored(['bright_magenta on_black'],$running)." OK: ".colored(['bright_green on_black'],$ok)."  Failed: ".colored(['bright_red on_black'],$failed)."  ";
		
		while ($myproc->poll()){
			my $ttt = Time::Seconds->new( time -$start_time );
			$status->subText(colored(['bright_yellow on_black'],$ttt->pretty));
			sleep 1;
		}
	my $ttt = Time::Seconds->new( time - $t1 );
	$hcmd->{elapse} = $ttt->pretty;
	$hcmd->{time} = time - $t1;
	
	print  $FH "\033[".$line2.";0H",colored(['bright_blue on_black'],"last job: ".$hcmd->{name}." : ".$hcmd->{elapse});
	
	if ($myproc->exit_status() == 0){
		$ok ++;
		$hcmd->{ok} ++;
	}
	else {
		$failed ++;
		$hcmd->{failed} ++;
	}
		 $status->update();
	}
undef $status;	
system("clear");
return $jobs;	

}



sub text_system {
	my ($title,$jobs) = @_;
	system("clear");
	my $nbo = 0;
	my $sum_time =0;
	my (@error) = grep {exists $_->{failed}} @$jobs;
	my (@ok) = grep {exists $_->{ok}} @$jobs;
	my $sum = 0;
	map {$sum+=$_->{time}} @ok;
	my $mean =0;
	$mean = int($sum/scalar(@ok)) if @ok; 
	my $ttt = Time::Seconds->new($mean);
	my $tt = Time::Seconds->new($sum);
	my @text;
	
	if (@error) {
		foreach my $job (@$jobs){
			print colored(['bright_red on_black'],$job->{name}." ".$job->{elapse})."\n";
			print colored(['bright_magenta on_black'],$job->{cmd})."\n";
			print "log : ".$dir_log."/".$job->{name}.".err";
		}
	die();
	}
	my $new = colored(['bright_cyan on_black'],"$title NONE ")."\n";
	if (@$jobs){
		$new =  colored(['bright_white'],"$title "). colored(['bright_green on_black'],"OK"). colored(['bright_white']," : ".$tt->pretty." (".$ttt->pretty).")\n";
	}
	push(@$sentences,$new);
	print join("",@$sentences)."\n";
}






#################
# RUN CLUSTER
##################
###############
# CLUSTER
###############
sub steps_cluster {
	my ($name,$jobs,$limit) = @_;
	
	$num_jobs ++;
	my $text = "$num_jobs- $name";
	running_text($text,$num_jobs);
	$jobs = run_cluster($jobs,$num_jobs,500);
	text_system($text,$jobs);
}

sub run_cluster { 
my ($commands,$row,$limit_jobs)	= @_;
# init 
 $runnings_jobs = {};
 $cluster_jobs = {};

my $t =0;

#STATUS BAR
my $status = new Term::StatusBar (
                    label => 'jobs Done : ',
                   showTime=>1,
                   subTextAlign =>"center",
                    totalItems => scalar(@$commands),  ## Equiv to $status->setItems(10)
                    startRow => $row,
                   
 );
 $status->{maxCol} = 200;
 $limit_jobs =10000 unless $limit_jobs;
 
 
 
 my $nb_jobs = scalar(@$commands);
while (@$commands) {
	
	last if ($t>=$limit_jobs); 
	my $cmd = shift(@$commands);
	$t++;
	run_cmd($cmd);
	
	
}

$status->subText("Submitting Jobs ...." );

sleep(1);
$|=1;
#system("clear");
$status->start();
my $completed ;
my $failed ;	
my $cancel;
my $start_time = time;
while (keys %$runnings_jobs ) {
	my $jobidline = join(",", keys %$runnings_jobs);	
	
	my @t = `/cm/shared/apps/slurm/16.05.8/bin/sacct -b -n -j $jobidline`;
	chomp(@t);
	
	my $hstatus;
	my $running = 0;
	my $pending =0;
	foreach my $l (@t) {
		my @case = split(" ",$l);
		my $jobid = $case[0];
		next unless exists $cluster_jobs->{$jobid};
		my $exit = pop(@case);
		my $stat = pop(@case);
		$cluster_jobs->{$jobid}->{status} = $stat;
		
		$hstatus->{$stat}->{$jobid} ++;
		if (lc($stat) eq "failed"){
			my $ttt = Time::Seconds->new( time - $cluster_jobs->{$jobid}->{t1} );
			$cluster_jobs->{$jobid}->{elapse} = $ttt->pretty;
			$cluster_jobs->{$jobid}->{time} = time - $cluster_jobs->{$jobid}->{t1};
			$cluster_jobs->{$jobid}->{failed} ++;
			#system("mv slurm-".$jobid."out"." slurm-".$jobid."failed")
		}
		if (lc($stat) eq "completed"){
			$cluster_jobs->{$jobid}->{t1} = time unless exists $cluster_jobs->{$jobid}->{t1};
			$cluster_jobs->{$jobid}->{time} = time - $cluster_jobs->{$jobid}->{t1};
			$cluster_jobs->{$jobid}->{time} = 1 if $cluster_jobs->{$jobid}->{time}  ==0;
			my $ttt = Time::Seconds->new( time - $cluster_jobs->{$jobid}->{t1} );
			$cluster_jobs->{$jobid}->{time} = time - $cluster_jobs->{$jobid}->{t1};
			$cluster_jobs->{$jobid}->{elapse} = $ttt->pretty;
			delete $cluster_jobs->{$jobid}->{failed};
			$cluster_jobs->{$jobid}->{ok} ++;
			delete $runnings_jobs->{$jobid};
#			unlink "slurm-".$jobid.".out";
			$completed->{$jobid} ++;
#			unlink "slurm-".$jobid."out";
			if (@$commands){
				my $cmd = shift(@$commands);
				run_cmd($cmd);
			}
			#$status->update();
		}
		elsif (lc($stat) eq "pending"){
			$pending ++;
		}
		elsif (lc($stat) eq "running"){
			$cluster_jobs->{$jobid}->{t1} = time  unless exists $cluster_jobs->{$jobid}->{t1};
			$running ++;
		}
		elsif (lc($stat) eq "cancel"){
			$cluster_jobs->{$jobid}->{failed} ++;
			delete $runnings_jobs->{$jobid};
			$cancel->{$jobid} ++;
#			unlink "slurm-".$jobid."out";
			$status->update()  ;
		}
		
		else {
			delete $runnings_jobs->{$jobid};
			if (@$commands){
				my $cmd = shift(@$commands);
				run_cmd($cmd);
			}
			$failed->{$jobid} ++;
			$status->update();
		}
		my $failed = scalar(keys %$failed) +scalar(keys %$cancel);
		my $ok = scalar(keys %$completed);
		my $ttt = Time::Seconds->new( time -$start_time );
		$status->subText(colored(['bright_yellow on_black'],$ttt->pretty));
		
		my $FH = $status->{fh};
		my $r = 4;
		$r = ($row -1) + $r if $row > 0;
		print $FH "\033[$r;0H", (' 'x(5))."total: ".colored(['bright_cyan on_black'],$nb_jobs)." Pending: ".colored(['bright_blue on_black'],$pending)." Running: ".colored(['bright_magenta on_black'],$running)." OK: ".colored(['bright_green on_black'],$ok)."  Failed: ".colored(['bright_red on_black'],$failed)."  ";
		
	}
	
	sleep(2) if keys %$runnings_jobs;
	}
	
	my @j = values %$cluster_jobs;
	return \@j;
	die() if keys %$runnings_jobs;
}







	
sub run_cmd {
	my ($cmd) = @_;
	my $out;
	my $cpus = $cmd->{cpus};
	my $script = qq{#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=$cpus
};
	my $run = $script."\n  ".$cmd->{cmd}."\n";
	open(TOTO,">$file_cluster") or die();
	print TOTO $run;
	close TOTO;
	my $n1  = $dir_log."/".$cmd->{name}."-slurm-%A.out";
	my $n2  = $dir_log."/".$cmd->{name}."-slurm-%A.out";
	
	my $jobid = `/cm/shared/apps/slurm/16.05.8/bin/sbatch  --parsable $file_cluster -o $n1 -e $n2`;
	unlink $file_cluster;
	chomp($jobid);
	#warn $file;
	#run3 ["sbatch  "],\$file,\$out ;
	$cmd->{jobid} = $jobid;
	$cluster_jobs->{$jobid}->{log} = $dir_log."/".$cmd->{name}."-slurm";
	$cluster_jobs->{$jobid} =  $cmd;
	#$cluster_jobs->{$jobid}->{cmd} = $cmd->{cmd};
	#$cluster_jobs->{$jobid}->{name} = $cmd->{name};
	#$cmd->{id} = $jobid;
	$cluster_jobs->{$jobid}->{status} = "submitted";
	$runnings_jobs->{$jobid} ++;
	
}








