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
use Text::Table;

my @running_jobs;
my $cluster_jobs;
my $runnings_jobs;
my $sentences = [];
my $file_cluster = "/tmp/tmp.".time."cmd";
my $script_perl = $Bin."/scripts/";
my $script_pipeline = $Bin."/../scripts/scripts_pipeline/";
my $num_jobs =0;
my $lims;
my $elapsed;

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
my $umi;
my $dude;
my $version ="";
my $dry;

##########
my $cmd_failed = [];
my $cmd_cancel = [];
my $step;
my $force;

GetOptions(
	'project=s' => \$project_name,
	'patients=s' => \$patients_name,
	'umi=s' => \$umi,
	'version=s' => \$version,
	'force=s' => \$force,
	'step=s'=> \$step,
	"dry=s" => \$dry,
	#'low_calling=s' => \$low_calling,
);

unless($step){
	$step = "align,gvcf,vcf,sv,cnv,lmdb";
}

my $hstep ;
map {$hstep->{$_}++} split(",",$step);

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
$patients_name ="" unless $patients_name;
my @apatients_name = split(":",$patients_name);
my $status_jobs; 
my $test_umi;
foreach my $pname (split(",",$project_name)){
	my $buffer = GBuffer->new();
	my $project = $buffer->newProject( -name => $pname , -version =>$version);
	#my $project = $buffer->newProject( -name => $pname );
	$project->isGenome;
	$project->get_only_list_patients($apatients_name[0]);
	 $test_umi=1 if grep{$_->umi} @{$project->getCaptures};
	push(@$projects,$project);
	push(@$buffers,$buffer);
}
 $project_name =~ s/\,/\./g ;
my $dir_log = $projects->[0]->buffer->config->{project_pipeline}->{bds}."/".$project_name.".dragen.".time;
system("mkdir $dir_log && chmod a+rwx $dir_log");

if ($test_umi && !($umi)){
		print colored::stabilo("red","Hey Sylvain, it seems to me that you didn't put the UMI=1 option but your project had UMIs  ", 1)."\n";
		my $choice = prompt("do you  want to use it anyway   (y/n) ? ");
		$umi=1 if ($choice eq "y"); 
	
}

 $|=1;
 
# unless $noprint;
system("clear") unless $dry;
my $start_time = time;
my $jobs =[];
my $pipeline_dragen_steps = ["align","gvcf","vcf","sv","cnv"];
####### Alignement
my $ppd  = patient_pipeline_dragen($projects);
run_command($ppd);
run_move($ppd);
if (exists $hstep->{lmdb}){
run_lmdb_depth_melt($ppd) unless $version =~ /HG38/;
}
run_genotype($projects);
#run_dude($projects) if $dude;
end_report($projects,$ppd);

exit(0);



sub running_text {
	my ($text,$row) = @_;
	print   "\033[".$row.";0H",colored(['bright_white'],"$text :").colored(['bright_cyan on_black']," RUNNING ")."\n";
}



########################################################################################################
#########################---------------------  PIPELINE CLUSTER --------------##########################
########################################################################################################

#### dragen_command
sub patient_pipeline_dragen {
my ($projects) = @_;
my $patients_jobs;
$patients_name = "all" unless $patients_name;
my @apatients_name = split(":",$patients_name);
my $dir_log ;;
foreach my $project (@$projects){
	my $projectName = $project->name;
	 $dir_log = $project->buffer->config->{project_pipeline}->{bds}."/".$project_name.".dragen.".time unless $dir_log;
	$project->isGenome;
	$project->get_only_list_patients($apatients_name[0]);
	foreach my $patient (@{$project->getPatients}){
		if ($version =~ /HG38/){
			my ($m) = grep{$_ eq "bwa" or $_ eq "dragen-align" } @{$patient->alignmentMethods};
			next unless $m;
		}
		$status_jobs->{$patient->name."_".$project->name}->{progress} ="waiting" ;
		my $h;
		$h->{run} = [];
		$h->{run_pipeline} = [];
		my $dir_pipeline = $patient->getDragenDirName("pipeline");
		my $prefix = $patient->name;
		$h->{align_dragen}->{file}  = $patient->getBamFileName("dragen-align");
		
		$h->{name} = $patient->name;
		
	 	$h->{project} = $patient->project->name;
	 	$h->{dir_pipeline} = $patient->getDragenDir("pipeline");
	 	$h->{command_option} = join(",",@{$h->{run}});
		
		#####  
		#####  BAM
		#####  
		$h->{pipeline}->{align}   = $dir_pipeline."/".$patient->name.".bam";
		$h->{pipeline}->{align}   = $dir_pipeline."/".$patient->name.".cram" if $version && $version =~ /HG38/;
		
		$h->{prod}->{align} = $patient->getBamFileName(); 
		$h->{prod}->{align} = $patient->getBamFileName("dragen-align") unless -e $patient->getBamFileName; 
		
		#$h->{prod}->{align_HG38} = $patient->getBamFileName("dragen-align"); 
		$h->{prod}->{file}   = $patient->getBamFileName("dragen-align") if $version;
		$h->{prod}->{align}   = $patient->getCramFileName("dragen-align") if $version && $version =~ /HG38/; 
		#####  
		#####  GVCF
		#####  
		
		 $h->{prod}->{gvcf} = $patient->gvcfFileName("dragen-calling");
		 $h->{pipeline}->{gvcf} = "$dir_pipeline/".$prefix.".hard-filtered.gvcf.gz";
		
		#####  
		#####  VCF
		#####  
		
		 $h->{prod}->{vcf} = $patient->vcfFileName("dragen-calling");
		 $h->{pipeline}->{vcf} = "$dir_pipeline/".$prefix.".hard-filtered.vcf.gz";
		
		#####  
		#####  CNV
		#####  
	
		my $dir = $patient->project->getTargetCountDir();
		$h->{prod}->{cnv}  = $dir."/".$patient->name.".target.counts.gc-corrected.gz";
		$h->{pipeline}->{cnv}  = "$dir_pipeline/".$prefix.".target.counts.gc-corrected.gz";
		
		#####  
		#####  SV
		#####  
		my $dir_prod2 = $project->getVariationsDir("dragen-sv");
	 	$h->{prod}->{sv}  =  $dir_prod2."/".$patient->name.".sv.vcf.gz";
	 	$h->{pipeline}->{sv}  =  $dir_pipeline."/".$prefix.".sv.vcf.gz";;
	 	
	 	#####  
		#####  lmdb
		#####  
	 	$h->{prod}->{lmdb} = $patient->fileNoSqlDepth;
	 	
	 	#####  
		#####  melt
		#####  
		
		$h->{prod}->{melt} = $project->getVariationsDir("melt")."/".$patient->name.".vcf.gz";	 	

	 

		push(@{$patients_jobs},$h);
	}
	}
	return $patients_jobs;
}

my @all_list ;



sub run_command {
	my ($patients_jobs) = @_;
	 @all_list = @$pipeline_dragen_steps;
	foreach my $hp (@$patients_jobs) {
	my $job;
	my $dir_pipeline = $hp->{dir_pipeline};
	my $pname =  $hp->{name}."_".$hp->{project};
	$hp->{run_pipeline} = [];
	
	foreach my $t (@$pipeline_dragen_steps){
		$lims->{$pname}->{$t} = "SKIP"; 
		next unless $hstep->{$t};
		next if (-e $hp->{prod}->{$t} && !($force));
		if ( $force  &&  -e $hp->{pipeline}->{$t} ){
			unlink $hp->{pipeline}->{$t};
		}
		next if -e $hp->{pipeline}->{$t} ;
		
		$lims->{$pname}->{$t} = "PLANNED"; 
		push(@{$hp->{run_pipeline}},$t);
		
	}
	next unless @{$hp->{run_pipeline}};
	
	$job->{name} = $hp->{name}.join("_",@{$hp->{run_pipeline}});
	$job->{step_name} = join(";",@{$hp->{run_pipeline}});
	$job->{patient} = $hp->{name}."_".$hp->{project};
	$job->{cmd} = "perl $script_perl/dragen_command.pl -project=".$hp->{project}." -patient=".$hp->{name} ." -command=".join(",",@{$hp->{run_pipeline}});;
	$job->{cmd} .= " -umi=1 " if $umi;
	$job->{cmd} .= " -version=$version " if $version;
	
	
	$job->{jobs_type_list} = join(",",@{$hp->{run_pipeline}});
	#die();
	my $first_cmd = $hp->{run_pipeline}->[0];
	$job->{out} =  $hp->{$first_cmd}->{pipeline};
	push(@$jobs,$job);
}
	my $text = "$num_jobs-  DRAGEN ALIGN";
	steps_system("Dragen :",$jobs);	
	
}



### MOVE 
sub run_move {
	my ($patients_jobs) = @_;
	my @tt;
	foreach my $a (@$pipeline_dragen_steps){
		push(@all_list,"move_".$a);
		push(@tt,"move_".$a);
	}
	my $jobs =[];
foreach my $hp (@$patients_jobs) {
	my $job;
	my $dir_pipeline = $hp->{dir_pipeline};
	my $nname = $job->{patient} = $hp->{name}."_".$hp->{project};
	next if $status_jobs->{$nname}->{progress} eq "failed";
	my $ppn = 20;
	$job->{name} = $hp->{name}.".move";
	
	$hp->{move_options} = [];
	
	foreach my $t (@$pipeline_dragen_steps){
		$lims->{$nname}->{"move_".$t} = "SKIP"; 
		unless (-e $hp->{prod}->{$t}){
			$lims->{$nname}->{"move_".$t} = "PENDING"; 
			push(@{$hp->{move_options}},$t);
			die("$t => ".$hp->{patient}) unless -e $hp->{pipeline}->{$t};
		}
		
	}
	
	next if scalar (@{$hp->{move_options}}) == 0;
	
	my $option = join(",",@{$hp->{move_options}});
	
	next unless $option;
	$job->{jobs_type_list} = join(",",@tt);
	$job->{cmd} = "perl $script_perl/dragen_move.pl -project=".$hp->{project}." -patient=".$hp->{name} ." -command=".$option;
	#$status_jobs->{patient}=$nname;
	$job->{cmd} .= " -version=$version " if $version;
	$job->{cpus} = $ppn;
	push(@$jobs,$job);
}
	my $text = "$num_jobs- MOVE BAM";
	my $limit = 5;
	steps_cluster("MOVE BAM  ",$jobs,$limit);

	
}



my @failed;
### LMDB 
sub run_lmdb_depth_melt {
	my ($patients_jobs) = @_;
	my $jobs =[];
	push(@all_list,"lmdb_depth");
	
foreach my $hp (@$patients_jobs) {
	my $ppn = 20;
	my $project_name = $hp->{project};
	my $job;
	my $nname  = $hp->{name}."_".$hp->{project};
	my $name = $hp->{name};
	$job->{patient} = $nname;
	$status_jobs->{$nname}->{progress} = "failed" unless -e $hp->{prod}->{align};
	next if $status_jobs->{$nname}->{progress} eq "failed";
	my $fileout = $hp->{prod}->{lmdb};
	$lims->{$nname}->{lmdb_depth} = "SKIP" if -e $fileout; 
	next if -e $fileout;
	my  $cmd = qq{perl $script_pipeline/coverage_genome.pl -patient=$name  -fork=$ppn  -project=$project_name  };
	$cmd .= qq{ -version=$version  } if $version;
	$cmd .= qq{ && perl $script_pipeline/coverage_statistics_genome.pl -patient=$name  -fork=$ppn  -project=$project_name};
	$cmd .= qq{ -version=$version  } if $version;
		$job->{name} =  $name.".lmdb";
		$job->{patient} = $nname;
		$job->{cmd} =$cmd;
		$job->{cpus} = $ppn;
		$job->{jobs_type} ="lmdb_depth";
		push(@$jobs,$job);
	}
unless ($umi){	
push(@all_list,"melt");		
foreach my $hp (@$patients_jobs) {
	my $project_name = $hp->{project};
	my $ppn = 5;
	my $job;
	my $nname  = $hp->{name}."_".$hp->{project};
	my $name = $hp->{name};
	$job->{patient} = $nname;
	$status_jobs->{$nname}->{progress} = "failed" unless -e $hp->{prod}->{align};
	next if $status_jobs->{$nname}->{progress} eq "failed";
	my $fileout = $hp->{prod}->{melt};
	$lims->{$nname}->{melt} = "SKIP" if -e $fileout; 
	next if -e $fileout;
	my $cmd1 = "perl $script_pipeline//melt/melt.pl -project=$project_name  -patient=$name -fork=$ppn ";
	$cmd1 .= qq{ -version=$version  } if $version;
	$job->{name} = $name.".melt";
	$job->{cmd} =$cmd1;
	$job->{cpus} = $ppn;
	$job->{jobs_type} ="melt";
	push(@$jobs,$job);
	
}	
steps_cluster("LMDBDepth+Melt ",$jobs);
return;
}
	
steps_cluster("LMDBDepth(-melt) ",$jobs);	
}

### DUDE  
sub run_dude {
	my ($projects) = @_;
	my $jobs =[];

foreach my $project (@$projects){
	
	my $ppn = 20;
	my $project_name = $project->name;
	my $patients = $project->getPatients($patients_name);
	 	my	 $cmd = qq{perl $script_pipeline/dude/dude.pl -fork=$ppn  -project=$project_name  };
	 	my $final_dir = $project->getVariationsDir("dude");
		my $fileout = $final_dir."/".$project->name.".dude";
		next if -e $fileout;
		my $job;
		$job->{name} = $project->name.".dude";
		$job->{cmd} =$cmd;
		$job->{cpus} = $ppn;
		push(@$jobs,$job);
	}
	steps_cluster("DUDE ",$jobs);
}



########################################################################################################
#########################---------------------  PIPELINE DRAGEN --------------##########################
########################################################################################################



##########################
### GENOTYPE 
##########################


sub run_genotype {
my ($projects,$force) = @_;
	my $jobs =[];
	push(@all_list,"genotype");
foreach my $project (@$projects){
	my $projectName = $project->name;
 my @ps;
 my $patients = $project->getPatients();
 my @ps1;
 foreach my $p (@$patients) {
 	my $nname = $p->name."_".$project_name;
 	$lims->{$nname}->{genotype} = "PLANNED" ;
 	my $dir_out= $project->getVariationsDir("dragen-calling");
 	my $f = $dir_out."/".$p->name.".vcf.gz";
 	$lims->{$nname}->{genotype} = "SKIP" if -e $f; 
 	next if -e $f;
 	next unless -e $p->gvcfFileName("dragen-calling");
 	push(@ps,$p->name);
 	push(@ps1,$p->name."_".$project_name);
 }
 if(@ps){
	my $cmd_genotype = "perl $script_perl/dragen_genotype.pl -project=$projectName -patient=".join(",",@ps);
	$cmd_genotype .= " -version=$version " if $version;
	push(@$jobs,{cmd=>$cmd_genotype,name=>$project->name.".genotype",patient=>join(",",@ps1)});
 }
 
}
 	steps_system("Dragen Genotype ",$jobs);
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
	if ($dry){
		foreach my $hcmd (@$jobs){
			print $hcmd->{cmd}."\n";
		}
		die();
	}
	running_text($text,$num_jobs);
	run_system($jobs,$num_jobs);
	text_system($text,$jobs);
}

##########################
#GENERIC SUB 
################
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
	my @tjs = split(",",$hcmd->{jobs_type_list});
	my @pnames = split(",",$hcmd->{patient}); 
	
	if ($myproc->exit_status() == 0){
		$ok ++;
		foreach my $pname (@pnames){
		map{$elapsed->{$pname}->{$_} = $ttt->pretty} @tjs; 	
		map{$lims->{$pname}->{$_} = "OK"} @tjs; 
		}
		$hcmd->{ok} ++;
	}
	else {
		$failed ++;
		foreach my $pname (@pnames){
			map{$lims->{$pname}->{$_} = "FAILED"} @tjs; 
		}
		push(@$cmd_failed,$hcmd->{cmd});
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
	my $et ="";
	
	if (@error) {
		$et ="-";
		
		foreach my $job (@$jobs){
			$status_jobs->{$job->{patient}}->{progress} = "failed";
			$status_jobs->{$job->{patient}}->{failed} = 1;
			push(@{$status_jobs->{$job->{patient}}->{cmd}},$job->{cmd});
		}
		foreach my $k (keys %{$status_jobs}) {
			next unless $status_jobs->{$k}->{progress} eq "failed";
			$et.= $k.";";
		}
		my $new =  colored(['bright_white'],"$title "). colored(['bright_magenta on_black'],"Error: "). colored(['bright_red'],$et)."\n";
		push(@$sentences,$new);
		print join("\n",@$sentences)."\n";
		return;
		#exit(1);	
	#die();
	}
	
	
	my $new = colored(['bright_cyan on_black'],"$title NONE ");
	
		
	$new =  colored(['bright_white'],"$title "). colored(['bright_green on_black'],"OK"). colored(['bright_white']," : ".$tt->pretty." (".$ttt->pretty).")" if (@$jobs);
	
	push(@$sentences,$new);
	print join("\n",@$sentences)."\n";
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
	$jobs = run_cluster($jobs,$num_jobs,$limit);
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
		my $pname = $cluster_jobs->{$jobid}->{patient};
		my $tj = $cluster_jobs->{$jobid}->{jobs_type};
		if (lc($stat) eq "failed"){
			my $ttt = Time::Seconds->new( time - $cluster_jobs->{$jobid}->{t1} );
			$cluster_jobs->{$jobid}->{elapse} = $ttt->pretty;
			$cluster_jobs->{$jobid}->{time} = time - $cluster_jobs->{$jobid}->{t1};
			$cluster_jobs->{$jobid}->{failed} ++;
			my $pname = $cluster_jobs->{$jobid}->{patient};
			my $tj = $cluster_jobs->{$jobid}->{jobs_type};
			push(@$cmd_failed,$cluster_jobs->{$jobid}->{cmd}->{cmd});
			$lims->{$pname}->{$tj} = "FAILED"; 
			#system("mv slurm-".$jobid."out"." slurm-".$jobid."failed")
		}
		if (lc($stat) eq "completed"){
			$cluster_jobs->{$jobid}->{t1} = time unless exists $cluster_jobs->{$jobid}->{t1};
			$cluster_jobs->{$jobid}->{time} = time - $cluster_jobs->{$jobid}->{t1};
			$cluster_jobs->{$jobid}->{time} = 1 if $cluster_jobs->{$jobid}->{time}  ==0;
			my $ttt = Time::Seconds->new( time - $cluster_jobs->{$jobid}->{t1} );
			$cluster_jobs->{$jobid}->{time} = time - $cluster_jobs->{$jobid}->{t1};
			$cluster_jobs->{$jobid}->{elapse} = $ttt->pretty;
			
			$elapsed->{$pname}->{$tj} = $cluster_jobs->{$jobid}->{elapse};
			$lims->{$pname}->{$tj} = "OK"; 
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
		elsif (lc($stat) =~ /cancel/){
			$cluster_jobs->{$jobid}->{failed} ++;
			push(@$cmd_cancel,$cluster_jobs->{$jobid}->{cmd}->{cmd});
			delete $runnings_jobs->{$jobid};
			$cancel->{$jobid} ++;
			my $tj = $cluster_jobs->{$jobid}->{jobs_type};
			$lims->{$pname}->{$tj} = "CANCEL"; 
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
	die();
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
	#run3 ["sbatch  "],\$file,\$out ;
	$cmd->{jobid} = $jobid;
	$cluster_jobs->{$jobid}->{log} = $dir_log."/".$cmd->{name}."-slurm";
	$cluster_jobs->{$jobid}->{cmd} =  $cmd;
	$cluster_jobs->{$jobid}->{jobs_type} = $cmd->{jobs_type};
	#$cluster_jobs->{$jobid}->{cmd} = $cmd->{cmd};
	#$cluster_jobs->{$jobid}->{name} = $cmd->{name};
	#$cmd->{id} = $jobid;
	$cluster_jobs->{$jobid}->{status} = "submitted";
	$cluster_jobs->{$jobid}->{patient} = $cmd->{patient};
	$runnings_jobs->{$jobid} ++;
	
}

sub end_report {
	my ($projects,$patients_jobs) = @_;
	 
	
	foreach my $project (@$projects){
		my @rows = ();
		my $tb = Text::Table->new( (colored::stabilo("blue", $project->name , 1),  @all_list) ) ; # if ($type == 1);
		foreach my $hp ( grep{$_->{project} eq $project->name} @$patients_jobs) {
		
			my @row;
			my $pname = $hp->{name}."_".$hp->{project};
			push(@row,colored::stabilo("blue",$hp->{name},1));
			foreach my $c (@all_list){
				$lims->{$pname}->{$c} = "-" unless exists $lims->{$pname}->{$c};
				my $text = $lims->{$pname}->{$c};
				my $color = "white";
				$color = "green" if $lims->{$pname}->{$c} eq "OK";
				$color = "red" if lc($lims->{$pname}->{$c}) eq "failed";
				$color = "magenta" if lc($lims->{$pname}->{$c}) eq "-";
				$color = "blue" if lc($lims->{$pname}->{$c}) eq "cancel";
				$text = $elapsed->{$pname}->{$c} if $lims->{$pname}->{$c} eq "OK";
				push(@row,colored::print_color("$color",$text,1));
			#print " $c : ".$lims->{$pname}->{$c}."\t";
			}
			push(@rows,\@row);
		}
		
		$tb->load(@rows);
		print $tb;
		print "\n ------------------------------\n";
		#print "\n";
	}
	if(@$cmd_failed){
	print "\n --------ERROR JOB ----------------\n";
	foreach my $c (@$cmd_failed){
		print colored::print_color("red",$c,1)."\n";
	}
	}
	if(@$cmd_cancel){
	print "\n --------CANCEL JOB ----------------\n";
	foreach my $c (@$cmd_cancel){
		print colored::print_color("blue",$c,1)."\n";
	}
	}
		
}







