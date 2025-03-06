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
my $myproc = Proc::Simple->new(); 
 $myproc->redirect_output ("/dev/null", "/dev/null");
my $fork = 10;
my $project_name;
my $debug;
my $patient_name;
GetOptions(
	'project=s' => \$project_name,
	'debug=s' => \$debug,
	'fork=s' => \$fork,
	#'patient=s' => \$patient_name
);


my $buffer = new GBuffer;
my $project = $buffer->newProject( -name 			=> $project_name );
my $file_done = $project->getCNVDir()."/".$project->name.".done";
unlink $file_done if -e $file_done;
my $dir;
my $dir_tmp= $project->getCallingPipelineDir("annot_sv");

$dir->{manta}->{dir}= $project->getVariationsDir("manta");
$dir->{canvas}->{dir}= $project->getVariationsDir("canvas");
$dir->{wisecondor}->{dir}= $project->getVariationsDir("wisecondor");
my $load = qq{module load tcltk/8.6.9};
my $annotsv = qq{/software/distrib/AnnotSV_2.0/bin/AnnotSV};
my $pm = new Parallel::ForkManager($fork);
my $hjobs;
my $job_id = time;
$pm->run_on_finish(
    	sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    		my $j = $data->{job};
    		delete $hjobs->{$j};
    		
    }
    );	

foreach my $patient (@{$project->getPatients}) {
	 $patient->isGenome();
	 my $listCallers = $patient->callingSVMethods();
}
$project->disconnect();

foreach my $patient (@{$project->getPatients}) {
	#next unless $patient->name eq "AO-LR-EA2";
	next unless $patient->isGenome();
	my $listCallers = $patient->callingSVMethods();
	$job_id ++;
	#$hjobs->{$job_id} ++;
	$hjobs->{$job_id}->{file}=   $dir_tmp."/".$patient->name;
	$project->disconnect();
	 my $pid = $pm->start and next;
	 	$listCallers = $patient->callingSVMethods();
	 	$dir_tmp= $dir_tmp."/".$patient->name;
	 	system("mkdir $dir_tmp && chmod a+rwx $dir_tmp") unless -e $dir_tmp;

		foreach my $type (@$listCallers){
		#next if $type ne "wisecondor";
		my $fileout = $patient->getAnnotSVFileName($type);
		my $file_tmp = $dir_tmp."/".$patient->name.".annotsv.tsv";
		#next if -e $fileout;;
		unlink $fileout if -e $fileout;
		
		my $filein = $patient->getSVFile($type);
		my $filebed = $filein;
		$filebed =~ s/\.gz//;
		my $opt = "";
		if ($filein =~ /\.bed/){
			my $t1 = `zcat $filein | wc -l`;
			chomp($t1);
			if($t1 ==1){
				system("touch $fileout");
				next;
			}
		}
		$opt = "-svtBEDcol 6" if $filein =~ /\.bed/;
		my $opt2 = "-genomeBuild GRCh37";
		$opt2 = "-genomeBuild GRCh38" if $project->genome_version_generic() =~/HG38/;
			my $cmd = qq{ export ANNOTSV=/software/distrib/AnnotSV_2.0 && $load && gunzip -c $filein > $filebed && $annotsv $opt2 -SVinputFile $filebed -outputFile $file_tmp $opt 2>/dev/null && mv $file_tmp $fileout };

		
		system("$cmd");
		unlink $filebed;
		die("problem ") unless -e $fileout;
		#my $dir_out = $dir->{$type}->{dir}."/".
		}
	
		$pm->finish(0,{job=>$job_id});
	
}
$pm->wait_all_children();
confess(Dumper %$hjobs ) if scalar(keys %$hjobs); 
#warn "end ";
exit(0);

