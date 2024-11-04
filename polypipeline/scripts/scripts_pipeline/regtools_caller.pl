#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use GBuffer;
use Data::Dumper;
use Getopt::Long;
use Carp;
use Bio::DB::Sam;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use List::Util qw(sum);
my $filein;
my $dir;
my $file_bed;
 
my $dp_limit = 5;
my $al = 3;
my $project_name;
my $fork;
$| =1;
my $log_file;
my $vcf_final;
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"log=s" =>\$log_file,
	"out=s"=>\$vcf_final
);

$fork = 2 unless $fork; 

my $date = `date`;
chomp($date);
system ("echo  ---- start running regtools $project_name $date");

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $chrs = $project->getChromosomes();
my $reference = $project->getGenomeFasta;
$project->getPatients();

system("rm $log_file") if ($log_file and -e $log_file);
foreach my $patient (@{$project->getPatients()}) {
	my $path_out = $project->getJunctionsDir('regtools').'/';
	if (not -d $path_out) {
		system("mkdir $path_out");
		system("chmod 777 $path_out ");
	}
	$path_out .= $patient->alignmentMethod().'/';
	if (not -d $path_out) {
		system("mkdir $path_out");
		system("chmod 777 $path_out ");
	}
	$path_out .= $patient->name().'/';
	if (not -d $path_out) {
		system("mkdir $path_out");
		system("chmod 777 $path_out ");
	}
}

my $this_fork = $fork / 2;
$this_fork = 1 if $this_fork < 2;
my $pm = new Parallel::ForkManager($this_fork);

foreach my $patient (@{$project->getPatients()}) {
	$pm->start and next;
	my $patient_name = $patient->name();
	print "echo  ---- start running regtools $patient_name $date\n";
	
	my $path_out = $project->getJunctionsDir('regtools').'/'.$patient->alignmentMethod().'/';
	my $final_file = $path_out.'/'.$patient->name().'.tsv';
	$path_out .= $patient->name().'/';
	
	my $out1 = "$path_out/extracted_$patient_name.bed";
	system("rm $out1") if -e $out1;
	my $cmd = "singularity run --writable-tmpfs";
	$cmd .= " --pwd $path_out -B /data-isilon:/data-isilon -B /data-beegfs:/data-beegfs";
	$cmd .= " /data-isilon/bipd-src/mhamici/Regtools/regtools_v_1_0_0.sif regtools junctions extract -o $out1";
	$cmd .= " -s FR ".$patient->getBamFiles->[0];
#	print "\necho  ---- CMD1\n$cmd\n";
	system($cmd);
	
	my $out2 = "$path_out/annotated_$patient_name.bed";
	system("rm $out2") if -e $out2;
	my $cmd2 = "singularity run --writable-tmpfs";
	$cmd2 .= " --pwd $path_out -B /data-isilon:/data-isilon -B /data-beegfs:/data-beegfs";
	$cmd2 .= " /data-isilon/bipd-src/mhamici/Regtools/regtools_v_1_0_0.sif regtools junctions annotate -o $out2";
	$cmd2 .= " $path_out/extracted_$patient_name.bed";
	$cmd2 .= " ".$project->genomeFasta();
	$cmd2 .= " ".$project->gtf_file();
#	print "\necho  ---- CMD2\n$cmd2\n";
	system($cmd2);
	
	my $cmd3 = "/bin/python $Bin/regtools_annotate_merge_polysplice.py -p $path_out -o $final_file";
	system($cmd3);
	
	my $cmd4 = "sort -k2,2 -k3,3n $final_file | bgzip > $final_file.gz";
	system($cmd4);
	
	my $cmd5 = "tabix -s 2 -b 3 $final_file.gz";
	system($cmd5);
	
	$pm->finish();
}
$pm->wait_all_children();

system ("echo  ---- end running regtools $project_name $date");

if ($log_file) {
	open (LOG , '>'.$log_file);
	my $date = `date`;
	chomp($date);
	print LOG "OK - end running regtools $project_name $date";
	close(LOG);
}

exit(0);




