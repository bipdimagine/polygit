#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
#use Set::IntSpan;
use GBuffer; 
use Data::Dumper;
use Getopt::Long;
use Carp;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use String::ProgressBar;
use List::Util qw(sum);
use Proc::Simple;



my $project_name;
 my $fork;
 my $callable_intspan_file;
 my $patient_name;
 #my $low_calling;
 my $method;
my $chemistry;

GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"chemistry=s" => \$chemistry,
);

my $buffer = GBuffer->new();

my $project = $buffer->newProject( -name => $project_name );

my $cell  = $project->getSoftware("cellranger");
my $ref_root = $buffer->config->{public_data}->{$project->getVersion()};
warn $ref_root;
my $index = $ref_root."/".$buffer->index("cellranger");
my @seqDir = $project->getAlignmentPipelineDir("cellranger_star");
count($project->getPatients,$index);
sub count{
	my ($patients,$index) = @_;
	foreach my $p (@$patients){
		my $patientName = $p->name; 
		warn $p->name();
		my $seqDir=$p->getSequencesDirectory();
		my $outputDir = $project->getCellRangerDir()."/";
		my $dir_tmp = $project->getCallingPipelineDir("cell_ranger");
		chop($seqDir)  unless -e $seqDir;
		$seqDir.="_SC.saved"  unless -e $seqDir;
		my $myproc = Proc::Simple->new();
		my $chem_string;
		$chem_string = "--chemistry=$chemistry" if $chemistry;
		my $count_command = "cd $dir_tmp;$cell count --id=$patientName --transcriptome=$index --fastqs=$seqDir --sample=$patientName  $chem_string --jobmode=slurm && mv $dir_tmp/$patientName $outputDir";
		print $count_command."\n";
		$myproc->start($count_command); 
	}

}