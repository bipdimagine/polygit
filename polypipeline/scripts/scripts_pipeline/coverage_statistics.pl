#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";

use Data::Dumper;
use GBuffer;
use Getopt::Long;
use Carp;

use Storable qw(store retrieve freeze);
use PBS::Client;
use Term::ANSIColor;
use pipeline_steps; 

my $buffer_1 = GBuffer->new();

my $projectName;
my $fileName;
my $patientName;
my $method_name;
my $output;
my $no_insert;
my $max_proc=1;
 use Set::IntSpan::Island;
$| =1;
GetOptions(
	'project=s'   => \$projectName,
);


unless ($projectName ){	
	confess("usage :\n $0 -project=projectname  -patient=patient_name -method=type_of_program(i.e maq) (-no_insert=1)" );
}
my $project = $buffer_1->newProject(-name=>$projectName);
my $build = $project->buffer->build();
die ("unknown project ".$projectName ) unless $project;
my $methods = $project->getAlignmentMethods();
my $bam_dir = $project->getAlignmentDir($methods->[0]);
my $coverage_dir = $bam_dir."./coverage_all/";
mkdir($coverage_dir);
my $patients = $project->getPatients();
my $client = PBS::Client->new();
my $text;
my $pname = $project->name();
my $bed = $project->getCaptureBedFile();
#my $bed = "/data-xfs/public-data/$build/capture/agilent/agilent.v50.bed";
#$bed =  qq{/data-xfs/FRED/SEQ/SOFT/Alport/alport_exons.bed};
my $dir = $project->getPipelineDir();
my $file = $dir."/".$project->name().".coverage.log";
system ("date > $file ");
my $pipeline = pipeline_steps->new(log_file=>$file);
foreach my $p (@$patients){
	by_patient($p,$pipeline);
}

print colored ['black ON_BRIGHT_GREEN'],"the log file is here : $file\n";
print color 'reset';
print "\n";
exit(0);
#	warn $p->name();
#	my $file = $p->getBamFile();
#	warn "WARNING can't find $file"  unless -e $file;;
#	next unless -e $file;
#	my $name = $p->name();
#	my $fileout= $coverage_dir."/".$p->name();
#	$text="-I $file ";
#	my $command = "/opt/java/jdk1.6.0_14/bin/java -Xmx10g -jar /bip-d/soft/distrib/GATK/gatk-latest/GenomeAnalysisTK.jar -R /data-xfs/public-data/$build/genome/fasta/all_nochr.fa   -T DepthOfCoverage  -o $fileout --omitDepthOutputAtEachBase  --omitLocusTable -L $bed $text -ct 1 -ct 5 -ct 15 -ct 30 -ct 50"";
#	my $jobs = PBS::Client::Job->new(
#         # Job declaration options
#         # Resources options
#       
#         ppn       => 4,             # process per node
#         nodes => 1,
#       	 script=>$project->name."_coverage",
#         cmd       => [$command],# command to be submitted
#     
# );
#my $commandOK = qq{echo "$pname	$name OK" >> OK.txt};
#	my $jobsOK = PBS::Client::Job->new(
#         # Job declaration options
#         # Resources options
#       
#         ppn       => 1,             # process per node
#         nodes => 1,
#       	 script=>$project->name."_OK",
#         cmd       => [$commandOK],# command to be submitted
#     
# );
#my $commandERROR = qq{echo "$pname	$name ERROR" >> ERROR.txt};
#my $jobsERROR = PBS::Client::Job->new(
#         # Job declaration options
#         # Resources options
#       
#         ppn       => 1,             # process per node
#         nodes => 1,
#       	 script=>$project->name."_ERR",
#         cmd       => [$commandERROR],# command to be submitted
#     
# );
#$jobs->next({ok => $jobsOK,fail => $jobsERROR});
#
#
#$client->qsub($jobs);

 
#}

sub by_patient {
	my ($p,$pipeline) = @_;
	$pipeline->patient($p);
	my $coverage_job;
	my $current_job;
	my $file = $p->getBamFile();
	warn "WARNING can't find $file"  unless -e $file;
		my $toto;
	unless (ischrornot($file)){
		warn $p->name();
		($current_job,$toto) = $pipeline->add_chr(filein=>$file,previous=>undef);
	}
	next unless -e $file;
	($coverage_job,$toto) = $pipeline->coverage(filein=>$file,previous=>$current_job);
	my $client = PBS::Client->new();
	$client->qsub($coverage_job);
}



sub ischrornot {
	my ($file) = @_;
	my @res = `samtools view -H $file`;

	my ($find) = grep {$_=~ /SN:chr/} @res;
	warn $file ;
	return $find;
}

