#!/usr/bin/perl
#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use colored;
#use Set::IntSpan;
use GBuffer;
use Data::Dumper;
use Getopt::Long;
use Carp;
use calling_target;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use String::ProgressBar;
use List::Util qw(sum);



my $project_name;
my $input_file;
my $dir;
my $patient_name;
GetOptions(
	'dir=s'   => \$dir,
	'project=s' => \$project_name,
);

my $buffer = GBuffer->new();

my $project = $buffer->newProject( -name => $project_name );
my $patients = $project->getPatients();

foreach my $patient (@$patients){
my $bam = $patient->getBamFileName("bwa",1);
next if -e $bam;

$patient_name = $patient->name;
my $bc =  $patient->barcode();

my $file = `ls $dir/*$bc*.bam`;
chomp($file);
warn "PROBLEM $bc ".$patient->name() unless -e $file;
next unless -e $file;
#system("mv $file.bai $dir/backup2/");
my $output_file = $patient->getBamFileName;#"$dir/final/$patient_name.bam";
print qq{\ntask ( "$output_file" <- "$file" , taskName := "$patient_name",cpus :=10)\{\n } ;


my $cmd = qq{sys  perl $Bin/run_cmd_change_bam.pl -file1=$file -file2=$output_file -bc=$bc -patient=$patient_name };
print $cmd;
print "\n}";
}
die();
