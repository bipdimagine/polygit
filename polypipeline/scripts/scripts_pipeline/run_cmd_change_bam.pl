
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
my $file2;
my $file1;
my $patient_name;
my $bc;
GetOptions(
	'file1=s'   => \$file1,
	'file2=s'   => \$file2,
	'patient=s' => \$patient_name,
	'bc=s' => \$bc,
);


my $cmd = qq{  samtools view -h $file1  | perl -lane '\$a =\$_;\$a=~s/chrMT/chrM/g; \$a=~s/:$bc/:$patient_name/;print \$a ' | samtools view -bS - > $file2 && samtools index $file2 };

system($cmd);