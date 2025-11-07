#!/usr/bin/perl
use FindBin qw($Bin $RealBin);
use strict;

use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
use GBuffer ;
use Data::Dumper;
use Getopt::Long;
use Bio::DB::HTS;
use Bio::DB::HTS::VCF ;
use List::Util qw(max);
### File Open ###

#!/usr/bin/perl
use FindBin qw($Bin $RealBin);
use strict;

use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
use GBuffer ;
use Data::Dumper;
use Getopt::Long;
use Bio::DB::HTS;
use Bio::DB::HTS::VCF ;
use List::Util qw(max);
### File Open ###



my $project_name;
my $fork;
my $patient_name;
my $version;
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"patients=s" => \$patient_name,
	"version=s" => \$version,
);
# test
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name , -version =>$version);


my $patients =  $project->getPatients();

my $dir_out= $project->getCountingDir("featureCounts");
my $fileout = $dir_out."/".$project_name.".count.genes.txt";
my $fileout2 = $dir_out."/".$project_name.".count.exons.txt";
	
	
	my @bams;
	my @sed_cmd;
	my @sed_cmd2;
	my $nb =0;
	my $align_method;
	my $strand ;
	foreach my $patient (@$patients){
		my $bam = $patient->getBamFile;
		my $name = $patient->name;
		push(@bams,$bam);
		$bam =~ s/\//\\\//g;
		$align_method = $patient->alignmentMethod();
		push(@sed_cmd,qq{sed -i "2s/$bam/$name/" $fileout} );
		push(@sed_cmd2,qq{sed -i "2s/$bam/$name/" $fileout2} );
	}
	my $ppn =16;
	my $gtf = $project->gtf_file();
	$gtf = $project->gtf_file_star() if $align_method eq "star";
	
	my $sed = join(" && ",@sed_cmd);
	my $sed2 = join(" && ",@sed_cmd2);
	my $featureCounts = $project->buffer->software("featureCounts");
	my $cmd = "$featureCounts -T $ppn   -a $gtf --ignoreDup -o $fileout -p -t exon  -s 1 ".join(" ",@bams)." && $sed";
	my $cmd2 = "$featureCounts -T $ppn -a $gtf -f   -t exon  -O --ignoreDup -p  -o $fileout2 -s 1 ".join(" ",@bams)." && $sed2";
	my $result = system($cmd);
	die() unless $result ne 0;
	$result = system($cmd2);
	die() unless $result ne 0;

