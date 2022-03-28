#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use GBuffer;
use Data::Dumper;
use Getopt::Long;
use Carp;
use Logfile::Rotate; 
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
my $filein;
my $dir;
my $file_bed;

my $fork = 1;
my $dp_limit = 7;
my $project_name;
my $patient_name;
my $only_snp;
my $log_file;
my $method;
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"log=s" =>\$log_file,
	"patient=s" => \$patient_name,
	"method=s"=>\$method
);
my $date = `date`;
chomp($date);

if ($log_file){
	system ("echo  ---- start running haplotypecaller $patient_name $date >> $log_file");
	open (STDOUT,">>".$log_file);
}
print "echo  ---- start running UnifiedGenotyper $patient_name $date\n";

my $buffer = GBuffer->new();
my $tabix = $buffer->software("tabix");
my $bgzip = $buffer->software("bgzip");
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatient($patient_name);
my $vcftools = qq{/bip-d/soft/distrib/vcftools/latest/bin/vcftools };
my $bcftools = $buffer->software("bcftools");


my $dir_in= $project->getCallingPipelineDir("unifiedgenotyper")."/".$patient_name."/";
warn $dir_in;
my $vcf_in = $dir_in."$patient_name.final.vcf";
die("$vcf_in not found") unless -e $vcf_in;
# SNP
my $dir_snp= $project->getVariationsDir($method);
warn $dir_snp;

backup($patient,$dir_snp);

my $vcf_snp = $dir_snp."/$patient_name.vcf";

 my $cmd1 = qq{$bcftools view -v snps,mnps $vcf_in  -o $vcf_snp};
system($cmd1);

system ("$bgzip $vcf_snp && $tabix -p vcf -f  $vcf_snp.gz");

warn "$vcf_snp";
if ($only_snp){
print  colored ['black ON_BRIGHT_MAGENTA']," \n ==> $patient_name MOVE ONLY SNP STOP HERE, FOR ION TORRENT  \n"; 
print  color 'reset';
print  "\n"; 
}
# INDEL
my $dir_indel= $project->getIndelsDir($method);
backup($patient,$dir_indel);

my $vcf_indel = $dir_indel."/$patient_name.vcf";

 my $cmd2 =  qq{$bcftools view -V snps,mnps $vcf_in  -o $vcf_indel}; #qq{$vcftools --vcf $vcf_final --keep-only-indels --out $vcf_indel --recode-INFO-all --recode >/dev/null};
 system($cmd2);

system ("$bgzip $vcf_indel && $tabix -p vcf -f  $vcf_indel.gz");








sub backup {
	my ($patient,$dir) = @_;
	my $final_gz = $dir.$patient->name().".vcf.gz";
	my $dirb = $dir."/backup";
	unless (-e $dirb){
		mkdir $dirb unless -e $dirb;
		system("chmod a+w $dirb");
	}
	if (-e $final_gz){
	my $log = new Logfile::Rotate( File => $final_gz,	
	 								Count => 7,
	 								Gzip  => 'lib',
	 								 Flock  => 'no',
	 								 Dir => $dirb
	 								);
		warn "start log";
		$log->rotate();	 								
		}
		unlink $final_gz;
		unlink $final_gz.".tbi";
}


sub move_vcf {
	my ($file_in,$final_out,$type,$dir) = @_;
	my $final_gz = $final_out.".gz";
	my $dir = $dir."/backup";
	mkdir $dir unless -e $dir;
	if (-e $final_gz){
	my $log = new Logfile::Rotate( File => $final_gz,	
	 								Count => 7,
	 								Gzip  => 'lib',
	 								 Flock  => 'no',
	 								 Dir => $dir
	 								);
		warn "start log";
		$log->rotate();	 								
		}
system("cp $file_in $final_out");		
#my $cmd = $cmd_select.qq{ --variant $file_in -o $final_out -selectType $type >/dev/null};
##system($cmd);
`bgzip -f $final_out`;
`$tabix -p vcf $final_gz`;
exit(1) unless -e $final_gz;
exit(1) if -s $final_gz ==0;

}
