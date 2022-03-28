#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use  File::Temp;
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use colored;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
 use IPC::Open2;
 use List::MoreUtils qw(natatime);
 
  
my $project_name;
my $final_vcf;
my $log_file;
my $list_patients;
my $fork;

GetOptions(
	'project=s'   => \$project_name,
	"log=s" =>\$log_file,
	"vcf=s" => \$final_vcf,
	"patient=s"=>\$list_patients,
	"fork=s" =>\$fork,
);
my $date = `date`;
chomp($date);

my $date = `date`;
chomp($date);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );

my $vcffirstheader = $buffer->software("vcffirstheader");
my $vcfstreamsort = $buffer->software("vcfstreamsort");
my $vcfuniq  = $buffer->software("vcfuniq");
my $bcftools = $buffer->software("bcftools");
my $bgzip = $buffer->software("bgzip");
my $tabix = $buffer->software("tabix");
my $javac = $project->getSoftware('java');
my $gatk  = $project->getSoftware('gatk');
my $reference = $project->genomeFasta();

my $dir_gvcf_out= $project->getCallingPipelineDir("gvcf");
my $dir_vcf_out= $project->getCallingPipelineDir("haplotypecaller");

print "GVCF \n";
my $patients =  $project->get_list_patients($list_patients);
my @gvcf;	
foreach my $patient (@$patients){
	my $gz = $patient->getGvcfFile();
	confess("no gvcf for ".$patient->name) unless -e $gz;
	#system("$tabix -p vcf $gz") unless -e $gz.".tbi";
	#confess($gz.".tbi") unless -e $gz.".tbi";
	push(@gvcf,$gz); 
}	
	
my $string_variant = join (" --variant ",@gvcf);

my $rregions;

foreach my $chr (@{$project->getChromosomes}){
#my $chr = $project->getChromosome(1);
#next if $chr->name ne "21";
my $intspan = $chr->getIntSpanCaptureForCalling(250);
my @regions = split(",",$intspan->as_string());
foreach my $t (@regions){
	$t =~s/-/\t/;
	my $h;
	$h->{bed_line} = $chr->fasta_name."\t".$t;
	push(@$rregions,$h);
}
}

warn scalar (@$rregions);
my $ns = int(scalar (@$rregions) / ($fork*4)) +1;
 my $slice_regions;
  my $iter = natatime $ns, @$rregions;
  my $n =0;
  while( my @tmp = $iter->() ){
  	warn "coucou $n ".scalar(@tmp);
  	$n++;
  	my $t1;
  	#$t1->{regions} =  \@tmp;
  	foreach my $r (@tmp){
	my @t = split("-",$r);
	push(@{$t1->{array_bed}},$r->{bed_line} );
	}
  	$t1->{vcf} =  getTmpFile($dir_vcf_out,"all","vcf");
  	$t1->{vcf1} =  getTmpFile($dir_vcf_out,"all","vcf");
  	
  	$t1->{bed} = getTmpFile($dir_vcf_out,"all","bed");
  	$t1->{ok} = getTmpFile($dir_vcf_out,"all","ok");
    push (@{$slice_regions},$t1);
  }

my $pm = new Parallel::ForkManager($fork/2);

my $index =0;
 
  foreach my $slice (@{$slice_regions}){
  	$index ++;

  	 my $pid = $pm->start and next;
  		genotype($slice);
  		$pm->finish();
}
$pm->wait_all_children();


warn "end gatk";
	my $vcffirstheader = $buffer->software("vcffirstheader");
	my $vcfstreamsort = $buffer->software("vcfstreamsort");
	my $vcfuniq  = $buffer->software("vcfuniq");

	open(GLOBAL, " | $vcffirstheader | $vcfstreamsort | $vcfuniq > $final_vcf");
   foreach my $slice (@{$slice_regions}){
   	die() unless -e $slice->{vcf};
   	
  
   		open (OUT,$slice->{vcf});
					while(<OUT>){
				
						print GLOBAL $_;
					}
					close OUT;
		
		unlink 	$slice->{vcf};
   }

   
   close GLOBAL;
exit(0);


sub genotype {
my ($slice) = @_;
my $bed = $slice->{bed};
my $vcf = $slice->{vcf};
my $vcf1 = $slice->{vcf1};
my $ok = $slice->{ok};
open (BED,">$bed");
print BED join("\n",@{$slice->{array_bed}});
close BED;
my $dir1 = "/data-beegfs/bug/tmp/";

my $cmd3 = qq{$bcftools filter $vcf1   -e "MAX(AD[*:1]+AD[*:0])<7 || MAX(AD[*:1])<5" >$vcf};
my $cmd2 = qq{$javac  -Xmx48g -jar $gatk -T GenotypeGVCFs -R $reference --out $vcf1 --variant $string_variant   -L $bed -nt 2 -l off && $cmd3 && date > $ok};
warn $cmd2;
system($cmd2);
unlink $vcf1;
unlink $bed;
return;
}




sub getTmpFile {
	my ($dir,$chr_name,$ext) = @_;
	die() unless $ext;
	die() unless $chr_name;
	#confess($chr_name) unless $chr_name !~/[0,1,2,3,4,5,6,7,8,9,10,11,12,13,15,15,16,17,18,19,20,21,22,X,Y,M]/;
	$dir .="/$chr_name";
	system ("mkdir -p $dir") unless -d $dir;
	 my $file_cat_tmp =  File::Temp->new( TEMPLATE =>"TMP.XXXXXXX",
                        DIR => $dir,
                        UNLINK=>1,
                        SUFFIX => ".$ext");
              $file_cat_tmp->unlink_on_destroy( 1 );         
  return $file_cat_tmp->filename();
}