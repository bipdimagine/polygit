#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages/";
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
my $gatk  = $project->getSoftware('gatk4');
my $reference = $project->genomeFasta();

my $dir_gvcf_out= $project->getCallingPipelineDir("gvcf");
my $dir_vcf_out= $project->getCallingPipelineDir("haplotypecaller4");

print "GVCF \n";
my $patients =  $project->get_list_patients($list_patients);
my @gvcf;	
foreach my $patient (@$patients){
	my $gz = $patient->getGvcfFile("haplotypecaller4");
	confess("no gvcf for ".$patient->name) unless -e $gz;
	#system("$tabix -p vcf $gz") unless -e $gz.".tbi";
	#confess($gz.".tbi") unless -e $gz.".tbi";
	push(@gvcf,$gz); 
}	


	
my $string_variant = join (" --variant ",@gvcf);

my $rregions;

foreach my $chr (@{$project->getChromosomes}){
#my $chr = $project->getChromosome(1);
#next if $chr->name eq "MT";
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
  my $idex =time;
  while( my @tmp = $iter->() ){
  	warn "coucou $n ".scalar(@tmp);
  	$n++;
  	$idex ++;
  	my $t1;
  	#$t1->{regions} =  \@tmp;
  	foreach my $r (@tmp){
	my @t = split("-",$r);
	push(@{$t1->{array_bed}},$r->{bed_line} );
	}
  	$t1->{vcf} =  getTmpFile($dir_vcf_out,"all","vcf");
  	$t1->{vcf1} =  getTmpFile($dir_vcf_out,"all","vcf");
  	$t1->{vcf1_bug} = "/data-beegfs/bug/tmp/1/".$idex.".vcf";
  	$t1->{vcf_bug} = "/data-beegfs/bug/tmp/2/".$idex.".vcf";
  	$t1->{bed} = getTmpFile($dir_vcf_out,"all","bed");
  	$t1->{ok} = getTmpFile($dir_vcf_out,"all","ok");
    push (@{$slice_regions},$t1);
  }

my $pm = new Parallel::ForkManager($fork/2);

my $index =0;
 
  foreach my $slice (@{$slice_regions}){
  	$index ++;
	
	$project->disconnect();
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
   	die($slice->{vcf}) unless -e $slice->{vcf};
   		open (OUT,$slice->{vcf});
					while(<OUT>){
				
						print GLOBAL $_;
					}
					close OUT;
		#unlink 	$slice->{vcf};
   }

   
   close GLOBAL;
exit(0);


sub genotype {
my ($slice,$index) = @_;
my $bed = $slice->{bed};
my $vcf = $slice->{vcf};
my $vcf1 = $slice->{vcf1};

my $gvcf = $vcf1.".db";
my $ok = $slice->{ok};
open (BED,">$bed");
print BED join("\n",@{$slice->{array_bed}});
close BED;

#my $cmd1 = qq{$gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport -R $reference --genomicsdb-workspace-path $gvcf  -V $string_variant   -L $bed  --verbosity ERROR };
##my $cmd3 = qq{$bcftools filter $vcf1   -e "MAX(AD[*:1]+AD[*:0])<7 || MAX(AD[*:1])<5" >$vcf};
##my $cmd3 = qq{/software/distrib/bcftools/bcftools-1.4/bcftools filter $vcf1   -e "MAX(AD[1]+AD[0])<7 || MAX(AD[1])<5" >$vcf};
#my $cmd3 ="";
##my $cmd2 = qq{$gatk   GenotypeGVCFs --java-options "-Xmx8g" -R $reference -O $vcf1 -V $gvcf   -L $bed  --verbosity ERROR && $cmd3 && date > $ok};
#my $cmd2 = qq{ $gatk   GenotypeGVCFs --java-options "-Xmx8g" -R $reference -O $vcf -V gendb://$gvcf   -L $bed  --verbosity ERROR  && date > $ok};

my $cmd1 = qq{$gatk --java-options "-Xmx4g -Xms4g" CombineGVCFs -R $reference -O $gvcf  -V $string_variant   -L $bed  --verbosity ERROR };
#my $cmd3 = qq{$bcftools filter $vcf1   -e "MAX(AD[*:1]+AD[*:0])<7 || MAX(AD[*:1])<5" >$vcf};
my $cmd3 = qq{$bcftools filter $vcf1  -i " N_PASS(FMT/GT = 'alt' & (FMT/DP > 5 | FMT/GQ > 95) )> 0"  >$vcf};
my $cmd2 = qq{ $gatk   GenotypeGVCFs --java-options "-Xmx8g" -R $reference -O $vcf1 -V $gvcf   -L $bed  --verbosity ERROR && $cmd3 && date > $ok };
#my $cmd2 = qq{ $gatk   GenotypeGVCFs --java-options "-Xmx8g" -R $reference -O $vcf -V $gvcf   -L $bed  --verbosity ERROR  && date > $ok};

warn $cmd2;
system($cmd1." && ".$cmd2." ");
unlink $gvcf;
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