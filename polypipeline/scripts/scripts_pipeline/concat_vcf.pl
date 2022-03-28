#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
my $filein;
my $dir;
my $file_bed;
my $project_name;
my $chr_name;
my $fork = 1;
use List::MoreUtils qw(firstidx);
use Vcf;

my @tmp_files;

my $dir_out;
my $log_file;
my $norecal;
my $file_out;
GetOptions(
	"log=s" =>\$log_file,
	"project=s" =>\$project_name,
	"fork=s" =>\$fork,
	"norecal=s"  =>\$norecal,
	"file_out=s"  =>\$file_out,
);
if ($log_file){
	
	open (STDOUT,">>".$log_file);
	print "echo  ---- start running concat-vcf.pl $project_name \n";
}

my ($can_launch_recal_snps, $can_launch_recal_indels);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
$norecal=1 if $project->isDiagnostic();
my $genomeFasta_file = $project->genomeFasta();
my $hapmap_file = $project->gatk_hapmap_file();
my $omni_file = $project->gatk_omni_file();
my $dbsnp_file = $project->gatk_dbsnp_file();
my $genomes1k_snps_phase1_file = $project->gatk_genomes1k_snps_phase1_file();
my $genome1k_indels = $project->gatk_genomes1k_indels_file();
my $mills_indels_file = $project->gatk_mills_indels_file();
my $vcfconcat =  $buffer->getSoftware('vcfconcat');
warn "\n\n";

my $tabix = $buffer->software("tabix")." -f";
my $bgzip = $buffer->software("bgzip")." -f";
my $gatk = $project->getSoftware('java')." -jar ".$project->getSoftware('gatk')." -R $genomeFasta_file -l off ";
my $gatk_recal = $gatk."-T VariantRecalibrator   -rf BadCigar  -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum  -tranche 100.0  -tranche 99.9 -tranche 99.0 -tranche 90.0 -nt $fork ";
my $gatk_recal_snp = $gatk_recal." -mode SNP";
$gatk_recal_snp .= " -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap_file" if ($hapmap_file);
$gatk_recal_snp .= " -resource:omni,known=false,training=true,truth=true,prior=12.0 $omni_file" if ($omni_file);
$gatk_recal_snp .= " -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp_file" if ($dbsnp_file);
$gatk_recal_snp .= " -resource:1000G,known=false,training=true,truth=false,prior=10.0 $genomes1k_snps_phase1_file" if ($genomes1k_snps_phase1_file);

my $gatk_recal_indel = $gatk_recal." --maxGaussians 4 -mode INDEL";
if ($mills_indels_file) { $gatk_recal_indel .= " -resource:mills,known=true,training=true,truth=true,prior=12.0 $mills_indels_file"; }
elsif ($genome1k_indels) { $gatk_recal_indel .= " -resource:1000G,known=false,training=true,truth=true,prior=10.0 $genome1k_indels"; }

$can_launch_recal_snps   = 1 if ($hapmap_file or $omni_file or $genomes1k_snps_phase1_file);
$can_launch_recal_indels = 1 if ($mills_indels_file or $genome1k_indels);
$norecal = 1 unless ($can_launch_recal_snps or $can_launch_recal_indels);

my $apply_recal = $gatk." -T ApplyRecalibration --ts_filter_level 99.0  -nt $fork";#--excludeFiltered



my $captures = $project->getCaptures();
my %chr_list_uniq ;
foreach my $capture (@$captures){
	my @list_chr_captured = $project->getCapturedChromosomes();
	@chr_list_uniq{@list_chr_captured} = ();
	}
my @chr_list_uniq = keys %chr_list_uniq;




my $error =0;
my %hpatients;
#map{$hpatients{$_->name} =0} @{$project->get_list_patients($patient_name)};

if ($log_file){
open (STDOUT,">>".$log_file);
}
my $dir_out= $project->getCallingPipelineDir("unifiedgenotyper");

my $recal_file_snp = $dir_out."/recalibrate_SNP.recal";
my $tranche_file_snp = $dir_out."/recalibrate_SNP.tranches";
my $rscript = $dir_out."/r.R";
$gatk_recal_snp .= " -recalFile $recal_file_snp -tranchesFile $tranche_file_snp -rscriptFile $rscript  ";

my $recal_file_indel = $dir_out."/recalibrate_INDEL.recal";
my $tranche_file_indel = $dir_out."/recalibrate_INDEL.tranches";

$gatk_recal_indel .= " -recalFile $recal_file_indel -tranchesFile $tranche_file_indel -rscriptFile $rscript ";

my $chrs = [1..22,'X','Y'];
my @files;
foreach my $chr_name (@$chrs){
	my $name = "chr".$chr_name;
	my $ext = "uni";
	my $fileout  = $dir_out."/".$name.".$ext.vcf";
	
	unless (-e $fileout){
		unless (grep{/$name/} @chr_list_uniq ){
			print  colored ['black ON_BRIGHT_RED']," no file for chromosome $chr_name . NOT OK because in capture"; 
			print  color 'reset';
			confess();
		}
		print  colored ['black ON_BRIGHT_RED']," no file for chromosome $chr_name .OK because not in capture"; 
		print  color 'reset';
		print  "\n";
		next;
		
	}
	if (-s $fileout == 0){
		print  colored ['black ON_BRIGHT_MAGENTA']," no variations for chromosome $chr_name"; 
		print  color 'reset';
		print  "\n";
		next;
	}
	push(@files,$fileout);
}

my  @temp_files;

my $fileout = $dir_out."/".$project->name.".uni.1.vcf";

$fileout = $dir_out."/".$project->name.".uni.vcf" if $norecal;
print $fileout."\n";
open (VCF,">$fileout");
my $cmd = "$vcfconcat ".join(" ",@files);
my @out = `$cmd`;


my @sall = sort par_chr grep {$_ !~ /^#/} @out;
my @s = grep {$_ =~ /^#CHR/} @out;
my $nb_p = scalar(@s) - 8;
print VCF grep {$_ =~ /^#/} @out;
print VCF @sall;
close(VCF);
exit(0) if $norecal;

system("$bgzip  $fileout ;$tabix -p vcf $fileout.gz");
print  colored ['black ON_BRIGHT_BLUE']," END MERGE VCF "; 
print  color 'reset';
print  "\n";
push(@tmp_files,$fileout."");
my $cmd1 = qq{$gatk_recal_snp -input $fileout.gz};
$cmd1 .= " --maxGaussians 4 " if ($nb_p < 8); 
warn $cmd1;
system($cmd1);
print  colored ['black ON_BRIGHT_BLUE']," END RECALIBRATION SNP"; 
print  color 'reset';
print  "\n";
my $cmd2 = qq{$gatk_recal_indel -input $fileout.gz};
warn $cmd2;
system($cmd2);
print  colored ['black ON_BRIGHT_BLUE']," END RECALIBRATION INDEL"; 
print  color 'reset';
print  "\n";
system ("$bgzip $recal_file_indel ;$tabix -p vcf $recal_file_indel.gz");
system ("$bgzip $recal_file_snp ;$tabix -p vcf $recal_file_snp.gz");
push(@tmp_files,$recal_file_indel."");
push(@tmp_files,$recal_file_snp."");
push(@tmp_files,$tranche_file_snp."");
push(@tmp_files,$tranche_file_indel."");
push(@tmp_files,$rscript."");
my $fileout2= $dir_out."/".$project->name.".uni.2.vcf";

my $cmd3 = qq{$apply_recal -input $fileout.gz -mode SNP -recalFile $recal_file_snp.gz -tranchesFile $tranche_file_snp -o $fileout2};

system($cmd3);
print  colored ['black ON_BRIGHT_BLUE']," END APPLY RECALIBRATION SNP"; 
print  color 'reset';
print  "\n";

unlink $fileout.".gz";
unlink $fileout.".gz.tbi";
system("$bgzip $fileout2 ;$tabix -p vcf $fileout2.gz");
 
my $fileout3= $dir_out."/".$project->name.".uni.vcf";
unlink $fileout3 if -e $fileout3;
my $cmd4 = qq{$apply_recal -input $fileout2.gz -mode INDEL -recalFile $recal_file_indel.gz -tranchesFile $tranche_file_indel -o $fileout3 };
system($cmd4);
print  colored ['black ON_BRIGHT_BLUE']," END APPLY RECALIBRATION INDEL"; 
print  color 'reset';
print  "\n";


unless (-e  $fileout3){
	print  colored ['black ON_BRIGHT_RED']," PROBLEM ON RECALIBRATION "; 
print  color 'reset';
print  "\n";
die();
}
push(@tmp_files,$fileout2);


foreach my $f (@tmp_files){
	if (-e $f){
		unlink $f;
	}
	my $f1 ="$f.idx";
	if (-e $f1){
		unlink $f1;
	}
	$f .=".gz";
	if (-e $f){
		unlink $f;
	}
	$f .=".gz.tbi";
	if (-e $f){
		unlink $f;
	}
}
	print  colored ['black ON_BRIGHT_GREEN']," END RECALIBRATION PROCESS "; 
print  color 'reset';
print  "\n";
sub par_chr ($$) {
	my ($left,$right) = @_;
	my @l = split(" ",$left);
	my @r = split(" ",$right);
	return give_chr_index($l[0]) <=> give_chr_index($r[0]) ||  $l[1] <=> $r[1];
	 
	
}

sub give_chr_index {
	my ($c) = @_;
	$c =~ s/chr//;
	if ($c eq  'X'){
		return 23;
	}
	if ($c eq 'Y'){
		return 24;
	}
	if ($c eq "MT" || $c eq "M"){
		return 25;
	}
	return $c;
}

