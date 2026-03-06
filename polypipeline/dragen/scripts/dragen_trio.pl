#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/";
use lib "$Bin/../../../GenBo/lib/GenBoDB";
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../../../GenBo/lib/kyoto/";
use lib "$Bin/../../../GenBo/lib/GenBoDB/writeDB";
use lib "$Bin/../../packages";
use Logfile::Rotate;
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moose;
use GBuffer;
use GenBoProject;
use colored; 
use Config::Std;
use Text::Table;
use file_util;
use File::Temp qw/ tempfile tempdir /;; 

use Term::Menus;
 use Proc::Simple;
 use Storable;
use JSON::XS;
use Net::SSH::Perl; 

 
my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);


my $project_name;
my $patient_name;
GetOptions(
	'project=s' => \$project_name,
);

my $buffer = GBuffer->new();
my $project = $buffer->newProjectCache( -name 			=> $project_name );
my $pipeline_dir = $project->project_dragen_pipeline_path_name();
my $bcftools = $buffer->software("bcftools");
my $tabix = $buffer->software("tabix");
foreach my  $f (@{$project->getFamilies}){
	my @ped;
	next unless $f->isTrio;
	my @ped;
	my @gvcf;
	foreach my $m (@{$f->getMembers}){
		push(@ped,$m->pedigreeLine);
		push(@gvcf,"--variant ".$m->getGvcfFile("dragen-calling"));
	}
	my $fdir = $project->project_dragen_pipeline_path_name()."/".$f->name;
	mkdir $fdir;
	my $fname = $f->name;
	my $ped_file = "$fdir/$fname.ped";
	open(PED, ">", $ped_file) or die "Cannot open $ped_file: $!";
	print PED join("\n",@ped)."\n";
	close PED;
	my $ref_dragen = $project->getGenomeIndex("dragen");
	my $ref_fasta = $project->getGenomeFasta();
	my $cmd_dragen = qq{dragen -f -r $ref_dragen --enable-joint-genotyping true --pedigree-file $ped_file --output-directory $fdir  }.join(" ",@gvcf);
	warn qq{$Bin/../run_dragen.pl -cmd=\"$cmd_dragen\"};
	#my $exit = system(qq{$Bin/../run_dragen.pl -cmd=\"$cmd_dragen\"}) ;#unless -e $f1;
	#die if $exit != 0;
	
	foreach my $m (@{$f->getMembers}){
		#bcftools view -s DYSPH_32_M  dragen.hard-filtered.vcf.gz -Oz -o DYSPH_32_M.sv.vcf.gz
		my $vcf = $m->vcfFileName();
		
		my $vcf1  = "$fdir/".$m->name.".1.vcf.gz";
		
		my $bcf = "$bcftools view -s ".$m->name." ".$fdir."/dragen.hard-filtered.vcf.gz -i \'GT!=\"mis\"\' -Oz -o $vcf1 --threads 10 && $tabix -p vcf $vcf1 ";
		warn $bcf;
		my $exit = system($bcf);
		die if $exit != 0;
		if (-e $vcf){
			my $vcf2 = "$fdir/".$m->name.".2.vcf.gz";
			system ("cp $vcf $vcf2   && tabix -p vcf $vcf2 ");
			my $vcf2_norm = "$fdir/".$m->name.".2.norm.vcf.gz";
			my $vcf1_norm  = "$fdir/".$m->name.".1.norm.vcf.gz";
			#bcftools norm -m -any -f reference.fa patient_run1.vcf.gz -Oz -o run1.norm.vcf.gz
			my $c1 = "$bcftools norm -m -any -f $ref_fasta $vcf1 -Oz -o $vcf1_norm --threads 10 && $tabix -p vcf $vcf1_norm";
			warn $c1;
		 	$exit = system($c1);
			die if $exit != 0;
			$c1 = "$bcftools norm -m -any -f $ref_fasta $vcf2 -Oz -o $vcf2_norm --threads 10 && $tabix -p vcf $vcf2_norm";	
			warn $c1;
			 $exit = system($c1);
			die if $exit != 0;
				warn $c1;
			$c1 = "$bcftools concat -a -D $vcf1 $vcf2 -Oz -o $vcf --threads 10 && $tabix -p vcf $vcf ";
			warn $c1;
			$exit = system($c1);
			die if $exit != 0;
			
		}
		else {
			system("cp $vcf1 $vcf && $tabix -p vcf $vcf");
		}
		
	}
	
}


