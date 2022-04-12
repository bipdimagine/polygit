#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
#use Set::IntSpan;
use GenBoNoSqlLmdb;
use Carp;
use strict;
use Set::IntSpan::Fast::XS ;
use Data::Dumper;
use GBuffer;
use Getopt::Long;
use Carp;
 use JSON::XS;
 use List::MoreUtils qw(natatime uniq);
 use Tabix;
 
 my $buffer = new GBuffer;
my $project_name= "NGS2017_1534";
my $fork;
my $fastq1;
my $fastq2;
my $dirin;
my $patient_name;
my $bamin;
my $bamout;
my $version;
GetOptions(
	'project=s' => \$project_name,
	'patient=s' => \$patient_name,
	'bamin=s' => \$bamin,
	'bamout=s' => \$bamout,
	'fork=s' => \$fork,
	'version=s' =>\$version,
);
die("hey man,  no fork ") unless $fork;


 my $project = $buffer->newProject( -name 			=> $project_name,version=> $version );
 my $patient = $project->getPatient($patient_name);
  my $name = $patient->name;
 my $cnvnator_singularity = $buffer->software("cnvnator-singularity");
 my $singularity  = $buffer->software("singularity");
  my $samtools  = $buffer->software("samtools");
  my $tabix  = $buffer->software("tabix");
  my $bgzip  = $buffer->software("bgzip");
 my $dir_out = $project->getCallingPipelineDir("cnvnator")."/$name";
 system("mkdir  $dir_out") unless -e $dir_out
; my @dd = split("/", $patient->getBamFile());
 my $bamfile = pop @dd;
 my $dirbam = join("/",@dd);

 my $rootfile;
 my $fasta = $project->getGenomeFasta();
 
 warn $dir_out;
# foreach my $chr (@{$project->getChromosomes}){
 #	system("mkdir $dir_out/fasta" );
 #	system("$samtools faidx $fasta ".$chr->fasta_name." >$dir_out/fasta/".$chr->fasta_name.".fa");
 #}
 my $pm = new Parallel::ForkManager($fork);
   	my $singularity_cmd  = qq{$singularity exec --bind $dirbam:/databam --bind $dir_out:/data $cnvnator_singularity };
   	system("mkdir $dir_out/fasta" ) unless -e "$dir_out/fasta";
   	
  foreach my $chr (@{$project->getChromosomes}){
  	my $chrname = $chr->fasta_name;

  	 my $pid = $pm->start and next;
  	
  
 	system("$samtools faidx $fasta ".$chr->fasta_name." >$dir_out/fasta/".$chr->fasta_name.".fa");
  	my $rootfile = $patient->name.".$chrname.root";
 	my $cmd = qq{$singularity_cmd /usr/local/bin/cnvnator -root /data/$rootfile  -tree /databam/$bamfile  -call 500   -chrom $chrname  -d /data/fasta/};
  	system($cmd) unless -e $dir_out."/".$rootfile;
  		my $cmd2 = qq{$singularity_cmd /usr/local/bin/cnvnator -root /data/$rootfile  -his 1000     -chrom $chrname  -d /data/fasta/};
  		my $cmd3 = qq{$singularity_cmd /usr/local/bin/cnvnator -root /data/$rootfile  -stat 1000     -chrom $chrname  -d /data/fasta/};
  		my $cmd4 = qq{$singularity_cmd /usr/local/bin/cnvnator -root /data/$rootfile  -partition 1000     -chrom $chrname  -d /data/fasta/};
  		my $sh = "$name.$chrname.sh";
  		open("TEST",">$dir_out/$sh");
  		print TEST qq{/usr/local/bin/cnvnator -root /data/$rootfile  -call 1000     -chrom $chrname  -d /data/fasta/ > /data/$name.$chrname.call\n};
  		close(TEST);
  		my $cmd5 = qq{$singularity_cmd /bin/bash /data/$sh};
  	#	my $cmd5 = qq{$singularity_cmd /usr/local/bin/cnvnator "-root /data/$rootfile  -call 1000     -chrom $chrname  -d /data/fasta/ > /data/$name.$chrname.call"};
  		system($cmd2."&& ".$cmd3."&& ".$cmd4."&& ".$cmd5);
  		unlink "$dir_out/fasta/".$chr->fasta_name.".fa";
  	$pm->finish();
  }
  $pm->wait_all_children();
  my $all_call = $dir_out."/$name.call";
  unlink $all_call if -e $all_call;
    foreach my $chr (@{$project->getChromosomes}){
    	my $chrname = $chr->fasta_name;
    	system ("cat $dir_out/$name.$chrname.call >> $all_call");
    }
    my $cnvnator2VCF = $buffer->software("cnvnator2VCF");
    my $final_vcf = $project->getVariationsDir("cnvnator")."/".$name.".vcf";
    system("$cnvnator2VCF -prefix $name  -reference HG19  $all_call >$final_vcf");
    system("$bgzip $final_vcf && $tabix -p vcf $final_vcf.gz");
    
  
  