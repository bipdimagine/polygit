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
GetOptions(
	'project=s' => \$project_name,
	'patient=s' => \$patient_name,
	'bamin=s' => \$bamin,
	'bamout=s' => \$bamout,
	'fork=s' => \$fork,
);
die("hey man,  no fork ") unless $fork;


 my $project = $buffer->newProject( -name 			=> $project_name );
 my $patient = $project->getPatient($patient_name);
 my $cnvnator_singularity = $buffer->software("cnvnator-singularity");
 my $singularity  = $buffer->software("singularity");
  my $samtools  = $buffer->software("samtools");
 my $dir_out = $project->getCallingPipelineDir("cnvnator");
 my @dd = split("/", $patient->getBamFile());
 my $bamfile = pop @dd;
 my $dirbam = join("/",@dd);
 my $name = $patient->name;
 my $rootfile;
 my $fasta = $project->getGenomeFasta();
 
 warn $dir_out;
# foreach my $chr (@{$project->getChromosomes}){
 #	system("mkdir $dir_out/fasta" );
 #	system("$samtools faidx $fasta ".$chr->fasta_name." >$dir_out/fasta/".$chr->fasta_name.".fa");
 #}
 my $pm = new Parallel::ForkManager($fork);
  foreach my $chr (@{$project->getChromosomes}){
  	my $chrname = $chr->fasta_name;
  	my $singularity_cmd  = qq{$singularity exec --bind $dirbam:/databam --bind $dir_out:/data $cnvnator_singularity /usr/local/bin/cnvnator};;
  	 my $pid = $pm->start and next;
  	
  	system("mkdir $dir_out/fasta" );
 	system("$samtools faidx $fasta ".$chr->fasta_name." >$dir_out/fasta/".$chr->fasta_name.".fa");
  	my $rootfile = $patient->name.".$chrname.root";
 	my $cmd = qq{$singularity_cmd -root /data/$rootfile  -tree /databam/$bamfile  -call 500   -chrom $chrname  -d /data/fasta/};
  	system($cmd);
  		my $cmd2 = qq{$singularity_cmd -root /data/$rootfile  -his 1000     -chrom $chrname  -d /data/fasta/};
  		my $cmd3 = qq{$singularity_cmd -root /data/$rootfile  -stat 1000     -chrom $chrname  -d /data/fasta/};
  		my $cmd4 = qq{$singularity_cmd -root /data/$rootfile  -partition 1000     -chrom $chrname  -d /data/fasta/};
  		my $cmd5 = qq{$singularity_cmd -root /data/$rootfile  -call 1000     -chrom $chrname  -d /data/fasta/ > $name.$chrname.call};
  			system($cmd2."&& ".$cmd3."&& ".$cmd4."&& ".$cmd5);
  	$pm->finish();
  }
  $pm->wait_all_children();
  