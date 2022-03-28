#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use Digest::MD5::File qw (file_md5_hex);

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
use Path::Tiny;

use Set::IntSpan::Island;
use Set::IntSpan;
#use Tree::R;
#use Bio::DB::HTS;
use Bio::DB::Sam;
use Set::IntSpan;
use Path::Tiny;
require "$Bin/lib/interval_util.pm"; 
 my $project_name;
 my $fork;
 my $fileout;
 my $bam;
 my $patient_name;
 my $print;
 my $force ;

GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"patient=s"=>\$patient_name,
	"force=s"=>\$force,
		
);
$fork =1 unless $fork;
my $bug;
my $buffer = GBuffer->new();
my  $project = $buffer->newProject( -name => $project_name );
my $patients = $project->get_only_list_patients($patient_name);
foreach my $patient (@$patients){
	my $bam = $patient->getBamFile();
}
my $nb =0;
 my $pm1 = new Parallel::ForkManager($fork);
 
  $pm1->run_on_finish(
    sub { 
    	my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    		unless (defined($data) or $exit_code>0){
				print qq|No message received from child process $exit_code $pid!\n|; 
				$bug =1;
				#die();
				return;
			}
			my $intspan = $data->{intspan_type};
			my $chr = $data->{chr};
			my $lmdb_file =  $data->{patient}.".ok.callable";
			my $lmdb = GenBoNoSqlLmdb->new(dir=>$project->getCoverageCallable(),mode=>"w",name=>$lmdb_file,is_compress=>1);
			
			  foreach my $k (keys %{$intspan}){
  						$lmdb->put($chr."_".$k,$intspan->{$k});
 				}
 				$lmdb->close();
    	}
    	);
foreach my $patient (@$patients){ 	
	$nb++;
	my $bam = $patient->getBamFile();
	
	
	
	my $lmdb_file =  $patient->name.".callable";
	unlink $project->getCoverageCallable()."/$lmdb_file" if -e  $project->getCoverageCallable()."/$lmdb_file";
	$lmdb_file =  $patient->name.".ok.callable";
	if (-e  $project->getCoverageCallable()."/$lmdb_file"){
		warn $project->getCoverageCallable()."/$lmdb_file";
		my $lmdb = GenBoNoSqlLmdb->new(dir=>$project->getCoverageCallable(),mode=>"w",name=>$lmdb_file,is_compress=>1);
		 my $res = $lmdb->get("timestamp");
		 $lmdb->close();
		 next if $res;
	}
 	my $gvcf = $patient->getGvcfFile();
 		
 
	my $md5 = file_md5_hex($bam.".bai");
	my $md5 = path($bam.".bai")->digest;
	
	my $lmdb;
	$lmdb = GenBoNoSqlLmdb->new(dir=>$project->getCoverageCallable(),mode=>"c",name=>$lmdb_file,is_compress=>1);
	$lmdb->put("start_md5",$md5);
 
 
 
    	warn $nb."/".scalar(@$patients);
 foreach my $chr (@{$project->getChromosomes}){
 	
 	 	my $pid = $pm1->start and next;
 	 	
 		my $data;
 		$data->{chr} = $chr->ucsc_name;
 		$data->{patient} = $patient->name;
		 $data->{intspan_type}  = interval_util::callable_patient_gvcf($project_name,$patient_name,$force,$chr->ucsc_name);
		 $pm1->finish(0,$data);
 }


}
 $pm1->wait_all_children();
  die() if $bug;
  
 foreach my $patient (@$patients){ 	
 	my $lmdb_file =  $patient->name.".ok.callable";
 my $lmdb = GenBoNoSqlLmdb->new(dir=>$project->getCoverageCallable(),mode=>"w",name=>$lmdb_file,is_compress=>1);
 $lmdb->put("timestamp",time);
  $lmdb->put("ok","ok");
  $lmdb->close;

 }

