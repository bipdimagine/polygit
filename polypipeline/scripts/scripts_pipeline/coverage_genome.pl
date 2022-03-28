#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use Set::IntSpan::Island;
use Set::IntSpan::Fast::XS;
use Array::IntSpan;
use lib $Bin;
 use GenBoBinaryFile;
use GBuffer;
#use Tie::IntegerArray;
use IPC::Open2;
use Data::Dumper;
use Getopt::Long;
use Carp;
use Scalar::Util qw(looks_like_number);
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use colored;
use Bio::DB::Sam;
use Parallel::ForkManager;
use String::ProgressBar;
use Set::IntSpan::Fast::XS;
#use Bio::DB::HTS;
 use JSON::XS;
use List::Util  qw(sum);
use IO::Handle;
use Fcntl 'SEEK_SET'; 
use IO::File ;
use Array::Diff;
#use File::Binary qw($BIG_ENDIAN $LITTLE_ENDIAN $NATIVE_ENDIAN);
use List::MoreUtils qw{ natatime };

#use DB_File ;

my $filein;
my $dir;
my $file_bed;
my $name;
my $fork = 1;
my $project_name;
my $patient_name;
my $verbose;
my $use_samtools;
my $log_file;
GetOptions(
	"fork=s"   => \$fork,
	"project=s" =>\$project_name,
	"patient=s" =>\$patient_name,
	"verbose=i" =>\$verbose,
);
unless($patient_name){
	$patient_name = $name;
}
unless ($name){
	$name = $patient_name;
}
my $buffer = GBuffer->new();
my $project = $buffer->newProject(-name=>$project_name);

my $patient = $project->getPatientOrControl($patient_name);


my $fbout =  $patient->getNoSqlDepth("c");#GenBoBinaryFile->new(name=>$patient->name,dir=>$dir,mode=>"c");
warn $fbout->no->filename;
$fbout->no->put("toto","titi");
$fbout->close();




my $dir = $project->getAlignmentPipelineDir($patient_name."_depth");#$buffer->config->{project_pipeline}->{tmp};

my $pm = new Parallel::ForkManager($fork);
$project->getChromosomes;
	$project->buffer->dbh_deconnect();
	my @tall = (0) x (50_000);
foreach my $chr (@{$project->getChromosomes}) {
	
	my $pid = $pm->start and next;
	
	$project->buffer->dbh_reconnect();
	my $intspan = $chr->getExtentedGenomicSpan(5000);
	my $regions = $chr->chunk(50_000);
	my $nb;
	my $fb =  GenBoBinaryFile->new(name=>$chr->name,dir=>$dir,mode=>"w");
	foreach my $r (@$regions){
		$nb ++;
		my $rspan =  Set::IntSpan::Fast::XS->new();
		$rspan->add_range($r->{start},$r->{end});
		my $aspan= $rspan->intersection($intspan);
		if ($aspan->is_empty){
			my $l = (($r->{end} -$r->{start})+1);
			if ($l == 500000){
				$fb->putDepth($chr->name,$r->{start},$r->{end},\@tall);
			}
			else {
				my @t = (0) x ($l);
			$fb->putDepth($chr->name,$r->{start},$r->{end},\@t);
			}
			
			next;
			
		}
		warn $chr->name." ".$nb."/". scalar(@$regions) ;#if $nb%100 ==0;
		my $gc =  GenBoCoverageSamtools->new(chromosome=>$chr, patient=>$patient, start=>$r->{start}, end=>$r->{end});
		$fb->putDepth($chr->name,$r->{start},$r->{end},$gc->array);
	}
	$fb->save_index();
	$fb->close();
	$pm->finish(0,{});
	}
$pm->wait_all_children();


my $fbout =  $patient->getNoSqlDepth("c");#GenBoBinaryFile->new(name=>$patient->name,dir=>$dir,mode=>"c");

foreach my $chr (@{$project->getChromosomes}){
	my $fb =  GenBoBinaryFile->new(name=>$chr->name,dir=>$dir,mode=>"r");
	my $list = $fb->no->get("index_".$chr->name);
	push(@{$fbout->{tree_array}},@$list);
	foreach my $line (@$list){
		my $z = $line->[0];
	#	my ($chr,$start,$end,$id,$type) = split(" ",$line);
		#my $v = $fb->getDepth($chr,$start,$end);
	#	print INDEX $line."\n";
		my $v = $fb->no->get($z->{name});	
		$fbout->no->put($z->{name},$v);
	}
	$fbout->no->put("index_".$chr->name,$list);
	#warn Dumper $fbout->no->get("index_".$chr->name);
	#$fb->unlink;
}

$fbout->close();

exit(0);	





