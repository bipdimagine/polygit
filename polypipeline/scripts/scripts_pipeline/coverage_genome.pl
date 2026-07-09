#!/usr/bin/env perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use lib $Bin;
 use GenBoBinaryFile;
use GBuffer;
use Getopt::Long;
use Carp;
use Term::ANSIColor;
use colored;
use Parallel::ForkManager;
#use Bio::DB::HTS;
use List::Util  qw(sum);

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
my $version;
my $align;
my $force;
GetOptions(
	"fork=s"   => \$fork,
	"project=s" =>\$project_name,
	"patient=s" =>\$patient_name,
	"version=s" =>\$version,
	"align=s" =>\$align,
	"verbose=i" =>\$verbose,
	"force=s" => \$force,
);
unless($patient_name){
	$patient_name = $name;
}
unless ($name){
	$name = $patient_name;
}
my $buffer = GBuffer->new();
my $project = $buffer->newProject(-name=>$project_name, -version =>$version);

my $patient = $project->getPatientOrControl($patient_name);
if($align){
 $patient->{alignmentMethods} = [$align];
}


my $dir_prod = $patient->NoSqlDepthDir."/".$patient->name . ".depth.lmdb";

if (-e $dir_prod){
	unless ($force){
		warn "already run add -force=1";
		exit(0);
	}
	else {
		if (-d $dir_prod."/1/"){
			system("rm $dir_prod/*/*");
		}
		else {
			system("rm $dir_prod/*");
		}
	}
	
}

my $dir = $project->getAlignmentPipelineDir($patient_name."_depth".time);
my $process;
my $pm = new Parallel::ForkManager($fork);
$pm->run_on_finish(
	sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hRes) = @_;
		unless (defined($hRes) or $exit_code > 0) {
			print qq|No message received from child process $exit_code $pid!\n|;
			die();
			return;
		}
#		warn "@@@@@ ".$hRes->{end}." ".$hRes->{ttime};
		delete $process->{$hRes->{end}};
#	
	}
);

$project->getChromosomes;
$project->preload();
	#$project->buffer->dbh_deconnect();
	$project->disconnect();
	
	my @tall = (0) x (50_000);
	my $id = time;

	
	
foreach my $chr (@{$project->getChromosomes}) {
	$id ++;
	$process->{$id} = 1;
	my $pid = $pm->start and next;
	
	my $intspan = $chr->getExtentedGenomicSpan(5000);
#	if ($chr->name eq "Y"){
#		my $gene = $project->newGene("SRY");
#		warn $gene->start." ".$gene->end;
#	}
	my $regions = $chr->chunk(50_000);
	my $nb;
	my $fb =  GenBoBinaryFile->new(name=>$chr->name,dir=>$dir,mode=>"w");
	foreach my $r (@$regions){
		#warn Dumper $r if $chr->name eq "Y";
		$nb ++;
		my $rspan =  Set::IntSpan::Fast->new();
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
		my $gc =  GenBoCoverageSamtools->new(chromosome=>$chr, patient=>$patient, start=>$r->{start}, end=>$r->{end});
		$fb->putDepth($chr->name,$r->{start},$r->{end},$gc->array);
	}
	$fb->save_index();
	$fb->close();
	$pm->finish(0,{end=>$id});
	}
$pm->wait_all_children();

if (keys %$process){
	warn Dumper $process;
	die();
}
my $fbout =  $patient->getNoSqlDepth("c");#GenBoBinaryFile->new(name=>$patient->name,dir=>$dir,mode=>"c");
warn $fbout->no->filename;
$fbout->no->put("toto","titi");
$fbout->close();
warn "END !!!!";
$fbout =  $patient->getNoSqlDepth("c");#GenBoBinaryFile->new(name=>$patient->name,dir=>$dir,mode=>"c");

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

unless ($project->isCapture){
	system("$Bin/coverage_statistics_genome.pl -project=$project_name -patient=$patient_name -fork=$fork");
}
exit(0);	





