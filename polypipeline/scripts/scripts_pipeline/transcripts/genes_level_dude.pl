#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use Parallel::ForkManager;
use Storable qw(store retrieve freeze thaw);
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use Cwd 'abs_path';
use Digest::MD5 qw(md5_hex);
use Statistics::Descriptive::Smoother;
use File::Temp;

use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
#use lib "/software/polyweb/poly-disk/poly-src/GenBo/lib/obj-nodb/packages";
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/packages";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use Set::IntSpan::Fast::XS;
use image_coverage;
 use List::Util qw( min sum max);
#use Cache_Commons;

my $fork = 1;
my $force;
my ($project_name, $patient_name);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'patient=s'    => \$patient_name,
	'force=s'  => \$force,
	
);


unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
#unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }
my $t = `hostname`;
my $nbErrors = 0;
my $buffer = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProjectCache( -name => $project_name);
 


my $lists = $project-> getListTranscripts();
$patient_name = "all" unless $patient_name;
$project->get_only_list_patients($patient_name);
$project->getCaptures();
$project->getRuns();

#if ($patient_name eq "all" ){
 foreach my $patient (sort {$a->name cmp $b->name} @{$project->getPatients} ){
 		by_gene($patient);
 }
exit(0);
	sub by_gene {
		my ($patient) = @_;
		my $no2 = $patient->getTranscriptsDude("r");
		warn $no2->filename;
		my $no3 = $patient->getGenesDude("c");
		my @levels =("high");
		foreach my $l (@levels){
			get_list_genes($l,$no2,$no3);
		}
		warn $no3->filename;
		warn Dumper $no3->get("ENSG00000120733_5");
		$no3->close();
		$no2->close();
	}
	sub get_list_genes {
		my ($level,$no2,$no3 ) = @_;
		my $array = $no2->get("$level");
		
		warn Dumper $array;#." ".$level;
		my $ts = $project->newTranscripts($array);
		my $h;
		foreach my $t (@$ts){
			my $debug;
			$debug = 1 if $t->getGene->external_name() eq "CEL";
			my $matrix = $no2->get($t->id);
			my $nb = $matrix->{nb};
			my $matrix2 = [map {$_/100} unpack("w".$nb,$matrix->{smooth_expo})];;
			my $dup; 
			my $del;
			foreach my $v (@$matrix2){
			 $dup ++ if $v > 1.4;
			 $del ++ if $v < 0.7;
			}
			$h->{$t->getGene->id}->{dup} += $dup; 
			$h->{$t->getGene->id}->{del} += $del; 
			warn $t->getGene->id." ".$h->{$t->getGene->id}->{dup}." ".$h->{$t->getGene->id}->{dup} if $debug;
		}
		foreach my $g (keys %{$h}){
		#	warn $level." ".$g;
			#warn $g;
			my $type = "dup";
			 if ($h->{$g}->{del} == 0 && $h->{$g}->{dup} == 0){
			 	delete $h->{$g} ;
			 	next;
			 }
			 
			$type = "del" if $h->{$g}->{del} >= $h->{$g}->{dup};
			warn "$level:$type";
			$no3->put($g,"$level:$type"); 
			
		}
		warn Dumper $h;
		$no3->put("genes_".$level,$h);
	}
	



	
	