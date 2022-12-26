#!/usr/bin/perl
use FindBin qw($Bin $RealBin);
use strict;

use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages/"; 
use GBuffer ;
use Data::Dumper;
use Getopt::Long;
use Carp;
use Bio::DB::Sam;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use List::Util qw(sum);
use File::Temp;
 use Time::Elapsed qw( elapsed );
 use Time::ETA;
use Storable qw(store retrieve freeze);
use JSON::XS;

my $buffer = GBuffer->new();

my $project_name;
my $p1name;
my $p2name;
GetOptions(
	'project=s' => \$project_name,
	'giab=s' => \$p2name,
	#'low_calling=s' => \$low_calling,
);
#NGS2021_3851


my $project = $buffer->newProject( -name => "$project_name" );
my $patient = $project->getPatient($p2name);
my $is_ucsc = undef;
$is_ucsc = 1 if $project->getChromosome(1)->fasta_name eq "chr1";
update_type($buffer->dbh,$project->id);
die();
#my $file_out = $patient->vcf_file;
my $GIAB_DIR;

if ($p2name =~ /001/){
	$GIAB_DIR = "/data-isilon/public-data/repository/HG19/GIAB/HG001";
}
elsif ($p2name =~ /002/){
	$GIAB_DIR = "/data-isilon/public-data/repository/HG19/GIAB/HG002";
}
else {
	die($p2name." not foud HG001 or HG002");
}
my $dir_pipeline = $project->getVariationsDir("bed");
#intspanToBed
#getIntSpanCapture
my $bed = $dir_pipeline."/".$patient->name.time.".bed";
my $list = $dir_pipeline."/".$patient->name.".list";
my $vcf1 = $dir_pipeline."/".$patient->name.".vcf";
unlink $bed if -e $bed;
open(BED,">$bed");
open(LIST,">$list");
foreach my $chr (@{$project->getChromosomes}){
	my $chr_name = $chr->name;
	#next unless $chr_name eq "18";
	print LIST $chr->name."\t".$chr->fasta_name."\n";
	#getIntSpanCapture
	#getIntSpanCaptureForCalling
	my @line = intspanToBed($chr_name,$chr->getIntSpanCapture(1000) );
	print BED join("\n",@line)."\n" if @line;
}

close(BED);
close (LIST);
warn $list;


my $vcf_giab = $GIAB_DIR."/vcf.gz";
my $bcftools = $buffer->software("bcftools");
my $bgzip = $buffer->software("bgzip");
my $tabix = $buffer->software("tabix");
my $fileout = $project->getVariationsDir("haplotypecaller4")."/".$patient->name.".vcf.gz";
system("$bgzip $bed ; $tabix $bed.gz");

$bed =$bed.".gz";
die() unless -e $bed.".tbi";

unlink  $fileout.".tbi" if -e $fileout.".tbi";
unlink  $fileout if -e $fileout;
warn "$bcftools view $vcf_giab -R $bed | $bcftools annotate --rename-chrs $list - -o $fileout -O z ";
system("$bcftools view $vcf_giab -R $bed | $bcftools annotate --rename-chrs $list - -o $fileout -O z ");#&& $tabix -f -p vcf $fileout");
#die();
#system("$bcftools view $vcf_giab | $bcftools annotate --rename-chrs $list - -o $fileout -O z ");
system("$tabix -f -p vcf $fileout  ");
die() unless -e $fileout.".tbi";
system("/software/polyweb/poly-disk/poly-src/polygit/polypipeline/scripts/scripts_cache/step/./step1.pl -project=NGS2022_6106 | run_cluster.pl -cpu=5");
warn "\n\n--------------------\n\n";
warn "run stats ";
system("$Bin/control_panel_giab.pl -project=$project_name -giab=$p2name");

warn "\n-----------------------------\n\n";
warn "-----------------------------\n\n";
warn "------ THAT's ALL FOLKS ------\n\n";
warn "-----------------------------\n\n";
warn "-----------------------------\n\n";
exit(0);

warn $fileout;
warn "$bcftools view $vcf_giab -R $bed | $bcftools annotate --rename-chrs $list - -o $fileout -O z";
#if ($is_ucsc){
#	open(BCF, "bcftools view $vcf_giab -R $bed|");
#	
#	while(BCF){
#		
#	}
#	
#}

sub update_type{
	my ($dbh,$project_id) = @_;
	my $query = qq{
		update  PolyprojectNGS.projects  set validation_db = 'control_giab' where project_id=$project_id;
	};
	$dbh->do($query) ;
}
sub intspanToBed{
	my ($chr_name,$intspan) = @_;
	my $iter = $intspan->iterate_runs();
	my @tt;
	my $size = 0;
    while (my ( $from, $to ) = $iter->()) {
    	#warn "+".$from."-".$to if ($from<=31263320);
    	warn $from."-".$to if ($from<=31263320 && $to >= 31263320);
    	$size += abs($from-$to);
    		push(@tt,$chr_name."\t".$from."\t".$to);
    	
    }
    warn $chr_name." ".$size;
	#	my @tt = map{$_ = $chr->ucsc_name."\t".$_} split(";",$intspan->as_string({ sep => ";", range => "\t" }) );
		return @tt;
}