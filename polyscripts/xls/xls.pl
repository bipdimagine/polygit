#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
#use lib "/data-isilon/bipd-src/pnitschk/git2-polyweb/polygit/GenBo/lib/obj-nodb";
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
#use lib "/software/polyweb/poly-disk/poly-src/GenBo/lib/obj-nodb/";
use GBuffer;
use Set::IntSpan::Fast::XS;
use Getopt::Long;
my $fork = 1;

my $file;
my $transcript = "all";
my $cache ;
GetOptions(
	'file=s'    => \$file,
	'transcript' => \$transcript,
	'cache=s' => \$cache,
);
##
my $root_cmd = "perl /software/polyweb/poly-disk/www/cgi-bin/polymorphism-cgi//cache_nodb/scripts/cache_polydiag.pl";
	
	
my $opt_transcript = "transcripts=$transcript ";
open (FILE,"$file");
my $cmd = qq{/software/polyweb/poly-disk/www/cgi-bin/polymorphism-cgi/validation_variation/patient_report.pl edit_mode=1 never=1 this=6 impact=2 frequence=2 allele_quality=5 report_mode=1 project_summary=1 graphic_coverage=1 sanger_variations=1 validated_variations=1 todo_variations=1 list=1 span=20 limit=30  user_name=pnitschk xls=1 };
while(<FILE>){
	chomp();
	my $buffer = new GBuffer;
	my $project = $buffer->newProjectCache( -name => $_);
	my $dir_out= " /data-beegfs/tmp/".$project->name;
	system ("mkdir $dir_out") unless -e $dir_out;
		my $ppn=1;
		foreach my $patient (@{$project->getPatients}){
			my $cmd_root = " $root_cmd -fork=$ppn -project=".$project->name." -patient=".$patient->name." && ";
			my $opt2 = " $opt_transcript patients=".$patient->name." project=".$project->name;
			my $cmd3.=$opt2." | tail -n +4 > $dir_out/".$patient->name.".xls";
			
			if ($cache){
				#print $cmd_root."\n";
				print $cmd_root.$cmd." $cmd3\n";	
				
			}
			else {
				print $cmd." $cmd3\n";
			}
			#
		}
	#delete $project->{objects};
	#$buffer = undef;
}

