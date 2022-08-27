#!/usr/bin/perl
use CGI qw/:standard :html3/;
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use colored;
use lib $Bin;
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use Set::IntSpan::Fast::XS;
use List::MoreUtils qw(part);
use List::Util qw(sum);
use Text::Table;
use List::MoreUtils qw(natatime);
use Statistics::Descriptive;
use CGI qw/:standard :html3/;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Parallel::ForkManager;

require "$Bin/../../../GenBo/lib/obj-nodb//packages/cache/polydiag/utility.pm";
require "$Bin/lib/quality_check.pm";


my $projectName;
my $end_ext = "uni";
my $details;
my $fork = 1;
my $vcf_file;
my $log_file;
my $noclean;
my $chr_names;
my $fork;
my $cache;
my $steps;
my $patients_name;

GetOptions(
	'project=s'		=> \$projectName,
	'fork=s'		=> \$fork,
	'cache=s'		=> \$cache,
	'step=s' => \$steps,
	'patients=s' => \$patients_name,
);

my $hresume;
$fork = 1 unless $fork;

my $buffer = GBuffer->new();

my $samtools = $buffer->software("samtools");
#colored::stabilo('blue',"CHECK FILE  : ");
my @types = ();
#@types =  ("files","mendelian","duplicate_regions","coverage_stats","bam_stats","coverage_transcripts","statistics_variations");
#@types =  ("identity","statistics_variations");
my $project = $buffer->newProject( -name => $projectName );
unless($steps){
	 @types = ("files","mendelian","duplicate_regions","coverage_stats","bam_stats","coverage_transcripts","statistics_variations");
	 if ($project->isGenome){
	 	@types = ("files","mendelian","coverage_stats","bam_stats","coverage_transcripts","statistics_variations");
	 	
	 }
}
else {
	 @types = split(",",$steps);
}

#warn $patients_name;

if ($patients_name){
	#warn $patients_name;
	my $z  = $project->get_only_list_patients($patients_name);
}
my( @selected_patient) = map {$_->name } @{$project->getPatients()};
die() unless @selected_patient;
my $no = $project->noSqlQuality("c");
$no->put($project->name,"timestamp",time);
$no->close;
my $steps = {
	"statistics_variations" => sub{
			warn "statistics_variations" ;
			my $buffer = GBuffer->new();
			my $project;


if ($cache == 1){
	#warn "cache";
	#die();
	 $project = $buffer->newProjectCache( -name => $projectName );
	 
}
else {
 $project = $buffer->newProject( -name => $projectName );
}
	
			my $no = $project->noSqlQuality("w");
				$project->get_only_list_patients(join(",",@selected_patient));
			my $a = quality_check::statistics_variations($project,$fork);
			$no->put($project->name,"statistics_variations",$a);
			$no->close();
			
	},
	"identity" => sub{
		warn "identity" ;
			my $buffer = GBuffer->new();
			my $project;


if ($cache == 1){
	 $project = $buffer->newProjectCache( -name => $projectName );
}
else {
 $project = $buffer->newProject( -name => $projectName );
}
			next if $project->isGenome();
			my $no = $project->noSqlQuality("w");
				$project->get_only_list_patients(join(",",@selected_patient));
			$no->put($project->name,"identity",quality_check::identity($project,$fork));
			$no->close();
	},
	"files" => sub{
			warn "files" ;
			my $buffer = GBuffer->new();
				$project->get_only_list_patients(join(",",@selected_patient));
			my $project = $buffer->newProject( -name => $projectName );
			my $no = $project->noSqlQuality("w");
			$no->put($project->name,"files",quality_check::files($project));
			$no->close();
	},
	"mendelian" =>  sub{
			warn "mendelian" ;
			my $buffer = GBuffer->new();
			my $project;


if ($cache == 1){
	warn "cache";
	 $project = $buffer->newProjectCache( -name => $projectName );
}
else {
 $project = $buffer->newProject( -name => $projectName );
}
		$project->get_only_list_patients(join(",",@selected_patient));
			my $no = $project->noSqlQuality("w");
			$no->put($project->name,"mendelian",quality_check::mendelian_statistics($project,$fork));
			$no->close();
	},
	"duplicate_regions" => sub{
				warn "duplicate_regions" ;
				
			my $buffer = GBuffer->new();
			my $project = $buffer->newProject( -name => $projectName );
				$project->get_only_list_patients(join(",",@selected_patient));
			my $no = $project->noSqlQuality("w");
			$no->put($project->name,"duplicate_regions",quality_check::duplicate_regions($project));
			$no->close();
	},
	"coverage_transcripts" =>sub{
		warn "coverage_transcripts" ;
			my $buffer = GBuffer->new();
			my $project;


if ($cache == 1){
	warn "cache";
	 $project = $buffer->newProjectCache( -name => $projectName );
}
else {
 $project = $buffer->newProject( -name => $projectName );
}
	$project->get_only_list_patients(join(",",@selected_patient));
			my $no = $project->noSqlQuality("w");
			$no->put($project->name,"coverage_transcripts",quality_check::coverage_transcripts($project,$fork));
			$no->close();
	},
	
	"bam_stats" => sub{
		warn "bam_stats" ;
			my $buffer = GBuffer->new();
			my $project = $buffer->newProject( -name => $projectName );
				$project->get_only_list_patients(join(",",@selected_patient));
			my $no = $project->noSqlQuality("w");
			$no->put($project->name,"bam_stats",quality_check::bam_stats($project));
			$no->close();
	},
	"coverage_stats" => sub{
		warn "coverage_stats" ;
			my $buffer = GBuffer->new();
				$project->get_only_list_patients(join(",",@selected_patient));
			my $project = $buffer->newProject( -name => $projectName );
			my $no = $project->noSqlQuality("w");
			$no->put($project->name,"coverage_stats",quality_check::coverage_stats($project));
			$no->close();
	},
	
};

foreach my $type (@types){
	$steps->{$type}();
}

exit(0);






 
 