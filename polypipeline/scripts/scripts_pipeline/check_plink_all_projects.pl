#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use Data::Dumper;
use Carp;
use GBuffer;

my $buffer = GBuffer->new();
my @lProjects = ('NGS2013_0300');
foreach my $projectName (@lProjects) {
	print "\n\n##### $projectName\n";
	my $project;
	my $vcfCatFile;
	my $hVcfFiles;
	eval { $project = $buffer->newProject( -name => $projectName ); };
	next if ($@);
	my $build = $project->buffer->build();
	
	my @lPatientsObj;
	eval { @lPatientsObj = @{$project->getPatients()}; };
	next if ($@);
	
	my $pedFile;
	eval {
		$pedFile = $project->getPedigreeFile();
		next unless (-e $pedFile);
	};
	
	foreach my $pat (@lPatientsObj) {
		my (@lSnpsFiles, @lIndelsFiles);
		eval { @lSnpsFiles = @{$pat->getVariationsFiles()}; };
		if (scalar(@lSnpsFiles) > 0) { $hVcfFiles->{$pat->name()}->{'snps'} = join(',', @lSnpsFiles); }
	}
	eval { unless (-d $project->getCallingPipelineDir('indiv_merge')) { makedir($project->getCallingPipelineDir('indiv_merge')); } };
	my $cmd1 = "/opt/java/latest/bin/java -jar /bip-d/soft/distrib/GATK/gatk-latest.2.7/GenomeAnalysisTK.jar -T CombineVariants -R /data-xfs/public-data/$build/genome/fasta/all.fa ";
	foreach my $k (keys(%$hVcfFiles)) {
		foreach my $file (split(',', $hVcfFiles->{$k}->{'snps'})) { $cmd1 .= " --variant:$k $file"; }
		#foreach my $file (split(',', $hVcfFiles->{$k}->{'indels'})) { $cmd1 .= " --variant:$k $file"; }
	}
	$vcfCatFile = $project->getCallingPipelineDir('indiv_merge').'/'.$projectName.'_indiv_combine.vcf';
	$cmd1 .= "  -o $vcfCatFile";
	$cmd1 .= "  -genotypeMergeOptions PRIORITIZE";
	$cmd1 .= "  -priority ".join(',', sort(keys(%$hVcfFiles)));
	
	warn "CMD1: ".$cmd1;
	if ($cmd1 =~ /variant/) { `$cmd1`; }

	if (-e $vcfCatFile) {
		if (-e "$vcfCatFile.gz") {
			my $cmdRm = "rm $vcfCatFile.gz";
			`$cmdRm`;
		}
		my $cmd2 = "bgzip $vcfCatFile";
		warn "CMD2: ".$cmd2;
		`$cmd2`;
		$vcfCatFile .= '.gz';
	}
	else { next; }
	if ($@) { next; }
	unless (-e $vcfCatFile) { next; }
	#my $cmd3 = "perl /bip-d/perl/GenBo/script/ngs_exome/pipeline/plink.pl -project=$projectName -vcf=$vcfCatFile -ped=$pedFile";
	my $cmd3 = "perl /data-xfs/dev/mbras/Workspace/GenBo/script/ngs_exome/pipeline/plink.pl -project=$projectName -vcf=$vcfCatFile -ped=$pedFile -indiv=true";
	warn "CMD3: ".$cmd3;
	my $res = `$cmd3`;
	warn Dumper $res;
	warn "\nPROJECT NAME: $projectName\n";
	$project = undef;
	print "\nLaunch NEXT ?? press any key....";
	my $next = <>;
}