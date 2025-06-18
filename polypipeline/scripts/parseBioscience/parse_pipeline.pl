#!/usr/bin/perl
use Data::Dumper;
use File::Find;
use Getopt::Long;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/";
use lib "$Bin/../../../GenBo/lib/GenBoDB";
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../packages";
use GBuffer;

my $projectName;
my $parent_project_name;
my $patients_name;
my @steps;
my $analysis;
my $run_name_option;
my $kit;
my $chem;
my $target;
my $reuse;
my $no_exec;
my $cpu = 40;
my $help;

GetOptions(
	'project=s'			=> \$projectName,
	'patients=s'		=> \$patients_name,
	'steps|mode=s{1,}'	=> \@steps,
	'run=s'				=> \$run_name_option,
	'kit=s'				=> \$kit,
	'chemistry=s'		=> \$chem,
	'bcr'				=> sub{ $analysis = 'bcr' },
	'tcr'				=> sub{ $analysis = 'tcr' },
	'parent_project=s'	=> \$parent_project_name,
	'target|panel=s'	=> \$target,
	'reuse'				=> \$reuse,
	'no_exec'			=> \$no_exec,
	'cpu=i'				=> \$cpu,
	'help'				=> \$help,
) || die("Error in command line arguments\n");

usage() if ($help);
warn ("--project option is mandatory\n") && usage() unless ($projectName);
#warn ("--analysis should be 'tcr' or 'bcr' or empty for WT analysis\n") && usage() unless ($analysis == undef or $analysis =~ /^bct|tcr$/);
unless ($analysis) {
	warn ("--kit option is mandatory except for TCR/BCR analysis\n") && usage() unless ($kit);
	my @possible_kits = qw(WT_mini WT WT_mega WT_mega_384 WT_penta WT_penta_384);
	warn ("--kit should be one of '".join("','",@possible_kits)."' , given '$chem'\n") && usage() unless( grep($kit, @possible_kits));
	$chem = 'v3' if (grep($kit, @possible_kits[-3..-1]));
}
warn ("--chemistry option is mandatory\n") && usage() unless ($chem);
warn ("--chemistry should be one of ('v1','v2','v3'), given '$chem'\n") && usage() unless( grep($chem, ('v1','v2','v3')) );
warn ("--parent_project option is mandatory for BCR/TCR analysis\n") && usage() if ($analysis =~ /^bcr|tcr$/i and not $parent_project_name);
@steps = qw{all comb} unless (scalar @steps);



my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
my $parent_project = $buffer->newProject( -name => $parent_project_name ) if ($analysis =~ /^bcr|tcr$/i);
my $dir_pipeline = $project->getAlignmentPipelineDir("split-pipe");
#my $dir = $project->getProjectRootPath;
my $dir = $project->getCountingDir("split-pipe");
my $index = $project->getGenomeIndex("split-pipe");
my %sublib;
my $release = $project->annotation_genome_version;

#my $sample_list = get_plate_des($project);

my $runs = $project->getRuns;
my $run;
if (scalar(@$runs) > 1){
	unless ($run_name_option){
		die ("You have ".scalar(@$runs)." runs in the project $projectName. You have to choose one and add --run on the command line\n"
		. map {$_->plateform_run_name." ".$run->date."\n"} @$runs);
	}
	($run) = grep{$_->plateform_run_name eq "$run_name_option"} @$runs;
	die("Unable to find $run_name_option ".$projectName) unless ($run);
}
else {
	$run = $runs->[0];
}

my $all_patients = $project->getPatients;
die("Project $projectName is empty: no patient in project $projectName") unless (scalar @$all_patients);
warn $patients_name;
my $patients = $project->get_only_list_patients($patients_name);
die("No patient $patients_name in project $projectName") unless $patients;



# ALL: process data from each sublibrary individually
if (grep {/all/} @steps) {
	foreach my $p (@$patients){
		my ($fastq1,$fastq2,$dirf) = get_fastq_file($p,$dir_pipeline);
	}
	warn "\n";

	mkdir("$dir/sublibraries") unless (-d "$dir/sublibraries");
	open(my $jobs_all, ">");
#	warn ("$dir/jobs_all.txt");
	foreach my $subl (sort keys(%sublib)){
#	warn $subl;
		my $subcmd = "singularity run";
		# si erreur ne trouve pas la librairie libxml2: "error while loading shared libraries: libxml2.so.2"
#		$subcmd .= ' --env LD_LIBRARY_PATH=/miniconda/lib:$LD_LIBRARY_PATH';
		$subcmd .= " -B $dir:/PROJECT/ -B $index:/REF/" unless ($analysis);
		if ($parent_project and $analysis =~ /^bcr|tcr$/) {
			my $parent_dir;
			$parent_dir = $parent_project->getCountingDir("split-pipe") if (-d $parent_project->getCountingDir("split-pipe").$1);
			$parent_dir = $parent_project->getProjectRootPath if (-d $parent_project->getProjectRootPath.$1);
			die ("No parent directory found for the corresponding WT sublibrary '$1' in project $parent_project_name") unless (-d $parent_dir);
			$subcmd .= " -B $parent_dir:/PARENT/";
		}
		$subcmd .= " -B $dir_pipeline:/FASTQ/ /data-beegfs/software/sif/splitpipe.1.5.1.sif split-pipe --mode all";
		$subcmd .= " --chemistry $chem ";
		if ($analysis =~ /^bcr|tcr$/) {
			$subcmd .= " --$analysis\_analysis";
			$subcmd .= ' --immune_genome human' if ($release =~ /^HG/);
			$subcmd .= ' --immune_genome mouse' if ($release =~ /^MM/);
			die ("BCR/TCR analysis only supported for human and mouse. Actual release: $release") unless ($release =~ /^HG|MM/);
			$subl =~ /(.*)_[bt]cr/i;
			$subcmd .= " --parent_dir /PARENT/$1" if ($parent_project);
			
		}
		else { # WT analysis
			$subcmd .= " --kit $kit";
			$subcmd .= " --genome_dir /REF " ;
			my $plate_des = $run->sample_sheet;
			die ("No SampleLoadingTable. Please upload the SampleLoadingTable excel file to the run document.") unless ( $plate_des);
			my $csv_tmp = $dir."plate.xlsm";
#			warn $csv_tmp;
			open(TOTO,">$csv_tmp");
			print TOTO $plate_des;
			close TOTO;
			$subcmd .= " --samp_sltab /PROJECT/plate.xlsm";
			$subcmd .= " --targeted_list /PROJECT/$target" if ($target);
		}
		$subcmd .= " --output_dir /FASTQ/sublibraries/$subl ";
		my $fastq1 = $sublib{$subl}{R1};
		my $fastq2 = $sublib{$subl}{R2};
		$subcmd .= " --fq1 /FASTQ/".$fastq1." --fq2 /FASTQ/".$fastq2;
		$subcmd .= " --reuse" if ($reuse);
		
		$subcmd .= " && cp -r $dir_pipeline/sublibraries/$subl $dir/sublibraries/";
		#	warn $csv_tmp;
		print $subcmd."\n";
		print {$jobs_all} $subcmd."\n";
	}
	close ($jobs_all);
	my $exit = system("cat $dir/jobs_all.txt | run_cluster.pl -cpu=$cpu") unless ($no_exec);
	die if ($exit);
}
warn "\n";



# COMB: combine the processed data from each sublibrary into a single dataset
if (grep {/comb(ine)?/} @steps) {
	my @names = map{ "/PROJECT/sublibraries/".$_->name()} @$all_patients;
	my $sublib = join (" ", sort @names);
#	warn $sublib;
	my $cmd2 = "singularity run -B  $dir:/PROJECT/ -B $index:/REF/ -B $dir_pipeline:/FASTQ/";
	if ($parent_project and $analysis =~ /^bcr|tcr$/) {
		my $parent_dir;
		$parent_dir = $parent_project->getCountingDir("split-pipe").'/comb/' if (-d $parent_project->getCountingDir("split-pipe").'/comb/');
		$parent_dir = $parent_project->getCountingDir("split-pipe") if (-d $parent_project->getCountingDir("split-pipe"));
		$parent_dir = $parent_project->getProjectRootPath if (-d $parent_project->getProjectRootPath);
		die ("No parent directory found for the corresponding WT sublibrary '$1' in project $parent_project_name") unless (-d $parent_dir);
		$cmd2 .= " -B $parent_dir:/PARENT/";
	}
	$cmd2 .= " /data-beegfs/software/sif/splitpipe.1.5.1.sif";
	$cmd2 .= " split-pipe --mode comb" ;
	if ($analysis =~ /^bcr|tcr$/) {
		$cmd2 .= ' --immune_genome human' if ($release =~ /^HG/);
		$cmd2 .= ' --immune_genome mouse' if ($release =~ /^MM/);
		$cmd2 .= " --parent_dir /PARENT" if ($parent_project);
	}
	$cmd2 .= " --output_dir /FASTQ/comb" ;
	$cmd2 .= " --sublibraries $sublib" ;
	$cmd2 .= " --reuse" if ($reuse);
	$cmd2 .= " && cp -r $dir_pipeline/comb $dir";
	print $cmd2."\n";
	
	open(my $jobs_comb, ">$dir/jobs_comb.txt");
	print {$jobs_comb} $cmd2;
	close($jobs_comb);
	
	my $exit = system("cat $dir/jobs_comb.txt | run_cluster.pl -cpu=$cpu") unless ($no_exec);
	die if ($exit);
}
warn "\n";


sub get_fastq_file {
	my ($patient,$dir_pipeline) = @_;
	my $name=$patient->name();
	my $dir_fastq = $patient->getSequencesDirectory() ;
#	warn $dir_fastq;
#	$dir_fastq = $patient->getSequencesDirectory() unless $dir_fastq;
	my $fastq1 = $name."_R1_L001.fastq.gz";
	my $fastq2 = $name."_R2_L001.fastq.gz";
	my @r1 = glob($dir_fastq."/*".$name."*_R1*.fastq.gz");
	my @r2 = glob($dir_fastq."/*".$name."*_R2*.fastq.gz");
#	warn substr($r1[0],rindex($r1[0],"/")+1,length($r1[0]));
#	my @r1_bis = map{substr($_,rindex($_,"/"),length($_))};
	$fastq1 = substr($r1[0],rindex($r1[0],"/")+1,length($r1[0])) if scalar(@r1) == 1;
	$fastq2 = substr($r2[0],rindex($r2[0],"/")+1,length($r2[0])) if scalar(@r2) == 1;
	my $cmd1 = join(" ",@r1);
	my $cmd2 = join(" ",@r2);
	
	$sublib{$name}{R1} = $fastq1;
	$sublib{$name}{R2} = $fastq2;

	return ($fastq1,$fastq2,$dir_fastq) if (-e $fastq1 && -e $fastq2);
	warn "cat $cmd1 > ".$dir_pipeline."/".$fastq1;
	
	unless ($no_exec) {
		system "cat $cmd1 > ".$dir_pipeline."/".$fastq1 unless -e $dir_pipeline."/".$fastq1;
		system "cat $cmd2 > ".$dir_pipeline."/".$fastq2 unless -e $dir_pipeline."/".$fastq2 ;
	}
	
	return  ($fastq1,$fastq2,$dir_pipeline);
}




sub get_plate_des{
	my ($project) = @_;
#	my $patients = $project->get_only_list_patients($patients_name);
	my @cond;
	my @names;
	foreach my $p (@$patients){
		
		my $group = uc($p->somatic_group());
		next() if $group eq "SUB";
		my $bc = $p->barcode();
		$bc =~ s/-/:/;
		my $xcond = "--sample ".$p->name." ".$bc;
		push(@cond, $xcond);
		push(@names, "/FASTQ/".$p->name());
	}
	my $list = join(" ",@cond);
	return $list;
}




sub usage {
	print "
$0
-------------------
Obligatoires:
	project <s>			Nom du projet
	chemistry <s>			Version de la chimie
					Valeurs possibles: 'v1', 'v2', 'v3'
	kit <s>				Version du kit. Obligatoire sauf pour l'analyse de BCR ou TCR.
					Valeurs possibles: 'WT_mini', 'WT', 'WT_mega', 'WT_mega_384', 'WT_penta', 'WT_penta_384'
	
Optionels:
	patients <s>			Noms de patients/échantillons, séparés par des virgules (utilisé seulement pour l'étape 'all'
	steps/mode <s>			Etapes à réaliser:
					all -> analyse des sublibrairies individuellement,
					comb -> combine les données traitées de chaque sublibrairie
					par défaut: fait les deux à la suite
	run <s>				Numéro du run si le projet en contient plusieurs (utilisé pour récupérer la plan de plaque)
	bcr/tcr				Analyse des BCR ou TCR
	parent_project <s>		Chemin des résultats de l'analyse whole transcriptome associée. Pour l'analyse de BCR/TCR uniquement
	target/panel <s>		Fichier csv spécifiant la liste des gènes cibles.
					Le fichier csv dois suivre le format suivant:
					gene_id,gene_name
					ENSG00000003096,KLHL13
	reuse				Re-use existing files if found (vs generate fresh)
	cpu <i>				Nombre de cpu à utiliser, défaut: 40
	no_exec				Ne pas exécuter les commandes
	help				Affiche ce message

";
	exit(1);

}


