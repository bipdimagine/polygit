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
use lib "$Bin/../../packages";
use lib "$Bin/../../../polypipeline/dragen/scripts/";
use dragen_util;
use file_util;
use GBuffer;

my $projectName;
my $parent_project_name;
my $patients_name;
my @steps;
my $analysis;
my $sample_list;
my $run_name_option;
my $kit;
my $chem;
my $target;
my $reuse;
my $no_exec;
my $dry_run;
my $other_opt;
my $cpu = 40;
my $help;

GetOptions(
	'project=s'					=> \$projectName,
	'patients=s'				=> \$patients_name,
	'steps|mode=s{1,}'			=> \@steps,
	'run=s'						=> \$run_name_option,
	'kit=s'						=> \$kit,
	'chemistry=s'				=> \$chem,
	'sample_list=s'				=> \$sample_list,
	'bcr'						=> sub{ $analysis = 'bcr' },
	'tcr'						=> sub{ $analysis = 'tcr' },
	'parent_project=s'			=> \$parent_project_name,
	'target|panel=s'			=> \$target,
	'reuse'						=> \$reuse,
	'no_exec'					=> \$no_exec,
	'dry_run'					=> \$dry_run,
	'other_options|options=s'	=> \$other_opt,
	'cpu=i'						=> \$cpu,
	'help'						=> \$help,
) || die("Error in command line arguments\n");

usage() if ($help);
warn ("--project option is mandatory\n") && usage() unless ($projectName);
#warn ("--analysis should be 'tcr' or 'bcr' or empty for WT analysis\n") && usage() unless ($analysis == undef or $analysis =~ /^bct|tcr$/);
unless ($analysis) {
	warn ("--kit option is mandatory except for TCR/BCR analysis\n") && usage() unless ($kit);
	my @possible_kits = qw(WT_mini WT WT_mega WT_mega_384 WT_penta WT_penta_384);
	warn ("--kit should be one of '".join("','",@possible_kits)."' , given '$chem'\n") && usage() unless ( grep($kit, @possible_kits));
	$chem = 'v3' if (grep($kit, @possible_kits[-3..-1]));
}
warn ("--chemistry option is mandatory\n") && usage() unless ($chem);
warn ("--chemistry should be one of ('v1','v2','v3'), given '$chem'\n") && usage() unless ( grep($chem, ('v1','v2','v3')) );
warn ("Sample list '$sample_list' doesn't exist\n") && usage() unless ( $sample_list and not -e $sample_list );
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
my $patients = $project->get_only_list_patients($patients_name);
die("No patient $patients_name in project $projectName") unless $patients;



# ALL: process data from each sublibrary individually
if (grep {/all/} @steps) {
	my ($fastq1,$fastq2,$dirf);
	foreach my $p (sort {$a->name cmp $b->name} @$patients){
	warn "cp fastq ".$p->name;
		($fastq1,$fastq2,$dirf) = dragen_util::get_fastq_file($p,$dir_pipeline);
	}
	warn "\n";

	mkdir("$dir/sublibraries") unless (-d "$dir/sublibraries");
	open(my $jobs_all, ">", $dir."jobs_all.txt") or die("Can't open '$dir/jobs_all.txt': $!");
#	warn ("$dir/jobs_all.txt");
#	warn sort keys(%sublib);
	foreach my $subl (sort {$a->name cmp $b->name} @$patients){
#	warn $subl;
		my $name = $subl->name;
		my $subcmd = "singularity run";
		# si erreur ne trouve pas la librairie libxml2: "error while loading shared libraries: libxml2.so.2"
#		$subcmd .= ' --env LD_LIBRARY_PATH=/miniconda/lib:$LD_LIBRARY_PATH';
		$subcmd .= " -B $dir -B $index" unless ($analysis);
		my $parent_dir;
		if ($parent_project and $analysis =~ /^bcr|tcr$/) {
			$name =~ /(.*)_[bt]cr/i;
			$parent_dir = $parent_project->getProjectRootPath if (-d $parent_project->getProjectRootPath.$1);
			$parent_dir = $parent_project->getCountingDir("split-pipe") if (-d $parent_project->getCountingDir("split-pipe").$1);
			die ("No parent directory found for the corresponding WT sublibrary '$1' in project $parent_project_name") unless (-d $parent_dir);
			$subcmd .= " -B $parent_dir";
		}
		$subcmd .= " -B $dir_pipeline /data-beegfs/software/sif/splitpipe.1.5.1.sif split-pipe --mode all";
		$subcmd .= " --chemistry $chem ";
		if ($analysis =~ /^bcr|tcr$/) {
			$subcmd .= " --$analysis\_analysis";
			$subcmd .= ' --immune_genome human' if ($release =~ /^HG/);
			$subcmd .= ' --immune_genome mouse' if ($release =~ /^MM/);
			die ("BCR/TCR analysis only supported for human and mouse. Actual release: $release") unless ($release =~ /^HG|MM/);
			$name =~ /(.*)_[bt]cr/i;
			$subcmd .= " --parent_dir $parent_dir$1" if ($parent_project);
			
		}
		else { # WT analysis
			$subcmd .= " --kit $kit";
			$subcmd .= " --genome_dir $index " ;
			unless ($sample_list) {
				my $plate_des = $run->sample_sheet;
				die ("No SampleLoadingTable. Please upload the SampleLoadingTable excel file to the run document.") unless ($plate_des);
				my $csv_tmp = $dir."SampleLoadingTable.xlsm";
#				warn $csv_tmp;
				open(TOTO,">$csv_tmp");
				print TOTO $plate_des;
				close TOTO;
				$subcmd .= " --samp_list $dir/sample-list.txt";
				$subcmd .= " --samp_sltab $dir/SampleLoadingTable.xlsm";
			}
			if ($sample_list) {
				$subcmd .= " --samp_list $dir/sample-list.txt";
			}
			$subcmd .= " --targeted_list $dir$target" if ($target);
		}
		$subcmd .= " --output_dir $dir_pipeline/sublibraries/$name ";
#		my $fastq1 = $sublib{$subl}{R1};
#		my $fastq2 = $sublib{$subl}{R2};
		$subcmd .= " --fq1 $fastq1 --fq2 $fastq2";
		$subcmd .= " --reuse" if ($reuse);
		$subcmd .= " --dryrun" if ($dry_run);
		$subcmd .= ' '.$other_opt if ($other_opt);
		
		$subcmd .= " && cp -r $dir_pipeline/sublibraries/$name $dir/sublibraries/" unless ($dry_run);
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
	my @names = map{ $dir."sublibraries/".$_->name()} @$all_patients;
	my $cmd2 = "singularity run -B  $dir -B $index -B $dir_pipeline";
	my $parent_dir;
	if ($parent_project and $analysis =~ /^bcr|tcr$/) {
		$parent_dir = $parent_project->getProjectRootPath if (-d $parent_project->getProjectRootPath);
		$parent_dir = $parent_project->getCountingDir("split-pipe") if (-d $parent_project->getCountingDir("split-pipe"));
		$parent_dir = $parent_project->getCountingDir("split-pipe").'/comb/' if (-d $parent_project->getCountingDir("split-pipe").'/comb/');
		die ("No parent directory found for the corresponding combined WT sublibrary in project $parent_project_name") unless (-d $parent_dir);
		$cmd2 .= " -B $parent_dir";
	}
	$cmd2 .= " /data-beegfs/software/sif/splitpipe.1.5.1.sif";
	$cmd2 .= " split-pipe --mode comb" ;
	if ($analysis =~ /^bcr|tcr$/) {
		$cmd2 .= ' --immune_genome human' if ($release =~ /^HG/);
		$cmd2 .= ' --immune_genome mouse' if ($release =~ /^MM/);
		$cmd2 .= " --parent_dir $parent_dir" if ($parent_project);
	}
	$cmd2 .= " --output_dir $dir_pipeline/comb" ;
	$cmd2 .= " --sublibraries ".join (" ", sort @names) ;
	$cmd2 .= " --reuse" if ($reuse);
	$cmd2 .= " --dryrun" if ($dry_run);
	$cmd2 .= ' '.$other_opt if ($other_opt);
	$cmd2 .= " && cp -r $dir_pipeline/comb $dir" unless ($dry_run);
	print $cmd2."\n";
	
	open(my $jobs_comb, ">$dir/jobs_comb.txt");
	print {$jobs_comb} $cmd2;
	close($jobs_comb);
	
	my $exit = system("cat $dir/jobs_comb.txt | run_cluster.pl -cpu=$cpu") unless ($no_exec);
	die if ($exit);
}
warn "\n";


#sub get_fastq_file {
#	my ($patient,$dir_pipeline) = @_;
#	my $name=$patient->name();
#	my $dir_fastq = $patient->getSequencesDirectory() ;
##	warn $dir_fastq;
##	$dir_fastq = $patient->getSequencesDirectory() unless $dir_fastq;
#	my $fastq1 = $name."_R1_L001.fastq.gz";
#	my $fastq2 = $name."_R2_L001.fastq.gz";
#	my @r1 = glob($dir_fastq."/*".$name."*_R1*.fastq.gz");
#	my @r2 = glob($dir_fastq."/*".$name."*_R2*.fastq.gz");
##	warn substr($r1[0],rindex($r1[0],"/")+1,length($r1[0]));
##	my @r1_bis = map{substr($_,rindex($_,"/"),length($_))};
#	$fastq1 = substr($r1[0],rindex($r1[0],"/")+1,length($r1[0])) if scalar(@r1) == 1;
#	$fastq2 = substr($r2[0],rindex($r2[0],"/")+1,length($r2[0])) if scalar(@r2) == 1;
#	my $cmd1 = join(" ",@r1);
#	my $cmd2 = join(" ",@r2);
#	
#	$sublib{$name}{R1} = $fastq1;
#	$sublib{$name}{R2} = $fastq2;
#
#	return ($fastq1,$fastq2,$dir_fastq) if (-e $fastq1 && -e $fastq2);
#	warn "cat $cmd1 > ".$dir_pipeline."/".$fastq1;
#	
#	unless ($no_exec) {
#		system "cat $cmd1 > ".$dir_pipeline."/".$fastq1 unless -e $dir_pipeline."/".$fastq1;
#		system "cat $cmd2 > ".$dir_pipeline."/".$fastq2 unless -e $dir_pipeline."/".$fastq2 ;
#	}
#	
#	return  ($fastq1,$fastq2,$dir_pipeline);
#}




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
		push(@names, $dir_pipeline.$p->name());
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
	sample_list <s>			Liste des échantillons et de leurs puis correspondants, séparés par un espace.
					Si omis, utilise le document dans la base de donnée comme plan de plaque
	bcr/tcr				Analyse des BCR ou TCR
	parent_project <s>		Chemin des résultats de l'analyse whole transcriptome associée. Pour l'analyse de BCR/TCR uniquement
	target/panel <s>		Fichier csv spécifiant la liste des gènes cibles.
					Le fichier csv dois suivre le format suivant:
					gene_id,gene_name
					ENSG00000003096,KLHL13
	cpu <i>				Nombre de cpu à utiliser, défaut: 40
	reuse				Re-utilise les fichiers existant si trouvés (vs génère de nouveau)
	no_exec				Ne pas exécuter les commandes
	dry_run				Initialise et reporte le status puis s'arrête
	other_options			Autres options pour split-pipe, à écrire entre guillemets
	help				Affiche ce message

";
	exit(1);

}


