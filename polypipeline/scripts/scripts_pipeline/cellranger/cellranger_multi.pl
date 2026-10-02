#!/usr/bin/env perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/";
use lib "$Bin/../../../../GenBo/lib/GenBoDB";
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages";
use Logfile::Rotate;
use Cwd;
use PBS::Client;
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moose;
use MooseX::Method::Signatures;
#use bds_steps;   
use file_util;
use Class::Inspector;
use Digest::MD5::File ;
use GBuffer;
use GenBoProject;
use colored; 
use Config::Std;
use Text::Table;
use File::Temp qw/ tempfile tempdir /;
use Term::Menus;
use Proc::Simple;
use Storable;
use JSON::XS;
use XML::Simple qw(:strict);
use Cwd 'abs_path';
use File::Path qw(make_path);
use Carp;
use autodie qw(system open);

  

my $projectName;
my $patients_name;
my @steps;
my $lane;
my $mismatch = 0;
my $feature_ref;
my $cmo_ref;
my $no_exec;
my $aggr_name;
my $choose_exec;
my $choose_transcriptome;
my $chemistry;
my $create_bam;
my $multi_meth;
my $cpu = 20;
$cpu = 64 if (`hostname` =~ /^master/);
my $force;
my $project_version;
my $help;

GetOptions(
	'project=s'										=> \$projectName,
	'patients=s'									=> \$patients_name,
	'steps=s{1,}'									=> \@steps,
	'multi_meth=s'									=> \$multi_meth,
	'mismatches=i'									=> \$mismatch,
	'create_bam!'									=> \$create_bam,
	'feature_ref|feature_csv=s'						=> \$feature_ref,
#	'cmo_ref|cmo_csv=s'								=> \$cmo_ref,
	'aggr_name=s'									=> \$aggr_name,
	'choose_exec|choose_version|version'			=> \$choose_exec,
	'choose_transcriptome|choose_reference'			=> \$choose_transcriptome,
	'chemistry=s'									=> \$chemistry,
	'cpu=i'											=> \$cpu,
	'force'											=> \$force,
	'no_exec'										=> \$no_exec,
	'version=s'										=> \$project_version,
	'help'											=> \$help,
) || confess("Error in command line arguments\n");

usage() if $help;
die("-project argument is mandatory") unless ($projectName);
my $buffer = GBuffer->new();
my $max_cpu = 40;
$max_cpu = 256 if ($buffer->biocluster);
die("cpu must be in [1;$max_cpu], given $cpu") unless ($cpu > 0 and $cpu <= $max_cpu);

my $project = $buffer->newProject( -name => $projectName , -version => $project_version);
my $all_patients = $project->getPatients;
my $patients = $project->get_only_list_patients($patients_name);
die("No patient in project $projectName") unless ($patients);

# Vérifie les caractères non acceptés
my @patient_names = map {$_->name} @$patients;
my @invalid_names = grep { $_ !~ /^[A-Za-z0-9_-]+$/ } @patient_names;
die ("Patient names can only contain letters, numbers, hyphens and underscores. Sample names ".join(@invalid_names, ', ')." are invalids.") if (@invalid_names);

# Vérifie que le projet est en somatic et que les groupes sont correctement remplis
unless ($project->isSomatic) {
	my $warn = "Project $projectName is not in somatic mode. "
		."Activate somatic mode and check that the groups have been filled in, so that they can be taken into account in the analysis.";
#	die ($warn) if (grep{/@steps/} ('count', 'aggr'));
	die ($warn);
}
my @groups = map {$_->somatic_group} @$patients;
die("Check that the groups are not empty and correctly filled. Accepted groups are exp, nuclei, adt, vdj") 
	unless (grep{/^(exp|nuclei|adt|vdj)$/i} @groups);

# Warn PolyProject
warn (qq{l
Les échantillons doivent être rentrés dans PolyProject comme suit:
Patient	Family	Group	BC	BC2
sample_name	sample_pool	exp/adt/vdj	SI-XX-XXX	BC001

Pour ce projet $projectName:
Patient	Family	Group	BC	BC2
} . join("\n", map{$_->name."\t".$_->family."\t".$_->somatic_group."\t".$_->barcode."\t".$_->barcode2} sort {$a->family cmp $b->family || $a->name cmp $b->name} @$patients)."\n");
#my %families = map {%families{$_->family} ++} @$patients;
#if (grep {/vdj|adt/} keys %families or min values %families <= 1) {
if (grep {/vdj|adt/} map {$_->family()} @$patients) {
	die("\nbye") unless (prompt("Are you sure the family/pool names are correctly filled ? ", -yes_no));
}
sleep(2);
die("\nbye") if (prompt("Tap any key to continue, 'q' to quit.", -w=>'q', -e=>'', -1));
print("\n");


# Steps
@steps = split(/,/, join(',',@steps));
my $list_steps = ['dragen_demultiplex', 'teleport', 'count', 'aggr', 'tar', 'cp', 'cp_web_summary', 'all'];
#my @correct_steps = grep{/@steps/} @$list_steps;
#undef @steps if (@steps and scalar @correct_steps == 0);
unless (@steps) {
	my %Menu_1 = (
		Item_1 => {
			Text   => "]Convey[",
			Convey => $list_steps,
		},
		Select => 'Many',
		Banner => "   Select steps:"
	);
	@steps = &Menu( \%Menu_1 );
	die if ( @steps eq ']quit[' );
}
warn 'steps='.join(',',@steps);

# Multi method
if (grep(/^count$/i, @steps)) {
	my @multi_meth_availables = qw{5P_singleplex flex_multiplex flex_singleplex ocm hash};
	unless (grep{/^$multi_meth$/} @multi_meth_availables) {
		$multi_meth = prompt("Choose a sample multiplexing method for project $projectName \"".$project->description.'":',
				-menu=>\@multi_meth_availables);
	}
	warn 'multi method: '.$multi_meth;
}


my $run = $project->getRun();
my $run_name = $run->plateform_run_name;
my $type = $run->infosRun->{method};
my $machine = $run->infosRun->{machine};

# Executable & version
my $exec = "cellranger";
$exec .= '-atac' if ($type eq 'atac');
$exec .= '-arc' if ($type eq 'arc');
$exec = 'spaceranger' if  ($type eq 'spatial');
if ($choose_exec) {
	my $exec_type = $exec;
	my $path_exec = "/data-bipd/data-pure/software/distrib/$exec_type/";
	opendir(my $dh, $path_exec) || die "Can't opendir '$path_exec': $!";
	my @exec = reverse grep {-x "$path_exec$_/$exec_type" && ! /^\./} readdir($dh);
	closedir ($dh);
	confess("No directory with executable '$exec_type' found in '$path_exec'") unless (scalar @exec);
	$exec = $path_exec.$exec[0] if (scalar @exec == 1);
	$exec = $path_exec.prompt("Choose a $exec_type version for project ".$projectName.':', -menu=>\@exec);
	$exec .= "/$exec_type";
	confess("'$exec' is not an executable") unless (-x $exec);
}
warn $exec;

qx/$exec --version/ =~ /^(?:cell|space)ranger-?(?:arc|atac)? ((?:(?:cell|space)ranger-?(?:arc|atac)?-)?(\d+\.\d+\.\d+))$/
	|| qx/realpath \$(which $exec)/ =~ /((?:cell|space)ranger-?(?:arc|atac)?-(\d+\.\d+\.\d+))/
	|| confess("Error determining the version");
my $version = $1;
my $version_nb = $2;
warn 'version: '.$version;
warn 'version nb: '.$version_nb;


#my $dir = $project->getProjectRootPath();
my $dir = $project->getCountingDir('cellranger');
$dir = $project->getCountingDir('spaceranger') if ($type eq 'spatial');
warn $dir;


my $hfamily;
my $prog;
foreach my $patient (@{$patients}) {
	my $hpatient;
	$hpatient->{name}= $patient->name();
	$hpatient->{group}= $patient->somatic_group();
	$hpatient->{pool}= $patient->family();
	$hpatient->{bc} = $patient->barcode;
	$hpatient->{bc2} = $patient->barcode2;
	$hpatient->{objet} = $patient;
	push(@{$hfamily->{$patient->barcode()}},$hpatient);
	$prog = $patient->alignmentMethod();
}



#------------------------------
# DRAGEN DEMULIPLEXAGE
#------------------------------
if (grep(/(dragen_)?demultiplex|^all$/i, @steps)){
	my $cmd_demultiplex = "$Bin/cellranger_samplesheet.pl -project=$projectName -mismatch=$mismatch -multi";
	$cmd_demultiplex .= "-no_exec " if ($no_exec);
	system($cmd_demultiplex);
}




#------------------------------
# DEMULIPLEXAGE OLD
#------------------------------
if (grep(/^demultiplex_old$/i, @steps)){
	
	my $bcl_dir = $run->bcl_dir;
	warn $bcl_dir;
	my $tmp = $project->getAlignmentPipelineDir("cellranger_demultiplex");
	warn $tmp;

	my $sampleSheet = $bcl_dir."/sampleSheet.csv";
	open (SAMPLESHEET,">$sampleSheet");
	print SAMPLESHEET "Lane,Sample,Index\n";
	
	for (my $i=1; $i<=$lane; $i++){
		foreach my $k (sort keys(%$hfamily)){
			my @pat = @{$hfamily->{$k}};
			print SAMPLESHEET $i.",".$pat[0]->{pool}.",".$k."\n";
		}
	}
	close(SAMPLESHEET);
	
	my $cmd = "cd $tmp; $Bin/../../demultiplex/demultiplex.pl -dir=$bcl_dir -run=$run_name -hiseq=10X -sample_sheet=$sampleSheet -cellranger_type=$exec -mismatch=$mismatch";
	warn $cmd;
	my $exit = system ($cmd) unless ($no_exec);
	exit($exit) if ($exit);
	unless ($@ or $no_exec) {
		print "\t------------------------------------------\n";
		print("Check the demultiplex stats\n");
		print("https://www.polyweb.fr/NGS/demultiplex/$run_name/$run_name\_laneBarcode.html\n\n");
		system ("firefox https://www.polyweb.fr/NGS/demultiplex/$run_name/$run_name\_laneBarcode.html &");
		print "\t------------------------------------------\n";
	}
}




#-----
# Choose transcriptome reference (COUNT/VELOCYTO)
#-----
# https://www.10xgenomics.com/support/software/cell-ranger/latest/release-notes/cr-reference-release-notes
my $htranscriptome;
if (grep(/count|^all$|velocyto/i, @steps) and $choose_transcriptome){
	my $transcriptome_dir = '/data-bipd/data-pure/public-data/10X/'.$project->genome_version.'/';
	if ($choose_transcriptome) {
		opendir(my $dh, $transcriptome_dir) || die "Can't opendir '$transcriptome_dir': $!";
		# gex/spatial
		my @transcriptomes = sort { -M $transcriptome_dir.$a <=> -M $transcriptome_dir.$b } grep {-d $transcriptome_dir.$_} readdir($dh);
		closedir ($dh);
		if (grep {/exp|nuclei|adt/i} @groups) {
			my @transcriptomes_gex = grep {/^refdata-(cellranger-|gex-)GRC/} @transcriptomes;
			confess("No transcriptome reference found in '$transcriptome_dir'") unless (scalar @transcriptomes_gex);
			$htranscriptome->{gex} = $transcriptome_dir.$transcriptomes_gex[0] if (scalar @transcriptomes_gex == 1);
			if (scalar @transcriptomes_gex > 1) {
				my $selected = prompt("Choose a transcritome reference for project ".$projectName.':', -menu=>\@transcriptomes_gex);
				die unless ($selected and -d $transcriptome_dir.$selected);
				$htranscriptome->{gex} = $transcriptome_dir.$selected;
			}
		}
		
		# vdj
		if (grep {/vdj/i} @groups) {
			my @transcriptomes_vdj = grep {/^refdata-cellranger-vdj-GRC/} @transcriptomes;
	#		my @transcriptomes = sort { -M $transcriptome_dir.$a <=> -M $transcriptome_dir.$b } grep {-d $transcriptome_dir.$_ && /^refdata-cellranger-vdj-GRC/} readdir($dh);
			confess("No V(D)J reference found in '$transcriptome_dir'") unless (scalar @transcriptomes_vdj);
			$htranscriptome->{vdj} = $transcriptome_dir.$transcriptomes_vdj[0] if (scalar @transcriptomes_vdj == 1);
			if (scalar @transcriptomes_vdj > 1) {
				my $selected = prompt("Choose a V(D)J reference for project ".$projectName.':', -menu=>\@transcriptomes_vdj);
				die unless ($selected and -d $transcriptome_dir.$selected);
				$htranscriptome->{vdj} = $transcriptome_dir.$selected;
			}
		}
		#warn Dumper $htranscriptome;
	}
}



#------------------------------
# COUNT
#------------------------------
if (grep(/count|^all$/i, @steps)){
	
	$create_bam = 'false' unless ($create_bam);
	$create_bam = 'true' if ($create_bam ne 'false');
	warn "create-bam=$create_bam";
	
	# Check chemistry if option provided
	 if ($chemistry) {
 		warn "chemistry=$chemistry";
		my @chemistries = qw/auto threeprime fiveprime 
		SC3Pv1 SC3Pv2 SC3Pv3 SC3Pv3-polyA SC3Pv4 SC3Pv4-polyA SC3Pv3HT SC3Pv3HT-polyA SC-FB SC3Pv3-polyA-OCM SC3Pv3-CS1-OCM SC3Pv4-polyA-OCM SC3Pv4-CS1-OCM 
		SC5P-PE SC5P-PE-v3 SC5P-R2 SC5P-R2-v3 SC5PHT SC5P-R1-OCM-v3 SC5P-R2-OCM SC5P-R2-OCM-v3 SC5P-PE-OCM-v3 SCVDJ-v3-OCM SCVDJ-R2-OCM-v3/;
		die ("Chemistry option '$chemistry' not valid, should be one of: ". join(', ', @chemistries)) unless (grep { $_ eq $chemistry } @chemistries);
	 }
	
	# Vérifie qu'il n'y ait pas de patient associé oublié
	my @pat_to_add;
	my %patients_id = map {$_->id => 1} @$patients;
	my %pat_to_add_id;
	foreach my $patient (@{$patients}) {
	    my $pname = $patient->name;
	    my $pfam = $patient->family;
	    my @oublis = grep {
	        ($_->name =~ /$pname/i or $_->family =~ /^$pfam$/i)
	        and $_->id() ne $patient->id
	        and !exists $patients_id{$_->id()}
	        and !exists $pat_to_add_id{$_->id()}
	    } @$all_patients;
#	    warn Dumper([map {$_->name} @oublis]) if @oublis;
	    foreach my $oubli (@oublis) {
	        if (prompt("'".$oubli->name."' in project but not in patients list. Do you want to add it? ", -yes_no, -1)) {
	            push(@pat_to_add, $oubli);
	            push(@groups,$oubli->somatic_group) unless (exists $pat_to_add_id{$oubli->id});
	            my $hpatient->{name}= $oubli->name();
				$hpatient->{group}= $oubli->somatic_group();
				$hpatient->{pool}= $oubli->family();
				$hpatient->{bc} = $oubli->barcode;
				$hpatient->{objet} = $oubli;
				push(@{$hfamily->{$oubli->barcode()}},$hpatient) unless (grep {$_->{objet}->id eq $oubli->id} @{$hfamily->{$patient->barcode()}});
	            $pat_to_add_id{$oubli->id} = 1;
	        }
	    }
	}
	push(@$patients, @pat_to_add) if @pat_to_add;
#	warn Dumper([map {$_->name} @$patients]);
	
	# Vérifie que les fastq sont bien nommés
#	foreach my $patient (@{$patients}) {
#		my $pname = $patient->name;
#		# todo: for flex multiplex
#		$pname = $patient->family if ($multi_meth =~ /^flex/);
#		warn 'fastq name: '.$pname;
#		my $fastq_files = $patient->fastqFiles();
#		my @fastq_files = map {values %$_} @$fastq_files;
#		confess ("Fastq file names (sample $pname) must follow the following naming convention: [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz"
#			 ."\n".Dumper \@fastq_files)
#			unless (scalar (grep {/$pname\_S\d*_L\d{3}_[IR][123]_\d{3,}\.fastq\.gz$/} @fastq_files) == scalar @fastq_files);
#	}
#	
#	foreach my $patient (@{$patients}) {
#	    my $pname = $patient->name;
#	    # todo: for flex multiplex
##	    $pname = $patient->family if ($multi_meth =~ /^flex/);
#	
#	    # Identifiants candidats pour un sample poolé (poolName et/ou barcode)
#	    my @candidates = grep { defined && length } ($patient->poolName, $patient->barcode, $patient->name);
#	    warn 'fastq name: '.$pname;
#	    my $fastq_files = $patient->fastqFiles();
#	    my @fastq_files = map {values %$_} @$fastq_files;
#	
#	    my $matched_id;
#	    for my $id (@candidates) {
#	        my $n_match = scalar grep {/^\Q$id\E_S\d*_L\d{3}_[IR][123]_\d{3,}\.fastq\.gz$/} @fastq_files;
#	        if ($n_match == scalar @fastq_files) {
#	            $matched_id = $id;
#	            last;
#	        }
#	    }
#	
#	    confess ("Fastq file names (sample $pname) must follow the following naming convention: [Sample Name or Pool Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz"
#	             ."\n".Dumper \@fastq_files)
#	        unless defined $matched_id;
#	}
	
	my $tmp = $project->getAlignmentPipelineDir("cellranger_multi");
	warn $tmp;
	
	# Si adt, vérifie feature_ref
	if (grep {$_ =~ /adt/i} @groups) {
		die("feature_ref csv required\n") unless ($feature_ref);
		die("'$feature_ref' not found") unless (-e $feature_ref);
		$feature_ref = abs_path($feature_ref);
	}

	
	# Choose the probe set
	my $probe_set_dir = "/data-bipd/data-pure/software/distrib/$exec/$exec-$version_nb/probe_sets/";
	$probe_set_dir = $exec =~ s/(cell|space)ranger(-(arc|atac))?$/probe_sets\//r unless (-e $probe_set_dir);
	my $probe_set;
	my $index = $project->getGenomeIndex($prog);
	my $htranscriptome;
	if ($multi_meth =~ /^flex/) {
		opendir(my $dh, $probe_set_dir) || die "Can't opendir '$probe_set_dir': $!";
		my @probe_sets = grep {-f $probe_set_dir.$_ && /^Chromium_\w*_Transcriptome_Probe_Set_v.*\.csv$/} readdir($dh);
		closedir ($dh);
		confess("No probe set found in $probe_set_dir") unless (scalar @probe_sets);
		
		my @sets;
		if ($project->getVersion() =~ /^HG38/) {
			@sets = grep{/^Chromium_Human_Transcriptome_Probe_Set_v[.0-9]*_GRCh38/} @probe_sets;
		}
		elsif ($project->getVersion() =~ /^MM38/) {
			@sets = grep{/^Chromium_Mouse_Transcriptome_Probe_Set_v[.0-9]*_mm10/} @probe_sets;
			
		}
		elsif ($project->getVersion() =~ /^MM39/) {
			@sets = grep{/^Chromium_Mouse_Transcriptome_Probe_Set_v[.0-9]*_GRCm39/} @probe_sets;
		}
		else {
			die("No probe set for release ".$project->getVersion().'. See https://www.10xgenomics.com/support/software/cell-ranger/downloads');
		}
		confess("No probe set found for release ".$project->getVersion()." in $probe_set_dir") unless (scalar @sets);
		$probe_set = $probe_set_dir.$sets[0] if (scalar @sets == 1);
		$probe_set = $probe_set_dir.prompt("Choose a probe set for project ".$projectName.':', -menu=>\@sets);
		die("No probe set '$probe_set'") unless (-f $probe_set);
	}
	
	
	# Vérifie que le pipeline n'ait pas déjà tourné : check si web_summary existe
	my @patients = @$patients;
	foreach my $patient (@patients) {
		my $pname = $patient->name;
		my $pool = $patient->family;
		my $group = $patient->somatic_group;
		next if ($group =~ /adt/i);
		my $file_out = "$dir$pool/$pool/" if ($multi_meth eq '5P_singleplex');
		$file_out = "$dir$pool/$pname/" if ($multi_meth eq '5P_singleplex');
		$file_out .= "web_summary.html" if ($group =~ /exp|nuclei/i);
		$file_out .= "vdj_.*/vloupe.vloupe" if ($group =~ /vdj/i);
#		warn $file_out;
		if (-e  $file_out and not $force) {
			warn "NEXT: $pname pipeline already completed: web_summary already exists" if ($group =~ /exp|nuclei/i);
			warn "NEXT: $pname pipeline already completed: vloupe already exists" if ($group =~ /vdj/i);
			@$patients = grep{$_->family ne $pool } @$patients;
			delete $hfamily->{$patient->barcode};
		}
	}
	undef @patients;
	confess("All done ! If you want to rerun samples, use -force") unless (scalar @$patients);
	
	# Vérifie qu'il n'y ait pas déjà un répertoire patient dans le $tmp
	foreach my $patient (@$patients) {
		my $pname = $patient->name;
		my $pool = $patient->family;
		my $group = $patient->somatic_group;
		next if ($group =~ /adt/i);
		if (-d $tmp.$pool) {
			die("'$tmp$pool/' already exists") unless(prompt("'$tmp$pool/' already exists. Continue ? ",-yes_no, -1));
		}
	}
	
	
	open (JOBS, ">$dir/jobs_multi.txt");
	foreach my $bc (keys(%$hfamily)){
		my @pat = @{$hfamily->{$bc}};
		next if (grep {$_->{group} =~ /adt/i} @pat);
		my $poolName = $pat[0]->{pool};
		my $name = $pat[0]->{name};
		if (-d $tmp.$poolName) {
			die("'$tmp$poolName/' already exists") unless (prompt("'$tmp$poolName' already exists.\nContinue anyway ? (y/n) ", -yes_no));
		}
		my @adt_flex = map { grep { $_->{group} =~ /adt/i } @{$hfamily->{$_}} } keys %$hfamily;
		warn scalar @adt_flex .' ADT found';
		#warn Dumper map {$_->{name}} @adt_flex;

		my $pobj=$pat[0]->{objet};
		my $seq_dir = $pobj->getSequencesDirectory();
		my $tmp_fastq = $tmp.'fastq/';
		make_path("$tmp_fastq", { mode => 0775 }) unless (-d $tmp_fastq);
		opendir(my $dh, $seq_dir) or die "Impossible d'ouvrir $seq_dir : $!";
		my @fastq = sort #map  { $seq_dir. $_}
		            grep { /^($poolName|$bc|$name).*\.fastq\.gz$/ }
		            readdir($dh);
		closedir($dh);
		$fastq[0] =~ /^($poolName|$bc|$name).*\.fastq\.gz$/;
		my $fastq_id = $1;
#		my $cmd_multi = "cp -u ".join(' ', map { $seq_dir.$_} @fastq )." $tmp_fastq && ";# unless ($multi_meth eq '5P_singleplex');
		my $cmd_multi = "cp $seq_dir$fastq_id*.fastq.gz $tmp_fastq && ";# unless ($multi_meth eq '5P_singleplex');
#		$cmd_multi = "cp -u $seq_dir$poolName*.fastq.gz $seq_dir".$adt_flex[0]->{pool}."*.fastq.gz $tmp_fastq && " if ($multi_meth ne '5P_singleplex' and scalar @adt_flex);
#		$cmd_multi = "cp -u ".join(' ', map{$seq_dir.$_->{name}.'*.fastq.gz'} @pat)." $tmp_fastq && " if ($multi_meth eq '5P_singleplex');
		
		my $config_csv = "$dir$poolName.csv";
		open (CSV,">$config_csv") or confess ("Can't open '$config_csv': $@");
		print CSV "[gene-expression]\n";
		my $transcriptome = $index;
		$transcriptome = qx/realpath $index/ if (-l $index);
		$transcriptome = $htranscriptome->{gex} if ($choose_transcriptome);
		chomp $transcriptome;
		print CSV "reference,".$transcriptome."\n";# unless ($multi_meth =~ /^flex/);
		print CSV "probe-set,".$probe_set."\n" if ($multi_meth =~ /^flex/);
		print CSV "create-bam,$create_bam\n";
		print CSV "chemistry,$chemistry\n" if ($chemistry);
		print CSV "\n";
		print CSV "[libraries]\n";
		print CSV "fastq_id,fastqs,feature_types\n";
		if ($multi_meth eq '5P_singleplex') {
			foreach my $p (sort { $a->{name} cmp $b->{name} } @pat){
				my $feature_type = 'Gene Expression' if ($p->{group} =~ /^exp$/i);
				$feature_type = 'VDJ' if ($p->{group} =~ /^vdj/i);
				$feature_type = 'Antibody Capture' if ($p->{group} =~ /^adt$/i);
				confess ("Error: could not attribute a feature type to sample ".$p->{name}." with group ".$p->{group}) unless ($feature_type);
				print CSV $p->{name}.','.$tmp_fastq.','.$feature_type."\n";
			}
		}
		else {
			print CSV $fastq_id.",".$tmp_fastq.",Gene Expression\n";
			print CSV $adt_flex[0]->{pool}.",".$tmp_fastq.",Antibody Capture\n" if (scalar @adt_flex);
		}
		print CSV "\n";
		if (grep {/vdj/i} @groups) {
			my $transcriptome = readlink $index."_vdj";
			$transcriptome = $htranscriptome->{vdj} if ($choose_transcriptome);
			print CSV "[vdj]\n";
			print CSV "reference,".$transcriptome."\n";
			print CSV "\n";
		}
		if (grep {/adt/i} @groups) {
			print CSV "[feature]\n";
			print CSV "reference,".$feature_ref."\n";
			print CSV "\n";
		}
		if ($multi_meth !~ /singleplex/) {
			print CSV "[samples]\n";
			print CSV "sample_id,probe_barcode_ids\n" if ($multi_meth =~ 'flex');
			print CSV "sample_id,ocm_barcode_ids\n" if ($multi_meth eq 'ocm');
			print CSV "sample_id,hashtag_ids\n" if ($multi_meth eq 'hastag');
			foreach my $p (sort { $a->{name} cmp $b->{name} } @pat){
				my $pname = $p->{name};
				confess ("$pname missing BC2") unless ($p->{bc2});
				print CSV $pname.','.$p->{bc2}."\n";
				if (scalar @adt_flex) {
					my ($adt) = grep {$_->{name} =~ /ADT/i and $_->{name} =~ /$pname/} @adt_flex;
					#print CSV $pname.','.$p->{bc2}.'+'.$adt->{bc2}."\n"; # Flex v1
				}
			}
		}
		close CSV;
		my $cmd_cellranger = "$exec multi --id=$poolName --csv=$config_csv --localcores=$cpu "; # --localcores=$cpu # --jobmode slurm
		warn $cmd_cellranger;
		unless ($no_exec) {
			foreach my $pat (map {$_->{objet}} @pat) {
				$pat->update_software_version('cellranger',$cmd_cellranger,$version_nb);
			}
		}
		$cmd_multi .= "cd $tmp && $cmd_cellranger";
		$cmd_multi .= "&& mkdir $dir$poolName --mode 775 && cp -r $tmp$poolName/outs/per_sample_outs/* $tmp$poolName/outs/qc_report.html $tmp$poolName/outs/config.csv $tmp$poolName/_versions $tmp$poolName/_cmdline $tmp$poolName/outs/config.csv $dir$poolName/ ";
		print JOBS $cmd_multi."\n";
	}
	close JOBS;
	
	
#	my $cmd2 = "cat $dir/jobs_multi.txt | $Bin/../../../../polyscripts/system_utility/run_cluster.pl -name cellranger -cpu=$cpu";
	my $cmd2 = "cat $dir/jobs_multi.txt | run_cluster.pl -cpu=$cpu";
	if (grep(/count|^all$/i, @steps)){
		warn $cmd2;
		sleep(5) unless ($no_exec);
		system ($cmd2) unless ($no_exec);
	}
	
	# Open web summaries
	my @error;
	my $web_summaries = "";
	foreach my $bc (sort keys %$hfamily) {
		my $pool = $hfamily->{$bc}->[0]->{pool};
		foreach my $hpatient (grep { $_->{group} =~ /exp|nuclei/i } @{$hfamily->{$bc}}) {
			my $file = $dir.$pool.'/'.$hpatient->{name}."/web_summary.html";
			$file = $dir.$pool.'/'.$pool."/web_summary.html" if ($multi_meth eq '5P_singleplex');
			$web_summaries .= $file.' ' if (-e $file);
			push(@error, $file) unless (-e $file or $no_exec);
		}
		unless($multi_meth eq '5P_singleplex' or $multi_meth eq 'flex_singleplex') {
			my $qc_report = $dir.$pool."/qc_report.html";
			$web_summaries .= $qc_report.' ' if (-e $qc_report);
			push(@error, $qc_report) unless (-e $qc_report or $no_exec);
		}
	}
	my $cmd3 = "firefox ".$web_summaries;
	$cmd3 = "google-chrome ".$web_summaries if (getpwuid($<) eq 'shanein');
	warn $cmd3 if ($web_summaries);
	system($cmd3.' &') if ($web_summaries and not $no_exec);
	die("Web summaries not found: ".join(', ', @error)) if (@error and not $no_exec);

	unless ($no_exec) {
		print "\t------------------------------------------\n";
		print("\tCheck the qc reports and/or the web summaries:\n");
		print("\t$dir*/qc_report.html $dir*/*/web_summary.html\n");
		print "\t------------------------------------------\n\n";
	}
	
}


#------------------------------
# AGGREGATION
#------------------------------
if (grep(/aggr/, @steps)){
	my $aggr_file = $dir."jobs_aggr.txt";
	my $id = $projectName.'_aggregation';
	$id = $aggr_name if $aggr_name;
	open (JOBS_AGGR, ">$aggr_file");
	my $type = $project->getRun->infosRun->{method};
	print JOBS_AGGR "sample_id,molecule_h5\n";
	foreach my $patient (sort {$a->name cmp $b->name} @$patients) {
		print JOBS_AGGR $patient->name().",".$dir."/".$patient->somatic_group.'/'.$patient->name."/molecule_info.h5\n";
	}
	close JOBS_AGGR;
	my $aggr_cmd = "cd $dir && $exec aggr --id=$id --csv=$aggr_file";
	warn $aggr_cmd;
	system ($aggr_cmd)  unless $no_exec;
}



#------------------------------
# COPY to SingleCell shared directory
#------------------------------
if (grep(/^cp(_web_summar(y|ies))?$|^all$/i, @steps)){
	my $cmd_cp = "$Bin/cellranger_copy.pl -project=$projectName ";
	$cmd_cp .= "-patients=$patients_name " if ($patients_name);
	$cmd_cp .= "-all_outs " if (grep(/^cp$/i, @steps));
	$cmd_cp .= "-no_exec " if ($no_exec);
	system($cmd_cp);
}



#------------------------------
# ARCHIVE / TAR
#------------------------------
if (grep(/tar|archive|^all$/i, @steps)){
	my $cmd_tar = "$Bin/cellranger_tar.pl -multi -project=$projectName ";
	$cmd_tar .= "-patients=$patients_name " if ($patients_name);
	$cmd_tar .= "-create_bam " if ($create_bam and $create_bam ne 'false');
	$cmd_tar .= "-no_exec " if ($no_exec);
	system($cmd_tar);
}



#------------------------------
# Velocyto
#------------------------------
if (grep(/velocyto/, @steps)) {
	my $cmd_velocyto = "$Bin/velocyto.pl -multi -project=$projectName ";
	$cmd_velocyto .= "-patients=$patients_name " if ($patients_name);
	$cmd_velocyto .= '-transcriptome='.$htranscriptome->{gex}.' '  if ($choose_transcriptome);
	$cmd_velocyto .= "-cpu=$cpu ";
	$cmd_velocyto .= "-no_exec " if ($no_exec);
	$cmd_velocyto .= "-version $project_version " if ($project_version);
	system($cmd_velocyto);
}



exit;



sub usage {
	print "
$0
-------------
Obligatoires:
	project <s>                nom du projet
	feature_ref	<s>            tableau des ADT, obligatoire seulement si step=count et qu'il y a des ADT
Optionels:
	steps <s>                  étape à réaliser: dragedemultiplex, demultiplex_old, count, tar, aggr, aggr_vdj, cp, cp_web_summary ou all (= demultiplex, count, cp_web_summary, tar)
	patients <s>               noms de patients/échantillons, séparés par des virgules
	cpu <i>                    nombre de cpu à utiliser, défaut: 20
	lane <i>                   nombre de lanes sur la flowcell, défaut: lit le RunInfo.xml
	mismatches <i>             nombre de mismatches autorisés lors du démultiplexage, défaut: 0
	create-bam/nocreate-bam    générer ou non les bams lors du count, défaut: nocreate-bam
	aggr_name <s>              nom de l'aggrégation, lors de step=aggr ou aggr_vdj
	chemistry                  chemistry , défaut: auto (pour librairies exp et adt)
	no_exec                    ne pas exécuter les commandes
	help                       affiche ce message

";
	exit(1);
}


