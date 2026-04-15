#!/usr/bin/perl

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
use Text::CSV qw( csv );
use Time::Local 'timelocal';
use File::Temp qw/ tempfile tempdir /;
use Term::Menus;
use Proc::Simple;
use Storable;
use JSON::XS;
use XML::Simple qw(:strict);
use Cwd 'abs_path';
use File::Path qw(make_path);
use Carp;


my $projectName;
my $patients_name;
my @steps;
#my $lane;
my $mismatch = 1;
my $feature_ref;
my $cmo_ref;
my $no_exec;
my $aggr_name;
my $choose_exec;
my $choose_transcriptome;
my $chemistry;
my $create_bam;
my $slide_infos_csv;
my $probe_set;
my $add_image;
my $loupe_alignment;
my $cpu = 20;
my $force;
my $help;

GetOptions(
	'project=s'										=> \$projectName,
	'patients=s'									=> \$patients_name,
	'steps=s{1,}'									=> \@steps,
#	'lane|nb_lane=i'								=> \$lane,
	'mismatches=i'									=> \$mismatch,
	'create_bam!'									=> \$create_bam,
	'feature_ref|feature_csv=s'						=> \$feature_ref,
	'cmo_ref|cmo_csv=s'								=> \$cmo_ref,
	'aggr_name=s'									=> \$aggr_name,
	'choose_exec|choose_version|version'			=> \$choose_exec,
	'choose_transcriptome|select_transcriptome'		=> \$choose_transcriptome,
	'chemistry=s'									=> \$chemistry,
	'slide_infos=s'									=> \$slide_infos_csv,
	'probe_set=s'									=> \$probe_set,
	'add_image'										=> \$add_image,
	'loupe_alignment|json'							=> \$loupe_alignment,
	'cpu=i'											=> \$cpu,
	'no_exec'										=> \$no_exec,
	'force'											=> \$force,
	'help'											=> \$help,
) || confess("Error in command line arguments\n");

usage() if $help;
die("-project argument is mandatory") unless ($projectName);
die("cpu must be in [1;40], given $cpu") unless ($cpu > 0 and $cpu <= 40);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
my $all_patients = $project->getPatients;
my $patients = $project->get_only_list_patients($patients_name);
die("No patient in project $projectName") unless ($patients);

# Vérifie les caractères non acceptés
my @patient_names = map {$_->name} @$patients;
my @invalid_names = grep { $_ !~ /^[A-Za-z0-9_-]+$/ } @patient_names;
die ("Patient names can only contain letters, numbers, hyphens and underscores. Sample names ".join(@invalid_names, ', ')." are invalids.") if (@invalid_names);

# Vérifie que le projet est en somatic et que les groupes sont correctement remplis
unless ($project->isSomatic) {
	my $warn = "Error: Project $projectName is not in somatic mode. "
		."Activate somatic mode and check that the groups have been filled in, so that they can be taken into account in the analysis.";
#	die ($warn) if (grep{/cont|aggr/} @steps);
	die ($warn);
}
my @groups = map {$_->somatic_group} @$patients;
die("Error: Check that the groups are not empty and correctly filled. Accepted groups are exp, adt, vdj, atac, spatial, cmo, arc") 
	unless (grep{/^(exp|adt|vdj|atac|spatial|cmo|arc)$/i} @groups);


# Steps
@steps = split(/,/, join(',',@steps));
my $list_steps = ['dragen_demultiplex', 'teleport', 'count', 'aggr', 'aggr_vdj', 'tar', 'cp', 'cp_web_summary', 'infos', 'velocyto', 'all'];
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

$create_bam = 1 if (grep(/velocyto/i, @steps) or grep (/velocyto/, map($_->getCallingMethods, @$patients)));


my $run = $project->getRun();
my $run_name = $run->plateform_run_name;
my $type = $run->infosRun->{method};
my $machine = $run->infosRun->{machine};

# Executable 
my $exec = "cellranger";
$exec .= '-atac' if ($type eq 'atac');
$exec .= '-arc' if ($type eq 'arc');
$exec = 'spaceranger' if  ($type eq 'spatial');
if ($choose_exec) {
	my $exec_type = $exec;
	my $path_exec = "/software/distrib/$exec_type/";
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

my $dir_project = $project->getProjectRootPath();
my $dir = $project->getCountingDir('cellranger');
$dir = $project->getCountingDir('spaceranger') if ($type eq 'spatial');
warn $dir;



#-----
# Choose transcriptome reference (COUNT/VELOCYTO)
#-----
# https://www.10xgenomics.com/support/software/cell-ranger/latest/release-notes/cr-reference-release-notes
my $htranscriptome;
if (grep(/count|velocyto/, @steps)) {
	if (grep(/count|^all$|velocyto/i, @steps) and $choose_transcriptome){
		my $transcriptome_dir = '/data-isilon/public-data/10X/'.$project->genome_version.'/';
		opendir(my $dh, $transcriptome_dir) || die "Can't opendir '$transcriptome_dir': $!";
		# gex/spatial
		my @transcriptomes = sort { -M $transcriptome_dir.$a <=> -M $transcriptome_dir.$b } grep {-d $transcriptome_dir.$_} readdir($dh);
		closedir ($dh);
		if (grep {/exp|adt/i} @groups) {
			my @transcriptomes_gex = grep {/^refdata-(cellranger-|gex-)GRC/} @transcriptomes;
			confess("No transcriptome reference found in '$transcriptome_dir'") unless (scalar @transcriptomes_gex);
			warn Dumper \@transcriptomes_gex;
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
			confess("No V(D)J reference found in '$transcriptome_dir'") unless (scalar @transcriptomes_vdj);
			$htranscriptome->{vdj} = $transcriptome_dir.$transcriptomes_vdj[0] if (scalar @transcriptomes_vdj == 1);
			if (scalar @transcriptomes_vdj > 1) {
				my $selected = prompt("Choose a V(D)J reference for project ".$projectName.':', -menu=>\@transcriptomes_vdj);
				die unless ($selected and -d $transcriptome_dir.$selected);
				$htranscriptome->{vdj} = $transcriptome_dir.$selected;
			}
		}
		
		# atac
		if (grep {/atac/i} @groups) {
			my @transcriptomes_atac = grep {/^refdata-cellranger-(arc|atac)-GRC/} @transcriptomes;
			confess("No ATAC reference found in '$transcriptome_dir'") unless (scalar @transcriptomes_atac);
			$htranscriptome->{atac} = $transcriptome_dir.$transcriptomes_atac[0] if (scalar @transcriptomes_atac == 1);
			if (scalar @transcriptomes_atac > 1) {
				my $selected = prompt("Choose a ATAC reference for project ".$projectName.':', -menu=>\@transcriptomes_atac);
				die unless ($selected and -d $transcriptome_dir.$selected);
				$htranscriptome->{atac} = $transcriptome_dir.$selected;
			}
		}
		warn Dumper $htranscriptome;
	}
}


	
#------------------------------
# DRAGEN DEMULIPLEXAGE
#------------------------------
if (grep(/dragen_demultiplex|^all$/i, @steps)){
	my $cmd_demultiplex = "$Bin/cellranger_samplesheet.pl -project=$projectName -mismatch=$mismatch ";
	$cmd_demultiplex .= "-no_exec " if ($no_exec);
	system($cmd_demultiplex);
}



#------------------------------
# OLD DEMULIPLEXAGE
#------------------------------
if (grep(/^demultiplex_old$/i, @steps)){
	
	die("-mismatch can be 0, 1, 2. Two entries, comma delimited, allowed for dual index.") unless ($mismatch =~ /[012](,[012])?/);
	warn 'mismatch(es)=$mismatch';
	
	my $bcl_dir = $run->bcl_dir;
	warn $bcl_dir;
	my $fastq_dir = $run->fastq_dir;
	warn $fastq_dir;
#	my $config = XMLin("$bcl_dir/RunInfo.xml", KeyAttr => { reads => 'Reads' }, ForceArray => [ 'reads', 'read' ]);
#	my $lane_config = $config->{Run}->{FlowcellLayout}->{LaneCount};
#	if ($lane and $lane != $lane_config) {
#		die unless prompt("$lane_config lanes found in the RunInfo.xml. You entered -lane=$lane. Continue anyway ? (y/n)  ", -yes_no);
#	}
#	unless ($lane) {
#		$lane = $lane_config;
#		warn 'LaneCount=$lane';
#	}
	
	my $sampleSheet = $bcl_dir."/sampleSheet.csv";
	open (SAMPLESHEET,">$sampleSheet");
	print SAMPLESHEET "Lane,Sample,Index\n";
	
	foreach my $patient (sort {$a->name cmp $b->name} @{$patients}) {
		my $name = $patient->name();
		my $bc = $patient->barcode();
		print SAMPLESHEET "*,$name,$bc"."\n";
#		for (my $i=1;  $i<=$lane;$i++){
#			print SAMPLESHEET $i.",".$name.",".$bc."\n";
#		}
	}
	close(SAMPLESHEET);
	
	my $cmd = "$Bin/../../demultiplex/demultiplex.pl -dir=$bcl_dir -fastq_dir=$fastq_dir -run=$run_name -hiseq=10X -sample_sheet=$sampleSheet -cellranger_type=$exec -mismatch=$mismatch";
	warn $cmd;
	my $exit = system ($cmd) unless ($no_exec);
	exit($exit) if ($exit);
	unless ($@ or $no_exec) {
		print "\t------------------------------------------\n";
		print("Check the demultiplex stats\n");
		print("https://www.polyweb.fr/NGS/demultiplex/$run_name/$run_name\_laneBarcode.html\n");
		system("firefox https://www.polyweb.fr/NGS/demultiplex/$run_name/$run_name\_laneBarcode.html &");
		print "\t------------------------------------------\n";
	}
}



#------------------------------
# TELEPORT
#------------------------------
if (grep(/teleport/, @steps)) {
	my $cmd = "$Bin/../../../teleport_patients.pl -project=$projectName -force=1";
	system($cmd) unless $no_exec;
}



#------------------------------
# COUNT
#------------------------------
if (grep(/count|^all$/i, @steps)){
	
	if (grep {$_ =~ /adt/i} @groups) {
		die("feature_ref csv required\n") unless ($feature_ref);
		die("'$feature_ref' not found") unless (-e $feature_ref);
		$feature_ref = abs_path($feature_ref);
	}
	if (grep {$_ =~ /cmo/i} @groups) {
		die("cmo_ref csv required\n") unless ($cmo_ref);
		die("'$cmo_ref' not found") unless (-e $cmo_ref);
		$feature_ref = abs_path($feature_ref);
	}
	
	$create_bam = 'false' unless ($create_bam);
	$create_bam = 'true' if ($create_bam ne 'false');
	warn "create-bam=$create_bam";
	
	# $exec version
	qx/$exec --version/ =~ /^(?:cell|space)ranger-?(?:arc|atac)? ((?:(?:cell|space)ranger-?(?:arc|atac)?-)?(\d+\.\d+\.\d+))$/
		|| qx/realpath \$(which $exec)/ =~ /((?:cell|space)ranger-?(?:arc|atac)?-(\d+\.\d+\.\d+))/
		|| confess("Error determining the version");
	my $version = $1;
	my $version_nb = $2;
	warn 'version: '.$version;
	warn 'version nb: '.$version_nb;
	
	# Check chemistry if option provided
	 if ($chemistry) {
 		warn "chemistry=$chemistry";
		my @chemistries = qw/auto threeprime fiveprime SC3Pv2 SC3Pv3 SC3Pv3-polyA SC3Pv4 SC3Pv4-polyA SC3Pv3HT SC5P-PE SC5P-PE-v3 SC5P-R2 SC3Pv1 ARC-v1/;
		die ("Chemistry option '$chemistry' not valid, should be one of: ". join(', ', @chemistries)) unless (grep { $_ eq $chemistry } @chemistries);
	 }
	
	my $prog =  $patients->[0]->alignmentMethod();
	my $index = $project->getGenomeIndex($prog);
	
	# Vérifie que les fastq sont bien nommés
	foreach my $patient (@{$patients}) {
		my $pname = $patient->name;
		my $fastq_files = $patient->fastqFiles();
		my @fastq_files = map {values %$_} @$fastq_files;
		confess ("Fastq file names (sample $pname) must follow the following naming convention: [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz"
			 ."\n".Dumper \@fastq_files)
			unless (scalar (grep {/$pname\_S\d*_L\d{3}_[IR][123]_\d{3,}\.fastq\.gz$/} @fastq_files) == scalar @fastq_files);
	}
	
	# Vérifie que le pipeline n'ait pas déjà tourné : check si web_summary existe
	my @patients = @$patients;
	foreach my $patient (@patients) {
		my $pname = $patient->name;
		my $group = $patient->somatic_group;
		next if ($group =~ /adt/i);
		my $file_out = "$dir$pname/web_summary.html"; # if ($group =~ /exp|nuclei/i);
#		$file_out = "$dir$pname/vloupe.vloupe" if ($group =~ /vdj/i);
		if (-e  $file_out and not $force) {
			warn "NEXT: $pname pipeline already completed: web_summary already exists"; # if ($group =~ /exp|nuclei|atac|arc/i);
#			warn "NEXT: $pname pipeline already completed: vloupe already exists" if ($group =~ /vdj/i);
			@$patients = grep{ $_->name ne $pname } @$patients;
		}
	}
	undef @patients;
	confess("All done ! If you want to rerun samples, use -force") unless (scalar @$patients);
	
	my $tmp = $project->getAlignmentPipelineDir("cellranger");
	$tmp = $project->getAlignmentPipelineDir("spaceranger") if ($type eq 'spatial');;
	warn $tmp;

	# Vérifie qu'il n'y ait pas déjà un répertoire patient dans le $tmp
	foreach my $patient (@$patients) {
		my $pname = $patient->name;
		my $pool = $patient->family;
		my $group = $patient->somatic_group;
		next if ($group =~ /adt/i);
		if (-d $tmp.$pname) {
			die("'$tmp$pname/' already exists") unless(prompt("'$tmp$pname/' already exists. Continue ? ",-yes_no, -1));
		}
	}

	
	# Vérifie qu'il n'y ait pas des adt oubliés
	my @pat_to_add;
	my %patients_id = map {$_->id => 1} @$patients;
	my %pat_to_add_id;
	foreach my $patient (@{$patients}) {
	    my $pname = $patient->name;
	    my $pfam = $patient->family;
	    my @oublis = grep {
	        ($_->name =~ /$pname/i or $_->family =~ /^$pfam$/i)
	        and $_->somatic_group =~ /^(adt|exp)$/i
	        and $_->id() ne $patient->id
	        and !exists $patients_id{$_->id()}
	        and !exists $pat_to_add_id{$_->id()}
	    } @$all_patients;
#	    warn Dumper([map {$_->name} @oublis]) if @oublis;
	    foreach my $oubli (@oublis) {
	        if (prompt("'".$oubli->name."' in project but not in patients list. Do you want to add it? ", -yes_no, -1)) {
	            push(@pat_to_add, $oubli);
	            $pat_to_add_id{$oubli->id} = 1;
	            push(@groups,$oubli->somatic_group) unless (exists $pat_to_add_id{$oubli->id});
	        }
	    }
	}
	push(@$patients, @pat_to_add) if @pat_to_add;
#	warn Dumper([map {$_->name} @$patients]);

	
	sub full_cmd {
		my ($pat, $cmd) = @_;
		chomp $cmd;
		my $name = $pat->name;
		my $tmp_fastq = $tmp.'fastq/';
		make_path($tmp_fastq, { mode => 0775 }) unless (-d "$tmp_fastq");
		my $seq_dir = $pat->getSequencesDirectory;
		my $cmd1 = "cp $seq_dir$name*.fastq.gz $tmp_fastq && ";
		my $adt = 'ADT_'.$name;
		$cmd1 = "cp $seq_dir$name*.fastq.gz $seq_dir/ADT_$name*.fastq.gz $tmp_fastq && " if (grep (/$adt/, @patient_names));
#		my $cmd2 = " --localcores=$cpu ";
		my $cmd2 = "&& mkdir $dir$name --mode=775 " unless (-d $dir.$name);
		$cmd2 .= "&& cp -r $tmp$name/outs/* $tmp$name/_versions $tmp$name/_cmdline $dir$name/ ";
#		make_path("$dir$name", { mode => 0775 }) unless (-d "$dir$name");
#		$cmd2 .= "&& rm $dir$name/possorted_bam.bam $dir$name/possorted_bam.bam.bai" if ($exec eq 'cellranger-atac' and $create_bam eq 'false');
		$pat->update_software_version($exec, $cmd, $version_nb) unless ($no_exec);
		$cmd =~ s/$seq_dir/$tmp_fastq/;
		return $cmd1.$cmd.$cmd2."\n";
	}
	
	# EXP
	my $type_exp = 1 if map {uc($_) =~ /(EXP|NUCLEI)/ } @groups;
	system("rm $dir/jobs_count.txt") if (-f "$dir/jobs_count.txt");
	if($type_exp){
		open (JOBS_EXP, ">$dir/jobs_count.txt");
		my @exp = grep { uc($_->somatic_group()) eq "EXP" || uc($_->somatic_group()) eq "NUCLEI" } @{$patients};
		warn "EXP/NUCLEI: ".join(',',map($_->name,@exp));
		foreach my $e (@exp){
			my $name = $e->name;
#			warn 'EXP: '.$name;
			my $group = $e->somatic_group();
			my $fastq = $e->getSequencesDirectory;
			my $transcriptome = $index;
			$transcriptome = (qx/realpath $index/) if (-l $index);
			chomp $transcriptome;
			$transcriptome = '/data-isilon/public-data/10X/HG38/refdata-gex-GRCh38-2020-A' if ($project->description =~ /at?traction/i);
			$transcriptome = $htranscriptome->{gex} if ($choose_transcriptome);
			my $cmd = "cd $tmp && $exec count --id=$name --sample=$name --fastqs=$fastq --create-bam=$create_bam --transcriptome=$transcriptome ";
			$cmd .= " --include-introns true " if ($type eq "nuclei" or lc($group) eq "nuclei"); # true par défaut
			$cmd .= " --chemistry $chemistry " if ($chemistry);
			$cmd .= "\n";
			$cmd = full_cmd($e, $cmd);
			my @cmd_cellranger = grep {/^$exec/} split (/ +&& +/, $cmd);
			warn $cmd_cellranger[0];
	#		warn $cmd;
			print JOBS_EXP $cmd;
		}
	}
	
	
	# ADT
	my $type_adt = 1 if map {uc($_) =~ /ADT/ } @groups;
	if ($type_adt){
		open (JOBS_ADT, ">$dir/jobs_count.txt");
		my @exp = grep { uc($_->somatic_group()) eq "EXP"} @{$patients};
		warn "EXP+ADT: ".join(',',map($_->name,@exp));
		foreach my $e(@exp){
			my $name = $e->name();
#			warn 'ADT: '.$ename;
			my $efam = $e->family();
#			my $egroup = uc($e->somatic_group());
			my $lib_file = $dir.$name."_library.csv";
			open(LIB,">$lib_file") or die "impossible $lib_file";
			my @adt = grep {$_->family() eq $efam && $_->name() =~ $name && uc($_->somatic_group()) eq "ADT"} @{$patients};
			warn ("no associated ADT library") if scalar(@adt)==0 ;
			next() if scalar(@adt)==0 ; 
			my $adt_name = $adt[0]->name();
			my $lib = "fastqs,sample,library_type\n".$tmp.'fastq/,'.$name.",Gene Expression\n";
			$lib .= $tmp.'fastq/'.",".$adt_name.",Antibody Capture\n";
			print LIB $lib;
			close(LIB);
			my $transcriptome = $index;
			$transcriptome = qx/realpath $index/ if (-l $index);
			chomp $transcriptome;
			$transcriptome = $htranscriptome->{gex} if ($choose_transcriptome);
			my $cmd = "cd $tmp && $exec count --id=$name --feature-ref=$feature_ref --transcriptome=$transcriptome  --libraries=$lib_file --create-bam=$create_bam ";
			$cmd .= " --chemistry $chemistry " if ($chemistry);
			$cmd .= "\n";
			$cmd = full_cmd($e, $cmd);
			my @cmd_cellranger = grep {/^$exec/} split (/ +&& +/, $cmd);
			warn $cmd_cellranger[0];
			print JOBS_ADT $cmd;
		}
	}

	
	# VDJ
	my $type_vdj = 1 if map {uc($_) =~ /VDJ/ } @groups;
	system("rm $dir/jobs_vdj.txt") if (-f "$dir/jobs_vdj.txt");
	if($type_vdj){
		open (JOBS_VDJ, ">$dir/jobs_vdj.txt");
		my @vdj = grep { uc($_->somatic_group()) eq "VDJ"} @{$patients};
		warn "VDJ: ".join(',',map($_->name,@vdj));
		foreach my $v (@vdj){
			my $name = $v->name(); 
#			my $vfam = $v->family();
#			my $vgroup = uc($v->somatic_group());
			my $fastq = $v->getSequencesDirectory;
			my $transcriptome = $index.'_vdj';
			$transcriptome = qx/realpath $index"_vdj"/ if (-l $index);
			chomp $transcriptome;
			$transcriptome = $htranscriptome->{vdj} if ($choose_transcriptome);
			my $cmd = "cd $tmp && cellranger vdj --sample=$name --id=$name --fastqs=$fastq --reference=$transcriptome ";
			$cmd = full_cmd($v, $cmd);
			my @cmd_cellranger = grep {/^$exec/} split (/ +&& +/, $cmd);
			warn $cmd_cellranger[0];
			print JOBS_VDJ $cmd;
		}
	}


	# SPATIAL
	my $type_spatial = 1 if map {uc($_) =~ /SPATIAL/ } @groups;
	system("rm $dir/jobs_spatial.txt") if (-f "$dir/jobs_spatial.txt");
 	if($type_spatial){
 		confess("Error: Project release must be HG38") if ($project->getVersion() =~ /^HG19/);
		open (JOBS_SPATIAL, ">$dir/jobs_spatial.txt");
		my @spatial = grep { uc($_->somatic_group()) eq "SPATIAL"} @{$patients};
		warn "SPATIAL: ".join(',',map($_->name,@spatial));
 		my $choice_image_type = {
 			'Brightfield image generated by the CytAssist instrument (the CytAssist image)' => 'cytaimage',
 			'Brightfield microscope image' => 'image',
 			'Dark background fluorescence microscope image' => 'darkimage',
 			'Composite colored fluorescence microscope image' => 'colorizedimage',
		};
		my $imagetype = prompt('Choose the type of image you have: ', -m=>$choice_image_type);
		die ("No image type selected") unless ($imagetype);
		my $imagetype2;
		if ($add_image) {
			my $choice_second_image_type = {
			    map { $_ => $choice_image_type->{$_} }
			    grep { $choice_image_type->{$_} ne $imagetype }
			    keys %$choice_image_type
			};
			$imagetype2 = prompt('Choose the second type of image you have: ', -m=>$choice_second_image_type);
			die ("No image type selected") unless ($imagetype2);
		}
		
		my $slide_infos;
		if ($slide_infos_csv and -f $slide_infos_csv) {
			my $a_slide_infos = csv (in => $slide_infos_csv, headers => "auto", filter => "not_empty") or confess("Can't open '$slide_infos_csv': $!") unless ($imagetype ne 'cytaimage');
	 		foreach my $line (@$a_slide_infos) {
				$slide_infos->{$line->{'sample'}}->{'slide'} = $line->{'slide id'};
				$slide_infos->{$line->{'sample'}}->{'area'} = $line->{'aera'};
			}
		}
		foreach my $s (@spatial){
			my $sname = $s->name();
			
			# Slide id and area
			my ($slide_id, $area);
			unless ($imagetype eq 'cytaimage') {
				$slide_id = $slide_infos->{$sname}->{slide} if (exists $slide_infos->{$sname} and exists $slide_infos->{$sname}->{slide});
		 		$slide_id = prompt("$sname Visium Slide ID: ") unless ($slide_id);
		 		confess("Error: Visuim slide ID must start with V1, V4, V5 or H1") unless ($slide_id =~ /^(V[145]|H1)-/);
		 		confess("Error: Visuim Spatial slide v1 and Visuim CytAssist Spatial slide v2 must be in MM38") if ($slide_id =~ /^V[145]/ and $project->getVersion() eq 'MM39');
				$area = $slide_infos->{$sname}->{area} if (exists $slide_infos->{$sname} and exists $slide_infos->{$sname}->{area});
		 		$area = uc(prompt("$sname area for Visium Slide $slide_id: ")) unless ($area);
		 		confess("Error: Area for Visuim Spatial slide, v1, 6.5 mm (V1) can be one of A1, B1, C1, D1") if ($slide_id =~ /^V1/ and grep {$_ ne $area} qw{A1 B1 C1 D1});
		 		confess("Error: Area for Visuim CytAssist Spatial slide, v2, 11 mm (V5) can A or B") if ($slide_id =~ /^V5/ and grep {$_ ne $area} qw{A B});
		 		confess("Error: Area for Visuim CytAssist Spatial slide, v2, 6.5 mm (V4), and Visium HD slides, 6.5 mm (H1) can be one of A1, D1") if ($slide_id =~ /^(V4|H1)/ and grep {$_ ne $area} qw{A1 D1});
				#my ($slide,$slide2,$area) = split("-",$bc2);
				confess("slide id and area are required") if (not $area or not $slide_id);
			}
		
			# Image
			opendir(my $dh, $dir_project) or confess ("Can't opendir '$dir_project': $!");
			my @all_images = grep {/$sname.*\.(tif{1,2}|jpe?g)$/} readdir($dh);
			close($dh);
			if (scalar @all_images < 1) {confess("No tiff/jpg image found in '$dir_project'")}
			
			my @images = @all_images;
			@images = grep {/$sname.*\.tif{1,2}$/} @all_images if ($imagetype eq 'cytaimage');
			if (scalar @images < 1) {confess("No tiff image found in '$dir_project'")}
			
			my $image1;
			if (scalar @images == 1) {$image1 = $images[0]}
			elsif (scalar @images > 1) {$image1 = prompt("Choose an image file:", -menu=>\@images)}
#			warn "$imagetype = $image1";
			confess("No image") unless($image1);
			$image1 = $dir_project.$image1;
			
			my $image2;
			if ($add_image) {
				@images = grep {$image1 !~ /$_$/} @all_images;
				if (scalar @images < 1) {confess("No other tiff/jpg image found in '$dir_project'")}
				@images = grep {$_ =~ /$sname.*\.tif{1,2}$/} @images if ($imagetype2 eq 'cytaimage');
				if (scalar @images < 1) {confess("No other tiff image found in '$dir_project'")}
				elsif (scalar @images == 1) {$image2 = $images[0]}
				elsif (scalar @images > 1) {$image2 = prompt("Choose an image file:", -menu=>\@images)}
	#			warn "$imagetype2 = $image2";
				confess("No image") unless($image2);
				$image2 = $dir_project.$image2;
			}
			
			# Image alignment
			my $json;
			if ($loupe_alignment) {
				opendir(my $dh, $dir_project) or confess ("Can't opendir '$dir_project': $!");
				my @json = grep {/$sname.*\.json$/} readdir($dh);
				close($dh);
				if (scalar @json < 1) {confess("No json loupe alignment file found in '$dir_project'")}
				elsif (scalar @json == 1) {$json = $json[0]}
				elsif (scalar @json > 1) {$json = prompt("Choose an image file:", -menu=>\@json)}
	#			warn "loupe-alignment = $json";
				confess("No json") unless($json);
				$json = $dir_project.$json;
				
			}
			
			my $fastq = $s->getSequencesDirectory;
			my $transcriptome = $index;
			$transcriptome = qx/realpath $index/ if (-l $index);
			chomp $transcriptome;
			$transcriptome = $htranscriptome->{gex} if ($choose_transcriptome);
			
			unless($probe_set and -f $probe_set) {
				$probe_set = "/software/distrib/$exec/$exec-$version_nb/probe_sets/";
				opendir(my $dh, $probe_set) || die "Can't opendir '$probe_set': $!";
				my @probe_sets;
				@probe_sets = sort grep {-f $probe_set.$_ && /^Visium_\w+_Transcriptome_Probe_Set_v.+\.csv$/} readdir($dh);
				close($dh);
				@probe_sets = grep {/^Visium_Human_Transcriptome_Probe_Set_v[0-9.]+_GRCh38-20\d{2}-A\.csv$/} @probe_sets if ($project->getVersion() =~ /^HG/);
				@probe_sets = grep {/^Visium_Mouse_Transcriptome_Probe_Set_v[0-9.]+_mm10-20\d{2}-A\.csv$/} @probe_sets if ($project->getVersion() =~ /^MM38/);
				@probe_sets = grep {/^Visium_Mouse_Transcriptome_Probe_Set_v[0-9.]+_GRCm39-20\d{2}-A\.csv$/} @probe_sets if ($project->getVersion() =~ /^MM39/);
				#Visium_Mouse_Transcriptome_Probe_Set_v2.1.0_GRCm39-2024-A.csv
				my $probe_set_compatibility = "";
				if ($slide_id) {
					@probe_sets = grep {/^Visium_Human_Transcriptome_Probe_Set_v1[0-9.]+_GRCh38-20\d{2}-A.*\.csv$/} @probe_sets if ($project->getVersion() =~ /^HG/ && $slide_id =~ /^(V1)/);
					@probe_sets = grep {/^Visium_Human_Transcriptome_Probe_Set_v2[0-9.]+_GRCh38-20\d{2}-A.*\.csv$/} @probe_sets if ($project->getVersion() =~ /^HG/ && $slide_id =~ /^(V[45]|H1)/);
					@probe_sets = grep {/^Visium_Mouse_Transcriptome_Probe_Set_v1[0-9.]+_mm10-20\d{2}-A.*\.csv$/} @probe_sets if ($project->getVersion() =~ /^MM38/ && $slide_id =~ /^(V[145])/);
					@probe_sets = grep {/^Visium_Mouse_Transcriptome_Probe_Set_v2[0-9.]+_mm10-20\d{2}-A.*\.csv$/} @probe_sets if ($project->getVersion() =~ /^MM38/ && $slide_id =~ /^H1/);
					@probe_sets = grep {/^Visium_Mouse_Transcriptome_Probe_Set_v2[0-9.]+_GRCm39-20\d{2}-A.*\.csv$/} @probe_sets if ($project->getVersion() =~ /^MM39/ && $slide_id =~ /^H1/);
				}
				else {
					$probe_set_compatibility = qq{Probe set	Compatible Assay
	Human Probe Set v2	Visium HD Spatial Gene Expression
	                 	Visium CytAssist Spatial Gene and Protein Expression
	                  	Visium CytAssist Spatial Gene Expression (FFPE, Fresh Frozen)
	Human Probe Set v1	Visium Spatial Gene Expression for FFPE\n} if ($project->getVersion() =~ /^HG/);
					$probe_set_compatibility = qq{Mouse Probe Set v2	Visium HD Spatial Gene Expression
	Mouse Probe Set v1	Visium CytAssist Spatial Gene Expression (FFPE, Fresh Frozen, Fixed Frozen)
	                  	Visium Spatial Gene Expression for FFPE\n} if ($project->getVersion() =~ /^MM/);
				}
				confess("No probe set found at '$probe_set'") unless (scalar @probe_sets);
				$probe_set .= $probe_sets[0] if (scalar @probe_sets == 1);
				if (scalar @probe_sets > 1) {
					print("Probe set compatibility:\n".$probe_set_compatibility) unless ($slide_id);
					$probe_set .= prompt("Choose a probe set for project ".$projectName.':', -menu=>\@probe_sets);
				}
				confess("'$probe_set' does not exist") unless (-f $probe_set);
			}
			
			my $cmd = "cd $tmp && $exec count --id=$sname --sample=$sname --fastqs=$fastq --transcriptome=$transcriptome --create-bam=$create_bam ";
			$cmd .= " --probe-set=$probe_set ";
			$cmd .= "--$imagetype=$image1 ";
			$cmd .= "--$imagetype2=$image2 " if ($add_image);
			$cmd .= "--area=$area --slide=$slide_id " if ($area and $slide_id);
			$cmd .= " --loupe-alignment=$json " if ($loupe_alignment);
			$cmd = full_cmd($s, $cmd);
			my @cmd_cellranger = grep {/^spaceranger/} split (/ +&& +/, $cmd);
			warn $cmd_cellranger[0];
			print JOBS_SPATIAL $cmd;
			
# The slide serial number is the unique identifier printed on the label of each Visium slide. The serial number starts with a prefix indicating the slide type:
# V1: Visium Spatial Gene Expression Slide (v1, 6.5 mm)
# V4: Visium CytAssist Spatial Gene Expression Slide (v2, 6.5 mm)
# V5: Visium CytAssist Spatial Gene Expression Slide (v2, 11 mm)
# H1: Visium HD and Visium HD 3' Slides (6.5 mm)

# Capture areas are active regions for capturing expression data on a Visium slide. Slides either have two or four capture areas. Slide areas are named consecutively from top to bottom:
# A1, B1, C1, D1 for Visium Spatial Gene Expression slides, v1, 6.5 mm (V1)
# A, B for Visium CytAssist Spatial Gene Expression slide, v2, 11 mm (V5)
# A1, D1 for Visium CytAssist Spatial Gene Expression slides, v2, 6.5 mm, and Visium HD slides, 6.5 mm (V4 and H1)
		}
	}

	
	# ATAC
	my $type_atac = 1 if map {uc($_) =~ /ATAC/ } @groups;
	system("rm $dir/jobs_atac.txt") if (-f "$dir/jobs_atac.txt");
	if($type_atac ){
		warn 'ATAC';
		open (JOBS_ATAC, ">$dir/jobs_atac.txt");
		my @atac = grep { uc($_->somatic_group()) eq "ATAC"} @{$patients};
		warn "ATAC: ".join(',',map($_->name,@atac));
		foreach my $a (@atac){
			my $vname = $a->name(); 
			my $vfam = $a->family();
			my $vgroup = uc($a->somatic_group());
			my $fastq = $a->getSequencesDirectory;
#			$exec = "cellranger-atac" if $type eq "atac";
			my $transcriptome = $index.'_atac';
			$transcriptome = qx/realpath $index"_atac"/ if (-l $index);
			chomp $transcriptome;
			$transcriptome = $htranscriptome->{atac} if ($choose_transcriptome);
			$transcriptome = '/data-isilon/public-data/10X/HG38/refdata-cellranger-arc-GRCh38-2020-A-2.0.0' if ($project->description =~ /at?traction/i);
			die("No transcriptome") unless ($transcriptome);
			my $cmd = "cd $tmp && $exec count --sample=$vname --id=$vname --fastqs=$fastq --reference=$transcriptome ";
			$cmd = full_cmd($a, $cmd);
			my @cmd_cellranger = grep {/^$exec/} split (/ +&& +/, $cmd);
			warn $cmd_cellranger[0];
			print JOBS_ATAC $cmd;
		}
	}

	
	# CMO
	my $type_cmo = 1 if map {uc($_) =~ /CMO/ } @groups;
	system("rm $dir/jobs_cmo.txt") if (-f "$dir/jobs_cmo.txt");
	if ($type =~ /cmo/ ){
		open (JOBS_CMO, ">$dir/jobs_cmo.txt");
#		warn $patient->somatic_group();
		my @cmo = grep { uc($_->somatic_group()) eq "CMO"} @{$patients};
		warn "CMO: ".join(',',map($_->name,@cmo));
		foreach my $e(@cmo){
			my $ename = $e->name(); 
			my $efam = $e->family();
			my $egroup = uc($e->somatic_group());
			my $transcriptome = $index;
			$transcriptome = qx/realpath $index/ if (-l $index);
			chomp $transcriptome;
			$transcriptome = $htranscriptome->{gex} if ($choose_transcriptome);
			my $lib_file = $dir."/".$ename."_multi.csv";
			open(LIB,">$lib_file") or die "impossible $lib_file";
			print LIB "[gene-expression]\nreference,".$transcriptome."\n$cmo_ref.csv\n";
			
			print LIB "[libraries]\nfastq_id,fastqs,feature_types\n";
			print LIB "$ename,".$tmp.'fastq/'.$ename.",Gene Expression\n";
			my @cmo = grep {$_->family() eq $efam && uc($_->somatic_group()) eq "CMO"} @{$patients};
			my $cmo_name = $cmo[0]->name();
			print LIB "$cmo_name,".$tmp.'fastq/'.$ename.",Antibody Capture\n";
			print LIB "[samples]\nsample_id,cmo_ids\n";
			print LIB $ename."_B251,B251\n";
			print LIB $ename."_B252,B252\n";
			print LIB $ename."_B253,B253\n";
			close(LIB);
			my $cmd = "cd $tmp && $exec multi --id=$ename --csv=$lib_file ";
			$cmd = full_cmd($e, $cmd);
			my @cmd_cellranger = grep {/^$exec/} split (/ +&& +/, $cmd);
			warn $cmd_cellranger[0];
			print JOBS_CMO $cmd;
		}
	}

	
	# ARC
	my $type_arc = 1 if map {uc($_) =~ /ARC/ } @groups;
	system("rm $dir/jobs_arc.txt") if (-f "$dir/jobs_arc.txt");
	if ($type_arc ){
		open (JOBS_ARC, ">$dir/jobs_arc.txt");
		$exec = "cellranger-arc";
		my @arc = grep { uc($_->somatic_group()) eq "ARC"} @{$patients};
		warn "ARC: ".join(',',map($_->name,@arc));
		foreach my $e (@arc){
			my $ename = $e->name(); 
			my $efam = $e->family();
			my $egroup = uc($e->somatic_group());
			my $lib_file = $dir."/".$ename."_library.csv";
			open(LIB,">$lib_file") or die "impossible $lib_file";
			my @atac = grep {$_->family() eq $efam && $_->name() =~ $ename && uc($_->somatic_group()) eq "ARC"} @{$patients};
			next() if scalar(@atac == 0);
			my $atac_name = $atac[0]->name() ;
			my $transcriptome = $index.'_arc';
			$transcriptome = qx/realpath $index"_arc"/ if (-l $index);
			chomp $transcriptome;
			$transcriptome = $htranscriptome->{atac} if ($choose_transcriptome);
			my $lib = "fastqs,sample,library_type\n".$tmp.'fastq/'.",".$ename.",Gene Expression\n";
			$lib .= $tmp.'fastq/'.",".$atac_name.",Chromatin Accessibility\n";
			print LIB $lib;
			close(LIB);
			my $cmd = "cd $tmp && $exec count --id=$ename --transcriptome=$transcriptome  --libraries=$lib_file ";
			$cmd = full_cmd($e, $cmd);
			my @cmd_cellranger = grep {/^$exec/} split (/ +&& +/, $cmd);
			warn $cmd_cellranger[0];
			print JOBS_ARC $cmd;
		}
	}
		
	close(JOBS_EXP);
	close(JOBS_ADT);
	close(JOBS_VDJ);
	close(JOBS_SPATIAL);
	close(JOBS_ATAC);
	close(JOBS_CMO);
	close(JOBS_ARC);


	my $cmd2 = "cat $dir/jobs*.txt | run_cluster.pl -cpu=$cpu";
	warn $cmd2;
	sleep(5) unless ($no_exec);
	my $exit = system ($cmd2) unless ($no_exec);
	exit($exit) if ($exit);
	
	# Open web summaries
	my @error;
	my $web_summaries = "";
	foreach my $patient (@{$patients}) {
		next if ($patient->somatic_group =~ /^ADT$/i);
		my $file = $dir.$patient->name."/web_summary.html";
		$web_summaries .= $file.' ' if (-e $file);
		push(@error, $file) unless (-e $file or $no_exec);
	}
	my $cmd3 = "firefox ".$web_summaries;
	$cmd3 = "google-chrome ".$web_summaries if (getpwuid($<) eq 'shanein');
	warn $cmd3 if ($web_summaries);
	system($cmd3.' &') if ($web_summaries and not $no_exec);
	die("Web summaries not found: ".join(',', @error)) if (@error and not $no_exec);

	unless ($no_exec) {
		print "\t------------------------------------------\n";
		print("\tCheck the web summaries:\n");
		print("\t$dir*/web_summary.html\n");
		print "\t------------------------------------------\n\n";
	}
	
}




#------------------------------
# AGGREGATION
#------------------------------
if (grep(/aggr/, @steps)){
	warn ("Can't aggregate gene expression and vdj librairies together") if (grep {'exp'} @groups and grep {'vdj'} @groups);
#	die("No 'exp' sample to aggregate") unless (grep {'exp'} @groups);
	my $id = $aggr_name if $aggr_name;
	my $aggr_file = $dir."/jobs_aggr.txt";
	open (JOBS_AGGR, ">$aggr_file");
	
	if (grep {'exp'} @groups) {
		$id = 'aggregation_exp' unless ($id);
		my $aggr_csv = $dir."/aggr.csv";
		open (AGGR_CSV, ">$aggr_csv");
		print AGGR_CSV "sample_id,molecule_h5\n";
		foreach my $patient (@{$patients}) {
			warn 
			my $group_type = lc($patient->somatic_group());
			if ($group_type eq "exp") {	# ne "adt" or $_ ne "vdj"
				print AGGR_CSV $patient->name().",".$dir."/".$patient->name()."/molecule_info.h5\n";
			}
		}
		close AGGR_CSV;
		print JOBS_AGGR "cd $dir ; $exec aggr --id=$id --csv=$aggr_csv";
		warn ("cat $aggr_file | run_cluster.pl -cpu=$cpu");
		system ("cat $aggr_file | run_cluster.pl -cpu=$cpu") unless $no_exec;
		die("Error while running cellranger aggr") unless (-d $id);
		if (-d "$dir/$id/" and not $no_exec) {
			print("--------------------\n");
			print("$dir/$id/\n\n");
		}
	}
	
	if (grep {'vdj'} @groups) {
		$id = 'VDJ_'.$id if ($id);
		$id = 'aggregation_vdj' unless ($id);
		my $aggr_csv_vdj = $dir."/aggr_vdj.csv";
		if (-e $aggr_csv_vdj) {
			my $overwrite = prompt("'aggr_vdj.csv' already exists. Overwrite the file ?  ",-yes);
			print ("If 'aggr_vdj.csv' is completed, you can run 'echo \"cd $dir ; $exec aggr --id=$id\_VDJ --csv=$aggr_csv_vdj\" | run_cluster.pl -cpu=$cpu'") unless ($overwrite);
			print JOBS_AGGR "cd $dir ; $exec aggr --id=$id --csv=$aggr_csv_vdj" unless ($overwrite);
			die("'$aggr_csv_vdj' already exists") unless ($overwrite);
		}
		open (AGGR_CSV_VDJ, ">$aggr_csv_vdj");
		print AGGR_CSV_VDJ "sample_id,vdj_contig_info,donor,origin\n" ;
		foreach my $patient (@{$patients}) {
			my $group_type = lc($patient->somatic_group());
			print AGGR_CSV_VDJ $patient->name().",$dir".$patient->name()."/vdj_contig_info.pb,,\n" if ($group_type eq "vdj");
		}
		close (AGGR_CSV_VDJ);
		print("--------------------\n");
		print("Fill the 'donor' and 'origin' columns for each vdj sample in '$aggr_csv_vdj'.\n");
		print("Then run 'echo \"cd $dir && $exec aggr --id=$id --csv=$aggr_csv_vdj\" | run_cluster.pl -cpu=$cpu'\n");
		print("Then make an archive: 'tar -czf $dir/$projectName\_aggr.tar.gz $dir/aggregation_*/web_summary.html $dir/aggregation_*/count/cloupe.cloupe $dir/aggregation_*/count/*_bc_matrix/* $dir/aggregation_*/*/vloupe.vloupe'\n\n");
	}
	close (JOBS_AGGR);
#	my $cmd_tar = "tar -czf $dir/$projectName\_aggr.tar.gz $dir/aggregation_*/web_summary.html $dir/aggregation_*/count/cloupe.cloupe $dir/aggregation_*/count/*_bc_matrix/* $dir/aggregation_*/*/vloupe.vloupe\n";
#	print("Then, you can make an archive: '$cmd_tar'\n");
#	system($cmd_tar) if (-d $dir.'aggregation_exp/' and not $no_exec);
#	print("Archive of aggragation : $dir/$projectName\_aggr.tar.gz\n") if (-d "$dir/$projectName\_aggr.tar.gz" and not $no_exec);

}




#------------------------------
# AGGREGATION VDJ
#------------------------------
if (grep(/aggr_vdj/, @steps)) {
#	my $type = $project->getRun->infosRun->{method};
#	die ("No vdj in project $projectName") if ($type !~ /vdj/);
	die("No 'vdj' sample to aggregate") unless (grep {'vdj'} @groups);
	my $id = $projectName.'_VDJ_aggregation';
	$id = $aggr_name if $aggr_name;
#	my $aggr_file_vdj = $dir."/jobs_aggr_vdj.txt";
#	open (JOBS_AGGR_VDJ, ">$aggr_file_vdj");
	my $aggr_csv_vdj = $dir."/aggr_vdj.csv";
	if (-e $aggr_csv_vdj) {
		my $overwrite = prompt("'aggr_vdj.csv' already exists. Overwrite the file ?  ",-yes);
		print ("If 'aggr_vdj.csv' is completed, you can run 'echo \"cd $dir ; $exec aggr --id=$id --csv=$aggr_csv_vdj\" | run_cluster.pl -cpu=$cpu'") unless ($overwrite);
		die("'$aggr_csv_vdj' already exists") unless ($overwrite);
	}
	open (AGGR_CSV_VDJ, ">$aggr_csv_vdj");
	print AGGR_CSV_VDJ "sample_id,vdj_contig_info,donor,origin\n" ;
	foreach my $patient (@{$patients}) {
		my $group_type = lc($patient->somatic_group());
		if ($group_type eq "vdj") {
			print AGGR_CSV_VDJ $patient->name().",".$dir."/".$patient->name()."/vdj_contig_info.pb,,\n";
		}
	}
	close (AGGR_CSV_VDJ);
	print "--------------------\n";
	print "Fill the 'donor' and 'origin' columns for each vdj sample in '$dir$aggr_csv_vdj'.\n";
	print "Then run 'echo \"cd $dir ; $exec aggr --id=$id\_VDJ --csv=$aggr_csv_vdj\" | run_cluster.pl -cpu=$cpu'";
	print "Then make an archive: 'tar -czf $dir/$projectName\_aggr.tar.gz $dir/aggregation_*/web_summary.html $dir/aggregation_*/count/cloupe.cloupe $dir/aggregation_*/count/*_bc_matrix/* $dir/aggregation_*/*/vloupe.vloupe'";
}



#------------------------------
# CSV with infos/paths for Francesco
#------------------------------
if (grep(/info/i, @steps)){
	system("$Bin/cellranger_infos.pl -project=$projectName -out_dir=$dir");
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
	my $cmd_tar = "$Bin/cellranger_tar.pl -project=$projectName ";
	$cmd_tar .= "-patients=$patients_name " if ($patients_name);
	$cmd_tar .= "-create_bam " if ($create_bam and $create_bam ne 'false');
	$cmd_tar .= "-no_exec " if ($no_exec);
	system($cmd_tar);
}




#------------------------------
# Velocyto
#------------------------------
if (grep(/velocyto/, @steps)) {
	my $cmd_velocyto = "$Bin/velocyto.pl -project=$projectName ";
	$cmd_velocyto .= "-patients=$patients_name " if ($patients_name);
	$cmd_velocyto .= '-transcriptome='.$htranscriptome->{gex}.' '  if ($choose_transcriptome);
	$cmd_velocyto .= "-cpu=$cpu ";
	$cmd_velocyto .= "-no_exec " if ($no_exec);
	system($cmd_velocyto);
}



sub usage {
	print "
$0
-------------
Obligatoires:
	project <s>                nom du projet
	feature_ref	<s>            tableau des ADT, obligatoire seulement si step=count et qu'il y a des ADT
	cmo_ref	<s>                tableau des CMO, obligatoire seulement si step=count et qu'il y a des CMO
Optionels:
	steps <s>                  étape(s) à réaliser: demultiplex, teleport, count, tar, aggr, aggr_vdj, cp, 
	                           cp_web_summary, velocyto, infos ou all (= demultiplex, count, cp_web_summary, tar)
	patients <s>               noms de patients/échantillons, séparés par des virgules
	mismatches <i>             nombre de mismatches autorisés lors du démultiplexage, défaut: 0
	create-bam/nocreate-bam    générer ou non les bams lors du count, défaut: nocreate-bam
	aggr_name <s>              nom de l'aggrégation, lors de step=aggr ou aggr_vdj
	choose_version|version     choisir la version de cellranger/spaceranger à executer
	choose_transcriptome       choisir le transcriptome à utiliser pour les comptages
	chemistry                  chemistry , défaut: auto (pour librairies exp et adt)
	slide_infos <s>            pour spaceranger, ficher csv avec header, contenant 3 colonnes séparées par des virgules: sample, area and slide id
	probe_set <s>              pour spaceranger, probe set à utiliser pour les comptages
	add_image                  pour spaceranger, ajouter une seconde image
	loupe_alignment            pour spaceranger, fichier json d'alignement loupe
	no_exec                    ne pas exécuter les commandes
	force                      relance le pipeline même s'il a déjà tourné
	help                       affiche ce message
	cpu <i>                    nombre de cpu à utiliser, défaut: 20

";
	exit(1);
}

