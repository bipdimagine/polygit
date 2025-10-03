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
use Time::Local 'timelocal';
use Cwd 'abs_path';

use File::Temp qw/ tempfile tempdir /;

use Term::Menus;
use Proc::Simple;
use Storable;
use JSON::XS;
use XML::Simple qw(:strict);


# Toutes ces étapes ne sont plus nécessaires avec cellranger v9
# https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-3p-multi#hashing
  
my $bin_cecile=qq{$Bin/scripts/scripts_db_polypipeline};
my $bin_script_pipeline = qq{$Bin/scripts/scripts_pipeline};


my $projectName;
my $associated_projectName;
my $patients_name;
my @steps;
my $feature_ref;
my $no_exec;
my $aggr_name;
#my $cmo_ref;
my $lane;
my $create_bam = 1;
my $cpu = 20;
my $help;

GetOptions(
	'project=s'				=> \$projectName,
	'associated_project=s'	=> \$associated_projectName,
	'patients=s'			=> \$patients_name,
	'steps=s{1,}'			=> \@steps,
	'create_bam!'			=> \$create_bam,
	'feature_ref=s'			=> \$feature_ref,
	'aggr_name=s'			=> \$aggr_name,
#	'cmo_ref|cmo_set=s'		=> \$cmo_ref,
	'cpu=i'					=> \$cpu,
	'no_exec'				=> \$no_exec,
	'help'					=> \$help,
) || die("Error in command line arguments\n");

#usage() if $help;
die ('-project option required') unless $projectName;
#die ("No -step option specified, can be set to 'demultiplex', 'count', 'aggr', 'cp', 'tar' or 'all'.") unless $step;
#usage() unless ($projectName and $step);
$feature_ref = abs_path($feature_ref) if ($feature_ref);


my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
$patients_name = $project->get_only_list_patients($patients_name);
die("No patient in project '$projectName'") unless $patients_name;
unless ($project->isSomatic) {
	my $warn = "/!\\ Project $projectName is not in somatic mode. Activate somatic mode and check that the groups have been filled in, "
		."so that they can be taken into account in the analysis.";
	die ($warn) unless (grep{/@steps/} ('demultiplex', 'teleport', 'tar', 'cp', 'cp_web_summaries'));
	warn ($warn);
}

@steps = split(/,/, join(',',@steps));
unless (@steps) {
	my $liste_steps = ['demultiplex', 'count_multiplex', 'bamtofastq', 'count', 'aggr', 'tar', 'cp', 'all'];
	my %Menu_1 = (
		Item_1 => {
			Text   => "]Convey[",
			Convey => $liste_steps,
		},
		Select => 'Many',
		Banner => "   Select steps:"
	);
	@steps = &Menu( \%Menu_1 );
	die if ( @steps eq ']quit[' );
}
warn 'steps='.join(',',@steps);

my $exec = "cellranger multi";
warn $exec;
my $run = $project->getRun();
my $run_name = $run->plateform_run_name;
my $type = $run->infosRun->{method};
#my $dir = $project->getProjectRootPath();
my $dir = $project->getCountingDir('cellranger');
warn $dir;
#my $tmp = $project->getAlignmentPipelineDir("cellranger");
#warn $tmp;


my $hfamily;
my $prog;
foreach my $patient (@{$patients_name}) {
	my $hpatient;
	$hpatient->{name}= $patient->name();
#	$hpatient->{family}= $patient->family();
	$hpatient->{group}= $patient->somatic_group();
	$hpatient->{bc} = $patient->barcode2;
	$hpatient->{objet} = $patient;
	my $sbc = $patient->barcode; 
	push(@{$hfamily->{$sbc}},$hpatient);
	$prog = $patient->alignmentMethod();
}

my @group;
foreach my $patient (@{$patients_name}) {
	my $group = $patient->somatic_group();
	push(@group,$group);
}



###############
# DEMULIPLEXAGE
###############
if (grep(/^demultiplex(ing|age)?$|^all$/i, @steps)){

	my $tmp = $project->getAlignmentPipelineDir("demultiplex");
	
	my $bcl_dir = $run->bcl_dir;
	warn $bcl_dir;
	unless ( $lane ){
		my $config = XMLin("$bcl_dir/RunInfo.xml", KeyAttr => { reads => 'Reads' }, ForceArray => [ 'reads', 'read' ]);
		$lane = $config->{Run}->{FlowcellLayout}->{LaneCount};
		warn "LaneCount=$lane";
	}
	
	my $sampleSheet = $bcl_dir."/sampleSheet.csv";
	open (SAMPLESHEET,">$sampleSheet");
	print SAMPLESHEET "Lane,Sample,Index\n";
	for (my $i=1;  $i<=$lane;$i++){
		foreach my $k (keys(%$hfamily)){
			my @pat = @{$hfamily->{$k}};
			print SAMPLESHEET $i.",".$pat[0]->{group}.",".$k."\n";
		}
	}
	close(SAMPLESHEET);
	
	my $cmd = "cd $tmp; /software/bin/demultiplex.pl -dir=$bcl_dir -run=$run_name -hiseq=10X -sample_sheet=$sampleSheet -cellranger_type=$exec";
	warn $cmd;
	my $exit = system ($cmd) unless $no_exec==1;
	exit($exit) if ($exit);
}




###############
# cellranger multi (on multiplexed samples)
###############
if (grep(/^demultiplex_hastag|^count_multiplex|^all$/i, @steps)){
	
	# Sur PolyProject, créer un nouveau run et projet avec les échantillons démultiplexés/dé-hashtagés,
	# dont le group correspond à l'échantillon multiplexé,
	# la famille sans VDJ_ ou ADT_ au début,
	# le BC2 au nom du CMO/hashtag associé
	
	die("-feature_ref option is mandatory with step 'demultiplex_hastags'") unless ($feature_ref);
	die("-associated_project option is mandatory with step 'demultiplex_hastags'") unless ($associated_projectName);
	my $associated_project = $buffer->newProject(-name=>$associated_projectName);
	my @associated_patients = @{$associated_project->getPatients};
	my @associated_pnames = map {$_->name} @associated_patients;
	
	my $adt_set;
	die("'$feature_ref' does not exists") unless (-e $feature_ref);
	open(ADT, '<', $feature_ref) or die("Can't open '$feature_ref'");
	while (my $line = <ADT>) {
		chomp $line;
		my $adt_id = (split /,/ , $line)[0];
		push (@$adt_set, $adt_id) unless ($adt_id eq 'id');
	}
	close(ADT);
	
	## Step 1.2: Config CSV setup for Cell Ranger
	# We will use cellranger multi to demultiplex samples. The multi pipeline accesses library information from the config CSV.
	# The demultiplexing config CSV file must have the following sections:
	# [gene-expression]
	# reference,/path/to/transcriptome
	# create-bam,true/false
	# [vdj]
	# reference,/path/to/vdj_reference
	# [feature]
	# reference,/path/to/feature_ref.csv
	# [libraries]
	# fastq_id,fastqs,feature_types
	# gex1,/path/to/gex_fastqs,Gene Expression
	# vdj,/path/to/vdj_fastqs,VDJ
	# ab1,/path/to/ab_fastqs,Antibody Capture
	# [samples]
	# sample_id,hashtag_ids
	# Sample1,TotalSeqC_Hashtag_1
	# Sample2,TotalSeqC_Hashtag_2
	# Sample3,TotalSeqC_Hashtag_3
	# Sample4,TotalSeqC_Hashtag_4

	system("mkdir $dir/multiplex") unless (-e $dir."multiplex");
	my $tmp = $project->getAlignmentPipelineDir("cellranger_multi");
	
	sub full_cmd {
		my ($pool, $cmd, $dir, $tmp) = @_;
		chomp $cmd;
		my $pool_name = $pool->name;
		my $tmp_fastq = $tmp.'fastq/';
		system("mkdir $tmp_fastq") unless (-d "$tmp_fastq");
		my $seq_dir = $pool->getSequencesDirectory;
		$cmd =~ s/$seq_dir/$tmp_fastq/;
		system("mkdir $dir/multiplex/$pool_name/") unless (-d "$dir/multiplex/$pool_name/");
#		my $cmd2 .= " && rsync -ra $tmp$pool_name/outs/per_sample_outs/* $tmp$pool_name/_versions $dir/multiplex/$pool_name/ ";
		my $cmd2 .= " && rsync -ra $tmp$pool_name $dir/multiplex/ ";
		system("mkdir $dir/multiplex/$pool_name") unless (-d "$dir/multiplex/$pool_name");
		return $cmd.$cmd2."\n";
	}
	
	my @fastq;
	foreach my $apat (@associated_patients) {
		push(@fastq, $apat->getSequencesDirectory.$apat->name.'*.fastq.gz');
	}
	my $cp = 'rsync -v '.join (' ', @fastq)." $tmp/fastq/";
	warn ($cp) if ($no_exec);
	system ($cp) unless ($no_exec);
	
	my $index = $project->getGenomeIndex($prog);
	open (JOBS, ">$dir/jobs_count_multiplex.txt");
	foreach my $k (keys(%$hfamily)){
		my @pat = @{$hfamily->{$k}};
		my $poolName = $pat[0]->{group};
		next if ($poolName =~ /^VDJ_|^ADT_/);
		die("'$poolName' does not exists in the associated project '$associated_projectName'") unless (grep {$poolName} @associated_pnames);
		next if ($poolName =~ /^(ADT|CMO)_/);
		my $demultiplex_csv = $dir.'multiplex/'.$poolName."_config.csv";
		open (CSV,">$demultiplex_csv") or die ("Can't open '$demultiplex_csv: $!");
		print CSV "[gene-expression]\n";
		print CSV "reference,".$index."\n";
		print CSV "create-bam,true\n" if ($create_bam);
		print CSV "create-bam,false\n" unless ($create_bam);
		print CSV "\n";
		print CSV "[vdj]\n";
		print CSV "reference,".$index."_vdj"."\n";
		print CSV "\n";
		print CSV "[feature]\n";
		print CSV "reference,".$feature_ref."\n";
		print CSV "\n";
		print CSV "[libraries]\n";
		my $pobj=$pat[0]->{objet};
		print CSV "fastq_id,fastqs,feature_types\n";
		my @pool = grep {$_->name =~ /^$poolName$/} @associated_patients;
		die("Found ".scalar @pool." $poolName in the associated project '$associated_projectName'") unless (scalar @pool == 1);
		print CSV $poolName.",$tmp/fastq/,Gene Expression\n";
		my @vdj = grep {$_->name =~ /^VDJ_$poolName$/} @associated_patients;
		die("Found ".scalar @pool." $poolName in the associated project '$associated_projectName'") unless (scalar @pool == 1);
		print CSV $vdj[0]->name.",$tmp/fastq/,VDJ\n";
		my @adt = grep {$_->name =~ /^ADT_$poolName$/} @associated_patients;
		die("No cmo/adt 'CMO_".$poolName."' or 'ADT_".$poolName."' found in associated project '$associated_projectName'") unless (scalar @adt);
		print CSV $adt[0]->name.",$tmp/fastq/,Antibody Capture\n";
		print CSV "\n";
		print CSV "[samples]\n";
		print CSV "sample_id,hashtag_ids\n";
		foreach my $p (@pat) {
			die("Can't find bc2 '".$p->{bc}."' in feature reference csv\n", Dumper $adt_set) unless (grep(/^$p->{bc}$/, @$adt_set));
			print CSV $p->{name}.','.$p->{bc}."\n";
		}
		close (CSV);
		my $cmd = "cd $tmp; cellranger multi --id=$poolName --csv=$demultiplex_csv";
		$cmd = full_cmd($pool[0], $cmd, $dir, $tmp);
		print JOBS $cmd;
	}
	close JOBS;
	
	## Step 1.3: Run the demultiplexing command
	my $exit = system("cat $dir/jobs_count_multiplex.txt | run_cluster.pl -cpu=$cpu") unless ($no_exec);
	exit($exit) if ($exit);
		
	# Open web summaries
	my @error;
	my $web_summaries = "";
	foreach my $patient (@$patients_name) {
		next if ($patient->somatic_group =~ /^ADT$/i);
		my $file = $dir.'multiplex/*/outs/per_sample_outs/'.$patient->name."/web_summary.html";
		$web_summaries .= $file.' ' if (-e $file);
		push(@error, $file) unless (-e $file or $no_exec);
	}
	my $cmd3 = "firefox ".$web_summaries;
	$cmd3 = "google-chrome ".$web_summaries if (getpwuid($<) eq 'shanein');
	warn $cmd3 if ($web_summaries);
	system($cmd3.' &') if ($web_summaries and not $no_exec);
	die("Web summaries not found: ".join(', ', @error)) if (@error and not $no_exec);
}



###############
# BamToFastq
###############
if (grep(/bamtofastq|^all$/i, @steps)) {
	my $tmp = $project->getAlignmentPipelineDir("cellranger_bamtofastq");
	system("mkdir $tmp") unless (-e "$tmp");
	system("mkdir $dir/bamtofastq") unless (-e "$dir/bamtofastq");
	
	## Step 1.5: Find and record the total number of reads sequenced per library
	# $dir/multiplex/$poolName/$sample/metrics_summary.csv 
	# Library,Gene Expression,Fastq ID,Phago_1,Number of reads,"545,239,728"
	my $nb_reads;
	open(JOBS, '>', "$dir/jobs_bamtofastq.txt") or die("Can't open '$dir/jobs_bamtofastq.txt': $!");
	foreach my $k (keys(%$hfamily)) {
		my @pat = @{$hfamily->{$k}};
		my $poolName = $pat[0]->{group};
		next if ($poolName =~ /^(ADT|VDJ)_/);
		my $pname = $pat[0]->{name};
		my $metrics_csv = "$dir/multiplex/$poolName/$pname/metrics_summary.csv";
		open(METRICS, '<', $metrics_csv) or die("Can not open '$metrics_csv': $!");
		while (my $line = <METRICS>) {
#			next unless ($line =~ /^Library,Gene Expression,Fastq ID,$poolName,Number of reads,/);
			chomp $line;
			$line =~ /^Library,(Gene Expression|VDJ [BT]),Fastq ID,(VDJ_)?$poolName,Number of reads,"([\d,])+"$/ || next;
			$nb_reads->{$poolName}->{'GEX'} = $3 if ($1 eq 'Gene Expression');
			$nb_reads->{$poolName}->{'VDJ'} = $3 if ($1 =~ /VDJ [BT]/);
#			$nb_reads->{$poolName} = $3;
		}
		close(METRICS);
#		warn $nb_reads->{$poolName};
	
		## Step 2: Convert per sample BAM files to FASTQs for the GEX data
		# Since cellranger multi requires FASTQ files as the input, you must convert the BAM files back to per individual FASTQ files.
		# Run bamtofastq
		# We recommend setting the `--reads-per-fastq` argument higher than the total number of reads recorded in Step 1.5
		# (so that the output FASTQ files are not split into chunks by bamtofastq).
		my $fastq_dir = $pat[0]->{objet}->getSequencesDirectory;
		
		my $rounded;
		foreach my $feature_type ('GEX', 'VDJ') {
			$rounded->{$feature_type} = int(($nb_reads->{$poolName}->{$feature_type} + 9_999_999) / 10_000_000) * 10_000_000; # Arrondi à la dizaine de million supérieure
			$rounded->{$feature_type} = 500_000_000 unless ($rounded->{$feature_type} >= 500_000_000);
		}
		foreach my $pat (@pat) {
			my $pname = $pat->{name};
			my $tmp_fastq = $tmp."$pname/";
			my $bam = $dir."multiplex/$poolName/$pname/count/sample_alignments.bam";
			die ("No bam found for sample '$pname', pool '$poolName': '$bam") unless (-e $bam);
			my $cmd = "cellranger bamtofastq --reads-per-fastq=".$rounded->{'GEX'}." $bam $tmp$pname --nthreads=$cpu";
			my $bam_vdj = $dir."multiplex/$poolName/$pname/vdj_b/concat_ref.bam" if (-d $dir."multiplex/$poolName/$pname/vdj_b/");
			$bam_vdj = $dir."multiplex/$poolName/$pname/vdj_t/concat_ref.bam" if (-d $dir."multiplex/$poolName/$pname/vdj_t/");
#			my $cmd_vdj = "cellranger bamtofastq --reads-per-fastq=".$rounded->{'VDJ'}." $bam_vdj $tmp/VDJ_$pname --nthreads=$cpu";
			system("mkdir $tmp/VDJ_$pname") unless (-d "$tmp/VDJ_$pname");
			my $cmd_vdj = "samtools fastq $bam_vdj -2 $tmp/VDJ_$pname/VDJ_$pname\_R2.fastq.gz --threads=$cpu";
			print JOBS $cmd."\n".$cmd_vdj."\n";
		}
#		warn $cmd if ($no_exec);
	}
	close(JOBS);
	my $exit = system("cat $dir/jobs_bamtofastq.txt | run_cluster.pl -cpu=$cpu") unless ($no_exec);
#	exit($exit) if ($exit);
die;
				
	## Step 2.2: Identify the FASTQ directory corresponding to GEX
	# Please note that the order of the libraries depends on the order that the libraries are listed in '$config.csv' (Step 1.3).
	my $cmd2;
	my %errors;
	open(JOBS, '>', "$dir/bamtofastq/jobs_mv_fastq.txt") or die("Can't open '$dir/bamtofastq/jobs_mv_fastq.txt': $!");
	foreach my $k (keys(%$hfamily)) {
		my @pat = @{$hfamily->{$k}};
		my $poolName = $pat[0]->{group};
		next if ($poolName =~ /^(ADT|CMO)_/);
		my $fastq_dir = $pat[0]->{objet}->getSequencesDirectory;
		foreach my $pat (@pat) {
			my $pname = $pat->{name};
			my $tmp_fastq = $tmp."$pname/";
			$errors{$pname} = $tmp_fastq unless ($tmp_fastq or $no_exec);
			next unless (-d $tmp_fastq);
			opendir(DIR, $tmp_fastq) or die "Can't opendir $tmp_fastq: $!";
			my @bamtofastq_dir = grep {/^$poolName\_[01]_1_[A-Z0-9]{9}$/} readdir(DIR);	# $poolName\_0_1_* pour les GEX, $poolName\_1_1_* pour les VDJ (dans l'ordre du csv config)
			close(DIR);
			warn Dumper \@bamtofastq_dir;
			warn ("No directory found in $tmp_fastq") unless (scalar @bamtofastq_dir or not $no_exec);
			die ("No directory found in $tmp_fastq") unless (scalar @bamtofastq_dir or $no_exec);
			die ("Expected 2 directories, found ".scalar @bamtofastq_dir) unless (scalar @bamtofastq_dir == 2 or $no_exec);
			my $gex_dir = $tmp_fastq.$bamtofastq_dir[0].'/';
			opendir(DIR, $gex_dir) or die "Can't opendir $gex_dir: $!";
			while (my $file = readdir(DIR)) {
				next unless ($file =~ /^bamtofastq_S1_L[0-9]{3}_R[1-2]_001.fastq.gz$/);
				my $tmp_cmd2 .= "cp $gex_dir$file $fastq_dir".$file=~s/^bamtofastq/$pname/r;
				$cmd2 .= $tmp_cmd2."\n";
				print JOBS $tmp_cmd2."\n";
			}
			my $vdj_dir = $tmp_fastq.$bamtofastq_dir[1].'/';
			close(DIR);
			opendir(DIR, $vdj_dir) or die "Can't opendir $vdj_dir: $!";
			while (my $file = readdir(DIR)) {
				next unless ($file =~ /^bamtofastq_S1_L[0-9]{3}_R[1-2]_001.fastq.gz$/);
				my $tmp_cmd2 .= "cp $vdj_dir$file $fastq_dir".$file=~s/^bamtofastq/VDJ_$pname/r;
				$cmd2 .= $tmp_cmd2."\n";
				print JOBS $tmp_cmd2."\n";
			}
			close(DIR);
		}
	}
	close(JOBS);
#	warn $cmd2 if ($no_exec);
	chomp $cmd2;
	warn "Copying fastq...\n" unless ($no_exec);
	my $exit = system($cmd2) unless ($no_exec);
	exit($exit) if ($exit);
#	die ("Could not copy fastq for samples '".join("', '", keys %errors)."'. GEX directories not found:\n".join("\n", values %errors)) if (keys %errors);
}



###############
# cellranger multi (on demultiplexed samples)
###############
if (grep(/^count(_demultiplex(ed)?)?$|^all$/i, @steps)){
	system("mkdir $dir/count") unless (-e "$dir/count");
	my $tmp = $project->getAlignmentPipelineDir("cellranger_count");
	
	my $index = $project->getGenomeIndex($prog);
	my $set = $index."/probe_sets/";
	$set .= "Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv" if ($project->getVersion() =~ /^HG/);
	$set .= "Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv" if ($project->getVersion() =~ /^MM/);
	
	## Step 1.4: Find and record the number of cells called per sample
	my $nb_cells;
	my $cmd;
	open (JOBS, ">$dir/jobs_count.txt");
	foreach my $k (keys(%$hfamily)) {
		my @pat = @{$hfamily->{$k}};
		my $poolName = $pat[0]->{group};
		foreach my $pat (@pat) {
			my $pname = $pat->{name};
			next if ($pname =~ /^(ADT|VDJ|CMO)_/);
			my $metrics_csv = "$dir/demultiplexing/$poolName/outs/per_sample_outs/$pname/metrics_summary.csv";
			open(METRICS, '<', $metrics_csv) or die("Can not open '$metrics_csv': $!");
			while (my $line = <METRICS>) {
				next unless ($line =~ /^Cells,Gene Expression,,,Cells,/);
				chomp $line;
				$nb_cells->{$pname} = $line =~ s/^Cells,Gene Expression,,,Cells,//rg;
				$nb_cells->{$pname} =~ s/"|,//g;
			}
			close(METRICS);
			
	## Step 3.3: Create sample-specific multi config CSV files
#			next unless ($pname eq $pat->{objet}->family); 
#			next if ($pname =~ /^VDJ_/);
#			warn $pname;
#			warn $pat->{objet}->family;
#			warn ($pname eq $pat->{objet}->family);
			my $csv_count = $dir."count/".$pname.".csv";
			open (CSV,">$csv_count") or die ("Can't open '$csv_count': $!");
			print CSV "[gene-expression]\n";
			print CSV "reference,".$index."\n";
#			print CSV "probe-set,".$set."\n";
			print CSV "create-bam,true\n" if ($create_bam);
			print CSV "create-bam,false\n" unless ($create_bam);
#			print CSV "force-cells,".$nb_cells->{$pname}."\n";		# donne une erreur
			print CSV "expect-cells,".$nb_cells->{$pname}."\n";
			print CSV "check-library-compatibility,false\n";
			print CSV "\n";
			print CSV "[vdj]\n";
			print CSV "reference,".$index."_vdj\n";
			print CSV "\n";
			# todo: ajouter ADT config csv
			#print CSV "[feature]\n";
			#print CSV "reference,".abs_path($feature_ref);
			#print CSV "\n";
			print CSV "[libraries]\n";
			my $pobj=$pat->{objet};
			my $fastq = $pobj->getSequencesDirectory();
			print CSV "fastq_id,fastqs,feature_types\n";
			print CSV $pname.",".$fastq.",Gene Expression\n";
			# todo: ajouter ADT et VDJ
			my @VDJ = grep {$_->name eq 'VDJ_'.$pat[0]->{name}} @$patients_name;
			print CSV $VDJ[0]->name.",".$VDJ[0]->getSequencesDirectory.",VDJ\n" if (scalar @VDJ);
#			my @ADT = grep {$_->name eq 'ADT_'.$pat[0]->{name}} @$patients_name;
#			print CSV $ADT[0]->name.",".$ADT[0]->getSequencesDirectory.",Antibody Capture\n" if (scalar @ADT);
			print CSV "\n";
#			print CSV "[samples]\n";
#			print CSV "sample_id,probe_barcode_ids\n";
#			foreach my $p (@pat){
#				print CSV $p->{name}.",";
#				print CSV $p->{bc}."\n";
#			}
			my $tmp_cmd = "cd $dir/count; cellranger multi --id=$pname --csv=$csv_count";
			$tmp_cmd = full_cmd($pat, $tmp_cmd, $dir, $tmp);
			$cmd .= $tmp_cmd."\n";
			print JOBS $tmp_cmd."\n";
			close CSV;
		}
	}
	close JOBS;
	
#	warn Dumper $nb_cells;
	
	warn $cmd if ($no_exec);
	chomp $cmd;
#	system("echo \"$cmd\" | run_cluster.pl -cpu=$cpu") unless ($no_exec);
	my $exit = system ("cat $dir/jobs_count.txt | run_cluster.pl -cpu=$cpu") unless ($no_exec);
	exit($exit) if ($exit);
}


# todo: count mais avec cellranger count ou cellranger vdj séparément ?


if (grep(/aggr|^all$/i, @steps)){
	my $aggr_file = $dir."/jobs_aggr.txt";
	my $id = $projectName;
	$id = $aggr_name if $aggr_name;
	open (JOBS_AGGR, ">$aggr_file");
	if ($type =~ /vdj/) {
		print JOBS_AGGR "sample_id,vdj_contig_info,donor,origin\n" ;
	}
	else{
		print JOBS_AGGR "sample_id,molecule_h5\n";
	}	
	foreach my $patient (@{$patients_name}) {
		my $group_type = uc($patient->somatic_group());
		if ($group_type eq "VDJ") {
			print JOBS_AGGR $patient->name().",".$dir."/".$patient->name()."/"."outs/vdj_contig_info.pb\n";
		}
		else{	
			print JOBS_AGGR $patient->name().",".$dir."/".$patient->name()."/"."outs/molecule_info.h5\n";
		}
	}
	close JOBS_AGGR;
	my $aggr_cmd = "cd $dir ; $exec aggr --id=$id --csv=$aggr_file";
	warn $aggr_cmd;
	die();
	system ($aggr_cmd) unless $no_exec==1;
}


#foreach my $patient (@{$patients_name}) {
#	my $name = $patient->name();
#	my $file = $dir."/".$name."/outs/web_summary.html";
#	print $file."\n" if -e $file ;
#}


if (grep(/^tar$|^all$/i, @steps)){
	my $tar_cmd = "tar -cvzf $dir/$projectName.tar.gz $dir/count/*/outs/per_sample_outs/*/web_summary.html \
	$dir/count/*/outs/per_sample_outs/*/count/sample_cloupe.cloupe \
	$dir/*/outs/per_sample_outs/*/vdj*/*vloupe.vloupe \
	$dir/*/outs/per_sample_outs/*/count/sample_filtered_bc_matrix/* ";
#	$tar_cmd .= "$dir/count/*/outs/multi/count/raw_feature_bc_matrix/* ";
#	count/PP_Phago1_1/outs/per_sample_outs/PP_Phago1_1/web_summary.html
#	count/PP_Phago1_1/outs/per_sample_outs/PP_Phago1_1/count/sample_cloupe.cloupe
#	count/PP_Phago1_1/outs/per_sample_outs/PP_Phago1_1/vdj_t/vloupe.vloupe
#	count/PP_Phago1_1/outs/per_sample_outs/PP_Phago1_1/count/sample_filtered_feature_bc_matrix/
#	count/PP_Phago1_1/outs/multi/count/raw_feature_bc_matrix/
	
	die ("archive $dir/$run.tar.gz already exists") if -e $dir."/".$run.".tar.gz";
	system ($tar_cmd) unless $no_exec==1;
	#or die "impossible $tar_cmd";
	print "\t##########################################\n";
	print "\tlink to send to the users : \n";
	print "\twww.polyweb.fr/NGS/$projectName/$run.tar.gz \n";
	print "\t##########################################\n";
}

# Commande pour copier sur /data-isilon/singleCell
if (grep(/^cp$|^all$/i, @steps)){
	my $dirout = "/data-isilon/SingleCell/$projectName";
	
	system("mkdir $dirout") unless -e $dirout;
	my $hgroups ;
	foreach my $patient (@{$patients_name}) {
		$hgroups->{$patient->somatic_group}++ 
	}
	foreach my $group_name (keys %$hgroups){
		warn  "$dir/$group_name";
		my $cmd= qq{rsync -rav  $dir/$group_name $dirout && cp $dir/$group_name.csv $dirout };
		warn $cmd;
		system($cmd);
	}
	print "\t##########################################\n";
	print "\tcp to $dirout \n";
	print "\t##########################################\n";
}

exit;
