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
use File::Temp qw/ tempfile tempdir /;
use Term::Menus;
use Proc::Simple;
use Storable;
use JSON::XS;
use XML::Simple qw(:strict);
use Cwd 'abs_path';
use File::Path qw(make_path);
use Text::CSV qw(csv);
use Carp;


my $projectNames;
my $patients_name;
#my $lane;
my $mismatch = 0;
my $multi;
my $no_exec;
my $force;
my $help;

GetOptions(
	'projects=s'				=> \$projectNames,
#	'lane|nb_lane=i'			=> \$lane,
	'mismatches=i'				=> \$mismatch,
	'multi'						=> \$multi,
	'no_exec'					=> \$no_exec,
	'force'						=> \$force,
	'help'						=> \$help,
) || die("Error in command line arguments\n");

usage() if $help;
die("-project argument is mandatory") unless ($projectNames);
die("-mismatches can be 0, 1, 2.") unless ($mismatch =~ /^[012]$/);
warn "mismatch(es)=$mismatch";
warn '-multi=' if ($multi);

my $run;
my $projects;
my $buffer = GBuffer->new();
# Récupère le run et vérifie qu'il n'y en ait qu'un
foreach my $projectName (split(",",$projectNames)){
	my $project = $buffer->newProject( -name => $projectName );
	push(@$projects, $project);
	my $runs = $project->getRuns();
#	warn scalar @{$run->getProjects};
	die("More than one run") if (scalar @$runs > 1 or ($run and $runs->[0]->id ne $run->id));
	$run = $runs->[0];
}
warn 'Run: '. $run->plateform_run_name;
my $machine = $run->infosRun->{machine};
my $bcl_dir = $run->bcl_dir;
warn 'BCL dir: '.$bcl_dir;

# Récupère tous les projets du run
my $dbh = $buffer->dbh();
my $sel = $dbh->selectall_hashref('SELECT project_id FROM PolyprojectNGS.patient WHERE run_id = '.$run->id.' ;', 'project_id');
my @project_ids = keys %$sel;
my $cmd_dbh = 'SELECT name FROM PolyprojectNGS.projects WHERE project_id = '.join(' or project_id = ', @project_ids).' ;';
my $project_names = $dbh->selectall_arrayref($cmd_dbh);
for my $project_name (@$project_names) {
	push(@$projects, $buffer->newProject( -name => $project_name->[0])) unless (grep{$_->name eq $project_name->[0]} @$projects);
}
warn 'Projects in this run:'."\n\t". join("\n\t", map {$_->name} @$projects);
# Récupère tous les patients du run
my @sel = $dbh->selectall_array('SELECT name FROM PolyprojectNGS.patient WHERE run_id = '.$run->id.' ;');
#warn 'Patients in this run: '."\n\t".. join("\n\t",  map {$_->[0]} @sel);
my $patient_names = join(',',map {$_->[0]} @sel);


	
my $config = XMLin("$bcl_dir/RunInfo.xml", KeyAttr => { reads => 'Reads' }, ForceArray => [ 'reads', 'read' ]);
my $lane_count = $config->{Run}->{FlowcellLayout}->{LaneCount};
	
my $samplesheet = $bcl_dir."SampleSheet10X.csv";
my $outcsv_headers;
#push(@$outcsv_headers, ["[Header]"],["FileFormatVersion","2"],["[BCLConvert_Settings]"],["CreateFastqForIndexReads","0"],["TrimUMI","0"],["[BCLConvert_Data]"]);
push(@$outcsv_headers, ["[Settings]"],["CreateFastqForIndexReads","0"],["[Data]"]) unless ($run->sequencing_method eq 'atac');
push(@$outcsv_headers, ["[Settings]"],["CreateFastqForIndexReads","1"],["TrimUMI","0"],["[Data]"]) if ($run->sequencing_method eq 'atac');


my $outcsv;
my $nb_index;
foreach my $project (@$projects){
	my $desc = $project->description =~ s/ /_/r;
#	$desc = $1 if ($desc =~ /^(SC\d+)/);
	my $patients = $project->get_only_list_patients($patient_names);
	@$patients = sort {$a->name cmp $b->name} @$patients;
	
	# pour multiplexage
	my $hpool;
	my $patpool = [];
	foreach my $patient (@$patients) {
		$multi = 1 if ($patient->barcode2 =~ /^BC00\d$|^OB[1-4]$|^CMO\d{3}$|^TotalSeq/ or $patient->somatic_group =~ /^pool/i);
		push(@$patpool,$patient) unless (exists $hpool->{$patient->barcode});
		$hpool->{$patient->barcode} ++;
	}
	warn '-multi' if ($multi);
	$patients = $patpool if ($multi);
	
	foreach my $pat (@$patients) {
		my $tproj = $buffer->getQuery->getProjectDestination($pat->id);
		$desc = $project->description =~ s/ /_/r if ($tproj);
#		$desc = $1 if ($desc =~ /^(SC\d+)/);
		my $pname = $pat->name;
		$pname = $pat->somatic_group if ($multi);
		my $bc = $pat->barcode;
		die ("ERROR in sample barcode for '$pname': $bc.\nShould be SI-XX-[A-H][1-12] or barcode sequence of 8 or 10 nt. without 'N'.")
			unless ($bc =~ /^SI-[GNT][ANST]-[A-H](?:[1-9]|1[0-2])$/ or $bc =~ /^[ATCG]{8,10}$/);
		
		# BC 10X SI-XX
		if ($bc =~ /^SI-([3GNPT][0ANST]3?)-[A-H](?:[1-9]|1[0-2])$/) {
			my $bc_name = $bc;
#			warn $bc_name;
			my $kit = $1;
			my $file_index = "/data-isilon/public-data/10X/sample_indexes_set_sequences/";
			
			# Single Index
			if ($kit =~ /^[NG]A$/) {
				die("Single and dual indexes mixed, run separately") if ($nb_index and $nb_index != 1);
				$nb_index = 1;
				$file_index .= "Single_Index_Kit_N_Set_A.csv" if ($kit eq 'NA');
				$file_index .= "Single_Index_Kit_T_Set_A.csv" if ($kit eq 'GA');
				my $indexes = csv (in => $file_index, sep => ",");
				while($indexes->[0]->[0] =~ /^#|^index_name$/) {shift @$indexes};
				my @barcodes = grep {$_->[0] eq $bc_name} @{$indexes};
				die("ERROR ".scalar @barcodes."barcodes '$bc_name' found in $file_index") unless (scalar @barcodes == 1);
				shift @{$barcodes[0]};
				foreach my $b (@{$barcodes[0]}) {
					push(@$outcsv, [$_,"$pname","$pname","$bc_name","$b","$desc"]) for (1 .. $lane_count);
				}
			}
			
			# Dual Index
			elsif ($kit =~ /^[3GNPT][0ANST]3?$/) {
				die("Single and dual indexes mixed, run separately") if ($nb_index and $nb_index != 2);
				$nb_index = 2;
				$file_index .= "Dual_Index_Kit_$kit\_Set_A.csv";
				my $indexes = csv (in => $file_index, sep => ",");
				while($indexes->[0]->[0] =~ /^#|^index_name$/) {shift @$indexes};
				my @barcodes = grep {$_->[0] eq $bc_name} @{$indexes};
				die("ERROR ".scalar @barcodes."barcodes '$bc_name' found in $file_index") unless (scalar @barcodes == 1);
				my $bc1 = $barcodes[0][1];
				my $bc2 = $barcodes[0][2];
				my $bc2_rc = $barcodes[0][3];
				for my $lane (1 .. $lane_count) {
					 if ($machine eq 'NOVASEQX') {
					 	push(@$outcsv, ["$lane","$pname","$pname","$bc_name","$bc1","$bc_name","$bc2","$desc"]);
						push(@$outcsv, ["$lane","$pname\_RC","$pname\_RC","$bc_name","$bc1","$bc_name","$bc2_rc","$desc"]);
					 }
					 else {
					 	push(@$outcsv, ["$lane","$pname","$pname","$bc_name","$bc1","$bc_name","$bc2_rc","$desc"]);
						push(@$outcsv, ["$lane","$pname\_RC","$pname\_RC","$bc_name","$bc1","$bc_name","$bc2","$desc"]);
					 }
				}
			}
			else {die("Error in $pname BC: $bc_name. No sample index kit set $kit: '$file_index'.")}
		}
		
		# BC sequence
		elsif ($bc =~ /^[ATCG]{8,10}$/ and not $multi) {
			my $bc2 = $pat->barcode2;
			die ("ERROR sample barcode2 for '$pname': $bc2.\nShould be barcode sequence of 8 or 10 nt. without 'N'.") if ($bc2 !~ /^[ATCG]{8,10}$/);
			my $bc2_rc = reverse($bc2) =~ tr/ATCG/TAGC/r;
			for my $lane (1 .. $lane_count) {
				push(@$outcsv, ["$lane","$pname","$pname","$bc","$bc2","$desc"]);
				push(@$outcsv, ["$lane","$pname\_RC","$pname\_RC","$bc","$bc2_rc","$desc"]);
			}
		}
		elsif ($bc =~ /^[ATCG]{8,10}$/ and not $multi) {
			die("ERROR for multiplexed samples, BC should be 10X barcode SI_XX and BC2 OCM/hastag/CMO barcode id.");
		}
		else {die("ERROR BC should be 10X barcode SI_XX or corresponding sequence barcode.")}
	}
}

push(@$outcsv_headers, ["Lane","Sample_ID","Sample_Name","I7_Index_ID","index","Sample_Project"]) if ($nb_index == 1);
push(@$outcsv_headers, ["Lane","Sample_ID","Sample_Name","I7_Index_ID","index","I5_Index_ID","index2","Sample_Project"]) if ($nb_index == 2);
unshift(@$outcsv, @$outcsv_headers);

csv (in => $outcsv, out => $samplesheet, sep_char=> ",");
warn 'Sample Sheet: '.$samplesheet;

my $cmd = "$Bin/../../../dragen/scripts/dragen_demultiplex.pl --project=$projectNames -mismatch=$mismatch";
$cmd .= " -keep_umi -fastq_index" if ($run->sequencing_method eq 'atac');
warn $cmd;
if ($run->sequencing_method eq 'atac') {
	print("For scATACseq, use mask ");
	print colored("Y50;I8;U16;Y50", 'bold'), "\n";
#	print colored("Y50;I8;Y16;Y50", 'bold'), "\n";
}
my $exit = system ($cmd) unless ($no_exec);
confess("ERROR $cmd") if ($exit);




sub usage {
	print "
$0
-------------
Obligatoires:
	project <s>                nom du projet
Optionels:
	patients <s>               noms de patients/échantillons, séparés par des virgules
	mismatches <i>             nombre de mismatches autorisés lors du démultiplexage, peut être 0, 1 ou 2, défaut: 0
	no_exec                    ne pas lancer les commandes dragen_demultiplex après avoir généré la sample sheet
	force                      relance le pipeline même s'il a déjà tourné
	help                       affiche ce message

";
	exit(1);
}

