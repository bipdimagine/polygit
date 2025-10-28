#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo/lib/";
use lib "$Bin/../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../GenBo/lib/obj-nodb/packages/";
use lib "$Bin/../../GenBo/lib/GenBoDB/writeDB/";
use lib "$Bin/../../polypipeline/packages/";
use Data::Dumper;
use Getopt::Long;
use Carp;
use IO::Prompt;
use Term::Menus;
use colored;
use Cwd 'abs_path';
use Net::SFTP;
use Archive::Tar;
use file_util;

use GBuffer;
my $buffer = new GBuffer;

my $sftp_host = "192.168.2.111";
my $sftp_user = "bioinfo";
my $sftp_key_path = "/home/masson/.ssh/known_hosts";


my $project_names;
my $patient_names;
my $file_types;
my $archive_dir = "/data-pure/workspace/downloads/transfert_file/";
my $archive_name;
my $transfert_sftp;
my $no_die;
my $no_archive;
my $help;

GetOptions(
	'projects=s'		=> \$project_names,    # Projet(s)
	'patients=s'		=> \$patient_names,    # Patient(s)
	'files|type=s'		=> \$file_types,       # fastq, bam, vcf, htlv1 or all
	'archive_dir=s'		=> \$archive_dir,      # directory path
	'archive_name=s'	=> \$archive_name,
	'sftp!'				=> \$transfert_sftp,
	'no_die'			=> \$no_die,
	'no_archive'		=> \$no_archive,
	'help'				=> \$help,
) || die("Error in command line arguments\n");

usage() if ($help);
die("-project is required") unless ($project_names);
$archive_dir .= '/' unless ($archive_dir =~ /\/$/);
confess("Directory \"$archive_dir\" does not exit") unless ( -d $archive_dir );

$transfert_sftp = prompt( "Do you want to tranfert files to the sftp ?  (y/n)  ", -yes_no ) unless ( defined($transfert_sftp) );

my $sftp;
my $dir_sftp;
if ($transfert_sftp) {

	# Connexion au serveur sftp, à faire sous le login de Cécile
	my $current_user = getpwuid($<);
	confess("You must be logged in under Cécile user ID: run 'su masson' before running the command line. Current user: '$current_user'") if ( $current_user ne "masson" );
	$sftp = Net::SFTP->new(
		$sftp_host,
		user     => $sftp_user,
		key_path => $sftp_key_path
	);

	# Liste les fichiers dans le répertoire distant
	my @ls           = $sftp->ls('.');
	my @sorted_files = sort { $a->{filename} cmp $b->{filename} } @ls;
	@sorted_files =
	  grep { $_->{filename} !~ /^\./ } @sorted_files;    # !~ /^\.{1,2}$/
	my @file_names = map { $_->{filename} } @sorted_files;
	$dir_sftp = prompt( "Select a directory to put the files: ", -menu => \@file_names );
	$sftp = undef;
	$dir_sftp .= '/filetransfer/' unless $dir_sftp =~ /filetransfer$/;
	$dir_sftp =~ s/^/\// unless ($dir_sftp =~ /^\//);
	print "\n";
}

my @file_types;
unless ($file_types) {
	my @types  = [ 'fastq', 'bam', 'vcf', 'htlv1' ];
	my %Menu_1 = (
		Item_1 => {
			Text   => "]Convey[",
			Convey => @types,
		},
		Select => 'Many',
		Banner => "   Select file types to transfert:"
	);
	@file_types = &Menu( \%Menu_1 );
	die if ( @file_types eq ']quit[' );
}
else {
	@file_types = split( ',', $file_types );
}

# Récupère les fichiers
my $files;
my @patient_names = split( ',', $patient_names ) if ( $patient_names ne "all" );
$project_names = [ split( ',', $project_names ) ];
foreach my $project_name (@$project_names) {
	my $project = $buffer->newProject( -name => $project_name );
	colored::stabilo( "white", "Project " . $project->name );

	my $patients = $project->get_only_list_patients($patient_names);
	warn( "No patient in project " . $project->name . "\n" ) unless $patients;

	foreach my $pat (@$patients) {
		my $pat_name = $pat->name;
		@patient_names = grep { $_ ne $pat_name } @patient_names;
		colored::stabilo( "yellow", $pat_name );

		#fastq
		if ( grep { $_ =~ /fastq/i or $_ =~ /all/i } @file_types ) {
			my $fastq_files;
			print "fastq:\n";
			eval {
				foreach my $fastq ( @{ file_util::find_file_pe($pat) } ) {	# $pat->fastqFiles
					foreach my $f ( values %$fastq ) {
						next unless ( $f =~ /_R(1|2)(_\d{3,})?.fastq.gz$/ );
						push( @$fastq_files, ( $f ) );
						my $index = $f =~ s/_R(1|2)(_\d{3,})\.fastq\.gz$/_I$1($2)\.fastq\.gz/r;
						push( @$fastq_files, ( $index ) ) if ( -e $index );
					}
				}
#				warn Dumper $fastq_files;
				print "\t" . join( "\n\t", @$fastq_files ) . "\n";
				push( @$files, @$fastq_files );
			};
			if ($@) {
				colored::stabilo( "red", "No fastq file found" );
				die($@) unless ($no_die);
				warn $@;
			}
		}

		# bam
		if ( grep { $_ =~ /bam/i or $_ =~ /all/i } @file_types ) {
			print "bam:\n";
			my $bam_files = $pat->getBamFiles;
			if ( scalar @$bam_files ) {
				print "\t" . join( "\n\t", @$bam_files ) . "\n";
				push( @$files, @$bam_files );

				# *.bam.bai
				foreach my $bam (@$bam_files) {
					$bam .= '.crai' if ($bam =~ /\.cram$/);
					$bam .= '.bai' if ($bam =~ /\.bam$/);
				}
				print "\t" . join( "\n\t", @$bam_files ) . "\n";
				push( @$files, @$bam_files );
			}
			else {
				colored::stabilo( "red", "No bam file found" );
				print "\n";
				die("No bam/cram file found for patient '$pat_name': ".$pat->getBamFileName) unless ($no_die);
			}
		}

		# vcf
		if ( grep { $_ =~ /vcf/i or $_ =~ /all/i } @file_types ) {
			print "vcf:\n";
			my $vcf_files = $pat->getVariationsFiles;
			if ( scalar @$vcf_files ) {
				print "\t" . join( "\n\t", @$vcf_files ) . "\n";
				push( @$files, @$vcf_files );

				# *.vcf.tbi
				foreach my $vcf (@$vcf_files) {
					$vcf .= '.tbi';
				}
				print "\t" . join( "\n\t", @$vcf_files ) . "\n";
				push( @$files, @$vcf_files );
			}
			else {
				colored::stabilo( "red", "No vcf file found" );
				print "\n";
				die("No vcf file found for patient '$pat_name'") unless ($no_die);
			}
		}

		# htlv1
		if ( grep { $_ =~ /htlv1/i } @file_types ) {
			print "htlv1:\n";
			my $htlv1_dir = $project->getVariationsDir("htlv1_calling");
			my $f = File::Util->new();
			my @htlv1_files = $f->list_dir($htlv1_dir => {
				files_only => 1,
				with_paths => 1,
				files_match => qr/^$pat_name-(clonalityResults.txt|mergedIS.txt|mergedIS.xls|SIMPLIFIED_mergedIS.txt)$/
			});
			if ( scalar @htlv1_files ) {
				print "\t" . join( "\n\t", @htlv1_files ) . "\n";
				push( @$files, @htlv1_files );
			}
			if ( scalar @htlv1_files < 4 ) {
				colored::stabilo( "magenta", 4 - scalar(@htlv1_files) . " missing file(s):" );
				unless ( grep { $_ =~ /-clonalityResults.txt$/ } @htlv1_files ) {
					print "\t$htlv1_dir$pat_name-clonalityResults.txt\n";
				}
				unless ( grep { $_ =~ /-mergedIS.txt$/ } @htlv1_files ) {
					print "\t$htlv1_dir$pat_name-mergedIS.txt\n";
				}
				unless ( grep { $_ =~ /-mergedIS.xls$/ } @htlv1_files ) {
					print "\t$htlv1_dir$pat_name-mergedIS.xls\n";
				}
				unless ( grep { $_ =~ /-SIMPLIFIED_mergedIS.txt$/ }	@htlv1_files ) {
					print "\t$htlv1_dir$pat_name-SIMPLIFIED_mergedIS.txt\n";
				}
			}
			elsif ( not scalar @htlv1_files ) {
				colored::stabilo( "red", "No htlv1 file found" ); # (*-clonalityResults.txt, *-mergedIS.txt, *-mergedIS.xls, *-SIMPLIFIED_mergedIS.txt)
				print "\n";
				die("No htlv1 file found for patient '$pat_name'") unless ($no_die);
			}
		}

		print "\n";
	}
	print "\n";
}

if (@patient_names) {
	print "Patient(s) not found in project(s) $project_names:\n";
	foreach my $name (@patient_names) {
		colored::stabilo( 'red', "\"$name\"" );
	}
	print "\n";
	die unless ($no_die);
}
print "Total: " . scalar @$files . " files\n\n" if ($files);

die("No file to copy\n") unless $files;



# todo: vérifier qu'il n'y ait pas déjà un fichier de ce nom dans le sftp ?

unless ($no_archive) {
	# Create the archive
	$archive_name = join( '_', @$project_names ) . '-' . join( '_', @file_types ) . '.tar' unless ($archive_name);
	$archive_name =~ s/\.tar\.gz$/\.tar/;
	$archive_name .= '.tar' unless ( $archive_name && $archive_name =~ /.tar$/ );
	$archive_name .= '.gz' if ( grep { $_ =~ /htlv1/i } @file_types );
	
	my $question = "Make an archive '$archive_name' of these files in '". abs_path($archive_dir) . "/'";
	$question .= " then put it on the sftp in '$dir_sftp'" if ($transfert_sftp);
	$question .= " ?  (y/n) ";
	my $choice = prompt($question);
	die("") if ( $choice ne "y" );
	
	my $cmd = "tar -cvf $archive_dir$archive_name " . join( " \\\n", @$files );
	$cmd =~ s/-cvf/-cvzf/ if ( grep { $_ =~ /htlv1/i } @file_types );
	print "\n";
	
	#warn $cmd;
	system $cmd;
	confess('Error while making the archive. The archive was not created.'."\n".$cmd) unless ( -e "$archive_dir$archive_name" );
	
	unless ($transfert_sftp) {
		colored::stabilo( 'white', "--------DONE--------" );
		colored::stabilo( 'white', "Archive created: '$archive_dir$archive_name'" );
		exit;
	}
	
	if ($transfert_sftp) {
	
		# Put the archive on the sftp
		$sftp = Net::SFTP->new(
			$sftp_host,
			user     => $sftp_user,
			key_path => $sftp_key_path
		);
		print "Putting the archive '$archive_name' on the sftp:\n";
		$sftp->put( "$archive_dir$archive_name", "/$dir_sftp$archive_name", \&callback ) || confess( "Could not transfer the archive on the sftp: " . $sftp->error );
		print "Putting the files on the sftp...\t100%\n";
	
		# Checks that the archive has been put on sftp
		my @ls = $sftp->ls("/$dir_sftp");
		@ls = map { $_->{filename} } @ls;
		if ( grep { $_ eq $archive_name } @ls ) {
			system("rm $archive_dir$archive_name");
			colored::stabilo( 'white', "--------DONE--------" );
			colored::stabilo( 'white', "Archive '$archive_name' put on the sftp in '$dir_sftp'" );
			exit;
		}
		else {
			colored::stabilo( 'red', 'ERROR' );
			confess("Archive '$archive_name' not put on the sftp");
		}
	}
	
}


if ($transfert_sftp and $no_archive ) {
	$archive_name = join( '_', @$project_names ) unless ($archive_name);
	$archive_name =~ s/\.tar\.gz$//;
	$archive_name =~ s/\.tar$//;
	$archive_name .= '/' unless ($archive_name =~ /\/$/ and $archive_name);
	
	my $question = "Put these files on the sftp in '$dir_sftp$archive_name'";
	$question .= " ?  (y/n) ";
	my $choice = prompt($question);
	die if ( $choice ne "y" );
	
	$sftp = Net::SFTP->new(
		$sftp_host,
		user     => $sftp_user,
		key_path => $sftp_key_path,
		debug => 1,
	);
	$sftp->do_mkdir("$archive_dir$archive_name") or confess("Impossible 'mkdir $archive_dir$archive_name' on sftp: ".$sftp->error);

	my $i = 0;
	for my $file (@$files) {
		printf( "%d / %d", $i, scalar @$files );
		$sftp->put( "$archive_dir$archive_name", $file ) or confess("Could not transfer file '$file' on the sftp: " . $sftp->error );
	}
	print "\n";

	# Checks that the files have been put on sftp
	my @ls = $sftp->ls("/$dir_sftp$archive_name");
	@ls = map { $_->{filename} } @ls;
	my $search = '('.join('$|',@$files).')$';
	my $nb_ok = grep ( /$search/, @ls );
	if ( $nb_ok == scalar @$files ) {
		colored::stabilo( 'white', "--------DONE--------" );
		colored::stabilo( 'white', "Files put on the sftp in '$dir_sftp$archive_name'" );
		exit;
	}
	else {
		colored::stabilo( 'red', 'ERROR' );
		confess("Missing files on the sftp");
	}

}

sub usage {
	print "
transfert_files
---------------
Mandatory arguments
	-project <s>		project names separated with a comma
	
Optional arguments
	-files <s>		file types separated with a comma (fastq, bam, vcf, htlv1 or all (all=fastq,bam,vcf))
	-patients <s>		patient names separated with a comma
	-archive_dir <s>	path to create the archive before putting it on the sftp, default '/data-pure/download/transfert_file/'
	-archive_name <s>	archive name
	-sftp / nosftp		whether to put or not the archive on the sftp
	-no_die			do not die if file(s) or patient(s) are not found
	-no_archive		with -sftp, do not create an archive but put the files individually in the sftp
	-help				

";
	exit(1);
}

sub callback {
	my ( $sftp, $data, $offset, $size ) = @_;
	my $percent = int( $offset / $size * 100 );
	printf( "Putting the file(s) on the sftp...\t%d%%\r", $percent );
}


