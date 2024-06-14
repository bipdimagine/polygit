use strict; 
use Data::Dumper;
use List::Util qw(max);
use Getopt::Long;
use Parallel::ForkManager;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo";
use lib "$Bin/../../GenBo/lib/obj-nodb";
use GBuffer;
use File::Basename;


#def variables
my $infile; 
my $force;
my $fork = 1;
my @patient_name; 
my @fastq_patient;
my $project_name;

GetOptions(
	'project|p=s' => \$project_name, 
	'infile|i=s' => \$infile, 
	'fork|f=s' => \$fork, 
	'force=s' => \$force, 
);

#die ("\n\nERROR: option -infile is mandatory. Die....\n\n") if (not $infile);
#die("\n\nERROR: -infile does not exist. Die....\n\n") unless (-e $infile);
 
my $buffer = new GBuffer;

my $h;
if ($project_name) {
	my $project = $buffer->newProject( -name => $project_name );
	foreach my $run (@{$project->getRuns()}) {
		my $out_dir = $run->demultiplex_dir."/fastq_screen/";
		if (not -d $out_dir) {
			my $cmd_dossier = "mkdir $out_dir"; 
			system($cmd_dossier);  
		}
		foreach my $patient (@{$run->getPatients()}) {
			my $patient_name = $patient->name();
			my $path_seq = $patient->getSequencesDirectory();
			my $fastq_file = select_fastq_file($patient_name, $path_seq);
			$h->{$patient_name}->{fastq} = $fastq_file;
			$h->{$patient_name}->{outdir} = $out_dir;
		}
	}
}
elsif ($infile)  {
	my $out_dir = $buffer->config->{project_data}->{demultiplex}.'/'.$infile.'/';
	if (not -d $out_dir) {
		die("\n\n$out_dir doesn't exists... DIE !!! \n\n");
	}
	#open(FILE, "$infile"); 
	#
	#while (<FILE>) {
	#	
	#	chomp($_); 
	#	if ($_ =~ "output_dir"){
	#		my @split = split(":", $_);
	#		$out_dir = $split[1]; 
	#	}
	#	
	#	else{
	#		my @split = split(":", $_);   
	#		push(@patient_name, $split[0]);
	#		push(@fastq_patient, $split[1]);		
	#	}
	#}
	#
	#close(FILE); 
	#
	#my $path = "$out_dir" . "fastq_screen_launch.sh";
	$out_dir = "$out_dir/fastq_screen/";
	if (not -d $out_dir) {
		my $cmd_dossier = "mkdir $out_dir"; 
		system($cmd_dossier);  
	}
	$h = fastq_file($infile);
	foreach my $pat_name(keys %$h) { $h->{$pat_name}->{outdir} = $out_dir; }
}


my $pm = new Parallel::ForkManager($fork);
foreach my $patient_name (keys %$h) {
	my $pid = $pm->start and next;
	my $fastq = $h->{$patient_name}->{fastq};
	my $out_dir = $h->{$patient_name}->{outdir};
	if (-e $fastq) {
		my $cmd = "perl $Bin/methodes_fqs.pl -i $fastq -o $out_dir -c $Bin/fastq_screen.conf -n $patient_name -t 1";
		system($cmd);
	}
	$pm->finish();
}
$pm->wait_all_children();





sub fastq_file {
	my ($run_name) = @_ ; #"240527_NB501645_0851_AHV33YAFX5.NGS2024_7774"
	$run_name =~ s/\.NGS20.+//;
	my $h = $buffer->get_demultiplex_run_infos($run_name);
	die("\n\nNo infos found for run_name $run_name in DB... DIE !!! \n\n") if not $h or scalar keys %$h == 0;
	my %hash; 
	my $last_patient_name;
	foreach my $patient_name (keys %{$h}) {
		$last_patient_name = $patient_name;
		my $proj = $buffer->newProject( -name => $h->{$patient_name}->{'project_name'} );
		my $pat = $proj->getPatient($patient_name);
		my $path_seq = $pat->getSequencesDirectory();
		my $fastq_file = select_fastq_file($patient_name,$path_seq );
		$hash{$patient_name}{fastq} = $fastq_file;
	}
	return \%hash; 
}

sub select_fastq_file {
	my ($patient_name, $path_seq) = @_;
	my $fastq_file;
	opendir my ($dir), $path_seq;
	my @found_files = readdir $dir;
	closedir $dir;
	my (@lFiles);
	foreach my $file (@found_files) {
		next if $file eq '.';
		next if $file eq '..';
		my $file_2 = $path_seq.'/'.$patient_name;
		my $regexp = $patient_name.'_S.+_R2_[L]?001.fastq.gz';
		if ($file =~ /$regexp/) {
			$fastq_file = $path_seq.'/'.$file;
			last if -e $fastq_file;
		}
		my $regexp_RC = $patient_name.'_RC_S.+_R2_[L]?001.fastq.gz';
		if ($file =~ /$regexp_RC/) {
			$fastq_file = $path_seq.'/'.$file;
			last if -e $fastq_file;
		}
	}
	return $fastq_file;
}





