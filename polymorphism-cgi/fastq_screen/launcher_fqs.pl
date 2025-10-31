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
use Carp;


#def variables
my $infile; 
my $force;
my $fork = 1;
my @patient_name; 
my @fastq_patient;
my $project_name;
my $patient_name_spec;

GetOptions(
	'project=s' => \$project_name, 
	'patient=s' => \$patient_name_spec, 
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
		foreach my $patient (@{$run->getPatients()}) {
			my $patient_name = $patient->name();
			next if $patient_name_spec and $patient_name ne $patient_name_spec; 
			my $fastq_file = get_fastq_file($patient);
			$h->{$patient_name}->{fastq} = $fastq_file;
			$h->{$patient_name}->{outdir} = $project->fastq_screen_path();
			if (not $fastq_file or not -e $fastq_file) {
				my $fileout = $patient->fastq_screen_path().'/'.$patient_name.'_screen_nom_espece.txt';
				open (F, ">$fileout");
				print F "no_fastq_found";
				close (F);
			}
		}
	}
}
elsif ($infile)  {
	$h = get_fastq_file_from_run($infile);
}


my $pm = new Parallel::ForkManager($fork);
foreach my $patient_name (keys %$h) {
	my $pid = $pm->start and next;
	my $fastq = $h->{$patient_name}->{fastq};
	my $out_dir = $h->{$patient_name}->{outdir};
	if (-e $fastq) {
		my $cmd = "perl $Bin/methodes_fqs.pl -i $fastq -o $out_dir -c $Bin/fastq_screen.conf -n $patient_name -t 1";
		system($cmd);
		if ($project_name) {
			my $b = new GBuffer;
			my $p = $b->newProject( -name => $project_name );
			my $patient = $p->getPatient($patient_name);
			my $run = $patient->getRun();
			my $fileout = $run->demultiplex_dir().'/fastq_screen/fastq_screen_'.$patient_name.'/'.$patient_name.'_screen_nom_espece.txt';
			my $fileout_error = $fileout.'.error';
			`rm $fileout_error` if -e $fileout_error;
			confess ("\n\nERROR: $fileout not found for $patient_name. Die...\n\n") if not -e $fileout;
			open (F, $fileout);
			my $specie_found = <F>;
			chomp($specie_found);
			close (F);
			my $h_db = $buffer->get_demultiplex_run_infos($run->run_name());
			my $specie_db = lc($h_db->{$patient_name}->{specie});
			if (lc($specie_db) ne lc($specie_found)) {
				`mv $fileout $fileout_error`;
				warn "\n\n";
				warn 'Patient: '.$project_name;
				warn 'Patient: '.$patient_name;
				warn Dumper $h_db;
				warn "\n\nERROR: no same specie found DB:$specie_db - FastqScreen:$specie_found. Die...\n\n";
			}
		}
	}
	$pm->finish();
}
$pm->wait_all_children();




sub get_fastq_file {
	my ($patient) = @_;
	foreach my $h_fastq (@{$patient->fastqFiles()}) {
		return $h_fastq->{R2} if exists $h_fastq->{R2} and -e $h_fastq->{R2};
	}
	return;
}


sub get_fastq_file_from_run {
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
		my $fastq_file = get_fastq_file($pat);
		$hash{$patient_name}{fastq} = $fastq_file;
		$hash{$patient_name}{outdir} = $proj->fastq_screen_path();
	}
	return \%hash; 
}






