package bds_steps;

use Term::ReadKey;
use IO::Prompt;
use Moo;

#use Data::Printer;
use FindBin qw($Bin);
use Time::Local;
use POSIX qw(strftime);
use Term::ANSIColor;
use Data::Dumper;
use job_bds;
use job_bds_tracking;
use sample;
use Carp;
use Digest::MD5 qw(md5 md5_hex md5_base64);
extends(qw(bds_root));

###########################################
# Variables globales pour les paths, id et versions des programmes #
###########################################

has "again" => (
	is      => 'rw',
	default => sub { undef },
);

has 'fastq_files' => (
	is      => 'rw',
	default => sub { {} },
);

#########
# Méthodes  #
#########

sub job_super_types () {
	confess();

#my (@types) = sort{$self->{jobs_super_type}->{$a} <=> $self->{jobs_super_type}->{$b}} keys %{$self->{jobs_super_type}};
#return \@types;
}

###########################
# Méthodes liées aux étapes  du pipeline #
###########################

###méthode d'aln mem utilisée pour le pipeline miseq, miseq_primer (avec trimming)

###méthode utilisée dans le pipeline illumina (sans trimming)

sub run_alignment_frag {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $name = $self->patient()->name();

	#confess("check the code man");
	my ($dirin) = $self->patient()->getSequencesDirectory();
	my $method  = $self->patient()->alignmentMethod();
	my $dirout  = $self->project->getAlignmentPipelineDir("$method");
	$filein =
	  file_util::find_fragment_file( $self->project, $name, $dirin, "gz" );

	#my $ppn = 16;
	my $ppn = $self->nproc;    #if $self->nocluster;

	my $project_name = $self->project->name();
	my $fileout      = $dirout . $name . ".align.bam";
	my $type         = "align-frag";
	my $stepname     = $self->patient->name . "@" . $type;
	my $bin_dev      = $self->script_dir();
	my $cmd =
qq{perl $bin_dev/align.pl -file1=$filein -method=$method   -mode=frag -project=$project_name -name=$name -bam=$fileout -fork=$ppn};

#my $job_bds = job_bds_tracking->new(uuid=>$self->bds_uuid,cmd=>["$cmd"],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds,software=>$synonym_program->{$method},sample_name=>$self->patient->name(),project_name=>$self->patient->getProject->name);

	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);

#
#	$self->hash_steps->{$stepname}->{id}= 2;
#	my ($job_hisat,$job_final) = $self->construct_jobs(stepname=>$stepname,cmd=>[$cmd],ppn=>$proc);
#	return ([$job_hisat],[$bam]);
}

sub alignment {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};
	my $name   = $self->patient()->name();
	my $method = $self->patient()->alignmentMethod();
	my $run    = $self->patient->getRun();
	if ( $run->infosRun->{method} eq "fragment" ) {
		return $self->run_alignment_frag({ filein => $filein} );
	}
	return $self->run_alignment_pe( {filein => $filein} );
}

sub run_alignment_pe {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $method       = $self->patient()->alignmentMethod();
	if ($method eq "star"){
		return $self->star_align($hash);
	}
	my $name         = $self->patient()->name();
	my ($dirin)      = $self->patient()->getSequencesDirectory();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	$dirout .= "/" . $self->patient()->name() . "/";
	my $run  = $self->patient->getRun();
	my $type = $run->infosRun->{method};
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;

	#my $ppn = 10;
	my $ppn = $self->nproc;    # if $self->nocluster;

	#	my $split = "_";
	#	my $ext1 = $split."1".$split."C";
	#	my $ext2 = $ext1;
	#	$ext2 =~ s/1/2/;
	#my $files_pe1 = file_util::test($self->patient);
	my $files_pe1 = file_util::find_file_pe( $self->patient, "" );

#my $files_pe1 = file_util::my $files_pe1 = file_util::find_file_pe($self->patient,$ext1);($self->patient->name,$dirin,$ext1);
#p $files_pe1;
	my $count_lane = scalar(@$files_pe1);
	my $nb_bam     = 1;
	my $already    = 0;
	my @cmds;
	my @bams;
	my $files;
	my @jobs;
	print $self->patient->name() . ":\n";

	foreach my $cp (@$files_pe1) {
		my $file1 = $cp->{R1};
		my $file2 = $cp->{R2};
		print "\t $nb_bam : " . $cp->{R1} . " " . $cp->{R2} . "\n";

		if (   exists $self->fastq_files->{$file1}
			or exists $self->fastq_files->{$file2} )
		{
			die(
"same fastq file present in two different patient : $name $file1 $file2"
			);
		}
		$self->fastq_files->{$file1}++;
		$self->fastq_files->{$file2}++;
		die( "problem $file1 $file2 $dirin :" . $name )
		  unless -e $dirin . $file2;
		die( "problem $file1 $file2 $dirin :" . $name )
		  unless -e $dirin . $file1;
		my $bam = $dirout . $name . ".F$nb_bam." . "bwa.bam";
		my $f1  = $dirin . $file1;
		my $f2  = $dirin . $file2;

		my $bin_dev = $self->script_dir();

		my $cmd =
qq{perl $bin_dev/align.pl -file1=$f1 -file2=$f2 -method=$method -lane=$nb_bam  -project=$project_name -name=$name -bam=$bam -fork=$ppn };
		my $type     = "$method#" . $nb_bam;
		my $stepname = $self->patient->name . "@" .$method;
		my $job_bds = job_bds_tracking->new(
			uuid         => $self->bds_uuid,
			cmd          => [$cmd],
			name         => $stepname,
			ppn          => $ppn,
			filein       => [ $f1, $f2 ],
			fileout      => $bam,
			type         => $type,
			dir_bds      => $self->dir_bds,
			sample_name  => $name,
			project_name => $project_name,
			software     => $method
		);
		$self->current_sample->add_job( { job => $job_bds } );
		push( @jobs, $job_bds );
		$nb_bam++;

		if ( $self->unforce()
			&& ( -e $bam or -e $self->patient()->getBamFileName() ) )
		{
			$job_bds->skip();
		}
		push( @bams, $bam );
	}

	my ($fileout) = $self->merge_bam( {filein => \@bams} );
	return ($fileout);

}

sub bwa2 {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $method       = "bwa2";
	my $name         = $self->patient()->name();
	my ($dirin)      = $self->patient()->getSequencesDirectory();
	my $project_name = $self->project->name();
	my $run          = $self->patient->getRun();
	my $ppn = 40;    # if $self->nocluster;

	my $files_pe1  = file_util::find_file_pe( $self->patient, "" );
	my $count_lane = scalar(@$files_pe1);
	my $nb_bam     = 1;
	my $already    = 0;
	my @cmds;
	my @bams;
	my $files;
	my @jobs;
	print $self->patient->name() . ":\n";
	my @af1;
	my @af2;

	foreach my $cp (@$files_pe1) {
		push( @af1, $dirin . "/" . $cp->{R1} );
		push( @af2, $dirin . "/" . $cp->{R2} );
	}
	my $fileout =
		$self->patient()->project->getAlignmentPipelineDir("bwa2") . "/"
	  . $name
	  . ".align.bam";
	my $bin_dev = $self->script_dir();
	my $f1      = join( ",", @af1 );
	my $f2      = join( ",", @af2 );
	my $version = $self->patient()->project->genome_version();

	my $cmd =
qq{cd $bin_dev;perl $bin_dev/bwa2/bwa2.pl -file1=$f1 -file2=$f2  -bam=$fileout -project=$project_name -patient=$name -fork=$ppn -version=$version};
	open( OUT, ">>/home/bds/list.txt" );
	print OUT $cmd . "\n";
	close OUT;
	my $dir_out =
	  $self->patient()->project->getCallingPipelineDir("genome") . "/$name";
	system("mkdir -p $dir_out");

	my $type     = "bwa2";
	my $stepname = $self->patient->name . "@" . $type;

	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [ $f1, $f2 ],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		sample_name  => $name,
		project_name => $project_name,
		software     => $method
	);
	$self->current_sample->add_job( { job => $job_bds } );
	push( @jobs, $job_bds );
	$nb_bam++;

	if ( $self->unforce()
		&& ( -e $fileout or -e $self->patient()->getBamFileName ) )
	{
		$job_bds->skip();
	}

	return ($fileout);

}

sub elprep5_genome {
	my ( $self, $hash ) = @_;
	my $filein  = $hash->{filein};
	my $name    = $self->patient()->name();
	my $method  = "elprep";
	my $fileout = $filein;
	$fileout =~ s/bam/sort\.rmdup\.bam/;

	my $project      = $self->project;
	my $project_name = $project->name;
	my $bin_dev      = $self->script_dir;
	my $ppn          = 40;
	my $cmd =
qq{perl $bin_dev/elprep/elprep5.pl  -project=$project_name -patient=$name -bamin=$filein  -bamout=$fileout -fork=$ppn };

	my $type     = "gelprep";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		cmd          => ["$cmd"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "elprep5",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce()
		&& ( -e $fileout or -e $self->patient()->getBamFileName() ) )
	{
		$job_bds->skip();
	}
	else {

		open( OUT, ">>/home/bds/$project_name.elprep.txt" );
		print OUT $cmd . "\n";
		close OUT;
	}
	return ($fileout);
}

#sub elprep5_target {
my ( $self, $hash ) = @_;
my $filein = $hash->{filein};

#	my $project = $self->project;
#	my $name = $project->name();
#	my $arg = $self->argument_patient();
#	my $patient = $self->patient();
#	my $dir_out_vcf= $project->getCallingPipelineDir("haplotypecaller4");
#	my $bamin= $patient->getBamFileName();
#	my $dir_out= $self->project->getCallingPipelineDir($self->method_calling);
#	my $bamout= $dir_out."/".$patient->name()."_elprep5.bam" ;
#	my $fileout = $project->getGvcfDir."/haplotypecaller4/".$patient->name.".g.vcf.gz";
#	my $type = "elprep5";
#	my $stepname =$patient->name()."@".$type;
#	my $fork =$self->nproc;
#	my $bin_dev = $self->script_dir();
#	my $cmd = "perl $bin_dev/elprep/elprep5_target.pl -project=".$name." -gvcf=".$fileout." -patients=".$patient->name()." -bamin=".$bamin." -bamout=".$bamout ;
#	 my $job_bds = job_bds_tracking->new(uuid=>$self->bds_uuid,software=>"",sample_name=>$self->argument_patient(),project_name=>$self->patient->getProject->name,cmd=>[$cmd],name=>$stepname,ppn=>$fork,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
#
#	$self->current_sample->add_job({job=>$job_bds});
#	if ($self->unforce() && -e $fileout){
#	  	$job_bds->skip();
#	}
#	return ($fileout);
#}

sub elprep5_target {
	my ( $self, $hash ) = @_;
	my $filein  = $hash->{filein};
	my $name    = $self->patient()->name();
	my $method  = "elprep";
	my $fileout = $filein;
	$fileout =~ s/bam/sort\.rmdup\.bam/;

	my $project      = $self->project;
	my $project_name = $project->name;
	my $bin_dev      = $self->script_dir;
	my $ppn          = 20;
	my $cmd =
qq{perl $bin_dev/elprep/elprep5_target.pl  -project=$project_name -patient=$name -bamin=$filein  -bamout=$fileout -fork=$ppn };
	my $type     = "telprep";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		cmd          => ["$cmd"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "elprep5",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce()
		&& ( -e $fileout or -e $self->patient()->getBamFileName() ) )
	{
		$job_bds->skip();
	}
	else {

		open( OUT, ">>/home/bds/$project_name.elprep.txt" );
		print OUT $cmd . "\n";
		close OUT;
	}
	return ($fileout);
}

sub sortdedup {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my $fileout      = $filein;
	$fileout =~ s/bam/sort\.bam/;
	my $ppn     = 40;                    # if $self->nocluster;
	my $bin_dev = $self->script_dir();
	$method = "bamsormadup";
	my $cmd =
qq{cd $bin_dev;perl $bin_dev/biobambam/sortdedup.pl -bamin=$filein -bamout=$fileout  -bam=$fileout -project=$project_name -patient=$name -fork=$ppn };
	my $type     = "sortdedup";
	my $stepname = $self->patient->name . "@" . $type;

	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		sample_name  => $name,
		project_name => $project_name,
		software     => $method
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout
		or -e $self->patient()->getBamFileName() )
	{
		$job_bds->skip();
	}
	return ($fileout);
}

sub sparkSort {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my $fileout      = $filein;
	$fileout =~ s/bam/sort\.bam/;
	my $ppn     = 40;                    # if $self->nocluster;
	my $bin_dev = $self->script_dir();

	my $cmd =
qq{cd $bin_dev;perl $bin_dev/spark_picard/sort.pl -bamin=$filein -bamout=$fileout  -bam=$fileout -project=$project_name -patient=$name -fork=$ppn };

	my $type     = "sparkSort";
	my $stepname = $self->patient->name . "@" . $type;

	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		sample_name  => $name,
		project_name => $project_name,
		software     => $method
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout
		or -e $self->patient()->getBamFileName() )
	{
		$job_bds->skip();
	}
	return ($fileout);
}

sub sparkDedup {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my $fileout      = $filein;
	$fileout =~ s/bam/dedup\.bam/;
	my $ppn     = 40;                    # if $self->nocluster;
	my $bin_dev = $self->script_dir();

	my $cmd =
qq{cd $bin_dev;perl $bin_dev/spark_picard/dedup.pl -bamin=$filein -bamout=$fileout  -bam=$fileout -project=$project_name -patient=$name -fork=$ppn };
	my $type     = "sparkdedup";
	my $stepname = $self->patient->name . "@" . $type;

	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		sample_name  => $name,
		project_name => $project_name,
		software     => $method
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout
		or -e $self->patient()->getBamFileName() )
	{
		$job_bds->skip();
	}
	return ($fileout);
}

sub pipeline_genome {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my ($dirin)      = $self->patient()->getSequencesDirectory();
	my $project_name = $self->project->name();
	my $run          = $self->patient->getRun();
	my $ppn = 40;    # if $self->nocluster;

	my $files_pe1  = file_util::find_file_pe( $self->patient, "" );
	my $count_lane = scalar(@$files_pe1);
	my $nb_bam     = 1;
	my $already    = 0;
	my @cmds;
	my @bams;
	my $files;
	my @jobs;
	print $self->patient->name() . ":\n";
	my @af1;
	my @af2;

	foreach my $cp (@$files_pe1) {
		push( @af1, $dirin . "/" . $cp->{R1} );
		push( @af2, $dirin . "/" . $cp->{R2} );
	}
	my $bin_dev = $self->script_dir();
	my $f1      = join( ",", @af1 );
	my $f2      = join( ",", @af2 );
	my $cmd =
qq{cd $bin_dev;perl $bin_dev/genome/bwa2.pl -file1=$f1 -file2=$f2  -dir=$dirin -project=$project_name -patient=$name -fork=$ppn };
	my $dir_out =
	  $self->patient()->project->getCallingPipelineDir("genome") . "/$name";
	system("mkdir -p $dir_out");

	my $fileout =
		$self->patient()->project->getAlignmentPipelineDir("bwa2") . "/"
	  . $name
	  . ".align.bam";
	my $type     = "genome";
	my $stepname = $self->patient->name . "@" . $type;

	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [ $f1, $f2 ],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		sample_name  => $name,
		project_name => $project_name,
		software     => $method
	);
	$self->current_sample->add_job( { job => $job_bds } );
	push( @jobs, $job_bds );
	$nb_bam++;

	if ( $self->unforce() && -e $fileout )
	{    # or -e  $self->patient()->getBamFileName())){
		$job_bds->skip();
	}
	return ($fileout);

}

sub run_alignment_umi {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my ($dirin)      = $self->patient()->getSequencesDirectory();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn = $self->nproc;    # if $self->nocluster;
	my $files_pe1  = file_util::find_file_pe_umi( $self->patient, "" );
	my $count_lane = scalar(@$files_pe1);
	my $nb_bam     = 1;
	my $already    = 0;
	my @cmds;
	my @bams;
	my $files;
	my @jobs;
	print $self->patient->name() . ":\n";
	my $file1;
	my $file2;
	my $R1;
	my $R3;

	foreach my $cp (@$files_pe1) {
		$file1 = $cp->{R1};
		$file2 = $cp->{R3};
		die( "problem $file1 $file2 $dirin :" . $name )
		  unless -e $dirin . $file2;
		die( "problem $file1 $file2 $dirin :" . $name )
		  unless -e $dirin . $file1;
		my $bam     = $dirout . $name . ".F$nb_bam." . "bwa.bam";
		my $f1      = $dirin . $file1;
		my $f2      = $dirin . $file2;
		my $bin_dev = $self->script_dir();

		my $cmd =
qq{perl $bin_dev/align.pl -file1=$f1 -file2=$f2 -method=$method -lane=$nb_bam  -project=$project_name -name=$name -bam=$bam -fork=$ppn };
		my $type     = "align_umi#" . $nb_bam;
		my $stepname = $self->patient->name . "@" . $type;
		my $job_bds  = job_bds_tracking->new(
			uuid         => $self->bds_uuid,
			cmd          => [$cmd],
			name         => $stepname,
			ppn          => $ppn,
			filein       => [ $f1, $f2 ],
			fileout      => $bam,
			type         => $type,
			dir_bds      => $self->dir_bds,
			sample_name  => $name,
			project_name => $project_name,
			software     => $method
		);
		$self->current_sample->add_job( { job => $job_bds } );
		push( @jobs, $job_bds );
		$nb_bam++;

		if ( $self->unforce() && -e $bam ) {
			$job_bds->skip();
		}
		push( @bams, $bam );
	}

	#	#my ($fileout) = $self->elprep_all(filein=>\@bams);
	my ($fileout) = $self->merge_bam( {filein => \@bams} );

	#my ($fileout) =$self->sort_sam_umi(filein=>\@bams);
	return ($fileout);
}

sub agent_trimmer {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my ($dirin)      = $self->patient()->getSequencesDirectory();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn = $self->nproc;    # if $self->nocluster;

	my @r1;
	my @r2;
	my $files_pe1 = file_util::find_file_pe( $self->patient, "" );
	foreach my $cp (@$files_pe1) {
		my $file1 = $cp->{R1};
		my $file2 = $cp->{R2};
		die( "problem $file1 $file2 $dirin :" . $name )
		  unless -e $dirin . $file2;
		die( "problem $file1 $file2 $dirin :" . $name )
		  unless -e $dirin . $file1;
		push( @r1, $dirin . $file1 );
		push( @r2, $dirin . $file2 );
	}

	my $f1 = join( ",", @r1 );
	my $f2 = join( ",", @r2 );

	my $fileout = $dirout . "/" . $name . "_RN.txt.gz";
	my $bin_dev = $self->script_dir();

	my $cmd =
qq{perl $bin_dev/xt_hs2.pl -file1=$f1 -file2=$f2 -method=$method -lane=1  -project=$project_name -name=$name -fileout=$fileout -fork=$ppn };
	my $type     = "agent_trimmer#";
	my $stepname = $self->patient->name . "@" . $type;

	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [ $f1, $f2 ],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		sample_name  => $name,
		project_name => $project_name,
		software     => $method
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);

}

sub flexbar {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my ($dirin)      = $self->patient()->getSequencesDirectory();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);

	my $illumina_adaptors = $self->project->getIlluminaAdaptors();
	my $tso_adaptors      = $self->project->getTsoAdaptors();
	my $fileout           = $dirout . "/" . $name . ".flexbarOut.log";

	#system("mkdir -p $dirout/metrics;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn    = $self->nproc;    # if $self->nocluster;
	my $nb_bam = 1;
	my @jobs;
	my @fastq_out;
	my $files_pe1 = file_util::find_file_pe( $self->patient, "" );
	foreach my $cp (@$files_pe1) {
		my $file1 = $cp->{R1};
		my $file2 = $cp->{R2};
		die( "problem $file1 $file2 $dirin :" . $name )
		  unless -e $dirin . $file2;
		die( "problem $file1 $file2 $dirin :" . $name )
		  unless -e $dirin . $file1;
		my $f1         = $dirin . $file1;
		my $f2         = $dirin . $file2;
		my $bin_dev    = $self->script_dir();
		my $fastq_out1 = $dirout . "/" . $name . "_F" . $nb_bam . "_1.fastq.gz";
		my $fastq_out2 = $dirout . "/" . $name . "_F" . $nb_bam . "_2.fastq.gz";

		my $cmd =
qq{perl $bin_dev/flexbar.pl -file1=$f1 -file2=$f2 -method=$method -lane=$nb_bam  -project=$project_name -name=$name -fork=$ppn -illumina=$illumina_adaptors -tso=$tso_adaptors };
		my $type     = "flexbar#" . $nb_bam;
		my $stepname = $self->patient->name . "@" . $type;
		my $job_bds  = job_bds_tracking->new(
			uuid         => $self->bds_uuid,
			cmd          => [$cmd],
			name         => $stepname,
			ppn          => $ppn,
			filein       => [ $f1, $f2 ],
			fileout      => $fastq_out1,
			type         => $type,
			dir_bds      => $self->dir_bds,
			software     => $method,
			sample_name  => $name,
			project_name => $project_name
		);

#my $job_bds = job_bds_tracking->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$f1,$f2],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds,sample_name=>$name,project_name=>$project_name,software=>$method);
		$self->current_sample->add_job( { job => $job_bds } );
		push( @jobs, $job_bds );
		$nb_bam++;
		if ( $self->unforce() && ( -e $fastq_out1 ) ) {
			$job_bds->skip();
		}

		#push(@fastq_out,$fastq_out1);
		#push(@fastq_out,$fastq_out2);
	}

	return ($fileout);
}

sub run_alignment_flexbar {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn = $self->nproc;    # if $self->nocluster;
	my $files_pe1  = file_util::find_file_in_dir( $self->patient, "", $dirout );
	my $count_lane = scalar(@$files_pe1);
	my $nb_bam     = 1;
	my $already    = 0;
	my @cmds;
	my @bams;
	my $files;
	my @jobs;
	print $self->patient->name() . ":\n";

	foreach my $cp (@$files_pe1) {
		my $file1 = $cp->{R1};
		my $file2 = $cp->{R2};
		print "\t $nb_bam : " . $cp->{R1} . " " . $cp->{R2} . "\n";

		if (   exists $self->fastq_files->{$file1}
			or exists $self->fastq_files->{$file2} )
		{
			die(
"same fastq file present in two different patient : $name $file1 $file2"
			);
		}
		$self->fastq_files->{$file1}++;
		$self->fastq_files->{$file2}++;
		die( "problem $file1 $file2 $dirout :" . $name )
		  unless -e $dirout . $file2;
		die( "problem $file1 $file2 $dirout :" . $name )
		  unless -e $dirout . $file1;
		my $bam = $dirout . $name . ".F$nb_bam." . "bwa.bam";
		my $f1  = $dirout . $file1;
		my $f2  = $dirout . $file2;

		my $bin_dev = $self->script_dir();

		my $cmd =
qq{perl $bin_dev/align.pl -file1=$f1 -file2=$f2 -method=$method -lane=$nb_bam  -project=$project_name -name=$name -bam=$bam -fork=$ppn };
		my $type     = "align#" . $nb_bam;
		my $stepname = $self->patient->name . "@" . $type;
		my $job_bds  = job_bds_tracking->new(
			uuid         => $self->bds_uuid,
			cmd          => [$cmd],
			name         => $stepname,
			ppn          => $ppn,
			filein       => [ $f1, $f2 ],
			fileout      => $bam,
			type         => $type,
			dir_bds      => $self->dir_bds,
			software     => $method,
			sample_name  => $name,
			project_name => $project_name
		);

#		my $job_bds = job_bds_tracking->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$f1,$f2],fileout=>$bam,type=>$type,dir_bds=>$self->dir_bds,sample_name=>$name,project_name=>$project_name,software=>$method);
		$self->current_sample->add_job( { job => $job_bds } );
		push( @jobs, $job_bds );
		$nb_bam++;
		if ( $self->unforce()
			&& ( -e $bam or -e $self->patient()->getBamFileName() ) )
		{
			$job_bds->skip();
		}
		push( @bams, $bam );
	}

	#my ($fileout) = $self->elprep_all(filein=>\@bams);
	my ($fileout) = $self->merge_bam( {filein => \@bams} );
	return ($fileout);

}

sub run_alignment_xths2 {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn      = $self->nproc;                          # if $self->nocluster;
	my $file1    = $dirout . "/" . $name . "_R1.fastq.gz";
	my $file2    = $dirout . "/" . $name . "_R2.fastq.gz";
	my $mbc_file = $dirout . "/" . $name . "_RN.txt.gz";

	my $bam     = $dirout . $name . ".xths2.bam";
	my $bin_dev = $self->script_dir();

	my $cmd =
qq{perl $bin_dev/align.pl -file1=$file1 -file2=$file2 -method=$method -lane=1 -mode=xths2  -project=$project_name -name=$name -bam=$bam -fork=$ppn };
	my $type     = "align_xths2#";
	my $stepname = $self->patient->name . "@" . $type;

	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [ $file1, $file2 ],
		fileout      => $bam,
		type         => $type,
		dir_bds      => $self->dir_bds,
		sample_name  => $name,
		project_name => $project_name,
		software     => $method
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $bam ) {
		$job_bds->skip();
	}
	return ($bam);

}

sub run_alignment_deepseq {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn   = $self->nproc;    # if $self->nocluster;
	my $file1 = $filein;
	my $file2 = $file1;
	$file2 =~ s/_R1/_R2/;
	my $bam     = $dirout . $name . "_umi_extracted_aligned.bam";
	my $bin_dev = $self->script_dir();
	my $cmd =
qq{perl $bin_dev/align.pl -file1=$file1 -file2=$file2 -method=$method -lane=1 -mode=pe  -project=$project_name -name=$name -bam=$bam -fork=$ppn };
	my $type     = "align_deepseq#";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [ $file1, $file2 ],
		fileout      => $bam,
		type         => $type,
		dir_bds      => $self->dir_bds,
		sample_name  => $name,
		project_name => $project_name,
		software     => $method
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $bam ) {
		$job_bds->skip();
	}
	return ($bam);

}

sub agent_locatit_duplex {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn = $self->nproc;    # if $self->nocluster;

	#	my @umi_files =  glob($dirout."/".$name."*_RN.txt.gz");
	#	my $mbc_file= $umi_files[0];
	my $agent_locatit =
	  "/software/distrib/AGeNT/AGeNT-latest/lib/locatit-2.0.5.jar";
	my $output_dir = $self->project->getAlignmentPipelineDir($method);
	my $fileout    = $dirout . $name . ".bam";
	my $java       = $self->project->getSoftware('java');
	$java = "java" unless -e $java;
	my $cmd =
qq{$java -Xmx12G -jar $agent_locatit -S -v2Duplex -d 1 -m 3 -q 25 -Q 25 -o $fileout $filein };
	my $type     = "agent_locatit";
	my $stepname = $name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		sample_name  => $name,
		project_name => $project_name,
		software     => $method
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}

	return ($fileout);
}

sub agent_locatit_single {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn      = $self->nproc;                          # if $self->nocluster;
	my $mbc_file = $dirout . "/" . $name . "_RN.txt.gz";
	my $agent_locatit =
	  "/software/distrib/AGeNT/AGeNT-latest/lib/locatit-2.0.5.jar";
	my $output_dir = $self->project->getAlignmentPipelineDir($method);
	my $fileout    = $dirout . $name . ".bam";
	my $java       = $self->project->getSoftware('java');
	$java = "java" unless -e $java;
	my $cmd =
qq {$java -Xmx12G -jar $agent_locatit -S  -v2 -d 1 -m 3 -q 25  -Q 25 PM:xm,Q:xq,q:nQ,r:nR   -IB -OB  -o $fileout  $filein };

#my $cmd = qq{$java -Xmx12G -jar $agent_locatit -S  -PM:xm,Q:xq,q:nQ,r:nR    -q 25 -m 1 -U -IS -OB -C -i -r  -o $fileout  $filein $mbc_file};
	my $type     = "agent_locatit";
	my $stepname = $name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		sample_name  => $name,
		project_name => $project_name,
		software     => $method
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}

	return ($fileout);
}

sub agent_locatit_hybrid {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn      = $self->nproc;                          # if $self->nocluster;
	my $mbc_file = $dirout . "/" . $name . "_RN.txt.gz";
	my $agent_locatit =
	  "/software/distrib/AGeNT/AGeNT-latest/lib/locatit-2.0.5.jar";
	my $output_dir = $self->project->getAlignmentPipelineDir($method);
	my $fileout    = $dirout . $name . ".bam";
	my $java       = $self->project->getSoftware('java');
	$java = "java" unless -e $java;

#my $cmd= qq {$java -Xmx12G -jar $agent_locatit -S  -v2 -d 1 -m 3 -q 25  -Q 25 PM:xm,Q:xq,q:nQ,r:nR   -IB -OB  -o $fileout  $filein };
#my $cmd = qq{$java -Xmx12G -jar $agent_locatit -S  -PM:xm,Q:xq,q:nQ,r:nR    -q 25 -m 1 -U -IS -OB -C -i -r  -o $fileout  $filein $mbc_file};
	my $cmd =
qq {$java -Xmx12G -jar $agent_locatit -S  -v2 -d 1 -m 3 -q 25  -Q 25 PM:xm,Q:xq,q:nQ,r:nR   -IB -OB  -o $fileout  $filein };
	my $type     = "agent_locatit";
	my $stepname = $name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		sample_name  => $name,
		project_name => $project_name,
		software     => $method
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}

	return ($fileout);
}

sub merge_picard {
	my ( $self, $hash ) = @_;
	my $filein   = $hash->{filein};
	my @allfiles = @$filein;
	my $name     = $self->patient->name();
	my $project  = $self->patient->getProject;

	#my ($f)  = File::Util->new();
	my $m      = $self->patient->alignmentMethod();
	my $suffix = $filein->[0];
	$suffix =~ s/$name//;
	$suffix = split( /\.sam/, $filein );

	my $fileout =
	  $project->getAlignmentPipelineDir($m) . "/" . $name . "_merge.bam";

#$fileout = $project->getAlignmentPipelineDir($m) . "/" . $name . "F[12345678].u.sort.bam" unless -e $fileout;
	my $merge_files = join( " I=", @allfiles );

	my $reference = $self->project->genomeFasta();
	my $java      = $self->project->getSoftware('java');
	$java = "java" unless -e $java;
	my $picard = $java . " -jar " . $self->project->getSoftware('picard_path');

	my $ppn = $self->nproc;
	$ppn = int( $self->nproc / 2 ) if $self->nocluster;
	my $cmd =
		$picard
	  . " MergeSamFiles CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate R="
	  . $reference . " I="
	  . $merge_files
	  . " OUTPUT="
	  . $fileout;
	my $type     = "merge-picard";
	my $stepname = $self->patient->name . "@" . $type;

	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => $filein,
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "picard",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub umi_ubam {
	my ( $self, $hash ) = @_;
	my $filein    = $hash->{filein};
	my $method    = $self->patient()->alignmentMethod();
	my $name      = $self->patient()->name();
	my ($dirin)   = $self->patient()->getSequencesDirectory();
	my $files_pe1 = file_util::find_file_pe_umi( $self->patient, "" );
	foreach my $cp (@$files_pe1) {

	}
}

sub fastq_to_bam {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my ($dirin)      = $self->patient()->getSequencesDirectory();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn  = $self->nproc;                          # if $self->nocluster;
	my $java = $self->project->getSoftware('java');
	$java = "java" unless -e $java;
	my $fgbio     = $java . " -jar " . $self->project->getSoftware('fgbio');
	my $files_pe1 = file_util::find_file_pe_umi( $self->patient, "" );
	my $nb_bam    = 1;
	my @bams;
	my @jobs;
	print $self->patient->name() . ":\n";

	foreach my $cp (@$files_pe1) {
		my $file1    = $cp->{R1};
		my $file2    = $cp->{R3};
		my $file_umi = $cp->{R2};
		print "\t $nb_bam : " . $cp->{R1} . " " . $cp->{R3} . "\n";

		if (   exists $self->fastq_files->{$file1}
			or exists $self->fastq_files->{$file2} )
		{
			die(
"same fastq file present in two different patient : $name $file1 $file2"
			);
		}
		$self->fastq_files->{$file1}++;
		$self->fastq_files->{$file2}++;
		$self->fastq_files->{$file_umi}++;
		die( "problem $file1 $file2 $dirin :" . $name )
		  unless -e $dirin . $file2;
		die( "problem $file1 $file2 $dirin :" . $name )
		  unless -e $dirin . $file1;
		my $bam = $dirout . $name . ".F$nb_bam." . "u.bam";
		my $f1  = $dirin . $file1;
		my $f2  = $dirin . $file2;

		my $cmd =
			$fgbio
		  . " FastqToBam -i  "
		  . $f1 . " "
		  . $f2 . " -o "
		  . $bam
		  . " --sample="
		  . $name
		  . " --library="
		  . $name
		  . " --read-group-id="
		  . $name;
		my $type     = "fastq_to_bam#" . $nb_bam;
		my $stepname = $self->patient->name . "@" . $type;

		my $job_bds = job_bds_tracking->new(
			uuid         => $self->bds_uuid,
			cmd          => [$cmd],
			name         => $stepname,
			ppn          => $ppn,
			filein       => [ $f1, $f2 ],
			fileout      => $bam,
			type         => $type,
			dir_bds      => $self->dir_bds,
			sample_name  => $name,
			project_name => $project_name,
			software     => $method
		);
		$self->current_sample->add_job( { job => $job_bds } );
		push( @jobs, $job_bds );
		$nb_bam++;

		if ( $self->unforce() && -e $bam ) {
			$job_bds->skip();
		}
		push( @bams, $bam );
	}
	my ($fileout) = $self->sort_sam_umi( filein => \@bams );
	return ($fileout);

	#	my ($fileout) = $self->merge_bam(filein=>\@bams);
	#	return ($fileout);
}

sub fastq_to_sam {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my ($dirin)      = $self->patient()->getSequencesDirectory();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn  = $self->nproc;                          # if $self->nocluster;
	my $java = $self->project->getSoftware('java');
	$java = "java" unless -e $java;
	my $picard = $java . " -jar " . $self->project->getSoftware('picard_path');

	#my $fgbio =  $java." -jar ".$self->project->getSoftware('fgbio');
	my $files_pe1 = file_util::find_file_pe( $self->patient, "" );
	my $nb_bam    = 1;
	my @bams;
	my @jobs;
	print $self->patient->name() . ":\n";
	foreach my $cp (@$files_pe1) {
		my $file1 = $cp->{R1};

		#my $file2 =  $cp->{R3};
		my $file2 = $cp->{R2};
		print "\t $nb_bam : " . $cp->{R1} . " " . $cp->{R2} . "\n";

		if (   exists $self->fastq_files->{$file1}
			or exists $self->fastq_files->{$file2} )
		{
			die(
"same fastq file present in two different patient : $name $file1 $file2"
			);
		}
		$self->fastq_files->{$file1}++;
		$self->fastq_files->{$file2}++;

		# $self->fastq_files->{$file_umi} ++;
		die( "problem $file1 $file2 $dirin :" . $name )
		  unless -e $dirin . $file2;
		die( "problem $file1 $file2 $dirin :" . $name )
		  unless -e $dirin . $file1;
		my $bam = $dirout . $name . ".F$nb_bam" . "_unmapped.sam";
		my $f1  = $dirin . $file1;
		my $f2  = $dirin . $file2;

		my $cmd =
			$picard
		  . " FastqToSam F1="
		  . $f1 . " F2="
		  . $f2 . " o="
		  . $bam
		  . " SM='umi_consensus'";
		my $type     = "fastq_to_sam#" . $nb_bam;
		my $stepname = $self->patient->name . "@" . $type;

		my $job_bds = job_bds_tracking->new(
			uuid         => $self->bds_uuid,
			cmd          => [$cmd],
			name         => $stepname,
			ppn          => $ppn,
			filein       => [ $f1, $f2 ],
			fileout      => $bam,
			type         => $type,
			dir_bds      => $self->dir_bds,
			sample_name  => $name,
			project_name => $project_name,
			software     => $method
		);
		$self->current_sample->add_job( { job => $job_bds } );
		push( @jobs, $job_bds );
		$nb_bam++;

		if ( $self->unforce() && -e $bam ) {
			$job_bds->skip();
		}
		push( @bams, $bam );
	}

	#	my ($fileout) = $self->sort_sam_umi(filein=>\@bams);
	#	return ($fileout);
	my ($fileout) = $self->extract_umi_from_bam( filein => \@bams );
	return ($fileout);
}

sub extract_umi_from_bam {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	my @bams;
	my @ubams;
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn  = $self->nproc;                          # if $self->nocluster;
	my $java = $self->project->getSoftware('java');
	$java = "java" unless -e $java;
	my $fgbio = $java . " -jar " . $self->project->getSoftware('fgbio');
	#
	my $nb_bam = 1;
	foreach my $f (@$filein) {
		my $bam = $f;
		$bam =~ s/sam/_umi_extracted\.sam/;
		my @jobs;
		print $self->patient->name() . ":\n";
		my $cmd =
			$fgbio
		  . " ExtractUmisFromBam -i "
		  . $f . " -o "
		  . $bam
		  . " -r 3M3S+T 3M3S+T -t RX -a ";
		my $type     = "extract_umi#" . $nb_bam;
		my $stepname = $self->patient->name . "@" . $type;
		my $job_bds  = job_bds_tracking->new(
			uuid         => $self->bds_uuid,
			cmd          => ["$cmd"],
			name         => $stepname,
			ppn          => $ppn,
			filein       => [$f],
			fileout      => $bam,
			type         => $type,
			dir_bds      => $self->dir_bds,
			software     => "picard",
			sample_name  => $self->patient->name(),
			project_name => $self->patient->getProject->name
		);

		$self->current_sample->add_job( { job => $job_bds } );
		push( @jobs, $job_bds );
		$nb_bam++;
		if ( $self->unforce() && -e $bam ) {
			$job_bds->skip();
		}
		push( @ubams, $f );
		push( @bams,  $bam );
	}

	#my ($ubam_merged) = $self->merge_picard(filein=>\@ubams);
	#my $file = $ubam_merged ;
	#$file =~ s/merge/unmapped.merge\.sam/;
	#system("mv $ubam_merged $file");
	my ($fileout) = $self->merge_picard( filein => \@bams );
	return ($fileout);
}

sub sam_to_fastq {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	my $f1           = $dirout . $name . "_umi_extracted_R1.fastq";
	my $f2           = $dirout . $name . "_umi_extracted_R2.fastq";

	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn  = $self->nproc;                          # if $self->nocluster;
	my $java = $self->project->getSoftware('java');
	$java = "java" unless -e $java;
	my $picard = $java . " -jar " . $self->project->getSoftware('picard_path');
	my $bedtools = $self->project->getSoftware('bedtools');
	my @jobs;
	print $self->patient->name() . ":\n";

	#my $cmd = $bedtools. " bamtofastq -i ". $filein ." -fq ".$f1."  -fq2 ".$f2;
	my $cmd =
		$picard
	  . " SamToFastq I="
	  . $filein
	  . " FASTQ="
	  . $f1
	  . " SECOND_END_FASTQ="
	  . $f2
	  . " CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2";
	my $fileout = $f1;

	my $type     = "sam_to_fastq";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "picard",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $f1 && -e $f2 ) {
		$job_bds->skip();
	}

	#		my @files = ($f1,$f2);
	#		($fileout) = $self->fastp(filein=>\@files);
	return ($fileout);
}

sub fastp {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	my $in1          = $filein;
	my $in2          = $in1;
	$in2 =~ s/_R1/_R2/;
	my $f1         = $dirout . $name . "_umi_extracted_trimmed_R1.fastq";
	my $f2         = $dirout . $name . "_umi_extracted_trimmed_R2.fastq";
	my $fastp_html = $dirout . "/" . $name . "_fastp.html";
	my $fastp_json = $dirout . "/" . $name . "_fastp.json";
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn  = $self->nproc;                          # if $self->nocluster;
	my $java = $self->project->getSoftware('java');
	$java = "java" unless -e $java;

	# my $picard =  $java." -jar ".$self->project->getSoftware('picard_path');
	my $fastp = $self->project->getSoftware('fastp');
	my @jobs;
	print $self->patient->name() . ":\n";
	my $cmd =
		$fastp . " -i "
	  . $in1 . " -o "
	  . $f1 . " -I "
	  . $in2 . " -O "
	  . $f2
	  . " -g -W 5 -q 20 -u 75 -x -3 -g -l 75 -c -j $fastp_json -h $fastp_html -w 12";
	my $fileout = $f1;

	my $type     = "fastp";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "picard",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $f1 && -e $f2 ) {
		$job_bds->skip();
	}

	#my @files = ($f1,$f2);
	#($fileout) = $self->fastp(filein=>\@files);
	return ($fileout);
}

sub merge_ubam {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	my @extracted_umi =
	  glob( $dirout . "/" . $name . "*_unmapped_umi_extracted.sam" );
	my $ppn = $self->nproc;    # if $self->nocluster;
	my $fileout;
	($fileout) = $self->merge_picard( filein => \@extracted_umi );

	return ($fileout);
}

sub filter_proper_pairs {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	my $fileout      = $filein;
	$fileout =~ s/bam/filtered\.bam/;

	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn = $self->nproc;    # if $self->nocluster;
	my $samtools = $self->project->getSoftware('samtools');
	my @jobs;
	print $self->patient->name() . ":\n";
	my $cmd = $samtools . " view -f 2 -bh " . $filein . "> " . $fileout;

	my $type     = "filter_proper_pairs";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => $filein,
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "picard",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}

	#($fileout) = $self->group_reads_by_umi(filein=>\@files);
	return ($fileout);
}

sub bam_to_fastq_umi {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	my $f1           = $dirout . $name . "_R1_L1.fastq";
	my $f2           = $dirout . $name . "_R3_L1.fastq";

	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn  = $self->nproc;                          # if $self->nocluster;
	my $java = $self->project->getSoftware('java');
	$java = "java" unless -e $java;
	my $picard = $java . " -jar " . $self->project->getSoftware('picard_path');
	my $bedtools = $self->project->getSoftware('bedtools');
	my @jobs;
	print $self->patient->name() . ":\n";

	#my $cmd = $bedtools. " bamtofastq -i ". $filein ." -fq ".$f1."  -fq2 ".$f2;
	my $cmd =
		$picard
	  . " SamToFastq I="
	  . $filein . " F="
	  . $f1 . " F2="
	  . $f2
	  . " CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2";
	my $fileout = $f1;

	my $type     = "bam_to_fastq";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "picard",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $f1 && -e $f2 ) {
		$job_bds->skip();
	}
	my @files = ( $f1, $f2 );
	($fileout) = $self->run_alignment_consensus( filein => \@files );
	return ($fileout);
}

sub run_alignment_consensus {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);

	# $dirout .= "/".$self->patient()->name()."/";
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn     = $self->nproc;                    # if $self->nocluster;
	my $fileout = $dirout . $name . ".cons.bam";
	my $nb_bam  = 1;
	my @cmds;
	my @bams;
	my $files;
	my @jobs;
	print $self->patient->name() . ":\n";
	my $file1 = $filein->[0];
	my $file2 = $filein->[1];

	my $bin_dev = $self->script_dir();
	my $f1gz    = $file1 . '.gz';
	my $f2gz    = $file2 . '.gz';
	my $cmd =
qq{gzip -f $file2 && gzip -f $file1 && perl $bin_dev/align.pl -file1=$f1gz -file2=$f2gz -method=$method -lane=$nb_bam  -project=$project_name -name=$name -bam=$fileout -fork=$ppn };
	my $type     = "align_cons#" . $nb_bam;
	my $stepname = $self->patient->name . "@" . $type;

	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [ $file1, $file2 ],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		sample_name  => $name,
		project_name => $project_name,
		software     => $method
	);
	$self->current_sample->add_job( { job => $job_bds } );
	push( @jobs, $job_bds );
	$nb_bam++;

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}

	return ($fileout);
}

#merge_bam_ubam
sub merge_bam_ubam {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	my $fileout      = $filein;
	$fileout =~ s/bam/merge_u_bam\.bam/;
	my $ubam = $dirout . "/" . $self->patient()->name . ".annotate.u.bam";
	$ubam = $dirout . "/" . $self->patient()->name . "_merge.bam"
	  unless -e ($ubam);

#$ubam = $dirout."/".$self->patient()->name."umi_extracted.sam" if -e ($dirout."/".$self->patient()->name.".annotate.u.bam") ;
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn  = $self->nproc;                          # if $self->nocluster;
	my $java = $self->project->getSoftware('java');
	$java = "java" unless -e $java;
	my $picard = $java . " -jar " . $self->project->getSoftware('picard_path');
	my @jobs;
	print $self->patient->name() . ":\n";
	my $reference = $self->project->genomeFasta();
	my $cmd =
		$picard
	  . " MergeBamAlignment CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT  ALIGNER_PROPER_PAIR_FLAGS=false R="
	  . $reference
	  . " ALIGNED="
	  . $filein
	  . " UNMAPPED="
	  . $ubam
	  . " OUTPUT="
	  . $fileout;
	my $type     = "merge_ubam_bam";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [ $filein, $ubam ],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "picard",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}

	return ($fileout);
}

sub group_reads_by_umi {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();

	#my ($dirin) = $self->patient()->getSequencesDirectory();
	my $dirout  = $self->project->getAlignmentPipelineDir($method);
	my $fileout = $filein;
	$fileout =~ s/bam/group\.bam/;
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn  = $self->nproc;                          # if $self->nocluster;
	my $java = $self->project->getSoftware('java');
	$java = "java" unless -e $java;
	my $fgbio = $java . " -jar " . $self->project->getSoftware('fgbio');
	my @jobs;
	my $file_hist = $dirout . "/" . $name . "_umi_group_data.tsv";
	print $name. ":\n";
	my $cmd =
		$fgbio
	  . " GroupReadsByUmi  -s adjacency -o "
	  . $fileout . " -i "
	  . $filein
	  . " -f $file_hist";
	my $type     = "group";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "fgbio",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}

	return ($fileout);
}

sub call_consensus_reads {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();

	#my ($dirin) = $self->patient()->getSequencesDirectory();
	my $dirout  = $self->project->getAlignmentPipelineDir($method);
	my $fileout = $filein;
	$fileout =~ s/bam/consensus\.bam/;
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn  = $self->nproc;                          # if $self->nocluster;
	my $java = $self->project->getSoftware('java');
	$java = "java" unless -e $java;
	my $fgbio = $java . " -jar " . $self->project->getSoftware('fgbio');
	my @jobs;
	print $self->patient->name() . ":\n";
	my $cmd =
		$fgbio
	  . " CallMolecularConsensusReads  -M 1 --read-group-id="
	  . $name . " -o "
	  . $fileout . " -i "
	  . $filein;
	my $type     = "consensus";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "fgbio",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}

	return ($fileout);
}

sub filter_consensus_read {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();

	#my ($dirin) = $self->patient()->getSequencesDirectory();
	my $dirout  = $self->project->getAlignmentPipelineDir($method);
	my $fileout = $filein;
	$fileout =~ s/bam/filter\.bam/;
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn  = $self->nproc;                          # if $self->nocluster;
	my $java = $self->project->getSoftware('java');
	$java = "java" unless -e $java;
	my $fgbio = $java . " -jar " . $self->project->getSoftware('fgbio');
	my $ref   = $self->project->genomeFasta();
	my @jobs;
	print $self->patient->name() . ":\n";
	my $cmd =
		$fgbio
	  . " FilterConsensusReads  -M 1  -N 30  -o "
	  . $fileout . " -i "
	  . $filein . " -r "
	  . $ref;
	my $type     = "filter";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "fgbio",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}

	return ($fileout);
}

sub concat_fastq_umi {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my ($dirin)      = $self->patient()->getSequencesDirectory();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn     = $self->nproc;                           # if $self->nocluster;
	my $fileout = $dirout . "/" . $name . ".umi.fastq.gz";
	my @seq = glob( $dirin . $name . "*_R2_*.fastq.gz" );
	my $in  = join( " ", @seq );
	my @jobs;
	print $self->patient->name() . ":\n";

	my $cmd = "cat " . $in . ">" . $fileout;

	#$cmd = $cmd."  && mv $filebai $fileout.bai ";

	my $type     = "concat";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "zcat",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}

	return ($fileout);
}

sub annotate_with_umi {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my ($dirin)      = $self->patient()->getSequencesDirectory();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	my @seq          = glob( $dirout . $name . "*umi*.fastq.gz" );

	die("several fastq with umi for $name") if scalar(@seq) > 1;
	my $fileout = $filein;
	$fileout =~ s/merge/annotate/;
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn  = $self->nproc;                          # if $self->nocluster;
	my $java = $self->project->getSoftware('java');
	$java = "java" unless -e $java;
	my $fgbio = $java . " -jar " . $self->project->getSoftware('fgbio');
	my @jobs;
	print $self->patient->name() . ":\n";

	my $file_umi = $seq[0];

	#	die("problem $file_umi $dirin :".$name) unless -e $file_umi;
	my $cmd =
		$fgbio
	  . " AnnotateBamWithUmis -i "
	  . $filein . " -f "
	  . $file_umi . " -o "
	  . $fileout;
	my $type     = "annotate";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "fgbio",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub sort_sam {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	my $fileout      = $filein;
	$fileout =~ s/bam/sort\.bam/;
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn  = $self->nproc;                          # if $self->nocluster;
	my $java = $self->project->getSoftware('java');
	$java = "java" unless -e $java;
	my $picard = $java . " -jar " . $self->project->getSoftware('picard_path');
	my @jobs;
	print $self->patient->name() . ":\n";

	my $cmd =
		$picard
	  . " SortSam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT  INPUT="
	  . $filein
	  . " OUTPUT="
	  . $fileout
	  . " SORT_ORDER=coordinate";

	#$cmd = $cmd."  && mv $filebai $fileout.bai ";

	my $type     = "sort_sam";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "picard",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}

	return ($fileout);
}

sub sort_sam_umi {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);
	my @bams;
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn  = $self->nproc;                          # if $self->nocluster;
	my $java = $self->project->getSoftware('java');
	$java = "java" unless -e $java;
	my $picard = $java . " -jar " . $self->project->getSoftware('picard_path');
	my $nb_bam = 1;

	foreach my $f (@$filein) {
		my $bam = $f;
		$bam =~ s/bam/sort\.bam/;
		my @jobs;
		print $self->patient->name() . ":\n";
		my $cmd =
			$picard
		  . " SortSam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT  INPUT="
		  . $f
		  . " OUTPUT="
		  . $bam
		  . " SORT_ORDER=coordinate";
		my $type     = "sort_sam#" . $nb_bam;
		my $stepname = $self->patient->name . "@" . $type;
		my $job_bds  = job_bds_tracking->new(
			uuid         => $self->bds_uuid,
			cmd          => ["$cmd"],
			name         => $stepname,
			ppn          => $ppn,
			filein       => [$f],
			fileout      => $bam,
			type         => $type,
			dir_bds      => $self->dir_bds,
			software     => "picard",
			sample_name  => $self->patient->name(),
			project_name => $self->patient->getProject->name
		);

		$self->current_sample->add_job( { job => $job_bds } );
		push( @jobs, $job_bds );
		$nb_bam++;
		if ( $self->unforce() && -e $bam ) {
			$job_bds->skip();
		}
		push( @bams, $bam );
	}
	my ($fileout) = $self->merge_picard( filein => \@bams );
	return ($fileout);
}

sub merge_bam {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};
	return $self->merge_bamba( {filein => $filein} );
}

sub merge_bamba {
	my ( $self, $hash ) = @_;
	my $filein  = $hash->{filein};
	my $name    = $self->patient->name();
	my $project = $self->patient->getProject;
	my ($f)     = File::Util->new();
	my $m       = $self->patient->alignmentMethod();
	my $fileout =
	  $project->getAlignmentPipelineDir($m) . "/" . $name . ".align.bam";

#$fileout = $project->getAlignmentPipelineDir($m) . "/" . $name . "F[12345678].u.sort.bam" unless -e $fileout;
	my $merge_files = join( " ", @$filein );
	my $bamba       = $self->project->getSoftware('sambamba');

	my $ppn = $self->nproc;
	$ppn = int( $self->nproc / 2 ) if $self->nocluster;

	my $cmd = " $bamba  merge -t $ppn $fileout $merge_files  ";

	if ( scalar(@$filein) == 1 ) {
		if ( -e $fileout ) {
			$cmd = "ln -s  $merge_files $fileout  ";  #if (scalar(@files) == 1);

		}
		else { $cmd = "ln -s  $merge_files $fileout " }

	}
	my $type     = "merge-bamba";
	my $stepname = $self->patient->name . "@" . $type;

	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => $filein,
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "sambamba",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce()
		&& ( -e $fileout or -e $self->patient()->getBamFileName() ) )
	{
		$job_bds->skip();
	}
	return ($fileout);
}

sub read_group_illumina {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $name = $self->patient()->name();
	my $fileout;
	unless ($filein) {
		$filein = $self->patient()->getBamFileName();

		#	 $filein =~s/bwa/bwa\/moitie/;
		die() unless -e $filein;
		my $outputdir = $self->project->getAlignmentPipelineDir("bwa");
		$fileout = $outputdir . "/" . $self->patient()->name . ".rg.bam";
	}
	else {
		$fileout = $filein;
		$fileout =~ s/bam/rg\.bam/;
	}
	my $filebai = $fileout;
	$filebai =~ s/bam$/bai/;

	my $run                = $self->patient()->getRun();
	my $machine            = $run->machine;
	my $run_name           = $run->plateform_run_name();
	my $run_date           = $run->date;
	my $constructor        = $run->machine_constructor();
	my $constructormachine = $constructor . "-" . $machine;
	my $plateform          = $run->plateform();
	my $bar_code           = $self->patient()->barcode();
	my $project            = $self->patient->getProject;
	my $project_name       = $project->name();
	my $patient_name       = $self->patient()->name();

	my $picard_path = $self->project->getSoftware('picard_path');

	my $ppn = 2;
	$ppn = 1 if $self->nocluster;
	my $java = $self->project->getSoftware('java');
	$java = "java" unless -e $java;

	my $picard = $java . " -jar " . $self->project->getSoftware('picard_path');
	my $cmd =
		$picard
	  . " AddOrReplaceReadGroups CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT  INPUT="
	  . $filein
	  . " OUTPUT="
	  . $fileout
	  . " RGDS="
	  . $project_name
	  . " RGLB="
	  . $run_name
	  . " RGSM="
	  . $name
	  . " RGCN="
	  . $plateform
	  . " RGID="
	  . $run_name
	  . " RGPL="
	  . $constructormachine
	  . " RGPU=1 ";
	$cmd = $cmd . "  && mv $filebai $fileout.bai ";

	my $type     = "read-group";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "picard",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}

	return ($fileout);
}

sub move_bam {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $m            = $self->patient->alignmentMethod();
	my $project_name = $self->patient()->getProject->name();

	my $method  = $self->patient()->alignmentMethod();
	my $dirout  = $self->patient()->getProject->getAlignmentDir($method);
	my $fileout = $self->patient()->getBamFileName()
	  ;    #$dirout."/".$self->patient->name.".bam";

	my $version = $self->patient()->project->genome_version();

	my $ppn = 4;
	$ppn = 1 if $self->nocluster;

	die() if $fileout eq $filein;
	my $bin_dev = $self->script_dir;

	my $cmd =
"perl $bin_dev/move_bam.pl -bam=$filein  -project=$project_name -patient=$name -fork=$ppn -version=$version && ln -s $fileout $filein  ";
	if ( $self->again ) {
		my $dir_again = $self->patient->project->getAlignmentDir("start_again");
		my $f         = $dir_again . "/" . $self->patient->name . ".bam";
		die("problem start again") unless -e $f;

		$cmd .= "&& perl $bin_dev/rm_bam.pl $f ";
	}
	my $type     = "move-bam";
	my $stepname = $self->patient->name . "@" . $type;

	# my $job_bds = job_bds_tracking->new();

	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}

	return ($fileout);
}

sub rnaseq_metrics {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	$filein = $self->patient->getBamFile() unless  $filein;
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->project;
	my $m            = $self->patient->alignmentMethod();
	my $project_name = $self->patient()->getProject->name();
	my $dir_out      = $project->getCountingDir("featureCounts") . "/metrics";
	system("mkdir $dir_out && chmod a+rwx $dir_out") unless -e $dir_out;

	my $method  = $self->patient()->alignmentMethod();
	my $fileout = $dir_out . "/$name.metrics";

	my $ppn = 2;
	$ppn = 1 if $self->nocluster;

	die() if $fileout eq $filein;
	$filein = $self->patient->getBamFileName() unless -e $filein;
	my $refFlat = $project->refFlat_file();
	warn $refFlat;
	
	#$refFlat = "/data-isilon/public-data/repository/HG38/refFlat/refFlat_no_chr.txt" if $method eq "star" ;
	
	#$refFlat = $project->refFlat_file_star() if $method eq "star" ;
	#$refFlat = $project->refFlat_file_dragen() if $method eq "dragen-align";
	my $rRNA_file = $project->rRNA_file();
	my $opt       = "";
	$opt = "RIBOSOMAL_INTERVALS=$rRNA_file" if -e $rRNA_file;
	unless ( -e $refFlat ) {
		die("can't find $refFlat $rRNA_file");
	}

	my $java   = $project->buffer->software("java");
	my $picard = $project->buffer->software("picard");

	my $cmd =
"$java -jar $picard  CollectRnaSeqMetrics I=$filein O=$fileout REF_FLAT=$refFlat STRAND=FIRST_READ_TRANSCRIPTION_STRAND $opt  ";
	my $type     = "metrics";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}

	return ($fileout);
}

sub fastqScreen {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->project;
	my $m            = $self->patient->alignmentMethod();
	my $project_name = $self->patient()->getProject->name();
	my ($dirin)      = $self->patient()->getSequencesDirectory();
	my @fastqfiles   = glob("$dirin/*");
	my $dir_out = $project->getCountingDir("featureCounts") . "/fastqScreen";
	system("mkdir $dir_out && chmod a+rwx $dir_out") unless -e $dir_out;
	my $fastqScreen = $project->buffer->software("fastqScreen");
	my $fileout;

	foreach my $f (@fastqfiles) {
		$fileout = $f;
		$fileout =~ s/\.fastq.gz/_screen\.html/;
		my $ppn = 20;
		$ppn = 1 if $self->nocluster;
		my $cmd      = "$fastqScreen $f --outdir $dir_out --threads $ppn";
		my $type     = "fastqScreen";
		my $stepname = $self->patient->name . "@" . $type;
		my $job_bds  = job_bds_tracking->new(
			uuid         => $self->bds_uuid,
			software     => "",
			sample_name  => $self->patient->name(),
			project_name => $self->patient->getProject->name,
			cmd          => [$cmd],
			name         => $stepname,
			ppn          => $ppn,
			filein       => [$filein],
			fileout      => $fileout,
			type         => $type,
			dir_bds      => $self->dir_bds
		);
		$self->current_sample->add_job( { job => $job_bds } );

		if ( $self->unforce() && -e $fileout ) {
			$job_bds->skip();
		}

	}
	return ($fileout);
}

sub bam_sort {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};
	return $self->bam_sort_bamba( filein => $filein );
}

sub bam_sort_bamba {
	my ( $self, $hash ) = @_;
	my $filein  = $hash->{filein};
	my $name    = $self->patient()->name();
	my $fileout = $filein;
	$fileout =~ s/bam/sort\.bam/;

	my $bamba = $self->project->getSoftware('sambamba');

	my $patient_name = $self->patient()->name();
	my $project      = $self->project->name();

	#	my $tmpdir = "--tmpdir=/scratch";

	my $ppn = $self->nproc;
	$ppn = int( $self->nproc / 2 ) if $self->nocluster;
	my $cmd =
		$bamba
	  . " sort   -o=$fileout  --nthreads=$ppn  --out=$fileout "
	  . $filein . " ";

	my $type     = "sort-bam";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}

	return ($fileout);

}

sub rmdup {
	my ( $self, $hash ) = @_;

	my $filein = $hash->{filein};
	return $self->rmdup_bamba( {filein => $filein });

}

sub rmdup_nudup {
	my ( $self, $hash ) = @_;
	my $filein  = $hash->{filein};
	my $name    = $self->patient()->name();
	my $fileout = $filein;
	my $m       = $self->patient->alignmentMethod();
	$fileout =~ s/bam/nudup\.bam\.sorted\.dedup\.bam/;
	my $dirout   = $self->project->getAlignmentPipelineDir($m);
	my $umi_file = $dirout . "/" . $name . ".umi.fastq.gz";
	my $nudup    = "/software/distrib/nudup-master/nudup.py";
	my $python   = "python";
	my $tmpdir   = "";
	$tmpdir = "--tmpdir=/tmp" if $self->host eq "morgan";
	my $ppn = $self->nproc;    # if $self->nocluster;
							   #$tmpdir = "" if $self->nocluster;
	my $cmd      = "$python $nudup -f $umi_file  -s 8 -l 8 $filein";
	my $type     = "rmdup_nudup";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "nudup",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub start_again {
	my ( $self, $hash ) = @_;
	my $filein  = $hash->{filein};
	my $project = $self->patient()->getProject();
	$filein = $self->patient()->getBamFileName();

	my $dir_out = $project->getAlignmentDir("start_again");
	my $fileout = $dir_out . "/" . $self->patient()->name . ".bam";
	if ( -e $fileout && -e $filein ) {
		die("problem file in prod and start_again");
	}
	unless ( -e $fileout ) {
		warn $fileout;
		die() unless -e $filein;
	}

	$self->again(1);
	if ( -e $fileout ) {
		warn "warning already exists $fileout";
		return ($fileout);
	}

	system("mkdir -p $dir_out ") unless -e $dir_out;
	$filein = $self->patient()->getBamFileName();
	my $cmd = "mv $filein* $dir_out/ ";
	system($cmd);
	die unless -e $fileout;
	return ($fileout);

}

sub replace_bam {
	my ( $self, $hash ) = @_;
	my $filein  = $hash->{filein};
	my $project = $self->patient()->getProject();
	my $fileout = $self->patient()->getBamFileName();
	warn( "can't have write access to your bam file " . $filein )
	  unless -w $fileout;
	die() unless -e $fileout;
	my $ppn     = 1;
	my $bin_dev = $self->script_dir;
	my $cmd =
" test -e $filein || exit 1 &&  mv $filein $fileout || exit 1 && mv $filein.bai $fileout.bai; test -e $fileout.bai"
	  ; #"mv $filein $fileout; samtools index $fileout && test -e $fileout.bai";
	my $type     = "replace_bam";
	my $stepname = $self->patient->name . "@" . $type;
	$filein .= ".bai";
	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout . ".out",
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );
	return ($fileout);
}

sub bazam {
	my ( $self, $hash ) = @_;
	my $filein  = $hash->{filein};
	my $project = $self->patient()->getProject();
	unless ($filein) {
		$filein = $self->patient()->getBamFileName();
	}
	my $name         = $self->patient->name;
	my $project_name = $project->name();
	my $dir_out      = $project->getAlignmentPipelineDir("bwa");
	system("mkdir -p $dir_out ") unless -e $dir_out;
	my $fileout = $dir_out . "/" . $name . ".bazam.bam";
	my $bc      = $self->patient()->barcode();
	my $bin_dev = $self->script_dir;

	my $ppn  = 5;
	my $type = "bazam";

	my $cmd =
qq{ $bin_dev/run_cmd_change_bazam.pl -project=$project_name -patient=$name;test -e $fileout};
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "bazam"
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub change_chrMT {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};
	my $name   = $self->patient()->name();
	unless ($filein) {
		$filein = $self->patient()->getBamFileName();
	}
	my $project = $self->patient()->getProject();
	my $dir_out = $project->getAlignmentPipelineDir("bwa");
	system("mkdir -p $dir_out ") unless -e $dir_out;
	my $fileout = $dir_out . "/" . $self->patient->name . ".cng.bam";
	my $bc      = $self->patient()->barcode();
	my $bin_dev = $self->script_dir;

	my $ppn  = 40;
	my $type = "chrMT";

	my $cmd =
qq{ $bin_dev/run_cmd_change_new.bam.pl -fork=20 -dir=$dir_out -file1=$filein -file2=$fileout -bc=$bc -patient=$name;test -e $fileout};
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub elprep5 {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};
	my $name   = $self->patient()->name();

	my $fileout = $filein;
	my $m       = $self->patient->alignmentMethod();
	$fileout =~ s/bam/elprep\.bam/;

#RGDS=".$project_name." RGLB=".$run_name." RGSM=".$name." RGCN=".$plateform." RGID=".$run_name." RGPL=".$constructormachine." RGPU=1 " ;
#@RG	ID:190529_NB501645_0220_AHFYCNBGXB	CN:LAVOISIER	DS:NGS2019_2521	LB:190529_NB501645_0220_AHFYCNBGXB	PL:ILLUMINA-NEXTSEQ500	PU:1	SM:OUT-DEL
	my $project = $self->project;
	my $run     = $self->patient->getRun();
	my $ID      = "ID:" . $run->name;
	my $CN      = "CN:" . $run->plateform;
	my $DS      = "DS:" . $project->name;
	my $PL      = "PL:" . $run->machine_constructor() . "-" . $run->machine;
	my $PU      = "PU:1";
	my $LB      = "LB:" . $run->name;
	my $SM      = "SM:" . $self->patient->name();
	my $rg_string =
		$ID . " "
	  . $CN . " "
	  . $DS . " "
	  . $PL . " "
	  . $PU . " "
	  . $LB . " "
	  . $SM;
	my $bed = $self->patient()->getCaptureBedFile();
	my $tmpdir =
	  "/tmp/";    #$self->patient->project->getAlignmentPipelineDir("elprep");
	$tmpdir .= "/$name";

	unless ( -e $tmpdir ) {
		system("mkdir -p $tmpdir; chmod a+rwx $tmpdir");
	}
	my $elprep = $project->buffer->software("elprep");

	my $samtools = $project->buffer->software("samtools");

	my $ref_root = $project->get_public_data_directory;
	my $ref      = $project->dirGenome() . $project->buffer->index("elprep");
	die($ref) unless -e $ref;
	my $known_sites = $ref_root
	  . "/elprep/dbsnp_137.hg19.elsites,$ref_root/elprep/Mills_and_1000G_gold_standard.indels.hg19.sites.elsite";
	die($known_sites) unless -e $ref_root . "/elprep/dbsnp_137.hg19.elsites";
	die($known_sites)
	  unless -e $ref_root
	  . "/elprep/Mills_and_1000G_gold_standard.indels.hg19.sites.elsite";
	my $recal =
		$self->patient->project->getRecalDir($m) . "/"
	  . $self->patient->name
	  . ".recal.table";
	my $vcfout   = $fileout . ".vcf.gz";
	my $fileout2 = $filein;
	$fileout2 =~ s/bam/elprep_tmp\.bam/;
	my $cmd2     = "filter ";
	my $tmp      = "";
	my $gvcf_arg = "";

	if ( $project->isGenome ) {
		$cmd2 = "sfm ";
		$tmp  = "--tmp-path $tmpdir";
		my $dir_gvcf_out =
		  $self->patient()->project->getGvcfDir("haplotypecaller4");
		my $gvcfout = $dir_gvcf_out . "/" . $name . ".g.vcf.gz";
		$gvcf_arg = " --haplotypecaller $gvcfout ";
		my $ref = $project->dirGenome() . $project->buffer->index("elprep");
		die($ref) unless -e $ref;
		my $known_sites = $ref_root
		  . "/elprep/dbsnp_137.hg19.elsites,$ref_root/elprep/Mills_and_1000G_gold_standard.indels.hg19.sites.elsite";
		die($known_sites)
		  unless -e $ref_root . "/elprep/dbsnp_137.hg19.elsites";
		die($known_sites)
		  unless -e $ref_root
		  . "/elprep/Mills_and_1000G_gold_standard.indels.hg19.sites.elsite";

		$gvcf_arg = " --known-sites $known_sites " . $gvcf_arg;

		#$cmd2 = "sfm "
	}

#--tmp-path /tmp --target-regions $bed --known-sites $known_sites --haplotypecaller $vcfout

	my $cmd =
qq{/home/pnitschk/go/bin/elprep  $cmd2  $filein $fileout $tmp  --replace-read-group "$rg_string" --mark-duplicates  --sorting-order coordinate --reference $ref  $gvcf_arg };
	die();
	my $ppn      = 40;
	my $type     = "elprep";
	my $stepname = $self->patient->name . "@" . $type;
	my $bin_dev  = $self->script_dir;
	my $cmd_verif =
"$samtools index -@ $ppn $fileout && perl $bin_dev/verif_bam.pl -file1=$filein -file2=$fileout ";
	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		cmd          => ["$cmd && $cmd_verif"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "elprep",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce()
		&& ( -e $fileout or -e $self->patient()->getBamFileName() ) )
	{
		$job_bds->skip();
	}
	return ($fileout);
}

sub elprep5_gvcf {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->project;
	my $project_name = $project->name;
	my $dir_gvcf_out =
	  $self->patient()->project->getGvcfDir("haplotypecaller4");
	my $fileout = $dir_gvcf_out . "/" . $name . ".g.vcf.gz";
	my $ppn     = 40;
	my $bin_dev = $self->script_dir;
	my $cmd =
"perl $bin_dev/elprep/elprep5_gvcf.pl /elprep5_gvcf.pl -project=$project_name -patient=$name -bam=$filein -padding=250 -fileout=$fileout";
	my $type      = "elprep_gvcf";
	my $stepname  = $self->patient->name . "@" . $type;
	my $cmd_verif = "";
	my $job_bds   = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		cmd          => ["$cmd"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "elprep",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && ( -e $fileout ) ) {
		$job_bds->skip();
	}
	return ($filein);
}

sub elprep {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};
	my $name   = $self->patient()->name();

	my $fileout = $filein;
	my $m       = $self->patient->alignmentMethod();
	$fileout =~ s/bam/elprep\.bam/;

#RGDS=".$project_name." RGLB=".$run_name." RGSM=".$name." RGCN=".$plateform." RGID=".$run_name." RGPL=".$constructormachine." RGPU=1 " ;
#@RG	ID:190529_NB501645_0220_AHFYCNBGXB	CN:LAVOISIER	DS:NGS2019_2521	LB:190529_NB501645_0220_AHFYCNBGXB	PL:ILLUMINA-NEXTSEQ500	PU:1	SM:OUT-DEL
	my $project = $self->project;
	my $run     = $self->patient->getRun();
	my $ID      = "ID:" . $run->name;
	my $CN      = "CN:" . $run->plateform;
	my $DS      = "DS:" . $project->name;
	my $PL      = "PL:" . $run->machine_constructor() . "-" . $run->machine;
	my $PU      = "PU:1";
	my $LB      = "LB:" . $run->name;
	my $SM      = "SM:" . $self->patient->name();
	my $rg_string =
		$ID . " "
	  . $CN . " "
	  . $DS . " "
	  . $PL . " "
	  . $PU . " "
	  . $LB . " "
	  . $SM;
	my $tmpdir = $self->patient->project->getAlignmentPipelineDir("elprep");
	$tmpdir .= "/$name";

	unless ( -e $tmpdir ) {
		system("mkdir -p $tmpdir; chmod a+rwx $tmpdir");
	}
	my $elprep = $project->buffer->software("elprep");

	my $samtools = $project->buffer->software("samtools");

	my $ref_root = $project->get_public_data_directory;
	my $ref      = $project->dirGenome() . $project->buffer->index("elprep");
	die($ref) unless -e $ref;
	my $known_sites = $ref_root
	  . "/elprep/dbsnp_137.hg19.elsites,$ref_root/elprep/Mills_and_1000G_gold_standard.indels.hg19.sites.elsite";
	die($known_sites) unless -e $ref_root . "/elprep/dbsnp_137.hg19.elsites";
	die($known_sites)
	  unless -e $ref_root
	  . "/elprep/Mills_and_1000G_gold_standard.indels.hg19.sites.elsite";
	my $recal =
		$self->patient->project->getRecalDir($m) . "/"
	  . $self->patient->name
	  . ".recal.table";

	my $cmd =
qq{$elprep  sfm $filein $fileout --tmp-path $tmpdir --replace-read-group "$rg_string" --mark-duplicates  --sorting-order coordinate --bqsr $recal --bqsr-reference $ref --known-sites $known_sites };
	my $ppn      = 40;
	my $type     = "elprep";
	my $stepname = $self->patient->name . "@" . $type;
	my $bin_dev  = $self->script_dir;
	my $cmd_verif =
"$samtools index $fileout && perl $bin_dev/verif_bam.pl -file1=$filein -file2=$fileout ";
	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		cmd          => ["$cmd && $cmd_verif"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "elprep",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce()
		&& ( -e $fileout or -e $self->patient()->getBamFileName() ) )
	{
		$job_bds->skip();
	}
	return ($fileout);
}

sub change_ref_bazam {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};
	my $name   = $self->patient()->name();
	unless ($filein) {
		$filein = $self->patient()->getBamFileName();
	}
	my $project = $self->patient()->getProject();
	my $dir_out = $project->getAlignmentPipelineDir("bwa");
	system("mkdir -p $dir_out ") unless -e $dir_out;
	my $fileout = $dir_out . "/" . $self->patient->name . ".bazam.bam";
	my $bc      = $self->patient()->barcode();
	my $bin_dev = $self->script_dir;

	my $ppn  = 40;
	my $type = "bazam";

	my $cmd =
qq{ $bin_dev/run_cmd_change_bazam.pl -fork=20 -dir=$dir_out -file1=$filein -file2=$fileout -bc=$bc -patient=$name;test -e $fileout};
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub rmdup_bamba {
	my ( $self, $hash ) = @_;
	my $filein  = $hash->{filein};
	my $name    = $self->patient()->name();
	my $fileout = $filein;
	my $m       = $self->patient->alignmentMethod();
	$fileout =~ s/bam/rmdup\.bam/;
	my $bamba  = $self->project->getSoftware('sambamba');
	my $tmpdir = "";
	$tmpdir = "--tmpdir=/tmp" if $self->host eq "morgan";
	my $ppn = $self->nproc;    # if $self->nocluster;
	$ppn = 40;

	#$tmpdir = "" if $self->nocluster;
	my $cmd =
		$bamba
	  . " markdup   $tmpdir  --nthreads=$ppn --overflow-list-size=2000000  "
	  . $filein
	  . " $fileout";
	my $type     = "rmdup";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub SplitNCigarReads {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};
	my $name   = $self->patient()->name();
	unless ($filein) {
		die();
	}
	my $fileout = $filein;
	$fileout =~ s/bam/splincigars\.bam/;

	my $capture   = $self->patient->getCapture();
	my $multiplex = $capture->multiplexFile();
	die("can't find multiplex file") unless -e $multiplex;
	my $reference = $self->reference();
	my $cmd       = $self->gatk()
	  . " -I  $filein -R $reference -T SplitNCigarReads -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -o $fileout ";
	my $ppn = 4;
	$ppn = 1 if $self->nocluster;

	my $type     = "splitncigar";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "gatk",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub mask_primer_start_end {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};
	my $name   = $self->patient()->name();
	unless ($filein) {
		die();
	}
	my $fileout = $filein;
	$fileout =~ s/bam/mask\.bam/;
	my $capture   = $self->patient->getCapture();
	my $multiplex = $capture->multiplexFile();
	die("can't find multiplex file") unless -e $multiplex;
	my $sam = $filein;
	$sam =~ s/bam/mask\.sam/;
	my $bin_dev = $self->script_dir;

	my $ppn = 4;
	$ppn = $self->nproc if $self->nocluster;

	my $cmd =
"$bin_dev/mask_primer_start_end.pl -file=$filein -primer=$multiplex -fork=$ppn | samtools view -Sb - > $fileout  ";

	my $type     = "mask-primer";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub covariate_illumina {
	my ( $self, $hash ) = @_;
	my $filein  = $hash->{filein};
	my $name    = $self->patient()->name();
	my $pname   = $self->project->name;
	my $fileout = $filein;
	$fileout =~ s/bam/recal\.bam/;
	my $csv = $filein . ".tmp1";

	# $fileout =~ s/bam/recal1\.bam/;
	unless ($filein) {
		confess();
	}
	my $filein_bai, my $ppn = $self->nproc;    #if $self->nocluster;

	my $real_ppn = 8;
	$real_ppn = int( $self->nproc / 2 ) if $self->nocluster;

	my $csv2 = $self->patient->getRecalFile();
	my $bed  = $self->patient()->getCaptureBedFile();
	die("not found bed file : $bed") unless -e $bed;
	my $bin_dev = $self->script_dir;
	my $cmd =
"$bin_dev/recalibration.pl -filein=$filein -fileout=$fileout -fork=$real_ppn -project=$pname -patient=$name ";

	my $type     = "recalibration";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "gatk",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub recalibration_table {
	my ( $self, $hash ) = @_;
	my $filein  = $hash->{filein};
	my $name    = $self->patient()->name();
	my $pname   = $self->project->name;
	my $fileout = $filein;
	my $csv     = $filein . ".tmp1";
	$fileout =~ s/bam/recal1\.bam/;
	unless ($filein) {
		confess();
	}
	my $filein_bai, my $ppn = $self->nproc;    #if $self->nocluster;

	my $real_ppn = 8;
	$real_ppn = int( $self->nproc / 2 ) if $self->nocluster;

	my $csv2 = $self->patient->getRecalFile();
	my $bed  = $self->patient()->getCaptureBedFile();
	die("not found bed file : $bed") unless -e $bed;
	my $bin_dev = $self->script_dir;
	my $cmd =
"$bin_dev/recalibration_table.pl -filein=$filein -fileout=$fileout -fork=$real_ppn -project=$pname -patient=$name && test -e $csv2 && mv $filein $fileout ";

	my $type     = "recalibration";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($filein);
}

sub realign_recal {
	my ( $self, $hash ) = @_;
	my $filein   = $hash->{filein};
	my $previous = $hash->{previous};
	my $name     = $self->patient()->name();
	my $pname    = $self->project->name;
	$filein = $self->current_sample->current_bam;
	unless ($filein) {
		confess() unless -e $filein;
	}

	my $outputdir = $self->project->getAlignmentPipelineDir("bwa");
	my $csv       = $outputdir . "/" . $name . ".interval_list";
	my $table     = $outputdir . "/" . $name . ".table";
	my $fileout   = $filein;
	$fileout =~ s/bam/realign\.bam/;

	my $ppn     = $self->nproc;        # if $self->nocluster;
									   #info_step
	my $bin_dev = $self->script_dir;

#my $cmd1 = "$bin_dev/realign_csv.pl -filein=$filein -fileout=$fileout -fork=$ppn -project=$pname -patient=$name ";
#my $cmd2 = "$bin_dev/recal_table.pl -filein=$filein -fileout=$fileout -fork=$ppn -project=$pname -patient=$name ";
	my $cmd =
"$bin_dev/realign_recal.pl -filein=$filein -fileout=$fileout -fork=$ppn -project=$pname -patient=$name ";
	my $type     = "realign-recal";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		uuid         => $self->bds_uuid,
		software     => "gatk",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub breakdancer {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	my $final_dir    = $project->getVariationsDir("breakdancer");
	my $fileout      = $final_dir . "/" . $name . ".breakdancer.txt.gz";
	$filein = $self->patient()->getBamFileName();    # unless $filein;

	my $ppn = $self->nproc;

	$ppn = int( $self->nproc / 2 ) if $self->nocluster;
	die( "-" . $filein ) unless $filein;

	my $bin_dev = $self->script_dir;

	my $cmd =
"perl $bin_dev/breakdancer.pl -project=$project_name  -patient=$name -fork=$ppn";
	my $type     = "breakdancer";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "breakdancer",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub cnvnator {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	my $fileout =
	  $project->getVariationsDir("cnvnator") . "/" . $name . ".vcf.gz";
	$filein = $self->patient()->getBamFileName();    # unless $filein;

	my $ppn = $self->nproc;

	$ppn = int( $self->nproc / 2 ) if $self->nocluster;
	die( "-" . $filein ) unless $filein;

	my $bin_dev = $self->script_dir;
	my $version = $self->patient()->project->genome_version();
	my $cmd =
"perl $bin_dev/cnvnator/cnvnator.pl -project=$project_name  -patient=$name -fork=$ppn -version=$version";
	my $type     = "cnvnator";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "cnvnator",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub canvas {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	my $fileout =
	  $project->getVariationsDir("canvas") . "/" . $name . ".vcf.gz";
	$filein = $self->patient()->getBamFileName();    # unless $filein;

	my $ppn = $self->nproc;

	$ppn = int( $self->nproc / 2 ) if $self->nocluster;
	die( "-" . $filein ) unless $filein;

	my $bin_dev = $self->script_dir;
	my $version = $self->patient()->project->genome_version();
	my $cmd =
"perl $bin_dev/canvas.pl -project=$project_name  -patient=$name -fork=$ppn -version=$version";
	my $type     = "canvas-calling";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "canvas",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub melt {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	my $fileout = $project->getVariationsDir("melt") . "/" . $name . ".vcf.gz";
	$filein = $self->patient()->getBamFileName();    # unless $filein;

	my $ppn = $self->nproc;

	$ppn = int( $self->nproc / 2 ) if $self->nocluster;
	die( "-" . $filein ) unless $filein;

	my $bin_dev = $self->script_dir;
	my $version = $self->patient()->project->genome_version();
	my $cmd =
"perl $bin_dev/melt/melt.pl -project=$project_name  -patient=$name -fork=$ppn -version=$version";
	my $type     = "melt-calling";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "melt",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub lumpy {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	my $fileout = $project->getVariationsDir("lumpy") . "/" . $name . ".vcf.gz";
	$filein = $self->patient()->getBamFileName();    # unless $filein;

	my $ppn = $self->nproc;

	$ppn = int( $self->nproc / 2 ) if $self->nocluster;
	die( "-" . $filein ) unless $filein;

	my $bin_dev = $self->script_dir;

	my $cmd =
"perl $bin_dev/lumpy.pl -project=$project_name  -patient=$name -fork=$ppn";
	my $type     = "lumpy-calling";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "lumpy",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub manta {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	my $fileout = $project->getVariationsDir("manta") . "/" . $name . ".vcf.gz";
	$filein = $self->patient()->getBamFileName();    # unless $filein;

	my $ppn = $self->nproc;

	$ppn = int( $self->nproc / 2 ) if $self->nocluster;
	die( "-" . $filein ) unless $filein;

	my $bin_dev = $self->script_dir;
	my $version = $self->patient()->project->genome_version();
	my $cmd =
"perl $bin_dev/manta.pl -project=$project_name  -patient=$name -fork=$ppn -version=$version";
	my $type     = "manta-calling";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "manta",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub callable_region_panel {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	my $fileout      = $filein;
	my $low_calling  = "";
	$filein = $self->patient()->getBamFileName();    # unless $filein;

	my $dirout = $project->getCallingPipelineDir("callable");
	$fileout = $dirout . "/" . $name . ".freeze";
	my $ppn       = $self->nproc;
	my $reference = $project->genomeFasta();

	$ppn = int( $self->nproc / 2 ) if $self->nocluster;

	#$ppn =5;
	die( "-" . $filein ) unless $filein;

	my $bin_dev = $self->script_dir;

	my $cmd =
"perl $bin_dev/callable_region.pl -project=$project_name  -patient=$name -fork=$ppn  -fileout=$fileout  -filein=$filein";
	my $type     = "callable-region";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub calling_panel {
	my ( $self, $hash ) = @_;
	my $filein      = $hash->{filein};
	my $low_calling = $hash->{low_calling};
	my $methods     = $self->patient()->getCallingMethods();

	foreach my $m (@$methods) {
		next if $m eq "seqnext";
		next if $m eq "SmCounter";
		next if $m eq "haplotypecaller4";
		next if $m eq "dude";
		next if $m eq "casava";
		next if $m eq "melt";

		#next unless $m eq "duplicate_region_calling";
		$self->calling_generic(
			{filein      => $filein,
			method      => $m,
			low_calling => $low_calling}
		);
	}
	return $filein;
}

my $synonym_program = {
	"unifiedgenotyper"         => "gatk",
	"haplotypecaller"          => "gatk",
	"haplotypecaller4"         => "gatk",
	"duplicate_region_calling" => "gatk",
	"freebayes"                => "freebayes",
	"p1_freebayes"             => "freebayes",
	"samtools"                 => "bcftools",
	"eif6_freebayes"           => "eif6_freebayes",
	"mutect2"                  => "mutect2",
	"lofreq"                   => "lofreq",
	"dude"                     => "dude",
};

sub calling_generic {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $method       = $hash->{method};
	my $low_calling  = $hash->{low_calling};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	$filein = $self->patient()->getBamFileName() unless $filein =~ /bam/;
	my $dirout  = $project->getVariationsDir($method);
	my $fileout = $dirout . "/$name.vcf.gz";
	my $dirout1 = $project->getCallingPipelineDir("callable");

	#$filein = $dirout1 ."/".$name.".freeze";
	my $ppn = $self->nproc;
	$ppn = int( $self->nproc / 2 ) if $self->nocluster;
	$ppn = 40;
	die( "-" . $filein ) unless $filein;

	#my $low_calling_string = $low_calling;
	#	 $low_calling_string = "-low_calling=1"  if $low_calling;

	my $m       = $self->patient()->alignmentMethod();
	my $dir_bam = $project->getAlignmentDir($m);
	my $bin_dev = $self->script_dir;

#my $cmd = "perl $bin_dev/calling_panel.pl -project=$project_name  -patient=$name -fork=$ppn  -fileout=$fileout -method=$method -filein=$filein $low_calling_string";
	my $cmd =
"perl $bin_dev/calling_panel.pl -project=$project_name  -patient=$name -fork=$ppn  -fileout=$fileout -method=$method -filein=$filein ";
	my $type = "calling-" . $method;

	#	$type = "lc-".$type  if $low_calling;
	my $stepname = $self->patient->name . "@" . $type;
	warn $method unless exists $synonym_program->{$method};
	return       unless exists $synonym_program->{$method};
	die($method) unless exists $synonym_program->{$method};
	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => $synonym_program->{$method},
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub calling_gvcf4 {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	my $low_calling  = "";
	$filein = $self->patient()->getBamFileName();    #unless $filein !~/bam/;

	my $dir_gvcf_out = $project->getGvcfDir("haplotypecaller4");
	my $dir_prod     = $project->getGvcfDir("haplotypecaller4");

	#my $fileout = $dir_prod."/".$name.".end";
	my $fileout = $dir_gvcf_out . "/" . $name . ".g.vcf.gz";
	print "$name \n " unless -e $dir_gvcf_out . "/" . $name . ".g.vcf.gz";

	#probleme du fichier de log non défini dans le pipeline lancé sans pbs
	#	my $cmd = "" ;
	my $ppn = $self->nproc;    # if $self->nocluster;
	$ppn = 40;
	my $real_ppn = $ppn;       #int($self->nproc / 2);
	$real_ppn = 20 if $self->host eq "morgan";
	die( "-" . $filein ) unless $filein;

	#	die($filein. " is empty") if (-z $filein);
	my $bin_dev = $self->script_dir;
	my $version = $self->patient()->project->genome_version();
	my $cmd =
"perl $bin_dev/gatk-4/calling_individual_gvcf.pl -version=$version -project=$project_name  -patient=$name -fork=$real_ppn -out=$fileout -window=5_000_000  ";
	my $type     = "gvcf4";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "gatk4",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub calling_gvcf {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	my $low_calling  = "";
	$filein = $self->patient()->getBamFileName();    #unless $filein !~/bam/;

	my $dir_gvcf_out = $project->getGvcfDir("haplotypecaller");
	my $fileout      = $dir_gvcf_out . "/" . $name . ".g.vcf.gz";
	print "$name \n " unless -e $dir_gvcf_out . "/" . $name . ".g.vcf.gz";

	#probleme du fichier de log non défini dans le pipeline lancé sans pbs
	#	my $cmd = "" ;
	my $ppn      = $self->nproc;              # if $self->nocluster;
	my $real_ppn = int( $self->nproc / 2 );
	$real_ppn = $ppn if $self->host eq "morgan";

	die( "-" . $filein ) unless $filein;

	#	die($filein. " is empty") if (-z $filein);
	my $bin_dev = $self->script_dir;

	my $cmd =
"perl $bin_dev/calling_individual_gvcf.pl -project=$project_name  -patient=$name -fork=$real_ppn -out=$fileout -window=1_000_000  ";
	my $type     = "gvcf";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "gatk",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub callable_regions {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	unless ($filein) {
		my $dir_gvcf_out = $project->getGvcfDir("haplotypecaller");
		$filein = $dir_gvcf_out . "/" . $name . ".g.vcf.gz";
	}
	die() unless $filein;

	my $ppn = 1;
	die( "-" . $filein ) unless $filein;

	#	die($filein. " is empty") if (-z $filein);
	my $bin_dev          = $self->script_dir;
	my $dir_tmp_callable = $project->getCallingPipelineDir("callable");
	my $lmdb_file        = $name . ".ok.callable";
	my $fileout          = $project->getCoverageCallable() . "/" . $lmdb_file;
	my $cmd =
"perl $bin_dev/callable_regions_patients.pl -project=$project_name  -patient=$name   ";
	my $type     = "callable";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub move_vcf {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();

	#my $fileout = $filein;
	my $ppn    = 1;
	my $method = "unifiedgenotyper";
	my $dirout = $project->getCallingPipelineDir($method);
	$filein = $dirout . "/$name/" . $name . ".final.vcf" unless $filein;
	my $dirin   = $project->getCallingPipelineDir($method);
	my $dir_snp = $project->getVariationsDir($method);
	my $fileout = $dir_snp . "/" . $project->name . ".vcf.gz";
	warn $fileout;
	die();
	my $bin_dev = $self->script_dir;
	my $cmd =
"perl $bin_dev/move_individual_vcf.pl -project=$project_name  -patient=$name -method=unifiedgenotyper ";
	my $type     = "move-vcf";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub picard_stats {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	my $fileout      = $self->patient->getMetricsFile();

	my $bin_dev = $self->script_dir;
	my $cmd =
	  "perl $bin_dev/picard_stats.pl -project=$project_name  -patient=$name ";
	my $type     = "stats";
	my $stepname = $self->patient->name . "@" . $type;

	my $ppn = $self->nproc;
	$ppn = int( $self->nproc / 2 ) if $self->nocluster;
	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub lmdb_depth {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	$filein = $self->patient()->getBamFileName();
	my $fileout = $self->patient()->fileNoSqlDepth;
	my $ppn     = $self->nproc;
	$ppn = int( $self->nproc / 2 ) if $self->nocluster;
	### Récupération de la version de samtools

	my $bin_dev = $self->script_dir;

	my $cmd =
qq{perl $bin_dev/coverage_genome.pl -patient=$name  -fork=$ppn  -project=$project_name  };
	if ( $project->isGenome ) {
		$cmd .=
qq{ && perl $bin_dev/coverage_statistics_genome.pl -patient=$name  -fork=$ppn  -project=$project_name};
	}
	my $type     = "lmdb_depth";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}

	return ($fileout);

}

sub transcripts_coverage {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();

	my $coverage_dir = $project->getRootDir() . "/align/coverage/depth/";

	$filein = $coverage_dir . "/" . $name . ".depth";

	my $no      = $self->patient()->getTranscriptsCoverageDepth();
	my $fileout = $no->filename();
	my $ppn     = $self->nproc;
	$ppn = int( $self->nproc / 2 ) if $self->nocluster;

	### Récupération de la version de samtools

	my $bin_dev = $self->script_dir;

	my $cmd =
qq{perl $bin_dev/transcripts/transcripts_coverage.pl -patient=$name  -fork=$ppn  -project=$project_name };

	my $type     = "transcripts_depth";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		,
		cmd     => [$cmd],
		name    => $stepname,
		ppn     => $ppn,
		filein  => [$filein],
		fileout => $fileout,
		type    => $type,
		dir_bds => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);

}

sub genes_dude {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	unless ($filein) {
		my $no = $self->patient()->getTranscriptsDude();
		$filein = $no->filename();
	}

	my $no      = $self->patient()->getGenesDude();
	my $fileout = $no->filename();
	my $ppn     = $self->nproc;
	$ppn = int( $self->nproc / 2 ) if $self->nocluster;

	### Récupération de la version de samtools

	my $bin_dev = $self->script_dir;

	my $cmd =
qq{perl $bin_dev/transcripts/genes_level_dude.pl -patient=$name  -fork=$ppn  -project=$project_name };
	my $type     = "genes_dude";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && $no->nb_keys > 10 ) {
		$job_bds->skip();
	}
	$no->close();
	return ($fileout);
}

sub transcripts_dude {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();

	my $dir_out = $project->getVariationsDir("dude");
	$filein = $dir_out . "/" . $name . ".dude.lid.gz";

	my $no      = $self->patient()->getTranscriptsDude();
	my $fileout = $no->filename();
	my $ppn     = $self->nproc;
	$ppn = int( $self->nproc / 2 ) if $self->nocluster;

	### Récupération de la version de samtools

	my $bin_dev = $self->script_dir;

	my $cmd =
qq{perl $bin_dev/transcripts/transcripts_dude.pl -patient=$name  -fork=$ppn  -project=$project_name };
	my $type     = "transcripts_dude";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && $no->nb_keys > 10 ) {
		$job_bds->skip();
	}
	$no->close();
	return ($fileout);
}

sub wisecondor {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	$filein = $self->patient()->getBamFileName();

	my $fileout = $self->patient()->fileWiseCondor;

	my $ppn = 10;
	$ppn = 2 if $self->nocluster;

	### Récupération de la version de samtools

	my $bin_dev = $self->script_dir;

	my $cmd =
qq{perl $bin_dev/wisecondor.pl -patient=$name    -project=$project_name -fork=$ppn};

	my $type     = "wisecondor";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "wisecondor",
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "wisecondor",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);

}

sub calling_wisecondor {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	unless ($filein) {
		$filein = $self->patient()->fileWiseCondor;
		die() unless -e $filein;
	}

	my $fileout =
		$project->getVariationsDir("wisecondor") . "/"
	  . $self->patient->name
	  . "_aberrations.bed.gz";
	my $ppn = 1;
	$ppn = 1 if $self->nocluster;

	### Récupération de la version de samtools

	my $bin_dev = $self->script_dir;

	my $cmd =
qq{perl $bin_dev/calling_wisecondor.pl -patient=$name    -project=$project_name };

	my $type     = "callingWise";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "wisecondor",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "wisecondor",
		sample_name  => $self->patient->name(),
		project      => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);

}

sub reorder_picard {
	my ( $self, $hash ) = @_;
	my $filein   = $hash->{filein};
	my $previous = $hash->{previous};
	my $ppn      = $self->nproc;
	my $name     = $self->patient()->name();
	my $fileout;
	my $project  = $self->patient()->project;
	my $m        = $self->patient->alignmentMethod();
	my $dir_prod = $project->getAlignmentDir($m);
	unless ($filein) {

		#redemarrage sans fichier bam dans dossier pipeline
		my $outputdir      = $self->project->getAlignmentPipelineDir("bwa");
		my $outputdirstart = $outputdir . "/start_again";

		my $fileprod = $self->patient()->getBamFileName();
		$filein = $outputdirstart . "\/" . $name . ".bam";
		die($fileprod) unless -e $fileprod;
		if ( -e $fileprod ) {

			#              $filein = $fileprod;
			#              $filein =~ s/bwa/bwa\/start_again/;
			$filein = $outputdirstart . "\/" . $name . ".bam";
			system("mkdir $outputdirstart && chmod a+rwx $outputdirstart ")
			  unless -e $outputdirstart;
			system(" chmod a+w $fileprod ");
			system("mv $fileprod $filein");
		}
		confess() unless -e $filein;

		$filein = $outputdirstart . "/" . $self->patient()->name . ".bam";
		$fileout =
			$self->project->getAlignmentPipelineDir("bwa") . "/"
		  . $self->patient()->name()
		  . ".reorder.bam";

	}
	else {
		$fileout = $filein;
		$fileout =~ s/bam/reord\.bam/;

	}
	my $filebai = $fileout;
	$filebai =~ s/bam/bai/;

	my $reference = $project->genomeFasta();

	my $picard_path = $self->project->getSoftware('picard_path');
	my $picard =
		$self->project->getSoftware('java')
	  . " -jar "
	  . $self->project->getSoftware('picard_path');
	my $dict = $reference;
	$dict =~ s/\.fa/\.dict/;
	my $cmd =
		$picard
	  . " ReorderSam   -CREATE_INDEX true -INPUT "
	  . $filein
	  . " -OUTPUT "
	  . $fileout
	  . "  -REFERENCE_SEQUENCE "
	  . $reference
	  . " -ALLOW_INCOMPLETE_DICT_CONCORDANCE true -ALLOW_CONTIG_LENGTH_DISCORDANCE true -SEQUENCE_DICTIONARY $dict";

	my $type     = "reorder_picard";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);

}

##methode remplacée par le script coverage.pl

sub coverage_samtools {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	$filein = $self->patient()->getBamFileName();

	my $coverage_dir = $project->getRootDir() . "/align/coverage/";
	mkdir $coverage_dir unless -e $coverage_dir;
	my $fileout = $coverage_dir . "/" . $name . ".cov.gz";
	my $bed     = $self->patient()->getCaptureFile();
	my $ppn     = $self->nproc;
	$ppn = int( $self->nproc / 2 ) if $self->nocluster;

	### Récupération de la version de samtools

	my $bin_dev = $self->script_dir;

	my $cmd =
qq{perl $bin_dev/coverage.pl -patient=$name -filein=$filein -dir=$coverage_dir -bed=$bed -fork=$ppn -name=$name -project=$project_name };
	my $type     = "coverage";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);

}

sub depthofcoverage {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	$filein = $self->patient()->getBamFileName();

	my $coverage_dir = $project->getRootDir() . "/align/coverage/depth/";
	mkdir $coverage_dir unless -e $coverage_dir;
	my $fileout = $coverage_dir . "/" . $name . ".depth";
	my $bed     = $self->patient()->getCaptureFile();
	my $ppn     = $self->nproc;
	$ppn = int( $self->nproc / 2 ) if $self->nocluster;

	#$ppn =;

	### Récupération de la version de samtools

	my $bin_dev = $self->script_dir;

	my $cmd =
qq{perl $bin_dev/depthofcoverage.pl -patient=$name -filein=$filein -fileout=$fileout -fork=$ppn  -project=$project_name };

	my $type     = "depthofcoverage";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);

}
############################
## Methodes pour les UMIs
############################

sub generate_ubam_umi {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $bin_dev      = $self->script_dir;
	my $method       = $self->patient->alignmentMethod();
	my $project_name = $project->name();
	my $dir     = $project->getAlignmentPipelineDir( $method . "/" . $name );
	my $fileout = "$dir/files.json";
	my $ppn     = 20;
	my $cmd =
qq{perl $bin_dev/umi/split_fast_ubam_dragen.pl -patient=$name -fork=$ppn  -project=$project_name  };
	warn "$cmd";

	#my $fileout = `$cmd -out=1`;
	#chomp($fileout);
	#die() unless $fileout;
	my $type     = "ubam";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub align_bam_combine_ubam_umi {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $name         = $self->patient()->name();
	my $patient      = $self->patient();
	my $project      = $self->patient()->getProject();
	my $bin_dev      = $self->script_dir;
	my $project_name = $project->name();
	my $ppn          = 40;
	my $cmd =
qq{perl $bin_dev/umi/align_combine.pl -patient=$name -fork=$ppn  -project=$project_name -filein=$filein };
	warn $cmd;
	my $method = $patient->alignmentMethod();
	my $dir =
	  $project->getAlignmentPipelineDir( $method . "/" . $patient->name );
	my $fileout = "$dir/list_file.txt";

	my $type     = "ubambam";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub merge_split_bam_umi {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $name         = $self->patient()->name();
	my $patient      = $self->patient();
	my $project      = $self->patient()->getProject();
	my $bin_dev      = $self->script_dir;
	my $project_name = $project->name();
	my $ppn          = $self->nproc;
	my $method       = $patient->alignmentMethod();
	my $dir =
	  $project->getAlignmentPipelineDir( $method . "/" . $patient->name );
	my $fileout = $dir . "/" . $patient->name . ".consensus.umi.bam";
	my $cmd =
qq{perl $bin_dev/umi/merge_bam.pl -patient=$name -fork=$ppn  -project=$project_name -filein=$filein -bamout=$fileout};
	warn $cmd;

	my $type     = "merge_umi";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	#$job_bds->skip();
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub consensus_bam_umi {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $bin_dev      = $self->script_dir;
	my $project_name = $project->name();

	my $ppn = $self->nproc;
	my $cmd =
qq{perl $bin_dev/umi/generateFastqFromConsensusBam.pl -patient=$name -fork=$ppn  -project=$project_name -bamin=$filein  };
	my $patient = $self->patient();
	my $method  = $patient->alignmentMethod();
	my $dir =
	  $project->getAlignmentPipelineDir( $method . "/" . $patient->name );
	my $fileout = "$dir/list_combined.txt";
	$cmd =
qq{perl $bin_dev/umi/generateFastqFromConsensusBam.pl -patient=$name -fork=$ppn  -project=$project_name -bamin=$filein -bamout=$fileout };
	warn $cmd;
	my $type     = "consensus_bam";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {

		$job_bds->skip();
	}
	return ($fileout);
}

sub merge_final_bam {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $bin_dev      = $self->script_dir;
	my $project_name = $project->name();
	my $ppn          = $self->nproc;
	my $method       = $self->patient()->alignmentMethod();
	my $fileout =
		$project->getAlignmentPipelineDir($method) . "/"
	  . $name
	  . ".consensus.umi.bam";
	my $cmd =
qq{perl $bin_dev/umi/merge_bam.pl -patient=$name -fork=$ppn  -project=$project_name -filein=$filein -bamout=$fileout};

	$cmd =
qq{perl $bin_dev/umi/merge_bam.pl -patient=$name -fork=$ppn  -project=$project_name -filein=$filein -bamout=$fileout};
	warn $cmd;
	my $type     = "merge_final_umi";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

####Méthodes générales de manipulation des bam

sub bam_to_fastq {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();

	my $ppn     = 8;
	my $m       = $self->patient()->alignmentMethod();
	my $dir_bam = $project->getAlignmentDir($m);         # où trouver les bams
	my $dir_seq =
	  $self->patient()->getSequencesDirectory();         # où ecrire les fastqs

	my $picard =
		$self->project->getSoftware('java')
	  . " -jar "
	  . $self->project->getSoftware('picard_path');

	$filein = $self->patient()->getBamFileName();
	my $fastq1 = $dir_seq . "/" . $name . "_R1_L001.fastq";
	my $fastq2 = $dir_seq . "/" . $name . "_R2_L001.fastq";

	my $fileout  = $fastq1;
	my $bedtools = $self->project->getSoftware('bedtools');

#my $cmd = $bedtools. " bamtofastq -i ". $filein ." -fq ".$fastq1."  -fq2 ".$fastq2;

	my $cmd =
		$picard
	  . " SamToFastq I="
	  . $filein
	  . " FASTQ="
	  . $fastq1
	  . " SECOND_END_FASTQ="
	  . $fastq2;

	my $type     = "bam-to-fastq";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}

	#my ($fileout) = $self->run_alignment_consensus(filein=>[$fastq1,$fastq2]);
	($fileout) = $self->zip_fastq( { filein => [ $fastq1, $fastq2 ], R => 0 } );

	#	my ($fileout1) = $sef->zip_fastq(filein=>[$fastq1,$fastq2],R=>1);
	#my ($fileout) = $fastq1;
	return ($fileout);

}

sub zip_fastq {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};
	my $R      = $hash->{R};
	my $name   = $self->patient()->name();

	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();

	my $ppn     = 8;
	my $fileout = $filein->[$R] . ".gz";

	my $cmd = "bgzip " . $filein->[$R];

	my $type     = "zip-fastq#" . $R;
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => $filein,
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub cellranger {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $name = $self->patient()->name();
	my $run  = $self->patient->getRun();

	my $project = $self->patient()->getProject();

	my $dirout = $self->project->getAlignmentPipelineDir("cellranger");

	my $ppn = 40;

	#$ppn = int($self->nproc/2) if $self->nocluster;
	#$ppn =;

	### Récupération de la version de samtools

	my $bin_dev    = $self->script_dir;
	my $cellranger = $self->patient->buffer->getSoftware("cellranger");

	my $dir_fastq = $self->patient()->getSequencesDirectory();

	#/public-data/cellranger/MM38-GFP-TOMATO/index
	#--chemistry=SC3Pv2
	my $dir_index = $project->getCellRangerIndex();
	my $cmd =
qq{cd $dirout && $cellranger  count --transcriptome=$dir_index --id=$name --fastqs=$dir_fastq --sample=$name --chemistry=SC3Pv3};
	my $fileout  = $dirout . "$name/outs/cloupe.cloupe";
	my $type     = "cellranger";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);

}

######################
# former bds_calling_steps
######################
has 'method_calling' => (
	is      => 'rw',
	lazy    => 1,
	default => sub {
		return "haplotypecaller4";
	}
);

sub test_genotype {
	my ( $self, $hash ) = @_;
	my $filein      = $hash->{filein};
	my $project     = $self->project;
	my $name        = $project->name();
	my $arg         = $self->argument_patient();
	my $dir_out_vcf = $project->getCallingPipelineDir("haplotypecaller4");
	my $fileout     = $dir_out_vcf . "/" . $name . ".hc.vcf";
	my $type        = "genotype";
	my $stepname    = $name . "@" . $type;
	my $fileout2    = $dir_out_vcf . "/" . $name . ".hc.vcf";

	#my $list_patients = $project->get_list_patients($arg);
	my $patients = $project->get_list_patients($arg);
	my $filesin;
	foreach my $patient (@$patients) {
		my $gvcf = $patient->gvcfFileName("haplotypecaller4")
		  ;    #$dir_out."/".$patient->name.".g.vcf.gz";
		push( @$filesin, $gvcf );
	}
	my $bin_dev = $self->script_dir();
	my $nb      = 1;
	my $files;
	foreach my $chr ( @{ $project->getChromosomes } ) {

		#my $chr = $project->getChromosome(1);
		next if $chr->name eq "MT";

		my $dir_vcf_out = $project->getCallingPipelineDir("haplotypecaller4");
		my @p           = sort { $a cmp $b } split( ",", $arg );
		$dir_vcf_out .= "/" . join( "@", @p );
		system("mkdir $dir_vcf_out ; chmod a+rwx $dir_vcf_out")
		  unless -e $dir_vcf_out;

		#$dir_vcf_out.= "/".$dir_name;
		my $windows = $chr->getWindowCaptureForCalling( 250, 1_000_000 );
		foreach my $w (@$windows) {
			my $as_string = $w->{intspan}->as_string();
			my $digest    = md5_hex($as_string);
			my $vcf =
				$dir_vcf_out . "/"
			  . $chr->name . "."
			  . $w->{start} . "."
			  . $w->{end}
			  . ".$digest.vcf";
			my $project_name = $project->name();
			my $ppn          = 1;
			my $cmd =
"perl $bin_dev/gatk-4//chunk_genotype_gvcf.pl -project=$project_name -patient=$arg -vcf=$vcf -intspan=\"$as_string\" -fork=1 -chr="
			  . $chr->name;
			my $type     = "geno#" . $nb;
			my $stepname = $project->name . "@" . $type;
			my $job_bds  = job_bds_tracking->new(
				uuid         => $self->bds_uuid,
				software     => "",
				sample_name  => $self->patient->name(),
				project_name => $self->patient->getProject->name,
				cmd          => [$cmd],
				name         => $stepname,
				ppn          => $ppn,
				filein       => $filesin,
				fileout      => $vcf,
				type         => $type,
				dir_bds      => $self->dir_bds
			);

#my $job_bds = job_bds_tracking->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$vcf,type=>$type,dir_bds=>$self->dir_bds,sample_name=>$name,project_name=>$project_name,software=>$method);
			$self->current_sample->add_job( { job => $job_bds } );
			if ( $self->unforce() && ( -e $vcf or -e $fileout2 ) ) {
				$job_bds->skip();
			}
			push( @$files, $vcf );
			$nb++;
		}
	}    #end chromosomes

	($fileout) = $self->merge_vcfs( filein => $files );
	return ($fileout);
}

sub merge_vcfs {
	my ( $self, $hash ) = @_;
	my $filein      = $hash->{filein};
	my $project     = $self->project;
	my $name        = $project->name();
	my $arg         = $self->argument_patient();
	my $dir_out_vcf = $project->getCallingPipelineDir("haplotypecaller4");
	my $fileout     = $dir_out_vcf . "/" . $name . ".hc.vcf";
	my $type        = "merge";
	my $stepname    = $name . "@" . $type;
	my $patients    = $project->get_list_patients($arg);

#	foreach my $patient (@$patients){
#		my $gvcf = $patient->getGvcfFile("haplotypecaller4");#$dir_out."/".$patient->name.".g.vcf.gz";
#		die("$gvcf you don't have gvcf for at least this patient restart after_lifescope on this project ") unless $gvcf ;#-e $gvcf;
#}
	my $fork        = $self->nproc;
	my $size_window = 5_000_000;
	my $bin_dev     = $self->script_dir();
	my $cmd =
		"perl $bin_dev/gatk-4/./fast_concat_sort_vcf.pl -project="
	  . $name
	  . " -vcf=$fileout -patient=$arg  ";
	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => 1,
		filein       => $filein,
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);

# my $job_bds = job_bds_tracking->new(uuid=>$self->bds_uuid,cmd=>["$cmd"],name=>$stepname,ppn=>$fork,filein=>$filein,fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds,software=>"gatk4",sample_name=>$arg,project_name=>$self->patient->getProject->name);

	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub genotype_gvcf4 {
	my ( $self, $hash ) = @_;
	my $filein      = $hash->{filein};
	my $project     = $self->project;
	my $name        = $project->name();
	my $arg         = $self->argument_patient();
	my $dir_out_vcf = $project->getCallingPipelineDir("haplotypecaller4");
	my $fileout     = $dir_out_vcf . "/" . $name . ".hc.vcf";
	my $type        = "genotype";
	my $stepname    = $name . "@" . $type;
	my $patients    = $project->get_list_patients($arg);

#	foreach my $patient (@$patients){
#		my $gvcf = $patient->getGvcfFile("haplotypecaller4");#$dir_out."/".$patient->name.".g.vcf.gz";
#		die("$gvcf you don't have gvcf for at least this patient restart after_lifescope on this project ") unless $gvcf ;#-e $gvcf;
#}
	my $fork        = $self->nproc;
	my $size_window = 5_000_000;
	my $bin_dev     = $self->script_dir();
	my $cmd =
		"perl $bin_dev/gatk-4/genotype_gvcf.pl -project="
	  . $name
	  . " -vcf=$fileout -window=$size_window -patient=$arg -fork=$fork ";

#my $job_bds = job_bds_tracking->new(uuid=>$self->bds_uuid,software=>"",sample_name=>$self->patient->name(),project_name=>$self->patient->getProject->name,cmd=>["$cmd"],name=>$stepname,ppn=>$fork,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds,software=>"gatk4",sample_name=>$arg,project_name=>$self->patient->getProject->name);
	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->argument_patient(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $fork,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);

	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub correct_vcf {
	my ( $self, $hash ) = @_;
	my $filein     = $hash->{filein};
	my $project    = $self->project;
	my $name       = $project->name();
	my $dir_out    = $project->getCallingPipelineDir( $self->method_calling );
	my $snp_out    = $dir_out . "/" . $project->name . ".snps.vcf";
	my $indels_out = $dir_out . "/" . $project->name . ".indels.vcf";

	my $type = "correct";

	#my $fileout = $filein;
	unless ($filein) {
		$filein = $dir_out . "/$name.uni.vcf";
		die("can't find $filein ") unless -e $filein;
	}
	my $bin_dev = $self->script_dir();
	my $fork    = 1;
	my $patient = $self->argument_patient();
	my $fileout = "$dir_out/$name.final.vcf";

	#@@@@@@@@@@@@@@@@@@@@@@@
	#todo
	#check if patient is in filein
	#@@@@@@@@@@@@@@@@@@@@@@
	my $stepname = $name . "@" . $type;
	my $cmd =
"perl $bin_dev/correct_gatk.pl -vcf=$filein -project=$name -dir=$dir_out -snp_out=$snp_out -indels_out=$indels_out -patient=$patient -method=haplotypecaller4";
	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->argument_patient(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $fork,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();

		#  	$self->add_skip_steps($stepname);
		#	return ($previous,$snp_out);
	}
	return ($fileout);
}

sub move_vcf_hc4 {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};
	my $name   = $self->project->name();
	my $dir_out =
	  $self->project->getCallingPipelineDir( $self->method_calling );
	my $arg      = $self->argument_patient();
	my $type     = "move-vcf";
	my $stepname = $name . "@" . $type;
	my $fileout =
		$self->project->getVariationsDir( $self->method_calling ) . "/"
	  . $name
	  . ".vcf.gz";
	my $bin_dev = $self->script_dir();

	my $cmd =
"perl $bin_dev/move_vcf.pl  -project=$name -vcf_dir=$dir_out -patient=$arg -method_calling=haplotypecaller4";
	my $ppn     = 1;
	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->argument_patient(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	#@@@@@@@@@@@@@@@@@@@@@
	#todo check for all individual file
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub move_and_split_vcf4 {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};
	my $name   = $self->project->name();
	my $dir_out =
	  $self->project->getCallingPipelineDir( $self->method_calling );
	unless ($filein) {
		$filein = $dir_out . "/" . $self->project->name() . ".final.vcf";
		die($filein) unless -e $filein;
	}
	my $bin_dev  = $self->script_dir();
	my $arg      = $self->argument_patient();
	my $type     = "move_split-vcf";
	my $stepname = $name . "@" . $type;
	my $fileout =
		$self->project->getVariationsDir( $self->method_calling ) . "/"
	  . $name
	  . ".vcf.gz";
	my $cmd =
"perl $bin_dev/split_haplotypecaller.pl  -project=$name -vcf=$filein -patient=$arg -method_calling=haplotypecaller4";
	my $ppn     = 1;
	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->argument_patient(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	#@@@@@@@@@@@@@@@@@@@@@
	#todo check for all individual file
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

####
# By project
####

sub dude {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $projectName = $self->project->name();
	my $logPlink =
	  $self->project->getProjectPath() . '../' . $projectName . '.plink.resume';
	my $ppn = $self->nproc;
	my $tfilein;
	my $coverage_dir = $self->project->getRootDir() . "/align/coverage/depth/";
	foreach my $patient ( @{ $self->project->getPatients } ) {
		my $file = $patient->fileNoSqlDepth;

		#	die("file coverage not found : $file") unless -e $file;
		push( @$tfilein, $file );
	}
	my $bin_dev   = $self->script_dir();
	my $type      = "dude";
	my $stepname  = $projectName . "@" . $type;
	my $dir       = $self->project->project_log();
	my $final_dir = $self->project->getVariationsDir("dude");
	my $fileout   = $final_dir . "/" . $self->project->name . ".dude";

	my $cmd =
"perl $bin_dev/dude/dude.pl -project=$projectName  -fork=$ppn && date > $fileout";
	my $job_bds = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $self->argument_patient,
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => \@$tfilein,
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	#@@@@@@@@@@@@@@@@@@@@@
	#todo check for all individual file
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($filein);

}

sub htlv1_insertion {
	my ( $self, $hash ) = @_;
	my $filein = $hash->{filein};

	my $method       = $self->patient()->alignmentMethod();
	my $name         = $self->patient()->name();
	my ($dirin)      = $self->patient()->getSequencesDirectory();
	my $project_name = $self->project->name();
	my $dirout       = $self->project->getAlignmentPipelineDir($method);

	#$dirout .= "/".$self->patient()->name()."/";
	system("mkdir -p $dirout;chmod a+rwx $dirout") unless -e $dirout;
	my $ppn = $self->nproc;    # if $self->nocluster;
	my $insertion = $self->project->getSoftware('insertion');
	my $bin_dev   = $self->script_dir;

#linker
#my @linkers =  ("1:CGCTCTTCCGATCT","2:GGCACATGCGTTCT","3:CGGTGAACCGTTCT","4:GGCTGTACGGATGT","5:CGGACTTGCGAAGT","6:CGGTGTTCGCATGT","7:GGCACTACCGTACT","8:CCCTCATGGCATCT");
#my $ln = $self->patient->barcode();
#my $ls = $linkers[$ln];
#my ($ln,$ls) = split(/:/,$self->patient->barcode());
	my @linker = (
		"1:CGCTCTTCCGATCT", "2:GGCACATGCGTTCT",
		"3:CGGTGAACCGTTCT", "4:GGCTGTACGGATGT",
		"5:CGGACTTGCGAAGT", "6:CGGTGTTCGCATGT",
		"7:GGCACTACCGTACT", "8:CCCTCATGGCATCT"
	);

	my $files_pe1  = file_util::find_file_pe( $self->patient, "" );
	my $count_lane = scalar(@$files_pe1);
	print $name. ":\n";
	my @jobs;
	my $fileout =
		$self->project->getVariationsDir("htlv1_calling") . "/"
	  . $name
	  . "-SIMPLIFIED_mergedIS.txt";

	#	 my $fileout = $dirout."/".$name;
	foreach my $cp (@$files_pe1) {
		my $file1 = $cp->{R1};
		my $file2 = $cp->{R2};

		if (   exists $self->fastq_files->{$file1}
			or exists $self->fastq_files->{$file2} )
		{
			die(
"same fastq file present in two different patient : $name $file1 $file2"
			);
		}
		$self->fastq_files->{$file1}++;
		$self->fastq_files->{$file2}++;
		die( "problem $file1 $file2 $dirin :" . $name )
		  unless -e $dirin . $file2;
		die( "problem $file1 $file2 $dirin :" . $name )
		  unless -e $dirin . $file1;

		my $f1 = $dirin . $file1;
		my $f2 = $dirin . $file2;

		# foreach my $l (@linker){
		my $bc   = $self->patient->barcode();
		my $ucbc = uc($bc);
		warn $ucbc;
		$ucbc =~ s/LIN//;
		warn $ucbc;
		my $l = $linker[ $ucbc - 1 ];

		my ( $ln, $ls ) = split( /:/, $l );
		my $out_dir = $name . "_" . $ln;

		my $fileout = $dirout . "/" . $name . "-SIMPLIFIED_mergedIS.txt";
		my $cmd =
"perl $bin_dev/htlv1.pl -project=$project_name  -patient=$name  -f1=$f1 -f2=$f2 $out_dir -ls=$ls -ln=$ln -fork=$ppn ";
		warn $cmd;
		my $type     = "insertion";
		my $stepname = $self->patient->name . "_" . $ln . "@" . $type;

		my $job_bds = job_bds_tracking->new(
			uuid         => $self->bds_uuid,
			cmd          => [$cmd],
			name         => $stepname,
			ppn          => $ppn,
			filein       => [ $f1, $f2 ],
			fileout      => $fileout,
			type         => $type,
			dir_bds      => $self->dir_bds,
			sample_name  => $name,
			project_name => $project_name,
			software     => $method
		);
		$self->current_sample->add_job( { job => $job_bds } );
		push( @jobs, $job_bds );

		if ( $self->unforce() && ( -e $fileout ) ) {
			$job_bds->skip();
		}

		#}
	}
	return ($fileout);
}

sub muc1 {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	my $dir_prod     = $project->getVariationsDir("vntyper") . "/muc1/";
	system("mkdir -p $dir_prod;chmod a+rwx $dir_prod") unless -e $dir_prod;
	my $fileout = $dir_prod . "/" . $name . ".vcf";

	$filein = $self->patient()->getBamFileName();    # unless $filein;

	my $ppn = 2;

	$ppn = int( $self->nproc / 2 ) if $self->nocluster;
	die( "-" . $filein ) unless $filein;

	my $bin_dev = $self->script_dir;
	my $version = $self->patient()->project->genome_version();
	my $cmd = "perl $bin_dev/muc1/vntyper.pl -project=$project_name  -patient=$name";
	my $type     = "muc1-calling";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "manta",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

sub advntr {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	my $dir_prod     = $project->getVariationsDir("advntr") ."/";
	#system("mkdir -p $dir_prod;chmod a+rwx $dir_prod") unless -e $dir_prod;
	my $fileout = $dir_prod . "/" . $name . ".vcf";

	$filein = $self->patient()->getBamFileName();    # unless $filein;

	my $ppn = 1;
	die( "-" . $filein ) unless $filein;

	my $bin_dev = $self->script_dir;
	my $version = $self->patient()->project->genome_version();
	my $cmd = "perl $bin_dev/muc1/advntr.pl -project=$project_name  -patient=$name";
	my $type     = "advntr-calling";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "advntr",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}
sub star_align {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	my $dir_prod     = $project->getVariationsDir("advntr") ."/";
	my $dir_pipeline = $project->getAlignmentPipelineDir("star-".$self->patient()->name);
	
	#dir_pipeline/$patient_name"."Aligned.sortedByCoord.out.bam
	#system("mkdir -p $dir_prod;chmod a+rwx $dir_prod") unless -e $dir_prod;
	my $fileout = "$dir_pipeline/$name".".Aligned.sortedByCoord.out.bam";


	my $ppn = 20 ;

	my $bin_dev = $self->script_dir;
	my $version = $self->patient()->project->genome_version();
	my $cmd = "perl $bin_dev/star/star_align.pl -project=$project_name  -patient=$name";
	my $type     = "star-align";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "star",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}


sub deepvariant {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $name         = $self->patient()->name();
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	my $low_calling  = "";
	$filein = $self->patient()->getBamFileName();    #unless $filein !~/bam/;
	my $fileout = $self->patient()->getVariationsFileName("deepvariant");
	#probleme du fichier de log non défini dans le pipeline lancé sans pbs
	#	my $cmd = "" ;
	my $ppn = $self->nproc;    # if $self->nocluster;
	$ppn = 20;
	$ppn =40 if $self->patient()->project->isGenome();
	my $real_ppn = $ppn;       #int($self->nproc / 2);
	$real_ppn = 40 if $self->host eq "morgan";
	die( "-" . $filein ) unless $filein;

	#	die($filein. " is empty") if (-z $filein);
	my $bin_dev = $self->script_dir;
	my $version = $self->patient()->project->genome_version();
	my $cmd = "perl $bin_dev/deepvariant/deepvariant.pl -version=$version -project=$project_name  -patient=$name -fork=$real_ppn   ";
	my $type     = "deepvariant";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		cmd          => ["$cmd"],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds,
		software     => "deepvariant",
		sample_name  => $self->patient->name(),
		project_name => $self->patient->getProject->name
	);
	$self->current_sample->add_job( { job => $job_bds } );

	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}


sub rnaseqsea_capture {
	my ( $self, $hash ) = @_;
	my $filein       = $hash->{filein};
	my $project      = $self->patient()->getProject();
	my $project_name = $project->name();
	my $name = $project->getPatients->[0]->name();
	my $ppn    = 40;
	my $method = "rnaseqsea_capture";
	my $dirout = $project->project_path . "/analysis/AllRes/";
	my $fileout = $dirout . "/allResRI.txt.gz";
	my $bin_dev = $self->script_dir;
	my $cmd_json = "$bin_dev/polyrnaseqsea/create_config_splices_analyse_file.pl -project=$project_name -force=1";
	my $json_file = `$cmd_json`;
	my $cmd = "Rscript $bin_dev/polyrnaseqsea/junctions/RNAseqSEA_capt_js_launch.r idprojet=$project_name fork=$ppn config_file=$json_file";
	$cmd .= " && $bin_dev/polyrnaseqsea/merge_all_junctions_files.pl -project=$project_name";
	my $type     = "rnaseqsea_capture";
	my $stepname = $self->patient->name . "@" . $type;
	my $job_bds  = job_bds_tracking->new(
		uuid         => $self->bds_uuid,
		software     => "",
		sample_name  => $name,
		project_name => $project_name,
		cmd          => [$cmd],
		name         => $stepname,
		ppn          => $ppn,
		filein       => [$filein],
		fileout      => $fileout,
		type         => $type,
		dir_bds      => $self->dir_bds
	);
	$self->current_sample->add_job( { job => $job_bds } );
	if ( $self->unforce() && -e $fileout ) {
		$job_bds->skip();
	}
	return ($fileout);
}

1;
