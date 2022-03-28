#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use Data::Dumper;
use Getopt::Long;
use JSON::XS;
use KyotoCabinet;
use GBuffer; 
use Carp;





my ($project_name, $patient_actual, $patient_new ) ;

GetOptions( 'project=s' => \$project_name,
			'patient_actual=s' => \$patient_actual,
			'patient_new=s' => \$patient_new,
);




die("\n\nERROR: -project missing... Die...\n\n") unless ($project_name);
die("\n\nERROR: -patient_actual missing... Die...\n\n") unless ($patient_actual);
die("\n\nERROR: -patient_new... Die...\n\n") unless ($patient_new);

#rename_one_patient_in_db($project_name, $patient_actual, $patient_new);

check_change_in_db($project_name,$patient_new);



### Methods ####


sub rename_one_patient_in_db {
	my ($project_name, $ancien, $newName) = @_ ;
	my $buffer = GBuffer->new();
	my $dbh = $buffer->dbh();
	$dbh->{AutoCommit} = 0;
	$dbh->do("use PolyprojectNGS;");

	my $project = $buffer->newProject( -name => $project_name );
	my $idProject = $project->id();
	my $patient = $project->getPatient($ancien) ; ###si patient inconnu, plante ici, et pas de test.
	unless ($patient ne undef){
		 warn "$ancien est inconnu dans le projet. \nListe des patients du projet : ". join(" ",map {$_->name} @{$project->getPatients()});
		 die();
	}
	my $idPat= $patient->id();
	my $patname = $patient->name();
	
	warn "ID project: ".$idProject." et ID patient:".$idPat."et name: ".$patname."\n" ;

		####check bam and vcf files exist
		my $bamFile = $patient-> getBamFile();
		die ("No bam file for patient ".$patient->name()) unless -e $bamFile ;

		warn "coucou";
		####rename cov, bam and vcf, and index files
		my $m = $patient->alignmentMethod();
		my $dir_bam = $project->getAlignmentDir($m);
		my $NbamFile = $dir_bam."/".$newName.".bam";
		warn $NbamFile ;
		my $cmd1 = system(`mv $bamFile $NbamFile`);
		my $cmd2 = system(`mv $bamFile".bai" $NbamFile".bai"`);
		my $recal_file = $dir_bam."/recal/".$ancien.".recal.table" ;
		my $Nrecal_file = $dir_bam."/recal/".$newName.".recal.table" ;
		my $cmd2bis = system(`mv $recal_file $Nrecal_file`);
	
		my $cov_dir = $project->getCoverageDir();	

		opendir(DIR, $cov_dir); 
		my @covFiles = readdir(DIR);
		closedir(DIR);
		foreach my $cov_file (@covFiles){

			if ($cov_file =~ $ancien){
				my $path_cov_file =  $cov_dir.'/'.$cov_file ;
				my $cov_file2 = $cov_file ;
				$cov_file2 =~ s/$ancien/$newName/g ;
				warn my $path_cov_file2 = $cov_dir.'/'.$cov_file2 ;
				my $cmd3 = system(`mv $path_cov_file $path_cov_file2`);
			} 
		}
	
		my $lMethodName = $project->getCallingMethods();
		foreach my $m2 (@$lMethodName) {
			my $dir_var = $project->getVariationsDir($m2);	
			my $VarFile = $dir_var."/".$ancien.".vcf.gz";
			my $VarFileTbi = $dir_var."/".$ancien.".vcf.gz.tbi";
			my $LogVarFile = $dir_var."/".$ancien.".log" ;
			my $NVarFile = $dir_var."/".$newName.".vcf.gz";
			my $NVarFileTbi = $dir_var."/".$newName.".vcf.gz.tbi";
			my $NLogVarFile = $dir_var."/".$newName.".log" ;
	#		warn $NVarFile;
			my $cmd4 = system(`mv $VarFile $NVarFile `);
			my $cmd5 = system(`mv $VarFileTbi $NVarFileTbi `);
			my $cmd6 = system(`mv  $LogVarFile $NLogVarFile`);
			my $dir_indels = $project->getIndelsDir($m2);	
			my $IndelFile = $dir_indels."/".$ancien.".vcf.gz";
			my $IndelFileTbi = $dir_indels."/".$ancien.".vcf.gz.tbi";
			my $LogIndelFile = $dir_indels."/".$ancien.".log" ;
			my $NIndelFile = $dir_indels."/".$newName.".vcf.gz";
			my $NIndelFileTbi = $dir_indels."/".$newName.".vcf.gz.tbi";
			my $NLogIndelFile = $dir_indels."/".$newName.".log" ;
	#		warn $NIndelFile;
			my $cmd7 = system(`mv $IndelFile $NIndelFile `);
			my $cmd8 = system(`mv  $IndelFileTbi $NIndelFileTbi`);
			my $cmd9 = system(`mv  $LogIndelFile $NLogIndelFile`);
			
	####replace name in db
	my $sql = qq{UPDATE `patient` SET `name`="$newName", `origin`="$newName", `family`="" WHERE `patient_id`="$idPat";};
	$dbh->do($sql);	
	$dbh->commit();

		
	}
}

sub check_change_in_db{
	my ($project_name, $patient_new) = @_ ;
	my $buffer2 = GBuffer->new();
	my $project = $buffer2->newProject( -name => $project_name );
	my $idProject = $project->id();
	my $dbh2 = $buffer2->dbh();
	$dbh2->do("use PolyprojectNGS;"); 
	
	my $sql = qq{SELECT `project_id`, `patient_id`,`name`, `origin`, `family` FROM `patient` where `project_id`="$idProject" and  `name`= "$patient_new" ;};
	my $sth = $dbh2->prepare($sql);
	$sth->execute();
	my $res = $sth->fetchall_arrayref();

	warn "my res".Dumper($res) ;
	
	my @db_patient_newName = keys $sth->fetchall_hashref(`name`);
	warn Dumper(@db_patient_newName);
	unless (scalar(values $res) != 0) {
		confess("New patient is not in database. Problem");
	}

	
};



