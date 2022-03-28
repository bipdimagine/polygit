#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo/lib/";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb/";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../GenBo/lib/kyoto/";
use lib "$Bin/../GenBo/lib/GenBoDB/writeDB";
use lib "$Bin/packages";


 
use PBS::Client;
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moose;
use MooseX::Method::Signatures;
use pipeline_steps;   
use file_util;
use Class::Inspector;
use Digest::MD5::File ;
use GBuffer;
use GenBoProject;
use colored; 

my $bin_cecile=qq{$Bin/../../cecile/polypipeline_scripts};
my $projectName ;

GetOptions(
	'project=s' => \$projectName
	);

unless ($projectName) {
	usage();
}

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );

my $dbh = $buffer->dbh();
$dbh->{AutoCommit} = 0;
$dbh->do("use Polypipeline;");
$ENV{'DATABASE'} = "";

my $dir = $project->getPipelineDir();
my $file = $dir."/".$project->name().".log";


my $steps_id = {
				"bwa_pe_clean"  => 1,
				"bwa_pe" => 2,
				"merge_bam_picard_illu" => 3,
				"sort_picard" => 4,
				"bamindex"  => 5,
				"readgroup_illumina" => 6,
				"covariate_illumina"  => 7,
				"realign1"  => 8,
				"realign2"  => 9,
				"move_bam" => 10,
				"rehead_bam" => 11,
				"coverage_samtools" => 12,
				"start_again"  => 13,
				"rmdup_picard" => 14,
				"mask_primer" => 15,
				"calling_merge" => 16,
				"calling_merge_varscan"  => 17,
				"move_vcf"  => 18,
				"stats"  => 19,
				"large_indels_finder"  => 20,
				"bam_sort_bamba" =>21,
				"rmdup_sambamba" => 22
};




fill_or_update_program_table($project, $dbh) ;
fill_or_update_reference_table($project, $dbh) ;
fill_or_update_versions_table ($project, $dbh) ;
fill_or_update_step_table($steps_id, $dbh);
fill_or_update_status_step_table($dbh);
fill_or_update_status_analysis_table($dbh);
$dbh->commit;
$dbh->disconnect;
print colored ['black ON_BRIGHT_BLUE'], "OK : Changes done in db !" ;


####Methods#####

sub fill_or_update_program_table {
	my ($project, $dbh)= @_ ;
	my @lprog = keys(%{$project->softwares_path});
	foreach my $prog (@lprog){
		my $md5 ;
		my $path= $buffer->getSoftware($prog);
		my @lprog = split('', $project->softwares_path->{$prog});
		if ($lprog[$#lprog] eq "/"){ 
			
			warn  colored::stabilo("red","program directory $buffer->getSoftware($prog) not found for $prog") unless -e $buffer->getSoftware($prog) ;
			$md5 = get_md5_dir($buffer->getSoftware($prog));
		} 
		else {
			warn colored::stabilo("red","program $prog not found : $buffer->getSoftware($prog) doesn't exist") unless -e $buffer->getSoftware($prog) ;
			$md5 = get_md5_string($buffer->getSoftware($prog));
		}
		#check md5 and path exist in Programs table and add program in Programs if not yet in table
		my $sql = qq{SELECT * from Polypipeline.`Programs` where `program_md5` = "$md5" ;};
		my $sth = $dbh->prepare($sql);
		$sth->execute();
		my $hprogram_id =  $sth->fetchall_hashref("program_id");
		my @tab_program_id = keys %$hprogram_id ;
		if (! defined $tab_program_id[0]) {
			print colored ['black ON_BRIGHT_RED'], "$prog md5 key is not yet in Programs table of Polypipeline database and will be added." ;
			print "\n";
			my $sql= qq{INSERT INTO Polypipeline.`Programs` ( `program_name`, `program_path`, `program_md5`) VALUES ("$prog" , "$path", "$md5");};
			$dbh->do($sql);
		}
	}
}

sub fill_or_update_versions_table {
	my ($project, $dbh)= @_ ;
	my @lpublic_data = keys(%{$project->buffer()->config->{'public_data'}});
	foreach my $public_data (@lpublic_data){
		my $path= $project->buffer()->config->{'public_data'}->{$public_data};
		#check path exist and add public_data in Versions if not in table
		my $sql = qq{SELECT * from Polypipeline.`Versions` where `version_name` = "$public_data" and `version_path` = "$path";};
		my $sth = $dbh->prepare($sql);
		$sth->execute();
		my $hpublic_data_id =  $sth->fetchall_hashref("version_id");
		my @tab_public_data_id = keys %$hpublic_data_id ;
		if (! defined $tab_public_data_id[0]) {
			print colored ['black ON_BRIGHT_RED'], "$public_data name and path is not yet in Versions table of Polypipeline database and will be added." ;
			print "\n";
			my $sql= qq{INSERT INTO Polypipeline.`Versions` ( `version_name`, `version_path`) VALUES ("$public_data" , "$path");};
			$dbh->do($sql);
		}
	}
	my @lensembl = grep{/ensembl_/} keys(%{$project->buffer()->config->{'ensembl_versions'}});
	foreach my $ensembl (@lensembl){
		my $ensembl_release= $project->buffer()->config->{'ensembl_versions'}->{$ensembl};
		#check path exist and add ensembl in ensembls if not in table
		my $sql = qq{SELECT * from Polypipeline.`Versions` where `version_name` = "$ensembl" and `version_path` = "$ensembl_release";};
		my $sth = $dbh->prepare($sql);
		$sth->execute();
		my $hensembl_id =  $sth->fetchall_hashref("version_id");
		my @tab_ensembl_id = keys %$hensembl_id ;
		if (! defined $tab_ensembl_id[0]) {
			print colored ['black ON_BRIGHT_RED'], "$ensembl name and path is not yet in Versions table of Polypipeline database and will be added." ;
			print "\n";
			my $sql= qq{INSERT INTO Polypipeline.`Versions` ( `version_name`, `version_path`) VALUES ("$ensembl" , "$ensembl_release");};
			$dbh->do($sql);
		}
	}
}


sub fill_or_update_reference_table {
	my ($project, $dbh)= @_ ;
	my $version = $project->getVersion();
	my $dir = $project->buffer()->config->{'public_data'}->{$version} . '/' ;
	my @lref_gatk = keys(%{$project->buffer->config->{'gatk_files'}});
	foreach my $file_ref (@lref_gatk){
		my $md5 ;
		my $path= $dir.$project->buffer->config->{'gatk_files'}->{$file_ref};
		my @lref_gatk = split('', $path);
		if ($lref_gatk[$#lref_gatk] eq "/"){ 
			warn colored::stabilo("red","file path not found : $path") unless -e $path ;
			$md5 = get_md5_dir($path);
		} 
		else {
			warn colored::stabilo("red","file not found : $path") unless -e $path ;
			$md5 = get_md5_string($path);
		}
		#check md5 and path exist in References table and add reference in Reference if not yet in table
		my $sql = qq{SELECT * from Polypipeline.`References` where `ref_md5` = "$md5" and `ref_path` = "$path";};
		my $sth = $dbh->prepare($sql);
		$sth->execute();
		my $href_id =  $sth->fetchall_hashref("ref_id");
		my @tab_ref_id = keys %$href_id ;
		if (! defined $tab_ref_id[0]) {
			print colored ['black ON_BRIGHT_RED'], "$file_ref md5 key and path is not yet in References table of Polypipeline database and will be added." ;
			print "\n";
			my $sql= qq{INSERT INTO Polypipeline.`References` ( `ref_name`, `ref_path`, `ref_md5`) VALUES ("$file_ref" , "$path", "$md5");};
			$dbh->do($sql);
		}
	}
}

sub fill_or_update_step_table {
	my ($steps_id, $dbh)= @_ ;
	foreach my $step (keys %$steps_id){
		my $step_id = $steps_id->{$step} ;
#		#check md5 and path exist in Programs table and add program in Programs if not yet in table
		my $sql = qq{SELECT * from Polypipeline.`Steps_infos` where `step_id` = "$step_id" ;};
		my $sth = $dbh->prepare($sql);
		$sth->execute();
		my $hstep_id =  $sth->fetchall_hashref("step_id");
		my @tab_step_id = keys %$hstep_id ;
		if (! defined $tab_step_id[0]) {
			print colored ['black ON_BRIGHT_RED'], "$step_id is not yet in Steps infos table of Polypipeline database and will be added." ;
			print "\n";
			my $sql= qq{INSERT INTO Polypipeline.`Steps_infos` ( `step_name`, `step_id`) VALUES ("$step" , "$step_id");};
			$dbh->do($sql);
		}
	}
}

sub fill_or_update_status_step_table {
	my  $dbh= shift ;
	my $steps_status = {
				"complete"  => 1,
				"in progress" => 2,
				"waiting" => 3,
				"error" => 4
	};
	foreach my $status (keys %$steps_status){
		my $status_id = $steps_status->{$status} ;
		my $sql = qq{SELECT * from Polypipeline.`Status_Steps` where `status_id` = "$status_id" ;};
		my $sth = $dbh->prepare($sql);
		$sth->execute();
		my $hstatus_id =  $sth->fetchall_hashref("status_id");
		my @tab_status_id = keys %$hstatus_id ;
		if (! defined $tab_status_id[0]) {
			print colored ['black ON_BRIGHT_RED'], "$status_id is not yet in Status_Steps table of Polypipeline database and will be added." ;
			print "\n";
			my $sql= qq{INSERT INTO Polypipeline.`Status_Steps` ( `status_name`, `status_id`) VALUES ("$status" , "$status_id");};
			$dbh->do($sql);
		}
	}
}


sub fill_or_update_status_analysis_table {
	my  $dbh= shift ;
	my $analysis_status = {
				"complete"  => 1,
				"in progress" => 2,
				"error" => 3
	};
	foreach my $status (keys %$analysis_status){
		my $status_id = $analysis_status->{$status} ;
		my $sql = qq{SELECT * from Polypipeline.`Status_Analysis` where `status_id` = "$status_id" ;};
		my $sth = $dbh->prepare($sql);
		$sth->execute();
		my $hstatus_id =  $sth->fetchall_hashref("status_id");
		my @tab_status_id = keys %$hstatus_id ;
		if (! defined $tab_status_id[0]) {
			print colored ['black ON_BRIGHT_RED'], "$status_id is not yet in Status_Analysis table of Polypipeline database and will be added." ;
			print "\n";
			my $sql= qq{INSERT INTO Polypipeline.`Status_Analysis` ( `status_name`, `status_id`) VALUES ("$status" , "$status_id");};
			$dbh->do($sql);
		}
	}
}

#calcul md5 key de la chaine de caractère date+taille
sub get_md5_string {
	my $f = shift ;
	my $cmd_date = qq{find $f -printf "%CH:%CM-%CD_" };
	$cmd_date =$cmd_date . "2>&1 ";
	my $cmd_size = qq{stat -c %s $f };
	$cmd_size= $cmd_size . "2>&1 ";

	my $string = `$cmd_date`;
	#my $string =system($cmd_date . "2>&1 "); #ne pas utiliser car accole le 0 du résultat de la commande system réussie à la sortie stout
	my $string2 = `$cmd_size`;
	my $stringtot = $string.$string2 ;
#	warn $stringtot ;

	my $md5 = Digest::MD5->new;
	$md5->add($stringtot );
	my $digest = $md5->hexdigest;
#	warn $digest ; # pour vérifier le résultat : faire echo $stringtot |md5sum (et pas echo -n )
	return $digest ;
	
}



#utile pour les fichiers de référence et exécutables uniques déclarés dans la section software du fichier de config
sub get_md5 {
	my $f = shift ;
	my $md5 = Digest::MD5->new;
    $md5->addpath($f) or die ($f."is not where you said: $!");
    my $digest = $md5->hexdigest;
	return $digest ;
}

#utile pour les dossiers déclarés dans la section software du fichier de config 
sub get_md5_dir {
	my $dir_path = shift ;
	my $md5 = Digest::MD5->new;
    $md5->adddir($dir_path);
    my $digest = $md5->hexdigest;
    return $digest ;
}


sub usage {
	print colored ['red'],  "\n======================= USAGE ========================\n";
	print colored ['blue'], $0." -project=project_name\n\n"; 
	die();
}