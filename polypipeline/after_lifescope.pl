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

my $bin_cecile=qq{$Bin/scripts/scripts_db_polypipeline};
my $bin_script_pipeline = qq{$Bin/scripts/scripts_pipeline};


my $projectName;
my $filename;

my $name;
my $patients_name;
my $steps_name;
my $force;
my $type;
my $fastq_ext;
my $exclude_patients;
my $max_cpu ;

 my @running_steps;
 my $predef_steps;
$predef_steps->{illumina} = ["alignment","rmdup","readgroup_illumina","covariate_illumina","move_bam","coverage","gvcf"];
$predef_steps->{bwa_frag} = ["bwa_frag","rmdup","realign","readgroup_illumina","covariate_illumina","move_bam","coverage","gvcf"];
$predef_steps->{hg38} = ["bwa_pe_hg38","merge_bamba","rmdup","readgroup_illumina","realign"];
#$predef_steps->{illumina} = ["bwa_pe","merge_bam_picard_illu","sort_bamba","rmdup_bamba","realign","readgroup_illumina","covariate_illumina","rehead_bam","move_bam","coverage","lifinder"];
$predef_steps->{diag} = ["alignment","rmdup","readgroup_illumina","realign","covariate_illumina","move_bam","coverage","calling_merge","move_vcf"];
$predef_steps->{full_genome} = ["alignment","rmdup","readgroup_illumina","realign","covariate_illumina","move_bam","coverage",];

$predef_steps->{diag_primer} = ["alignment","mask_primer","bam_sort","readgroup_illumina","realign","covariate_illumina","move_bam","coverage","calling_merge","move_vcf"];
$predef_steps->{diag_pcr} = ["alignment","readgroup_illumina","realign","covariate_illumina","move_bam","coverage","calling_merge","move_vcf"];
$predef_steps->{miseq_low_calling} = ["alignment","merge_bam_picard_illu","bam_sort","readgroup_illumina","covariate_illumina","bamindex","move_bam","coverage","calling_merge_low","move_vcf"];

$predef_steps->{hegp_pcr} = ["alignment","mask_primer_se","bam_sort","covariate_illumina","readgroup_illumina","move_bam","coverage", "calling_merge","move_vcf"];
$predef_steps->{hegp_capture} = ["alignment","rmdup","readgroup_illumina","realign","covariate_illumina","move_bam","coverage", "calling_merge","move_vcf"];
$predef_steps->{calling} =["calling_merge","move_vcf"];
$predef_steps->{rna_seq} =["alignment","rmdup","move_bam"];
$predef_steps->{calling_rna_seq} =["alignment", "rmdup","splitntrim","covariate_illumina","move_bam","gvcf"];
GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'exclude=s' => \$exclude_patients,
	'steps=s' => \$steps_name,
	'force=s' => \$force,
	'type=s' => \$type,
	'max_cpu=s' => \$max_cpu,
);


my $buffer = GBuffer->new();

unless ($projectName) {
	usage();
}
my $project = $buffer->newProject( -name => $projectName );

my $dbh = $buffer->dbh();
$dbh->{AutoCommit} = 0;
$dbh->do("use Polypipeline;");
$ENV{'DATABASE'} = "";


my $jobs;
my $dir = $project->getPipelineDir();
my $file = $dir."/".$project->name().".log";
my $pipeline = pipeline_steps->new(log_file=>$file, project=>$project, dbh=>$dbh,patient_argument=>$patients_name);
$pipeline->fastq_extend($fastq_ext) if $fastq_ext;

$pipeline->max_cpu($max_cpu) if $max_cpu ;


my $steps = {
				"readgroup_illumina"=> sub {$pipeline->read_group_illumina(@_)},
				"move_bam"=> sub {$pipeline->move_bam(@_)},
				"bwa_pe_hg38"=> sub {$pipeline->bwa_pe_hg38(@_)},
				"merge_bam"=> sub {$pipeline->merge_bam(@_)},
				"bam_sort"=> sub {$pipeline->bam_sort(@_)},
				"bamindex" =>sub {$pipeline->bamindex_samtools(@_)},
				"rmdup"=> sub {$pipeline->rmdup(@_)},
				"realign" =>sub {$pipeline->realign_gatk(@_)},
				"covariate_illumina" => sub {$pipeline->covariate_illumina(@_)},
				"coverage" => sub {$pipeline->coverage_samtools(@_)},
				"bcf" => sub {$pipeline->fast_bcf(@_)},
				"purge_fastq" => sub {$pipeline->purge_fastq(@_)},
				"cutadapt" => sub {$pipeline->cutadapt(@_)},
				"align_tmap"=> sub {$pipeline->align_tmap(@_)},
				"add_chr"=> sub {$pipeline->add_chr(@_)},
				"start_again" => sub {$pipeline->start_again(@_)},
				"chrM_chrMT"=> sub {$pipeline->chrM_chrMT(@_)},
				"reorder_picard"=> sub {$pipeline->reorder_picard(@_)},
				"calling_merge" => sub {$pipeline->calling_merge(@_)},
				"calling_merge_low" => sub {$pipeline->calling_merge_low(@_)},
				"calling_merge_varscan" => sub {$pipeline->calling_merge_varscan(@_)},
				"move_vcf" => sub {$pipeline->move_vcf(@_)},
				"move_snp_vcf" => sub {$pipeline->move_snp_vcf(@_)},
				"move_bam_covariate" => sub {$pipeline->move_bam_after_covariate(@_)},
				"mask_primer" => sub {$pipeline->mask_primer(@_)},
				"mask_primer_se" => sub {$pipeline->mask_primer_start_end(@_)},
				"haplotypeCaller" => sub {$pipeline->haplotypeCaller(@_)},
				"delete_chrM" => sub {$pipeline->delete_chrM(@_)},
				"stats" => sub {$pipeline->picard_stats(@_)},
				"rehead_bam" => sub {$pipeline->rehead_bam(@_)},
				"launch_analysis" => sub {$pipeline->launch_analysis(@_)},
				"merge_bam" => sub {$pipeline->merge_bam(@_)},
				"gvcf" => sub {$pipeline->calling_gvcf(@_)},
				"splitntrim"=>  sub {$pipeline->SplitNCigarReads(@_)},
				"test"=>  sub {$pipeline->test(@_)},
				"alignment"=>  sub {$pipeline->alignment(@_)},
			};




if ($type && $steps_name eq "all") {
	if (exists $predef_steps->{$type}){
			 @running_steps =@{$predef_steps->{$type}};
	}
	else {
		usage();
	}		
}
unless ($steps_name) {
	usage();
}

usage() unless $patients_name;

my $patients = file_util::return_patients( $project, $patients_name );
confess("no patients") unless scalar(@$patients);



$pipeline->unforce(0) if $force;



my $ind_steps  = {
	"coverage" => sub {$pipeline->coverage(@_)},
	"bcf" => sub {$pipeline->bcf(@_)},
	 "stats" => sub {$pipeline->picard_stats(@_)},
};


if ($steps_name ne "all") {
	@running_steps = split(",",$steps_name);
}

if (defined $steps_name && $steps_name ne "all") {
	my @list_steps = split(",",$steps_name);
	foreach my $s (@list_steps){
		unless (grep{/$s/} keys %$steps ){
			usage() ;
		}
	}
}


foreach my $s (keys %$steps){
	my ($f) = grep {$_ =~ /$s/} @running_steps;
	delete $steps->{$s} unless $f;	
	my ($f2) = grep {$_ =~ /$s/} @running_steps;
	delete $ind_steps->{$s} unless $f;	
}


	

my @jobs;
my @exclude;
@exclude = split(",",$exclude_patients) if ($exclude_patients);
my @lpatients ;
my @lp_ids ;
foreach my $p (@$patients) {
	my $find = grep {$p->name eq $_ } @exclude;
	next if $find;
	warn $p->name;
	push(@lpatients, $p->name);
	push(@lp_ids, $p->id);
}




#construction de l'analyse id
my $analysis_id =add_analysis_in_db($project, $dbh);
warn "Analysis id : ".$analysis_id ;

#verif prog et fichiers de références sont conformes à ceux rentrés dans la db
read_config_prog_check_md5_from_db($project, $dbh, $analysis_id) ;
read_version_project_and_add_version_id_analysis_id($project, $dbh, $analysis_id) ;
read_config_ref_files_check_md5_from_db($project, $dbh, $analysis_id) ;
#die();

##############################################################
# Lancement du pipeline : cas 1 sans cluster (option max_cpu)  cas 2 avec cluster et protocole PBS #
##############################################################
print "\n";
print "\n";
print "\n";

print colored ['black '], "################################################################################################" ;
			print "\n";
print colored ['black '], "################################################################################################" ;
print "\n";
print "\n";
if ($max_cpu){
	foreach my $p (@$patients) {
	my $find = grep {$p->name eq $_ } @exclude;
	next if $find;
	push(@jobs,prepare_jobs_by_patient($p,$pipeline));
	}

	$pipeline->print_all_steps;
	
#	warn Dumper $pipeline->hash_cmds() ;
#	warn Dumper $pipeline->hash_steps() ;
print "\n";
print colored ['black '], "################################################################################################" ;
			print "\n";
	warn "nb patient : ".scalar(@lpatients)."\n";
	warn "nb step : ".scalar(keys %{$pipeline->hash_steps()})."\n";
print colored ['black '], "################################################################################################" ;
print "\n";
	

	my $choice = prompt("run this/these step(s)   (y/n) ? ");
	die() if ($choice ne "y"); 
	$dbh->commit;

	#Patient loop
	foreach my $p (@$patients) {
		#log_file 
		my $patient_name = $p->name;
		colored::stabilo("blue",$p->name);
		my $log_file = $patient_name."_log.txt";
		#order list of step
		my $h_patient_steps ;
		foreach my $stepname (keys %{$pipeline->hash_cmds()->{$p->name}}){
		
#			warn $pipeline->hash_steps()->{$stepname}->{$p->name}->{order};
			my $order = $pipeline->hash_steps()->{$stepname}->{$p->name}->{order};
			$h_patient_steps->{$order}=$stepname ;
		} ;

		foreach my $order (sort {$a<=>$b} keys %{$h_patient_steps}){
			my $stepname = $h_patient_steps->{$order} ;
			colored::stabilo("yellow", $patient_name." : ".$stepname) ;
			my $cmd= $pipeline->hash_cmds->{$p->name}->{$stepname};
			my $c;
			my $step_id=$pipeline->hash_steps()->{$stepname}->{id};
			eval {
				foreach  $c (split("&&",$cmd)){
				print "####  ".$c." ####1\n";	
				my $r = system($c);
				my $r;
				#warn $c;
				warn $c unless $r==0;
				my $cmd_log = qq{echo "\nerror $patient_name  :  $stepname" >>$log_file ; echo ++++++++++++>>$log_file ; };
				#redirection de l'erreur dans le log
				system($cmd_log.$c." 2>>$log_file" ) unless $r==0;  #lignes d'erreur : pb c'est du binaire, comment le récupérer pour le mettre dans le fichier de log ? avec 2>>$log_file
				
#				my $db_end = qq{perl $bin_cecile/change_step_status_end_in_db.pl -analysis_id=$analysis_id -step_id=$step_id -status=4 -patient=$patient_name;};
#				system($db_end) unless $r ==0;
				die()  unless $r==0;
				};
				colored::stabilo("green","ok");
			};

			if ($@){
				colored::stabilo("red",$p->name." : $stepname not ok => see error in ".$log_file) ;
				my $db_end = qq{perl $bin_cecile/change_step_status_end_in_db.pl -analysis_id=$analysis_id -step_id=$step_id -status=4 -patient=$patient_name;};
				system($db_end);
				last;
			};
		};
	};
	#changer statut anaylyse quand terminé. Faire requette sur les statuts des étapes et si toutes OK : complete, sinon, ERROR
	my $command_analysis_final_status = qq{perl $bin_cecile/set_analysis_status_end_in_db.pl -analysis_id=$analysis_id ;};
	system($command_analysis_final_status);
	my $quality = "$bin_script_pipeline/quality_check_target_gene.pl -project=$projectName";

	$dbh->commit;
	$dbh->disconnect;
	warn $quality;
	system($quality);
	
}

else {
	warn "Launching pipeline all patient together, jobs are submitted using PBS protocol" ;
	foreach my $p (@$patients) {
		my $find = grep {$p->name eq $_ } @exclude;
		next if $find;
		push(@jobs,prepare_jobs_by_patient($p,$pipeline));
	}

	$pipeline->print_all_steps;
	warn "analysis id : ".$analysis_id ;
	warn "Total patients : ".scalar(@lpatients);
	warn "Total steps : ".scalar(keys %{$pipeline->hash_steps()});
	
	my $choice = prompt("run this/these step(s)   (y/n) ? ");
	die() if ($choice ne "y"); 
	$dbh->commit;
	
	
	my $client = PBS::Client->new();
	my $final_jobs = $pipeline->complete(previous=>\@jobs);
	$client->qsub($final_jobs);
	

	print colored ['black ON_BRIGHT_GREEN'],"the log file is here : $file";
	print color 'reset';
	print "\n";
	
	$dbh->commit;
	$dbh->disconnect;
};



exit(0);


########
# Methods #
########

sub prepare_jobs_by_patient {
	my ($p)  = @_;
	my $name = $p->name();
	my $projecName = $project->name ;
	$pipeline->patient($p);
	$pipeline->{$name}->{analysis_id}=$analysis_id ;
		
	
	#dès que $running_steps->patient existe, écrire dans la table des étapes, statut waiting
	my $current_job = $pipeline->start_job();
	my $next_file="";
	foreach my $step (@running_steps){
		warn "c est quoi cette etape ?  $step " unless exists $steps->{$step};
		($current_job,$next_file) = $steps->{$step}->((filein=>$next_file,previous=>$current_job));
	}
	
	#remplissage de la table Steps
	my $c =1 ;
	foreach my $step (@{$pipeline->running_steps->{$name}}){
		my $step_name = $step ;
		my $step_id = $pipeline->hash_steps->{$step}->{"id"};
		my $step_order = $c ;
		$pipeline->hash_steps->{$step}->{$name}->{order} = $step_order ;
		add_steps_in_db($analysis_id, $step_id, $step_order,$name, $dbh);
		$c ++ ;
	}
	
	my $final_job = $pipeline->end_job();
	$final_job->prev({ ok => [$current_job]} );
	return $final_job;
}


sub add_analysis_in_db {
	my ($project, $dbh) = @_ ;
	my $project_id = $project->id;
	my $project_name = $project->name;
	my $user =$ENV{'LOGNAME'};
	my $date = localtime();
#	my $status = "in progress" ;
	my $status = 2 ;
	my $list_steps = join(",", @running_steps);
	my $list_patients = join(",", @lpatients);
	my $list_p_ids = join(",", @lp_ids);
	
	my $sql= qq{INSERT INTO Polypipeline.`Analysis` (`project_name`,`project_id`, `user`, `status`, `date`, type, pipeline_steps, patients_names, patients_ids) VALUES ("$project_name","$project_id", "$user", "$status", "$date","$type", "$list_steps", "$list_patients", "$list_p_ids");};
	$dbh->do($sql);

	#get analysis_id
	my $sql2 = qq{SELECT * FROM Polypipeline.`Analysis` where `project_id`="$project_id" and  `date`= "$date" and `user`= "$user";};
	
	my $sth = $dbh->prepare($sql2);
	$sth->execute();
	my $analysis_id =  $sth->fetchall_hashref("analysis_id");
	my @tab_analysis_id = keys %$analysis_id ;
	warn "Not a unique analysis_id for this patient/project/date/user analysis" unless (scalar(@tab_analysis_id) ==1 );
	
	return $tab_analysis_id[0] ;
}

sub add_steps_in_db {
	my ($analysis_id, $step_id, $step_order,$patient_name, $dbh) = @_ ;
#	my $step_status="waiting" ;
	my $step_status=3 ;

	my $sql5= qq{INSERT INTO Polypipeline.`Analysis_Steps` (`analysis_id`,`step_id`, `step_order`, `patient`,`status`) VALUES ("$analysis_id","$step_id", "$step_order","$patient_name", "$step_status");};
	$dbh->do($sql5);
	
}

sub read_config_prog_check_md5_from_db {
	my ($project, $dbh, $analysis_id) = @_ ;
	my @lprog = keys(%{$project->softwares_path});
	foreach my $prog (@lprog){
		my $md5 ;
		my @lprog = split('', $project->softwares_path->{$prog});
		if ($lprog[$#lprog] eq "/"){ 
			$md5 = get_md5_dir($buffer->getSoftware($prog));
		} 
		else {
			$md5 = get_md5_string($buffer->getSoftware($prog));
		}
		#check md5 exist in Programs and get id
		my $sql = qq{SELECT * from Polypipeline.`Programs` where `program_md5` = "$md5" ;};
		my $sth = $dbh->prepare($sql);
		$sth->execute();
		my $hprogram_id =  $sth->fetchall_hashref("program_id");
		my @tab_program_id = keys %$hprogram_id ;
		if (! defined $tab_program_id[0] ) {
			print colored ['black ON_BRIGHT_CYAN'], "$prog md5 key is not recognized in program database : if it is a new file corresponding to a new version, you should update db with init_db_polypipeline.pl ; if not, original file has been changed and is not conform. Check the file." ;
			print "\n";
		}
		else{
		my $program_id = $tab_program_id[0] ;
		#add program in Analysis_programs 
		my $sql2= qq{INSERT INTO Polypipeline.`Analysis_programs` (`analysis_id`, `program_id`) VALUES ("$analysis_id","$program_id");};
	
		$dbh->do($sql2);
		}
	}
}

sub read_config_ref_files_check_md5_from_db {
	my ($project, $dbh, $analysis_id) = @_ ;
	my $version = $project->getVersion();
	my $dir = $project->buffer()->config->{'public_data'}->{$version} . '/' ;
	my @lref_gatk = keys(%{$project->buffer->config->{'gatk_files'}});
	foreach my $file_ref (@lref_gatk){
		my $md5 ;
		my $path= $dir.$project->buffer->config->{'gatk_files'}->{$file_ref};
		my @lref_gatk = split('', $path);
		if ($lref_gatk[$#lref_gatk] eq "/"){ 
			$md5 = get_md5_dir($path);
		} 
		else {
			$md5 = get_md5_string($path);
		}
		#check md5 exist in References table and get id
		my $sql = qq{SELECT * from Polypipeline.`References` where `ref_md5` = "$md5" and `ref_path` = "$path";}; #md5 not unique=> md5 + path unique
		my $sth = $dbh->prepare($sql);
		$sth->execute();
		my $href_id =  $sth->fetchall_hashref("ref_id");
		
		my @tab_ref_id = keys %$href_id ;
		
		if (! defined $tab_ref_id[0] ) {
			print colored ['black ON_BRIGHT_CYAN'], "$file_ref md5 key is not recognized in References table of Polypipeline database : \n if it is a new file corresponding to a new version, you should update db with init_db_polypipeline.pl ; \n if not, original file has been changed and is not conform. Check the file." ;
			print "\n";
		}
		else{
			my $ref_id = $tab_ref_id[0] ;
			#add program in Analysis_programs 
			my $sql2= qq{INSERT INTO Polypipeline.`Analysis_references` (`analysis_id`, `ref_id`) VALUES ("$analysis_id","$ref_id");};
			$dbh->do($sql2);
		}
	}
}

sub read_version_project_and_add_version_id_analysis_id {
	my ($project, $dbh, $analysis_id) = @_ ;
	my $version = $project->getVersion();
	#check version exist in Versions table and get id
		my $sql = qq{SELECT * from Polypipeline.`Versions` where `version_name` = "$version" ;};
		my $sth = $dbh->prepare($sql);
		$sth->execute();
		my $hversion_id =  $sth->fetchall_hashref("version_id");
		my @tab_version_id = keys %$hversion_id ;
		if (! defined $tab_version_id[0] ) {
			print colored ['black ON_BRIGHT_RED'], "$version is not in Versions table of Polypipeline database : \n if it is a new one , you should update db with init_db_polypipeline.pl." ;
			print "\n";
		}
		else{
			my $version_id = $tab_version_id[0] ;
			#add program in Analysis_programs 
			my $sql2= qq{INSERT INTO Polypipeline.`Analysis_version` (`analysis_id`, `version_id`) VALUES ("$analysis_id","$version_id");};
			$dbh->do($sql2);
		}	
	my $ensembl = $project->ensembl_version();
	
	#check ensembl exist in Versions table and get id
		my $sql3 = qq{SELECT * from Polypipeline.`Versions` where `version_path` = "$ensembl" ;};
		my $sth3 = $dbh->prepare($sql3);
		$sth3->execute();
		my $hensembl_id =  $sth3->fetchall_hashref("version_id");
		my @tab_ensembl_id = keys %$hensembl_id ;
		if (! defined $tab_ensembl_id[0] ) {
			print colored ['black ON_BRIGHT_CYAN'], "$ensembl is not in Versions table of Polypipeline database : \n if it is a new one , you should update db with init_db_polypipeline.pl." ;
			print "\n";
		}
		else{
			my $ensembl_id = $tab_ensembl_id[0] ;
			#add program in Analysis_programs 
			my $sql4= qq{INSERT INTO Polypipeline.`Analysis_version` (`analysis_id`, `version_id`) VALUES ("$analysis_id","$ensembl_id");};
			$dbh->do($sql4);
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

#utile pour les dossiers déclarés dans la section software du fichier de config (picard et vcftools qui n'ont pas un exécutable unique et qu'on récupère par une méthode qui va chercher le path du dossier)
sub get_md5_dir {
	my $dir_path = shift ;
	my $md5 = Digest::MD5->new;
    $md5->adddir($dir_path);
    my $digest = $md5->hexdigest;
    return $digest ;
}



sub add_run_and_skipped_steps_in_analysis {
	my ($analysis_id,$run_step_p,$skip_step_p, $dbh) = @_ ;
	my $sql4= qq{UPDATE Polypipeline.`Analysis` SET `run_steps`="$run_step_p",`skip_steps`="$skip_step_p" WHERE `analysis_id` = $analysis_id ;};
	$dbh->do($sql4);
}



sub usage {
	print colored ['red'],  "\n======================= USAGE ========================\n";
	print colored ['blue'], $0." -project=project_name -steps=(step or all) -patients=(patient or all) -type =(illumina, miseq, miseq_primer or hegp_pcr, hegp_capture) -force=(1 : force restart step) -cpu_max = (available CPU number if run without PBS protocol) \n";
	print colored ['green'], "type : illumina ===> steps : " .join(",",@{$predef_steps->{illumina}})."\n";
	print colored ['green'], "type : illumina diag ===> steps : " .join(",",@{$predef_steps->{illumina_diag}})."\n";
	print colored ['green'], "type : miseq ===> steps : " .join(",",@{$predef_steps->{miseq}})."\n";
	print colored ['green'], "type : miseq primer ===> steps : " .join(",",@{$predef_steps->{miseq_primer}})."\n";
	print colored ['green'], "type : hegp_pcr ===> steps : " .join(",",@{$predef_steps->{hegp_pcr}})."\n";
	print colored ['green'], "type : hegp_capture ===> steps : " .join(",",@{$predef_steps->{hegp_capture}})."\n";
	print colored ['green'], "type : calling ===> steps : " .join(",",@{$predef_steps->{calling}})."\n";
	print colored ['blue'], "================== steps list ================\n";
	print colored ['green'], join(", ",keys %$steps)."\n";
	print colored ['blue'], "===================================================\n";
	if (defined $steps_name && $steps_name ne "all") {
	my @list_steps = split(",",$steps_name);
	foreach my $s (@list_steps){
		unless (grep{/$s/} keys %$steps ){
			print colored ['red'], $s." is not a valid step.\n";
		}
	}
}
print colored ['red'],"=================================================\n";
	die();
}

