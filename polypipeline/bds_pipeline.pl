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
use Logfile::Rotate;
 use Cwd;
use PBS::Client;
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moo;

use bds_steps;   
use file_util;
use Class::Inspector;
use Digest::MD5::File ;
use GBuffer;
use GenBoProject;
use colored; 
use Config::Std;
use Text::Table;

use File::Temp qw/ tempfile tempdir /;; 
use YAML::Syck;
use YAML::Dumper;
use Term::Menus;
 use Proc::Simple;
 use Storable;
use JSON::XS;
 
 

 
  
my $bin_cecile=qq{$Bin/scripts/scripts_db_polypipeline};
my $bin_script_pipeline = qq{$Bin/scripts/scripts_pipeline};


my $projectName;
my $filename;

my $name;
my $patients_name;
my $steps_name;
my $force;
my $type ="";
my $fastq_ext;
my $exclude_patients;
my $max_cpu ;
my $bds;
 my @running_steps;
 my $predef_steps;
 my $nocluster = 0;
 my $myproc;
my $low_calling;
my $predef_type;
my $define_steps;
my $yes =0;
my $cnv =0;


#$define_steps->{pipeline}->{all} = ["alignment","elprep","move_bam","coverage","gvcf4","callable_regions","binary_depth"];
#$define_steps->{pipeline}->{all_genome} = ["alignment","elprep","move_bam","gvcf4","callable_regions","binary_depth","breakdancer","manta","canvas","wisecondor","calling_wisecondor"];
#$define_steps->{pipeline}->{cnv} = ["breakdancer","manta","canvas","wisecondor","calling_wisecondor"];#,"manta","canvas","wisecondor","calling_wisecondor"];
#
#$define_steps->{pipeline}->{diag_capture} = ["alignment","elprep","move_bam","coverage","gvcf4","calling_panel","binary_depth"];
#$define_steps->{pipeline}->{diag_primer} = ["alignment","mask_primer","bam_sort","readgroup","realign_recal","move_bam","coverage","gvcf4","calling_panel","binary_depth"];
#$define_steps->{pipeline}->{diag_pcr} = ["alignment","readgroup","realign_recal","move_bam","coverage",'gvcf4',"calling_panel","binary_depth"];
#$define_steps->{pipeline}->{picard_stat} = ["stats"];
#$define_steps->{pipeline}->{diag_mito} = ["alignment","rmdup","readgroup","move_bam","coverage","calling_panel","binary_depth"];
#$define_steps->{pipeline}->{just_alignement} = ["alignment","elprep","move_bam"];
#$define_steps->{pipeline}->{gvcf_binary_depth} = ["gvcf4","binary_depth"];
#$define_steps->{pipeline}->{rna_seq} =["alignment","rmdup","move_bam","rnaseq_metrics","binary_depth"];
#$define_steps->{pipeline}->{rnaseq_umi} =["alignment","rmdup_nudup","move_bam","rnaseq_metrics","binary_depth"];
#$define_steps->{pipeline}->{exome_umi} =["concat_fastq_umi","fastq_to_bam","annotate_with_umi","run_alignment_umi","merge_bam_ubam","group_reads_by_umi","call_consensus_reads","filter_consensus_read","bam_to_fastq_umi","move_bam","calling_panel","coverage","binary_depth"];
#$define_steps->{pipeline}->{qiagen} =["reorder_picard","readgroup","move_bam","coverage","binary_depth"];
#
#
#$define_steps->{calling}->{all} = ["genotype_gvcf4","correct_vcf","move_vcf_hc4","dude"];
#$define_steps->{calling}->{genotype_and_move} = ["genotype_gvcf4","correct_vcf","move_vcf_hc4"];
#$define_steps->{calling}->{dude} = ["dude"];
#$define_steps->{transcripts}->{all} = ["transcripts_dude","genes_dude","transcripts_coverage"];
#$define_steps->{pipeline_diag} = ["alignment","elprep","move_bam","coverage","gvcf4","callable_regions","binary_depth"];
$predef_steps->{getFastqFromBam}=["bam_to_fastq"];



#$predef_steps->{exome_umi} =["concat_fastq_umi","fastq_to_bam","annotate_with_umi","run_alignment_umi","merge_bam_ubam","group_reads_by_umi","call_consensus_reads","filter_consensus_read","bam_to_fastq_umi","move_bam","calling_panel","coverage",binary_depth"];
#$predef_steps->{calling_rna_seq} =["alignment", "rmdup","splitntrim","covariate","move_bam","gvcf4","callable_regions"];
my $filename_cfg;
my $limit;
my $version;
my $arg_steps;
my $pipeline_name;
my $secret;
my $pad;

GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'exclude=s' => \$exclude_patients,
	'steps=s' => \$arg_steps,
	'force=s' => \$force,
	'type=s' => \$type,
	'max_cpu=s' => \$max_cpu,
	'bds=s' => \$bds,
	'nocluster=s' => \$nocluster,
	'config=s' => \$filename_cfg,
	'limit=s' => \$limit,
	'version=s' => \$version,
	'pipeline=s' =>\$pipeline_name,
	'yes=s' =>\$yes,
	'cnv=s' =>\$cnv,
	'nolimit=s' =>\$secret,
	#'low_calling=s' => \$low_calling,
	'padding=s' =>\$pad,
);
$patients_name = "all" unless $patients_name;
my $report;
my $buffer = GBuffer->new();

unless ($filename_cfg) {
 $filename_cfg = $buffer->config_directory()."/pipeline/pipeline.cfg";
}
if ($cnv == 1) {
	$define_steps = {}; 
	$define_steps->{pipeline}->{short_read} = "binary_depth,manta,cnvnator,canvas,wisecondor,calling_wisecondor,hificnv";
	$define_steps->{pipeline}->{pacbio} = "binary_depth,hificnv";
}
else {
	read_config $filename_cfg =>  %{$define_steps};
}
unless ($projectName) {
	usage();
}
my $project = $buffer->newProject( -name => $projectName, -version=>$version );

#if($HG38) {
#	$project->genome_version("HG38_CNG");
#	$project->version("HG38_CNG");
#}
my $pipeline = bds_steps->new( project=>$project,argument_patient=>$patients_name,nocluster=>$nocluster);
if ($limit) {
	$pipeline->limit($limit);
}
$pipeline->yes($yes);
if ($pipeline->bipd){
	push(@{$predef_steps->{diag_capture}},"gvcf","callable_regions");
	push(@{$predef_steps->{diag_pcr}},"gvcf","callable_regions");
	push(@{$predef_steps->{diag_calling}},"gvcf","callable_regions");
	push(@{$predef_steps->{diag_primer}},"gvcf","callable_regions");
	push(@{$predef_steps->{calling_merge}},"gvcf","callable_regions");
}
my $steps = {
				"alignment"=>  sub {$pipeline->alignment(@_)},
				"fsgs"=>  sub {$pipeline->alignment_fsgs(@_)},
				"readgroup"=> sub {$pipeline->read_group_illumina(@_)},
				"move_bam"=> sub {$pipeline->move_bam(@_)},
				"merge_bam"=> sub {$pipeline->merge_bam(@_)},
				"bam_sort"=> sub {$pipeline->bam_sort(@_)},
				"bamindex" =>sub {$pipeline->bamindex_samtools(@_)},
				"rmdup"=> sub {$pipeline->rmdup(@_)},
				"realign" =>sub {$pipeline->realign_gatk(@_)},
				"covariate" => sub {$pipeline->covariate_illumina(@_)},
				"recalibration_table" => sub {$pipeline->recalibration_table(@_)},
				"coverage" => sub {$pipeline->coverage_samtools(@_)},
				"bcf" => sub {$pipeline->fast_bcf(@_)},
				"purge_fastq" => sub {$pipeline->purge_fastq(@_)},
				"add_chr"=> sub {$pipeline->add_chr(@_)},
				"start_again" => sub {$pipeline->start_again(@_)},
				"chrM_chrMT"=> sub {$pipeline->chrM_chrMT(@_)},
				"reorder_picard"=> sub {$pipeline->reorder_picard(@_)},
				#"calling_merge" => sub {$pipeline->calling_merge(@_)},
				"calling_merge_low" => sub {$pipeline->calling_merge_low(@_)},
			#	"move_vcf" => sub {$pipeline->move_vcf(@_)},
				"mask_primer" => sub {$pipeline->mask_primer_start_end(@_)},
				"delete_chrM" => sub {$pipeline->delete_chrM(@_)},
				"stats" => sub {$pipeline->picard_stats(@_)},
				"rehead_bam" => sub {$pipeline->rehead_bam(@_)},
				"merge_bam" => sub {$pipeline->merge_bam(@_)},
				"gvcf" => sub {$pipeline->calling_gvcf(@_)},
				"gvcf4" => sub {$pipeline->calling_gvcf4(@_)},
				"splitntrim"=>  sub {$pipeline->SplitNCigarReads(@_)},
				"callable_region"=>  sub {$pipeline->callable_region(@_)},
				#"calling_panel"=>  sub {$pipeline->calling_panel(@_,low_calling=>$low_calling)},
				"calling_panel"=>  sub {$pipeline->calling_panel(@_,$pad)},
				"depthofcoverage"=>  sub {$pipeline->depthofcoverage(@_)},
				"duplicate_region_calling"=>  sub {$pipeline->calling_generic(@_,method=>"duplicate_region_calling")},
				"realign_recal"=>  sub {$pipeline->realign_recal(@_)},
				"bam_to_fastq"=>  sub {$pipeline->bam_to_fastq(@_)},
				"zip_fastq"=>  sub {$pipeline->zip_fastq(@_)},
				"rnaseq_metrics"=>  sub {$pipeline->rnaseq_metrics(@_)},
				"callable_regions"=>  sub {$pipeline->callable_regions(@_)},
				"dude"=>  sub {$pipeline->dude(@_)},
				"rmdup_nudup"=>  sub {$pipeline->rmdup_nudup(@_)},
				"binary_depth" =>  sub {$pipeline->lmdb_depth(@_)},
				"elprep" =>  sub {$pipeline->elprep(@_)},
				"manta" =>  sub {$pipeline->manta(@_)},
				"canvas" =>  sub {$pipeline->canvas(@_)},
				"breakdancer" =>  sub {$pipeline->breakdancer(@_)},
				"lumpy" =>  sub {$pipeline->lumpy(@_)},
				"change_chrMT" =>  sub {$pipeline->change_chrMT(@_)},
				"replace_bam" =>  sub {$pipeline->replace_bam(@_)},
				"wisecondor" =>  sub {$pipeline->wisecondor(@_)},
				"calling_wisecondor" =>  sub {$pipeline->calling_wisecondor(@_)},
				"cellranger" =>  sub {$pipeline->cellranger(@_)},
				"annotate_with_umi" =>  sub {$pipeline->annotate_with_umi(@_)},
				"run_alignment_umi" =>  sub {$pipeline->run_alignment_umi(@_)},
				"merge_bam_ubam"  =>  sub {$pipeline->merge_bam_ubam(@_)},
				"bam_to_fastq_umi"  =>  sub {$pipeline->bam_to_fastq_umi(@_)},
				"group_reads_by_umi"  =>  sub {$pipeline-> group_reads_by_umi(@_)},
				"call_consensus_reads"  =>  sub {$pipeline-> call_consensus_reads(@_)},
				"filter_consensus_read"  =>  sub {$pipeline-> filter_consensus_read(@_)},
				"run_alignment_consensus"  =>  sub {$pipeline-> run_alignment_consensus(@_)},
				"split_bam"  =>  sub {$pipeline-> split_bam(@_)},
				"concat_fastq_umi" => sub {$pipeline ->concat_fastq_umi(@_)},
				"fastq_to_bam" => sub {$pipeline ->fastq_to_bam(@_)},
				"fastqScreen" => sub {$pipeline ->fastqScreen(@_)},
				#calling
				"dude" =>  sub {$pipeline->dude(@_)},
				"genotype_gvcf4" =>  sub {$pipeline->genotype_gvcf4(@_)},
				"correct_vcf4" =>  sub {$pipeline->correct_vcf4(@_)},
				"move_and_split_vcf4" =>  sub {$pipeline->move_and_split_vcf4(@_)},
				"correct_vcf"					=> sub{$pipeline->correct_vcf(@_)},
				"move_vcf"						=> sub{$pipeline->move_vcf(@_)},
				"elprep5_target"  => sub{$pipeline->elprep5_target(@_)},
				#by project
				"transcripts_depth" =>  sub {$pipeline->transcripts_depth(@_)},
				"transcripts_coverage" =>  sub {$pipeline->transcripts_coverage(@_)},
				"transcripts_dude" =>  sub {$pipeline->transcripts_dude(@_)},
				"genes_dude" =>  sub {$pipeline->genes_dude(@_)},
				"test_genotype" =>  sub {$pipeline->test_genotype(@_)},
				"move_vcf_hc4" =>  sub {$pipeline->move_vcf_hc4(@_)},
				"cellranger_count" =>  sub {$pipeline->cellranger_count(@_)},
				"cellranger_count_jules" =>  sub {$pipeline->cellranger_count_jules(@_)},
				#umi patrick
				"generate_ubam_umi" => sub {$pipeline->generate_ubam_umi(@_)},
				"align_bam_combine_ubam_umi" => sub {$pipeline->align_bam_combine_ubam_umi(@_)},
				"merge_split_bam_umi" => sub {$pipeline->merge_split_bam_umi(@_)},
				"consensus_bam_umi"=> sub {$pipeline->consensus_bam_umi(@_)},
				"merge_final_bam" => sub {$pipeline->merge_final_bam(@_)},
				"fastq_to_bam" => sub {$pipeline ->fastq_to_bam(@_)},
				"agent_trimmer" => sub {$pipeline->agent_trimmer(@_)},
				"mv_agent_files" => sub {$pipeline->mv_agent_files(@_)},
				"agent_locatit_duplex" => sub {$pipeline->agent_locatit_duplex(@_)},
				"agent_locatit_single" => sub {$pipeline->agent_locatit_single(@_)},
				"agent_locatit_hybrid" => sub {$pipeline->agent_locatit_duplex(@_)},
				"run_alignment_xths2" => sub {$pipeline->run_alignment_xths2(@_)},
				#rnaseq NEB
				"flexbar" => sub {$pipeline->flexbar(@_)},
				"run_alignment_flexbar" => sub {$pipeline->run_alignment_flexbar(@_)},
				"bazam" => sub {$pipeline->bazam(@_)},
				#insertion villarese
				"htlv1_insertion" => sub {$pipeline->htlv1_insertion(@_)},
				#deepseq umi bichat
				"fastq_to_sam" => sub {$pipeline->fastq_to_sam(@_)},
				"sam_to_fastq" => sub {$pipeline->sam_to_fastq(@_)},
				"extract_umi_from_bam" => sub {$pipeline->extract_umi_from_bam(@_)},
				"fastp" => sub {$pipeline->fastp(@_)},
				"run_alignment_deepseq" => sub {$pipeline->run_alignment_deepseq(@_)},
				"filter_proper_pairs" => sub {$pipeline->filter_proper_pairs(@_)},
				"merge_bam_ubam" => sub {$pipeline->merge_ubam(@_)},
				#umi tools
				"nudup" => sub {$pipeline->nudup(@_)},
				"pipeline_genome" => sub {$pipeline->pipeline_genome(@_)},
				"cnvnator" => sub {$pipeline->cnvnator(@_)},
				"melt" => sub {$pipeline->melt(@_)},
				"sortdedup" => sub {$pipeline->sortdedup(@_)},
				"bwa2" => sub {$pipeline->bwa2(@_)},
				"elprep5_genome" => sub {$pipeline->elprep5_genome(@_)},
				"muc1" => sub {$pipeline->muc1(@_)},
				"advntr" => sub {$pipeline->advntr(@_)},
				"star_align" => sub {$pipeline->star_align(@_)},
				"deepvariant" => sub {$pipeline->deepvariant(@_)},
				"rnaseqsea_capture" => sub {$pipeline->rnaseqsea_capture(@_)},
				"hificnv" => sub {$pipeline->hificnv(@_)},	
				"rnaseqsea_rnaseq" => sub {$pipeline->rnaseqsea_rnaseq(@_)},
				"specie_contaminant_check" => sub {$pipeline->check_specie_contaminant(@_)},

			};
			
my @types_steps = ('pipeline','calling');
my $list_steps;
my $list_steps_types;
if ($arg_steps){
	my $type = "pipeline";
	my $x;
	foreach my $n (split(",",$arg_steps)){
		die() unless exists $steps->{$n};
		push (@$x,$n);
		
	}
	 push(@$list_steps,$x);
	 push(@$list_steps_types,$type);
	
}
else {
foreach  my $type (@types_steps){
	my @list = sort {$a cmp $b} keys %{$define_steps->{$type}};
	if ($pipeline_name){
		 push(@$list_steps,[split(",",$define_steps->{$type}->{$pipeline_name})]);
   			push(@$list_steps_types,$type);
   			next;
	}
	push(@list,'none');
	my $x =  colored ['black ON_BRIGHT_GREEN'];
	my $banner=" $x Please Pick a $type  for $projectName :";
	
   	my $choice = &pick(\@list,$banner,20);
   	my $text;
   ($steps_name,$text)=split(":",$choice);
   next if $steps_name eq 'none';
   push(@$list_steps,[split(",",$define_steps->{$type}->{$steps_name})]);
   push(@$list_steps_types,$type);
}
}
$pipeline->priority_name($list_steps_types);
$pipeline->queue("-q testq");
# if $project->isDiag;
#$pipeline->queue("") if $secret;
my $dir_bds =$pipeline->dir_bds();
#$pipeline->fastq_extend($fastq_ext) if $fastq_ext;
$pipeline->max_cpu($max_cpu) if $max_cpu ;

$SIG{'INT'} = sub {
	if (defined $pipeline->daemon){
		$pipeline->daemon->Kill_Daemon($pipeline->pid,15) if $pipeline->daemon->Status($pipeline->pid);
		sleep(5);
	 	system("$Bin//scripts/scripts_pipeline/parse_qstat.pl ".$pipeline->timestamp());
	 	sleep(2);
	 	
	}
	 if (defined $pipeline->process) {
	 	$pipeline->process->kill();
	 	sleep(5);
	 	system("$Bin//scripts/scripts_pipeline/parse_qstat.pl ".$pipeline->timestamp());
	 	sleep(2);
	 };
	 $pipeline->clean_error;
	 exit(1);
};


#@running_steps =@$list_steps;



my $patients = file_util::return_patients( $project, $patients_name );
confess("no patients") unless scalar(@$patients);



$pipeline->unforce(0) if $force;


print colored ['black '], "################################################################################################" ;

print "\n";
print  join("->",@running_steps)."\n";
print colored ['black '], "################################################################################################" ;
print "\n";
print "\n";

foreach my $p (@$patients) {
		my $file = $p->trackingFile().".lock";
		unlink $file if -e $file;
		if (-e $p->trackingFile()){
			my $backup_dir = $project->getPipelineTrackingDir()."/backup";
			system ("mkdir $backup_dir ; chmod a+rwx $backup_dir") unless (-d $backup_dir);
			my $log = new Logfile::Rotate( File => $p->trackingFile(),	
		 								Count => 7,
		 								Gzip  => 'lib',
		 								 Flock  => 'no',
		 								 Dir => $backup_dir,
		 								);
			warn "start log";
			$log->rotate();
			unlink $p->trackingFile; 								
		}

my $file2 = $p->trackingFile();
next if -e "$file2.lite";
my $dbh = DBI->connect("dbi:SQLite:dbname=$file2.lite","","");
$dbh->{AutoCommit} = 0;
my $sql_create = qq{CREATE TABLE IF NOT EXISTS `steps` (
  `step_id` INTEGER PRIMARY KEY,
  `run_id` TEXT ,
  `project_id` TEXT ,
  `software` TEXT ,
  `software_version` TEXT ,
  `machine` TEXT ,
  `start` INTEGER ,
  `status` TEXT ,
  `cmd` TEXT ,
  `step_name` TEXT ,
  `end` INTEGER ,
  UNIQUE(run_id,step_name)
);
};
$dbh->do($sql_create) ;
$dbh->commit();
system("chmod a+rw $file2.lite");
		
}

my $priority = 0;
my $n = 0;
my $end_files;
my $nb_type = 0;
warn Dumper $list_steps;

foreach my $list_requests (@{$list_steps}) {
	
	my $type_objects = $list_steps_types->[$nb_type];
	warn $type_objects;
	$nb_type ++;
	if ($type_objects eq 'calling') {
		$pipeline->add_sample_with_priority($project, $priority);
		push(@$end_files,prepare_jobs($list_requests, $steps));
	}
	elsif ($type_objects eq 'chromosomes') {
		foreach my $chr (@{$project->getChromosomes()}) {
			#warn $chr->name();
			$pipeline->add_sample_with_priority($chr, $priority);
			push(@$end_files,prepare_jobs($list_requests, $steps));
		}
	}
	elsif ($type_objects eq 'pipeline') {
		foreach my $patient (@{$patients}) {
			$pipeline->add_sample_with_priority($patient, $priority);
			push(@$end_files,prepare_jobs($list_requests, $steps));
		}
	}
	elsif ($type_objects eq 'transcripts') {
		foreach my $patient (@{$patients}) {
			$pipeline->add_sample_with_priority($patient, $priority);
			push(@$end_files,prepare_jobs($list_requests, $steps));
		}
	}
	$priority++;
	$n++; 
}
$pipeline->samples();
$pipeline->print_all_steps_by_prority;

print "\n";
print colored ['black '], "################################################################################################" ;
			print "\n";

print colored ['black '], "################################################################################################" ;
print "\n";

unless ($yes){
	my $choice = prompt("run this/these step(s)   (y/n) ? ");
	die() if ($choice ne "y"); 
}

$pipeline->clean();



#$pipeline->list_cmd();

$pipeline->launch_bds_daemon_by_priority();

if ($pipeline->error > 0){
	die();
}


exit(0);

sub prepare_jobs {
	my ($running_steps,$steps) = @_;
	my $next_file = "";
	foreach my $step (@$running_steps){
		($next_file) = $steps->{$step}->({filein=>$next_file});
	}
	return $next_file;
}



sub DESTROY {
      my $self = shift;
     if (defined $pipeline->daemon){
		$pipeline->daemon->Kill_Daemon($pipeline->pid,15) if $pipeline->daemon->Status($pipeline->pid);
		sleep(5);
	 	system("$Bin//scripts/scripts_pipeline/parse_qstat.pl ".$pipeline->timestamp());
	 	sleep(2);
	 	
	}
	 if (defined $pipeline->process) {
	 	$pipeline->process->kill();
	 	sleep(5);
	 	system("$Bin//scripts/scripts_pipeline/parse_qstat.pl ".$pipeline->timestamp());
	 	sleep(2);
	 };
	 $pipeline->clean_error;
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

