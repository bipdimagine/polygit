#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo/lib/obj-nodb/";
use lib "$Bin/packages";
use GBuffer;
#use GenBoProject;
#use PBS::Client;
use Getopt::Long;
use Data::Dumper;
use IO::Prompt;
use Sys::Hostname;
use Parallel::ForkManager;
use Term::ANSIColor;
use Moo;

#use bds_cache_steps;    
use bds_cache_rocks; 
#use file_util;
use Class::Inspector;
#use check_utils;
use Text::Table;
use colored; 
use Term::Menus;
use Config::Std;
 
my ($projectName, $filename, $name, $patients_name, $steps_name, $force, $type, $fastq_ext, $somatic, $method, $no_cluster, $stdout, $help, $annot_version);
my $fork = 1;
my $yes = 0;
my $nocluster = 0;
my $menu= 0;
my $filename_cfg;
my $secret;
my $analyse_type;
my $nobackup;
my $giab;
my $clean;

GetOptions(
	'project=s' => \$projectName,
	'patients=s' => \$patients_name,
	'steps=s' => \$steps_name,
	'force=s' => \$force,
	'yes=s' => \$yes,
	'annot_version=s' => \$annot_version,
	'nocluster=s' => \$nocluster,
	'menu=s' => \$menu,
	'nolimit=s' => \$secret,
	'config=s' => \$filename_cfg,
	'type=s' => \$analyse_type,
	'h!' => \$help,
	"nobackup=s" => \$nobackup,
	"control=s" => \$giab,
	"clean=s" => \$clean,
);

#$steps_name = "all" unless $steps_name;

check_version($projectName) unless ($steps_name eq 'update' );
clean($projectName) if $clean;
my $report;
my $buffer = GBuffer->new();
my $define_steps;
unless ($filename_cfg) {
 $filename_cfg = $buffer->config_directory()."/pipeline/cache.cfg";
}

read_config $filename_cfg =>  %{$define_steps};

my $project = $buffer->newProject( -name => $projectName );

if ($annot_version) {
	$project->changeAnnotationVersion($annot_version);
}

unless ($steps_name){
	if ($project->isDiagnostic()){
		$steps_name = "diag";
	}
	else {
		$steps_name = "query";
	}
}



#$predef_steps->{update} = [["update_variants","update_chromosomes","update_genes","polyviewer"]];
#$predef_steps->{update_type} = ["chromosomes"];
my $pipeline;
 $pipeline = bds_cache_rocks->new( project=>$project, nocluster=>$nocluster,cache=>1 );
$pipeline->queue("-q testq") unless $secret;
my $steps = {
				"store_ids"=>  sub {$pipeline->store_ids(@_)},
				"store_annotations"=> sub {$pipeline->store_annotations(@_)},
				"strict_denovo"=> sub {$pipeline->strict_denovo(@_)},
				"loh"=> sub {$pipeline->loh(@_)},
				"global_infos"=> sub {$pipeline->global_infos(@_)},
				"coverage" =>sub {$pipeline->coverage(@_)},
				"dejavu" =>sub {$pipeline->dejavu(@_)},
				"check_store_ids"=>  sub {$pipeline->check_store_ids(@_)},
				"check_store_annotations"=> sub {$pipeline->check_store_annotations(@_)},
				"check_strict_denovo"=> sub {$pipeline->check_strict_denovo(@_)},
				"check_loh"=> sub {$pipeline->check_loh(@_)},
				"check_global_infos"=> sub {$pipeline->check_loh(@_)},
				"check_coverage" =>sub {$pipeline->check_coverage(@_)},
				"polydiag" =>sub {$pipeline->polydiag(@_)},
				"quality_check" =>sub {$pipeline->quality_check(@_)},
				"update_annotation" =>sub {$pipeline->update_annotation(@_)},
				"tree_cache" =>sub {$pipeline->tree_cache(@_)},
				"update_score" =>sub {$pipeline->update_score(@_)},
				"update_coverage" =>sub {$pipeline->update_coverage(@_)},
				"diagHash"=>  sub {$pipeline->diagHash(@_)},
				"polyviewer"=>  sub {$pipeline->polyviewer(@_)},
				"polydiag" =>sub {$pipeline->polydiag(@_)},
				"dude_bed" =>sub {$pipeline->dude_bed(@_)},
				"transcripts_dude" =>sub {$pipeline->transcripts_dude(@_)},
				"genes_dude" =>sub {$pipeline->genes_dude(@_)},
				"update_variants" =>sub {$pipeline->update_variants(@_)},
				"update_chromosomes" =>sub {$pipeline->update_chromosomes(@_)},
				"update_genes" =>sub {$pipeline->update_genes(@_)},
				"cnv" =>sub {$pipeline->cnv_manue(@_)},
				"identito_vigilence" =>sub {$pipeline->identito_vigilence(@_)}, 
				"cache_html_polyviewer" =>sub {$pipeline->html_cache_polyviewer(@_)}, 
				"cache_html_polycyto" =>sub {$pipeline->html_cache_polycyto(@_)}, 
				"polydude" =>sub {$pipeline->polydude(@_)}, 
				"sashimi_plots"=>sub {$pipeline->sashimi_plots(@_)}, 
				"store_rna_junction_ids"=>sub {$pipeline->store_rna_junction_ids(@_)}, 
				"merge_patients" => sub {$pipeline->merge_patients(@_)}, 
				"merge_objects" => sub {$pipeline->merge_objects(@_)}, 
			};



if ($help) { usage(); }
unless ($steps_name) { usage(); }

$pipeline->yes($yes);
$pipeline->unforce(0) if $force;
#$pipeline->add_sample(patient=>$project);

my @types_steps = ('dude','chromosomes','project','polydiag','html_cache');
 @types_steps = ('chromosomes','rocks','html_cache','project') if $project->isGenome or $project->isExome ;
@types_steps = ('chromosomes') if $giab ;
my $list_steps;
my $list_steps_types;

unless ($menu) {
	foreach  my $type (@types_steps){

		my $ht = "exome";
		if ($project->isDiagnostic){
			$ht = "diag";
		}
		if ($project->isGenome){
			$ht = "genome";
		}
		if ($analyse_type){
			$ht = "$analyse_type";
		}

	
		next unless $define_steps->{$type}->{$ht};

		push(@$list_steps,[split(",",$define_steps->{$type}->{$ht})]);
   		push(@$list_steps_types,$type);
	}
}
else {
foreach  my $type (@types_steps){
	my @list = sort {$a cmp $b} keys %{$define_steps->{$type}};
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

#warn Dumper $list_steps_types;
$pipeline->priority_name($list_steps_types);





my $bin_dev = qq{$Bin/scripts/scripts_pipeline/};

my $cmd_first = qq{/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/scripts/cache_global_infos.pl -project=$projectName };
system("$cmd_first");

my $end_files;

#confess ("\nERROR: $steps_name doesn't exists. Die.\n") unless (exists $predef_steps->{$steps_name});
my $n = 0;
my $priority = 0;
my $nb_type = 0;
foreach my $list_requests (@{$list_steps}) {
	my $type_objects = $list_steps_types->[$nb_type];
	$nb_type ++;
	$priority++;
	if ($type_objects eq 'dude') {
		foreach my $patient (@{$project->getPatients()}) {
			$pipeline->add_sample_with_priority($patient, $priority);
			push(@$end_files,prepare_calling_jobs($list_requests, $steps));
		}
	}
	elsif ($type_objects eq 'project') {
		$pipeline->add_sample_with_priority($project, $priority);
		push(@$end_files,prepare_calling_jobs($list_requests, $steps));
	}
	elsif ($type_objects eq 'chromosomes') {
		foreach my $chr (@{$project->getChromosomes()}) {
			$pipeline->add_sample_with_priority($chr, $priority, 'chr');
			
			push(@$end_files,prepare_calling_jobs($list_requests, $steps));
		}
	}
	elsif ($type_objects eq 'polydiag' or $type_objects eq 'html_cache' ) {
		foreach my $patient (@{$project->getPatients()}) {
			$pipeline->add_sample_with_priority($patient, $priority);
			push(@$end_files,prepare_calling_jobs($list_requests, $steps));
		}
	}
	$n++; 
}




$pipeline->samples();
$pipeline->print_all_steps_by_prority;

unless ($yes){
	my $choice = prompt("run this/these step(s)   (y/n) ? ");
	die() if ($choice ne "y"); 
}
$pipeline->clean();
$buffer->getQuery()->insertHistoryCacheVersion($project->id,$project->annotation_version);
if ($pipeline->nocluster ne 2){
	$SIG{'INT'} = sub {
		 if (defined $pipeline->daemon){
			$pipeline->daemon->Kill_Daemon($pipeline->pid,15) if $pipeline->daemon->Status($pipeline->pid);
			sleep(2);
		 	system("$Bin//scripts/scripts_pipeline/parse_qstat.pl ".$pipeline->timestamp());
		 	sleep(2);
		}
		elsif (defined $pipeline->process()) {
			warn "process";
		 	$pipeline->process->kill();
		 	sleep(5);
		 	system("$Bin//scripts/scripts_pipeline/parse_qstat.pl ".$pipeline->timestamp());
		 };
		 $pipeline->clean_error;
		 exit(1);
	};
}

$pipeline->launch_bds_daemon_by_priority();

exit(0);



##### METHODS ######

sub prepare_calling_jobs {
	my ($running_steps,$steps) = @_;
	
	my $next_file = "";
	foreach my $step (@$running_steps){
		($next_file) = $steps->{$step}->({filein=>$next_file});
	}
	return $next_file;
}



sub usage {
	print colored ['red'],  "\n======================= USAGE ========================\n";
	print colored ['blue'], $0." -project=project_name -steps=(check or polyquery or polydiag or all) -force=(1 : force restart step) -cpu_max = (available CPU number if run without PBS protocol) \n";
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

sub clean {
		my ($project_name) = @_;
		my $buffer = GBuffer->new();
		my $project = $buffer->newProject( -name => $project_name );
		my $cache_directory_actual = $project->getCacheDir();
		my $choice = prompt(colored ['black ON_BRIGHT_YELLOW'],"delete  $cache_directory_actual (y/n) ? ") unless $yes;;
		die() if ($choice ne "y"); 
		$choice = prompt(colored ['black ON_BRIGHT_CYAN'],"no regret :  $cache_directory_actual (y/n) ? ") unless $yes;
		die() if ($choice ne "y");
		colored::stabilo("cyan","As you wish ....");
		system ("rm -r $cache_directory_actual/*");
		system("rmdir $cache_directory_actual");
		
}


sub check_version {
	my ($project_name) = @_;
	my $tb = Text::Table->new( (" ", colored::stabilo("magenta","Actual", 1), colored::stabilo("magenta","Latest", 1),, colored::stabilo("magenta","Status", 1) ) );
	
	my $hchange ={};
	my $buffer = GBuffer->new();
	my $project = $buffer->newProject( -name => $project_name );
	my $query =  $buffer->getQuery();
	my $gencode_actual = $project->gencode_version();
	my $change;
	my $gencode_available = $query->getMaxGencodeVersion($project->id);;
	my $cache_directory_actual = $project->getCacheDir();
	colored::stabilo("magenta","======================= Check Version   =========================");
	#colored::stabilo("red","Too bad you are not eligible to a genCode update.") unless $project->is_gencode_upgradable();
	my @lines;
	push(@{$lines[0]},"GENCODE");
	push(@{$lines[0]},colored::stabilo("cyan","$gencode_actual", 1));
	push(@{$lines[0]},colored::stabilo("white","N/A", 1)) unless $project->is_gencode_upgradable();
	if ($project->is_gencode_upgradable()){
	if ($gencode_actual < $gencode_available) {
		push(@{$lines[0]},colored::stabilo("yellow","$gencode_available", 1)); 
		push(@{$lines[0]},colored::stabilo("red","Update", 1)); 
	#	colored::stabilo("yellow","You are Lucky a new GENCODE version is avaliable ");
	#	colored::stabilo("yellow","Project $project_name GENCODE actual version:  * $gencode_actual *   vs latest version : * $gencode_available *");
		$hchange->{gencode} ++;
	
	}
	else {
		
		push(@{$lines[0]},colored::stabilo("green","$gencode_available", 1)); 
		push(@{$lines[0]},colored::stabilo("green","OK", 1)); 
	#	colored::stabilo("green","Project $project_name GENCODE actual version: $gencode_actual  vs latest version : $gencode_available");
		print  "\n ";
	}
	}
	
#	colored::stabilo("magenta","======================= Public Database  =========================");
	
	my $version_actual = $project->public_database_version();
	my $version_available = $query->getMaxPublicDatabaseVersion($project->id);
	push(@{$lines[1]},"Public Database");
	push(@{$lines[1]},colored::stabilo("cyan","$version_actual", 1));
	
	if ( $version_actual < $version_available) {
		push(@{$lines[1]},colored::stabilo("yellow","$version_available", 1));
		push(@{$lines[1]},colored::stabilo("red","Update", 1)); 
		$hchange->{public_database} ++;
	}
	else {
		push(@{$lines[1]},colored::stabilo("green","$version_available", 1)); 
		push(@{$lines[1]},colored::stabilo("green","OK", 1)); 
	}
		$tb->load(@lines);
		print $tb;
	
	if (exists $hchange->{gencode}){
          my $choice = "y";
			 $choice = prompt("\nDo you want to upgrade your GENCODE Version   (y/n) ? ") unless $yes;
			if ($choice eq "y"){
				$change =1;
			} 
			else {
				delete $hchange->{gencode};
			}
	}

	if (exists $hchange->{public_database}){
			print "\n-------------------************----------------\n";
          	my $choice = "y";
			 $choice = prompt("\nDo you want to upgrade your PUBLIC  RELEASE   (y/n) ? ") unless $yes;
			if ($choice eq "y"){
				$change = 1;
			} 
			else {
				delete $hchange->{public_database};
			}
	}
	
	
	
	if (keys %$hchange){
		if ($project->isGenome or $project->isExome){
			my $choice = "y";
			  $choice = prompt("\nDo you want to backup  your cache  (y/n) ? ") unless $yes;
			if ($choice ne "y"){
				$nobackup = 1;
			} 
		}
		colored::stabilo("cyan","please wait , archiving old cache  : $nobackup");
		my $archives_dir = $buffer->getDataDirectory("archives");
		my $archives_name = $archives_dir."/".$project->name.".".$project->genome_version.".$gencode_actual.$version_actual".".tar";
		unless ($nobackup){
		system("tar -cvf $archives_name $cache_directory_actual >/dev/null && pigz -p 10 $archives_name");
		unless (-e $archives_name.".gz"){
			colored::stabilo("red","Opps I'm not able to backup your files "."$archives_name");
			die();
		}
		}
		print colored ['green'],"\n backup ok ... "."\n";
		my $choice = "y";
		 $choice = prompt(colored ['black ON_BRIGHT_YELLOW'],"delete  $cache_directory_actual (y/n) ? ") unless $yes;;
		die() if ($choice ne "y"); 
		$choice = prompt(colored ['black ON_BRIGHT_CYAN'],"no regret :  $cache_directory_actual (y/n) ? ") unless $yes;
		die() if ($choice ne "y");
		colored::stabilo("cyan","As you wish ....");
		system ("rm -r $cache_directory_actual/*");
		system("rmdir $cache_directory_actual");
		colored::stabilo("yellow"," proceed to database "."\n");
		if (exists $hchange->{gencode}){
			$project->update_gencode_version();
			delete $project->{gencode_version};
			print colored ['red'],"\nHmmm .... it's look like your upgrade is OK ...\n".$project->gencode_version();
		}
		if (exists $hchange->{public_database}){
			
			$project->update_public_database_version();
			delete $project->{public_database_version};
			colored::stabilo("green","\nAnd now the current Public Release is  :".$project->public_database_version()."\n");
		}
	}
	print "\n";
	colored::stabilo("green"," LET'S GO   !!!! " );
	print "Press <Enter> or <Return> to continue: ";
	my $resp = <STDIN> unless $yes;
	$buffer = undef;
	$project = undef;
	
}