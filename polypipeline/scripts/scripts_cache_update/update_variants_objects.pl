#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use Parallel::ForkManager;
use Storable qw(store retrieve freeze thaw);
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use Cwd 'abs_path';
use Digest::MD5 qw(md5_hex);
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use Set::IntSpan::Fast::XS;
use Sys::Hostname;

my $host = hostname();
#warn "*_*_*_*_*_ ".$host." _*_*_*_*_*_";

my $fork = 1;
my $force;
my ($project_name, $chr_name);
GetOptions(
	'fork=s'    => \$fork,
	'project=s' => \$project_name,
	'chr=s'     => \$chr_name,
	'force=s' 	=> \$force,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }

my @methods_to_delete = (
	"annotation",
);

my @methods_to_run = (
	"gnomad",
	"dejaVuInfosForDiag",
	"name",
	"cosmic",
	"getGenes",
	"getTranscripts",
	"score_clinical_local",
	"score_clinvar",
	"text_clinvar",
	"comment_clinical_local",
	"cadd_score",
	"dbscsnv_ada",
	"dbscsnv_rf",
	"revel_score",
	"hgmd_id",
#	"get_codon_text(1)",
#	"get_codon_text(-1)",
	"annotation",
	"getSequence",
	"ref_allele",
	"in_this_run_ratio",
	"ncboost_score",
	"spliceAI",
);

my $nbErrors = 0;
my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name => $project_name, -cache => '1', -typeFilters => 'familial' );
my $chromosomes = $project->getChromosomes;

my $log_file = $project->lmdb_cache_variations_dir()."/".$chr_name.'.update';
if (-e $log_file) {
	my $cmd = 'rm '.$log_file;
	`$cmd`;
}

foreach my $chr (@{$chromosomes}){
	if ($chr_name){
		next if $chr->name ne "$chr_name";
	}
	run_update($chr->name);
}
exit(0) if ($chr_name);


sub run_update {
	my ($chr_name) = @_;
	my $project = $buffer->newProjectCache( -name => $project_name, -cache => '1', -typeFilters => 'familial' );
	my $chr = $project->getChromosome($chr_name);
	my $no = $chr->lmdb_variations("r");

#	my $no3 = $chr->lmdb_score_impact("c");
#	$no3->put("date",time);
#	$no3->close;
	exit(0) unless $no->exists_db;

	my $ranges = $no->ranges($fork);
	$no->close();
	my $pm = new Parallel::ForkManager($fork);
	my $dir_out = $project->lmdb_cache_variations_dir()."/test/";
	$dir_out = "/tmp/";
	my $f1 = $dir_out."/".$chr_name;
	unlink $f1."_index" if -e $f1."_index";
	unlink $f1 if -e $f1;

	$pm->run_on_finish(
	    sub { 
	    	my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$h)=@_;
	    	unless (defined($h) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
		 	my $no  = GenBoNoSqlLmdb->new(dir=>$dir_out,mode=>"w",is_index=>1,name=>$chr_name,is_compress=>1);
		 	foreach my $o (@{$h->{obj}}){
		 		$no->put_with_index($o->{id},$o);
		 	}
		 	$no->close();
			$no = undef;
	    }
	);
   
	$project->getPatients;
	$project->buffer->dbh_deconnect();
	foreach my $r (@$ranges){
		my $pid = $pm->start and next;
		my $nb =0;
		
		$project->buffer->dbh_reconnect();
		my $vquery = validationQuery->new(dbh=>$buffer->dbh,capture_name=>"idefix");
		my $no = $chr->lmdb_variations("r");
		my $h;
		for (my $i=$r->[0];$i<$r->[1]+1;$i++){
			$nb ++;
#			warn $nb if $nb%10000 == 0;
			my $v = $no->get_index($i);
			$v->{project} =  $project;
			$v->{buffer} = $buffer;
			#purge data
			if ($force){
				$v->purge_deja_vu();
				$v->purge_public_data();
			}
			my $debug ;
			$debug =1 if $v->id eq "21_44589921_G_GC";
			
			my $chr_name = $v->getChromosome()->name;
			
			foreach my $cat (@methods_to_delete) {
				delete $v->{$cat};
			}
			foreach my $method_name (@methods_to_run) {
				$v->$method_name;
			}
			
		  	foreach my $tr ( @{$v->getTranscripts}){
				$v->getNomenclature($tr);
				foreach my $p (@{$tr->getProteins}){
					$v->polyphenScore($p);
					$v->siftScore($p);
		  	 	}
		  	}
			
			delete $v->{project} ;
			delete $v->{buffer};
			push(@{$h->{obj}},$v);
		}
		$no->close();
		$no = undef;
		$pm->finish(0,$h);
		die();
	}
	$pm->wait_all_children();
	$project->buffer->dbh_reconnect();
	
#	warn "end variants";
	my $no1 = $chr->lmdb_variations("r");
	
	my $no2 = GenBoNoSqlLmdb->new(dir=>$dir_out,mode=>"r",is_index=>1,name=>$chr_name,is_compress=>1);
	foreach my $r (@$ranges){
		for (my $i=$r->[0];$i<$r->[1]+1;$i++){
			my $v = $no1->get_index($i);
			my $v2  = $no2->get_index($i);
			die()  if $v->{id} ne $v2->{id};
		}
	}
	$no1->close();
	$no2->close();
	
	my $chr = $project->getChromosome($chr_name);
	my $f1 = $dir_out."/".$chr_name;
	my $f2 = $project->lmdb_cache_variations_dir()."/".$chr_name;
	system("mv $f1 $f2 ");
	die() unless -e $f2;
	$f1 = $dir_out."/".$chr_name."_index";
	$f2 = $project->lmdb_cache_variations_dir()."/".$chr_name."_index";
	system("mv $f1 $f2 ");
	
	open (LOG, ">>".$log_file);
	print LOG "\n";
	print LOG "# Methods deleted\n";
	foreach my $cat (@methods_to_delete) {
		print LOG $cat."\n";
	}
	print LOG "\n";
	print LOG "# Methods launch\n";
	foreach my $method_name (@methods_to_run) {
		print LOG $method_name."\n";
	}
	print LOG "Polyphen\n";
	print LOG "Sift\n";
	close(LOG);
	return 1;
}