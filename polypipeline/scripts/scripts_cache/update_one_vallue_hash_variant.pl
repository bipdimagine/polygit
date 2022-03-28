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
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
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
use File::Temp qw/ tempfile tempdir /;
#use Mojo::DOM58; 

require "$RealBin/../../../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
 my $host = hostname();


warn "*_*_*_*_*_".$host."*_*_*_*_*_";
#use Cache_Commons;

my $fork = 1;
my $force;
my ($project_name, $chr_name, $annot_version);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'chr_name=s'        => \$chr_name,
	'force=s'  => \$force,
);




unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
#unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }
my $t = `hostname`;
warn $t;
my $nbErrors = 0;
my $buffer = GBuffer->new();
my $project = $buffer->newProjectCache( -name => $project_name );

foreach my $chr (@{$project->getChromosomes} ){
	my $tt =time;
	warn $chr->name. '----';
	run_update_chr($chr->name);
	warn abs(time - $tt);
}

exit(0);

sub run_update_chr {
	my ($chr_name) = @_;
my $buffer = GBuffer->new();
my $project = $buffer->newProjectCache( -name => $project_name );
warn '**__START '.$chr_name;
my $chr = $project->getChromosome($chr_name);

my $pm = new Parallel::ForkManager($fork);

my $dir_tmp = tempdir( CLEANUP => 1,DIR => "/tmp/" );

$pm->run_on_finish(
    sub { 
    	my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$h)=@_;
  
    	unless (defined($h) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
	
		return;
    }
    
    );
    
$project->getPatients;
$project->buffer->dbh_deconnect();
my $list =listVariants($chr);
warn scalar(@$list);
my $nb = int(scalar(@$list)/($fork) +0.5);
warn $nb;
my $iter = natatime($nb, @$list);

my $part =0;
my @files ;
	my $javascript_id = time + int(rand(10000));
	warn $a;
  while( my @tmp = $iter->() ){
	$part ++;
	warn $part;
	my $f = $chr->name."_".$part;
	push(@files,$f);
	my $pid = $pm->start and next;
	my $nb =0;
	$project->buffer->dbh_reconnect();
	my $array;
	$project->setListVariants(\@tmp);
	my $max =scalar(@tmp);

	my $no  = $chr->lmdb_hash_variants("r");
	my $nb =0;
	while(my $v = $project->nextVariant){
		$v->spliceAI();
		$nb ++;
		warn $nb if $nb %10000 ==0;
		#die() if $nb %100000 ==0;
		foreach my $patient (@{$v->getPatients}) {
			my $z = $no->get($v->id."@".$patient->name);
			#warn Dumper keys %$z;
		}
		
		
	}#end for ecah patient
	$no->close();
	$pm->finish(0,{});
	
	#die();
}
warn 'wait';
$pm->wait_all_children();
$project->buffer->dbh_reconnect();
#warn "create ";
#foreach my $f (@files){
#	unlink $dir_tmp."/".$f;
#	unlink $dir_tmp."/".$f."_index";
#}
#die();

return 1;
}


sub listVariants {
		my($chr) =@_;
		my $vector = $chr->getVectorScore("gnomad_ho_ac_all");
		my $vector = $chr->getVectorVariations;#("gnomad_ho_ac_all");
		 $vector += $chr->getVectorScore("gnomad_ho_ac_50");
		my $v1 = $chr->vectorDM();
	 	$v1  |=  $chr->vectorClinvarPathogenic();
	
		 $vector  |= $v1;
		
		my @list_variants;
		my $already;
	
		foreach my $id (@{to_array($vector,$chr->name)}) {
				push (@list_variants,$id);
			}
		return( \@list_variants);
}
sub to_array {
	my ($v,$name) = @_;
	my $set = Set::IntSpan::Fast::XS->new($v->to_Enum);
	my $iter = $set->iterate_runs();
	my @t;
	while (my ( $from, $to ) = $iter->()) {
   		for my $member ($from .. $to) {
   			push(@t,$name."!".$member);
   		}
    }
    return \@t;
}
