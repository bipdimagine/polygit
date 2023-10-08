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
use lib "$RealBin/../../../GenBo/lib/obj-nodb/packages";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use Set::IntSpan::Fast::XS;
use GBuffer;

#use Cache_Commons;

my $fork = 1;
my $force;
my ($project_name, $chr_name,$annot_version);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'chr=s'        => \$chr_name,
	'annot_version=s'    => \$annot_version,
	'force=s'  => \$force,
);
my $hrun;
unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
#unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }
my $t = `hostname`;

my $nbErrors = 0;
my $buffer = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProjectCache( -name => $project_name, -cache => '1', -typeFilters => 'familial' );
if ($annot_version) {
	$project->changeAnnotationVersion($annot_version);
}

my $scaled_score_ncboost =  $project->buffer->config->{scaled_score_ncboost};
my $scaled_score_frequence_public =  $project->buffer->config->{scaled_score_frequence_public};
my $scaled_score_gnomad_ho =  $project->buffer->config->{scaled_score_gnomad_ho_ac};
my $scaled_score_gnomad_ac =  $project->buffer->config->{scaled_score_gnomad_ac};
my $scaled_score_ratio =  $project->buffer->config->{scaled_score_ratio};
my $scaled_score_project_dejavu =  $project->buffer->config->{project_dejavu};
my $scaled_score_sample_dejavu =  $project->buffer->config->{sample_dejavu};
my $scaled_score_sample_ho_dejavu =  $project->buffer->config->{sample_ho_dejavu};
my $functional_categorie = $project->buffer->config->{functional_annotations};
my $scaled_score_ncboost =  $project->buffer->config->{scaled_score_ncboost};

 my $pathogenic_categories = ["hgmd","clinvar_pathogenic","dm"];

my @all_categories;
foreach my $cat (keys %$scaled_score_ncboost) {
	push(@all_categories,'ncboost_'.$cat);
}
push(@all_categories,keys %$scaled_score_frequence_public);
push(@all_categories,keys %$scaled_score_gnomad_ac);
push(@all_categories,keys %$scaled_score_sample_ho_dejavu);

push(@all_categories,keys %$scaled_score_project_dejavu);
push(@all_categories,keys %$scaled_score_sample_dejavu);
push(@all_categories,keys %$scaled_score_gnomad_ho);

foreach my $p (@{$project->getPatients}){
	foreach my $k (keys %$scaled_score_ratio) {
		push(@all_categories,$p->name."+".$k);
	}
}

push(@all_categories,keys %$functional_categorie);


push(@all_categories,@$pathogenic_categories);
push(@all_categories,"genes");
push(@all_categories,"intergenic");
push(@all_categories,"intronic");
push(@all_categories,"acmg");

my $final_vector_categories_chromosomes ={};

my $chr = $project->getChromosome($chr_name);
my $empty_vector = $chr->getNewVector;	
	
my $list;
my $hv;
my $no      = $chr->get_rocks_variations("r");
 my $ranges = $no->ranges($fork);
$no->close;
delete $chr->{rocks};
#foreach my $chr (@{$project->getChromosomes}){


my $pm = new Parallel::ForkManager($fork);



$pm->run_on_finish(
    sub { 
    	my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$h)=@_;
  
    	unless (defined($h) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
	
		delete $hrun->{$h->{run_id}};
		my $hintspan = $h->{vector_categories_chromosomes};
		
		foreach my $chr (keys %{$h->{vector_categories_chromosomes}}){
			foreach my $c (keys %{$h->{vector_categories_chromosomes}->{$chr}}){
					#warn $c."::".$h->{vector_categories_chromosomes}->{$chr}->{$c}->to_Enum;
					$final_vector_categories_chromosomes->{$chr}->{$c} = $empty_vector->Clone unless exists $final_vector_categories_chromosomes->{$chr}->{$c};
					$final_vector_categories_chromosomes->{$chr}->{$c} |= $h->{vector_categories_chromosomes}->{$chr}->{$c};
			}
		}
    }
    
    );



	my $id = time;
$project->preload_patients();
	$project->buffer->disconnect();
	
		#$project->buffer->{dbh} ="-";
		
	foreach my $r (@$ranges){
		$id ++;
		$hrun->{$id} ++;
		my $pid = $pm->start and next;
		my $res = construct_vector($project,$r);
		$res->{run_id} = $id;
		$pm->finish(0,$res);
	}
$pm->wait_all_children();
die() if keys %$hrun;



 $chr = $project->getChromosome($chr_name);
	#my $no3 = $chr->lmdb_score_impact("c");
	my $rocks3 = $chr->rocks_vector("w");
	
	foreach my $c (@all_categories){
		unless (exists $final_vector_categories_chromosomes->{$chr->name}->{$c}){
			#warn "vide $c " if $c eq "intergenic";
			
			$rocks3->put_batch($c,$empty_vector->Clone);
		}
		else{
		#	warn $final_vector_categories_chromosomes->{$chr_name}->{$c} if $c eq "intergenic";
			$rocks3->put_batch($c,$final_vector_categories_chromosomes->{$chr->name}->{$c});
		}
	}
	$rocks3->write_batch();
	$rocks3->close();

exit(0);



sub construct_vector {
	my ($project,$r) =@_;
#	$project->buffer->dbh_reconnect();
	#$project->buffer->disconnect();
	#$project->buffer->{debug} =1;
	#my $variations =   $project->myflushobjects($list,"variants");
	my $vector_categories_chromosomes ={};
	for (my $id=$r->[0];$id<=$r->[1];$id++){
	#foreach my $id (@$list){
			my $v = $project->returnVariants($chr->name."!".$id);
		
			my $nb_ho = $v->getGnomadHO();
	 		my $nb_ac = $v->getGnomadAC();
	
	 	
	 		 delete $v->{project} ;
			delete $v->{buffer};
			$v = undef;
	}
	return {vector_categories_chromosomes=>$vector_categories_chromosomes};
} 

sub addBit {
	my ($hash,$chr,$cat,$index) = @_;
	my $chr_name = $chr->name;
	unless (exists $hash->{$chr_name}->{$cat}){
		$hash->{$chr_name}->{$cat} = $chr->getVariantsVector->Shadow();
	}
	$hash->{$chr_name}->{$cat}->Bit_On($index);
}

sub listVariants {
		my($vectors) =@_;
		my @list_variants;
		foreach my $k (keys %$vectors){
			push(@list_variants,@{to_array($vectors->{$k},$k)});
		}
		return \@list_variants;
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