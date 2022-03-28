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
		push(@all_categories,$p->name."_".$k);
	}
}

push(@all_categories,keys %$functional_categorie);


push(@all_categories,@$pathogenic_categories);
push(@all_categories,"genes");
push(@all_categories,"intergenic");
push(@all_categories,"intronic");
push(@all_categories,"acmg");

my $final_vector_categories_chromosomes ={};
#foreach my $c (@all_categories){
#	$g_vector_categories_chromosomes->{$c} = $chr->getVariantsVector->Shadow();
#}


my $list;
my $hv;
#foreach my $chr (@{$project->getChromosomes}){
	my $chr = $project->getChromosome($chr_name);
	$hv->{$chr->name} =  $chr->getVectorVariations();
	

#}
 $list = listVariants($hv);
if  (scalar(@$list) == 0 ){
	#die("zero");
}


my $nb = int(scalar(@$list)/$fork+1);
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
					#warn $h->{vector_categories_chromosomes}->{$c}->to_Enum;
					$final_vector_categories_chromosomes->{$chr}->{$c} = $h->{vector_categories_chromosomes}->{$chr}->{$c}->Shadow unless exists $final_vector_categories_chromosomes->{$chr}->{$c};
					$final_vector_categories_chromosomes->{$chr}->{$c} |= $h->{vector_categories_chromosomes}->{$chr}->{$c};
			}
			#$final_vector_categories_chromosomes->{$chr}->{all} |= $h->{vector_categories_chromosomes}->{$chr}->{all};
		}
		
    }
    
    );





my $iter = natatime($nb, @$list);
	my $id = time;
$project->preload_patients();
	$project->buffer->disconnect();
	
		#$project->buffer->{dbh} ="-";
		
	  while( my @tmp = $iter->() ){
		$id ++;
		$hrun->{$id} ++;
		my $pid = $pm->start and next;
		my $res = construct_vector($project,\@tmp);
		$res->{run_id} = $id;
		$pm->finish(0,$res);
	}
$pm->wait_all_children();



 $chr = $project->getChromosome($chr_name);
#foreach my $chr (@{$project->getChromosome}){ 
	my $no3 = $chr->lmdb_score_impact("c");
	my $chr_name = $chr->name();
	
	foreach my $c (@all_categories){
		unless (exists $final_vector_categories_chromosomes->{$chr_name}->{$c}){
			#warn "vide $c " if $c eq "intergenic";
			$no3->put($c,$chr->getVariantsVector->Shadow());
		}
		else{
			#warn $final_vector_categories_chromosomes->{$chr_name}->{$c} if $c eq "intergenic";
			$no3->put($c,$final_vector_categories_chromosomes->{$chr_name}->{$c});
		}
	}
	$no3->close();
#}
die() if keys %$hrun;
 warn $project->dir_lmdb_score_impact();
#exit(0);

warn "OK";
#die();

sub construct_vector {
	my ($project,$list) =@_;
#	$project->buffer->dbh_reconnect();
	#$project->buffer->disconnect();
	#$project->buffer->{debug} =1;
	#my $variations =   $project->myflushobjects($list,"variants");
	my $vector_categories_chromosomes ={};
	foreach my $id (@$list){
			my $v = $project->returnVariants($id);
			my $nb_ho = $v->getGnomadHO();
	 		my $nb_ac = $v->getGnomadAC();
	 		my $ncboost = $v->ncboost_score();
	 		my $pdv = $v->other_projects;
	 		my $sdv = $v->other_patients ; 
	 		my $sdv_ho = $v->other_patients_ho; 
	 		
	 		$nb_ac = 0 unless $nb_ac;
	 		$nb_ho = 0 unless defined $nb_ho;
	 		my $chr = $v->getChromosome;
	 		my $chr_name = $chr->name;
	 		#ratio allelique
	 		foreach my $p (@{$v->getPatients}){
	 			my $r = $v->getRatio($p);
	 				foreach my $c (keys %$scaled_score_ratio){
	 					my $c1 = $p->name."_".$c;
	 					if ($r >= $scaled_score_ratio->{$c} ){
	 						addBit($vector_categories_chromosomes,$chr,$c1,$v->{index_lmdb});
	 					}
	 				}
	 			
	 		}
	 		
	 		
	 		
			foreach my $c (keys %$scaled_score_gnomad_ho){
	 			
	 				if ($nb_ho <= $scaled_score_gnomad_ho->{$c} ){
	 					addBit($vector_categories_chromosomes,$chr,$c,$v->{index_lmdb});
	 					#$vector_categories_chromosomes->{$c}->Bit_On($v->{index_lmdb});
	 				}
	 				
	 			}
	 			
	 		foreach my $c (sort {$a <=>$b} keys %$scaled_score_project_dejavu){
	 			
	 				if ($pdv <= $scaled_score_project_dejavu->{$c} ){
	 					addBit($vector_categories_chromosomes,$chr,$c,$v->{index_lmdb});
	 					#$vector_categories_chromosomes->{$c}->Bit_On($v->{index_lmdb});
	 				}

	 				
	 			}
	 			
	 			
	 			foreach my $c (sort {$a <=>$b}keys %$scaled_score_sample_dejavu){
	 			
	 				if ($sdv <= $scaled_score_sample_dejavu->{$c} ){
	 					addBit($vector_categories_chromosomes,$chr,$c,$v->{index_lmdb});
	 					#$vector_categories_chromosomes->{$c}->Bit_On($v->{index_lmdb});
	 				}

	 				
	 			}
	 			
	 			foreach my $c (sort {$a <=>$b} keys %$scaled_score_sample_ho_dejavu){
	 				if ($sdv_ho <= $scaled_score_sample_ho_dejavu->{$c} ){
	 					addBit($vector_categories_chromosomes,$chr,$c,$v->{index_lmdb});
	 					#$vector_categories_chromosomes->{$c}->Bit_On($v->{index_lmdb});
	 				}

	 			}	
	 			
	 			
	 			
	 	if ($ncboost and $ncboost ne '-') {
	 		$ncboost = $ncboost * 100;
		 	foreach my $c (keys %$scaled_score_ncboost){
		 		my $value_cat = $scaled_score_ncboost->{$c} * 100;
		 		if ($ncboost >= $value_cat){
		 			addBit($vector_categories_chromosomes,$chr,'ncboost_'.$c,$v->{index_lmdb});
		 		}
		 	}
	 	}	
	 	foreach my $c (keys %$scaled_score_gnomad_ac){
	 		if ($nb_ac <= $scaled_score_gnomad_ac->{$c} ){
	 			addBit($vector_categories_chromosomes,$chr,$c,$v->{index_lmdb});
	 		}
	 	}
	 		 if ($v->isDM){
	 		 	addBit($vector_categories_chromosomes,$chr,"dm",$v->{index_lmdb});
	 		 }
	 		 
	 		  if ($v->hgmd_id){
	 		  	addBit($vector_categories_chromosomes,$chr,"hgmd",$v->{index_lmdb});
	 		 }
	 		 if ($v->score_clinvar == 5) {
	 		 	 addBit($vector_categories_chromosomes,$chr,"clinvar_pathogenic",$v->{index_lmdb});
	 		 }
	 		 my $mask = $v->{annotation}->{all}->{mask} ;
	 		 foreach my $f ( keys %{$project->buffer->config->{functional_annotations}}){
	 		 	 if ($mask & $project->getMaskCoding($f)){
	 		 		addBit($vector_categories_chromosomes,$chr,$f,$v->{index_lmdb});
	 		 	}
	 		 }
	 		 
	 		
	 		 if (scalar(@{$v->getGenes}) >0 ) {
	 		 	 addBit($vector_categories_chromosomes,$chr,"genes",$v->{index_lmdb});
	 		 	  foreach my $g (@{$v->getGenes}) {
	 		 		next unless exists $g->panels_name->{'ACMG-Actionable'};
	 		 		addBit($vector_categories_chromosomes,$chr,"acmg",$v->{index_lmdb});
	 		 	}
	 		 	$project->buffer->disconnect();
	 		 }
	 		 else {
	 		 	
	 		 	addBit($vector_categories_chromosomes,$chr,"intergenic",$v->{index_lmdb});
	 		 }
	 		 
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