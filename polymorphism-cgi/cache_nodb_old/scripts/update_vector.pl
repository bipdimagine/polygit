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

#use Cache_Commons;

my $fork = 1;
my ($project_name, $chr_name);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'chr=s'        => \$chr_name,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }
my $t = `hostname`;

my $nbErrors = 0;
my $buffer = new GBuffer;

my $project = $buffer->newProjectCache( -name => $project_name, -cache => '1', -typeFilters => 'familial' );


my $scaled_score_ratio =   $buffer->config->{scaled_score_ratio};
my $scaled_score_impact = $project->buffer->config->{scaled_score_impact};
my $scaled_score_frequence_public =  $project->buffer->config->{scaled_score_frequence_public};
my $scaled_score_in_this_run =  $project->buffer->config->{scaled_score_in_this_run};
my $chr = $project->getChromosome($chr_name);



my $no3 = $chr->lmdb_score_impact("c");
 $no3->put("date",time);
 #my $z = $no3->get("dm");
 #warn $z;
 #die();
 #$no3->close();
my $no = $chr->lmdb_variations("r");
my $vector_ratio;
foreach my $p (@{$project->getPatients}){
	foreach my $r (keys %$scaled_score_ratio){
		$vector_ratio->{$p->id."_ratio_".$scaled_score_ratio->{$r}} = $chr->getNewVector();
	}
}

my $vector_impact; 
my $vector_frequences;
my $vector_impact;
my $vector_clinical;
my $vector_in_this_run;
my $vector_genes;
my $vector_transcripts;
my $nb =0;
while (my $h = $no->get_next()){
	my ($v) = values %$h;


	
#	$project->buffer->dbh_reconnect();
	

		$nb ++;
		warn $nb if $nb%1000 == 0;
		$v->{project} =  $project;
		$v->{buffer} = $buffer;
		my $i = $v->{index_lmdb};
		#vector ratio coverage 
		foreach my $p (@{$v->getPatients}){
			foreach my $c (keys %$scaled_score_ratio){
			if ($v->getRatio($p) >= $scaled_score_ratio->{$c}) {
						$vector_ratio->{$p->id."_ratio_".$scaled_score_ratio->{$c}}->Bit_On($i);
				}
			}
		}
		
		#vector 
		foreach my $tr ( @{$v->getTranscripts}){
			 $vector_transcripts->{$tr->id} = $chr->getNewVector() unless exists $vector_transcripts->{$tr->id};
			 $vector_transcripts->{$tr->id}->Bit_On($i);
			 
			foreach my $c (keys %$scaled_score_impact){
						$vector_impact->{$tr->id."_".$c} =  $chr->getNewVector() unless exists $vector_impact->{$tr->id."_".$c};
						$vector_impact->{$tr->id."_".$c}->Bit_On($i)  if $v->isEffectImpact($tr,$c);
		  	}
		}
		
		foreach my $gene ( @{$v->getGenes}){
			 $vector_genes->{$gene->id} = $chr->getNewVector() unless exists $vector_genes->{$gene->id};
			 $vector_genes->{$gene->id}->Bit_On($i);
			$vector_impact->{$gene->id."_all"} =  $chr->getNewVector() unless exists $vector_impact->{$gene->id."_all"};
			$vector_impact->{$gene->id."_all"}->Bit_On($i);
			foreach my $c (keys %$scaled_score_impact){
						$vector_impact->{$gene->id."_".$c} =  $chr->getNewVector() unless exists $vector_impact->{$gene->id."_".$c};
						$vector_impact->{$gene->id."_".$c}->Bit_On($i)  if $v->isEffectImpact($gene,$c);
		  	}
		}
		
			#vector frequency
			
			my $f = $v->frequency;
	 		$f =0 unless $f;
	 		foreach my $c (keys %$scaled_score_frequence_public){
	 			
	 			if ($f <= $scaled_score_frequence_public->{$c} ){
	 				$vector_frequences->{"freq_".$c}=  $chr->getNewVector() unless exists $vector_frequences->{"freq_".$c};
	 				$vector_frequences->{"freq_".$c}->Bit_On($v->{index_lmdb});
	 			}
	 		}
	 		
	 		#vector clinical
	 		
	 		 if ($v->isDM){
	 		 	$vector_clinical->{dm}=  $chr->getNewVector() unless exists $vector_clinical->{dm};
	 		 	$vector_clinical->{dm}->Bit_On($v->{index_lmdb});
	 		 }
	 		 
	 		  if ($v->hgmd_id){
	 		  	 	$vector_clinical->{hgmd}=  $chr->getNewVector() unless exists $vector_clinical->{hgmd};
	 		 		$vector_clinical->{hgmd}->Bit_On($v->{index_lmdb});
	 		 }
	 		 
	 		 if ($v->score_clinvar == 5) {
	 		 	$vector_clinical->{clinvar_pathogenic} =  $chr->getNewVector() unless exists $vector_clinical->{clinvar_pathogenic};
	 		 	warn $vector_clinical->{clinvar_pathogenic}->Size();
	 		 	$vector_clinical->{clinvar_pathogenic}->Bit_On($v->{index_lmdb})
	 		 }
	 		 
	 	#vector run_ratio ?	 
	 	
	 	my $r = $v->in_this_run_ratio;
	 	foreach my $c (keys %$scaled_score_in_this_run){
	 			if ($r <= $scaled_score_in_this_run->{$c} ) {
	 					$vector_in_this_run->{"this_".$c}=  $chr->getNewVector() unless exists $vector_in_this_run->{"this_".$c};
	 					$vector_in_this_run->{"this_".$c}->Bit_On($v->{index_lmdb});
	 			}
	 	}
		
		#warn Dumper $v->{getRatio}; 
		delete $v->{project} ;
		delete $v->{buffer};
		$v = undef;
}


#my $no3 = $chr->lmdb_score_impact("w");
	foreach my $r (keys %$vector_ratio){
      $no3->put($r,$vector_ratio->{$r});	
	}

foreach my $r (keys %$vector_impact){
      $no3->put($r,$vector_impact->{$r});	
}


foreach my $r (keys %$vector_frequences){
      $no3->put($r,$vector_frequences->{$r});	
}

foreach my $r (keys %$vector_clinical){
      $no3->put($r,$vector_clinical->{$r});	
}

foreach my $r (keys %$vector_in_this_run){
      $no3->put($r,$vector_in_this_run->{$r});	
}

#my $no4 = GenBoNoSqlIntervalTree->new(dir=>$project->dir_lmdb_score_impact(),mode=>"c");
my $array;
foreach my $tr (keys %$vector_transcripts){
	#warn $vector_transcripts->{$tr}->to_Enum();
	$no3->put($tr,$vector_transcripts->{$tr});
	my ($start,$end) = split("-", $vector_transcripts->{$tr}->to_Enum());
	push(@$array,[$tr,$start,$end+1]);
} 
$no3->put("transcritps_vector_annotations",$array);
#$no4->put("vector_annotations","transcripts_".$chr->name,$array);
#my $no3 = GenBoNoSqlIntervalTree->new(dir=>$sqliteDir,mode=>"c");
my $array_genes;
foreach my $gene (keys %$vector_genes){
	#warn $vector_transcripts->{$tr}->to_Enum();
	$no3->put($gene,$vector_genes->{$gene});	
	my ($start,$end) = split("-", $vector_genes->{$gene}->to_Enum());
	push(@$array_genes,[$gene,$start,$end+1]);
}
$no3->put("genes_vector_annotations",$array_genes);
#$no3->close();
#$no4->close();