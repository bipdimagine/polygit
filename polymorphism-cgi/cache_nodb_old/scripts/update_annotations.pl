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

my $chr = $project->getChromosome($chr_name);
#my $no = $chr->lmdb_variations("r");
my $no = $chr->lmdb_variations("r");
my $no3 = $chr->lmdb_score_impact("c");
$no3->put("date",time);
exit(0) unless $no;
my $ranges = $no->ranges($fork*2);
$no->close();
my $pm = new Parallel::ForkManager($fork);
my $dir_out = $project->lmdb_cache_variations_dir()."/test/";
$dir_out = "/tmp/";
my $f1 = $dir_out."/".$chr_name;
unlink $f1."_index" if -e $f1."_index";
unlink $f1 if -e $f1;
my $intspan_categories_transcripts;
my $intspan_patients;

my $impact_factors = $project->buffer->config->{impact_factors};
my $scaled_score_impact = $project->buffer->config->{scaled_score_impact};
my $scaled_score_frequence_public =  $project->buffer->config->{scaled_score_frequence_public};
my $hgmd_public_dm =  $project->buffer->config->{scaled_score_frequence_public};
my $scaled_score_ratio =   $project->buffer->config->{scaled_score_ratio};
my $scaled_score_in_this_run =  $project->buffer->config->{scaled_score_in_this_run};
my $vector_categories_chromosomes;
my $vector_categories_chromosomes;
 
 my @pathogenic_categories =("hgmd","clinvar_pathogenic","dm");

foreach my $c (keys %$scaled_score_frequence_public){
	$vector_categories_chromosomes->{$c} = $chr->getVariantsVector->Shadow()
}

foreach my $c (keys %$scaled_score_in_this_run){
	$vector_categories_chromosomes->{$c} = $chr->getVariantsVector->Shadow()
}

foreach my $c (@pathogenic_categories){
	$vector_categories_chromosomes->{$c} = $chr->getVariantsVector->Shadow()
}



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
	 		my $f = $o->frequency;
	 		my $r = $o->in_this_run_ratio;
	 		$f =0 unless $f;
	 		my $fp = $o->exome_patients()/$o->total_exome_patients();
	 		foreach my $c (keys %$scaled_score_frequence_public){
	 			
	 			if ($f <= $scaled_score_frequence_public->{$c} && $fp <= $scaled_score_frequence_public->{$c}  ){
	 				$vector_categories_chromosomes->{$c}->Bit_On($o->{index_lmdb});
	 			}
	 		}
	 		 if ($o->isDM){
	 		 	$vector_categories_chromosomes->{dm}->Bit_On($o->{index_lmdb});
	 		 }
	 		 
	 		  if ($o->hgmd_id){
	 		 	$vector_categories_chromosomes->{hgmd}->Bit_On($o->{index_lmdb});
	 		 }
	 		 if ($o->score_clinvar == 5) {
	 		 	$vector_categories_chromosomes->{clinvar_pathogenic}->Bit_On($o->{index_lmdb});
	 		 }
	 	
	 	
	 	foreach my $c (keys %$scaled_score_in_this_run){
	 			if ($r <= $scaled_score_in_this_run->{$c} ) {
	 					
	 			}
	 	}
	 	

	 	}
	 	$no->close();
		$no = undef;

	 	foreach my $tr (keys %{$h->{intspan_transcripts}}) {
	 		foreach my $c (keys %{$h->{intspan_transcripts}->{$tr} } ) {
	 			$intspan_categories_transcripts->{$c}->{$tr} =Set::IntSpan::Fast::XS->new()  unless exists $intspan_categories_transcripts->{$c}->{$tr};
	 			$intspan_categories_transcripts->{$c}->{$tr} = $intspan_categories_transcripts->{$c} ->{$tr}->union($h->{intspan_transcripts}->{$tr}->{$c});
	 		}
	 	}
	 	
	 	foreach my $p (keys %{$h->{intspan_patients}}) {
	 		foreach my $c (keys %{$h->{intspan_patients}->{$p} } ) {
	 			$intspan_patients->{$p}->{$c} =Set::IntSpan::Fast::XS->new()  unless exists $intspan_patients->{$p}->{$c};
	 			$intspan_patients->{$p}->{$c} = $intspan_patients->{$p} ->{$c}->union($h->{intspan_patients}->{$p}->{$c});
	 		}
	 	}
		
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
	my $intspan_transcripts;
	my $intspan_chromosome;
	my $intspan_patients;
	for (my $i=$r->[0];$i<$r->[1]+1;$i++){
		$nb ++;
		warn $nb if $nb%1000 == 0;
		my $v = $no->get_index($i);
		$v->{project} =  $project;
		$v->{buffer} = $buffer;
		#purge data
		
		$v->purge_deja_vu();
		$v->purge_public_data();
		# delete $v->{scaled_score_frequence_public};
		foreach my $patient (@{$v->getPatients}){
				  		unless (exists $intspan_transcripts->{$patient->id}){
		  		
		  					foreach my $c (keys %$scaled_score_ratio){
		  						$intspan_patients->{$patient->id}->{"ratio_".$c} = Set::IntSpan::Fast::XS->new();
		  					}
		  					$intspan_patients->{$patient->id}->{"ratio_all"} = Set::IntSpan::Fast::XS->new();
				  		}
				  		
				  		foreach my $c (keys %$scaled_score_ratio){
				  			if ($v->getRatio($patient) >= $scaled_score_ratio->{$c}) {
				  				$intspan_patients->{$patient->id}->{"ratio_".$c}->add($i);
				  			}
				  		}
		}
	
	
	
		my $chr_name = $v->getChromosome()->name;
		
 	
		  delete $v->{isPublic};
	#	  delete $v->{scaled_score_frequence_public};
		  
		$v->gnomad;
		$v->dejaVuInfosForDiag();
		 $v->name();
		  $v->cosmic();
		  $v->getGenes();
		#  delete $v->{isPublic};
		 # $v->frequency();
		  $v->score_clinical_local();
		  $v->score_clinvar();
		  $v->text_clinvar();
		  $v->comment_clinical_local();
		  $v->cadd_score();
		 $v->hgmd_id();
		 $v->get_codon_text(1);
		 $v->get_codon_text(-1);
	
		 #warn $i."/".$r->[1];
		  $v->annotation();
		 delete $v-> {scaled_score_frequence_public};
		  
		  	foreach my $tr ( @{$v->getTranscripts}){
		  		unless (exists $intspan_transcripts->{$tr->id}){
		  		
		  			foreach my $c (keys %$scaled_score_impact){
		  				$intspan_transcripts->{$tr->id}->{$c} = Set::IntSpan::Fast::XS->new();
		  			}
		  			$intspan_transcripts->{$tr->id}->{all} = Set::IntSpan::Fast::XS->new();
		  		}
		  		
		  		 foreach my $patient (@{$v->getPatients}){
		  	 			$v->scaledScoreVariant($tr,$patient,$vquery);
		  		 }
		  		 
		  		 	foreach my $c (keys %$scaled_score_impact){
		  		 		$intspan_transcripts->{$tr->id}->{$c}->add($i) if $v->isEffectImpact($tr,$c);
		  		 	}
		  		 $intspan_transcripts->{$tr->id}->{low}->add($i);
		  		  $intspan_transcripts->{$tr->id}->{all}->add($i);
		  		 $v->getNomenclature($tr);
		  	 	foreach my $p (@{$tr->getProteins}){
		  	 		$v->polyphenScore($p);
		  	 		 $v->siftScore($p);
		  	 	}
		  	 	
		  	}
		  
		  
		
		
		$v->getSequence();
		$v->ref_allele();
		$v->in_this_run_ratio();
		$v->dejaVuInfosForDiag();
		
		### 
		# deja vu 
		####
		
		
		
	
		
		my $codon;
		delete $v->{project} ;
		delete $v->{buffer};
		push(@{$h->{obj}},$v);
		$h->{intspan_transcripts} = $intspan_transcripts;
		$h->{intspan_patients} = $intspan_patients;
	}
	$no->close();
	$no = undef;
	$pm->finish(0,$h);
	die();
	#die();
}
$pm->wait_all_children();
$project->buffer->dbh_reconnect();

warn "end variants";

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

warn "save vector";
my $chr = $project->getChromosome($chr_name);

my $no3 = $chr->lmdb_score_impact("c");

my $h1;
my $h2;
foreach my $c (keys %{$intspan_categories_transcripts}){
	$h2->{$c} ++;
	foreach my $tr (keys %{$intspan_categories_transcripts->{$c}}){
			my $debug ;
			$debug =1 if $tr eq "ENST00000369037_6";
			warn "coucou " if $debug;
		my ($id,$chr_name) = split("_",$tr);
		$h2->{$c} ++;
		#next if $intspan_categories_transcripts->{$c}->{$tr}->is_empty();
		my $v;
		unless ($intspan_categories_transcripts->{$c}->{$tr}->is_empty){
			 $v =  Bit::Vector->new_Enum($chr->size_vector(), join(',', $intspan_categories_transcripts->{$c}->{$tr}->as_string));
			# $no3->put($tr."_".$c,$v);
			
			 $h1->{$tr}->{$c} = $v;
		}
		else {
			$v = $chr->getNewVector();
		}
		  $no3->put($tr."_".$c,$v);
		  if ($c eq "all"){
		  	 $no3->put($tr,$v);
		  }
		  
		#$no3->put($chr_name,$tr."_".$c,$v);
		#$no3->put($tr."_".$c,$intspan_categories_transcripts->{$c}->{$tr});
		
	}
}
warn Dumper $no3->get("ENST00000369037_6_moderate");

foreach my $p (keys %{$intspan_patients}){
	
	foreach my $c (keys %{$intspan_patients->{$p}}){
		#next if $intspan_categories_transcripts->{$p}->{$c}->is_empty();
		my $v;
		unless ($intspan_patients->{$p}->{$c}->is_empty){
			 $v =  Bit::Vector->new_Enum($chr->size_vector(), join(',', $intspan_patients->{$p}->{$c}->as_array));
		}
		else {
			$v = $chr->getNewVector();
			
			
		}
	
		 $no3->put($p."_".$c,$v);
		  if ($c eq "all"){
		  	 $no3->put($p,$v);
		  }
		#$no3->put($chr_name,$tr."_".$c,$v);
		#$no3->put($tr."_".$c,$intspan_categories_transcripts->{$c}->{$tr});
		
	}

}


foreach my $c (keys %$scaled_score_frequence_public){

	$no3->put($c,$vector_categories_chromosomes->{$c});
}

foreach my $c (@pathogenic_categories){
	#warn $c." ".$vector_categories_chromosomes->{$c};
		warn $c;
		$no3->put($c,$vector_categories_chromosomes->{$c});
	#$vector_categories_chromosomes->{$c} = $chr->getVariantsVector->Shadow()
}
$no3->close();
my $f1 = $dir_out."/".$chr_name;
my $f2 = $project->lmdb_cache_variations_dir()."/".$chr_name;
system("mv $f1 $f2 ");
die() unless -e $f2;
$f1 = $dir_out."/".$chr_name."_index";
$f2 = $project->lmdb_cache_variations_dir()."/".$chr_name."_index";
system("mv $f1 $f2 ");







exit(0);


#verif 


my $v = $no->get_index(1);

$v->{project} =  $project;
$v->{buffer} = $buffer;

warn $v->name();
warn $v->getGenes->[0]->name;


my $cursor = $no->cursor( 0, 10);

while (my $v = $cursor->next_key()){
		
	
}


warn $no->get_key_index(1);
die();

foreach my $var (@{$chr->getListVarObjects($chr->getVariantsVector)}) {
	warn $var;
}
foreach my $v (@{$chr->getStructuralVariations}){
	
	warn $v->name;
	warn $v->{index_lmdb};
	
}



