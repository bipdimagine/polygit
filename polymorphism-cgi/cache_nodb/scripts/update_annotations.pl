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
 my $host = hostname();


warn "*_*_*_*_*_".$host."*_*_*_*_*_";
#use Cache_Commons;

my $fork = 1;
my $force;
my ($project_name, $chr_name);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'chr=s'        => \$chr_name,
	'force=s'  => \$force,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
#unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }
my $t = `hostname`;

my $nbErrors = 0;
my $buffer = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProjectCache( -name => $project_name, -cache => '1', -typeFilters => 'familial' );



#
my $chromosomes = $project->getChromosomes;

foreach my $chr (@{$chromosomes}){
	
	if ($chr_name){
		next if $chr->name ne "$chr_name";
	}
	run_update($chr->name);
}
exit(0) if ($chr_name);
#system("$RealBin/tree_cache.pl -project=$project_name -fork=$fork");


sub run_update {
	my ($chr_name) = @_;
	
my $project = $buffer->newProjectCache( -name => $project_name, -cache => '1', -typeFilters => 'familial' );
	my $chr = $project->getChromosome($chr_name);
#my $no = $chr->lmdb_variations("r");
my $no = $chr->lmdb_variations("r");

my $no3 = $chr->lmdb_score_impact("c");
$no3->put("date",time);
$no3->close;
exit(0) unless $no->exists_db;

my $ranges = $no->ranges($fork);
#warn Dumper $ranges;
$no->close();
my $pm = new Parallel::ForkManager($fork);
my $dir_out = $project->lmdb_cache_variations_dir()."/test/";
$dir_out = "/tmp/";
my $f1 = $dir_out."/".$chr_name;
unlink $f1."_index" if -e $f1."_index";
unlink $f1 if -e $f1;
my $intspan_categories_transcripts;
my $intspan_patients;


my $scaled_score_impact = $project->buffer->config->{scaled_score_impact};
#warn Dumper $scaled_score_impact;
#die();
my $scaled_score_frequence_public =  $project->buffer->config->{scaled_score_frequence_public};
my $scaled_score_gnomad_ho =  $project->buffer->config->{scaled_score_gnomad_ho_ac};
my $scaled_score_gnomad_ac =  $project->buffer->config->{scaled_score_gnomad_ac};
my $scaled_score_ratio =   $project->buffer->config->{scaled_score_ratio};
my $scaled_score_in_this_run =  $project->buffer->config->{scaled_score_in_this_run};
my $functional_categorie = $project->buffer->config->{functional_annotations};

my $g_vector_categories_chromosomes;
 
 my @pathogenic_categories =("hgmd","clinvar_pathogenic","dm");
my @all_categories;
push(@all_categories,keys %$scaled_score_frequence_public);
push(@all_categories,keys %$scaled_score_impact);
push(@all_categories,keys %$scaled_score_gnomad_ac);
push(@all_categories,keys %$scaled_score_gnomad_ho);
push(@all_categories,keys %$scaled_score_in_this_run);
push(@all_categories,keys %$functional_categorie);
#warn $functional_categorie;
push(@all_categories,@pathogenic_categories);

foreach my $c (@all_categories){
	$g_vector_categories_chromosomes->{$c} = $chr->getVariantsVector->Shadow();
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
	 	


	 	}
	 	$no->close();
		$no = undef;
		my $hintspan = $h->{vector_categories_chromosomes};
		foreach my $c (@all_categories){
				#warn $h->{vector_categories_chromosomes}->{$c}->to_Enum;
				$g_vector_categories_chromosomes->{$c} |= $h->{vector_categories_chromosomes}->{$c};
		}
	 	foreach my $tr (keys %{$h->{intspan_transcripts}}) {
	 		my ($id,$chr) = split("_",$tr);
	 		 
	 		foreach my $c (keys %{$h->{intspan_transcripts}->{$tr} } ) {
	 			$intspan_categories_transcripts->{$c}->{$tr} =Set::IntSpan::Fast::XS->new()  unless exists $intspan_categories_transcripts->{$c}->{$tr};
	 			$intspan_categories_transcripts->{$c}->{$tr} = $intspan_categories_transcripts->{$c} ->{$tr}->union($h->{intspan_transcripts}->{$tr}->{$c});
	 			my @z = $intspan_categories_transcripts->{$c}->{$tr}->as_array;
	 			foreach my $x (@z){
	 				$g_vector_categories_chromosomes->{$c}->Bit_On($x);
	 			}
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
my $vector_categories_chromosomes;
	 		foreach my $c (@all_categories){
	 			warn $c;
				$vector_categories_chromosomes->{$c} = $chr->getVariantsVector->Shadow();
			}
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
		warn $nb if $nb%10000 == 0;
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
 	
		
		 $v->gnomad;
		 #warn Dumper $v->gnomad	if $v->id eq "7_151433401_G_A";
		$v->dejaVuInfosForDiag();
		 $v->name();
		  $v->cosmic();
		  $v->getGenes();
		  $v->score_clinical_local();
		  $v->score_clinvar();
		  $v->text_clinvar();
		  $v->comment_clinical_local();
		  $v->cadd_score();
		  $v->dbscsnv_ada;
		  $v->dbscsnv_rf;
		   $v->revel_score;
		 $v->hgmd_id();
		 $v->get_codon_text(1);
		 $v->get_codon_text(-1);
	
		 #warn $i."/".$r->[1];
		  $v->annotation();
	#	 delete $v-> {scaled_score_frequence_public};


		  	foreach my $tr ( @{$v->getTranscripts}){
		  		warn $tr->id." ".$tr->getGene->external_name if $debug;
		  		unless (exists $intspan_transcripts->{$tr->id}){
		  		
		  			foreach my $c (keys %$scaled_score_impact){
		  				$intspan_transcripts->{$tr->id}->{$c} = Set::IntSpan::Fast::XS->new();
		  			}
		  			 foreach my $c (keys %$functional_categorie){
		  			 	$intspan_transcripts->{$tr->id}->{$c} = Set::IntSpan::Fast::XS->new();
		  			 }
		  			$intspan_transcripts->{$tr->id}->{all} = Set::IntSpan::Fast::XS->new();
		  		}
		  		
		  		 #foreach my $patient (@{$v->getPatients}){
		  	 		#	$v->scaledScoreVariant($tr,$patient,$vquery);
		  		 #}
		  		 
		  		 	foreach my $c (keys %$scaled_score_impact){
		  		 		$intspan_transcripts->{$tr->id}->{$c}->add($i) if $v->isEffectImpact($tr,$c);
		  		 	}
		  		 foreach my $c (keys %$functional_categorie){
		  		 		$intspan_transcripts->{$tr->id}->{$c}->add($i) if $v->isConsequence($tr,$c);
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
		# chromosome vecteurs 
		####
		 my $o = $v;
			my $f = $o->frequency;
	 		$f = 0 unless $f;
	 		my $r = $o->in_this_run_ratio;
	 		$o->gnomad;
	 		
	 		
	 		
	 		$f =0 unless $f;
	 		my $fp = $o->exome_patients()/$o->total_exome_patients();
	 		my $nb_ho = $o->getGnomadHO();
	 		my $nb_ac = $o->getGnomadAC();
	 		$nb_ac = 0 unless $nb_ac;
	 		$nb_ho = 0 unless defined $nb_ho;
	 		#warn $nb_ho." ".$nb_ac;
	 		
	 		#chromosome categorie
	 		
	 		
	 		foreach my $c (keys %$scaled_score_frequence_public){
	 			
	 			if ($f <= $scaled_score_frequence_public->{$c} && $fp <= $scaled_score_frequence_public->{$c}  ){
	 				$vector_categories_chromosomes->{$c}->Bit_On($o->{index_lmdb});
	 			}
	 		}
	 		
	 		 foreach my $c (keys %$functional_categorie){
		  		 	$vector_categories_chromosomes->{$c}->Bit_On($o->{index_lmdb})  if $v->isConsequence("all",$c);
		  	 }
		  			
	 			foreach my $c (keys %$scaled_score_gnomad_ho){
	 			
	 				if ($nb_ho <= $scaled_score_gnomad_ho->{$c} ){
	 			
	 					$vector_categories_chromosomes->{$c}->Bit_On($o->{index_lmdb});
	 				}
	 				
	 			}
	 			
	 	foreach my $c (keys %$scaled_score_gnomad_ac){
	 		if ($nb_ac <= $scaled_score_gnomad_ac->{$c} ){
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
		
	
		######
		# END
		#######
		my $codon;
		delete $v->{project} ;
		delete $v->{buffer};
		push(@{$h->{obj}},$v);
		$h->{intspan_transcripts} = $intspan_transcripts;
		$h->{intspan_patients} = $intspan_patients;
		
	}
	$no->close();
	$no = undef;
	$h->{vector_categories_chromosomes} = $vector_categories_chromosomes;
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


my $chr = $project->getChromosome($chr_name);

my $no3 = $chr->lmdb_score_impact("w");

my $tree_transcripts;

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
		  	 my @z = $intspan_categories_transcripts->{$c}->{$tr}->as_array;
		  	 push(@{$tree_transcripts->{$c}},[$tr,$z[0],$z[-1]]);
		  }
		  
		#$no3->put($chr_name,$tr."_".$c,$v);
		#$no3->put($tr."_".$c,$intspan_categories_transcripts->{$c}->{$tr});
		
	}
}

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
foreach my $c (@all_categories){
	$no3->put($c,$g_vector_categories_chromosomes->{$c});
}


my $f1 = $dir_out."/".$chr_name;
my $f2 = $project->lmdb_cache_variations_dir()."/".$chr_name;
system("mv $f1 $f2 ");
die() unless -e $f2;
$f1 = $dir_out."/".$chr_name."_index";
$f2 = $project->lmdb_cache_variations_dir()."/".$chr_name."_index";
system("mv $f1 $f2 ");



foreach my $chr (@{$project->getChromosomes}) {
	my $tree;

	my @keys = grep {$_ =~/_all/} @{$no3->get_keys()};
	foreach my $k (@keys){
		next unless $k =~/ENST/;
		my $v = $no3->get($k);
		my @pos =  $v->Index_List_Read();
		$k =~s/_all//;
		push(@$tree,[$k,$pos[0],$pos[-1]+1]);
	}
		
		$no3->put("intervaltree_transcripts".$chr->name,$tree);
}
$no3->close();

return 1;
}