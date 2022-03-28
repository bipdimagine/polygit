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
use lib "$RealBin/../../../GenBo/lib/obj-nodb/packages/";
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
require "$RealBin/../../../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
my $host = hostname();


#warn "*_*_*_*_*_".$host."*_*_*_*_*_";
#use Cache_Commons;

my $fork = 1;
my $force;
my ($project_name, $chr_name, $annot_version);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'chr_name=s'        => \$chr_name,
	'annot_version=s'    => \$annot_version,
	'force=s'  => \$force,
);



my @header_transcripts = ("consequence","enst","nm","ccds","appris","exon","nomenclature","codons","codons_AA", "polyphen","sift","cadd","revel","dbscsnv");
unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
#unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }
my $t = `hostname`;
warn $t;
my $nbErrors = 0;
my $buffer = GBuffer->new();
my $project = $buffer->newProjectCache( -name => $project_name, -cache => '1', -typeFilters => 'familial' );

if ($annot_version) {
	$project->changeAnnotationVersion($annot_version);
}



my $file;
foreach my $chr (@{$project->getChromosomes} ){
	if ($chr_name ne "all"){
		next if  $chr->name ne $chr_name;
	}
	
	 $file = run_update_chr($chr->name);
	warn "OK $file" if -e $file;
	warn "NOT OK " unless -e $file;
}

warn "END seems prefect exit now ...";
warn "OK $file" if -e $file;
warn "NOT OK " unless -e $file;

exit(0);

sub run_update_chr {
	my ($chr_name) = @_;
my $chr = $project->getChromosome($chr_name);
$project->validations;
my $a = "s:146487:ENSG00000227855_7-9802";

my $pm = new Parallel::ForkManager($fork);

my $dir_tmp = tempdir( CLEANUP => 1,DIR => "/tmp/" );
my $proc;
my $hgenes;
my $zh;
$pm->run_on_finish(
    sub { 
    	my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$h)=@_;
    	unless (defined($h) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
		warn "end " .$h->{index};
		delete $proc->{$h->{index}};
		foreach my $k (keys %{$h->{ids}}){
			push(@{$zh->{$k}},@{$h->{ids}->{$k}});
		}
		my $global_genes = $h->{genes};
		
		foreach my $g (keys %$global_genes){
			foreach my $p (keys %{$global_genes->{$g} } ){
				unless (exists $hgenes->{$g}->{$p}){
					$hgenes->{$g}->{$p} = $global_genes->{$g}->{$p};
					next;
				}
				else {
					if($hgenes->{$g}->{$p}->{father}->{score} < $global_genes->{$g}->{$p}->{father}->{score}){
						$hgenes->{$g}->{$p}->{father} = $global_genes->{$g}->{$p}->{father};
					}
					if($hgenes->{$g}->{$p}->{mother}->{score} < $global_genes->{$g}->{$p}->{mother}->{score}){
						$hgenes->{$g}->{$p}->{mother} = $global_genes->{$g}->{$p}->{mother};
					}
				}
			}	
		}
		
		
    }
    
    );
    
$project->getPatients;

my $list =listVariants($chr);
#$project->newVariant($list->[0]) if @$list;
my $nb = int(scalar(@$list)/($fork));
$nb = 20 if $nb == 0;
my $iter = natatime($nb, @$list);

my $part =0;
my $index = 0;
my @files ;
	my $javascript_id = time + int(rand(10000));
	warn $a;
	$project->preload_patients();
	$project->buffer->dbh_deconnect();
	
	#die();
  while( my @tmp = $iter->() ){
	$part ++;
	warn $part;
	my $f = $chr->name."_".$part;
	push(@files,$f);
	$index ++;
	$proc->{$index} ++;
	my $pid = $pm->start and next;
	warn "start ".$proc->{$index};
	my $nb =0;
	$project->buffer->dbh_reconnect();
	
	#$project->buffer->{dbh} = "-";
	my $no  = GenBoNoSqlLmdb->new(dir=>$dir_tmp,mode=>"c",name=>$f,is_compress=>1);
	my $array;
	$project->setListVariants(\@tmp);
	my $max =scalar(@tmp);
	my $nb =1 ;
	my $global_genes;
	my $zh;
	#sleep(2);
	while(my $v = $project->nextVariant){
		#warn $v->name;
		$nb++;
		if ($nb%1000 ==0){
		#	$no->close();
		#	$no  = GenBoNoSqlLmdb->new(dir=>$dir_tmp,mode=>"w",name=>$f,is_compress=>1);
			print "$nb/$max \n" ;
		}
		
		
	#	warn $v->id;
		my $vmask;
		foreach my $patient (@{$v->getPatients}) {
			my $fam = $patient->getFamily();
		 	my $h =  update_variant_editor::construct_hash_variant ( $project, $v,undef,$patient,1);
		 	$h->{vector_id} = $v->vector_id;
		 	$h->{global_vector_id} = $v->global_vector_id;
		 	update_variant_editor::vvarsome($h,$patient);
		
			my $mother = $patient->getFamily->getMother();
			my $father = $patient->getFamily->getFather();
		 
		 	
		 	my $agenes =[];
		 	my $score_genes;
		 	foreach my $g (@{$v->getGenes}){
		 		my $debug;
		 		$debug = 1 if $g->id eq "ENSG00000167522_16";
		 		#warn $v->id if $debug;
		 		#warn $g->name if $debug;
		 		#warn $patient->name if $debug;
		 		my $hgenes;
				$javascript_id ++;
				my $score_variant = $v->scaledScoreVariant($g,$patient);
				#print "".$v->id."\t".$score_variant."\n" if $debug;
				my $parent_score = 0;
				$parent_score = $score_variant  if $score_variant >0 &&  $v->isConsequence($g,"medium");
				my $nid = "s:".$v->{vector_id}.":".$g->id."@".$patient->name;
				$h->{scaled_score}->{$g->id} = $score_variant;
				
				$no->put($nid,{id=>$v->id,score=>$score_variant,parent_score=>$parent_score,id1=>"s:".$v->{vector_id}.":".$g->id,patient=>$patient->name,type=>"genes",pp=>$patient->name});
			#	warn $g->name." ".$score_variant if $debug;
				#my $score_variant = 0;
				#$res->{variants}->{$v->id}->{score}->{$g->id} = $score_variant;
				
				unless (exists $global_genes->{$patient->name}->{$g->id}) {
				#	$buffer->dbh_reconnect();
					$global_genes->{$patient->name}->{$g->id}->{id} = $g->id;
					$global_genes->{$patient->name}->{$g->id}->{name} = $g->external_name;
					$global_genes->{$patient->name}->{$g->id}->{score} = $g->score;
					$global_genes->{$patient->name}->{$g->id}->{omim_inheritance} = $g->omim_inheritance;
					$global_genes->{$patient->name}->{$g->id}->{omim_inheritance} = $g->omim_inheritance;
					$global_genes->{$patient->name}->{$g->id}->{external_name} = $g->external_name;
					$global_genes->{$patient->name}->{$g->id}->{pLI} = $g->pLI;
					$global_genes->{$patient->name}->{$g->id}->{omim_id} = $g->omim_id;
					$global_genes->{$patient->name}->{$g->id}->{panels} = $buffer->queryPanel()->getPanelsForGeneName($g->external_name);
					$global_genes->{$patient->name}->{$g->id}->{js_id} = $javascript_id."_".$g->id;
					$global_genes->{$patient->name}->{$g->id}->{father}->{score} = -100;
					$global_genes->{$patient->name}->{$g->id}->{mother}->{score} = -100;
					$global_genes->{$patient->name}->{$g->id}->{description} = $g->description;
					$global_genes->{$patient->name}->{$g->id}->{phenotypes} = $g->phenotypes;
					$global_genes->{$patient->name}->{$g->id}->{vector} = $g->getCurrentVector();
					if ($patient->isMale && $chr->name eq "X" && $chr->isAutosomal($g->start,$g->end) ){
							$global_genes->{$patient->name}->{$g->id}->{no_compound} ++;
					}
					unless($mother){
						$global_genes->{$patient->name}->{$g->id}->{no_compound} ++;
					}
					unless($father){
							$global_genes->{$patient->name}->{$g->id}->{no_compound} ++;
					}
				# $project->buffer->disconnect();
				}
				unless (exists $hgenes->{$g->id}){
					$hgenes->{score_mother} = -100;
					$hgenes->{score_father} = -100;
					$hgenes->{score_biallelic} = -10;
					$hgenes->{id} = $g->id;
			}
			$hgenes->{score_variant}->{$v->id} = $score_variant;
			die() if exists $hgenes->{variant}->{$v->id};
			
			$hgenes->{variant}->{$v->id}->{score} = $score_variant;
			$hgenes->{pathogenic} ++  if  $v->isDM_for_gene($g) or $v->isClinvarPathogenic_for_gene($g);
			$hgenes->{clinvar_hgmd} ++ if  $v->hgmd or $v->clinvar;
			
			$hgenes->{cnv_del} ++ if  $v->isLargeDeletion();
			$hgenes->{cnv_dup} ++ if  $v->isLargeDuplication();
#			warn $v->id." :: ". $v->isDenovoTransmission($fam,$patient)." :: ".$v->getGnomadAC." denovo_rare : ".$patient->name.":".$hgenes->{denovo_rare}." ".$global_genes->{$patient->name}->{$g->id}->{denovo_rare} if $debug;
#			warn "\t\t ------- ".$hgenes->{cnv_dup}." ".$hgenes->{cnv_del} if $debug;

			$hgenes->{denovo_rare} ++ if $v->getGnomadAC< 10 && $v->isDenovoTransmission($fam,$patient) ;#&& $v->effectImpact($g) eq "moderate";
			$global_genes->{$patient->name}->{$g->id}->{denovo_rare} ++ if $v->getGnomadAC< 10 && $v->isDenovoTransmission($fam,$patient);
			$global_genes->{$patient->name}->{$g->id}->{cnv_del} ++ if  $v->isLargeDeletion();;
			$global_genes->{$patient->name}->{$g->id}->{cnv_dup} ++ if  $v->isLargeDuplication();
			
			
			my $val = $v->score_validations($g);
			#my $val = 0;
			if ($val){
				$hgenes->{validations} = $val->{validation} ;#if $val->{validation} > $hgenes->{$g->id}->{validations};
			}
			 
			push (@{$hgenes->{variants}}, $v->id);
			my $debug;
			#$debug =1 if $v->id eq "5_289406_A_T"; 
				
			if ($patient->isChild && $patient->getFamily->isTrio()){
				
				if  ($v->isUniparentalDisomyTransmission($fam,$patient)){
					$hgenes->{score_biallelic} = $score_variant;
					$hgenes->{biallelic} ++;
					$hgenes->{variant}->{$v->id}->{biallelic} ++;
				}
				
				if  ($v->isMosaicTransmission($fam,$patient)){
					warn "mos" if $debug;
					$hgenes->{score_biallelic} = $score_variant;
					$hgenes->{biallelic} ++;
					$hgenes->{variant}->{$v->id}->{biallelic} ++;
				}
				elsif  ($v->isDominantTransmission($fam,$patient)){
					$hgenes->{score_biallelic} = $score_variant;
					$hgenes->{biallelic} ++;
					$hgenes->{variant}->{$v->id}->{biallelic} ++;
				}
				 elsif  ($v->isMotherTransmission($fam,$patient)){
				 #	warn "mot" if $debug;
				 	$hgenes->{score_mother} = $score_variant;#>=$hgenes->{score_mother}?$score_variant:$hgenes->->{score_mother});
				 	$hgenes->{mother} ++;
				 	$hgenes->{variant}->{$v->id}->{mother} ++;
				 	#push(@)
				 	if ($global_genes->{$patient->name}->{$g->id}->{mother}->{score} < $score_variant  ){
				 		$global_genes->{$patient->name}->{$g->id}->{mother}->{id} = $v->id;
				 		$global_genes->{$patient->name}->{$g->id}->{mother}->{score} = $score_variant;
				 		$global_genes->{$patient->name}->{$g->id}->{mother}->{vector_id} .= $v->{vector_id};
				 		$global_genes->{$patient->name}->{$g->id}->{mother}->{global_vector_id} .= $v->{global_vector_id};
				 	}
				 	
				 }
				 elsif  ($v->isFatherTransmission($fam,$patient)){
				 #	warn "father ".$v->id if $debug;
				 	$hgenes->{score_father} = $score_variant;#>=$hgenes->{score_father}?$score_variant:$hgenes->{score_father});
				 	$hgenes->{father} ++;
				 	$hgenes->{variant}->{$v->id}->{father} ++;
				 	if ($global_genes->{$patient->name}->{$g->id}->{father}->{score} < $score_variant){
				 		$global_genes->{$patient->name}->{$g->id}->{father}->{id} = $v->id;
				 		$global_genes->{$patient->name}->{$g->id}->{father}->{score} = $score_variant;
				 		$global_genes->{$patient->name}->{$g->id}->{father}->{vector_id} .= $v->{vector_id};
				 		$global_genes->{$patient->name}->{$g->id}->{father}->{global_vector_id} .= $v->{global_vector_id};
				 	}
				 }
				 
				 else {
				 	warn "else" if $debug;
				 	$hgenes->{variant}->{$v->id}->{biallelic} ++;
				 	$hgenes->{score_biallelic} = $score_variant;#>=$hgenes->{score_biallelic}?$score_variant:$hgenes->{score_biallelic});
				 }
				
			}
			else {
				$hgenes->{biallelic} ++;	
				$hgenes->{score_biallelic} = $score_variant;#>=$hgenes->{score_biallelic}?$score_variant:$hgenes->{score_biallelic});
				$hgenes->{variant}->{$v->id}->{biallelic} ++;
				#next unless  $v->annotation()->{$g->id}->{mask} & $maskcoding;
			}
			#warn $hgenes->{variant}->{$v->id}->{biallelic}."---------------".$patient->name." ".$hgenes->{id} if $debug;
			#warn $hgenes->{variant}->{$v->id}->{father}."---------------" if $debug;
			my $mask = $v->annotation()->{$g->id}->{mask};
			$hgenes->{mask} = $mask;
			push(@{$vmask->{$mask}},$g->id);
			 		
			push(@$agenes,$hgenes);
		}#end Genes 
		 	$h->{patient} = $patient->name;
		 	$h->{id1} = $v->id;
		 	$h->{type} = "variants";
		 	$no->put($v->{global_vector_id}."-hashvariants@".$patient->name,$h);
		 	
		 	#$no->put($v->id."-html@".$patient->name,$h->{html});
		 	
		 	$no->put($v->{global_vector_id}."-mask@".$patient->name,{value=>$vmask,patient=>$patient->name,id1=>$v->vector_id,type=>"mask"}) if $vmask;
		 	$no->put($v->{global_vector_id}."-genes@".$patient->name,{array=>$agenes,patient=>$patient->name,id1=>$v->{global_vector_id},type=>"genes"}) if $agenes;
		 	$no->put($v->{global_vector_id}."@".$patient->name,{id=>$v->id,id1=>$v->{global_vector_id},patient=>$patient->name,type=>"variants",});
		 	push(@{$zh->{$patient->name}},$v->{global_vector_id});
		 	
		}
		
	}#end for each patient
	foreach my $p (keys %$zh){
		$no->put($p,$zh->{$p});
	}
	
	warn "stop 1";
	$no->close();
	warn "stop 2";
	my $h= {};
	
	$h->{index} = $index;
	
	
	$h->{genes} = $global_genes;
	$h->{ids} = $zh;
	#foreach my $g (keys %$global_genes){
	#	foreach my $p (keys %{$global_genes->{$g} } ){
	#		my $hdata = $global_genes->{$g}->{$p};
	#		$no->put($g."-genes_global@".$p,{value=>$hdata,patient=>$p,id1=>$g,type=>"genes_global"}) ;
	#	}
	#}
	warn "finish";
	$pm->finish(0,$h);
	
	#die();
}
warn "++++ end children";
$pm->wait_all_children();
#die();
######## END STEP 1 
######### NOW SAVING DATA
##########################



$project->buffer->dbh_reconnect();



if (keys %$proc) {
	warn "argggggggg ";
	foreach my $f (@files){
			unlink $dir_tmp."/".$f;
			unlink $dir_tmp."/".$f."_index";
	}
}

die("problem") if keys %$proc;

warn "create";

$pm = new Parallel::ForkManager($fork);


$pm->run_on_finish (
    sub { 
    	my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$h)=@_;
    	unless (defined($h) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
		warn "end " .$h->{index};
		delete $proc->{$h->{index}};
		
    }
    
    );
   warn "file 1";
	foreach my $patient (@{$project->getPatients}){	
	my $no_p = $chr->lmdb_polyviewer_variants_genes($patient,"c");
	$no_p->put("date",time);
	$chr->lmdb_polyviewer_variants_genes($patient,"close");
	
	my $no_v = $chr->lmdb_polyviewer_variants($patient,"c");
	$no_v->put("date",time);
	$chr->lmdb_polyviewer_variants($patient,"close");
	$no_v->close();
	$no_p->close();
	}


 $project->buffer->disconnect();
 #my $buffer = GBuffer->new();
#my $project = $buffer->newProjectCache( -name => $project_name, -cache => '1', -typeFilters => 'familial' );
 
 
 warn "ile2";

		foreach my $patient (@{$project->getPatients}){
			$proc->{$patient->name} ++;
			
			my $pid = $pm->start and next;
			 $project->buffer->dbh_reconnect();
			my $no_p = $chr->lmdb_polyviewer_variants_genes($patient,"w");
			my $no_v = $chr->lmdb_polyviewer_variants($patient,"w");
			foreach my $f (@files){
				my $noin = GenBoNoSqlLmdb->new(dir=>$dir_tmp,mode=>"r",name=>$f,is_compress=>1);
				save_patient($chr,$noin,$patient,$no_p,$no_v);
				
				$noin->close();
				$noin = undef;
				
			}
			my $h= {};
		$h->{index} = $patient->name;	
		 $chr->lmdb_polyviewer_variants_genes($patient,"close");
		$chr->lmdb_polyviewer_variants($patient,"close");
		$no_p->close;
		$no_v->close;
		
		$pm->finish(0,$h);
	#die();

			
		}
		$pm->wait_all_children();
	 $project->buffer->dbh_reconnect();

#$pm->wait_all_children();

warn "end step 2";

if (keys %$proc) {
	warn "argggggggg ";
	foreach my $f (@files){
			unlink $dir_tmp."/".$f;
			unlink $dir_tmp."/".$f."_index";
	}
	foreach my $patient (@{$project->getPatients}){
		my $f = $chr->lmdb_polyviewer_variants($patient,"r")->filename;
		unlink $f if -e $f;
		 $f = $chr->lmdb_polyviewer_variants_genes($patient,"r")->filename;
		unlink $f if -e $f;
		
	}
	die("problem phase 2");
}
warn "start phase 3 ";
	foreach my $p (keys %$hgenes){
		$chr->lmdb_polyviewer_genes($p,"c")->put("date",time);
		foreach my $g (keys %{$hgenes->{$p} } ){
			$hgenes->{$p}->{$g}->{penality} = 0;
			$hgenes->{$p}->{$g}->{penality} = $hgenes->{$p}->{$g}->{denovo_rare}* -0.3  if $hgenes->{$p}->{$g}->{denovo_rare} > 2;
			$chr->lmdb_polyviewer_genes($p)->put($g,$hgenes->{$p}->{$g});
		}
		$chr->lmdb_polyviewer_genes($p,"close");
	}
foreach my $f (@files){
	unlink $dir_tmp."/".$f;
	unlink $dir_tmp."/".$f."_index";
}


my $dir = $chr->lmdb_cache_dir();
warn "date >$dir/lmdb.ok && chmod a+rw $dir/lmdb.ok";
system ("date >$dir/lmdb.ok && chmod a+rw $dir/lmdb.ok");
warn "$dir/lmdb.ok";
$dir_tmp = undef;
return ("$dir/lmdb.ok");
}

sub save_patient {
	my ($chr,$noin,$patient,$no_p,$no_v) = @_;
	
	
	#my $noin = GenBoNoSqlLmdb->new(dir=>$dir_tmp,mode=>"r",name=>$f,is_compress=>1);
		my $list = $noin->get($patient->name);
	
			return  unless $list;
#			warn scalar (@$list);
			foreach my $id (@$list){
					#variatns 
					my $value =  $noin->get($id."@".$patient->name);
					my $a = $value->{id1};
					delete $value->{id1};
					delete $value->{patient};
					delete $value->{type};
					die() unless $value;
					$no_v->put($a,$value);
					
					 $value =  $noin->get($id."-hashvariants@".$patient->name);
					 $a = $value->{id1};
					delete $value->{id1};
					delete $value->{patient};
					delete $value->{type};
					die() unless $value;
					$no_v->put($a,$value);
					
					#genes
					
					$value =  $noin->get($id."-genes@".$patient->name);
					 $a = $value->{id1};
					delete $value->{id1};
					delete $value->{patient};
					delete $value->{type};
					$no_p->put($a,$value);
					#mask
					
					#$chr->lmdb_polyviewer_mask("w")->put($a,$value->{value});
					
					
			}
#			warn $no_v->filename;
		
			#$no_v->close();
			#$no_p->close();
#			warn $no_p->nb_keys()." ".$patient->name;
}



sub listVariants {
		my($chr) =@_;
		my $vector;
		if ($project->isDiagnostic || $project->isExome ) {
			 $vector = $chr->getVectorScore("gnomad_ac_all");
		}
		else {
		if ($chr->name eq "MT"){
			 $vector = $chr->getVectorScore("gnomad_ac_all");
		}
		else {
		 $vector = $chr->getVectorScore("gnomad_ac_5000");
		$vector += $chr->getVectorScore("gnomad_ho_ac_1000");
		 #$vector += $chr->getVariantsVector();
		}
		}
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
