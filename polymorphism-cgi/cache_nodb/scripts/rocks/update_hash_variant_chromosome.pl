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
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/packages/";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Getopt::Long;
use Carp;
use GBuffer;
use Set::IntSpan::Fast::XS;
use Sys::Hostname;
use File::Temp qw/ tempfile tempdir /;
use Storable qw(dclone);
use Clone qw(clone);
require "$RealBin/../../../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
#use Try::Tiny;


my $host = hostname();


#warn "*_*_*_*_*_".$host."*_*_*_*_*_";
#use Cache_Commons;
my $fork = 1;
my $force;
my ($project_name, $chr_name, $annot_version);
my $ok_file;
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'chr_name=s'        => \$chr_name,
	'annot_version=s'    => \$annot_version,
	'force=s'  => \$force,
	'file=s'  => \$ok_file,
);

 if ($ok_file && -e $ok_file) {
 	system("rm $ok_file");
 }

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
#try {
foreach my $chr (@{$project->getChromosomes} ){
	if ($chr_name ne "all"){
		next if  $chr->name ne $chr_name;
	}
	warn "coucou";
	 $file = run_update_chr($chr->name);
}


warn "OK";
system("date > $ok_file") if $ok_file;
exit(0);
#}
#catch {
#	$buffer = undef;
#	die();
#};

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
my %htansmission = {};
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
		%htansmission = (%{$h->{transmission}},%htansmission);
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

my $no_var = $chr->get_rocks_variations("r");
my $ranges = $no_var->ranges($fork);
$no_var->close;
$no_var = undef;


my @headers = ( "varsome", "igv", "alamut", "var_name", "trio", "gnomad", "deja_vu" );

my $part =0;
my $index = 0;
my @files ;

	my $javascript_id = time + int(rand(10000));
	$project->preload();
	
	$project->disconnect();
	foreach my $range (@$ranges){	
	$part ++;
	warn $part;
	my $f = $chr->name."_".$part."_".$range->[0]."_".$range->[-1];
	push(@files,$f);
	$index ++;
	$proc->{$index} ++;
	my $pid = $pm->start and next;

	warn "start ".$proc->{$index};
	my $nb =0;
	
	#$project->buffer->{dbh} = "-";
	
	my $no  = GenBoNoSqlLmdb->new(dir=>$dir_tmp,mode=>"c",name=>$f,is_compress=>1);
	my $array;

	my $nb =1 ;
	my $global_genes;
	my $zh;
	#sleep(2);
	my $chr = $project->getChromosome($chr_name);
	my $transmission;
	my $no_var = $chr->get_rocks_variations("r");
	my $max = $range->[1] - $range->[0];
	my @tmp;
	my $prints;
	my $vtime = time;
	for (my $aid=$range->[0];$aid <= $range->[1];$aid++){
		push(@tmp,$aid);
	#while(my $v = $project->nextVariant){
		my $v = $no_var->get_index($aid);
		#next if $v->isCnv;
		$v->{buffer} = $buffer;
		$v->{project} = $project;
		#warn $v->name;
		$nb++;
		if ($nb%3000 ==0){
		#	$no->close();
		#	$no  = GenBoNoSqlLmdb->new(dir=>$dir_tmp,mode=>"w",name=>$f,is_compress=>1);
			print "$nb/$max \n" ;
			warn "$nb/$max ".abs(time - $vtime); 
			$vtime = time;
		}
		my $hvariant = {};
		# update_variant_editor::construct_hash_variant_global ( $project, $v,undef,1);
		
	#	warn $v->id;
		my $vmask;
		foreach my $patient (@{$v->getPatients}) {
			
			my $fam = $patient->getFamily();
			#my $h = {};
		 	#update_variant_editor::construct_hash_variant_patient( $project, $v,$patient,$h);
		 	#$h->{vector_id} = $v->vector_id;
		 	#update_variant_editor::vvarsome($h,$patient);
			my $mother = $patient->getFamily->getMother();
			my $father = $patient->getFamily->getFather();
		  	$transmission->{ $v->global_vector_id}->{$patient->id} = $v->getTransmissionModelType($patient->getFamily(),$patient);
		 	my $agenes =[];
		 	my $score_genes;
		 	foreach my $g (@{$v->getGenes}){
		 		my $debug;
		 		$debug = 1 if $g->id eq "ENSG00000167522_16";
		 		my $hgenes;
				$javascript_id ++;
				my $score_variant = $v->scaledScoreVariant($g,$patient);
				my $parent_score = 0;
				$parent_score = $score_variant  if $score_variant >0 &&  $v->isConsequence($g,"medium");
				my $nid = "s:".$v->{vector_id}.":".$g->id."@".$patient->name;
				#$h->{scaled_score}->{$g->id} = $score_variant;
				
				$no->put($nid,{id=>$v->id,score=>$score_variant,parent_score=>$parent_score,id1=>"s:".$v->{vector_id}.":".$g->id,patient=>$patient->name,type=>"genes",pp=>$patient->name});
		
				
				unless (exists $global_genes->{$patient->name}->{$g->id}) {
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
					$global_genes->{$patient->name}->{$g->id}->{vector} = 0;
					# $g->getCurrentVector();;
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
			$hgenes->{vector_ids}->{$v->id} = $v->vector_id;
			$hgenes->{variant}->{$v->id}->{score} = $score_variant;
			$hgenes->{pathogenic} ++  if  $v->isDM_for_gene($g) or $v->is_clinvar_pathogenic_for_gene($g);
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
			$hgenes->{vector_ids}->{$v->id} = $v->vector_id;
			my $debug;
			#$debug =1 if $v->id eq "5_289406_A_T"; 
				
			if ($patient->isChild && $patient->getFamily->isTrio()){
				
				if  ($v->isUniparentalDisomyTransmission($fam,$patient)){
					$hgenes->{score_biallelic} = $score_variant;
					$hgenes->{biallelic} ++;
					$hgenes->{variant}->{$v->id}->{biallelic} ++;
				}
				
				if  ($v->isMosaicTransmission($fam,$patient)){
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
				 	$hgenes->{vector_ids}->{$v->id} = $v->vector_id;
				 	#push(@)
				 	if ($global_genes->{$patient->name}->{$g->id}->{mother}->{score} < $score_variant  ){
				 		$global_genes->{$patient->name}->{$g->id}->{mother}->{id} = $v->id;
				 		$global_genes->{$patient->name}->{$g->id}->{mother}->{vector_id} = $v->vector_id;
				 		$global_genes->{$patient->name}->{$g->id}->{mother}->{score} = $score_variant;
				 	}
				 	
				 }
				 elsif  ($v->isFatherTransmission($fam,$patient)){
				 #	warn "father ".$v->id if $debug;
				 	$hgenes->{father} ++;
				 	$hgenes->{variant}->{$v->id}->{father} ++;
				 	$hgenes->{vector_ids}->{$v->id} = $v->vector_id;
				 	if ($global_genes->{$patient->name}->{$g->id}->{father}->{score} < $score_variant){
				 		$global_genes->{$patient->name}->{$g->id}->{father}->{id} = $v->id;
				 		$global_genes->{$patient->name}->{$g->id}->{father}->{vector_id} = $v->vector_id;
				 		$global_genes->{$patient->name}->{$g->id}->{father}->{score} = $score_variant;
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
			}
			my $mask = $v->annotation()->{$g->id}->{mask};
			$hgenes->{mask} = $mask;
			push(@{$vmask->{$mask}},$g->id);
			 		
			push(@$agenes,$hgenes);
		}#end Genes 
		 	$no->put($v->global_vector_id."-mask@".$patient->name,{value=>$vmask,patient=>$patient->name,id1=>$v->vector_id,type=>"mask"}) if $vmask;
		 	$no->put($v->global_vector_id."-genes@".$patient->name,{array=>$agenes,patient=>$patient->name,id1=>$v->global_vector_id,type=>"genes"}) if $agenes;
		  	delete $v->{genes};
		 	$no->put($v->global_vector_id."@".$patient->name,{id=>$v->id,id1=>$v->global_vector_id,patient=>$patient->name,type=>"variants",});
		 	push(@{$zh->{$patient->name}},$v->global_vector_id);
		 	
		}
		delete $v->{buffer};
		delete $v->{project};
	}#end for each patient
	foreach my $p (keys %$zh){
		$no->put($p,$zh->{$p});
	}
	$no->put("list",\@tmp);
	$no->close();
	my $h= {};
	$h->{transmission} = $transmission;
	$h->{index} = $index;
	$h->{genes} = $global_genes;
	$h->{ids} = $zh;

	$pm->finish(0,$h);
	
	#die();
}

$pm->wait_all_children();
warn "++++ end children";
warn "OKkkkkk";
#die();
######## END STEP 1 
######### NOW SAVING DATA
##########################






if (keys %$proc) {
	warn "argggggggg ";
	foreach my $f (@files){
			unlink $dir_tmp."/".$f;
			unlink $dir_tmp."/".$f."_index";
	}
}

die("problem") if keys %$proc;

warn "create";

 
 warn "************************************ile2";
 
  warn "ile2";



my $final_polyviewer = GenBoNoSqlRocks->new(dir=>$project->rocks_pipeline_directory("polyviewer_raw"),mode=>"w",name=>$chr->name);
#$hh->{model} = "-";#$vh->getTransmissionModelType($p->getFamily(),$p);
foreach my $k (keys %htansmission){
	
	my $a = $final_polyviewer->get($k);
	# 'patients_calling' => {
     #                                    '55661' => {
	foreach my $pid (keys %{$htansmission{$k}} ){
		die() unless exists $a->{patients_calling}->{$pid};
		$a->{patients_calling}->{$pid}->{model} = $htansmission{$k}->{$pid};
	}
	$final_polyviewer->put_batch($k,$a);
}
$final_polyviewer->write_batch();
$final_polyviewer->close();
warn " end transmission ";
my $dir_pipeline = $project->rocks_pipeline_directory("patients");

my $no_p = {};

		foreach my $patient (@{$project->getPatients}){
			 $no_p->{$patient->id} =  GenBoNoSqlRocks->new(dir=>$dir_pipeline."/".$patient->name,mode=>"c",name=>$chr->name);
			$proc->{$patient->name} ++;
			
			 $project->buffer->dbh_reconnect();
			my $no_v = undef;
			foreach my $f (@files){
				my $noin = GenBoNoSqlLmdb->new(dir=>$dir_tmp,mode=>"r",name=>$f,is_compress=>1);
				save_patient($chr,$noin,$patient, $no_p->{$patient->id});
				
				$noin->close();
				$noin = undef;
				
			}
			my $h= {};
		$h->{index} = $patient->name;	
		}
	

warn "--------------------";
warn "start phase 3 ";

	foreach my $p (keys %$hgenes){
		my $patient = $project->getPatient($p);
		
		foreach my $g (keys %{$hgenes->{$p} } ){
			$hgenes->{$p}->{$g}->{penality} = 0;
			$hgenes->{$p}->{$g}->{penality} = $hgenes->{$p}->{$g}->{denovo_rare}* -0.3  if $hgenes->{$p}->{$g}->{denovo_rare} > 2;
			 $no_p->{$patient->id}->put_batch($g,$hgenes->{$p}->{$g});
		}
		 $no_p->{$patient->id}->write_batch();
		 $no_p->{$patient->id}->close();
	}


 

foreach my $f (@files){
	unlink $dir_tmp."/".$f;
	unlink $dir_tmp."/".$f."_index";
}


my $dir = $project->rocks_cache_dir();
my $file = "$dir/lmdb".$chr->name.".ok";

warn "date >$dir/lmdb.ok && chmod a+rw $dir/lmdb.ok";
system ("date >$dir/$file && chmod a+rw $dir/$file");
warn "$dir/lmdb.ok";
$dir_tmp = undef;
return ("$dir/$file");
}





sub save_patient {
	my ($chr,$noin,$patient,$no_p,$no_v) = @_;
	
		my $list = $noin->get($patient->name);
			return  unless $list;
			foreach my $id (@$list){
				
					my $value =  $noin->get($id."-genes@".$patient->name);
					$a = $value->{id1};
					$no_p->put_batch($a,$value);
				
					
					
			}
			$no_p->write_batch;

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
		 $vector = $chr->getVectorScore("gnomad_ac_all");
		 #$vector += $chr->getVariantsVector();
		}
		}
		# $vector = $chr->getVectorScore("gnomad_ac_all");
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
