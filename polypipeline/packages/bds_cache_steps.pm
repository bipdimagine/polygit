package bds_cache_steps;
use strict;
use FindBin qw($Bin);
use lib "$Bin";
#use root_steps;
use Moose;  

use job_bds;
use sample;
use Data::Dumper;
extends (qw(bds_root));
my $bin_dev = qq{$Bin/scripts/scripts_pipeline/};

	
sub global_infos   {
	my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $project = $self->patient();	
	my $project_name = $project->name();
	my $ppn =1;
	warn $project->name();;
	my $fileout = $project->getCacheBitVectorDir()."/global_infos.freeze";
	my $log_error = $project->getCacheBitVectorDir()."/log/global_infos.log";
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/scripts/cache_global_infos.pl -project=$project_name ";
 	my $type = "global_infos";
	 my $stepname = $self->patient->name."@".$type;
	 my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}

	return ($fileout);
}
sub define_ppn {
	my ($self,$kid) = @_;
	my $ppn = $self->nproc;

	return 20;

	if (between($kid, 8, 15)) {
    	$ppn = 10;
   }
   if (between($kid, 15, 30)) {
    	$ppn = 5;
   }
   if ($kid == 23){
   	$ppn = 20;
   }
	
	if ($self->nocluster){
		$ppn = 10;
	}
	return 20;
}

sub update_variants  {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $chr = $self->patient();
	my $chr_name = $chr->name();
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->define_ppn($chr->karyotypeId);
	$ppn = 1 if $ppn < 1;
	$filein ="";
	my $fileout = $project->lmdb_cache_variations_dir()."/".$chr_name.'.update';
	my $cmd = "/usr/bin/perl $Bin/scripts/scripts_cache_update/update_variants_objects.pl -project=$project_name -chr=$chr_name -fork=$ppn";
	my $type = "update_variants";
	my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
		$job_bds->skip();
	}
	return ($fileout);
}

sub update_chromosomes  {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $chr = $self->patient();
	my $chr_name = $chr->name();
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->define_ppn($chr->karyotypeId);
	$ppn = 1 if $ppn < 1;
	$filein ="";
	my $fileout = $project->lmdb_cache_variations_dir()."/".$chr_name.'.chr_update';
	my $cmd = "/usr/bin/perl $Bin/scripts/scripts_cache_update/update_chromosomes_vectors.pl -project=$project_name -chr=$chr_name -fork=$ppn && date > $fileout";
	my $type = "update_chromosomes";
	my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
		$job_bds->skip();
	}
	return ($fileout);
	
}

sub update_genes  {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $chr = $self->patient();
	my $chr_name = $chr->name();
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->define_ppn($chr->karyotypeId);
	$ppn = 1 if $ppn < 1;
#	$filein ="";
	my $fileout = $project->lmdb_cache_variations_dir()."/".$chr_name.'.genes_update';
	my $cmd = "/usr/bin/perl $Bin/scripts/scripts_cache_update/update_genes_vectors.pl -project=$project_name -chr=$chr_name -fork=$ppn && date > $fileout";
	my $type = "update_genes";
	my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
		$job_bds->skip();
	}
	return ($fileout);
	
}

sub store_ids {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $chr = $self->patient();
	my $chr_name = $chr->name();
	
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->define_ppn($chr->karyotypeId);
	#warn $chr->karyotypeId." ".$ppn;
	
	$ppn = 1 if $ppn < 1;
	$ppn = 15;
	$filein ="";
	my $fileout = $project->getCacheBitVectorDir()."/lmdb_cache/".$chr_name.".dv.freeze";
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/scripts/cache_store_ids.pl -project=$project_name -chr=$chr_name -fork=$ppn";
	
	 my $type = "cache_store_ids";
	 my $stepname = $self->patient->name."@".$type;
	 my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>20,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
}

sub store_rna_junction_ids  {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $chr = $self->patient();
	my $chr_name = $chr->name();
	
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->define_ppn($chr->karyotypeId);
	#warn $chr->karyotypeId." ".$ppn;
	
	$ppn = 1 if $ppn < 1;
	$ppn = 1;
	$filein ="";
	my $fileout = $project->getCacheBitVectorDir()."/lmdb_cache/".$chr_name.".dv.freeze";
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/scripts/cache_rna_junctions_store_ids.pl -project=$project_name -chr=$chr_name -fork=$ppn";
	
	 my $type = "cache_store_junctions_ids";
	 my $stepname = $self->patient->name."@".$type;
	 my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>20,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
}

sub store_annotations  {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $chr = $self->patient();
	
	my $chr_name = $chr->name();
	
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->define_ppn($chr->karyotypeId);
	$ppn = 1 if $ppn < 1;
	
	my $fileout = $project->getCacheBitVectorDir()."/lmdb_cache/".$chr_name."/genes_index";
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/scripts/cache_store_annotations.pl -project=$project_name -chr=$chr_name -fork=$ppn";
	
	my $log_error = $project->getCacheBitVectorDir()."/log/store_anotations_$chr_name.log";
	 my $type = "store_annotations";
	 my $stepname = $self->patient->name."@".$type;
	 my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
	
}

sub polyviewer  {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $chr = $self->patient();
	my $chr_name = $chr->name();
	
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->define_ppn($chr->karyotypeId);
	$ppn = 1 if $ppn < 1;
	my $fileout = $chr->lmdb_cache_dir()."/lmdb.ok";

	#my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/scripts/update_vector_score.pl -project=$project_name -chr=$chr_name -fork=$ppn";
	# $cmd .= " && /usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/scripts/update_hash_variant_chromosome.pl -project=$project_name -chr=$chr_name -fork=$ppn && echo ok && sleep 2";
	
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/scripts/polyviewer.pl -project=$project_name -chr=$chr_name -fork=$ppn";
	
	
	my $log_error = $project->getCacheBitVectorDir()."/log/store_anotations_$chr_name.log";
	 my $type = "polyviewer";
	 my $stepname = $self->patient->name."@".$type;
	 my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
}

sub  diagHash {
	my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $chr = $self->patient();
	my $chr_name = $chr->name();
	
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->define_ppn($chr->karyotypeId);
	$ppn = 1 if $ppn < 1;
	my $fileout = $project->getCacheBitVectorDir()."/lmdb_cache/".$chr_name."/hashes_variants";
	
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/scripts/update_hash_variant_chromosome.pl -project=$project_name -chr=$chr_name -fork=$ppn";
	
	my $log_error = $project->getCacheBitVectorDir()."/log/store_anotations_$chr_name.log";
	 my $type = "diagHash";
	 my $stepname = $self->patient->name."@".$type;
	 my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
}

sub update_score  {
	my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $chr = $self->patient();
	
	my $chr_name = $chr->name();
	
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->define_ppn($chr->karyotypeId);
	$ppn = 1 if $ppn < 1;
	my $fileout =   $project->dir_lmdb_score_impact()."/".$chr->name;
	
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/scripts/update_vector_score.pl -project=$project_name -chr=$chr_name -fork=$ppn";
	
	
	my $log_error = $project->getCacheBitVectorDir()."/log/store_anotations_$chr_name.log";
	 my $type = "update_score";
	 my $stepname = $self->patient->name."@".$type;
	 my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
		my $no = $chr->lmdb_variations("r");
			
		 		$job_bds->skip() if $no->get("gnomad_ac_10") && $no->get("gnomad_ac_100");
	}
	return ($fileout);
}


sub update_coverage {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $chr = $self->patient();
	
	my $chr_name = $chr->name();
	
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->define_ppn($chr->karyotypeId);
	$ppn = 1 if $ppn < 1;
	
	my $fileout =  $project->getCoverageDir()."/lmdb_images_uri/".$chr_name;
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/scripts/update_coverage.pl -project=$project_name -chr=$chr_name -fork=$ppn";
	
	
	my $log_error = $project->getCacheBitVectorDir()."/log/store_anotations_$chr_name.log";
	 my $type = "update_coverage";
	 my $stepname = $self->patient->name."@".$type;
	 my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
}



sub transcripts_coverage {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $p = $self->patient();
	
	my $patient_name = $p->name();
	
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = 10;
	$ppn = 1 if $ppn < 1;
	my $no2 = $p->getTranscriptsCoverageDepth("w");
	my $fileout =  $no2->filename();
	 my $bin_dev = $self->script_dir();
	
	my $cmd = "/usr/bin/perl $bin_dev/transcripts/transcripts_coverage.pl  -project=$project_name -patient=$patient_name -fork=$ppn";
	
	my $log_error = $project->getCacheBitVectorDir()."/log/transcript_coverage__$patient_name.log";
	 my $type = "transcripts_coverage";
	 my $stepname = $self->patient->name."@".$type;
	 my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
}

sub identito_vigilence {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();
	my $project_name = $project->name();
	
	
	my $no = $self->patient()->getGenesDude();
	my $fileout = $no->filename();
	my $ppn =$self->nproc;
	$ppn = int($self->nproc/2) if $self->nocluster;
	
	
	### Récupération de la version de samtools
	
	my $bin_dev = $self->script_dir;
	
	my $cmd = qq{perl $bin_dev/identito_vigilence.pl -patient=$name  -fork=$ppn  -project=$project_name };
	my $type = "genes_dude";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	my $no_cache = $self->patient()->get_lmdb_cache("r");
	if ($self->unforce() &&  -e $no_cache->filename ){
		 		$job_bds->skip();
	}
	$no->close();
	return ($fileout);
}

sub html_cache_polyviewer {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();
	my $project_name = $project->name();
	
	
	my $no = $self->patient()->getGenesDude();
	my $fileout = $no->filename();
	my $ppn =$self->nproc;
	$ppn = int($self->nproc/2) if $self->nocluster;
	
	
	### Récupération de la version de samtools
	
	my $bin_dev = $self->script_dir;

	#my $cmd = qq{/$Bin/variations_editor.pl  $args allele_quality=- annot="splicing+essential_splicing+nonsynonymous+stop+phase+maturemirna+frameshift+non-frameshift+predicted_splice_site" denovo=1 dv=3 dv_ho=2 edit_mode=1 in_this_run=6 keep_pathogenic=1 never=1 patients=$name $ph project=$project_name recessive=1 report_mode=1 strict_denovo=1 user_name= xor=1 >/dev/null && date > $totofile};
	my $no_cache = $self->patient()->get_lmdb_cache("r");
	$fileout = $no_cache->filename.".polyviewer";
	my $cmd = qq{perl $bin_dev/../scripts_cache/polyviewer/variations_editor_cache.pl -patient=$name  -project=$project_name && touch $fileout};
	my $type = "html_polyviewer";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() &&  -e $no_cache->filename && -e $fileout ){
		 		$job_bds->skip();
	}
	$no->close();
	return ($fileout);
}

sub html_cache_polycyto{
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();
	my $project_name = $project->name();
	
	
	my $no = $self->patient()->getGenesDude();
	my $fileout = $no->filename();
	my $ppn =$self->nproc;
	$ppn = int($self->nproc/2) if $self->nocluster;
	
	
	### Récupération de la version de samtools
	
	my $bin_dev = $self->script_dir;
	$ppn = 2;
	#my $cmd = qq{/$Bin/variations_editor.pl  $args allele_quality=- annot="splicing+essential_splicing+nonsynonymous+stop+phase+maturemirna+frameshift+non-frameshift+predicted_splice_site" denovo=1 dv=3 dv_ho=2 edit_mode=1 in_this_run=6 keep_pathogenic=1 never=1 patients=$name $ph project=$project_name recessive=1 report_mode=1 strict_denovo=1 user_name= xor=1 >/dev/null && date > $totofile};
	my $no_cache = $self->patient()->get_lmdb_cache("r");
	$fileout = $no_cache->filename.".polycyto";
	my $cmd = qq{perl $bin_dev/../scripts_cache/polycyto/polycyto_cache_html.pl -patient=$name -project=$project_name && touch $fileout};
	my $type = "html_polycyto";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() &&  -e $no_cache->filename && -e $fileout){
		 		$job_bds->skip();
	}
	
	return ($fileout);
}
sub genes_dude {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();
	my $project_name = $project->name();
	unless ($filein){
 		my $no = $self->patient()->getTranscriptsDude();
	 	$filein = $no->filename();
	}
	
	my $no = $self->patient()->getGenesDude();
	my $fileout = $no->filename();
	my $ppn =$self->nproc;
	$ppn = int($self->nproc/2) if $self->nocluster;
	
	
	### Récupération de la version de samtools
	
	my $bin_dev = $self->script_dir;
	
	#my $cmd = qq{perl $bin_dev/transcripts/genes_level_dude.pl -patient=$name  -fork=$ppn  -project=$project_name && perl $bin_dev/transcripts/update_level.pl -patient=$name  -fork=$ppn  -project=$project_name};
	my $cmd = qq{perl $bin_dev/transcripts/genes_level_dude.pl -patient=$name  -fork=$ppn  -project=$project_name};
	my $type = "genes_dude";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() &&  $no->nb_keys >10 ){
		 		$job_bds->skip();
	}
	$no->close();
	return ($fileout);
}

sub polydude {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();
	my $project_name = $project->name();
	unless ($filein){
 		my $no = $self->patient()->getTranscriptsDude();
	 	$filein = $no->filename();
	}
	
	my $no = $self->patient()->getGenesDude();
	my $fileout = $no->filename()."update_log";
	my $ppn =$self->nproc;
	$ppn = int($self->nproc/2) if $self->nocluster;
	
	
	### Récupération de la version de samtools
	
	my $bin_dev = $self->script_dir;
	
	my $cmd = qq{perl $bin_dev/transcripts/update_level.pl -patient=$name  -fork=$ppn  -project=$project_name && touch $fileout};
	my $type = "polydude";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce()  && -e $fileout ){
		 		$job_bds->skip();
	}
	$no->close();
	return ($fileout);
}

#sub level_dude (Str :$filein!) {
#	
#	my $patient = $self->patient();
#	my $patient_name = $patient->name();
#	
#	my $project = $self->project();
#	my $project_name = $project->name();
#	#warn $chr->karyotypeId." ".$ppn;
#	my $ppn = 20;
#	my $no3 = $patient->getGenesDude("r");
#	my $fileout = $no3->filename.".done"; 
#	my $cmd = qq{perl $bin_dev/transcripts/update_level.pl -patient=$patient_name  -fork=$ppn  -project=$project_name && $fileout };
#	
#	 my $type = "level_dude";
#	 my $stepname = $self->patient->name."@".$type;
#	 my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>20,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
#	$self->current_sample->add_job({job=>$job_bds});
#	
#	if ($self->unforce() && -e $fileout){
#			if($no3->get($no3->get("low_update_date") ) ) {
#		 			$job_bds->skip();
#			}
#	}
#	return ($filein);
#}

sub dude_bed {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();
	my $project_name = $project->name();
 	#warn $filein;
 	my $no = $self->patient()->getTranscriptsDude();
	 $filein = $no->filename();
	my $dir_out= $project->getVariationsDir("dude");
	#$filein = $dir_out."/".$name.".dude.lid.gz";
	#my $no = $project->noSqlCnvs("r");
	#$filein =  $no->filename;
	my $fileout = $dir_out."/".$name.".dude.lid.gz";
	
	my $ppn =int($self->nproc/4);
	$ppn = int($self->nproc/2) if $self->nocluster;
	
	
	### Récupération de la version de samtools
	
	my $bin_dev = $self->script_dir;
	
	my $cmd = qq{perl $bin_dev/transcripts/dude_bed.pl -patient=$name  -fork=$ppn  -project=$project_name };
	my $type = "dude_bed";
	 my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() &&  -e $fileout ){
		 		$job_bds->skip();
	}
	$no->close();
	return ($fileout);
}

sub sashimi_plots {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();
	my $project_name = $project->name();
	my $log_error = $project->project_path."/align/sashimi_plots/sashimi_$name.log";
	my $fileout = $log_error.'.ok';
	my $ppn = $self->nproc;
	$ppn = int($self->nproc/2) if $self->nocluster;
	my $bin_dev = $self->script_dir;
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/scripts/cache_rna_sashimi_plot.pl -project=$project_name -patient=$name -fork=$ppn -fileout=$fileout";
	my $type = "sashimi_plots";
	my $stepname = $self->patient->name."@".$type;
	$ppn =20;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() &&  -e $fileout){
		$job_bds->skip();
	}
	return ($fileout);
}

sub transcripts_dude {
	my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();
	my $project_name = $project->name();
 	my $no = $project->noSqlCnvs("r");
	$filein =  $no->dir."/raw_data.lite";
	#my $dir_out= $project->getVariationsDir("dude");
	#$filein = $dir_out."/".$name.".dude.lid.gz";
	my $no1 = $self->patient()->getTranscriptsDude();
	my $fileout = $no1->filename();
	my $ppn =$self->nproc;
	$ppn = int($self->nproc/2) if $self->nocluster;
	
	
	### Récupération de la version de samtools
	
	my $bin_dev = $self->script_dir;
	
	my $cmd = qq{perl $bin_dev/transcripts/transcripts_cache.pl -patient=$name  -fork=$ppn  -project=$project_name };
	my $type = "transcripts_dude";
	 my $stepname = $self->patient->name."@".$type;
	 $ppn =20;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() &&  $no1->nb_keys >10 ){
		 		$job_bds->skip();
	}
		
	warn Dumper $self->current_sample->jobs();
	$no->close();
	return ($fileout);
}


sub update_annotation{
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	
	my $chr = $self->patient();
	
	my $chr_name = $chr->name();
	
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->nproc;
	if ($self->nocluster){
		$ppn = $self->nproc;
	}
	$ppn = int($ppn/4);
	$ppn = 1 if $ppn < 1;
	
	my $dir = $project->dir_lmdb_score_impact();
	my $fileout = $dir."/".$chr_name.".update";
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/scripts/update_annotations.pl -project=$project_name -chr=$chr_name -fork=$ppn && date > $fileout";
	
	

	 my $type = "update_annotations";
	 my $stepname = $self->patient->name."@".$type;
	 my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
	
}


sub check_store_annotations {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	
	my $chr = $self->patient();
	
	my $chr_name = $chr->name();
	
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->nproc;
	if ($self->nocluster){
		$ppn = $self->nproc;
	}
	$ppn = int($ppn/4);
	$ppn = 1 if $ppn < 1;
#	/data-xfs/polycache//HG19.28/NGS2019_2335/vector//log//check_store_annotations.17.ok
	my $fileout = $project->getCacheBitVectorDir()."/log/check_store_annotations.$chr_name.ok";
	my $log_error = $fileout;
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/scripts/cache_check_step.pl -step=store_annotations  -project=$project_name -chr=$chr_name -fork=$ppn  ";
	
	 my $type = "check_store_annot";
	 my $stepname = $self->patient->name."@".$type;
	 my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($filein);
	
}

sub  strict_denovo{
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	
	my $chr = $self->patient();
	
	my $chr_name = $chr->name();
	
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->nproc;
	if ($self->nocluster){
		$ppn = $self->nproc;
	}
	$ppn = int($ppn/4);
	$ppn = 1 if $ppn < 1;
	
	my $fileout = $project->getCacheBitVectorDir()."/strict-denovo/$chr_name.lite";
	my $fileout2 = $project->getCacheBitVectorDir()."/log/check_strict_denovo.$chr_name.ok";
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/scripts/cache_strict_denovo.pl  -project=$project_name -chr=$chr_name -fork=$ppn";
	
	$cmd .= " && test $fileout && perl $Bin/../polymorphism-cgi/cache_nodb/scripts/cache_check_step.pl -project=$project_name -step=strict_denovo -project=$project_name -chr=$chr_name -fork=$ppn";
	
	
	 my $type = "strict_denovo";
	 my $stepname = $self->patient->name."@".$type;
	 my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout2,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout2);
}

sub  loh {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $chr = $self->patient();
	
	my $chr_name = $chr->name();
	
	my $project = $self->project();
	return ($filein) unless $project->isSomaticStudy();
	my $project_name = $project->name();
	my $ppn = $self->nproc;
	if ($self->nocluster){
		$ppn = $self->nproc;
	}
	$ppn = int($ppn/4);
	$ppn = 1 if $ppn < 1;
	
	my $fileout = $project->getCacheBitVectorDir()."/somatic_loh/$chr_name.lite";
	my $fileout2 = $project->getCacheBitVectorDir()."/log/check_loh.$chr_name.ok";
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/scripts/cache_loh.pl  -project=$project_name -chr=$chr_name -fork=$ppn";
	
	$cmd .= " && test $fileout  && perl $Bin/../polymorphism-cgi/cache_nodb/scripts/cache_check_step.pl -project=$project_name -step=loh -chr=$chr_name -fork=$ppn";
	
	my $type = "loh";
	my $stepname = $self->patient->name."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($filein);
}

sub between {
  my($test ,$fom, $tom)=@_;
  no warnings;
  $fom<$tom ? $test>=$fom && $test<=$tom
            : $test>=$tom && $test<=$fom;
}

sub polydiag {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $project = $self->project();
	my $dir_out = $project->getCacheDir() . "/polydiag_lite/";
	my $p = $self->patient();
	my $project_name = $self->patient()->getProject->name();
	my $type = "polydiag";
	my $root_cmd = "perl $Bin/../polymorphism-cgi//cache_nodb_old/scripts/cache_polydiag.pl";
	
	my $fout;
	my $ppn =1;
	my $cmd = " $root_cmd -fork=$ppn -project=$project_name -patient=".$p->name;
	my $fileout = $dir_out."/".$p->name.".lite";
	my $stepname = $project->name.".".$p->name."@".$type;
	push(@$fout,$fileout);
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
		$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
  		$job_bds->skip();
	}
	return $filein;
}

sub cnv_manue {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $projectName = $self->project->name();
	
	my $ppn = $self->nproc;
	my $type = "cnv-manue";
	my $stepname = $projectName."@".$type;
	my $dir = $self->project->project_log();
	my $methods = $self->project->callingSVMethods();
	return $filein unless @$methods;
	my $fileout   = $self->project->getCNVDir()."/".$self->project->name.".done";

	
	my $cmd = "perl $bin_dev/manue_cnv/SV_global.pl -project=$projectName -fork=$ppn > $fileout";
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	$job_bds->isLogging(1);
	if ($self->unforce() && -e $fileout){
		$job_bds->skip();
	}
	return ($filein);
}

sub identito_vigilence {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $name = $self->patient()->name();
	my $project = $self->patient()->getProject();
	my $project_name = $project->name();
	
	my $projectName = $self->project->name();
	my $ppn = 1;
	$ppn = int($self->nproc/2) if $self->nocluster;
	my $output   = $self->project->quality_dir;#$self->getCacheDir() . "/check_quality";
	my $fileout = "$output/".$self->project->name.".identito_vigilence.txt";
	
	### Récupération de la version de samtools
	
	my $bin_dev = $self->script_dir;
	
	my $cmd = qq{perl $bin_dev/identito_vigilence.pl -project=$projectName > $fileout };
	my $type = "identito";
	my $stepname = $projectName."@".$type;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
		my $nofound;
		foreach my $patient (@{$project->getPatients()}) {
			my $identity_vigilance= $patient->identity_vigilance();
			$nofound ++ unless $identity_vigilance;
		}
		
		 		$job_bds->skip() unless $nofound;
	}
	return ($fileout);
}

sub quality_check {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	
	my $projectName = $self->project->name();
	my $logPlink = $self->project->getProjectPath().'../'.$projectName.'.qualitycheck.resume';
	my $ppn = $self->nproc;
	my $type = "qualtity-check";
	my $stepname = $projectName."@".$type;
	my $dir = $self->project->project_log();
	my $output   = $self->project->quality_dir;#$self->getCacheDir() . "/check_quality";
	my $fileout = "$output/".$self->project->name.".lite";
	
	my $cmd = "perl $bin_dev/quality_check.pl -project=$projectName -fork=$ppn -cache=1>$fileout";
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	$job_bds->isLogging(1);
	if ($self->unforce() && -e $fileout){
		my $no = $self->project->noSqlQuality("r");
		my $data = $no->get($self->project->name,"mendelian");
		my $quality_patients;
		foreach my $line (@{$data->{data}}){
			my $name = $line->{sample}->{text};
			$quality_patients->{$name} ++;
		}
		my $error;
		foreach my $p (@{$self->project->getPatients}){
			warn $p->name unless exists $quality_patients->{$p->name};
			$error =1 unless exists $quality_patients->{$p->name};
		}
		
  		$job_bds->skip() unless $error;
	}
	return ($filein);
}

sub dejavu {
	my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $projectName = $self->project->name();
	my $logPlink = $self->project->getProjectPath().'../'.$projectName.'.dejavu.resume';
	my $ppn = $self->nproc;
	my $type = "dejavu";
	my $stepname = $projectName."@".$type;
	my $dir = $self->project->project_log();
	my $fileout = $dir."/dejavu.log";
	my $cmd = "perl $Bin/../polymorphism-cgi/cache_nodb/scripts/cache_lite_dejavu.pl -project=$projectName >$fileout";
	
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	$job_bds->isLogging(1);
	if ($self->unforce() && -e $fileout){
  		$job_bds->skip();
	}
	return ($filein);
}

sub coverage {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $project = $self->project();
	my $patients = $project->getPatients();
	my $type = "coverage";
	my $fileout = $project->getCacheDir().'/coverage_lite/'.$project->name.'.lite';
	my $stepname = $project->name."@"."$type";
	my $root_cmd = "perl $Bin/../polymorphism-cgi//cache_nodb/scripts/cache_coverage.pl";
	
	my $ppn = "10";
	 $ppn = scalar(@$patients) if  scalar(@$patients) < 40;
	if ($self->nocluster){
		$ppn = $self->nproc;
	}
	my $cmd = " $root_cmd -fork=$ppn -project=".$project->name;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>40,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
  		$job_bds->skip();
	}
	return $filein;
	
}

1;