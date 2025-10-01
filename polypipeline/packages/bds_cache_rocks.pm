package bds_cache_rocks;
use strict;
use FindBin qw($Bin);
use lib "$Bin";
#use root_steps;
use Moose;  
use job_bds;
use sample;
use Data::Dumper;
extends (qw(bds_cache_steps));
my $bin_dev = qq{$Bin/scripts/scripts_pipeline/rocks/};
my $bin_cache = qq{$Bin/../polymorphism-cgi/cache_nodb/scripts/rocks/};
	
sub isrocksOK{
	my ($dirs) = @_;
	
	foreach my $dir (@$dirs){
		return undef unless -e ($dir."/CURRENT/");
	}
	return 1;
} 



sub cache_store_duck {
	my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $projectName = $self->project->name();
	my $ppn = $self->nproc * 2 ;
	my $type = "duck_cache_store";
	my $stepname = $projectName."@".$type;
#	my $fileout = $self->project->project_log()."/dejavu_parquet.log";
	my $dir_parquet = $self->project->buffer->dejavu_parquet_dir();
	my $fileout = $self->project->parquet_cache_variants; ;
	if (not $self->project->infosProject->{dejavu}) { $fileout .= '.no_dejavu'; }
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/scripts/rocks/duck_cache_store_annotations.pl  -project=$projectName ";
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	$job_bds->isLogging(1);
	if ($self->unforce() && -e $fileout){
  		$job_bds->skip();
	}
	return ($filein);
}

sub tiny_rocks {
	my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $projectName = $self->project->name();
	my $ppn = $self->nproc * 2 ;
	my $type = "tiny_rocks";
	my $stepname = $projectName."@".$type;
#	my $fileout = $self->project->project_log()."/dejavu_parquet.log";
	# my $dir_parquet = $self->project->rocks_directory_beegfs()."/tiny_rocks/".$self->project->name.".store";;
	my $fileout = $self->project->rocks_directory_beegfs()."/tiny_rocks/".$self->project->name.".store";
	if (not $self->project->infosProject->{dejavu}) { $fileout .= '.no_dejavu'; }
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/scripts/rocks/tiny_rocks.pl  -project=$projectName ";
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	$job_bds->isLogging(1);
	if ($self->unforce() && -e $fileout){
  		$job_bds->skip();
	}
	return ($filein);
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
	my $fileout = $project->rocks_directory("logs")."/cache_store_ids.".$chr_name.".ok";
	my $cmd = "/usr/bin/perl $bin_cache/cache_store_ids.pl -project=$project_name -chr=$chr_name -fork=$ppn -file=$fileout &&  test -e  $fileout";
	 my $type = "cache_store_ids";
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
	my $fileout =  $project->rocks_directory("logs")."/cache_store_annotations.".$chr_name.".ok";
	my $cmd = "/usr/bin/perl $bin_cache/cache_store_annotations.pl -project=$project_name -chr=$chr_name -fork=$ppn -file=$fileout &&  test -e  $fileout";
	
	 my $type = "store_annotations";
	 my $stepname = $self->patient->name."@".$type;
	 my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
		 		$job_bds->skip();
	}
	return ($fileout);
	
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
	
	my $fileout =  $project->rocks_directory("logs")."/cache_strict_denovo.".$chr_name.".ok";
	my $cmd = "/usr/bin/perl $bin_cache/cache_strict_denovo.pl  -project=$project_name -chr=$chr_name -fork=$ppn -file=$fileout &&  test -e  $fileout";

	
	
	 my $type = "strict_denovo";
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
 
	my $fileout =  $project->rocks_directory("logs")."/update_hash_variant_chromosome.".$chr_name.".ok";
	my $cmd = "/usr/bin/perl $bin_cache/update_hash_variant_chromosome.pl -project=$project_name -chr=$chr_name -fork=$ppn -file=$fileout &&  test -e  $fileout";
	
	
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

sub merge_objects {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $ppn = 5;
	my $project = $self->project;
	my $project_name = $self->project->name();
	my $type = "merge-objs";
	
	my $stepname = $project_name."@".$type;
	my $dir = $self->project->project_log();
	my $output   = $self->project->quality_dir;#$self->getCacheDir() . "/check_quality";
	my $fileout =  $project->rocks_directory("logs")."/merge_objects".".ok";
	my $cmd = "/usr/bin/perl $bin_cache/merge_polyviewer_db.pl -project=$project_name  -file=$fileout &&  test -e  $fileout";
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	$job_bds->isLogging(1);
	if ($self->unforce() && -e $fileout){
#		my $no = $self->project->noSqlQuality("r");
#		my $data = $no->get($self->project->name,"mendelian");
#		my $quality_patients;
#		foreach my $line (@{$data->{data}}){
#			my $name = $line->{sample}->{text};
#			$quality_patients->{$name} ++;
#		}
#		my $error;
#		foreach my $p (@{$self->project->getPatients}){
#			warn $p->name unless exists $quality_patients->{$p->name};
#			$error =1 unless exists $quality_patients->{$p->name};
#		}
#  		$job_bds->skip() unless $error;
		$job_bds->skip();
	}
	return ($filein);
}
sub merge_patients {
		my ($self,$hash) = @_;
	my $filein = $hash->{filein};
	my $ppn = 10;
	my $project = $self->project;
	my $project_name = $self->project->name();
	my $type = "merge-patients";
	
	my $stepname = $project_name."@".$type;
	my $dir = $self->project->project_log();
	my $output   = $self->project->quality_dir;#$self->getCacheDir() . "/check_quality";
	my $fileout =  $project->rocks_directory("logs")."/merge_objects".".ok";

	my $cmd = "/usr/bin/perl $bin_cache/global_annotation_by_patient.pl -project=$project_name -file=$fileout -fork=$ppn &&  test -e  $fileout";
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	$job_bds->isLogging(1);
	if ($self->unforce() && -e $fileout){
#		my $no = $self->project->noSqlQuality("r");
#		my $data = $no->get($self->project->name,"mendelian");
#		my $quality_patients;
#		foreach my $line (@{$data->{data}}){
#			my $name = $line->{sample}->{text};
#			$quality_patients->{$name} ++;
#		}
#		my $error;
#		foreach my $p (@{$self->project->getPatients}){
#			warn $p->name unless exists $quality_patients->{$p->name};
#			$error =1 unless exists $quality_patients->{$p->name};
#		}
#  		$job_bds->skip() unless $error;
		$job_bds->skip();
	}
	return ($filein);
}
1;