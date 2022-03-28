package bds_cache_steps;
use strict;
use FindBin qw($Bin);
use lib "$Bin";
#use root_steps;
use Moose;  
use MooseX::Method::Signatures;
use job_bds;
use sample;
use Data::Dumper;
extends (qw(bds_root));
my $bin_dev = qq{$Bin/scripts/scripts_pipeline/};

	
method get_list_cmds_by_chr (Int :$fork, Str :$cmd, Str :$fileout, Str :$log_error) {
	my $project =  $self->project();
	my $project_name = $project->name();
	mkdir $project->getCacheBitVectorDir().'/log' unless (-d $project->getCacheBitVectorDir().'/log');
	my $cmds;
	#my @lChr = (1..22, 'X', 'Y', 'MT');
	foreach my $chr (@{$project->getChromosomes}) {
		my $chr_name = $chr->name();
		my $hcmds;
		my $this_fileout = $fileout;
		$this_fileout =~ s/CHR_NAME/$chr_name/;
		my $this_log_error = $log_error;
		$this_log_error =~ s/CHR_NAME/$chr_name/;
		$hcmds->{fileout} = $this_fileout;
		$hcmds->{cmd} = "$cmd -fork=$fork -project=$project_name -chr=$chr_name 2>$this_log_error";
		$hcmds->{chr} =  $chr_name;
		$hcmds->{ppn} =  $fork;
		push(@$cmds,$hcmds)
	}
	return ($cmds);
}

method store_ids {
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->nproc;
	if ($self->nocluster){
		$ppn = $self->nproc;
	}
	$ppn = int($ppn/4);
	$ppn = 1 if $ppn < 1;
	my $filein = $project->getCacheBitVectorDir().'/global_infos.freeze';
	my $fileout = $project->getCacheBitVectorDir().'/lmdb_cache/CHR_NAME.dv.freeze';
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi//cache_nodb_old/scripts/cache_store_ids.pl";
	my $log_error = $project->getCacheBitVectorDir()."/log/store_ids_CHR_NAME.log";
	my $cmds = $self->get_list_cmds_by_chr(fork => $ppn, cmd => $cmd, fileout => $fileout, log_error => $log_error);
 	my $type = "s_ids";
 	my $fout;
	foreach my $hcmd (@$cmds) {
		my $stepname = $project_name.".".$hcmd->{chr}."@".$type;
		my $this_cmd = $hcmd->{cmd};
		my $fileout = $hcmd->{fileout};
		push(@$fout,$fileout);
		push(@$fout, $fileout);
		my $job_bds = job_bds->new (
			cmd     => [$this_cmd],
			name    => $stepname,
			ppn     => $ppn,
			filein  => [$filein],
			fileout => $fileout,
			type    => $type,
			dir_bds => $self->dir_bds
		);
		$self->current_sample->add_job(job=>$job_bds);
		if ($self->unforce() && -e $fileout){
	  		$job_bds->skip();
		}
		
	}
	my $cmd2 = "/usr/bin/perl $Bin/../polymorphism-cgi//cache_nodb_old/scripts/cache_lite_dejavu.pl  -project=$project_name";
	my $type2 = "dejavu";
	my $stepname = $project_name."@".$type2;
	my $fileout_gv = $project->deja_vu_lite_dir() . "/projects/" . $project_name . ".lite";
	my $job_bds = job_bds->new(cmd=>[$cmd2],name=>$stepname,ppn=>1,filein=>$fout,fileout=>$fileout_gv,type=>$type2,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout_gv){
	  		$job_bds->skip();
		}
	
	return ($fileout_gv);
}

method local_config {
	my $project = $self->project();
	my $project_name = $project->name();
	my $filein = $project->getCacheBitVectorDir().'/global_infos.freeze';
	my $cmd2 = "/usr/bin/perl $Bin/../polymorphism-cgi//cache_nodb_old/scripts/cache_local_config.pl -project=$project_name";
	my $type2 = "local_config";
	my $stepname = $project_name."@".$type2;
	my $fileout = $project->getCacheBitVectorDir().'/genbo_'.$project_name.'.cfg';
	my $job_bds = job_bds->new(cmd=>[$cmd2],name=>$stepname,ppn=>1,filein=>[$filein],fileout=>$fileout,type=>$type2,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout){
  		$job_bds->skip();
	}
	return ($fileout);
}

method dejavu {
	my $project = $self->project();
	my $project_name = $project->name();
	my $filesin;
	foreach my $chr (@{$project->getChromosomes()}) {
		my $filein = $project->getCacheBitVectorDir().'/lmdb_cache/'.$chr->name().'.dv.freeze';
		push(@$filesin, $filein) if (-e $filein);
	}
	my $cmd2 = "/usr/bin/perl $Bin/../polymorphism-cgi//cache_nodb_old/scripts/cache_lite_dejavu.pl  -project=$project_name";
	my $type2 = "dejavu";
	my $stepname = $project_name."@".$type2;
	my $fileout_gv = $project->deja_vu_lite_dir() . "/projects/" . $project_name . ".lite";
	my $job_bds = job_bds->new(cmd=>[$cmd2],name=>$stepname,ppn=>1,filein=>$filesin,fileout=>$fileout_gv,type=>$type2,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	if ($self->unforce() && -e $fileout_gv){
  		$job_bds->skip();
	}
	return ($fileout_gv);
}

method store_annotations {
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->nproc;
	if ($self->nocluster){
		$ppn = $self->nproc;
	}
	$ppn = int($ppn/4);
	$ppn =1 if $ppn <1;
	my $filein = $project->getCacheBitVectorDir().'/lmdb_cache/CHR_NAME.dv.freeze';
	my $fileout = $project->getCacheBitVectorDir().'/lmdb_cache/CHR_NAME/genes_index';
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi//cache_nodb_old/scripts/cache_store_annotations.pl";
	my $log_error = $project->getCacheBitVectorDir()."/log/store_anotations_CHR_NAME.log";
	my $cmds = $self->get_list_cmds_by_chr(fork => $ppn, cmd => $cmd, fileout => $fileout, log_error => $log_error);
	my $fout;	
 	my $type = "s_annot";
	foreach my $hcmd (@$cmds) {
		my $chr_name = $hcmd->{chr};
		my $this_filein = $filein;
		$this_filein =~ s/CHR_NAME/$chr_name/;
		my $stepname = $project_name.".".$hcmd->{chr}."@".$type;
		my $this_cmd = $hcmd->{cmd};
		my $fileout = $hcmd->{fileout};
		push(@$fout, $fileout);
		my $job_bds = job_bds->new (
			cmd     => [$this_cmd],
			name    => $stepname,
			ppn     => $ppn,
			filein  => [$this_filein],
			fileout => $fileout,
			type    => $type,
			dir_bds => $self->dir_bds
		);
		$self->current_sample->add_job(job=>$job_bds);
		if ($self->unforce() && -e $fileout){
	  		$job_bds->skip();
		}
	}
	return ($filein);
}

method strict_denovo {
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->nproc;
	if ($self->nocluster){
		$ppn = $self->nproc;
	}
	$ppn = int($ppn/4);
	my $filein = $project->getCacheBitVectorDir().'/lmdb_cache/CHR_NAME/genes_index';
	my $fileout = $project->getCacheBitVectorDir().'/strict-denovo/CHR_NAME.lite';
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi//cache_nodb_old/scripts/cache_strict_denovo.pl";
	my $log_error = $project->getCacheBitVectorDir()."/log/strict_denovo_CHR_NAME.log";
	my $cmds = $self->get_list_cmds_by_chr(fork => $ppn, cmd => $cmd, fileout => $fileout, log_error => $log_error);
	my $fout;	
 	my $type = "s_denovo";
	foreach my $hcmd (@$cmds) {
		my $chr_name = $hcmd->{chr};
		my $this_filein = $filein;
		$this_filein =~ s/CHR_NAME/$chr_name/;
		my $stepname = $project_name.".".$hcmd->{chr}."@".$type;
		my $this_cmd = $hcmd->{cmd};
		my $fileout = $hcmd->{fileout};
		push(@$fout, $fileout);
		my $job_bds = job_bds->new (
			cmd     => [$this_cmd],
			name    => $stepname,
			ppn     => $ppn,
			filein  => [$this_filein],
			fileout => $fileout,
			type    => $type,
			dir_bds => $self->dir_bds
		);
		$self->current_sample->add_job(job=>$job_bds);
		if ($self->unforce() && -e $fileout){
	  		$job_bds->skip();
		}
	}
	return ($filein);
}

method loh {
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->nproc;
	if ($self->nocluster){
		$ppn = $self->nproc;
	}
	$ppn = int($ppn/4);
	my $filein = $project->getCacheBitVectorDir().'/lmdb_cache/CHR_NAME/genes_index';
	my $fileout = $project->getCacheBitVectorDir().'/somatic_loh/CHR_NAME.lite';
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi//cache_nodb_old/scripts/cache_loh.pl";
	my $log_error = $project->getCacheBitVectorDir()."/log/loh_CHR_NAME.log";
	unless ($project->isSomaticStudy()) {
		return ($filein);
	}
	my $cmds = $self->get_list_cmds_by_chr(fork => $ppn, cmd => $cmd, fileout => $fileout, log_error => $log_error);
	my $fout;	
 	my $type = "loh";
	foreach my $hcmd (@$cmds) {
		my $chr_name = $hcmd->{chr};
		my $this_filein = $filein;
		$this_filein =~ s/CHR_NAME/$chr_name/;
		my $stepname = $project_name.".".$hcmd->{chr}."@".$type;
		my $this_cmd = $hcmd->{cmd};
		my $fileout = $hcmd->{fileout};
		push(@$fout, $fileout);
		my $job_bds = job_bds->new (
			cmd     => [$this_cmd],
			name    => $stepname,
			ppn     => $ppn,
			filein  => [$this_filein],
			fileout => $fileout,
			type    => $type,
			dir_bds => $self->dir_bds
		);
		$self->current_sample->add_job(job=>$job_bds);
		if ($self->unforce() && -e $fileout){
	  		$job_bds->skip();
		}
	}
	return ($filein);
}

method global_infos {
	my $project = $self->project();
	my $project_name = $project->name();
	my $filein = '';
	my $fileout = $project->getCacheBitVectorDir().'/global_infos.freeze';
	my $log_error = $project->getCacheBitVectorDir()."/log/global_infos.log";
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi//cache_nodb_old/scripts/cache_global_infos.pl -project=$project_name 2>$log_error";
 	my $type = "global_infos";
	my $stepname = $project_name."@".$type;
	my $job_bds = job_bds->new (
		cmd     => [$cmd],
		name    => $stepname,
		ppn     => 1,
		filein  => [$filein],
		fileout => $fileout,
		type    => $type,
		dir_bds => $self->dir_bds
	);
	$self->current_sample->add_job(job=>$job_bds);
	return ($filein);
}



###### CHECKS METHODS ######



method check_store_ids {
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->nproc;
	if ($self->nocluster){
		$ppn = $self->nproc;
	}
	$ppn = int($ppn/4);
	my $filein = $project->getCacheBitVectorDir().'/lmdb_cache/CHR_NAME.dv.freeze';
	my $fileout = $project->getCacheBitVectorDir()."/log/check_store_ids_CHR_NAME.log";
	my $log_error = $fileout;
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi//cache_nodb_old/scripts/cache_check_query_vcf.pl";
	my $cmds = $self->get_list_cmds_by_chr(fork => $ppn, cmd => $cmd, fileout => $fileout, log_error => $log_error);
	my $fout;	
 	my $type = "check_ids";
	foreach my $hcmd (@$cmds) {
		my $stepname = $project_name.".".$hcmd->{chr}."@".$type;
		my $chr_name = $hcmd->{chr};
		$filein =~ s/CHR_NAME/$chr_name/;
		my $this_cmd = $hcmd->{cmd};
		my $fileout = $hcmd->{fileout};
		push(@$fout, $fileout);
		my $job_bds = job_bds->new (
			cmd     => [$this_cmd],
			name    => $stepname,
			ppn     => $ppn,
			filein  => [$filein],
			fileout => $fileout,
			type    => $type,
			dir_bds => $self->dir_bds
		);
		$self->current_sample->add_job(job=>$job_bds);
	}
	return ($fileout);
}

method check_store_annotations {
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->nproc;
	if ($self->nocluster){
		$ppn = $self->nproc;
	}
	$ppn = int($ppn/4);
	my $filein = $project->getCacheBitVectorDir().'/lmdb_cache/CHR_NAME/genes_index';
	my $fileout = $project->getCacheBitVectorDir()."/log/check_store_anotations_CHR_NAME.log";
	my $log_error = $fileout;
	my $cmd = "/usr/bin/perl $Bin/../polymorphism-cgi//cache_nodb_old/scripts/cache_check_step.pl -step=store_annotations";
	my $cmds = $self->get_list_cmds_by_chr(fork => $ppn, cmd => $cmd, fileout => $fileout, log_error => $log_error);
	my $fout;	
 	my $type = "check_annot";
	foreach my $hcmd (@$cmds) {
		my $chr_name = $hcmd->{chr};
		my $this_filein = $filein;
		$this_filein =~ s/CHR_NAME/$chr_name/;
		my $stepname = $project_name.".".$hcmd->{chr}."@".$type;
		my $this_cmd = $hcmd->{cmd};
		my $fileout = $hcmd->{fileout};
		push(@$fout, $fileout);
		my $job_bds = job_bds->new (
			cmd     => [$this_cmd],
			name    => $stepname,
			ppn     => $ppn,
			filein  => [$this_filein],
			fileout => $fileout,
			type    => $type,
			dir_bds => $self->dir_bds
		);
		$self->current_sample->add_job(job=>$job_bds);
	}
	return ($filein);
}

method check_strict_denovo {
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->nproc;
	if ($self->nocluster){
		$ppn = $self->nproc;
	}
	$ppn = int($ppn/4);
	my $filein = $project->getCacheBitVectorDir().'/strict-denovo/CHR_NAME.lite';
	my $fileout = $project->getCacheBitVectorDir()."/log/check_strict_denovo_CHR_NAME.log";
	my $log_error = $fileout;
	my $cmd = "perl $Bin/../polymorphism-cgi//cache_nodb_old/scripts/cache_check_step.pl -step=strict_denovo";
	my $cmds = $self->get_list_cmds_by_chr(fork => $ppn, cmd => $cmd, fileout => $fileout, log_error => $log_error);
	my $fout;	
 	my $type = "check_s_denovo";
	foreach my $hcmd (@$cmds) {
		my $chr_name = $hcmd->{chr};
		my $this_filein = $filein;
		$this_filein =~ s/CHR_NAME/$chr_name/;
		my $stepname = $project_name.".".$hcmd->{chr}."@".$type;
		my $this_cmd = $hcmd->{cmd};
		my $fileout = $hcmd->{fileout};
		push(@$fout, $fileout);
		my $job_bds = job_bds->new (
			cmd     => [$this_cmd],
			name    => $stepname,
			ppn     => $ppn,
			filein  => [$this_filein],
			fileout => $fileout,
			type    => $type,
			dir_bds => $self->dir_bds
		);
		$self->current_sample->add_job(job=>$job_bds);
	}
	return ($filein);
}

method check_loh {
	my $project = $self->project();
	my $project_name = $project->name();
	my $ppn = $self->nproc;
	if ($self->nocluster){
		$ppn = $self->nproc;
	}
	$ppn = int($ppn/4);
	my $filein = $project->getCacheBitVectorDir().'/somatic_loh/CHR_NAME.lite';
	my $fileout = $project->getCacheBitVectorDir()."/log/check_loh_CHR_NAME.log";
	my $log_error = $fileout;
	my $cmd = "perl $Bin/../polymorphism-cgi//cache_nodb_old/scripts/cache_check_step.pl -step=loh";
	unless ($project->isSomaticStudy()) {
		return ($filein);
	}
	my $cmds = $self->get_list_cmds_by_chr(fork => $ppn, cmd => $cmd, fileout => $fileout, log_error => $log_error);
	my $fout;	
 	my $type = "check_loh";
	foreach my $hcmd (@$cmds) {
		my $chr_name = $hcmd->{chr};
		my $this_filein = $filein;
		$this_filein =~ s/CHR_NAME/$chr_name/;
		my $stepname = $project_name.".".$hcmd->{chr}."@".$type;
		my $this_cmd = $hcmd->{cmd};
		my $fileout = $hcmd->{fileout};
		push(@$fout, $fileout);
		my $job_bds = job_bds->new (
			cmd     => [$this_cmd],
			name    => $stepname,
			ppn     => $ppn,
			filein  => [$this_filein],
			fileout => $fileout,
			type    => $type,
			dir_bds => $self->dir_bds
		);
		$self->current_sample->add_job(job=>$job_bds);
	}
	return ($filein);
}




method polydiag (Str :$filein){
	my $project = $self->project();
	my $dir_out = $project->getCacheDir() . "/polydiag_lite/";
	my $patients =  $project->getPatients();
	my $project_name = $self->patient()->getProject->name();
	 my $type = "cacheDiag";
	my $root_cmd = "perl $Bin/../polymorphism-cgi//cache_nodb_old/scripts/cache_polydiag.pl";
	my $fout;
	my $ppn =1;
	foreach my $p (@$patients) {
		my $cmd = " $root_cmd -fork=$ppn -project=$project_name -patient=".$p->name;
		my $fileout = $dir_out."/".$p->name.".lite";
		my $stepname = $project->name.".".$p->name."@".$type;
		push(@$fout,$fileout);
		my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
			$self->current_sample->add_job(job=>$job_bds);
		if ($self->unforce() && -e $fileout){
	  		$job_bds->skip();
		}
		
		
	}
	return $filein;
	
}


method quality_check (Str :$filein){
	my $projectName = $self->project->name();
	my $logPlink = $self->project->getProjectPath().'../'.$projectName.'.plink.resume';
	my $ppn = $self->nproc;
#	unless ($filein){
#		$filein = $self->project->getVariationsDir($self->method_calling)."/".$projectName.".vcf.gz";
#		die() unless -e $filein;
#	}
	my $type = "qualtity-check";
	my $stepname = $projectName."@".$type;
	my $dir = $self->project->project_log();
	my $fileout = $dir."/quality_check.log";
#	my $cmd = "perl $bin_dev/quality_check.pl -project=$projectName  -vcf_file=$filein -fork=$ppn -cache=1>$fileout";
	my $cmd = "perl $bin_dev/quality_check.pl -project=$projectName -fork=$ppn -cache=1>$fileout";
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job(job=>$job_bds);
	$job_bds->isLogging(1);
		#@@@@@@@@@@@@@@@@@@@@@
		#todo check for all individual file 
	if ($self->unforce() && -e $fileout){
	  		$job_bds->skip();
	}
	return ($filein);
	
}



method coverage (Str :$filein){
	my $project = $self->project();
	my $patients = $project->getPatients();
	my $type = "coverage";
	my $fileout = $project->getCacheDir().'/coverage_lite/'.$project->name.'.lite';
	my $stepname = $project->name."@"."$type";
	my $root_cmd = "perl $Bin/../polymorphism-cgi//cache_nodb_old/scripts/cache_coverage.pl";
	my $ppn = "40";
	 $ppn = scalar(@$patients) if  scalar(@$patients) < 40;
	if ($self->nocluster){
		$ppn = $self->nproc;
	}
	#$ppn = int($ppn/4);
	
	#warn $ppn;
	#die();

		my $cmd = " $root_cmd -fork=$ppn -project=".$project->name;
		my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
		$self->current_sample->add_job(job=>$job_bds);
		if ($self->unforce() && -e $fileout){
	  		$job_bds->skip();
		}
	return $filein;
	
}







#
#method check_global_infos (Str :$filein!) {
#	
#}
#
#method check_coverage (Str :$filein!) {
#	
#}

1;