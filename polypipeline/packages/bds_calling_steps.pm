package bds_calling_steps;
use strict;
use FindBin qw($Bin);
use lib "$Bin";
#use root_steps;
use Moose;  
use MooseX::Method::Signatures;
use job_bds;
use sample;
use Data::Dumper;
use Text::CSV qw( csv );
extends (qw(bds_root));
my $bin_dev = qq{$Bin/scripts/scripts_pipeline/};





has 'method_calling' =>(
is =>'rw',
lazy=>1,
		default => sub {
			return "haplotypecaller";
		}
);






method correct_vcf(Str :$filein!){
	my $project = $self->project;
	my $name = $project->name();
	my $dir_out= $project->getCallingPipelineDir($self->method_calling);
	my $snp_out = $dir_out . "/" . $project->name . ".snps.vcf";
	my $indels_out = $dir_out . "/" . $project->name . ".indels.vcf";
	
	my $type = "correct";	
	
	#my $fileout = $filein;
	unless ($filein){
		$filein = $dir_out."/$name.uni.vcf";
		die("can't find $filein ") unless -e $filein;
	}
	my $fork = 1;
	my $patient =  $self->argument_patient(); 
	my $fileout = "$dir_out/$name.final.vcf";
	#@@@@@@@@@@@@@@@@@@@@@@@
	#todo 
	#check if patient is in filein 
	#@@@@@@@@@@@@@@@@@@@@@@
	my $stepname = $name."@".$type;
	 my $cmd = "perl $bin_dev/correct_gatk.pl -vcf=$filein -project=$name -dir=$dir_out -snp_out=$snp_out -indels_out=$indels_out -patient=$patient -method=haplotypecaller";
	  my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$fork,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
		$self->current_sample->add_job({job=>$job_bds});
	 if ($self->unforce() && -e $fileout){
	 	$job_bds->skip();
	#  	$self->add_skip_steps($stepname);
	#	return ($previous,$snp_out);
	}
	  return ($fileout);	
}

method filter_vcf4(Str :$filein!){
	my $project = $self->project;
	my $name = $project->name();
	my $dir_out= $project->getCallingPipelineDir($self->method_calling);

	
	my $type = "correct";	
	
	#my $fileout = $filein;
	unless ($filein){
		$filein = $dir_out."/$name.final.vcf";
		die("can't find $filein ") unless -e $filein;
	}
	my $fork = 1;
	my $patient =  $self->argument_patient(); 
	my $fileout = "$dir_out/$name.final.vcf";
	#@@@@@@@@@@@@@@@@@@@@@@@
	#todo 
	#check if patient is in filein 
	#@@@@@@@@@@@@@@@@@@@@@@
	my $stepname = $name."@".$type;
	 my $cmd = "perl $bin_dev/filter_haplotypecaller.pl -vcf=$filein -project=$name -dir=$dir_out  -patient=$patient -method=haplotypecaller4";

	  my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$fork,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
		$self->current_sample->add_job({job=>$job_bds});
	 if ($self->unforce() && -e $fileout){
	 	$job_bds->skip();
	#  	$self->add_skip_steps($stepname);
	#	return ($previous,$snp_out);
	}
	  return ($fileout);	
}


method correct_vcf4(Str :$filein!){
	my $project = $self->project;
	my $name = $project->name();
	my $dir_out= $project->getCallingPipelineDir($self->method_calling);
	my $snp_out = $dir_out . "/" . $project->name . ".snps.vcf";
	my $indels_out = $dir_out . "/" . $project->name . ".indels.vcf";
	
	my $type = "correct";	
	
	#my $fileout = $filein;
	unless ($filein){
		$filein = $dir_out."/$name.uni.vcf";
		die("can't find $filein ") unless -e $filein;
	}
	my $fork = 1;
	my $patient =  $self->argument_patient(); 
	my $fileout = "$dir_out/$name.final.vcf";
	#@@@@@@@@@@@@@@@@@@@@@@@
	#todo 
	#check if patient is in filein 
	#@@@@@@@@@@@@@@@@@@@@@@
	my $stepname = $name."@".$type;
	 my $cmd = "perl $bin_dev/correct_gatk.pl -vcf=$filein -project=$name -dir=$dir_out -snp_out=$snp_out -indels_out=$indels_out -patient=$patient -method=haplotypecaller4";
	 warn $cmd;
	  my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$fork,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
		$self->current_sample->add_job({job=>$job_bds});
	 if ($self->unforce() && -e $fileout){
	 	$job_bds->skip();
	#  	$self->add_skip_steps($stepname);
	#	return ($previous,$snp_out);
	}
	  return ($fileout);	
}

method genotype_gvcf4 (Str :$filein! ){
	my $project = $self->project;
	my $name = $project->name();
	my $arg = $self->argument_patient(); 
	my $dir_out_vcf= $project->getCallingPipelineDir("haplotypecaller4");
	my $fileout = $dir_out_vcf."/".$name.".hc.vcf";
	my $type = "genotype";
	my $stepname =$name."@".$type;	
	my $patients = $project->get_list_patients($arg);
	foreach my $patient (@$patients){
			my $gvcf = $patient->getGvcfFile("haplotypecaller4");#$dir_out."/".$patient->name.".g.vcf.gz";
			die("$gvcf you don't have gvcf for at least this patient restart after_lifescope on this project ") unless $gvcf ;#-e $gvcf;
	}
	my $fork =$self->nproc;
	my $size_window = 5_000_000;
	my $cmd = "perl $bin_dev/gatk-4/genotype_gvcf.pl -project=".$name." -vcf=$fileout -window=$size_window -patient=$arg -fork=$fork ";;
	warn $cmd;
	 my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>40,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
			$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
	  	$job_bds->skip();
	}
	return ($fileout);
}


method genotype_gvcf (Str :$filein! ){
	my $project = $self->project;
	my $name = $project->name();
	my $arg = $self->argument_patient(); 
	my $dir_out_vcf= $project->getCallingPipelineDir("haplotypecaller");
	my $fileout = $dir_out_vcf."/".$name.".hc.vcf";
	my $type = "genotype";
	my $stepname =$name."@".$type;
	
	
	my $patients = $project->get_list_patients($arg);
	foreach my $patient (@$patients){
			my $gvcf = $patient->getGvcfFile();#$dir_out."/".$patient->name.".g.vcf.gz";
			my $patient_name = $patient->name();
			die("\nYou don't have gvcf for at least this patient restart after_lifescope on this project. Die.\nPatient Name: $patient_name\nExpected GVCF: $gvcf\n\nERROR: die.\n\n") unless $gvcf ;#-e $gvcf;
	}
	my $fork =$self->nproc;
	my $size_window = 1_000_000;
	my $cmd = "perl $bin_dev/genotype_gvcf.pl -project=".$name." -vcf=$fileout -window=$size_window -patient=$arg -fork=$fork ";;
	 my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>40,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
			$self->current_sample->add_job({job=>$job_bds});
	if ($self->unforce() && -e $fileout){
	  	$job_bds->skip();
	}
	return ($fileout);
}
method move_vcf (Str :$filein! ){
	my $name = $self->project->name();
	my $dir_out= $self->project->getCallingPipelineDir($self->method_calling);
	my $arg = $self->argument_patient(); 
	my $type = "move-vcf";
	my $stepname = $name."@".$type;
	my $fileout = $self->project->getVariationsDir($self->method_calling)."/".$name.".vcf.gz";
	
	my $cmd = "perl $bin_dev/move_vcf.pl  -project=$name -vcf_dir=$dir_out -patient=$arg -method_calling=haplotypecaller";
	my $ppn = 1;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
		$self->current_sample->add_job({job=>$job_bds});
	#@@@@@@@@@@@@@@@@@@@@@
	#todo check for all individual file 
	if ($self->unforce() && -e $fileout){
	  	$job_bds->skip();
	}
	return ($fileout);
}

method move_and_split_vcf4 (Str :$filein! ){
	my $name = $self->project->name();
	my $dir_out= $self->project->getCallingPipelineDir($self->method_calling);
	unless ($filein){
		$filein = $dir_out."/".$self->project->name().".final.vcf";
		die($filein) unless -e $filein;
	}
	warn $filein;
	my $arg = $self->argument_patient(); 
	my $type = "move_split-vcf";
	my $stepname = $name."@".$type;
	my $fileout = $self->project->getVariationsDir($self->method_calling)."/".$name.".vcf.gz";
	warn $fileout;
	my $cmd = "perl $bin_dev/split_haplotypecaller.pl  -project=$name -vcf=$filein -patient=$arg -method_calling=haplotypecaller4";
	my $ppn = 1;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
		$self->current_sample->add_job({job=>$job_bds});
	#@@@@@@@@@@@@@@@@@@@@@
	#todo check for all individual file 
	if ($self->unforce() && -e $fileout){
	  	$job_bds->skip();
	}
	return ($fileout);
}

method move_vcf4 (Str :$filein! ){
	my $name = $self->project->name();
	my $dir_out= $self->project->getCallingPipelineDir($self->method_calling);
	my $arg = $self->argument_patient(); 
	my $type = "move-vcf";
	my $stepname = $name."@".$type;
	my $fileout = $self->project->getVariationsDir($self->method_calling)."/".$name.".vcf.gz";
	warn $fileout;
	my $cmd = "perl $bin_dev/move_vcf.pl  -project=$name -vcf_dir=$dir_out -patient=$arg -method_calling=haplotypecaller4";
	my $ppn = 1;
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
		$self->current_sample->add_job({job=>$job_bds});
	#@@@@@@@@@@@@@@@@@@@@@
	#todo check for all individual file 
	if ($self->unforce() && -e $fileout){
	  	$job_bds->skip();
	}
	return ($fileout);
}
 

method calling_larges_indels  (Str :$filein! ){
	my $project = $self->project;
	my $project_name = $project->name();
	my @jobs;
	my @files;
	
	my $dir_out= $project->getCallingPipelineDir("gvcf");
	my $arg = $self->argument_patient(); 
	my $patients = $project->get_list_patients($arg);
	my $fork =16;
	my $nb =1;
	foreach my $patient (@$patients){
		my $fileout =  $patient->getCnvsFile();
		my $filein =  $patient->getBamFile();
		next   if ($self->unforce() && -e $fileout);
		my $name = $patient->name;
	
	
		my $cmd = qq{perl $bin_dev/calculate_large_indels.pl -patient=$name -fork=16 -project=$project_name };
		
	
		my $ppn = 16;
		my $type = "lifinder";
		my $stepname = $name."@".$type."#$nb";
			$nb++;
		my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
		$self->current_sample->add_job({job=>$job_bds});
		#@@@@@@@@@@@@@@@@@@@@@
		#todo check for all individual file 
		if ($self->unforce() && -e $fileout){
	  		$job_bds->skip();
		}
		
		
	}
	
	return ("");
}


method count_featureCounts  (Str :$filein! ){
	my $project = $self->project;
	my $project_name = $project->name();
	
	my $dir_out= $project->getCountingDir("featureCounts");
	my $fileout = $dir_out."/".$project_name.".count.genes.txt";
	my $fileout2 = $dir_out."/".$project_name.".count.exons.txt";
	my $arg = $self->argument_patient(); 
	my $patients = $project->get_list_patients($arg);
	my @bams;
	my @sed_cmd;
	my @sed_cmd2;
	my $align_method;
	my $profile;
	my %strands;
	foreach my $patient (@$patients){
		my $run = $patient->getRun();
		my $type = $run->infosRun->{method};
		my $bam = $patient->getBamFile;
		my $name = $patient->name;
		my $metrics = $project->getCountingDir("featureCounts") . "/metrics/$name.metrics";
		die("Can't find '$metrics'") unless (-e $metrics);
		push(@bams,$bam);
		$bam =~ s/\//\\\//g;
		push(@sed_cmd,qq{sed -i "2s/$bam/$name/" $fileout} );
		push(@sed_cmd2,qq{sed -i "2s/$bam/$name/" $fileout2} );
		
		my $aoa = csv (in => $metrics, sep => "\t");
		my $pct_r1 = $aoa->[7]->[13] if ($aoa->[6]->[13] eq 'PCT_R1_TRANSCRIPT_STRAND_READS');
		my $pct_r2 = $aoa->[7]->[14] if ($aoa->[6]->[14] eq 'PCT_R2_TRANSCRIPT_STRAND_READS');
		die("ERROR parsing '$metrics': no 'PCT_R1_TRANSCRIPT_STRAND_READS' found: ".$aoa->[6]->[13].' -> '.$aoa->[7]->[13]) unless (defined $pct_r1);
		die("ERROR parsing '$metrics': no 'PCT_R1_TRANSCRIPT_STRAND_READS' found: ".$aoa->[6]->[14].' -> '.$aoa->[7]->[14]) unless (defined $pct_r2);
		die("ERROR pct R1 and R2 transcript strand reads are both zero / anormal for '$name': R1 = $pct_r1\tR2 = $pct_r2\n$metrics") if ($pct_r1 + $pct_r2 != 1);
		warn "$name\tR1 = $pct_r1";
		$strands{'-s 1 '} ++ if ($pct_r1 >= 0.9);
		$strands{'-s 2 '} ++ if ($pct_r1 <= 0.1);
		$strands{'-s 0 '} ++ if ($pct_r1 >= 0.4 and $pct_r1 <= 0.6);
		$strands{'error'}->{"$name"} = $pct_r1 if (($pct_r1 > 0.1 and $pct_r1 < 0.4) or ($pct_r1 > 0.6 and $pct_r1 < 0.9));
		$align_method = $patient->alignmentMethod();
		$profile = $patient->getSampleProfile();
	}
	warn 'Strands'.Dumper \%strands;
	my @strands = keys %strands;
	die("Error: pct R1 transcript strand reads:\n".Dumper \%strands) if (grep{/error/} @strands);
	die("More than one strand for the ".scalar @$patients." patients in project $project_name:\n".Dumper \%strands) unless (scalar @strands);

	my $ppn = 16;
	my $gtf = $project->gtf_file();
	#$gtf = $project->gtf_file_dragen() if $align_method eq "star" || $align_method eq "dragen-align";
	
	my $featureCounts = $project->buffer->software("featureCounts");
	my $strand = "-s 1 ";
	$strand = "-s 2 " if $profile eq "bulk illumina pcr-free" or $profile eq "bulk ribozero pcr-free" or $profile eq "bulk NEB-directional pcr-free" or $profile eq "bulk watchmaker pcr-free";
	$strand = "-s 0 " if $profile eq "bulk neb pcr-free" ;
	die("Strands from metrics (".$strands[0].") and from profile ($strand) don't match") unless ($strand eq $strands[0]);
	my $strand = $strands[0];
	
	my $sed = join(" && ",@sed_cmd);
	my $cmd = "$featureCounts -T $ppn   -a $gtf --ignoreDup -o $fileout -p -t exon  $strand ".join(" ",@bams)." && $sed";
	my $type = "featureCounts-genes";
	my $stepname = $project_name."@".$type;
	my $job_bds = job_bds->new(
		cmd=>[$cmd],
		name=>$stepname,
		ppn=>$ppn,
		filein=>[@bams],
		fileout=>$fileout,
		type=>$type,
		dir_bds=>$self->dir_bds
		);
	$self->current_sample->add_job( {job=>$job_bds} );
	if ($self->unforce() && -e $fileout){
  		$job_bds->skip();
	}
		
	my $sed2 = join(" && ",@sed_cmd2);
	my $cmd2 = "$featureCounts -T $ppn -a $gtf -f   -t exon  -O --ignoreDup -p  -o $fileout2 $strand ".join(" ",@bams)." && $sed2";
	my $type2 = "featureCounts-exons";
	my $stepname2 = $project_name."@".$type2;
	my $job_bds2 = job_bds->new(
		cmd=>[$cmd2],
		name=>$stepname2,
		ppn=>$ppn,
		filein=>[@bams],
		fileout=>$fileout2,
		type=>$type2,
		dir_bds=>$self->dir_bds
		);
	$self->current_sample->add_job( {job=>$job_bds2} );
	if ($self->unforce() && -e $fileout2){
  		$job_bds2->skip();
	}
	return ($filein);
}


method xhmm (Str :$filein){
	
	my $projectName = $self->project->name();
	my $logPlink = $self->project->getProjectPath().'../'.$projectName.'.plink.resume';
	my $ppn = 16;
	my $tfilein;
	my $coverage_dir = $self->project->getRootDir() . "/align/coverage/depth/";
	foreach my $patient (@{$self->project->getPatients}){
		my $file = $coverage_dir."/".$patient->name.".depth";
		die("file depth not found : $file") unless -e $file;
		push(@$tfilein,$file);
	}
	unless ($filein){
		$filein = $self->project->getVariationsDir($self->method_calling)."/".$projectName.".vcf.gz";
		#die() unless -e $filein;
	}
	
	my $type = "xhmm";
	my $stepname = $projectName."@".$type;
	my $dir = $self->project->project_log();
	my $final_dir = $self->project->getVariationsDir("xhmm");
	my $fileout = $final_dir."/".$self->project->name.".xcnv";
	
	my $cmd = "perl $bin_dev/xhmm.pl -project=$projectName  ";
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>\@$tfilein,fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
		#@@@@@@@@@@@@@@@@@@@@@
		#todo check for all individual file 
	if ($self->unforce() && -e $fileout){
	  		$job_bds->skip();
	}
	return ($filein);
	
}

method dude (Str :$filein){
	
	my $projectName = $self->project->name();
	my $logPlink = $self->project->getProjectPath().'../'.$projectName.'.plink.resume';
	my $ppn =$self->nproc;
	my $tfilein;
	my $coverage_dir = $self->project->getRootDir() . "/align/coverage/depth/";
	foreach my $patient (@{$self->project->getPatients}){
		my $file   = $patient->fileNoSqlDepth;
		die("file coverage not found : $file") unless -e $file;
		push(@$tfilein,$file);
	}
	
	my $type = "dude";
	my $stepname = $projectName."@".$type;
	my $dir = $self->project->project_log();
	my $final_dir = $self->project->getVariationsDir("dude");
	my $fileout = $final_dir."/".$self->project->name.".dude";
	
	my $cmd = "perl $bin_dev/dude/dude.pl -project=$projectName  -fork=$ppn && date > $fileout";
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>\@$tfilein,fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
		#@@@@@@@@@@@@@@@@@@@@@
		#todo check for all individual file 
	if ($self->unforce() && -e $fileout){
	  		$job_bds->skip();
	}
	return ($filein);
	
}


method quality_check (Str :$filein){
	
	my $projectName = $self->project->name();
	my $logPlink = $self->project->getProjectPath().'../'.$projectName.'.plink.resume';
	my $ppn = $self->nproc;
	unless ($filein){
		$filein = $self->project->getVariationsDir($self->method_calling)."/".$projectName.".vcf.gz";
		die() unless -e $filein;
	}
	my $type = "qualtity-check";
	my $stepname = $projectName."@".$type;
	my $dir = $self->project->project_log();
	my $fileout = $dir."/quality_check.log";
		my $cmd = "perl $bin_dev/quality_check.pl -project=$projectName  -vcf_file=$filein -fork=$ppn >$fileout";
	my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	$job_bds->isLogging(1);
		#@@@@@@@@@@@@@@@@@@@@@
		#todo check for all individual file 
	if ($self->unforce() && -e $fileout){
	  		$job_bds->skip();
	}
	return ($filein);
	
}

method cache (Str :$filein){
	
	my $projectName = $self->project->name();
	my $logPlink = $self->project->getProjectPath().'../'.$projectName.'.plink.resume';
	my $ppn = $self->nproc;
	if ($self->nocluster){
		$ppn = $self->nproc;
	}
	
	$ppn = int($ppn/4);
#	unless ($filein){
#		$filein = $self->project->getVariationsDir($self->method_calling)."/".$projectName.".vcf.gz";
#		die() unless -e $filein;
#	}

	my $cmds= $self->get_list_cmds_vector(fork=>$ppn);
	
 my $fout;	
 my $type = "cacheV";
	foreach my $hcmd (@$cmds){
			
		my $stepname = $projectName.".".$hcmd->{chr}."@".$type;
		my $cmd = $hcmd->{cmd};
		my $fileout = $hcmd->{fileout};
		push(@$fout,$fileout);
		my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
		$self->current_sample->add_job({job=>$job_bds});
		if ($self->unforce() && -e $fileout){
	  		$job_bds->skip();
		}
		
	}
	
	my $fileout_gv = $self->project->getCacheBitVectorDir()."/log/$projectName.global.log";
	my $ppn2 =  $self->ppn;
		if ($self->nocluster){
		$ppn2 = $self->nproc;
	}
    my $cmd2 = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/cache.pl -fork=$ppn2 -project=$projectName -type=global -force=1 -bds" ;
	my $type2 = "global_infos";
	my $stepname = $projectName."@".$type2;
	my $job_bds = job_bds->new(cmd=>[$cmd2],name=>$stepname,ppn=>$ppn2,filein=>$fout,fileout=>$fileout_gv,type=>$type2,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	

	
	return ($filein);
}

method cache_test (Str :$filein){
	
	my $projectName = $self->project->name();
	my $logPlink = $self->project->getProjectPath().'../'.$projectName.'.plink.resume';
	my $ppn = $self->ppn;
	if ($self->nocluster){
		$ppn = $self->nproc;
	}
	
	$ppn = int($ppn/4);
#	unless ($filein){
#		$filein = $self->project->getVariationsDir($self->method_calling)."/".$projectName.".vcf.gz";
#		die() unless -e $filein;
#	}

	my $cmds= $self->get_list_cmds_vector(fork=>$ppn);
	
 my $fout;	
 my $type = "cacheV";
	foreach my $hcmd (@$cmds){
			
		my $stepname = $projectName.".".$hcmd->{chr}."@".$type;
		my $cmd = $hcmd->{cmd_test};
		my $fileout = $hcmd->{fileout};
		push(@$fout,$fileout);
		my $job_bds = job_bds->new(cmd=>[$cmd],name=>$stepname,ppn=>$ppn,filein=>[$filein],fileout=>$fileout,type=>$type,dir_bds=>$self->dir_bds);
		$self->current_sample->add_job({job=>$job_bds});
		if ($self->unforce() && -e $fileout){
	  		$job_bds->skip();
		}
		
	}
	
	my $fileout_gv = $self->project->getCacheBitVectorDir()."/log/$projectName.global.log";
	my $ppn2 =  $self->ppn;
		if ($self->nocluster){
		$ppn2 = $self->nproc;
	}
    my $cmd2 = "echo toto" ;
	my $type2 = "global_infos";
	my $stepname = $projectName."@".$type2;
	my $job_bds = job_bds->new(cmd=>[$cmd2],name=>$stepname,ppn=>$ppn2,filein=>$fout,fileout=>$fileout_gv,type=>$type2,dir_bds=>$self->dir_bds);
	$self->current_sample->add_job({job=>$job_bds});
	

	
	return ($filein);
}

method get_list_cmds_vector (Int :$fork) {
	my $project =  $self->project;
	my $project_name = $project->name();
	mkdir $project->getCacheBitVectorDir().'/log' unless (-d $project->getCacheBitVectorDir().'/log');
	my  @lChr = (1..22, 'X', 'Y', 'MT'); 
	my $cmds;
	foreach my $chr_name (@lChr) {
		my $hcmds;
		$hcmds->{fileout} = $project->getCacheBitVectorDir().'/log/cache_'.$chr_name."_P1.log";
		$hcmds->{fileout2} = $project->getCacheBitVectorDir().'/log/cache_'.$chr_name."_P2.log";
			$hcmds->{cmd} = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/cache.pl -fork=$fork -project=$project_name -chr=$chr_name -type=part1 -force=1 -bds " ;
		$hcmds->{cmd_test} = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/cache.pl -fork=$fork -project=$project_name -chr=$chr_name -type=part1 -force=1 -test  " ;
		$hcmds->{cmd2} = "/usr/bin/perl $Bin/../polymorphism-cgi/cache_nodb/cache.pl -fork=$fork -project=$project_name -chr=$chr_name -type=denovo -force=1 -bds  " ;
	#	warn $hcmds->{cmd};
		$hcmds->{chr} =  $chr_name;
		push(@$cmds,$hcmds)
	}
	return ($cmds);
}


1;
