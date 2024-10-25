package calling_steps;
use strict;
use FindBin qw($Bin);
use lib "$Bin";
#use root_steps;
use Moose;  
use MooseX::Method::Signatures;
use Data::Dumper;

extends (qw(root_steps));


my $bin_dev = qq{$Bin/scripts/scripts_pipeline/};

my $ext_recal = {
	snp=> {
		  tranches=>"snp.tranches",
		  recal=>"snp.recal",
	},
	indel=> {
		  tranches=>"indel.tranches",
		  recal=>"indel.recal",
	}
};

has 'method_calling' => (
	is =>'rw',
		required =>1,
);


has 'project' => (
	is =>'rw',
	isa =>"Object",
);
has 'somatic' => (
	is =>'rw',
	default => undef,
);
has 'bam_files'  =>(
	is => 'rw',
	isa => 'Str',
	lazy => 1,
	default => \&init_bam_files,
);

has 'vcftools_path'  =>(
	is => 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		return $self->project->getSoftware('vcftools_path');
	},
);

has 'javac' => (
	is => 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		return $self->project->getSoftware('java').' -Xmx20g -XX:ParallelGCThreads=8';
	}
);

has 'gatk' => (
	is => 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		return $self->javac.' -jar '.$self->project->getSoftware('gatk');
	}
);

has 'public_data' => (
	is => 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		my $build = $self->project->getVersion();
		return $self->project->buffer->config->{public_data}->{$build}.'/';
	}
);

has 'reference'  => (
	is        => 'ro',
	isa       => 'Str',
	lazy =>1,
	default   => sub {
		my $self = shift;
		return $self->public_data()."/genome/fasta/all.fa";
	},
); 



method correct_vcf(Str :$filein! ,Object :$previous! ){
	my $project = $self->project;
	my $name = $project->name();
	my $dir_out= $project->getCallingPipelineDir($self->method_calling);
	my $snp_out = $dir_out . "/" . $project->name . ".snps.vcf";
	my $indels_out = $dir_out . "/" . $project->name . ".indels.vcf";
		my $stepname = $name."_correct";
	if ($self->unforce() && -e $snp_out && -e $indels_out){
	  	$self->add_skip_steps($stepname);
		return ($previous,$snp_out);
	}
	#my $fileout = $filein;
	my $fileout =$dir_out."/".$name.".snps.vcf";
	warn $filein;
	unless ($filein){
		$filein = $dir_out."/$name.uni.vcf";
		die("can't find $filein ") unless -e $filein;
	}
	my $fork = 1;
	my $patient =  $self->argument_patient(); 
	 my $cmd = "perl $bin_dev/correct_gatk.pl -vcf=$filein -project=$name -dir=$dir_out -snp_out=$snp_out -indels_out=$indels_out -patient=$patient -log=".$self->log_file." -method=".$self->method_calling;
	my ($job,$job_next) =  $self->construct_jobs_bds(stepname=>$stepname,cmd=>$cmd,ppn=>$fork,type=>"wait_exit");
	$self->add_running_steps($stepname);
	  if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	  return ($job,$fileout);	
}

method split_vcf (Str :$filein! ,Object :$previous! ){
	my $project = $self->project;
	my $name = $project->name();
	my $dir_out= $project->getCallingPipelineDir("unifiedgenotyper");
	my $fileout = $dir_out."/".$name.".snps.vcf";
	my $tmpout = $dir_out."/".$name.".tmp.snps.vcf";
	my $tmpout_indels = $dir_out."/".$name.".tmp.indels.vcf";
	my $fileout_indels= $dir_out."/".$name.".indels.vcf";
	my $stepname = $name."_split";

	my @patients_name = map {$_->name()} @{$self->project->get_list_patients($self->argument_patient())};
	my $sn;
	my @query = ("QD > 5.0 ");
	my @query2;
	my @query3 ;
	foreach my $pname (@patients_name){
		$sn .= "-sn $pname ";
		push(@query, 
		qq{vc.getGenotype("$pname").getAD().1 > 7 });
		push(@query2, 
		qq{ vc.getGenotype("$pname").getAD().1 > 3});
		push(@query3, 
		qq{ vc.getGenotype("$pname").getAD().1 > 5});
	}
	
	my $st_query = join(" || ",@query);
	my $st_query2 = join(" || ",@query2);
	my $st_query3 = join(" || ",@query3);
	#my $cmd = $self->gatk()."-T SelectVariants -R ".$self->reference()." -selectType SNP -o $fileout --variant $filein $sn  -select '((QD>3.0 && MQ >40.0 && FS < 60.0 && HaplotypeScore < 13.0 && MQRankSum > -12.5 && ReadPosRankSum > -8.0) || (".$st_query."))&&($st_query2)'";
	my $cmd = $self->gatk()."-T SelectVariants -R ".$self->reference()." -selectType SNP -o $tmpout --variant $filein $sn  -select '".$st_query."'";
	my $cmd3 = $self->gatk()."-T SelectVariants -R ".$self->reference()." -selectType SNP -o $fileout --variant $tmpout -select '".$st_query2."'";
	$cmd = $cmd.";".$cmd3;
	#QD < 2.0
	#ReadPosRankSum < -20.0
	#InbreedingCoeff < -0.8
	#FS > 200.0
	my $select_indels = "QD > 2.0 || ReadPosRankSum > -20.0 || InbreedingCoeff > -0.8 || FS < 200.0";
	my $cmd2 = $self->gatk()."-T SelectVariants -R ".$self->reference()." -selectType INDEL -o $tmpout_indels --variant $filein $sn -select '".$st_query."'";;
	my $cmd4 = $self->gatk()."-T SelectVariants -R ".$self->reference()." -selectType INDEL -o $fileout_indels --variant $tmpout_indels -select '".$st_query3."'";
	$cmd2 = $cmd2.";".$cmd4;
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd.";".$cmd2,ppn=>1);
	$self->add_running_steps($stepname);
	  if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	  return ($job,$fileout);	
}

method calling_unifiedgenotyper (Str :$filein! ,Object :$previous! ){
	my $project = $self->project;
	my $name = $project->name();
	my @jobs;
	my @files;
	mkdir ( $project->getCallingPipelineDir("unifiedgenotyper")) unless -e  $project->getCallingPipelineDir("unifiedgenotyper");
	my $dir_out= $project->getCallingPipelineDir("unifiedgenotyper");
	my $fileout = $dir_out."/".$name.".uni.vcf";
	my $stepname = $name."_catvcf";
	if ($self->unforce() && -e $fileout ){
	  	$self->add_skip_steps($stepname);
		return ($previous,$fileout);
	}
	my $chrs = [1..22,'X','Y'];
	foreach my $chr_name (@$chrs){
		#my $chr = $self->project->getChromosome($chr_name);
		#next unless $chr;
		my ($job,$file) = $self->calling_unifiedgenotyper_chr(chr=>$chr_name,previous=>$previous);
		if ($job){
			push(@jobs,$job);
		}
		push(@files,$file);
	}
	my $cmd = "perl $bin_dev/concat_vcf.pl -fork=16 -project=".$name." -file_out=$fileout -log=".$self->log_file;#." -norecal=1";
	my ($job,$job_next) = $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>16);
	if (scalar @jobs){
		foreach my $pre (@jobs){
				$self->prev(job=>$job,prev=>$pre);
		}
	}
	else {
		$self->prev(job=>$job,prev=>$previous);
	}
	$self->add_running_steps($stepname);
	return ($job,$fileout);
}

method unifiedgenotyper(Any :$chr!,Any :$previous){
	my $project = $self->project;
	my $ped_file = $project->getPedigreeFile();
	
	
}


#method calling_unifiedgenotyper_chr (Any :$chr!,Any :$previous){ 
#	my $name = "chr".$chr;
#	my $project = $self->project;
#	my $pname = $project->name();
#	mkdir ( $project->getCallingPipelineDir("unifiedgenotyper")) unless -e  $project->getCallingPipelineDir("unifiedgenotyper");
#	my $dir_out= $project->getCallingPipelineDir("unifiedgenotyper");
#	my $ext = "uni";
#	my $fileout  = $dir_out."/".$name.".$ext.vcf";
#	my $ppn =16;
#	my $zz = $chr;
#	$zz = 23 if $zz eq 'X';
#	$zz = 24 if $zz eq 'Y';
#	$zz = 25 if $zz eq 'MT';
#	$zz = 25 if $zz eq 'M';
#	if ($zz < 4 ){ 
#		$ppn = 24;
#	}
#	my $arg = $self->argument_patient(); 
#	my $arg_somatic =  "";
#	 $arg_somatic =  " -somatic=1 " if $self->somatic(); 
#	 my $cmd = "perl $bin_dev/fast_unifiedGenotyper.pl -patient=$arg -project=$pname -chr=$name -fork=$ppn -ext=$ext $arg_somatic -log=".$self->log_file;
#	# my $cmd = "perl $bin_dev/fastest_haplotypecaller.pl   -patient=$arg -project=$pname -chr=$name -fork=$ppn -ext=$ext -log=".$self->log_file;
#	 my $stepname = $name."_unified";
#	 if ($self->unforce() && -e $fileout && -s $fileout > 0){
#	  	$self->add_skip_steps($stepname);
#		return (undef,$fileout);
#	}
#	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$ppn);
#		$self->add_running_steps($stepname);
#	  if (defined $previous) {
#		$self->prev(job=>$job,prev=>$previous);
#	}	
#	  return ($job,$fileout);
#	 
#}


method calling_larges_indels  (Str :$filein! ,Object :$previous! ){
	my $project = $self->project;
	my $project_name = $project->name();
	my @jobs;
	my @files;
	
	my $dir_out= $project->getCallingPipelineDir("gvcf");
	my $arg = $self->argument_patient(); 
	my $patients = $project->get_list_patients($arg);
	my $fork =16;

		warn "coucou";
	foreach my $patient (@$patients){
		my $fileout =  $patient->getCnvsFile();
		next   if ($self->unforce() && -e $fileout);
		my $name = $patient->name;
		my $log_file = $self->log_file;
		my $stepname = $name."_lifinder";
		
		my $cmd = qq{perl $bin_dev/calculate_large_indels.pl -patient=$name -fork=16 -project=$project_name -log=$log_file };
			my ($job,$job_next) = $self->construct_jobs_bds(stepname=>$stepname,cmd=>$cmd,ppn=>$fork,type=>"run");
			$self->add_running_steps($stepname);
	}
	
	return ($previous," ");
}

method count_featureCounts  (Str :$filein! ,Object :$previous! ){
	my $project = $self->project;
	my $project_name = $project->name();
	my @jobs;
	my @files;
	
	my $dir_out= $project->getCountingDir("featureCounts");
	my $fileout = $dir_out."/".$project_name.".count.genes.txt";
	my $fileout2 = $dir_out."/".$project_name.".count.exons.txt";
	my $arg = $self->argument_patient(); 
	my $patients = $project->get_list_patients($arg);
	my @bams;
	my @sed_cmd;
	my @sed_cmd2;
	my $align_method;
	
	foreach my $patient (@$patients){
		my $bam = $patient->getBamFile;
		my $name = $patient->name;
		push(@bams,$bam);
		$bam =~ s/\//\\\//g;
		$align_method = $patient->alignmentMethod();
		push(@sed_cmd,qq{sed -i "2s/$bam/$name/" $fileout} );
		push(@sed_cmd2,qq{sed -i "2s/$bam/$name/" $fileout2} );
	}
	my $fork =16;
	my $gtf = $project->gtf_file();
	#$gtf = $project->gtf_file_star() if $align_method eq "star";

	my $sed = join(" && ",@sed_cmd);
	my $featureCounts = $project->buffer->software("featureCounts");
	my $cmd = "$featureCounts -T $fork -a $gtf  --ignoreDup  -o $fileout -p -t exon -s 1".join(" ",@bams)."&& $sed";
#	warn $cmd;
	my $stepname = "featureCounts_genes";
		my ($job,$job_next) = $self->construct_jobs_bds(stepname=>$stepname,cmd=>$cmd,ppn=>$fork,type=>"run");
		$self->add_running_steps($stepname);
			my $sed2 = join(" && ",@sed_cmd2);
	my $cmd2 = "$featureCounts -T $fork -a $gtf -f  -p  -t exon  --ignoreDup  --minOverlap 10 -o $fileout2 -s 1  ".join(" ",@bams)." && $sed2";
		my $stepname2 = "featureCounts_exons";
		my ($job1,$job_next1) = $self->construct_jobs_bds(stepname=>$stepname2,cmd=>$cmd2,ppn=>$fork,type=>"run");
		$self->add_running_steps($stepname2);
	return ($previous," ");
}



method calling_haplotypecaller_gvcf (Str :$filein! ,Object :$previous! ){
	warn $filein."-".$previous;
	my $project = $self->project;
	my $name = $project->name();
	my @jobs;
	my @files;
	
	my $dir_out= $project->getCallingPipelineDir("gvcf");
	my $arg = $self->argument_patient(); 
	my $patients = $project->get_list_patients($arg);
	my $list_patients= [];
	# test if global vcf exists !!!!
	 
	my $dir_out_vcf= $project->getCallingPipelineDir("haplotypecaller");
	
	my $stepname = $name."_join.gvcf";
	
	if ($self->unforce()){
		foreach my $patient (@$patients){
			my $gvcf = $patient->getGvcfFile(); #$dir_out."/".$patient->name.".g.vcf.gz";
			warn "skip : ".$patient->name." g.vcf exists " if  $gvcf;
			$self->add_skip_steps("gvcf_".$patient->name) if $gvcf;
			push(@$list_patients,$patient->name) unless  $gvcf;
		}
		
	}
	else {
		foreach my $patient (@$patients){
			my $gvcf = $dir_out."/".$patient->name.".g.vcf.gz";
			unlink $gvcf if -e $gvcf;
				push(@$list_patients,$patient->name);
		}
		
	}
	unless (scalar (@$list_patients)){
	 	#$self->add_skip_steps($stepname);
		return ($previous,"undef");
	}
	
	my $size_window = 1_000_000;
	my $total_jobs =0;
	foreach my $chr (@{$project->getChromosomes()}){
		next unless scalar(@$list_patients);
		my $windows = $chr->getWindowCaptureForCalling(250,$size_window);
	
		foreach my $window (@$windows){
				next if $window->{intspan}->is_empty();
				my ($job,$file) = $self->calling_hc_region(chr=>$chr,previous=>$previous,start=> $window->{start},end=> $window->{end},window=>$window,list=>$list_patients);
				
				if ($job){
				
					push(@jobs,@$job);
				}
		}
		warn "end chr ".$chr->name." :: ".scalar(@jobs);
	
		#last;
	}
		$total_jobs = scalar(@jobs);
    warn "total jobs : $total_jobs";
	my @jobs_concat;
	my $cmd = qq{echo 'END CALLING'};

	if (scalar @jobs){
			my ($job_echo,$job_next) =  $self->construct_jobs_bds(stepname=>"END_CALLING",cmd=>$cmd,ppn=>1,type=>"wait");
		foreach my $pre (@jobs){
				$self->prev(job=>$job_echo,prev=>$pre);
		}
		return ($job_echo," ");	
	}
	else {
		return ($previous,"");
		#$self->prev(job=>$job_echo,prev=>$previous);
	}
	
 return ($previous," ");	

	
	
}


method join_gvcf (Str :$filein! ,Object :$previous! ){
	my $project = $self->project;
	my $name = $project->name();
	my @jobs;
	my @files;
	my $size_window = 1_000_000;
	my $dir_out= $project->getCallingPipelineDir("gvcf");
	my $arg = $self->argument_patient(); 
	my $patients = $project->get_list_patients($arg);
	my $list_patients= [];
	my $fileout = 
	# test if global vcf exists !!!!
	 
	my $dir_out_vcf= $project->getCallingPipelineDir("haplotypecaller");
	
	my $stepname = $name."_join.gvcf";
	
	if ($self->unforce()){
		foreach my $patient (@$patients){
			my $gvcf = $patient->getGvcfFile(); #$dir_out."/".$patient->name.".g.vcf.gz";
			warn "skip : ".$patient->name." g.vcf exists " if  $gvcf;
			$self->add_skip_steps("gvcf_".$patient->name) if $gvcf;
			push(@$list_patients,$patient->name) unless  $gvcf;
		}
		
	}
	else {
		foreach my $patient (@$patients){
			my $gvcf = $dir_out."/".$patient->name.".g.vcf.gz";
			unlink $gvcf if -e $gvcf;
				push(@$list_patients,$patient->name);
		}
		
	}
	#unless (scalar (@$list_patients)){
	 	#$self->add_skip_steps($stepname);
	#	return ($previous,"undef");
	#}
	my @jobs_concat;

	foreach my $patient_name (@$list_patients){
		my $fork = 8;
		my $cmd = "perl $bin_dev/join_gvcf.pl -project=".$name." -vcf=$fileout -window=$size_window -patient=$patient_name -fork=$fork";
		
		my ($job,$job_next) = $self->construct_jobs_bds(stepname=>"join_gvcf_".$patient_name,cmd=>$cmd,ppn=>$fork,type=>"run");
			$self->add_running_steps("join_gvcf_".$patient_name);
		if ($job){
					push(@jobs_concat,$job);
				}
		
				$self->prev(job=>$job,prev=>$previous);
	}
	my $cmdf = qq{echo "END GVCF"};
	
	
	
	

	if (scalar @jobs_concat){
		my ($final_gvcf_job,$job_next2) =   $self->construct_jobs_bds(stepname=>"END JOIN GVCF ",cmd=>$cmdf,ppn=>1,type=>"wait");
#			$final_gvcf_job->prev( { ok => \@jobs_concat} );
			#$previous->next({ok=>[@jobs_concat]});
		#	warn Dumper $final_gvcf_job;
		foreach my $pre (@jobs_concat){
				$self->prev(job=>$final_gvcf_job,prev=>$pre);
		}
		return ($final_gvcf_job,"toto");
	}
	else {
		return ($previous,"");
		#$final_gvcf_job->prev( { ok => [$previous]} );
		
		#$self->prev(job=>$final_gvcf_job,prev=>$previous);
	}
	$self->add_running_steps("END JOIN GVCF");
	
	
}
method genotype_gvcf4 (Str :$filein! ,Object :$previous! ){
	my $project = $self->project;
	my $name = $project->name();
	my $arg = $self->argument_patient(); 
	my $dir_out_vcf= $project->getCallingPipelineDir("haplotypecaller");
	my $fileout = $dir_out_vcf."/".$name.".hc.vcf";
	my $stepname = "genotype_gvcf";
	if ($self->unforce() && -e $fileout){
	  	$self->add_skip_steps($stepname);
		return ($previous,$fileout);
	}
	my $patients = $project->get_list_patients($arg);
	foreach my $patient (@$patients){
			my $gvcf = $patient->getGvcfFile();#$dir_out."/".$patient->name.".g.vcf.gz";
			die("you don't have gvcf for at least this patient restart after_lifescope on this project ") unless $gvcf ;#-e $gvcf;
	}
	my $fork =20;
	my $size_window = 5_000_000;
	my $cmd = "perl $bin_dev/gatk4/genotype_gvcf.pl -project=".$name." -vcf=$fileout -window=$size_window -patient=$arg -fork=$fork -log=".$self->log_file;#." -norecal=1";
	my ($job,$job_next) = $self->construct_jobs_bds(stepname=>$stepname,cmd=>$cmd,ppn=>$fork,type=>"wait_exit");

#
		$self->prev(job=>$job,prev=>$previous);
#
	$self->add_running_steps($stepname);
	return ($job,$fileout);
	
}
method genotype_gvcf (Str :$filein! ,Object :$previous! ){
	my $project = $self->project;
	my $name = $project->name();
	my $arg = $self->argument_patient(); 
	my $dir_out_vcf= $project->getCallingPipelineDir("haplotypecaller");
	my $fileout = $dir_out_vcf."/".$name.".hc.vcf";
	my $stepname = "genotype_gvcf";
	if ($self->unforce() && -e $fileout){
	  	$self->add_skip_steps($stepname);
		return ($previous,$fileout);
	}
	my $patients = $project->get_list_patients($arg);
	foreach my $patient (@$patients){
			my $gvcf = $patient->getGvcfFile();#$dir_out."/".$patient->name.".g.vcf.gz";
			die("you don't have gvcf for at least this patient restart after_lifescope on this project ") unless $gvcf ;#-e $gvcf;
	}
	my $fork =16;
	my $size_window = 1_000_000;
	my $cmd = "perl $bin_dev/genotype_gvcf.pl -project=".$name." -vcf=$fileout -window=$size_window -patient=$arg -fork=$fork -log=".$self->log_file;#." -norecal=1";
	my ($job,$job_next) = $self->construct_jobs_bds(stepname=>$stepname,cmd=>$cmd,ppn=>$fork,type=>"wait_exit");

#
		$self->prev(job=>$job,prev=>$previous);
#
	$self->add_running_steps($stepname);
	return ($job,$fileout);
	
}



method calling_hc_region (Any :$chr!,Any :$previous,Any :$start!,Any :$end!,Any :$window!,Any :$list!){ 
	my $name = $chr->ucsc_name."-$start-$end";
	my $project = $self->project;
	my $pname = $project->name();
	mkdir ( $project->getCallingPipelineDir("haplotypecaller")) unless -e  $project->getCallingPipelineDir("haplotypecaller");
	my $dir_out= $project->getCallingPipelineDir("gvcf");
	#mkdir $dir_out."/vcf/" unless -e $dir_out."/vcf/";
	#my $patients = $project->get_list_patients(join(",",@$list);
	 my $arg = join(",",@$list);

	
	
	

	
	 my $chr_name = $chr->ucsc_name;
	my $ext ="";
	 #my $cmd = "perl $bin_dev/haplotypeCaller_gvcf.pl -patient=$arg -project=$pname -chr=$chr_name -start=$start -end=$end -fork=16  -log=".$self->log_file;
	# my $cmd = "perl $bin_dev/haplotypeCaller.pl   -patient=$arg -project=$pname -chr=$name -fork=$ppn -ext=$ext -log=".$self->log_file;
	 my $stepname = $name."_hc";
	 my $pp;
	 my $outfiles = [];
	 if ($self->unforce()){
	 			my $patients = $project->get_list_patients($arg);
	 			my $no_file;
	 			my @names;
	 			foreach my $patient (@$patients){
	 				my $outfile = $patient->getWindowGvcf($window);
	 				
	 			
	 				if (-e $outfile && -s $outfile > 0){
	 					next;
	 				}
	 				else {
	 					unlink $outfile if -e $outfile;
	 					$no_file =1;
	 					push(@$pp,$patient);
	 					push(@names,$patient->name());
	 					push(@$outfiles,$outfile);
	 					#last;
	 				}
	 				
	 			}
#	 			
	 			unless($no_file){
	 					#$self->add_skip_steps($stepname);
	 					return (undef,undef);
	 			}
	 			 $arg = join(",",@names);
	 	
	 }
	# if ($self->unforce() && -e $fileout && -s $fileout > 0){
	 	
	  #	$self->add_skip_steps($stepname);
	#	return (undef,$fileout);
	#}
	
	my @all_jobs;
	my $nb_p = scalar(split(",",$arg));
	return (undef,undef) if $nb_p == 0;
	my $fork = 16;
	if ($nb_p < 16){
		$fork = $nb_p *2;
	}
	$fork = 16 if $fork >16 ;
#	$fork =2;
	my $real_fork = int($fork /2);
	#foreach my $p (@$pp){
	#	my $pname = $p->name();
	# my $stepname = $name."_hc.$pname";
	 my $cmd = "perl $bin_dev/haplotypeCaller_gvcf.pl -patient=$arg -project=$pname -chr=$chr_name -start=$start -end=$end -fork=$real_fork ";
		
	my ($job,$job_next) =  $self->construct_jobs_bds(stepname=>$stepname,cmd=>$cmd,ppn=>$fork,type=>"run");
	
		$self->add_running_steps($stepname);
	  if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}
	push(@all_jobs,$job);
	#}
	  return (\@all_jobs,undef);
	 
}



method snp_recalibrator (Str :$filein! ,Object :$previous! ){
	my $name = $self->project->name();
	my $dir_out= $self->project->getPipelineCallingDir("unifiedgenotyper");
	my $tranches  = $dir_out."/".$name.$ext_recal->{snp}->{tranches};
	my $recal  =  $dir_out."/".$name..$ext_recal->{snp}->{recal};
	my $fileout = $filein;
	$fileout =~ s/vcf/recal_snp\.vcf/;
	my $cmd = $self->gatk()." -T VariantRecalibrator -R ".$self->reference()." -input $filein -recalFile $recal -tranchesFile $tranches ";
	$cmd .= " --maxGaussians 6 ";
	$cmd .= " -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ".$self->project->gatk_hapmap_file() if ($self->project->gatk_hapmap_file());
	$cmd .= " -resource:omni,known=false,training=true,truth=false,prior=12.0 ".$self->project->gatk_omni_file() if ($self->project->gatk_omni_file());
	$cmd .= " -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 ".$self->project->gatk_dbsnp_file() if ($self->project->gatk_dbsnp_file());
	$cmd .= " -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ  -an DP -mode SNP ";
	#-an InbreedingCoeff
	my $stepname = $name."_snp_recal1";
	my $ppn =2;
	if ($self->unforce() && -e $fileout){
	  	$self->add_skip_steps($stepname);
		return ($previous,$filein);
	}
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$ppn);
		$self->add_running_steps($stepname);
	  if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	
	  return ($job,$filein);
	 
}

method snp_apply_recalibration (Str :$filein! ,Object :$previous! ){
	my $name = $self->project->name();
	my $dir_out= $self->project->getPipelineCallingDir("unifiedgenotyper");
	
	my $tranches  = $dir_out."/".$name.$ext_recal->{snp}->{tranches};
	my $recal  =  $dir_out."/".$name..$ext_recal->{snp}->{recal};
	my $fileout = $filein;
	$fileout =~ s/vcf/recal_snp\.vcf/;
	#my $fileout = $dir_out."/".$name.".final.vcf";
	my $cmd = $self->gatk()." -T ApplyRecalibration -R ".$self->reference()." -input $filein -tranchesFile $tranches -recalFile $recal -o $fileout "; 
	 $cmd .= qq{--ts_filter_level 95.0 -mode SNP};
	 
	 my $stepname = $name."snp_recal2";
	my $ppn =2;
	 if ($self->unforce() && -e $fileout){
	  	$self->add_skip_steps($stepname);
		return ($previous,$fileout);
	}
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$ppn);
		$self->add_running_steps($stepname);
	  if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	  return ($job,$fileout);
	 
	
}

method indels_recalibrator (Str :$filein! ,Object :$previous! ){
	
	my $name = $self->project->name();
	my $dir_out= $self->project->getPipelineCallingDir("unifiedgenotyper");
	my $tranches  = $dir_out."/".$name.$ext_recal->{indel}->{tranches};
	my $recal  =  $dir_out."/".$name.$ext_recal->{indel}->{recal};

	my $cmd = $self->gatk()." -T VariantRecalibrator -R ".$self->reference()." -input $filein -recalFile $recal -tranchesFile $tranches";
	 $cmd .= " --maxGaussians 4 -std 10.0 -percentBad 0.12";
	 $cmd .= " -resource:mills,known=true,training=true,truth=true,prior=12.0 ".$self->project->gatk_mills_indels_file() if ($self->project->gatk_mills_indels_file());
	 $cmd .= " -an QD -an FS -an HaplotypeScore -an ReadPosRankSum -an InbreedingCoeff";
	 $cmd .= " -mode INDEL ";
	#-an InbreedingCoeff
	my $stepname = $name."indel_recal1";
	
	my $ppn =2;
	my $fileout = $filein;
	 $fileout =~ s/vcf/recal_indel\.vcf/; 
	 warn $fileout;
	 if ($self->unforce()  && -e $tranches && -e $recal){
	  	$self->add_skip_steps($stepname);
		return ($previous,$filein);
	}
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$ppn);
		$self->add_running_steps($stepname);
	  if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	  return ($job,$filein);
	 
}



method indels_apply_recalibration (Str :$filein! ,Object :$previous! ){
	my $name = $self->project->name();
	
	my $dir_out= $self->project->getPipelineCallingDir("unifiedgenotyper");
	my $tranches  = $dir_out."/".$name.$ext_recal->{indel}->{tranches};
	my $recal  =  $dir_out."/".$name.$ext_recal->{indel}->{recal};
	
	my $fileout = $filein;
	 $fileout =~ s/vcf/recal_indel\.vcf/;
	my $cmd = $self->gatk()." -T ApplyRecalibration -R ".$self->reference()." -input $filein -tranchesFile $tranches -recalFile $recal -o $fileout "; 
	 $cmd .= qq{-ts_filter_level 95.0 -mode INDEL};
	 my $stepname = $name."indel_recal2";
	my $ppn =2;
	 if ($self->unforce() && -e $fileout){
	  	$self->add_skip_steps($stepname);
		return ($previous,$fileout);
	}
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$ppn);
		$self->add_running_steps($stepname);
	  if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	  return ($job,$filein);
}

method mendelian_verification (Str :$filein! ,Object :$previous! ){
	my $name = $self->project->name();
	
	my $dir_out= $self->project->getCallingPipelineDir("unifiedgenotyper");
	
	my $fileout = $filein;
	
	 my $cmd = "perl $bin_dev/mendelian_error.pl -vcf=$filein -project=$name -log=".$self->log_file;
	 my $ppn = 1;
	 my $stepname = $name."_mendel";
	
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$ppn);
		$self->add_running_steps($stepname);
	  if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	  return ($job,$fileout);

}


method sex_verification (Str :$filein! ,Object :$previous! ){
	my $name = $self->project->name();
	
	
	
	my $fileout = $filein;
	
	 my $cmd = "perl $bin_dev/check_sex.pl -project=$name -log=".$self->log_file;
	 my $ppn = 1;
	 my $stepname = $name."_sex";
	
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$ppn);
		$self->add_running_steps($stepname);
	  if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	  return ($job,$fileout);

}


method coverage_verification (Str :$filein! ,Object :$previous! ){
	my $name = $self->project->name();
	my $fileout = $filein;
	my $cmd = "perl $bin_dev/check_coverage.pl -project=$name -log=".$self->log_file;
	my $ppn = 1;
	my $stepname = $name."_coverage_verif";
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$ppn);
	$self->add_running_steps($stepname);
	if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	return ($job,$fileout);
}

method check_dbsnp (Str :$filein! ,Object :$previous! ){
	my $name = $self->project->name();
	
	my $dir_out= $self->project->getCallingPipelineDir("unifiedgenotyper");
	
	my $fileout = $filein;
	my $arg = $self->argument_patient(); 
	
	 my $cmd = "perl $bin_dev/check_dbsnp.pl -patient=$arg -project=$name -vcf=$filein -log=".$self->log_file;
	 my $ppn = 1;
	 my $stepname = $name."_dbsnp";
	
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$ppn);
		$self->add_running_steps($stepname);
	  if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	  return ($job,$fileout);
}

method move_vcf (Str :$filein! ,Object :$previous! ){
	my $name = $self->project->name();
	my $dir_out= $self->project->getCallingPipelineDir($self->method_calling);
	my $arg = $self->argument_patient(); 
	my $cmd = "perl $bin_dev/move_vcf.pl  -project=$name -vcf_dir=$dir_out -patient=$arg -method_calling=".$self->method_calling." -log=".$self->log_file;
	my $stepname = $name."_move_vcf";
	my $ppn = 1;
	my ($job,$job_next) =  $self->construct_jobs_bds(stepname=>$stepname,cmd=>$cmd,ppn=>$ppn,type=>"wait_exit");
	$self->add_running_steps($stepname);
	if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	return ($job,$filein);
}

method quality_check (Str :$filein! ,Object :$previous! ){
	
	my $projectName = $self->project->name();
	my $logPlink = $self->project->getProjectPath().'../'.$projectName.'.plink.resume';
	my $ppn = 16;
	my $cmd = "perl $bin_dev/quality_check_project.pl -project=$projectName  -vcf_file=$filein -log=".$self->log_file;
	my $stepname = $projectName."_quality_check";
	my ($job, $job_next) = $self->construct_jobs_bds(stepname=>$stepname, cmd=>$cmd, ppn=>$ppn, type=>"run");
	$self->add_running_steps($stepname);
	if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	return ($job, $filein);
}

method plink (Str :$filein! ,Object :$previous! ){
	my $projectName = $self->project->name();
	my $logPlink = $self->project->getProjectPath().'../'.$projectName.'.plink.resume';
	my $ppn = 16;
	my $cmd = "perl $bin_dev/plink.pl -project=$projectName -fork=$ppn -vcf_file=$filein -log=".$self->log_file;
	my $stepname = $projectName."_plink";
	my ($job, $job_next) = $self->construct_jobs(stepname=>$stepname, cmd=>$cmd, ppn=>$ppn);
	$self->add_running_steps($stepname);
	if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	return ($job, $filein);
}


method move_snp_vcf (Str :$filein! ,Object :$previous! ){
	my $name = $self->project->name();
	my $dir_out= $self->project->getCallingPipelineDir("unifiedgenotyper");
	my $arg = $self->argument_patient(); 
	my $cmd = "perl $bin_dev/move_vcf.pl  -type=snp -project=$name -vcf_dir=$dir_out -patient=$arg -log=".$self->log_file;
	my $stepname = $name."_move_vcf";
	 my $ppn = 1;
	my ($job,$job_next) =  $self->construct_jobs(stepname=>$stepname,cmd=>$cmd,ppn=>$ppn);
		$self->add_running_steps($stepname);
	  if (defined $previous) {
		$self->prev(job=>$job,prev=>$previous);
	}	
	  return ($job,$filein);
}

1;
