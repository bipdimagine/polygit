package calling;
use strict;
use FindBin qw($Bin);
use Moose;  
use MooseX::Method::Signatures;
use Data::Dumper;
use Data::Printer;
use colored;
use threads;
use Thread::Queue;
use List::MoreUtils qw(part natatime);

my $bin_dev = qq{$Bin/scripts/scripts_pipeline};

 #(Str :$filein! ,Object :$previous! ){
#method calling_merge1  (Int :$fork! ){
	
	my @tmp_files;
	
	sub calling_merge {
		my ($project_name,$chr_name,$patient_names,$fork,$methods,$dir_out,$low_calling) = @_;
	
	# Object :$chr! ,ArrayRef :$patients! ,Int :$fork!
	# Object :$chr!
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $patients =$project->get_list_patients($patient_names);
my $chr = $project->getChromosome($chr_name);
my $reference = $project->genomeFasta();
unless (-e $dir_out){
system("mkdir -p $dir_out ") ;
system("chmod a+rwx $dir_out");
}
my $freebayes_min_alternate = 0.2;
my $unified_min_alternate = 0.1;

#$low_calling =1;
my $low_calling_filter = "" ;
if ($low_calling){
	colored::stabilo('red',"  ----SOMATIC PARAMETER  ----");
	$freebayes_min_alternate = 0.02;
	$unified_min_alternate = 0.02;
	$low_calling_filter = qq{| $bin_dev/correct_tumoral.pl -limit $freebayes_min_alternate };
}
#my $dir_out= $project->getCallingPipelineDir("unifiedgenotyper"); 
my $javac = $project->getSoftware('java');
my $gatk  = $project->getSoftware('gatk');
my $freebayes = $buffer->software("freebayes");
my $bcftools = $buffer->software("bcftools");
my $vcfutil = $buffer->software("vcfutils");
my $samtools  = $buffer->software("samtools");
my $vcffilter =  $buffer->software("vcffilter");
my $vcffirstheader = $buffer->software("vcffirstheader");
	my $vcfstreamsort = $buffer->software("vcfstreamsort");
	my $vcfuniq  = $buffer->software("vcfuniq");
	
my $hvcf;
$hvcf->{freebayes}  = getTmpFile($dir_out,$chr_name,"final.freebayes.vcf");
$hvcf->{samtools}   = getTmpFile($dir_out,$chr_name,"final.samtools.vcf");
 $hvcf->{unifiedgenotyper}  =getTmpFile($dir_out,$chr_name,"final.unifiedgenotyper.vcf");
 $hvcf->{haplotypecaller}  = getTmpFile($dir_out,$chr_name,"final.haplotypecaller.vcf");
 
	foreach my $type (keys %$methods){
		die($type) unless exists  $hvcf->{$type};
		unlink $hvcf->{$type} if -e $hvcf->{$type};
	}



my $span = $chr->getIntSpanCaptureForCalling(300);
return undef if $span->is_empty;

my @t = $span->as_array;

my $nb_bases = scalar(@t);

my $limit = 10000;
my $it = natatime($limit, @t);
my @bed2;

my $bed_array;
my $zz =0;
  while (my @vals = $it->())
  {
  	my $set = Set::IntSpan::Fast::XS->new();
  	$set->add(@vals);
  	my @tt = map{$_ = $chr->ucsc_name."\t".$_} split(";",$set->as_string({ sep => ";", range => "\t" }) );
  	if (scalar(@vals) > $limit/2){
  	push(@{$bed2[$zz]},@tt);
  	$zz ++;
  	}
  	else {
  		push(@{$bed2[$zz]},@tt);
  	}
  #	push(@$bed_array,\@tt );
  }
my $limit2 = int(scalar(@bed2) / ($fork*5)); 
$limit2 = scalar(@bed2) if $limit2 == 0;
my $it2 = natatime($limit2, @bed2);

  while (my @vals = $it2->()){
  		my $at;
  		foreach my $a (@vals){
  				push(@$at,@$a);
  		}
  	
  		push(@$bed_array,$at );
  }
 

my @bed = map{$_ = $chr->ucsc_name."\t".$_} split(";",$span->as_string({ sep => ";", range => "\t" }));


###################
#### #RECAL FILE
#################
colored::stabilo('blue',"  ---- start recalibration $chr_name  ----");
my $bam_files_recal = printReadsFork($patients,$chr_name,\@bed,$fork,$dir_out);
 my $file_list_bam = getTmpFile($dir_out,$chr_name,"list");  
 my $bam_file_string_hc ="";
 
open(LIST,">$file_list_bam");

foreach my $f (@$bam_files_recal){
	print LIST $f."\n";
	$bam_file_string_hc  .= " -I ".$f;
}

close(LIST);
###########
 my $hfiles;
#############
colored::stabilo('green',"  ---- end  recalibration $chr_name  ----");
#my @toto = `cat $file_list_bam`;

foreach my $f (@$bam_files_recal){
     		die $f unless -e $f;
     	}
#warn $bam_file_string_hc;
my $dbsnp = $project->gatk_dbsnp_file();
my $onefile;
	if (scalar ( @$bam_files_recal)==1){
		$onefile = $bam_files_recal->[0];
	}
	
 my $nb = scalar(@$bed_array);
 my $cpt= 0;
 colored::stabilo('blue',"  ---- start  calling :  ".join(",",keys %$methods)." $chr_name  fork : $fork ----");
#my $limit = 10;
my $time_start = time;
my $run = 0;
#my $bed_array;

my $l =0;
my $index =scalar(@$bed_array);

	
my $pm = new Parallel::ForkManager(int($fork));

my $finish = 0;
$pm->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    	
    		$finish ++;
    		my $complete = int(($finish/$index) *100);
    		#colored::stabilo('cyan'," $chr_name : $finish $index $complete  ".$complete%10);
    		if ($complete%10  ==0 && $complete >0){
    			my $elapsed =  int (abs(time -$time_start)/60);	
    		my $still = abs(int(( 1-(100 / $complete)) * $elapsed));
    	
    		colored::stabilo('white',"  ---- end  $finish/$index freebayes   $chr_name   $elapsed mn  remained : $still ----");
    		}
		    		
    		
    }
  );

 foreach my $bed (@$bed_array){
 	foreach my $type (keys %$methods){
 		next if $type eq "unifiedgenotyper";
 		push(@{$hfiles->{$type}->{array_files}},getTmpFile($dir_out,$chr_name,"vcf"));
 	}  
	my $pid = $pm->start and next;
       my $bed_temp = getTmpFile($dir_out,$chr_name,"bed");
          my $region_temp = getTmpFile($dir_out,$chr_name,"region");
 	open (BED,">$bed_temp" ) or die();
 	open (REGION,">$region_temp" ) or die();
 	my ($titi,$m_start,$toto) =  split(" ",$bed->[0]);
 	my $m_end;
 	foreach my $bed_region (@$bed){
 			my ($chr,$start,$end) = split(" ",$bed_region);
 			$end = $start +1 unless $end;
 			$m_end = $end;
 			my $max = abs($end -$start);
 				print BED $bed_region."\n";
 				print REGION "$chr:$start-$end\n";
 		
 	}
 	close BED; 
 	close REGION;
	my $time1 = time;
		
	if (exists $methods->{freebayes}){
		my $file = $hfiles->{freebayes}->{array_files}->[-1];
		my $cmd_free;
			
		if ($onefile){
			 $cmd_free = qq{ cat $region_temp | xargs -I \{\} $freebayes -r \{\} -b $onefile -f   $reference --min-coverage 20  -0  -F $freebayes_min_alternate  --max-complex-gap 0 | $vcffirstheader | $vcfstreamsort -w 1000 | $vcfuniq | $vcffilter -f 'QUAL > 50' $low_calling_filter>$file};
			# $cmd_free = qq{ $samtools view $bam_file_recal -bh -L $bed_temp | $freebayes  --stdin -f  $reference  -F 0.1  | $vcffilter -f 'QUAL > 50' >>$file};
		}
		else {
			$cmd_free = qq{ for i in `cat $region_temp`;do $samtools merge  - -b $file_list_bam -R \$i | $freebayes --stdin -f   $reference --min-coverage 15 -0  -F $freebayes_min_alternate --max-complex-gap 0 -r \$i >>$file;done};
	
		
		 #$cmd_free = qq{$samtools merge  - -b $file_list_bam -R $chr_name:$m_start-$m_end | $samtools view - -bh -L $bed_temp | $freebayes  --stdin -f  $reference -F 0.1 --standard-filters | $vcffilter -f 'QUAL > 20' >$file};
		}
		my $tt = time;
		#colored::stabilo('cyan',"  ---- start  $run freebayes $run - $chr_name  ----");
		
     	system($cmd_free);
     	#colored::stabilo('magenta',"  ---- end  $run freebayes ".int (abs($tt -time)/60)." $run - $chr_name  ----");
     
    
	}
	
   if (exists $methods->{haplotypecaller}){
   	my $file = $hfiles->{haplotypecaller}->{array_files}->[-1];
     	my $cmd = qq{ for i in `cat $region_temp`;do $javac -jar $gatk -T HaplotypeCaller -R $reference --dbsnp $dbsnp $bam_file_string_hc  -stand_call_conf 30 -stand_emit_conf 10 -l off -L  \$i;done > $file };          
     	my $tt = time;
     	system($cmd);
    }
    
     if (exists $methods->{samtools}){
   
	
     	my $file = $hfiles->{samtools}->{array_files}->[-1];
     	
     	
     	my $param_mp = "";#$methods->{samtools}->{param};
     	my $cmd_mpileup;
     	if ($onefile){
			 $cmd_mpileup =qq{$samtools mpileup $onefile -ugf $reference  $param_mp -l $bed_temp -L 8000 | $bcftools call -mv - |  $vcfutil  varFilter -a 2 -d 10  | grep -v "ID=GL" > $file};
		}
		else {
     	 	$cmd_mpileup = qq{$samtools merge - -b $file_list_bam -R $chr_name:$m_start-$m_end | $samtools mpileup -ugf $reference - $param_mp -l $bed_temp | $bcftools call -mv - |  $vcfutil  varFilter -a 2 -d 10  | grep -v "ID=GL" > $file};
		}		
		
		system ("$cmd_mpileup");
    
    } 
#        if (exists $methods->{unifiedgenotyper}){
#        	  $time1 = time;
#        		my $file = $hfiles->{unifiedgenotyper}->{array_files}->[-1];
#        		my $cmd_uni = qq{cat $bed_temp | perl -lane 'print "\$F[0]:\$F[1]-\$F[2]";' | xargs -I \{\} $javac -jar $gatk -T UnifiedGenotyper -L  \{\} --min_indel_fraction_per_sample 0.1 -rf BadCigar -R $reference --dbsnp $dbsnp $bam_file_string_hc  --genotype_likelihoods_model BOTH  -l off | $vcffirstheader | $vcfstreamsort -w 1000 | $vcfuniq  >$file };
#    				#my $cmd_uni = qq{$javac -jar $gatk -T UnifiedGenotyper    --min_indel_fraction_per_sample 0.1 -rf BadCigar -R $reference  -L  $bed_temp --dbsnp $dbsnp $bam_file_string_hc  --genotype_likelihoods_model BOTH  -l off -o $file };
#      			warn $cmd_uni;
#      			system ("$cmd_uni  ");
#      					warn "end time1 $run ". int (abs(time -$time1)/60);
#        }
    
	$pm->finish();

}
$pm->wait_all_children();
	colored::stabilo('yellow',"  ---- end  calling freebayes  ".int (abs(time -$time_start)/60)." $chr_name  ----");
  if (exists $methods->{unifiedgenotyper}){
  	  my $bed_temp = getTmpFile($dir_out,$chr_name,"bed");
 		open (BED,">$bed_temp" ) or die();
 		foreach my $bed_region (@bed){
 			print BED $bed_region."\n";
 	}
 	close BED; 
      my  $time1 = time;
		push(@{$hfiles->{unifiedgenotyper}->{array_files}},getTmpFile($dir_out,$chr_name,"vcf"));
		my $file = $hfiles->{unifiedgenotyper}->{array_files}->[-1];
		my $cmd_uni = qq{$javac -jar $gatk -T UnifiedGenotyper  -dcov 10000  --min_indel_fraction_per_sample $unified_min_alternate -rf BadCigar -R $reference  -nt $fork -L  $bed_temp --dbsnp $dbsnp $bam_file_string_hc  --genotype_likelihoods_model BOTH  -l off -o $file };
 		system ("$cmd_uni ");
      	unlink $bed_temp;
  }


my $t1 = abs(time - $time_start);
 		my $elapse = int($t1/60);
 		
 colored::stabilo('magenta'," @@@  -+-+-+-+-+ end  calling ".join(",",keys %$methods)." $chr_name : $elapse -+-+-+-+ @@@");
foreach my $type (keys %$methods){
	concat_vcf($hfiles->{$type}->{array_files},$hvcf->{$type},$project) ;
	warn $hvcf->{$type};
}
my $param_merge ="";
my $priority ="-priority ";
foreach my $type (sort {$methods->{$a} <=>$methods->{$b} }keys %$methods){
	$param_merge .=" --variant:".$type." ".$hvcf->{$type};
	$priority .="$type,";
}

my $end_ext = "uni";
my $output1  = getTmpFile($dir_out,$chr_name,"vcf");#$dir_out."/".$chr_name.".$end_ext.vcf";


my $cmd_merge =qq{$javac -jar $gatk -T CombineVariants -R $reference  $param_merge  -o $output1 -genotypeMergeOptions PRIORITIZE $priority };
		system($cmd_merge."   " );
 colored::stabilo('green',"  ---- end  calling".join(",",keys %$methods)." $chr_name  ----");
warn "end";
return $output1;


}

sub calling_merge_one_patient {
my ($project_name,$chr_name,$patient_names,$fork,$methods,$dir_out,$somatic) = @_;
	# Object :$chr! ,ArrayRef :$patients! ,Int :$fork!
	# Object :$chr!
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $patient =$project->getPatient($patient_names);

my $chr = $project->getChromosome($chr_name);
my $reference = $project->genomeFasta();
unless (-e $dir_out){
system("mkdir -p $dir_out ") ;
system("chmod a+rwx $dir_out");
}
my $freebayes_min_alternate = 0.1;
my $unified_min_alternate = 0.1;
#$somatic =1;
if ($somatic){
	$freebayes_min_alternate = 0.02;
	$unified_min_alternate = 0.02;
}
#my $dir_out= $project->getCallingPipelineDir("unifiedgenotyper"); 
my $javac = $project->getSoftware('java');
my $gatk  = $project->getSoftware('gatk');
my $freebayes = $buffer->software("freebayes");
my $bcftools = $buffer->software("bcftools");
my $vcfutil = $buffer->software("vcfutils");
my $samtools  = $buffer->software("samtools");
my $vcffilter =  $buffer->software("vcffilter");
my $vcffirstheader = $buffer->software("vcffirstheader");
	my $vcfstreamsort = $buffer->software("vcfstreamsort");
	my $vcfuniq  = $buffer->software("vcfuniq");
	
my $hvcf;
$hvcf->{freebayes}  = getTmpFile($dir_out,$chr_name,"final.freebayes.vcf");
$hvcf->{samtools}   = getTmpFile($dir_out,$chr_name,"final.samtools.vcf");
 $hvcf->{unifiedgenotyper}  =getTmpFile($dir_out,$chr_name,"final.unifiedgenotyper.vcf");
 $hvcf->{haplotypecaller}  = getTmpFile($dir_out,$chr_name,"final.haplotypecaller.vcf");
 
	foreach my $type (keys %$methods){
		die($type) unless exists  $hvcf->{$type};
		unlink $hvcf->{$type} if -e $hvcf->{$type};
	}



my $span = $chr->getIntSpanCaptureForCalling(300);
return undef if $span->is_empty;


my @bed = map{$_ = $chr->ucsc_name."\t".$_} split(";",$span->as_string({ sep => ";", range => "\t" }));

#colored::stabilo('blue',"  ---- start  calling :  ".join(",",keys %$methods)." chromosome $chr_name $patient_names----");
###################
#### #RECAL FILE
#################
my $onefile = printReads($patient,$chr_name,\@bed,$fork,$dir_out);
 my $file_list_bam = getTmpFile($dir_out,$chr_name,"list");  
 
 my $bam_file_string_hc =" -I ".$onefile;
 
die() unless -e $onefile;

###########
#############

 
my $time_start = time;
my $run = 0;
my $l =0;

my $dbsnp = $project->gatk_dbsnp_file();	

my $hfiles;

 	foreach my $type (keys %$methods){
 		next if $type eq "unifiedgenotyper";
 		push(@{$hfiles->{$type}->{array_files}},getTmpFile($dir_out,$chr_name,"vcf"));
 	}  
     my $bed_temp = getTmpFile($dir_out,$chr_name,"bed");
    my $region_temp = getTmpFile($dir_out,$chr_name,"region");
 	open (BED,">$bed_temp" ) or die();
 	open (REGION,">$region_temp" ) or die();
 	my ($titi,$m_start,$toto) =  split(" ",$bed[0]);
 	my $m_end;
 	foreach my $bed_region (@bed){
 			my ($chr,$start,$end) = split(" ",$bed_region);
 			$end = $start +1 unless $end;
 			$m_end = $end;
 			my $max = abs($end -$start);
 				print BED $bed_region."\n";
 				print REGION "$chr:$start-$end\n";
 		
 	}
 	close BED; 
 	close REGION;
	my $time1 = time;
		
	if (exists $methods->{freebayes}){
		my $file = $hfiles->{freebayes}->{array_files}->[-1];
		my $cmd_free;
			 $cmd_free = qq{ cat $region_temp | xargs -I \{\} $freebayes -r \{\} -b $onefile -f   $reference --min-coverage 5 -0  -F $freebayes_min_alternate --max-complex-gap 0 | $vcffirstheader | $vcfstreamsort -w 1000 | $vcfuniq | $vcffilter -f 'QUAL > 50' >>$file};
		my $tt = time;
		warn $cmd_free;
     	system($cmd_free);
    
	}
	
   if (exists $methods->{haplotypecaller}){
   		my $file = $hfiles->{haplotypecaller}->{array_files}->[-1];
     	my $cmd = qq{ for i in `cat $region_temp`;do $javac -jar $gatk -T HaplotypeCaller -R $reference --dbsnp $dbsnp $bam_file_string_hc  -stand_call_conf 30 -stand_emit_conf 10 -l off -L  \$i;done > $file };          
     	my $tt = time;
     	warn "start haplo";
     	system($cmd);
     	warn "end haplo";
    }
    
     if (exists $methods->{samtools}){
   
	
     	my $file = $hfiles->{samtools}->{array_files}->[-1];
     	
     	
     	my $param_mp = "";#$methods->{samtools}->{param};
     	my $cmd_mpileup;
		$cmd_mpileup =qq{$samtools mpileup $onefile -ugf $reference  $param_mp -l $bed_temp -L 8000 | $bcftools call -mv - |  $vcfutil  varFilter -a 2 -d 10  | grep -v "ID=GL" > $file};
		system ("$cmd_mpileup");
    
    } 


  if (exists $methods->{unifiedgenotyper}){
  	  my $bed_temp = getTmpFile($dir_out,$chr_name,"bed");
 		open (BED,">$bed_temp" ) or die();
 		foreach my $bed_region (@bed){
 			print BED $bed_region."\n";
 	}
 	close BED; 
      my  $time1 = time;
		push(@{$hfiles->{unifiedgenotyper}->{array_files}},getTmpFile($dir_out,$chr_name,"vcf"));
		my $file = $hfiles->{unifiedgenotyper}->{array_files}->[-1];
		my $cmd_uni = qq{$javac -jar $gatk -T UnifiedGenotyper  -dcov 10000  --min_indel_fraction_per_sample $unified_min_alternate -rf BadCigar -R $reference  -nt 1 -L  $bed_temp --dbsnp $dbsnp $bam_file_string_hc  --genotype_likelihoods_model BOTH  -l off -o $file };
 		system ("$cmd_uni ");
      	unlink $bed_temp;
  }


my $t1 = abs(time - $time_start);
 		my $elapse = int($t1/60);
 		

foreach my $type (keys %$methods){
	concat_vcf($hfiles->{$type}->{array_files},$hvcf->{$type},$project) ;
}
my $param_merge ="";
my $priority ="-priority ";
foreach my $type (sort {$methods->{$a} <=>$methods->{$b} }keys %$methods){
	$param_merge .=" --variant:".$type." ".$hvcf->{$type};
	$priority .="$type,";
}

my $end_ext = "uni";
my $output1  = getTmpFile($dir_out,$chr_name,"vcf");#$dir_out."/".$chr_name.".$end_ext.vcf";


my $cmd_merge =qq{$javac -jar $gatk -T CombineVariants -R $reference  $param_merge  -o $output1 -genotypeMergeOptions PRIORITIZE $priority -l off};
#warn $cmd_merge;
		system($cmd_merge."   " );
 #colored::stabilo('magenta'," @@@  -+-+-+-+-+ end  calling $patient_names ".join(",",keys %$methods)." chromosome : $chr_name time:$elapse -+-+-+-+ @@@");
return $output1;


}

sub printReads {
	my ($patient,$chr_name,$beds,$fork,$dir_out) = @_;
	
		my @files;
		my $project = $patient->project;
		my $reference = $project->genomeFasta();

#my $dir_out= $project->getCallingPipelineDir("unifiedgenotyper"); 
my $javac = $project->getSoftware('java');
my $gatk  = $project->getSoftware('gatk');
my $samtools  = $project->getSoftware('samtools');
 my $bed_file = getTmpFile($dir_out,$chr_name,"bed");
 open(BED,">$bed_file");map{print BED $_."\n";} @$beds;close BED;
 

		
         my $bam =  $patient->getBamFile();
		my $recal = $patient->getRecalFile();
		die($bam) unless -e $bam;
		my $out = "$dir_out/".$patient->name.".$chr_name.bam";
		 unlink $out if	-e $out;
		if (-e $recal){
		my $cmd = qq{$javac -jar $gatk -T PrintReads   -o $out  -rf BadCigar -R $reference -I $bam -BQSR $recal -L $bed_file  --interval_padding 50 -l off >/dev/null};
		system($cmd);
		}
		else {
	       
         	system("ln -s $bam $out");
         	#warn $bam;
         	#system ("$samtools index $out");
         	unlink $out.".bai" if	-e $out.".bai";
         	system("ln -s $bam.bai $out.bai");
		}

	
	
    return $out;
}



sub printReadsFork {
	my ($patients,$chr_name,$beds,$fork,$dir_out) = @_;
	
		my @files;
		my $project = $patients->[0]->project;
		my $reference = $project->genomeFasta();

#my $dir_out= $project->getCallingPipelineDir("unifiedgenotyper"); 
my $javac = $project->getSoftware('java');
my $gatk  = $project->getSoftware('gatk');
my $samtools  = $project->getSoftware('samtools');
 my $bed_file = getTmpFile($dir_out,$chr_name,"bed");
 open(BED,">$bed_file");map{print BED $_."\n";} @$beds;close BED;
 
my %hlist;
my $tlist;
my $z =0;

foreach my $p (@$patients){
         my $bam =  $p->getBamFile();
		my $recal = $p->getRecalFile();
		die($bam) unless -e $bam;
			my $out = "$dir_out/".$p->name.".$chr_name.bam";
			  push(@files,$out);
			my $t;
			$t->{bam} = $bam;
			$t->{recal} = $recal;
			$t->{out} = $out;
			$t->{bed} = $bed_file;
			$t->{java} = $javac;
			$t->{gatk} = $gatk;
			$t->{reference} = $reference;
			
		  
		  if (-e $out){
			my $od = -M $out;
			my $id = -M $bam;
			my ($size1) =`$samtools idxstats $out | grep $chr_name | cut -f 3`;
			my ($size2) =`$samtools idxstats $bam | grep $chr_name | cut -f 3`;
			chomp($size1);
			
			chomp($size2);
		#	print $od." ".$id." ".$size1. " ".$size2."  $bam \n";
			next if $od < $id && $size1 > ($size2*0.75);
		}
		push(@$tlist,$t);
		$hlist{$p->name}->{bam} = $bam;
		 $hlist{$p->name}->{recal} = $recal;
		  $hlist{$p->name}->{out} = $out;
		  
		  $hlist{$p->name}->{bed} = $bed_file;
		  
		
}


my $fork1 = int ($fork/2);
#my $fork1 = $fork -4;
#$fork1= 2;
$fork1 =6 if $fork1>6;
#$fork1 = 1 if scalar(@$tlist) < 6;
my $pm = new Parallel::ForkManager($fork1);
my $thr;
print "FORK $fork1\n";

	foreach my $t (@$tlist){
		my $pid = $pm->start and next;
		my $f1 = $t->{recal};
		sleep(5) if $f1 eq "none";
		my $bam = $t->{bam};
		my $out = $t->{out};
		my $file_bed = $t->{bed};
		my $javac =  $t->{java};
		my $gatk =  $t->{gatk};
		my $reference =  $t->{reference};

		#next if -e $out;
		if (-e $f1){
		my $cmd = qq{$javac -jar $gatk -T PrintReads   -o $out  -rf BadCigar -R $reference -I $bam -BQSR $f1 -L $file_bed  --interval_padding 50 -l off >/dev/null};
		system($cmd);
		}
		else {
	        unlink $out if	-e $out;
         	system("ln -s $bam $out");
         	#warn $bam;
         	#system ("$samtools index $out");
         	unlink $out.".bai" if	-e $out.".bai";
         	system("ln -s $bam.bai $out.bai");
		}

         $pm->finish();
	}
	
	
	$pm->wait_all_children();
	
    return \@files;
}



sub getTmpFile {
	my ($dir,$chr_name,$ext) = @_;
	die() unless $ext;
	die() unless $chr_name;
	#confess($chr_name) unless $chr_name !~/[0,1,2,3,4,5,6,7,8,9,10,11,12,13,15,15,16,17,18,19,20,21,22,X,Y,M]/;
	$dir .="/$chr_name";
	system ("mkdir -p $dir") unless -d $dir;
	 my $file_cat_tmp =  File::Temp->new( TEMPLATE =>"TMP.XXXXXXX",
                        DIR => $dir,
                        SUFFIX => ".$ext");
  return $file_cat_tmp->filename();
}

sub concat_vcf {
	my ($files,$vcf,$project) = @_;
	my $buffer = $project->buffer;
	my $dir= $project->getCallingPipelineDir("unifiedgenotyper");
	my $vcffirstheader = $buffer->software("vcffirstheader");
	my $vcfstreamsort = $buffer->software("vcfstreamsort");
	my $vcfuniq  = $buffer->software("vcfuniq");

  my $cat_name = getTmpFile($dir,"vcf_cat","vcf");
  for my $list_tmp_file (@$files){
	next unless -e $list_tmp_file;
	system "cat $list_tmp_file>> $cat_name";
	unlink $list_tmp_file;
}

my $join_cmd = "cat $cat_name | $vcffirstheader | $vcfstreamsort -w 1000 | $vcfuniq > $vcf" if -e $cat_name;
system($join_cmd);
unlink $cat_name;
die() unless -e $vcf;	
}


1;

