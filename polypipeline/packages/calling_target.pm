package calling_target;
use strict;
use FindBin qw($Bin);
use Moose;
use MooseX::Method::Signatures;
use Data::Dumper;
use Data::Printer;

#use colored;
use threads;
use Thread::Queue;
use List::MoreUtils qw(part natatime);

my @tmp_files;


sub calling_merge {
	my ( $project_name, $chr_name, $start,   $end,$patient_name, $fork,     $methods, $dir_out,$bin_dev) = @_;

	# Object :$chr! ,ArrayRef :$patients! ,Int :$fork!
	# Object :$chr!
	my $buffer    = GBuffer->new();
	my $project   = $buffer->newProject( -name => $project_name );
	my $patient   = $project->getPatient($patient_name);
	my $vcffilter = $buffer->software("vcffilter");
	my $chr   = $project->getChromosome($chr_name);
	my $reference = $project->genomeFasta();

	unless ( -e $dir_out ) {
		system("mkdir -p $dir_out ");
		system("chmod a+rwx $dir_out");
	}
	my $freebayes_min_alternate = 0.1;
	my $unified_min_alternate   = 0.1;
	my $low_calling_filter      = "";
	my $limit_qual              = " | $vcffilter -f 'QUAL > 50' ";
#	if ($low_calling) {
#		$freebayes_min_alternate = 0.01;
#		$unified_min_alternate   = 0.01;
#		$limit_qual              = "";
#		$low_calling_filter =
#		  qq{| $bin_dev/correct_tumoral.pl -limit $freebayes_min_alternate };
#	}

	#my $dir_out= $project->getCallingPipelineDir("unifiedgenotyper");
	my $javac = $project->getSoftware('java');

	#$javac ="java";
	my $gatk      = $project->getSoftware('gatk');
	my $freebayes = $buffer->software("freebayes");
	my $bcftools  = $buffer->software("bcftools");
	my $vcfutil   = $buffer->software("vcfutils");
	my $samtools  = $buffer->software("samtools");
	my $lofreq  = $buffer->software("lofreq");

	my $vcffirstheader = $buffer->software("vcffirstheader");
	my $vcfstreamsort  = $buffer->software("vcfstreamsort");
	my $vcfuniq        = $buffer->software("vcfuniq");
	my $sub            = "$chr_name.$start.$end";
	my $hvcf;
	$hvcf->{freebayes} =	  getTmpFile( $dir_out, $chr_name, "$sub.final.freebayes.vcf" );
	$hvcf->{samtools} =   getTmpFile( $dir_out, $chr_name, "$sub.final.samtools.vcf" );
	$hvcf->{unifiedgenotyper} =  getTmpFile( $dir_out, $chr_name, "$sub.final.unifiedgenotyper.vcf" );
	$hvcf->{haplotypecaller} = 	  getTmpFile( $dir_out, $chr_name, "$sub.final.haplotypecaller.g.vcf" );
	$hvcf->{haplotypecaller4} = 	  getTmpFile( $dir_out, $chr_name, "$sub.final.haplotypecaller4.g.vcf" );
	$hvcf->{mutect2} =	  getTmpFile( $dir_out, $chr_name, "$sub.final.mutect2.vcf" );
	$hvcf->{lofreq} =	  getTmpFile( $dir_out, $chr_name, "$sub.final.lofreq.vcf" );
	
	foreach my $type ( keys %$methods ) {
		die($type) unless exists $hvcf->{$type};
		unlink $hvcf->{$type} if -e $hvcf->{$type};
	}

###################
#### #RECAL FILE
#################
	my $onefile = $patient->getBamFile();
	die() unless -e $onefile;

	my $bam_file_string_hc = " -I " . $onefile;

	die() unless -e $onefile;

###########
#############

	my $time_start = time;
	my $run        = 0;
	my $l          = 0;

	my $dbsnp = $project->gatk_dbsnp_file();

	my $hfiles;

	my $gatk_region = "-L " . $chr->ucsc_name . ":$start-$end";
	my $free_region = "-r " . $chr->ucsc_name . ":$start-$end";

	my $time1 = time;

	if ( exists $methods->{freebayes} ) {

		#die();

		my $file = $hvcf->{freebayes};

#	my	 $cmd_free = qq{  $freebayes -b $onefile -f   $reference --min-coverage 20 -0  -F $freebayes_min_alternate  $free_region | $vcffirstheader | $vcfstreamsort -w 1000 | $vcfuniq $limit_qual $low_calling_filter>>$file};
		my $cmd_free =qq{  $freebayes -b $onefile -f   $reference --min-coverage 20 -0  -F $freebayes_min_alternate  $free_region | $vcfuniq $limit_qual $low_calling_filter >$file};
		warn $cmd_free;

#	my  $cmd_free = qq{  $freebayes --target $bed_temp -b $onefile -f   $reference --min-coverage 15 -0  -F $freebayes_min_alternate --max-complex-gap 0  | $vcffilter -f 'QUAL > 50' >$file};
#$cmd_free = qq{ for i in `cat $region_temp`;do $samtools view   -b $file_list_bam  \$i | $freebayes --stdin -f   $reference --min-coverage 15 -0  -F $freebayes_min_alternate --max-complex-gap 0  >>$file;done};
#	warn $cmd_free;

		system($cmd_free);

	}

	if ( exists $methods->{haplotypecaller} ) {
		my $file         = $hvcf->{haplotypecaller};
		my $recal        = $patient->getRecalFile();
		my $recal_string = "";
		$recal_string = "--BQSR $recal" if -e $recal;
		my $cmd =
qq{$javac -jar $gatk -T HaplotypeCaller -R $reference  $recal_string $bam_file_string_hc  -o $file --emitRefConfidence GVCF -pairHMM VECTOR_LOGLESS_CACHING $gatk_region --disable_indel_quals -nct 1 -variant_index_type LINEAR -variant_index_parameter 128000  -l off };
		system($cmd);
		return $file;

#die();
#   		my $file = $hfiles->{haplotypecaller}->{array_files}->[-1];
#     	my $cmd = qq{ for i in `cat $region_temp`;do $javac -jar $gatk -T HaplotypeCaller -R $reference --dbsnp $dbsnp $bam_file_string_hc  -stand_call_conf 30 -stand_emit_conf 10 -l off -L  \$i;done > $file };
#     	my $tt = time;
#     	warn "start haplo";
#     	system($cmd);
#     	warn "end haplo";
	}

	if ( exists $methods->{samtools} ) {

		#	die();

		#     	my $file = $hfiles->{samtools}->{array_files}->[-1];
		#
		#
		my $param_mp = "";                  #$methods->{samtools}->{param};
		my $file     = $hvcf->{samtools};
		my $cmd_mpileup;
		$cmd_mpileup =qq{$samtools mpileup $onefile -ugf $reference  $param_mp $free_region -L 8000 2>/dev/null | $bcftools call -mv - |  perl $vcfutil  varFilter -a 5 -d 15  |   $vcfutil  varFilter -a 5 | grep -v "ID=GL" > $file 2>/dev/null};
		warn $cmd_mpileup;
		system("$cmd_mpileup");

	}

	if ( exists $methods->{unifiedgenotyper} ) {

		my $time1        = time;
		my $recal        = $patient->getRecalFile();
		my $recal_string = "";
		$recal_string = "--BQSR $recal" if -e $recal;

		my $file = $hvcf->{unifiedgenotyper};
		my $cmd_uni =
qq{$javac -jar $gatk -T UnifiedGenotyper  $recal_string -dcov 10000   --min_indel_fraction_per_sample $unified_min_alternate -rf BadCigar -R $reference $gatk_region  --dbsnp $dbsnp $bam_file_string_hc  --genotype_likelihoods_model BOTH  -l off -o $file };

# $cmd_uni = qq{$javac -jar $gatk -T HaplotypeCaller $recal_string -R $reference --dbsnp $dbsnp $bam_file_string_hc  -stand_call_conf 30 -stand_emit_conf 10 -l off $gatk_region   -o $file};
#warn $cmd_uni;

		system("$cmd_uni ");

		if ( $? == -1 ) {
			warn "failed to execute: $! \n $cmd_uni \n ";
			die();
		}

		#warn "$cmd_uni";
		#unlink $bed_temp;
	}
	
	if ( exists $methods->{mutect2} ) {
		my $name = $patient->name();

		#	die();

		#     	my $file = $hfiles->{samtools}->{array_files}->[-1];
		#
		#
		my $param_mp = "";                  #$methods->{samtools}->{param};
		my $file     = $hvcf->{mutect2};
		my $cmd_mutect2;
		$cmd_mutect2 =qq{$javac -jar $gatk -T MuTect2  --callable-depth 10 --base-quality-score-threshold 25 --dont-use-soft-clipped-bases -tumor $name -R $reference $gatk_region   $bam_file_string_hc   -o $file -l off  };
		warn $cmd_mutect2;
		#system("$cmd_mutect2");
	}
		my $t1     = abs( time - $time_start );
		my $elapse = int( $t1 / 60 );

	my $param_merge = "";
	my $priority    = "-priority ";
	foreach
	  my $type ( sort { $methods->{$a} <=> $methods->{$b} } keys %$methods )
	{
		next unless -e $hvcf->{$type};
		$param_merge .= " --variant:" . $type . " " . $hvcf->{$type};
		$priority .= "$type,";
	}

	my $end_ext = "uni";
	my $output1 = getTmpFile( $dir_out, $chr_name, "$sub.vcf" )
	  ;    #$dir_out."/".$chr_name.".$end_ext.vcf";

	#$output1 = "$dir_out/$chr_name/"

	my $cmd_merge =
qq{$javac -jar $gatk -T CombineVariants -R $reference  $param_merge  -o $output1 -genotypeMergeOptions PRIORITIZE $priority -l off};

	#warn $cmd_merge;

	system( $cmd_merge. "   " );

#colored::stabilo('magenta'," @@@  -+-+-+-+-+ end  calling $patient_names ".join(",",keys %$methods)." chromosome : $chr_name time:$elapse -+-+-+-+ @@@");
	foreach my $type ( keys %$methods ) {
		unlink $hvcf->{$type} if -e $hvcf->{$type};
		unlink $hvcf->{$type} . ".idx" if -e $hvcf->{$type} . ".idx";
	}
	return $output1;

}

sub calling_merge2 {
	my ( $project_name, $chr_name, $patient_names, $fork, $methods, $dir_out,
		$bin_dev )
	  = @_;

	# Object :$chr! ,ArrayRef :$patients! ,Int :$fork!
	# Object :$chr!
	my $buffer    = GBuffer->new();
	my $project   = $buffer->newProject( -name => $project_name );
	my $patient   = $project->getPatient($patient_names);
	my $vcffilter = $buffer->software("vcffilter");
	my $chr       = $project->getChromosome($chr_name);
	my $reference = $project->genomeFasta();
	unless ( -e $dir_out ) {
		system("mkdir -p $dir_out ");
		system("chmod a+rwx $dir_out");
	}
	my $freebayes_min_alternate = 0.1;
	my $unified_min_alternate   = 0.1;
	my $low_calling_filter      = "";
	my $limit_qual              = " | $vcffilter -f 'QUAL > 50' ";
#	if ($low_calling) {
#		$freebayes_min_alternate = 0.02;
#		$unified_min_alternate   = 0.02;
#
#		$limit_qual = "";
#		$low_calling_filter = 		  qq{| $bin_dev/correct_tumoral.pl -limit $freebayes_min_alternate };
#
#	}

	#my $dir_out= $project->getCallingPipelineDir("unifiedgenotyper");
	my $javac     = $project->getSoftware('java');
	my $gatk      = $project->getSoftware('gatk');
	my $freebayes = $buffer->software("freebayes");
	my $bcftools  = $buffer->software("bcftools");
	my $vcfutil   = $buffer->software("vcfutils");
	my $samtools  = $buffer->software("samtools");
	my $lofreq  = $buffer->software("lofreq");

	my $vcffirstheader = $buffer->software("vcffirstheader");
	my $vcfstreamsort  = $buffer->software("vcfstreamsort");
	my $vcfuniq        = $buffer->software("vcfuniq");

	my $hvcf;
	$hvcf->{freebayes} =
	  getTmpFile( $dir_out, $chr_name, "final.freebayes.vcf" );
	$hvcf->{samtools} = getTmpFile( $dir_out, $chr_name, "final.samtools.vcf" );
	$hvcf->{unifiedgenotyper} =
	  getTmpFile( $dir_out, $chr_name, "final.unifiedgenotyper.vcf" );
	$hvcf->{haplotypecaller} =
	  getTmpFile( $dir_out, $chr_name, "final.haplotypecaller.vcf" );
	$hvcf->{mutect2} =
	  getTmpFile( $dir_out, $chr_name, "final.mutect2.vcf" );
	  $hvcf->{lofreq} =
	  getTmpFile( $dir_out, $chr_name, "final.lofreq.vcf" );
	  

	foreach my $type ( keys %$methods ) {
		die($type) unless exists $hvcf->{$type};
		unlink $hvcf->{$type} if -e $hvcf->{$type};
	}

	my $recal = $patient->getRecalFile();
	die( $recal . " is empty" ) if ( -z $recal );

	my $span = $chr->getIntSpanCaptureForCalling(150);
	return undef if $span->is_empty;

	my @bed =
	  map { $_ = $chr->ucsc_name . "\t" . $_ }
	  split( ";", $span->as_string( { sep => ";", range => "\t" } ) );

###################
#### #RECAL FILE
#################
	my $onefile = $patient->getBamFile();

	#printReads($patient,$chr_name,\@bed,$fork,$dir_out);
	my $file_list_bam = getTmpFile( $dir_out, $chr_name, "list" );

	my $bam_file_string_hc = " -I " . $onefile;

	#die() unless -e $onefile;
	die() if ( -e $onefile && -z $onefile );
###########
#############

	my $time_start = time;
	my $run        = 0;
	my $l          = 0;

	my $dbsnp = $project->gatk_dbsnp_file();

	my $hfiles;

	my $time1 = time;

	if ( exists $methods->{freebayes} ) {

		#die();
		my $region_temp = getTmpFile( $dir_out, $chr_name, "region" );

		my $bed_temp = getTmpFile( $dir_out, $chr_name, "bed" );
		open( REGION, ">$region_temp" ) or die();
		open( BED,    ">$bed_temp" )    or die();
		foreach my $bed_region (@bed) {
			my ( $chr, $start, $end ) = split( " ", $bed_region );
			print BED $bed_region . "\n";
			print REGION "$chr:$start-$end\n";

		}
		close BED;
		close REGION;

		my $file = $hvcf->{freebayes};
		my $cmd_free =
qq{  $freebayes -b $onefile -f   $reference --min-coverage 10 -0  -F $freebayes_min_alternate   | $vcffirstheader | $vcfstreamsort -w 1000 | $vcfuniq $limit_qual $low_calling_filter>>$file};

		warn $cmd_free;
		sleep(5000);

#	my  $cmd_free = qq{  $freebayes --target $bed_temp -b $onefile -f   $reference --min-coverage 5 -0  -F $freebayes_min_alternate --max-complex-gap 0  | $vcffilter -f 'QUAL > 50' >$file};

#$cmd_free = qq{ for i in `cat $region_temp`;do $samtools view   -b $file_list_bam  \$i | $freebayes --stdin -f   $reference --min-coverage 15 -0  -F $freebayes_min_alternate --max-complex-gap 0  >>$file;done};

		my $tt = time;
		warn $cmd_free;
		system($cmd_free);

	}

	if ( exists $methods->{haplotypecaller} ) {
		die();

#   		my $file = $hfiles->{haplotypecaller}->{array_files}->[-1];
#     	my $cmd = qq{ for i in `cat $region_temp`;do $javac -jar $gatk -T HaplotypeCaller -R $reference --dbsnp $dbsnp $bam_file_string_hc  -stand_call_conf 30 -stand_emit_conf 10 -l off -L  \$i;done > $file };
#     	my $tt = time;
#     	warn "start haplo";
#     	system($cmd);
#     	warn "end haplo";
	}

	if ( exists $methods->{samtools} ) {
		die();

#     	my $file = $hfiles->{samtools}->{array_files}->[-1];
#
#
#     	my $param_mp = "";#$methods->{samtools}->{param};
#     	my $cmd_mpileup;
#		$cmd_mpileup =qq{$samtools mpileup $onefile -ugf $reference  $param_mp -l $bed_temp -L 8000 | $bcftools call -mv - |  $vcfutil  varFilter -a 2 -d 10  | grep -v "ID=GL" > $file};
#		system ("$cmd_mpileup");

	}

	if ( exists $methods->{unifiedgenotyper} ) {
		my $bed_temp = getTmpFile( $dir_out, $chr_name, "bed" );
		open( BED, ">$bed_temp" ) or die();
		foreach my $bed_region (@bed) {
			print BED $bed_region . "\n";
		}
		close BED;
		my $time1        = time;
		my $recal        = $patient->getRecalFile();
		my $recal_string = "";
		$recal_string = "--BQSR $recal" if -e $recal;

		my $file = $hvcf->{unifiedgenotyper};
		my $cmd_uni =
qq{$javac -jar $gatk -T UnifiedGenotyper  -dcov 10000 $recal_string  --min_indel_fraction_per_sample $unified_min_alternate -rf BadCigar -R $reference  -nt $fork  --dbsnp $dbsnp $bam_file_string_hc  --genotype_likelihoods_model BOTH  -l off -o $file };
		$cmd_uni =
qq{$javac -jar $gatk -T HaplotypeCaller -R $reference --dbsnp $dbsnp $bam_file_string_hc  -stand_call_conf 30 -stand_emit_conf 10 -l off  -nct $fork  -o $file};
		warn $cmd_uni;
		system("$cmd_uni ");

		#unlink $bed_temp;
	}

	my $t1     = abs( time - $time_start );
	my $elapse = int( $t1 / 60 );

	my $param_merge = "";
	my $priority    = "-priority ";
	foreach
	  my $type ( sort { $methods->{$a} <=> $methods->{$b} } keys %$methods )
	{
		$param_merge .= " --variant:" . $type . " " . $hvcf->{$type};
		$priority .= "$type,";
	}

	my $end_ext = "uni";
	my $output1 = getTmpFile( $dir_out, $chr_name, "vcf" )
	  ;    #$dir_out."/".$chr_name.".$end_ext.vcf";

	my $cmd_merge =
qq{$javac -jar $gatk -T CombineVariants -R $reference  $param_merge  -o $output1 -genotypeMergeOptions PRIORITIZE $priority -l off};

	#warn $cmd_merge;
	system( $cmd_merge. "   " );

#colored::stabilo('magenta'," @@@  -+-+-+-+-+ end  calling $patient_names ".join(",",keys %$methods)." chromosome : $chr_name time:$elapse -+-+-+-+ @@@");
	return $output1;

}

sub printReads {
	my ( $patient, $chr_name, $beds, $fork, $dir_out ) = @_;

	my @files;
	my $project   = $patient->project;
	my $reference = $project->genomeFasta();

	#my $dir_out= $project->getCallingPipelineDir("unifiedgenotyper");
	my $javac    = $project->getSoftware('java');
	my $gatk     = $project->getSoftware('gatk');
	my $samtools = $project->getSoftware('samtools');
	my $bed_file = getTmpFile( $dir_out, $chr_name, "bed" );
	open( BED, ">$bed_file" );
	map { print BED $_ . "\n"; } @$beds;
	close BED;

	my $bam   = $patient->getBamFile();
	my $recal = $patient->getRecalFile();

	#		die($bam) unless -e $bam;
	die($bam) if -z $bam;
	my $out = "$dir_out/" . $patient->name . ".$chr_name.bam";
	unlink $out if -e $out;

	#		if (-e $recal){
	die( $recal . " is absent" ) unless ( -e $recal );
	die( $recal . " is empty" ) if ( -z $recal );
	unless ( -z $recal ) {

		my $cmd =
qq{$javac -Xmx2048m -Xms256m -jar $gatk -T PrintReads   -o $out  -rf BadCigar -R $reference -I $bam -BQSR $recal -L $bed_file  --interval_padding 50 -l off >/dev/null};
		system($cmd);
	}
	else {
		die( $recal . " is empty" ) if ( -z $recal );
		system("ln -s $bam $out");

		#warn $bam;
		#system ("$samtools index $out");
		unlink $out . ".bai" if -e $out . ".bai";
		system("ln -s $bam.bai $out.bai");
	}

	return $out;
}


sub getTmpFile {
	my ( $dir, $chr_name, $ext ) = @_;
	die() unless $ext;
	die() unless $chr_name;

#confess($chr_name) unless $chr_name !~/[0,1,2,3,4,5,6,7,8,9,10,11,12,13,15,15,16,17,18,19,20,21,22,X,Y,M]/;
	system("mkdir  $dir && chmod a+rwx $dir") unless -d $dir;
	$dir .= "/$chr_name";
	system("mkdir -p $dir && chmod  a+rwx $dir") unless -d $dir;

	my $file_cat_tmp = File::Temp->new(
		TEMPLATE => "TMP.XXXXXXX",
		DIR      => $dir,
		SUFFIX   => ".$ext"
	);
	return $file_cat_tmp->filename();
}

sub concat_vcf {
	my ( $files, $vcf, $project ) = @_;
	my $buffer         = $project->buffer;
	my $dir            = $project->getCallingPipelineDir("unifiedgenotyper");
	my $vcffirstheader = $buffer->software("vcffirstheader");
	my $vcfstreamsort  = $buffer->software("vcfstreamsort");
	my $vcfuniq        = $buffer->software("vcfuniq");
	$project->disconnect();
	my $cat_name = getTmpFile( $dir, "vcf_cat", "vcf" );
	unlink $cat_name if -e $cat_name;
	for my $list_tmp_file (@$files) {
		next unless -e $list_tmp_file;
		system "cat $list_tmp_file>> $cat_name";

		#unlink $list_tmp_file;
	}

	my $vcffilter = $buffer->software("vcffilter");
	my $join_cmd =
"cat $cat_name | $vcffirstheader | $vcfstreamsort -w 1000 | $vcfuniq | $vcffilter -f  'DP > 10'  > $vcf"
	  if -e $cat_name;

	system($join_cmd);

	unlink $cat_name;
	die() unless -e $vcf;
}


my $filter_bcftools                = "(AD[0:0]+AD[0:1])>=10 && AD[0:1]/(AD[0:0]+AD[0:1])>=0.07";

sub haplotypecaller {
	my ( $project, $patient, $intspans, $fork, $verbose ) = @_;
	my $buffer    = $project->buffer;
	my $dir_out   = $project->getCallingPipelineDir("haplotypecaller");
	my $reference = $project->genomeFasta();

	my $javac = $buffer->getSoftware('java');

	#$javac ="java";
	my $gatk      = $buffer->getSoftware('gatk');
	my $freebayes = $buffer->software("freebayes");
	my $bcftools  = $buffer->software("bcftools");
	my $vcfutil   = $buffer->software("vcfutils");
	my $samtools  = $buffer->software("samtools");
	my $out = calling_target::getTmpFile( $dir_out, $project->name, "hc.vcf" );
	my $out1 =
	  calling_target::getTmpFile( $dir_out, $project->name, "hc.1.vcf" );
	my $low_calling_filter = "";
	my $beds = return_all_beds( $project, $intspans );

	if ( scalar(@$beds) == 0 ) {
		print_empty_vcf( $out, $patient );
		return $out;
	}
	my $bed1 = calling_target::getTmpFile( $dir_out, $project->name, "bed" );
	open( BED, ">$bed1" );
	print BED join( "\n", @$beds );
	print BED "\n";
	close BED;
	my $unified_min_alternate = 0.07;

#	if ($low_calling) {
#		$unified_min_alternate = 0.02;
#		$filter = "(AD[0]+AD[1])>=10";    # && AD[1]/(AD[0]+AD[1])>0.1
#	}

	my $recal              = $patient->getRecalFile();
	my $recal_string       = "";
	my $bam_file_string_hc = " -I " . $patient->getBamFile();
	my $gatk_region        = "-L $bed1";
	my $gatk4  = $project->getSoftware('gatk4');
	my $cmd_uni =
qq{$gatk4  HaplotypeCaller   $recal_string   -R $reference $gatk_region   $bam_file_string_hc   -o $out1 -l off };
	system( $cmd_uni   . qq{&& $bcftools  filter -i " $filter_bcftools " $out1 > $out && rm $out1} );
	return $out;

}

sub unifiedGenotyper_old {
	my ( $project, $patient, $intspans, $fork, $verbose ) = @_;
	my $buffer    = $project->buffer;
	my $dir_out   = $project->getCallingPipelineDir("unifiedgenotyper");
	my $reference = $project->genomeFasta();
	my $name = $patient->name;
	my $javac = $buffer->getSoftware('java');

	#$javac ="java";
	my $gatk      = $buffer->getSoftware('gatk');
	my $freebayes = $buffer->software("freebayes");
	my $bcftools  = $buffer->software("bcftools");
	my $vcfutil   = $buffer->software("vcfutils");
	my $samtools  = $buffer->software("samtools");
	my $out = calling_target::getTmpFile( $dir_out, $project->name, "uni.vcf" );
	my $out1 = 	  calling_target::getTmpFile( $dir_out, $project->name, "uni.1.vcf" );
	my $low_calling_filter = "";
	my $beds = return_all_beds( $project, $intspans );

	if ( scalar(@$beds) == 0 ) {
		print_empty_vcf( $out, $patient );
		return $out;
	}
	my $bed1 = calling_target::getTmpFile( $dir_out, $project->name, "bed" );
	open( BED, ">$bed1" );
	print BED join( "\n", @$beds );
	print BED "\n";
	close BED;
	my $unified_min_alternate = 0.1;

#	if ($low_calling) {
#		$unified_min_alternate = 0.02;
#		$filter = "(AD[0]+AD[1])>=10";    # && AD[1]/(AD[0]+AD[1])>0.1
#	}

	my $recal              = $patient->getRecalFile();
	my $recal_string       = "";
	my $bam_file_string_hc = " -I " . $patient->getBamFile();
	my $gatk_region        = "-L $bed1";
	my $tmp_dir = "";
	$tmp_dir = "-Djava.io.tmpdir=/data-beegfs/tmp/" if -e "/data-beegfs/tmp/";
	$gatk = "/software/distrib/GATK/3.5/GenomeAnalysisTK.jar";
	my $cmd_uni = qq{$javac $tmp_dir -jar $gatk  -T UnifiedGenotyper  -nt $fork $recal_string    --min_indel_fraction_per_sample $unified_min_alternate -rf BadCigar -R $reference $gatk_region   $bam_file_string_hc  --genotype_likelihoods_model BOTH   -o $out1 };
	warn $cmd_uni;
	system( $cmd_uni
		  . qq{&& $bcftools  filter -i " $filter_bcftools " $out1 > $out && rm $out1} );
	warn $cmd_uni
	  . qq{&& $bcftools  filter -i " $filter_bcftools " $out1 > $out && rm $out1};
	return $out;

}

sub cp {
	my ( $project, $patient, $intspans, $fork, $verbose ) = @_;
	my $buffer    = $project->buffer;
	my $dir_out   = $project->getCallingPipelineDir("unifiedgenotyper_cp");
	my $reference = $project->genomeFasta();

	my $javac = $buffer->getSoftware('java');

	#$javac ="java";
	my $gatk      = $buffer->getSoftware('gatk');
	my $freebayes = $buffer->software("freebayes");
	my $bcftools  = $buffer->software("bcftools");
	my $vcfutil   = $buffer->software("vcfutils");
	my $samtools  = $buffer->software("samtools");
	my $bgzip     = $buffer->software("bgzip");
	my $tabix     = $buffer->software("tabix");

	my $out1 =
	  calling_target::getTmpFile( $dir_out, $project->name, ".uni1.vcf" );
	my $out =
	  calling_target::getTmpFile( $dir_out, $project->name, ".uni.vcf" );
	my $low_calling_filter = "";
	my $beds = return_all_beds( $project, $intspans );
	if ( scalar(@$beds) == 0 ) {
		print_empty_vcf( $out, $patient );
		return $out;
	}
	my $bed1 = calling_target::getTmpFile( $dir_out, $project->name, "bed" );
	open( BED, ">$bed1" );
	print BED join( "\n", @$beds );
	print BED "\n";
	close BED;
	my $unified_min_alternate = 0.1;
#	if ($low_calling) {
#		$unified_min_alternate = 0.02;
#	}

	my $vcf_uni = $patient->getVariationsFile("unifiedgenotyper");

	my $recal              = $patient->getRecalFile();
	my $recal_string       = "";
	my $bam_file_string_hc = " -I " . $patient->getBamFile();
	my $gatk_region        = "-L $bed1";
	my $cmd_uni =
qq{$javac -jar $gatk -T UnifiedGenotyper  -nt $fork $recal_string -dcov 10000   --min_indel_fraction_per_sample $unified_min_alternate -rf BadCigar -R $reference $gatk_region   $bam_file_string_hc  --genotype_likelihoods_model BOTH   -o $out1 -l off };
	$cmd_uni .=
qq{&& $bgzip $out1 && $tabix -p vcf $out1.gz && $bcftools  isec -C $out1.gz $vcf_uni -e "MIN(DP)<20" -O v -o $out -w 1 };
	system($cmd_uni);
	return $out;

}

sub print_empty_vcf {
	my ( $file, $patient ) = @_;
	confess() unless $patient;
	open( OUT, ">$file" );
	print OUT
"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
	  . $patient->name . "\n";
	close OUT;
}

sub freebayes {
	my ( $project, $patient, $intspans, $fork, $verbose ) = @_;

	my $buffer      = $project->buffer;
	my $freebayes   = $buffer->software("freebayes");
	my $pm_samtools = new Parallel::ForkManager($fork);
	my $dir_out     = $project->getCallingPipelineDir("freebayes");
	my $vcf_final =
	  calling_target::getTmpFile( $dir_out, "freebayes", "final.vcf" );
	my @res_samtools;
	my @beds;
	my $c = 0;

	my $reference               = $project->genomeFasta();
	my $onefile                 = $patient->getBamFile();
	my $samtools                = $buffer->software("samtools");
	my $freebayes_min_alternate = 0.1;
	my $correct_calling         = "";

		$correct_calling = qq{|  $Bin/correct_calling_freebayes.pl -bam=$onefile -samtools=$samtools };

	my $beds = return_all_beds( $project, $intspans );
	if ( scalar(@$beds) == 0 ) {
		print_empty_vcf( $vcf_final, $patient );
		return $vcf_final;
	}

#	my $hotspots =[] ;#= $patient->getCapture()->hotspot;
#	my $chr_intspan_hotspots;
#	foreach my $h (@$hotspots){
#		my $chr = $project->getChromosome($h->{chromosome});
#		$chr_intspan_hotspots->{$chr->name} = new Set::IntSpan::Fast::XS->new() unless exists $chr_intspan_hotspots->{$chr->name};
#		$chr_intspan_hotspots->{$chr->name}->add_range($h->{start},$h->{end});
#	}
#	my $pr = String::ProgressBar->new( max => scalar(@$beds) );
	$pm_samtools->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data ) =
			  @_;

			#  	  $pr->update($c++);

			#	$pr->write() if $verbose;
			push( @res_samtools, $data ) if $data->{file};

		}
	);
	my $time     = time;
	my $nb       = 0;
	my $bcftools = $buffer->software("bcftools");

	my $it = natatime $fork , @$beds;
	while ( my @vals = $it->() ) {

		#foreach my $bed (@$beds){

		$nb++;
		my $pid = $pm_samtools->start and next;
		my $bed = calling_target::getTmpFile( $dir_out, "samtools", "$nb.bed" );
		open( BED, ">$bed" );
		foreach my $v (@vals) {
			print BED $v . "\n";
		}
		close BED;

		$nb++;
		my $region;

	 # ($region->{chr},$region->{start},$region->{end}) = split(" ",$bed);
	 #my $chr = $project->getChromosome($region->{chr});
		my $out =
		  calling_target::getTmpFile( $dir_out, "freebayes", "$nb.vcf" );
		my $cmd_free = qq{  $freebayes -b $onefile -f   $reference --min-coverage 20 -0  -F $freebayes_min_alternate  -t $bed $correct_calling >$out 2>/dev/null};


		system($cmd_free);
		unlink $bed;
		my %h;
		$h{file} = $out;
		$h{nb}   = $nb;
		$pm_samtools->finish( 0, \%h );

	}

	$pm_samtools->wait_all_children();

	warn "samtools time :" . abs( time - $time );
	my @all_vcfs =
	  map { $_->{file} } sort { $a->{nb} <=> $b->{nb} } @res_samtools;

	calling_target::concat_vcf( \@all_vcfs, $vcf_final, $project );
	foreach my $f (@all_vcfs) {
		unlink $f;
	}

	warn "freebayes " . abs( $time - time );
	return $vcf_final;
}
sub deepvariant {
	my ( $project, $patient, $intspans, $fork, $verbose ) = @_;
	my $name = $patient->name();
	my $buffer      = $project->buffer;
	my $gatk   = $buffer->software("gatk");
	my $javac   = $buffer->software("java");
	my $dir_out     = $project->getCallingPipelineDir("deepvariant");
	my $vcf_out = calling_target::getTmpFile( $dir_out, "$name", "deep.1.vcf" );
	my $vcf_final = calling_target::getTmpFile( $dir_out, "$name", "deep.vcf" );
	my $beds = return_all_beds( $project, $intspans );
	my $capture = calling_target::getTmpFile( $dir_out, "$name", ".bed" );
	my $ref               = $project->genomeFasta();
	my $bam                 = $patient->getBamFile();
	open( BED, ">$capture" );
	print BED join("\n",@$beds);
	close BED;
	my $cmd = qq{singularity run  --bind /data-isilon:/data-isilon --bind /data-beegfs:/data-beegfs  /software/distrib/deepvariant/deepvariant.sif /opt/deepvariant/bin/run_deepvariant  --num_shards=$fork --model_type=WES --ref=$ref --reads=$bam --regions=$capture --output_vcf=$vcf_out};
	system($cmd);
	my $bcftools = $buffer->software("bcftools");
	my $tabix = $buffer->software("tabix");
	system(qq{$bcftools filter -e 'FILTER="RefCall" || (FORMAT/DP[0]<5) || (FORMAT/DP[0]<10 && QUAL<15)  $vcf_out -o $vcf_final });
	return $vcf_final;
}
sub unifiedGenotyper {
	my ( $project, $patient, $intspans, $fork, $verbose ) = @_;
	my $name = $patient->name();
	my $buffer      = $project->buffer;
	my $gatk   = $buffer->software("gatk");
	my $javac   = $buffer->software("java");
	my $pm_samtools = new Parallel::ForkManager($fork);
	my $dir_out     = $project->getCallingPipelineDir("unifiedgenotyper");
	my $vcf_final = calling_target::getTmpFile( $dir_out, "$name", "uni.vcf" );
	my @res_samtools;
	my @beds;
	my $c = 0;

	my $reference               = $project->genomeFasta();
	my $onefile                 = $patient->getBamFile();
	my $samtools                = $buffer->software("samtools");
	
	my $beds = return_all_beds( $project, $intspans );
	if ( scalar(@$beds) == 0 ) {
		print_empty_vcf( $vcf_final, $patient );
		return $vcf_final;
	}


	$pm_samtools->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data ) =  @_;
			push( @res_samtools, $data ) if $data->{file};
		}
	);
	my $unified_min_alternate = 0.1;
	my $time     = time;
	my $nb       = 0;
	my $bcftools = $buffer->software("bcftools");

	my $it = natatime $fork , @$beds;
	while ( my @vals = $it->() ) {
		$nb++;
		my $out =
		  calling_target::getTmpFile( $dir_out, "$name", "$nb.vcf" );
		  
		my $pid = $pm_samtools->start and next;
		my $bed = calling_target::getTmpFile( $dir_out, "$name", "$nb.bed" );
		open( BED, ">$bed" );
		foreach my $v (@vals) {
			print BED $v . "\n";
		}
		close BED;
		$nb++;
		my $region;
		
		my $out1 = calling_target::getTmpFile( $dir_out, "$name", "$nb.tmp.vcf" );  
#		my $cmd_free = qq{  $freebayes -b $onefile -f   $reference --min-coverage 20 -0  -F $freebayes_min_alternate  -t $bed $correct_calling >$out 2>/dev/null};
		$javac = "/opt/jdk1.8.0_121/bin/java";
		my $cmd_uni = qq{$javac  -jar $gatk  -T UnifiedGenotyper  --min_indel_fraction_per_sample $unified_min_alternate -rf BadCigar -R $reference -L $bed   -I $onefile  --genotype_likelihoods_model BOTH   -o $out1 };
		warn $cmd_uni;
	system( $cmd_uni
		  . qq{&& $bcftools  filter -i " $filter_bcftools " $out1 > $out && rm $out1} );
		unlink $bed;
		my %h;
		$h{file} = $out;
		$h{nb}   = $nb;
		$pm_samtools->finish( 0, \%h );

	}

	$pm_samtools->wait_all_children();

	warn "samtools time :" . abs( time - $time );
	my @all_vcfs =
	  map { $_->{file} } sort { $a->{nb} <=> $b->{nb} } @res_samtools;

	calling_target::concat_vcf( \@all_vcfs, $vcf_final, $project );
	foreach my $f (@all_vcfs) {
		unlink $f;
	}

	warn "unifiedgenotyper " . abs( $time - time );
	return $vcf_final;
}


sub p1_freebayes {
	my ( $project, $patient, $intspans, $fork, $verbose ) = @_;
	my $buffer      = $project->buffer;
	my $freebayes   = $buffer->software("freebayes");
	my $pm_samtools = new Parallel::ForkManager($fork);
	my $dir_out     = $project->getCallingPipelineDir("freebayes");
	my $vcf_final = calling_target::getTmpFile( $dir_out, "freebayes", "final.vcf" );
	my @res_samtools;
	my @beds;
	my $c   = 0;
	my $reference    = $project->genomeFasta();
	my $onefile    = $patient->getBamFile();
	my $samtools    = $buffer->software("samtools");
	my $freebayes_min_alternate = 0.01;
	my $correct_calling         = "";
	$correct_calling = 	  qq{| $Bin/correct_tumoral.pl -limit $freebayes_min_alternate };
	my $beds = return_all_beds( $project, $intspans );

	if ( scalar(@$beds) == 0 ) {
		print_empty_vcf( $vcf_final, $patient );
		return $vcf_final;
	}
	$pm_samtools->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data ) =
			  @_;
			push( @res_samtools, $data ) if $data->{file};
		}
	);
	my $time     = time;
	my $nb       = 0;
	my $bcftools = $buffer->software("bcftools");

	$fork = scalar(@$beds) if scalar(@$beds) < $fork;
	my $it = natatime $fork , @$beds;
	
	$project->disconnect();
	while ( my @vals = $it->() ) {
		$nb++;
		my $pid = $pm_samtools->start and next;
		my $bed = calling_target::getTmpFile( $dir_out, "samtools", "$nb.bed" );
		open( BED, ">$bed" );
		foreach my $v (@vals) {
			print BED $v . "\n";
		}
		close BED;

		$nb++;
		my $region;

		my $out =   calling_target::getTmpFile( $dir_out, "freebayes", "$nb.vcf" );
		my $cmd_free = qq{  $freebayes -b $onefile -f   $reference  -0  -F $freebayes_min_alternate  -t $bed $correct_calling >$out 2>/dev/null};
		#my $cmd_free = qq{  $freebayes -b $onefile -f   $reference --min-coverage 20 -0  -F $freebayes_min_alternate  -t $bed $correct_calling >$out 2>/dev/null};

		system($cmd_free);
		unlink $bed;
		my %h;
		$h{file} = $out;
		$h{nb}   = $nb;
		$pm_samtools->finish( 0, \%h );

		#}
	}
	sleep(3);
	$pm_samtools->wait_all_children();

	$project->disconnect();
	
	warn "samtools time :" . abs( time - $time );
	my @all_vcfs =
	  map { $_->{file} } sort { $a->{nb} <=> $b->{nb} } @res_samtools;

	calling_target::concat_vcf( \@all_vcfs, $vcf_final, $project );
	foreach my $f (@all_vcfs) {
		unlink $f;
	}

	warn "1p freebayes " . abs( $time - time );
	return $vcf_final;
}

sub eif6_freebayes {
	my ( $project, $patient, $intspans, $fork, $verbose ) = @_;
	my $buffer      = $project->buffer;
	my $freebayes   = $buffer->software("freebayes");
	my $pm_samtools = new Parallel::ForkManager($fork);
	my $dir_out     = $project->getCallingPipelineDir("freebayes");
	my $vcf_final =
	  calling_target::getTmpFile( $dir_out, "freebayes", "final.vcf" );
	my @res_samtools;
	my @beds;
	my $c   = 0;
	my $reference    = $project->genomeFasta();
	my $onefile    = $patient->getBamFile();
	my $samtools    = $buffer->software("samtools");
	my $freebayes_min_alternate = 0.0005;
	my $correct_calling         = "";
	$correct_calling = 	  qq{| $Bin/correct_tumoral.pl -limit $freebayes_min_alternate };
	my $beds = return_all_beds( $project, $intspans );

	if ( scalar(@$beds) == 0 ) {
		print_empty_vcf( $vcf_final, $patient );
		return $vcf_final;
	}
	$pm_samtools->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data ) =
			  @_;
			push( @res_samtools, $data ) if $data->{file};
		}
	);
	my $time     = time;
	my $nb       = 0;
	my $bcftools = $buffer->software("bcftools");

	my $it = natatime $fork , @$beds;
	while ( my @vals = $it->() ) {
		$nb++;
		my $pid = $pm_samtools->start and next;
		my $bed = calling_target::getTmpFile( $dir_out, "samtools", "$nb.bed" );
		open( BED, ">$bed" );
		foreach my $v (@vals) {
			print BED $v . "\n";
		}
		close BED;

		$nb++;
		my $region;

		my $out =   calling_target::getTmpFile( $dir_out, "freebayes", "$nb.vcf" );
		my $cmd_free = qq{  $freebayes -b $onefile -f   $reference --min-coverage 20 -0  -F $freebayes_min_alternate  -t $bed $correct_calling >$out 2>/dev/null};
		system($cmd_free);
		unlink $bed;
		my %h;
		$h{file} = $out;
		$h{nb}   = $nb;
		$pm_samtools->finish( 0, \%h );
	}


	$pm_samtools->wait_all_children();

	warn "samtools time :" . abs( time - $time );
	my @all_vcfs =
	  map { $_->{file} } sort { $a->{nb} <=> $b->{nb} } @res_samtools;

	calling_target::concat_vcf( \@all_vcfs, $vcf_final, $project );
	foreach my $f (@all_vcfs) {
		unlink $f;
	}

	warn "freebayes eif6" . abs( $time - time );
	return $vcf_final;
}

sub freebayes_eif6 {
	my ($project,$patient,$intspans,$low_calling,$fork,$verbose) = @_;
	my $buffer = $project->buffer;
	my $freebayes  = $buffer->software("freebayes");
	$low_calling =1;
	my $pm_samtools = new Parallel::ForkManager($fork);
		my $dir_out= $project->getCallingPipelineDir("freebayes_lc");
	my $vcf_final =  calling_target::getTmpFile($dir_out,"freebayes_lc","final.vcf");
	my @res_samtools;
	my @beds ;
	my $c =0;

	my $reference = $project->genomeFasta();
	
	 my $freebayes_min_alternate = 0.1;
	 my $low_calling_filter="";
	if ($low_calling){
		$freebayes_min_alternate = 0.0005;
		warn $freebayes_min_alternate ;
	#$unified_min_alternate = 0.02;
	#$limit_qual = "";
	$low_calling_filter = qq{| $Bin/correct_tumoral.pl -limit $freebayes_min_alternate };
	}

	my $beds = return_all_beds($project,$intspans);
	if (scalar(@$beds)==0){
		print_empty_vcf($vcf_final,$patient);
		return $vcf_final;
	}
	
	my $pr = String::ProgressBar->new( max => scalar(@$beds) );
 $pm_samtools->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    	  $pr->update($c++);
	
 		$pr->write() if $verbose;
    		push (@res_samtools,$data) if  $data->{file};
		    		
    		
    }
  );
 my $time = time;
 my $nb =0;
 my $onefile = $patient->getBamFile();
 my $bcftools = $buffer->software("bcftools");

	
 	foreach my $bed (@$beds){
 		
	$nb ++;
	 my $pid = $pm_samtools->start and next;
	 my $region;
	 ($region->{chr},$region->{start},$region->{end}) = split(" ",$bed);
	my $chr = $project->getChromosome($region->{chr});
	my $free_region = $chr->fasta_name.":".$region->{start}."-".$region->{end};
	my $name_region = $chr->fasta_name.".".$region->{start}.".".$region->{end};
	my $out =  calling_target::getTmpFile($dir_out,"freebayes","$name_region.vcf");
	
	my	$cmd_free = qq{  $freebayes -b $onefile -f   $reference --min-coverage 20 -0  -F $freebayes_min_alternate  -r $free_region $low_calling_filter >$out};
	
	#system ($cmd_mpileup);
	system ($cmd_free);
	my %h;
	$h{file} = $out;
	$h{nb} = $nb;
	$pm_samtools->finish(0,\%h);
	
 	#}
 	}
   

	
	
$pm_samtools->wait_all_children();

warn "samtools time :".abs(time - $time);
my @all_vcfs = map {$_->{file}} sort{$a->{nb}<=> $b->{nb}} @res_samtools;


calling_target::concat_vcf(\@all_vcfs,$vcf_final,$project) ;
foreach my $f (@all_vcfs){
	unlink $f;
}

warn "freebayes ".abs($time-time);
warn 	$vcf_final;
sleep(50);
return $vcf_final;
}


sub p02_freebayes {
	my ( $project, $patient, $intspans, $fork, $verbose ) = @_;
	my $buffer      = $project->buffer;
	my $freebayes   = $buffer->software("freebayes");
	my $pm_samtools = new Parallel::ForkManager($fork);
	my $dir_out     = $project->getCallingPipelineDir("freebayes");
	my $vcf_final = calling_target::getTmpFile( $dir_out, "freebayes", "final.vcf" );
	my @res_samtools;
	my @beds;
	my $c   = 0;
	my $reference    = $project->genomeFasta();
	my $onefile    = $patient->getBamFile();
	my $samtools    = $buffer->software("samtools");
	my $freebayes_min_alternate = 0.005;
	my $correct_calling         = "";
	$correct_calling = 	  qq{| $Bin/correct_tumoral.pl -limit $freebayes_min_alternate };
	my $beds = return_all_beds( $project, $intspans );

	if ( scalar(@$beds) == 0 ) {
		print_empty_vcf( $vcf_final, $patient );
		return $vcf_final;
	}
	$pm_samtools->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data ) =
			  @_;
			push( @res_samtools, $data ) if $data->{file};
		}
	);
	my $time     = time;
	my $nb       = 0;
	my $bcftools = $buffer->software("bcftools");

	my $it = natatime $fork , @$beds;
	while ( my @vals = $it->() ) {
		$nb++;
		my $pid = $pm_samtools->start and next;
		my $bed = calling_target::getTmpFile( $dir_out, "samtools", "$nb.bed" );
		open( BED, ">$bed" );
		foreach my $v (@vals) {
			print BED $v . "\n";
		}
		close BED;

		$nb++;
		my $region;

		my $out =   calling_target::getTmpFile( $dir_out, "freebayes", "$nb.vcf" );
		my $cmd_free = qq{  $freebayes -b $onefile -f   $reference  -0  -F $freebayes_min_alternate  -t $bed $correct_calling >$out 2>/dev/null};
		#my $cmd_free = qq{  $freebayes -b $onefile -f   $reference --min-coverage 20 -0  -F $freebayes_min_alternate  -t $bed $correct_calling >$out 2>/dev/null};

		system($cmd_free);
#		unlink $bed;
		my %h;
		$h{file} = $out;
		$h{nb}   = $nb;
		$pm_samtools->finish( 0, \%h );

		#}
	}

	$pm_samtools->wait_all_children();

	warn "samtools time :" . abs( time - $time );
	my @all_vcfs =
	  map { $_->{file} } sort { $a->{nb} <=> $b->{nb} } @res_samtools;

	calling_target::concat_vcf( \@all_vcfs, $vcf_final, $project );
	foreach my $f (@all_vcfs) {
		unlink $f;
	}

	warn "02p freebayes " . abs( $time - time );
	return $vcf_final;
}


sub return_all_beds {
	my ( $project, $intspans ) = @_;
	my @beds;
	foreach my $chr ( @{ $project->getChromosomes() } ) {

		next unless exists $intspans->{ $chr->fasta_name };
		push(
			@beds,
			$project->buffer->intspanToBed(
				$chr, $intspans->{ $chr->fasta_name }
			)
		);
	}
	return \@beds;
}

sub return_all_beds_without_primers {
	my ( $project, $intspans ) = @_;
	my @beds;

	foreach my $chr ( @{ $project->getChromosomes() } ) {
		next unless exists $intspans->{ $chr->fasta_name };
		my $primers         = $project->getPrimers();
		my $intspan_primers = Set::IntSpan::Fast->new();
		foreach my $p (@$primers) {
			$intspan_primers = $intspan_primers->union( $p->intspan_pcr );
		}
		my $intspan_temp =
		$intspans->{ $chr->fasta_name }->diff($intspan_primers);

		push( @beds, $project->buffer->intspanToBed( $chr, $intspan_temp ) );
	}
	return \@beds;
}

sub samtools {
	my ( $project, $patient, $intspans, $fork, $verbose ) = @_;

	$fork = 10;
	my $buffer    = $project->buffer;
	my $dir_out   = $project->getCallingPipelineDir("samtools");
	my $reference = $project->genomeFasta();

	#$reference = "/data-xfs/public-data/HG19/genome/fasta/all.fa";
	#warn $reference
	my $vcf_final = 	  calling_target::getTmpFile( $dir_out, "samtools", "final.vcf" );
	my $samtools    = $buffer->software("samtools");
	my $vcfutil     = $buffer->software("vcfutils");
	my $pm_samtools = new Parallel::ForkManager($fork);
	my $beds        = return_all_beds_without_primers( $project, $intspans );
	my $filter_samtools = "$Bin/filter_samtools.pl";
	if ( scalar(@$beds) == 0 ) {
		print_empty_vcf( $vcf_final, $patient );
		return $vcf_final;
	}

	#my $pr = String::ProgressBar->new( max => scalar(@$beds) );
	my @res_samtools;
	my $c = 0;
	$pm_samtools->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data ) =
			  @_;

			#  $pr->update($c++);

			#$pr->write() if $verbose;
			push( @res_samtools, $data ) if $data->{file};

		}
	);
	my $time     = time;
	my $nb       = 0;
	my $onefile  = $patient->getBamFile();
	my $bcftools = $buffer->software("bcftools");
	my $it       = natatime $fork , @$beds;
	while ( my @vals = $it->() ) {

		#foreach my $bed (@$beds){

		$nb++;
		my $pid = $pm_samtools->start and next;
		my $bed = calling_target::getTmpFile( $dir_out, "samtools", "$nb.bed" );
		open( BED, ">$bed" );
		foreach my $v (@vals) {
			print BED $v . "\n";
		}
		close BED;

	# my $region;
	#($region->{chr},$region->{start},$region->{end}) = split(" ",$bed);
	#	my $chr = $project->getChromosome($region->{chr});
		my $out = calling_target::getTmpFile( $dir_out, "samtools", "$nb.vcf" );

#my $cmd_mpileup =qq{$samtools mpileup $onefile -ugf $reference   -r $free_region -L 8000 2>/dev/null | $bcftools call -mv - |  perl $vcfutil  varFilter -a 5 -d 15  |   $vcfutil  varFilter -a 5 | grep -v "ID=GL" > $out 2>/dev/null};
		#my $cmd_mpileup = qq{$samtools mpileup $onefile -ugf $reference   -l $bed  2>/dev/null | $bcftools call -mv - 2>/dev/null   | grep -v "ID=GL" |  $vcfutil  varFilter -a 5 -d 15 | $bcftools norm -f $reference  - 2>/dev/null |  $filter_samtools > $out 2>/dev/null};
		my $cmd_mpileup = qq{$bcftools mpileup -Ou -f $reference $onefile -R $bed | bcftools call -vmO v | grep -v "ID=GL" |  $vcfutil  varFilter -a 5 -d 15   | $bcftools norm -f $reference  - | $filter_samtools > $out 2>/dev/null};
		warn $cmd_mpileup;
		
		#system ($cmd_mpileup);
		#warn $cmd_mpileup;
		system($cmd_mpileup);
		unlink $bed;
		my %h;
		$h{file} = $out;
		$h{nb}   = $nb;
		$pm_samtools->finish( 0, \%h );

		#}
	}

	$pm_samtools->wait_all_children();

	warn "samtools time :" . abs( time - $time );
	my @all_vcfs =
	  map { $_->{file} } sort { $a->{nb} <=> $b->{nb} } @res_samtools;

	calling_target::concat_vcf( \@all_vcfs, $vcf_final, $project );
	foreach my $f (@all_vcfs) {
		unlink $f;
	}

	warn "samtools " . abs( $time - time );

	return $vcf_final;
}

sub duplicate_regions {
	my ($patient,$noforce) = @_;

	my $project = $patient->project();
	my $bam     = $patient->getBamFile();
	my @bams;
	my $regions;
	my $dir =
	  $project->getVariationsDir("duplicate_region_calling") . "/regions";
	system("mkdir $dir && chmod a+rwx $dir") unless -e $dir;
	my $file = "$dir/" . $patient->name() . ".dup.bed";
#	if ($noforce && -e $file){
#		
#		my @res = `cat $file`;
#		my @zz;
#		chomp(@res);
#		foreach my $r (@res){
#			my ($a,$b,$c) = split(" ",$r);
#			push(@zz,"$a:$b-$c");
#		}
#		return \@zz;
#	}
	open( REG, ">$dir/" . $patient->name() . ".dup.bed" );

	foreach my $chr ( @{ $project->getChromosomes } ) {
		my $cn      = $chr->name();
		my $intspan = $chr->getIntSpanCaptureForCalling(200);
		next unless $intspan;
		next unless $intspan->as_array();
		my $intspan_res->{$cn} = Set::IntSpan::Fast->new();

		my $callback = sub {
			my ( $seqid, $pos, $pileup ) = @_;
			my $nb  = 0;
			my $nbd = 0;
			return if scalar(@$pileup) < 50;
			for my $p (@$pileup) {
				my $alignment = $p->alignment;
				$nbd++  if $alignment->qual ==0;    #&&  $alignment->get_tag_values("XA");
				$nb++;
			}

			#	 warn $pos." $nb :: $nbd " if $nbd > 50;
			if ( $nbd > 0.75 * $nb ) {
				$intspan_res->{$cn}->add_range( $pos, ( $pos + 150 ) );
			}

		};

		my $sam = Bio::DB::Sam->new(
			-fasta => $project->getGenomeFasta()
			,    #"/data-xfs/public-data/HG19//genome/fasta/all.fa",
			-bam => $bam
		);

		my $iter = $intspan->iterate_runs();
		while ( my ( $from, $to ) = $iter->() ) {
			$sam->fast_pileup( $chr->fasta_name . ":$from-$to", $callback );
		}
		push( @$regions,	map { $chr->fasta_name . ":" . $_ }   split( ",", $intspan_res->{$cn}->as_string ) );
		my $iter2 = $intspan_res->{$cn}->iterate_runs();
		while ( my ( $from, $to ) = $iter2->() ) {
			print REG $chr->fasta_name . "\t" . $from . "\t" . $to . "\n";
		}
	}

	close(REG);
	return $regions;
}

sub duplicate_region_calling {
	my ( $project, $patient, $intspans,  $fork, $verbose ) = @_;
	my $dir_out = $project->getCallingPipelineDir("duplicate_region_calling");
	my $vcf_final =	  calling_target::getTmpFile( $dir_out, $patient->name, "dup.vcf" );
	my $regions = duplicate_regions( $patient, $intspans );
	unless (@$regions) {

		print_empty_vcf( $vcf_final, $patient );
		return $vcf_final;
	}
	return unless @$regions;

#my $regions = ["chr11:118895063-118898583","chrX:153770668-153793248"];
#my $regions = ["chr7:74773962-74789315","chrX:153770668-153793248", "chr15:43891596-43910998","chrX:101096262-101096273","chrX:101093768-101093792","chrX:101093344-101093352","chrX:101093197-101093202","chr6:32012128-32012523","chr6:32011754-32011936","chr6:32011184-32011363","chr6:32010922-32011125","chr6:32010698-32010888","chr6:32010446-32010638","chr6:32010013-32010169","chr6:32009759-32009980","chr6:32009518-32009741","chr6:32008902-32009257"];
	my $buffer = $project->buffer;

	my $reference = $project->genomeFasta();

	my $javac = $buffer->getSoftware('java');

	#$javac ="java";
	my $gatk     = $buffer->getSoftware('gatk');
	my $bcftools = $buffer->software("bcftools");
	my $vcfutil  = $buffer->software("vcfutils");
	my $samtools = $buffer->software("samtools");
	my $gatk4  = $project->getSoftware('gatk4');
	my $vcf_uni =	  calling_target::getTmpFile( $dir_out, $patient->name, "dup.vcf" );
	my $bamtmp =	  calling_target::getTmpFile( $dir_out, $patient->name, "tmp.bam" );
	my $bed = calling_target::getTmpFile( $dir_out, $patient->name, "tmp.bed" );
	my $bam = $patient->getBamFile();
	my $cmd2 =qq{ | perl -lane 'if(\$_=~/^@/){print \$_} else{\@t=split(" ",\$_); \$t[4] = 30 if \$t[4] eq 0;print join("\t",\@t);}' | samtools view -bS - | samtools sort - > $bamtmp && samtools index $bamtmp};
	my $cmd1 = qq{samtools view -h $bam };

	my $cmd = $cmd1 . join( " ", @$regions ) . $cmd2;
	warn $cmd;

	system($cmd);
	die() unless -e $bamtmp . ".bai";

	my $low_calling_filter = "";

	open( BED, ">$bed" );
	foreach my $r (@$regions) {
		my ( $chr,   $pos ) = split( ":", $r );
		my ( $start, $end ) = split( "-", $pos );
		print BED $chr . "\t" . "$start\t$end\n";
	}
	close BED;

	my $recal              = $patient->getRecalFile();
	my $recal_string       = "";
	my $bam_file_string_hc = " -I " . $patient->getBamFile();
	my $gatk_region        = "";
	
	my $cmd_uni =qq{$gatk4  HaplotypeCaller   -R $reference         -I $bamtmp   -O $vcf_uni -L $bed };

	system($cmd_uni);
	unlink $bamtmp;
	my @res = `cat $vcf_uni`;
	warn $res[0];
	$res[0] .= "##REGION_DUPLICATE=" . join( ",", @$regions ) . "\n";
	open( VCF, ">", $vcf_final );
	print VCF @res;
	close VCF;
	#my @res2 = `cat $vcf_final`;
	#print @res2;
	warn $vcf_final;
	return $vcf_final;
}

sub duplicate_region_calling2 {
	my ( $project, $patient, $intspans,  $fork, $verbose ) = @_;
	my $dir_out = $project->getCallingPipelineDir("duplicate_region_calling");
	my $vcf_final =	  calling_target::getTmpFile( $dir_out, $patient->name, "dup.vcf" );
	my $regions = duplicate_regions( $patient, $intspans,1);
	unless (@$regions) {

		print_empty_vcf( $vcf_final, $patient );
		return $vcf_final;
	}
	return unless @$regions;

#my $regions = ["chr11:118895063-118898583","chrX:153770668-153793248"];
#my $regions = ["chr7:74773962-74789315","chrX:153770668-153793248", "chr15:43891596-43910998","chrX:101096262-101096273","chrX:101093768-101093792","chrX:101093344-101093352","chrX:101093197-101093202","chr6:32012128-32012523","chr6:32011754-32011936","chr6:32011184-32011363","chr6:32010922-32011125","chr6:32010698-32010888","chr6:32010446-32010638","chr6:32010013-32010169","chr6:32009759-32009980","chr6:32009518-32009741","chr6:32008902-32009257"];
	my $buffer = $project->buffer;

	my $reference = $project->genomeFasta();

	my $javac = $buffer->getSoftware('java');

	#$javac ="java";
	my $gatk     = $buffer->getSoftware('gatk');
	my $bcftools = $buffer->software("bcftools");
	my $vcfutil  = $buffer->software("vcfutils");
	my $samtools = $buffer->software("samtools");
	my $vcf_uni =	  calling_target::getTmpFile( $dir_out, $patient->name, "dup.vcf" );
	my $bamtmp =	  calling_target::getTmpFile( $dir_out, $patient->name, "tmp.bam" );
	my $bed = calling_target::getTmpFile( $dir_out, $patient->name, "tmp.bed" );
	my $bam = $patient->getBamFile();
	my $cmd2 =qq{ | perl -lane 'if(\$_=~/^@/){print \$_} else{\@t=split(" ",\$_); \$t[4] = 30 if \$t[4] eq 0;print join("\t",\@t);}' | samtools view -bS - | samtools sort - > $bamtmp && samtools index $bamtmp};
	my $cmd1 = qq{samtools view -h $bam };

	my $cmd = $cmd1 . join( " ", @$regions ) . $cmd2;
	warn $cmd;

	system($cmd);
	die() unless -e $bamtmp . ".bai";

	my $low_calling_filter = "";

	open( BED, ">$bed" );
	foreach my $r (@$regions) {
		my ( $chr,   $pos ) = split( ":", $r );
		my ( $start, $end ) = split( "-", $pos );
		print BED $chr . "\t" . "$start\t$end\n";
	}
	close BED;

	my $recal              = $patient->getRecalFile();
	my $recal_string       = "";
	my $bam_file_string_hc = " -I " . $patient->getBamFile();
	my $gatk_region        = "";
	my $cmd_uni =qq{$javac -jar $gatk -T HaplotypeCaller   -R $reference     -rf BadCigar    -I $bamtmp   -o $vcf_uni -L $bed};

	warn $cmd_uni;
	system($cmd_uni);
	unlink $bamtmp;
	my @res = `cat $vcf_uni`;
	warn $res[0];
	$res[0] .= "##REGION_DUPLICATE=" . join( ",", @$regions ) . "\n";
	open( VCF, ">", $vcf_final );
	print VCF @res;
	close VCF;
	#my @res2 = `cat $vcf_final`;
	#print @res2;
	warn $vcf_final;
	return $vcf_final;
}





#sub mutect2 {
#	my ( $project, $patient, $intspans, $fork, $verbose ) = @_;
#	my $buffer    = $project->buffer;
#	my $dir_out   = $project->getCallingPipelineDir("mutect2");
#	my $reference = $project->genomeFasta();
#	my $name = $patient->name();
#	#my $javac = $buffer->getSoftware('java');
#	
#	#$javac ="java";
#	my $gatk      = $buffer->getSoftware('gatk4');
#	#my $freebayes = $buffer->software("freebayes");
#	my $bcftools  = $buffer->software("bcftools");
#	#my $vcfutil   = $buffer->software("vcfutils");
#	#my $samtools  = $buffer->software("samtools");
#	my $mutect_out = calling_target::getTmpFile( $dir_out, $project->name, "mutect.vcf" );
#	my $out = calling_target::getTmpFile( $dir_out, $project->name, "hc.vcf" );
#	#my $out1 = 	  calling_target::getTmpFile( $dir_out, $project->name, "hc.1.vcf" );
#	#my $low_calling_filter = "";
#	my $beds = return_all_beds( $project, $intspans );
#
#	if ( scalar(@$beds) == 0 ) {
#		print_empty_vcf( $out, $patient );
#		return $out;
#	}
#	my $bed1 = calling_target::getTmpFile( $dir_out, $project->name, "bed" );
#	open( BED, ">$bed1" );
#	print BED join( "\n", @$beds );
#	print BED "\n";
#	close BED;
#	#my $unified_min_alternate = 0.07;
#
#
##	my $recal              = $patient->getRecalFile();
##	my $recal_string       = "";
#	my $bam_file_string_hc = " -I " . $patient->getBamFile();
#	my $gatk_region        = "-L $bed1";
#	
#	my $cmd_mutect2 = qq{ $gatk Mutect2 --base-quality-score-threshold 25 --dont-use-soft-clipped-bases -tumor $name -R $reference $gatk_region   $bam_file_string_hc   -O $mutect_out ;  $bcftools norm -f $reference -m -any $mutect_out -o $out  };
#	warn $cmd_mutect2;
#	system($cmd_mutect2);
#	#system( $cmd_uni   . qq{&& $bcftools  filter -i " $filter_bcftools " $out1 > $out && rm $out1} );
#	return $out;
#
#}

sub mutect2 {
	my ( $project, $patient, $intspans, $fork, $verbose ) = @_;

	$fork = 10;
	my $buffer    = $project->buffer;
	my $dir_out   = $project->getCallingPipelineDir("mutect2");
	my $reference = $project->genomeFasta();
	my $name = $patient->name();
	#$reference = "/data-xfs/public-data/HG19/genome/fasta/all.fa";
	#warn $reference
	my $vcf_final = 	 calling_target::getTmpFile( $dir_out, $project->name, "mutect.vcf" );
	my $gatk      = $buffer->getSoftware('gatk4');
	
	my $beds = return_all_beds( $project, $intspans );
	
	if ( scalar(@$beds) == 0 ) {
		print_empty_vcf( $vcf_final, $patient );
		return $vcf_final;
	}

	#my $pr = String::ProgressBar->new( max => scalar(@$beds) );
	my @res_samtools;
	my $c = 0;
	my $pm = new Parallel::ForkManager($fork);
	$pm->run_on_finish(
		sub {
			my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data) = @_;
			push(@res_samtools, $data);
		}
	);
	my $time     = time;
	my $nb       = 0;
	my $onefile  = $patient->getBamFile();
	my $bcftools = $buffer->software("bcftools");
	
	$fork = scalar(@$beds) if scalar(@$beds) < $fork;
	my $it       = natatime $fork , @$beds;
	my $bam_file_string_hc = " -I " . $patient->getBamFile();
	
	$project->disconnect();
	while ( my @vals = $it->() ) {
 	 	my $pid = $pm->start and next;
		my $res;
		$res->{start}   = 'started';
		#foreach my $bed (@$beds){

		$nb++;
		my $bed = calling_target::getTmpFile( $dir_out, "mutect", "$nb.bed" );
		open( BED, ">$bed" );
		foreach my $v (@vals) {
			print BED $v . "\n";
		}
		close BED;
		my $mutect_out = calling_target::getTmpFile( $dir_out, "mutect.tmp.", "$nb.vcf" );
		my $vcf_out = calling_target::getTmpFile( $dir_out, "mutect", "$nb.vcf" );
		my $gatk_region        = "-L $bed";
		my $cmd_mutect2 = qq{ $gatk Mutect2 --base-quality-score-threshold 25 --dont-use-soft-clipped-bases -tumor $name -R $reference $gatk_region   $bam_file_string_hc   -O $mutect_out ;  $bcftools norm -f $reference -m -any $mutect_out -o $vcf_out  };
		system($cmd_mutect2);
		unlink $bed;
		$res->{file} = $vcf_out;
		$res->{nb}   = $nb;
	 	$pm->finish(0, $res);
	}
	sleep(3); 
	$pm->wait_all_children();
	
	$project->disconnect;

	warn "samtools time :" . abs( time - $time );
	my @all_vcfs =
	  map { $_->{file} } sort { $a->{nb} <=> $b->{nb} } @res_samtools;

	calling_target::concat_vcf( \@all_vcfs, $vcf_final, $project );
	foreach my $f (@all_vcfs) {
		unlink $f;
	}

	warn "samtools " . abs( $time - time );

	return $vcf_final;
}




sub lofreq {
	my ( $project, $patient, $intspans, $fork, $verbose ) = @_;
	my $buffer    = $project->buffer;
	my $dir_out   = $project->getCallingPipelineDir("lofreq");
	my $reference = $project->genomeFasta();
	my $name = $patient->name();
	#my $javac = $buffer->getSoftware('java');
	
	#$javac ="java";
	my $lofreq      = $buffer->getSoftware('lofreq');
	my $out = calling_target::getTmpFile( $dir_out, $project->name, "hc.vcf" );
	my $beds = return_all_beds( $project, $intspans );

	if ( scalar(@$beds) == 0 ) {
		print_empty_vcf( $out, $patient );
		return $out;
	}
	my $bed1 = calling_target::getTmpFile( $dir_out, $project->name, "bed" );
	open( BED, ">$bed1" );
	print BED join( "\n", @$beds );
	print BED "\n";
	close BED;

	my $bam_file =  $patient->getBamFile();
	my $gatk_region        = "-l $bed1";
	my $cmd_lofreq = qq{ $lofreq call-parallel -f $reference $gatk_region  --pp-threads 20 $bam_file  -o $out  };
	warn $cmd_lofreq;
	system($cmd_lofreq);
	#system( $cmd_uni   . qq{&& $bcftools  filter -i " $filter_bcftools " $out1 > $out && rm $out1} );
	return $out;

}



1;

