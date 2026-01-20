package dragen_util;
use strict;
use Data::Dumper;

sub get_fastq_file {
	my ($patient,$dir_pipeline, $dir_fastq) = @_;
	my $name=$patient->name();
	$dir_fastq = $patient->getSequencesDirectory() unless $dir_fastq;
	my $files_pe1 = file_util::find_file_pe($patient,"",$dir_fastq,1);
	return unless @$files_pe1;
	my $cmd;
	my @r1;
	my @r2;
	my $find_I1;
	my $indice2 = "R2";
	#if exists 
	my $r2 = "R2";
	foreach my $cp (@$files_pe1) {
		
		if (exists $cp->{R3}){
			$r2 = "R3";
		}
		$find_I1 ++ if $cp->{R1} =~ /_I1_0/;;
		next if $cp->{R1} =~ /_I1_0/;
		die($dir_fastq." ".Dumper $cp) if ($cp->{R2} !~ /_R2_/ and $cp->{R2} !~ /_R2./);
		my $file1 = $dir_fastq."/".$cp->{R1};
		die($file1) unless -e  $file1;
		my $file2 = $dir_fastq."/".$cp->{$r2};
		die($file2) unless -e  $file1;
		push(@r1,$file1);
		push(@r2,$file2);##
	}
	my $cmd1 = join(" ",@r1);
	my $cmd2 = join(" ",@r2);
	my $fastq1 = $dir_pipeline."/".$patient->name."_S1_L001_R1_001.fastq.gz";
	my $fastq2 = $dir_pipeline."/".$patient->name."_S1_L001_R2_001.fastq.gz";
	#if ($step eq "align"){
		system "cat $cmd1 > $fastq1";# unless -e $fastq1;
		system "cat $cmd2 > $fastq2";# unless -e $fastq2;
	#}
	
	if ($find_I1){
		
		my @r1;
		my @r2;
		foreach my $cp (@$files_pe1) {
			next if $cp->{R1} =~ /_R1_/;
			die() if $cp->{R2} !~ /_I2_/;
			my $file1 = $dir_fastq."/".$cp->{R1};
			die($file1) unless -e  $file1;
			my $file2 = $dir_fastq."/".$cp->{R2};
			die($file2) unless -e  $file1;
			push(@r1,$file1);
			push(@r2,$file2);##
		}
		my $cmd1 = join(" ",@r1);
		my $cmd2 = join(" ",@r2);
		my $fastq1 = $dir_pipeline."/".$patient->name."_S1_L001_I1_001.fastq.gz";
		my $fastq2 = $dir_pipeline."/".$patient->name."_S1_L001_I2_001.fastq.gz";
		
		#if ($step eq "align"){
			system "cat $cmd1 > $fastq1";# unless -e $fastq1;
			system "cat $cmd2 > $fastq2";# unless -e $fastq2;
		#}	
	
	}
	return  ($fastq1,$fastq2);
}


sub get_capture_file {
	my ($patient,$capture_file) = @_;
	my $project = $patient->project;
	my $capture_file_prod  = $patient->getCapture->gzFileName();
	open (ZCAT,"zcat $capture_file_prod |");
	open (CAP,">$capture_file");
	while(<ZCAT>){
		chomp();
		my ($a,$b,$c) = split(" ",$_);
		my $chr = $project->getChromosome($a);
		print CAP $chr->fasta_name."\t".$b."\t".$c."\n";
	}
	close ZCAT;
	close CAP;
	my $bgzip = $project->buffer->software("bgzip");
	my $tabix = $project->buffer->software("tabix");
	system("bgzip $capture_file && tabix -p bed $capture_file.gz");
	die() unless -e "$capture_file.gz.tbi";
	return "$capture_file.gz";
	
}
1;