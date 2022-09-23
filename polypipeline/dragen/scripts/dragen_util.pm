package dragen_util;
use strict;

sub get_fastq_file {
	my ($patient,$dir_pipeline) = @_;
	
	my $name=$patient->name();
	
	my $files_pe1 = file_util::find_file_pe($patient,"");
	my $cmd;
	my @r1;
	my @r2;
	foreach my $cp (@$files_pe1) {
		my $file1 = $patient->getSequencesDirectory()."/".$cp->{R1};
		my $file2 = $patient->getSequencesDirectory()."/".$cp->{R2};
		push(@r1,$file1);
		push(@r2,$file2);
	}
	my $cmd1 = join(" ",@r1);
	my $cmd2 = join(" ",@r2);
	my $fastq1 = $dir_pipeline."/".$patient->name.".R1.fastq.gz";
	my $fastq2 = $dir_pipeline."/".$patient->name.".R2.fastq.gz";
	#if ($step eq "align"){
		system "cat $cmd1 > $fastq1";# unless -e $fastq1;
		system "cat $cmd2 > $fastq2";# unless -e $fastq2;
	#}
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