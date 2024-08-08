package epi2me_util;
use strict;
use Data::Dumper;

sub get_fastq_file {
	my ($patient,$dir_pipeline,$dir_fastq) = @_;
	my $name=$patient->name();
	$dir_fastq = $patient->getSequencesDirectory() unless $dir_fastq;
	if($patient->isGenome()){
		return $dir_fastq;
	}
	my $files_pe1 = file_util::find_file_epi2me($patient,$dir_fastq,1);
	my $path2 = $dir_fastq."/".lc($patient->barcode()."/");
	$dir_fastq = $path2  unless defined $files_pe1;
	$files_pe1 =  file_util::find_file_epi2me($patient,$path2,1) unless defined $files_pe1;
	return unless @$files_pe1;
	my @files = map { $dir_fastq."/".$_} @$files_pe1;
	my $cmd1 = join(" ",@files);
	my $fastq1 = $dir_pipeline."/".$patient->name.".fastq.gz";
	
	system "cat $cmd1 > $fastq1";# unless -e $fastq1;
	return  ($fastq1);
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