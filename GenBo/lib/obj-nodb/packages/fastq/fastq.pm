package fastq;
use strict;
use Data::Printer;
use File::Util;
use Carp;
use Data::Dumper;



sub find_paired_files_umi {
	my ($files,$dir) = @_;
	my @r1 = grep{$_  =~ /_R1_/} @$files;
	my @all;
	my @couple;
	foreach my $r1 (@r1){
		my $d;
		$d->{R1} = $dir."/".$r1;
		
		my $r2 = $r1;
		if (-e $dir."/".$r2){
			$r2 =~ s/_R1_/_R2_/;
			$d->{R2} = $dir."/".$r2;
		}
		else {
			confess($r1 ." don t find R2" );
		}
		my $r3 = $r1;
		$r3 =~ s/_R1_/_R3_/;
		if ( -e $dir."/".$r3) {
			$d->{R3} = $dir."/".$r3;
		}
		push(@couple,$d);
	}
	
	return \@couple;
}


my %cached_dir;

sub find_file_pe {
	my ($patient,$separator) = @_;
	my $dir = $patient->getSequencesDirectory();
	my @names;
	my $couple;
	push(@names,$patient->name);
	push(@names,$patient->barcode) if length($patient->barcode)>1 ;
	unless (exists $cached_dir{$dir}){
		$cached_dir{$dir}= [];
		opendir(DIR,$dir);
		my @allFiles= readdir(DIR);
		$cached_dir{$dir} = \@allFiles;
	}
	NAME: foreach my $name (@names){
	my @pattern = ("^".$name."_[ATGC][ATGC][ATGC]","^".$name."_S[1-9]+","^".$name."_","$name");
	foreach my $find (@pattern){
		my (@titi) = grep { /$find/} grep { /fastq/} @{$cached_dir{$dir}} ;
		
		if (@titi) {
			$couple = find_paired_files_umi(\@titi,$dir);
			last NAME if ($couple);
		}
	}
	}
	confess("NO fastq file for : -".$patient->name()."- ".$patient->barcode." ".$dir) unless $couple;	
	return $couple if scalar(@$couple)>0;
	die();
}





1;
