#!/usr/bin/perl
use strict;
use Getopt::Long;
use Data::Dumper;
use Parallel::ForkManager;

my $projectName;
my $input_file;
my $fork;

GetOptions(
	'file=s'   => \$input_file,
	'fork=s'   => \$fork,
	'project=s' => \$projectName,
);

my @chromosomes = `samtools idxstats $input_file | cut -f 1`;
chomp(@chromosomes);

warn Dumper(@chromosomes);

my $chr_syno ={
		1=> "chr1",
		2=>"chr2",
		3=>"chr3",
		4=>"chr4",
		5=>"chr5",
		6=>"chr6",
		7=>"chr7",
		8=>"chr8",
		9=>"chr9",
		10=>"chr10",
		11=>"chr11",
		12=>"chr12",
		13=>"chr13",
		14=>"chr14",
		15=>"chr15",
		16=>"chr16",
		17=>"chr17",
		18=>"chr18",
		19=>"chr19",
		20=>"chr20",
		21=>"chr21",
		22=>"chr22",
		X=>"chrX",
		Y=>"chrY",
		MT=>"chrM",
};



	die() unless -e $input_file;
	my $output_file = $input_file;
	my $backup_file = $input_file;
	$output_file =~ s/bam/new.bam/;
	$backup_file =~ s/bam/old.bam/;
my $samtools = "/software/bin/samtools"; 	
;
 my $pm2 = new Parallel::ForkManager($fork);
 my @tmp_file;
 
 foreach my $chr (@chromosomes){
my $bamout = "/data-beegfs/tmp/toto.$chr.bam";
push(@tmp_file,$bamout);
my $pid      = $pm2->start() and next;	
open (BAM,"$samtools  view -h $input_file $chr|");
my $pid = open(WRITEME, "| $samtools view -bS -   > $bamout") or die "Couldn't find file $input_file \n $!\n";
while (<BAM>){
	chomp();
	my $line = $_;
	my @t = split(" ",$line);
	
	if ($line =~ /^\@SQ/){
		
		my ($sn,$ch) = split(":",$t[1]);
		if (exists $chr_syno->{$ch}){
			$t[1] = "SN:".$chr_syno->{$ch};
		}
		$line =  join("\t",@t);
		
	}

	else {
		
		my $c1 = $t[2];
			if (exists $chr_syno->{$c1}){
			$t[2] = $chr_syno->{$c1};
		}
		my $c2 = $t[6];
			if (exists $chr_syno->{$c2}){
			$t[6] = $chr_syno->{$c2};
		}
	
	}
	#warn $line;
	$line =  join("\t",@t)."\n";
	print WRITEME $line;
}
close(BAM);
close(WRITEME);

 	$pm2->finish(0);
	}
$pm2->wait_all_children;
 my $split = join(" ",@tmp_file);
 my $final = "/data-beegfs/tmp/final.bam";
 
  system("sambamba  merge -t $fork $final $split");
 
#system ("mv $input_file $backup_file");
#system ("mv $output_file $input_file");
warn "start index !!!";
#system ("$samtools index $input_file");
#unlink ("$backup_file");
exit(0);

#

sub ischrornot {
	my ($file) = @_;
	my @res = 'samtools view -H $file';
	my ($find) = grep {$_=~ /SN:chr/} @res;
	return $find;
}