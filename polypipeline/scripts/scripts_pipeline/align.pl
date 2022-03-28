#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use File::Basename;
use  File::Temp;
use Sys::Hostname;

my $filein;
my $dir;
my $file_bed;
my $project_name;
my $chr_name;
my $fork = 1;
use List::MoreUtils qw(firstidx);
use Vcf;

my @tmp_files;

my $dir_out;
my $log_file;
my $norecal;
my $file1;
my $file2;
my $aligner;
my $fileout;
my $name;
my $prog;
my $lane ="";
my $mode ;
my $verbose;
my $nosort;
my $notest;
GetOptions(
	"log=s" =>\$log_file,
	"project=s" =>\$project_name,
	"fork=s" =>\$fork,
	"bam=s"  =>\$fileout,
	"file1=s"  =>\$file1,
	"file2=s"  =>\$file2,
	"method=s" => \$prog,
	"name=s" => \$name,
	"lane=s" => \$lane,
	"mode=s" => \$mode,
	"verbose=s" => \$mode,
	"nosrot=s" => \$nosort,
	"nosort=s" => \$nosort, 
	"notest=s" => \$notest,
);
 
if ($log_file){
	open (STDOUT,">>".$log_file);
}
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
 $mode = "pe" unless $mode;
if ($prog =~/_UMI/){
	my @dir_out  = split("/",$fileout);
	pop(@dir_out);
	my $dd =join("/",@dir_out);
	my($f1,$f2) = mark_UMI($file1,$file2,$dd);
	$file1 = $f1;
	$file2 = $f2;
	$prog =~s/_UMI//;
	
}
warn hostname;

$prog = "bwa" if $prog eq "umi_dedup";
$prog = "star" if $prog eq "star_umi";
$prog = "bwa" if $prog eq "bwa_xths2";
my $cmd = $buffer->software($prog);
warn $cmd;

die("problem with method:$prog (bwa,star,hisat2,bwa_xths2") if ($prog ne "star" && $prog ne "hisat2" && $prog ne "bwa" && $prog ne "bowtie2" && $prog ne "umi_dedup" && $prog ne "bwa_xths2");

die("file1 not found : -$file1-") unless -e $file1;
if ($file2){
	die("file2 not found : $file2") unless -e $file2;
}
else {
	$mode ="frag";
}

my $bgzip = $buffer->software("bgzip");
if($file1 !~ /fastq.gz/){
	system("$bgzip -f -@ $fork $file1");
	system("$bgzip -f -@ $fork $file2");
	$file1 .= ".gz";
	$file2 .= ".gz";
} 

my $index;
my $ref_root =  $project->dirGenome();#$buffer->config->{public_data}->{$project->getVersion()};



my $index = $project->getGenomeIndex($prog);#
$index .= "/all.fa" if $prog eq "bwa";#$ref_root.$buffer->index($prog);




my $file;
my $path;
my $suffix;

 ($file,$path,$suffix) = fileparse($fileout);
warn $file." ".$path;

#my $pathtmp  = "/data-ibrixfs";                        

 my $dir_tmp_cat = File::Temp->newdir( "$file.XXXXXXX", 
                        DIR => $path,
                        CLEANUP => 1 );
my $dir_tmp = $dir_tmp_cat->dirname();

 system("chmod a+w $dir_tmp");
my $file_out_tmp = $dir_tmp."/".$file;
my $file_out_tmp2 = $dir_tmp."/".$file.".tmp.bam";
my $samtools =$buffer->software('samtools');
my $sambamba =$buffer->software('sambamba');
my $pgzip = $buffer->software('pigz');

my $cmd_test = "";
$cmd_test = "pigz -t $file1 $file2 -p $fork && " if -e $pgzip;

my $cmd_sort = " $samtools view -Sb - > $file_out_tmp2 && $samtools  sort -@ $fork $file_out_tmp2 -o  $file_out_tmp  && rm $file_out_tmp2  ";
$cmd_sort = " $samtools view -Sb - > $file_out_tmp" if $nosort; 

my $readgroup = qq{\@RG\\tID:$name\\tSM:$name\\tPL:ILLUMINA\\tDS:$name\\tPG:$prog};
my $readgroupstar = "ID:".$name." SM:".$name." PL:ILLUMINA DS:".$name." PG:".$prog;
my $cmd_move = "mv $file_out_tmp $fileout && test -e $file_out_tmp.bai &&  mv $file_out_tmp.bai $fileout.bai";
my $cmd_index = "$samtools index $file_out_tmp";                     
#$sammem $reference_bwa $file1 $file2 -R readgroup -M

my $align_command ={
	pe=>{
		bwa =>  "cd $dir_tmp && $cmd mem $index  $file1 $file2 -t $fork  -R '$readgroup' -M  | $cmd_sort  && $cmd_move",
		hisat2 => "cd $dir_tmp && $cmd --dta-cufflinks -x $index -1 $file1  -2 $file2 --rg-id $name --rg SM:$name  -p $fork | $cmd_sort && $cmd_move",
		star =>  "cd $dir_tmp && $cmd --runThreadN $fork --genomeDir $index  --readFilesIn  $file1 $file2 --outSAMunmapped Within --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate  --outSAMattrRGline $readgroupstar && mv Aligned.sortedByCoord.out.bam $file_out_tmp && $cmd_index && $cmd_move",
		bowtie2 =>  "cd $dir_tmp && $cmd  -x  $index -1 $file1 -2 $file2 --threads $fork| $cmd_sort  && $cmd_move",
	},
	frag =>{
		bwa =>  "cd $dir_tmp && $cmd mem $index  $file1  -t $fork  -R '$readgroup' -M  | $cmd_sort  && $cmd_move",
		hisat2 => "cd $dir_tmp && $cmd --dta-cufflinks -x $index  -U $file1 -p $fork | $cmd_sort && $cmd_move",
		star =>  "cd $dir_tmp && $cmd --runThreadN $fork --genomeDir $index  --readFilesIn  $file1   --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMattrRGline $readgroupstar && mv Aligned.sortedByCoord.out.bam $file_out_tmp && $cmd_index && $cmd_move"
	},
	xths2 => {
		bwa=>  "cd $dir_tmp && $cmd mem -C $index  $file1 $file2 -t $fork  -R '$readgroup' -M  | $cmd_sort  && $cmd_move",
	}

	} ;

warn "@@@@@@@@@@@";
warn $align_command->{$mode}->{$prog};
warn $cmd_test.$align_command->{$mode}->{$prog};
my $cmd_align = $cmd_test.$align_command->{$mode}->{$prog};
$cmd_align = $align_command->{$mode}->{$prog} if $notest;
system( $cmd_align);


sub mark_UMI{
	my ($read1,$read2,$dir_out) = @_;
	my %humi;
my @part = split("/",$read1);
my $name = pop(@part);
warn $name;
 $name =~ s/_R[1,2,3]_/_R\*_/;
warn $name;
my $path = join("/",@part);
#warn $path.$name;
my @l = `ls $path/index/$name`;
chomp(@l);
die() if scalar(@l) ne 1;
my $umi = $l[0];


#open my $fumi, " zcat ", $umi." " or die "Unable to open $umi for reading : $!";
  open (FUMI," zcat $umi | ") or die();;
while(1){
	my $header = <FUMI>;
	last unless $header;
	chomp($header);
	my $seq = <FUMI>;
	my $plus = <FUMI>;
	my $quality = <FUMI>;
	chomp($seq);
	my ($h1,$h2) = split(" ",$header) ;
		if (length($seq) ne 8 ) {
		$seq = "NNNNNNNN";
	}
	$humi{$h1} = $h1.":".$seq;
	
}

close (FUMI);
warn scalar(keys %humi);
my @fs = ($read1,$read2);
	my @fouts;
foreach my $f (@fs){
	my @part = split("/",$f);
	my $name = pop(@part);
	$name =~ s/\.gz//;
	my $fout = $name;
	open (FUMI," zcat $f | ") or die();;
	warn $fout;
	open (OUT," >$dir_out/$fout ") ;;

while(1){
	my $header = <FUMI>;
	last unless $header;
	chomp($header);
	my $seq = <FUMI>;
	my $plus = <FUMI>;
	my $quality = <FUMI>;
	#chomp($seq);
	my ($h1,$h2) = split(" ",$header) ;
	next unless exists $humi{$h1};
	#warn $header unless exists $humi{$h1};
	if ($h2 =~/^3:/){
		$h2=~s/3:/2:/;
	}
	print OUT $humi{$h1}." ".$h2."\n";
	print OUT $seq;
	print OUT $plus;
	print OUT $quality;
	#warn $humi{$h1};
}
	close FUMI;
	my $bgzip = $buffer->software("bgzip");
	system("$bgzip -f -@ $fork $dir_out/$fout");
	push(@fouts, "$dir_out/$fout.gz");
}
	warn "END UMI";
	return (@fouts);
}


