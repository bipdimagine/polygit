#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages/";

use GBuffer;
use GenBoProject;
use Getopt::Long;
use Carp;
use strict;
use Set::IntSpan::Fast::XS;
use Data::Dumper;
use Data::Printer;
use Storable qw(store retrieve freeze);
use file_util;
use Term::Menus;
use IO::Prompt;
use GenBoNoSqlLmdb;
use Bio::Seq::Quality;
use Bio::SeqIO::fastq;
use GenOO::Data::File::FASTQ;
use Parallel::ForkManager;
use JSON::XS;
my $buffer = GBuffer->new();
my $patient_name;
my $project_name;
my $fork;
my $out;

GetOptions(
	'project=s' => \$project_name,
	'patient=s' => \$patient_name,
	'fork=s' => \$fork,
	'out=s'  => \$out,
);



#my $patient_name = "2006N080737_LEBMI";
my $project = $buffer->newProject(-name=>"$project_name");
my $patient = $project->getPatient($patient_name);
my $method = $patient->alignmentMethod();
my $dir = $project->getAlignmentPipelineDir($method."/".$patient->name);
my $final2 = "$dir/list_file.txt";
$fork=40 unless $fork;

 my $json_file = $dir."/files.json";
if($out){
	print $json_file;
	exit(0);
}

	my $picard = $buffer->software("picard");
	my $java = $buffer->software("java");
	my $pigz = $buffer->software("pigz");
	my $dir_fastq = $patient->getRun->fastq_dir();
	my $files_pe1 = file_util::find_file_pe($patient,$dir_fastq);

#java -jar picard.jar FastqToSam \
#      F1=file_1.fastq \
#      O=fastq_to_bam.bam \
#      SM=for_tool_testing 

#system(qq{$picard FastqToSam F1=$fastq1  F2=$fastq2 SM="$rg" ID="$rg" PL=illumina "})
#my $files_pe1 = file_util::find_file_pe_umi($patient,"");

my $pm = new Parallel::ForkManager(scalar(@$files_pe1));

my $generic = "$dir".$patient->name;

my @global_files;
my $rg= $patient->name."_1";
my ($dirin) = $patient->getSequencesDirectory();

	$pm->run_on_finish(
    	sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    		push(@global_files,@{$data->{files}});
    }
    );


my $fq =0;

foreach my $cp (@$files_pe1){
$fq ++;
my $pid = $pm->start and next;
$| =1;
my @files;
my $nb_file =0;
my $fastq1 = $generic.".$nb_file.$fq"."_R1.fastq.gz";
my $fastq2 = $generic.".$nb_file.$fq"."_R2.fastq.gz";
my $ubam = $generic.".$nb_file.$fq".".u.bam";
open(my $ofastq1,"| $pigz -p 5  > $fastq1");
open(my $ofastq2,"| $pigz -p 5 > $fastq2");
#open(my $ofastq1,"> $fastq1");
#open(my $ofastq2,">$fastq2");
open(my $obam,"|-","$java -Djava.io.tmpdir=/tmp -jar $picard SortSam -I /dev/stdin  -VERBOSITY ERROR -SORT_ORDER queryname -O $ubam -QUIET true");
#open(my $obam,"|-","samtools view -S -b - >$ubam");
print $obam qq{\@HD	VN:1.6	SO:queryname	GO:none\n\@RG	ID:$rg	SM:$rg	LB:$rg	PL:illumina\n};

push(@files,{fastq1=>$fastq1,fastq2=>$fastq2,ubam=>$ubam});

my $R1 = $dirin."/".$cp->{R1};#"/data-isilon/data/sequences/ILLUMINA/NOVASEQ/IMAGINE/run103.saved/2006N080737_LEBMI_S24_L001_R1_001.fastq.gz";
my $R2 = $dirin."/".$cp->{R2};#"/data-isilon/data/sequences/ILLUMINA/NOVASEQ/IMAGINE/run103.saved/2006N080737_LEBMI_S24_L001_R3_001.fastq.gz";
#my $R3 = $dirin."/".$cp->{R2};#"/data-isilon/data/sequences/ILLUMINA/NOVASEQ/IMAGINE/run103.saved/2006N080737_LEBMI_S24_L001_R2_001.fastq.gz";

open(my $fR1,"zcat $R1 |");
open(my $fR2,"zcat $R2 |");

my @desc = ($fR1,$fR2);
my @types = ("name","seq","plus","qual");
my @line2 = (141,"*",0,0,"*","*",0,0);
my @line1 = (77,"*",0,0,"*","*",0,0);
my $line1 = join("\t",@line1);
my $line2 = join("\t",@line2);


my $nb_paired;

my $nb =0;
my $t1 = time;
my $finish;


while ( ! eof($desc[0]) ) {
 	my $name  = undef;
 	my $z=0;
 	my $umi;
	for (my $i = 0;$i<2;$i++){
		my $x;
		my $y;
		my $xl = readline($desc[$i]);
		($y,$x) = split(" ",$xl);
		
		substr($y, 0, 1) = "";
		
		$name = $y unless $name;
		die() if $name ne $y; 
	} 
	my @seq;
	for (my $i = 0;$i<2;$i++){
		 $seq[$i] = readline($desc[$i]);
		 chomp($seq[$i]);
	} 
	
	for (my $i = 0;$i<2;$i++){
		 readline($desc[$i]);
	} 
	
	my @qual;
	for (my $i = 0;$i<2;$i++){
		 $qual[$i] = readline($desc[$i]);
		  chomp($qual[$i]);
	} 

	my @t1 = split(":",$name);
	my $umi = $t1[-1];
	print $obam join("\t",($name,$line1,$seq[0],$qual[0],"RG:Z:$rg\tRX:Z:".$umi))."\n";
	print $obam join("\t",($name,$line2,$seq[1],$qual[1],"RG:Z:$rg\tRX:Z:".$umi))."\n";
	print $ofastq1 join("\n","@".$name,$seq[0],"+",$qual[0])."\n";
	print $ofastq2 join("\n","@".$name,$seq[1],"+",$qual[1])."\n";
	
	#print  join("\t",($name,$line1,$seq[0],$qual[0],"RG:Z:$rg\tRX:Z:".$seq[2]))."\n";
	#print  join("\t",($name,$line2,$name,$seq[1],$qual[1],"RG:Z:$rg\tRX:Z:".$seq[2]))."\n";
	
 	
 	$nb++;
 	if ($nb % 1_000_000 == 0 ){
 		warn "elapsed ".abs($t1-time)." ".$nb;
 	}
	if ($nb % 5_000_000 == 0 ){
	
		close($obam);
		close($ofastq1);
		close($ofastq2);
		#push(@files,{fastq1=>$fastq1,fastq2=>$fastq2,ubam=>$ubam});
		$nb_file ++ ;
		$fastq1 = $generic.".$nb_file.$fq"."_R1.fastq.gz";
		warn $fastq1;
		$fastq2 = $generic.".$nb_file.$fq"."_R2.fastq.gz";
		$ubam = $generic.".$nb_file.$fq".".u.bam";
		push(@files,{fastq1=>$fastq1,fastq2=>$fastq2,ubam=>$ubam});
		open( $ofastq1,"| bgzip > $fastq1");
		open( $ofastq2,"| bgzip > $fastq2");
		open($obam,"|-","java -Djava.io.tmpdir=/tmp -jar /software/distrib/picard/2.23.3/picard.jar SortSam -I /dev/stdin  -VERBOSITY ERROR -SORT_ORDER queryname -O $ubam -QUIET true 2>/dev/null");
		print $obam qq{\@HD	VN:1.6	SO:queryname	GO:none
\@RG	ID:$rg	SM:$rg	LB:$rg	PL:illumina
};
		
		
		warn "elapsed ".abs($t1-time)." ".$nb;
		$t1 =time;
		

		
	}
 }
 
close($obam);
close($ofastq1);
close($ofastq2);
 $pm->finish(0,{files=>\@files});
}

	$pm->wait_all_children();
	
	
	warn $json_file;
	
	open (JSON,">$json_file");
	print JSON encode_json \@global_files;
	close JSON;
	exit(0);
	
