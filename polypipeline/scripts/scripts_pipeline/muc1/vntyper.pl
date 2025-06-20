#!/usr/bin/perl
use FindBin qw($Bin);
use strict;

use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages/";
#use Set::IntSpan;
use GBuffer; 
use Data::Dumper;
use Getopt::Long;
use Carp;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use String::ProgressBar;
use List::Util qw(sum);
use JSON::XS;
use file_util;


 my $project_name;
 my $fork;
 my $callable_intspan_file;
 my $patient_name2;
 #my $low_calling;
 my $version;
 my $method;
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"patient=s" => \$patient_name2,
	"version=s" => \$version,
);


my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name,-version=>$version );

if ($patient_name2 =~ / /){
	my $query = $buffer->getQuery();
	my $res   = $query->getPatients( $project->id );
	my ($a,$b) = split(" ",$patient_name2);
	my (@find) = grep {$_->{name} =~/$a/ && $_->{name} =~ /$b/} @$res;
	die() if scalar(@find) ne 1;
	$patient_name2 = $find[0]->{name};
	warn $find[0]->{name};
}

#exit(0) if -e "/data-isilon/sequencing/muc1/renome/".$patient_name2.".vcf";
my $patients = $project->get_only_list_patients($patient_name2);
die() unless scalar(@$patients);
warn "coucou ".scalar(@$patients);
 my $dir_pipeline = "/tmp/pipeline/$project_name.".time."/";#

 system("mkdir -p $dir_pipeline") unless -e $dir_pipeline;
foreach my $patient (@$patients) {
my $patient_name = $patient->name;
my $bam = $patient->getBamFile();
  my $tmp_dir = $project->getCallingPipelineDir("muc1-vntyper.".$patient->name);
  $dir_pipeline = $tmp_dir;
my $fileoutx=  $project->getVariationsDir("muc1_vntyper6")."/".$patient_name.".vcf";
 my $dir_pipeline2 = $dir_pipeline ."/$patient_name/";
 system("rm -r $dir_pipeline2") if -e $dir_pipeline2;

 #exit(0) if -e $fileout;

 my $dir_working = "$dir_pipeline/$patient_name"  ;
  if (-e $dir_working){
 	system("rm $dir_pipeline/*.log*");
 }
 else {
 system ("mkdir $dir_working");
 }
my $singularity = $buffer->software("singularity");
#my $image = "/software/distrib/ADVNTR/SINGLARITY/vntyper_v7_20221209.sif";
my $image = "/software/distrib/ADVNTR/SINGLARITY/vntyper.sif";
my $db = "/data-isilon/public-data/repository/HG19/vntr/";
system ("mkdir $dir_pipeline/temp") unless -e "$dir_pipeline/temp";

my $release = 'hg19';
$release = 'hg38' if $project->getVersion() =~ /HG38/;

my $cmd = "$singularity run --pwd /opt/vntyper -B /data-isilon:/data-isilon -B /data-beegfs/:/data-beegfs /data-beegfs/software/sif/vntyper_main.sif vntyper   pipeline --bam $bam -o $tmp_dir   --reference-assembly $release"; 
warn $cmd;
#my $cmd = qq{$singularity run --pwd /DATA/adVNTR/ -B /data-isilon:/data-isilon -B /tmp/:/tmp/ -H $tmp_dir  $image $scommand};

#system("samtools index $bam -\@2");
my $fad = $tmp_dir."/advntr/output_adVNTR.vcf" ;
my $fad2 = $tmp_dir."/advntr/output_adVNTR_result.tsv" ;
 my $fk = $tmp_dir."/kestrel/kestrel_result.tsv";
system($cmd." && touch $dir_pipeline/$patient_name.ok") unless -e $fad2;


unless (-e "$dir_pipeline/$patient_name.ok"){
 #	system("rm -r $dir_pipeline");
 	die();
 }
 
 my $dir_prod = $project->getVariationsDir("vntyper");

  die() unless -e $fk;
  my $kestrel = parse_tsv($fk);
  returnkestrel_hash($kestrel);

  #die() unless -e $fad2;
  
# my $fad2 = $tmp_dir."/advntr/output_adVNTR_result.tsv" ;
#  my $advntr = parse_adVNTR($fad2,$fad);
 my $json_text = encode_json {kestrel=>$kestrel};
 my $file_json = $dir_prod."/".$patient->name.".json";
open(my $fh, ">", $file_json) or die "Impossible d'ouvrir le fichier $file_json: $!";
	print $fh $json_text;
	close($fh);
	warn $file_json;
}
 exit(0);

sub returnkestrel_hash{
	my($hash) = @_;
	return if $hash->{data}->{Confidence} eq "Negative";
	my $ref = $hash->{data}->{REF} ;
	my $alt = $hash->{data}->{ALT};
	my $motif = $hash->{data}->{Motif_sequence};
	
	my $l = length ($motif);
	my $pos = $hash->{data}->{POS};
	my $position = ($l - $pos);
	$hash->{data}->{Motif_sequence_origin} =$hash->{data}->{Motif_sequence};
	if (length($alt)> length($ref)){
			#insertion 
			my $seq_alt = substr($motif, 0, $pos)."XXX".substr($motif, $pos);
			my $reverse_alt = BioTools::complement_sequence($seq_alt);
			 my $ralt = BioTools::complement_sequence(substr($alt, 1));
			my $ins = qq{<span style="color:#0F79B2;">[<span style="text-emphasis: double-circle #0F79B2; ">$ralt</span/>]</span/>};
			#  my $ins = qq{<span style="color:blue;text-emphasis: double-circle blue; ">$ralt</span/>};
			$reverse_alt =~ s/XXX/$ins/;
			$hash->{data}->{Motif_sequence_origin} =$hash->{data}->{Motif_sequence};
			 $hash->{data}->{Motif_sequence} = $reverse_alt;
			$hash->{data}->{Motif_sequence} = $reverse_alt;
		}
		elsif (length($alt)< length($ref)){
			$pos ++;
			 $ref  = substr($ref, 1);
			my $end = $pos + length($ref);
			my $seq_alt = $motif;
			my $lref = length($motif);
			if ($end > $lref ){
				my $add = substr($ref,($lref-$pos)+1);
				$seq_alt .= lc($add);
			}
			 my $seq_alt2 = substr($seq_alt, 0, $pos-1)."XXX".substr($seq_alt, ($pos+length($ref)-1));
			my $alt = substr($seq_alt,$pos-1,length($ref));
			my $ralt = BioTools::complement_sequence($alt);
			my $ins = qq{<span style="color:red;">[<span style="text-emphasis: double-circle red; ">$ralt</span/>]</span/>};
			$seq_alt2 =  BioTools::complement_sequence($seq_alt2);
			$seq_alt2 =~ s/XXX/$ins/;
			$hash->{data}->{Motif_sequence}  =$seq_alt2;
		}
		elsif (length($alt) == length($ref)){
			my $reverse = BioTools::complement_sequence($hash->{data}->{Motif_sequence});
			my $ralt = BioTools::complement_sequence($alt);
			 my $ins = qq{<span style="color:red">[$ref/$ralt]</span/>};
			 my $left = substr($reverse, 0, $position);
			my $right = substr($reverse, $position+2);
			$hash->{data}->{Motif_sequence}  =$left.$ins.$right;
		}
		
	
	
}

sub parse_adVNTR {
	my ($file,$file2) = @_;
	my $res;
	my $date = POSIX::strftime( 
             "%d/%m/%y", 
             localtime( 
               		(stat $file2 )[10]
                 )
             );
	unless( -e $file){
		$res->{header} = ["Caller","date","Confidence"];
		push( @{$res->{data}},["adVNTR",$date,"Negative"]);
		return $res;
	}
	my @lines = `grep -v "#" $file`;  
	chomp(@lines);
	
	$res->{header} = ["date","State","NumberOfSupportingReads","MeanCoverage","Pvalue"];
	
	my @hs = split(" ",shift @lines);
	$res->{header} = ["caller","date",@hs];
	foreach my $l (@lines) {
		my @tt = split(" ",$l);
		my @aa = split("&",$tt[1]);
		my @bb = split("_",$aa[0]);
		my $repeat = $bb[1];
		$tt[0] = "Insertion" ;
		$tt[0] = "Deletetion" if $aa[0] =~ /^D/;
		push( @{$res->{data}},["adVNTR",$date,$repeat,$tt[0],$tt[1],"-",$tt[2],$tt[3],$tt[4]] ) ;
	}
	push(@{$res->{data}},["adVNTR",$date]) unless @lines; 
	return $res;
}



sub parse_tsv {
	my ($fk) =@_;
	open (FH,"$fk");
 my @lines;
 while (<FH>) {
 	next if $_ =~/^#/;
 	chomp($_);
 	push(@lines,$_);
 }
 close (FH);
 my @header = split(" ",shift @lines);
 
 my @data;
 my $res;
  $res->{date} = POSIX::strftime( 
             "%d/%m/%y", 
             localtime( 
               		(stat $fk )[10]
                 )
             );
 foreach my $l (@lines) {
 	my @dataa = split(" ",$l);
 	
 	for (my $i =0;$i<@header;$i++){
 		$res->{$header[$i]} = $dataa[$i];
 	}
 	push(@data,$res);
 	last;
 }
 my $t;
 $t->{header} = \@header;
  $t->{data} = $res;
 return $t;
}
