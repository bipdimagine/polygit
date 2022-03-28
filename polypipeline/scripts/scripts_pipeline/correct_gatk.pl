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
my $filein;
my $dir;
my $file_bed;
my $project_name;
my $chr_name;
my $fork = 1;
use List::Util qw(max min sum);
use List::MoreUtils qw(firstidx);
use calling;

my $DP_MIN_FOR_CORRECTION = 3;

my $AC_CUT   = 5;
my $LIMIT_AC = 3;
my $LIMIT_PC = 0.25;

my $vcf_file;
my $dir_out;
my $log_file;
my $project_name;
my $patient_name;
my $snp_out;
my $indels_out;
my $method;

GetOptions(
	"vcf=s"          => \$vcf_file,
	"dir=s"          => \$dir_out,
	"snp_out=s" => \$snp_out,
	"indels_out=s" => \$indels_out,
	"log=s"          => \$log_file,
	"project=s"      => \$project_name,
	"patient_name=s" => \$patient_name,
	"method=s"=> \$method,
);

my $buffer  = GBuffer->new();

my $project = $buffer->newProject( -name => $project_name );
my $error   = 0;
my %hpatients;
map { $hpatients{ $_->name } = 0 }
  @{ $project->get_list_patients($patient_name) };

my $patients = $project->get_list_patients($patient_name);
my %hbarcode;
foreach my $p (@$patients){
	 $hbarcode{$p->barcode} = $p->name();
}

open( VCF, " bcftools view $vcf_file | " );
#my $snp_out = $dir_out . "/" . $project->name . ".snps.vcf";
#warn $snp_out;
my $dir_out =  $project->getCallingPipelineDir($method);
my $output1  = calling::getTmpFile($dir_out,"correct","vcf");
open( OUT, ">" . $output1 );

#open( OUT_SNP, ">" . $snp_out );
#my $indels_out = $dir_out . "/" . $project->name . ".indels.vcf";
#open( OUT_INDEL, ">" . $indels_out );

if ($log_file) {
	open( STDOUT, ">>" . $log_file );
}

die($vcf_file) unless -e $vcf_file;

print colored ['black ON_BRIGHT_MAGENTA'],
  "======= Correct and filter GATK ==== ";
print color 'reset';
print "\n";

my $format = 0;
my $nline;
my $cline;
my $ref = 0;
my $alt = 0;
my $info;
my $filter;
my $dbsnp_before;
my $dbsnp_after;
my $total_snp;

while ( my $line = <VCF> ) {
	if ( $line =~ /^##/ ) {
		print OUT $line;
		next;
	}
	chomp($line);
	if ( $line =~ /^#CHROM/ ) {
		my @array = split( " ", $line );

		for ( my $i = 0 ; $i < @array ; $i++ ) {
			$format = $i if $array[$i] eq "FORMAT";
			$ref    = $i if $array[$i] eq "REF";
			$alt    = $i if $array[$i] eq "ALT";
			$info   = $i if $array[$i] eq "INFO";
			if ($format) {
				my $n = $array[$i];
				if (exists $hbarcode{$n}){
					$array[$i] = $hbarcode{$n};
				}
				$hpatients{ $array[$i] }++;

			}
		}
		my @toto = grep { $_ ne 1 } values %hpatients;
		if ( scalar(@toto) > 1 ) {
			warn Dumper %hpatients;
			error( "problem in vcf with sample  " . join( ":", @toto ) );
			
		}
		print OUT join("\t",@array) . "\n";
		if ( $format == 0 || $ref == 0 || $alt == 0 || $info == 0 ) {
			error("problem vcf format $line ");

		}
		next;
	}    #end if CHROM
	my @array = split( " ", $line );
	my $debug;
	$debug = 1 if $array[1] == 130914965;
	warn $line if $debug;
	my @infos =   split(";",$array[$info]);
	my $set_idx =  firstidx { $_ =~/set/ } @infos;
	if ($infos[$set_idx] =~/free/ || $infos[$set_idx] =~/Inter/){
		
			print OUT $line."\n";
			next;
	}
	my $indel     = is_indel( join( ",", ( $array[$ref], $array[$alt] ) ) );
	my $hinfo     = parse_info("$array[$info]");
	my @format_gt = split( ":", $array[$format] );
	my $gt        = firstidx { $_ eq "GT" } @format_gt;
	my $ad        = firstidx { $_ eq "AD" } @format_gt;

	die("GT NO FIRST POSITION ".$line ) if $gt == -1;
	die() unless $ad;
	my @newline;
	for ( my $i = 0 ; $i < $format + 1 ; $i++ ) {
		push( @newline, $array[$i] );
	}
	my $change;
	my $max_all = 0;
	my $max_dp;

	for ( my $ii = $format + 1 ; $ii < @array ; $ii++ ) {
		my @array_sample = split( ":", $array[$ii] );
		my $sample_gt = $array_sample[$gt];
		my $sample_ad = $array_sample[$ad];
		$sample_gt = "./." if ($sample_ad eq ".");
		if ( $sample_gt eq "./." ) {
			push( @newline, "./." );
			next;
		}
		my $sample_ad = $array_sample[$ad];
		
		my @geno = split( "/", $sample_gt );

		my @array_ad = split( ",", $sample_ad );
		my $max_index = -1;
		my $max       = -999;
		my $max_dp_local = 0;
		for ( my $i = 1 ; $i < @array_ad ; $i++ ) {
			if ( $i > 0 && $max_all < $array_ad[$i] ) {
				$max_all = $array_ad[$i];
				$max_dp  = sum @array_ad;
			}
			if ( $max < $array_ad[$i] ) {
				$max_index = $i;
				$max       = $array_ad[$i];
				$max_dp_local = sum @array_ad;
			}    #end if max
		}
		die($line) if $max_index == -1;
		my $pc = -1;
		my $dp = sum @array_ad;
		if ( $max > 0 || $array_ad[0] > 0 ) {
			$pc = int( $max / ($dp) * 100 );

		}    #end if
		my $limit_ho     = 90;
		my $limit_he_inf = 15;
		my $limit_he_sup = 60;
		#########################
		### CORRECT WEIRD CALLING
		##########################
	
		if ( $max_dp_local > $DP_MIN_FOR_CORRECTION ) {
			my $string = $array_sample[$gt];
			if ( $pc > $limit_ho && $geno[1] ne $geno[0] ) {
				$array_sample[$gt] = "$max_index/$max_index";
			}
			if (   $pc >= $limit_he_inf
				&& $pc <= $limit_he_sup
				&& $geno[1] eq $geno[0] )
			{

				$array_sample[$gt] = "0/$max_index";
			}

		}    #end if
		
		push( @newline, join( ":", @array_sample ) );
	}
	$total_snp ++ unless $indel ;
	next if $max_dp ==0;
	my $max_pc = ($max_all/$max_dp);
	unless ($indel){
		my $st = join("\t",@newline);
		$dbsnp_before ++ if $st =~ /rs/;
	}
	warn "coucou " if $debug;
	warn $max_all." ".$max_pc if $debug;
	
	
	#if ($max_all >= $LIMIT_AC || ($max_all >= 10 && $max_pc > 0.3)) {
	if ($max_all >= $LIMIT_AC || ($max_all >= 5 && $max_pc > 0.25)) {
		die("1") if $debug;
			if ($indel){
				#--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \ 
			#	if ($hinfo->{QD} < 2 || $hinfo->{FS} > 200.0 || $hinfo->{ReadPosRankSum} < -20 ){
					#$filter ++;
				#	next;
			#	}
				print OUT join("\t",@newline)."\n";
			}
			else {
				
				#"QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || ReadPosRankSum < -8.0"
				$hinfo->{ReadPosRankSum} =0 unless exists $hinfo->{ReadPosRankSum};
				$hinfo->{MappingQualityRankSum} = 0 unless exists $hinfo->{MappingQualityRankSum};
			#	if ($hinfo->{QD} < 2 || $hinfo->{FS} > 60.0 || $hinfo->{MQ} < 40.0 || $hinfo->{HaplotypeScore} > 13.0  || $hinfo->{ReadPosRankSum} < -8.0){# || $hinfo->{MappingQualityRankSum} < 12.5 ){
					#|| $hinfo->{HaplotypeScore} > 13.0 || $hinfo->{ReadPosRankSum} < -8.0 || $hinfo->{MappingQualityRankSum} < 12.5 ) {
				#		$filter ++;
				#		next;
				#	}
				print OUT join("\t",@newline)."\n";
			}
		unless ($indel){
			my $st = join("\t",@newline);
			$dbsnp_after ++ if $st =~ /rs/;
		}
			next;
	}
	

	if ($max_pc < 0.25  && $max_all <= 4){
		die("2") if $debug;
			$filter ++;
				next;
	}
	my $filter_text ="";
	if ($indel){
		if ($hinfo->{QD} < 2){
			$filter_text .="QD2;";
		}
		if ($hinfo->{QUALITY} < 30){
			$filter_text .="QUALITY;";
		}
		if ($hinfo->{FS} > 200){
			$filter_text .="FS;";
		}
		
		if ($filter_text){
			$filter ++;
			next;
		}
		
		print OUT join("\t",@newline)."\n";# if $change;
	}
	else{
#		die("3") if $debug;
#		#AC=1;AF=0.167;AN=6;BaseQRankSum=0.519;DP=110;ExcessHet=3.0103;FS=1.527;MLEAC=1;MLEAF=0.167;
#		#MQ=60;MQRankSum=0;QD=5.63;
#		#ReadPosRankSum=0.839;SOR=1.061	GT:AD:DP:GQ:PL	0/1:23,9:32:99:189,0,580
		my $filter_text ="";
		if ($hinfo->{QD} < 2){
			$filter_text .="QD;";
		}
		if ($hinfo->{MQ} < 40){
			$filter_text .="MQ;";
		}
		if ($hinfo->{FS} > 60){
			$filter_text .="FS;";
		}
		if ($hinfo->{MQRankSum} < -12.5){
			$filter_text .="MQR;";
		}
		if ( $hinfo->{ReadPosRankSum} < -8.0){
			$filter_text .="RPR;";
		}
		if ($filter_text){
			$filter ++;
			next;
		}
		print OUT join("\t",@newline)."\n";
	}

	unless ($indel){
		my $st = join("\t",@newline);
		$dbsnp_after ++ if $st =~ /rs/;
	}
}
warn "$filter :: $total_snp :: $dbsnp_before :: $dbsnp_after";
close(OUT);
close(VCF);


warn "$dir_out/$project_name.final.vcf $output1";
system("mv $output1 $dir_out/$project_name.final.vcf");
warn "ok";
exit(0);
#open( OUT_SNP, ">" . $snp_out );
#my $indels_out = $dir_out . "/" . $project->name . ".indels.vcf";
#open( OUT_INDEL, ">" . $indels_out );

exit(0);

sub parse_info {
	my ($st) = @_;
	my %data;
	foreach my $s ( split( ";", $st ) ) {
		my ( $k, $v ) = split( "=", $s );
		$data{$k} = $v;

	}
	return \%data;
}

sub error {
	my ($st) = @_;
	print colored ['black ON_BRIGHT_RED'], $st;
	print color 'reset';
	print "\n";
	exit(1);
}

sub is_indel {
	my ($st) = @_;
	foreach my $s ( split( ",", $st ) ) {
		return 1 if length($s) > 1;
	}
	return;
}
