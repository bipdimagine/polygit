#!/usr/bin/perl
use FindBin qw($Bin $RealBin);
use strict;

use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
use GBuffer ;
use Data::Dumper;
use Getopt::Long;
use Bio::DB::HTS;
use Bio::DB::HTS::VCF ;
use List::Util qw(max);
use Text::CSV qw(csv);
use List::MoreUtils qw(firstidx);


my $project_name;
my $fork;
my $patient_names;
my $version;
my $exons;
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"patients=s" => \$patient_names,
	"version=s" => \$version,
	"genes" => sub { $exons = 0 },
	"exons" => \$exons,
) || die ("Error in command line arguments\n");

die("-project mandatory") unless ($project_name);
die("One of -genes or -exons is mandatory") unless (defined $exons);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name , -version =>$version);
my $patients =  $project->getPatients();
#my $patients = $project->get_list_patients($patient_names);

my $dir_out= $project->getCountingDir("featureCounts");
my $fileout = $dir_out."/".$project_name.".count.genes.txt" unless ($exons);
$fileout = $dir_out."/".$project_name.".count.exons.txt" if ($exons);
	
	
	my @bams;
	my @sed_cmd;
	my $align_method;
	my $profile;
	my %strands;
	foreach my $patient (@$patients){
		my $bam = $patient->getBamFile;
		my $name = $patient->name;
		my $metrics_file = $project->getCountingDir("featureCounts") . "/metrics/$name.metrics";
		die("Can't find '$metrics_file'") unless (-e $metrics_file);
		push(@bams,$bam);
		$bam =~ s/\//\\\//g;
		push(@sed_cmd,qq{sed -i "2s/$bam/$name/" $fileout} );
		my $aoa = csv (in => $metrics_file, sep => "\t");
		my $metrics;
		for my $metric (qw {PF_ALIGNED_BASES PF_BASES PCT_CODING_BASES PCT_MRNA_BASES PCT_USABLE_BASES PCT_R1_TRANSCRIPT_STRAND_READS PCT_R2_TRANSCRIPT_STRAND_READS}) {
			my $index = firstidx{ $_ eq $metric } @{$aoa->[6]};
			$metrics->{$metric} = $aoa->[7]->[$index];
			die("ERROR parsing '$metrics_file': no '$metric' found: ".join("\t",$aoa->[6])) if ($index == -1);
		}
		#warn $name.": ". Dumper $metrics;
		die ("ERROR $name aligned bases too low: ".sprintf("%.2f%%",$metrics->{PF_ALIGNED_BASES}/$metrics->{PF_BASES} *100)
			."\nCheck release (current = ".$project->getVersion().") or contaminations") if ($metrics->{PF_ALIGNED_BASES}/$metrics->{PF_BASES} < 0.4);
		die ("ERROR $name mRNA bases too low: ".sprintf("%.2f%%",$metrics->{PCT_MRNA_BASES} *100)) if ($metrics->{PCT_MRNA_BASES} < 0.2); # (UTR_BASES + CODING_BASES)/PF_ALIGNED_BASES
		die ("ERROR $name usable bases too low: ".sprintf("%.2f%%",$metrics->{PCT_USABLE_BASES} *100)) if ($metrics->{PCT_USABLE_BASES} < 0.2); # (CODING_BASES + UTR_BASES)/PF_BASES
		warn ("$name aligned bases low: ".sprintf("%.2f%%",$metrics->{PF_ALIGNED_BASES}/$metrics->{PF_BASES} *100)) and sleep(1) if ($metrics->{PF_ALIGNED_BASES}/$metrics->{PF_BASES} < 0.6);
		warn ("$name mRNA bases low: ".sprintf("%.2f%%",$metrics->{PCT_MRNA_BASES} *100)) if ($metrics->{PCT_MRNA_BASES} < 0.5); # (UTR_BASES + CODING_BASES)/PF_ALIGNED_BASES
		warn ("$name usable bases low: ".sprintf("%.2f%%",$metrics->{PCT_USABLE_BASES} *100)) if ($metrics->{PCT_USABLE_BASES} < 0.5); # (CODING_BASES + UTR_BASES)/PF_BASES
		my $pct_r1 = $metrics->{PCT_R1_TRANSCRIPT_STRAND_READS};
		my $pct_r2 = $metrics->{PCT_R2_TRANSCRIPT_STRAND_READS};
		die("ERROR pct R1 and R2 transcript strand reads are anormal for '$name': R1 = $pct_r1\tR2 = $pct_r2\n$metrics_file") if ($pct_r1 + $pct_r2 != 1);
		warn "$name\tR1 = $pct_r1";
		$strands{'-s 1 '} ++ if ($pct_r1 >= 0.9);
		$strands{'-s 2 '} ++ if ($pct_r1 <= 0.1);
		$strands{'-s 0 '} ++ if ($pct_r1 >= 0.4 and $pct_r1 <= 0.6);
		$strands{'error'}->{"$name"} = $pct_r1 if (($pct_r1 > 0.1 and $pct_r1 < 0.4) or ($pct_r1 > 0.6 and $pct_r1 < 0.9));
		$align_method = $patient->alignmentMethod();
		$profile = $patient->getSampleProfile();
	}
	warn 'Strands'.Dumper \%strands;
	my @strands = keys %strands;
	die("Error: pct R1 transcript strand reads:\n".Dumper \%strands) if (grep{/error/} @strands);
	die("More than one strand for the ".scalar @$patients." patients in project $project_name:\n".Dumper \%strands) unless (scalar @strands);

	my $ppn = 16;
	my $gtf = $project->gtf_file();
	
	my $featureCounts = $project->buffer->software("featureCounts");
	my $strand = " -s 1 ";
	$strand = " -s 2 " if $profile eq "bulk illumina pcr-free" or $profile eq "bulk ribozero pcr-free" or $profile eq "bulk NEB-directional pcr-free" or $profile eq "bulk watchmaker pcr-free";
	$strand = " -s 0 " if $profile eq "bulk neb pcr-free" ;
	die("Strands from metrics (".$strands[0].") and from profile ($strand) don't match") unless ($strand eq $strands[0]);
	my $strand = $strands[0];
	
	my $sed = join(" && ",@sed_cmd);
	my $cmd = "$featureCounts -T $ppn -a $gtf --ignoreDup -o $fileout -p -t exon $strand ";
	$cmd .= "-f -O " if ($exons);
	$cmd .= join(" ",@bams)." && $sed";
#	my $cmd1 = "$featureCounts -T $ppn -a $gtf --ignoreDup -o $fileout -p -t exon  $strand ".join(" ",@bams)." && $sed";
#	my $cmd2 = "$featureCounts -T $ppn -a $gtf --ignoreDup -o $fileout2 -p -t exon -f -O $strand ".join(" ",@bams)." && $sed2";
	my $result = system($cmd);
	die() unless $result ne 0;

