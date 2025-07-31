#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/kyoto";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";

#use lib "/bip-d/soft/distrib/tabix/latest/perl";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/coverage";
use lib "$Bin/../packages/validation_variation";
require
  "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";
use html;

#use Set::;
use Carp;
use export_data;
use strict;
use Set::IntSpan::Fast::XS;
use Data::Dumper;
use GBuffer;
use Getopt::Long;
use Bio::DB::HTS;







my $buffer = GBuffer->new();

my $fsize        = "font-size:10px";
my $cgi          = new CGI();
my $project_name = $cgi->param('project');

my $project = $buffer->newProjectCache(
	-name        => $project_name,
	-typeFilters => 'individual',
	-cgi_object  => 1
);

my $patient_name = $cgi->param('patient');

my $patient = $project->getPatient($patient_name);
my $bam = $patient->getAlignmentFile();
my $bam_obj = Bio::DB::HTS->new(-bam => $bam,-fasta => $project->genomeFasta,);


my $l1 = qq{chr21	46283735	46283735	G	A	rs843345
chr3	184188727	184188727	T	C	rs1058018
chr17	48922889	48922889	C	T	rs8017
chr16	2771572	2771572	C	T	rs3738494
chr1	42659188	42659188	C	T	rs1065483
chr17	5381475	5381475	G	A	rs2839181
chr21	46266025	46266025	A	G	rs11059924
chr12	128808801	128808801	C	T	rs2075144
chr19	46354029	46354029	G	A	rs6795772
chr3	49327836	49327836	C	T	rs456261
chr6	33290666	33290666	G	A	rs1131620
chr19	40611963	40611963	A	G	rs2231926
chr3	73062658	73062658	A	G	rs352169
chr2	105038258	105038258	C	T	rs3739160
};
my $l1 = qq{chr21	46283735	46283735	G	A	rs843345
chr3	184188727	184188727	T	C	rs1058018
chr17	48922889	48922889	C	T	rs8017
chr16	2771572	2771572	C	T	rs3738494
chr1	42659188	42659188	C	T	rs1065483
chr17	5381475	5381475	G	A	rs2839181
chr21	46266025	46266025	A	G	rs11059924
chr12	128808801	128808801	C	T	rs2075144
chr19	46354029	46354029	G	A	rs6795772
chr3	49327836	49327836	C	T	rs456261
chr6	33290666	33290666	G	A	rs1131620
chr19	40611963	40611963	A	G	rs2231926
chr3	73062658	73062658	A	G	rs352169
chr3	52202746	52202746	G	A	rs3739160
chr2	105038258	105038258	C	T	rs3739161
} if $project->genome_version_generic =~/HG38/;




my @snps;

foreach my $l (split("\n",$l1)){
	chomp($l);
	my @z = split(" ",$l);
	$z[0]=~ s/chr//;
	push(@snps,{chr=>$z[0],start=>$z[1],ref=>$z[3],alt=>$z[4],rs=>$z[5]});
}
my @iv = split( "", $patient->identity_vigilance() );
my @iv_vcf = split( "", $patient->identity_vigilance_vcf() );
$iv_vcf[5] = 4;

print $cgi->header('text/html; charset=utf-8');

print <<'HTML';
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>Tableau des SNPs</title>
  <style>
    table { border-collapse: collapse; margin-top: 20px; }
    th, td { border: 1px solid #999; padding: 5px 10px; text-align: center; }
    th { background-color: #DADaff; }
    caption { font-weight: bold; margin-bottom: 10px; }
  </style>
</head>
<body>
<h1>Profil SNP</h1>
<table>
  <caption>Identito Vigilence snpflex</caption>
  <tr>
HTML

print $patient->identity_vigilance()." ->". $patient->identity_vigilance_vcf()."<br>";

my %error;
for (my$i=0;$i<@iv;$i++){
	$error{$i} ++ if $iv[$i] ne $iv_vcf[$i];
}
my $i =0;
foreach my $snp (@snps) {
	my $n = $snp->{rs} ." ".$snp->{chr}.":".$snp->{start};
	my $style = qq{};
	if (exists $error{$i}){
		$style = qq{style="background-color:tomato";};
	}
    print "<th $style >$n</th>";
    $i++;
}
print "</tr>\n<tr>";

for (my $i=0;$i<@snps;$i++){
	my $n = return_text($snps[$i],$iv[$i]);
	my$n2 = return_text($snps[$i],$iv_vcf[$i]);
	my $style = qq{style="background-color:#f8f8f8";};
	if (exists $error{$i}){
		$style = qq{style="background-color:tomato";};
	}
	print "<td $style>$n-$n2</td>";
}
print "</tr>\n<tr>";


for (my $i=0;$i<@snps;$i++){
	my $h = return_pileup($snps[$i],$bam_obj);
		my $v1 = $h->{A};
		my $v2 = $h->{T};
		my $v3 = $h->{C};
		my $v4 = $h->{G};
		
	my $n = "A($v1) T ($v2) C($v3) G($v4)";
	my $style = qq{style="background-color:white";};
	if (exists $error{$i}){
		$style = qq{style="background-color:tomato";};
	}
	print "<td $style >$n</td>";
}

print "</tr>\n</table>";


#for (my $i=0;$i<@snps;$i++){
#	#if ($iv[$i] == $iv_vcf[$i]){
#	#	print return_text($snps[$i],$iv[$i])."\n";
#	#}
#	#else {
#		my $h = return_pileup($snps[$i],$bam_obj);
#		my $t = return_text($snps[$i],$iv_vcf[$i]);
#	
#		my @a = split("",$t);
#		my $v1 = $h->{A};
#		my $v2 = $h->{T};
#		my $v3 = $h->{C};
#		my $v4 = $h->{G};
#		print $snps[$i]->{rs}."::".return_text($snps[$i],$iv[$i])." <=> A($v1) T ($v2) C($v3) G($v4) \n";
#	#}
#}

sub return_pileup {
	my ($snp,$bam_obj) =@_;
	my $hash = {};
		$hash->{'A'} = 0;
		$hash->{'T'} = 0;
		$hash->{'G'} = 0;
		$hash->{'C'} = 0;
		$hash->{'INS'} = 0;
		$hash->{'DEL'} = 0;
		my $chr = $project->getChromosome($snp->{chr});
		my $start = $snp->{start};
		#ATTENTION: le start presente un decalage d une base (+1) apr rapport au BED. Samtools + bed = demarage a 0 donc decalage d une base
		my $end = $snp->{start}+1;
		my $region = $chr->fasta_name.":".($snp->{start}-1)."-".$end; 
		$bam_obj->pileup($region, sub {
    		my ($seqid, $pos, $pileups) = @_;
			return if $pos != $start;
    		foreach my $p (@$pileups) {
        	next if $p->is_refskip;  # ignorer les sauts de type 'N' (RNA-seq)

       		if ($p->is_del) {
            	$hash->{'DEL'}++;
        	} elsif ($p->indel > 0) {
            	$hash->{'INS'}++;
        	} else {
        		my $aln = $p->alignment;
				my $qpos = $p->qpos;
				my $base = substr($aln->qseq, $qpos, 1);
            	$hash->{$base}++ if $base =~ /^[ATCG]$/;
        	}
    		}
		});
		return $hash;
}

sub return_text {
	my ($s,$c) = @_;
	return $s->{ref}. $s->{ref} if ($c == 1);
	return $s->{ref}. $s->{alt} if ($c == 3);
	return $s->{alt}. $s->{alt} if ($c == 2);
	return "NC";
}
