#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;

use strict;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";
use GBuffer;
use JSON;
use Set::IntSpan::Fast;


my $cgi = new CGI();
my $project_name	= $cgi->param('project');
my $enst_name		= $cgi->param('transcript');
my $out				= $cgi->param('out');
my $check_only_position = $cgi->param('only_position');
my $gencode_v1		= $cgi->param('gencode_v1');
my $gencode_v2		= $cgi->param('gencode_v2');
my $skip_intronic	= $cgi->param('skip_intronic');
my $skip_utr		= $cgi->param('skip_utr');
my $skip_pseudogene	= $cgi->param('skip_pseudogene');
my $type_simulation	= $cgi->param('type');

#NGS2014_0566
#UNC93B1
#ENST00000227471_11

die() unless $project_name;
die() unless $enst_name;

$gencode_v1 = '19' unless ($gencode_v1);
$gencode_v2 = '43' unless ($gencode_v2);
my $gencode_v1_text = 'v'.$gencode_v1;
my $gencode_v2_text = 'v'.$gencode_v2;

my $h_var;
my $h_cores_nt;
if ($type_simulation eq 'ins') {
	$h_cores_nt->{'A'}->{'CGTCGA'} = undef;
	$h_cores_nt->{'A'}->{'CGTCGAC'} = undef;
	$h_cores_nt->{'T'}->{'GGCTACGAT'} = undef;
	$h_cores_nt->{'T'}->{'GGCTACGATC'} = undef;
	$h_cores_nt->{'G'}->{'AAACTAGC'} = undef;
	$h_cores_nt->{'G'}->{'AAACTAGCC'} = undef;
	$h_cores_nt->{'C'}->{'ACGACGCATCGTAGCTAC'} = undef;
	$h_cores_nt->{'C'}->{'ACGACGCATCGTAGCTACG'} = undef;
}
elsif ($type_simulation eq 'del') {
	$h_cores_nt->{'A'}->{'ACGTCGA'} = undef;
	$h_cores_nt->{'A'}->{'ACGTCGAC'} = undef;
	$h_cores_nt->{'T'}->{'TGGCTACGAT'} = undef;
	$h_cores_nt->{'T'}->{'TGGCTACGATC'} = undef;
	$h_cores_nt->{'G'}->{'GAAACTAGC'} = undef;
	$h_cores_nt->{'G'}->{'GAAACTAGCC'} = undef;
	$h_cores_nt->{'C'}->{'CACGACGCATCGTAGCTAC'} = undef;
	$h_cores_nt->{'C'}->{'CACGACGCATCGTAGCTACG'} = undef;
}
else {
	$h_cores_nt->{'A'}->{'T'} = undef;
	$h_cores_nt->{'A'}->{'G'} = undef;
	$h_cores_nt->{'A'}->{'C'} = undef;
	$h_cores_nt->{'T'}->{'A'} = undef;
	$h_cores_nt->{'T'}->{'G'} = undef;
	$h_cores_nt->{'T'}->{'C'} = undef;
	$h_cores_nt->{'G'}->{'A'} = undef;
	$h_cores_nt->{'G'}->{'T'} = undef;
	$h_cores_nt->{'G'}->{'C'} = undef;
	$h_cores_nt->{'C'}->{'A'} = undef;
	$h_cores_nt->{'C'}->{'T'} = undef;
	$h_cores_nt->{'C'}->{'G'} = undef;
}

my $buffer  = new GBuffer;


#my $project = $buffer->newProject( -name => 'NGS2014_0437', -verbose => 1 );
my $project = $buffer->newProject( -name => $project_name, -verbose => 1 );
$project->changeAnnotationVersion($gencode_v1.'.18', 1);
#my $gene = $project->newGene($gene_name);

print $cgi->header('text/json-comment-filtered');
print "{\"progress\":\".";

my $intspan_pos = Set::IntSpan::Fast->new();
my $t = $project->newTranscript($enst_name);

my $tid = $t->id();

my $h_exons_introns;
foreach my $e (@{$t->getExons()}) {
	my $eiid = $e->id();
	$eiid =~ s/$tid//;
	$h_exons_introns->{$gencode_v1_text}->{$eiid} = Set::IntSpan::Fast->new();
	$h_exons_introns->{$gencode_v1_text}->{$eiid}->add_range($e->start(), $e->end());
}
foreach my $i (@{$t->getIntrons()}) {
	my $eiid = $i->id();
	$eiid =~ s/$tid//;
	$h_exons_introns->{$gencode_v1_text}->{$eiid} = Set::IntSpan::Fast->new();
	$h_exons_introns->{$gencode_v1_text}->{$eiid}->add_range($i->start(), $i->end());
}

my $seq = $t->getChromosome->sequence($t->start, $t->end);
my @lSeq = split('', $seq);
if ($check_only_position) {
	if ($check_only_position =~ '-') {
		my @lTmp = split('-', $check_only_position);
		$intspan_pos->add_range($lTmp[0], $lTmp[1]);
	}
	else {
		$intspan_pos->add($check_only_position);
	}
}
else {
	$intspan_pos->add_range($t->start(), $t->end());
}
my $start = $t->start();
my $end = $t->end();
my $i = 0;
foreach my $nt (@lSeq) {
	my ($ref_all, $alt_all);
	if ($intspan_pos->contains($start + $i)) {
		foreach my $nt2 (sort keys %{$h_cores_nt->{$nt}}) {
			if ($type_simulation eq 'del') {
				$h_var->{'chr'.$t->getChromosome->id().'_'.($start + $i).'_'.$nt2.'_'.$nt} = undef;
			}
			else {
				$h_var->{'chr'.$t->getChromosome->id().'_'.($start + $i).'_'.$nt.'_'.$nt2} = undef;
			}
		}
	}
	$i++;
}

my $j = 0;
foreach my $varid (keys %$h_var) {
	$j++;
	print '.' if $j % 25 == 0;
	my $var = $project->_newVariant($varid);
#	$h_var->{$varid}->{$gencode_v1_text}->{annot} = $var->variationTypeInterface();
	$h_var->{$varid}->{$gencode_v1_text}->{annot} = $var->variationTypeInterface($t);
	eval {
		if ( $var->isCoding($t) ) {
			my $prot = $t->getProtein();
			if ($prot) {
				$h_var->{$varid}->{$gencode_v1_text}->{'nomenclature'} = $var->getNomenclature($t);
				my $cds_pos = $var->getOrfPosition($prot);
				$cds_pos = '-' if (not $cds_pos or $cds_pos eq '.');
				$h_var->{$varid}->{$gencode_v1_text}->{'cds_position'} = $cds_pos;
				$h_var->{$varid}->{$gencode_v1_text}->{'protein_position'} = $var->getProteinPosition($prot);
				my $protAA = $var->getProteinAA($prot);
				my $chanAA = $var->changeAA($prot);
				if ( $protAA and $chanAA ) {
					$h_var->{$varid}->{$gencode_v1_text}->{'aa'} = $protAA . '/' . $chanAA;
				}
				$h_var->{$varid}->{$gencode_v1_text}->{'coding_sequence'} = $t->coding_sequence();
			}
		}
	};
	if ($@) {
		$h_var->{$varid}->{$gencode_v1_text}->{'nomenclature'} = '-';
		$h_var->{$varid}->{$gencode_v1_text}->{'cds_position'} = '-';
		$h_var->{$varid}->{$gencode_v1_text}->{'protein_position'} = '-';
		$h_var->{$varid}->{$gencode_v1_text}->{'aa'} = '-';
	}
	$h_var->{$varid}->{$gencode_v1_text}->{'exon_intron'} = '-';
	foreach my $eiid (keys %{$h_exons_introns->{$gencode_v1_text}}) {
		if ($h_exons_introns->{$gencode_v1_text}->{$eiid}->contains($var->start())) {
			$h_var->{$varid}->{$gencode_v1_text}->{'exon_intron'} = $eiid;
			last;
		}
	}
}



$project = undef;
$buffer = undef;
$t = undef;



$buffer  = new GBuffer;
#my $project = $buffer->newProject( -name => 'NGS2014_0437', -verbose => 1 );
$project = $buffer->newProject( -name => 'NGS2014_0566', -verbose => 1 );
$project->changeAnnotationVersion($gencode_v2.'.18', 1);
$t = $project->newTranscript($enst_name);

foreach my $e (@{$t->getExons()}) {
	my $eiid = $e->id();
	$eiid =~ s/$tid//;
	$h_exons_introns->{$gencode_v2_text}->{$eiid} = Set::IntSpan::Fast->new();
	$h_exons_introns->{$gencode_v2_text}->{$eiid}->add_range($e->start(), $e->end());
}
foreach my $i (@{$t->getIntrons()}) {
	my $eiid = $i->id();
	$eiid =~ s/$tid//;
	$h_exons_introns->{$gencode_v2_text}->{$eiid} = Set::IntSpan::Fast->new();
	$h_exons_introns->{$gencode_v2_text}->{$eiid}->add_range($i->start(), $i->end());
}

foreach my $varid (sort keys %$h_var) {
	$j++;
	print '.' if $j % 25 == 0;
	my $var = $project->_newVariant($varid);
	eval {
		#$h_var->{$varid}->{$gencode_v2_text}->{annot} = $var->variationTypeInterface();
		$h_var->{$varid}->{$gencode_v2_text}->{annot} = $var->variationTypeInterface($t);
	};
	if ($@) {
		$h_var->{$varid}->{$gencode_v2_text}->{annot} = 'error';
	}
	eval {
		if ( $var->isCoding($t) ) {
			my $prot = $t->getProtein();
			if ($prot) {
				$h_var->{$varid}->{$gencode_v2_text}->{'nomenclature'} = $var->getNomenclature($t);
				my $cds_pos = $var->getOrfPosition($prot);
				$cds_pos = '-' if (not $cds_pos or $cds_pos eq '.');
				$h_var->{$varid}->{$gencode_v2_text}->{'cds_position'} = $cds_pos;
				$h_var->{$varid}->{$gencode_v2_text}->{'protein_position'} = $var->getProteinPosition($prot);
				my $protAA = $var->getProteinAA($prot);
				my $chanAA = $var->changeAA($prot);
				if ( $protAA and $chanAA ) {
					$h_var->{$varid}->{$gencode_v2_text}->{'aa'} = $protAA . '/' . $chanAA;
				}
				$h_var->{$varid}->{$gencode_v2_text}->{'coding_sequence'} = $t->coding_sequence();
			}
		}
	};
	if ($@) {
		$h_var->{$varid}->{$gencode_v2_text}->{'nomenclature'} = $@;
		$h_var->{$varid}->{$gencode_v2_text}->{'cds_position'} = '-';
		$h_var->{$varid}->{$gencode_v2_text}->{'protein_position'} = '-';
		$h_var->{$varid}->{$gencode_v2_text}->{'aa'} = '-';
	}
	$h_var->{$varid}->{$gencode_v2_text}->{'exon_intron'} = '-';
	foreach my $eiid (keys %{$h_exons_introns->{$gencode_v2_text}}) {
		if ($h_exons_introns->{$gencode_v2_text}->{$eiid}->contains($var->start())) {
			$h_var->{$varid}->{$gencode_v2_text}->{'exon_intron'} = $eiid;
			last;
		}
	}
	
	#print ' : '.$h_var->{$varid}->{$gencode_v1_text}->{annot}.' / '.$h_var->{$varid}->{$gencode_v2_text}->{annot}."\n";
	#die if $h_var->{$varid}->{$gencode_v1_text} ne $h_var->{$varid}->{$gencode_v2_text};
}


my $hash;
$i = 0;
foreach my $varid (sort keys %$h_var) {
	
	my $annot1 = $h_var->{$varid}->{$gencode_v1_text}->{annot};
	my $annot2 = $h_var->{$varid}->{$gencode_v1_text}->{annot};
	if ($skip_intronic) {
		next if ($annot1 eq $annot2 and lc($annot1) eq 'intronic');
	}
	if ($skip_utr) {
		next if ($annot1 eq $annot2 and lc($annot1) eq 'utr');
	}
	if ($skip_pseudogene) {
		next if ($annot1 eq $annot2 and lc($annot1) eq 'pseudogene');
	}
	
	$i++;
	my $line;
	
	$line .=  $varid;
	$line .=  ';'.$h_var->{$varid}->{$gencode_v1_text}->{annot};
	$line .=  ';'.$h_var->{$varid}->{$gencode_v1_text}->{exon_intron};
	
	$line .=  ';';
	if (exists $h_var->{$varid}->{$gencode_v1_text}->{'nomenclature'}) { $line .=  $h_var->{$varid}->{$gencode_v1_text}->{'nomenclature'}; }
	else { $line .=  '-'; }
	
	$line .=  ';';
	if (exists $h_var->{$varid}->{$gencode_v1_text}->{'cds_position'}) { $line .=  $h_var->{$varid}->{$gencode_v1_text}->{'cds_position'}; }
	else { $line .=  '-'; }
	
	$line .=  ';';
	if (exists $h_var->{$varid}->{$gencode_v1_text}->{'protein_position'}) { $line .=  $h_var->{$varid}->{$gencode_v1_text}->{'protein_position'}; }
	else { $line .=  '-'; }
	
	$line .=  ';';
	if (exists $h_var->{$varid}->{$gencode_v1_text}->{'aa'}) { $line .=  $h_var->{$varid}->{$gencode_v1_text}->{'aa'}; }
	else { $line .=  '-'; }
	
	
	$line .=  ';'.$h_var->{$varid}->{$gencode_v2_text}->{annot};
	$line .=  ';'.$h_var->{$varid}->{$gencode_v2_text}->{exon_intron};
	
	$line .=  ';';
	if (exists $h_var->{$varid}->{$gencode_v2_text}->{'nomenclature'}) { $line .=  $h_var->{$varid}->{$gencode_v2_text}->{'nomenclature'}; }
	else { $line .=  '-'; }
	
	$line .=  ';';
	if (exists $h_var->{$varid}->{$gencode_v2_text}->{'cds_position'}) { $line .=  $h_var->{$varid}->{$gencode_v2_text}->{'cds_position'}; }
	else { $line .=  '-'; }
	
	$line .=  ';';
	if (exists $h_var->{$varid}->{$gencode_v2_text}->{'protein_position'}) { $line .=  $h_var->{$varid}->{$gencode_v2_text}->{'protein_position'}; }
	else { $line .=  '-'; }
	
	$line .=  ';';
	if (exists $h_var->{$varid}->{$gencode_v2_text}->{'aa'}) { $line .=  $h_var->{$varid}->{$gencode_v2_text}->{'aa'}; }
	else { $line .=  '-'; }
	$hash->{$i} = $line;
}

my $json_encode = encode_json $hash;
print ".\",";
$json_encode =~ s/{//;
print $json_encode;
exit(0);

