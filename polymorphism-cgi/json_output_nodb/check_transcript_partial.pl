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

#NGS2014_0566
#UNC93B1
#ENST00000227471_11

die() unless $project_name;
die() unless $enst_name;


my $h_var;
my $h_cores_nt;
$h_cores_nt->{'A'} = 'C';
$h_cores_nt->{'T'} = 'G';
$h_cores_nt->{'G'} = 'A';
$h_cores_nt->{'C'} = 'A';

my $buffer  = new GBuffer;


#my $project = $buffer->newProject( -name => 'NGS2014_0437', -verbose => 1 );
my $project = $buffer->newProject( -name => $project_name, -verbose => 1 );
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
	$h_exons_introns->{v19}->{$eiid} = Set::IntSpan::Fast->new();
	$h_exons_introns->{v19}->{$eiid}->add_range($e->start(), $e->end());
}
foreach my $i (@{$t->getIntrons()}) {
	my $eiid = $i->id();
	$eiid =~ s/$tid//;
	$h_exons_introns->{v19}->{$eiid} = Set::IntSpan::Fast->new();
	$h_exons_introns->{v19}->{$eiid}->add_range($i->start(), $i->end());
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
my $seq = $t->getChromosome->sequence($start, $end);
my @lSeq = split('', $seq);
my $i = 0;
foreach my $nt (@lSeq) {
	my $id = 'chr'.$t->getChromosome->id().'_'.($start + $i).'_'.$nt.'_'.$h_cores_nt->{$nt};
	if ($intspan_pos->contains($start + $i)) { 
		$h_var->{$id} = undef;
	}
	$i++;
}

my $j = 0;
foreach my $varid (keys %$h_var) {
	$j++;
	print '.' if $j % 25 == 0;
	my $var = $project->_newVariant($varid);
#	$h_var->{$varid}->{v19}->{annot} = $var->variationTypeInterface();
	$h_var->{$varid}->{v19}->{annot} = $var->variationTypeInterface($t);
	eval {
		if ( $var->isCoding($t) ) {
			my $prot = $t->getProtein();
			if ($prot) {
				$h_var->{$varid}->{v19}->{'nomenclature'} = $var->getNomenclature($t);
				my $cds_pos = $var->getOrfPosition($prot);
				$cds_pos = '-' if (not $cds_pos or $cds_pos eq '.');
				$h_var->{$varid}->{v19}->{'cds_position'} = $cds_pos;
				$h_var->{$varid}->{v19}->{'protein_position'} = $var->getProteinPosition($prot);
				my $protAA = $var->getProteinAA($prot);
				my $chanAA = $var->changeAA($prot);
				if ( $protAA and $chanAA ) {
					$h_var->{$varid}->{v19}->{'aa'} = $protAA . '/' . $chanAA;
				}
				$h_var->{$varid}->{v19}->{'coding_sequence'} = $t->coding_sequence();
			}
		}
	};
	if ($@) {
		$h_var->{$varid}->{v19}->{'nomenclature'} = '-';
		$h_var->{$varid}->{v19}->{'cds_position'} = '-';
		$h_var->{$varid}->{v19}->{'protein_position'} = '-';
		$h_var->{$varid}->{v19}->{'aa'} = '-';
	}
	$h_var->{$varid}->{v19}->{'exon_intron'} = '-';
	foreach my $eiid (keys %{$h_exons_introns->{v19}}) {
		if ($h_exons_introns->{v19}->{$eiid}->contains($var->start())) {
			$h_var->{$varid}->{v19}->{'exon_intron'} = $eiid;
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
$project->changeAnnotationVersion('43.18', 1);
$t = $project->newTranscript($enst_name);

foreach my $e (@{$t->getExons()}) {
	my $eiid = $e->id();
	$eiid =~ s/$tid//;
	$h_exons_introns->{v43}->{$eiid} = Set::IntSpan::Fast->new();
	$h_exons_introns->{v43}->{$eiid}->add_range($e->start(), $e->end());
}
foreach my $i (@{$t->getIntrons()}) {
	my $eiid = $i->id();
	$eiid =~ s/$tid//;
	$h_exons_introns->{v43}->{$eiid} = Set::IntSpan::Fast->new();
	$h_exons_introns->{v43}->{$eiid}->add_range($i->start(), $i->end());
}

foreach my $varid (sort keys %$h_var) {
	$j++;
	print '.' if $j % 25 == 0;
	my $var = $project->_newVariant($varid);
	eval {
		#$h_var->{$varid}->{v43}->{annot} = $var->variationTypeInterface();
		$h_var->{$varid}->{v43}->{annot} = $var->variationTypeInterface($t);
	};
	if ($@) {
		$h_var->{$varid}->{v43}->{annot} = 'error';
	}
	eval {
		if ( $var->isCoding($t) ) {
			my $prot = $t->getProtein();
			if ($prot) {
				$h_var->{$varid}->{v43}->{'nomenclature'} = $var->getNomenclature($t);
				my $cds_pos = $var->getOrfPosition($prot);
				$cds_pos = '-' if (not $cds_pos or $cds_pos eq '.');
				$h_var->{$varid}->{v43}->{'cds_position'} = $cds_pos;
				$h_var->{$varid}->{v43}->{'protein_position'} = $var->getProteinPosition($prot);
				my $protAA = $var->getProteinAA($prot);
				my $chanAA = $var->changeAA($prot);
				if ( $protAA and $chanAA ) {
					$h_var->{$varid}->{v43}->{'aa'} = $protAA . '/' . $chanAA;
				}
				$h_var->{$varid}->{v43}->{'coding_sequence'} = $t->coding_sequence();
			}
		}
	};
	if ($@) {
		$h_var->{$varid}->{v43}->{'nomenclature'} = $@;
		$h_var->{$varid}->{v43}->{'cds_position'} = '-';
		$h_var->{$varid}->{v43}->{'protein_position'} = '-';
		$h_var->{$varid}->{v43}->{'aa'} = '-';
	}
	$h_var->{$varid}->{v43}->{'exon_intron'} = '-';
	foreach my $eiid (keys %{$h_exons_introns->{v43}}) {
		if ($h_exons_introns->{v43}->{$eiid}->contains($var->start())) {
			$h_var->{$varid}->{v43}->{'exon_intron'} = $eiid;
			last;
		}
	}
	
	#print ' : '.$h_var->{$varid}->{v19}->{annot}.' / '.$h_var->{$varid}->{v43}->{annot}."\n";
	#die if $h_var->{$varid}->{v19} ne $h_var->{$varid}->{v43};
}


my $hash;
my $i = 0;
foreach my $varid (sort keys %$h_var) {
	$i++;
	my $line;
	
	$line .=  $varid;
	$line .=  ';'.$h_var->{$varid}->{v19}->{annot};
	$line .=  ';'.$h_var->{$varid}->{v19}->{exon_intron};
	
	$line .=  ';';
	if (exists $h_var->{$varid}->{v19}->{'nomenclature'}) { $line .=  $h_var->{$varid}->{v19}->{'nomenclature'}; }
	else { $line .=  '-'; }
	
	$line .=  ';';
	if (exists $h_var->{$varid}->{v19}->{'cds_position'}) { $line .=  $h_var->{$varid}->{v19}->{'cds_position'}; }
	else { $line .=  '-'; }
	
	$line .=  ';';
	if (exists $h_var->{$varid}->{v19}->{'protein_position'}) { $line .=  $h_var->{$varid}->{v19}->{'protein_position'}; }
	else { $line .=  '-'; }
	
	$line .=  ';';
	if (exists $h_var->{$varid}->{v19}->{'aa'}) { $line .=  $h_var->{$varid}->{v19}->{'aa'}; }
	else { $line .=  '-'; }
	
	
	$line .=  ';'.$h_var->{$varid}->{v43}->{annot};
	$line .=  ';'.$h_var->{$varid}->{v43}->{exon_intron};
	
	$line .=  ';';
	if (exists $h_var->{$varid}->{v43}->{'nomenclature'}) { $line .=  $h_var->{$varid}->{v43}->{'nomenclature'}; }
	else { $line .=  '-'; }
	
	$line .=  ';';
	if (exists $h_var->{$varid}->{v43}->{'cds_position'}) { $line .=  $h_var->{$varid}->{v43}->{'cds_position'}; }
	else { $line .=  '-'; }
	
	$line .=  ';';
	if (exists $h_var->{$varid}->{v43}->{'protein_position'}) { $line .=  $h_var->{$varid}->{v43}->{'protein_position'}; }
	else { $line .=  '-'; }
	
	$line .=  ';';
	if (exists $h_var->{$varid}->{v43}->{'aa'}) { $line .=  $h_var->{$varid}->{v43}->{'aa'}; }
	else { $line .=  '-'; }
	$hash->{$i} = $line;
}

my $json_encode = encode_json $hash;
print ".\",";
$json_encode =~ s/{//;
print $json_encode;
exit(0);

