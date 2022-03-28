#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";
use GBuffer;
use Storable qw/freeze thaw/;
use JSON;

my $cgi = new CGI();
my $project_name	= $cgi->param('project');
my $gene_id			= $cgi->param('gene');
my $chr_name		= $cgi->param('chromosome');
my $var_ids			= $cgi->param('ids');
my $used_filters	= $cgi->param('used_filters');
my $dejavu			= $cgi->param('dejavu');
my $dejavu_ho		= $cgi->param('dejavu_ho');

my @lIds = split(',', $var_ids);
my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $project_name );

my $hash;
if ($gene_id eq 'all' and $chr_name eq 'all') {
	my $h_tmp;
	foreach my $ids (split('::', $var_ids)) {
		my @lTmp = split(':', $ids);
		foreach my $id (split(',', $lTmp[-1])) {
			$h_tmp->{$lTmp[0]}->{$id} = undef;
		}
	}
	my $nb = 0;
	foreach my $chr_id (keys %$h_tmp) {
		$nb += scalar keys %{$h_tmp->{$chr_id}};
	}
	$hash->{'id'}				= 'ALL';
	$hash->{'complete_name'}	= 'ALL';
	$hash->{'position'} 		= 'ALL';
	$hash->{'strand'} 			= 'Forward';
	$hash->{'description'} 		= 'ALL';
	$hash->{'variations'} 		= $nb;
}
elsif ($gene_id =~ /intergenic/) {
	my ($inter, $chr_name, $start, $end) = split('_', $gene_id);
	$hash->{'id'}				= $gene_id;
	$hash->{'complete_name'}	= $gene_id;
	$hash->{'position'} 		= $start.'Reverse'.$end;
	$hash->{'strand'} 			= 'Forward';
	$hash->{'description'} 		= "Intergenic from $start to $end";
	$hash->{'variations'} 		= scalar(@lIds);
}
else {
	my $chr = $project->getChromosome( $chr_name );
	my $gene;
	if ($gene_id =~ /ENSG/) { $gene = $chr->getGene( $gene_id ); }
	else {
		foreach my $g (@{$chr->getGenes()}) {
			if ($g->external_name() eq $gene_id) {
				$gene = $g;
				last;
			}
		}
	}
	$hash->{'id'} = $gene_id;
	$hash->{'complete_name'} = '['.$gene_id.'] '.$gene->external_name();
	$hash->{'position'}	= $gene->start().'-'.$gene->end();
	if ($gene->strand() eq '1') { $hash->{'strand'} = 'Forward'; }
	else { $hash->{'strand'} = 'Reverse'; }
	my $description = $gene->description();
	$description =~ s/ \[.+\]//;
	$hash->{'description'} = $description;
	$hash->{'variations'} = scalar(@lIds);
}

my @lPatNames;
foreach my $patient (@{$project->getPatients()}) {
	push(@lPatNames, $patient->name());	
}

my @lPat;
foreach my $name (sort @lPatNames) {
	my $patient = $project->getPatient($name);
	my $align = join(',', @{$patient->alignmentMethods()});
	#my $calling = join(',', @{$patient->getCallingMethods()});
	push(@lPat, $name.':'.$align);
	#push(@lPat, $name.':'.$align.':'.$calling);
}
$hash->{'patients'} = join(' ', @lPat);

my (@lFilters, $hFilters, $hFilters2, $hFreq);
foreach my $filter (keys %{$buffer->config->{ensembl_annotations}}) {
	$hFilters->{$filter} = undef;
}
foreach my $filter (sort split(' ', $used_filters)) {
	next unless ($filter);
	if (exists $hFilters->{$filter}) {
		delete $hFilters->{$filter};
		if ($filter eq 'coding' or $filter eq 'nonsynonymous' or $filter eq 'non-synonymous') {
			delete $hFilters->{'coding'} if (exists $hFilters->{'coding'});
			delete $hFilters->{'nonsynonymous'} if (exists $hFilters->{'nonsynonymous'});
			delete $hFilters->{'non-synonymous'} if (exists $hFilters->{'non-synonymous'});
		}
	}
	elsif ($filter =~ /freq/) { $hFreq->{$filter} = undef; }
}

foreach my $filter (sort keys %$hFilters) {
	my @lTmp = split(';', $buffer->config->{ensembl_annotations}->{$filter});
	$hFilters2->{$lTmp[-1]} = undef;
}

foreach my $impact (keys %{$project->impacts_ensembl_annotations()}) {
	my $is_ok = 1;
	my @lTmp;
	foreach my $annot (keys %{$project->impacts_ensembl_annotations->{$impact}}) {
		if (exists $hFilters2->{$annot}) { push(@lTmp, $annot); }
		else { $is_ok = undef; }
	}
	if ($is_ok) {
		if ($impact eq 'intergenic') { push(@lFilters, 'Intergenic'); }
		else { push(@lFilters, 'ALL '.$impact.' Impact'); }
	}
	else {
		push(@lFilters, @lTmp);
	}
}

$hash->{'used_filters'} = join(', ', sort @lFilters);
if (exists $hFreq->{'freq_001'})   { $hash->{'freq'} = '<1&#8241;'; }
elsif (exists $hFreq->{'freq_01'}) { $hash->{'freq'} = '<1&#8240;'; }
elsif (exists $hFreq->{'freq_05'}) { $hash->{'freq'} = '<1&#37;'; }
elsif (exists $hFreq->{'freq_1'})  { $hash->{'freq'} = '<5&#37;'; }
else { $hash->{'freq'} = 'ALL'; }

$hash->{'dejavu'} = 'No';
if ($dejavu) {
	$hash->{'dejavu'} = $dejavu;
	if ($dejavu_ho) {
		$hash->{'dejavu'} .= ' (Ho)';
	}
}
	
my @list;
push(@list, $hash);
my $hashInfos;
$hashInfos->{label} = 'name';
$hashInfos->{items} = \@list;

print $cgi->header('text/json-comment-filtered');
print encode_json $hashInfos;
exit(0);