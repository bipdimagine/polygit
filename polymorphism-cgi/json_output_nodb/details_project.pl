#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";
use GBuffer;
use Digest::MD5 qw(md5_hex);
use JSON;

my $cgi = new CGI();
my $projectName	= $cgi->param('project');
my $type	    = $cgi->param('type');
my $user	    = $cgi->param('user');

# analyse  global_infos  check->vcf->variations/indels


my $buffer = new GBuffer;
my $project = $buffer->newProjectCache( -name 			=> $projectName,
										-cache 			=> '1', );
										
my @lItems;
my $isGenome;

if ($type eq 'annot_infos') {
	@lItems = @{$project->json_ensembl_annotations()};
	my $h;
	$h->{id} = 'can_use_hgmd';
	if ($project->hasHgmdAccess($user)) { $h->{value} = 'YES'; }
	else { $h->{value} = 'NO'; }
	push(@lItems, $h);
	
	$h = undef;
	$h->{id} = 'is_cache_updated';
	if ($project->isUpdate()) { $h->{value} = 'YES'; }
	else { $h->{value} = 'NO'; }
	push(@lItems, $h);
	
	$h = undef;
	$h->{id} = 'annotation_version';
	my @lAnnots;
	
	foreach my $annot_version (@{$project->annotation_version_history()}) {
		next unless ($project->isCacheDone($annot_version));
		push(@lAnnots, $annot_version);
	}
#	push(@lAnnots, $project->annotation_version());
	$h->{value} = join(',', @lAnnots);
#	if ($project->annotation_version_current() ne $project->annotation_version()) {
#		$h->{value} .= ' <span style="color:red"><b> -> new release '.$project->annotation_version_current().' available ! </b></span>';
#		$h->{value} .= ' '.$project->annotation_version_current();
#	}
	
	push(@lItems, $h);
	
	foreach my $annot_version (@{$project->annotation_version_history()}) {
		next unless ($annot_version =~ /\./);
		$h = undef;
		my $version_hgmd = $project->buffer->get_public_data_version("hgmd",$annot_version);
		my $version_clinvar = $project->buffer->get_public_data_version("clinvar",$annot_version);
		my $version_cosmic = $project->buffer->get_public_data_version("cosmic",$annot_version);
		$h->{id} = 'annotation_version_details_'.$annot_version;
		my @lValues;
		my @lTmp = split('\.', $annot_version);
		foreach my $cat (sort keys %{$project->buffer->public_data->{$lTmp[1]}}) {

			next if ($cat eq 'cldb');
			next if ($cat eq 'cytoband');
			next if ($cat eq 'gnomad-genome');
			next if ($cat eq 'gnomad-rsname');
			my $cat_name = $cat;
			$cat_name = 'gnomad' if ($cat eq 'gnomad-exome');
			my $version = $lTmp[1];
			if (int($version) > 0) {
				my $previous_version = int($version) - 1;
				unless ($project->buffer->public_data->{$previous_version}->{$cat}->{'version'} eq $project->buffer->public_data->{$version}->{$cat}->{'version'}) {
					push(@lValues, "<span style='color:orange;'><b>$cat_name </b>".$project->buffer->public_data->{$version}->{$cat}->{'version'}.'</span>');
				}
				else {
					push(@lValues, "<b>$cat_name </b>".$project->buffer->public_data->{$version}->{$cat}->{'version'});
				}
			}
			else {
				push(@lValues, "<b>$cat_name </b>".$project->buffer->public_data->{$version}->{$cat}->{'version'});
			}
		}
		$h->{value} = "<b>genecode </b> ".$lTmp[0].', '.join(', ', @lValues);
		push(@lItems, $h);
	}
}
else {
	my $hInfos = $project->global_infos->{$type};
	my $i = 0;
	foreach my $key (@{listOrderKeyToCheckFromType($type, $hInfos)}) {
		if ($key eq 'genome') {
			unless (exists $hInfos->{genome}) {
				my $hash;
				$hash->{line} = $i;
				$hash->{key} = uc($key);
				$hash->{value} = 'NO';
				push(@lItems, $hash);
				$i++;
				next;
			}
		}
		my $value = $hInfos->{$key};
		if (ref($value) eq 'ARRAY') {
			if (scalar @$value > 1) {
				my $j = 1;
				foreach my $val (sort @$value) {
					my $hash;
					$hash->{line} = $i;
					$hash->{key} = uc($key).'_'.$j;
					$hash->{value} = $val;
					$i++;
					$j++;
					push(@lItems, $hash);
				}
			}
			else {
				my $hash;
				$hash->{line} = $i;
				$hash->{key} = uc($key);
				$hash->{value} = $value->[0];
				$i++;
				push(@lItems, $hash);
			}
		}
		elsif (ref($value) eq 'HASH') {
			foreach my $second_key (sort keys %$value) {
				my $hash;
				my $second_value = $value->{$second_key};
				$hash->{line} = $i;
				$hash->{key} = uc($key).' '.uc($second_key);
				$hash->{value} = $second_value;
				$i++;
				push(@lItems, $hash);
			}
		}
		else {
			my $hash;
			$hash->{line} = $i;
			$hash->{key} = uc($key);
			if ($key eq 'diagnostic' or $key eq 'exome' or $key eq 'genome') {
				if ($value) {
					$hash->{value} = 'YES';
					$isGenome = 1 if ($key eq 'genome');
				}
				else { $hash->{value} = 'NO'; }
			}
			else { $hash->{value} = $value;	}
			$i++;
			push(@lItems, $hash);
		}
	}
	foreach my $key (keys %{$project->maskImpactText()}) {
		my $hash;
		my $value = $project->maskImpactText->{$key};
		$hash->{key} = 'impact_'.$key;
		$hash->{value} = join(', ', @$value);
		push(@lItems, $hash);
	}
	my $hash;
	$hash->{key} = 'IS_BIG_DATA';
	$hash->{value} = '0';
	if ($isGenome) { $hash->{value} = '1'; }
	elsif (scalar @{$project->getPatients()} > 80) { $hash->{value} = '1'; }
	elsif (scalar @{$project->getPatients()} > 50) { $hash->{value} = '2'; }
	my @lChr = @{$project->getChromosomes()};
	if ($lChr[0]->getVariantsVector->Size() > 500000) { $hash->{value} = '4'; }
	push(@lItems, $hash);
	
	my $h1;
	$h1->{key} = 'ho_regions_short';
	$h1->{value} = $project->ho_regions_short_value();
	push(@lItems, $h1);
	my $h2;
	$h2->{key} = 'ho_regions_medium';
	$h2->{value} = $project->ho_regions_medium_value();
	push(@lItems, $h2);
	my $h3;
	$h3->{key} = 'ho_regions_large';
	$h3->{value} = $project->ho_regions_large_value();
	push(@lItems, $h3);
	my $h4;
	$h4->{key} = 'ho_regions_extra_large';
	$h4->{value} = $project->ho_regions_extra_large_value();
	push(@lItems, $h4);
}

my $hashRes;
$hashRes->{'label'} = 'line';
$hashRes->{'items'} = \@lItems;
print $cgi->header('text/json-comment-filtered');
print encode_json $hashRes;
exit(0);


sub listOrderKeyToCheckFromType {
	my ($type, $hash) = @_;
	my @list;
	if ($type eq 'global_infos') { @list = ('creation_date', 'name', 'description', 'project_type'); }
	elsif ($type eq 'analyse') { @list = ('build', 'genome', 'exome', 'diagnostic', 'capture', 'version', 'cache', 'alignment', 'calling'); }
	else { @list = sort keys %$hash; }
	return \@list;
}
