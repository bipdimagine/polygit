#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/kyoto";
use lib "$Bin/../packages/export";
use Carp;
use export_data;
use strict;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use JSON::XS;
use validationQuery;
my $buffer = GBuffer->new();
my $cgi    = new CGI();
my $project_name = $cgi->param('project');
my $genes_names = $cgi->param('genes');
my @genes = split(",",$genes_names);
my $project = $buffer->newProject(-name=>$project_name);
my @out;


my $hTransFromCaptures;
my @lObjs;
foreach my $genes_name (sort @genes){
	if ($genes_name =~ /capture/) {
		my @lTmp = split(':', $genes_name);
		my $capture_name = $lTmp[-1];
		foreach my $parent_capture_name (keys %{$buffer->getAllGenesNamesInAllBundle()}) {
			next unless (exists $buffer->getAllGenesNamesInAllBundle->{$parent_capture_name}->{$capture_name});
			foreach my $tr_id (@{$buffer->getAllGenesNamesInAllBundle->{$parent_capture_name}->{$capture_name}->{transcripts}}) {
				unless (exists $hTransFromCaptures->{$tr_id}) {
					push(@lObjs, $tr_id);
					$hTransFromCaptures->{$tr_id} = undef;
				}
			}
		}
	}
	else {
		push(@lObjs, $genes_name);
	}
}

foreach my $genes_name (sort @lObjs){
	my $z = $project->liteAnnotations->get_like("annotations","$genes_name");
	next unless defined $z;
	foreach my $k (keys %$z){
 		my ($n,$g,$gene_name,$id) = split(" ",$k);
		my $hgene;
 		unless ( $gene_name =~/ENSG/){
			$genes_name = $z->{$k}->{gene_kyoto_id};
			my $zz = $project->liteAnnotations->get("annotations","$genes_name");
			my @tr;
			my ($nt,$chr) = split("_",$z->{$k}->{genbo_id});
			foreach my $t (@{$zz->{transcripts}}){
				next if $t ne $nt;
				push (@tr,$t);
			}
			$zz->{transcripts} = \@tr;
			$hgene = $zz;
		}
		else{
			$hgene = $z->{$k};
		}
		unless ($hgene){
			my $item;
			$item->{name} = "not found";
			$item->{gene} = $gene_name;
			push(@out,$item);
			next;
		}
		my $found_default;
		my @out1;
		foreach my $t (@{$hgene->{transcripts}}){
			unless ($t =~ /_/) {
				$t = $t."_".$hgene->{chromosome};
			}
			my $tr = $project->liteObject($t);
			my (@exons)  = split(",",$tr->{genomic_span}->as_string());
			my $item;
			$item->{name} = $tr->{id};
			$item->{transcript} = $t;
			$item->{'chr'} = $tr->{chromosome};
			$item->{default} =0;
			$item->{default} = 1 if $tr->{ccds_name};
			$found_default = 1 if $tr->{ccds_name};
			$item->{ccds} = $tr->{ccds_name};
			$item->{ccds} =~ s/\..*//;
			$item->{refseq} = $tr->{external_name};
			$item->{refseq} =~ s/\..*//;
			$item->{protein} = $tr->{external_protein_name};
			$item->{start} = $tr->{start};
			$item->{end} = $tr->{end};
			$item->{gene} = $hgene->{external_name};
			$item->{label} = $hgene->{external_name};
			$item->{description} = $hgene->{description};
			$item->{description} =~ s/\[.*\]//;
			$item->{bundle} = "";#$hgene->{external_name};
			$item->{nb_exons} = "".scalar(@exons)."";
			push(@out1,$item);
		}
		@out1 = sort{$a->{nb_exons} <=> $b->{nb_exons}} @out1;
		$out1[-1]->{default} =1;
		push(@out,@out1);
	}
}	

export_data::print_json($cgi,\@out);
exit(0);

