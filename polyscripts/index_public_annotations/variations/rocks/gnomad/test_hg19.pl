#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
#use lib "/data-isilon/bipd-src/pnitschk/git2-polyweb/polygit/GenBo/lib/obj-nodb";
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
use GBuffer;
use Set::IntSpan::Fast::XS;
use Getopt::Long;
use GenBoNoSqlRocksAnnotation;
use Sereal qw(sereal_encode_with_object sereal_decode_with_object write_sereal_file read_sereal_file);
my $fork = 1;
my $chr = "21";

my $dir_out = "/data-isilon/public-data/repository/HG38/gnomad/4.0"."/rocksdb/";
my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => "NGS2019_2653");#NGS2021_3659 NGS2023_6384
my $rg = GenBoNoSqlRocksAnnotation->new(dir=>$dir_out,mode=>"r",index=>"genomic",name=>$chr);
my $vs = $project->getChromosome($chr)->getVariations();
warn "end get";

my $t = time;
my @vals;
foreach my $v (sort{$a->start <=> $b->start} @$vs){
	warn $v->rocksdb_id;
	warn $v->name;
	warn Dumper  $v->getChromosome->rocksdb("gnomad")->value($v->rocksdb_id);
	next unless $v->getGnomadHO();
	my $k1 = sprintf("%010d", $v->start);
	my $allele = $v->alternate_allele;
	push(@vals,$k1."!".$allele);
	my $h = $rg->gnomad($v->gnomad_id);
	warn Dumper $h;
	#warn $v->gnomad_id;
	#my $h;
	#die();
#	my $r = $v->buffer->get_gnomad($v->getChromosome->name,$v->type_public_db,$v->start,$v->alternate_allele);
#	print_erreur($v,$h,"AN") if  $v->cadd_score() ne $h->{an};
##	warn Dumper $h;
#	#my (@s) = unpack("w4 f2 A3 A3 ",$h);
	#print_erreur($v,$h,"AN") if  $v->getGnomadAN() ne $h->{an};
	#print_erreur($v,$h,"AC") if  $v->getGnomadAC() ne $h->{ac};
	#print_erreur($v,$h,"XY") if  $v->getGnomadAC_Male() ne $h->{xy};
#	
	#print_erreur($v,$h,"MAX ".abs($v->max_pop_freq() - $h->{max})." ".$h->{max_pop}) if  abs($v->max_pop_freq() - $h->{max})> 0.01;
	#print_erreur($v,$h,"MIN ". abs($v->min_pop_freq() - $h->{min})) if  abs($v->min_pop_freq() - $h->{min})> 0.01;
	#delete $v->{buffer};
	#delete $v->{project};
	#die() if $v->getGnomadAN();
}
#my $enc = Sereal::Encoder->new({compress=>Sereal::SRL_ZSTD,compress_threshold=>0});
#$enc->sereal_encoder->encode_to_file("NGS2023_6384.sereal",\@vals);
warn "OK => ". abs(time -$t);

$t =time;
undef $project->{objects};
warn "delete :".abs(time -$t);
undef $rg;

warn abs(time -$t);
warn time;
exit(0);

sub print_erreur {
	my ($v,$h,$text) =@_;
	warn $v->getGnomadAN()." ".$h->{an};
	warn $v->getGnomadAC()." ".$h->{ac};
	warn $v->getGnomadAC_Male()." ".$h->{xy};
	warn $v->getGnomadHO()." ".$h->{ho};
	warn $v->max_pop_freq()." ".$v->max_pop_name." => ".$h->{max}." ".$h->{max_pop};
	warn $v->min_pop_freq()." "." ".$v->min_pop_name." => ".$h->{min}." ".$h->{min_pop};
	warn $v->gnomad_id;
	warn Dumper $h;
	warn $text;
	warn $v->getPatients->[0]->name;
	warn "gnomad :".$v->gnomad_id;
	warn "vcf :".$v->vcf_id;
	warn $v->id;
	warn "=========================================================================================\n";
}
