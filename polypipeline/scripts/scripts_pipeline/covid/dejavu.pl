#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use GBuffer;
#use lib "/bip-d/perl/ensembl64/ensembl/modules";
use Bio::SearchIO;
use strict;
use Getopt::Long;
use Carp;
use JSON::XS;
use Getopt::Long; 
use List::MoreUtils qw(uniq);
#use ensembl_buffer;
use Storable qw/freeze thaw nfreeze nstore_fd nstore/;
use Set::IntSpan::Fast::XS ;
use Sys::Hostname;
use Storable qw/retrieve store/;
use  Set::IntSpan::Island;
use Parallel::ForkManager;
use Array::IntSpan;
use Data::Dumper;
use Array::IntSpan;
use GenBoNoSqlLmdb;
use GenBoNoSqlAnnotation;
use GenBoNoSqlDejaVu;
#use Bio::Ensembl::Registry;
use Storable qw(dclone);
use List::MoreUtils qw(uniq natatime);
use Set::IntervalTree;
use GenBoGene;
use GenBoTranscript;
use GenBoProtein;
use GenBoExon;
use GenBoIntron;
 use Storable;
 
 my $covid = {
 	NGS2020_2994 => 1,
 	NGS2020_2978 => 1,
 	NGS2020_2976 => 1,
 	NGS2020_3002 => 1,
 	NGS2020_3016 => 1,
 	NGS2020_3023 => 1,
 	NGS2020_3030 => 1,
 	NGS2020_3045 => 1,
 	NGS2020_3049 => 1,
 	NGS2020_3067 => 1,
 };
 my $chr_name;
GetOptions(

	'chr=s' =>\$chr_name,
);
 my @chromosomes = (1..22,'X','Y','MT');
 #foreach my $chr (@chromosomes){
 	 my $total ={};
 	foreach my $p (keys %$covid){
 		one_project($p,$total,$chr_name);
 	
 	}
 	open (DV,">$chr_name.dv.txt");
 	foreach my $k (sort{$a <=> $b} keys %$total){
	foreach my $id ( sort {$total->{$k}->{$a}->{pos} <=> $total->{$k}->{$b}->{pos} } keys %{$total->{$k}}){
		print DV $total->{$k}->{$id}->{string};
	}
	}
	close DV;


sub one_project {
my ($pname,$total,$chr_name) = @_;

my $buffer = new GBuffer;
my $project_name= "$pname";
warn "start ".$pname;


my $project = $buffer->newProjectCache( -name 			=> $project_name );
my $capture =  $project->getCapture();
 my $query = $project->buffer->getQuery->listAllProjectsNameByCaptureId($capture->id());
 my $nb_twist = 0;
 my$htwist;
 foreach my $pn ( @$query){
 	next if exists $covid->{$pn};
 	my $b = new GBuffer;
 	my $p = $b->newProject( -name 			=> $pn );
 	
 	foreach my $patient (@{$p->getPatients}){
 		my $c = $patient->getCapture;
 		if ($c->id eq $capture->id()){
 			$nb_twist ++;
 			$htwist->{$pn}->{$patient->name} = undef;
 			#push(@{$htwist->{$pn}},$patient->name);
 		}
 		
 	}
 }
 
 my $fork = 40;
 my $pm = new Parallel::ForkManager($fork);

	my $id = 0;
	my $hrun;
	$pm->run_on_finish(
    sub { 
    	my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$h)=@_;
  
    	unless (defined($h) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				die();
				return;
			}
			foreach my $k (keys %{$h->{total}}){
				foreach my $id ( keys %{$h->{total}->{$k}}){
		 		$total->{$k}->{$id} = $h->{total}->{$k}->{$id};
				}
			}
			
    }
	);
 	my $chr_query = $project->getChromosome($chr_name);
	my $list = $project->getVariantsList($chr_query);
	my $nb = int(scalar(@$list)/$fork+1);
	my $iter = natatime($nb, @$list);
 #warn Dumper $query;
 my %h  = @$query;
 my $old_chr;
 warn "start";
 while( my @tmp = $iter->() ){
 	my $pid = $pm->start and next;
		my $t =time;
		my $res;
		 $res->{total} = compute_list($project,\@tmp,$htwist,$nb_twist);
#		warn "end ".abs($t-time);
		
		$pm->finish(0,$res);
	}
	$pm->wait_all_children();
 }
 
sub compute_list{
	my($project,$list,$htwist,$nb_twist) = @_;
	$project->setListVariants($list);
	my $total;
#	#foreach my $id (@$list){
	while (my $v = $project->nextVariant){
		warn scalar(@{$project->{list_variant}}) if scalar(@{$project->{list_variant}})%1000 == 0;
		my $dv =  $project->getDejaVuInfos($v->id);
		my $twist;
		my $other;
		$twist->{he} = 0;
		$twist->{ho} = 0;
		foreach my $p (keys %$dv){
			next if exists $covid->{$p};
		
			if (exists $htwist->{$p}){
				my @ps = split(";",$dv->{$p}->{string});
				foreach my $tt (@ps){
					my ($pp,$type) = split(":",$tt);
					if (exists $htwist->{$p}->{$pp}){
						
						$twist->{ho} ++ if $type ==1;
						$twist->{he} ++ if $type ==2;
					}
				}

			}
			else {
				$other->{he} += $dv->{$p}->{he};
				$other->{ho} += $dv->{$p}->{ho};
			}
			
		}
		my ($chr,$pos,$a1,$a2) = split("_",$v->vcf_id);

		my $o = $project->getChromosome($chr);
		$total->{$o->karyotypeId}->{$v->vcf_id}->{pos} = $pos;
		$total->{$o->karyotypeId}->{$v->vcf_id}->{string} = $chr."\t".$pos."\t".$a1."\t".$a2."\t"."\t".$twist->{he}."\t".$twist->{ho}."\t".$nb_twist."\n";
		#print $chr."\t".$pos."\t".$a1."\t".$a2."\t"."\t".$twist->{he}."\t".$twist->{ho}."\t".$nb_twist."\n";
	}
	return $total;
}
				
