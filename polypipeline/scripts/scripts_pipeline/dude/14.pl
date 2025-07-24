#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use GenBoNoSqlLmdb;

use Carp;
use strict;
use Set::IntSpan::Fast::XS ;
use Data::Dumper;
use GBuffer;
use Getopt::Long;
use Carp;
 use JSON::XS;
 use List::MoreUtils qw(natatime uniq);
  use List::Util qw(sum);
 use Statistics::Descriptive;
 
 my $buffer = new GBuffer;
my $project_name= "";
my $max_length_event = 3;
my $fork;
GetOptions(
	'project=s' => \$project_name,
	'fork=s' => \$fork,
);
die("man !!!  no fork ") unless $fork;
my $project = $buffer->newProjectCache( -name 			=> $project_name );
if ($project->isGenome){
	$max_length_event = 9;
} 
my $jobs;

my $pm2 = new Parallel::ForkManager($fork);
my $all_tab;
my $dir_out= $project->getVariationsDir("dude");
$project->getCaptures();
foreach my $patient (@{$project->getPatients}){
	#my $files = $patient->getVariationsFiles();
	my $file = $dir_out."/".$patient->name.".dude.lid";
	my $filegz = $dir_out."/".$patient->name.".dude.lid.gz";

	if (-e $filegz){
		unlink $filegz;
		unlink $filegz.".tbi";
	}
	if (-e $file){
		unlink $file;
		unlink $file;
	}
	
	
}
$pm2->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data ) =
			  @_;
			 if ($data) {
			 	my $id = $data->{job};
#			 	warn $id;
			 	delete $jobs->{$id};
			 	my $patient_name = $data->{patient};
			 	
			 	my $file = $dir_out."/".$patient_name.".dude.lid";
			 	open (BED,">>$file");
				print BED join("\n",@{$data->{lines}})."\n" if $data->{lines};		
				close(BED);		  
			 #	$all_tab->{$data->{patient}}->{$data->{chromosome}} = ;
			 #	warn Dumper $all_tab->{$data->{patient}}->{$data->{chromosome}};
#			 	warn scalar keys %$jobs;
			 }
		}
	);
	
my $id = time;	
my $no =  $project->noSqlCnvs("r");
my $control = $no->get("raw_data","controls");
my $cai_count = $no->get("raw_data","cai_count");
$project->getPatients();
	$project->preload();
$project->buffer->dbh_deconnect();
my $bcftools = $buffer->software("bcftools");
foreach my $patient (@{$project->getPatients}){
	#next unless $patient->name eq "ALL_RAF";
	$id ++;
	$jobs->{$id} ++; 
	 
	 my @pcontrol = map{$_->{id}} @{$control->{$patient->getRun->id}};
	 die() unless @pcontrol;
	

#	my $variations = $patient->getStructuralVariations();
	 my $tree = Set::IntervalTree->new;
	 my $htree;
	 
	 foreach my $chr (@{$project->getChromosomes}){
	 
	 	my $files = $patient->getVariationsFiles();
	 	my $chr_name = $chr->ucsc_name();
	 	$htree->{$chr->name} = Set::IntervalTree->new;;
	 	foreach my $file (@$files){
	 		next unless $file =~/vcf/;
	 		open("VCF","$bcftools  view -g het $file -r $chr_name -H | cut -f 2 |") or die;
	 		
	 		while(<VCF>){
	 			chomp($_);
	 			$htree->{$chr->name}->insert($_,$_-1,$_+1);
	 		}
	 		
	 		close VCF;
	 	}
	 	#$htree->{$chr->name} = $chr->getTreeVariants($patient);
	 }
	 

#	warn $file;
	my $nb_genes;
	my @beds;
	
	foreach my $chr (@{$project->getChromosomes}){
	#	next if $chr->name ne "1";
		 my $pid      = $pm2->start() and next;	
	   	#$project->buffer->dbh_reconnect();
		#next if $chr->name ne "21";
		 	my %hgenes;
			delete $chr->{natatime};
			my $vtree = $htree->{$chr->name};
			my @primer_he;
			#my $primers = [sort{$a->start <=> $b->start} @{$chr->getPrimers()}];
			my $first ;
			my $i = 0;
			my $primers;
			#my $previous;
			my $index_primer = 0;
			my $intspan_dup_ho = Set::IntSpan::Fast::XS->new;
			my $intspan_dup_he = Set::IntSpan::Fast::XS->new;
			my $intspan_del_ho = Set::IntSpan::Fast::XS->new;
			my $intspan_del_ho_true = Set::IntSpan::Fast::XS->new;
			my $intspan_del_he = Set::IntSpan::Fast::XS->new;
			my $intspan_deletion = Set::IntSpan::Fast::XS->new;
			my $intspan_duplication = Set::IntSpan::Fast::XS->new;
			my $res;
			my $end_primers;
			while (my $primer = $chr->next_primer){
#				warn $primer->start." ".$primer->end()." ".$primer->id;
				$index_primer ++;
				my $debug;
				my $genes = $primer->getGenes();
				my $name ;
				$name = $genes->[0]->external_name if @$genes;
				
				
				push(@$primers,$primer->id);
				my $level = $primer->level($patient);
				
				my $score = $primer->cnv_score($patient);
				warn $level." ".$score;
				next if $level ne  2 and $level ne 1;
				next if $score == -1;
				
				my $nb1 = scalar(@pcontrol)	;
			
			
				my $nb2 =  scalar( grep {$_ >0.6  and $_ < 1.4 } @{$primer->{zscore_data}->{$patient->getRun->id}});
				next if ($nb2 < $nb1*0.7 && $nb1 > 10);
				warn " $index_primer ".$score ." nb2:$nb2 $nb1 ".$level." ".$patient->name if $debug;
				my $status ;
				my $type;
				if ($level == 1 ){
					## deletion 	
						
					my $results = $vtree->fetch($primer->start-50,$primer->end+50);
					next if (@$results);
					
					if ($score <= 0.1  ) {
						$status = "1/1";
						$type = "del";
						$intspan_del_ho->add($index_primer);
						
						#$primer->{$patient->id}->{type} = "del";
						#$primer->{$patient->id}->{status} = "1/1";
						#push(@$end_primers,$primer);
					}
					else {
						$status = "0/1";
						$type = "del";
						$intspan_del_he->add($index_primer);
						warn "\t hetero  $index_primer" if $debug;
					}
				
				}
				elsif ($level == 2 ) {
					
					if ($score >= 1.9  ) {
						$status = "1/1";
						$type = "dup";
						#$primer->{$patient->id}->{type} = "dup";
						#$primer->{$patient->id}->{status} = "1/1";
						#push(@$end_primers,$primer);
						
						$intspan_dup_ho->add($index_primer);
					}
					else {
							$status = "0/1";
							$type = "del";
							$intspan_dup_he->add($index_primer);
						#	warn "coucou";
					}
				}
				else {
					warn $level;
					confess();
				}
				
				$res->{$index_primer}->{id} = $primer->id;
				$res->{$index_primer}->{start} = $primer->start;
				$res->{$index_primer}->{end} = $primer->end;
				$res->{$index_primer}->{type} = $type;
				$res->{$index_primer}->{status} = $type;
				next unless $type;
				}
				
			
				#my ($intspan,$res,$chr,$limit) =@_;
				
				my $genes_del_ho = get_genes($intspan_del_ho,$res,$chr,2);
				
				push(@$end_primers,@{select_primers($patient,$chr,$genes_del_ho,"del",1)});
				my $genes_del_he = get_genes($intspan_del_he,$res,$chr,2);
				push(@$end_primers,@{select_primers($patient,$chr,$genes_del_he,"del",undef,1)});
				my $genes_dup_ho = get_genes($intspan_dup_ho,$res,$chr,2);
																		
				push(@$end_primers,@{select_primers($patient,$chr,$genes_dup_ho,"dup",1)});
#				warn $intspan_dup_he->as_string;
				my $genes_dup_he = get_genes($intspan_dup_he,$res,$chr,2);
					warn Dumper $genes_dup_he;
				#warn scalar @{select_primers($patient,$chr,$genes_dup_he,"dup",1)};
					
				push(@$end_primers,@{select_primers($patient,$chr,$genes_dup_he,"dup",undef,1)});
				my %dj;
				my $lines;
				foreach my $primer (sort{$a->start <=> $b->start} @$end_primers){
					next if exists $dj{$primer->id};
					 $dj{$primer->id} ++; 
					my @genes;
					foreach my $gene (@{$chr->getGenesByPosition($primer->start,$primer->end) }){
		 		 		push(@genes,$gene->name()."-".$gene->external_name());
		 		 	}
		 		 	my $type = $primer->{$patient->id}->{type};
		 		 	die() unless exists $primer->{$patient->id}->{type};
		 		 	my $status = $primer->{$patient->id}->{status};
					my $line =  $chr->name."\t".$primer->start."\t".$primer->end."\t".$status."\t".$type."\t".join(";",sort{$a cmp $b} @genes)."\t".$primer->cnv_score($patient)."\t".join(";",map{sprintf("%.1f",$_)} @{$primer->{zscore_data}->{$patient->getRun->id}})."\t".int($cai_count->{$primer->id}->{$patient->id})."\t".$primer->sd($patient);
					push(@$lines,$line);
					#print BED $line."\n";
				}
				my $hres;
				 $hres->{patient} = $patient->name;
				 $hres->{chromosome} = $chr->name;
				$hres->{lines} = $lines;
				 $hres->{job} = $id;			  
				$pm2->finish(0,$hres);
				
	}#end chromoseome

		
}	
	
	$pm2->wait_all_children;
	warn "******************************************************************** ";
	warn "*************************** end fork ******************************* ";
	warn "******************************************************************** ";
	warn Dumper $all_tab;
	#
			#warn "end patient";

	if  (scalar keys %$jobs ne 0){
		die("problem ");
	}
		my $bgzip = $buffer->software("bgzip");
	my $tabix = $buffer->software("tabix");
	
	foreach my $patient (@{$project->getPatients}){
		my $file = $dir_out."/".$patient->name.".dude.lid";
		my $filegz = $dir_out."/".$patient->name.".dude.lid.gz";
		system("$bgzip $file && $tabix -p bed $file.gz ");
	}
	
sub get_genes {
	my ($intspan,$res,$chr,$limit) =@_;
my $iter = $intspan->iterate_runs();
my %hgenes;
			while (my ( $from, $to ) = $iter->()) {
			 	for my $index_primer ($from .. $to) {
			 		my $start = $res->{$index_primer}->{start};
			 		my $end = $res->{$index_primer}->{end};
			 		foreach my $gene (@{$chr->getGenesByPosition($start,$end) }){
			 			$hgenes{$gene->id} ++;
		 			}
				 }
			}
		my @genes = grep {$hgenes{$_}>= $limit} keys %hgenes;	
		return \@genes;
		
}	

sub select_primers {
	my ($patient,$chr,$genes,$type,$ho,$debug) = @_;
	my $return_primers =[];
	my $limit_del = 0.6;
	my $limit_dup = 1.4;
	my $status = "0/1";
	if ($ho == 1){
		 $limit_del = 0.1;
		 $limit_dup = 1.9;
		 $status = "1/1";
	}
	my $run_id = $patient->getRun->id;
	foreach my $g (@$genes){
					my $hgene =   $project->newGene($g);
					my $gene_intspan = Set::IntSpan::Fast::XS->new; 
					foreach my $t (@{$hgene->getTranscripts}){
						#warn $t->getGenomicSpan->as_string();
						$gene_intspan = $gene_intspan->union( $t->getSpanCoding)
					}
					
					my $nb_bad = 0;
					my @tprimers ;
					foreach my $pt (@{$project->getPrimersByPosition($chr,$hgene->{start},$hgene->{end})}){
						#warn $pt if $debug;
					 	next unless exists $pt->runs_object->{$run_id};
						if ($pt->level($patient) == -1){
							$nb_bad ++ ;
						}
						else {
							push(@tprimers,$pt);
						}
					}
					
					next if $nb_bad > 1;
					next  if scalar @tprimers ==0; 
					my @v ;
					my $coding = Set::IntSpan::Fast::XS->new; 
					my $i =0;
					foreach my $p (@tprimers) {
						push(@v,$p->cnv_score($patient));
						my $z = $gene_intspan->intersection($p->getGenomicSpan());
						next if $z->is_empty();
							#$coding->add($i);
						#if (@{$p->getExons()} > 0){
							$coding->add($i);
						#}
						$i++;
					}
					my $s = smooth_data(\@v);
					my $in = Set::IntSpan::Fast::XS->new; 
					my $ct =0;
#					warn $type." ".$limit_dup;
					for (my $i=0;$i<@$s;$i++){
						 if ($type eq "del") {
						 	
						 	$in->add($i) if $s->[$i] <= $limit_del;
						 	$ct ++ if $s->[$i] >= $limit_dup;
						 }
						 elsif($type eq "dup") {
						 	$in->add($i) if $s->[$i] >= $limit_dup; 
						 	$ct ++if $s->[$i] <= $limit_del;
						 }
					}
					next if $ct > 2;
				
					my $iter = $in->iterate_runs();
					my $find;
    				while (my ( $from, $to ) = $iter->()) {
    					next if abs($from-$to) < 1;
    					$find =1;
    				}
    				next  unless $find;
					my @ind = $in->as_array;
					my $z = $coding->intersection($in);
					#next if $z->is_empty;
			#		warn $z->as_string."\n ==>".$coding->as_string()."++".$in->as_string();
					my $iter = $in->iterate_runs();
					my $max_len = 0 ;
					while (my ( $from, $to ) = $iter->()) {
						$max_len =  abs($from-$to) if abs($from-$to) > $max_len;
					}
				#	warn Dumper @ind;
					my $nb = scalar($in->as_array);
#					warn "----- $nb ".scalar(@tprimers) if $debug;
					unless ($ho){
						if ( $max_len < $max_length_event && (($nb/scalar(@tprimers))< 0.25)){
						#if ( $max_len < $max_length_event && (($nb/scalar(@tprimers))< 0.25)){
							#warn "coucou";
							next;
						}
					}
					for my $index (@ind){
						$tprimers[$index]->{$patient->id}->{type} = $type;
						$tprimers[$index]->{$patient->id}->{status} = $status;
						push(@$return_primers,$tprimers[$index]);
					}
					
				}
			return $return_primers;
}
sub smooth_data {
	my ($data) = shift;
	my @data2;
		for (my $i = 0;$i<@$data;$i++){
		my @trio;
		 push (@trio,$data->[$i-2]) if ($i>1);
		 push (@trio,$data->[$i+1]) if ($i+1<@$data);
		 push (@trio,$data->[$i-1]) if ($i>0);
		 push (@trio,$data->[$i+1]) if ($i+1<@$data);
		  push (@trio,$data->[$i+2]) if ($i+2<@$data);
		  push (@trio,$data->[$i]*2) ;#if ($i+1<@data);
		  my $z = sum(@trio);
		  push(@data2,$z/(scalar(@trio)+1));
		 
	}
	return \@data2;
}	
exit(0);