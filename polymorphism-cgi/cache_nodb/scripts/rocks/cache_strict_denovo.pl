#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use Parallel::ForkManager;
use Storable qw(store retrieve freeze thaw);
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use Cwd 'abs_path';
use Digest::MD5 qw(md5_hex);
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use Bit::Vector::Overload;
use Sys::Hostname;

my $fork = 1;
my ($project_name, $chr_name,$annot_version);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'chr=s'        => \$chr_name,
	'annot_version=s'    => \$annot_version,
);
warn "******";
warn hostname;
warn "******";
unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }



warn "\n### CACHE: strict-denovo model step\n";
my $nbErrors = 0;
my $buffer = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProjectCache( -name => $project_name, -cache => '1', -typeFilters => 'familial' );
my $limit = 2;
if ($project->isGenome){
	$limit = 10;
}
if ($annot_version) {
	$project->changeAnnotationVersion($annot_version);
}
$project->getPatients();
my $chr = $project->getChromosome($chr_name);

my $hResults;
my $hErrors;

my $nbErrors = 0;

my $pm = new Parallel::ForkManager($fork);

my $vector_denovo;
my $total_job;
foreach my $patient (@{$project->getPatients()}) {
	 	$vector_denovo->{$patient->name} =  $chr->getNewVector();
}
$pm->run_on_finish (
	sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $res) = @_;
		my $pname = $res->{patient};
		unless (exists $vector_denovo->{$pname} ) {
			die();
		}
		$vector_denovo->{$pname} +=   $res->{vector};
		delete $total_job->{$res->{run_id}};
		
	}
);
my $id =time;
$project->preload_patients();


my $no      = $chr->get_rocks_variations("r");
 my $ranges = $no->ranges($fork);
$no->close;

	$project->disconnect();
	$buffer->{dbh} = "-";
	
foreach my $family (@{$project->getFamilies()}) {
	foreach my $children  (@{$family->getChildren}){
		#my $vector_denovo = $family->getVector_individual_denovo($chr,$children)->Clone();
		#my $vector_ratio_name = $children->name . "_ratio_" . 20;

		#my @lVarIds = @{$chr->getListVarVectorIds($vector_denovo)};
		#my $nb        = int( scalar(@lVarIds) / $fork + 1 );
		#my $iter      = natatime( $nb, @lVarIds );
		my $sam;
		#warn Dumper  $chr->{rocks};
		#delete $chr->{rocks}->{"vector-r"};
		#delete $chr->{rocks};
		my @tmp;
		foreach my $r (@$ranges){
		#while ( my @tmp = $iter->() ) {
			my $run_id = $id ++;
			$total_job->{$run_id} ++;
			
			$pm->start() and next;
			my $vector_denovo = $family->getVector_individual_denovo($chr,$children)->Clone();
			my $no = $project->getChromosome($chr_name)->get_rocks_variations("r");
			delete $no->{rocks};
			foreach my $patient (@{$family->getParents()}) {
				
				$sam->{$patient->name} =  Bio::DB::Sam->new(-bam=>$patient->getBamFile, -fasta=>$project->getGenomeFasta());
			}
			my $res;
			$res->{run_id} = $run_id;
			my @strict_denovo;
			my $hVarDeleted;
			
			for (my $vector_id=$r->[0];$vector_id<= $r->[1];$vector_id ++){
					next unless $vector_denovo->contains($vector_id);
					warn $vector_id;
			#foreach my $vector_id (@tmp){ 
				my $debug;
				push(@tmp,$vector_id);
				my $var = $no->get_index($vector_id);
		
				unless ($var) {
				confess();
				}
			confess() unless $var->id;
				$var->{buffer}  = $buffer;
				$var->{project} = $project;
				my $r = $var->getRatio($children);
				my $percent;
				if ($r > 40 ){
					$limit = 5;
					$percent =0.08;
					}
				
				
			
				$debug =1 if $var->id eq "7_70231236_C_G";
				
				#die() if $debug;
				my $to_keep;
				foreach my $parent (@{$family->getParents()}) {
					my $to_keep;
					if ($var->isLarge()) {
						
						#push(@strict_denovo,$vector_id);
						next;
					}
					elsif ($var->isFoundBySVCaller($children)) {
						#warn "coucou SV****" if $debug;
						#warn "coucou SV ++++ ".$vector_id unless $debug;
						if ($var->length > 60 ) {
							$hVarDeleted->{$vector_id} ++;
						}
						elsif ($var->length == 1 &&  length($var->var_allele()) > 60) {
							$hVarDeleted->{$vector_id} ++;
						}
						#$hVarDeleted->{$vector_id} ++;
						next;
					}
					elsif ($var->isVariation) {
						$to_keep = check_substitution($sam->{$parent->name} , $var, $limit,$percent,$parent);
						warn "keep var ".$to_keep." ".$parent->name if $debug;
						#$hVarDeleted->{$vector_id} ++;
					}
					elsif ($var->isDeletion) {
						my $chr = $var->getChromosome->fasta_name();
						my $start = $var->start() - 2;
						my $end = $var->end() + 2;
						my $locus = $chr.':'.$start.'-'.$end;
						$to_keep = check_del($parent->getBamFile, $var->type_public_db(), $locus, $limit);
						warn $to_keep if $debug;
					}
					elsif ($var->isInsertion) {
						my $chr = $var->getChromosome->fasta_name();
						my $length = 3;
						if ($length > 3) { $length = $var->length(); }
						my $start = $var->start() - $length - 1;
						my $end = $var->start() + $length;
						my $locus = $chr.':'.$start.'-'.$end;
						$limit = 1;
						$to_keep = check_ins($parent->getBamFile, $var->type_public_db(), $locus, $limit);
						warn $to_keep." ".$parent->name if $debug;
					}
					else {
						die();
					}
					if ($to_keep) { $hVarDeleted->{$vector_id} ++; }
				#	else {
				#		last;
				#	}
				#die()  if $debug;
#				warn $hVarDeleted->{$vector_id}."+++"
				}
				
			}
			my $vector = $chr->getNewVector();
			
			foreach my $vector_id (@tmp){ 
				if ( $hVarDeleted->{$vector_id} == scalar(@{$family->getParents()}) ){
					$vector->Bit_On($vector_id);
				}
			}
		
			$res->{patient} = $children->name;
			$res->{vector}  = $vector;
			delete $chr->{rocks};
			sleep(5);
			$pm->finish(0, $res);
		} #end @tmp 
		
	}
	
	
}
#warn "wait";
$pm->wait_all_children();
confess() if keys %$total_job;
my $nosql = GenBoNoSql->new(dir => $chr->project->getCacheBitVectorDir().'/strict-denovo', mode => 'c');
my $rocks4 = $chr->rocks_vector("w");
foreach my $family (@{$project->getFamilies}){
	foreach my $child  (@{$family->getChildren}){
		 $rocks4->put_batch_vector_transmission($child,"ind_strict_denovo",$vector_denovo->{$child->name});
	}
}
$rocks4->write_batch();
$rocks4->close();
#$nosql->put_bulk($chr->id(), $vector_denovo);
#$nosql->close();
exit(0);



sub check_substitution {
	my ($sam, $var, $limit,$pourcent,$parent) = @_;
	my $debug = "19_13373567_T_G";
	$debug = undef if $var->id ne $debug; 
	my $chr_name = $chr->fasta_name();
	my $start = $var->start();
	my $end = $var->end();
	my $all_var = $var->var_allele();
	my $count = 0;
	my $nb_mut = 0;
	 my $nb_reads =0;
	my $callback = sub {
		my ($seqid, $pos1, $pileups) = @_;
		return if ($pos1 ne $start);
		if (scalar(@$pileups) < 3) {
			$nb_mut = 99;
			return;
		}
		 $nb_reads = scalar(@$pileups);
		
		foreach my $pileup (@$pileups){
			my $b     = $pileup->alignment;
			my $qbase  = substr($b->qseq,$pileup->qpos,1);
			warn $qbase." ".$all_var if $debug;
			$nb_mut ++ if ($qbase eq $all_var);
			last if $nb_mut/$nb_reads >= $pourcent;
			last if ($nb_mut >= $limit);
		}
	};
	
	my $callback2 = sub {
		my ($seqid, $pos1, $pileups) = @_;
		return if ($pos1 ne $start);
		if (scalar(@$pileups) < 3) {
			$nb_mut = 99;
			return;
		}
		my $nb_reads = scalar(@$pileups);
		
		foreach my $pileup (@$pileups){
			my $b     = $pileup->alignment;
			my $qbase  = substr($b->qseq,$pileup->qpos,1);
			warn $qbase." ".$all_var;
		}
	};
	$sam->fast_pileup("$chr_name:$start-$end", $callback);
	warn $nb_mut."****" if $debug;#." ".$all_var;
	if ($nb_reads == 0){
		return;
	}
	
	return if (($nb_mut/$nb_reads > $pourcent) or ($nb_mut >= $limit));
	return 1;
}

sub check_del {
	my ($bam_file, $type, $locus, $limit) = @_;
	my $cmd = $buffer->getSoftware('sambamba')." depth base $bam_file -L $locus 2>/dev/null";
	my $res = `$cmd`;
	my @lRes = split("\n", $res);
	my $name_to_check = 'DEL';
	my $col_to_check;
	my $i = 0;
	foreach my $name (split(' ', shift(@lRes))) {
		if ($name eq $name_to_check) {
			$col_to_check = $i;
			last;
		}
		else { $i++; }
	}
	foreach my $res_pos (@lRes) {
		my @lCol = split(' ', $res_pos);
		my $value = $lCol[$col_to_check];
		return if ($value >= $limit)
	}
	return 1;
}

sub check_ins {
	my ($bam_file, $type, $locus, $limit) = @_;
	#return 1;
	my $cmd = $buffer->getSoftware('samtools')." mpileup $bam_file -r $locus 2>/dev/null";
	my $res = `$cmd`;
	my @lRes = split("\n", $res);
	my $count = 0;
	foreach my $line (@lRes) {
		my @lCol = split(' ', $line);
		my $this_count = ($lCol[4] =~ tr/\+//);
		$count += $this_count;
		return if ($count >= $limit);
	}
	return 1;
}