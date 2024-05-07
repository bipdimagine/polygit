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
use Bio::DB::HTS;
use GBuffer;
use Bit::Vector::Overload;
use Sys::Hostname;


my $fork = 1;
my ($project_name, $chr_name,$annot_version);
my $ok_file;
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'chr=s'        => \$chr_name,
	'annot_version=s'    => \$annot_version,
	'file=s' => \$ok_file,
);

 if ($ok_file && -e $ok_file) {
 	system("rm $ok_file");
 }
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
	
		my $sam;
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
				warn $patient->getBamFile;
				$sam->{$patient->name} =  Bio::DB::HTS->new(-bam=>$patient->getBamFile, -fasta=>$project->getGenomeFasta());
			}
			my $res;
			$res->{run_id} = $run_id;
			my @strict_denovo;
			my $hVarDeleted;
			
			for (my $vector_id=$r->[0];$vector_id<= $r->[1];$vector_id ++){
				next unless $vector_denovo->contains($vector_id);
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
				if ($r > 40 ) {
					$limit = 5;
					$percent =0.05;
				}
				
				
			
				$debug =1 if $var->name eq '17-78064045-del-187';
				my $to_keep;
				warn "coucou " if $var->gnomad_id eq '17-78064045-del-187';
				foreach my $parent (@{$family->getParents()}) {
					my $to_keep;
						if ($var->isCnv && $var->isSrPr){
								$to_keep =check_cnv($var,$children,$parent);	
								warn "res ".$to_keep if $debug;
								warn $vector_id if $debug;
						}
					elsif ($var->isSrPr){
								$to_keep = check_srpr($var,$children,$parent);	
						}
					elsif ($var->isVariation) {
						next if ($parent->depth($var->getChromosome->name,$var->start,$var->start) < $limit);
						$to_keep = check_substitution($sam->{$parent->name},$var->getChromosome, $var->start,$var->sequence,$percent);
						#$hVarDeleted->{$vector_id} ++;
					}
					elsif ($var->isDeletion) {
						my $chr = $var->getChromosome->fasta_name();
						my $start = $var->start() - 2;
						my $end = $var->end() + 2;
						my $locus = $chr.':'.$start.'-'.$end;
						$to_keep = check_del($sam->{$parent->name},$var->getChromosome, $start,$var->delete_sequence,$limit);
					}
					elsif ($var->isInsertion) {
						my $chr = $var->getChromosome->fasta_name();
						my $length = 3;
						if ($length > 3) { $length = $var->length(); }
						my $start = $var->start() - $length - 1;
						my $end = $var->start() + $length;
						my $locus = $chr.':'.$start.'-'.$end;
						$to_keep = check_ins($sam->{$parent->name},$var->getChromosome, $start,$var->sequence,$limit);
					}
					else {
						die();
					}
					if ($to_keep) { $hVarDeleted->{$vector_id} ++; }
					last unless $to_keep;
				}
				
			}
			my $vector = $chr->getNewVector();
			warn $hVarDeleted->{168908};
			
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
system("date > $ok_file") if $ok_file;
exit(0);



sub check_srpr {
	my ($var,$children,$parent) = @_;
	my @v1 =  $parent->sr_raw($var->getChromosome,$var->start);
	return if ($v1[0] < 5 or $v1[1]>5 or  $v1[2]>5);
	@v1 =  $parent->sr_raw($var->getChromosome,$var->end);
	return if ($v1[0] < 5 or $v1[1]>5 or  $v1[2]>5);
	return 1;
}


sub check_cnv {
	my ($var,$children,$parent) = @_;
	my $debug;
	$debug =1 if $var->gnomad_id eq '17-78064045-del-187';
	$debug =1 if $var->name eq '17-78064045-del-187';
	my $dp = $var->getNormDP($parent);
	my $dpc = $var->getNormDP($children);
	warn " dp parent ".$dp." ".$dpc if $debug == 1; 
	return if $var->getMeanDP($parent) < 5;
	return if $dp < 5; 
	
	$dpc = 0.00000001 if $dpc == 0; 
	
	my $pc = int((($dpc-$dp)/$dpc)*100);
	$pc = abs($pc);
	return if $pc < 30; 
	
	warn $pc if $debug;
	if ($var->isDeletion) {
		return if $pc > -30; 
		return 1;
		#return 1 if 
	}
	elsif ($var->isLargeDuplication) {
		return 1 if $pc > 30; 
		return;
		#return 1 if 
	}
	else {
		confess($var->name." ".$var->type)
	}
	
	
	
}

sub check_substitution {
	my ($sam,$chr,$pos, $sequence_alt,$limit) = @_;
	
	my ($start,$end) = get_start_end($chr,$pos,$chr->sequence($pos,$pos));
	
	my $res = pileup($sam,$chr,$start,$end);
	my $count = 0;
	foreach my $pos (sort {$a <=> $b} keys %$res) {
		if (exists $res->{$pos}->{$sequence_alt}){
			my ($d) = values %{$res->{$pos}->{ref}};
			
			my $v =  ( $res->{$pos}->{$sequence_alt} / ($d + $res->{$pos}->{$sequence_alt}));
			$count ++ if $v > $limit;
		}
	}
	return if $count;
	return 1;
	
}

sub check_del {
	my ($sam,$chr,$pos, $sequence_alt,$limit) = @_;
	my ($start,$end) = get_start_end($chr,$pos+1,$sequence_alt,1);
	my ($start1,$end1) = get_start_end($chr,$pos-1,$sequence_alt,1);
	$start = $start1 if ($start1 < $start);
	$end = $end1 if ($end1 > $end);
	my $res = pileup($sam,$chr,$start,$end);

	my $count = 0;
	foreach my $pos (sort {$a <=> $b} keys %$res) {
		$count += $res->{$pos}->{del} if exists $res->{$pos}->{del};
	}
	return if ($count >= $limit);
	return 1;
	
}


sub pileup {
	my ($sam,$chr,$start,$end, $locus) = @_;
	my %res;
	my $callback = sub {
		my ($seqid, $pos1, $pileups) = @_;
		return if ($pos1 < $start);
		return if ($pos1 > $end);
	
		 my $nb_reads = scalar(@$pileups);
		
		foreach my $pileup (@$pileups){
			
			
			if ($pileup->indel > 0){
				$res{$pos1}->{ins} ++;
			}
			elsif ($pileup->indel < 0){
					$res{$pos1}->{del} ++;
			}
			else {
				my $b     = $pileup->alignment;
				my $ref = $chr->sequence($pos1,$pos1);
				my $qbase  = substr($b->qseq,$pileup->qpos,1);
				if ($ref eq $qbase){
					$res{$pos1}->{ref}->{$ref} ++;
				}else {
					$res{$pos1}->{$qbase} ++;
				}
			}
		}
	};
	$sam->fast_pileup($chr->fasta_name.":$start-$end", $callback);
	return \%res;
}


sub get_start_end {
	my ($chromosome,$pos,$ref1,$debug) = @_;
	my $start =$pos;
	my $ref2 = $chromosome->sequence($start,$start);
	do {
		$start --;
		$ref2 = $chromosome->sequence($start,$start);
	} while ($ref2 eq $ref1);
	$start ++ if $start < $pos;
	my $end = $pos;
	
	$ref2 = $chromosome->sequence($end,$end);
	
	 do {
		$end ++;
		$ref2 = $chromosome->sequence($end,$end);
	} while ($ref2 eq $ref1);
	$end -- if $end > $pos;
	return($start,$end);
	
	
}

sub check_ins {
	my ($sam,$chr,$pos, $sequence_alt,$limit) = @_;
	my ($start,$end) = get_start_end($chr,$pos,$chr->sequence($pos,$pos));
	my $res = pileup($sam,$chr,$start,$end);
	my $count = 0;
	foreach my $pos (sort {$a <=> $b} keys %$res) {
		$count += $res->{$pos}->{ins} if exists $res->{$pos}->{ins};
	}
	return if ($count >= $limit);
	return 1;
	#return 1;
#	my $cmd = $buffer->getSoftware('samtools')." mpileup $bam_file -r $locus 2>/dev/null";
#	warn $cmd;
#	die($locus);
#	my $res = `$cmd`;
#	my @lRes = split("\n", $res);
#	warn Dumper @lRes;
#	my $count = 0;
#	foreach my $line (@lRes) {
#		my @lCol = split(' ', $line);
#		
#		my $this_count = ($lCol[4] =~ tr/\+//);
#		$count += $this_count;
#		return if ($count >= $limit);
#	}
#	return 1;
}