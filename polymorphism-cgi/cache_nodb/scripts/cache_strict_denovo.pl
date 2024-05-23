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
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use Cache_Commons;
use Bit::Vector::Overload;
use Sys::Hostname;
use Bio::DB::HTS;

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
	$limit = 2;
}
if ($annot_version) {
	$project->changeAnnotationVersion($annot_version);
}
$project->getPatients();
my $chr = $project->getChromosome($chr_name);
mkdir ($project->getCacheBitVectorDir().'/strict-denovo') unless (-d $project->getCacheBitVectorDir().'/strict-denovo');
my $freeze = $project->getCacheBitVectorDir()."/lmdb_cache/$chr_name.dv.freeze";
unless (-e $freeze) {
	die ("\n\nERROR: $freeze doesn't exist... Die\n\n");
}
my $hAllVar = retrieve $freeze;
if (scalar keys %{$hAllVar} == 0) {
	my $cmd = "touch ".$project->getCacheBitVectorDir()."/strict-denovo/".$chr->id().".lite";
	`$cmd`;
	exit();
}
my $fasta_file = $project->getGenomeFasta();
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

#$project->getFamilies();
#foreach my $p (@{$project->getPatients()}) {
#		$p->callingSVMethods();
#		$p->callingMethods();
#		$p->getBamFile();
#}
	my @t = (31346,54941,64616,84747,86289,108838,116069,125388,220948);
	$project->disconnect();
	$buffer->{dbh} = "-";
	
foreach my $family (@{$project->getFamilies()}) {
	warn "..";
	
	foreach my $children  (@{$family->getChildren}) {
		my $vector_denovo = $family->getVector_individual_denovo($chr,$children)->Clone();
		my $vector_ratio_name = $children->name . "_ratio_" . 20;
#		 my $vquality2 = $chr->getVectorScore($vector_ratio_name);
#		 $vector_ratio_name = $children->name . "_ratio_" . 40;
#		  my $vquality3 = $chr->getVectorScore($vector_ratio_name);
#		my $vector_denovo = $family->getVector_individual_denovo($chr,$children);
		my @lVarIds = @{$chr->getListVarVectorIds($vector_denovo)};
		my $nb        = int( scalar(@lVarIds) / $fork + 1 );
		my $iter      = natatime( $nb, @lVarIds );
		my $sam;
		
		while ( my @tmp = $iter->() ) {
			my $run_id = $id ++;
			$total_job->{$run_id} ++;
			$pm->start() and next;
			#$project->disconnect();
			#$buffer->{dbh} = "-";
			foreach my $patient (@{$family->getParents()}) {
				warn $patient->getBamFile;
				$sam->{$patient->name} =  Bio::DB::HTS->new(-bam=>$patient->getBamFile, -fasta=>$project->getGenomeFasta());
			}
			my $res;
			$res->{run_id} = $run_id;
			my @strict_denovo;
			my $hVarDeleted;
			
			
			foreach my $vector_id (@tmp){ 
				
				my $debug;
				# 31346,54941,64616,84747,86289,108838,116069,125388,220948
				my ($find) = grep {$vector_id == $_} @t;
				
				
				my $var = $chr->getVarObject($vector_id);
				my $r = $var->getRatio($children);
				#warn "ccuiiui" if $var->id eq "7_6042930_A_G";
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
						$to_keep = check_substitution($sam->{$parent->name} , $var, $limit,$percent);
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
			$pm->finish(0, $res);
		} #end @tmp 
		
	}
	
	
}
#warn "wait";
$pm->wait_all_children();
warn "end";
confess() if keys %$total_job;

my $nosql = GenBoNoSql->new(dir => $chr->project->getCacheBitVectorDir().'/strict-denovo', mode => 'c');
$nosql->put_bulk($chr->id(), $vector_denovo);
$nosql->close();
exit(0);
#die();
#die();
#$buffer->dbh_reconnect();
#foreach my $family (@{$project->getFamilies()}) {
#	foreach my $children  (@{$family->getChildren}){
#		warn  $vector_denovo->{$children->name};
#	}
#}
# $pm = new Parallel::ForkManager(1);
#$pm->run_on_finish (
#	sub {
#		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hResGlobal) = @_;
#		if (defined($hResGlobal)) {
#			foreach my $fam_name (keys %{$hResGlobal}) {
#				delete $hErrors->{$fam_name};
#				foreach my $pat_name (keys %{$hResGlobal->{$fam_name}}) {
#					$hResults->{$pat_name} = $hResGlobal->{$fam_name}->{$pat_name};
#				}
#			}
#		}
#		else {
#			$nbErrors++;
#			print qq|No message received from child process $pid!\n|;
#		}
#	}
#);
#
#foreach my $family (@{$project->getFamilies()}) {
#	#next unless $family->name eq '19';
#	warn $family->name;
#	$pm->start() and next;
#	$buffer->dbh_reconnect();
#	my $hTmp;
#	$hTmp->{$family->name()} = undef;
#	foreach my $patient (@{$family->getPatients()}) {
#		$hTmp->{$family->name()}->{$patient->name()} = $chr->getNewVector();
#	}
#	if (scalar (@{$family->getChildrenIll()}) == 0 or scalar (@{$family->getHealthy()}) == 0) {
#		$pm->finish(0, $hTmp);
#		next;
#	}
#	else {
#		my $nb_init = $chr->countThisVariants( $family->getVariantsVector($chr) );
#		my $vector_denovo = $family->getModelVector_fam_denovo($chr)->Clone();
#		$hResults->{$family->name()}->{global}->{denovo} = $vector_denovo->Clone();
#		my $nb_denovo = $chr->countThisVariants( $hResults->{$family->name()}->{global}->{denovo} );
#		my @lVarIds = @{$chr->getListVarVectorIds($vector_denovo)};
#		foreach my $patient (@{$family->getHealthy()}) {
#			my $hVarDeleted;
#			my @lBamFiles = @{$patient->getBamFiles()};
#			foreach my $bam_file (@lBamFiles) {
#				#TODO: mettre HTS plus tard
#				my $sam = Bio::DB::Sam->new(-bam=>$bam_file, -fasta=>$project->getGenomeFasta());
#				foreach my $vector_id (@lVarIds) {
#					next if (exists $hVarDeleted->{$vector_id});
#					my $debug;
#					my ($find) = grep {$vector_id == $_} @t;
#					$debug =1 if $find;
#					my $var = $chr->getVarObject($vector_id);
#					
#					my $to_keep;
#			
#					if ($var->isLarge()) {
#						warn "large " if $debug;
#						$hVarDeleted->{$vector_id} = undef;
#						next;
#					}
#					elsif ($var->isFoundBySVCaller($patient)) {
#						warn "sv <-------> " if $debug;
#						$hVarDeleted->{$vector_id} = undef;
#						next;
#					}
#					elsif ($var->isVariation) {
#						$to_keep = check_substitution($sam, $var, $limit);
#					}
#					elsif ($var->isDeletion) {
#						my $chr = $var->getChromosome->fasta_name();
#						my $start = $var->start() - 2;
#						my $end = $var->end() + 2;
#						my $locus = $chr.':'.$start.'-'.$end;
#						$to_keep = check_del($bam_file, $var->type_public_db(), $locus, $limit);
#					}
#					elsif ($var->isInsertion) {
#						my $chr = $var->getChromosome->fasta_name();
#						my $length = 3;
#						if ($length > 3) { $length = $var->length(); }
#						my $start = $var->start() - $length - 1;
#						my $end = $var->start() + $length;
#						my $locus = $chr.':'.$start.'-'.$end;
#						$limit = 1;
#						$to_keep = check_ins($bam_file, $var->type_public_db(), $locus, $limit);
#					}
#					warn "pppppp ".$vector_id if $debug;
#					unless ($to_keep) { $hVarDeleted->{$vector_id} = undef; }
#				}
#			}
#			my $vector_to_del = Bit::Vector->new_Enum($vector_denovo->Size(), join(',', keys %{$hVarDeleted}));
#			$vector_denovo -= $vector_to_del;
#		}
#		my $vector_strict_denovo = $vector_denovo->Clone();
#		foreach my $patient (@{$family->getPatients()}) {
#			$hTmp->{$family->name()}->{$patient->name()} = $patient->getVariantsVector($chr)->Clone();
#			$hTmp->{$family->name()}->{$patient->name()}->Intersection($hTmp->{$family->name()}->{$patient->name()}, $vector_strict_denovo);
#		}
#		$pm->finish(0, $hTmp);
#	}
#}
#$pm->wait_all_children();
#$buffer->dbh_reconnect();
#
#if ($nbErrors > 0) {
#	warn ("\n\nERROR: $nbErrors found in fork... DIE...\n\n");
#	warn Dumper $hErrors;
#	die();
#}
#
#foreach my $family (@{$project->getFamilies()}) {
#	
#	foreach my $children  (@{$family->getChildren}){
#		#warn $children->name;
#		#warn  $vector_denovo->{$children->name};
#		#warn $hResults->{$children->name};
#		warn "denovo : " .$chr->countThisVariants($family->getVector_individual_denovo($chr,$children));
#		warn $chr->countThisVariants( $hResults->{$children->name} );
#		warn $chr->countThisVariants(  $vector_denovo->{$children->name} );
#		my $set  = Set::IntSpan::Fast::XS->new( $hResults->{$children->name}->to_Enum );
#	
#		#$set  = Set::IntSpan::Fast::XS->new( $vector_denovo->{$children->name}->to_Enum );
#		#warn $set->as_string();
#		my $z = $hResults->{$children->name} -  $vector_denovo->{$children->name};
#		
#		$vector_denovo->{$children->name} -=  $hResults->{$children->name};
#		warn "----";
#		warn $chr->countThisVariants(  $vector_denovo->{$children->name} );
#		warn $chr->countThisVariants( $z );
#		$set  = Set::IntSpan::Fast::XS->new( $z->to_Enum );
#		warn $set->as_string();
#	}
#	
#	foreach my $parents  (@{$family->getParents}){
#		warn $chr->countThisVariants( $hResults->{$parents->name} );
#	}
#}
#
#
#	warn 'store 1/1: nosql strict-denovo';
#	#warn Dumper $hResults;
#	#my $nosql = GenBoNoSql->new(dir => $chr->project->getCacheBitVectorDir().'/strict-denovo', mode => 'c');
#	#$nosql->put_bulk($chr->id(), $hResults);
#	#$nosql->close();



sub check_substitution {
	my ($sam, $var, $limit,$pourcent) = @_;
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