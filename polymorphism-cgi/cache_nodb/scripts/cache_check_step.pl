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
use Bit::Vector;
use Bit::Vector::Overload;
use Carp;
use GBuffer;
use Cache_Commons;
use Sys::Hostname;
my $host = hostname();




my $fork = 1;
my ($project_name,$chr, $chr_name, $step_name, $no_verbose, $annot_version);
GetOptions(
	'fork=s'     => \$fork,
	'project=s'  => \$project_name,
	'chr=s'      => \$chr,
	'step=s'     => \$step_name,
	'no_verbose=s' => \$no_verbose,
	'annot_version=s'    => \$annot_version,
);
warn "*_*_*_*_*_ $0 ".$host."*_*_*_*_*_".Dumper \@ARGV;
warn "$chr ";
unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
unless ($step_name) { confess("\n\nERROR: -step option missing... confess...\n\n"); }
unless ($chr) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }

my $buffer = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProject( -name => $project_name );
if ($annot_version) {
	$project->changeAnnotationVersion($annot_version);
}
my $hOutputs = Cache_Commons::output_files();
my $hSpecificCheck = Cache_Commons::specific_check_steps();

my $dir_log = $project->getCacheBitVectorDir()."/log/";
unless (-d $dir_log) { mkdir($dir_log); }

my @lChrNames;
if ($chr eq 'all') { @lChrNames = (1..22, 'X', 'Y', 'MT'); }
else { push(@lChrNames, $chr); }

foreach my $chr_name (@lChrNames) {
	if (-e "$dir_log/check_$step_name.$chr_name.ok") {
		`rm $dir_log/check_$step_name.$chr_name.ok`;
	}
	
	unless (exists $hOutputs->{'cache_'.$step_name}) {
		confess("\n\nERROR: '$step_name' step doesn't exists... confess.\n\n");
	}
	
	warn "\n### CHECK $step_name step\n" unless ($no_verbose);
	if ($step_name eq 'loh') {
		if (not $project->isSomaticStudy()) {
			my $cmd_touch = "date > $dir_log/check_$step_name.$chr_name.ok";
			system($cmd_touch);
			warn "-> NOT necessary (not a somatic study)\n" unless ($no_verbose);
			exit(0);
			if ($chr eq 'all') { next; }
			else { exit(0); }
		}
		else {
			my $cmd_touch = "date > $dir_log/check_$step_name.$chr_name.ok";
			system($cmd_touch);
		}
	}
	
	if ($step_name eq 'store_ids') {
		my ($hPos) = estim_nb_varids_from_vcf($project_name, $chr_name, $fork);
		my $estim_sub_vcf = $hPos->{$chr_name};
		my $freeze = $project->getCacheBitVectorDir().'/'.$hSpecificCheck->{'cache_store_ids'}->{'estim_nb_var'};
		$freeze =~ s/CHR_NAME/$chr_name/g;
		my $estim_sub_freeze = 0;
		my $hAllVar = retrieve $freeze;
		foreach my $var_id (keys %{$hAllVar}) {
			my @lFields = split('_', $var_id);
			if (length($lFields[2]) == 1 and length($lFields[3]) == 1) {
				$estim_sub_freeze++;
			}
		}
		if ($estim_sub_freeze >= $estim_sub_vcf) {
			warn "  -> Estimation nb substitution(s) in VCF $estim_sub_vcf and in FREEZE $estim_sub_freeze: OK\n" unless ($no_verbose);
		}
		else {
			warn "  -> Estimation nb substitution(s) in VCF $estim_sub_vcf and in FREEZE $estim_sub_freeze: ERROR\n" unless ($no_verbose);
			confess("\n\nERROR: missing variations...\nERROR: $step_name in error...\nERROR: confess\n\n");
		}
		if ($estim_sub_freeze == 0 and -e $project->getCacheBitVectorDir()."/lmdb_cache/$chr_name.empty") {
			my $cmd_touch = "date > $dir_log/check_$step_name.$chr_name.ok";
			warn '-> '.$cmd_touch unless ($no_verbose);
			`$cmd_touch`;
			my $cmd_chmod = "chmod 777 $dir_log/check_$step_name.$chr_name.ok";
			`$cmd_chmod`;
			if ($chr eq 'all') { next; }
			else { exit(0); }
		}
	}
	
	if ($step_name eq 'store_annotations') {
		my $freeze = $project->getCacheBitVectorDir().'/'.$hSpecificCheck->{'cache_store_ids'}->{'estim_nb_var'};
		$freeze =~ s/CHR_NAME/$chr_name/g;
		my $estim_sub_freeze = 0;
		my $estim_ins_freeze = 0;
		my $estim_del_freeze = 0;
		warn $freeze;
		my $hAllVar = retrieve $freeze;
		foreach my $var_id (keys %{$hAllVar}) {
			next if ($var_id =~ /-/);
			#next if ($var_id =~ /N/);
			my @lFields = split('_', $var_id);
			if (length($lFields[2]) == 1 and length($lFields[3]) == 1) {
				$estim_sub_freeze++;
			}
			elsif (length($lFields[2]) < length($lFields[3])) {
				$estim_ins_freeze++;
			}
			elsif (length($lFields[2]) > length($lFields[3])) {
				$estim_del_freeze++;
			}
			else {
				confess();
			}
		}
		system("date> ".$project->getCacheBitVectorDir()."/log/check_store_annotations.$chr_name.ok")  if  scalar keys %{$hAllVar} == 0;
		exit(0) if  scalar keys %{$hAllVar} == 0;
		
		my $buffer_cache  = new GBuffer;
		my $project_cache = $buffer_cache->newProjectCache( -name => $project_name, -cache => '1', -typeFilters => 'individual' );
		$project_cache->getPatients();
		my $chr = $project_cache->getChromosome($chr_name);
		
		my $var = $chr->getNewVector();
		
#		warn $var;
		if (exists $chr->global_categories->{substitution}) {
			$var += $chr->global_categories->{substitution};
		}
#		warn $var;
		my $estim_sub_obj = $chr->countThisVariants($var);
#		warn $estim_sub_obj;
		
		if ($estim_sub_freeze == $estim_sub_obj) {
			warn "  -> Nb substitution(s) in FREEZE $estim_sub_freeze and in CACHE $estim_sub_obj: OK\n" unless ($no_verbose);
		}
		else {
			warn "  -> Nb substitution(s) in FREEZE $estim_sub_freeze and in CACHE $estim_sub_obj: ERROR\n" unless ($no_verbose);
			
			confess("\n\nERROR: missing substitutions...\nERROR: $step_name in error...\nERROR: confess\n\n");
		}
		$var->Empty();
		if (exists $chr->global_categories->{insertion}) {
			$var += $chr->global_categories->{insertion};
		}
		if (exists $chr->global_categories->{large_duplication}) {
			$var += $chr->global_categories->{large_duplication};
		}
		my $estim_ins_obj = $chr->countThisVariants($var);
		if ($estim_ins_freeze <= $estim_ins_obj) {
			warn "  -> Nb insertion(s) in FREEZE $estim_ins_freeze and in CACHE $estim_ins_obj: OK\n" unless ($no_verbose);
		}
		else {
			warn "   -> Nb insertion(s) in FREEZE $estim_ins_freeze and in CACHE $estim_ins_obj: ERROR\n" unless ($no_verbose);
			confess("\n\nERROR: missing insertions...\nERROR: $step_name in error...\nERROR: confess\n\n");
		}
		$var->Empty();
		if (exists $chr->global_categories->{deletion}) {
			$var += $chr->global_categories->{deletion};
		}
		if (exists $chr->global_categories->{large_deletion}) {
			$var += $chr->global_categories->{large_deletion};
		}
		my $estim_del_obj = $chr->countThisVariants($var);
		if ($estim_del_freeze <= $estim_del_obj) {
			warn "  -> Nb deletion(s) in FREEZE $estim_del_freeze and in CACHE $estim_del_obj: OK\n" unless ($no_verbose);
		}
		else {
			warn "  -> Nb deletion(s) in FREEZE $estim_del_freeze and in CACHE $estim_del_obj: ERROR\n" unless ($no_verbose);
			#confess("\n\nERROR: missing deletions...\nERROR: $step_name in error...\nERROR: confess\n\n");
		}
		
		# Check si genes OK seulement sil y a des variants...
		my $total_var = $estim_sub_obj + $estim_ins_obj + $estim_del_obj;
		if ($total_var > 0) {
			my @lGeneIndexStored_lmdb_genes;
			foreach my $hGene (@{$chr->values_lmdb_genes()}) {
				push(@lGeneIndexStored_lmdb_genes, $hGene->{index_lmdb});
			}
			my @lGeneIndexStored_lmdb_genes_sorted = sort {$a <=> $b} @lGeneIndexStored_lmdb_genes;
			my $nb_index = scalar(@lGeneIndexStored_lmdb_genes_sorted);
			my $first_elem = $lGeneIndexStored_lmdb_genes_sorted[0];
			my $last_elem = $lGeneIndexStored_lmdb_genes_sorted[-1];
			if (scalar(@lGeneIndexStored_lmdb_genes_sorted) == $lGeneIndexStored_lmdb_genes_sorted[-1] + 1) {
				warn "  -> Nb Genes in LMDB_GENES = $nb_index (index $first_elem to $last_elem): OK\n" unless ($no_verbose);
			}
			else {
				warn "  -> Nb Genes in LMDB_GENES = $nb_index (index $first_elem to $last_elem): ERROR\n" unless ($no_verbose);
				confess("\n\nERROR: pb nb genes stored [NB]...\nERROR: $step_name in error...\nERROR: confess\n\n");
			}
			my $nb_index_ok = 0;
			my $no_genes = $chr->get_lmdb_genes("r");
			foreach my $index (@lGeneIndexStored_lmdb_genes_sorted) {
				my $res;
				eval { $res = $no_genes->get_index($index); };
				if ($@) {
					warn "\n\n  -> Problem genes indexation (ex with index $index)\n\n";
					confess("\n\nERROR: pb nb genes stored [INDEX]...\nERROR: $step_name in error...\nERROR: confess\n\n");
				}
				$nb_index_ok ++ if ($res); 
			}
			$no_genes->close();
			if (scalar(@lGeneIndexStored_lmdb_genes_sorted) == $nb_index_ok) {
				warn "  -> Nb Genes in LMDB_GENES = $nb_index (found $nb_index_ok index): OK\n" unless ($no_verbose);
			}
			else {
				warn "  -> Nb Genes in LMDB_GENES = $nb_index (found only $nb_index_ok index): ERROR\n" unless ($no_verbose);
				confess("\n\nERROR: pb nb genes stored [INDEX]...\nERROR: $step_name in error...\nERROR: confess\n\n");
			}
			
			my $hCommons = Cache_Commons::categories->{global};
			foreach my $type (keys %{$hCommons}){
				next if ($type ne 'variation_type');
				foreach my $cat (keys %{$hCommons->{$type}}) {
					#next if ($cat eq 'large_deletion');
					#next if ($cat eq 'large_duplication');
					unless (exists $chr->global_categories->{$cat}) {
						warn "\n\n";
						warn Dumper $chr->global_categories();
						confess("\n\nERROR: [CHR".$chr->id()."] global_categories->{$cat} not found. Pb cache ?? Die.\n\n");
					}
				}
			}
			warn "  -> All global_categories found: OK\n" unless ($no_verbose);
		}
	}
	
	
	if (-e $project->getCacheBitVectorDir()."/lmdb_cache/$chr_name.empty" and  $step_name ne 'global_infos') {
		warn "-> empty chromosome. SKIP\n" unless ($no_verbose);
		my $cmd_touch = "date > $dir_log/check_$step_name.$chr_name.ok";
		`$cmd_touch`;
		my $cmd_chmod = "chmod 777 $dir_log/check_$step_name.$chr_name.ok";
		`$cmd_chmod`;
		if ($chr eq 'all') { next; }
		else { exit(0); }
	}
	
	foreach my $key_dir (keys %{$hOutputs->{'cache_'.$step_name}}) {
		my $dir = $project->getCacheBitVectorDir().'/'.$key_dir;
		$dir =~ s/CHR_NAME/$chr_name/g;
		unless (-d $dir) {
			warn "-> $dir: ERROR\n" unless ($no_verbose);
			confess("\n\nERROR: directory $dir doesn't exist...\nERROR: $step_name in error...\nERROR: confess\n\n");
		}
		else {
			warn "-> $dir: OK\n" unless ($no_verbose);
		}
		foreach my $key_file (keys %{$hOutputs->{'cache_'.$step_name}->{$key_dir}}) {
			my $file = $dir.'/'.$key_file;
			$file =~ s/CHR_NAME/$chr_name/g;
			$file =~ s/\/\//\//g;
			unless (-e $file) {
				warn "  -> $file: ERROR\n" unless ($no_verbose);
				confess("\n\nERROR: file $file doesn't exist...\nERROR: $step_name in error...\nERROR: confess\n\n");
			}
			else {
				warn "  -> $file: OK\n" unless ($no_verbose);
			}
		}
	}
	
	my $cmd_touch = "date > $dir_log/check_$step_name.$chr_name.ok";
	warn '-> '.$cmd_touch unless ($no_verbose);
	`$cmd_touch`;
	my $cmd_chmod = "chmod 777 $dir_log/check_$step_name.$chr_name.ok";
	`$cmd_chmod`;
}


sub estim_nb_varids_from_vcf {
	my ($project_name, $chr_name, $fork) = @_;
	my @lChrNames;
	my $buffer = new GBuffer;
	$buffer->vmtouch(1);
	my $project = $buffer->newProject( -name => $project_name );
	my (@lVcf, $hPos);
	my @lPatients = @{$project->getPatients()};
	my $nb_pat = 0;
	my ($h, $hPos);
	foreach my $patient (@lPatients) {
		my @lVcf_Pat;
		foreach my $vcf (@{$patient->getVariationsFiles()}) {
			if (scalar(keys %{$patient->callingFiles()}) == 1) {
				push(@lVcf_Pat, $vcf);
			}
			elsif ($project->isDiagnostic()) { push(@lVcf_Pat, $vcf); }
			elsif (exists $patient->callingFiles()->{haplotypecaller} and $vcf =~ /haplotypecaller/) {
				push(@lVcf_Pat, $vcf);
			}
			elsif (not exists $patient->callingFiles()->{haplotypecaller} and $vcf =~ /unifiedgenotyper/) {
				push(@lVcf_Pat, $vcf);
			}
		}
		foreach my $vcf (@lVcf_Pat) {
			next unless (-e $vcf);
			next unless (-e $vcf.'.tbi');
			my $tabix = new Tabix(-data => $vcf);
			my $res = $tabix->query("chr$chr_name");
			eval { my $line = $tabix->read($res); };
			next if ($@);
			while(my $line = $tabix->read($res)){
				my (@lCol) = split(' ', $line);
				next if (length($lCol[3]) > 1);
				next if (length($lCol[4]) > 1);  
				next if ($line =~ /0\/0/);
				next if ($line =~ /\.\/\./);
				my $id = $chr_name.'_'.$lCol[1].'_'.$lCol[3].'_'.$lCol[4];
				unless (exists $hPos->{$id}) {
					$hPos->{$id} = undef;
					$h->{$chr_name}++;
				}
			}
		}
		$nb_pat++;
		last if ($nb_pat == 10);
	}
	return $h;
}