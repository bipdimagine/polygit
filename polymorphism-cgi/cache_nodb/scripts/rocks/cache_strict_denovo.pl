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
use Net::FTP;
use Carp;
use Bio::DB::HTS;
use Devel::Cycle;
use GBuffer;
use Bit::Vector::Overload;
use Sys::Hostname;
use List::Util qw(min max sum);
use Parallel::ForkManager;



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
$project->getFamilies();
my $limit = 5;
if ($project->isGenome){
	$limit = 3;
}

if ($annot_version) {
	$project->changeAnnotationVersion($annot_version);
}
$project->getPatients();
my $chr = $project->getChromosome($chr_name);

my $hResults;
my $hErrors;

my $nbErrors = 0;
$project->preload_patients();

my $vector_denovo;
my $total_job;
my $vector_total = $chr->getNewVector();

warn 'here';
foreach my $patient (@{$project->getPatients()}) {
	 $vector_denovo->{$patient->name} =  $chr->getNewVector();
	 $vector_total += $patient->getVariantsVector($chr);
}

if ($vector_total->is_empty()) {
	warn 'empty chromosome';
	my $rocks4 = $chr->rocks_vector("w");
	foreach my $family (@{$project->getFamilies()}) {
		foreach my $children  (@{$family->getChildren}){
			my $vector_denovo = $chr->getNewVector();
			$rocks4->put_batch_vector_transmission($children,"ind_strict_denovo",$vector_denovo);
			warn "save";
		}
	}
	$rocks4->write_batch();
	$rocks4->close();
	
	system("date > $ok_file") if $ok_file;
	exit(0);
}


# AO-F2

### Generate BAM file
my $hfile;
my $tall = time;
my $fork_samtools = 10;
my $fm_fork = int($fork/$fork_samtools);
$fm_fork = 1 if $fm_fork == 0;
$project->getPatients();
$project->disconnect();
my $pm = new Parallel::ForkManager($fm_fork);

foreach my $family (@{$project->getFamilies()}) {
	$family->{tmp_dir}  = $project->getCallingPipelineDir("strict_denovo_".$family->name);
	my $tmp_dir = $family->{tmp_dir};
	foreach my $parent  (@{$family->getParents}){
		
		$hfile->{$parent->name} = "$tmp_dir/".$parent->name.".".$chr->name.".bam";
		next if -e $hfile->{$parent->name}.".bai";
		$project->disconnect();
		$pm->start() and next;
		warn "---> start ".$parent->name;
		my $cram =  $parent->getBamFile;
		my $chra = $chr->fasta_name();
		my $t =time;
		system("samtools view -T /data-isilon/public-data/genome/HG38_CNG/fasta/all.fa -b $cram $chra -@ $fork_samtools >".$hfile->{$parent->name});
	
		system("samtools index ".$hfile->{$parent->name}." -@ $fork_samtools");
		die($hfile->{$parent->name}) unless -e $hfile->{$parent->name}.".bai";
		warn "\t\t ++ ".$parent->name." ".abs(time-$t);
		$pm->finish(0, {});
	}
	
}
$pm->wait_all_children();

warn "--------------------------------------------------------";
warn "--- END BAM  ".abs(time-$tall);
warn "--------------------------------------------------------";

$tall = time ;
my $hbed ;
$pm = new Parallel::ForkManager($fork);
my $all_var;

my $no = $chr->flush_rocks_vector("r");
foreach my $family (@{$project->getFamilies()}) {
		my $tmp_dir = $family->{tmp_dir};
	foreach my $children  (@{$family->getChildren}){
		$hbed->{$children->id}= $tmp_dir."/".$children->name.'.'.$chr->name.".bed";
		#next if -e $hbed->{$children->id};
		#$pm->start() and next;
		warn "\n";
		warn $children->name;
		my $vector_denovo =  $no->get_vector_transmission($children,"ind_denovo");#$family->getVector_individual_denovo($chr,$children)->Clone();
		my @bits = $vector_denovo->Index_List_Read();
		open (BED , ">".$hbed->{$children->id});		
		my $nov = $project->getChromosome($chr_name)->get_rocks_variations("r");
		foreach my $vector_id (@bits){
				
				my $var = $nov->get_index($vector_id);
				
				$all_var->{$vector_id}->{start} = $var->start;
				$all_var->{$vector_id}->{end} = $var->end;
				 $all_var->{$vector_id}->{seq} = uc ($var->sequence);
				 
				print BED $chr->fasta_name."\t".($var->start-2)."\t".($var->end+2)."\n";
		}
		close BED;
		warn $hbed->{$children->id};
		#$pm->finish(0, {});
	}
}
$project->close_rocks();
$no->close;
$no = undef;

warn "----";
#$pm->wait_all_children();
warn "--------------------------------------------------------";
warn "--- END BED  ".abs(time-$tall);
warn "--------------------------------------------------------";
$tall = time ;
my $res_sambamba;
my $file_store;
my $hash_pileup;
 my $pm2 = new Parallel::ForkManager($fork);
#my $pm2 = new Parallel::ForkManager(10);
$pm2->run_on_finish (
	sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $res) = @_;
	
		warn "\t\t\t ====> END ".$res->{patient}." :: ".$res->{nb};
		
	}
);
my $nbp =0;
foreach my $family (@{$project->getFamilies()}) {
		my $tmp_dir = $family->{tmp_dir};
	foreach my $children  (@{$family->getChildren}){
		my $bed = $hbed->{$children->id};
		my $hbamba;
		foreach my $parent (@{$family->getParents()}) {
				my $file = $hfile->{$parent->name};
				$res_sambamba->{$parent->id} = $file.".".$children->name.".pileup";
				$nbp ++;
				$project->disconnect();
				$pm2->start() and next;
				warn "\t\t ---> start ".$parent->name;
			 	my $t =time;
			 	my $cmd = qq{samtools mpileup  $file -l $bed > }.$res_sambamba->{$parent->id};
			 	system($cmd);
			 	
			 	#`$cmd`;
				warn " \t\t  ++ pileup ".$parent->name." ".abs(time -$t)." ".$res_sambamba->{$parent->id};	
				
				$pm2->finish(0, {nb=>$nbp,patient=>$res_sambamba->{$parent->id}});
				
			}
	}
}
warn "wait";
$pm2->wait_all_children();

warn "--------------------------------------------------------";
warn "--- END PILEUP  ".abs(time-$tall);
warn "--------------------------------------------------------";
$tall = time ;
print join("\t","pos","A","C",,"G","T","del","ins")."\n";
foreach my $family (@{$project->getFamilies()}) {
		my $tmp_dir = $family->{tmp_dir};
	foreach my $children  (@{$family->getChildren}){
		my $bed = $hbed->{$children->id};
		my $hbamba;
		foreach my $parent (@{$family->getParents()})  {
	
			open (BAMBA,$res_sambamba->{$parent->id});
				while (my $res = <BAMBA>) {
					next if $res =~ /^RES/;
					my @tab = split(" ",$res);
				
					chomp($res);
					my $sequence = uc($tab[4]);
					my $count_A = () = $sequence =~ /A/g;
					my $count_T = () = $sequence =~ /T/g;
					my $count_C = () = $sequence =~ /C/g;
					my $count_G = () = $sequence =~ /G/g;
					my @deletions = ($sequence =~ /-(\d+)/g);
					my @insertions = ($sequence =~ /\+(\d+)/g);
					print join("\t",$tab[1],$count_A,$count_C,$count_G,$count_T,scalar(@deletions),scalar(@insertions))."\n";
					
					push(@{$hbamba->{$tab[1]}->{COV}}, $tab[3]);
					push(@{$hbamba->{$tab[1]}->{A}}, $count_A);
					push(@{$hbamba->{$tab[1]}->{C}}, $count_C);
					push(@{$hbamba->{$tab[1]}->{G} }, $count_G);
					push(@{$hbamba->{$tab[1]}->{T}} , $count_T);
					push(@{$hbamba->{$tab[1]}->{DEL}} , sum(@deletions)+0);
					push(@{$hbamba->{$tab[1]}->{INS}} , sum(@insertions)+0);
					warn Dumper $hbamba->{$tab[1]} if ($tab[1] == 61190108);
					warn Dumper $hbamba->{$tab[1]} if ($tab[1] == 61190108);
					warn Dumper $hbamba->{$tab[1]} if ($tab[1] == 61190108);
				}
			close (BAMBA);	
			
			
		}
		$hash_pileup->{$children->id} = $hbamba;
	
	}
}

warn "--------------------------------------------------------";
warn "--- END READ PILEUP  ".abs(time-$tall);
warn "--------------------------------------------------------";
$tall = time ;
#my $no = $chr->flush_rocks_vector("r");
my $rocks4 = $chr->rocks_vector("w");
warn $project;
foreach my $family (@{$project->getFamilies()}) {
	warn $family->name;
	warn $family->project;
	foreach my $children  (@{$family->getChildren}){
		warn "\t".$children->name();
		my $vector_denovo =  $rocks4->get_vector_transmission($children,"ind_denovo");#$family->getVector_individual_denovo($chr,$children)->Clone();
		my @bits = $vector_denovo->Index_List_Read();
		warn "\t".$children->name();
		my $vdenovo =construct_strict_denovo(\@bits,$children,$chr->getNewVector(),$hash_pileup->{$children->id},$chr);
		warn "\t end".$children->name();
		 $rocks4->put_batch_vector_transmission($children,"ind_strict_denovo",$vdenovo);
	}
	warn "denf ".$family->name;;
}

warn "--------------------------------------------------------";
warn "--- END   CONSTRUCT STRICT DENOVO ".abs(time-$tall);
warn "--------------------------------------------------------";
$rocks4->write_batch();
$rocks4->close();


#foreach my $family (@{$project->getFamilies}){
#	foreach my $child  (@{$family->getChildren}){
#		warn "**** ".$child->name." " .$vector_denovo->{$child->name}->Norm;
#		 $rocks4->put_batch_vector_transmission($child,"ind_strict_denovo",$vector_denovo->{$child->name});
#	}
#}



system("date > $ok_file") if $ok_file;
exit(0);
sub construct_strict_denovo {
	my ($bits , $children,$vdenovo,$hbamba,$chr) =@_;
	my $family = $children->getFamily();
		my $no = $project->getChromosome($chr_name)->get_rocks_variations("r");
			delete $no->{rocks};
				foreach my $vector_id (@$bits) {
				my $local_limit = $limit ;
				my $var = $no->get_index($vector_id);
				my $debug;
				$debug =1 if $var ->start  == 61190109;
				#next if $var->name =~/dup/;
				warn $var->name ;
				my $percent;
				my $r = $var->getRatio($children);
				if ($r > 40 ) {
					$local_limit = 5;
					$percent =0.05;
				}
				if ($var->getNbAlleleAlt($children) < 10){
					$local_limit = 3;
				}
				if ($var->isCnv && $var->isSrPr){
					warn "coucou";
					my $to_keep  = 0;
						foreach my $parent (@{$family->getParents()}) {
							warn "\t\t toto";
								 $to_keep += check_cnv($var,$children,$parent,$chr);	
									last  if ($to_keep == 0 );
								}
							$vdenovo->Bit_On($vector_id) if $to_keep == 2;
					warn "end";		
						}
					elsif ($var->isSrPr){
						my $to_keep = 0 ;
						
								foreach my $parent (@{$family->getParents()}) {
									 $to_keep += check_srpr($var,$children,$parent,$chr);	
									last  if ($to_keep == 0 );
								}
							$vdenovo->Bit_On($vector_id) if $to_keep == 2;
						}
				elsif ($var->isVariation) {
					my $alt = uc ($var->sequence);
					my $start = $var->start;
					my $min_cov = min($hbamba->{$start}->{COV}->[0],$hbamba->{$start}->{COV}->[1]);
					#next if $min_cov <5;
					
					my $nb_alt = $hbamba->{$start}->{$alt}->[0] + $hbamba->{$start}->{$alt}->[1];
					#warn $nb_alt." ".$min_cov." ".$var->name if $nb_alt >0;
					if ($nb_alt < $local_limit && $min_cov >= 5){
						$vdenovo->Bit_On($vector_id);
					}
				}
				elsif ($var->isInsertion) {
					my $start = $var->start;
					my $nb_alt = $hbamba->{$start}->{INS}->[0] + $hbamba->{$start}->{INS}->[1]+$hbamba->{$start+1}->{INS}->[0] + $hbamba->{$start+1}->{INS}->[1];
					# + $hbamba->{$start+2}->{INS}->[0] + $hbamba->{$start+2}->{INS}->[1];
					my $min_cov = min($hbamba->{$start}->{COV}->[0] + $hbamba->{$start}->{COV}->[1]);
					
					#if $nb_alt >0;
					$vdenovo->Bit_On($vector_id) if $nb_alt < $local_limit && $min_cov >= 5 ;
				}
				elsif ($var->isDeletion) {
					my $start = $var->start;
					my $nb_alt = $hbamba->{$start}->{DEL}->[0] + $hbamba->{$start}->{DEL}->[1]+$hbamba->{$start+1}->{DEL}->[0] + $hbamba->{$start+1}->{DEL}->[1]+$hbamba->{$start-1}->{DEL}->[0] + $hbamba->{$start-1}->{DEL}->[1];
				
					my $min_cov = max (min($hbamba->{$start-1}->{COV}->[0] + $hbamba->{$start-1}->{COV}->[1]),min($hbamba->{$start}->{COV}->[0] + $hbamba->{$start}->{COV}->[1]),min($hbamba->{$start+1}->{COV}->[0] + $hbamba->{$start+1}->{COV}->[1])) ;
						warn $nb_alt." ".$min_cov." limit ".$local_limit if $debug;
					warn Dumper $hbamba->{$start}->{DEL} if $debug;
					warn Dumper $hbamba->{$start-1}->{DEL} if $debug;
					$vdenovo->Bit_On($vector_id) if $nb_alt < $local_limit && $min_cov >= 5 ;
				}
				
				else {
					
				}
				
		}
		return $vdenovo;
	
	
}


sub check_cnv {
	my ($var,$children,$parent) = @_;
	my $debug =1;
	$debug =1 if $var->gnomad_id eq '17-78064045-del-187';
	$debug =1 if $var->name eq '17-78064045-del-187';
	warn "xxx";
	my $dp = $var->getNormDP($parent);
	my $dpc = $var->getNormDP($children);
	warn " dp parent ".$dp." ".$dpc if $debug == 1; 
	return 0 if $parent->meanDepth($chr->name,$var->start,$var->end) < 5;
	return 0 if $dp < 5; 
	warn "yyy";
	$dpc = 0.00000001 if $dpc == 0; 
	
	my $pc = int((($dpc-$dp)/$dpc)*100);
	$pc = abs($pc);
	return if $pc < 30; 
	
	warn $pc if $debug;
	if ($var->isDeletion) {
		return 0 if $pc > -30; 
		return 1;
		#return 1 if 
	}
	elsif ($var->isLargeDuplication) {
		return 1 if $pc > 30; 
		return 0 ;
		#return 1 if 
	}
	else {
		confess($var->name." ".$var->type)
	}
}

sub check_srpr {
	my ($var,$children,$parent,$chr) = @_;
	my @v1 =  $parent->sr_raw($chr,$var->start);
	return 0  if ($v1[0] < 5 or $v1[1]>5 or  $v1[2]>5);
	@v1 =  $parent->sr_raw($chr,$var->end);
	return 0 if ($v1[0] < 3 or $v1[1]>3 or  $v1[2]>3);
	return 1;
}


