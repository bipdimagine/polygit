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
use lib "$RealBin/../../../../../../../GenBo/lib/obj-nodb/";
use lib "$RealBin/../../../../../../../GenBo/lib/obj-nodb/packages";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use DBI;
use Compress::Snappy;
use Storable qw/thaw freeze/;
use GenBoNoSqlRocksGenome;
use GenBoNoSqlRocksAnnotation;
use File::Slurp qw(write_file);
use List::Util qw/shuffle/;
use lib "$RealBin/../../../utility";
use liftOverRegions;
use chunks;
use Sereal qw(sereal_encode_with_object sereal_decode_with_object write_sereal_file read_sereal_file);

my $fork = 1;
my $project_name;
GetOptions(
	'fork=s' => \$fork,
	'project=s' => \$project_name,
);







my $buffer = new GBuffer;

my $project = $buffer->newProjectCache( -name => "$project_name");
my $f2 = "/tmp/".$project->name.".lite";

$SIG{INT} = sub { cleanup_and_exit() };  # CTRL+C
$SIG{TERM} = sub { cleanup_and_exit() };
sub cleanup_and_exit {
    unlink $f2 if -e $f2; # Supprime le fichier si il existe
    exit 0;
}
my $hp;

my $nb;
my $exists_current = chunks::exists_project($project,$project->current_genome_version);
my $exists_lift  = chunks::exists_project($project,$project->lift_genome_version);
warn "ok ".$project->lift_genome_version;
#exit(0) if $exists_lift;

#exit(0) if chunks::exists_project($project,$project->current_genome_version);

#exit(0) if 	$project->current_genome_version ne "HG19";
 my $pm = new Parallel::ForkManager($fork);
 warn chunks::exists_project($project,$project->lift_genome_version);
 $project->getPatients();
 
my $lift = liftOverRegions->new(project=>$project,version=>$project->lift_genome_version);  
my $data_lift;	
$data_lift = HG38_HG19() if $project->current_genome_version eq "HG38"  && !($project->isGenome);
$data_lift = HG19_HG38() if $project->current_genome_version eq "HG19"; ;


my ($chunks_lift,$tree_lift) = chunks::chunks_and_tree($project,$project->lift_genome_version);
#### save liftover
exit(0) unless $chunks_lift;
 my $pm2 = new Parallel::ForkManager(5);
		my $ht = time;
	 
		foreach my $chr ( keys %$data_lift){
			next unless $project->isChromosomeName($chr);
			my $pid = $pm2->start and next;
#			warn $chr;
			my $chr_obj = $project->getChromosome($chr);
			my $chr_name = $chr_obj->name;
			my $current_regions = $chunks_lift->{$chr_name}->[0];
			
			foreach my $vhh (sort {$a->{LIFT}->{start} <=> $b->{LIFT}->{start}} @{$data_lift->{$chr}}) {
				
				if ($vhh->{LIFT}->{start} >= $current_regions->{end}){
					my $oid = $current_regions->{id};
					 $current_regions = chunks::get_region($tree_lift->{$chr_name},$vhh->{LIFT}->{start});
					 die($oid." : : ". $vhh->{LIFT}->{start}." :: ".$current_regions->{id} ) if $oid eq $current_regions->{id};
				}
				
				my $genoboid = $vhh->{LIFT}->{chromosome}."_".$vhh->{LIFT}->{start}."_".$vhh->{ref}."_".$vhh->{alt};
				my $rocksid = chunks::return_rocks_id_from_genbo_id($genoboid);
				push(@{$current_regions->{variants}},{id=>$rocksid,value=>$vhh->{value}});
			}
			chunks::save_chromosome_chunks($project,$chr_obj,$chunks_lift->{$chr_name},$project->lift_genome_version);
			$pm2->finish( 0, {data=>[],jobid=>0,snps=>[]} );
 			}
		 $pm2->wait_all_children();	
		chunks::save_final($project,$project->lift_genome_version);	
	






sub HG38_HG19 {

my $process;
 $pm->run_on_finish(
	sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hRes) = @_;
		
		unless (defined($hRes) or $exit_code > 0) {
			print qq|No message received from child process $exit_code $pid!\n|;
			return;
		}
		delete $process->{$hRes->{jobid}};
		foreach my $s (@{$hRes->{snps}}){
			$lift->add_region($s);
		}

	}
);

my ($chunks,$tree) = chunks::chunks_and_tree($project,$project->current_genome_version);

my $jobid = time;
my @chromosomes = shuffle (@{$project->getChromosomes});

foreach my $chr (@chromosomes){		
	#my $list =listVariants($chr);
	#$project->setListVariants($list);
	 $jobid ++;
	 $process->{$jobid};
	my $pid = $pm->start and next;
	my $chr_name = $chr->name;
	my $current_regions = $chunks->{$chr->name}->[0];
	
	my $no = $chr->get_rocks_variations("r");
#	warn $chr_name;
	my $snps =[];
	for (my $i =0;$i<$chr->size_vector();$i++) {
		my $v = $no->get_index( $i);
		$v->{buffer} = $buffer;
		$v->{project} = $project;
			if ($v->start >= $current_regions->{end}){
					my $oid = $current_regions->{id};
					 $current_regions = chunks::get_region($tree->{$chr_name},$v->start);
					 die($oid." ". $v->start." ".$current_regions->{id} ) if $oid eq $current_regions->{id};
				}
			my @he;
			my @ho;
			foreach my $p (@{$v->getPatients}){
				push(@ho,$p->id) if $v->isHomozygote($p);
				push(@he,$p->id) if $v->isHeterozygote($p);
			}
			my $c = compress1($project->id,\@he,\@ho);
			my $vhh;
			my ($c,$p,$a,$b) = split("_",$v->id);
			$vhh->{chromosome} = $chr->ucsc_name;
			$vhh->{start} = $p;
			$vhh->{end} = $p+length($b);
			$vhh->{rocksid} = $v->rocksdb_id;
			push(@$snps,$vhh);
			push(@{$current_regions->{variants}},{id=>$v->rocksdb_id,value=>$c});
	}
	
	chunks::save_chromosome_chunks($project,$chr,$chunks->{$chr->name},$project->current_genome_version);
	$pm->finish( 0, {data=>[],jobid=>$jobid,snps=>$snps} );
 }
 $pm->wait_all_children();
 
 return (undef,undef) if chunks::exists_project($project,$project->lift_genome_version);	
chunks::save_final($project,$project->current_genome_version);
return $lift->liftOver_regions();;
}

sub HG19_HG38 {
	 my $pm = new Parallel::ForkManager($fork);
	 my $process;
		 $pm->run_on_finish(
		sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hRes) = @_;
		
			unless (defined($hRes) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
			delete $process->{$hRes->{jobid}};
			foreach my $s (@{$hRes->{snps}}){
				$lift->add_region($s);
			}

		}
		);
		my $dir1 = $buffer->deja_vu_project_sqlite_dir("HG19","variations");
	
		my $f1 = "/data-isilon/DejaVu/HG19/variations/projects/".$project->name.".lite";

		my $no = GenBoNoSql->new( dir => $dir1 , mode => "r" );
		my $patients = $no->get($project_name, "patients");
		my $hpat ={};
		my $lt = time;
		my ($chunks,$tree) = chunks::chunks_and_tree($project,$project->current_genome_version);
		foreach my $name (keys %$patients){
			$hpat->{$patients->{$name}} = $project->getPatient($name)->id;
			}
		$nb ++;
		my @chromosomes = shuffle( @{$project->getChromosomes});
		my $jobid = time;
		foreach my $chr (@chromosomes){
			 $jobid ++;
	 		$process->{$jobid};
			my $pid = $pm->start and next;
			my $snps;
			warn $chr->name;
			my $h_hg19 = $no->get($project_name, $chr->name);
			my $variations;
			foreach my $vid (keys %{$h_hg19}){
				my $vh;
				my ($c,$p,$a,$b) = split("_",$vid);
				my $vhh;
				$vhh->{chromosome} = $chr->ucsc_name;
				$vhh->{start} = $p;
				$vhh->{end} = $p+length($b);
				$vhh->{end} = $p+length($b);
				
				my @hohe = split(" ",$h_hg19->{$vid});
				my @he = split(",",$hohe[0]); 
				my @ho = split(",",$hohe[1]); 
				
				my $list;
					#$list = join(",",@ho);
					my @uniq;
					my $huniq ={};
					my @nho;
					foreach my $a (@ho){
						my $id = $hpat->{$a};
						$huniq->{$id}++;
						push(@nho,$id);
					}
					my @nhe;
					foreach my $a (@he){
						my $id = $hpat->{$a};
						next if exists $huniq->{$id};
						push(@nhe,$id);
					}
				
				my $c = compress1($project->id,\@nhe,\@nho);
				$vhh->{value} = $c;
				$vhh->{he} = scalar(@he);
				$vhh->{ho} = scalar(@ho);
				$vhh->{ref} = $a;
				$vhh->{alt} = $b;
				push(@$variations,$vhh);
				push(@$snps,$vhh);
			}
			my $current_regions = $chunks->{$chr->name}->[0];
			warn "save";
			my $t =time;
			my @vorder = sort{$a->{start} <=> $b->{start}} @$variations;
			warn "order ".abs(time-$t);
			$t =time;
				foreach my $vhh (@vorder){
					if ($vhh->{start} >= $current_regions->{end}){
						my $oid = $current_regions->{id};
					 	$current_regions = chunks::get_region($tree->{$chr->name},$vhh->{start});
					 	die($oid." : : ". $vhh->{LIFT}->{start}." :: ".$current_regions->{id} ) if $oid eq $current_regions->{id};
					}
					my $genoboid = $chr->name."_".$vhh->{start}."_".$vhh->{ref}."_".$vhh->{alt};
					my $rocksid = chunks::return_rocks_id_from_genbo_id($genoboid);
					push(@{$current_regions->{variants}},{id=>$rocksid,value=>$vhh->{value}});
					
				}
			warn "before save :".abs(time-$t);		
				$t = time;
				chunks::save_chromosome_chunks($project,$chr,$chunks->{$chr->name},$project->current_genome_version);		
				warn "after save :".abs(time-$t);	
				$pm->finish( 0, {data=>[],jobid=>$jobid,snps=>$snps} );
				
		}#END CHROMOSOME
		$no->close();
		warn "end";
		 $pm->wait_all_children();
		chunks::save_final($project,$project->current_genome_version);	
		return $lift->liftOver_regions();
	
}



exit(0);

sub compress1 {
	my ($pid,$list1,$list2) =@_;
	my $compressed = pack("w w w w* w*",$pid, scalar(@$list1), scalar(@$list2), @$list1, @$list2);
	return $compressed;
}  
sub decompress1{
	my ($c1) =@_;
	my ($a,$b,@t) = unpack("w w w*",$c1);
	my @decompressed_list1 = @t[0..$a-1];
	my @decompressed_list2 = @t[$a..$a+$b-1];
	warn $a." ".$b." ".join(";",@decompressed_list1)." ++ ".join(";",@decompressed_list2);
}

END {
    unlink $f2 if -e $f2; # Supprime le fichier si il existe
}






