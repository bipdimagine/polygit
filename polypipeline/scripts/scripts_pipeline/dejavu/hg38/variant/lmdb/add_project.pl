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
my $exit = system("$RealBin/del_project.pl -project=$project_name -version=".$project->current_genome_version) if $exists_current;
die("$RealBin/del_project.pl -project=$project_name -version=".$project->current_genome_version) unless $exit == 0;
my $exists_lift  = chunks::exists_project($project,$project->lift_genome_version);
$exit = system("$RealBin/del_project.pl -project=$project_name -version=".$project->lift_genome_version) if $exists_lift;
die("$RealBin/del_project.pl -project=$project_name -version=".$project->lift_genome_version) unless $exit == 0;
warn "ok ".$project->lift_genome_version." ".$project->current_genome_version;
$project->isGenome();
#exit(0) if $exists_lift;

#exit(0) if chunks::exists_project($project,$project->current_genome_version);

#exit(0) if 	$project->current_genome_version ne "HG19";
 my $pm = new Parallel::ForkManager($fork);
 $project->getPatients();
 
 
 
my $lift = liftOverRegions->new(project=>$project,version=>$project->lift_genome_version);  
my $data_lift;	
$data_lift = HG38_HG19() if $project->current_genome_version eq "HG38"  ;
$data_lift = HG19_HG38() if $project->current_genome_version eq "HG19"; ;

#exit(0) if $project->isGenome && $project->current_genome_version eq "HG38";

my ($chunks_lift,$tree_lift) = chunks::chunks_and_tree($project->buffer,$project->lift_genome_version);
#### save liftover
exit(0) unless $chunks_lift;


 my $pm2 = new Parallel::ForkManager(5);
 
 my @final_regions;
 
 $pm2->run_on_finish(
	sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hRes) = @_;
		
		unless (defined($hRes) or $exit_code > 0) {
			print qq|No message received from child process $exit_code $pid!\n|;
			return;
		}
		foreach my $k (keys %{$hRes->{regions}}){
			push(@final_regions,$k) if $hRes->{regions};
		}
		

	}
);
 
 
 
		my $ht = time;
	 
		foreach my $chr ( keys %$data_lift){
			next unless $project->isChromosomeName($chr);
			my $pid = $pm2->start and next;
			my $chr_obj = $project->getChromosome($chr);
			warn "lift : ".$chr_obj->name;
			my $chr_name = $chr_obj->name;
			my $current_regions = $chunks_lift->{$chr_name}->[0];
			my $hregions ={};
			foreach my $vhh (sort {$a->{LIFT}->{start} <=> $b->{LIFT}->{start}} @{$data_lift->{$chr}}) {
				
				if ($vhh->{LIFT}->{start} >= $current_regions->{end}){
					my $oid = $current_regions->{id};
					 $current_regions = chunks::get_region($tree_lift->{$chr_name},$vhh->{LIFT}->{start});
					 die($oid." : : ". $vhh->{LIFT}->{start}." :: ".$current_regions->{id} ) if $oid eq $current_regions->{id};
				}
				
				my $genoboid = $vhh->{LIFT}->{chromosome}."_".$vhh->{LIFT}->{start}."_".$vhh->{ref}."_".$vhh->{alt};
				my $rocksid = chunks::return_rocks_id_from_genbo_id($genoboid);
				push(@{$current_regions->{variants}},{id=>$rocksid,value=>$vhh->{value}});
				$hregions->{$current_regions->{id}} ++ ;
			}
			chunks::save_chromosome_chunks_lmdb($project,$chr_obj,$chunks_lift->{$chr_name},$project->lift_genome_version);
			$pm2->finish( 0, {data=>[],jobid=>0,snps=>[],regions=>$hregions} );
 			}
		 $pm2->wait_all_children();	
	 warn Dumper @final_regions;
	chunks::save_index($project,$project->lift_genome_version,\@final_regions);





sub HG38_HG19 {
warn "here";
my $process;
my @final_regions;
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
		foreach my $k (keys %{$hRes->{regions}}){
			push(@final_regions,$k) if $hRes->{regions};
		}
		

	}
);

my ($chunks,$tree) = chunks::chunks_and_tree($project->buffer,$project->current_genome_version,);

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
	my $t1 =time;
	my $hho;
	my $hhe;
	my $arrayHo;
	my $arrayHe;
	foreach my $p (@{$project->getPatients}){
		my $v =  $p->getVectorOriginHo($chr);
		my $set = Set::IntSpan::Fast::XS->new($v->to_Enum());
		foreach my $a ($set->as_array){
			push(@{$arrayHo->[$a]},$p->id)
		}
		my $set = Set::IntSpan::Fast::XS->new($p->getVectorOriginHe($chr)->to_Enum());
		foreach my $a ($set->as_array){
			push(@{$arrayHe->[$a]},$p->id)
		}
	}
	my $hregions;
	for (my $i =0;$i<$chr->size_vector();$i++) {
		my $v = $no->get_index( $i);
		
			
		$v->{buffer} = $buffer;
		$v->{project} = $project;
			if ($v->start >= $current_regions->{end}){
					my $oid = $current_regions->{id};
					 $current_regions = chunks::get_region($tree->{$chr_name},$v->start);
					 die($oid." ". $v->start." ".$current_regions->{id} ) if $oid eq $current_regions->{id};
				}
			die($i ) if $arrayHe->[$i] == undef  && $arrayHo->[$i] == undef;	
			
			my $value = compress1($project->id,$arrayHe->[$i],$arrayHo->[$i]);
			my ($c,$p,$a,$b) = split("_",$v->id);
			my $vhh;
			$vhh->{chromosome} = $chr->ucsc_name;
			$vhh->{chromosome} = $chr->ucsc_name;
			$vhh->{start} = $p;
			$vhh->{end} = $p+length($b);
			$vhh->{rocksid} = $v->rocksdb_id;
			$vhh->{value} = $value;
			push(@$snps,$vhh);
			push(@{$current_regions->{variants}},{id=>$v->rocksdb_id,value=>$value});
			$hregions->{$current_regions->{id}} ++ ;

	}
	warn $chr->name." ::  end : ". abs (time - $t1);
	$t1 =time;
	chunks::save_chromosome_chunks_lmdb($project,$chr,$chunks->{$chr->name},$project->current_genome_version);
	warn "\t".$chr_name." save : ". abs (time - $t1);
	$pm->finish( 0, {data=>[],jobid=>$jobid,snps=>$snps,regions=>$hregions} );
 }
 $pm->wait_all_children();
 
 chunks::save_index($project,$project->current_genome_version,\@final_regions);
return $lift->liftOver_regions() ;
}



sub HG19_HG38 {
	 my $pm = new Parallel::ForkManager($fork);
	 my $process;
	 my $dir1 = $buffer->deja_vu_public_dir("HG19","projects_lite");
	my @final_regions;
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
			foreach my $k (keys %{$hRes->{regions}}){
				push(@final_regions,$k);
				}
		}
		);

		my $f1 = "$dir1".$project->name.".lite";

		my $no = GenBoNoSql->new( dir => $dir1 , mode => "r" );
		my $patients = $no->get($project_name, "patients");
		$no->close;
		my $hpat ={};
		my $lt = time;
		my ($chunks,$tree) = chunks::chunks_and_tree($project->buffer,$project->current_genome_version);
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
			my $t =time;
			my $no = GenBoNoSql->new( dir => $dir1 , mode => "r" );
			my $snps;
			my $h_hg19 = $no->get($project_name, $chr->name);
			my $variations;
			my $hregions;
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
			
			my @vorder = sort{$a->{start} <=> $b->{start}} @$variations;
		
				foreach my $vhh (@vorder){
					if ($vhh->{start} >= $current_regions->{end}){
						my $oid = $current_regions->{id};
					 	$current_regions = chunks::get_region($tree->{$chr->name},$vhh->{start});
					 	die($oid." : : ". $vhh->{LIFT}->{start}." :: ".$current_regions->{id} ) if $oid eq $current_regions->{id};
					}
					my $genoboid = $chr->name."_".$vhh->{start}."_".$vhh->{ref}."_".$vhh->{alt};
					my $rocksid = chunks::return_rocks_id_from_genbo_id($genoboid);
					push(@{$current_regions->{variants}},{id=>$rocksid,value=>$vhh->{value}});
					$hregions->{$current_regions->{id}} ++ ;
				}
				my $tt1 = abs(time-$t);
				$t = time;
				chunks::save_chromosome_chunks_lmdb($project,$chr,$chunks->{$chr->name},$project->current_genome_version);		
				warn "\t\t  :".$chr->name." ".$tt1." = ".abs(time-$t);	
				$no->close();
				$pm->finish( 0, {data=>[],jobid=>$jobid,snps=>$snps,regions=>$hregions} );
				
		}#END CHROMOSOME
		
		 $pm->wait_all_children();
 		chunks::save_index($project,$project->current_genome_version,\@final_regions);
		return $lift->liftOver_regions();
}


sub HG19_CHR_MCE {
	
}

exit(0);

sub compress1 {
	my ($pid,$list1,$list2) =@_;
	$list1 = [] if $list1 ==undef;
	$list2 = [] if $list2 ==undef;
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






