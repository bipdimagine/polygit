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
use lib "$RealBin/../../../../../../GenBo/lib/obj-nodb/";
use lib "$RealBin/../../../../../../GenBo/lib/obj-nodb/packages";
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
use liftOver;
use lib "$RealBin/../../utility";
use List::Util qw( shuffle);
use chunks;
use File::NFSLock qw(uncache);
use Fcntl qw(LOCK_EX LOCK_NB);
use Fcntl ':flock';
use Sereal qw(sereal_encode_with_object sereal_decode_with_object write_sereal_file read_sereal_file);
my $fork = 1;
my $project_name;
GetOptions(
	'fork=s' => \$fork,
	'project=s' => \$project_name,
);


my $dir_hg19_dv = "/data-isilon/DejaVu/HG19/variations/projects/";


system ("mkdir $dir38") unless -e $dir38;

 my $encoder = Sereal::Encoder->new({compress=>Sereal::SRL_ZSTD,compress_threshold=>0});
my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => "$project_name");
my $dir38 = $project->deja_vu_rocks_project_dir;
$project->getChromosomes();
my $hp;

my $lift = liftOver->new(project=>$project,version=>$project->lift_genome_version);
my $nb;



my ($chunks,$tree) = chunks::chunks_and_tree($project,$project->current_genome_version);
my ($chunks_lift,$tree_lift) = chunks::chunks_and_tree($project,$project->lift_genome_version);



		my $no = GenBoNoSql->new( dir => $dir_hg19_dv, mode => "r" );
		my $patients = $no->get($project_name, "patients");
		my $hpat ={};
		my $lt = time;
		foreach my $name (keys %$patients){
			$hpat->{$patients->{$name}} = $project->getPatient($name)->id;
			}
		$nb ++;
		foreach my $chr (@{$project->getChromosomes}){
			my $dd= "$dir38".$chr->name."/";
			system ("mkdir $dd") unless -e $dir38;
			
			my $h_hg19 = $no->get($project_name, $chr->name);
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
				$lift->add_region($vhh);
				#push(@$vh,$vhh)
			}
			
		}
		$no->close();
		warn "lift write ".abs($lt-time);
		my $nb1 = 0;
		my $t =time;
		my $lift = $lift->liftOver_regions();
		
		$t =time;

		my $ht = time;
		
		foreach my $chr ( keys %$lift){
			next unless $project->isChromosomeName($chr);
			warn "ok";
			my $chr_name = $project->getChromosome($chr)->name;
			my $current_regions = $chunks->{$chr_name}->[0];
			foreach my $vhh (@{$lift->{$chr}}) {
				if ($vhh->{HG38}->{start} >= $current_regions->{end}){
					my $oid = $current_regions->{id};
					 $current_regions = chunks::get_region($tree->{$chr_name},$vhh->{HG38}->{start});
					 die($oid." ". $vhh->{HG38}->{start}." ".$current_regions->{id} ) if $oid eq $current_regions->{id};
				}
				$nb1 ++;
				warn "$nb1 " if $nb1 % 500000 == 0;
				
				my $genoboid = $vhh->{HG38}->{chromosome}."_".$vhh->{HG38}->{start}."_".$vhh->{ref}."_".$vhh->{alt};
				my $rocksid = chunks::return_rocks_id($vhh->{HG38}->{start},$vhh->{ref},$vhh->{alt});
				push(@{$current_regions->{variants}},{id=>$rocksid,value=>$vhh->{value}});
			}
			}
			
		warn "end regions  ".abs(time-$ht);
		save_chunks($project,$chunks);

		warn $dir38."-";
	
		warn "end write :  ".abs(time -$t);
		exit(0) ;




#}
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
