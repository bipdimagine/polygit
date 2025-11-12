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
use lib "$RealBin/../";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime); 
use JSON::XS;
use Compress::Snappy;
use Getopt::Long;
use Carp; 
use GBuffer;
use Cache_Commons;
use Sys::Hostname;
use GenBoNoSqlRocksVector;
use Text::CSV;
use Devel::Size qw(size total_size);
use MCE::Loop;
use Storable qw(dclone);  # pour cloner le project par process si nécessaire
use Devel::Size qw(total_size);
use Storable qw(dclone);
 my $host = hostname();
use Scalar::Util qw(looks_like_number);
use lib "$RealBin/../../../GenBo/lib/obj-nodb/polyviewer/";
use PolyviewerVariant;
use File::Slurp;
use Digest::MD5 qw(md5 md5_hex md5_base64);
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone store);
use List::Util qw(sum);
use strict;
use warnings;
use RocksDB;
use GenBoNoSqlRocksTinyPolyviewerVariant;
# Ouvrir la DB avec options


# Forker plusieurs process qui écrivent en parallèle


warn "*_*_*_*_*_".$host."*_*_*_*_*_";

my $fork = 1;
my ($project_name, $chr_name, $annot_version);
my $ok_file;
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'annot_version=s'    => \$annot_version,
	'chr=s'        => \$chr_name,
	'file=s' => \$ok_file,
);

 if ($ok_file && -e $ok_file) {
 	system("rm $ok_file");
 }

my $parquets =[];



unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }


my $buffer  = new GBuffer;
$buffer->vmtouch(1);

my $project = $buffer->newProjectCache( -name => $project_name );

my $chunk_size = 3;
my @groups;
my @data = @{$project->getPatients};
while (@data) {
    push @groups, [ splice(@data, 0, $chunk_size) ];
}
my $ph;
my $gh;
foreach my $g (@groups){
	my $id =  join("-",map{$_->id} @$g);
	$gh->{$id}->{patients} = $g;
}

#
my  $tiny = GenBoNoSqlRocksTinyPolyviewerVariant->new(mode=>"c",pipeline=>1,project=>$project);
 $tiny->columns("");
my $data_final ;
my $chr = $project->getChromosome(1);
my @chroms =  @{$project->getChromosomes};
my $rocks_all;

foreach my $chr (@chroms){	
my $rocks;# = $tiny->rocksdb($chr);

#foreach my $id (keys %$gh){
#	 $rocks->{$id}->{tiny} = $tiny->rocksdb($chr,$id);
#	  $rocks->{$id}->{vector} = $chr->getNewVector();
#	  foreach my $p (@{$gh->{$id}->{patients}}) {
#	  	$rocks->{$id}->{vector} |= $p->getVectorOrigin($chr);
#	  }
#}
foreach my $f (@{$project->getFamilies}){
	$rocks->{$f->id} =  $tiny->rocksdb($chr,$f->id);
}

MCE::Loop::init {
    chunk_size => 'auto',  # chaque worker recevra un seul élément
    max_workers => 'auto', 
     gather => sub {
        my ($mce, $data) = @_;
        foreach my $id (keys %$data){
        	foreach my $fid (keys %{$data->{$id}}){
        		$rocks->{$fid}->put_batch_raw($id,$data->{$id}->{$fid});
        	}
		}
		foreach my $f (@{$project->getFamilies}){
 			$rocks->{$f->id}->write_batch();
		}
		#$rocks->write_batch();
    },
};

 my $no = $chr->get_rocks_variations("r");
my $size =$no->size;
$no->close();
$project->disconnect();
mce_loop {
    my ($mce, $chra,$chunk_id) = @_;
   warn $chunk_id;
   # my $chr = $chra->[0];
	my $hash = compute_chr($chr, $chra->[0], $chra->[-1]+1);
	 MCE->gather($chunk_id,$hash);
}(0..$size);
#$rocks->close();
MCE::Loop->finish;

warn "------->".$chr->name." end";
foreach my $f (@{$project->getFamilies}){
			warn $f->id;
 			$rocks->{$f->id}->close();
	}
#foreach my $id (keys %$gh){
#			 $rocks->{$id}->{tiny}->close();
#			 delete $rocks->{$id}->{vector};
#	}
}

warn $tiny->dir;
exit(0);




#warn "start  write ";
sub compute_chr {
	my ($chr,$start,$end) =@_;
  #  my $no = $chr->get_rocks_variations("r");
	my $res;
	my $final_polyviewer_all = GenBoNoSqlRocks->new(dir=>$project->rocks_directory(),mode=>"r",name=>"polyviewer_objects",cache=>1);

    my $nb = 0;
   	my $N ;
	my $hash_final;
    for (my $i = $start; $i < $end; $i++) {
        my $index= $chr->name."!".$i;
        my $vp = $final_polyviewer_all->get($index,1);
        
        warn $i unless $vp;
        next unless $vp;
		
		$hash_final->{$index} = $tiny->transform_polyviewer_variant($chr,$index,$vp,"encode");
		
#		$rocks_all->put($index,$array);
    }
    
 #   $rocks_all->close();
	$final_polyviewer_all->close();
	delete $project->{rocks};
	return $hash_final;
} 
