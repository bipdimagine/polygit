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
use GenBoNoSqlRocksTinyPolyviewerVariant;
use PolyviewerVariant;
use MCE::Shared;

my $htr = {
	 
                      'mane' => -99,
                      'end' => '-',
                      'spliceAI_cat' => '-',
                      'appris' => '-',
                      'codons_AA' => '-',
                      'name' => '-',
                      'spliceAI' => -99,
                      'dbscsnv' => -99,
                      'impact_score_text' => '-',
                      'consequence' => '-',
                      'polyphen' => -99,
                      'start' => '-',
                      'codons' => '-',
                      'alphamissense' => -99,
                      'sift' => -99,
                      'prot' => '-',
                      'nomenclature' => 'c.-1--07dupT',
                      'impact_score' => -99,
                      'cadd' => -99,
                      'nm' => "-",
                      'exon' => '-',
                      'revel' => -99,
                      'ccds' => "-",
                      'main' => -99
}; 


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
my $dir_tmp_cvs =  $project->getCallingPipelineDir($project->name.".parquet.".time);


##################
#construct_sql 
##################
$project->preload();

my  $tiny = GenBoNoSqlRocksTinyPolyviewerVariant->new(mode=>"c",pipeline=>1,project=>$project);
 $tiny->columns("");


#	
my $tume = time;
my $diro = $project->rocks_directory();
my $nb = 0;
foreach my $chr (@{$project->getChromosomes}) {
		my $rocks;
		my $error;
	my @batches;
	foreach my $f (@{$project->getFamilies}){
		$rocks->{$f->id} =  $tiny->rocksdb($chr,$f->id);
		warn $f->id;
	}
	#next if $chr->name ne "21";
	my $no =  GenBoNoSqlRocksVariation->new(dir=>$chr->project->rocks_directory("genbo"),mode=>"r",name=>$chr->name.".genbo.rocks");
   	my $nproc = 10;	
	my $size = $no->size;
 	$no->close;
 

$project->disconnect();
delete $project->{rocks};
foreach my $pat (@{$project->getPatients}){
	 $chr->rocks_vector("r")->cache_memory_patient($pat);
}
delete $chr->rocks_vector("r")->{rocks};
my $shared_hash = MCE::Shared->hash();
MCE::Loop::init {
    chunk_size => 'auto',  # chaque worker recevra un seul élément
    max_workers => 'auto', 
     gather => sub {
        my ($mce, $data,$a,$b) = @_;
        die unless exists $shared_hash->{$mce};
        delete $shared_hash->{$mce};
        foreach my $id (keys %$data){
        	
        	foreach my $fid (keys %{$data->{$id}}){
        		$rocks->{$fid}->put_batch_raw($id,$data->{$id}->{$fid});
        	}
		}
		foreach my $f (@{$project->getFamilies}){
 			$rocks->{$f->id}->write_batch();
		}
#		$rocks->write_batch();
    },
    on_post_exit => sub {
        my ($mce, $pid, $exit_code, $ident) = @_;
        if ($exit_code != 0) {
            warn "?? Worker $pid (ident=$ident) exited with error $exit_code\n";
            $error ++;
        } else {
            print "? Worker $pid (ident=$ident) exited normally\n";
        }
    }
};


mce_loop {
    my ($mce, $chra,$chunk_id) = @_;
   # my $chr = $chra->[0];
   	$shared_hash->{$chunk_id} ++;
	my $hash = compute($chr, $chra->[0], $chra->[-1]+1);
	#$project->disconnect(1);
	#delete $project->{rocks};

	 MCE->gather($chunk_id,$hash,$chra->[0], $chra->[-1]);
	 #$project->disconnect(1);
}(0..$size);
$project->disconnect();
delete $project->{rocks};
MCE::Loop->finish;
confess()  if %$shared_hash;
warn "end chromosome ".$chr->name;
foreach my $f (@{$project->getFamilies}){
 			$rocks->{$f->id}->close();
}	
$project->disconnect();
delete $project->{rocks};
$rocks = undef;
}
exit(0);


sub compute {
	my ($chr,$start,$end) =@_;
    my $no = GenBoNoSqlRocksVariation->new(dir=>$chr->project->rocks_directory("genbo"),mode=>"r",name=>$chr->name.".genbo.rocks");
# $chr->get_rocks_variations("r");
    	my $hash_final;
    for (my $i = $start; $i < $end; $i++) {
        my $variation = $no->get_index($i);
        next unless $variation;
        $variation->{buffer}  = $buffer;
        $variation->{project} = $project;
         my $index= $chr->name."!".$i;
       my $vp =  PolyviewerVariant->new();
		$vp->setLmdbVariant($variation);
		$vp->{hgenes} = {};
		$vp->{genes_id} = [];
		my $code =0;
		foreach my $g (@{$variation->getGenes}){
			my $h = $vp->set_gene($variation,$g);
			$h->{code} = $code;
			$vp->{hgenes}->{$g->{id}} = $h;
			
			push(@{$vp->{genes_id}},$g->{id});
			$code ++;
		}
		##############
		#	next;
		##############0
		$vp->{hpatients} ={};
		$vp->{patients_id} = [];
		my $dvp;
		
		
##		warn " before setPatients - time: ".abs(time - $tloop) if $nb %5000 == 0; 
		foreach my $pat (@{$variation->getPatients}){
			
		foreach my $p (@{$pat->getFamily()->getMembers}){
#				
				next if exists $dvp->{$p->id};
				$dvp->{$p->id} ++;
				my $h = $vp->set_patient_cache($variation,$p);
				$vp->{patients_calling}->{$p->id} =$h; 
				#$chr->rocks_vector("r")->get_vector_transmission($p,"ind_recessive") if $p->isChild; 
				#warn  $variation->getTransmissionModelType($p->getFamily(),$p);
				$vp->{patients_calling}->{$p->id}->{model} = $variation->getTransmissionModelType($p->getFamily(),$p);
			}
		}
		
		$hash_final->{$index} = $tiny->transform_polyviewer_variant($chr,$index,$vp,"encode");
	
	
     
	
 	delete $variation->{project};
 	delete $variation->{buffer};
 	$variation = undef;
			
    }
   # $chr->rocks_vector("close");
    #->close();
	$no->close();

	return $hash_final;
   # warn "FIN ".$chr->name." ".$region->{start}."-".$region->{end};

	
} 


