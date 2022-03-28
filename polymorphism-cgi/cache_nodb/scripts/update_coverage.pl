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
use lib "$RealBin/../../../GenBo/lib/obj-nodb/packages";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use Set::IntSpan::Fast::XS;
use image_coverage;
#use Cache_Commons;

my $fork = 1;
my $force;
my ($project_name, $chr_name,$annot_version);
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'chr=s'        => \$chr_name,
	'force=s'  => \$force,
	'annot_version=s'    => \$annot_version,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }
#unless ($chr_name) { confess("\n\nERROR: -chr option missing... confess...\n\n"); }
my $t = `hostname`;

my $nbErrors = 0;
my $buffer = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProjectCache( -name => $project_name);
if ($annot_version) {
	$project->changeAnnotationVersion($annot_version);
}

my $panels = $project->getPanels();
my $transcripts;
#foreach my $panel (@$panels){
#	foreach my $gene (@{$panel->getGenes}){
#		
#		next if $gene->getChromosome()->name ne $chr_name;
#		push(@$transcripts,@{$gene->getMainTranscripts});
#
#	}
#}
my $synonym = $project->liteAnnotations;

 my $z = $synonym->get_like("synonyms", "gene");
 my $annot4 =  GenBoNoSqlAnnotation->new(dir=>".",mode=>"c");
 

# warn Dumper $z;
 #die();
 $transcripts = [];
 
 my $genes = $project->newGenes([grep {$_=~/_$chr_name$/} values %$z]);
 	foreach my $gene (@$genes){
		
		next if $gene->getChromosome()->name ne $chr_name;
		push(@$transcripts,@{$gene->getMainTranscripts});

	}

my $images = uri_image($transcripts);
my $chr = $project->getChromosome("$chr_name");
my $no = $chr->lmdb_image_transcripts_uri("c");
$no->put("date",time);

foreach my $k (keys %$images){
	warn $k;
	$no->put($k,$images->{$k});
}

foreach my $k (keys %$images){
	#warn Dumper $no->get($k);
}


$no->close();
	
	sub uri_image {
	my ($transcripts) = @_;
	my $patients = $project->getPatients();
	my $fork = 10;
	my $nb = int(scalar(@$transcripts)/($fork*2))+1;
	warn $nb;
	my $genes ;
	foreach my $t (@$transcripts){
		my $gene = $t->getGene;
		next if exists $genes->{$gene->id};
		$genes->{$gene->id} = $gene;
	}

	#$transcripts = [values %$genes];
	
	
	my $pm = new Parallel::ForkManager($fork);
	my $iter = natatime $nb, @$transcripts;
	my @t_final;
	print qq{<div style="visibility: hidden">};
	my $images;
	$pm->run_on_finish(
    sub { 
    	my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$h)=@_;
  
    	unless (defined($h) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
		#	my $no = $chr->lmdb_image_transcripts_uri("w");
			
		foreach my $k (keys %{$h}){
			$images->{$k} = $h->{$k};
			
		}
		#$no->close();
    }
    );
  
	
  	$project->buffer->dbh_deconnect();
  	$|=1;
  	my $t =time;
  	
 	 while( my @tmp = $iter->() ){
 	 		my $pid = $pm->start and next;
 	 	
			$project->buffer->dbh_reconnect();
			my $himages ={};
			my $znb =0;
			my $dj;
 	 	foreach my $tr1  ( @tmp){ 
 	 		
 	 			$znb ++;
 	 			#my $gene = $tr1->getGene;
 	 			#next if exists $dj->{$gene->id};
 	 			 #$dj->{$gene->id} ++;
 	 			my $res;
 	 		if ($project->isNoSqlDepth){	
				 $res  = image_coverage::image_depth_lmdb ($patients, $tr1,0,0, 10, 10 );
			
 	 		}
 	 		else {
 	 			 $res  = image_coverage::image ($patients, $tr1,0,0, 10, 10 );
 	 		}
 	 		#next unless exists $res->{alert};
			my $uri = URI->new("data:");
			$uri->media_type("image/png");
			$uri->data($res->{image}->png);
			$himages->{$tr1->id}->{uri} = $uri;
			$himages->{$tr1->id}->{data} = $res->{data};
			$himages->{$tr1->id}->{alert} = $res->{alert} if exists $res->{alert} ;
 	 	}
 	 	
 	 	$pm->finish(0,$himages);
	}
	$pm->wait_all_children();
	$project->buffer->dbh_reconnect();
	print qq{</div>};
	return $images;
}
	
	