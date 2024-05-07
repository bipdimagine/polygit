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
#use lib "/software/polyweb/poly-disk/poly-src/GenBo/lib/obj-nodb/";
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use GenBoNoSqlRocksGenome;
use GenBoNoSqlRocksVariation;
use List::Util qw(shuffle);
use Sys::Hostname;
use GenBoNoSqlRocksPolyviewerVariant;
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/packages/";
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/polyviewer/";
use Deep::Hash::Utils qw(reach slurp nest deepvalue);
use Carp;
use JSON::XS;

my ($project_name, $chr_name, $no_verbose, $skip_pseudo_autosomal,$version,$annot_version);
my $ok_file;
GetOptions(
	'project=s'    => \$project_name,
	'annot_version=s'    => \$annot_version,
	'no_verbose=s' => \$no_verbose,
	'skip_pseudo_autosomal=s' => \$skip_pseudo_autosomal,
	'version=s' => \$version,
'file=s' => \$ok_file,
);

 if ($ok_file && -e $ok_file) {
 	system("rm $ok_file");
 }
unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }

my $buffer = new GBuffer;


#my $color = $colors[ rand @colors ];
my $project = $buffer->newProject( -name => $project_name );
warn $project->rocks_directory();
my $final_polyviewer_all = GenBoNoSqlRocks->new(dir=>$project->rocks_directory(),mode=>"c",name=>"polyviewer_objects",pipeline=>1);
my $nb =0;
my $hpatients;
foreach my $patient (sort { $a->{name} cmp $b->{name} } @{ $project->getPatients }){
	$hpatients->{$patient->id} = $nb;
	$nb ++;
}


my @patient_names = sort { $a cmp $b } map { $_->name } @{ $project->getPatients };
my $hpatients2;
for ( my $i = 0 ; $i < @patient_names ; $i++ ) {
	$hpatients2->{ $patient_names[$i] } = $i;
}
my $root_dir = $project->rocks_directory(). "/deja_vu/";
mkdir $root_dir unless -e $root_dir;
#my $no = GenBoNoSql->new( dir => $root_dir, mode => "c" );
#$no->put( $project_name, "patients", $hpatients2 );
my $no_dv_rocks = GenBoNoSqlRocks->new( dir => $root_dir,name=>"dejavu", mode => "c" );
foreach my $chr (@{$project->getChromosomes} ){
	my $dv;
	my $final_polyviewer = GenBoNoSqlRocks->new(dir=>$project->rocks_pipeline_directory("polyviewer_raw"),mode=>"r",name=>$chr->name);
	my $iter = $final_polyviewer->rocks->new_iterator->seek_to_first;
		my $nb = 0;
		warn "start chromosome ".$chr->name;
		while (my ($var_id, $value) = $iter->each) {
			$nb ++;
			my $c = $chr->name."!";
			next unless $var_id =~/^$c/;
			if  (ref ($var_id) =~ /HASH/) {
				warn Dumper $var_id;
				die();
			}
			warn $var_id." ".$nb if $nb%10000 == 0;
			$final_polyviewer_all->put_batch_raw($var_id,$value);
			next;
			
			
			my $v = $final_polyviewer_all->decode($value);
			my $res;
			my $ap =[];
			my $aho =[];
			foreach my $pid (keys %{$v->{patients_calling}}){
				
				next unless exists $v->{patients_calling}->{$pid}->{array_text_calling};
				push(@{$res->[0]},$pid);
				push(@{$res->[1]},$pid)  if   lc($v->{patients_calling}->{$pid}->{gt}) eq "ho";;
				push( @$ap, $hpatients->{$pid} );
				push( @$aho, $hpatients->{$pid} ) if   lc($v->{patients_calling}->{$pid}->{gt}) eq "ho";
			}
			my $heho_string =  join( ",", sort { $a <=> $b } @$ap );
			if ( scalar(@$aho) ) {
				$heho_string .= " ".join( ",", sort { $a <=> $b } @$aho ) . " HO";
			}
			$dv->{$v->{id}} = $heho_string;
			#my $rocks_id = $no_dv_rocks->return_rocks_id_from_gnomad_id($v->{gnomad_id});
			$no_dv_rocks->put_batch_raw($v->{rocksdb_id},$heho_string);
			#$final_polyviewer_all->put_batch_raw($var_id,$value);
			
			
		#$finalrg->write_batch();
	}
	
	
	#$no->put( $project->name, $chr->name, $dv );
	$final_polyviewer->close();
	#$final_polyviewer_all->write_batch();
	warn "end chromosome ".$chr->name;
}
	warn "end all";
	my $t =time;
	#$final_polyviewer_all->write_batch();
	$final_polyviewer_all->close();
	warn "end all".abs(time - $t);
	#$no->close();
	warn $root_dir;
	system("date > $ok_file") if $ok_file;
	