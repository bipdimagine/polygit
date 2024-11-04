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



my $fork = 1;
my ($project_name, $chr_name);
my $version;
GetOptions(
	'fork=s'       => \$fork,
	'project=s'    => \$project_name,
	'version=s'    => \$version,
);

unless ($project_name) { confess("\n\nERROR: -project option missing... confess...\n\n"); }


warn "\n### Cache For Deja Vu\n";
my $buffer1 = new GBuffer;
$buffer1->vmtouch(1);
my $project = $buffer1->newProjectCache( -name => $project_name, -verbose => 1 ,-version=>$version);
my $root_dir = $project->deja_vu_lite_dir() . "/projects/";
warn $root_dir;
#mkdir $root_dir unless -e $root_dir;
#unlink $root_dir . "/" . $project_name . ".lite" if -e $root_dir . "/" . $project_name . ".lite";
my @chr_names = map { $_->name } @{ $project->getChromosomes };
my @patient_names = sort { $a cmp $b } map { $_->name } @{ $project->getPatients };
my $dir_out = $project->rocks_cache_dir() . "/dejavu/";
mkdir $dir_out if not -d $dir_out;
warn $dir_out;
my $hpatients;
for ( my $i = 0 ; $i < @patient_names ; $i++ ) {
	$hpatients->{ $patient_names[$i] } = $i;
}


my $pm = new Parallel::ForkManager($fork);
$project->disconnect();
foreach my $chr (@{$project->getChromosomes} ){
	my $chr_name = $chr->name();
	my $pid = $pm->start and next;
#	next if $chr->id ne 21;
	my $hh;
	foreach my $variation (@{$chr->getVariations()}) {
		my $hv;
		my $array_patients;
		my $aho = [];
		my $ap  = [];
		foreach my $pat ( @{ $variation->getPatients() } ) {
			my $pn = $pat->name();
			my $hp;
			$hp->{name} = $pn;
			my $patient_id = $hpatients->{$pn};
			push( @$ap, $patient_id );
			push( @$aho, $patient_id ) if ( $variation->isHomozygote($pat) );
			$variation->{patients_details}->{$pn}->{vcf_infos} = $variation->{check_id};
			$variation->{patients_details}->{$pn}->{he} = $variation->{annex}->{ $pat->id() }->{he};
			$variation->{patients_details}->{$pn}->{ho} = $variation->{annex}->{ $pat->id() }->{ho};
			$variation->{patients_details}->{$pn}->{he_ho_details} = "he";
			$variation->{patients_details}->{$pn}->{he_ho_details} = "ho" if ( $variation->{annex}->{ $pat->id() }->{ho} eq '1' );
			$variation->{patients_details}->{$pn}->{he_ho_details} .= ':'.$variation->{annex}->{ $pat->id() }->{nb_all_ref}.':'.$variation->{annex}->{ $pat->id() }->{nb_all_mut};
			$variation->{patients_details}->{$pn}->{type} = "he";
			$variation->{patients_details}->{$pn}->{type} = "ho" if ( $variation->isHomozygote($pat) );
		}
		delete $variation->{buffer};
		delete $variation->{project};
		$variation->{heho_string} = join( ",", sort { $a <=> $b } @$ap );
		if ( scalar(@$aho) ) {
			$variation->{heho_string} = $variation->{heho_string}." ".join( ",", sort { $a <=> $b } @$aho ) . " HO";
		}
		$hh->{$variation->id} = $variation->{heho_string};
	}
	store( $hh, $dir_out."/$chr_name.dv.freeze" ) if $hh;
	warn $dir_out."/$chr_name.dv.freeze STORED!";
	
	$pm->finish();
}
$pm->wait_all_children();



exit(0);

my $chr = $project->getChromosome("Y");
#foreach my $chr (@{$project->getChromosomes} ){
	
	warn "!! start ".$chr->name;
	my $final_polyviewer = GenBoNoSqlRocks->new(dir=>$project->getRocksCacheDir("polyviewer_raw"),mode=>"r",name=>$chr->name);
		warn $project->rocks_pipeline_directory("polyviewer_raw");
		my $iter = $final_polyviewer->rocks->new_iterator->seek_to_first;
		my $nb = 0;
		warn "start chromosome ".$chr->name;
		while (my ($var_id, $value) = $iter->each) {
			$nb ++;
			my $c = $chr->name."!";
			next unless $var_id =~/^$c/;
			if  (ref ($var_id) =~ /HASH/) {
				warn Dumper $var_id;
				
			}
			warn $var_id." ".$nb if $nb%1000 == 0;
			my $v = $final_polyviewer->decode($value);
			warn Dumper $v;
			die();
		
	}
	$final_polyviewer->close();
#}

