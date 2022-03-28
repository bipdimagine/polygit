#!/usr/bin/perl

use FindBin qw($RealBin);
use FindBin qw($Bin);
use FindBin qw($RealBin);
use lib "$RealBin/../GenBo";
use lib "$RealBin/../GenBo/lib/GenBoDB";
use lib "$RealBin/../GenBo/lib/obj-nodb";
use lib "$RealBin/../GenBo/lib/obj-nodb/packages";
use GenBoNoSql;
use lib "$RealBin/packages";
use Carp;
use strict;
use Data::Dumper;
use GBuffer;
use Getopt::Long;
use Carp;
use Digest::MD5 qw(md5 md5_hex md5_base64);
use Vcf;
use Storable qw/thaw freeze/;
use JSON::XS;
use Parallel::ForkManager;
use Logfile::Rotate;
use String::ProgressBar;
use File::Temp;
use Digest::MD5::File qw( file_md5_hex );
use File::stat;
use Devel::Size qw(size total_size);

#use NoSQL;
use Tie::LevelDB;
use Compress::Snappy;
use DBI;

#use  List::Compare;
use GenBoNoSqlText;
use GenBoNoSqlDejaVu;
use Storable qw/thaw freeze/;
use List::MoreUtils qw(natatime);
use Logfile::Rotate;
use File::Util;
use File::Temp qw/ tempfile tempdir /;
use Digest::MD5 qw(md5 md5_hex md5_base64);

#use RocksDB;
use GenBoNoSqlLmdb;
use LMDB_File qw(:flags :cursor_op :error);
use List::MoreUtils qw{ natatime uniq};
use List::Util qw(max);

my $buffer = GBuffer->new();

my $vquery = $buffer->validations_query();
my $zen    = $vquery->getAllValidations(1);
my $dir    = $buffer->config->{deja_vu}->{path};
my $no     = GenBoNoSqlDejaVu->new( dir => $dir, mode => "r" );
my $hprojects;
my $dir_public_data = $buffer->config->{deja_vu}->{path};

# $buffer->config->{'public_data'}->{root} . "/HG19/snp/deja_vu/lite/";
my $dir_projects = $dir . "/" . "projects/";

#2_37319312_C_T NGS2021_3816
#1_179526214_C_T
#ENSG00000116218_1!1_179526214_C_T
my $notodo = GenBoNoSql->new( dir => $dir_projects, mode => "r" );
foreach my $zid ( keys %$zen ) {

	#	next if $zid ne "ENSG00000116218_1!1_179526214_C_T";
	my ( $g, $vid ) = split( "!", $zid );
	my ( $chr, $pos, $b, $c ) = split( "_", $vid );
	my $k = $no->get( $chr, $vid );

	unless ($k) {
		next;
		my $pr      = $zen->{$zid}->[0]->{project_name};
		my $patient = $notodo->get( $pr, "patients" );
		my $vh      = $notodo->get( "$pr", $chr, 1 );

	}
	my $idacmg = $zen->{$zid}->[-1]->{idacmg};
	foreach my $l ( split( "!", $k ) ) {
		my ( $p, $nho, $nhe, $info ) = split( ":", $l );

		$hprojects->{$p}->{$chr}->{$vid} = $idacmg;

		#push(@{$hprojects->{$p}->{$chr}},$vid);
	}

	#next unless exists $vh->{$vid};
	#	my $no = $project->lite_deja_vu2();
	#my $k = $no->get($chr->name,$vid);
}
$notodo->close();
$no->close();

#my $chrs = [1..22,'X','Y','MT'];
#foreach my $c (@$chrs){
#warn "end 1";
my $fork         = 10;
my $pm           = new Parallel::ForkManager($fork);
my $novalidation = GenBoNoSqlLmdb->new(
	dir         => $dir,
	name        => "local.new",
	mode        => "c",
	is_compress => 1
);
$novalidation->put("date",time);
$novalidation->close();
my $hrun;
$pm->run_on_finish(
	sub {
		my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $h ) = @_;

		unless ( defined($h) or $exit_code > 0 ) {
			print qq|No message received from child process $exit_code $pid!\n|;
			die();
			return;
		}
		my $novalidation = GenBoNoSqlLmdb->new(
			dir         => $dir,
			name        => "local.new",
			mode        => "w",
			is_compress => 1
		);
		warn "end";
		#warn Dumper  $h->{data};
		delete $hrun->{$h->{run_id}};
		foreach my $s ( @{ $h->{data} } ) {
			$novalidation->put( $s->{id}, $s->{vector} );
		}
		$novalidation->close();
	}
);

my $nb = 0;
warn scalar( keys %$hprojects );

 $nb = int( scalar( keys %$hprojects ) / 50 + 1 );

my $iter = natatime($nb, (keys %$hprojects) );
my $id =0;
	while ( my @tmp = $iter->() ) {
		$id++;
		
		$hrun->{$id}++;
		my $pid = $pm->start and next;
		my @tres;
		foreach my $pr1 (@tmp) {
			#warn $pr1;
			my $buffer1 = GBuffer->new();
			#2_37319312_C_T NGS2021_3816
			my $project1 = $buffer1->newProjectCache( -name => "NGS20" . $pr1 );

			next unless $project1->genome_version() =~ /HG19/;
			$nb++;

			foreach my $chrn ( keys %{ $hprojects->{$pr1} } ) {
				my $chr    = $project1->getChromosome($chrn);
				my $vLocal = $chr->getNewVector();
				foreach my $vid ( keys %{ $hprojects->{$pr1}->{$chrn} } ) {

					next if $hprojects->{$pr1}->{$chrn}->{$vid} < 1;
					eval {
						my $k = $chr->cache_lmdb_variations->get($vid);
						next unless $k;
						$vLocal->Bit_On( $k->{index_lmdb} );
					};
					
				}
				push(@tres,{id=> $project1->name . "_" . $chrn,vector=>$vLocal});
				$chr->close_lmdb();
			}
			$project1 = undef;
			$buffer1  = undef;
		}
	#	warn "end";
		my $res;
		$res->{data} = \@tres;
		$res->{run_id} = $id;
		$pm->finish( 0, $res );
	}
	$pm->wait_all_children();
	die(Dumper $hrun) if keys %$hrun;
	#$novalidation->close();
	system("mv $dir/local.new $dir/local && chmod a+rw $dir/local");
	warn "$dir/local";

	#warn Dumper $hprojects;

