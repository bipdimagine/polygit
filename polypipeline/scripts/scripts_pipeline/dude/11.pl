#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
#use Set::IntSpan;
use GenBoNoSqlLmdb;

use Carp;
use strict;
use Set::IntSpan::Fast::XS ;
use Data::Dumper;
use GBuffer;
use Getopt::Long;
use Carp;
 use JSON::XS;
 use List::MoreUtils qw(natatime uniq);
 use Tabix;
#use GenBoNoSql;
#use GenBoNoSqlDejaVu;
#use Array::Diff;
#use UnQLite;
 # use Tie::LevelDB; 
 # use Devel::Size qw(size total_size);
  #use Devel::Size::Report qw/report_size/;
  #use Statistics::Descriptive;
  

 $| =1;               
#ENSG00000234585
my $buffer = new GBuffer;
my $project_name= "NGS2017_1534";
my $fork;
GetOptions(
	'project=s' => \$project_name,
	'fork=s' => \$fork,
);
die("hey man,  no fork ") unless $fork;

my $project = $buffer->newProject( -name 			=> $project_name );

my $f1 = $project->getCacheDir() . "/coverage_lite/primers.lite";
warn $f1;
unlink $f1  if -e $f1;


my $pm2 = new Parallel::ForkManager($fork);

	my $total;
	$pm2->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data ) =
			  @_;
			  warn "end".$data->{chr};
			$total->{ $data->{chr} } = $data->{primers};
		}
	);
	my $dir_out   =$project->noSqlCnvsDir;
	foreach my $chr ( @{ $project->getChromosomes } ) {
		#next if $chr->name() ne "4";
		my $chr_name = $chr->name();
		
		my $pid      = $pm2->start() and next;
		my $buffer   = new GBuffer;
		my $project =
		  $buffer->newProject( -name => $project_name, -verbose => 1 );
		my $chr1    = $project->getChromosome($chr_name);
		warn $chr_name." read";
		my $primers = $chr1->getPrimers();
			warn $chr_name." end ".scalar(@$primers);
	#	map { delete $_->{project}; delete $_->{buffer} } @$primers;
	#	my $dir_out   =$project->getCacheDir() . "/cnv_lite";
		my $no2 = GenBoNoSqlLmdb->new(dir=>$dir_out,mode=>"c",name=>$chr->name,is_compress=>1);
			warn $chr_name." create ".scalar(@$primers)." ".$dir_out;
		$no2->put("toto","tutu");
		warn $chr_name." put ".scalar(@$primers);
		my $captures;
		my $lists;
		warn "--->".scalar(@$primers);
		my $nbp =0;
		foreach my $primer (@$primers){
			warn $nbp."/". scalar(@$primers) if $nbp %100==0;
			$nbp++;
		#warn $primer->start()." ".$primer->end;
			foreach my $c (@{$primer->getCaptures}){
				push(@{$lists->{"capture_".$c->id}},[$primer->id,$primer->start]);
			}
			
			foreach my $c (@{$primer->getRuns}){
				push(@{$lists->{"run_".$c->id}},[$primer->id,$primer->start]);
			}
				foreach my $c (@{$primer->getPatients}){
				push(@{$lists->{"patient_".$c->id}},[$primer->id,$primer->start]);
			}
				delete $primer->{project}; delete $primer->{buffer};
			$no2->put($primer->id,$primer);
			}
		map {$no2->put($_,$lists->{$_})} keys %$lists;
			
		$no2->close();
		my $res;
		$res->{primers} = $primers;
		$res->{chr}     = $chr->name;
		
		$pm2->finish( 0, $res );
	}
	$pm2->wait_all_children;
	
	my $no  = $project->noSqlCnvs("c");
	$no->put_bulk( "primers", $total );
	$no->close();
 my $file = "$dir_out/primers.bed";
 warn $file;
 my $nb = 0;
 open (BED,">$file") or die();
 
 foreach my $chr ( @{ $project->getChromosomes } ) {
 	
 	my @primers = sort{$a->{start} <=> $b->{start}} @{$total->{$chr->name}} if exists $total->{$chr->name};
 	foreach my $p (@primers){
 		$nb ++;
 		print BED $chr->ucsc_name."\t".$p->start."\t".$p->end."\t".$p->{id}."\n";
 	}
 }
 close(BED);
 my $tabix = $buffer->software("tabix");
 my $bgzip  = $buffer->software("bgzip");
 unlink $file.".gz" if -e $file.".gz";
 system("$bgzip -f $file && tabix -f -p bed $file.gz ");
 die() unless -e $file.".gz.tbi";
 warn "nb de primers = ".$nb;
 exit(0);
 
 
 
 
 