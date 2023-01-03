#!/usr/bin/perl
use FindBin qw($Bin $RealBin);
use strict;

use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
use GBuffer ;
use Data::Dumper;
use Getopt::Long;
use Bio::DB::HTS;
use Bio::DB::HTS::VCF ;
use List::Util qw(max);
use GenBoNoSqlRocks;
### File Open ###



my $project_name;
my $fork;
my $patient_name;
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"patients=s" => \$patient_name,
);


 my $rocks = GenBoNoSqlRocks->new(
					dir         => ".",
					mode        => "w",
					is_index    => 1,
					name        => "test",
					is_compress => 1,
				);
				
					

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );

foreach my $patient (@{$project->getPatients}){
 my $sv_file =  $project->getVariationsDir("dragen-sv")."/".$patient->name.".sv.vcf.gz";
 next unless -e $sv_file;
 my $v = Bio::DB::HTS::VCF->new( filename => "$sv_file" );
 #my $v2 = Bio::DB::HTS::Tabix->new( filename => "$manta_file" );
	my $h = $v->header;
while (my $row = $v->next()){
	my $c = $row->chromosome($h) ;
	
	my $chr = $project->getChromosome($c);
	my $p = $row->position() ;
	my $z = $rocks->get($chr->name.":".$p);
	warn $z.":".$chr->name.":".$p;
	unless ($z){
		$rocks->put($chr->name.":".$p,1);
	}
	else {
		$z++;
		$rocks->put($chr->name.":".$p,$z);
	}
		
}
 
}