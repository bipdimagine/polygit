#!/usr/bin/perl
use FindBin qw($Bin);
use Carp;
use strict;
use JSON;
use Data::Dumper;
use CGI qw/:standard :html3/;
use Set::IntSpan::Fast::XS;
use Set::IntervalTree;
use List::Util qw[min max];
use FindBin qw($Bin);
use Storable qw(store retrieve freeze);
use Clone qw(clone);
use Parallel::ForkManager;
use strict;
use Text::CSV qw( csv );
use lib "$Bin/../../../../../GenBo/lib/obj-nodb/";
use Getopt::Long;
use lib "$Bin";
require  "$Bin/../SVParser.pm";
require  "$Bin/../parser/parse_pbsv.pm";
require  "$Bin/../parser/parse_hificnv.pm";
require  "$Bin/../parser/parse_sniffles2.pm";
require  "$Bin/../parser/parse_wisecondor.pm";
use lib "$Bin/../../dejavu/utility/";
use liftOverRegions;
use GBuffer;
use Text::CSV;
use Bio::DB::HTS::Tabix;
#
#je vais filtrer les cnvs en utilsant wisecondor le ratio 
#

my $fork = 5;
my $cgi = new CGI;


my $limit;
my $project_name;
my $patient_name;
my $fork =1 ;

GetOptions(
	'project=s' => \$project_name,
	'patient=s' => \$patient_name,
	'fork=s' => \$fork,
);
my $caller_type_flag = {
	"caller_sr" => 1,
	"caller_depth" => 2,
	"caller_coverage" => 4,
};
my $buffer = GBuffer->new();	
my $project = $buffer->newProject( -name => $project_name);
my $dir = $project->getCacheCNV(). "/rocks/";
my $rocks = GenBoNoSqlRocks->new(dir=>"$dir",mode=>"r",name=>"cnv");
my $parquet_file_quality = $project->getCacheCNV()."/".$project->name.".".$project->id.".cnv_quality.parquet";
my $parquet_file = $project->getCacheCNV()."/".$project->name.".".$project->id.".parquet";

my $dir_tmp = "/data-beegfs/tmp/";
my $filename = "$dir_tmp".$project->name.time.".quality.csv";
my $fh;
open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
	my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
	#my $col = ["id","patient_id","zscore","ratio","na","nr","pna","len","nb_caller","caller_type","dejavu"];
	my $col = ["id","type","chr","start","end","patient","len","nb_dejavu_patients","nb_dejavu_projects","genes","mz","mr","na","nr","pa","nb_caller","caller_sr","caller_depth","caller_coverage","blacklist"];

	$csv->print($fh, $col); 
	my $tab_pid;
foreach my $patient (@{$project->getPatients}){
	my $pid= $patient->id;
	
my $sql =qq{select * from '$parquet_file' where id <> 'Z' and patient=$pid ; };
warn $sql;
#my $sql =qq{select id from '$parquet_file' where  len > $minlength and patient=$patient_id; };
my $file = $patient->rawDataWiseCondor();
my $tabix = Bio::DB::HTS::Tabix->new( filename => $file );
my $cmd = qq{duckdb -json -c "$sql"};
my $res =`$cmd`;
my $array_ref = [];
$array_ref  = decode_json $res if $res;
my $bl_name = "spectre.blacklist.bed.gz";
$bl_name = "encode.blacklist.bed.gz";

my $tabix_blacklist = Bio::DB::HTS::Tabix->new( filename => "/data-pure/public-data/repository/HG38/blacklist/$bl_name"); 
	
#my $sth = $dbh->prepare($sql);
#$sth->execute();
my $nb = 0;
foreach my $row (@$array_ref) {
	$nb++;
	warn $nb."/".scalar(@$array_ref) if $nb %10 ==0;
	 my $cnv = $rocks->get($row->{id});
	  my @t =  keys %{$cnv->{score}};
	 #next  if @t > 1;
	 #next unless exists $cnv->{score}->{score_caller_sr};
	 my $l = $cnv->{len};
	 #if ($l > 10_000) {
		my $start = $cnv->{start};
		my $end = $cnv->{end};
		
		my $region = $cnv->{chromosome}.":$start-$end";
		my $res_blacklist = $tabix_blacklist->query("chr".$region);
		
		my $cnv_span = Set::IntSpan::Fast::XS->new("$start-$end");
		my $cnv_len = $end - $start + 1;
		my $bl_span = Set::IntSpan::Fast::XS->new();
		if ($res_blacklist){
			while (my $line = $res_blacklist->next) {
    		my ($b_chr, $b_start, $b_end) = split /\t/, $line;
    			$bl_span->add_range($b_start, $b_end);
			}
		}
		my $intersect = $cnv_span->intersection($bl_span);
		my $overlap_bp = scalar($intersect->as_array);   # taille totale des bases recouvertes
		my $overlap_pct = $overlap_bp / $cnv_len * 100;
		$row->{blacklist} = int($overlap_pct);
		#return $tree if not $res;
			my $s;
		my $sumr;
		my $sumz;
		my $nr =0;
		my $nz;
		
		my $na =0 ;
		
		my $res = $tabix->query( $region );
		if ($res){
		while ( my $line = $res->next ) {
			chomp($line);
		my @z = split(" ",$line);
		#warn "tabix $file $region;";
		my $a = $z[4];
		my $b = $z[5];
	#	warn $a." == ".$b;
	#	die();
		#	my ($a,$b) = split(" ",$p);
			$na ++ if $a eq "NaN";
			next if $a eq "NaN";
			$nr++;
			$sumr += $a;
			$sumz += $b;
			
		}
		}
		#warn "----";
		delete $cnv->{genes};
		my $mz = 0;
		$mz = $sumz/$nr if $nr > 0;
		my $mr = 0;
		$mr  = $sumr/$nr if $nr > 0;
		my $pa =0; 
		$pa = $na /($nr+$na) if ($nr+$na) > 0;
		#my $col = ["id","type","chr","start","end","patient","len","nb_dejavu_patients","nb_dejavu_projects","genes","mz","mr","na","nr","pa","nb_caller","caller_sr","caller_depth","caller_coverage"];
		$row->{mz} = $mz;
		$row->{mr} = $mr;
		$row->{na} = $na;
		$row->{nr} = $nr;
		$row->{pa} = $pa;
		$row->{nb_caller} = scalar(@t);
		$row->{caller_sr} = 0;
		$row->{caller_depth} = 0;
		$row->{caller_coverage} = 0;
		
		$row->{caller_sr} = 1 if exists $cnv->{score}->{score_caller_sr};
		$row->{caller_depth}= 1 if exists $cnv->{score}->{score_caller_depth};
		$row->{caller_coverage} = 1 if exists $cnv->{score}->{score_caller_coverage};
		
		my @t;
		foreach my $k (@$col){
			push(@t,$row->{$k});
			
		}
		#my $colv = [$cnv->{id},$patient->id,$mz,$mr,$na,$nr,$pa,$l,scalar(@t)];
		$csv->print($fh, \@t); 
	 #}
}


}
$fh->close();
	 my $query = "
	COPY (
        SELECT * from read_csv_auto(['$filename']) order by patient,type,chr,type,start
    )
    TO '$parquet_file_quality'  (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE);";
    print "\n# $query\n";
     
     system("duckdb :memory: -c \"$query \"");	 
	 
sub test_type {
	my($hash,$flag) = @_;
	$flag = "caller_".$flag;
	confess() unless $caller_type_flag->{$flag};
	return $hash->{caller_type_flag} & $caller_type_flag->{$flag};
}
