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
my $dir = $project->getCacheSV(). "/rocks/";
my $rocks = GenBoNoSqlRocks->new(dir=>"$dir",mode=>"r",name=>"sv");
my $parquet_file_quality = $project->getCacheSV()."/".$project->name.".".$project->id.".sv_quality.parquet";
my $parquet_file = $project->getCacheSV()."/".$project->name.".".$project->id.".parquet";

my $dir_tmp = "/data-beegfs/tmp/";
my $filename = "$dir_tmp".$project->name.time.".quality.csv";
my $fh;
my $bl_name = "spectre.blacklist.bed.gz";
$bl_name = "encode.blacklist.bed.gz";
warn "/data-pure/public-data/repository/HG38/blacklist/$bl_name";
my $tabix_blacklist = Bio::DB::HTS::Tabix->new( filename => "/data-pure/public-data/repository/HG38/blacklist/$bl_name"); 

open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
	my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
	my $col = ["id","type","chr1","pos1","chr2","pos2","patient","nb_dejavu_patients","nb_dejavu_projects","genes","dv1","dv2","bl1","bl2"];
	$csv->print($fh, $col); 
	
	my $tab_pid;
#my $tabix_blacklist = Bio::DB::HTS::Tabix->new( filename => "/data-pure/public-data/repository/HG38/blacklist/$bl_name"); 
my 	$tabix_inversion = Bio::DB::HTS::Tabix->new( filename => "/data-beegfs/sawfish/sawfish.inv.gz"); 
my 	$tabix_bnd = Bio::DB::HTS::Tabix->new( filename => "/data-beegfs/sawfish/sawfish.bnd.gz"); 
foreach my $patient (@{$project->getPatients}){
	my $pid= $patient->id;
	warn $patient->name;
my $sql =qq{select * from '$parquet_file' where id <> 'Z' and patient=$pid ; };
warn $sql;
#my $sql =qq{select id from '$parquet_file' where  len > $minlength and patient=$patient_id; };

my $cmd = qq{duckdb -json -c "$sql"};
my $res =`$cmd`;
my $array_ref = [];
$array_ref  = decode_json $res if $res;


#my $sth = $dbh->prepare($sql);
#$sth->execute();
my $nb = 0;
 my $bam = $patient->getAlignmentFile();;
foreach my $row (@$array_ref) {
	$nb++;
	warn $nb."/".scalar(@$array_ref) if $nb %100 ==0;
	 my $sv = $rocks->get($row->{id});
	 my $chr1 = $project->getChromosome($sv->{chrom1});
	  my $chr2 = $project->getChromosome($sv->{chrom2});
	 my $region =  $sv->{chrom1}.":".($sv->{pos1}-5)."-".($sv->{pos1}+5);
	 my $region2 =  $sv->{chrom2}.":".($sv->{pos2}-5)."-".($sv->{pos2}+5);
	 my $bl1 = blacklist( $sv->{chrom1},$sv->{pos1}-20,$sv->{pos1}+20,$tabix_blacklist);
	  my $bl2 = blacklist( $sv->{chrom2},$sv->{pos2}-20,$sv->{pos2}+20,$tabix_blacklist);
	  
	  #chr7:107417608-107418608
	 #  next unless $row->{pos1} == 13391895;
	 my $tabixf = $tabix_bnd;
	 my $a =0;
	 my $b=0;
	if ($sv->{type} eq "INV"){
		$tabixf = $tabix_inversion;
		 $a = tabix_inversion($tabixf, $sv->{chrom1}, $sv->{pos1},$sv->{pos2});
		  my $region =  $sv->{chrom1}.":".$sv->{pos1}."-".$sv->{pos2};
		  $bl1  = blacklist( $sv->{chrom1},$sv->{pos1},$sv->{pos2},$tabix_blacklist,1);
		  warn $bl1;
		  $bl2 = $bl2;
		 $b = $a;
		#next;
	 }
	 if ($a == 0) {
	# else {next;}
	  $a = tabix($tabixf,$region);
	  $b = tabix($tabixf,$region2);
	 }
	 warn $a." ".$b." ".$region." ".$region2 if $row->{pos1} == 13391895;
	 warn $bl1." ".$bl2  if $row->{pos1} == 13391895;;
	 $row->{dv1} = $a+0;
	 $row->{dv2} = $b+0;
	  $row->{bl1} = $bl1+0;
	 $row->{bl2} = $bl2+0;
		  my @t;
		foreach my $k (@$col){
			push(@t,$row->{$k});
			
		}
		$csv->print($fh, \@t); 
	
	 #my $col = [,"patient","nb_dejavu_patients","nb_dejavu_projects","genes","dv1","dv2"];

	#my $chr = $project->getChromosome($row->{chr1});
	#my $chr2 = $project->getChromosome($row->{chr2});
	#next if $row->{type} eq "INV";
	#my $res = {};
	#count_rnext($bam, $chr1,($row->{pos1}-50),$row->{pos1},$chr2,$res);
	#warn Dumper $res;
	#die();
	#count_rnext($bam,$chr2,$row->{pos2},($row->{pos2}+50),$chr,$res);
	#warn Dumper $res;
	#warn Dumper $sv;
	 #die();
}
warn "--->".$nb;
}

$fh->close();
	 my $query = "
	COPY (
        SELECT * from read_csv_auto(['$filename']) order by patient,type,nb_dejavu_patients
    )
    TO '$parquet_file_quality'  (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE);";
    print "\n# $query\n";
     
     system("duckdb :memory: -c \"$query \"");	 

sub count_rnext{
	my ($bam,$chr,$start,$end,$chr2  , $res ) =@_;
	my $region =  $chr->fasta_name.":".$start."-".$end;
	open(my $SAM, "samtools view $bam $region |") or die $!;
while (<$SAM>) {
    chomp;
    my @f = split /\t/;
    my ($qname, $flag, $rname, $pos, undef, $cigar,
        $rnext, $pnext, $tlen) = @f[0..8];
    # mate sur un autre chromosome
    warn $rnext." ".$chr2->fasta_name;
    if ($rnext eq '=' || $rnext eq $rname) {
		warn $rnext;
		$res->{ok} ++;
		
		}
		elsif ($rnext eq $chr2->fasta_name ){
			$res->{tl} ++;
		}
		
   	 else {
		$res->{bad} ++;
    }
}
	
}

sub tabix_inversion {
	my ($tabix,$chr1,$start1,$end1) = @_;
	my $res = $tabix->query("chr".$chr1.":".$start1."-".$end1);
	my $l = abs($start1-$end1);
	return 0 unless $res;
	my @val;
		while (my $line = $res->next) {
			my ($b_chr, $start2, $end21, $a,$b,$nb) = split /\t/, $line;
			my ($mate_chr, $end2) = $a =~ /([\w\.]+):(\d+)/;
			my $ov = min($end1,$end2) - max($start1,$start2) + 1;
			next if $ov <=0;
			my $p = $ov/$l;
			push(@val,$nb) if $p > 0.75;
			
			
			}
	return max(@val);
	
	



	
}

sub blacklist {
	my ($chr,$start,$end,$tabix,$debug) = @_;
		my $cnv_len = $end - $start + 1;
	my $cnv_span = Set::IntSpan::Fast::XS->new("$start-$end");
	my $res_blacklist = $tabix->query("chr".$chr.":".$start."-".$end) if $debug;
	return 0 unless $res_blacklist;

		my $bl_span = Set::IntSpan::Fast::XS->new();
		
		while (my $line = $res_blacklist->next) {
    		my ($b_chr, $b_start, $b_end) = split /\t/, $line;
    			$bl_span->add_range($b_start, $b_end);
		}
			
		my $intersect = $cnv_span->intersection($bl_span);
		my $overlap_bp = scalar($intersect->as_array);   # taille totale des bases recouvertes
		my $overlap_pct = int($overlap_bp / $cnv_len * 100);
	return $overlap_pct;
}

sub tabix {
	my ($tabixf,$region) = @_;
	my $res = $tabixf->query("chr".$region);
	my @val;
	if ($res){
			while (my $line = $res->next) {
    		my ($b_chr, $b_start, $b_end, $a,$b,$nb) = split /\t/, $line;
    			push(@val,$nb);
    			#$bl_span->add_range($b_start, $b_end);
			}
			return max @val;
		}
		return 0;
}
