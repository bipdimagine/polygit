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
use lib "$RealBin/../../../../../../GenBo/lib/obj-nodb/";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use Getopt::Long;
use Carp;
use GBuffer;
use DBI;
use Compress::Snappy;
use Storable qw/thaw freeze/;
use GenBoNoSqlRocksGenome;
use File::Slurp qw(write_file);

my $fork = 1;
my ($project_name, $chr_name);
my $version;
my $chr;
GetOptions(
	'chr=s'		   => \$chr,
	'fork=s'		   => \$fork,
);
  my $tmp_dir = "/data-beegfs/tmp/djv/";
system("mkdir $tmp_dir") unless -e $tmp_dir;
warn "\n### Cache For Deja Vu\n";
my $buffer = new GBuffer;
my $ucsc_chr = "chr".$chr;
$ucsc_chr = "chrM" if $ucsc_chr eq "chrMT";


my $dir = $buffer->config->{deja_vu}->{path_rocks}."/HG19/".$buffer->config->{deja_vu}->{variations};
my $dir_out = $dir."/rocks/";
warn $dir_out;
my $dir38 = $buffer->config->{deja_vu}->{path_rocks}."/HG38/".$buffer->config->{deja_vu}->{variations}."/rocks/";
my $dir_put_38 = $dir38."/rocks/";
 warn $buffer->config->{deja_vu}->{path_rocks}."/XXXXX" . "/".$buffer->config->{deja_vu}->{variations};
warn  $dir38;
warn $chr;
warn "dbi:SQLite:dbname=$dir/$chr.dejavu.lite";
my $dbh = DBI->connect( "dbi:SQLite:dbname=$dir/$chr.dejavu.lite", "", "",{ sqlite_use_immediate_transaction => 0, } );
warn "$dir"."/1.dejavu.lite";
warn $dbh;
 my $sth = $dbh->prepare( 
      'SELECT _key,_value from __DATA__ '
      );			
 warn $sth;     
 my $rc = $sth->execute();
 my $rg = GenBoNoSqlRocksGenome->new(dir=>$dir_out,mode=>"c",index=>"genomic",chromosome=>$chr,genome=>"HG19",pack=>"",description=>[],pipeline=>1);
 my $rg38 = GenBoNoSqlRocksGenome->new(dir=>$dir38,mode=>"c",index=>"genomic",chromosome=>$chr,genome=>"HG38",pack=>"",description=>[],pipeline=>1);
 my $nb =0;
 warn $dir38;
   warn "- end -";
 #  $fp->write_batch();
   $rg->close();
   my $pm = new Parallel::ForkManager($fork);
    my $f4 = "$tmp_dir/"."$chr.later.bed";
    unlink $f4 if -e $f4;
    my %jobs;
   $pm->run_on_finish(
	sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hRes) = @_;
		unless (defined($hRes) or $exit_code > 0) {
			print qq|No message received from child process $exit_code $pid!\n|;
			die();
			return;
		}
		warn "write";
		 write_file($f4, {append => 1}, @{$hRes->{data}});
		warn "end write";
		delete $jobs{$hRes->{job}};
	}
);
   my $jobid= 0;
   $rg38->regions();
   foreach my $r (@{$rg->regions}) {
   	$jobid ++;
   	$jobs{$jobid} ++;
   	my $pid = $pm->start and next;
   	warn "start";
   	my $data = by_region($rg,$r,$rg38);
   	warn "end";
   	$pm->finish( 0, {data=>$data,job=>$jobid} );
	}
$pm->wait_all_children();
if (scalar(keys %jobs)){
	warn Dumper %jobs;
	die();
}

warn "---------------------------------------- STEP 2 --------------------------------------------";

system ("sort -u -k 2,2n $f4 > $f4.sort && bgzip -f $f4.sort && tabix -f -p bed $f4.sort.gz");

#unlink $f4;
warn "$RealBin/dejavu_rocks2.pl -chr=$chr -file=$f4.sort.gz";
system("$RealBin/dejavu_rocks2.pl -chr=$chr -file=$f4.sort.gz");
warn "$RealBin/dejavu_rocks2.pl -chr=$chr -file=$f4.sort.gz";
exit(0);
$rg38 = GenBoNoSqlRocksGenome->new(dir=>$dir38,mode=>"w",index=>"genomic",chromosome=>$chr,genome=>"HG38",pack=>"",description=>[]);

    my $v1 = Bio::DB::HTS::Tabix->new( filename => "$f4.sort.gz");
   foreach my $region (@{$rg38->regions}) {
   	warn "*";
	my $no38 =  $rg38->nosql_rocks($region);
	$no38->put_raw("date",time);
	my $start = $region->{start};
	my $end = $region->{end};
	warn $ucsc_chr.":".$start."-".$end;
	my $iter = $v1->query($ucsc_chr.":".$start."-".$end);
	
	if ($iter){
	my $xs;
	while (my $line = $iter->next) {
		chomp($line);
   		my ($chr,$pos38,$end,$id,$data) = split(" ",$line);
   		warn $pos38;
 		my ($a,$b,$c,$d) = split("_",$id);
   		my $id =  join("_",$a,$pos38,$c,$d);
   		$no38->put_raw($id,$data);
	}
	
	warn "close ==> ".$ucsc_chr.":".$start."-".$end;
 	warn "end ".$ucsc_chr.":".$start."-".$end;
   }
   }
   $rg38->close();
    exit(0);
   
   sub by_region {
   	my ($rg,$region,$rg38) = @_;
	
   	my $no = $rg->nosql_rocks($region);
   	my $no38 =  $rg38->nosql_rocks($region);
   	my $dbh = DBI->connect( "dbi:SQLite:dbname=$dir"."/$chr.dejavu.lite", "", "",{ sqlite_use_immediate_transaction => 0, } );
	my $start = $region->{start};
	my $end = $region->{end};
	#return if $end != 45000001;
	
	warn $start." ".$end;
 	my $sth = $dbh->prepare( 
      "SELECT _key,_value,start,end from __DATA__ where start > $start and end <= $end"
      );			
 my $rc = $sth->execute();
 my @pos;
  while (my @s = $sth->fetchrow()) {
  	warn $nb if $nb%1000000 ==0;
  	# $no->write_batch() if $nb%2000000 ==0; 
  
  	$nb ++;
  	  my $obj = thaw (decompress($s[1]));
  	  push(@pos,$ucsc_chr."\t".$s[2]."\t".$s[3]."\t".$s[0]."\t".$obj->{data}."\n");
  	  my $rid = $no->return_rocks_id_from_genbo_id($s[0]);
  	  next unless $rid;
  	 $no->put_batch_raw($rid,$obj->{data});
 
   }


   my $f1 = "$tmp_dir/"."$chr.$start.$end.bed";
   my $f2 = "$tmp_dir/"."$chr.$start.$end.out.bed";
   my $f3 = "$tmp_dir/"."$chr.$start.$end.error.bed";
  
   write_file("$tmp_dir/"."$chr.$start.$end.bed", @pos);
   warn "$tmp_dir/"."$chr.$start.$end.bed";
   system("/software/distrib/ucsc_util/liftOver $f1 /data-isilon/bipd-src/pnitschk/git2-polyweb/polygit/polymorphism-cgi/cache_nodb/scripts/rocks/hg19ToHg38.over.chain $f2 $f3");
   open (BED,"$f2");
  
   my @unsave;
   while(my $line = <BED>){
   	my $debug;
   	$debug =1 if $line =~ /42565845/;
   	chomp($line);
   	my ($chr,$pos38,$end,$id,$data) = split(" ",$line);
   	warn "coucou $chr,$pos38,$end,$id " if $debug;
	next if $chr ne $ucsc_chr;   	
	warn "1 " if $debug;
   	if ($pos38>$start && $pos38 <= $end){
   		my ($a,$b,$c,$d) = split("_",$id);
   		my $id =  join("_",$a,$pos38,$c,$d);
   		my $rockid = $no38->return_rocks_id_from_genbo_id($id);
   		$no38->put_batch_raw($rockid,$data);
   		warn "**********".$id." ".$d if $debug;
   	}
   	else {
   		push(@unsave,$line."\n");
   	}
   }
   	warn "END ".$start." ".$end;
   #	unlink $f1;
  # 	unlink $f2;
  # 	unlink $f3;
   	warn "$f1 unlink";
   $no38->write_batch;
   $no->write_batch;
  $no->close();
  $no38->close();
  return \@unsave;
   }