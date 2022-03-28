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

my $fork;
my $force;
my $verbose = 1;
my $chr;
GetOptions(
	'fork=i'    => \$fork,
	'verbose=i' => \$verbose,
	'force=i'   => \$force,
	'chr=s'     => \$chr,
);
$| = 1;
die() unless $chr;
 my @he =(1,2,5,7);
 my @ho =(9,11);
				my $nbhe = scalar(@he);
				my $nbho = scalar(@ho);
				my $test = pack("C* x C*",@he,chr(0),chr(0),chr(0),@ho);
				my ($a,$b) = split(chr(0).chr(0).chr(0),$test);
				my (@tz1) = unpack("(C4)",$a);
				my (@tz2) = unpack("(C2)",$b);
				#my (@tz1,$a,@tz2) = unpack("(C4) x (C2)",$test);
				warn Dumper @tz1;
				warn Dumper @tz2;

my $buffer    = GBuffer->new();

my $db_final;
$SIG{'INT'} = sub {
	print "Caught One!\n";
	$db_final->close if $db_final;
	exit(0);
	die();
};
my $dir_prod = "/data-isilon/tmp/dejavu/";

my %exclude;
$exclude{NGS2014_0001} = 1;
sleep(1);

#$buffer->config->{public_data}->{HG19} = "/data-isilon/DejaVu/HG19/snp/";

my $dir_public_data =  $buffer->config->{deja_vu}->{path};
 # $buffer->config->{'public_data'}->{root} . "/HG19/snp/deja_vu/lite/";
my $dir_projects = $dir_public_data . "/" . "projects/";
unless ( -d $dir_projects ) {
	confess($dir_projects);
}


#
my $no_reload = undef;
#my $dir_tmp = tempdir ( DIR => "/data-beegfs/tmp/",CLEANUP => 1 );
my $dir_tmp = "/tmp/pipeline";
system("mkdir -p $dir_tmp && chmod a+rwx $dir_tmp") unless -e $dir_tmp;

my @projects = grep { !( exists $exclude{$_} ) } @{ $buffer->listProjectsForDejaVu() };
  my $id = time;
  my $dbfile = "$dir_tmp/$chr.$id.sqlite";
  create_projects_sqlite() if $chr == 1;
create_first_sqlite($dbfile) unless -e $dbfile;
#exit(0);
create_dejavu_lite($dbfile);






exit(0);

sub create_projects_sqlite {
	
	my $pr = String::ProgressBar->new( max => scalar(@projects) );

my $noprojects = GenBoNoSqlDejaVu->new( dir => $dir_prod, mode => "c" );
warn $dir_prod;

foreach my $project_name ( sort @projects ) {
	next unless -e $dir_projects . "/$project_name.lite";
	next if -z $dir_projects . "/$project_name.lite";
	my $patient;
	my %h2;
	eval {
	my $notodo  = GenBoNoSql->new( dir => $dir_projects, mode => "r" );
	 $patient = $notodo->get( $project_name, "patients" );
	next unless $patient;
	 %h2 = reverse %$patient;
	$noprojects->put( "projects", $project_name, \%h2 );
	};
	#if ($@) {
	#	$noprojects->put( "projects", $project_name, \%h2 );
	#}

}

warn "first step";
$noprojects->close();
#die();	
}


sub create_first_sqlite {
		my ($dbfile) = @_;
	my $pm = new Parallel::ForkManager($fork);
	warn warn "---------------------------------------\n";
	warn "--------------Compute lite database ---------------------\n";
	warn "---------------------------------------\n";

	my $hash;
my $nbp = scalar(@projects);

#last;
warn "---------------------------------------\n";
$| = 1;

my $id1 = time;
my $running_jobs;
my $jid = time;

my $htotal;
my @files;

my $dbh = DBI->connect( "dbi:SQLite:dbname=$dbfile", "", "",
		{ sqlite_use_immediate_transaction => 0, } );
		$dbh->do("PRAGMA journal_mode = MEMORY") or die $DBI::errstr;;
		$dbh->do("PRAGMA synchronous = OFF") or die $DBI::errstr;;
	$dbh->do("create table tabLines(idx varchar(10000) , start INTEGER, end INTEGER,ho INTEGER,line varchar(5000), variation_type varchar(2));")
	  or die $DBI::errstr;
	 $dbh->disconnect;
	 
$pm->run_on_finish(
    sub { 
    	my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$h)=@_;
  
    	unless (defined($h) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
			delete $running_jobs->{$h->{jobid}};
    		my $c = $h->{chr};
    		push(@files,$h->{file});
    		my $f = $h->{file};
    		my $t = time;
    		
    		system(qq{sqlite3 $dbfile ".import $f tablines "});
    		warn "\t end import $dbfile $f ".abs(time -$t);
			unlink $f;
    	
    }
    
    );


	my $nbj =-1;
	my $nfork = 20;
	my $nb = int((scalar(@projects)+1)/($nfork));
	
	$nb = 20 if $nb == 0;
	my $iter = natatime($nb, @projects);
	my $tloop = time;
	 while( my @tmp = $iter->() ){
	 	$jid++;
		$running_jobs->{$jid} =1;
		$nbj ++;
		my $pid = $pm->start and next;
		my $tpart = time;
		my $h;
		my $text_file = $dir_tmp."/$chr.$nbj.$id.txt";
		open (TMP, ">".$text_file);
		foreach my $project_name (  @tmp) {
			next unless -e $dir_projects . "/$project_name.lite";
		
			my $notodo = GenBoNoSql->new( dir => $dir_projects, mode => "r" );
			my $vh = $notodo->get( "$project_name", $chr, 1 );
			$notodo->close();
			my ($y,$n) = split("_",$project_name);
			$y =~s/NGS//;
			
			foreach my $vid ( keys %$vh ) {
				next if $vid =~ /\|/;
				my ($chr,$pos,$a,$b) = split("_",$vid);
				my $line = $project_name . ";" . $vh->{$vid};
				my($chr,$start,$seq1,$seq2) = split("_",$vid);
				my $end = $start+1;
				my $type = "S";
				
				if (length($seq1) > length($seq2)){
					$end = $start+length($seq1);
					$type = "D";
				}
				if (length($seq1) < length($seq2)){
					$end = $start+length($seq2);
					$type = "I";
				}
				if (length($seq1)>1000 or length($seq2)>1000){
						$type = $type."L";
				}
									
				#warn "COUOCU $vid" if abs(length($seq1) - length($seq2)) > 1500; 	
				my @hohe = split(" ",$vh->{$vid});
				my @he = split(",",$hohe[0]); 
				my @ho = split(",",$hohe[1]); 
				my $pn = $project_name;
				$pn =~ s/NGS20//; 
				my $list;
				if (@ho){
					#$list = join(",",@ho);
					my @uniq;
					my $huniq;
					map{$huniq->{$_}++} @ho;
					map{push(@uniq,$_) unless exists $huniq->{$_}}  @he;
					$list = join(",",(@ho,@uniq));
				}
				else {
					$list = join(",",(@he));
				}
				my $string2 =  $pn . ":" .scalar(@ho).":".scalar(@he).":".$list;
#				#decompose 
#				my ($p1,$nb,$v) = split(":",$string2);
#				my @he1 = split(",",$v);
#				my @ho1 = $he1[0..$nb];
#				my $ts = "NGS20".$p1.";" .join(",",@he);
#				if (scalar(@ho)){
#					$ts .=" ".join(",",@ho)." HO";
#				}
#				chomp($line);
#				die("\n;$line;\n;$ts;\n\n") if $line ne $ts;		
#				my $nbhe = scalar(@he);
#				my $nbho = scalar(@ho);
					
				print TMP "$vid|$start|$end|$nbho|$string2|$type\n";
			}#end vid  
			
		}#end project_name
		warn "finish";close(TMP);
		warn "end par $nbj ".abs($tpart - time);
		#$dbh->disconnect();
		$pm->finish( 0, {data=>$h,jobid=>$jid,file=>$text_file} );

	 	
	 }#end while temp
	#}
$pm->wait_all_children();

warn "end step 1 ".abs(time - $tloop);
my $tindex = time;
 $dbh = DBI->connect( "dbi:SQLite:dbname=$dbfile", "", "",
		{ sqlite_use_immediate_transaction => 0, } );
		$dbh->do("PRAGMA journal_mode = MEMORY") or die $DBI::errstr;;
		$dbh->do("PRAGMA synchronous = OFF") or die $DBI::errstr;;
		
$dbh->do("create index idxx1 on tabLines(idx)") or die $DBI::errstr;
$dbh->disconnect();
warn "*** end creation index $chr :".abs(time -$tindex);
return $dbfile;
}


sub create_dejavu_lite {
my ($dbfile) = @_;
my $dbh = DBI->connect( "dbi:SQLite:dbname=$dbfile", "", "",
		{ sqlite_use_immediate_transaction => 0, } );
		$dbh->do("PRAGMA journal_mode = MEMORY") or die $DBI::errstr;;
		$dbh->do("PRAGMA synchronous = OFF") or die $DBI::errstr;;
		
my $sth = $dbh->prepare(
		'select idx,group_concat(line,"!") from tabLines group by idx  ;')
	 or die $DBI::errstr;

	$sth->execute();
my 	 $tinsertglobal =time;
my $nodejavu = GenBoNoSqlDejaVu->new(dir=>$dir_tmp,mode=>"c");
$nodejavu->create_table($chr);
$nodejavu->dbh($chr)->do(qq{attach database "$dbfile" as db;});
#,variation_type,count(idx) from tabLines group by idx limit 10;

$nodejavu->dbh($chr)->do(qq{insert into __DATA__(_key,_value,start,end,variation_type,ho,projects) select idx,group_concat(line,"!"),start,end,variation_type,sum(ho),count(idx) from db.tabLines group by idx ;});
$nodejavu->dbh($chr)->do(qq{CREATE UNIQUE INDEX if not exists _key_idx  on __DATA__ (_key);});
$nodejavu->dbh($chr)->do(qq{CREATE  INDEX if not exists _start_idx  on __DATA__ (start);});
$nodejavu->dbh($chr)->do(qq{CREATE  INDEX if not exists _end_idx  on __DATA__ (end);});
$nodejavu->dbh($chr)->do(qq{CREATE  INDEX if not exists _type_idx  on __DATA__ (variation_type);});
my $filedb = $nodejavu->db_name($chr);
	#$nodejavu->close();
	unlink $dbfile;
	#	$lmdb->close();
	$nodejavu->dbh($chr)->do(qq{ PRAGMA auto_vacuum = FULL;});
	warn "end $chr insert :".abs(time-$tinsertglobal);
my $tcompress = time ;
my $sth2 = $nodejavu->dbh($chr)->prepare(
		'select rowid,_key,_value,start,end from __DATA__;');
		$sth2->execute();
my $sth3 = $nodejavu->dbh($chr)->prepare(
		'update __DATA__   set _value=? where _key=? ;')
	 or die $DBI::errstr;
	 my $sth4 = $nodejavu->dbh($chr)->prepare(' insert  into __DATA__POSITION (id,start,end) values(?,?,?)') or die $DBI::errstr;
	
	 
while ( my @row = $sth2->fetchrow_array ) {
	#warn Dumper(@row);
	my $rowid = $row[0]; 
	my $id = $row[1];
	my $text = $row[2];
	my $start = $row[3];
	my $end = $row[4];
	my $tp;
	
	#foreach my $l (split("!",$text)) {
	#		my($p,$info) = split(";",$l);
	#		$tp->{$p} = $info;
	#	}
	#my $z = compress($text);
	$sth4->execute($rowid,$start,$end);
	$sth3->execute($nodejavu->encode($text),$id);
}
warn "end compress ".($tcompress -time);
	$nodejavu->dbh($chr)->do(qq{ VACUUM;});
$nodejavu->close();
system("mv $filedb $dir_prod");


#$nodejavu->dbh($chr)->do("DROP  INDEX _key_idx; ")  or die $DBI::errstr;;
#my $lmdb = GenBoNoSqlLmdb->new(dir=>$dir_tmp,mode=>"c",name=>$chr,is_compress=>1,is_index=>1) ;
###  INSERT
#$nodejavu->close();
#warn "end select $chr ";#.abs(time-$t);
#	my $tp;
#	my $i =0;
# my 	 $tinsertglobal =time;
#	 my $t1= time;
#	while ( my @row = $sth->fetchrow_array ) {
#		$i++;
#		my $tem;
#		my $id= $row[0];
#		my $nho = 0;
#		my $nall =0;
#		foreach my $l (split("!",$row[1])) {
#			my($p,$info) = split(";",$l);
#			$tp->{$id}->{data}->{$p} = $info;
#			$nho ++ if $info =~/HO/;
#			$nall ++;
#		}
#		my($chr,$start,$seq1,$seq2) = split("_",$id);
#		my $end = $start+1;
#		
#		if (length($seq1) > length($seq2)){
#			$end = $start+length($seq1);
#		}
#		$tp->{$id}->{start} = $start;
#		$tp->{$id}->{end} = $end;
#		$tp->{$id}->{ho} =  $nho;
#		$tp->{$id}->{all} =  $nall;
#		#$no->put( $id, $nall);
#		#$no->put( $id . "_ho",$nho  ) if $nho > 0;
##		$lmdb->put($id,$tp->{$id});
#		if ($i%100000 == 0){
#			
#			print "==> $i  $chr ".(abs(time-$t1)."\n");
#			$t1 = time;
#		
#			#$noglobal->put_bulk($chr,$tt);
#			$nodejavu->put_dejavu_bulk($chr,$tp);
#			print "\t==> $i  $chr ".(abs(time-$t1)."\n");
#			
#			$t1 = time;
#			$tp = {};
#			#$tt ={};
#		}
#		
#	}
#	
#	$nodejavu->put_dejavu_bulk($chr,$tp); 
#

}
	