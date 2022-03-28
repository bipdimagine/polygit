#!/usr/bin/perl

use FindBin qw($RealBin);
use FindBin qw($Bin);
use FindBin qw($RealBin);
use lib "$RealBin/../GenBo";
use lib "$RealBin/../GenBo/lib/GenBoDB";
use lib "$RealBin/../GenBo/lib/obj-nodb";
use lib "$RealBin/../GenBo/lib/obj-nodb/packages";

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
 
  use LMDB_File qw(:flags :cursor_op);
  
my $fork;
my $force ;
my $verbose=1;
my $only_chr_name;
GetOptions(
	'fork=i'  	=> \$fork,
	'verbose=i' =>\$verbose,
	'force=i' =>\$force,
	'chr=s' =>\$only_chr_name,
);
$|=1; 
my $CHANGE_SOMETHING;
my %exclude;
my @chromosomes = (1..22,'X','Y','MT');
@chromosomes = ($only_chr_name) if $only_chr_name;
#my @chromosomes = (21..22);#,'X','Y','MT');
my @sub_files  = ("merge_count_ho","merge_count","merge");
my $buffer = GBuffer->new();




my $db_final;
$SIG{'INT'} = sub {print "Caught One!\n";
	$db_final->close if $db_final;
	exit(0);
	die();
};

#warn "### Step 1: exclude this project";
#my $dir_root = $buffer->config->{public_data}->{HG19}."/".$buffer->config->{kyoto}->{deja_vu}."/lite";
#my $exclude_project_file = $dir_in."/exclude_project.txt";
#warn $exclude_project_file;
#if (-e $exclude_project_file){
#	my @lproject = `cat $exclude_project_file`;
#	chomp(@lproject);
#	warn "exclude project : \n".join("\n",@lproject);
#	chomp(@lproject);
#	map{$exclude{$_} = 1} @lproject;
#} 
$exclude{NGS2014_0001} = 1;
sleep(1);
my $dir_root = $buffer->config->{public_data}->{HG19}."/".$buffer->config->{kyoto}->{deja_vu}."/lite";
unless (-d $dir_root){
	system(" mkdir  $dir_root");
	system("chmod a+rwx $dir_root");
}


my $dir_projects = $dir_root."/projects/";
unless (-d $dir_projects){
	system(" mkdir -p $dir_projects");
	system("chmod a+rwx $dir_projects");
}



#
my $no_reload = undef;
my $dir_tmp = tempdir ( DIR => $dir_root,CLEANUP => 1 );
my @projects = grep { !(exists $exclude{$_})}  @{$buffer->listProjectsForDejaVu()}; 
my $pm = new Parallel::ForkManager($fork);
#$no_reload =1;
unless ( defined $no_reload){
my %list_files;
 opendir (DIR, $dir_projects) or die $!;
  while (my $file = readdir(DIR)) {
		my ($pname,$ext) = split(/\./,$file);
		$list_files{$pname} = $dir_projects.$file;
       # print "$file\n";
    }    
  close (DIR);


my @restart_projects;
  my $notodo = GenBoNoSql->new(dir=>$dir_projects,mode=>"r");
foreach my $pname (@projects){
		next if $pname =~/2010/;
	#	warn $pname;
		my $b =   GBuffer->new();
	my $p = $buffer->newProject(-name=>"$pname");
		my $d;
		eval{
		 	$d = $p->getRootDir()."/variations/";
		};
		if ($@){
			next;
		}
			my @files = grep {$_ !~/tbi/} `find $d  -name "*.vcf*" `;
			chomp(@files);
			$d = $p->getRootDir()."/indels/";
			my @files2 =  grep {$_ !~/tbi/} `find $d  -name "*.vcf*" `;
			chomp(@files2);
			push(@files,@files2);
			#warn Dumper(@files);
			next unless @files;
			unless (exists $list_files{$pname}){
		push(@restart_projects,$pname);
		next;
	}
			unless ( $notodo->exists_db($pname)){
				warn "$pname nodatabase";
				push(@restart_projects,$pname);
				next;
			}
			push(@files,$list_files{$pname});
			my @sorted = sort{ -M $a <=> -M $b} @files;
			if ($sorted[0] ne $list_files{$pname}){
					push(@restart_projects,$pname);
					next;
			}
			
				unless ( $notodo->exists_db($pname)){
					push(@restart_projects,$pname);
					next;
			}
			
			my $patients;
			
			 map{$patients->{$_->name} = 1} @{$p->getPatients()};
			 my $patients_lite = $notodo->get($pname,"patients");
			 map{$patients->{$_} ++}  keys %{$patients_lite} ;
			 my @diff = grep{$patients->{$_} ne 2} keys %$patients;
			 
			
			
	
			if (@diff){
				#warn Dumper @diff;
			
				push(@restart_projects,$pname);
				next;
			}

#		warn "ok";
		
	
	
}
$notodo->close();
$notodo = undef;
warn Dumper @restart_projects;
#die();

#exit(0) unless @restart_projects;
my $npp = scalar(@restart_projects);
my $cpp = 0;
foreach my $pname(@restart_projects){
	$cpp++;
	#warn "cache $pname $npp/$cpp";
	eval {
		warn "$RealBin/cache.pl -project=$pname -type=dejavu -fork=$fork";
		system("$RealBin/cache.pl -project=$pname -type=dejavu -fork=$fork");
	};
	
	#cache_dejavu2($pname,$fork,$dir_projects);
}

#exit(0);
warn "### Step 1bis: prepare project  ";#my $projects = $buffer->listProjects();

#warn Dumper(@projects);
my $nbp=0;

my $pr = String::ProgressBar->new( max => scalar(@projects) );


foreach my $project_name ( @projects){
 	next unless -e $dir_projects."/$project_name.lite";
   my $notodo = GenBoNoSql->new(dir=>$dir_projects,mode=>"r");
   my $noprojects = GenBoNoSqlDejaVu->new(dir=>$dir_tmp,mode=>"w");
   
   my $patient = $notodo->get($project_name,"patients");
 next unless $patient;
   my %h2 = reverse %$patient;
   $noprojects->put("projects",$project_name,\%h2);
   
   
 }
}

# exit(0);
my $c =0;
$| = 1;
  my $tt =0;
  
  unless (-e $dir_tmp){
  	mkdir $dir_tmp;
  	system("chmod a+rwx $dir_tmp");
  }
  #system("rm $dir_tmp/*");
  warn 
  warn "---------------------------------------\n";
    warn "--------------Compute lite database ---------------------\n";
      warn "---------------------------------------\n";
      
   my $hash;
   my $nbp = scalar(@projects);

  #last;
   warn "---------------------------------------\n";
   
foreach my $chr (@chromosomes){
	  my $pid = $pm->start and next;
	   my $tp =0;
	   my $temp_hash;
  		my $nodejavu = GenBoNoSqlDejaVu->new(dir=>$dir_tmp,mode=>"w");
    	  
 foreach my $project_name (sort @projects){
 	next unless -e $dir_projects."/$project_name.lite";
 	 my $count;
 	  my $count_ho={};
 	$tp ++;
 	$tt++;
# 	warn "chr: $chr ". $tt."/$nbp" if $tt%100==0;
 	
 	my $list;
   my $notodo = GenBoNoSql->new(dir=>$dir_projects,mode=>"r");
  
  
    my $vh = $notodo->get("$project_name",$chr,1);
    my $last;

    	foreach my $vid (keys %$vh){
    		$last = $vid;
    		$temp_hash->{$vid}->{$project_name} = $vh->{$vid};
    		
    		
     	}
     	$count = {};
     	$notodo->close();
		$notodo = undef;
}
 #$notext->dbh($chr."_count")->do("INSERT INTO __DATA__(__DATA__) VALUES('optimize');");
my $t = time;
warn "start : insert  $chr";
my $i=0;
my $tt;
my $tp;
my $max = scalar( keys %$temp_hash);
#$noglobal->create_table($chr);
#$noglobal->dbh($chr)->do("DROP  INDEX  if  exists  __DATA__._key_idx ")  or die $DBI::errstr;;
$nodejavu->create_table($chr);
$nodejavu->dbh($chr)->do("DROP  INDEX  if  exists  __DATA__._key_idx ")  or die $DBI::errstr;;

foreach my $id (keys %$temp_hash){
			$i++;
		my $text = encode_json $temp_hash->{$id};
		my($chr,$start,$seq1,$seq2) = split("_",$id);
		my $end = $start+1;
#		$tp->{$id}->{start}= $start;
#		$tp->{$id}->{end}= $start;
		if (length($seq1) > length($seq2)){
			$end = $start+length($seq1);
#		
#		for (my $i =1; $i<length($seq1);$i++){
# 					$tp->{$start+$i}->{$id} ++;
# 			}
#
		}
		$tp->{$id}->{start} = $start;
		$tp->{$id}->{end} = $end;
		
		$tp->{$id}->{data} =  $temp_hash->{$id};
		my @ho = grep {$_=~/HO/} values %{$temp_hash->{$id}};
		my $nho = scalar(@ho);
		my @all = values %{$temp_hash->{$id}};	
		my $nall = scalar(@all);
		$tp->{$id}->{ho} =  $nho;
		$tp->{$id}->{all} =  $nall;
		#$tt->{$id} = $temp_hash->{$id};
	
	#	$notext->put($chr,$id." ".$text);
		if ($i%100000 == 0){
				print "==> $i/$max  $chr \n";
			#$noglobal->put_bulk($chr,$tt);
			$nodejavu->put_dejavu_bulk($chr,$tp);
			$tp = {};
			#$tt ={};
		}
     	#warn encode_json $temp_hash->{$last};
   }
  $nodejavu->put_dejavu_bulk($chr,$tp); 
#$noglobal->put_bulk($chr,$tt);#
#$noglobal->dbh($chr)->do("CREATE UNIQUE INDEX if not exists _key_idx  on __DATA__ (_key); ")  or die $DBI::errstr;;
$nodejavu->dbh($chr)->do("CREATE UNIQUE INDEX if not exists _key_idx  on __DATA__ (_key); ")  or die $DBI::errstr;;
#warn "end global ---> $chr :: ".abs(time -$t);
#warn "end $chr ".abs(time- $t);
$pm->finish();

  }
  warn "wait -----";
$pm->wait_all_children();

warn "*-*-*-* start backup sqlite *-*-*-*";
my $dir_backup = $dir_root."/backup";

unless (-e $dir_backup){
mkdir $dir_backup;
system("chmod a+rwx $dir_backup");
}



 my $nodejavu = GenBoNoSqlDejaVu->new(dir=>$dir_root,mode=>"w");
  my $lmdb_prod_name = $nodejavu->lmdb_extension;
  

warn "wait ----- extra homozygote/heterozygote  ";

  my $nodejavu = GenBoNoSqlDejaVu->new(dir=>$dir_root,mode=>"r");
 my $ext =  $nodejavu->extension();
 warn $dir_tmp;
 die("problem with dir_temp : $dir_tmp ") unless -e $dir_tmp;
 sleep(50);
my $root_files = get_files($dir_tmp);
warn "*/*/*/*/ start backup $dir_root /*/*/*/* ";
# my $root_files = get_files($dir_root);
  my @backup_files;
  my $pm2 = new Parallel::ForkManager($fork);
  push(@backup_files,"projects.$ext");
   foreach my $chr (@chromosomes){
   	push(@backup_files,"$chr.$ext");
   }
   
 foreach my $file_name (@backup_files){
 	my $f = $dir_root."/".$file_name; 
 	my $f1 = $dir_tmp."/".$file_name; 
 	my $pid = $pm2->start and next;
 	warn "\t backup $file_name";
 		logrotate_File($f,$dir_backup) if  -e $f && -s $f;
 		rename $f1 , $f if -e $f1 && -s $f1;
 		warn "\t Done $file_name";
	$pm2->finish();	
	}	
$pm2->wait_all_children;

warn "*/*/*/*/ END --> backup & move /*/*/*/* \n";


 ###
 ####
 ###
 
 warn "*/*/*/*/ create He/Ho database /*/*/*/* \n";
 
 
  my $nodejavu = GenBoNoSqlDejaVu->new(dir=>$dir_root,mode=>"w");
  my $lmdb_prod_name = $nodejavu->lmdb_extension;
  my $lmdb_temp_name = $lmdb_prod_name."_temp";
	my $dir_lmdb_temp = $dir_root."/".$lmdb_temp_name;
	my $dir_lmdb_prod = $dir_root."/".$lmdb_prod_name;
	
  $nodejavu->lmdb_extension($lmdb_temp_name);
if (-e $dir_lmdb_temp){
	system ("rm -r $dir_lmdb_temp");
}

  my %hash;

my $dir_out2 = $buffer->config->{public_data}->{HG19}."/".$buffer->config->{kyoto}->{deja_vu}."lite/lmdb_chr";
unless (-e $dir_out2){
	mkdir $dir_out2;
	system("chmod a+rwx $dir_out2");
}

   
 foreach my $chr (@chromosomes){
 	 my $sth = $nodejavu->dbh($chr)->prepare("SELECT _key,ho,projects FROM __DATA__;") or die $DBI::errstr;
	$sth->execute();
	my $no = GenBoNoSqlLmdb->new(dir=>$dir_out2,mode=>"c",name=>$chr);
  my $h;
  my $nb =0;
  while (my @row = $sth->fetchrow_array ) {
  	$nb ++;
  	$row[1] +=0;
  	$row[2] +=0;
  	if ($nb%100000 == 0){
  			warn $chr." ".$nb;
  	}
     $nodejavu->put_lmdb($row[0],$row[2]); 
      $nodejavu->put_lmdb($row[0]."_ho",$row[1]) if $row[1] > 0; 
      $no->put($row[0],$row[2]); 
    	$no->put($row[0]."_ho",$row[1]) if $row[1] > 0; 
  }
  $no->close();
 }
 warn "----\n";
 my $nb_total_temp = 	$nodejavu->lmdb()->stat->{entries};
 warn $nb_total_temp;
  warn "----\n";
 $nodejavu->close();
 warn "tar -cvf $dir_lmdb_prod.tar $dir_lmdb_prod \n";
 system("tar -cvf $dir_lmdb_prod.tar $dir_lmdb_prod 2>/dev/null");
 logrotate_File($dir_lmdb_prod.".tar",$dir_backup) if -e $dir_lmdb_prod.".tar"; 
 system("rm $dir_lmdb_prod.tar");

 system("rm -r $dir_lmdb_prod");
 system ("mv $dir_lmdb_temp $dir_lmdb_prod");
 system "chmod a+rwx $dir_lmdb_prod";
  system "chmod -R a+rw $dir_lmdb_prod";
  $nodejavu = GenBoNoSqlDejaVu->new(dir=>$dir_root,mode=>"r");
 my $nb_total_prod = 	$nodejavu->lmdb()->stat->{entries};
 warn $nb_total_prod;
 warn "----\n";

 warn "*/*/*/*/ That's all folks !!!! /*/*/*/* \n";

exit(0);



sub logrotate_File {
	my ($file,$dir_backup) = @_;
	
	if (-e $file) {
		my $merge_rotate = new Logfile::Rotate (
				File  => $file,
				Count => 3,
				Gzip  => 'lib',
				Dir	  => $dir_backup,
				Flock => 'no' );		
		$merge_rotate->rotate();
	}
	
}

sub get_files {
	my ($dir) =@_;
	warn $dir;
	my @toto = `ls $dir/*.lite`;
	chomp(@toto);
	return \@toto;
	warn Dumper @toto;
	my $fu = File::Util->new;
	my @contents = grep {!(-d $_) }$fu->list_dir( $dir=>{no_fsdots =>1});
	return \@contents;
}




sub cache_dejavu2 {
	my ($project_name,$fork,$root_dir) = @_;
	mkdir $root_dir unless -e $root_dir;
	
	my $buffer1 = new GBuffer;
	my $project = $buffer1->newProject( -name => $project_name, -verbose =>1 );
	my @chr_names = map{$_->name} @{$project->getChromosomes};
	my $pm = new Parallel::ForkManager($fork);
	my @patient_names =  sort{$a cmp $b} map{$_->name} @{$project->getPatients};
	my $hpatients;
	for (my $i=0;$i<@patient_names;$i++){
		$hpatients->{$patient_names[$i]} = $i;
	}
	warn $root_dir;
	my $no = GenBoNoSql->new(dir=>$root_dir,mode=>"w");
	$no->put($project_name,"patients",$hpatients);
	#unlink $root_dir."$project_name.lite";
	$no->close;
	$no = undef;
	
	foreach my $chr_name (@chr_names) {
		 #next if $chr ne "15";
			my $pid = $pm->start and next;
			my $buf = new GBuffer;
			my $pr = 	$buf->newProject( -name => $project_name, -verbose =>1 );
			my $chr = $pr->getChromosome($chr_name);
			my $references = [$chr];#->getReferences(undef,undef,50);
			my $results2;
			my $nr =0;
			foreach my $r (@$references){
			my $vs = $r->getStructuralVariations;
			$nr ++;
			warn "----> $chr_name $nr/50 " ;#if $nr%5 ==0;
			foreach my $v (@$vs){
				my $key = $v->id;
				my $aho= [];
				my $ap=[];
				foreach my $p (@{$v->getPatients()}){
					my $pn = $p->name();
					my $patient_id = $hpatients->{$pn};
					push(@$ap,$patient_id);
					
					#$results2->{data}->{$key}->{patients}->{$pn} ++; 
					push(@$aho,$patient_id) if ($v->isHomozygote($p));
			}
			$results2->{data}->{$key} = join(",",sort{$a <=> $b} @$ap);
			if (scalar(@$aho)) {
				$results2->{data}->{$key} .=" ".join(",",sort{$a <=> $b} @$aho)." HO";
			}
			}#enbd var
			$pr->purge_memory_reference($r->id);
			}#end reference
			
			
				my $no = GenBoNoSql->new(dir=>$root_dir,mode=>"w");
				  	$no->put($project_name,$chr_name,$results2->{data});
      				$no->close;
      				$no = undef;
      				delete $results2->{data};
			$pm->finish();	
	}	
		$pm->wait_all_children;
	
	$project->buffer->dbh->disconnect();
		return 1;
}



sub lmdb {
	my ($path,$name) = @_;
	  my $env = LMDB::Env->new($path, {
      mapsize => 100 * 1024 * 1024 * 1024, # Plenty space, don't worry
      maxdbs => 20, # Some databases
      mode   => 0777,
      #flags => LMDB_File::MDB_RDONLY,
      # More options
  });
 my $txn = $env->BeginTxn(); # Open a new transactionif ()
 my $db;
 $db = $txn->OpenDB( {    # Create a new database
      dbname => $name,
      flags => MDB_CREATE
  });
return $db;	
	
}

sub lmdb_key {
	my ($key) = @_;
	  if (length($key) > 500){
     	my (@s) = split ("_",$key);
     	my $m = md5_hex($s[2]."_".$s[3]);
     	$key = $s[0]."_".$s[1]."_".$m;
     }
     return $key;
}

