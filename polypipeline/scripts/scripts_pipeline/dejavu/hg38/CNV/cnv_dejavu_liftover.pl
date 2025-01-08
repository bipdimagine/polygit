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
use File::Basename;
use File::Slurp qw(write_file);
use GenBoNoSqlDejaVuCNV;


 

my $buffer = GBuffer->new();
my @releases = ("HG38");

my $dir_HG19 = $buffer()->buffer()->deja_vu_public_dir("HG19","CNV")#;
my $dir_HG38 = $buffer()->buffer()->deja_vu_public_dir("HG38","CNV")#;
my @HG19_files = `ls $dir_HG19/projects/*.dejavu`;

chomp(@HG19_files);
my %HG19_hash_files = @HG19_files;
my @HG38_files = `ls $dir_HG38/projects/*.dejavu`;
chomp(@HG38_files);

my %HG38_hash_files;
map{$HG38_hash_files{basename($_)} ++} @HG38_files;

my $id_HG19;
my $n;
warn "start";

my $total;
foreach my $CNVFile (keys %HG19_hash_files){
	
	my $filename = basename($CNVFile);
	my ($projectname,$rien) = split(/\./,$filename);
	#warn $projectname;
	next if exists $HG38_hash_files{$filename};
	
	my $hdejavu = retrieve($CNVFile) or die "Can't retrieve datas from ".$CNVFile." !\n";
	
	foreach my $type (keys %{$hdejavu})
	{
		foreach my $num (keys %{$hdejavu->{$type}})
		{
			
				foreach my $id (keys %{$hdejavu->{$type}->{$num}})
				{		
				
					$total->{$id}->{$projectname} = $hdejavu->{$type}->{$num}->{$id};
					my ( $t, $c, $d, $f ) = split( /_/,$id);
					 $id_HG19->{$c}->{$id}->{start} = $d;
					 $id_HG19->{$c}->{$id}->{end} = $f;
					 $id_HG19->{$c}->{$id}->{type} = $t;
					 $id_HG19->{$c}->{$id}->{data}->{$projectname} = $hdejavu->{$type}->{$num}->{$id};
					#$halldejavu->{$type}->{$num}->{$id}->{$projectname} = $hdejavu->{$type}->{$num}->{$id};
					#$nbCNV++;
				}
		}
	}	
	$n ++;
	warn $n;
	#last if $n >50;
	
}

#save($dir_HG19,$id_HG19);
#die();

my @pos;
foreach my $c (keys %$id_HG19){
	my $hash = $id_HG19->{$c};
	foreach my $id (sort {$hash->{$a}->{start} <=> $hash->{$b}->{start}} keys %$hash){
		my $chr = $c;
		$chr = "chr".$c;
		$chr = "chrM" if $chr eq "chrMT";
		push(@pos,$chr."\t".$hash->{$id}->{start}."\t".$hash->{$id}->{end}."\t".$hash->{$id}->{type}."\t".$id."\n");
	}
}

my $tmp_dir = ".";
my $name = "tmp.".time.".".rand(time);
my $f1 = "$tmp_dir/"."$name.bed";
my $f2 = "$tmp_dir/"."$name.out.bed";
 my $f3 = "$tmp_dir/"."$name.error.bed";
 
  
 write_file("$tmp_dir/"."$name.bed", @pos);
 my $chain = $buffer->liftover_chain_file("HG19","HG38");
 system("/software/distrib/ucsc_util/liftOver $f1 $chain $f2 $f3 ");
 
 
   system("/software/distrib/ucsc_util/liftOver $f1 /data-isilon/bipd-src/pnitschk/git2-polyweb/polygit/polymorphism-cgi/cache_nodb/scripts/rocks/hg19ToHg38.over.chain $f2 $f3");
   my $nn ;
   open (BED,"$f2");
  my $ids_HG38;
   my @unsave;
   my $newTotal;
   while (my $line = <BED>){
   		chomp($line);
   		my ($chr,$start,$end,$type,$id) = split(" ",$line);
   		$chr =~ s/chr//;
   		$chr = "MT" if ($chr eq "M");
   		my $new_id = join("_",$type,$chr,$start,$end);
   		$nn ++;
   		$ids_HG38->{$chr}->{$new_id}->{data} = $id_HG19->{$chr}->{$id}->{data};
   }
   warn scalar(@pos)." ".$nn;
   
   	unlink $f1;
   	unlink $f2;
   	unlink $f3;
 
foreach my $filename (keys %HG38_hash_files){
	my $CNVFile = $dir_HG38."/projects/".$filename;
	my $filename = basename($CNVFile);
	my ($projectname,$rien) = split(/\./,$filename);
	#warn $projectname;
	
	my $hdejavu = retrieve($CNVFile) or die "Can't retrieve datas from ".$CNVFile." !\n";
	foreach my $type (keys %{$hdejavu})
	{
		foreach my $num (keys %{$hdejavu->{$type}})
		{
			
				foreach my $id (keys %{$hdejavu->{$type}->{$num}})
				{		
					$newTotal->{$id}->{$projectname} = $hdejavu->{$type}->{$num}->{$id};
					my ( $type, $chr, $d, $f ) = split( /_/,$id);
					$ids_HG38->{$chr}->{$id}->{data}->{$projectname} = $hdejavu->{$type}->{$num}->{$id} ++;
				}
		}
	}	
	$n ++;
}
warn " ==> 1 ".scalar(keys %{$id_HG19->{1}})." ".scalar(keys %{$ids_HG38->{1}});
save($dir_HG38,$ids_HG38);

 
exit(0);


sub save {
	my ($dir,$ids) = @_;
	my $nodejavu = GenBoNoSqlDejaVuCNV->new( dir => $dir, mode => "c" );

	 warn $dir;
	 
	foreach my $chr (keys %{$ids})
	{
	next unless $chr;
	next if $chr =~ /KI/;
	next if $chr =~ /GL/;
	next if length($chr) > 4;
	warn $chr;
	$nodejavu->create_table($chr);
	my $xx =0;
	my $nbid =0;
	my @string;
	
	foreach my $id  (keys %{$ids->{$chr}})
				{	
					$nbid ++;
					my @prjs = keys %{$ids->{$chr}->{$id}->{data}};
					my $hs; 
					$hs->{caller}->{manta} =0;
					$hs->{caller}->{canvas} =0;
					$hs->{caller}->{wisecondor} =0;
					$hs->{projects} = \@prjs;
					$hs->{project_sample} =[];
					$hs->{projects_stats};
					my $ps ={};
					my $st ="";
					
					foreach my $pr (@prjs){
						foreach my $sample (keys %{$ids->{$chr}->{$id}->{data}->{$pr}}){
							$ps->{$pr}->{patients} ++;
							push(@{$hs->{patients}},$sample);
							$st.="$sample";
							my $res ="";
							$ps->{$pr}->{sr} = 0;
							$ps->{$pr}->{depth} = 0;
							$ps->{$pr}->{coverage} = 0;
							foreach my $caller (keys %{$ids->{$chr}->{$id}->{data}->{$pr}->{$sample}}){
								$ps->{$pr}->{sr} ++ if $caller eq "manta";
								$ps->{$pr}->{depth} ++ if $caller eq "canvas";
								$ps->{$pr}->{coverage} ++ if $caller eq "wisecondor";
								$hs->{caller}->{$caller} ++;
								$res .= "!X!%_".$caller." ";
							}
							push(@{$hs->{manue_string}->{$pr}},$pr.":".$sample.":".$res);
						}
						$ps->{$pr}->{string} = join("",@{$hs->{manue_string}->{$pr}});
						
					}
					
					my $dv_project = scalar(@prjs);
					next if $dv_project == 0;
					my $dv_sample = scalar(@{$hs->{patients}});
					my $dv_sr = 0+$hs->{caller}->{manta};
					my $dv_depth = 0+$hs->{caller}->{canvas};
					my $dv_cov = 0+$hs->{caller}->{wisecondor};
					$nodejavu->insert_cnv($id,$ids->{$chr}->{$id}->{data},$ps,$dv_project,$dv_sample,$dv_sr,$dv_depth,$dv_cov);
				}
	$nodejavu->create_index($chr);
	}
$nodejavu->close();
}

#sub save2 {
#	my ($dir,$ids) = @_;
#	my $nodejavu = GenBoNoSqlDejaVu->new( dir => $dir, mode => "c" );
#
#	 
#foreach my $chr (keys %{$ids})
#{
#	next if $chr =~ /KI/;
#	next if $chr =~ /GL/;
#	next if length($chr) > 4;
#	$nodejavu->create_table($chr);
#	my $sth = $nodejavu->dbh($chr)->prepare(
#		'insert into  __DATA__(_key,_value,start,end,variation_type,ho,projects)  values(?,?,?,?,?,?,?) ;') or die $DBI::errstr;
#		$sth->execute();
#				my $tree;
#				
#				foreach my $id  (keys %{$ids->{$chr}})
#				{
#								my ( $t, $c, $d, $f ) = split( /_/,$id);
#								my $text;
#								$sth->execute($id,$nodejavu->encode($ids->{$c}->{$id}->{data}),$d,$f,$t,0,0);
#								my $id2 = $nodejavu->dbh($chr)->sqlite_last_insert_rowid();
#								#warn Dumper $total->{$id};
#	  							$nodejavu->sth_insert_position_cached($chr)->execute($id2,$d,$f) ;
#				}
#$nodejavu->dbh($chr)->do(qq{CREATE UNIQUE INDEX if not exists _key_idx  on __DATA__ (_key);});
#$nodejavu->dbh($chr)->do(qq{CREATE  INDEX if not exists _start_idx  on __DATA__ (start);});
#$nodejavu->dbh($chr)->do(qq{CREATE  INDEX if not exists _end_idx  on __DATA__ (end);});
#$nodejavu->dbh($chr)->do(qq{CREATE  INDEX if not exists _type_idx  on __DATA__ (variation_type);});
#
#}
#
#$nodejavu->close();
#}
