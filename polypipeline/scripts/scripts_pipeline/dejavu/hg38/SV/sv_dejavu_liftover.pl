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
#use GenBoNoSqlRocksGenome;
use File::Basename;
use File::Slurp qw(write_file);

#######################################################################
#  1) liste  l'ensemble des fichiers  project.dejavu
#  2) retrieve chaque table de hash correspondante et l'integre dans une table de hash globale 
#  3) pour tout les project_dejavu plus récents que le dejavuglobal : retrieve la hash et l'integre dans la hash globale
#  4) freeze de la table de hash gllobale
########################################################################

#my $cgi = new CGI;

# lister les fichiers allSV
#my $cmd = "ls /data-xfs/Manue/Test_SV/DejaVu/TransLoc/*.store";
#my @res = `$cmd`;

my $chrs =[1..22,'X','Y','MT'];
my $ucsc_chr;
my $hash_num_ucsc;
my $hash_ucsc_num;
my $nc = 0;
foreach my $chr (@$chrs){
	$nc ++;
	 my $u = "chr".$chr;
 	$u = "chrM" if $chr eq "chrMT";
 	push(@$ucsc_chr,$u);
	$hash_num_ucsc->{$nc} = $u;
	$hash_ucsc_num->{$u} = $nc;
	
}
my $halldejavu;
my $nbPatient;
my $buffer = GBuffer->new();
my @releases = ("HG19","HG38");	

# pour acceder aux dejavu de chaque projet
my $dir = $buffer()->buffer()->deja_vu_public_dir($releases[0],"SVeq")#; 
warn $dir;
my $cmd = "ls $dir/projects/*.SVeqDejavu";
my @res = `$cmd`;
chomp(@res);
 my $dir_HG38 = $buffer()->buffer()->deja_vu_public_dir($releases[1],"SVeq");
 
my @HG38_files = `ls $dir_HG38/projects/*.dejavu`;
chomp(@HG38_files);
my %HG38_hash_files;
map{$HG38_hash_files{basename($_)} ++} @HG38_files;
my $nb;
my $nbo = 0;
foreach my $TransLocFile (@res)
{
	chomp($TransLocFile);
	my $filename = basename($TransLocFile);
	my ($projectname,$rien) = split(/\./,$filename);
	next if exists $HG38_hash_files{$filename};
	warn $projectname;
	
	my $hdejavu = retrieve($TransLocFile) or die "Can't retrieve datas from ".$TransLocFile." !\n";

	foreach my $id (keys %{$hdejavu})
	{
			$nbo ++ unless exists $halldejavu->{$id};
			$halldejavu->{$id}->{$projectname} = $hdejavu->{$id};
	}
	
	$nb ++;
}
warn $nbo ;
#my $file_alldejavu = "/data-xfs/Manue/Test_SV/DejaVu/TransLoc/TranslocDejavu.all";
#store(\ %{$halldejavu}, $file_alldejavu) or die "Can't store $file_alldejavu!\n";

#Pour stocker le dejavu global
#save_dejavu($dir,$halldejavu);

my $file_alldejavu = $dir."/SVeqDejavu.all";
my @pos;
foreach my $id (keys %$halldejavu) {
	my ($c1,$p1,$c2,$p2) = split("_",$id);
	
	push(@pos,join("\t",$hash_num_ucsc->{$c1},$p1,$p1+1,$id)."\n");
	push(@pos,join("\t",$hash_num_ucsc->{$c2},$p2,$p2+1,$id)."\n");
}


my $tmp_dir = ".";
my $name = "tmp.".time.".".rand(time);
my $f1 = "$tmp_dir/"."$name.bed";
my $f2 = "$tmp_dir/"."$name.out.bed";
 my $f3 = "$tmp_dir/"."$name.error.bed";
 
  
 write_file("$tmp_dir/"."$name.bed", @pos);
 warn `wc -l $tmp_dir/$name.bed`;
 system("/software/distrib/ucsc_util/liftOver $f1 $RealBin/../hg19ToHg38.over.chain $f2 $f3 ");
   open (BED,"$f2");
  
   my @unsave;
   my $newTotal;
   my $boundary;

   while (my $line = <BED>){
   	chomp($line);
   	my ($chr,$start,$end,$id) = split(" ",$line);
   	$chr = $hash_ucsc_num->{$chr};
   	push(@{$boundary->{$id}->{chr}},$chr);
   	push(@{$boundary->{$id}->{start}},$start);
   
   }
   
 unlink $f1;
 unlink $f2;
 unlink $f3;
 

 
   
foreach my $id (keys %$boundary) {
   	next if scalar @{$boundary->{$id}->{start}} < 2;
   	my $new_id = $boundary->{$id}->{chr}->[0]."_".$boundary->{$id}->{start}->[0]."_".$boundary->{$id}->{chr}->[1]."_".$boundary->{$id}->{start}->[1];
   		$newTotal->{$new_id} = $halldejavu->{$id};
}
 $dir = $buffer()->buffer()->deja_vu_public_dir($releases[1],"SVeq");
 $cmd = "ls $dir/projects/*.SVeqDejavu";
 @res = `$cmd`;

foreach my $TransLocFile (@res)
{
	chomp($TransLocFile);
	my $filename = basename($TransLocFile);
	my ($projectname,$rien) = split(/\./,$filename);
	
	
	my $hdejavu = retrieve($TransLocFile) or die "Can't retrieve datas from ".$TransLocFile." !\n";

	foreach my $id (keys %{$hdejavu})
	{
			$newTotal->{$id}->{$projectname} = $hdejavu->{$id};
	}
	
}
   
 my $dir_HG38 = $buffer->config_path("root","dejavu_SV")."/HG38/SVeq/";
 warn scalar(keys %$newTotal);
 save_dejavu($dir_HG38,$newTotal);

exit(0);


sub save_dejavu {
	my ($dir,$data) = @_;
	my $file_alldejavu = $dir."/SVeqDejavu.all";
	
	store(\ %{$data}, $file_alldejavu) or die "Can't store $file_alldejavu!\n";
	
	my $nodejavu = GenBoNoSqlDejaVuSV->new( dir => $dir, mode => "c" );


	 

	$nodejavu->create_table();

				foreach my $id  (keys %$data)
				{
					$nodejavu->insert_sv($id,$data->{$id},"manta");
						
				}
	$nodejavu->create_index();	

$nodejavu->close();

	warn $dir;
}



 	