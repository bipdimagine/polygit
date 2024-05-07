#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin/../../../../../GenBo/lib/obj-nodb/";
use lib "$Bin";
use strict; 
#use GBuffer;
use Data::Printer;
use Getopt::Long;
use Carp;
use Data::Dumper;
use Getopt::Long; 
use Storable qw/retrieve freeze thaw nfreeze nstore_fd nstore/;
use Set::IntSpan::Fast::XS ;
use Sys::Hostname;
 use JSON::XS;
 use bytes;
use Bio::DB::HTS::VCF ;
use Parallel::ForkManager;
use Date::Tiny;
use Vcf;
use GenBoNoSqlRocksGenome;
use Sereal qw(sereal_encode_with_object sereal_decode_with_object write_sereal_file read_sereal_file);
use Time::ETA;
use Bio::DB::HTS::Tabix;
my $chr_name;
my $version;
my $fork;
my $genome_version;
my $merge;
my $type = "genome";

GetOptions(
	'chr=s' => \$chr_name,
	'version=s' => \$version,
	'genome=s' => \$genome_version,
	'fork=s'  => \$fork,
	'merge=s'  => \$merge,
	'type=s'  => \$type,
);
die("fork") unless $fork;
die("genome") unless $genome_version; 

#die("-type=genome or exome") unless $type_vcf;
die("-version= (/data-xfs/public-data/HG19/vcf/gnomad/(version)/") unless $version;

my $dir_root = "/data-isilon/public-data/repository/$genome_version";
my $config = {
		exome=>{
			dir=>$dir_root."/gnomad/$version/vcf/",
			file=>"gnomad.exomes.v$version.sites.chr$chr_name.vcf.bgz",
			ext=>".vcf.bgz",
			name=>"gnomad-exome",
		},
		genome =>{
			dir=>$dir_root."/gnomad/$version/vcf/",
			file=>"gnomad.genomes.v$version.sites.chr$chr_name.vcf.bgz",
			ext=>".vcf.bgz",
			name=>"gnomad-genome",
		}
	
};
my @vv = split(".",$version);

#my $vcf1 = $config->{genome}->{dir}. $config->{genome}->{file};
my $dir_out;
my $mode;

	$dir_out = "/data-isilon/public-data/repository/$genome_version/gnomad/".$version."/rocksdb_split/";
	$mode ="c";
$dir_out = "/data-isilon/public-data/repository/$genome_version/gnomad/".$version."/rocksdb_split/";
my $vcf1 = $config->{genome}->{dir}. $config->{genome}->{file};
system("mkdir -p $dir_out") unless -e $dir_out;
die($vcf1) unless -e $vcf1;
my $vcf2;
if ($vv[0] < 3 ){
	 $vcf2 = $config->{exome}->{dir}. $config->{exome}->{file};
	die($vcf2) unless -e $vcf2;
}	
# require "$Bin/save.pm";
my @pops =       ("afr","amr","asj","eas","fin","nfe","oth","sas","remaining");
my @order_pops = ("afr","ami","amr","asj","eas","fin","mid","nfe","oth","remaining","sas");
my @types= ("AC","AN","HO"); #HO synonym in vcf nhomalt
#my @sex_types= ("male","female");
my @sex_types= ("XY","XX");
my @filter_type = ("AC0","RF","InbreedingCoeff","LCR","SEGDUP");

my $keys =["rs","ac","an","xy","ho","min","max","min_pop","max_pop"];
	#rs ac an xy ho min max max_pop min pop	
my $pack= "w5 f2 A3 A3 ";

my $rg = GenBoNoSqlRocksGenome->new(dir=>$dir_out,mode=>"c",chromosome=>$chr_name,genome=>$genome_version,pack=>$pack,description=>$keys,pipeline=>1);#,pack=>"w7A3A3A*");
my $nb_regions = scalar(@{$rg->regions});
 	my $eta = Time::ETA->new(
    milestones =>$nb_regions,
);





my $pm = new Parallel::ForkManager($fork);
my $hash_proc;

$pm->run_on_finish (
		sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hres) = @_;
			delete $hash_proc->{$hres->{process}};
			 $eta->pass_milestone();
 			my $remaining = $eta->get_remaining_time();
			my $elapsed = $eta->get_elapsed_time();
			my $p = int($eta->get_completed_percent());
			print $chr_name." complete $p % => elapsed :".$elapsed." remaining : ".$remaining."\n";
		});

my $idp = time;
#21-9999984-T-C
#0050697561!2
foreach my $r (@{$rg->regions}){
	$idp ++;
	$hash_proc->{$idp} ++;
	#next unless $r->{start} < 9999984 && $r->{end} >9999984 ;
	
	my $pid = $pm->start and next;
	my $rz = GenBoNoSqlRocks->new(dir=>"/data-beegfs/tmp/",mode=>"c",name=>"TMP_".$r->{id},temporary=>1);#,pack=>"w7A3A3A*");
	
	my $hmemory ={};
	warn $vcf1." 1 ";
	run_region($r,$vcf1,$rz);
	warn "$vcf2  2";
	run_region($r,$vcf2,$rz);
	 my $no = $rg->nosql_rocks_tmp($r);
	 $no->put_raw("date",time);
	my $iter = $rz->rocks->new_iterator->seek_to_first;

	while (my ($key, $value) = $iter->each) {
		my $h = $no->decode($value);
		
		my $hfinal ={};
		#warn Dumper ($hmemory->{$key});
		my @oc;
		next if $h->{an} < 1;
		next unless $h->{an};
		 @oc = sort{($h->{min_max}->{$a}->[0]/$h->{min_max}->{$a}->[1]) <=>( $h->{min_max}->{$b}->[0]/$h->{min_max}->{$b}->[1])} keys %{$h->{min_max}};
		my $min_pop = $oc[0];
		my $max_pop = $oc[-1];
		my $values;
		#if ($key == 9999984){

			my $max = ($h->{min_max}->{$max_pop}->[0]/$h->{min_max}->{$max_pop}->[1]);
			my $min = ($h->{min_max}->{$min_pop}->[0]/$h->{min_max}->{$min_pop}->[1]);
			#"rs","ac","an","xy","ho","min","max","min_pop","max_pop"
			my @toto = ($h->{rs},$h->{ac},$h->{an},$h->{xy},$h->{ho});
			$values  = pack("w5",@toto);
			
			#$values = pack($pack,($h->{rs},$h->{ac},$h->{an},$h->{xy},$h->{ho},$min,$max,$min_pop,$max_pop));
			warn $h->{rs} if $key eq "0050697561!2";
			my $test = pack("w ",$h->{rs});
			# if $key eq "0050697561!2";
			$test.= pack ("w ",int($h->{ac}));
			$test.= pack("w ",$h->{an});
			$test.= pack("w ",$h->{xy});
			$test.= pack("w ",$h->{ho});
			$test.= pack("f ",$min);
			$test.= pack("f ",$max);
			$test.= pack("A3 ",$min_pop);
			$test.= pack("A3 ",$max_pop);
			warn Dumper unpack ("w5 f2 A3 A3",$test) if $key eq "0050697561!2" ;
			warn "-----" if $key eq "0050697561!2" ;
			
			#warn sprintf("%010d", $key)."!".$k;
				warn Dumper unpack("w5 f2 A3 A3",$test) if $key eq "0050697561!2";
				warn Dumper ($h->{rs},$h->{ac},$h->{an},$h->{xy},$h->{ho},$min,$max,$min_pop,$max_pop) if $key eq "0050697561!2";
				warn $pack if $key eq "0050697561!2";
			$no->put_raw($key,$test);
		
	}
	$no->close();
	warn "end --------";
	$rz->close();
	my $hres;
	$hres->{process} = $idp;
	$pm->finish(0,$hres);
}

$pm->wait_all_children;

die(Dumper $hash_proc) if scalar(keys %$hash_proc);
warn "END";

  
exit(0);

sub run_region {
	my ($region,$file,$no) = @_;
	
	my $hvariant = {};
	
	#my $no = $rg->nosql_rocks($region);
	my $v = Bio::DB::HTS::Tabix->new( filename => "$file" );
	
	
	
	
	my $h = $v->header();
	my $chr ="";
	#my ($test) = grep {$_ eq "chr1"} @{$h->get_seqnames()};
	$chr= "chr" ;#if $test;
	warn $chr.$region->{tabix};
	
	my $iter  = $v->query($chr.$region->{tabix}); 
	#my $iter  = $v->query($chr."22:50697559-50697562"); 
	return unless $iter;
	my $nb =0;
	while (my $result = $iter->next) {
		my $vh = parse_vcf_line($result,$no);
		
		my $hvariant;
		$hvariant = $no->get($vh->{zid});
		$hvariant= {} unless $hvariant;
		foreach my $pop (@order_pops){
			my $info = "AC_".$pop;
			$vh->{infos}->{$info} = 0 unless exists $vh->{infos}->{$info};
			my $value = $vh->{infos}->{$info};
			$info = "AN_".$pop; 
			$vh->{infos}->{$info} = 0 unless exists $vh->{infos}->{$info};
			my $value_AN = $vh->{infos}->{$info};
			$hvariant->{ac} += $value;
			$hvariant->{an} += $value_AN;
			if (defined $value_AN  && $value_AN > 0){
				$hvariant->{min_max}->{$pop}->[0]  += $value;
				$hvariant->{min_max}->{$pop}->[1]  += $value_AN;
			}
			
			$info = "nhomalt_".$pop;
			$vh->{infos}->{$info} = 0 unless exists $vh->{infos}->{$info};
			$hvariant->{ho} += $vh->{infos}->{$info};
			
			 $info = "AC_".$pop."_XY";
			 unless (exists $vh->{infos}->{$info}){
			 	$info = "AC_".$pop."_male";
			 }
			 $vh->{infos}->{$info} = 0 unless exists $vh->{infos}->{$info};
			 $hvariant->{xy} += $vh->{infos}->{$info};
			 
			 
			
		}
		$hvariant->{rs} = $vh->{rs};
	if ($vh->{zid} eq "0050697561!2"){
	 	warn Dumper ($hvariant);
	 	delete $vh->{infos};
		warn Dumper $vh;
	 	
	 }
	$no->put($vh->{zid},$hvariant);
	
#
#		foreach my $pop (@order_pops){
#			my $info = "nhomalt_".$pop;
#			my $value = return_value($h,$result,$info);
#			next unless defined $value;
#			$hvariant->{ho} += $value;
#		
#		}
#	
#
#	foreach my $pop (@order_pops){
#			my $info = "AC_".$pop."_XY";
#			my $value = return_value($h,$result,$info);
#			unless (defined $value){
#				my $info = "AC_".$pop."_male";
#				$value = return_value($h,$result,$info);
#			}
#			next unless defined $value;
#			$hvariant->{xy} += $value;
#	}
#	 if ($zid eq "0050697561!2"){
#	 	warn Dumper ($hvariant);
#	 	warn  $result->reference()." :: ".$varAllele;
#	 	
#	 }
	#$no->put($zid,$hvariant);
	
	}
	return 1;
}

sub parse_vcf_line {
	my ($line,$no) = @_;
	chomp($line);
	my @tab = split("\t",$line);
	my $vh = {};
	my $chr =  $tab[0];
	my $pos =  $tab[1];
	my $ref =  $tab[3];
	my $alt =  $tab[4];
	my $rs = $tab[2];
	$rs =~ s/rs//;
	$rs = 0 if ($rs eq ".");
	$vh->{rs} = $rs;
	$vh->{id} = join("-",$chr,$pos,$ref,$alt);
	die() if $ref =~ /,/;
	die() if $ref =~ /;/;
	die() if $alt =~ /,/;
	die() if $alt =~ /;/;
	$vh->{ref} = $ref;
	$vh->{alt} = $alt;
	$vh->{pos} = $pos;
	$vh->{zid} = $no->return_rocks_id($pos,$ref,$alt);
#	$vh->{chr} = ;
#	$vh->{ref} = ;
#	$vh->{alt} = ;
	
	my $info = $tab[7];
	my @pairs = split /;/, $info;
	foreach my $pair (@pairs) {
    		my ($key, $value) = split /=/, $pair;
    		$vh->{infos}->{$key} = $value;
		}
	return $vh;
}


sub compress_vcf_position {
	my ($ref,$alt) = @_;
	my $l1 = length($ref);
	my $l2 = length($alt);
	if ($l1 ==1 && $l2 ==1){
		return $alt;
	}
	if ($l1 ==1 && $l2 > 1){
		return "+".substr($alt, 1);
	}
	elsif ($l1 >1 && $l2 == 1){
		$ref = substr($ref, 1);
		return ($l1 -1);
	}
	elsif ($l1 >1 && $l2 > 1){
		#confess();
		#$ref = substr($ref, 1);
		return $ref."*".$alt;
	}
	
}



sub return_value {
	my ($h,$result,$info,$hinfo) = @_;
	
			  if  ($hinfo->{$info} eq "ID_NOT_FOUND"){
			  	return undef;
			  
			   }
			  else {
			  my $val = $hinfo->{$info}->[0];
			 die() if scalar(@{$hinfo->{$info}}) > 1;
			 return undef unless defined $val;
			 return $val;
			  	 }
}

sub return_rocks_id_from_vcf {
	my ($pos,$ref,$alt) = @_;
	my $l1 = length($ref);
	my $l2 = length($alt);
	return  (sprintf("%010d", $pos)."!".$alt) if ($l1 == 1 && $l2 ==1);
	my $seqid = $alt;
	
	if ($l1 ==1 && $l2 > 1){
		
		$seqid = "+".substr($alt, 1);
		return  (sprintf("%010d", $pos+1)."!".$seqid);
	}
	elsif ($l1 >1 && $l2 ==1){
		$ref = substr($ref, 1);
		$seqid = ($l1 -1);
		return  (sprintf("%010d", $pos+1)."!".$seqid);
	}
	 elsif ($l1 >1 && $l2 == $l1 && $l2>1 ){
	 	
		$ref = substr($ref, 1);
		$alt = substr($alt, 1);
		$seqid = "$ref*$alt";
		return  (sprintf("%010d", $pos+1)."!".$seqid);
	}
	
	else {
		#return  ($pos."!".$ref."*".$alt);
		confess($l1." ".$l2);
	}
	
}
