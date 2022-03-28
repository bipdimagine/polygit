#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use Digest::MD5::File qw (file_md5_hex);

use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages/";
use colored;
#use Set::IntSpan;
use GBuffer;
use Data::Dumper;
use Getopt::Long;
use Carp;
use calling_target;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use String::ProgressBar;
use List::Util qw(sum);
use Path::Tiny;
use Digest::MD5::File qw (file_md5_hex);
use Set::IntSpan::Island;
use Set::IntSpan;
#use Tree::R;
#use Bio::DB::HTS;
use Bio::DB::Sam;
use Set::IntSpan;
use Path::Tiny;
require "$Bin/../lib/interval_util.pm"; 
 my $project_name;
 my $fork;
 my $fileout;
 my $bam;
 my $patient_name;
my $limit;
my $type;
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"limit=s"=>\$limit,
	"type=s"=>\$type,
);
die("type=count,sample") unless $type;
die() unless $limit;
#my $types = ["1","5","10","15"];
my $buffer = GBuffer->new();

my $project = $buffer->newProject( -name => $project_name );
my $patients = $project->getPatients();
my @cmds;
my $pm1 = new Parallel::ForkManager($fork);
my $lmdb_patients = [];
my $lmdb_patients = list_patients();
warn Dumper $lmdb_patients;
 project_count($project_name) if $type eq "count";
 project_samples($project_name) if $type eq "sample";
 
sub list_patients {
	open(FILE,"list_patient.txt");
	my $lmdb_patients;
	while(<FILE>){
		warn $_;
	chomp();
	my ($patient_name,$project_name) = split(" ",$_);
		my $root = "/data-isilon/sequencing/ngs/"."$project_name/HG19/align/callable";
		die($root."/".$patient_name.".ok.callable") unless -e $root."/".$patient_name.".ok.callable";
#		warn $project->getCoverageCallable()."/".$patient->name.".callable";
		push(@$lmdb_patients,{name=>$patient_name,id=>$patient_name,dir=>$root});
	}
	return $lmdb_patients;
}

sub project_samples {
	
	my $project_name ="toto";
	my $dir =".";
	my @ids = sort{$a <=> $b} map{$_->{name}} @$lmdb_patients; 
	
   my $full_ids = join(";" , sort{$a <=> $b} @ids);
   
	my $pm1 = new Parallel::ForkManager($fork);
	## finish concat intspan
	my $first =1;
#	 $pm1->run_on_finish(
#    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
#		unless (defined($data) or $exit_code>0){
#				print qq|No message received from child process $exit_code $pid!\n|; 
#				#die();
#				return;
#			}
#			my $mode ="w";
#			$mode = "c" if $first;
#			
#			$first = undef;
#			
#			my $lmdb2= GenBoNoSqlLmdb->new(dir=>$dir,mode=>$mode,name=>$project_name.".lmdb.callable",is_compress=>1);
#			my $key = $data->{key};
#			my $v = $data->{value};
#			foreach my $k (keys %$v){
#				$lmdb2->put($k,$v->{$k});
#			}
# 			$lmdb2->close();
#    }
#  );
  
  ### principal loop 

  	
my $chrs = [1..22,'X','Y','MT'];
#my $chrs = [1];
#$chrs = [20,21];

# $types = ["15"];
#$chrs = [1];
my $t = time;

	my $t = $limit;
foreach my $chr (@$chrs){

#$type =["1"];
	
	 	my $key = 	"chr".$chr."_".$t;
	  	 my $pid = $pm1->start and next;
	 	  	my $res = callable_project($lmdb_patients,$key,$chr);
	 		my $data;
	 		$data->{key} = $key; 
	 		#$data->{value} = $res; 
	  		$pm1->finish(0,$data);
	}
	 $pm1->wait_all_children();
	
my $f = "$dir/all.bed";
 foreach my $chr (@$chrs){
 	my $key = 	"chr".$chr."_".$t;
 	system("cat /tmp/$key.bed >>/tmp/$limit.bed ; rm   /tmp/$key.bed;") 
 }
  my $bgzip = $buffer->software("bgzip");
 my $tabix = $buffer->software("tabix");
  system("$bgzip -@ $fork  -f /tmp/$limit.bed  && mv  /tmp/$limit.bed.gz $limit.gz");
 
 exit(0);


  
}


sub project_count {
	
	my $project_name ="toto";
	my $dir =".";
	my @ids = sort{$a <=> $b} map{$_->{name}} @$lmdb_patients; 
	
   my $full_ids = join(";" , sort{$a <=> $b} @ids);
   
	my $pm1 = new Parallel::ForkManager($fork);
	## finish concat intspan
	my $first =1;

  ### principal loop 

  	
my $chrs = [1..22,'X','Y','MT'];
#my $chrs = [1];
#$chrs = [20,21];

# $types = ["15"];
#$chrs = [1];
my $t = time;

	my $t = $limit;
foreach my $chr (@$chrs){

#$type =["1"];
	
	 	my $key = 	"chr".$chr."_".$t;
	  	 my $pid = $pm1->start and next;
	 	  	my $res = count_patient($lmdb_patients,$key,$chr);

	  		$pm1->finish(0,[]);
	}
	 $pm1->wait_all_children();
	
my $f = "$dir/all.bed";
 foreach my $chr (@$chrs){
 	my $key = 	"chr".$chr."_".$t;
 	system("sort -k2,2n -k3,3n /tmp/$key.bed >>/tmp/$limit.count.bed ; rm   /tmp/$key.bed;") 
 }
  my $bgzip = $buffer->software("bgzip");
 my $tabix = $buffer->software("tabix");
  system("$bgzip -@ $fork  -f /tmp/$limit.count.bed  && mv  /tmp/$limit.count.bed.gz $limit.count.gz");
 
 exit(0);


  
}

sub count_patient {
	my ($patients,$key,$chr) = @_;
	warn $key;
	my $regions;
	my $union = Set::IntSpan::Fast::XS->new();
	foreach my $patient (@$patients){
		my $lmdb = GenBoNoSqlLmdb->new(dir=>$patient->{dir},mode=>"r",name=>$patient->{name}.".ok.callable",is_compress=>1);
		my $a = $lmdb->get($key);
		next unless $a;
		$union = $union->union($a);
		my $region = interval_util::construct_regions_from_intspan($a,$patient->{id});
		push(@$regions, @$region) if $region;
		#warn $a;
	}
	
	$regions =  interval_util::coverages_regions3($regions,$key,$chr);
	return ;
}



sub callable_project {
	my ($patients,$key,$chr) = @_;
	warn $key;
	my $regions;
	my $union = Set::IntSpan::Fast::XS->new();
	foreach my $patient (@$patients){
		my $lmdb = GenBoNoSqlLmdb->new(dir=>$patient->{dir},mode=>"r",name=>$patient->{name}.".ok.callable",is_compress=>1);
		my $a = $lmdb->get($key);
		next unless $a;
		$union = $union->union($a);
		my $region = interval_util::construct_regions_from_intspan($a,$patient->{id});
		push(@$regions, @$region) if $region;
		#warn $a;
	}
	warn "step 1";
	$regions =  interval_util::coverages_regions2($regions,$key,$chr);
	warn "step 2";
	
	#my $aunions = interval_util::construct_regions_from_intspan($union);
	#warn "step 3";
	
	
	return ({$key=>$regions,$key."_union"=>$union});
}


