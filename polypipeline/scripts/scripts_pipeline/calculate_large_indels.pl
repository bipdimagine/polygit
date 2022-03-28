#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
#use lib "$RealBin/../../../lib/obj-nodb/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use Bio::DB::Sam; 
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
#use QueryMoose;
use Tabix;
use Parallel::ForkManager;
use IO::Handle;
use List::Util qw( shuffle sum max min);

my $patient_name;
my $project_name;
my $fork = 1;
use List::MoreUtils qw(firstidx);
use Vcf;

my @tmp_files;
my $dir_out;
my $log_file;
my $vcftools = "/bip-d/soft/distrib/vcftools/latest/bin/";
my $cov;

my $file_out;
GetOptions(
	"log=s" =>\$log_file,
	"patient=s" =>\$patient_name,
	"project=s" =>\$project_name,
	"fork=s" =>\$fork,
);
if ($log_file){
	
	open (STDOUT,">>".$log_file);
	print "echo  ---- start running calculate_large_indels.pl \n";
}



#my $dbh = $buffer->dbh();
#my $sql = qq{select pa.name as patient, p.name as project from PolyprojectNGS.patient pa, PolyprojectNGS.projects p where pa.run_id =? and pa.project_id=p.project_id ;};
#my $sth = $dbh->prepare($sql);
#$sth->execute($run_id);
#my @s = @{$sth->fetchall_arrayref({})};
#foreach my $toto (@s) {
#	my $project_name = $toto->{project};
#	my $patient_name = $toto->{patient};

my $buffer = GBuffer->new();
my $dir_temp =  $buffer->config->{project_pipeline}->{tmp};
my $pm = new Parallel::ForkManager($fork);
my $project = $buffer->newProject( -name => $project_name );
my $patient = $project->getPatient($patient_name);
my $bgzip = $buffer->software("bgzip");
my $tabix = $buffer->software("tabix");
my $samtools = $buffer->software("samtools");
my $family = $patient->family;

mkdir ( $project->getCallingPipelineDir("lifinder")) unless -e  $project->getCallingPipelineDir("lifinder");
my @bams;
my $run = $patient->getRun();
	my $infos = $run->getAllPatientsInfos();
	my @covFiles;
	foreach my $k (@$infos){
		next if $k->{patient} eq $patient_name;
		my $buffertmp = GBuffer->new();
		my $projecttmp = $buffertmp->newProject( -name => $k->{project} );
		my $patienttmp = $projecttmp->getPatient($k->{patient});
		next if ($k->{project} eq $project_name && $patienttmp->family eq $family);
		push(@covFiles,$patienttmp->getCoverageFile());
	#	warn $patienttmp->getCoverageFile();
		push(@bams,$patienttmp->getBamFile());
		$buffertmp->dbh_deconnect();
		$projecttmp = undef;
		$buffertmp = undef;
	}
my @covs = @covFiles;


if (scalar(@covFiles) > 10){
my @covFileRandom = shuffle(@covFiles);

 @covs = @covFileRandom[0..10];
}
my $limit = int(scalar(@covs)/2);
$limit = 1 if $limit == 0;

my $resFile = $patient->getCnvsFile();
warn $resFile;
unlink $resFile if -e $resFile;
unlink $resFile.".tbi" if -e $resFile.".tbi";
my $bed_file = $resFile;
$bed_file =~s/\.gz//;

#my %intspan_chr;
#
# $pm -> run_on_finish (
#    sub {
#      my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $intspan) = @_;
#		$intspan_chr{}
#     
#    }
#  );
 my %results;
$pm->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    		
    				$results{$data->{chr}} =  Set::IntSpan::Fast::XS->new()  unless exists $results{$data->{chr}}  ;
    			$results{$data->{chr}}  = $results{$data->{chr}}->union( $data->{span});
    		
    }
  );
	my $size_window = 1_000_000;
#	$size_window *= 100 if $project->isExome(); 
#	warn $size_window;
my @chr_names = map {$_->ucsc_name} @{$project->getChromosomes()};
	my $real_intspan;
	my $total_jobs =0;
	foreach my $chr (@{$project->getChromosomes()}){
			next if $chr->ucsc_name eq "chrY";
				next if $chr->ucsc_name eq "chrM";
			
		my $windows = $chr->getWindowCaptureForCalling(100,$size_window);
		

		

		warn "compute ".scalar(@$windows);
	foreach my $w (@$windows){

			my $pid = $pm->start and next;
			#my $span = large_indels_by_patient_by_region ($patient->getCoverageFile(),$chr->ucsc_name,$w,\@covFiles) ;
			my $span = large_indels_by_patient_by_region ($patient->getCoverageFile(),$chr->ucsc_name,$w,\@covs,\@bams) ;
			$pm->finish(0,{span=>$span,chr=>$chr->ucsc_name});
    		}
	}
$pm->wait_all_children();

warn "end fork";


open(BED,">$bed_file") or die();
foreach my $chr_name  (@chr_names){
	next unless exists $results{$chr_name};
		
	#next if $chr_name eq 
	my $span = $results{$chr_name};
	next if $span->is_empty;
	
	warn $chr_name;
	my @bed = map{$_ = $chr_name."\t".$_} split(";",$span->as_string({ sep => ";", range => "\t" }));
	my $blimit = int(scalar(@bams)*0.8+1);
	my $nlimit = scalar(@bams) - $blimit;
	foreach my $b (@bed){
		my ($chr,$start,$end) = split(" ",$b);
		$end +=0;
		next if abs($start-$end) < 50;
		next unless $b;
		next if $end ==0;
		my $score = 0;
		my $nb =0;
		my $no= 0;
		warn "$chr:$start-$end";
		foreach my $bam (@bams){
			$nb ++;
			my @res = `$samtools mpileup $bam -r $chr:$start-$end -Q 30 -q 30 | cut -f 4 `;
			my @t = grep {$_ < 7} @res;
			if (scalar(@t) > 5){
				$score ++;
			}
			else {
				$no ++;
			}
			
			last if $no > $nlimit;
		}
		next if $no > $nlimit;
		next if abs($start-$end) < 50;
		next unless $end;
		print BED $b."\t $score \t $nb \n";
	}
}
close BED;
warn $bed_file;
system ("$bgzip $bed_file;$tabix -p bed $bed_file.gz");
exit(0);
#my ($unzipResFile,$ext) = split("\.gz",$resFile);
#warn $unzipResFile;
#
#open (LIFILE,">>$unzipResFile");
#
#foreach my $chr (@{$project->getChromosomes()}){
#	my $chr_name = $chr->ucsc_name;
#	#my $file =$temp_file{$chr}->filename;
#	#my $ucsc = $chr->name();
#	my $file = $dir_temp."/".$patient_name."_".$chr_name;
#	warn $file;
#	next if $chr_name eq "chrM" ;
#
#	open( TITI, $file );
#	while ( my $line = <TITI> ) {
#		
#		print LIFILE $line;
#	}
#	close TITI;
#	#unlink $file;
#}
#run("$bgzip -f $unzipResFile;");
#system("chmod a+w $resFile");
#
#run(" $tabix -f  -s 1 -b 2 -e 2  $resFile");

#sub large_indels_by_patient{
#	my ($project_name,$patient_name,$chr,$dir_temp) = @_;
#	my $buffer = GBuffer->new();
#	my $project = $buffer->newProject( -name => $project_name );
#	my $patient = $project->getPatient($patient_name);
#	my $run = $patient->getRun();
#	my $infos = $run->getAllPatientsInfos();
#	my @covFiles;
#	foreach my $k (@$infos){
#		my $buffertmp = GBuffer->new();
#		my $projecttmp = $buffertmp->newProject( -name => $k->{project} );
#		my $patienttmp = $projecttmp->getPatient($k->{patient});
#		push(@covFiles,$patienttmp->getCoverageFile());
#	}
#	my @otherCovFiles = grep{$_ ne $patient->getCoverageFile()} @covFiles;
#	my $c = $project->getChromosome($chr);
#	my $chr_name = $c->ucsc_name();
#
#	my $all = Set::IntSpan::Fast::XS->new();
#	my $allbis = Set::IntSpan::Fast::XS->new();
#	return $all if $chr_name eq "chrM";
#	my $test =0;
#	my $start = -1;
#	my $end;
#	my %data;
#	my $captureChr = $c->getExtentedGenomicSpan(20);
#	my %dejavu;
#	my @bed = map{$_ = $chr_name.":".$_} split(";",$captureChr->as_string({ sep => ";", range => ":" }));
#	foreach my $of (@otherCovFiles){
#		$test++;
#		#last if $test >5;
#		my $tabixFile = new Tabix(-data => $of);
#		
#	
#		foreach my $b (@bed){
#			next if exists $dejavu{$b} ;
#			my ($ch,$st,$en) = split(":",$b);
#			
#		my $res = $tabixFile->query($ch,$st,$en) ;
#		next unless $res->{_};
#		my $nb=0;
#		my $size = abs($en-$st);
#		while(my $line = $tabixFile->read($res)){
#			#warn $line;
#			my($a,$b,$c) = split(" ",$line);
#			$nb++  if exists $data{$b};
#			next if exists $data{$b};
#			$end = $b;
#			if ($c>=5){
#				$nb++;
#				$data{$b} ++;
#			}
#		
#		
#		}
#		$dejavu{$b} ++ if $nb eq $size;
#		}	
#		 warn $chr_name." patient : $test";
#	}
#	
##	$all->add(keys %data);
#	warn "end $chr_name";
#	
#	my %data_patients;
#	my $start = -1;
#	my $end;
#	my $intSpanPat = Set::IntSpan::Fast::XS->new();
#	my $tabixFile = new Tabix(-data => $patient->getCoverageFile());
#	foreach my $b (@bed){
#		my ($ch,$st,$en) = split(":",$b);
#	my $res = $tabixFile->query($ch,$st,$en) ;
#	next unless $res->{_};
#	while(my $line = $tabixFile->read($res)){
#		
#		my($a,$b,$c) = split("\t",$line);
#		if (exists $data{$b}){
#				$intSpanPat->add($b) if $c < 3;
#		}
#		
#	}
#	}
#	mkdir ( $project->getCallingPipelineDir("lifinder")) ;
#
#	return $intSpanPat;
#	
##	my $final;
##	$final = Set::IntSpan::Fast::XS->new();
##	$final= $all->diff($intSpanPat);
#
#	
#	
#}

sub large_indels_by_patient_by_region {
	my ($patient_cov,$chr_name,$window,$otherCovFiles) = @_;
	my $all = Set::IntSpan::Fast::XS->new();
	my $allbis = Set::IntSpan::Fast::XS->new();
		
	return $all if $chr_name eq "chrM";
	my $htabix;
		

	 
    		
    		 my $iter = $window->{intspan}->iterate_runs();
	 
    while (my ( $from, $to ) = $iter->()) {
    	my %tt;
    	my %xx;
    		my $region = Set::IntSpan::Fast::XS->new("$from-$to");
    		my $nb =0;
    		my $all_region = Set::IntSpan::Fast::XS->new();
    		my $bad_region = Set::IntSpan::Fast::XS->new();
    		my $viewed_region = Set::IntSpan::Fast::XS->new();
    		foreach my $of (@$otherCovFiles){
		#last if $test >5;
		
			
			$htabix->{$of} = new Tabix(-data => $of) unless exists $htabix->{$of};

			my $tabixFile  = $htabix->{$of};	
    			
    			my $tintspan = $region->diff($viewed_region);
			last if $tintspan->is_empty;
			if (scalar($tintspan->as_array) < 50 ){
				#last;
			}
			
			#warn "\t".scalar($tintspan->as_array)." ".$nb++." ".$tintspan->as_string." ".scalar($all_region->as_array);
			#die() if scalar($tintspan->as_array) < 10;
			#warn "start :".scalar($tintspan->as_array);
		
		 my $iter2 = $tintspan->iterate_runs();
	 		my $string ="";
   			 while (my ( $from2, $to2 ) = $iter2->()) {
   			 #	next if abs($from2-$to2) <100;
   			 	#$string .= $chr_name.":$from2-$to2 ";
   			 	my $res = $htabix->{$of}->query($chr_name,$from2,$to2) ;
   			 	while(my $line = $tabixFile->read($res)){
   			 			my($a,$b,$c) = split(" ",$line);
   			 				if ($c >= 10 ) {
			 					$tt{$b} ++;
			 		}
			 		if ($c <= 2){
			 			$xx{$b} ++;
			 			
			 		}
   			 	}
   			 }
   			 
   			# last if $string eq "";
#   			 warn $string." ".$nb++;
#			my $start = $from;
#			my $end = $to ;
#			#$chr_name:$start-$end;
#			
#			open(TABIX,"tabix $of $string | ");
#				my $nb=0;
#				while(<TABIX>) {
#					my($a,$b,$c) = split(" ",$_);
#					#next if  $tt{$b} >2;
#			 		if ($c >= 10 ) {
#			 			$tt{$b} ++;
#			 		}
#			 	
#				}
#				close TABIX;
				
			my @t = sort { $a <=> $b} grep {$tt{$_} >$limit} keys %tt;
			my @x = sort { $a <=> $b} grep {$xx{$_} >$limit} keys %xx;
			$all_region->add(@t) ;#if $first;
			$bad_region->add(@x) ;
			$viewed_region =  $viewed_region->union($all_region);
			$viewed_region =  $viewed_region->union($bad_region);
			#$all_region->add(@x) ;#if $first;
			
		 }
		# warn "end";
		 $all = $all->union($all_region);
		 $all = $all->diff($bad_region);
    }
		 
    		#warn $all_region->as_string();
  
    		
   	#	die();
#	$all->add(keys %data);
	 my $iter = $all->iterate_runs();
	 my $intSpanPat = Set::IntSpan::Fast::XS->new();
	my $tabixFile = new Tabix(-data => $patient->getCoverageFile());
	 
    while (my ( $from, $to ) = $iter->()) {
    			my $res = $tabixFile->query($chr_name,$from,$to) ;
    			unless ($res->{_}) {
    					$intSpanPat->add_range($from,$to);
    					next;
    			}
    			while(my $line = $tabixFile->read($res)){
    					my($a,$b,$c) = split(" ",$line);
    					next if $c>=2;
    					#	if ($all->contains($b)){
    							$intSpanPat->add($b);
    					#	}
    			
    }
	

	}

	
#warn scalar(keys %tt);
	return $intSpanPat;
	
#	my $final;
#	$final = Set::IntSpan::Fast::XS->new();
#	$final= $all->diff($intSpanPat);

	
	
}


sub run {
	my ($cmd) = @_;
	my $return = system($cmd);
	if ($return ne 0){
		confess("error : $cmd");
	}
}
