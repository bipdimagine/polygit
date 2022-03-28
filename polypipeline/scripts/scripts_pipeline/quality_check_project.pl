#!/usr/bin/perl
use CGI qw/:standard :html3/;
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use colored;
use lib $Bin;
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use List::MoreUtils qw(part);
use File::Temp;
use JSON::XS;
use check_utils; 
use Statistics::Zscore; 

my $projectName;
my $end_ext = "uni";
my $details;
my $fork = 1;
my $vcf_file;
my $log_file;
my $noclean;
my $chr_names;
my $family_test;
GetOptions(
	'project=s'		=> \$projectName,
	"details=s"		=> \$details,
	"fork=s"		=> \$fork,
	"vcf_file=s" => \$vcf_file,
	 "log=s" => \$log_file,
	  "noclean=s" => \$noclean,
	  'chr=s' =>\$chr_names,
	  'family=s' => \$family_test,
);


my $date = `date`;

chomp($date);
if ($log_file){
	open (STDOUT,">>".$log_file);
}

colored::stabilo('blue',"START QUALITY CHECK ");
die("\n\nERROR: -project option mandatory. Exit...\n\n") unless ($projectName);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
my $ped_file = $project->getPedigreeFile();

my $gatk  = $buffer->getSoftware('java')." -jar ".$buffer->getSoftware('gatk');
my $plink = $buffer->getSoftware('plink');


my $patients = $project->getPatients();
foreach my $p (@$patients){
	print $p->name();
	print "\n";
}

my @chrs_plink = (21..22,'X','Y');
if ($chr_names eq "all"){
	
	 @chrs_plink = (1..22,'X','Y');
}

my @lChr = ('21', 'X', 'Y');
my $children;
my $family;
 $family = $project->families();
eval { $family = $project->families(); };
if ($@) { die("\n\nERROR: can't construct family relations... Maybe a problem in PED file. Die...\n\n"); }
my $statistics= {};
foreach my $p (@{$project->getPatients}){
	next() unless $p->family eq $family_test;
	my $name = $p->name;
	$statistics->{$name}->{indb} = 1;
	$statistics->{$name}->{family} = $p->family;
	$statistics->{$name}->{relation_code} = 0;
	$statistics->{$name}->{relation_code} = 1 if $p->isMother;
	$statistics->{$name}->{relation_code} = 2 if $p->isChild;
	$statistics->{$name}->{sex} =  $p->sex;
	
	
}
if ($project->isFamilialStudy){
	#die("\n\nERROR: no pedigree file found. Exit...\n\n") unless (-e $ped_file);
	colored::stabilo("blue","check pedigree");
	#die() unless check_utils::printErrorPedigree($project);
	colored::stabilo("blue","check mendelian ");
}
else {
	 colored::stabilo("magenta", " YOU DON'T PROVIDE ANY PEDIGREE FILE SKIP PLINK ....");
	 
}
if (-e $vcf_file){
	colored::stabilo("blue","check dbsnp");
	 check_dbsnp_with_no_split($project,$vcf_file,$statistics,$family_test) ;
	 }
else {
	$vcf_file = concatVcf_by_patient();
	colored::stabilo("blue","check dbsnp");
	#check_dbsnp($project,$statistics);
}	 
if ($project->isFamilialStudy){
	fast_plink($vcf_file,$statistics);
}

colored::stabilo("blue","check sex ");
check_sex($project,$statistics);
colored::stabilo("blue","check coverage ");
check_coverage($project,$statistics);

	
foreach my $k (keys %$statistics){
	delete $statistics->{$k} unless exists $statistics->{$k}->{indb} ;
}
colored::stabilo("blue","check cache");
checkCache($project,$statistics);
check_utils::printTable($statistics,$project->name());


 open( JSON, ">".$project->getStatsFile );
 print JSON encode_json( $statistics );
 close(JSON);
 warn "coucou";
#warn $project->getStatsFile;
my $dir_out = $project->getCallingPipelineDir("unifiedgenotyper").'/plink';
unless ($noclean){
	system("rm $dir_out/*"); 
	system ("rmdir $dir_out");
}
else {
	warn $dir_out;
}

#warn $dir_out if $noclean; 
sub check_sex {
		my ($project,$statistics) = @_;
		
	foreach my $patient (@{$project->getPatients()}){
	my $csex = $patient->compute_sex();
#	my $bam = $patient->getBamFile();
#	
#	my @apos = ("chrY:6736371-6736371","chrY:6736443-6736443","chrY:6737845-6737845");
#	my $samtools = $project->buffer()->software("samtools");
#	my $cmd = qq{$samtools depth $bam };
#	
#	foreach my $pos (@apos){
#		my $cmd2 = $cmd." -r $pos | cut -f 3 ";
#		my ($depth) = `$cmd2`;
#		chomp($depth);
#		$depth = 0 if ($depth eq '');
#		my $sex;
#		$sex= 1 if $depth > 10;
#		$sex =2 if $depth <=10;
		$statistics->{$patient->name}->{check_sex}->{bipd}->{sex} = $csex;
		$statistics->{$patient->name}->{check_sex}->{bipd}->{score} = "-";
#	}
	}
}
sub check_coverage {
		my ($project,$statistics) = @_;
foreach my $patient (@{$project->getPatients}){
	my $patName = $patient->name();
	my $hashCov = $patient->coverage();
	$statistics->{$patient->name}->{coverage}->{mean} = $hashCov->{'mean'};
	$statistics->{$patient->name}->{coverage}->{5} = $hashCov->{'5x'};;
	$statistics->{$patient->name}->{coverage}->{15} = $hashCov->{'15x'};
	$statistics->{$patient->name}->{coverage}->{30} = $hashCov->{'30x'};
	}
	
}


sub check_dbsnp {
	my ($project,$statistics) = @_;
	my @indels;
	my @snps;
	
	foreach my $patient (@{$project->getPatients()}){
		my $name = $patient->{name};
	  	my $vcf_snp = $patient->getVariationsFile("unifiedgenotyper");
		die() unless -e $vcf_snp;	
	
	my $debug;
	 $debug =1  if $patient->name eq "PEL_P";
	my $cmd1 = qq{ zgrep -v "#" $vcf_snp |grep -v "0/0" | grep -v "\\./\\." | };
	open ("TOTO", $cmd1 );
	my $nb = 0;
	my $rs = 0;
	my $last_chr;
	my $real_nb;
	while (my $line = <TOTO>){
		
		$rs ++ if $line =~ /rs/;
		
	 	my @t = split(" ",$line);
	 	my $s = $t[4];
	 	my @sp = split(":",$t[9]);
	 	my %all;
	 	map {$all{$_}++ if $_>0} sort{$a <=> $b} split("/",$sp[0]);
	  $nb += scalar(keys %all);
	  
	 	#$real_nb ++ if scalar(split(",",$s)) >1;

	}

close (TOTO);
warn $rs if $debug;

my $pourcent =0;
$pourcent = int( ($rs/$nb) *10000) if $nb >0; 
 $pourcent =  int($pourcent/100);
 $statistics->{$name}->{snp}->{nb} = $nb;
  $statistics->{$name}->{snp}->{real_nb} = $real_nb;
 $statistics->{$name}->{snp}->{dbsnp} = $rs;
 $statistics->{$name}->{snp}->{percent} = $pourcent;
 push(@snps,$nb);
 my $vcf_indel = $patient->getIndelsFile();
die($patient->name) unless -e $vcf_indel;	
	
	
	
	my $cmd1 = qq{ zgrep -v "#" $vcf_indel |grep -v "0/0" | grep -v "\\./\\." | };
	open ("TOTO", $cmd1 );
	 $nb = 0;
	 $rs = 0;
	while (my $line = <TOTO>){
		$rs ++ if $line =~ /rs/;
			my @t = split(" ",$line);
	 	my $s = $t[4];
		my @sp = split(":",$t[9]);
	 	my %all;
	map {$all{$_}++ if $_>0} sort{$a <=> $b} split("/",$sp[0]);
	  $nb += scalar(keys %all);

	}

close (TOTO);


 $pourcent =0;
$pourcent = int( ($rs/$nb) *10000) if $nb >0; 
 $pourcent =  int($pourcent/100);
 $statistics->{$name}->{indel}->{nb} = $nb;
 $statistics->{$name}->{indel}->{dbsnp} = $rs;
 $statistics->{$name}->{indel}->{percent} = $pourcent;
 push(@indels,$nb);

 
 
 
}
	my $z = Statistics::Zscore->new;
	my $zscore_snps = $z->standardize( \@snps );
	my $zscore_indels = $z->standardize( \@indels );
	my $i=0;
	foreach my $patient (@{$project->getPatients()}){
		my $name = $patient->name();
		  $statistics->{$name}->{snp}->{zscore} = int($zscore_snps->[$i]*100)/100;
		   $statistics->{$name}->{indel}->{zscore} = int($zscore_indels->[$i]*100)/100;
		   $i++;
		  
	}
	
}


sub check_dbsnp_with_no_split {
	my ($project,$vcf_file,$statistics,$family_test) = @_;
	my @indels;
	my @snps;
	my $no = $project->lite_public_snps();
	my $hrs;
	
	my @lines = ` zgrep -v "#" $vcf_file `;
	
	my $hindex;
	my $hhindex;
	my $index;
		my $line_header = `zgrep "#CHR" $vcf_file`;
	foreach my $patient (@{$project->getPatients()}){
		next() unless $patient->family eq $family_test;
		my $name = $patient->name;
		#warn $name;
	
		chomp($line_header);
		my @theader =  split(" ",$line_header);
		my $index=0;
		for (my $i=8;$i< @theader;$i++){
			$index = $i if $theader[$i] eq $name;
			
		}
		$hindex->{ $patient->name} = $index;
		$hhindex->{$index} = $patient->name;
		if  ($index ==0){
			my $string;
			for (my $i=8;$i< @theader;$i++){
				$string.=$theader[$i].";"; 
			}
			print " ERROR NoT FOUND $name in vcf\n";
			print $string."\n";
			warn $string;
	}
		confess() unless $index;
		$index ++;
	}

	
#	exit(0);
my $stat;
my $n=0;
foreach my $line (@lines){
	chomp($line);
	my @t = split(" ",$line);
	my $rs = 0 ;
	my $snp = 0 ;
	my $indel = 0;
	#warn $n."/".scalar(@lines);
	#last if $n > 50000;
	$n++;
	if (length($t[3]) eq length($t[4])){
		my $id = $t[0]."_".$t[1]."_".$t[3]."_".$t[4];
		$snp =1;
		my $chr = $project->getChromosome($t[0]);
		my $z = $no->get($chr->name,$t[1]);
		$rs = 1  if exists $z->{$t[4]};
	}
	else {
		$indel = 1;
	}
	for  (my $i=9;$i<@t;$i++){
		
		my $name = $hhindex->{$i};
		my $code = $t[$i];
		next if $code =~ /0\/0/;
		next if $code=~/\.\/\./;
		$statistics->{$name}->{snp}->{nb} += $snp;
		$statistics->{$name}->{snp}->{dbsnp} += $rs;
		$statistics->{$name}->{indel}->{nb} += $indel;
		
	}
	
	
	
}
my @snps;
my @indels;

foreach my $patient (@{$project->getPatients()}){
	my $name = $patient->name();
	push(@snps,$statistics->{$name}->{snp}->{nb});
	push(@indels,$statistics->{$name}->{indel}->{nb});
	my $pourcent =0;
	 $pourcent = int( ($statistics->{$name}->{snp}->{dbsnp}/$statistics->{$name}->{snp}->{nb}) *10000) if $statistics->{$name}->{snp}->{nb} >0; 
	$statistics->{$name}->{snp}->{percent} =  int($pourcent/100);
}
	my $z = Statistics::Zscore->new;
	my $zscore_snps = $z->standardize( \@snps );
	my $zscore_indels = $z->standardize( \@indels );
	my $i=0;
	foreach my $patient (@{$project->getPatients()}){
		my $name = $patient->name();
		  $statistics->{$name}->{snp}->{zscore} = int($zscore_snps->[$i]*100)/100;
		   $statistics->{$name}->{indel}->{zscore} = int($zscore_indels->[$i]*100)/100;
		   $i++;
	}



#	
#		my $cmd1 = qq{ zgrep -v "#" $vcf_file | cut -f 1-5,$index |grep -v "0/0" | grep -v "\\./\\." | };
#	
#		open ("TOTO", $cmd1 );
#		my $nb = 0;
#		my $rs = 0;
#		my $last_chr;
#		my $nb_snp = 0;
#		my $nb_indel =0;
#		
#		while (my $line = <TOTO>){
#			my @t = split(" ",$line);
#			if (length($t[3]) eq length($t[4])){
#				$nb_snp ++;
#				my $id = $t[0]."_".$t[1]."_".$t[3]."_".$t[4];
#			
#				if (exists $hrs->{$id}){
#					$rs += $hrs->{$id};
#					warn "*";
#				}
#				else {
#				warn $nb;
#				my $chr = $project->getChromosome($t[0]);
#				my $z = $no->get($chr->name,$t[1]);
#				$hrs->{$id} =0;
#				if ($z){
#				$rs ++ if exists $z->{$t[4]};
#				$hrs->{$id} =1;
#				}
#				}
#
#			#	my $n = $project->getNbDejaVu($id);
#				#$id =~ s/chr//;
#				 #$n = $project->getNbDejaVu($id);
#				 #	$rs ++ if $n > 0;
#			}
#			else {
#				$nb_indel ++;
#			}
#			
#			$nb++;
#			
#			#$rs ++ if $line =~ /rs/;
#			#last if $nb > 8000;
#			#last if $line =~ /chr2/;
#		}
#	close (TOTO);
#
#	my $pourcent =0;
#	$pourcent = int( ($rs/$nb) *10000) if $nb >0; 
#	
#	$pourcent =  int($pourcent/100);
# 	$statistics->{$name}->{snp}->{nb} = $nb_snp;
# 	$statistics->{$name}->{snp}->{dbsnp} = $rs;
# 	$statistics->{$name}->{snp}->{percent} = $pourcent;
# 	$statistics->{$name}->{indel}->{nb} = $nb_indel;
# 	$statistics->{$name}->{indel}->{dbsnp} = $rs;
# 	$statistics->{$name}->{indel}->{percent} = $pourcent;
# 	push(@indels,$nb);
#	push(@snps,$nb);
#	}#end_patient

#	my $z = Statistics::Zscore->new;
#	my $zscore_snps = $z->standardize( \@snps );
#	my $zscore_indels = $z->standardize( \@indels );
#	my $i=0;
#	foreach my $patient (@{$project->getPatients()}){
#		my $name = $patient->name();
#		  $statistics->{$name}->{snp}->{zscore} = int($zscore_snps->[$i]*100)/100;
#		   $statistics->{$name}->{indel}->{zscore} = int($zscore_indels->[$i]*100)/100;
#		   $i++;
#		  
#	}
}



sub checkPedFile {
	my $ped_file;
	unless ($ped_file) { $ped_file = $project->getPedigreeFile(); }
	warn $ped_file;
	die() unless -e $ped_file;
	open (PED, $ped_file);
	
	while (<PED>) {
		chomp($_);
		my $line = $_;
		$line =~ s/\s+/\t/g;
		unless  ($line =~ /\t/) { warn "======= WARN: no tabulation found in this line. Maybe a problem ???\n\nLINE: $line\n\n"; }
		my @lFields = split("\t", $line);
		my $nbCol = scalar(@lFields);
		if ($nbCol < 6) {
			warn  "======= ERROR: line with only $nbCol columns in ped file ! 6 columns are mandatory... Die !\n\nLINE: $line\n\n";
			warn Dumper @lFields;
			die();
		}
	}
	close(PED);
}
 



sub concatVcf_by_patient {
	print "\n\n### Step 2: concat VCF files.\n";
	my $dir_out = $project->getCallingPipelineDir("unifiedgenotyper").'/plink';
	if (-e $dir_out){
		system("rm $dir_out/*");
	}
	else {
	 	system("mkdir $dir_out");
	}
	my $fileout = $dir_out."/".$project->name.".merge.vcf";
	my $fileout1 = $dir_out."/".$project->name.".merge1.vcf";
	my $bed = $dir_out."/".$project->name.".bed";
	unlink $fileout if -e $fileout;
	unlink $bed  if -e $bed;
	my @vcfs;
	foreach my $p (@{$project->getPatients}){
		next() unless $p->family eq $family_test;
		my $vcf = $p->getVariationsFile("unifiedgenotyper");
		my $vcf = $p->getVariationsFile("haplotypecaller")  unless -e $vcf;
		die($vcf."-".$p->name) unless -e $vcf;
		push(@vcfs,$vcf);
	}
	## create bed file
	open(BED,">$bed"); 
	foreach my $chr_name (@chrs_plink){
		my $chr = $project->getChromosome($chr_name);
		my $span = $chr->getIntSpanCapture();

		# LINK 1 PAR REGIONS: http://www.ensembl.org/info/genome/genebuild/assembly.html
		#	# LINK 2 PAR REGIONS: http://blog.kokocinski.net/index.php/par-regions?blog=2
	if ($chr->name() eq 'X') {
		my $spanX = Set::IntSpan::Fast::XS->new();
			$spanX->add_from_string('60001-2699520, 154931044-155270560');
			$span = $span->diff($spanX);
		}
			elsif ($chr->name() eq 'Y') {
			my $spanY = Set::IntSpan::Fast::XS->new();
			$spanY->add_from_string('10001-2649520, 59034050-59034050');
			$span = $span->diff($spanY);
			}

		my @bed = map{$_ = $chr->ucsc_name."\t".$_} split(";",$span->as_string({ sep => ";", range => "\t" }));
		print BED join("\n",@bed);
		print BED "\n";
	}
	
	
	close BED;
	
	my $variant_string = join(" --variant ",@vcfs);
	#my $cmd = $gatk." -R ".$project->getGenomeFasta()." -nt 8 -T CombineVariants  --variant $variant_string -o $fileout "." -L ". $bed." 2>/dev/null >/dev/null" ;
	my $cmd = $gatk." -R ".$project->getGenomeFasta()." -nt 8 -T CombineVariants  --variant $variant_string -o $fileout  2>/dev/null >/dev/null" ;
	my $variant_string = join(" ",@vcfs);
	my $cmd = " bcftools merge   $variant_string -o $fileout  2>/dev/null >$fileout " ;
	system($cmd);
	#system ("grep -v './.'  $fileout1 >$fileout;rm $fileout1" );
	
	return $fileout;
}

sub fast_plink {
	my ($vcf_file,$statistics) = @_;
	my $dir_out = $project->getCallingPipelineDir("unifiedgenotyper").'/plink';
	

	
	print "\n\n### Step 3: launching PLINK.\n";
	my $dir = $project->getCallingPipelineDir("unifiedgenotyper");
	
	my $logPlink = $dir.'/plink/'.$projectName.'.plink.resume';
	
	my @snps;
	if ($vcf_file =~/gz/){
	open (VCF, "zcat $vcf_file | ");
	}
	else {
		open (VCF, "cat $vcf_file | ");
	}
	my @samples;
	my $DP;
	while (my $line = <VCF>){
		next if $line =~ /##/;
		chomp($line);
		if ($line =~ /#CHR/){
			my @tdata = split(" ",$line);
			@samples = splice(@tdata,9);#$tdata[9..-1];
			next;
		}
		my ($chrom,$pos,$pid,$ref,$alt,$qual,$filter,$info,$format,@gsamples) = split(" ",$line); 
		next if $ref =~ /,/;
		next if $alt =~ /,/;
		my $snp;
		my $id = $chrom."_".$pos."_".$alt;
		$snp->{id} = $id;
		$snp->{position} = $pos;
		$chrom =~ s/chr//;
		$chrom = "MT" if $chrom eq "M";
		$snp->{chr} = $chrom;
		unless ($DP){
			my @t = split(":",$format);
			for (my $i =0;$i<@t;$i++){
				$DP=$i if $t[$i] eq "DP";
			}
			$DP =-1 unless $DP;
		}
		for (my $i=0;$i<@samples;$i++){
			my $name = $samples[$i];
			my $string = $gsamples[$i];
			my $geno;
			if ($string =~ /0\/0/) { $geno ="$ref $ref"; }
			elsif ($string =~ /0\/1/) { $geno ="$ref $alt"; }
			elsif ($string =~ /1\/0/) { $geno ="$ref $alt"; }
			elsif ($string =~ /1\/1/) { $geno ="$alt $alt"; }
			elsif ($string =~ /.\/./) { $geno ="0 0"; }
			elsif ($string =~ /.:./) { $geno ="0 0"; }
			else { die($string); }
			if ($DP > 0){
				my @tt = split(":",$string);
				my $v = $tt[$DP];
				$geno = "0 0" if $v < 10;
								
			}
			$snp->{samples}->{$name} = $geno;
		}
		
		foreach my $f (values %{$project->families()}) {
			my ($find) = grep {$snp->{samples}->{$_} eq "0 0" } @{$f->{members}};
			if ($find){
				map {$snp->{samples}->{$_} =  "0 0" } @{$f->{members}};
			}	
			else {
				$f->{nb_snp} ++;
			}
		}
		#my ($find) = grep {$snp->{samples}->{$_}  eq "0 0"} @samples;
		#next if $find;
		push(@snps,$snp);
	}
	
	$dir .= "/plink/";
	mkdir $dir unless -e $dir;
	my $ped_file = $dir."/".$project->name.".ped";
	my @samples_name;
	open (PED,">$ped_file");
	
		foreach my $p (@{$project->getPatients()}){
			print PED $p->pedigreeLine()."\n";
			push(@samples_name,$p->name);
			
		}
	close PED;	
	
	
	die() unless -e $ped_file;
	
	
	
	my $tped_file = $dir."/".$project->name.".tped";
	open (TPED,">".$tped_file);
	foreach my $snp (@snps){
		my @data;
		push(@data,$snp->{chr});
		push(@data,$snp->{id});
		push(@data,0);
		push(@data,$snp->{position});
		foreach my $name  (@samples_name){
			
			if (exists $snp->{samples}->{$name} ) { push(@data,$snp->{samples}->{$name}); }
			else { push(@data,$snp->{samples}->{$name}, "0 0");	}
		}
		print TPED join("\t",@data)."\n";
	}
	close TPED;
	my $cmd2 = "$plink --tped $tped_file --tfam $ped_file --noweb --mendel  --out $dir/$projectName";
	warn $cmd2;
	my @log = `$cmd2`;
	my $file1 = "$dir/".$project->name();
	open(MENDEL,"$file1.imendel") || die("problem with $file1.imendel \n$cmd2");
	while(my $line =<MENDEL> ){	
		chomp($line);
		next if $line =~ /FID/;
		my ($f, $name,$nb) = split(" ",$line);
		
		my $nb_snp = $project->families()->{$f}->{nb_snp};
		next if $nb_snp == 0;
		my $p = int($nb/$nb_snp*10000)/100;
		$statistics->{$name}->{mendelian_errors}->{percent} = $p;
			$statistics->{$name}->{mendelian_errors}->{nb_snp} = $project->families()->{$f}->{nb_snp};
	}
	close(MENDEL);
	my $cmd3 = "$plink --tped $tped_file --tfam $ped_file  --noweb  --check-sex --out $dir/$projectName 1>$logPlink 2>>$logPlink";
	warn $cmd3;
	`$cmd3`;
	my $sexFile = $project->getCallingPipelineDir("unifiedgenotyper").'/plink/'.$projectName.'.sexcheck';
	warn $sexFile;
	unless (-e $sexFile) { die("\n\nERROR: PLINK sex file doesn't exist... Die...\n\n") };
	open(SEX, "$sexFile");
	while(<SEX>) {
		chomp($_);
		my $line = $_;
		next if($line =~ /PEDSEX/);
		$line =~ s/\s+/ /g;
		my @lFields = split(" ", $line);
		my $name = $lFields[1];
		$statistics->{$name}->{check_sex}->{plink}->{sex} = $lFields[3];
		#$statistics->{$name}->{check_sex}->{plink}->{snp} = $lFields[2];
		$statistics->{$name}->{check_sex}->{plink}->{status} = $lFields[4];
		$statistics->{$name}->{check_sex}->{plink}->{score} = int($lFields[5]*100)/100;
	}
	close(SEX);
}





sub checkCache {
	my ($project,$statistics) = @_;
	 my $dir_cache  = $project->getCacheDir();
	 foreach my $patient (@{$project->getPatients}){
	 		my $name = $patient->name();
	 	 	my $file_kct = $dir_cache."/".$patient->name().".kct";
	 	 	my $mtime_kct = (stat $file_kct)[9];
	 	 	unless (-e $file_kct){
	 	 	 $statistics->{$name}->{cache}->{date} = 999;
	 	 	next;
	 	 	}
	 	#my $cpt = `/bip-d/soft/bin/get.pl $file_kct | wc -l`;
	 	my $cmd= qq{get.pl $file_kct | cut -f 1 | perl -lane ' \@t = split("_",\$_);print \$_ if (length(\$t[2]) eq length(\$t[3]));' | wc -l };
	 	my $cpt = `$cmd`;
	 	chomp($cpt);
	 	my $cmd_indel= qq{get.pl $file_kct | cut -f 1 | perl -lane ' \@t = split("_",\$_);print \$_ if (length(\$t[2]) ne length(\$t[3]));' | wc -l};
	 	my $cpt_indel = `$cmd_indel`;
	 	chomp($cpt_indel);
	 	my $mtime_vcf = (stat $patient->getVariationsFile())[9];
	 	my $delta = $mtime_kct - $mtime_vcf;
	 	$statistics->{$name}->{cache}->{snps} = $cpt - $statistics->{$name}->{snp}->{nb};
	 	$statistics->{$name}->{cache}->{indels} = $cpt_indel - $statistics->{$name}->{indel}->{nb};
	 	 $statistics->{$name}->{cache}->{date} = 1;
	 	 $statistics->{$name}->{cache}->{date} = 0 if $delta <0;
	 	 $statistics->{$name}->{variations}->{cache} = $cpt;
	 } 
	
	 
	
}
