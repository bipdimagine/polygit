#!/usr/bin/perl
use CGI qw/:standard :html3/;
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
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
use JSON;
use check_utils;

my $projectName;
my $end_ext = "uni";
my $details;
my $fork = 1;
my $vcf_file;
my $log_file;

GetOptions(
	'project=s'		=> \$projectName,
	"details=s"		=> \$details,
	"fork=s"		=> \$fork,
	"vcf_file=s" => \$vcf_file,
	 "log=s" => \$log_file,
);

my $date = `date`;

chomp($date);
if ($log_file){
	open (STDOUT,">>".$log_file);
}

print "echo  ---- start plink $date\n";
die("\n\nERROR: -project option mandatory. Exit...\n\n") unless ($projectName);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
my $ped_file = $project->getPedigreeFile();
die("\n\nERROR: no pedigree file found. Exit...\n\n") unless (-e $ped_file);
my @lChr = ('21', 'X', 'Y');
my $children;
my $family;
eval { $family = $project->families(); };
if ($@) { die("\n\nERROR: can't construct family relations... Maybe a problem in PED file. Die...\n\n"); }
checkPedFile();
checkPatientsName();
unless (-e $vcf_file) {
callingChr21XYAllPatients();
$vcf_file = concatVcf();
}

launchPlink($vcf_file);
checkPlink();

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
 
sub checkPatientsName {
	my $hPatNames;
	foreach my $pat (@{$project->getPatients()}) { $hPatNames->{$pat->name()} = 'false'; }
	my $ped_file;
	unless ($ped_file) { $ped_file = $project->getPedigreeFile(); }
	die() unless -e $ped_file;
	open (PED, $ped_file);
	while (<PED>) {
		chomp($_);
		my $line = $_;
		my @lFields = split("\t", $line);
		my $patName = $lFields[1];
		unless (exists ($hPatNames->{$patName})) { warn "======= WARN: $patName in ped file not declared in this project !"; }
		$hPatNames->{$patName} = 'true';
	}
	close(PED);
	foreach my $patName (keys(%$hPatNames)) {
		if ($hPatNames->{$patName} eq 'false') { warn "======= WARN: $patName is declared in this project but not included in ped file !"; }
	}
}

sub callingChr21XYAllPatients {
	print "\n\n### Step 1: calling chr 21, X and Y for each patient.\n";
	my $ext = "uni";
	my $dir_out = $project->getCallingPipelineDir("unifiedgenotyper").'/plink';
	unless (-d $dir_out) { `mkdir $dir_out`; }
	my $log_file = $dir_out."/callingChr_21_X_Y.log";
	foreach my $chrName (@lChr) {
		my $cmd = "perl $Bin/fast_unifiedGenotyper.pl -patient=all -project=$projectName -keep=1 -chr=$chrName -remove_par_regions=true -fork=$fork -ext=$ext -dir_out=$dir_out -log=$log_file";
		`$cmd`;
	}
}

sub concatVcf {
	print "\n\n### Step 2: concat VCF files.\n";
	my $dir_out = $project->getCallingPipelineDir("unifiedgenotyper").'/plink';
	my @lVcfFiles;
	foreach my $chrName (@lChr) { push(@lVcfFiles, $dir_out.'/'.$chrName.'.uni.vcf') }
	my $fileOut = $dir_out."/".$projectName.".uni.vcf";
	my $vcftools = "/bip-d/soft/distrib/vcftools/latest/bin/";
	my $cmd = "$vcftools/vcf-concat ".join(" ", @lVcfFiles)." >".$fileOut;
	warn $cmd;
	`$cmd`;
	return $fileOut;
}

sub launchPlink {
	my ($vcf_file) = @_;
	
	print "\n\n### Step 3: launching PLINK.\n";
	my $vcftools = "/bip-d/soft/distrib/vcftools/latest/bin/vcftools";
	my $dir = $project->getCallingPipelineDir("unifiedgenotyper");
	
	my $logPlink = $dir.'/plink/'.$projectName.'.plink.resume';
	
	my @snps;

	warn $vcf_file;
	open (VCF, "<$vcf_file");
	my @samples;
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
		for (my $i=0;$i<@samples;$i++){
			my $name = $samples[$i];
			my $string = $gsamples[$i];
			my $geno;
			if ($string =~ /0\/0/) { $geno ="$ref $ref"; }
			elsif ($string =~ /0\/1/) { $geno ="$ref $alt"; }
			elsif ($string =~ /1\/0/) { $geno ="$ref $alt"; }
			elsif ($string =~ /1\/1/) { $geno ="$alt $alt"; }
			elsif ($string =~ /.\/./) { $geno ="0 0"; }
			else { die($string); }
			$snp->{samples}->{$name} = $geno;
		}
		push(@snps,$snp);
	}
	
	$dir .= "/plink/";
	mkdir $dir unless -e $dir;
	unless ($ped_file) { $ped_file = $project->getPedigreeFile(); }
	die() unless -e $ped_file;
	my @samples_name;
	open(PED,$ped_file);
	while (my $line = <PED>){
		chomp($line);
		my @t = split(" ",$line);
		push(@samples_name,$t[1]);
	}
	close PED;
	
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
	my $cmd2 = "/bip-d/soft/distrib/plink/plink-1.07-x86_64/plink --tped $tped_file --tfam $ped_file --noweb --check-sex --out $dir/$projectName";
	warn $cmd2;
	`$cmd2`;
	my $cmd3 = "/bip-d/soft/distrib/plink/plink-1.07-x86_64/plink --tped $tped_file --tfam $ped_file --mendel --noweb --map3 --check-sex --out $dir/$projectName 1>$logPlink 2>>$logPlink";
	warn $cmd3;
	`$cmd3`;
	my $file1 = "$dir/".$project->name();
	my @titi = `cat  $file1.imendel`;
	foreach my $line (@titi){
		chomp($line);
		next if $line =~ /FID/;
		my ($f, $name,$nb) = split(" ",$line);
		my $p = int($nb/scalar(@snps)*1000)/10;
		print $line."\t".$p."\n";
		$family->{$f}->{$name}->{'mendelian_errors'} = $p.'%';
		if (exists($children->{$name})) { $children->{$name}->{'mendelian_errors'} = $p.'%'; }
	}
}

sub checkPlink {
	print "\n\n### Step 4: checking PLINK.\n";
	my $pedFile = $project->getPedigreeFile();
	die("No PED file found... Die.") unless ($pedFile);
	my $hashSex = getHashVerifSex();
	my $listHash = createListHashes($hashSex);
	my $jsonFile = $project->getProjectPath().'../'.$projectName.'.resume.json';
	my $jsonRes = check_utils::createJson($listHash, $jsonFile);
	checkSex($listHash);
	checkMendelianErrors();
	open(JSON_FILE, ">$jsonFile");
	print JSON_FILE $jsonRes;
	close(JSON_FILE);
	my ($nbWarn, $nbErr) = check_utils::getTableFromJsonFile($jsonFile);
}

sub getHashVerifSex {
	my $hashSex;
	my $sexFile = $project->getCallingPipelineDir("unifiedgenotyper").'/plink/'.$projectName.'.sexcheck';
	unless (-e $sexFile) { die("\n\nERROR: PLINK sex file doesn't exist... Die...\n\n") }
	open(SEX, "<$sexFile");
	while(<SEX>) {
		chomp($_);
		my $line = $_;
		next if($line =~ /PEDSEX/);
		$line =~ s/\s+/ /g;
		my @lFields = split(" ", $line);
		$hashSex->{$lFields[0]}->{$lFields[1]}->{'ped_sex'} = $lFields[2];
		$hashSex->{$lFields[0]}->{$lFields[1]}->{'snp_sex'} = $lFields[3];
		$hashSex->{$lFields[0]}->{$lFields[1]}->{'status_sex'} = $lFields[4];
		$hashSex->{$lFields[0]}->{$lFields[1]}->{'f_sex'} = $lFields[5];
	}
	close(SEX);
	return $hashSex;
}

sub createListHashes {
	my ($hashSex) = @_;
	my $thisFamily;
	foreach my $famName (%{$project->families()}) {
		foreach my $type (%{$project->families()->{$famName}}) {
			if (exists $project->families()->{$famName}->{'father'}) {
				my $fatherName = $project->families()->{$famName}->{'father'};
				$thisFamily->{$famName}->{$fatherName}->{'relation'} = 'father';
					$thisFamily->{$famName}->{$fatherName}->{'relation_code'} = '0';
			}
			if (exists $project->families()->{$famName}->{'mother'}) {
				my $motherName = $project->families()->{$famName}->{'mother'};
				$thisFamily->{$famName}->{$motherName}->{'relation'} = 'mother';
				$thisFamily->{$famName}->{$motherName}->{'relation_code'} = '1';
			}
			if (exists $project->families()->{$famName}->{'child'}) {
				foreach my $childName (@{$project->families()->{$famName}->{'child'}}) {
					$thisFamily->{$famName}->{$childName}->{'relation'} = 'child';
						$thisFamily->{$famName}->{$childName}->{'relation_code'} = '2';
				}
			}
		}
	}
	my @lHashFamily;
	foreach my $famName (sort(keys(%$thisFamily))) {
		foreach my $patName (keys(%{$thisFamily->{$famName}})) {
			my $hash;
			$hash->{'family'} = $famName;
			$hash->{'patient'} = $patName;
			$hash->{'relation'} = $thisFamily->{$famName}->{$patName}->{'relation'};
			$hash->{'relation_code'} = $thisFamily->{$famName}->{$patName}->{'relation_code'};
			if (exists($family->{$famName}->{$patName}->{'mendelian_errors'})) { $hash->{'plink_mendelian_errors'} = $family->{$famName}->{$patName}->{'mendelian_errors'}; }
			else { $hash->{'plink_mendelian_errors'} = 'PLINK_not_used'; }
			$hash->{'ped_sex'} = $hashSex->{$famName}->{$patName}->{'ped_sex'};
			$hash->{'plink_f_sex'} = $hashSex->{$famName}->{$patName}->{'f_sex'};
			$hash->{'plink_snp_sex'} = $hashSex->{$famName}->{$patName}->{'snp_sex'};
			if (($hash->{'plink_snp_sex'} eq '0') and ($hash->{'plink_f_sex'} ne 'nan')) { $hash->{'plink_status_sex'} = "Not Sure"; }
			else { $hash->{'plink_status_sex'}  = $hashSex->{$famName}->{$patName}->{'status_sex'}; }
			push(@lHashFamily, $hash);
		}
	}
	return \@lHashFamily;
}

sub checkSex {
	my ($listHash) = shift;
	if ($details) { warn Dumper $listHash; }
	my $hasProblem = 0;
	foreach my $hash (@$listHash) {
		if ($hash->{'plink_status_sex'} eq 'PROBLEM') { $hasProblem++; }
		if ($hash->{'plink_status_sex'} eq 'Not Sure') { $hasProblem++; }
	}
	print "\n\n";
	if ($hasProblem > 0) {
		print  colored ['black ON_BRIGHT_RED'],"Check Sex: ERROR !!!"; 
		print  color 'reset';
		print  "\n";
	}
	else {
		print  colored ['black ON_BRIGHT_GREEN'],"Check Sex: OK"; 
		print  color 'reset';
		print  "\n";
	}
}

sub checkMendelianErrors {
	my $hasProblem = 0;
	my @patErrors;
	foreach my $childName (keys(%$children)) {
		my $perc = $children->{$childName}->{'plink_mendelian_errors'};
		$perc =~ s/%//;
		if (int($perc) >= 10) {
			my $familyName = $children->{$childName}->{'family'};
			push(@patErrors, @{$family->{$familyName}->{'child'}});
			push(@patErrors, @{$family->{$familyName}->{'parents'}});
		}
	}
	if (scalar(@patErrors) > 0) { foreach my $patName (@patErrors) { $hasProblem++; }	}
	if ($hasProblem > 0) {
		print  colored ['black ON_BRIGHT_RED'],"Check Mendelian Errors: ERROR !!!"; 
		print  color 'reset';
		print  "\n\n";
	}
	else {
		print  colored ['black ON_BRIGHT_GREEN'],"Check Mendelian Errors: OK"; 
		print  color 'reset';
		print  "\n\n";
	}
}


