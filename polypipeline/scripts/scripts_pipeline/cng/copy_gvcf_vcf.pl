#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
#use Set::IntSpan;
use GenBoNoSqlLmdb;

use Carp;
use strict;
use Set::IntSpan::Fast::XS ;
use Data::Dumper;
use GBuffer;
use Proc::Simple;
use Getopt::Long;
use File::Util;
my $myproc = Proc::Simple->new(); 

 $myproc->redirect_output ("/dev/null", "/dev/null");
my $fork;
my $project_name;
my $debug;
my $type;
my $set;
GetOptions(
	'project=s' => \$project_name,
	'debug=s' => \$debug,
	'name=s' => \$type,
	'set=s' => \$set,
);

my @list = `cat /software/polyweb/poly-disk/poly-src/defidiag/project/download/list.txt`;
chomp(@list);
my ($l) = grep{$_ =~ /$type/} @list;
die() unless $l;

my($name,$password,$type2) = split(" ",$l);

my $buffer = new GBuffer;
my $project = $buffer->newProject( -name 			=> $project_name );

my $dir_in = "/data-isilon/download/$type2/www.cnrgh.fr/data/".$type;
$dir_in = "/data-isilon/download/$type2/www.cnrgh.fr/data/".$type if lc($type) =~/devo/;
my $bcftools = $buffer->software("bcftools");
foreach my $patient (@{$project->getPatients}){
	my $bc = $patient->barcode();
	my @files = `find $dir_in -name *$bc*`;
	warn "find $dir_in -name *$bc*";
	chomp(@files);
	
	my ($tar) = grep{$_=~ /tar.gz$/} @files;
	warn "\t".$tar;
	next if -e $tar;
		
	 ($tar) = grep{$_=~ /tar$/} @files;
	next if -e $tar;
	
	my $dir = "/data-isilon/download/$type2/";
	my $cmd_var="wget  -N  --recursive --no-parent --user $name  --no-check-certificate --password  $password https://www.cnrgh.fr/data/$name/variants/$set/ -A \"*$bc*\" -R cram";  
	system("cd $dir && $cmd_var");
	 @files = `find $dir_in -name *$bc*`;
	chomp(@files);
	 ($tar) = grep{$_=~ /tar.gz$/} @files;
	die($tar) unless -e $tar;
	#die() unless -e 
}
warn "untar";
untar();
warn "bcftools";
bcftools();
exit(0);


foreach my $patient (@{$project->getPatients}){
	my $bc = $patient->barcode();
	warn $bc." ".$patient->name;
	next if -e $patient->gvcfFileName("haplotypecaller4");
	my $files;
	my @files = `find $dir_in -name *$bc*`;
	chomp(@files);
	my ($tar) = grep{$_=~ /tar.gz$/} @files;
	my ($tar2) = grep{$_=~ /tar$/} @files;
	my $cgzip="";
	if ($tar) {
		 $tar2 = $tar;
		$tar2 =~ s/\.gz//;
		$cgzip = "gunzip $tar";
	
	}
	else {
		die();
		
	#	die($bc." ==> ".$patient->name) unless $tar;
	}
	
	
	
	my ($manta) = grep{$_=~ /diploidSV.vcf.gz$/} @files;
	warn "cd $dir_in;$cgzip;tar -xvf ".$tar2 unless -e $manta;;
	#die();
	system("cd $dir_in;$cgzip;tar -xvf ".$tar2) unless -e $manta;
	warn "end ";
	@files = `find $dir_in -name *$bc*`;
	chomp(@files);
	 ($manta) = grep{$_=~ /diploidSV.vcf.gz$/} @files;
	die(Dumper @files." ".$manta) unless $manta;
	my $ff = $dir_in."/".$bc.".txt";
	system("echo ".$patient->name.">".$ff);
	reheader($manta,$patient->getVariationsFileName("manta"),$ff);
	my ($canvas) = grep{$_=~ /canvas_CNV.vcf.gz$/} @files;
	reheader($canvas,$patient->getVariationsFileName("canvas"),$ff);
	die("--") unless $canvas;
	my ($gvcf) = grep{$_=~ /g.vcf.gz$/} @files;
	die("gvcf ".$bc) unless $gvcf;
	reheader($gvcf,$patient->gvcfFileName("haplotypecaller4"),$ff);
	
	
#	warn "tar -xvf $dir_in/".$gz[0];
#	system("cd $dir_in;tar -xvf ".$gz[0]);
#	my $dir1 = $gz[0];
#	die() unless -e $dir_in."/".$gz[0];
#	$dir1 =~ s/\.tar\.gz//;
#	$dir1 = $dir_in."/".$dir1;
#
#	
#	my @gvcf = $f->list_dir($dir1=> { 
#			no_fsdots =>1,
#			files_only => 1,
#			files_match => {or => [qr/$bc/],and =>[".g.vcf.gz\$"]}
#	});
#	die() if scalar (@gvcf) != 1;
#	$files->{gvcf} = $dir1."/".$gvcf[0];
#	
#	my $out = $patient->gvcfFileName("haplotypecaller4");
#	system("echo ".$patient->name.">".$dir1."/p.txt");
#	my $cmd = "$bcftools reheader -s ".$dir1."/p.txt ".$files->{gvcf}.">".$out ." && tabix -f -p vcf ".$out;
#	system("$cmd");
#	
#	warn $cmd; 

}

sub bcftools {
	my @cmd;
	warn 'bcftools -----------';
	foreach my $patient (@{$project->getPatients}){
		my $bc = $patient->barcode();
		warn $bc." ".$patient->name;
		#next if -e $patient->gvcfFileName("haplotypecaller4");
		warn $bc." ".$patient->name;
		warn "find $dir_in -name '*$bc*'";
		my @files = `find $dir_in -name '*$bc*'`;
		chomp(@files);
		
	my ($manta) = grep{$_=~ /diploidSV.vcf.gz$/} @files;
	die(Dumper @files." ----- ".$manta) unless $manta;
	my $ff = $dir_in."/".$bc.".txt";
	#push(@cmd,"echo ".$patient->name.">".$ff);
	system("echo ".$patient->name.">".$ff);
	push(@cmd,reheader($manta,$patient->getVariationsFileName("manta"),$ff) ) unless -e $patient->getVariationsFileName("manta").".tbi";
	my ($canvas) = grep{$_=~ /canvas_CNV.vcf.gz$/} @files;
	push(@cmd , reheader($canvas,$patient->getVariationsFileName("canvas"),$ff)) unless -e $patient->getVariationsFileName("canvas").".tbi" ;
	die("--") unless $canvas;
	my ($gvcf) = grep{$_=~ /g.vcf.gz$/} @files;
	die("gvcf ".$bc) unless $gvcf;
	push(@cmd,reheader($gvcf,$patient->gvcfFileName("haplotypecaller4"),$ff)) unless -e $patient->gvcfFileName("haplotypecaller4").".tbi" ;
	
	my ($vcf_gatk4) = grep{$_=~ /gatk4HC_annotated.vcf.gz$/} @files;
	push(@cmd,reheader2($vcf_gatk4,$patient->getVariationsFileName("haplotypecaller4"),$ff)) unless -e $patient->getVariationsFileName("haplotypecaller4").".tbi" ;
	warn  $patient->getVariationsFileName("haplotypecaller4").".tbi";
	
	}
	warn Dumper @cmd;
	open(CMD,'>cmd.test.txt');
	print CMD join("\n",@cmd)."\n";
	close CMD;
	system("cat cmd.test.txt | parallel -j 8 "); 
}
sub untar {
	foreach my $patient (@{$project->getPatients}){
	my $bc = $patient->barcode();
	warn $bc." ".$patient->name;
	warn $patient->gvcfFileName("haplotypecaller4");
	next if -e $patient->gvcfFileName("haplotypecaller4");
	my $files;
	warn qq{find $dir_in -name \'*$bc*\'};
	my @files = `find $dir_in -name \'*$bc*\'`;
	chomp(@files);
	my ($targz) = grep{$_=~ /tar.gz$/} @files;
	#my ($tar2) = grep{$_=~ /.tar$/} @files;
	my $cgzip="";
	if ($targz) {
		system("pigz -d -p 5  $targz");
		
	}
	@files = `find $dir_in -name \'*$bc*\'`;
	chomp(@files);
	my ($tar2) = grep{$_=~ /.tar$/} @files;
	die($bc) unless -e $tar2;
	
	my ($manta) = grep{$_=~ /diploidSV.vcf.gz$/} @files;
	warn "cd $dir_in;tar -xvf ".$tar2 unless -e $manta;;
	#die();
	system("cd $dir_in;tar -xvf ".$tar2) unless -e $manta;
	warn 'next';
	}
	warn "end ";
}
sub reheader2 {
	my($f1,$f2,$ff) = @_;
	my $cmd = "$bcftools annotate -x INFO/ANN $f1 | $bcftools reheader - -s $ff | bgzip -c >".$f2." && tabix -f -p vcf ".$f2;
	return $cmd;
	warn $cmd;
	system($cmd);
}
sub reheader {
	my($f1,$f2,$ff) = @_;
	my $cmd = "$bcftools reheader -s $ff $f1>".$f2 ." && tabix -f -p vcf ".$f2;
	return $cmd;
	warn $cmd;
	system($cmd);
}