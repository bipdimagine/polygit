#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../../GenBo/lib/obj-nodb/packages";
use GBuffer;
use GenBoProject;
use Getopt::Long;
use Data::Dumper;
#use illumina_util;
use IO::Prompt;
use Sys::Hostname;
 use File::Find::Rule ;
use Text::Table;
use Term::Twiddle;
use colored; 
my $project_name_trio;
my $project_name_solo;
my $project_name;
my $chr_name;
my $ppn;
my $set;
my $name;
my $solo;
GetOptions(
	'project=s' => \$project_name,
	'fork=s' => \$ppn,
	'set=s' => \$set,
	'name=s' => \$name,
	'solo=s' => \$solo,
);
die() unless $solo;
$project_name_solo = $solo;
my @list;
my $dir2 = "";

$set = "set".$set unless $set=~/set/;
if($name){
 @list = `cat ../../../../defidiag/project/$name/$set.txt`;
 chomp(@list);
}
elsif($project_name){
	push(@list,$project_name);
}
die() unless @list;





foreach my $project_name_trio (@list) {
my $buffer_solo = GBuffer->new();
my $project_solo = $buffer_solo->newProject( -name => $project_name_solo );
my $patients_solo = $project_solo->getPatients();
warn "patient solo =>".scalar(@$patients_solo);
my $buffer_trio = GBuffer->new();
warn $project_name_trio;
my $project_trio = $buffer_trio->newProject( -name => $project_name_trio );
my $patients_trio = $project_trio->getPatients();
warn scalar(@$patients_trio);

my @common ;
my $find;
foreach my $pt (@$patients_trio){
		next unless $pt->isChild();
	foreach my $ps (@$patients_solo){
		
		warn $pt->name ." ==>  ". $ps->name() if $pt->name eq $ps->name();
		push(@common,$pt->name) if $pt->name eq $ps->name();
		$find ++ if $pt->name eq $ps->name();
	}
} 

#die("not found ") unless $find ==1;
my @allcmd;

foreach my $pname (@common){
	push(@allcmd,move_patient($project_trio,$project_solo,$project_trio->getPatient($pname)));
}





my $primer_trio  = $project_trio->noSqlCnvsDir."/primers.bed.gz";
my $primer_solo  = $project_solo->noSqlCnvsDir."/primers.bed.gz";

#unless (-e $primer_solo){
#		my $cmd = "rsync -rav".$project_trio->noSqlCnvsDir."/ ".$project_solo->noSqlCnvsDir."/";
#		push(@allcmd,$cmd);
#}
foreach my $cmd (@allcmd){
	system($cmd);
	
}

}
system("$Bin/add_defidiag_phentoype.pl -project=$solo");
warn "ADD PHENOTYPE AND CALLING METHODS !!!!!!";
exit(0);





sub move_patient{
	my ($project,$newProject,$patient) =@_;
	my $project_in = $project->name;
	warn $project_in;
	my $projectName = $newProject->name;	
	my $pn = $patient->name();
 	my $dir1 = $project->getProjectPath();
 	my $dir2 = $newProject->getProjectPath();
 	my $dir3 = $project->getSequencesRootDir();

 #my @files = `find $dir1 -name $pn.*`;
  my @files = File::Find::Rule->file()
                                  ->name( "$pn.*" )
                                  ->in( $dir1);
	warn Dumper @files;
	warn $dir1;                                  
  my (@bam) = grep {/bam$/} @files;
   my (@vcf) = grep {/vcf.gz$/} @files;
   warn scalar @vcf;
  
 if (scalar(@bam) ne 1 || scalar(@vcf) < 2){
 	colored::stabilo("red","----------------$pn problem with files---------------");
	colored::stabilo("red","----------------".join(",",@bam)." bam file  ---------------");
	colored::stabilo("red","---------------- ".join(",",@vcf)." vcf  file  ---------------");
	colored::stabilo("red","-----------------------------------------------------"); 
	warn Dumper @files;
 	
 }
 
my @files2 = File::Find::Rule->file()
                                  ->name( "$pn"."_*" )
                                  ->in( $dir1);

warn Dumper @files2;
push(@files,@files2);
my $bam;
my $vcf;
my @cmd;
foreach my $f (@files){
	my $f2 = $f;
	$bam ++ if $f2 =~ /bam$/;
	$vcf ++ if $f2 =~ /vcf.gz$/;
	$f2 =~   s/$project_in/$projectName/;
	warn ("file exists $f2") if -e $f2;
	next if -e $f2;
	
	die ("file exists $f2 $f") if -e $f2;
	my @dir = split("/",$f2);
	pop(@dir);
	
	my $dir_out = join("/",@dir);
	system("mkdir -p $dir_out") unless -e $dir_out;
	
	die("no directory  $dir_out") unless -e $dir_out;
	push(@cmd,"ln -s $f $f2");
	
} 

if ($bam ne 1 || $vcf < 2){
	colored::stabilo("red","----------------$pn problem with files---------------");
	colored::stabilo("red","----------------$bam bam file  ---------------") if $bam ne 1; 
	colored::stabilo("red","---------------- $vcf vcf  file  ---------------") if $vcf < 2; 
	colored::stabilo("red","-----------------------------------------------------"); 
	warn Dumper @cmd;
	die();
}
	
	
 return @cmd;
}
exit(0);