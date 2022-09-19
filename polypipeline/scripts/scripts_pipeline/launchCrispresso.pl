#!/usr/bin/perl
use Data::Dumper;
use File::Find;
use Getopt::Long;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/";
use lib "$Bin/../../../GenBo/lib/GenBoDB";
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../../../GenBo/lib/kyoto/";
use lib "$Bin/../../../GenBo/lib/GenBoDB/writeDB";
use lib "$Bin/../packages";
use GBuffer;

my $editing;
my $file;
my $projectName;
my $no_exec;
my $fork =20;

GetOptions(
	'project=s' => \$projectName,
	'file=s' => \$file,
	'editing=s' => \$editing,
	'no_exec=s' => \$no_exec,
	'fork=s' => \$fork,
 );





my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
my $dir = $project->getProjectRootPath();

open(SUMMARY, $file) or die ("$file fichier inconnu ");
my %hash;
while(<SUMMARY>){
	
	chomp($_);
	my($patient,$protocol,$size,$udi,$index1,$index2,$gRNA,$ngRNA,$amplicon,$scaffold,$extension)=split(",",$_);
	$hash{$patient}->{gRNA}=$gRNA;
	$hash{$patient}->{amplicon}=$amplicon;
	$hash{$patient}->{scaffold}=$scaffold;
	$hash{$patient}->{ngRNA}=$ngRNA;
	$hash{$patient}->{extension}=$extension;
	
}

warn Dumper \%hash;
close SUMMARY;

my $dirin;
my $script_file = $dir."/jobs.txt";
open(SCRIPT,">".$script_file) or die "impossible d'ouvrir le fichier de donnees $script_file";
foreach my $k (keys(%hash)){
	warn $k;
		my $pat = $project->getPatient($k);
		my $bc = $pat->barcode();
		$dirin = $pat->getSequencesDirectory();
		#$dirin = $dir;
		my $name = $pat->name();
		my $r1 = wanted($name,$bc,$dirin,"R1");
		#my $r1 = $bc."S_L003_R1_001.fastq.gz";
		warn $r1;
		my $r2 = wanted($name,$bc,$dirin,"R2");
		#my $r2 =  $bc."S_L003_R2_001.fastq.gz";
		my $run = $pat->getRun();
		my $type = $run->infosRun->{method};
		warn $type;
		my $enzyme = $pat->somatic_group();
		my $wnuc;
		my $enuc;
		if($enzyme eq "CBE"){
			$wnuc = "C";
			$enuc =  "T";
		}
		elsif($enzyme eq "ABE"){
			$wnuc = "A";
			$enuc = "G";
		}
		my $cmd = "docker run --rm -v $dir:/DATA -v $dirin:/SEQ -w /DATA  pinellolab/crispresso2 CRISPResso  -r1 /DATA/$r1 -a $hash{$name}->{amplicon} -g $hash{$name}->{gRNA} --no_rerun --output_folder /DATA/$name ";
		$cmd .= "-r2 /DATA/$r2" if $type eq "paired-end";
		$cmd .= " --quantification_window_size 10 --quantification_window_center -10 --base_editor_output --conversion_nuc_from $wnuc --conversion_nuc_to $enuc" if $editing==1;
		$cmd .= " --prime_editing_pegRNA_spacer_seq ".$hash{$name}->{gRNA};
		$cmd .= " --prime_editing_nicking_guide_seq ".$hash{$name}->{ngRNA};
		$cmd .= " --prime_editing_pegRNA_extension_seq ".$hash{$name}->{extension};
		
		print SCRIPT $cmd."\n\n";
		warn $cmd;
		 
}
close SCRIPT;
my $cmd2 = "cat $script_file \| run_cluster.pl -cpu=$fork";
system ($cmd2) unless $no_exec==1; 


#exec($script_file,$fork,$dirin,$dir) ;
#system ("mv $dirin/CRISPR* $dir");
#if ($step eq "tar" or $step eq "all"){
	
my $tar_cmd = "tar -cvzf $dir/$projectName.tar.gz $dir/CRISPR*";
die ("archive $dir/$projectName.tar.gz already exists") if -e $dir."/".$projectName.".tar.gz";
system ($tar_cmd)  unless $no_exec==1;
#or die "impossible $tar_cmd";
print "\t#########################################\n";
print "\t  link to send to the users : \n";
print "\t www.polyweb.fr/NGS/$projectName/$projectName.tar.gz \n";
print "\t#########################################\n";



#sub wanted {
#	my ($patient, $nb) =@_;
#	return unless -f;
#	if ($_ =~ /^$patient*bigfile1*.fastq.gz/) {return $_};
#}

sub wanted {
	my ($patient, $bc,$dirin,$pair) =@_;
	warn $dirin."/".$patient."*_".$pair."_*.fastq.gz";
	my @seq = glob($dirin."/".$patient."*_".$pair."_*.fastq.gz");
	my $list = join (" ",@seq);
	#my $fileout =  $project->getProjectRootPath()."/".$patient."_".$pair.".fastq.gz";
	my $fileout =$patient."_".$pair.".fastq.gz";
	my $path_fileout = $project->getProjectRootPath()."/".$fileout;
	my $cmd = "cat $list > $path_fileout";
	warn $cmd;
	system($cmd) unless $no_exec==1;
	chmod (0777,$fileout);
	return ($fileout);
}

sub exec {
	my ($script,$fork,$dirin,$dirout) = @_;
	system ("cat $script | run_cluster.pl -cpu=$fork");
	system ("mv $dirin/CRISPR* $dirout")
	
}



