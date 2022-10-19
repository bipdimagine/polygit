#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/";
use lib "$Bin/../../../GenBo/lib/GenBoDB";
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../../../GenBo/lib/kyoto/";
use lib "$Bin/../../../GenBo/lib/GenBoDB/writeDB";
use lib "$Bin/../../packages";
use Text::CSV qw( csv );
use Getopt::Long;
use GBuffer;
use Data::Dumper;
use List::MoreUtils qw(firstidx);
use colored; 
use IO::Prompt;
use Sys::Hostname;
use Term::ANSIColor;
use Moose;
use GBuffer;
use GenBoProject;
use colored; 
use Text::Table;
use file_util;
use File::Temp qw/ tempfile tempdir /;; 
use Term::Menus;
use Proc::Simple;
use Text::Table;
use dragen_util;
use XML::Simple qw(:strict);
use IO::Prompt;
my $project_names;
my $sample_sheet;
my $mask = undef;
my $l2;
my $force;
GetOptions(
	'project=s' => \$project_names,
	'l2=s' => \$l2,
);

my $bcl_dir;
my $aoa;
my $dir_out;
my %patients ;
my $dir_fastq;
my $run_name;

foreach my $project_name (split(",",$project_names)){
	my $buffer = GBuffer->new();
	my $project = $buffer->newProject( -name 			=> $project_name );
	
	foreach my $capture (@{$project->getCaptures}){
		
		 $mask = $capture->umi->{mask} unless defined $mask;
		 die("problem more than one mask ") if $mask ne $capture->umi->{mask};
		 
	}
	map {$patients{$_->name}++} @{$project->getPatients};
	my $run = $project->getRun();
	
	$run_name = $run->run_name;
	
	die("can't find bcl directory : ".$run->bcl_dir()) unless -e $run->bcl_dir();

	if ($bcl_dir && $run->bcl_dir() ne $bcl_dir){
		die("problem different bcl dir : $bcl_dir ".$run->bcl_dir);
	}
	$bcl_dir = $run->bcl_dir;
	next if $aoa;
	my $csv_tmp = $bcl_dir."/SP.".time.".csv";
		die("no sample sheet in database ") unless ($run->sample_sheet);
	my $toto = $run->sample_sheet;
	#$run->sample_sheet =~ s/;/,/g;	
	$toto  =~ s/;/,/g;
	open(TOTO,">$csv_tmp");
	print TOTO $toto;
	close TOTO;
	$aoa = csv (in => $csv_tmp); 
	$dir_out = $project->project_dragen_demultiplex_path();
	$dir_fastq = $project->dragen_fastq;
}
if ($mask){
	warn $bcl_dir;
	my $config = XMLin("$bcl_dir/RunInfo.xml", KeyAttr => { reads => 'Reads' }, ForceArray => [ 'reads', 'read']);
	my $read = $config->{Run}->{Reads}->{Read};
	my $y = "Y".$read->[0]->{NumCycles};
	 $mask = $y.";".$mask.";".$y; 
	
}


my $choice = prompt("use \"$mask\" for demultipexing  (y/n) ? ");
if ($choice ne "y") {
	$mask =  prompt("enter your mask  ? ");
	die($mask);
}

#die() unless -e $dir_in;




my $lines;
my $titles;
my $current_title;


foreach my $line(@$aoa){
	
		if ($line->[0] =~ /^\[/) {
			$current_title = $line->[0];
			push(@$titles,$current_title);
		}
		else {
			push(@{$lines->{$current_title}},$line);
		}
	
	
}

#change setting
foreach my $set (@{$lines->{"[Settings]"}}){
	if ($set->[0] eq "Adapter") {
		$set->[0] = "AdapterRead1"
	}
	
}

### read mask ;
my @amask = split(";",$mask);
my $pos_umi = firstidx { $_ =~ /U/ } @amask;

$lines->{"[Settings]"} = [];
push(@{$lines->{"[Settings]"}},["BarcodeMismatchesIndex1",0]);
push(@{$lines->{"[Settings]"}},["BarcodeMismatchesIndex2",0]) unless $pos_umi;
push(@{$lines->{"[Settings]"}},["OverrideCycles",$mask]) if $mask;

my $lheader_data = shift @{$lines->{"[Data]"}};
if (scalar (@$lheader_data) ne  scalar (@{$lines->{"[Data]"}->[0]})){
	my $tb = Text::Table->new(
	{is_sep => 1},
	@$lheader_data,
    );
    $tb->load($lines->{"[Data]"}->[0]);
    print $tb;
	print "\n";
	die(scalar (@$lheader_data)."vs" .scalar (@{$lines->{"[Data]"}->[0]}));
}

my $error_not_in_project = {};
my $ok;
my $pos_sample = firstidx { $_ eq "Sample_ID" } @$lheader_data;
die("no sample id in header ") if $pos_sample eq -1;
my $pos_sample_name = firstidx { $_ eq "Sample_Name" } @$lheader_data;
my $pos_cb1 = firstidx { $_ eq "index" } @$lheader_data;
my $pos_cb2 = firstidx { $_ eq "index2" } @$lheader_data;

### read mask ;
#my @amask = split(";",$mask);
#my $pos_umi = firstidx { $_ =~ /U/ } @amask;
if($pos_umi == 2){
	splice(@$lheader_data, $pos_cb2, 1);
	foreach my $data (@{$lines->{"[Data]"}}){
	if($pos_umi == 2){
		splice(@$data, $pos_cb2, 1);
		#my $pos_cb2 = firstidx { $_ eq "index2" } @$lheader_data;
	}
	}
}

my $dj;
foreach my $data (@{$lines->{"[Data]"}}){
	
	if($l2){
 		my $index2 =  $data->[$pos_cb2];
 		$data->[$pos_cb2] = substr($index2, 0, $l2);
 	}
	my $name = $data->[$pos_sample];
	$data->[$pos_sample_name] = $name;
	next if $name =~ /_RC/; 
	next if exists $dj->{$name};
	next unless $name;
	$error_not_in_project->{$name} ++ unless exists $patients{$name};
	$ok->{$name} ++;
 	delete $patients{$name};
 	$dj->{$name}++;
 	
	
}

my $error;
if (keys %$ok) {
	print colored(['bright_green on_black']," SAMPLE OK : ".scalar(keys %$ok) )."\n";
}
if (keys %patients) {
	print colored(['bright_red on_black']," SAMPLE IN PROJECT NOT IN SAMPLE SHEET :")."\n";
	map {print $_."\n"} keys %patients;
	$error =1;
}

if (keys %$error_not_in_project) {
		print colored(['bright_red on_black']," SAMPLE IN SAMPLE SHEET  NOT IN PROJECT :")."\n";
		
		map {print $_."\n"} keys %$error_not_in_project;
		print " nb error : ".scalar (keys %$error_not_in_project)."\n";
		my $choice = prompt("continue anyway  (y/n) ? ");
		die() if ($choice ne "y"); 
		#$lines->{"[Data]"} = change_sample_sheet($pos_sample,$lines->{"[Data]"},$error_not_in_project);
		
}

die() if $error;
unshift( @{$lines->{"[Data]"}},$lheader_data);
my $outcsv;

foreach my $title (@{$titles}){
	push(@$outcsv,[$title]);
	push(@$outcsv,@{$lines->{$title}});
}
my $ss = $bcl_dir."/file".time.".csv";
csv (in => $outcsv, out => $ss, sep_char=> ",");

sleep(1);



my $cmd = qq{dragen --bcl-conversion-only=true --bcl-input-directory $bcl_dir --output-directory $dir_out --sample-sheet $ss --force  };
my $exit = 0;
warn qq{$Bin/../run_dragen.pl -cmd="$cmd"};
$exit = system(qq{$Bin/../run_dragen.pl -cmd="$cmd"});


die() if $exit ne 0;

warn "END DEMULTIPEX \n let's copy ";
my $fork =6;
my $pm   = new Parallel::ForkManager($fork);

foreach my $project_name (split(",",$project_names)){
	my $buffer = GBuffer->new();
	my $project = $buffer->newProject( -name 			=> $project_name );
	my $run = $project->getRun;
	my $out_fastq = $run->fastq_dir();
	system("mkdir $out_fastq ; chmod g+rwx $out_fastq ");
	
	foreach my $p (@{$project->getPatients}){
		my $pid = $pm->start and next;
		my ($fastq1,$fastq2) = dragen_util::get_fastq_file($p,$out_fastq,$dir_out);
		
	#	create_3_fastq($fastq1,$fastq2,$p);
	#	warn "end ".$p->name;
	#system ("rsync -rav $dir_out/".$p->name."_S* $out_fastq/");
		$pm->finish( 0, {});
	}

}

$pm->wait_all_children();

my $dir_stats = "/data-isilon/sequencing/ngs/demultiplex/".$run_name;


system("mkdir -p $dir_stats && rsync -rav --remove-source-files ".$dir_out."/Reports/ $dir_stats/ && rm $dir_out/Reports/* && rmdir $dir_out/Reports/ && chmod -R a+rwx $dir_stats  ");

exit(0);
###

sub create_3_fastq {
	my ($fastq1,$fastq2,$patient) = @_;
	my $fastq3_prod  = $fastq2;
	$fastq3_prod =~ s/_R2/_R3/;
	#system("mv $fastq2 $fastq3");
	open(FASTQR,"zcat $fastq2 | ");
	
	my $fastq3 = "/data-beegfs/tmp/".$patient->name."_R3.fastq";
	open(FASTQ3,"> $fastq3 ");
	my $fastq_umi ="/data-beegfs/tmp/".$patient->name."_R2.fastq";
	open(FASTQU,"> $fastq_umi ");
	my $line;
	while($line = <FASTQR>){
		chomp($line);
		if ($line =~ /^@/) {
			my $line3 = $line;
			$line3 =~ s/2:N:0:/3:N:0:/;
			print FASTQ3 $line3."\n";
			my @t = split(" ",$line);
			my @us = split(":",$t[0]);
			
			print FASTQU "$line\n".$us[-1]."\n+\nFFFFFFFFFF\n";
				
			}
		else {
			print FASTQ3 $line."\n";
		}
	}
	close(FASTQ3);
	close(FASTQU);
	system("gzip $fastq3");
	system("gzip $fastq_umi");
	warn "mv $fastq3.gz $fastq3_prod";
	system("mv $fastq3.gz $fastq3_prod");
	system("mv $fastq_umi.gz $fastq2");
}

sub change_sample_sheet {
	my ($pos,$array,$patients) = @_;
	my $new_array;
	foreach my $l (@$array){
		my $name = $l->[$pos];
		if (exists $patients->{$name}){
			next;
		}
		push(@$new_array,$l);
	}
	return $new_array;
	
}


