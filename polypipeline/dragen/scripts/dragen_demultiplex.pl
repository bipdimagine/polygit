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
use Moo;
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
use Statistics::Descriptive;

my $project_names;
my $sample_sheet;
my $mask = undef;
my $l2;
my $force;
my $run_name_option;
my $mismatch; 

GetOptions(
	'project=s' => \$project_names,
	'l2=s' => \$l2,
	'run=s' => \$run_name_option,
	'mismatch=s' => \$mismatch,
);
 
my $bcl_dir;
my $aoa;
my $dir_out;
my %patients ;
my $dir_fastq;
my $run_name;
my $umi_name;
my $dir_bcl_tmp;
foreach my $project_name (split(",",$project_names)){
	my $buffer = GBuffer->new();
	my $project = $buffer->newProject( -name 			=> $project_name );
	my $runs = $project->getRuns;
	my $run;
	if (scalar(@$runs)> 1){
		unless ($run_name_option){
			print colored(['bright_red on_black'],"Hey Manue you have ".scalar(@$runs)." runs in the projects $project_name you have to choose one and add -run on the command line" )."\n";
			foreach my $run (@$runs){
				print colored(['bright_green on_black'],$run->plateform_run_name)." ".colored(['bright_blue on_black'],$run->date)."\n";
			}
			die();
		}
		($run) = grep{$_->plateform_run_name eq "$run_name_option"} @$runs;
		die("unable to find $run_name_option ".$project_name) unless $run;
	}
	else {
		$run = $runs->[0];
	}
	warn $run->name();
	#die() if scalar(@$runs)> 1;
	#warn 
	foreach my $capture (@{$project->getCaptures}){
		if ($capture->umi){
		 $mask = $capture->umi->{mask} unless defined $mask;
		 $umi_name = $capture->umi->{name};
		 die("problem more than one mask ") if $mask ne $capture->umi->{mask};
		}
	}
	map {$patients{$_->name}= $project->name} @{$project->getPatients};
	#my $run = $project->getRun();
	
	$run_name = $run->run_name;
	
	die("can't find bcl directory : ".$run->bcl_dir()) unless -e $run->bcl_dir();

	if ($bcl_dir && $run->bcl_dir() ne $bcl_dir){
		die("problem different bcl dir : $bcl_dir ".$run->bcl_dir);
	}
	$bcl_dir = $run->bcl_dir;
	 $dir_bcl_tmp = "/data-dragen/bcl/".$project->getRun->run_name()."/";
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


my $reads;

###############
#SAMPLE SHEET 
###############

my $lines;
my $titles;
my $current_title;


foreach my $line(@$aoa){
		if ($line->[0] =~ /^\[/) {
			$current_title = $line->[0];
			push(@$titles,$current_title);
		}
		else {
			next unless $current_title;
			push(@{$lines->{$current_title}},$line);
		}	
}
unshift (@$titles, "[Settings]") unless grep(/^\[Settings\]$/, @$titles);
# todo: faire pareil avec [Header] et [Reads] ?

my $cb1_len;
my $cb2_len;

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
my $pos_sample = firstidx{ $_ eq "Sample_ID" } @$lheader_data;
die("no sample id in header ") if $pos_sample eq -1;
my $pos_sample_name = firstidx { $_ eq "Sample_Name" } @$lheader_data;

my $pos_cb1 = firstidx{ $_ eq "index" } @$lheader_data;
my $pos_cb2 = firstidx{ $_ eq "index2" } @$lheader_data;
my $len_cb;
 $len_cb->[0] = length($lines->{"[Data]"}->[0]->[$pos_cb1]);
 $len_cb->[1] = length($lines->{"[Data]"}->[0]->[$pos_cb2]); ;
 
foreach my $data (@{$lines->{"[Data]"}}){
	next unless  $data->[$pos_cb1];
	die($len_cb->[0]." ::  ".$data->[$pos_cb1]) if  $len_cb->[0] ne length($data->[$pos_cb1]);
	die("CB de taille differente") if  ($pos_cb2 ne '-1') and ($len_cb->[1] ne length($data->[$pos_cb2]));
}





my $guess_mask;
#if ($mask){

	my $config = XMLin("$bcl_dir/RunInfo.xml", KeyAttr => { reads => 'Reads' }, ForceArray => [ 'reads', 'read']);
	$reads = $config->{Run}->{Reads}->{Read};
	my $i_index =0;

	foreach  my $read (@$reads){
		if ($read->{IsIndexedRead} eq 'N'){
			push(@$guess_mask,"Y".$read->{NumCycles});
		}
		elsif ($read->{IsIndexedRead} eq 'Y'){
			my $l = $read->{NumCycles};
			if ($l == $len_cb->[$i_index]) {
				push(@$guess_mask,"I".$read->{NumCycles});
			}
			elsif ($l> $len_cb->[$i_index]){
				push(@$guess_mask,"I". $len_cb->[$i_index]."N".($l-$len_cb->[$i_index]));
			}
			else {
				die();
			}
			$i_index ++;
		}
		else {die;}
	}
	#warn Dumper $reads;
	
	#my $pos_sample_name = firstidx { $_ eq "Sample_Name" } @$lheader_data;
	
	#my $y = "Y".$read->[0]->{NumCycles};
	 #$mask = $y.";".$mask.";".$y; 
	
#}
my @real_mask;
my @index = (1,2);
@index = ();
if ($umi_name ){
	
	my @vmask = split(";",$mask);
	die() if scalar(@vmask) ne scalar(@$reads);
	for (my $i=0;$i< @vmask;$i++){
		my $l = $reads->[$i]->{NumCycles};
	
		
		my @mread = split("-",$vmask[$i]);
		my $nr = 0;
		foreach my $m1 (@mread){
			my ($v,$c) =split("",$m1);
			if ($v ne "Y"){
				$nr += $c;
			}
			if ($v eq "I"){
				push(@index,$i);
			}
		} 
		$l -= $nr;
		$vmask[$i] = join("",@mread) ;
		$vmask[$i]  =~ s/\*/$l/;
	}
	$mask = join(";",@vmask);	
}
else {
		foreach my $read (@$reads){
			if ($read->{IsIndexedRead} eq  "Y"){
				push(@index,$read->{Number});
			}
		}
}


$mask = join(";",@$guess_mask) unless $mask;


my $choice = prompt("use - ".colored(['bright_red on_black'],"$mask")." - for demultipexing  (y/n) ? ");
if ($choice ne "y") {
	$mask =  prompt("enter your mask  ? ");
	warn "use this mask $mask";
	#die($mask);
}

#die() unless -e $dir_in;


#change setting
foreach my $set (@{$lines->{"[Settings]"}}){
	if ($set->[0] eq "Adapter") {
		$set->[0] = "AdapterRead1"
	}
	
}


### read mask ;
my @amask = split(";",$mask);
#my $pos_umi = firstidx { $_ =~ /U/ } @amask;



### read mask ;
#my @amask = split(";",$mask);
#my $pos_umi = firstidx { $_ =~ /U/ } @amask;



if(scalar(@index) == 1){
	splice(@$lheader_data, $pos_cb2, 1);
	foreach my $data (@{$lines->{"[Data]"}}){
	if(scalar(@index) == 1){

		splice(@$data, $pos_cb2, 1);
		#my $pos_cb2 = firstidx { $_ eq "index2" } @$lheader_data;
	}
	}
}



$lines->{"[Settings]"} = [];
my $nb_mis = 0;
$nb_mis= $mismatch if $mismatch;
push(@{$lines->{"[Settings]"}},["BarcodeMismatchesIndex1",$nb_mis]);
push(@{$lines->{"[Settings]"}},["BarcodeMismatchesIndex2",$nb_mis]) if scalar(@index) == 2;
push(@{$lines->{"[Settings]"}},["OverrideCycles",$mask]) if $mask; 

my $dj;
my $ok_in_project;
foreach my $data (@{$lines->{"[Data]"}}){
	
	if($l2){
 		my $index2 =  $data->[$pos_cb2];
 		$data->[$pos_cb2] = substr($index2, 0, $l2);
 	}
	my $name = $data->[$pos_sample];
	$data->[$pos_sample_name] = $name;
	next if $name =~ /_RC$/; 
	next if exists $dj->{$name};
	next unless $name;
	$error_not_in_project->{$name} = $patients{$name}  unless exists $patients{$name};
	$ok_in_project->{$name} ++ unless not $patients{$name};
	$ok->{$name} ++;
 	delete $patients{$name};
 	$dj->{$name}++;
 	
	
}

my $error;
if (keys %$ok) {
	print colored(['bright_green on_black']," SAMPLES OK : ".scalar(keys %$ok) )."\n";
}
if (keys %patients) {
	print colored(['bright_red on_black']," SAMPLES IN PROJECT NOT IN SAMPLE SHEET :".scalar(keys %patients))."\n";
	map {print $_."\t".$patients{$_}."\n"} keys %patients;
	if ($run_name_option) {
		my $choice = prompt("continue anyway  (y/n) ? ");
		die() if ($choice ne "y"); 
	}
	else {
		$error =1;
	}
}

if (keys %$error_not_in_project) {
		print colored(['bright_red on_black']," SAMPLES IN SAMPLE SHEET NOT IN PROJECT : ".scalar(keys %$error_not_in_project) )."\n";
		map {print $_."\n"} keys %$error_not_in_project;
		print " nb error : ".scalar (keys %$error_not_in_project)."\n";
		my $choice = prompt("continue anyway  (y/n) ? ");
		die() if ($choice ne "y"); 
		#$lines->{"[Data]"} = change_sample_sheet($pos_sample,$lines->{"[Data]"},$error_not_in_project);
		
}

if ($error && $run_name_option) {
		my $choice = prompt("continue anyway  (y/n) ? ");
		die() if ($choice ne "y"); 
}
elsif  ($error) {die();}

unshift( @{$lines->{"[Data]"}},$lheader_data);
my $outcsv;

foreach my $title (@{$titles}){
	push(@$outcsv,[$title]);
	push(@$outcsv,@{$lines->{$title}});
}
my $ss = $bcl_dir."/file".time.".csv";
csv (in => $outcsv, out => $ss, sep_char=> ",");

sleep(1);

# sleep tant que le run n'est pas fini
my $complete = $bcl_dir.'RTAComplete.txt'; # MISEQ
$complete = $bcl_dir.'CopyComplete.txt' if ($bcl_dir =~ m{/(10X|ISEQ|NEXTSEQ500|NOVASEQ)/}); # 10X, ISEQ, NEXTSEQ500, NOVASEQ
warn $complete;
my $checkComplete = 1;
$checkComplete = 0 if -f $complete;
while($checkComplete == 1){
	warn "Run not complete, sleep 1h";
	sleep(3600);
	$checkComplete = 0 if (-f $complete);
}
system("mkdir $dir_bcl_tmp") unless -e $dir_bcl_tmp;
my $exit2 = system("rsync -rav --no-times $bcl_dir  $dir_bcl_tmp ");

my $cmd = qq{dragen --bcl-conversion-only=true --bcl-input-directory $dir_bcl_tmp --output-directory $dir_out --sample-sheet $ss --force --bcl-num-parallel-tiles 4   --bcl-num-conversion-threads 4   --bcl-num-compression-threads 4   --bcl-num-decompression-threads 4 };

my $exit = 0;
warn qq{$Bin/../run_dragen.pl -cmd="$cmd"};

$exit = system(qq{$Bin/../run_dragen.pl -cmd="$cmd"});
die() if $exit ne 0;
#exit(0);
warn "END DEMULTIPEX \n let's copy ";
my $fork =6;
my $pm   = new Parallel::ForkManager($fork);

foreach my $project_name (split(",",$project_names)){
	my $buffer = GBuffer->new();
	my $project = $buffer->newProject( -name 			=> $project_name );
	my $runs = $project->getRuns;
	my $run;
	if($run_name_option){
		($run) = grep{$_->plateform_run_name eq "$run_name_option"} @$runs;
	}
	else {
		$run = $runs->[0];
	}
	my $out_fastq = $run->fastq_dir();
	system("mkdir $out_fastq ; chmod g+rwx $out_fastq ");
	
	foreach my $p (@{$project->getPatients}){

		my $pid = $pm->start and next;

		my ($fastq1,$fastq2) = dragen_util::get_fastq_file($p,$out_fastq,$dir_out);
			warn $fastq1." ".$fastq2;



		#system ("rsync -rav $dir_out/".$p->name."_S* $out_fastq/");

		$pm->finish( 0, {});
	}

}

$pm->wait_all_children();
my $pr = $project_names;
$pr =~ s/,/_/g;
my $dir_stats = "/data-isilon/sequencing/ngs/demultiplex/".$run_name.".".$pr;


system("mkdir -p $dir_stats ;chmod -R a+rwx $dir_stats; rsync -rav ".$dir_out."/Reports/ $dir_stats/ ; chmod -R a+rwx $dir_stats;");

report($dir_stats."/Demultiplex_Stats.csv");

exit(0);
###

sub report {
	my ($csv) = @_;
	confess($csv) unless -e $csv;
	my $aoa = csv (in => $csv); 
	my $header = shift @$aoa;
	
	my $byline;
	my $bypatient;
	
	foreach my $line (@$aoa){
		my $l = $line->[0];
		$byline->{$l} += $line->[3];
		my $p = $line->[1];
		$bypatient->{$p} += $line->[3];;
	}
	
print "-- BY LINES : --\n";	
my $tb = Text::Table->new( (colored::stabilo("blue", "Lane" , 1),  "# read") ) ; # if ($type == 1);
my @l ;
my @rows;
my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(values %$byline);
my $mean = $stat->mean();
my $sd = $stat->standard_deviation();
foreach my $l (keys %$byline){
		my @row;
			push(@row,colored::stabilo("blue",$l,1));
			
			if ($byline->{$l} < 100){
				push(@row,colored::stabilo("red",$byline->{$l},1)) ;
			}
			elsif ($byline->{$l} < (3*$sd) + $mean){
				push(@row,colored::stabilo("red",$byline->{$l},1)) ;
			}
			elsif ($byline->{$l} < $sd + $mean){
				push(@row,colored::stabilo("orange",$byline->{$l},1)) ;
			}
			
			else {
					push(@row,colored::stabilo("green",$byline->{$l},1)) ;
			}
			push(@rows,\@row);
		}
		
		$tb->load(@rows);
		print $tb;
		print "\n ------------------------------\n";

print "-- By Sample : --\n";	
my $tb2 = Text::Table->new( (colored::stabilo("blue", "Sample" , 1),  "# read") ) ; # if ($type == 1);
 @rows = ();
foreach my $l (keys %$bypatient){
		my @row;
			push(@row,colored::stabilo("blue",$l,1));
			if ($byline->{$l} > 1000){
				push(@row,colored::stabilo("green",$bypatient->{$l},1)) ;
			}
			else {
					push(@row,colored::stabilo("red",$bypatient->{$l},1)) ;
			}
			push(@rows,\@row);
		}
		
		$tb2->load(@rows);
		print $tb2;
		print "\n ------------------------------\n";
}

	
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


