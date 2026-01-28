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
use lib "$Bin/../../dragen/scripts/";
use Text::CSV qw( csv );
use Getopt::Long;
use Data::Dumper;
use List::MoreUtils qw(firstidx);
use colored; 
use IO::Prompt;
use Sys::Hostname;
use Term::ANSIColor;
use Moo;
use GBuffer;
use GenBoProject;
use Text::Table;
use file_util;
use File::Temp qw/ tempfile tempdir /;; 
use Term::Menus;
use Proc::Simple;
use dragen_util;
use XML::Simple qw(:strict);
use Statistics::Descriptive;
use Bio::SeqIO;
use Carp;
use JSON::XS;
use File::Path qw(make_path);
use Pod::Usage;

my $project_names;
my $run_name_option;
my $mismatch = 0;
my $create_fastq_umi;
my $no_demux_only;

GetOptions(
	'projects=s'			=> \$project_names,
	'run=s'					=> \$run_name_option,
	'mismatches=i'			=> \$mismatch,
	'no_demux_only'			=> \$no_demux_only,
	'create_fastq_umi'		=> \$create_fastq_umi,
	'help|?'				=> sub {pod2usage(-verbose => 2,-noperldoc=>1)},
) || confess("Possibles options are: projects, run, mimsatches, create_fastq_umi, help.
Use $0 --help to see the full documentation.\n");

# Vérifications options
pod2usage(
    -message => "Error : --projects option is mandatory\n",
    -exitval => 2
) unless $project_names;

if ( $mismatch !~ /^[012]$/ ) {
    pod2usage(
        -message => "Error : --mismatches must be 0, 1 ou 2\n",
        -exitval => 2
    );
}
 
my $run_name;
my $umi_mask;
my %patients;
my $neb;
my $adaptors;
my $bcl_dir;
my $dir_tmp;
my $bcl_tmp;
my $dir_out;
my $aoa;
foreach my $project_name (split(",",$project_names)){
	my $buffer = GBuffer->new();
	my $project = $buffer->newProject( -name => $project_name );
	my $runs = $project->getRuns;
	my $run;
	if (scalar(@$runs)> 1){
		unless ($run_name_option){
			print colored(['bright_red on_black'],"You have ".scalar(@$runs)." runs in the project $project_name. You have to choose one and add -run on the command line" )."\n";
			foreach my $run (sort {$a->plateform_run_name cmp $b->plateform_run_name} @$runs){
				print colored(['bright_green on_black'],$run->plateform_run_name)." ".colored(['bright_blue on_black'],$run->date)."\n";
			}
			die();
		}
		($run) = grep{$_->plateform_run_name eq "$run_name_option"} @$runs;
		die("Unable to find $run_name_option ".$project_name.":\n\t".join("\n\t", sort map{$_->plateform_run_name} @$runs)) unless ($run);
	}
	else {
		$run = $runs->[0];
	}
	confess("This run is not from an Aviti, but ".$run->machine) unless ($run->machine eq 'AVITI');
	$run_name = $run->run_name();
	warn $run_name;

	
	# UMI mask
	foreach my $capture (@{$project->getCaptures}){
		if ($capture->umi){
		 $umi_mask = $capture->umi->{mask} unless (defined $umi_mask);
		 die("Problem more than one mask ") if ($umi_mask ne $capture->umi->{mask});
		}
	}
	
	# RNAseq NEB: récupère les adaptateurs à trimmer
	my @profiles = map {$_->getSampleProfile} @{$project->getPatients};
	$neb = 1 if (grep{/neb/i} @profiles);
	if ($neb) {
		undef $adaptors;
		my $illuminaAdaptors = Bio::SeqIO->new(-file => $project->getIlluminaAdaptors, '-format' => 'Fasta');
		my $tsoAdaptors = Bio::SeqIO->new(-file => $project->getTsoAdaptors, '-format' => 'Fasta');
		foreach my $fasta ($illuminaAdaptors,$tsoAdaptors) {
			while ( my $seq = $fasta->next_seq() ) {
				push(@$adaptors,$seq->seq);
			}
		}
	}
	
	
	map {$patients{$_->name} = $project->name} @{$project->getPatients};

	$bcl_dir = $run->bcl_dir;
	die("Can't find bcl directory : ".$bcl_dir) unless (-d $bcl_dir);
	die("Problem different bcl dir : $bcl_dir ".$run->bcl_dir) if ($bcl_dir && $run->bcl_dir() ne $bcl_dir);
	$dir_tmp = $buffer->config_path("root","tmp");
	$bcl_tmp = $dir_tmp.'bcl/'.$run_name.'/';
	make_path($bcl_tmp) unless (-d $bcl_tmp);
	$dir_out = $dir_tmp.'fastq/'.$run_name.'/';
	make_path($dir_out) unless (-d $dir_out);
	
	# Run Manifest
#	if ($project->getCapture->name =~ /^transcriptome_10X/) {
#		die("No SampleSheet10X.csv") unless (-f $bcl_dir.'SampleSheet10X.csv');
#		$aoa = csv (in => $bcl_dir.'SampleSheet10X.csv');
#		next;
#	}
	my $csv_tmp = $bcl_dir."/RunManifest.tmp.csv";
	my $sp = $run->sample_sheet;
	if (not $sp and -e $bcl_dir."/RunManifest.csv") {
		my $prompt = prompt("No RunManifest in database. Use '$bcl_dir/RunManifest.csv' instead ? (y/n) ", -yes_no);
		$aoa = csv (in => $bcl_dir."/RunManifest.csv" ) if ($prompt);
		next
	}
	die("No RunManifest in database ") unless ($sp or -e $bcl_dir.'SampleSheet10X.csv');
	$sp =~ s/;/,/g;	
	open(my $SP,">$csv_tmp") or die("Can't open '$csv_tmp': $!");
	print {$SP} $sp;
	close ($SP);
	$aoa = csv (in => $csv_tmp);
	system("rm $csv_tmp") if (-e $csv_tmp);
}



###############
# RUN MANIFEST 
###############


my $current_title;
my $titles;
my $lines;
foreach my $line(@$aoa){
	if ($line->[0] =~ /^\[/) {
		$current_title = $line->[0];
		push(@$titles,uc($current_title));
	}
	else {
		next if ($line->[0] =~ /^\s*#/);
		next unless $current_title;
		push(@{$lines->{$current_title}},$line);
	}	
}
unshift (@$titles, "[SETTINGS]", ["SettingName","Value"]) unless grep(/^\[SETTINGS\]/i, @$titles);

my $cb1_len;
my $cb2_len;


# Samples section
die("Missing [SAMPLES] header/section in RunManifest") unless (exists $lines->{"[SAMPLES]"});
my $lheader_data = shift @{$lines->{"[SAMPLES]"}};
if (scalar (@$lheader_data) ne  scalar (@{$lines->{"[SAMPLES]"}->[0]})){
	my $tb = Text::Table->new(
		{is_sep => 1},
		@$lheader_data,
    );
    $tb->load($lines->{"[SAMPLES]"}->[0]);
    print $tb;
	print "\n";
	die(scalar (@$lheader_data)." fields in the Samples header vs " .scalar (@{$lines->{"[SAMPLES]"}->[0]})." in the 1rst line of Samples section");
}

my $pos_sample_name = firstidx { $_ eq "SampleName" } @$lheader_data;
die("No SampleName in [SAMPLES] section header:\n".join(',',@$lheader_data)) if ($pos_sample_name eq -1);

my $pos_cb1 = firstidx{ $_ eq "Index1" } @$lheader_data;
my $pos_cb2 = firstidx{ $_ eq "Index2" } @$lheader_data;
my $len_cb;
$len_cb->[0] = length($lines->{"[SAMPLES]"}->[0]->[$pos_cb1]) unless ($pos_cb1 eq '-1');
$len_cb->[1] = length($lines->{"[SAMPLES]"}->[0]->[$pos_cb2]) unless ($pos_cb2 eq '-1');

foreach my $data (@{$lines->{"[SAMPLES]"}}){
	next unless ($data->[$pos_cb1]);
	die("CB1 de tailles differentes: ".$len_cb->[0]." :: ".$data->[$pos_sample_name].' '..$data->[$pos_cb1]) if  ($pos_cb1 ne '-1') and ($len_cb->[0] ne length($data->[$pos_cb1]));
	die("CB2 de tailles differentes: ".$len_cb->[1].' :: '.$data->[$pos_sample_name].' '.$data->[$pos_cb2]) if  ($pos_cb2 ne '-1') and ($len_cb->[1] ne length($data->[$pos_cb2]));
}



# Mask
open(my $runParam, '<', "$bcl_dir/RunParameters.json") or die ("Can't open '$bcl_dir/RunParameters.json': $!");
my $json;
while (my $line = <$runParam>) {
	chomp $line;
	$json .= $line;
}
close($runParam);
my $config = decode_json $json;
my $cycles = $config->{"Cycles"};
my $lanes = $config->{"AnalysisLanes"};
warn 'Lanes: '.$lanes; 
warn 'Cycles:'.Dumper $cycles;

my $mask;
my $i_index = 0;
foreach my $read (keys %$cycles){
	if ($read =~ /^R/){
		$mask->{$read} = "Y".$cycles->{$read};
	}
	elsif ($read =~ /^I/){
		my $l = $cycles->{$read};
		if ($l == $len_cb->[$i_index]) {
			$mask->{$read} = "Y".$cycles->{$read};
		}
		elsif ($l > $len_cb->[$i_index]){
			$mask->{$read} = "Y". $len_cb->[$i_index]."N".($l-$len_cb->[$i_index]);
		}
		else {
			die("Error: length of BC $read (".$len_cb->[$i_index].") is larger than the number of cycles ($l)");
		}
		$i_index ++;
	}
	else {
		die "Error guessing the mask:\n". Dumper $cycles;
	}
}

if ($umi_mask){
	warn $umi_mask;
	my @vmask = split(";",$umi_mask);
	die("Error: mask $mask doesn't fit the run cycles ".join(',', keys %$cycles ) ) if scalar(@vmask) ne scalar(keys %$cycles);
	my ($r,$i) = (0,0);
	undef $mask;
	foreach my $n (0 .. $#vmask) {
		if ($n == 0 or $n == $#vmask) {
#		if ($vmask[$n] =~ /Y/) {
			$r ++;
#			$mask->{'R'.$r} = 'R'.$r.':';
			if ($vmask[$n] =~ /U/) {
				$mask->{'Umi'} .= '-' if ($mask->{'Umi'});
				$mask->{'Umi'} .= 'R'.$r.':';
			}
			foreach my $m1 (split("-",$vmask[$n])) {
				my ($v,$c) = split("",$m1);
				$mask->{'R'.$r} .= 'Y'.$c if ($v eq "Y");
				$mask->{'R'.$r} .= 'N'.$c if ($v ne "Y");
				$mask->{'Umi'} .= 'Y'.$c if ($v eq "U");
				$mask->{'Umi'} .= 'N'.$c if ($v ne "U");
			}
		}
		else {
#		elsif ($vmask[$n] =~ /I/) {
			$i ++;
#			$mask->{'I'.$i} = 'I'.$i.':';
			if ($vmask[$n] =~ /U/) {
				$mask->{'Umi'} .= '-' if ($mask->{'Umi'});
				$mask->{'Umi'} .= 'I'.$i.':';
			}
			foreach my $m1 (split("-",$vmask[$n])) {
				my ($v,$c) = split("",$m1);
				$mask->{'I'.$i} .= 'Y'.$c if ($v eq "I");
				$mask->{'I'.$i} .= 'N'.$c if ($v ne "I");
				$mask->{'Umi'} .= 'Y'.$c if ($v eq "U");
				$mask->{'Umi'} .= 'N'.$c if ($v ne "U");
			}
		}
	}
	foreach my $key (keys %$mask) {
		while ($mask->{$key} =~ /([YN])(\d|\*)\g1(\d|\*)/) {
			my $s = '*';
			$s = $2 + $3 unless ($2 eq '*' or $3 eq '*');
			$mask->{$key} =~ s/([YN])(\d|\*)\g1(\d|\*)/$1$s/g;
		}
	}
}
warn 'Mask:'.Dumper $mask;


unless ( prompt( "Use - " . colored(['bright_red on_black'], join(",", map {$_.':'.$mask->{$_}} sort keys %$mask)) . " - for demultipexing  (y/n) ? ", -yes_no) ) {
	undef $mask;
	my $new_mask = uc(prompt("Enter your mask: "));
	$new_mask =~ s/ //;
#	warn "Use this mask: ".$new_mask;
	my $motif_mask = '([RI][12]):([YN](?:\d+|\*))+';
	my $motif_mask_umi = '(Umi):([RI][12]:(?:[YN](?:\d+|\*))+(?:-[RI][12]:(?:[YN](?:\d+|\*))+)*)';
	my $re = qr{^(?:$motif_mask|$motif_mask_umi)(?:,(?:$motif_mask|$motif_mask_umi))*$};
	die("Error in mask entered") unless ($new_mask =~ $re);
	$mask->{$1} = $2 while ($new_mask =~ /([RI][12]|Umi):((?:[YN](?:\d+|\*))+)/g);
	warn "Using this mask: ".colored(['bright_red on_black'], join(",", map {$_.':'.$mask->{$_}} sort keys %$mask) );
	sleep(2);
}

# todo: check that the mask matches with the cycles ?


my $settings;
shift @{$lines->{"[SETTINGS]"}};
foreach my $line (@{$lines->{"[SETTINGS]"}}) {
	$settings->{$line->[0]} = $line->[1];
}

$lines->{"[SETTINGS]"} = [['SettingName','Value','Lane']];
$settings->{"I1MismatchThreshold"} = $mismatch;
$settings->{"I2MismatchThreshold"} = $mismatch if (grep (/Index2/,@$lheader_data));
$settings->{"I1Mask"} = $mask->{'I1'};
$settings->{"I2Mask"} = $mask->{'I2'};
$settings->{"R1FastQMask"} = $mask->{'R1'};
$settings->{"R2FastQMask"} = $mask->{'R2'};
$settings->{"UmiMask"} = $mask->{'Umi'} if (exists $mask->{'Umi'});
$settings->{"UmiFastQ"} = 'True' if ($create_fastq_umi);
if ($neb) {	# RNAseq NEB: Trim les adaptateurs
	$settings->{"R1Adapter"} .= '+' if (exists $settings->{"R1Adapter"} and $settings->{"R1Adapter"});
	$settings->{"R2Adapter"} .= '+' if (exists $settings->{"R2Adapter"} and $settings->{"R2Adapter"});
	$settings->{"R1Adapter"} .= join('+',@$adaptors);
	$settings->{"R2Adapter"} .= join('+',@$adaptors);
	$settings->{"R1AdapterTrim"} = 'True';
	$settings->{"R1AdapterTrim"} = 'True';
}
foreach my $param (sort keys %$settings) {
	push(@{$lines->{"[SETTINGS]"}},[$param,$settings->{$param},$lanes]);
}


my $dj;
my $ok_in_project;
my $error_not_in_project = {};
my $ok_in_runmanifest;
foreach my $data (@{$lines->{"[SAMPLES]"}}){
	my $name = $data->[$pos_sample_name];
	next if ($name =~ /_RC$/); 
	next if ($name =~ /^PhiX$/);
	next if (exists $dj->{$name});
	next unless ($name);
	$error_not_in_project->{$name} = $patients{$name} unless (exists $patients{$name});
	$ok_in_project->{$name} ++ if ($patients{$name});
	$ok_in_runmanifest->{$name} ++;
 	delete $patients{$name};
 	$dj->{$name}++;
}

if (keys %$ok_in_project) {
	print colored(['bright_green on_black']," Samples OK in project(s) and RunManifest: ".scalar(keys %$ok_in_project) )."\n";
}
if (keys %$ok_in_runmanifest) {
	print colored(['bright_green on_black']," Samples OK in RunManifest: ".scalar(keys %$ok_in_runmanifest) )."\n";
}
if (keys %patients) {
	print colored(['bright_red on_black']," Samples OK in project(s) NOT in RunManifest:".scalar(keys %patients))."\n";
	map {print $_."\t".$patients{$_}."\n"} keys %patients;
	print "nb error: ".scalar (keys %patients)."\n";
	die() unless (prompt("continue anyway ? (y/n)  ", -yes_no));
}
if (keys %$error_not_in_project) {
	print colored(['bright_yellow on_black']," Samples in RunManifest NOT in project: ".scalar(keys %$error_not_in_project) )."\n";
	map {print $_."\n"} keys %$error_not_in_project;
	print " nb error: ".scalar (keys %$error_not_in_project)."\n";
	die() unless (prompt("continue anyway ? (y/n)  ", -yes_no));
		
}

unshift(@{$lines->{"[SAMPLES]"}},$lheader_data);
my $outcsv;
foreach my $title (@{$titles}){
	push(@$outcsv,[$title]);
	push(@$outcsv,@{$lines->{$title}});
}
my $runm_name = "RunManifest".time.".csv";
my $ss = $bcl_dir.$runm_name;
csv (in => $outcsv, out => $ss, sep_char=> ",");
warn $ss;


###############
# Demultipexage
###############

# sleep tant que le run n'est pas fini
my $complete = $bcl_dir.'RunUploaded.json';
warn $complete;
my $checkComplete = 1;
$checkComplete = 0 if -f $complete;
while($checkComplete == 1){
	my ($sec, $min, $hour) = (localtime)[0,1,2];
	my $time = sprintf("%02d:%02d:%02d", $hour, $min, $sec);
	warn "Run not complete, sleep 1h (3600 s) at $time";
	sleep(3600);
	$checkComplete = 0 if (-f $complete);
}

# Copy bcl
#my $cmd_rsync = "rsync -rah --no-times --size-only $bcl_dir $bcl_tmp" =~ s/\/\//\//rg;
#warn $cmd_rsync;
#my $exit_rsync = system($cmd_rsync);
#die("Rsync error, please retry") if ($exit_rsync);
#my $ss1 = $bcl_tmp.$runm_name;

# Demultiplex command
my $cmd = "singularity run -B $bcl_dir -B $dir_tmp /software/distrib/BASE2FASTQ/bases2fastq.2.2.sif bases2fastq ";
$cmd .= "--error-on-missing ";
#$cmd .= "--r2-cycles 58 " if (getpwuid($<) eq 'mperin'); # si le run n'est pas fini
#$cmd .= "--settings 'I1MismatchThreshold,$mismatch' --settings 'I2MismatchThreshold,$mismatch' ";
#$cmd .= "--settings 'I1Mask,".$mask->{'I1'}."' --settings 'I2Mask,".$mask->{'I2'}."' ";
#$cmd .= "--settings 'R1FastQMask,".$mask->{'R1'}."' --settings 'R2FastQMask,".$mask->{'R2'}."' ";
#$cmd .= "--settings 'UmiMask,".$mask->{'Umi'}."' ";
$cmd .= "--run-manifest $ss --num-unassigned 500 --num-threads 40 "; # $ss1
$cmd .= "--group-fastq --no-projects ";
$cmd .= "$bcl_dir $dir_out"; # bcl_tmp

# Demux only
unless ($no_demux_only) {
	my $out_demux_only = $dir_out.'demux_only/'.time.'/';
	make_path($out_demux_only) unless (-d $out_demux_only);
	my $cmd_demux_only = $cmd;
	$cmd_demux_only =~ s/--num-threads 40/--num-threads 5 --demux-only/;
	$cmd_demux_only =~ s/$dir_out/$out_demux_only/;
	warn $cmd_demux_only;
	my $exit = system('time '.$cmd_demux_only);
	confess ("Error while demux only") if ($exit);
	
	# Open demutiplex stats
	opendir(my $dh, $out_demux_only) || die ("Can't opendir '$out_demux_only': $!");
	my $demux_qc_html;
	while (my $file = readdir($dh)) {
		$demux_qc_html = $out_demux_only.$file if ($file =~ /_QC\.html$/);
	}
	close ($dh);
	system("firefox $demux_qc_html &") unless (getpwuid($<) eq 'shanein');
	system("google-chrome $demux_qc_html &") if (getpwuid($<) eq 'shanein');
	die("\nbye") if (prompt("Check the demux only stats. Then, tap any key to continue, 'q' to quit.", -w=>'q'));
	print "\n";
}

# Complete demultiplex
warn $cmd;
sleep(2);
my $exit = system("echo $cmd | run_cluster.pl -cpu=40");
die("Error while demultiplexing") if ($exit);

# Copy fastq
warn "END DEMULTIPEX \n let's copy ";
my $fork = 6;
my $pm = new Parallel::ForkManager($fork);

foreach my $project_name (sort split(",",$project_names)){
	my $buffer = GBuffer->new();
	my $project = $buffer->newProject( -name => $project_name );
	my $runs = $project->getRuns;
	my $run;
	if($run_name_option){
		($run) = grep{$_->plateform_run_name eq "$run_name_option"} @$runs;
	}
	else {
		$run = $runs->[0];
	}
	my $out_fastq = $run->fastq_dir();
	make_path($out_fastq, {mode => 0774}) unless (-d $out_fastq);
	
	foreach my $p (@{$project->getPatients}){
		my $pid = $pm->start and next;
		my $tmp_fastq = $dir_out.'Samples/';
#		if (grep (/Project/, @$lheader_data)) { # si pas --no-projects --group-fastq
#			opendir(my $dh, $dir_out.'Samples/') || die ("Can't open '$dir_out/Samples/': $!");
#			while (my $dir = readdir($dh)) {
#				next if ($dir =~ /^(\.|PhiX|Unassigned$)/);
#				next unless (-d $dir);
#				$tmp_fastq = $dir_out.'Samples/'.$dir.'/'.$p->name.'/' if (-d $dir_out.'Samples/'.$dir.'/'.$p->name.'/');
#			}
#			close($dh);
#		}
		my ($fastq1,$fastq2) = dragen_util::get_fastq_file($p,$out_fastq,$tmp_fastq);
		warn $fastq1." ".$fastq2;
		$pm->finish( 0, {});
	}

}
$pm->wait_all_children();

# Demultiplex stats
my $pr = $project_names;
$pr =~ s/,/_/g;
my $dir_stats = "/data-isilon/sequencing/ngs/demultiplex/".$run_name.".".$pr.'/';
make_path($dir_stats,{mode=>0777}) unless (-d $dir_stats);
system("rsync -a $dir_out $dir_stats ; chmod -R a+rwx $dir_stats");
warn("rsync $dir_out* $dir_stats ; chmod -R a+rwx $dir_stats");
system("rsync $dir_out* $dir_stats ; chmod -R a+rwx $dir_stats");

report($dir_stats."/IndexAssignment.csv");

# Open demutiplex stats
opendir(my $dh, $dir_stats) || die ("Can't opendir '$dir_stats': $!");
my $demux_qc_html;
while (my $file = readdir($dh)) {
	$demux_qc_html = $dir_stats.$file if ($file =~ /_QC\.html$/);
}
close ($dh);
system("firefox $demux_qc_html") unless (getpwuid($<) eq 'shanein');
system("google-chrome $demux_qc_html") if (getpwuid($<) eq 'shanein');


exit(0);

###

sub report {
	my ($csv) = @_;
	confess("No '$csv'") unless (-e $csv);
	my $aoa = csv (in => $csv); 
	my $header = shift @$aoa;
	
	my $byline;
	my $bypatient;
	
	# read csv
	foreach my $line (@$aoa){
		# SampleNumber,SampleName,I1,I2,NumPoloniesAssigned,PercentPoloniesAssigned,Yield(Gb),Lane
		my $l = $line->[7];
		$byline->{$l}->{'count'} += $line->[4];
		$byline->{$l}->{'percent'} += $line->[5];
		my $p = $line->[1];
		unless ($l eq '1+2') {
			$bypatient->{$p}->{'count'} += $line->[4];
			$bypatient->{$p}->{'percent'} += $line->[5];
		}
	}
	
	print "-- By Lines : --\n";	
	my $tb = Text::Table->new( (colored::stabilo("blue", "Lane" , 1),  "# read") ) ; # if ($type == 1);
	my @l ;
	my @rows;
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(values %$byline);
	my $mean = $stat->mean();
	my $sd = $stat->standard_deviation();
	foreach my $l (sort keys %$byline){
		my @row;
		push(@row,colored::stabilo("blue",$l,1));
		
		if ($byline->{$l} < 100){
			push(@row,colored::stabilo("red",$byline->{$l}->{'count'}.' ('.sprintf("%.2f",$byline->{$l}->{'percent'}).'%)',1)) ;
		}
		elsif (abs( $byline->{$l} - $mean ) > 3*$sd){
			push(@row,colored::stabilo("red",$byline->{$l}->{'count'}.' ('.sprintf("%.2f",$byline->{$l}->{'percent'}).'%)',1)) ;
		}
		elsif (abs($byline->{$l} - $mean) > $sd){
			push(@row,colored::stabilo("yellow",$byline->{$l}->{'count'}.' ('.sprintf("%.2f",$byline->{$l}->{'percent'}).'%)',1)) ;
		}
		else {
			push(@row,colored::stabilo("green",$byline->{$l}->{'count'}.' ('.sprintf("%.2f",$byline->{$l}->{'percent'}).'%)',1)) ;
		}
		push(@rows,\@row);
	}
	$tb->load(@rows);
	print $tb;
	print "\n ------------------------------\n";
	
	print "-- By Samples : --\n";	
	my $tb2 = Text::Table->new( (colored::stabilo("blue", "Sample" , 1),  "# read") ) ; # if ($type == 1);
	@rows = ();
	foreach my $p (sort keys %$bypatient){
		my @row;
		push(@row,colored::stabilo("blue",$p,1));
		if ($bypatient->{$p} < 1000000){
			push(@row,colored::stabilo("red",$bypatient->{$p}->{'count'}.' ('.sprintf("%.2f",$bypatient->{$p}->{'percent'}).'%)',1)) ;
		}
		elsif (abs( $bypatient->{$p} - $mean ) > 3*$sd){
			push(@row,colored::stabilo("red",$bypatient->{$p}->{'count'}.' ('.sprintf("%.2f",$bypatient->{$p}->{'percent'}).'%)',1)) ;
		}
		elsif (abs( $bypatient->{$p} - $mean ) > $sd){
			push(@row,colored::stabilo("yellow",$bypatient->{$p}->{'count'}.' ('.sprintf("%.2f",$bypatient->{$p}->{'percent'}).'%)',1)) ;
		}
		else {
			push(@row,colored::stabilo("green",$bypatient->{$p}->{'count'}.' ('.sprintf("%.2f",$bypatient->{$p}->{'percent'}).'%)',1)) ;
		}
		push(@rows,\@row);
	}
	$tb2->load(@rows);
	print $tb2;
	print "\n ------------------------------\n";
}

	
sub usage {
	print "
$0
-------------
Obligatoires:
	projects <s>               nom de projet(s) séparés par une virgule
Optionels:
	run <s>                    nom du run à démultiplexer si un projet est sur plusieurs runs
	mismatches <i>             nombre de mismatch(es) à autoriser,
	                           valeurs possibles: 0,1,2, défaut: 0
	no_demux_only              ne pas faire l'étape demux only, lancer directement le démultiplexage complet
	create_fastq_umi           générer des fastq pour les UMI.
	help                       affiche un message d'aide détaillé

";
	exit(1);
}


__END__

=pod

=encoding UTF-8

=head1 NAME

aviti_demultiplex.pl - Lance le démultiplexage des runs AVITI via bases2fastq

=head1 SYNOPSIS

  perl aviti_demultiplex.pl --projects <projet1,projet2> [options]

=head1 DESCRIPTION

Ce script permet de lancer le démultiplexage de runs issus d'un séquenceur
AVITI en utilisant l'outil B<bases2fastq>.

=head1 REQUIRED OPTIONS

=over 4

=item B<--projects> I<string>

Nom ou liste de noms de projets à démultiplexer, séparés par des virgules.

Exemple :

  --projects NGS20XX_XXXXX,NGS20YY_YYYYY

=back

=head1 OPTIONAL OPTIONS

=over 4

=item B<--run> I<string>

Nom du run à démultiplexer.  
À utiliser lorsque le projet est présent sur plusieurs runs.

=item B<--mismatches> I<int>

Nombre de mismatches autorisés pour l'index.

Valeurs possibles : 0, 1, 2  
Valeur par défaut : 0

=item B<--no_demux_only>

Désactive l'étape demux only, et lance directement le démultiplexage complet.

=item B<--create_fastq_umi>

Active la génération de fichiers FASTQ pour les UMI.

=back

=head1 EXAMPLES

Démultiplexage simple d'un projet :

  perl aviti_demultiplex.pl --projects NGS20XX_XXXXX

Démultiplexage d'un projet avec 0 mismatch :

  perl aviti_demultiplex.pl --projects NGS20XX_XXXXX --mismatches 0

Démultiplexage avec génération des FASTQ UMI :

  perl aviti_demultiplex.pl --projects NGS20XX_XXXXX,NGS20YY_YYYYY --create_fastq_umi

=cut

