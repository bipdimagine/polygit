#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../../../GenBo/lib/GenBoDB/writeDB";
use lib "$Bin/../../packages";
use lib "$Bin/../../../polyprojectNGS/cgi-bin/";
use GBuffer;
use GenBoProject;
use Getopt::Long;
use Data::Dumper;
#use illumina_util;
use IO::Prompt;
use Sys::Hostname;
use File::Find::Rule ;
use queryPolyproject; 
 

use colored; 
use file_util;
use check_utils;
use Text::Table;
use Term::Twiddle;
my $projectName;
my $filename;

my $name;
my $patients_name;
my $steps_name;
my $force;
my $type;
my $fastq_ext;
GetOptions(
	'project=s' => \$projectName,

);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
my $ped_ori =  $project->getPedigreeFile();
my @no_moove;
my $moove;
my $dbh = $buffer->dbh;
$dbh->{AutoCommit} = 0;

ExtractPedigreeSection($project->name());


sub ExtractPedigreeSection {
	my ($project_name) =@_;
	my $buffer = GBuffer->new;
	my $project = $buffer->newProject(-name => $project_name);
#	$project->check_ped_db_only(1);
	die( "unknown project" . $project_name ) unless $project->id();
	my $file_out = $project->getPedigreeFile();
  	my $patientList = queryPolyproject::getPatientProjectInfo($buffer->dbh,$project->id());
  	warn Dumper ($patientList);
	my %fat;
	foreach my $f (@$patientList) {
		warn Dumper $f;
		$f->{family}=$f->{name} unless $f->{family};
		$fat{$f->{family}."_".$f->{father}}=1 if $f->{father};
	}
	my %mot;
	foreach my $m (sort {$b->{family} cmp $a->{family}}@$patientList) {
		$m->{family}=$m->{name} unless $m->{family};
		$mot{$m->{family}."_".$m->{mother}}=1 if $m->{mother};
	}
#	seek $FHO,0,0;
	my @data;
	my %seen;
	foreach my $p (@$patientList) {
		my %s;
		my $fatnew=0;
		my $motnew=0;
		$fatnew=1 unless $p->{father};
		$motnew=1 unless $p->{mother};
		$p->{father}=$p->{family}."_Pere" unless $p->{father};
		$p->{mother}=$p->{family}."_Mere" unless $p->{mother};		
		$p->{father}=0 if (($fat{$p->{family}."_".$p->{name}}));
		$p->{mother}=0 if (($mot{$p->{family}."_".$p->{name}}));
		$p->{father}=0 unless  $p->{mother};
		$p->{mother}=0 unless  $p->{father};

		$s{family} = $p->{family};
		$s{name} = $p->{name};
		$s{father}=$p->{father};
		$s{mother}=$p->{mother};
		$s{sex}=$p->{sex};
		$s{sex}=1 unless $p->{sex};
		$s{status}=$p->{status};
		$s{status}=2 unless $p->{status};		
		my $father=$s{father};
		my $mother=$s{mother};
		push(@data,\%s);

		my %a;
		if ($father eq $p->{family}."_Pere" && $fatnew) {
			$a{family} = $p->{family};
			$a{name} = $p->{family}."_Pere";
			$a{father}=0;
			$a{mother}=0;
			$a{sex}=1;
			$a{status}=1;
			push(@data,\%a) unless $seen{$a{name}}++;
		}
		my %b;
		if ($mother eq $p->{family}."_Mere" && $motnew) {
			$b{family} = $p->{family};
			$b{name} = $p->{family}."_Mere";
			$b{father}=0;
			$b{mother}=0;
			$b{sex}=2;
			$b{status}=1;
			push(@data,\%b) unless $seen{$b{name}}++;
		}
	}
	my @result_sorted=sort {$b->{family} cmp $a->{family}||$b->{father} cmp $a->{father}||$b->{mother} cmp $a->{mother}} @data;
	my $dlo="\t";
	my $dir = $project->getRootDir;
	warn $file_out;
	logrotate_File("$file_out",$dir);
	open( my $FHO, '>', $file_out ) or die("Can't create file: $file_out\n");
	foreach my $l (@result_sorted) {
		print $l->{family}.$dlo.$l->{name}.$dlo.$l->{father}.$dlo.$l->{mother}.$dlo.$l->{sex}.$dlo.$l->{status}."\n";
		print $FHO $l->{family}.$dlo.$l->{name}.$dlo.$l->{father}.$dlo.$l->{mother}.$dlo.$l->{sex}.$dlo.$l->{status}."\n";
	}
	chmod 0777, $FHO;
	close($FHO);
#	sendOK("Successful validation for Project Patient Ped file: <b>". $project_name.".ped</b>" );
}

sub logrotate_File {
	my ($sub,$dir_backup) = @_;
	my $file_out = $dir_backup.'/'.$sub;
	if (-e $file_out) {
		my $merge_rotate = new Logfile::Rotate (
				File  => $file_out,
				Count => 9,
				Dir	  => $dir_backup,
				Gzip  => 'no',
				 );		
		$merge_rotate->rotate();
		unlink $file_out;
	}
}