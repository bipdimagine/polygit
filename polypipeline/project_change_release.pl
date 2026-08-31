#!/usr/bin/env perl

use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use lib "$RealBin/../GenBo/lib/obj-nodb/";
use lib "$RealBin/../GenBo/lib/obj-nodb/";
#use lib "/data-isilon/bipd-src/pnitschk/git-master/polygit/GenBo/lib/obj-nodb";

use GBuffer;
use Set::IntSpan::Fast::XS;
use Getopt::Long;
use Digest::MD5 qw(md5 md5_hex md5_base64);
my ($project_name, $patient_name,$vid);
my $file;
my $release;
GetOptions(
	'project=s'    => \$project_name,
	'release=s' => \$release,
);

my $buffer = new GBuffer;
my @list = split(",",$project_name);
foreach my $pname (@list){
my $project = $buffer->newProject( -name => $pname);
my $query = $buffer->getQuery();
$buffer->dbh->{'AutoCommit'} = 0;
foreach my $p (@{$project->getPatients}){
	my $c = $p->getCapture();
	
	 
	next if  $c->genome_version() =~ /HG38/;
	my $cid;
	if ($p->isgenome){
		 $cid =  $query->getCaptureId("genome_hg38");
		 	my $cmd2 ="del_calling_method.sh -project=$project_name  -method=haplotypecaller,haplotypecaller4 ";
		system($cmd2);
		$cmd2 ="add_calling_method.sh -project=$project_name  -method=dragen-calling,dragen-sv,dragen-cnv,wisecondor";
		system($cmd2);
		 
	}
	else {
	 my $capNameHG38=$c->name."_HG38";
	 $cid =  $query->getCaptureId($capNameHG38);
	unless  ($cid){
		die("no cpature HG38 for ".$capNameHG38);
	}
		my $cmd2 ="del_calling_method.sh -project=$project_name  -method=haplotypecaller,haplotypecaller4 ";
		system($cmd2);
		$cmd2 ="add_calling_method.sh -project=$project_name  -method=dragen-calling,dragen-sv,deepvariant,melt ";
		system($cmd2);
	
	
	}
	
	my $capInfo = $query->getCaptureInfos( $cid );
	update_table_patient($buffer->dbh,$capInfo->{capture_id},$p->id);
}

$buffer->dbh->commit();

my $dir_backup;

my $cmd1 ="change_genome_release.sh -project=$project_name -release=$release";

system($cmd1);


my $cmd2 ="del_calling_method.sh -project=$project_name  -method=haplotypecaller,haplotypecaller4 ";
 

system($cmd2);


$cmd2 ="add_calling_method.sh -project=$project_name  -method=dragen-calling,dragen-sv,deepvariant,melt ";
 
system($cmd2);

system("$RealBin/change_genome_release.pl -project=$pname -release=$release");

system("$RealBin/backup_hg19.pl -project=$pname -release=$release");
}
sub update_table_patient {
my ($dbh, $capture_id, $patient_id) = @_;
warn "UPDATE PolyprojectNGS.patient SET capture_id=$capture_id WHERE patient_id=$patient_id";
$dbh->do("UPDATE PolyprojectNGS.patient SET capture_id=$capture_id WHERE patient_id=$patient_id");
}