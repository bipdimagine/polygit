#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/kyoto";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use Carp;
use export_data;
use strict;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use JSON::XS;
use validationQuery;
use Digest::MD5::File qw(dir_md5_hex file_md5_hex url_md5_hex file_md5);
use Capture::Tiny ':all';
my $buffer = GBuffer->new();
my $cgi    = new CGI();
my $project_name = $cgi->param('project');
my $project = $buffer->newProjectCache(-name=>$project_name);
my $server =  $ENV{HTTP_HOST};
$server = "darwin.bipd.fr" if $server eq "bipd";
$server = "darwin.bipd.fr" if $server =~/10\.200\.27/;
my $user = $cgi->param('user_name');

my $is_genome_project = 0;
my $is_diag_project = 0;
if ($project->isGenome()) { $is_genome_project = 1; }
elsif (not $project->isExome()) { $is_diag_project = 1; }

my $sid = join(";",map{$_->id} sort {$a->name cmp $b->name() } @{$project->getPatients});
$sid .= ";".file_md5_hex($Bin."/fast_patients_json.pl");	
my 	$no_cache = $project->get_lmdb_cache_summary("r");
my  $out = $no_cache->get_cache($sid); 
if($out){
	print $out;
	exit(0);
} 
 $no_cache->close();
# $out = undef;
unless ($out){
foreach my $p (sort {  $a->name cmp $b->name  or  $a->getFamily->name cmp $b->getFamily->name or $b->isChild <=> $a->isChild }@{$project->getPatients()}){
	#next unless $p->isChild();
	my $item;
	my $test;

	
	$item->{label} = $p->name();
	$item->{is_genome} = $is_genome_project;
	$item->{is_diag} = $is_diag_project;
	
	$item->{sex} = "M";
	$item->{sex} = "F" if $p->sex == 2;	
	#$item->{sex} = $p->return_icon;	
	
	
	$item->{child} = 0;	
	$item->{child} = 1 if $p->isChild;
	
	my  $class ="icon_babyboy_d";
	if($p->isChild &&  $p->sex == 1 && $p->status == 1) {
		$class ="icon_babyboy_s";
	}
	elsif($p->isChild &&  $p->sex == 1 && $p->status == 2) {
		$class ="icon_babyboy_d";
	}
	elsif($p->isChild &&  $p->sex == 2 && $p->status == 1) {
		$class ="icon_babygirl_s";
	}
	elsif($p->isChild &&  $p->sex == 2 && $p->status == 2) {
		$class ="icon_babygirl_d";
	}
	elsif(!($p->isChild) &&  $p->sex == 1 && $p->status == 1) {
		$class ="icon_male_s";
	}
	elsif(!($p->isChild) &&  $p->sex == 1 && $p->status == 2) {
		$class ="icon_male_d";
	}
	elsif(!($p->isChild) &&  $p->sex == 2 && $p->status == 1) {
		$class ="icon_female_s";
	}
	elsif(!($p->isChild) &&  $p->sex == 2 && $p->status == 2) {
		$class ="icon_female_d";
	}
	$item->{icon} =  $class;
	$item->{status} =  $p->status;
	$item->{id} = $p->name();
	my $vstatus = $p->getLatestValidationStatus($user);
	my $hstatus = $p->getLatestValidationStatus($user);
	
	$item->{validation_term} = "-";
	$item->{validation_date} = "-";
	$item->{validation_user} = "-";
	if  ($hstatus){
	my $xx;
	my $date2;
	($date2,$xx) = split(" ",$hstatus->{modification_date});
	my $u = $hstatus->{user_name}; 
	$item->{validation_term} = $hstatus->{term};
	$item->{validation_date} = $date2;
	$item->{validation_user} = $u;
	}
	$item->{validation_status} = "0";
	if($vstatus){
		$item->{validation_status} = $vstatus->{status};	
	}
	my $fam = $p->getFamily();
	$item->{fam} = "-";
	$item->{fam} = $p->getFamily()->name if $fam &&  $p->getFamily()->name ne $p->name or scalar(@{ $p->getFamily()->getParents})>0 ;
	
	$item->{hasJonction} = 0;
	$item->{hasJonction} = 1 if $p->hasJunctions(0);
	
	$item->{bam} = $p->bamUrl();
	
	push(@$out,$item);
	}
	
}

my $stdout2 = tee_stdout {
export_data::print_simpleJson($cgi,$out);
};
 eval {
	my $no_cache = $project->get_lmdb_cache_summary("w");
   	$no_cache->put_cache($sid,$stdout2,2400); 
 	$no_cache->close();
 };
exit(0);
