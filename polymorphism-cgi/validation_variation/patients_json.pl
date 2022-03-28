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
my $buffer = GBuffer->new();
my $cgi    = new CGI();
my $project_name = $cgi->param('project');
my $project = $buffer->newProject(-name=>$project_name);
my $server =  $ENV{HTTP_HOST};
$server = "darwin.bipd.fr" if $server eq "bipd";
$server = "darwin.bipd.fr" if $server =~/10\.200\.27/;

my @out;

foreach my $p (sort {$a->name cmp $b->name() }@{$project->getPatients()}){
	#next unless $p->isChild();
	my $item;
	my $test;
	if ($project->isDiagnostic){
	my $capture = $p->getCapture();
	my $vquery = validationQuery->new(dbh=>$buffer->dbh,capture_name=>$project->validation_db());
	my $variations_todo = $vquery->get_variations_todo(project_name=>$project_name,sample_name=>$p->{name});
	my $variations_ions = $vquery->get_variations_ions(project_name=>$project_name,sample_name=>$p->{name});
	my $variations_sanger = $vquery->get_variations_sanger(project_name=>$project_name,sample_name=>$p->{name});
	my $exons_todo = $vquery->get_exons(project_name=>$project_name,sample_name=>$p->{name});
	my $item;
	if ($variations_ions ){
		$item->{ion} = scalar(grep {/chr/} keys %$variations_ions);
	}
	if ($variations_todo ){
		
		$item->{todo} = scalar(grep {/chr/}keys %$variations_todo);
	}
	if ($variations_sanger ){
		$item->{sanger} = scalar(grep {/chr/} keys %$variations_sanger);
	}
	 $test = $vquery->get_report(project=>$project->name,sample=>$p->name);
	}
	
	$item->{label} = $p->name();
	
	$item->{save} = "-" ;
	$item->{save} = "OK"  if $test;
	my $res2 = "- - -1";
	eval {
			my $tabix = $p->tabix_coverage();
			my $res = $tabix->query("mean_chrY",14,16);
			 $res2 = $tabix->read($res);
	};
			my($a,$b,$val) = split(" ",$res2);
			my $sex ='F';
			$sex = 'M' if $val>50;
			$sex = 'NC' if $val<0;
	$item->{sex} = "M";
	$item->{sex} = "F" if $p->sex == 2;	
	$item->{child} = 0;	
	$item->{child} = 1 if $p->isChild;
	#$p->{child}= 1 unless $project->isFamilial;
	$item->{compute_sex} = $p->compute_sex;
	$item->{id} = $p->name();
	$item->{machine} = $p->getRun()->machine();
	$item->{date} = $p->getRun()->date();
	$item->{active} = -1;
	my $bam_file = $p->getBamFiles->[0];
	$bam_file =~ s/\/\//\//g;
	my @lTmp = split('/', $bam_file);
	my $max = scalar(@lTmp);
	my $bam_url = "http://$server/NGS/".$project->name()."/".$project->getVersion()."/align/".$lTmp[$max-2].'/'.$lTmp[$max-1];
	$item->{bam} = $bam_url;
		push(@out,$item);
	}

export_data::print_simpleJson($cgi,\@out);
exit(0);
