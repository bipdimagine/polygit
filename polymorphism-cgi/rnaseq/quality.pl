#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/kyoto";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
#use lib "/bip-d/soft/distrib/tabix/latest/perl";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/coverage";
use lib "$Bin/../packages/validation_variation"; 
use html; 
#use Set::;
use Carp;
use export_data;
use strict;
use Set::IntSpan::Fast::XS;
use Data::Dumper;
use GBuffer;
use Getopt::Long;
use Carp;
use Set::Intersection;
use Tabix;
use Storable qw/thaw/;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use coverage;
use Spreadsheet::WriteExcel;
use POSIX;
use validationQuery;
use Date::Tiny;
use JSON::XS;
use List::MoreUtils qw(firstidx uniq);
  use List::Util qw(sum);
#use Time::Piece;

my $w = "500px";

my @stats_name = ("total_read","assigned","ambiguity","nofeature","duplicate","unmapped");
my $cgi          = new CGI();
html::print_header_polydiag($cgi);

my $buffer1 = GBuffer->new();
my $list = $buffer1->listProjectsByAnalyse("rnaseq");


my $stats_projects;
foreach my $project_name (sort {$b cmp $a} @$list){
	my $buffer = GBuffer->new();

my $project = $buffer->newProject(-name=>$project_name);

my $file = $project->getCountingDir("featureCounts")."$project_name.count.genes.txt.summary";
next unless -e $file;
print $cgi->start_div({class=>"panel panel-primary",style=>"font-size: 11px;font-family:  Verdana;"});
print $cgi->div({class=>"panel-heading"},'<i class="fa fa-flag fa-lg" style="text-align: center;font-color: white;"></i> &nbsp'. $project_name." &nbsp ".$project->description."&nbsp".$project->getVersion);
	
print $cgi->start_table({class=>"table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 10px;font-family:  Verdana;", 'data-click-to-select'=>"true",'data-toggle'=>"table"});

my @lines = `cat $file`;
chomp(@lines);
my $header = shift(@lines);
my @files = split(" ",$header);

my $stats;

foreach my $patient (@{$project->getPatients}){
	my $bam = $patient->getBamFile();
	
	my $idx = firstidx {$_ eq $bam } @files;
	my %stat_patient;
	foreach my $line (@lines){
		my @t = split(" ",$line);
		
			$stat_patient{$t[0]} = $t[$idx] ;
		}
		
		$stats->{$patient->name} = \%stat_patient;
		
	}
	
my $type = {
	aligned => ["Assigned",""]
	
};
my $stat_project;
foreach my $n (keys %{$stats}){
	my $st = $stats->{$n};
	my $sum = sum values %$st;
	$stat_project->{global}->{total_read} += $sum;

	$stat_project->{sample}->{$n}->{total_read} = $sum;
	$_ = reverse $stat_project->{sample}->{$n}->{total_read};
	s/(\d{3})(?=\d)(?!\d*\.)/$1,/g;
	$stat_project->{sample}->{$n}->{total_read} = reverse $_;
	
	$stat_project->{sample}->{$n}->{unmapped} = sprintf("%.2f", ($st->{Unassigned_Unmapped} / $sum *100));
	$stat_project->{global}->{unmapped} += $stat_project->{sample}->{$n}->{unmapped};
	$stat_project->{sample}->{$n}->{duplicate} = sprintf("%.2f", ($st->{Unassigned_Duplicate} / $sum *100));
	$stat_project->{global}->{duplicate} += $stat_project->{sample}->{$n}->{duplicate};
	
	$stat_project->{sample}->{$n}->{assigned} = sprintf("%.2f", ($st->{Assigned} / $sum *100));
	$stat_project->{global}->{assigned} += 	$stat_project->{sample}->{$n}->{assigned};
		
	$stat_project->{sample}->{$n}->{ambiguity} = sprintf("%.2f", ($st->{Unassigned_Ambiguity} / $sum *100));	
	$stat_project->{global}->{ambiguity} += 	$stat_project->{sample}->{$n}->{ambiguity};
	
	$stat_project->{sample}->{$n}->{nofeature} = sprintf("%.2f", ($st->{Unassigned_NoFeatures} / $sum *100));	
	$stat_project->{global}->{nofeature} += 	$stat_project->{sample}->{$n}->{nofeature};
}


	
print qq{<thead>};
		print $cgi->start_Tr({class=>"success"});
my @header = ("patient",@stats_name);
print $cgi->th( \@header);
#print join("\t",@header)."\n";
	print $cgi->end_Tr();
print qq{</thead>};

foreach my $s (keys %{$stat_project->{sample}} ){
		print $cgi->start_Tr();
	my @line;
	push(@line,$s);
	foreach my $stat (@stats_name){
		push(@line,$stat_project->{sample}->{$s}->{$stat});
	}
	
	print $cgi->td(\@line);
	
		print $cgi->end_Tr();
	
}
#
my @line;
	push(@line,$project_name);

	my $nb = scalar(keys %{$stat_project->{sample}} );
#	warn $nb;
	foreach my $stat (@stats_name){
		$stat_project->{global}->{$stat} /= $nb;
			$stat_project->{global}->{$stat} = int ($stat_project->{global}->{$stat}*100) / 100;
		$_ = reverse 	$stat_project->{global}->{$stat};
	s/(\d{3})(?=\d)(?!\d*\.)/$1,/g;
		$stat_project->{global}->{$stat} = reverse $_;
		push(@line,$stat_project->{global}->{$stat} );
		#push(@line,sprintf("%.2f",$stat_project->{global}->{$stat}));
	}
#	print join("\t",@line)."\n";
print $cgi->start_Tr({class=>"warning"});
print $cgi->td(\@line);
	
		print $cgi->end_Tr();
		print $cgi->end_table;
print $cgi->end_div();
print $cgi->end_div();
	print "<br>";
		print "<br>";
} #end list;


	



