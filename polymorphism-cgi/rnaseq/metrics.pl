#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo";
use lib "$Bin/../../GenBo/lib/GenBoDB";
use lib "$Bin/../../GenBo/lib/obj-nodb";
use lib "$Bin/../../GenBo/lib/kyoto";
use lib "$Bin/../../GenBo/lib/obj-nodb/packages";
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
my $stats = {
		"PCT_USABLE_BASES" => "Usable",
		"PCT_CODING_BASES" => "Coding ",
		"PCT_UTR_BASES" => "UTR ",
		"PCT_INTRONIC_BASES" => "Intronic",
		"PCT_INTERGENIC_BASES" => "Intergenic",
		"PCT_MRNA_BASES" => "Mrna",
		"PF_BASES" => "Bases",
		"PF_ALIGNED_BASES"=>"Aligned",
		"PCT_RIBOSOMAL_BASES" => "ribosomal",
};

my @sheaders = (
		"PF_BASES" ,
		"PF_ALIGNED_BASES",
		"PCT_USABLE_BASES",
		"PCT_MRNA_BASES",
		"PCT_CODING_BASES" ,
		"PCT_UTR_BASES" ,
		"PCT_INTRONIC_BASES" ,
		"PCT_INTERGENIC_BASES" ,
		"PCT_RIBOSOMAL_BASES"
	
		);

html::print_header_polydiag($cgi);

my $buffer1 = GBuffer->new();
my $list = $buffer1->listProjectsByAnalyse("rnaseq");
foreach my $l (@$list){
	warn $l;
}


my $stats_projects;
foreach my $project_name (sort {$b cmp $a} @$list){
	my $buffer = GBuffer->new();

my $project = $buffer->newProject(-name=>$project_name);
my $file = $project->getCountingDir("featureCounts")."$project_name.count.exons.txt.summary";
next unless -e $file;

	my $dir_out= $project->getCountingDir("featureCounts")."/metrics";
	system("mkdir $dir_out && chmod a+rwx $dir_out") unless -e $dir_out;
	#my $fileout  =$dir_out."/$name.metrics";
my $dir_out = $project->getCountingDir("featureCounts")."/metrics/";
next unless -e $dir_out;


my $global;
my $alert = {};
foreach my $patient (@{$project->getPatients}){
	my $name = $patient->name();
	my $fileout  =$dir_out."/$name.metrics";
	open (OUT,$fileout);
	my @title;
	my $res;
	while(<OUT>){
		next if $_ =~/^#/;
		chomp();
		warn $_   if $name eq "aa2";
		if ($_ =~ /PF_BASES/){
			
			@title = split(" ",$_);
			
			next;
		}
		
		if (@title){
			my @temp = split("\t",$_);
			my $zz; 
			my $xx;
			for (my$i=0;$i<@title;$i++){
				next unless exists $stats->{$title[$i]};
				if ($title[$i] =~/PCT/){
				$global->{$name}->{$title[$i]} = sprintf("%.1f",($temp[$i]*100))  ;
				$zz->{$title[$i]} =$temp[$i] ;
				if ($title[$i]  eq "PCT_USABLE_BASES" && $temp[$i] <0.10) {
					$alert->{$name} ="danger";
				}
				elsif ($title[$i]  eq "PCT_USABLE_BASES" && $temp[$i] <0.2) {
					$alert->{$name} = "warning";
				}
				}
				#elsif ($title[$i] =~/BASE/){
					#$xx->{$title[$i]} =$temp[$i] ;
				#}
				else {
					if ($title[$i]  eq "PF_ALIGNED_BASES"){
						$global->{$name}->{$title[$i]} = sprintf("%.1f",($temp[$i] / $global->{$name}->{PF_BASES})*100);
					}
					else {
						$global->{$name}->{$title[$i]} =$temp[$i] ;
					}
					
				}
					
			}
			#warn Dumper $zz if $name eq "aa2";
			#warn Dumper $xx if $name eq "aa2";
			#warn Dumper @title if $name eq "aa2";
			#warn Dumper @temp if $name eq "aa2";
			#warn $fileout if $name eq "aa2";
				last;			 
			}
		}
		 close(OUT);
}	
 	print $cgi->start_div({class=>"panel panel-primary",style=>"font-size: 11px;font-family:  Verdana;"}) ;
my $danger =  grep {$_ eq "danger"} values %$alert;
 my $warning =  grep {$_ eq "warning"} values %$alert;
 if ($danger){

 	print $cgi->div({class=>"panel-heading"},'<i class="fa fa-exclamation-triangle fa-lg" style="text-align: center;font-color: red;"></i> &nbsp'. $project_name." &nbsp ".$project->description."&nbsp".$project->getVersion);
 }
 elsif ($warning){
 	print $cgi->div({class=>"panel-heading"},'<i class="fa fa-bell fa-lg" style="text-align: center;font-color: white;"></i> &nbsp'. $project_name." &nbsp ".$project->description."&nbsp".$project->getVersion);
 }
 else {
 		print $cgi->div({class=>"panel-heading"},'<i class="fa fa-flag fa-lg" style="text-align: center;font-color: white;"></i> &nbsp'. $project_name." &nbsp ".$project->description."&nbsp".$project->getVersion);
 }
 

	
print $cgi->start_table({class=>"table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 10px;font-family:  Verdana;", 'data-click-to-select'=>"true",'data-toggle'=>"table"});

print $cgi->start_Tr({class=>"success"});
			my @header = ("patient",map{$stats->{$_}}@sheaders );
			print $cgi->th( \@header);
			print $cgi->end_Tr();
			print qq{</thead>};
my $sum;			
foreach my $patient (sort{$a->name cmp $b->name}  @{$project->getPatients} ){
		my $name = $patient->name();
		my @line;
		if (exists $alert->{$name}){
			print $cgi->start_Tr({class=>$alert->{$name}});
		}
		else {
			print $cgi->start_Tr();
		}	
			
			push(@line,$name);
			my $i=0;
			foreach my $stat (@sheaders){
				push(@line,$global->{$name}->{$stat});
				$sum->[$i] += $global->{$name}->{$stat};
				$i++;
				
			}
	
				print $cgi->td(\@line);
		
	
		print $cgi->end_Tr();
		}
print $cgi->start_Tr({class=>"info"});
my @line ="project";
my $nb_patient = scalar(@{$project->getPatients});

foreach my $v (@$sum){
	push(@line,int($v/$nb_patient));
}
print $cgi->td(\@line);		
print $cgi->end_Tr();
print $cgi->end_table;
	
print $cgi->end_div();
print $cgi->end_div();
	print "<br>";
		print "<br>";		
			
}

exit(0);

