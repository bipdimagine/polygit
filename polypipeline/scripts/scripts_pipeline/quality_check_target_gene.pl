#!/usr/bin/perl
use CGI qw/:standard :html3/;
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use colored;
use lib $Bin;
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use Set::IntSpan::Fast::XS;
use List::MoreUtils qw(part);
use List::Util qw(sum);
use Text::Table;

my $projectName;
my $end_ext = "uni";
my $details;
my $fork = 1;
my $vcf_file;
my $log_file;
my $noclean;
my $chr_names;
GetOptions(
	'project=s'		=> \$projectName,
);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
my $samtools = $buffer->software("samtools");
colored::stabilo('blue',"CHECK FILE  : ");
my $patients = $project->getPatients();

my $bam_stats;
my $size;
my $errors;

foreach my $p (@$patients){
	my $bam = $p->getBamFile();
	push(@{$errors->{file}},"NO BAM FOR : ".$p->name()) unless -e $bam;
	my $vcf_indels = $p->getIndelsFiles();
	foreach my $vcf_indel (@$vcf_indels){
	push(@{$errors->{file}},"NO VCF INDEL FOR : ".$p->name()) unless -e $vcf_indel;
	}
	
	my $vcf_vars = $p->getVariationsFiles();
	foreach my $vcf_var (@$vcf_vars){
	push(@{$errors->{file}},"NO VCF VARIATION FOR : ".$p->name()) unless -e $vcf_var;
	}
}

if (exists $errors->{file}){
	foreach my $e (@{$errors->{bam}}){
		colored::stabilo('red',$e);
	}
	die();
	
}
colored::stabilo('green'," YOUR FILE SEEMS TO BE OK ");

colored::stabilo('blue',"CHECK STATISTICS BAM  : ");
foreach my $p (@$patients){
	my $bam = $p->getBamFile();
	push(@{$errors->{bam}},"NO BAM FOR : ".$p->name()) unless -e $bam;
	my $fsize = -s  $p->getBamFile();
	my $cmd = qq{$samtools idxstats $bam};
	my @t = `$cmd`;
	chomp(@t);
	foreach my $l (@t){
	
		my($chr,$size,$amount) = split(" ",$l);
			next if $chr eq "*";
		push(@{$bam_stats->{$chr}->{data}},$amount);
		
	}
}

my $table;


my @header =("patient");
foreach my $chr (sort {$project->getChromosome($a)->karyotypeId <=> $project->getChromosome($b)->karyotypeId} keys %$bam_stats){
	next if $chr eq "*";
	push(@header,$chr);
	my $t = sum @{$bam_stats->{$chr}->{data}};
	my $mean = int($t / scalar @{$bam_stats->{$chr}->{data}});
	$bam_stats->{$chr}->{mean} = $mean;
}
my $tb = Text::Table->new(
      @header
    );
    $table->{stat_chromsome}->{header} = \@header;
my @lines;
for (my $i=0;$i<@$patients;$i++){
	my @line;
	push(@line,$patients->[$i]->name);
	
	foreach my $chr (sort {$project->getChromosome($a)->karyotypeId <=> $project->getChromosome($b)->karyotypeId}  keys %$bam_stats){
	#push(@line,$bam_stats->{$chr}->{mean});
	my $v = $bam_stats->{$chr}->{data}->[$i];
	my $m =  $bam_stats->{$chr}->{mean};
	
	if ($v eq 0){
			push(@line,colored::stabilo("red","NONE",1));
	}
	elsif ($v < abs($m*0.3)){
		push(@line,colored::stabilo("cyan","BAD",1));
	} 
	elsif ($v < abs($m*0.5)){
		push(@line,colored::stabilo("orange","POOR",1));
	} 
	else {
		push(@line,colored::stabilo("green","OK",1));
	}
	
	}
		push(@lines,\@line);
	}

	$tb->load(@lines);
	print $tb;
	colored::stabilo('blue',"---------VARIATIONS ------------");
	
	my @transcripts_cgi = @{$project->bundle_transcripts() } ;
	
	my @header2 =("patient");
	my $transcripts;
	foreach my $t (@transcripts_cgi){
		 my $tr1 = $project->newTranscript($t);
		 push(@$transcripts,$tr1);
		 push(@header2,$tr1->getGene->external_name);
	}
	my @lines2;
	foreach my $p (@$patients){
		my @line;
		my $vs = $p->getStructuralVariations();
		push(@line,$p->name);
		foreach my $t (@$transcripts){
				my $n =0;
				foreach my $v (@$vs){
					$n++ if $v->start>$t->start && $v->end < $t->end;
				}
				push(@line,colored::stabilo("green",$n,1)) if $n > 0;
					push(@line,colored::stabilo("red",$n,1)) if $n <= 0;		
		}
		push(@lines2,\@line);
	}
my $tb2 = Text::Table->new(
      @header2
    );
    
    $tb2->load(@lines2);
	print $tb2;
	colored::stabilo('blue',"----------COVERAGE-----------");
	
	my @lines2;
	foreach my $p (@$patients){
		my @line;
		push(@line,$p->name);
		foreach my $t (@$transcripts){
				my $n = int($t->mean_coding_coverage($p));
			
				if ($n <10){
					push(@line,colored::stabilo("red",$n,1));
				}
				if ($n <50){
					push(@line,colored::stabilo("orange",$n,1));
				}
				else{
					push(@line,colored::stabilo("green",$n,1));
				}
				
		}
		push(@lines2,\@line);
	}
my $tb2 = Text::Table->new(
      @header2
    );
    
    $tb2->load(@lines2);
	print $tb2;
my $stat;

foreach my $p (@{$project->getPatients()}){
	my $name = $p->name();
foreach my $v (@{$p->getVariations()}){
	next if $v->frequency > 0.1;
	$stat->{$name}->{all} ++;
	foreach my $p1 (@{$v->getPatients}){
		$stat->{$name}->{$p1->name} ++;
	}
}
}	
	colored::stabilo('blue',"----------% IDENTITY -----------");
my @lines3;
my @header3;
push(@header3,"sample");
foreach my $p (@{$project->getPatients()}){
	my $name = $p->name;
	push(@header3,$name);
	my $nb = $stat->{$name}->{$name};
	my @lines;
	push(@lines,$name);
	foreach my $p1 (@{$project->getPatients()}){
		my $name1 = $p1->name();
		
		my $nb1 = $stat->{$name}->{$name1};
		my $p= "-";
		if ($name ne $name1 && $nb >0 ){
			$p = int(($nb1/$nb)*100);
		}
		push(@lines,colored::stabilo("red",$p,1))	if ($p>90);
		push(@lines,colored::stabilo("green",$p,1))	if ($p<=90);
	}
	push(@lines3,\@lines);
}

my $tb3 = Text::Table->new(
      @header3
    );
        $tb3->load(@lines3);
        print $tb3;
#	colored::stabilo('blue',"---------- CONTROLE PARSING -----------");
#	my $bcftools = $buffer->software("bcftools");
#	foreach my $p (@{$project->getPatients()}){
#		my $files = $p->getVariationsFiles();
#		push(@$files,@{$p->getIndelsFiles()});
#		
#		my $nbs =0;
#		my $nbi =0;
#			foreach my $f (@$files){
#				warn $f;
#			my @res = `zgrep -v "0\/0" $f | $bcftools stats - `;
#			
#				foreach my $l (@res){
#						chomp($l);
#					next unless $l =~ /^SN/;
#					my @t = split(" ",$l);
#				
#					
#					if ( $l =~/SNPs|MNPs/){
#						$nbs+= $t[-1];
#					}
#					if ( $l =~/indels/){
#						$nbi+= $t[-1];
#					}
#				#	next unless $l =~/SNPs/;
#					next unless $l =~/SNPs|MNPs|indels|others/;
#										
#					
#				}
#			}
#			my $vs = $p->getVariations();
#			my $vi = $p->getIndels();
#			warn scalar(@$vs)." ".$nbs."   ".scalar(@$vi)." ".$nbi;
#			
#		
#	}
	


