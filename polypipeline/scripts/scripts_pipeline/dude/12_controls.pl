#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";

#use Set::IntSpan;
use GenBoNoSqlLmdb;

use Carp;
use strict;
use Set::IntSpan::Fast::XS ;
use Data::Dumper;
use GBuffer;
use Getopt::Long;
use Carp;
 use JSON::XS;
 use List::MoreUtils qw(natatime uniq samples  );
  use List::Util qw(sum shuffle );
 use Tabix;
 
#use GenBoNoSql;
#use GenBoNoSqlDejaVu;
#use Array::Diff;
#use UnQLite;
 # use Tie::LevelDB; 
 # use Devel::Size qw(size total_size);
  #use Devel::Size::Report qw/report_size/;
  #use Statistics::Descriptive;
  

 $| =1;               
#ENSG00000234585
my $buffer = new GBuffer;
my $project_name;
my $fork;
GetOptions(
	'project=s' => \$project_name,
	'fork=s' => \$fork,
);
die("he no fork ") unless $fork;
my $jobs;
my $project = $buffer->newProjectCache( -name 			=> $project_name );
my $runs = $project->getRuns();

 my $sambamba = $buffer->software("sambamba");

#choose control by run .
my $controls;
my $max_controls = 50;
foreach my $r (@$runs) {
	#find_other_patient($r);
	#die();
	my $hps =  $r->getAllPatientsInfos();
	#warn Dumper($hps) if $r->id() eq 982;
	
	#only control
#	if (scalar(@$hps) < 8){
#		push(@{$controls->{$r->id}},@$hps);
#		next;
#	}
	my  @hps2 =  grep{$_->{project} ne $project_name && $_->{status} eq 1 && $_->{patient} =~ /GIAB/ } @$hps;
warn "coucou";
	if (@hps2>=$max_controls){
		@hps2 = samples $max_controls,@hps2;
		#push(@{$controls->{$r->id}},samples 6,@hps2);
	}
	
	if (@hps2<=$max_controls) {
		
		foreach my $h  (grep{$_->{project} ne $project_name && $_->{status} ne 1} @$hps){
			push(@hps2,$h);
			last if scalar(@hps2) >= $max_controls;
		}
			
		}
	
	if (@hps2<=$max_controls) {
		foreach my $h  (grep{$_->{project} eq $project_name && $_->{status} eq 1} @$hps){
			push(@hps2,$h);
			last if scalar(@hps2) >= $max_controls;
			
		}
	}
	
		if (@hps2<=$max_controls) {
		foreach my $h  (grep{$_->{project} eq $project_name && $_->{status} ne 1} @$hps){
			push(@hps2,$h);
			last if scalar(@hps2) >= $max_controls;
			
		}
	}
	#push(@{$controls->{$r->id}},@hps2);

	
	warn "\t*******************".scalar (@hps2);
	#die();
	 $max_controls = 10 if $project->isExome;
	if (scalar (@hps2) < $max_controls ) {
		warn scalar (@hps2);
		find_other_patient($r,\@hps2);
	}
	warn "\t*******************".scalar (@hps2);
	my %contr_projects;
	map {$contr_projects{$_->{project}} ++} @hps2;
warn "\t*******************".scalar (@hps2);
	#warn Dumper @hps2;


	foreach my $pr (keys %contr_projects){
		my $buffer1 = new GBuffer;
		my $project1 = $buffer1->newProjectCache( -name 			=> $pr );
		#next() if $project1->name() eq "NGS2018_2286";
		foreach my $p (grep{$_->{project} eq $pr} @hps2){	
		
			my $patient = $project1->getPatient($p->{patient});
			my $b = $patient->getBamFileName();
			next unless -e $b;		
			my $scompute = &compute_sex($patient);
			my $cov_file = $patient->fileNoSqlDepth;
			if ($project1->name eq $project_name && !( -e $cov_file) ){
				warn "coverage"; 
				warn "$Bin/../coverage_genome.pl -project=$project_name -fork=$fork -patient=".$patient->name;
				system("$Bin/../coverage_genome.pl -project=$project_name -fork=$fork -patient=".$patient->name);
				die($cov_file) unless -e $cov_file;
			}
			warn $patient->name." ".$pr  unless -e $cov_file;
			next unless -e $cov_file;
			$p->{file_depth} = $cov_file;
			$p->{dir_depth} = $project1->getCoverageDir()."/lmdb_depth";
			$p->{bam} = $patient->getBamFile();
			$p->{compute_sex} = compute_sex($patient);
			$p->{bd_sex} = $patient->sex;
			if($p->{compute_sex} ne $p->{bd_sex}){
				#warn ("--> pb sexe patient : " .$patient->name."\n compute sex : ".$p->{compute_sex}."\n bd sex : ".$p->{bd_sex}) ;
				warn ("pb sexe patient : " .$patient->name."\n compute sex : ".$p->{compute_sex}."\n bd sex : ".$p->{bd_sex});# if $patient->project->name eq $project_name ;
				$p->{bd_sex} = $p->{compute_sex};
#				die();
				#next;
			}
			push(@{$controls->{$r->id}},$p) ;
		} 
		
	}
}

#my $chr = $project->getChromosome("Y");
foreach my $r (keys %$controls){
	foreach my $c (@{$controls->{$r}}){
		next unless $r ==1453;
		warn "run : ".$r. " c : ".Dumper $c;
	}
}


foreach my $r (@$runs){
	warn $r->id." ".$r->name." ". scalar(@{$controls->{$r->id}});
	#die("not enough control ". $r->id." ".$r->name) if scalar(@{$controls->{$r->id}}) <= 2;
}
my $no =  $project->noSqlCnvs("c");
warn Dumper  $controls;
$no->put("raw_data","controls",$controls);

warn "OK END";
sub compute_sex {
	my $patient = shift;
	#return 1;
	my $bam = $patient->getBamFile();
	my $sambamba = $buffer->software("sambamba");
	my $chry = $project->getChromosome("Y")->fasta_name;
	my @vv =`$sambamba depth region $bam -L $chry:2654896-2655740 2>/dev/null | cut -f 5`; #`tabix $f mean_chrY:1-15| cut -f 3`;
	my $v = $vv[-1];
	my $sex;
	if ($v < 6){
		$sex = 2;
	}
	else{
		$sex =1;
	}
	return $sex;
}
	
my %patients_captures;

sub find_other_patient {
	my ($run,$controls) = @_;
	 my $capture = $run->project->getCaptures()->[0];
	 if (exists $patients_captures{$capture->name}){
	 	 push(@$controls, @{$patients_captures{$capture->name}});
	 	 return; 
	 	
	 }
	 my $mean_norm;
	 my @apos ;
	 my $file_bed = $run->getCapture->gzFileName();
	 open(SCHUF," zcat $file_bed | shuf -n 50 |");
	 while(<SCHUF>){
	 	chomp();
	 	my ($a,$b,$c) = split(" ",$_);
	 	my $chr = $project->getChromosome($a);
	 	push(@apos,{chr=>$chr->name,start=>$b-50,end=>$c+50});
	 	
	 }
	 close(SCHUF);
	 my $nv;
	 
	 foreach my $c (@$controls){
	 		my $buffer2 = GBuffer->new();
			my $project2 =  $buffer2->newProject( -name 			=> $c->{project});
			my $patient = $project2->getPatient($c->{patient});
			my $hp = $patient->nb_reads;
			warn Dumper $hp;
			foreach my $ps (@apos){
				warn Dumper $ps;
				warn $hp->{$ps->{chr}};
				$mean_norm += ($patient->maxDepth($ps->{chr},$ps->{start},$ps->{end})/$hp->{$ps->{chr}});
				$nv ++;
			}
	 }
	 $mean_norm /= $nv;
	 
	 my $query = $project->buffer->getQuery->listAllProjectsNameByCaptureId($capture->id());
	my $x;
	my $res;
	my $limit = 10 - scalar(@$controls);
	
	 foreach my $project_name2 (@$query){
	 		next if $project_name2 =~ /NGS2010/;
	 		my $buffer2 = GBuffer->new();
	 	
			my $project2 =  $buffer2->newProject( -name 			=> $project_name2);
			my $nbx = 0;
			foreach my $p (@{$project2->getPatients}){
				my $capture2  = $p->getCapture();
				next if $p->status == 1;
				my $bam; 
				eval {
				
				next if $capture->name ne $capture2->name;
				$bam =  $p->getBamFileName();
				
				};
				next unless -e $bam;
				next unless -e $p->NoSqlDepthDir()."/".$p->name . ".depth.lmdb";
			#	warn Dumper $p->nb_reads();
				$nbx++;
				my $hp;
				$hp->{family} = $p->getFamily()->name;
				$hp->{id} = $p->id;
				$hp->{patient} = $p->name;
				$hp->{status} = $p->status;
				$hp->{sex} = $p->sex;
				$hp->{father} = "";
				$hp->{mother} = "";
				my $nb_reads = $p->nb_reads();
				$hp->{norm} =  $nb_reads->{norm};
				$hp->{norm_c} =  abs($mean_norm - $nb_reads->{norm});
				$hp->{project} = $project2->name;#->getFamily()->getMother->name;
				 my $nv1;
				 my $mean_norm1 =0;
				 my $hpp = $p->nb_reads;
				foreach my $ps (@apos){
					
					$mean_norm1 += ($p->maxDepth($ps->{chr},$ps->{start},$ps->{end})/$hpp->{$ps->{chr}});
					$nv1 ++;
				}
				 $mean_norm1 /= $nv1;
				 $hp->{mean_max} = $mean_norm1;
				 $hp->{mean_sort} = abs($mean_norm1-$mean_norm);
				push(@$res,$hp);
				last if $nbx > 3;
		
			}
			last if scalar(@$res) > 50;
	 }
	  @$res = sort {$a->{mean_sort} <=> $b->{mean_sort}} (@$res);
	 my @toto = splice (@$res,0,$limit);
	# warn Dumper(@toto);
	 # die();
#	warn Dumper $res;
	#die();
	$patients_captures{$capture->name} = \@toto;
	 push(@$controls, @toto ); 
	
}
	