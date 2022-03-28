#!/usr/bin/perl
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
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
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use List::MoreUtils qw(part);
use File::Temp;
use JSON::XS;
use check_utils; 
use Statistics::Zscore; 
use Math::Combinatorics ;

my $projectName;
my $end_ext = "uni";
my $details;
my $fork = 1;
my $vcf_file;
my $log_file;
my $noclean;
my $chr_names;
my $patient_name;
my $run;
GetOptions(
	'project=s'		=> \$projectName,
	"vcf_file=s" => \$vcf_file,
	"patients=s" =>\$patient_name,
	"run=s" =>\$run,
);

 if ($run eq 'all' || $run eq ''){
		$run ="duo,trio" ;
 }

#die("run=all,duo,trio") unless $run ~= /duotrioall/;

my $date = `date`;

chomp($date);
if ($log_file){
	open (STDOUT,">>".$log_file);
}

colored::stabilo('blue',"START QUALITY CHECK ");
die("\n\nERROR: -project option mandatory. Exit...\n\n") unless ($projectName);

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $projectName );
$patient_name ="all" unless $patient_name;
my $patients = $project->get_list_patients($patient_name,",");



my $gatk  = $buffer->getSoftware('java')." -jar ".$buffer->getSoftware('gatk');
my $plink = $buffer->getSoftware('plink');


foreach my $p (@$patients){
	warn $p->name()."\n";
}

my @chrs_plink = (1..22,'X','Y');
if ($chr_names eq "all"){
	
	 @chrs_plink = (1..22,'X','Y');
}

my @lChr = ('21', 'X', 'Y');

fast_plink3($vcf_file,$patients) if $run =~/trio/;
find_duo($vcf_file,$patients) if $run =~ /duo/; 

exit(0);

#warn $dir_out if $noclean; 



sub fast_plink2 {
	my ($vcf_file,$patients) = @_;

	my $dir_out = $project->getCallingPipelineDir("unifiedgenotyper").'/plink';
	

	
	my $dir = $project->getCallingPipelineDir("unifiedgenotyper");
	
	my $logPlink = $dir.'/plink/'.$projectName.'.plink.resume';
	
	my @snps;

	open (VCF, "<$vcf_file");
	my @samples;
	my $DP;
	while (my $line = <VCF>){
		next if $line =~ /##/;
		chomp($line);
		if ($line =~ /#CHR/){
			my @tdata = split(" ",$line);
			@samples = splice(@tdata,9);#$tdata[9..-1];
			next;
		}
		
		my ($chrom,$pos,$pid,$ref,$alt,$qual,$filter,$info,$format,@gsamples) = split(" ",$line); 
		next if $ref =~ /,/;
		next if $alt =~ /,/;
		my $snp;
		my $id = $chrom."_".$pos."_".$alt;
		$snp->{id} = $id;
		$snp->{position} = $pos;
		$chrom =~ s/chr//;
		$chrom = "MT" if $chrom eq "M";
		$snp->{chr} = $chrom;
		unless ($DP){
			my @t = split(":",$format);
			for (my $i =0;$i<@t;$i++){
				$DP=$i if $t[$i] eq "DP";
			}
			$DP =-1 unless $DP;
		}
		for (my $i=0;$i<@samples;$i++){
			my $name = $samples[$i];
			my $string = $gsamples[$i];
			my $geno;
			if ($string =~ /0\/0/) { $geno ="$ref $ref"; }
			elsif ($string =~ /0\/1/) { $geno ="$ref $alt"; }
			elsif ($string =~ /1\/0/) { $geno ="$ref $alt"; }
			elsif ($string =~ /1\/1/) { $geno ="$alt $alt"; }
			elsif ($string =~ /.\/./) { $geno ="0 0"; }
			elsif ($string =~ /./) { $geno ="0 0"; }
			else { die($string); }
			if ($DP > 0){
				my @tt = split(":",$string);
				my $v = $tt[$DP];
				$geno = "0 0" if $v < 10;
								
			}
			$snp->{samples}->{$name} = $geno;
		}
		
	
		#my ($find) = grep {$snp->{samples}->{$_}  eq "0 0"} @samples;
		#next if $find;
		push(@snps,$snp);
	}
	
	$dir .= "/plink_trio/";
	mkdir $dir unless -e $dir;
	my @samples_name;
	foreach my $p (sort {$a->name cmp $b->name} @{$patients}){
			push(@samples_name,$p->name());
	}
	

	 my $combinat = Math::Combinatorics->new(count => 3,
                                          data => [@samples_name],
                                         );
    my $res;      
    	my %dejavu;  
    	my @comb = combine(3,@samples_name) ;   
    	my $nb =  scalar(@comb);   
    	my $cpt =0;
    	foreach my $c (@comb){                     
    		$cpt ++;
    		warn "$cpt/$nb\n"; 
#	 while(my @combo = $combinat->next_combination){
	 		my @toto = permute(@$c);
	 		my $dp;
	 	
    			foreach my $fam (@toto){
    				my $child = $fam->[-1];
    				next if exists $dejavu{$child};
    				next if exists $dejavu{$fam->[0]};
    				next if exists $dejavu{$fam->[1]};
    				my $s = join("-",sort{$a cmp $b} ($fam->[0],$fam->[1]));
    				#warn $s;
    				next if exists $dp->{$s};
    				$dp->{$s} ++;
    				$res->{$child}->{score} = 999 unless exists $res->{$child};
    				my $tped_file = $dir."/".$project->name.".1.tped";
				my $ped_file = $dir."/".$project->name.".ped";
				open (PED,">".$ped_file);
				print PED join("\t",("F",$fam->[0],0,0,1,1))."\n";
				print PED join("\t",("F",$fam->[1],0,0,2,1))."\n";
				print PED join("\t",("F",$fam->[2],$fam->[0],$fam->[1],2,1))."\n";
				close PED;
				
			open (TPED,">".$tped_file);
			my $nb_snp =0;
		
			foreach my $snp (@snps){
					my @data;
					push(@data,$snp->{chr});
					push(@data,$snp->{id});
					push(@data,0);
					push(@data,$snp->{position});
				 my $found =0;
				 my %tt;
				 foreach my $name  (@$fam){
				 	$found ++ if (exists $snp->{samples}->{$name}  && $snp->{samples}->{$name} ne "0 0") ;
				 	$tt{$snp->{samples}->{$name}} ++;
				 	
				 }
				 next if $found <3;
				 next if scalar(keys %tt) == 1;
				 
				 $nb_snp ++;
				 last if $nb_snp > 10000;
				foreach my $name  (@$fam){
			
					if (exists $snp->{samples}->{$name} ) { push(@data,$snp->{samples}->{$name}); }
						else { warn "coucou";push(@data,$snp->{samples}->{$name}, "0 0");	}
					}
						print TPED join("\t",@data)."\n";
					}
					close TPED;
					my $cmd2 = "$plink --tped $tped_file --tfam $ped_file --noweb --mendel  --out $dir/$projectName.1";		
					my @log = `$cmd2`;
					open(MENDEL,"$dir/$projectName.1.imendel") || die("problem with $dir/$projectName.imendel \n$cmd2");
					while(my $line =<MENDEL> ){	
						chomp($line);
						next if $line =~ /FID/;
						my ($f, $name,$nb) = split(" ",$line);
						next if $name ne $child;
						#warn  $nb_snp;
						next if $nb_snp == 0;
						my $p = int($nb/$nb_snp*10000)/100;
						#warn $p."\n";
						if ($p < $res->{$child}->{score}){
							$res->{$child}->{score} = $p;
							$res->{$child}->{fam} = $fam;
							warn join("->",@{$fam})."\n"  if $p < 2;
							$dejavu{$child} ++ if $p < 1;
						}
					
			}
					close(MENDEL);
    				
    			}
    			
    			
    			
  }
  
  my %f_dejavu;
  foreach my $k (sort {$res->{$a}->{score} <=> $res->{$b}->{score} } keys %$res){
  	next if $res->{$k}->{score} > 3;
  	map {$f_dejavu{$_} ++} @{$res->{$k}->{fam}};
  	print $k."\t". $res->{$k}->{score}."\t".join("->",@{$res->{$k}->{fam}})."\n";
  }
  print "------------------------\n";
    foreach my $k (sort {$res->{$a}->{score} <=> $res->{$b}->{score} } keys %$res){
  	next if exists $f_dejavu{$res->{$k}->{fam}->[-1]};
  	print $k."\t". $res->{$k}->{score}."\t".join("->",@{$res->{$k}->{fam}})."\n";
  }
	

}


sub fast_plink3 {
	my ($vcf_file,$patients) = @_;

	my $dir_out = $project->getCallingPipelineDir("unifiedgenotyper").'/plink';
	

	
	my $dir = $project->getCallingPipelineDir("unifiedgenotyper");
	
	my $logPlink = $dir.'/plink/'.$projectName.'.plink.resume';
	
	my @snps;

	open (VCF, "<$vcf_file");
	my @samples;
	my $DP;
	while (my $line = <VCF>){
		next if $line =~ /##/;
		chomp($line);
		if ($line =~ /#CHR/){
			my @tdata = split(" ",$line);
			@samples = splice(@tdata,9);#$tdata[9..-1];
			next;
		}
		
		my ($chrom,$pos,$pid,$ref,$alt,$qual,$filter,$info,$format,@gsamples) = split(" ",$line); 
		next if $ref =~ /,/;
		next if $alt =~ /,/;
		my $snp;
		my $id = $chrom."_".$pos."_".$alt;
		$snp->{id} = $id;
		$snp->{position} = $pos;
		$chrom =~ s/chr//;
		$chrom = "MT" if $chrom eq "M";
		$snp->{chr} = $chrom;
		unless ($DP){
			my @t = split(":",$format);
			for (my $i =0;$i<@t;$i++){
				$DP=$i if $t[$i] eq "DP";
			}
			$DP =-1 unless $DP;
		}
		for (my $i=0;$i<@samples;$i++){
			my $name = $samples[$i];
			my $string = $gsamples[$i];
			my $geno;
			if ($string =~ /0\/0/) { $geno ="$ref $ref"; }
			elsif ($string =~ /0\/1/) { $geno ="$ref $alt"; }
			elsif ($string =~ /1\/0/) { $geno ="$ref $alt"; }
			elsif ($string =~ /1\/1/) { $geno ="$alt $alt"; }
			elsif ($string =~ /.\/./) { $geno ="0 0"; }
			elsif ($string =~ /./) { $geno ="0 0"; }
			else { die($string); }
			if ($DP > 0){
				my @tt = split(":",$string);
				my $v = $tt[$DP];
				$geno = "0 0" if $v < 10;
								
			}
			$snp->{samples}->{$name} = $geno;
		}
		
	
		#my ($find) = grep {$snp->{samples}->{$_}  eq "0 0"} @samples;
		#next if $find;
		push(@snps,$snp);
	}
	
	$dir .= "/plink_trio/";
	mkdir $dir unless -e $dir;
	  my $res;      

    	my %dejavu;  
    	 	my $nb =  @{$patients};   
    	my $cpt =0;
	foreach my $p1 (@{$patients}){
		my $res_p =[];	 
		my @samples_name;
		$cpt++;
		warn $p1->name()." $cpt/$nb\n";
		foreach my $p (sort {$a->name cmp $b->name} @{$project->getPatients()}){
			next if $p->name() eq $p1->name();
				push(@samples_name,$p->name());
		}
	
warn "start compute"; 
	 my $combinat = Math::Combinatorics->new(count => 2,
                                          data => [@samples_name],
                                         );
  
    	my @comb = combine(2,@samples_name) ;   
   	my $pr = String::ProgressBar->new( max => scalar(@comb) );
   	my $cc =0;
   	
    	foreach my $c (@comb){                     
  		print STDERR $pr->string()."\r";
		$cc++;
     		$pr->update($cc);
    		push(@$c,$p1->name());
#	 while(my @combo = $combinat->next_combination){
	 		my @toto = permute(@$c);
	 		my $dp;
	 	
    			foreach my $fam (@toto){
    				my $child = $fam->[-1];
    				next if exists $dejavu{$p1->name};
    				next if exists $dejavu{$fam->[0]};
    				next if exists $dejavu{$fam->[1]};
    				my $s = join("-",sort{$a cmp $b} ($fam->[0],$fam->[1]));
    				#warn $s;
    				next if exists $dp->{$s};
    				$dp->{$s} ++;
    				$res->{$child}->{score} = 999 unless exists $res->{$child};
    				my $tped_file = $dir."/".$project->name.".1.tped";
				my $ped_file = $dir."/".$project->name.".ped";
				open (PED,">".$ped_file);
				print PED join("\t",("F",$fam->[0],0,0,1,1))."\n";
				print PED join("\t",("F",$fam->[1],0,0,2,1))."\n";
				print PED join("\t",("F",$fam->[2],$fam->[0],$fam->[1],2,1))."\n";
				close PED;
				
			open (TPED,">".$tped_file);
			my $nb_snp =0;
		
			foreach my $snp (@snps){
					my @data;
					push(@data,$snp->{chr});
					push(@data,$snp->{id});
					push(@data,0);
					push(@data,$snp->{position});
				 my $found =0;
				 my %tt;
				 foreach my $name  (@$fam){
				 	$found ++ if (exists $snp->{samples}->{$name}  && $snp->{samples}->{$name} ne "0 0") ;
				 	$tt{$snp->{samples}->{$name}} ++;
				 	
				 }
				 next if $found <3;
				 next if scalar(keys %tt) == 1;
				 
				 $nb_snp ++;
				 last if $nb_snp > 1000;
				foreach my $name  (@$fam){
			
					if (exists $snp->{samples}->{$name} ) { push(@data,$snp->{samples}->{$name}); }
						else { warn "coucou";push(@data,$snp->{samples}->{$name}, "0 0");	}
					}
						print TPED join("\t",@data)."\n";
					}
					close TPED;
					my $cmd2 = "$plink --tped $tped_file --tfam $ped_file --noweb --mendel  --out $dir/$projectName.1";		
					my @log = `$cmd2`;
					open(MENDEL,"$dir/$projectName.1.imendel") || die("problem with $dir/$projectName.imendel \n$cmd2");
					while(my $line =<MENDEL> ){	
						chomp($line);
						next if $line =~ /FID/;
						my ($f, $name,$nb) = split(" ",$line);
						next if $name ne $child;
						#warn  $nb_snp;
						next if $nb_snp == 0;
						my $p = int($nb/$nb_snp*10000)/100;
						#warn $p."\n";
						if ($p < $res->{$child}->{score}){
							$res->{$child}->{score} = $p;
							$res->{$child}->{fam} = $fam;
							warn join("->",@{$fam})."\n"  if $p < 2;
							$dejavu{$child} ++ if $p < 1;
							if ($p<4){
								my $t;
								push(@$res_p, {score=>$p,fam=>$fam});
							}
							#$dejavu{$p1->name} ++ if $p < 1;
							#$res->{$p1->name}->{score} = $p;
							#$res->{$p1->name}->{fam} = $fam; 
						}
					
			}
					close(MENDEL);
					
    				
    			}
    		
    			
  	}
  	print STDERR "\n";
  		warn "NONE \n" unless scalar(@$res_p);
  		foreach my $r (@$res_p){
  			warn $r->{score}."\t".join("->",@{$r->{fam}})."\n";
  		}
  		
	}#end patients
  
  my %f_dejavu;
  foreach my $k (sort {$res->{$a}->{score} <=> $res->{$b}->{score} } keys %$res){
  	next if $res->{$k}->{score} > 3;
  	map {$f_dejavu{$_} ++} @{$res->{$k}->{fam}};
  	print $k."\t". $res->{$k}->{score}."\t".join("->",@{$res->{$k}->{fam}})."\n";
  }
  print "------------------------\n";
    foreach my $k (sort {$res->{$a}->{score} <=> $res->{$b}->{score} } keys %$res){
  	next if exists $f_dejavu{$res->{$k}->{fam}->[-1]};
  	print $k."\t". $res->{$k}->{score}."\t".join("->",@{$res->{$k}->{fam}})."\n";
  }
	

}



sub find_duo {
	my ($vcf_file,$patients) = @_;

	my $dir_out = $project->getCallingPipelineDir("unifiedgenotyper").'/plink';
	

	warn "read VCF";

	my $dir = $project->getCallingPipelineDir("unifiedgenotyper");
	
	my $logPlink = $dir.'/plink/'.$projectName.'.plink.resume';
	
	my @snps;

	open (VCF, "<$vcf_file");
	my @samples;
	my $DP;
	while (my $line = <VCF>){
		next if $line =~ /##/;
		chomp($line);
		if ($line =~ /#CHR/){
			my @tdata = split(" ",$line);
			@samples = splice(@tdata,9);#$tdata[9..-1];
			next;
		}
		
		my ($chrom,$pos,$pid,$ref,$alt,$qual,$filter,$info,$format,@gsamples) = split(" ",$line); 
		next if $ref =~ /,/;
		next if $alt =~ /,/;
		my $snp;
		my $id = $chrom."_".$pos."_".$alt;
		$snp->{id} = $id;
		$snp->{position} = $pos;
		$chrom =~ s/chr//;
		$chrom = "MT" if $chrom eq "M";
		$snp->{chr} = $chrom;
		unless ($DP){
			my @t = split(":",$format);
			for (my $i =0;$i<@t;$i++){
				$DP=$i if $t[$i] eq "DP";
			}
			$DP =-1 unless $DP;
		}
		for (my $i=0;$i<@samples;$i++){
			my $name = $samples[$i];
			my $string = $gsamples[$i];
			my $geno;
			if ($string =~ /0\/0/) { $geno ="$ref $ref"; }
			elsif ($string =~ /0\/1/) { $geno ="$ref $alt"; }
			elsif ($string =~ /1\/0/) { $geno ="$ref $alt"; }
			elsif ($string =~ /1\/1/) { $geno ="$alt $alt"; }
			elsif ($string =~ /.\/./) { $geno ="0 0"; }
			elsif ($string =~ /./) { $geno ="0 0"; }
			else { die($string); }
			if ($DP > 0){
				my @tt = split(":",$string);
				my $v = $tt[$DP];
				$geno = "0 0" if $v < 10;
								
			}
			$snp->{samples}->{$name} = $geno;
		}
		
	
		#my ($find) = grep {$snp->{samples}->{$_}  eq "0 0"} @samples;
		#next if $find;
		push(@snps,$snp);
	}
	close VCF;
	$dir .= "/plink_trio/";
	mkdir $dir unless -e $dir;
	  my $res;      
	$| =1;
    	my %dejavu;  
  	my $res;
    	my $cpt =0;
    	my @samples_name ;
    	map {push(@samples_name,$_->name)} @{$project->getPatients};
    	my %pt;
    	map{$pt{$_->name} ++} @$patients;
    	warn "start compute";
    	 	my @comb = combine(2,@samples_name) ;  
    	 my $nb =  scalar(@comb);   
    	my $pr = String::ProgressBar->new( max => scalar(@comb) );
   	my $cc =0;
	foreach my $fam (@comb){
		my $find = grep{exists $pt{$_}} @$fam;
		next unless $find;
		print STDERR $pr->string()."\r";
		$cc++;
     	$pr->update($cc);
		 my $score;
		 my $nb_snp =0;
			foreach my $snp (@snps){
					my @data;
					push(@data,$snp->{chr});
					push(@data,$snp->{id});
					push(@data,0);
					push(@data,$snp->{position});
				 my $found =0;
				 my %tt;
				 foreach my $name  (@$fam){
				 	$found ++ if (exists $snp->{samples}->{$name}  && $snp->{samples}->{$name} ne "0 0") ;
				 	$tt{$snp->{samples}->{$name}} ++;
				 	
				 }
				 next if $found ne 2;
				 $nb_snp ++;
				 $score ++ if scalar(keys %tt) == 1;
				 last if $nb_snp == 10000;
    			}
    			next if $nb_snp ==0;
    			my $p = int($score/$nb_snp*10000)/100;
    			
 
    			my $n1 = $fam->[0];
    			my $n2 = $fam->[1];
    			$res->{$n1}->{$n2} = $p;
    			$res->{$n2}->{$n1} = $p;
    			
    			
	}#end combi

  
  my %f_dejavu;
  
  foreach my $p (@$patients){
  	my $k = $p->name();
  	 my $th = $res->{$k};
  	 print $k." : ";
  	 my $cpt;
  	  foreach my $k2 (sort {$th->{$b} <=> $th->{$a} } keys %$th){
  	  	$cpt ++ if $th->{$k2} < 75;
  	  	last if $cpt > 2;
  	  	 print $k2.":".$th->{$k2}."% ";
  	  }
  	print "\n";
  	  	print "\n";
  }


}


