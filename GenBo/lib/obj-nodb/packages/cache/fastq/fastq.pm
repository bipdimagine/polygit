package file_util;
use strict;
use Data::Printer;
use File::Util;
use Carp;
use Data::Dumper;



sub find_paired_files {
	my ($files,$dir) = @_;
	my $couple;
	my $associated;
	for (my $i=0 ;$i<@$files;$i++){
		my @f1 = split("",$files->[$i]);
		
		for (my $j=$i+1 ;$j<@$files;$j++){
			my $adif =[];
			my @f2 = split("",$files->[$j]);
			next if scalar(@f2) ne scalar(@f1);
			my $error =0;
			for (my $z =0 ;$z<@f1;$z++){
				if ($f1[$z] eq $f2[$z] ){
					next;
				}
				#last if $error > 1;
				my $aa = join("",sort ($f2[$z],$f1[$z] ) );
				$error ++;
				$error ++ if $aa ne "12" ;
				last if $error > 1;
		
			#	$error ++ if $aa ne "12" ;
			
			#	last  if $aa ne "12" ;
				my $dif;
				my $ad = $f2[$z];
				$ad =  $f1[$z] if ( $f1[$z] < $f2[$z]);
			
				$dif->{pos} = $z;
				$dif->{dif} = $ad;
				if ($f1[$z] eq "1"){
					$dif->{R1} = $files->[$i] ;
					$dif->{R2} = $files->[$j] ;
				}
				else {
					$dif->{R2} = $files->[$i] ;
					$dif->{R1} = $files->[$j] ;
				}
				
				push(@$adif,$dif);
				
			}
			if ($error == 1){ # &&  scalar(@$adif) ==1 ){
				foreach my $d (@$adif){
				$associated->{$d->{R1}} ++;
				$associated->{$d->{R2}} ++;
				push(@{$couple->{$d->{pos}}},$d);
				}
			}
			
			
		}
		
	}
	
 	#simple case 
	my @pos = sort{scalar($couple->{$a}) <=> scalar($couple->{$b}) } keys %{$couple};
	my (@z) = grep {scalar(@{$couple->{$_}}) *2 eq  scalar(@$files)} keys %{$couple};
	return [] unless @z;
#	die()  if (scalar(@z)>1);

	my $hcouple;
	my %bon;
	foreach my $zz (@z){
		$bon{$zz} ++;
		$hcouple = $couple->{$zz};
		return [] if (scalar(@$hcouple)*2 ne scalar(@$files));
		foreach my $cp (@$hcouple){
			my $cmd1 = "zcat $dir/".$cp->{R1} ." 2>/dev/null  | head -1";
			my ($h1) = `$cmd1 2>/dev/null`;
			chomp($h1);
			$h1 =~s/\// /g;
			my ($h11,$h12) = split(" ",$h1);
			my $cmd2 = "zcat $dir/".$cp->{R2} ."  2>/dev/null | head -1";
			my ($h2) = `$cmd2`;
		
			chomp($h2);
				$h2 =~s/\// /g;
			my ($h21,$h22) = split(" ",$h2);
			my @h23 = split(":",$h22);
			my @h13 = split(":",$h12);

		#	if($h11 ne $h21 && $h23[-1] ne $h13[-1]){
			if($h11 ne $h21){
				delete $bon{$zz};
				last;
			}
			
		}
	
	}
	my @nb =  keys %bon;
	die("pb couple seq") if scalar(@nb) >1;
	return $couple->{$nb[0]};
	die("probleme couple ")  ;
	
}

sub find_paired_files_umi {
	my ($files,$dir) = @_;
	my @r1 = grep{$_  =~ /_R1_/} @$files;
	my @all;
	my @couple;
	foreach my $r1 (@r1){
		my $d;
		$d->{R1} = $r1;
		
		my $r2 = $r1;
		if (-e $r2){
			$r2 =~ s/_R1_/_R2_/;
			$d->{R2} = $r2;
		}
		else {
			confess($r1 ." don t find R2" );
		}
		my $r3 = $r1;
		$r3 =~ s/_R1_/_R3_/;
		if ( -e $r3) {
			$d->{R3} = $r3;
		}
		push(@couple,$d);
	}
	
	return \@couple;
}


my %cached_dir;

sub find_file_pe {
	my ($patient,$separator) = @_;
	my $dir = $patient->getSequencesDirectory();
	my @names;
	my $couple;
	push(@names,$patient->name);
	push(@names,$patient->barcode) if length($patient->barcode)>1 ;
	unless (exists $cached_dir{$dir}){
		$cached_dir{$dir}= [];
		opendir(DIR,$dir);
		my @allFiles= readdir(DIR);
		$cached_dir{$dir} = \@allFiles;
	}
	NAME: foreach my $name (@names){
	my @pattern = ("^".$name."_[ATGC][ATGC][ATGC]","^".$name."_S[1-9]+","^".$name."_","$name");
	foreach my $find (@pattern){
		my (@titi) = grep { /$find/} grep { /fastq/} @{$cached_dir{$dir}} ;
		
		if (@titi) {
			$couple = find_paired_files_umi(\@titi,$dir);
			last NAME if ($couple);
		}
	}
	}
	confess("NO fastq file for : -".$patient->name()."- ".$patient->barcode." ".$dir) unless $couple;	
	return $couple if scalar(@$couple)>0;
	die();
}





1;
