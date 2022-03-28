package file_util;
use strict;
use Data::Printer;
use File::Util;
use Carp;
use Data::Dumper;

sub return_patients {
	my ($project,$patients_name) = @_;
	my $patients;
if ($patients_name eq 'all'){
	$patients = $project->getPatients();
}
else {
	my @names = split(",",$patients_name);
	foreach my $name (@names){
		 unless ($project->getPatient($name)){
		 	warn "$name =>". join(" ",map {$_->name} @{$project->getPatients()});
		 	die();
		 }
	}
	map{push(@$patients,$project->getPatient($_))} split(",",$patients_name);
}
	return $patients;
}

sub count_files{
my ($patient,$ext) = @_;

my($f) = File::Util->new();

my(@files);
my $patern = "--pattern=".$patient."*".$ext."\$";
@files = $f->list_dir("/",'--files-only',$patern);

return scalar(@files);	
}

sub find_files {
my ($patients,$dir,$ext,$lane_number,$ext2) = @_;
$ext = "xsq" unless $ext;
$ext2 = "" unless $ext2;
my %dir_patients;
my $nb =$lane_number;
my($f) = File::Util->new();
foreach my $p (@$patients) {
	my $name = $p->name();
	
	my(@files);
	my $patern = "--pattern=\\.".$ext."\$";

	@files = grep {/^$name/} $f->list_dir($dir,'--files-only',$patern);
#	p @files;
	$dir_patients{$p->name()} = \@files;
	if ($lane_number>-1){
		warn " !!!! problem for $patern $dir  *".$p->name(). "* only ".scalar(@files) ." expected $nb" if scalar(@files) != $nb;
	}	
	confess("no files for ".$p->name()) unless  scalar(@files);  
}
return \%dir_patients;
}


sub find_fragment_file{


	my ($project,$name,$dir,$ext) = @_;
	my $patient = $project->getPatient($name);
	die("no patient with name : ".$name) unless $patient; 
	my $bc = $patient->barcode;
	my($f) = File::Util->new();
	my (@titi) = $f->list_dir($dir=> { 
			no_fsdots =>1,
			files_only => 1,
			files_match => {or => [qr/$name/ ]}
	});
	die() if (scalar(@titi) > 1);
	return $dir.$titi[0];
}


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
			my ($h11,$h12) = split(" ",$h1);
			my $cmd2 = "zcat $dir/".$cp->{R2} ."  2>/dev/null | head -1";
			my ($h2) = `$cmd2`;
			chomp($h2);
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


sub find_paired_files_pn {
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
	die()  if (scalar(@z)>1);

	my $good_pos = $z[0];
	my $hcouple = $couple->{$good_pos};
	return [] if (scalar(@$hcouple)*2 ne scalar(@$files));
	
	foreach my $cp (@$hcouple){
		my $cmd1 = "zcat $dir/".$cp->{R1} ." 2>/dev/null  | head -1";
		my ($h1) = `$cmd1 2>/dev/null`;
		chomp($h1);
		$h1 =~s/\/1//;
		$h1 =~s/\/2//;
		$h1 =~s/\s{1}1//;
		$h1 =~s/\s{1}2//;
		my $cmd2 = "zcat $dir/".$cp->{R2} ."  2>/dev/null | head -1";
		my ($h2) = `$cmd2`;
		chomp($h2);
		$h2 =~s/\/2//;
		$h2 =~s/\/1//;
		$h2 =~s/\s{1}2//;
		$h2 =~s/\s{1}1//;
		die() if $h1 ne $h2;
	}
	
	

	return $hcouple;
	
	
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
			warn "ici";
			$couple = find_paired_files(\@titi,$dir);
			last NAME if ($couple);
		}
	}
	warn $name;
	}
	warn "NO fastq file for : -".$patient->name()."- ".$patient->barcode." ".$dir unless $couple;	
	return $couple if scalar(@$couple)>0;
	warn "NO fastq file for : -".$patient->name()."- ".$patient->barcode." ".$dir;
	die();
}


sub find_file_pe_old {
	my ($patient,$ext) = @_;
	my $dir = $patient->getSequencesDirectory();
	my @names;
	push(@names,$patient->name);
	push(@names,$patient->barcode) if length($patient->barcode)>1 ;
	
	unless (exists $cached_dir{$dir}){
		$cached_dir{$dir}= [];
		opendir(DIR,$dir);
		my @allFiles= readdir(DIR);
		$cached_dir{$dir} = \@allFiles;
	}
	my @titi;
	foreach my $name (@names){
	
	my $find = "^".$name."_[ATGC][ATGC][ATGC][ATGC][ATGC][ATGC]";
	my (@titi) = grep { /^$find/} grep { /$ext/} @{$cached_dir{$dir}} ;
	return \@titi if (scalar(@titi));
 	$find = $name."_";
	my (@titi) = grep { /^$find/} grep { /$ext/} @{$cached_dir{$dir}} ;
	return \@titi if (scalar(@titi));
	 $find = "_".$name."_";
	my (@titi) = grep { /$find/} grep { /$ext/} @{$cached_dir{$dir}} ;
	return \@titi if (scalar(@titi));
	
		}
	confess("no files for ".$patient->name." in : ".$dir) ;  
#	my($f) = File::Util->new();
#my	$name = "^".$patient."_[ATGC][ATGC][ATGC][ATGC][ATGC][ATGC]";
#	
#	
#	
#	
#	
#	my (@titi) = $f->list_dir($dir=> { 
#			no_fsdots =>1,
#			files_only => 1,
#			files_match => {and => [qr/$name/ , qr/$ext/]}
#	});
#	
#	p (@titi);
#	unless  (scalar(@titi)){
#		 $name = "_".$patient."_";
#		(@titi) = $f->list_dir($dir=> { 
#			no_fsdots =>1,
#			files_only => 1,
#			files_match => {and => [qr/$name/ , qr/$ext/]}
#			
#	});
#	} 
#	unless  (scalar(@titi)){
#		 $name = $patient."_";
#		(@titi) = $f->list_dir($dir=> { 
#			no_fsdots =>1,
#			files_only => 1,
#			files_match => {and => [qr/$name/ , qr/$ext/]}
#			
#	});
#	} 

	
#	return \@titi;
}
sub find_patient_casava_pe_1 {
	my ($patient,$dir,$ext) = @_;
	
	
	unless (exists $cached_dir{$dir}){
		$cached_dir{$dir}= [];
		opendir(DIR,$dir);
		my @allFiles= readdir(DIR);
		$cached_dir{$dir} = \@allFiles;
	}
	my $find = "^".$patient."_[ATGC][ATGC][ATGC][ATGC][ATGC][ATGC]";
	my (@titi) = grep { /^$find/} grep { /$ext/} @{$cached_dir{$dir}} ;
	return \@titi if (scalar(@titi));

	 $find = "_".$patient."_";
	 
	my (@titi) = grep { /^$find/} grep { /$ext/} @{$cached_dir{$dir}} ;
	return \@titi if (scalar(@titi));
	 $find = $patient."_";
	my (@titi) = grep { /^$find/} grep { /$ext/} @{$cached_dir{$dir}} ;
	return \@titi if (scalar(@titi));
	confess("no files for ".$patient." in : ".$dir) ;  
#	my($f) = File::Util->new();
#my	$name = "^".$patient."_[ATGC][ATGC][ATGC][ATGC][ATGC][ATGC]";
#	
#	
#	
#	
#	
#	my (@titi) = $f->list_dir($dir=> { 
#			no_fsdots =>1,
#			files_only => 1,
#			files_match => {and => [qr/$name/ , qr/$ext/]}
#	});
#	
#	p (@titi);
#	unless  (scalar(@titi)){
#		 $name = "_".$patient."_";
#		(@titi) = $f->list_dir($dir=> { 
#			no_fsdots =>1,
#			files_only => 1,
#			files_match => {and => [qr/$name/ , qr/$ext/]}
#			
#	});
#	} 
#	unless  (scalar(@titi)){
#		 $name = $patient."_";
#		(@titi) = $f->list_dir($dir=> { 
#			no_fsdots =>1,
#			files_only => 1,
#			files_match => {and => [qr/$name/ , qr/$ext/]}
#			
#	});
#	} 

	
#	return \@titi;
}
sub find_patient_casava_file {
my ($patients,$dir) = @_;
warn "new";
my $ext = "*";
my %dir_patients;

foreach my $p (@$patients) {
	my $name = $p->name()."_";
	
	my(@files);
	my $patern = "--pattern=\\.".$ext."\$";

	my($f) = File::Util->new();

	@files =  $f->list_dir($dir=>{no_fsdots =>1,
			files_only => 1,
			files_match => qr/$name/, 
			recurse=>0,
	});
	p @files;
#warn @files;
	$dir_patients{$p->name()} = \@files;
	confess("no files for ".$p->name()) unless  scalar(@files);  
}
warn "end";
return \%dir_patients;
}

1;
