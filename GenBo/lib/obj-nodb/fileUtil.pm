package file_util;
use strict;
use Data::Printer;
use File::Util;
use Carp;


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
#warn("ICIIII ".$ext." .. ".$ext2."..");
my %dir_patients;
my $nb =$lane_number;
my($f) = File::Util->new();
foreach my $p (@$patients) {
	my $name = $p->name().$ext2;
	
	my(@files);
	my $patern = "--pattern=\\.".$ext."\$";
#	warn($patern);
#	my @toto = $f->list_dir($dir,'--files-only',$patern);
#	p @toto;
#	p  $f->list_dir($dir,'--files-only',$patern);
#	warn("LAAAAA");
#	warn($ext2);
	@files = grep {/$name/} $f->list_dir($dir,'--files-only',$patern);
#	p @files;
	$dir_patients{$p->name()} = \@files;
	if ($lane_number>-1){
		warn " !!!! problem for $patern $dir  *".$p->name(). "* only ".scalar(@files) ." expected $nb" if scalar(@files) != $nb;
	}	
	confess("no files for ".$p->name()) unless  scalar(@files);  
}
return \%dir_patients;
}



1;