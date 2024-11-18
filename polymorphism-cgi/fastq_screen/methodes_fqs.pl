use strict; 
use Data::Dumper;
use List::Util qw(max);
use Getopt::Long;

#def de variables
my $nom; 
my $loc; 
my $infile;
my $conf; 
my $subset = 100000; 
my $thread = 1; 
my $espece; 



#options de la fonction getopt::Long 
GetOptions(
	'infile|i=s' => \$infile,
	'outfile|o=s' => \$loc,
	'conf|c=s' => \$conf,
	'espece|e=s' => \$espece, 
	'subset|s=s' => \$subset,
	'threads|t=s' => \$thread, 
	'name|n=s' => \$nom, 
);


#Warnings si certaines options ne sont pas presentes
die ("\n\nERROR: option -infile is mandatory. Die....\n\n") if (not $infile);
die ("\n\nERROR: option -outfile is mandatory. Die....\n\n") if (not $loc);
die ("\n\nERROR: option -config file is mandatory. Die....\n\n") if (not $conf);
#die ("\n\nERROR: option -espece is mandatory. Die....\n\n") if (not $espece);
die ("\n\nERROR: option -name is mandatory. Die....\n\n") if (not $nom);




#Warnings si les fichiers pour certaines options n'existent pas 
die("\n\nERROR: -infile does not exist. Die....\n\n") unless (-e $infile);
die("\n\nERROR: -outfile does not exist. Die....\n\n") unless(-e $loc);
die("\n\nERROR: -config file does not exist. Die....\n\n") unless(-e $conf); 



#enleve les derniers caractères du nom du fichiers fastq --> on enleve ".fastq.gz"
my @split = split("/", $infile); 
my $nom_fichier_temporaire = $split[-1] ;
my $nom_fichier; 
if ($infile =~ ".gz" and $infile =~ ".fastq"){
	$nom_fichier = substr($nom_fichier_temporaire, 0, length($nom_fichier_temporaire)- 9);
}
elsif ($infile =~ ".fastq"){
	$nom_fichier = substr($nom_fichier_temporaire, 0, length($nom_fichier_temporaire)- 6);  
}
else {
	$nom_fichier = $nom ;
}


my $nom_dossier = $loc . "fastq_screen_". $nom ; 
if (not -d $nom_dossier){
	my $cmd_fichier = "mkdir $nom_dossier"; 
	system($cmd_fichier);
}


#fastq_screen
my $cmd = "fastq_screen --outdir $loc --conf $conf $infile -force --subset=$subset --threads=$thread"; 
system($cmd);



#recupere le nom du fichier fastq pour nommer les fichiers

#warn "\n\n";

my $nom_screen = $nom_fichier . "_screen"; 

#warn "nom_screen: $nom_screen";


my $outfile;
my $nom_dossier_sortie = "$nom_screen" . "_fichier_sortie"; 

#warn "nom_dossier_sortie: $nom_dossier_sortie";


my $base_files = $infile;
$base_files =~ s/.fastq.gz//;
$base_files =~ s/.+\///;
$base_files = $loc.'/'.$base_files;

my $cmd_mv = "mv $base_files* $nom_dossier/.";
system($cmd_mv);

my $txt = $infile;
$txt =~ s/.fastq.gz/_screen.txt/;
$txt =~ s/.+\///;
$txt = $nom_dossier.'/'.$txt;

#warn "nom_dossier fin: $txt";


die ("\n\nERROR PATH $txt not exists\n\n") if not -e $txt;

$outfile = "$nom_dossier/$nom_dossier" . "_screen.tab"; 

#warn "outfile: $outfile";

my $nom_espece = "$nom_dossier/$nom". "_screen" . "_nom_espece.txt"; 


#warn "nom_espece: $nom_espece";

my $genome = "ERROR: possible contamination"; #si impossible de déterminer le génome

$genome = det_espece($txt);
print "\n\n$genome\n\n"; 

#creer fichier contenant uniquement le nom de l'espece 
open(FILE, ">$nom_espece"); 
print FILE $genome;
close(FILE);



if ($genome ne $espece and $espece ){
	die( "\n\n\nERROR: wrong species given \n\n\n ");
}

sub lecture {
	my ($lecture_file) = @_; 
	my %h;
	open(FILE, "$lecture_file");
	#warn $lecture_file;
	while(<FILE>){
		chomp($_);
		if ($_ =~ /%Hit_no_genomes:/) {
			my @split = split(" ", $_); 
			$h{'%NoGenome'} = $split[-1];
		}
		elsif ($_ !~ "#" and $_ !~ "Genome" and $_ !~ "%" and $_ =~ /[A-Za-z]/ ){
			my @split = split("\t", $_); 
			$h{$split[0]}{"%unmapped"} = $split[2];
			$h{$split[0]}{"%one hit one genome"} = $split[5];
			$h{$split[0]}{"%multiple hits one genome"} = $split[7];
		}
	}
	close(FILE);
	return %h; 
}

sub det_espece {
	my ($filein) = @_; 
	my %h = lecture($filein);
	my @liste = (100, 0 x 2);
	my $perc_nogenome = $h{'%NoGenome'};
	return 'no_genome' if int($perc_nogenome) >= 10;
	my $geno; 
	foreach my $gen (keys %h){
		next if $gen eq '%NoGenome';
		if ($h{$gen}{"unmapped"} <= $liste[0]){
			if ($h{$gen}{"%multiple hits one genome"} >= $liste[1]){
				if ($h{$gen}{"%one hit one genome"} >= $liste[2]){
					$geno = $gen; 
					$liste[0] = $h{$gen}{"unmapped"}; 
					$liste[1] = $h{$gen}{"%multiple hits one genome"}; 
					$liste[2] = $h{$gen}{"%one hit one genome"}; 
				}
			}
		}
	}
	return $geno; 
}






##subs
#sub parsing {
#	 
#	my ($filein) = @_;
#	my %hash;
#	my @list = (0, 0, 0, 100, 0 x 11); 
#	open(FILE, $filein);
#	while (<FILE>) {
#		chomp($_); 
#		if ($_ =~ /[A-Za-z]/ and $_ !~ "#" and $_ !~ "%" and $_ !~ "Genome" ) { 
#			my @l = split(" ", $_); 
#			$hash{$l[0]}{$l[1]}{'%unmapped'} = $l[3]; 
#			$hash{$l[0]}{$l[1]}{'%one hit one genome'} = $l[5]; 
#			$hash{$l[0]}{$l[1]}{'%multiple hits one genome'} = $l[7]; 
#			$hash{$l[0]}{$l[1]}{'%one hit multiple genomes'} = $l[9]; 
#			$hash{$l[0]}{$l[1]}{'%multiple hits multiple genomes'} = $l[10]; 
#			
#			
#			if ($l[3] <= $list[3]) {  #unmapped
#				if ($l[5] >= $list[5]) {  #one hit one genome
#					if ($l[7] >= $list[7]) {	#multiple hits one genome 
#						$genome = $l[0]; 
#						
#						@list = (); 
#						push(@list, $_);
#			 
#		}}}}}
#	close(FILE);
#	
#	open(FILE, ">$outfile"); 
#	
#	foreach my $key_gen (keys %hash) {
#		my %hash_readstot = %{$hash{$key_gen}};
#	
#		foreach my $key_reads_tot (keys %hash_readstot){
#			my %hash_hits= %{$hash_readstot{$key_reads_tot}}; 
#				
#			my @keys_hash = keys %hash_hits;
#			my @values_hash_hits = values %hash_hits; 
#				
#			my $ligne2 = "$key_gen	\t $key_reads_tot	";
#				
#			print FILE  $ligne2 ;
#				
#			foreach my $key_hits (@keys_hash){
#				my $lignedude22 = "\t $key_hits--> $hash_hits{$key_hits} :	"; 
#				print FILE ($lignedude22); }	
#				print FILE "\n"; 	
#		}
#	}
#	close(FILE);
#	return $genome; 
#}
#
#

 




