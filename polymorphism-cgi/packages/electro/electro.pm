package electro;

#use ABI; 
#use GenBoCacheQuery;
use strict;
#use connect;
use Data::Dumper;
use Carp;
use JSON::XS; 
#use util_file; 
#use Bio::SCF; 
my $cgi    = new CGI;

#ce script est utilisé pour tracer les electrophoregrammes des séquences pour le sequencage classique ou les histogrammes pour les puces de resequencage
my @bases = ( "A", "T", "C", "G" );
#fonction pour creer le cache permettant de tracer la variation graphiquement
sub construct_cache {
	my ($project,$contig,$variation,$trace,$position) = @_;	
	#chargement du buffer
	my $buffer = $project->buffer();
	#on recupère l'id du cache
	#my $json_id = GenBoCacheQuery::getid($buffer->dbh,$project->id,$variation->id,$trace->id);
	return unless $trace; 
	#si il s'agit d'une varaition trouvée par la methode gseq
	my $data;
	if ($project->projectType()->name() eq "array"){
		my $patient = $trace->getPatient();
		my $lecture = util_file::get_gseq_file({project=>$project,patient_name=>$patient->name(),method=>"affymetrix"});
		#my $lecture =  "/temporary/Sequence/analyses/". $project->name . "/intensity/" . $patient->name . "\.txt";
		#on recupère les données d'intensité correspondantes
		$data = gseqData($lecture,$project,$variation->position($variation->getContig())->start(),$trace);
	}
	#si projet de type ngs on ne fait rien
	elsif ($project->projectType()->name() eq "ngs"){
		return [];
	}
	#sinon on recupère les infos dans fichier phd de la trace possèdant la variation
	else{
		my @fichiers = glob("/temporary/Sequence/analyses/". $project->name . "/" . $contig->name . "/chromat_dir/*");
		my $extension; 
		foreach my $fic(@fichiers){
			#warn $fic;
			if($fic=~ /$\.ab./){
				my ($nom,$ext)=split(/\./,$fic);
				$extension=$ext;
				last();
			}
			if($fic=~ /$\.SC./){
				my ($nom,$ext)=split(/\./,$fic);
				$extension=$ext;
				last();
			}
		}
		#$extension = "SCF";
		my $lecture     =   "/temporary/Sequence/analyses/". $project->name() . "/" . $contig->name() . "/chromat_dir/" . $trace->name() . "\.".$extension;
		my $pos_var_contig;
		if ($variation) { 
			 $pos_var_contig = $variation->position($contig)->start;
		}
		else {
			 $pos_var_contig = $position;
		}
		#$data = electroData($lecture,$extension,$project,$contig,$variation,$trace,$pos_var_contig);
		$data = electroData($lecture,$extension,$project,$contig,$trace,$pos_var_contig);
	}
	return $data;
}

#fonction permettant de recupérer les données d'intensité pour une variation
sub gseqData{
	my ($lecture,$project,$position,$trace) = @_;
	warn $lecture;
	open (IN,$lecture);
	my @data = <IN>;
	close(IN);
	chomp(@data);
	my %all;
	my %intensity;
	my $pos =  $position;
	#on parcourt le fichier d'intensité et on recherche la ligne commencant par le nom du contig auquel appartient la variation
	my @elt; 
	foreach my $d (@data){
		my($premier) = split(/\t/,$d);
		my $contigName=substr($premier,0,rindex($premier,"-"));
		my $arrayPos = substr($premier,rindex($premier,"-")+1);
		$arrayPos++ if $arrayPos<3108;
		$arrayPos =$arrayPos +2 if $arrayPos>3107;
		$contigName =~ s/\+/-/;
		#on stcoke toutes les lignes correspondantes dans la table @donnee
		if($contigName eq "human_mtDNA_RCRS" && $arrayPos == $pos){
			warn $d;
			@elt = split(" ",$d);
		}
	}
	
	#on parse donc les données du tableau @donnee dont la ligne correspond à la position de la variation sur la trace
	my @values =@elt;
	shift(@values);
	my %h;
	$all{identifier} = "name";
	$all{label}      = "name";
	$h{name} = "toto";
	$h{type} = 2;
	$h{data} = \@values;
	push(@{$all{items}},\%h);
	my @mut; 
	push( @mut, \%h ); 
	#on renvoit les données correspondantes au format json directement dans le cache
	#print printJson2( $project->buffer,\@mut,$project->id,$variation->id,$trace->id);
	return \@mut;
}


sub ABI_Electro {
	my ($file) = @_;
	#on utilise la module ABI pour lire les données du fuchiers phd contenant les infos pour tracer les chromatogrammes
	my $abi = ABI->new( -file => "$file" )
 	 or die "impossible d'ouvrir le fichier abd $file";
  my %traces;
#pour chacune des bases on recupère le tracé de la base en fonction des positions
	foreach my $b (@bases) {
		@{ $traces{$b} } = $abi->get_trace($b);   # Get the raw traces for each base
	}
	return \%traces;
}
sub SCF_Electro2 {
	my ($hash) = @_;
	my @bases = ( "A", "T", "C", "G" );

	my $traces;
	
	
	my $samples = $hash->{samples};
	
 for (my $i = 0; $i<$hash->{samples_length}; $i++) {
 	foreach my $b (@bases){
 		$traces->{$b} =  $samples->{$b};
 
 	}
      
   }

 return $traces;
}
sub SCF_Electro {
	my ($scf) = @_;
	
	my $traces;
	
	
	#my $samples = $hash->{samples};

 for (my $i = 0; $i<$scf->samples_length; $i++) {
 	foreach my $b (@bases){
 		push(@{$traces->{$b}},$scf->sample($b,$i));
 
 	}
      
   }

 return $traces;
}

#fonction appelée pour recupérer les données pour tracer l'electrophoregramme de la variation
sub electroData {
my ($lecture,$extension,$project,$contig,$trace,$contig_position) = @_;

my $position_trace       =  $contig->getPositionOnTrace($trace,$contig_position);

my %tab;
#on recherche la trace correspondante à celle passée en paramètre
#my ($tv)  = grep {$_->id eq $trace->id} @{$variation->getTraces()};

my %tab;
my @mut;
my %traces;
if ($extension =~ /ab/){
	
	require "ABI.pm";
	(%traces) = %{ABI_Electro($lecture)};
}
else {
	my %hash;	
	require "Bio::SCF.pm";
	my $scf = Bio::SCF->new($lecture);
	(%traces) = %{SCF_Electro($scf)};
}
my $data_electro = readPhd($project,$trace,$contig,$extension);

my $seqPhd =  join("",map{uc($_->{lettre})} @$data_electro);

#on recherche la position de la varaition sur la trace

if ($position_trace-1 > scalar(@$data_electro)) {
	#confess();
	#$pos = $trace->position($contig)->length -1;
}

my $base_calls = $data_electro->[$position_trace-1]->{electro_pos};



#on stoke les infos sur la trace dans la tableau tab
my $direction = $trace->position($contig)->strand();
$tab{type} = 1;
$tab{strand} = $direction;
#$tab{allele} = $contig->sequence($variation)."/".$variation->sequence();
$tab{name} = $trace->name();


my %data;
my $z      = 1;




#longueur des traces
my $max;

#pos5 = valeur en position $base (position sur electrophregramme de la base 25 dans la séquence)
# - ecart entre position sur graph de la base et de la base précédente divisé par 2
my $limit = 100;
my $pos5 = $base_calls - $limit;

#pos 3 = = valeur en position $base (position sur electrophregramme de la base 25 dans la séquence)
# + ecart entre position sur graph de la base suivante et de la base divisé par 2
my $pos3    = $base_calls + $limit;
my $nb_plot = $pos3 - $pos5;
my $size    = $nb_plot;
my $echx;
my $echy;
#pour chacune des positions des traces et pour chaque des traces (A,C,G et T), la valeur du maximum des traces est assigné à $max,
for ( my $i = $pos5 ; $i < $pos3 ; $i++ ) {
	foreach my $b (@bases) {
		$max = $traces{$b}->[$i] if ( $max < $traces{$b}->[$i] );
	}
}

my $height = 100;
my $echy = ($height-5) / $max;
foreach $b (@bases) {
	my $posx = 0;
	my @array2;
	my $x = 0;

	#pour chaque valeur entre pos5 et pos3, ajouter points de pos5 à pos3
	if ( $direction == 1 ) {
		for ( my $i = $pos5 ; $i < $pos3 ; $i++ ) {
			my %base;
			$base{x} = $x;
			$x++;
			$base{y} = int($height - ( $traces{$b}->[$i] * $echy ) +20);
		
			push( @array2, \%base );

			$posx++;
		}
		$data{$b} = \@array2;
	}
	#complement des bases si trace est en -1
	elsif ( $direction == -1 ) {
		my $bc;
		$bc = uc($b);
		if ( $bc eq "A" ) {
			$bc = "T";
		}
		elsif ( $bc eq "T" ) {
			$bc = "A";
		}
		elsif ( $bc eq "G" ) {
			$bc = "C";
		}
		elsif ( $bc eq "C" ) {
			$bc = "G";
		}
		
		for ( my $i = $pos3 ; $i > $pos5 ; $i-- ) {
			my %base;
			$base{x} = $x;
			$x++;
			$base{y} = int(100 - ( $traces{$b}->[$i] * $echy ) + 20);
			push( @array2, \%base );

			$posx++;
		}
		$data{$bc} = \@array2;
	}
}
my $posontrace=0;
my @keepBases;



my $toto;




#pour associer les lettres à chacune des positions du graph
for my $d (@$data_electro){
	$posontrace ++;
	next if $d->{electro_pos} < $pos5;
	last if $d->{electro_pos} > $pos3;

	if ($direction ==1){
		$d->{electro_pos} = ($d->{electro_pos} - $pos5);
		 $d->{lettre} = uc($d->{lettre});
	}
	else {
		$d->{electro_pos} = ($pos3 -$d->{electro_pos});
		 $d->{lettre} = reverseAllele(uc($d->{lettre}));
		
	}
	$d->{variation} = -1;
	#warn $posontrace ." ".$posVariation;
	$d->{variation} = 1 if ($posontrace == $position_trace);

	$d->{lettreRef} =  $trace->getContigBase($contig,$posontrace,$d->{lettre});
	$toto .=  $d->{lettreRef};

	push(@keepBases,$d);
	
}
#die();
my %base;
$base{x1} = $limit;
$base{y1} = 120;
$base{x2} = $limit;
$base{y2} = 20;
####
#
####
$data{found} = -1;
$data{variation_id} = undef;
my ($find_var) = grep {$_->position($contig)->start == $contig_position} @{$trace->getVariations};

if ($find_var){	
	$data{found} = 1;
	$data{variation_id} = $find_var->id;	
}

######################################
$data{strand} = $direction;
$data{toto} = \%base;#array2;
$data{bases} = \@keepBases;#array2;
$data{trace_id} = $trace->id;

$data{project} = $project->name;
$tab{data}  = \%data;
push( @mut, \%tab );
return \@mut;
}

#fonction pour inserer les données en cache dans la BD (pour accelerer l'affichage)
sub printJson2 {
	my ($buffer,$data,$pid,$vid,$tid) = @_;
	
	my ($z) = encode_json($data);
	
	GenBoCacheQuery::insertJson($buffer->dbh,$pid,$vid,$tid,$z);
}

#fonction pour obtenir le complement d'une base
sub reverseAllele{
	my $bc =shift;
	
	if ( $bc eq "A" ) {
			$bc = "T";
		}
		elsif ( $bc eq "T" ) {
			$bc = "A";
		}
		elsif ( $bc eq "G" ) {
			$bc = "C";
		}
		elsif ( $bc eq "C" ) {
			$bc = "G";
		}
		return $bc;
}

sub readPhd {
 my ($project,$trace,$contig,$extension) = @_;	
	#################à revoir
my $fichier_phd =
    "/temporary/Sequence/analyses/"
  . $project->name() . "/"
  . $contig->name()
  . "/phd_dir/".$trace->name() . "\.".$extension."\.phd\.1";
 
###############
my @data_electro;
my $start2;

#on ouvre les fichier phd
open( PHD, "$fichier_phd" );

while ( my $lines = <PHD> ) {

	if ( $lines =~ /^BEGIN_DNA/ ) {
		$start2 = 1;
		$lines  = <PHD>;
	}
	if ( $lines =~ /^END_DNA/ ) {
		$start2 = 0;
	}
	#on recupère les lignes qui nous interressent dans la fichiers phd
	if ( $start2 == 1 ) {
		my %st;
		
		( $st{lettre}, $st{qualite}, $st{electro_pos} ) = split( /\s+/, $lines );
	
		push(@data_electro,\%st);
	}
}
close(PHD);

return \@data_electro;
}

1;