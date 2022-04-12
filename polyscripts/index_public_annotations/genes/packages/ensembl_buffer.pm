package ensembl_buffer;
use strict;
use Config::Std; 
use Data::Dumper;
use lib "/software/distrib/ensembl-variation/modules";
use lib "/software/distrib//ensembl/modules";
use lib "/software/distrib//bioperl-live/";
use lib "/software/distrib//ensembl-variation/modules";
use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix;
my $species;


sub new {
	my ($class, @hash) = @_;
	
	my $self ={};
	bless ($self);
	$self->getConfig();
	warn "WARNING WORKING ON ".$self->getBuild();
	return $self;
}

sub getBuild {
	my $self = shift;
	return $self->getConfig->{global}->{current_build};
}
sub getdbh_variations{ 
	my ($self,$version) = @_;
	return $self->{getdbh_variations} if exists  $self->{getdbh_variations};
	my $ip =  $self->{config}->{ensembl}->{ip};
	my $db_user_name = $self->{config}->{ensembl}->{user};
	my $db_password = "";
	 my $genome_ensembl_release = $self->{config}->{global}->{current_ensembl_build};
	my $genome_release = $self->{config}->{global}->{current_release};
	my $ensembl_release = $self->{config}->{ensembl}->{$genome_release};
	warn $self->get_ensembl_version;
	my $dsn = "DBI:mysql:homo_sapiens_variation_".$self->get_ensembl_version."_".$genome_ensembl_release.":$ip:5306\n";
	my $dbh = DBI->connect($dsn, $db_user_name, $db_password)|| die "Database connection not made: $DBI::errstr";
	return $dbh;
}  

sub get_genome_release {
	my ($self,$version) = @_;
	return $self->{config}->{global}->{current_build};
}

sub get_ensembl_version{
		my ($self,$version) = @_;
	return$self->{config}->{ensembl}->{$self->get_genome_release};
	
}

sub get_genome_ensembl_release{
	my ($self,$version) = @_;
	return $self->{config}->{global}->{current_ensembl_build};
}
sub setRegistry {
	my ($self,$version) = @_;
	#push(@INC,"/bip-d/perl/ensembl59/ensembl/");
	my $genome_release = $self->{config}->{global}->{current_build};
	my $version = $self->{config}->{ensembl}->{$genome_release};

	my $path_ensembl = $self->{config}->{ensembl}->{base}.$self->{config}->{ensembl}->{base_end};
	warn $path_ensembl;
	$species = "Human" if $version > 61;
	push(@INC,$self->{config}->{ensembl}->{base}."/ensembl/modules");
	#push(@INC,$self->{config}->{ensembl}->{base}.$self->{config}->{ensembl}->{$gp}."/ensembl-variation/modules");	
	
	require($path_ensembl."/Registry.pm");
	
	$self->{registry} = "Bio::EnsEMBL::Registry";
	
	#require($self->{registry});
		$self->{registry}->load_registry_from_db(-host => $self->{config}->{ensembl}->{ip},
									  -user => $self->{config}->{ensembl}->{user},
									   -db_version => $version,
									 
	 											);
	 											
					 					  
}

sub getMirnomeFile {
	my $file =  $_[0]->{config}->{public_data}->{HG19}."/mirbase/".$_[0]->{config}->{mirbase}->{version}."/".$_[0]->{config}->{mirbase}->{file};
	die() unless -e $file;
	return $file;
}


sub getMatrixAdaptor {
	return $_[0]->{matrix_adaptor} if exists  $_[0]->{matrix_adaptor};
	
	my $db = $_[0]->getRegistry()->get_DBAdaptor( $species, "variation" );
	$_[0]->{matrix_adaptor} = $db->get_ProteinFunctionPredictionMatrixAdaptor;
	return $_[0]->{matrix_adaptor};
}


sub getSliceAdaptor {
#	warn "WARNING WORKING ON ".$_[0]->getBuild();
	return $_[0]->{slice_adaptor} if exists  $_[0]->{slice_adaptor};
	my $db = $_[0]->getRegistry()->get_DBAdaptor( $species, "core" );

	$_[0]->{slice_adaptor} = $db->get_SliceAdaptor();
	return $_[0]->{slice_adaptor};
}

sub getGeneAdaptor {
	#warn "WARNING WORKING ON ".$_[0]->getBuild();
	return $_[0]->{gene_adaptor} if exists  $_[0]->{gene_adaptor};
	my $db = $_[0]->getRegistry()->get_DBAdaptor( $species, "core" );
	
	$_[0]->{gene_adaptor} = $db->get_GeneAdaptor();
	return $_[0]->{gene_adaptor};
}

sub getRegulationAdaptor {
	return $_[0]->{reg_adaptor} if exists  $_[0]->{reg_adaptor};
	#my $db = $_[0]->getRegistry()->get_DBAdaptor( $species, "Funcgen","RegulatoryFeature" );
	
	$_[0]->{reg_adaptor} = $_[0]->getRegistry()->get_adaptor($species, 'Funcgen', 'RegulatoryFeature');
	return $_[0]->{reg_adaptor};
}
sub getTranslationAdaptor {
	#warn "WARNING WORKING ON ".$_[0]->getBuild();
	return $_[0]->{translation_adaptor} if exists  $_[0]->{translation_adaptor};
	my $db = $_[0]->getRegistry()->get_DBAdaptor($species , "core");
	
	
	$_[0]->{translation_adaptor} = $db->get_TranslationAdaptor();
	return $_[0]->{translation_adaptor};
}

sub getTranscriptAdaptor {
	#warn "WARNING WORKING ON ".$_[0]->getBuild();
	return $_[0]->{transcript_adaptor} if exists  $_[0]->{transcript_adaptor};
	my $db = $_[0]->getRegistry()->get_DBAdaptor($species , "core");
	
	
	$_[0]->{transcript_adaptor} = $db->get_TranscriptAdaptor();
	return $_[0]->{transcript_adaptor};
}

sub getRegistry {
	my ($self) = shift;
	$self->setRegistry() unless exists  $self->{registry};
	return $self->{registry};
}


sub getConfig{
	my ($self) = @_;
	read_config "../packages/ensembl.cfg" => my %config;
	$self->{config} = \%config;
	
}





sub getDescription {
	my ($self,$description) = @_;
	my @dd = split (";",$description);
	my %desc;
	foreach my $d (@dd){
		my ($a,$b) = split("=",$d);
		$desc{lc($a)} = $b;
	}
	return \%desc;
}
sub getMirData {
	my $self = shift;
	my $file_mir = $self->getMirnomeFile();
	open (MIR,"zcat $file_mir |");
	my $mirs_primary;
	my $mirs;
	while (my $line = <MIR>){
	chomp($line);
	next if $line =~ /^#/;
	my ($chr_ucsc,$point,$type,$start,$end,$point2,$strand,$point3,$description) = split(" ",$line);
	my $hmir;
	$hmir->{start} = $start;
	$hmir->{end} = $end;
	$hmir->{strand} = 1;
	$hmir->{strand} = -1 if $strand ne "+";
	
	my $hdesc = $self->getDescription($description);
	$hdesc->{id} =~ s/_/-/;
	$hdesc->{derives_from}=~ s/_/-/;
	$hmir->{id} = $hdesc->{id};
	$hmir->{name} = $hdesc->{name};
	$hmir->{description} = $description;
	$hmir->{type} = $type;
	
	my $chr = $chr_ucsc;
	$chr =~ s/chr//;
	$chr = "MT" if $chr eq "M";
	if ($hmir->{type} eq "miRNA_primary_transcript"){
		$mirs_primary->{$chr}->{id}->{$hmir->{id}} = $hmir;
		$mirs_primary->{$chr}->{pos}->{$start."-".$end} = $hmir;
	}
	else {
		
		die($hdesc->{derives_from}) unless exists $mirs_primary->{$chr}->{id}->{$hdesc->{derives_from}};
		my $mir1 = $mirs_primary->{$chr}->{id}->{$hdesc->{derives_from}};
		$mir1->{span_mature} =  Set::IntSpan::Fast::XS->new() unless exists $mir1->{span_mature}; 
		$mir1->{span_mature}->add_from_string($start."-".$end); 
		#warn $hdesc->{derives_from} if exists $mirs_primary->{$chr}->{$hdesc->{derives_from}};
	}
	#push(@{$mirs->{$chr}},$hmir);
	 
}

 return $mirs_primary;
}
1;