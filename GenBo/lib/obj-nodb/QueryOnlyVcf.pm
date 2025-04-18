package QueryOnlyVcf;

use strict;
use GBuffer;
use GenBoChromosome;
use GenBoVariation;
use GenBoInsertion;
use GenBoDeletion;
use Moo;
use Data::Dumper;
use Compress::Snappy;
use Carp;
use List::MoreUtils qw{ natatime };
use Storable qw/thaw freeze/;
extends 'QueryVcf';



has buffer => (
	is	 => 'ro',
	lazy => 1,
	default => sub{ 
		my $self = shift;
		my $buffer = GBuffer->new();
		return $buffer;
	}
);

has project => (
	is	 => 'rw',
	lazy => 1,
	default => sub { 
		my $self = shift;
		my $project_name;
		if ($self->release()) {
			if ($self->release() =~ /HG19/) {
				$project_name = $self->buffer->getRandomProjectName($self->release(), '46');
			}
			else {
				$project_name = $self->buffer->getRandomProjectName($self->release());
			}
		}
		else { $project_name = $self->buffer->getRandomProjectName(); }
		my $project = $self->buffer->newProject( -name => $project_name );
		return $project;
	}
);

has release => (
	is		=> 'rw',
	default => undef,
);

has patient => (
	is		=> 'rw',
	default => undef,
);

has method => (
	is		=> 'rw',
	default => 'no_method',
);

has parse_hgmd => (
	is		=> 'rw',
	default => undef,
);

has file => (
	is		=> 'rw',
	default => undef,
);

has version => (
	is		=> 'rw',
	default => 'HG19',
);

has annotation_version => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $last_annot_version = $self->max_gencode_version().'.'.$self->max_public_data_annot();
		return $last_annot_version;
	},
);

has max_public_data_annot => (
	is		=> 'rw',
	default => sub {
		my $self = shift;
		my $buffer = new GBuffer;
		my $annotdb = $buffer->getQuery->getMaxPublicDatabaseVersion();
		$buffer = undef;
		return $annotdb;
	},
);

has max_gencode_version => (
	is		=> 'rw',
	default => sub {
		my $self = shift;
		my $buffer = new GBuffer;
		my $genecode = $buffer->getQuery->getMaxGencodeVersion();
		$buffer = undef;
		return $genecode;
	},
);

has getObjects => (
	is		=> 'rw',
	default => undef,
);

has genomeFasta => (
	is		=> 'rw',
	lazy 	=> 1,
	reader  => 'getGenomeFasta',
	default => undef,
);

has genomeFai => (
	is		=> 'rw',
	lazy 	=> 1,
	reader	=> 'getGenomeFai',
	default => undef,
);

has hashChromosomes => (
	is		=> 'rw',
	lazy 	=> 1,
	default => sub {
		my $self = shift;
		my $hash;
		foreach my $hchr (@{$self->getGenomeFai()}){
			#next unless ($chr_name eq $hchr->{'name'});
			
			$hchr->{fasta_name} = $self->project->dirGenome()."/genome/ucsc/chr".$hchr->{id}.".fa";
			
			my $obj = GenBoChromosome->new($hchr);
			$obj->{genomeFasta} = $self->getGenomeFasta();
			$hash->{$obj->id()} = $obj
		}
		return $hash;
	},
);

has passVarAlreadyDone => (
	is		=> 'rw',
	default => undef,
);

has hashVarIdsDone => (
	is		=> 'rw',
	default => undef,
);

has noPatient => (
	is		=> 'rw',
	default => undef,
);

has force_all_gt_he => (
	is		=> 'rw',
	default => undef,
);



########## METHODS ##########



#sub newPseudoProject {
#	my ($self, $project_name) = @_;
#	my $test = 0;
#	
#	confess();
#	
#	my $annot_version = $self->annotation_version();
#	my $max_gencode = $self->max_gencode_version();
#	my $max_db_annot = $self->max_public_data_annot();
#	my @lTmp = split('\.', $annot_version);
#	my $project = GenBoProject -> new ( name	=> $project_name,
#										version => $self->version(),
#										genome_version_generic => $self->version(),
#										genome_version => $self->version(),
#										gencode_version => $max_gencode,
#										annotation_version => $max_db_annot,
#										public_database_version => $max_db_annot,
#										release => undef,
#										buffer	=> $self->buffer() );
#	$self->buffer->genome_version($project->gencode_version());	
#	$self->buffer->annotation_version($project->public_database_version());		
#	#$self->buffer->lmdb_public_dir($project->annotation_public_path);
#	my $buffer2 = GBuffer->new();
#	my $project2 = $buffer2->newProject ( -name	=> $buffer2->getRandomProjectName() );
#	#$project2->changeAnnotationVersion($annot_version, '1');
#	
#	
#	$project->public_data_root($project2->public_data_root());
##	$self->buffer->public_data_version($project2->public_data_version());
##	$self->buffer->lmdb_public_dir($project2->annotation_public_path);
#	$self->buffer->annotation_version($project2->annotation_version());
#	$self->buffer->annotation_genome_version($project2->annotation_genome_version());
#	$self->buffer->public_data_version($project2->buffer->public_data_version());
#	$self->buffer->{public_data} = $project2->buffer->public_data();
#	
#	$project2 = undef;
#	$buffer2 = undef;
#	return $project;
#}

sub getReferences {
	my $self = shift;
	return $self->project->getReferences();
}

sub get_regions {
	my ($self, $chr, $fork ) = @_;
	if ( $fork == 1 ) {
		my $regions;
		my $hregions;
		$hregions->{part_id} = 1;
		$hregions->{chromosome} = $chr->name;
		$hregions->{start}      = 1;
		$hregions->{end}        = $chr->length;
		push( @$regions, $hregions );
		return \@$regions;
	}
	my $tabix = $chr->buffer->software("tabix");
	my %hv;
				
	my $file = $self->file();
#	warn $file;
	die() unless -e $file;
	next unless -e $file;
#	warn $file." ok";
#	warn "$tabix $file " . $chr->fasta_name() . " | cut -f 2 |" ;
	open( TOTO, "$tabix $file " . $chr->fasta_name() . " | cut -f 2 |" );
	while ( my $pos = <TOTO> ) {
		chomp($pos);
		$hv{$pos}++;
	}
	close TOTO;
	
	my @snps = sort { $a <=> $b } keys %hv;
	warn "Found snps : ".scalar(@snps)."\n";
	if ( scalar(@snps) == 0 ) {
		my $hregions;
		$hregions->{part_id} = 1;
		$hregions->{chromosome} = $chr->name;
		$hregions->{start}      = 1;
		$hregions->{end}        = $chr->length;
		return [$hregions];
	}
	my $nb;
	if ( $fork == 1 ) {
		die();
		$nb = scalar(@snps);
	}
	else {
		$nb = int( scalar(@snps) / ( $fork - 1 ) );
	}
	$nb = 2_000 if $nb > 2_000;
	my $regions;
	my $iter = natatime $nb, @snps;
	my $old_end;
	my $i = 1;
	while ( my @tmp = $iter->() ) {
		my $hregions;
		$hregions->{part_id} = $i;
		$hregions->{chromosome} = $chr->name;
		if ($old_end) {
			$hregions->{start} = $old_end + 1;
		}
		else {
			$hregions->{start} = 1;
		}
		$hregions->{end} = $tmp[-1] + 100;
		$old_end = $hregions->{end};
		push( @$regions, $hregions );
		$i++;
	}
	$regions->[0]->{start} = 1;
	$regions->[-1]->{end}  = $chr->length + 1000;
	return $regions;
}

sub checkPatientNameInVcf {
	my $self = shift;
	die("\n\nERROR: FILE ".$self->file()." doesn't exists... Die...\n\n") unless (-e $self->file());
	if ($self->file() =~ /.gz/) { open (FILE, "zcat ".$self->file()." |"); }
	else { open (FILE, $self->file()); }
	while (<FILE>) {
		my $line = $_;
		chomp($line);
		next unless ($line =~ /#CHROM/);
		my @lTmp = split(' ', $line);
		if (scalar(@lTmp) > 10) {
			close(FILE);
			die("\n\nERROR: multiple patients in VCF file. Please defined you wanted patient name. Die.\n\n");
		}
		close(FILE);
		return $lTmp[-1];
	}
	close(FILE);
	return;
}

sub getThisChromosome {
	my ($self, $chr_name) = @_;
	return $self->hashChromosomes->{$chr_name} if (exists $self->hashChromosomes->{$chr_name});
	die("\n\nERROR: can't construct chromosome object with name $chr_name. Die.\n\n");	
}

sub parseVcfFile {
	my ($self, $reference) = @_;
	unless ($self->patient()) {
		unless ($self->noPatient()) {
			my $patientName = $self->checkPatientNameInVcf();
			my $hPat;
			$hPat->{id}   = $patientName;
			$hPat->{name} = $patientName;
			$self->patient($hPat);
			warn '-> Patient name found in VCF file: '.$self->patient->{name};
		}
	}
	my $hashVcf = $self->parseVcfFileForReference($reference);
	return $self->getVariantsObjects($hashVcf, $self->noPatient()) if ($self->getObjects());
	return $self->getHashDecompress($hashVcf);
}

sub getVariantsObjects {
	my ($self, $hashVcf, $is_noPatient) = @_;
	my $hashRes;
	foreach my $structType (keys %$hashVcf) {
		my @lObjs;
		while (my ($keyId, $valHash) = each(%{$hashVcf->{$structType}})) {
			$valHash = thaw(decompress($valHash));
			#my $h 
			#warn Dumper $valHash; die;
			#chromosome_object
			my $typeObj;
			if    ($structType eq 'snp') { $typeObj = 'GenBoVariation'; }
			elsif ($structType eq 'ins') { $typeObj = 'GenBoInsertion'; }
			elsif ($structType eq 'del') { $typeObj = 'GenBoDeletion'; }
			my $object = $typeObj->new($valHash);
			unless ($is_noPatient) {
				$object->patients_object->{$self->patient->{id}} = undef;
				$object->annex->{$self->patient->{id}} = $valHash->{annex}->{$self->patient->{id}};
			}
			my $chr = $self->getThisChromosome($object->chromosome_name());
			$chr->project( $self->project() );
			$object->chromosome_object( $chr );
			$object->project( $self->project() );
			$valHash = undef;
			delete $hashVcf->{$structType}->{$keyId};	
			push(@lObjs, $object);
		}
		$hashRes->{$structType} = \@lObjs;
	}
	return $hashRes;
}

sub getHashDecompress {
	my ($self, $hashVcf) = @_;
	my $hashRes;
	foreach my $type (keys %$hashVcf) {
		foreach my $id (keys %{$hashVcf->{$type}}) {
			my $h2 = thaw(decompress($hashVcf->{$type}->{$id}));
			if (exists $hashRes->{$type}->{$id}){
				my $h1 = thaw(decompress($hashRes->{$type}->{$id}));
				$h1->{annex}->{$self->patient->{name}}->{method_calling}->{$self->method()} = $h2->{annex}->{$self->patient->{name}}->{method_calling}->{$self->method()};
			}
			else {
				$hashRes->{$type}->{$id} = $h2;
			}
		}
	}
	return $hashRes;
}

1;