package GenBoCapture;

use Moo;

use Data::Dumper;
use Config::Std;
use FindBin qw($Bin);
use lib "$Bin/";
use lib "$Bin/../GenBoDB";
use Set::IntSpan::Fast::XS;
use Bio::DB::HTS::Tabix;
use Carp;
use Storable qw(store retrieve freeze thaw);
use List::Util qw( shuffle sum max min);
extends "GenBo";

has name => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		$self->infos()->{name};
		
	},
);

has type_object => (
	is		=> 'ro',
	default	=> "captures_object",
);

has isPcr =>(
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		return undef;
	},

);
has isCapture =>(
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		return 1;
	},

);
	
has gencode_version =>(
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->infos()->{gencode_version};
		
	},
);

has analyse => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		$self->infos()->{analyse};
		
	},
);
has validation_db => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		$self->infos()->{validation_db};
		
	},
);
has umi => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
				my $query = $self->getProject->buffer->getQuery();
		my $hash1 = $query->getUmiInfos($self->id);
		return $hash1;
		
	},
);


has hash_intspan => (
	is		=> 'ro',
	reader	=> 'getHashGenomicSpan',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		confess();
		my $hashIntSpan;
		my $listChrom = $self->getChromosomes();
		foreach my $chromObj (@{$self->getChromosomes}) {
			my $intSpan = $self->getIntSpanForChromosome($chromObj);
			$hashIntSpan->{$chromObj->id()} = $intSpan;
		}
		return $hashIntSpan;
	},
);

sub genomic_span {
		my ($self,$chr) = @_;
		die() unless $chr;
		
		return $self->{gs}->{$chr->name} if exists  $self->{gs}->{$chr->name};
		
		my $iter = $self->tabix->query($chr->name);
		unless ($iter) {
			$iter = $self->tabix->query($chr->ucsc_name);
		}
		$self->{gs}->{$chr->name} =  Set::IntSpan::Fast::XS->new();
		while (my $line = $iter->next){
			chomp($line);
			my($a,$b,$c,@d) = split(" ",$line);
			$self->{gs}->{$chr->name}->add_range($b,$c);
		}
		return $self->{gs}->{$chr->name};
}

has hash_intspan_extended => (
	is		=> 'ro',
	reader	=> 'getHashGenomicSpan_extended',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $hashIntSpan;
		my $listChrom = $self->getChromosomes();
		foreach my $chromObj (@{$self->getChromosomes}) {
			my $intSpan = $self->getIntSpanForChromosome_extended($chromObj);
			$hashIntSpan->{$chromObj->id()} = $intSpan;
		}
		return $hashIntSpan;
	},
);

has hash_intspan_referenceExtended => (
	is		=> 'ro',
	reader	=> 'getHashGenomicSpan_referenceExtended',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $hashIntSpan;
		my $listChrom = $self->getChromosomes();
		foreach my $chromObj (@{$self->getChromosomes}) {
			my $intSpan = $self->getIntSpanForChromosome_referenceExtended($chromObj);
			$hashIntSpan->{$chromObj->id()} = $intSpan;
		}
		return $hashIntSpan;
	},
);



sub parsePrimers{
	my ($self,$list) = @_;
	my $primers;
	foreach my $line (@$list){
		my $hpos;
		chomp($line);
		my (@toto) = split(" ",$line);
		my ($chrname,$startf,$endf,$startr,$endr,$plex) ;
		my $start_id;
		if (scalar(@toto) == 6){
		 ($chrname,$startf,$endf,$startr,$endr,$plex) = split(" ",$line);
		 $start_id = $startf;
		}
		elsif (scalar(@toto) == 4){ 
			 $chrname = $toto[0];
			 $startf =  $toto[1]-15;
			  $endf=  $toto[1];
			  $startr = $toto[2];
			  $endr = $toto[2]+15;
			  $plex = $toto[3];
			   $start_id = $startf;
			
		}
		elsif (scalar(@toto) == 3){ 
			 $chrname = $toto[0];
			 $startf =  $toto[1]-15;
			  $endf=  $toto[1];
			  $startr = $toto[2];
			  $endr = $toto[2]+15;
			  $plex = 1;
			   $start_id = $startf;
			
		}
		else{
			confess("you defined Pcr amplification method but I don't have enought column in primer file : $line ;")
		}
		my $chromosome = $self->getProject->getChromosome($chrname);
		$hpos->{chromosomes_object}->{$chromosome->id} =undef;
		$hpos->{gstart} = $startf;
		$hpos->{gend} = $endr;
		$hpos->{start} = $endf+1;
		$hpos->{end}   = $startr-1;
		$hpos->{start_forward}   = $startf;
		$hpos->{end_forward}   = $endf;
		$hpos->{id}="primer".$chrname."_$start_id";
		$hpos->{name}= $plex."_".$chrname."_$startf";
		$hpos->{multiplex}= $plex;
		$hpos->{length}   = abs($endf-$startr)+1;
		$hpos->{intspan} =  Set::IntSpan::Fast::XS->new($hpos->{start}."-".$hpos->{end});
		$hpos->{cnv} ={};
		$self->{primer_size}->{$hpos->{id}} = abs($endf-$startr)-1; 
		#$self->{primer_size}->{$hpos->{id}} = 1; 
		push(@$primers,$hpos);
		
	
	}
	return $primers;
} 
###### SET OBJECTS #####
sub parsePrimersForMultiplex {
	my ($self,$multiplex) = @_;
	my $data;
	foreach my $line (@{$self->primers_lines()}){
	
		
		push(@$data,$line);
	}
	
	my $primers = $self->parsePrimers($data);
	my $objs = $self->getProject()->flushObjects("primers",$primers);
	foreach my $o (@$objs){
		$o->{$self->type_object}->{$self->id} = undef;
	}
	return $objs;
}

sub getPrimersByMultiplex {
	my ($self,$multi) = @_;
	return $self->{multiplexes}->{$multi}->{primers} if exists $self->{multiplexes}->{$multi}->{primers};
	
	$self->{multiplexes}->{$multi}->{primers}  = $self->parsePrimersForMultiplex($multi);
	return $self->{multiplexes}->{$multi}->{primers} ;
}

sub setPatients {
	my $self = shift;
	my $project = $self->getProject();
	my %patientsId;
	foreach my $thisPatient (@{$project->getPatients()}) {
		if ($thisPatient->getCaptureId() eq $self->id()) {
			$patientsId{$thisPatient->id()} = undef;
		}
	}
	return \%patientsId;
}
has tabix => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $captureFileName = $self->gzFileName();
		my $captureTabixFileName = $captureFileName . '.tbi';
		if (not -e $captureTabixFileName) { die("\n\nERROR: $captureTabixFileName doesn't exist !! die...\n\n"); }
		return  Bio::DB::HTS::Tabix->new( filename => $captureFileName );
		
	},
);

has seqnames => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->tabix->seqnames();#$self->infos()->{name};
		
	},
);


sub setChromosomes {
	my $self = shift;
	my %chrIds;
	my $project = $self->getProject();
	confess();
	foreach my $name (@{self->seqnames}) {
		my $chr = $project->getChromosome($name);
		$chrIds{$chr->id} = undef;
	}
	return \%chrIds;
}

sub setVariants {
	my ($self, $methodName) = @_;
	my $hashVarIds;
	my $listPatients = $self->getPatients();
	my $hashIntSpan = $self->getHashGenomicSpan();
	my $listChromosomes = $self->getChromosomes();
	foreach my $chrom (@$listChromosomes) {
		my $chromIntSpan = $chrom->getGenomicSpan();
		foreach my $patient (@$listPatients) {
			my $listVar = $patient->$methodName($chrom->id(), $chromIntSpan);
			foreach my $var (@$listVar) {
				my $interIntSpan = $chromIntSpan->intersection($var->getGenomicSpan());
				unless($interIntSpan->is_empty()) { $hashVarIds->{$var->id()} = undef; }
			}
		}
	}
	return $hashVarIds;
}

sub setVariations {
	my $self = shift;
	return $self->setVariants('getVariations');
}

sub setDeletions {
	my $self = shift;
	return $self->setVariants('getDeletions');
}

sub setInsertions {
	my $self = shift;
	return $self->setVariants('getInsertions');
}



###### METHODS #####

sub contains {
	my ($self,$chr,$pos) = @_;
	my $intspan = $self->getIntSpanForChromosome($chr,100);
	return  $intspan->contains($pos);
}



sub getCaptureChromosomesName {
	my ($self) = @_;
	return $self->{chr_name} if exists  $self->{chr_name};
	my @l_tabix_chr_names = @{$self->seqnames};
	my @t;
	if ($l_tabix_chr_names[0] =~ /chr/) { @t = map {$self->project->getChromosome($_)->ucsc_name}  @l_tabix_chr_names; }
	else { @t = map {$self->project->getChromosome($_)->id}  @l_tabix_chr_names; }
	$self->{chr_name} = \@t;
	return $self->{chr_name};
}

sub getIntSpanForChromosome {
	my ($self, $GenBoChrom, $nbExtendNt) = @_;
	unless ($nbExtendNt) { $nbExtendNt = 0; }

	my $intSpan = Set::IntSpan::Fast::XS->new();
	my $chrName = $GenBoChrom->fasta_name();
	my $start = $GenBoChrom->start();
	my $end = $GenBoChrom->end();
	
	my ($find) = grep {$_ eq $GenBoChrom->name or $_ eq $GenBoChrom->ucsc_name} @{$self->getCaptureChromosomesName};
	return $intSpan unless $find;
	my $res = $self->tabix->query($find);
	
	while(my $line = $res->next()){
    	my @lFields = split("\t", $line);
    	my $start = int($lFields[1]) - $nbExtendNt;
    	$start =0 if $start <0;
    	my $end = int($lFields[2]) + $nbExtendNt;
    	$intSpan->add_range($start, $end);
  	}
  	$res->close();
	return $intSpan;
}

sub getIntSpanForChromosome_extended {
	my ($self, $GenBoChrom) = @_;
	return $self->getIntSpanForChromosome($GenBoChrom, 100);
}

sub getIntSpanForChromosome_referenceExtended {
	my ($self, $GenBoChrom) = @_;
	my $nbNt = 100_000;
	if ($self->getProject()->test()) { $nbNt = 2000; }
	return $self->getIntSpanForChromosome($GenBoChrom, $nbNt);
}

has infos =>(
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		warn "coucou";
		confess() unless $self->getProject;
		my $query = $self->getProject->buffer->getQuery();
		my $hash1 = $query->getCaptureInfos($self->id);
		warn Dumper $hash1;
		die();
		my $hash2 = $query->getCaptureTranscripts($self->id);
		
		my %newHash = (%$hash1, %$hash2);
		$self->{gencode_version} = $newHash{gencode_version};
		$self->{version} = $newHash{version};
		return \%newHash;
	}
);

has transcripts_name =>(
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		#my @ls =   keys %{$self->infos()->{transcripts}};
		my @ltrs;
		foreach my $tr_id (@{$self->infos()->{transcripts_name}}) {
			if ($tr_id =~ /_[0-9]+/) {
				push(@ltrs, $tr_id);
				next;
			}
			else {
				foreach my $chr (@{$self->getProject->getChromosomes()}) {
					my $tr_id_2 = $tr_id.'_'.$chr->id();
					if ($self->getProject->rocksGenBo->exists($tr_id_2)) {
						push(@ltrs, $tr_id_2);
						last;
					}
				}
			}
		}
		return \@ltrs;
		#return $self->infos()->{transcripts_name};
	}
);

sub return_bundle {
	my ($self,$tr) = @_;
	return $self->infos()->{transcripts}->{$tr};
}


#
has type =>(
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->infos()->{type};
	}
);

has version =>(
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->infos()->{version};
	}
);

has description => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->infos()->{description};
	}
);


has file_name => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->infos()->{filename};
	}
);
#https://www.polyweb.fr/NGS/NGS2023_6960/HG19_MT/align/dragen-align/CLE_CHR_6123GM005310_FLP.bam
has files =>(
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $project = $self->getProject();
		my $version = $project->getVersion();
		
		$version = "HG19" if $version =~/HG19/;
		my  $dir = $self->project->capture_dir;
		my $return = {};
		my $files;
		my $file = $dir . "/" . $self->type . "/" . $self->file_name . ".gz";
		confess("Unable to find capture file :" .$file."\n") unless -e $file;
		$files->{gz} = $file;
		my $file2 = $dir . "/" . $self->type . "/" . $self->file_name;
		$files->{bed} = $file2;
		#warn $files->{bed};
		$files->{hotspot} = $dir . "/hotspot/" .lc($self->analyse) . ".bed";;
		my $fp =  $dir . "/" . $self->type . "/" . $self->primers_filename;
		
	#	if (-s $fp){
			 
	#		$files->{primers} =  $dir . "/" . $self->type . "/" .$self->primers_filename;
	#		warn $files->{primers};
	#	}
	#	else {
			$files->{primers} = $files->{gz};
	#	}
		return $files;
	}
);

has dude_bed =>(
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $project = $self->getProject();
		my $version = $project->getVersion();
		
		$version = "HG19" if $version =~/HG19/;
		my  $dir = $self->project->capture_dir;;
		my $return = {};
		my $files;
		my $file = $dir . "/dude/dude.bed.gz";#. $self->type . "/" . $self->file_name . ".gz";
		confess("Unable to find capture file :" .$file."\n") unless -e $file;
		return $file;
	}
);

has gzFileName => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->files()->{gz}; 
	},
);
has primers_filename => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		$self->infos()->{primers_filename};
		
	},
);

has hotspots_filename => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		#return "/data-isilon/bipd-src/pnitschk/git2-polyweb/polygit/polyscripts/hotspot/cloves.bed"; 
		return $self->files()->{hotspot}; 
		
	},
); 

has hotspots => (
	is		=> 'ro',
	lazy	=> 1,
	
	default => sub {
		my $self = shift;
		my $file = $self->files()->{hotspot}; 
		return [] unless -e $self->hotspots_filename();
		open (PLEX,$self->hotspots_filename()) ;
		my $hotspots;
		#my $motif ={
#	chromosome=>'chr7',
#	start =>117232265,
#	end=>117232275,
#	seq=>'ACAAAAAAACA',
#} ;
 	my $h;
	while (my $line = <PLEX>){
		chomp($line);
		my @t = split(" ",$line);
		 my $id =$t[0].":".$t[1];
		 $h->{$id}->{gene} = $t[3];
		 my $gene = $self->project->newGene($t[3]);
		
		 ($h->{$id}->{ref},$h->{$id}->{alt}) = split("/",$t[4]);
		 $h->{$id}->{name} = $t[5];
		 if ($gene->strand == -1){
			$h->{$id}->{ref} = BioTools::complement($h->{$id}->{ref});
			$h->{$id}->{alt} = BioTools::complement($h->{$id}->{alt});
			$h->{$id}->{name} = "-".$t[5];
			
			}
		my $chr = $self->project->getChromosome($t[0],$t[0]);
		if (length($h->{$id}->{ref})==1 && length($h->{$id}->{alt})==1){
		$h->{$id}->{genbo_id} = $chr->name."_".($t[1]+1)."_".$h->{$id}->{ref}."_".$h->{$id}->{alt};
		}
		elsif(length($h->{$id}->{ref})>1){
			my $before = $chr->getSequence($t[1],$t[1]);
			$h->{$id}->{genbo_id} = $chr->name."_".($t[1]+1)."_".$before.$h->{$id}->{ref}."_".$before;
		}
		elsif(length($h->{$id}->{alt})>1){
			my $before = $chr->getSequence($t[1],$t[1]);
			$h->{$id}->{genbo_id} = $chr->name."_".($t[1]+1)."_".$before."_".$before.$h->{$id}->{alt};
		}
		else {
			confess();
		}
		$h->{$id}->{protid} = $t[7];
	}
	return $h;
	},
);
has multiplexFile => (
	is		=> 'ro',
	lazy	=> 1,
	
	default => sub {
		my $self = shift;
		return $self->files()->{primers}; 
		
	},
);

has multiplex =>(
	is		=> 'ro',
	lazy	=> 1,
	reader	=> 'getMultiplex',
	default => sub {
		my $self = shift;
		return[1];
		confess();

		my $primers = $self->getPrimers();

		my %t;
		foreach my $p (@$primers){
			$t{$p->multiplex} = undef;
		}
		my @tt = keys %t;
		return \@tt;
	},

);



has primers_lines =>(
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
#		warn "*********";
		if ($self->analyse eq "genome"){
			my @lines;
			foreach my $chr (@{$self->project->getChromosomes()}){
				
				my $from =1;
    			my $to = $chr->length;
  	
    	my $window =1000;
    	my $z =0;
    	
    	while ($from < $to){
    		$z++;
        my $start = $from;
        my $end = $from + $window;
        if ($end > $to){
            $end = $to;
        }
        push(@lines,$chr->ucsc_name."\t".$from."\t".$end."\n"."1000");
        $from = $end+1;
       }
	
		
			}
			return \@lines;	
		}
		
		my $capture_file = $self->gzFileName();
		die("unable to find $capture_file , I'm dyiiiiiinnnng." ) unless -e $capture_file;
			my @lines;
			if ($capture_file =~/\.gz/){
				my $bedtools = $self->buffer->software("bedtools");
				
		 		@lines = `zcat $capture_file | cut -f 1-3 | $bedtools makewindows -b /dev/stdin -w 500 `;
		 		confess() if $? ne 0;
		 		
			}
			else {
				confess();
			}
		chomp(@lines);
		return \@lines;
	},
);

has nb_primers => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		scalar(@{$self->primers_lines});	
	},
);

sub setListPrimers {
		my ($self,$multi) = @_;
		warn "cuicui";
		warn $self->nb_primers();
		die();
}
sub getListPrimers {
	my ($self,$multi) = @_;
	return $self->{list_primers} if exists $self->{list_primers};
		my $no = $self->project->noSqlCoverage();
		my $list = $no->get($self->project->name,$self->id);
		if ($list){
			$self->{list_primers} = $list;
			return $self->{list_primers};
		}
		if ($self->nb_primers() > 200000){
		 $self->{list_primers} =   $self->getRandomPrimers(500);
		}
		else {
	
			 $self->{list_primers} =   $self->getRandomPrimers(500);
		}
		
	return $self->{list_primers};
}

#sub getPrimersByMultiplex {
#	my ($self,$multi) = @_;
#	return $self->{multiplexes}->{$multi}->{primers} if exists $self->{multiplexes}->{$multi}->{primers};
#	my @primers = grep {$_->multiplex eq $multi } @{$self->getPrimers};
#	#$self->{multiplexes}->{$multi}->{primers} =  [@primers[0..100]];
#	$self->{multiplexes}->{$multi}->{primers}  = \@primers;
#	$self->{multiplexes}->{$multi}->{random_primers} =  [shuffle(@primers)];
#	
#	return $self->{multiplexes}->{$multi}->{primers} ;
#}




sub getRandomPrimers {
	my ($self,$nb) = @_;
	my $z = $self->nb_primers;
	

	my @t = shuffle @{$self->primers_lines()};
	my %rprimers;
	
	
	my $score =0;
	my $runs = $self->getProject->getRuns();
	
	my $zzz =0;	
	foreach my $p (@t){
		my ($chr,$start,$end,@s) = split(" ",$p);
		next if ($chr =~/[X,Y,M]/);
		$zzz++;
	
	#	$chr =~s/chr//;
	#	$chr = 'MT' if $chr eq 'M';
		my $id = "primer".$chr."_".($start-15);
		$self->{primer_size}->{$id} = abs($start-$end)-1; 
		
		$rprimers{$id} = $p; 
		last if $zzz>=$nb;
	
	}
	
	
	return \%rprimers;
}


sub primer_size{
	my($self,$id) = @_;
	return $self->{primer_size}->{$id}  if exists $self->{primer_size}->{$id}; 
	foreach my $t (@{$self->primers_lines()}){
			my ($chr,$start,$end,@s) = split(" ",$t);
			next if ($chr =~/[X,Y,M]/);
			
			#$zzz++;
	
	#	$chr =~s/chr//;
	#	$chr = 'MT' if $chr eq 'M';
			#next if $chr ne "chr5";
		my $id2 = "primer".$chr."_".($start-15);
		my $id3 = "primer".$chr."_".$start;
		#next unless $id2 =~/primerchr8/;
		#warn $id;
		$self->{primer_size}->{$id2} = abs($start-$end)-1;
		$self->{primer_size}->{$id3} = abs($start-$end)-1; 
	#	warn  abs($start-$end)-1; 
	}
	
	return $self->{primer_size}->{$id}  if exists $self->{primer_size}->{$id};
	confess($self->name." ".$id);
}


sub count_multiplex {
	my ($self,$multi,$gene_id,$patient) = @_;

	my $run_id = $patient->getRun->id;
	return $self->{count}->{$run_id} if exists $self->{count}->{$run_id};
	#return $self->{multiplex}->{$multi}->{$gene_id}->{sum} if exists  $self->{multiplex}->{$multi}->{$gene_id}->{sum};
	my $primers = $self->getListPrimers($multi,$gene_id);
	$self->{count}->{$run_id} = 0;
	my $nn =0;
	 my @patients = grep {$_->getRun->id eq $run_id} @{$self->getProject->getPatients()};
		my @data;
		my $nbp =0;
		foreach my $p (@patients){
			next if $p->getRun->id ne $run_id;
			$nbp ++;
		}

	my $nbh = scalar(keys %$primers );
	my $z =0;
		foreach my $primer (keys %$primers){
			my $sum1;
				$nbp =0;
			$z ++;
			foreach my $p (@patients){
			next if $p->getRun->id ne $run_id;
			my $count =  $p->count($primer,$self);
			$nbp++;
		#	last if $nbp > 20;
			
			$self->{count}->{$run_id}  +=  $count;
			push(@data,$count);
			
		}
	}
	 my $sum = sum @data;
	 if (@data>10){
	 my $min = min(@data);
	 my $max = max(@data);
	 $sum = $sum -($min+$max);
	 }
	 $self->{count}->{$run_id} = $sum;
		return $self->{count}->{$run_id}  ;
		
}


has count_all_multiplex => (
	is		=> 'ro',
	lazy=>1,
	default => sub {
		my $self = shift;
		my $capture = $self->getProject->getCapture();
		my $primers = $capture->getPrimersByMultiplex($self->multiplex());
		my $sum =0;
		#my $t1 = $self->getTranscripts->[0];
		
		foreach my $primer (@$primers){
			next if $primer->getChromosome()->name eq "X";
			#my $t2 = $primer->getTranscripts->[0];
			#next if $t2 && $t1 && $t2->id eq $t1->id; 
			$sum += $primer->count_all();
		}
		return $sum;
	} 
);

sub restore_save_primers {
	my ($self,$no,$chromosome) = @_;
	
		my $objs = $no->get("primers",$chromosome->name);
			foreach my $o (@$objs){
					$self->project->{objects}->{primers}->{$o->{id}} = $o;
					$o->{project} =  $self->project;
					$o->{buffer} = $self->buffer;
			}
		return $objs;
}


has intspan_exon => (
	is		=> 'ro',
	lazy=>1,
	default => sub {
		my $self = shift;
		warn "start";
		my $transcripts = $self->project->getListTranscripts();
		
		warn scalar(@$transcripts);
		my $nb;
		my $intspan;
		foreach my $tr (@$transcripts){
			$nb ++;
			warn $nb if $nb % 1000 ==0;
			my $t = $self->project->newTranscript($tr);
			foreach my $e (@{$t->getExons()}){
				next  if $e->intspan_no_utr->is_empty;
				
				$intspan->{$e->getChromosome->name} = Set::IntSpan::Fast::XS->new() unless exists $intspan->{$e->getChromosome->name};
				$intspan->{$e->getChromosome->name}->add_range($e->start_utr,$e->end_utr);
			}
		}
		return $intspan;
	} 
);

sub constructPrimersForExome {
		my ($self,$chromosome) = @_;
		my $file = $self->dude_bed;
		my $tabix = Bio::DB::HTS::Tabix->new( filename => $file );
		my $iter = $tabix->query($chromosome->name);
		#warn $file;
		my $hreturn;
		my $primers;
		while (my $line = $iter->next){
				chomp($line);
				my($chr,$s,$e) = split("\t",$line);
				my $hpos;
				$hpos->{chromosomes_object}->{$chromosome->id} =undef;
				my ($chrname,$startf,$endf,$startr,$endr,$plex) ;
				 $chrname = $chromosome->name();
				my $start_id;
				$start_id = $s;
				$startf =  $start_id-15;
				$startf =1 if $startf <= 0;
			 	$endf=  $s;
			 	$startr =  $e;
			 	$endr = $e+15;
			 	
			 	$plex = 1;
				$hpos->{gstart} = $startf;
				$hpos->{gend} = $endr;
				$hpos->{start} = $endf-1;
				$hpos->{start} = 1 if $hpos->{start} <=0; 
				$hpos->{end}   = $startr+1;
				$hpos->{start_forward}   = $startf;
				$hpos->{end_forward}   = $endf;
				$hpos->{id}="primer".$chrname."_$startf";
				
				$hpos->{id2}="primer".$chrname."_$start_id";
				$hpos->{name}= $plex."_".$chrname."_$startf";
				$hpos->{multiplex}= $plex;
				$hpos->{length}   = ($endf-$startr)+1;
				#warn $hpos->{start}."-".$hpos->{end};
				$hpos->{intspan} =  Set::IntSpan::Fast::XS->new($hpos->{start}."-".$hpos->{end} );
				$hpos->{cnv} ={};
				$self->{primer_size}->{$hpos->{id}} = abs($startf-$startr)-1; 
				$hpos->{$self->type_object}->{$self->id} = undef;
				$hreturn->{$hpos->{id}} = undef;
#				warn $hpos->{start}." ".$hpos->{end} if $hpos->{id} eq "primer1_6854";;
				push(@$primers,$hpos);
			}
	$iter->close();	
	#warn scalar(@$primers);	
	my $objs = $self->getProject()->flushObjects("primers",$primers);
	
	foreach my $o (@$objs){
		$o->{$self->type_object}->{$self->id} = undef;
	}
	return $objs;
}
sub constructPrimersForGenome {
		my ($self,$chromosome) = @_;
		my $intspan = $chromosome->project->liteIntervalTree->get_intspan("genes_padding",$chromosome->name);
		
		#my $intspan = $chromosome->project->liteIntervalTree->get_intspan("transcripts_padding",$chromosome->name);
		my $iter = $intspan->iterate_runs();
	
		my $hreturn;
			my $primers;
		while (my ( $from, $to ) = $iter->()) {
			my $intervals = $self->buffer->divide_by_chunks($from,$to,500);
			#warn $from." ".$to;
			#warn Dumper $intervals;
			#die();
			foreach my $interval (@$intervals){
				my $hpos;
				$hpos->{chromosomes_object}->{$chromosome->id} =undef;
				my ($chrname,$startf,$endf,$startr,$endr,$plex) ;
				 $chrname = $chromosome->name();
				my $start_id;
				$start_id = $interval->[0];
				$startf =  $start_id-15;
				$startf =1 if $startf <= 0;
			 	$endf=  $start_id;
			 	$startr =  $interval->[1];
			 	$endr =$interval->[1]+15;
			 	
			 	$plex = 1;
				$hpos->{gstart} = $startf;
				$hpos->{gend} = $endr;
				$hpos->{start} = $endf-1;
				$hpos->{start} = 1 if $hpos->{start} <=0; 
				$hpos->{end}   = $startr+1;
				$hpos->{start_forward}   = $startf;
				$hpos->{end_forward}   = $endf;
				$hpos->{id}="primer".$chrname."_$startf";
				
				$hpos->{id2}="primer".$chrname."_$start_id";
				$hpos->{name}= $plex."_".$chrname."_$startf";
				$hpos->{multiplex}= $plex;
				$hpos->{length}   = ($endf-$startr)+1;
				#warn $hpos->{start}."-".$hpos->{end};
				$hpos->{intspan} =  Set::IntSpan::Fast::XS->new($hpos->{start}."-".$hpos->{end} );
				$hpos->{cnv} ={};
				$self->{primer_size}->{$hpos->{id}} = abs($startf-$startr)-1; 
				$hpos->{$self->type_object}->{$self->id} = undef;
				$hreturn->{$hpos->{id}} = undef;
#				warn $hpos->{start}." ".$hpos->{end} if $hpos->{id} eq "primer1_6854";;
				push(@$primers,$hpos);
			}
		}
		
	#warn scalar(@$primers);	
	my $objs = $self->getProject()->flushObjects("primers",$primers);
	foreach my $o (@$objs){
		warn $o->start." ".$o->end if  $o->id eq "primer1_6854";
		last  if  $o->id eq "primer1_6854";
	}
	#die();
	foreach my $o (@$objs){
		$o->{$self->type_object}->{$self->id} = undef;
	}
	return $objs;
}

#sub constructPrimersForGenome {
#		my ($self,$chromosome) = @_;
#		my $regions = $chromosome->chunk(1000);
#		my $primers;
#			my $hreturn;
#		foreach my $region (@$regions){
#			my $v = $chromosome->genesIntervalTree->fetch($region->{start},$region->{end});
#			next unless $v;
#				my $hpos;
#				$hpos->{chromosomes_object}->{$chromosome->id} =undef;
#				 my ($chrname,$startf,$endf,$startr,$endr,$plex) ;
#				  $chrname = $chromosome->name();
#				my $start_id;
#				  $start_id = $region->{start};
#					$startf =  $start_id-15;
#			 		$endf=  $start_id;
#			 		$startr = $region->{end};
#			 		$endr =$region->{end}+15;
#			 		$plex = 1;
#				$hpos->{gstart} = $startf;
#				$hpos->{gend} = $endr;
#				$hpos->{start} = $endf+1;
#				$hpos->{end}   = $startr-1;
#				$hpos->{start_forward}   = $startf;
#				$hpos->{end_forward}   = $endf;
#				$hpos->{id}="primer".$chrname."_$startf";
#				$hpos->{id2}="primer".$chrname."_$start_id";
#				$hpos->{name}= $plex."_".$chrname."_$startf";
#				$hpos->{multiplex}= $plex;
#				$hpos->{length}   = ($endf-$startr)+1;
#				$hpos->{intspan} =  Set::IntSpan::Fast::XS->new($hpos->{start}."-".$hpos->{end} );
#				$hpos->{cnv} ={};
#				$self->{primer_size}->{$hpos->{id}} = abs($startf-$startr)-1; 
#				$hpos->{$self->type_object}->{$self->id} = undef;
#				$hreturn->{$hpos->{id}} = undef;
#				push(@$primers,$hpos);
#		}
#		warn scalar(@$primers);	
#		my $objs = $self->getProject()->flushObjects("primers",$primers);
#	foreach my $o (@$objs){
#		$o->{$self->type_object}->{$self->id} = undef;
#	}
#	return $objs;
#}

sub parsePrimersForChromosome {
	my ($self,$chromosome) = @_;

	my $project = $self->project;
	my $print = $self->project->print_waiting();
	if ($self->analyse eq "genome"){
		return $self->constructPrimersForGenome($chromosome);
	}
##	if ($self->analyse eq "exome"){
#		return $self->constructPrimersForExome($chromosome);
#	}
	my $objs;
	my $multiplex;
	my $all ;
	my $primers;
	my $n = $chromosome->name;
	my $u = $chromosome->ucsc_name;
	my $hreturn;
	my @lines = grep{$_=~/^$n|^$u\s/} @{$self->primers_lines()};
	my $intspan =Set::IntSpan::Fast->new();
	my $nb =1;
foreach my $line (@lines){
	$nb ++;
		my $hpos;
		chomp($line);
		my (@toto) = split(" ",$line);
			
		#die($toto[1]) if $toto[1] eq 70181901;
		if (abs($toto[1]-$toto[2]) < 5 ){
				$toto[1] -= 15;
				$toto[2] +=15;
				#warn "coucou";
			}
			
		#	die($toto[1]) if $toto[1] eq 70181901;
			my ($chrname,$startf,$endf,$startr,$endr,$plex) ;
			my $start_id;
			 $chrname = $toto[0];
			 $start_id = $toto[1];
			 $startf =  $toto[1]-15;
			 $endf=  $toto[1];
			 $startr = $toto[2];
			 $endr = $toto[2]+15;
			 $plex = 1;
		$hpos->{chromosomes_object}->{$chromosome->id} =undef;
		print "." if $print && $nb%2000==0;
		$hpos->{gstart} = $startf;
		$hpos->{gend} = $endr;
		$hpos->{start} = $endf+1;
		$hpos->{end}   = $startr-1;
		confess() if ($hpos->{start} > $hpos->{end});
		($hpos->{end}, $hpos->{start}) = ($hpos->{start},$hpos->{end})  if ($hpos->{start} > $hpos->{end});
		$hpos->{start_forward}   = $startf;
		$hpos->{end_forward}   = $endf;
		$hpos->{id}="primer".$chrname."_$startf";
		$hpos->{id2}="primer".$chrname."_$start_id";
		$hpos->{name}= $plex."_".$chrname."_$startf";
		$hpos->{multiplex}= $plex;
		$hpos->{length}   = ($endf-$startr)+1;
		$hpos->{intspan} =  Set::IntSpan::Fast::XS->new($hpos->{start}."-".$hpos->{end} );
		
		$hpos->{cnv} ={};
		$self->{primer_size}->{$hpos->{id}} = abs($startf-$startr)-1; 
		$hpos->{$self->type_object}->{$self->id} = undef;
		
		$hreturn->{$hpos->{id}} = undef;
		push(@$primers,$hpos);
		
		
	}
	print "@" if $print;
	 $objs = $self->getProject()->flushObjects("primers",$primers);
	foreach my $o (@$objs){
		$o->{$self->type_object}->{$self->id} = undef;
	}
	return $objs;
}
sub setPrimersForOneChromosome {
	
	
}
sub setPrimers {
	my ($self) = @_;
#	die();
	my $chromosomes = $self->project->getChromosomes();
	my %hchrs;
	my @objs;
	my $no = $self->project->noSqlCoverage();
	foreach my $chr (@$chromosomes){
	
		push(@objs,@{$self->parsePrimersForChromosome($chr)});
	}
	my %hash;
	  map{$hash{$_->id}++} @objs;
	return \%hash;;
	
}


has bedFileName => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->files->{bed}; 
	},
);

has dir_capture_controls_dude => (
	is		=> 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		my $dir = $self->buffer->public_data_root."/".$self->project->annotation_genome_version."/dude/NOVASEQ/".$self->name;
		#$dir ="/data-xfs/dev/pnitschk/svn-genbo/polypipeline/scripts/scripts_pipeline/dude/genome/NOVASEQ//genome_hg19_cng/";
		confess("no dude for cpature :".$self->name." ".$dir) unless -e $dir;
		return $dir;
	},
);
has sd_controls_dude => (
	is		=> 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		 
	 	$self->{control_nosql} =  GenBoBinaryFile->new(name=>"controls_sd",dir=>$self->dir_capture_controls_dude,mode=>"r");
		return $self->{control_nosql};
	},
);

has depth_controls_dude => (
	is		=> 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		
	 	$self->{control_nosql} =  GenBoBinaryFile->new(name=>"controls_depth",dir=>$self->dir_capture_controls_dude,mode=>"r");
		return $self->{control_nosql};
	},
);

has sd_controls_dude_dev => (
	is		=> 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		
	 	$self->{control_nosql} =  GenBoBinaryFile->new(name=>"controls_sd",dir=>"/data-xfs/dev/pnitschk/svn-genbo/polypipeline/scripts/scripts_pipeline/dude/genome/NOVASEQ/genome_hg19_cng",mode=>"r");
		return $self->{control_nosql};
	},
);



has ratio_controls_dude_dev => (
	is		=> 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		return  GenBoBinaryFile->new(name=>"controls_ratio_16",dir=>"/data-xfs/dev/pnitschk/svn-genbo/polypipeline/scripts/scripts_pipeline/dude/genome/NOVASEQ/genome_hg19_cng",mode=>"r");
	},
);

has infos_controls_dude => (
	is		=> 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		#my $lmdb2 = GenBoNoSqlLmdb->new(dir=>$dir_control,mode=>"r",name=>"",is_compress=>1);
		return  GenBoNoSqlLmdb->new(name=>"control.infos",dir=>$self->dir_capture_controls_dude,mode=>"r",is_compress=>1);
	},
);
has ratio_controls_dude => (
	is		=> 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		return  GenBoBinaryFile->new(name=>"controls_ratio",dir=>$self->dir_capture_controls_dude,mode=>"r");
	},
);
has ratio_controls_dude_new => (
	is		=> 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		return  GenBoBinaryFile->new(name=>"controls_ratio.new",dir=>$self->dir_capture_controls_dude,mode=>"r");
	},
);


1;
