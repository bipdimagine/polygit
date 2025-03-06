package GenBoProjectTests;
use strict;
use FindBin qw($Bin);
use Storable qw(retrieve);
use Moo;
use Data::Dumper;
use GenBoProject;
use GenBoProjectCache;
use Storable qw(store retrieve freeze dclone thaw);

extends 'GenBoProject';




has path_tmp => (
	is		 => 'ro',
	required => 1,
);

has test => (
	is 		=> 'ro',
	default	=> 1,
);

has release => (
	is		=> 'rw',
	default => 'HG38',
);

has test_name => (
	is		=> 'rw',
	lazy    => 1,
	default => undef,
);

has name => (
	is		=> 'rw',
	lazy    => 1,
	default => undef,
);

has name_tests_genbo_annovar  => (
	is 		=> 'ro',
	lazy    => 1,
	default	=> 'NGSTEST_FREQ',
);

has path_annovar => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return '/software/distrib/annovar/latest/';
	},
);

has project_root_path => (
	is		=> 'rw',
	lazy 	=> 1,
	reader 	=> 'getProjectRootPath',
	default => sub {
		my $self = shift;
		#my $path = $Bin.'/datas/tests_frequences/'.$self->name().'/';
		my $path = $Bin.'/../../../../PolyDatasTests/datas/compare_genbo_annovar/'.$self->name().'/';
		return $path;
	},
);

has year => (
	is		=> 'rw',
	lazy 	=> 1,
	default => '2025',
);

has in_this_run_patients => (
	is => 'rw',
	default => sub {
		my $self = shift;
		my $h;
		$h->{nb_patients} = 1;
		return $h;
	},
);

has annotation_version => (
	is => 'rw',
	default => undef,
);

has infosProject => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $res = {
          'version' => 'HG38',
          'name' => 'NGSTEST_FREQ',
          'description' => 'Project Tests',
          'dbname' => 'Polyexome',
          'somatic' => '0',
          'creation_date' => '2017-05-07 20:00:00',
          'projectTypeId' => '3',
          'projectType' => 'ngs',
          'id' => '1',
          'dejavu' => '1'
        };
		return $res;
	},
);



##### METHODES GLOBALES #####



sub setPatients {
	my $self = shift;
	my $res;
	if ($self->test_name() eq 'tests_frequences') {
		$res = $self->setPatients_tests_genbo_annovar();
	}
	elsif ($self->test_name() eq 'tests_annotations') {
		$res = $self->setPatients_tests_genbo_annovar();
	}
	my %names;
	foreach my $h (@{$res}){
		$h->{id} = $h->{patient_id};
		$names{$h->{id}} = undef;
		$h->{project} = $self;
		next if exists $self->{objects}->{patients}->{$h->{id}};
		$self->{objects}->{patients}->{$h->{id}} = new GenBoPatient($h);
	}
	return \%names;
}

sub setPatients_tests_genbo_annovar {
	my $self = shift;
	my @lH;
	my $h = {
            'status' => '2',
            'flowcell' => 'A',
            'genbo_id' => '0',
            'patient_id' => '28865',
            'capture_id' => '146',
            'project_id' => '2724',
            'name' => 'patient',
            'origin' => 'patient',
            'mother' => 'mother',
            'sex' => '2',
            'description' => undef,
            'run_id' => '1462',
            'bar_code' => 'AAGGTACA',
            'creation_date' => '2017-02-23 10:40:23',
            'father' => 'father',
            'project_id_dest' => '2724',
            'family' => 'fam'
          };
	push(@lH, $h);
	return \@lH;
}

sub init_genbo_annovar {
	my $self = shift;
	my $project_name = $self->name_tests_genbo_annovar();
	$self->name( $project_name );
	my $obj_patient = $self->getPatient('patient');
	my @callingMethods;
	push (@callingMethods, 'haplotypecaller');
	$obj_patient->callingMethods(\@callingMethods);
	return $obj_patient;
}

sub hash_ccds {
	my $ccds_file_path = '/data-xfs/dev/qbrandao/mart_ccds.txt';
	my $ccds_trans;
	open (CCDS, "<$ccds_file_path");
	foreach my $ligne (<CCDS>) {
		next if ($ligne =~ /Gene stable ID/); 	# ici on passe l'en-tete du fichier
		my @lCCDS = split (",", $ligne);		# on split le csv
		my $ENSG = $lCCDS[0];					# le premier element est l'ENSG
		my $ENST = $lCCDS[1];					# le deuxieme element est l'ENST
		my $CCDS = $lCCDS[2];					# le troisieme element est l'ID CCDS
		if (length($CCDS) != 1) {				# Si pas d'id ccds, le troisieme element contient un caractere
			$ccds_trans->{$ENST}++;				# On ajoute au hash que les ENST qui possedent un id ccds
		}
	}
	print "\nParsing of CCDS file Complete !\n";
	close (CCDS);
	return $ccds_trans;
}


##### TESTS ANNOTATIONS ######


sub launch_tests_annotations {
	my $self = shift;
	$self->test_name('tests_annotations');
	my ($obj_patient) = $self->init_genbo_annovar();
	print "\n##### TEST ANNOTATIONS #####\n\n";
	my $ccds_trans = hash_ccds();
	print "\n### LAUNCHING annot_genbo\n";
	$self->t_annotation_annot_genbo ( $obj_patient);
	print "\n### LAUNCHING annot_annovar\n";
	$self->t_annotation_annot_annovar ($obj_patient);
	print "\n### LAUNCHING parsing_annovar\n";
	$self->t_annotation_annovar_multianno($ccds_trans);
	$self->t_annotation_annovar_exonic($ccds_trans);
	print "\n### LAUNCHING comparaison\n";
	my $path_output = $self->path_tmp() . '/t_annotation';
	$self->t_annotation_compar ($path_output, 1, 1);
	$self->t_annotation_compar ($path_output, 1, 0);
	$self->t_annotation_compar ($path_output, 0, 1);
	
	#TODO:
	
	
	print "\n### LAUNCHING comptage\n";
	$self->t_annotation_comptage ($path_output);
}

sub t_annotation_comptage {
	my ($self, $path_output ) = @_;
	### Fonction de comptage permettant d'obtenir un resume du pipeline
	open( RES, ">$path_output.comptage" );
	open( COM, "<$path_output.common" );
	open( OTH, "<$path_output.others" );
	
	### <-- Partie de comparaison des features obtenues --> ###
	
	my ( $l, $t, $n, $r, $enst, $c, $o, $nb_ok, $nb_not, $nb_total );
	my $h;
	
	foreach my $ligne (<COM>) {
		chomp($ligne);
		my @tab = split( "\t", $ligne );
		if ($tab[11] eq 'same_trans') { #on compte les same_trans
			$t++;
		}
		if ($tab[11] eq 'same_annot') { #on compte les same_annot
			$n++;
		}
		if ($tab[4] eq $tab[7]) { #on compte les same_rs
			$r++;
		}
		$l++;   #compte le nombre de variants, sachant qu'il n'y a pas de header
		$c++;	#compte le nombre de ligne du .common
		###on met les infos du .common dans un hash
		$h->{"common"}->{$tab[0]}++;
		$h->{"annot_ok"}->{"$tab[5];$tab[8]::$tab[10]"}->{'ligne'}++;
		if (not exists $h->{"annot_ok"}->{"$tab[5];$tab[8]::$tab[10]"}->{'variants'}->{$tab[0]}) {
			$h->{"annot_ok"}->{"$tab[5];$tab[8]::$tab[10]"}->{'variants'}->{$tab[0]}++;
			$h->{"annot_ok"}->{"$tab[5];$tab[8]::$tab[10]"}->{'nb_var'}++;
			$nb_ok++;
			$nb_total++;
		}
	}
	foreach my $ligne (<OTH>) {
		chomp($ligne);
		my @tab = split( "\t", $ligne );
		if ($tab[4] eq $tab[7]) {
			$r++;
		}
		$l++;   #compte le nombre de variants, sachant qu'il n'y a pas de header
		$o++;	#compte le nombre d ligne du .others
		###on met les infos du .others dans un hash
		$h->{"others"}->{$tab[0]}++;
		$h->{"annot_non_ok"}->{"$tab[5];$tab[8]::$tab[10]"}->{'ligne'}++;
		if (not exists $h->{"annot_non_ok"}->{"$tab[5];$tab[8]::$tab[10]"}->{'variants'}->{$tab[0]}) {
			$h->{"annot_non_ok"}->{"$tab[5];$tab[8]::$tab[10]"}->{'variants'}->{$tab[0]}++;
			$h->{"annot_non_ok"}->{"$tab[5];$tab[8]::$tab[10]"}->{'nb_var'}++;
			$nb_not++;
			$nb_total++;
		}
	}
	###on print les infos du comptage general dans [features]
	print ("[features]\n");
	print ("Total_variants:$nb_total\n");
	print ("Total_common:$nb_ok\n");
	print ("Total_others:$nb_not\n");
	print ("Total_lines_variants:$l\n");
	print ("Total_lines_common:$c\n");
	print ("Total_lines_others:$o\n");
	print ("Same_transcripts:$t\n");
	print ("Same_annotation:$n\n");
	print ("Same_rs:$r\n");
	print ("\n");
	print ("[common]\n");
	foreach my $k ( reverse sort { $h->{"annot_ok"}->{$a}->{'nb_var'} <=> $h->{"annot_ok"}->{$b}->{'nb_var'} } keys %{$h->{"annot_ok"}}) {
		print  $k.":".$h->{"annot_ok"}->{$k}->{'ligne'}."lines".";".$h->{"annot_ok"}->{$k}->{'nb_var'}."variants"."\n";
	}
	print ("\n");
	print ("[others]\n");
	foreach my $k ( reverse sort { $h->{"annot_non_ok"}->{$a}->{'nb_var'} <=> $h->{"annot_non_ok"}->{$b}->{'nb_var'} } keys %{$h->{"annot_non_ok"}}) {
		print  $k.":".$h->{"annot_non_ok"}->{$k}->{'ligne'}."lines".";".$h->{"annot_non_ok"}->{$k}->{'nb_var'}."variants"."\n";
	}
	print  ("\n");
	my $com_tot = scalar (keys %{$h->{"common"}});
	my $oth_tot = scalar (keys %{$h->{"others"}});
	my $tot = $com_tot + $oth_tot;
	my $var_com = $com_tot * 100 / $tot;
	my $var_oth = $oth_tot * 100 / $tot;
	print 'VAR COMMON: '.($com_tot)."\t-> ".sprintf("%.2f",$var_com)."%\n";
	print 'VAR OTHERS: '.($oth_tot)."\t\t-> ".sprintf("%.2f",$var_oth)."%\n";
	close (COM);
	close (OTH);
	close(RES);
}

sub t_annotation_parser {
	my ($ligne) = @_;
	### Fonction permettant le parsing ds fichiers formates
	my @liste = split( '\t', $ligne );
	my $chr = $liste[0];
	$chr =~ s/chr//;
	my $pos            = $liste[1];
	my $ref            = $liste[2];
	my $alt            = $liste[3];
	my $gene_name      = $liste[4];
	my $trans_name     = $liste[5];
	my $annot          = $liste[6];
	my $rs_name   	   = '';
	if ($liste[7]) { $rs_name = $liste[7]; }
	my $var_id     = "$chr\_$pos\_$ref\_$alt";
	if ( $rs_name eq ' ' or $rs_name eq '.' ) {
		$rs_name = undef;
	}
	return (
		$chr,        $pos,   $ref,    $alt, $gene_name,
		$trans_name, $annot, $var_id, $rs_name
	);
}

sub t_annotation_compar {
	### ici on compare les resultats de genbo et d'annovar
	my ($self, $path_output, $levier_genbo, $levier_annovar ) = @_;
	my $projectName = $self->name();
	my $path_tmp = $self->path_tmp();
	my $genbo    = "$path_tmp/$projectName\_tab-annotation.tsv";
	my $annovar1 = "$path_tmp/$projectName.hg19_multianno.txt.tsv";
	my $annovar2 = "$path_tmp/$projectName.avinput.exonic_variant_function.tsv";
	###On commence par genbo
	open( TSV1, "<$genbo" ) or die "$!";
	my $hash;
	foreach my $ligne1 (<TSV1>) {
		chomp($ligne1);
		next if ( $ligne1 =~ /annotation/ );
		my (
			$chr1,   $pos1,       $ref1,
			$alt1,   $gene_name1, $trans_name1,
			$annot1, $id1,        $rs_name1
		) = t_annotation_parser($ligne1); #on utilise le parser de fichiers formates
		#on remplit le hash 
		$hash->{$chr1}->{$id1}->{'gene'}->{$gene_name1}->{'GenBo'}
		  ->{'annotation'}->{$trans_name1}->{$annot1} = undef;
		$hash->{$chr1}->{$id1}->{'ref'} = $ref1;
		$hash->{$chr1}->{$id1}->{'alt'} = $alt1;
		$hash->{$chr1}->{$id1}->{'gene'}->{$gene_name1}->{'GenBo'}->{'rs'} =
		  $rs_name1;
	}
	close(TSV1);
	
	### pour chacun des deux fichiers d'annovar
	foreach my $annovar ( $annovar1, $annovar2 ) {
		open( TSV2, "<$annovar" ) or die "$!";
		foreach my $ligne2 (<TSV2>) {
			chomp($ligne2);
			my (
				$chr2,   $pos2,       $ref2,
				$alt2,   $gene_name2, $trans_name2,
				$annot2, $id2,        $rs_name2
			) = t_annotation_parser($ligne2);
			if ($annot2) {
				my @lTmp = split( ',', $annot2 );
				$annot2 = $lTmp[0];
			}
			else { $annot2 = ''; } 

			$hash->{$chr2}->{$id2}->{'gene'}->{$gene_name2}->{'Annovar'}->{'annotation'}->{$trans_name2}->{$annot2} = undef;
			$hash->{$chr2}->{$id2}->{'ref'} = $ref2;
			$hash->{$chr2}->{$id2}->{'alt'} = $alt2;
			if ($rs_name2) {
				$hash->{$chr2}->{$id2}->{'gene'}->{$gene_name2}->{'Annovar'}->{'rs'} = $rs_name2;  
			}
		}
		close(TSV2);
	}

	my $hAnnot;
	my $total;
	### Table d'association des annotations annovar en genbo
	my $hAnnot_keys = {
		'UTR3'              	 => 'utr',
		'ncRNA_exonic'       	 => 'ncrna,exonic',
		'upstream,downstream'	 => 'intronic',
		'stoploss'           	 => 'phase',
		'splicing'           	 => 'splicing',
		'UTR3,UTR5'         	 => 'utr',
		'exonic,splicing'        => 'splicing,coding',
		'upstream'               => 'intronic',
		'frameshift deletion'    => 'frameshift',
		'frameshift_deletion'    => 'frameshift',
		'nonframeshift deletion' => 'non-frameshift',
		'nonframeshift_deletion' => 'non-frameshift',
		'exonic'                 => 'coding',
		'nonsynonymous SNV'      => 'coding',
		'nonsynonymous_SNV'      => 'coding',
		'intergenic'             => 'intergenic',
		'frameshift insertion'   => 'frameshift',
		'frameshift_insertion'   => 'frameshift',
		'synonymous SNV'         => 'silent',
		'synonymous_SNV'         => 'silent',
		'ncRNA_intronic'         => 'ncrna,intronic',
		'stopgain'               => 'stop',
		'ncRNA_splicing'         => 'ncrna,splicing',
		'ncRNA_exonic,splicing'  => 'ncrna,coding,splicing',
		'UTR5'                   => 'utr',
		'downstream'             => 'intronic',
		'intronic'               => 'intronic'
	};
	### On choisit de separer les variants en deux categories, les cas ou les infos sont identiques ou semblables, et les cas ou les infos sont discordantes
	open( COM, ">$path_output.common" ); #pour les infos identiques ou semblables
	open( DIF, ">$path_output.others" ); #pour les infos discordantes
	foreach my $chr ( 1 .. 22, 'X', 'Y', 'MT' ) { #on parcourt tous les chromosomes
		foreach my $id ( keys( %{ $hash->{$chr} } ) ) { #pour chaque variant
			my $ref = $hash->{$chr}->{$id}->{'ref'};
			my $alt = $hash->{$chr}->{$id}->{'alt'};
			my ( $is_common_trans, @lTmp, $new_annot_annovar );
			my $is_common_annot = 0;
			foreach my $gene ( keys( %{ $hash->{$chr}->{$id}->{'gene'} } ) ) { #pour chaque gene de chaque variant
				foreach my $trans_genbo (keys(%{$hash->{$chr}->{$id}->{'gene'}->{$gene}->{'GenBo'}->{'annotation'}})) { #pour chaque transcrit genbo de chaque gene de chaque variant
					my ( $rs_genbo, $annot_genbo );
					next if ( $levier_genbo and not exists $hash->{$chr}->{$id}->{'gene'}->{$gene}->{'GenBo'} );
					next  if ( $levier_annovar and not exists $hash->{$chr}->{$id}->{'gene'}->{$gene}->{'Annovar'} );
					if ( exists $hash->{$chr}->{$id}->{'gene'}->{$gene} ->{'GenBo'} ) { #si les infos de genbo existent
						$rs_genbo = $hash->{$chr}->{$id}->{'gene'}->{$gene}->{'GenBo'}->{'rs'}; #on prend le rs
						if ( not $rs_genbo ) { $rs_genbo = '.'; } #si pas de rs on prend un point
						$annot_genbo = join( ',', keys %{$hash->{$chr}->{$id}->{'gene'}->{$gene}->{'GenBo'}->{'annotation'}->{$trans_genbo}}); #on prend toutes les annotations concatene par des virgules
					}
					else {
						$rs_genbo    = '.';
						$annot_genbo = '.';
					}
					foreach my $trans_annovar ( keys (%{$hash->{$chr}->{$id}->{'gene'}->{$gene}->{'Annovar'}->{'annotation'}})) { #pour chaque transcrit annovar de chaque gene de chaque variant
						my ( $rs_annovar, $annot_annovar);
						if ( exists $hash->{$chr}->{$id}->{'gene'}->{$gene}->{'Annovar'} ) {
							$rs_annovar = $hash->{$chr}->{$id}->{'gene'}->{$gene}->{'Annovar'}->{'rs'};
							if ( not $rs_annovar ) { $rs_annovar = '.'; }
							$annot_annovar = $hash->{$chr}->{$id}->{'gene'}->{$gene}->{'Annovar'}->{'annotation'}->{$trans_annovar};
							$annot_annovar = join( ',', keys %{$hash->{$chr}->{$id}->{'gene'}->{$gene}->{'Annovar'}->{'annotation'}->{$trans_annovar}});
							$new_annot_annovar = $hAnnot_keys->{$annot_annovar}
							  if ( exists $hAnnot_keys->{$annot_annovar} ); #ici on remplace l'annotation annovar par une annotation genbo
						}
						else {
							$rs_annovar    = '.';
							$annot_annovar = '.';
						}

						$hAnnot->{ $annot_genbo . '-' . $annot_annovar }++;
						$total++;
						### a ce niveau on a obtenu toutes les annotations de genbo du transcrit, du gene, du variant dans $annot_genbo et idem pour annovar dans $new_annot_annovar
						my @lAnnG = split (",", $annot_genbo);
						my @lAnnA;
						if ($new_annot_annovar) { @lAnnA = split (",", $new_annot_annovar); }
						else { $new_annot_annovar = ''; }
						foreach my $annG (@lAnnG) { #pour chaque annotation genbo
							foreach my $annA (@lAnnA) { #pour chaque annotation annovar
								if ( $annG eq $annA ) { #si elles sont equivalentes
									my $lastCol;
									if (    $trans_genbo eq $trans_annovar and $trans_annovar ne '.' and $trans_genbo   ne '.' ) { #si meme transcrit mais pas un point
										$is_common_trans = 1;
										$is_common_annot = 1;
										$lastCol = 'same_trans'; #on certifie que le transcrit est bien identique
									}
									else {
										$is_common_annot = 1;
										$lastCol = 'same_annot'; #sinon on certifie que l'annotation est identique
									}
									### et on print la ligne dans le .common
									print COM "$id\t$ref\t$alt\t$gene\t$rs_genbo\t$annot_genbo\t$trans_genbo\t$rs_annovar\t$annot_annovar\t$trans_annovar\t$annot_genbo;$new_annot_annovar\t$lastCol\n";
								}
							}
						}
						if (not $is_common_annot) { #si les annotation genbo et annovar ne sont pas identiques on les sauvegarde dans une liste
							push (@lTmp, "$id\t$ref\t$alt\t$gene\t$rs_genbo\t$annot_genbo\t$trans_genbo\t$rs_annovar\t$annot_annovar\t$trans_annovar\t$annot_genbo;$new_annot_annovar\tnot_the_same\n");
						}
					}
				}
			}
			if (not $is_common_annot) { #pour le gene, si aucune annot n'est commune alors on print dans .others les divergeances
				foreach my $diff (@lTmp) {
					print DIF $diff;
				}
			}
		}
	}
	close(COM);
	close(DIF);
}

sub t_annotation_annovar_exonic {
	my ($self, $ccds_trans) = @_;	
	my $projectName = $self->name();
	my $path_tmp = $self->path_tmp();
	my $file = "$projectName.avinput.exonic_variant_function";
	open( FILE, "<$path_tmp/$file" )     or die "$!";
	open( TSV,  ">$path_tmp/$file.tsv" ) or die "$!";
	foreach my $ligne (<FILE>) {
		my @liste = split( '\t', $ligne );
		my $chr   = $liste[3];
		my $pos   = $liste[4];
		my $ref   = $liste[6];
		my $alt   = $liste[7];
		my $names = $liste[2];
		my @lTmp  = split( ":", $names );
		my $gene  = $lTmp[0];
		my $c;
		my @lTrans = ();
		foreach $c (@lTmp) {
			if ( $c =~ /(ENST[0-9]*)/ ) {
				push( @lTrans, $c );
			}
		}
		my $annot = $liste[1];
		foreach my $t (@lTrans) {
			if (exists $ccds_trans->{$t}) {
				print TSV "$chr\t$pos\t$ref\t$alt\t$gene\t$t\t$annot\n";
			}
		}
	}
}

sub t_annotation_annovar_multianno {
	my ($self, $ccds_trans) = @_;
	my $projectName = $self->name();
	my $path_tmp = $self->path_tmp();
	my $in = "$projectName.hg19_multianno.txt";
	open( TXT, "<$path_tmp/$in" )     or die "$!"; #on va chercher le fichier multianno d'annovar dans le dossier tmp
	open( TSV, ">$path_tmp/$in.tsv" ) or die "$!"; #on ouvre sa version tsv pour y entrer le formatage
	print TSV "chr\tstart\tref\tmut\tgene\ttranscript\tannot\trs_name\n";
	my $nb = 0;
	foreach my $ligne (<TXT>) {
		$nb++;
		next if ( $nb == 1 );#permet de passer l'en-tete
		chomp($ligne);
		my @lFields            = split( ' ', $ligne );
		my $chr                = $lFields[0];
		my $start              = $lFields[1];
		my $ref                = $lFields[3];
		my $mut                = $lFields[4];
		my $func_refgene       = $lFields[5];
		my $genes_name         = $lFields[6];
		my $geneDetail_refGene = $lFields[7];
		my $rs_name            = $lFields[-1];
		my @l_func_refgene     = split( ';', $func_refgene );
		my @lGene_byFct        = split( ';', $genes_name );
		my $nb_g               = 0;

		foreach my $this_funct (@l_func_refgene) {
			if ( $this_funct eq 'UTR5' or $this_funct eq 'UTR3' ) {
				$geneDetail_refGene = '.';
			}
			if ( $this_funct eq 'intergenic' ) {
				$geneDetail_refGene = '.';
			}
			my $this_genes = $lGene_byFct[$nb_g];
			my @lGenes = split( ',', $this_genes );

			my @lAnnot;
			my @lTrans;
			if ( $this_funct ne '.' ) { push( @lAnnot, $this_funct ); }
			if ( $geneDetail_refGene ne '.' ) {
				my @lTmp = split( ",", $geneDetail_refGene );
				foreach my $c (@lTmp) {
					if ( $c =~ /(ENST[0-9]*)/ ) {
						push( @lTrans, $1 );
					}
				}
			}
			my $annot = join( ',', @lAnnot );

			foreach my $gene (@lGenes) {
				if ( $gene ne 'NONE' ) { #pour chaque gene qui n'est pas NONE
					if ( $geneDetail_refGene ne '.' ) { #et qui possede une annotation /!\ position 7 et 9 des transcrits
						foreach my $transcript (@lTrans) { #pour chacun de ses transcrits
							if (exists $ccds_trans->{$transcript}) { #si le transcrit est dans la liste de tous les ccds humains
								print TSV "$chr\t$start\t$ref\t$mut\t$gene\t$transcript\t$annot\t$rs_name\n"; #on le print dans le fichier
							}
						}
					}
				}
			}
			$nb_g++;
		}
	}
	close(TXT);
	close(TSV);
}

sub t_annotation_annot_annovar {
	my ($self, $obj_patient) = @_;
	my $project = $self;
	my $projectName = $self->name();
	my $path_tmp = $self->path_tmp();
	my $path_annovar = $self->path_annovar();
	if (-e "$path_tmp/$projectName.avinput") {
		print "SKIP...\n";
		return;
	}
	### Annovar	
	# Genere un avinput a partir d'un vcf
	my $path_vcf = "$Bin/datas/tests_frequences/NGSTEST_FREQ/HG19/variations/haplotypecaller/patient.vcf.gz";
	
	my $cmd1 =
	"$path_annovar/convert2annovar.pl -format vcf4 -allsample -withfreq $path_vcf > $path_tmp/$projectName.avinput";
	warn $cmd1;
	system("$cmd1");	
	# Parse les DB d'annovar et donne un fichier tab
	my $cmd2 =
	"$path_annovar/table_annovar.pl $path_tmp/$projectName.avinput $path_annovar/humandb/ -buildver hg19 -out $path_tmp/$projectName -remove -protocol ensGene,genomicSuperDups,avsnp147 -operation g,r,f -nastring .";
	warn $cmd2;
	system("$cmd2");	
	# Genere un exonic_variant_function
	my $cmd3 =
	"$path_annovar/annotate_variation.pl -geneanno -buildver hg19 $path_tmp/$projectName.avinput $path_annovar/humandb/ -dbtype ensGene";
	warn $cmd3;
	system("$cmd3");
}

sub t_annotation_annot_genbo {
	my ($self, $obj_patient) = @_;
	my $project = $self;
	my $projectName = $self->name();
	my $path_tmp = $self->path_tmp();
	### Ouverture du tsv
	if ( -e "$path_tmp/$projectName\_tab-annotation.tsv" )
	{    #si le fichier existe deja on passe cette etape chronophage
		print "SKIP...\n";
		return;
	}
	open( TSV, ">$path_tmp/$projectName\_tab-annotation.tsv" ) or die "$!";
	## Creation de l'en-tete
	print TSV "chr\tstart\tref\talt\tgene\ttranscrit\tannotation\trs_name\n";
	## Choix chromosomes et appel de la fonction
	foreach my $this_chr_name ( 1 .. 22, 'X', 'Y', 'MT' ) {
		my ($listLignes) = $self->t_annotation_launch_genbo($this_chr_name, $obj_patient);
		foreach my $ligne (@{$listLignes}) {
			print TSV $ligne."\n";
		}
	}
	## Fermeture du tsv
	close(TSV);
}

sub t_annotation_launch_genbo {
	my ($self, $chr_name, $obj_patient) = @_;
	my $project = $self;
	my @lLignes;
	warn ' -> checking chr' . $chr_name;
	my $chr = $project->getChromosome($chr_name);  #on charge l'objet chromosome
	foreach my $var ( @{ $chr->getStructuralVariations() } ){
		if ($obj_patient) {
			next unless ( exists $var->annex->{ $obj_patient->id() } );
		}
		my $gene;
		my @lGene = @{ $var->getGenes() };
		foreach $gene (@lGene) { #pour chaque gene de la liste
			foreach my $t ( @{ $gene->getTranscripts() } ) { #pour chaque transcrit de ce gene
				my $chr_id = $chr->id(); #on prend son id chromosome
				my $start = $var->start(); #on prend la position start de la variation
				my $ref = $var->ref_allele(); #on prend l'allele de ref
				my $alt = $var->var_allele(); #on prend l'allele altere
				my ( $gene_id, $chrname1 ) = split( '_', $gene->id() ); #on prend le gene id EnsEMBL
				my ( $trans_id, $chrname2 ) = split( '_', $t->kyotoId() ); #on prend l'id du transcrit EnsEMBL
				my $annot = $var->variationType($t); #on prend l'annotation de la variation sur le transcrit
				my $rs = '';
				if ($var->rs_name()) {$rs = $var->rs_name();} #on recupere le rs
				my $ligne = "$chr_id\t$start\t$ref\t$alt\t$gene_id\t$trans_id\t$annot\t$rs";
				if ($t->ccds_name) { #on test si le transcrit est un ccds
					push(@lLignes, $ligne);
				}
			}
		}
	}
	return \@lLignes;
}

##### TESTS FREQUENCES ######


sub launch_tests_frequences {
	my $self = shift;
	$self->test_name('tests_frequences');
	my ($obj_patient) = $self->init_genbo_annovar();
	print "\n##### TEST FREQUENCES #####\n\n";
	
	print "\n### LAUNCHING annot_genbo\n";
	$self->t_frequences_annot_genbo($obj_patient);
	
	print "\n### LAUNCHING annot_annovar\n";
	$self->t_frequences_annot_annovar($obj_patient);
	$self->t_frequences_annovar2tsv();
	
	print "\n### LAUNCHING comparaison\n";
	$self->t_frequences_comparaison();
} 

sub t_frequences_comparaison {
	my $self = shift;
	my $projectName = $self->name();
	my $path_tmp = $self->path_tmp();
	open (GB, "<$path_tmp/$projectName\_tab-freq.tsv");
	open (AN, "<$path_tmp/$projectName\_freq.hg19_multianno.txt.tsv");
	my $h;
	my $hCat;
	foreach my $ligne (<GB>) {
		chomp($ligne);
		my @lCol = split ("\t", $ligne);
		my $var_id	= $lCol[0];
		my $rs_name = $lCol[1];
		my $freq 	= $lCol[2];
		$h->{$var_id}->{'GenBo'}->{'freq'} = $freq;
		$h->{$var_id}->{'GenBo'}->{'rs'} = $rs_name;
	}
	foreach my $ligne (<AN>) {
		chomp($ligne);
		my @lCol = split ("\t", $ligne);
		my $var_id	= $lCol[0];
		my $rs_name = $lCol[1];
		my $milleG 	= $lCol[2];
		my $ExAC 	= $lCol[3];
		$h->{$var_id}->{'Annovar'}->{'1000G'} = $milleG;
		$h->{$var_id}->{'Annovar'}->{'ExAC'} = $ExAC;
		$h->{$var_id}->{'Annovar'}->{'rs'} = $rs_name;
		
		#warn Dumper $h->{$var_id}; 
	}
	close(GB);
	close(AN);
	open (COM, ">$path_tmp/$projectName\_comparaison_freq_common.txt");
	open (OTH, ">$path_tmp/$projectName\_comparaison_freq_others.txt");
	foreach my $var_id (sort keys %{$h}) {
		next if (not exists $h->{$var_id}->{'GenBo'});
		next if (not exists $h->{$var_id}->{'Annovar'});
		my $rs_genbo = $h->{$var_id}->{'GenBo'}->{'rs'};
		my $rs_annovar = $h->{$var_id}->{'Annovar'}->{'rs'};
		my $rs;
		if ($rs_genbo eq $rs_annovar) {$rs = $rs_genbo;}
		else {$rs = "$rs_genbo\/$rs_annovar";}
		my $this_freq_Genbo = $h->{$var_id}->{'GenBo'}->{'freq'};
		my $this_freq_1KG = $h->{$var_id}->{'Annovar'}->{'1000G'};
		my $this_freq_ExAC = $h->{$var_id}->{'Annovar'}->{'ExAC'};
		my ($cat_Genbo, $freq_Genbo) = $self->t_frequences_categorie($this_freq_Genbo);
		my ($cat_1KG, $freq_1KG) = $self->t_frequences_categorie($this_freq_1KG);
		my ($cat_ExAC, $freq_ExAC) = $self->t_frequences_categorie($this_freq_ExAC);
		my $triplet = "$cat_Genbo,$cat_1KG,$cat_ExAC";
		if ($cat_Genbo eq $cat_1KG) {
			print COM "$var_id,$rs;$freq_Genbo,$freq_1KG,$freq_ExAC;$cat_Genbo,$cat_1KG,$cat_ExAC\n";
			$hCat->{'common'}->{$triplet}++;
		}
		elsif ($cat_Genbo eq $cat_ExAC) {
			print COM "$var_id,$rs;$freq_Genbo,$freq_1KG,$freq_ExAC;$cat_Genbo,$cat_1KG,$cat_ExAC\n";
			$hCat->{'common'}->{$triplet}++;
		}
		else {
			print OTH "$var_id,$rs;$freq_Genbo,$freq_1KG,$freq_ExAC;$cat_Genbo,$cat_1KG,$cat_ExAC\n";
			$hCat->{'others'}->{$triplet}++;
		}
		$hCat->{$triplet}++;
	}
	close (COM);
	close (OTH);
	my $ok  = 0;
	my $not = 0;
	my (@lCommon, @lOthers) = ();
	foreach my $triplets (sort keys %{$hCat}) {
		if (exists $hCat->{'common'}->{$triplets}){
			my $occ = $hCat->{'common'}->{$triplets};
			push (@lCommon, "$triplets:$occ\n");
			$ok += $occ;
		}
		if (exists $hCat->{'others'}->{$triplets}){
			my $occ = $hCat->{'others'}->{$triplets};
			push (@lOthers, "$triplets:$occ\n");
			$not += $occ;
		}
	}
	print "\n";
	print "[common]\n";
	foreach my $ligne (@lCommon) {
		print $ligne;
	}
	print "\n";
	print "[others]\n";
	foreach my $ligne (@lOthers) {
		print $ligne;
	}
	print "\n";
	print "[resume]\n";
	my $p_ok = ($ok / ($ok + $not)) * 100;
	my $p_not = 100 - $p_ok;
	$p_ok  = sprintf("%.2f", $p_ok);
	$p_not = sprintf("%.2f", $p_not);
	print 'Freq OK: '.$ok." -> $p_ok%\n";
	print 'Freq PB: '.$not." -> $p_not%\n\n";
}

sub t_frequences_categorie {
	my ($self, $freq) = @_;
	if (not defined($freq)) {
		return ('no_value', $freq);
		die;
	}
	if ( $freq eq "." or $freq == -1) {
		return ("no_freq", $freq);
	}
	$freq = sprintf ("%f", $freq);
	if ( $freq >= 0.05 ) {
		return ("more_than_5", $freq);
	}
	elsif ( $freq < 0.05 and $freq >= 0.01 ) {
		return ("less_than_5", $freq);
	}
	elsif ( $freq < 0.01 and $freq >= 0.001 ) {
		return ("less_than_1", $freq);
	}
	elsif ( $freq < 0.001 and $freq >= 0.0001 ) {
		return ("less_than_01", $freq);
	}
	elsif ( $freq < 0.0001 ) {
		return ("less_than_0001", $freq);
	}
}

sub t_frequences_annovar2tsv {
	my $self = shift;
	my $projectName = $self->name();
	my $path_tmp = $self->path_tmp();
	my $in = "$projectName\_freq.hg19_multianno.txt";
	open( TXT, "<$path_tmp/$in" )     or die "$!"; #on va chercher le fichier multianno d'annovar dans le dossier tmp
	open( TSV, ">$path_tmp/$in.tsv" ) or die "$!"; #on ouvre sa version tsv pour y entrer le formatage
	my $nb = 0;
	foreach my $ligne (<TXT>) {
		$nb++;
		next if ( $nb == 1 );#permet de passer l'en-tete
		chomp($ligne);
		my @lFields	 = split( ' ', $ligne );
		my $chr		 = $lFields[0];
		$chr 		 =~ s/chr//;
		my $start	 = $lFields[1];
		my $ref		 = $lFields[3];
		my $alt		 = $lFields[4];
		my $var_id	 = "$chr\_$start\_$ref\_$alt";
		my $milleG	 = $lFields[5];
		my $ExAC_ALL = $lFields[6];
		my $ExAC_AFR = $lFields[7];
		my $ExAC_AMR = $lFields[8];
		my $ExAC_EAS = $lFields[9];
		my $ExAC_FIN = $lFields[10];
		my $ExAC_NFE = $lFields[11];
		my $ExAC_OTH = $lFields[12];
		my $ExAC_SAS = $lFields[13];
		my $rs_name	 = $lFields[-1];
		print TSV "$var_id\t$rs_name\t$milleG\t$ExAC_ALL\n"; #on le print dans le fichier
	}
	close(TXT);
	close(TSV);
}

sub t_frequences_annot_annovar {
	my ($self, $obj_patient) = @_;
	my $project = $self;
	my $projectName = $self->name();
	my $path_tmp = $self->path_tmp();
	my $path_annovar = $self->path_annovar();
	if (-e "$path_tmp/$projectName\_freq.avinput") {
		print "SKIP...\n";
		return;
	}
	### Annovar	
	# Genere un avinput a partir d'un vcf
	my $path_vcf = $obj_patient->getVariationsFiles->[0];
	my $cmd1 =
	"$path_annovar/convert2annovar.pl -format vcf4 -allsample -withfreq $path_vcf > $path_tmp/$projectName\_freq.avinput";
	warn $cmd1;
	system("$cmd1");	
	# Parse les DB d'annovar et donne un fichier tab
	my $cmd2 =
	"$path_annovar/table_annovar.pl $path_tmp/$projectName\_freq.avinput $path_annovar/humandb/ -buildver hg19 -out $path_tmp/$projectName\_freq -remove -protocol 1000g2015aug_all,exac03,avsnp147 -operation f,f,f -nastring .";
	warn $cmd2;
	system("$cmd2");
} 

sub t_frequences_annot_genbo {
	my ($self, $obj_patient) = @_;
	my $project = $self;
	my $projectName = $self->name();
	my $path_tmp = $self->path_tmp();
	### Ouverture du tsv
	if ( -e "$path_tmp/$projectName\_tab-freq.tsv" )
	{    #si le fichier existe deja on passe cette etape chronophage
		print "SKIP...\n";
		return;
	}
	open( TSV, ">$path_tmp/$projectName\_tab-freq.tsv" ) or die "$!";
	## Creation de l'en-tete
	print TSV "var_id\trs_name\tfreq\n";
	## Choix chromosomes et appel de la fonction
	foreach my $this_chr_name ( 1 .. 22, 'X', 'Y', 'MT' ) {
		my ($listLignes) = $self->t_frequences_launch_genbo($this_chr_name, $obj_patient);
		foreach my $ligne (@{$listLignes}) {
			print TSV $ligne."\n";
		}
	}
	## Fermeture du tsv
	close(TSV);
}

sub t_frequences_launch_genbo {
	my ($self, $chr_name, $obj_patient) = @_;
	my $project = $self;
	my @lLignes;
	warn ' -> checking chr' . $chr_name;
	my $chr = $project->getChromosome($chr_name);  #on charge l'objet chromosome
	foreach my $var ( @{ $chr->getStructuralVariations() } ){
		if ($obj_patient) {
			next unless ( exists $var->annex->{ $obj_patient->id() } );
		}
		my $chr_id 	= $chr->id(); #on prend son id chromosome
		my $start 	= $var->start(); #on prend la position start de la variation
		my $ref 	= $var->ref_allele(); #on prend l'allele de ref
		my $alt 	= $var->var_allele(); #on prend l'allele altere
		my $var_id 	= "$chr_id\_$start\_$ref\_$alt";
		my $rs    	= '';
		$rs    		= $var->rs_name() if ($var->rs_name()); #on recupere le rs
		my $freq  	= $var->frequency();
		my $ligne 	= "$var_id\t$rs\t$freq";
		push(@lLignes, $ligne);
	}
	return \@lLignes;
}


1;