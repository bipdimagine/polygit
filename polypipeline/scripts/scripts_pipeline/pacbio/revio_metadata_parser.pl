#!/usr/bin/perl

# revio_metadata_parser.pl
#
# Description : Extrait le nom du run, les noms de patients (BioSamples),
#               les barcodes et construit les identifiants Lane depuis les
#               fichiers metadata XML des SMRT Cells PacBio Revio.
#
# Usage   : revio_metadata_parser.pl [--run <run_name> [--run <run_name2> ...]]
#           revio_metadata_parser.pl [--runs <run_name>[,<run_name2> ...]]
#           Sans option : menu interactif de sélection des runs disponibles
#

use strict;
use warnings;
use XML::LibXML;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../packages/";
use GBuffer;
use Data::Dumper;
use Term::Menus;

my $base_dir = '/data-pure/raw-data/PACBIO/REVIO/IMAGINE/';

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
my @run_dir_names;

GetOptions(
	'runs=s{1,}'	=> \@run_dir_names,
    'help|h'		=> sub {pod2usage(-verbose => 2, -noperldoc=>1)},
) or pod2usage(-message => "Error in command line arguments\n", -verbose => 2, -exitval => 2, -noperldoc=>1);


@run_dir_names = split(/,/, join(',',@run_dir_names));

# ---------------------------------------------------------------------------
# Si aucun run fourni : menu interactif de sélection
# ---------------------------------------------------------------------------
unless (@run_dir_names) {
    my @available = sort {$b cmp $a} grep { -d "$base_dir/$_" } do {
        opendir(my $dh, $base_dir) or die "Impossible d'ouvrir $base_dir : $!\n";
        grep { /^r\d+_/ } readdir($dh);
    };
#    warn Dumper \@available;

    die "Aucun run trouvé dans $base_dir\n" unless scalar @available;

#    my @selected = pick(
#        \@available,
#        " Sélectionnez un ou plusieurs runs (Espace = sélection, Entrée = valider) :"
#    );
    my %Menu_1 = (
		Item_1 => {
			Text   => "]Convey[",
			Convey => \@available,
		},
		Select => 'Many',
		Banner => " Sélectionnez un ou plusieurs runs :",
		Display => 20,
	);
	my @selected = &Menu( \%Menu_1 );
    

    unless (@selected && defined $selected[0]) {
        print "Aucun run sélectionné. Abandon.\n";
        exit 0;
    }
    die("bye") if $selected[0] eq ']quit[';

    @run_dir_names = @selected;
}

# ---------------------------------------------------------------------------
# Structures de données
#
# %patients = (
#   patient_name => {
#     bc   => 'bc2076',                          # barcode court (sans '--...')
#     runs => { run_name => [ cell_id, ... ] }   # ex: ['1_C01', '1_D01']
#   }
# )
# %run_info = (
#   run_name => { name => '...', instrument => '...', cells => [...] }
# )
# ---------------------------------------------------------------------------
my %patients;
my @run_order;   # conserve l'ordre d'insertion pour l'affichage
my %run_info;

my $parser = XML::LibXML->new();

# ---------------------------------------------------------------------------
# Parsing de chaque run
# ---------------------------------------------------------------------------
foreach my $run_dir_name (@run_dir_names) {
    $run_dir_name =~ s/\/$//;
    my $run_dir = $base_dir . $run_dir_name;

    die "Répertoire introuvable : $run_dir\n" unless -d $run_dir;

    push @run_order, $run_dir_name;
    $run_info{$run_dir_name} = { name => 'N/A', instrument => 'N/A', cells => [] };

    # Exclure les fichiers *preview* (*_basic_preview_*, *_full_preview_*, etc.)
    my @metadata_files =
      sort grep { !/preview/ } glob("$run_dir/*/metadata/*.metadata.xml");

    warn "Aucun fichier metadata XML trouvé dans : $run_dir\n"
      unless @metadata_files;

    foreach my $xml_file (@metadata_files) {

        # Extraction du cell_id depuis le chemin (ex: '1_C01')
        my $cell_id = 'unknown';
        $cell_id = $1 if $xml_file =~ m{/(\d+_[A-Z]\d+)/metadata/};

        # Enregistrement du cell_id (sans doublon)
        unless (grep { $_ eq $cell_id } @{ $run_info{$run_dir_name}{cells} }) {
            push @{ $run_info{$run_dir_name}{cells} }, $cell_id;
        }

        my $doc;
        eval { $doc = $parser->parse_file($xml_file) };
        if ($@) {
            warn "Erreur parsing $xml_file : $@\n";
            next;
        }

        my $xpc = XML::LibXML::XPathContext->new($doc);
        $xpc->registerNs('pbdm', 'http://pacificbiosciences.com/PacBioDataModel.xsd');
        $xpc->registerNs('pbcm', 'http://pacificbiosciences.com/PacBioCollectionMetadata.xsd');
        $xpc->registerNs('pbsi', 'http://pacificbiosciences.com/PacBioSampleInfo.xsd');

        # Nom du run ? collecté une seule fois par run
        if ($run_info{$run_dir_name}{name} eq 'N/A') {
            my ($run_node) = $xpc->findnodes('//pbdm:Run');
            if ($run_node) {
                my $rn = $run_node->getAttribute('Name') // 'N/A';
                $rn =~ s/^\s+|\s+$//g;
                $run_info{$run_dir_name}{name} = $rn;
            }
        }

        my @collections = $xpc->findnodes('//pbcm:CollectionMetadata');

        foreach my $coll (@collections) {

            # Instrument - collecté une seule fois par run
            if ($run_info{$run_dir_name}{instrument} eq 'N/A') {
                $run_info{$run_dir_name}{instrument} =
                  $coll->getAttribute('InstrumentName') // 'N/A';
            }

            my @biosamples = $xpc->findnodes(
                'pbcm:WellSample/pbsi:BioSamples/pbsi:BioSample', $coll
            );

            foreach my $bs (@biosamples) {
                my $bs_name = $bs->getAttribute('Name') // 'N/A';
                $bs_name =~ s/^\s+|\s+$//g;

                my ($bc_node) =
                  $xpc->findnodes('pbsi:DNABarcodes/pbsi:DNABarcode', $bs);
                next unless $bc_node;

                # Barcode court : 'bc2076--bc2076' ? 'bc2076'
                my $bc_full  = $bc_node->getAttribute('Name') // 'N/A';
                my $bc_short = $bc_full;
                $bc_short =~ s/--.*$//;

                $patients{$bs_name}{bc} //= $bc_short;

                # Ajout du cell_id sans doublon pour ce patient/run
                unless (grep { $_ eq $cell_id }
                    @{ $patients{$bs_name}{runs}{$run_dir_name} // [] })
                {
                    push @{ $patients{$bs_name}{runs}{$run_dir_name} }, $cell_id;
                }
            }
        }
    }
}

# ---------------------------------------------------------------------------
# Affichage des infos de run ? une seule fois par run
# ---------------------------------------------------------------------------
foreach my $run_dir_name (@run_order) {
    my $info = $run_info{$run_dir_name};
    printf "Run dir    : %s\n", $run_dir_name;
    printf "Run name   : %s\n", $info->{name};
#    printf "Instrument : %s\n", $info->{instrument};
    printf "SMRT Cells : %s\n", join(', ', @{ $info->{cells} });
    print "\n";
}

# ---------------------------------------------------------------------------
# Tableau TSV : Patient / BC / Lane
#
# Lane par run = run:wells:bc
#   run   :  ex: r84301,20260605,125505
#   wells : cell_ids, joinés par ','
#               ex: 1_C01 + 1_D01 -> 1_C01,1_D01
#   bc        : barcode court ex: bc2076
#
# Si plusieurs runs : Lane = lane_run1;lane_run2;...
# ---------------------------------------------------------------------------
print join("\t", "Patient", "BC", "Lane") . "\n";

foreach my $patient (sort keys %patients) {
    my $bc = $patients{$patient}{bc} // 'N/A';

    my @lane_parts;
    foreach my $run_dir_name (@run_order) {
        next unless exists $patients{$patient}{runs}{$run_dir_name};

        # run_fmt : underscores ? commas
        (my $run_fmt = $run_dir_name);

        # wells_fmt : chaque cell_id underscore ? comma, puis joinés par ','
        my @wells     = @{ $patients{$patient}{runs}{$run_dir_name} };
        my $wells_fmt = join(',', @wells);

        push @lane_parts, "${run_fmt}:${wells_fmt}:${bc}";
    }

    print join("\t", $patient, $bc, join(';', @lane_parts)) . "\n";
}


__END__

=pod

=encoding UTF-8

=head1 NOM

revio_metadata_parser.pl - Extraction des métadonnées de runs PacBio Revio

=head1 SYNOPSIS

    revio_metadata_parser.pl [options]

    Options :
      --runs <run_name>   Nom du répertoire de run (répétable)
      --help             Affiche la documentation

    Sans option --run : menu interactif de sélection des runs disponibles.

=head1 OPTIONS

=over 4

=item B<--run> I<run_name>

Nom du répertoire de run à analyser, tel qu'il apparaît dans :

    /data-pure/raw-data/PACBIO/REVIO/IMAGINE/

Exemple : C<r84301_20260605_125505>

L'option peut être répétée pour combiner plusieurs runs :

    --run r84301_20260605_125505 --run r84302_20260610_130000

ou les runs peuvent être séparés par C<,>:

    --run r84301_20260605_125505,r84302_20260610_130000

=item B<--help>, B<-h>

Affiche l'aide et quitte.

=back

=head1 DESCRIPTION

Ce script parcourt les fichiers metadata XML des SMRT Cells d'un ou plusieurs
runs PacBio Revio et produit :

=over 4

=item * Un résumé par run (nom, instrument, SMRT Cells détectées)

=item * Un tableau TSV (Patient / BC / Lane) sur la sortie standard

=back

Les fichiers C<*_preview_*.metadata.xml> (basic_preview, full_preview) sont
automatiquement exclus.

=head2 Format de la colonne Lane

    <run>:<wells>:<bc>

=over 4

=item * B<run>   : nom du répertoire de run 
                       ex : C<r84301_20260605_125505>

=item * B<wells> : liste des SMRT Cell IDs, séparées par C<,>
                       ex : C<1_C01,1_D01>

=item * B<bc>        : barcode court (partie gauche du barcode symétrique)
                       ex : C<bc2076>

=back

Si plusieurs runs sont fournis, les lanes sont concaténées avec C<;> :

    r84301_20260605_125505:1_C01:bc2076;r84302_20260610_130000:1_B01:bc2076

=head1 EXEMPLES

    # Menu interactif
    ./revio_metadata_parser.pl

    # Un seul run en ligne de commande
    ./revio_metadata_parser.pl --run r84301_20260605_125505

    # Deux runs combinés
    ./revio_metadata_parser.pl --run r84301_20260605_125505 --run r84302_20260610_130000
    ./revio_metadata_parser.pl --runs r84301_20260605_125505,r84302_20260610_130000

    # Redirection de la sortie TSV
    ./revio_metadata_parser.pl --run r84301_20260605_125505 > samples.tsv

=head1 AUTEUR

Plateforme IMAGINE - Bioinformatique

=cut
