#!/usr/bin/env perl
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use  File::Temp;
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use CGI qw/:standard :html3/;
use Storable qw(store retrieve freeze);
use List::MoreUtils qw(natatime);
use JSON::XS;
use GBuffer;  
use MCE::Loop;

 my $project_name;
my $final_vcf;
my $log_file;
my $patient_name;
my $fork;
my $cgi = new CGI();

my $project_name = $cgi->param("project");
my $force = $cgi->param("force");

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );
my $dir = $project->getVariationsDir("amplicon");
	my $json_file = "$dir/pileup.json";
	unlink $json_file if $force;
if (-e $json_file){


open(my $fh, "<", $json_file)

    or die "Impossible d'ouvrir $json_file : $!";

local $/;

my $json = <$fh>;

close($fh);

print $cgi->header(

    -type => 'application/json',

    -charset => 'utf-8'

);
 print $json;
 exit(0);
 
}

my $file =  $project->project_path."/hotspot.bed";
my $positions = read_bed($file);
$project->preload();
$project->disconnect();
my $pileup = {};

MCE::Loop::init {

    chunk_size => 1,

    max_workers => 8,

    gather => sub {

        my ($result) = @_;

        my $patient = $result->{patient};

        foreach my $region (keys %{$result->{data}}) {

            $pileup->{$region}{$patient} = $result->{data}{$region};

        }

    }

};

mce_loop {

    my ($mce,$patients) = @_;
	my $patient = $patients->[0];
    my $patient_name = $patient->name;

    my $bam = $patient->getAlignmentFile();

    warn "START $patient_name : $bam\n";

    my $hts = Bio::DB::HTS->new(

        -bam => $bam,

    );
	$hts->max_pileup_cnt(100000);
    my %patient_pileup;

    foreach my $p (@$positions) {

        my $region =  $p->{chr} . ":" . $p->{start} . "-" . $p->{end}.":".$p->{ref}."/".$p->{alt};

        my $r = _pileup(

            $p->{chr},

            $p->{start},

            $p->{end},

            $hts

        );

        $patient_pileup{$region} = {

            pileup => $r

        };

    }

    MCE->gather({

        patient => $patient_name,

        data    => \%patient_pileup

    });

} @{$project->getPatients};

MCE::Loop::finish;

#
#my $pileup;
#
#foreach my $patient (@{$project->getPatients()}){
#	my $bam = $patient->getAlignmentFile();
#	    my $hts = Bio::DB::HTS->new(
#
#        -bam => $bam,
#
#    );
#    
#	foreach my $p (@$positions) {
#	 my $r = _pileup ($p->{chr},$p->{start},$p->{end},$hts); 
#	 $pileup->{$p->{chr}.":".$p->{start}."-".$p->{end}}->{$patient->name}->{pileup} = $r;
#	 
#		
#	}
#}

print $cgi->header(

    -type => 'application/json',

    -charset => 'utf-8'

);

open(my $fh, ">", $json_file)

    or die "Impossible d'ouvrir $json_file : $!";

print $fh encode_json($pileup);

close($fh) or die "Erreur fermeture $json_file : $!";
print encode_json($pileup);

exit(0);


sub _pileup {

    my ($chr, $start, $end, $hts) = @_;

    my %res;

    # BED = 0-based / end exclusive

    #

    # pileup() utilise une région 1-based.

    #

    # BED chr 100 101

    # correspond à la position 101 en coordonnées 1-based.

    my $region = $chr . ":" . ($start ) . "-" . ($end);
	warn $region;
    # --------------------------------------------------------

    # Callback appelé pour chaque position

    # --------------------------------------------------------

    my $callback = sub {

        my ($seqid, $pos, $pileups) = @_;

        # Sécurité : pileup peut légèrement déborder

        # de la région demandée.
		
        return if $pos < $start;

        return if $pos > $end;
        my $r = {

            A     => 0,

            T     => 0,

            C     => 0,

            G     => 0,

            N     => 0,

            INS   => 0,

            DEL   => 0,

            depth => 0,

        };

        foreach my $p (@$pileups) {
            $r->{depth}++;

            # ------------------------------------------------

            # Délétion

            # ------------------------------------------------

            if ($p->is_del) {

                $r->{DEL}++;

                next;

            }

            # ------------------------------------------------

            # Base du read

            # ------------------------------------------------

            my $alignment = $p->alignment;

            my $qpos = $p->qpos;

            my $seq = $alignment->qseq;

            my $base = substr($seq, $qpos, 1);

            $base = uc($base);

            if ($base eq 'A') {

                $r->{A}++;

            }

            elsif ($base eq 'T') {

                $r->{T}++;

            }

            elsif ($base eq 'C') {

                $r->{C}++;

            }

            elsif ($base eq 'G') {

                $r->{G}++;

            }

            else {

                $r->{N}++;

            }

            # ------------------------------------------------

            # Insertion

            #

            # indel > 0 = insertion après cette position

            # ------------------------------------------------

            my $indel = $p->indel;

            if ($indel > 0) {

                $r->{INS}++;

            }

        }

        # On ne stocke que les positions réellement couvertes.

        #

        # Si tu veux aussi les positions avec profondeur 0,

        # il faudra les initialiser avant le pileup.

        $res{$pos} = $r;

    };

    # --------------------------------------------------------

    # Exécution du pileup

    # --------------------------------------------------------

    $hts->pileup(

        $region,

        $callback

    );

    return \%res;

}

sub read_bed {
	my ($bed_file) = @_;
	my @positions;

open(my $BED, '<', $bed_file)

    or die "Cannot open $bed_file: $!\n";

while (<$BED>) {

    chomp;

    next if /^\s*$/;

    next if /^\s*#/;

    my @line = split(" ");

    die "BED invalide : $_ \n".scalar(@line) if @line < 3;

    my ($chr, $start, $end,$all) = @line[0,1,2,3];
	my ($ref,$alt) = split("/",$all);

    push @positions, {
        chr   => $chr,
        start => $start,
        end   => $end,
        ref=> $ref,
        alt=> $alt
    };

}

close($BED);
return \@positions;
	
}