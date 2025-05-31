package parse_wisecondor;
use strict;
use Data::Dumper;
use Storable qw(dclone);
sub parse_cnv {
	 my ($patient) = @_;
	 my $caller = "wisecondor";
	  my $fichierPatient = $patient->getSVFile($caller);
	 		my $res;
			my $project = $patient->project;
			my $name = $patient->name;
			# ouverture du fichier wisecondor et parsing
			my $fd;
			open($fd," zcat $fichierPatient | ") or die("open: $!");
			
			# lecture de la première ligne
			my $ligne = <$fd>;
			my $Caller = "wisecondor";
			# traitement des lignes suivantes pour le fichier en entier
			while( defined( $ligne = <$fd> ) )
			{
						my @champs = split(/\t/,$ligne);
						my $chr = $project->getChromosome($champs[0],1);
						my $WC_ratio=$champs[3];
						my $QUAL=$champs[4];
				  		 
						
						my $num = $chr->name();
						
						my $h;
						
						
						$h->{'SVTYPE'}=$champs[5];
						$h->{'SVTYPE'} = "DEL"  if  ($h->{'SVTYPE'} =~ m/loss/);
						$h->{'SVTYPE'} = "DUP" if ($h->{'SVTYPE'} =~ m/gain/);
						$h->{'CHROM'}=$num;
						$h->{'START'}=$champs[1];
						$h->{'END'}=$champs[2];
						$h->{'SVLEN'}= abs($champs[1]-$champs[2])+1;
						$h->{'GT'}="-";
						$h->{'CN'}="-";
						$h->{'RATIO'}=$champs[3];
						$h->{'QUAL'}=$champs[4];
						$h->{KARYOTYPE_ID}= $chr->karyotypeId;
						$h->{'INFOS'}->{zscore} = $champs[4];
						my $id = $h->{'SVTYPE'}."_".$chr->name."_".$h->{'START'}."_".$h->{'END'};
						$h->{'ELEMENTARY'}= [$id];
						$h->{'REAL_CALLER'} = $caller;
						$h->{'CALLER'} = $caller;
						$h->{id}=$id;
						my $type = $h->{'SVTYPE'};
						$res->{$id} = $h;
			}
			return $res;
}
1;
