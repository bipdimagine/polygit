package SVParser;
use strict;
use Data::Dumper;
use Storable qw(dclone);
sub parse_wisecondor {
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
						my $id = $h->{'SVTYPE'}."_".$chr->name."_".$h->{'START'}."_".$h->{'END'};
						$h->{'ELEMENTARY'}= [$id];
						$h->{'INFOS'}->{zscore} = $champs[4];
						$h->{'REAL_CALLER'} = $caller;
						$h->{id}=$id;
						my $type = $h->{'SVTYPE'};
						control_object($h);
						$res->{$id} = $h;
			}
			return $res;
}



sub parse_vcf {
	my ($patient,$caller) =@_;
	 my $vcf_file = $patient->getSVFile($caller);
	 my $project = $patient->project;
	# ouverture du fichier manta zippé
	my $res;
	my $vcf = Bio::DB::HTS::VCF->new( filename => "$vcf_file" );
	my $header = $vcf->header();
	
	while (my $row = $vcf->next){
		
		my $chr = $project->getChromosome($row->chromosome($header),1);
		
		next unless $chr;
	
		my $h;
		#champ infos 
		next unless $row->has_filter($header,".");
		my $svtype = 	$row->get_info($header, "SVTYPE");
	
		$h->{'SVTYPE'} = get_value($row->get_info($header, "SVTYPE"));
		if ($h->{'SVTYPE'}  eq "CNV"){
			$h->{'SVTYPE'} = "DEL"  if  ($row->id()  =~ m/loss/i);
			$h->{'SVTYPE'} = "DUP" if ($row->id()  =~ m/gain/i);
		}
		next unless ( ($h->{'SVTYPE'} eq "DUP") || ($h->{'SVTYPE'} eq "DEL") );
	
		$h->{'CHROM'}=$chr->name;
		$h->{"CALLER"} = $caller;
		$h->{'END'} = $row->get_info($header, "END")->[0];
		$h->{'START'} = $row->position();
		$h->{'SVLEN'} =  abs($h->{'END'} - $h->{'START'});
		$h->{'KARYOTYPE_ID'}= $chr->karyotypeId;
		$h->{'QUAL'} = $row->quality() ;
		if ($caller eq "pbsv"){
			$h->{'QUAL'} =  1;
		}
		$h->{'GT_a'} = $row->get_format($header, "GT");
		
		$h->{'GT'} = "0/1";
		$h->{'GT'} = "1/1" if ($h->{'GT_a'}->[0] == $h->{'GT_a'}->[1]);  
		$h->{'CN'} ="-" ;
		$h->{'RATIO'} ="-" ;
		$h->{'CN'} = get_value($row->get_format($header, "CN")); 
		$h->{'CN'} = 0   unless $h->{'CN'};
		my $id = $h->{'SVTYPE'}."_".$chr->name."_".$h->{'START'}."_".$h->{'END'};
		$h->{'ELEMENTARY'}= [$id];
		$h->{'INFOS'} = $row->get_format($header);
		$h->{id}=$id;
		$h->{'REAL_CALLER'} = $caller;
		$h->{'CALLER'} = $caller;
		my $type = $h->{'SVTYPE'};
		my $num = $chr->name();
		$res->{$id} = $h;
		control_object($h);
		
	}
	
return $res;
	
}

sub get_value {
	my($res) = @_;

	if (ref($res) eq 'ARRAY') {
			return  $res->[0];	
	}
	else {
		return undef;
	}
}

sub gatherCNV_from_samecaller{
	my ($name,$hcnv) = @_;
	my $htree;
	my $hintspan;

	# 1)  detecter les SV chevauchants 
	
	# creer les arbres
	

foreach my $id  (keys %{$hcnv}) {
						my ($t,$c,$d,$f) = split( /_/, $id );
						
						$htree->{$t}->{$c}= Set::IntervalTree->new unless exists $htree->{$t}->{$c} ;
						$hintspan->{$t}->{$c}= Set::IntSpan::Fast->new() unless exists $hintspan->{$t}->{$c} ;;
						$htree->{$t}->{$c}->insert($id,$d,$f);
	}
	

	# 2) associer a chaque SV trouve avec le même caller ceux qui lui sont proches	
foreach my $id  (keys %{$hcnv}) 
				{
						#on cherche les regroupements
						my ($t,$c,$dtheSV,$ftheSV) = split( /_/, $id );
		
						my $padding = 0.1 * abs($ftheSV-$dtheSV);
		
						my $tab_id = $htree->{$t}->{$c}->fetch($dtheSV-$padding,$ftheSV+$padding);

						# on regroupe les id dans un intspan
						my $gdeb=0;
						my $gend=0;
			
						foreach my $ind_id ( @$tab_id ) 
						{
							my ($t,$c,$d,$f) = split(/_/, $ind_id);
				
		 				   	$gdeb = $d if ( ($d < $gdeb) || ($gdeb==0));
		 					$gend = $f if ( ($f > $gend) );
 						}
						$hintspan->{$t}->{$c}->add_range($gdeb,$gend);
				}

	return gather_id($htree,$hintspan,$hcnv);
}



sub gather_id {
	
		my ($htree,$hintspan,$hcnv) = @_;
		my $nb = 0;
		my $hcnv_gather;
		foreach my $type (keys %{$hintspan})
		{
			foreach my $num (keys %{$hintspan->{$type}})
			{
				my $liste_of_bornes = $hintspan->{$type}->{$num}->as_string() if defined( $hintspan->{$type}->{$num});
				$liste_of_bornes .= "," if ($liste_of_bornes !~ m/,/);
				my @tab_bornes = split(/,/,$liste_of_bornes);
				
				foreach my $bornes (@tab_bornes)
				{
					my ($start,$end) = split(/-/,$bornes);
					
					my $tab_ind_id = $htree->{$type}->{$num}->fetch($start,$end);
					
					my $global_id = $type."_".$num."_".$start."_". $end;
		
					# enregistrer les infos dans la table de hash 
					$hcnv_gather->{$global_id}->{'id'} = $global_id;			
					$hcnv_gather->{$global_id}->{'SVTYPE'} = $type;
					$hcnv_gather->{$global_id}->{'CHROM'} = $num;
					$hcnv_gather->{$global_id}->{'START'} = $start;
					$hcnv_gather->{$global_id}->{'END'} = $end;
					$hcnv_gather->{$global_id}->{'SVLEN'}  = abs( $start - $end );
		
					# retrouver les informations correspondant aux differents id regroupes 
					# ici il faut recuperer les valeurs max ou les concatenation de liste
					
					my $ggfreq = -1;
					my $glfreq = -1;
					my $mg ="no";
					my $dbVar_status = "";
					my $rankannot = 0;
					my $gt = "";
					my $elementary_ids = [];
					my $cn;
					my $ratio;
					my $qual;
					my $index;
					my $cn = -99;
					
					foreach my $id ( @$tab_ind_id )
					{		
						$ggfreq = $hcnv->{$id}->{'GOLD_G_FREQ'} if ($hcnv->{$id}->{'GOLD_G_FREQ'} > $ggfreq);
						$glfreq = $hcnv->{$id}->{'GOLD_L_FREQ'} if ($hcnv->{$id}->{'GOLD_L_FREQ'} > $glfreq);
						$mg = "yes" if ($hcnv->{$id}->{'OMIN_MG'} eq "yes");
						$dbVar_status .= $hcnv->{$id}->{'dbVar_status'}." ";
						$rankannot = $hcnv->{$id}->{'RANKAnnot'} if ($hcnv->{$id}->{'RANKAnnot'} > $rankannot);
						my $text = $hcnv->{$id}->{'GT'};
						$gt .= $hcnv->{$id}->{'GT'};#." " unless ($gt =~ m/$text/) ;
						$cn = $hcnv->{$id}->{'CN'} if $cn < $hcnv->{$id}->{'CN'};
						$ratio += $hcnv->{$id}->{'RATIO'};
						$qual += $hcnv->{$id}->{'QUAL'};
						$hcnv_gather->{$global_id}->{'INFOS'}  = $hcnv->{$id}->{'INFOS'};
						$hcnv_gather->{$global_id}->{REAL_CALLER} = $hcnv->{$id}->{'REAL_CALLER'};
						$hcnv_gather->{$global_id}->{CALLER} = $hcnv->{$id}->{'CALLER'};
						push(@$elementary_ids,$id);
					}
					$hcnv_gather->{$global_id}->{'ELEMENTARY'} = $elementary_ids;
					
					$hcnv_gather->{$global_id}->{'GOLD_G_FREQ'} = $ggfreq;
					$hcnv_gather->{$global_id}->{'GOLD_L_FREQ'} = $glfreq;
					$hcnv_gather->{$global_id}->{'OMIN_MG'} = $mg;
					$hcnv_gather->{$global_id}->{'dbVar_status'} = $dbVar_status;
					$hcnv_gather->{$global_id}->{'RANKAnnot'} = $rankannot;
					$hcnv_gather->{$global_id}->{'GT'}=$gt;
					$hcnv_gather->{$global_id}->{'CN'}=$cn;
					$hcnv_gather->{$global_id}->{'RATIO'}=$ratio/scalar(@$tab_ind_id);			#on prend la moyenne
					$hcnv_gather->{$global_id}->{'QUAL'}=$qual/scalar(@$tab_ind_id);			#on prend la moyenne
					
			}
		}
	}
	return $hcnv_gather;
}


sub control_object {
	my ($h) = @_;
	my @keys = ('id','SVTYPE','CHROM','START','END','SVLEN','GT','CN','RATIO','QUAL','ELEMENTARY');
	foreach my $k (@keys){
		die($k) unless exists $h->{$k};
	}
	
}


sub parse_sniffles_bnd {
	 my ($hBND) = @_;

	foreach my $gid (keys %{$hBND}){
		my @mate = values %{$hBND->{$gid}};
		confess(Dumper $hBND->{$gid}) if scalar(@mate) ne 1;
		my $h1 = dclone $mate[0];
		my $hchr;
		my $alt_field = $h1->{alt};
		if ($alt_field =~ /[\[\]](?:chr)?([A-Za-z0-9]+):(\d+)[\[\]]/i) {
   		 $hchr->{chr} = $1;  # Capture le chromosome
    	$hchr->{pos} = $2;    # Capture la position
    	 $hchr->{chr} =~ s/chr//;
   		  $hchr->{chr} ="MT" if $hchr->{chr} eq "M";
		} else {
			
 	   die()
	}
	$hchr->{infos} = $h1->{infos};
	$hchr->{qual} = $h1->{qual};
	$hchr->{alt} = "[".$h1->{chr}.":".$h1->{pos}."[";
	$hBND->{$gid}->{$gid."_2"} = $hchr; 
	}
	return;
}

sub parse_vcf_bnd {
 my ($patient,$caller) = @_;
  my $project = $patient->project;
	 my $vcf_file = $patient->getSVFile($caller);
	my $vcf = Bio::DB::HTS::VCF->new( filename => "$vcf_file" );
	my $header = $vcf->header();
	my $hBND;
	while (my $row = $vcf->next){
			
			
			my $type =  get_value($row->get_info($header, "SVTYPE"));
			next if  $type ne "BND";
			my $event_id = $row->id();
			my $mate_id = get_value($row->get_info($header, "MATEID"));
			my @tt =  sort($event_id,$mate_id);
			my $g_id = join(":",@tt);
			my $chr = $project->getChromosome($row->chromosome($header),1);
		
			
			next unless $chr; 
			next  unless  $row->has_filter($header,"PASS");
			$hBND->{$g_id}->{$event_id}->{qual} = $row->quality;
			$hBND->{$g_id}->{$event_id}->{qual} = 0 if  $row->quality eq "NaN";
			$hBND->{$g_id}->{$event_id}->{chr} = $chr->name ;
			$hBND->{$g_id}->{$event_id}->{alt} = $row->get_alleles()->[0];
			$hBND->{$g_id}->{$event_id}->{pos} = $row->position + 1;
			$hBND->{$g_id}->{$event_id}->{infos} = $row->get_format($header);
			
	}
	
	
	$vcf->close();
	parse_sniffles_bnd($hBND) if $caller eq "Sniffles2" or $caller eq "Spectre";
	
	return create_transloc_hash($hBND,$patient);
}

sub create_transloc_hash {
	my ($hBND,$patient) = @_;
		my $hTransLoc;
	foreach my $gid (keys %{$hBND}){
		my @mate = values %{$hBND->{$gid}};
		next if scalar(@mate) ne 2;
		my ($b1,$b2) = sort {$a->{chr} <=> $b->{chr} or $a->{pos} <=> $b->{pos}} @mate;
		my $event_id = $b1->{chr}."_".$b1->{pos}."_".$b2->{chr}."_".$b2->{pos};
			$hTransLoc->{$event_id}->{"ID"}= $event_id;
			$hTransLoc->{$event_id}->{"TRANSLOC"}= $b1->{chr}."to".$b2->{chr};
			$hTransLoc->{$event_id}->{"QUAL"}= $b1->{qual};
			$hTransLoc->{$event_id}->{"TYPE"}= "TL";
			$hTransLoc->{$event_id}->{"TYPE"}= "INV" if $b1->{chr} eq $b2->{chr};
			# info du chromosome de départ
			$hTransLoc->{$event_id}->{"CHROM1"}= $b1->{chr};
			$hTransLoc->{$event_id}->{"POS1"}= $b1->{pos};
			#$hTransLoc->{$event_id}->{"CYTOBAND1"}= $b1->{cytoband};
			$hTransLoc->{$event_id}->{"ALT1"} = $b1->{alt};
			$hTransLoc->{$event_id}->{INFOS} = delete $b1->{infos};
			# info du chromosome d'arrivée
			$hTransLoc->{$event_id}->{"CHROM2"}=  $b2->{chr};
			$hTransLoc->{$event_id}->{"POS2"}=  $b2->{pos};
			$hTransLoc->{$event_id}->{"ALT2"}= $b2->{alt};
			$hTransLoc->{$event_id}->{PATIENT_ID}= $patient->id;
			$hTransLoc->{$event_id}->{PATIENT_NAME}= $patient->name;
			#push(@$transloc,$hTransLoc);
	}
	return ($hTransLoc);
}

1;
