package parse_pbsv;
use strict;
use Data::Dumper;
use Storable qw(dclone);


sub factorial {
    my ($n) = @_;
    return 1 if $n == 0;
    my $fact = 1;
    $fact *= $_ for (1 .. $n);
    return $fact;
}

# Fonction pour le coefficient binomial "n choose k"
sub binomial_coeff {
    my ($n, $k) = @_;
    return 0 if $k > $n;
    return factorial($n) / (factorial($k) * factorial($n - $k));
}

# Fonction pour la probabilité binomiale
sub binomial_prob {
    my ($n, $k, $p) = @_;
    my $q = 1 - $p;
    return binomial_coeff($n, $k) * ($p ** $k) * ($q ** ($n - $k));
}

# Fonction pour la p-value (P(X ? k))
sub p_value_binomial_ge {
    my ($n, $k, $p) = @_;
    my $sum = 0;
    for my $i ($k .. $n) {
        $sum += binomial_prob($n, $i, $p);
    }
    return $sum;
}

sub parse_cnv {
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
		next unless ( ($h->{'SVTYPE'} eq "DUP") || ($h->{'SVTYPE'} eq "DEL") );
	
		$h->{'CHROM'}=$chr->name;
		$h->{"CALLER"} = $caller;
		$h->{'END'} = $row->get_info($header, "END")->[0];
		$h->{'START'} = $row->position();
		$h->{'SVLEN'} =  abs($h->{'END'} - $h->{'START'});
		$h->{'KARYOTYPE_ID'}= $chr->karyotypeId;
		
		
		$h->{'GT_a'} = $row->get_format($header, "GT");
		
		$h->{'GT'} = "0/1";
		$h->{'GT'} = "1/1" if ($h->{'GT_a'}->[0] == $h->{'GT_a'}->[1]);  
		$h->{'CN'} ="-" ;
		$h->{'RATIO'} ="-" ;
		if ($h->{'SVTYPE'} eq "DUP"){
		$h->{'CN'} = 3;
		$h->{'CN'} = 4 if $h->{'GT'} eq "1/1"; 
		}
		else {
			$h->{'CN'} = 1;
			$h->{'CN'} = 0 if $h->{'GT'} eq "1/1"; 
		}
		my $id = $h->{'SVTYPE'}."_".$chr->name."_".$h->{'START'}."_".$h->{'END'};
		$h->{'ELEMENTARY'}= [$id];
		$h->{'INFOS'} = $row->get_format($header);
		$h->{'INFOS'}->{GT} = $h->{'GT'};
		my $dp = $h->{'INFOS'}->{DP}->[0];
		my $dp_alt = $h->{'INFOS'}->{AD}->[1];
		$h->{'INFOS'}->{SR} = [$h->{'INFOS'}->{AD}->[0],$h->{'INFOS'}->{AD}->[1]];
		$h->{'INFOS'}->{PR} =  [$h->{'INFOS'}->{AD}->[0],$h->{'INFOS'}->{AD}->[1]];
		my $error_rate = 0.01;
		my $pval = p_value_binomial_ge($dp, $dp_alt, $error_rate);
		my $pval1 =  p_value_binomial_ge(4, 2, $error_rate);
		# Calcul Phred
		my $phred = -10 * log($pval) / log(10);
		$phred = 999 if ($phred eq "NaN");
		$h->{id}=$id;
		$h->{'REAL_CALLER'} = $caller;
		my $type = $h->{'SVTYPE'};
		
		$h->{'QUAL'} = $phred ;
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
sub control_object {
	my ($h) = @_;
	my @keys = ('id','SVTYPE','CHROM','START','END','SVLEN','GT','CN','RATIO','QUAL','ELEMENTARY');
	foreach my $k (@keys){
		die($k) unless exists $h->{$k};
	}
	
}


1;
