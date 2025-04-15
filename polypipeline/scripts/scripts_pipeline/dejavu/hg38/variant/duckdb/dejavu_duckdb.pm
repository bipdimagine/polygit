package dejavu_duckdb;

use FindBin qw($Bin);
use strict;
use JSON;
use Archive::Tar;
use Fcntl ':flock';
use File::NFSLock qw(uncache);
use List::Util qw(shuffle);
use Data::Dumper;
use MIME::Base64;




sub save_csv {
	my ($chr,$snps,$dir_tmp) = @_;
	my $project = $chr->project;
	my $filename = "$dir_tmp/".$chr->project->name."_".$chr->name.".csv";
	
	my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
	
	my $fh;		
	open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
	$csv->print($fh, ["project","chr38","pos38","chr19","pos19","allele","max_ratio","max_dp","transmissions","he","ho","patients","dp_ratios"]); 
	$csv->print($fh, [0,"NONE",-1,"NONE",-1,"W",0,0,"z",0,0,"NONE","NONE"]); 
	$csv->print($fh, [0,$chr->name,-1,$chr->name,-1,"W",0,0,"z",0,0,"NONE","NONE"]); 
	
	unless ($snps){
		close($fh);
		return "\'".$filename."\'"  ;
	}
	my $mt;
	if ($chr->name eq "MT"){
		$mt=1 if $project->getChromosome("MT")->length() == 16571 or $project->current_genome_version eq "HG38";
	}
	foreach my $vhh (values %$snps){
		my $he = $vhh->{he};
		my $ho = $vhh->{ho};
		my $chr0 = $chr->name;
		my $pos0 = $vhh->{start};
		my $poslift ="0";
		my $chrlift ="NONE";
		if (exists $vhh->{LIFT} ) {
			 if ($project->isChromosomeName($vhh->{LIFT}->{chromosome})){
			 	
				$poslift = $vhh->{LIFT}->{start};
				$chrlift = $project->getChromosome($vhh->{LIFT}->{chromosome})->name;
				
			 }
			
		}
		
		my $pos38 = $vhh->{start};
		my $chr38 = $chr->name;
		my $pos19 = $poslift ;
		my $chr19 = $chrlift ;
		if ($project->current_genome_version eq "HG19"){
			$pos19 = $vhh->{start};
			$chr19 = $chr->name;
			$pos38 = $poslift ;
			$chr38 = $chrlift ;
			
		}
		
		if ($mt == 1 ) {
				$pos19 = $pos38;
		}
		
		my $encoded_data = encode_base64($vhh->{value},""); 
		my $encoded_data_patients = encode_base64($vhh->{patients},""); 
		
		my ($c1,$p1) = split("!",$pos38);
		$p1 += 0;
		my ($c2,$p2) = split("!",$pos19);
		$p2 += 0;
		$csv->print($fh, [$project->id,$chr38,$pos38,$chr19,$pos19,$vhh->{allele},$vhh->{max_ratio},$vhh->{max_dp},$vhh->{transmissions},$vhh->{he},$vhh->{ho},$encoded_data_patients,$encoded_data]);
	}
	close($fh);
	
	return "\'".$filename."\'";
	
}	


my $h_models_ids = {
    solo          => 1 << 0,  # 2^0 = 1
    father        => 1 << 1,  # 2^1 = 2
    mother        => 1 << 2,  # 2^2 = 4
    both          => 1 << 3,  # 2^3 = 8
    is_parent     => 1 << 4,  # 2^4 = 16
    recessif      => 1 << 5,  # 2^5 = 32
    dominant      => 1 << 6,  # 2^6 = 64
    denovo        => 1 << 7,  # 2^7 = 128
    strict_denovo => 1 << 8,  # 2^8 = 256
    error         => 1 << 9,  # 2^9 = 512
};

sub find_variant_model {
	my ($h_models, $vid, $patient_id) = @_;
	return $h_models_ids->{solo} if exists $h_models->{$patient_id}->{solo};
	if (exists $h_models->{$patient_id}->{'dominant'}) {
		return $h_models_ids->{dominant} if $h_models->{$patient_id}->{dominant}->contains(int($vid));
	}
	foreach my $model_name ('strict_denovo', 'denovo', 'recessif', 'father', 'mother', 'both') {
		return $h_models_ids->{$model_name} if $h_models->{$patient_id}->{$model_name}->contains(int($vid));
	}
	return $h_models_ids->{error};
}

sub get_hash_model_variant {
	my ($chr, $vector_id) = @_;
	my $project = $chr->project;
	my $hvector;
	foreach my $patient (@{$project->getPatients()}) {
		my $fam = $patient->getFamily();
		my $patient_id = $patient->id();
		if ($fam->isTrio()) {
			if ($fam->isDominant()) {
				$hvector->{$patient_id}->{dominant} = $fam->getVector_individual_dominant($chr, $patient);
			}
			$hvector->{$patient_id}->{strict_denovo} = $fam->getVector_individual_strict_denovo($chr,$patient);
			$hvector->{$patient_id}->{denovo} = $fam->getVector_individual_denovo($chr,$patient);
			$hvector->{$patient_id}->{recessif} = $fam->getVector_individual_recessive($chr,$patient);
			$hvector->{$patient_id}->{father} = $fam->getFatherVector($chr);
			$hvector->{$patient_id}->{mother} = $fam->getMotherVector($chr);
			$hvector->{$patient_id}->{both}  = $hvector->{$patient_id}->{father} & $hvector->{$patient_id}->{father};
			$hvector->{$patient_id}->{both} -= $hvector->{$patient_id}->{recessif};
		}
		else {
			$hvector->{$patient_id}->{solo} = 1;
		}
	}
	return $hvector;
}		


sub compress1 {
	my ($pid,$list1,$list2) =@_;
	$list1 = [] if $list1 ==undef;
	$list2 = [] if $list2 ==undef;
	my $compressed = pack("w*",$pid, scalar(@$list1)/3, scalar(@$list2)/3, @$list1, @$list2);
	return $compressed;
}
	