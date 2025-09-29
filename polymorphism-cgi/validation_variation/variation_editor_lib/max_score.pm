package max_score;
use strict;
use Data::Dumper;
use Time::HiRes qw( time);
use List::Util qw( max);





sub calculate {
	my ($project,$patient, $list,$final_polyviewer_all,$rocksdb_pv) = @_;
	$project->buffer->close_lmdb();
	my $mother = $patient->getFamily->getMother();
	my $father = $patient->getFamily->getFather();


	my $hno;
	my $tsum = 0;
	my $t    = time;
	my $xp   = 0;

	my $total_time = 0;
	my $ztotal     = 0;

	my $h_all_variants_validations = $patient->getProject->validations();
	my $no_dude                    = $patient->getGenesDude();
	$no_dude = undef if $project->isGenome();
	my @ids = map{$_->{id}} @{$list};
	my $xx =time;
	$final_polyviewer_all->prepare(\@ids);
	warn abs($xx - time)." final polyviewer";
	#warn "\t score :".abs(time -$t)." ".scalar(@ids);
	my $st =0;
	my $st2;
	

	foreach my $hgene (@$list) {
		my $gid = $hgene->{id};
		my $debug ;
		
		#my $gene = $project->newGenes( [$gid] )->[0];
		my ( $n, $chr_name ) = split( "_", $gid );
		#my $chr = $gene->getChromosome();
		$hgene->{chr_name} = $chr_name;
		
		if ( $no_dude && -e $no_dude->filename ) {
			$hgene->{level_dude} = $no_dude->get($gid);
		}
		else {
			$hgene->{level_dude} = -1;
		}
		my $t1 = time;
		my $global_gene = $final_polyviewer_all->get_cached($gid); #$chr->get_polyviewer_genes($patient,$gid);
		$st += abs(time-$t1);
		$t1 = time;
		foreach my $k ( keys %{$global_gene} ) {
			next if $k eq "penality";
			next if $k eq "denovo_rare";
			$hgene->{$k} = $global_gene->{$k};
		}
		
		#$hgene->{score} = $gene->score;
		my $class;
		$class->{biallelic} = [];
		$class->{mother}    = [];
		$class->{father}    = [];
		my $chromosome = $hgene->{chr_name};
		foreach my $k ( keys %{ $hgene->{all_variants} } ) {
		#	if ($version_db){
		#		$hgene->{score} += update_score_clinvar($k);
		#	}
		#   my $pub = $db->get_with_sequence($self->start,$self->alternate_allele);
			if ( exists $h_all_variants_validations->{ $gid . '!' . $k } ) {
				my $score_validation = $h_all_variants_validations->{ $gid . '!' . $k }->[0]->{validation};
				$hgene->{score} += 0.5 if ( $score_validation == 3 );
				$hgene->{score} += 2   if ( $score_validation == 4 );
				$hgene->{score} += 3   if ( $score_validation == 5 );
			}
			
			push(@{ $class->{biallelic} },$hgene->{all_variants}->{$k}->{score}) if exists $hgene->{all_variants}->{$k}->{biallelic};
			push( @{ $class->{mother} }, $hgene->{all_variants}->{$k}->{score} ) if exists $hgene->{all_variants}->{$k}->{mother};
			push( @{ $class->{father} }, $hgene->{all_variants}->{$k}->{score} ) if exists $hgene->{all_variants}->{$k}->{father};
		}
		
		if (scalar( @{ $class->{mother} } ) > 0 && scalar( @{ $class->{father} } ) == 0 && exists $hgene->{father}->{id} ) {
			
			my $nid = $hgene->{father}->{id};
#			$hgene->{all_variants}->{$nid}->{father} = 1;
			$hgene->{all_variants}->{$nid}->{score} = $hgene->{father}->{score};
			$hgene->{all_variants}->{$nid}->{added}++;
			my $id = $chromosome."!".$hgene->{father}->{vector_id};
			$rocksdb_pv->add_index($chromosome,$id);
			$hgene->{all_vector_ids}->{$nid} = $hgene->{father}->{vector_id};
			push( @{ $class->{father} }, ( $hgene->{father}->{score} - 2 ) );

		}
		elsif (scalar( @{ $class->{father} } ) > 0 && scalar( @{ $class->{mother} } ) == 0 && exists $hgene->{mother}->{id} )
		{
			my $nid = $hgene->{mother}->{id};
			$hgene->{all_variants}->{$nid}->{mother} = 1;
			$hgene->{all_variants}->{$nid}->{score} = $hgene->{mother}->{score};
			$hgene->{all_vector_ids}->{$nid} = $hgene->{mother}->{vector_id};
			push( @{ $class->{mother} }, ( $hgene->{mother}->{score} - 2 ) );
			my $id = $chromosome."!".$hgene->{mother}->{vector_id};
			$rocksdb_pv->add_index($chromosome,$id);
			$hgene->{all_variants}->{$nid}->{added}++;
		}

		my $score_father    = max( @{ $class->{father} } );
		my $score_mother    = max( @{ $class->{mother} } );
		$score_mother = 0 if $score_mother <0;
		$score_father = 0 if $score_father <0;
		my $score_biallelic = max( @{ $class->{biallelic} } );
		if ( $score_father + $score_mother > $score_biallelic && $patient->getFamily->isTrio()) {
			$hgene->{max_score} = $hgene->{score} + $score_father + $score_mother;
		}
		else {
			warn $hgene->{score} ." ". $score_biallelic." ".$hgene->{penality} if $hgene->{name} eq "ABCC8";
			#die()  if $hgene->{name} eq "ABCC8";
			$hgene->{max_score} = $hgene->{score} + $score_biallelic + $hgene->{penality} ;
		}
		$st2 += abs(time -$t1);
	}
	#warn "cache :=> ".$st."  ".$st2;
	return $list;
#	$project->buffer->close_lmdb();
	
}




1;