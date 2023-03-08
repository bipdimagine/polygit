package quality_check;
use FindBin qw($Bin);
use strict;
use File::stat;
use Time::localtime;
use Data::Dumper;
use colored;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);

sub mendelian_statistics {
	my ( $project, $fork ) = @_;
	my $res;
	$res  = {};
	$fork = 1 unless $fork;
	#$fork=1;
	my $pm = new Parallel::ForkManager($fork);
	my $results;
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $hres ) =
			  @_;
			unless ( defined($hres) or $exit_code > 0 ) {
				print
				  qq|No message received from child process $exit_code $pid!\n|;
				warn Dumper $hres;
				return;
			}

			#warn $hres->{data};
			die() unless $hres->{data};
			push( @{ $results->{data} }, @{ $hres->{data} } );

		}
	);
	map { $_->getMembers } @{ $project->getFamilies };
	$project->buffer->dbh_deconnect();

	my $pid;
	foreach my $f ( @{ $project->getFamilies } ) {
		$f->getMembers();

		$pid = $pm->start and next;

		$project->buffer->dbh_reconnect();
		#my $ps = $f->getMembers();
		#	next if scalar(@$ps) > 1;
		my $vcf = concatVcf( $project, $f );

		my $hres;
		$hres->{data} = fast_plink( $project, $vcf, $f );

		#warn Dumper $hres->{data};
		$pm->finish( 0, $hres );
	}    #end for range range
	$pm->wait_all_children();

	#	die();
	add_columns(
		$results,
		[
			"results", "familly",   'sample', "sex",
			"SRY",     "plink_sex", "mendelian"
		]
	) if $project->isFamilial;
	add_columns( $results,
		[ "results", "familly", 'sample', "sex", "SRY", "plink_sex" ] )
	  unless $project->isFamilial;
	return $results;
}

sub mendelian_statistics2 {
	my ($project) = @_;

	my $vcf = concatVcf($project);
	warn ":::>".$vcf;
	return fast_plink( $project, $vcf );
}

sub concatVcf {
	my ( $project, $fam ) = @_;

	my $bcftools = $project->buffer->software("bcftools");
	print "\n\n### Step 2: concat VCF files.\n";
	my $dir_out =
	  $project->getCallingPipelineDir("unifiedgenotyper") . '/plink';
	if ( -e $dir_out ) {

		#system("rm $dir_out/*");
	}
	else {
		system("mkdir -p $dir_out");
	}

	my $fileout =  $dir_out . "/" . $project->name . "." . $fam->name . ".merge.vcf";
	my $fileout1 =
	  $dir_out . "/" . $project->name . "." . $fam->name . ".merge1.vcf";
	my $bed = $dir_out . "/" . $project->name . ".bed";
	#return $fileout if -e $fileout;
	unlink $fileout if -e $fileout;
	unlink $bed     if -e $bed;
	my @vcfs;

	foreach my $p ( @{ $fam->getMembers } ) {
		my $methods = $p->getCallingMethods();
		my $m;
		if ( scalar(@$methods) == 1 ) {
			($m) = $methods->[0];

		}
		else {
			($m) = grep { $_ eq "unifiedgenotyper" } @$methods;
			($m) = grep { $_ =~ /haplotypecaller4/ } @$methods unless $m;
			($m) = grep { $_ =~ /dragen_calling/ } @$methods unless $m;
			($m) = grep { $_ =~ /dragen-calling/ } @$methods unless $m;
			($m) = grep { $_ =~ /haplotypecaller/ } @$methods unless $m;
			($m) = grep { $_ =~ /mpileup/ } @$methods unless $m;
			($m) = grep { $_ =~ /dibayes/ } @$methods unless $m;
			($m) = grep { $_ =~ /gatk/ } @$methods unless $m;
			($m) = grep { $_ =~ /p1_freebayes/ } @$methods unless $m;
		}
		unless ($m) {
			warn "\n\n";
			warn "No methods calling for patient " . $p->name();
			warn Dumper $methods;
			die;
		}
		my $vcfs = $p->getVariationsFiles($m);
		my $vcf;
		$vcf = $p->getVariationsFile($m);
		die( $vcf . "-" . $p->name ) unless -e $vcf;
		push( @vcfs, $vcf );
	}
	if ( scalar(@vcfs) == 1 ) {
		return $vcfs[0];
	}
	## create bed file
#my $variant_string = join(" --variant ",@vcfs);
#my $cmd = $gatk." -R ".$project->getGenomeFasta()." -nt 8 -T CombineVariants  --variant $variant_string -o $fileout "." -L ". $bed." 2>/dev/null >/dev/null" ;
	my $variant_string = join( " ", @vcfs );

#my $cmd = " $bcftools merge   $variant_string -o $fileout --force-samples >$fileout " ;
	my $cmd = " $bcftools merge   $variant_string -o $fileout  >$fileout 2>/dev/null";
	system($cmd) unless -e $fileout;
	warn "end";
	return $fileout;
}

sub fast_plink {
	my ( $project, $vcf_file, $fam ) = @_;
	my $statistics;
	my $plink       = $project->buffer->software("plink");
	my $projectName = $fam->name();

	print "\n\n### Step 3: launching PLINK.\n";
	my $dir = $project->getCallingPipelineDir("plink");

	my $logPlink = $dir . $fam->name . '.plink.resume';
	my @snps;

	if ( $vcf_file =~ /gz/ ) {
		open( VCF, "zcat $vcf_file | " );
	}
	else {
		open( VCF, "cat $vcf_file | " );
	}
	my @samples;
	my $DP;
	my $AD;
	my $GT;
	my %max_chr;
	while ( my $line = <VCF> ) {
		next if $line =~ /##/;
		chomp($line);
		if ( $line =~ /#CHR/ ) {
			my @tdata = split( " ", $line );
			@samples = splice( @tdata, 9 );    #$tdata[9..-1];
			next;
		}
		my (
			$chrom, $pos,    $pid,  $ref,    $alt,
			$qual,  $filter, $info, $format, @gsamples
		) = split( " ", $line );
		next if length($ref) > 1;
		next if length($alt) > 1;
		$max_chr{$chrom}++;
		next if ( $max_chr{$chrom} > 80_000 );
		next if $ref =~ /,/;
		next if $alt =~ /,/;
		next if scalar(@gsamples) ne scalar(@samples);

		my $snp;
		my $id = $chrom . "_" . $pos . "_" . $alt;
		$snp->{id}       = $id;
		$snp->{position} = $pos;
		$chrom =~ s/chr//;
		my $debug;
		$debug = 1    if $id eq "1_754964_T";
		$chrom = "MT" if $chrom eq "M";
		$snp->{chr} = $chrom;

		unless ($DP) {
			my @t = split( ":", $format );
			for ( my $i = 0 ; $i < @t ; $i++ ) {
				$DP = $i if $t[$i] eq "DP";
			}
			$DP = -1 unless $DP;
		}
		unless ( $AD or $GT ) {
			my @t = split( ":", $format );
			for ( my $i = 0 ; $i < @t ; $i++ ) {
				$AD = $i if $t[$i] eq "AD";
				$GT = $i if $t[$i] eq "GT";
			}
			$AD = -1 unless $AD;
		}
		for ( my $i = 0 ; $i < @samples ; $i++ ) {
			my $name = $samples[$i];

			#warn Dumper @string;
			#die();
			my @tt = split( ":", $gsamples[$i] );
			next unless $tt[$GT];
			my @string = split( "", $tt[$GT] );
			if ( @string == 1 ) {
				$string[1] = "/";
				$string[2] = $string[0];
			}

			#else { 	warn scalar(@gsamples) ;die($line." ".@string." ".$name); }
			my $a;
			my $b;
			if ( $AD > 0 ) {

				#my @tt = split(":",$string);
				my $v = $tt[$AD];
				( $a, $b ) = split( ",", $v );

				#$geno = "0 0" if $v < 10;
			}

			#my $v = $tt[$DP];
			#$geno = "0 0" if $v < 10;
			my $geno;
			if ( $a + $b < 10 ) {
				$geno = "0 0";
			}
			elsif ( $string[0] eq "0" && $string[2] eq "0" ) {
				$geno = "$ref $ref";
			}
			elsif ( $string[0] eq "0" && $string[2] eq "1" ) {
				$geno = "$ref $alt";
				$geno = "0 0" if ( $b / ( $a + $b ) ) < 0.30;
			}
			elsif ( $string[0] eq "1" && $string[2] eq "0" ) {
				$geno = "$ref $alt";
				$geno = "0 0" if ( $b / ( $a + $b ) ) < 0.30;
			}
			elsif ( $string[0] eq "1" && $string[2] eq "1" ) {
				$geno = "$alt $alt";
				$geno = "0 0" if ( $a + $b ) < 20;
			}
			elsif ( $string[0] eq "." && $string[2] eq "." ) {
				$geno = "0 0";
				my $patient = $project->getPatient($name);
				my $d       = $patient->depth( $chrom, $pos, $pos );
				$geno = "0 0";
				$geno = "$ref $ref" if ( $d->[0] > 30 );
				warn "$name " . $geno . " " . $d->[0] if $debug;

				#warn $d->[0] if ($d->[0] < 20);
			}

			#elsif ($string =~ /.:./) { $geno ="$ref $ref"; }
			else {
				warn scalar(@gsamples);
				die( $line . " " . @string . " " . $name );
			}
			if ( $DP > 0 ) {

				#my @tt = split(":",$string);
				#my $v = $tt[$DP];
				#$geno = "0 0" if $v < 10;

			}
			warn $name . " " . $geno if $debug;
			$snp->{samples}->{$name} = $geno;
		}

		#die() if $debug;
		#foreach my $f (values %{$project->families()}) {
		my ($find) =
		  grep { $snp->{samples}->{$_} eq "0 0" } @{ $fam->{members} };

		if ($find) {

			map { $snp->{samples}->{$_} = "0 0" } @{ $fam->{members} };
		}
		else {

			$fam->{nb_snp}++;
		}

		#}

		#my ($find) = grep {$snp->{samples}->{$_}  eq "0 0"} @samples;
		#next if $find;
		push( @snps, $snp );
	}
	$dir .= "/plink/";
	mkdir $dir unless -e $dir;
	my $ped_file = $dir . "/" . $fam->name . ".ped";
	my @samples_name;
	open( PED, ">$ped_file" );

	foreach my $p ( @{ $fam->getMembers } ) {
		print PED $p->pedigreeLine() . "\n";
		push( @samples_name, $p->name );

	}
	close PED;

	die() unless -e $ped_file;

	my $tped_file = $dir . "/" . $fam->name . ".tped";
	open( TPED, ">" . $tped_file );
	foreach my $snp (@snps) {
		my @data;
		push( @data, $snp->{chr} );
		push( @data, $snp->{id} );
		push( @data, 0 );
		push( @data, $snp->{position} );
		foreach my $name (@samples_name) {
			if ( exists $snp->{samples}->{$name} ) {
				push( @data, $snp->{samples}->{$name} );
			}
			else { push( @data, $snp->{samples}->{$name}, "0 0" ); }
		}
		print TPED join( "\t", @data ) . "\n";
	}
	close TPED;
	warn $tped_file;
	if ( @{ $fam->getParents } ) {
		my $cmd2 =
"$plink --tped $tped_file --tfam $ped_file --noweb --mendel  --mendel-duos --out $dir/$projectName";

		#	warn $cmd2;
		my @log   = `$cmd2`;
		my $file1 = "$dir/" . $fam->name();
		open( MENDEL, "$file1.imendel" )
		  || die("problem with $file1.imendel \n$cmd2");

		while ( my $line = <MENDEL> ) {
			chomp($line);

			#next if $line =~ /FID\s*IID\s*N/;
			my ( $f, $name, $nb ) = split( " ", $line );
			next if $f eq "FID" && $name eq "IID" && $nb eq "N";
			my $nb_snp = $fam->{nb_snp};
			next if $nb_snp == 0;
			my $patient = $project->getPatient($name);
			my $fam     = $patient->getFamily();
			my $p       = int( $nb / $nb_snp * 10000 ) / 100;
			warn "------------>" . $fam->name . " " . $p;
			$statistics->{ $fam->name }->{$name}->{mendelian_errors}->{percent}
			  = $p;
			$statistics->{ $fam->name }->{$name}->{mendelian_errors}->{nb_snp}
			  = $project->families()->{$f}->{nb_snp};
		}
		close(MENDEL);
	}

	my $cmd3 =
"$plink --tped $tped_file --tfam $ped_file  --noweb  --check-sex --out $dir/$projectName 1>$logPlink 2>>$logPlink";
	system($cmd3);

	#	warn $cmd3;
	#`$cmd3`;
	my $sexFile = $dir . '/' . $projectName . '.sexcheck';

#	unless (-e $sexFile) { die("\n\nERROR: PLINK sex file doesn't exist... Die...\n\n") };

	if ( -e $sexFile ) {
		open( SEX, "$sexFile" );
		warn $sexFile;
		while (<SEX>) {
			chomp($_);
			my $line = $_;
			next if ( $line =~ /PEDSEX/ );
			$line =~ s/\s+/ /g;
			my @lFields = split( " ", $line );
			my $name    = $lFields[1];
			my $patient = $project->getPatient($name);
			my $fam     = $patient->getFamily();
			$statistics->{ $fam->name }->{$name}->{check_sex}->{plink}->{sex} =
			  $lFields[3];

			#$statistics->{$name}->{check_sex}->{plink}->{snp} = $lFields[2];
			$statistics->{ $fam->name }->{$name}->{check_sex}->{plink}->{status}
			  = $lFields[4];
			$statistics->{ $fam->name }->{$name}->{check_sex}->{plink}->{score}
			  = int( $lFields[5] * 100 ) / 100;
		}
		close(SEX);
	}
	else {
		foreach my $p ( @{ $fam->getMembers } ) {
			my $name = $p->name();
			my $fam  = $p->getFamily();
			$statistics->{ $fam->name }->{$name}->{check_sex}->{plink}->{sex} =
			  0;
			$statistics->{ $fam->name }->{$name}->{check_sex}->{plink}->{status}
			  = "?";
			$statistics->{ $fam->name }->{$name}->{check_sex}->{plink}->{score}
			  = "?";
		}
	}
	warn "end";
	return ( table_json_mendelian( $project, $statistics, $fam ) );
}

sub add_columns {
	my ( $h, $names ) = @_;
	foreach my $name (@$names) {
		push( @{ $h->{columns} }, { field => $name, title => $name } );
	}
}

sub transform_array_to_json_like {
	my ($resume) = @_;
	my $hresume = {};
	add_columns( $hresume, $resume->{header} );
	foreach my $line ( @{ $resume->{lines} } ) {
		my $z = 0;
		my $hline;
		foreach my $id ( @{ $resume->{header} } ) {
			$hline->{$id} = $line->[$z];
			$z++;
		}
		push( @{ $hresume->{data} }, $hline );
	}
	return $hresume;

}

sub files {
	my ($project) = @_;
	my $patients = $project->getPatients();
	my $resume;
	my $hresume;
	$hresume->{header} = [ "patients", "bam" ];
	my $lines;
	my %hmethods;
	my $methods = $project->getCallingMethods();
	push( @{ $hresume->{header} }, @$methods );

	foreach my $p (@$patients) {
		my $methods = $p->getCallingMethods();
		map { $hmethods{$_}++ } @$methods;
		my $line;
		my $hline;
		my $bam       = $p->getBamFile();
		my $timestamp = ctime( stat($bam)->mtime );
		$hline->{"patients"} = { text => $p->name, type => "default" };
		push( @$line, { text => $p->name, type => "default" } );
		my $col_bam = { text => "NONE", type => "error" };

		if ( -e $bam ) {
			my $t = utility::return_date_from_file($bam);
			$col_bam = { text => "$t", type => "success" };
		}
		push( @$line, $col_bam );
		$hline->{"bam"} = $col_bam;

		#push(@{$errors->{file}},"NO BAM FOR : ".$p->name()) unless -e $bam;

		my $vcf_vars = $p->getVariationsFiles();
		my $files    = $p->callingFiles();
		foreach my $method (@$methods) {
			my $col_method = { text => "NONE", type => "danger" };
			if ( exists $files->{variations}->{$method} ) {

				my $t = utility::return_date_from_file(
					$files->{variations}->{$method} );
				$col_method = { text => "$t", type => "success" };

			}
			push( @$line, $col_method );
			$hline->{$method} = $col_method;
		}
		push( @{ $resume->{data} }, $hline );

		push( @{ $hresume->{lines} }, $line );
	}
	return transform_array_to_json_like($hresume);
	return $hresume;
}

sub statistics_variations_patient {
	my ($p) = @_;
	my $hnb;
	$hnb->{ho}     = $p->countHomozygote();
	$hnb->{he}     = $p->countHeterozygote();
	$hnb->{snp}    = $p->countSubstitutions();
	$hnb->{indel}  = $p->countIndels();
	$hnb->{total}  = $p->countVariations();
	$hnb->{public} = $p->countPublicVariations();
	return $hnb;

}

sub statistics_variations {
	my ( $project, $fork ) = @_;
	my $resume;
	my $patients = $project->getPatients();
	warn $project;
	my $h = fork_patients( $project, \&statistics_variations_patient, $fork );
	push(
		@{ $resume->{header} },
		( "patients", "snp", "indel", "%he", "%public" )
	);
	my $sum;
	foreach my $k ( keys %$h ) {
		$sum->{snp}   += $h->{$k}->{snp};
		$sum->{indel} += $h->{$k}->{indel};
	}
	my $nbp = scalar( @{$patients} );
	my $mean;
	$mean->{indel} = int( $sum->{indel} / $nbp );
	$mean->{snp}   = int( $sum->{snp} / $nbp );
	foreach my $p (@$patients) {
		my $n = $p->name();
		my $line2;
		push( @$line2, { text => $p->name, type => "default" } );
		my $color = "success";

		#my $d = $hnb->{snp}->{$p->name}/abs();
		if (   $h->{$n}->{snp} > $mean->{snp} * 1.5
			or $h->{$n}->{snp} < $mean->{snp} * 0.75 )
		{
			$color = "warning";
		}
		if (   $h->{$n}->{snp} > $mean->{snp} * 2
			or $h->{$n}->{snp} < $mean->{snp} * 0.5 )
		{
			$color = "danger";
		}
		push( @$line2, { text => $h->{$n}->{snp}, type => "$color" } );
		$color = "success";

		#my $d = $hnb->{snp}->{$p->name}/abs();
		if (   $h->{$n}->{indel} > $mean->{indel} * 1.25
			or $h->{$n}->{indel} < $mean->{indel} * 0.75 )
		{
			$color = "warning";
		}
		if (   $h->{$n}->{indel} > $mean->{indel} * 1.5
			or $h->{$n}->{indel} < $mean->{indel} * 0.5 )
		{
			$color = "danger";
		}
		push( @$line2, { text => $h->{$n}->{indel}, type => "$color" } );
		$h->{$n}->{total} = -1 if $h->{$n}->{total} == 0;
		$h->{$n}->{total} = -1 unless $h->{$n}->{total};
		my $z    = int( ( $h->{$n}->{he} / $h->{$n}->{total} ) * 100 );
		my $type = "danger";

		if ( $z < 75 ) {
			$type = "success";
		}
		elsif ( $z <= 50 ) {
			$type = "danger";
		}
		push( @$line2, { text => "$z%", type => "$type" } );

		#CAS des temoins negatifs
		my $p;
		if ( int( $h->{$n}->{total} ) == 0 ) { $p = 0; }
		else { $p = int( ( $h->{$n}->{public} / $h->{$n}->{total} ) * 100 ); }
		if ( $p >= 90 ) {
			push( @$line2, { text => "$p%", type => "success" } );
		}
		elsif ( $p >= 85 ) {
			push( @$line2, { text => "$p%", type => "success" } );
		}
		elsif ( $p >= 75 ) {
			push( @$line2, { text => "$p%", type => "warning" } );
		}
		else {
			push( @$line2, { text => "$p%", type => "danger" } );
		}

		push( @{ $resume->{lines} }, $line2 );
	}
	warn "end statistics variations";
	return ( transform_array_to_json_like($resume) );

}

sub statistics_variations2 {
	my ($project) = @_;
	my $patients = $project->getPatients();
	my $resume;
	my $hnb;
	my $sum;

	foreach my $p (@$patients) {
		$hnb->{ho}->{ $p->name }     = $p->countHomozygote();
		$hnb->{he}->{ $p->name }     = $p->countHeterozygote();
		$hnb->{snp}->{ $p->name }    = $p->countSubstitutions();
		$hnb->{indel}->{ $p->name }  = $p->countIndels();
		$hnb->{total}->{ $p->name }  = $p->countVariations();
		$hnb->{public}->{ $p->name } = $p->countPublicVariations();
		$sum->{snp}   += $hnb->{snp}->{ $p->name };
		$sum->{indel} += $hnb->{indel}->{ $p->name };

		#warn $p->countPublicVariations()." ".$hnb->{total}->{$p->name};
	}
	push(
		@{ $resume->{header} },
		( "patients", "snp", "indel", "%he", "%public" )
	);
	my $nbp = scalar( @{$patients} );
	my $mean;
	$mean->{indel} = int( $sum->{indel} / $nbp );

	$mean->{snp} = int( $sum->{snp} / $nbp );

	foreach my $p (@$patients) {
		my $line2;
		push( @$line2, { text => $p->name, type => "default" } );
		my $color = "success";

		#my $d = $hnb->{snp}->{$p->name}/abs();
		if (   $hnb->{snp}->{ $p->name } > $mean->{snp} * 1.25
			or $hnb->{snp}->{ $p->name } < $mean->{snp} * 0.75 )
		{
			$color = "warning";
		}
		if (   $hnb->{snp}->{ $p->name } > $mean->{snp} * 1.5
			or $hnb->{snp}->{ $p->name } < $mean->{snp} * 0.5 )
		{
			$color = "danger";
		}
		push( @$line2,
			{ text => $hnb->{snp}->{ $p->name }, type => "$color" } );
		$color = "success";

		#my $d = $hnb->{snp}->{$p->name}/abs();
		if (   $hnb->{indel}->{ $p->name } > $mean->{indel} * 1.25
			or $hnb->{indel}->{ $p->name } < $mean->{indel} * 0.75 )
		{
			$color = "warning";
		}
		if (   $hnb->{indel}->{ $p->name } > $mean->{indel} * 1.5
			or $hnb->{indel}->{ $p->name } < $mean->{indel} * 0.5 )
		{
			$color = "danger";
		}
		push( @$line2,
			{ text => $hnb->{indel}->{ $p->name }, type => "$color" } );
		my $z = int(
			( $hnb->{he}->{ $p->name } / $hnb->{total}->{ $p->name } ) * 100 );
		my $type = "success";
		if ( $z < 60 ) {
			$type = "success";
		}
		elsif ( $z <= 50 ) {
			$type = "danger";
		}
		push( @$line2, { text => "$z%", type => "$type" } );

		my $p =
		  int( ( $hnb->{public}->{ $p->name } / $hnb->{total}->{ $p->name } ) *
			  100 );
		if ( $p >= 90 ) {
			push( @$line2, { text => "$p%", type => "success" } );
		}
		elsif ( $p >= 85 ) {
			push( @$line2, { text => "$p%", type => "success" } );
		}
		elsif ( $p >= 75 ) {
			push( @$line2, { text => "$p%", type => "warning" } );
		}
		else {
			push( @$line2, { text => "$p%", type => "danger" } );
		}

		push( @{ $resume->{lines} }, $line2 );
	}
	return ( transform_array_to_json_like($resume) );

}

sub statistics_variations2 {
	my ($project) = @_;
	die();
	my $patients = $project->getPatients();
	my $resume;
	my $hnb;

	foreach my $p (@$patients) {

		foreach my $v ( @{ $p->getStructuralVariations() } ) {
			$hnb->{total}->{ $p->name }++;
			if ( $v->isVariation ) {
				$hnb->{snp}->{ $p->name }++;
				if ( $v->isPublic ) {
					$hnb->{public}->{ $p->name }++;
				}
			}
			else {
				$hnb->{indel}->{ $p->name }++;
			}
			if ( $v->isHomozygote($p) ) {
				$hnb->{ho}->{ $p->name }++;
			}
			else {
				$hnb->{he}->{ $p->name }++;
			}
		}
	}

	push(
		@{ $resume->{header} },
		( "patients", "snp", "indel", "%he", "%public" )
	);
	foreach my $p (@$patients) {
		my $line2;
		push( @$line2, { text => $p->name, type => "default" } );
		push( @$line2,
			{ text => $hnb->{snp}->{ $p->name }, type => "default" } );
		push( @$line2,
			{ text => $hnb->{indel}->{ $p->name }, type => "default" } );
		my $z = int(
			( $hnb->{he}->{ $p->name } / $hnb->{total}->{ $p->name } ) * 100 );
		my $type = "success";
		if ( $z < 60 ) {
			$type = "success";
		}
		elsif ( $z <= 50 ) {
			$type = "danger";
		}
		push( @$line2, { text => "$z%", type => "$type" } );

		my $p =
		  int( ( $hnb->{public}->{ $p->name } / $hnb->{snp}->{ $p->name } ) *
			  100 );
		if ( $p > 90 ) {
			push( @$line2, { text => "$p%", type => "success" } );
		}
		elsif ( $p > 85 ) {
			push( @$line2, { text => "$p%", type => "default" } );
		}
		elsif ( $p > 75 ) {
			push( @$line2, { text => "$p%", type => "warning" } );
		}
		else {
			push( @$line2, { text => "$p%", type => "danger" } );
		}

		push( @{ $resume->{lines} }, $line2 );
	}
	return ( transform_array_to_json_like($resume) );

}


sub check_calling_methods_in_cache {
	my ($project) = @_;
	foreach my $patient (@{$project->getPatients()}) {
		my $pname = $patient->name();
		print "\n";
		print colored::stabilo('white',"# Patient $pname",1);
		print "\n";
		my ($h_methods_patient, $h_methods_variants);
		foreach my $method (@{$patient->getCallingMethods()}) { $h_methods_patient->{$method} = undef; }
		my @lChr = @{$project->getChromosomes()};
		foreach my $chr (reverse @lChr) {
			next if $chr->not_used();
			#print "-> checking chr".$chr->id()."\n";
			my $no = $chr->get_lmdb_variations("r");
			my $cursor = $no->cursor( 1, $chr->end() );
			while ( my $var_id = $cursor->next_key ) {
				my $lmdb_index = $cursor->current_index();
				my $variation = $no->get($var_id);
				my $h_details = $variation->sequencing_details($patient);
				foreach my $method (keys %{$h_details}) {
					$h_methods_variants->{$method}++;
				}
			}
			$no->close();
			last if (scalar(keys %$h_methods_patient) == scalar(keys %$h_methods_variants));
		}
		
		
		foreach my $method (sort keys %$h_methods_patient) {
			my @lLetters = split('', $method);
			my $method_3letters = $lLetters[0].$lLetters[1].$lLetters[2];
			if (exists $h_methods_variants->{$method_3letters}) {
				my $text = $method;
				print colored::stabilo('green',$text,1);
				print ' => found '.$method_3letters.' => found nb '.$h_methods_variants->{$method_3letters};
				print "\n";
				delete $h_methods_variants->{$method_3letters};
			}
			else {
				my $text = $method;
				print colored::stabilo('red',$text,1);
				print ' => not found '.$method_3letters;
				print "\n";
			}
		}
		if (scalar keys %$h_methods_variants > 0) {
			print "\n";
			print colored::stabilo('red','=> found this method(s) in cache: '.sort join(', ', keys %$h_methods_variants),1); 
			print "\n";
		}
	}
}

sub transcripts_variations {
	my ($project) = @_;
}

sub public_variations {
	my ($project) = @_;
	my $patients = $project->getPatients();
	my $resume;
	my $resume2;
	my @transcripts_cgi = @{ $project->bundle_transcripts() };
	my $hvp;
	my $hv;
	my $ht;
	my $hnbv;
	my $hnbindel;
	my $hnbhe;
	my $hnbtotal;
	my $hnb;

	foreach my $t (@transcripts_cgi) {
		$ht->{$t}++;
		push( @{ $resume2->{header} }, $t );

	}
	my $hres;

	foreach my $p (@$patients) {

		foreach my $v ( @{ $p->getStructuralVariations() } ) {
			$hnb->{total}->{ $p->name }++;
			if ( $v->isVariation ) {
				$hnb->{snp}->{ $p->name }++;
				if ( $v->isPublic ) {
					$hnb->{public}->{ $p->name }++;
				}
			}
			else {
				$hnb->{indel}->{ $p->name }++;
			}
			if ( $v->isHomozygote($p) ) {
				$hnb->{ho}->{ $p->name }++;
			}
			else {
				$hnb->{he}->{ $p->name }++;
			}

			$hnbtotal->{ $p->name }++;
			foreach my $t ( @{ $v->getTranscripts } ) {
				next unless exists $ht->{ $t->name };
				$hv->{ $t->name }->{ $p->name }++;
				$hvp->{ $t->name }->{ $p->name }++ if $v->isPublic();
			}

		}

	}

	#compute stats ;
	my $tstat;
	foreach my $t (@transcripts_cgi) {
		my @values = values %{ $hv->{$t} };
		my $stat   = Statistics::Descriptive::Full->new();
		$stat->add_data(@values);
		$tstat->{$t}->{mean}               = $stat->mean;
		$tstat->{$t}->{standard_deviation} = $stat->standard_deviation;
	}
	foreach my $p (@$patients) {
		my $line;
		my $line2;

		push( @$line,  { text => $p->name, type => "default" } );
		push( @$line2, { text => $p->name, type => "default" } );
		my $nb;
		my $nbp;

		foreach my $t (@transcripts_cgi) {
			$nb  += $hv->{$t}->{ $p->name };
			$nbp += $hvp->{$t}->{ $p->name };
			my $dif  = abs( $hv->{$t}->{ $p->name } - $tstat->{$t}->{mean} );
			my $text = $hv->{$t}->{ $p->name };
			if ( $dif <= $tstat->{$t}->{standard_deviation} * 2 ) {
				push( @$line, { text => $text, type => "default" } );
			}
			elsif ( $dif <= $tstat->{$t}->{standard_deviation} * 3 ) {
				push( @$line, { text => $text, type => "warning" } );
			}
			else {
				push( @$line, { text => $text, type => "danger" } );
			}
		}

		my $p = 0;
		$p = int( ( $nbp / $nb ) * 100 ) if $nb > 0;

		push( @{ $resume2->{lines} }, $line );

	}    #for patients

	push(
		@{ $resume->{header} },
		( "patients", "snp", "indel", "%he", "%public" )
	);
	foreach my $p (@$patients) {
		my $line2;
		push( @$line2, { text => $p->name, type => "default" } );
		push( @$line2,
			{ text => $hnb->{snp}->{ $p->name }, type => "default" } );
		push( @$line2,
			{ text => $hnb->{indel}->{ $p->name }, type => "default" } );
		my $z = int(
			( $hnb->{he}->{ $p->name } / $hnb->{total}->{ $p->name } ) * 100 );
		my $type = "success";
		if ( $z < 60 ) {
			$type = "success";
		}
		elsif ( $z <= 50 ) {
			$type = "danger";
		}
		push( @$line2, { text => "$z%", type => "$type" } );

		my $p =
		  int( ( $hnb->{public}->{ $p->name } / $hnb->{snp}->{ $p->name } ) *
			  100 );
		if ( $p > 90 ) {
			push( @$line2, { text => "$p%", type => "success" } );
		}
		elsif ( $p > 85 ) {
			push( @$line2, { text => "$p%", type => "default" } );
		}
		elsif ( $p > 75 ) {
			push( @$line2, { text => "$p%", type => "warning" } );
		}
		else {
			push( @$line2, { text => "$p%", type => "danger" } );
		}

		push( @{ $resume->{lines} }, $line2 );
	}
	return (
		transform_array_to_json_like($resume),
		transform_array_to_json_like($resume2)
	);
}

sub identity_patient {
	my ($p) = @_;
	die();
	my $project = $p->getProject();
	my $line;
	my $hthis;
	push( @$line, { text => $p->name, type => "default" } );
	my $name = $p->name();
	my $h;
	my $nb = 0;
	my $hn;

	foreach my $v ( @{ $p->getVariations() } ) {
		last if $nb > 1000;
		next if $v->getPourcentAllele($p) < 10;
		next if $v->frequency > 0.05;

		my $hres = $project->getDejaVuInfos( $v->id );
		foreach my $p ( keys %$hres ) {

			next if $p eq $project->name;
			my @ps = split( ";", $hres->{$p}->{patients} );
			foreach my $pname (@ps) {
				next if $pname eq $name;
				$h->{ "<td>" . $p . "</td>" . "<td>" . $pname . "</td>" }++;
			}
		}
		$nb++;
		foreach my $p1 ( @{ $v->getPatients } ) {
			next if $p1->name eq $name;
			$hn->{ $p1->name }++;
			$h->{   "<td>this</td>" . "<td>"
				  . $p1->name . " "
				  . $p->getFamilly->name
				  . "</td>" }++;
		}
	}

	#	warn Dumper $hn;
	#	die();
	my $max = -1;
	my $namemax;
	my $hstats;
	my @t    = sort { $h->{$b} <=> $h->{$a} } keys %$h;
	my $text = "<table class='table table-sm table-striped table-bordered '>";
	my $max_pct = 0;
	for ( my $i = 0 ; $i < 3 ; $i++ ) {
		my $id  = $t[$i];
		my $pct = int( ( $h->{$id} / $nb ) * 100 );
		$max_pct = $pct if $pct > $max_pct;
		$text .= $id . " <td>" . $pct . "%</td></tr>";
	}
	$text .= "</table>";
	my $color = "success";
	if ( $max_pct > 90 ) {
		$color = "danger";
	}
	elsif ( $max_pct > 70 ) {
		$color = "warning";
	}

	push( @$line, { text => "$text", type => "$color" } );
	return $line;
}

sub identity {
	my ( $project, $fork ) = @_;
	my $patients = $project->getPatients();
	my $chrs;
	unless ( $project->isDiagnostic ) {
		$chrs = $project->get_only_list_chromosomes("1,5,10,20,21,7,9")
		  if $project->isExome;
		$chrs = $project->get_only_list_chromosomes("1,20,21")
		  if $project->isGenome;

	}

	my $vs = $project->getVariations();

	my $h;
	my $nb;
	my $hn;
	my $limit = 0;
	foreach my $v ( shuffle @$vs ) {
		next unless $v->frequency;
		next if $v->frequency > 0.2;

		#  		warn $limit;
		$limit++;
		last if $limit > 10000;

		#next unless $v->isCoding;
		my $patients = $v->getPatients();
		my @ps;
		my $hres = $project->getDejaVuInfos( $v->id );
		my $z    = 0;
		foreach my $p ( keys %$hres ) {
			foreach my $pname ( split( ";", $hres->{$p}->{patients} ) ) {
				$z++;
			}
		}

		#next if $z > 30;
		my $found;
		foreach my $p ( keys %$hres ) {
			next if $p eq $project->name;

			foreach my $pname ( split( ";", $hres->{$p}->{patients} ) ) {
				push( @ps, $p . " " . $pname );

				#push(@ps, "<td>".$p."</td>"."<td>".$pname."</td>");

			}

		}
		foreach my $p ( @{$patients} ) {
			my $pname = $p->name();
			my $p     = $project->name();
			push( @ps, $p . " " . $pname );

			#push(@ps, "<td>".$p."</td>"."<td>".$pname."</td>");
		}
		foreach my $p ( @{$patients} ) {
			next if $v->getDepth($p) < 10;
			next if $v->getPourcentAllele($p) < 40;

			#  			 warn $v->getPourcentAllele($p);
			my $name = $p->name();
			map { $h->{$name}->{$_}++ } @ps;
			$nb->{$name}++;
		}
	}
	my $results;
	foreach my $p ( @{ $project->getPatients } ) {
		my $pname = $p->name();
		my $line;
		push( @$line, { text => $p->name, type => "default" } );
		my @t1 = sort { $h->{$pname}->{$b} <=> $h->{$pname}->{$a} }
		  keys %{ $h->{$pname} };
		my @t;
		foreach my $a (@t1) {
			my ( $pr, $na ) = split( " ", $a );
			push( @t, $a ) if $na ne $pname;
		}

		#shift(@t);
		my $text =
		  "<table class='table table-sm table-striped table-bordered '>";
		my $max_pct = 0;
		for ( my $i = 0 ; $i < 5 ; $i++ ) {

			my $id = $t[$i];
			warn $pname if $nb->{$pname} == 0;
			$nb->{$pname} = 1 if $nb->{$pname} == 0;
			$nb->{$pname} = 1 unless $nb->{$pname};
			my $pct = int( ( $h->{$pname}->{$id} / $nb->{$pname} ) * 100 );
			$max_pct = $pct if $pct > $max_pct;
			last if $i > 1 && $pct < 10;
			my $class = "";
			my ( $pr, $na ) = split( " ", $id );
			my $sf;
			my $text_f = "";

			if ( $pr eq $p->project->name ) {
				my $p1 = $project->getPatient($na);
				if ( $p1->getFamily->name eq $p->getFamily->name ) {
					$sf     = 1;
					$text_f = " (child" if $p->isChild;
					$text_f = " (mother" if $p->isMother;
					$text_f = " (father" if $p->isFather;
					$text_f .= "-child)"  if $p1->isChild;
					$text_f .= "-mother)" if $p1->isMother;
					$text_f .= "-father)" if $p1->isFather;

				}
			}

			my $aid = "<td>" . $pr . "</td>" . "<td>" . $na . $text_f . "</td>";
			if ($sf) {
				$class = qq{class='bg-info'};
			}
			else {
				if ( $pct >= 85 ) {
					$class = qq{class='bg-danger'};
				}
				elsif ( $pct >= 70 ) {
					$class = qq{class='bg-warning'};
				}
			}
			$text .=
				"<tr $class>"
			  . $aid
			  . " <td >"
			  . $pct . "%"
			  . $h->{$pname}->{$id} . " "
			  . $nb->{$pname}
			  . " </td></tr>";

		}
		$text .= "</table>";
		my $color = "success";
		if ( $max_pct >= 85 ) {
			$color = "danger";
		}
		elsif ( $max_pct > 70 ) {
			$color = "warning";
		}
		push( @$line, { text => "$text", type => "$color" } );
		$results->{$pname} = $line;
	}
	my $resume = {};
	$resume->{header} = [ "patients", "identity" ];
	foreach my $k ( keys %$results ) {

		push( @{ $resume->{lines} }, $results->{$k} );
	}

	#$resume->{lines} = values %$h;
	return ( transform_array_to_json_like($resume) );

	#warn Dumper $h;
	#die();

}

sub identity2 {
	my ( $project, $fork ) = @_;
	my $patients = $project->getPatients();
	unless ( $project->isDiagnostic ) {
		my $chrs = $project->get_only_list_chromosomes("1,20,21");
	}
	warn "coucou";
	my $vs = $project->getVariations();
	$project->buffer->dbh_deconnect();
	my $h;
	my $nb;
	my $hn;
	foreach my $v (@$vs) {
		next if $v->frequency > 0.01;
		my $patients = $v->getPatients();
		my $hres     = $project->getDejaVuInfos( $v->id );
		my @ps;

		foreach my $p ( keys %$hres ) {
			next if $p eq $project->name;

			foreach my $pname ( split( ";", $hres->{$p}->{patients} ) ) {
				push( @ps, "<td>" . $p . "</td>" . "<td>" . $pname . "</td>" );
			}

		}
		foreach my $p ( @{$patients} ) {
			my $pname = $p->name();
			my $p     = $project->name();
			push( @ps, "<td>" . $p . "</td>" . "<td>" . $pname . "</td>" );
		}
		foreach my $p ( @{$patients} ) {

			my $name = $p->name();
			map { $h->{$name}->{$_}++ } @ps;
			$nb->{$name}++;
		}
	}
	warn Dumper $nb;
	my $results;
	foreach my $p ( @{ $project->getPatients } ) {
		my $pname = $p->name();
		my $line;
		push( @$line, { text => $p->name, type => "default" } );
		my @t = sort { $h->{$pname}->{$b} <=> $h->{$pname}->{$a} }
		  keys %{ $h->{$pname} };
		shift(@t);
		my $text =
		  "<table class='table table-sm table-striped table-bordered '>";
		my $max_pct = 0;
		for ( my $i = 0 ; $i < 3 ; $i++ ) {
			my $id = $t[$i];

			#CAS des temoins negatifs
			my $pct;
			if ( int( $nb->{$pname} ) == 0 ) { $pct = 0; }
			else {
				$pct = int( ( $h->{$pname}->{$id} / $nb->{$pname} ) * 100 );
			}
			$max_pct = $pct if $pct > $max_pct;
			$text .= $id . " <td>" . $pct . "%</td></tr>";
		}
		$text .= "</table>";
		my $color = "success";
		if ( $max_pct > 90 ) {
			$color = "danger";
		}
		elsif ( $max_pct > 70 ) {
			$color = "warning";
		}
		push( @$line, { text => "$text", type => "$color" } );
		$results->{$pname} = $line;
	}
	my $resume = {};
	$resume->{header} = [ "patients", "identity" ];
	foreach my $k ( keys %$results ) {

		push( @{ $resume->{lines} }, $results->{$k} );
	}

	#$resume->{lines} = values %$h;
	return ( transform_array_to_json_like($resume) );

	#warn Dumper $h;
	#die();

}

sub coverage_stats {
	my ($project) = @_;
	my $patients = $project->getPatients();
	my $resume;
	my $bam_stats;
	$resume->{title} = "Global Coverage Statistics  : ";
	my $samtools = $project->buffer->software("samtools");
	my $mean_cov;
	foreach my $p (@$patients) {
		my $bam = $p->getBamFile();
		my $cov = $p->coverage();
		$mean_cov += $cov->{mean};
	}
	$mean_cov = $mean_cov / scalar(@$patients);
	my $limit_cov;
	$limit_cov->{danger}  = $mean_cov * 0.5;
	$limit_cov->{warning} = $mean_cov * 0.75;

	my @header = ("patient");
	$resume->{header} = [ "patients", "mean", "15X", "30X", "100X" ];

	# my $it = natatime , @tchromosomes;

	my @lines;
	for ( my $i = 0 ; $i < @$patients ; $i++ ) {
		my @line;
		my $hline;
		push( @$hline, { text => $patients->[$i]->name, type => "default" } );
		my $cov   = $patients->[$i]->coverage();
		my $color = "success";
		$color = "warning" if $cov->{mean} < $limit_cov->{warning};
		$color = "danger"  if $cov->{mean} < $limit_cov->{danger};

		push( @$hline, { text => $cov->{mean}, type => "$color" } );
		$color = "success";
		$color = "warning" if $cov->{"15x"} < 95;
		$color = "danger" if $cov->{"15x"} < 85;

		push( @$hline, { text => $cov->{"15x"}, type => "$color" } );
		$color = "success";
		$color = "warning" if $cov->{"30x"} < 92;
		$color = "danger" if $cov->{"30x"} < 80;

		push( @$hline, { text => $cov->{"30x"}, type => "$color" } );
		$color = "success";
		$color = "warning" if $cov->{"30x"} < 90;
		$color = "danger" if $cov->{"30x"} < 75;

		push( @$hline, { text => $cov->{"100x"}, type => "$color" } );
		push( @{ $resume->{lines} }, $hline );
	}

	return ( transform_array_to_json_like($resume) );

}

sub bam_stats {
	my ($project) = @_;
	my $patients = $project->getPatients();
	my $resume;
	my $bam_stats;
	$resume->{title} = "Global Coverage Statistics  : ";
	my $samtools = $project->buffer->software("samtools");
	my $mean_cov;
	foreach my $p (@$patients) {
		my $bam = $p->getBamFile();
		my $cov = $p->coverage();
		$mean_cov += $cov->{mean};
		my $fsize = -s $p->getBamFile();
		my $cmd   = qq{$samtools idxstats $bam};
		my @t     = `$cmd`;
		chomp(@t);
		foreach my $l (@t) {

			my ( $chr, $size, $amount ) = split( " ", $l );
			next if $chr eq "*";
			next unless $project->isChromosomeName($chr);

			#next if uc($chr) =~ /GL/;
			#next if uc($chr) =~ /HS/;
			$amount *= 2 if ( $chr eq "chrX" && $p->isMale() );

			#warn "coucou"  if ($chr eq "chrX" && $p->isMale());

			push( @{ $bam_stats->{$chr}->{data} }, $amount );

		}
	}
	$mean_cov = $mean_cov / scalar(@$patients);
	my $limit_cov;
	$limit_cov->{danger}  = $mean_cov * 0.5;
	$limit_cov->{warning} = $mean_cov * 0.75;

	my @header = ("patient");

	my @tchromosomes = sort {
		$project->getChromosome($a)
		  ->karyotypeId <=> $project->getChromosome($b)->karyotypeId
	} keys %$bam_stats;
	$resume->{header} =
	  [ "patients", "mean", "15X", "30X", "100X", @tchromosomes ];

	# my $it = natatime , @tchromosomes;

	foreach my $chr (@tchromosomes) {

		next if $chr eq "*";
		push( @header, $chr );
		my $t    = sum @{ $bam_stats->{$chr}->{data} };
		my $mean = int( $t / scalar @{ $bam_stats->{$chr}->{data} } );
		$bam_stats->{$chr}->{mean} = $mean;
	}

	my @lines;
	for ( my $i = 0 ; $i < @$patients ; $i++ ) {
		my @line;
		my $hline;
		push( @$hline, { text => $patients->[$i]->name, type => "default" } );
		my $cov   = $patients->[$i]->coverage();
		my $color = "success";
		$color = "danger"  if $cov->{mean} < $limit_cov->{danger};
		$color = "warning" if $cov->{mean} < $limit_cov->{warning};
		push( @$hline, { text => $cov->{mean}, type => "$color" } );
		$color = "success";
		$color = "danger" if $cov->{"15x"} < 85;
		$color = "warning" if $cov->{"15x"} < 95;
		push( @$hline, { text => $cov->{"15x"}, type => "$color" } );
		$color = "success";
		$color = "danger" if $cov->{"30x"} < 80;
		$color = "warning" if $cov->{"30x"} < 92;
		push( @$hline, { text => $cov->{"30x"}, type => "$color" } );
		$color = "success";
		$color = "danger" if $cov->{"30x"} < 75;
		$color = "warning" if $cov->{"30x"} < 90;
		push( @$hline, { text => $cov->{"100x"}, type => "$color" } );

		foreach my $chr (@tchromosomes) {

			#push(@line,$bam_stats->{$chr}->{mean});

			my $v = $bam_stats->{$chr}->{data}->[$i];
			my $m = $bam_stats->{$chr}->{mean};
			if ( $chr eq "chrY" && $patients->[$i]->isFemale ) {
				push( @$hline, { text => "-", type => "info" } );
				next;
			}

			if ( $v eq 0 ) {
				push( @$hline, { text => "NONE", type => "danger" } );
			}
			elsif ( $v < abs( $m * 0.3 ) ) {
				push( @$hline, { text => "POORLY", type => "warning" } );
			}
			elsif ( $v < abs( $m * 0.5 ) ) {
				push( @$hline, { text => "LOW", type => "info" } );
			}
			else {
				push( @$hline, { text => "OK", type => "success" } );
			}

		}
		push( @{ $resume->{lines} }, $hline );
	}

	return ( transform_array_to_json_like($resume) );

}

sub fork_patients {
	my ( $project, $sub, $fork, $chr ) = @_;

	my $pm = new Parallel::ForkManager($fork);
	my $results;
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $hres ) =
			  @_;
			unless ( defined($hres) or $exit_code > 0 ) {
				print
				  qq|No message received from child process $exit_code $pid!\n|;
				warn Dumper $hres;
				return;
			}
			foreach my $p ( keys %$hres ) {
				$results->{$p} = $hres->{$p};
			}

		}
	);
	my @patients_name = map { $_->name } @{ $project->getPatients };
	$project->buffer->dbh_deconnect();
	foreach my $p (@patients_name) {
		my $pid;
		$pid = $pm->start and next;
		$project->buffer->dbh_reconnect();
		my $patient = $project->getPatient($p);

		my $tab = $sub->($patient);
		my $hres;
		$hres->{$p} = $tab;
		$pm->finish( 0, $hres );
	}    #end for range range
	$pm->wait_all_children();

	$project->buffer->dbh_reconnect();
	return $results;
}

sub coverage_transcripts_patient {
	my ($patient)       = @_;
	my $project         = $patient->getProject();
	my @transcripts_cgi = @{ $project->bundle_transcripts() };
	my $transcripts;

	foreach my $t (@transcripts_cgi) {
		push( @$transcripts, $project->newTranscript($t) );
	}
	my $cols = [];
	push( @$cols, { text => $patient->name, type => "default" } );
	my $tsum = 0;
	foreach my $t (@$transcripts) {
		my $toto = $t->getGene()->get_coverage($patient)
		  ->coverage_intspan( $t->getSpanCoding );
		my $n    = int( $toto->{mean} );
		my $min  = $toto->{min};
		my $type = "success";
		$tsum += $n;
		if ( $patient->isFemale() and $t->getChromosome()->name eq 'Y' ) {
			push( @$cols, { text => "-", type => "default" } );
		}
		my $problem;
		if ( $n < 10 ) {
			$type    = "danger";
			$problem = 1;
		}
		elsif ( $n < 30 ) {
			$type    = "warning";
			$problem = 1;
		}
		elsif ( $n < 50 ) {

			$type = "info";
		}
		my $text = $n . " ($min)";
		push( @$cols, { text => $text, type => $type } );
	}
	my $n       = int( $tsum / scalar(@$transcripts) );
	my $problem = 1;
	my $type    = "success";
	if ( $n < 10 ) {
		$type = "danger";
	}
	elsif ( $n < 30 ) {
		$type = "warning";
	}
	elsif ( $n < 50 ) {
		$problem = undef;
		$type    = "info";
	}

	push( @$cols, { text => $n, type => $type } );
	return $cols;
}

sub coverage_transcripts {
	my ( $project, $fork ) = @_;
	my $resume;
	return unless $project->bundle_transcripts();
	my @transcripts_cgi = @{ $project->bundle_transcripts() };
	$resume->{header} = [ "patients", @transcripts_cgi, "all" ];
	my $h = fork_patients( $project, \&coverage_transcripts_patient, $fork );
	foreach my $k ( keys %$h ) {
		push( @{ $resume->{lines} }, $h->{$k} );
	}

	#$resume->{lines} = values %$h;
	return ( transform_array_to_json_like($resume) );

}

sub coverage_transcripts2 {
	my ($project) = @_;
	my $patients = $project->getPatients();
	my $resume2;

	my $resume;
	my @transcripts_cgi = @{ $project->bundle_transcripts() };
	$resume2->{header} = ["patients"];
	my $index = 0;
	foreach my $t (@transcripts_cgi) {
		$index++;
		push( @{ $resume2->{header} }, $index );
	}
	$resume->{header} = [ "patients", @transcripts_cgi, "all" ];

	$resume->{title} = " Coverage Statistics by Transcripts : ";
	my $transcripts;

	foreach my $t (@transcripts_cgi) {
		push( @$transcripts, $project->newTranscript($t) );
	}
	my $nb;
	my $tr_problem;
	foreach my $p (@$patients) {
		warn $p->name;
		$nb++;
		my $cols;
		my $cols2;
		push( @$cols, { text => $p->name, type => "default" } );
		my $tsum  = 0;
		my $index = 0;
		foreach my $t (@$transcripts) {
			my $toto = $t->getGene()->get_coverage($p)
			  ->coverage_intspan( $t->getSpanCoding );
			my $n    = int( $toto->{mean} );
			my $min  = $toto->{min};
			my $type = "success";
			$tsum += $n;
			if ( $p->isFemale() and $t->getChromosome()->name eq 'Y' ) {
				push( @$cols,  { text => "-", type => "default" } );
				push( @$cols2, { text => "-", type => "default" } );
			}
			my $problem;
			if ( $n < 10 ) {
				$type    = "danger";
				$problem = 1;
			}
			elsif ( $n < 30 ) {
				$type    = "warning";
				$problem = 1;
			}
			elsif ( $n < 50 ) {

				$type = "info";
			}
			my $text = $n . " ($min)";
			$tr_problem->{$index}++;    # if $problem;
			push( @$cols,  { text => $text,    type => $type } );
			push( @$cols2, { text => "&nbsp.", type => $type } );
			$index++;
		}
		my $n       = int( $tsum / scalar(@$transcripts) );
		my $problem = 1;
		my $type    = "success";
		if ( $n < 10 ) {
			$type = "danger";
		}
		elsif ( $n < 30 ) {
			$type = "warning";
		}
		elsif ( $n < 50 ) {
			$problem = undef;
			$type    = "info";
		}

		push( @$cols, { text => $n, type => $type } );
		push( @{ $resume2->{lines} }, $cols2 );
		push( @{ $resume->{lines} },  $cols );
	}
	return ( transform_array_to_json_like($resume) );
	my $resume3;
	$resume3->{header} = ["patients"];
	foreach my $index ( keys %{$tr_problem} ) {
		push( @{ $resume3->{header} }, $transcripts->[$index]->name );
	}

	foreach my $l ( @{ $resume->{lines} } ) {
		my $l1;
		push( @{$l1}, $l->[0] );
		foreach my $index ( keys %{$tr_problem} ) {

			push( @{$l1}, $l->[ $index + 1 ] );
		}
		push( @{ $resume3->{lines} }, $l1 );
	}

	$resume->{legend} = [
		{ text => "<10",   type => "danger" },
		{ text => "<10",   type => "danger" },
		{ text => "<30",   type => "warning" },
		{ text => "<50",   type => "info" },
		{ text => ">= 50", type => "success" }
	];
	return ( transform_array_to_json_like($resume3) );

}

sub duplicate_regions {
	my ($project) = @_;

	my $patients = $project->getPatients();
	my $resume;
	$resume->{title} = "Looking For Duplicate Regions ";

	my @bams;

	$resume->{header} = [ "Genes", "transcript", "exon", "% dup", "position" ];
	my $hash_dup;
	my $htr;
	my $he;
	my $ws;
	return unless $project->bundle_transcripts();
	map { $htr->{$_}++ } @{ $project->bundle_transcripts() };
	foreach my $patient ( @{ $project->getPatients } ) {

		my $filebed =
			$patient->project->getVariationsDir("duplicate_region_calling")
		  . "/regions/"
		  . $patient->name()
		  . ".dup.bed";

		#system("mkdir $dir && chmod a+rwx $dir" ) unless -e $dir;
		return unless -e $filebed;
		if ( -e $filebed ) {
			open( BED, $filebed );
			while (<BED>) {
				chomp();
				my ( $chr, $start, $end ) = split(" ");
				unless ( exists $hash_dup->{$chr} ) {
					$hash_dup->{$chr} = Set::IntSpan::Fast::XS->new();
				}

				$hash_dup->{$chr}->add_range( $start, $end );
			}
		}
	}

	foreach my $cn ( keys %$hash_dup ) {
		my $chr  = $project->getChromosome($cn);
		my $iter = $hash_dup->{$cn}->iterate_runs();
		while ( my ( $from, $to ) = $iter->() ) {

			my $ts = $chr->getTranscriptsByPosition( $from, $to );
			foreach my $t (@$ts) {
				next unless exists $htr->{ $t->name };
				foreach my $e ( @{ $t->getExons } ) {
					next if exists $he->{ $e->name };

					my $span =
					  $hash_dup->{$cn}->intersection( $e->getGenomicSpan );
					my $l1 = scalar( $e->getGenomicSpan->as_array );
					my $l2 = scalar( $span->as_array );
					my $p  = int( ( $l2 / $l1 * 100 ) );
					unless ( $span->is_empty ) {
						my $i1 = $span->iterate_runs();
						my @st;
						while ( my ( $from, $to ) = $i1->() ) {
							push( @st,
									"["
								  . abs( $from - $e->start ) . "-"
								  . abs( $to - $e->start )
								  . "]" );
						}
						$he->{ $e->name }++;
						$ws->{ $t->getGene->external_name() }->{tr}
						  ->{ $t->name }->{ $e->name }->{pos} =
						  join( " ", @st );
						$ws->{ $t->getGene->external_name() }->{tr}
						  ->{ $t->name }->{ $e->name }->{percent} = $p . "%";
						$ws->{ $t->getGene->external_name() }->{row}++;

					}

				}
			}
		}
	}
	return $resume unless ( keys %{$ws} );
	foreach my $g ( keys %{$ws} ) {

		my $nbc = scalar( keys %{ $ws->{$g}->{tr} } );

		foreach my $t ( keys %{ $ws->{$g}->{tr} } ) {
			my $nbe = scalar( keys %{ $ws->{$g}->{tr}->{$t} } );
			foreach my $e ( keys %{ $ws->{$g}->{tr}->{$t} } ) {
				my $hline;
				push( @$hline,
					{ text => "$g", type => "default", rowspan => $nbe } );
				push( @$hline,
					{ text => "$t", type => "default", rowspan => $nbe } );
				push( @$hline, { text => "$e", type => "default" } );
				push(
					@$hline,
					{
						text => $ws->{$g}->{tr}->{$t}->{$e}->{percent},
						type => "danger"
					}
				);
				push(
					@$hline,
					{
						text => $ws->{$g}->{tr}->{$t}->{$e}->{pos},
						type => "danger"
					}
				);
				push( @{ $resume->{lines} }, $hline );

			}
		}
	}
	return ( transform_array_to_json_like($resume) );

}

sub return_mendelian_color {
	my ($value) = @_;
	if ( $value < 1 ) {
		return ( "success", "OK" );
	}
	elsif ( $value < 2 ) {
		return ( "warning", "WARNING" );
	}
	else {
		return ( "danger", "ERROR" );
	}

}

sub return_sex_color {
	my ( $sex1, $sex2, $sex3, $score ) = @_;

	my $color = "danger";
	$color = "success" if $sex1 eq $sex2;
	my $text = "-";
	if ( $sex2 eq "2" ) {
		$text = qq{<i class="fa fa-venus "></i>};
	}
	elsif ( $sex2 eq "1" ) {
		$text = qq{<i class="fa fa-mars  "></i>};
	}
	else {
		$text  = qq{<i class="fa fa-question  "></i>};
		$color = "warning";
	}
	if ( $sex3 && $sex3 eq $sex1 && $sex2 ne $sex1 ) {
		$color = "warning";
	}
	if ( $sex3 && $sex3 eq $sex1 && $sex2 == 0 ) {
		$color = "default";
	}

	#  	return ($color,$text) unless $sex3;
	#
	#  	if ($sex3 eq $sex1 && $sex1 ne $sex2) {
	#  		$color= "warning";
	#  	}

	return ( $color, $text );

}

sub table_json_mendelian {
	my ( $project, $h, $fam ) = @_;
	my $color = {
		"danger"  => "danger",
		"warning" => "orange",
		"info"    => "cyan",
		"success" => "green",
	};
	my $resume = {};
	add_columns(
		$resume,
		[
			"results", "familly",   'sample', "sex",
			"SRY",     "plink_sex", "mendelian"
		]
	) if $project->isFamilial;
	add_columns( $resume,
		[ "results", "familly", 'sample', "sex", "SRY", "plink_sex" ] )
	  unless $project->isFamilial;

	#	foreach my $fam (@{$project->getFamilies}){
	my $fname  = $fam->name();
	my $ps     = $fam->getPatients();
	my @values = map { $h->{$fname}->{$_}->{mendelian_errors}->{percent} }
	  keys %{ $h->{$fname} };
	my $max = max(@values);
	my ( $color, $text ) = return_mendelian_color($max);
	my $hfirst_cols;
	my $hcol1;

	if ( $text eq "OK" ) {
		$hcol1 = {
			type    => "default",
			rowspan => scalar(@$ps),
			text =>
			  qq{<i class="fa fa-check-circle fa-2x " style="color:green"></i>}
		};
	}
	elsif ( $text eq "WARNING" ) {
		$hcol1 = {
			type    => "warning",
			rowspan => scalar(@$ps),
			text    => qq{<i class="fa fa-exclamation fa-2x " ></i>}
		};
	}
	else {
		$hcol1 = {
			type    => "danger",
			rowspan => scalar(@$ps),
			text    => qq{<i class="fa fa-exclamation-triangle fa-2x" ></i>}
		};
	}
	$hcol1->{ index_line => $fname };
	$hfirst_cols->{results} = $hcol1;
	$hfirst_cols->{familly} = {
		type       => "default",
		rowspan    => scalar(@$ps),
		text       => $fname,
		index_line => $fname
	};

	#print $cgi->td({class=>"$color",rowspan=>scalar(@$ps)},$fname);
	#

	my $children = $fam->getChildren();
	my $mother   = $fam->getMother;
	my $father   = $fam->getFather;
	my $pc       = $h->{$fname};          #mendelian_errors}->{percent};

	if ($mother) {
		my $hline = dclone $hfirst_cols;

		my $value =
		  $h->{$fname}->{ $mother->name }->{mendelian_errors}->{percent};
		my $value2 =
		  $h->{$fname}->{ $mother->name }->{mendelian_errors}->{nb_snp};
		my ( $color, $text ) = return_mendelian_color($value);
		$hline->{sample} = { type => "default", text => $mother->name };
		$hline->{sex}    = {
			type => "default",
			text =>
qq{<i class="fa fa-female fa-2x" aria-hidden="true" style="color:pink"></i>}
		};
		my $name1 = $mother->name;
		my ( $color_s2, $text2 ) =
		  return_sex_color( $mother->sex,, $mother->compute_sex );
		$hline->{SRY} = { type => "$color_s2", text => $text2 };

		#print $cgi->td({class=>"$color_s2"},[$text2]);
		my ( $color_s1, $text1 ) =
		  return_sex_color( $mother->sex,
			$h->{$fname}->{$name1}->{check_sex}->{plink}->{sex},
			$mother->compute_sex );
		$hline->{plink_sex} = {
			type => "$color_s1",
			text => $text1 . "["
			  . $h->{$fname}->{$name1}->{check_sex}->{plink}->{score} . "]"
		};
		$hline->{mendelian} =
		  { type => "$color", text => "$text [$value% $value2]" };
		push( @{ $resume->{data} }, $hline );
	}
	if ($father) {
		my $hline = dclone $hfirst_cols;

		my $value =
		  $h->{$fname}->{ $father->name }->{mendelian_errors}->{percent};
		my $value2 =
		  $h->{$fname}->{ $father->name }->{mendelian_errors}->{nb_snp};

		my ( $color, $text ) = return_mendelian_color($value);
		$hline->{sample} = { type => "default", text => $father->name };
		$hline->{sex}    = {
			type => "default",
			text =>
qq{<i class="fa fa-male fa-2x" aria-hidden="true" style="color:blue"></i>}
		};
		my ( $color_s2, $text2 ) =
		  return_sex_color( $father->sex, $father->compute_sex );
		$hline->{SRY} = { type => "$color_s2", text => $text2 };
		my $name1 = $father->name;
		my ( $color_s1, $text1 ) =
		  return_sex_color( $father->sex,
			$h->{$fname}->{$name1}->{check_sex}->{plink}->{sex},
			$father->compute_sex );
		$hline->{plink_sex} = {
			type => "$color_s1",
			text => $text1 . "["
			  . $h->{$fname}->{$name1}->{check_sex}->{plink}->{score} . "]"
		};
		$hline->{mendelian} =
		  { type => "$color", text => "$text [$value% $value2]" };
		push( @{ $resume->{data} }, $hline );
	}
	foreach my $p (@$children) {
		my $hline  = dclone $hfirst_cols;
		my $value  = $h->{$fname}->{ $p->name }->{mendelian_errors}->{percent};
		my $value2 = $h->{$fname}->{ $p->name }->{mendelian_errors}->{nb_snp};
		my ( $color, $text ) = return_mendelian_color($value);
		$hline->{sample} = { type => "default", text => $p->name };
		$hline->{sex}    = {
			type => "default",
			text =>
qq{<i class="fa fa-child fa-2x" aria-hidden="true" style="color:black"></i> }
		};
		my $name1 = $p->name;
		my ( $color_s2, $text2 ) = return_sex_color( $p->sex, $p->compute_sex );
		$hline->{SRY} = { type => "$color_s2", text => $text2 };
		my ( $color_s1, $text1 ) =
		  return_sex_color( $p->sex,
			$h->{$fname}->{$name1}->{check_sex}->{plink}->{sex},
			$p->compute_sex );
		$hline->{plink_sex} = {
			type => "$color_s1",
			text => $text1 . "["
			  . $h->{$fname}->{$name1}->{check_sex}->{plink}->{score} . "]"
		};
		$hline->{mendelian} =
		  { type => "$color", text => "$text [$value% $value2]" };

		#return $hline;
		push( @{ $resume->{data} }, $hline );
	}

	#}
	#return $hline;
	return $resume->{data};
}

1;
