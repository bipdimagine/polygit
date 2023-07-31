package parse_gff;
use strict;
use Data::Dumper;

sub read_gff_genes {
	my ($file,$translate_id,$gene_code_version,$id) = @_;
	#open (GFF ,qq{zcat $file  | grep PAR_Y  |  perl -lane 'print \$_ if \$F[2] eq "gene"' | });
	warn $file;
	open (GFF ,qq{zcat $file   |  perl -lane 'print \$_ if \$F[2] eq "gene"' | });
	my $genes;
	while (my $line = <GFF>){
	next if $line =~/^#/;

	chomp($line);
	my @line = split("\t",$line);
	my $type = $line[2];
	my $chr = get_chr($line[0]);
	my $start = $line[3];
	my $end = $line[4];
	my $strand =get_strand($line[6]);
	my $infos = parse_infos($line[8]);
	if ($id) {
		#next if $id ne getId($infos->{ID});
	}
	my $genbo_id =  getId($infos->{ID})."_".$chr;

	push(@{$translate_id->{getId($infos->{ID})}},$genbo_id);
		
    push(@{$translate_id->{$infos->{ID}}},$genbo_id);
		
	my $gene = getId($infos->{gene_id});
	my $genboid_gene = $gene."_".$chr;
	
	#next if $infos->{gene_name} ne "$ngene" and $ngene ne "all" ;
	my $parent = getId($infos->{Parent});
	#next if $chr ne "Y";
	my $genboid_parent = $parent."_".$chr;

		my $id = getId($infos->{ID})."_".$chr;
			$genes->{$genbo_id}->{id} = $id;
		$genes->{$genbo_id}->{strand} = $strand;
		$genes->{$genbo_id}->{external_name} = $infos->{gene_name};
		$genes->{$genbo_id}->{end} = $end;
		$genes->{$genbo_id}->{start} = $start;
		$genes->{$genbo_id}->{transcripts} = [];
		$genes->{$genbo_id}->{name} =  getId($infos->{ID});
		$genes->{$genbo_id}->{genbo_id} =  $id;
		$genes->{$genbo_id}->{stable_id} =  $infos->{ID};
		$genes->{$genbo_id}->{chromosome} =  $chr;
		 $genes->{$genbo_id}->{gene_code}  =$gene_code_version;
		 $genes->{$genbo_id}->{genbo_type} = "gene";
		 #$genes->{$genbo_id}->{chr_name} = $chr;
		  $genes->{$genbo_id}->{chromosome} = $chr;
		  $genes->{$genbo_id}->{description} =   $infos->{description};
		   $genes->{$genbo_id}->{description} =~s/-/ /g if $genes->{$genbo_id}->{description};
		  #objects 
		  $genes->{$genbo_id}->{remap_status} = $infos->{remap_status};
		   $genes->{$genbo_id}->{chromosomes_object}->{$chr} = undef;
		  
		 
}

return $genes;
}


sub read_gff_transcripts {
	my ($file,$genes,$translate_id,$gene_code_version,$gid) = @_;
	my $proteins;
	my $transcripts;
	my $unprocess;
	open (GFF ,qq{zcat $file  |  perl -lane 'print \$_ if \$F[2] ne "gene"' | });
	#open (GFF ,qq{zcat $file | grep PAR_Y |  perl -lane 'print \$_ if \$F[2] ne "gene"' | });
	while (my $line = <GFF>){
	next if $line =~/^#/;
	chomp($line);
	my @line = split(" ",$line);
	my $type = $line[2];
	my $chr = get_chr($line[0]);
	my $start = $line[3];
	my $end = $line[4];
	my $strand =get_strand($line[6]);
	my $infos = parse_infos($line[8]);
	my $genbo_id =  getId($infos->{ID})."_".$chr;

	my $parent = getId($infos->{Parent});
	my $genboid_parent = $parent."_".$chr;
	my $gene_id = getId($infos->{gene_id});
	next if exists $unprocess->{$genboid_parent};
	if($gid){
    	#	next if $gid ne $gene_id;
    	}
	if ($type eq "transcript") {
		
    	push(@{$translate_id->{$infos->{ID}}},$genbo_id);
   
 		

    	
    #	warn $genboid_parent;
    
		die($line) unless exists $genes->{$genboid_parent};
		my $hgene = $genes->{$genboid_parent};
		my $new_transcript;
		#$VAR2 = 'gene_description';
		$new_transcript->{gene_description} = $genes->{$genboid_parent}->{description};
		#$VAR2 = 'name';
		my $tl = $infos->{transcript_support_level};
		$tl = 6 if $tl eq "NA";
		$tl = 6 if $tl eq "Missing";
		$new_transcript->{transcript_support_level} =  $tl;
		foreach my $t (split(",",$infos->{tag})){
	 		
	 		$new_transcript->{tag}->{$t} ++;
	 		$t = uc($t);
	 		if (uc($t) =~/PRINCIPAL/){
	 			$new_transcript->{tag}->{PRINCIPAL} ++;
	 		}
	 		if ($t =~/APPRIS/ && $t !~/CANDIDATE/) {
	 			$new_transcript->{tag}->{APPRIS} ++;
	 		}
	 		if ($t =~/MANE/) {
	 			$new_transcript->{tag}->{MANE} ++;
	 		}
	 		if ($t =~/_NF/ ) {
	 			$new_transcript->{tag}->{NF} ++;
	 		}
	 	}
		$new_transcript->{name} =  getId($infos->{ID});
		$new_transcript->{stable_id} = $infos->{ID};
		#$VAR3 = 'ccds_name';
		$new_transcript->{ccds_name} = getId($infos->{ccdsid});
		#$VAR4 = 'chromosome';
		$new_transcript->{chromosome} = $chr;
		#$VAR7 = 'start';
		 $new_transcript->{start} = $start;
		 #$VAR6 = 'end';
	 	$new_transcript->{end} = $end;
	 	#$VAR10 = 'strand';
	 	$new_transcript->{strand} = $strand;
	 	#$VAR8 = 'gene_external_name';
	 	$new_transcript->{gene_external_name}  = $infos->{gene_name};
	    $new_transcript->{gene_code}  =$gene_code_version;
		#$VAR9 = 'genbo_id';
		$new_transcript->{id} =  $genbo_id;
		$new_transcript->{genbo_id} = getId($infos->{ID})."_".$chr;
		$new_transcript->{chr_name} = $chr;
		$new_transcript->{chromosome} = $chr;
		#$VAR12 = 'gene_kyoto_id';
		$new_transcript->{gene_kyoto_id} = getId($infos->{gene_id})."_".$chr;
		#$new_transcript->{gene} = getId($infos->{gene_id})."_".$chr;
		#$VAR15 = 'gene';
		$new_transcript->{gene} = $hgene->{genbo_id};#getId($infos->{gene_id});
		$new_transcript->{gene} = $hgene->{genbo_id};
		
	 	#$VAR19 = 'biotype_ensembl';
	 	
	 	
	 	
	 	
	 	$new_transcript->{biotype_ensembl} =$infos->{transcript_type};
	 	$new_transcript->{transcript_name} =$infos->{transcript_name};
	 	#$new_transcript->{gencode_tag} =$infos->{tag};
		#$VAR23 = 'biotype';
		$new_transcript->{biotype} =$infos->{transcript_type};
		$new_transcript->{genbo_type} = "transcript";
		#intspan
		init_span($new_transcript);
		$new_transcript->{span_transcript}->add_range($start,$end);
		
		#_objects
		$new_transcript->{genes_object} ->{$hgene->{genbo_id}} = undef; #getId($infos->{gene_id});
		$new_transcript->{chromosomes_object} ->{$hgene->{chromosome}} = undef; #getId($infos->{gene_id});
		$hgene->{transcripts_object}->{$new_transcript->{id}} = undef;
		
		
		#PROTEIN 
		if ( $infos->{protein_id}){
		#$VAR17 = 'protein';
		
		$new_transcript->{protein_id} = getId($infos->{protein_id})."_".$chr;
		$new_transcript->{protein} = $new_transcript->{protein_id};
		$new_transcript->{proteins_object}->{$new_transcript->{protein}} = undef;
		my $prot_id = $new_transcript->{protein_id} ;
		die( $prot_id) if $prot_id =~/ENST/;
 		$new_transcript->{protein_stable_id} = $infos->{protein_id};
		$new_transcript->{protein_genbo_id} = $infos->{protein_id};
		
	 	$proteins->{$prot_id}->{id} = $prot_id;
	 	$proteins->{$prot_id}->{genbo_id} = $prot_id;
	 	$proteins->{$prot_id}->{name} =  getId($infos->{protein_id});
	 	$proteins->{$prot_id}->{transcript} =$new_transcript->{genbo_id};
	 	$proteins->{$prot_id}->{gene} = $new_transcript->{gene};
	 	$proteins->{$prot_id}->{gene_kyoto_id} = $new_transcript->{gene};
	 	
	 	
	 	$proteins->{$prot_id}->{strand} = $new_transcript->{strand};
	 	$proteins->{$prot_id}->{start} = 1;
	 	$proteins->{$prot_id}->{genbo_type} = "protein";
	 	$proteins->{$prot_id}->{chr_name} = $chr;
	 	$proteins->{$prot_id}->{chromosome} = $chr;
	 	#objects 
	 	$new_transcript->{proteins_object}->{$new_transcript->{protein}} = undef;
	 	$proteins->{$prot_id}->{transcripts_object}->{$new_transcript->{genbo_id}} = undef;
	 	$proteins->{$prot_id}->{genes_object} ->{$new_transcript->{gene}} = undef;
	 	$proteins->{$prot_id}->{chromosomes_object} ->{$new_transcript->{chromosome}} = undef;
	 	
	 	push(@{$translate_id->{$infos->{protein_id}}} ,$prot_id);
	 	
		}

		$transcripts->{$new_transcript->{genbo_id}} = $new_transcript;
		push(@{$hgene->{transcripts}},$new_transcript->{genbo_id});

	} #end if transcripts
	
	if ($type eq "stop_codon" or $type eq "start_codon"){
		
		die() unless exists  $transcripts->{$genboid_parent};
		my $prot_id = $transcripts->{$genboid_parent}->{protein_id} ;
		die() unless $prot_id;
		if ($type eq "start_codon"){
		 $transcripts->{$genboid_parent}->{start_codon}->{start} = $start;
		 $transcripts->{$genboid_parent}->{start_codon}->{end} = $end; 
		 $proteins->{$prot_id}->{start} = $start;
		 
		}
		 else {
		 		$transcripts->{$genboid_parent}->{end_codon}->{start} = $start;
		 		$transcripts->{$genboid_parent}->{end_codon}->{end} = $end; 
		 		$proteins->{$prot_id}->{end} = $end;
		 }
		 $transcripts->{$genboid_parent}->{spanCodonStartEnd}->add_range($start,$end);
		 #warn $transcripts->{$genboid_parent}->{spanCodonStartEnd}->as_string();
	} #end start stop
	
	if ($type eq "five_prime_UTR") {
			die() unless exists  $transcripts->{$genboid_parent};
		#warn "\t".$type." ".$genboid_parent." $start $end";
		$transcripts->{$genboid_parent}->{five_prime_UTR}->add_range($start,$end);
		
		#warn Dumper $infos;
		
	}
		if ($type eq "three_prime_UTR") {
			die() unless exists  $transcripts->{$genboid_parent};
		#warn "\t".$type." ".$genboid_parent." $start $end";
		$transcripts->{$genboid_parent}->{three_prime_UTR}->add_range($start,$end);
		
		#warn Dumper $infos;
		
	}
	if ($type eq "CDS") {
			die() unless exists  $transcripts->{$genboid_parent};
			my $prot_id = $transcripts->{$genboid_parent}->{protein_id} ;
			die($genboid_parent." ".$line) unless $prot_id;
		#warn "\t".$type." ".$genboid_parent." $start $end";
		$transcripts->{$genboid_parent}->{cds} = 1;
		$transcripts->{$genboid_parent}->{span_coding}->add_range($start,$end);
		$proteins->{$prot_id}->{genomic_span} = $transcripts->{$genboid_parent}->{span_coding}->copy;
		
		
		#warn Dumper $infos;
		
	}
		if ($type eq "exon") {
			
			die($genboid_parent) unless exists  $transcripts->{$genboid_parent};
			 $transcripts->{$genboid_parent}->{genomic_span}->add_range($start,$end);
			 my $prot_id = $transcripts->{$genboid_parent}->{protein_id} ;
			 $proteins->{$prot_id}->{transcript_genomic_span} = $transcripts->{$genboid_parent}->{genomic_span}->copy if $prot_id;
			 
			$transcripts->{$genboid_parent}->{essential_splice_site_span}->add_range($start-2,$start-1);
			$transcripts->{$genboid_parent}->{essential_splice_site_span}->add_range($end+1,$end+2);
			$transcripts->{$genboid_parent}->{splice_site_span}->add_range($start-8,$start-2);
			$transcripts->{$genboid_parent}->{splice_site_span}->add_range($start,$start+2);
			$transcripts->{$genboid_parent}->{splice_site_span}->add_range($end+2,$end+8);
			$transcripts->{$genboid_parent}->{splice_site_span}->add_range($end-2,$end);
		
	}
	
	
	}
	
	foreach my $prot (values %$proteins){
		my $tr_id = $prot->{transcript};
		my $transcript  = $transcripts->{$tr_id};
		die() unless $transcript;
		$prot->{genomic_span} = $transcript->{span_coding}->copy;
		my $cstart =	 scalar($transcript->{five_prime_UTR}->as_array) +1;
		my $cend =  scalar($transcript->{genomic_span}->as_array)- scalar($transcript->{three_prime_UTR}->as_array) ;
		$prot->{cdna_start} = $cstart;
		$prot->{cdna_end} = $cend;
		#$VAR18 = 'orf_start';relatif au transcript
		$transcript->{orf_start} = $cstart;
		$transcript->{orf_end} = $cend;
		#$VAR13 = 'external_name';
		
	}
	
#	warn 'FILTER transcript !!!!!! ';
#	foreach my $tid (keys %$transcripts) {
#		
#		my $new_transcript = $transcripts->{$tid};
#		my $tl = $new_transcript->{transcript_support_level};
#		next if $tl >= 1 ;
#		next if exists  $new_transcript->{tag}->{basic};
#		
#		next if exists $new_transcript->{tag}->{CCDS};
#		
#		next if exists  $new_transcript->{tag}->{APPRIS};
#		if ($tl >=4 or $new_transcript->{biotype_ensembl} eq 'processed_transcript' or exists $new_transcript->{tag}->{NF})  {
#		my $protein_id = $new_transcript->{protein_id};
#		delete $proteins->{$protein_id};
#		my @gene_ids = keys %{$new_transcript->{genes_object}};
##		warn $gene_ids[0];
#		die("ici") if scalar(@gene_ids)>1;
#		die("or here") unless exists $genes->{$gene_ids[0]}->{transcripts_object}->{$new_transcript->{id}};
#		delete $genes->{$gene_ids[0]}->{transcripts_object}->{$new_transcript->{id}};
#		my @ts = keys %{$genes->{$gene_ids[0]}->{transcripts_object}};
#		delete $genes->{$gene_ids[0]} if scalar(@ts) == 0;
#		warn $gene_ids[0] if scalar(@ts) == 0;
#		next ;
#		}
#
#	} 
	
	
	return ($transcripts,$proteins);
}

sub init_span {
	my ($new_transcript,$start,$end) = @_;
	my @types = ("span_mature","span_transcript","spanCodonStartEnd","genomic_span","span_coding","splice_site_span","essential_splice_site_span","five_prime_UTR","three_prime_UTR");
	foreach my $type (@types){
		$new_transcript->{$type} = Set::IntSpan::Fast::XS->new();
	}
	
		
}


sub get_strand{
my ($strand) = @_;
return 1 if $strand eq "+";
return -1 if $strand eq "-";
}

sub get_chr {
my ($chr) = @_;
$chr =~s/chr//;
return "MT" if $chr eq  "M";
return $chr;
}


sub getId{
	my ($id) = @_;
	my ($n,$v) = split(/\./,$id);
	return $n;
	$id =~s/\.*//g;
	return $id;
}

sub parse_infos {
	my ($string) = @_;
	my $infos ={};
	foreach my$cp (split(";",$string)){
		my $code = " ";
		$code ="=" if $cp =~ /=/;
		my ($a,$b) = split("$code",$cp);
		$infos->{$a} = $b;
		
		
	}
	return $infos;
}

sub save_sqlite {
	my ($no,$obj,$type,$replace,$debug) = @_;
	 foreach my $gid (keys %$obj){
	 	next unless $obj->{$gid}->{id};
	 	if ($replace == 1){
	 	my $o = $no->get("annotations",$gid);
	 	warn "skip $gid";
	 	next if $o;
	 	}
	 	if ($replace == 2){
	 		my $o = $no->get("annotations",$gid);
	 		my ($name,$z1)=  $no->get_key_value("annotations",$gid);	
			#die($gid) unless $name;
			$no->delete_bulk("annotations",$name) if $name;
			
	 	#next if $o;
	 	}
	 
	 	#warn $gid." :: ".$obj->{$gid}->{id}." ".$obj->{$gid}->{sequence} if $debug && !($obj->{$gid}->{sequence});
	 	$obj->{$gid}->{kyotoId} = 	$obj->{$gid}->{id};
	 	warn "save $gid" if $replace ==2;
	 	
	 	my @c =  keys %{$obj->{$gid}->{chromosomes_object}};
	 	warn $gid if @c ne 1;
	 	die( Dumper $obj->{$gid}) if @c ne 1;
 		my @synonym = ($obj->{$gid}->{external_name},$obj->{$gid}->{name},$obj->{$gid}->{genbo_id}," genboid:".$obj->{$gid}->{genbo_id},$type,"chr:".$obj->{$gid}->{chromosome});
 		if (exists $obj->{$gid}->{target_genes}){
 			push(@synonym,@{$obj->{$gid}->{target_genes}});
 			
 		}
 		#warn Dumper $obj->{$gid}->{tag} if ($type  eq "transcript");
 	#	warn join(" ",@synonym);
# 	warn $obj->{$gid}->{ccds_name} if $obj->{$gid}->{ccds_name};
 	push(@synonym,$obj->{$gid}->{ccds_name}) if $obj->{$gid}->{ccds_name};
 	$no->put("annotations",join(" ",@synonym),$obj->{$gid});
 }
 
}


1;