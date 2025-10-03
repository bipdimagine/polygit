#!/bin/bash

#echo $pad &&samtools bedcov -H -d 30 RENOME_V5_Hg19.padding_regions_$pad.bed align/dragen-align/GIAB_V5.bam > align/coverage/GIAB_V5.padding_regions_$pad.bed.cov


pad=50
echo "pad=$pad"
for project in NGS2025_08659 NGS2025_08660 NGS2025_08679 NGS2025_08680 NGS2025_08964 NGS2025_09282
  do
	echo "$project"
	dir="/data-isilon/sequencing/ngs/$project/HG19_MT/variations/"
#	bedB="/data-isilon/sequencing/ngs/NGS2025_08659/HG19_MT/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed"
	bedB="/data-isilon/sequencing/ngs/NGS2025_08659/HG19_MT/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.with_chr.bed"
	if [[ ! -e "$bedB" ]]; then
		sed 's/^/chr/; s/chrMT/chrM/' "/data-isilon/sequencing/ngs/NGS2025_08659/HG19_MT/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed" > "$bedB"
	fi
#	echo "$bedB"
	
	# RENOME_V5_Hg19
	if [[ "$project" == "NGS2025_08659" || "$project" == "NGS2025_08660" ]]; then
		patient="GIAB_V5"
		bedA="/data-isilon/sequencing/ngs/NGS2025_08659/HG19_MT/RENOME_V5_Hg19.bed"
		bed_intersect="/data-isilon/sequencing/ngs/NGS2025_08659/HG19_MT/intersect_renome_giab.bed"
		
	# Rapid_TubeF_V4_UMI
	elif [[ "$project" == "NGS2025_08679" || "$project" == "NGS2025_08680" || "$project" == "NGS2025_08964" ]]; then
		patient="GIA_HG0_6123GM002954_STB"
#		bedA="/data-isilon/sequencing/ngs/NGS2025_08679/HG19_MT/Rapid_V4_TubeF_UMI.no_chr.bed"
		bedA="/data-isilon/public-data/capture/HG19/agilent/Rapid_V4_TubeF_UMI.bed"
		bed_intersect="/data-isilon/sequencing/ngs/NGS2025_08679/HG19_MT/intersect_tubeF_giab.bed"

	# Exome agilent_50_V5
	elif [[ "$project" == "NGS2025_09282" ]]; then
		patient="AJ_SON_CTRL_dragen_rep1"
#		bedA="/data-isilon/sequencing/ngs/NGS2025_09282/HG19_MT/agilent.50.v5.no_chr.bed"
		bedA="/data-isilon/public-data/capture/HG19/agilent/agilent.50.v5.bed"
		bed_intersect="/data-isilon/sequencing/ngs/NGS2025_09282/HG19_MT/intersect_agilentV5_giab.bed"
	fi
	
	echo "$patient"
	ls -lh "$bedA"
	ls -lh "$bedB"
	ls -lh "$bed_intersect"
	
	bed_extended="$(echo "$bed_intersect" | sed s/\.bed$/\.extended_$pad\.bed/)"
	if [[ "$project" == "NGS2025_08659" || "$project" == "NGS2025_08679"  || "$project" == "NGS2025_09082" ]]; then
		bedtools intersect -a "$bedA" -b "$bedB" > "$bed_intersect"
		while read -r chr start stop name; do
			echo -e "$chr\t$(($start-$pad))\t$(($stop+$pad))"
		done < "$bed_intersect" > "$bed_extended"
		head -n 3 "$bed_extended"
		link="$(echo $bed_extended | sed s/NGS2025_08659/NGS2025_08660/ | sed s/NGS2025_08679/NGS2025_08680/ )"
		if [[ -L "$link" ]]; then rm "$link" ; fi
		ln -s "$bed_extended" "$link"
		if [[ "$project" == "NGS2025_08679" ]]; then
			link="$(echo $bed_extended | sed s/NGS2025_08679/NGS2025_08964/ )"
			if [[ -L "$link" ]]; then rm "$link" ; fi
			ln -s "$bed_extended" "$link"
		fi
		
		bed_no_chr="$(echo "$bed_extended" | sed s/\.bed$/\.no_chr\.bed/)"
		sed 's/^chrM/MT/;s/^chr//' "$bed_extended" > "$bed_no_chr"
		if [[  ! -e "$dir""giab/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz" ]]; then
			ln -s "/data-isilon/sequencing/ngs/NGS2025_08659/HG19_MT/variations/giab/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz" "$dir""giab/"
		fi
		bcftools view -T "$bed_no_chr" "$dir""giab/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz" -O z -o "$dir""giab/HG002.vcf.gz"
		bcftools index -ft "$dir""giab/HG002.vcf.gz"
		link="$(echo "$dir""giab/HG002.vcf.gz" | sed s/NGS2025_08659/NGS2025_08660/ | sed s/NGS2025_08679/NGS2025_08680/ )"
		if [[ -L "$link" ]]; then rm "$link" ; fi
		ln -s "$dir""giab/HG002.vcf.gz" "$link"
		if [[ -L "$link"".tbi" ]]; then rm "$link"".tbi" ; fi
		ln -s "$dir""giab/HG002.vcf.gz"".tbi" "$link"".tbi"
	fi
	ls -lh "$bed_extended"
	
#	echo "$dir"
	cd "$dir"
	
	for caller in dragen-calling haplotypecaller4 samtools freebayes unifiedgenotyper melt ; #deepvariant
	do
		vcf_original="$dir$caller/$patient.original.vcf.gz"
		if [[ ! -f "$vcf_original" ]]; then
			mv "$dir$caller/$patient.vcf.gz" "$vcf_original"
			mv "$dir$caller/$patient.vcf.gz.tbi" "$vcf_original"".tbi"
		fi
#		echo $vcf_original
		vcf_restreint="$(echo $vcf_original | sed s/\.original\.vcf\.gz$/\.restreint\.vcf\.gz/)"
		echo $vcf_restreint
		
#		echo "bcftools view -T $bed_extended $vcf_original -O z -o $vcf_restreint"
		bcftools view "$vcf_original" -T "$bed_extended" -O z -o "$vcf_restreint" 
		
#		echo "bcftools index -ft $vcf_restreint"
		bcftools index -ft "$vcf_restreint"
		
		if [[ -L "$dir$caller/$patient.vcf.gz" ]]; then
			rm "$dir$caller/$patient.vcf.gz"
		fi
		ln -s "$vcf_restreint" "$caller/$patient.vcf.gz"
	
		if [[ -L "$dir$caller/$patient.vcf.gz.tbi" ]]; then
			rm "$dir$caller/$patient.vcf.gz.tbi"
		fi
		ln -s "$vcf_restreint.tbi" "$caller/$patient.vcf.gz.tbi"
		
#		zgrep -E '^##bcftools_viewCommand=view -T ' "$vcf_restreint"
	 done
	echo
  done

