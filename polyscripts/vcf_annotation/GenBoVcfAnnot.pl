#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$Bin/packages";
use lib "$Bin/../../GenBo/lib/";
use lib "$Bin/../../GenBo/lib/obj-nodb";
use strict;
use QueryOnlyVcf;
use Getopt::Long;
use Data::Dumper;
use String::ProgressBar;
use Compress::Snappy;
use Parallel::ForkManager;
use Carp;
use Storable qw/thaw freeze/;
use File::Temp qw/ :mktemp  /;

#perl ./GenBoVcfAnnot.pl -vcf_file=FILE.vcf.gz -tab_out=OUTFILE.tab -parse_all -fork=10

my ($help, $format, $fork, $fileName, $tab_outfile, $join_characters, $method, $getObjects, $addAnnotVcf, $filterTrans, @lFiltersAnnot, $filter_score_polyphen, $filter_score_sift, $filter_freq_5P, $filter_freq_1P, $filter_freq_01P, $filter_freq_001P,$no_dejavu );
my ($use_main_transcripts, $use_ccds_transcripts);
my $patientName = 'none';
my $force_all_gt_he;
my $release;
GetOptions(
	'help|h!'          => \$help,
	'vcf_file|f=s'     => \$fileName,
	'patient|p=s'      => \$patientName,
	'method|m=s'       => \$method,
	'getobj!'        => \$getObjects,
	'addAnnotVcf!'     => \$addAnnotVcf,
	'filterTrans=s'    => \$filterTrans,
	'annot=s'          => \@lFiltersAnnot,
	'score_polyphen=s' => \$filter_score_polyphen,
	'score_sift=s'     => \$filter_score_sift,
	'freq_5P!'         => \$filter_freq_5P,
	'freq_1P!'         => \$filter_freq_1P,
	'freq_01P!'        => \$filter_freq_01P,
	'freq_001P!'       => \$filter_freq_001P,
	'force_all_gt_he!' => \$force_all_gt_he,
	'parse_all!'	   => \$force_all_gt_he,
	'no_dejavu!'	   => \$no_dejavu,
	'tab_out=s'		   => \$tab_outfile,
	'fork=s'		   => \$fork,
	'join_characters|j=s' => \$join_characters,
	'use_main_transcripts=s' => \$use_main_transcripts,
	'use_ccds_transcripts=s' => \$use_ccds_transcripts,
	'format=s'	=> \$format,
	'release=s'	=> \$release,
);
die "\n\nNo -file or -f option... Die...\n\n" unless ($fileName);
die "\n\nNeed -tab_outfile option... Die...\n\n" unless ($tab_outfile);

confess("\n\nFile $fileName not found ! Die\n\n") if not -e $fileName;
confess("\n\nINDEX File $fileName not found ! Die\n\n") if not -e $fileName.'.tbi';
confess("\n\nFile $fileName not zip ! Use bgzip + tabix -p vcf. Die,\n\n") if not $fileName =~ /\.gz/;

my $buffer = new GBuffer;
my $genecode = $buffer->getQuery->getMaxGencodeVersion();
my $annotdb = $buffer->getQuery->getMaxPublicDatabaseVersion();
my $last_annot_version = $genecode.'.'.$annotdb;
my $proj_tmp_name;
if ($release) { $proj_tmp_name = $buffer->getRandomProjectName($release); }
else { $proj_tmp_name = $buffer->getRandomProjectName($release); }
my $proj_tmp = $buffer->newProject( -name => $proj_tmp_name );
my $getGenomeFasta = $proj_tmp->getGenomeFasta();
my $getGenomeFai = $proj_tmp->getGenomeFai(); 

if ($filter_score_polyphen) { $filter_score_polyphen =~ s/,/./; }
if ($filter_score_sift) { $filter_score_sift =~ s/,/./; }

unless ($join_characters) { $join_characters = "\t"; }
unless ($fork) { $fork = 1; }
unless ($format) { $format= 'by_transcript'; }


my ($hAnnot_authorized, $hFilters);
$hAnnot_authorized->{'Intergenic'} = undef;
$hAnnot_authorized->{'Intronic'} = undef;
$hAnnot_authorized->{'Pseudogene'} = undef;
$hAnnot_authorized->{'ncRNA'} = undef;
$hAnnot_authorized->{'mature-miRNA'} = undef;
$hAnnot_authorized->{'Utr'} = undef;
$hAnnot_authorized->{'Splice_Region'} = undef;
$hAnnot_authorized->{'Splice_Acceptor/Donor'} = undef;
$hAnnot_authorized->{'Missense'} = undef;
$hAnnot_authorized->{'Frameshift'} = undef;
$hAnnot_authorized->{'Stop-gained'} = undef;
$hAnnot_authorized->{'(Start/Stop)-lost'} = undef;
$hAnnot_authorized->{'Synonymous'} = undef;
$hAnnot_authorized->{'No-frameshift'} = undef;
foreach my $annot (@lFiltersAnnot) {
	unless (exists $hAnnot_authorized->{$annot}) {
		warn "\n\nERROR: $annot not a listed annotation... you can use theses annotations:\n\n";
		warn join("\n", keys %$hAnnot_authorized)."\n\nDie...\n\n";
		die;
	}
	$hFilters->{'annot'}->{$annot} = undef;
}

$getObjects = 1 if ($addAnnotVcf);
$getObjects = 1 if ($tab_outfile);

my %args;
$args{annotation_version} = $last_annot_version;
$args{file}        = $fileName;
$args{method}      = $method if ($method);
$args{getObjects}  = $getObjects if ($getObjects);
$args{genomeFasta} = $getGenomeFasta;
$args{genomeFai}   = $getGenomeFai;
$args{release}     = $release if ($release);
if ($patientName) {
	if ($patientName eq 'none') {
		$args{noPatient} = 1;
	}
	else {
		$args{patient}->{id} = $patientName;
		$args{patient}->{name} = $patientName
	}
}
if ($force_all_gt_he) { $args{force_all_gt_he} = 1; }

my ($hOnlyGenesTrans, $hOnlyGenesTransLog);
if ($filterTrans) {
	die("\n\nERROR: $filterTrans file doesn't exists... Die...\n\n") unless (-e $filterTrans);
	open(FILE, $filterTrans);
	my $hGenesTrans;
	my $i = 0;
	while(<FILE>) {
		chomp($_);
		my $line = $_;
		$i++;
		next if ($i == 1);
		my @lTmp = split(',', $line);
		$hOnlyGenesTrans->{$lTmp[0]}->{$lTmp[1]} = undef;
	}
	close(FILE);
}
warn "\n";
warn "Fork used: $fork\n";

my $hashRes;
my $nb_ref_done = 1;
my $nb_var_total = 0;
my $hNbVarTotal_byRef;
#my @lRef = @{$queryVcf->getReferences()};

my $tmpdir_name = mkdtemp('/tmp/tmp_genbovcfannot_XXXX');
print "TMP dir: ".$tmpdir_name."\n";

my $hChr_found;
my @lChrVcf = `tabix -l $fileName`;
foreach my $chr_id (@lChrVcf) {
	chomp($chr_id);
	$hChr_found->{$chr_id} = undef;
	$chr_id =~ s/chr//;
	$hChr_found->{$chr_id} = undef;
	if ($chr_id eq 'M') { $hChr_found->{'MT'} = undef; }
	if ($chr_id eq 'MT') { $hChr_found->{'M'} = undef; }
}

my (@lRef, @lPartFiles);
$proj_tmp->getChromosomes();
my $h_tmp_ref;
my $nb = 0;
foreach my $chr_id (1..22, 'X', 'Y', 'MT', 'M') {
	next if not exists $hChr_found->{$chr_id};
	my $chr = $proj_tmp->getChromosome($chr_id);
	foreach my $ref (@{$chr->getReferences($chr->start(), $chr->end(), 50)}) {
		push(@lRef, $ref);
		$h_tmp_ref->{$ref->id()} = $nb;
		my $tmp = $tmpdir_name.'/vcf_annot.part.'.$nb.'.tab';
		push(@lPartFiles, $tmp);
		$nb++;
	}
}

print "NB References: ".scalar(@lRef)."\n";

warn "\n";
warn "# Parsing / Annoting VCF\n";
warn '-> File: '.$fileName."\n";
my $queryVcf = QueryOnlyVcf -> new(\%args);
my $gencode_version = $queryVcf->project->gencode_version();
print "# Build: ".$queryVcf->project->getVersion()."\n";
print "# GenCode Version: $gencode_version\n";
my $version = $queryVcf->project->annotation_version();
print "# GenBo Annotation Version: $version\n";

$queryVcf->project();
$queryVcf->buffer->hash_genes_omim_morbid();

my @lTmp = split('\.', $version);
my $annot_version = $lTmp[1];


my (@lGlobal, $hok);
my $pm = new Parallel::ForkManager($fork);

foreach my $reference (@lRef) {
	my $pid = $pm->start() and next;
	print 'parsing '.$reference->id()."\n";
	$queryVcf->project->disconnect();
	$queryVcf->project->buffer->dbh_deconnect();
	my $res; 
	$res->{ok} = 1;
	$res->{data} = $queryVcf->parseVcfFile($reference);
	my $hok = annote_output_file($res->{data});
	my $nb = int($h_tmp_ref->{$reference->id()});
	open (TMP, '>'.$lPartFiles[$nb]);
	foreach my $chr_id (1..22, 'X', 'Y', 'MT', 'M') {
		next unless exists $hok->{$chr_id};
		
		my $hPos;
		foreach my $var_id (sort keys %{$hok->{$chr_id}}) {
			my @lTmp = split('_', $var_id);
			$hPos->{$lTmp[1]}->{$var_id} = $hok->{$chr_id}->{$var_id};
		}
		foreach my $pos (sort {$a <=> $b} keys %{$hPos}) {
			foreach my $var_id (sort keys %{$hPos->{$pos}}) {
				if ($tab_outfile) {
					print TMP $hPos->{$pos}->{$var_id}."\n";
				}
			}
		}
	}
	close (TMP);
	
	$pm->finish();
}
$pm->wait_all_children;

sleep(3);


print "\n\n-> Parsing / Annoting Done !\n\n";


open (FILE, '>'.$tab_outfile);
print FILE "# Build: ".$queryVcf->project->getVersion()."\n";
print FILE "# GenBo Annotation Version: $version\n";
foreach my $cat (sort keys %{$queryVcf->buffer->public_data->{$annot_version}}) {
	print FILE "# $cat: ".$queryVcf->buffer->public_data->{$annot_version}->{$cat}->{'version'}."\n";
}
close (FILE);


print 'cat tmp files ';
my $cmd_cat = "cat ".join(' ', @lPartFiles)." >>$tab_outfile";
`$cmd_cat`;
print "-> Ok!\n\n";

print 'bgzip file '; 
my $cmd_bgzip = qq{bgzip -f $tab_outfile};
`$cmd_bgzip`;
print "-> Ok!\n\n";

print 'tabix file '; 
my $cmd_tabix = qq{tabix -f -b 4 -e 5 -s 3 $tab_outfile.gz};
`$cmd_tabix`;
print "-> Ok!\n\n";

print 'rm tmpDir ';
`rm -r $tmpdir_name`;
print "-> Ok!\n\n";

print "\n\nJOB Done: $tab_outfile.gz\n\n";

$proj_tmp = undef;
$buffer = undef;

exit(0);


############################


sub getAllPatientsNameInVcf {
	my $file = shift;
	my @lPatNames;
	die("\n\nERROR: FILE ".$file." doesn't exists... Die...\n\n") unless (-e $file);
	if ($file =~ /.gz/) { open (FILE, "zcat ".$file." |"); }
	else { open (FILE, $file); }
	while (<FILE>) {
		my $line = $_;
		chomp($line);
		next unless ($line =~ /#CHROM/);
		my @lTmp = split(' ', $line);
		my $i = 9;
		while ($i < scalar(@lTmp)) {
			push(@lPatNames, $lTmp[$i]);
			$i++;
		}
		close(FILE);
		return \@lPatNames;
	}
	close(FILE);
	return;
}

sub annote_output_file {
	my ($hashRes) = @_;
	my $hres;
	my $buffer = $queryVcf->buffer();
	$queryVcf->project->{buffer} = $buffer;
	my $project = $queryVcf->project();
	foreach my $type (keys %{$hashRes}) {
		foreach my $var (@{$hashRes->{$type}}) {
			$var->{buffer} = undef; 
			$var->{project} = undef; 
			$var->{buffer} = $buffer;
			$var->{project} = $project;
			$var->getChromosome->{buffer} = $buffer;
			$var->getChromosome->{project} = $project;
			my @lFields;
			my @lTmp = split('_', $var->id());
			my $vcf_pos = $var->vcf_position();
			my $h;
			if ($tab_outfile) {
				$h->{id} = $var->id();
				$h->{type} = $var->type_public_db();
				$h->{chr} = 'chr'.$var->getChromosome->id();
				$h->{start} = $var->start();
				$h->{end} = $var->end();
				$h->{ref_allele} = $var->ref_allele();
				$h->{var_allele} = $var->var_allele();
				
				$h->{rs_name} = $var->rs_name();
				$h->{rs_name} = '-' unless ($h->{rs_name});
				
				$h->{hgmd_id} = $var->hgmd_id();
				$h->{hgmd_id} = '-' unless ($h->{hgmd_id});
				
				$h->{clinvar_id} = $var->clinvar_id();
				$h->{clinvar_id} = '-' unless ($h->{clinvar_id});
				
				$h->{rs_name} = $var->rs_name();
				$h->{rs_name} = '-' unless ($h->{rs_name});
				
				$h->{frequency} = sprintf("%.3f", $var->percent()).'%' if ($var->frequency());
				$h->{frequency} = '-' unless ($h->{frequency});
				
#				$h->{dejavu} = 'NA';
#				$h->{dejavu} = $var->other_projects().'/'.$var->other_patients() if not ($no_dejavu);
#				
#				warn "\n";
#				warn "\n";
#				warn Dumper $var-> dejaVuInfosForDiag2;
				my $h_dv = $var->getChromosome->rocks_dejavu->dejavu($var->rocksdb_id);
				my $dv_proj = 0;
				my $dv_samples_he = 0;
				my $dv_samples_ho = 0;
				foreach my $proj_id (keys %$h_dv) {
					$dv_proj++;
					$dv_samples_he += $h_dv->{$proj_id}->{he};
					$dv_samples_ho += $h_dv->{$proj_id}->{ho};
				}
				my $dv_samples = $dv_samples_he + $dv_samples_ho;
				$h->{dejavu} = $dv_proj.'/'.$dv_samples.'('.$dv_samples_ho.'ho)';
				
				
				
				$h->{hgmd_class} = $var->hgmd_class();
				$h->{hgmd_class} = '-' unless ($h->{hgmd_class});
				
				$h->{cadd_score} = $var->cadd_score();
				$h->{cadd_score} = '-' if not defined ($h->{cadd_score});
				
				$h->{ncboost_score} = $var->ncboost_score();
				$h->{ncboost_score} = '-' if not defined ($h->{ncboost_score});
				
				$h->{class_clinvar} = '-';
				$h->{class_clinvar} = $var->text_clinvar() if $var->text_clinvar();
				
				$h->{score_clinvar} = $var->score_clinvar();
				$h->{score_clinvar} = '-' if not defined ($h->{score_clinvar});
				
				$h->{isClinvarPathogenic} = $var->isClinvarPathogenic();
				$h->{isClinvarPathogenic} = '-' unless ($h->{isClinvarPathogenic});
				
				my $hTrIds;
				foreach my $tr (@{$var->getTranscripts()}) {
					next if ($use_main_transcripts and not $tr->isMain());
					next if ($use_ccds_transcripts and not $tr->ccds_name());
					$hTrIds->{$tr->id()} = undef;
				}
				
				foreach my $gene (@{$var->getGenes()}) {
					$h->{genes}->{$gene->id()}->{external_name} = $gene->external_name();
					foreach my $tr (@{$gene->getTranscripts()}) {
						next unless (exists $hTrIds->{$tr->id});
						$h->{genes}->{$gene->id()}->{transcripts}->{$tr->id} = undef;
					}
					# SPLICE_AI
					my $h_score_spliceAI = $var->spliceAI_score($gene);
					if ($h_score_spliceAI and exists $h_score_spliceAI->{'AG'}) {
						my @l_score_spliceAI;
						foreach my $cat (sort keys %$h_score_spliceAI) {
							push(@l_score_spliceAI, $cat.'='.$h_score_spliceAI->{$cat});
						}
						$h->{genes}->{$gene->id()}->{spliceAI_score} = join(',', @l_score_spliceAI);
					}
					else { $h->{genes}->{$gene->id()}->{spliceAI_score} = '-'; }
					# OMIM
					if ($gene->omim_id()) {
						$h->{genes}->{$gene->id()}->{omim_id} = $gene->omim_id();
						$h->{genes}->{$gene->id()}->{omim_inheritance} = $gene->omim_inheritance();
						$h->{genes}->{$gene->id()}->{pLI} = $gene->pLI();
					}
					else {
						$h->{genes}->{$gene->id()}->{omim_id} = '-';
						$h->{genes}->{$gene->id()}->{omim_inheritance} = '-';
						$h->{genes}->{$gene->id()}->{pLI} = '-';
					}
				}
				
				foreach my $tr (@{$var->getTranscripts()}) {
					next unless (exists $hTrIds->{$tr->id});
					next if ($use_main_transcripts and not $tr->isMain());
					next if ($use_ccds_transcripts and not $tr->ccds_name());
					$h->{transcripts}->{$tr->id()}->{annotations} = $var->variationTypeInterface($tr);
					
					if ($tr->getProtein) {
						$h->{transcripts}->{$tr->id()}->{alphamissense} = $var->alphamissense($tr);
						
#						my $polyphen_status = $var->polyphenStatus($tr);
#						if ($polyphen_status) { $polyphen_status =~ s/ /_/g; }
#						else { $polyphen_status = '-'; }
#						$h->{transcripts}->{$tr->id()}->{polyphen_status} = $polyphen_status;
						
						
#						my $sift_status = $var->siftStatus($tr);
#						if ($sift_status) { $sift_status =~ s/ /_/g; }
#						else { $sift_status = '-'; }
#						$h->{transcripts}->{$tr->id()}->{sift_status} = $sift_status;
						
						my $polyphen_score = '-';
						my $sift_score = '-';
						if ($var->isVariation()) {
							$polyphen_score = $var->polyphenScore($tr);
							$h->{transcripts}->{$tr->id()}->{polyphen_score} = $polyphen_score;
							$sift_score = $var->siftScore($tr);
							$h->{transcripts}->{$tr->id()}->{sift_score} = $sift_score;
						}
						
						
						
						
						
					}
					else {
						$h->{transcripts}->{$tr->id()}->{polyphen_status} = '-';
						$h->{transcripts}->{$tr->id()}->{polyphen_score} = '-';
						$h->{transcripts}->{$tr->id()}->{sift_status} = '-';
						$h->{transcripts}->{$tr->id()}->{sift_score} = '-';
						$h->{transcripts}->{$tr->id()}->{alphamissense} = '-';
					}
				}
				$h->{id} = $var->id();
				$h->{type} = $var->type_public_db();
				$h->{chr} = 'chr'.$var->getChromosome->id();
				$h->{start} = $var->start();
				$h->{end} = $var->end();
				$h->{ref_allele} = $var->ref_allele();
				$h->{var_allele} = $var->var_allele();
				
				if ($format eq 'by_transcript') {
					my @lCol;
					push(@lCol, $h->{id});
					push(@lCol, $h->{type});
					push(@lCol, $h->{chr});
					push(@lCol, $h->{start});
					push(@lCol, $h->{end});
					push(@lCol, $h->{ref_allele});
					push(@lCol, $h->{var_allele});
					push(@lCol, 'id:'.$h->{id});
					push(@lCol, 'rs:'.$h->{rs_name});
					push(@lCol, 'hgmd_id:'.$h->{hgmd_id});
					push(@lCol, 'clinvar_id:'.$h->{clinvar_id});
					push(@lCol, 'db_freq:'.$h->{frequency});
					push(@lCol, 'dejavu:'.$h->{dejavu});
					push(@lCol, 'hgmd_class:'.$h->{hgmd_class});
					push(@lCol, 'cadd:'.$h->{cadd_score});
					push(@lCol, 'ncboost:'.$h->{ncboost_score});
					push(@lCol, 'clinvar_class:'.$h->{class_clinvar});
					push(@lCol, 'clinvar_score:'.$h->{score_clinvar});
					push(@lCol, 'clinvar_pathogenic:'.$h->{isClinvarPathogenic});
					if ($h->{genes}) {
						foreach my $gene_id (sort keys %{$h->{genes}}) {
							my @lCol_g;
							push(@lCol_g, 'gene:'.$gene_id.'|'.$h->{genes}->{$gene_id}->{external_name});
							push(@lCol_g, 'omim_id:'.$h->{genes}->{$gene_id}->{omim_id});
							push(@lCol_g, 'omim_inheritance:'.$h->{genes}->{$gene_id}->{omim_inheritance});
							push(@lCol_g, 'pLI:'.$h->{genes}->{$gene_id}->{pLI});
							push(@lCol_g, 'spliceAI_score:'.$h->{genes}->{$gene_id}->{spliceAI_score});
							if ($h->{genes}->{$gene_id}->{transcripts}) {
								foreach my $tr_id (sort keys %{$h->{genes}->{$gene_id}->{transcripts}}) {
									my @lCol_tr;
									push(@lCol_tr, $join_characters);
									push(@lCol_tr, $tr_id.':'.$h->{transcripts}->{$tr_id}->{annotations});
									#push(@lCol_tr, 'polyphen_status:'.$h->{transcripts}->{$tr_id}->{polyphen_status});
									push(@lCol_tr, 'polyphen_score:'.$h->{transcripts}->{$tr_id}->{polyphen_score});
									#push(@lCol_tr, 'sift_status:'.$h->{transcripts}->{$tr_id}->{sift_status});
									push(@lCol_tr, 'sift_score:'.$h->{transcripts}->{$tr_id}->{sift_score});
									push(@lCol_tr, 'alphamissense:'.$h->{transcripts}->{$tr_id}->{alphamissense});
									$hres->{$var->getChromosome->id()}->{$var->id().'_'.$gene_id.'_'.$tr_id}  = join($join_characters, @lCol).$join_characters.join($join_characters, @lCol_g).$join_characters.join($join_characters, @lCol_tr);
								}
							}
							else {
								$hres->{$var->getChromosome->id()}->{$var->id().'_'.$gene_id}  = join($join_characters, @lCol).$join_characters.join($join_characters, @lCol_g).$join_characters.'no_transcript_found';
							}
						}
					}
					else {
						$hres->{$var->getChromosome->id()}->{$var->id().'_intronic'}  = join($join_characters, @lCol).$join_characters.'intronic';
					}
				}
			}
			$h = undef;
			$buffer = $var->{buffer};
			$project = $var->{project};
			$var->{buffer} = undef;
			$var->{project} = undef;
			$var = undef;
		}
	}
	return $hres;
}
