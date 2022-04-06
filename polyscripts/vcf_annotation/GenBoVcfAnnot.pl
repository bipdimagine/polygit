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
use Storable qw/thaw freeze/;

my ($help, $format, $fork, $fileName, $tab_outfile, $join_characters, $method, $getObjects, $addAnnotVcf, $filterTrans, @lFiltersAnnot, $filter_score_polyphen, $filter_score_sift, $filter_freq_5P, $filter_freq_1P, $filter_freq_01P, $filter_freq_001P,$no_dejavu );
my ($use_main_transcripts, $use_ccds_transcripts);
my $patientName = 'none';
my $force_all_gt_he;
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
);
die "\n\nNo -file or -f option... Die...\n\n" unless ($fileName);

my $buffer = new GBuffer;
my $genecode = $buffer->getQuery->getMaxGencodeVersion();
my $annotdb = $buffer->getQuery->getMaxPublicDatabaseVersion();
my $last_annot_version = $genecode.'.'.$annotdb;
my $proj_tmp = $buffer->newProject( -name => $buffer->get_random_project_name_with_this_annotations_and_genecode() );
my $getGenomeFasta = $proj_tmp->getGenomeFasta();
my $getGenomeFai = $proj_tmp->getGenomeFai();
$proj_tmp = undef;
$buffer = undef; 

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

warn "\n";
warn "# Parsing / Annoting VCF\n";
warn '-> File: '.$fileName."\n";
my $queryVcf = QueryOnlyVcf -> new(\%args);
my $gencode_version = $queryVcf->project->gencode_version();
print "# GenCode Version: $gencode_version\n";
my $version = $queryVcf->project->annotation_version();
print "# GenBo Annotation Version: $version\n";


$queryVcf->project();
$queryVcf->buffer->hash_genes_omim_morbid();

my $hashRes;
my $nb_ref_done = 1;
my $nb_var_total = 0;
my $hNbVarTotal_byRef;
my @lRef = @{$queryVcf->getReferences()};

my (@lGlobal, $hok);
my $pm = new Parallel::ForkManager($fork);
$pm->run_on_finish(
	sub {
		my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data ) = @_;
		foreach my $chr_id (keys %$data) {
			$hok->{$chr_id} = $data->{$chr_id};
			warn "chr$chr_id annotation done.\n";
		}
	}
);

foreach my $reference (@lRef) {
	my $pid = $pm->start() and next;
	my $res; 
	$res->{ok} = 1;
	$res->{data} = $queryVcf->parseVcfFile($reference);
	$res->{id} = $reference->getChromosome->id();
	my $chr_id = $reference->getChromosome->id();
	my $hres = annote_output_file($res->{id}, $res->{data});
	$pm->finish( 0, $hres );
}
$pm->wait_all_children;

sleep(3);

if ($tab_outfile) {
	open (FILE, '>'.$tab_outfile);
	my $version = $queryVcf->project->annotation_version();
	print FILE "# GenBo Annotation Version: $version\n";
	my @lTmp = split('\.', $version);
	my $annot_version = $lTmp[1];
	foreach my $cat (sort keys %{$queryVcf->buffer->public_data->{$annot_version}}) {
		print FILE "# $cat: ".$queryVcf->buffer->public_data->{$annot_version}->{$cat}->{'version'}."\n";
	}
	foreach my $chr_id (1..22, 'X', 'Y', 'MT', 'M') {
		next unless exists $hok->{$chr_id};
		
		my $hPos;
		foreach my $var_id (sort keys %{$hok->{$chr_id}}) {
			my @lTmp = split('_', $var_id);
			$hPos->{$lTmp[1]}->{$var_id} = $hok->{$chr_id}->{$var_id};
		}
		foreach my $pos (sort {$a <=> $b} keys %{$hPos}) {
			foreach my $var_id (sort keys %{$hPos->{$pos}}) {
				print FILE $hPos->{$pos}->{$var_id}."\n";
			}
		}
	}
	close (FILE);
}

warn "\n\n-> Parsing / Annoting Done !\n\n";

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
	my ($part_id, $hashRes) = @_;
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
				
				$h->{dejavu} = 'NA';
				$h->{dejavu} = $var->nb_deja_vu_projects().'/'.$var->nb_deja_vu_samples() unless ($no_dejavu);
				
				$h->{hgmd_class} = $var->hgmd_class();
				$h->{hgmd_class} = '-' unless ($h->{hgmd_class});
				
				$h->{cadd_score} = $var->cadd_score();
				$h->{cadd_score} = '-' unless ($h->{cadd_score});
				
				$h->{ncboost_score} = $var->ncboost_score();
				$h->{ncboost_score} = '-' unless ($h->{ncboost_score});
				
				$h->{score_clinvar} = $var->score_clinvar();
				$h->{score_clinvar} = '-' unless ($h->{score_clinvar});
				
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
					if ($h_score_spliceAI) {
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
						$h->{transcripts}->{$tr->id()}->{polyphen_status} = $var->polyphenStatusText($tr);
						$h->{transcripts}->{$tr->id()}->{polyphen_status} = '-' unless ($h->{transcripts}->{$tr->id()}->{polyphen_status});
						$h->{transcripts}->{$tr->id()}->{polyphen_score} = $var->polyphenScore($tr);
						$h->{transcripts}->{$tr->id()}->{polyphen_score} = '-' unless ($h->{transcripts}->{$tr->id()}->{polyphen_score});
						$h->{transcripts}->{$tr->id()}->{sift_status} = $var->siftStatusText($tr);
						$h->{transcripts}->{$tr->id()}->{sift_status} = '-' unless ($h->{transcripts}->{$tr->id()}->{sift_status});
						$h->{transcripts}->{$tr->id()}->{sift_score} = $var->siftScore($tr);
						$h->{transcripts}->{$tr->id()}->{sift_score} = '-' unless ($h->{transcripts}->{$tr->id()}->{sift_score});
					}
					else {
						$h->{transcripts}->{$tr->id()}->{polyphen_status} = '-';
						$h->{transcripts}->{$tr->id()}->{polyphen_score} = '-';
						$h->{transcripts}->{$tr->id()}->{sift_status} = '-';
						$h->{transcripts}->{$tr->id()}->{sift_score} = '-';
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
									push(@lCol_tr, 'polyphen_status:'.$h->{transcripts}->{$tr_id}->{polyphen_status});
									push(@lCol_tr, 'polyphen_score:'.$h->{transcripts}->{$tr_id}->{polyphen_score});
									push(@lCol_tr, 'sift_status:'.$h->{transcripts}->{$tr_id}->{sift_status});
									push(@lCol_tr, 'sift_score:'.$h->{transcripts}->{$tr_id}->{sift_score});
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
