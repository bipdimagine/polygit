#!/usr/bin/perl
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/pedigree";

use GBuffer;
use Data::Dumper;
use export_excel;  
use export_data;
use Sys::Hostname;
use get_variations;
use pedigree;  
use validationQuery; 
use Set::IntSpan::Fast::XS;
use Storable qw/freeze thaw nfreeze nstore_fd nstore retrieve/;
use JSON;


my %types = ( variations => "variations" );

my $cgi    = new CGI();
my $buffer = new GBuffer;



my $host = hostname;

my $project_name = $cgi->param('project');
my $user = $cgi->param('user');
my $chr_id = $cgi->param('chromosome');
my $vector_ids = $cgi->param('vector_ids');
my $v_resume = $cgi->param('vector_resume');
my $type_cache = $cgi->param('type_cache');
my $type_name = $cgi->param('type');
$type_name = "variations";
my $type2    = $cgi->param('type2');
my $select_ids   = [split(",",$cgi->param('ids'))];
my $famind = $cgi->param('famind');
my $only_genes = $cgi->param('only_genes');
my $filter_nbvar_regionho = $cgi->param('filter_nbvar_regionho');

my $patients = $cgi->param('patients');

my $project = $buffer->newProjectCache( -name => $project_name );
die( "unknown project" . $project_name ) unless $project;
my @final_data;
my $stream;
$filter_nbvar_regionho = $project->ho_regions_short_value() unless ($filter_nbvar_regionho);

if ($cgi->param('gene') and $patients) {
	
	my $hGenes;
	my @lGenes = split('\+', $cgi->param('gene'));
	foreach my $g_name (@lGenes) { $hGenes->{$g_name} = undef; }
	
	my $hPat;
	my @lPat = split(',', $patients);
	foreach my $pat_name (@lPat) { $hPat->{$pat_name} = undef; }
	my $projectCache = $buffer->newProjectCache( -name 	=> $project_name );
	my $chr = $projectCache->getChromosome($chr_id);
	
	my $vector_pat;
	foreach my $patient (@{$project->getPatients()}) {
		next unless (exists $hPat->{$patient->name()});
		unless ($vector_pat) { $vector_pat = $patient->getVariantsVector($chr)->Clone(); }
		else { $vector_pat += $patient->getVariantsVector($chr); }
	}
	
	foreach my $gene (@{$chr->getGenes()}) {
		next unless (exists $hGenes->{$gene->external_name()});
		my $vector = $gene->getVariantsVector();
		$vector->Intersection($vector, $vector_pat);
		$vector_ids = join(',', @{$chr->getListVarVectorIds($vector)});
		
		#TODO: multi gene -> pb il faut accumuler les ids var pour prendre 
		#last;
	}
}

#if (($vector_ids or $v_resume) and $type_cache) {
if (defined($vector_ids) or $v_resume) {
	unless ( $cgi->param('xls') == 1 ) {
		$stream =1;
		$| =1;
		export_data::print_header_stream($cgi);
	}

	my $typeFilters = 'individual';
	$typeFilters = 'familial' if ($famind eq 'Fam');
	$typeFilters = 'somatic' if ($famind eq 'Som');
	my $projectCache = $buffer->newProjectCache( -name 	=> $project_name,
												 -cache 	=> '1',
												 -typeFilters => $typeFilters, );
	die( "unknown project" . $project_name ) unless $projectCache;
	#$projectCache->subdata_filter($type_cache);
	my (@lChr, @lIds, $h_ids);
	if ($cgi->param('chromosome') eq 'all') {
		foreach my $tmp (split('::', $v_resume)) {
			my @lTmp = split(':', $tmp);
			foreach my $id (split(',', $lTmp[-1])) { 
				$h_ids->{$lTmp[0]}->{$id} = undef;
			}
		}
		@lChr = sort {$a <=> $b} keys %$h_ids;
	}
	else {
		@lChr = split(" ", $chr_id);
		@lIds = split(',', $vector_ids);
	}
	$projectCache->getChromosomes();
	foreach my $chr_name (@lChr) {
		my $chr_cache = $projectCache->getChromosome($chr_name);
		my $hRegion;
		foreach my $type ('short', 'medium', 'large', 'extra_large') {
			my $nb_var = $projectCache->getNbVarHoRegionsFromType($type);
			$hRegion->{$type} = $chr_cache->getNewVector();
			foreach my $patient (@{$chr_cache->getPatients()}) {
				$hRegion->{$type} += $patient->getRegionHo($chr_cache, $nb_var);
			}
		}
		my $hIds;
		if ($cgi->param('chromosome') eq 'all') {
			@lIds = sort {$a <=> $b} keys %{$h_ids->{$chr_name}};
		}
		my $vector_all = Bit::Vector->new_Enum($chr_cache->getVariantsVector->Size(), join(',', @lIds));
		$chr_cache->setVariantsVector($vector_all);
		
		foreach my $v_id (@lIds) {
			my $var = $chr_cache->getVarObject($v_id);
			my $id = $var->id();
			#warn ref($var);
			$hIds->{$id}->{v_id} = $v_id;
			$hIds->{$id}->{chr_name} = $chr_cache->name;
			$hIds->{$id}->{he} = 0;
			$hIds->{$id}->{ho} = 0;
			$hIds->{$id}->{ho_region} = 0;
			foreach my $type ('extra_large', 'large', 'medium', 'short') {
				if ($hRegion->{$type}->contains($v_id)) {
					$hIds->{$id}->{ho_region} = $type;
					last;
				}
			}
			$hIds->{$id}->{nb_patients} = 0;
			foreach my $patient (@{$chr_cache->getPatients()}) {
				$patient->getVariantsVector($chr_cache);
				if ($patient->getHe($chr_cache)->contains($v_id)) {
					$hIds->{$id}->{he} = 1;
					$hIds->{$id}->{nb_patients}++;
				}
				if ($patient->getHo($chr_cache)->contains($v_id)) {
					$hIds->{$id}->{ho} = 1;
					$hIds->{$id}->{nb_patients}++;
				}
				next if (exists $hIds->{$id}->{json_ped}->{$patient->family()});
				next unless ($patient->family()); 
				my $h;
				my $ok = 0;
				foreach my $family (@{$chr_cache->getFamilies()}) {
					my $nb_parent = 0;
					my $nb_child = 0;
					next unless ($family->name() eq $patient->family());
					$h->{familly}->{name} = $family->name();
					$h->{familly}->{he_ho} = 'fam';
					$h->{familly}->{nb_ref} = 'fam';
					$h->{familly}->{nb_mut} = 'fam';
					foreach my $patient (@{$family->getPatients()}) {
						if (exists $var->sequencing_infos->{$patient->id()}) {
							$ok = 1;
							last;
						}
					}
					next unless ($ok);
					foreach my $patient_name (keys %{$family->parents()}) {
						$nb_parent++;
						my $patient = $chr_cache->getPatient($patient_name);
						my $h_parent;
						$h_parent->{name} = $patient_name;
						$h_parent->{nb_parent} = $nb_parent;
						$h_parent->{is_parent_child} = 'parent';
						if ($var->isVariation()) { $h_parent->{type} = 'snp'; }
						else { $h_parent->{type} = 'indel'; }
						$h_parent->{ref_allele} = $var->ref_allele();
						$h_parent->{var_allele} = $var->var_allele();
						$h_parent->{sex} = $patient->sex();
						$h_parent->{status} = 'healthy' if ($patient->status() == 1);
						$h_parent->{status} = 'ill' if ($patient->status() == 2);				
						if ($patient->getVariantsVector($chr_cache)->contains($v_id)) {
							my $he_ho ;
							$he_ho = 'He' if ($var->isHeterozygote($patient));
							$he_ho = 'Ho' if ($var->isHomozygote($patient));
							$h_parent->{he_ho} = $he_ho;
							if (exists $var->{seq}) {
								$h_parent->{nb_ref} = $var->{seq}->{$patient->id()}->{nref};
								$h_parent->{nb_mut} = $var->{seq}->{$patient->id()}->{nalt};
							}
							else {
								$h_parent->{nb_ref} = $var->getNbAlleleRef($patient);
								$h_parent->{nb_mut} = $var->getNbAlleleAlt($patient);
							}
						}
						else {
							my $cov;
							my $tabixFile = $patient->tabix_coverage();
							if ($tabixFile) {
								eval {
									my $res = $tabixFile->query('chr'.$chr_cache->id(), $var->start()-1, ($var->start() +1));
									while (my $line = $tabixFile->read($res)) {
										my ($thisChr, $thisPos, $thisCov) = split("\t", $line);
										if ($thisPos == $var->start()) {
											$cov = $thisCov;
										}
									}
									if ($cov) {
										$h_parent->{he_ho} = 'ref';
										$h_parent->{dp} = $cov;
										$h_parent->{nb_ref} = 'not';
										$h_parent->{nb_mut} = 'not';
									}
								};
								if ($@) {
									$h_parent->{he_ho} = 'pb_cov';
									$h_parent->{nb_ref} = 'pb_cov';
									$h_parent->{nb_mut} = 'pb_cov';
								}
							}
							else {
								$h_parent->{he_ho} = 'not';
								$h_parent->{nb_ref} = 'not';
								$h_parent->{nb_mut} = 'not';
							}
						}
						push(@{$h->{parents}}, $h_parent);
					}
					foreach my $patient_name (keys %{$family->children()}) {
						$nb_child++;
						my $patient = $chr_cache->getPatient($patient_name);
						my $h_child;
						$h_child->{name} = $patient_name;
						$h_child->{nb_child} = $nb_child;
						$h_child->{max_child} = scalar(keys %{$family->children()});
						$h_child->{is_parent_child} = 'child';
						if ($var->isVariation()) { $h_child->{type} = 'snp'; }
						else { $h_child->{type} = 'indel'; }
						$h_child->{ref_allele} = $var->ref_allele();
						$h_child->{var_allele} = $var->var_allele();
						$h_child->{sex} = $patient->sex();
						$h_child->{status} = 'healthy' if ($patient->status() == 1);
						$h_child->{status} = 'ill' if ($patient->status() == 2);
						
						my $transmission = $var->getTransmissionModel($patient->getFamily(), $patient);
						if ($transmission eq 'Father') { $h_child->{transmission} = $transmission; }
						elsif ($transmission eq 'Mother') { $h_child->{transmission} = $transmission; }
						elsif ($transmission eq 'Both') { $h_child->{transmission} = $transmission; }
						elsif ($transmission eq '?') { $h_child->{transmission} = $transmission; }
						else { $h_child->{transmission} = 'Model'; }
						
						if ($patient->getVariantsVector($chr_cache)->contains($v_id)) {
							my $he_ho ;
							$he_ho = 'He' if ($var->isHeterozygote($patient));
							$he_ho = 'Ho' if ($var->isHomozygote($patient));
							$h_child->{he_ho} = $he_ho;
							if (exists $var->{seq}) {
								$h_child->{nb_ref} = $var->{seq}->{$patient->id()}->{nref};
								$h_child->{nb_mut} = $var->{seq}->{$patient->id()}->{nalt};
							}
							else {
								$h_child->{nb_ref} = $var->getNbAlleleRef($patient);
								$h_child->{nb_mut} = $var->getNbAlleleAlt($patient);
							}
						}
						else {
							my ($cov, $nb);
							my $tabixFile = $patient->tabix_coverage();
							if ($tabixFile) {
								eval {
									my $res = $tabixFile->query('chr'.$chr_cache->id(), $var->start()-1, $var->end()+1);
									while (my $line = $tabixFile->read($res)) {
										my ($thisChr, $thisPos, $thisCov) = split("\t", $line);
										$cov += $thisCov;
										$nb++;
									}
									if ($cov and $nb) {
										$h_child->{he_ho} = 'ho_ref';
										$h_child->{dp} = int($cov/$nb);
										$h_child->{nb_ref} = 'not';
										$h_child->{nb_mut} = 'not';
									}
								};
								if ($@) {
									$h_child->{he_ho} = 'pb_cov';
									$h_child->{nb_ref} = 'pb_cov';
									$h_child->{nb_mut} = 'pb_cov';
								}
							}
							else {
								$h_child->{he_ho} = 'not';
								$h_child->{nb_ref} = 'not';
								$h_child->{nb_mut} = 'not';
							}
						}
						push(@{$h->{children}}, $h_child);
					}
				}
				next unless ($ok);
				my $nb_parent = 0;
				$nb_parent = scalar(@{$h->{parents}}) if (exists $h->{parents});
				my ($has_mother, $has_father);
				foreach my $h_parent (@{$h->{parents}}) {
					$has_father = 1 if ($h_parent->{sex} eq '1');
					$has_mother = 1 if ($h_parent->{sex} eq '2');
				}
				while ($nb_parent < 2) {
					$nb_parent++;
					my $h_parent;
					$h_parent->{name} = 'Undefined';
					$h_parent->{nb_parent} = $nb_parent;
					$h_parent->{is_parent_child} = 'parent';
					if ($var->isVariation()) { $h_parent->{type} = 'snp'; }
					else { $h_parent->{type} = 'indel'; }
					$h_parent->{ref_allele} = '?';
					$h_parent->{var_allele} = '?';
					if (not $has_father) {
						$h_parent->{sex} = '1';
						$has_father = 1;
					}
					elsif (not $has_mother) {
						$h_parent->{sex} = '2';
						$has_mother = 1;
					}
					$h_parent->{status} = 'healthy';
					$h_parent->{he_ho} = 'not';
					$h_parent->{nb_ref} = 'not';
					$h_parent->{nb_mut} = 'not';
					push(@{$h->{parents}}, $h_parent);
				}
				my (@lParents, @lChildren);
				if (exists $h->{parents}) {
					foreach my $h_parent (@{$h->{parents}}) { $h_parent->{len_cat} = scalar(@{$h->{parents}}); }
					@lParents = @{$h->{parents}};
				}
				if (exists $h->{children}) {
					foreach my $h_child (@{$h->{children}}) { $h_child->{len_cat} = scalar(@{$h->{children}}); }
					@lChildren = @{$h->{children}};
				}
				if (scalar(@lParents) > 0 and scalar(@lChildren) > 0) {
					$lParents[0]->{'_children'} = \@lChildren;
					my $hRes;
					$hRes->{name} = $h->{familly}->{name};
					$hRes->{he_ho} = 'fam';
					$hRes->{nb_ref} = 'fam';
					$hRes->{nb_mut} = 'fam';
					$hRes->{'_children'} = \@lParents;
					$hIds->{$id}->{json_ped}->{$patient->family()} = encode_json $hRes;
				}
				elsif (scalar(@lParents) == 0 and scalar(@lChildren) > 0) {
					my $hRes;
					$hRes->{name} = $h->{familly}->{name};
					$hRes->{he_ho} = 'fam';
					$hRes->{nb_ref} = 'fam';
					$hRes->{nb_mut} = 'fam';
					$hRes->{'_children'} = \@lChildren;
					$hIds->{$id}->{json_ped}->{$patient->family()} = encode_json $hRes;
				}
			}
		}
		
		my ($data) = get_variations::getIds_byCache_onlive($buffer, $project, $chr_cache, $hIds, $user);
		
		export_data::update_deja_vu($project, $data, $user);
		export_data::print_stream($data) if $stream;
		push(@final_data,@$data);
	}
	if ($stream){
		export_data::print_end_stream() ;
		#export_data::print($project,$cgi,\@final_data,"id");
		exit(0);
	}
}

else {
	my $pedigree = {};
	foreach my $patient (@{$project->getPatients()}) {
		$pedigree->{$patient->name()}->{sex} = $patient->sex();
		$pedigree->{$patient->name()}->{status} = $patient->status();
		if ($patient->getFamily()) {
			$pedigree->{$patient->name()}->{fam} = $patient->getFamily->name();
			$pedigree->{$patient->name()}->{father} = $patient->getFamily->getFather->name() if ($patient->getFamily->getFather());
			$pedigree->{$patient->name()}->{mother} = $patient->getFamily->getMother->name() if ($patient->getFamily->getMother());
		}
		if ($patient->isFather())    { $pedigree->{$patient->name()}->{type} = 'F'; }
		elsif ($patient->isMother()) { $pedigree->{$patient->name()}->{type} = 'M'; }
		elsif ($patient->sex() == 1) { $pedigree->{$patient->name()}->{type} = 'C1'; }
		else                         { $pedigree->{$patient->name()}->{type} = 'C2'; }
	}
	my $temp_type_name = $type_name;
	my @tchr = split(" ",$cgi->param('chromosome'));
	foreach my $nchr (@tchr) {
		my $chr = $project->getChromosome($nchr);
		die() unless $chr;
		
		my @lRegions;
		my $vector_resume_ids = Bit::Vector->new_Enum($chr->getVariantsVector->Size(), join(',', @{$select_ids}));
		#TODO: supress regionHO?
		foreach my $patient (@{$chr->getPatients()}) {
			my $vector_region_ho = $patient->getRegionHo($chr, 25)->Clone();
			my @lReg = split(',', $vector_region_ho->to_Enum());
			foreach my $region (@lReg) {
				my $vector_region = Bit::Vector->new_Enum($chr->getVariantsVector->Size(), $region);
				$vector_region->Intersection($vector_region, $vector_resume_ids);
				next if ($vector_region->is_empty());
				my ($v_id_start, $v_id_end) = split('-', $region);
				my $start = $chr->getVarObject($v_id_start)->start();
				my $end   = $chr->getVarObject($v_id_end)->end();
				my $region_id = $patient->name().';chr'.$chr->id().':'.$start.'-'.$end;
				my $length = $end - $start + 1;
				my $hThisRegion;
				$hThisRegion->{id} = $region_id;
				$hThisRegion->{patient_name} = $patient->name();
				$hThisRegion->{start} = $start;
				$hThisRegion->{end} = $end;
				$hThisRegion->{length} = $length;
				push (@lRegions, $hThisRegion);
			}
		}
		purge($chr);
		
		
		my $data;
		#($data) = get_variations_new::getIds($buffer,$project,$chr,$temp_type_name,$select_ids);
		($data) = get_variations::getIds_onlive($buffer,$project,$chr,$select_ids,\@lRegions);
		update_pedigree($project,$data,$pedigree);
		export_data::update_deja_vu($project,$data,$user);
		#update_deja_vu($project,$data);
		update_genes_names($project,$data);
		push(@final_data,@$data);
	}
}
my $type_label = "id";
if ($type_name eq "patients"){
	by_patients($project,\@final_data);
}
if ( $cgi->param('xls') == 1 ) {
	export_data::variation_report_xls($project,\@final_data);
	exit(0);
}
else {
	export_data::print($project,$cgi,\@final_data,$type_label);
}
exit(0);


sub update_pedigree {
	my ($project,$data,$pedigree) = @_;
	foreach my $d (@$data) {
		
		foreach my $p (@{$d->{patient_name}}){
			my $ped = $pedigree->{$p};
			push(@{$d->{pedigree_type}},$ped->{type}.$ped->{status});
			push(@{$d->{pedigree_status}},$ped->{status});
			push(@{$d->{pedigree_fam}},$ped->{fam});
			#$d->{pedigree}->{status}->{$p} = $ped->{status};
			
		}
	}

}


sub update_genes_names {
	my ($project,$data) = @_;
	foreach my $d (@$data) {
		my $z=0;
		my %name2;
			foreach my $c (@{$d->{tab_consequences}}){
				my $a = $c->{gene};
					$a =~ /(.*) \((.*)\)/;
				my $ensg = $1;
				$name2{$1} = $2;
			} 
			$z++;
		$d->{external_names} = join(";",map{$name2{$_}} @{$d->{genes}});
		$d->{genes_name} = join(";",@{$d->{genes}});
	}
}

sub by_patients {
	my ($project,$data) = @_;
	my %patients;
	foreach my $d (@$data) {
		map {$patients{$_} ++}  @{ $d->{patient_name} };
	}
	
	my $out;
	foreach my $name (sort{$a cmp $b} keys %patients){
		my $item;
		$item->{name} = $name;
		push(@$out,$item);
	}
	export_data::print($project,$cgi,$out);
	exit(0);
}


sub update_valid {
	my ($project, $data ,$user) = @_;
	
	my $captures = $project->getCaptures();
	
	my $capture = @$captures[-1];
	warn $capture->validation_db;
	return unless $capture->validation_db;
	my $vquery = validationQuery->new(dbh=>$project->buffer->dbh,capture_name=>$capture->validation_db);

	return unless $vquery->exists_db();
	
	foreach my $d (@$data) {
		my $valid =0;
		my $invalid = 0;
		my $notseq=0;
		my $ho = 0;
		my $he = 0 ;
		my $globalValid = 0;
		my $nbpatients = scalar (@{ $d->{patient_name} });
		my $id;
		
		 my $validation_vid = $vquery->getVariationByGenBoId(id=>$d->{id});
		 
		 warn $validation_vid." ::  ".$d->{id};
		if ($validation_vid){
			foreach my $patient_name (@{ $d->{patient_name} }){
				my ($value_valid,$validation_sanger) = $vquery->getValidations(id=>$validation_vid,project=>$project->name(),sample=>$patient_name);
				$value_valid = $validation_sanger if $validation_sanger;
				$value_valid += 100 if $validation_sanger;
				$d->{ "valid!" .$patient_name} = $value_valid if $value_valid;
				$d->{valid} = "1" if $value_valid;
			}
		}
	}

}

sub purge {
	my $chr = shift;
	$project->{objects}->{proteins}    = {};
	$project->{objects}->{genes}       = {};
	$project->{objects}->{deletions}   = {};
	$project->{objects}->{insertions}  = {};
	$project->{objects}->{variations}  = {};
	$project->{objects}->{exons}       = {};
	$project->{objects}->{transcripts} = {};
	$chr->available_genes_ids(1);
	$chr->hash_filters_deleted(1);
	$chr->hash_freeze_file_genes(1);
	$chr->hash_freeze_file_all_genes(1);
	$chr->hash_freeze_file_patients(1);
	$chr->hash_filters_keeped(1);
	$chr->genes_object(1);
	delete $project->{objects}->{chromosomes}->{$chr->id()}->{fastGenes};
	delete $project->{objects}->{chromosomes}->{$chr->id()}->{fastTranscripts};
#	$chr->cache_lmdb_variations->close();
#	$chr->get_lmdb_patients->close();
#	$chr->get_lmdb_categories->close();
#	$chr->get_lmdb_genes->close() if ($chr->get_lmdb_genes());
#	$chr->get_lmdb_variations->close();
}

