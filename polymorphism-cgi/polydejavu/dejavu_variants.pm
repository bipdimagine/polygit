package dejavu_variants;

use strict;
use Moo;
use Bit::Vector;
use Bit::Vector::Overload;
use Set::IntSpan::Fast;
use Data::Dumper;
use List::MoreUtils qw(natatime);
use Compress::Snappy;
use MIME::Base64;
use Storable qw(store retrieve freeze dclone thaw);
use JSON;
use IPC::Open2;
use Text::CSV;
use MCE::Loop;


#use polyviewer_html;

#use html_polygenescout;

my @headers = ("Viewer", "var_name", "trio", "gnomad", "deja_vu");


has fork => (
	is		=> 'rw',
	lazy    => 1,
	default => 1
);

has buffer => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $buffer = new GBuffer;
		return $buffer;
	}
);

has use_phenotype => (
	is		=> 'rw',
	lazy    => 1,
	#default => sub { return 'intellectual disability'; }
	default => sub { return 'heart'; }
);

has xls_export_session => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $xls_export = new xls_export();
		return $xls_export;
	}
);

has project => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		return $self->buffer->newProject( -name => $self->project_name() );
	}
);

has project_name => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $buffer = new GBuffer;
		return $buffer->getRandomProjectName();
	}
);

has hash_projects_ids_names => (
	is		=> 'rw',
	lazy    => 1,
);

has hash_filters_cons => (
	is		=> 'rw',
	lazy    => 1,
);

has is_magic_user => (
	is		=> 'rw',
	lazy    => 1,
);

has user_name => (
	is		=> 'rw',
	lazy    => 1,
);

has pwd => (
	is		=> 'rw',
	lazy    => 1,
);

has min_reads => (
	is		=> 'rw',
	lazy    => 1,
);

has min_ratio => (
	is		=> 'rw',
	lazy    => 1,
);

has min_promoter_ai => (
	is		=> 'rw',
	lazy    => 1,
);

has min_ncboost => (
	is		=> 'rw',
	lazy    => 1,
);

has hash_ncboost_values => (
	is		=> 'rw',
	lazy    => 1,
);

has max_dejavu => (
	is		=> 'rw',
	lazy    => 1,
);

has max_dejavu_ho => (
	is		=> 'rw',
	lazy    => 1,
);

has max_gnomad_ac => (
	is		=> 'rw',
	lazy    => 1,
);

has max_gnomad_ac_ho => (
	is		=> 'rw',
	lazy    => 1,
);

has only_ill_patients => (
	is		=> 'rw',
	lazy    => 1,
);

has only_strict_ill_patients => (
	is		=> 'rw',
	lazy    => 1,
);

has models => (
	is		=> 'rw',
	lazy    => 1,
);

has view_all_projects => (
	is		=> 'rw',
	lazy    => 1,
);

has sql_projects_parquet => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my $sql;
		if ($self->is_magic_user() or $self->view_all_projects()) {
			my $dir_parquets = $self->buffer->dejavu_parquet_dir();
			$sql = "read_parquet('".$dir_parquets."/NGS20*.parquet')";
		}
		else {
			my @lParquets;
			foreach my $file (keys %{$self->hash_users_projects_parquet()}) {
				next if not -e $file;
				push(@lParquets, "'$file'");
			}
			$sql = "read_parquet([".join(', ', @lParquets)."])";
		}
		return $sql;
	}
);

has hash_users_projects_parquet =>  (
	is		=> 'rw',
	lazy    => 1,
);

has hash_users_projects => (
	is		=> 'rw',
	lazy    => 1,
	default => sub {
		my $self = shift;
		my ($h_projects, $h_projectsNames, $h_parquets, @list_hash);
		if (not $self->hash_exclude_projects() and ($self->view_all_projects() or $self->buffer->getQuery->isUserMagic($self->user_name(), $self->pwd()))) {
			$self->{is_magic_user} = 1;
			my $dir_parquets = $self->buffer->dejavu_parquet_dir();
			opendir my $dir, $dir_parquets or die "Cannot open directory: $!";
			my @files = readdir $dir;
			closedir $dir;
			foreach my $file (@files) {
				next if $file eq '.';
				next if $file eq '..';
				my @ltmp = split('\.', $file);
				next if scalar(@ltmp) != 3;
				next if $ltmp[-1] eq 'no_dejavu';
				my $proj_name = $ltmp[0];
				my $proj_id = $ltmp[1];
				$h_projects->{$proj_name}->{description} = '-';
				$h_projects->{$proj_name}->{name} = $proj_name;
				$h_projects->{$proj_name}->{id} = $proj_id;
				my $list_patients = $self->buffer->getQuery->getPatients($proj_id);
				foreach my $h_pat (@$list_patients) {
					$h_projects->{$proj_name}->{patients}->{$h_pat->{name}}->{family} = $h_pat->{family};
					$h_projects->{$proj_name}->{patients}->{$h_pat->{name}}->{name} = $h_pat->{name};
					$h_projects->{$proj_name}->{patients}->{$h_pat->{name}}->{status} = $h_pat->{status};
					$h_projects->{$proj_name}->{patients}->{$h_pat->{name}}->{sex} = $h_pat->{sex};
					$h_projects->{$proj_name}->{patients}->{$h_pat->{name}}->{father} = $h_pat->{father};
					$h_projects->{$proj_name}->{patients}->{$h_pat->{name}}->{mother} = $h_pat->{mother};
					$h_projects->{$proj_name}->{patients}->{$h_pat->{patient_id}} = $h_projects->{$proj_name}->{patients}->{$h_pat->{name}};
				}
				$h_projects->{$proj_id} = $h_projects->{$proj_name};
			}
		}
		else {
			$self->{is_magic_user} = undef;
			@list_hash = @{$self->buffer->getQuery()->getProjectListForUser($self->user_name(), $self->pwd())};
			foreach my $hash (@list_hash) {
				my $proj_name = $hash->{name};
				next unless ($proj_name =~ /NGS20/);
				next if $self->hash_exclude_projects() and exists $self->hash_exclude_projects->{$proj_name};
				next if $self->hash_only_project() and not exists $self->hash_only_project->{$proj_name};
				$h_projects->{$proj_name}->{description} = $hash->{description};
				$h_projects->{$proj_name}->{name} = $proj_name;
				$h_projects->{$proj_name}->{id} = $hash->{id};
				my $list_patients = $self->buffer->getQuery->getPatients($hash->{id});
				foreach my $h_pat (@$list_patients) {
					$h_projects->{$proj_name}->{patients}->{$h_pat->{name}}->{family} = $h_pat->{family};
					$h_projects->{$proj_name}->{patients}->{$h_pat->{name}}->{name} = $h_pat->{name};
					$h_projects->{$proj_name}->{patients}->{$h_pat->{name}}->{status} = $h_pat->{status};
					$h_projects->{$proj_name}->{patients}->{$h_pat->{name}}->{sex} = $h_pat->{sex};
					$h_projects->{$proj_name}->{patients}->{$h_pat->{name}}->{father} = $h_pat->{father};
					$h_projects->{$proj_name}->{patients}->{$h_pat->{name}}->{mother} = $h_pat->{mother};
					$h_projects->{$proj_name}->{patients}->{$h_pat->{patient_id}} = $h_projects->{$proj_name}->{patients}->{$h_pat->{name}};
				}
				
				if ($self->hash_only_project()) {
					my $only_patient = $self->hash_only_project->{$proj_name};
					if ($only_patient) {
						my $fam_name = $h_projects->{$proj_name}->{patients}->{$only_patient}->{family};
						$self->{hash_only_patients}->{$proj_name}->{$only_patient} = undef;
						foreach my $pat (keys %{$h_projects->{$proj_name}->{patients}}) {
							my $this_fam_name = $h_projects->{$proj_name}->{patients}->{$pat}->{family};
							$self->{hash_only_patients}->{$proj_name}->{$pat} = undef if ($this_fam_name eq $fam_name);
						}
					}
				}
				
				$h_projects->{$hash->{id}} = $h_projects->{$proj_name};
				
				my $parquet = $self->buffer->dejavu_parquet_dir().'/'.$proj_name.'.'.$hash->{id}.'.parquet';
				$h_parquets->{$parquet} = undef if -e $parquet;
			}
		}
		$self->{hash_users_projects_parquet} = $h_parquets;
		return $h_projects;
	}
);

has hash_only_project =>  (
	is		=> 'rw',
	lazy    => 1,
);

has hash_only_patients =>  (
	is		=> 'rw',
	lazy    => 1,
);

has hash_exclude_projects => (
	is		=> 'rw',
	lazy    => 1,
);

has hash_lift_variants => (
	is		=> 'rw',
	lazy    => 1
);

has alert_too_much_results  => (
	is		=> 'rw',
	lazy    => 1,
	default => undef, 
);

has alert_ncboost_min_cadd_25  => (
	is		=> 'rw',
	lazy    => 1,
	default => undef, 
);

has only_chromosome  => (
	is		=> 'rw',
	lazy    => 1,
);

has only_genes  => (
	is		=> 'rw',
	lazy    => 1,
);

has hash_only_genes  => (
	is		=> 'rw',
	lazy    => 1,
);

has hash_variant_pathogenic  => (
	is		=> 'rw',
	lazy    => 1,
);


sub check_variants_from_gene {
	my ($self, $h_dv_rocks_ids) = @_;
	print '.checkvar.!';
	print '!';
	my $ii = 0;
	my ($h_dv_var_ids, @lVarIds, @lVar, $h_dv_var_ids_erros_lift);
	my ($total, $total_pass);
	$self->buffer->dbh_deconnect();
	$self->project->disconnect();
	$self->project->getChromosomes();

	my $fork = $self->fork();
	my $h_cons = $self->hash_filters_cons();
	my $project_name = $self->project_name();
	my $h_dv = $h_dv_rocks_ids;   # alias simple
	
	my ($h_var_pos, $hGenes, $h_var_ids);
	my $nb_var = 0;
	foreach my $chr_id (keys %$h_dv) {
		print '.chr'.$chr_id.'.';
        my $chr = $self->project->getChromosome($chr_id);
		my $nodv = $chr->rocks_dejavu();
		my $fork2 = $fork;
		my $nb_part = 0;
		my $iter = natatime(100000, keys %{$h_dv->{$chr_id}});
		while( my @tmp = $iter->() ){
			my $nb_var = scalar (@tmp);
			my $chunk_size = int($nb_var/$fork2)+1;
			$nb_part++;
			print '.part.'.$nb_part.'.';
			MCE::Loop->init(
			    max_workers => $fork2,
			    chunk_size  => $chunk_size,
			    gather      => sub {
			        my ($data) = @_;
			        return unless $data;
			        print '|';
			        if ($data->{var_pos}) {
			        	#$hres->{var_pos}->{$chr_id}->{$start}->{$var_allele}->{var_id} = $var_id;
			            foreach my $chr_id (keys %{$data->{var_pos}}) {
			            	foreach my $start (keys %{$data->{var_pos}->{$chr_id}}) {
			            		foreach my $var_allele (keys %{$data->{var_pos}->{$chr_id}->{$start}}) {
			            			$h_var_pos->{$chr_id}->{$start}->{$var_allele} = $data->{var_pos}->{$chr_id}->{$start}->{$var_allele};
#			            			if ($self->is_magic_user()) {
#			            				my $var_id = $data->{var_pos}->{$chr_id}->{$start}->{$var_allele}->{var_id};
#			            				$h_var_ids->{$var_id} = undef;
#			            			}
			            			$nb_var++;
			            		}
			            	}
			            }
			        }
			        if ($data->{genes}) {
			            foreach my $gene_id (keys %{$data->{genes}}) {
			            	foreach my $var_id (keys %{$data->{genes}->{$gene_id}}) {
			                	$hGenes->{$gene_id}->{$var_id} = $data->{genes}->{$gene_id}->{$var_id};
			            	}
			            }
			        }
			        if ($data->{hash_projects_ids_names}) {
			        	foreach my $proj_id (keys %{$data->{hash_projects_ids_names}}) {
							$self->{hash_projects_ids_names}->{$proj_id} = $data->{hash_projects_ids_names}->{$proj_id};
			        	}
			        }
			    }
			);
			
			my $worker = sub {
			    my ($mce, $chunk_ref, $chunk_id) = @_;
			    my $b = GBuffer->new();
			    my $p = $b->newProject(-name => $project_name);
			    my $hres = {
			        variants => {},
			        genes    => {}
			    };
		        my $ii = 0;
		        foreach my $rocks_id (@$chunk_ref) {
		        	$ii++;
		        	print '.' if ($ii % 5000 == 0);
		            my $var_id = $chr->transform_rocksid_to_varid($rocks_id);
		            my $var    = $p->_newVariant($var_id);
		            
		            #TODO: HERE
#		            if (not $var->isVariation()) {
#		            	my @ltmp = split('_', $var_id);
#		            	my @ltmp2 = split('-', $ltmp[-1]);
#		            	if (lc($ltmp2[0]) eq 'inv') {
#							bless $var , 'GenBoInversion';
#							$var->{start} = $ltmp[1]; 
#							$var->{length} = $ltmp2[1];
##							delete $var->{annotation};
##							$var->annotation();
##							warn "\n\n";
##							warn Dumper $var->{annotation};
#		            	}
#		            }
					
					my $is_pathognic;
					$is_pathognic = 1 if $self->hash_variant_pathogenic() and exists $self->{hash_variant_pathogenic}->{$rocks_id};
					next if not $is_pathognic and $var->other_patients() > $self->max_dejavu();
					next if not $is_pathognic and $var->other_patients_ho() > $self->max_dejavu_ho();
					
#					my $h_dv = $nodv->dejavu($rocks_id);
#					next if not $h_dv;
                	my $is_ok = 0;
		            if ($chr->intergenic_intspan->contains($var->start())) {
		                if (exists $h_cons->{intergenic}) {
		                    $hres->{genes}->{intergenic}->{$var_id} = 1;
		                    $is_ok = 1;
		                }
		            }
		            else {
		                my $genes = $var->getGenes();
		                foreach my $gene (@$genes) {
		                	my $var_annot;
		                	eval {
			                    $var_annot = $var->variationTypeInterface($gene);
#			                    warn ref($var).' - '.$var->id.' -> '.$var_annot;
			                    foreach my $annot (split(',', $var_annot)) {
			                        $annot =~ s/ /_/g;
			                        if (exists $h_cons->{lc($annot)}) {
			                            $is_ok = 1;
			                            last;
			                        }
			                    }
		                	};
		                	if ($@) { $is_ok = 0; }
		                    next unless $is_ok;
		                    $hres->{genes}->{$gene->id}->{$var_id} = 1;
		                }
		            }
		            next if not $is_ok;
		            
					my $start = $var->start();
					my $var_allele = $var->var_allele();
					my $ref_allele = $var->ref_allele();
					my $gnomad_id = $var->gnomad_id();
				
					my $found;
					foreach my $project_id (keys %{$h_dv_rocks_ids->{$chr_id}->{$rocks_id}}) {
						next if not exists $self->hash_users_projects->{all} and not exists $self->hash_users_projects->{$project_id};
						my $project_name = $self->hash_users_projects->{$project_id}->{name};
						my $project_type_cache = $self->{hash_projects_ids_names}->{$project_id}->{type};
						$hres->{hash_projects_ids_names}->{$project_id}->{name} = $project_name;
						$hres->{hash_projects_ids_names}->{$project_id}->{type} = $project_type_cache;
						$found++;
					}
					if ($found) {
						if (not $var->isVariation()) {
							$var_allele = 'indel';
							my @ltmp = split('-', $gnomad_id);
							if (lc($ltmp[-2]) eq 'cnv') {
								$gnomad_id = $var_id;
#								$var_allele = 'cnv';
							}
							else { $start = $ltmp[1]; }
						}
						$hres->{var_pos}->{$chr_id}->{$start}->{$var_allele}->{var_id} = $var_id;
						$hres->{var_pos}->{$chr_id}->{$start}->{$var_allele}->{ref_all} = $ref_allele;
						$hres->{var_pos}->{$chr_id}->{$start}->{$var_allele}->{gnomad_id} = $gnomad_id;
					}
		            $var = undef;
		        }
			    $p = undef;
			    $b = undef;
			    MCE->gather($hres);
			};
			MCE::Loop->run($worker, [ @tmp ]);
			MCE::Loop->finish();
		}
		$nodv->close();
	}
	$self->buffer->dbh_deconnect();
	$self->project->disconnect();
	
	#TODO: detail des patients;
	my ($h_projects_patients, $h_gnomadid);
#	if ($self->is_magic_user()) {
#		$h_projects_patients = $h_var_ids;
#	}
#	else {
		print '.before_duck_projects.nbvar.';
		($h_projects_patients, $h_gnomadid) = $self->get_from_duckdb_project_patients_infos_global($h_var_pos);
		print '.after_duck_projects.';
#	}
	
	my $hVariants_ok;
	my $h_models = $self->models();
	my $min_ratio = $self->min_ratio();
	my $fork2 = 1;
	
	my ($B, $P, $H_intergenic);
	MCE::Loop->init(
	   max_workers => $fork2,
	   chunk_size => 'auto',
	   user_begin => sub {
	       $B = new GBuffer;
	       $P = $B->newProject(-name => $project_name);
	       foreach my $chr (@{$P->getChromosomes()}) {
				$H_intergenic->{$chr->id()} = $chr->intergenic_intspan();
	       }
	   },
	   gather => sub {
	        my ($data) = @_;
	        print '|';
			my $iii = 0;
			if (exists $data->{lift}) {
				foreach my $gid (keys %{$data->{lift}}) {
					$self->{hash_lift_variants}->{$gid} = $data->{lift}->{$gid};
				}
			}
			print '.';
			if (exists $data->{variants}) {
				foreach my $var_id (keys %{$data->{variants}}) {
					$hVariants_ok->{$var_id} = $data->{variants}->{$var_id};
				}
			}
			print '.';
	    }
	);
	mce_loop {
		my ($mce, $chunk_ref, $chunk_id) = @_;
		print '.';
		my ($list, $hres);
		my $p  = $P;
		my $ip = 0;
		my $can_construct;
		my $ii = 0;
		my $h_genes_trans;
		
		foreach my $var_id (@$chunk_ref) {
			my $var = $p->_newVariant($var_id);
#			if ($self->is_magic_user()) {
#				$can_construct = 1;
#			}
#			else {	
				foreach my $project_name (keys %{$h_projects_patients->{$var_id}}) {
					my ($h_pat_done, $h_pat_filtred);
					next if ($self->hash_only_project() and not exists $self->{hash_only_project}->{$project_name});
#					next if (not $self->is_magic_user() and $self->hash_only_project() and not exists $self->{hash_only_project}->{$project_name});
					foreach my $patient_name (keys %{$h_projects_patients->{$var_id}->{$project_name}}) {
						next if ($self->hash_only_patients() and not exists $self->{hash_only_patients}->{$project_name}->{$patient_name});
#						next if (not $self->is_magic_user() and $self->hash_only_patients() and not exists $self->{hash_only_patients}->{$project_name}->{$patient_name});
						push(@{$hres->{variants}->{$var_id}->{polyviewer_html}}, $h_projects_patients->{$var_id}->{$project_name}->{$patient_name}->{print_html});
						my $hh;
						$hh->{project_name} = $project_name;
						$hh->{patient_name} = $patient_name;
						$hh->{heho} = $h_projects_patients->{$var_id}->{$project_name}->{$patient_name}->{heho};
						$hh->{ratio} = $h_projects_patients->{$var_id}->{$project_name}->{$patient_name}->{ratio};
						$hh->{dp} = $h_projects_patients->{$var_id}->{$project_name}->{$patient_name}->{dp};
						$hh->{model} = $h_projects_patients->{$var_id}->{$project_name}->{$patient_name}->{model};
						if ($h_models) {
							next if ($self->hash_only_patients() and not exists $h_models->{$hh->{model}});
#							next if (not $self->is_magic_user() and $self->hash_only_patients() and not exists $h_models->{$hh->{model}});
							$h_pat_filtred->{$patient_name}->{$hh->{model}} = 1 if (not exists $h_models->{$hh->{model}});
						}
						if ($min_ratio) {
							my $this_ratio = $h_projects_patients->{$var_id}->{$project_name}->{$patient_name}->{ratio_value};
							$h_pat_filtred->{$patient_name}->{'ratio_'.$this_ratio} = 1 if ($this_ratio ne '?' and $this_ratio < $min_ratio);
						}
						$h_pat_done->{$patient_name} = $hh;
					}
					if (keys %$h_pat_done > keys %$h_pat_filtred) {
						foreach my $patient_name (keys %{$h_pat_done}) {
							push(@{$hres->{variants}->{$var_id}->{polyviewer_html_details_proj_pat}}, $h_pat_done->{$patient_name});
							$can_construct = 1;
						}
					}
				}
				foreach my $gid (keys %$h_gnomadid) {
					$hres->{lift}->{$gid} = $h_gnomadid->{$gid};
				}
#			}
			
			if ($can_construct) {
				my $vp = PolyviewerVariant->new();
				$vp->setLmdbVariant($var);
				$vp->{hgenes} = {};
				$vp->{genes_id} = [];
				my $code = 0;
				my @ltmp = split('_', $var->id);
				if ($H_intergenic->{$ltmp[0]}->contains($var->start())) {
					my $h = $vp->set_intergenic($var);
					$h->{code} = $code;
					$vp->{hgenes}->{intergenic} = $h;
					push(@{$vp->{intergenic}},'intergenic');
					$code ++;
				}
				else {
					foreach my $g (@{$var->getGenes}){
						my $h = $vp->set_gene($var,$g);
						$h->{code} = $code;
						$vp->{hgenes}->{$g->id} = $h;
						push(@{$vp->{genes_id}},$g->id);
						$code ++;
					}
				}
				$vp->{hpatients} ={};
				$vp->{patients_id} = [];
				$hres->{variants}->{$var_id}->{polyviewer_variant} = $vp;
				$hres->{variants}->{$var_id}->{can_export} = 1;
			}
		}
		MCE->gather($hres);
	} keys %{$h_projects_patients};	
	MCE::Loop->finish();
	$self->project->disconnect();
	print '._end_MCE_check_.';
	return ($hGenes, $hVariants_ok);
}


sub get_table_project_patients_infos {
	my ($self, $project_name, $duck_he, $duck_patients, $duck_ratios) = @_;
	my $nb_he = $duck_he;
	my ($h_infos_patients, $h_tmp_pat);
	my $is_heho = 'He';
	my $nb_pat = 0;
	foreach my $pat_id (unpack("w*",decode_base64($duck_patients))) {
		$nb_pat++;
		$h_infos_patients->{$nb_pat}->{id} = $pat_id;
		$h_tmp_pat->{$pat_id} = $nb_pat;
		if (int($nb_pat) <= $nb_he) {
			$is_heho = 'He';
		}
		else { $is_heho = 'Ho'; }
		$h_infos_patients->{$nb_pat}->{heho} = $is_heho;
	}
	
	my $i = 0;
	$nb_pat = 1;
	foreach my $info (unpack("w*",decode_base64($duck_ratios))) {
		$i++;
		if ($i == 1) { $h_infos_patients->{$nb_pat}->{dp} = $info; }
		elsif ($i == 2) {
			my $ratio = '?';
			$ratio = ($info / $h_infos_patients->{$nb_pat}->{dp}) * 100 if ($h_infos_patients->{$nb_pat}->{dp});
			my $text = 'Reads:'.$info;
			$h_infos_patients->{$nb_pat}->{nb_reads} = $info;
			if ($ratio and $ratio eq '?') {
				$text .= ', Ratio: ?';
				$h_infos_patients->{$nb_pat}->{ratio_value} = '?';
			}
			elsif ($ratio) {
				$text .= ', Ratio:'.int($ratio);
				$h_infos_patients->{$nb_pat}->{ratio_value} = int($ratio);
			}
			$h_infos_patients->{$nb_pat}->{ratio} = $text;
		}
		elsif ($i == 3) {
    		my $model;
    		if ($info == 1) { $model = 'solo'; }
    		elsif ($info == 2) { $model = 'father'; }
    		elsif ($info == 4) { $model = 'mother'; }
    		elsif ($info == 8) { $model = 'both'; }
    		elsif ($info == 16) { $model = 'is_parent'; }
    		elsif ($info == 32) { $model = 'recessif'; }
    		elsif ($info == 64) { $model = 'dominant'; }
    		elsif ($info == 128) { $model = 'denovo'; }
    		elsif ($info == 256) { $model = 'strict_denovo'; }
    		elsif ($info == 512) { $model = 'error'; }
    		else { $model = 'error2'; }
			$h_infos_patients->{$nb_pat}->{model} = $model;
			$i = 0;
			$nb_pat++;
		}
	}
	my $hres;
	
	my $found_healthy_patient = 0;
	my $found_ill_patient = 0;
	my $ok_nb_reads = 0;
	foreach my $pat_id (keys %{$self->hash_users_projects->{$project_name}->{patients}}) {
		next if not exists $h_tmp_pat->{$pat_id};
		my $h_pat = $self->hash_users_projects->{$project_name}->{patients}->{$pat_id};
		$h_infos_patients->{$h_tmp_pat->{$pat_id}}->{name} = $h_pat->{name};
		$h_infos_patients->{$h_tmp_pat->{$pat_id}}->{sex} = '-';
		$h_infos_patients->{$h_tmp_pat->{$pat_id}}->{sex} = 'male' if $h_pat->{sex} eq '1';
		$h_infos_patients->{$h_tmp_pat->{$pat_id}}->{sex} = 'female' if $h_pat->{sex} eq '2';
		$h_infos_patients->{$h_tmp_pat->{$pat_id}}->{status_txt} = '-';
		if ($h_pat->{status} eq '1') {
			$h_infos_patients->{$h_tmp_pat->{$pat_id}}->{status_txt} = 'healthy';
			$found_healthy_patient++;
			if ($self->min_reads() and not $self->only_ill_patients() and not $self->only_strict_ill_patients()) {
				$ok_nb_reads = 1 if $h_infos_patients->{$h_tmp_pat->{$pat_id}}->{nb_reads} >= $self->min_reads();
			}
		}
		if ($h_pat->{status} eq '2') {
			$h_infos_patients->{$h_tmp_pat->{$pat_id}}->{status_txt} = 'ill';
			$found_ill_patient++;
			if ($self->min_reads()) {
				$ok_nb_reads = 1 if $h_infos_patients->{$h_tmp_pat->{$pat_id}}->{nb_reads} >= $self->min_reads();
			}
		}
		$h_infos_patients->{$h_tmp_pat->{$pat_id}}->{family} = $h_pat->{family};
		$h_infos_patients->{$h_tmp_pat->{$pat_id}}->{description} = $self->hash_users_projects->{$project_name}->{description};
	}
	$ok_nb_reads = 1 if not $self->min_reads();
	
	
	my $ok_ill = 1;
	if ($self->only_ill_patients()) {
		$ok_ill = undef if $found_ill_patient == 0;
	}
	elsif ($self->only_strict_ill_patients()) {
		$ok_ill = undef if $found_healthy_patient > 0;
	}
	if ($ok_ill and $ok_nb_reads) {
		foreach my $id (keys %{$h_infos_patients}) {
			my $pat_name = $h_infos_patients->{$id}->{name};
			next if not $pat_name;
			$hres->{$pat_name} = $h_infos_patients->{$id};
		}
	}
	return $hres;
}

sub get_from_duckdb_project_patients_infos_global {
	my ($self, $h_var_pos) = @_;
	print '.';
	return if not $h_var_pos;
	
	my $fork = $self->fork();
	my $sql_parquets = $self->sql_projects_parquet();
	my ($h_gnomadid, $h_projects_patients);
	
	foreach my $chr_id (sort keys %$h_var_pos) {
		my $sql = qq{
			PRAGMA threads=$fork;
			CREATE TEMP TABLE positions( pos38 INT );
			INSERT INTO positions VALUES
		};
		my $pid = open2(my $out, my $in, "duckdb -csv");
		print $in $sql;
		my $first = 1;
		foreach my $chr_id (sort keys %$h_var_pos) {
			foreach my $pos (sort keys %{$h_var_pos->{$chr_id}}) {
			    my $line = "($pos)";
			    if ($first) {
			        print $in $line;
			        $first = 0;
			    } else {
			        print $in ",\n$line";
			    }
			}
		}
		print $in ";\n";
		print '.temp.done.';
		
		print $in qq{
			WITH filtered AS ( SELECT * FROM $sql_parquets where chr38='$chr_id' )
			SELECT f.project, f.chr38, f.chr19, f.pos38, f.pos19, f.he, f.allele, f.patients, f.dp_ratios
			FROM filtered f JOIN positions p ON f.pos38 = p.pos38;
		};
		close($in);
		
		my $iii = 0;
		my $csv_in = Text::CSV->new({ binary => 1 });
		my $header = $csv_in->getline($out);
		while (my $row = $csv_in->getline($out)) {
		    my ($project_id,$this_chr38,$this_chr19,$this_pos38,$this_pos19,$he,$var_all,$patients,$dp_ratios) = @$row;
		    $var_all = 'indel' if ($var_all ne 'T' and $var_all ne 'G' and $var_all ne 'A' and $var_all ne 'C');
		    next if (not exists $h_var_pos->{$this_chr38}->{$this_pos38}->{$var_all});
		    
		    my $var_id    = $h_var_pos->{$this_chr38}->{$this_pos38}->{$var_all}->{var_id};
		    my $gnomad_id = $h_var_pos->{$this_chr38}->{$this_pos38}->{$var_all}->{gnomad_id};
		    my $ref       = $h_var_pos->{$this_chr38}->{$this_pos38}->{$var_all}->{ref_all};
		    $h_gnomadid->{$gnomad_id}->{chr19} = $this_chr19;
		    $h_gnomadid->{$gnomad_id}->{pos19} = $this_pos19;
		    $h_gnomadid->{$gnomad_id}->{ref_all} = $ref;
		    $h_gnomadid->{$gnomad_id}->{var_all} = $var_all;
	    
		    my $project_name;
		    if (exists $self->{hash_projects_ids_names}->{$project_id}) {
		    	$project_name = $self->{hash_projects_ids_names}->{$project_id}->{name};
		    }
		    else {
		    	$project_name = $self->buffer->getQuery->getProjectNameFromId($project_id);
		    	$self->{hash_projects_ids_names}->{$project_id}->{name} = $project_name;
				my $list_patients = $self->buffer->getQuery->getPatients($project_id);
				foreach my $h_pat (@$list_patients) {
					$self->{hash_projects_ids_names}->{$project_id}->{patients}->{$h_pat->{name}}->{family} = $h_pat->{family};
					$self->{hash_projects_ids_names}->{$project_id}->{patients}->{$h_pat->{name}}->{name} = $h_pat->{name};
					$self->{hash_projects_ids_names}->{$project_id}->{patients}->{$h_pat->{name}}->{status} = $h_pat->{status};
					$self->{hash_projects_ids_names}->{$project_id}->{patients}->{$h_pat->{name}}->{sex} = $h_pat->{sex};
					$self->{hash_projects_ids_names}->{$project_id}->{patients}->{$h_pat->{name}}->{father} = $h_pat->{father};
					$self->{hash_projects_ids_names}->{$project_id}->{patients}->{$h_pat->{name}}->{mother} = $h_pat->{mother};
					$self->{hash_projects_ids_names}->{$project_id}->{patients}->{$h_pat->{patient_id}} = $self->{hash_projects_ids_names}->{$project_id}->{patients}->{$h_pat->{name}};
				}
				$self->{hash_projects_ids_names}->{$project_name} = $self->{hash_projects_ids_names}->{$project_id};
		    }
		    $h_projects_patients->{$var_id}->{$project_name} = $self->get_table_project_patients_infos($project_name, $he, $patients, $dp_ratios);
		    $iii++;
		    print '.' if $iii % 500 == 0;
		}
		print '.';
		close($out);
		waitpid($pid, 0);
	}
	return ($h_projects_patients, $h_gnomadid);
}

sub print_html_gene {
	my ($self, $gene_id, $list_variants, $hVariantsDetails) = @_;
	print '|';
	my $nb_ok = 0;
	my @lVar_ok;
	foreach my $var_id (sort @$list_variants) {
		next if not exists $hVariantsDetails->{$var_id};
#		next if not $self->is_magic_user() and not $hVariantsDetails->{$var_id}->{polyviewer_html_details_proj_pat};
		next if not $hVariantsDetails->{$var_id}->{polyviewer_html_details_proj_pat};
		$nb_ok++;
		push(@lVar_ok, $var_id);
	}
	return if $nb_ok == 0;
	
#	warn "\n\n\nNB VAR: ".scalar(@$list_variants).' -> '.$nb_ok;
	
	my $pr = $self->project;
	my $pat = $self->project->getPatients->[0];
	my $print_html = polyviewer_html->new( project=>$pr, patient=>$pat,header=>\@headers);

	my ($g, $gene, $max_gene_score, $h_var_scores);
	if ($gene_id eq 'intergenic') {
		$g->{id} = $gene_id;
		$g->{uid} = $gene_id;
		$g->{name} = 'Intergenic';
		$g->{external_name} = 'Intergenic';
		$g->{chr_name} = 'Intergenic';
		$g->{nb} = $nb_ok;
		$g->{omim_id} = undef;
		$g->{pLI} = undef;
		$g->{phenotypes} = 'Intergenic';
		$g->{variants} = \@lVar_ok;
		$g->{panels} = undef;
		($max_gene_score, $h_var_scores) = $self->get_score_variant_from_gene_without_patient(\@lVar_ok, $hVariantsDetails);
		$g->{max_score} = $max_gene_score;
	}
	else {
		$gene = $self->project->newGene($gene_id);
		$g->{id} = $gene_id;
		$g->{uid} = $gene_id;
		$g->{name} = $gene->external_name();
		$g->{external_name} = $gene->external_name();
		$g->{chr_name} = $gene->getChromosome->id();
		$g->{nb} = $nb_ok;
		$g->{omim_id} = $gene->omim_id();
		$g->{pLI} = $gene->pLI();
		$g->{phenotypes} = $gene->description().' '.$gene->polyquery_phenotypes();
		$g->{variants} = \@lVar_ok;
		$g->{panels} = $self->buffer->queryPanel()->getPanelsForGeneName($gene->external_name);
		($max_gene_score, $h_var_scores) = $self->get_score_variant_from_gene_without_patient(\@lVar_ok, $hVariantsDetails, $gene);
		$g->{max_score} = $max_gene_score;
	}
	
	next if ($self->only_genes() and not exists $self->{hash_only_genes}->{$g->{name}});
	
	my $panel_id = "panel_".$gene_id;
	$panel_id .= "_he_comp" if $self->only_genes() and $self->hash_only_project() and $self->hash_only_patients();
	my ($out, $h_phenos);
	my $bg_color = $print_html->bgcolor;
	$out .= qq{<div class="panel-heading panel-face panel-grey"	style="$bg_color;height:43px;padding:10px;border:0px;width:100%;">};
	my $out_g = update_variant_editor::panel_gene($g, $panel_id, $self->project->name);
	$out_g =~ s/>CNV/ hidden>CNV/;
	$out .= $out_g;
	my $table_id = 'table_'.$panel_id;
	$out .= qq{</div>};
	$out .= qq{<div style="height:3px;" loading="lazy"></div>};
	$out .= qq{<div loading="lazy" class="panel-body panel-collapse collapse" style="font-size: 09px;font-family:Verdana;" id="$panel_id" loading="lazy">};
	
	if ($self->only_genes() and $self->hash_only_project() and $self->hash_only_patients() and $nb_ok >= 5) {
		$out .= qq{<table loading="lazy" id='$table_id' data-filter-control='true' data-virtual-scroll="true" data-toggle="table" data-show-extended-pagination="false" data-cache="false" data-pagination-loop="true" data-pagination-v-align="both" data-pagination-h-align="left" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="5" data-page-list="[5, 10, 20]" data-resizable='true' class='table table-striped table-borderless table-bootstraptable' style='font-size:9px;'>};
	}
	elsif ($nb_ok >= 20) {
		$out .= qq{<table loading="lazy" id='$table_id' data-filter-control='true' data-virtual-scroll="true" data-toggle="table" data-show-extended-pagination="false" data-cache="false" data-pagination-loop="true" data-pagination-v-align="top" data-pagination-h-align="left" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="10" data-page-list="[5, 10, 20]" data-resizable='true' class='table table-striped table-borderless table-bootstraptable' style='font-size:9px;'>};
	}
	else {
		$out .= qq{<table class="table table-striped table-borderless" style="vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;line-height: 25px;min-height: 25px;height: 25px;box-shadow: 3px 3px 5px #555;">};
	}
	
	
	$out .= "<thead><tr style='background-color:aliceblue;color:black'>";
	$out .= qq{<th data-field="1"></th>};
	$out .= qq{<th data-field="2" data-sortable="true" data-filter-control="input" data-filter-control-placeholder="">Position</th>};
	$out .= qq{<th data-field="3" data-filter-control="input" data-filter-control-placeholder="">Project / Patient</th>};
	$out .= qq{<th data-field="4">gnomAD</th>};
	$out .= qq{<th data-field="5">DejaVu</th>};
	$out .= qq{<th data-field="6" data-filter-control="input" data-filter-control-placeholder="">Clinvar</th>};
	$out .= qq{<th data-field="7" data-filter-control="input" data-filter-control-placeholder="">Annotation</th>};
	$out .= "</tr></thead>";
	$out .= "<tbody>";
	
	#my $h_var_done;
	my $color_validation = "grey";
	my $i = 0;
	foreach my $var_id (sort @lVar_ok) {
		next if not exists $hVariantsDetails->{$var_id};
#		next if not $self->is_magic_user() and not $hVariantsDetails->{$var_id}->{polyviewer_html_details_proj_pat};
		next if not $hVariantsDetails->{$var_id}->{polyviewer_html_details_proj_pat};
		$i++;
		if ($i == 50) {
			$i = 0;
			print '.';
		}
		my $opacity = undef;
		my $var = $self->project->_newVariant($var_id);
		my $polyviewer_variant = $hVariantsDetails->{$var_id}->{polyviewer_variant};
		$polyviewer_variant->{gene} = $g;
		$polyviewer_variant->{ncboost_value} = $hVariantsDetails->{$var_id}->{ncboost};
		$polyviewer_variant->{transcripts} = $polyviewer_variant->{hgenes}->{$g->{id}}->{tr};
		
#		my @list_polyviewer_h_details;
#		if ($self->is_magic_user()) { @list_polyviewer_h_details = []; }
#		else { @list_polyviewer_h_details = @{$hVariantsDetails->{$var_id}->{polyviewer_html_details_proj_pat}}; }
		my @list_polyviewer_h_details = @{$hVariantsDetails->{$var_id}->{polyviewer_html_details_proj_pat}};
		
		$print_html->variant($polyviewer_variant);
		my ($this_out, $this_h_pheno) = $self->print_line_variant_all_patients($polyviewer_variant, $print_html, \@list_polyviewer_h_details, $opacity);
		#$h_var_done->{$var_id};
		$out .= $this_out;
		foreach my $pheno_name (keys %$this_h_pheno) { $h_phenos->{$pheno_name} = $this_h_pheno->{$pheno_name}; }
	}
	
	$out .= "</tbody>";
	$out .= qq{</table>};
	$out .= qq{</div>};
	
	$gene = undef;
	return ($out, $h_phenos, $max_gene_score);
}


sub print_line_variant_all_patients {
	my ($self,$polyviewer_variant,$print_html,$list_h_details,$opacity) = @_;
	my $style = { style => "max-height:180px;" };
	$style = { style => "background-color: #DAEEED;opacity:0.5;max-height:180px;" }  if $opacity;
	my ($out, $h_phenos);
	my $hpatients;
	my $i = 0;
	my $h_proj_pat_list_print_html;
#	if (not $self->is_magic_user()) {
		foreach my $h_proj_pat (@{$list_h_details}) {
			my $b = new GBuffer;
			my $pr = $b->newProject(-name => $h_proj_pat->{project_name});
			$pr->getPatients();
			$pr->getFamilies();
			$pr->project_root_path();
			if (not exists $self->{hash_users_projects}->{$h_proj_pat->{project_name}}->{phenotypes}) {
				$self->{hash_users_projects}->{$h_proj_pat->{project_name}}->{phenotypes} = join(', ', sort @{$pr->phenotypes()});
				foreach my $pheno (@{$pr->phenotypes()}) {
					$self->{hash_users_projects}->{$h_proj_pat->{project_name}}->{phenotypes_tags} .=  ' phenotype '.$pheno;
					$h_phenos->{$pheno} = ' phenotype '.$pheno;
				}
			}
			my $pat = $pr->getPatient($h_proj_pat->{patient_name});
			my $fam_name = $pat->getFamily->name();
			foreach my $this_pat (@{$pat->getFamily->getPatients()}) { $this_pat->alignmentMethods(); }
			my $this_print_html = polyviewer_html->new( project=>$pr, patient=>$pat,header=>\@headers, bgcolor=>"background-color:#607D8B" );
			$this_print_html->variant($polyviewer_variant);
			$this_print_html->variant->{patients_calling}->{$pat->id()}->{gt} = $h_proj_pat->{heho};
			$this_print_html->variant->{patients_calling}->{$pat->id()}->{pc} = $h_proj_pat->{ratio};
			$this_print_html->variant->{patients_calling}->{$pat->id()}->{dp} = 'DP:'.$h_proj_pat->{dp};
			$this_print_html->variant->{patients_calling}->{$pat->id()}->{model} = $h_proj_pat->{model};
			if ($pat->isParent) {
				push (@{$h_proj_pat_list_print_html->{$h_proj_pat->{project_name}}->{$fam_name}->{parents}}, $this_print_html);
			}
			else {
				push (@{$h_proj_pat_list_print_html->{$h_proj_pat->{project_name}}->{$fam_name}->{children}}, $this_print_html);
			}
		}
#	}
	my $cgi = $print_html->cgi;
	my $icon = qq{<img width="32" height="32" src="https://img.icons8.com/external-gliphyline-royyan-wijaya/32/external-laptop-laptop-collection-glyphyline-gliphyline-royyan-wijaya-15.png" alt="external-laptop-laptop-collection-glyphyline-gliphyline-royyan-wijaya-15"/>};
	$icon   = qq{<img width="24" height="24" src="https://img.icons8.com/external-tal-revivo-filled-tal-revivo/24/external-live-preview-of-a-smart-class-education-school-filled-tal-revivo.png" alt="external-live-preview-of-a-smart-class-education-school-filled-tal-revivo"/>};
	my $dropdown = qq{
		<div class="dropdown">
		<button class="btn btn-primary btn-xs dropdown-toggle " type="button" id="dropdownMenuButton" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false" style="font-size:10px;background-color:#C67FAE">
		$icon
		</button>
		<div class="dropdown-menu" aria-labelledby="dropdownMenuButton" style="font-size:12px;background-color:beige;color:black">
	};
	$dropdown .= "<li>".$print_html->mobidetails()."</li>";
  	$dropdown .= "<li>".$print_html->gnomadurl()."</li>";
	$dropdown .= "<li>".$print_html->alamuturl()."</li>";
	$dropdown .= "<li>".$print_html->varsome()."</li>";
	$dropdown .= qq{</div></div>};
	$out .= $cgi->td( $style, $dropdown );
	$out .= "\n";
	
	my $var_text = $print_html->var_name();
#	if (not $self->is_magic_user()) {
		$var_text .= "<br><i><b>HG19: ".$self->{hash_lift_variants}->{$polyviewer_variant->gnomad_id()}->{chr19}.'-'.$self->{hash_lift_variants}->{$polyviewer_variant->gnomad_id()}->{pos19}.'</b></i>';
#	}
	$out .= $cgi->td($style, $var_text);
	$out .= "\n";
	
	my (@l_html_calling, $h_fam_done);
	foreach my $proj_name (sort keys %{$h_proj_pat_list_print_html}) {
		foreach my $fam_name (sort keys %{$h_proj_pat_list_print_html->{$proj_name}}) {
			my $list_print_html;
			if (exists $h_proj_pat_list_print_html->{$proj_name}->{$fam_name}->{children}) {
				push(@$list_print_html, $h_proj_pat_list_print_html->{$proj_name}->{$fam_name}->{children}->[0]);
			}
			else {
				@$list_print_html = @{$h_proj_pat_list_print_html->{$proj_name}->{$fam_name}->{parents}};
			}
			foreach my $this_print_html (@$list_print_html) {
				my $html_pat = $this_print_html->calling();
				
				my (@l_names, @l_bam);
				foreach my $pat (@{$this_print_html->patient->getFamily->getPatients()}) {
					push(@l_names, $pat->name());
					push(@l_bam, $pat->bamUrl());
				}
				my $pnames = join(';', @l_names);
				my $f = join(';', @l_bam);
				my $gn = $this_print_html->patient->getProject->getVersion();
				my $locus;
				if ($gn =~ /HG19/) { $locus = $self->{hash_lift_variants}->{$polyviewer_variant->gnomad_id()}->{chr19}.':'.$self->{hash_lift_variants}->{$polyviewer_variant->gnomad_id()}->{pos19}.'-'.$self->{hash_lift_variants}->{$polyviewer_variant->gnomad_id()}->{pos19}; }
				else { $locus = $polyviewer_variant->locus(); }
				my $igv_b = qq{<button class='igvIcon2 rounded-circle p-2 lh-1' onclick='launch_web_igv_js("$proj_name","$pnames","$f","$locus","/","$gn")' style="color:black;padding-top:2px;"></button>};
				my $project_phenotypes = '';
				if ($self->{hash_users_projects}->{$proj_name}->{phenotypes}) {
					$project_phenotypes = qq{<span hidden>}.$self->{hash_users_projects}->{$proj_name}->{phenotypes_tags}.qq{</span>};
					$project_phenotypes .= qq{<span style='color:green;'><i>}.$self->{hash_users_projects}->{$proj_name}->{phenotypes}.qq{</i></span><br>};
				}
				my $project_description = $this_print_html->patient->getProject->name().', Description: '.$this_print_html->patient->getProject->description();
				push(@l_html_calling, "<tr style='padding:6px;'><td style='padding-right:10px;width:120px;'><center>".$project_phenotypes."<button class='btn btn-outline-secondary d-inline-flex align-items-left' onClick='alert(\"$project_description\")' style='padding-bottom:2px;font-size:10px;'><b>".$this_print_html->patient->getProject->name().'</b></button><br>'.$igv_b.'&nbsp;&nbsp;&nbsp;<b>'.$this_print_html->patient->getFamily->name()."</b></td></center></td><td><nobr>".$html_pat."</nobr></td></tr>");
				push(@l_html_calling, "<tr style='padding:6px;'><td><br></td><td><br></td></tr>");
			}
		}
	}
	
	my $out_calling = qq{<center><div style='width:100%;max-height:200px;overflow-y: scroll;'><table style='width:95%;max-height:200px;'>};
	$out_calling .= join('', @l_html_calling);
	$out_calling .= qq{</table></div></center>};
	
	$out .= $cgi->td( $style, $out_calling) ;
	$out .= "\n";
	$out .= $cgi->td( $style, $print_html->gnomad() );
	$out .= "\n";
	$out .= $cgi->td( $style, $print_html->dejavu() );
	$out .= "\n";
	$out .= $cgi->td( $style, $print_html->validations );
#	$out .= $cgi->td( $style, '' );
	$out .= "\n";
	$out .= $cgi->td( $style, $print_html->transcripts() );
	$out .= $cgi->end_Tr();
	$out .= "\n";
	return ($out, $h_phenos);
}

sub get_score_variant_from_gene_without_patient {
	my ($self, $list_variants, $hVariantsDetails, $gene) = @_;
	my $max_gene_score = -999;
	my $h_var_scores;
	foreach my $var_id (@$list_variants) {
		next if not exists $hVariantsDetails->{$var_id};
		next if not $hVariantsDetails->{$var_id}->{polyviewer_html_details_proj_pat};# and not $self->is_magic_user();
		my $max_score = -99;
		my $var = $self->project->_newVariant($var_id);
		
		my $max_scaled_score = -999;
		if ($gene) {
			foreach my $tr (@{$gene->getTranscripts()}) {
				my $this_score = $var->score_transcript_consequence($tr, undef);
				$this_score +=  $var->score_frequence_public();
				$this_score += $var->scaled_score_cadd();
				my $scaled_score = 1;
				$scaled_score = 2 if $this_score > 50;
				$scaled_score = 3 if $this_score >= 150;
				$scaled_score = 4 if $this_score >= 200;
				
				my $val = $var->score_validations($tr);
				if ($val) {
					$scaled_score = $val->{validation};
				}
				if ($var->is_clinvar_pathogenic_for_gene($gene)){
					my $gac  = $var->getGnomadAC ;
					$gac = 0 unless $gac;
					$scaled_score ++;
					$scaled_score += 0.5 if ($var->is_clinvar_pathogenic_for_gene($gene) );
					$scaled_score += 0.2 if ($var->is_clinvar_pathogenic_for_gene($gene) && $gac < 1000);
					$scaled_score += 0.5 if $gac < 100 ;
					$scaled_score ++ if $gac < 30;
				}
				$scaled_score += $gene->score($self->use_phenotype());
				$max_scaled_score = $scaled_score if $max_scaled_score < $scaled_score;
			}
		}
		else {
			my $this_score;
			$this_score +=  $var->score_frequence_public();
			$this_score += $var->scaled_score_cadd();
			my $scaled_score = 1;
			$scaled_score = 2 if $this_score > 50;
			$scaled_score = 3 if $this_score >= 150;
			$scaled_score = 4 if $this_score >= 200;
			$max_scaled_score = $scaled_score if $max_scaled_score < $scaled_score;
		}
		$max_gene_score = $max_scaled_score;
	}
	return ($max_gene_score, $h_var_scores);
}

1;