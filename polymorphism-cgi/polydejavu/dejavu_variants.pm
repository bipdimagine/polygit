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

has user_name => (
	is		=> 'rw',
	lazy    => 1,
);

has pwd => (
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
#		if ($self->buffer->getQuery->isUserMagic($self->user_name(), $self->pwd())) {
#			@list_hash = @{$self->buffer->getQuery()->getAllProjects()};
#			$h_projects->{all}->{description} = 'all';
#			$h_projects->{all}->{name} = 'all';
#			$h_projects->{all}->{id} = 'all';
#		}
#		else {
			@list_hash = @{$self->buffer->getQuery()->getProjectListForUser($self->user_name(), $self->pwd())};
			foreach my $hash (@list_hash) {
				my $proj_name = $hash->{name};
				next unless ($proj_name =~ /NGS20/);
				$h_projects->{$proj_name}->{description} = $hash->{description};
				$h_projects->{$proj_name}->{name} = $proj_name;
				$h_projects->{$proj_name}->{id} = $hash->{id};
				$h_projects->{$hash->{id}} = $h_projects->{$proj_name};
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
				
				my $parquet = $self->buffer->dejavu_parquet_dir().'/'.$proj_name.'.'.$hash->{id}.'.parquet';
				$h_parquets->{$parquet} = undef if -e $parquet;
			}
#		}
		$self->{hash_users_projects_parquet} = $h_parquets;
		return $h_projects;
	}
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
	
	my ($h_variants, $hGenes);
	
	foreach my $chr_id (keys %$h_dv) {
        my $chr = $self->project->getChromosome($chr_id);
		my $nodv = $chr->rocks_dejavu();
		my @l_var_chr = keys %{$h_dv->{$chr_id}};
		my $nb_var = scalar (@l_var_chr);
		my $fork2 = $fork;
		my $chunk_size = 'auto';
		if ($nb_var > 100000) {
			#$fork2 = $fork/2;
			$chunk_size = 10000;
		}
		my $nb_part = 0;
		my $iter = natatime(100000, @l_var_chr);
		while( my @tmp = $iter->() ){
			$nb_part++;
			print '.part.'.$nb_part.'.';
			MCE::Loop->init(
			    max_workers => $fork2,
			    chunk_size  => $chunk_size,
			    gather      => sub {
			        my ($data) = @_;
			        return unless $data;
			        print '|';
			        if ($data->{variants}) {
			            foreach my $var_id (keys %{$data->{variants}}) {
			                $h_variants->{$var_id} = 1;
			            }
			        }
			        if ($data->{genes}) {
			            foreach my $gene_id (keys %{$data->{genes}}) {
			                $hGenes->{$gene_id} = $data->{genes}->{$gene_id};
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
					my $h_dv = $nodv->dejavu($rocks_id);
					next if not $h_dv;
		            if ($chr->intergenic_intspan->contains($var->start())) {
		                if (exists $h_cons->{intergenic}) {
		                    $hres->{genes}->{intergenic}->{$var_id} = 1;
		                    $hres->{variants}->{$var_id} = 1;
		                }
		            }
		            else {
		                my $genes = $var->getGenes();
		                foreach my $gene (@$genes) {
		                    my $var_annot = $var->variationTypeInterface($gene);
		                    my $is_ok = 0;
		                    foreach my $annot (split(',', $var_annot)) {
		                        $annot =~ s/ /_/g;
		                        if (exists $h_cons->{lc($annot)}) {
		                            $is_ok = 1;
		                            last;
		                        }
		                    }
		                    next unless $is_ok;
		                    $hres->{genes}->{$gene->id}->{$var_id} = 1;
		                    $hres->{variants}->{$var_id} = 1;
		                }
		            }
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
	
	#TODO: here 
#	my $limit_check = 50000;
#	if (scalar (keys %$var_id) > $limit_check) {
#		my @list_var_2;
#		my $j = 0;
#		foreach my $var (values %$var_id) {
#			push(@list_var_2, $var);
#			$j++;
#			$self->{alert_too_much_results} = 1 if $j == $limit_check;
#			last if $j == $limit_check;
#		}
#		@list_var = undef;
#		@list_var = @list_var_2;
#	}
	
	my $nb_var = scalar keys %$h_variants;
	print ".var.after.gad.dv.".$nb_var;

	my $h_parquets;
	my ($h_var_pos, @list_var_pos);
	foreach my $var_id (keys %$h_variants) {
		my $var = $self->project->_newVariant($var_id);
		my $chr = $var->getChromosome();
		my $chr_id = $chr->id();
		my $rocks_id = $var->rocksdb_id();
		$h_var_pos->{$chr_id}->{$var->start()}->{$var->var_allele()}->{ref_all} = $var->ref_allele();
		$h_var_pos->{$chr_id}->{$var->start()}->{$var->var_allele()}->{var_id} = $var->id();
		$h_var_pos->{$chr_id}->{$var->start()}->{$var->var_allele()}->{gnomad_id} = $var->gnomad_id();
		
		my $var_id = $var->id();
		push(@list_var_pos, $var->start());
		
		foreach my $project_id (sort keys %{$h_dv_rocks_ids->{$chr_id}->{$rocks_id}}) {
			my $project_name = $self->hash_users_projects->{$project_id}->{name};
			next if not exists $self->hash_users_projects->{all} and not exists $self->hash_users_projects->{$project_name};
			my $project_type_cache = $self->{hash_projects_ids_names}->{$project_id}->{type};
			$self->{hash_projects_ids_names}->{$project_id}->{name} = $project_name;
			$self->{hash_projects_ids_names}->{$project_id}->{type} = $project_type_cache;
			my $parquet = $self->buffer->dejavu_parquet_dir().'/'.$project_name.'.'.$project_id.'.parquet';
			$h_parquets->{$parquet} = undef if -e $parquet;
		}
		$var = undef;
	}
	
	my @list_parquet;
	foreach my $file (keys %$h_parquets) {
		push(@list_parquet, "'".$file."'");
	}
	
	print '.before_duck_projects.';
	my $hVariants_ok;
	
	my $nb_part = 0;
	my $iter = natatime(5000, @list_var_pos);
	while( my @tmp = $iter->() ){
		$nb_part++;
		print '.part.'.$nb_part.'.';
		my ($h_projects_patients, $h_gnomadid) = $self->get_from_duckdb_project_patients_infos_global(\@tmp, \@list_parquet, $h_var_pos);
		my $nb_2 = scalar (keys %{$h_projects_patients});
		print '.nb.'.$nb_2.'.';
		print '.after_duck_projects.';
		my $h_models = $self->models();
		my $min_ratio = $self->min_ratio();
		my $fork2 = $fork;
#		my $fork2 = $fork/2;
		$fork2 = 1 if ($nb_2 < 20);
		
		my $chunk_size;
		$chunk_size = int($nb_2/$fork2);
		print '.chunk.'.$chunk_size.'fork'.$fork2.'.';
		my ($B, $P);
		MCE::Loop->init(
		   max_workers => $fork2,
		   chunk_size => $chunk_size,
		   user_begin => sub {
		       $B = new GBuffer;
		       $P = $B->newProject(-name => $project_name);
		   },
#		   user_end => sub {
#		       $P->disconnect if $P;
#		   },
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
#				warn $var_id;
				foreach my $project_name (keys %{$h_projects_patients->{$var_id}}) {
					my ($h_pat_done, $h_pat_filtred);
					foreach my $patient_name (keys %{$h_projects_patients->{$var_id}->{$project_name}}) {
						push(@{$hres->{variants}->{$var_id}->{polyviewer_html}}, $h_projects_patients->{$var_id}->{$project_name}->{$patient_name}->{print_html});
						my $hh;
						$hh->{project_name} = $project_name;
						$hh->{patient_name} = $patient_name;
						$hh->{ratio} = $h_projects_patients->{$var_id}->{$project_name}->{$patient_name}->{ratio};
						$hh->{dp} = $h_projects_patients->{$var_id}->{$project_name}->{$patient_name}->{dp};
						$hh->{model} = $h_projects_patients->{$var_id}->{$project_name}->{$patient_name}->{model};
						if ($h_models) {
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
				if ($can_construct) {
					my $vp = PolyviewerVariant->new();
					$vp->setLmdbVariant($var);
					$vp->{hgenes} = {};
					$vp->{genes_id} = [];
					my $code = 0;
					if ($p->getChromosome($var->getChromosome->id())->intergenic_intspan->contains($var->start())) {
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
				}
			}
			MCE->gather($hres);
		} keys %{$h_projects_patients};	
		MCE::Loop->finish();
		$self->project->disconnect();
	}
	print '._end_MCE_check_.';
	return ($hGenes, $hVariants_ok);
}


sub get_table_project_patients_infos {
	my ($self, $project_name, $duck_he, $duck_patients, $duck_ratios) = @_;
	my $nb_he = $duck_he;
	my ($h_infos_patients, $h_tmp_pat);
	my $nb_pat = 0;
	foreach my $pat_id (unpack("w*",decode_base64($duck_patients))) {
		$nb_pat++;
		$h_infos_patients->{$nb_pat}->{id} = $pat_id;
		$h_tmp_pat->{$pat_id} = $nb_pat;
	}
	my $i = 0;
	$nb_pat = 1;
	foreach my $info (unpack("w*",decode_base64($duck_ratios))) {
		$i++;
		if ($i == 1) { $h_infos_patients->{$nb_pat}->{dp} = $info; }
		elsif ($i == 2) {
			my $ratio = '?';
			$ratio = ($info / $h_infos_patients->{$nb_pat}->{dp}) * 100 if ($h_infos_patients->{$nb_pat}->{dp});
			my $text = 'AC:'.$info;
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
	
	#TODO: here
	my $found_healthy_patient = 0;
	my $found_ill_patient = 0;
	foreach my $pat_id (keys %{$self->hash_users_projects->{$project_name}->{patients}}) {
		next if not exists $h_tmp_pat->{$pat_id};
		my $h_pat = $self->hash_users_projects->{$project_name}->{patients}->{$pat_id};
		my $icon = 'TODO';
		$h_infos_patients->{$h_tmp_pat->{$pat_id}}->{name} = $h_pat->{name};
		$h_infos_patients->{$h_tmp_pat->{$pat_id}}->{status} = $icon;
		$h_infos_patients->{$h_tmp_pat->{$pat_id}}->{sex} = '-';
		$h_infos_patients->{$h_tmp_pat->{$pat_id}}->{sex} = 'male' if $h_pat->{sex} eq '1';
		$h_infos_patients->{$h_tmp_pat->{$pat_id}}->{sex} = 'female' if $h_pat->{sex} eq '2';
		$h_infos_patients->{$h_tmp_pat->{$pat_id}}->{status_txt} = '-';
		if ($h_pat->{status} eq '1') {
			$h_infos_patients->{$h_tmp_pat->{$pat_id}}->{status_txt} = 'healthy';
			$found_healthy_patient++;
		}
		if ($h_pat->{status} eq '2') {
			$h_infos_patients->{$h_tmp_pat->{$pat_id}}->{status_txt} = 'ill';
			$found_ill_patient++;
		}
		$h_infos_patients->{$h_tmp_pat->{$pat_id}}->{family} = $h_pat->{family};
		$h_infos_patients->{$h_tmp_pat->{$pat_id}}->{description} = $self->hash_users_projects->{$project_name}->{description};
	}
	
	
#	my $b = new GBuffer;
#	my $p = $b->newProject( -name => $project_name);
#	my @lPat = @{$p->getPatients()};
#	return undef if scalar(@lPat) == 0;
#	my $found_healthy_patient = 0;
#	my $found_ill_patient = 0;
#	foreach my $pat (@lPat) {
#		next if not exists $h_tmp_pat->{$pat->id};
#		if ($pat->isParent) {
#			$h_infos_patients->{$nb_pat}->{model} = 'parent';
#		}
#		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{name} = $pat->name;
#		my $icon = $pat->small_icon();
#		$icon =~ s/"/'/g;
#		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{status} = $icon;
#		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{sex} = '-';
#		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{sex} = 'male' if $pat->sex() eq '1';
#		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{sex} = 'female' if $pat->sex() eq '2';
#		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{status_txt} = '-';
#		if ($pat->status() eq '1') {
#			$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{status_txt} = 'healthy';
#			$found_healthy_patient++;
#		}
#		if ($pat->status() eq '2') {
#			$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{status_txt} = 'ill';
#			$found_ill_patient++;
#		}
#		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{family} = $pat->getFamily->name();
#		$h_infos_patients->{$h_tmp_pat->{$pat->id}}->{description} = $p->description();
#	}
#	$p = undef;
#	$b =undef;
	
	#TODO: idem, je dois garder l'info des parents avec l'option only_ill si au moins un pat malade (faire h filters ici aussi pour le faire)
	my $ok_ill = 1;
	if ($self->only_ill_patients()) {
		$ok_ill = undef if $found_ill_patient == 0;
	}
	elsif ($self->only_strict_ill_patients()) {
		$ok_ill = undef if $found_healthy_patient > 0;
	}
	if ($ok_ill) {
		foreach my $id (keys %{$h_infos_patients}) {
			my $pat_name = $h_infos_patients->{$id}->{name};
			next if not $pat_name;
			$hres->{$pat_name} = $h_infos_patients->{$id};
		}
	}
	return $hres;
}

sub get_from_duckdb_project_patients_infos_global {
	my ($self, $list_var_pos, $list_files, $h_var_pos) = @_;
	print '.';	
	return if scalar(@$list_var_pos) == 0;
	return if scalar(@$list_files) == 0;
	my $text_pos = join(',', @$list_var_pos);
	my $fork = $self->fork();
#	my $sql = "PRAGMA threads=$fork; SELECT project,chr38,chr19,pos38,pos19,he,allele,patients,dp_ratios FROM read_parquet([".join(', ', @$list_files)."])
#					WHERE (allele='A' or allele='T' or allele='G' or allele='C') and pos38 IN ($text_pos)";
	my $sql = "PRAGMA threads=$fork; SELECT project,chr38,chr19,pos38,pos19,he,allele,patients,dp_ratios FROM read_parquet([".join(', ', @$list_files)."])
					WHERE pos38 IN ($text_pos)";
	my ($h_gnomadid, $h_projects_patients);
	my $duckdb = $self->buffer->software('duckdb');
	open(my $fh, "-|", "$duckdb -csv -c \"$sql\"") or die "duckdb failed";
	while (my $line = <$fh>) {
	    chomp $line;
	    my ($project_id,$this_chr38,$this_chr19,$this_pos38,$this_pos19,$he,$var_all,$patients,$dp_ratios) = split(/,/, $line);
    	next if $project_id eq 'project';
    	next if not exists $h_var_pos->{$this_chr38}->{$this_pos38}->{$var_all};
    	my $var_id = $h_var_pos->{$this_chr38}->{$this_pos38}->{$var_all}->{var_id};
    	my $gnomad_id = $h_var_pos->{$this_chr38}->{$this_pos38}->{$var_all}->{gnomad_id};
    	my $ref = $h_var_pos->{$this_chr38}->{$this_pos38}->{$var_all}->{ref_all};
		$h_gnomadid->{$gnomad_id}->{chr19} = $this_chr19;
		$h_gnomadid->{$gnomad_id}->{pos19} = $this_pos19;
		$h_gnomadid->{$gnomad_id}->{ref_all} = $ref;
		$h_gnomadid->{$gnomad_id}->{var_all} = $var_all;
		my $project_name = $self->{hash_projects_ids_names}->{$project_id}->{name};
		$h_projects_patients->{$var_id}->{$project_name} = $self->get_table_project_patients_infos($project_name, $he, $patients, $dp_ratios);
	}
	close($fh);
	print '.';
	return ($h_projects_patients, $h_gnomadid);
}

#sub get_from_duckdb_project_patients_infos {
#	my ($self, $var, $list_files) = @_;
#	return if scalar(@$list_files) == 0;
#	my $sql = "PRAGMA threads=3; SELECT project,chr38,chr19,pos38,pos19,he,allele,patients,dp_ratios FROM read_parquet([".join(', ', @$list_files)."])";
#	
#	my ($pos38, $alt38) = split('!', $var->rocksdb_id());
#	my ($h_gnomadid, $h_projects_patients);
#	
#	if ($var->getProject->current_genome_version() eq 'HG38') {
#		$sql .= " WHERE chr38='".$var->getChromosome->id()."' and pos38=".$var->start() ;
#		my $duckdb = $var->buffer->software('duckdb');
#		open(my $fh, "-|", "$duckdb -csv -c \"$sql\"") or die "duckdb failed";
#		while (my $line = <$fh>) {
#		    chomp $line;
#		    my ($project_id,$this_chr38,$this_chr19,$this_pos38,$this_pos19,$he,$var_all,$patients,$dp_ratios) = split(/,/, $line);
#	    	next if $project_id eq 'project';
#			next if $this_pos38 ne int($pos38);
#			$var_all =~ s/\+//;
#			if (not $var->isDeletion) {
#				next if $var_all ne $var->var_allele();
#			}
#			$h_gnomadid->{$var->gnomad_id()}->{chr19} = $this_chr19;
#			$h_gnomadid->{$var->gnomad_id()}->{pos19} = $this_pos19;
#			$h_gnomadid->{$var->gnomad_id()}->{ref_all} = $var->ref_allele();
#			$h_gnomadid->{$var->gnomad_id()}->{var_all} = $var->var_allele();
#			my $project_name = $self->{hash_projects_ids_names}->{$project_id}->{name};
#			$h_projects_patients->{$project_name} = $self->get_table_project_patients_infos($project_name, $he, $patients, $dp_ratios);
#		}
#		close($fh);
#	}
#	else {
#		warn "HG19";
#		warn $sql if $self->debug;
#		confesss("HG19!!");
#	}
#	return ($h_projects_patients, $h_gnomadid);
#}

sub print_html_gene {
	my ($self, $gene_id, $list_variants, $hVariantsDetails) = @_;
	print '|';
	my $found;
	foreach my $var_id (sort @$list_variants) {
		next if not exists $hVariantsDetails->{$var_id};
		next if not $hVariantsDetails->{$var_id}->{polyviewer_html_details_proj_pat};
		$found = 1;
	}
	return if not $found;
	
	my @l_var_ids = keys %{$hVariantsDetails};
	
	my $pr = $self->project;
	my $pat = $self->project->getPatients->[0];
	my $print_html = polyviewer_html->new( project=>$pr, patient=>$pat,header=>\@headers, bgcolor=>"background-color:#607D8B" );

	my ($g, $gene, $max_gene_score, $h_var_scores);
	if ($gene_id eq 'intergenic') {
		$g->{id} = $gene_id;
		$g->{uid} = $gene_id;
		$g->{name} = 'Intergenic';
		$g->{external_name} = 'Intergenic';
		$g->{chr_name} = 'Intergenic';
		$g->{nb} = scalar(@$list_variants);
		$g->{omim_id} = undef;
		$g->{pLI} = undef;
		$g->{phenotypes} = 'Intergenic';
		$g->{variants} = $list_variants;
		$g->{panels} = undef;
		($max_gene_score, $h_var_scores) = $self->get_score_variant_from_gene_without_patient($list_variants, $hVariantsDetails);
		$g->{max_score} = $max_gene_score;
	}
	else {
		$gene = $self->project->newGene($gene_id);
		$g->{id} = $gene_id;
		$g->{uid} = $gene_id;
		$g->{name} = $gene->external_name();
		$g->{external_name} = $gene->external_name();
		$g->{chr_name} = $gene->getChromosome->id();
		$g->{nb} = scalar(@$list_variants);
		$g->{omim_id} = $gene->omim_id();
		$g->{pLI} = $gene->pLI();
		$g->{phenotypes} = $gene->description();
		$g->{variants} = $list_variants;
		$g->{panels} = $self->buffer->queryPanel()->getPanelsForGeneName($gene->external_name);
		($max_gene_score, $h_var_scores) = $self->get_score_variant_from_gene_without_patient($list_variants, $hVariantsDetails, $gene);
		$g->{max_score} = $max_gene_score;
	}
	
	my $panel_id = "panel_".$gene_id;
	my ($out, $h_phenos);
	my $bg_color = $print_html->bgcolor;
	$out .= qq{<div class="panel-heading panel-face panel-grey"	style="$bg_color;height:43px;padding:10px;border:0px;width:100%;">};
	$out .= update_variant_editor::panel_gene($g, $panel_id, $self->project->name);
	$out .= qq{</div>};
	$out .= qq{<div style="height:3px;"></div>};
	$out .= qq{<div class="panel-body panel-collapse collapse" style="font-size: 09px;font-family:Verdana;" id="$panel_id">};
	$out .= qq{<table class="table table-striped table-condensed table-bordered table-hover table-mybordered" style="vertical-align:middle;text-align: center;font-size: 8px;font-family:  Verdana;line-height: 25px;min-height: 25px;height: 25px;box-shadow: 3px 3px 5px #555;">};
	$out .= $print_html->print_header("background-color:aliceblue;color:black");
	
	my $color_validation = "grey";
	my $i = 0;
	foreach my $var_id (sort @$list_variants) {
		next if not exists $hVariantsDetails->{$var_id};
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
		
		my @list_polyviewer_h_details = @{$hVariantsDetails->{$var_id}->{polyviewer_html_details_proj_pat}};
		$print_html->variant($polyviewer_variant);
		my ($this_out, $this_h_pheno) = $self->print_line_variant_all_patients($polyviewer_variant, $print_html, \@list_polyviewer_h_details, $opacity);
		$out .= $this_out;
		foreach my $pheno_name (keys %$this_h_pheno) { $h_phenos->{$pheno_name} = $this_h_pheno->{$pheno_name}; }
	}
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
		#$this_print_html->variant->{patients_calling}->{$pat->id()}->{gt} = '';
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
	$var_text .= "<br><i><b>HG19: ".$self->{hash_lift_variants}->{$polyviewer_variant->gnomad_id()}->{chr19}.'-'.$self->{hash_lift_variants}->{$polyviewer_variant->gnomad_id()}->{pos19}.'</b></i>';
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
				my $igv_b = qq{<button class='igvIcon2' onclick='launch_web_igv_js("$proj_name","$pnames","$f","$locus","/","$gn")' style="color:black"></button>};
				my $project_phenotypes = '';
				if ($self->{hash_users_projects}->{$proj_name}->{phenotypes}) {
					$project_phenotypes = qq{<span hidden>}.$self->{hash_users_projects}->{$proj_name}->{phenotypes_tags}.qq{</span>};
					$project_phenotypes .= qq{<span style='color:green;'><i>}.$self->{hash_users_projects}->{$proj_name}->{phenotypes}.qq{</i></span><br>};
				}
				my $project_description = $this_print_html->patient->getProject->name().', Description: '.$this_print_html->patient->getProject->description();
				push(@l_html_calling, "<tr style='padding:6px;'><td style='padding-right:10px;width:120px;'><center>".$project_phenotypes."<button onClick='alert(\"$project_description\")'>".$this_print_html->patient->getProject->name().'</button><br><b>'.$this_print_html->patient->getFamily->name()."</b><td style='padding-right:5px;'>".$igv_b."</td></center></td><td>".$html_pat."</td></tr>");
				push(@l_html_calling, "<tr style='padding:6px;'><td><br></td><td><br></td></tr>");
			}
		}
	}
	
	my $out_calling = qq{<center><table style='width:98%;max-height:200px;'>};
	$out_calling .= join('', @l_html_calling);
	$out_calling .= qq{</table></center>};
	
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
		next if not $hVariantsDetails->{$var_id}->{polyviewer_html_details_proj_pat};
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
			
		my $max_score_pat = -999;
		my @list_polyviewer_h_details = @{$hVariantsDetails->{$var_id}->{polyviewer_html_details_proj_pat}};
		foreach my $hpat (@list_polyviewer_h_details) {
			my $score_pat = $max_scaled_score;
			my $pc = $hpat->{ratio};
			$pc =~ s/.*Ratio://g;
			my $dp = $hpat->{dp};
			
			if (lc($hpat->{model}) eq 'solo') {
				$score_pat += 0.4 if $pc > 70;
				$score_pat += 0.3 if  ($pc > 90);
				$score_pat += 0.3 if $dp > 10;
				my $n = $var->getGnomadHO();
				$n = 0 unless $n;
				$n = 0 if $n eq "-";
				$score_pat += 0.5 if $n == 0;
				$score_pat -= 0.4 if $n > 3;
			}
			
			if (lc($hpat->{model}) eq 'both' or lc($hpat->{model}) eq 'parent') {
				$score_pat -= 2;
			}
				
			if (lc($hpat->{model}) eq 'recessif' or lc($hpat->{model}) eq 'recessive') {
				$score_pat += 0.4 if $pc > 70;
				$score_pat += 0.3 if  ($pc > 90);
				$score_pat += 0.3 if $dp > 10;
				my $n = $var->getGnomadHO();
				$n = 0 unless $n;
				$score_pat += 0.5 if $n == 0;
				$score_pat -= 0.4 if $n > 3;
				$score_pat =  $score_pat*2 ;
			}
			
			if (lc($hpat->{model}) eq 'mother' or lc($hpat->{model}) eq 'father') {
			 	$score_pat += 0.3 if $pc > 15;
			 	$score_pat += 0.4 if $pc > 20;
				$score_pat += 0.3 if $dp > 10;
				my $n = $var->getGnomadHO();
				$n = 0 unless $n;
				$n = 0 if $n eq "-";
				$score_pat += 0.25 if $n == 0;
				$score_pat -= 0.4 if $n > 3;
			}
			
			if (lc($hpat->{model}) eq 'denovo' or lc($hpat->{model}) eq 'strict-denovo' or lc($hpat->{model}) eq 'strict_denovo') {
				my $n = $var->getGnomadAC();
				$n = 0 unless $n;
				$score_pat += 0.5 if $n == 0;
				$score_pat -= 0.4 if $n > 30;
				if (lc($hpat->{model}) eq 'strict-denovo' or lc($hpat->{model}) eq 'strict_denovo'){
					$score_pat += 0.3 if $pc > 20  ;
					$score_pat += 0.4 if $pc > 30  ;
					$score_pat += 0.3 if $dp > 10;
				}
				else {
					$score_pat += 0.1 if   $pc > 20  ;
					$score_pat += 0.2 if   $pc > 30  ;
					$score_pat -= 0.5;# if   $pc > 20  ;
				}
				$score_pat += 0.4 if ($pc > 20 && ($pc < 75 && $var->getChromosome->name ne "X"));
				$score_pat += 0.3 if $dp > 10;
				$score_pat =  $score_pat*2 ;
			}
			$max_score_pat = $score_pat if $max_score_pat < $score_pat;
		}
		if ($max_score_pat > $max_score) {
			$max_score = $max_score_pat;
		}
		$h_var_scores->{$var->id()} = $max_score;
		$max_gene_score = $max_score if $max_gene_score < $max_score;
	}
	return ($max_gene_score, $h_var_scores);
}

1;