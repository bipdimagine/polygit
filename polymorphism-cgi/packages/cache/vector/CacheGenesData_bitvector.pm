package CacheGenesData_bitvector;
use strict;

use FindBin qw($RealBin);
use Data::Dumper;
use Parallel::ForkManager;
use Bio::DB::Sam;
use Storable qw(store retrieve freeze dclone thaw);
use Time::Duration;
use Math::Combinatorics;
use Set::IntSpan::Fast::XS;
use Bit::Vector;
use Bit::Vector::Overload;
use Set::IntRange; 
use CacheGenesData_nodb;
use Tabix;
use JSON::XS;
use Bio::DB::Sam;
use POSIX qw(strftime);
use lib "$RealBin";
use lib "$RealBin/../../GenBo/lib/obj-nodb/";
use CacheDiag;
use validationQuery;
use GenBoNoSql;
use String::ProgressBar;
use Tabix;
use GenBoNoSqlLmdb;

my @lSubdata_filters = ('all', 'exome-intronic', 'exome');
my ($hset, $hCatUsed);
my $himpact_sorted = {
	"high" 		=> "4",
	"moderate" 	=> "3",
	"low" 		=> "1",
	"unknown"	=> "0",
};
my $fork = 1;
my $verbose = 1;
my $hLargeDelIds;




sub create_cache_genes {
	my ($project, $chr_name, $fork_wanted, $limit) = @_;
	$fork = $fork_wanted if ($fork_wanted);
	my $buffer = $project->buffer();
	$buffer->dbh_reconnect();
	$verbose = undef unless ($project->cache_verbose());
	my $project_name = $project->name();
	my $chr = $project->getChromosome($chr_name);
	
	my ($hRefCoord, $limit) = getHashReferencesCoord($project, $chr, $limit);
	my $nb_max = scalar(keys %$hRefCoord);
	my $dir_temp = $project->getCacheBitVectorDir().'/tmp_chr'.$chr_name.'/';
	`rm -r $dir_temp` if (-d $dir_temp);
	`mkdir $dir_temp`;
	my @lPat = @{$project->getPatients()};

	print "\n### STEP 1: get variants / genes ids and dejavu (limit=".$limit."pb, fork=$fork)\n" if ($verbose);
	print strftime "Start at %H:%M:%S\n", localtime if ($verbose);
    my $fileout = $dir_temp.'/../'.$chr_name.'.1.freeze';
    my ($hVarIds_all, $hGenesIds_all, $nb_subpart) ;
  #  unless (-e $fileout){
		 ($hVarIds_all, $hGenesIds_all, $nb_subpart) = getIdsAndDejaVu($dir_temp, $hRefCoord, $project, $chr, $fork);
	
		updateDejaVuPn($project,$chr_name,$fork);

	
# 		my $h;
#    		$h->{var}=$hVarIds_all;
#    		$h->{gene}=$hGenesIds_all;
#    		$h->{subpart}=$nb_subpart;
#    		store($h, $fileout);
#	#}
#	#else {
#		warn "*-*retrieve*-*";
#		my $h = retrieve($fileout);
#		$hVarIds_all = $h->{var};
#		$hGenesIds_all= $h->{gene};
#		$nb_subpart = $h->{subpart};
#	#}


	print strftime "\nEnd at %H:%M:%S\n", localtime if ($verbose);
	unless ($hVarIds_all) {
		print "\n\n-> NO VARIANTS !\n" if ($verbose);
		`rm -r $dir_temp`;
		return 1;
	}
	my $nb_var = scalar(keys %{$hVarIds_all});
	print "\n\n### STEP 2: ALL subdata (nb_var=".$nb_var.", fork=$fork)\n" if ($verbose);
	print strftime "Start at %H:%M:%S\n", localtime if ($verbose);
	my ($hSet_all, $hCatUsed_all, $size_genes) = doCache_all($dir_temp, $project, $chr_name, $hVarIds_all, $hGenesIds_all, 'all', $nb_subpart);
	print strftime "End at %H:%M:%S\n", localtime if ($verbose);
 	my $h2;
 	$h2->{set}=$hSet_all;
 	$h2->{cat}=$hCatUsed_all;
 	$h2->{size}=$size_genes;
	print "\n\n### STEP 3: EXOME-INTRONIC subdata\n" if ($verbose);
	print strftime "Start at %H:%M:%S\n", localtime if ($verbose);
	doCache_other($dir_temp, $project, $chr_name, $hSet_all, $hCatUsed_all, $size_genes, 'exome-intronic', 'intergenic');
	print strftime "End at %H:%M:%S\n", localtime if ($verbose);
	`rm -r $dir_temp`;
	return;
}

sub getHashBitVector_other {
	my ($hCatUsed, $hset, $size_var, $size_genes) = @_;
	my $hVector;
	foreach my $key (keys %$hCatUsed) {
		my @lFields = split(';', $key);
		my @lIds;
		if  (scalar @lFields == 3) {
			@lIds = sort {$a <=> $b} keys %{$hset->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}};
			unless (exists $hVector->{$lFields[1]}->{$lFields[2]}) {
				$hVector->{$lFields[1]}->{$lFields[2]} = Bit::Vector->new_Enum( $size_genes, join(',', @lIds) );
			}
			else {
				my $v = Bit::Vector->new_Enum( $size_genes, join(',', @lIds) );
				$hVector->{$lFields[1]}->{$lFields[2]} += $v;
			}
		}
		elsif (scalar @lFields == 5) {
			@lIds = sort {$a <=> $b} keys %{$hset->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}};
			if ($lFields[1] eq 'genes') {
				unless (exists $hVector->{$lFields[1]}->{$lFields[2]}->{start_subpart}) {
					my @lTmp = sort {$a <=> $b} keys %{$hset->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{'global_categories'}->{'all'}};
					$hVector->{$lFields[1]}->{$lFields[2]}->{start_subpart} = $lTmp[0];
					$hVector->{$lFields[1]}->{$lFields[2]}->{end_subpart} = $lTmp[-1];
					$hVector->{$lFields[1]}->{$lFields[2]}->{size} = $lTmp[-1] - $lTmp[0] + 1;
				}
				my $this_size = $hVector->{$lFields[1]}->{$lFields[2]}->{size};
				my @lNewIds;
				my $i = $hVector->{$lFields[1]}->{$lFields[2]}->{start_subpart};
				my $j = 0;
				while ($i <= $hVector->{$lFields[1]}->{$lFields[2]}->{end_subpart}) {
					push(@lNewIds, $j) if (exists $hset->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->{$i});
					$i++;
					$j++;
				}
				unless (exists $hVector->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}) {
					$hVector->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]} = Bit::Vector->new_Enum( $this_size, join(',', @lNewIds) );
				}
				else {
					my $v = Bit::Vector->new_Enum( $this_size, join(',', @lNewIds) );
					$hVector->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]} += $v;
				}
			}
			else {
				unless (exists $hVector->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}) {
					$hVector->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]} = Bit::Vector->new_Enum( $size_var, join(',', @lIds) );
				}
				else {
					my $v = Bit::Vector->new_Enum( $size_var, join(',', @lIds) );
					$hVector->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]} += $v;
				}
			}
		}
		elsif (scalar @lFields == 6) {
			@lIds = sort {$a <=> $b} keys %{$hset->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->{$lFields[5]}};
			unless (exists $hVector->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->{$lFields[5]}) {
				$hVector->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->{$lFields[5]} = Bit::Vector->new_Enum( $size_var, join(',', @lIds) );
			}
			else {
				my $v = Bit::Vector->new_Enum( $size_var, join(',', @lIds) );
				$hVector->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->{$lFields[5]} += $v;
			}
		}
		else { die; }
	}
	return $hVector;
}

sub getHashBitVector_all {
	my ($dir_temp, $chr_name, $size_var, $size_genes, $nb_max) = @_;
	my $i_progress = 1;				
	my $hVector;
	my ($hSet_all, $hCatUsed_all);
	my $nb_key = 0;
	while ($nb_key < $nb_max) {
		$nb_key++;
		my ($hset, $hCatUsed) = prepareHashesAfterForkAnnot($dir_temp, $chr_name, $nb_key);
		foreach my $key (keys %$hCatUsed) {
			$hCatUsed_all->{$key} = undef;
			my @lFields = split(';', $key);
			my @lIds;
			if  (scalar @lFields == 3) {
				@lIds = sort {$a <=> $b} keys %{$hset->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}};
				foreach my $id (@lIds) { $hSet_all->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$id} = undef; }
				unless (exists $hVector->{$lFields[1]}->{$lFields[2]}) {
					$hVector->{$lFields[1]}->{$lFields[2]} = Bit::Vector->new_Enum( $size_genes, join(',', @lIds) );
				}
				elsif(scalar(@lIds) > 0) {
					my $v = Bit::Vector->new_Enum( $size_genes, join(',', @lIds) );
					$hVector->{$lFields[1]}->{$lFields[2]} += $v;
				}
			}
			elsif (scalar @lFields == 5) {
				@lIds = sort {$a <=> $b} keys %{$hset->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}};
				foreach my $id (@lIds) { $hSet_all->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->{$id} = undef; }
				if ($lFields[1] eq 'genes') {
					unless (exists $hVector->{$lFields[1]}->{$lFields[2]}->{start_subpart}) {
						my @lTmp = sort {$a <=> $b} keys %{$hset->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{'global_categories'}->{'all'}};
						$hVector->{$lFields[1]}->{$lFields[2]}->{start_subpart} = $lTmp[0];
						$hVector->{$lFields[1]}->{$lFields[2]}->{end_subpart} = $lTmp[-1];
						$hVector->{$lFields[1]}->{$lFields[2]}->{size} = $lTmp[-1] - $lTmp[0] + 1;
					}
					else {
						my @lTmp = sort {$a <=> $b} keys %{$hset->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{'global_categories'}->{'all'}};
						$hVector->{$lFields[1]}->{$lFields[2]}->{end_subpart} = $lTmp[-1] if ($hVector->{$lFields[1]}->{$lFields[2]}->{end_subpart} < $lTmp[-1]);
						$hVector->{$lFields[1]}->{$lFields[2]}->{size} = $lTmp[-1] - $hVector->{$lFields[1]}->{$lFields[2]}->{start_subpart} + 1;
						my $this_size = $hVector->{$lFields[1]}->{$lFields[2]}->{size};
						foreach my $cat_name (keys %{$hVector->{$lFields[1]}->{$lFields[2]}->{global_categories}}) {
							$hVector->{$lFields[1]}->{$lFields[2]}->{global_categories}->{$cat_name}->Resize($this_size);
						}
						foreach my $cat_name (keys %{$hVector->{$lFields[1]}->{$lFields[2]}->{categories}}) {
							$hVector->{$lFields[1]}->{$lFields[2]}->{categories}->{$cat_name}->Resize($this_size);
						}
					}
					my $this_size = $hVector->{$lFields[1]}->{$lFields[2]}->{size};
					my @lNewIds;
					my $i = $hVector->{$lFields[1]}->{$lFields[2]}->{start_subpart};
					my $j = 0;
					while ($i <= $hVector->{$lFields[1]}->{$lFields[2]}->{end_subpart}) {
						push(@lNewIds, $j) if (exists $hset->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->{$i});
						$i++;
						$j++;
					}
					unless (exists $hVector->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}) {
						$hVector->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]} = Bit::Vector->new_Enum( $this_size, join(',', @lNewIds) );
					}
					elsif(scalar(@lIds) > 0) {
						$hVector->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->Resize($this_size);
						my $v = Bit::Vector->new_Enum( $this_size, join(',', @lNewIds) );
						$hVector->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]} += $v;
					}
				}
				else {
					unless (exists $hVector->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}) {
						$hVector->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]} = Bit::Vector->new_Enum( $size_var, join(',', @lIds) );
					}
					elsif(scalar(@lIds) > 0) {
						my $v = Bit::Vector->new_Enum( $size_var, join(',', @lIds) );
						$hVector->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]} += $v;
					}
				}
			}
			elsif (scalar @lFields == 6) {
				@lIds = sort {$a <=> $b} keys %{$hset->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->{$lFields[5]}};
				foreach my $id (@lIds) { $hSet_all->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->{$lFields[5]}->{$id} = undef; }
				unless (exists $hVector->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->{$lFields[5]}) {
					$hVector->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->{$lFields[5]} = Bit::Vector->new_Enum( $size_var, join(',', @lIds) );
				}
				elsif(scalar(@lIds) > 0) {
					my $v = Bit::Vector->new_Enum( $size_var, join(',', @lIds) );
					$hVector->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->{$lFields[5]} += $v;
				}
			}
			else { die; }
		}
		$hset = undef;
		$hCatUsed = undef;
	}
	return ($hVector, $hSet_all, $hCatUsed_all);
}

sub changeIds {
	my $hIds = shift;
	my $hIds_new;
	my $new_id = 1;
	foreach my $old_id (sort {$a <=> $b} keys %{$hIds->{ok}}) {
		$hIds_new->{$old_id} = $new_id;
		$new_id++;
	}
	return ($hIds_new, $new_id);
}

sub doCache_other {
	my ($dir_temp, $project, $chr_name, $hSet_all, $hCatUsed_all, $size_genes, $subdata_filter, $filters) = @_;

	print "Filters: $filters\n" if ($verbose);
	my $project_name = $project->name();
	print "Filter global informations\n" if ($verbose);
	my @lFilters = split(';', $filters);
	my ($hSet_other, $hCatUsed_other, $hIds);
	foreach my $key (keys %$hCatUsed_all) {
	
		my @lFields = split(';', $key);
		my $ok = 1;
		foreach my $filter (@lFilters) {
			if ($filter eq $lFields[-1]) {
				$ok = undef;
			}
		}
		
		next unless ($ok);
		#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		#!!!!!!!!!! ca sert a quoi et pourquoi certaine fois all_global et pas d'autre ?????
		#next unless ($lFields[1] eq 'patients');
		#
		$hCatUsed_other->{$key} = undef;
		my @lIds;
#		if (scalar @lFields == 3) {
#			@lIds = keys %{$hSet_all->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}};
#		}
		if (scalar @lFields == 5) {
			@lIds = keys %{$hSet_all->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}};
		}
		elsif (scalar @lFields == 6) {
			@lIds = keys %{$hSet_all->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->{$lFields[5]}};
		}
		foreach my $id (@lIds) { $hIds->{ok}->{int($id)} = undef; }
	}
	
	print "Create new ids\n" if ($verbose);
	my ($hIds_new, $size_var) = changeIds($hIds);
	
	print "Create new bit_vectors\n" if ($verbose);
	foreach my $key (keys %$hCatUsed_other) {
		my @lIds;
		my @lFields = split(';', $key);
		if (scalar @lFields == 3) {
			@lIds = keys %{$hSet_all->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}};
			foreach my $id (@lIds) {
				next unless (exists $hIds_new->{$id});
				my $id_new = $hIds_new->{$id};
				$hSet_other->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$id_new} = undef;
			}
		}
		elsif (scalar @lFields == 5) {
			@lIds = keys %{$hSet_all->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}};
			foreach my $id (@lIds) {
				next unless (exists $hIds_new->{$id});
				my $id_new = $hIds_new->{$id};
				$hSet_other->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->{$id_new} = undef;
			}
		}
		elsif (scalar @lFields == 6) {
			@lIds = keys %{$hSet_all->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->{$lFields[5]}};
			foreach my $id (@lIds) {
				next unless (exists $hIds_new->{$id});
				my $id_new = $hIds_new->{$id};
				$hSet_other->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->{$lFields[5]}->{$id_new} = undef;
			}
		}
		
	}

	my $nosql_var_ids_all = GenBoNoSql->new(dir => $dir_temp.'../all_variants_ids', mode => 'r');
	my $nosql_var_ids_new = GenBoNoSql->new(dir => $dir_temp.'../'.$subdata_filter.'_variants_ids', mode => 'c');
	#$nosql_var_ids_new->clear($chr_name);
	my $h_var_ids_all = $nosql_var_ids_all->get_bulk($chr_name);
	my $h_var_ids_new;
	foreach my $id (keys %$hIds_new) {
		my $var_id = $h_var_ids_all->{$id};
		$h_var_ids_new->{$hIds_new->{$id}} = $var_id;
		$h_var_ids_new->{$var_id} = $hIds_new->{$id};
	}
	
	$nosql_var_ids_new->put_bulk($chr_name, $h_var_ids_new);
	$nosql_var_ids_all->close();
	$nosql_var_ids_new->close();
	
	delete $hSet_other->{all}->{all_genes};
	my ($hSet_tmp, $hCat_tmp) = createAllGeneKey($hSet_other);
	$hSet_other->{all}->{all_genes} = $hSet_tmp;
	foreach my $key (keys %$hCat_tmp) { $hCatUsed->{$key} = undef; }
	
	my $hVector;
	$hVector->{$subdata_filter} = getHashBitVector_other($hCatUsed_other, $hSet_other, $size_var, $size_genes);
	
	foreach my $patient (@{$project->getPatients()}) {
		my $pat_name = $patient->name();
		unless (exists $hVector->{$subdata_filter}->{patients}->{$pat_name}) {
			$hVector->{$subdata_filter}->{patients}->{$pat_name}->{categories}->{intergenic} = Bit::Vector->new( $size_var );
		}
		foreach my $cat_name ('substitution', 'insertion', 'deletion', 'ho', 'he') {
			unless (exists $hVector->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{$cat_name}) {
				$hVector->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{$cat_name} = Bit::Vector->new( $size_var );
			}
		}
	}
	
	my $exp_var = scalar(keys %{$hIds->{ok}});
	my ($is_ok) = checkResults($exp_var, $subdata_filter, $hVector, $hIds_new) if ($exp_var > 0);
	unless ($is_ok) {
		my $cache_dir = $project->getCacheBitVectorDir();
		warn "\n\nERROR in cache process... Supress $cache_dir. Die.\n\n";
		if (-d $cache_dir) {
			my $new_dir = $cache_dir;
			$new_dir =~ s/\/vector\//\/vector.pb_cache\//;
			`mv $cache_dir $new_dir`;
		}
		die();
	}
	
	mkdir $dir_temp.'../'.$subdata_filter.'_global/' unless (-d $dir_temp.'../'.$subdata_filter.'_global/');
	my $fileout = $dir_temp.'../'.$subdata_filter.'_global/'.$chr_name.'.freeze';
	store($hVector->{$subdata_filter}, $fileout);
}

sub createAllGeneKey {
	my $hSet = shift;
	my ($hSet_tmp, $hCat_tmp);
	foreach my $gene_id (keys %{$hSet->{all}->{genes}}) {
		foreach my $type_cat (keys %{$hSet->{all}->{genes}->{$gene_id}}) {
			foreach my $cat (keys %{$hSet->{all}->{genes}->{$gene_id}->{$type_cat}}) {
				$hSet_tmp->{$cat}->{$gene_id} = undef;
				$hCat_tmp->{"all;all_genes;$cat"} = undef;
			}
		}
	}
	return ($hSet_tmp, $hCat_tmp);
}

sub prepareHashesAfterForkAnnot {
	my ($dir_temp, $chr_name, $nb_key) = @_;
	my $nosql_annot = GenBoNoSql->new(dir => $dir_temp.'/annot', mode => 'r');
	my $hset_all = $nosql_annot->get_bulk($nb_key);
	$nosql_annot->close();
	my $nosql_cat_used = GenBoNoSql->new(dir => $dir_temp.'/cat_used', mode => 'r');
	my $hCatUsed_all = $nosql_cat_used->get_bulk($nb_key);
	$nosql_cat_used->close();
	my ($hset, $hCatUsed);
	foreach my $nb (keys %$hCatUsed_all) {
		foreach my $key (keys %{$hCatUsed_all->{$nb}}) {
			$hCatUsed->{$key} = undef;
			my @lFields = split(';', $key);
			if (scalar @lFields == 3) {
				foreach my $id (keys %{$hset_all->{$nb}->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}}) {
					$hset->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$id} = undef;
				}
			}
			elsif (scalar @lFields == 5) {
				foreach my $id (keys %{$hset_all->{$nb}->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}}) {
					$hset->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->{$id} = undef;
				}
			}
			elsif (scalar @lFields == 6) {
				foreach my $id (keys %{$hset_all->{$nb}->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->{$lFields[5]}}) {
					$hset->{$lFields[0]}->{$lFields[1]}->{$lFields[2]}->{$lFields[3]}->{$lFields[4]}->{$lFields[5]}->{$id} = undef;
				}
			}
		}
	}
	
	my ($hSet_tmp, $hCat_tmp) = createAllGeneKey($hset);
	$hset->{all}->{all_genes} = $hSet_tmp;
	foreach my $key (keys %$hCat_tmp) { $hCatUsed->{$key} = undef; }
	return ($hset, $hCatUsed);
}
sub openNoSql {
	my ($dir_temp, $chr_name, $mode) = @_;
	my $hConnect;
	$hConnect->{var_infos} = GenBoNoSql->new(dir => $dir_temp.'../variants_infos', mode => $mode);
	$hConnect->{var_ids} = GenBoNoSql->new(dir => $dir_temp.'../all_variants_ids', mode => $mode);
	$hConnect->{genes_ids} = GenBoNoSql->new(dir => $dir_temp.'../genes_ids', mode => $mode);
	$hConnect->{genes_infos} = GenBoNoSql->new(dir => $dir_temp.'../genes_infos', mode => $mode);
	$hConnect->{cat_used} = GenBoNoSql->new(dir => $dir_temp.'cat_used', mode => $mode);
	$hConnect->{annot} = GenBoNoSql->new(dir => $dir_temp.'annot', mode => $mode);
	foreach my $name (keys %$hConnect) {
		$hConnect->{$name}->cache_limit(2000);
	}
	return $hConnect;
}

sub closeNoSql {
	my ($hConnect) = @_;
	foreach my $key (keys %$hConnect) {
		$hConnect->{$key}->close();
	}
}

sub getVarFrequencies {
	my ($variation) = shift;
	my ($freq_text, $freq_ho_text, $freq_he_text);
	$freq_text = "freq_none";
	$freq_ho_text = "freq_ho_none";
	$freq_he_text = "freq_he_none";
	unless ($variation->isNoFreqVariant()) {
		my $freq = $variation->frequency();
		if    ($freq == -1)     { $freq_text = "freq_none"; }
		elsif ($freq <= 0.0001) { $freq_text = "freq_0001"; }
		elsif ($freq <= 0.001)  { $freq_text = "freq_001"; }
		elsif ($freq <= 0.01)   { $freq_text = "freq_01"; }
		elsif ($freq <= 0.05)   { $freq_text = "freq_05"; }
		else                    { $freq_text = "freq_1"; }
		my $freq_ho = $variation->homozygous_frequency();
		if    ($freq_ho == -1)     { $freq_ho_text = "freq_ho_none"; }
		elsif ($freq_ho <= 0.0001) { $freq_ho_text = "freq_ho_0001"; }
		elsif ($freq_ho <= 0.001)  { $freq_ho_text = "freq_ho_001"; }
		elsif ($freq_ho <= 0.01)   { $freq_ho_text = "freq_ho_01"; }
		elsif ($freq_ho <= 0.05)   { $freq_ho_text = "freq_ho_05"; }
		else                       { $freq_ho_text = "freq_ho_1"; }
		my $freq_he = $variation->heterozygous_frequency();
		if    ($freq_he == -1)     { $freq_he_text = "freq_he_none"; }
		elsif ($freq_he <= 0.0001) { $freq_he_text = "freq_he_0001"; }
		elsif ($freq_he <= 0.001)  { $freq_he_text = "freq_he_001"; }
		elsif ($freq_he <= 0.01)   { $freq_he_text = "freq_he_01"; }
		elsif ($freq_he <= 0.05)   { $freq_he_text = "freq_he_05"; }
		else                       { $freq_he_text = "freq_he_1"; }
	}
	return ($freq_text, $freq_ho_text, $freq_he_text);
}

sub getVarTypes {
	my $variation = shift;
	my @lType;
	if    ($variation->unknown())       { push(@lType, "new"); } 
	if    ($variation->isClinical())    { push(@lType, "pheno_snp"); } 
	if    ($variation->isVariation()) 	{ push(@lType, "substitution"); }
	elsif ($variation->isInsertion()) 	{ push(@lType, "insertion"); }
	elsif ($variation->isDeletion())  	{
		if (exists $hLargeDelIds->{$variation->id()}) {
			push(@lType, "large_deletion");
		}
		else { push(@lType, "deletion"); }
	}
	else {
		warn "\n\nERROR: no type found for ".$variation->id()." ! Die...\n\n";
		die;
	}
	return \@lType;
}

sub doCache_all {
	my ($dir_temp, $project, $chr_name, $hVarIds, $hGenesIds, $subdata_filter, $max) = @_;
	my $project_name = $project->name();
	my $hConnect = openNoSql($dir_temp, $chr_name, 'c');
	my $h_fork;
	my $nb_key = 0;
	my $hGlobal;
	
	my $i_progress = 1;
	my $pr;
	if ($verbose) {		
		$pr = String::ProgressBar->new( max => $max, info => 'Check each variant annotations' );
		$pr->write();
	}
	my $pm = new Parallel::ForkManager($fork, $dir_temp);
	$pm->run_on_finish (
		sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hRes) = @_;
			if (defined($hRes)) {
				foreach my $key (keys %$hRes) {
					if ($key eq 'var_infos') {
						my $hTmp;
						foreach my $i (keys %{$hRes->{$key}}) {
							foreach my $pat_name (keys %{$hRes->{$key}->{$i}}) {
								foreach my $var_id (keys %{$hRes->{$key}->{$i}->{$pat_name}}) {
									$hTmp->{$var_id}->{$pat_name} =  $hRes->{$key}->{$i}->{$pat_name}->{$var_id};
								}
							}
						}
						$hConnect->{$key}->put_bulk($chr_name, $hTmp);
						$hTmp = undef;
						$hRes->{$key} = undef;
					}
					elsif ($key eq 'var_ids') {
						$hConnect->{'var_ids'}->put_bulk($chr_name, $hRes->{'var_ids'});
						$hRes->{$key} = undef;
					}
					elsif ($key eq 'genes_ids') {
						$hConnect->{'genes_ids'}->put_bulk($chr_name, $hRes->{'genes_ids'});
						$hRes->{$key} = undef;
					}
					elsif ($key eq 'genes_infos') {
						$hConnect->{'genes_infos'}->put_bulk($chr_name, $hRes->{'genes_infos'});
						$hRes->{$key} = undef;
					}
					elsif ($key eq 'cat_used') {
						$hConnect->{'cat_used'}->put_bulk($hRes->{'sub_part'}, $hRes->{'cat_used'});
						$hRes->{$key} = undef;
					}
					elsif ($key eq 'annot') {
						$hConnect->{'annot'}->put_bulk($hRes->{'sub_part'}, $hRes->{'annot'});
						$hRes->{$key} = undef;
					}
				}
				$hRes = undef;
				resetHset();
				$hRes = undef;
				$i_progress++;
				if ($verbose) {
					$pr->update($i_progress);
					$pr->write();
				}
			}
			else { print qq|No message received from child process $pid!\n|; }
		}
	);
	
	my $this_part = 0;
	while ($this_part <= $max) {
		$this_part++;
		my $pid = $pm->start and next;
		my $nosql_patients = GenBoNoSql->new(dir => $dir_temp.'tmp_var_patients', mode => 'r');
		my $hVarAnnex = $nosql_patients->get($chr_name, $this_part);
		$nosql_patients->close();
		my ($hFork, $hFork_catUsed);
		my ($h_var_ids, $h_var_infos, $h_genes_ids, $h_genes_infos);
		$project->getChromosomes();
	    foreach my $var_id (keys %$hVarAnnex) {
	    	die($var_id) unless exists  $hVarIds->{$var_id};
			my $vector_var_ids = $hVarIds->{$var_id};
			unless ($vector_var_ids) {
				warn ("\n\n### ERROR: $var_id don't have $vector_var_ids\n\n");
				die;
			}
			unless ($var_id) {
				warn ("\n\n### ERROR: non var_id...\n\n");
				die;
			}
			my $variation = $project->_newVariant($var_id);
			die unless($variation);
			$h_var_ids->{$vector_var_ids} = $var_id;
			$h_var_ids->{$var_id} = $vector_var_ids;
			foreach my $pat_name (keys %{$hVarAnnex->{$var_id}->{'patients'}}) {
				my ($chr_name, $pos_vcf, $ref_vcf, $var_vcf) = split('_', $hVarAnnex->{$var_id}->{'vcf_infos'});
				$h_var_infos->{$pat_name}->{$var_id}->{chr_vcf} = $chr_name;
				$h_var_infos->{$pat_name}->{$var_id}->{pos_vcf} = $pos_vcf;
				$h_var_infos->{$pat_name}->{$var_id}->{ref_vcf} = $ref_vcf;
				$h_var_infos->{$pat_name}->{$var_id}->{var_vcf} = $var_vcf;
				$h_var_infos->{$pat_name}->{$var_id}->{he_ho_details} = $hVarAnnex->{$var_id}->{'patients'}->{$pat_name}->{'he_ho_details'};
			}
			# calcul des scores de frequence pour la variation (PolyDiag)
			my $bank = $variation->database();
			my (@lType) = @{getVarTypes($variation)};
			my ($freq_text, $freq_ho_text, $freq_he_text) = getVarFrequencies($variation);
			
			foreach my $pat_name (keys %{$hVarAnnex->{$var_id}->{'patients'}}) {
				$hFork_catUsed->{"$subdata_filter;patients;$pat_name;global_categories;$freq_text"} = undef;
				$hFork->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{$freq_text}->{$vector_var_ids} = undef;
				$hFork_catUsed->{"$subdata_filter;patients;$pat_name;global_categories;$freq_ho_text"} = undef;
				$hFork->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{$freq_ho_text}->{$vector_var_ids} = undef;
				$hFork_catUsed->{"$subdata_filter;patients;$pat_name;global_categories;$freq_he_text"} = undef;
				$hFork->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{$freq_he_text}->{$vector_var_ids} = undef;
				$hFork_catUsed->{"$subdata_filter;patients;$pat_name;global_categories;$bank"} = undef;
				$hFork->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{$bank}->{$vector_var_ids} = undef;
		 		foreach my $type (@lType) {
					$hFork_catUsed->{"$subdata_filter;patients;$pat_name;global_categories;$type"} = undef;
					$hFork->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{$type}->{$vector_var_ids} = undef;
		 		}
			 	if ($hVarAnnex->{$var_id}->{'patients'}->{$pat_name}->{'he'} == 1){
					$hFork_catUsed->{"$subdata_filter;patients;$pat_name;global_categories;he"} = undef;
					$hFork->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{he}->{$vector_var_ids} = undef;
			 	}
			 	elsif ($hVarAnnex->{$var_id}->{'patients'}->{$pat_name}->{'ho'} == 1){
					$hFork_catUsed->{"$subdata_filter;patients;$pat_name;global_categories;ho"} = undef;
					$hFork->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{ho}->{$vector_var_ids} = undef;
			 	}
				# check infos Cosmic for Cancer projects
				if ($project->isSomaticStudy) {
					if ($variation->isCosmic()) {
						$hFork_catUsed->{"$subdata_filter;patients;$pat_name;global_categories;cosmic"} = undef;
						$hFork->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{cosmic}->{$vector_var_ids} = undef;
					}
					else {
						$hFork_catUsed->{"$subdata_filter;patients;$pat_name;global_categories;not_cosmic"} = undef;
						$hFork->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{not_cosmic}->{$vector_var_ids} = undef;
					}
				}
			}
			
			# consequences par genes (AVEC des genes)
			my @lVarGenes = @{$variation->getGenes()};
			if ( scalar(@lVarGenes) > 0 ) {
			 	foreach my $g (@lVarGenes) {
			 		my $id_gene = int($hGenesIds->{by_obj_id}->{$g->id()});
					$h_genes_ids->{$id_gene} = $g->id();
					$h_genes_ids->{$g->id()} = $id_gene;
					unless ($id_gene) {
						warn '  -> ERROR '.$g->id();
					}
			 		die unless ($id_gene);
	 	 			unless (exists $h_genes_infos->{$id_gene}) {
						$h_genes_infos->{$id_gene} = basic_gene_infos($g);
	 				}
		 	 		my $cons_text = $variation->variationType($g);
		 	 		my @gene_variation_type =  split(",", lc($cons_text));
		 	 		die()  if (scalar (@gene_variation_type) > 2);
		 	 		foreach my $typev (@gene_variation_type) {
			 			foreach my $pat_name (keys %{$hVarAnnex->{$var_id}->{'patients'}}) {
							$hFork_catUsed->{"$subdata_filter;patients;$pat_name;categories;$typev"} = undef;
							$hFork->{$subdata_filter}->{patients}->{$pat_name}->{categories}->{$typev}->{$vector_var_ids} = undef;
			 			}
						$hFork_catUsed->{"$subdata_filter;genes;$id_gene;categories;$typev"} = undef;
						$hFork->{$subdata_filter}->{genes}->{int($id_gene)}->{categories}->{$typev}->{$vector_var_ids} = undef;
		 	 		}
					$hFork_catUsed->{"$subdata_filter;genes;$id_gene;global_categories;$bank"} = undef;
					$hFork->{$subdata_filter}->{genes}->{int($id_gene)}->{global_categories}->{$bank}->{$vector_var_ids} = undef;
					foreach my $type (@lType) {
						$hFork_catUsed->{"$subdata_filter;genes;$id_gene;global_categories;$type"} = undef;
						$hFork->{$subdata_filter}->{genes}->{int($id_gene)}->{global_categories}->{$type}->{$vector_var_ids} = undef;
					
					}
					$hFork_catUsed->{"$subdata_filter;genes;$id_gene;global_categories;all"} = undef;
					$hFork->{$subdata_filter}->{genes}->{int($id_gene)}->{global_categories}->{all}->{$vector_var_ids} = undef;
	 				# calcul des scores polyphen et sift ainsi que l effectImpact (text et score) pour PolyDiag.
		 	 		my $key_polyphen = "polyphen".$variation->polyphenStatus($g);
					my $key_sift = "sift".$variation->siftStatus($g);
					my $key_polyphen_sift;
					if (($key_polyphen eq 'polyphen3') and ($key_sift eq 'sift2')) {
						$key_polyphen_sift = 'prediction_3';
					}
					elsif (($key_polyphen eq 'polyphen3') or ($key_sift eq 'sift2')) {
						$key_polyphen_sift = 'prediction_2';
					}
					elsif ($key_polyphen eq 'polyphen2') {
						$key_polyphen_sift = 'prediction_1';
					}
					else {
						$key_polyphen_sift = 'prediction_0';
					}
					$hFork_catUsed->{"$subdata_filter;genes;$id_gene;global_categories;$key_polyphen_sift"} = undef;
					$hFork->{$subdata_filter}->{genes}->{int($id_gene)}->{global_categories}->{$key_polyphen_sift}->{$vector_var_ids} = undef;
					
					my (@lEffectImpactText, @lEffectImpactScore);
					foreach my $t (@{$variation->getTranscripts()}) {
			 	 		my $effectImpactText = 'impact_text_'.$variation->effectImpact($t);
			 	 		my $effectImpactScore = 'impact_score_'.$himpact_sorted->{$variation->effectImpact($t)};
						push(@lEffectImpactText, $effectImpactText);
						push(@lEffectImpactScore, $effectImpactScore);
					}
		 			foreach my $pat_name (keys %{$hVarAnnex->{$var_id}->{'patients'}}) {
						$hFork_catUsed->{"$subdata_filter;patients;$pat_name;global_categories;$key_polyphen_sift"} = undef;
						$hFork->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{$key_polyphen_sift}->{$vector_var_ids} = undef;
			 			foreach my $effectImpactText (@lEffectImpactText) {
							$hFork_catUsed->{"$subdata_filter;patients;$pat_name;global_categories;$effectImpactText"} = undef;
							$hFork->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{$effectImpactText}->{$vector_var_ids} = undef;
			 			}
			 			foreach my $effectImpactScore (@lEffectImpactScore) {
							$hFork_catUsed->{"$subdata_filter;patients;$pat_name;global_categories;$effectImpactScore"} = undef;
							$hFork->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{$effectImpactScore}->{$vector_var_ids} = undef;
			 			}
		 			}
		 			$g = undef;
		 	 	}
		 	 	@lVarGenes = undef;
			}
		 	# consequences par genes (SANS gene)
			else {
				foreach my $pat_name (keys %{$hVarAnnex->{$var_id}->{'patients'}}) {
					$hFork_catUsed->{"$subdata_filter;patients;$pat_name;categories;intergenic"} = undef;
					$hFork->{$subdata_filter}->{patients}->{$pat_name}->{categories}->{intergenic}->{$vector_var_ids} = undef;
		 			
					$hFork_catUsed->{"$subdata_filter;patients;$pat_name;global_categories;polyphen0"} = undef;
					$hFork->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{polyphen0}->{$vector_var_ids} = undef;
		 			
					$hFork_catUsed->{"$subdata_filter;patients;$pat_name;global_categories;sift0"} = undef;
					$hFork->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{sift0}->{$vector_var_ids} = undef;
				}
	 			my ($intergenic_name, $id_gene, $hFork, $hFork_catUsed) = addVarIntergenic($hGenesIds, $vector_var_ids, $bank, \@lType, $variation->start(), $subdata_filter, $dir_temp, $hFork, $hFork_catUsed);
		 		unless ($id_gene) {
	 				warn "\n\n### ERROR\n\n";
	 				warn ref($variation).' -> '.$variation->id().' (start: '.$variation->start().' - end: '.$variation->end().')';
	 				warn '$vector_var_ids: '.$vector_var_ids;
	 				warn '$subdata_filter: '.$subdata_filter;
	 				warn '$bank: '.$bank;
	 				warn '@lType: '.join(', ', @lType);
	 				warn 'Positions start intergenic: ';
	 				foreach my $start (sort {$a <=> $b} keys %{$hGenesIds->{intergenic}}) {
	 					warn '  -> INTERGENIC start: '.$start.' - end: '.$hGenesIds->{intergenic}->{$start}->{end};
	 				}
	 				die;
		 		}
				$h_genes_ids->{$id_gene} = $intergenic_name;
				$h_genes_ids->{$intergenic_name} = $id_gene;
				$hset->{$subdata_filter}->{$intergenic_name} = int($id_gene);
			}
			$variation = undef;
	    }
		my $hRes;
		$hRes->{var_ids} = $h_var_ids;
		$hRes->{var_infos}->{$this_part} = $h_var_infos;
		$hRes->{genes_ids} = $h_genes_ids;
		$hRes->{genes_infos} = $h_genes_infos;
		$hRes->{annot}->{$this_part} = $hFork;
		$hRes->{cat_used}->{$this_part} = $hFork_catUsed;
		$hRes->{sub_part} = $this_part;
	    $hFork_catUsed = undef;
	    $hFork = undef;
		$pm->finish(0, $hRes);
	}
	$pm->wait_all_children();
	sleep(1);
	warn "\n" if ($verbose);
	my $obs_var = int($hConnect->{'var_ids'}->count_bulk($chr_name)) / 2;
	my $obs_genes = int($hConnect->{'genes_ids'}->count_bulk($chr_name)) / 2;
	my $exp_var = scalar(keys %{$hVarIds});
	my $exp_genes = scalar(keys %{$hGenesIds->{by_obj_id}});
	my $size_var   = $exp_var + 10;
	my $size_genes = $exp_genes + 10;
	if ($verbose) {
		warn 'Obs Nb Var: '.$obs_var;
		warn 'Exp Nb Var: '.$exp_var.' -> vector size: '.$size_var;
		warn "\n\nWARNING: expected $exp_var variants but found $obs_var variations...\n\n" unless ($obs_var == $exp_var);
		warn 'Obs Genes: '.$obs_genes;
		warn 'Exp Genes: '.$exp_genes.' -> vector size: '.$size_genes;
		warn "\n\nWARNING: expected $exp_genes genes but found $obs_genes genes...\n\n" unless ($obs_genes == $exp_genes);
	}
	closeNoSql($hConnect);
	$hConnect = undef;
	
	my $hVector;
	my ($hSet_all, $hCatUsed_all);
	($hVector->{$subdata_filter}, $hSet_all, $hCatUsed_all) = getHashBitVector_all($dir_temp, $chr_name, $size_var, $size_genes, $max);
	
	my $vTMP = Bit::Vector->new( $size_var );
	foreach my $patient (@{$project->getPatients()}) {
		my $pat_name = $patient->name();
		unless (exists $hVector->{$subdata_filter}->{patients}->{$pat_name}) {
			$hVector->{$subdata_filter}->{patients}->{$pat_name}->{categories}->{intergenic} = Bit::Vector->new( $size_var );
		}
		foreach my $cat_name ('substitution', 'insertion', 'deletion', 'ho', 'he') {
			unless (exists $hVector->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{$cat_name}) {
				$hVector->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{$cat_name} = Bit::Vector->new( $size_var );
			}
			else { $vTMP += $hVector->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{$cat_name}; }
		}
	}
	
	my ($is_ok) = checkResults($exp_var, $subdata_filter, $hVector, $hVarIds) if ($exp_var > 0);
	unless ($is_ok) {
		my $cache_dir = $project->getCacheBitVectorDir();
		warn "\n\nERROR in cache process... Supress $cache_dir. Die.\n\n";
		if (-d $cache_dir) {
			my $new_dir = $cache_dir;
			$new_dir =~ s/\/vector\//\/vector.pb_cache\//;
			`mv $cache_dir $new_dir`;
		}
		die();
	}
	mkdir $dir_temp.'../all_global/' unless (-d $dir_temp.'../all_global/');
	my $fileout = $dir_temp.'../all_global/'.$chr_name.'.freeze';
	store($hVector->{$subdata_filter}, $fileout);
	return ($hSet_all, $hCatUsed_all, $size_genes);
}

sub checkResults {
	my ($exp_var, $subdata_filter, $hVector, $hVarIds) = @_;
	my @lValues = sort {$a <=> $b} values %$hVarIds;
	my $min = $lValues[0];
	my $max = $lValues[-1];
	my $len = scalar(@lValues);
	unless ($min == 1) {
		warn "\n\nERROR: first vector_id expected is 1 and I found $min...\n\n";
		return;
	}
	unless ($len == $max) {
		warn "\n\nERROR: expected $max vector_id and I found $len vector_ids...\n\n";
		return;
	}
	
	print "Checking BitVector created:\n" if ($verbose);
	my $v_patients;
	foreach my $pat_name (keys %{$hVector->{$subdata_filter}->{patients}}) {
		foreach my $cat_name (keys %{$hVector->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}}) {
			unless ($v_patients) { $v_patients = $hVector->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{$cat_name}; }
			else { $v_patients += $hVector->{$subdata_filter}->{patients}->{$pat_name}->{global_categories}->{$cat_name}; }
		}
	}
	my @lVar = sort {$a <=> $b} @{transformBitVectorToList($v_patients)};
	my $nbvarInPatients = scalar(@lVar);
	if ($nbvarInPatients == $exp_var) { print "  -> from patients ($nbvarInPatients variants) ok !\n" if ($verbose); }
	else {
		warn "\n\nERROR: expected $exp_var var results but found $nbvarInPatients var from patients vectors...\n\n";
		return;
	}
	if ($lVar[0] == $min) { print "  -> from patients first vector_id $min found ok !\n" if ($verbose); }
	else {
		warn "\n\nERROR: vector_id $min not found from patients...\n\n";
		return;
	}
	if ($lVar[-1] == $max) { print "  -> from patients last vector_id $max found ok !\n" if ($verbose); }
	else {
		warn "\n\nERROR: vector_id $max not found from patients...\n\n";
		return;
	}
	@lVar = undef;
	if ($subdata_filter eq 'all') {
		my $v_genes = Bit::Vector->new( ($exp_var + 10) );
		foreach my $g_name (keys %{$hVector->{$subdata_filter}->{genes}}) {
			my $start = int($hVector->{$subdata_filter}->{genes}->{$g_name}->{start_subpart});
			foreach my $id (@{transformBitVectorToList($hVector->{$subdata_filter}->{genes}->{$g_name}->{global_categories}->{all})}) {
				$v_genes->Bit_On( ($id + $start) );
			}
		}
		my $nbvarInGenes = scalar(@{transformBitVectorToList($v_genes)});
		if ($nbvarInGenes == $exp_var) { print "  -> from genes ($nbvarInGenes variants) ok !\n" if ($verbose); }
		else {
			warn "\n\n[Subdata: $subdata_filter] ERROR: expected $exp_var va results but found $nbvarInGenes var from genes vectors...\n\n";
			return;
		}
	}
	return 1;
}

sub transformBitVectorToList {
	my $bitvector = shift;
	my @lIds;
	foreach my $id (split(',', $bitvector->to_Enum())) {
		my @lTmp = split('-', $id);
		if (scalar(@lTmp) > 1) {
			my $i = int($lTmp[0]);
			while ($i <= int($lTmp[1])) {
				push(@lIds, $i);
				$i++;
			}
		}
		else { push(@lIds, $id); }
	};
	return \@lIds;
}

sub addVarIntergenic {
	my ($hGenesIds, $id, $bank, $listType, $pos, $subdata_filter, $dir_temp, $hFork, $hFork_catUsed) = @_;
	my ($name, $id_gene_bitvector);
	foreach my $start (sort {$a <=> $b} keys %{$hGenesIds->{intergenic}}) {
		if ($start <= int($pos)) {
			if (int($pos) <= $hGenesIds->{intergenic}->{$start}->{end}) {
				$name = $hGenesIds->{intergenic}->{$start}->{id};
				$id_gene_bitvector = int($hGenesIds->{by_obj_id}->{$name});
			}
		}
		last if ($id_gene_bitvector);
	}
	$hFork_catUsed->{"$subdata_filter;genes;$id_gene_bitvector;global_categories;$bank"} = undef;
	$hFork->{$subdata_filter}->{genes}->{$id_gene_bitvector}->{global_categories}->{$bank}->{$id} = undef;
	foreach my $type (@$listType) {
		$hFork_catUsed->{"$subdata_filter;genes;$id_gene_bitvector;global_categories;$type"} = undef;
		$hFork->{$subdata_filter}->{genes}->{$id_gene_bitvector}->{global_categories}->{$type}->{$id} = undef;
	}
	$hFork_catUsed->{"$subdata_filter;genes;$id_gene_bitvector;categories;intergenic"} = undef;
	$hFork->{$subdata_filter}->{genes}->{$id_gene_bitvector}->{categories}->{intergenic}->{$id} = undef;
	$hFork_catUsed->{"$subdata_filter;genes;$id_gene_bitvector;global_categories;polyphen0"} = undef;
	$hFork->{$subdata_filter}->{genes}->{$id_gene_bitvector}->{global_categories}->{polyphen0}->{$id} = undef;
	$hFork_catUsed->{"$subdata_filter;genes;$id_gene_bitvector;global_categories;sift0"} = undef;
	$hFork->{$subdata_filter}->{genes}->{$id_gene_bitvector}->{global_categories}->{sift0}->{$id} = undef;
	$hFork_catUsed->{"$subdata_filter;genes;$id_gene_bitvector;global_categories;all"} = undef;
	$hFork->{$subdata_filter}->{genes}->{$id_gene_bitvector}->{global_categories}->{all}->{$id} = undef;
	return ($name, $id_gene_bitvector, $hFork, $hFork_catUsed);
}

sub resetHset {
	my @lKeys = keys %$hset;
	foreach my $key (@lKeys) { delete $hset->{$key}; }
	$hCatUsed = undef;
}

sub getHashReferencesCoord {
	my ($project, $chr, $this_limit) = @_;
	unless ($this_limit) {
		my $len = $chr->end() + 100;
		#$this_limit = int($len / (2 * $fork));
		$this_limit = int($len /  $fork);
	}
	my $hRefCoord;
	my $id = 1;
	my $start = 1;
	my $next  = $this_limit;
	my $end   = $next;
	while ($end <= $chr->end()) {
		$hRefCoord->{$id}->{start} = $start;
		$hRefCoord->{$id}->{end} = $end;
		$id++;
		$start = $end;
		$end  += $next;
	}
	$hRefCoord->{$id}->{start} = $start;
	$hRefCoord->{$id}->{end} = $chr->end();
	return ($hRefCoord, $this_limit);
}

sub basic_gene_infos {
	my $g = shift;
	my $hash;
	my $pos = $g->position($g->getChromosome);
	$hash->{start}  	= $pos->start();
	$hash->{start}  	= 1 if ($hash->{start} <= 1);
	$hash->{end}    	= $pos->end();
	$hash->{chromosome}	= $g->getChromosome->name();
	$hash->{strand} 	= $pos->strand();
	$hash->{name}   	= $g->name();
	$hash->{id}   		= $g->id();
	$hash->{cover} 		= 0; 
	if ($g->isGene()){
		$hash->{gene} = 1;
		$hash->{xref}      		= $g->external_name();
		$hash->{description} 	= $g->description();
		$hash->{transcripts} 	= join (";",map {$_->name} @{$g->getTranscripts()} );
	}
	return $hash;
}

sub IdsByRegion {
	my ($project_name,$chr_name,$start,$end,$key,$hpatients) = @_;
		my $buffer1 = new GBuffer;
		$fork = 1 unless $fork;
		my $project1 = $buffer1->newProject( -name => $project_name );
		my $chr = $project1->getChromosome($chr_name);
		#$buffer->dbh_reconnect();
		my $intspan_intergenic = Set::IntSpan::Fast::XS->new();
		my $hRes;
		$hRes->{variants} = [];
		$hRes->{genes} = [];
		$hRes->{patients} = {};
		$hRes->{dejavu} = {};
		$hRes->{dejavu_ho} = {};
		$hRes->{intergenic} =Set::IntSpan::Fast::XS->new();
		$hRes->{part} = $key;
#		die if ($key == 2);
		my $use_low_memory = 1;
		my $refObj = $project1->getChromosome($chr_name)->getReferences($start, $end)->[0];
		my @vars = @{$refObj->getStructuralVariations()};
		
		push(@vars, @{$refObj->getLargeDeletions()});
		if (scalar(@vars) == 0) {
			$hRes->{null} = 1;
			$project1->lite_deja_vu->close();
			#confess();
			return $hRes;
			#next;
		}
		my ($hGeneDone, $hIntergenic, $nbIntergenic);
		my $i = 1;
		foreach my $variation (sort {$a->start() <=> $b->start() } @vars) {
			$hRes->{large_deletions}->{$variation->id()} = undef if ($variation->isLargeDeletion);
			my $hPatients;
			
			my $aho= [];
			my $ap=[];
			my $hv;
			$hv->{id} = $variation->id();
			$hv->{start} = $variation->start();
			$hv->{end} = $variation->end();
			
			foreach my $pat (@{$variation->getPatients()}) {
				my $pn = $pat->name();
				my $hp;
				$hp->{name} = $pn;
				
				
				my $patient_id = $hpatients->{$pn};
				push(@$ap,$patient_id);
				push(@$aho,$patient_id) if ($variation->isHomozygote($pat));
				$hRes->{patients}->{$variation->id()}->{'vcf_infos'} = $variation->{check_id};
				$hp->{vcf_infos} =  $variation->{check_id};
				$hp->{he} =  $variation->{annex}->{$pat->id()}->{he};
				$hp->{ho} =  $variation->{annex}->{$pat->id()}->{ho};
				$hp->{he_ho_details} =  "he";
				$hp->{he_ho_details} =  "ho"  if ($variation->{annex}->{$pat->id()}->{ho} eq '1');
				$hp->{he_ho_details} .= ':'.$variation->{annex}->{$pat->id()}->{nb_all_ref}.':'.$variation->{annex}->{$pat->id()}->{nb_all_mut};
				
				$hRes->{patients}->{$variation->id()}->{'patients'}->{$pat->name()}->{'he'} = $variation->{annex}->{$pat->id()}->{he};
				$hRes->{patients}->{$variation->id()}->{'patients'}->{$pat->name()}->{'ho'} = $variation->{annex}->{$pat->id()}->{ho};
				$hRes->{patients}->{$variation->id()}->{'patients'}->{$pat->name()}->{'he_ho_details'} = 'he';
				$hRes->{patients}->{$variation->id()}->{'patients'}->{$pat->name()}->{'he_ho_details'} = 'ho' if ($variation->{annex}->{$pat->id()}->{ho} eq '1');
				$hRes->{patients}->{$variation->id()}->{'patients'}->{$pat->name()}->{'he_ho_details'} .= ':'.$variation->{annex}->{$pat->id()}->{nb_all_ref}.':'.$variation->{annex}->{$pat->id()}->{nb_all_mut};
				$hv->{patients}->{$pn} = $hp;
				#push(@{$hv->{patients}},$hp);
		
			}
			
			# line to prepare dejavu global;
			 $hv->{heho_string} =  join(",",sort{$a <=> $b} @$ap);
			 if (scalar(@$aho)) {
			 	$hv->{heho_string} =  $hv->{heho_string}."=".join(",",sort{$a <=> $b} @$aho)."=HO";
			}
			 push(@{$hRes->{variants}},$hv);
			my @lGenes = @{$variation->getGenes()};
			
		
		 	if ( scalar(@lGenes) > 0 ) {
		 		
		 		if ($nbIntergenic > 0) {
		 			my $t = $intspan_intergenic->as_string();
 					my @s = $t =~/(\d+)/g;
 					my $sstart = $s[0];
 					my $send = $s[-1];
 					my $hg ;
 					$hg->{id} = 'intergenic_'.$chr_name.'_'.$sstart.'_'.$send;
 					$hg->{start} = $sstart;
 					$hg->{end} = $send;
 					push(@{$hRes->{genes}},$hg);
					$intspan_intergenic->empty();
		 			$nbIntergenic = 0;
		 			delete $hIntergenic->{start};
		 			delete $hIntergenic->{end};
		 		}
		 		foreach my $g (@lGenes) {
	 				next if (exists $hGeneDone->{$g->id()});
	 				my $hg;
	 				$hg->{id} = $g->id;
	 				$hg->{start} = $g->start;
	 				$hg->{end} = $g->end;
	 				push(@{$hRes->{genes}},$hg);
					#$hRes->{genes}->{int($g->position($g->getChromosome)->start())}->{int($g->position($g->getChromosome)->end())}->{$g->id()} = undef;
					$hGeneDone->{$g->id()} = undef;
		 		}
		 	}
		 	else {
		 		
		 		$nbIntergenic++;
		 		my $var_start = $variation->start();
		 		my $sstart = $variation->start()-50000;
		 		$sstart =1 if $sstart <1;
		 		my $send = $variation->end+50000;
		 		$send = $chr->length() if $send>$chr->length() ;
		 		$hRes->{intergenic}->add_range($sstart,$send);
	 			if ($variation->isInsertion()) {
	 				my @lTmp = split('_', $variation->id());
	 				$var_start = int($variation->start()) - length($lTmp[-1]) + 1;
	 			}
	 			elsif ($variation->isDeletion()) {
	 				my @lTmp = split('_', $variation->id());
	 				$var_start = int($variation->start()) - length($lTmp[-2]) + 1;
	 			}
		 		unless (exists $hIntergenic->{start}) { $hIntergenic->{start} = $var_start; }
		 		if ($var_start < $hIntergenic->{start}) {
		 			$hIntergenic->{start} = $var_start;
		 		}
		 		$intspan_intergenic->add_range($var_start,$variation->end());
		 		$hIntergenic->{end} = $variation->end();
		 	}
		 	#$variation->frequency();
		 	
		 	#delete $variation->{project} ;
			#delete $variation ->{buffer};
			if ($i == 10000) {
				$variation->getProject->purge_memory($variation->end());
				foreach my $id (keys %{$variation->getProject->{objects}->{chromosomes}}) {
					delete $variation->getProject->{objects}->{chromosomes}->{$id} unless ($id eq $chr_name);
				}
				$i = 1;
			}
			$variation = undef;
			$i++;
		}
	#	warn $nbIntergenic;
#				warn "end  variations ***  $key ***";
 		if ($nbIntergenic > 0) {
 			#my @t = $hIntergenic->as_array();
 			my $t = $intspan_intergenic->as_string();
 			my @s = $t =~/(\d+)/g;
 			my $sstart = $s[0];
 			my $send = $s[-1];
 			my $hg ;
 			$hg->{id} = 'intergenic_'.$chr_name.'_'.$sstart.'_'.$send;
 			$hg->{start} = $sstart;
 			$hg->{end} = $send;
 			push(@{$hRes->{genes}},$hg);
 			$nbIntergenic = 0;
 			delete $hIntergenic->{start};
 			delete $hIntergenic->{end};
 		}
	
	return $hRes;
}


sub getIdsAndDejaVu {
	my ($dir_temp, $hRefCoord, $project, $chr, $fork) = @_;
	my $buffer = $project->buffer();
	my $project_name = $project->name();
		my $chr_name = $chr->id();



	mkdir $dir_temp unless (-d $dir_temp);
	my $hGenes = [];	
	my $txt_file = $dir_temp.'list.txt';
	my $intergenic_span =Set::IntSpan::Fast::XS->new();
	my $nosql_var       = GenBoNoSql->new(dir => $dir_temp.'tmp_var_ids', mode => 'c');
	
	my $nosql_patients  = GenBoNoSql->new(dir => $dir_temp.'tmp_var_patients', mode => 'c');
	
	my $nosql_dejavu    = GenBoNoSql->new(dir => $dir_temp.'../dejavu', mode => 'c');
	$nosql_dejavu->clear($chr_name);
	my $nosql_dejavu_ho = GenBoNoSql->new(dir => $dir_temp.'../dejavu_ho', mode => 'c');
	$nosql_dejavu_ho->clear($chr_name);
	#open database for dejavu global
	
	my $root_dir =$project->deja_vu_lite_dir()."/projects/";
	mkdir $root_dir unless -e $root_dir;

	my @patient_names =  sort{$a cmp $b} map{$_->name} @{$project->getPatients};
	my $hpatients;
	for (my $i=0;$i<@patient_names;$i++){
		$hpatients->{$patient_names[$i]} = $i;
	}
	

	
	my $pm = new Parallel::ForkManager($fork);
	
	my (@lCoord_var, @lCoord_genes, $hTmpPatients);
	my $nb_subpart = 1;
	my $tmp_i = 0;
	
	my $i_progress = 1;	
	my $nb_key = 0;
	my $nb_max = scalar(keys %$hRefCoord);

	my $pr;
	if ($verbose) {
		$pr = String::ProgressBar->new( max => $nb_max, info => 'Get all variants IDS' );
		$pr->write();
	}
	my $h_part_done;
	my $pm = new Parallel::ForkManager($fork, $dir_temp);
	my %amoi;
	my %vdejavu;
	unlink $txt_file if -e $txt_file;
	open(TXT,">".$txt_file);
	$pm->run_on_finish (
		sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hRes) = @_;
			unless (defined($hRes)){
				print qq|No message received from child process $exit_code $pid!\n|; 
				#die();
				return;
			}
			$h_part_done->{$hRes->{part}} = undef;
			delete $hRefCoord->{$hRes->{part}};
			$i_progress++;
			if ($verbose) {
					$pr->update($i_progress);
					$pr->write();
			}
		
			if (exists $hRes->{null}){
				$hRes = undef;
				return;
			}
			
					if (exists $hRes->{'large_deletions'}) {
						foreach my $id (keys %{$hRes->{'large_deletions'}}) {
							$hLargeDelIds->{$id} = undef;
						}
					}
				foreach my $hv (@{$hRes->{'variants'}}){
										print TXT  $hv->{start}."\t".$hv->{end}."\t".$hv->{id}."\t".$hv->{heho_string}."\t".encode_json($hv)."\n";
				}
				$intergenic_span = $intergenic_span->union($hRes->{intergenic});
				$nosql_patients->put($chr_name, $hRes->{part}, $hRes->{'patients'});
				push(@$hGenes,@{$hRes->{genes}}) if $hRes->{genes};
				$hRes = undef;
		}
	);
	
	
	my $toto;
	my $turn =0;
	
	#####
	#start fork annotation list var
	#####
	
	while ($turn <2){
	foreach my $key (sort {$a <=> $b} keys %$hRefCoord) {
		$nb_key++;
		
		my $pid = $pm->start and next;
		
		my $hRes =  IdsByRegion($project_name,$chr_name,$hRefCoord->{$key}->{start},$hRefCoord->{$key}->{end},$key,$hpatients);

		$pm->finish(0, $hRes);
	}
	
	$pm->wait_all_children();
	
	######
	# end fork annotation list var
	######
	
	
	my $exp = scalar(keys %$hRefCoord);
	if ($exp == 0) {
		last;
	}
	else {
			foreach my $part (keys %$hRefCoord) {
			unless (exists $h_part_done->{$part}) {
				my $this_start = $hRefCoord->{$part}->{start};
				my $this_end   = $hRefCoord->{$part}->{end};
				warn "  -> subpart $part (start:$this_start - end:$this_end) missing...\n";
			}
			}
		$turn ++;
	}
	
	}
	
	close (TXT);
	
	my $exp = scalar(keys %$hRefCoord);
	
	
	if ($exp ne 0){
			foreach my $part (keys %$hRefCoord) {
			unless (exists $h_part_done->{$part}) {
				my $this_start = $hRefCoord->{$part}->{start};
				my $this_end   = $hRefCoord->{$part}->{end};
				warn "  -> subpart $part (start:$this_start - end:$this_end) missing...\n";
			}
			die();
	}
	} 
	
	
	
	if (scalar(keys %{$hTmpPatients}) > 0) {
		$nosql_patients->put($chr_name, $nb_subpart, $hTmpPatients);
		$hTmpPatients = undef;
	}

my $real_intergenic_span = $chr->intronic_intspan->intersection( $intergenic_span);
my $iter = $real_intergenic_span->iterate_runs();
my $array_intergenic_intspan = Array::IntSpan->new();
my $nb1 =0;
 while (my ( $from, $to ) = $iter->()) {
 	my $id = 'intergenic_'.$chr_name.'_'.$from.'_'.$to;
 	$array_intergenic_intspan->set_range($from,$to,$id);
 }
 

	my $hVarIds;
	my $vector_id_var  = 1;
	my $hintergenic;
	my $no2 = $project->getChromosome($chr_name)->get_lmdb_variations("c");
	
	my $nolmdb_patients;
	
#open database for dejavu global
	

	my @patient_names =  sort{$a cmp $b} map{$_->name} @{$project->getPatients};
	my $hpatients;
	for (my $i=0;$i<@patient_names;$i++){
		my $pname = $patient_names[$i] ;
		my $patient = $project->getPatient($pname);
		$nolmdb_patients->{$patient_names[$i]}->{intspan}->{all}   = Set::IntSpan::Fast::XS->new();#GenBoNoSqlLmdb->new(dir=>$dir_out,mode=>"c",name=>$patient_names[$i]);
		$nolmdb_patients->{$patient_names[$i]}->{intspan}->{he}   = Set::IntSpan::Fast::XS->new();
		$nolmdb_patients->{$patient_names[$i]}->{intspan}->{ho}   = Set::IntSpan::Fast::XS->new();
		#$nolmdb_patients->{$patient_names[$i]}->{lmdb} = $project->getChromosome($chr_name)->get_lmdb_patients_variations("c",$patient );
		
	}
#	$no3->close();
	
	#my $no = GenBoNoSql->new(dir=>$root_dir,mode=>"c");
	
my $nn = `wc  -l $txt_file`;
my $hh = {};
	open(SORT," sort -n  -k1,1 -k2,2 $txt_file |  ");
	while( my$line = <SORT>){
		chomp($line);
		my ($start,$end,$var_id,$heho_string,$json) = split(" ",$line);
			next if exists $hVarIds->{$var_id};
			my $hv =  decode_json($json);
			  my $index_lmdb = $no2->put($var_id, $hv);
			  
			foreach my $p (keys %{$hv->{patients}}){
				$nolmdb_patients->{$p}->{intspan}->{all}->add($index_lmdb);
				my $type = "he";
				$type = "ho" if $hv->{patients}->{$p}->{ho} eq 1; #he or ho
				$nolmdb_patients->{$p}->{intspan}->{$type}->add($index_lmdb);
				#$nolmdb_patients->{$p}->{lmdb}->put($index_lmdb,$type);
			}
			$hv->{id} = $vector_id_var;
			#$hv->{patients} = decode_json($json);
			$heho_string =~s/=/ /g;
			#$hv->{heho_string} = $heho_string;
		  
		     $hh->{$var_id} = $heho_string;
		     $hVarIds->{$var_id}  = $vector_id_var;
		     #my $hv = decode_json($json);
		     #warn Dumper $hv;
#		     foreach my $patient (@{$hv->{patients}}){
#		     		$nolmdb_patients->{$patient_names[$i]}->put()
#		     }
		     
			$vector_id_var ++;
	}
my $no4 = $project->getChromosome($chr_name)->get_lmdb_patients("c");
	foreach my $p (keys %{$nolmdb_patients}){
		$no4->put($p,$nolmdb_patients->{$p});
		#warn $nolmdb_patients->{$p}->as_string();
	}
	
	
close SORT;

store($hh, $project->lmdb_cache_dir."/$chr_name.dv.freeze") if $hh;

###
# add verification on size here 
####
die() if $no2->nb_keys+1 ne $vector_id_var;
$no2->close();
$nb1 = scalar(keys %$hintergenic);

	my $nb2=0;
	my $hGenesIds;
	my $vector_id_gene = 1;	
	foreach my $gene (sort {$a->{start} <=>  $b->{start} or $a->{end} <=>  $b->{end}} @$hGenes){
		my $gene_id = $gene->{id};
		next if exists $hGenesIds->{by_obj_id}->{$gene_id};
		$hGenesIds->{by_obj_id}->{$gene_id} = $vector_id_gene;
			my $start = $gene->{start};
			my $end = $gene->{end};
			if ($gene_id =~ /intergenic/) {
					$nb2++;
					$hGenesIds->{intergenic}->{$start}->{end} = $end;
					$hGenesIds->{intergenic}->{$start}->{id}  = $gene_id;
				}
				$vector_id_gene++;
	}
	
	$nosql_var->close();
	$nosql_patients->close();
	$nosql_dejavu->close();
	$nosql_dejavu_ho->close();
	
	warn "Launching $exp subparts and received $nb_max subparts\n" if ($verbose);
	return ($hVarIds, $hGenesIds,$nb_max);
}

sub cache_lite_for_dejavu {
	my ($project_name,$fork) = @_;

	#my $root_dir = "/data-isilon/dejavu/projects/";
	warn "*** CAche For Deja Vu *****";
	my $buffer1 = new GBuffer;
	my $project = $buffer1->newProject( -name => $project_name, -verbose =>1 );
	my $root_dir =$project->deja_vu_lite_dir()."/projects/";

		mkdir $root_dir unless -e $root_dir;
	unlink $root_dir."/".$project_name.".lite" if -e $root_dir."/".$project_name.".lite";
	my @chr_names = map{$_->name} @{$project->getChromosomes};
	#my $pm = new Parallel::ForkManager(1);
	my @patient_names =  sort{$a cmp $b} map{$_->name} @{$project->getPatients};
	my $dir_out = 	 $project->getCacheBitVectorDir()."/lmdb_cache";
	my $hpatients;
	for (my $i=0;$i<@patient_names;$i++){
		$hpatients->{$patient_names[$i]} = $i;
	}
	my $no = GenBoNoSql->new(dir=>$root_dir,mode=>"c");
	$no->put($project_name,"patients",$hpatients);
	foreach my $chr (@{$project->getChromosomes}) {
		my $fileout = $dir_out."/".$chr->name.".dv.freeze";; 
		warn "miss $fileout " unless -e $fileout;
		next unless -e $fileout;
		my $h =  retrieve $fileout;
		$no->put($project_name,$chr->name,$h);
	}
	

}



sub updateDejaVuPn {
	my ($project, $chr_name, $fork) = @_;
	my $buffer = $project->buffer;
	my $pm = new Parallel::ForkManager($fork);
	my $nosql_dejavu    = GenBoNoSql->new(dir => $project->getCacheBitVectorDir().'dejavu', mode => 'c');
	$nosql_dejavu->clear($chr_name);
	my $nosql_dejavu_ho = GenBoNoSql->new(dir => $project->getCacheBitVectorDir().'dejavu_ho', mode => 'c');
	$nosql_dejavu_ho->clear($chr_name);
		my $tranche;
	$pm->run_on_finish (
		sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hRes) = @_;

			unless (defined($hRes)){
				print qq|No message received from child process $exit_code $pid!\n|; 
				#die();
				return;
			}
				my $part = $hRes->{part};
				delete $tranche->{$hRes->{part}};
				$nosql_dejavu->put_bulk($chr_name, $hRes->{'dejavu'});
				$nosql_dejavu_ho->put_bulk($chr_name, $hRes->{'dejavu_ho'});
		
		}
	);
	
	
	
	my $dir_dejavu = $buffer->config->{public_data}->{HG19}."/".$buffer->config->{kyoto}->{deja_vu}."/lite/lmdb_chr/";
	system ("mkdir -p $dir_dejavu && chmod a+rwx $dir_dejavu ") unless -e $dir_dejavu ;
	my $no = $project->getChromosome($chr_name)->get_lmdb_variations("r");
	my $ranges = $no->ranges($fork);
	
	$no->close();
	
	my $dir_out = $project->lmdb_cache_variations_dir();

	#my $dir_dejavu = $project->deja_vu_lite_dir();
		
	$buffer->dbh_deconnect();
	my $turn =0;
	my $nbt = 0;
while($turn < 2){
	foreach my $r (@$ranges){
			$tranche->{$r->[0]."-".$r->[1]} ++;
			$nbt ++;
			my $pid = $pm->start and next;
			my $hRes;
			warn $r->[0]."-".$r->[1] if $turn > 0;
			
			my $no = GenBoNoSqlLmdb->new(dir=>$dir_out,mode=>"r",is_index=>1,name=>$chr_name);
			my $no2 = GenBoNoSqlLmdb->new(dir=>$dir_dejavu,mode=>"r",name=>$chr_name);
			#my $no2 = GenBoNoSqlDejaVu->new(dir=>$dir_dejavu,mode=>"r");
			my $cursor =  $no->cursor($r->[0],$r->[1]);
			while (my $vid = $cursor->next_key) {
				my $nb =0;
				my $nb1 =0;		
 			 	$nb += $no2->get($vid."_ho");
 			 
 			 	$nb1 += $no2->get($vid);
 				$hRes->{dejavu}->{$vid} = $nb1;
				$hRes->{dejavu_ho}->{$vid} = $nb;
				$hRes->{part} = $r->[0]."-".$r->[1];
			}
			#warn $r->[0]."-".$r->[1];
		#$project->lite_deja_vu->close();
		#$buffer->dbh_deconnect();
		$no->close();
		$no2->close();
		#warn "\t end ".$r->[0]."-".$r->[1];
		$pm->finish(0,$hRes);
		 }
		 $pm->wait_all_children();
			my $exp = scalar(keys %$tranche);
		if ($exp == 0) {
			last;
		}
		else {
			 $ranges =[];
			foreach my $part (keys %$tranche) {
				my ($start,$end) = split("-",$part);
				push(@$ranges,[$start,$end]);
				warn "POUET  on ". $part;
			}
			$tranche={};
		$turn ++;
	}
}
	
		 die() if scalar(keys %$tranche);
		 
	}



sub updateDejaVu {
	my ($project, $chr_name, $fork, $limit) = @_;
	my $buffer = $project->buffer();
	my $chr = $project->getChromosome($chr_name);
	foreach my $patient (@{$project->getPatients()}) {
		$patient->alignmentMethods();
		$patient->callingMethods();
	}
	unless ($limit) {
		my $len = $chr->end() + 100;
		$limit = int($len / (3 * $fork));
	}
	warn "\n### DEJAVU: update\n" if ($verbose);
	warn strftime "Start at %H:%M:%S\n", localtime if ($verbose);
	my $dir_temp = $project->getCacheBitVectorDir().'/tmp_dejavu_chr'.$chr_name.'/';
	mkdir $dir_temp unless (-d $dir_temp);
	my ($hRefCoord, $limit) = getHashReferencesCoord($project, $chr, $limit);
	my $file1 = $dir_temp.'../dejavu/'.$chr_name.'.lite';
	my $file2 = $dir_temp.'../dejavu_ho/'.$chr_name.'.lite';
	unless (-d $dir_temp.'../dejavu/') { mkdir($dir_temp.'../dejavu/'); }
	unless (-d $dir_temp.'../dejavu_ho/') { mkdir($dir_temp.'../dejavu_ho/'); }
	if (-e $file1) {
		my $cmd = 'rm '.$file1;
		`$cmd`;
	}
	if (-e $file2) {
		my $cmd = 'rm '.$file2;
		`$cmd`;
	}
	my $nosql_dejavu    = GenBoNoSql->new(dir => $dir_temp.'../dejavu', mode => 'c');
	my $nosql_dejavu_ho = GenBoNoSql->new(dir => $dir_temp.'../dejavu_ho', mode => 'c');
	my $i_progress = 1;				
	my $nb_key = 0;
	my $nb_max = scalar(keys %$hRefCoord);
	my $pr;
	if ($verbose) {
		$pr = String::ProgressBar->new( max => $nb_max, info => 'Update Deja Vu' );
		$pr->write();
	}
	my $pm = new Parallel::ForkManager($fork, $dir_temp);
	$pm->run_on_finish (
		sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hRes) = @_;
			if (defined($hRes)) {
				$nosql_dejavu->put_bulk($chr_name, $hRes->{'dejavu'});
				$nosql_dejavu_ho->put_bulk($chr_name, $hRes->{'dejavu_ho'});
				$hRes = undef;
				$i_progress++;
				if ($verbose) {
					$pr->update($i_progress);
					$pr->write();
				}
			}
			else { print qq|No message received from child process $pid!\n|; }
		}
	);
	
	foreach my $key (sort {$a <=> $b} keys %$hRefCoord) {
		$nb_key++;
		my $pid = $pm->start and next;
		$buffer->dbh_reconnect();
		my $hRes;
		$hRes->{dejavu} = {};
		$hRes->{dejavu_ho} = {};
		my $start = int($hRefCoord->{$key}->{start});
		my $end   = int($hRefCoord->{$key}->{end});
		my $refObj = $chr->getReferences($start, $end)->[0];
		my @vars = sort {$a->start <=> $b->start } @{$refObj->getStructuralVariations()};
		my $i = 1;
		foreach my $variation (sort {$a->start() <=> $b->start() } @vars) {
			$hRes->{dejavu}->{$variation->id()} = $variation->nb_dejavu();
			$hRes->{dejavu_ho}->{$variation->id()} = $variation->nb_dejavu_ho();
			if ($i == 10000) {
				$variation->getProject->purge_memory($variation->end());
				foreach my $id (keys %{$variation->getProject->{objects}->{chromosomes}}) {
					delete $variation->getProject->{objects}->{chromosomes}->{$id} unless ($id eq $chr_name);
				}
				$i = 1;
			}
			$variation = undef;
			$i++;
		}
		$project->lite_deja_vu->close();
		$buffer->dbh_deconnect();
		$pm->finish(0, $hRes);
	}
	$pm->wait_all_children();
	$nosql_dejavu->close();
	$nosql_dejavu_ho->close();
	`rm -r $dir_temp`;
	
	my $freeze_infos = $project->getCacheBitVectorDir().'/global_infos.freeze';
	return unless -e $freeze_infos;
	my $hGlobalInfos = retrieve $freeze_infos;
	$hGlobalInfos->{analyse}->{cache}->{dejavu} = strftime '%Y-%m-%d', localtime;
	`rm $freeze_infos` if (-e $freeze_infos);
	store($hGlobalInfos, $freeze_infos);
	`chmod 777 $freeze_infos`;
	warn strftime "\nEnd at %H:%M:%S\n", localtime if ($verbose);
}

1;