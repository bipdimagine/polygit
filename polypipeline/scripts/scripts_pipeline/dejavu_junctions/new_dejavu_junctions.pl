#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use JSON;
use Data::Dumper;
use CGI qw/:standard :html3/;
use Set::IntSpan::Fast::XS;
use Set::IntervalTree;
use List::Util qw[min max];
use FindBin qw($Bin);
use Storable qw(store retrieve freeze);
use Number::Format qw(:subs);
use File::Basename;
use Getopt::Long;
use JSON;

use GBuffer;
use GenBoProject;
use GenBoCache;
use GenBoNoSqlIntervalTree;


#######################################################################
#  1) liste  l'ensemble des fichiers  project.dejavu
#  2) retrieve chaque table de hash correspondante et l'integre dans une table de hash globale 
#  3) pour tout les project_dejavu plus r??cents que le dejavuglobal : retrieve la hash et l'integre dans la hash globale
#  4) freeze de la table de hash gllobale
########################################################################

my $cgi = new CGI;
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => 'NGS2022_5612' );
$project->getChromosomes();

my $nbok;
my $nbSV;
my $halldejavu;
my $htree_dejavu;
my $ids;
my $total;
my $xx =0;

my $fork = 1;
my $release = 'HG19';

GetOptions(
	'release=s' => \$release,
	'fork=s'    => \$fork,
);

my $dv_dir_path = $project->DejaVuJunction_path();
if ($release =~ /HG19/) { $release = 'HG19'; }
else { $dv_dir_path =~ s/HG19/$release/; }

my $dv_dir_projects_path = $dv_dir_path.'/projects/';
print "\n# Checking RNA-Projects Junctions";

my $test = 0;
my ($hash_projects, $hash_project_phenotypes);
if (-d $dv_dir_projects_path) {
	opendir DIR, $dv_dir_projects_path;
	my @dir = readdir(DIR);
	close DIR;
	
	foreach my $project_name (@dir) {
#		next if $test == 10;
		next if not $project_name =~ /NGS/;
		$hash_projects->{$project_name} = undef;
#		$test++;
		my $b1 = GBuffer->new;
		my $p1 = $b1->newProject( -name => $project_name );
		if ($p1->is_human_genome()) {
			foreach my $pheno_name (@{$p1->phenotypes()}) {
				$pheno_name =~ s/ /_/g;
				$hash_project_phenotypes->{lc($pheno_name)}->{$project_name} = undef;
			}
		}
	}
}

#warn Dumper $hash_project_phenotypes;
#die;
print " -> Done!\n";
print " -> Found ".scalar(keys %{$hash_projects})." projects\n\n";


#foreach my $this_project_name (keys %$hash_projects) {
#	my $dir_dv_proj = $dv_dir_path.'/projects/'.$this_project_name.'/';
#	my $file = $dir_dv_proj.'/'.$this_project_name.'.json';
#	my $file_gz = $dir_dv_proj.'/'.$this_project_name.'.json.gz';
#	if (-e $file_gz) {
#		my $cmd_zip = "gunzip $dir_dv_proj/*.gz";
#		my $cmd_rm_tbi = "rm $dir_dv_proj/*.tbi";
#		`$cmd_zip`;
#		`$cmd_rm_tbi`;
#	}
#}


my $h_not_found;
foreach my $this_project_name (keys %$hash_projects) {
	my $dir_dv_proj = $dv_dir_path.'/projects/'.$this_project_name.'/';
	my $file = $dir_dv_proj.'/'.$this_project_name.'.json';
	my $file_gz = $dir_dv_proj.'/'.$this_project_name.'.json.gz';
	unless (-d $dir_dv_proj) {
		warn "WARN DejaVu $this_project_name not done yet....";
		$h_not_found->{$this_project_name} = undef;
	}
	if (not -e $file and not -e $file_gz) {
		warn "WARN DejaVu $this_project_name not done yet....";
		$h_not_found->{$this_project_name} = undef;
	}
}

my $dir = $dv_dir_path;
unless (-d $dir){
	`mkdir $dir`;
	`chmod 777 $dir`;
}

my $dir_canoniques = $dv_dir_path.'/canoniques/';
unless (-d $dir_canoniques){
	`mkdir $dir_canoniques`;
	`chmod 777 $dir_canoniques`;
}

my @lTypes;
push(@lTypes, 'all');
push(@lTypes, 'canoniques');
#foreach my $pheno_name (keys %{$hash_project_phenotypes}) {
#	push(@lTypes, $pheno_name);
#}

foreach my $type (@lTypes) {
	if ($type eq 'all') { insert_all_junctions(); }
#	elsif ($type eq 'canoniques') { insert_all_junctions_canoniques(); }
#	else { insert_all_junctions_type($type); }
}


print "\n\nALL DONE!\n\n";


sub insert_all_junctions {
	my $h_all_junctions;
	foreach my $this_project_name (keys %$hash_projects) {
		
		#next if $this_project_name ne 'NGS2021_4218';
		
		print "-> $this_project_name ";
		my $dir_dv_proj = $dv_dir_path.'/projects/'.$this_project_name.'/';
		my $file = $dir_dv_proj.'/'.$this_project_name.'.json';
		my $file_gz = $dir_dv_proj.'/'.$this_project_name.'.json.gz';
		my $file_canonique = $dir_dv_proj.'/'.$this_project_name.'.canoniques.json';
		my $file_canonique_gz = $dir_dv_proj.'/'.$this_project_name.'.canoniques.json.gz';
		next if exists $h_not_found->{$this_project_name};
		my ($h_proj_junctions, $h_proj_junctions_canonique);
		
		
		if (-e $file) {
			open (FILE, $file);
			my $json = <FILE>;
			$h_proj_junctions = decode_json $json;
			close (FILE);
		}
		elsif (-e $file_gz) {
			open (FILE, "gunzip -c $file_gz |");
			my $json = <FILE>;
			$h_proj_junctions = decode_json $json;
			close (FILE);
		}
		else {
			die("\n\nNot exists: $file  - DIE\n");
		}
		
		my $hres;
		foreach my $chr_id (1..22, 'X', 'Y', 'MT') {
			my $todo;
			if (exists $h_proj_junctions->{$chr_id}) {
				foreach my $junction_id (keys %{$h_proj_junctions->{$chr_id}}) {
					my $start = $h_proj_junctions->{$chr_id}->{$junction_id}->{start};
					my $end = $h_proj_junctions->{$chr_id}->{$junction_id}->{end};
					foreach my $pat_name (keys %{$h_proj_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$this_project_name}}) {
						$hres->{$chr_id}->{$start}->{$end}->{$pat_name} = $h_proj_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$this_project_name}->{$pat_name};
						$hres->{$chr_id}->{$start}->{$end}->{$pat_name}->{junction_id} = $junction_id;
					}
					my $new_id = $h_proj_junctions->{$chr_id}->{$junction_id}->{start}.'-'.$h_proj_junctions->{$chr_id}->{$junction_id}->{end};
					if (exists $h_all_junctions->{$chr_id}->{$new_id}) {
						$h_all_junctions->{$chr_id}->{$new_id}->{dejavu}->{$this_project_name} = $h_proj_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$this_project_name};
					}
					else {
						$h_all_junctions->{$chr_id}->{$new_id} = $h_proj_junctions->{$chr_id}->{$junction_id};
					}
				}
			}
		}
		$h_proj_junctions = undef;
			
		
		if (-e $file_canonique) {
			open (FILE, $file_canonique);
			my $json = <FILE>;
			eval { $h_proj_junctions_canonique = decode_json $json; };
			if ($@) { $h_proj_junctions_canonique = {}; }
			close (FILE);
		}
		elsif (-e $file_canonique_gz) {
			open (FILE, "gunzip -c $file_gz |");
			my $json = <FILE>;
			$h_proj_junctions_canonique = decode_json $json;
			close (FILE);
		}	
		foreach my $chr_id (1..22, 'X', 'Y', 'MT') {
			if (exists $h_proj_junctions_canonique->{$chr_id}) {
				foreach my $junction_id (keys %{$h_proj_junctions_canonique->{$chr_id}}) {
					my $start = $h_proj_junctions_canonique->{$chr_id}->{$junction_id}->{start};
					my $end = $h_proj_junctions_canonique->{$chr_id}->{$junction_id}->{end};
					foreach my $pat_name (keys %{$h_proj_junctions_canonique->{$chr_id}->{$junction_id}->{dejavu}->{$this_project_name}}) {
						$hres->{$chr_id}->{$start}->{$end}->{$pat_name} = $h_proj_junctions_canonique->{$chr_id}->{$junction_id}->{dejavu}->{$this_project_name}->{$pat_name};
						$hres->{$chr_id}->{$start}->{$end}->{$pat_name}->{junction_id} = $junction_id;
						my $a = $hres->{$chr_id}->{$start}->{$end}->{$pat_name}->{count_junctions};
						my $b = $hres->{$chr_id}->{$start}->{$end}->{$pat_name}->{count_normal};
						if (($a == 0 and $b > 0) or $a == $b) {
							$hres->{$chr_id}->{$start}->{$end}->{$pat_name}->{count_junctions} = 0;
						}
						else {
							$hres->{$chr_id}->{$start}->{$end}->{$pat_name}->{count_junctions} = 0;
							$hres->{$chr_id}->{$start}->{$end}->{$pat_name}->{count_normal} = $a;
						}
						$hres->{$chr_id}->{$start}->{$end}->{$pat_name}->{score} = 0;
					}
					
					my $new_id = $h_proj_junctions_canonique->{$chr_id}->{$junction_id}->{start}.'-'.$h_proj_junctions_canonique->{$chr_id}->{$junction_id}->{end};
					if (exists $h_all_junctions->{$chr_id}->{$new_id}) {
						$h_all_junctions->{$chr_id}->{$new_id}->{dejavu}->{$this_project_name} = $h_proj_junctions_canonique->{$chr_id}->{$junction_id}->{dejavu}->{$this_project_name};
					}
					else {
						$h_all_junctions->{$chr_id}->{$new_id} = $h_proj_junctions_canonique->{$chr_id}->{$junction_id};
					}
				}
			}
		}
		$h_proj_junctions_canonique = undef;
			
		my $file_tab = $dir_dv_proj.'/'.$this_project_name.'.tab';
		open (FILE, ">$file_tab");
		foreach my $chr_id (1..22, 'X', 'Y', 'MT') {
			foreach my $start (sort {$a <=> $b} keys %{$hres->{$chr_id}}) {
				foreach my $end (sort {$a <=> $b} keys %{$hres->{$chr_id}->{$start}}) {
					foreach my $pat_name (sort keys %{$hres->{$chr_id}->{$start}->{$end}}) {
						my $junction_id = $hres->{$chr_id}->{$start}->{$end}->{$pat_name}->{junction_id};
						my $type = $hres->{$chr_id}->{$start}->{$end}->{$pat_name}->{type};
						my ($count_junctions, $count_normal, $ratio);
						$count_junctions = $hres->{$chr_id}->{$start}->{$end}->{$pat_name}->{count_junctions};
						$count_normal = $hres->{$chr_id}->{$start}->{$end}->{$pat_name}->{count_normal};
						$ratio = $hres->{$chr_id}->{$start}->{$end}->{$pat_name}->{score};
						my @line;
						push(@line, $chr_id);
						push(@line, $start);
						push(@line, $end);
						push(@line, $pat_name);
						push(@line, $type);
						push(@line, $count_junctions);
						push(@line, $count_normal);
						push(@line, $ratio);
						push(@line, $junction_id);
						print FILE join("\t", @line)."\n";
					}
				}
			}
		}
		close (FILE);
		
		$hres = undef;
		
		
		if (-e "$file_tab.gz") {
			`rm $file_tab.gz`;
		}
		my $cmd1 = "bgzip $file_tab";
		`$cmd1`;
		my $cmd2 = "tabix -p bed $file_tab.gz";
		`$cmd2`;
		if (not -e $file_gz) {
			my $cmd3 = "bgzip $file";
			`$cmd3`;
		}
		if (not -e $file_canonique_gz) {
			my $cmd4 = "bgzip $file_canonique";
			`$cmd4`;
		}
		
		print "-> ok\n";
	}

#	my $nodejavu = GenBoNoSqlDejaVuJunctions->new( dir => $dir, mode => "c", is_compress => 1 );
#	insert_in_dejavu_jonctions($nodejavu, $h_all_junctions);
#	$nodejavu->close();
	print "# CREATED TABIX Junctions -> DONE!\n";
	
	my $h_all_junctions_resume;
	foreach my $chr_id (keys %{$h_all_junctions}) {
		foreach my $junction_id (keys %{$h_all_junctions->{$chr_id}}) {
			
			my $h_resume;
			$h_resume->{start} = $h_all_junctions->{$chr_id}->{$junction_id}->{start};
			$h_resume->{end} = $h_all_junctions->{$chr_id}->{$junction_id}->{end};
			
			$h_resume->{pat} = 0;
			$h_resume->{pat_r10} = 0;
			$h_resume->{pat_r20} = 0;
			my (@lDetails,@lDetailsR10, @lDetailsR20);
			foreach my $proj (keys %{$h_all_junctions->{$chr_id}->{$junction_id}->{dejavu}}) {
				foreach my $pat_name (keys %{$h_all_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$proj}}) {
					$h_resume->{pat}++;
					my $new = $h_all_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$proj}->{$pat_name}->{count_junctions};
					my $normal = $h_all_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$proj}->{$pat_name}->{count_normal};
					my $dp = ($new + $normal);
					my $ratio = 0;
					$ratio = ($new / $dp) * 100 if $new > 0;
					
					my $patname_text = $proj.'_'.$pat_name;
					$patname_text =~ s/NGS20//;
					
					if ($ratio >= 10)  {
						$h_resume->{pat_r10}++;
						push(@lDetailsR10, $patname_text);
					}
					if ($ratio >= 20)  {
						$h_resume->{pat_r20}++;
						push(@lDetailsR20, $patname_text);
					}
					push(@lDetails, $patname_text);
				}
			}
			$h_resume->{details} = join(';', @lDetails);
			$h_resume->{details_r10} = join(';', @lDetailsR10);
			$h_resume->{details_r20} = join(';', @lDetailsR20);
			$h_all_junctions_resume->{$chr_id}->{$junction_id} = $h_resume;
		}
	}
	
	my $nodejavu_resume = GenBoNoSqlDejaVuJunctionsResume->new( dir => $dir, mode => "c", is_compress => 1 );
	insert_in_dejavu_jonctions_resume($nodejavu_resume, $h_all_junctions_resume);
	$nodejavu_resume->close();
	print "# INSERT in DejaVuLMDB Junctions Resume -> DONE!\n";
}

#sub insert_all_junctions_canoniques {
#	my $h_all_junctions_canoniques;
#	foreach my $this_project_name (keys %$hash_projects) {
#		my $dir_dv_proj = $dv_dir_path.'/projects/'.$this_project_name.'/';
#		my $file = $dir_dv_proj.'/'.$this_project_name.'.canoniques.json';
#		next if exists $h_not_found->{$this_project_name};
#		die("\n\nNot exists: $file  - DIE\n") if (not -e $file);
#		open (FILE, $file);
#		my $json = <FILE>;
#		my $h_proj_junctions_canoniques; 
#		$h_proj_junctions_canoniques = decode_json $json if ($json);
#		close (FILE);
#		next if (not $h_proj_junctions_canoniques); 
#		foreach my $chr_id (keys %{$h_proj_junctions_canoniques}) {
#			foreach my $junction_id (keys %{$h_proj_junctions_canoniques->{$chr_id}}) {
#				my $jid = $junction_id;
#				$jid =~ s/_junction//;
#				if (exists $h_all_junctions_canoniques->{$chr_id}->{$junction_id}) {
#					$h_all_junctions_canoniques->{$chr_id}->{$jid}->{dejavu}->{$this_project_name} = $h_proj_junctions_canoniques->{$chr_id}->{$junction_id}->{dejavu}->{$this_project_name};
#				}
#				else {
#					$h_all_junctions_canoniques->{$chr_id}->{$jid} = $h_proj_junctions_canoniques->{$chr_id}->{$junction_id};
#				}
#			}	
#		}
#	}
#	
#	my $nodejavu_resume = GenBoNoSqlDejaVuJunctionsResume->new( dir => $dir, mode => "w", is_compress => 1 );
#	insert_in_dejavu_jonctions_resume($nodejavu_resume, $h_all_junctions_canoniques);
#	$nodejavu_resume->close();
#	print "# INSERT in DejaVuLMDB Junctions canoniques -> DONE!\n";
#}

#sub insert_all_junctions_type {
#	my ($pheno_name) = @_;
#	my $h_all_junctions_pheno;
#	foreach my $this_project_name (keys %{$hash_project_phenotypes->{$pheno_name}}) {
#		my $dir_dv_proj = $dv_dir_path.'/projects/'.$this_project_name.'/';
#		my $file = $dir_dv_proj.'/'.$this_project_name.'.json';
#		next if exists $h_not_found->{$this_project_name};
#		die("\n\nNot exists: $file  - DIE\n") if (not -e $file);
#		open (FILE, $file);
#		my $json = <FILE>;
#		my $h_proj_junctions = decode_json $json;
#		close (FILE);
#		foreach my $chr_id (keys %{$h_proj_junctions}) {
#			foreach my $junction_id (keys %{$h_proj_junctions->{$chr_id}}) {
#				my $jid = $junction_id;
#				$jid =~ s/_junction//;
#				if (exists $h_all_junctions_pheno->{$chr_id}->{$junction_id}) {
#					$h_all_junctions_pheno->{$chr_id}->{$jid}->{dejavu}->{$this_project_name} = $h_proj_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$this_project_name};
#				}
#				else {
#					$h_all_junctions_pheno->{$chr_id}->{$jid} = $h_proj_junctions->{$chr_id}->{$junction_id};
#				}
#			}	
#		}
#	} 
#	
#	my $dir_pheno = $dv_dir_path.'/'.$pheno_name.'/';
#	unless (-d $dir_pheno){
#		`mkdir $dir_pheno`;
#		`chmod 777 $dir_pheno`;
#	}
#	my $nodejavupheno = GenBoNoSqlDejaVuJunctionsPhenotype->new( phenotype_name => $pheno_name, dir => $dir_pheno, mode => "c", is_compress => 1 );
#	insert_in_dejavu_jonctions($nodejavupheno, $h_all_junctions_pheno);
#	$nodejavupheno->close();
#	print "# INSERT in DejaVuLMDB Junctions $pheno_name -> DONE!\n";
#}


#sub insert_in_dejavu_jonctions {
#	my ($nodejavu, $h_junctions_to_insert) = @_;
##	print "\n-> DIR: ".$nodejavu->dir()."\n";
#	foreach my $chr_id (keys %{$h_junctions_to_insert}) {
#		$nodejavu->create_table($chr_id);
#		my $sth = $nodejavu->dbh($chr_id)->prepare('insert into  __DATA__(_key,_value,start,end,variation_type,patients,projects,ratios,gene_name)  values(?,?,?,?,?,?,?,?,?) ;') or die $DBI::errstr;
#		$sth->execute();
#		my $tree;
#		foreach my $junction_id (keys %{$h_junctions_to_insert->{$chr_id}}) {
#			my $type = 'all';
#			my $start = $h_junctions_to_insert->{$chr_id}->{$junction_id}->{start};
#			my $end = $h_junctions_to_insert->{$chr_id}->{$junction_id}->{end};
#			
#			my $value = $nodejavu->encode($h_junctions_to_insert->{$chr_id}->{$junction_id}->{dejavu});
#			my (@l_proj, @l_pat, @l_ratio, @l_gene_name);
#			foreach my $proj (sort keys %{$h_junctions_to_insert->{$chr_id}->{$junction_id}->{dejavu}}) {
#				#my $patients = join(',', keys %{$h_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$proj}});
#				
#				my (@local_pat, @local_ratios, @local_gene);
#				foreach my $pat_name (sort keys %{$h_junctions_to_insert->{$chr_id}->{$junction_id}->{dejavu}->{$proj}}) {
#					push(@local_pat, $pat_name);
#					push(@local_ratios, $h_junctions_to_insert->{$chr_id}->{$junction_id}->{dejavu}->{$proj}->{$pat_name}->{score});
#					push(@local_gene, $h_junctions_to_insert->{$chr_id}->{$junction_id}->{dejavu}->{$proj}->{$pat_name}->{gene_name});
#				}
#				my $patients = join(',',@local_pat);
#				my $ratios = join(',',@local_ratios);
#				my $genes_names = join(',',@local_gene);
#				$proj =~ s/NGS20//;
#				push(@l_proj, $proj);
#				push(@l_pat, $patients);
#				push(@l_ratio, $ratios);
#				push(@l_gene_name, $genes_names);
#			}
#			my $pr = join(';', @l_proj);
#			my $pt = join(';', @l_pat);
#			my $ra = join(';', @l_ratio);
#			my $gn = join(';', @l_gene_name);
#			$sth->execute($junction_id, $value, $start, $end, $type, $pt, $pr, $ra, $gn);
#			push(@$tree, [$junction_id, $start, $end]);
#		}
#		$nodejavu->dbh($chr_id)->do(qq{CREATE UNIQUE INDEX if not exists _key_idx  on __DATA__ (_key);});
#		$nodejavu->dbh($chr_id)->do(qq{CREATE  INDEX if not exists _start_idx  on __DATA__ (start);});
#		$nodejavu->dbh($chr_id)->do(qq{CREATE  INDEX if not exists _end_idx  on __DATA__ (end);});
#		$nodejavu->dbh($chr_id)->do(qq{CREATE  INDEX if not exists _type_idx  on __DATA__ (variation_type);});
#		$nodejavu->dbh($chr_id)->do(qq{CREATE  INDEX if not exists _value_idx  on __DATA__ (_value);});
#		$nodejavu->dbh($chr_id)->do(qq{CREATE  INDEX if not exists _gene_name_idx  on __DATA__ (gene_name);});
#	}
#}

sub insert_in_dejavu_jonctions_resume {
	my ($nodejavu, $h_junctions_to_insert) = @_;
#	print "\n-> DIR: ".$nodejavu->dir()."\n";
	foreach my $chr_id (keys %{$h_junctions_to_insert}) {
		$nodejavu->create_table($chr_id);
		my $sth = $nodejavu->dbh($chr_id)->prepare('insert into  __DATA__(_key,start,end,patients,patients_ratios10,patients_ratios20,details,details_ratios10,details_ratios20)  values(?,?,?,?,?,?,?,?,?) ;') or die $DBI::errstr;
		$sth->execute();
		my $tree;
		foreach my $junction_id (keys %{$h_junctions_to_insert->{$chr_id}}) {
			my $start = $h_junctions_to_insert->{$chr_id}->{$junction_id}->{start};
			my $end = $h_junctions_to_insert->{$chr_id}->{$junction_id}->{end};
			my $nb_pat = $h_junctions_to_insert->{$chr_id}->{$junction_id}->{pat};
			my $nb_pat_r10 = $h_junctions_to_insert->{$chr_id}->{$junction_id}->{pat_r10};
			my $nb_pat_r20 = $h_junctions_to_insert->{$chr_id}->{$junction_id}->{pat_r20};
			my $details = $h_junctions_to_insert->{$chr_id}->{$junction_id}->{details};
			my $details_ratios10 = $h_junctions_to_insert->{$chr_id}->{$junction_id}->{details_r10};
			my $details_ratios20 = $h_junctions_to_insert->{$chr_id}->{$junction_id}->{details_r20};
			
			#my $value = $nodejavu->encode($h_junctions_to_insert->{$chr_id}->{$junction_id});
			my (@l_proj, @l_pat, @l_ratio, @l_gene_name);
			$sth->execute($junction_id, $start, $end, $nb_pat, $nb_pat_r10, $nb_pat_r20, $details,$details_ratios10, $details_ratios20);
			push(@$tree, [$junction_id, $start, $end]);
		}
		$nodejavu->dbh($chr_id)->do(qq{CREATE UNIQUE INDEX if not exists _key_idx  on __DATA__ (_key);});
		$nodejavu->dbh($chr_id)->do(qq{CREATE  INDEX if not exists _start_idx  on __DATA__ (start);});
		$nodejavu->dbh($chr_id)->do(qq{CREATE  INDEX if not exists _end_idx  on __DATA__ (end);});
	}
}

 	