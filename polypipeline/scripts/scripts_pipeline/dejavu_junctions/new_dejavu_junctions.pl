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

print "\n# Checking RNA-Projects Junctions";

my $toto = 0;
my ($hash_projects, $hash_project_phenotypes);
foreach my $h (@{$buffer->getQuery->getListProjectsRnaSeq()}) {
	next if not $h->{relname} =~ /$release/;
	my $h_users;
	foreach my $this_user_name (split(',', $h->{username})) { $h_users->{$this_user_name} = undef; }
	my $name = $h->{name};
	
	my $b1 = GBuffer->new;
	my $p1 = $b1->newProjectCache( -name => $name );
	$h->{button} = '';
	if (-d $p1->get_path_rna_seq_junctions_analyse_all_res()) {
		my $ok;
		my $path = $p1->get_path_rna_seq_junctions_analyse_all_res();
		my $se_file = $path.'/allResSE.txt';
		my $ri_file = $path.'/allResRI.txt';
		$ok = 1 if (-e $se_file);
		$ok = 1 if (-e $ri_file);
		if ($ok) {
			$hash_projects->{$name} = $h->{relname};
			if ($release =~ /HG/) {
				foreach my $pheno_name (@{$p1->phenotypes()}) {
					$pheno_name =~ s/ /_/g;
					$hash_project_phenotypes->{lc($pheno_name)}->{$name} = undef;
				}
			}
			$toto++;
		}
	}
	else {
		my $ok_dragen;
		foreach my $patient (@{$p1->getPatients()}) { 
			my $dragen_file = $p1->getVariationsDir('dragen-sj').'/'.$patient->name().'.SJ.out.tab.gz';
			$ok_dragen = 1 if (-e $dragen_file);
		}
		if ($ok_dragen) {
			$hash_projects->{$name} = $h->{relname};
			if ($release =~ /HG/) {
				foreach my $pheno_name (@{$p1->phenotypes()}) {
					$pheno_name =~ s/ /_/g;
					$hash_project_phenotypes->{lc($pheno_name)}->{$name} = undef;
				}
			}
			$toto++;
		}
	}
}
print " -> Done!\n";
print " -> Found ".scalar(keys %{$hash_projects})." projects\n\n";



foreach my $this_project_name (keys %$hash_projects) {
	my $dir_dv_proj = $dv_dir_path.'/projects/'.$this_project_name.'/';
	unless (-d $dir_dv_proj) {
		warn "WARN DejaVu $this_project_name not done yet....";
	}
	unless (-e $dir_dv_proj.'/'.$this_project_name.'.json') {
		warn "WARN DejaVu $this_project_name not done yet....";
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
foreach my $pheno_name (keys %{$hash_project_phenotypes}) {
	push(@lTypes, $pheno_name);
}

my $pm2 = new Parallel::ForkManager($fork);
foreach my $type (@lTypes) {
	my $pid = $pm2->start and next;
	if ($type eq 'all') { insert_all_junctions(); }
	elsif ($type eq 'canoniques') { insert_all_junctions_canoniques(); }
	else { insert_all_junctions_type($type); }
	$pm2->finish();
}
$pm2->wait_all_children();


print "\n\nALL DONE!\n\n";


sub insert_all_junctions {
	my $h_all_junctions;
	foreach my $this_project_name (keys %$hash_projects) {
		my $dir_dv_proj = $dv_dir_path.'/projects/'.$this_project_name.'/';
		my $file = $dir_dv_proj.'/'.$this_project_name.'.json';
		die("\n\nNot exists: $file  - DIE\n") if (not -e $file);
		open (FILE, $file);
		my $json = <FILE>;
		my $h_proj_junctions = decode_json $json;
		close (FILE);
		foreach my $chr_id (keys %{$h_proj_junctions}) {
			foreach my $junction_id (keys %{$h_proj_junctions->{$chr_id}}) {
				if (exists $h_all_junctions->{$chr_id}->{$junction_id}) {
					$h_all_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$this_project_name} = $h_proj_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$this_project_name};
				}
				else {
					$h_all_junctions->{$chr_id}->{$junction_id} = $h_proj_junctions->{$chr_id}->{$junction_id};
				}
			}	
		}
	}
	my $nodejavu = GenBoNoSqlDejaVuJunctions->new( dir => $dir, mode => "c" );
	insert_in_dejavu_jonctions($nodejavu, $h_all_junctions);
	$nodejavu->close();
	print "# INSERT in DejaVuLMDB Junctions -> DONE!\n";
}

sub insert_all_junctions_canoniques {
	my $h_all_junctions_canoniques;
	foreach my $this_project_name (keys %$hash_projects) {
		my $dir_dv_proj = $dv_dir_path.'/projects/'.$this_project_name.'/';
		my $file = $dir_dv_proj.'/'.$this_project_name.'.canoniques.json';
		die("\n\nNot exists: $file  - DIE\n") if (not -e $file);
		open (FILE, $file);
		my $json = <FILE>;
		my $h_proj_junctions_canoniques; 
		$h_proj_junctions_canoniques = decode_json $json if ($json);
		close (FILE);
		next if (not $h_proj_junctions_canoniques); 
		foreach my $chr_id (keys %{$h_proj_junctions_canoniques}) {
			foreach my $junction_id (keys %{$h_proj_junctions_canoniques->{$chr_id}}) {
				if (exists $h_all_junctions_canoniques->{$chr_id}->{$junction_id}) {
					$h_all_junctions_canoniques->{$chr_id}->{$junction_id}->{dejavu}->{$this_project_name} = $h_proj_junctions_canoniques->{$chr_id}->{$junction_id}->{dejavu}->{$this_project_name};
				}
				else {
					$h_all_junctions_canoniques->{$chr_id}->{$junction_id} = $h_proj_junctions_canoniques->{$chr_id}->{$junction_id};
				}
			}	
		}
	}
	my $nodejavucanoniques = GenBoNoSqlDejaVuJunctionsCanoniques->new( dir => $dir_canoniques, mode => "c" );
	insert_in_dejavu_jonctions($nodejavucanoniques, $h_all_junctions_canoniques);
	$nodejavucanoniques->close();
	print "# INSERT in DejaVuLMDB Junctions canoniques -> DONE!\n";
}

sub insert_all_junctions_type {
	my ($pheno_name) = @_;
	my $h_all_junctions_pheno;
	foreach my $this_project_name (keys %{$hash_project_phenotypes->{$pheno_name}}) {
		my $dir_dv_proj = $dv_dir_path.'/projects/'.$this_project_name.'/';
		my $file = $dir_dv_proj.'/'.$this_project_name.'.json';
		die("\n\nNot exists: $file  - DIE\n") if (not -e $file);
		open (FILE, $file);
		my $json = <FILE>;
		my $h_proj_junctions = decode_json $json;
		close (FILE);
		foreach my $chr_id (keys %{$h_proj_junctions}) {
			foreach my $junction_id (keys %{$h_proj_junctions->{$chr_id}}) {
				if (exists $h_all_junctions_pheno->{$chr_id}->{$junction_id}) {
					$h_all_junctions_pheno->{$chr_id}->{$junction_id}->{dejavu}->{$this_project_name} = $h_proj_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$this_project_name};
				}
				else {
					$h_all_junctions_pheno->{$chr_id}->{$junction_id} = $h_proj_junctions->{$chr_id}->{$junction_id};
				}
			}	
		}
	} 
	
	my $dir_pheno = $dv_dir_path.'/'.$pheno_name.'/';
	unless (-d $dir_pheno){
		`mkdir $dir_pheno`;
		`chmod 777 $dir_pheno`;
	}
	my $nodejavupheno = GenBoNoSqlDejaVuJunctionsPhenotype->new( phenotype_name => $pheno_name, dir => $dir_pheno, mode => "c" );
	insert_in_dejavu_jonctions($nodejavupheno, $h_all_junctions_pheno);
	$nodejavupheno->close();
	print "# INSERT in DejaVuLMDB Junctions $pheno_name -> DONE!\n";
}

sub insert_in_dejavu_jonctions {
	my ($nodejavu, $h_junctions_to_insert) = @_;
#	print "\n-> DIR: ".$nodejavu->dir()."\n";
	foreach my $chr_id (keys %{$h_junctions_to_insert}) {
		$nodejavu->create_table($chr_id);
		my $sth = $nodejavu->dbh($chr_id)->prepare('insert into  __DATA__(_key,_value,start,end,variation_type,patients,projects,ratios,gene_name)  values(?,?,?,?,?,?,?,?,?) ;') or die $DBI::errstr;
		$sth->execute();
		my $tree;
		foreach my $junction_id (keys %{$h_junctions_to_insert->{$chr_id}}) {
			my $type = 'all';
			my $start = $h_junctions_to_insert->{$chr_id}->{$junction_id}->{start};
			my $end = $h_junctions_to_insert->{$chr_id}->{$junction_id}->{end};
			
			my $value = $nodejavu->encode($h_junctions_to_insert->{$chr_id}->{$junction_id}->{dejavu});
			my (@l_proj, @l_pat, @l_ratio, @l_gene_name);
			foreach my $proj (sort keys %{$h_junctions_to_insert->{$chr_id}->{$junction_id}->{dejavu}}) {
				#my $patients = join(',', keys %{$h_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$proj}});
				
				my (@local_pat, @local_ratios, @local_gene);
				foreach my $pat_name (sort keys %{$h_junctions_to_insert->{$chr_id}->{$junction_id}->{dejavu}->{$proj}}) {
					push(@local_pat, $pat_name);
					push(@local_ratios, $h_junctions_to_insert->{$chr_id}->{$junction_id}->{dejavu}->{$proj}->{$pat_name}->{score});
					push(@local_gene, $h_junctions_to_insert->{$chr_id}->{$junction_id}->{dejavu}->{$proj}->{$pat_name}->{gene_name});
				}
				my $patients = join(',',@local_pat);
				my $ratios = join(',',@local_ratios);
				my $genes_names = join(',',@local_gene);
				$proj =~ s/NGS20//;
				push(@l_proj, $proj);
				push(@l_pat, $patients);
				push(@l_ratio, $ratios);
				push(@l_gene_name, $genes_names);
			}
			my $pr = join(';', @l_proj);
			my $pt = join(';', @l_pat);
			my $ra = join(';', @l_ratio);
			my $gn = join(';', @l_gene_name);
			$sth->execute($junction_id, $value, $start, $end, $type, $pt, $pr, $ra, $gn);
			push(@$tree, [$junction_id, $start, $end]);
		}
		$nodejavu->dbh($chr_id)->do(qq{CREATE UNIQUE INDEX if not exists _key_idx  on __DATA__ (_key);});
		$nodejavu->dbh($chr_id)->do(qq{CREATE  INDEX if not exists _start_idx  on __DATA__ (start);});
		$nodejavu->dbh($chr_id)->do(qq{CREATE  INDEX if not exists _end_idx  on __DATA__ (end);});
		$nodejavu->dbh($chr_id)->do(qq{CREATE  INDEX if not exists _type_idx  on __DATA__ (variation_type);});
		$nodejavu->dbh($chr_id)->do(qq{CREATE  INDEX if not exists _value_idx  on __DATA__ (_value);});
		$nodejavu->dbh($chr_id)->do(qq{CREATE  INDEX if not exists _gene_name_idx  on __DATA__ (gene_name);});
	}
}

 	