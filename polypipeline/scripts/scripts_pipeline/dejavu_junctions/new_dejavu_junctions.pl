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

my $release = $cgi->param('release');
if (not $release) { $release = 'HG19'; }
elsif ($release =~ /HG19/) { $release = 'HG19'; }


my $dir = $project->DejaVuJunction_path();

unless (-d $dir){
	`mkdir $dir`;
	`chmod 777 $dir`;
}


print "\n# Checking RNA-Projects Junctions";

#my $nb = 0;


foreach my $h (@{$buffer->getQuery->getListProjectsRnaSeq()}) {
#	warn $h->{name} if $h->{relname} eq 'MM38';
#	next if not $h->{relname} eq 'MM38';
	
	next if not $h->{relname} =~ /$release/;
	
	my $b1 = GBuffer->new;
	my $p1 = $b1->newProjectCache( -name => $h->{name} );
	if (-d $p1->get_path_rna_seq_junctions_analyse_all_res()) {
		warn $h->{name}.' - '.$h->{relname}.' ---> OK !!';
		my $ok;
		my $path = $p1->get_path_rna_seq_junctions_analyse_all_res();
		my $se_file = $path.'/allResSE.txt';
		my $ri_file = $path.'/allResRI.txt';
		$ok = 1 if (-e $se_file);
		$ok = 1 if (-e $ri_file);
		warn $h->{name}.' - '.$h->{relname}.' ------> OK 2 !!!!!' if $ok;
	}
}
#die;

my $hash_projects;
foreach my $h (@{$buffer->getQuery->getListProjectsRnaSeq()}) {
	my $h_users;
	foreach my $this_user_name (split(',', $h->{username})) { $h_users->{$this_user_name} = undef; }
	my $name = $h->{name};
	my $b1 = GBuffer->new;
	my $p1 = $b1->newProjectCache( -name => $name );
	$h->{button} = '';
	if (-d $p1->get_path_rna_seq_junctions_analyse_all_res()) {
#		$nb++;
#		last if $nb == 4;
#		warn $name;
		my $ok;
		my $path = $p1->get_path_rna_seq_junctions_analyse_all_res();
		my $se_file = $path.'/allResSE.txt';
		my $ri_file = $path.'/allResRI.txt';
		$ok = 1 if (-e $se_file);
		$ok = 1 if (-e $ri_file);
		$hash_projects->{$name} = $h->{relname} if $ok;
	}
}
print " -> Done!\n";
print "   -> Found ".scalar(keys %{$hash_projects})." projects\n\n";

#warn Dumper $hash_projects;

#die;

my ($h_junctions, $h_junctions_canonique);
foreach my $this_project_name (keys %$hash_projects) {
	print "$this_project_name "; 
	my $buffer_tmp = GBuffer->new();
	my $project_tmp = $buffer_tmp->newProject( -name => $this_project_name );
	if (not $project_tmp->annotation_genome_version() =~ /$release/) {
		print "-> SKIPPED (".$project_tmp->annotation_genome_version().")\n";
		$project_tmp = undef;
		$buffer_tmp = undef;
		next;
	}
	print " -> OK release ".$project_tmp->annotation_genome_version();
	$project_tmp->getChromosomes();
	
	my $hType_patients;
	$hType_patients = $project_tmp->get_hash_patients_description_rna_seq_junction_analyse() if (-d $project_tmp->get_path_rna_seq_junctions_analyse_description_root());
	
	foreach my $this_patient (@{$project_tmp->getPatients()}) {
		#$this_patient->use_not_filtred_junction_files(0);
		if (($hType_patients and exists $hType_patients->{$this_patient->name()}->{pat}) or not $hType_patients) {
			my @lJunctions;
			eval { @lJunctions = @{$this_patient->getJunctions()}; };
			if ($@) {
				print "-> ERROR ".$this_patient->name();
				next;
			}
			foreach my $junction (@lJunctions) {
				#next if ($junction->isCanonique($this_patient));
				next if ($junction->get_ratio_new_count($this_patient) == 1);
				$junction->getPatients();
				my $type = $junction->getTypeDescription($this_patient);
				my $chr_id = $junction->getChromosome->id();
				my $start = $junction->start();
				my $end = $junction->end();
				my $gene_name = $junction->annex->{$this_patient->name()}->{ensid};
				my $count_new_junction = $junction->get_nb_new_count($this_patient);
				my $count_normal_junction = $junction->get_nb_normal_count($this_patient);
				my $score = int($junction->get_percent_new_count($this_patient));
				my $junction_id = $chr_id.'_'.$start.'_'.$end.'_junction';
				
				if ($junction->isCanonique($this_patient)) {
					if (not exists $h_junctions_canonique->{$chr_id}->{$junction_id}) {
						$h_junctions_canonique->{$chr_id}->{$junction_id}->{start} = $start;
						$h_junctions_canonique->{$chr_id}->{$junction_id}->{end} = $end;
					}
					$h_junctions_canonique->{$chr_id}->{$junction_id}->{dejavu}->{$project_tmp->name()}->{$this_patient->name()}->{count_junctions} = $count_new_junction;
					$h_junctions_canonique->{$chr_id}->{$junction_id}->{dejavu}->{$project_tmp->name()}->{$this_patient->name()}->{count_normal} = $count_normal_junction;
					$h_junctions_canonique->{$chr_id}->{$junction_id}->{dejavu}->{$project_tmp->name()}->{$this_patient->name()}->{score} = $score;
					$h_junctions_canonique->{$chr_id}->{$junction_id}->{dejavu}->{$project_tmp->name()}->{$this_patient->name()}->{type} = $type;
					$h_junctions_canonique->{$chr_id}->{$junction_id}->{dejavu}->{$project_tmp->name()}->{$this_patient->name()}->{gene_name} = $gene_name;
				}
				else {
					if (not exists $h_junctions->{$chr_id}->{$junction_id}) {
						$h_junctions->{$chr_id}->{$junction_id}->{start} = $start;
						$h_junctions->{$chr_id}->{$junction_id}->{end} = $end;
					}
					$h_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$project_tmp->name()}->{$this_patient->name()}->{count_junctions} = $count_new_junction;
					$h_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$project_tmp->name()}->{$this_patient->name()}->{count_normal} = $count_normal_junction;
					$h_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$project_tmp->name()}->{$this_patient->name()}->{score} = $score;
					$h_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$project_tmp->name()}->{$this_patient->name()}->{type} = $type;
					$h_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$project_tmp->name()}->{$this_patient->name()}->{gene_name} = $gene_name;
				}
			}
		}
	}
	print "-> Done!\n";
	$project_tmp = undef;
	$buffer_tmp = undef;
}


print "\n# INSERT in DejaVuLMDB Junctions\n";
my $nodejavu = GenBoNoSqlDejaVuJunctions->new( dir => $dir, mode => "c" );
insert_in_dejavu_jonctions($nodejavu, $h_junctions);
$nodejavu->close();
print "-> DONE!\n\n";

print "\n# INSERT in DejaVuLMDB Junctions Canoniques\n";
my $nodejavucanoniques = GenBoNoSqlDejaVuJunctionsCanoniques->new( dir => $dir, mode => "c" );
insert_in_dejavu_jonctions($nodejavucanoniques, $h_junctions_canonique);
$nodejavucanoniques->close();
print "-> DONE!\n\n";



sub insert_in_dejavu_jonctions {
	my ($nodejavu, $h_junctions_to_insert) = @_;
	print "\n-> DIR: $dir\n";
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

 	