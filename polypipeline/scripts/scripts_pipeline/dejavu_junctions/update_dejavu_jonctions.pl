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

my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => 'NGS2022_5612' );
$project->getChromosomes();


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
print "\n# Checking RNA-Projects Junctions\n";

my ($hash_projects, $hash_project_phenotypes, $hash_projects_already_done);
if (-d $dv_dir_projects_path) {
	opendir DIR, $dv_dir_projects_path;
	my @dir = readdir(DIR);
	close DIR;
	
	foreach my $project_name (@dir) {
		next if not $project_name =~ /NGS/;
		
		my $json1 = $dv_dir_projects_path.'/'.$project_name.'/'.$project_name.'.json';
		my $json2 = $dv_dir_projects_path.'/'.$project_name.'/'.$project_name.'.canoniques.json';
		$hash_projects_already_done->{$project_name} = undef if (-e $json1.'.gz');
		next if not -e $json1;
		next if not -e $json2;
		$hash_projects->{$project_name} = undef;
		
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

print " -> Done!\n";
print " -> Found ALREADY DONE: ".scalar(keys %{$hash_projects_already_done})." projects\n";
print " -> Found TO DO: ".scalar(keys %{$hash_projects})." projects\n\n";


my $h_not_found;
foreach my $this_project_name (keys %$hash_projects) {
	my $dir_dv_proj = $dv_dir_path.'/projects/'.$this_project_name.'/';
	unless (-d $dir_dv_proj) {
		warn "WARN DejaVu $this_project_name not done yet....";
		$h_not_found->{$this_project_name} = undef;
	}
	unless (-e $dir_dv_proj.'/'.$this_project_name.'.json') {
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
foreach my $pheno_name (keys %{$hash_project_phenotypes}) {
	push(@lTypes, $pheno_name);
}

insert_all_junctions();


zip_files($hash_projects);

#foreach my $type (@lTypes) {
#	if ($type eq 'all') { insert_all_junctions(); }
##	elsif ($type eq 'canoniques') { insert_all_junctions_canoniques(); }
##	else { insert_all_junctions_type($type); }
#}


print "\n\nALL DONE!\n\n";


sub zip_files {
	my ($hash_projects) = @_;
	foreach my $this_project_name (keys %$hash_projects) {
		my $dir_dv_proj = $dv_dir_path.'/projects/'.$this_project_name.'/';
		my $file1 = $dir_dv_proj.'/'.$this_project_name.'.json';
		my $file2 = $dir_dv_proj.'/'.$this_project_name.'.canoniques.json';
		if (-e $file1) {
			my $cmd = "bgzip -f $file1";
			`$cmd`;
		}
		if (-e $file2) {
			my $cmd = "bgzip -f $file2";
			`$cmd`;
		}
	}
	
}


sub insert_all_junctions {
	my ($h_all_junctions, @lCmds);
	foreach my $this_project_name (keys %$hash_projects) {
		my $dir_dv_proj = $dv_dir_path.'/projects/'.$this_project_name.'/';
		my $file = $dir_dv_proj.'/'.$this_project_name.'.json';
		next if exists $h_not_found->{$this_project_name};
		die("\n\nNot exists: $file  - DIE\n") if (not -e $file);
		open (FILE, $file);
		my $json = <FILE>;
		my $h_proj_junctions;
		eval {
			$h_proj_junctions = decode_json $json;	
			close (FILE);
		};
		if ($@) {
			close (FILE);
			warn 'pb with project '.$this_project_name;
			next;
		}
		
		my $hres;
		my $file_tab = $dir_dv_proj.'/'.$this_project_name.'.tab';
		open (FILE, ">$file_tab");
		foreach my $chr_id (1..22, 'X', 'Y', 'MT') {
			next unless exists $h_proj_junctions->{$chr_id};
			foreach my $junction_id (keys %{$h_proj_junctions->{$chr_id}}) {
				my $start = $h_proj_junctions->{$chr_id}->{$junction_id}->{start};
				my $end = $h_proj_junctions->{$chr_id}->{$junction_id}->{end};
				foreach my $pat_name (keys %{$h_proj_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$this_project_name}}) {
					$hres->{$chr_id}->{$start}->{$end}->{$pat_name} = $h_proj_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$this_project_name}->{$pat_name};
					$hres->{$chr_id}->{$start}->{$end}->{$pat_name}->{junction_id} = $junction_id;
				}
			}
			
			foreach my $start (sort {$a <=> $b} keys %{$hres->{$chr_id}}) {
				foreach my $end (sort {$a <=> $b} keys %{$hres->{$chr_id}->{$start}}) {
					foreach my $pat_name (sort keys %{$hres->{$chr_id}->{$start}->{$end}}) {
						my $junction_id = $hres->{$chr_id}->{$start}->{$end}->{$pat_name}->{junction_id};
						my $type = $hres->{$chr_id}->{$start}->{$end}->{$pat_name}->{type};
						my $count_junctions = $hres->{$chr_id}->{$start}->{$end}->{$pat_name}->{count_junctions};
						my $count_normal = $hres->{$chr_id}->{$start}->{$end}->{$pat_name}->{count_normal};
						my $ratio = $hres->{$chr_id}->{$start}->{$end}->{$pat_name}->{score};
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
		my $cmd1 = "bgzip -f $file_tab";
		my $cmd2 = "tabix -p bed $file_tab.gz";
		push(@lCmds, $cmd1);
		push(@lCmds, $cmd2);
		
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
#	update_in_dejavu_jonctions($h_all_junctions);
	
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
					my $patname_text = $proj.'_'.$pat_name;
					$patname_text =~ s/NGS20//;
					$h_resume->{pat}++;
					my $new = $h_all_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$proj}->{$pat_name}->{count_junctions};
					my $normal = $h_all_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$proj}->{$pat_name}->{count_normal};
					my $dp = ($new + $normal);
					my $ratio =  ($new / $dp) * 100;
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
	update_in_dejavu_jonctions_resume($h_all_junctions_resume);
	
	foreach my $cmd (@lCmds) { `$cmd`; }
	
	print "\n\nUPDATE DEJAVU JUNCTIONS DONE\n\n";
	
}

sub update_in_dejavu_jonctions_resume {
	my ($h_all_junctions_resume) = @_;
	my $pm = new Parallel::ForkManager($fork);
	foreach my $chr_id (keys %{$h_all_junctions_resume}) {
		#next if $chr_id ne '4';
		my $pid = $pm->start and next;
		my $nodejavu = GenBoNoSqlDejaVuJunctionsResume->new( dir => $dir, name => $chr_id, mode => "w" );
		foreach my $junction_id (keys %{$h_all_junctions_resume->{$chr_id}}) {
			my $start = $h_all_junctions_resume->{$chr_id}->{$junction_id}->{start};
			my $end = $h_all_junctions_resume->{$chr_id}->{$junction_id}->{end};
			my $hres = $nodejavu->get_position($chr_id, $start, $end);
			
			my ($h_tmp, $jid);
			if ($hres) {
				my @lPat_this_all = split(';', $h_all_junctions_resume->{$chr_id}->{$junction_id}->{details});
				my @lPat_this_r10 = split(';', $h_all_junctions_resume->{$chr_id}->{$junction_id}->{details_r10});
				my @lPat_this_r20 = split(';', $h_all_junctions_resume->{$chr_id}->{$junction_id}->{details_r20});

				my @lkeys = keys %$hres;
				$jid = $lkeys[0];
				
				my @lPat_all = split(';', $hres->{$jid}->{details});
				my @lPat_r10 = split(';', $hres->{$jid}->{details_r10});
				my @lPat_r20 = split(';', $hres->{$jid}->{details_r20});
				
				foreach my $pat_name (@lPat_all) { $h_tmp->{details}->{$pat_name} = undef;	}
				foreach my $pat_name (@lPat_this_all) { $h_tmp->{details}->{$pat_name} = undef;	} 
				foreach my $pat_name (@lPat_r10) { $h_tmp->{details_r10}->{$pat_name} = undef;	}
				foreach my $pat_name (@lPat_this_r10) { $h_tmp->{details_r10}->{$pat_name} = undef;	} 
				foreach my $pat_name (@lPat_r20) { $h_tmp->{details_r20}->{$pat_name} = undef;	}
				foreach my $pat_name (@lPat_this_r20) { $h_tmp->{details_r20}->{$pat_name} = undef;	} 
				
				my $nb_pat = scalar (keys %{$h_tmp->{details}});
				my $nb_pat_r10 = scalar (keys %{$h_tmp->{details_r10}});
				my $nb_pat_r20 = scalar (keys %{$h_tmp->{details_r20}});
				my $details = join(';', keys %{$h_tmp->{details}});
				my $details_ratios10 = join(';', keys %{$h_tmp->{details_r10}});
				my $details_ratios20 = join(';', keys %{$h_tmp->{details_r20}});
				
				my $sth_del = $nodejavu->dbh($chr_id)->prepare('DELETE FROM __DATA__  WHERE _key=?');
				$sth_del->execute($jid);	

				my $sth = $nodejavu->dbh($chr_id)->prepare('insert into  __DATA__(_key,start,end,patients,patients_ratios10,patients_ratios20,details,details_ratios10,details_ratios20)  values(?,?,?,?,?,?,?,?,?) ;') or die $DBI::errstr;
				$sth->execute($jid, $start, $end, $nb_pat, $nb_pat_r10, $nb_pat_r20, $details,$details_ratios10, $details_ratios20);
			}
			else {
				$jid = $start.'-'.$end;
				my @lPat_this_all = split(';', $h_all_junctions_resume->{$chr_id}->{$junction_id}->{details});
				my @lPat_this_r10 = split(';', $h_all_junctions_resume->{$chr_id}->{$junction_id}->{details_r10});
				my @lPat_this_r20 = split(';', $h_all_junctions_resume->{$chr_id}->{$junction_id}->{details_r20});
				foreach my $pat_name (@lPat_this_all) { $h_tmp->{details}->{$pat_name} = undef;	} 
				foreach my $pat_name (@lPat_this_r10) { $h_tmp->{details_r10}->{$pat_name} = undef;	} 
				foreach my $pat_name (@lPat_this_r20) { $h_tmp->{details_r20}->{$pat_name} = undef;	} 
			
				my $nb_pat = scalar (keys %{$h_tmp->{details}});
				my $nb_pat_r10 = scalar (keys %{$h_tmp->{details_r10}});
				my $nb_pat_r20 = scalar (keys %{$h_tmp->{details_r20}});
				my $details = join(';', keys %{$h_tmp->{details}});
				my $details_ratios10 = join(';', keys %{$h_tmp->{details_r10}});
				my $details_ratios20 = join(';', keys %{$h_tmp->{details_r20}});
				
				my $sth = $nodejavu->dbh($chr_id)->prepare('insert into  __DATA__(_key,start,end,patients,patients_ratios10,patients_ratios20,details,details_ratios10,details_ratios20)  values(?,?,?,?,?,?,?,?,?) ;') or die $DBI::errstr;
				$sth->execute($jid, $start, $end, $nb_pat, $nb_pat_r10, $nb_pat_r20, $details,$details_ratios10, $details_ratios20);
			}
			
		}
		$nodejavu->dbh($chr_id)->do(qq{CREATE UNIQUE INDEX if not exists _key_idx  on __DATA__ (_key);});
		$nodejavu->dbh($chr_id)->do(qq{CREATE  INDEX if not exists _start_idx  on __DATA__ (start);});
		$nodejavu->dbh($chr_id)->do(qq{CREATE  INDEX if not exists _end_idx  on __DATA__ (end);});
		$nodejavu->close();
		
		print "-> chr$chr_id Done\n";
		$pm->finish();
	}
	$pm->wait_all_children();
	print "# INSERT in DejaVuLMDB Junctions Resume -> DONE!\n";
}


#sub update_in_dejavu_jonctions {
#	my ($h_junctions_to_insert) = @_;
#	
#	my $pm = new Parallel::ForkManager($fork);
#	foreach my $chr_id (keys %{$h_junctions_to_insert}) {
#		my $pid = $pm->start and next;
#		my $nodejavu = GenBoNoSqlDejaVuJunctions->new( dir => $dir, name => $chr_id, mode => "w" );
##		my $tree;
#		foreach my $junction_id (keys %{$h_junctions_to_insert->{$chr_id}}) {
#			my $type = 'all';
#			my $start = $h_junctions_to_insert->{$chr_id}->{$junction_id}->{start};
#			my $end = $h_junctions_to_insert->{$chr_id}->{$junction_id}->{end};
#			
#			my $hres = $nodejavu->get_position($chr_id, $start, $end, $type);
#			
#			if ($hres) {
#				my $h_geneid;
#				foreach my $id (keys %{$hres}) {
#					foreach my $project_name (keys %{$h_junctions_to_insert->{$chr_id}->{$junction_id}->{dejavu}}) {
#						$hres->{$id}->{details}->{$project_name} = $h_junctions_to_insert->{$chr_id}->{$junction_id}->{dejavu}->{$project_name};
#						foreach my $pat_name (keys %{$h_junctions_to_insert->{$chr_id}->{$junction_id}->{dejavu}->{$project_name}}) {
#							my $gene_name = $h_junctions_to_insert->{$chr_id}->{$junction_id}->{dejavu}->{$project_name}->{$pat_name}->{gene_name};
#							$h_geneid->{$gene_name} = undef;
#						}
#					}
#				}
#				foreach my $id (keys %{$hres}) {
#					my $gn = $hres->{$id}->{gene_id};
#					$gn =~ s/,/;/g;
#					foreach my $gene_name (split(';', $gn)) {
#						$h_geneid->{$gene_name} = undef;
#					}
#					$hres->{$id}->{gene_id} = join(',', sort keys %{$h_geneid});
#					my $value = $nodejavu->encode($hres->{$id}->{details});
#					my (@l_proj, @l_pat, @l_ratio);
#					foreach my $proj (sort keys %{$hres->{$id}->{details}}) {
#						my (@local_pat, @local_ratios);
#						foreach my $pat_name (sort keys %{$hres->{$id}->{details}->{$proj}}) {
#							push(@local_pat, $pat_name);
#							push(@local_ratios, $hres->{$id}->{details}->{$proj}->{$pat_name}->{score});
#						}
#						my $patients = join(',',@local_pat);
#						my $ratios = join(',',@local_ratios);
#						$proj =~ s/NGS20//;
#						push(@l_proj, $proj);
#						push(@l_pat, $patients);
#						push(@l_ratio, $ratios);
#					}
#					my $pr = join(';', @l_proj);
#					my $pt = join(';', @l_pat);
#					my $ra = join(';', @l_ratio);
#					my $gn = join(';', sort keys %{$h_geneid});
#					my $sth = $nodejavu->dbh($chr_id)->prepare('UPDATE __DATA__ SET _value=?, patients=?, projects=?, ratios=?, gene_name=? WHERE start=? and end=? and variation_type=?');
#					$sth->execute($value, $pt, $pr, $ra, $gn, $start, $end, $type);
#				}
##				push(@$tree, [$junction_id, $start, $end]);
#			}
#			else {
#				my $value = $nodejavu->encode($h_junctions_to_insert->{$chr_id}->{$junction_id}->{dejavu});
#				my (@l_proj, @l_pat, @l_ratio, @l_gene_name);
#				foreach my $proj (sort keys %{$h_junctions_to_insert->{$chr_id}->{$junction_id}->{dejavu}}) {
#					#my $patients = join(',', keys %{$h_junctions->{$chr_id}->{$junction_id}->{dejavu}->{$proj}});
#					
#					my (@local_pat, @local_ratios, @local_gene);
#					foreach my $pat_name (sort keys %{$h_junctions_to_insert->{$chr_id}->{$junction_id}->{dejavu}->{$proj}}) {
#						push(@local_pat, $pat_name);
#						push(@local_ratios, $h_junctions_to_insert->{$chr_id}->{$junction_id}->{dejavu}->{$proj}->{$pat_name}->{score});
#						push(@local_gene, $h_junctions_to_insert->{$chr_id}->{$junction_id}->{dejavu}->{$proj}->{$pat_name}->{gene_name});
#					}
#					my $patients = join(',',@local_pat);
#					my $ratios = join(',',@local_ratios);
#					my $genes_names = join(',',@local_gene);
#					$proj =~ s/NGS20//;
#					push(@l_proj, $proj);
#					push(@l_pat, $patients);
#					push(@l_ratio, $ratios);
#					push(@l_gene_name, $genes_names);
#				}
#				my $pr = join(';', @l_proj);
#				my $pt = join(';', @l_pat);
#				my $ra = join(';', @l_ratio);
#				my $gn = join(';', @l_gene_name);
#				my $sth = $nodejavu->dbh($chr_id)->prepare('insert into  __DATA__(_key,_value,start,end,variation_type,patients,projects,ratios,gene_name)  values(?,?,?,?,?,?,?,?,?) ;') or die $DBI::errstr;
#				$sth->execute($junction_id, $value, $start, $end, $type, $pt, $pr, $ra, $gn);
##				push(@$tree, [$junction_id, $start, $end]);
#			}
#		}
#		
#		$nodejavu->dbh($chr_id)->do(qq{CREATE UNIQUE INDEX if not exists _key_idx  on __DATA__ (_key);});
#		$nodejavu->dbh($chr_id)->do(qq{CREATE  INDEX if not exists _start_idx  on __DATA__ (start);});
#		$nodejavu->dbh($chr_id)->do(qq{CREATE  INDEX if not exists _end_idx  on __DATA__ (end);});
#		$nodejavu->dbh($chr_id)->do(qq{CREATE  INDEX if not exists _type_idx  on __DATA__ (variation_type);});
#		$nodejavu->dbh($chr_id)->do(qq{CREATE  INDEX if not exists _value_idx  on __DATA__ (_value);});
#		$nodejavu->dbh($chr_id)->do(qq{CREATE  INDEX if not exists _gene_name_idx  on __DATA__ (gene_name);});
#		$pm->finish();
#	}
#	$pm->wait_all_children();
#	
#}

