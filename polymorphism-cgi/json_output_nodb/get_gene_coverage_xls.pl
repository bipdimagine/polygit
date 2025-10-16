#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";

use Data::Dumper;
use GBuffer;
use xls_export;
use Carp;
use CGI;

my $fork = 1;
my $cgi = new CGI();
my $dir_out			= $cgi->param('dir_out');
my $project_name	= $cgi->param('project');
my $gene_name		= $cgi->param('gene');
my $only_coding		= $cgi->param('only_coding');
my $intronic		= $cgi->param('intronic');

my $details_by_pos	= $cgi->param('details_by_pos');
my $session_id	= $cgi->param('session_id');

#confess("\n\nOption -dir_out mandatory. Die\n\n") unless $dir_out;
confess("\n\nOption -project mandatory. Die\n\n") unless $project_name;
confess("\n\nOption -gene mandatory. Die\n\n") unless $gene_name;


my $buffer = new GBuffer;
my $project = $buffer->newProject( -name => $project_name);
my $gene = $project->newGene($gene_name);

my @lPatients = @{$project->getPatients()};
my $h_patients;
foreach my $p (@lPatients) {
	$h_patients->{$p->name()} = undef;
}
#print "\n# NB PATIENTS: ".scalar(@lPatients)."\n";



my $h_res;

if ($session_id) {
	my $xls_export = get_xls_handler();
	$h_res = loadSessionXLS($session_id);
	if ($details_by_pos) {
		write_xls_pos($xls_export, $gene, $h_res);
		exit(0);
	}
	else {
		foreach my $t (sort @{$gene->getTranscripts()}) {
			write_xls_transcript($xls_export, $t, $h_res);
			exit(0);
		}
	}
	exit(0);
}
else {
	if ($details_by_pos) {
		print $cgi->header('text/json-comment-filtered');
		print "{\"progress\":\".";
		my $start = $gene->start();
		my $end = $gene->end();
		my $pos = $start;
		while ($pos <= $end) {
			$project->print_dot(250);
			foreach my $p (@{$project->getPatients()}) {
				my $cov = get_coverage_at_positions_start($p, $gene->getChromosome->name, $pos);
				$h_res->{$pos}->{'all'}->{$p->name()} = $cov;
			}
			$pos++;
		}
		saveSessionXLS($project, $h_res);
	}
	
	else {
		print $cgi->header('text/json-comment-filtered');
		print "{\"progress\":\".";
		foreach my $t (sort @{$gene->getTranscripts()}) {
			$project->print_dot(1);
			#warn Dumper $t->getSpanCoding();
			
			my $t_id = $t->id();
			foreach my $e (@{$t->getExons()}) {
				my $e_id = $e->id();
				$e_id =~ s/$t_id//;
				foreach my $p (@lPatients) {
					my $h = get_coverage_at_positions($p, $gene->getChromosome->name, $t, $e);
					my $pos = $e_id;
					$pos =~ s/ex//;
					$h_res->{$pos}->{$e_id}->{$p->name()} = $h;
				}
			}
			if ($intronic) {
				foreach my $i (@{$t->getIntrons()}) {
					my $i_id = $i->id();
					$i_id =~ s/$t_id//;
					foreach my $p (@{$project->getPatients()}) {
						my $h = get_coverage_at_positions($p, $gene->getChromosome->name, $t, $i);
						my $pos = $i_id;
						$pos =~ s/intron//;
						$h_res->{$pos}->{$i_id}->{$p->name()} = $h;
					}
				}
			}
		}
		saveSessionXLS($project, $h_res);
		exit(0);
	}
}






sub saveSessionXLS {
	my ($project, $hashRes) = @_;
	print ".";
	my $xls_export = new xls_export();
	$xls_export->{hash_specific_infos} = $hashRes;
	my $session_id = $xls_export->save();
	print '@@@';
	print "\",\"session_id\":\"";
	print $session_id;
	print "\"}";
    exit(0);
}

sub loadSessionXLS {
	my ($session_id) = @_;
	my $session = new xls_export();
	$session->load($session_id);
	my $hashRes = $session->{hash_specific_infos};
	return $hashRes;
}

sub get_xls_handler {
	my $xls_export = new xls_export();
	my $title_xls = 'Coverage_'.$project->name().'_'.$gene_name;
	$title_xls .= '_only_coding' if ($only_coding);
	$title_xls .= '.xls';
	$xls_export->title_page($title_xls);
	if ($dir_out) { $xls_export->output_dir($dir_out); }
	$xls_export->open_xls_file();
	return $xls_export;
}

sub write_xls_transcript {
	my ($xls_export, $t, $h_res) = @_;
	my $xls_tr = $xls_export->workbook->add_worksheet($t->id());
	my @lHeaders;
	push(@lHeaders, 'Exon/Intron');
	foreach my $p_name (sort keys %{$h_patients}) { push(@lHeaders, $p_name); }
	my $i = 0;
	my $j = 0;
	foreach my $cat (@lHeaders) {
		$xls_tr->write( $i, $j, $cat, $xls_export->get_format_header() );
		$j++;
	}
	$xls_tr->write( $i, $j, 'Mean Project', $xls_export->get_format_header() );
	
	my $h_lines;
	my $nb_line = 1;
	
	my @lPos = sort {$a <=> $b} keys %$h_res;
	foreach my $pos (@lPos) {
		foreach my $e_i_id (sort keys %{$h_res->{$pos}}) {
			my @_this_line;
			push(@_this_line, $e_i_id);
			
			$h_lines->{$nb_line}->{format}->{0} = $xls_export->get_format_header();
			
			my $nb_pat = 0;
			my $total_pat = 0;
			foreach my $p_name (sort keys %{$h_res->{$pos}->{$e_i_id}}) {
				$nb_pat++;
				$total_pat += $h_res->{$pos}->{$e_i_id}->{$p_name}->{mean};
			}
			my $mean_pat = $total_pat / $nb_pat;
			
			$nb_pat = 0;
			foreach my $p_name (sort keys %{$h_res->{$pos}->{$e_i_id}}) {
				$nb_pat++;
				my $min = $h_res->{$pos}->{$e_i_id}->{$p_name}->{min};
				my $max = $h_res->{$pos}->{$e_i_id}->{$p_name}->{max};
				my $mean = $h_res->{$pos}->{$e_i_id}->{$p_name}->{mean};
				my $text = "$mean ($min->$max)";
				push(@_this_line, $text);
				
				my $limit_down = $mean_pat * 0.6;
				my $limit_up = $mean_pat * 1.4;
				if ($min == 0) {
					$h_lines->{$nb_line}->{format}->{$nb_pat} = $xls_export->get_format_red();
				}
				elsif ($mean < $limit_down) {
					$h_lines->{$nb_line}->{format}->{$nb_pat} = $xls_export->get_format_orange();
				}
				elsif ($mean > $limit_up) {
					$h_lines->{$nb_line}->{format}->{$nb_pat} = $xls_export->get_format_blue();
				}
				else {
					$h_lines->{$nb_line}->{format}->{$nb_pat} = undef;
				}
			}
			my $text = sprintf("%.2f", $mean_pat);
			push(@_this_line, $text);
			$h_lines->{$nb_line}->{line} = \@_this_line;
			$nb_line++;
		}
	}
	
	foreach my $nb_line (sort {$a <=> $b} keys %{$h_lines}) {
		$i = $nb_line;
		$j = 0;
		foreach my $col (@{$h_lines->{$nb_line}->{line}}) {
			my $format = $h_lines->{$nb_line}->{format}->{$j};
			$xls_tr->write( $i, $j, $col, $format );
			$j++;
		}
	}
	
	$i += 3;
	$xls_tr->write( $i, 0, 'Legend Value', $xls_export->get_format_header() );
	$xls_tr->write( $i, 1, 'Mean Cov (Min Cov -> Max Cov)');
	
	$i += 2;
	$xls_tr->write( $i, 0, 'RED Value', $xls_export->get_format_red() );
	$xls_tr->write( $i, 1, 'Min coverage 0');
	
	$i++;
	$xls_tr->write( $i, 0, 'ORANGE Value', $xls_export->get_format_orange() );
	$xls_tr->write( $i, 1, 'Mean coverage of this patient < Mean coverage project * 0.6');
	
	$i++;
	$xls_tr->write( $i, 0, 'BLUE Value', $xls_export->get_format_blue() );
	$xls_tr->write( $i, 1, 'Mean coverage of this patient > Mean coverage project * 1.4');
}

sub write_xls_pos {
	my ($xls_export, $g, $h_res) = @_;
	my $xls_tr = $xls_export->workbook->add_worksheet('ALL POS');
	my $xls_cov_low = $xls_export->workbook->add_worksheet('COV < 30');
	my $xls_cov_null = $xls_export->workbook->add_worksheet('COV 0');
	my @lHeaders;
	push(@lHeaders, 'chr'.$g->getChromosome->id().' pos');
	foreach my $p_name (sort keys %{$h_patients}) { push(@lHeaders, $p_name); }
	my $i = 0;
	my $j = 0;
	foreach my $cat (@lHeaders) {
		$xls_tr->write( $i, $j, $cat, $xls_export->get_format_header() );
		$j++;
	}
	$j = 0;
	foreach my $cat (@lHeaders) {
		$xls_cov_low->write( $i, $j, $cat, $xls_export->get_format_header() );
		$j++;
	}
	$j = 0;
	foreach my $cat (@lHeaders) {
		$xls_cov_null->write( $i, $j, $cat, $xls_export->get_format_header() );
		$j++;
	}
	
	my $h_lines;
	my $nb_line = 1;
	my $h_lines_low_cov;
	my $nb_line_low_cov = 1;
	my $h_lines_null_cov;
	my $nb_line_null_cov = 1;
	my $limit_down = 15;
	my $limit_up = 30;
	
	my @lPos = sort {$a <=> $b} keys %$h_res;
	foreach my $pos (@lPos) {
		my @_this_line;
		push(@_this_line, $pos);
		
		$h_lines->{$nb_line}->{format}->{0} = $xls_export->get_format_header();
		
		my $nb_pat = 0;
		my $total_pat = 0;
		foreach my $p_name (sort keys %{$h_res->{$pos}->{'all'}}) {
			$nb_pat++;
			$total_pat += $h_res->{$pos}->{'all'}->{$p_name};
		}
		my $mean_pat = $total_pat / $nb_pat;
		
		$nb_pat = 0;
		my ($has_low_cov, $has_null_cov);
		foreach my $p_name (sort keys %{$h_res->{$pos}->{'all'}}) {
			$nb_pat++;
			my $cov = $h_res->{$pos}->{'all'}->{$p_name};
			my $text = "$cov";
			push(@_this_line, $text);
			
			if ($cov == 0) {
				$h_lines->{$nb_line}->{format}->{$nb_pat} = $xls_export->get_format_red();
				$h_lines_low_cov->{$nb_line_low_cov}->{format}->{$nb_pat} = $xls_export->get_format_red();
				$h_lines_null_cov->{$nb_line_null_cov}->{format}->{$nb_pat} = $xls_export->get_format_red();
				$has_low_cov = 1;
				$has_null_cov = 1;
			}
			elsif ($cov < $limit_down) {
				$h_lines->{$nb_line}->{format}->{$nb_pat} = $xls_export->get_format_orange();
				$h_lines_low_cov->{$nb_line_low_cov}->{format}->{$nb_pat} = $xls_export->get_format_orange();
				$h_lines_null_cov->{$nb_line_null_cov}->{format}->{$nb_pat} = $xls_export->get_format_orange();
				$has_low_cov = 1;
			}
			elsif ($cov < $limit_up) {
				$h_lines->{$nb_line}->{format}->{$nb_pat} = $xls_export->get_format_blue();
				$h_lines_low_cov->{$nb_line_low_cov}->{format}->{$nb_pat} = $xls_export->get_format_blue();
				$h_lines_null_cov->{$nb_line_null_cov}->{format}->{$nb_pat} = $xls_export->get_format_blue();
			}
			else {
				$h_lines->{$nb_line}->{format}->{$nb_pat} = undef;
				$h_lines_low_cov->{$nb_line_low_cov}->{format}->{$nb_pat} = undef;
				$h_lines_null_cov->{$nb_line_null_cov}->{format}->{$nb_pat} = undef;
			}
		}
		
		my $text2 = sprintf("%.2f", $mean_pat);
		push(@_this_line, $text2);
		$h_lines->{$nb_line}->{line} = \@_this_line;
		$nb_line++;
		
		if ($has_low_cov) {
			$h_lines_low_cov->{$nb_line_low_cov}->{line} = \@_this_line;
			$nb_line_low_cov++;
		}
		if ($has_null_cov) {
			$h_lines_null_cov->{$nb_line_null_cov}->{line} = \@_this_line;
			$nb_line_null_cov++;
		}
		$has_low_cov = undef;
		$has_null_cov = undef;
	}
	
	foreach my $nb_line (sort {$a <=> $b} keys %{$h_lines}) {
		$i = $nb_line;
		$j = 0;
		foreach my $col (@{$h_lines->{$nb_line}->{line}}) {
			my $format = $h_lines->{$nb_line}->{format}->{$j};
			$xls_tr->write( $i, $j, $col, $format );
			$j++;
		}
	}
	
	foreach my $nb_line (sort {$a <=> $b} keys %{$h_lines_low_cov}) {
		$i = $nb_line;
		$j = 0;
		foreach my $col (@{$h_lines_low_cov->{$nb_line}->{line}}) {
			my $format = $h_lines_low_cov->{$nb_line}->{format}->{$j};
			$xls_cov_low->write( $i, $j, $col, $format );
			$j++;
		}
	}
	
	foreach my $nb_line (sort {$a <=> $b} keys %{$h_lines_null_cov}) {
		$i = $nb_line;
		$j = 0;
		foreach my $col (@{$h_lines_null_cov->{$nb_line}->{line}}) {
			my $format = $h_lines_null_cov->{$nb_line}->{format}->{$j};
			$xls_cov_null->write( $i, $j, $col, $format );
			$j++;
		}
	}
}


sub get_coverage_at_positions {
	my ($patient, $chr_name, $t, $ei) = @_;
	my $start = $ei->start();
	my $end = $ei->end();
	
	my $array = $patient->getNoSqlDepth->getDepth($chr_name, $start, $end);
	my $min = 99999;
	my $max = 0;
	my $total = 0;
	my $nb = 0;
	my $i = 0;
	foreach my $val (@$array) {
		my $this_pos = $start + $i;
		$i++;
		if ($only_coding) {
			next if not $t->getSpanCoding->contains($this_pos);
		}
		$min = $val if $val < $min;
		$max = $val if $val > $max;
		$total += $val;
		$nb++;
	}
	my $h;
	$h->{min} = $min;
	$h->{max} = $max;
	if ($total == 0 or $nb ==0) {
		$h->{min} = 0;
		$h->{max} = 0;
		$h->{mean} = 0;
	}
	else { $h->{mean} = sprintf("%.2f", ($total / $nb)); }
	return $h;
}

sub get_coverage_at_positions_start {
	my ($patient, $chr_name, $start) = @_;
	my $end = $start;
	
	my $array = $patient->getNoSqlDepth->getDepth($chr_name, $start, $end);
	my $min = 99999;
	my $max = 0;
	my $total = 0;
	my $nb = 0;
	my $i = 0;
	return $$array[0];
}




