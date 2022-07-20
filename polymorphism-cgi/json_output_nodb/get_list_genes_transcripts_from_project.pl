#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use Data::Dumper;
use JSON;
use Parallel::ForkManager;

use lib "$Bin/../GenBo/lib/obj-nodb/";
use GBuffer;
use GenBoProject;


my $cgi = new CGI();
print $cgi->header('text/json-comment-filtered');
print "{\"progress\":\".";

my $project_name	= $cgi->param('project');

my $buffer = GBuffer->new();
my $project = $buffer->newProject(-name=> $project_name);


my $h;
my $pm = new Parallel::ForkManager(3);
$pm->run_on_finish(
	sub {
		my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $hres ) = @_;
		unless ( defined($hres) or $exit_code > 0 ) {
			print qq|No message received from child process $exit_code $pid!\n|;
			warn Dumper $hres;
			return;
		}
		
		my $chr_id = $hres->{'chr'};
		foreach my $gene_id (sort keys %{$hres->{genes}}) {
			$h->{$chr_id}->{$gene_id} = $hres->{genes}->{$gene_id};
		}
	}
);
my @lCaptures = @{$project->getCaptures()};


my ($h_chr_found, $h_tr_captures_list);
foreach my $capture (@lCaptures) {
	foreach my $chr_id (keys %{$capture->getHashGenomicSpan()}) {
		$h_chr_found->{$chr_id} = undef;
	}
	if ($capture->transcripts_name()) {
		foreach my $tr_id (@{$capture->transcripts_name()}) {
			$h_tr_captures_list->{$tr_id} = undef;
		}
	}
}

foreach my $chr_id (keys %$h_chr_found) {
	my $chr = $project->getChromosome($chr_id);
	my $pid = $pm->start and next;
	my $hrestmp;
	$hrestmp->{'chr'} = $chr->id();
	my $i = 0;
	foreach my $gene (@{$chr->getGenes()}) {
		my $ok;
		foreach my $capture (@lCaptures) {
			my $intspan_gene = Set::IntSpan::Fast->new( $gene->start().'-'.$gene->end() );
			$intspan_gene = $intspan_gene->intersection($capture->getHashGenomicSpan->{$chr->id()});
			next if $intspan_gene->is_empty();
			$ok = 1;
			last;
		}
		next unless $ok;
		my @ltmp = split('_', $gene->id());
		my $ensg = $ltmp[0];
		my @lTr;
		foreach my $tr (@{$gene->getTranscripts()}) {
			my $ok_t;
			my @ltmp2 = split('_', $tr->id());
			my $ensgtr = $ltmp2[0];
			if ($h_tr_captures_list) {
				$ok_t = 1 if (exists $h_tr_captures_list->{$tr->id()});
				$ok_t = 1 if (exists $h_tr_captures_list->{$ensgtr});
			}
			else { $ok_t = 1; }
			next unless $ok_t;
			push(@lTr, $ensgtr);
		}
		next unless @lTr;
		$hrestmp->{genes}->{$ensg}->{external_name} = $gene->external_name();
		$hrestmp->{genes}->{$ensg}->{transcripts} = join(', ', sort @lTr);
		$hrestmp->{genes}->{$ensg}->{locus} = 'chr'.$chr->id().':'.$gene->start().'-'.$gene->end();
		$i++;
		print '.' if $i%10000; 
	}
	$pm->finish( 0, $hrestmp );
}
$pm->wait_all_children();

my $table_id = 'table_genes_tr';
my $html;
$html .= qq{<table id='$table_id' data-sort-name='name' data-filter-control='true' data-toggle="table" data-show-extended-pagination="true" data-cache="false" data-pagination-loop="false" data-total-not-filtered-field="totalNotFiltered" data-virtual-scroll="true" data-pagination-pre-text="Previous" data-pagination-next-text="Next" data-pagination="true" data-page-size="10" data-page-list="[10, 20, 50, 100, 200, 300]" data-resizable='true' class='table table-striped' style='font-size:9px;'>};
$html .= qq{<thead>};
$html .= qq{<th data-filter-control='input' data-sortable="true" data-field="locus"><b>LOCUS</b></th>};
$html .= qq{<th data-filter-control='input' data-sortable="true" data-field="ensg"><b>ENSG</b></th>};
$html .= qq{<th data-filter-control='input' data-sortable="true" data-field="name"><b>EXTERNAL NAME</b></th>};
$html .= qq{<th data-filter-control='input' data-sortable="true" data-field="transcript"><b>TRANSCRIPTS</b></th>};
$html .= qq{</thead>};
$html .= qq{<tbody>};

my @lExportCsv;
my @header;
push(@header, '#LOCUS');
push(@header, 'GENE_ENSG');
push(@header, 'GENE_NAME');
push(@header, 'TRANSCRIPTS');
push(@lExportCsv, \@header);
foreach my $chr_id (1..22, 'X', 'Y', 'M', 'MT') {
	next if not exists $h->{$chr_id};
	foreach my $ensg (sort keys %{$h->{$chr_id}}) {
		my $name = $h->{$chr_id}->{$ensg}->{external_name};
		my $transcripts = $h->{$chr_id}->{$ensg}->{transcripts};
		my $locus = $h->{$chr_id}->{$ensg}->{locus};
		$html .= qq{<tr style="font-size:11px;">};
		$html .= qq{<td>$locus</td>};
		$html .= qq{<td>$ensg</td>};
		$html .= qq{<td>$name</td>};
		$html .= qq{<td>$transcripts</td>};
		$html .= qq{</tr>};
		
		my @line;
		push(@line, $locus);
		push(@line, $ensg);
		push(@line, $name);
		push(@line, $transcripts);
		push(@lExportCsv, \@line);
	}
}
$html .= qq{</tbody>};
$html .= qq{</table>};

my $hash;
$hash->{html} = $html;
$hash->{table_id} = $table_id;
$hash->{description} = $project->description();
my @lCapturesNames;
foreach my $c (@{$project->getCaptures()}) {
	push(@lCapturesNames, $c->name());
}
$hash->{capture} = join(', ', sort @lCapturesNames);
$hash->{export_csv} = \@lExportCsv;
$hash->{gencode} = $project->gencode_version();

printJson($hash);
exit(0);

sub printJson {
	my ($hashRes) = @_;
	my $json_encode = encode_json $hashRes;
	print ".\",";
	$json_encode =~ s/{//;
	print $json_encode;
	exit(0);
}
