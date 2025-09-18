#!/usr/bin/perl
$| = 1;
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../GenBo";
use lib "$Bin/../../GenBo/lib/obj-nodb";
use lib "$Bin/../../GenBo/lib/obj-nodb/polyviewer/";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
require
  "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/update_variant_editor.pm";

use connect;
use GBuffer;
use Getopt::Long;
use Data::Dumper;
use JSON;
use xls_export;
use session_export;
use List::MoreUtils qw{ natatime };
use Parallel::ForkManager;
use PolyviewerJunction;
use polysplice_html;

my $cgi                  = new CGI;
my $project_name         = $cgi->param('project');
my $patient_name         = $cgi->param('patient');
my $fork         = $cgi->param('fork');
#$only_positions = '10:17000000-17080000';


my $buffer  = GBuffer->new;
my $project = $buffer->newProjectCache( -name => $project_name );
$project->getChromosomes();
$project->getFamilies();
my $patient = $project->getPatient($patient_name);

$fork =6 unless $fork;


my ( @lJunctions, $h_var_linked_ids );

add_linked_hash_in_cache ($patient );

my $cmd = "$Bin/rna_junctions_patient.pl project=$project_name patient=$patient_name dejavu=10 dejavu_percent=100 min_score=10 only_dejavu_ratio_10=1 only_junctions_NDA=1 only_junctions_A=1 only_junctions_D=1";
`$cmd`;

exit(0);




sub add_linked_hash_in_cache {
	my ( $patient ) = @_;
	my $h_vector;
	$patient->getProject->disconnect();
	my $fork     =  $fork;
my $pm       = new Parallel::ForkManager($fork);
my $nbErrors = 0;
my $global ={};
my $global_linked = {};
$pm->run_on_finish(
	sub {
		my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $hres ) = @_;
		unless ( defined($hres) or $exit_code > 0 ) {
			$nbErrors++;
			print qq|No message received from child process $exit_code $pid!\n|;
			return;
		}
		$global->{$hres->{chr}} = $hres->{vector};
		foreach my $jid (keys %{$hres->{linked}} ){
			$global_linked->{$jid} = $hres->{linked}->{$jid};
		}
		}
);
	
	
	
	foreach my $chr ( @{ $patient->getProject->getChromosomes() } ) {
		$pm->start and next;
		my $h_var_linked_ids = {};
		my $hres;
		my $h_vector_chr;
		
		my $vector_patient = $patient->getJunctionsVector($chr);
		
		$h_vector_chr->{min0} = $vector_patient->Clone();
		$h_vector_chr->{min2} = $vector_patient->Clone();
		$h_vector_chr->{min4} = $vector_patient->Clone();
		$h_vector_chr->{min6} = $vector_patient->Clone();
		$h_vector_chr->{min8} = $vector_patient->Clone();
		foreach my $junction ( @{ $chr->getListVarObjects($vector_patient) } ) {
			
			if ( $junction->is_junctions_linked($patient) ) {
				if ( exists $h_var_linked_ids->{ $junction->id() } ) {
					$h_var_linked_ids->{ $junction->id() }->{vector_id} = $chr->id() . '-' . $junction->vector_id();
				}
				else {
					$h_var_linked_ids->{ $junction->id() }->{vector_id} =
					  $chr->id() . '-' . $junction->vector_id();
					my $hdone;
					my @lOthers = (
						keys %{
							$junction->get_hash_junctions_linked_to_me->{ $patient->name() }
						}
					);
					my $nb_others = scalar(@lOthers);
					while ( $nb_others > 0 ) {
						my $other_id = $lOthers[0];

						#if (not exists $hdone->{$other_id}) {
						foreach my $id (keys %{ $junction->get_hash_junctions_linked_to_me->{ $patient->name()} }){
							$h_var_linked_ids->{ $junction->id() }->{linked_to}->{$id} = undef;
							$h_var_linked_ids->{$id}->{linked_to}->{ $junction->id() } = undef;
							push( @lOthers, $id ) if not( exists $hdone->{$other_id} );
						}
						$hdone->{$other_id} = undef;
						my $supress = shift(@lOthers);
						$nb_others = scalar(@lOthers);
					}

					my $hallids = $h_var_linked_ids->{ $junction->id() }->{linked_to};
					$hallids->{ $junction->id() } = undef;
					foreach my $jid ( keys %$hallids ) {
						next unless exists $h_var_linked_ids->{$jid}->{linked_to};
						$h_var_linked_ids->{$jid}->{linked_to} = $hallids;
					}
				}
			}
			my $score_this_j = $junction->junction_score_without_dejavu_global($patient);
			
			
			if ( $score_this_j < 0 ) {
				$h_vector_chr->{min0}->Bit_Off( $junction->vector_id() );
				$h_vector_chr->{min2}->Bit_Off( $junction->vector_id() );
				$h_vector_chr->{min4}->Bit_Off( $junction->vector_id() );
				$h_vector_chr->{min6}->Bit_Off( $junction->vector_id() );
				$h_vector_chr->{min8}->Bit_Off( $junction->vector_id() );
			}
			if ( $score_this_j < 2 ) {
				$h_vector_chr->{min2}->Bit_Off( $junction->vector_id() );
				$h_vector_chr->{min4}->Bit_Off( $junction->vector_id() );
				$h_vector_chr->{min6}->Bit_Off( $junction->vector_id() );
				$h_vector_chr->{min8}->Bit_Off( $junction->vector_id() );
			}
			if ( $score_this_j < 4 ) {
				$h_vector_chr->{min4}->Bit_Off( $junction->vector_id() );
				$h_vector_chr->{min6}->Bit_Off( $junction->vector_id() );
				$h_vector_chr->{min8}->Bit_Off( $junction->vector_id() );
			}
			if ( $score_this_j < 6 ) {
				$h_vector_chr->{min6}->Bit_Off( $junction->vector_id() );
				$h_vector_chr->{min8}->Bit_Off( $junction->vector_id() );
			}
			if ( $score_this_j < 8 ) {
				$h_vector_chr->{min8}->Bit_Off( $junction->vector_id() );
			}
		}
		$h_vector_chr->{min0} = $h_vector_chr->{min0}->to_Enum();
		$h_vector_chr->{min2} = $h_vector_chr->{min2}->to_Enum();
		$h_vector_chr->{min4} = $h_vector_chr->{min4}->to_Enum();
		$h_vector_chr->{min6} = $h_vector_chr->{min6}->to_Enum();
		$h_vector_chr->{min8} = $h_vector_chr->{min8}->to_Enum();
		$hres->{done} = 1;
		$hres->{chr} = $chr->id;
		$hres->{vector} = $h_vector_chr;
		$hres->{linked} = $h_var_linked_ids;
		$pm->finish( 0, $hres );
		
	}
	$pm->wait_all_children();
	my $no_cache    = $patient->get_lmdb_cache("c");
	warn $no_cache->filename();
	my $outfile_log = $no_cache->filename() . '.ok';
	my $cache_id = 'splices_linked_' . $patient->name();
	$no_cache->put_cache_hash( $cache_id, $global_linked ) if $global_linked;
	$no_cache->put_cache_hash(  $patient->name() . '_chr_vectors_enum', $global );
	$no_cache->close();
		my $cmd = 'touch ' . $outfile_log;
		`$cmd`;
}


sub printJson {
	my ($hashRes) = @_;
	my $json_encode = encode_json $hashRes;
	print ".\",";
	$json_encode =~ s/{//;
	print $json_encode;
	exit(0);
}
