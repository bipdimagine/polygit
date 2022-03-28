package Cache_PolyDiag;
use strict;
use FindBin qw($RealBin);

use lib "$RealBin";
use lib "$RealBin/../";
use lib "$RealBin/../../cache/";
use Data::Dumper;
use Parallel::ForkManager;
use Bio::DB::Sam;
use JSON;
#use CacheDiag;
use Storable qw(store retrieve freeze thaw);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Cwd 'abs_path';
use Digest::MD5 qw(md5_hex);
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use lib "$RealBin/../../../GenBo/lib/obj-nodb/packages/";
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use POSIX qw(strftime);
use JSON;
use Compress::Snappy;
use preload_coverage;
use Digest::MD5::File qw(dir_md5_hex file_md5_hex url_md5_hex);



my $himpact_sorted = {
	"high" => "4",
	"moderate" =>"3",
	"low" =>"1",
};

sub cache_coverage {
	my ( $project_name, $fork ) = @_;
	warn "\n\n";
	warn "###################################################\n";
	warn "Prepare Cache for  COVERAGE !!!\n";
	warn "###################################################\n";
	cache_delete_coverage_file( $project_name, $fork );
	cache_coverage_primers( $project_name, $fork );
	cache_coverage_list_primers( $project_name, $fork );
	my $buffer1 = new GBuffer;
	my $projectP = $buffer1->newProject( -name => $project_name, -verbose => 1 );
	my $pm2           = new Parallel::ForkManager($fork);
	my @patients_name = map { $_->name() } @{ $projectP->getPatients };
	my $nbt           = scalar(@patients_name);
	return if $nbt eq 0;
	my $pr            = String::ProgressBar->new( "max" => $nbt );
	my $nbb           = 1;
	$pm2->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data ) =
			  @_;
			$pr->update( $nbb++ );
			$pr->write();
		}
	);
	my $diag;
	$diag = 1 if $projectP->isDiagnostic();
	unless ($diag) {
		foreach my $patient_name (@patients_name) {
			my $pid = $pm2->start() and next;
			compute_coverage_exome( $project_name, $patient_name );
			$pm2->finish(0);
		}
		$pm2->wait_all_children;
		return 1;
	}
	my $nb = 0;
	warn "coverage by patient  *******************\n";
	$pr = String::ProgressBar->new( max => scalar(@patients_name) );
	$nbb = 0;
	$pr->update(0);
	$pr->write();
	foreach my $patient_name (@patients_name) {
		my $pid = $pm2->start() and next;
		compute_coverage_diagnostic1( $project_name, $patient_name );
		$pm2->finish(0);
	}
	$pm2->wait_all_children;
	warn "now Primers *******************";
	$pr = String::ProgressBar->new( max => $nbt );
	$nbb = 0;
	$pr->update(0);
	$pr->write();
	foreach my $patient_name (@patients_name) {
		my $pid = $pm2->start() and next;
		compute_coverage_diagnostic2( $project_name, $patient_name );
		$pm2->finish(0);
	}
	$pm2->wait_all_children;
	warn "and finally ..  exons *****************************";
	$pr = String::ProgressBar->new( max => $nbt );
	$nbb = 0;
	$pr->update(0);
	$pr->write();
	my @utrs = (1);

	foreach my $patient_name (@patients_name) {
		$nb++;
		my $pid = $pm2->start() and next;
		compute_coverage_diagnostic3( $project_name, $patient_name, 0 );
		$pm2->finish(0);
	}
	$pm2->wait_all_children;
	print "\n";
}

sub compute_coverage_exome{
	my ($project_name,$patient_name) = @_;
		my $buffer = new GBuffer;
	my $project = $buffer->newProject( -name => $project_name, -verbose =>1 );
	my $patient = $project->getPatient($patient_name);
	my $no =  $project->noSqlCoverage();
	$no->clear($patient->name);
	$no->clear($patient->name."_cnv");
	my $patient = $project->getPatient($patient_name);
	preload_coverage::load_coverage_list_primers($project,[$patient],[],undef);
	$no->close();
	return;
}

sub cache_delete_coverage_file {
	my ( $project_name, $fork ) = @_;
	my $buffer1 = new GBuffer;
	my $projectP =
	  $buffer1->newProject( -name => $project_name, -verbose => 1 );
	my $no     = $projectP->noSqlCoverage();
	my $output = $projectP->getCacheDir() . "/coverage_lite";
	unlink $output . "/" . $projectP->name() . ".lite";
	unlink $output . "/primers.lite" if -e $output . "/primers.lite";
	my $no = $projectP->noSqlCoverage();
	$no->put( $projectP->name(), "test", "tutu" );
	$no->close();
}

sub cache_coverage_primers {
	my ( $project_name, $fork ) = @_;
	my $buffer1 = new GBuffer;
	my $projectP =
	  $buffer1->newProject( -name => $project_name, -verbose => 1 );
	my $pm2 = new Parallel::ForkManager($fork);
	my $no  = $projectP->noSqlCoverage();
	my $total;
	$pm2->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data ) =
			  @_;
			$total->{ $data->{chr} } = $data->{primers};
		}
	);
	foreach my $chr ( @{ $projectP->getChromosomes } ) {
		my $chr_name = $chr->name();
		my $pid      = $pm2->start() and next;
		my $buffer   = new GBuffer;
		my $project =
		  $buffer->newProject( -name => $project_name, -verbose => 1 );
		my $chr1    = $project->getChromosome($chr_name);
		my $primers = $chr1->getPrimers();
		map { delete $_->{project}; delete $_->{buffer} } @$primers;
		my $res;
		$res->{primers} = $primers;
		$res->{chr}     = $chr->name;
		$pm2->finish( 0, $res );
	}
	$pm2->wait_all_children;
	$no->put_bulk( "primers", $total );
	$no->close();
}

sub cache_coverage_list_primers {
	my ( $project_name, $fork ) = @_;
	my $buffer1 = new GBuffer;
	my $projectP =
	  $buffer1->newProject( -name => $project_name, -verbose => 1 );
	my $no = $projectP->noSqlCoverage();
	foreach my $capture ( @{ $projectP->getCaptures } ) {
		next if $capture->isPcr();
		my $list = $capture->getListPrimers();
		$no->put( $projectP->name(), $capture->id, $list );
	}
	$no->close();
}

sub compute_coverage_diagnostic1 {
	my ( $project_name, $patient_name ) = @_;
	my $buffer  = new GBuffer;
	my $project = $buffer->newProject( -name => $project_name, -verbose => 1 );
	my $patient = $project->getPatient($patient_name);
	my $no      = $project->noSqlCoverage();
	$no->clear( $patient->name );
	$no->clear( $patient->name . "_cnv" );
	my @trs =
	  map { $project->newTranscript($_) } @{ $project->bundle_transcripts() };
	preload_coverage::load_coverage_transcripts( $project, [$patient], \@trs );
	$no->close();
	undef($no);
}

sub compute_coverage_diagnostic2 {
	my ( $project_name, $patient_name ) = @_;
	my $buffer  = new GBuffer;
	my $project = $buffer->newProject( -name => $project_name, -verbose => 1 );
	my $patient = $project->getPatient($patient_name);
	my $no      = $project->noSqlCoverage();
	my @trs =
	  map { $project->newTranscript($_) } @{ $project->bundle_transcripts() };
	preload_coverage::load_coverage_primers( $project, [$patient], \@trs );
	$no->close();
	undef($no);
	return;
}

sub compute_coverage_diagnostic3 {
	my ( $project_name, $patient_name, $utr ) = @_;
	my $buffer   = new GBuffer;
	my $project  = $buffer->newProject( -name => $project_name, -verbose => 1 );
	my $patient  = $project->getPatient($patient_name);
	my $no       = $project->noSqlCoverage();
	my @paddings = ( 0, 5, 10, 15, 20, 30 );
	my @utrs     = ( 0, 1 );
	my @trs =
	  map { $project->newTranscript($_) } @{ $project->bundle_transcripts() };
	my $utr     = 0;
	my $padding = 0;
	preload_coverage::load_coverage_for_cache1( $project, [$patient], \@trs,
		$padding, $utr );
	$no->close();
	undef($no);
	return;
}

sub cache_polydiag {
	my ( $project_name, $fork ) = @_;
	my $buffer = new GBuffer;
	my $project = $buffer->newProject( -name => $project_name, -verbose => 1 );
	my $tbundle;
	foreach my $t ( @{ $project->bundle_transcripts() } ) {
		$tbundle->{$t}++;
	}
	my $nbt = scalar( @{ $project->getPatients } );
	my $pr  = String::ProgressBar->new( max => $nbt );
	my $nbb = 1;
	warn "###################################################\n";
	warn "Prepare Cache for  Polydiag !!!\n";
	warn "###################################################\n";
	$| = 1;
	my $pm = new Parallel::ForkManager($fork);
	my $error;
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump,
				$data_structure_reference )
			  = @_;
			$pr->update( $nbb++ );
			$pr->write();
			$error = 1 unless defined $data_structure_reference;
			die("\npb cache_polydiag") unless defined $data_structure_reference;
		}
	);

	my $project_name = $project->name();

	$pr->write();
	foreach my $p ( @{ $project->getPatients } ) {
		my $pname = $p->name();

		#next if $pname ne "AS1502552";
		my $pid = $pm->start and next;
		my $resp = run_cache_polydiag( $project_name, $p->name, $tbundle );
		$resp = $resp + 0;
		$pm->finish( 0, \$resp );
	}

	$pm->wait_all_children;
	$project->buffer->dbh->disconnect();
	#warn Dumper $project
	die() if $error;

}

sub run_cache_polydiag {
	my ( $project_name, $patient_name, $tbundle ) = @_;
	my $buffer1 = new GBuffer;
	my $project = $buffer1->newProject( -name => $project_name, -verbose => 1 );

	my $vquery = validationQuery->new(
		dbh          => $buffer1->dbh,
		capture_name => $project->validation_db()
	);
	my $pa      = $project->get_only_list_patients($patient_name);
	my $p       = $pa->[0];
	my $db_lite = $project->noSqlPolydiag("c");
	my $vtr;
	my %th;

	foreach my $chr ( @{ $project->getChromosomes } ) {

		#	my $db = $project->buffer->open_kyoto_db($file_out,'c');

		my $variations = $chr->getStructuralVariations();

		#	warn scalar(@$variations);
		#	die();
		my $ii = 0;
		my $dd = 0;
		foreach my $v ( @{$variations} ) {
			my $debug;

			#	$dd++;
			#warn $dd."/".scalar(@$variations);
			#$debug = 1 if $v->id eq "14_93670213_A_AT";
			#die() if $debug;
			#warn $ii++;

			my $ok;
			my $transcripts = $v->getTranscripts;
			foreach my $tr ( @{$transcripts} ) {

				next unless exists $tbundle->{ $tr->name };
				$ok = 1;
				my $h = construct_variant( $project, $v, $tr, $p,
					$vquery );
					CacheDiag::update_deja_vu($project,$tr,$h); 
					CacheDiag::update_cadd($project,$h); 
				my $id = join( ";", $tr->id, $v->id );
				$db_lite->put( $patient_name, $id, $h );

				#		$db->set($id,freeze $h );
				push( @{ $vtr->{ $tr->id } }, $v->id );
				$th{ $tr->id }++;
			}

			if ( $v->isIntergenic ) {
				my $h =
				  construct_intergenic_variant( $project, $v, $p,
					$vquery );
					CacheDiag::update_cadd($project,$h); 
					CacheDiag::update_deja_vu( $project, $v, $h );
				
				my $id = join( ";", "intergenic", $v->id );
				$db_lite->put( $patient_name, $id, $h );

				#	$db->set($id,freeze $h );
				push( @{ $vtr->{intergenic} }, $v->id );
				$th{"intergenic"}++;

			}
		}
		$project->purge_memory( $chr->length );
	}    #end chromosome

	foreach my $t ( keys %$vtr ) {
		$db_lite->put( $patient_name, "list_$t", join( ";", @{ $vtr->{$t} } ) );

		#$db->set("list_$t",join(";",@{$vtr->{$t}}));
	}
	$db_lite->put( $patient_name, "transcripts", join( ";", keys %th ) );

	#$db->set("transcripts",join(";",keys %th));
	my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
	my @days   = qw(Sun Mon Tue Wed Thu Fri Sat Sun);

	my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) =
	  localtime();
	my $date = $mday . "/" . ( $mon + 1 ) . "/" . ( 1900 + $year );

	#$db->set("date",$date);
	$db_lite->put( $patient_name, "date", $date );
	my $f1s = $p->getVariationsFiles();
	my $tf;
	foreach my $f1 (@$f1s) {
		if ( -e $f1 ) {
			$tf->{$f1} = file_md5_hex($f1);

			#push(@$tf,file_md5_hex($f1));

			#$db->set("variations_vcf_md5",file_md5_hex($f1)) ;
		}
	}
	$db_lite->put( $patient_name, "variations_vcf_md5", $tf );
	$tf = {};
	my $f2s = $p->getIndelsFiles();
	foreach my $f1 (@$f2s) {
		$tf->{$f1} = file_md5_hex($f1);

		# 			push(@$tf,file_md5_hex($f1)) if -e $f1;;
	}

	$db_lite->put( $patient_name, "indels_vcf_md5", $tf );

	#$db->set("indels_vcf_md5",file_md5_hex($f2)) if -e $f2;
	#$db->close();
	#warn "end run cache";
	$db_lite->close();
	$project->buffer->dbh->disconnect();
	if ( exists $project->{cosmic_db} ) {
		$project->{cosmic_db}->close();
	}
	
#	warn Dumper keys %{$project};
	
	$project = undef;

	#$project->{cosmic_db}->close();

	return 1;
}

sub cache_cnv {
	my ( $project_name, $fork ) = @_;
	my $buffer1 = new GBuffer;
	my $projectP =
	  $buffer1->newProject( -name => $project_name, -verbose => 1 );
	my @transcripts_cgi = @{ $projectP->bundle_transcripts() };
	my @transcripts =
	  sort { $a->getGene->external_name cmp $b->getGene->external_name }
	  map  { $projectP->newTranscript($_) } @transcripts_cgi;
	warn "###################################################\n";
	warn "Prepare Cache for  CNV  !!!\n";
	warn "###################################################\n";
	my $primers = $projectP->getPrimers();

	foreach my $patient ( @{ $projectP->getPatients } ) {
		foreach my $primer (@$primers) {
			$primer->cached_cnv($patient);
		}
	}
	preload_coverage::load_cnv_score( $projectP, $projectP->getPatients,
		\@transcripts );
	warn "###################################################\n";
	warn "END CNV  !!!\n";
	warn "###################################################\n";
	warn "\n";
}

sub construct_variant {
	my ($project,$v,$tr1,$patient,$vquery) = @_;
	my $hvariation;
	$hvariation->{id} = $v->id;
	if ($project->isSomatic){
		
		$hvariation->{cosmic} = $v->cosmic();
		my $id  = $v->cosmic();
		if ($id){
			my ($id1,$nb) = split(":",$id);
			my $text = $id1;
			my $url = "http://cancer.sanger.ac.uk/cosmic/mutation/overview?id="; 
			 $url = "http://grch37-cancer.sanger.ac.uk/cosmic/ncv/overview?id="  if ($id1 =~/COSN/);
			  $text = "$id1 [$nb]" if $nb;
			$id1 =~s/COSM//;
			$id1 =~s/COSN//;
		
			
			
			
			
			$hvariation->{cosmic} = qq{<a href=\"$url$id1\" target="_blank">$text</a>};
		}
	}
		
		$hvariation->{impact_text} = $v->effectImpact($tr1);
		$hvariation->{impact_score} = $himpact_sorted->{$v->effectImpact($tr1)};
		$hvariation->{gene} = $tr1->getGene->external_name();
		$hvariation->{var_name} = $v->name();
		if ($v->name() =~ /rs/){
			my $vn = $v->name();
			$hvariation->{var_name} = qq{<a href='http://www.ncbi.nlm.nih.gov/snp/?term=$vn' target = '_blank'>$vn</a> };	
		}
		my @asequence_info;
		my @apc;
		my @methods;
		my $nb_methods;
		foreach my $method (@{$patient->callingMethods}){
			next unless exists $v->annex()->{$patient->id}->{method_calling}->{$method}->{nb_all_ref};
			
			my $all_annex = $v->annex()->{$patient->id}->{method_calling}->{$method};
			my $nb_ref =$all_annex->{nb_all_ref};
			my $nb_alt =  $all_annex->{nb_all_mut};
			
		my $method_name = substr $method,0,3;
		push(@methods,$method_name);
		my $sequence_info = "he("; 
		my $pc ="-";		
		if ($v->annex()->{$patient->id}->{nb_all_ref} eq "?"){
			$sequence_info = "??";
		}
		else {
		$sequence_info = "ho(" if $all_annex->{ho};
		
		my $sum = $nb_ref + $nb_alt;
		if ($sum >0){
		 $pc = int ($nb_alt *100/($sum));
		}
		$sequence_info .= $nb_ref."/".$nb_alt.")";
	
		}
		$sequence_info = $method_name.":".$sequence_info;
		$pc = $method_name.":".$pc."%";
		push(@apc,$pc);
		push(@asequence_info,$sequence_info);
		$nb_methods ++;
		}
		
		 if ($v->validation_method eq "sanger" ) {
		 	#$sequence_info = "-";
		 	push(@asequence_info,"-");
		 }
		
	
		$hvariation->{ngs} = join("<br>",@asequence_info);
	
		$hvariation->{ratio} =  join("<br>",@apc);
		$hvariation->{caller} =  join("<br>",@methods);
		$hvariation->{allele} = $v->getSequence();
		$hvariation->{ref_allele} = $v->ref_allele();
		$hvariation->{genomique} = $v->getChromosome()->name.":".$v->start;
		$hvariation->{start} = $v->start;
		$hvariation->{transcript} = $tr1->name;
		$hvariation->{chromosome} = $v->getChromosome()->name;
		$hvariation->{trans} = $v->start * $tr1->strand; 
		
		$hvariation->{cds} = "";
		$hvariation->{prot} ="";
		$hvariation->{codons_AA} = "";
		$hvariation->{polyphen} ="-";
		$hvariation->{sift} ="-";

		my $debug;
		$debug = 1   if $v->name eq "rs113993960";
		#warn  $v->delete_sequence if $debug;
	
		if ($v->isDeletion){
	
			$hvariation->{codons}  =  $v->delete_sequence."/".$v->sequence();
		}
		else {
			$hvariation->{codons}  =  $v->getChromosome()->sequence($v->start,$v->end)."/".$v->sequence();
		}
		if ($tr1->strand() == -1 ){
			$hvariation->{codons}  =  BioTools::complement_sequence($v->getChromosome()->sequence($v->start,$v->end))."/".BioTools::complement_sequence($v->sequence());
			}
			
		my $start = $v->start;
		my $chr = $v->getChromosome()->name();
		my $vid = $hvariation->{id};
		my $a0 = $v->ref_allele;
		my $a1 = $v->var_allele;
	
		my $pname = $patient->name;
		
			my $qq4 = qq{	<div  data-dojo-type="dijit/Toolbar"><button dojoType="dijit.form.Button"   iconClass="igvIcon" onClick ="displayInIGV('$chr',$start,$start);"></button></div>};
			my $qq4 = qq{	<button  class="igvIcon2" onClick ="displayInIGV('$chr',$start,$start);">&nbsp;&nbsp;&nbsp;&nbsp;</button>};
			$hvariation->{igv} = $qq4; 		
			
			my $qq5 = qq{	<div  data-dojo-type="dijit/Toolbar"><button dojoType="dijit.form.Button"   iconClass="alamutView" onClick ="displayInAlamut('$chr',$start,['$a0','$a1']);"></button></div>};
			my $qq5 = qq{	<button    class="alamutView3" onClick ="displayInAlamut('$chr',$start,['$a0','$a1']);"></button>};
			$hvariation->{alamut} = $qq5; 
						
						my $qq3 = qq{
      				<div  data-dojo-type="dijit/Toolbar">
					<button dojoType="dijit.form.Button"   iconClass="dijitEditorIcon dijitEditorIconFullScreen" onClick = viewElectro('$pname','$vid')></button>
					</div>
      					};
      					my $qq3 = qq{
      				
					<button  class="alignIcon" onClick = viewElectro('$pname','$vid')></button>
			
      					};	
			$hvariation->{align} = $qq3; 	
		
		
		if ($v->isCoding($tr1)){
			my $prot = $tr1->getProtein();
		
			$hvariation->{cds} = $v->getOrfPosition($prot);
			$hvariation->{prot} = $v->getProteinPosition($prot);
			$hvariation->{codons_AA} =   $v->protein_nomenclature($prot);#$v->getProteinAA($prot).$hvariation->{prot}.$v->changeAA($prot);
			$hvariation->{polyphen} = $v->polyphenScore($tr1->getProtein);
			$hvariation->{sift} = $v->siftScore($tr1->getProtein);
			$hvariation->{codons} =   $v->getCodons($tr1);
		}
		#warn $hvariation->{codons} if $debug;
			#			die if $debug;
		#warn $hvariation->{prot};
		$hvariation->{exon} = $tr1->findExonNumber($v->start);
		$hvariation->{exon} = $tr1->findNearestExon($v->start) if $hvariation->{exon} == -1;
		$hvariation->{nomenclature} =  $v->getNomenclature($tr1);
		$hvariation->{consequence} =  $v->variationType($tr1);
		
		$hvariation->{freq}  =  $v->frequency;
		$hvariation->{freq}  =  0 if $hvariation->{freq} == -999;
		#$hvariation->{freq} = $hvariation->{freq}/100;
		$hvariation->{scaled_score} = $v->scaledScoreVariant($tr1,$patient,$vquery);
		if($nb_methods == 1 && $hvariation->{ngs} =~/dup/){
			$hvariation->{dup} = 1;
		}
		$hvariation->{score} = $v->scoreVariant($tr1,$patient,$vquery);
		
				
		if ($hvariation->{freq} < 0 ){
							$hvariation->{freq_score} = 4; 
							$hvariation->{freq_level} = 1; 
							$hvariation->{freq} = 0;
						
		}
		
		if ($hvariation->{freq} <= 0.01){
							$hvariation->{freq_score} = 4; 
							$hvariation->{freq_level} = 1; 
		}
		
		elsif ($hvariation->{freq} <= 0.05){
				$hvariation->{freq_score} = 2; 
					$hvariation->{freq_level} = 3; 
		}
		else {
					$hvariation->{freq_score} = 1; 
					$hvariation->{freq_level} = 4; 
			}
		$hvariation->{freq} = "-" unless $v->frequency();
		return $hvariation;
}

sub update_deja_vu{
	 my ($project,$tr,$hvariation) = @_;
	#my $dejavu =$tr->getChromosome->hash_kyoto_dejavu;

	my $vid = $hvariation->{id};
	my $debug;
	$debug=1 if $vid eq "14_93670213_A_AT";
	
	warn $vid  if $debug;;

	my $similar = $project->similarProjects();
	my $hres = $project->getDejaVuInfosForDiag($vid);

	my $nb_pat =0;
	my $pname = $project->name();
	my $proj_dejavu =0;
	$hvariation->{sim_deja_vu} = $hres->{similar_patient};
	$hvariation->{sim_proj_deja_vu} = $hres->{similar_project};
	
	
		my $nb_dejavu = $project->getDejaVuThisProject($vid);
		die() if $nb_dejavu eq 0;
		$hvariation->{this_deja_vu} =$nb_dejavu;
		$hvariation->{diff_project_deja_vu} =$hres->{other_project}  ; 
		$hvariation->{project_deja_vu} = $hres->{other_project} + $hres->{similar_project}; 
		$hvariation->{deja_vu} =$hres->{other_project}.":".$hres->{other_patient}; 
		$hvariation->{in_this_run} = $nb_dejavu."/". scalar(@{$project->getPatients});
		$hvariation->{similar_projects} = $hres->{similar_project}.":". $hres->{similar_patient};
		
		my $freq_score;
		my $freq_level;
		
		if ($hvariation->{freq} < 0 ){
							$hvariation->{freq_score} = 4; 
							$hvariation->{freq_level} = 1; 
							$hvariation->{freq} = 0;
						
		}
		
		if ($hvariation->{freq} <= 0.01){
							$hvariation->{freq_score} = 4; 
							$hvariation->{freq_level} = 1; 
		}
		
		elsif ($hvariation->{freq} <= 0.05){
				$hvariation->{freq_score} = 2; 
					$hvariation->{freq_level} = 3; 
		}
		else {
					$hvariation->{freq_score} = 1; 
					$hvariation->{freq_level} = 4; 
			}
	
			
		if ($hvariation->{freq_score} == 4) {
				if ($hvariation->{diff_project_deja_vu} == 0){
					$freq_score = 4; 
					$freq_level =1;
					
				}
				elsif ($hvariation->{diff_project_deja_vu} <= 20 ){
							$freq_score = 3; 
							$freq_level =2;
					
				}
				elsif ($hvariation->{diff_project_deja_vu} <= 50){
						$freq_score = 2; 
						$freq_level =3;
					
				}
				else {
						$freq_score = 1;
						$freq_level =4; 
				}
		}	
		
		$hvariation->{freq_score} = $freq_score if $freq_score < $hvariation->{freq_score} && $freq_score;
		$hvariation->{freq_level} = $freq_level if $freq_level > $hvariation->{freq_level} && $freq_level ;
		
}

sub construct_intergenic_variant {
	my ($project,$v,$patient,$vquery) = @_;

	my $hvariation;
	$hvariation->{id} = $v->id;


		
		$hvariation->{impact_text} = 0;
		$hvariation->{impact_score} = 0;
		$hvariation->{gene} = "intergenic";
		$hvariation->{var_name} = $v->name();
		if ($v->name() =~ /rs/){
			my $vn = $v->name();
			$hvariation->{var_name} = qq{<a href='http://www.ncbi.nlm.nih.gov/snp/?term=$vn' target = '_blank'>$vn</a> };	
		}
		
			my $sequence_info = "he("; 
		my $pc ="-";		
		if ($v->annex()->{$patient->id}->{nb_all_ref} eq "?"){
			$sequence_info = "??";
		}
		else {

		$sequence_info = "ho(" if $v->annex()->{$patient->id}->{ho};
		
		my $sum = $v->annex()->{$patient->id}->{nb_all_ref}+$v->annex()->{$patient->id}->{nb_all_mut};
		if ($sum >0){
		 $pc = int ($v->annex()->{$patient->id}->{nb_all_mut} *100/($sum));
		}
		$sequence_info .= $v->annex()->{$patient->id}->{nb_all_ref}."/".$v->annex()->{$patient->id}->{nb_all_mut}.")<br>";
	
		}
		
		 if ($v->validation_method eq "sanger" ) {
		 	$sequence_info = "-";
		 }
		$hvariation->{ngs} = $sequence_info;
		$hvariation->{ratio} = $pc."%";
		$hvariation->{allele} = $v->getSequence();
		$hvariation->{ref_allele} = $v->ref_allele();
		$hvariation->{genomique} = $v->getChromosome()->name.":".$v->start;
		$hvariation->{start} = $v->start;
		$hvariation->{transcript} = "-";
		$hvariation->{chromosome} = $v->getChromosome()->name;
		$hvariation->{trans} = $v->start ;
		$hvariation->{cds} = "";
		$hvariation->{prot} ="";
		$hvariation->{codons_AA} = "";
		$hvariation->{polyphen} ="-";
		$hvariation->{sift} ="-";

		$hvariation->{codons}  =  $v->getChromosome()->sequence($v->start,$v->end)."/".$v->sequence();
	
		
		my $start = $v->start;
		my $chr = $v->getChromosome()->name();
		my $vid = $hvariation->{id};
		my $a0 = $v->ref_allele;
		my $a1 = $v->var_allele;

		my $pname = $patient->name;
		
			my $qq4 = qq{	<div  data-dojo-type="dijit/Toolbar"><button dojoType="dijit.form.Button"   iconClass="igvIcon" onClick ="displayInIGV('$chr',$start,$start);"></button></div>};
			my $qq4 = qq{	<button  class="igvIcon2" onClick ="displayInIGV('$chr',$start,$start);">&nbsp;&nbsp;&nbsp;&nbsp;</button>};
			$hvariation->{igv} = $qq4; 		
			
			my $qq5 = qq{	<div  data-dojo-type="dijit/Toolbar"><button dojoType="dijit.form.Button"   iconClass="alamutView" onClick ="displayInAlamut('$chr',$start,['$a0','$a1']);"></button></div>};
			my $qq5 = qq{	<button    class="alamutView3" onClick ="displayInAlamut('$chr',$start,['$a0','$a1']);"></button>};
			$hvariation->{alamut} = $qq5; 
						
						my $qq3 = qq{
      				<div  data-dojo-type="dijit/Toolbar">
					<button dojoType="dijit.form.Button"   iconClass="dijitEditorIcon dijitEditorIconFullScreen" onClick = viewElectro('$pname','$vid')></button>
					</div>
      					};
      					my $qq3 = qq{
      				
					<button  class="alignIcon" onClick = viewElectro('$pname','$vid')></button>
			
      					};	
			$hvariation->{align} = $qq3; 	
		
		
		
		
		#warn $hvariation->{prot};
		$hvariation->{exon} = "-";
		$hvariation->{exon} = "-";
		$hvariation->{nomenclature} =  "-";
		$hvariation->{consequence} =  "-";
		
		$hvariation->{freq}  =  $v->frequency;
		$hvariation->{freq}  =  0 if $hvariation->{freq} == -999;
		#$hvariation->{freq} = $hvariation->{freq}/100;
		$hvariation->{scaled_score} = 1;
		$hvariation->{impact_score} =1;
		$hvariation->{scaled_score} =1;
		my $debug ;
		$debug =1 if $v->name eq "4_86916249_C_T";
		warn $v->frequency if $debug;
		if ($hvariation->{freq} < 0 ){
							$hvariation->{freq_score} = 0; 
							$hvariation->{freq_level} = 1; 
							$hvariation->{freq} = 0;
						
		}
		
		if ($hvariation->{freq} <= 0.01){
							$hvariation->{freq_score} = 4; 
							$hvariation->{freq_level} = 1; 
		}
		
		elsif ($hvariation->{freq} <= 0.05){
				$hvariation->{freq_score} = 2; 
					$hvariation->{freq_level} = 3; 
		}
		else {
					$hvariation->{freq_score} = 1; 
					$hvariation->{freq_level} = 4; 
			}
				
		$hvariation->{freq} = "-" unless $v->frequency;
		
		return $hvariation;
}



1;
