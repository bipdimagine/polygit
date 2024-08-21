package Cache_nodb;
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use Data::Dumper;
use Data::Printer;
use Parallel::ForkManager;
use Bio::DB::Sam;
use CacheGenesData_nodb; 
use CacheVariationsData;
use Storable qw(store retrieve freeze thaw);
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use List::MoreUtils qw{ natatime };
use String::ProgressBar;
use Digest::MD5::File qw(dir_md5_hex file_md5_hex url_md5_hex);
#use List::Compare::Functional qw(:main :mainrefs);
#use Hash::Merge  qw( merge );
use CacheDiag;
use Math::Combinatorics;
use GenBoNoSql;

sub forkCache {
	my ($project,$maxProc) = @_;

	
	my $project_id = $project->id();
	
	my $pm = new Parallel::ForkManager($maxProc);

	
	my $project_name = 	$project->name();
	#my $database = $project->getDataBase;
	my $chrs = $project->getChromosomes();
	
	my @chr_names = sort{$a <=> $b} map {$_->name} @{$project->getChromosomes()};
	
	$project->buffer->dbh->disconnect();
	
	foreach my $chr (@chr_names) {
		 #next if $chr ne "1";
			my $pid = $pm->start and next;
			my $resp = start_thread($project_name,$chr);
			$pm->finish();	
	}	
	
	warn "###################################################\n";
	warn "WAIT !!!\n";
	warn "###################################################\n";	
	$pm->wait_all_children;
	warn "###################################################\n";
	warn "Saving in database  !!!\n";
	warn "###################################################\n";	
	###
	# saving in database 
	###
	#exit(0);
	my $buffer2 = new GBuffer;
	my $project_end = $buffer2->newProject( -name => $project_name );
	my $dir_store = $project_end->getRootDir()."/";
	my $data;
	$data->{type} = "nodb";
	my $homo;
	my $hetero;
	my $homo_file;
	my $hetero_file;
	 my $dir_cache  = $project_end->getCacheDir();
	
	foreach my $p (@{$project_end->getPatients}){
		push(@{$data->{patients}},$p->name());  
		my $file1 = $dir_cache."/".$p->name.".homo.kct";
		my $file2 = $dir_cache."/".$p->name.".hetero.kct";
		 $homo_file->{$p->name()} = $file1;
		 $hetero_file->{$p->name()} = $file2;
		my $db = new KyotoCabinet::DB;
		if (!$db->open($file1, $db->OWRITER|$db->ONOLOCK | $db->OCREATE|$db->OTRUNCATE )) {
     			printf STDERR ("open error: %s\n", $db->error);
     			die();
 		}
 		system("chmod a+w $file1");
 		$homo->{$p->name()} = $db;
 		$db->close();
 		my $db2 = new KyotoCabinet::DB;
 		if (!$db2->open($file2, $db2->OWRITER|$db2->ONOLOCK | $db2->OCREATE|$db2->OTRUNCATE )) {
     			printf STDERR ("open error: %s\n", $db2->error);
     			die();
 		}
 		
 		$hetero->{$p->name()} = $db2;
 		$db2->close();
 		system("chmod a+w $file2");
 		 
		#tie(%{$homo->{$p->name()}}, 'KyotoCabinet::DB', $file1 , KyotoCabinet::DB::ONOLOCK | KyotoCabinet::DB::OWRITER | KyotoCabinet::DB::OCREATE| KyotoCabinet::DB::OTRUNCATE) || die($file1);
	#	tie(%{$hetero->{$p->name()}}, 'KyotoCabinet::DB', $file2 , KyotoCabinet::DB::ONOLOCK | KyotoCabinet::DB::OWRITER | KyotoCabinet::DB::OCREATE| KyotoCabinet::DB::OTRUNCATE) || die($file2);
	}
	
	 my $db;
	
	 
	foreach my $chr (@{$project_end->getChromosomes()}) {	
		#next unless $chr->getVariations();
		my $file_store = $dir_store."/".$chr->name.".".$project_end->getDataBase.".store";
		unless (-e $file_store) {
			warn $file_store;
			die("!!!!! miss $file_store") if $chr->getVariations();
			next;
		}
	
		my $dref = retrieve($file_store);	
		my $toto = $dref->{data}->{coverage};
      	my $titi =  $dref->{data}->{patients};
      	delete ($dref->{data}->{coverage});
      	delete ($dref->{data}->{patients});
      	my $kyoto_cahe_genes = $project->getCacheGenesKyotoFile($chr->name);
    
   	
		
      	$data->{chromosomes}->{$dref->{name}}->{genes}=$dref->{data};
      	$data->{chromosomes}->{$dref->{name}}->{coverage} = $toto;
      	#$data->{chromosomes}->{$dref->{name}}->{patients} = $titi;

      	foreach my $patient (keys %$titi){
      		my $db_homo = new KyotoCabinet::DB;
      		if (!$db_homo->open($homo_file->{$patient} , $db_homo->OWRITER|$db_homo->ONOLOCK )) {
     			printf STDERR ("open error: %s\n", $db_homo->error);
     			die();
 		}
 		my $db_hetero = new KyotoCabinet::DB;
      		if (!$db_hetero->open($hetero_file->{$patient}, $db_hetero->OWRITER|$db_hetero->ONOLOCK )) {
     			printf STDERR ("open error: %s\n", $db_hetero->error);
     			die();
 		}
      			if (exists $titi->{$patient}->{homozygote}){
      			foreach my $v (@{$titi->{$patient}->{homozygote}}){
      				$db_homo->set($v,1);
      				
      			}
      			}
      			if (exists $titi->{$patient}->{heterozygote}){
      			foreach my $v (@{$titi->{$patient}->{heterozygote}}){
      				$db_hetero->set($v,1);
      				
      			}
      			}
      			$db_hetero->close();
      			$db_homo->close();
      	}
  	
     
		
	}
	warn  $project_end->getCacheDir();
	
	$data->{type_cache} = 1;
	my $file_cache_genes = $project_end->getCacheGenesFile();
	$project_end->buffer->saveStore($data,$file_cache_genes);
	system("chmod a+rw $file_cache_genes");
	 #GenBoStorable::insertStorable( $project_end->buffer->dbh, $project_end->id, $project_end->id,"genes", $data);	
	
	warn "\t => end save  cache ";
	foreach my $chr (@{$project_end->getChromosomes()}) {
		my $store_variations_file = $dir_store."/".$chr->name.".".$project_end->getDataBase.".var.kct";
		my $kyoto_gene_file =   $dir_store."/".$chr->name.".genes.kct";
		my $kyoto_filter_genes =   $dir_store."/".$chr->name.".genes.filter.kct";
			
		unless (-e $store_variations_file) {
			die("miss $store_variations_file") if $chr->getVariations();
			next;
		}
		
	
		my $kyoto_cache_variation = $project->getCacheVariationsKyotoFile($chr->name);
		my $kyoto_cache_genes = $project->getCacheGenesKyotoFile($chr->name);
		my $kyoto_cache_filter_genes = $project->getCacheFilterGenesKyotoFile($chr->name);
		
		#warn "mv $store_variations_file $kyoto_cache_variation";
		system ("mv $store_variations_file $kyoto_cache_variation");
		system("chmod a+rw $kyoto_cache_variation");
		
#		warn "mv $kyoto_gene_file $kyoto_cache_genes";
			system ("mv $kyoto_gene_file $kyoto_cache_genes");
		system("chmod a+rw $kyoto_cache_genes");
		
		
#		warn "mv $kyoto_filter_genes $kyoto_cache_filter_genes";
		system ("mv $kyoto_filter_genes $kyoto_cache_filter_genes");
		system("chmod a+rw $kyoto_cache_filter_genes");

			
			
			
			}
		
	
		#GenBoStorable::insertStorable( $project_end->buffer->dbh, $project_end->id, $chr->id,"variations", $dref);
	
	 
	warn "###################################################\n";
	warn "Prepare Sample kct  file !!!\n";
	warn "###################################################\n";	
	
	foreach my $patient (@{$project_end->getPatients()}){
		warn "\t patient: ".$patient->name();
	my $kyoto_filter_patients =  $project_end->getCachePatientKyotoFile($patient->name);# $dir_store."/".$patient->name().".kct";
	system("rm -f $kyoto_filter_patients ") if -e $kyoto_filter_patients;
	 my $db = new KyotoCabinet::DB;
	if (!$db->open($kyoto_filter_patients, $db->OWRITER|$db->ONOLOCK | $db->OCREATE|$db->OTRUNCATE )) {
     	printf STDERR ("open error: %s\n", $db->error);
     	die();
 	}
 	system ("chmod a+w $kyoto_filter_patients");
	foreach my $chr (@{$project_end->getChromosomes()}) {
			my $kyoto_filter_patients =   $dir_store."/".$chr->name.".".$patient->name().".kct";
			my $kh;
			die();
#			tie(%{$kh}, 'KyotoCabinet::DB', $kyoto_filter_patients ,KyotoCabinet::DB::ONOLOCK | KyotoCabinet::DB::OREADER ) || die($kyoto_filter_patients);
			$db->set_bulk($kh);
			untie %$kh;
			unlink $kyoto_filter_patients;
		}
		$db->close();;
	}
		
		
		
	warn "###################################################\n";
	warn "Purging store file !!!\n";
	warn "###################################################\n";	
	
	
	
	 
	 foreach my $chr (@{$project_end->getChromosomes()}) {	
		my $file_store = $dir_store."/".$chr->name.".".$project_end->getDataBase.".store";
		unlink($file_store);
		my $file_cache_variation = $project_end->getCacheVariationsFile($chr->name());
		if (-e $file_cache_variation){
			unlink($file_cache_variation);
		}
		#my $store_variations_file = $dir_store."/".$chr->name.".".$project_end->getDataBase.".var.store";
	
		#unlink($store_variations_file);
	 }

	warn "\n\n###################################################\n";
	warn "Done $project_name\n";
	warn "###################################################\n";	
	#return;

#
#	  foreach my $p (@{$project_end->getPatients}){
#	  	my $pname = $p->name();
#	 			 my $output   =  $project_end->getCoverageDir."/$pname.tmp.kct";
#	 
#	 	my $db1 = new KyotoCabinet::DB;
#	 
#		if (!$db1->open($output, $db1->ONOLOCK | $db1->OCREATE | $db1->OREADER )){
#				printf STDERR ("open error: %s   => %s\n", $output,$db1->error);
#				confess();
#		}
#	 	
#	 		my $cur = $db1->cursor;
# 			$cur->jump;
# 			while (my ($key, $value) = $cur->get(1)) {
# 					$p->transcriptsCoverage->set($key,$value);
# 			}
# 			unlink $output;
#	}
#	 
#
#
##my $pm3 = new Parallel::ForkManager(1);  
# 
# foreach my $p (@{$project_end->getPatients}){
#	 	warn "close ".$p->name();
#		$p->transcriptsCoverage->close();
#	}
#	


return;
}




sub cache_polydiag {
	my ($project_name,$fork) = @_;
	my $buffer = new GBuffer;
	my $project = $buffer->newProject( -name => $project_name, -verbose =>1 );
	

	my $tbundle;
	foreach my $t (@{$project->bundle_transcripts()}){
		$tbundle->{$t} ++;
	}
	warn "end 1";
	 my $nbt = scalar (@{$project->getPatients});
	 my $pr = String::ProgressBar->new( max => $nbt );
  	my $nbb =1;
  	warn "###################################################\n";
	warn "Prepare Cache for  Polydiag !!!\n";
	warn "###################################################\n";
		$| =1;
	 
		my $pm = new Parallel::ForkManager($fork);
		my $error;
	 $pm->run_on_finish( sub {
       my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
       $pr->update($nbb++);
    	 	$pr->write();
      $error =1 unless defined $data_structure_reference ;
      die("fdjkdfjk") unless defined $data_structure_reference ;
  });
	
	my $project_name = 	$project->name();
	
	$pr->write();
	foreach my $p (@{$project->getPatients}){
			my $pname = $p->name();
			#next if $pname ne "AS1502552";
			my $pid = $pm->start and next;
			my $resp = run_cache_polydiag($project_name,$p->name,$tbundle);
			$resp = $resp+0;
			$pm->finish(0,\$resp);	
	}	
	

	$pm->wait_all_children;
	$project->buffer->dbh->disconnect();
	die() if $error;

}
#sub cache_dejavu {
#	my ($project_name,$fork) = @_;
#
#	#my $root_dir = "/data-isilon/dejavu/projects/";
#	warn "*** CAche For Deja Vu *****";
#	my $buffer1 = new GBuffer;
#	my $project = $buffer1->newProject( -name => $project_name, -verbose =>1 );
#	my $root_dir =$project->deja_vu_lite_dir()."/projects/";
#
#		mkdir $root_dir unless -e $root_dir;
#	unlink $root_dir."/".$project_name.".lite" if -e $root_dir."/".$project_name.".lite";
#	my @chr_names = map{$_->name} @{$project->getChromosomes};
#	my $pm = new Parallel::ForkManager($fork);
#	my @patient_names =  sort{$a cmp $b} map{$_->name} @{$project->getPatients};
#	my $hpatients;
#	for (my $i=0;$i<@patient_names;$i++){
#		$hpatients->{$patient_names[$i]} = $i;
#	}
#	my $no = GenBoNoSql->new(dir=>$root_dir,mode=>"w");
#	$no->put($project_name,"patients",$hpatients);
#	
#	#unlink $root_dir."$project_name.lite";
#	$no->close;
#	$no = undef;
#	$pm->run_on_finish(
#    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
#  
#      	my $no = GenBoNoSql->new(dir=>$root_dir,mode=>"w");
#      	my $chr_name = $data->{chromosome};
#      	delete $data->{chromosome};
#      	delete $data->{project};
#      	
#      	$no->put($project_name,$chr_name,$data->{data});
#      	$no->close;
#      	$no = undef;
#    }
#    );
#	
#	foreach my $chr_name (@chr_names) {
#		#warn $chr_name;
#		 #next if $chr_name ne "X";
#		
#			my $pid = $pm->start and next;
#			my $buf = new GBuffer;
#			my $pr = 	$buf->newProject( -name => $project_name, -verbose =>1 );
#			my $chr = $pr->getChromosome($chr_name);
#			my $references = $chr->getReferences(undef,undef,10);
#			my $results2;
#			foreach my $r (@$references){
#			my $vs = $r->getStructuralVariations;
#			foreach my $v (@$vs){
#				my $key = $v->id;
#				
#				my $debug;
#			#	$debug =1 if $v->id eq "X_153696898_C_T";
#				my $aho= [];
#				my $ap=[];
#				foreach my $p (@{$v->getPatients()}){
#				#	warn $p->name if $debug;
#					my $pn = $p->name();
#					my $patient_id = $hpatients->{$pn};
#					push(@$ap,$patient_id);
#					
#					#$results2->{data}->{$key}->{patients}->{$pn} ++; 
#					push(@$aho,$patient_id) if ($v->isHomozygote($p));
#			}
#			$results2->{data}->{$key} = join(",",sort{$a <=> $b} @$ap);
#			if (scalar(@$aho)) {
#				$results2->{data}->{$key} = $results2->{data}->{$key}." ".join(",",sort{$a <=> $b} @$aho)." HO";
#			}
#			}
#				$pr->purge_memory_reference($r->id);
#			}
#			
#			$results2->{chromosome} = $chr_name;
#			$results2->{project} = $project_name;
#			# warn "$chr_name \n";
#			$pm->finish(0,$results2);	
#	}	
#		$pm->wait_all_children;
#	
#	$project->buffer->dbh->disconnect();
#		return 1;
#}

sub run_cache_polydiag1 {
	my ($project_name,$patient_name,$tbundle) = @_;
	my $buffer1 = new GBuffer;
	my $project = $buffer1->newProject( -name => $project_name, -verbose =>1 );

	my $vquery = validationQuery->new(dbh=>$buffer1->dbh,capture_name=>$project->validation_db());
	my $p = $project->getPatient($patient_name);

			my $db_lite = $project->noSqlPolydiag("c");
			
				
		#	my $db = $project->buffer->open_kyoto_db($file_out,'c');
			my $vtr;
			my $variations = $p->getStructuralVariations();
			my %th;
		#	warn scalar(@$variations);
			#	die();
			my $ii =0;
			my $dd =0;
		foreach my $v (@{$variations}){
				my $debug;
			#	$dd++;
				#warn $dd."/".scalar(@$variations);
				#$debug = 1 if $v->id eq "14_93670213_A_AT";
				#die() if $debug;
				#warn $ii++;
				
			my $ok;
			my $transcripts = $v->getTranscripts;
			foreach my $tr (@{$transcripts}){
		
			#	next unless exists $tbundle->{$tr->name};
					$ok =1;
					my $h = CacheDiag::construct_variant($project,$v,$tr,$p,$vquery);
					CacheDiag::update_deja_vu($project,$tr,$h); 
					my $id = join(";",$tr->id,$v->id);
					$db_lite->put($patient_name,$id,$h);
			#		$db->set($id,freeze $h );
					push(@{$vtr->{$tr->id}},$v->id);
					$th{$tr->id}++;
			}
			
			if ($v->isIntergenic){
				my $h = CacheDiag::construct_intergenic_variant($project,$v,$p,$vquery);
					CacheDiag::update_deja_vu($project,$v,$h); 
					my $id = join(";","intergenic",$v->id);
					$db_lite->put($patient_name,$id,$h);
				#	$db->set($id,freeze $h );
					push(@{$vtr->{intergenic}},$v->id);
					$th{"intergenic"}++;
				
			}

			
		}
		foreach my $t (keys %$vtr){
			$db_lite->put($patient_name,"list_$t",join(";",@{$vtr->{$t}}));
			 #$db->set("list_$t",join(";",@{$vtr->{$t}}));
		}
		$db_lite->put($patient_name,"transcripts",join(";",keys %th));
		 #$db->set("transcripts",join(";",keys %th));
		my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
		my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);

		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
		my $date = $mday."/".($mon+1)."/".(1900+$year);
 		#$db->set("date",$date);
 		$db_lite->put($patient_name,"date",$date);
 		my $f1s = $p->getVariationsFiles();
 		my $tf;
 		foreach my $f1 (@$f1s){
 		if (-e $f1){
 			$tf->{$f1} = file_md5_hex($f1);
 			#push(@$tf,file_md5_hex($f1));
 		
 			#$db->set("variations_vcf_md5",file_md5_hex($f1)) ;
 		}
 		}
 			$db_lite->put($patient_name,"variations_vcf_md5",$tf) ;
 			$tf ={};
 		my $f2s = $p->getIndelsFiles();
 		foreach my $f1 (@$f2s){
 			$tf->{$f1} = file_md5_hex($f1);
# 			push(@$tf,file_md5_hex($f1)) if -e $f1;;
 	 		}
 	 		
 	 	$db_lite->put($patient_name,"indels_vcf_md5",$tf);
 	 		
 		#$db->set("indels_vcf_md5",file_md5_hex($f2)) if -e $f2;
		#$db->close();
		#warn "end run cache";
		$db_lite->close();
		$project->buffer->dbh->disconnect();
			if (exists $project->{cosmic_db}){
			$project->{cosmic_db}->close();
		}
		$project = undef;
		#$project->{cosmic_db}->close();
		
		return 1;
}
sub run_cache_polydiag {
	my ($project_name,$patient_name,$tbundle) = @_;
	my $buffer1 = new GBuffer;
	my $project = $buffer1->newProject( -name => $project_name, -verbose =>1 );

	my $vquery = validationQuery->new(dbh=>$buffer1->dbh,capture_name=>$project->validation_db());
	my $pa = $project->get_only_list_patients($patient_name);
	my $p = $pa->[0];
		my $db_lite = $project->noSqlPolydiag("c");
		my $vtr;
		my %th;
		foreach my $chr (@{$project->getChromosomes} ){
				warn $chr->name();
		#	my $db = $project->buffer->open_kyoto_db($file_out,'c');
		
			my $variations = $chr->getStructuralVariations();
			warn "end get variations ".$chr->name;
			
		#	warn scalar(@$variations);
			#	die();
			my $ii =0;
			my $dd =0;
		foreach my $v (@{$variations}){
				my $debug;
			#	$dd++;
				#warn $dd."/".scalar(@$variations);
				#$debug = 1 if $v->id eq "14_93670213_A_AT";
				#die() if $debug;
				#warn $ii++;
				
			my $ok;
			my $transcripts = $v->getTranscripts;
			foreach my $tr (@{$transcripts}){
		
				next unless exists $tbundle->{$tr->name};
					$ok =1;
					my $h = CacheDiag::construct_variant($project,$v,$tr,$p,$vquery);
					CacheDiag::update_deja_vu($project,$tr,$h); 
					my $id = join(";",$tr->id,$v->id);
					$db_lite->put($patient_name,$id,$h);
			#		$db->set($id,freeze $h );
					push(@{$vtr->{$tr->id}},$v->id);
					$th{$tr->id}++;
			}
			
			if ($v->isIntergenic){
				my $h = CacheDiag::construct_intergenic_variant($project,$v,$p,$vquery);
					CacheDiag::update_deja_vu($project,$v,$h); 
					my $id = join(";","intergenic",$v->id);
					$db_lite->put($patient_name,$id,$h);
				#	$db->set($id,freeze $h );
					push(@{$vtr->{intergenic}},$v->id);
					$th{"intergenic"}++;
				
			}
		}
		$project->purge_memory($chr->length);
		}#end chromosome
		
		foreach my $t (keys %$vtr){
			$db_lite->put($patient_name,"list_$t",join(";",@{$vtr->{$t}}));
			 #$db->set("list_$t",join(";",@{$vtr->{$t}}));
		}
		$db_lite->put($patient_name,"transcripts",join(";",keys %th));
		 #$db->set("transcripts",join(";",keys %th));
		my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
		my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);

		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
		my $date = $mday."/".($mon+1)."/".(1900+$year);
 		#$db->set("date",$date);
 		$db_lite->put($patient_name,"date",$date);
 		my $f1s = $p->getVariationsFiles();
 		my $tf;
 		foreach my $f1 (@$f1s){
 		if (-e $f1){
 			$tf->{$f1} = file_md5_hex($f1);
 			#push(@$tf,file_md5_hex($f1));
 		
 			#$db->set("variations_vcf_md5",file_md5_hex($f1)) ;
 		}
 		}
 			$db_lite->put($patient_name,"variations_vcf_md5",$tf) ;
 			$tf ={};
 		my $f2s = $p->getIndelsFiles();
 		foreach my $f1 (@$f2s){
 			$tf->{$f1} = file_md5_hex($f1);
# 			push(@$tf,file_md5_hex($f1)) if -e $f1;;
 	 		}
 	 		
 	 	$db_lite->put($patient_name,"indels_vcf_md5",$tf);
 	 		
 		#$db->set("indels_vcf_md5",file_md5_hex($f2)) if -e $f2;
		#$db->close();
		#warn "end run cache";
		$db_lite->close();
		$project->buffer->dbh->disconnect();
			if (exists $project->{cosmic_db}){
			$project->{cosmic_db}->close();
		}
		$project = undef;
		#$project->{cosmic_db}->close();
		
		return 1;
}
sub cache_cnv {
	my ($project_name,$fork) = @_;
	my $buffer1 = new GBuffer;
	my $projectP = $buffer1->newProject( -name => $project_name, -verbose =>1 );
	my @transcripts_cgi = @{$projectP->bundle_transcripts() } ;
	my @transcripts = sort{$a->getGene->external_name cmp $b->getGene->external_name} map{$projectP->newTranscript($_)} @transcripts_cgi ;
	warn "###################################################\n";
	warn "Prepare Cache for  CNV  !!!\n";
	warn "###################################################\n";
	my $primers = $projectP->getPrimers();
	foreach my $patient (@{$projectP->getPatients}){
	foreach my $primer (@$primers){
		$primer->cached_cnv($patient);
		#print ".";
	}
	}
	#print "\n";
	preload_coverage::load_cnv_score($projectP,$projectP->getPatients,\@transcripts);
	
	warn "###################################################\n";
	warn "END CNV  !!!\n";
	warn "###################################################\n";
	warn "\n";
}

sub cache_coverage_primers {
	my ($project_name,$fork) = @_;
	my $buffer1 = new GBuffer;
	my $projectP = $buffer1->newProject( -name => $project_name, -verbose =>1 );
	my $pm2 = new Parallel::ForkManager($fork);  
	
	my $no = $projectP->noSqlCoverage();
	my $total ;
	
	$pm2->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
   		$total->{$data->{chr}} = $data->{primers};
    }
  );
  
  foreach my $chr (@{$projectP->getChromosomes}){
  		my $chr_name = $chr->name();
  		my $pid = $pm2->start() and next;
			my $buffer = new GBuffer;
			my $project = $buffer->newProject( -name => $project_name, -verbose =>1 );
			my $chr1 = $project->getChromosome($chr_name);
			
			my $primers = $chr1->getPrimers();
			map {delete $_->{project}; delete $_->{buffer}} @$primers;
			my $res;
			$res->{primers} = $primers;
			$res->{chr} = $chr->name;
			$pm2->finish(0,$res);	
	}
	$pm2->wait_all_children;	
	$no->put_bulk("primers",$total);
  
	$no->close();
}


#list primers compute a random list of primers/capture => if you have to much primers/capture
 
 sub cache_coverage_list_primers {
	my ($project_name,$fork) = @_;
	my $buffer1 = new GBuffer;
	my $projectP = $buffer1->newProject( -name => $project_name, -verbose =>1 );
	my $no = $projectP->noSqlCoverage();
	foreach my $capture (@{$projectP->getCaptures}){
	next if $capture->isPcr();	
		my $list = $capture->getListPrimers();
		$no->put($projectP->name(),$capture->id,$list);
	}
	$no->close();
 }
 sub cache_delete_coverage_file {
	my ($project_name,$fork) = @_;
	my $buffer1 = new GBuffer;
	my $projectP = $buffer1->newProject( -name => $project_name, -verbose =>1 );
	my $no = $projectP->noSqlCoverage();
	my $output   =$projectP->getCacheDir() . "/coverage_lite";
unlink $output."/".$projectP->name().".lite";
unlink $output."/primers.lite" if -e $output."/primers.lite" ;

my $no = $projectP->noSqlCoverage();
# #my $no = $projectP->noSqlCoverage();
# warn "clean";
$no->put($projectP->name(),"test","tutu");
$no->close();
 }
 
sub cache_coverage {
	my ($project_name,$fork) = @_;
	warn "###################################################\n";
	warn "Prepare Cache for  COVERAGE !!!\n";
	warn "###################################################\n";
	
	cache_delete_coverage_file ($project_name,$fork);
	cache_coverage_primers ($project_name,$fork);
	cache_coverage_list_primers($project_name,$fork);
	
	#############
	# 
	#############
	
	my $buffer1 = new GBuffer;
	my $projectP = $buffer1->newProject( -name => $project_name, -verbose =>1 );
	my $pm2 = new Parallel::ForkManager($fork);  
	warn $fork;
	my @patients_name = map{$_->name()} @{$projectP->getPatients};
	my $nbt = scalar(@patients_name);
  	my $pr = String::ProgressBar->new( max => $nbt );
  	my $nbb =1;
  	$pm2->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    	 $pr->update($nbb++);
    	 	$pr->write();
 	
    }
  );
  my $diag;
  $diag=1 if $projectP->isDiagnostic();
  
# unless ($diag){
	

#$no->clear($projectP->name());
 unless ($diag){
  	 foreach my $patient_name (@patients_name){
  	 
  		my $pid = $pm2->start() and next;
  				compute_coverage_exome($project_name,$patient_name);
  		
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
foreach my $patient_name (@patients_name){
  	 #	last if $nb == 1;
  	 #	$nb++;
  	 	
  		my $pid = $pm2->start() and next;
  		
  		compute_coverage_diagnostic1($project_name,$patient_name);
  		$pm2->finish(0);	
 	#	last; 		
	}

	$pm2->wait_all_children;
	warn "now Primers *******************";

	 $pr = String::ProgressBar->new( max => $nbt );
	 $nbb = 0;
	$pr->update(0);
		$pr->write();
  	 
  	 foreach my $patient_name (@patients_name){
  	 
  		my $pid = $pm2->start() and next;
  		compute_coverage_diagnostic2($project_name,$patient_name);
  		
  			$pm2->finish(0);	
	}

	$pm2->wait_all_children;
	#$fork =1;
	#$pm2 = new Parallel::ForkManager($fork);  
	warn "and finally ..  exons *****************************";
	
		
	 $pr = String::ProgressBar->new( max => $nbt );
	 warn $nbt;
	 $nbb = 0;
	$pr->update(0);
		$pr->write();
  	 my @utrs = (1);
  	 #	my $pm3 = new Parallel::ForkManager($fork);  
  	 foreach my $patient_name (@patients_name){
  	# 	last if $nb == 2;
  	 	$nb++;
  		my $pid = $pm2->start() and next;
  		#foreach my $utr (@utrs){
  			
  			compute_coverage_diagnostic3($project_name,$patient_name,0);
  			
		#}
	$pm2->finish(0);	

  	}
	$pm2->wait_all_children;


	print "\n";
	
#	preload_coverage::load_cnv_score($projectP,$project->getPatients,\@transcripts,1);
}
sub compute_coverage_diagnostic1{
	my ($project_name,$patient_name) = @_;
		my $buffer = new GBuffer;
	my $project = $buffer->newProject( -name => $project_name, -verbose =>1 );
	my $patient = $project->getPatient($patient_name);
	my $no =  $project->noSqlCoverage();
	$no->clear($patient->name);
	$no->clear($patient->name."_cnv");
	my @trs = map{$project->newTranscript($_)} @{$project->bundle_transcripts()} ;
	preload_coverage::load_coverage_transcripts($project,[$patient],\@trs);
	$no->close();
	undef($no);

}
sub compute_coverage_diagnostic2{
	my ($project_name,$patient_name) = @_;
		my $buffer = new GBuffer;
	my $project = $buffer->newProject( -name => $project_name, -verbose =>1 );
	my $patient = $project->getPatient($patient_name);
	my $no =  $project->noSqlCoverage();
	
	my @trs = map{$project->newTranscript($_)} @{$project->bundle_transcripts()} ;
	preload_coverage::load_coverage_primers($project,[$patient],\@trs) ;
	$no->close();
	undef($no);
		return;
}
sub compute_coverage_diagnostic3{
	my ($project_name,$patient_name,$utr) = @_;
		my $buffer = new GBuffer;
	my $project = $buffer->newProject( -name => $project_name, -verbose =>1 );
	my $patient = $project->getPatient($patient_name);
	my $no =  $project->noSqlCoverage();
	my @paddings= (0,5,10,15,20,30);
	my @utrs = (0,1);
	my @trs = map{$project->newTranscript($_)} @{$project->bundle_transcripts()} ;
	warn scalar(@trs);
	my $utr = 0;
	my $padding =0;
	#	foreach my $utr (@utrs){
					#foreach my $padding (@paddings){
						warn "preload";
						preload_coverage::load_coverage_for_cache1($project,[$patient],\@trs ,$padding,$utr);
						
					
					#}
			#	preload_coverage::load_coverage($project,[$patient],\@trs ,20,1);	
	#	}
		$no->close();
		undef($no);
		return;
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

sub cache_exons_coverage {
	my ($project_name,$fork) = @_;
	my $buffer1 = new GBuffer;
	my $projectP = $buffer1->newProject( -name => $project_name, -verbose =>1 );
	my @paddings= (0,5,10,15,20,30);
	my @utrs = (0,1);
	my $pm2 = new Parallel::ForkManager($fork);  

	my $nb =0;
		my @transcripts_cgi = @{$projectP->bundle_transcripts() } ;
	#@transcripts_cgi = splice(@transcripts_cgi,0,24);
	my $nbt = scalar(@transcripts_cgi);
	
  my @patients_name = map{$_->name()} @{$projectP->getPatients};
  my $nbt = scalar(@patients_name);
  	my $pr = String::ProgressBar->new( max => $nbt );
  	my $c =0;
	$pm2->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    	 $pr->update($c++);
 	
    }
  );
  foreach my $patient_name (@patients_name){
  		$pr->write();
  	my $pid = $pm2->start() and next;
		my $buffer_fork = new GBuffer;
		my $project_fork = $buffer_fork->newProject( -name => $project_name, -verbose =>1 );
		my $patient = $project_fork->getPatient($patient_name);
		my $no =  $project_fork->noSqlCoverage();
		foreach my $transcript_name (@transcripts_cgi){
			my $transcript =  $project_fork->newTranscript($transcript_name);
			my $exons  = $transcript->getAllGenomicsParts();
			foreach my $exon  (@$exons) {
				foreach my $utr (@utrs){
					foreach my $padding (@paddings){
						 my ($mean,$intspan,$min)  = $exon->statistic_coverage_coding(patient=>$patient,padding=>$padding,limit=>1,utr=>$utr);
						 $no->put($patient_name,$exon->id."_".$padding."_".$utr,{mean=>$mean,min=>$min});
						 
					 	#warn $mean." ".$min;
					 
					}
				}
			}
		}
		my $resp;
	 	$pm2->finish(0,$resp);	
  }
		$pm2->wait_all_children;
				$pr->write();
}


sub cache_images{
	my ($project_name,$fork) = @_;
		warn "\n\n###################################################\n";
	warn "IMAGES : COVERAGE \n";
	warn "###################################################\n";		
		my $buffer1 = new GBuffer;
	my $projectP = $buffer1->newProject( -name => $project_name, -verbose =>1 );
	
	my @paddings= (0,5,10,15,20,30);
	my @limits= (0,5,10,15,20,30);
	my @utrs = (0,1);
	my @intronics = (0);
	my @params;
	my @transcripts_cgi = @{$projectP->bundle_transcripts() } ;
	
  
	#foreach my $tr (@transcripts_cgi){
		foreach my $intronic (@intronics){
			foreach my $utr (@utrs) {
				foreach my $limit (@limits) {
					foreach my $padding (@paddings){
						my %t;
				#		$t{transcript} = $tr;
						$t{utr} = $utr;
						$t{limit} = $limit;
						$t{padding} = $padding;
						$t{intronic} = $intronic;
						push(@params,\%t);
					}
				}
			}
		}
#	}
	

	my $patients = $projectP->getPatients();
	

	
	$| =1;
	
	my $pm2 = new Parallel::ForkManager($fork);  
	my $nb =0;
	#@transcripts_cgi = splice(@transcripts_cgi,0,24);
	my $nbt = scalar(@transcripts_cgi);
	#my $pr = String::ProgressBar->new( max => scalar(@params) );
	my $pr = String::ProgressBar->new( max => $nbt );
	my $c =0;
	 my %pids;	
	
	my $total;
	
	$pm2->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
     $pr->update($c++);
  	foreach my $k (keys %$data){
    				$total->{$k} = $data->{$k};
    			}

 	$pr->write();
    }
  );
  	$pr->write();
foreach my $t (@transcripts_cgi){
	my $pid = $pm2->start() and next;
	my $resp = process_cache_images($project_name,\@params,$t);
	 $pm2->finish(0,$resp);	
}

#
#foreach my $p (@params){
#	my $pid = $pm2->start() and next;
#	
#	my $resp = process_cache_images($project_name,$p);
#	
#  $pm2->finish();	
#}

	
	$pm2->wait_all_children;
	
 $pr->write();
 warn "###################################################\n";
warn "END  IMAGES : COVERAGE !!!\n";
warn "###################################################\n";		
warn "\n";
 warn "###################################################\n";
warn "START  IMAGES : CNV !!!\n";
warn "###################################################\n";	
my $pm3 = new Parallel::ForkManager($fork);  
	my $nb =0;
	#@transcripts_cgi = splice(@transcripts_cgi,0,24);
	my $nbt = scalar(@transcripts_cgi);
	my $pr2 = String::ProgressBar->new( max => $nbt);
	 $c =0;
	 my %pids;	
	
	
	
	$pm3->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
     $pr2->update($c++);
  	foreach my $k (keys %$data){
    				$total->{$k} = $data->{$k};
    			}

 	$pr2->write();
    }
  );


foreach my $t (@transcripts_cgi){
	my $pid = $pm3->start() and next;
	
	my $resp = process_cache_images_cnv($project_name,$t);
	
  $pm3->finish(0,$resp);	
}

	
	$pm3->wait_all_children;

warn "###################################################\n";
warn "END  IMAGES : CNV !!!\n";
warn "###################################################\n";	
 			 			
 	
#	 my $db1 = new KyotoCabinet::DB;
  #my $db1 = tie(my %db, 'KyotoCabinet::DB', '/scratch/toto.kch');

 # open the database
 

 my $file_project = $projectP->coverageImagesCoverageFile;

 
 my $db_out = $projectP->buffer->open_kyoto_db($file_project,'c');
  foreach my $k (keys %$total){
  	$db_out->set($k,$total->{$k});
  }
 

		$db_out->close();






warn "###################################################\n";
warn "END IMAGES !!!\n";
warn "###################################################\n";	
warn $projectP->coverageImagesCoverageFile;	
	
	
}

sub process_cache_images_cnv{
	my ($project_name,$tr) = @_;
	
	my $buffer1 = new GBuffer;
	my $projectP = $buffer1->newProject( -name => $project_name, -verbose =>1 );
	my $patients = $projectP->getPatients();
	my $tr1 = $projectP->newTranscript($tr);
	
	my $ret  = image_coverage::image_cnv ($patients, $tr1 );
	
	my $kyoto_id = join("_",("cnv",$tr));
	my %res;
	$res{$kyoto_id}	= $ret->{image}->png;
	$res{$kyoto_id."_data"} = freeze $ret->{data};

	return \%res;
	
}

sub process_cache_images {
	my ($project_name,$params,$tr) = @_;
	
	my $buffer1 = new GBuffer;
	my $projectP = $buffer1->newProject( -name => $project_name, -verbose =>1 );
	my $patients = $projectP->getPatients();
	
	
	
	
	#my $tr =  $p->{transcript};
	my $tr1 = $projectP->newTranscript($tr);
			
	
	my %res;
	foreach my $p (@$params){
		my $utr = $p->{utr};
		my $intronic = $p->{intronic};
		my $limit = $p->{limit};
		my $padding = $p->{limit};
		my $kyoto_id = join("_",("all",$tr,$utr,$intronic,$limit,$padding));
		
		my $ret = image_coverage::image ($patients, $tr1,$intronic,$utr, $padding, $limit,1);
		$res{$kyoto_id} = $ret->{image}->png;
		$res{$kyoto_id."_data"} = freeze $ret->{data};
	}
		
	return \%res;
}



sub start_primer {
	my ($project_name,$chr_name) = @_;
	my $buffer1 = new GBuffer;
	my $projectP = $buffer1->newProject( -name => $project_name, -verbose =>1 );
	my $chr = $projectP->getChromosome($chr_name);
	my $primers = $chr->getPrimers();
	warn "chr".$chr_name." : ".scalar(@$primers);
	my $res;
	my $nb = scalar(@$primers);
	my $z =0;
	foreach my $primer (@$primers){
		warn "$z/$nb";
		$z++;
		foreach my $p (@{$projectP->getPatients}){
			warn $p->name();
			my $id = $primer->id."_".$p->id;
			my $cnv = $primer->compute_cnv_score($p);
			$res->{$p->name}->{$primer->id}->{cnv} = $cnv;
			$res->{$p->name}->{$primer->id}->{level_cnv} = $primer->level($p);
			my $min = $primer->minimum($p);
			warn $min;
			my $mean = $primer->mean($p);
			$res->{$p->name}->{$primer->id}->{mean} = $mean;
			$res->{$p->name}->{$primer->id}->{min} = $min;
		}
		
	}
	return $res;	
}




sub start_thread {
	my ($project_name,$chr_name) = @_;
	die();
	my $buffer1 = new GBuffer;
	my $data;
	my $projectP = $buffer1->newProject( -name => $project_name, -verbose =>1 );

	my $chr = $projectP->getChromosome($chr_name);
	my $dir_store = $projectP->getRootDir()."/";

	
	my $data_variations;
	my $store_filename = $dir_store."/".$chr->name.".".$projectP->getDataBase.".store";
	my $store_variations_file = $dir_store."/".$chr->name.".".$projectP->getDataBase.".var.store";
	my $file_variations =  $dir_store."/".$chr->name.".".$projectP->getDataBase.".var.kct";
	my $kyoto_gene_file =   $dir_store."/".$chr->name.".genes.kct";
	my $kyoto_filter_genes =   $dir_store."/".$chr->name.".genes.filter.kct";
	
	unlink ($store_filename)  if -e $store_filename;
	unlink ($file_variations) if -e $file_variations;
	unlink ($kyoto_gene_file) if -e $kyoto_gene_file;
	
	unless (-e $store_filename && -e $file_variations){
		warn "\t\t Running :-- ". $chr->name." --\n";
	
		#$data->{$chr->name} = @{$chr->getVariations};
		$data->{name} = $chr->name;
		my $hvars;
#		tie(%{$hvars}, 'KyotoCabinet::DB', $file_variations , KyotoCabinet::DB::ONOLOCK | KyotoCabinet::DB::OWRITER | KyotoCabinet::DB::OCREATE| KyotoCabinet::DB::OTRUNCATE) || die($file_variations);
		my $hgenes;
#		tie(%{$hgenes}, 'KyotoCabinet::DB', $kyoto_gene_file , KyotoCabinet::DB::ONOLOCK | KyotoCabinet::DB::OWRITER | KyotoCabinet::DB::OCREATE| KyotoCabinet::DB::OTRUNCATE) || die($kyoto_gene_file);
		my $fgenes;
#		tie(%{$fgenes}, 'KyotoCabinet::DB', $kyoto_filter_genes , KyotoCabinet::DB::ONOLOCK | KyotoCabinet::DB::OWRITER | KyotoCabinet::DB::OCREATE| KyotoCabinet::DB::OTRUNCATE) || die($kyoto_filter_genes);
#	
		my $hpatients;
		foreach my $patient (@{$projectP->getPatients()}){
			my $kyoto_filter_patients =   $dir_store."/".$chr->name.".".$patient->name().".kct";
#			tie(%{$hpatients->{$patient->name}}, 'KyotoCabinet::DB', $kyoto_filter_patients , KyotoCabinet::DB::ONOLOCK | KyotoCabinet::DB::OWRITER | KyotoCabinet::DB::OCREATE| KyotoCabinet::DB::OTRUNCATE) || die($kyoto_filter_patients);
		}
	
		
		($data->{data},$data_variations) = CacheGenesData_nodb::create_cache_genes($projectP,$chr,$hvars,$hgenes,$fgenes,$hpatients);
			
		$hvars->{type} = "nodb";
		foreach my $patient (@{$projectP->getPatients()}){
		untie  $hpatients->{$patient->name};
		#	tie(%{$hpatients->{$patient->name}}, 'KyotoCabinet::DB', $kyoto_filter_patients , KyotoCabinet::DB::ONOLOCK | KyotoCabinet::DB::OWRITER | KyotoCabinet::DB::OCREATE| KyotoCabinet::DB::OTRUNCATE) || die($kyoto_filter_patients);
		}
		
		store($data,$store_filename);
		#store($data_variations,$store_variations_file);
		$data = undef;
	} 
	else {
		warn "======== GENE : nothing to do for ".$chr->name ."$store_filename use cache ======\n"; 
		
	}
	
	warn "==== End GENE CACHE FOR ".$chr->name." =====\n";
	warn "$kyoto_filter_genes\n";
	return 1;
	
	
}  


sub exons_coverage{
	my ($trs,$dir,$project_name) =@_;
	my @paddings= (0,5,10,15,20,30);
	my $buffer2 = GBuffer->new();
	my $project2 = $buffer2->newProject(-name=>$project_name);
	my $patients = $project2->getPatients();
	
		
	foreach my $tr (@$trs){
		warn $tr;
		my $transcript = $project2->newTranscript($tr);
		my $capture_intspan = $transcript->getChromosome->getIntSpanCapture();
		my $htr;
			$htr->{name} = $transcript->name();
			$htr->{start} = $transcript->start();	
			$htr->{end} = $transcript->end();	
			$htr->{strand} = $transcript->strand();	
			$htr->{chromosome} = $transcript->getChromosome()->name;
		
		foreach my $pad (@paddings){
				$htr->{nb_exons} =0;
				$htr->{nb_introns} =0;
			
	 	foreach my $exon (@{$transcript->getAllGenomicsParts}){
	 		my $hexon;
	 		$hexon->{name} = $exon->name;
	 		$hexon->{id} = $exon->id;
	 		$hexon->{start} = $exon->start;
	 		$hexon->{end} = $exon->end;
	 		$hexon->{exon} = $exon->isExon();
	 		$hexon->{coding}  = 1;
	 		$hexon->{coding}  = undef if $exon->is_noncoding() ;
	 		if ($exon->isExon()){
	 				 $htr->{nb_exons}  ++;
	 				 $hexon->{exon} = 1;
	 		}
	 		else {
	 			$htr->{nb_introns}  ++;
	 			 $hexon->{exon} = undef;
	 		}

	 		my $s1 = $exon->getGenomicSpan()->intersection($capture_intspan);
	 		my $covered = 1;
			if ($exon->isExon()){
				$covered =undef if $s1->is_empty;
			}
			else {
				$covered =undef if scalar($s1->as_array)<20;
			}
			$hexon->{covered} = $covered;
			
		
	 		foreach my $patient (@$patients) {
	 				my $pid = $patient->id;
	 				$hexon->{$pid}->{stat}->{mean} = -1;
	 				$hexon->{$pid}->{stat}->{min} = -1;
	 				$hexon->{$pid}->{stat}->{coding_stats} = -1;
	 				$hexon->{$pid}->{stat}->{coding_stats} = -1;
	 				if ($covered){
	 					my @z = $exon->mean_intspan_coverage(patient=>$patient,padding=>$pad,limit=>5);
	 					$hexon->{$pid}->{stat}->{mean} = $z[0];
	 					$hexon->{$pid}->{stat}->{min} = $z[2];
	 					my @z2 = $exon->mean_intspan_coverage_coding(patient=>$patient,padding=>$pad,limit=>5);
	 					$hexon->{$pid}->{coding_stats}->{mean} = $z2[0];
	 					$hexon->{$pid}->{coding_stats}->{min} = $z2[2];
	 				}
	 			}
	 			
	 			push(@{$htr->{exons}},$hexon);
	 		
	 		}#end exons
				my $fileout = $dir."/".$transcript->name.".stats.$pad.freeze";
				store($htr,$fileout);
				
	 		delete $htr->{exons};
	 	
		}#end padding
	}
}

 
 1;