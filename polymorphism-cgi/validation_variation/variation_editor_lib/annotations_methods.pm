package annotations_methods;
use strict;
use Data::Dumper;
use Time::HiRes qw( time);
use Carp;

sub annotations {
	my ( $project,$patient, $tmp, $list_saved, $maskcoding,$final_polyviewer_all,$hash_genes_panel,$hash_variants_DM ) =@_;

	 if ($project->isRocks){
	 	return annotations_rocks($project,$patient, $tmp, $list_saved, $maskcoding,$final_polyviewer_all,$hash_genes_panel,$hash_variants_DM);
	 }
	 else {
	 	return annotations_lmdb($project,$patient, $tmp, $list_saved, $maskcoding,$final_polyviewer_all,$hash_genes_panel,$hash_variants_DM);
	 }
}




sub annotations_rocks {
	my ( $project,$patient, $list, $list_saved, $maskcoding,$final_polyviewer_all,$hash_genes_panel,$hash_variants_DM ) = @_;
	my $cgi = new CGI();
	my $final_polyviewer_all = GenBoNoSqlRocks->new(dir=>$project->rocks_directory."/patients/",mode=>"r",name=>$patient->name);
	
	#delete $project->{rocks};

	print "." unless ($cgi->param('export_xls'));
	my $tsum = 0;

	my $tglobal = time;
	my $res     = {};
	my $agenes  = [];
	my $e;
	my $nb = 1;
	my $javascript_id = int( time + rand(400000) );
	my $list2 = [];
	my $dchr;
	my $t = time;
	my $chr_view;
	my $current;
	my $values;
	my $tsumz = 0;
	my $tt1 = time;
	my $no;
#	warn "start prepare ".$final_polyviewer_all->dir;
	my ($all_hash) = $final_polyviewer_all->prepare($list);
	
	foreach my $id (@$list) {
		my $debug;
		$nb++;
		unless ($cgi->param('export_xls')) {
			print "." if $nb % 100 == 0;
		}

		my $hg = $final_polyviewer_all->get($id);
		if ($hash_genes_panel) {
			foreach my $ag ( @{ $hg->{array} } ) {
				next unless exists $hash_genes_panel->{ $ag->{id} };
				push( @$agenes, $ag );    # if $ag->{mask} & $maskcoding;
			}
		}
		elsif ( exists $hash_variants_DM->{$id} ) {
			foreach my $a ( @{ $hg->{array} } ) {
				$a->{DM}++;
				push( @$agenes, $a );
			}

		}
		else {
			push( @$agenes,grep { $_->{mask} & $maskcoding } @{ $hg->{array} } );
		}
	}
	$res->{genes} = $agenes;
	return $res ;
}


sub annotations_lmdb {
		my ( $project,$patient, $list, $list_saved, $maskcoding,$final_polyviewer_all,$hash_genes_panel,$hash_variants_DM ) = @_;
	my $cgi = new CGI();
	#$project->buffer->dbh_reconnect();
	print "." unless ($cgi->param('export_xls'));
	my $tsum = 0;

	#my $vs =  $project->myflushobjects($list,"variants");;
	my $tglobal = time;
	my $res     = {};
	my $agenes  = [];
	my $e;
	my $nb = 1;
	my $mtime_flush;
	my $mtime_gene;
	my $mtime;
	my $javascript_id = int( time + rand(400000) );
	my $list2 = [];
	my $dchr;
	foreach my $id (@$list) {
		my $debug;
		my ( $cname, $xid ) = split( "!", $id );
		$nb++;
		unless ($cgi->param('export_xls')) {
			print "." if $nb % 100 == 0;
		}
		my $t = time;

		$dchr->{$cname} = $project->getChromosome($cname)->lmdb_polyviewer_variants_genes( $patient, "r" )  unless exists $dchr->{$cname};

		$tsum += abs( $t - time() );

		my $hg = $dchr->{$cname}->get($id);
		unless ($hg){
			$hg = hvariant::hash_variant($patient,$id);
		}
		

		if ($hg) {
			if ($hash_genes_panel) {
				foreach my $ag ( @{ $hg->{array} } ) {
					next unless exists $hash_genes_panel->{ $ag->{id} };

					push( @$agenes, $ag );   
				}
			}
			elsif ( exists $hash_variants_DM->{$id} ) {
				foreach my $a ( @{ $hg->{array} } ) {
					$a->{DM}++;
					push( @$agenes, $a );
				}

			}
			else {
				push( @$agenes,grep { $_->{mask} & $maskcoding } @{ $hg->{array} } );
			}

			next;
		}
		push( @$list2, $id );
	}
	$res->{genes} = $agenes;
	
	return $res unless @$list2;

	###RETURN ####
	warn Dumper @$list2;
	foreach my $z (@$list2){
		my $vobj = $project->returnVariants($z);
		#warn $vobj->name;
		warn Dumper $vobj->gnomad;
	}
	confess();

}	
	
	
	
1;
