package parse_sniffles2;
use strict;
use Data::Dumper;
use Storable qw(dclone);




sub parse_cnv {
	my ($patient,$caller) =@_;
	 my $vcf_file = $patient->getSVFile($caller);
	 my $project = $patient->project;
	# ouverture du fichier manta zippé
	my $res;
	my $vcf = Bio::DB::HTS::VCF->new( filename => "$vcf_file" );
	my $header = $vcf->header();
	
	while (my $row = $vcf->next){
		my $chr = $project->getChromosome($row->chromosome($header),1);
		
		next unless $chr;
	
		my $h;
		#champ infos 
		next unless $row->has_filter($header,".");
		my $svtype = 	$row->get_info($header, "SVTYPE");
		next unless $row->has_filter($header,"PASS");
		
		$h->{'SVTYPE'} = get_value($row->get_info($header, "SVTYPE"));
		next unless ( ($h->{'SVTYPE'} eq "DUP") || ($h->{'SVTYPE'} eq "DEL") );
	
		$h->{'CHROM'}=$chr->name;
		$h->{"CALLER"} = $caller;
		$h->{'END'} = $row->get_info($header, "END")->[0];
		$h->{'START'} = $row->position();
		$h->{'SVLEN'} =  abs($h->{'END'} - $h->{'START'});
		$h->{'KARYOTYPE_ID'}= $chr->karyotypeId;
		
		
		$h->{'CN'} ="-" ;
		$h->{'RATIO'} ="-" ;
		
		my $id = $h->{'SVTYPE'}."_".$chr->name."_".$h->{'START'}."_".$h->{'END'};
		$h->{'ELEMENTARY'}= [$id];
		$h->{'INFOS'} = $row->get_format($header);
		
		#GT:GQ:DR:DV	0/1:9:18:6
		
		$h->{'INFOS'}->{sr_qual} = "0";
		$h->{'INFOS'}->{SR} = ["0","0"];
		$h->{'INFOS'}->{PR} = ["0","0"];
		if (exists $h->{'INFOS'}->{DR} and $h->{'INFOS'}->{DV}) {
			$h->{'INFOS'}->{SR} = [$h->{'INFOS'}->{DR}->[0],$h->{'INFOS'}->{DV}->[0]];
			$h->{'INFOS'}->{PR} =  [$h->{'INFOS'}->{DR}->[0],$h->{'INFOS'}->{DV}->[0]];
		}
		$h->{'CN'} = $h->{'INFOS'}->{GT}->[0];
		$h->{'GT'} = "0/1";
		$h->{'GT'} = "1/1" if ($h->{'INFOS'}->{GT}->[0] == {'INFOS'}->{GT}->[1]);  
		
		$h->{id}=$id;
		$h->{'REAL_CALLER'} = $caller;
		my $type = $h->{'SVTYPE'};
		$h->{'QUAL'} = $row->quality;;
		$h->{'QUAL'} = "0" if $h->{'QUAL'} eq 'NaN';
		my $num = $chr->name();
		$res->{$id} = $h;
		control_object($h);
		
	}
	
return $res;
}

sub get_value {
	my($res) = @_;

	if (ref($res) eq 'ARRAY') {
			return  $res->[0];	
	}
	else {
		return undef;
	}
}
sub control_object {
	my ($h) = @_;
	my @keys = ('id','SVTYPE','CHROM','START','END','SVLEN','GT','CN','RATIO','QUAL','ELEMENTARY');
	foreach my $k (@keys){
		die($k) unless exists $h->{$k};
	}
	
}


1;
