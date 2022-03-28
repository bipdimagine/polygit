package formatter;
use strict;


sub formatSnp {
	my $value = shift;
	if ($value =~ /^rs/) {
	return 	{type=>"url", url=>"http://www.ncbi.nlm.nih.gov/snp/?term="
			  .$value,value=>$value};
	}
	 
	
			  
	return {type=>"value", value=>$value};
}

sub formatMethods {
	my $value = shift;
	return {type=>"value", value=>$value};
}

sub formatReferenceGene {
	my $value = shift;
	
	return {type=>"value", value=>join( ";", @$value )};
}

sub formatFilter {
	my $value = shift;
	return {type=>"value", value=>$value};
}

sub formatValid {
	my $value = shift;
	return {type=>"value", value=>"No validation"} if ( $value == 0 );
	return {type=>"value", value=>"Unvalidated" }  if ( $value == -1 );
#	if ( $value == -2 ) { return {type=>"value", value=>"Uncertain"} }
	if ( $value == 1 )  { return {type=>"value", value=>"Validated"} }
	if ( $value == 2 )  { return {type=>"value", value=>"Validated (ho)"} }
	if ( $value == 3 )  { return {type=>"value", value=>"Validated (h?) (he)" }}
#	if ( $value == -8 ) { return {type=>"value", value=>"Ambiguities"} }
	return {type=>"value", value=>$value};
}

sub formatPolyphen {
	my $value = shift;
	
	return {type=>"value", value=>$value};
}

sub formatHapmap {
	my $value = shift;

	if ($value) {
		return {type=>"url", url=>"http://hapmap.ncbi.nlm.nih.gov/cgi-perl/gbrowse/hapmap27_B36/?name=SNP:$value+",value=>"hapmap"};
	
	
	}
	
	return {type=>"value", value=>$value};
}

sub formatDejavu {
	my $value = shift;
	return $value;
	
	
}

sub formatHomozygote {
	my $value = shift;
	return $value;
	
	
}

sub formatValidFile {
	my $value = shift;
	return {type=>"value", value=>"file"} if ( $value == 1 ) ;
}












1;