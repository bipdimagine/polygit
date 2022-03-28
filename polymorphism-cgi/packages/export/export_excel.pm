package export_excel;
use strict;
use FindBin qw($Bin);
#use lib "$Bin/layout/";
#use lib "$Bin/../layout/";
#use layout;
use Spreadsheet::WriteExcel; 
#use formatter;
use Data::Dumper;



my %formatter = (
	"formatSnp"         => \&formatter::formatSnp,
	formatMethods       => \&formatter::formatMethods,
	formatReferenceGene => \&formatter::formatReferenceGene,
	formatFilter        => \&formatter::formatFilter,
	formatValid    => \&formatter::formatValid,
	formatPolyphen => \&formatter::formatPolyphen,
	formatHapmap  => \&formatter::formatHapmap
);


sub exportTable {
	my ($array,$type_name,$project) = @_;
	
	if ( $type_name =~ /variation/ ) {
		exportFormatTable( $array, "variations",$project );
		return;
	}
	elsif ( $type_name =~ /indel/ ) {
		exportFormatTable( $array, "indels",$project );
		return;
	}
	elsif ( $type_name =~ /gene/ ) {
		exportFormatTable( $array, "genes",$project );
		return;
	}

}

sub exportFormatTable {
	my ( $data, $type_layout,$project ) = @_;

	my $layoutVariation = layout::returnLayout($type_layout,$project);
	
	
	die("problem layout") unless $layoutVariation;
		
	my $workbook        = Spreadsheet::WriteExcel->new( \*STDOUT );
	my $worksheet       = $workbook->add_worksheet();
	
	my $col = $layoutVariation->[0];

	#
	# line header
	#
	my $nbline = 0;
	my $formatTh;
	my $nbcol = 0;
	foreach my $col (@$layoutVariation) {
		next unless exists $col->{field};
		$worksheet->write_string( $nbline, $nbcol, $col->{name}, $formatTh );

		$nbcol++;
	}
	foreach my $line (@$data) {
		$nbcol = 0;
	#	next unless $line->{checkGoodBipd};
		foreach my $col (@$layoutVariation) {
			next unless exists $col->{field};
			my $funct = $col->{formatter};
			my $value = $line->{ $col->{field} };
			
			my $old_value= $value;
			my $new_value;
			if ($funct && exists $formatter{$funct}) {
				 $new_value = $formatter{$funct}->($value);					
			}
			else {
				$new_value = {type=>"value",value=>$value};
				#$worksheet->write( $nbline + 1, $nbcol, $value, $formatTh );
			}
		
			printXlsLine($worksheet,$nbline+1,$nbcol,$new_value,$formatTh);
			
			$nbcol++;
		}
		$nbline++;

	}
	exit(0);
}

sub printXlsLine {
	my ($worksheet,$line,$col,$value,$format) = @_;
	if ($value->{type} eq "value") {
		$worksheet->write( $line, $col, $value->{value}, $format );
	} 
	elsif ($value->{type} eq "url") {
		$worksheet->write( $line, $col, $value->{url},$value->{value}, $format );
	}
	
}



1;