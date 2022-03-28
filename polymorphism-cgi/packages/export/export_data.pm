package export_data;
use FindBin qw($Bin);
use strict;
use Data::Printer;
use Data::Dumper;
use JSON::XS;
use export_excel;
use Spreadsheet::WriteExcel; 
use Spreadsheet::WriteExcel::Utility;               # Import everything


use Spreadsheet::WriteExcel;
sub print {
	my ($project,$cgi,$data,$type_identifier) = @_;
	
	$type_identifier = "name" unless defined $type_identifier;
	print_bed($cgi,$project,$data) if ($cgi->param('ucsc') == 1);
	my %all;
#	$all{identifier} = "name";
	$all{label}      = $type_identifier;
	#$all{identifier} ="id";
	
	if ($type_identifier eq "uniq") {
		my $zz =1111;
		map {$_->{uniq} = $zz++;} @$data;+
	   
		$all{identifier} = $type_identifier;
	}
	elsif ($type_identifier eq "forest"){
		
		$all{label}      = "name";
		$all{identifier} = "id";#$id;
	}
	
	$all{items}      = $data;
	if ( $cgi->param('xls') == 1 ) {
		print "Content-type: application/msexcel\n";
		print "Content-Disposition: attachment;filename=".$project->name().".xls\n\n";
		my $type = $cgi->param('type');
		export_excel::exportTable($data,$type,$project);
	}
	else {
		
		print $cgi->header('text/json-comment-filtered');
		print encode_json \%all;
		print "\n";
	}
	
	exit(0);
}

my $first;
sub print_header_stream {
		my ($cgi) = @_;
		print $cgi->header('text/json-comment-filtered');
		print qq{\{"label":"id","items":[};
		
}

sub print_stream {
		my ($data) = @_;
	
		foreach my $d (@$data){ 
				if ($first){
			print ",";
		}
		print encode_json $d;
		$first ++;
		}
}
sub print_end_stream {
		
		print qq{]\}};
		
}

sub print_scalar {
	my ($project,$cgi,$data,$type_identifier) = @_;
	my $pp;
	$type_identifier = "name" unless defined $type_identifier;
	
	my %all;
	#$all{identifier} = "name";
	$all{label}      = $type_identifier;
	$all{identifier} ="id";
	if ($type_identifier eq "uniq") {
		my $zz =1111;
		map {$_->{uniq} = $zz++;} @$data;
	
		$all{identifier} = $type_identifier;
	}
	
	$all{items}      = $data;
	$pp.=  $cgi->header('text/json-comment-filtered');
	$pp.=  encode_json \%all;
	$pp.=  "\n";
	
	
	return $pp;
	
}

sub print_json {
	my ($cgi,$data) = @_;
		print $cgi->header('text/json-comment-filtered');
		print encode_json $data;
		print "\n";
	
	exit(0);
} 

sub print_simpleJson {
	my ($cgi,$data,$type_identifier) = @_;
	
	$type_identifier = "name" unless defined $type_identifier;
	my %all;
	#$all{identifier} = "name";
	$all{label}      = $type_identifier;

	
	$all{items}      = $data;		
		print $cgi->header('text/json-comment-filtered');
		print encode_json \%all;
		print "\n";
	
	exit(0);
}

sub print_bed {
	my ($cgi,$project,$data) = @_;
 
	print $cgi->header();
	my $name = $project->name();
	my $desc = "description=\"$name\"";
	if ($project->is_ngs){
		print "track name=".$name." ".$desc."color=0,0,200 \n";
	}
	elsif ($project->is_cnv){
		print "track name=".$name." ".$desc." color=200,0,200  \n";
	}
	else{
			print "track name=".$name." ".$desc." color=200,10,0 \n";
	}
	 
	
	
	#print "browser position chr$chr:$start-$end\n";
	foreach my $v  (@{$data}){
	
		my $chr_name = "chr".$v->{chromosome};
		$chr_name = "chrM" if $v->{chromosome} eq "MT";
		
		
		print $chr_name."\t".($v->{start})."\t".($v->{end})."\t".$v->{name}."\t".$v->{bed_score}."\n";
	}	
		
	

	exit(0);
	
}


sub printWithLabel {
	my ($project,$cgi,$data) = @_;
	my %all;
	$all{identifier} = "name";
	$all{label}      = "label";
	$all{items}      = $data;
	if ( $cgi->param('xls') == 1 ) {
		print "Content-type: application/msexcel\n";
		print "Content-Disposition: attachment;filename=".$project->name().".xls\n\n";
		exportTable($data);
	}
	else {
		print $cgi->header('text/json-comment-filtered');
		print encode_json \%all;
	}
	exit(0);
}


sub update_deja_vu {
	my ( $project,$data,$user) = @_;
	my $t;
	# map {push(@{$t->{$_->{chromosome}}},$_->{id})} @$data;

	
	foreach my $d (@$data) {
		my $nb_project= 0;
		my $nb_patient= 0;
		my $he =0;
		my $ho =0;
		my $list = $project->getDejaVuInfos($d->{id});
		delete $list->{$project->name};
		my $st_project;
		foreach my $oproject (keys %$list){
			$nb_project ++;
			$st_project=$oproject.":";
			$st_project .= $list->{$oproject}->{patients};
			$he += $list->{$oproject}->{he};
			$ho += $list->{$oproject}->{ho};
			
			
		}
		$d->{bipd_db} = "-";
		$nb_patient = $he + $ho;
		if ($nb_project > 0) {
			$d->{bipd_db} = $nb_project."/".$nb_patient."\nho:$ho,he:$he|".$d->{id};
			if ($user eq "hgidnyc" || $user eq "pnitschk" || $user eq "masson"){
				$d->{bipd_db} .= ":".$st_project if $nb_project <  5;
				
			}
		}

		
	}
	

}
sub update_valid {
	my ($project, $data ,$user) = @_;

	my $vquery = validationQuery->new(dbh=>$project->buffer->dbh,capture_name=>$project->capture->[0]->{name});
	
	foreach my $d (@$data) {
		my $valid =0;
		my $invalid = 0;
		my $notseq=0;
		my $ho = 0;
		my $he = 0 ;
		my $globalValid = 0;
		my $nbpatients = scalar (@{ $d->{patient_name} });
		my $validation_vid = $vquery->getVariationByGenBoId(id=>$d->{id});
		if ($validation_vid){
			foreach my $patient_name (@{ $d->{patient_name} }){
				my $value_valid = $vquery->getValidations(id=>$validation_vid,project=>$project->name(),sample=>$patient_name);
				
				$d->{ "valid!" .$patient_name} = $value_valid if $value_valid;
				$d->{valid} = "1" if $value_valid;
			}
		}
		
	}
	
}



sub variation_report_xls{
my ($project,$full_data) = @_;
$| = 1;
print "Content-type: application/msexcel\n";
print "Content-Disposition: attachment;filename=".$project->name().".xls\n\n";
my $line_sheet2;

my $calling_methods = $project->getCallingMethods();
my $workbook  = Spreadsheet::WriteExcel->new(\*STDOUT);
#my $workbook  = Spreadsheet::WriteExcel->new("/data-xfs/dev/pnitschk/svn-genbo/polymorphism-cgi/json_output/toto.xls");
my $worksheet = $workbook->add_worksheet();
my $worksheet2 =  $workbook->add_worksheet();
my $row_sheet2 = 1;
my @patients = map {$_->name} @{$project->getPatients};
my $nb_patients = scalar(@patients);
my @header = ("variation","type","consequence","dejavu","chr","position","allele","sequence",@patients,"he","ho",
"gene","consequence","transcript","transcript xref","description","exon","cdna pos","cds pos",
"protein","protein_xref","AA","protein pos","nomenclature","Polyphen","sift"); 
#my @col_name = ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","W","X","Y","Z"AAABACAD);
my $nb_col = scalar(@header);
my $nb_line = scalar(@$full_data);

my $bold       = $workbook->add_format(bold => 1);
my $row = 0;
#$worksheet->set_row(0,$nb_col,$bold);
my $bold_merge ;
 $bold_merge->[0]      = $workbook->add_format(valign      => 'vcentre',
                                        align       => 'centre',
                                        );
  $bold_merge->[1]      = $workbook->add_format(valign      => 'vcentre',
                                        align       => 'centre',
                                        );   
my $bg_color;

 $bg_color->{red} = $workbook->add_format(valign      => 'vcentre',
                                        align       => 'centre',
                                        bg_color => 'red',
                     );                                 
   $bg_color->{orange} = $workbook->add_format(valign      => 'vcentre',
                                        align       => 'centre',
                                        bg_color => 'orange',
                     );       
   $bg_color->{blue} = $workbook->add_format(valign      => 'vcentre',
                                        align       => 'centre',
                                        bg_color => 'blue',
                   );                      
   $bg_color->{green} = $workbook->add_format(valign      => 'vcentre',
                                        align       => 'centre',
                                        bg_color => 'green',
                   );   
   $bg_color->{gray} = $workbook->add_format(valign      => 'vcentre',
                                        align       => 'centre',
                                        bg_color => 'gray',
                   ); 
      $bg_color->{silver} = $workbook->add_format(valign      => 'vcentre',
                                        align       => 'centre',
                                        bg_color => 'silver',
                   );  
       $bg_color->{cyan} = $workbook->add_format(valign      => 'vcentre',
                                        align       => 'centre',
                                        bg_color => 'cyan',
                   );                 
                                                                                                                                   
                                        
my ($col1,$line1) = split("",xl_rowcol_to_cell(0,0)); 
my ($col2,$line2) = split("",xl_rowcol_to_cell(5+$nb_patients,0)); 
$worksheet->set_column($col1.":".$col2, 12);
$worksheet2->set_column($col1.":".$col2, 12);
($col1,$line1) = split("",xl_rowcol_to_cell(6+$nb_patients,0)); 
($col2,$line2) = split("",xl_rowcol_to_cell(7+$nb_patients,0));
 $worksheet->set_column($col1.":".$col2, 5);
 $worksheet2->set_column($col1.":".$col2, 12);
 
($col1,$line1) = split("",xl_rowcol_to_cell(8+$nb_patients,0));
($col2,$line2) = split("",xl_rowcol_to_cell(13+$nb_patients,0));

 $worksheet->set_column($col1.":".$col2, 30);
 $worksheet2->set_column($col1.":".$col2, 12);


$worksheet->write(0,0, \@header);
$worksheet2->write(0,0, \@header);
$row++;
my %test;
my $total = scalar(@$full_data);
my $nbv =0;
my $seq_in;
$seq_in = 1 if $total < 7000;
foreach my $data (@$full_data){
	$nbv++;
	my $col = 0;
	my $nb_merge = 0;
	my @first_part;
	foreach my $gene (@{$data->{genes}}){
		
		my @trans = grep {$_->{gene}=~/$gene/} @{$data->{tab_consequences}};
		if  (scalar(@trans)){
			$nb_merge += scalar(@trans);
		}
		else {
			$nb_merge++;
		} 
	}
	$nb_merge= 1 if $nb_merge == 0;
	push(@first_part,$data->{name});
	printWorksheet($worksheet,$row,$col,$nb_merge,$data->{name},$bold_merge);
	$test{$row." ".$col}++;
	$col++;
	printWorksheet($worksheet,$row,$col,$nb_merge,$data->{type},$bold_merge);
	$test{$row." ".$col}++;
	$col++;
	push(@first_part,$data->{type});
	my $h = "consequence!all";
	printWorksheet($worksheet,$row,$col,$nb_merge,$data->{$h},$bold_merge);
	$test{$row." ".$col}++;
	push(@first_part,$data->{$h});
	$col++;
	printWorksheet($worksheet,$row,$col,$nb_merge,$data->{bipd_db},$bold_merge);
	push(@first_part,$data->{bipd_db});
	$test{$row." ".$col}++;
	$col++;
	printWorksheet($worksheet,$row,$col,$nb_merge,$data->{chromosome},$bold_merge);
	push(@first_part,$data->{chromosome});
	my $chr = $project->getChromosome($data->{chromosome});
	$test{$row." ".$col}++;
	$col++;
	printWorksheet($worksheet,$row,$col,$nb_merge,$data->{start},$bold_merge);
	push(@first_part,$data->{start});
	$test{$row." ".$col}++;
	$col++;
	my $seq_around="";
	$seq_around = $chr->getSequence($data->{start}-20,$data->{start}+20);
	
	my $z = substr $seq_around,20,1, "[".$data->{text}."]";
	
	
	printWorksheet($worksheet,$row,$col,$nb_merge,$data->{text},$bold_merge);
	push(@first_part,$data->{text});
	$test{$row." ".$col}++;
	$col++;
	
	printWorksheet($worksheet,$row,$col,$nb_merge,$seq_around,$bold_merge);
	push(@first_part,$seq_around);
	$test{$row." ".$col}++;
	$col++;
	
	my $string = $data->{allpatients};

	foreach my $p (@patients){
		my $find =" ";
		
	for (my $i=0;$i<@{$data->{patient_name}};$i++){  
		
		next if $data->{patient_name}->[$i] ne $p;
		my $seq = $data->{patient_base}->[$i];
		#$find ="-";
		#$find = "He";
		#$find = "Ho" if lc($seq) eq "t" || lc($seq) eq "a" || lc($seq) eq "c" || lc($seq) eq "g" ; #if lc($seq) =~ /atcg/;
		
			foreach my $cl1  (@$calling_methods){
			my	$cl="unifiedgenotyper";
			$find .= "$cl ?? - " unless exists $data->{"patient_".$cl."heho"};
			next unless exists $data->{"patient_".$cl."heho"};
			$find .=  $data->{"patient_".$cl."heho"}->[$i];
			$find .="(".$data->{"patient_".$cl."2"}->[$i]."/".$data->{"patient_".$cl."3"}->[$i].")";
		}
		}
		push(@first_part,$find);
		printWorksheet($worksheet,$row,$col,$nb_merge,$find,$bold_merge);
		$test{$row." ".$col}++;
			$col++;
	}
		printWorksheet($worksheet,$row,$col,$nb_merge,$data->{heterozygote},$bold_merge);
			push(@first_part,$data->{heterozygote});
		$test{$row." ".$col}++;
	
	$col++;
	printWorksheet($worksheet,$row,$col,$nb_merge,$data->{homozygote},$bold_merge);
	push(@first_part,$data->{homozygote});
	$test{$row." ".$col}++;
	$col++;
	unless (@{$data->{genes}}){
			$worksheet2->write($row_sheet2,0, \@first_part);
			$row_sheet2 ++;
	}
	foreach my $gene (@{$data->{genes}}){
		
		my @trans = grep {$_->{gene}=~/$gene/} @{$data->{tab_consequences}};
		my @genes_part = @first_part;
		unless (scalar @trans){
			printWorksheet($worksheet,$row,$col,1,$gene,$bold_merge);
			push(@genes_part,$gene);
			$worksheet2->write($row_sheet2,0, \@genes_part);
			$row_sheet2 ++;
			$row++;
			next;
		}
		my $nb_merge_trans = scalar(@trans) ;
		$nb_merge_trans =1 unless $nb_merge_trans;
		printWorksheet($worksheet,$row,$col,$nb_merge_trans,$trans[0]->{gene},$bold_merge);
		push(@genes_part,$trans[0]->{gene});
		$test{$row." ".$col}++;
		
		
	
		
		foreach my $tt (@trans){
			my @trans_part = @genes_part;
			my $col_start = $col+1;
			my ($ens,$xref) = split(/\+/,$tt->{transcript});
					#	warn "--".$ens;
		#	my $vid = $project->getKyotoGenes($ens);

			#warn  $vid->[0];
			#die();
			my $chromosome = $data->{chromosome};
			
			my  $htr=  $project->liteObject($ens."_".$chromosome);
			
			my $h = "consequence!".$ens;
			my $color = "silver";
			if ($data->{$h} eq "utr" ){
				$color = "blue";
			}
			if ($data->{$h} =~ /splicing/){
				$color = "cyan";
			}
			if ( $data->{$h} eq "stop" || $data->{$h} eq "phase"){
				$color = "red";
			}
			if ($data->{$h} eq "coding" || $data->{$h} eq "frameshift"){
				$color = "orange";
			}
			
			my $bg =  $bg_color->{$color};
			$worksheet->write($row,$col_start ++,$data->{$h},$bg);
			push(@trans_part,$data->{$h});
			# p $tt;
			
			$test{$row." ".$col_start}++;
			$worksheet->write($row,$col_start ++,$ens,$bg);
				push(@trans_part,$ens);
			$test{$row." ".$col_start}++;
			$worksheet->write($row,$col_start ++,$htr->{external_name},$bg);
			push(@trans_part,$htr->{external_name});
			$test{$row." ".$col_start}++;
			$worksheet->write($row,$col_start ++,$htr->{gene_description});
			push(@trans_part,$htr->{gene_description});
			$test{$row." ".$col_start}++;
			##"exon","cdna pos","cds pos"
			
			$worksheet->write($row,$col_start ++,$tt->{exon});
			push(@trans_part,$tt->{exon});
			$test{$row." ".$col_start}++;
			
			$worksheet->write($row,$col_start ++,$tt->{cdna_position});
			push(@trans_part,$tt->{cdna_position});
			$test{$row." ".$col_start}++;
			
			$worksheet->write($row,$col_start ++,$tt->{cds_position});
			push(@trans_part,$tt->{cds_position});
			$test{$row." ".$col_start}++;
			## Protein ....
			($ens,$xref) = split(/\+/,$tt->{protein},$bg);
			$worksheet->write($row,$col_start ++,$ens);
				push(@trans_part,$ens);
				$test{$row." ".$col_start}++;
			$worksheet->write($row,$col_start ++,$htr->{external_protein_name});
				$test{$row." ".$col_start}++;
				push(@trans_part,$htr->{external_protein_name});
			$worksheet->write($row,$col_start ++,$tt->{AA_protein}."/".$tt->{AA_variation});
			push(@trans_part,$tt->{AA_protein}."/".$tt->{AA_variation});
			$test{$row." ".$col_start}++;
			
			$worksheet->write($row,$col_start ++,$tt->{protein_position});
			push(@trans_part,$tt->{protein_position});
			#nomenclature
			$worksheet->write($row,$col_start ++,$tt->{nomenclature});
			push(@trans_part,$tt->{nomenclature});
			## polyphen sift
			my ($status,$score) = split(/\+/,$tt->{polyphen_status});
			$color = "silver";
			$color="red" if $status ==3;
			$color = "orange" if $status == 2;
			$color = "green" if $status == 1;
			$worksheet->write($row,$col_start ++,$score,$bg_color->{$color});
			push(@trans_part,$score);
			#status
			($status,$score) = split(/\+/,$tt->{sift_status});
			 $color = "silver";
			#$color="red" if $status ==3;
			$color = "red" if $status >= 2;
			$color = "green" if $status == 1;
			$worksheet->write($row,$col_start ++,$score,$bg_color->{$color});
			push(@trans_part,$score);
			$test{$row." ".$col_start}++;
				
			$test{$row." ".$col_start}++;
			
			
			
			
			$row++;
			$worksheet2->write($row_sheet2,0, \@trans_part);
			$row_sheet2 ++;
		}
		#$row ++;
	}
	#$row ++;
}
	
 $worksheet->autofilter('E1:E'.$row);
  $workbook->close();

exit(0);
}

sub printWorksheet{
	my ($worksheet,$row,$col,$merge,$text,$bold_merge) = @_;
	
	if ($merge >1) {
		$merge --;
		$worksheet->merge_range($row,$col,$row+$merge,$col,$text,$bold_merge->[0]);
	}
	else {
		$worksheet->write($row,$col,$text,$bold_merge->[1]);
	}
}

sub printVcfFile {
	my ($projectName, $patientName, $lVcfLines) = @_;
	print "Content-type: application/text\n";
	print "Content-Disposition: attachment;filename=$projectName\_$patientName.vcf\n\n";
	my $nbResPrinted = 0;
	foreach my $line (@$lVcfLines) {
		print $line;
		$nbResPrinted++;
	}
	if ($nbResPrinted == 0) { print "NO RESULT for $patientName and this filters...\n"; }
	exit(0);
}

1;