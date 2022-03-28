#!/usr/bin/perl
use CGI qw/:standard :html3/;

use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/kyoto";

#use lib "/bip-d/soft/distrib/tabix/latest/perl";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/coverage";

#use Set::;
use Carp;
use export_data; 
use strict;
use Set::IntSpan::Fast::XS;
use Data::Dumper;
use GBuffer;
use Getopt::Long;
use Carp;
use Set::Intersection;
use Tabix;
use Storable qw/thaw/;
use Data::Dumper;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use coverage; 
use Spreadsheet::WriteExcel;
use POSIX;
use validationQuery;

my $nb_gene_by_patient = 3;
my $nb_exon_by_genes = 10;

my $buffer = GBuffer->new();

my $cgi          = new CGI();
my $project_name = $cgi->param('project');

my $project = $buffer->newProject(-name=>$project_name);
my $vquery = validationQuery->new(dbh=>$buffer->dbh);
my $patient_name = $cgi->param('patients');
$patient_name ="all" unless $patient_name;

my $sql = qq{SELECT validation_id ,sample_name as patient_name ,vcfid as variation_name ,validation  FROM validation.validations v, validation.variations vr where v.project_name="NGS2013_0201" and v.variation_id=vr.variation_id; };

my $sth = $buffer->dbh->prepare($sql);
$sth->execute();
my $hres = $sth->fetchall_hashref("validation_id");

my $validations;
foreach my $id (keys %$hres){
	my $var;
	($var->{chr},$var->{start},$var->{ref},$var->{alt}) = split("_",$hres->{$id}->{variation_name});
	$var->{chr_name} =~ s/chr//;
	$var->{validation} = $hres->{$id}->{validation};
	push(@{$validations->{$hres->{$id}->{patient_name}}},$var);
	
}


my $patients = $project->get_list_patients($patient_name," ");

my $hpatient;

my $coverage_dir = $project->getRootDir()."/align/coverage/";
my $hjson;
my $patient_name =  $cgi->param('patient');
#2:227867427-228028829
my $tr_index = "ENST00000396625";
my @transcripts_cgi = split(" ",$cgi->param("transcripts"));


my $exon_type = "genomic_span";
#
#my $workbook  = Spreadsheet::WriteExcel->new(\*STDOUT);
#my $worksheet = $workbook->add_worksheet();
#my $mergef      = $workbook->add_format(valign      => 'vcentre',
#                                        align       => 'centre',
#                                        );
#my $format      = $workbook->add_format(valign      => 'vcentre',
#                                        align       => 'centre',
#                                        );
my $splice_5 = 10;
my $limit    = 0;

$limit = $cgi->param('limit');
$splice_5 = $cgi->param('span');
$splice_5 = 1 unless $splice_5;
my $apos;
my $coding;



my $row = 1;
my $sum;
my $nb_exons;
my $patients_data;
my $transcripts_exons;
my $capture = coverage::get_capture($project);
#foreach my $gene_name (@genes_name){
	
#	my $gene = thaw  $chr->getKyotoHash("ensembl","gene")->get($gene_name);
#	next unless $gene->{transcripts}; 

my @chromosomes = (1..22,'X','Y','MT');
my $chromosomes = $project->getChromosomes();

my @transcripts;
foreach my $chr (@$chromosomes){
	
	
	foreach my $tr_name  (@transcripts_cgi){
			my $tr1 =   $chr->getKyotoHash("ensembl","transcript")->get($tr_name);
			next unless $tr1;
			
			my $tr = thaw $tr1;
			 $tr->{ochromosome} = $chr;
			$tr->{name} = $tr_name;
			push(@transcripts,$tr);
	}
}


	foreach my $htr  (@transcripts){
		my $exons;
		my $tr;
		
		 ($exons,$tr) = coverage::get_exons($htr->{ochromosome},$htr->{name},$splice_5);
		$transcripts_exons->{$htr->{name}} = $exons;
		coverage::get_statistics( $htr->{name}, $exons, $patients,$splice_5 );
		coverage::get_all_span( $htr->{name}, $exons, $patients,
					$splice_5, $limit, $capture );
	}
	
	my $hresults;
	
	foreach my $patient (@{$patients}){
	my $hpatient;
	my $validations_patient = [];
	 $validations_patient = $validations->{$patient->name()} if exists $validations->{$patient->name()};
	$hpatient->{name} = $patient->name();
	
	foreach my $tr (@transcripts){
		my $htranscript;
		
		my $tr_name = $tr->{name};
		my $variations = [];
		
		foreach my $v (@$validations_patient){
			next unless $v->{chr_name} ne $tr->{chromosome};
			next unless $tr->{genomic_span}->contains($v->{start});
			push(@$variations,$v);
		}
		
		my $exons = $transcripts_exons->{$tr_name};
		$htranscript->{lname} = $tr->{gene}." $tr_name :".$tr->{gene_external_name};
		$htranscript->{name} = $tr_name;
		my $exons = $transcripts_exons->{$tr_name};
		my $nb_exons =0;
		
		foreach my $exon (@$exons){
			my $hexons;
			$hexons->{name} = $exon->{name};
			$hexons->{id} = $exon->{id};
			#$hexons->{label} = $tr->{gene_external_name}.":".$tr->{name}.":".$exon->{id}.":".$tr->{chromosome}.":".$exon->{gstart}.":".$exon->{gend};
			$nb_exons++;
			$hexons->{label}->{$patient->name()} = $patient->name().":".$tr->{gene_external_name}.":".$tr->{name}.":".$exon->{id}.":".$tr->{chromosome}.":".$exon->{gstart}.":".$exon->{gend};
			
			 $hexons->{todo} = 1 if  $vquery->is_todo(project_name=>$project_name,id=>$hexons->{label}->{$patient->name()});
			
			my $ename  = $exon->{name};
			my $variation_inside;
			foreach my $v (@$variations){
				$variation_inside= 1 if $v->{start}>= $exon->{start} && $v->{start} <= $exon->{end};
			}
			my @nb = ();
			my $debug;
		
			
			 @nb = $patient->{_span}->{$tr_name}->{$ename}->as_array() if $patient->{_span}->{$tr_name}->{$ename};
			 $hexons->{color} = "green";
			 
			my $color = "tdgreen";
			if (scalar(@nb)> 0 ){
				 $color = "tdred";
				 $hexons->{color} = "red";
			}
			$color = "tdviolet" if $variation_inside;
			
			$hexons->{color} = "violet" if $variation_inside;
			$hexons->{color} = "strike" if $exon->{utr};
			$hexons->{color} = "blue" if exists $hexons->{todo};
			
			push(@{$htranscript->{exons}},$hexons);
		}
		push(@{$hpatient->{transcripts}},$htranscript);

	}
	push(@{$hresults->{patients}},$hpatient);
}

unless ( $cgi->param('xls')){
	html();
	exit(0);
}

xls();	
sub html {	
	
my $nb_case;
print $cgi -> header;
my $w = "0px";
my $CSS = <<END;
<style type="text/css"> 
.table1 {
	align:top;
	font-size:12px;
	height:100px;
	
	
}
.table2 {
	vertical-align:bottom;
	
}
.thtop{
	
}
.tdgreen {
	padding:$w;
	border: 1px solid black;
	background-color:#009900;
}
.tdblue {
	padding:$w;
	background-color:#5CB3FF;
}

.tdorange {
	padding:$w;
	background-color:#FFA500;
}
.tdviolet {
	padding:$w;
	background-color:#8D38C9;
}
 	
.tdred{
	padding:$w;
	background-color:#990000;
}
.tdstrike{
	padding:$w;
	background-color:#101010;
}
.check{
	
}
.label{
	font-size:9px;
}
</style>

END

print $CSS;
my $nb_col =5;
my $color1 = "#009966";
my $color2 = " #990000";

my $modulo_patient =5;
my $nb_patient =0; 



print $cgi->start_table({class=>"table2"});
my $jsid = 0;
foreach my $patient (@{$hresults->{patients}}){
	print $cgi->start_Tr();
	print $cgi->td({colspan=>$nb_gene_by_patient,align=>"center"},$patient->{name});
	print $cgi->end_Tr();
	print $cgi->start_Tr();
	foreach my $tr (@{$patient->{transcripts}}){
		print $cgi->start_td();
		print $cgi->start_table({caption=>$tr->{name}, class=>"table1"});
		print $cgi->th({colspan=>$nb_exon_by_genes},$tr->{lname});
		print $cgi->start_Tr();
		my $nb_exons = 0;
			foreach my $exon (sort{$a->{id} <=> $b->{id} }@{$tr->{exons}}){
			my $color = "td".$exon->{color};
			$nb_exons++;
			print $cgi->start_td({class=>$color});
			my $value = $exon->{label}->{$patient->{name}} ;#$patient->{name}.":".$exon->{label};
			my $label =$exon->{name}; 
			$jsid++;
			my $check = qq{
				<input id="$jsid" name ="$jsid" class ="check" data-dojo-type="dijit/form/CheckBox" value="$value"  onChange="test(this)" /> 
				<label for="mycheck" class="label">$label</label>
			};
			#print $check."\n";
#			my $value = $patient->{name}.":".$exon->{name};
#			my $toto = $cgi->checkbox(-name=>'checkbox_name',
#			   -checked=>0,
#			   -value=>$exon->{name},
#			   -label=>$exon->{name},
#			   -class=>"check",
#			   -onClick=>"test(this)");
			if ($exon->{color} eq "red"){
			print $cgi->div({style=>"color:white"},$check) ;
			}
			else {
			print $cgi->div({style=>"color:white"},$label) ;
			}
			print $cgi->end_td;
			if ($nb_exons%$nb_exon_by_genes == 0) {
				print $cgi->end_Tr();
				print $cgi->start_Tr();
				
			}
		}
		print $cgi->end_Tr();
		print $cgi->end_table;
		print $cgi->end_td(); 
		
	}
	
	print "\n";
	
	print $cgi->end_Tr();
	
}

print $cgi->end_table();
exit(0);
}

sub xls{

	print "Content-type: application/msexcel\n";
	print "Content-Disposition: attachment;filename=".$project->name().".xls\n\n";
	my $workbook  = Spreadsheet::WriteExcel->new(\*STDOUT);
	my $worksheet = $workbook->add_worksheet();
	my $bg_color;
	my @colors = ("red","orange","blue","green","cyan","gray");

                   
	foreach my $c (@colors){
		$bg_color->{$c} = $workbook->add_format(valign      => 'vcentre',
                                        align       => 'centre',
                                        #bg_color => 'white',
                                        color => $c,
                                        bold=>1,
                     );                                 
	}
	$bg_color->{strike} = $workbook->add_format(valign      => 'vcentre',
                                        align       => 'centre',
                                        #bg_color => 'white',
                                        color => "gray",
                                        bold=>1,
                                        font_strikeout=>1,
                     );    
                     
	my $row_patient = 1;
	my $col = 0;
	my $row =0;
	$worksheet->write($row,$col,"Coverage :".$limit." span : ".$splice_5);
	$row++;
foreach my $patient (@{$hresults->{patients}}){
	$col =0;
	$worksheet->write($row,$col,$patient->{name});
	#$row++;
	$col++;
	my $start_row = $row;
	my $row_next_patient;
	my $bold_merge;
	$bold_merge->[0]      = $workbook->add_format(valign      => 'vcentre',
                                        align       => 'centre',
                                        );
	foreach my $tr (@{$patient->{transcripts}}){
		
		#$worksheet->write($start_row,$col,$tr->{name});
		$worksheet->merge_range($start_row,$col,$start_row,$col+$nb_exon_by_genes -1,$tr->{lname},$bold_merge->[0]);
		my $row_exon = $start_row + 1;
	
		
		my $nb_exons = 0;
		my $start_col = $col;
		foreach my $exon (sort{$a->{id} <=> $b->{id} }@{$tr->{exons}}){
			my $color = $exon->{color};
			$nb_exons++;
			$color ="cyan" if $color eq "violet";
			$worksheet->write($row_exon,$col++,$exon->{name},$bg_color->{$color});
			if ($nb_exons % $nb_exon_by_genes == 0) {
				$row_exon ++;
				$col = $start_col;
				$row_next_patient= $row_exon if $row_next_patient < $row_exon;
			}
		}
		
		$col = $start_col + $nb_exon_by_genes +1;
	}
	
	$row = $row_next_patient +1;
}
            
                   
} 
