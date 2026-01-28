#!/usr/bin/perl

use Carp;
use strict;
use JSON;
use Data::Dumper;
use CGI qw/:standard :html3/;
use Set::IntSpan::Fast::XS;
use Set::IntervalTree;
use List::Util qw[min max];
use FindBin qw($Bin);
use Storable qw(store retrieve freeze);
 
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";

use GBuffer;
use GenBoProject;

my $compteur;
my $cgi = new CGI;
my $projectname = $cgi->param('project');
my $patientname = $cgi->param('patient');
my $chr         = $cgi->param('chrom');
#my $qual = $cgi->param('qual');
my $dejavu      = $cgi->param('dejavu');
my $listOfGenes = $cgi->param('genes');
my $omim      = $cgi->param('omim');
my $transmis = $cgi->param('transmission');
my $scoreEvent = $cgi->param('score_event');
#

my $type;

my $hfilteredEvent;

my @listHashRes;
my $hbreakpoint;
my $hallPatient;


# creer un objet projet
my $buffer = GBuffer->new();	
my $project = $buffer->newProjectCache( -name => $projectname);

# pour distinguer parmi les autres patients du projet
# ceux qui sont de la même famille
my $thePatient = $project->getPatient($patientname);
my $patient = $thePatient;
my $thePatientFamily = $thePatient->getFamily();

my $dir = $project->getCacheSV(). "/rocks/";
my $rocks = GenBoNoSqlRocks->new(dir=>"$dir",mode=>"r",name=>"sv");

#my $file_in = $thePatient->getSRFile;
#my $tabix = Bio::DB::HTS::Tabix->new( filename => $file_in );
	
my $trio;
my $mother;
my $father;
if ($patient->isChild){
	 $mother = $patient->getFamily->getMother();
	 $father = $patient->getFamily->getFather();
	
}

my $hintpsan = {};
warn $listOfGenes;
if ($listOfGenes ne "all"){
	foreach my $chr (@{$project->getChromosomes}){
		$hintpsan->{$chr->name} = Set::IntSpan::Fast::XS->new();
	}
	Set::IntSpan::Fast::XS->new();
	foreach my $gname (split(",",$listOfGenes)){
		my $g = $project->newGene($gname);
		my $chr = $g->getChromosome();
		$hintpsan->{$chr->name}->add_range($g->start,$g->end);
		
	}
}

#my $dbh = DBI->connect("dbi:ODBC:Driver=DuckDB;Database=:memory:", "", "", { RaiseError => 1 , AutoCommit => 1});
my $parquet_file = $project->getCacheSV()."/".$project->name.".".$project->id.".parquet";
my $dir = $project->getCacheSV(). "/rocks/";
my $rocks = GenBoNoSqlRocks->new(dir=>"$dir",mode=>"r",name=>"sv");



#my $patient = $project->getPatients()->[0];
my $patient = $project->getPatient($patientname);
my $patient_id = $patient->id;
my $colpatient = "p".$patient->id;
my $query = qq{CREATE TABLE SV  AS
                           SELECT * 
                           FROM '$parquet_file'
                           WHERE pos1 > -1  and patient = $patient_id  ;
	};
	
#$dbh->do($query);
my $sql =qq{select * from '$parquet_file' where patient = $patient_id and nb_dejavu_patients <= $dejavu ; };
my $cmd = qq{duckdb -json -c "$sql"};
my $res =`$cmd`;
my $array_ref = [];
$array_ref  = decode_json $res if $res;
	
#my $sth = $dbh->prepare($sql);
#$sth->execute();
$dejavu =20;
#	while (my $row = @{$array_ref}) {
	foreach my $row ( @{$array_ref}) {
	 my $sv = $rocks->get($row->{id});
	next unless $sv;
	 next if $sv->{dejavu}->{nb_patients} > $dejavu;
	 delete $sv->{score};
	 getScoreEvent($sv);
	  my $hash;
	my $name = "t(".$sv->{chrom1}.",".$sv->{chrom2}.")(".$sv->{cytoband1}.",".$sv->{cytoband2}.")";
	 $hash->{"LENGTH"}=  "-";
	 if ($sv->{type} eq "INV") {
	 	$name = "inv(".$sv->{chrom1}.")(".$sv->{cytoband1}.",".$sv->{cytoband2}.")";
	 	$hash->{"LENGTH"}=  $sv->{type}.";".abs($sv->{pos1}  -$sv->{pos2});
	 	 if (exists $hintpsan->{$sv->{chrom1}}){
	 	 	next if  $hintpsan->{$sv->{chrom1}}->is_empty;
	 		my $intspan = Set::IntSpan::Fast::XS->new($sv->{pos1}."-".$sv->{pos2});
	 		my $ai = $hintpsan->{$sv->{chrom1}}->intersection($intspan);
	 		next if $ai->is_empty;
		 }

	 }
	 else {
	 	
	 	my @chr1 = ($sv->{chrom1},$sv->{chrom2});
	 	my @pos = ($sv->{pos1},$sv->{pos2});
	 	my $found;
	 	for (my $i=0;$i< @chr1;$i++ ){
	 		next unless exists  $hintpsan->{$chr1[$i]};
	 		my $intspan = Set::IntSpan::Fast::XS->new(($pos[$i]-2000)."-".($pos[$i]+2000));
	 		my $ai = $hintpsan->{$chr1[$i]}->intersection($intspan);
	 		$found ++ unless  $ai->is_empty;
	 	}
	 	next unless $found;
	 }
		#$name = "INV(".$sv->{chrom1}."[".$sv->{pos1}."-".$sv->{pos2}."])" if ($sv->{type} eq "INV");
	 
	
	
	$hash->{TRANSLOC} = $sv->{type}.";".$name;
	$hash->{"id"}= $sv->{id};
	$hash->{"QUAL"}=$sv->{qual};
	
	
	
	$hash->{"POS1"}= $sv->{pos1};
	$hash->{"CYTOBAND1"}=  $sv->{type}.";".$sv->{cytoband1};
	
	
#	$hash->{"REF/ALT1"}= "-";
	$hash->{"TYPE"}= $sv->{type};
	$hash->{"POS2"}= $sv->{pos2};
	$hash->{"CYTOBAND2"}=  $sv->{type}.";".$sv->{cytoband2};
	##################
	## GENES 
	##################
	$hash->{"ALLGENES"}= "-";	
	my $before = $sv->{chromosome}.":".$sv->{pos1}."_".$sv->{pos2}."##";
	my $debug;
	$debug=1 if $sv->{pos1} == 90072610;
	$hash->{"GENES1"}= listGenes($sv->{genes},$sv->{id},$debug);"-";#$sv->{genes1};
	##################
	## OMIM
	##################
	$hash->{"OMIM2"}= $sv->{omim2};
		$hash->{"OMIM1"}= $sv->{omim1};
	if ($omim){
		next if $sv->{omim1} + $sv->{omim2} == 0;# and $sv->{omim1};
	}
	
	#IGV 
	my $members = $thePatientFamily->getMembers();
	my $bamNames = $patientname;
	my $bamFiles = $patient->bamUrl(); 
		foreach my $m (@$members)		{
			my $membername = $m->name(); 
			$bamFiles .= ",".$m->bamUrl() unless ($membername eq $patientname);
			$bamNames .= ",".$membername unless ($membername eq $patientname);
		}
	$hash->{'IGV'} = $bamFiles.";".$bamNames."et".$sv->{chrom1}."_".$sv->{pos1}."_".$sv->{chrom2}."_".$sv->{pos2};
	
	$hash->{"dejavuBP1"} = "-";
	$hash->{"dejavuBP2"} = "-";
	$sv->{dejavu}->{id} = $sv->{id};
	$sv->{dejavu}->{project} = $project->name();
	$sv->{dejavu}->{nb_projects} +=0;
	$sv->{dejavu}->{nb_patients} +=0;
	$sv->{dejavu}->{nb_itp} +=0;
	$sv->{dejavu}->{nb_itp} = scalar(keys %{$sv->{dejavu}->{this_project}});
	$sv->{dejavu}->{total_itp} = scalar(@{$project->getPatients});
	$hash->{"dejavu"} = encode_json $sv->{dejavu};#$sv->{DEJAVU}->{STRING};
	$hash->{'TRANSMISSION'} =  getTransmission($sv,$mother,$father); 
	$hash->{'CNV'}= "-";#$sv->{CNV};
	$hash->{'CNV_mother'}= "#";#$sv->{CNV_mother};
	$hash->{'CNV_father'}= "#";#$sv->{CNV_father};
#	
	$hash->{"SCORE_EVENT"} = $sv->{score};
	$hash->{"PR"} = $sv->{pr1}."/".$sv->{pr2};
	$hash->{"SR"} =  $sv->{sr1}."/".$sv->{sr2};
	push( @listHashRes, $hash);


}	
#<span onclick="update_grid_gene_phenotypes('ENSG00000006377_7', 'NGS2019_2653')" ;"="">Autism | Autism spectrum disorder | Cochlear incomplete partition, typ... <span style="color:#5cf0d3">+ 16 terms</span></span>
#
#warn "coucou";

if ( (scalar(@listHashRes) == 0) )
{
	my $hash;
	$hash->{"id"}= "-";
	$hash->{"TRANSLOC"}= "-";
	$hash->{"IGV"}= "-";
	$hash->{"QUAL"}= "-";
	$hash->{"ALLGENES"}= "-";	
	$hash->{"LENGTH"}= "-";
	$hash->{"POS1"}= "-";
	$hash->{"CYTOBAND1"}= "-";
	$hash->{"GENES1"}= "-";
	$hash->{"OMIM1"}= "-";
	$hash->{"REF/ALT1"}= "-";
	$hash->{"TYPE"}= "-";
	$hash->{"POS2"}= "-";
	$hash->{"CYTOBAND2"}= "-";
	$hash->{"GENES2"}= "-";
	$hash->{"OMIM2"}= "-";
	$hash->{"REF/ALT2"}= "-";
	$hash->{"dejavuBP1"} = "-";
	$hash->{"dejavuBP2"} = "-";
	$hash->{"dejavu"} = "-";
	$hash->{'TRANSMISSION'} = "-";
	$hash->{'CNV'}= "-";
	$hash->{'CNV_mother'}= "-";
	$hash->{'CNV_father'}= "-";
	$hash->{"SCORE_EVENT"} = "-";
	$hash->{"PR"} = "-";
	$hash->{"SR"} = "-";
	push( @listHashRes, { %{$hash} } );
}

#close($fd);
printJson( \@listHashRes );

exit(0);


############
#          #
# methodes #
#          #
############


sub printJson {
	my ($listHash) = @_;
	my $hash;
	$hash->{'identifier'} = 'id';
	$hash->{'label'} = 'id';

	my @t = sort {$b->{SCORE_EVENT} <=> $a->{SCORE_EVENT} or $a->{CHROM1} <=> $b->{CHROM1}or $a->{CHROM2} <=> $b->{CHROM2}} @$listHash;
	
	$hash->{'items'} = \@t;
	print $cgi->header('text/json-comment-filtered');
	print encode_json $hash;
	print "\n";
}

sub getScoreEvent
 {
	my ($sv) = @_;
	
	my $score = 0;
	
	# break point dans un gene omim morbid
	
	#27374152
	my $debug;
	$debug=1 if $sv->{pos1} == 90072610;
	
	$score += 1  if $sv->{"max_score_gene"} >= 5 ;	
	$score += 1  if $sv->{"max_score_gene"} >= 4 ;
	$score += 0.5  if $sv->{"max_score_gene"} > 3 ;	
	if ($sv->{type} eq "INV" && scalar(@{$sv->{genes}}) == 0){
		$score -= 5;
	}
	# score sur le dejavu
	my $nbdejavu = $sv->{dejavu}->{nb_patients};#$sv->{"nbdejavu"};
	warn $nbdejavu if $debug;
	if ($nbdejavu <= 10)
	{
		$score +=0.5;
	}
	if ($nbdejavu <= 5)
	{
		$score +=0.5;
	}
	if ($nbdejavu ==0)
	{
		$score ++;
	}
	warn $score if $debug;
	# score sur la QUAL
	$score += int($sv->{qual} /100)/10;
	warn $score if $debug;
	my $dsr = ($sv->{sr2}/($sv->{sr2}+$sv->{sr1}+1));
	my $dpr = ($sv->{pr2}/($sv->{pr2}+$sv->{pr1}+1));
	$score += $dsr;
	$score += $dpr;
	#die() if $debug;
	$score = sprintf("%.1f", $score);
	#pour présenter en premier les CNV denovo
#	if ( ($sv->{PATIENTS}->{$patient->id}->{'TRANSMISSION'}  =~ m/denovo/) && !($sv->{PATIENTS}->{$patient->id}->{'TRANSMISSION'} =~ m/mother/ ) && !($sv->{PATIENTS}->{$patient->id}->{'TRANSMISSION'}  =~ m/father/ ) && !($sv->{PATIENTS}->{$patient->id}  =~ m/both/) )
#	{
#			$score += 1.5; 
#	}
$sv->{score} = $score;
	return $score;
 }
  sub listGenes {
 	my ($agenes,$id,$debug) = @_;
 	my $debug ;
 	my @merged;

	my $res;
	$res->{genes} = [];
	$res->{nb_genes} = scalar(@$agenes);
	my $n =0;
	foreach my $g (@$agenes){
		delete $g->{phenotypes};
		push(@{$res->{genes}},$g);
		last if $n >= 4;
		$n ++;
	}
	$res->{project} = $project->name;
	$res->{id} = $id;
	$res->{patient} = $patient->name;
	$res->{nbh} = 0;
	$res->{nbm} = 0;
	$res->{nbl} = 0;
	
 	return encode_json $res;
 }
# sub listGenes1 {
# 	my ($agenes,$id,$debug) = @_;
# 	warn $id;
#
# 	my @genes_list_bp0 = sort{$b->{score} <=> $a->{score}}  @{$agenes};
# 	my @merged1;
#	while (@genes_list_bp1 or @genes_list_bp2 ) {
# 		push @merged, shift @genes_list_bp1 if @genes_list_bp1;
#    	push @merged, shift @genes_list_bp2 if @genes_list_bp2;
#	}
#	push (@merged, @genes_list_bp0);
#	@merged = sort{$b->{score} <=> $a->{score}} @merged;
#	my $res;
#	$res->{nb} = scalar(@merged);
#	my $n = 0;
#	foreach my $g (@merged){
#		delete $g->{phenotypes};
#		push(@{$res->{genes}},$g);
#		last if $n >= 4;
#		$n ++;
#	}
#	$res->{project} = $project->name;
#	$res->{id} = $id;
#	$res->{patient} = $patient->name;
# 	
#	$res->{table_id} = "table_".time."_".$id;
# 	return encode_json $res;
# }


sub getTransmission
{
	my ($sv,$mother,$father) = @_;
	
	return "?" unless $mother && $father;
		
	my $trm;
	$trm = "?";
	my $hash =  $sv->{dejavu}->{this_project};
	
	if ($mother){
		$trm = "-";
		$trm = "+"  if (exists $hash->{$mother->id});
	}
	my $trf;
	$trf = "?";
	if ($father){
		$trf = "-";
		$trf = "+"  if (exists $hash->{$father->id});
	}
	
	if ($trf eq "?" && $trm eq "?") {
		return "?";
	}
	if ($trf eq "?" && $trm eq "m") {
		return "m?";
	}
	if ($trf eq "+" && $trm eq "?") {
		return "f?";
	}
	if ($trf eq "+" && $trm eq "+") {
		return "fm";
	}
	if ($trf eq "+" && $trm eq "-") {
		return "f";
	}
	if ($trf eq "-" && $trm eq "+") {
		return "m";
	}
	if ($trf eq "-" && $trm eq "-") {
		return "d";
	}
	return "?";
}

 
# sub getLinkedCNV{
#	
#	my ($event_id) = @_;
#	my $out = ";";
#	my $djv=0;
#	
#	my ($chr1, $bp1, $chr2, $bp2) = split(/_/, $event_id);
#	
#	$chr1 = "X" if ($chr1 == 23);
#	$chr1 = "Y" if ($chr1 == 24);
#	$chr2 = "X" if ($chr2 == 23);
#	$chr2 = "Y" if ($chr2 == 24);
#	
#	# on recherche dans le dejavu des  CNV chevauchants avec bp1 et bp2
#	foreach my $t (keys %$htree_dejavuCNVProject)
#	{
#		my $htab_id;
#		
#		if ( defined ($htree_dejavuCNVProject->{$t}->{$chr1}))
#		{
#			my $tab_id1 = $htree_dejavuCNVProject->{$t}->{$chr1}->fetch($bp1-5000,$bp1+5000);
#			foreach my $cnv_id (@$tab_id1)
#			{
#				$htab_id->{$cnv_id}=1;
#			}
#		}
#		
#		if ( defined ($htree_dejavuCNVProject->{$t}->{$chr2}))
#		{
#			my $tab_id2 = $htree_dejavuCNVProject->{$t}->{$chr2}->fetch($bp2-5000,$bp2+5000);
#			foreach my $cnv_id (@$tab_id2)
#			{
#				$htab_id->{$cnv_id}=1;
#			}
#		}
#		
#		foreach my $cnv_id ( keys %{$htab_id})
#		{
#			$djv = getDejavuCNV($cnv_id);
#			next if ($djv > 20);
#			
#			$hTransLoc->{$event_id}->{'CNV'} .= $cnv_id.";"  if $hCNV->{$patientname}->{$cnv_id} == 1;
#			$hTransLoc->{$event_id}->{'CNV_mother'} .= $cnv_id.";"  if $hCNV->{$mothername}->{$cnv_id} == 1;
#			$hTransLoc->{$event_id}->{'CNV_father'} .= $cnv_id.";"  if $hCNV->{$fathername}->{$cnv_id} == 1;
#			
#		}
#	}
#}
#
#sub get_local_CNV {
#my ($sv,$chr_name,$start,$type,$pid,$mid,$fid) = @_;
#return;
#my $dir = $project->getCacheDir(). "/CNV/";
#my $nodejavu = GenBoNoSqlDejaVuCNV->new( dir => $dir, mode => "r" );
#
#my $chr = $patient->project->getChromosome($chr_name);
#my $l = abs ($start- $chr->length);
#my $a = $start;
#my $b =  $chr->length;
#warn "coucou ".$chr_name;
#my $del1 = $nodejavu->get_cnv_project($type,$chr_name,$a,$b,30);
#$a = 1;
#$b =  $start;
#
#my $del2 = $nodejavu->get_cnv_project($type,$chr_name,$a,$b,30);
#
#my @objs;
#push(@objs,@$del1);
#push(@objs,@$del2);
#
#return  unless @objs;
#
#		foreach my $cnv (@objs)
#		{
#			next if $cnv->{dejavu}->{nb_patients} > 30;
#			$sv->{'CNV'}.=$cnv->{ID}.";" if exists $cnv->{PATIENTS}->{$pid};
#			$sv->{'CNV_mother'}.=$cnv->{ID}.";" if exists $cnv->{PATIENTS}->{$mid};
#			$sv->{'CNV_father'}.=$cnv->{ID}.";" if exists $cnv->{PATIENTS}->{$fid};
#		}
#		
#return $objs[-1];
#}
#
#
#sub getDejavuCNV
#{
#	my ($event_id) = @_;
#	my $nbproject=0;
#	my $nbpatient=0;
#	my $scorecaller_evt=0;
#	my $list_of_other_patient;
#	
#	my ($t,$c,$start1,$end1) = split(/_/,$event_id);
#
#	my $no = $project->dejavuCNV();
#	my $hashdv = $no->get_cnv($c,$start1,$end1,$t,$dejavu,90);
#
#	foreach my $project_name (keys %{$hashdv})
#	{
#			$nbproject++;
#			foreach my $pname (keys %{$hashdv->{$project_name}})
#			{
#					$nbpatient++;
#			}
#	}
#
#	return $nbpatient;
#}
