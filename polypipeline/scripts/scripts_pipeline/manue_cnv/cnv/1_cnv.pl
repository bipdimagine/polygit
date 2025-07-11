#!/usr/bin/perl
use FindBin qw($Bin);
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
use Clone qw(clone);
use Parallel::ForkManager;
use strict;
use Text::CSV qw( csv );
use lib "$Bin/../../../../../GenBo/lib/obj-nodb/";
use Getopt::Long;
use lib "$Bin";
require  "$Bin/../SVParser.pm";
require  "$Bin/../parser/parse_pbsv.pm";
require  "$Bin/../parser/parse_hificnv.pm";
require  "$Bin/../parser/parse_sniffles2.pm";
require  "$Bin/../parser/parse_wisecondor.pm";
use lib "$Bin/../../dejavu/utility/";
use liftOverRegions;
use GBuffer;
use Text::CSV;



###################################################################
# Cherche tous les variants structuraux d'un projet pour chaque patient et tous callers confondus.
# Reconstruit les CNV fragmentes (même caller et bornes a moins de 10% de la longeur)
# Construit une table de hash par patient et freeze ces tables  : nom du fihier = patient.allSV
# Pour le dejavu construit egalement une table de hash qui conserve tous les CNV du projet 
# et garde pour chaque CNV l'info du patient et du caller qui l'a détecté.
# Pour chaque CNV la localisation chr_start_end donne accès à la liste des gènes compris dans cet interval
####################################################################

my $fork = 5;
my $cgi = new CGI;


my $limit;
my $projectname;
my $patient_name;
my $fork =1 ;
GetOptions(
	'project=s' => \$projectname,
	'patient=s' => \$patient_name,
	'fork=s' => \$fork,
);

 #= $cgi->param('projectname');
#my $patient_name = $cgi->param('patient');
#$fork = $cgi->param('fork');

#$fork=1;

# pour récupérer les objets project et patient
my $buffer = GBuffer->new();

my $project = $buffer->newProjectCache( -name => $projectname);

my $caller_type_flag = {
	"caller_sr" => 1,
	"caller_depth" => 2,
	"caller_coverage" => 4,
};
my $value_caller_type_flag = {
	1=>"caller_sr" ,
	2=>"caller_depth",
	4=>"caller_coverage" ,
};
my $type_by_caller_name = {
	"wisecondor" => "caller_coverage",
	"manta" => "caller_sr",
	"pbsv" => "caller_sr",
	"dragen-sv" =>"caller_sr",
	"Sniffles2" =>"caller_sr",
	"Spectre" =>"caller_sr",
	"hificnv" =>"caller_depth",
	"canvas" =>"caller_depth",
	"dragen-cnv" =>"caller_depth",
	"cnvnator" =>"caller_depth",
};


my $listPatients = $project->getPatients();
if ($patient_name) {
	$listPatients = [$project->getPatient($patient_name)];
}



# on gardera que les DUP DEL de taille supérieur à 1kb
#my $length = 1000;


# pour le fichier annoSV
my $hannot;
my $fdannot;

# pour lire les vcfs
my $ligne;
my $fd;
my $compteur=0;
	
# pour stocker les SV

my $hPat_CNV; # pour les id_globaux

# pour le dejavu
my $hdejavu;



#########################################
#  Pour les duplications segmentales
#########################################
my $cytodir = $project->get_cytology_directory;	
my $fichier_dupseg = $cytodir."/segmental_duplication.bed";
my $fdds;

#  Lecture du fichier dupseg
	my $hdupseg;
	open($fdds,$fichier_dupseg) or die("open: $!");
	
	# lecture ligne par ligne
	while( defined( $ligne = <$fdds> ) )
	{
		next unless ($ligne =~ m/^chr/ );
		my @champs = split(/\t/,$ligne);
		
		my $chr = $champs[0];
		my $start = $champs[1];
		my $end = $champs[2];
		
		my $idLocation = $chr.":".$start."_".$end;
		# version avec IntervalTree
		if (!exists $hdupseg->{$chr})
		{
		 	my $DupSegTree = Set::IntervalTree->new;
		 	$hdupseg->{$chr}=$DupSegTree;
		}		
		$hdupseg->{$chr}->insert($idLocation,$start,$end);		
	}	
	
	
############
#  lecture des vcf 
############


my $SVtype;
my $SVend ;
my $SVlength;
my $SVchr;
my $SVdeb;


# pour la sauvegarde des fichiers de sortie
#my $variationsDir = "/data-xfs/Manue/Test_SV/".$projectname."/newVersion/";
#my $variationsDir= $project->getVariationsDir();
my $variationsDir = $project->getCNVDir();

my $pm = new Parallel::ForkManager($fork);
$project->getChromosomes();		
my $job_id = time;
my $hjobs ;
my $chrs;
foreach my $c (@{$project->getChromosomes}){
	next if $c->name eq "MT";
	$chrs->{$c->name} ++;
	$chrs->{$c->ucsc_name} ++;
}
my $type_caller = {
	"caller_sr" => 1,
	"caller_depth" => 2,
	"caller_coverage" => 4,
};

#"manta" or $caller  eq "pbsv" or $caller  eq "dragen-sv"
#$caller  eq "hificnv" or $caller  eq "canvas" or $caller  eq "dragen-cnv"
my $type_by_caller = {
	"wisecondor" => 4,
	"manta" => 1,
	"pbsv" => 1,
	"dragen-sv" =>1,
	"hificnv" =>2,
	"canvas" =>2,
	"dragen-cnv" =>2,
};

my $type_by_caller = {
	"wisecondor" => 4,
	"manta" => 1,
	"pbsv" => 1,
	"dragen-sv" =>1,
	"Sniffles2" =>1,
	"Spectre" =>1,
	"hificnv" =>2,
	"canvas" =>2,
	"dragen-cnv" =>2,
	"cnvnator" =>2,
};

my $cnv_callers = {
    "wisecondor"          => 1 << 0,  # 2^0 = 1
    "canvas"        => 1 << 1,  # 2^1 = 2
    "manta"        => 1 << 2,  # 2^2 = 4
    "pbsv"        => 1 << 3,  
    "dragen-sv"        => 1 << 4,  
    "hificnv"        => 1 << 5,  
    "dragen-cnv"        => 1 << 6, 
     "cnvnator"        => 1 << 7, 
      "Sniffles2"        => 1 << 8, 
      "Spectre"        => 1 << 9, 
};



$pm->run_on_finish(
    	sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    		my $j = $data->{job};
    		delete $hjobs->{$j};
    		
    }
    );
  $project->buffer->hash_genes_omim_morbid();
  $project->setPhenotypes();
    
  $project->disconnect();
  $buffer->dbh_deconnect();
  
  
  foreach my $patobj (@$listPatients)
{
	my $listCallers = $patobj->callingSVMethods();
	 $patobj->isGenome;
}
  $project->disconnect();
  $buffer->dbh_deconnect();
my $total_cnvs;
foreach my $patobj (@$listPatients)
{
	warn "------------";
	warn $patobj->name;
	warn "---------------";
	my $listCallers = $patobj->callingSVMethods();
	#or $patobj->name ne "dl-2-E-sg-A";
	#next unless $patobj->name =~ /short/;
	next unless $patobj->isGenome;
	$job_id ++;
	$hjobs->{$job_id} ++;
	#my $pid = $pm->start and next;
	$listCallers = $patobj->callingSVMethods();
	$project->buffer->dbh_reconnect();
	my $patname= $patobj->name();
	my $file_out = $variationsDir.$patname.".allCNV";
	unlink  $file_out if -e $file_out;
	#next if (-e $file_out);
	my $patient = $patobj;
	my $hPat_elementaryCNV;  # tous les CNV du vcf  (élémentaires = avant regroupement)
	warn Dumper @$listCallers;
	$hPat_CNV ={};
	foreach my $caller (@$listCallers)
	{
		confess() unless exists $type_by_caller->{$caller};
		warn "+++".$caller;
		
		my $dir = $project->getVariationsDir($caller);
	
		#    Lecture de la première ligne du fichier Annot 
		#    pour recuperer le format

		
		my $fichierPatient = $patobj->getSVFile($caller);
		my $res;	
		if ($type_by_caller->{$caller}== $type_caller->{caller_coverage})
		{
	
				my $hash = parse_wisecondor::parse_cnv($patient);
				warn "wisecondor".$patient->name;
				$hPat_CNV->{'caller_coverage'} = SVParser::gatherCNV_from_samecaller($patname,$hash);
		}
		if ($type_by_caller->{$caller} == $type_caller->{caller_sr}){
			if ($caller eq "pbsv"){
				$hPat_CNV->{'caller_sr'}  = parse_pbsv::parse_cnv($patient,$caller);
			}
			elsif (lc($caller) eq "sniffles2"){
				warn "coucou";
				$hPat_CNV->{'caller_sr'}  = parse_sniffles2::parse_cnv($patient,$caller);
			}
			else {
				$hPat_CNV->{'caller_sr'}  = SVParser::parse_vcf($patient,$caller);
			}
		}
		elsif ($type_by_caller->{$caller}== $type_caller->{caller_depth}) {
			my $hash;
			if ($caller eq "hificnv"){
				$hash  = parse_hificnv::parse_cnv($patient,$caller);
			}
			else {
				 $hash  = SVParser::parse_vcf($patient,$caller);
			}
			$hPat_CNV->{'caller_depth'} = SVParser::gatherCNV_from_samecaller($patname,$hash);
		}
	}
	my $cnvs = gather_identical_CNV($hPat_CNV,$patient);
	push(@$total_cnvs,@$cnvs);

} # fin boucle sur les patients
#$pm->wait_all_children();
liftover_and_dejavu_parquet($total_cnvs);
exit(0);



warn "END PATIENT !!!";
confess() if scalar(keys %{$hjobs});
$project->buffer->dbh_reconnect();
exit(0);

sub liftover_and_dejavu_parquet {
	my ($all_hash) = @_;
		my $lift = liftOverRegions->new(project=>$project,version=>$project->lift_genome_version);
foreach my $hcnv1 (@{$all_hash})
				{	
#					my $hcnv;
#					my $id = $hcnv1->{ID};
#					my ($type,$chr,$start,$end) = split("_",$id);	
#					
#					$hcnv->{chromosome} = $project->getChromosome($chr)->ucsc_name;
#					$hcnv->{type} = $type;
#					$hcnv->{id} = $id;
#					my ($t,$c,$s,$end) = split("_",$id);
#					my $debug;
#					$debug = 1  if  $start == 34955437;
#					$hcnv->{start} = $start;
#					$hcnv->{end} = $end;
#					$hcnv->{patient} = $hcnv->{patient};
#					my $callers = 0;
#					foreach my $caller (keys %{$hcnv1->{PATIENT}->{CALLERS}}){
#						$callers = $callers | $cnv_callers->{$caller};
#					}
#					
#					$hcnv->{callers} = $callers;
					$lift->add_region_id($hcnv1);
				}
warn "lift";
	my $lift =  $lift->liftOver_regions_cnv($project->name);
	warn "end";
	 save_csv($project,$lift);
	
		
}
sub gather_identical_CNV
{
	
	my ($hPat_CNV,$patient) = @_;
	my $hCNV;
	
	foreach my $caller_flag (keys %{$hPat_CNV})
	{
				foreach my $id  (keys %{$hPat_CNV->{$caller_flag}})
				{
					
					my $current_cnv = $hPat_CNV->{$caller_flag}->{$id};
					$hCNV->{$id}->{"elementary_caller_sr"} = "" unless exists $hCNV->{$id}->{"elementary_caller_sr"};
					$hCNV->{$id}->{"elementary_caller_depth"} = "" unless exists $hCNV->{$id}->{"elementary_caller_depth"};
					$hCNV->{$id}->{"elementary_caller_coverage"} = "" unless exists $hCNV->{$id}->{"elementary_caller_coverage"}; 
					$hCNV->{$id}->{'id'} = $id."!".$patient->id;
								
					$hCNV->{$id}->{'type'} = $current_cnv->{'SVTYPE'};		
					$hCNV->{$id}->{'chromosome'} = $project->getChromosome($current_cnv->{'CHROM'})->ucsc_name;
					$hCNV->{$id}->{'start'} = $current_cnv->{'START'};
					$hCNV->{$id}->{'end'} = $current_cnv->{'END'};
					$hCNV->{$id}->{'len'} = abs ($hCNV->{$id}->{'END'} - $hCNV->{$id}->{'START'}) +1;				
					# pour sauvegarder l'info propre aux differents callers
					$hCNV->{$id}->{callers} = 0 unless exists $hCNV->{$id}->{caller_flag};
					 my $vc = $cnv_callers->{$current_cnv->{CALLER}};
#					warn $current_cnv->{REAL_CALLER}." ".$vc;
					 die() unless $vc; 
					
					$hCNV->{$id}->{callers} =  $hCNV->{$id}->{callers} | $vc;
					
					$hCNV->{$id}->{caller_type_flag} = 0 unless exists $hCNV->{$id}->{caller_type_flag};
					$hCNV->{$id}->{caller_type_flag} =  $hCNV->{$id}->{caller_type_flag} | $caller_type_flag->{$caller_flag};
					$hCNV->{$id}->{sr1} = -1; 
					$hCNV->{$id}->{sr2} = -1; 
					$hCNV->{$id}->{pr1} = -1; 
					$hCNV->{$id}->{pr2} = -1;
					$hCNV->{$id}->{sr_qual} = -1;
					my $caller = $hCNV->{$id}->{CALLER};
					if ($caller_flag eq "caller_sr") {
					$hCNV->{$id}->{sr1} = $current_cnv->{'INFOS'}->{SR}->[0]; 
					$hCNV->{$id}->{sr2} = $current_cnv->{'INFOS'}->{SR}->[1];; 
					$hCNV->{$id}->{pr1} =  $current_cnv->{'INFOS'}->{PR}->[0]; ; 
					$hCNV->{$id}->{pr2} =  $current_cnv->{'INFOS'}->{PR}->[1]; 
					$hCNV->{$id}->{sr_qual} =  $current_cnv->{'QUAL'};
					}
					$hCNV->{$id}->{depth_CN} = -1; 
					$hCNV->{$id}->{depth_qual} = -1;
					$hCNV->{$id}->{"elementary_".$caller_flag} = join(";",@{$current_cnv->{'ELEMENTARY'}});
					$hCNV->{$id}->{gt} =  $current_cnv->{'GT'};
					
					if ($caller_flag eq "caller_depth")  {
						$hCNV->{$id}->{depth_CN} = $current_cnv->{'CN'};
						$hCNV->{$id}->{depth_qual} = $current_cnv->{'QUAL'};
					}
					$hCNV->{$id}->{coverage_zscore} = -1;
					$hCNV->{$id}->{coverage_ratio} = -1;
					if ($caller_flag eq "caller_coverage")  {
						$hCNV->{$id}->{coverage_ratio} = $current_cnv->{'RATIO'};
						$hCNV->{$id}->{coverage_zscore} =$current_cnv->{'INFOS'}->{zscore}; 
					}
					$hCNV->{$id}->{patient} = $patient->id;
					getScoreCallers($hCNV->{$id});
				}
					
	}
	return [values %$hCNV];
}
sub test_type {
	my($value,$flag) = @_;
	confess() unless $caller_type_flag->{$flag};
	return $value & $caller_type_flag->{$flag};
}





sub getScoreCallers {
		my ($icnv) = @_;
		my $score = 0;
		my $sc_total = 0;
				my $iscore =0;
				$icnv->{score_callers}->{nb} = 0;
				$icnv->{score_caller_coverage} =  $icnv->{score_caller_depth} = $icnv->{score_caller_sr} = 0;
				if (test_type($icnv->{caller_type_flag},"caller_coverage")) {
					$icnv->{score_callers}->{nb} ++;
					$icnv->{score_caller_coverage} = 0.75;
					$icnv->{score_caller_coverage} += 0.25 if abs($icnv->{coverage_zscore}) > 40;
				}
				if (test_type($icnv->{caller_type_flag},"caller_depth")){
					$icnv->{score_callers}->{nb} ++;
					$icnv->{score_caller_depth} = 0.75;
					$icnv->{score_caller_depth} += 0.25 if abs($icnv->{depth_qual}) > 30;
				}
				if (test_type($icnv->{caller_type_flag},"caller_sr")){
					$icnv->{score_callers}->{nb} ++;
					$icnv->{score_caller_sr} = 0.50;
					my $limit = 400;
					if ($icnv->{callers} & $cnv_callers->{pbsv}){
						$limit = 100;
					}
					elsif ($icnv->{callers} & $cnv_callers->{Sniffles2}){
						$limit = 100;
					}
					elsif ($icnv->{callers} & $cnv_callers->{Spectre}){
						$limit = 100;
					}
					elsif ($icnv->{callers} & $cnv_callers->{manta}){
						$limit = 400;
						
					}
					$icnv->{score_caller_sr} = 0.25 if $icnv->{sr2} > 10;
					$icnv->{score_caller_sr} = 0.25 if $icnv->{pr2} > 7;
					$icnv->{score_caller_sr} += 0.25 if abs($icnv->{sr_qual}) > $limit;
				}
}
sub save_csv {
	my ($project,$snps) = @_;
	my $dir_tmp = "/data-beegfs/tmp/";
	my $filename = "$dir_tmp/".$project->name.".csv";
	my $fh;		
	open( $fh, ">", $filename) or die "Impossible d'ouvrir $filename: $!";
	my $csv = Text::CSV->new({ binary => 1, eol => "\n" });
	my $col = ["project","patient","chr38","start38","end38","chr19","start19","end19","type","callers","caller_type_flag",'sr1','sr2','sr_qual','pr1','pr2',"depth_qual",'depth_CN','coverage_zscore','coverage_ratio','elementary_caller_sr','elementary_caller_depth','elementary_caller_coverage','score_caller_sr','score_caller_depth','score_caller_coverage','gt'];
	$csv->print($fh, $col); 
	$csv->print($fh, [0,0,"Z",-1,-1,"Z",-1,-1,"Z",0,0,0,0,0,0,0,0,0,0,0,"Z","Z","Z",0,0,0,"Z"]); 
	foreach my $vhh (values %$snps){
		my $chr = $project->getChromosome($vhh->{chromosome});
		my $chr0 = $chr->name;
		my $pos0 = $vhh->{start};
		my $startlift =-1;
		my $endlift =-1;
		my $chrlift ="Z";
		my $mt;
		if ($chr->name eq "MT"){
			$mt=1 if $project->getChromosome("MT")->length() == 16571 or $project->current_genome_version eq "HG38";
		}
		if (exists $vhh->{LIFT} && $vhh->{LIFT}->{MULTI}  == 1 ) {
			 if ($project->isChromosomeName($vhh->{LIFT}->{chromosome})){
				$startlift = $vhh->{LIFT}->{start};
				$endlift = $vhh->{LIFT}->{end};
				$chrlift = $project->getChromosome($vhh->{LIFT}->{chromosome})->name;
			 }
			
		}
		my $start38 = $vhh->{start};
		my $end38 = $vhh->{end};
		my $chr38 = $chr->name;
		my $start19 = $startlift ;
		my $end19 = $endlift;
		my $chr19 = $chrlift ;
		if ($project->current_genome_version eq "HG19"){
			$start19 = $vhh->{start};
			$end19 = $vhh->{end};
			$chr19 = $chr->name;
			$start38 = $startlift ;
			$end38 = $endlift ;
			$chr38 = $chrlift ;
		}
		
		if ($mt == 1 ) {
				$start19 = $start38;
				$end19 = $end38;
		}
		$vhh->{chr38}=$chr38;
		$vhh->{start38}=$start38;
		$vhh->{end38}=$end38;
		$vhh->{chr19}=$chr19;
		$vhh->{start19}=$start19;
		$vhh->{end19}=$end19;
		$vhh->{project}=$project->id;
		my $line;
		foreach my $c (@$col){
			die(Dumper $vhh." ".$c) unless exists $vhh->{$c};
			push(@$line,$vhh->{$c});
		}
			$csv->print($fh, $line);
		
	}
	
	close($fh);
	warn $filename;
	my $txt_filename =  "\'".$filename."\'";
	my $dbh = DBI->connect("dbi:ODBC:Driver=DuckDB;Database=:memory:", "", "", { RaiseError => 1 , AutoCommit => 1});
	my $parquet_file = $project->buffer->config_path("dejavu","cnv").$project->name.".".$project->id.".parquet";
	my $query = "
	COPY (
        SELECT * from read_csv_auto([$txt_filename]) order by type,chr38,end38,chr19,end19
    )
    TO '$parquet_file' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE TRUE,ROW_GROUP_SIZE 1000000);";
    warn $query;
	$dbh->do($query);
	$dbh->disconnect;
	
	
}


