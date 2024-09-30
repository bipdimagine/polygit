#!/usr/bin/perl

use FindBin qw($Bin);
use strict;

use lib $Bin;
use strict;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";



use Data::Dumper;
use Getopt::Long;
use GBuffer;
use Term::ANSIColor;
use colored;

my $project_name ;
my $SNP_id_list ;
my $filename ;

GetOptions(
	'project=s' => \$project_name,
	'file=s' => \$filename,
	);
	


my $buffer = GBuffer -> new();
my $project = $buffer->newProjectCache( -name => $project_name );
my $project_id= $project->id() ;

my $hPatbarrecode ;
my $hBarrecodePat ;
my @listbarrecodeGenotypage ;
my @list = ('21_47703649_G_A', '3_183906515_T_C','17_47000251_C_T', '16_2821573_C_T', '1_43124859_C_T', '17_5284770_G_A','21_47685939_A_G', '12_129293346_C_T','19_46857286_G_A' ,'3_49365269_C_T','6_33258443_G_A','19_41117869_A_G',  '3_73111809_A_G',  '3_52236762_G_A', '2_105654716_C_T'  );



my $pos =0;
my $id_bc;
my $limit = 10;
foreach my $var_id (@list) {  
	#eval {
	my $v = $project->_newVariant($var_id);
	#}
	foreach my $patient (@{$project->getPatients()}) {
	my $depth = $patient->depth($v->getChromosome()->id(),$v->start,$v->start+1);
	if ($depth->[0] > $limit){
		push(@{$id_bc->{$patient->name}->{coverage}},"c");
	}
	else {
		push(@{$id_bc->{$patient->name}->{coverage}},"mc");
	}
	}
#	my $ps = $v->getPatients();
#	foreach my $p (@{$ps}){
#		warn $p->name;
#	}
	my $chr = $v->getChromosome();
	#next;
	unless (defined $v->vector_id){
		foreach my $patient (@{$project->getPatients()}) {
			push(@{$id_bc->{$patient->name}->{code}},1);
		}
		next;
	}
	
	foreach my $patient (@{$project->getPatients()}) {
		if ($patient->getVectorHe($chr)->contains($v->vector_id)){
			push(@{$id_bc->{$patient->name}->{code}},3);
		}
		elsif ($patient->getVectorHo($chr)->contains($v->vector_id)){
			push(@{$id_bc->{$patient->name}->{code}},2);
		}
		else {
			my $depth= $patient->depth($v->getChromosome()->id(),$v->start,$v->start+1);
			if ($depth > 10){
				push(@{$id_bc->{$patient->name}->{code}},1) ;
			}
			else {
					push(@{$id_bc->{$patient->name}->{code}},4) ;
			}
		#	warn $depth->[0];
			
		}
		
	}
}

my $dbh = $buffer->dbh();
$dbh->{AutoCommit} = 0;
$dbh->do("use PolyprojectNGS;");

foreach my $patient (@{$project->getPatients()}) {
	my $identity_vigilance= $patient->identity_vigilance();
	my $iv_vcf=join("",@{$id_bc->{$patient->name}->{code}});
	my $iv_vcf_cov=join(" ",@{$id_bc->{$patient->name}->{coverage}});
	warn "barre code final -".$patient->name()." : ".$iv_vcf."---couvertures : ".$iv_vcf_cov; 
#	
	
#	#comparaison avec bcgenotypage et ajout d'une lettre de code couleur d'interprétation 
#	my $identity_vigilance= $patient->identity_vigilance();
#	###############################
		warn $identity_vigilance." --> ". $iv_vcf;
		if ($identity_vigilance ne "" && $iv_vcf ne "") {
			#die();
		warn "pas vide : ".$identity_vigilance."__".$iv_vcf ;
			
		##### comparaison : retour ok ou pas #####
		#my @bc = () ;
		my @bc = split("",$iv_vcf );
		warn @bc[0];

		#recherche du taux d'erreur le plus bas en comparant le barre code de NGS avec chaque barre code de génotypage du run
		#print "\n\n*****************************************************************\n".$patient->name()." : Recherche à l'aveugle du barrecode de genotypage le plus proche\n*****************************************************************\n";
		my $hbc_erreur ;
		foreach my $test_barre_code (@listbarrecodeGenotypage){
			my $ccbis=0 ; #compteur nb erreurs
			my $cbis=1; #compteur parallele chaine
			my @codebarrelist = split(//, $test_barre_code);
			#warn $test_barre_code;
			my @num_color ;
			my $global_color ="vert"; #vert par défaut, orange si il y a un 4 parmi les 15 valeurs
				foreach my $j (@codebarrelist){
				if ($j==4){
					push(@num_color, "orange");
					$global_color= "orange";
				}
				else{push(@num_color,"vert");}
			}

			
			foreach my $i (@bc){
				
				if ($i ne $codebarrelist[$cbis-1] ){ 
							
								if ($codebarrelist[$cbis-1]== 4){
									}
									else {$ccbis ++;}
				}
				$cbis ++;
				
			}
#			print "barrecode cache -".$patient->name()." : ".join("",@bc). "--resultat genotypage : ".$test_barre_code." : ".$ccbis. "/15 erreurs\n";

			$hbc_erreur->{$test_barre_code}=$ccbis;
		}

#		warn Dumper $hbc_erreur;
					
	my $h_percent ; 
	my $ctoplist =1 ;
		foreach my $k ( sort { $hbc_erreur->{$a} <=> $hbc_erreur->{$b} } keys % $hbc_erreur){
			my $percent = $hbc_erreur->{$k}/15;
			$h_percent->{$ctoplist}=$percent ; 

			if ($ctoplist ==1  ){
#			print $patient->name()."-".$k." avec : ".$hbc_erreur->{$k}." erreurs, soit : ".$percent." %, correspondant au barre code genotype de : ".$hBarrecodePat->{$k}."\n";
			if ($patient->name() eq $hBarrecodePat->{$k} && $hbc_erreur->{$k} <2) { 
				print $patient->name()." OK-vert\n";
				push(@bc,'v');
			
			}
			elsif ($patient->name() eq $hBarrecodePat->{$k} && $hbc_erreur->{$k} >1 && $hbc_erreur->{$k} <5) { 
				print $patient->name()." OK-orange\n";
				push(@bc,'w');
			}
			elsif ($patient->name() eq $hBarrecodePat->{$k} && $hbc_erreur->{$k} >4 ) { 
				print $patient->name()." NOK-rouge\n";
				push(@bc,'e');
			}
			elsif ($patient->name() ne $hBarrecodePat->{$k}  ) { 
				print $patient->name()." NOK-rouge\n";
				push(@bc,'e');
			}
			
			}
			$ctoplist ++;
			
		}
	
#warn Dumper $h_percent ;
##si % erreur du 2ème barre code est significativement différent du % d'erreur du 1er barre code trouvé, alors ok et afficher le résultat ; sinon, pas ok
##faire une table de hash [{1=>err} {2=>err {3=>err}...}] ou travailler sur la liste

	##systeme de couleurs vert orange rouge
	
 	$iv_vcf=join("",@bc);
		warn $iv_vcf ;
		
		
		
		}

		#else {warn "pas d'identitovigilance sur ce patient"; die();}
	else {warn "pas d'identitovigilance sur ce patient"; }
		####################################################
	#mettre des côtes cat $iv_vcf est une chaine de caractères
	warn "insert";
	my $patient_name = $patient->name;
	
	my $sql = qq{UPDATE  patient set identity_vigilance_vcf="$iv_vcf" where project_id="$project_id" and  name= "$patient_name" ;};
	my $sth = $dbh->prepare($sql);
	warn $sql;
	$sth->execute();
	
	
}
$dbh->commit();
