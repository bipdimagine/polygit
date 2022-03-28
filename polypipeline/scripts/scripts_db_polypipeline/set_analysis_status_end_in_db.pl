#!/usr/bin/perl

use strict;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";

use Getopt::Long;
use GBuffer;
use colored; 

my $analysis_id =0;




GetOptions(
	'analysis_id=s' => \$analysis_id,
	);


	
change_step_status($analysis_id);
	
sub change_step_status {
	my ($analysis_id) = @_ ;
	my $buffer = GBuffer -> new();
	my $dbh = $buffer->dbh();
	$dbh->{AutoCommit} = 0;
	$dbh->do("use Polypipeline;");
	
	$ENV{'DATABASE'} = "";
	my $sql = qq{SELECT distinct status from Polypipeline.`Analysis_Steps` where analysis_id="$analysis_id";};
	
	my $sth = $dbh->prepare($sql);
	$sth->execute() ;
	my $res = $sth->fetchall_hashref("status");
	warn Dumper $res ;
	#tri numérique décroissant
	my @lstep_status = sort({$b <=> $a} keys %$res );
	#tri lexical décroissant
#	my @lstep_status = sort({$b cmp $a} keys %$res );
	#statut le plus élevé = 1er élément du tableau
#	warn $lstep_status[0] ;
	#correspondance avec statut final de l'analyse
	my $analysis_final_status ;
	if ( $lstep_status[0] == 1)  {
		$analysis_final_status = 1;
		warn "statut final : ".colored::stabilo("green","COMPLETE") ;
   	}
   	elsif ( $lstep_status[0] == 2 || $lstep_status[0] == 3)  {
		$analysis_final_status = 2;
		warn "statut final : ".colored::stabilo("yellow","in progress" );
   	}
   	elsif ( $lstep_status[0] == 4 )  {
		$analysis_final_status = 4;
		warn colored::stabilo("red","statut final : ERROR") ;
   	}
   	warn $analysis_final_status ;

	
	

	my $sql2= qq{UPDATE Polypipeline.`Analysis` SET status ="$analysis_final_status", date=NOW() WHERE analysis_id="$analysis_id";};

	$dbh->do($sql2);
	$dbh->commit;
	$dbh->disconnect;
}


exit(0);






	

	