#!/usr/bin/perl
use FindBin qw($Bin);
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
#use lib "$Bin/../../../../../lib/GenBoDB";
#use lib "$Bin/../../../../../lib/GenBoDB/writeDB";

use Bio::SearchIO;
use strict;
use GBuffer;
use Data::Dumper;
use Getopt::Long;
use Carp;
use JSON::XS;
use Getopt::Long; 
use decode_prediction_matrix;
use Storable qw/freeze thaw nfreeze nstore_fd nstore/;
use Set::IntSpan::Fast::XS ;
use Sys::Hostname;
use Storable qw/retrieve/;
use decode_prediction_matrix;
use POSIX qw(ceil);


 	
use Parallel::ForkManager;
use Digest::MD5 qw(md5 md5_hex md5_base64);
require("$Bin/ensembl_buffer.pm");
my $sqliteDir =  "/tmp/lmdb/annotation";

my $fork =5;
my $pm = new Parallel::ForkManager($fork);
my @aas = qw(A C D E F G H I K L M N P Q R S T V W Y);
my $pack="s".(scalar(@aas)*2);
my $version;

GetOptions(
	'version=s' => \$version,
);
die("hep version ?") unless $version;

my $factor = [];
my $description = [];
foreach my $a (@aas){
	push(@$factor,"0.001");
	push(@$description,"polyphen!".$a);
}
foreach my $a (@aas){
	push(@$factor,"0.001");
	push(@$description,"sift!".$a);
}
my $chromosomes = [1..22,'X','Y','MT'];
#my $chromosomes = [1];

#my $no = GenBoNoSql->new(dir=>$sqliteDir,mode=>"c");
$pm->run_on_finish(
		sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hash) = @_;
			unless (defined($hash) or $exit_code > 0) {
				print qq|No message received from child process $exit_code $pid!\n|;
				return;
			}
			
#			foreach my $k (keys %$hash){
#				$no->put("prediction_matrix",$k,$hash->{$k});
#			}
			#$no->close();
			}
	);



foreach my $chr (@$chromosomes){
	
	my $pid = $pm->start and next;
	my $resp = save_matrix($chr);
	$pm->finish(0,$resp);	
}
warn "waiting";
$pm->wait_all_children;

sub save_matrix {
	my ($chr) = @_;
my $nop = GenBoNoSqlRocksAnnotation->new(dir=>"/data-isilon/public-data/repository/HG19/prediction_matrix/$version/rocks/",mode=>"w",name=>$chr);	
my $buffer = ensembl_buffer->new();
$buffer->getConfig->{ensembl}->{HG38} = 75;
my $chromosome = $buffer->getSliceAdaptor->fetch_by_region('chromosome', $chr);

# Récupérer les gènes sur le chromosome 21
my @genes = @{ $chromosome->get_all_Genes };
my $nb = 0;	
# Parcourir les gènes et récupérer les protéines
foreach my $gene (@genes) {
    my $gene_id = $gene->stable_id;
    my @transcripts = @{ $gene->get_all_Transcripts };
	my $matrix_adaptor = $buffer->getMatrixAdaptor;
    foreach my $transcript (@transcripts) {
        my $protein = $transcript->translation;
        if ($protein) {
        	$nb ++;
        	warn $nb if $nb%200 ==0;
        	
        	my $all_results = [];
        	my $protein_id = $protein->stable_id;
        	my $genbo_id = $protein_id."_".$chr;
        	my $md5 =  md5_hex($transcript->translate()->seq());
        	my $test = $nop->get_raw($md5);
        	next if $test;
        
			my $matrix_polyphen = $matrix_adaptor->fetch_polyphen_predictions_by_translation_md5($md5);	
        	my $matrix_sift = $matrix_adaptor->fetch_sift_predictions_by_translation_md5($md5);
        	unless ($matrix_polyphen){
				 	next;
			}
			unless ($matrix_sift){
				 	next;
			}
		#	warn $protein_id;
        #	next;
			$matrix_polyphen->expand_matrix ;
			save_values($transcript->translate()->seq(),$matrix_polyphen->{matrix},$nop,$protein_id,$all_results);
			$matrix_sift->expand_matrix ;
			save_values($transcript->translate()->seq(),$matrix_sift->{matrix},$nop,$protein_id,$all_results);
			 for (my $i=0;$i<@$all_results;$i++){
				$nop->put_batch_raw($md5."!".($i+1),$all_results->[$i]);
			}
			$nop->put_batch_raw($md5,$genbo_id."!hg19");
			$nop->put_batch_raw($genbo_id."!hg19",$md5);
			warn $genbo_id."!hg19";
        }
       
    }
}
 $nop->rocks->write($nop->batch);
 $nop->rocks->compact_range();
 $nop->close;
 return {1=>1};
}


sub save_values {
	my ($seq,$matrix,$no,$id,$all_results) =@_;
	
	for (my $i=0;$i<length($seq);$i++){
		 my $values =[];
		foreach my $aa (@aas){
	  	 	my ($v) = decode_prediction_matrix::compact_value_from_matrix($matrix,$i+1, $aa );
	  	 	my ($p,$score) = decode_prediction_matrix::prediction_from_matrix($matrix,$i+1, $aa );
	  		$score = int($score*1000);
	  	 	push(@$values,$score);
		}
		$all_results->[$i].= pack("s".scalar(@aas),@$values);
	}
}
