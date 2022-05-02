#!/usr/bin/perl
use FindBin qw($Bin);
use lib $Bin;
use lib "$Bin/../../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../packages/";
#use lib "/bip-d/perl/ensembl64/ensembl/modules";
use Bio::SearchIO;
use strict;
use Getopt::Long;
use Carp;
use JSON::XS;
use Getopt::Long; 
use List::MoreUtils qw(uniq);
#use ensembl_buffer;
use Storable qw/freeze thaw nfreeze nstore_fd nstore/;
use Set::IntSpan::Fast::XS ;
use Sys::Hostname;
use Storable qw/retrieve store/;
use  Set::IntSpan::Island;
use Parallel::ForkManager;
use Array::IntSpan;
use Data::Dumper;
use Array::IntSpan;
use GenBoNoSql;
use GenBoNoSqlAnnotation;
#use Bio::Ensembl::Registry;
use Storable qw(dclone);
use List::MoreUtils qw(uniq);
use Set::IntervalTree;
use Date::Tiny;
use Digest::MD5::File qw( file_md5_hex );
use Text::CSV;
use Storable;


require "$Bin/../packages/parse_gff.pm";

my $csv = Text::CSV->new({ sep_char => ',' });
require "$Bin/../packages/ensembl_buffer.pm";
#my $sqliteDir =  "/data-xfs/public-data/HG19/sqlite/75/annotation_test";
my $sqliteDir =  "/tmp/lmdb/annotation";
my $gene_code_version = "30";
#my $gff = "/data-xfs/public-data/HG19/gencode/v28/gencode.v28lift37.annotation.gff3.gz";

my $version;
my $genome_version = "HG19";
GetOptions(
	'version=s' => \$version,
);
die("add version ") unless $version; 
my $sqliteDir =  "/tmp/lmdb/$version/annotations";
system("mkdir -p $sqliteDir") unless -e $sqliteDir;
my $dir_genecode =  "/data-isilon/public-data/repository/$genome_version/annotations/gencode.v$version/tabix/";
my @files = `ls $dir_genecode/*.gz`;
chomp(@files);
my $ftype = {
	#gff =>{name => "all"},
	gff =>{name => "annotation.gff3"},
	fasta_transcripts =>{name => ".pc_transcripts.fa"},
	fasta_proteins =>{name => ".pc_translations.fa"},
	refseq =>{name => "RefSeq"},
	swissprot =>{name => "SwissProt"},
	mart =>{name => "mart"},
	
};
my $add_genes = {};
my $add_genes = {
	miss =>{file => "$Bin/missed/genes.gff3.gz"},
	regulatory => {file => "$Bin/regulatory/regulatory.gff.gz"},
	fasta_transcripts=> {file => "$Bin/missed/transcripts.fa.gz"},
	fasta_proteins=> {file => "$Bin/missed/translations.fa.gz"},
};

foreach my $ft (keys %$add_genes){
	
	die("problem with file $ft  ".Dumper $add_genes) unless -e $add_genes->{$ft}->{file};
}

my $hversion;
$hversion->{name} = "gencode";
$hversion->{version} = "$version";

my @tfiles;
my @md5;
foreach my $ft (keys %$ftype){
	my $fname = $ftype->{$ft}->{name};
	my $find;
	foreach my $f (@files){
		if ($f =~/$fname/){
			$find =1;
			warn 'OK '.$fname;
			$ftype->{$ft}->{file} = $f;
			push(@tfiles,$f);
			push(@md5,file_md5_hex( $f));
			last;
		}
	}
	die("problem with file $ft  ".Dumper @files) unless $find;
}
#concat  gff 
$ftype->{gff_all}->{file} = "$sqliteDir/all.gff3.gz";

system("cat ".$ftype->{gff}->{file}." ".$add_genes->{miss}->{file}." ".$add_genes->{regulatory}->{file}.">".$ftype->{gff_all}->{file});

$ftype->{fasta_transcripts_all}->{file} = "$sqliteDir/transcript.fa.gz";
die($ftype->{fasta_transcripts}->{file}) unless -e $ftype->{fasta_transcripts}->{file};
system("cat ".$ftype->{fasta_transcripts}->{file}." ".$add_genes->{fasta_transcripts}->{file}.">".$ftype->{fasta_transcripts_all}->{file} );
warn "cat ".$ftype->{fasta_transcripts}->{file}." ".$add_genes->{fasta_transcripts}->{file}.">".$ftype->{fasta_transcripts_all}->{file};

$ftype->{fasta_proteins_all}->{file} = "$sqliteDir/translations.fa.gz";
die() unless -e $ftype->{fasta_proteins}->{file};
system("cat ".$ftype->{fasta_proteins}->{file}." ".$add_genes->{fasta_proteins}->{file}.">".$ftype->{fasta_proteins_all}->{file} );
warn "cat ".$ftype->{fasta_proteins}->{file}." ".$add_genes->{fasta_proteins}->{file}.">".$ftype->{fasta_proteins_all}->{file};

$hversion->{file} = join(";",@tfiles);
$hversion->{md5sum} = join(";",@md5);
warn Dumper @tfiles;
my $gene_name;

my $translate_id ={};
warn 'parse gene';
my $genes = parse_gff::read_gff_genes($ftype->{gff_all}->{file},$translate_id,$gene_code_version);
warn "description";
warn $ftype->{mart}->{file};
readMetaDataMartGenes($ftype->{mart}->{file},$translate_id,$genes,"description");
my $transcripts;
my $proteins;
if (-e $sqliteDir."/translate_id.freeze"){
	$proteins = retrieve($sqliteDir."/proteins.freeze");
	$transcripts = retrieve($sqliteDir."/transcripts.freeze");
	$genes = retrieve($sqliteDir."/genes.freeze");
	$translate_id = retrieve($sqliteDir."/translate_id.freeze");
#	warn "coucou";
}
else {
	
 ($transcripts,$proteins) = parse_gff::read_gff_transcripts($ftype->{gff_all}->{file},$genes,$translate_id,$gene_code_version);

store $transcripts, $sqliteDir."/transcripts.freeze";
store $proteins, $sqliteDir."/proteins.freeze";
store $genes, $sqliteDir."/genes.freeze";
store $translate_id, $sqliteDir."/translate_id.freeze";
}
foreach my $prot_id (keys %$proteins) {
	next if $prot_id =~/ENSP/;
	warn Dumper $proteins->{$prot_id};
	warn $prot_id;
	die();
}

 warn "parse metadata +++++++ ";
 readMetaDataProteins($ftype->{swissprot}->{file},$translate_id,$transcripts,$proteins,"external_name");
 foreach my $prot_id (keys %$proteins) {
	next if $prot_id =~/ENSP/;
	warn Dumper $proteins->{$prot_id};
	warn $prot_id;
	die();
}
 my $proteins_seq = readFastaProtein($ftype->{fasta_proteins_all}->{file},$translate_id,$transcripts,$proteins);
 foreach my $prot_id (keys %$proteins) {
	next if $prot_id =~/ENSP/;
	warn Dumper $proteins->{$prot_id};
	warn $prot_id;
	die();
}
my $transripts_seq = readFasta($ftype->{fasta_transcripts_all}->{file},$translate_id,$transcripts);


 #warn Dumper $proteins->{ENSP00000498926_1};

readMetaDataTranscripts($ftype->{refseq}->{file},$translate_id,$transcripts,"external_name");


warn 'filering transcripts';


my %delete_transcripts;
	
	foreach my $tid (keys %$transcripts) {
		my $new_transcript = $transcripts->{$tid};
		my $debug;
		$debug =1 if $tid eq "ENST00000341947_6";
		warn Dumper  $new_transcript->{tag} if $debug;
		my $tl = $new_transcript->{transcript_support_level};
		warn "delete $tl" if $debug;
		next if $tl < 3 ;
		
		next if exists  $new_transcript->{tag}->{basic};
		warn "deleet basic" if $debug;
		next if exists $new_transcript->{tag}->{CCDS};
		warn "deleet CCDS" if $debug;
		warn Dumper  $new_transcript->{tag} if $debug;
		next if exists  $new_transcript->{tag}->{APPRIS};
		next if exists $new_transcript->{tag}->{MANE};
		if ($tl >=4 or $new_transcript->{biotype_ensembl} eq 'processed_transcript' or exists $new_transcript->{tag}->{NF})  {
			$delete_transcripts{$tid} ++;
			warn "deleet " if $debug;
#			
#			my $protein_id = $new_transcript->{protein_id};
#			delete $proteins->{$protein_id};
#			my @gene_ids = keys %{$new_transcript->{genes_object}};
##		warn $gene_ids[0];
#			die("ici") if scalar(@gene_ids)>1;
#			die("or here") unless exists $genes->{$gene_ids[0]}->{transcripts_object}->{$new_transcript->{id}};
#			delete $genes->{$gene_ids[0]}->{transcripts_object}->{$new_transcript->{id}};
#			my @ts = keys %{$genes->{$gene_ids[0]}->{transcripts_object}};
#			delete $genes->{$gene_ids[0]} if scalar(@ts) == 0;
#			warn $gene_ids[0] if scalar(@ts) == 0;
			next ;
		}

	} 
	foreach my $tid (keys %delete_transcripts) {
		my $protein_id = $transcripts->{$tid}->{protein_id};
		delete $proteins->{$protein_id} if $protein_id;
		my $gene_id = $transcripts->{$tid}->{gene};
		
		die($tid) unless exists $genes->{$gene_id};
		unless (exists $genes->{$gene_id}->{transcripts_object}->{$tid}){
			warn $tid;
			warn $gene_id;
			warn Dumper $genes->{$gene_id};
		}
		die() unless exists $genes->{$gene_id}->{transcripts_object}->{$tid};
		delete $genes->{$gene_id}->{transcripts_object}->{$tid};
		delete $transcripts->{$tid};
	}
foreach my $prot_id (keys %$proteins) {
		next if $proteins->{$prot_id}->{sequence};
		#warn "coucou ".$prot_id;
		die( $prot_id) if $prot_id =~/ENST/;
		#warn $prot_id;
#		warn $prot_id;
		my $tid = $proteins->{$prot_id}->{transcript};
		delete $transcripts->{$tid}->{protein_id};
		delete $transcripts->{$tid}->{protein};
		delete $transcripts->{$tid}->{proteins_object}->{$prot_id};
		delete $transcripts->{$tid}->{protein_stable_id};# = $infos->{protein_id};
		delete $transcripts->{$tid}->{protein_genbo_id} ;#= $infos->{protein_id};
		delete $proteins->{$prot_id};
	
}
 warn scalar keys %$proteins;
 die() if (scalar keys %$proteins) eq 0;
 foreach my $gid (keys %$genes){
 	my $g = $genes->{$gid};
 	#warn $gid unless  $g->{transcripts_object};
 	 delete $genes->{$gid} unless  $g->{transcripts_object};
 	  delete $genes->{$gid} unless  keys %{$g->{transcripts_object}};
 	  $genes->{$gid}->{transcripts} = [keys %{$g->{transcripts_object}}];
 }


warn "SAVE";
 my $no2 = GenBoNoSqlAnnotation->new(dir=>$sqliteDir,mode=>"c");
 warn "\t genes";
 parse_gff::save_sqlite ($no2,$genes,"gene");
 warn "\t transcripts";
 parse_gff::save_sqlite ($no2,$transcripts,"transcript",undef,1);
  warn "\t protein";
   warn "\t protein";
   
 parse_gff::save_sqlite ($no2,$proteins,"protein",undef,1);
$no2->close();


my $d = Date::Tiny->now;
$hversion->{name} = "genecode";
$hversion->{version} = "$version";
$hversion->{date} =  $d->as_string;
warn "$sqliteDir//version.json";
open(JSON,">$sqliteDir/version.json") or die();
print JSON encode_json $hversion;
close JSON;
warn "END !!!! ";

exit(0);

sub readMetaDataMartGenes{
		my ($file,$translate_id,$genes,$type) = @_;
warn "mart ".$file;
open (READ, "zcat $file | ") || die "Cannot open $file: $!.\n";
	while (my $line = <READ>){
		chomp($line);
		
		my($zid,$description,@t) = split(";",$line) ;
		next unless exists $translate_id->{$zid};
		warn $zid unless exists $translate_id->{$zid};
		foreach my $rid (@{$translate_id->{$zid}}) {
		if ($description){
		 $genes->{$rid}->{description} = $description unless $genes->{$rid}->{description};
		}
	
		}
	

	}
close (READ);
}

sub readMetaDataMartGenes_old{
		my ($file,$translate_id,$genes,$type) = @_;
warn "mart ".$file;
open (READ, "zcat $file | ") || die "Cannot open $file: $!.\n";
	while (my $line = <READ>){
		chomp($line);
		
		my($zid,$description,$pheno,$source,@t);
		if ($line =~/"/){
			if ($csv->parse($line)) {
				
				($zid,$pheno,$source,$description,@t) =  $csv->fields();
			}
		}
		else {
			($zid,$pheno,$source,$description,@t) =  split(",",$line);
		#	warn $zid;
		#	warn $description.' '.$pheno;
		#	die();
			#next;
		}
		#warn $zid;
		
		next unless exists $translate_id->{$zid};
		warn $zid unless exists $translate_id->{$zid};
		foreach my $rid (@{$translate_id->{$zid}}){
			
		
			
#		warn $rid." -->".$description;
		if ($description){
		 $genes->{$rid}->{description} = $description unless $genes->{$rid}->{description};
		}
		
		 next unless $pheno;
		 next if $pheno eq "";
		 $genes->{$rid}->{hphenotype}->{ensembl}->{lc($pheno)} ++;
		$source = "omim" if $source =~/MIM/;
		$source = lc($source);
		$genes->{$rid}->{hphenotype}->{$source}->{lc($pheno)} ++;
	 
	
		}
	

	}
close (READ);
foreach my $rid  (keys %$genes){
	next unless exists $genes->{$rid}->{hphenotype};
	foreach my $source (keys %{$genes->{$rid}->{hphenotype}}){
		$genes->{$rid}->{phenotype}->{$source} = join(';',keys %{$genes->{$rid}->{hphenotype}->{$source}});
		delete $genes->{$rid}->{hphenotype}->{$source};
#		warn $source.' '.$genes->{$rid}->{phenotype}->{$source};
	}
#	warn $genes->{$rid}->{description};
	
}
		
}




sub readMetaDataTranscripts{
		my ($file,$translate,$transcripts,$type) = @_;

open (READ, "zcat $file | ") || die "Cannot open $file: $!.\n";
	while (my $line = <READ>){
		chomp($line);
		my($rid,$a,@t) = split(" ",$line);
		next unless exists $translate_id->{$rid};
		
		foreach my $id (@{$translate_id->{$rid}}){
			next unless exists $transcripts->{$id}->{genbo_id};
			my $sid = parse_gff::getId($a);
			$transcripts->{$id}->{$type} = $sid	;#unless exists $transcripts->{$id}->{$type};
			$transcripts->{$id}->{refseq}->{$sid} ++; 
		
		} 
		}
		
close (READ);
		
}
sub readMetaDataGenes{
		my ($file,$translate,$genes,$type) = @_;

open (READ, "zcat $file | ") || die "Cannot open $file: $!.\n";
	while (my $line = <READ>){
		chomp($line);
		my($rid,$a,@t) = split(",",$line);
		next unless exists $translate_id->{$rid};
		foreach my $id (@{$translate_id->{$rid}}){
		$genes->{$id}->{$type} = $a; 
		}
	}
close (READ);
		
}
sub readMetaDataProteins{
		my ($file,$translate,$transcripts,$proteins,$type) = @_;

open (READ, "zcat $file | ") || die "Cannot open $file: $!.\n";
	while (my $line = <READ>){
		chomp($line);
		my($rid,$a,@t) = split(" ",$line);
		
		
		next unless exists $translate_id->{$rid};
		foreach my $id (@{$translate_id->{$rid}}){
		next unless exists $translate_id->{$id};
		next unless exists $translate_id->{$id}->{genbo_id};
		next unless $transcripts->{$id}->{protein_id};
		my $pid = $transcripts->{$id}->{protein_id};
		$proteins->{$pid}->{external_name} = $a;
		}
	}
close (READ);
		
}

sub readFastaProtein {
	my ($file,$translate,$tr,$prot) = @_;
	warn $file;
	my %notfound;
	open (READ, "zcat $file | ") || die "Cannot open $file: $!.\n";
	my %seqs = ();
	my $header = '';
	my $ids;
	while (my $line = <READ>){
    	chomp $line;
    if($line =~ /^>(.+)/){
    		$ids =[];
            $header = $1;
            my @t = split(/\|/,$header);
           
            my $pary ;
        #    $pary ==1 if ($header =~ /PAR_Y/);
            my $id1 = $t[0];
            my $tt;
            my $rid;
            ($rid,$tt) = split("_",$t[0]);
#            warn $rid;
            $rid =~s/_.*//;
           
            $ids = $translate_id->{$rid};
            #warn $header;
        	}
        	else {
        		foreach my $id (@$ids){
        			
        			my $protid;
        			if (exists $tr->{$id}){
        				$protid = $tr->{$id}->{protein_id};
        			}
        			elsif (exists $prot->{$id}) {
        				$protid = $id;
        			}
        			else {
        				$notfound{$id} ++ unless $protid;
        				next;
        			}
        			die($protid) if $protid =~/ENST/;
#        			next unless 
				
            		$prot->{$protid}->{sequence} .= $line;
        		}
        	}
}
close (READ);
}

sub readFasta {
	my ($file,$translate,$obj) = @_;
	my %notfound;
	open (READ, "zcat $file | ") || die "Cannot open $file: $!.\n";

	my %seqs = ();
	my $header = '';
	my $ids;
	while (my $line = <READ>){
    	chomp $line;
    if($line =~ /^>(.+)/){
    		$ids =[];
            $header = $1;
            my @t = split(/\|/,$header);
            my $pary ;
        #    $pary ==1 if ($header =~ /PAR_Y/);
            my $id1 = $t[0];
            my $tt;
            my $rid;
            ($rid,$tt) = split("_",$t[0]);
            $rid =~s/_.*//;
           
            $ids = $translate_id->{$rid};
            
            #warn $header;
        	}else{
        		foreach my $id (@$ids){
        		 	$notfound{$id} ++ unless exists $obj->{$id}->{genbo_id};
        		 	warn $id if $id eq "ENST00000291688_21";
        			next unless exists $obj->{$id}->{genbo_id};
#        			next unless 
            		$obj->{$id}->{sequence} .= $line;
            		warn $obj->{$id}->{sequence}  if $id eq "ENST00000291688_21";
        		}
        	}
}
close (READ);
 
foreach my $transcript (values %$obj){
	next if $transcript->{genbo_type} ne "transcript";
	next unless exists $transcript->{orf_start};
	my $cstart =	$transcript->{orf_start};
	my $l3 = abs ($transcript->{orf_end}-$transcript->{orf_start})+1;
	$transcript->{coding_sequence} = substr($transcript->{sequence},$cstart-1,$l3)	;
}

}
	


