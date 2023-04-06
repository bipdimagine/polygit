#!/usr/bin/perl
use strict;
use Getopt::Long;
#include the module
use Bio::DB::HTS::Faidx;
 use Data::Dumper;
#create the index object
use FindBin qw($RealBin);
use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
use GBuffer;
use Bio::DB::HTS::Tabix;
use List::MoreUtils qw{ natatime uniq};
use Parallel::ForkManager;
use Try::Catch;
use Set::IntervalTree;
use GenBoNoSqlRocks;

my $dir;
my $version;
my $fork;
GetOptions(
	'version=s'  => \$version,
	'fork=s'  => \$fork,
);
$fork = 1 unless $fork;
my $file = "transcripts.fa.gz";
my $fileprot = "proteins.fa.gz";
my $dir1 =  "/data-isilon/public-data/repository/HG19/annotations/gencode.v$version/";
my $dir_HG19 = "$dir1/fasta/HG19/";
my $dir_HG38 = "$dir1/fasta/HG38/";
my $file19 = $dir_HG19.$file;
my $file38 = $dir_HG38.$file;
my $file_gff_HG38 = "$dir_HG38/gencode.v".$version.".annotation.gff3.gz";
my $file_gff = "$dir1/tabix/gencode.v".$version."lift37.annotation.gff3.gz";
die($file_gff) unless -e $file_gff;

#open (FILE,"zgrep \"remap_status=partial\" $file_gff | ");
my $annot2 =  GenBoNoSqlRocks->new(name=>"partial_transcripts",dir=>$dir1."/fasta/",mode=>"c",is_compress=>1);

$annot2->put("date",{time1=>time});
$annot2->close();
#my $annot2 =  GenBoNoSqlLmdb->new(name=>"partial_transcripts",dir=>$dir1."/fasta/",mode=>"c",is_compress=>1);
my @all =  `zgrep \"remap_status=partial\" $file_gff | awk \'\$3 \~ \/transcript\|exon\/\'`;

chomp(@all);
my (@genes);
my (@all_exons);
foreach my $l (@all){
	my @f = split(" ",$l);
	if($f[2] =~ /transcript/){
		push(@genes,$l);
	}
	elsif ($f[2] =~ /exon/){
			push(@all_exons,$l);
	}
	else{
		die($l);
	}
}

my $max =  scalar(@genes);
my $nb =0;
my $sid= 0;
my $hrun ={};
my $pm           = new Parallel::ForkManager($fork);


#$annot2->close();

$pm->run_on_finish(
	sub {
		my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $h ) = @_;

		unless ( defined($h) or $exit_code > 0 ) {
			print qq|No message received from child process $exit_code $pid!\n|;
			die();
			return;
		}
		#
		#warn Dumper  $h->{data};
		
		delete $hrun->{$h->{run_id}};
		my $annot2 =  GenBoNoSqlRocks->new(name=>"partial_transcripts",dir=>$dir1."/fasta/",mode=>"w",is_compress=>1);
		foreach my $s ( @{ $h->{data} } ) {
			my $tid = delete $s->{id};
			$annot2->put( $tid, $s );
			warn $tid;
		}
		
		$annot2->close();
	}
);


my $nbt = int(scalar(@genes)/($fork*1)) +1;
my $it = natatime $nbt, @genes;

while (my @vals = $it->())
{
		$sid++;
		
		$hrun->{$sid}++;
		warn "NEXT";
		my $pid = $pm->start and next;	
		my $allres;
		my $index19 = Bio::DB::HTS::Faidx->new($file19);
		my $index38 = Bio::DB::HTS::Faidx->new($file38);
		#my $indexprot19 = Bio::DB::HTS::Faidx->new($dir_HG19.$fileprot);
		#my $indexprot38 = Bio::DB::HTS::Faidx->new($dir_HG38.$fileprot);
		
	$nb = 0;
	 my $buffer = new GBuffer;
	foreach my $l (@vals){
		
	next if $l =~ /ENST00000674361/;
	#$l = "ENST00000227471";
	#next unless $l =~ /ENST00000433992/;
	
		$nb++;
	#	next if $nb > 5;
		my @F = split(" ",$l);
		next  if $F[2] ne "transcript";
	#next unless ($F[8] =~ /CCDS/);
		my $chr = $F[0];
		#next if  $F[6] eq "-";
		#next if $nb < 100;
		$chr =~s/chr//;
		$chr = "MT" if $chr eq "M";
		my $coding;
		$coding = 1 if ($F[8] =~ /protein_id/);
		warn $nb."/".$nbt;
		my @z =split(";",$F[8]);
	
		my ($enst) = grep{$_ =~ /transcript_id/} @z;
	
		my ($gene) = grep{$_ =~ /gene_name/} @z;
		$enst =~ s/transcript_id=//;
		$gene =~ s/gene_name=//;
		$enst =~ s/\..*//;
		my $enst_ori = $enst;
		my $y ="";
		$y="_Y"	if $z[0] =~ /PAR_Y/;
    	$enst = $enst.$y;
 		#my ($seq1, $length) = $indexprot19->get_sequence("$enst");
		#my ($seq2, $length2) = $indexprot38->get_sequence("$enst");
		#next if $seq1 && ($seq1 eq $seq2) ;
		
		my $res  = return_hash($enst,$index19,$index38,$buffer,$coding);
		next unless $res;
		
		my @exons =  grep {$_ =~/$enst_ori/} @all_exons;#`zgrep \"$enst_ori\" $file_gff | awk \'\$3 \~ \/exon\/\'`;
		chomp(@exons);
		my $span = Set::IntSpan::Fast::XS->new;
		foreach my $exonl (@exons) {
			my @l = split(" ",$exonl);
			$span->add_range($l[3],$l[4]);
		}
		my $iter = $span->iterate_runs();
		
		$res->{splice_site_span} = Set::IntSpan::Fast::XS->new;
		$res->{essential_splice_site_span} = Set::IntSpan::Fast::XS->new;
    	while (my ( $start, $end ) = $iter->()) {
    	$res->{essential_splice_site_span}->add_range($start-2,$start-1);
		$res->{essential_splice_site_span}->add_range($end+1,$end+2);
		$res->{splice_site_span}->add_range($start-8,$start-2);
		$res->{splice_site_span}->add_range($start,$start+2);
		$res->{splice_site_span}->add_range($end+2,$end+8);
		$res->{splice_site_span}->add_range($end-2,$end);
    		
    	}
		
		$res->{id} = $enst_ori."_".$chr;
	
		push(@{$allres->{data}},$res);
}
	$allres->{run_id} = $sid;
	warn "END $sid";
	$pm->finish( 0, $allres );
}

$pm->wait_all_children();
my $iter = $annot2->rocks->new_iterator->seek_to_first;
while (my ($key, $value) = $iter->each) {
    printf "%s: %s\n", $key, $value;
}
die(Dumper ($hrun)) if keys %$hrun;
#@seq_ids = ("ENST00000227471");
 exit(0);

 sub return_hash {
 	my ($id,$index19,$index38,$buffer,$coding) = @_;
	
	my ($seq1, $length) = $index19->get_sequence("$id");
	my ($seq2, $length2) = $index38->get_sequence("$id");	
	return unless $seq2;
	


 my $res = {};
 
try{

#warn $id unless $aln->{align1};
my ($align1,$align2) = return_alignment($buffer,$seq1,$seq2);

my @correspondence_table;

my $pos2 = 1;
my $pos1 = 1;
my $find;

for my $i (0 .. length($align1)-1) {
    my $char1 = substr($align1, $i, 1);
    my $char2 = substr($align2, $i , 1);
     if ($char1 eq "-" ) {
     	$pos2 ++;
     	next;
     }
     if ($char2 eq "-" ){
     		$res->{miss}->{$pos1} = $char1;
     		$correspondence_table[$pos1] = $pos2;
     		$pos1 ++;
     		next;
     } 
     	$correspondence_table[$pos1] = $pos2;
     	 $pos1 ++;
     	 $pos2++;

}

# Find the position on sequence 2 corresponding to position 2 on sequence 1
#my $pos2 = $correspondence_table[$pos1-1];

for (my $i=0;$i<@correspondence_table;$i++){
	my $dec = $correspondence_table[$i] - $i;
	#warn $dec." i ".$i." c:".$correspondence_table[$i] if $find && $dec ne 0;
	$res->{intspan}->{$dec} =  Set::IntSpan::Fast::XS->new unless exists $res->{intspan}->{$dec};
	$res->{intspan}->{$dec}->add($i);	
#	my $span = Set::IntSpan::Fast::XS->new;
#	warn $i." =>  ".$correspondence_table[$i];
}

$res->{seq38} = $seq2;
}
catch{
	print $id."\n";
	#print $seq1."\n".$seq2."\n";
	die("$length $length2");
};
if($coding){
$res->{protein_coding} = $coding;
my @t = `zgrep $id $dir_HG38/transcripts.pc.fa.gz`;
if (@t){
my @l = split(/\|/,$t[0]);
my ($cds) = grep {$_ =~ /CDS:/} @l;
$cds =~ /CDS:(\d+)-(\d+)/;
$res->{HG38}->{cds}->{start} = $1;
$res->{HG38}->{cds}->{end} = $2;
$res->{HG38}->{cds}->{length} = $length2;
}
#unless (exists $res->{HG38}->{cds}){
	
unless (exists $res->{HG38} ){
		
	my @lines =  `zgrep \"$id\" $file_gff_HG38 | awk \'\$3 \~ \/exon\|CDS\/\'`;
	chomp(@lines);
	readGtf(\@lines,$res);
}
}
#	}
return $res;
 }
 
 sub readGtf {
 	my ($lines,$res) = @_;
# Initialisation des variables
my %transcripts;
 my $tree = Set::IntervalTree->new;
 my $span_cds = Set::IntSpan::Fast::XS->new;
  my $span_exons = Set::IntSpan::Fast::XS->new;
	my $exons = [];
	my $cds;
    # Créer un objet pour chaque gène, transcript et exon
    my $cds_length =0;
    my $start_cds = -1;
    my $end_cds = -1 ;
# Boucle sur chaque ligne du fichier GTF
 my $clen = 0;
 my $strand = 1;
foreach my $line (@$lines){

    # Ignorer les lignes de commentaires commençant par #
    next if $line =~ /^#/;

    # Séparer les colonnes en un tableau
    my @columns = split("\t", $line);

    # Récupérer les informations d'intérêt
    my $seqname   = $columns[0];
    my $source    = $columns[1];
    my $feature   = $columns[2];
    my $start     = $columns[3];
    my $end       = $columns[4];
    my $score     = $columns[5];
    my $cstrand    = $columns[6];
    my $frame     = $columns[7];
    my $attribute = $columns[8];
	$strand = -1 if $cstrand eq "-";
    # Extraire l'identifiant du gène, du transcript et l'ordre de l'exon
    my ($gene_id)      = $attribute =~ /gene_id "(\S+)"/;
    my ($transcript_id) = $attribute =~ /transcript_id "(\S+)"/;
    my ($exon_number)   = $attribute =~ /exon_number "(\d+)"/;
	
  	if ($feature eq 'exon') {
     
		$span_exons->add_range($start,$end);
        
  	}
     elsif ($feature eq 'CDS') {
        	$start_cds = $start if $start_cds == -1;
        	$end_cds = $end if  $end_cds == -1;
         	$start_cds = $start if $start < $start_cds;
         	$end_cds = $end if $end > $end_cds;
    		$cds_length += abs($start-$end) +1;
    		$span_cds->add_range($start,$end);
    }
}
    my $values; 
	my $len =0;
	my $cds_codon;
    foreach my $e (sort {$a->{end} <=> $b->{end}} @$exons){
    		if ($start_cds>= $e->{start} && $start_cds <= $e->{end} ){
    			$cds_codon = abs($start_cds - $e->{start}) + $len +1;
    			last;
    		}
    		$len += abs($e->{start}-$e->{end}) +1;
    }
 
    ($tree,$clen) = return_tree($span_exons);
    my @array = $span_cds->as_array();
    my @texons = $span_exons->as_array();
    my $ltranscript = scalar(@texons);
	if(@array){
		my $cds_start;
		my $cds_end;
	if ($strand == -1 ){
		
			$res->{HG38}->{cds}->{start} = translate_position($array[-1],$tree,$ltranscript,$strand);
			$res->{HG38}->{cds}->{end} = translate_position($array[0],$tree,$ltranscript,$strand);
	}
	else {
			$res->{HG38}->{cds}->{start} = translate_position($array[0],$tree,$ltranscript,$strand);
			$res->{HG38}->{cds}->{end} = translate_position($array[-1],$tree,$ltranscript,$strand);
		}
		$res->{HG38}->{cds}->{length} = abs($res->{HG38}->{cds}->{end} - $res->{HG38}->{cds}->{start})+1;
	}

 	
 }
 
 sub return_tree {
 	my ($span) = @_;
 	my $iter = $span->iterate_runs();
		 my $tree = Set::IntervalTree->new;
		 my $index = 0;
		 my $len=0;
			while (my ( $from, $to ) = $iter->()) {
				
				$tree->insert([$from,$len],$from,$to+1);
				$len += abs($to-$from) +1;
				#push(@t,[$from, $to]);
			}
			return ($tree,$len);
 }
 sub translate_position{
	my ($pos,$tree,$cds_length,$strand) = @_;
	#return -1 unless ($self->getSpanCoding->contains($pos));
	#return -1 if $self->getGenomicSpan->is_empty;
	my $r1 = $tree->fetch($pos,$pos+1);
	my $r = $r1->[0];
	my $rpos = abs($r->[0]-$pos)+1 + $r->[1];
 	$rpos = abs($rpos - $cds_length) +1 if $strand == -1;
 	return $rpos;
}
 
 
 
 sub return_alignment {
 	my ($buffer,$seq1,$seq2) = @_;
 	#if (length($seq1)>50){
 		return return_long_alignment($seq1,$seq2);
 	#}
 	#else {
 	#	return return_normal_alignment($buffer,$seq1,$seq2);
 	#}
 }
 
sub return_normal_alignment {
 	my ($buffer,$seq1,$seq2) = @_;
 	 my @opt = ('timeout'=>100000);
 	my $sw;	
   $sw = $buffer->getNeedlemanWunsch(\@opt);
 	$sw->do_alignment($seq1,$seq2);
	my $aln =  $sw->get_next_hit();
	my $a = $aln->{align1};
	my $b = $aln->{align2};
 	return ($a,$b);	
 }
 sub return_long_alignment {
 	my ($seq1,$seq2) = @_;
 my $suf = time.substr($seq1,1,int(rand(20))+1)."_".substr($seq2,int(rand(20)+2),int(rand(20)+3));	
open (A,">a.$suf.fasta");
print A ">A\n";
print A $seq1."\n";
close(A);
open (A,">b.$suf.fasta");
print A ">B\n";
print A $seq2."\n";
close(A);
my $ff = "AB".$suf.".fasta";
system(qq{/software/distrib/emboss/EMBOSS-6.6.0/emboss/stretcher   -asequence a.$suf.fasta -bsequence b.$suf.fasta -aformat markx3 -outfile /dev/stdout 2>/dev/null | grep -v "#" > $ff});
#warn qq{/software/distrib/emboss/EMBOSS-6.6.0/emboss/stretcher   -asequence a.$suf.fasta -bsequence b.$suf.fasta -aformat markx3 -outfile /dev/stdout 2>/dev/null | grep -v "#" > $ff};
#die();
unlink "a.$suf.fasta";
unlink "b.$suf.fasta";

my $a = `samtools faidx $ff A -n 10000000 | grep -v '>'`;
chomp($a);
my $b = `samtools faidx $ff B -n 10000000 | grep -v '>'`;
chomp($b);
system("rm $ff");
system("rm $ff.fai");
die() unless $a;
die() unless $b;

 return ($a,$b);	
 }
 