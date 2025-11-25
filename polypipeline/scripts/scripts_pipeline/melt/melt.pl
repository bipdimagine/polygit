#!/usr/bin/perl
use FindBin qw($Bin $RealBin);
use strict;

use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
use GBuffer ;
use Data::Dumper;
use Getopt::Long;
use Bio::DB::HTS;
use Bio::DB::HTS::VCF ;
use List::Util qw(max);
### File Open ###

#!/usr/bin/perl
use FindBin qw($Bin $RealBin);
use strict;

use lib "$RealBin/../../../../GenBo/lib/obj-nodb/";
use GBuffer ;
use Data::Dumper;
use Getopt::Long;
use Bio::DB::HTS;
use Bio::DB::HTS::VCF ;
use List::Util qw(max);
### File Open ###



my $project_name;
my $fork;
my $patient_name;
my $version;
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"patients=s" => \$patient_name,
	"version=s" => \$version,
);
#test
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name , -version =>$version);

my $bcftools = $buffer->software("bcftools");
my $bgzip = $buffer->software("bgzip");
my $tabix = $buffer->software("tabix");
my $gatk  = $buffer->software("gatk4");

my $patient =  $project->getPatientOrControl($patient_name);
my $bam_prod = $patient->getBamFile();
my $ref = $project->genomeFasta($bam_prod);
my $fileout = $project->getVariationsDir("melt")."/".$patient->name.".vcf.gz";
unlink $fileout.".tbi" if -e $fileout;
unlink $fileout if -e $fileout;

my $melt_dir= $project->getCallingPipelineDir("melt-".$patient->name.".".time);
my $samtools = $buffer->software("samtools");
my $dir_out = $melt_dir;
my $bam_tmp = $dir_out."/".$patient->name.".bam";

my $done = $dir_out."melt.done";

warn $bam_prod;
my @header = `$samtools view -H  $bam_prod`;
chomp(@header);
my ($umi) = grep{$_ =~ /umi-correction-scheme/} @header;
if ($umi){
	
		my $cmd = qq{samtools view -h $bam_prod }.q{| perl -pe 'if (!/^@/) { @fields = split("\t"); die() if length($fields[9]) ne length($fields[10]); $fields[10] = "I" x length($fields[9]); $_ = join("\t", @fields); }' | samtools view -Sb - >}.$bam_tmp.q{ && samtools index }.$bam_tmp.' -@ 5';
		system($cmd);
		
	die() unless -e $bam_tmp.".bai";
}
elsif ($bam_prod =~/.cram/){
		my $list;
		foreach my $chr (@{$project->getChromosomes}){
			push(@$list,$chr->fasta_name);
		} 
	
		$bam_tmp = $dir_out."/".$patient->name.".cram2bam.bam";
		my $samtools = $buffer->software("samtools");
		warn "$samtools view -@ $fork  -T $ref  $bam_prod ".join(" ",@$list) ." -o $bam_tmp --write-index";
		system("$samtools view -@ $fork  -T $ref  $bam_prod ".join(" ",@$list) ." -o $bam_tmp --write-index");
}
else {
system("ln -s $bam_prod $bam_tmp");
system("ln -s $bam_prod.bai $bam_tmp.bai");
}

my $bam_tmp2 = $bam_tmp;
if ($bam_prod =~/.cram/){
	
}
$dir_out."/".$patient->name.".cram.bam";
my $bed = $dir_out."/".$patient->name.time.".bed";

#my $list = $dir_out."/".$patient->name.".list";

unlink $bed if -e $bed;
open(BED,">$bed");
#open(LIST,">$list");
foreach my $chr (@{$project->getChromosomes}){
	my $chr_name = $chr->fasta_name;
	#next unless $chr_name eq "18";
#	print LIST $chr->name."\t".$chr->fasta_name."\n";
	#getIntSpanCapture
	#getIntSpanCaptureForCalling
	my @line = intspanToBed($chr_name,$chr->getIntSpanCapture(50));
	print BED join("\n",@line)."\n" if @line;
}
##
close(BED);


my $meltd = "/software/distrib/MELT/MELTv2.2.2";
my $melt_soft = $buffer->software("melt");
my $java = $buffer->software("java");
my $melt = "$java -jar $melt_soft Single -a -c 8 ";
my $dir_melt = $buffer->config_path("public_data") . '/repository/'.$project->annotation_genome_version  . '/mei/';
#my $dir_melt = $buffer->public_data_root . '/'.$project->annotation_genome_version  . '/mei/';
my $bedg = $dir_melt."/bed/".lc($project->annotation_genome_version).".genes.bed";

my @files = `ls $dir_melt/me_refs/*.zip`;
chomp @files;
my $list = $dir_out."/list.txt";

#update_method($buffer->dbh,$patient->id);
system("add_calling_method.sh -project=$project_name -patient=$patient_name -method=melt");
open (LIST,">$list");
print LIST join("\n",@files);
close LIST;

#java -jar MELT.jar Single -a -c 8 -h /data-isilon/public-data/genome/HG19/fasta/all.fa -bamfile /data-isilon/sequencing/ngs/NGS2018_2224/HG19/align/bwa/1806245.bam -n ./add_bed_files/1KGP_Hg19/hg19.genes.bed  -w /data-isilon/sequencing/ngs/NGS2018_2224/HG19/align/test -t ./me_refs/list.txt

	system("mkdir $dir_out && chmod a+rwx $dir_out ") unless -e $dir_out;
	
	#system("sambamba slice $bam ".$chr->fasta_name." >$bout && samtools index $bout");

	unless (-e $done){
	if ($patient->isGenome){
		#cat gencode.v43.annotation.gtf |  awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$10,$16,$7}}' | tr -d '";'
		my $bed = $dir_out."/".$patient->name.".genes.bed";
		my $gtf ="";
		#system qq{cat gencode.v43.annotation.gtf |  awk 'OFS="\t" {if ($3=="gene") {print \$1,\$4-1,\$5,\$10,\$16,\$7}}' | tr -d '";'" > $bed};
		warn "$melt -h $ref -bamfile $bam_tmp  -w $dir_out -t $list  -n /data-isilon/software/distrib/centos/MELTv2.2.2/add_bed_files/1KGP_Hg19/hg19.genes.bed -b hs37d5/NC_007605 && touch $done" ;
		system("$melt -h $ref -bamfile $bam_tmp  -w $dir_out -t $list  -n /data-isilon/software/distrib/centos/MELTv2.2.2/add_bed_files/1KGP_Hg19/hg19.genes.bed -b hs37d5/NC_007605 && touch $done");

		#system("$melt -h $ref -bamfile $bam_tmp  -w $dir_out -t $list  -n /data-isilon/software/distrib/centos/MELTv2.2.2/add_bed_files/1KGP_Hg19/hg19.genes.bed -b hs37d5/NC_007605 && touch $done") unless -e $done;
		
	}
	else {
		warn "$melt -h $ref -bamfile $bam_tmp -n $bedg  -w $dir_out -t $list  -exome 1 ";
		my $log = $dir_out."melt.log";
		warn $log;
		system("$melt -h $ref -bamfile $bam_tmp -n $bedg  -w $dir_out -t $list  -exome 1 1>$log 2>>$log && touch $done");
		unless (-e $done) {
			my $no_results_found;
			open (LOG, $log);
			while (<LOG>) {
				my $line = $_;
				if ($line =~ /MELT likely did not find any mobile elements during/) { $no_results_found = 1; }
			}
			close (LOG);
			warn $no_results_found;
			
			if ($no_results_found) {
				warn "cuicui";
				print "\n\nExiting no mobile elements found\n\n";
				warn $done;
				system("touch $done");
				my $fileout = $project->getVariationsDir("melt")."/".$patient->name.".vcf";
				print_empty_vcf($patient,$fileout);
				exit(0);
			}
			die();
		}
	}
	}
	my $files = {ALU=>"$dir_out/ALU.final_comp.vcf",LINE1=>"$dir_out/LINE1.final_comp.vcf",SVA=>"$dir_out/SVA.final_comp.vcf"};
	warn Dumper $files;
	foreach my $f (keys %$files){
		system("gunzip ".$files->{$f}.".gz") if -e $files->{$f}.".gz";
		unless (-e $files->{$f}){
			delete $files->{$f};
			next;
		}
		$files->{$f} = $buffer->gzip_tabix($files->{$f},"vcf");
		my $ff =  $files->{$f};
		delete  $files->{$f} unless -s $ff;
		my $res = ` bcftools view  $ff -U -c 1 2>/dev/null | grep -v "#" | wc -l `;
		chomp($res);
		warn $res;
		delete  $files->{$f} if $res == 0 ;
		
	}
	
	unless (keys %$files){
		my $fileout = $project->getVariationsDir("melt")."/".$patient->name.".vcf";
		print_empty_vcf($fileout,$patient);
		system("$bgzip -f $fileout");
		$fileout .= ".gz";
		system("$tabix -f -p vcf $fileout");
		exit(0);
		#write empty vcf;
	}
	my $tvcf = $dir_out."/".$patient->name.".".time.".vcf";
	my $tvcf2 = $dir_out."/".$patient->name.".".time.".2.vcf.gz";
	my $list_file = join(" ",values %$files);
	

#close (LIST);

my $dict = $ref;
$dict =~ s/\.fa/\.dict/;
$bed = $buffer->gzip_tabix($bed,"bed");
my $fileout = $project->getVariationsDir("melt")."/".$patient->name.".vcf.gz";

#warn "$bcftools concat $list_file | $bcftools view - -R $bed | $bcftools view  - -U -c 1  > $tvcf; $gatk UpdateVCFSequenceDictionary -V $tvcf --source-dictionary /data-isilon/public-data/genome/HG19/fasta/all.dict  --output $tvcf2 --replace; $bcftools sort $tvcf2 -O z -o $fileout; tabix -f -p vcf $fileout";
my $cmd = qq{$bcftools concat -a $list_file  | perl -lane 's/GL,Number=\\d/GL,Number=G/;print \$_' | $bcftools view  - -U -c 1  > $tvcf;$gatk UpdateVCFSequenceDictionary -V $tvcf --source-dictionary $dict  --output $tvcf2 --replace 2>/dev/null;};
system ($cmd);
warn $cmd;
	
	my $tvcf3 = $dir_out."/".$patient->name.".".time.".3.vcf.gz";
	my $cmd2 = qq{$bcftools sort $tvcf2 -O z -o $tvcf3 ;$tabix -p vcf $tvcf3; $bcftools view $tvcf3 -R $bed -O z -o $fileout};
	system($cmd2);
	warn $cmd2;
	$buffer->gzip_tabix($fileout,"vcf");
	
	warn "OK " ;
	
	exit(0);
	
	
	sub print_empty_vcf {
	my ( $file, $patient ) = @_;
	confess() unless $patient;
	open( OUT, ">$file" );
	print OUT
"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
	  . $patient->name . "\n";
	close OUT;
	
}





sub intspanToBed{
	my ($chr_name,$intspan) = @_;
	my $iter = $intspan->iterate_runs();
	my @tt;
	my $size = 0;
    while (my ( $from, $to ) = $iter->()) {
    	#warn "+".$from."-".$to if ($from<=31263320);
    	warn $from."-".$to if ($from<=31263320 && $to >= 31263320);
    	$size += abs($from-$to);
    		push(@tt,$chr_name."\t".$from."\t".$to);
    	
    }
	#	my @tt = map{$_ = $chr->ucsc_name."\t".$_} split(";",$intspan->as_string({ sep => ";", range => "\t" }) );
	return @tt;
}
sub update_method {
	my ($dbh,$patient_id) = @_;
	my $query = qq{
		insert into PolyprojectNGS.patient_methods (patient_id,method_id) SELECT $patient_id ,method_id as methodId FROM PolyprojectNGS.methods where methods.name="melt" ;
	};
	warn $query;
	$dbh->do($query) ;

}
	
