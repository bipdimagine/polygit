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
GetOptions(
	'project=s'   => \$project_name,
	"fork=s"  => \$fork,
	"patients=s" => \$patient_name,
);
#test
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );

my $bcftools = $buffer->software("bcftools");
my $bgzip = $buffer->software("bgzip");
my $tabix = $buffer->software("tabix");
my $gatk  = $buffer->software("gatk4");

my $patient =  $project->getPatient($patient_name);
my $dir_in = $project->getAlignmentDir("bwa");
my $ref = $project->genomeFasta();
my $bam = $patient->getBamFile();

my $melt_dir= $project->getCallingPipelineDir("melt-".$patient->name."-".time);

my $samtools = $buffer->software("samtools");
my $dir_out = $melt_dir;
my $bam_tmp = $dir_out."/".$patient->name.time.".bam";
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
#system("sambamba slice $bam -L $bed -o $bam_tmp;samtools index $bam_tmp ");
#warn $bam_tmp;
#die();

my $meltd = "/software/distrib/MELT/MELTv2.2.2";
my $melt = $buffer->software("melt");
my $java = $buffer->software("java");
my $melt = "$java -jar $melt Single -a -c 8 ";
my $dir_melt = $buffer->config->{'public_data'}->{root} . '/repository/'.$project->annotation_genome_version  . '/mei/';
my $bed = $dir_melt."/bed/hg19.genes.bed";
my @files = `ls $dir_melt/me_refs/*.zip`;
my $bcftools = $buffer->software("bcftools");
my $bgzip =$buffer->software("bgzip"); 
my $tabix =$buffer->software("tabix"); 
my $gatk=$buffer->software("gatk4");
chomp @files;
my $list = $dir_out."/list.txt";
update_method($buffer->dbh,$patient->id);
system("add_calling_method.sh -project=$project_name -patient=$patient_name -method=melt");
open (LIST,">$list");
print LIST join("\n",@files);
close LIST;

#java -jar MELT.jar Single -a -c 8 -h /data-isilon/public-data/genome/HG19/fasta/all.fa -bamfile /data-isilon/sequencing/ngs/NGS2018_2224/HG19/align/bwa/1806245.bam -n ./add_bed_files/1KGP_Hg19/hg19.genes.bed  -w /data-isilon/sequencing/ngs/NGS2018_2224/HG19/align/test -t ./me_refs/list.txt

	system("mkdir $dir_out && chmpd a+rwx $dir_out ") unless -e $dir_out;
	
	#system("sambamba slice $bam ".$chr->fasta_name." >$bout && samtools index $bout");
	warn "$melt -h $ref -bamfile $bam -n $bed  -w $dir_out -t $list  -exome 1";
	system("$melt -h $ref -bamfile $bam -n $bed  -w $dir_out -t $list  -exome 1");
	my $files = {ALU=>"$dir_out/ALU.final_comp.vcf",LINE1=>"$dir_out/LINE1.final_comp.vcf",SVA=>"$dir_out/SVA.final_comp.vcf"};
	
	foreach my $f (keys %$files){
		warn $files->{$f};
		my $ff =  $files->{$f};
		delete  $files->{$f} unless -s $ff;
		my $res = ` bcftools view  $ff -U -c 1 2>/dev/null | grep -v "#" | wc -l `;
		chomp($res);
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



system("$bgzip $bed ; $tabix $bed.gz");
$bed = $bed.".gz";
	
	
 	my $fileout = $project->getVariationsDir("melt")."/".$patient->name.".vcf.gz";

	#warn "$bcftools concat $list_file | $bcftools view - -T $bed | $bcftools view  - -U -c 1  > $tvcf; $gatk UpdateVCFSequenceDictionary -V $tvcf --source-dictionary /data-isilon/public-data/genome/HG19/fasta/all.dict  --output $tvcf2 --replace; $bcftools sort $tvcf2 -O z -o $fileout; tabix -f -p vcf $fileout";
	my $cmd = qq{$bcftools concat $list_file | $bcftools view - -T $bed | perl -lane 's/GL,Number=\\d/GL,Number=G/;print \$_' | $bcftools view  - -U -c 1  > $tvcf; $gatk UpdateVCFSequenceDictionary -V $tvcf --source-dictionary /data-isilon/public-data/genome/HG19/fasta/all.dict  --output $tvcf2 --replace; $bcftools sort $tvcf2 -O z -o $fileout; $tabix -f -p vcf $fileout};
	warn $cmd;
	system ($cmd);
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
    warn $chr_name." ".$size;
	#	my @tt = map{$_ = $chr->ucsc_name."\t".$_} split(";",$intspan->as_string({ sep => ";", range => "\t" }) );

sub update_method{
	my ($dbh,$patient_id) = @_;
	my $query = qq{
		insert into PolyprojectNGS.patient_methods (patient_id,method_id) SELECT $patient_id ,method_id as methodId FROM PolyprojectNGS.methods where methods.name="melt" ;
	};
	warn $query;
	$dbh->do($query) ;

}
	