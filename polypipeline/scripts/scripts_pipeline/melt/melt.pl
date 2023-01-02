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
my $patient =  $project->getPatient($patient_name);
my $dir_in = $project->getAlignmentDir("bwa");
my $ref = $project->genomeFasta();
my $bam = $patient->getBamFile();
my $dir_out= $project->getCallingPipelineDir("melt-".$patient->name);
my $samtools = $buffer->software("samtools");
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
		system("bgzip $fileout");
		$fileout .= ".gz";
		system("tabix -p vcf $fileout");
		exit(0);
		#write empty vcf;
	}
	my $tvcf = $dir_out."/".$patient->name.".".time.".vcf";
	my $tvcf2 = $dir_out."/".$patient->name.".".time.".2.vcf.gz";
	my $list_file = join(" ",values %$files);
	
 	my $fileout = $project->getVariationsDir("melt")."/".$patient->name.".vcf.gz";
	warn "bcftools concat $list_file | bcftools view  - -U -c 1  > $tvcf; gatk UpdateVCFSequenceDictionary -V $tvcf --source-dictionary /data-isilon/public-data/genome/HG19/fasta/all.dict  --output $tvcf2 --replace; bcftools sort $tvcf2";
	system ("$bcftools concat $list_file | $bcftools view  - -U -c 1  > $tvcf; $gatk UpdateVCFSequenceDictionary -V $tvcf --source-dictionary /data-isilon/public-data/genome/HG19/fasta/all.dict  --output $tvcf2 --replace; $bcftools sort $tvcf2 -O z -o $fileout; $tabix -f -p vcf $fileout");
	exit(0);
	my $f2 = "/data-isilon/bipd-src/pnitschk/git/repository/polypipeline/scripts/scripts_pipeline/melt/1806245.vcf.gz";
	my $v = Bio::DB::HTS::VCF->new( filename => $f2 );
	my $header = $v->header();
	warn $header->get_sample_names()->[0];
	#die if $header->num_samples() ne 1;
	
	my %hashRes;
	my $iter = $v->query("chr10:1-89685296");
	
		while (my $row = $iter->next) {
			my $structType = "mei";
			my $alleles = $row->get_alleles();
			my $allele_index = 1;
			warn $row->reference;
			
			#warn Dumper $row->get_format($header, "DP") ;
			
			#warn Dumper $row->get_format_type($header, "GT") ;
			#warn Dumper $header;
			warn Dumper $row->get_format($header, "GT") ;
			warn Dumper $row->get_format($header, "AD") ;
			warn Dumper $row->get_genotypes($header);
			
			next ;
			die();
			#warn Dumper $row->get_format($header, "AD") ;
			next ;
			for my $a (@$alleles) {
				
				warn Dumper $row;
				
  				printf( "(%s, %s)\n", $a, $row->get_variant_type($allele_index++) ) ;
			}
		}
	die();
	
	sub print_empty_vcf {
	my ( $file, $patient ) = @_;
	confess() unless $patient;
	open( OUT, ">$file" );
	print OUT
"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
	  . $patient->name . "\n";
	close OUT;
	
}
	