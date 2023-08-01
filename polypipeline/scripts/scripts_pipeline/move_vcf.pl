#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use Logfile::Rotate; 
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
use FileHandle;
use List::Util qw(sum);
my $filein;
my $dir;
my $file_bed;
my $project_name;
my $chr_name;
my $fork = 1;
use List::MoreUtils qw(firstidx);


my $log_file;
my $end_ext = "uni";
my $vcf_dir;
my $patient_name;
my $type;
my $merge;
my $calling_type = "unifiedgenotyper";
my $version;
GetOptions(
	'project=s'   => \$project_name,
	"vcf_dir=s"=>\$vcf_dir,
	"patient=s"=>\$patient_name,
	"log=s" =>\$log_file,
	"type=s" =>\$type,
	"merge=s" =>\$merge,
	"method_calling=s" =>\$calling_type,
	"version=s" =>\$version,
);

#my $calling_type = "unifiedgenotyper";
$calling_type = "ion_merge" if $merge;
$calling_type =~  s/_vcf//;

my $date = `date`;
chomp($date);
#$log_file = "toto.log";
if ($log_file){
open (STDOUT,">>".$log_file);
open (STDERR,">>".$log_file);
}

print  colored ['black ON_BRIGHT_MAGENTA'],"======= Move VCF file to prod directory ==== "; 
print  color 'reset';
print  "\n";


my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name,-version=> $version );

my $bgzip = $buffer->software("bgzip");
my $tabix = $buffer->software("tabix");
my $javac = $project->getSoftware('java');
my $reference = $project->genomeFasta();

my $patients = $project->get_list_patients($patient_name);

my $dir_out= $project->getVariationsDir($calling_type);
my $variation_out = $dir_out."/".$project->name().".vcf";
warn $variation_out;
my $total =  $vcf_dir."/".$project_name.".final.vcf";
warn $project->getVariationsDir($calling_type);
if (-e $total){
print  colored ['black ON_BRIGHT_MAGENTA'],"======= Start with BACKUP SNP+INDELS ==== "; 
print  color 'reset';
print  "\n";

backup_vcf($patients);
print  colored ['black ON_BRIGHT_MAGENTA'],"======= Start with MOVE SNP+INDELS ==== "; 
print  color 'reset';
print  "\n";
move_vcf($total,$variation_out,"SNP",$project->getVariationsDir($calling_type));
split2_vcf($variation_out.".gz",$dir_out,$patients);



exit(0);	
}



print  colored ['black ON_BRIGHT_MAGENTA'],"======= Start with SNPs ==== "; 
print  color 'reset';
print  "\n";
my $variations_in = $vcf_dir."/".$project_name.".snps.vcf";
warn $variations_in;
#die();
warn $variations_in;
die($variations_in."--") unless -e $variations_in;


my $indels_in = $vcf_dir."/".$project_name.".indels.vcf";
unless (-e $variations_in ){#&& -e $indels_in){
print  colored ['black ON_BRIGHT_RED'], "file not found $variations_in or $indels_in";
print  color 'reset';
print "\n";
die();
}

move_vcf($variations_in,$variation_out,"SNP",$project->getVariationsDir($calling_type));
warn "test";
split2_vcf($variation_out.".gz",$dir_out,$patients);
exit(0) if $type eq "snp";
print  colored ['black ON_BRIGHT_MAGENTA'],"======= and now the Indels ==== "; 
print  color 'reset';
print  "\n";
my $indel_out = $project->getIndelsDir($calling_type)."/".$project->name().".vcf";
warn $indel_out;

move_vcf($indels_in,$indel_out,"INDEL",$project->getIndelsDir($calling_type));
warn "split\n";
split2_vcf($indel_out.".gz",$project->getIndelsDir($calling_type),$patients);
exit(0);

sub read_header{
	my ($in) = @_;
	my $bcftools = $buffer->software("bcftools");

	open (VCF,"zcat $in |") or die();
	my $format;
	my $vcf_patient;
	
	my $headers;
	my $chr_line;
while (my $line = <VCF>){

	if ($line =~ /^##/){
		push (@$headers,$line);
		next;
	};
	last;
	}
	close VCF;
	return $headers;
}


sub split2_vcf {
	my ($in,$dir_out,$patients) = @_;
	#my $nb = `zgrep -cv "#" $in`;
	#warn $nb;  
	my %hpatients;
	map{$hpatients{$_->name} =0} @{$patients};
	my %families;
	foreach my $p (@$patients){
		$families{$p->name} = $p->family();
	}
	my %col_patients;
	#system ("tabix -f -p vcf $in");
	my @chrs = `tabix -l $in`;
	chomp(@chrs);
	my $headers = read_header($in);
	my $file;
	foreach my $p (@$patients){
	my $nb =0;
	my $rs;
		my $final_out = $dir_out."/".$p->name.".vcf";
		my $final_gz = $final_out.".gz";
		$file->{filename}->{$p->name} =$final_out;
		$file->{filegz}->{$p->name} =$final_gz;
		#open my $fh,">".$final_out or die();
		#$file->{fh}->{$p->name} = $fh;
		$file->{fh}->{$p->name} = FileHandle->new(">".$final_out);
		my $fh = $file->{fh}->{$p->name};
		foreach my $l (@$headers){
			 $fh->print($l);
		}
	}
		
;
	
	open (VCF,"zcat $in  |") or die();
	my $format;
	my $vcf_patient;
	
	my $headers;
	my $chr_line;
	my $chr_current;
while (my $line = <VCF>){
	warn $line;
	if ($line =~ /^##/){
		push (@$headers,$line);
		next;
	};
	chomp($line);
	if ($line =~ /^#CHROM/){
		
		my @array = split(" ",$line);
		
		for (my $i=0 ;$i< @array;$i++){
			unless ($format) {
				push(@$chr_line,$array[$i]);
			}
			$format =$i if $array[$i] eq "FORMAT";
			if ($format){
				$hpatients{$array[$i]} ++;
				$col_patients{$i} = $array[$i];

			} 
		}
		my @toto = grep {$_ ne 1 } values %hpatients;
		if (scalar (@toto) > 1){
			error("problem in vcf with sample  ".join(":",@toto)); 
		
		}
		
		print OUT_SNP $line."\n";
		print OUT_INDEL $line."\n";
		if ($format == 0){
			error ("problem vcf format $line "); 
			
		};
		next;
		}
		
		my @array = split(" ",$line);

	
		if ($array[0] ne $chr_current){
			$chr_current = $array[0];
				 print  colored ['black ON_BRIGHT_MAGENTA'],"\t working $chr_current";  
				 print  color 'reset';
				print  "\n";
			#warn "start $chr_current";
			#last unless $chr_current eq "chr1";
		}
		my @newline;
		for( my $i =0;$i<$format+1;$i++){
		 	push(@newline,$array[$i]); 
		}
		my @format_gt = split(":",$array[$format]);
		my $gt = firstidx{$_ eq "GT"} @format_gt;
		my $ad = firstidx{$_ eq "AD"} @format_gt;
		my $a0 = firstidx{$_ eq "AO"} @format_gt;
		my $dp = firstidx{$_ eq "RO"} @format_gt;
		my %genotype_families;
		my %ad_families;
		for( my $ii =$format +1;$ii<@array;$ii++){
			my @array_sample = split(":",$array[$ii]);
		
			my $sample_gt = $array_sample[$gt];
				
			if ($sample_gt eq "."){
				#$ad_families{$fam} = -1 unless exists $ad_families{$fam};
				next;
			}
			if ($sample_gt eq "./."){
				#$ad_families{$fam} = -1 unless exists $ad_families{$fam};
				next;
			}
			my $name = $col_patients{$ii};
			unless (exists $file->{first}->{$name}){
				my $fh =  $file->{fh}->{$name};
				 $file->{first}->{$name} =1;
				$fh->print(join("\t",@$chr_line)."\t".$name."\n") if $fh ;
			}
			my $fam = $families{$name};

			die() unless $name;
			my $read_all;
			
			if ($ad ne -1){
				my @sample_ad = split(",",$array_sample[$ad]);
				my @ttT = @sample_ad;
				shift(@ttT);
				my $max = sum(@ttT);
				die() unless defined $max;
				
				$ad_families{$fam} =$max  if $max >= $ad_families{$fam};
			}
			else {
				$ad_families{$fam} = $array_sample[$a0] if  $array_sample[$a0] >= $ad_families{$fam};
			}
			
			if ($sample_gt ne "0/0"){
				$genotype_families{$fam} = 1;
				#next;
			}
			
	
		}
		for( my $ii =$format +1;$ii<@array;$ii++){
			my $name = $col_patients{$ii};
			my $fam = $families{$name};
			next unless exists $genotype_families{$families{$name}};
			next if $ad_families{$fam} < 3;
			my @array_sample = split(":",$array[$ii]);
			
			my $sample_gt = $array_sample[$gt];
			my $name = $col_patients{$ii};
			die() unless $name;
			if ($sample_gt eq "./."){
				next;
			}
			if ($sample_gt eq "."){
				next;
			}
			
			my $line2 = join("\t",@newline);
			$line2 .= "\t".$array[$ii]."\n";
			if (exists $file->{fh}->{$name} ) {
			my $fh =  $file->{fh}->{$name};
			#warn "$line2";
		
			  $fh->print($line2);
		}
			#push(@{$vcf_patient->{$name}},$line2);
		}
	}
	
	foreach my $p (@$patients){
		my $name = $p->name;
		$file->{fh}->{$name}->close;
	 my $f1 = $file->{filename}->{$p->name};
	 my $fz = $file->{filegz}->{$p->name};
	system ("$bgzip -f $f1 ");
	system("$tabix -f -p vcf $fz");
	my $nb = `zgrep -cv "#" $fz`;
	my $rs = `zgrep -c "rs" $fz`;
	chomp($nb);
	chomp($rs);
	$nb++;
	my $z = $rs/$nb;
		$z = $z*100;
	 warn  colored ['black ON_BRIGHT_GREEN'],"\t ".$p->name()." ".$nb." ".$z.color 'reset'."\n";  
	 #warn  color 'reset';
	#warn  "\n";
	}
#die();	
#warn "end read ...";
#	foreach my $p (@$patients){
#	my $nb =0;
#	my $rs;
#		die() unless exists $vcf_patient->{$p->name};
#		my $final_out = $dir_out."/".$p->name.".vcf";
#		my $final_gz = $final_out.".gz";
#		open(VCF,">".$final_out) or die();
#		foreach my $l (@$headers){
#			print VCF $l;
#		}
#		print VCF join("\t",@$chr_line);
#		print VCF "\t".$p->name."\n";
#		
#		foreach my $l2 (sort par_chr @{$vcf_patient->{$p->name}}){
#			print VCF $l2;
#			$nb++;
#			$rs++ if $l2 =~/rs/;
#		
#		}
#		close VCF;
#		
#		warn $final_out." ".$final_gz;
#		`bgzip -f $final_out`;
#		`$tabix -p vcf $final_gz -f`;
#		my $z = $rs/$nb;
#		$z = $z*100;
#			print  colored ['black ON_BRIGHT_GREEN'],"\t ".$p->name()." ".$nb." ".$z; 
#		print  color 'reset';
#		print  "\n";
#	}
	
}

sub backup_vcf {
	my ($patients) = @_;
	foreach my $p (@$patients){
	my $final_out = $dir_out."/".$p->name.".vcf";
	my $f = $p->getVariationsFile($calling_type,1);
	next unless $f;
	next unless -e $f;
	system("$Bin/rm_vcf.pl $f")
	}
}

sub move_vcf {
	my ($file_in,$final_out,$type,$dir) = @_;
	warn $file_in;
	my $final_gz = $final_out.".gz";
	my $dir = $dir."/backup";
	mkdir $dir unless -e $dir;
	
	if (-e $final_gz){
		#return;
	my $log = new Logfile::Rotate( File => $final_gz,	
	 								Count => 7,
	 								Gzip  => 'lib',
	 								 Flock  => 'no',
	 								 Dir => $dir
	 								);
		warn "start log";
		$log->rotate();	 								
		}
	warn "cp";	
#system("cp $file_in $final_out");		
#my $cmd = $cmd_select.qq{ --variant $file_in -o $final_out -selectType $type >/dev/null};
##system($cmd);
warn "bgzip $file_in $final_gz";
system("$bgzip -f -c $file_in > $final_gz");
#`bgzip -f $final_out`;
warn "$tabix";
`$tabix -p vcf $final_gz`;
warn "end";
exit(1) unless -e $final_gz;
exit(1) if -s $final_gz ==0;

}

sub par_chr ($$) {
	my ($left,$right) = @_;
	my @l = split(" ",$left);
	my @r = split(" ",$right);
	return give_chr_index($l[0]) <=> give_chr_index($r[0]) ||  $l[1] <=> $r[1];
	 
	
}

sub give_chr_index {
	my ($c) = @_;
	$c =~ s/chr//;
	if ($c eq  'X'){
		return 23;
	}
	if ($c eq 'Y'){
		return 24;
	}
	if ($c eq "MT" || $c eq "M"){
		return 25;
	}
	return $c;
}