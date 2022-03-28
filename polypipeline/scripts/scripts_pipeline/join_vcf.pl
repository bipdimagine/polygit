#!/usr/bin/perl
use FindBin qw($Bin);
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../packages/";
use  File::Temp;
use Data::Dumper;
use Getopt::Long;
use Carp;
use GBuffer;
use colored;
use Storable qw(store retrieve freeze);
use Term::ANSIColor;
use threads;
use Thread::Queue;
use Set::IntSpan::Fast::XS;
  use IPC::Open2;
  
my $project_name;
my $final_vcf;
my $log_file;
my $window_length;
my $list_patients;
my $fork;

GetOptions(
	'project=s'   => \$project_name,
	"log=s" =>\$log_file,
	"vcf=s" => \$final_vcf,
	"window=s" => \$window_length,
	"patient=s"=>\$list_patients,
	"fork=s" =>\$fork,
);
my $date = `date`;
chomp($date);
confess();
my $buffer = GBuffer->new();
my $project = $buffer->newProject( -name => $project_name );

my $vcffirstheader = $buffer->software("vcffirstheader");
my $vcfstreamsort = $buffer->software("vcfstreamsort");
my $vcfuniq  = $buffer->software("vcfuniq");
my $bcftools = $buffer->software("bcftools");
my $bgzip = $buffer->software("bgzip");
my $tabix = $buffer->software("tabix");

if ($log_file){
	open (STDOUT,">>".$log_file);
}
#$SIG{INT} = \&interrupt;




my $dir_gvcf_out= $project->getCallingPipelineDir("gvcf");

#open(OUT,">$bed_file");
#open(OUT2,">$gvcf_file");
#foreach my $chr (@{$project->getChromosomes()}){
#		foreach my $window (@{$chr->getWindowCaptureForCalling(100,$window_length)}){
#			my $filename = $dir_out."/vcf/".$chr->ucsc_name.".".$window->{start}.".".$window->{end}.".vcf";
#			my $filename2 = $dir_out."/vcf/".$chr->ucsc_name.".".$window->{start}.".".$window->{end}.".none";
#			next if  -e $filename2;
#		#	my $gfilename = $dir_out."/gvcf/".$chr->ucsc_name.".".$window->{start}.".".$window->{end}.".g.vcf";
#			die($filename) unless -e $filename;
#			print OUT "$filename\n";
#			#if (-e $gfilename) {
#			#		print OUT2 "$gfilename\n";
#			#}
#		}
#}	
#close(OUT);
#close(OUT2);


#my $vcf_temp = getTmpFile($dir_out,$project_name,"vcf");
#my $cmd = qq{for i in `cat $bed_file`;do cat \$i ;done | $vcffirstheader | $vcfstreamsort -w 1000 | $vcfuniq > $vcf_temp};
#system ($cmd);
##	my $final_gvcf = $dir_out."/".$patient->name.".g.vcf";
##my $cmd2 = qq{for i in `cat $gvcf_file`;do cat \$i ;done | $vcffirstheader | $vcfstreamsort -w 1000 | $vcfuniq > $final_gvcf};
##system ($cmd2." && gzip $final_gvcf");
##warn $cmd2;

#my @bed;
#
#system("$bgzip $vcf_temp && $tabix -p vcf $vcf_temp.gz");
#foreach my$chr (@{$project->getChromosomes()}){
#	my $intspan = $chr->getIntSpanCaptureForCalling(250);
#	map{push(@bed,$chr->ucsc_name."\t".$_."\n")} split(";",$intspan->as_string({ sep => ";", range => "\t" }));	
#	
#}
#
#my $bed_file = getTmpFile($dir_out,$project_name,"bed");
#
#open(BED,">$bed_file");
#foreach my $b (@bed){
#	print BED $b;
#}
#close BED;
#system (qq{$bcftools filter $vcf_temp.gz -R $bed_file max -e "MAX(AD[1]+AD[0])<7" >$final_vcf});


#unlink "$vcf_temp.gz" if -e "$vcf_temp.gz";
#unlink "$vcf_temp.gz.tbi" if -e "$vcf_temp.gz.tbi";
#unlink $bed_file if -e $bed_file;
print "GVCF \n";
my $patients =  $project->get_list_patients($list_patients);
	my $dir_out_gvcf = $dir_out."/gvcf/";
	my @gvcf;
	 my $pm1 = new Parallel::ForkManager($fork);
	 
	 my $chr_gvcf_files;
	 
	 $pm1->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    	 # $pr->update($c++);
	
 		#$pr->write() if $verbose;
 		my $p = $data->{patient};
 		my $c = $data->{chr};
 		
 		$chr_gvcf_files->{$p}->{$c} = $data->{file};
 		#push(@gvcf,$data->{file});
    		
    }
  );
	 
	 my $vcffirstheader = $buffer->software("vcffirstheader");
	my $vcfstreamsort = $buffer->software("vcfstreamsort");
	my $vcfuniq  = $buffer->software("vcfuniq");
	
	foreach my $patient (@$patients){
			my $first;
			my $final_gvcf = $dir_out."/".$patient->name.".g.vcf";
				#push(@gvcf,$final_gvcf);
				my $gz = $final_gvcf.".gz";	
				 if (-e $gz){
				 	push(@gvcf,$gz);
				 	next;
				 }
				
				 warn $patient->name;
		#	next if -e $gz;
	
		my $final_gvcf = $dir_out."/".$patient->name.".g.vcf";
		
				#push(@gvcf,$final_gvcf);
		my $gz = $final_gvcf.".gz";
		#unlink $gz if -e $gz;
		unlink $final_gvcf if -e $final_gvcf;
		my $temp_gvcf = 	$dir_out."/".$patient->name.".g.vcf.tmp";		
	#	my $pid = pen3(\*CHLD_IN, \*CHLD_OUT, \*CHLD_ERR, '  | $vcffirstheader | $vcfstreamsort -w 1000 | $vcfuniq ')
		foreach my $chr (@{$project->getChromosomes()}){
				my $windows = $chr->getWindowCaptureForCalling(250,$window_length);
				#warn $chr->name;
				my $chr_gvcf = $final_gvcf.".".$chr->name;
				my $files;
				foreach my $window (@{$windows}){
				next if $window->{intspan}->is_empty();
				my $outfile = $patient->getWindowGvcf($window);
				confess($outfile) unless -e $outfile;
				push(@$files,$outfiles);
				}
				my $pid = $pm1->start and next;
				unless (-e $chr_gvcf){
				unlink $chr_gvcf if -e $chr_gvcf;
				open(GLOBAL, " | $vcffirstheader | $vcfstreamsort | $vcfuniq > $chr_gvcf");
				
	
	 			my @all_vcfs;
				foreach my $outfile (@$files){
					open (OUT,"$outfile");
					while(<OUT>){
						print GLOBAL $_;
					}
					close OUT;

				unlink $outfile.".idx" if -e  $outfile.".idx";
			
				
		}
		close GLOBAL;
				}
		my $f;
		$f->{file} = $chr_gvcf;
		$f->{chr} = $chr->name;
		$f->{patient} = $patient->name;
		$pm1->finish(0,$f);
	}
	$pm1->wait_all_children();
	 
	}
 
 my @gvcf;
	 my $pm2 = new Parallel::ForkManager($fork);
	 
	 
	 $pm2->run_on_finish(
    sub { my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;

 		push(@gvcf,$data->{file});
    		
    }
  );
  
	foreach my $patient (@$patients){
		my $first;
		my $final_gvcf = $dir_out."/".$patient->name.".g.vcf";
		my $gz = $final_gvcf.".gz";
		
		my $files;
		my $pid = $pm2->start and next;
		if (-e $gz){
		system(" $tabix -p vcf $gz"); unless -e $gz.".tbi";
		my $f;
		$f->{file} = $gz;		
		$pm2->finish(0,$f);
		
		}
		foreach my $chr (@{$project->getChromosomes()}){
			my $file = $chr_gvcf_files->{$patient->name}->{$chr->name};
			confess($file) unless -e $file;
			push(@$files,$file);
		
		}
		my $string = join(" ",@$files);
		system ("cat $string | $vcffirstheader > $final_gvcf");
		warn "cat $string | $vcffirstheader > $final_gvcf";
		for my $f (@$files){
			unlink $f;
		}
		system("$bgzip $final_gvcf && $tabix -p vcf $gz");		
		my $f;
		$f->{file} = $gz;		
		$pm2->finish(0,$f);
	}
$pm2->wait_all_children();
	
#	push(@gvcf,$final_gvcf);
#	warn "end concat ".$patient->name;
#	system("$bgzip $final_gvcf && $tabix -p vcf $gz");
#	warn $gz;
#	my $f;
#	$f->{file} = $gz;		
#	
#	
#	#system("$bgzip $final_gvcf");
#	my $final_gvcf = $dir_out."/".$patient->name.".g.vcf";
#		
#		
#		#warn $final_gvcf;
#		#my $cmd2 = qq{for i in `cat $list_file`;do cat \$i ;done | $vcffirstheader > $final_gvcf};
#
#		
#		
#}
# $pm1->wait_all_children();
my $string_variant = join (" --variant ",@gvcf);

my $javac = $project->getSoftware('java');
my $gatk  = $project->getSoftware('gatk');
my $reference = $project->genomeFasta();
my $vcf1 = getTmpFile($dir_out,"all","vcf");
my $cmd2 = qq{$javac -jar $gatk -T GenotypeGVCFs -R $reference --out $vcf1 --variant $string_variant -nt $fork };
warn $cmd2;
system($cmd2);
system("$bgzip $vcf1");

system (qq{$bcftools filter $vcf1.gz   -e "MAX(AD[1]+AD[0])<7 || MAX(AD[1])<5" >$final_vcf});

exit(0);
foreach my $patient (@{$project->getPatients}){
	my $bed_file = $dir_out."/file.".$patient->name."list";
	my $final_gvcf = $dir_out."/".$patient->name.".g.vcf";
	open(OUT,">$bed_file");
	foreach my $chr (@{$project->getChromosomes()}) { 
		foreach my $window (@{$chr->getWindowCaptureForCalling(100,$window_length)}){
			my $filename = $dir_out."/gvcf/".$patient->name.".".$chr->ucsc_name.".".$window->{start}.".".$window->{end}.".g.vcf";
			die($filename) unless -e $filename;
			
			print OUT "$filename\n";
			
		}
		
		
	}
	close(OUT);
	warn $bed_file;
	my $cmd = qq{for i in `cat $bed_file`;do cat \$i ;done | $vcffirstheader  > $final_gvcf};
	system ($cmd);
	warn 	$final_gvcf;
	die();
}


sub getTmpFile {
	my ($dir,$chr_name,$ext) = @_;
	die() unless $ext;
	die() unless $chr_name;
	#confess($chr_name) unless $chr_name !~/[0,1,2,3,4,5,6,7,8,9,10,11,12,13,15,15,16,17,18,19,20,21,22,X,Y,M]/;
	$dir .="/$chr_name";
	system ("mkdir -p $dir") unless -d $dir;
	 my $file_cat_tmp =  File::Temp->new( TEMPLATE =>"TMP.XXXXXXX",
                        DIR => $dir,
                        SUFFIX => ".$ext");
  return $file_cat_tmp->filename();
}