#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/../../../GenBo/lib/obj-nodb/packages/";
use lib "$Bin/../../../GenBo/lib/GenBoDB/writeDB/";
use lib "$Bin/packages";
use Data::Dumper;
use Getopt::Long;
use Carp;
use IO::Prompt;
use colored;
use Cwd 'abs_path';
use Vcf;
use LWP::UserAgent;
use URI;
use JSON;
use HTTP::Request::Common qw(POST);
use List::Util qw( min max );


my $project_name = 'NGS2024_8164'; # 'NGS2024_8164','NGS2022_6048';
my $patient_name = 'OHO_Emm'; # 'DCAS'
my $no_exec;
my $cpu = 8;

GetOptions(
	'project=s'		=> \$project_name,		# Projet
	'patient=s'		=> \$patient_name,
	'no_exec'		=> \$no_exec,
	'cpu=i'			=> \$cpu,
);


use GBuffer;
my $buffer = new GBuffer;

my $project = $buffer->newProject(-name=>$project_name) or confess ("Can't open project '$project_name': $!");
warn $project_name;

my $release = $project->annotation_genome_version;
my $hmmfile;
$hmmfile = 'homo_sapiens_dfam.hmm' if ($release =~ /^HG/);
$hmmfile = 'mus_musculus_dfam.hmm' if ($release =~ /^MM/);
$hmmfile = 'rattus_norvegicus_dfam.hmm' if ($release =~ /^RN/);
$hmmfile = prompt("Choose a dfam.hmm file :") unless $hmmfile;
confess ("$hmmfile does not exists or is not a dfamm.hmm file: $!") unless (-e '/home/mperin/git/polygit/polyscripts/melo/dfamscan/'.$hmmfile and $hmmfile =~ /_dfam.hmm$/);


my $pat = $project->getPatient($patient_name) or confess ("Can't find patient '$patient_name' in project $project_name: $!");
my $pat_name = $pat->name;
warn $pat_name;
#my $vcf_file = $pat->getVariationsFile("dragen-sv");
my $vcf_file = $project->getVariationsDir("dragen-sv")."$pat_name.sv.vcf.gz";
confess ("Can't find dragen-sv variation file for '$pat_name' : $vcf_file. $!") unless (-e $vcf_file);
warn $vcf_file;
my $vcf = Vcf->new(file=>$vcf_file) or confess ("Can\'t open vcf file $vcf_file: $!");
$vcf->parse_header();

my $fasta_path = "/home/mperin/git/polygit/polyscripts/melo/dfamscan/$pat_name\_ins.fasta";
warn $fasta_path;
open(my $fasta_ins, '>', $fasta_path);
while (my $x=$vcf->next_data_hash) {
	if ($x->{'INFO'}->{'SVTYPE'} eq 'INS') {
		my $lsvinsseq = $x->{'INFO'}->{'LEFT_SVINSSEQ'};
		my $rsvinsseq = $x->{'INFO'}->{'RIGHT_SVINSSEQ'};
		my $alt_svinsseq = $x->{'ALT'}[0];
				
		# Récupère les séquences des insertions (si >= 10 nt)
		# todo: est-ce qu'on impose un taille min de la seq pour faire la recherche Dfam ? sur le site seq min = 50 nt
		print $fasta_ins '>'.$x->{'ID'}."-ALT\n$alt_svinsseq\n" if ($alt_svinsseq ne '<INS>'); # && length $alt_svinsseq >= 50);
		print $fasta_ins '>'.$x->{'ID'}."-LEFT_SVINSSEQ\n$lsvinsseq\n" if $lsvinsseq; # if (length $lsvinsseq >= 50);
		print $fasta_ins '>'.$x->{'ID'}."-RIGHT_SVINSSEQ\n$rsvinsseq\n" if $rsvinsseq; # if (length $rsvinsseq >= 50);
	}
}
close($fasta_ins);
$vcf->close();


# Lance dfamscan
my $dfamscan_outfile = "/home/mperin/git/polygit/polyscripts/melo/dfamscan/$pat_name.dfamscan.out";
#warn $dfamscan_outfile;
my $cmd = "cd /home/mperin/git/polygit/polyscripts/melo/dfamscan/; perl dfamscan.pl --fastafile=$fasta_path --hmmfile=$hmmfile --dfam_outfile=$dfamscan_outfile -cpu=$cpu";
warn $cmd;
# todo: lancer sur cluster ? (option dfamscan --cpu)
system($cmd) unless ($no_exec);
my $dfamscan_results = parse_dfamscan ($dfamscan_outfile) unless ($no_exec);
#print Dumper $dfamscan_results;

my $dfamscan_vcf = "/home/mperin/git/polygit/polyscripts/melo/dfamscan/$pat_name.dfamscan.vcf";
warn $dfamscan_vcf;
write_results ($dfamscan_results, $vcf_file, $dfamscan_vcf) unless ($no_exec);





sub parse_dfamscan {
	my $dfam_out_file = shift;
	
	open(my $dfam_out, '<', $dfam_out_file) or confess("Can't open $dfam_out_file: $!");
	my @columns = ('target name', 'acc', 'query name', 'bits', 'e-value', 'bias', 'hmm-st', 'hmm-en', 'strand', 'ali-st', 'ali-en', 'env-st', 'env-en', 'modlen', 'description of target');
	my $dfamscan_results;
	while (my $line = readline($dfam_out)) {
		next if $line =~ /^#/;
		my @hit = (split(/ {2,}/, $line));
		chomp @hit;
		confess(scalar @columns.' columns expected, got '.scalar @hit.":\n".join("\t",@hit)) if (scalar @hit != scalar @columns);
		my %hit = map { $columns[$_] => $hit[$_] } 0..$#columns;
		push (@$dfamscan_results, \%hit);
	}
	close($dfam_out);
	return $dfamscan_results;
}



sub write_results {
	my ($dfamscan_results, $vcf_file, $out_path) = @_;

	open( my $out_file, '>', $out_path ) or confess("Can't write in $out_path: $!");
	
	# copy header and add the corresponding meta-informations to describe the INFO entries to be added
	open (my $vcf_fh, "zcat $vcf_file |");
	while ( my $x = <$vcf_fh> ) {
		next unless $x =~ /^#/;
		if ($x =~ /^##fileDate=/){
			print $out_file "##fileDate=".localtime."\n";
			next;
		}
		if ($x =~ /^##source=/){
			print $out_file "##source=DFAMSCAN\n";
			next;
		}
		if ($x =~ /^##ALT=<ID=INS,Description="Insertion">$/){
			print $out_file q{##ALT=<ID=INS:DFAM,Description="Insertion of transposable element">
};
#			print $out_file q{##ALT=<ID=INS:ME,Description="Insertion of mobile element">
#};
#			print $out_file q{##ALT=<ID=INS:DFAM:ALU,Description="Insertion of ALU element">
###ALT=<ID=INS:DFAM:LINE1,Description="Insertion of LINE1 element">
###ALT=<ID=INS:DFAM:SVA,Description="Insertion of SVA element">
#};
			next;
		}
		next if ($x =~ /^##ALT=/);
		print $out_file $x;
		if ($x =~ /^##INFO=<ID=JUNCTION_QUAL/){
#			print $out_file q{##INFO=<ID=DFAMSCAN_TARGET_NAME,Number=.,Type=String,Description="Target name">
###INFO=<ID=DFAMSCAN_ACC,Number=.,Type=String,Description="Accesion">
###INFO=<ID=DFAMSCAN_QUERY,Number=.,Type=String,Description="Name of the query sequence">
###INFO=<ID=DFAMSCAN_TARGET_DESCRIPTION,Number=.,Type=String,Description="Description of the target">
###INFO=<ID=DFAMSCAN_BITS,Number=.,Type=Float,Description="Bit score">
###INFO=<ID=DFAMSCAN_E-VALUE,Number=.,Type=Float,Description="E-value">
###INFO=<ID=DFAMSCAN_BIAS,Number=.,Type=Float,Description="Bias">
###INFO=<ID=DFAMSCAN_HMM-ST,Number=.,Type=Integer,Description="Start of the hmm alignement">
###INFO=<ID=DFAMSCAN_HMM-EN,Number=.,Type=Integer,Description="End of the hmm alignement">
###INFO=<ID=DFAMSCAN_STRAND,Number=.,Type=String,Description="Strand">
###INFO=<ID=DFAMSCAN_ENV-ST,Number=.,Type=Integer,Description="Start of the alignement on the entry sequence">
###INFO=<ID=DFAMSCAN_ENV-EN,Number=.,Type=Integer,Description="End of the alignment on the entry sequence">
###INFO=<ID=DFAMSCAN_ALI-ST,Number=.,Type=Integer,Description="Start of alignment on the model sequence">
###INFO=<ID=DFAMSCAN_ALI-EN,Number=.,Type=Integer,Description=>"End of alignment on the model sequence">
###INFO=<ID=DFAMSCAN_MODLEN,Number=.,Type=Integer,Description="Model length">
#};
			# todo: mettre name, start, end, strand, score ?, 
			print $out_file q{##INFO=<ID=DFAMINFO,Number=4,Type=String,Description="Mobile Element info of the form NAME,START,END,POLARITY; If START or END is unknown, will be -1; If POLARITY is unknown, will be 'null'">
};
#			print $out_file q{##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile Element info of the form NAME,START,END,POLARITY; If START or END is unknown, will be -1; If POLARITY is unknown, will be 'null'">
#};
			next;
		}
	}
	close($vcf_fh);
	
	# add INFO entries for the Dfamscan results
	my $vcf = Vcf->new( file => $vcf_file );
	$vcf->parse_header();
	while ( my $x = $vcf->next_data_array ) {
		my @dfamscan_query = grep { $_->{'query name'} =~ $$x[2] } @$dfamscan_results;
		while (my $hit = shift @dfamscan_query) {
#			my $TE_type = 'OTHER';
#			$TE_type = 'ALU' if ($hit->{'target name'} =~ /^Alu/);
#			my %added_info = 
			$$x[7] = $vcf->add_info_field( $$x[7], 
#				'DFAMSCAN_TARGET_NAME'			=> $hit->{'target name'},
#				'DFAMSCAN_ACC'					=> $hit->{'acc'},
#				'DFAMSCAN_QUERY'				=> $hit->{'query name'}, # split('-', $hit->{'query name'})[1]; 
#				'DFAMSCAN_TARGET_DESCRIPTION'	=> '"'.$hit->{'description of target'}.'"',
#				'DFAMSCAN_BITS'					=> $hit->{'bits'},
#				'DFAMSCAN_E-VALUE'				=> $hit->{'e-value'},
#				'DFAMSCAN_BIAS'					=> $hit->{'bias'},
#				'DFAMSCAN_HMM-ST'				=> $hit->{'hmm-st'},
#				'DFAMSCAN_HMM-EN'				=> $hit->{'hmm-en'},
#				'DFAMSCAN_STRAND'				=> $hit->{'strand'},
#				'DFAMSCAN_ENV-ST'				=> $hit->{'env-st'},
#				'DFAMSCAN_ENV-EN'				=> $hit->{'env-en'},
#				'DFAMSCAN_ALI-ST'				=> $hit->{'ali-st'},
#				'DFAMSCAN_ALI-EN'				=> $hit->{'ali-en'},
#				'DFAMSCAN_MODLEN'				=> $hit->{'modlen'},
				'SVTYPE'				=> "TE",
				'DFAMINFO'				=> $hit->{'target name'}.','.$hit->{'ali-st'}.','.$hit->{'ali-en'}.','.$hit->{'strand'}.';',
#				'MEINFO'				=> $hit->{'target name'}.','.$hit->{'ali-st'}.','.$hit->{'ali-en'}.','.$hit->{'strand'}.';',
			);
#			$$x[7] = $vcf->add_info_field( $$x[7], %added_info );
			print $out_file join( "\t", @$x ) . "\n";
		}
	}
	close($out_file);
	$vcf->close();
	
	return;
}






#	# Réécrire les résultats dans un csv (séparé par \t plus facile à parser)
#	my $dfamscan_csv_path = "/home/mperin/git/polygit/polyscripts/melo/dfamscan/$pat_name\_dfamscan.csv";
#	warn $dfamscan_csv_path;
#	my @dfamscan_keys = ('target name', 'acc', 'query name', 'bits', 'e-value', 'bias', 'hmm-st', 'hmm-en', 'strand', 'ali-st', 'ali-en', 'env-st', 'env-en', 'modlen', 'description of target');
##	print '#'.join("\t", @dfamscan_keys)."\n";
#	print $dfamscan_csv '#'.join("\t", @dfamscan_keys)."\n";
#	foreach my $hit (@$dfamscan_results) {
##		print join("\t", map {$hit->{$_}} @dfamscan_keys)."\n";
#		print $dfamscan_csv join("\t", map {$hit->{$_}} @dfamscan_keys)."\n";
#	}
#	close($dfamscan_csv);




