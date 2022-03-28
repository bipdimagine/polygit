#!/usr/bin/perl
use CGI qw/:standard :html3/;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/GenBoDB";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/pedigree";
use GBuffer;
use GenBoStorable;
use Data::Dumper;
use export_excel;  
use export_data;
use Sys::Hostname;
use get_variations;
use kyoto::public_data; 
use kyoto::db_public;
use kyoto::kyoto_database;
use pedigree;  
use validationQuery; 
use Storable qw/freeze thaw nfreeze nstore_fd nstore retrieve/;


my %types = ( variations => "variations" );

my $cgi    = new CGI();
my $buffer = new GBuffer;



my $host = hostname;



my $project_name = $cgi->param('project');
my $user = $cgi->param('user');
my $type_name = $cgi->param('type');
$type_name = "variations";
my $type2    = $cgi->param('type2');
my $select_ids   = [split(",",$cgi->param('ids'))];
#creation du projet  partir du buffer
my $project = $buffer->newProject( -name => $project_name );
die( "unknown project" . $project_name ) unless $project;
#$project->getMethods($buffer->getType("variation"));


my $ped_file = $project->getRootDir()."/../".$project->name().".ped";
my $pedigree ={};
my $ped_fam =[];

if (-e $ped_file){
	my $patients = $project->getPatients();
	my @names = map{$_->name()} @$patients;
	($pedigree,$ped_fam) = pedigree::parse_ped($ped_file);
}

my $temp_type_name = $type_name;
my @tchr = split(" ",$cgi->param('chromosome'));
my @final_data;

foreach my $nchr (@tchr) {
my $chr = $project->getChromosome($nchr);
die() unless $chr;
my $data;

($data) = get_variations::getIds($buffer,$project,$chr,$temp_type_name,$select_ids);
	update_pedigree($project,$data,$pedigree);
	export_data::update_deja_vu($project,$data,$user);
	#update_deja_vu($project,$data);
	update_genes_names($project,$data);
	eval{
	update_valid( $project,$data,$user);
	};
push(@final_data,@$data);
}
	my $type_label = "id";
if ($type_name eq "patients"){
	by_patients($project,\@final_data);
}
if ( $cgi->param('xls') == 1 ) {
	export_data::variation_report_xls($project,\@final_data);
	exit(0);
}
else {
export_data::print($project,$cgi,\@final_data,$type_label);
}

exit(0);

sub update_pedigree {
	my ($project,$data,$pedigree) = @_;
	foreach my $d (@$data) {
		
		foreach my $p (@{$d->{patient_name}}){
			my $ped = $pedigree->{$p};
			push(@{$d->{pedigree_type}},$ped->{type}.$ped->{status});
			push(@{$d->{pedigree_status}},$ped->{status});
			push(@{$d->{pedigree_fam}},$ped->{fam});
			#$d->{pedigree}->{status}->{$p} = $ped->{status};
			
		}
	}

}


sub update_genes_names {
	my ($project,$data) = @_;
	foreach my $d (@$data) {
		my $z=0;
		my %name2;
			foreach my $c (@{$d->{tab_consequences}}){
				my $a = $c->{gene};
					$a =~ /(.*) \((.*)\)/;
				my $ensg = $1;
				$name2{$1} = $2;
			} 
			$z++;
		$d->{external_names} = join(";",map{$name2{$_}} @{$d->{genes}});
		$d->{genes_name} = join(";",@{$d->{genes}});
	}
}

sub by_patients {
	my ($project,$data) = @_;
	my %patients;
	foreach my $d (@$data) {
		map {$patients{$_} ++}  @{ $d->{patient_name} };
	}
	
	my $out;
	foreach my $name (sort{$a cmp $b} keys %patients){
		my $item;
		$item->{name} = $name;
		push(@$out,$item);
	}
	export_data::print($project,$cgi,$out);
	exit(0);
}


sub update_valid {
	my ($project, $data ,$user) = @_;
	
	my $captures = $project->getCaptures();
	
	my $capture = @$captures[-1];
	return unless $capture->validation_db;
	my $vquery = validationQuery->new(dbh=>$project->buffer->dbh,capture_name=>$capture->validation_db);

	return unless $vquery->exists_db();
	
	foreach my $d (@$data) {
		my $valid =0;
		my $invalid = 0;
		my $notseq=0;
		my $ho = 0;
		my $he = 0 ;
		my $globalValid = 0;
		my $nbpatients = scalar (@{ $d->{patient_name} });
		my $id;
		
		 my $validation_vid = $vquery->getVariationByGenBoId(id=>$d->{id});
		 
		 #warn $validation_vid." ::  ".$d->{id};
		if ($validation_vid){
			foreach my $patient_name (@{ $d->{patient_name} }){
				my ($value_valid,$validation_sanger) = $vquery->getValidations(id=>$validation_vid,project=>$project->name(),sample=>$patient_name);
				$value_valid = $validation_sanger if $validation_sanger;
				$value_valid += 100 if $validation_sanger;
				$d->{ "valid!" .$patient_name} = $value_valid if $value_valid;
				$d->{valid} = "1" if $value_valid;
			}
		}
	}

}


