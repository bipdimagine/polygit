package root_filter;
use strict;
use Set::Intersection;
use Moose;
#use MooseX::Method::Signatures;
use Data::Dumper;
use Carp;

has 'data' => (
	is =>'rw',
	isa =>"HashRef",
	required => 1,
);

has 'patients' => (
	is =>'rw',
	isa =>"ArrayRef",
	required => 1,
);
has 'pedigree' => (
	is =>'rw',
	isa =>"ArrayRef",
);

has 'individual_pedigree' => (
	is =>'rw',
	isa =>"HashRef",
);
has 'attic' => (
	is =>'rw',
	isa =>"Str",
	default=>"",
);
has 'heho_patients'  => (
	is =>'rw',
	isa =>"HashRef",
);


has 'chr'  => (
	is =>'rw',
	isa =>"Str",
	trigger   => \&set_chr,
);

has 'gid'  => (
	is =>'rw',
	isa =>"Str",
	trigger => \&set_gid,
);

has 'hashvar' => (
	is =>'rw',
	isa =>"HashRef",
	weak_ref =>1,
	
);

has 'gene' =>(
	is =>'rw',
	isa =>"HashRef",
	weak_ref =>1,
);

has 'gene_data' =>(
	is =>'rw',
	isa =>"HashRef",
	weak_ref =>1,
);

has 'attic' =>(
	is =>'rw',
	isa =>"HashRef",
);
has 'special_filters' =>(
	is =>'ro',
	isa =>'HashRef',
	default => \&set_special_filter,
);
has 'need_refresh'=>(
	is =>'rw',
	isa =>'HashRef',
	#default => \&set_special_filter,
);
 
sub set_special_filter {
	my $self = shift;
	$self->{specialfilter} = {
		"essential_splicing" => 1,
		"splicing" => 1,
		"phase" => 1,	
		"coding" => 1,
		"silent" => 1,
		"phase"=> 1,
		"non-frameshift" => 1,
		"frameshift" => 1,
		"stop" => 1
	};
	
}

sub set_gid {
	my $self =shift;

	
	$self->gene($self->data->{$self->chr}->{genes}->{$self->gid});
	
	unless (exists  $self->gene->{hash_all_variations}){
	
		my $vars =$self->gene->{all_variations};
		my %hash;
		@hash{@$vars} = ();
		$self->gene->{hash_all_variations} = \%hash ;
	}
	$self->hashvar($self->gene->{hash_all_variations});
	$self->gene_data($self->gene->{data});
	
#	my @filter = keys %{$self->{gene_data}};
#	my %filters;
#	@filters{@filter} = ();
#	$self->{need_refresh} = \%filters;
}


sub set_chr {
	my $self = shift;
	 my $hash_heho_patients;
	 my $data = $self->data->{$self->chr}->{patients};
	 $self->heho_patients({});
#	foreach my $patient (@{$self->patients}) { 
#		$data->{$patient}->{homozygote} = [] unless exists  $data->{$patient}->{homozygote};
#		my %hash_ho_patients;	
#		@hash_ho_patients{@{$data->{$patient}->{homozygote}}} = ();
#		$self->heho_patients->{$patient}->{homozygote} = \%hash_ho_patients;
#	}
	
}


sub make_he_ho () {
	my ($self,$lists,$type) = @_;
	warn "coucou";
	die();
	#my ($lists,$type,$data,$hash_heho_patients) = @_;
	my $data = $self->gene_data();
	foreach my $p (@$lists) {
		$data->{$type."_".$p} = [];
		
		foreach my $t (@{$data->{$p}}){
			
			if ($type eq "homozygote"){
				push(@{$data->{$type."_".$p}},$t)  if exists $self->heho_patients->{$p}->{homozygote}->{$t};
			}
			else {
				push(@{$data->{$type."_".$p}},$t)  unless exists $self->heho_patients->{$p}->{homozygote}->{$t};
			}
		}
	}
	
	
}



sub refresh_data_filter {
	my ($self,$filter) =@_;
#	warn $filter;
	my $res = $self->new_intersection($self->{hashvar},$self->{gene_data}->{$filter});
	$self->{gene_data}->{$filter} = $res; 
	#$self->put($filter,$res);
		
	 $self->{need_refresh}->{$filter} = undef;
	
}



sub refresh_data {
	my $self =shift;
	#my $data = $self->get_root_gene();
	
	my @filters = grep {$_ !~ /zygote_/} keys %{$self->gene_data};
	foreach my $name (@filters) {
		my $res = $self->new_intersection($self->hashvar,$self->get($name));
		$self->put($name,$res);
	
	}
}

sub refresh_data2 {
	my $self =shift;
	$self->{need_refresh} = {};
	return;
	#my $data = $self->get_root_gene();
	my @filter = keys %{$self->gene_data};
	my %filters;
	@filters{@filter} = ();
	$self->{need_refresh} = \%filters;

#	my @filters = grep {$_ !~ /zygote_/} keys %{$self->gene_data};
#	foreach my $name (@filters) {
#		my $res = $self->new_intersection($self->hashvar,$self->get($name));
#		$self->put($name,$res);
#	
#	}
}

sub new_intersection  {
	my ($self,$hin,$left) = @_;
	return [] unless %$hin;
	my @toto=  grep exists $hin->{$_}, @$left;
	return (\@toto);
	
  
}	



sub delete_variations {
	my ($self,$filters,$first) = @_;
	return;
	foreach my $filter_name (@$filters) {
			#last if $self->isEmpty;
			map {delete  $self->hashvar->{$_}}  @{$self->get($filter_name)};
			 $self->delete_key($filter_name);
			
	}
	
	return;
}


sub delete_special_filters {
	my ($self,$filters) = @_;
	
	return 1 unless scalar(@$filters);
	
						  
						  
	my %hfilters;
	 @hfilters{@$filters} = undef;
	
	my %delete_var;
	my %all_vars;
	
	my @ff = keys %{$self->special_filters};
	foreach my $filter_name (@ff) {
		next if exists $hfilters{$filter_name};
		
		 @all_vars{@{$self->get($filter_name)}} =undef;
	}
	
	
	foreach my $filter_name (@$filters) {
			#next unless (exists $data->{$filter_name});
			last if $self->isEmpty();
			map {delete  $self->hashvar->{$_}} grep {! exists $all_vars{$_}} @{$self->get($filter_name)};
			 $self->delete_key($filter_name);
	}
	return;

	
}

sub change_hashvar {
	my ($self,$ref) = @_;
	$self->hashvar($ref);
}

sub delete {
	my $self =shift;
	delete $self->data->{$self->chr}->{genes}->{$self->gid};
}

sub isEmpty {
	my $self = shift;
	unless (%{$self->hashvar}) {
		#$self->delete_gene();
		return 1 ;
	}
	return undef;
}


sub diseases {
	my ($self,$disease) = @_;
	
}

sub get {
	my ($self,$filter) = @_;
	
	if (exists  $self->{gene_data}->{$filter}){
		
		#$self->refresh_data_filter($filter) unless exists $self->{need_refresh}->{$filter};
		return $self->{gene_data}->{$filter} ;
	}
	return [];
}
sub put {
	my ($self,$filter,$data) = @_;
	unless (@$data){
		$self->delete_key($filter);
	}
	$self->gene_data->{$filter} = $data;
}

sub has_key {
	my ($self,$filter) = @_;
	return scalar(@{$self->gene_data->{$filter}}) > 0;
}

sub delete_key {
	my ($self,$filter) = @_;
	delete $self->gene_data->{$filter};
	
}

sub prepare_sift_polyphen {
	my ($self,$list_polyphen,$list_sift) = shift;
	my %var_tmp_polyphen;
	foreach my $l (@$list_polyphen){
		
		my $vars = $self->get($l);
		foreach my $vv (@$vars){
							$var_tmp_polyphen{$vv} = 1;
				}
		}	 
			
		my $all_vars;
		my @var_tmp_sift;
		foreach my $l (@$list_sift){
					#next unless exists $self->gene_data->{$l};
					push(@var_tmp_sift, @{$self->get($l)});
			} 
		$self->put("polyphen_sift",$self->new_intersection(\%var_tmp_polyphen,\@var_tmp_sift));
}

sub set_attic {
	my ($self,$attic_string) = @_;
	my %hattic;
	@hattic{split(" ",$attic_string)} = ();
	$self->attic(\%hattic);
}

sub in_the_attic  {
	my ($self,$attic_string) = @_;
	

	
	my $shadow_variations;
	my $vars =[];
	my $vars_homo = [];
	foreach my $p (@{$self->patients}) {
		next if exists $self->attic->{$p};
		push(@$vars_homo,@{$self->get("homozygote_".$p)});
		push(@$vars,@{$self->get($p)});
				
		}
		
		unless (scalar(@$vars)){
			$self->change_hashvar({});
			#$self->delete();
			return;
		}
			
			my $common_variations;
			my %hash;
			@hash{@$vars} = ();
			$self->change_hashvar(\%hash);
			my %hash2;
			@hash2{@$vars_homo} = ();
			$self->put("homozygote",[keys %hash2]);
		
#			}
			foreach my $p (keys %{$self->attic}) {
				$self->delete_key($p);
				$self->delete_key("homozygote_".$p);
				$self->delete_key("heterozygote_".$p);

			}
}

1;