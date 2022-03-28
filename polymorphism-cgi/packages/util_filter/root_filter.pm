package root_filter;
use strict;
use Set::Intersection;
use Moose;
use Data::Dumper;
use Storable qw(store retrieve freeze thaw);
use KyotoCabinet;
#use MooseX::Method::Signatures;
   use Carp qw(cluck longmess shortmess);

has 'data' => (
	is =>'rw',
	required => 1,
);
has 'type' => (
	is =>'ro',
	required => 1,
);
has 'project' => (
	is =>'rw',
	required => 1,
);

has 'patients' => (
	is =>'rw',
	required => 1,
);

has 'pedigree' => (
	is =>'rw',
	default=>sub{ return [];},
);

has 'individual_pedigree' => (
	is =>'rw',
);
has somatic_file => (
	is =>'rw',
	default => "",
);

has 'patients_groups' => (
	is =>'ro',
	lazy =>1,
	default =>
		sub {
		my $self = shift;
		my $file = $self->somatic_file();
		return {} unless -e $file;
		open(FILE,$file);
		my $al;
		while (my $line = <FILE>){
			chomp($line);
			next if $line eq "";
			my $p;
			($p->{group},$p->{name},$p->{status}) = split(" ",$line);
			$al->{$p->{name}} = $p;
			#push(@al,$p);
		
		}
	

	return $al;
	
		}
);


has 'groups' => (
	is =>'ro',
	lazy =>1,
	default =>
		sub {
		my $self = shift;
		my $groups;
		foreach my $p (values %{$self->patients_groups}){
			push(@{$groups->{$p->{group}}->{samples}},$p);  
				if ($p->{status} == 1){
					push(@{$groups->{$p->{group}}->{germinal}},$p->{name});  
				}
				elsif   ($p->{status} == 2){
					push(@{$groups->{$p->{group}}->{somatic}},$p->{name});  
				}
				else {
					push(@{$groups->{$p->{group}}->{unknown}},$p->{name});  
				}
				$groups->{$p->{group}}->{name} = $p->{group};
		}
		return $groups;
	
		}
);


has 'individual_group' => (
	is =>'rw',
);
has 'attic' => (
	is =>'rw',
	isa =>"Str",
	default=>"",
);
has 'heho_patients'  => (
	is =>'rw',
);


has 'chr'  => (
	is =>'rw',
	trigger   => \&set_chr,
);

has 'gid'  => (
	is =>'rw',
	trigger => \&set_gid,
);

has 'hashvar' => (
	is =>'rw',
	
);

has 'gene' =>(
	is =>'rw',
);

has 'gene_data' =>(
	is =>'rw',
);

has 'attic' =>(
	is =>'rw',
);
has 'special_filters' =>(
	is =>'ro',
	default => \&set_special_filter,
);
has 'need_refresh'=>(
	is =>'rw',
	#default => \&set_special_filter,
);



#
#sub special_filters {
#	my $self = shift;
#	$self->set_special_filter() unless exists $self->{specialfilter};
#	return $self->{specialfilter};
#}


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
	my ($self,$array) =@_;


	$self->{need_refresh} = {};
	if ($self->type == 1){
		$self->data->{$self->chr}->{genes}->{$self->{gid}} = $self->get_data_from_kyoto($self->{gid});
	}
		$self->{gene} = $self->{data}->{$self->chr}->{genes}->{$self->{gid}};
		
	
	unless 	(keys %{$self->{gene}}){
		if ($self->chr eq 'X' && $self->{gid} =~ /_Y/){
			return;
		}
		warn $self->{gid}." ".$self->chr;
		warn Dumper $self->data->{$self->chr}->{genes}->{$self->{gid}};
		 confess();
	} 
	unless (exists  $self->{gene}->{hash_all_variations}){
#		die();
		my $vars =$self->{gene}->{all_variations};
	#	warn Dumper($self->{data}->{$self->chr}->{genes}->{$self->{gid}});
		confess($self->{gid});
		my %hash;
		@hash{@$vars} = ();
		$self->{gene}->{hash_all_variations} = \%hash ;
	}
	$self->{hashvar} = $self->{gene}->{hash_all_variations};
	
	$self->{gene_data} = $self->{gene}->{data};
	$self->{first} = 1;
	return;
}

sub _heho {
	my ($self,$name,$type,$data) =@_;
	if ($data){
		$self->put($type."_$name",$data);
		return 1;
	}
	else {
		return $self->get($type."_".$name,$data);
	}
	
}
sub homozygote {
	my ($self,$name,$data) =@_;
	return $self->_heho($name,"homozygote",$data);
}
sub heterozygote {
	my ($self,$name,$data) =@_;
	return $self->_heho($name,"heterozygote",$data);
}

sub set_chr {
	my $self = shift;
	 my $hash_heho_patients;
	 
	 my $data = $self->{data}->{$self->{chr}}->{patients};
	 $self->{heho_patients} = {};
	foreach my $patient (@{$self->patients}) { 
		my $t = $self->get_db_homo($patient);
		unless ($t){
		$data->{$patient}->{homozygote} = [] unless exists  $data->{$patient}->{homozygote};
		my %hash_ho_patients;	
		@hash_ho_patients{@{$data->{$patient}->{homozygote}}} = ();
		$self->{heho_patients}->{$patient}->{homozygote} = \%hash_ho_patients;
		}
	}

		if ($self->type ==1) {
			$self->set_kyoto_gene();
		}
}

sub set_kyoto_gene{
		my ($self,$chr) = @_;
		if (exists $self->{db_genes}){
			$self->{db_genes}->close();
		}
		my $db = new KyotoCabinet::DB;
		my $file =  $self->project->getCacheGenesKyotoFile($self->{chr});
		if (!$db->open($file, $db->OREADER | $db->ONOLOCK)) {
     		printf STDERR ("open error: %s\n", $db->error);
 		}
 		
		$self->{db_genes} = $db;
} 
sub get_data_from_kyoto {
	my ($self,$id) = @_;
	return thaw $self->{db_genes}->get($id);
	
}

sub make_he_ho () {
	my ($self,$lists,$type) = @_;
	#my ($lists,$type,$data,$hash_heho_patients) = @_;

	my $data = $self->{gene_data};
	#$lists = ["homozygote","heterozygote"];
	
	foreach my $p (@$lists) {
		$data->{$type."_".$p} = [];
		foreach my $t (@{$self->get($p)}){
			if ($type eq "homozygote"){
				push(@{$data->{$type."_".$p}},$t)  if exists $self->heho_patients->{$p}->{homozygote}->{$t};
				#push(@{$data->{$type."_".$p}},$t)  if exists $h->{$t};
			}
			else {
				push(@{$data->{$type."_".$p}},$t)  unless exists $self->heho_patients->{$p}->{homozygote}->{$t};
				#push(@{$data->{$type."_".$p}},$t)  unless exists $h->{$t};
			}
		}
	}
	
	
}

sub get_db_homo{
	my ( $self,$name) = @_;
	my $file = $self->project->getCacheDir()."/$name.homo.kct";
	return undef unless -e $file;
	return $self->{heho_patients}->{$name}->{homozygote}  if exists $self->{heho_patients}->{$name}->{homozygote};
	$self->{heho_patients}->{$name} = {};
	if (-e $file){
			tie(%{$self->{heho_patients}->{$name}->{homozygote} }, 'KyotoCabinet::DB', $file , KyotoCabinet::DB::ONOLOCK | KyotoCabinet::DB::OREADER) || die($file);
	}
	else {
		die($file);
		#my %hash_ho_patients;	
		#@hash_ho_patients{@{$data->{$patient}->{homozygote}}} = ();
		#$self->{heho_patients}->{$name}->{homozygote} = \%hash_ho_patients;
	}
	
	return $self->{heho_patients}->{$name}->{homozygote} ;
}


sub refresh_data_filter {
	my ($self,$filter) =@_;
	
		if ($filter =~ /gote_/ ) {
		
			# this is for homo or heterozygote 	
			my ($type,$name) = split("_",$filter);
			#	cluck "$type $name";
			$self->make_he_ho([$name],$type);  
		}
		my $res = $self->new_intersection($self->{hashvar},$self->{gene_data}->{$filter});
		if ($res) {
			$self->{gene_data}->{$filter} = $res;
			#$self->put($filter,$res);
		}
		else {
			delete $self->{gene_data}->{$filter};
		}
		#warn Dumper  $self->need_refresh." $filter";
		$self->{need_refresh}->{$filter} =();
}

sub refresh_data {
	my $self =shift;
	delete $self->{first};
	$self->{need_refresh} = {};
	return;
	my @filter = keys %{$self->{gene_data}};
	my %filters;
	@filters{@filter} = ();
	$self->{need_refresh} = \%filters;
}

sub new_intersection  {
	my ($self,$hin,$left) = @_;
	return [] unless %$hin;
	my @toto=  grep exists $hin->{$_}, @$left;
	return (\@toto);
	
  
}	


sub delete_variations {
	my ($self,$filters,$first) = @_;
	my $dejavu;
	my $nb;
		
	foreach my $filter_name (@$filters) {
		next unless (exists $self->{gene_data}->{$filter_name});
		

		last unless %{$self->{hashvar}};
			map {delete  $self->{hashvar}->{$_}}  @{$self->get($filter_name)};
			$self->{debug}++;
			delete $self->{gene_data}->{$filter_name};
			
			
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
		next unless exists $self->{gene_data}->{$filter_name};
		 @all_vars{@{$self->{gene_data}->{$filter_name}}} =undef;
	}
	
	
	foreach my $filter_name (@$filters) {
			#next unless (exists $data->{$filter_name});
			last if $self->isEmpty();
			next unless exists $self->{gene_data}->{$filter_name};
			map {delete  $self->{hashvar}->{$_}} grep {! exists $all_vars{$_}} @{$self->{gene_data}->{$filter_name}};
			
			$self->delete_key($filter_name);
	}
	return;

	
}

sub change_hashvar {
	my ($self,$ref) = @_;
	$self->{gene}->{hash_all_variations} = $ref;
	$self->{hashvar} = $self->{gene}->{hash_all_variations};
	$self->refresh_data();
}

sub delete {
	my ($self,$debug) =@_;
	delete $self->{data}->{$self->{chr}}->{genes}->{$self->{gid}};
}

sub isEmpty {
	my $self = shift;
	unless (%{$self->{hashvar}}) {
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
	return [] unless %{$self->{hashvar}};
	
	if (exists  $self->{gene_data}->{$filter}){
		return $self->{gene_data}->{$filter} if exists $self->{first}; 
		$self->refresh_data_filter($filter) unless exists $self->{need_refresh}->{$filter};
		return $self->{gene_data}->{$filter} ;
	}
	return [];
}

sub final_stats {
	my ($self,$type,$stats,$filter) =@_;
	warn $filter;
	#my $res = $self->new_intersection($self->{hashvar},$self->{gene_data}->{$filter});
	#$self->put($filter,$res);
	#warn Dumper  $self->need_refresh." $filter";
	delete $self->{need_refresh}->{$filter};
	
}
sub put {
	my ($self,$filter,$data) = @_;
	unless (@$data){
		$self->delete_key($filter);
	}
	$self->{gene_data}->{$filter} = $data;
}

sub has_key {
	my ($self,$filter) = @_;
	return scalar(@{$self->{gene_data}->{$filter}}) > 0;
}

sub delete_key {
	my ($self,$filter) = @_;
	delete $self->{gene_data}->{$filter};
	
}
sub delete_patient {
	my ($self,$name) = @_;
	delete $self->{gene_data}->{$name};
	delete $self->{gene_data}->{"homozygote_".$name};
	delete $self->{gene_data}->{"heterozygote_".$name};
}

sub prepare_sift_polyphen {
	my ($self,$list_polyphen,$list_sift) = @_;
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
#
#sub attic {
#	my $self = shift;
#	$self->set_attic() unless exists $self->{attic};
#	return $self->{attic};
#}
#


sub set_attic {
	my ($self,$attic_string) = @_;
	my %hattic;
	@hattic{split(" ",$attic_string)} = ();
	$self->{attic} = \%hattic;
}

sub in_the_attic  {
	my ($self,$hash) = @_;
	my $filter_patients = $self->attic;
	if ($hash) {
		die();
		$filter_patients = $hash;
	} 
	my $shadow_variations;
	my $vars =[];
	my $vars_homo = [];
	foreach my $p (@{$self->patients}) {
		next if exists $filter_patients->{$p};
		#push(@$vars_homo,@{$self->get("homozygote_".$p)});
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
			#my %hash2;
			#@hash2{@$vars_homo} = ();
			#$self->put("homozygote",[keys %hash2]);
#			}
			foreach my $p (keys %{$self->attic}) {
				$self->delete_patient($p);
			}
}


sub exclude_genes {
	my ($self,$names) = @_;
	foreach my $name (@$names){
		next if scalar (@{$self->get($name)}) == 0;
		 $self->change_hashvar({});
		 last;
	}
	
}

sub exclude_variations{
	my ($self,$names) = @_; 
	$self->delete_variations($names);
}

sub and_genes {
	my ($self,$names) = @_; 
	my $nb =0;
	foreach my $name (@$names){
		next if scalar (@{$self->get($name)}) == 0;
		$nb++;
	}
	$self->change_hashvar({}) if $nb ne scalar(@$names);
}

sub and_variations {
	my ($self,$names) = @_;
	my %var_ids;
		foreach my $name (@$names) {
				next if scalar (@{$self->get($name)}) == 0;
				#next unless exists $data->{$name};
				foreach my $vid (@{$self->get($name)}){
						$var_ids{$vid} ++;
				}
			}
			my $nb = scalar(@$names);
			my @kk = keys %var_ids;
			foreach my $k (@kk){
				next if $var_ids{$k}== $nb;
				delete $var_ids{$k};
			}
			
		$self->change_hashvar(\%var_ids);		
	
}
sub sam {
	my ($self,$name) = @_;
	return $self->{sam}->{$name} if exists $self->{sam}->{$name};
	my $patient = $self->project->getPatient($name);
	my $bam = $patient->getBamFile();
	$self->{sam}->{$name} = Bio::DB::Sam->new(-bam  =>$bam,
                           #  -fasta=>"/data-xfs/public-data/HG19/genome/fasta/all.fa",
                             );
    return $self->{sam}->{$name};
}

sub get_sample_variations {
	my ($self,$name) = @_;
	return $self->{kct}->{$name} if exists $self->{kct}->{$name};
	my $kyoto_filter_patients =  $self->project->getCachePatientKyotoFile($name);
	tie(%{$self->{kct}->{$name}}, 'KyotoCabinet::DB', $kyoto_filter_patients ,KyotoCabinet::DB::ONOLOCK | KyotoCabinet::DB::OREADER ) || die($kyoto_filter_patients);
	return  $self->{kct}->{$name};
}


1;