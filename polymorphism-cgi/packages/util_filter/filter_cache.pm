package filter_cache;
use Set::Intersection;
use strict;
use Data::Dumper;
require Exporter;
use KyotoCabinet;
use Storable qw/freeze thaw nfreeze nstore_fd nstore retrieve/;
use Math::Combinatorics;
use List::MoreUtils qw{uniq};
use Carp;
our @ISA = qw/Exporter/;

our @EXPORT = qw/get_genes_all get_combination_param_for_cache get_data delete_variations  new_intersection  get_variations_cache/;



sub get_genes_all {
	return get_data(@_);
}
sub get_data {
	my ($project,$key) = @_;
	my $data;

	return get_kyoto_raw_cache($project,undef);
}
sub get_kyoto_delete_genes_cache {
	my ($project,$key) = @_;
	
	my $file = "/data-xfs/sequencing/ngs/Polycache/www/".$project->name().".kct";
	#get_kyoto_raw_cache();
	my $hcache;
	 tie(%{$hcache}, 'KyotoCabinet::DB', $file , KyotoCabinet::DB::OREADER | KyotoCabinet::DB::ONOLOCK) || die($file);
	
	 my $d = $hcache->{$key};
	  return thaw $d if $d;
	  return;
	
}

sub get_kyoto_param_cache{
	my ($project,$key) = @_;
	confess();
	my $file = "/data-xfs/sequencing/ngs/Polycache/www/".$project->name().".kct";
	#get_kyoto_raw_cache();
	get_kyoto_raw_cache($project) unless -e $file;
	warn "kyoto";
	  my $hcache;
	 tie(%{$hcache}, 'KyotoCabinet::DB', $file , KyotoCabinet::DB::OREADER | KyotoCabinet::DB::ONOLOCK) || die($file);

	 my $d = $hcache->{$key};
	warn $key if $d;
	  return thaw $d if $d;
	  return;
	
}

sub get_kyoto_raw_cache {
	my ($project) = @_;
	
	# my $file_cache_genes = $project->getCacheGenesFile().".mask";
 	my $file_cache_genes = $project->getCacheGenesFile();
	 if (-e $file_cache_genes){
	 	
	 	my $data = $project->buffer->getStore($file_cache_genes);
	 	return $data;
	 }
	confess($file_cache_genes);
} 

sub get_combination_param_for_cache {
	my @param_for_cache = ("dbsnp","evs","1000genomes","silent","dbsnp_1p","evs_1p","1000genomes_1p","intergenic");
@param_for_cache = sort {$a cmp $b} @param_for_cache;
my @all_combination;
for (my $i =@param_for_cache; $i > 2 ;$i-- ){
 push(@all_combination,combine($i,@param_for_cache));
}
return \@all_combination;
}





sub delete_variations {
	my ($filters,$data,$first) = @_;
	my %delete_var;
	my $current_variations = $data->{all_variations};
	foreach my $filter_name (@$filters) {
			
			next unless (exists $data->{data}->{$filter_name});
		#	warn "coucou" unless %{$data->{hash_all_variations}};
			last unless %{$data->{hash_all_variations}};
			map {delete  $data->{hash_all_variations}->{$_}}  @{$data->{data}->{$filter_name}};
			delete $data->{data}->{$filter_name};
			$first->{debug}++;
	}
	
	return scalar( keys %{$data->{hash_all_variations}});
}




sub delete_intersection {
	my ($data,$hash_delete,$debug) =@_;
	my @vars = diff($data->{all_variations},$hash_delete);
	
	$data->{all_variations} = \@vars;	
	
	return unless (@vars);
	refresh_data($data);
	return scalar  @vars;
}



	
sub diff {
  my ($left,$hin) = @_;
  my %n;
  %n = map { $_ => undef } grep !(exists $hin->{$_}), @$left;
  return (keys %n);
  my @t = (keys %n);
  return \@t;
}

sub new_intersection {
	 my ($hin,$left) = @_;
		my @toto=  grep exists $hin->{$_}, @$left;
  	return (\@toto);
}	


1;