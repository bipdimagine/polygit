package filter_cache;
use Set::Intersection;
use strict;
require Exporter;
our @ISA = qw/Exporter/;

our @EXPORT = qw/get_genes_all delete_variations delete_intersection new_intersection keep_intersection get_variations_cache/;



sub get_genes_all {
	my ($project,$flag_dbsnp) = @_;
	 my $file_cache_genes = $project->getCacheGenesFile();
	 my $md5 = `md5sum $file_cache_genes`;
	 my $file2 = $project->getCacheDbsnpGenesFile();
	 
	
	 if ($flag_dbsnp && -e $file2 && -e $file_cache_genes){
	 	my $data = $project->buffer->getStore($file2);
	 	
	 	
	 	chomp($md5);
		if ($md5 eq $data->{md5}){
			warn "dbsnp cache file !!!!";
			delete $data->{md5};
	 		return $data;
		}
	 }
	 if (-e $file_cache_genes){
	 	my $data = $project->buffer->getStore($file_cache_genes);
	 	warn "file cache !!!!!!!!!!!!!!";
	 	return $data;
	 }
	die();
	my $id =  GenBoStorable::getStoreId( $project->buffer->dbh, $project->id, $project->id,"genes" );
	return unless $id;
	return (GenBoStorable::getStore( $project->buffer->dbh, $id ));
} 


sub delete_variations {
	my ($filters,$data) = @_;
	my %delete_var;
	my $current_variations = $data->{all_variations};
	foreach my $filter_name (@$filters) {
			next unless (exists $data->{data}->{$filter_name});
			map {$delete_var{$_}++}  @{$data->{data}->{$filter_name}};
			delete $data->{data}->{$filter_name};
	}

		return delete_intersection($data,\%delete_var);
}

sub delete_intersection {
	my ($data,$hash_delete,$debug) =@_;
	my @vars = diff($data->{all_variations},$hash_delete);
	$data->{all_variations} = \@vars;	
	return unless (@vars);
	refresh_data($data);
	return 1;
}

sub refresh_data {
	my ($data) = @_;
	my $vars = $data->{all_variations};
	
	my %hash;
	 @hash{@$vars} = ();
	 
	foreach my $name (keys %{$data->{data}}) {
		my $debug;
		my @res = new_intersection(\%hash,$data->{data}->{$name});
		#my @res = get_intersection($vars,$data->{data}->{$name});
		
		if (@res) {
			$data->{data}->{$name} = \@res;
		}
		else {
			delete $data->{data}->{$name} ;
		}
	}
	
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
	 
	my %n = map { $_ => undef } grep exists $hin->{$_}, @$left;
  return (keys %n);
}	

sub keep_intersection {
	my ($data,$array,$debug) =@_;
	
	return  unless @$array;
	

	my (@vars) = get_intersection($data->{all_variations},$array);
	
	return unless (@vars);
	$data->{all_variations} = \@vars;	
	#$data->{h_all_variations} = $hvars;
	refresh_data($data);
	return 1;
}

1;