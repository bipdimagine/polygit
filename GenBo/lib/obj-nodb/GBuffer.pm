package GBuffer;
use Moo;
use Carp;
use Data::Dumper;
use Config::Std;
use GenBoProject; 
use GenBoProjectCache;
use GenBo;
use QueryMoosePolyProjectNgs;
use QueryPolyPanel;
use QueryHgmd;
use QueryPolyPhenotype;
use QueryClinvarPathogenic; 
use Storable qw(store retrieve freeze thaw fd_retrieve);
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use QueryMoosePolyProject; 
use QueryFilter;
use validationQuery;
use connect;
use FindBin;
use Parallel::ForkManager;
use Bio::DB::HTS::Tabix;
use packages::SmithWaterman;
use packages::NeedlemanWunsch;
use Storable qw(dclone);
use JSON::XS;
use Sys::Hostname;
use POSIX qw(strftime);
use DateTime;
use List::Util qw[min max];
use Cwd qw(abs_path);
use Auth::GoogleAuth;

#use Sereal qw(sereal_encode_with_object sereal_decode_with_object);


has maxProc => (
	is 		=> 'rw',
	default	=> '1',
);

has type_db => (
	is 		=> 'ro',
	default	=> 'PolyprojectNGS',
);

has verbose => (
	is		=> 'rw',
	default	=> sub {
		my $self = shift;
		return 1;
		my $user ="";
		 $user .= getlogin();
		 
		if ($user =~ /apache/) { return 1; }
		return 0;
	},
);

has vmtouch => (
	is		=> 'rw',
	default	=> sub {
		if (hostname =~ /img/ ){
			return 1;
		}
		return undef;
	},
);

has genbo_dir =>(
		is		=> 'ro',
		lazy	=> 1,
	default	=> sub {
			my $dir = $INC{"GBuffer.pm"};
			
			 $dir =~ s/GBuffer\.pm//;
			return $dir;
	}

);

has google_auth_issuer =>(
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return 'POLYWEB';
	}
);

has google_auth_key_id =>(
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return 'STAFF';
	}
);

has google_auth =>(
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		my $date_now = DateTime->now;
		my $auth = Auth::GoogleAuth->new({issuer => $self->google_auth_issuer(), key_id => $self->google_auth_key_id()});
		return $auth;
	}
);

sub use_otp_for_login {
	my ($self, $login) = @_;
	my $h_otp = $self->getQuery->getSecretOtpKeyForUserId($login);
	return $h_otp->{'uKey'};
	return;
}

sub google_auth_secret_pwd {
	my ($self, $login) = @_;
	my $dbh = $self->getQuery();
	return if not $self->use_otp_for_login($login);
	my $h_otp = $dbh->getSecretOtpKeyForUserId($login);
	return $h_otp->{'Key'};
	confess();
}

has google_auth_qr_code =>(
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		$self->google_auth->secret32( $self->google_auth_secret_pwd() );
		return $self->google_auth->qr_code;
	}
);

sub hasHgmdAccess {
	my ($self, $user) = @_;
	return 1 if ($self->queryHgmd->getHGMD($user) == 1);
	return;
}

has config => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift; 
		my $dir = $self->genbo_dir;
		my $filename = $dir."genbo.cfg";
		
		my $filename2 = $dir."genbo-vector-filter.cfg";
		confess($filename) unless -e $filename;
		confess($filename2) unless -e $filename2;
	#	warn $filename2;
	#	my @t = `cat $filename2`;
	#	warn Dumper @t;
		read_config $filename => my %config1;
		read_config $filename2 => my %config2;
		my $hConfig;
		foreach my $k (keys %config1){
			foreach my $k2 (keys %{$config1{$k}}) {
				$hConfig->{$k}->{$k2} = $config1{$k}{$k2};
			}
		}
		foreach my $k (keys %config2){
			foreach my $k2 (keys %{$config2{$k}}) {
				$hConfig->{$k}->{$k2} = $config2{$k}{$k2};
			}
		}
		return $hConfig;
	},
);
has public_data_versions => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift; 
		return [sort{$a <=> $b} keys %{$self->public_data}];
	}
	);
has config_directory => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift; 
		my $dir = $self->genbo_dir."../../../GenBoConfig/";
		confess("problem with config directory :".$dir ) unless -e $dir;
		return $dir; 
	}
	);
	
has polybtf_infos => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift; 
		my $dir = $self->config_directory."/public_data/";
		my $filename = $dir."public-data.version.cfg";
		confess($filename) unless -e $filename;
		read_config $filename => my %config;
		my @list_versions = (sort{$a <=> $b} keys %config);
		my $last_version = $list_versions[-1];
		my $hash;
		$hash->{last_version} = $list_versions[-1];
		my ($date, $since);
		foreach my $cat (sort keys %{$config{$last_version}}) {
			push(@{$hash->{news}}, ucfirst($cat).': '.$config{$last_version}->{$cat});
			my $dir = $self->public_data_root().'/'.$self->build().'/'.$cat.'/'.$config{$last_version}->{$cat}.'/lmdb/';
			my $version_json_file = $dir.'/version.json';
			next unless -e $version_json_file;
			open( JSON, $version_json_file );
			my $desc = decode_json <JSON>;
			close(JSON);
			$date = $desc->{date};
			my $t = (stat $dir )[9];
			my $t1 = localtime();
			$since = int((time - $t) / 86400);
			last;
		}
		$hash->{date_release} = $date;
		$hash->{date_now} = strftime "%d/%m/%Y", localtime();
		$hash->{date_last_days} = $since;
		return $hash;
	},
);
	
has public_data => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift; 
		my $dir = $self->config_directory."/public_data/";
		my $filename = $dir."public-data.cfg";
		confess($filename) unless -e $filename;
		read_config $filename => my %config1;
		my $filename2 = $dir."public-data.version.cfg";
		confess($filename2) unless -e $filename2;
		read_config $filename2 => my %config2;
		confess($filename2) unless -e $filename2;
		my $public_data;
		my $previous;
		foreach my $v (sort{$a <=> $b} keys %config2){
			if ($previous){
				$public_data->{$v} = dclone $public_data->{$previous};
				foreach my $k (keys %{$config2{$v}}){
					$public_data->{$v}->{$k}->{version} = $config2{$v}{$k};
				}
			}
			else {
				foreach my $k (keys %{$config2{$v}}){
					$public_data->{$v}->{$k}->{version} = $config2{$v}{$k};
				}
			}
			$previous = $v;
			
			foreach my $name (keys %{$public_data->{$v}}) {
				die($name) unless $config1{$name};
				$public_data->{$v}->{$name}->{config} = dclone $config1{$name} unless exists $public_data->{$v}->{$name}->{config};
				#warn  $public_data->{$v}->{$name}->{version};
				$public_data->{$v}->{$name}->{config}->{version} = $public_data->{$v}->{$name}->{version};
				$public_data->{$v}->{$name}->{config}->{directory} =  $public_data->{$v}->{$name}->{config}->{name}."/".
																$public_data->{$v}->{$name}->{config}->{version}."/".
																$public_data->{$v}->{$name}->{config}->{dir}."/";
				#delete $public_data->{$v}->{$name}->{version};
			}
		}

		return $public_data;
		 
		
		},
);	

has gencode => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift; 
		my $dir = $self->config_directory."/public_data/";
		my $filename = $dir."gencode.cfg";
		confess($filename) unless -e $filename;
		read_config $filename => my %config1;
		my $filename2 = $dir."public-data.version.cfg";
		read_config $filename2 => my %config2;
		confess($filename2) unless -e $filename2;
		my $public_data;
		my $previous;

		foreach my $v ( keys %config1){
			$config1{$v}->{directory} = $config1{$v}->{name}."/".$config1{$v}->{version}."/".$config1{$v}->{dir};
		}
		return \%config1;
		},
);	

sub deja_vu_public_dir {
	my ($self,$version)= @_;
	confess() unless $version;
	return $self->{dj_pub_dir}->{$version} if exists $self->{dj_pub_dir}->{$version};
	 $self->{dj_pub_dir}->{$version} =  $self->config->{deja_vu}->{path_rocks}."/".$version . "/".$self->config->{deja_vu}->{variations} if (exists $self->config->{deja_vu}->{path});
	return $self->{dj_pub_dir}->{$version}  if -e $self->{dj_pub_dir}->{$version};
	confess("\n\nERROR: path dejavu not found in genbo.cfg  -> $version Die\n\n".$self->{dj_pub_dir}->{$version} );
}

sub getNeedlemanWunsch {
	my ($self,$opt,$type) = @_;
	$type ++;
	return $self->{nw}->{$type} if exists $self->{nw}->{$type};
	$self->closeNeedlemanWunsch() if scalar (keys %{$self->{nw}}) > 2;
	$self->{nw}->{$type} =  new NeedlemanWunsch($self,@$opt);
	return $self->{nw}->{$type} ;
}

sub closeNeedlemanWunsch {
	my ($self,$opt,$type) = @_;
	foreach my $type (keys %{$self->{nw}}){
		$self->{nw}->{$type}->destructor;
		$self->{nw}->{$type} = undef;
		delete $self->{nw}->{$type};
	}
	
	
	
}

sub get_genbo_annotation_name {
	my ($self, $annot) = @_;
	confess("\n\nERROR: $a annotation not found in config->{genbo_annotations_names} file... Die...\n\n") unless (exists $self->config->{genbo_annotations_names}->{lc($annot)});
	return $self->config->{genbo_annotations_names}->{lc($annot)};
}

sub get_annotations_files {
	my ($self,$version,$type,$file_type) = @_;
	die($file_type." ".$type." ".Dumper $self->public_data_config) unless exists $self->public_data_config->{$file_type}->{$type};
	my $name = $self->config->{'public_data'}->{$version}."/".$self->public_data_config->{general}->{root_dir}."/".$type."/".$self->public_data_config->{version}->{$type}."/".$self->public_data_config->{$file_type}->{$type};
	
	die(Dumper( $self->public_data_config)."\n$type $file_type") unless exists $self->public_data_config->{$file_type}->{$type};
	die($name) unless -e $name;
	return $name;	
	
}

sub get_public_data_version {
	my($self,$database,$version) = @_;
	if ($version =~ /\./) {
		my @lTmp = split('\.', $version);
		$version = $lTmp[1];
	}
	return $self->public_data->{$version}->{$database}->{version} if (exists $self->public_data->{$version}->{$database});
	confess("\n\nERROR: dir:$database no version found for release $version\n\n");
} 

sub get_alamut_api_key_from_user_name {
	my ($self, $user_name) = @_;
	return $self->getQuery->getAlamutApiKeyFromUserName($user_name);
}

has query => (
	is		=> 'ro',
	lazy	=> 1,
	reader => 'getQuery',
	default	=> sub { 
		my $self = shift;
		my $config = $self->config;
#		warn "\t\t\t ***************   Query ";
		my $query;
		if (exists $config->{projectdb}->{polyprojectNGS}){ 
			$main::my_project_db = "PolyprojectNGS";
		
			 $query = QueryMoosePolyProjectNgs -> new ( all_config => $config, dbh=>$self->dbh,buffer=>$self );
		}
		else {
			$main::my_project_db = "Polyproject";
			 $query = QueryMoosePolyProject -> new ( all_config => $config, dbh=>$self->dbh );
		}
		return $query;
	},
);



has queryPanel => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub { 
		my $self = shift;
		my $config = $self->config;
		my $query;
		$query = QueryPolyPanel -> new ( all_config => $config, dbh=>$self->dbh );
		return $query;
	},
);
has queryPhenotype => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub { 
		my $self = shift;
		my $config = $self->config;
		my $query;
		$query = QueryPolyPhenotype -> new ( all_config => $config, dbh=>$self->dbh );
		return $query;
	},
);
has queryHgmd => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub { 
		my $self = shift;
		my $config = $self->config;
		my $query;
		$query = QueryHgmd -> new ( all_config => $config, dbh=>$self->dbh );
		return $query;
	},
);

has queryClinvarPathogenic => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub { 
		my $self = shift;
		my $config = $self->config;
		my $query;
		$query = QueryClinvarPathogenic -> new ( all_config => $config, dbh=>$self->dbh );
		return $query;
	},
);

has queryFilter => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub { 
		my $self = shift;
		my $config = $self->config;
		my $query;
		my $project_db = "Polyproject";
		if (exists $config->{projectdb}->{polyprojectNGS}){ 
			$project_db = "PolyprojectNGS";
		}
		
		return QueryFilter->new(dbh=>$self->dbh,database=>"$project_db");
	},
);



has dbh =>(
	is		=> 'ro',
	lazy =>1,
	default => sub {
	my $self = shift;
	confess() if (exists $self->{debug} );
	my $dbh = connect::getdbh($self->config->{polyprod}) ;
	return $dbh;
	}

);


has getAllGenesNamesInAllBundle =>(
	is		=> 'rw',
	lazy =>1,
	default => sub {
		my $self = shift;
		my $query = $self->getQuery();
		my $res = $query->getAllGenesNamesInAllBundle();
		return $res;
	}

);

has getOmimTranscriptsNames =>(
	is		=> 'rw',
	lazy =>1,
	default => sub {
		my $self = shift;
		my $h;
		my $query = $self->getQuery();
		my $res = $query->getOmimGenesTranscriptsNames();
		foreach my $capture (keys %$res) {
			foreach my $panel (keys %{$res->{$capture}}) {
				foreach my $tr_name (@{$res->{$capture}->{$panel}->{transcripts}}) {
					$h->{$tr_name} = undef;
				}
			}
		}
		return $h;
	}

);


sub ucsc2ensembl {
	my ($shift,$chr) = @_;
	if ($chr =~ /chr/){
		$chr =~ s/chr//;
		$chr = 'MT' if $chr eq 'M';
	}
	return $chr;
}

sub newProjectCache {
	my $self = shift;
	my $args = _checkArguments(@_);
	my $name = 'undef';
	my $release = 'undef';
	my $typeFilters = 'individual';
	my $test = 0;
	if (exists $args->{-test}){ $test = 1; }
	if (exists $args->{-name}){ $name = $args->{-name}; }
	if (exists $args->{-release}){ $release = $args->{-release}; }
	if (exists $args->{-typeFilters}){ $typeFilters = $args->{-typeFilters}; }

	my $project = GenBoProjectCache -> new ( 	name => $name,
												buffer => $self,
												release => $release,
												cache => 1,
												typeFilters => $typeFilters );
	$self->genome_version($project->genome_version);	
	$self->annotation_genome_version($project->annotation_genome_version);
	$self->annotation_version($project->annotation_version);
	 $self->public_data_version($project->public_database_version);			
#	if ($project->annotation_version) { $self->lmdb_public_dir($project->annotation_public_path); }	
	 $self->public_data_version($project->public_database_version);							
	return $project;
}
has genome_version => (
	is		=> 'rw',
);
has annotation_genome_version => (
	is		=> 'rw',
);
has annotation_version => (
	is		=> 'rw',
);
has public_data_version => (
	is		=> 'rw',
);
has public_data_root => (
	is		=> 'ro',
	lazy => 1,
	default => sub {
		my $self = shift;
		my $d = $self->config->{'public_data_annotation'}->{root} ."/".$self->config->{'public_data_annotation'}->{repository};
		confess($d) unless -e $d;
		return $self->config->{'public_data_annotation'}->{root} ."/".$self->config->{'public_data_annotation'}->{repository};
	},
);
sub newProject {
	my $self = shift;
	my $args = _checkArguments(@_);
	my $name = 'undef';
	my $release = 'undef';
	my $test = 0;
	if (exists $args->{-test}){ $test = 1; }
	if (exists $args->{-name}){ $name = $args->{-name}; }
	if (exists $args->{-release}){ $release = $args->{-release}; }
	my $project = GenBoProject -> new ( name	=> $name,
										release	=> $release,
										test	=> $test,
										buffer	=> $self );
	$self->genome_version($project->genome_version);	
	$self->annotation_genome_version($project->annotation_genome_version);
	if ($project->gencode_version =~ /M/) {
		$self->annotation_version($project->annotation_version);
	 	$self->public_data_version($project->public_database_version);	
	}
	elsif ($project->gencode_version > -1){
		$self->annotation_version($project->annotation_version);
	 	$self->public_data_version($project->public_database_version);	
	}
	if (exists $args->{-version} && $args->{-version}){ 
 			$project->genome_version( $args->{-version});
 			$project->version( $args->{-version});
	}	
	#if ($project->annotation_version) { $self->lmdb_public_dir($project->annotation_public_path); }	
	return $project;
}

sub _checkArguments {
	my $index;
    for ($index=0; $index < @_; $index += 2) {
        my $key = $_[$index];
        unless ($key =~ /^\-/o) {
            confess ("Please, could you be so kind as to check your arguments for method \'construct\'? I have the impression you wrote \'$key\' instead of \'$key\' -- didn\'t you?\n");
            die();
            return undef;
        }
    }
    my %args = @_;
    # use lowercase in %args's keys -- usefull in case nameSpace was written instead of namespace!
    foreach my $key (keys %args) {
        $args{lc($key)} = $args{$key};
    }
   return \%args;
}

sub config_database {
	my ($self,$db) = @_;
	$db = lc($db);
	$self->config->{server}->{type_db} = 1;
	delete $self->config->{server}->{name};
	return;
}

sub getDataDirectory {
	my ($self, $type) = @_;
	confess unless exists  $self->config->{'project_data'}->{$type};
	return $self->config->{'project_data'}->{$type};
} 

sub saveStore {
	my ($self,$data,$file) = @_;
	my $serialized = freeze $data;	
	my $status = gzip \$serialized => $file or die $GzipError;
}

sub getStore {
	my ($self,$file) = @_;
	my $buffer;
	if ($file =~ /\.gz/){
	
	my $status = gunzip $file => \$buffer or die $GunzipError;
	
	my $code = thaw $buffer;
	return $code;
	}
	confess();
}
has get_go_db => (
	is 		=> 'ro',
	lazy =>1,
	default	=> sub{
		my $self =shift;
		eval {
	require "GO/AppHandle.pm";
	};
	if ($@) {
#		die();
		return undef;
	} 
#	confess();
eval {
	my $dbname =$self->config->{gene_ontology}->{dbname};
	my $mysqlhost =$self->config->{gene_ontology}->{ip};
	$mysqlhost=$ENV{POLY_DB_IP} if exists $ENV{POLY_DB_IP};
	my $login =$self->config->{gene_ontology}->{user};
	my $pwd =$self->config->{gene_ontology}->{pwd};
	$pwd =$self->config->{gene_ontology}->{pw} unless $pwd;
	my $aph;
	if ($pwd){
	 $aph = GO::AppHandle->connect(-dbname=>$dbname, -dbhost=>$mysqlhost, -dbuser=>$login,-dbauth=>$pwd);
	}
	else {
		$aph = GO::AppHandle->connect(-dbname=>$dbname, -dbhost=>$mysqlhost, -dbuser=>$login);
	}
	$aph->filters({evcodes=>["!IEA"], taxid => ["9606"]});
	return $aph;
};
	if ($@) {
#		die();
		return undef;
	} 
	} ,
);

sub listProjects {
		my $self = shift;
		my $query = $self->getQuery();
		my $list = $query->listProjects;
		return $list;
}

sub listProjectsByAnalyse {
		my ($self,$analyse )= @_;
		my $query = $self->getQuery();
		my $list = $query->listAllProjectsNameByAnalyse($analyse);
		return $list;
}

sub listProjectsForDejaVu{
		my $self = shift;
		my $query = $self->getQuery();
		my $list = $query->listProjectsForDejaVu;
		return $list;
}
sub listProjectsExomeForDejaVu{
		my $self = shift;
		my $query = $self->getQuery();
		my $list = $query->listProjectsExomeForDejaVu;
		return $list;
}

sub listAllPatientsId {
		my $self = shift;
		my $query = $self->getQuery();
		my $list = $query->listAllPatientsId;
		return $list;
}

sub hashAllPatientsInfos {
	my ($self, $by) = @_;
	my $query = $self->getQuery();
	my $list = $query->hashAllPatientsInfos($by);
	return $list;
	
}

sub getFindPatientDescription {
	my ($self, $patient_name) = @_;
	my @lHash;
	return $self->hashAllPatientsInfos('name')->{$patient_name};
}

sub getProjectNameFromId {
	my ($self,$poject_id) = @_;
	my $query = $self->getQuery();
	return $query->getProjectNameFromId($poject_id);
}

sub software {
		my ($self,$name) = @_;
		return $self->getSoftware($name);
}

sub software_version {
		my ($self,$name,$nodie) = @_;
		my $prog =  $self->getSoftware($name,1);
		return {"name" => $prog,"not_avalaible"=>1} unless $prog;
		my $rp = abs_path($prog);
		my @p = split("/",$rp);
		pop @p;
		my $version_json = join("/",@p)."/version.json";
		if ($nodie){
		return {} unless  -e $version_json;
		}
		else {
			confess("$name no version $version_json") unless -e $version_json;
		}
		open(JSON, $version_json);
		#warn $self->lmdb_public_dir."$database/"."description.json";
 		my $desc = decode_json <JSON>;
 		close (JSON);
		return $desc;
		#chomp($rp);
}

sub index {
		my ($self,$name) = @_;
		confess($name) unless exists $self->config->{index}->{$name};
		return $self->config->{index}->{$name};
}
sub getSoftware {
	my ($self,$name,$nodie) = @_;
	if($nodie){
		return unless exists $self->config->{software}->{$name};
	}
	confess($name) unless exists $self->config->{software}->{$name};
	return $self->config->{software}->{$name};
}

sub samtools {
	my ($self) = @_;
	return $self->getSoftware("samtools");
}

sub getListAllProjectName {
	my $self = shift;
	my @lProjectName;
	my $dir = $self->config()->{project_data}->{root}.'/'.$self->config()->{project_data}->{'ngs'}.'/';
	opendir(PROJECTS_PATH, "$dir");
	### TODO: change method (use DB and not readdir... can have some dir (not deleted) from deleted project)
	my @lProjectsPath = sort(readdir(PROJECTS_PATH));
	closedir(PROJECTS_PATH);
	foreach my $projectName (@lProjectsPath) { if ($projectName =~ /NGS20/) { push(@lProjectName, $projectName); } }
	return \@lProjectName;
}

has cacheFile => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		### TODO: to define path sqlite result...
		return '/home/mbras/sqlite/test.sqlite';
	}
);

has getHashTransIdWithCaptureDiag => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $query = $self->getQuery();
		return $query->hashTransIdWithCaptureDiag();
	}
);


has build => (
	is => 'rw',
	lazy => 1,
	default => 'HG19',
);

has ensemblVersion => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
	
		return $self->config->{ensembl}->{$self->build()};
	},
);





sub cache {
	my ($self, $args) = @_;
	my @lProjectName = split(',', $args->{projectName});
	my @lBatchFiles;
	my $pm = new Parallel::ForkManager(int($args->{maxProc}));
	foreach my $projectName (@lProjectName) {
		$pm->start and next;
		my $project = $self->newProject( -name=>$projectName );
		$project->createTmpCache();
		$project = undef;
		$pm->finish;
	}
	$pm->wait_all_children;
	foreach my $projectName (@lProjectName) {
		my $project = $self->newProject( -name=>$projectName );
		#$project->getQuerySqlLite()->insertSqlFile();
		#if ($projectName eq $lProjectName[-1]) {
		#	$project->getQuerySqlLite()->createIndex();
		#	$project->getQuerySqlLite()->vacuumFreeSpace();
		#}
		$project = undef;
	}
}

sub open_kyoto_db {
	my ($self,$file,$type) = @_;
	confess() unless $file;
	confess() unless $type;
	confess();
# my $db1 = new KyotoCabinet::DB;
#	if ($type eq 'r'){
#	 		 if (!$db1->open($file, $db1->ONOLOCK | $db1->OREADER )){
#						printf STDERR ("open error: %s\n", $db1->error);
#						confess();
#			}
#	}
#	elsif ($type eq 'c'){
#			unlink $file if -e $file;
#				if (!$db1->open($file, $db1->ONOLOCK |  $db1->OCREATE |$db1->OWRITER)){
#						printf STDERR ("open error: %s   => %s\n", $file,$db1->error);
#						confess();
#				}
#					system("chmod a+w $file");
#					$db1->clear();	
#					$db1->set("coucou","coucou");
#	}
#	elsif ($type eq 'w'){
#			my $new;
#			$new = 1 unless -e $file;
#				if (!$db1->open($file, $db1->ONOLOCK  |  $db1->OCREATE | $db1->OWRITER| $db1->OREADER )){
#						printf STDERR ("open error: %s   => %s\n", $file,$db1->error);
#						confess();
#				}
#				if ($new){
#						system("chmod a+w $file");
#				}
#					#$db1->clear();	
#	}
#	else {
#		die();
#	}
#		
#	return $db1;
}


has lmdb_public_dir => (
	is		=> 'rw',
);

has lmdb_public_dir1 => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
	my $self = shift;
	confess();
	}
);


has value_mask_database => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
	my $self = shift;
	return {
	'dbsnp'			=> 1,
 	'gnomad-genome'	=> 2,
 	'gnomad-exome'		=> 4,
 	'evs'			=> 8,	
 	'exac'			=> 16,
 	'clinvar'			=> 32,
 	'clinvar_precious'			=> 64,
};		
	
	}
);

sub description_gnomad {
	my ($self,$db) = @_;
	return $self->{description}->{$db} if exists $self->{description}->{$db};
	$self->{description}->{$db}= $self->description_public_lmdb_database($db);
	my @pops =  @{$self->{description}->{$db}->{array}->{populations}};
	$self->{description}->{$db}->{array}->{populations} = \@pops;
	return 	$self->{description}->{$db};
}

sub close_gnomad {
	my ($self,$chr) = @_;
	foreach my $database (keys %{$self->{lmdb}->{$chr}}){
		
		foreach my $type (keys %{$self->{lmdb}->{$chr}->{$database}}){
			
			$self->{lmdb}->{$chr}->{$database}->{$type}->close();
			
			delete $self->{lmdb}->{$chr}->{$database}->{$type};
		}
		delete $self->{lmdb}->{$chr}->{$database};
	}
	
}

sub get_gnomad {
	my ($self,$chr,$type,$start,$allele) = @_;
	my @databases = ("gnomad-genome","gnomad-exome");
	my $res ={};
	my $rsname ="";
	my $public;
	my $h;
			foreach my $db  (@databases) {
				my $desc = $self->description_gnomad($db);
				#warn Dumper $desc;
				my $pops = $desc->{array}->{populations};
				my $infos = $desc->{array}->{infos};
				my $hash = $self->get_lmdb_database($db,$chr,$type)->get($start);
				
				
				next unless $hash;
				warn $start if $hash =~ /-/;
				next unless exists $hash->{$allele};
				$public ++;
				 $h = $hash->{$allele};
					
				foreach my $info (@{$infos}){
					my @array = unpack($h->{$info}->[0],$h->{$info}->[1]);
				
					my $all;
					my $i =0;
					foreach my $p (@$pops){
						if (lc($p) eq "all"){
							$i++;
							next;
						}
						$res->{$p}->{$info} += $array[$i];
						$res->{all}->{$info} += $array[$i];
						$i++;
					}
				}
					$rsname = $h->{rs} if exists $h->{rs};
					$rsname = $h->{rsname} if exists $h->{rsname};
			}
			return {} unless $public;
			my $desc = $self->description_gnomad("gnomad-exome");
			my $pops = $desc->{array}->{populations};
			my $fmin;
			$fmin = 1000;
			my $minpop ="-";
			
			my $fmax= -2;
			my $maxpop="-";
			
			foreach my $p (@$pops){
				next unless exists  $res->{$p}->{AN};
				$res->{$p}->{F}  = "-1";
				warn Dumper $res unless exists  $res->{$p}->{AN};
				warn Dumper $h unless exists  $res->{$p}->{AN};
				die($p)  unless exists  $res->{$p}->{AN};
				if ($res->{$p}->{AN}>0){
					$res->{$p}->{F} = $res->{$p}->{AC}/$res->{$p}->{AN} if $res->{$p}->{AN}>0;
					if ($res->{$p}->{F} < $fmin){
						$fmin = $res->{$p}->{F};
						$minpop = $p;
					}
					if ($res->{$p}->{F} > $fmax){
						$fmax = $res->{$p}->{F};
						$maxpop = $p;
					}
				}
				
			}
	my $gnomad = {};
	$gnomad->{public} =1 if $public ==1;
	
	$gnomad->{populations} = $res;
	$gnomad->{rsname} = 	$rsname;
	$gnomad->{minpop} = 	$minpop;
	$gnomad->{maxpop} = 	$maxpop;

	return 	$gnomad;

	}
	
sub get_lmdb_database_directory{
	my ($self,$database)= @_;
	my $version = $self->public_data_version;
	return $self->public_data_root."/".$self->annotation_genome_version."/".$self->public_data->{$version}->{$database}->{config}->{directory};
}


sub get_version_database{
	my ($self,$database)= @_;
	my $version = $self->public_data_version;
	return $self->public_data->{$version}->{$database}->{config}->{version};
}


sub get_index_database_directory{
	my ($self,$database)= @_;
	my $version = $self->public_data_version;
	return $self->public_data_root."/".$self->annotation_genome_version."/$database/".$self->public_data->{$version}->{$database}->{config}->{version}."/".$self->public_data->{$version}->{$database}->{config}->{dir};
}


sub description_public_lmdb_database {
	my ($self,$database)= @_;

	 return $self->{config}->{$database} if exists $self->{config}->{$database};
	my $version = $self->public_data_version;
	my $dir = $self->public_data_root."/".$self->annotation_genome_version."/".$self->public_data->{$version}->{$database}->{config}->{directory};
	
	my $f ="$dir/description.json";
	 $f  =  "$dir/version.json"  unless -e $f;
	 confess($f) unless -e $f;
	  open(JSON, $f);
#warn $self->lmdb_public_dir."$database/"."description.json";
 	my $desc = decode_json <JSON>;
	#my $pop = $desc->{array}->{populations};
	$self->{config}->{$database} = $desc;	
	return $self->{config}->{$database};
}

sub mask_database {
	my ($self,$db)= @_;
	confess() unless exists $self->value_mask_database->{$db};
	return  $self->value_mask_database->{$db};
}

sub get_lmdb_database {
	my ($self,$database,$chr,$type) = @_;
	confess() unless $type;
	return $self->{lmdb}->{$chr}->{$database}->{$type} if exists  $self->{lmdb}->{$chr}->{$database}->{$type};
	my $dir = $self->get_lmdb_database_directory($database).'/'.$type;
	confess($dir) unless -e $dir;
	unless (-e $dir."/".$chr ){
		if ($chr eq "Y"){
				 my $lmdb = GenBoNoSqlLmdb->new(dir=>$dir,mode=>"w",name=>$chr,is_compress=>1,is_integer=>1,vmtouch=>$self->vmtouch);
				 $lmdb->put(0,"toto");
				 $lmdb = undef;
		}
		if ($chr eq "MT"){
		 my $lmdb = GenBoNoSqlLmdb->new(dir=>$dir,mode=>"c",name=>$chr,is_compress=>1,is_integer=>1,vmtouch=>$self->vmtouch);
		 $lmdb->put(0,"toto");
		 $lmdb = undef;
		} 
		else {
		warn $database;
		warn $dir."/".$chr;
		confess();
		die();
		}
	#	die();
	}
	$self->{lmdb}->{$chr}->{$database}->{$type} = GenBoNoSqlLmdb->new(dir=>$dir,mode=>"r",name=>$chr,is_compress=>1,vmtouch=>$self->vmtouch);
	return $self->{lmdb}->{$chr}->{$database}->{$type};
}

sub get_lmdb_clinical_local_db {
	my ($self,$name) = @_;	
	my $database = "clinical_local_db";
	return $self->{lmdb}->{$name}->{$database}->{all} if exists  $self->{lmdb}->{$name}->{$database}->{all};
	my $dir = $self->get_lmdb_database_directory($database);
	confess($self->lmdb_public_dir."/$database/$name") unless -e $self->lmdb_public_dir."/$database/$name";
	 $self->{lmdb}->{$name}->{$database}->{all}  = GenBoNoSqlLmdb->new(dir=>$dir,mode=>"r",name=>$name,is_compress=>1,vmtouch=>$self->vmtouch);
	return  $self->{lmdb}->{$name}->{$database}->{all} ;
	
}



sub close_lmdb {
		my ($self) = @_;
	foreach my $chr (keys %{$self->{lmdb}}){
			foreach my $db (keys %{$self->{lmdb}->{$chr}}){
					foreach my $type (keys %{$self->{lmdb}->{$chr}->{$db}}){
							$self->{lmdb}->{$chr}->{$db}->{$type}->close() if $self->{lmdb}->{$chr}->{$db}->{$type};
					}
			}
		
	}
	
	foreach my $a (values %{$self->{lmdb_hash}}) {
		foreach my $c (values %{$a}) {
		#	warn $c;
			$c->close();
			
		}
	}
	delete $self->{lmdb_hash};
	delete $self->{lmdb};
	
}

  
sub DESTROY {
	my ($self) = @_;
	#delete $self->{lmdb};
	my $t =time;
	#swarn "detroy buffer ";
	foreach my $chr (keys %{$self->{lmdb}}){
			foreach my $db (keys %{$self->{lmdb}->{$chr}}){
					foreach my $type (keys %{$self->{lmdb}->{$chr}->{$db}}){
							$self->{lmdb}->{$chr}->{$db}->{$type}->close() if $self->{lmdb}->{$chr}->{$db}->{$type};
					}
			}
		
	}
foreach my $database (keys %{$self->{lmdb_score}}){
		foreach my $chr (keys %{$self->{lmdb_score}->{$database}}){
			$self->{lmdb_score}->{$database}->{$chr}->close  if $self->{lmdb_score}->{$database}->{$chr};
			delete $self->{lmdb_score}->{$database}->{$chr};
		}
}
delete 	$self->{lmdb_score};
delete $self->{lmdb};
}

sub intersection {
   my ($self,$hin,$left) = @_;
	my %n = map { $_ => undef } grep exists $hin->{$_}, @$left;
  return (keys %n);
}
sub coverage_samtools {
	my($self,$bam,$chr,$start,$end) = @_;
	my $samtools = $self->samtools;
		my $len = abs($start - $end) +1;
		my @v = ((0) x $len);
		my $region = $chr.":".$start."-".$end;
		my @t = `$samtools depth -Q 1 $bam -r $region | cut -f 2-3`;
		chomp(@t);
		foreach my $tt (@t ){
			my ($p,$c) = split(" ",$tt);
			my $pv = $p - $start ;
			next if $pv <0;
			$v[$pv] = $c;
			last if $pv > $end;
		
		}
		my $gc  = GenBoCoverage->new(start=>$start,end=>$end,array=>\@v);
		return ($gc);
}

sub coverage_tabix {
	my($self,$coverage_file,$chr,$start,$end) = @_;
	confess() unless -e $coverage_file;
	my $tabix = Bio::DB::HTS::Tabix->new( filename => $coverage_file );
	my $region = $chr.":".$start."-".$end;
	my $len = abs($start - $end) +1;
	my @v = ((0) x $len);
	 my $res = $tabix->query($region) if $start;
		 my @data;
		
		 while(my $line = $res->next){
				my($a,$p,$c) = split(" ",$line);
				confess() if $a ne $chr;
				my $pv = $p - $start ;
				 $v[$pv] = $c;
			}

		my $gc  = GenBoCoverage->new(start=>$start,end=>$end,array=>\@v);
		return ($gc);
}

sub getListCaptureDiagFromTransId {
	my ($self, $t_id) = @_;
	my @lTmp = split('_', $t_id);
	my $id = $lTmp[0];
	my $query = $self->getQuery();
	my $list = $query->listCaptureDiagnosticFromTranscriptId($id);
	return $list;
}


sub disconnect {
	my($self,$project) = @_;
	$self->dbh_deconnect($project);
}
sub dbh_deconnect {
	my($self,$project) = @_;
	$self->{dbh}->disconnect if exists $self->{dbh};
	$self->{query}->{dbh}->disconnect if exists $self->{query}->{dbh};
	$self->{queryPanel}->{dbh}->disconnect if (exists $self->{queryPanel}->{dbh});
	delete $self->{validations_query};
	delete  $self->{queryPanel};
	delete $self->{query};
	delete $self->{dbh};
	#delete $self->{query};
	delete $self->{queryHgmd};
	delete $self->{queryPhenotype};
	
	if ($project){
		 delete  $project->{noSqlCoverage};
		 
	}
	
}

sub validations_query {
	my $self = shift;
	 return $self->{validations_query} if exists  $self->{validations_query};
	 $self->{validations_query} = QueryValidationAcmg->new(
				dbh      => $self->dbh,
				database => "ACMG"
			);
	 return $self->{validations_query};# if exists   $self->{validations_query1};
}

sub dbh_reconnect {
	my $self = shift;
	$self->dbh;
	#delete $self->{dbh};
	#delete $self->{query}->{dbh};
	#delete $self->{queryPanel}->{dbh};
	#$self->{dbh} = connect::getdbh($self->config->{polyprod});
	#$self->{query}->{dbh} = connect::getdbh($self->config->{polyprod});
	#$self->{queryPanel}->{dbh} = connect::getdbh($self->config->{polyprod}) if (exists $self->{queryPanel});
}
 
sub intspanToBed{
	my ($self,$chr,$intspan) = @_;
	my $iter = $intspan->iterate_runs();
	my @tt;
    while (my ( $from, $to ) = $iter->()) {
    		push(@tt,$chr->fasta_name."\t".$from."\t".$to);
    	
    }
	#	my @tt = map{$_ = $chr->ucsc_name."\t".$_} split(";",$intspan->as_string({ sep => ";", range => "\t" }) );
		return @tt;
}

sub gzip_tabix {
	my ($self,$file,$type) = @_;
	my $bgzip = $self->software('bgzip');
	my $tabix = $self->software('tabix');

	if ($file =~ /gz$/){
		system(" $tabix -f -p $type $file");
	}
	else {
		system("$bgzip -f $file && $tabix -p $type $file.gz");
		$file .= ".gz";
	}
	
	confess($file) unless -e "$file.tbi";
	return "$file"
}


sub public_data_annotation_root {
		my ($self) = @_;
		return $self->config->{'public_data_annotation'}->{root};
}

sub Intspan_length{
	my ($self,$intspan) = @_;
	my $iter = $intspan->iterate_runs();
	my @tt;
	my $l =0;
    while (my ( $from, $to ) = $iter->()) {
    		$l += ($to-$from) +1;
    	
    }
	#	my @tt = map{$_ = $chr->ucsc_name."\t".$_} split(";",$intspan->as_string({ sep => ";", range => "\t" }) );
		return $l;
}

sub Intspan_start_end{
	my ($self,$intspan) = @_;
	my $iter = $intspan->iterate_runs();
	my @tt;
	my $l =0;
	my $start = -1;
	my $end =0;
	
    while (my ( $from, $to ) = $iter->()) {
    	$start = $from if $start == -1;
    	$end = $to;
    	
    }
	#	my @tt = map{$_ = $chr->ucsc_name."\t".$_} split(";",$intspan->as_string({ sep => ";", range => "\t" }) );
		return ($start,$end);
}

sub divide_by_chunks {
	my ($self,$start,$end,$chunksize) = @_;
	
	my $number = ($end-$start)-1;
	die if $number <0 ;
	
	if ($chunksize >= $number){
		my $r = [$start,$end];
		return [$r];
	}
	my $interval;
	
	#my $chunksize = int($number/$parts)+1;             # size of each interval
    my $chunkstart = 1;                        # start of interval
    my $chunkend = $chunkstart + $chunksize -1;  # end of that interval
    $start = $start -1;
    while ($chunkstart < $number){            # don't go beyond the range
       push(@$interval,[$chunkstart+$start,$chunkend+$start]);
        $chunkstart += $chunksize;           # take me to beginning of next interval
        $chunkend += $chunksize;             # also tell me where to end that
	 if ($chunkend >= $number)  {          # interval end is beyond the range
             push(@$interval,[$chunkstart+$start,$end]);
            last;                         # we are beyond the range now
        }
}
	return $interval;
}

sub get_annotation_terms{
	my ($self,$term) = @_;
	return $term unless exists $self->config->{ensembl_annotations}->{$term};
	confess($term) unless exists $self->config->{ensembl_annotations}->{$term};
	my @terms =  split(";", $self->config->{ensembl_annotations}->{$term}); 
	return lc($terms[0]);
}

sub getListFiles {
	my ($self, $dir) = @_;
	return $self->{dir}->{$dir} if exists $self->{dir}->{$dir};
	opendir (DIR, $dir) or die $!;
	$self->{dir}->{$dir} = [];
	while (my $file = readdir(DIR)) {

       push(@{$self->{dir}->{$dir}},$file);

    }
	return $self->{dir}->{$dir};
}


sub appendTreeFromVector {
        my ($self,$tree,$vector,$name,$dec) = @_;
        $dec =0 unless $dec;
         my @enum =  split(",",$vector->to_Enum());
                #warn Dumper @enum;
                        foreach my $en (@enum){
                                my ($start,$end) = split("-",$en);
                                $end = $start unless $end;
                                $tree->insert($name,$start+$dec,$end+$dec+1);
                        }
        return ;
}

sub bundle_infos_all_proj_users {
	my ($self, $user) = @_;
	confess() unless ($user);
	my $query = $self->getQuery();
	my $res = $query->getBundleGenesAllProjectsUsers($user);
	return $res;
}

sub getGenesNamesInBundle {
	my ($self, $bundle_id) = @_;
	confess() unless ($bundle_id);
	my $query = $self->getQuery();
	my $res = $query->getGenesNamesInBundle($bundle_id);
	return $res;
}

has hash_genes_omim_morbid =>(
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		my $self = shift;
		return $self->queryPanel->getGenesOmimMorbid();
	}
);

sub test_fisher {
	my ($self, $sample_var_1, $sample_var_2, $min_p) = @_;
	my ($t,$a,$c)  = split(":", $sample_var_1);
	my ($t1,$b,$d) = split(":", $sample_var_2);
	$a = 0 unless ($a);
	$b = 0 unless ($b);
	$c = 0 unless ($c);
	$d = 0 unless ($d);
	my $p = fisher::fishers_exact( $a, $b, $c, $d, 1);
	return ($p, 0) if $p > $min_p;
	return ($p, 1);
}


		my $validation =  {
		 pathogenic => 5,
  		"likely pathogenic" => 4,
  		"Uncertain significance" =>3,
  		"Likely benign" => 2,
  		"benign" => 1,
  		"False Positive" => -1,
  		"ToDo" => -3
	};


has value_validation => (
	is 		=> 'rw',
	default	=> sub {
		my $self = shift;
		my $z ;
		foreach my $k (keys %{$validation} ){
			my $v = $validation->{$k};
			$z->{$v} = $k; 
		}
		return $z;
},
	
);

sub getValidationTerm {
	my ($self,$value) = @_;

	return "unknown"  unless exists $self->value_validation->{$value};
	return $self->value_validation->{$value};
	
}


has gene_annotation_version_current => (
	is		=> 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		return $self->getQuery->getCurrentGenesReleasesAnntotations();
	},
);

has list_all_annotation_version => (
	is		=> 'rw',
	lazy => 1,
	default => sub {
		my $self = shift;
		my $query = $self->getQuery();
		my @lReleases = $query->getListAllReleasesAnntotations();
		return \@lReleases;
	},
);

sub get_polybtf_default_release {
	my ($self) = shift;
	return $self->config->{polybtf_default_releases}->{default};
}

sub get_polybtf_path {
	my ($self, $release) = @_;
	unless ($release) {
		$release = $self->get_polybtf_default_release();
	}
	my $dir = $self->config->{project_data}->{btf}.'/'.$release.'/';
	return $dir;
}

sub get_polybtf_project_resume {
	my ($self, $project_name, $release) = @_;
	my $json_file = $self->get_polybtf_path($release).'/'.$project_name.'/'.$project_name.'_new_public_db.resume.json';
	return unless (-e $json_file);
	open (FILE, "$json_file");
	my $json = <FILE>;
	close (FILE);
	return decode_json $json;
}

sub areSameCNV {
	my ( $self, $start1, $end1, $start2, $end2 ) = @_;

	my $len1 = abs( $start1 - $end1 );
	my $len2 = abs( $start2 - $end2 );

	# les CNV doivent se chevaucher à 60% dans les deux sens
	my $overlap = min( $end1, $end2 ) - max( $start1, $start2 );

	return 0 if $overlap < 0;
	return 0 if ( $overlap < ( 0.6 * $len1 ) );
	return 0 if ( $overlap < ( 0.6 * $len2 ) );

	return 1;
}

sub color_model {
	my ($self,$model) = @_;
	my $color = "#555";
	return $color unless $model;
	if ($model eq 'mother'){
		$color = "#F7E4E4";
	}
	elsif ($model eq 'mother_c'){
		$color = "#779ECB";
	}
	elsif ($model eq 'father'){
		$color = "#A9BCD1";
	}
	elsif ($model eq 'father_c'){
		$color = "#779ECB";
	}
	elsif ($model =~ /denovo/){
		$color = "#e74c3c";
	}
	elsif ($model =~ /rece/){
		$color = "#EE82EE";
	}
	elsif ($model =~ /mosa/){
		$color = "#F9885C";
		$color = "#FDE6B0";
	}
	elsif ($model =~ /uni/){
		$color = "#F9885C";
		$color = "#45B8AC";
	}
	return $color;
}

sub get_lmdb_cache_btf_view {
	my ( $self, $release, $user_name, $mode) = @_;
	$mode = "r" unless $mode;
	my $dir_out = $self->get_polybtf_path($release).'/users/';
	my $file = $user_name.'.polybtf.view.cache';
	unless (-e $dir_out.'/'.$file){
		$mode = "c";
	}
	my $no2  = GenBoNoSqlLmdbCache->new(
		dir     => $dir_out,
		mode    => $mode,
		name    => $file,
		is_compress => 1,
		vmtouch => $self->vmtouch
	);
	if ( $mode eq "c"){
		$no2->put("cdate",time);
		system("chmod a+w ".$no2->filename);
	}
	return $no2;
}

sub get_lmdb_cache_btf_resume {
	my ( $self, $release, $user_name, $mode) = @_;
	$mode = "r" unless $mode;
	my $dir_out = $self->get_polybtf_path($release).'/users/';
	my $file = $user_name.'.polybtf.resume.cache';
	unless (-e $dir_out.'/'.$file){
		$mode = "c";
	}
	my $no2  = GenBoNoSqlLmdbCache->new(
		dir     => $dir_out,
		mode    => $mode,
		name    => $file,
		is_compress => 1,
		vmtouch => $self->vmtouch
	);
	if ( $mode eq "c"){
		$no2->put("cdate",time);
		system("chmod a+w ".$no2->filename);
	}
	return $no2;
}

has get_hash_projects_ids_genes_databases_version => (
	is		=> 'rw',
	lazy => 1,
	default	=> sub {
		my $self = shift;
		my $h_proj = $self->getQuery->get_hash_projects_ids_genes_databases_version();
		return $h_proj;
	},
);

sub get_random_project_name_with_this_annotations_and_genecode {
	my ($self, $genecode, $annotdb) = @_;
	$genecode = $self->getQuery->getMaxGencodeVersion() unless ($genecode);
	$annotdb = $self->getQuery->getMaxPublicDatabaseVersion() unless ($annotdb);
	my $h_proj = $self->get_hash_projects_ids_genes_databases_version();
	foreach my $proj_id (sort {$a <=> $b} keys %$h_proj) {
		my $this_genecode = $h_proj->{$proj_id}->{rel_gene_id};
		my $this_db = $h_proj->{$proj_id}->{version_id};
		if ($genecode eq $this_genecode and $annotdb eq $this_db) {
			return $h_proj->{$proj_id}->{name};
		}
	}
	confess("\n\nNo project fount with GENECODE $genecode and ANNOT DB $annotdb... DIE...\n\n");
}

sub log2 {
	 my ($buffer,$n) = @_;
	 return -2 if $n <0.01;
	 my $v = log($n)/log(2);
	 $v =-2 if $v < -2; 
    return $v;
}
########
# SEREAL 
#######
#has sereal_encoder => (
#	is      => 'rw',
#	lazy    => 1,
#	default => sub {
#		my $self = shift;
#		#return Sereal::Encoder->new();
#		return Sereal::Encoder->new({compress=>Sereal::SRL_ZSTD,compress_threshold=>0});
#		return 0;
#	},
#);
#
#sub sereal_encode {
#	my ($self,$value)  =@_;
#	return sereal_encode_with_object($self->sereal_encoder, $value);
#}
#
#
#has sereal_decoder => (
#	is      => 'rw',
#	lazy    => 1,
#	default => sub {
#		my $self = shift;
#		return Sereal::Decoder->new({compress=>Sereal::SRL_ZSTD,compress_threshold=>0});
#		return 0;
#	},
#);
#sub sereal_decode {
#	my ($self,$value)  = @_;
#	return unless $value;
#	return sereal_decode_with_object($self->sereal_decoder, $value);
#}


sub get_url_polyrna {
	my ($self) = shift;
	return $self->get_base_url_polyrna().':'.$self->get_port_url_polyrna();
}

sub get_base_url_polyrna {
	my ($self) = shift;
	return $self->config->{polyrna}->{base_url_polyrna};
}

sub get_port_url_polyrna {
	my ($self) = shift;
	return $self->config->{polyrna}->{port_url_polyrna};
}

sub get_polyrna_file_server_to_docker {
	my ($self) = shift;
	return $self->config->{polyrna}->{server_to_docker};
}

sub get_polyrna_file_docker_to_server {
	my ($self) = shift;
	return $self->config->{polyrna}->{docker_to_server};
}

#####################
# COMPRESS HASH 
#####################
my $hash_dejavu = {
          'similar_patients_ho' => 0,
          'total_in_this_run_patients' => 1,
          'total_similar_projects' => 2,
          'total_exome_projects' => 3,
          'other_patients_ho' => 4,
          'in_this_run_patients' => 5,
          'other_patients' => 6,
          'exome_patients_ho' => 7,
          'exome_patients' => 8,
          'other_projects' => 9,
          'total_exome_patients' => 10,
          'similar_projects' => 11,
          'similar_patients' => 12,
          'total_similar_patients' => 13,
          'exome_projects' => 14,
          'other_projects_ho' => 15,
        };

sub index_dejavu {
	my ($self,$key) =@_; 
	confess() unless exists $hash_dejavu->{$key};
	return $hash_dejavu->{$key};
}

sub hash_to_array_dejavu{
	my ($self,$hash) =@_;
	my $array = [];
	foreach my $key (keys %$hash){
		my $index = $self->index_dejavu($key);
		$array->[$index] = $hash->{$key};
	}
	return $array;
}


sub get_demultiplex_run_infos {
	my ($self, $run_name) = @_;
	my $h_db = $self->getQuery->getInfosFromRunMachineId($run_name);
	$h_db = $self->getQuery->getInfosFromRunName($run_name) if scalar keys %$h_db == 0;
	return $h_db;
} 


1;