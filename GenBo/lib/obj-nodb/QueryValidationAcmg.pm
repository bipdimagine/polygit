package QueryValidationAcmg;
use strict;
use Data::Dumper;
use Carp;
use Moo;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use JSON::XS;


has 'dbh' => (
	is =>'ro',
	required => 1,
);
 
has 'db' =>(
	is =>'ro',	
	lazy=>1,
	default =>sub {
		my $self = shift;
		return "validation_OI" unless $self->database();
		my $dbname = "validation_".$self->database();
		my $query = qq{SHOW DATABASES LIKE '$dbname'};
		my $sth = $self->dbh->prepare($query);
		$sth->execute();
		my $s = $sth->fetchrow_hashref();
		return "validation_".$self->database() if $s;
		confess("no database validation for this capture validation_".$self->database());
		return "validation_OI" ;
	}
);

has 'database' =>(
	is =>'ro',	
	required => 1,
);

has 'exists_db'=>(
	is =>'ro',	
	lazy=>1,
	default =>sub {
		my $self = shift;
		my $dbname =  "validation_".$self->capture_name;
		my $query = qq{SHOW DATABASES LIKE '$dbname'};
		my $sth = $self->dbh->prepare($query);
		$sth->execute();
		my $s = $sth->fetchrow_hashref();
		return 1 if $s;
		return undef;
	}
);



sub getIdFromPositionSequence{
		my ($self,$chr,$start,$end,$seq,$type) = @_; 
	my $db = $self->db;
	my $query = qq{
		select variation_id as id  from $db.variations where chromosome = "$chr" and start=$start and end = $end and sequence = "$seq"  and type = "$type" and version = 'HG19';
	};
	my $sth = $self->dbh->prepare($query);
	$sth->execute();
	my $s = $sth->fetchrow_hashref();	
	return $s->{id};
}
		
sub getIdFromVariation{
		my ($self,$v) = @_; 
		return $self->getIdFromPositionSequence($v->getChromosome->name,$v->start,$v->end,$v->alternate_allele,$v->type_public_db);
}

sub insertVariation {
		my ($self,$v) = @_; 
	my $id = $v->id;
	 my $vcf_id = $v->vcf_id;
	 my $chr = $v->getChromosome->name();
	 my $start = $v->start;
	 my $end = $v->end;
	 my $seq = $v->alternate_allele();
	 my $type = $v->type_public_db();
	 my $name = $v->name;
	 my $db = $self->db;
	my $query = qq{
		insert into $db.variations (polyid,vcfid,version,chromosome,start,end,sequence,type,rs) values ("$id","$vcf_id",'HG19',"$chr",$start,$end,"$seq","$type","$name");
	};

	$self->dbh->do($query) || die($query);
	
	return $self->getIdFromVariation($v);
}

sub createVariation {
	my ($self,$v) = @_;
	my $id = $self->getIdFromVariation($v);

	unless ($id){
		
		$id = $self->insertVariation($v);
	}

	confess() unless $id;
	return $id;
}


sub createValidation {
	my ($self,$vid,$v,$patient,$gene,$user_name,$validation_type) = @_;
	my $project_name =$v->project->name();
	my $project_id = $v->project->id();
	my $phenotypes_id = join(';', map {$_->id()} @{$v->project->getPhenotypes()} );
	my $phenotypes_name = join(';', map {$_->name()} @{$v->project->getPhenotypes()} );
	my $patient_name = $patient->name;
	my $patient_id = $patient->id;
	my $gene_name = $gene->external_name;
	my $gene_id = $gene->id;
	my $i =  $v->annex()->{$patient->id};
	 my $infos = encode_json $i;
	 my $db = $self->db;
	my $query = qq{
		insert into $db.validations (variation_id,project_name,project_id,sample_name,sample_id,gene_name,gene_id, user_name,vcf_line,creation_date,modification_date,validation,phenotypes_id,phenotypes_name) values
		                                             ($vid, "$project_name",$project_id,"$patient_name",$patient_id,"$gene_name","$gene_id","$user_name",'$infos',NOW(),NOW(),$validation_type,"$phenotypes_id","$phenotypes_name") ON DUPLICATE KEY UPDATE modification_date = NOW()  ;
	};

	$self->dbh->do($query) || die($query);
	$self->addPatientsStatus($patient,$user_name,$validation_type);
	
}



sub getValidationPatientVariation {
	my ($self,$v,$patient) = @_;
	my $id = $self->getIdFromVariation($v);
		 my $db = $self->db;
	return unless $id;
	my $patient_name = $patient->name();
	my $project_name = $v->project->name();
	my $query = qq{select * from $db.validations where variation_id=$id and sample_name = "$patient_name" and project_name = "$project_name" order by validation_id}; 
	my $sth = $self->dbh->prepare($query) ;
	$sth->execute() || die();
	my $res = [];
	while (my $s =$sth->fetchrow_hashref()){
		my $gid = $s->{gene_id};
		push(@{$res->{$gid}},$s);
	}
	return $res;
}

sub getValidationPatient {
	my ($self,$patient) = @_;
	return $self->getValidationPatientName($patient->name(),$patient->project->name());
}

sub getValidationPatientName {
	my ($self,$patient_name,$project_name) = @_;
	my $db = $self->db;
	my $query = qq{select * from $db.validations as va ,$db.variations as v   where va.variation_id=v.variation_id and sample_name = "$patient_name" and project_name = "$project_name" order by  validation_id desc}; 
	my $sth = $self->dbh->prepare($query) ;
	
	$sth->execute() || die();
#	confess();
#	warn Dumper $sth->fetchrow_hashref();
return $self->hresults($sth);
}

sub getValidationPatientId {
	my ($self,$patient_id,$project_name) = @_;
	my $db = $self->db;
	my $query = qq{select * from $db.validations as va ,$db.variations as v   where va.variation_id=v.variation_id and sample_id = "$patient_id" and project_name = "$project_name" order by  validation_id desc}; 
	my $sth = $self->dbh->prepare($query) ;
	
	$sth->execute() || die();
	
	return $self->hresults($sth);
}

sub getValidationProjectName {
	my ($self,$project_name) = @_;
	my $db = $self->db;
	my $query = qq{select * from $db.validations as va   where  project_name = ? order by  modification_date asc}; 
	my $sth = $self->dbh->prepare($query) ;
	$sth->execute($project_name) || die();
	return $self->mhresults($sth);
	#return $self->hresults($sth);
}

sub  mhresults {
	my ($self,$sth) = @_;
	my $res =[];
	while (my $s =$sth->fetchrow_hashref()){
		
		push(@{$res},$s);
	}
	return $res;
}
sub  hresults {
	my ($self,$sth) = @_;
	my $res = {};
	while (my $s =$sth->fetchrow_hashref()){
		my $gid = $s->{gene_id};
		my $vid = $s->{polyid};
		$gid="+" unless $gid;
		push(@{$res->{$gid."!".$vid}},$s);
	}
	return $res;
}

sub  hByPatients{
	my ($self,$sth) = @_;
	my $res = {};
	while (my $s =$sth->fetchrow_hashref()){
		my $vid = $s->{sample_id};
	#	$gid="+" unless $gid;
		push(@{$res->{$vid}},$s);
	}
	return $res;
}

sub getValidationsProject {
	my ($self,$project) = @_;
	my $db = $self->db;
	my $project_id = $project->id();
	my $query = qq{select * from $db.validations as va ,$db.variations as v   where va.variation_id=v.variation_id and project_id = "$project_id" order by validation_id desc}; 
	my $sth = $self->dbh->prepare($query) ;
	$sth->execute() || die();
	return $self->hresults($sth);
}

sub getAllValidations {
	my ($self,$validation) = @_;
	my $db = $self->db;
	$validation = -100 unless $validation;
	my $query = qq{select *,UNIX_TIMESTAMP(va.modification_date) as unix_time from $db.validations as va ,$db.variations as v,$db.acmg     where va.variation_id=v.variation_id and idacmg=validation and validation >= $validation order by unix_time desc}; 
	my $sth = $self->dbh->prepare($query) ;
	$sth->execute() || die();
	return $self->hresults($sth);
}

sub getAllValidations_last_months {
	my ($self,$validation,$nb_month) = @_;
	my $db = $self->db;
	$validation = -100 unless $validation;
	my $query = qq{select * , UNIX_TIMESTAMP(va.modification_date) as unix_time 
		from validation_ACMG.validations as va ,validation_ACMG.variations as v,validation_ACMG.acmg 
			where va.variation_id=v.variation_id and idacmg=validation and validation >= $validation and va.modification_date > DATE_ADD(NOW(), INTERVAL -$nb_month MONTH)
				 order by unix_time desc;};
	my $sth = $self->dbh->prepare($query) ;
	$sth->execute() || die();
	return $self->hresults($sth);
}

sub getAllValidationsForVariation {
	my ($self,$vid) = @_;
	my $db = $self->db;
	my $query = qq{select *,UNIX_TIMESTAMP(va.modification_date) as unix_time from $db.validations as va ,$db.variations as v,$db.acmg     where v.polyid= "$vid" and va.variation_id=v.variation_id and idacmg=validation order by unix_time desc}; 
	my $sth = $self->dbh->prepare($query) ;
	$sth->execute() || die();
	return $self->hresults($sth);
}

sub getAllValidationsForPatient {
	my ($self,$patient,$user_name) = @_;
	my $db = $self->db;
	my $pid = $patient->id;
	my $query = qq{select *,term from $db.validations as va ,$db.variations as v,$db.acmg   where va.variation_id=v.variation_id  and va.sample_id=$pid and idacmg=validation order by validation_id desc}; 
	if ($user_name){
		$query = qq{select *,term from $db.validations as va ,$db.variations as v,$db.acmg   where va.variation_id=v.variation_id  and va.sample_id=$pid and user_name="$user_name" and idacmg=validation order by validation_id desc}; 
	}
	my $sth = $self->dbh->prepare($query) or confess();
	
	$sth->execute() || die();
	return $self->hresults($sth);
}

###################
# CNV METHODS
###################

sub get_exons {
	return {};
}
sub  interval_tree_cnv {
	my ($self,$sth) = @_;
	my $res = {};
	 my $tree = Set::IntervalTree->new;
	 my $hh;
	while (my $s =$sth->fetchrow_hashref()){
		my $id = $s->{type}."_".$s->{chromosome}."_".$s->{start}."_".$s->{end};
		$s->{$id} = $id;
		$hh->{$id}->{patient}->{$s->{patient_id}} ++;
		
		push(@{$hh->{$id}->{patient}->{$s->{patient_id}}},$s);
		push(@{$hh->{patient}->{$s->{patient_id}}},$s);
	}
	
	foreach my $k (keys %$hh){
		my ($type,$chr,$start,$end) =split("_",$k);
		$tree->insert($start,$end,$hh->{$k});
	}
	return $tree;
}

sub  hresults_cnv {
	my ($self,$sth) = @_;
	my $res = {};
	while (my $s =$sth->fetchrow_hashref()){
		my $id = $s->{type}."_".$s->{chromosome}."_".$s->{start}."_".$s->{end};
		push(@{$res->{$id}},$s);
	}
	return $res;
}


sub getAllCnvsValidationsForPatient {
	my ($self,$patient,$user_name) = @_;
	my $db = $self->db;
	my $pid = $patient->id;
	my $query = qq{select * from  $db.cnv_validations as va ,$db.cnvs as v,$db.acmg as acmg     where v.id=va.cnv_id and va.validation=acmg.idacmg and va.sample_id=$pid order by validation_id desc}; 
	if ($user_name){
		 $query = qq{select * from $db.cnv_validations as va ,$db.cnvs as v ,$db.acmg as acmg     where v.id=va.cnv_id and va.validation=acmg.idacmg and va.sample_id=$pid and user_name="$user_name" order by validation_id desc}; 
	}
	my $sth = $self->dbh->prepare($query) ;
	
	$sth->execute() || die();
	return $self->hresults_cnv($sth);
}

sub getCNVId{
		my ($self,$hcnv) = @_; 
	my $db = $self->db;
	my $chr = $hcnv->{chromosome};
	my $start = $hcnv->{start};
	my $end = $hcnv->{end};
	my $type = $hcnv->{type};
	my $query = qq{
		select id   from $db.cnvs where chromosome = "$chr" and start=$start and end = $end   and type = "$type" and version = 'HG19';
	};
	my $sth = $self->dbh->prepare($query);
	$sth->execute();
	my $s = $sth->fetchrow_hashref();	
	return $s->{id};
}


sub insertCNV {
	my ($self,$hcnv) = @_; 
	my $chr = $hcnv->{chromosome};
	my $start = $hcnv->{start};
	my $end = $hcnv->{end};
	my $type = $hcnv->{type};
	 my $db = $self->db;
	my $query = qq{
		insert into $db.cnvs (chromosome,start,end,type) values ("$chr",$start,$end,"$type");
	};
	$self->dbh->do($query) || die($query);
	
	return $self->getCNVId($hcnv);
}

sub createCNV {
	my ($self,$hcnv) = @_;
	my $id = $self->getCNVId($hcnv);
	unless ($id){
		$id = $self->insertCNV($hcnv);
	}
	confess() unless $id;
	return $id;
}


sub createValidationCNV {
	my ($self,$hcnv,$patient,$user_name,$validation_type) = @_;
	my $cnv_id = $self->createCNV($hcnv);
	#$hcnv = {chr_name} {start} {end} {type}
	my $project_name =$patient->project->name();
	my $project_id = $patient->project->id();
	my $phenotypes_id = join(';', map {$_->id()} @{$patient->project->getPhenotypes()} );
	my $phenotypes_name = join(';', map {$_->name()} @{$patient->project->getPhenotypes()} );
	my $patient_name = $patient->name;
	my $patient_id = $patient->id;
	my $i =  $hcnv->{infos}->{$patient->id};
	$i = {} unless $i;
	 my $infos = encode_json $i;
	 my $db = $self->db;
	my $query = qq{
		insert into $db.cnv_validations (cnv_id,project_name,project_id,sample_name,sample_id, user_name,infos,creation_date,modification_date,validation,phenotypes_id,phenotypes_name) values
		                                             ($cnv_id, "$project_name",$project_id,"$patient_name",$patient_id,"$user_name",'$infos',NOW(),NOW(),$validation_type,"$phenotypes_id","$phenotypes_name") ON DUPLICATE KEY UPDATE modification_date = NOW()  ;
	};

	$self->dbh->do($query) || die($query);
	$self->addPatientsStatus($patient,$user_name,$validation_type);
}

##########################
# Patient status Methods
#############################

sub addPatientsStatus {
	my($self,$patient,$user_name,$value) =@_;
	 my $patient_name = $patient->name;
	 my $patient_id = $patient->id;
	 my $project_name =$patient->project->name();
		my $project_id = $patient->project->id();
	 my $db = $self->db;
	 my $query = qq{
		insert into $db.samples_status (project_name,project_id,sample_name,sample_id, user_name,creation_date,modification_date,status) values
		                                             ("$project_name",$project_id,"$patient_name",$patient_id,"$user_name",NOW(),NOW(),$value) ON DUPLICATE KEY UPDATE modification_date = NOW() ;
		};
	$self->dbh->do($query) || die($query);
}

sub getLatestStatusPatientsForUser {
	my ($self,$patient,$user_name) = @_;
	confess();
	 my $db = $self->db;
	# my $pid = $patient->project->id;
	 my $sid = $patient->id;
	 my $query = qq{
	 	SELECT sample_id,term,creation_date as date  FROM $db.samples_status,$db.acmg
	 	 where sample_id=$sid and user_name="$user_name" and idacmg=status order by modification_date desc limit 1 ;
				};
	my $sth = $self->dbh->prepare($query) ;
	$sth->execute() || die();
	return $sth->fetchrow_hashref();
	
}

sub getLatestStatusPatients {
	my ($self,$patient,$user_name) = @_;
	 my $db = $self->db;
	# my $pid = $patient->project->id;
	 my $sid = $patient->id;
	 my $query = qq{
	 	SELECT sample_id,term,creation_date as date,id,modification_date,UNIX_TIMESTAMP(modification_date) as unix_time, status,user_name   FROM $db.samples_status,$db.acmg
	 	 where sample_id=$sid  and idacmg=status order by unix_time desc ;
				};
		if ($user_name){
			$query = qq{
	 	SELECT sample_id,term,creation_date as date,id,modification_date,UNIX_TIMESTAMP(modification_date) as unix_time, status,user_name   FROM $db.samples_status,$db.acmg
	 	 where sample_id=$sid  and idacmg=status and user_name="$user_name" order by unix_time desc ;
				};
		}		
	my $sth = $self->dbh->prepare($query) ;
	$sth->execute() || die();
	return $self->hByPatients($sth);
}

sub getLatestStatusProject {
	my ($self,$project_id) = @_;
	 my $db = $self->db;
	# my $pid = $patient->project->id;
	 my $query = qq{
	 	SELECT sample_id,term,creation_date as date,id,modification_date,UNIX_TIMESTAMP(modification_date) as unix_time, status,user_name   FROM $db.samples_status,$db.acmg
	 	 where project_id=$project_id  and idacmg=status order by unix_time desc ;
				};
	my $sth = $self->dbh->prepare($query) ;
	$sth->execute() || die();
	return $self->hByPatients($sth);
}

1;
