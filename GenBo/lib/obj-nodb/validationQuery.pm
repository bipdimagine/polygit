package validationQuery;
use strict;
use Data::Dumper;
use Carp;
use Moo;

use IO::Compress::Gzip qw(gzip $GzipError) ;

has 'dbh' => (
	is =>'ro',
	
	#weaken=>1,
	required => 1,
);
 
has 'db' =>(
	is =>'ro',	
	lazy=>1,
	default =>sub {
		my $self = shift;
		return "validation_OI" unless $self->capture_name();
		my $dbname = "validation_".$self->capture_name();
		my $query = qq{SHOW DATABASES LIKE '$dbname'};
		my $sth = $self->dbh->prepare($query);
		$sth->execute();
		my $s = $sth->fetchrow_hashref();
		return "validation_".$self->capture_name() if $s;
		confess("no database validation for this capture validation_".$self->capture_name());
		return "validation_OI" ;
	}
);
has 'db2' =>(
	is =>'ro',	
	lazy=>1,
	default =>sub {
		my $self = shift;
		return "validation_OI" unless $self->capture_name();
		my $dbname = "validation_".$self->capture_name();
		my $query = qq{SHOW DATABASES LIKE '$dbname'};
		my $sth = $self->dbh->prepare($query);
		$sth->execute();
		my $s = $sth->fetchrow_hashref();
		return "validation_".$self->capture_name() if $s;
		confess("no database validation for this capture validation_".$self->capture_name());
		return "validation_OI" ;
	}
);
has 'capture_name' =>(
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

has variations=>(
is =>'ro',	
	lazy=>1,
	default =>sub {
	my ($self) =@_;
	my $db = $self->db;
	warn $db;
	my $query = qq{
		select variation_id as query_id, polyid as genbo_id,vcfid as vcf_id from $db.variations ;
	};
	my $sth = $self->dbh->prepare($query);
	$sth->execute();
	my $s = $sth->fetchrow_hashref("genbo_id");	
	my $vs;
	while ( my $v  = $sth->fetchrow_hashref() ){
		$vs->{$v->{genbo_id}} = $v->{query_id};
		$vs->{$v->{vcf_id}} = $v->{query_id};
	}
	return $vs;
	}
);


sub isValidated {
	my ($self,$id) =@_;
	return $self->variations->{$id}  if exists $self->variations->{$id};
	return;
#	return undef:
}


sub getVariation {
	my ($self,$id) =@_;
	my $db = $self->db;
	my $query = qq{
		select variation_id as id from $db.variations where polyid='$id';
	};
	my $sth = $self->dbh->prepare($query);
	$sth->execute();

	my $s = $sth->fetchrow_hashref();	
	return $s->{id};
}

sub getVariationByVcfId  {
	my ($self,$id) =@_;
	my $db = $self->db;
	my $query = qq{
		select variation_id as id  from $db.variations where vcfid='$id';
	};
	my $sth = $self->dbh->prepare($query);
	$sth->execute();
	my $s = $sth->fetchrow_hashref();	
	return $s->{id};
}

sub getVariationByGenBoId {
	my ($self,%arg) =@_;
	my $id = $arg{id};
	my $db = $self->db;
	my $id2 = "chr".$id;
	my $query = qq{
		select variation_id as id  from $db.variations where polyid='$id' or vcfid ='$id' or vcfid = "$id2";
	};
	my $sth = $self->dbh->prepare($query);
	$sth->execute();
	my $s = $sth->fetchrow_hashref();	
	return $s->{id};
}
sub getValidations{
	my ($self,%arg) =@_;
	my $db = $self->db;
	my $id = $arg{id};
	my $project = $arg{project};
	my $sample = $arg{sample};
	my $query = qq{select validation as validation ,validation_sanger as sanger from $db.validations where variation_id=$id and project_name='$project'  and sample_name='$sample' order by modification_date; };
	my $sth = $self->dbh->prepare($query);
	$sth->execute();
	my $s = $sth->fetchrow_hashref();	
	
	return ($s->{validation},$s->{sanger});
	
}


sub getIdFromPositionSequence{
		my ($self,$chr,$start,$end,$seq,$type) = @_; 
	my $db = $self->db;
	my $query = qq{
		select variation_id as id  from $db.variations where chromosome = $chr and start=$start and end = $end and sequence = "$seq"  and type = "$type" and version = 'HG19';
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
	
	return $self->hresults($sth);
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
sub get_variations_validated{
	my ($self,%arg) = @_;
	my $project_name =$arg{project_name};
	my $sample_name = $arg{sample_name};
	my $db = $self->db;
	my $query = qq{SELECT v.validation_id as validation_id,  vr.vcfid as vcfid, validation, validation_sanger FROM $db.validations v, $db.variations vr where
v.variation_id=vr.variation_id and v.sample_name="$sample_name" and v.project_name="$project_name";};
	my $sth = $self->dbh->prepare($query);
	$sth->execute();
	#my $s = $sth->fetchall_arrayref();	
	my $variations;
	while ( my $v  = $sth->fetchrow_hashref() ){
	#foreach my $v (@{$sth->fetchall_arrayref()}){
		
		my $var;
		($var->{chromosome},$var->{start},$var->{ref},$var->{alt}) = split("_",$v->{vcfid});
		$var->{validation} = $v->{validation};
		$var->{validation_id} = $v->{validation_id};
		$var->{validation_sanger} = $v->{validation_sanger} if $v->{validation_sanger};
		
		push(@$variations,$var);
	}
	return $variations;
}

sub get_variations_ions {
	confess();
#	method get_variations_ions(Str :$project_name!, Str :$sample_name!){
#	my $db = $self->db;
#	my $query = qq{SELECT v.validation_id as validation_id,  vr.vcfid as vcfid, validation, validation_sanger FROM $db.validations v, $db.variations vr where
#v.variation_id=vr.variation_id and v.sample_name="$sample_name" and v.project_name="$project_name" and v.validation_sanger=0 and v.validation!=-3;};
#	my $sth = $self->dbh->prepare($query);
#	$sth->execute();
#	#my $s = $sth->fetchall_arrayref();	
#	my $variations;
#	while ( my $v  = $sth->fetchrow_hashref() ){
#		my $var;
#		($var->{chromosome},$var->{start},$var->{ref},$var->{alt}) = split("_",$v->{vcfid});
#		$var->{vcf_id} = $v->{vcfid};
#		$var->{validation} = $v->{validation};
#		$var->{validation_id} = $v->{validation_id};
#		$var->{validation_sanger} = $v->{validation_sanger} if $v->{validation_sanger};
#		$variations->{$v->{vcfid}} = $var;
#		#push(@$variations,$var);
#	}
#	return $variations;
#}
}

sub get_variations_in_validation_table {
	my ($self,%arg) = @_;
	my $project_name = $arg{project_name};
	my $sample_name= $arg{sample_name};
	my $db = $self->db;
	my $query = qq{SELECT * FROM $db.validations v, $db.variations vr where v.variation_id=vr.variation_id 
		and v.sample_name="$sample_name" and v.project_name="$project_name" order by modification_date;};
    #my $query = qq{SELECT * FROM $db.validations v, $db.variations vr where v.variation_id=vr.variation_id and v.sample_name="$sample_name" and v.project_name="$project_name" and  v.validation_sanger>0 ;};
	my $sth = $self->dbh->prepare($query);
	$sth->execute();
	
	#my $s = $sth->fetchall_arrayref();	
	my $variations;
	while ( my $v  = $sth->fetchrow_hashref() ){
	#foreach my $v (@{$sth->fetchall_arrayref()}){
		delete $v->{sam_lines};
		$variations->{$v->{vcfid}} = $v; 
		$variations->{$v->{polyid}} = $v; 
		
	}
	return $variations;
}

sub get_validations_users {
		my ($self,%arg) = @_;
	my $project_name = $arg{project_name};
	my $sample_name = $arg{sample_name};
	my $db = $self->db;
	my $query = qq{SELECT user_name as user ,modification_date as date FROM $db.validations v, $db.variations vr where
		 v.variation_id=vr.variation_id and v.sample_name="$sample_name" and v.project_name="$project_name" 
		 group by user_name order by modification_date;};
		 my $sth = $self->dbh->prepare($query);
		$sth->execute();
		my $variations;
		my $last_date;
		my %users;
		my $last_user;
	while ( my $v  = $sth->fetchrow_hashref() ){
	#foreach my $v (@{$sth->fetchall_arrayref()}){
		$users{$v->{user}} ++;
		$last_user = $v->{user};
		$last_date = $v->{date};
	
		
	}
	return (\%users,$last_user,$last_date);
		 
}
sub get_variations_todo{
	my ($self,%arg) = @_;
	my $project_name = $arg{project_name};
	my $sample_name = $arg{sample_name};
	my $uniq  = $arg{uniq};
	my $db = $self->db;
	my $query = qq{SELECT * FROM $db.validations v, $db.variations vr where v.variation_id=vr.variation_id and v.sample_name="$sample_name" and v.project_name="$project_name" and v.validation_sanger=0 and v.validation =-3;};
    #my $query = qq{SELECT * FROM $db.validations v, $db.variations vr where v.variation_id=vr.variation_id and v.sample_name="$sample_name" and v.project_name="$project_name" and  v.validation_sanger>0 ;};
	my $sth = $self->dbh->prepare($query);
	$sth->execute();
	#my $s = $sth->fetchall_arrayref();	
	my $variations;
	while ( my $v  = $sth->fetchrow_hashref() ){
	#foreach my $v (@{$sth->fetchall_arrayref()}){
		delete $v->{sam_lines};
			($v->{chromosome},$v->{start},$v->{ref},$v->{alt}) = split("_",$v->{vcfid});
		$variations->{$v->{vcfid}} = $v unless $uniq; 
		$variations->{$v->{polyid}} = $v; 
		
	}
	return $variations;
}
sub get_variations_sanger {
	my ($self,%arg) = @_;
	my $project_name = $arg{project_name};
	my $sample_name = $arg{sample_name};
	my $db = $self->db;
	my $query = qq{SELECT v.validation_id as validation_id, v.user_name as user,  vr.vcfid as vcfid, validation, validation_sanger FROM $db.validations v, $db.variations vr where
v.variation_id=vr.variation_id and v.sample_name="$sample_name" and v.project_name="$project_name" and v.validation_sanger != 0;};
	my $sth = $self->dbh->prepare($query);
	$sth->execute();
	#my $s = $sth->fetchall_arrayref();	
	my $variations;
	while ( my $v  = $sth->fetchrow_hashref() ){
		#foreach my $v (@{$sth->fetchall_arrayref()}){
		
		my $var;
		($var->{chromosome},$var->{start},$var->{ref},$var->{alt}) = split("_",$v->{vcfid});
		$var->{vcf_id} = $v->{vcfid};
		$var->{validation} = $v->{validation};
			$var->{user} = $v->{user};
		$var->{validation_id} = $v->{validation_id};
		$var->{validation_sanger} = $v->{validation_sanger} if $v->{validation_sanger};
		$variations->{$v->{vcfid}} = $var;
		
	}
	return $variations;
}
sub get_exons {
	my ($self,%arg) = @_;
	my $project_name = $arg{project_name};
	my $sample_name = $arg{sample_name};
	my $db = $self->db;
	my $query = qq{select * from $db.exons where  project_name="$project_name" and sample_name ="$sample_name";};
	my $sth = $self->dbh->prepare($query) ;
	
	$sth->execute() || die();
	my $exons;
	while ( my $v  = $sth->fetchrow_hashref() ){
	#foreach my $v (@{$sth->fetchall_arrayref()}){
		my $id = join("-",$v->{chromosome},$v->{transcript},$v->{start},$v->{end});
		$exons->{$id} = $v;
	
	}
	return $exons;
	}
	
	
sub getTodoVariations {
	my ($self,%arg) = @_;
	my $project_name = $arg{project_name};
	my $sample_name = $arg{sample_name};
	my $db = $self->db;
my $sql = qq{SELECT validation_id ,sample_name as patient_name ,vcfid as variation_name ,validation  FROM $db.validations v, $db.variations vr where v.project_name="NGS2013_0201" and v.variation_id=vr.variation_id; };

my $sth = $self->dbh->prepare($sql);
$sth->execute();
my $hres = $sth->fetchall_hashref("validation_id");
my $validations;
foreach my $id (keys %$hres){
	my $var;
	($var->{chr},$var->{start},$var->{ref},$var->{alt}) = split("_",$hres->{$id}->{variation_name});
	$var->{chr_name} =~ s/chr//;
	$var->{validation} = $hres->{$id}->{validation};
	push(@{$validations->{$hres->{$id}->{patient_name}}},$var);
	
}

}
sub existsValidation {
	my ($self,%arg) = @_;
	my $project = $arg{project};
	my $sample = $arg{sample};
	my $variation_id = $arg{variation_id};
	my $user =  $arg{user};
	
	my $db = $self->db;
	my $query = qq{
		select validation_id as id  from $db.validations where variation_id=$variation_id and project_name='$project' and user_name='$user' and sample_name='$sample' ;
	};

	my $sth = $self->dbh->prepare($query);
	$sth->execute();
	my $s = $sth->fetchrow_hashref();
	return $s->{id};
}
sub createValidation {
		my ($self,%arg) = @_;
	my $project = $arg{project};
	my $sample = $arg{sample};
	my $variation_id = $arg{variation_id};
	my $user =  $arg{user};
	my $vcf_line = $arg{vcf_line};
	my $validation_status = $arg{validation_status};
	my $bam_line = $arg{bam_line};
	my $method = $arg{method};

	my $db = $self->db;
	my $validation_id = $self->existsValidation(variation_id=>$variation_id,project=>$project,sample=>$sample,user=>$user);

	if ($validation_id){
		my $query = qq{
		update $db.validations set vcf_line=?, modification_date=NOW(),sam_lines=COMPRESS(?),validation=?,method=?,validation_sanger=? where validation_id = $validation_id;
	};
	my $sth= $self->dbh->prepare($query);
	$sth->execute($vcf_line,$bam_line,$validation_status,$method,0);
	$sth->finish;
	return 1;
	}
	else {								 
	my $query = qq{
		insert into $db.validations (variation_id,project_name,sample_name,user_name,validation,vcf_line,creation_date,modification_date,sam_lines,method,project_id,validation_sanger) values (?,?,?,?,?,?,NOW(),NOW(),COMPRESS(?),?,?,?);
	};
	my $sth= $self->dbh->prepare($query);
	$sth->execute($variation_id,$project,$sample,$user,$validation_status,$vcf_line,$bam_line,$method,0,0);
	$sth->finish;
	#$self->dbh->do($query) || return undef;
	return 1;
	}
}
sub update_variation_validation {
		my ($self,%arg) = @_;
	my $heho = $arg{heho};
	my $validation_id = $arg{validation_id};
	my $method = $arg{method};
	my $db = $self->db;
	my $query = qq{update $db.validations set validation_sanger=$heho , modification_date=NOW() where validation_id=$validation_id};
	#my $query = qq{insert into $db.sanger_validations (validation_id,validation,method,creation_date) values(?,?,?,NOW()) };
	warn $query;
	$self->dbh->do($query);
	
	#$self->dbh->do($query) || return undef;
	return 1;
}
sub update_exon_validation {
	my ($self,%arg) = @_;
my $exon_id = $arg{exon_id};
my $value = $arg{value};

	my $db = $self->db;
	
	my $query = qq{update $db.exons set done=$value , modification_date=NOW() where exon_id="$exon_id"};
	$self->dbh->do($query);
	
	
	return 1;
}
sub compressData{
	my ($data) = @_;
	my $data2;
	 gzip \$data => \$data2
        or die "gzip failed: $GzipError\n";
    return $data2;
}
sub exon_todo {
	confess();
#method exon_todo (Str :$project_name!, Str :$id!,Str :$sample_name, Str :$chromosome!, Str :$start!, Str :$end!, Str :$transcript!, Str :$user_name, Int :$todo!, Str :$name, Str :$gene! ){ 
#	#$bam_line = "toto";
#	my $db = $self->db;
#	my $query = qq{
#		insert into $db.exons (exon_id , project_name,sample_name,gene,chromosome,start,end,transcript,user_name,todo,creation_date,modification_date,done) 
#		                       values("$id","$project_name","$sample_name","$gene","$chromosome",$start,$end,"$transcript","$user_name",$todo,NOW(),NOW(),0 )
#		on DUPLICATE KEY UPDATE todo=$todo, modification_date= NOW();
#	};
#	warn $query;
#	$self->dbh->do($query) || die($query);
#	return 1;
	}
sub is_todo{
	confess();
#method is_todo (Str :$project_name!, Str :$id!){ 
#	#$bam_line = "toto";
#	my $db = $self->db;
#	my $query = qq{select * from $db.exons where exon_id="$id" and project_name="$project_name";};
#	my $sth = $self->dbh->prepare($query) ;
#	
#	$sth->execute() || die();
#	my $s = $sth->fetchrow_hashref();
#	 if (exists $s->{exon_id} && !$s->{done}){
#	 	return (1);
#	 }
#	return;
	}
sub save_report{
	confess();
#method save_report (Str :$project!, Str :$sample!,Str :$conclusion!, Str :$json, Str :$user_name){ 
#	
#	my $db = $self->db;
#	my $query =qq{insert into $db.reports(project,sample,creation_date,json,conclusion,user_name)  values (?,?,NOW(),?,?,?)
#		on DUPLICATE KEY UPDATE creation_date=NOW(),json=?,conclusion=?;
#	};
#	my $sth = $self->dbh->prepare($query) ;
#	$sth->bind_param(1,$project);
#	$sth->bind_param(2,$sample);
#	$sth->bind_param(3,$json);
#	$sth->bind_param(4,$conclusion);
#	$sth->bind_param(5,$user_name);
#	$sth->bind_param(6,$json);
#	$sth->bind_param(7,$conclusion);
#	
#	$sth->execute() || die();
#	#$self->dbh->do($query) || die($query);
#	return 1;
#}
}
sub get_report{
	confess();
#method get_report (Str :$project!, Str :$sample!){ 
#	my $db = $self->db;
#	my $query =qq{select * from $db.reports where project="$project" and sample="$sample"};
#	my $sth = $self->dbh->prepare($query) ;
#	$sth->execute() || die();
#	my $s = $sth->fetchrow_hashref();
#	
#	return $s;
#	return 1;
}

sub getVariantValidated{
	my ($self,$project,$varid) =@_;
	my $db = $self->db;
	my $query = qq{select * from $db.variations var, $db.validations val 
		where var.variation_id = val.variation_id and val.project_name != '$project' and var.polyid = '$varid' and (val.validation>0 or val.validation_sanger>0);
	};
	my $sth = $self->dbh->prepare($query);
	$sth->execute();
	my $s = $sth->fetchrow_hashref();
	my $row = $sth->rows;	
	return $row;
}

sub getVariantRejected{
	my ($self,$project,$varid) =@_;
	my $db = $self->db;
	my $query = qq{select * from $db.variations var, $db.validations val 
		where var.variation_id = val.variation_id and val.project_name != '$project' and var.polyid = '$varid' and (val.validation_sanger < 0 and val.validation_sanger != -3);
	};
	my $sth = $self->dbh->prepare($query);
	$sth->execute();
	my $s = $sth->fetchrow_hashref();
	my $row = $sth->rows;	
	return $row;
}
sub createVariation {
	my ($self,%arg) = @_;
	my $polyid = $arg{polyid};
	my $vcfid = $arg{vcfid};
	my $genboid = $arg{genboid};
	my $version = $arg{version};
	unless ($vcfid =~ /chr/){
		$vcfid = "chr".$vcfid;
		$vcfid="chrM" if $vcfid eq "chrMT";
		
	}
	my $id = $self->getVariationByVcfId($vcfid);
	return $id if $id;
	my $test_vcf_id = $vcfid;
	if ($test_vcf_id =~ /chr/){
		$test_vcf_id =~ s/chr//;
		$test_vcf_id= "MT" if $test_vcf_id eq "M";
	}
	else {
		die();
	}
	 $id = $self->getVariationByVcfId($test_vcf_id);
	return $id if $id;
	my $db = $self->db;
	my $query = qq{
		insert into $db.variations (polyid,vcfid,genbo_id,version) values ('$polyid','$vcfid','1','$version');
	};
	
	$self->dbh->do($query) || return undef;
	return $self->getVariationByVcfId($vcfid);
	
}
sub getAllCnvsValidationsForPatient {
	
	return{};
}
sub getAllValidationsForPatient {
	my ($self,$patient,$user_name) = @_;
	my $db = $self->db;
	my $pid = $patient->name;
	my $query = qq{select * from $db.validations as va ,$db.variations as v   where va.variation_id=v.variation_id  and va.sample_name="$pid" order by validation_id desc}; 
	if ($user_name){
		$query = qq{select * from $db.validations as va ,$db.variations as v   where va.variation_id=v.variation_id  and va.sample_name="$pid" and user_name="$user_name" order by validation_id desc}; 
	}
	my $sth = $self->dbh->prepare($query) ;
	$sth->execute() || die();
	return $self->hresults($sth);
}

sub getAllValidations {
	my ($self) = @_;
	my $db = $self->db;
	my $query = qq{select * from $db.validations as va ,$db.variations as v   where va.variation_id=v.variation_id  order by validation_id desc}; 
	my $sth = $self->dbh->prepare($query) ;
	$sth->execute() || die();
	return $self->hresults($sth);
}

sub getLatestStatusPatientsForUser {
	my ($self,$patient,$user_name) = @_;
	return {};
	
}
sub getLatestStatusPatients {
	return {};
}
1;