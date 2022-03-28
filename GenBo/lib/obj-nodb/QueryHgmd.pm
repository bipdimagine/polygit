package QueryHgmd;
use strict;
use Vcf;
use Moose;
use MooseX::Method::Signatures;
use Data::Dumper;
use Config::Std;
use JSON;
use GBuffer;

has config => (
	is		=> 'ro',
	isa		=> 'HashRef',
	reader	=> 'getConfig',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $config = $self->all_config;
		return $$config{polyprod};
	},
);

has build => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return 'HG19';
	}
);

has database => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		return $self->last_database_name();
	}
);

has last_database_name => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $buffer = GBuffer->new();
		my $h = $buffer->config->{hgmd_db_current_version};
		my @lKeys = keys %$h;
		confess("\n\nERROR: hgmd_db_current_version must have one key / value only in config file. Die...\n\n") if (scalar(@lKeys) > 1);
		return $lKeys[0];
	}
);

has previous_database => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $buffer = GBuffer->new();
		my $h = $buffer->config->{hgmd_db_previous_version};
		my @lKeys = keys %$h;
		confess("\n\nERROR: hgmd_db_previous_version must have one key / value only in config file. Die...\n\n") if (scalar(@lKeys) > 1);
		return $lKeys[0];
	}
);

has current_database_name => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $buffer = GBuffer->new();
		return $buffer->config->{hgmd_db_current_version}->{$self->database()};
	}
);

has getHashOldDatabases  => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $buffer = GBuffer->new();
		my $h = $buffer->config->{hgmd_db_old_versions};
		return $h;
	}
);

has all_config => (
	is		=> 'ro',
	#isa		=> 'HashRef',
	required=> 1,
);

has dbh => (
	is		=> 'ro',
	isa		=> 'DBI::db',
	required=> 1,
	reader	=> 'getDbh',
);

has hashAccNum => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $buffer = GBuffer->new();
		my $h;
		my @list_acc_num = $self->getAllAccNum($self->database(), 'DM');
		foreach my $acc_num (@list_acc_num) {
			$h->{$acc_num}->{$self->database()}++;
		}
		return $h;
	}
);

has hashOldAccNum => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $buffer = GBuffer->new();
		my $h;
		foreach my $hgmd_version (keys %{$self->getHashOldDatabases()}) {
			my $value = $self->getHashOldDatabases->{$hgmd_version};
			my @list_acc_num = $self->getAllAccNum($hgmd_version, 'DM');
			foreach my $acc_num (@list_acc_num) {
				$h->{$acc_num}->{$value}++;
			}
		}
		return $h;
	}
);

sub getAllAccNum {
        my ($self, $database, $class) = @_;
        my $dbh =  $self->getDbh;
        my $sql = qq{ SELECT acc_num FROM `$database`.allmut };
        if ($class) { $sql .= qq{ where tag='$class'; }; }
        else { $sql .= ';'; }
        my $sth = $dbh->prepare($sql);
        $sth->execute();
        return keys %{$sth->fetchall_hashref("acc_num")};
}

sub getAllAccNum_with_hg19_coord {
        my ($self, $database, $class) = @_;
        my $dbh =  $self->getDbh;
        my $sql = qq{ SELECT a.acc_num FROM `$database`.allmut as a, `$database`.hg19_coords as c where a.acc_num=c.acc_num };
        if ($class) { $sql .= qq{ and a.tag='$class'; }; }
        else { $sql .= ';'; }
        my $sth = $dbh->prepare($sql);
        $sth->execute();
        return keys %{$sth->fetchall_hashref("acc_num")};
}

sub getDataHGMDPro {
	my ($self,$id) = @_;
	my $database = $self->database();
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT * FROM `$database`.allmut where acc_num = ?;};
	my $sth = $dbh->prepare($sql);
	$sth->execute($id);
	my $toto = $sth->fetchall_hashref("acc_num");
	$sql = qq{ SELECT * FROM `$database`.extrarefs where acc_num = ?;};
	$sth = $dbh->prepare($sql);
	$sth->execute($id);
	my $titi = $sth->fetchall_hashref("pmid");
	$toto->{$id}->{pubmed} = $titi; 
	return $toto->{$id};
}

sub getDataHGMDPro_gene_infos {
	my ($self, $gene_name) = @_;
	my $database = $self->database();
	confess("\n\nERROR: need a gene_name for QueryMoose::getDataHGMDPro_gene_infos. Die.\n\n") unless ($gene_name);
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT * FROM `$database`.allgenes where gene=?; };
	my $sth = $dbh->prepare($sql);
	$sth->execute("$gene_name");
	my $toto = $sth->fetchall_hashref("gene");
	return $toto;
}

sub getDataHGMDPro_from_gene_name {
	my ($self, $gene_name, $class_name) = @_;
	my $database = $self->database();
	confess("\n\nERROR: need a gene_name for QueryMoose::getDataHGMDPro_from_gene_name. Die.\n\n") unless ($gene_name);
	confess("\n\nERROR: need a class_name for QueryMoose::getDataHGMDPro_from_gene_name. Die.\n\n") unless ($class_name);
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT gene FROM `$database`.allmut where gene=? and tag=?; };
	my $sth = $dbh->prepare($sql);
	$sth->execute($gene_name, $class_name);
	my $toto = $sth->fetchall_hashref("gene");
	return $toto;
}

sub getDataHGMDPro_positions_for_acc_num {
	my ($self, $acc_num) = @_;
	my $database = $self->database();
	confess("\n\nERROR: need a acc_num for QueryMoose::getDataHGMDPro_positions_for_acc_num. Die.\n\n") unless ($acc_num);
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT * FROM `$database`.hg19_coords where acc_num=?; };
	my $sth = $dbh->prepare($sql);
	$sth->execute($acc_num);
	my $h = $sth->fetchall_hashref("acc_num");
	return $h;
}

sub getDataHGMDPro_positions_for_class {
	my ($self, $chr_name, $class_name) = @_;
	my $database = $self->database();
	confess("\n\nERROR: need a chr_name for QueryMoose::getDataHGMDPro_positions_for_class. Die.\n\n") unless ($chr_name);
	confess("\n\nERROR: need a class_name for QueryMoose::getDataHGMDPro_positions_for_class. Die.\n\n") unless ($class_name);
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT hg19.coordSTART FROM `$database`.allmut as allmut, `$database`.hg19_coords as hg19 where allmut.tag='DM' and allmut.chromosome=? and allmut.acc_num=hg19.acc_num and hg19.coordSTART; };
	my $sth = $dbh->prepare($sql);
	$sth->execute($chr_name);
	my $toto = $sth->fetchall_hashref("coordSTART");
	my @lPos = keys %$toto;
	return \@lPos;
}

sub getDataHGMDPro_genes_for_class {
	my ($self, $chr_name, $class_name) = @_;
	my $database = $self->database();
	confess("\n\nERROR: need a chr_name for QueryMoose::getDataHGMDPro_genes_for_class. Die.\n\n") unless ($chr_name);
	confess("\n\nERROR: need a class_name for QueryMoose::getDataHGMDPro_genes_for_class. Die.\n\n") unless ($class_name);
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT gene FROM `$database`.allmut where chromosome=? and tag=?; };
	my $sth = $dbh->prepare($sql);
	$sth->execute($chr_name, $class_name);
	return $sth->fetchall_hashref("gene");
}

sub getHGMD {
	my ($self,$login)=@_;
	my $dbh = $self->getDbh;
	
		my $sql = qq{	  	
			select BU.hgmd as hmd
			from   bipd_users.USER BU 
			  where 
			  BU.LOGIN= ? ;
	};
my $sth = $dbh->prepare($sql);
	$sth->execute($login);
	my $toto = $sth->fetchall_arrayref();
	return  $toto->[0]->[0] if $toto;
	return 0;
}

sub search_gene_name {
	my ($self, $gene_name) = @_;
	my $database = $self->database();
	confess("\n\nERROR: need a gene_name for QueryMoose::search_gene_name. Die.\n\n") unless ($gene_name);
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT * FROM `$database`.allgenes where gene like ?; };
	my $sth = $dbh->prepare($sql);
	$gene_name = '%'.$gene_name.'%';
	$sth->execute($gene_name);
	return $sth->fetchall_hashref("gene");
}

sub search_acc_num {
	my ($self, $acc_num) = @_;
	my $database = $self->database();
	confess("\n\nERROR: need a acc_num for QueryMoose::search_acc_num. Die.\n\n") unless ($acc_num);
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT * FROM `$database`.allmut as a, `$database`.hg19_coords as c where a.acc_num=c.acc_num and a.acc_num like ?; };
	my $sth = $dbh->prepare($sql);
	$acc_num = '%'.$acc_num.'%';
	$sth->execute($acc_num);
	return $sth->fetchall_hashref("acc_num");
}

sub search_rsname {
	my ($self, $rsname) = @_;
	my $database = $self->database();
	confess("\n\nERROR: need a rsname for QueryMoose::search_rsname. Die.\n\n") unless ($rsname);
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT * FROM `$database`.allmut as a, `$database`.hg19_coords as c where a.acc_num=c.acc_num and a.dbsnp like ?; };
	my $sth = $dbh->prepare($sql);
	$rsname = '%'.$rsname.'%';
	$sth->execute($rsname);
	return $sth->fetchall_hashref("dbsnp");
}

sub search_refseq {
	my ($self, $refseq) = @_;
	my $database = $self->database();
	confess("\n\nERROR: need a refseq for QueryMoose::search_refseq. Die.\n\n") unless ($refseq);
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT * FROM `$database`.allmut where refseq like ?; };
	my $sth = $dbh->prepare($sql);
	$refseq = '%'.$refseq.'%';
	$sth->execute($refseq);
	return $sth->fetchall_hashref("refseq");
}

sub search_disease {
	my ($self, $disease) = @_;
	my $database = $self->database();
	confess("\n\nERROR: need a disease for QueryMoose::search_disease. Die.\n\n") unless ($disease);
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT disease, acc_num, tag FROM `$database`.allmut where disease like ?; };
	my $sth = $dbh->prepare($sql);
	$disease = '%'.$disease.'%';
	$sth->execute($disease);
	my $h = $sth->fetchall_hashref("acc_num");
	my $h_diseases;
	foreach my $acc_num (keys %$h) {
		my $disease = $h->{$acc_num}->{disease};
		$h_diseases->{$disease}->{$acc_num} = $h->{$acc_num};
	}
	return $h_diseases;
}

sub search_variant_for_description {
	my ($self, $description) = @_;
	my $database = $self->database();
	confess("\n\nERROR: need a description for QueryMoose::search_variant_for_description. Die.\n\n") unless ($description);
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT * FROM `$database`.allmut as a, `$database`.hg19_coords as c where a.acc_num=c.acc_num and (a.dbsnp like ? or a.disease like ? or a.gene like ? or a.acc_num like ?); };
	my $sth = $dbh->prepare($sql);
	$description = '%'.$description.'%';
	$sth->execute($description, $description, $description, $description);
	return $sth->fetchall_hashref("acc_num");
}

sub search_variant_for_disease {
	my ($self, $disease) = @_;
	my $database = $self->database();
	confess("\n\nERROR: need a disease for QueryMoose::search_variant_for_disease. Die.\n\n") unless ($disease);
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT * FROM `$database`.allmut as a, `$database`.hg19_coords as c where a.acc_num=c.acc_num and a.disease like ?; };
	my $sth = $dbh->prepare($sql);
	$sth->execute($disease);
	return $sth->fetchall_hashref("acc_num");
}

sub search_variant_for_gene_without_coord {
	my ($self, $gene_name) = @_;
	my $database = $self->database();
	confess("\n\nERROR: need a gene_name for QueryMoose::search_variant_for_gene. Die.\n\n") unless ($gene_name);
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT * FROM `$database`.allmut where gene=? and startCoord is null; };
	my $sth = $dbh->prepare($sql);
	$sth->execute($gene_name);
	return $sth->fetchall_hashref("acc_num");
}

sub search_variant_for_gene {
	my ($self, $gene_name) = @_;
	my $database = $self->database();
	confess("\n\nERROR: need a gene_name for QueryMoose::search_variant_for_gene. Die.\n\n") unless ($gene_name);
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT * FROM `$database`.allmut as a, `$database`.hgmd_hg19_vcf as vcf, `$database`.hg19_coords as c where a.acc_num=c.acc_num and a.gene=? and vcf.id=a.acc_num;; };
	my $sth = $dbh->prepare($sql);
	$sth->execute($gene_name);
	return $sth->fetchall_hashref("acc_num");
}

sub get_score_concept {
	my ($self, $gene_name,$concept) = @_;
	my $database = $self->database();
	confess("\n\nERROR: need a gene_name for QueryMoose::search_variant_for_gene. Die.\n\n") unless ($gene_name);
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT *  FROM `$database`.gene_concept_scores where gene=? and concept=?;};
	my $sth = $dbh->prepare($sql);
	$sth->execute($gene_name,$concept);
	return $sth->fetchall_hashref("gene");
}
sub get_allele_insertion_for_acc_num {
	my ($self, $acc_num) = @_;
	my $database = $self->database();
	confess("\n\nERROR: need a acc_num for QueryMoose::get_allele_insertion_for_acc_num. Die.\n\n") unless ($acc_num);
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT d.acc_num, d.inserted FROM `$database`.allmut as a, `$database`.delins as d where a.acc_num=d.acc_num and a.acc_num=?; };
	my $sth = $dbh->prepare($sql);
	$sth->execute($acc_num);
	my $h = $sth->fetchall_hashref("acc_num");
	return $h->{$acc_num}->{inserted};
}

sub get_allele_deletion_for_acc_num {
	my ($self, $acc_num) = @_;
	my $database = $self->database();
	confess("\n\nERROR: need a acc_num for QueryMoose::get_allele_deletion_for_acc_num. Die.\n\n") unless ($acc_num);
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT d.acc_num, d.deleted FROM `$database`.allmut as a, `$database`.delins as d where a.acc_num=d.acc_num and a.acc_num=?; };
	my $sth = $dbh->prepare($sql);
	$sth->execute($acc_num);
	my $h = $sth->fetchall_hashref("acc_num");
	return $h->{$acc_num}->{deleted};
}

sub get_hash_last_released_DM {
	my ($self) = @_;
	my $new_database = $self->database();
	my $buffer = GBuffer->new();
	my $version_db = $buffer->config->{hgmd_db_current_version}->{$self->database()};
	my $json_new = $self->json_hgmd_new_dm($version_db);
	
	# Comme cette requete est trop longue, le resultat est stocke lors de sa premiere execution dans son dossier annotation (derniere release connu)
	if (-e $json_new) {
		open (JSON, $json_new);
		my $json = <JSON>;
		close (JSON);
		my $h = decode_json $json;
		return $h;
	}
	my $old_database = $self->previous_database();
	my $dbh =  $self->getDbh;
	
	# Requete 1: HGMD DM dans la derniere release
	my $sql = qq{ SELECT vcf.id, vcf.chrom, vcf.pos, vcf.ref, vcf.alt, new.tag, new.gene, new.disease, new.expected_inheritance, new.hgvs, new.rankscore FROM `$new_database`.hgmd_hg19_vcf as vcf, `$new_database`.allmut as new 
	    LEFT JOIN `$old_database`.allmut as old 
	        ON old.acc_num=new.acc_num 
	            WHERE old.acc_num IS NULL and new.tag='DM' and vcf.id=new.acc_num; };
	my $sth = $dbh->prepare($sql);
	$sth->execute();
	my $h = $sth->fetchall_hashref("id");
	
	# Requete 2: Variants note nouveau HGMD DM dans historique derniere release.
	# ATTENTION: historique pas seulement depuis la release precedente MAIS historique global de TOUTES les releases
	my $sql2 = qq{ SELECT vcf.id, vcf.chrom, vcf.pos, vcf.ref, vcf.alt, new.tag, new.mutype, new.gene, new.new_date, new.disease, new.expected_inheritance, new.hgvs, new.rankscore  FROM `$new_database`.hgmd_hg19_vcf as vcf, `$new_database`.allmut as new,`$new_database`.history as h 
    					WHERE new.acc_num=h.acc_num and vcf.id=new.acc_num and h.afterUpd ="DM"; };
	my $sth2 = $dbh->prepare($sql2);
	$sth2->execute();
	my $h2 = $sth2->fetchall_hashref("id");
	foreach my $id (keys %$h2) { $h->{$id} = $h2->{$id}; }
	
	# Requete 3: Variants DM de la release precedente que l on pourra retirer de notre liste de nouveau DM
	my $sql3 = qq{ SELECT acc_num, tag FROM `$old_database`.allmut where tag="DM"; };
	my $sth3 = $dbh->prepare($sql3);
	$sth3->execute();
	my $h3 = $sth3->fetchall_hashref("acc_num");
	foreach my $id (keys %$h3) {
		if (exists $h->{$id}) {
			delete $h->{$id};
		}
	}
	
	open (JSON, '>'.$json_new);
	my $encode_json = encode_json $h;
	print JSON $encode_json;
	close (JSON);
	my $cmd = 'chmod 755 '.$json_new;
	`$cmd`;
	return $h;
}

sub get_hash_concepts {
	my $self = shift;
	
	my $new_database = $self->database();
	my $buffer = GBuffer->new();
	my $version_db = $buffer->config->{hgmd_db_current_version}->{$self->database()};
	my $old_database = $self->previous_database();
	my $dbh =  $self->getDbh;
	
	my $sql = qq{
		SELECT * FROM `$new_database`.gene_concept_scores where num_matching > 0;
	};
	my $sth = $dbh->prepare($sql);
	$sth->execute();
	return $sth->fetchall_hashref("concept");
}

sub get_hash_new_genes_with_concept {
	my $self = shift;
	
	my $new_database = $self->database();
	my $buffer = GBuffer->new();
	my $version_db = $buffer->config->{hgmd_db_current_version}->{$self->database()};
	my $old_database = $self->previous_database();
	my $dbh =  $self->getDbh;
	
	my $sql = qq{
		SELECT * FROM `$new_database`.gene_concept_scores where num_matching > 0;
	};
	my $sth = $dbh->prepare($sql);
	$sth->execute();
	my $h_now = $sth->fetchall_hashref("gene");
	
	
	
	my $sql2 = qq{
		SELECT * FROM `$old_database`.gene_concept_scores where num_matching > 0;
	};
	my $sth2 = $dbh->prepare($sql2);
	$sth2->execute();
	my $h_old = $sth2->fetchall_hashref("gene");
	
	my $h_new;
	foreach my $gene_name (keys %{$h_now}) {
		my $h_concepts_new = $self->get_gene_concept($gene_name, $new_database);
		if (exists $h_old->{$gene_name}) {
			my $h_concepts_old = $self->get_gene_concept($gene_name, $old_database);
			my $is_ok;
			foreach my $concept (keys %{$h_concepts_new}) {
				if (exists $h_concepts_old->{$concept}) {
					$h_concepts_new->{$concept}->{'new'} = 0;
				}
				elsif (exists $h_concepts_old->{'"'.$concept.'"'}) {
					$h_concepts_new->{$concept}->{'new'} = 0;
				}
				else {
					$h_concepts_new->{$concept}->{'new'} = 1;
				}
			}
			$h_new->{$gene_name} = $h_concepts_new;
		}
		else {
			$h_new->{$gene_name} = $h_concepts_new;
			foreach my $concept_name (keys %{$h_new->{$gene_name}}) {
				$h_new->{$gene_name}->{$concept_name}->{'new'} = 1;
			}
		}
	}
	return $h_new;
}

sub json_hgmd_new_dm {
	my ($self, $version) = @_;
	my $build = $self->build();
	my $buffer = GBuffer->new();
	unless ($version) {
		$version = $buffer->config->{hgmd_db_current_version}->{$self->database()};
	}
	my $json_file = $buffer->config->{public_data_annotation}->{root}.'/repository/'.$build.'/hgmd/'.$version.'/new_dm.json';
	return $json_file;
}

sub get_gene_infos {
	my ($self, $gene_name) = @_;
	return $self->{genes_infos}->{$gene_name} if (exists $self->{genes_infos} and $self->{genes_infos}->{$gene_name});
	my $database = $self->database();
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT * FROM `$database`.allgenes where gene=?; };
	my $sth = $dbh->prepare($sql);
	$sth->execute($gene_name);
	my $h = $sth->fetchall_hashref("gene");
	$self->{genes_infos}->{$gene_name} = $h->{$gene_name};
	return $h->{$gene_name};
}

sub get_gene_disease {
	my ($self, $gene_name) = @_;
	return $self->get_gene_infos($gene_name)->{disease};
}

sub get_gene_expected_inheritance {
	my ($self, $gene_name) = @_;
	return $self->get_gene_infos($gene_name)->{expected_inheritance};
}

sub get_gene_mut_total {
	my ($self, $gene_name) = @_;
	return $self->get_gene_infos($gene_name)->{mut_total};
}

sub get_gene_new_mut_total {
	my ($self, $gene_name) = @_;
	return $self->get_gene_infos($gene_name)->{new_mut_total};
}

sub get_gene_go_terms_acc {
	my ($self, $gene_name) = @_;
	return $self->get_gene_infos($gene_name)->{go_terms_acc};
}

sub get_gene_go_terms_name {
	my ($self, $gene_name) = @_;
	return $self->get_gene_infos($gene_name)->{go_terms_name};
}

sub get_gene_concept {
	my ($self, $gene_name, $database) = @_;
	$database = $self->database() unless ($database);
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT * FROM `$database`.gene_concept_scores where gene=?; };
	my $sth = $dbh->prepare($sql);
	$sth->execute($gene_name);
	my $h = $sth->fetchall_hashref("concept");
	return $h;
}

sub get_gene_complex_infos {
	my ($self, $gene_name) = @_;
	my $database = $self->database();
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT * FROM `$database`.complex where gene=?; };
	my $sth = $dbh->prepare($sql);
	$sth->execute($gene_name);
	my $h = $sth->fetchall_hashref("disease");
	return $h;
}

sub get_genes_from_concept_name {
	my ($self, $concept_name) = @_;
	my $database = $self->database();
	my $dbh =  $self->getDbh;
	my $sql = qq{
		SELECT gene FROM `$database`.gene_concept_scores
			where concept='?' and num_matching>0;
	};
	my $sth = $dbh->prepare($sql);
	$sth->execute($concept_name);
	my $h = $sth->fetchall_hashref('gene');
	return [keys %$h];
}

sub get_gene_infos_from_gene_name {
	my ($self, $gene_name) = @_;
	my $database = $self->database();
	my $dbh =  $self->getDbh;
	my $sql = qq{ SELECT * FROM `$database`.allgenes where gene=?; };
	my $sth = $dbh->prepare($sql);
	$sth->execute($gene_name);
	my $h = $sth->fetchall_hashref("gene");
	return $h->{$gene_name};
}

sub get_refseq_from_gene_name {
	my ($self, $gene_name) = @_;
	return $self->get_gene_infos_from_gene_name($gene_name)->{refseq};
}

sub get_disease_from_gene_name {
	my ($self, $gene_name) = @_;
	return $self->get_gene_infos_from_gene_name($gene_name)->{disease};
}

1;