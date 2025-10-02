=head1 NAME

connect : Specific GenBo for connection to the database 

=head1 SYNOPSIS


=head1 DESCRIPTION

connect provides a set of functions to get basics informations on the connection to the database

=head1 METHODS

=cut
 
package connect;

use strict;
use DBI; 
use DBD::mysql;


use Getopt::Long;
use Data::Dumper;
use Carp;
use IO::Uncompress::Gunzip qw(gunzip  $GunzipError) ;
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError) ;
  use Sys::Hostname;
use Carp qw(cluck longmess shortmess);
=head2 getdbh
	Title   : getdbh
 	Usage   : connect::getdbh($ip,$user,$pwd);
 	Function: Establish the connection to the database for the host $ip, the user $user with the password $pwd
 	Returns : A connection to the databse
 	Args	: The ip adress of the host, the name of the user, the password of the user (string)
	Note    : 
=cut

sub getdbh{
	my ($config,$database) = @_;
	#cluck "here";
	my $ip = $config->{ip};
	my $db_user_name = $config->{user};
	my $db_password = $config->{pw}."";
	$ip=$ENV{POLY_DB_IP} if exists $ENV{POLY_DB_IP};
	my $port = $config->{port};
	$port = 3306 unless $port;

	confess unless exists $config->{user};
	my $dbh;
	my $host =  hostname;
	my $ip2 = $config->{ip2};
	if ($ip2){
		($ip2,$ip) = ($ip,$ip2);
	}
	else {
		$ip2 =$ip;
	}
	my $dsn = "DBI:mysql::$ip;port=$port\n";
	#	my $dsn = "DBI:MariaDB::$ip;hostname=$ip;port=$port";
	eval{
		
		require "DBD/MariaDB.pm";
		my $dsn = "DBI:MariaDB::$ip;port=$port";
	 	$dbh = DBI->connect($dsn, $db_user_name, "$db_password")|| die "Database connection not made: $DBI::errstr";
	};
	if ( $@){
		my $dsn = "DBI:mysql::$ip;port=$port\n";
		$dbh = DBI->connect($dsn, $db_user_name, "$db_password")|| die "Database connection not made: $DBI::errstr";
	}
	if ( $@){
		warn "connect 3" ;
		$dsn = "DBI:mysql::$ip2;port=$port\n";
		 $dbh = DBI->connect($dsn, $db_user_name, "$db_password")|| die "Database connection not made: $DBI::errstr";
	}
	return $dbh;
}  


=head2 returnOneVal
	Title   : returnOneVal
 	Usage   : connect::returnOneVal($dbh, $sql);
 	Function: To return the mysql query results where there is only one value
 	Returns : A value
 	Args	: The mysql query
	Note    : 
=cut

sub returnOneVal{
	my ($dbh,$sql) = @_;
	my $sth = $dbh->prepare( $sql );
	$sth->execute() || confess();
	my ($id) = $sth->fetchrow_array() ;#or confess($sql);
	return $id;
} 



=head2 compressData
	Title   : compressData
 	Usage   : connect::compressData($data);
 	Function: Compress the data to store them in the database
 	Returns : The compress data
 	Args	: Tha data to compress
	Note    : 
=cut

sub compressData{
	my ($data) = @_;
	my $data2;
	 gzip \$data => \$data2
        or die "gzip failed: $GzipError\n";
    return $data2;
        
	my $filename = "cp_".time;
	my $tempfile = "/temporary/tmp/".$filename;
	open UPLOADFILE, ">".$tempfile;
	binmode UPLOADFILE;
	print UPLOADFILE $data;
	close UPLOADFILE;
	`gzip  $tempfile`;
	
	open UPLOADFILE, $tempfile.".gz";
	binmode UPLOADFILE;
	read(UPLOADFILE, $data2, -s UPLOADFILE ) ;
	unlink($tempfile.".gz");
	unlink($tempfile);
	return $data2;	
}

=head2 uncompressData
	Title   : uncompressData
 	Usage   : connect::uncompressData($data);
 	Function: Uncompress the data to retrieve them from the database
 	Returns : The uncompress data
 	Args	: The data to uncompress
	Note    : 
=cut

sub uncompressData{
	my ($data) = @_;
	my $test;	
	 gunzip \$data => \$test
         || die "gunzip failed \n";
 
      return $test;   
}

sub return_arrayref {
	my ($dbh,$sql) = @_;
	my $sth = $dbh->prepare($sql);
	 $sth->execute();
	my $res = $sth->fetchall_arrayref();
	return $res;
}

sub return_hashref {
	my ($dbh,$sql,$key) = @_;
	my $sth = $dbh->prepare($sql);
	 $sth->execute();
	my $res = $sth->fetchall_hashref($key);
	return $res;
}

1;