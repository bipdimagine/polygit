#!/usr/bin/perl
# permet de renvoyer petit a petit les print et non pas de tout mettre en buffer et tout sortir a la fin du script
$|=1;
use CGI qw/:standard :html3/;

use strict;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../GenBo";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../GenBo/lib/obj-nodb";
use lib "$Bin/../packages/export";
use lib "$Bin/../packages/layout";
use lib "$Bin/../packages/validation_variation"; 
use lib "$Bin/../cache_nodb/scripts/";
use lib "$Bin/../GenBo/lib/obj-nodb/packages";
use lib "$Bin/../GenBo/lib/obj-nodb/packages/cache/polydiag/";


use GBuffer;
use export_data;
use JSON;
use polyweb_dude;
use VcfMerge;
use GenBoNoSql;
use Set::IntervalTree;
use Spreadsheet::WriteExcel;
use Compress::Snappy;
use Storable qw(store retrieve freeze dclone thaw);
use POSIX qw(strftime);
use CGI::Session;
use html; 
use Carp;
use Cache_Commons;
use QueryVectorFilter;
use IO::Handle;
use update_variant_editor;
use Getopt::Long;
use File::Basename;


my $io = IO::Handle->new();
$io->autoflush(1);

my $fork;
my $release;
my $only_users;
GetOptions(
	'fork=s' => \$fork,
);


$fork = 1 unless ($fork);
my $buffer = new GBuffer;

my @lCmds;
my $dirname = dirname(__FILE__);


my $forkA = int($fork / 2) if $fork > 1;
my $forkB = 1;
$forkB = 2 if $fork > 1;

my $i = 0;
my $pm = new Parallel::ForkManager($forkA);
my $h_login = get_hash_login_pwd($buffer);
foreach my $user (sort keys %$h_login) {
	my $pid = $pm->start and next;
	my $pwd = $h_login->{$user}->{PW};
	my $cmd1 = "$dirname/user_polybtf.pl login=$user pwd=$pwd";
	`$cmd1`;
	warn "-> login: $user done.\n";
	$pm->finish();
}
$pm->wait_all_children();

#my $pm2 = new Parallel::ForkManager($fork);
#foreach my $user (sort keys %$h_login) {
#	my $pid = $pm2->start and next;
#	my $pwd = $h_login->{$user}->{PW};
#	my $cmd2 = "$dirname/check_new_hgmd_clinvar.pl user_name=$user pwd=$pwd print=1 view_others=1 fork=1";
#	$cmd2 .= " force_db_annot=$release" if $release;
#	`$cmd2`;
#	warn "-> login: $user done.\n";
#	$pm2->finish();
#}
#$pm2->wait_all_children();



sub get_hash_login_pwd {
	my ($buffer) = @_;
	my $sql = "SELECT LOGIN, PW FROM bipd_users.USER;";
	my $sth = $buffer->dbh->prepare($sql);
	$sth->execute();
	my $h = $sth->fetchall_hashref("LOGIN");
	return $h;
}
