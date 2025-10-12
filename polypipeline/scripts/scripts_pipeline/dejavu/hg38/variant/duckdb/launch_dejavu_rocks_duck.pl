#!/usr/bin/perl
use CGI qw/:standard :html3/;
use FindBin qw($Bin); 
use strict;
use Getopt::Long;
use List::Util qw(shuffle); 
use Email::Send;
use Email::Send::Gmail;
use Email::Simple::Creator;
 

 my @chromosomes = (1..22,'X','Y','MT');
  my $dir_tmp = "/data-beegfs/tmp/";
 my $filetmp =  "/data-beegfs/tmp/cmd.".time;


 open(CMD,">$filetmp");
 foreach my $chr (@chromosomes){
 	print CMD "$Bin/dejavu_rocks.pl -chr=$chr -fork=5 -version=HG38\n"; 
 	print CMD "$Bin/dejavu_rocks.pl -chr=$chr -fork=5 -version=HG19\n"; 
 	
 }
 
 system("cd /data-beegfs/tmp/ ; cat $filetmp | /software/bin/run_cluster.pl -cpu=20 -limit=10 -noprint=1 && touch $filetmp.ok");
 
 if (-e "$filetmp.ok"){
	exit(0);
 	system("$Bin/check_dejavu.pl  && touch $filetmp.check ");
 	if (-e "$filetmp.check"){
 		unlink $filetmp;
 		unlink  "$filetmp.ok";
 		unlink  "$filetmp.check";
 		my $d = scalar localtime time;
 		error ("DEJAVU :::  OK ".$d);
 		exit(0);
 	}
 	else {
 		error("DEJAVU ::: problem check dejavu ");
 	}
 }
 
 error ("dejavu_lite_memory => $dir_tmp ");
 
 
 
 sub error {
	my ($message) = @_;
exit(1);
   my $email = Email::Simple->create(
      header => [
          From    => 'bipd.imagine@gmail.com',
          To      => 'pnitschke@gmail.com',
          Subject => $message,
      ],
      body => "Only one  thing to say : ===> * ".$message."*",
  );
  my $sender = Email::Send->new(
      {   mailer      => 'Gmail',
          mailer_args => [
              username => 'backup.bipd@gmail.com',
              password => 'backupnecobi',
          ]
      }
  );
  eval { $sender->send($email);  };
  die "Error sending email: $@" if $@;                                               
 
    exit(0);
}