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
 
 
 my $dir_tmp = "/data-isilon/tmp/dejavu/";
 my $filetmp =  "/data-isilon/tmp/dejavu/cmd.".time;
 open(CMD,">$filetmp");
 
 foreach my $chr (@chromosomes){
 	print CMD "$Bin//dejavu_lite_memory.pl -chr=$chr -fork=10 \n"; 
 	
 }
 
 system("cat $filetmp | run_cluster.pl -cpu=10 && touch $filetmp.ok");
 
 if (-e "$filetmp.ok"){
 	system("$Bin/check_dejavu.pl  && touch $filetmp.check ");
 	if (-e "$filetmp.check"){
 		unlink $filetmp;
 		unlink  "$filetmp.ok";
 		unlink  "$filetmp.check";
 		error ("dejavu OK");
 		exit(0);
 	}
 	else {
 		error("problem check dejavu ");
 	}
 }
 
 error ("dejavu_lite_memory => $dir_tmp ");
 
 
 
 sub error {
	my ($message) = @_;

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