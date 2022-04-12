#!/usr/bin/perl
use strict;



my ($cmd,$opt) = split(/\?/,$ARGV[0]);

 $opt =~ s/&/ /g;
 $cmd =~s/\/\//\//g;	

my $cmd2 = "/poly-disk/www/";
my @pcmds = split("/",$cmd);
shift(@pcmds);
shift(@pcmds);


my $cmd2 = "/software/polyweb/poly-disk/www/".join("/",@pcmds)." ".$opt;
warn "$cmd2 \n";

warn "-*-*-*-*-*-*-*-*-*-*\n";
system($cmd2);


warn "-*-*-*-*-*-*-*-*-*-*\n";
warn "$cmd2 \n";
