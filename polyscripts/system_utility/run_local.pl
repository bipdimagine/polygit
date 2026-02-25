#!/usr/bin/perl
use strict;
use warnings;
use Time::HiRes qw(time);
use Cwd;

my $file = shift or die "Usage: $0 <fichier_commandes>\n";

die "Fichier non trouvé: $file\n" unless -f $file;

open my $fh, '<', $file or die "Erreur lecture: $!";
my @commands = grep { /\S/ && !/^#/ } <$fh>;
close $fh;

chomp @commands;

my $total = scalar @commands;
die "Aucune commande trouvée\n" unless $total > 0;

my $success = 0;
my $error = 0;
my $start_time = time();

open my $log_fh, '>', '/tmp/jobs_errors.log' or die "Erreur log: $!";



print "=" x 60 . "\n";
print "? Jobs : $total\n";
print "=" x 60 . "\n\n";
my @jobe;
my @joberror;
my $nb =0;
for my $i (0..$#commands) {
    my $cmd = $commands[$i];
    my $count = $i + 1;
    my $percent = int($count * 100 / $total);
    my $filled = int($percent / 5);
    my $bar = "[" . ("=" x $filled) . (" " x (20-$filled)) . "]";
    
    my $cmd_display = substr($cmd, 0, 55);
    $cmd_display .= "..." if length($cmd) > 55;
    
    printf "%s %3d%% (%2d/%d) %s", $bar, $percent, $count, $total, $cmd_display;
    
    my $cmd_start = time();
    my $ret = system($cmd . " > /tmp/job_$count.log 2>&1");
    my $cmd_time = time() - $cmd_start;
    
    if ($ret == 0) {
        $success++;
        printf " ? (%.1fs)\n", $cmd_time;
    } else {
		push(@jobe,$i);
		push(@joberror,$cmd);
        $error++;
        printf " ?\n";
        print $log_fh "JOB $count: ERREUR\nCmd: $cmd\n\n";
    }
}

close $log_fh;

my $total_time = time() - $start_time;

print "\n" . "=" x 60 . "\n";
print "? OK : $success/$total\n";
print "? ERROR: $error/$total\n";
print "??  TIME : " . sprintf("%.1f", $total_time) . "s\n";
print "=" x 60 . "\n";
print "=" x 60 . "\n" if @joberror;
foreach my $c (@joberror){
	print "$c.\n";
}
exit($error > 0 ? 1 : 0);