#!/usr/bin/perl
use strict;
use warnings;
use Net::SFTP;
use File::Basename;
use Time::Local;
use POSIX qw(strftime);
use Getopt::Long;
use IO::Prompt;
use Data::Dumper;

# === OPTIONS & FLAGS ===
my $dry_run     = 1;
my $verbose     = 0;
my $date_str    = '2025-01-01';
my $remote_path;
my $show_help   = 0;

GetOptions(
    'delete'     => sub { $dry_run = 0 },
    'verbose'    => \$verbose,
    'before=s'   => \$date_str,
    'remote=s'   => \$remote_path,
    'help'       => \$show_help,
) or usage();

usage() if $show_help;

# === AIDE / USAGE ===
sub usage {
    print << "EOF";
Usage: $0 [options]

Options :
  --delete                 Exécute réellement la suppression (sinon dry-run par défaut)
  --verbose                Affiche les étapes de suppression
  --before=YYYY-MM-DD      Date limite de modification (par défaut : 2025-01-01)
  --remote=PATH            Répertoire distant à traiter (ex: /inem_asnafi) (par défaut : liste les répertoires)
  --help                   Affiche cette aide

Exemple :
  $0 --verbose
  $0 --delete --remote=/inem_asnafi --before=2024-12-31
EOF
    exit;
}

# === PARSE DATE LIMITE ===
my ($y, $m, $d) = $date_str =~ /^(\d{4})-(\d{2})-(\d{2})$/
    or die "Format de date invalide (attendu : YYYY-MM-DD)\n";

my $limit_epoch = timelocal(0, 0, 0, $d, $m - 1, $y);

# === CONFIGURATION CONNEXION ===
my $host        = '192.168.2.111';
my $user        = 'bioinfo';
my $key_path    = '/home/masson/.ssh/known_hosts';    # <-- adapte si nécessaire

# Vérifier utilisateur (optionnel)
my $current_user = getpwuid($<);
die("You must be logged in as user 'masson' before running this script. Current user: '$current_user'\n")
    if ( $current_user ne "masson" );

# === CONNEXION SFTP ===
print "Connexion à $user\@$host...\n" if $verbose;

my $sftp = Net::SFTP->new(
    $host,
    user      => $user,
    key_path  => $key_path,
    more      => ['-oBatchMode=yes'],
);

# === DÉTERMINER LE RÉPERTOIRE DISTANT ===
$remote_path =~ s{/$}{};                       # enlever slash terminal
my @ls = $sftp->ls('.', wanted => sub { my $n = shift; return $n->{longname} =~ /^d/; });
my @sorted_dir = sort { $a->{filename} cmp $b->{filename} } @ls;
@sorted_dir = grep { $_->{filename} !~ /^\./ } @sorted_dir;
if (grep {$remote_path} @sorted_dir) {
    $remote_path .= '/filetransfer' unless ($remote_path =~ /filetransfer\/?$/);
} else {
    # Sinon on affiche un menu pour sélectionner un dossier
    my @sorted_dir = sort { $a->{filename} cmp $b->{filename} } @ls;
    @sorted_dir = grep { $_->{filename} !~ /^\./ } @sorted_dir;

    my @dir_names = map { $_->{filename} } @sorted_dir;
    $remote_path = prompt("Select a directory to put the files: ", -menu => \@dir_names);

    $remote_path .= "/filetransfer" unless ($remote_path =~ /filetransfer\/?$/);
    print "\n";
}
$remote_path =~ s{/$}{};

my $total_deleted = 0;   # cumul des tailles supprimées

# === LISTER LES ENTRÉES (fichiers + dossiers) DANS $remote_path ===
my @entries = $sftp->ls($remote_path)
    or die "Impossible de lister $remote_path : " . join(' ', $sftp->status) . "\n";

foreach my $entry (@entries) {
    next unless defined $entry && defined $entry->{filename};
    my $name = $entry->{filename};
    next if $name eq '.' || $name eq '..';

    my $attrs = $entry->{a} || {};
    my $mtime = $attrs->{mtime} || 0;
    my $date_fmt = strftime('%Y-%m-%d', localtime($mtime));
    my $full_path = "$remote_path/$name";

    # Détecter répertoire via bits de permission si possible, sinon fallback au longname
    my $is_dir = 0;
    if (defined $attrs->{perm}) {
        $is_dir = ($attrs->{perm} & 0040000) ? 1 : 0;
    } else {
        $is_dir = ($entry->{longname} && $entry->{longname} =~ /^d/) ? 1 : 0;
    }

    if ($is_dir) {
        # dossier : suppression récursive si date antérieure au seuil
        if ($mtime <= $limit_epoch) {
		    my $size_h = human_size(dir_size($full_path));    # Taille calculée du dossier avant suppression
            print ">> [$date_fmt] Dossier éligible : $full_path ($size_h)\n";
            delete_remote_dir($full_path);
        } else {
		    my $size_h = human_size(dir_size($full_path));	# Taille calculée du dossier avant suppression
            print ">> [$date_fmt] Dossier ignoré (trop récent) : $full_path ($size_h)\n" if $verbose;
        }
    } else {
        # fichier directement sous $remote_path : traiter aussi
        if ($mtime <= $limit_epoch) {
		    my $bytes = dir_size($full_path,1);
		    my $size_h = human_size($bytes);	# Taille calculée du dossier avant suppression
            print ">> [$date_fmt] Fichier éligible : $full_path ($size_h)\n";
            $total_deleted += $bytes;
            if ($dry_run) {
                print "   (dry-run) suppression simulée : $full_path\n" if $verbose;
            } else {
                my $exit = $sftp->do_remove($full_path);
                if ($exit) {
                    warn "Échec suppression fichier '$full_path': " . join(' ', $sftp->status) . "\n";
                } else {
                    print "   -> Fichier supprimé : $full_path\n" if $verbose;
                }
            }
        } else {
            print ">> [$date_fmt] Fichier ignoré (trop récent) : $full_path\n" if $verbose;
        }
    }
}

# === CALCUL TAILLE D'UN DOSSIER (récursif) ===
sub dir_size {
    my ($path, $is_file) = @_;
    my $total = 0;
	
	if ($is_file) {
		my @path = split('/', $path);
		my $file = pop(@path);
		my $previous_dir = join('/',@path);
		my @items = $sftp->ls($previous_dir) || ();
		die() if (scalar @items > 1);
		my @file = grep{$_->{filename} eq $file} @{$items[0]};
		my $attrs = $file[0]->{a} || {};
		$total += $attrs->{size} || 0;
	}
	else {
	    my @items = $sftp->ls($path) || ();
	    die() if (scalar @items > 1);
	    foreach my $item (@{$items[0]}) {
	        my $name = $item->{filename};
	        next if !$name || $name eq '.' || $name eq '..';
	
	        my $full = "$path/$name";
	        my $attrs = $item->{a} || {};
	        my $is_dir = (defined $attrs->{perm} && ($attrs->{perm} & 0040000)) ? 1 :
	                     ($item->{longname} && $item->{longname} =~ /^d/) ? 1 : 0;
	
	        if ($is_dir) {
	            $total += dir_size($full);
	        } else {
	            $total += $attrs->{size} || 0;
	        }
	}
   }
    return $total;
}

# === SUPPRESSION RÉCURSIVE ===
sub delete_remote_dir {
    my ($path) = @_;

    # Taille calculée du dossier avant suppression
    my $bytes  = dir_size($path);
    my $size_h = human_size($bytes);

#    print ">> [$date_str] Dossier éligible : $path ($size_h)\n";

    $total_deleted += $bytes;# unless $dry_run;

    my @items = $sftp->ls($path) || ();
    die() if (scalar @items > 1);
    foreach my $item (@{$items[0]}) {
        my $name = $item->{filename};
        next if !$name || $name eq '.' || $name eq '..';

        my $full = "$path/$name";
        my $attrs = $item->{a} || {};
        my $is_dir = (defined $attrs->{perm} && ($attrs->{perm} & 0040000)) ? 1 :
                     ($item->{longname} && $item->{longname} =~ /^d/) ? 1 : 0;

        if ($is_dir) {
            delete_remote_dir($full);
        } else {
            print "   -> Suppression fichier : $full\n" if $verbose;
            unless ($dry_run) {
                my $exit = $sftp->do_remove($full);
                warn "Échec suppression fichier '$full': ".join(' ', $sftp->status)."\n" if ($exit);
            }
        }
    }

    print "   -> Suppression dossier : $path\n" if $verbose;
    unless ($dry_run) {
        my $exit = $sftp->do_rmdir($path);
        warn "Échec suppression dossier '$path': ".join(' ', $sftp->status)."\n" if ($exit);
    }
}

# === FORMATAGE TAILLE LISIBLE ===
sub human_size {
    my ($bytes) = @_;
    return "0 B" unless $bytes;
    my @units = ('B','KB','MB','GB','TB');
    my $i = 0;
    while ($bytes > 1024 && $i < @units-1) {
        $bytes /= 1024;
        $i++;
    }
    return sprintf("%.2f %s", $bytes, $units[$i]);
}

# === FIN ===
print $dry_run
    ? sprintf("\n[MODE TEST] Aucun fichier n'a été supprimé. Eligible à supprimer : %s\n", human_size($total_deleted))
    : sprintf("\n[SUPPRESSION ACTIVE] Total supprimé : %s\n", human_size($total_deleted));
