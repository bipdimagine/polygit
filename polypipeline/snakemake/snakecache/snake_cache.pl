#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/packages";
use Getopt::Long;
use GBuffer;
use YAML::XS qw(DumpFile);
use Cwd qw(abs_path);
use Proc::Simple;
use Term::ANSIColor qw(colored);
use Text::Table;
use POSIX qw(strftime :sys_wait_h);
use utf8;
use open ':std', ':encoding(UTF-8)';

my $project_names;
my $cclean;
my $poll_interval = 10;
my $retry;
GetOptions(
    'project=s'  => \$project_names,
    'clean=s'    => \$cclean,
    'interval=i' => \$poll_interval,
    'retry=i' => \$retry,
);
die "Usage: $0 -project=proj1,proj2,...\n" unless $project_names;

my $buffer   = GBuffer->new();
my $pathRoot = $buffer->config_path("root", "project_pipeline");
my $id       =  time;
my $path;
my $fork    = 20;
my $version = "";
my @projects   = split(",", $project_names);
my @chrs       = (1..22, 'X', 'Y', 'MT');
my $total_chrs = scalar @chrs;
my @chr_steps = qw(
    store_ids store_annot denovo
    update_vector_score update_hash_variant_chromosome tiny_rocks
);

my @project_steps = qw(
    all_chromosomes merge_polyviewer_db global_annotation_by_patient
    duck_cache_store_annotations add_project_parquet hotspot
    quality_check cnv1 cnv2 sv1 sv2 project_done
);

my @display_steps;
push @display_steps, map { { name => $_, type => "chr" } } @chr_steps;
push @display_steps, map { { name => $_, type => "proj" } }
    grep { $_ ne "all_chromosomes" } @project_steps;

if ($retry){
	$path     = $pathRoot . "/tmp.snake_" . $retry . "/";
	die("no dir $path") unless -e $path;
	clean_error_steps();
}
else {
 $path     = $pathRoot . "/tmp.snake_" . $id . "/";
system("mkdir -p $path/log");



clean(\@projects) if $cclean;

#  Config YAML 
my $config = {
    projects         => \@projects,
    fork             => $fork,
    version          => $version,
    chromosomes      => \@chrs,
    base_dir         => abs_path("$Bin/../../../"),
    rocks_scripts    => abs_path("$Bin/../../../polymorphism-cgi/cache_nodb/scripts/rocks"),
    pipeline_scripts => abs_path("$Bin/../../../polypipeline/scripts/scripts_pipeline"),
    perl             => "/usr/bin/perl",
};
DumpFile("$path/config.yaml", $config);
}
# ??????????????????????????????????????????????????????????
#  LANCEMENT SNAKEMAKE
# ??????????????????????????????????????????????????????????

my $snakemake_cmd = join(" ",
    "snakemake3.13",
    "--snakefile $Bin/Snakefile",
    "-j 100 --keep-going",
    "--executor cluster-generic",
    '--cluster-generic-submit-cmd "sbatch --cpus-per-task=20 --partition=defq --parsable"',
    "--cluster-generic-cancel-cmd scancel",
    "--show-failed-logs -p --quiet progress",
);

my $proc = Proc::Simple->new();
$proc->redirect_output("$path/snake_stdout.log", "$path/snake_stderr.log");
$proc->start("bash", "-c", "cd $path && $snakemake_cmd");
sleep(2);

if (!$proc->poll()) {
    print colored("ERREUR: Snakemake n'a pas démarré\n", "red bold");
    exit 1;
}

my $snake_pid = $proc->pid();

# ??????????????????????????????????????????????????????????
#  GESTION DES SIGNAUX (Ctrl+C, kill, etc.)
# ??????????????????????????????????????????????????????????
#
#  Sans ça : Ctrl+C tue le Perl, mais les jobs SLURM
#  continuent à tourner en arrière-plan sur le cluster !
#
#  Avec ça :
#    1. On intercepte SIGINT / SIGTERM
#    2. On envoie SIGINT à snakemake
#    3. Snakemake utilise --cluster-generic-cancel-cmd scancel
#       pour annuler tous les jobs SLURM qu'il a soumis
#    4. On attend que snakemake finisse son ménage
#    5. On fait un scancel de sécurité au cas où
#    6. On affiche le dernier état + rapport
#    7. On quitte proprement

my $interrupted = 0;

$SIG{INT} = $SIG{TERM} = sub {
    my ($sig) = @_;
    return if $interrupted;    # éviter double-entrée
    $interrupted = 1;

    print "\n";
    print colored("----------------------------------------------------\n", "red bold");
    print colored("+      Emergency stop detected (SIG$sig)                  +\n", "red bold");
    print colored("+  please be patient during clean up. don't touch anything !!!!           +\n", "red bold");
    print colored("----------------------------------------------------\n", "red bold");

    # ?? Étape 1 : Envoyer SIGINT à snakemake ??
    # Snakemake intercepte SIGINT et appelle scancel sur chaque job
    if ($proc->poll()) {
        print colored("  => Envoi SIGINT à Snakemake (PID $snake_pid)...\n", "yellow");
        kill('INT', $snake_pid);

        # ?? Étape 2 : Attendre que snakemake fasse le ménage ??
        # (il appelle scancel pour chaque job soumis)
        my $wait_max = 60;    # secondes max d'attente
        my $waited   = 0;
        while ($proc->poll() && $waited < $wait_max) {
            print colored("  => wait for snakemake ... (${waited}s/${wait_max}s)\n", "yellow")
                if $waited % 10 == 0;
            sleep(2);
            $waited += 2;
        }

        # ?? Étape 3 : Si snakemake ne répond plus, forcer ??
        if ($proc->poll()) {
            print colored("  => snake make problem,  SIGKILL...\n", "red");
            kill('KILL', $snake_pid);
            sleep(1);
        }
    }

    # ?? Étape 4 : Scancel de sécurité ??
    # Chercher les jobs SLURM orphelins de cet utilisateur
    cleanup_slurm_jobs();

    # ?? Étape 5 : Dernier état + rapport ??
    scan_all();
    display_table();
    display_summary();
    save_report();

    print colored("\npipeline exit peacefully .\n", "red bold");
    print colored("to restart -retry=$id\n", "dark");
    exit 130;
};

# ?? Aussi intercepter la mort du process fils ??
$SIG{CHLD} = 'IGNORE';



my $start_time = time;
my %data;

while (1) {
    last if $interrupted;

    scan_all();
    display_table();

    last unless $proc->poll();

    # Sleep interruptible (pour réagir vite au Ctrl+C)
    for (1 .. $poll_interval) {
        last if $interrupted;
        sleep(1);
    }
}

unless ($interrupted) {
    scan_all();
    display_table();
    display_summary();
    save_report();
}

my $has_fail = grep { project_failed($_) } @projects;
exit($has_fail ? 1 : 0);


# ???????????????????????????????????????????????????????????
#  NETTOYAGE JOBS SLURM ORPHELINS
# ???????????????????????????????????????????????????????????

sub cleanup_slurm_jobs {
    # Récupérer les job IDs depuis le log snakemake
    my $log_file = "$path/snake_stdout.log";
    my @job_ids;

    if (-e $log_file) {
        open(my $fh, "<", $log_file) or return;
        while (<$fh>) {
            # Snakemake avec --parsable affiche les job IDs
            if (/Submitted.*?(\d{5,})/ || /^(\d{5,})$/) {
                push @job_ids, $1;
            }
        }
        close($fh);
    }

    if (@job_ids) {
        print colored("  ? Annulation de " . scalar(@job_ids) . " jobs SLURM...\n", "yellow");
        my $ids = join(",", @job_ids);
        system("scancel $ids 2>/dev/null");
        print colored("  ? scancel $ids\n", "dark");
    }

    # Aussi scanner par nom de job au cas où
    my $user = $ENV{USER} || (getpwuid($<))[0];
    my @running = `squeue -u $user -h -o '%i %j' 2>/dev/null`;
    my @snake_jobs;
    for my $line (@running) {
        chomp $line;
        if ($line =~ /^(\d+)\s+snakejob/) {
            push @snake_jobs, $1;
        }
    }

    if (@snake_jobs) {
        print colored("  ? " . scalar(@snake_jobs) . " jobs snakemake encore actifs, annulation...\n", "yellow");
        system("scancel " . join(" ", @snake_jobs) . " 2>/dev/null");
    }

    print colored("  ? Nettoyage SLURM terminé.\n", "green");
}


# ???????????????????????????????????????????????????????????
#  SCAN
# ???????????????????????????????????????????????????????????
sub clean_error_steps {
	  my $logs = "$path/logs";
	   for my $proj (@projects) {
        my $pd = "$logs/$proj";
        next unless -d $pd;

        for my $step (@chr_steps) {
            my ($ok, $err) = (0, 0);
            for my $chr (@chrs) {
                unlink "$pd/$chr/$step.err" if    (-e "$pd/$chr/$step.err");
            }
        }

        for my $step (@project_steps) {
            next if $step eq "all_chromosomes";
            my ($ok, $err) = (0, 0);
            unlink "$pd/$step.err"  if    (-e "$pd/$step.err");
        }
    }
}
sub scan_all {
    my $logs = "$path/logs";
    for my $proj (@projects) {
        my $pd = "$logs/$proj";
        next unless -d $pd;

        for my $step (@chr_steps) {
            my ($ok, $err) = (0, 0);
            for my $chr (@chrs) {
                if    (-e "$pd/$chr/$step.err") { $err++ }
                elsif (-e "$pd/$chr/$step.ok")  { $ok++  }
            }
            $data{$proj}{$step} = { ok => $ok, err => $err, total => $total_chrs };
        }

        for my $step (@project_steps) {
            next if $step eq "all_chromosomes";
            my ($ok, $err) = (0, 0);
            if    (-e "$pd/$step.err") { $err = 1 }
            elsif (-e "$pd/$step.ok")  { $ok  = 1 }
            $data{$proj}{$step} = { ok => $ok, err => $err, total => 1 };
        }
    }
}


# ???????????????????????????????????????????????????????????
#  AFFICHAGE
# ???????????????????????????????????????????????????????????
sub display_table {
    print "\033[2J\033[H";

    my $elapsed = time - $start_time;
    my $elapsed_str = sprintf("%02d:%02d:%02d",
        int($elapsed/3600), int(($elapsed%3600)/60), $elapsed%60);

    if ($interrupted) {
        print colored(" ARGGGGG !!!  $elapsed_str\n\n", "red bold");
    } else {
        print colored(" Waiting Room  $elapsed_str   $id" .
            strftime("%H:%M:%S", localtime) . "\n\n", "cyan bold");
    }

    my @table_def;
    push @table_def, \" | ";
    push @table_def, " Step ";
    for my $p (@projects) {
        push @table_def, \" | ";
        push @table_def, " $p ";
    }
    push @table_def, \" | ";

    my $table = Text::Table->new(@table_def);

    my @body_rows;
    my %sep_before;
    my $row_idx   = 0;
    my $prev_type = "";

    for my $sd (@display_steps) {
        my $step = $sd->{name};
        my $type = $sd->{type};

        if ($type ne $prev_type && $prev_type ne "") {
            $sep_before{$row_idx} = 1;
        }
        $prev_type = $type;

        my @row = ($step);
        for my $proj (@projects) {
            push @row, format_cell($data{$proj}{$step});
        }
        push @body_rows, \@row;
        $row_idx++;
    }

    $sep_before{$row_idx} = 1;
    my @state_row = ("State");
    push @state_row, project_state($_) for @projects;
    push @body_rows, \@state_row;

    $table->load(@body_rows);

    # ?? Règles ASCII simples (pas d'UTF-8 box-drawing) ??
    my $rule_thick = $table->rule('=', '+');
    my $rule_thin  = $table->rule('-', '+');

    print $rule_thick;
    print colorize($table->title);    # title retourne toujours quelque chose
    print $rule_thick;

    my $n = scalar @body_rows;
    for my $i (0 .. $n - 1) {
        print $rule_thin if $sep_before{$i};

        # ?? FIX : body($start, $nb_lignes) ??
        # body($i, $i) retourne $i lignes ? 0 lignes quand $i=0 ? undef !
        # body($i, 1)  retourne 1 ligne à partir de la position $i
        my $line = $table->body($i, 1);
        print colorize($line) if defined $line;
    }
    print $rule_thick;
	print colored("\nImportant Pipeline Id = $id \n", "dark")
}




sub format_cell {
    my ($d) = @_;
    my $ICON_OK      = "OK ";  # ?
my $ICON_ERR     = "ERR";  # ?
my $ICON_RUN     = "RUN";  # ?
my $ICON_PENDING = "PEN";  # ?
    return "$ICON_PENDING  0/0"                                       unless $d;
    return "$ICON_ERR ERR"                                            if $d->{err} > 0;
    return sprintf("$ICON_OK %2d/%-2d", $d->{ok}, $d->{total})       if $d->{ok} == $d->{total};
    return sprintf("$ICON_RUN %2d/%-2d", $d->{ok}, $d->{total})      if $d->{ok} > 0;
    return sprintf("$ICON_PENDING  0/%-2d", $d->{total});
}

sub colorize {
    my ($text) = @_;
my $ICON_OK      = "OK ";  # ?
my $ICON_ERR     = "X";  # ?
my $ICON_RUN     = ".";  # ?
my $ICON_PENDING = "?";  # ?
    # Utiliser quotemeta ou \Q...\E pour échapper les caractères Unicode
    $text =~ s/(\Q$ICON_OK\E *\d+\/\d+)/colored($1, "green")/ge;
    $text =~ s/(\Q$ICON_ERR\E *ERR)/colored($1, "red bold")/ge;
    $text =~ s/(\Q$ICON_RUN\E *\d+\/\d+)/colored($1, "yellow")/ge;
    $text =~ s/(\Q$ICON_PENDING\E *\d+\/\d+)/colored($1, "dark")/ge;
    $text =~ s/\b(DONE)\b/colored($1, "green bold")/ge;
    $text =~ s/\b(FAILED)\b/colored($1, "red bold")/ge;
    $text =~ s/\b(RUNNING)\b/colored($1, "yellow bold")/ge;
    $text =~ s/\b(INTERRUPTED)\b/colored($1, "red")/ge;
    return $text;
}

sub project_state {
    my ($proj) = @_;
    my $d = $data{$proj}{"project_done"};
    return "DONE" if $d && $d->{ok} == 1;
    return "INTERRUPTED" if $interrupted;
    for my $step (@chr_steps, @project_steps) {
        my $dd = $data{$proj}{$step};
        return "FAILED" if $dd && $dd->{err} > 0;
    }
    return "RUNNING";
}

sub project_failed {
    my ($proj) = @_;
    my $s = project_state($proj);
    return $s eq "FAILED" || $s eq "INTERRUPTED";
}


# ???????????????????????????????????????????????????????????
#  RÉSUMÉ FINAL
# ???????????????????????????????????????????????????????????

sub display_summary {
    my $elapsed = time - $start_time;
    print "\n";

    my @cols_def = (\" | ", " Projet ", \" | ", " État ", \" | ", " Détail ", \" | ");
    my $st = Text::Table->new(@cols_def);

    for my $proj (@projects) {
        my $state = project_state($proj);
        my @info;
        if ($state eq "FAILED") {
            for my $step (@chr_steps, @project_steps) {
                my $d = $data{$proj}{$step};
                push @info, $step if $d && $d->{err} > 0;
            }
        } elsif ($state eq "INTERRUPTED") {
            # Trouver la dernière étape complétée
            my $last_done = "-";
            for my $sd (@display_steps) {
                my $d = $data{$proj}{$sd->{name}};
                if ($d && $d->{ok} == $d->{total} && $d->{total} > 0) {
                    $last_done = $sd->{name};
                }
            }
            push @info, "dernière ok: $last_done";
        }
        $st->add($proj, $state, join(", ", @info) || "-");
    }

    print colored("  RÉSUMÉ FINAL\n", "cyan bold");
    print $st->rule('?', '?');
    print $st->title;
    print $st->rule('?', '?');
    print colorize($st->body);
    print $st->rule('?', '?');

    printf("\n  Durée: %02d:%02d:%02d\n",
        int($elapsed/3600), int(($elapsed%3600)/60), $elapsed%60);
}


# ???????????????????????????????????????????????????????????
#  RAPPORT FICHIER
# ???????????????????????????????????????????????????????????

sub save_report {
    my $f = "$path/pipeline_report.txt";
    open(my $fh, ">", $f) or return;

    my $elapsed = time - $start_time;
    print $fh "RAPPORT ? " . strftime("%Y-%m-%d %H:%M:%S", localtime) . "\n";
    print $fh "INTERRUPTION MANUELLE\n" if $interrupted;
    printf $fh "Durée: %02d:%02d:%02d\n\n",
        int($elapsed/3600), int(($elapsed%3600)/60), $elapsed%60;

    for my $proj (@projects) {
        print $fh "[$proj] " . project_state($proj) . "\n";
        for my $sd (@display_steps) {
            my $step = $sd->{name};
            my $d = $data{$proj}{$step} // { ok=>0, err=>0, total=>0 };
            my $status = $d->{err} > 0           ? "FAILED"
                       : $d->{ok} == $d->{total} ? "DONE"
                       : $d->{ok} > 0            ? "PARTIAL $d->{ok}/$d->{total}"
                       :                           "PENDING";
            printf $fh "  %-40s %s\n", $step, $status;
        }
        print $fh "\n";
    }
    close($fh);
    print colored("Rapport: $f\n", "dark");
}


# ???????????????????????????????????????????????????????????
#  CLEAN
# ???????????????????????????????????????????????????????????

sub clean {
    my ($list) = @_;
    for my $project_name (@$list) {
        my $buffer  = GBuffer->new();
        my $project = $buffer->newProject(-name => $project_name);
        my $cache   = $project->getCacheDir();
        warn "Nettoyage $cache";
        system("rm -rf $cache/*") if -e $cache;
        $cache = $project->getCacheDir();
        my $pf = $project->parquet_cache_variants;
        system("rm $pf") if -e $pf;
        my $tr  = $project->rocks_cache_2_root_dir() . "/" . $project->name;
        my $tr2 = $project->tiny_rocks_cache_dir();
        system("rm -rf $tr/*") if -e $tr2;
        system("rmdir $cache");
        system("rmdir $tr");
    }
}