#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../../../GenBo/lib/obj-nodb/";
use lib "$Bin/packages";
use Getopt::Long;
use GBuffer;
use YAML qw(DumpFile);
use Cwd qw(abs_path);
use Proc::Simple;
use Term::ANSIColor qw(colored);
use Text::Table;
use POSIX qw(strftime :sys_wait_h);
use utf8;
use open ':std', ':encoding(UTF-8)';
use Digest::MD5 qw(md5_hex);
use Data::Dumper;
use IO::Prompt;
# ???????????????????????????????????????????????????????????
#  OPTIONS
# ???????????????????????????????????????????????????????????

my $project_names;
my $cclean;
my $poll_interval = 10;
my $retry;

GetOptions(
    'project=s'  => \$project_names,
    'clean=s'    => \$cclean,
    'interval=i' => \$poll_interval,
    'retry=i'    => \$retry,
);
die "Usage: $0 -project=proj1,proj2,...\n" unless $project_names;

# ???????????????????????????????????????????????????????????
#  VARIABLES GLOBALES
# ???????????????????????????????????????????????????????????

my $buffer   = GBuffer->new();
my $pathRoot = $buffer->config_path("root", "project_pipeline");

my $path;
my $fork    = 64;
my $version = "";
my @projects   = sort split(",", $project_names);
my $id = substr(md5_hex(join(",", @projects)), 0, 12);  # 12 chars suffisent
my $id_job = $id."_".time;
$ENV{PIPELINE_ID} = $id_job;
$path = $pathRoot . "/tmp.snake_" . $id . "/";
my @chrs       = (1..22, 'X', 'Y', 'MT');
my $total_chrs = scalar @chrs;

my @chr_steps = qw(
    store_ids store_annot denovo
    update_vector_score update_hash_variant_chromosome tiny_rocks
);

my @patient_steps = qw(
    transcripts_dude genes_dude
);

my @project_steps_genome = qw(
    all_patients all_chromosomes merge_polyviewer_db global_annotation_by_patient
    duck_cache_store_annotations add_project_parquet hotspot 
    quality_check cnv1 cnv2 sv1 sv2 
);
my @project_steps_exome = qw(
    all_patients all_chromosomes merge_polyviewer_db global_annotation_by_patient
    duck_cache_store_annotations add_project_parquet hotspot polydude
    project_done
);

my @project_steps;
# ???????????????????????????????????????????????????????????
#  DÉTECTION TYPE (exome/genome) + HASH PROJET
# ???????????????????????????????????????????????????????????



# ???????????????????????????????????????????????????????????
#  CONSTRUCTION DES DISPLAY STEPS
# ???????????????????????????????????????????????????????????
my ($hash_project, $type) = return_hash_project(\@projects);
my @display_steps;
push @display_steps, map { { name => $_, type => "chr" } } @chr_steps;
if ($type eq "exome"){
push @display_steps, map { { name => $_, type => "patient" } } @patient_steps ;
push @display_steps,();
push @display_steps, map { { name => $_, type => "proj" } }  @project_steps_exome;
    unshift(@project_steps_exome, "all_patients","all_chromosomes");
  @project_steps = @project_steps_exome;
}
else {
push @display_steps, map { { name => $_, type => "proj" } } @project_steps_genome;
      unshift(@project_steps_genome,"all_chromosomes");
      @project_steps = @project_steps_genome;
 }
 
 push(@project_steps,"project_done");
# ???????????????????????????????????????????????????????????
#  SETUP : RETRY OU NOUVEAU RUN
# ???????????????????????????????????????????????????????????
my $restart_from_scratch = 1;
if (-e $path){
	my $choice = prompt(colored("** Cache detected, resuming from existing cache (y/n)", "cyan"));
	 if ($choice eq "n") {
		print colored("\n  Cleanup complete - fresh start initiated", "cyan");
        system("rm -rf $path");
        system("mkdir -p $path/logs");
        print colored("\n  => Done \n\n", "blue");
	 }
	 else {
		print colored("\n   Nothing changed ? really ?  maybe this time it'll magically work.", "green ");
	$restart_from_scratch = 0;	
    clean_error_steps();
    }
}
if  ($restart_from_scratch ) {
    $path = $pathRoot . "/tmp.snake_" . $id . "/";
    system("mkdir -p $path/logs");

    clean(\@projects) if $cclean;
 	my ($hash_project, $type) = return_hash_project(\@projects);
 
    # Config YAML
    my $config = {
        projects         => $hash_project,
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




# ???????????????????????????????????????????????????????????
#  LANCEMENT SNAKEMAKE
# ???????????????????????????????????????????????????????????

my $snfile = "Snakefile";
$snfile = "Snakefile.exome" if ($type eq "exome");

my $cmd = qq{run_singularity.sh snakemake.sif snakemake};
my $snakemake_cmd = join(" ",
    "snakemake3.13",
    "--snakefile $Bin/$snfile",
    "-j 100 --keep-going",
    "--executor cluster-generic",
     "--cluster-generic-submit-cmd '$Bin/submit_job.sh'",
   # '--cluster-generic-submit-cmd "sbatch --cpus-per-task=20 --partition=defq --parsable"',
    "--cluster-generic-cancel-cmd scancel",
    "--show-failed-logs -p --quiet progress",
);

my $proc = Proc::Simple->new();
$proc->redirect_output("$path/snake_stdout.log", "$path/snake_stderr.log");
#warn "cd $path && $snakemake_cmd";
#die();
#warn  "cd $path && $snakemake_cmd";
#die();
$proc->start("bash", "-c", "cd $path && $snakemake_cmd");
sleep(2);

if (!$proc->poll()) {
    print colored("ERREUR: Snakemake n'a pas démarré\n", "red bold");
    exit 1;
}

my $snake_pid = $proc->pid();

# ???????????????????????????????????????????????????????????
#  GESTION DES SIGNAUX (Ctrl+C, kill, etc.)
# ???????????????????????????????????????????????????????????

my $interrupted = 0;

$SIG{INT} = $SIG{TERM} = sub {
    my ($sig) = @_;
    return if $interrupted;
    $interrupted = 1;

    print "\n";
    print colored("----------------------------------------------------\n", "red bold");
    print colored("+      Emergency stop detected (SIG$sig)                  +\n", "red bold");
    print colored("+  please be patient during clean up. don't touch anything !!!!           +\n", "red bold");
    print colored("----------------------------------------------------\n", "red bold");

    if ($proc->poll()) {
        print colored("  => Envoi SIGINT à Snakemake (PID $snake_pid)...\n", "yellow");
        $proc->kill();
        my $wait_max = 20;
        my $waited   = 0;
        while ($proc->poll() && $waited < $wait_max) {
            print colored("  => wait for snakemake ... (${waited}s/${wait_max}s)\n", "yellow")
                if $waited % 5 == 0;
            sleep(2);
            $waited += 2;
        }
         $proc->kill();
		 cleanup_slurm_jobs();
        if ($proc->poll()) {
            print colored("  => snake make problem, SIGKILL...\n", "red");
            kill('KILL', $snake_pid);
            sleep(1);
        }
    }

   

    scan_all();
    display_table();
    display_summary();
    save_report();

    print colored("\npipeline exit peacefully.\n", "red bold");
    print colored("to restart -retry=$id\n", "dark");
    exit 130;
};

$SIG{CHLD} = 'IGNORE';

# ???????????????????????????????????????????????????????????
#  BOUCLE PRINCIPALE
# ???????????????????????????????????????????????????????????

my $start_time = time;
my %data;

# Pré-initialiser avec les bons totaux
for my $proj (@projects) {
    for my $step (@chr_steps) {
        $data{$proj}{$step} = { ok => 0, err => 0, total => $total_chrs };
    }
    if ($type eq "exome") {
        my @patients = @{$hash_project->{$proj}{patients}};
        my $total_patients = scalar @patients;
        for my $step (@patient_steps) {
            $data{$proj}{$step} = { ok => 0, err => 0, total => $total_patients };
        }
    }
    for my $step (@project_steps) {
        next if $step =~ /all_/;
        $data{$proj}{$step} = { ok => 0, err => 0, total => 1 };
    }
}
while (1) {
    last if $interrupted;

    scan_all();
    display_table();
    display_slurm_bar();    # ? une ligne en plus

    last unless $proc->poll();
	 sleep($poll_interval);
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

#sub cleanup_slurm_jobs {
#    my $log_file = "$path/snake_stdout.log";
#    my @job_ids;
#
#    if (-e $log_file) {
#        open(my $fh, "<", $log_file) or return;
#        while (<$fh>) {
#            if (/Submitted.*?(\d{5,})/ || /^(\d{5,})$/) {
#                push @job_ids, $1;
#            }
#        }
#        close($fh);
#    }
#
#    if (@job_ids) {
#        print colored("  => Annulation de " . scalar(@job_ids) . " jobs SLURM...\n", "yellow");
#        my $ids = join(",", @job_ids);
#        system("scancel $ids 2>/dev/null");
#        print colored("  => scancel $ids\n", "dark");
#    }
#
#    my $user = $ENV{USER} || (getpwuid($<))[0];
#    my @running = `squeue -u $user -h -o '%i %j' 2>/dev/null`;
#    my @snake_jobs;
#    for my $line (@running) {
#        chomp $line;
#        if ($line =~ /^(\d+)\s+snakejob/) {
#            push @snake_jobs, $1;
#        }
#    }
#
#    if (@snake_jobs) {
#        print colored("  => " . scalar(@snake_jobs) . " jobs snakemake encore actifs, annulation...\n", "yellow");
#        system("scancel " . join(" ", @snake_jobs) . " 2>/dev/null");
#    }
#
#    print colored("  => Nettoyage SLURM terminé.\n", "green");
#}


# ???????????????????????????????????????????????????????????
#  SCAN
# ???????????????????????????????????????????????????????????

sub clean_error_steps {
    my $logs = "$path/logs";
    for my $proj (@projects) {
        my $pd = "$logs/$proj";
        next unless -d $pd;

        # Per chromosome
        for my $step (@chr_steps) {
            for my $chr (@chrs) {
                unlink "$pd/$chr/$step.err" if (-e "$pd/$chr/$step.err");
            }
        }

        # Per patient (exome uniquement)
        if ($type eq "exome") {
            my @patients = @{$hash_project->{$proj}{patients}};
            for my $step (@patient_steps) {
                for my $pat (@patients) {
                    unlink "$pd/$pat/$step.err" if (-e "$pd/$pat/$step.err");
                }
            }
        }

        # Per project
        for my $step (@project_steps) {
            next if $step eq "all_chromosomes" or $step eq "all_patients";
            unlink "$pd/$step.err" if (-e "$pd/$step.err");
        }
    }
}

sub scan_all {
    my $logs = "$path/logs";
    for my $proj (@projects) {
        my $pd = "$logs/$proj";
        next unless -d $pd;

        # Per chromosome
        for my $step (@chr_steps) {
            my ($ok, $err) = (0, 0);
            for my $chr (@chrs) {
                if    (-e "$pd/$chr/$step.err") { $err++ }
                elsif (-e "$pd/$chr/$step.ok")  { $ok++  }
            }
            $data{$proj}{$step} = { ok => $ok, err => $err, total => $total_chrs };
        }

        # Per patient (exome uniquement)
        if ($type eq "exome") {
            my @patients = @{$hash_project->{$proj}{patients}};
            my $total_patients = scalar @patients;
            for my $step (@patient_steps) {
                my ($ok, $err) = (0, 0);
                for my $pat (@patients) {
                    if    (-e "$pd/$pat/$step.err") { $err++ }
                    elsif (-e "$pd/$pat/$step.ok")  { $ok++  }
                }
                $data{$proj}{$step} = { ok => $ok, err => $err, total => $total_patients };
            }
        }

        # Per project
        for my $step (@project_steps) {
            next if $step eq "all_chromosomes" or $step eq "all_patients";
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
        print colored(" Waiting Room  $elapsed_str   $id  " .
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

    my @type_steps  = ("chr", "patient", "proj");
    my $row_idx     = 0;
    my $first_group = 1;

    foreach my $ts (@type_steps) {
        my @asteps = grep { $_->{type} eq $ts } @display_steps;

        # Pas de patient steps en mode genome (déjà filtré dans @display_steps,
        # mais double sécurité)
        next if $ts eq "patient" && $type eq "genome";
        next unless @asteps;

        # Séparateur entre groupes
        if (!$first_group) {
            $sep_before{$row_idx} = 1;
        }
        $first_group = 0;

        foreach my $as (@asteps) {
            my $step = $as->{name};
            my @row  = ($step);
            for my $proj (@projects) {
                push @row, format_cell($data{$proj}{$step});
            }
            push @body_rows, \@row;
            $row_idx++;
        }
    }

    # Ligne d'état par projet
    $sep_before{$row_idx} = 1;
    my @state_row = ("State");
    push @state_row, project_state($_) for @projects;
    push @body_rows, \@state_row;

    $table->load(@body_rows);

    my $rule_thick = $table->rule('=', '+');
    my $rule_thin  = $table->rule('-', '+');

    print $rule_thick;
    print colorize($table->title);
    print $rule_thick;

    my $n = scalar @body_rows;
    for my $i (0 .. $n - 1) {
        print $rule_thin if $sep_before{$i};
        my $line = $table->body($i, 1);
        print colorize($line) if defined $line;
    }
    print $rule_thick;
    print colored("\ Pipeline Id = $id Job_id= $id_job \n", "blue");
}


sub format_cell {
    my ($d) = @_;
    my $ICON_OK      = "OK ";
    my $ICON_ERR     = "ERR";
    my $ICON_RUN     = "RUN";
    my $ICON_PENDING = "PEN";

    return "$ICON_PENDING  0/0"                                       unless $d;
    return "$ICON_ERR ERR"                                            if $d->{err} > 0;
    return sprintf("$ICON_OK %2d/%-2d", $d->{ok}, $d->{total})       if $d->{ok} == $d->{total};
    return sprintf("$ICON_RUN %2d/%-2d", $d->{ok}, $d->{total})      if $d->{ok} > 0;
    return sprintf("$ICON_PENDING  0/%-2d", $d->{total});
}

sub colorize {
    my ($text) = @_;
    my $ICON_OK      = "OK ";
    my $ICON_ERR     = "X";
    my $ICON_RUN     = ".";
    my $ICON_PENDING = "?";

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
    for my $step (@chr_steps, @patient_steps, @project_steps) {
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
            for my $step (@chr_steps, @patient_steps, @project_steps) {
                my $d = $data{$proj}{$step};
                push @info, $step if $d && $d->{err} > 0;
            }
        }
        elsif ($state eq "INTERRUPTED") {
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
    print $st->rule('=', '+');
    print $st->title;
    print $st->rule('=', '+');
    print colorize($st->body);
    print $st->rule('=', '+');

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
    print $fh "RAPPORT " . strftime("%Y-%m-%d %H:%M:%S", localtime) . "\n";
    print $fh "INTERRUPTION MANUELLE\n" if $interrupted;
    printf $fh "Durée: %02d:%02d:%02d\n\n",
        int($elapsed/3600), int(($elapsed%3600)/60), $elapsed%60;

    for my $proj (@projects) {
        print $fh "[$proj] " . project_state($proj) . "\n";
        for my $sd (@display_steps) {
            my $step = $sd->{name};
            my $d = $data{$proj}{$step} // { ok => 0, err => 0, total => 0 };
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
#  RETURN HASH PROJECT + DETECTION TYPE
# ???????????????????????????????????????????????????????????

sub return_hash_project {
    my ($list) = @_;
    my $hash;
    my $type_count = {};

    for my $project_name (@$list) {
        my $buffer  = GBuffer->new();
        my $project = $buffer->newProject(-name => $project_name);
	#	 my $log=   = $project->getCacheDir()."/logs";
# 		system("mkdir $log && chmod a+rwx $log") unless -e $log;
        my @patient_names;
        foreach my $patient (@{$project->getPatients}) {
            push @patient_names, $patient->name;
            
            $type_count->{exome}++  if $patient->isExome();
            $type_count->{genome}++ if $patient->isGenome();
        }
        $hash->{$project_name} = {
            patients => \@patient_names,
            #log => $log,
        };
    }
    my @types = keys %$type_count;
    die("Erreur : mélange exome/genome interdit (trouvé: " . join(", ", @types) . ")\n")
        if scalar @types != 1;

    return ($hash, $types[0]);
}

sub get_slurm_stats {
    my $prefix = "pp_${id_job}";
    my %stats = (
        RUNNING   => 0,
        PENDING   => 0,
        COMPLETED => 0,
        FAILED    => 0,
        TIMEOUT   => 0,
        CANCELLED => 0,
        OTHER     => 0,
        total     => 0,
    );

    # --allocations = uniquement les jobs principaux, pas les .batch/.extern
    my @lines = `sacct --name=$prefix --noheader -P --allocations --format=State 2>/dev/null`;
    for my $line (@lines) {
        chomp $line;
        $line =~ s/\s+//g;
        $line =~ s/\s+by.*//;    # "CANCELLED by 1000" ? "CANCELLED"

        $stats{total}++;
        if (exists $stats{$line}) { $stats{$line}++ }
        else                      { $stats{OTHER}++ }
    }

    return \%stats;
}
sub cleanup_slurm_jobs {
    my $prefix = "pp_${id}";
    print colored("  => scancel --name=$prefix\n", "yellow");
    system("scancel --name=$prefix 2>/dev/null");
    print colored("  => Nettoyage SLURM terminé.\n", "green");
}
sub display_slurm_bar {
    my $stats = get_slurm_stats();
    return unless $stats->{total} > 0;

    my $line = sprintf(
        "  SLURM: %d jobs | %s %s %s %s",
        $stats->{total},
        colored(sprintf("RUN:%d",  $stats->{RUNNING}),   "green"),
        colored(sprintf("PEN:%d",  $stats->{PENDING}),   "yellow"),
        colored(sprintf("OK:%d",   $stats->{COMPLETED}), "cyan"),
        colored(sprintf("FAIL:%d", $stats->{FAILED}),    "red"),
    );

    # Barre de progression visuelle
    my $done = $stats->{COMPLETED} + $stats->{FAILED};
    my $pct  = $stats->{total} > 0 ? int(100 * $done / $stats->{total}) : 0;
    my $bar_len = 30;
    my $filled  = int($bar_len * $done / ($stats->{total} || 1));
    my $bar = "[" . "#" x $filled . "." x ($bar_len - $filled) . "]";

    print "$line  $bar $pct%\n";
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