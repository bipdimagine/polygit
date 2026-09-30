package slurm;

use Moo;

use Data::Dumper;
use Config::Std;
use FindBin qw($Bin);
use Storable qw(store retrieve freeze thaw);
use List::Util qw( shuffle sum max min);
use Text::Table;
use Term::Table;
use Term::ANSIColor qw(color colored);
use Carp qw(confess);
use IO::Prompt;
has jobs => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		return [];
	},
);

has submited_jobs => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		return {};
	},
);

has jobs_uniq_id => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		return {};
	},
);



has names => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		return {};
	},
);


has types => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		return {};
	},
);


has current_id => (
	is		=> 'rw',
	lazy	=> 1,
	default	=> sub {
		return time;
	},
);
has uniq_id => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		return {};
	},
);


sub return_uniq_id {
	my ($self,$job) = @_;
	my $st = $job->{name}."!".$job->{type};
	
	if (exists $self->uniq_id->{$st}){
		return  $self->uniq_id->{$st};
	}
	my $v = $self->current_id + 1 ;
	$self->current_id($v);
	return $self->current_id;
}



sub add_job {
	my ($self,$job) = @_;
	my $uid = $self->return_uniq_id($job);
	$job->{uid} = $uid;
	$self->jobs_uniq_id->{$uid} = $job;
	push(@{$self->types->{$job->{type}}},$uid);
	push(@{$self->names->{$job->{name}}},$uid);
	push(@{$self->jobs},$uid);
	return $uid;
	
}


sub get_job {
	my ($self,$uid) = @_;
	confess("no uid") unless $uid;
	die("no job with $uid") unless exists $self->{jobs_uniq_id}->{$uid};
	
	return $self->{jobs_uniq_id}->{$uid};
}
sub run_slurm {
	my ($self,$cmds,$cpu) = @_;
	my $jobs = $self->submit_jobs();
	warn Dumper $jobs;
 	my $error = $self->run($jobs);
 	#die("eeror on jobs") if $error > 0;
	
}

sub signal_handlers {
    my ($self) = @_;

    $SIG{INT} = sub {
        $self->cleanup();
        die "Interrupted by Ctrl-C\n";
    };

    $SIG{TERM} = sub {
        $self->cleanup();
        die "Terminated\n";
    };

    $SIG{HUP} = sub {
        $self->cleanup();
        die "Terminated\n";
    };
}



has dir_log => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		 my $user = getpwuid($<);
 		my $dir_log = "/data-bipd/data-pure/workspace/slurm/${user}";
 		system("mkdir -p $dir_log");
		return $dir_log;
	},
);

has running_jobs => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		return {} ;	
	},
);
sub cleanup {
	my ($self) = @_;
	print colored("\n--------------------------------------------------------------------\n",'magenta');
	print colored("\n I'm on it! Cleaning up the mess you left when you went away...\n",'magenta');
	print colored("\n--------------------------------------------------------------------\n",'magenta');
	sleep(1);
	my @ids = keys %{$self->running_jobs};
	system("scancel", @ids) if @ids;
}

sub run {
my ($self) = @_;
 my %current_states;
 $self->signal_handlers();
 sleep(5);
 my %done;
  foreach my $uid (@{$self->jobs}){
	my $job = $self->get_job($uid);
	$self->running_jobs->{$job->{id}} ++ if $job->{status} eq "SUBMITED";
  }

  
while (keys %{$self->running_jobs}) {
    my $still_running = 0;
    
    # On interroge sacct pour TOUS les jobs en une seule commande (très rapide)
    # -P (pipe-separated), -X (no sub-steps), --format=JobID,State
    my $ids_string = join(",",keys %{$self->running_jobs});
     
    my $sacct_out = `sacct -j $ids_string --format=JobID,State,Elapsed -P -X 2>/dev/null`;
    
    foreach my $line (split(/\n/, $sacct_out)) {
        next if $line =~ /^JobID/; # Ignorer l'en-tête
        my ($jid, $state,$elapse) = split(/\|/, $line);
        # Si le statut contient des détails (ex: "CANCELLED by 1234"), on ne garde que le premier mot
        $state =~ s/\s+.*// if $state; 
        $current_states{$jid} = $state if $jid;
        $self->submited_jobs->{$jid}->{elapsed} = $elapse;
         if ($state =~ /^(PENDING|COMPLETING|CONFIGURING)$/) {
         	$self->submited_jobs->{$jid}->{status} = "PENDING";
         	 
         }
         elsif ($state =~ /^(RUNNING)$/) {
         	$self->submited_jobs->{$jid}->{status} = "RUNNING";
         }
         elsif ($state =~ /^(COMPLETED)$/){
				if (exists $self->running_jobs->{$jid}){
				#unlink $self->submited_jobs->{$jid}->{log};
         	 	delete $self->running_jobs->{$jid};
         	 	warn "OK ".$self->submited_jobs->{$jid}->{type};
         	 		$self->submited_jobs->{$jid}->{status} = "OK";
				}
			
		}
		else {
			if (exists $self->running_jobs->{$jid}){
				delete $self->running_jobs->{$jid};
         	 	$self->submited_jobs->{$jid}->{status} = "ERROR";
         	 	$self->submited_jobs->{$jid}->{message} = "$state";
			}
		}  
			
    $self->display_jobs_status();
    print "sacct -j $ids_string --format=JobID,State,Elapsed -P -X \n";
	}
	sleep(20) if keys %{$self->running_jobs};
}
 $self->report();

}


sub report {
	my ($self) = @_;
	$self->{jobs_by_type} = {};
 	my $jtypes = $self->types;
    my $pending =0;
    my $error=0;
    my $running=0;
    my $done = 0;
    my @headers = sort {$a cmp $b} keys %{$self->names};
    
    my $tb = Text::Table->new(
        'Type',@headers
    );
   my @rows;
    my $st_error = {};
    foreach my $type (sort{$a cmp $b} keys %{$self->types}){
		 my @row;
		push(@row,$type);
		foreach my $h ( @headers){
			my $job = $self->jobs_by_type_name->{$type}->{$h};
			unless ($job){
				push(@row,colored("-", 'magenta'));
			}
			elsif ($job->{status} eq "ERROR"){
				push(@row,colored("ERROR", 'magenta'));
				push(@{$st_error->{$h}}, colored($type, 'magenta')." -> job: ".$job->{cmd}." log: ".$job->{log}." ".$job->{message});

			}
			else {
				push(@row,colored($job->{elapsed}, 'green'));
				}
			
		}
		push(@rows,\@row);
		
	}
	$tb->load(@rows);
	print "\033[H\033[J";
	print $tb;
	if (keys %{$st_error}){
	foreach my $h ( keys %{$st_error} ){
		print $h."\n";
		print join("\n",@{$st_error->{$h}});
		print "\n";
	}
	}
	else {
		print colored("Job?s done!  Nothing to report! \n","green");
	}
   
		
}
has status_job => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> sub {
		return ["PENDING","RUNNING","OK","ERROR"];
	},
);

sub jobs_by_type {
	my ($self) = @_;
	my $type ={};
	foreach my $uid (@{$self->jobs}){
		my $value  = $self->get_job($uid);
		push (@{$type->{$value->{type}}},$value);
	}
	return $type;
}


sub jobs_by_type_name {
	my ($self) = @_;
	my $type ={};
	my $hash;
	foreach my $uid (@{$self->jobs}){
		my $value  = $self->get_job($uid);
		$hash->{$value->{type}}->{$value->{name}} = $value;
		#push (@{$type->{$value->{type}}->{$value->{name}}},$value->{id});
	}
	return $hash;
	}

sub progress_bar {

    my ($self,$done, $total, $width) = @_;

    $width //= 15;

    return "[".colored("-" x $width, 'white')."]" ." 0 %" if $total == 0;

    my $n = int($done / $total * $width);
	my $pc = int(($done / $total) * 100);
    return "[".colored("*" x $n, 'green')

         . colored("-" x ($width - $n), 'white')."]" ." $pc %";

}

sub display_jobs_status {
	
    my ($self) = @_;
    my $count;
    
    
    
    my $jtypes = $self->jobs_by_type;
    
    my @rows;
    my $pending =0;
    my $error=0;
    my $running=0;
    my $done = 0;
    foreach my $type (sort{$a cmp $b} keys %$jtypes){
		  my $count;
		  my $global_status = "PENDING";
		  my $total = 0;
		  my $nb_ok = 0 ;
		 foreach my $st (@{$self->status_job}){
			$count->{$type}->{$st} = 0;
		 }
		 my $st_error= "";
		 my $st_running= "";
		foreach my $value1 ( @{$jtypes->{$type}}) {
			$total ++;
			#my $value1 = $self->submited_jobs->{$jid};
			my $st =  $value1->{status};
			$count->{$type}->{$st} ++;
			$nb_ok ++ if $st eq "OK";
			my ($sname,@other) = split("!",$value1->{name});
			$st_error .= $sname." " if ($st eq "ERROR");
			$st_running .= $sname."[".$value1->{elapsed}."]-" if ($st eq "RUNNING");
			$error ++ if $st eq "ERROR";
			$running ++ if $st eq "RUNNING";
			$pending ++ if $st eq "PENDING";
		}
		
		if ($count->{$type}->{ERROR} >0 ) {
			$global_status = colored("ERROR : ".$st_error, 'red');
		}
		elsif ($count->{$type}->{RUNNING} >0 ) {
			$global_status = colored($st_running, 'cyan');
		}
		elsif ($nb_ok == $total ) {
			$global_status = colored("DONE", 'green');
		}
		elsif (($nb_ok+$count->{$type}->{ERROR}) == $total ) {
			$global_status = colored("FINISHED ERROR : ".$count->{$type}->{ERROR}, 'red'); 
		}
		
		$done += $nb_ok;
		my $line = [$type];
		push(@$line,$self->progress_bar($nb_ok,$total,50));
		
    
    	push(@$line,$global_status);
    	push (@rows,$line);
		
	}
	
	my $tb = Text::Table->new(
        'Type',
		'Progress',
        'Status'
    );



#my @text_table = (["delete","running"],["prepare primers","waiting"],["Coverage Patients","waiting"],["Coverage Primers","waiting"],["Coverage Exons","waiting"],["CNV","waiting"],["transcripts low","waiting"]);

#push(@text_table, ["delete","running"]);


$tb->load(@rows);
print "\033[H\033[J";
print $tb;

print "Pending : $pending Running:$running Error:$error Done: $done  \n";	


}

sub print_jobs {
	 my ($self) = @_;
	 my $h ;
	foreach my $k (keys %{$self->jobs_by_type_name}){
		foreach my $name (keys %{$self->jobs_by_type_name->{$k}}){
			my ($step,$project) = split("!",$name);
			$h->{$project}->{$k}->{$step} ++;
		}
	}
	
	my @rows;

foreach my $p (sort keys %$h) {

	foreach my $n (sort keys %{$h->{$p}}) {
		my $steps = join(", ", sort { $a cmp $b } keys %{$h->{$p}->{$n}});
		push @rows, [$p, $n, $steps];
	}

}
	# Largeurs des colonnes
	print "\033[H\033[J";
my @width = (

	length("PROJECT"),
	length("TYPE"),
	length("STEPS"),

);

foreach my $row (@rows) {

	for my $i (0..2) {
		$width[$i] = length($row->[$i]) if length($row->[$i]) > $width[$i];
	}

}

my $sep = "+-" .

          join("-+-", map { "-" x $_ } @width) .

          "-+";

print color('bold cyan');

print "$sep\n";
print color('bold cyan');
printf  "| %-*s | %-*s | %-*s |\n",

	$width[0], "PROJECT",
	$width[1], "TYPE",
	$width[2], "STEPS";

print color('reset');

print "$sep\n";

foreach my $row (@rows) {

	printf "| %-*s | %-*s | %-*s |\n",
		$width[0], $row->[0],
		$width[1], $row->[1],
		$width[2], $row->[2];

}

print "$sep\n";
print "\n";
print colored ['bold cyan '], "################################################################################################" ;
			print "\n";


my $yes;
unless ($yes){
	my $choice = prompt("run this/these step(s)   (y/n) ? ");
	die() if ($choice ne "y"); 
}
	#warn Dumper  $self->jobs_by_type_name;
	# die();
}

sub submit_jobs {
 my ($self) = @_;
  my @job_ids =();
foreach my $uid (@{$self->jobs}) {
	warn $uid;
	my $job = $self->get_job($uid);
    my $command = $job->{cmd};
    my $name  = $job->{type}."-".$job->{name};
    my $cpu = $job->{cpu};
    my $dir_log = $self->dir_log;
    my $dependency="";
    if (exists $job->{previous}){
		my $id_previous;
		#$dependency="--dependency=afterok";
		foreach my $uid (@{$job->{previous}}){
			next unless $uid;
			$dependency  ="--dependency=afterok" if $dependency eq "";
			my $previous_jobs = $self->get_job($uid);
			die() unless exists $previous_jobs->{id};
			$dependency.=":".$previous_jobs->{id};
				
			}
	}
    my $sbatch_cmd = "sbatch --parsable "
                   . "--job-name=$name "
                   . "--output=$dir_log/$name-%j.out "
                  # . "---nodes=1 "
                   . "--partition=bipd " 
                   . "--ntasks-per-node=$cpu "
                   . "--mem-per-cpu=3G "
                   .$dependency." "
                   . "--wrap=\"$command\"";
    # Exécution de sbatch et récupération de l'ID du job
    my $job_id = `$sbatch_cmd`;
    chomp($job_id);
    $job->{log} = "$dir_log/$name-$job_id.out";


    # Nettoyage au cas où Slurm renvoie "ID;nom_cluster"
    $job_id =~ s/;.*//; 
	$job->{id} = $job_id;
	$self->submited_jobs->{$job_id} = $job;
    if ($job_id && $job_id =~ /^\d+$/) {
		
        print " $name -> Job soumis avec succès (ID: $job_id)\n";
        $job->{status} = "SUBMITED";
        push(@job_ids, $job_id);
    } else {
		  $job->{status} = "ERROR";
        die "ERREUR FATALE: Impossible de soumettre le job : $command\nRetour Slurm: $job_id\n";
    }
}
 #return \@job_ids;
}


1;

