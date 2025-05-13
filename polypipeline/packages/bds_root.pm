package bds_root;
use Moo;  

use strict;
use FindBin qw($Bin);
use Term::StatusBar;
use lib "$Bin";
use lib ".";
#use root_steps;
use Term::ANSIColor;
use job_bds;
use sample;
use Data::Dumper;
use Proc::Simple;
use Proc::Background; 
use Proc::Daemon;
use Term::ReadKey;
#use Tree::Simple;
#use Graph::Directed;
use POSIX qw{strftime};
require Term::Screen;
use Term::StatusBar;
use UUID 'uuid';
use Carp;

my $bin_dev = "";
has 'patient' => (
	is =>'rw',
);

has "bds_uuid" =>(
	is        => 'rw',
	lazy=>1,
	default   => sub{  
	my $self = shift;
	my $uuid;
 UUID::generate($uuid); 
 my $string_uuid;
 UUID::unparse( $uuid, $string_uuid );
 return $string_uuid."-".$self->project->name."-".time; },
);
has "error" =>(
	is        => 'rw',
	lazy=>1,
	default   => sub{  
	return 0 },
);
has   'yes' =>(
	is =>'rw',
	default => sub {
		return 0 ;
	}
);
has   'limit' =>(
	is =>'rw',
	default => sub {
		return 0 ;
	}
);
has   'cache' =>(
	is =>'rw',
	default => sub {
		return 0 ;
	}
);
has   'queue' =>(
	is =>'rw',
	default => sub {
		return "" ;
	}
);
sub command_after {
	my ($self,$command) = @_;
	$self->{command_after} = [] unless exists $self->{command_after};
	unless ($command){
		return $self->{command_after};
	}
	else {
		push(@{$self->{command_after}},$command);
	}
}
has   'nocluster' =>(
	is =>'rw',
	default => sub {
		return 0 ;
	}
);

has 'report'  =>(
	is        => 'rw',
	default   => sub{ return undef},
);

has   'script_dir' =>(
	is =>'ro',
	default => sub {
		return qq{$Bin/scripts/scripts_pipeline};
	}
);


has 'project' => (
	is =>'ro',
	required=>1,
);

has 'max_cpu' => (
	is	=> 'rw',
); 

has 'argument_patient' => (
	is =>'rw',
);

has 'host' =>(
		is =>'rw',
		lazy=>1,
	default => sub {
		my ($host) = `hostname`;
		chomp($host);
		return $host;
	}
);

has 'bipd' =>(
		is =>'rw',
		lazy=>1,
	default => sub {
		my ($host) = `hostname`;
		my $self = shift;
		if ($self->host =~/morgan/) {
			return 1;
		}
		if ($self->host =~/sanger/) {
			return 1;
		}
		return undef;
	}
);

has 'nproc' => (
		is =>'ro',
		lazy=>1,
	default => sub {
		my $self = shift;
		if ($self->host =~/crick/){
			return 16;
		}
		if ($self->host =~/morgan/) {
			return 20;
		}
		warn 'HOST: '.$self->host."\n\n";
		if ($self->host =~/sanger/) {
			return 20;
		}
		#die($self->host."--");
		my $proc = `cat /proc/cpuinfo | grep processor | wc -l`;
		chomp($proc);
		$proc -= 2;
		$proc =1 if $proc <= 0;
		
		return $proc;
	}
);
has 'local_cpu' => (
		is =>'ro',
		lazy=>1,
	default => sub {
		my $self = shift;

		my $proc = `cat /proc/cpuinfo | grep processor | wc -l`;
		chomp($proc);
		$proc -= 4;
		$proc =1 if $proc <= 0;
		
		return $proc;
	}
);

has 'ppn' => (
		is =>'ro',
	default => sub {
		return 20;
	}
);
has 'root_dir' =>(
	is =>'ro',
	lazy=>1,
	default => sub {
		my $self = shift;
		confess("") unless -e $self->project->buffer->config_path("root","bds");
		return $self->project->buffer->config_path("root","bds");
	}
);
has 'timestamp' => (
	is =>'ro',
	lazy=>1,
	default  => sub {
		return time;
	}
);
has 'unforce' => (
	is        => 'rw',
	default   => 1,
); 
has 'dir_bds' =>(
		is =>'ro',
		lazy=>1,
		default => sub {
		my $self = shift;
		my $dir_root = 
		my $dir = $self->root_dir."/".$self->project->name."_calling_".$self->timestamp;
		system("mkdir -p $dir");
		#mkdir $dir;
		system ("chmod a+rwx $dir");
		return $dir;	
	}
);

has process =>(
		is =>'rw',
		lazy=>1,
		default => sub {
		my $self = shift;
		my $myproc = Proc::Simple->new();
	 	#$myproc->kill_on_destroy(1);
	 	#$myproc->signal_on_destroy("KILL");
	 	return $myproc;
	}
);
has daemon =>(
		is =>'rw',
);
has pid =>(
		is =>'rw',
);

has yaml =>(
		is =>'rw',
		lazy=>1,
		default => sub {
		return undef;
		}
);

has 'jobs' =>(
is =>'rw',
lazy=>1,
		default => sub {
			return [];
		}
);

has 'current_sample' =>(
is =>'rw',
lazy=>1,
		default => sub {
			return undef;
		}
);
sub get_index {
	my ($self) = @_;
	 $self->{index} = 0 unless exists $self->{index}; 
	 return $self->{index} ++;
}
has 'index' =>(
is =>'rw',
lazy=>1,
		default => sub {
			return undef;
		}
);


sub add_sample_with_priority {
	my ($self, $object, $priority, $type) = @_;
	$self->patient($object);
	if ($type) {
		$self->{samples}->{$type.'_'.$object->name."_".$priority} = sample->new(name=>$type.'_'.$object->name, priority=>$priority);
		$self->current_sample($self->{samples}->{$type.'_'.$object->name."_".$priority});
	}
	else {
		$self->{samples}->{$object->name."_".$priority} = sample->new(name=>$object->name, priority=>$priority);
		$self->current_sample($self->{samples}->{$object->name."_".$priority});
	}
}

sub add_sample{
	my ($self,$hash) = @_;
	confess();
	my $patient = $hash->{patient};
	$self->patient($patient);
	
	$self->{samples}->{$patient->name} = sample->new(name=>$patient->name);
	#$self->{samples}->{$patient->name}->index($self->get_index);
	$self->current_sample($self->{samples}->{$patient->name});
	
}
sub job_super_types() {
	confess();
}

sub print_bds{
		my ($self) = @_;
	my $dir = $self->dir_bds;
	open(BDS,">$dir/pipeline.bds") or die($dir);
	
	my @samples = grep { $_->is_pending_jobs()}  $self->samples;
	foreach my $s (@samples){
#			warn $s->name();
		print  BDS $s->bds_sub();
		
	}	
	my $xx = 0 ;
	foreach my $s (@samples){
		$xx =0;
		my $name = $s->task_name();
		print   BDS "par $name();\n sleep(1)\n";
		$xx++;
		
		
	}

	if ($self->command_after ){
		print BDS "wait\n";
		
		foreach my $c (@{$self->command_after}){
			print BDS "task (cpus := 20) { " ;
			print BDS qq{sys $c \n};
			 
				print BDS "\n}\n";
		}
	
	}
}

sub shell2 {
		my ($self) = @_;
	my @samples = grep { $_->is_pending_jobs()}  $self->samples;

	my $hSampleByPriority;
	my @commands;
	foreach my $s (@samples) {
		push(@{$hSampleByPriority->{$s->priority()}}, $s);
	}
	my $xx =0;
	foreach my $priority (sort keys %$hSampleByPriority) {
		
		foreach my $s (@{$hSampleByPriority->{$priority}}) {
			my $jobs = $s->jobs;
			
			foreach my $j (@$jobs){
				warn $j->type." ".$s->name." ".$priority;
				
				next if $j->is_skip();
					foreach my $c (@{$j->cmd}){
						push(@commands,{cmd=>$c,log=>$j->bds_log,name=>$j->type,sample=>$s->name});
					}
			}
			
	}
	
			#die($priority);
	}
	die();
	$self->run(\@commands);
	
	#warn Dumper $hSampleByPriority; die;

}
sub run_one_priority {
	my ($self,$hash) = @_;
	my $list_cmd;
	my $ppn;
	my $nproc = 40;
	my $nb_process = int($nproc/$ppn);
	

}

sub run {
	my ($self,$commands) = @_;
	my $total = scalar @$commands;

# Création de la barre de progression
my $progress = new Term::StatusBar (
                    label => 'jobs Done : ',
                   showTime=>1,
                   subTextAlign =>"center",
                    totalItems => scalar(@$commands),  ## Equiv to $status->setItems(10)
                    #startRow => $row,
                   
 );
# Exécution des commandes séquentiellement
my $index =0;
system("clear");
$progress->start();
foreach my $job (@$commands) {
    # Exécuter la commande
    my $output_file = $job->{log};
    print "\n\033[2K"; # Effacer la ligne précédente
 	print colored::stabilo("green",$job->{name}.":".$job->{sample});
 	my @cmd1 = split("&&",$job->{cmd});
 	my @cs;
 	foreach my $c1 (@cmd1){
 			my $cmd = $c1." >>$output_file 2>&1";
 			push(@cs,$cmd);
 	}
 
 	my $cmd = join (" && ",@cs);
    my $status = system($cmd);
  
    # Vérification des erreurs
    if ($status != 0) {
        my $exit_code = $status >> 8;
        my $command = $job->{cmd};
        print "\nErreur détectée lors de l'exécution de la commande '$command' (code de sortie : $exit_code).\n";
        print "log file is here : $output_file \n";
        exit 1;
    }
    # Mettre à jour la barre de progression
    $progress->update();
}

}
sub print_bds_by_priority{
		my ($self) = @_;
	my $dir = $self->dir_bds;
	open(BDS,">$dir/pipeline.bds") or die($dir);
	my @samples = grep { $_->is_pending_jobs()}  $self->samples;
	foreach my $s (@samples){
		print  BDS $s->bds_sub();
		warn $s->bds_sub();
		
	}
	my $hSampleByPriority;
	foreach my $s (@samples) {
		push(@{$hSampleByPriority->{$s->priority()}}, $s);
	}
	my $xx =0;
	foreach my $priority (sort keys %$hSampleByPriority) {
		
		foreach my $s (@{$hSampleByPriority->{$priority}}) {
			$xx++;
			my $name = $s->task_name();
			print   BDS "par $name();\n sleep(1)\n";
			if ($self->limit > 0 && $xx%$self->limit == 0 ){
			print BDS "wait\n";
		}
		}
		print BDS "wait\n";
	}

	#warn Dumper $hSampleByPriority; die;

}

sub launch_bds_daemon{
		my ($self) = @_;
	#return $self->launch_shell()  if $self->nocluster;
	$self->print_bds();
	$self->launch_bds_daemon_common();
}

sub launch_bds_daemon_by_priority{
		my ($self) = @_;
	#return $self->launch_shell()  if $self->nocluster;
	return $self->shell2()  if $self->nocluster;
	$self->print_bds_by_priority();
	$self->launch_bds_daemon_common();
}



sub launch_bds_daemon_common{
		my ($self) = @_;
	my $dir = $self->dir_bds;
	$| = 1;
	print colored ['black ON_BRIGHT_GREEN'],"RUNNING BDS .... ";
	print color 'reset';
	print "\n";
	#system("cd $dir; bds -reportYaml -quiet -s cluster  pipeline.bds");
	my $bds_exe = $self->project->buffer->software("bds-cluster")."  ";
	
	$bds_exe = $self->project->buffer->software("bds") if $self->nocluster > 0;
	return $self->shell2 if $self->nocluster == 2;
	return $self->launch_shell3 if $self->nocluster == 3;
	
	warn "mode  $bds_exe =>  ".$self->nocluster;#  if $self->nocluster == 2;
	my $libe = "$bds_exe  -reportYaml    pipeline.bds";
	$libe= "$bds_exe ".$self->queue." -reportYaml    pipeline.bds" if $self->cache ==1;
	my $daemon = Proc::Daemon->new(
        work_dir => "$dir",
   		exec_command =>$libe,
    );
	
	
	my $pid = $daemon->init;
	$self->daemon($daemon);
	$self->pid($pid);
	
	
	#$self->process->start(" cd $dir ; $bds_exe -reportYaml  $cluster  pipeline.bds");
	my $nb;
	while ( $daemon->Status($pid)){
		$nb ++;
		my $z =0;
		while ($z< 30){
			$z++;
			#$s->next();
			print ".";
			sleep(1);
		}
		print "\n";
		$self->print_status_bds();# if $nb %3 == 0;
		print "pid : $pid \n";
#		 if ($self->nocluster){
#		
#		my @ps =  `ps -edf`;
#		chomp(@ps);
#		warn Dumper (@ps);
#		my $results;
#		foreach my $l1 (@ps){
#			if ($l1 =~/$pid/){
#				push(@$results,"$l1");
#			}
#		}
#		if ($results){
#			print join("\n",$results);
#		}
#		else {print "problem with bds process\n;"}
#		 }
	}
	$daemon->Kill_Daemon($pid);
	print colored::stabilo("CYAN","----------------------------------\n",1);
	print colored::stabilo("CYAN","-----------FINISHED---------------\n",1);
	print colored::stabilo("CYAN","----------------------------------\n",1);
	print "\n";
	print "\n";	
	
	

	#$self->process(undef);
	my ($logfiles) = $self->report_final();
	
	 if ($self->report ){
	 	my $report_file = $self->report_target();
	 	push(@$logfiles,$report_file);
	 }
	 if (@$logfiles){
	 	print colored::stabilo("CYAN","-----------FINISHED---------------\n",1);
	 	print colored::stabilo("GREEN","Report file :  ".join(" ",@$logfiles),1);
		print colored::stabilo("CYAN","----------------------------------\n",1);
		print "\n";
		print "\n";	
	 }
	
	 
}
has   'sample_error' =>(
	is =>'rw',
	default => sub {
		return {} ;
	}
);

has   'start_jobs' =>(
	is =>'rw',
);

sub launch_shell2{
		my ($self) = @_;
		die();
#	my $max_proc = $self->nproc;
#	my @samples = grep { $_->is_pending_jobs()} $self->samples;
#	my @categories =  $samples[0]->list_categories();
#	foreach my $c (@categories){
#	
#		foreach my $s (@samples) {
#		
#			my $error = 0;
#			
#			my $jobs = $s->get_jobs_category({category=>$ctaegory});
#			my $nbr = scalar(grep{$_->is_run} @$jobs);
#			next if $nbr ==0;
#			
#				warn "-*-*-*-*-*-*-*-*-*-*-*--*\n";
#				warn "$c : ".scalar(@{$h->{$c}->{jobs}}) ."\n";
#				warn "-*-*-*-*-*-*-*-*-*-*-*--*\n";
#				
#				my $cpu = $self->nproc;
#				my $real_fork = int($cpu/$ppn);
#				#warn $real_fork;
#				$real_fork =1 if $real_fork ==0;
#		 		$real_fork = $self->nproc if $real_fork>$self->nproc;
#
#				
#				foreach my $j (@{$jobs}){
#					next if $error == 1;
#					next if $j->is_skip;
#					$self->print_status_bds();
#					my $pid = $pm->start and next;
#			  			foreach my $c ($j->command_bds){
#  							system ($c);
#  							#warn $c;
#  							
#  					}
#  					sleep 1;
#  					my $status = 0;
#  					$status = 1 if $j->is_ok;
#  					$status = 1 ;
#  						$self->print_status_bds();
#				}
#			} #end categories
#			
#	}#end samples
	
}


sub launch_shell3{
		my ($self) = @_;
	my $max_proc = $self->nproc;
	my @samples = grep { $_->is_pending_jobs()} $self->samples;
	my $scr = Term::Screen->new();
	my $status_samples;
	my $pos = 1;
	if (scalar(@samples) > 1){
	  $status_samples = new Term::StatusBar (
                    label => 'Samples : ',
                   showTime=>1,
                   subTextAlign =>"left",
                   totalItems => scalar(@samples),  ## Equiv to $status->setItems(10)
                    startRow => $pos,
    );
    $pos +=2;
}
	$scr->clrscr();
	
	$status_samples->start() if $status_samples;	
	my $error_jobs = [];
	my $ok_jobs ;
	my $htime;

	
	foreach my $s (@samples) {
			my @categories1 =  $s->list_categories();
			my $h = $samples[0]->list_hashes_categories();
			my @categories;
			foreach my $c (@categories1){
				my @jobs = grep {!($_->is_skip)} @{$h->{$c}->{jobs}};
				next unless @jobs;
				push(@categories,$c);
			}
			my $error = 0;
		
					my $status_cat = new Term::StatusBar (
                     label => "Steps   : ",
                   subTextAlign =>"left",
                     startRow => $pos,
                   totalItems => scalar(@categories),  ## Equiv to $status->setItems(10)
                  
    			);
    			$status_cat->start();
    			$pos +=2;
    			
    				
			foreach my $c (@categories){
			$htime->{$s->{name}}->{$c}->{start} =  time; # same as ( epoch => time() );
			
			#	warn "-*-*-*-*-*-*-*-*-*-*-*--*\n";
			#	warn "$c : ".scalar(@{$h->{$c}->{jobs}}) ."\n";
			#	warn "-*-*-*-*-*-*-*-*-*-*-*--*\n";
				my $ppn= $h->{$c}->{ppn};
				my $cpu = $self->nproc;
				my $real_fork = int($cpu/$ppn);
				#warn $real_fork;
				$real_fork =1 if $real_fork ==0;
		 		$real_fork = $self->nproc if $real_fork>$self->nproc;
		 			my $pm = new Parallel::ForkManager($real_fork);
		 		my @jobs = grep {!($_->is_skip)} @{$h->{$c}->{jobs}};
					my $status_jobs = new Term::StatusBar (
                     label => "Jobs       ",
                   subTextAlign =>"left",
                     startRow => $pos,
                   totalItems => scalar(@jobs),  ## Equiv to $status->setItems(10)
                   subText =>"$c",
                   showTime=>1,
    			);
    			$status_jobs->start();
    			my $running_jobs =0;
    			
				$pm->run_on_finish(
				sub {
					my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data ) = @_;
						my $merror =0;
						$merror ++ unless (exists $data->{status});
						$merror ++ if $data->{status} ne 1;
						#$scr->at($pos+2,1)->bold()->puts("hi!")->normal();
						#$self->print_status_bds();
						$running_jobs --;
						 $ok_jobs ++ if $merror == 0;
  						$scr->at($pos+2,1)->bold()->puts("Running Jobs : ".$running_jobs."    ")->normal();
  						$scr->at($pos+3,1)->bold()->puts("Finished Jobs : ".$ok_jobs."    ")->normal();
  						$scr->at($pos+4,1)->bold()->puts("Error Jobs : ".scalar(@$error_jobs)."    ")->normal();
  							if ($merror){
  						push(@$error_jobs,{sample=>$s->name,cat=>$c,name =>$data->{job}->name});
  						my $nb =5;
  						foreach my $e (@$error_jobs){
  							my $string =  $e->{sample}." ".$e->{cat}." ".$e->{name};
  							$scr->at($pos+$nb,1)->bold()->puts(colored::stabilo("GREEN",$string,1))->normal();
  							$nb ++;
  							last if $nb >10;
  						}
  							}
						$error += $merror;
						$status_jobs->update();
						
				}
				);
				
				foreach my $j (@jobs){
					next if $error == 1;
				#	$self->print_status_bds();
				$running_jobs ++;
				$scr->at($pos+2,1)->bold()->puts("Running Jobs : ".$running_jobs."    ")->normal();
					my $pid = $pm->start and next;
			  			foreach my $c ($j->command_bds){
  							system ($c);
  						
  							#warn $c;
  							
  					}
  				#	sleep 2;
  					my $status = 0;
  					$status = 1 if $j->is_ok;
  					$status = 1 ;
  					$pm->finish(0,{status=>$status,job=>$j});
  					
				}
				$pm->wait_all_children();
				#$htime->{$s->{name}}->{$c}->{end} =  DateTime->now; # same as ( epoch => time() );;
				$htime->{$s->{name}}->{$c}->{elapsed} = ( time - $htime->{$s->{name}}->{$c}->{start});# - time; # same as ( epoch => time() );;
				$status_cat->update();
				$htime->{$s->{name}}->{$c}->{elapsed_text} = strftime("\%H:\%M:\%S", gmtime($htime->{$s->{name}}->{$c}->{elapsed}));
			} #end categories
		$status_samples->update if $status_samples;	
		
	}#end samples
	$scr->clrscr();
		my $tb;
		my $lines;
		foreach my $s (@samples) {
			my @line;
			push(@line,$s->name);
			my @categories =  $s->list_categories();
			
			unless ($tb){
					 $tb = Text::Table->new(
						("Sample",@categories,"total")
    					);
			}
			my $htotal;
			foreach my $c (@categories){
				if (exists $htime->{$s->{name}}->{$c}->{elapsed_text}){
					push(@line,colored::stabilo("white",$htime->{$s->{name}}->{$c}->{elapsed_text},1)) ;
					$htotal += $htime->{$s->{name}}->{$c}->{elapsed};
				}
				else {
					push(@line,"-") 
				}
				}
				push(@line, strftime("\%H:\%M:\%S", gmtime($htotal)));
				push(@$lines,\@line);
			}
		 $tb->load(@$lines);
	print $tb;
	print "\n";
	#warn Dumper $htime;
	#$scr->at($pos+20,1);
	
}


sub launch_shell{
		my ($self) = @_;
	
	my $max_proc = $self->nproc;
	
	my @samples = grep { $_->is_pending_jobs()} $self->samples;
	 my $list_running_cmd;
	$| =1;
	my @categories =  $samples[0]->list_categories();
		
		
	foreach my $category (@categories) {
		my $run;
		$run->{ppn} =0;
		$run->{jobs} =[];
		foreach my $s (@samples){
			my $jobs = $s->get_jobs_category({category=>$category});
		
			my $nbr = scalar(grep{$_->is_run} @$jobs);
			next if $nbr ==0;
			foreach my $job (@$jobs) {
				next if $job->is_skip();
				$job->{sample} = $s;
				if ($run->{ppn} + $job->ppn >$max_proc){
					push(@$list_running_cmd,$run);
					$run = {};
					$run->{ppn} =0;
					$run->{jobs} =[];
				}
				push(@{$run->{jobs}},$job);
				$run->{ppn} += $job->ppn;
				$run->{sample} = $s;
				
			} 
		}
		if ($run->{jobs} ){
			push(@$list_running_cmd,$run);
		}
		
		
		
	}

	$self->start_jobs(time);		
	foreach my $run (@{$list_running_cmd}){
			my $dir = $self->dir_bds;
			my @process;
			foreach my $j (@{$run->{jobs}}) {
				my $sname = $j->{sample}->name;
			
				next if exists $self->sample_error->{$sname};
				my @com = $j->command_bds();
				my $proc;
				 $proc->{file} = "$dir/".$j->name.".sh";
				open(OUT,">". $proc->{file}) or die ;
				print OUT join("\n",@com);
				close (OUT);
				system("chmod a+x ".$proc->{file});
				$proc->{process} = Proc::Simple->new();
				$proc->{process}->kill_on_destroy(1);
				$j->{shell_start} = time;
				$proc->{process}->start($proc->{file});
				$proc->{job}= $j;
				push(@process,$proc);
			}
		while ($self->still_runnning(\@process)){
			my $z =0;
			while ($z< 30){
				$z++;
				print ".";
				sleep(1);
				$z =50 unless $self->still_runnning(\@process);
			}
			print "\n";
			
			foreach my $p (@process){
				
				if (exists  $p->{error}){
					print colored ['black ON_BRIGHT_RED']," ERROR !!!!!  ".$p->{job}->name();
					print color 'reset';
					print "\n";
					$self->sample_error->{$p->{job}->{sample}->name} ++;
				}
		}
		$self->print_status_bds();
		}		
		
	
	}#end for each command
	
	print colored::stabilo("CYAN","----------------------------------\n",1);
	print colored::stabilo("CYAN","-----------FINISHED---------------\n",1);
	print colored::stabilo("CYAN","----------------------------------\n",1);
	print "\n";
	print "\n";
	$self->	report_final_shell;
	my $logfiles = [];
	if ($self->report ){
	 	my $report_file = $self->report_target();
	 	push(@$logfiles,$report_file);
	 }
	 if (@$logfiles){
	 	print colored::stabilo("CYAN","-----------FINISHED---------------\n",1);
	 	print colored::stabilo("GREEN","Report file :  ".join(" ",@$logfiles),1);
		print colored::stabilo("CYAN","----------------------------------\n",1);
		print "\n";
		print "\n";	
	 }
	 
	 
	
	}
	
sub report_final_shell{
		my ($self) = @_;
	my @samples = grep { $_->is_pending_jobs()} $self->samples;
	 my $list_running_cmd;
	my @categories =  $samples[0]->list_categories();
	my @header =("patients",@categories,"total");
	my @lines;
	my %cat_time;
		foreach my $s (@samples){
			my @line;
			push(@line,$s->name);
			my $sample_time =0;
			foreach my $category (@categories) {
				my $jobs = $s->get_jobs_category({category=>$category});
				my $time = 0;
				my $status = 0;
				my $text = "NC";
				foreach my $j (@$jobs){
						
						if (exists $j->{shell_end})	{
							$time+=abs($j->{shell_end} - $j->{shell_start});
							
						}
						$status += $j->{shell_error} if exists  $j->{shell_error}	;	
				}
				
					$sample_time += $time;
					if ($time == 0) {
							push(@line,colored::stabilo("cyan","NC",1));
					}
					elsif ($status == 0){
							push(@line,colored::stabilo("green","ERROR",1));
						push(@line,colored::stabilo("green",strftime("\%H:\%M:\%S", gmtime($time)),1));
					}
					else {
						push(@line,colored::stabilo("red",strftime("\%H:\%M:\%S", gmtime($time)),1));
					}
					$cat_time{$category} += $time;
			}
			push(@line,colored::stabilo("white",strftime("\%H:\%M:\%S", gmtime($sample_time)),1));
			push(@lines,\@line);
		}
		my @line =(" ");
		my $total;
	foreach my $category (@categories) {
		$total += $cat_time{$category};
			push(@line,colored::stabilo("white",strftime("\%H:\%M:\%S", gmtime($cat_time{$category})),1));
	}
	
	push(@line,colored::stabilo("white",strftime("\%H:\%M:\%S", gmtime($total)),1));
		push(@lines,\@line);
	my $tb = Text::Table->new(
	@header,
    );
    	$tb->load(@lines);
	print $tb;
	print "\n";
	
}

	
	sub still_runnning {
		my ($self,$process) = @_;
		my $r;
		foreach my $p (@$process) {
			 $r = 1 if $p->{process}->poll();
			 next if $p->{process}->poll();
			 
			 next if exists  $p->{job}->{shell_end};
			 
			 $p->{job}->{shell_end} = time ;
			  $p->{error} = 1 if $p->{process}->exit_status() ne 0;
			  $p->{job}->{shell_error} = 1 if $p->{process}->exit_status() ne 0;
		}
	
		return $r;
	}


sub launch_bds{
		my ($self) = @_;
	$self->print_bds();
	my $dir = $self->dir_bds;
	$| = 1;
	die();
	print colored ['black ON_BRIGHT_GREEN'],"RUNNING BDS .... ";
	print color 'reset';
	print "\n";

	#system("cd $dir; bds -reportYaml -quiet -s cluster  pipeline.bds");
	my $cluster = "-s cluster";
	
	$cluster ="" if $self->nocluster;
	my $bds_exe = $self->project->buffer->software("bds");

	$self->process->start(" cd $dir ; $bds_exe -reportYaml -quiet $cluster  pipeline.bds");
	my $pid = $self->process->pid;
	#$self->process->start(" cd $dir ; $bds_exe -reportYaml  $cluster  pipeline.bds");
	while ( $self->process->poll()){
		my $z =0;
		while ($z< 30){
			$z++;
			#$s->next();
			print ".";
			sleep(1);
		}
		print "\n";
		$self->print_status_bds();
		print "pid : $pid \n";
		 if ($self->nocluster){
		
		my @ps =  `ps -edf`;
		chomp(@ps);
		my $results;
		foreach my $l1 (@ps){
			if ($l1 =~/$pid/){
				push(@$results,"$l1");
			}
		}
		if ($results){
			print join("\n",$results);
		}
		else {print "problem with bds process\n;"}
		 }
	}
	
	print colored::stabilo("CYAN","----------------------------------\n",1);
	print colored::stabilo("CYAN","-----------FINISHED---------------\n",1);
	print colored::stabilo("CYAN","----------------------------------\n",1);
	print "\n";
	print "\n";	
	
	

	$self->process(undef);
	my ($logfiles) = $self->report_final();
	
	 if ($self->report ){
	 	#my $report_file = $self->report_target();
	 	#push(@$logfiles,$report_file);
	 }
	 if (@$logfiles){
	 	print colored::stabilo("CYAN","-----------FINISHED---------------\n",1);
	 	print colored::stabilo("GREEN","Report file :  ".join(" ",@$logfiles),1);
		print colored::stabilo("CYAN","----------------------------------\n",1);
		print "\n";
		print "\n";	
	 }
	
	 
}



sub report_final {
		my ($self) = @_;
	my @yamlfiles;
	my $dir = $self->dir_bds;;
	 @yamlfiles = `ls -r $dir/*.yaml`;
	chomp(@yamlfiles);
	my $file_yaml = $yamlfiles[-1];
	if (-e  "$file_yaml"){
		$self->yaml($file_yaml);
    	system("$Bin//scripts/scripts_pipeline/parse_bds_report.pl $file_yaml");
	}
	my @samples = grep { $_->is_pending_jobs()} $self->samples;
	my @job_log;
	foreach my $s ($self->samples){
		push(@job_log, grep{$_->isLogging} @{$s->jobs});
	}
	
	foreach my $j (@job_log){
		warn $j->fileout;
		#system("cat ".$j->fileout);
	}
	return \@job_log;
}

sub report_target{
		my ($self) = @_;
	my @yamlfiles;
	 my $project_name = $self->project()->name();
	 my $dir = $self->project->project_log();
	 my $file_log = $dir."/quality_check.html";
	  my $file_log2 = $dir."/mendelian.log";
	   my $c = "$Bin//scripts/scripts_pipeline/quality_check.pl -project=$project_name -fork=".$self->local_cpu;
#	 my $c = "$Bin//scripts/scripts_pipeline/quality_check_target_gene.pl -project=$project_name > $file_log";
	system($c);	 
#	if ($p->isFamilial)	{
#		 my $c = "$Bin//scripts/scripts_pipeline/quality_check_project.pl -project=$project_name  ";
#	}
	
	return $file_log;
  
    	
}

sub print_status_bds {
		my ($self) = @_;
	print "\n";
	#[$status,$tcurrent,$terror,$nbok,$remaining];
	my @header = ("sample","status","Running","OK","error","remaining");
	my @lines;
	$self->start_jobs(time) unless defined $self->start_jobs();
	my $tt = abs(time-$self->start_jobs);
		foreach my $s ($self->samples){
		
			my @line;
			push(@line,$s->name);
			my $total =0;
			my $totalr= 0;
			push(@line,@{$s->status_jobs});
			
			push(@lines,\@line);
			
		}
	 	my $tb = Text::Table->new(
			@header,
    	);
   
    	$tb->load(@lines);
	print $tb;
	
	print "-------------------------------- ".strftime("\%H:\%M:\%S", gmtime($tt)). "--------------------------------\n";
	print "\n";
	print "\n";
	my @files;
	foreach my $s ($self->samples){
		
		next unless @{$s->log_for_error_jobs()};
		$self->error(1);
		push(@files,@{$s->log_for_error_jobs()});
	}
	if (@files){
		print "you have error in  jobs.  log file are here : \n";
		print join("\n",@files)."\n";
	}
	print "\n-----\n";
	
}


sub clean {
		my ($self) = @_;
	my @samples = grep { $_->is_pending_jobs()} values %{$self->{samples}};
		
		my $delete_files;
		my $delete_prod_files;
		foreach my $s (@samples){
		my @dp;
		my @dt;
		 @dt = grep {-e  $_ or -l $_ } @{$s->delete_tmp_files};
		 @dp = grep {-e  $_} @{$s->delete_prod_files};
		 push(@dt,grep {-l  $_} @{$s->delete_prod_files});
		 chomp(@dt);
		 chomp(@dp);
		 $delete_files->{pipeline}->{$s->name} = \@dt if @dt;
		 $delete_files->{prod}->{$s->name} = \@dp if @dp;
		}
		
		return if scalar(keys %$delete_files) == 0;#&& scalar (@dp)==0; 
	my @header = ("sample","pipeline","prod");
	my @lines;
	foreach my $s (@samples){
		my @line;
		push(@line,$s->name);
		
		#my @dt = grep {-e  $_ or -l $_} @{$s->delete_tmp_files};
		my $text = colored::stabilo("green", "0",1);
		if (exists  $delete_files->{pipeline}->{$s->name}){
			$text = colored::stabilo("cyan",scalar(@{ $delete_files->{pipeline}->{$s->name}}). " Files ",1);
		}
		push(@line,$text);
		 $text = colored::stabilo("green", "0",1);
		
		if (exists $delete_files->{prod}->{$s->name}){
			$text="";
			my @type_file= ("bam","cov","vcf.gz","gvcf.gz","log");
			my @res;
			foreach my $t (@type_file){
				my $f = grep {$_=~/$t/} @{$delete_files->{prod}->{$s->name}};
				if ($f) {
					push(@res,colored::stabilo("red",uc($t),1));
				}	
			}
			$text = join("-",@res);
			}
			push(@line,$text);
		
		
		
		push(@lines,\@line)
	}	
	 	my $tb = Text::Table->new(
			( @header),
    	);
   
   	$tb->load(@lines);
	print $tb;
	print "\n";
	print "\n";
	unless ($self->yes){
	print "delete or backup all these files   (y/n) ? ";
	my $key = $self->key;
	print "\n";
	die() if ($key ne "y"); 
	}
foreach my $s (@samples){
		#log_file 
		my $name = $s->name;
		print colored ['black ON_BRIGHT_YELLOW'], "delete $name."   ;
		
	if (exists  $delete_files->{pipeline}->{$s->name}){
		 
			foreach my $f (@{ $delete_files->{pipeline}->{$s->name}}){
				$self->clean_temporary_file($f);
			}
				
		 }
		 
		  print  "...";
		if (exists  $delete_files->{prod}->{$name}){			
		 	my $arg = join(" ",@{ $delete_files->{prod}->{$s->name}});
		 	my $bin_dev = $self->script_dir;
		 			my $cmd = "$bin_dev/rm_prod.pl $arg";
				system($cmd);
				foreach my $f (@{ $delete_files->{prod}->{$s->name}})
					{
							if (-e $f ){
								warn "unable to delete  : $f   !!!!!!!!!!!!!! \n check permissions" ;
								die();
							}
					}			
		 }	
		 print "Done ";
		 print color 'reset';
		print "\n";
}
}

sub clean_temporary_file {
	my ($self,$f) = @_;
	my @tt = `ls -d $f*`;
	chomp(@tt);
	foreach my $ff (@tt){
		if (-d $ff){
						
						system("rm $ff/* ; rmdir $ff ");
					}
	}
		system("rm ".$f."* 2>/dev/null");
}

sub clean_error{
		my ($self) = @_;
	my @samples = values %{$self->{samples}};
	
	foreach my $s  (@samples){
		foreach my $j (@{$s->jobs}){
			next unless $j->is_run;
			next unless -e $j->bds_start;
			unless (-e $j->bds_ok){
				warn "delete temporary file : ".$j->name;		
				# if (-e $j->fileout){
				#warn $j->fileout;
						$self->clean_temporary_file( $j->fileout) if $j->fileout =~/pipeline/;
						unlink $j->fileout  if -e $j->fileout;
				 #}
				
			}
		}
	}
}

sub samples {
		my ($self) = @_;
	return (sort {$a->{index} <=> $b->{index}}values %{$self->{samples}});
}


has   'priority_name' =>(
	is =>'rw',
	default => sub {
		
		return ["Chromosomes","Project","Patients"]
	}
);

sub print_all_steps_by_prority{
		my ($self) = @_;
	my @sams = $self->samples;
	
	my $hCat;
	my $types;
	my %jdj;
	foreach my $sam (@sams) {
		my $type =  $sam->priority();
		foreach my $j (@{$sam->jobs}){ 
			my $cat = $j->category;
			unless (exists $jdj{$cat}){
				push(@{$hCat->{$type}},$cat);
				$jdj{$cat}++;
			}
			#$hCat->{$type}->{$j->category}->{pos} = $njob unless exists ;
			
			}
		
	}
	
	my $hSampleByPriority;
	foreach my $s (@sams) {
		push(@{$hSampleByPriority->{$s->priority()}}, $s);
	}
	my @types = sort {$a <=> $b} keys %$hSampleByPriority;
	
	my $nb_type= 0;
	foreach my $type (@types){
		$nb_type ++;
		my $samples = $hSampleByPriority->{$type};
		my $tb;
		my @lines;
		foreach my $s (sort @$samples){
			my @line;
			push(@line,$s->name);
			my $total =0;
			my $totalr= 0;
			my $categories = $hCat->{$type};
			my $text = $self->priority_name->[($type-1)];
	#		$text = "Project" if ($type == 2);
	#		$text = "Patient" if ($type == 3);
			
		 	$tb = Text::Table->new( (colored::stabilo("cyan", $text, 1),  @$categories, "total") ) unless $tb; # if ($type == 1);
		foreach my $cat ( @$categories) {
				my $jobs = $s->get_jobs_category({category=>$cat});
			my $nb = scalar(@$jobs);
		$total += $nb;
		my $nbp =0;
			unless ($jobs){
					push(@line, colored::stabilo("white","NONE",1) );
					next;
				}
				my $nbj = scalar(@$jobs);
				my $nbr = scalar(grep{$_->is_run} @$jobs);
				$totalr +=$nbr;
				if ($nbr ==0){
					push(@line, colored::stabilo("white","SKIP",1) );
				}
				elsif ($nbj eq $nbr){
					my $text = "$nbr" ;
					$nbp++;
					$text = "Planned" if $nbr== 1;
					push(@line, colored::stabilo("green","$text",1));
				}
				else {
					my $text = "$nbr/$nbj" ;
					push(@line, colored::stabilo("green","$text",1));
				}
			}
			push(@line, colored::stabilo("cyan","$totalr/$total",1) );
			push(@lines,\@line);
		}
		
		$tb->load(@lines);
		print $tb;
	
		print "\n ";
		print "\n";
		if ($nb_type< scalar(@types)){
		unless ($self->yes){
		  print "Press <Enter> or <Return> to continue: ";
		 my $resp = <STDIN>;
		}
		}
	}
	
	
	
}


sub print_all_steps{
		my ($self) = @_;
	my @sams = $self->samples;#keys %{$self->samples};
	my $hCat;
	foreach my $sam (@sams) {
		foreach my $j (@{$sam->jobs}){ $hCat->{$j->category}++; }
	}
	my @categories = keys %{$hCat};
	
	my $n = 1;
	foreach my $cat (sort keys %{$hCat}) {
		$hCat->{$cat} = $n;
		$n++;
	}

	my @lines;
	my @types;
	foreach my $s (sort @sams){
		my @line;
		push(@line,$s->name);
		my $total =0;
		my $totalr= 0;
		#$s->reconstruct_dependancy();
		my $jobs = $s->jobs();
		
		my $nb = scalar(@$jobs);
		$total += $nb;
		my $nbp =0;
		foreach my $cat (sort @categories){
			my $jobs = $s->get_jobs_category({category=>$cat});
			unless ($jobs){
				push(@line, colored::stabilo("white","NONE",1) );
				next;
			}
			my $nbj = scalar(@$jobs);
			my $nbr = scalar(grep{$_->is_run} @$jobs);
			$totalr +=$nbr;
			if ($nbr ==0){
				push(@line, colored::stabilo("white","SKIP",1) );
			}
			elsif ($nbj eq $nbr){
				my $text = "$nbr" ;
				$nbp++;
				$text = "Planned" if $nbr== 1;
				push(@line, colored::stabilo("green","$text",1));
			}
			else {
				my $text = "$nbr/$nbj" ;
				push(@line, colored::stabilo("green","$text",1));
			}
		}
		#map{push(@types,$_->type)} @$jobs unless @types;
		push(@line, colored::stabilo("cyan","$totalr/$total",1) );
		push(@lines,\@line);
	}
 	my $tb = Text::Table->new(
		("patients", sort @categories,"total"),
   	);
   	$tb->load(@lines);
	print $tb;
	print "\n";
	print "\n";
}

sub key {
	my ($self) = @_;
	my $key;
	 ReadMode 0;
	 while (not defined ($key = ReadKey(-1))) {
        # No key yet
    }
    ReadMode 1;
    return $key;
}
1;