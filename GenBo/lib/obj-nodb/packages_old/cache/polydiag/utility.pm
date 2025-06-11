package utility;
use strict;
use FindBin qw($RealBin);
use Data::Dumper;


sub return_list_all_transcripts {
	my ($project,$patient,$tr) = @_;
	my $d = $project->getCacheDir();
	my $project_name= $project->name;
	my $patient_name = $patient->name();
	my $string = "";
	$string = $project->noSqlPolydiag()->get($patient->name,"transcripts");

	if ($string){
	return [split(";",$string)];
	}
	else {
		confess();
	}
}
sub return_date_from_file {
	my ($file_out) = @_;
	my $t = (stat $file_out )[9];
	
	my $date = POSIX::strftime( 
             "%d/%m/%y", 
             localtime( 
               		$t
                 )
             );
      return $date;       
}

sub return_last_days_from_file {
	my ($file_out) = @_;
	confess() unless (-e $file_out);
	my $t = (stat $file_out )[9];
	my $t1 = localtime();
	my $date = POSIX::strftime( "%d/%m/%y", localtime( $t ) );
	my $days_difference = int((time - $t) / 86400);
	return ($date, $days_difference);
}

sub get_date {
		my ($patient) = @_;
		my $file_out;
		if ($patient->project->noSqlPolydiag()->exists_db($patient->name)){
			 $file_out = $patient->project->noSqlPolydiag()->dir."/".$patient->name.".".$patient->project->noSqlPolydiag()->extension;
		}
		else {
			#confess();
		}
		
		my $t = (stat $file_out )[9];
		my $t1 = localtime();
		my $date = POSIX::strftime( 
             "%d/%m/%y", 
             localtime( 
               		$t
                 )
             );
             
          
		my $days_difference = int((time - $t) / 86400);
			
           return ($date,$days_difference);
}


sub return_list_variants {
	my ($project,$patient,$tr_id) = @_;
	my $d = $project->getCacheDir();
	my $project_name= $project->name;
	my $patient_name = $patient->name();
	my @vars;
	my $key =$tr_id;
	my $string ="";
	eval {
		$string = $project->noSqlPolydiag()->get($patient->name,"list_$key")."";
	};
	if ($@) {
		my @l;
		return \@l;
	}
	#warn $patient->name." ".$string;
	#confess() unless $string;
	return [split(";",$string)];

	return \@vars;
	
}

sub return_hash_variant {
	my ($project,$vid,$tr_id,$patient,$vquery) = @_;
	my $project_name= $project->name;
	my $patient_name = $patient->name();
		#my $db = return_db($patient);
	my $id = join(";",$tr_id,$vid);
	#my $id = $patient->name()."_".$tr1->id."_".$v->id;
	my $z;
	
		$z = $project->noSqlPolydiag()->get($patient->name,$id);
	unless ($z) {
		return undef;
		die($id."-");
#	
	}
	
	return $z;
}


1;