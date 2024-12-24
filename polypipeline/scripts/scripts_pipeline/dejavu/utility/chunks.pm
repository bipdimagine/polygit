package chunks;

use FindBin qw($Bin);
use strict;
use JSON;
use Archive::Tar;
use Fcntl ':flock';
use File::NFSLock qw(uncache);
use List::Util qw(shuffle);
use Data::Dumper;
use Archive::SevenZip;
sub hg38_tar {
	return "/data-beegfs/tmp/sereal";
	return  "/data-isilon/DejaVu/HG38/variations/projects/rocks.sereal.tar/";
	
}
sub hg38_rocks_projects {
	return "/data-beegfs/tmp/rocks_sereal";
}
sub chromosomes {
return [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y','MT'];
}
sub divide_by_chunks{
	my ($chr,$version) = @_;
	$version = "HG38" unless $version;
	my $json_chr_length = qq{{"HG19":{"6":"171115067","11":"135006516","9":"141213431","15":"102531392","14":"107349540","1":"249250621","8":"146364022","17":"81195210","7":"159138663","13":"115169878","MT":"16569","22":"51304566","Y":"59373566","3":"198022430","21":"48129895","18":"78077248","5":"180915260","X":"155270560","20":"63025520","16":"90354753","2":"243199373","12":"133851895","19":"59128983","4":"191154276","10":"135534747"},
"HG38":{"21":"46709983","Y":"57227415","18":"80373285","3":"198295559","22":"50818468","MT":"16569","20":"64444167","16":"90338345","X":"156040895","5":"181538259","2":"242193529","10":"133797422","4":"190214555","19":"58617616","12":"133275309","11":"135086622","6":"170805979","14":"107043718","15":"101991189","9":"138394717","17":"83257441","8":"145138636","1":"248956422","13":"114364328","7":"159345973"}}};

		my $region;
		my $size = 2_500_000;
		my $h = decode_json($json_chr_length);
		my $from =1;
    	my $to = $h->{$version}->{$chr};
    	my $start;
    	my $end;
    	
    	while ($from < $to){
        $start = $from;
        $end = $from + $size;
        if ($end > $to){
            $end = $from + $size;
        }
        my $id_chunk = $chr.".".$from.".".$end;
        $region->{$id_chunk}->{start} = $from;
        $region->{$id_chunk}->{end} = $end;
        $from = $end;
          $region->{$id_chunk}->{variants} = ["toto"];
       }
       my $r;
		foreach my $id (sort{$region->{$a}->{start} <=> $region->{$b}->{start}} keys %{$region}){
			push(@$r,{id=>$id,start=>$region->{$id}->{start},end=>$region->{$id}->{end},tabix=>$chr.":".$region->{$id}->{start}."-".$region->{$id}->{end}});
		}
 		return $r; 
}
sub chunks_and_tree {
	my ($project,$version) = @_;
	 $version = $project->genome_version_generic unless $version;
	 my $chunks;
	 my $tree;
	foreach my $chr (@{$project->getChromosomes}){
			my @rs = sort {$a->{start} <=> $b->{start}} @{chunks::divide_by_chunks($chr->name,$version)};
			$chunks->{$chr->name} = \@rs;
			$tree->{$chr->name} = chunks::construct_tree($chunks->{$chr->name});
}
	return ($chunks,$tree);
}
sub get_region {
	my ($tree,$pos) = @_;
	my $results = $tree->fetch($pos,$pos+1);
	confess(@$results) if scalar(@$results) > 1;
	return $results->[0];
}
	
	
sub construct_tree {
	my ($regions) = @_;
	my $tree =  Set::IntervalTree->new;
 	for my $r (@$regions) {
 		 $tree->insert($r, $r->{start}, $r->{end});
 	}
 		return $tree;
	}


sub return_rocks_id_from_genbo_id {
	my ($id) = @_;
	my ($chr,$pos,$ref,$alt) = split("_",$id);

	
	return return_rocks_id($pos,$ref,$alt);
}

sub stringify_pos {
	my ($pos) = @_;
	return ($pos,sprintf("%010d", $pos));
}

sub return_rocks_id {
	my ($pos,$ref,$alt) = @_;
	my $l1 = length($ref);
	my $l2 = length($alt);
	return  (stringify_pos($pos)."!".$alt) if ($l1 == 1 && $l2 ==1);
	my $seqid = $alt;
	if ($alt =~ /del/ ){
		return  (stringify_pos($pos)."!".$alt);
	}
	elsif ($ref =~ /del/ ){
		return  (stringify_pos($pos)."!".$ref);
	}
	elsif ($alt=~ /inv/ ){
		return  (stringify_pos($pos)."!".$alt);
	}
	elsif ($alt=~ /ins/ ){
		return  (stringify_pos($pos)."!".$alt);
	}
	elsif ($alt=~ /dup/ ){
		return  (stringify_pos($pos)."!".$alt);
	}
	
	
	die($ref." ".$alt) if $alt =~ /INV/;
	#die(ref." ".$alt) if $alt =~ /DEL/;s
	die(ref." ".$alt) if $alt =~ /DUP/;
	
	
	if ($l1 ==1 && $l2 > 1){
		
		$seqid = "+".substr($alt, 1);
		return  (stringify_pos($pos)."!".$seqid);
	}
	elsif ($l1 >1 && $l2 ==1){
		$ref = substr($ref, 1);
		$seqid = ($l1 -1);
		return  (stringify_pos($pos)."!".$seqid);
	}
	 elsif ($l1 >1 && $l2 == $l1 && $l2>1 ){
	 	
		$ref = substr($ref, 1);
		$alt = substr($alt, 1);
		$seqid = "$ref*$alt";
		return  (stringify_pos($pos)."!".$seqid);
	}
	
	else {
		return  (stringify_pos($pos)."!".$ref."*".$alt);
		confess($l1." ".$l2);
	}
	
}


#    while (-e $lock_file) {
#        # Attendre que le fichier de verrouillage disparaisse
#        sleep(10);
#        warn "wait $lock_file";
#        $x++;
#        die($project_name." ".$lock_file) if $x > 3;
#    }
#    # Créer un fichier de verrouillage
#    open my $lock_fh, '>', $lock_file or die "Impossible de créer le fichier de verrouillage : $! $lock_file";
#    print $lock_fh "Verrouillage pour $tar_file\n";
#    close $lock_fh;
sub create_table {
	my ($dbh) = @_;

# Créer une table pour stocker les BLOBs si elle n'existe pas
$dbh->do("CREATE TABLE IF NOT EXISTS results_blob (project_name TEXT PRIMARY KEY, data BLOB)");
$dbh->do("CREATE INDEX IF NOT EXISTS idx_process_id ON results_blob (project_name)");

}
sub put_rocks_with_lock_file {
	 my ($no,$project_name,$data) = @_;
	 my $lock_file = $no->dir.".".$no->name.".lock";
    warn $lock_file;
       if (my $lock = new File::NFSLock {
  file      => $lock_file,
  lock_type => LOCK_EX,
  #blocking_timeout   => 600,      # 10 sec
  stale_lock_timeout => 5 * 60, # 30 min
}) {
	
	$no->put($project_name,$data);
	$no->rocks->compact_range();
	$lock->unlock();
	}else{
  		die "I couldn't lock the file [$File::NFSLock::errstr]";
	}
}



sub add_to_tar_with_lock_file {
    my ($dir_tar,$region_id,$project_name,$data,$htime,$debug) = @_;
    system("mkdir $dir_tar") unless -e  $dir_tar;
    my $tar_file = $dir_tar.$region_id.".tar";
    my $lock_file = $tar_file.".lock";
    if (my $lock = new File::NFSLock {
  file      => $lock_file,
  lock_type => LOCK_EX,
  #blocking_timeout   => 600,      # 10 sec
  stale_lock_timeout => 5 * 60, # 30 min
}) {
    uncache($tar_file);
    my $x;
	my $tar_file2 = "/tmp/pipeline/".$region_id.".tar.new";
	my $t = time;
	system("cp $tar_file $tar_file2");
	 $htime->{cp} += abs(time -$t);
    # Ouvrir et modifier le fichier TAR
    my $tar = Archive::Tar->new;
     $t = time;
    
    $tar->read($tar_file2) if -e $tar_file;
    $htime->{read} += abs(time -$t);
    $t = time;
    if ($tar->contains_file($project_name)){
			 $tar->replace_content ($project_name, $data);
	}
	else {
		
		
    	$tar->add_data ( $project_name, $data);
    }
    $htime->{add} += abs(time -$t);
    $t = time;
    $tar->write($tar_file2);
     $htime->{write} += abs(time -$t);
      $t = time;
    system("mv $tar_file2 $tar_file");
     $htime->{mv} += abs(time -$t);
   
    $lock->unlock();
}else{
  die "I couldn't lock the file [$File::NFSLock::errstr]";
}
}
sub save_chromosome_chunks {
	my ($project,$chr,$chunks,$version) = @_;
	confess() unless $version;
	my $htime = {};
	my $dirout = $project->deja_vu_rocks_project_dir($version);
	 my $encoder = Sereal::Encoder->new({compress=>Sereal::SRL_ZSTD,compress_threshold=>0});
	system("mkdir $dirout/".$chr->name) unless -e "$dirout/".$chr->name;
	my $tar_dir = $dirout."/".$chr->name."/";
	$tar_dir =~ s/\/\//\//g;
	my @vht =  shuffle @{$chunks} ;
	
	foreach my $region (@vht){
		next  unless exists $region->{variants};
		next unless @{$region->{variants}};
		
			my $t = time;
			my $data= $encoder->encode($region->{variants});
			my $debug;
			 $htime->{encode} += abs(time -$t);
			 $t =time;
			 
			add_to_tar_with_lock_file($tar_dir,$region->{id},$project->name,$data,$htime,$debug);
			 $htime->{total} += abs(time -$t);
	}
	warn Dumper $htime;
}

sub exists_project {
	 my ($project,$version) = @_;
	 my $dirout = $project->deja_vu_rocks_project_dir($version);
	 my $tar_file = $dirout."/"."projects.tar";
	 return unless  -e $tar_file;
	 my $tar = Archive::Tar->new;
	  $tar->read($tar_file);
	   return $tar->contains_file($project->name);
	  
}
sub save_final {
	 my ($project,$version) = @_;
	 confess() unless $version;
	 my $dirout = $project->deja_vu_rocks_project_dir($version);
	  add_to_tar_with_lock_file($dirout."/","projects",$project->name,time);
}
sub save_chunks {
	 my ($project,$chunks) = @_;
	
	 my $dirout = $project->deja_vu_rocks_project_dir();
	my @chrt = shuffle @{$project->getChromosomes};	
		foreach my $chr (@chrt){
			save_chromosome_chunks($project,$chr,$chunks->{$chr->name});
		}
		save_final($project,$project->deja_vu_rocks_project_dir());
		#chunks::add_to_sqlite($dir38."/","projects",$project->name,"!!!");
}


1;