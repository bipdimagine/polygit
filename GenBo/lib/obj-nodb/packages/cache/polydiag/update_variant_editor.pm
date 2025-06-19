package update_variant_editor;
use strict;
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../";
use lib "$Bin/../../";
use Storable qw/thaw freeze/;
#use BioTools;
use POSIX;
use Time::Piece;
use List::Util qw( max min sum);
use Carp qw(confess croak);
use JSON::XS;
use Data::Dumper;
#use Bio::DB::HTS::VCF;
#use Bio::DB::HTS::Tabix;
use Tabix;
use polyweb_dude;
use Number::Format qw(:subs);
use List::MoreUtils qw{ natatime };
use CGI qw/:standard :html3/;
use Number::Format qw(:subs);
sub init_bundle_infos{
my $self = shift;
#
		
}
 my $himpact_sorted = {
	"high" => "4",
	"moderate" =>"3",
	"low" =>"1",
};

my $minus = qq{<span class="glyphicon glyphicon-minus" aria-hidden="true"></span>};
my $plus  = qq{<span class="glyphicon glyphicon-plus" aria-hidden="true"></span>};
 my $boy = qq{<img src="https://img.icons8.com/color/24/000000/boy.png">};
 my $girl 	= 		qq{<img src="https://img.icons8.com/color/24/000000/girl.png">};
 my $female = qq{<i class="fa fa-venus fa-2x" aria-hidden="true" style="color:pink"></i>};
  my $male = qq{<i class="fa fa-mars fa-2x" aria-hidden="true" style="color:blue"></i>};
  
sub printBadge {
	my ($value,$types) = @_;
	my $color = "#4CAF50";
	 #$color = "#D62C1A" if $value > $type[0] ;
	 $color = "#FF8800" if $value > $types->[0] ;
	 $color = "#FF0025" if $value > $types->[1] ;
	 return qq{<span class="badge badge-success badge-xs" style="border-color:$color;background-color:#FFFFFF;color:$color;font-size :8px;">$value</span>} ;
}

sub printBadgeWithDesc{
	my ($description,$value,$types) = @_;
	my $color = "#4CAF50";
	 #$color = "#D62C1A" if $value > $type[0] ;
	 $color = "#FF8800" if $value > $types->[0] ;
	 $color = "#FF0025" if $value > $types->[1] ;
	 return qq{<span class="badge badge-success badge-xs" style="border-color:$color;background-color:#FFFFFF;color:$color;font-size :8px;">$description:$value</span>} ;
}

sub printSimpleBadge {
	my ($value) = @_;
	my $color = "black";
	 return qq{<span class="badge badge-success badge-xs" style="border-color:black;background-color:#FFFFFF;color:$color;font-size :8px;">$value</span>} ;
}

sub printBlueSimpleBadge {
	my ($value) = @_;
	my $color = "blue";
	 return qq{<span class="badge badge-success badge-xs" style="border-color:black;background-color:#FFFFFF;color:$color;font-size :8px;">$value</span>} ;
}

sub printInvBadge {
        my ($value,$types) = @_;

        my $color = "#4CAF50";
         #$color = "#D62C1A" if $value > $type[0] ;
         $color = "#FF8800" if $value < $types->[0] ;
         $color = "#FF0025" if $value < $types->[1] ;

         $value ="-" unless defined $value;
         $color = "#FF0025" if $value eq "-" ;

         return qq{<span class="badge badge-success badge-xs" style="border-color:$color;background-color:#FFFFFF;color:$color;font-size :8px;">$value</span>} ;
}


sub printInvBadgeWithDesc {
	my ($description,$value,$types) = @_;
	my $color = "#4CAF50";
	 #$color = "#D62C1A" if $value > $type[0] ;
	 $color = "#FF8800" if $value < $types->[0] ;
	 $color = "#FF0025" if $value < $types->[1] ;
	 $value ="-" unless defined $value;
	 $color = "#FF0025" if $value eq "-" ;
	
	 return qq{<span class="badge badge-success badge-xs" style="border-color:$color;background-color:#FFFFFF;color:$color;font-size :8px;">$description:$value</span>} ;
}
sub printInvButton {
	my ($value,$types,$text,$othercss) = @_;
	my $btn_class  = qq{class= "btn btn-xs btn-primary " style="background-color: #D0D0D0;font-size: 7px;font-family:  Verdana;color:black"};
	$btn_class = qq{class= "btn btn-xs btn-primary " style="background-color: #FF8800;font-size: 7px;font-family:  Verdana;;color:black"} if $value <  $types->[0] ;
	$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: #e74c3c;font-size: 7px;font-family:  Verdana;color:white"}  if $value < $types->[1] ;
	 $value ="-" unless defined $value;
	$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: #e74c3c;font-size: 7px;font-family:  Verdana;"}  if $value eq "-" ;
	$text  = $value unless $text;
	return  qq{<button type="button" $btn_class $othercss>$text</button>};
}
sub printButton {
	my ($value,$types,$text,$othercss,$text_alert) = @_;

	my $btn_class  = qq{class= "btn btn-xs btn-primary " style="background-color: #D0D0D0;font-size: 7px;font-family:  Verdana;color:black"};
	$btn_class = qq{class= "btn btn-xs btn-primary " style="background-color: #FF8800;font-size: 7px;font-family:  Verdana;;color:white"} if $value >=  $types->[0] ;
	$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: #e74c3c;font-size: 7px;font-family:  Verdana;color:white"}  if $value >=  $types->[1] ;
	 $value ="-" unless defined $value;
	 
	$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: red;font-size: 7px;font-family:  Verdana;"}  if $value eq "-" ;
	$text  = $value unless $text;
	return  qq{<button onClick="alert('$text_alert')" "type="button" $btn_class $othercss>$text</button>} if ($text_alert);
	return  qq{<button type="button" $btn_class $othercss>$text</button>};
}
##################
#value for hash variant 
###################



sub vcosmic {
	my ($v,$hvariation) = @_;
	 	$hvariation->{cosmic} = $v->cosmic();
		$hvariation->{value}->{cosmic} = $v->cosmic();
		my $id  = $v->cosmic();
		if ($id){
			my ($id1,$nb) = split(":",$id);
			my $text = $id1;
			my $url = "http://cancer.sanger.ac.uk/cosmic/mutation/overview?id="; 
			 $url = "http://grch37-cancer.sanger.ac.uk/cosmic/ncv/overview?id="  if ($id1 =~/COSN/);
			  $text = "$id1 [$nb]" if $nb;
			$id1 =~s/COSM//;
			$id1 =~s/COSN//;
		
			
			
			
			
			$hvariation->{cosmic} = qq{<a href=\"$url$id1\" target="_blank">$text</a>};
			$hvariation->{html}->{cosmic} = qq{<a href=\"$url$id1\" target="_blank">$text</a>};
		}
	
}

sub vname2 {
	my ($v,$hvariation) = @_;
	my $vn = $v->vcf_id;
	my $vclinvar_text = '';
#	$vclinvar_text = "<br><span style='color:red'>clinvar:".$v->clinvar_id()."</span>" if ($v->clinvar_id());
	my $vt = 'SNP' if ($v->isVariation());
	$vt = 'CNV' if ($v->isCnv());
	$vt = 'DEL' if ($v->isDeletion());
	$vt = 'INS' if ($v->isInsertion());
	$vn =~ s/_/-/g;
	$vn=~ s/chr//;
	my $project_name = $v->project->name();
	my $gnames = join(';',map{$_->name} @{$v->getGenes});
	my $dataset = "?dataset=gnomad_r2_1";
	$dataset = "?dataset=gnomad_r4" if $v->getProject->getVersion() =~ /HG38/;
	my $url_gnomad = qq{https://gnomad.broadinstitute.org};
	if ($v->isCnv()){
		my $len = $v->length();
		my $chr_id = $v->getChromosome->id();
		my $start = $v->start();
		my $end = $v->end();
		$vn = $v->name;
		my $pp = $v->getChromosome->name."-".$v->start."-".$v->end;
		$hvariation->{var_name} = printSimpleBadge(qq{<a href='$url_gnomad/region/$pp$dataset' target = '_blank' style="color:black;vertical-align:middle;">style="color:black;vertical-align:middle;"><div>$vn$vclinvar_text<br><small>$vt</small></div></a> });;
		$hvariation->{html}->{var_name} = printSimpleBadge(qq{<a href='$url_gnomad/region/$pp$dataset' target = '_blank' style="color:black;vertical-align:middle;" style="color:black;vertical-align:middle;"><div>$vn$vclinvar_text<br><small>$vt</small></div></a> });
	}
	elsif ($v->name() =~ m/rs/) {
		my $vname = $v->name();
		my $vname2;
		if ($v->isVariation) {
			$vname2 = $vn;
		}
		else {
			$vname2 = $v->getChromosome->id().'-'.$v->start().'-'.$v->end();
		}
		$hvariation->{var_name} = printSimpleBadge(qq{<a href='$url_gnomad/variant/$vn$dataset' target = '_blank' style="color:black;vertical-align:middle;"><i class="fa fa-users fa-2x" style="color:coral;padding-bottom:3px;"></i><div>$vname2<br><span style='color:red;'>$vname</span>$vclinvar_text<br><small>$vt</small></div></a> });	;#qq{<a href='http://www.ncbi.nlm.nih.gov/snp/?term=$vn' target = '_blank'>$vn</a> };	
		$hvariation->{html}->{var_name}  = printSimpleBadge(qq{<a href='$url_gnomad/variant/$vn$dataset' target = '_blank' style="color:black;vertical-align:middle;"><i class="fa fa-users fa-2x" style="color:coral;padding-bottom:3px;"></i><div>&nbsp$vname2<br><span style='color:red;'>$vname</span>$vclinvar_text<br><small>$vt</small></div></a> });	;
	}
	elsif ($v->name() =~ m/del|ins|dup/ or $v->isInsertion() or $v->isDeletion()) {
		my $vname = $v->name();
		if ($v->isDeletion() and $v->length() > 20) {
			$vn = $v->getChromosome->id().'-'.$v->start().'-del-'.$v->length();
		}
		my $len = length($v->var_allele())-1;
		if ($v->isInsertion() and $len > 20) {
			$vn = $v->getChromosome->id().'-'.$v->start().'-ins-'.$len;
		}
		my $pp = $v->getChromosome->name."-".$v->start."-".$v->end;
		$hvariation->{var_name} = printSimpleBadge(qq{<a href='$url_gnomad/region/$pp$dataset' target = '_blank' style="color:black;vertical-align:middle;"><div>$vn$vclinvar_text<br><small>$vt</small></div></a> });;
		$hvariation->{html}->{var_name} = printSimpleBadge(qq{<a href='$url_gnomad/region/$pp$dataset' target = '_blank' style="color:black;vertical-align:middle;"><div>$vn$vclinvar_text<br><small>$vt</small></div></a> });
	}
	elsif ($v->isVariation()) {
		my $vname = $v->name();
		$hvariation->{var_name} = printSimpleBadge(qq{<a href='$url_gnomad/variant/$vn$dataset' target = '_blank' style="color:black;vertical-align:middle;"></i><div><span style='color:black;'>$vname</span>$vclinvar_text<br><small>$vt</small></div></a> });	;#qq{<a href='http://www.ncbi.nlm.nih.gov/snp/?term=$vn' target = '_blank'>$vn</a> };	
		$hvariation->{html}->{var_name}  = printSimpleBadge(qq{<a href='$url_gnomad/variant/$vn$dataset' target = '_blank' style="color:black;vertical-align:middle;"><div><span style='color:black;'>$vname</span>$vclinvar_text<br><small>$vt</small></div></a> });	;
	}
	else {
		my $pp = $v->getChromosome->name."-".$v->start;
		$hvariation->{var_name} = printSimpleBadge(qq{<a href='$url_gnomad/region/$pp$dataset' target = '_blank' style="color:black;vertical-align:middle;"><div>$vn$vclinvar_text<br><small>$vt</small></div></a> });;
		#$hvariation->{value}->{var_name} = $vn;
		$hvariation->{html}->{var_name} = printSimpleBadge(qq{<a href='$url_gnomad/region/$pp$dataset' target = '_blank' style="color:black;vertical-align:middle;"><div>$vn$vclinvar_text<br><small>$vt</small></div></a> });
	}
}

sub vname {
	my ($v,$hvariation) = @_;
	my $dataset = "?dataset=gnomad_r2_1";
	$dataset = "?dataset=gnomad_r4" if $v->getProject->getVersion() =~ /HG38/;
	my $vn=$v->vcf_id;
	$vn =~ s/_/-/g;
	$vn=~ s/chr//;
	my $project_name = $v->project->name();
	my $gnames = join(';',map{$_->name} @{$v->getGenes});
	if ($v->isCnv()){
		my $len = $v->length();
		my $chr_id = $v->getChromosome->id();
		my $start = $v->start();
		my $end = $v->end();
		$vn = $v->name;
		$hvariation->{var_name} = printSimpleBadge(qq{<a>style="color:black">$vn</a> });;
		$hvariation->{html}->{var_name} = printSimpleBadge(qq{<a style="color:black">$vn</a> });
	}
	elsif ($v->name() =~ m/del|ins|dup/) {
		my $vname = $v->name();
		$hvariation->{var_name} = printSimpleBadge(qq{<a href='https://gnomad.broadinstitute.org/variant/$vn$dataset' target = '_blank' style="color:black"><i class="fa fa-users fa-2x" style="color:coral"></i>&nbsp$vname</a> });	;#qq{<a href='http://www.ncbi.nlm.nih.gov/snp/?term=$vn' target = '_blank'>$vn</a> };	
		$hvariation->{html}->{var_name}  = printSimpleBadge(qq{<a href='https://gnomad.broadinstitute.org/variant/$vn$dataset' target = '_blank' style="color:black"><i class="fa fa-users fa-2x" style="color:coral"></i>&nbsp$vname</a> });	;
	}
	elsif ($v->name() =~ m/rs/) {
		my $vname = $vn;
		$hvariation->{var_name} = printSimpleBadge(qq{<a href='https://gnomad.broadinstitute.org/variant/$vn$dataset' target = '_blank' style="color:black"><i class="fa fa-users fa-2x" style="color:coral"></i>&nbsp$vname</a> });	;#qq{<a href='http://www.ncbi.nlm.nih.gov/snp/?term=$vn' target = '_blank'>$vn</a> };	
		$hvariation->{html}->{var_name}  = printSimpleBadge(qq{<a href='https://gnomad.broadinstitute.org/variant/$vn$dataset' target = '_blank' style="color:black"><i class="fa fa-users fa-2x" style="color:coral"></i>&nbsp$vname</a> });	;
	}
	else {
		my $pp = $v->getChromosome->name."-".$v->start;
		$hvariation->{var_name} = printSimpleBadge(qq{<a href='https://gnomad.broadinstitute.org/region/$pp$dataset' target = '_blank' style="color:black">$vn</a> });;
		$hvariation->{value}->{var_name} = $vn;
		$hvariation->{html}->{var_name} = printSimpleBadge(qq{<a href='https://gnomad.broadinstitute.org/region/$pp$dataset' target = '_blank' style="color:black">$vn</a> });
	}
}

sub vsequencing  {
	my ($v,$hvariation,$patient,$debug) = @_;
	my @asequence_info;
	my @apc;
	my @methods;
	my $nb_methods;
	my $max_pc =-1;
	my $max_dp =  -1;
	 my $pid = $patient->id;
	 
	my $isCachedPatient;
	$isCachedPatient = 1 if (ref($patient) =~ /Cache/);
	
	foreach my $method (keys %{$v->sequencing_infos->{$pid}}) {
			next if $method eq "ok";
			next if $method eq "max";
			next if $method eq "values";
			#next if $method eq "dude";

		push(@methods,$method);
		
		if ($isCachedPatient) {
			my $array = $v->sequencing_infos->{$pid}->{$method};
			 
			my $sequence_info ; 
			my $pc ="-";		
			if ($v->getNbAlleleRef($patient,$method) eq "?"){
				$sequence_info = "??";
			}
			else {
				$sequence_info .=	$array->[2]."(";
				#$sequence_info = "ho(" if $all_annex->{ho};
				$sequence_info .= $v->getNbAlleleRef($patient,$method)."/".$v->getNbAlleleAlt($patient,$method).")";
			
				}
				$sequence_info = $method.":".$sequence_info;
				 $pc = $method.":".$v->getRatio($patient,$method)."%";
				push(@apc,$pc);
				push(@asequence_info,$sequence_info);
				$nb_methods ++;
			}
			
			 if ($v->validation_method eq "sanger" ) {
			 	#$sequence_info = "-";
			 	push(@asequence_info,"-");
			 }
		}
		
		$hvariation->{max_dp} = $v->getDP($patient);
		$hvariation->{value}->{max_dp} = $hvariation->{max_dp};
		$hvariation->{html}->{max_dp} = $hvariation->{max_dp};
		
		if ($isCachedPatient) {
			$hvariation->{max_pc} = $v->getRatio($patient);;
			$hvariation->{value}->{max_pc} =$hvariation->{max_pc};
			$hvariation->{html}->{max_pc} = $hvariation->{max_pc};
			
			$hvariation->{value}->{ngs} = \@asequence_info;
		
			$hvariation->{ngs} = printSimpleBadge(join("<br>",@asequence_info));
			$hvariation->{value}->{ngs} = \@asequence_info;
			$hvariation->{html}->{ngs} = printSimpleBadge(join("<br>",@asequence_info));
			$hvariation->{ratio} =  printSimpleBadge(join("<br>",@apc));
			$hvariation->{value}->{ratio} =  \@apc;
			$hvariation->{html}->{ratio} =  printSimpleBadge(join("<br>",@apc));
		}
		
		$hvariation->{caller} =  printSimpleBadge(join("<br>",@methods));
		$hvariation->{value}->{caller} =  \@methods;
		$hvariation->{html}->{caller} =  printSimpleBadge(join("<br>",@methods));
	
}

sub value_html {
	my ($hvariation,$type,$v1,$v2) = @_;
	$v2 = $v1 unless defined $v2;
	
	$hvariation->{value}->{$type} = $v1;
	$hvariation->{html}->{$type} = $v2 ; 
}

sub vdivers {
	my ($v,$hvariation) = @_;
		value_html($hvariation,"allele",$v->getSequence());
		value_html($hvariation,"ref_allele", $v->ref_allele());
		value_html($hvariation,"genomique", $v->getChromosome()->name.":".$v->start,printSimpleBadge($v->getChromosome()->name.":".$v->start));
		value_html($hvariation,"genomique_value", $v->getChromosome()->name.":".$v->start);
		value_html($hvariation,"start",$v->start);
		value_html($hvariation,"end",$v->end);
		value_html($hvariation,"chromosome",$v->getChromosome()->name);
	
}

sub vvarsome {
	my ($hvariation,$debug) = @_;
	return if exists $hvariation->{html}->{varsome}; 
	my $url = qq{https://varsome.com/variant/hg38/}.$hvariation->{value}->{gnomad_id};
	my $text =qq{<button dojoType="dijit.form.Button"   iconClass='https://img.icons8.com/ios-filled/24/000000/vimeo.png' onclick='window.open($url,"_blank")' style="color:black"></button>};
	my $text =qq{<a  type="button" class="btn btn-primary btn-xs" href="$url" target="_blank">V</a>};
	value_html($hvariation,"varsome",$url,$text);
}

sub alamut_link_js {
	my ($hvariation,$patient,$debug) = @_;
	my $start = $hvariation->{value}->{start};
	my $a0 = $hvariation->{value}->{ref_allele};
	
	my $a1 = $hvariation->{value}->{allele};
	my $chr_name = $hvariation->{value}->{chromosome};

	my $qq5 = qq{	<button    class="alamutView3" onClick ="displayInAlamut('$chr_name',$start,['$a0','$a1']);"></button>};
		value_html($hvariation,"alamut",$chr_name.":".$start.$a0.">".$a1,$qq5);
#	die();
}
sub valamut_igv {
	my ($v,$hvariation,$patient,$debug) = @_;
		
		my $bam;
		eval { $bam = $patient->getBamFileName(); };
		if ($@) { $bam = undef; }
		my $start = $v->start();
		my $chr = $v->getChromosome();
		my $chr_name = $v->getChromosome();
			
		if (not $bam or not -e $bam) {
			value_html($hvariation,"igv",$chr.":".$start,qq{<img src="https://img.icons8.com/ios/24/000000/select-none.png">});
			value_html($hvariation,"alamut",$chr.":".$start,qq{<img src="https://img.icons8.com/ios/24/000000/select-none.png">});
			return;
		}
	
		my $qq4 = qq{	<div  data-dojo-type="dijit/Toolbar"><button dojoType="dijit.form.Button"   iconClass="igvIcon" onClick ="displayInIGV('$chr_name',$start,$start);"></button></div>};
		my $qq4 = qq{	<button  class="igvIcon2" onClick ="displayInIGV('$chr',$start,$start);">&nbsp;&nbsp;&nbsp;&nbsp;</button>};
		#value_html($hvariation,"igv",$chr.":".$start,$qq4);
		
		
		#$hvariation->{igv} = $qq4; 		
				my @bams;
				my @names;
					foreach my $p (@{$patient->getFamily->getPatients()}){
						next unless -e $p->getBamFileName;
						push(@bams,$p->bamUrl);
						push(@names,$p->name());
					}
					
				my $f =  join(";",@bams);#$patient->{obj}->bamUrl;;
					my $l = $v->getChromosome()->name.":".$v->start;
					my $v1 = $hvariation->{ref_allele}."/".$hvariation->{allele};	
					my $gn = $patient->project->getVersion();
					my $project_name = $patient->project->name;
					my $pnames = join(";",@names);
					
					#$text =qq{<button dojoType="dijit.form.Button"   iconClass='igvIcon' onclick='launch_web_igv("$project_name","$pnames","$f","$l","$v","$gn")' style="color:black">toto</button>};
					#$text =qq{<button onclick='alert("coucou");' style="color:black">toto</button>};
					my $text =qq{<button class='igvIcon2' onclick='launch_web_igv("$project_name","$pnames","$f","$l","$v1","$gn")' style="color:black"></button>};
					#value_html($hvariation,"igv",$chr->name.":".$start,$text);
					#value_html($hvariation,"igv",$chr->name.":".$start,$text);
					$pnames = join(",",@names);
					my $text =qq{<button class='igvIcon2' onclick='launch_web_igv_js("$project_name","$pnames","$f","$l","$v1","$gn")' style="color:black"></button>};
					value_html($hvariation,"igv",$chr->name.":".$start,$text);
						
			my $a0 = $v->ref_allele;
			my $a1 = $v->var_allele;
			my $chr_name = $chr->fasta_name();		
			my $qq5 = qq{	<button    class="alamutView3" onClick ="displayInAlamut('$chr_name',$start,['$a0','$a1']);"></button>};
			value_html($hvariation,"alamut",$chr->name.":".$start.$a0."_".$a1,$qq5);
			
		
		
			
} 

sub put_text_minus {
	my ($value) = @_;
	return $value  if $value ;
	return   qq{<span class="glyphicon glyphicon-minus" aria-hidden="true"></span>};
	
}



sub table_gnomad {
	my ($v) = @_;
	my $cgi          = new CGI();
	my $value = $v->getGnomadAC;
	
	my $color = "grey";
	if ($value <5){
		$color = "red";
	}
	elsif ($value < 20){
		$color = "orange";
	} 
	elsif ($value < 50){
		$color = "blue";
	} 
	
	my $pp = $v->getChromosome->name."-".$v->start;
	my $dataset = "?dataset=gnomad_r2_1";
	$dataset = "?dataset=gnomad_r4" if $v->getProject->getVersion() =~ /HG38/;
	my $href = qq{https://gnomad.broadinstitute.org/region/$pp$dataset};
	my $vn=$v->vcf_id;
	$vn =~ s/_/-/g;
	$vn=~ s/chr//;

	if ($v->getGnomadAC && $v->getGnomadAC > 0){
			$href = qq{https://gnomad.broadinstitute.org/variant/$vn$dataset};
	}

	my $html= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary table_gnomad",style=>"box-shadow: 1px 1px 6px $color ;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
	$html.= $cgi->start_Tr();
	
	$html.= $cgi->th("AC");
	$html.= $cgi->th("Ho");
	$html.= $cgi->th('<i class="fa fa-mars" > </i> ') if  ($v->getChromosome->name eq "X" or $v->getChromosome->name eq "Y");
	$html.= $cgi->th("Max");
	$html.= $cgi->th("Min");
	$html.= $cgi->th("AN");
	$html.= $cgi->end_Tr();
	
	$html.= $cgi->start_Tr();
		
	
	$html.= $cgi->td(abutton($href,put_text_minus($v->getGnomadAC)));
	$html.= $cgi->td(abutton($href,put_text_minus($v->getGnomadHO)));
	$html.= $cgi->td(abutton($href,put_text_minus($v->getGnomadAC_Male))) if  ($v->getChromosome->name eq "X" or $v->getChromosome->name eq "Y");
	my $text =  qq{<span class="glyphicon glyphicon-minus" aria-hidden="true"></span>};
	$text  = $v->max_pop_name."<br>".sprintf("%.4f", $v->max_pop_freq)  if $v->max_pop_name;
	$html.= $cgi->td(abutton($href,$text));
	
	$text =qq{<span class="glyphicon glyphicon-minus" aria-hidden="true"></span>}; 
	$text  = $v->min_pop_name."<br>".sprintf("%.4f", $v->min_pop_freq) if $v->min_pop_name;
	$html.= $cgi->td(abutton($href,$text));
	$html.= $cgi->td(abutton($href,put_text_minus($v->getGnomadAN)));
	$html.= $cgi->end_Tr();
	$html.=$cgi->end_table();

	return $html;
}


sub vgnomad {
	my ($v,$hvariation) = @_;
	 
		value_html($hvariation,"gnomad",undef,table_gnomad($v));
		return if  exists $hvariation->{value}->{ac};
		#return;
	 	my $max  ="-";
	 	$max = $v->max_pop_name.":".sprintf("%.4f", $v->max_pop_freq ) if $v->max_pop_name;
	 	value_html($hvariation,"max_pop",$max, printSimpleBadge($max));
	 	
	 	my $min = "-";
		$min = $v->min_pop_name.":".sprintf("%.4f", $v->min_pop_freq ) if $v->min_pop_name;
		value_html($hvariation,"min_pop",$min, printSimpleBadge($min));
		
		my $freq_ho = "-"; 
		$freq_ho = sprintf("%.4f", $v->frequency_homozygote ) if $v->frequency_homozygote ;
		value_html($hvariation,"freq_ho",$freq_ho, printSimpleBadge($freq_ho));
		
		value_html($hvariation,"ac",$v->getGnomadAC+0, printInvButton($v->getGnomadAC,[200,10]));
		value_html($hvariation,"an",$v->getGnomadAN+0,  printInvBadge($v->getGnomadAN,[0,0]));
		
	
		 
		value_html($hvariation,"ac_ho",$v->getGnomadHO+0, printInvButton($v->getGnomadHO,[50,5]));
		
		
		if  ($v->getChromosome->name eq "X" or $v->getChromosome->name eq "Y") {
			$hvariation->{ac_ho} = printInvButton($v->getGnomadHO,[50,5], $v->getGnomadHO." -  ".qq{&nbsp<i class="fa fa-mars" > </i> &nbsp;} .$v->getGnomadAC_Male);#printInvButton($v->getGnomadAC_Male,[50,5]) unless ($v->is_in_pseudoautosomal_region );  
			value_html($hvariation,"ac_ho",$v->getGnomadHO.":".$v->getGnomadAC_Male+0,printInvButton($v->getGnomadHO,[50,5], $v->getGnomadHO." -  ".qq{&nbsp<i class="fa fa-mars" > </i> &nbsp;} .$v->getGnomadAC_Male));
		}
}


sub vclinical_local {
		my ($v,$hvariation) = @_;
	 my $v1 = $v->score_clinical_local();
	 $hvariation->{clinical_local}  = "" ;
	 if ($v1){
	 	 	$hvariation->{clinical_local}  ++ ;
	 	 #	$hvariation->{clinvar_alert}  ++ ;
	 	 	my $cm = $v->comment_clinical_local();
	 	 	my $txt =  qq{<span class="badge badge-warning" style="font-size:8px">}.$hvariation->{clinvar}."</span>".qq{<span class="badge badge-warning"  onClick='alert ($cm)' style="font-size:8px">Local Clinical</span>};
	 	 	#$hvariation->{clinvar}  = " $txt ";
	 		$hvariation->{scaled_score} = 4 if $hvariation->{freq_level} <= 2 ;	
	 		
	 }
}

sub check_is_hgmd_dm_for_gene {
	my ($project,$hvariation,$gene) = @_;
	return $hvariation->{value}->{dm_for_this_gene} if (exists $hvariation->{value}->{dm_for_this_gene});
	if ($hvariation->{value}->{dm}) {
		my $chr = $project->getChromosome($hvariation->{value}->{chromosome});
		my $hgmd_id = $hvariation->{value}->{hgmd_id};
		my $g = $project->newGene($gene->{id});
		if ($chr->is_hgmd_DM_for_gene($hgmd_id, $g)) {
			$hvariation->{value}->{dm_for_this_gene} = 1;
			return 1;
		}
		else {
			$hvariation->{value}->{dm_for_this_gene} = undef;
			$hvariation->{value}->{dm} = undef;
			$hvariation->{value}->{hgmd} = '';
			$hvariation->{html}->{hgmd} = '';
		}
	}
	return;
}

sub check_is_clinvar_pathogenic_for_gene  {
	my ($project,$hvariation,$gene) = @_;
	return $hvariation->{value}->{clinvar_pathogenic_for_this_gene} if (exists $hvariation->{value}->{clinvar_pathogenic_for_this_gene});
	if ($hvariation->{value}->{clinvar_pathogenic}) {
		my $chr = $project->getChromosome($hvariation->{value}->{chromosome});
		my $clinvar_id = $hvariation->{value}->{clinvar_id};
		my $g = $project->newGene($gene->{id});
		if ($chr->is_clinvar_pathogenic_for_gene($clinvar_id, $g)) {
			$hvariation->{value}->{clinvar_pathogenic_for_this_gene} = 1;
			return 1;
		}
		else {
			$hvariation->{value}->{clinvar_pathogenic_for_this_gene} = undef;
			$hvariation->{value}->{clinvar_pathogenic} = undef;
			$hvariation->{value}->{clinvar} = '';
			$hvariation->{html}->{clinvar} = '';
		}
	}
	return;
}

sub table_validation_without_local {
	my ($project, $hvariation, $gene) = @_;
	my $cgi = new CGI();
	my $color = "#555";
	check_is_hgmd_dm_for_gene($project, $hvariation, $gene);
	check_is_clinvar_pathogenic_for_gene($project, $hvariation, $gene);
	if ($hvariation->{value}->{dm} or $hvariation->{value}->{clinvar_pathogenic}){
		$color = "red";
	}
	elsif  ($hvariation->{value}->{hgmd_id} or $hvariation->{value}->{clinvar_id}) {
		$color = "orange";
	}
	my $html= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"max-width:300px;box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
	$html .= $cgi->start_Tr().$cgi->th(["HGMD","Clinvar"]).$cgi->end_Tr();
	$html .= $cgi->start_Tr();
		$hvariation->{html}->{table_validation} = $html;
		$hvariation->{html}->{table_validation} .= $cgi->td([$hvariation->{html}->{hgmd},$hvariation->{html}->{clinvar}]);
		
		my $hgmd_no_access = qq{<span class="glyphicon glyphicon-ban-circle" aria-hidden="true" style='font-size:12px;color:black;'></span>};
		$hvariation->{html}->{table_validation_hgmd_no_access} = $html;
		$hvariation->{html}->{table_validation_hgmd_no_access} .= $cgi->td([$hgmd_no_access,$hvariation->{html}->{clinvar}]);
		
		my $v_phen = ' ';
   		if (exists $hvariation->{hgmd_phenotype}) {
   			$v_phen = $hvariation->{hgmd_phenotype};
   			$v_phen =~ s/"//g;
			$hvariation->{html}->{table_validation} .= qq{<tr><td  colspan='3' style='color:#4CAF50;'>$v_phen</td></tr>};
			$hvariation->{html}->{table_validation_hgmd_no_access} .= qq{<tr><td  colspan='3' style='color:#4CAF50;'>$v_phen</td></tr>};
   		}
		$hvariation->{html}->{table_validation} .= $cgi->end_Tr().$cgi->end_table;
		$hvariation->{html}->{table_validation_hgmd_no_access} .= $cgi->end_Tr().$cgi->end_table;
		return $hvariation->{html}->{table_validation};
	
}

sub table_validation {
		my ($patient,$hvariation,$gene) = @_;
		my $cgi = new CGI();
		my $color = "#555";
		#check_is_hgmd_dm_for_gene($patient->getProject(), $hvariation, $gene);
		check_is_clinvar_pathogenic_for_gene($patient->getProject(), $hvariation, $gene);
		if ($hvariation->{value}->{dm} or $hvariation->{value}->{clinvar_pathogenic}){
			$color = "red";
		}
		elsif  ($hvariation->{value}->{hgmd_id} or $hvariation->{value}->{clinvar_id}) {
			$color = "orange";
		}
		my $html= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"max-width:300px;box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
		$html .= $cgi->start_Tr().$cgi->th(["HGMD","Clinvar","Local"]).$cgi->end_Tr();
		$html .= $cgi->start_Tr();
		my $all_validations = $patient->project->validations;
		my $val_id = $gene->{id}."!".$hvariation->{value}->{id};
		my $local_validation = $patient->project->getValidationVariation($val_id,$patient);
		
		my $polyweb = $minus ; 
	
		if ($local_validation){
				my $saved = $local_validation->{validation};
				#$hvariatio printButton(4,[3,4],$v->hgmd->{class},qq{onClick="$cmd"}); 
				my $id = $hvariation->{value}->{id};
				my $gene_id = $gene->{id};
				my $cmd = qq{view_variation_validation(\'$id\', \'$gene_id\')};
				$polyweb = printButton($saved,[4,5],$patient->buffer->value_validation->{$saved},qq{onClick="$cmd"}) ;
		}
		$hvariation->{html}->{table_validation} = $html;
		$hvariation->{html}->{table_validation} .= $cgi->td([$hvariation->{html}->{hgmd},$hvariation->{html}->{clinvar},$polyweb]);
		
		my $hgmd_no_access = qq{<span class="glyphicon glyphicon-ban-circle" aria-hidden="true" style='font-size:12px;color:black;'></span>};
		$hvariation->{html}->{table_validation_hgmd_no_access} = $html;
		$hvariation->{html}->{table_validation_hgmd_no_access} .= $cgi->td([$hgmd_no_access,$hvariation->{html}->{clinvar},$polyweb]);
		
		my $v_phen = ' ';
   		if (exists $hvariation->{hgmd_phenotype}) {
   			$v_phen = $hvariation->{hgmd_phenotype};
   			$v_phen =~ s/"//g;
			$hvariation->{html}->{table_validation} .= qq{<tr><td  colspan='3' style='color:#4CAF50;'>$v_phen</td></tr>};
			$hvariation->{html}->{table_validation_hgmd_no_access} .= qq{<tr><td  colspan='3' style='color:#4CAF50;'>$v_phen</td></tr>};
   		}
		$hvariation->{html}->{table_validation} .= $cgi->end_Tr().$cgi->end_table;
		$hvariation->{html}->{table_validation_hgmd_no_access} .= $cgi->end_Tr().$cgi->end_table;
		return $hvariation->{html}->{table_validation};
		
}

sub vclinvar {
		my ($v,$hvariation) = @_;
		 my $cl  = "" ;
		 my $alert = 0;
	 my $v1 = $v->score_clinvar();
	  value_html($hvariation,"clinvar","-");
	  $hvariation->{value}->{clinvar_pathogenic} = 1 if  $v->isClinvarPathogenic ;
	  
	   #  $hvariation->{clinvar_alert}  = 0 ;
	 if ($v1){
	 	 $hvariation->{value}->{clinvar_id} = $v->clinvar->{id};
	 	my $uc = qq{https://www.ncbi.nlm.nih.gov/clinvar/?term=}.$v->clinvar->{id}."[alleleid]";
	 	my $a = qq{<a href="$uc" target="_blank" style="color:white">}.$v->text_clinvar()."</a>"; 
	 	my $oc = qq{onClick='window.open("$uc")'};
	 	value_html($hvariation,"clinvar",$v->text_clinvar(),  printButton($v,[1,4],$v->text_clinvar(),$oc));
	 
	 	if (($v == 4 || $v==5)    ){
	 		$alert  = 4  if $hvariation->{value}->{freq_level} <= 2 ;
	 		$alert ++;
	 	}
	 }
	 	 value_html($hvariation,"clinvar_alert",$alert,$alert);
}

sub vhgmd {
		my ($v,$hvariation) = @_;
		 
		 value_html($hvariation,"hgmd",$minus);
	 if ($v->hgmd_id()){
	 	value_html($hvariation,"hgmd",1,1);
	 	if  ($v->isDM){
	 		my $dm =1;
	 		value_html($hvariation,"dm",1,1);
			 		
	 	}
	 	$hvariation->{value}->{hgmd_id} = $v->hgmd_id;
	 	my $txt = $v->hgmd->{phen}." - ".$v->hgmd->{class}." - ".$v->hgmd_id;
	 	$txt =~s/\"//g;
	 	my $n1 = $v->project->name;
	 	my $n2 = $v->hgmd_id;
	 	my $n3 = $v->id;
	 	my $cmd;
	 	if (exists $hvariation->{html}->{no_css_polydiag}) { $cmd = qq{zoomHgmdWithoutCss(\'$n1\',\'$n2\',\'$n3\')}; }
	 	else { $cmd = qq{zoomHgmd(\'$n1\',\'$n2\',\'$n3\')}; }
	 	my $nb = 1;
		if ($v->hgmd->{class} eq 'DM') { $nb = 4; }
		elsif ($v->hgmd->{class} eq 'DM?') { $nb = 3; }
		elsif ($v->hgmd->{class}) { $nb = 2; }
		$hvariation->{hgmd}    = printButton($nb,[3,4],$v->hgmd->{class},qq{onClick="$cmd"}); 
		value_html($hvariation,"hgmd",$v->hgmd->{class}.":".$v->hgmd->{id},printButton($nb,[3,4],$v->hgmd->{class},qq{onClick="$cmd"}));
	 }
	 
}


sub vfreq_level {
	my ($v,$hvariation) = @_;
		my $scaled =  $v->scaled_score_frequence(); 
		 value_html($hvariation,"freq_level",$scaled->{freq_level},$scaled->{freq_level});
		  value_html($hvariation,"freq_score",$scaled->{freq_score},$scaled->{freq_score});
}

my $server =  $ENV{HTTP_HOST};
my $variation_script = $ENV{SCRIPT_NAME};
$variation_script =~s/patient_/variation_/;
$server = "darwin.bipd.fr" if $server eq "bipd";
$server = "www.polyweb.fr" if $server =~/10\.200\.27/;
my $deja_vu_url = "http://$server//polyweb/polydejavu/dejavu.html?input=";

sub abutton {
	my ($href,$value) = @_;
	return qq{<a class="btn btn-xs btn-primary" href="$href" target="_blank" style="background-color: #D0D0D0;font-size: 7px;font-family:  Verdana;color:black" role="button">}.$value."</a>";
}
sub obutton {
	my ($onclick,$value) = @_;
	return qq{<a class="btn btn-xs btn-primary" onclick=\'$onclick\' target="_blank" style="background-color: #D0D0D0;font-size: 7px;font-family:  Verdana;color:black" role="button">}.$value."</a>";
}

sub vdejavu {
	my ($v,$hvariation) = @_;
		#delete $v->{dejaVuInfosForDiag};
		my $project = $v->project;
		$hvariation->{value}->{other_project} = $v->other_projects();
		$hvariation->{value}->{other_patients} = $v->other_patients();
		$hvariation->{value}->{other_patients_ho} = $v->other_patients_ho();
		$hvariation->{value}->{similar_projects} = $v->similar_projects();
		$hvariation->{value}->{similar_patients} = $v->similar_patients();
		$hvariation->{value}->{similar_patients_ho} = $v->similar_patients_ho();
		$hvariation->{value}->{this_run_project} = '-';
		$hvariation->{value}->{this_run_patients} = $v->in_this_run_patients()."/".scalar(@{$project->getPatients});
		#$v->{value}->{this_run__patients_ho} = $v->similar_patients_ho();
		value_html($hvariation,"deja_vu",$v->other_projects(),table_dejavu($v));
	
}

sub vcnv {
	my ($v, $hvariation, $patient) = @_;
	$hvariation->{value}->{manta}->{is_imprecise} = $v->is_imprecise();
	my $chr_id = $v->getChromosome->id();
	my $start = $v->start();
	my $end = $v->start + $v->length() + 1;
	foreach my $p (@{$patient->getFamily->getMembers()}) {
		$hvariation->{value}->{dp_details}->{mean}->{$p->name()} = $p->meanDepth($chr_id, $start, $end);
		$hvariation->{value}->{dp_details}->{min}->{$p->name()} = $p->minDepth($chr_id, $start, $end);
		$hvariation->{value}->{dp_details}->{max}->{$p->name()} = $p->maxDepth($chr_id, $start, $end);
		$hvariation->{value}->{dp_details}->{norm}->{$p->name()} = $p->cnv_region_ratio_norm($chr_id, $start, $end);
	}
	if ($v->isCnv()) {
		$hvariation->{value}->{is_cnv} = 1;
		#$hvariation->{value}->{sd_value_controls} = sprintf("%.2f",$patient->sd_value_dude($chr_id,$start,$end));
		$hvariation->{value}->{cnv_details_genes} = $v->get_genes_transcripts_details_dup_del();
		if ($v->getProject->isGenome()) {
			$hvariation->{value}->{cnv_confidence}->{$patient->name()} = $v->cnv_confidence($patient);
		}
	}
	
	#TODO: a pouvoir ajouter
	#ajout sd_value_dude de patient (= les controls) OK
	#ajout sd_value_patient de patien (= le pat) OK
	# enlever transcript + ajout des regions introns / exons de perdus OK
	# modifier le tableau trio avc les champs de MANTA PE, SR, etc... + score DUDE OK
	# voir clinSV -> results_test_data -> XLS -> tag HIGH / PASS / LOW OK
	# changer header violet en autre couleur si CNV ? OK
	
}

sub correct_table_dejavu {
		my ($v) = @_;
		
		#problem avec l'url du dejavu 
		my $z = $v->{html}->{deja_vu};
		my $x1 = 'href="http:////polyweb';
		my $server =  $ENV{HTTP_HOST};
		my $variation_script = $ENV{SCRIPT_NAME};
		$variation_script =~s/patient_/variation_/;
		$server = "darwin.bipd.fr" if $server eq "bipd";
		$server = "www.polyweb.fr" if $server =~/10\.200\.27/;
	
		my $x2 = "href=\"http://$server/polyweb";
		my $deja_vu_url = "http://$server//polyweb/polydejavu/dejavu.html?input=";
		$v->{html}->{deja_vu} =~ s/$x1/$x2/g;
		
}

sub table_dejavu_live {
	my ($v,$project,$patient,$g) = @_;
	#TODO: pour les dejavu genome interoger la base db dejavu pour mettre a jour les valeurs +  IGV dejavu
	if ($patient){
		my $vid1= $v->{id};
		my $onclick = qq{goto_dejavu("$vid1")};
		my $html =  $v->{html}->{deja_vu};
		my $name = $patient->name();
	 	$v->{html}->{deja_vu} =~ s/goto_dejavu\("$vid1"\)/goto_dejavu\("$vid1","$name"\)/g;
	}
	return if $project->isGenome();
	my $vobj = $project->_newVariant($v->{id});
	#die();
	my $cgi          = new CGI();
	my $html =  $v->{html}->{deja_vu};
	my $vobj = $project->_newVariant($v->{id});
		my $ph ="?";
		$ph =  $project->getPhenotypes->[0]->short_name() if @{$project->getPhenotypes};
		$html =~s/ID/$ph/;
		$v->{value}->{this_run_patients} = $vobj->in_this_run_patients()."/".scalar(@{$project->getPatients});
		my $vid1= $v->{id};
		my $name = $patient->name();
		my $gene = $project->newGene($g->{id});
		my @lTr = @{$gene->getMainTranscripts()};
		if (scalar(@lTr) > 0) {
			my $tr_id = $lTr[0]->id();
			my $proj_name = $project->name();
			my $text = qq{<tr onClick="goto_dejavu_in_this_run('$proj_name','$vid1','$tr_id')">};
			#my $text = qq{<tr onClick="goto_dejavu_in_this_run('$vid1','$name')">};
			$text .= $cgi->td("Run").$cgi->td({colspan=>3},obutton("",$v->{value}->{this_run_patients})).$cgi->end_Tr();
			$html =~s/<\/table>/$text <\/table>/;
		}
	$v->{html}->{deja_vu} =  $html;
	return ;
}

sub table_dejavu {
	my ($v, $no_phenotype) = @_;
	my $cgi          = new CGI();
	my $color = "#555";
		
	my $po = $v->other_projects();
	if ($po <5){
		$color = "red";
	}
	elsif ($po < 20){
		$color = "orange";
	} 
	elsif ($po < 50){
		$color = "blue";
	} 
#	my $deja_vu_url = "http://$server//polyweb/polydejavu/dejavu.html?input=";
	my $html= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
	
	$html.= $cgi->start_Tr();
		$html.= $cgi->th("");
	$html.= $cgi->th("Pr");
	$html.= $cgi->th("Sa");
	$html.= $cgi->th("Ho");
	$html.= $cgi->end_Tr();
	$html.= $cgi->start_Tr();
	$html.= $cgi->td("other");
	my $href = $deja_vu_url.$v->id;
	my $vid1= $v->id;
	my $onclick = qq{goto_dejavu("$vid1")};
	$html.= $cgi->td(obutton($onclick,$po));
	$html.= $cgi->td(obutton($onclick,$v->other_patients));
	$html.= $cgi->td(obutton($onclick,$v->other_patients_ho));
	$html.= $cgi->end_Tr();
	if (not $no_phenotype) {
		$html.= $cgi->start_Tr();
		my $ps = $v->similar_projects();
		my $ph;
		$ph =  $v->project->getPhenotypes->[0]->short_name() if @{$v->project->getPhenotypes};
		$html.= $cgi->td("$ph");
		$html.= $cgi->td(obutton($onclick,$ps));
		$html.= $cgi->td(obutton($onclick,$v->similar_patients));
		$html.= $cgi->td(obutton($onclick,$v->similar_patients_ho));
		$html.= $cgi->end_Tr();
	}
	$html.=$cgi->end_table();
	return $html;
}




################
# construct hash variant
#################

my @header_transcripts = ("consequence","enst","nm","ccds","appris","exon","nomenclature","codons","codons_AA", "polyphen","sift","cadd","revel","dbscsnv",'spliceAI');
sub construct_hash_variant_patient {
	my ($project,$v,$patient,$hvariation) = @_;
	my $debug;
	valamut_igv($v,$hvariation,$patient,$debug);
	vcnv($v, $hvariation, $patient) if ($v->isCnv());
	trio($v,$hvariation,$patient);
	vsequencing($v,$hvariation,$patient);
}

sub construct_hash_variant_global {
	my ($project,$v,$vquery,$debug) = @_;
	my $hvariation;
	$hvariation->{value}->{id} =  $v->id;
	$hvariation->{html}->{id} =  $v->id;
	
	
	$hvariation->{value}->{type} = $v->type;
	$hvariation->{html}->{type} = $v->type;
	foreach my $p (@{$v->getPatients}){
		$hvariation->{html}->{patients} .= $p->name." ";
		
		
	}
	
	$hvariation->{html}->{infos} .= $v->{infos};
	my $vn=$v->vcf_id;
	$vn =~ s/_/-/g;
	$vn=~ s/chr//;
	$hvariation->{value}->{gnomad_id} = $vn;
	$hvariation->{html}->{gnomad_id} = $vn;	
	$hvariation->{value}->{is_cnv} = 0;	
	vcosmic($v,$hvariation)	if ($project->isSomatic);
	vgnomad($v,$hvariation);
	vname($v,$hvariation);
	vspliceAI($v,$hvariation);
	vvarsome($hvariation,$debug);
	vdivers($v,$hvariation);

	vclinvar($v,$hvariation);
	vhgmd($v,$hvariation);
	vdejavu($v,$hvariation);

	

	my $gs;
	my $cgi          = new CGI();
	foreach my $g (@{$v->getGenes}){
		push(@$gs,$g->id);
		$hvariation->{genes}->{$g->id}= construct_hash_transcript($v,$cgi,\@header_transcripts,2,$g);
		#$hvariation->{html}->{genes}->{$g->id} = construct_table_transcript($v,$cgi,\@header_transcripts,2,$g);
	}

	return $hvariation;
	
}

sub href {
	my ($url,$name) = @_;
	$url.=$name;
	return qq{<a href="$url" target="_blank">}.$name."</a>";
}

sub value_html_badge {
	my($h,$type,$value) = @_;
	value_html($h,$type,$value,printBadge($value));
}

################
# construct table transcripts
#################

sub construct_hash_transcript {
	my ($v,$cgi,$header_transcripts,$level,$gene) = @_;
	#confess() if $v->start ==45745275;
	my $debug ;
	$debug = 1  if $v->id eq "1_8865904_CTTTTTTTTTTT_C";
	my  $transcripts = $v->getTranscripts();
	my $all_transcripts;
	my $max_ai = $v->max_spliceAI_score($gene);
	my $max_cat = $v->max_spliceAI_categorie($gene);
	my ($refseq_hgmd_1, $refseq_hgmd, $tmprefseq);
	if ($v->hgmd_details) {
		$refseq_hgmd_1 = $v->hgmd_details->{refseq};
		($refseq_hgmd, $tmprefseq) = split('\.', $refseq_hgmd_1);
	}
	foreach my $tr1 (sort { ($himpact_sorted->{$v->effectImpact($b)} <=>  $himpact_sorted->{$v->effectImpact($a)}) or ($a->appris_level <=> $b->appris_level)} @$transcripts) {
		
		next if $tr1->getGene->id ne $gene->id;
		my $htr = {};
		# score
		my $cadd = $v->cadd_score;
		
		value_html($htr,"uniq_id",$tr1->id."-".$v->id);
		value_html($htr,"cadd",$cadd,printBadge($cadd,[20,30]));
		my $rf = $v->dbscsnv_rf();
		my $ada = $v->dbscsnv_rf();
		value_html($htr,"dbscsnv",$rf.":".$ada,printBadge($v->dbscsnv_rf,[0.6,0.9]).printBadge($v->dbscsnv_ada,[0.6,0.9]));
		
		#revel 
		
		#TODO: a corriger
		my $revel = 0;
		#my $revel =  $v->revel_score;
		value_html($htr,"revel",$revel,printBadge($revel,[0.5,0.9]));
		
		
		
		#id , name synonym
			#name
			my $enst = $tr1->name;
		$htr->{value}->{name} = $tr1->name();
		value_html($htr,"name",$enst);
		#enst 
		
		my $u = qq{https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=}.$tr1->name;
		
		value_html($htr,"enst",$enst,href(qq{https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=},$enst));
		value_html($htr,"trid",$tr1->id);
		
		#nm
		
		#TODO: a corriger
		my $nm = '-';
		#my $nm =$tr1->refseq_names;
		my $is_same_nm_as_hgmd;
		if ($refseq_hgmd) {
			$is_same_nm_as_hgmd = 1 if $nm =~ /$refseq_hgmd/;
			if ($is_same_nm_as_hgmd) { value_html($htr,"nm",$htr,printBlueSimpleBadge($nm)); }
			else { value_html_badge($htr,"nm",$nm); }
		}
		else { value_html_badge($htr,"nm",$nm); }
		
		my $ccds = $tr1->ccds_name;
		value_html($htr,"ccds",$ccds,href(qq{https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=},$ccds));
		
		value_html_badge($htr,"appris",$tr1->appris_type);
		
		
		
		#init coding info 
		my @coding_infos = ("sift","polyphen","prot","codons","codons_AA");
		
		foreach my $val (@coding_infos) {
			value_html($htr,$val,"-");
		}
			
		#Protein Information
		if ($v->isCoding($tr1)){
			my $prot = $tr1->getProtein();
			#my $cds = $v->getOrfPosition($prot);
			if ($v->isLargeDeletion or $v->isLargeDuplication){
			value_html_badge($htr,"prot","-");
			value_html_badge($htr,"codons","-");
			value_html_badge($htr,"codons_AA","-");
			}
			else {
			value_html_badge($htr,"prot",$v->getProteinPosition($prot));
			value_html_badge($htr,"codons",$v->getCodons($tr1));
	#		warn "****".$v->protein_nomenclature($prot) if $v->start ==45745275;
			value_html_badge($htr,"codons_AA",$v->protein_nomenclature($prot));
			}
			#sift ;
			my $sift = $v->siftScore($tr1->getProtein);
			value_html($htr,"sift",$sift,printInvBadge($sift,[0,1,0,05]));
			
			my$polyphen = $v->polyphenScore($tr1->getProtein);
			value_html($htr,"polyphen",$polyphen,printBadge($polyphen,[0.446,0.908]));
		}
		
		value_html($htr,"spliceAI",$max_ai,printButton($max_ai,[0.5, 0.9],$max_ai,undef,$max_ai."/".$max_cat) );
		
		my $alphamissense = $v->alphamissense($tr1);
		value_html($htr,"alphamissense",$alphamissense,printBadge($alphamissense,[0.34, 0.564]));
		
	#	my $text_alert = 'SpliceAI values - '.join(', ', @l_score_spliceAI);
			#$htr->{html}->{spliceAI} = printButton($max_ai,[0.5, 0.9],$max_cat,undef,$max_cat);
		#exons information
		my $te = $tr1->findExonNumber($v->start, $v->end);
		 if ($te == -1){
		 	my $tc = $tr1->findNearestExon($v->start, $v->end);
		 	value_html_badge($htr,"exon",$tc);
		 }
		 else {
		 value_html_badge($htr,"exon",'exon'.$te);
		 }
		#nomenclature 
			if ($v->isLargeDeletion or $v->isLargeDuplication){
					$htr->{nomenclature} =  "-";
					value_html_badge($htr,"nomenclature","-");
			}
			else {
				$htr->{nomenclature} =  $v->getNomenclature($tr1);
				#warn "******".$htr->{nomenclature}." ".$v->getNomenclature($tr1) if $v->start ==45745275;
				#value_html_badge($htr,"nomenclature","XX");
				value_html_badge($htr,"nomenclature",$htr->{nomenclature});
			}
	
		value_html($htr,"impact_score_text",$v->effectImpact($tr1));
		value_html($htr,"impact_score",$himpact_sorted->{$v->effectImpact($tr1)});
		
		
		$htr->{impact_score}  =  $himpact_sorted->{$v->effectImpact($tr1)};
		my $cons = $v->variationTypeInterface($tr1);
		my $x = $himpact_sorted->{$v->effectImpact($tr1)};
		value_html($htr,"consequence",$cons,printButton($x,[3,4],$cons));
		my $main =  0 ;
		 $main  = 1 if $tr1->isMain();
		value_html($htr,"main",$main,$main);
#		delete $htr->{html};
		push(@$all_transcripts,$htr)
	#	if ($is_same_nm_as_hgmd){
	#		unshift(@$all_transcripts,$htr);
	#		next;
	#	} 
	#	push(@$all_transcripts,$htr);
	}#end for transcript
	return $all_transcripts;
}

sub table_transcripts_cnv {
	my ($v, $atr, $gene_id, $header_transcripts_cnv) =@_;
#	$h->{$g->id()}->{$t_id}->{exons_introns}->{$e->start()}->{$e_id}->{overlap} = $overlap.' nt ('.$perc.'%)';
	my $color_header = "white";
	if ($v->{value}->{type} eq 'large_deletion') {
		$color_header = "#fd8253";
	}
	elsif ($v->{value}->{type} eq 'large_duplication') {
		$color_header = "#5393fd";
	}
	my $cgi = new CGI();
	my $html;
	$html .= qq{<div style="max-height:250px;overflow-y:auto;border:solid 1px grey;">};
	$html .= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color_header;font-size: 7px;font-family: Verdana;margin-bottom:0px;max-height:150px;"});
	$html .=  $cgi->start_Tr({style=>"position:sticky;z-index:8;top:0;background-color:$color_header;color:white;"}).$cgi->th({style=>"text-align: center;"}, $header_transcripts_cnv)."".$cgi->end_Tr();
	
	#$hvariation->{value}->{cnv_details_genes}
	
	foreach my $htr (@$atr){
		my $t_id = $htr->{html}->{trid};
		#next unless (exists $v->{value}->{cnv_details_genes}->{$gene_id}->{$t_id});
		my @l_pos_e_i = sort {$a <=> $b} keys %{$v->{value}->{cnv_details_genes}->{$gene_id}->{$t_id}->{positions}};
		my $first = $v->{value}->{cnv_details_genes}->{$gene_id}->{$t_id}->{positions}->{$l_pos_e_i[0]};
		my $last;
		if (scalar(@l_pos_e_i) > 1) {
			$last = $v->{value}->{cnv_details_genes}->{$gene_id}->{$t_id}->{positions}->{$l_pos_e_i[-1]};
		}
		$html .=  $cgi->start_Tr();
		$html.= "<td>".$htr->{html}->{consequence}."</td>";
		$html.= "<td>".$htr->{html}->{enst}."</td>";
		$html.= "<td>".$htr->{html}->{nm}."</td>";
		$html.= "<td>".$htr->{html}->{ccds}."</td>";
		$html.= "<td>".$htr->{html}->{appris}."</td>";
		if ($last) {
			$html.= "<td><b>".$first."</b></td>";
			$html.= "<td><b>".$last."</b></td>";
		}
		else {
			$html.= "<td colspan='2'><b>".$first."</b></td>";
		}
		
		$html .=  $cgi->end_Tr();
	}
	$html.=$cgi->end_table();
	return $html;
}

sub table_transcripts {
	my ($atr,$header_transcripts,$compact_table) =@_;
	my $level = 2;
	my $cgi          = new CGI();
	my $value = $atr->[0]->{value}->{impact_score};
	my $color = "#555";
	if ($value > 3){
		$color = "red";
	}
	if ($value >= 2){
		$color = "#FF8800";
	}
	
	my $html;
	$html .= qq{<div style="max-height:250px;overflow-y:auto;border:solid 1px grey;">};
	
	$html .= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family: Verdana;margin-bottom:0px;max-height:150px;"});
	my @colors = ("#F9F6FF","#F9F9F9");
	
	my $nb =0;
	my $nb_skip = 0;
	 
	$html .=  $cgi->start_Tr({style=>"position:sticky;z-index:8;top:0;background-color:grey;color:white;"}).$cgi->th({style=>"text-align: center;"}, $header_transcripts)."".$cgi->end_Tr();
	 
	my $rids = [];
	my $hhide;
	foreach my $htr (@$atr){
			my $zid = $htr->{value}->{uniq_id} ;
			my $hide;
			
			#{style=>"border: 1px solid black;color:lightgrey;$style"}).$cgi->td([$mother->name,qq{<i class="fa fa-female fa-2x" aria-hidden="true" style="color:ligthgrey"></i>},$type,$ps."%",""]);
		 	#$htr->{html}->{spliceAI} = "-" ;
		 	my $hide_id= $htr->{value}->{ccds};
		 	$hide_id = $htr->{value}->{nomenclature} if $hide_id eq "";
		 	
		 	$hide = "display:none;"  if ($htr->{value}->{impact_score}  <  $level and $nb >0) or (exists $hhide->{$hide_id})  ;
		 	$hide = "display:none;"  if ($compact_table and not $htr->{value}->{ccds} and $nb != $nb_skip)  ;
		 	$hide = "display:none;"  if ($compact_table and $nb - $nb_skip >= 3)  ;
		 	$hhide->{$hide_id} ++ ;
		 	$hhide->{$htr->{value}->{nomenclature}} ++ ;
			$nb_skip ++ if $hide;
			my $c = $colors[$nb%2];
			$nb ++;
			my $rid = "row_".$htr->{value}->{trid}."_".$zid;
			push(@$rids,$rid) if $hide;
			
		 $html .=  $cgi->start_Tr({id=>$rid,style=>"border: 1px solid;background-color:$c ;".$hide});
		 foreach my $c (@$header_transcripts){
		 	$html.= $cgi->td("<center>".$htr->{html}->{$c}."</center>");
		 }
		 $html.= $cgi->end_Tr();
	}
	 if ($nb_skip ==0){
		$html.=$cgi->end_table();
		return $html
	 }
	my $js = encode_json $rids;
	my $za = "hide_tr_".time."_".int(rand(50000));
	$html .=  $cgi->start_Tr({id=>$za});

	$html.= $cgi->td({style=>"box-shadow: 1px 1px 2px #555;background-color:#CECFCE;",colspan=>scalar(@$header_transcripts),onClick=>qq{showTranscripts($js,"$za");}},qq{<span class="glyphicon glyphicon-plus"></span> }."view $nb_skip Transcripts");
	$html.= $cgi->end_Tr();
	
	$html.=$cgi->end_table();
	$html.=qq{</div>};
	#$atr->{html} = $html;
	return $html;
	
}
sub color_model {
	my ($model) = @_;
	my $color = "#555";
	if ($model eq 'mother'){
		$color = "#F7E4E4";
	}
	elsif ($model eq 'mother_c'){
		$color = "#779ECB";
	}
	elsif ($model eq 'father'){
		$color = "#A9BCD1";
	}
	elsif ($model eq 'father_c'){
		$color = "#779ECB";
	}
	elsif ($model =~ /denovo/){
		$color = "#e74c3c";
	}
	elsif ($model =~ /rece/){
		$color = "#EE82EE";
	}
	elsif ($model =~ /mosa/){
		$color = "#F9885C";
		$color = "#FDE6B0";
	}
	elsif ($model =~ /uni/){
		$color = "#F9885C";
		$color = "#45B8AC";
	}
	return $color;
}

sub text_model {
	my ($model) = @_;
	my $text = "-";
	if (lc($model) =~ /mother/){
		$text = "<img src='https://img.icons8.com/offices/24/000000/guest-female.png'>";
	}
	elsif ($model eq 'mother_c'){
		$text = "<img src='/icons/Polyicons/female-d.png'>";
	}
	elsif ($model =~ /father/){
	$text = "<img src='https://img.icons8.com/offices/24/000000/person-male.png'>";
	}
	elsif ($model eq 'father_c'){
		$text = "<img src='/icons/Polyicons/male-d.png'>";
	}
	elsif ($model =~ /denovo/){
		 $text = '<button class="btn btn-primary btn-xs" type="button" style="background-color:#e74c3c;color:white">'.$model.'</button>';
	}
	elsif ($model =~ /rece/){
		 $text = '<button class="btn btn-primary btn-xs" type="button" style="background-color:#EE82EE;color:white">'.$model.'</button>';
	}
	elsif ($model =~ /mosa/){
		$text = "mosaic";
	}
	elsif ($model =~ /uni/){
		$text = $model;
	}
	return $text;
}
sub solo {
	 my ($var,$hvariation,$patient) = @_;
	 	my $cgi          = new CGI();
	my $color = "green";
	
	my $value  = $var->getPourcentAllele($patient);
	if ($value< 10){
		$color = "red";
	}
	elsif ($value < 25){
		$color = "orange";
	}
	

	 my $html =$cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
	 my $type =$minus;
		$type = "ho" if $var->isHomozygote($patient);
		$type = "he" if $var->isHeterozygote($patient);
		my $depth = $var->getDepth($patient);
		my $icon = $male;
		$icon = $female if $patient->isFemale;
	 	my $tab = [$patient->name,$icon,$type,$var->getPourcentAllele($patient)."%","$depth"]; 	
	 	$html .=  $cgi->start_Tr({style=>"border: 1px solid black;color:black;"}).$cgi->td($tab).$cgi->end_Tr();
	 	$html .=$cgi->end_table
}

sub trio {
	 my ($var,$hvariation,$patient,$project) = @_;
#die();
	my $debug;
	my $cgi          = new CGI();
	 return $hvariation->{triom} if exists $hvariation->{triom} && exists $hvariation->{trio};
#	 
	 
	 $hvariation->{triom} ++;
	  $project = $var->project unless $project;
	  unless ($var) {
	  	$var = $project->_newVariant($hvariation->{id}); 
	  }
	  $debug = 1 if $hvariation->{id} eq "5_254512_G_A";
	 my $fam = $patient->getFamily();
	#die();
	my $patient_input_name = $patient->name();
	my $var_id = $var->id();
	 unless ($fam->isTrio){
	 	 $hvariation->{html}->{trio} = solo($var,$hvariation,$patient);
	 	 return;
	 }
	# ENST00000301067 
	  return if exists $hvariation->{trio};
	  $hvariation->{trio} = 0 unless $project->isFamilial();
	  
	  return $hvariation->{trio} unless $project->isFamilial();

	   
	  my $children;
	  my $father;
	  my $mother;
	  my $all_tab;
	
	  my $fam = $patient->getFamily();
	  if ($patient->isChild){
	  	$children = $fam->getChildren();
	  	$mother = $fam->getMother;
	  	$father = $fam->getFather;
	  }
	  else {
	  		$children = [$patient];
	  	$mother = undef;
	  	$father = undef;
	  }
#	  	unless ($father ){
#	  		$hvariation->{trio} ="-";
#	  		return ;
#	  	7
#	  	unless ($mother){
#	  		$hvariation->{trio} ="-";
#	  		return ;
#	  	}
	my $gene_used;
	if (exists $hvariation->{scaled_score_gene_used}) {
		foreach my $this_g (@{$var->getGenes()}) {
			if (($this_g->id() eq $hvariation->{scaled_score_gene_used}) or ($this_g->external_name() eq $hvariation->{scaled_score_gene_used})) {
				$gene_used = $this_g;
				last;
			}
		}
	}
	my $model =  lc ( $var->getTransmissionModelType($fam,$patient,$gene_used));
	#if ($model eq "denovo" && $var->isLargeDeletion){
		#14-105417183-cnv-del-1017 22-45745275-del-926
#	if ($var->start == 72088484){
#		warn "***". $patient->name();
#		warn $model." ".$var->getChromosome().":".$var->start."-".($var->start+$var->length)." =+>". ($var->other_projects + $var->similar_projects());
#		$var->score_refined($patient,1);
#		warn $var;
#		confess();
#	}
	
	my $color = "#555";
	$color = color_model($model);
		
	  	my $html;
		if ($var->isCnv) {
			my $cnv_confidence = $hvariation->{value}->{cnv_confidence}->{$patient->name()};
			my $box_color = "grey";
			if ($cnv_confidence eq 'high') { $box_color = "red"; }
			elsif ($cnv_confidence eq 'medium') { $box_color = "orange"; }
			elsif ($cnv_confidence eq 'low') { $box_color = "blue"; }
			$html = $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $box_color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
			my $pr_text = qq{<button onclick="alert('PR: Number of spanning read pairs which strongly (Q30) support the REF or ALT alleles')">PR</button>};
			my $sr_text = qq{<button onclick="alert('SR: Number of split-reads which strongly (Q30) support the REF or ALT alleles')">SR</button>};
			$html .=  $cgi->start_Tr({style=>"border: 1px solid black;color:black;text-align:center;font-weight: bold;"}).$cgi->td('pat').$cgi->td('sex/status').$cgi->td('he/ho').$cgi->td($pr_text).$cgi->td($sr_text).$cgi->td('norm dp').$cgi->td('cnv score').$cgi->td('trans').$cgi->end_Tr();
		}
		else {
			$html = $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
		}
		my $html_line;
		my $type = "-";
		my $sex_status;
		my $project_name = $project->name();
		my $gene_tr_name;
		if (exists $hvariation->{scaled_score_gene_used} and exists $hvariation->{scaled_score_transcript_used}) {
			$gene_tr_name = $hvariation->{scaled_score_gene_used}.':'.$hvariation->{scaled_score_transcript_used};
		}
		else {
			 $gene_tr_name = '';
			 $gene_tr_name = $var->getGenes->[0]->id() if scalar(@{$var->getGenes});
			  }
		
		my $chr_id = $var->getChromosome->id();
		my $start = $var->start();
		my $end = $var->start + $var->length();
		if ($mother){
			my $mean_dp_m = int($mother->meanDepth($chr_id, $start, $end));
			my $text_dp_m = -$mean_dp_m;
			my $norm_dp_m = $mother->cnv_region_ratio_norm($chr_id, $start, $end);
			my $text_norm_dp_m = int($norm_dp_m);
			if ($mother->isIll()) { $sex_status = "<center><img src='/icons/Polyicons/female-d.png'></center>"; }
			else { $sex_status = "<center><img src='/icons/Polyicons/female-s.png'></center>"; }
			my $patient_name = $mother->name();
			my $name = qq{<button style='color:black;' onClick="view_var_from_proj_gene_pat('$project_name', '$gene_tr_name', '$patient_input_name', '$var_id', 'mother')">$patient_name</button>};
			if (lc($model) eq 'mother' or lc($model) =~ /denovo/) {
				$name = qq{<button style='color:lightgrey;' onClick="view_var_from_proj_gene_pat('$project_name', '$gene_tr_name', '$patient_input_name', '$var_id', 'mother')">$patient_name</button>};
			}
			$type = "ho" if $var->isHomozygote($mother);
			$type = "he" if $var->isHeterozygote($mother);
			my $style = "";
			my $tab;
			if ($var->getProject->isGenome() && $var->isCnv) {
				my $cnv_score = sprintf("%.2f", log2($mother->cnv_value_dude($var->getChromosome->name,$var->start,$var->start+$var->length)));
				my $pr = $var->pr($mother);
				my $sr = $var->sr($mother);
				if ($pr eq '-' and $sr eq '-') {
					my $tab = [$name,$sex_status,$type,$minus,$minus,$text_norm_dp_m,$cnv_score,$minus];
					$html .=  $cgi->start_Tr({style=>"border: 1px solid black;color:lightgrey;$style"}).$cgi->td($tab).$cgi->end_Tr();
					$html_line = join(":",@$tab);
					$html_line .=" - &nbsp;";
				}
				else {
					
					my $tab = [$name,$sex_status,$type,$pr,$sr,$text_norm_dp_m,$cnv_score,$plus];
					
					$html .=  $cgi->start_Tr({style=>"border: 1px solid black;color:black;$style"}).$cgi->td($tab).$cgi->end_Tr();
					$html_line .= join(":",@$tab);
					$html_line .=" - &nbsp;";
				}
			}
			elsif ($var->getPourcentAllele($mother) eq "-"){
				my $ps = "-";
				$tab = [$name,$sex_status,$type,$minus,$text_dp_m,$minus];
				$html .=  $cgi->start_Tr({style=>"border: 1px solid black;color:lightgrey;$style"}).$cgi->td($tab).$cgi->end_Tr();
				$html_line = join(":",@$tab);
				$html_line .=" - &nbsp;";
			}
			else {
				my $tab = [$name,$sex_status,$type,$var->getPourcentAllele($mother)."%",$var->getDepth($mother),"$plus"];
				$html .=  $cgi->start_Tr({style=>"border: 1px solid black;color:black;$style"}).$cgi->td($tab).$cgi->end_Tr();
				$html_line .= join(":",@$tab);
				$html_line .=" - &nbsp;";
			}
		}
		if ($father){
			my $mean_dp_f = int($father->meanDepth($chr_id, $start, $end));
			my $text_dp_f =  $mean_dp_f;
			my $norm_dp_f = $father->cnv_region_ratio_norm($chr_id, $start, $end);
			my $text_norm_dp_f = int($norm_dp_f);
			if ($father->isIll()) { $sex_status = "<center><img src='/icons/Polyicons/male-d.png'></center>"; }
			else { $sex_status = "<center><img src='/icons/Polyicons/male-s.png'></center>"; }
			my $patient_name = $father->name();
			my $name = qq{<button style='color:black;' onClick="view_var_from_proj_gene_pat('$project_name', '$gene_tr_name', '$patient_input_name', '$var_id', 'father')">$patient_name</button>};
			if (lc($model) eq 'father' or lc($model) =~ /denovo/) {
				$name = qq{<center><span style='text-align:center'>$patient_name</span></center>};
			}
			my $style = "";
			#$style = "#F7C9C9" ;
			$type = "-";
			$type = "ho" if $var->isHomozygote($father);
			$type = "he" if $var->isHeterozygote($father);
			#$style = "#779ECB" ;
			if ($var->getProject->isGenome() && $var->isCnv) {
				my $cnv_score = sprintf("%.2f", log2($father->cnv_value_dude($var->getChromosome->name,$var->start,$var->start+$var->length)));
				my $pr = $var->pr($father);
				my $sr = $var->sr($father);
				if ($pr eq '-' and $sr eq '-') {
					my $tab = [$name,$sex_status,$type,$minus,$minus,$text_norm_dp_f,$cnv_score,$minus];
					$html .=  $cgi->start_Tr({style=>"border: 1px solid black;color:lightgrey;$style"}).$cgi->td($tab).$cgi->end_Tr();
					$html_line = join(":",@$tab);
					$html_line .=" - &nbsp;";
				}
				else {
					my $tab = [$name,$sex_status,$type,$pr,$sr,$text_norm_dp_f,$cnv_score,"$plus"];
					$html .=  $cgi->start_Tr({style=>"border: 1px solid black;color:black;$style"}).$cgi->td($tab).$cgi->end_Tr();
					$html_line .= join(":",@$tab);
					$html_line .=" - &nbsp;";
				}
			}
			elsif ($var->getPourcentAllele($father) eq "-"){
				my $tab = [$name,$sex_status,$type,$minus,$text_dp_f,$minus];
				$html .=  $cgi->start_Tr({style=>"border: 1px solid black;color:lightgrey;$style"}).$cgi->td($tab).$cgi->end_Tr();
				push(@$all_tab,@$tab);
				#splice @$tab, 1, 1;
				$html_line .= join(":",@$tab);
				$html_line .="- &nbsp;";
		
			}
			else {
				my $tab = [$name,$sex_status,$type,$var->getPourcentAllele($father)."%",$var->getDepth($father),"$plus"];
				$html .=  $cgi->start_Tr({style=>"border: 1px solid black;color:black;$style"}).$cgi->td($tab).$cgi->end_Tr();
				push(@$all_tab,@$tab);
				#splice @$tab, 1, 1;
				$html_line .= join(":",@$tab);
				$html_line .=" - &nbsp;";
			}
		}
		foreach my $child (@$children) {
			my $mean_dp_c = $child->meanDepth($chr_id, $start, $end);
			my $text_dp_c = sprintf("%.2f", $mean_dp_c);
			my $norm_dp_c = $child->cnv_region_ratio_norm($chr_id, $start, $end);
			my $text_norm_dp_c = int($norm_dp_c);
			my $model = $var->getTransmissionModel($fam,$child,$gene_used);
			my $isMotherTransmission = $var->isMotherTransmission($fam,$child);
			my $isFatherTransmission = $var->isFatherTransmission($fam,$child);
			my $patient_name = $child->name();
			my $name = qq{<button onClick="view_var_from_proj_gene_pat('$project_name', '$gene_tr_name', '$patient_name', '$var_id')">$patient_name</button>};
			
			if ($child->isIll()) {
				if ($child->sex() eq '1') { $sex_status = "<center><img src='/icons/Polyicons/baby-boy-d.png'></center>"; }
				else { $sex_status = "<center><img src='/icons/Polyicons/baby-girl-d.png'></center>"; }
			}
			else {
				if (lc($model) =~ /denovo/) {
					$name = qq{<center><span style='text-align:center'>$patient_name</span></center>};
				}
				if ($child->sex() eq '1') { $sex_status = "<center><img src='/icons/Polyicons/baby-boy-s.png'></center>"; }
				else { $sex_status = "<center><img src='/icons/Polyicons/baby-girl-s.png'></center>"; }
			}
		 $type = "-";
#		 warn $child->name;
		$type = "ho" if $var->isHomozygote($child);
		$type = "he" if $var->isHeterozygote($child);
			if ($type eq "-"){
				my $tab = [$name,$sex_status,$type,$minus,$minus,$minus];
				$html .=  $cgi->start_Tr({class=>"table-info",style=>"border: 1px solid black;color:black;background-color:lightgrey"}).$cgi->td({class=>"table-primary"},$tab).$cgi->end_Tr();
					push(@$all_tab,@$tab);
					#splice @$tab, 1, 1;
				
				$html_line .= join(":",@$tab);
				$html_line .=" - &nbsp;";
				next;
			}
		
		if ($model =~ /denovo/ && scalar (keys%{$fam->parents()}) ==1){
			$model = "denovo/?";
		}
		$hvariation->{transmission_model} = lc ($model);
		$hvariation->{transmission_model_m} = 1 if ($hvariation->{transmission_model} eq "mother" or ($hvariation->{transmission_model} eq "mother_c"));
		$hvariation->{transmission_model_p} = 1 if ($hvariation->{transmission_model} eq "father" or ($hvariation->{transmission_model} eq "father_c"));
		$hvariation->{transmission_model} = "xor" if $hvariation->{transmission_model} eq "mother" or $hvariation->{transmission_model} eq "father";#lc ($model);
		#my $color = "white";
		
		
		my $model2 = qq{<i class="fa fa-male  fa-2x" style="color:lightgrey"></i><i class="fa fa-female  fa-2x" style="color:lightgrey"></i>};
		if (lc($model) eq "father" or (lc($model) eq 'father_c' and $isFatherTransmission)) {
			$color = "#779ECB";
			my $text_compound = ' ';
			#$text_compound = ' father_c' if(lc($model) eq 'father_c');
			if ($child->getFamily->getFather && $child->getFamily->getFather->isIll() ) {
				$model2 = "<center><img src='/icons/Polyicons/male-d.png'>$text_compound</center>";;
			}
			else {
				$model2 = "<center><img src='/icons/Polyicons/male-s.png'>$text_compound</center>";
			}
		}
		if (lc($model) eq "mother" or (lc($model) eq 'mother_c')) {
			$color = "#F7C9C9";
			my $text_compound = ' ';
			#$text_compound = ' mother_c' if(lc($model) eq 'mother_c');
			if ($child->getFamily->getMother && $child->getFamily->getMother->isIll()) {
				$model2 = "<center><img src='/icons/Polyicons/female-d.png'>$text_compound</center>";
			}
			else { $model2 = "<center><img src='/icons/Polyicons/female-s.png'>$text_compound</center>"; }
		}
		#$color = "#779ECB" if (lc($model) eq "father");
		#$model2 = qq{<i class="fa fa-male  fa-2x" style="color:$color"></i>} if (lc($model) eq "father");
		
		#$color = "#F7C9C9" if (lc($model) eq "mother");
		#$model2 = qq{<i class="fa fa-female  fa-2x" style="color:$color"></i>} if (lc($model) eq "mother");
		#$color = "#e74c3c" if (lc($model) =~ "denovo");
		$model2 = qq{Denovo} if (lc($model) =~ "denovo");
		$model2 = qq{Strict Denovo} if (lc($model) =~ "strict");
		#$color = "#E55137" if (lc($model) eq "denovo/?");	
		$model2 = qq{Denovo/?} if (lc($model) eq "denovo/?");
	#	$color = "violet" if (lc($model) eq "recessive");
		$model2 = qq{Recessive} if (lc($model) eq "recessive");
		$color = "#F9885C" if (lc($model) eq "mosaic");
		$model2 = qq{mosaic }.$var->isMosaicTransmission($fam,$patient)  if (lc($model) eq "mosaic");
		$model2 = qq{Uniparental}  if (lc($model) =~ /uniparental/);
		
		my ($pr, $sr, $sd_value, $cnv_score);
		if ($var->getProject->isGenome() && $var->isCnv) {
			$pr = $var->pr($child);
			$sr = $var->sr($child);
			$cnv_score = sprintf("%.2f", log2($child->cnv_value_dude($var->getChromosome->name,$var->start,$var->start+$var->length)));
			#$cnv_score .= '<br>SD:'.sprintf("%.2f",$child->sd_value_patient($var->getChromosome->name,$var->start,$var->start+$var->length));
		}
		
		
		my $perc_allele = $var->getPourcentAllele($child);
		$perc_allele .= "%" if ($perc_allele ne '-');
		
		my $tab;
		if ($var->isCnv) {
			$tab = [$name,$sex_status,$type,$pr,$sr,$text_norm_dp_c,$cnv_score,$model2];
		}
		else {
			$tab = [$name,$sex_status,$type,$perc_allele,$var->getDepth($child),$model2];
		}
	
		
		$html .=  $cgi->start_Tr({class=>"table-info",style=>"border: 1px solid black;color:black;background-color:$color"}).$cgi->td({class=>"table-primary ",style=>"background-color:$color"},$tab).$cgi->end_Tr();
		push(@$all_tab,@$tab);
		#splice @$tab, 1, 1;
		$html_line .= join(":",@$tab);
		$html_line .=" - &nbsp;";

		}
		$html.= $cgi->end_table();
	#	die();
		$hvariation->{html}->{trio} = $html;
		
	#	if ($print ==1){

		#$hvariation->{trio} =$html_line;
		#}
		
	 # }
	
	 
}

sub log2 {
	 my $n = shift;
	 return -2 if $n <0.01;
	 my $v = log($n)/log(2);
	 $v =-2 if $v < -2; 
    return $v;
}


sub return_coverage {
	my ($p,$chr,$start) = @_;
	return -1 unless $p->isNoSqlDepth;
	return $p->depth($chr,$start,$start+1)->[0];
}


##########
#panel gene 
#@@@@@@@@


sub panel_gene_short_button {
	my ($hgene, $project_name,$list_html_cons) = @_;
	my $cgi          = new CGI();
	my $out;
	my $gene_id = $hgene->{id};
	my $homim;
	my $vval ={};
	my $label_id = "label_".$hgene->{js_id};
	my $max_score = $hgene->{max_score};
	my $bgcolor2 = "background-color:#F3F3F3;border-color:#607D8B";
	my $glyph = "";
	$glyph = qq{<span class="glyphicon glyphicon-star-empty text-default" aria-hidden="true"></span>} if $hgene->{nb_clinvar} > 0;
	$glyph = qq{<span class="glyphicon  glyphicon-alert text-alert" aria-hidden="true" style="color:red"></span>} if $hgene->{nb_clinvar_alert} > 0 or $hgene->{dm};
	$glyph = qq{<img src="https://img.icons8.com/color/24/000000/treatment-plan.png" style="float:right">} if exists $vval->{validation};# if $hgene->{validations} >2;
	my $astyle = "";
	my $gene_name = $hgene->{external_name};
	$out .=  $cgi->start_div({class=>" btn-group btn-xs ",style=>'position:relative;float:left;'});
	my $bcolor = "grey";
	$bcolor = "yellow" if $max_score >= 5;
	$bcolor = "orange" if $max_score >= 8;
	$bcolor = "coral" if $max_score >= 12;
	$bcolor = "red" if $max_score >= 14;
	my $div_id_gene = "div_".$gene_name;
	$div_id_gene =~ s/\./_/g;
	$div_id_gene =~ s/ //g;
	
	my $panel_name1 = join("-",keys %{$hgene->{panels}});
   	my $hPanels;
   	foreach my $panel_name (keys %{$hgene->{panels}}) {
   		my $pheno_name = $hgene->{panels}->{$panel_name}->{phenotype};
   		$hPanels->{$pheno_name}->{$panel_name} = undef;
  	}
   	my $panel_list;
   	foreach my $pheno_name (sort keys %$hPanels) {
   		$panel_list .= "<b><u>Phenotype: ".$pheno_name."</b></u><br>";
   		foreach my $panel_name (sort keys %{$hPanels->{$pheno_name}}) {
   			$panel_list .= '  - '.$panel_name."<br>";
   		}
   		$panel_list .= "<br>";
   	}
   	$panel_name1 = "Panel: ".scalar (keys %{$hgene->{panels}});
	my $b_panels = qq{<div class="btn btn-brown btn-xs" style="color:white;background-color:#4A4F53;$astyle;font-family: Verdana,Arial,sans-serif; text-shadow:1px 1px 2px black;" onclick="document.getElementById('span_list_panels').innerHTML='$panel_list';dijit.byId('dialog_list_panels').show();"><span style="font-size:10px;">$panel_name1</span></div>} if $panel_name1;
	
	my $pheno =$hgene->{phenotypes};
	my @t = split(";",$pheno);
	$pheno = $t[0];
	my $nb_other_terms = scalar(@t) - 1;
	my $b_desc = qq{<span style='color:white;'>}.$hgene->{description}.qq{</span>};
	my $to;
	($pheno,$to) = split(/\[/,$hgene->{description}) unless $pheno;
	my $color ;
	$color = qq{ style = "color:#E74C3C"} if $pheno =~/intellectual/ or $pheno =~/mental/ or $pheno =~/retar/;
	my $jtid = 'zz'.time.rand(500);
	my $b_pheno = qq{<div class="btn btn-brown btn-xs $bcolor" style="font-size:10px;color:white;background-color:#4A4F53;$astyle;font-family: Verdana,Arial,sans-serif; text-shadow:1px 1px 2px black;position:relative;">};
	if ($pheno) {
		if (length($pheno) > 70) { $pheno = substr($pheno, 0, 70).'...'; }
		if ($nb_other_terms > 0) { $pheno .= " <span style='color:#5cf0d3'>+ $nb_other_terms terms</span>"; }
#		if ($project_name) {
#			$b_pheno .= qq{<span onclick="update_grid_gene_phenotypes(\'$gene_id\', \'$project_name\')">$pheno</span>};
#		}
#		else {
			$b_pheno .= qq{<span onclick="update_grid_gene_phenotypes(\'$gene_id\')">$pheno</span>};
#		}
	}
	$b_pheno .= qq{</div>};
	
	my $in = $hgene->{omim_inheritance};
	$in ="" if $in eq "-";
	$in = "X-linked " if $in =~/X-linked/;
	
	$out .= "<table>";
	$out .= "<tr>";
	$out .= "<td>";
	if (exists $hgene->{collapse_with_id}) {
		my $this_collapse_id = $hgene->{collapse_with_id};
		$out .= qq{<div id="$div_id_gene" class="btn btn-brown btn-xs $bcolor" data-toggle='collapse' data-target="#$this_collapse_id" aria-expanded='false' aria-controls='$this_collapse_id' style="font-size:12px;color:white;background-color:#4A4F53;border-right: 4px solid $bcolor;border-left: 4px solid $bcolor;$astyle;font-family: Verdana,Arial,sans-serif; text-shadow:1px 1px 2px black;position:relative;bottom:0px;min-width:150px;"> </span> <span id= "$label_id" class="glyphicon glyphicon-triangle-right" aria-hidden="true" style="padding-top:2px;float:left;"></span> $gene_name<sup>&nbsp;$in</b></sup>&nbsp;&nbsp;&nbsp;</div></div>};
	}
	else {
		$out .= qq{<div id="$div_id_gene" class="btn btn-brown btn-xs $bcolor" style="font-size:12px;color:white;background-color:#4A4F53;border-right: 4px solid $bcolor;border-left: 4px solid $bcolor;$astyle;font-family: Verdana,Arial,sans-serif; text-shadow:1px 1px 2px black;position:relative;bottom:0px;min-width:150px;"> </span> <span id= "$label_id" class="glyphicon glyphicon-triangle-right" aria-hidden="true" style="padding-top:2px;float:left;"></span> $gene_name<sup>&nbsp;$in</b></sup>&nbsp;&nbsp;&nbsp;</div></div>};
	}

	$out .= "</td>";
	if ($list_html_cons) {
		$out .= '<td>&nbsp;&nbsp;'.join('&nbsp;', @$list_html_cons).'</td>';
	}
	$out .= "<td>".$b_panels."</td>";
	$out .= "<td>".$b_pheno."</td>";
	$out .= "<td>".$b_desc."</td>";
	$out .= "</table>";
	return $out;
} 
sub panel_gene { 
	my ($hgene,$panel_id,$project_name,$patient) = @_;
	my $cgi          = new CGI();
	my $out;
	my $gene_id = $hgene->{id};
	my $gene;
	$gene = $patient->project->newGene($gene_id) if ($patient);
	
	#my $homim;
	#$homim = $patient->project->lmdbOmim->get($gene_id) if ($patient);
	
	#warn $omim->{phenotype}->{omim};
	my $vval ={};
	my $all_validations;
	($all_validations) = grep {$_ =~ /$gene_id/ } keys %{$patient->validations} if ($patient);
	if ($all_validations) {
		 $vval =  $patient->validations->{$all_validations}->[0];
		#warn Dumper  $patient->validations->{$all_validations} if $all_validations;
	}
	
	
	#$out .=  $cgi->start_div({class=>"panel panel-success" });
	 #panel heading
	
	#$out .=  $cgi->start_div({class=>"panel-heading panel-warning warning ",style=>" min-height:13px;max-height:13px;padding:10px;border:1px"});
	#$out .=  $cgi->start_div({class=>" btn-group btn-xs "});
	my $label_id = "label_".$hgene->{uid};
	my $max_score = $hgene->{max_score};
	
my $bgcolor2 = "background-color:#607D8B;border-color:#607D8B";
	 #panel heading
	#$out .=  $cgi->start_div({class=>"panel-heading panel-warning warning ",style=>" min-height:13px;max-height:13px;padding:1px;border:1px"});
	 #$out.="coucou";
			my $glyph = "";
				$glyph = qq{<span class="glyphicon glyphicon-star-empty text-default" aria-hidden="true"></span>} if $hgene->{nb_clinvar} > 0;
				$glyph = qq{<span class="glyphicon  glyphicon-alert text-alert" aria-hidden="true" style="color:red"></span>} if $hgene->{nb_clinvar_alert} > 0 or $hgene->{dm};
				$glyph = qq{<img src="https://img.icons8.com/color/24/000000/treatment-plan.png" style="float:right">} if exists $vval->{validation};# if $hgene->{validations} >2;
				my $astyle = "";
				
				#$astyle = "background-color:#FF8800" ;# if $vval->{validation}>5;#if $hgene->{validations} >2 ;
				$astyle = "background-color:#E74C3C"  if $vval->{validation}>4;
				my $gene_name = $hgene->{external_name};
				
				$out .=  $cgi->start_div({class=>" btn-group btn-xs ",style=>'position:relative;float:left;bottom:5px;'});
				my $in;
				if ($gene) {
					if ($gene->getProject->getVersion() =~ /HG/) { $in = $gene->omim_inheritance; }
					else {$in = '-';}

				}
				else { $in = $hgene->{omim_inheritance}; }
				$in ="" if $in eq "-";
				$in = "X-linked " if $in =~/X-linked/;
				my $pli = 	$hgene->{pLI}*1.0;
				#$pli = $pli.;#"/".$hgene->{max_score};
				
				
				my $bcolor = "grey";
				$bcolor = "green" if $max_score >= 0;
				$bcolor = "yellow" if $max_score >= 5;
				$bcolor = "orange" if $max_score >= 8;
				$bcolor = "coral" if $max_score >= 12;
				$bcolor = "red" if $max_score >= 14;
			
				my $div_id_gene = "div_".$gene_name;
				$div_id_gene =~ s/\./_/g;
				$div_id_gene =~ s/ //g;
				
				my $cnv_status = "cnv_none";
				if ($hgene->{level_dude} ne '-1') {
					my ($l,$t) = split(":",$hgene->{level_dude});
					if ($t eq "del"){
						$cnv_status = "cnv_del_medium" if $l eq "medium";
						$cnv_status = "cnv_del_high" if $l eq "high";
						$cnv_status = "cnv_del_low" if $l eq "low";
					}
					if ($t eq "dup"){ 
						$cnv_status = "cnv_dup_medium" if $l eq "medium";
						$cnv_status = "cnv_dup_high" if $l eq "high";
						$cnv_status = "cnv_dup_low" if $l eq "low";
					}
				}
				#
				my $this_b_cmd = qq{collapse("$panel_id","$label_id")};
				if (exists $hgene->{specific_cmd}) {
					$this_b_cmd = $hgene->{specific_cmd};
				}
				
				
				if (exists $hgene->{collapse_with_id}) {
					my $this_collapse_id = $hgene->{collapse_with_id};
					$out .= qq{<div id="$div_id_gene" class="btn btn-brown btn-xs $bcolor $cnv_status" data-toggle='collapse' data-target="#$this_collapse_id" aria-expanded='false' aria-controls='$this_collapse_id' style="background-color:#4A4F53;border-top: 2px solid #4A4F53;border-bottom: 2px solid #4A4F53;border-right: 4px solid $bcolor;border-left: 4px solid $bcolor;$astyle;font-family: Verdana,Arial,sans-serif; text-shadow:1px 1px 2px black;position:relative;bottom:0px;min-width:150px;" onClick='$this_b_cmd'>  <font style='color:white;'><span id= "$label_id" class="glyphicon glyphicon-triangle-right" style="" aria-hidden="true"  style="float:left;top:4px;"></span> $gene_name<sup>&nbsp;$in</b></font></sup> $glyph}.$cgi->span({class=>"badge1 $bcolor"},$hgene->{max_score}).qq{</div>};
				}
				else {
					$out .= qq{<div id="$div_id_gene" class="btn btn-brown btn-xs $bcolor $cnv_status" style="background-color:#4A4F53;border-top: 2px solid #4A4F53;border-bottom: 2px solid #4A4F53;border-right: 4px solid $bcolor;border-left: 4px solid $bcolor;$astyle;font-family: Verdana,Arial,sans-serif; text-shadow:1px 1px 2px black;position:relative;bottom:0px;min-width:150px;" onClick='$this_b_cmd'>  <font style='color:white;'><span id= "$label_id" class="glyphicon glyphicon-triangle-right" style="" aria-hidden="true"  style="float:left;top:4px;"></span> $gene_name<sup>&nbsp;$in</b></font></sup> $glyph}.$cgi->span({class=>"badge1 $bcolor"},$hgene->{max_score}).qq{</div>};
				}
				
				
				
	   			my $nbv = $hgene->{nb};

				my $omim = $hgene->{omim_id};
				$out .=qq{<a class="btn btn-primary btn-xs" href="http://www.omim.org/entry/$omim" role="button" target="_blank" style="min-width:40px;$bgcolor2;text-shadow:1px 1px 1px black;color:white">Omim</a>} if $omim ne "";
				$out .=qq{<a class="btn btn-primary btn-xs"   style="$bgcolor2;min-width:40px"><i class="fa fa-minus"></i></a>} if $omim eq "";
				
				#gtex Portal
					my ($gid,$t) = split("_",$hgene->{id});
					
					$out .=qq{<a class="btn btn-primary btn-xs" href="https://gtexportal.org/home/gene/$gid" role="button" target="_blank" style="min-width:40px;$bgcolor2;text-shadow:1px 1px 1px black;color:white">Gtex</a>};# if $omim ne "";
				#ENSG00000124155
			my $oid = $hgene->{name};
				my $type ="green";
				#$type = "default" if $pli <= 0.1;
				$type = "orange" if $pli >= 0.75;
				$type = "red" if $pli >= 0.9;
				my $m = $hgene->{max_score};
				#$out .=qq{<a class="btn btn-primary btn-xs" href="https://gnomad.broadinstitute.org/gene/$oid" target="_blank" style="$bgcolor2;min-width:30px;height:22px;padding-top:3px;"><span class="badge" style="color:$type">$pli</span></a>};
 				
 				
				my $dataset = "?dataset=gnomad_r4";
				$dataset = "?dataset=gnomad_r2_1" if $gene and $gene->getProject->getVersion() =~ /HG19/;
				$dataset = "?dataset=gnomad_r2_1" if $patient and $patient->getProject->getVersion() =~ /HG19/;
 				my $b_id_pli = 'b_pli_'.$oid.'_'.$type;
 				my $popup_pli = qq{<div data-dojo-type="dijit/Tooltip" data-dojo-props="connectId:'$b_id_pli',position:['above']"><span><b>pLI</b> Score</span></div>};
 				if ($gene) {
 					my ($gidtmp,$gtmp) = split('_',$gene->id());
 					$out .=qq{<a class="btn btn-primary btn-xs" href="https://gnomad.broadinstitute.org/gene/$gidtmp$dataset" target="_blank" style="$bgcolor2;min-width:30px"><span id="$b_id_pli" class="badge" style="color:$type">$pli</span>$popup_pli</a>};
 				}
 				else {
 					my ($gidtmp,$gtmp) = split('_',$hgene->{id});
 					$out .=qq{<a class="btn btn-primary btn-xs" href="https://gnomad.broadinstitute.org/gene/$gidtmp$dataset" target="_blank" style="$bgcolor2;min-width:30px"><span id="$b_id_pli" class="badge" style="color:$type">$pli</span>$popup_pli</a>};
 				}
 				
 				
 				
	   			my $panel_name1 = join("-",keys %{$hgene->{panels}});
	   			my $hPanels;
	   			foreach my $panel_name (keys %{$hgene->{panels}}) {
	   				my $pheno_name = $hgene->{panels}->{$panel_name}->{phenotype};
	   				$hPanels->{$pheno_name}->{$panel_name} = undef;
	   			}
	   			my $panel_list;
	   			foreach my $pheno_name (sort keys %$hPanels) {
	   				$panel_list .= "<b><u>Phenotype: ".$pheno_name."</b></u><br>";
	   				foreach my $panel_name (sort keys %{$hPanels->{$pheno_name}}) {
	   					$panel_list .= '  - '.$panel_name."<br>";
	   				}
	   				$panel_list .= "<br>";
	   			}
	   			$panel_name1 = "Panel: ".scalar (keys %{$hgene->{panels}});
				$out .=qq{<a class="btn btn-primary btn-xs" href="#" role="button" style="top:-4px;$bgcolor2" onclick="document.getElementById('span_list_panels').innerHTML='$panel_list';dijit.byId('dialog_list_panels').show();"><p style="font-size:10px;text-shadow:0px 1px 1px #000;position:relative;bottom:-4px">$panel_name1</p></a>} if $panel_name1;

	   		my ($pheno,$nb_other_terms);
	   		if ($gene) {
	   			if ($gene->getProject->getVersion() =~ /HG/) { 
	   				($pheno,$nb_other_terms) = $gene->polyviewer_phentotypes();
	   			}
	   			else {
	   				$pheno = $gene->description();
	   			}
	   			
	   		}
	   		
	   		elsif (exists $hgene->{phenotypes} and exists $hgene->{phenotypes}) {
				$pheno = $hgene->{phenotypes};
				$nb_other_terms = scalar(split(";",$pheno));
	   		}
	   		
   			my $color ;
   			$color = qq{ style = "color:#E74C3C"} if $pheno =~/intellectual/ or $pheno =~/mental/ or $pheno =~/retar/;
   			my $jtid = 'zz'.time.rand(500);
   			my $div_pheno = qq{<a aria-disabled="true" class="btn btn-primary btn-xs" href="#" role="button" style="text-align:left;font-family: proxima-nova, sans-serif;font-style:normal;text-shadow:0px 1px 1px #000000;position:relative;bottom:4px;font-size: 0.9em;color:white;$bgcolor2;">};
   			if ($pheno) {
		   		if (length($pheno) > 70) { $pheno = substr($pheno, 0, 70).'...'; }
   				if ($nb_other_terms > 0) { $pheno .= " <span style='color:#5cf0d3'>+ $nb_other_terms terms</span>"; }
   				if ($project_name) {
   					$div_pheno .= qq{<i class="fa fa-circle fa-xs" $color ></i> <span onclick="update_grid_gene_phenotypes(\'$gene_id\', \'$project_name\')";">$pheno</span>};
   				}
   				else {
   					$div_pheno .= qq{<i class="fa fa-circle fa-xs" $color ></i> <span onclick="update_grid_gene_phenotypes(\'$gene_id\')";">$pheno</span>};
   				}
   			}
   			$div_pheno .= qq{</a>};
	   		$out .= $div_pheno;
	   		$out .= $cgi->end_div();
	   		
				$out .=  $cgi->start_div({class=>" btn-group btn  ",style=>'position:relative;float:right;bottom:5px;'});
			
						my $tlabel = "label-grey";
						$tlabel = "label-warning" if exists $hgene->{clinvar_hgmd};
						$tlabel = "label-danger" if exists $hgene->{pathogenic};
						$out .= $cgi->span({class=>"label $tlabel" },qq{<span class="glyphicon glyphicon-alert text-alert " aria-hidden="true" ></span>});
					
						my $nbv = scalar(@{$hgene->{variants}});
						$out .=$cgi->span({class=>"label label-grey"},qq{<span class='badge badge-infos badge-xs ' style="color:#00C851"  >$nbv </span> });
						
					
						$out .=$cgi->span({class=>"label label-grey"},qq{<span class="glyphicon glyphicon glyphicon-menu-hamburger " aria-hidden="true" style="color:orange" ></span>});
						my $clabel = " label-default";
						$clabel = "label-danger" if $all_validations;#if exists $hgene->{validations};
						$out .=$cgi->span({class=>"label $clabel"},qq{<span class="glyphicon glyphicon-pencil " aria-hidden="true" ></span>});
						my $style ="";
						$clabel = " label-grey";
						if ($hgene->{level_dude} eq '-1') {
							$out .=$cgi->span({class=>"label $clabel",style=>$style,onclick=>"alert('Not Available - Please contact BioInformatic Platform (genes_level_dude)');"},qq{<s><span>-</span></s> <span class="glyphicon glyphicon-ban-circle"></span>});
						}
						else {
							my ($l,$t) = split(":",$hgene->{level_dude});
							if ($t eq "del"){
								$style = "background-color:#FF6961" if $l eq "medium";
								$style = "background-color:#DD4124" if $l eq "high";
								$style = "background-color:#FFB2AE" if $l eq "low";
							}
							if ($t eq "dup"){ 
								$style = "background-color:#0080FF" if $l eq "medium";
								$style = "background-color:#0066FF" if $l eq "high";
								$style = "background-color:#9DFFFF" if $l eq "low";
							}
							#$clabel = "label-warning" if  $hgene->{level_dude} eq "medium";
							#$clabel = "label-danger" if  $hgene->{level_dude} eq "high";
							#onclick="zoomDude(\'$gnames\','50')"
							my $gnames = $hgene->{name};
							my $pname;
							if ($patient) {
								$pname = $patient->name();
								$out .=$cgi->span({class=>"label $clabel",style=>$style,onclick=>"zoomDude(\'$project_name\', \'$gnames\','$pname', 'force_visualisation')"},qq{<span>CNV</span>});
							}
							else {
								$out .=$cgi->span({class=>"label $clabel",style=>$style,onclick=>"zoomDude(\'$project_name\', \'$gnames\','', 'force_visualisation')"},qq{<span>CNV</span>});
							}
						}
				
				
			 	$out.= $cgi->end_div(); # end div lavel right 
}

sub table_line_gene {
	my ($hgene, $project_name) = @_;
	my $cgi = new CGI();
	my $out;
	my $gene_id = $hgene->{id};
	my $max_score = $hgene->{max_score};
	
	my $glyph = "";
	$glyph = qq{<span class="glyphicon glyphicon-star-empty text-default" aria-hidden="true"></span>} if $hgene->{nb_clinvar} > 0;
	$glyph = qq{<span class="glyphicon  glyphicon-alert text-alert" aria-hidden="true" style="color:red"></span>} if $hgene->{nb_clinvar_alert} > 0 or $hgene->{dm};
	
	my $gene_name = $hgene->{external_name};
	my $in = $hgene->{omim_inheritance};
	$in ="" if $in eq "-";
	$in = "X-linked " if $in =~/X-linked/;
	
	my $pli = 	$hgene->{pLI}*1.0;
				
	my $bcolor = "grey";
	$bcolor = "green" if $max_score >= 0;
	$bcolor = "yellow" if $max_score >= 5;
	$bcolor = "orange" if $max_score >= 8;
	$bcolor = "coral" if $max_score >= 12;
	$bcolor = "red" if $max_score >= 14;
			
	my $cnv_status = "cnv_none";
	if ($hgene->{level_dude} ne '-1') {
		my ($l,$t) = split(":",$hgene->{level_dude});
		if ($t eq "del"){
			$cnv_status = "cnv_del_medium" if $l eq "medium";
			$cnv_status = "cnv_del_high" if $l eq "high";
			$cnv_status = "cnv_del_low" if $l eq "low";
		}
		if ($t eq "dup"){ 
			$cnv_status = "cnv_dup_medium" if $l eq "medium";
			$cnv_status = "cnv_dup_high" if $l eq "high";
			$cnv_status = "cnv_dup_low" if $l eq "low";
		}
	}
	my $this_b_cmd;
	if (exists $hgene->{specific_cmd}) {
		$this_b_cmd = $hgene->{specific_cmd};
	}
	if ($this_b_cmd) {
		$out .= qq{<td><button class="btn btn-primary btn-xs" style="border-color:transparent;background-color:transparent;color:black;font-size:12px;" onClick='$this_b_cmd'><u>$gene_name</u><sup>&nbsp;$in</b></sup> $glyph}.$cgi->span({class=>"badge1 $bcolor"},$hgene->{max_score}).qq{</button></td>};
	}
	else {
		$out .= qq{<td><a class="btn btn-primary btn-xs" style="border-color:transparent;background-color:transparent;color:black;font-size:12px;"  href="https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=$gene_name" target="_blank"><u>$gene_name</u><sup>&nbsp;$in</b></sup> $glyph}.$cgi->span({class=>"badge1 $bcolor"},$hgene->{max_score}).qq{</a></td>};
	}

	my $nbv = $hgene->{nb};
	my $omim = $hgene->{omim_id};
	
	$out .= qq{<td><a class="btn btn-primary btn-xs" href="http://www.omim.org/entry/$omim" role="button" target="_blank" style="min-width:40px;border-color:transparent;background-color:transparent;color:black;"><u>Omim</u></a></td>} if $omim ne "";
	$out .= qq{<td><i class="fa fa-minus"></i></td>} if $omim eq "";
				
	my ($gid,$t) = split("_",$hgene->{id});
	$out .=qq{<td><a class="btn btn-primary btn-xs" href="https://gtexportal.org/home/gene/$gid" role="button" target="_blank" style="min-width:40px;border-color:transparent;background-color:transparent;color:black;"><u>Gtex</u></a></td>};
	
	my $oid = $hgene->{external_name};
	my $type ="green";
	$type = "orange" if $pli >= 0.75;
	$type = "red" if $pli >= 0.9;
	my $m = $hgene->{max_score};
	$out .=qq{<td><center><a class="btn btn-xs" href="https://gnomad.broadinstitute.org/gene/$oid" target="_blank" style="min-width:30px"><span class="badge" style="background-color:$type;color:white;font-size:11px;">$pli</span></a></center></td>};
	
	my $panel_name1 = join("-",keys %{$hgene->{panels}});
	my $hPanels;
	foreach my $panel_name (keys %{$hgene->{panels}}) {
		my $pheno_name = $hgene->{panels}->{$panel_name}->{phenotype};
		$hPanels->{$pheno_name}->{$panel_name} = undef;
	}
	my $panel_list;
	foreach my $pheno_name (sort keys %$hPanels) {
		$panel_list .= "<br><center><b><u>Phenotype: ".$pheno_name."</b></u></center><br>";
		foreach my $panel_name (sort keys %{$hPanels->{$pheno_name}}) {
			$panel_list .= '<center>'.$panel_name."</center>";
		}
	}
	$panel_list .= "<br>";
	$panel_name1 = scalar (keys %{$hgene->{panels}});
	if ($panel_name1) {
		$out .=qq{<td><center><a class="btn btn-primary btn-xs" href="#" role="button" style="border-color:transparent;background-color:transparent;color:black;font-size:11px;" onclick="document.getElementById('span_list_panels').innerHTML='$panel_list';dijit.byId('dialog_list_panels').show();"><span class="badge" style="background-color:green;color:white;font-size:11px;">$panel_name1</span></a></center></td>};
	}
	else { $out .= qq{<td><center><i class="fa fa-minus"></i></center></td>}; }

	my $pheno = $hgene->{phenotypes}->{pheno};
	my $nb_other_terms = $hgene->{phenotypes}->{nb_other_terms};
	my $color ;
	$color = qq{ style = "color:#E74C3C"} if $pheno =~/intellectual/ or $pheno =~/mental/ or $pheno =~/retar/;
	my $jtid = 'zz'.time.rand(500);
   	my $div_pheno = qq{<a aria-disabled="true" class="btn btn-primary btn-xs" href="#" role="button" style="text-align:left;font-family: proxima-nova, sans-serif;font-style:normal;border-color:transparent;background-color:transparent;color:black;">};
   	if ($pheno) {
		if (length($pheno) > 70) { $pheno = substr($pheno, 0, 70).'...'; }
   		if ($nb_other_terms > 0) {
	   		if ($project_name) {
	   			$pheno .= qq{</u><span style='color:blue' onclick="update_grid_gene_phenotypes(\'$gene_id\', \'$project_name\')";">+ $nb_other_terms terms</span>};
	   		}
	   		else {
	   			$pheno .= qq{</u><span style='color:blue' onclick="update_grid_gene_phenotypes(\'$gene_id\')";">+ $nb_other_terms terms</span>};
	   		}
   		}
   		if ($project_name) {
   			$div_pheno .= qq{<span onclick="update_grid_gene_phenotypes(\'$gene_id\', \'$project_name\')";">$pheno</span>};
   		}
   		else {
   			$div_pheno .= qq{<span onclick="update_grid_gene_phenotypes(\'$gene_id\')";">$pheno</span>};
   		}
   	}
   	$div_pheno .= qq{</a>};
	$out .= '<td><u>'.$div_pheno.'</u></td>';

#	$out .= '</tr></table>';
	return $out;	   		
}
  

sub print_table_nav_solo {
	my ($patient, $stamp,$my_style) =@_;
	my $fam = $patient->getFamily();
	my $name = $patient->name;
	if ($patient->getProject->validation_db() and $patient->getProject->validation_db() eq 'defidiag') {
		$name = $name." - ".$fam->name;
	}
	 my $icon_sex = qq{<mark style="border: 2px solid black;background-color:#D9EEF2"><img src="https://img.icons8.com/color/24/000000/male.png">}.$name."</mark>";
	 $icon_sex =  qq{<mark style="background-color:#FFD1E3"><img src="https://img.icons8.com/office/24/000000/female.png">}.$name."</mark>" unless $patient->isMale();
	 my $icon_warning = q{<mark style="background-color:#CCFF00"><img src="https://img.icons8.com/color/50/000000/break.png">}."Only ACMG Actionable Genes</mark>";#<img src="https://img.icons8.com/color/50/000000/break.png">
	#FFD1E3 pink
	my $text;
	my $n ="";
	if ($patient->project->getPhenotypes()){
		
		#my $pheno = $patient->project->getPhenotypes()->[0];
	 	#$n =  $pheno->name if $pheno;
	}
	my $text = $stamp;
	if ($n) { $text .= qq{Project: <span style="color:rgb(214, 44, 26);">$n</span>}; }
	if ($patient->isMother or $patient->isFather) {
		if ($patient->getProject->validation_db() and $patient->getProject->validation_db() eq 'defidiag') {
			$text .= qq{
			<table  style=" font: Verdana;position:relative;left:10%;width:80%; color: rgba(0,0,0,0.7);font-size:20px">
				<tr style="span:3px"><td style="padding-right:15px" align="center">$icon_warning</td></tr>
				<tr style="span:3px"><td style="padding-right:15px" align="center">$icon_sex</td></tr>
				
				</table>
			};
			return $text;
		}
		else {
			$text .= qq{
			<table  style=" font: Verdana;position:relative;left:10%;width:80%; color: rgba(0,0,0,0.7);font-size:20px">
				<tr style="span:3px"><td style="padding-right:15px" align="center">$icon_sex</td></tr>
				</table>
			};
			return $text;
			
		}
	}
	$my_style = "left:10%;width:80%;" unless ($my_style);
	 $text .= qq{
	<table  style=" font:  Verdana;position:relative;$my_style;color: rgba(0,0,0,0.7);font-size:20px">
		<tr style="span:3px"><td style="padding-right:15px" align="center">$icon_sex</td></tr>
		</table>
	};
	$text.= "<br><br>";
	return $text;
}
sub print_table_nav_trio {
	my ($patient,$stamp,$my_style,$string_filtering) =@_;
#	warn $my_style;
	my $fam = $patient->getFamily();
	my $fam_name = $fam->name();
	my $name = $patient->name;
	my $space= "&nbsp;";
	my $border_s = "box-shadow: 2px 2px 5px 0 rgba(0,0,0,0.75);";
	my $border_d = "border: 1px solid #E51400;box-shadow: 2px 2px 5px 0 rgba(255,0,0,0.75);";
	my $border = $border_s;
	 $border = $border_d if  $patient->status ne "1";
	 my $icon_sex = qq{<mark style="$border;background-color:#D9EEF2">}.$patient->return_icon.$space.$name."</mark>";
	# $icon_sex =  qq{<mark style="border: 2px solid black;background-color:#FFD1E3"><img src="https://img.icons8.com/color/32/000000/girl.png">}.$name."</mark>" unless $patient->isMale();
	#FFD1E3 pink
	
	 	my $bmother ="";
		 $border = $border_s;
		 
		 $border = $border_d if $fam->getMother && $fam->getMother->status ne "1";
		
	 $bmother = qq{<mark style="$border;background-color:#FFD1E3">}.$fam->getMother->return_icon.$space.$fam->getMother->name.qq{</mark>} if $fam->getMother; ;
		my	$bfather = "";
		$border = $border_s;
		$border = $border_d  if $fam->getFather && $fam->getFather->status ne "1";
		$bfather = qq{<mark style="$border;background-color:#D9EEF2">}.$fam->getFather->return_icon.$space.$fam->getFather->name().qq{</mark>}  if $fam->getFather; 
	$my_style = "left:10%;width:80%;" unless ($my_style);
	my $n ="";
	if ($patient->project->getPhenotypes()){
		my $pheno = $patient->project->getPhenotypes()->[0];
	 	$n =  $pheno->name if $pheno;
	}
	my $text;

	$text .= $stamp;
	
	if ($n) { $text .= qq{ <button style="letter-spacing: 2px;position: relative;top:-10px;font-size:13px">Family: <b>$fam_name</b><button>  }; }
	else { $text .= qq{Family: <b>$fam_name</b>}; }
	
	if ($fam->getMother && $fam->getMother->isHealthy && $fam->getFather && $fam->getFather->isHealthy ){
			$text .= qq{<button type ="button" class="btn btn-success" style="letter-spacing: 2px;position: relative;right:-200px;top:-10px;;font-size:13px;color:black;padding:3px;"> $n : <span style="color:white;"><b>Scoring Model  : (Recessive/Denovo)</b></span></button>};
		}
		elsif (scalar(@{$fam->getParents}) == 1  && $fam->getParents->[0]->isHealthy) {
					$text .= qq{<button type ="button" class="btn btn-success" style="letter-spacing: 2px;position: relative;right:-200px;top:-10px;;font-size:13px;color:black;padding:3px;"> $n : <span style="color:white;"><b>Scoring Model  : Recessive + Denovo</b></span></button>};
		}
	else {
				$text .= qq{<button type ="button" class="btn btn-danger" style="letter-spacing: 2px;position: relative;right:-200px;top:-10px;font-size:13px;color:black;padding:3px;"> $n : <span style="color:white;"><b>Score Model  : Dominant inherited</b></span></button>};
				#$text .= qq{&nbsp;-&nbsp;Model: Model:<span style="color:#FFCA3A;">Dominant</span>};
	}
	$text .= qq{
		
	<table  style=" font:  Verdana;position:relative;$my_style;color: rgba(0,0,0,0.7);font-size:20px">
		<tr style="border-bottom: 2px solid black;span:3px"><td style="padding-right:15px" align="center">$bmother</td><td style="padding:15px;" align="center">$bfather</td></tr>
			<tr><td style="border-right: 2px solid black;">&nbsp;</td> <td></td></tr>
			<tr ><td colspan="2" align="center" style="" align="center">$icon_sex</td></tr>
		</table>
		
	};
	return $text;
}

sub compose_string_filtering{
	my ($cgi) = @_;
	my @t;
	my $filter_transmission;
	push(@$filter_transmission,"denovo") if $cgi->param('denovo');	
	#$filter_transmission->{"strict_denovo"} = 1 if $cgi->param('denovo');
	push(@$filter_transmission,"denovo")   if $cgi->param('strict_denovo');
	push(@$filter_transmission,"recessive")  if $cgi->param('recessive');
	push(@$filter_transmission,"composite")  if $cgi->param('xor');
	push(@$filter_transmission,"Mother exclusive")  if $cgi->param('xor_mother');
	push(@$filter_transmission,"Father exclusive")  if $cgi->param('xor_father');
	$filter_transmission = ["All Variants"] if $cgi->param('both');
	return join("+",@$filter_transmission);
	
}

sub printNavBar {
	my ($patient,$hgenes,$statistics,$version,$date,$user,$ztime,$string_filtering) = @_;
	my $nb_var = $statistics->{variations};
	my $nb_dm = $statistics->{DM};
	my $nb_genes = $statistics->{genes};
	my $cgi = new CGI;
	my $name = $patient->name;
	my $all_panel=[];
	my $all_label=[];
	foreach my $gene ( @$hgenes){
		$gene->{uid} = $gene->{id}."_".int(rand(time)) unless exists $gene->{uid};
		push(@$all_panel,"panel_".$gene->{uid});
		push(@$all_label,"label_".$gene->{uid});
		$nb_genes++;
	}

#  if ($fam->isTrio && $patient->isChild()){
#		$small_panel = small_panel_trio($patient,$cgi);
#	}
#	else {
#		$small_panel = small_panel($patient,$cgi);
#	}
my $table ="";
my $validation = $patient->getLatestValidationStatus($user);

	my $term = "";
	my $display = "visibility:hidden;";
	my $date2= "";
	
	if ($validation){
		$term = $validation->{term};
		 $date2 = join("-",return_date($validation->{modification_date}));
		 $display = "";
	}
	my $patient_name = $patient->name;
	my $stamp  = qq{ <span id="stamp_$patient_name" class="stamp is-approved_red" style="$display"><span>$term</span><br><small>$date2</small></span>};
	


 if ($patient->getFamily->isTrio && $patient->isChild()){
 	
 	$table = print_table_nav_trio($patient,$stamp,undef,$string_filtering);
 }
 else {
 	$table =print_table_nav_solo($patient,$stamp);
 }
 my $hidden = "visible";
 
 $hidden = "hidden" if $validation;
		#$table =print_table_nav_solo($patient,$cgi);
	my $string_panel = join(";",@$all_panel);
my $string_label = join(";",@$all_label);
	my $pname = $patient->name();
    my $cmd = qq{validation_acmg('$pname','99',this,'xxx','xxx')};
	  my $icon = $cgi->span({class=>"glyphicon glyphicon-eye-open  pull-left",'aria-hidden'=>"true"});
 	  my $icon_help = $cgi->span({class=>"glyphicon glyphicon-question-sign pull-left",'aria-hidden'=>"true"});
	my $icon_save = $cgi->span({class=>"glyphicon glyphicon glyphicon-copy pull-left",'aria-hidden'=>"true" });
	my $styleb = qq{style="font-size:12px"};
	my $stylec = qq{style="font-size:10px"};
	#'onclick=byPatientDude($pname,"high")';
	my $pname = $patient->name();
	print qq{
		<center>
		<div class="mydiv" style="width:99%">
			<table style="width:100%;">
			<td>
	};
	
	if ($ztime) {
		print qq{<div id="group_buttons_time_exec" class="btn-group-vertical  btn-sm clearfix pull-left" style="padding:20px;min-width:200px;">};
		
		my (@lButtons_times, @lButtons_times_id);
		foreach my $cat (sort split(' ', $ztime)) {
			my ($cat_name, $time) = split(':', $cat);
			my $b_id = 'b_time_exec_'.$cat_name.'_'.$patient->name();
			push(@lButtons_times_id, $b_id);
			$cat_name =~ s/_/ /;
			$time = sprintf("%.4f", $time);
			push(@lButtons_times, qq{<button type="button" id="$b_id" class="btn btn-secondary btn-xs" style="font-size:10px;display:none;"><span style="float:left">$cat_name</span> <span class="badge badge-alert" style="font-size: 10px;float:right;">$time </span></button>});
		}
		
		my $cmd_b_time = "view_time_exec_buttons('". join(';', @lButtons_times_id) ."')";
		print qq{<button type="button" onclick="$cmd_b_time" class="btn btn-secondary btn-xs" style="font-size: 10px;"><span style="float:left"><span class="glyphicon glyphicon-plus" aria-hidden="true"></span> Time Execution</span></span></button>};
		foreach my $b (@lButtons_times) {
			print qq{$b};
		}
		
		print qq{</div>};
	}
	
	my (@lbam_alamut, @lPatientsNames);
	foreach my $p (@{$patient->getFamily->getPatients()}) {
		push(@lbam_alamut, $p->bamUrl());
		push(@lPatientsNames, $p->name());
	}
	my $string_url_bam = join(',', @lbam_alamut);
	my $string_url_names = join(',', @lPatientsNames);
	my $project_name = $patient->project->name;
	print qq{		
			</td>
			<td>
				<section class="notepad" style='width:100%'>
					<div class="notepad-heading">
						<hh1>
							<div class="btn-group  btn-sm clearfix pull-right" style="position:relative;padding-block:15px;float:right;top:-10px">
								<div class=" btn btn-info btn-sm " aria-label="Left Align"  onClick='collapse_all("$string_panel","$string_label");'  >
									$icon
								</div>
								
								<div class=" btn btn-primary btn-sm " aria-label="Left Align" onClick='load_polyviewer_export_xls(1,2);'  style="background-color: #009B77;" >
									<b><i class="fa fa-file-excel-o pull-left"></i></b>
								</div>
									<div id="negative_$pname" class=" btn btn-danger btn-sm" aria-label="Left Align" onClick='validation_negative("$pname",99,this);' style="visibility:$hidden" >
									$icon_save&nbsp;<span>negative</span>
									
								</div>
							</div>	
							<div class="btn-group  btn-sm clearfix pull-right" style="color:black;position:relative;padding-block:15px;float:left;top:-10px;right:10%">
								<button type="button" class="btn btn-alert btn-xs" $styleb>Genes <span class="badge badge-success" $stylec>$nb_genes </span></button>	
								<button type="button" class="btn btn-alert btn-xs" $styleb>Variations <span class="badge" $stylec>$nb_var </span></button>
								<button type="button" class="btn btn-danger btn-xs" onclick='load_polyviewer(1,0,1)' $styleb >DM <span class="badge" $stylec>$nb_dm</span></button>
								<button type="button" class="btn  btn-xs" style="font-size:12px;margin-left:8px;margin-right:2px;">CNV</button>
								<button type="button" class="btn btn-info btn-xs" onclick='byPatientDude("$pname","1")' style="background-color:#BC243C;font-size:8px;margin:2px" >High</button>
								<button type="button" class="btn btn-info btn-xs" onclick='byPatientDude("$pname","2")' style="background-color:#EFC050;font-size:8px;margin:2px" >Med</button>
								<button type="button" class="btn btn-info btn-xs" onclick='byPatientDude("$pname","3")' style="background-color:#98B4D4;font-size:8px;margin:2px">Low</button>
								<button type="button" class="btn  btn-xs" style="font-size:12px;margin-left:8px;margin-right:2px;" onClick="dijit.byId('dialog_igv_alamut_download_infos').show();"><span class="glyphicon glyphicon-info-sign" style="color:blue;"></span> LOAD BAM</button>
								<button type="button" class="btn btn-info btn-xs" style="background-color:white;font-size:8px;" onClick='httpGetLoadOnlyListBam("$string_url_bam");'><img style="width:16px;height:16px;" src="images/polyicons/alamut_visual.png"></img></button>
							</div>
						</hh1>
		      		</div>
		      		$table
		 	 </section>
	 	 </td>
	 	 <td>
	};
	
	print qq{<div class="btn-group-vertical  btn-sm clearfix pull-right" style="color:black;padding:20px;width:200px;">};
	foreach my $name (sort keys %$version){
		my $v = $version->{$name}->{version};
		print qq{<button type="button" class="btn btn-success btn-xs" style="font-size: 10px;"><span style="float:left">$name</span> <span class="badge badge-alert" style="font-size: 10px;float:right;">$v </span></button>	};
	}
	foreach my $name (sort keys %$date){
		my $v = $date->{$name};
		print qq{<button type="button" class="btn btn-primary btn-xs" style="font-size: 10px;"><span style="float:left">$name</span> <span class="badge badge-alert" style="font-size: 10px;float:right;">$v </span></button>	};
	}
	print qq{</div>};
	
	print qq{ 	 
	 	 </td>
	 </table>
	 </div>
	 </center>
	};
	
	return;
}

sub printNavBar1 {
	my ($patient,$hgenes,$statistics) = @_;
	my $cgi = new CGI;
	my $project = $patient->project;
	my $all_panel=[];
	my $all_label=[];
	my $nb_var = $statistics->{variations};
	my $nb_dm = $statistics->{DM};
	my $nb_genes = $statistics->{genes};

	

foreach my $gene ( @$hgenes){
	$gene->{js_id}.= rand(500);
	push(@$all_panel,"panel_".$gene->{js_id});
	push(@$all_label,"label_".$gene->{js_id});
	$nb_genes++;
}

my $string_panel = join(";",@$all_panel);
my $string_label = join(";",@$all_label);
	
	

	
	#info filter 
	my $small_panel ;
	my $fam = $patient->getFamily();
	if ($fam->isTrio && $patient->isChild()){
		$small_panel = small_panel_trio($patient,$cgi);
	}
	else {
		$small_panel = small_panel($patient,$cgi);
	}
	my ($st, $text1);
	
	
	  my $icon = $cgi->span({class=>"glyphicon glyphicon-eye-open  pull-left",'aria-hidden'=>"true"});
 	  my $icon_help = $cgi->span({class=>"glyphicon glyphicon-question-sign pull-left",'aria-hidden'=>"true"});
  
  
  my $label  = "label-default";
		$label  = "label-danger"  if $cgi->param('denovo');
		my $text2 =  qq{&nbsp;<h4 style="position:relative;top:-10px;"float:right><span class="label $label">denovo</span>};
		my $text2 =  qq{<span class="label $label">denovo</span>};
		$label  = "label-default";
		$label  = "label-danger"   if $cgi->param('recessive');
		$text2 .=  qq{&nbsp;<span class="label $label">recessive</span>};
		$label  = "label-default";
		$label  = "label-danger"  if $cgi->param('xor');
		$text2 .=  qq{&nbsp; <span class="label $label">exclusive (M or F)</span>};
		$label  = "label-default";
		$label  = "label-danger"  if $cgi->param('both');
		$text2 .=  qq{&nbsp;<span class="label $label">both</span></h4>};
  
  	my $text_panel_name = "ALL";#$panel_name;
  	#if ($vectors_hex) { $text_panel_name = 'PolyQuery Filters'; }
  
	my $text = qq{<nav class="navbar navbar-default">
  <div class="container-fluid">
      	$small_panel
     
		<div class="btn-group  btn-sm clearfix pull-right" style="position:relative;padding-block:15px;float:right">
		
			<div class=" btn btn-info btn-sm " aria-label="Left Align"  onClick='collapse_all("$string_panel","$string_label");' >
				$icon
			</div>
			 		
			<div class=" btn btn-info btn-sm" aria-label="Left Align" onClick='dijit.byId("dialog_help_1").show();' >
				$icon_help
			</div>
			<div class=" btn btn-danger btn-sm " aria-label="Left Align" onClick='load_polyviewer_export_xls(1);' >
				<b><i class="fa fa-file-excel-o pull-left"></i></b>
			</div>
		</div>
		<div class="btn-group  btn-sm clearfix pull-right" style="position:relative;padding-block:15px;float:left">
			<button type="button" class="btn btn-alert btn-xs">Genes <span class="badge badge-success">$nb_genes </span></button>	
			<button type="button" class="btn btn-alert btn-xs">Variations <span class="badge">$nb_var </span></button>
			<button type="button" class="btn btn-danger btn-xs">DM <span class="badge">$nb_dm</span></button>
		</div>
	    <div style="position:relative;float:right;">
	    	<span style="position:relative;padding-top:7px;padding-right: 6px;color:white"><span style="color:orange;">Used Filters [ </span>$text1<span style="color:orange;"> ] </span></span>
	    	<br><br>
		</div>
  </div>
  </div>
  
</nav>
	};
	return;
}


sub cnv_select{
	my ($patient,$cnv_id,$cgi) = @_;
	$cgi =  new CGI unless $cgi;
	my $buffer = $patient->buffer;
	my $project = $patient->project;
	my $out;
	my $pname = $patient->{name};
	my $tt = $patient->name."_".$cnv_id;
	my $menu = $tt."_menu";
	my $sp = $tt."_span";
	my $tdid = $cnv_id."_".rand(500);
	my $bgcolor = "info";

my $saved  ;
my $all_validations = $patient->validations_cnv;
my $validation_term;
my $validation_value;
 if (exists $all_validations->{$cnv_id}){
 	$saved =  $all_validations->{$cnv_id};
 	$validation_term = $saved->[0]->{term};
 	$validation_value = $saved->[0]->{validation};
 }
 
# warn Dumper $all_validations;
 #die();
##$saved = $all_validations->{$val_id}->[0]->{validation}  if exists $all_validations->{$val_id};

my $option;

foreach my $val (sort {$b <=> $a }keys %{$buffer->value_validation}){
	my $term = $buffer->value_validation->{$val};
	my $sel ="";
	if ($validation_term eq $term ){
		$sel = "selected";
	}
	$option .=  qq{<option $sel value="$val">$term</option>\n};
}
unless ($validation_value){
	$option .=  qq{<option selected value="0">-</option>\n};
}

	$option .="</select></div>";
	$bgcolor = "info";
	$bgcolor = "secondary" if  $validation_value >= 2;
	$bgcolor = "warning" if  $validation_value >= 3;
 	$bgcolor = "danger" if  $validation_value >= 5;
 	
	my $uniq_id ="$pname"."+"."_$tdid";
	my $label_id= "div_$uniq_id";
	my $select_id =  "select_$uniq_id";
	my $force_text = qq{
		</select>
		<label id ="$label_id">
		<input type="checkbox" onChange="document.getElementById('select_$uniq_id').disabled = false;"> force
		</label>
		</div>
	};
	my $disabled ="";
	if ($saved){
		$disabled ="disabled";
	}
	else {
 		$force_text=qq{</select><label id ="$label_id"></label></div>};
	}
	my $select_text = qq{
		<div  style="vertical-align:middle;padding:1px;border-radius: 5px; -moz-border-radius: 5px; -webkit-border-radius: 5px; border: 2px solid #FFFFFF;" >
		<select id="select_$uniq_id" style="padding:1px" onchange ="validation_acmg_cnv('$pname','$cnv_id',this,'$uniq_id');" $disabled>
	};
	my $td_text = $select_text.$option.$force_text;
	$out .= $cgi->td({class=>$bgcolor,style=>"color:#000000;text-align: center;vertical-align:middle",id=>"td_$uniq_id"},$td_text); 
	return $out;
}


sub validation_select{
	my ($patient,$variation,$gene,$cgi) = @_;
	$cgi =  new CGI unless $cgi;
	my $buffer = $patient->buffer;
	my $project = $patient->project;
	my $out;
	my $pname = $patient->{name};
	my $vid = $variation->{id};
	my $tt = $patient->name."_".$gene->{js_id};
	my $menu = $tt."_menu";
	my $sp = $tt."_span";
	my $tdid = $gene->{js_id}."_".rand(50);
			
#	my @tds = ( qq{<a href="#" id="$bi_todo" class="button clean" onclick="this.className='button todo';validation('$pname','$vid',this,-3,'todo','$sp','$menu',$bi);">todo</a>},
#						qq{<a href="#"  id="$bi_ho" class="button clean" onclick="this.className='button todo';validation('$pname','$vid',this,1,'todo','$sp','$menu',$bi);">Ho</a>},
#						qq{<a href="#" id="$bi_he"  class="button clean" onclick="this.className='button todo';validation('$pname','$vid',this,2,'todo','$sp','$menu',$bi);">He</a>},
#						qq{<a href="#"  id="$bi_fp" class="button clean" onclick="this.className='button todo';validation('$pname','$vid',this,-1,'todo','$sp','$menu',$bi);">FP</a>}		
#	
my $bgcolor = "info";

my $val_id = $gene->{id}."!".$variation->{html}->{id};

my $saved =[] ;
my $all_validations = $patient->validations;
my $validation_term;
my $validation_value = 0;
 if (exists $all_validations->{$val_id}){
 #	my @found = grep {$_->{sample_id} eq $patient->id} @{$all_validations->{$val_id}};
 	$saved =  $all_validations->{$val_id};
 	$validation_term = $saved->[0]->{term};
 	$validation_value = $saved->[0]->{validation};
 }
 
#$saved = $all_validations->{$val_id}->[0]->{validation}  if exists $all_validations->{$val_id};

my $option;
foreach my $val (sort {$b <=> $a }keys %{$buffer->value_validation}){
	my $term = $buffer->value_validation->{$val};
	my $sel ="";

	if (lc($validation_term) eq lc($term) ){
		$sel = "selected";
	}
	$option .=  qq{<option $sel value="$val">$term</option>\n};
}
unless ($validation_term){
	$option .=  qq{<option selected value="0">-</option>\n};
}

	$option .="</select></div>";
	$bgcolor = "info";
	$bgcolor = "secondary" if  $validation_value >= 3;
	$bgcolor = "warning" if  $validation_value >= 4;
 	$bgcolor = "danger" if  $validation_value >= 5;
 	
	my $uniq_id ="$pname"."+"."_$tdid";
	my $label_id= "div_$uniq_id";
	my $select_id =  "select_$uniq_id";
	my $force_text = qq{
		</select>
		<label id ="$label_id">
		<input type="checkbox" onChange="document.getElementById('select_$uniq_id').disabled = false;"> force
		</label>
		</div>
	};
	my $disabled ="";
	if ($validation_term){
		$disabled ="disabled";
	}
	else {
 		$force_text=qq{</select><label id ="$label_id"></label></div>};
	}
	my $gene_id = $gene->{id};
	my $select_text = qq{
		<div  style="vertical-align:middle;padding:1px;border-radius: 5px; -moz-border-radius: 5px; -webkit-border-radius: 5px; border: 2px solid #FFFFFF;" >
		<select id="select_$uniq_id" style="padding:1px" onchange ="validation_acmg('$pname','$vid',this,'$uniq_id','$gene_id');" $disabled>
	};
	my $td_text = $select_text.$option.$force_text;
	$out .= $cgi->td({class=>$bgcolor,style=>"color:#000000;",id=>"td_$uniq_id"},$td_text); 
	$variation->{html}->{validation_select}->{$gene->{id}} = $out;
	#value_html($variation,"validation_select",$tt,$out);
	return $out;
}




sub table_cnv_genes_transcripts {
	my ($patient, $genes, $hGenes_dude) =@_;
	
	my $h_chr_cnv;
	my $project = $patient->getProject();
	#return "" unless $genes;
	$genes = [] unless $genes;
	my $header_cnv_genes;
	if ($patient->isGenome()) {
		$header_cnv_genes = ["<center><u>Name</u></center>","<center><u>Omim</u></center>","<center><u>pLI</u></center>","<center><u>Panels</u></center>","<center><u>Phenotypes</u></center>","<center><u>View VAR</u></center>"];
		#$header_cnv_genes = ["<center><u>Name</u></center>","<center><u>Omim</u></center>","<center><u>pLI</u></center>","<center><u>Panels</u></center>","<center><u>Phenotypes</u></center>","<center><u>View CNV</u></center>","<center><u>View VAR</u></center>"];
	}
	else {
		#$header_cnv_genes = ["<center><u>Name</u></center>","<center><u>Type</u></center>","<center><u>Nb Exons</u></center>","<center><u>Locus</u></center>","<center><u>Main Trasncripts<br>Coverage Graph</u></center>","<center><u>Phenotypes</u></center>","<center><u>Description</u></center>","<center><u>Omim</u></center>","<center><u>pLI</u></center>","<center><u>View CNV</u></center>","<center><u>View VAR</u></center>"];
		$header_cnv_genes = ["<center><u>Name</u></center>","<center><u>Type</u></center>","<center><u>Nb Exons</u></center>","<center><u>Locus</u></center>","<center><u>Phenotypes</u></center>","<center><u>Description</u></center>","<center><u>Omim</u></center>","<center><u>pLI</u></center>","<center><u>View CNV</u></center>","<center><u>View VAR</u></center>"];
	}
	my $level = 2;
	my $cgi          = new CGI();
	#my $value = $atr->[0]->{value}->{impact_score};
	my $color = "#EEE";
	#if ($value > 3){
		$color = "red";
	#}
	#if ($value >= 2){
		$color = "#FF8800";
	#}
	my $html;
	$html .= qq{<div style="max-height:300px;overflow-y:auto;border:solid 1px grey;">};
	$html .= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 8px;font-family:  Verdana;margin-bottom:0px;"});
	my @colors = ("#F9F6FF","#F9F9F9");
	
	my $nb =0;
	my $nb_skip = 0;
	$html .=  $cgi->start_Tr({style=>"position:sticky;z-index:9;top:0;background-color:grey;color:white;"}).$cgi->th({style=>"text-align: center;"}, $header_cnv_genes)."".$cgi->end_Tr();
	my $rids = [];
	my $nb_gene = scalar(@$genes);
	foreach my $gene (sort {$b->score <=> $a->score} @$genes) {
		$nb ++;
		my $hide;
		if ($nb > 3) { $hide = "display:none;"; }
		if ($gene->score >= 1) { $hide = undef; }
		$nb_skip ++ if $hide;
		my $c = $colors[$nb%2];
		my $rid = "row_gene_".$gene->{id}."_".time."_".$nb;
		push(@$rids,$rid) if $hide;
		#GENE NAME
		my $in = "";#$gene->{omim_inheritance};
			
		$in ="" if $in eq "-";
		$in = "X-linked " if $in =~/X-linked/;
		my $uc = qq{https://gnomad.broadinstitute.org/gene/}.$gene->name;
	 	my $oc = qq{onClick='window.open("$uc")'};	
		
		#TODO: here
		my $gene_id = $gene->{id};
		
		#PHENOTYPE
		my ($pheno,$nb_other_terms) =$gene->polyviewer_phentotypes();
		
	   	my $b_span_pheno = '';
   			if ($pheno) {
		   		if (length($pheno) > 70) { $pheno = substr($pheno, 0, 70).'...'; }
   				if ($nb_other_terms > 0) { $pheno .= " <span style='color:blue'>+ $nb_other_terms terms</span>"; }
   				$b_span_pheno .= qq{<a class="btn btn-xs" role="button" style="font-size:8px;background-color:#EEE;color:black;border:solid 1px black;" onclick="update_grid_gene_phenotypes(\'$gene_id\')";"><span>$pheno</span></a>};
   			}
		 
		 #OMIM
		 my $b_omim;
		 my $omim = $gene->omim_id();
		 if ($omim ne "") { $b_omim .=qq{<a class="btn btn-xs" style="font-size:8px;background-color:#EEE;color:black;border:solid 1px black;" href="http://www.omim.org/entry/$omim" role="button" target="_blank">Omim</a>}; }
		 else { $b_omim .=qq{<a class="btn btn-xs" style="font-size:8px;background-color:#EEE;color:black;border:solid 1px black;" role="button">-</a>} if $omim eq ""; }
		 
		 #PANELS
		 my $panel_name1 = join("-",keys %{$gene->{panels}});
		 my $hPanels;
		 foreach my $panel_name (keys %{$gene->{panels}}) {
   			my $pheno_name = $gene->{panels}->{$panel_name}->{phenotype};
   			$hPanels->{$pheno_name}->{$panel_name} = undef;
   		 }
		 my $panel_list;
   		 foreach my $pheno_name (sort keys %$hPanels) {
   			$panel_list .= "<b><u>Phenotype: ".$pheno_name."</b></u><br>";
   			foreach my $panel_name (sort keys %{$hPanels->{$pheno_name}}) {
   				$panel_list .= '  - '.$panel_name."<br>";
   			}
   			$panel_list .= "<br>";
   		 }
   		 my $nb_panels = scalar(keys %{$gene->{panels}});
   		 my $b_panels;
		 if ($nb_panels > 0) { $b_panels .=qq{<a class="btn btn-xs" style="font-size:8px;background-color:#EEE;color:black;border:solid 1px black;" role="button" onclick="document.getElementById('span_list_panels').innerHTML='$panel_list';dijit.byId('dialog_list_panels').show();">$nb_panels</a>}; }
		 else { $b_panels .=qq{<a class="btn btn-xs" style="font-size:8px;background-color:#EEE;color:black;border:solid 1px black;" role="button" onclick="document.getElementById('span_list_panels').innerHTML='$panel_list';dijit.byId('dialog_list_panels').show();">-</a>}; }
		 
		 #GTEX
		my ($gid,$t) = split("_",$gene_id);
		my $b_gtex = qq{<a class="btn btn-xs" href="https://gtexportal.org/home/gene/$gid" role="button" target="_blank" style="font-size:8px;background-color:#EEE;color:black;border:solid 1px black;">Gtex</a>};
		 
		 #PLI
		my $pli = $gene->{pLI}*1.0;
		my $type ="green";
		$type = "orange" if $pli >= 0.75;
		$type = "red" if $pli >= 0.9;
		
 		my $b_pli =qq{<a class="btn btn-xs" role="button" href='https://gnomad.broadinstitute.org/gene/$gid' target='_blank' style="background-color:#EEE;color:black;border:solid 1px black;"><span style="font-size:8px;background-color:#EEE;">$pli</span></a>};
 		
 		#VARIANTS
		my ($nb_dude, $nb_var);
		$nb_var = 'View';
		
 		my $patient_name = $patient->name();
 		my $project_name = $project->name();
 		my $cmd_var = qq{dijit.byId('dialog_hgmd').show();view_var_from_proj_gene_pat('$project_name', '$gene_id', '$patient_name', '', 'all', 'nocnv');};
 		my $b_var = qq{<a class="btn btn-xs" role="button" onclick="$cmd_var" style="background-color:#3AB795;color:white;border:solid 1px black;" disabled><span style="font-size:8px;" class="glyphicon glyphicon-ban-circle"></span></a>};
		
		my $level_dude_text;
		my $level_dude = 'Not';
		if (exists $hGenes_dude->{$gene->external_name()}->{'high'}->{'all'} or exists $hGenes_dude->{$gene_id}->{'high'}->{'all'}) {
			$level_dude = 'High';
			$level_dude_text = "high";
		}
		elsif (exists $hGenes_dude->{$gene->external_name()}->{'medium'}->{'all'} or exists $hGenes_dude->{$gene_id}->{'medium'}->{'all'}) {
			$level_dude = 'Med';
			$level_dude_text = "medium";
		}
		elsif (exists $hGenes_dude->{$gene->external_name()}->{'low'}->{'all'} or exists $hGenes_dude->{$gene_id}->{'low'}->{'all'}) {
			$level_dude = 'Low';
			$level_dude_text = "low";
		}
		
		$html .=  $cgi->start_Tr({id=>$rid,style=>"border: 1px solid;background-color:$c;".$hide});
		my $bname = printButton($gene->score,[1,4],$gene->external_name.' <b><u>'.$in.'</b></u>',$oc) ;
		$gene->{name} = $bname;
		$html.= $cgi->td("<center>".$bname."</center>");
		
		my ($noise, $del_ho);
		if (exists $hGenes_dude->{$gene->id()}) {
			$noise = int($hGenes_dude->{$gene->id()}->{$level_dude_text}->{perc_noise});
			$del_ho = $hGenes_dude->{$gene->id()}->{$level_dude_text}->{del_ho};
		}
		elsif (exists $hGenes_dude->{$gene->external_name()}) {
			$noise = int($hGenes_dude->{$gene->external_name()}->{$level_dude_text}->{perc_noise});
			$del_ho = $hGenes_dude->{$gene->external_name()}->{$level_dude_text}->{del_ho};
		}
		
		my $level_dude_noise = '';
		if ($noise == 0) {
			$level_dude_noise = qq{<button type="button" class="btn btn-xs" style="background-color:#BC243C;color:white;border:solid 1px black;font-size:8px;">0%</button>};
		}
		else {
			$level_dude_noise = qq{<button type="button" class="btn btn-xs" style="background-color:white;color:black;border:solid 1px black;font-size:8px;">$noise%</button>};
		}

		my $color_dude = 'grey';
		$color_dude = '#BC243C' if (lc($level_dude) eq 'high');
		$color_dude = '#EFC050' if (lc($level_dude) eq 'med');
		$color_dude = '#98B4D4' if (lc($level_dude) eq 'low');
		
		my $b_cnv = qq{<button type="button" class="btn btn-xs" onclick="zoomDude(\'$project_name\','$gene_id','$patient_name', 'force_visualisation')" style="background-color:$color_dude;color:white;border:solid 1px black;font-size:8px;">$level_dude</button>};
 		if ($nb_var eq 'View' or $nb_var > 0)  { $b_var =qq{<a class="btn btn-xs" role="button" onclick="$cmd_var" style="background-color:#3AB795;color:white;border:solid 1px black;"><span style="font-size:8px;">$nb_var</span></a>}; }
		
		# LOCUS (EXOMES)
 		my $chr_id = $gene->getChromosome->id();
 		my $gene_start = $gene->start();
 		my $gene_end = $gene->end();
		my $locus = $chr_id.":".$gene_start."-".$gene_end;
		my $ucsc = qq{https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=$locus};
		my $b_locus = qq{<button class= "btn btn-xs btn-primary " style="background-color: #EEE;font-size: 1em;font-family:  Verdana;color:black"><a href=$ucsc target='_blank' style="color:black;">}.$locus.qq{</a></button>};
		
		if ($patient->isGenome) {
			$html.= $cgi->td("<center>".$b_omim."</center>");
			#$html.= $cgi->td("<center>".$b_gtex."</center>");
			$html.= $cgi->td("<center>".$b_pli."</center>");
			$html.= $cgi->td("<center>".$b_panels."</center>");
			$html.= $cgi->td("<center>".$b_span_pheno."</center>");
#			$html.= $cgi->td("<center>".$b_cnv."</center>");
			$html.= $cgi->td("<center>".$b_var."</center>");
		}
		else {
			my ($nb_dup, $nb_del, $nb_all);
			if (exists $hGenes_dude->{$gene->external_name()} and exists $hGenes_dude->{$gene->external_name()}->{$level_dude_text}->{dup}) {
				$nb_dup = $hGenes_dude->{$gene->external_name()}->{$level_dude_text}->{dup};
				$nb_del = $hGenes_dude->{$gene->external_name()}->{$level_dude_text}->{del};
				$nb_all = $hGenes_dude->{$gene->external_name()}->{$level_dude_text}->{all};
			}
			elsif (exists $hGenes_dude->{$gene_id} and exists $hGenes_dude->{$gene_id}->{$level_dude_text}->{dup}) {
				$nb_dup = $hGenes_dude->{$gene_id}->{$level_dude_text}->{dup};
				$nb_del = $hGenes_dude->{$gene_id}->{$level_dude_text}->{del};
				$nb_all = $hGenes_dude->{$gene_id}->{$level_dude_text}->{all};
			}
			else {
				warn "\n\n\n\n";
				warn $gene_id;
				warn Dumper $hGenes_dude->{$gene_id};
				confess();
			}
			
			#TYPE
	 		my ($b_type, $b_exons);
	 		$nb_del = $nb_del + $del_ho;
			if ($nb_dup > $nb_del) {
				$b_type = qq{<button type="button" class= "btn btn-xs btn-primary " style="background-color: #217DBB;font-size: 1em;font-family:Verdana;color:white">gain</button>};
				$b_exons = qq{<button type="button" class= "btn btn-xs btn-primary " style="background-color: #217DBB;font-size: 1em;font-family:Verdana;color:white">$nb_dup/$nb_all</button>};
	 		}
	 		elsif ($nb_dup < $nb_del) {
	 			$b_type = qq{<button type="button" class= "btn btn-xs btn-primary " style="background-color: #E74C3C;font-size: 1em;font-family:Verdana;color:white">del</button>};
	 			if ($del_ho > 0) {
					$b_type = qq{<button type="button" class="btn btn-xs" style="background-color:#E74C3C;color:white;border:solid 1px black;font-size:8px;">del HO</button>};
				}
	 			$b_exons = qq{<button type="button" class= "btn btn-xs btn-primary " style="background-color: #E74C3C;font-size: 1em;font-family:Verdana;color:white">$nb_del/$nb_all</button>};
	 		}
	 		else {
	 			my $nb_dup_del = $nb_dup + $nb_del;
				$b_type = qq{<button type="button" class= "btn btn-xs btn-primary " style="background-color: grey;font-size: 1em;font-family:Verdana;color:white">gain/del ?</button>};
				$b_exons = qq{<button type="button" class= "btn btn-xs btn-primary " style="background-color: grey;font-size: 1em;font-family:Verdana;color:white">$nb_dup_del/$nb_all</button>};
	 		}
	 		
			my @lbam_alamut;
			foreach my $p (@{$patient->getFamily->getPatients()}) {
				push(@lbam_alamut, 'https://www.polyweb.fr/'.$p->bamUrl());
			}
			my $string_url_bam = join(',', @lbam_alamut);
	 		my $b_igv = qq{<button type="button" class="btn btn-info btn-xs" style="background-color:white;font-size:8px;" onClick="displayInIGV('$chr_id', '$gene_start', '$gene_end');"><img style="width:16px;height:16px;" src="images/polyicons/igv_logo_32.png"></img></button>};
			my $b_alamut = qq{<button type="button" class="btn btn-danger btn-xs" style="background-color:white;font-size:8px;" onClick="displayLocusInAlamut('$chr_id', '$gene_start', '$gene_end');"><img style="width:16px;height:16px;" src="images/polyicons/alamut_visual.png"></img></button>};
			
#			my $i_tr = 0;
#			my $tab_graph;
#			foreach my $tr (@{$gene->getMainTranscripts()}) {
#				$i_tr++;
#               my ($tr_id, $chr_id) = split('_', $tr->id());
#				$tab_graph .= qq{<button class= "btn btn-primary btn-xs" style="background-color:#EEE;font-size: 1em;font-family: Verdana;color:black" onClick="load_graph_transcript('$patient_name', '$tr_id');">$tr_id</button>};
#			}
	 		
			$html.= $cgi->td("<center>".$b_type."</center>");
			$html.= $cgi->td("<center>".$b_exons."</center>");
			$html.= $cgi->td("<center>".$b_locus."</center>");
#			$html.= $cgi->td("<center>".$tab_graph."</center>");
			$html.= $cgi->td("<center>".$b_span_pheno."</center>");
			$html.= $cgi->td("<center>".$gene->description()."</center>");
			#$html.= $cgi->td("<center>".$b_panels."</center>");
			$html.= $cgi->td("<center>".$b_omim."</center>");
			#$html.= $cgi->td("<center>".$b_gtex."</center>");
			$html.= $cgi->td("<center>".$b_pli."</center>");
			$html.= $cgi->td("<center>".$b_cnv."</center>");
			$html.= $cgi->td("<center>".$b_var."</center>");
		}
		$html.= $cgi->end_Tr();
		
		#last if ($nb == 500);
	}
#	if ($nb == 500 and $nb_gene >= 500) {
#		$html .=  $cgi->start_Tr();
#		$html.= $cgi->td({style=>"box-shadow: 1px 1px 2px #555;background-color:#CECFCE;",colspan=>scalar(@$header_cnv_genes)},qq{<span class="glyphicon glyphicon-plus"></span> }."LIMITED to 500 Genes but found $nb_gene");
#		$html.= $cgi->end_Tr();
#	 	$html.=$cgi->end_table();
#		return $html
#	}
	if ($nb_skip == 0){
	 	$html.=$cgi->end_table();
		return $html
	}
	my $js = encode_json $rids;
	my $za = "hide_tr_".time."_".int(rand(50000));
	$html .=  $cgi->start_Tr({id=>$za});
	$html.= $cgi->td({style=>"box-shadow: 1px 1px 2px #555;background-color:#CECFCE;",colspan=>scalar(@$header_cnv_genes),onClick=>qq{showTranscripts($js,"$za");}},qq{<span class="glyphicon glyphicon-plus"></span> }."view $nb_skip Genes");
	$html.= $cgi->end_Tr();
	
	$html.=$cgi->end_table();
	$html.=qq{</div>};
	#$atr->{html} = $html;
	return $html;
	
}



sub get_hash_genes_dude {
	my ($patient, $by_names_or_ids, $list_type, $panel,$noprint) = @_;
	my ($h_type_levels, $h_panels_tr);
	#unless ($list_type) { push(@$list_type, 'high'); }
	#push(@$list_type, 'medium');
	#die();
	foreach my $type (@$list_type) {
		$h_type_levels->{$type} = undef;
	}
	#$h_type_levels->{medium} = 1;
	if ($panel) { $h_panels_tr = $panel->genes_id(); }

	$by_names_or_ids = 'names' unless ($by_names_or_ids);
	my $hGenes_dude = {};
	my $no = $patient->getTranscriptsDude("r");
	my @lPatients = @{$patient->getProject->getPatients()};
 	my @selected_patients;
 	push(@selected_patients, $patient);
 	my $nb_genes_done = 0;
	foreach my $type (@$list_type) {
		my $list_tr_high;
		eval { $list_tr_high = $no->get($type); };
		if($@) {}
#		warn $type;
#		warn Dumper $list_tr_high;
		#ENST00000631057
#		my ($toto) = grep {$_ =~/ENST/ } @$list_tr_high;
	#	warn Dumper $list_tr_high;
	#	die();
		unless ($list_tr_high) {
			$no->close();
			next;
		}
		if (scalar(@$list_tr_high) >= 250 and scalar keys %$h_type_levels == 1) {
			$hGenes_dude->{too_much} = 1;
			return $hGenes_dude;
		}
		if ($list_tr_high && scalar(@$list_tr_high) > 5000) {
			$hGenes_dude->{too_much} = 1;
			return $hGenes_dude;
		}
		my (@ltr, $hGenesToDo, $nbg);
		foreach my $t (@{$patient->getProject->newTranscripts($list_tr_high)}) {
			my $gene_external_name = $t->gene_external_name;
#			die() if $gene_external_name eq "GPC3";
			unless (exists $hGenesToDo->{$gene_external_name}) {
				$hGenesToDo->{$gene_external_name} = undef;
				push(@ltr, $t);
				$nbg++;
			}
			last if ($nbg == 250 and not exists $h_type_levels->{'medium'} and not $h_type_levels->{'low'})
		}
		my $fork = 5;
		$fork = 10 if (scalar (@ltr)>100);
		$fork = 15 if (scalar (@ltr)>150);
		$fork = 200 if (scalar (@ltr)>300);
		my $nb = int(scalar(@ltr)/($fork*2))+1;
		my $pm = new Parallel::ForkManager($fork);
		my $iter = natatime $nb, @ltr;
		$pm->run_on_finish(
		    sub { 
		    	my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$h)=@_;
		    	unless (defined($h) or $exit_code > 0) {
					print qq|No message received from child process $exit_code $pid!\n|;
					return;
				}
				foreach my $g (keys %{$h}){
#					warn $g;
					next if ($g eq 'not');
					foreach my $t (keys %{$h->{$g}}){
						next unless (exists $h_type_levels->{$t});
						$hGenes_dude->{$g}->{$t}->{dup} = $h->{$g}->{$t}->{dup};
						$hGenes_dude->{$g}->{$t}->{del} = $h->{$g}->{$t}->{del};
						$hGenes_dude->{$g}->{$t}->{del_ho} += $h->{$g}->{$t}->{del_ho};
						$hGenes_dude->{$g}->{$t}->{all} = $h->{$g}->{$t}->{all};
						$hGenes_dude->{$g}->{$t}->{all_others} = $h->{$g}->{$t}->{all_others};
						$hGenes_dude->{$g}->{$t}->{grey_others} = $h->{$g}->{$t}->{grey_others};
						$hGenes_dude->{$g}->{$t}->{noise} = $h->{$g}->{$t}->{noise};
						$hGenes_dude->{$g}->{$t}->{perc_noise} = $h->{$g}->{$t}->{perc_noise};
						$hGenes_dude->{$g}->{$t}->{perc_grey} = $h->{$g}->{$t}->{perc_grey};
					}
				}
		    }
	    );
		
		$patient->getProject->buffer->dbh_deconnect();
		my $xp = 0;
 	 	while( my @ltr_tmp = $iter->() ){
 	 		my $pid = $pm->start and next;
			$patient->getProject->buffer->dbh_reconnect();
			my $hGenes_dude_tmp;
			$hGenes_dude_tmp->{'not'}++;
		 	foreach my $t (@ltr_tmp){
		 		$xp++;
		 		unless ($noprint){
				print "." if ($xp % 20 == 0 and $patient->getProject->cgi_object());
				$patient->getProject->print_dot(50);
		 		}
		 		my $g_name_id;
		 		$g_name_id = $t->gene_external_name() if ($by_names_or_ids eq 'names');
		 		$g_name_id = $t->getGene->id() if ($by_names_or_ids eq 'ids');
		 		
		 		next if (exists $hGenes_dude->{$g_name_id});
		 		
	 	 		# FILTRE TRANSCRIPTS HIGH DUDE
				my $coverage = polyweb_dude->new(patients=>\@lPatients,transcript=>$t,limit=>undef,selected_patients=>\@selected_patients);
				$coverage->init_matrices();
				my $hcov = $coverage->quality;
				
				$hGenes_dude_tmp->{$g_name_id}->{$type}->{dup} += $hcov->{$patient->name()}->{dup};
				$hGenes_dude_tmp->{$g_name_id}->{$type}->{del} += $hcov->{$patient->name()}->{del};
				$hGenes_dude_tmp->{$g_name_id}->{$type}->{del_ho} += $hcov->{$patient->name()}->{del_ho};
				$hGenes_dude_tmp->{$g_name_id}->{$type}->{all} += $hcov->{$patient->name()}->{all};
				foreach my $other_patient (@lPatients) {
					next if ($other_patient->getFamily->name() eq $patient->getFamily->name());
					$hGenes_dude_tmp->{$g_name_id}->{$type}->{grey_others} += $hcov->{$other_patient->name()}->{grey};
					$hGenes_dude_tmp->{$g_name_id}->{$type}->{noise} += $hcov->{$other_patient->name()}->{dup};
					$hGenes_dude_tmp->{$g_name_id}->{$type}->{noise} += $hcov->{$other_patient->name()}->{del};
					$hGenes_dude_tmp->{$g_name_id}->{$type}->{all_others} += $hcov->{$other_patient->name()}->{all};
				}
				
				
			}
 	 		$pm->finish(0,$hGenes_dude_tmp);
 	 	}
		$pm->wait_all_children();
		$patient->getProject->buffer->dbh_reconnect();
	}
	$no->close();
	if (not $hGenes_dude  or scalar keys %{$hGenes_dude} == 0) {
		$hGenes_dude->{no_result} = 1;
		return $hGenes_dude;
	}
	#warn Dumper $hGenes_dude->{ENSG00000147257_X};
	#warn Dumper $hGenes_dude;
	#die();
	unless ($patient->isGenome()) {
			foreach my $g_name_id (keys %{$hGenes_dude}) {
				foreach my $type (keys %{$hGenes_dude->{$g_name_id}}) {
					#next unless (exists $hGenes_dude->{$g_name_id}->{$type});
					my $nb_grey_others = $hGenes_dude->{$g_name_id}->{$type}->{grey_others};
					my $nb_all_others = $hGenes_dude->{$g_name_id}->{$type}->{all_others};
					my $perc_grey = 0;
					$perc_grey = ($nb_grey_others / $nb_all_others) * 100 if ($nb_all_others > 0);
					my $nb_noise = $hGenes_dude->{$g_name_id}->{$type}->{noise};
					my $perc = 0;
					$perc = ($nb_noise / $nb_all_others) * 100 if ($nb_all_others > 0);
					$hGenes_dude->{$g_name_id}->{$type}->{'perc_noise'} = $perc;
					$hGenes_dude->{$g_name_id}->{$type}->{'perc_grey'} = $perc_grey;
				}
				if (exists $hGenes_dude->{$g_name_id}->{'high'} and $hGenes_dude->{$g_name_id}->{'high'}->{noise} > 0) {
					my $perc_noise = $hGenes_dude->{$g_name_id}->{'high'}->{'perc_noise'};
					my $nb_del = $hGenes_dude->{$g_name_id}->{'high'}->{del};
					my $nb_del_ho = $hGenes_dude->{$g_name_id}->{'high'}->{del_ho};
					$nb_del += $nb_del_ho;
					my $nb_dup = $hGenes_dude->{$g_name_id}->{'high'}->{dup};
					my $nb_all = $hGenes_dude->{$g_name_id}->{'high'}->{all};
					my $perc_gene_dup = ($nb_dup / $nb_all) * 100;
					my $perc_gene_del = ($nb_del / $nb_all) * 100;
					
					#CAS majorite exons concernes mais 
					if ($nb_del_ho) {
						
					}
					elsif ($nb_dup == $nb_del) {
						if ($perc_noise >= 20) {
							$hGenes_dude->{$g_name_id}->{'low'} = $hGenes_dude->{$g_name_id}->{'high'} if (exists $h_type_levels->{'low'});
							delete $hGenes_dude->{$g_name_id}->{'high'};
						}
						else {
							$hGenes_dude->{$g_name_id}->{'medium'} = $hGenes_dude->{$g_name_id}->{'high'} if (exists $h_type_levels->{'medium'});
							delete $hGenes_dude->{$g_name_id}->{'high'};
						}
					}
					#CAS majorite exons concernes mais 
					elsif ($perc_gene_dup >= 60 or $perc_gene_del >= 60) {
						if ($perc_noise >= 20) {
							$hGenes_dude->{$g_name_id}->{'medium'} = $hGenes_dude->{$g_name_id}->{'high'} if (exists $h_type_levels->{'medium'});
							delete $hGenes_dude->{$g_name_id}->{'high'};
						}
					}
					#CAS baisse si BCP TROP bruit
					elsif ($perc_noise >= 40) {
						$hGenes_dude->{$g_name_id}->{'low'} = $hGenes_dude->{$g_name_id}->{'high'} if (exists $h_type_levels->{'low'});
						delete $hGenes_dude->{$g_name_id}->{'high'};
					}
					#CAS baisse si TROP bruit avec DEL HO
					elsif ($perc_noise >= 10 and $nb_del_ho > 0) {
						$hGenes_dude->{$g_name_id}->{'medium'} = $hGenes_dude->{$g_name_id}->{'high'} if (exists $h_type_levels->{'medium'});
						delete $hGenes_dude->{$g_name_id}->{'high'};
					}
					#CAS baisse si TROP bruit sans DEL HO
					elsif ($perc_noise >= 5 and $nb_del_ho == 0) {
						$hGenes_dude->{$g_name_id}->{'medium'} = $hGenes_dude->{$g_name_id}->{'high'} if (exists $h_type_levels->{'medium'});
						delete $hGenes_dude->{$g_name_id}->{'high'};
					}
				}
				
				#CAS rattrapage MEDIUM DEL HO et peu NOISE
				elsif (exists $hGenes_dude->{$g_name_id}->{'medium'}) {
				
					my $nb_del = $hGenes_dude->{$g_name_id}->{'medium'}->{del};
					my $nb_del_ho = $hGenes_dude->{$g_name_id}->{'medium'}->{del_ho};
					$nb_del += $nb_del_ho;
					my $nb_dup = $hGenes_dude->{$g_name_id}->{'medium'}->{dup};
					my $perc_noise = $hGenes_dude->{$g_name_id}->{'medium'}->{perc_noise};
					if ($nb_del_ho > 0) {
					#if ($nb_dup != $nb_del and $perc_noise < 10 and $nb_del_ho > 0) {
						$hGenes_dude->{$g_name_id}->{'high'} = $hGenes_dude->{$g_name_id}->{'medium'};
						delete $hGenes_dude->{$g_name_id}->{'medium'};
					}
					elsif ($nb_dup != $nb_del and $perc_noise < 1) {
						$hGenes_dude->{$g_name_id}->{'high'} = $hGenes_dude->{$g_name_id}->{'medium'};
						delete $hGenes_dude->{$g_name_id}->{'medium'};
					}
				}
				
				#CAS rattrapage LOW  et peu NOISE
				elsif (exists $hGenes_dude->{$g_name_id}->{'low'}) {
					if ($hGenes_dude->{$g_name_id}->{'low'}->{perc_noise} < 1) {
						$hGenes_dude->{$g_name_id}->{'high'} = $hGenes_dude->{$g_name_id}->{'low'};
						delete $hGenes_dude->{$g_name_id}->{'low'};
					}
					elsif ($hGenes_dude->{$g_name_id}->{'low'}->{perc_noise} < 2) {
						$hGenes_dude->{$g_name_id}->{'medium'} = $hGenes_dude->{$g_name_id}->{'low'};
						delete $hGenes_dude->{$g_name_id}->{'low'};
					}
				}
				delete $hGenes_dude->{$g_name_id}->{'low'} unless (exists $h_type_levels->{'low'});
				delete $hGenes_dude->{$g_name_id}->{'medium'} unless (exists $h_type_levels->{'medium'});
				delete $hGenes_dude->{$g_name_id} if (scalar keys %{$hGenes_dude->{$g_name_id}} == 0);
			}
 	 }
	return $hGenes_dude;
}

sub print_dude_button {
		my ($patient) = @_;
		my $project_name = $patient->project->name;
		my $patient_name = $patient->name;
		my $cmd = qq{byPatientDude(\'$patient_name\')};
		
		return qq{<a type="button" class= "btn btn-xs  btn-primary"  role="button" target="_blank" onClick="$cmd"> <img src="https://img.icons8.com/color/30/000000/biotech.png">&nbspView CNV&nbsp</a>};
}




sub print_hotspot {
	my ($patient, $panel) = @_;
	
	my $hotspots = $patient->hotspot;
	return "" unless $hotspots;
	
	my @ltmp = split('/', $patient->getCapture->hotspots_filename);
	my $fileName = $ltmp[-1];
	
	my $project = $patient->project;
	my $pname = $patient->name();
	my $cgi          = new CGI();
	print $cgi->start_div({class=>"panel panel-warning ",style=>"border-color:white;-webkit-border-radius: 3px;-moz-border-radius: 3px;border-radius: 3px;border: 1px solid black;max-height:400px;overflow-y:auto;" });
	print $cgi->div({class=>"panel-heading"},qq{HOTSPOT - $fileName});
	my $div_alert;	
	my $s_id = $patient->{name};

my $t = time;

	my $out;

#$out .=  $cgi->start_div({class=>"panel-heading panel-warning warning ",style=>" min-height:13px;max-height:13px;padding:1px;border:1px"});
#	$out .= qq{<div class="btn  btn-success btn-xs " style="position:relative;bottom:1px;min-width:150px;" onClick='collapse("$panel_id","$label_id")'>  <span id= "$label_id" class="glyphicon glyphicon-triangle-right  "   style="float:left;"></span> $text &nbsp</div>};
	   		#	$out .=$cgi->span({class=>"label label-success"},qq{<span class='badge badge-primary badge-xs'  >$nbv</span>});
#		my $nbv = scalar (keys %{$hotspot->{results}->{$s_id}});
#		$out .=$cgi->span({class=>"label label-success"},qq{<span class='badge badge-primary badge-xs'  >$nbv</span>});	
	my @header = ("ID","NAME","PROT","A","C","G","T","DEL","INS","COV");	 
	$out .= $cgi->start_table({class=>"table table-striped table-condensed table-bordered table-hover table-mybordered",style=>"font-size: 9px;font-family:  Verdana;"});
foreach my $g (keys %$hotspots){
	#$out .=  $cgi->start_div({class=>"panel panel-info" });
	 #panel heading
	 
	 # $out.= $cgi->end_div();
		#REF	POS	COV	A	C	G	T	DEL	REFSKIP	SAMPLE
	#my $var_obj = $self->cache_lmdb_variations->get($vid);
	#  panel table
	$out.= $cgi->start_Tr();
	$out.=$cgi->th({colspan=>(scalar(@header)+1),style=>"background-color:#217DBB;color:white;font-size:12px"},$g);
	$out.= $cgi->end_Tr();
	$out.= $cgi->start_Tr();
	$out.=$cgi->th({style=>"background:#E0E0FF"},["igv",@header]);
	$out.= $cgi->end_Tr();
	
	my @bams;
	my @names;
	foreach my $p (@{$patient->getFamily->getPatients()}){
		push(@bams,$patient->bamUrl);
		push(@names,$patient->name());
	}
					
	my $f =  join(";",@bams);#$patient->{obj}->bamUrl;;
	 my $pnames = join(";",@names);
	foreach my $hotspot (@{$hotspots->{$g}}){
		my @td;
		my $chr = $project->getChromosome($hotspot->{REF});
		#my $var_obj = $chr->cache_lmdb_variations->get($hotspot->{GENBO_ID});
		my $var_obj = $project->_newVariant($hotspot->{GENBO_ID});
		#warn $chr->cache_lmdb_variations->get(0);
		my $style ={};
		 $style = {style=>"background-color:#D2386C;color:white"} if $var_obj && defined $var_obj->vector_id() && $var_obj->existsPatient($patient);
		 $out.= $cgi->start_Tr($style);
		 my $nba;
		 my $chrn = $chr->name;
		 my $start = $hotspot->{POS};
		 my $l = $chr->name.":".$start;
		 my $gn = $project->getVersion();
		 my $project_name = $project->name; 
		 my $v1 = "?/?";#.$hvariation->{allele};	
		# launch_web_igv_js
		my $text =qq{<button class='igvIcon2' onclick='launch_web_igv_js("$project_name","$pnames","$f","$l","$v1","$gn")' style="color:black"></button>};
		#my $text =qq{<button dojoType="dijit.form.Button"   iconClass='igvIcon' onclick='view_web_igv_bam("dialog_igv", "div_igv", "$l", "$f", "$pnames")' style="color:black"></button>};
		 
		 
		$out.=$cgi->td($text);
		foreach my $h (@header){
			if ($h eq  $hotspot->{A_ALT}){
				my $pc = int(($hotspot->{$h}/$hotspot->{COV})*1000)/10;
				my $color = "#f2dedc";
				$color = "#F7BFB9" if $pc>2;
				$color = "#E9897E" if $pc>5;
				$out.=$cgi->td({style=>"background-color:$color;color:black"},"$pc% (".$hotspot->{$h}.")");
			}
			elsif ($h eq  $hotspot->{A_REF}){
				my $pc = int(($hotspot->{$h}/$hotspot->{COV})*1000)/10;
				$out.=$cgi->td({style=>"background-color:#c7eadd;color:black"},"$pc% (".$hotspot->{$h}.")");
			}
			elsif ($h eq  "DEL" && $hotspot->{A_ALT} eq "-"){
				my $pc = int(($hotspot->{$h}/$hotspot->{COV})*1000)/10;
				my $color = "#f2dedc";
				$color = "#F7BFB9" if $pc>2;
				$color = "#E9897E" if $pc>5;
				$out.=$cgi->td({style=>"background-color:#F7BFB9;color:black"},"$pc% (".$hotspot->{$h}.")");
			}
			else {
				my $pc = int(($hotspot->{$h}/$hotspot->{COV})*1000)/10;
				$out .= $cgi->td($hotspot->{$h});
			}
			#push(@td, $hotspot->{$h});
		}
	
		#$out.=$cgi->td(\@td);	
		$out.= $cgi->end_Tr();
	}
	
	}
	$out.= $cgi->end_table();
	print $out;
	print "</div></div>";
	return;
}



sub print_cnv_exome {
	my ($patient, $types, $panel) = @_;
	if ($patient->isParent()) {
		if ($patient->getProject->validation_db() and $patient->getProject->validation_db() eq 'defidiag') {
			return;
		}
	}
	my $pname = $patient->name();
	my $cgi          = new CGI();
	print $cgi->start_div({class=>"panel panel-danger ",style=>"border-color:white;-webkit-border-radius: 3px;-moz-border-radius: 3px;border-radius: 3px;border: 1px solid black;max-height:400px;overflow-y:auto;" });
	if ($patient->getProject->isDiagnostic()) {
		print $cgi->div({class=>"panel-heading"},print_dude_button($patient));
	}
	else {
		print $cgi->div({class=>"panel-heading"},qq{Dude - HIGH Most Probable CNV Event <div class="btn-group  btn-sm clearfix " style="color:black;position:relative;"><button type="button" class="btn btn-info btn-xs" onclick="load_polyviewer(1, 'high')" style="background-color:#BC243C;font-size:8px;margin:2px">High</button><button type="button" class="btn btn-info btn-xs" onclick="load_polyviewer(1, 'high,medium')" style="background-color:#EFC050;font-size:8px;margin:2px">Med</button><button type="button" class="btn btn-info btn-xs" onclick="load_polyviewer(1, 'high,medium,low')" style="background-color:#98B4D4;font-size:8px;margin:2px">Low</button></div>});
	}
	my $hTypes;
	if ($types) {
		foreach my $type (split(',', $types)) { $hTypes->{$type} = undef; };
	}
	
#	$patient->getProject->cgi_object(1);
#	print "<div hidden>";
#	my $use_type = 'high';
#	$use_type = 'medium' if (exists $hTypes->{'medium'});
#	$use_type = 'low' if (exists $hTypes->{'low'});
#	my $hGenes_dude = get_hash_genes_dude_NEW($patient, $use_type);
#	print "</div>";
	
	my @ltypes = ('high');
	push(@ltypes, 'medium') if (exists $hTypes->{'medium'});
	push(@ltypes, 'low') if (exists $hTypes->{'low'});
	print qq{<div style="display: none">};
	$patient->getProject->cgi_object(1);
	my $hGenes_dude = get_hash_genes_dude($patient, 'ids', \@ltypes, '');
	print qq{</div>};
	my @lGenes_ids = sort keys %{$hGenes_dude};
	if ($lGenes_ids[0] eq 'too_much' or $lGenes_ids[-1] eq 'too_much') {
		print $cgi->start_table({class=>"table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 10px;font-family:  Verdana;", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
		print $cgi->td({style=>"text-align: center;vertical-align:middle"},'');
		print $cgi->td({style=>"text-align: center;vertical-align:middle"},"<span><b>TOO MANY results...</b></span>");
		print $cgi->end_table();
		print "</div></div>";
		return;
	}
	if ($lGenes_ids[0] eq 'no_result' or $lGenes_ids[-1] eq 'no_result') {
		print $cgi->start_table({class=>"table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 10px;font-family:  Verdana;", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
		print $cgi->td({style=>"text-align: center;vertical-align:middle"},'');
		print $cgi->td({style=>"text-align: center;vertical-align:middle"},"<span><b>No result... Sorry..</b></span>");
		print $cgi->end_table();
		print "</div></div>";
		return;
	}
	my @lGenes = @{$patient->getProject->newGenes(\@lGenes_ids)};
	if (@lGenes) {
		my $html = table_cnv_genes_transcripts($patient, \@lGenes, $hGenes_dude);
		print $html;
	}
	else {
		print $cgi->start_table({class=>"table table-striped table-bordered table-hover",style=>"text-align: center;vertical-align:middle;font-size: 10px;font-family:  Verdana;", 'data-click-to-select'=>"true",'data-toggle'=>"table"});
		print $cgi->td({style=>"text-align: center;vertical-align:middle"},'');
		#print $cgi->td({style=>"text-align: center;vertical-align:middle"},"<span><b>NO RESULT <u>OR</u> ANALYSE NOT FOUND</b></span>");
		#print $cgi->td({style=>"text-align: center;vertical-align:middle"},'');
		print $cgi->end_table();
	}
	print "</div></div>";
	return;
}

sub print_all_cnv {
	my ($patient, $panel,$test) = @_;
	return;
	print print_cnv($patient, $panel,$patient->getBestOne());
}
sub print_validated_cnv {
	my ($patient, $validated) = @_;
	my $bestone = $patient->getBestOne();
	my $b;
	foreach my $event (@$bestone){
		
			my $start = $event->{START};
			$start =~ s/,//g;
			my $end = $event->{END};
			$end =~ s/,//g;
			my $chr_name = $event->{CHROM};
			$chr_name = "X" if $event->{CHROM} == 23;
			$chr_name = "Y" if $event->{CHROM} == 24;
			$chr_name = "MT" if $event->{CHROM} == 25;
			my $type =  lc($event->{TYPE});# eq "dup" ;
			my $id = $type."_".$chr_name."_".$start."_".$end;
			
			next unless exists $validated->{$id};
			$event->{validation} =  $validated->{$id}->[0];
			push(@$b,$event);
	}
	return print_cnv($patient, undef,$b,1);
}

	sub return_date {
		my ($dd) = @_;
		 my @amonths = ('Jan', 'Feb', 'Mar', 'Apr','May',"Jun","Jul","Aug","Sep","Oct","Nov","Dec");
		my ($date,$time) = split(" ",$dd);
	    my ($year,$month,$day) =  split("-",$date);
		return ($day,$amonths[$month-1],$year);
	}
	
	
sub print_cnv {
	my ($patient, $panel,$bestone,$only_validation) = @_;
	return;
}

sub getNbVarInRegion {
	my ($chr, $vector, $start, $end) = @_;
	my $nb = 0;
	return $nb unless $vector;
	foreach my $var_obj (@{$chr->getListVarObjects($vector)}) {
		$nb++ if ($var_obj->start() >= $start and $var_obj->start() <= $end);
	}
	return $nb;
}

sub print_polyCyto_button {
		my ($patient,$style) = @_;
		my $project = $patient->project;
		my $u1  = qq{https://www.polyweb.fr/polyweb//html/manta/PolyCyto.html};
		# $u1 = qq{https://www.polyweb.fr//polyweb/html/manta/Patient_CNV_Editor.html};
		my $url = "$u1?projectname=".$project->name."&filename=".$patient->name."&transloc=yes";
		unless($style){
		return qq{<a type="button" class= "btn btn-xs  btn-primary" href="$url" role="button" target="_blank"><img src="https://img.icons8.com/color/30/000000/biotech.png">&nbspPolyCyto</a>};
		}
		else {
				return qq{<a type="button" class= "btn btn-xs  btn-primary" href="$url" role="button" target="_blank" style="font-size:10px"><img src="https://img.icons8.com/color/24/000000/biotech.png">&nbspPolyCyto</a>};
			
		}
}

sub vspliceAI {
	my ($v,$hvariation) = @_;
	foreach my $g (@{$v->getGenes}){
		my $max_value;
		my $max_cat = '-';
		
		my $h_score_spliceAI = $v->spliceAI_score($g);
		#warn Dumper  $h_score_spliceAI;
		if ($h_score_spliceAI) {
			foreach my $cat (sort keys %$h_score_spliceAI) {
				if (not $max_value or $h_score_spliceAI->{$cat} > $max_value) {
					$max_value = $h_score_spliceAI->{$cat};
					$max_cat = $cat;
				}
			}
			$hvariation->{value}->{spliceAI}->{$g->id}= "$max_cat:$max_value";
			my $text = $max_cat.':'.$max_value;
			if ($max_value == 0) { $text = '0'; }
			unless (defined $max_value ) { $text = '-'; $max_value =0; }
			my $text_alert = "";
			$hvariation->{html}->{spliceAI}->{$g->id} = printButton($max_value,[0.5, 0.9],$text,undef,$text_alert);
			
		}
	}
} 
sub spliceAI {
	my ($project,$hvariation,$hgene) = @_;
	my $gid = $hgene->{id};
	#warn Dumper $hvariation->{value}->{spliceAI};
	return if (exists $hvariation->{value}->{spliceAI}->{$gid});
#	die();
	my @l_score_spliceAI;
	my $v_id = $hvariation->{value}->{var_name};
	return unless ($v_id);
	$v_id =~ s/-/_/g;
	my $v = $project->_newVariant($v_id);
	my $genes = $project->newGenes([$gid]);
	my $g =  $genes->[0];
	my $max_value;
	my $max_cat = '-';
	my $h_score_spliceAI = $v->spliceAI_score($g);
	if ($h_score_spliceAI) {
		foreach my $cat (sort keys %$h_score_spliceAI) {
			if (not $max_value or $h_score_spliceAI->{$cat} > $max_value) {
				$max_value = $h_score_spliceAI->{$cat};
				$max_cat = $cat;
			}
			if ($cat eq 'DG') { push(@l_score_spliceAI, 'Donor Gain: '.$h_score_spliceAI->{$cat}); }
			if ($cat eq 'DL') { push(@l_score_spliceAI, 'Donor Lost: '.$h_score_spliceAI->{$cat}); }
			if ($cat eq 'AG') { push(@l_score_spliceAI, 'Acceptor Gain: '.$h_score_spliceAI->{$cat}); }
			if ($cat eq 'AL') { push(@l_score_spliceAI, 'Acceptor Lost: '.$h_score_spliceAI->{$cat}); }
		}
	}
	foreach my $htr (@{$hvariation->{genes}->{$gid}}){
		value_html($htr,'spliceAI',join(",", @l_score_spliceAI),join(",", @l_score_spliceAI));
		if (defined($max_value)) {
			my $text = $max_cat.':'.$max_value;
			if ($max_value == 0) { $text = 'All:'.$max_value; }
			my $text_alert = 'SpliceAI values - '.join(', ', @l_score_spliceAI);
			$htr->{html}->{spliceAI} = printButton($max_value,[0.5, 0.9],$text,undef,$text_alert);
		}
		else {
			$htr->{html}->{spliceAI} = '-';
		}
	}
}

1;


