package polyviewer_html;
use strict;
use FindBin qw($Bin);
use Moose;
use Data::Dumper;
use JSON::XS;

has project => (
	is		=> 'ro',
	required=> 1,
);

has patient => (
	is		=> 'ro',
	required=> 1,
);
has patient_name => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		
		return $self->patient->name;
	},
);
has patient_id => (
	is		=> 'ro',
	lazy	=> 1,
	default => sub {
		my $self = shift;return $self->patient->id;
	},
);

has bam_files => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		$self->init_fam;
	},

);

has names => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		$self->init_fam;
	},
);

has cgi => (
	is		=> 'rw',
	default => sub {
		return  new CGI();
	},
);


has short_phenotype => (
	is		=> 'rw',
	lazy	=> 1,
	default => sub {
		my $self = shift;
		my $ph = "-";
		$ph =  $self->project->getPhenotypes->[0]->short_name() if $self->project->getPhenotypes && @{$self->project->getPhenotypes};
		return $ph;
	},
);
has variant => (
	is		=> 'rw',
);
	
my $minus = qq{<span class="glyphicon glyphicon-minus" aria-hidden="true"></span>};
my $plus  = qq{<span class="glyphicon glyphicon-plus" aria-hidden="true"></span>};
 my $boy = qq{<img src="https://img.icons8.com/color/24/000000/boy.png">};
 my $girl 	= 		qq{<img src="https://img.icons8.com/color/24/000000/girl.png">};
 my $female = qq{<i class="fa fa-venus fa-2x" aria-hidden="true" style="color:pink"></i>};
 my $male = qq{<i class="fa fa-mars fa-2x" aria-hidden="true" style="color:blue"></i>};
 my $mother_icon_ill =   qq{<center><img src='/icons/Polyicons/female-d.png'><center>};
 my $mother_icon_healthy =  qq{<center><img src='/icons/Polyicons/female-s.png'></center>};
 

sub color_genetic_model {
	my ($self,$model,$patient) = @_;
	return "" unless $model;
	my $color = "";
	#FFF7F8 mother line
	if ($model eq "+m"){
		$color = "#fff4f4";
	}
	elsif ($model eq "-m" && $patient->isMother){
		$color = ";opacity:0.5";
	}
	elsif ($model eq "+f"){
		$color = "#E2EEEE";
	}
	 elsif ($model eq "-f" && $patient->isFather){ 
		$color = ";opacity:0.5";
	}
	elsif (lc ($model) =~ /mother/i){
		$color = "#f2cdcd";
	}
	elsif ($model eq 'mother_c'){
		$color = "#779ECB";
	}
	elsif (lc($model)=~ /father/i){
		$color = "#A9BCD1";
	}
	elsif (lc($model) eq 'father_c'){
		$color = "#779ECB";
	}
	elsif ($model =~ /denovo/i){
		$color = "#e74c3c";
	}
	elsif ($model =~ /rece/i){
		$color = "#EE82EE";
	}
	elsif ($model =~ /mosa/i){
		$color = "#F9885C";
		$color = "#FDE6B0";
	}
	elsif ($model =~ /uni/i){
		$color = "#F9885C";
		$color = "#45B8AC";
	}
	return $color;
}

sub return_genetic_model {
	my ($self,$model) = @_;
	return "?" unless $model;
	if($model =~ /\+/){
		return  $plus;
	}
	elsif ($model  =~ /\-/) {
			return  $minus;
	}
	elsif (lc($model) =~ /mother/i) {
		return $self->patient->getFamily->getMother()->small_icon();
	}
	elsif (lc($model) =~ /father/i) { 
		return $self->patient->getFamily->getFather()->small_icon();
	}
	elsif  (lc($model) =~ /strict/i) {
			return  qq{Strict Denovo};
	}
	elsif  (lc($model) eq "denovo/?") {
			return  qq{Denovo/?};
	}
	elsif  (lc($model) =~ /denovo/i) {
			return  qq{Denovo};
	}
	elsif  (lc($model) =~ /recessive/i) {
			return  qq{Recessive};
	}
	elsif  (lc($model) eq "mosaic") {
			return  qq{mosaic};
	}	
	elsif  (lc($model) =~ /uniparental/i) {
			return  qq{Uniparental};
	}
	else {
		return $model;
	}

	
}

sub return_onclick {
	my ($self,$gene_name,$variation_id,$type) = @_;
	
	return qq{onClick="view_var_from_proj_gene_pat('}.join(qq{', '},$self->project->name,$gene_name,$self->patient_name,$variation_id,$type).qq{')"};
}

sub return_button_name {
	my ($self,$color,$sample,$onclick) = @_;
	#box-shadow: 1px 1px 2px $color;
	my $size = "9px";
	if (length($sample->name)>10){
		$size = "8px";
	}
	if (length($sample->name)>15){
		$size = "7px";
	}
	return qq{<button type="button" class="btn btn-primary btn-xs" style="font-size:$size;color:black;background-color:aliceblue;border-color:$color;" $onclick>}.$sample->name.qq{</button>};
}
sub html_parent_line_variant {
	my ($self,$hsample,$sample,$gene_name) = @_;
	my $color = "black";
	my $pid = $self->patient_id;
	#if ($hvariation->{patients}->{$pid}->{model} =~ /denovo/ && ($sample->isMother or $sample->isFather))  {
	#	$color = "lightgrey";
	#}
	#$hsample->{model} = $v->getTransmissionModelType($patient->getFamily(),$p,undef)  unless $hsample->{model};
#		die() unless  $hsample->{model};
	$hsample->{model} = "-" unless $hsample->{model};
	my $model_text = $self->return_genetic_model($hsample->{model});
	
	
	my $bgcolor = $self->color_genetic_model($hsample->{model},$sample);
	#[$name,$sample->small_icon,$hvariation->{patients}->{$sid}->{gt},$hvariation->{patients}->{$sid}->{pc}."%", $hvariation->{patients}->{$sid}->{dp},$model_text];
	my $pc_text = '-';
	$pc_text = $hsample->{pc}.'%' if ($hsample->{pc});
	my $tab = [$hsample->{button},$sample->small_icon,$hsample->{gt},$pc_text,$hsample->{dp},$model_text];
	my $style = "color:$color;";
	$style .= "background-color:$bgcolor" if $bgcolor; 
	my $html;
	$html =  $self->cgi->start_Tr({style=>"border: 1px solid black;color:black;$style"}).$self->cgi->td($tab).$self->cgi->end_Tr();
	return $html;
}

sub calling_variation {
	my ($self) = @_;
	my $gene = $self->variant->gene;
	my $gene_name =  $gene->{name};
	my $patient_name  = $self->patient->name;
	my $text_caller = join( "<br>", @{ $self->variant->text_caller } ); 
	
	my $samples = $self->variant->patients_calling();
	my $model = $samples->{$self->patient_id}->{infos_text};
	my $variation_id = $self->variant->id;
	$text_caller ="" unless $text_caller;
	$model ="denovo" unless $model;
	my $color = "#555";
	die() unless $model;
	$color = $self->color_genetic_model($model,$self->patient);
	my $id_info = 'ti_'.$variation_id.'_'.$gene->{id}.'_'.$self->patient_id;
	my $html = $self->cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",id=>$id_info,style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
	  my $fam = $self->patient->getFamily();
	  if ($self->patient->isChild) {
	  	if ($fam->getMother){
	  		my $hsample = $samples->{$fam->getMother->id};
	  		$hsample->{button} =  qq{<button style='color:black;' onClick="view_var_from_proj_gene_pat('}.join(qq{', '},$self->project->name,$gene_name,$self->patient_name,$variation_id,"mother").qq{')">}.$fam->getMother->name.qq{</button>};
	  		#$hsample->{button} = qq{<span class="badge badge-primary">}.$fam->getMother->name." ".$fam->getMother->small_icon.qq{</span>};
	  		#
			my $onclick = $self->return_onclick($gene_name,$variation_id,"mother");
	  		$hsample->{button} =  $self->return_button_name("pink",$fam->getMother,$onclick);
	  		
	  		$html.= "\n".$self->html_parent_line_variant($hsample,$fam->getMother,$gene_name,undef);
	  	}
	  	if ($fam->getFather){
	  		my $hsample = $samples->{$fam->getFather->id};
	  		$hsample->{button} =  qq{<button style='color:black;' onClick="view_var_from_proj_gene_pat('}.join(qq{', '},$self->project->name,$gene_name,$self->patient_name,$variation_id,"father").qq{')">}.$fam->getFather->name.qq{</button>};
	  		my $onclick = $self->return_onclick($gene_name,$variation_id,"father");
	  		$hsample->{button} =  $self->return_button_name("#1BA1E2",$fam->getFather,$onclick);
	  		$html.= "\n".$self->html_parent_line_variant($hsample,$fam->getFather,$gene_name,undef);
	  	}
	  	my $ch = 0;
	  	foreach my $child (@{$fam->getChildren()}){
	  			my $hsample = $samples->{$child->id};
	  			$ch++;
	  			next unless $hsample;
	  			
	  			$hsample->{button} =  qq{<button style='color:black;' onClick="view_var_from_proj_gene_pat('}.join(qq{', '},$self->project->name,$gene_name,$self->patient_name,$variation_id,"").qq{')">}.$child->name.qq{</button>};
	  			my $onclick = $self->return_onclick($gene_name,$variation_id,"");
	  			$hsample->{button} =  $self->return_button_name("#333333",$child,$onclick);
	  			$html.= "\n".$self->html_parent_line_variant($hsample,$child,$gene_name,$color,undef);
	  	}
	  }
	  else {
	  	my $hsample = $samples->{$self->patient->id};
	  	$html.= "\n".$self->html_parent_line_variant($hsample,$self->patient,$gene_name);
	  }
	  $html .= $self->cgi->end_table();
	  $html .= qq{<div data-dojo-type="dijit/Tooltip" data-dojo-props="connectId:'$id_info',position:['after']"><div style="min-width:300px;width:auto;height:auto;"><center><b><u>Patient $patient_name</b></u><br><br>};
	  $html .= qq{<b><u>Calling Method(s):</b></u><br> $text_caller </center></div></div>};
	return $html;	
}

sub calling {
	my ($self,$gene) = @_;
	if ($self->variant->isDude) {
		$self->calling_dude();
	}
	elsif ($self->variant->isCnv) {
		$self->calling_cnv();
	}
	else {
		$self->calling_variation();
	}
}

sub calling_dude {
	my($self) = @_;
	my $gene = $self->variant->gene;
	my $variation_id = $self->variant->id;
	my $gene_name =  $gene->{name};
	my $patient_name  = $self->patient->name;
	my $samples = $self->variant->patients_calling();
	my $text_caller = join( "<br>", @{ $self->variant->text_caller } ); 
	$text_caller ="" unless $text_caller;
	my $color = "#555";
		#$samples->{$self->patient_id}->{model} = "denovo";#$hsample->{model} = "-";
	$color = $self->color_genetic_model($samples->{$self->patient_id}->{model});
	my $html;
	my $id_info = 'ti_'.$variation_id.'_'.$gene->{id}.'_'.$self->patient_id;

		my $cnv_confidence = "";#$hvariation->{value}->{cnv_confidence}->{$patient_name};
		
		my $box_color = "grey";
		if ($cnv_confidence eq 'high') { $box_color = "red"; }
		elsif ($cnv_confidence eq 'medium') { $box_color = "orange"; }
		elsif ($cnv_confidence eq 'low') { $box_color = "blue"; }
		$html .= $self->cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",id=>$id_info,style=>"box-shadow: 1px 1px 6px $box_color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
		$html .=  $self->cgi->start_Tr({style=>"border: 1px solid black;color:black;text-align:center;font-weight: bold;"}).$self->cgi->td(['pat','status','he/ho','norm dp','cnv score','trans']).$self->cgi->end_Tr();	

	  my $fam = $self->patient->getFamily();
	  if ($self->patient->isChild) {
	  	if ($fam->getMother){
	  		my $hsample = $samples->{$fam->getMother->id};
	  	#	$hsample->{button} =  qq{<button style='color:black;' onClick="view_var_from_proj_gene_pat('}.join(qq{', '},$self->project->name,$gene_name,$self->patient_name,$variation_id,"mother").qq{')">}.$fam->getMother->name.qq{</button>};
	  		my $onclick = $self->return_onclick($gene_name,$variation_id,"mother");
	  		$hsample->{button} =  $self->return_button_name("pink",$fam->getMother,$onclick);
	  		$html.= "\n".$self->html_parent_line_variant_dude($hsample,$fam->getMother,$gene_name,undef);
	  	}
	  	if ($fam->getFather){
	  		my $hsample = $samples->{$fam->getFather->id};
	  			$hsample->{button} =  qq{<button style='color:black;' onClick="view_var_from_proj_gene_pat('}.join(qq{', '},$self->project->name,$gene_name,$self->patient_name,$variation_id,"father").qq{')">}.$fam->getFather->name.qq{</button>};
	  			my $onclick = $self->return_onclick($gene_name,$variation_id,"father");
	  			$hsample->{button} =  $self->return_button_name("#1BA1E2",$fam->getFather,$onclick);
	  			$html.= "\n".$self->html_parent_line_variant_dude($hsample,$fam->getFather,$gene_name,undef);
	  	}
	  	foreach my $child (@{$fam->getChildren()}){
	  		my $hsample = $samples->{$child->id};
	  		next unless $hsample;
	 # 		$hsample->{button} =  qq{<button type="button" class="btn btn-primary btn-xs" style="font-size:8px;color:black;font-weight:bold;background-color:aliceblue;border-color:black;box-shadow: 1px 1px 2px black" onClick="view_var_from_proj_gene_pat('}.join(qq{', '},$self->project->name,$gene_name,$self->patient_name,$variation_id,"father").qq{')">}.$child->name." ".$child->small_icon.qq{</button>};
	  		my $onclick = $self->return_onclick($gene_name,$variation_id,"");
	  		$hsample->{button} =  $self->return_button_name("#333333",$child,$onclick);
	  		$html.= "\n".$self->html_parent_line_variant_dude($hsample,$child,$gene_name,undef);
	  	}
	  }
	  else {
	  		$html.= "\n".$self->html_parent_line_variant_dude($self->patient,$gene_name);
	  }
	  $html .= $self->cgi->end_table();
	  $html .= qq{<div data-dojo-type="dijit/Tooltip" data-dojo-props="connectId:'$id_info',position:['after']"><div style="min-width:300px;width:auto;height:auto;"><center><b><u>Patient $patient_name</b></u><br><br>};
	  $html .= qq{<b><u>Calling Method(s):</b></u><br> $text_caller </center></div></div>};
	return $html;	
}

sub calling_cnv {
	my($self) = @_;
	my $gene = $self->variant->gene;
	my $variation_id = $self->variant->id;
	my $gene_name =  $gene->{name};
	my $patient_name  = $self->patient->name;
	my $samples = $self->variant->patients_calling();
	my $text_caller = join( "<br>", @{ $self->variant->text_caller } ); 
	$text_caller ="" unless $text_caller;
	my $color = "#555";
		#$samples->{$self->patient_id}->{model} = "denovo";#$hsample->{model} = "-";
	$color = $self->color_genetic_model($samples->{$self->patient_id}->{model});
	my $html;
	my $id_info = 'ti_'.$variation_id.'_'.$gene->{id}.'_'.$self->patient_id;

		my $cnv_confidence = "";#$hvariation->{value}->{cnv_confidence}->{$patient_name};
		
		my $box_color = "grey";
		if ($cnv_confidence eq 'high') { $box_color = "red"; }
		elsif ($cnv_confidence eq 'medium') { $box_color = "orange"; }
		elsif ($cnv_confidence eq 'low') { $box_color = "blue"; }
		$html .= $self->cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",id=>$id_info,style=>"box-shadow: 1px 1px 6px $box_color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
		my $pr_text = qq{<button onclick="alert('PR: Number of spanning read pairs which strongly (Q30) support the REF or ALT alleles')">PR</button>};
		my $sr_text = qq{<button onclick="alert('SR: Number of split-reads which strongly (Q30) support the REF or ALT alleles')">SR</button>};
		
		$html .=  $self->cgi->start_Tr({style=>"border: 1px solid black;color:black;text-align:center;font-weight: bold;"}).$self->cgi->td(['pat','status','he/ho',$pr_text,$sr_text,'norm dp','cnv score','trans']).$self->cgi->end_Tr();	

	  my $fam = $self->patient->getFamily();
	  if ($self->patient->isChild) {
	  	if ($fam->getMother){
	  		my $hsample = $samples->{$fam->getMother->id};
	  	#	$hsample->{button} =  qq{<button style='color:black;' onClick="view_var_from_proj_gene_pat('}.join(qq{', '},$self->project->name,$gene_name,$self->patient_name,$variation_id,"mother").qq{')">}.$fam->getMother->name.qq{</button>};
	  		my $onclick = $self->return_onclick($gene_name,$variation_id,"mother");
	  		$hsample->{button} =  $self->return_button_name("pink",$fam->getMother,$onclick);
	  		$html.= "\n".$self->html_parent_line_variant_cnv($hsample,$fam->getMother,$gene_name,undef);
	  	}
	  	if ($fam->getFather){
	  		my $hsample = $samples->{$fam->getFather->id};
	  			$hsample->{button} =  qq{<button style='color:black;' onClick="view_var_from_proj_gene_pat('}.join(qq{', '},$self->project->name,$gene_name,$self->patient_name,$variation_id,"father").qq{')">}.$fam->getFather->name.qq{</button>};
	  			my $onclick = $self->return_onclick($gene_name,$variation_id,"father");
	  			$hsample->{button} =  $self->return_button_name("#1BA1E2",$fam->getFather,$onclick);
	  			$html.= "\n".$self->html_parent_line_variant_cnv($hsample,$fam->getFather,$gene_name,undef);
	  	}
	  	foreach my $child (@{$fam->getChildren()}){
	  		my $hsample = $samples->{$child->id};
	  		next unless $hsample;
	 # 		$hsample->{button} =  qq{<button type="button" class="btn btn-primary btn-xs" style="font-size:8px;color:black;font-weight:bold;background-color:aliceblue;border-color:black;box-shadow: 1px 1px 2px black" onClick="view_var_from_proj_gene_pat('}.join(qq{', '},$self->project->name,$gene_name,$self->patient_name,$variation_id,"father").qq{')">}.$child->name." ".$child->small_icon.qq{</button>};
	  		my $onclick = $self->return_onclick($gene_name,$variation_id,"");
	  		$hsample->{button} =  $self->return_button_name("#333333",$child,$onclick);
	  		$html.= "\n".$self->html_parent_line_variant_cnv($hsample,$child,$gene_name,undef);
	  	}
	  }
	  else {
	  	$html.= "\n".$self->html_parent_line_variant_cnv($self->patient,$gene_name);
	  }
	  if ($self->variant->isMantaImprecise) {
	  	$html .=  $self->cgi->start_Tr({style=>"font-size:09px;position:sticky;z-index:8;top:0;background-color:#fd8253;color:white;"}).$self->cgi->th({style=>"text-align: center;",colspan=>8}, 'IS IMPRECISE')."".$self->cgi->end_Tr();
	  }
	  $html .= $self->cgi->end_table();
	  $html .= qq{<div data-dojo-type="dijit/Tooltip" data-dojo-props="connectId:'$id_info',position:['after']"><div style="min-width:300px;width:auto;height:auto;"><center><b><u>Patient $patient_name</b></u><br><br>};
	  $html .= qq{<b><u>Calling Method(s):</b></u><br> $text_caller </center></div></div>};
	return $html;	
}

sub html_parent_line_variant_dude {
	my ($self,$hsample,$sample,$gene_name) = @_;
	my $color = "black";
	my $pid = $self->patient_id;
	my $parent_type = "" ;
	my $sid = $sample->id;
	$parent_type = "mother" if $sample->isMother;
	$parent_type = "father" if $sample->isFather;
	$hsample->{model} = "-" unless $hsample->{model};
	my $model_text = $self->return_genetic_model($hsample->{model});
	my $bgcolor = $self->color_genetic_model($hsample->{model},$sample);
	$hsample->{dude} =0 unless exists  $hsample->{dude};
	my $cnv_score = $hsample->{dude_score};
	my $text_norm_dp_m = $hsample->{norm_depth}; 
	my $tab = [$hsample->{button},$sample->small_icon,$hsample->{gt},$text_norm_dp_m,$cnv_score,$model_text];
	my $style = "color:$color;";
	$style .= "background-color:$bgcolor" if $bgcolor; 
	my $html;
	$html =  $self->cgi->start_Tr({style=>"border: 1px solid black;color:black;$style"}).$self->cgi->td($tab).$self->cgi->end_Tr();
	return $html;
}

sub html_parent_line_variant_cnv {
	my ($self,$hsample,$sample,$gene_name) = @_;
	my $color = "black";
	my $pid = $self->patient_id;
	
	#warn Dumper $hvariation->{patients};
	#die();
	#if ($hvariation->{patients}->{$pid}->{model} =~ /denovo/ && ($sample->isMother or $sample->isFather))  {
	#	$color = "lightgrey";
	#}
	my $parent_type = "" ;
	my $sid = $sample->id;
	$parent_type = "mother" if $sample->isMother;
	$parent_type = "father" if $sample->isFather;
	#$hsample->{model} = "-";
	$hsample->{model} = "-" unless $hsample->{model};
	my $model_text = $self->return_genetic_model($hsample->{model});
	my $bgcolor = $self->color_genetic_model($hsample->{model},$sample);
	my $pr = $hsample->{pr} ;
	my $sr =  $hsample->{sr};
	
	$hsample->{dude} =0 unless exists  $hsample->{dude};
	my $cnv_score = $hsample->{dude_score};
	my $text_norm_dp_m = $hsample->{norm_depth}; 
	
	#my $name =  qq{<button style='color:$color;' onClick="view_var_from_proj_gene_pat('}.join(qq{', '},$self->project->name,$gene_name,$self->patient_name,$hvariation->{value}->{id},$parent_type).qq{')">}.$sample->name.qq{</button>};
	#[$name,$sample->small_icon,$hvariation->{patients}->{$sid}->{gt},$hvariation->{patients}->{$sid}->{pc}."%", $hvariation->{patients}->{$sid}->{dp},$model_text];
	#$hvariation->{patients}->{$sid}->{pc} = "-" unless $hvariation->{patients}->{$sid}->{pc};
	my $tab = [$hsample->{button},$sample->small_icon,$hsample->{gt},$pr, $sr,$text_norm_dp_m,$cnv_score,$model_text];
	my $style = "color:$color;";
	$style .= "background-color:$bgcolor" if $bgcolor; 
	my $html;
	$html =  $self->cgi->start_Tr({style=>"border: 1px solid black;color:black;$style"}).$self->cgi->td($tab).$self->cgi->end_Tr();
	return $html;
}

sub printSimpleBadge {
	my ($self,$value) = @_;
	my $color = "black";
	 return qq{<span class="badge badge-success badge-xs" style="border-color:black;background-color:#FFFFFF;color:$color;font-size :8px;">$value</span>} ;
}

sub init {
	my($self) = @_;
	$self->init_fam();
	$self->short_phenotype;
	$self->patient_name;
	$self->project->validations;
}

sub init_fam {
	my($self) = @_;
	my @bams;
	my @names;
	foreach my $p (@{$self->patient->getFamily->getPatients()}){
			next unless -e $p->getBamFileName;
			push(@bams,$p->bamUrl);
			push(@names,$p->name());
		 }
		 $self->names(\@names);
		  $self->bam_files(\@bams);
}



sub printInvButton {
	my ($self,$value,$types,$text,$othercss) = @_;
	 $value ="-" unless defined $value;
	my $btn_class  = qq{class= "btn btn-xs btn-primary " style="background-color: #D0D0D0;font-size: 7px;font-family:  Verdana;color:black"};
	if ($value eq "-" ){
			$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: #e74c3c;font-size: 7px;font-family:  Verdana;"} 
	}
	elsif ($value <=  $types->[0]){
		$btn_class = qq{class= "btn btn-xs btn-primary " style="background-color: #FF8800;font-size: 7px;font-family:  Verdana;;color:black"};
	}
	elsif( $value < $types->[1] ){
		$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: #e74c3c;font-size: 7px;font-family:  Verdana;color:white"} ;
	}
	return  qq{<button type="button" $btn_class $othercss>$text</button>};
}


sub abutton {
	my ($self,$href,$value) = @_;
	return qq{<a class="btn btn-xs btn-primary" href="$href" target="_blank" style="background-color: #D0D0D0;font-size: 7px;font-family:  Verdana;color:black" role="button">}.$value."</a>";
}

sub obutton {
	my ($self,$onclick,$value) = @_;
	
	#confess($value) unless defined $value;
	$value ="!!!!" unless defined $value;
	return qq{<a class="btn btn-xs btn-primary" onclick=\'$onclick\' target="_blank" style="background-color: #D0D0D0;font-size: 7px;font-family:  Verdana;color:black" role="button">}.$value."</a>";
}

sub put_text_minus {
	my ($self,$value) = @_;
	return $value  if $value ;
	return   qq{<span class="glyphicon glyphicon-minus" aria-hidden="true"></span>};
	
}


sub varsome {
my ($self,$debug) = @_;
my $url = qq{https://varsome.com/variant/hg19/}.$self->variant->gnomad_id;
my $text =qq{<a  type="button" class="btn btn-primary btn-xs" href="$url" target="_blank">V</a>};
return $text;
}

sub igv {
	my ($self,$debug) = @_;
	my $gn = $self->project->getVersion();
	my $project_name = $self->project->name;
	my $f =  join(";",@{$self->bam_files});#$patient->{obj}->bamUrl;;
	my $pnames = join(";",@{$self->names});
	my $locus = $self->variant->locus;
	my $v1 = "/";
	return qq{<button class='igvIcon2' onclick='launch_web_igv_js("$project_name","$pnames","$f","$locus","$v1","$gn")' style="color:black"></button>};
}




sub alamut {
	my ($self,$debug) = @_;
	
	my $start = $self->variant->start;
	my $a0 = $self->variant->ref_allele;
	my $a1 = $self->variant->allele;
	my $chr_name = $self->variant->chromosome;
	#war "test"
	my $qq5 = qq{	<button    class="alamutView3" onClick ="displayInAlamut('$chr_name',$start,['$a0','$a1']);"></button>};
	return $qq5;
#	die();
}

sub var_name {
	#my ($self,$hvariation,$debug) = @_;
	my ($self,$debug) = @_;
	my $project_name = $self->project->name();
	my $patient_name = $self->patient_name;
	my @genes = ();#split(";",$hvariation->{text}->{genes_name});
	my $type = $self->variant->type;
	my $name = $self->variant->name;
	my $gnomad_id = $self->variant->gnomad_id;
	my $is_gnomad = $self->variant->gnomad_id;
	my $html ="";
	
#		if (exists $hvariation->{value}->{cnv_details_genes} and scalar keys %{$hvariation->{value}->{cnv_details_genes}} > 1) {
#			my @lGenesNames;
#			push(@lGenesNames, "<table class='table table-striped table-condensed table-bordered table-hover table-mybordered'>");
#			push(@lGenesNames, "<tr>");
#			push(@lGenesNames, "<td style='text-align:center;'><b>Gene Name</b></td>");
#			push(@lGenesNames, "<td style='text-align:center;'><b>Nb Panel(s)</b></td>");
#			push(@lGenesNames, "<td style='text-align:center;'><b>Description</b></td>");
#			push(@lGenesNames, "</tr>");
#			foreach my $g_id (keys %{$hvariation->{value}->{cnv_details_genes}}) {
#				 push(@lGenesNames, "<tr>");
#				my $cmd = qq{"view_var_from_proj_gene_pat('$project_name', '$g_id', '$patient_name', '', '', '');"};
#			#	my $disabled = ''; 
#			#	$disabled = 'disabled' if ($g_id eq $g->{id});
#				my $this_g = $self->project->newGene($g_id);
#				my $nb_panels = scalar(@{$this_g->getPanels()});
#				my $panels_text = '';
#				foreach my $p (@{$this_g->getPanels()}) {
#						$panels_text .= $p->name()."\n";
#				}
#				my $cmd_alert_panel = qq{alert("$panels_text");};
#				push(@lGenesNames, "<td style='text-align:center;'><button onclick=$cmd >".$this_g->external_name()."</button></td>");
#				push(@lGenesNames, "<td style='text-align:center;'>".$nb_panels."</td>");
#				push(@lGenesNames, "<td style='text-align:center;'>".$this_g->description()."</td>");
#				push(@lGenesNames, "</tr>");
#			}
#			push(@lGenesNames, "</table>");
#			my $genes_text = join("", @lGenesNames);
#			my $id_info = 'b_multi_genes'.$vid.'_'.rand(50000).'_'.$self->patient->name();
#			$html .= qq{<br><br><button id="$id_info" type="button" class="btn btn-xs  btn-primary" style="background-color: #9796C4;font-size: 7px;font-family:  Verdana;color:white">Multi Genes !</button>};
#			$html .= qq{<div data-dojo-type="dijit/Tooltip" data-dojo-props="connectId:'$id_info',position:['after']"><div style="width:450px;height:auto;text-align:center;">$genes_text</div></div>};
#		}
	
	#$hvariation->{value}->{ac}  = -1 unless defined $hvariation->{value}->{ac} ;
	
	
	if ($type =~ m/large/) {
		#$vn = $hvariation->{value}->{name};
		my $text="";
		#if (exists $hvariation->{value}->{manta}->{is_imprecise} && defined $hvariation->{value}->{manta}->{is_imprecise}){
		 #	if ($hvariation->{value}->{manta}->{is_imprecise} == 1) {
			#			$text .=  qq{<br><span style='font-size:7px'><b><u>IS IMPRECISE</b></u></span>};
			#		}
		#}
					
		return $self->printSimpleBadge(qq{<a style="color:black">$name</a>}).$text.$html;
		
	}
	elsif ($type =~ m/del|ins|dup/) {
		return $self->printSimpleBadge(qq{<a href='https://gnomad.broadinstitute.org/variant/$gnomad_id' target = '_blank' style="color:black"><i class="fa fa-users fa-2x" style="color:coral"></i>&nbsp$name</a> });	;
	}
	elsif ($is_gnomad) {
		return $self->printSimpleBadge(qq{<a href='https://gnomad.broadinstitute.org/variant/$gnomad_id' target = '_blank' style="color:black"><i class="fa fa-users fa-2x" style="color:coral"></i>&nbsp$name</a> }).$html;
	}
	else {
		my @z = split("-",$gnomad_id);
		my $pp = $z[0]."-".$z[1];
		return  $self->printSimpleBadge(qq{<a href='https://gnomad.broadinstitute.org/region/$pp' target = '_blank' style="color:black">$name</a> }).$html;
	}
}

sub gnomad {
	my ($self) = @_;
	my $cgi=$self->cgi;
	my $color = "grey";
	my $ac = $self->variant->gnomad_ac;
	my $an = $self->variant->gnomad_an;
	unless ($an){
	my $html;
	 $html= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px red ;font-size: 7px;font-family:  Verdana;margin-bottom:0px;opacity:0.5"});
	$html.= $cgi->start_Tr();
	
	$html.= $cgi->th("AC");
	$html.= $cgi->th("Ho");
#$html.= $cgi->th('<i class="fa fa-mars" > </i> ') if  ($chr_name eq "X" or $chr_name eq "Y");
	$html.= $cgi->th("Max");
	$html.= $cgi->th("Min");
	$html.= $cgi->th("AN");
	$html.= $cgi->end_Tr();
	$html.= $cgi->start_Tr();
	$html.= $cgi->td(["-","-","-","-","-"]);
	$html.=$cgi->end_table();
	return $html;
	
	}
	my $min_pop = "-";
	 $min_pop = $self->variant->gnomad_min_pop_name."<br>".$self->variant->gnomad_min_pop if  $self->variant->gnomad_min_pop; # !!!!!! 	$text  = $v->min_pop_name."<br>".sprintf("%.4f", $min_pop) if $min_pop;
	
	my $max_pop = $self->variant->gnomad_max_pop_name."<br>".$self->variant->gnomad_max_pop;
	
	my $ac_ho = $self->variant->gnomad_ho;
	my $ac_ho_male = $self->variant->gnomad_ho_male;
	my $chr_name = $self->variant->chromosome;	
	my $start = $self->variant->start;	
	
 	$ac = 0 unless $ac;
	
	
	if ($ac <5){
		$color = "red";
	}
	elsif ($ac < 20){
		$color = "orange";
	} 
	elsif ($ac < 50){
		$color = "blue";
	} 
	
	my $pp = $chr_name."-".$start;
	
	my $href = qq{https://gnomad.broadinstitute.org/region/$pp};
#	my $vn= $hvariation->{value}->{gnomad_id};
	my $html= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color ;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
	$html.= $cgi->start_Tr();
	
	$html.= $cgi->th("AC");
	$html.= $cgi->th("Ho");
	$html.= $cgi->th('<i class="fa fa-mars" > </i> ') if  ($chr_name eq "X" or $chr_name eq "Y");
	$html.= $cgi->th("Max");
	$html.= $cgi->th("Min");
	$html.= $cgi->th("AN");
	$html.= $cgi->end_Tr();
	
	$html.= $cgi->start_Tr();
		
	
	$html.= $cgi->td($self->abutton($href,$self->put_text_minus($ac)));
	#$html.= $cgi->td($self-abutton($href,));
	$html.= $cgi->td($self->abutton($href,$self->put_text_minus($ac_ho)));
	$html.= $cgi->td($self->abutton($href,$self->put_text_minus($ac_ho_male))) if  ($chr_name eq "X" or $chr_name eq "Y");
	my $text =  qq{<span class="glyphicon glyphicon-minus" aria-hidden="true"></span>};
	$text = $max_pop;
	$html.= $cgi->td($self->abutton($href,$text));
	
	$text =qq{<span class="glyphicon glyphicon-minus" aria-hidden="true"></span>}; 
	$text = $min_pop;

	$html.= $cgi->td($self->abutton($href,$text));
	$html.= $cgi->td($self->abutton($href,$self->put_text_minus($an)));
	$html.= $cgi->end_Tr();
	$html.=$cgi->end_table();

	return $html;
}

sub dejavu {
		my ($self,$no_phenotype) = @_;
		my $opr = $self->variant->dejavu_other_projects;# = $v->other_projects();
		my $opa =0;
		 $opa =  $self->variant->dejavu_other_patients if defined  $self->variant->dejavu_other_patients ;;# = $v->other_patients();
		my $opho =  $self->variant->dejavu_other_patients_ho;# = $v->other_patients_ho();
		my $smpr =  $self->variant->dejavu_similar_projects;;#= $v->similar_projects();
		my $smpa = $self->variant->dejavu_similar_patients;
		my $smho = $self->variant->dejavu_similar_patients_ho;# = $v->similar_patients_ho();
		my $this = $self->variant->dejavu_this_run_patients;# = '-';
		my $vid =  $self->variant->id;
		#$hvariation->{value}->{this_run_patients} = $v->in_this_run_patients()."/".scalar(@{$project->getPatients});
		my $cgi = $self->cgi;
		
		my $color = "#555";
		if ($opa <5){
			$color = "red";
		}
		elsif ($opa < 10){
			$color = "orange";
		} 
		elsif ($opa < 50){
			$color = "blue";
		} 
		
#	my $deja_vu_url = "http://$server//polyweb/polydejavu/dejavu.html?input=";

	my $html = $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
	
	$html.= $cgi->start_Tr();
	$html.= $cgi->th("");
	$html.= $cgi->th("Pr");
	$html.= $cgi->th("Sa");
	$html.= $cgi->th("Ho");
	$html.= $cgi->end_Tr();
	$html.= $cgi->start_Tr();
	$html.= $cgi->td("other");
	
	#my $href = $deja_vu_url.$vid;
	my $onclick = qq{goto_dejavu("$vid")};
	$html.= $cgi->td($self->obutton($onclick,$opr));
	$html.= $cgi->td($self->obutton($onclick,$opa));
	$html.= $cgi->td($self->obutton($onclick,$opho));
	$html.= $cgi->end_Tr();
	
	if (not $no_phenotype) {
		$html.= $self->cgi->start_Tr();
		my $ps = $smpr;
		my $ph;
		$ph =  $self->short_phenotype;
		$html.= $cgi->td("$ph");
		$html.= $cgi->td($self->obutton($onclick,$smpr));
		$html.= $cgi->td($self->obutton($onclick,$smpa));
		$html.= $cgi->td($self->obutton($onclick,$smho));
		$html.= $cgi->end_Tr();
	}
	unless ($self->project->isGenome){
			$html .= $cgi->td("Run").$cgi->td({colspan=>3},$self->obutton($onclick,$this)).$cgi->end_Tr();
	}
	$html.=$cgi->end_table();
	return $html;	
}




sub printButtonWithAlert {
	my ($self,$value,$types,$text,$text_alert) = @_;

	my $btn_class  = qq{class= "btn btn-xs btn-primary " style="background-color: #D0D0D0;font-size: 7px;font-family:  Verdana;color:black"};
	$btn_class = qq{class= "btn btn-xs btn-primary " style="background-color: #FF8800;font-size: 7px;font-family:  Verdana;;color:white"} if $value >=  $types->[0] ;
	$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: #e74c3c;font-size: 7px;font-family:  Verdana;color:white"}  if $value >=  $types->[1] ;
	 $value ="-" unless defined $value;
	 
	$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: red;font-size: 7px;font-family:  Verdana;"}  if $value eq "-" ;
	$text  = $value unless $text;
	return  qq{<button onClick="alert('$text_alert')" "type="button" $btn_class>$text</button>} if ($text_alert);
	return  qq{<button type="button" $btn_class >$text</button>};
}
sub printButtonWithCmd {
	my ($self,$value,$types,$text,$cmd) = @_;

	my $btn_class  = qq{class= "btn btn-xs btn-primary " style="background-color: #D0D0D0;font-size: 7px;font-family:  Verdana;color:black"};
	$btn_class = qq{class= "btn btn-xs btn-primary " style="background-color: #FF8800;font-size: 7px;font-family:  Verdana;;color:white"} if $value >=  $types->[0] ;
	$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: #e74c3c;font-size: 7px;font-family:  Verdana;color:white"}  if $value >=  $types->[1] ;
	 $value ="-" unless defined $value;
	 
	$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: red;font-size: 7px;font-family:  Verdana;"}  if $value eq "-" ;
	$text  = $value unless $text;
	return  qq{<button  $cmd "type="button" $btn_class>$text</button>};# if ($text_alert);
}
sub printButton {
	my ($self,$value,$types,$text) = @_;

	my $btn_class  = qq{class= "btn btn-xs btn-primary " style="background-color: #D0D0D0;font-size: 7px;font-family:  Verdana;color:black"};
	$btn_class = qq{class= "btn btn-xs btn-primary " style="background-color: #FF8800;font-size: 7px;font-family:  Verdana;;color:white"} if $value >=  $types->[0] ;
	$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: #e74c3c;font-size: 7px;font-family:  Verdana;color:white"}  if $value >=  $types->[1] ;
	 $value ="-" unless defined $value;
	 
	$btn_class = qq{class= "btn btn-xs  btn-primary" style="background-color: red;font-size: 7px;font-family:  Verdana;"}  if $value eq "-" ;
	$text  = $value unless $text;
	return  qq{<button type="button" $btn_class >$text</button>};
}

sub validations {
		my ($self) = @_;
		my $cgi = $self->cgi;
		my $color = "#555";
		my $project = $self->project;
		
	#	confess() unless (exists $hvariation->{value}->{dm_for_this_gene});
	#	confess() unless (exists $hvariation->{value}->{clinvar_pathogenic_for_this_gene});

		if ($self->variant->dm or $self->variant->clinvar_pathogenic){
			$color = "red";
		}
		elsif  ($self->variant->hgmd_id or $self->variant->clinvar_id) {
			$color = "orange";
		}
		
		my $html= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"max-width:300px;box-shadow: 1px 1px 6px $color;font-size: 7px;font-family:  Verdana;margin-bottom:0px"});
		$html .= $cgi->start_Tr().$cgi->th(["HGMD","Clinvar","Local"]).$cgi->end_Tr();
		$html .= $cgi->start_Tr();
		my $all_validations = $project->validations;
		my $val_id = $self->variant->gene->{id}."!".$self->variant->id;
		my $local_validation = $project->getValidationVariation($val_id,$self->patient);
		
		my $polyweb = $minus ; 
	
		if ($local_validation){
				my $saved = $local_validation->{validation};
				#$hvariatio printButton(4,[3,4],$v->hgmd->{class},qq{onClick="$cmd"}); 
				my $id = $self->variant->id;
				my $gene_id = $self->variant->gene->{id};
				my $cmd = qq{view_variation_validation(\'$id\', \'$gene_id\')};
				$polyweb = $self->printButtonWithAlert($saved,[4,5],$project->buffer->value_validation->{$saved},qq{onClick="$cmd"}) ;
		}
		if ($self->variant->hgmd_id){
	 	my $n1 = $self->project->name;
	 	my $n2 = $self->variant->hgmd_id;
	 	my $n3 = $self->variant->id;
	 	my $cmd = "";
	 	#if (exists $hvariation->{html}->{no_css_polydiag}) { $cmd = qq{zoomHgmdWithoutCss(\'$n1\',\'$n2\',\'$n3\')}; }
	 	#else { 
	 	$cmd = qq{zoomHgmd(\'$n1\',\'$n2\',\'$n3\')}; 
	 	
	 	my $nb =  $self->variant->hgmd_value;

		$html   .= $cgi->td($self->printButtonWithCmd($nb,[3,4],$self->variant->hgmd,qq{onClick="$cmd"}) ); 
			
		}
		else {
			$html   .= $cgi->td($minus); 
		}
		if ($self->variant->clinvar){
			my $nb =  $self->variant->clinvar_value;
			my $cmd ="";
			my $uc = qq{https://www.ncbi.nlm.nih.gov/clinvar/?term=}.$self->variant->clinvar_id."[alleleid]";
	 		my $oc = qq{onClick='window.open("$uc")'};
			$html   .= $cgi->td($self->printButtonWithCmd($nb,[4,5],$self->variant->clinvar,$oc ) ); 
		}
		else {
			$html   .= $cgi->td($minus); 
		}
		
		$html .= $cgi->td($polyweb);
		
		my $hgmd_no_access = qq{<span class="glyphicon glyphicon-ban-circle" aria-hidden="true" style='font-size:12px;color:black;'></span>};
		#$hvariation->{html}->{table_validation_hgmd_no_access} = $html;
		#$hvariation->{html}->{table_validation_hgmd_no_access} .= $cgi->td([$hgmd_no_access,$hvariation->{html}->{clinvar},$polyweb]);
		
		my $v_phen = ' ';
   		if ( $self->variant->hgmd_phenotype) {
   			$v_phen = $self->variant->hgmd_phenotype;
   			$v_phen =~ s/"//g;
			$html .= qq{<tr><td  colspan='3' style='color:#4CAF50;'>$v_phen</td></tr>};
			#$hvariation->{html}->{table_validation_hgmd_no_access} .= qq{<tr><td  colspan='3' style='color:#4CAF50;'>$v_phen</td></tr>};
   		}
		$html .= $cgi->end_Tr().$cgi->end_table;
		#$hvariation->{html}->{table_validation_hgmd_no_access} .= $cgi->end_Tr().$cgi->end_table;
		return $html;
		
}
my @header_transcripts = ("consequence","enst","nm","ccds","appris","exon","nomenclature","codons","codons_AA", "polyphen","sift","cadd","revel","dbscsnv",'spliceAI');
my @header_transcripts_cnv = ( "consequence", "enst", "nm", "ccds", "appris", "start", "end" );
sub href {
	my ($url,$name) = @_;
	$url.=$name;
	return qq{<a href="$url" target="_blank">}.$name."</a>";
}

sub printInvBadge {
        my ($self,$value,$types) = @_;
         return $minus unless $value;
		return $minus if $value eq "-";
	
        my $color = "#4CAF50";
         #$color = "#D62C1A" if $value > $type[0] ;
         $color = "#FF8800" if $value < $types->[0] ;
         $color = "#FF0025" if $value < $types->[1] ;

         $value ="-" unless defined $value;
         $color = "#FF0025" if $value eq "-" ;

         return qq{<span class="badge badge-success badge-xs" style="border-color:$color;background-color:#FFFFFF;color:$color;font-size :8px;">$value</span>} ;
}

sub transcripts_cnv {
	my ($self) =@_;

	my $color_header = "white";
	if ($self->variant->type eq 'large_deletion') {
		$color_header = "#fd8253";
	}
	elsif ($self->variant->type eq 'large_duplication') {
		$color_header = "#5393fd";
	}
	my $cgi = new CGI();
	my $html;
	$html .= qq{<div style="max-height:250px;overflow-y:auto;border:solid 1px grey;">};
	$html .= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color_header;font-size: 7px;font-family: Verdana;margin-bottom:0px;max-height:150px;"});
	$html .=  $cgi->start_Tr({style=>"position:sticky;z-index:8;top:0;background-color:$color_header;color:white;"}).$cgi->th({style=>"text-align: center;"}, \@header_transcripts_cnv)."".$cgi->end_Tr();
	
	my $atr = $self->variant->transcripts;
	my $gene = $self->variant->gene;
	my $gene_id = $gene->{id};
	foreach my $htr (@$atr){
		my $t_id = $htr->{trid};
		next if not $self->variant->cnv_details_genes;
		next if not exists $self->variant->cnv_details_genes->{$gene_id};
		next if not exists $self->variant->cnv_details_genes->{$gene_id}->{$t_id};
		my @l_pos_e_i = sort {$a <=> $b} keys %{$self->variant->cnv_details_genes->{$gene_id}->{$t_id}->{positions}};
		next unless @l_pos_e_i;
		my $first = $self->variant->cnv_details_genes->{$gene_id}->{$t_id}->{positions}->{$l_pos_e_i[0]};
		my $last;
		if (scalar(@l_pos_e_i) > 1) {
			$last = $self->variant->cnv_details_genes->{$gene_id}->{$t_id}->{positions}->{$l_pos_e_i[-1]};
		}
		$html .=  $cgi->start_Tr();
		$html.= "<td>".$htr->{consequence}."</td>";
		$html.= "<td>".$htr->{enst}."</td>";
		if ($htr->{nm}) { $html.= "<td>".$htr->{nm}."</td>"; }
		else { $html.= "<td>-</td>"; }
		if ($htr->{ccds}) { $html.= "<td>".$htr->{ccds}."</td>"; }
		else { $html.= "<td>-</td>"; }
		$html.= "<td>".$htr->{appris}."</td>";
		if ($last) {
			$html.= "<td><b>".$first."</b></td>";
			$html.= "<td><b>".$last."</b></td>";
		}
		else {
			$html.= "<td colspan='2'><b>".$first."</b></td>";
		}
		$html .=  $cgi->end_Tr();
	}
	
	if (@{$self->variant->other_genes}){
	my $text = "Other Genes :".join("-",@{$self->variant->other_genes});
	
	$html .=  $cgi->start_Tr({style=>"font-size:09px;position:sticky;z-index:8;top:0;background-color:#363945;color:white;"}).$cgi->th({style=>"text-align: center;",colspan=>7}, $text)."".$cgi->end_Tr();
	}
	$html.=$cgi->end_table();
	return $html;
}

sub transcripts {
		my ($self) = @_;
		if ($self->variant->isCnv){
			return $self->transcripts_cnv();
		}
		return $self->transcripts_variants();
}

sub transcripts_variants {
	my ($self) =@_;
	my $gene = $self->variant->gene;
	
	my $atr = $self->variant->transcripts;
	
	 my $revel = $self->variant->revel;
	
	  unless( $revel eq "-"){
	 	$revel = $self->printBadge($revel,[0.5,0.9]);
	 }
	 else {
	 	$revel = $minus;
	 }
	 my $dbscsnv = $self->variant->rf;
	 	my $ada =  $self->variant->ada;
	  unless( $dbscsnv eq "-"){
	 	$dbscsnv = $self->printBadge($dbscsnv,[0.6,0.9]).$self->printBadge($ada,[0.6,0.9]);
	 }
	 else {
	 	$dbscsnv = $minus;
	 }

	 my $cadd =  $self->variant->cadd;
	 $cadd = "-" unless $cadd;
	 unless( $cadd eq "-"){
	 	$cadd = $self->printBadge($cadd,[20,30]);
	 }
	 else {
	 	$cadd = $minus;
	 }
	 
	 my $spliceAI = $self->variant->spliceAI;
	 
	 my $cat = $self->variant->spliceAI_cat;
	 unless ($spliceAI){
	 	$spliceAI = -1;
	 }
	 unless ($spliceAI == -1){
			my $text = $cat.':'.$spliceAI;
			if ($spliceAI == 0) { $text = '0'; }
			my $text_alert = "";
			
			$spliceAI = $self->printButtonWithAlert($spliceAI,[0.5, 0.9],$text,$text_alert);
	 }
	 else {
	 	$spliceAI = $minus;
	 }
	 
	my $level = 2;
	my $cgi          = $self->cgi;
	#my $value = $atr->[0]->{value}->{impact_score};
	my $color = "#555";
	#if ($value > 3){
	#	$color = "red";
	#}
	#if ($value >= 2){
	#	$color = "#FF8800";
	#}
	
	my $html;
	$html .= qq{<div style="max-height:250px;overflow-y:auto;border:solid 1px grey;">};
	
	$html .= $cgi->start_table({class=>"table table-sm table-striped table-condensed table-bordered table-primary ",style=>"box-shadow: 1px 1px 6px $color;font-size: 7px;font-family: Verdana;margin-bottom:0px;max-height:150px;"});
	my @colors = ("#F9F6FF","#F9F9F9");
	
	my $nb =0;
	my $nb_skip = 0;
	 
	$html .=  $cgi->start_Tr({style=>"position:sticky;z-index:8;top:0;background-color:grey;color:white;"}).$cgi->th({style=>"text-align: center;"}, \@header_transcripts)."".$cgi->end_Tr();
	 
	my $rids = [];
	my $hhide;
	my $compact_table;
	foreach my $htr (@$atr){
			my $zid = $htr->{enst}."-".$self->variant->gnomad_id;
			my $hide ="";
			
			#{style=>"border: 1px solid black;color:lightgrey;$style"}).$cgi->td([$mother->name,qq{<i class="fa fa-female fa-2x" aria-hidden="true" style="color:ligthgrey"></i>},$type,$ps."%",""]);
		 	#$htr->{html}->{spliceAI} = "-" ;
		 	my $hide_id= "-";
		 	 $hide_id = $htr->{ccds} if exists $htr->{ccds} && $htr->{ccds};
		 #	my $hide_id = 
		 	$hide_id .= $htr->{nomenclature} if exists $htr->{nomenclature} && $htr->{nomenclature};# eq "" ;
		 	
		 	$hide = "display:none;"  if ($htr->{impact_score}  <  $level and $nb >0) or (exists $hhide->{$hide_id})  ;
		 	$hide = "display:none;"  if ($compact_table and not $htr->{value}->{ccds} and $nb != $nb_skip)  ;
		 	$hide = "display:none;"  if ($compact_table and $nb - $nb_skip >= 3)  ;
		 	$hhide->{$hide_id} ++ ;
		 	$hhide->{$htr->{nomenclature}} ++ if exists $htr->{nomenclature};
			$nb_skip ++ if $hide;
			my $c = $colors[$nb%2];
			$nb ++;
			my $rid = "rowTR_"."_".$zid;
			push(@$rids,$rid) if $hide;
			
		 $html .=  $cgi->start_Tr({id=>$rid,style=>"border: 1px solid;background-color:$c ;".$hide});
		 ##cosnequence 
		 $html.= $cgi->td( $self->printButton($htr->{impact_score},[3,4],$htr->{consequence}) );
		 $html.= $cgi->td($self->cgi->a({href=>qq{https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=}.$htr->{enst},target=>"_blank"},$htr->{enst}));
		 $html.= $cgi->td($htr->{nm});
		 $htr->{ccds} = "-" unless $htr->{ccds};
 		 $html.= $cgi->td($self->cgi->a({href=>qq{https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=}.$htr->{ccds},target=>"_blank"},$htr->{ccds}));
		 $html.=  $cgi->td($self->printBadge($htr->{appris}));
		 $html.=  $cgi->td($self->printBadge($htr->{exon}));
		 $html.=  $cgi->td($self->printBadge($htr->{nomenclature}));
		 $html.=  $cgi->td($self->printBadge($htr->{codons}));
		 $html.=  $cgi->td($self->printBadge($htr->{codons_AA})); 
		 
		 $html.= $cgi->td($self->printBadge($htr->{polyphen},[0.446,0.908]));
		 $html.=  $cgi->td($self->printInvBadge($htr->{sift},[0,1,0,05]));
		 $html.= $cgi->td($cadd);
		 $html.= $cgi->td($revel);
		 $html.= $cgi->td($dbscsnv);
		 $html.= $cgi->td($spliceAI);
		 
		 
		 $html.= $cgi->end_Tr();
	}
	 if ($nb_skip ==0){
		$html.=$cgi->end_table();
		return $html
	 }
	my $js = encode_json $rids;
	my $za = "hide_tr_".time."_".int(rand(50000));
	$html .=  $cgi->start_Tr({id=>$za});

	$html.= $cgi->td({style=>"box-shadow: 1px 1px 2px #555;background-color:#CECFCE;",colspan=>scalar(@header_transcripts),onClick=>qq{showTranscripts($js,"$za");}},qq{<span class="glyphicon glyphicon-plus"></span> }."view $nb_skip Transcripts");
	$html.= $cgi->end_Tr();
	
	$html.=$cgi->end_table();
	$html.=qq{</div>};
	#$atr->{html} = $html;
	return $html;
	
}

sub printBadge {
	my ($self,$value,$types) = @_;
	return $minus unless defined $value;
	return $minus if  $value eq "-";
	my $color = "#4CAF50";
	if ($types){
	 #$color = "#D62C1A" if $value > $type[0] ;
	 $color = "#FF8800" if $value > $types->[0] ;
	 $color = "#FF0025" if $value > $types->[1] ;
	}
	 return qq{<span class="badge badge-success badge-xs" style="border-color:$color;background-color:#FFFFFF;color:$color;font-size :8px;">$value</span>} ;

}


sub validation_select{
	my ($self) = @_;
	
	my $cgi = $self->cgi;
	my $patient = $self->patient;

	$cgi =  new CGI unless $cgi;
	my $buffer = $patient->buffer;
	my $project = $patient->project;
	my $out;
	my $pname = $patient->name;
	my $vid = $self->variant->id;
	my $gene = $self->variant->gene;
	my $tt = $patient->name."_".$self->variant->gene->{js_id};
	my $menu = $tt."_menu";
	my $sp = $tt."_span";
	my $tdid = $gene->{id}."_".rand(50);
			

	my $bgcolor = "info";

my $val_id = $gene->{id}."!".$self->variant->id;

my $saved =[] ;
my $all_validations = $patient->validations;
my $validation_term ="";
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
	return $out;
}




1;
