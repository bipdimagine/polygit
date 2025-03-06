/**
 * @author pnitschk and mbras 
 */


var formatters = new Array();
formatters["formatMethods"] = formatMethods;
formatters["formatReferenceGene"] = formatReferenceGene;
formatters["formatSnp"] = formatSnp;
formatters["formatCnv"] = formatCnv;
formatters["formatFilter"] = formatFilter;
formatters["formatValidWithFilter"] = formatValidWithFilter;
formatters["formatValid"] = formatValid;
formatters["formatPolyphen1"] = formatPolyphen1;
formatters["formatAlign"] = formatAlign;
formatters["formatHapmap"] = formatHapmap;
formatters["formatHomozygote"] = formatHomozygote;
formatters["formatPolyphen"] = formatPolyphen;
formatters["formatSift"] = formatSift;
formatters["formatCachePolyPhen"] = formatCachePolyPhen;
formatters["formatPolyphenUrl"] = formatPolyphenUrl;
formatters["formatSimple"] = formatSimple;
formatters["formatProteinName"] = formatProteinName;
formatters["formatTranscriptName"] = formatTranscriptName;
formatters["formatDgv"] = formatDgv;
formatters["formatChild"] = formatChild2;
formatters["formatDisease"] = formatDisease;
formatters["formatConsequence"] = formatConsequence;
formatters["formatVarId"] = formatVarId;
formatters["formatCellGeneTrans"] = formatCellGeneTrans;
formatters["formatCellDejaVU"] = formatCellDejaVU;
formatters["formatCellPatDetails"] = formatCellPatDetails;
formatters["formatDejavu"] = formatDejavu;
formatters["formatRsName"] = formatRsName;
formatters["formatEmail"] = formatEmail;
formatters["statusIcon"] = statusIcon;
formatters["dejavu_capture_diag"] = dejavu_capture_diag;
formatters["dejavu_genes"] = dejavu_genes;
formatters["formatRegionHo"] = formatRegionHo;
formatters["formatRegionHo2"] = formatRegionHo2;
formatters["formatColorCell"] = formatColorCell;
formatters["formatHoRegionUcsc"] = formatHoRegionUcsc;
formatters["formatHoRegionEnsembl"] = formatHoRegionEnsembl;
formatters["formatHgmd"] = formatHgmd;
formatters["has_hgmd"] = has_hgmd;
formatters["write_bold_purple"] = write_bold_purple;
formatters["write_bold_lightpurple"] = write_bold_lightpurple;
formatters["write_bold_red"] = write_bold_red;
formatters["write_bold_lightred"] = write_bold_lightred;
formatters["write_bold_orange"] = write_bold_orange;
formatters["write_bold_lightorange"] = write_bold_lightorange;
formatters["write_bold_green"] = write_bold_green;
formatters["write_bold_lightgreen"] = write_bold_lightgreen;
formatters["write_bold_grey"] = write_bold_grey;
formatters["write_bold_lightgrey"] = write_bold_lightgrey;
formatters["format_is_omim"] = format_is_omim;
formatters["format_is_omim_morbid"] = format_is_omim_morbid;
formatters["format_phenotype"] = format_phenotype;
formatters["format_displayButtonWithCmd"] = format_displayButtonWithCmd;
formatters["formatNcboost"] = formatNcboost;
formatters["formatGeneMaxScoreVariant"] = formatGeneMaxScoreVariant;
formatters["formatPolywebScoreDetailProject"] = formatPolywebScoreDetailProject;
formatters["formatSpliceAI"] = formatSpliceAI;
formatters["formatColorPhenotypes"] = formatColorPhenotypes;
formatters["formatColorProjectAndPhenotypes"] = formatColorProjectAndPhenotypes;
formatters["formatColorCadd"] = formatColorCadd;
formatters["formatGnomadAC"] = formatGnomadAC;
formatters["formatGnomadLink"] = formatGnomadLink;




function formatGnomadLink(value) {
	if (value == '-') { return '-'; }
	var list = value.split(';');
	var res = list[0];
	var link = list[1];
	return "<a href='"+link+"' target='_blank'>"+res+"</a>";
}

function formatSpliceAI(value) {
	if (value == '-') { return value; }
	var max_value = 0;
	var max_cat = 'All';
	var list_tmp = value.split(' ');
	var list_values = [];
	for (i=0;i<list_tmp.length;i++) {
		var list_tmp_2 = list_tmp[i].split('=');
		var cat = list_tmp_2[0];
		var cat_value = list_tmp_2[1];
		if (cat_value > max_value) {
			max_value = cat_value;
			max_cat = cat;
		}
		if (cat == 'DG') { list_values.push('Donor Gain: ' + cat_value); }
		if (cat == 'DL') { list_values.push('Donor Lost: ' + cat_value); }
		if (cat == 'AG') { list_values.push('Acceptor Gain: ' + cat_value); }
		if (cat == 'AL') { list_values.push('Acceptor Lost: ' + cat_value); }
	}
	var text = list_values.join(", ");
	var btn_class = "class='btn btn-xs btn-primary' style='border-radius:10px;background-color:#D0D0D0;font-size:7px;font-family:Verdana;color:black'";
	if (max_value >= 0.5 ) { btn_class = "class='btn btn-xs btn-primary' style='border-radius:10px;background-color:#FF8800;font-size:7px;font-family:Verdana;color:white'"; }
	if (max_value >= 0.8 ) { btn_class = "class='btn btn-xs btn-primary' style='border-radius:10px;background-color:#e74c3c;font-size:7px;font-family:Verdana;color:white'"; }
	var btn = "<button onClick='alert(\"" + text + "\")' type='button' " + btn_class + ">" + max_cat + " " + max_value + "</button>";
	return btn;
}

function formatPolywebScoreDetailProject(value) {
	var res = '';
	var lValues = value.split(";");
	for (i=0;i<lValues.length;i++) {
		var this_value = lValues[i];
		var new_value = '';
		if (this_value == '') {}
		else if (this_value >= 10) { new_value = "<img  width=16 height=16 src='/icons/Polyicons/difficulty_2_1.png'>"; }
		else if (this_value >= 7)  { new_value =  "<img  width=16 height=16 src='/icons/Polyicons/difficulty_1.png'>"; }
		else if (this_value >= 5)  { new_value = "<img  width=16 height=16 src='/icons/Polyicons/difficulty_0.png'>"; }
		else { new_value = '.'; }
		if (i > 0) { res += '<br>'; }
		res += new_value;
	}
	return res;
}

function formatGeneMaxScoreVariant(value) {
	if (value) {
		if (value == -999) { return ' '; }
		if (value >= 10) { return "<img  width=16 height=16 src='/icons/Polyicons/difficulty_2_1.png'> " + value; }
		if (value >= 7) { return "<img  width=16 height=16 src='/icons/Polyicons/difficulty_1.png'> " + value; }
		if (value >= 5) { return "<img  width=16 height=16 src='/icons/Polyicons/difficulty_0.png'> " + value; }
//		if (value >= 10) { return "<img  width=16 height=16 src='/icons/Polyicons/difficulty_2_1.png'> "+value; }
//		if (value >= 7) { return "<img  width=16 height=16 src='/icons/Polyicons/difficulty_1.png'> "+value; }
//		if (value >= 5) { return "<img  width=16 height=16 src='/icons/Polyicons/difficulty_0.png'> "+value; }
	}
	return ' ';
}

function formatNcboost(value) {
	if (value) {
		var lValues = value.split(";");
		if      (String(lValues[0]).toLowerCase() == 'high')     { return "<span style='color:#EA3434;'>High " + lValues[1] + "</span>"; }
		else if (String(lValues[0]).toLowerCase() == 'medium')   { return "<span style='color:#EAA134;'>Medium " + lValues[1] + "</span>"; }
		else if (String(lValues[0]).toLowerCase() == 'moderate') { return "<span style='color:#EACC34;'>Moderate " + lValues[1] + "</span>"; }
		else if (String(lValues[0]).toLowerCase() == 'all') 	 { return "<span style='color:black;'>" + lValues[1] + "</span>"; }
	}
	return '-';
}

function format_phenotype(value) {
	if (value) {
		value = value.replace("Hpo: ", "<span style='color:red;'>Hpo: </span>");
		value = value.replace("Hgmd: ", "<span style='color:red;'>Hgmd: </span>");
		value = value.replace("Ens: ", "<span style='color:red;'>Ens: </span>");
		value = value.replace("Ddg2p: ", "<span style='color:red;'>Ddg2p: </span>");
		value = value.replace("Omim: ", "<span style='color:red;'>Omim: </span>");
		return "<span style='color: #29b6f6;'>" + value + "</span>";
	}
	return ' ';
}

function format_is_omim(value) {
	if (value) {
		var res = '';
		var lValues = value.split(";");
		for (i=0;i<lValues.length;i++) {
			if (i > 0) { res += "<br>"; }
			if (lValues[i] == 'new') {
				res += " <span style='color:red;'><i>New!</i></span>";
			}
			else if (lValues[i]) {
				var url = "https://www.omim.org/entry/" + lValues[i] + "?search=" + lValues[i] + "&highlight=" + lValues[i];
				res += "<a href='"+url+"' target='_blank'>Omim</a>";
			}
			else { res += ' '; }
		}
		return res;
	}
	return ' ';
}

function format_is_omim_morbid(value) {
	//if (value == 1) { return "<img src='/icons/Polyicons/12-em-check.png'>"; }
	if (value == 1) { return "<span style='color:red;'><i>morbid</i></span>"; }
	if (value == '1;new') { return "<span style='color:red;'><i>morbid</i></span> <span style='color:red;'><i>New!</i></span>"; }
	return ' ';
}

function write_bold_lightpurple(value) {
	if (value == 0) { return "<button class='no-click' style='background-color:white;border-radius:10px;border:none;'><span style='color:grey;'>0</span></button>"; }
	return "<button class='no-click' style='background-color:#cca7d6;border-radius:10px;border:none;'><span style='color:white;'><b>" + value + "</b></span></button>";
}

function write_bold_purple(value) {
	if (value == 0) { return "<button class='no-click' style='background-color:white;border-radius:10px;border:none;'><span style='color:grey;'>0</span></button>"; }
	return "<button class='no-click' style='background-color:mediumorchid;border-radius:10px;border:none;'><span style='color:white;'><b>" + value + "</b></span></button>";
}

function write_bold_lightred(value) {
	if (value == 0) { return "<button class='no-click' style='background-color:white;border-radius:10px;border:none;'><span style='color:grey;'>0</span></button>"; }
	return "<button class='no-click' style='background-color:#fe8181;border-radius:10px;border:none;'><span style='color:white;'><b>" + value + "</b></span></button>";
}

function write_bold_red(value) {
	if (value == 0) { return "<button class='no-click' style='background-color:white;border-radius:10px;border:none;'><span style='color:grey;'>0</span></button>"; }
	return "<button class='no-click' style='background-color:red;border-radius:10px;border:none;'><span style='color:white;'><b>" + value + "</b></span></button>";
}

function write_bold_lightorange(value) {
	if (value == 0) { return "<button class='no-click' style='background-color:white;border-radius:10px;border:none;'><span style='color:grey;'>0</span></button>"; }
	return "<button class='no-click' style='background-color:#fed68d;border-radius:10px;border:none;'><span style='color:white;'><b>" + value + "</b></span></button>";
}

function write_bold_orange(value) {
	if (value == 0) { return "<button class='no-click' style='background-color:white;border-radius:10px;border:none;'><span style='color:grey;'>0</span></button>"; }
	return "<button class='no-click' style='background-color:orange;border-radius:10px;border:none;'><span style='color:white;'><b>" + value + "</b></span></button>";
}

function write_bold_lightgreen(value) {
	if (value == 0) { return "<button class='no-click' style='background-color:white;border-radius:10px;border:none;'><span style='color:grey;'>0</span></button>"; }
	return "<button class='no-click' style='background-color:#72bf72;border-radius:10px;border:none;'><span style='color:white;'><b>" + value + "</b></span></button>";
}

function write_bold_green(value) {
	if (value == 0) { return "<button class='no-click' style='background-color:white;border-radius:10px;border:none;'><span style='color:grey;'>0</span></button>"; }
	return "<button class='no-click' style='background-color:green;border-radius:10px;border:none;'><span style='color:white;'><b>" + value + "</b></span></button>";
}

function write_bold_grey(value) {
	if (value == 0) { return "<button class='no-click' style='background-color:white;border-radius:10px;border:none;'><span style='color:grey;'>0</span></button>"; }
	return "<button class='no-click' style='background-color:grey;border-radius:10px;border:none;'><span style='color:white;'><b>" + value + "</b></span></button>";
}

function write_bold_lightgrey(value) {
	if (value == 0) { return "<button class='no-click' style='background-color:white;border-radius:10px;border:none;'><span style='color:grey;'>0</span></button>"; }
	return "<button class='no-click' style='background-color:lightgrey;border-radius:10px;border:none;'><span style='color:white;'><b>" + value + "</b></span></button>";
}

function has_hgmd(value) {
	if (value) {
		var list = value.split(";");
		var b = ' ';
		if (list[0] == '2') {
			b = "<button onClick='viewGeneHgmdInfos(\"" + list[1] + "\")' style=\"font-size:8px\"><img src='/icons/Polyicons/12-em-check.png'></button>";
		}
		else if (list[0] == '1') {
			b  = "<button onClick='viewGeneHgmdInfos(\"" + list[1] + "\")' style=\"font-size:8px\"><img src='/icons/Polyicons/12-em-check.png'><span style='color:red;'><i>New!</i></span></button>";
			//b += "<br>";
			//b = "<button onClick='viewGeneHgmdInfos(\"" + list[1] + "\")' style=\"font-size:8px\"><span style='color:red;'><i>New Var!</i></span></button>";
		}
		return b;
	}
	return ' ';
}

function formatHgmd(value) {
	if (value) {
		if (value == 'N.A.') {
			return '<i>'+value+'</i>';
		}
		var list = value.split(";");
		if (list.length == 3 && list[2] == 'New!') {
			return "<button onClick=\"" + list[1] + "\" style=\"font-size:8px\">" + list[0] + " <span style='color:red;'><i>New!</i></span></button>";
		}
		return "<button onClick=\"" + list[1] + "\" style=\"font-size:8px\">" + list[0] + "</button>";
	}
	return ' ';
}

//fonction pour connaitre le numero de la ligne
function getRow(inRowIndex){
	return ' ' + (inRowIndex+1);
}
					
function formatMethods(value){
       
	if (!value) value =0 ;

	switch(value){
		case 0:
		return "<img src='/icons/Polyicons/bullet_red.png'>";
	}
	return "<img src='/icons/Polyicons/bullet_green.png'>"+value;
}  

function formatHoRegionUcsc(url) {
	if (url) return "<a href='"+url+"' target='_blank'>UCSC</a>";
	return ' ';
}

function formatHoRegionEnsembl(url) {
	if (url) return "<a href='"+url+"' target='_blank'>Ensembl</a>";
	return ' ';
}

function formatDejavu(value){
	if (value =="-") return value;
	var mySplitResult = value.split("|");
	var varId;
	if (mySplitResult[1].match(/:/g)) {
		var tmp = mySplitResult[1].split(":");
		varId = tmp[0];
	}
	else { varId = mySplitResult[1]; }
	var url = "polydejavu/dejavu.html?input="+varId;
	return "<a href='"+url+"' target='_blank'>"+mySplitResult[0]+"</a>";
}

function formatReferenceGene(value){
	if (!value) return "-";
	var mySplitResult = value.split(";");
	var string_return="";
	for (i=0;i<mySplitResult.length;i++){
		string_return += '<a href="http://www.ensembl.org/Homo_sapiens/geneview?gene='
			+ mySplitResult[i]
			+ '" target="_blank">'
			+ mySplitResult[i] + '</a>'+"<br>";
	}
	return string_return;
} 

function formatCnv(value){
	if (value.indexOf("Var") == -1) {
		return value;
	}
	return '<a href="http://projects.tcag.ca/cgi-bin/variation/xview?source=hg18&view=variation&id='
	  + value
	  + '" target="_blank">'
	  + value + '</a>';
} 

function formatDgv(value){
	if (!value) return "";
	var ref = "<a href='http://projects.tcag.ca/cgi-bin/variation/xview?source=hg18&view=variation&id="+value+"' target=_blank>";
	ref += "<img  width=8 height=8 src='/icons/Polyicons/8-em-check.png'></a>";
	return ref;
} 

function formatSnp(value){
	if (value.indexOf("rs") == -1) { return value; }
	return '<a href="http://www.ncbi.nlm.nih.gov/snp/?term='
	  + value
	  + '" target="_blank">'
	  + value + '</a>';
} 

function formatFilter(value){
	var img = "<img width=10 height=10 src='/icons/Polyicons/services.png'>";
	if (value == "medium") value = 1;
	else if (value == "strong") value = 2;
	else if (value == "weak") value = 0;
	switch(value){
		case 0:
			return img;
		case 1:
			return img+img;	
		case 2:
			return img+""+img+""+img;
	}
}

function formatValidWithFilter (tvalue) {
	if (tab_id_patientFilter >= 0) {
		value = tvalue[tab_id_patientFilter];
	}
	else { value = tvalue; }
	if (!value) value =0;
	console.log(value);
	return formatValid(value);
}

 function formatValid(value){
 	if (!value) value ="0";
	var text ="";
	if (value > 50){
		value = value - 100;
		text = "sanger";
	}
	value += "";
 	var tab = value.split("+");
	var resume ="";
	if (tab.length > 1) {
		resume = "<br>("+tab[1]+")";
	}
	if (tab[0] == 0 )return "<img src='/icons/Polyicons//bullet_black.png'>";
	if (tab[0] == -1) return "<img src='/icons/Polyicons/bullet_red.png'>";	
	if (tab[0] == -5) return "<img src='/icons/Polyicons/bullet_red.png'>"+text;	
	if (tab[0] == -2) return "<img src='/icons/Polyicons/bullet_purple.png'>"+resume;
	if (tab[0] == 1) return "<img src='/icons/Polyicons/bullet_green.png'>"+text;
	if (tab[0] == 2) return "<img src='/icons/Polyicons/bullet_green.png'>(he)"+text;
	if (tab[0] == 3) return "<img src='/icons/Polyicons/bullet_green.png'>(ho)"+text;
	if (value == -2) return "<img src='/icons/Polyicons/bullet_orange.png'>";
	if (value == -3) return "<img src='/icons/Polyicons/bullet_orange.png'>";
	if (value == -99) return "<img src='/icons/Polyicons/System-Security-Question-icon-16.png'>";
	return "-";
}

function formatPolyphen1(value){
	//value= parseInt(value);
	alert("cocuou");
	return "1";
	value = value *1.0;
	if ( value == 6 ) return "<img src='/icons/Polyicons/bullet_orange_grey.png'>"+"<img src='http://mendel.necker.fr/icons/myicons2/bullet_orange.png'>"+"<img src='http://mendel.necker.fr/icons/myicons2/bullet_orange_grey.png'>";
	if ( value == 7 ) return "<img src='/icons/Polyicons/bullet_red_grey.png'>"+"<img src='http://mendel.necker.fr/icons/myicons2/bullet_red_grey.png'>"+"<img src='http://mendel.necker.fr/icons/myicons2/bullet_red.png'>";;
	if ( value == 5 ) return "<img src='/icons/Polyicons/bullet_green.png'>"+"<img src='http://mendel.necker.fr/icons/myicons2/bullet_green_grey.png'>"+"<img src='http://mendel.necker.fr/icons/myicons2/bullet_green_grey.png'>";;
	if ( value == 3 ) return "<img src='/icons/Polyicons/bullet_grey.png'>"+"<img src='http://mendel.necker.fr/icons/myicons2/bullet_grey.png'>"+"<img src='http://mendel.necker.fr/icons/myicons2/bullet_grey.png'>";
	if ( value == 9 )  return "<img  width=16 height=16 src='/icons/Polyicons//stop_hunabkuc_software.png'>";
	if ( value == 0 )  return "<img  width=16 height=16 src='/icons/Polyicons/cancel.png'>";
	return "-";
}

function formatAlign (value){
	if (value == -1 ) return "<img src='/icons/Polyicons/bullet_red.png'>";
	return "<img src='/icons/Polyicons/bullet_green.png'>";
}

function formatHomozygote (value){
	 
	if ( value <0 ) return "<img src='/icons/Polyicons//bullet_grey.png'>"+"<img src='/icons/Polyicons/bullet_grey.png'>";
	if ( value == 1 ) return  "<img src='/icons/Polyicons//bullet_green.png'>"+"<img src='/icons/Polyicons//bullet_green.png'><br>ho";
	if ( value == 2 ) return "<img src='/icons/Polyicons//bullet_red.png'>"+"<img src='/icons/Polyicons//bullet_green.png'><br>he";
	if ( value == 3 ) return "<img src='/icons/Polyicons//bullet_purple.png'>"+"<img src='/icons/Polyicons//bullet_purple.png'><br>both";
    return  "<img src='/icons/Polyicons//bullet_grey.png'>"+"<img src='/icons/Polyicons//bullet_grey.png'>";
	
}

function formatHapmap (value){
	if ( value == 0 ) return "";
	
	var ref = "<a href='http://hapmap.ncbi.nlm.nih.gov/cgi-perl/gbrowse/hapmap27_B36/?name=SNP:"+value+"' target=_blank>";
	ref += "<img  width=8 height=8 src='/icons/Polyicons//8-em-check.png'></a>";
	
	 return ref;
}

function formatSift(value){
	if (!value) value ="0";
	
	value += "";
	var tab = value.split("+");
	var score ="";
	if (tab.length > 1) {
		score = "<br>("+tab[1]+")";
	}
	if ( tab[0] == 5 ) return  "<img  width=16 height=16 src='/icons/Polyicons/exclamation.png'>";
	if ( tab[0] == 4 ) return  "<img  width=16 height=16 src='/icons/Polyicons/error.png'>";
	if ( tab[0] == 2 ) return "<img  width=16 height=16 src='/icons/Polyicons/difficulty_2_1.png'>"+score;
	if (  tab[0] == 1 ) return "<img  width=16 height=16 src='/icons/Polyicons/difficulty_0.png'>"+score;	
	if (  tab[0] == 0 )   return "<img  width=16 height=16 src='/icons/Polyicons/difficulty_3_b.png'>";
	return "<img  width=16 height=16 src='/icons/Polyicons/System-Security-Question-icon-16.png'>";
	return "-";//"<img src='/icons/myicons/bullet_black.png'>";
}

function formatPolyphen(value){
	if (!value) value = "0";
	value += "";
	var tab = value.split("+");
	var score = "";
	if (tab.length > 1) { score = "<br>("+tab[1]+")"; }
	if ( tab[0] == 2 ) return "<img  width=16 height=16 src='/icons/Polyicons/difficulty_1.png'>"+score;
	if ( tab[0] == 3 ) return "<img  width=16 height=16 src='/icons/Polyicons/difficulty_2_1.png'>"+score;
	if ( tab[0] == 1 ) return "<img  width=16 height=16 src='/icons/Polyicons/difficulty_0.png'>"+score;	
	if ( tab[0] == 0 ) return "<img  width=16 height=16 src='/icons/Polyicons/difficulty_3_b.png'>";
	if ( tab[0] == 4 ) return  "<img  width=16 height=16 src='/icons/Polyicons/error.png'>";
	if ( tab[0] == 5 ) return  "<img  width=24 height=24 src='/icons/Polyicons/exclamation.png'>";
	//if ( value == 0 )  return "<img  width=16 height=16 src='/icons/myicons/warning.png'>";
	//if ( value == 0 )  
	return "<img  width=16 height=16 src='/icons/Polyicons/System-Security-Question-icon-16.png'>";
	return "-";//"<img src='/icons/myicons/bullet_black.png'>";
}

function formatCachePolyPhen (value){
	return "-";
	if (value == "-") return value;
	if (value == "") return value;
	return "<a href='"+url_query+"?"+value+"' target='_blank'>view</a>";
}
 
function formatPolyphenUrl (value){
	return "-";
	if (value == "-") return value;
	var reg=new RegExp("&lt;", "g");
	var g=value.replace(reg,"<");
	return "-";
	return '<form method="post" id ="harvard" action="http://genetics.bwh.harvard.edu/cgi-bin/ggi/ggi.cgi" target="_blank" enctype="application/x-www-form-urlencoded">'+g+'</form>';
}

function formatSimple (value){
	return value;
}

function formatProteinName (value){
	var tab = value.split("+");
	var eurl = '<a href="http://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p='
	  + tab[0]
	  + '" target="_blank">'
	  + tab[0] + '</a>';
	return eurl+"<br>("+tab[1]+")";
}

function formatTranscriptName (value){
	var tab = value.split("+");
	var eurl = '<a href="http://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t='
	  + tab[0]
	  + '" target="_blank">'
	  + tab[0] + '</a>';
	return eurl+"<br>("+tab[1]+")";
}

function formatDisease2(val) {
	if (val == 2) return "<img src='/icons/Polyicons/pill2.png'>";
	if (val == 1) return "<img src='/icons/Polyicons/bullet_green.png'>";
	return ' ';
}

function formatDisease(val) {
	if (val == 2)  return "<img src='/icons/Polyicons/pill2.png'>";
	return "<img src='/icons/Polyicons/bullet_green.png'>";
}
	
function formatChild2(val) {
	if (val == "C21")  return "<img src='/icons/Polyicons/baby-girl-s.png'>";
	if (val == "C22")  return "<img src='/icons/Polyicons/baby-girl-d.png'>";
	if (val == "C11")  return "<img src='/icons/Polyicons/baby-boy-s.png'>";
	if (val == "C12")  return "<img src='/icons/Polyicons/baby-boy-d.png'>";
	if (val == "F1")  return "<img src='/icons/Polyicons/male-s.png'>";
	if (val == "F2")  return "<img src='/icons/Polyicons/male-d.png'>";
	if (val == "M1")  return "<img src='/icons/Polyicons/female-s.png'>";
	if (val == "M2")  return "<img src='/icons/Polyicons/female-d.png'>";
	return "<img src='/icons/Polyicons/12-em-check.png'>";
}	
	
function formatChild(val) {
	if (val == "C2")  return "<img src='/icons/Polyicons/baby-girl.png'>";
	if (val == "C1")  return "<img src='/icons/Polyicons/baby-boy.png'>";
	if (val == "F")  return "<img src='/icons/Polyicons/male.png'>";
	if (val == "M")  return "<img src='/icons/Polyicons/female.png'>";
	if (val == "1")  return "<img src='/icons/Polyicons/male.png'>";
	if (val == "2")  return "<img src='/icons/Polyicons/female.png'>";
	return "<img src='/icons/Polyicons/12-em-check.png'>";
}	

function formathe(val) {
	if (val == 0)  return "<img src='/icons/Polyicons/12-em-cross.png'>";
	return "<img src='/icons/Polyicons/12-em-check.png'>";
}
	
function formatho(val) {	
	if (val == 1) return "<img src='/icons/Polyicons/12-em-cross.png'>";;	
	return "<img src='/icons/Polyicons/12-em-check.png'>";	
}
	
function formatInclude(value){
	value= value*1.0;
	if (!value) value =0 ;
	switch(value){
		case 2:{
			return "<center><img src='/icons/Polyicons/trashcan-full-icon.png'></center>";
		}
		case -2:{
			return "<center><img src='../images/polyicons/del_reg.png'></center>";
		}
		case -1:{
			return "<center><img src='/icons/Polyicons//delete.png'></center>";
		}
		case 0: {
			return "<center><img src='/icons/Polyicons/add.png'></center>";
		}
		case 3: {
			return "<img src='/icons/Polyicons/12-em-check.png'>";
		}
	}
	return "<center><img src='/icons/Polyicons/bullet_grey.png'></center>";
}   

function icon_server(value) {
	if (value == "dev"){
		return "<img src='/icons/Polyicons/bullet_purple.png'>";
	}
	return "<img src='/icons/Polyicons/bullet_green.png'>";
	
}

function icon_follow(value) {
	if (value == 1){
		return "<img src='/icons/Polyicons/bullet_green.png'>";
	}
	return value;
}

function formatColorCadd(value) {
	var new_text;
	if (parseInt(value) >= 30) { new_text = "<font color='red'>" + value + "</font>"; }
	else if (parseInt(value) >= 20) { new_text = "<font color='orange'>" + value + "</font>"; }
	else if (parseInt(value) >= 10) { new_text = "<font color='green'>" + value + "</font>"; }
	return new_text;
}

function formatColorPhenotypes(value) {
	var new_text = "<font color='red'>" + value + "</font>";
	return new_text;
}

function formatColorProjectAndPhenotypes(value) {
	var tab = value.split(";");
	var new_text = tab[0];
	if (tab.length > 1) { 
		new_text += "<br><font color='red'>" + tab[1] + "</font>";
	}
	return new_text;
}

function formatGnomadAC(value) {
	if (value == '.') { return value; }
	var new_text;
	if (parseInt(value) <= 10) { new_text = "<font color='red'>" + value + "</font>"; }
	else if (parseInt(value) <= 100) {  new_text = "<font color='orange'>" + value + "</font>";}
	else if (parseInt(value) <= 1000) {  new_text = "<font color='green'>" + value + "</font>";}
	else { new_text = value; }
	return new_text;
}

function formatConsequence(values){
	var new_tab = [];
	var tab = values.split(";");
	var new_text;
	var text = tab.join("<br>");
	if (tab[0].match(/Stop/)) { new_text = "<font color='red'>" + text + "</font>"; }
	else if (tab[0].match(/Acc/)) { new_text = "<font color='red'>" + text + "</font>"; }
	else if (tab[0].match(/Mature/)) { new_text = "<font color='red'>" + text + "</font>"; }
	else if (tab[0].match(/ncRNA/)) { new_text = "<font color='orange'>" + text + "</font>"; }
	else if (tab[0].match(/Region/)) { new_text = "<font color='orange'>" + text + "</font>"; }
	else if (tab[0].match(/Missense/)) { new_text = "<font color='orange'>" + text + "</font>"; }
	else if (tab[0].match(/No-/)) { new_text = "<font color='orange'>" + text + "</font>"; }
	else if (tab[0].match(/Frameshift/)) { new_text = "<font color='red'>" + text + "</font>"; }
	else { new_text = "<font color='green'>" + text + "</font>"; }
	return new_text;
}

function formatVarId(value){
	if (value.match(/ERROR/)) { return "<img src='http://mendel.necker.fr/icons/myicons2/cross.png'>" + " " + value; }
	else if (value == "No result...") { return "<img src='http://mendel.necker.fr/icons/myicons2/error.png'>" + " " + value; }
	return value;
}

function formatCellGeneTrans(value){
	if (/; /.test(value)) {
		var tab = value.replaceAll(";", ",");
		return tab;
	}
	return value;
}

function formatCellPat(value){
	var tab = value.split("|");
	return tab[0];
}

function formatCellDejaVU(value){
	if (value == '0') { return '0'; }
	var tab = value.split("; ");
	return tab.length;
}

function formatCellPatDetails(value){
	//alert("VALUE: " + value);
	return hashPatient[value];
}

function formatRsName(value){
	var eurl = value;
	if (value.match(/;/)) {
		var tab = value.split(";");
		var id = tab[0];
		var tabIdFields = id.split("_");
		var alleles = tabIdFields[2] + '/' + tabIdFields[3];
		if (tab[1].match(/TMP/)) { eurl = '<a href="http://www.ensembl.org/Homo_sapiens/Variation/Summary?v=' + tab[1] + '" target="_blank">' + tab[1] + '</a> (' + alleles + ')'; }
		else if (tab[1].match(/rs/)) { eurl = tab[0] + ' / <a href="http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=' + tab[1] + '" target="_blank">' + tab[1] + '</a> (' + alleles + ')'; }
		else { eurl = tab[0]; }
	}
	return eurl;
}

function formatShowYourVariations(value) {
	if (value == "X") { return "<img src='/icons/Polyicons/accept-icon-16.png'>"; }
	else { return "<img src='/icons/Polyicons/Cancel-icon-16.png'>"; }
	return value;
}

function formatEmail(value) {
	if (value.match(/;/)) {
		var tab = value.split(";");
		return tab.join('<br>');
	}
	return value;
}

function statusIcon(value) {
	if (value == 'OK') { return "<img src='./images/polyicons/accept-icon-16.png'>"; }
	else if (value == '?') { return "<img src='./images/polyicons/interrogation.png'>"; }
	else if (value == 'ERROR') { return "<img src='./images/polyicons/Cancel-icon-16.png'>"; }
	return "<img src='./images/Polyicons/problem.png'>";
}

function dejavu_capture_diag (value) {
	//if (value == 'yes') { return "<img src='/icons/Polyicons/0169_star2.png'>"; };
	if (value) {
		return "<button onClick=\"showCapturesFromTranscriptsIds('" + value + "')\" style=\"font-size:8px\"><img src='/icons/Polyicons/12-em-check.png'></button>";
	};
	return '';
}

function dejavu_genes (value) {
	if (value == 'yes') { return "<font color='green'>&#8730;</font>"; }
	else if (value == 'NEVER') { return "<font color='orange'>NOT IN REF</font>"; }
	else { return "<font color='red'>X</font>"; }
}

function formatRegionHo(this_value){
	if (this_value == null) {
		return "";
	}
	var value = this_value.split(';');
	var listPat = value[1].split('|');//value.length
	return "<button data-dojo-type=\"dojox/mobile/Button\" class=\"blueButton\" style=\"width:auto\" onclick=\"selectPatients_horegion('" + listPat.join(',') + "')\"><span style=\"color:white\">" + value[0] + "</span></button>";
}

function formatRegionRec(this_value){
	if (this_value == null) {
		return span = "<button type='button' style=\"font-size:08px;\">0</button>";
	}
	var value = this_value.split(';');
	var listFam = value[1].split('|');//value.length
	return "<button data-dojo-type=\"dojox/mobile/Button\" class=\"blueButton\" style=\"width:auto\" onclick=\"selectPatients_horegion('" + listFam.join(',') + "')\"><span style=\"color:white\">" + value[0] + "</span></button>";
}

function formatRegionHo2(val) {
	if (val == '')  return " ";
	if (val == 'false')  return " ";
	return "<button onClick='show_region_ho(\"" + String(val) + "\");' showLabel='true' style='width:40px;float:left;'><img src='/icons/Polyicons/12-em-check.png'></button>";
}

function formatColorCell(value) {
	var lFields= String(value).split(';');
	return '<span style="color:' + lFields[1] + '">' + lFields[0] + '</span>';
}

function format_displayButtonWithCmd(value) {
	var lTmp = value.split(';');
	var button_name = lTmp[0];
	var cmd_name = lTmp[1];
	return "<button type='button' style=\"font-size:08px;\" onclick=\"" + cmd_name + "\">" + button_name + "</button>";
}

function format_displayButtonWithCmd_v2(value) {
	var lTmp = value.split(';');
	var button_name = lTmp[0];
	var cmd_name = lTmp[1];
	cmd_name = cmd_name.replace("view_web_igv_bam('dialog_igv', ", "view_web_igv_bam_simple(");
	return "<button type='button' style=\"font-size:08px;\" onclick=\"" + cmd_name + "\">" + button_name + "</button>";
}
