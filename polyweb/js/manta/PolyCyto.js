
var dataStorePlotValue;
var gridPlot;
var tabPos ="";
var tabRatio = "";

var nb=0;
require([ "dojo/_base/xhr","dojox/widget/DialogSimple", "dojox/widget/Dialog","dijit/MenuBar", "dijit/MenuBarItem", "dijit/PopupMenuBarItem","dijit/DropDownMenu", "dijit/MenuItem"]);


function LaunchPolycyto(project,patient)
{
	var url_plot = "PolyCyto.html?projectname=" + project + "&filename=" + patient +"&trans=yes";
	var myWindow = window.open(url_plot,"_blank",""); 
}


function onload_polycyto_body() {
	setTitrePolyCyto();
	GetAllCNV();
	GetSVeq();
	GetChrsPloidy();
}

function launchAnalyse(transloc)
{
	setTitrePolyCyto();
	
	if(transloc == "yes")
	{
		document.getElementById("id_transloc").style.visibility = 'visible';
		GetSVeq();
	}
	GetAllCNV();
	GetChrsPloidy();
	GetChrsDisomie();
	refresh_grids();
}

function setTitrePolyCyto()
{
	//document.getElementById("waiting").hide;
	var parameters = location.search.split("&");
   	var project;
   	var patient;
   	var in_frame = 0;
	for (var i=0; i<parameters.length; i++) {
		parameters[i] = parameters[i].replace('?','');
		var par = parameters[i].split("=");
		if (par[0] == 'projectname') { project = par[1]; }
		if (par[0] == 'filename') { patient = par[1]; }
		if (par[0] == 'iframe') { in_frame = par[1]; }
	}
	if (in_frame == 1) 
	{ 
		titre = '<br>'; 
	}
	else {
		titre = '<br><b><font size="3" style="color:DarkOrange">' + project + '   /   ' + patient + '</font></b><br>';
		document.getElementById('titre').innerHTML=titre;
	}
}

//---------------------------------------
// pour les CNV
//---------------------------------------
var dataStore_SVCompare;
var gridSVCompare;
var url_CNV_Patient = url_path + "/manta/PolyCytoRetrieveCNV.pl";


var h_gridSVCompare = {};
function GetAllCNV() 
{
		document.getElementById("span_cnv_dup_del").innerHTML = "<img src='../../images/polyicons/wait18trans.gif'> CNV (DUP/DEL)";

		var parameters = location.search.split("&");
		var arg1 = parameters[0].split("=");
   		var arg2 = parameters[1].split("=");
   			
    	var project = arg1[1].split(",");
    	var patient = arg2[1].split(",");
    	
		var minlength = document.getElementById('minlengthSelect').value;
		var transmission = document.getElementById('transmissionSelect').value;
		var chr = document.getElementById('chromSelect').value;
		var dejavu = document.getElementById('dejavu').value;
		var thevalue = document.getElementById('scoreEvent').value;
		
		var genes = document.getElementById('geneSelect').value;
		if (genes == "") { genes = "all"};
		
		var omim=0;
		if (document.getElementById('OMIMSelect').checked == true)
		{
			omim=1;
		}
		
		
		// var genotype = document.getElementById('genotypeSelect').value;
		//var thevalue= 	document.querySelector('input[name=cnvSelect]:checked').value;
		//var thevalue = document.getElementById('scoreCNVSelect').value;
		//var maxlength = document.getElementById('maxlengthSelect').value;	
		//var maxfreq = document.getElementById('maxfreqSelect').value;
		//var cytoband = document.getElementById('cytobSelect').value;
		//if (cytoband == "") { cytoband = "all"};
		//if (thevalue == 1)
		//{
		//	accessCNVFilters("Off");
		//}
		
		
		// ceux que l'on pourra supprimer
		
		var genotype="both";
		var maxlength="nomax";
		var cytoband="all";
		var maxfreq="nomax";
	
			
		dataStore_SVCompare = new dojo.data.ItemFileReadStore({ url: url_CNV_Patient+ "?project=" + project + "&patient=" + patient + "&minlength=" + minlength + "&maxlength=" + maxlength + "&transmission=" + transmission + "&maxfreq=" + maxfreq + "&dejavu=" + dejavu + "&genotype=" + genotype  + "&select_best=" + thevalue + "&chrom=" +chr +"&cytoband=" +cytoband + "&genes=" +genes + "&omim=" +omim + "&print=" + 0 });
		
		h_gridSVCompare[patient] = dijit.byId("SVCompare_gridDetails");
		h_gridSVCompare[patient].setStore(dataStore_SVCompare);
		h_gridSVCompare[patient].store.fetch({
			onComplete: function(items){
				document.getElementById("span_cnv_dup_del").innerHTML = "CNV (DUP/DEL)";
			}
		});
	
		return;
}

function formaterNB(value)
{
	// pour activer le bouton a la fin du script CNV
	if (value == 1)
	{
			//document.getElementById("btnAbstract").style.visibility = 'visible'; 
			//document.getElementById("btnAbstract").disabled = false; 
	}
	return (value);
}
var color_del = "#D62C1A";
var color_dup = "#1079B2";

function formaterDUPDEL(value)
{

	var eurl;
	if (value == "-")
	{   
		return value;
	}
	var text ='';
	var color = return_color_by_type(value);
	 if ( /DEL/.test(value)) {
	 	//eurl = '<b><font style="color:red"> '+ value +' </font></b>';
	 	text = "DEL";
	 	
	 	//eurl=  '<button type="button" class="btn btn-xs btn-primary " style="border-color:black;padding:5px;background-color:'+color_del+';color:white;font-size: 14px;font-family:  Verdana;border-color:white;">DEL'+'</button>';
	 	//eurl = '<p style="font-size:10px;color:white;background-color:rgb(255,74,42);padding:8px">' + value + ' XXX</p>';
	 }
	 
	 if ( /DUP/.test(value)) {
		 text = "DUP";
	 	//eurl = '<b><font style="color:blue"> '+ value +' </font></b>';
	 	eurl=  '<button type="button" class="btn btn-xs btn-primary " style="padding:5px;background-color:'+color_dup+';color:white;font-size: 14px;font-family:  Verdana;border-color:white;">DUP'+'</button>';
	 	//eurl = '<p style="font-size:10px;color:white;background-color:rgb(12,148,230);padding:8px">' + value + ' </p>';
	 }
	 
	 if ( /TRI/.test(value)) {
			//eurl = '<b><font style="color:blue"> '+ value +' </font></b>';	 
				 text = "TRI";
			eurl = '<p style="font-size:10px;color:white;background-color:rgb(12,98,228);padding:8px">' + value + ' </p>';	
	}
	 
	 if ( /QUA/.test(value)) {
		 	 text = "QUA";
	 	//eurl = '<b><font style="color:indigo"> '+ value +' </font></b>';
	 	eurl = '<p style="font-size:10px;color:white;background-color:rgb(1,26,255);padding:8px">' + value + ' </p>';	
	 }
	 return return_style_value(color,value,11);
	 return '<button type="button" class="btn btn-xs btn-primary " style="background-color:'+color+';color:white;font-size: 8px;font-family:  Verdana;border-color:white;">'+text+'</button>';
	return eurl;
}


function formaterRatio(value)
{
		var out;
		var val;
		var eurl="-";
		
		
			if (value == "-")
		{
			return "-";
		}
		
	var color;
	var val = parseInt(value);

	if (val == 0 ){
		color="green";
	}
	else  if (val < 0){
		color="blue";
	}
	else {
		 color="tomato";
	}
	
	return '<button "type="button" class="btn btn-xs btn-primary " style="border-color:'+color+';background-color:white;color:'+color+';font-size: 11px;font-family:  Verdana">'+val.toFixed(2)+'</button>';
	
		
		
}

function formaterScoreDude(value)
{
		var eurl="-";
		
		if (value == "-")
		{
			return value;
		}
				
		if (value == -1)
		{
			return eurl;
		}
		
		eurl = '<b><font style="color:DarkOrange">' + value.toFixed(2) + '</font></b>';
		return eurl;	
}	
	
function formaterDUPSEG(value)
{
	var color;
	value = value +0;
	if (value == "-")
	{
			return "-";
	}
	var val = parseInt(value);
	
	if (val == 0 ){
		color="green";
	}
	else  if (val < 40){
		color="orange";
	}
	else {
		 color="tomato";
	}
	 return return_style_value (color,val,8);
	return '<button type="button" class="btn btn-xs btn-primary " style="background-color:'+color+';color:white;font-size: 8px;font-family:  Verdana;border-color:white;">'+val+'%</button>';
	
	return '<button "type="button" class="btn btn-xs btn-primary " style="border-color:'+color+';background-color:white;color:'+color+';font-size: 11px;font-family:  Verdana">'+val+'%</button>';
	
	
}

	
function formaterCN(value)
{
		var out = "";
		var val;
		
		if (value == "-")
		{
			return value;
		}
		
	
		if (! /[CN]/.test(value))
		{
			return "-";
		}
		else
		{
			var tval = value.split(" ");
			var len = tval.length;
			
			for ( i=0 ; i <len-1 ; i++ )
			{
				if ( /[CN]/.test(tval[i]) )
				{
						out = out + tval[i] +" ";
				}
			}
		}
		
		out = out.replace(/CN0/g, '<font style="color:red">  CN0 </font>');
		out = out.replace(/CN2/g, '<font style="color:blue">  CN2 </font>');
		out = out.replace(/CN3/g, '<font style="color:blue">  CN3 </font>');
		out = out.replace(/CN4/g, '<font style="color:indigo">  CN4 </font>');
		
		out = '<b>'+out+'</b>';
		return out;	
}

function formaterQUAL(value)
{
		var out = "";
		var val;

		var tvalue = value.split(" ");
		var val = tvalue[0];
		var qual = tvalue[1];
	
		var tval = val.split("/");
		var len = tval.length;
		
		var qwc = 0;
		var qc = 0;
		var qm = 0;
			
		for ( i=0 ; i <len-1 ; i++ )
		{
			var tval_caller = tval[i].split(":");
			val = tval_caller[1];
			val = val*1;
			
			if(tval_caller[0] == "wisecondor")
			{
				qwc = val.toFixed(2);
			}
			
			if(tval_caller[0] == "canvas")
			{
				qc = val;
			}
			 
			 if(tval_caller[0] == "manta")
			{
				qm = val;
			}
		}
		
		qual = qual * 1;
	
		
		if ( qual == 1 )
		{
				qual = '<font size = "2" color="green"><span class="glyphicon glyphicon-plus-sign"></span></font>';
		}
		else
		{
			if (qual >= 0.5)
			{
				qual = '<font size = "2" color="orange"><span class="glyphicon glyphicon-plus-sign"></span></font>';
			}
			else
			{
				qual = '<font size = "2" color="red"><span class="glyphicon glyphicon-minus-sign"></span></font>';
			}
		}
		
		out =  qwc + "/ " + qc + "/ " + qm;
		out = '<a href="#" data-toggle="tooltip" title="' + out + '">' +  qual  +'  </a>';
		return out;	
}

function formaterCHROM_int(value)
{
	var eurl = value;

	if (value == "-")
		{
			return value;
	}
	var v1 = value.split(";");
	var color = return_color_by_type(v1[0]); 
	value = v1[1];
	if (value == 23)
	{
		 eurl ="X";
	}
	if (value == 24)
	{
		 eurl="Y";
	}
	
	eurl = '<b>' + eurl + '</b>';
	return return_style_value(color,value);
}

function return_color_by_type(value){
	var color = "grey";
	 if ( /DEL/.test(value)) {
	 	//eurl = '<b><font style="color:red"> '+ value +' </font></b>';
	 	color ="#D62C1A";
	 }
	 
	 if ( /DUP/.test(value)) {
	 	color ="#1079B2";
	 }
	 
	 if ( /TRI/.test(value)) {
	 		color ="rgb(12,98,228)";
	}
	 
	 if ( /QUA/.test(value)) {
	 		color ="rgb(1,26,255)";
	 }
	 if ( /TL/.test(value)) {
	 		color ="#00AB78";
	 }
	  if ( /INV/.test(value)) {
	 		color ="#42BDCB";
	 }
	 return color;
	
}

function formaterCytoBand(value)
{
	var out=" ";
	var outshort;
	if (value == "-")
	{
		return value;
	}
	
	var v1 = value.split(";");
	var color = return_color_by_type(v1[0]); 
	var names = v1[1].split(",");
	var len = names.length;

	
	for ( i = 0 ; i < len; i++)
	{
				if (names[i] != "-")
				{	
					out = out +"  " + names[i] ;
				}
	}
	
	outshort = names[0]+" ... "+names[len-1];
	
	if (len < 3 )
	{
			out = ''+ out + '';
	}
	else
	{
			out = '<a href="#" data-toggle="tooltip" title="' + out + '"> ' + outshort + '</a>';
	}
	return return_style_value(color,out);
	//return out;
}	

function string_length(value) {
	var eurl;
	
	 if (value > 1000000 ) {
	 	value = value / 1000000;
	 	eurl = '<b>' + value.toFixed(2)  + ' Mb </b>';
	 }
	 
	 else if (value > 1000 ) {
	 	value = value / 1000;
	 	eurl = value.toFixed(2)  + ' kb';
	 }	
	 
	 else if (value < 1000 ) {
	 	eurl = value  + " pb";
	 }	
	 return eurl;
}

function formaterLocus(value)
{
	let obj = JSON.parse(value);

	var color = return_color_by_type(obj.type); 
	var out=" ";
	var l = obj.len;
	
	 
	out += '<font style="color:white">' + obj.start + '</font><font style="color:white">-</font><font style="color:white">' + obj.end + '</font><br>'+obj.cytoband+" ("+l+")";
	
	var text =  return_style_value(color,out);
	var out2 = "<div>"+formater_viewers1(obj.viewers)+' <span style=\"white-space: nowrap;position: relative; top: -5px; left: -5px;\">'+text+"</span><div>";
	
	return out2;
}	

function formater_viewers1(obj)
{
   			nb++;
   			var boutonID = 'boutonIGVCNV'+ nb;		
   				
   			var igv_value = obj.igv.locus_start + ';' + obj.igv.locus_end + ';' + obj.igv.bam_files + ";" + obj.igv.bam_names + ";"+ obj.igv.l;	
			var onvIGV = ' onclick="launch_web_igv_start_and_end(\'' + obj.igv.id + "et" + igv_value + '\')"';		
			var bigv  = return_simple_button_with_onclick("IGV","#F0A30A",onvIGV,"white");
			var onclick = 'onclick="launch_plot(\''+ obj.plot.transmission + '\',' + obj.plot.type + ',\'' + obj.plot.chromosome + '\',' + obj.plot.start + ',' + obj.plot.end + ' )" ';
			var bplot   = return_simple_button_with_onclick("SNP Array","#A20025",onclick,"white");
			
			//var bigv   = '<button class="igvIcon2"  id="' + boutonID + '" onclick="launch_web_igv_start_and_end(\'' + obj.igv.id + "et" + igv_value + '\')"></button>';
			//var bdgv = '<a href="' + obj.dgv.url  + '" target="_blank"><button class="btn btn-default btn-xs" ><b><i><font style="color:red;font-size:10px;padding:0px; ">  DGV  </font></i></b></button></a>';
			var bdgv  =  button_url("DGV",obj.dgv.url,"#3AB795","white");
			var bgnomad = button_url("GNOMAD",obj.gnomad.url,"#217DBB","white");
			out = ' <span style="white-space: nowrap;position: relative; top: -10px; left:0px;">'+bigv+bplot+bgnomad+bdgv+'</span></div>';
	return out;
}

function button_url(text,url,color1,color2){
	var onclickdgv='onclick=window.open(\"' + url + '\",\"_blank\")';
	return return_simple_button_with_onclick(text,color1,onclickdgv,color2);
}

function formaterPlot(value)
{
	if (value == "-")
	{
		return value;
	}

	var val1 = value.split(";");
	
	var type = val1[0]
	var valcnv = val1[1].split(":");
	
	var chr = valcnv[0];
	var pos = valcnv[1].split("-");
	var debcnv = pos[0];
	var fincnv = pos[1];
	
	var transmission = val1[2];

	var out=" - ";
	
	if (type == "DEL") { type=1; }
	if (type == "DUP") { type=2; }
	if (type == "TRI") 	 { type=3; }
	if (type == "QUA") { type=4; }
	
	var boutonID = 'boutonPlot'+ nb;				
   	var out   = '<button class="btn btn-classic btn-m" style="border: 1px solid black;padding:3px" id="' + boutonID + '" onclick="launch_plot(\''+ transmission + '\',' + type + ',\'' + chr + '\',' + debcnv + ',' + fincnv + ' )" style="font-size:16px;">&#128200;</button>';
  
	return out;
}	


function formaterLength(val)
{
	var eurl;
	
	if (val== "-")
	{
		return val;
	}
	var v1 = val.split(";");
	var color = return_color_by_type(v1[0]); 
	var value = v1[1];
	 if (value > 1000000 ) {
	 	value = value / 1000000;
	 	eurl = '<b>' + value.toFixed(2)  + ' Mb </b>';
	 }
	 
	 else if (value > 1000 ) {
	 	value = value / 1000;
	 	eurl = value.toFixed(2)  + ' kb';
	 }	
	 
	 else if (value < 1000 ) {
	 	eurl = value  + " pb";
	 }	
	 return return_style_value(color,eurl);
	 
}


function formater_viewers(value)
{
		let obj = JSON.parse(value);
   			nb++;
   			var boutonID = 'boutonIGVCNV'+ nb;		
   				
   			var igv_value = obj.igv.locus_start + ';' + obj.igv.locus_end + ';' + obj.igv.bam_files + ";" + obj.igv.bam_names + ";"+ obj.igv.l;	
			var onvIGV = ' onclick="launch_web_igv_start_and_end(\'' + obj.igv.id + "et" + igv_value + '\')"';		
			var bigv  = return_simple_button_with_onclick("igv","#2E5283",onvIGV,"#700CBC");
			//var bigv   = '<button class="igvIcon2"  id="' + boutonID + '" onclick="launch_web_igv_start_and_end(\'' + obj.igv.id + "et" + igv_value + '\')"></button>';
			var url_dgv ="";
			var bdgv = '<a href="' + obj.dgv.url  + '" target="_blank"><button class="btn btn-default btn-xs" ><b><i><font style="color:red;font-size:10px;padding:0px; ">  DGV  </font></i></b></button></a>';
			var url_gnomad ="";
			//onclick="window.open(\'' + url + '\', \'_blank\')"
			var bgnomad = '<a href="' + obj.gnomad.url  + '" target="_blank">'+return_simple_button("gnomAD","#C67FAE")+'</a>';
			var bdgv = '<a href="' + obj.dgv.url  + '" target="_blank">'+return_simple_button("DGV","#C67FAE")+'</a>';
			var bgnomad2 = '<a href="' + obj.gnomad.url  + '" target="_blank"><button class="btn btn-default btn-xs" ><b><font style="color:blue;font-size:8px; ">gnom</font><font style="color:black;font-size:10px; ">AD</font></b></button></a>';
			var onclick = 'onclick="launch_plot(\''+ obj.plot.transmission + '\',' + obj.plot.type + ',' + obj.plot.chromosome + ',' + obj.plot.start + ',' + obj.plot.end + ' )" ';
			
			var bplot   = return_simple_button_with_onclick("snp array","#2E5283",onclick,"#700CBC");
			//'<button type="button" class="btn btn-xs btn-primary " style="background-color:blue'+';color:white;font-size: 8px;font-family:  Verdana;border-color:white;"'+onclick+'>Snp Array'+'</button>';

			
			//var bplot   = '<button class="btn btn-classic btn-m" style="border: 1px solid black;padding:3px" " style="font-size:16px;">Snp array</button>';
  			//var bplot   = '<button class="btn btn-classic btn-m" style="border: 1px solid black;padding:3px"></button> '
  			var out = "";
  
			out = ' <span style="white-space: nowrap;">'+bigv+bplot+'</span><span style="white-space: nowrap;"'+bgnomad+bdgv+'</span></div>';
			
	return out;
}		

function return_simple_button_with_onclick(value,color,onclick,color2="white") {
	return '<button type="button" class="btn btn-xs btn-primary " style="background-color:'+color+';color:'+color2+';font-size: 8px;font-family:  Verdana;border-color:white;box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.4);"'+onclick+'>'+value+'</button>';
	}


function formater_infoPATIENTIGV(value)
{
		if (value == "-")
		{
			return value;
		}
		
		var parameters = location.search.split("&");
		var arg1 = parameters[0].split("=");
  	 	var arg2 = parameters[1].split("=");
		
  		 var projectname = arg1[1];
  	  	var patname = arg2[1];
    	
		var out = " ";
	
		var tv =  value.split(';');
		var bamnames = tv[0];
		var bamfiles = tv[1];
		var idCNV = tv[2];
	
		var w = idCNV.split("_");
		var chr = w[1];
				
			// 1  kb autour du start
			var s1 = parseInt(w[2]);
			var start1 = s1 - 500;
			var start2   = s1 + 500;
			var locusstart = 'chr' + chr  + ':' + start1 + '-' + start2;
				
			// 1  kb autour du end
			var s2 = parseInt(w[3]);
			var end1 = s2 - 500;
			var end2   = s2 + 500;
			var locusend = 'chr' + chr  + ':' + end1 + '-' + end2;
   				
   			var t = (s2-s1);
			var taille = 1;
	
	 		if (t > 1000000 ) {
	 				t = Math.round(t / 1000000);
	 				taille = t  + "Mb";
	 		}
	 		else
	 		{
				 if (t > 1000 ) {
	 					t = Math.round(t / 1000);
	 					taille = t  + "kb";
	 			}
	 			else
	 			{
	 					taille = t + "pb"
	 			}
	 		}	
   				
   			nb++;
   			var boutonID = 'boutonIGVCNV'+ nb;					
   			var igv_value = locusstart + ';' + locusend + ';' + bamfiles + ";" + bamnames + ";"+ taille;	
			var etiquette   = '<button class="igvIcon2"  id="' + boutonID + '" onclick="launch_web_igv_start_and_end(\'' + idCNV + "et" + igv_value + '\')"></button>';
			var url_dgv ="";
			var bdgv = '<a href="' + url_dgv  + '" target="_blank"><button class="btn btn-default btn-xs" ><b><i><font style="color:red;font-size:10px;padding:0px; ">  DGV  </font></i></b></button></a>';
			var url_gnomad ="";
			var bgnomad = '<a href="' + url_gnomad  + '" target="_blank"><button class="btn btn-default btn-xs" ><b><font style="color:blue;font-size:8px; ">gnom</font><font style="color:black;font-size:10px; ">AD</font></b></button></a>';
		//	var bplot   = '<button class="btn btn-classic btn-m" style="border: 1px solid black;padding:3px" id="' + boutonID + '" onclick="launch_plot(\''+ transmission + '\',' + type + ',' + chr + ',' + debcnv + ',' + fincnv + ' )" style="font-size:16px;">&#128200;</button>';
  				var bplot   = '<button class="btn btn-classic btn-m" style="border: 1px solid black;padding:3px"></button> '
  			var out = "";
  			out =  '<table class= "table table-sm table-striped table-condensed table-bordered table-primary table_gnomad" style="box-shadow: 1px 1px 6px #337AB7 ;font-size: 7px;font-family:  Verdana;margin-bottom:0px;opacity:0.2" >';
			out += '<tr>';
			out += '<td>'+etiquette+'</td>'+'<td>'+bdgv+'</td>'+'<td>'+bgnomad+'</td>'+'<td>'+bplot+'</td>';
			out += '</tr></table>';
	return out;
}		

function formaterInfoCallerCoverage(value) {
	var color= "#8D71B4";
	if (value === "-")
	{
	var out = '<div style ="opacity:0.2"> wisecondor';
	
	out =  '<table class= "table table-sm table-striped table-condensed table-bordered table-primary table_gnomad" style="box-shadow: 1px 1px 6px #337AB7 ;font-size: 7px;font-family:  Verdana;margin-bottom:0px;opacity:0.2" >';
	out = out + '<caption style="caption-side: top; text-align: center; font-weight: bold;opacity:0.2">wisecondor</caption>';
	out = out + '<tr>';
	out = out + '<th  '+th_style+'>Ratio</th>';
	out = out + '<th '+th_style+'>Score</th>';
	out = out + '</tr>';
	out = out + '<tr>';
	var b = return_simple_button ("-",color);
	out +=  '<td style="vertical-align: middle;">' + b + '</td>';
	b = return_simple_button ("-",color);
	out +=  '<td >' + b + '</td>';
	out = out + '</tr>';
	out = out + '</table></div>';
	return out;
	}
	var tvalue = value.split(";");
	 color = return_color_by_type(tvalue[3]);
	 
	var name = tvalue[0];
//	var qual = tvalue[1];
//	var text= name
	
	var hcolor = "#E8FFF8";
	var th_style = 'style="padding:1px;border: none;text-align: center;"';
	var out = name;
	out =   '<table class= "table table-sm table-striped table-condensed table-bordered table-primary table_gnomad" style="box-shadow: 1px 1px 6px #337AB7 ;font-size: 7px;font-family:  Verdana;margin-bottom:0px" >';
		out = out + '<caption style="caption-side: top; text-align: center; font-weight: bold;">wisecondor</caption>';
	out = out + '<tr>';
	out = out + '<th  '+th_style+'>Ratio</th>';
	out = out + '<th '+th_style+'>Score</th>';
	out = out + '</tr>';
	out = out + '<tr>';
	var b = return_simple_button (tvalue[1],color);
	out +=  '<td style="vertical-align: middle;">' + b + '</td>';
	b = return_simple_button (tvalue[2],color);
	out +=  '<td >' + b + '</td>';
	out = out + '</tr>';
	out = out + '</table>';
	return out;
	
}

function formaterInfoCallerSr(value) {
	if (value === "-")
	{
		var out = "manta";
		out =  '<table class= "table table-sm table-striped table-condensed table-bordered table-primary table_gnomad" style="box-shadow: 1px 1px 6px #337AB7 ;font-size: 7px;font-family:  Verdana;margin-bottom:0px;opacity:0.2" >';
		out = out + '<caption style="caption-side: top; text-align: center; font-weight: bold;opacity:0.2">Manta</caption>';
		out = out + '<tr>';
		out = out + '<th  '+th_style+'>sr</th>';
		out = out + '<th  '+th_style+'>Pr</th>';
		out = out + '<th '+th_style+'>Score</th>';
		out = out + '</tr>';
		out = out + '<tr>';
		var b = return_simple_button ("-",color);
		out +=  '<td style="vertical-align: middle;">' + b + '</td>';
		b = return_simple_button ("-",color);
		out +=  '<td >' + b + '</td>';
		b = return_simple_button ("-","red");
		out +=  '<td >' + b + '</td>';
		out = out + '</tr>';
		out = out + '</table>';
	return out;
	}
		
	
	var tvalue = value.split(";");
	var name = tvalue[0];
//	var qual = tvalue[1];
//	var text= name
	var color= "#8D71B4";
	 color = return_color_by_type(tvalue[4]);
	var hcolor = "#E8FFF8";
var th_style = 'style="padding:1px;border: none;text-align: center;"';
	var out = name;
	out =  '<table class= "table table-sm table-striped table-condensed table-bordered table-primary table_gnomad" style="box-shadow: 1px 1px 6px #337AB7 ;font-size: 7px;font-family:  Verdana;margin-bottom:0px" >';
	out = out + '<caption style="caption-side: top; text-align: center; font-weight: bold">Manta</caption>';
	out = out + '<tr>';
	out = out + '<th  '+th_style+'>sr</th>';
	out = out + '<th  '+th_style+'>Pr</th>';
	out = out + '<th '+th_style+'>Score</th>';
	out = out + '</tr>';
	out = out + '<tr>';
	var b = return_simple_button (tvalue[1],color);
	out +=  '<td style="vertical-align: middle;">' + b + '</td>';
	b = return_simple_button (tvalue[2],color);
	out +=  '<td >' + b + '</td>';
	b = return_simple_button (tvalue[3],color);
	out +=  '<td >' + b + '</td>';
	out = out + '</tr>';
	out = out + '</table>';
	return out;
	
}
function formaterInfoCallerDepth(value) {
	
	
	if (value === "-")
	{
		var out = "manta";
		out =  '<table class= "table table-sm table-striped table-condensed table-bordered table-primary table_gnomad" style="box-shadow: 1px 1px 6px #337AB7 ;font-size: 7px;font-family:  Verdana;margin-bottom:0px;opacity:0.2" >';
		out = out + '<caption style="caption-side: top; text-align: center; font-weight: bold;opacity:0.2">Canvas</caption>';
		out = out + '<tr>';
		out = out + '<th  '+th_style+'>Ratio</th>';
		out = out + '<th '+th_style+'>Score</th>';
		out = out + '</tr>';
		out = out + '<tr>';
		var b = return_simple_button ("-",color);
		out +=  '<td style="vertical-align: middle;">' + b + '</td>';
		b = return_simple_button ("-",color);
		out +=  '<td >' + b + '</td>';
		out = out + '</tr>';
		out = out + '</table>';
	return out;
		return  value;
	}
	var tvalue = value.split(";");
	var name = tvalue[0];
//	var qual = tvalue[1];
//	var text= name
	var color= "#8D71B4";
	var hcolor = "#E8FFF8";
	 color = return_color_by_type(tvalue[3]);
	var th_style = 'style="padding:1px;border: none;text-align: center;"';
	var out = name;
	out =  '<table class= "table table-sm table-striped table-condensed table-bordered table-primary table_gnomad" style="box-shadow: 1px 1px 6px #337AB7 ;font-size: 7px;font-family:  Verdana;margin-bottom:0px" >';
	out = out + '<caption style="caption-side: top; text-align: center; font-weight: bold;">Canvas</caption>';
	out = out + '<tr>';
	out = out + '<th  '+th_style+'>Ratio</th>';
	out = out + '<th '+th_style+'>Score</th>';
	out = out + '</tr>';
	out = out + '<tr>';
	var b = return_simple_button (tvalue[1],color);
	out +=  '<td style="vertical-align: middle;">' + b + '</td>';
	b = return_simple_button (tvalue[2],color);
	out +=  '<td >' + b + '</td>';
	out = out + '</tr>';
	out = out + '</table>';
	return out;
	
}


function formaterBP(value)
{
	
	var out="";
	var chr1;
	var chr2;
	var pos1;
	var pos2;
	
	if ((value == ";") || (value =="-") || (value ==" "))
	{
		out = "-";
	}
	else
	{
		if ( value == "X" )
		{
			out ='<b><font  style="color:grey"> X </font></b>';
		}
		else
		{
			var tval = value.split(";");
			var len = tval.length;
		
			for ( i = 1 ; i < len-1 ; i++)
			{
				tab = tval[i].split("_");
				
				chr1 = tab[0];
				chr2 = tab[2];
				
				if (chr1 == 23)
				{
					chr1 ="X";
				}
				if (chr1 == 24)
				{
					chr1 ="Y";
				}
				if (chr2 == 23)
				{
					chr2 ="X";
				}
				if (chr2 == 24)
				{
					chr2 ="Y";
				}
				
				pos1 = formatNumber(tab[1]);
				pos2 = formatNumber(tab[3]);
				
				var type = '<font style="color:DarkOrange;">t';
				if (chr2 == chr1) {type = '<font style="color:indigo">inv';}
				
				var event_id = '<b>' + type + '(' + chr1 + ',' + chr2 + ')<br></font></b>' + pos1 + '-' + pos2;
				out = out + "<div>" + event_id + "</div>";
			}
		}
	}
	return out;
}

function return_style_value (color,value,font = "9" ){
	return '<button type="button" class="btn btn-xs btn-primary " style="box-shadow: 1px 1px 2px rgba(255, 0, 0, 0.2) ;background-color:'+color+';color:white;font-size: '+font+'px;font-family:  Verdana;border-color:white;">'+value+'</button>';
}

function formaterScoreCNV(value)
{

	if (value <= -2 )
	{
		return "-";
	}

	var eurl;
	var val = Number(value);
	if (val >= 4  ){
		color="tomato";
	}
	else  if (val >= 2){
		color="orange";
	}
	else {
		 color="blue";
	}
	return_style_value(color,val);

	
}	

function formaterScoreGene(value)
{

	if (value <= -2 )
	{
		return '<button "type="button" class="btn btn-xs btn-primary " style="border-color:grey;background-color:white;color:grey;font-size: 12px;font-family:  Verdana">-</button>';
	}
	var val = Number(value);
	if (val >= 4  ){
		color="tomato";
	}
	else  if (val >= 3){
		color="#cd5c5c";
	}
	else  if (val >= 2){
		color="orange";
	}
	else {
		 color="blue";
	}
	return return_style_value(color,val);

}	


function formaterTransmission(value)
{

	
	if (value === "?")
	{
		 return return_style_value('grey','<span class="glyphicon glyphicon-question-sign" aria-hidden="true"></span>');;
		 return return_style_value('','<span class="glyphicon glyphicon-question-sign" aria-hidden="true"></span>');;
	}
	
	if ( value == "X" )
	{
			out ='<b><label  style="color:grey"> X </label></b>';
			return out;
	}
	var img = "?";
	var color = "white";
	
	if ( value === "m"){
		img = '<img src="../../images/polyicons/icons8-person-female-24.png" width="20" height="20">';
	}
	else if ( value === "f"){
		img = '<img src="../../images/polyicons/icons8-person-24.png" width="20" height="20">';
	}
	else if ( value === "fm") {
		img = '<img src="../../images/polyicons/icons8-person-24.png" width="16" height="16">'+'<img src="../../images/polyicons/icons8-person-female-24.png" width="16" height="16">';
	}
	else if ( value === "d"){
		img = 'denovo';
		color = " #FF5733"
	}
	return return_style_value(color,img);
	
//	var TransM="";
//	var TransF="";
//	var valM="";
//	var valF="";
//	var res;
//	
//	var tval = value.split(" ");
//	var len = tval.length;
//		
//	for ( i = 1 ; i < len-1 ; i++)
//	{
//		var val = tval[i];
//		
//		if ( val == "mother")
//		{
//			TransM =  '<img src="../../images/polyicons/icons8-person-female-24.png" >';
//		}
//		if (/maybeM/.test(val))
//		{
//			valM = val.replace("maybeM","mother");
//		}
//		
//		if ( val == "father")
//		{
//			TransF = '<img src="../../images/polyicons/icons8-person-24.png" >';
//		}
//		if (/maybeF/.test(val))
//		{
//			valF = val.replace("maybeF","father");
//		}
//	}
//	
//	if (TransM=="" && valM != "" )
//	{ 
//		TransM = '<b><i><label style="font-size:10px;color:pink">' + valM  + '</label></i></b>';
//	}
//	
//	if (TransF=="" && valF != "" )
//	{ 
//		TransF = '<b><i><label style="font-size:10px;color:lightblue">' + valF  + '</label></i></b>';
//	}
//	
//	if (TransM=="" && TransF == "" )
//	{
//		res = '<b><i><label style="font-size:10px;color:red"> strict-denovo  </label></i></b>';;
//	}
//	else
//	{
//		res = TransM+TransF;
//	}
//	
//	var eurl = '<div>' + res +  '</div>';
//	
//	
//	 return eurl;
}	


function formaterGenotype(value)
{
		var out = "";
		var val;

		if (! /[012]/.test(value))
		{
			return "-";
		}
		else
		{
			var tval = value.split(" ");
			var len = tval.length;
			
			for ( i=0 ; i <len-1 ; i++)
			{
				if ( /[012]/.test(tval[i]) )
				{	

					if (tval[i] == "1/2")
					{
						tval[i] = "0/1";
					}

					if ( tval[i]+" " != out )
					{
						out = out + tval[i] + " ";
					}
					
				}
			}
		}
		
		if (/1\/1/.test(out))
		{
			out = '<font style="color:red">' + out  + '</font>';
		}
		return out;	
}

function formaterRANK(value)
{
	var eurl;
	
		if (value == "-")
		{
			return value;
		}
	
	
	 if (value > 5) {
	 	eurl = '<label style="color:red">' + value + ' </label>';
	 }
	 else {
	 	eurl = '<label style="color:black">'+ value + '</label>';
	 }
	 
	return eurl;
	
}	


function formaterDGV(value)
{
		if (value == "-")
		{
			return value;
		}

	var eurl = '<a href="' + value  + '" target="_blank"><button class="btn btn-default btn-xs" ><b><i><font style="color:red;font-size:10px; ">  DGV  </font></i></b></button></a>';
	return eurl;
}

function formatergnomAD(value)
{
		if (value == "-")
		{
			return value;
		}

		var eurl = '<a href="' + value  + '" target="_blank"><button class="btn btn-default btn-xs" ><b><font style="color:blue;font-size:8px; ">gnom</font><font style="color:black;font-size:10px; ">AD</font></b></button></a>';
	return eurl;
}

function formaterFreq(value)
{
		if (value == "-")
		{
			return value;
		}

	var eurl;

	if (value > 0)
	{
		value = value*100;
		eurl = value.toFixed(3);
		
	 		if (value < 0.05)
			{
	 			eurl = '<font style="color:red">' + eurl + '</font>';
	 		}
	 		else
	 		{
	 			if (value < 0.1)
				{
	 				eurl = '<font style="color:orange">' + eurl + '</font>';
	 			}
	 			else
	 			{
	 				eurl = eurl ;
	 			}
	 		}	 	
	 		return eurl;
	 }
	 return "-";
}	


function formaterdbvarStatus(value)
{
	
	if (value == "-")
	{
			return value;
	}
	
	var eurl=value;
	var out=" ";
	
	var status = value.split(" ");
	var len = status.length;
	
	for ( i = 0 ; i < len; i++)
	{
				if (status[i] != "-")
				{	
							out = out + "  " + status[i] ;
				}
	}
	
	if (( /pathogenic/.test(out) ) || ( /Pathogenic/.test(out) ) )
	{
		out = '<a href="#" data-toggle="tooltip" title="' + out + '"><i class="fa fa-asterisk" style="font-size:14px;color:DarkOrange"> </i></a>';
	}
	else 
	{
		if (( /benign/.test(out) ) || ( /Benign/.test(out) ) )
		{
			out = '<a href="#" data-toggle="tooltip" title="' + out + '"><i class="fa fa-asterisk" style="font-size:14px;color:green"> </i></a>';
		}
		else
		{
			out = "-";
		}
	}
	return out;
}		



//--------------------------------------------------------------------
// pour les Translocations et les inversions
//---------------------------------------------------------------------

var dataStore_TRANS;
var url_BND_Project = url_path + "/manta/PolyCytoRetrieveBND.pl";

function setTitreBND()
{

	var parameters = location.search.split("&");
	var par0 = parameters[0].split("=");
	var par1 = parameters[1].split("=");
   	
   	var project =par0[1];
   	var patient = par1[1];
   	
	titre = '<br>Translocations BreakPoint Viewer  <br> ' + project + '   /   ' + patient ;
	document.getElementById('titre').innerHTML=titre;
}

var h_gridTRANS = {};
function GetSVeq() 
{
	document.getElementById("span_transloc_inv").innerHTML = "<img src='../../images/polyicons/wait18trans.gif'> TRANSLOC/INV";
	
		var parameters = location.search.split("&");
		var arg1 = parameters[0].split("=");
   		var arg2 = parameters[1].split("=");
		
    	var project = arg1[1];
    	var patient = arg2[1];
    	
		var chr = document.getElementById('chromSelect').value;
		var dejavu = document.getElementById('dejavu').value;
		var transmission = document.getElementById('transmissionSelect').value;
		var genes = document.getElementById('geneSelect').value;
		var score_event=document.getElementById('scoreEvent').value;
		var omim=0;
		
		var qual = 1;
		var type = "all";
		
		var this_grid = dijit.byId("TRANS_gridDetails");
		var emptyCells = { items: "" };
		var emptyStore = new dojo.data.ItemFileWriteStore({data: emptyCells});
		this_grid.setStore(emptyStore);
		this_grid.showMessage("Chargement...");
	

		if (genes == "") { genes = "all"};
		
		if (document.getElementById('OMIMSelect').checked == true)
		{
			omim=1;
		}
		
		var url = url_BND_Project;
		//var args = "project=" + project + "&patient=" + patient + "&dejavu=" + dejavu + "&chrom=" + chr + "&qual=" + qual;
		var args = "project=" + project + "&patient=" + patient + "&dejavu=" + dejavu + "&chrom=" + chr + "&transmission=" + transmission + "&genes=" + genes + "&omim=" + omim + "&score_event=" + score_event;
		var xhr = new XMLHttpRequest();
		xhr.responseType = 'json';
		xhr.open('POST', url, true);
		xhr.setRequestHeader('Content-type', 'application/x-www-form-urlencoded');
		xhr.onreadystatechange = function () {
			var json_reponse = xhr.response;
			dataStore_TRANS = new dojo.data.ItemFileReadStore({ data: json_reponse });
			this_grid.setStore(dataStore_TRANS);
			this_grid.store.fetch({ onComplete: function(items){
				document.getElementById("span_transloc_inv").innerHTML = 'TRANSLOC/INV';
				h_gridTRANS[patient] = this_grid;
			} });
		}
		xhr.send(args);
	return;
}

var nb_check_force_refresh = 0;
function force_refresh_grids() {
	if (nb_check_force_refresh >= 1) { return; }
	nb_check_force_refresh++;
	$('c1polycyto'+name).prevObject[0].links[1].click();
	$('c1polycyto'+name).prevObject[0].links[0].click();
	refresh_grids();
}

function refresh_grids() {
	var parameters = location.search.split("&");
	var arg1 = parameters[0].split("=");
	var arg2 = parameters[1].split("=");
	var projectname = arg1[1];
	var filename = arg2[1];
	$('a[data-toggle="tab"]').on('shown.bs.tab', function (e) {
		var target = $(e.target).attr("href") // activated tab
		if (filename in h_gridSVCompare) {
			h_gridSVCompare[filename].render();
		}
		if (filename in h_gridTRANS) {
			h_gridTRANS[filename].render();
		}
	});
	return;
}
//
function formaterPosBND(value)
{
	//alert('posbnd'+value);
	if (value=="-")
	{
		return(value);
	}
	return formatNumber(value);
}	



function formaterCHROM_BND1(value)
{
	//alert('chr1 = '+value);
	if (value=="-")
	{
		return(value);
	}
	
	var eurl = value;
	
	if (value == 23)
	{
		 eurl ="X";
	}
	if (value == 24)
	{
		 eurl="Y";
	}

	var out = '<center> <table><tr><td style="border: 2px solid indigo;width:40px; vertical-align:middle;font-size: 10px;text-align:center;">' + eurl + '</td></tr></table></center>';
	return out;
}

function formaterCHROM_BND2(value)
{
	//alert('chr2 = '+value);
	if (value=="-")
	{
		return(value);
	}
	
	var eurl = value;
	
	if (value == 23)
	{
		 eurl ="X";
	}
	if (value == 24)
	{
		 eurl="Y";
	}

	var out = '<center> <table><tr><td style="border: 2px solid DarkOrange;width:40px; vertical-align:middle;font-size: 10px;text-align:center;">' + eurl + '</td></tr></table></center>';
	return out;
}

function formaterGENES_SV(value) {
	
	let obj = JSON.parse(value);
	if  (obj.genes.length == 0) {return "-";}
	//unless  (obj.genes) {return "-";}
		
	var tid = obj.table_id;
	var hcolor = "#E8FFF8";
var th_style = 'style="padding:1px;border: none;text-align: center;"';
th_style = '';
	var out ="";
	out = out + '<table id="'+tid+'" class= "table table-sm table-striped table-condensed table-bordered table-primary table_gnomad" style="box-shadow: 1px 1px 6px #337AB7 ;font-size: 7px;font-family:  Verdana;margin-bottom:0px" >';
	out = out + '<tr>';

		var data =obj.genes.slice(0, 4); 
    
    	  data.forEach(gene => {
			  var color = return_color_by_score_gene(gene.score);
			  var url = 'https:\/\/www.ensembl.org\/Homo_sapiens\/Gene\/Summary\?db=core;g='+gene.name;
			  var b_djv_itp=  '<button type="button" class="btn btn-xs btn-primary " style="box-shadow: 1px 1px 2px rgba(0, 0, 0, 0.2);background-color:'+color+';color:white;font-size: 8px;font-family:  Verdana;border-color:white;" onclick="window.open(\'' + url + '\', \'_blank\')">'+gene.name+'</button>';
				out = out +'<td>'+ b_djv_itp+'</td>'
			});
		out = out + '</tr>';


	var stid="";
	out = out + '</tr>';

	out = out + '<tr>';
	

	out = out + '<td colspan=4>';
	var onclick = 'onclick=get_all_genes(\"'+obj.project+'\",\"'+obj.id+'\",\"'+obj.patient+'\")';
	var tname = "all";
	
	//var text = "+ high("+obj.nbh+") medium("+obj.nbm+") low("+obj.nbl+")";
	text =  obj.nb_genes;
	                                                                                  
	var ball=  '<button type="button" class="btn btn-xs btn-primary " style="  box-shadow: 1px 1px 4px rgba(0, 0, 255, 0.3);background-color:#0468BF;color:white;font-size: 8px;font-family:  Verdana;border-color:white;" '+onclick+'> View  details  ('+text+' genes)</button>';
	out = out + ball+'</td></tr>';
	
	out = out + '</table>';
	return out;
	
	
}



function get_all_genes(this_proj_name,id,patient) {
	
	var this_url = url_path + "/manta/getGenesForEvent.pl";
	var this_args;
	this_args = "project=" + this_proj_name + "&id=" + id+"&patient="+patient; 
		dijit.byId("genes_popup").show();
	 document.getElementById('genes_div').innerHTML = "Loading ...";
	var xhr = new XMLHttpRequest();
	xhr.open('POST', this_url, true);
	xhr.setRequestHeader('Content-type', 'application/x-www-form-urlencoded');
	xhr.onreadystatechange = function() {
		if (this.readyState === XMLHttpRequest.DONE && this.status === 200) {
			var data = JSON.parse(this.responseText);
			var out ="";
			var out ="";
	out = out + '<table class= "table table-sm table-striped table-condensed table-bordered table-primary table_gnomad" style="box-shadow: 1px 1px 6px #337AB7 ;font-size: 7px;font-family:  Verdana;margin-bottom:0px" >';
	out = out + '<tr>';

			var project = data.project;
			 data.genes.forEach(gene => {
			  var color = return_color_by_score_gene(gene.score);
			  var name = gene.name+"<sup>"+gene.score+"</sup>";
			  var url = 'https:\/\/www.genecards.org\/cgi-bin\/carddisp.pl?gene='+gene.name;
			  var url2 = 'https:\/\/www.ensembl.org\/Homo_sapiens\/Gene\/Summary\?db=core;g='+gene.name;
			  var url3 = 'https:\/\/gnomad.broadinstitute.org/gene/'+gene.name+'?dataset=gnomad_r4';
			  
			  https://gnomad.broadinstitute.org/gene/ENSG00000171634?dataset=gnomad_r4
			  var b_djv_itp=  '<button type="button" class="btn btn-xs btn-primary " style="box-shadow: 1px 1px 2px rgba(0, 0, 0, 0.2);background-color:'+color+';color:white;font-size: 12px;font-family:  Verdana;border-color:white;" onclick="window.open(\'' + url2 + '\', \'_blank\')">'+name+'</button>';
			   var b_ensembl=  '<button type="button" class="btn btn-xs btn-primary " style="box-shadow: 1px 1px 2px rgba(0, 0, 0, 0.2);background-color:'+color+';color:white;font-size: 8px;font-family:  Verdana;border-color:white;" onclick="window.open(\'' + url2 + '\', \'_blank\')">Ensembl</button>';
		 	 var b_card=  '<button type="button" class="btn btn-xs btn-primary " style="box-shadow: 1px 1px 2px rgba(0, 0, 0, 0.2);background-color:'+color+';color:white;font-size: 8px;font-family:  Verdana;border-color:white;" onclick="window.open(\'' + url + '\', \'_blank\')">gene cards</button>';
			var b_gnomad=  '<button type="button" class="btn btn-xs btn-primary " style="box-shadow: 1px 1px 2px rgba(0, 0, 0, 0.2);background-color:'+color+';color:white;font-size: 8px;font-family:  Verdana;border-color:white;" onclick="window.open(\'' + url3 + '\', \'_blank\')">Gnomad</button>';
		
			  
			  var button_view =  '<button type="button" class="btn btn-xs btn-primary " style="box-shadow: 1px 1px 2px rgba(0, 0, 0, 0.2);background-color:'+color+';color:white;font-size: 10px;font-family:  Verdana;border-color:white;\" ';
			
			  var phenotypes = '<ul style="list-style-type: disc; padding-left: 20px; font-family: Verdana; font-size: 14px; color: #333;">';
			  var ap  = gene.phenotypes.slice(0, 3); 
			  for (let i = 0; i < ap.length; i ++) {
				  phenotypes+='<li>'+gene.phenotypes[i]+"</li>";
				 }
				  button_view += 'onclick=get_gene_phenotypes2(\"'+project+'\",\"'+gene.id+ '\")>View All Phenotypes ('+gene.phenotypes.length+')</button>';
				  
				  phenotypes+='<li>'+button_view+'</li>';
			   phenotypes += '</ul>';
			 out = out +'<td>'+'<span style="white-space: nowrap;">'+b_ensembl+b_card+b_gnomad+'</span><br><center>'+b_djv_itp+'</center></td><td><span style="font-size: 10px;font-family:  Verdana;color:black">'+phenotypes+"</span></td>";
			out = out + '</tr>';
			});
				out = out + '</table>';

	document.getElementById('genes_div').innerHTML = out;
			}
		}
	
	xhr.send(this_args);
}

function return_color_by_score_gene(value){
	if(value >= 5){
		return "#CC202A";
	}
	else if(value >= 4){
		return "#E2542D";
	}
	else if(value >= 3){
		return "#D1805E";
	}
	else if(value >= 2){
		return "#E3BC33";
	}
	else {return "#BDBEBF"}
}
function formaterGENES(value)
{
	//alert("genes =" + value);
	var parameters = location.search.split("&");
	var par0 = parameters[0].split("=");
	var par1 = parameters[1].split("=");
   	
   	var project =par0[1];
   	var patient = par1[1];
   	
  
	if (value == "")
	{
		return "-";
	}
	
	if (value=="-")
	{
		return(value);
	}
	
	var out ="";
	var name = "";
	var score;
	var locus;
	var phen = "";
	
	var outshort;
	var values = value.split("##");
	var len = values.length-1;
	
	var liste1="";
	var liste2="";
	var liste3="";
	var nbghigh=0;
	var nbglow=0;
	var nbplus=0;
	
	locus = values[0];
	
	for ( i = 1 ; i < len; i++)
	{
		if (values[i] != "-")
		{	
					var geneinfos = values[i].split(";");
					name = geneinfos[0];
					score = geneinfos[1];
										
					if (score < 3 )
					{
						if(nbglow < 3)
						{
							if (score >= 2)
							{
								liste1 = liste1 + '<font style="font-size:10px;color:white;background-color:rgb(255,165,0);padding:1px">' + name +'</font><br>';
							}
							else
							{
								if (score >= 1)
								{
									liste1 = liste1 + '<font style="font-size:10px;color:black;background-color:rgb(255,255,2);padding:1px">' + name +'</font><br>';
								}
								else
								{
									liste1 = liste1 + '<font style="font-size:10px;color:white;background-color:rgb(128,128,128);padding:1px">' + name +'</font><br>';
								}
							}
							nbglow++;
						}
						else
						{
							nbplus++;
						}
					}
					else
					{
						if(nbghigh<10)
						{
							if (score >= 4)
							{
								liste2 = liste2 + '<font style="font-size:10px;color:white;background-color:rgb(255,0,0);padding:2px;">'+name+'</font><br>';
							}
							else
							{
								liste2 = liste2 + '<font style="font-size:10px;color:white;background-color:rgb(255,127,79);padding:2px">'+name+'</font><br>';							
							}
						}
						else
						{
							nbplus++;
						}
						nbghigh++;
						var flag = 1;
					}
		}
	}
	
	if (flag ==1)
	{
		if (nbplus > 0)
		{
			out = '<button class="btn btn-classic btn-xs"  style="border: 1px solid black;padding:2px" onclick="viewGeneByName(\'' + locus + '\', \''+ project + '\',\'' + patient + '\')">' + liste2 + '<font style="font-size:9px;"> (+' + nbplus + ')</button>';
		}
		else
		{
			out = '<button class="btn btn-classic btn-xs"   style="border: 1px solid black;padding:2px" onclick="viewGeneByName(\'' + locus + '\', \''+ project + '\',\'' + patient + '\')">' + liste2 + '</button>';
		}
		
	}
	else
	{
		if (nbplus > 0)
		{
			out = '<button class="btn btn-classic btn-xs"  style="border: 1px solid black;padding:2px" onclick="viewGeneByName(\'' + locus + '\', \''+ project + '\',\'' + patient + '\')">' + liste1  + '<font style="font-size:9px;"> (+' + nbplus + ')</button>';
		}
		else
		{
			out = '<button class="btn btn-classic btn-xs" style="border: 1px solid black;padding:2px"  onclick="viewGeneByName(\'' + locus + '\', \''+ project + '\',\'' + patient + '\')">' + liste1  + '</button>';
		}
	}
	return out;
}	


function formaterOMIM(value)
{
	var out = ""; 
	
	if (value == "-")
	{
		return("-");
	}
	
	if (value == "0")
	{
		return("-");
	}
	
	if (value == 1)
	{
		out =  '<i class="fa fa-asterisk" style="font-size:14px;color:DarkOrange"> </i>';
	}
	return out;
}



function formaterFreqBND(value)
{	
	//alert('freqbnd'+value);
	if (value=="-")
	{
		return(value);
	}
	
	out = value.toFixed(2);
	return out;
}


function formaterTransloc(value)
{
	//alert('transloc'+value);
	
	if (value=="-")
	{
		return(value);
	}
	var infos = value.split(";");
	
	var color = return_color_by_type(infos[0]);
	 	var text = infos[1];
	 	eurl=  '<button type="button" class="btn btn-xs btn-primary " style="border-color:black;padding:5px;background-color:'+color+';color:white;font-size: 14px;font-family:  Verdana;border-color:white;">'+text+'</button>';
	 	//eurl = '<p style="font-size:10px;color:white;background-color:rgb(255,74,42);padding:8px">' + value + ' XXX</p>';
	 return return_style_value(color,text,11);
	return eurl;
	
}

function formaterInv(value)
{
	if (value=="-")
	{
		return(value);
	}
	

	var chrs = value.split("to");
	
	if (chrs[0] == 23)
	{
		 chrs[0] ="X";
	}
	if (chrs[0] == 24)
	{
		 chrs[0]="Y";
	}
	
	
	var out = '<center><table><tr><td style="border: 2px solid indigo;width:30px; vertical-align:middle;text-align:center;font-size: 10px;">' + chrs[0] + '</td></tr></table></center>';
	return out;
}


function formaterTypeInv(value)
{
	if (value=="-")
	{
		return(value);
	}
	
	var eurl;
	
	if (value == "Paracentrique")
	{
		eurl = '<font style="color:black"> Paracent. </font>';
	}
	else
	{
		eurl = '<font style="color:red"> Pericent. </font>';
	}
	return eurl;
}

 function formaterREFALT(value) 
{
	//alert('refalt'+value);
	if (value=="-")
	{
		return(value);
	}
	
	var vals = value.split("/");
	
	out = '<a href="#" data-toggle="tooltip" title="' + vals[1] + '"> ' + vals[0] + ' / <u><b> Alt </b></u></a>';
	return out;
	
 }
 

function formater_BND_IGV(value)
{
		if (value=="-")
		{
			return(value);
		}
	
		var parameters = location.search.split("&");
		var arg1 = parameters[0].split("=");
  	 	var arg2 = parameters[1].split("=");
		
  		 var projectname = arg1[1];
  	  	var patientname = arg2[1];
    	
    
		var out = " ";
		
		var tv =  value.split('et');
		var baminfo = tv[0];
		var tbam = baminfo.split(';');
		
		var bamfiles = tbam[0];
		var bamnames = patientname;
		
		if(tbam.length ==2)
		{
			bamnames = tbam[1];
		}
		
		var real_value = tv[1];
	
		var w =  real_value.split("_");
		var chr1 = w[0];
		var chr2 = w[2];
		
		if (chr1 == 23) { chr1 = "X"; }
		if (chr2 == 23) { chr2 = "X"; }
		if (chr1 == 24) { chr1 = "Y"; }
		if (chr2 == 24) { chr2 = "Y"; }
		 
			
		// 1  kb autour du start
		var s1 = parseInt(w[1]);
		
		var start1 = s1 - 500;
		var end1   = s1 + 500;
		var locus1 = 'chr' + chr1  + ':' + start1 + '-' + end1;
				
		// 1  kb autour du end
		var s2 = parseInt(w[3]);
		var start2 = s2 - 500;
		var end2   = s2 + 500;
		var locus2 = 'chr' + chr2  + ':' + start2 + '-' + end2;
   			
   		
   		nb++;
   		
   		var idCNV= patientname + '_t(' + chr1 + ':' + chr2 + ')';
   		
   		// cas particulier des inversions
   		if(chr1==chr2)
   		{
   			idCNV= patientname + '_inv(' + chr1 + ':' + chr2 + ')';
   		}
   		
   		var boutonID = 'boutonIGVTI'+ nb;					
   		var igv_value = locus1 + ';' + locus2 + ';' + bamfiles + ";" +bamnames;	
   		
   		// pour lancer web_igv
		var etiquette   = '<button class="igvIcon2"  id="' + boutonID + '" onclick="launch_web_igv_start_and_end(\'' + idCNV + "et" + igv_value + '\')"></button>';
		
		
		return etiquette;
}

function formaterTransmissionBND(value)
{
	if ((value == "-") || (value == "X"))
	{
		return value;
	}
	
	var eurl;
	
	if (/strictdenovo/.test(value) )
	 {
	 	eurl = '<b><i><label style="font-size:10px;color:red"> strict denovo  </label></i></b>';
	 }
	
	if (/denovo/.test(value) )
	 {
	 		eurl = '<b><i><label style="font-size:10px;color:orange"> denovo  </label></i></b>';
	}
	
	 
	if ( /mother/.test(value) & !/father/.test(value) )
	{
	 	eurl = '<div class=circle2F><img src="../../images/polyicons/icons8-person-female-24.png" > </div>';
	 }
	 
	 if  ( /father/.test(value) & !/mother/.test(value) )
	 {
	 	eurl = '<div><img src="../../images/polyicons/icons8-person-24.png" ></div>';
	 }
	 
	 if  ( /both/.test(value) | ( /father/.test(value) & /mother/.test(value) ) ) {
	 	eurl = '<div><img src="../../images/polyicons/icons8-person-female-24.png" ><img src="../../images/polyicons/icons8-person-24.png"></div>';
	 }
	 
	 return eurl;
}	


function formaterCNVlinked(value)
{
	var out="";
	var chr1;
	var chr2;
	var pos1;
	var pos2;
	
	if ((value == ";") || (value =="-") || (value ==" "))
	{
		out = "-";
	}
	else
	{
		if ( value == "X" )
		{
			out ='<b><font  style="color:grey"> X </font></b>';
		}
		else
		{		
			var tab = value.split(";");
			var len = tab.length;
	
			for (i=0;i<len;i++)
			{
		 		if ( /DUP/.test(tab[i]) )
		 		{
					out = out + '<font style="color:blue;font-size:8px;"> ' + tab[i] + '</font>';	 
				}
				else
				{
					out =out +  '<font style="color:red;font-size:8px"> ' + tab[i]  + '</font>';	 
				}
			}
		}
	}
	return out;
}

//--------------------------------------------------------------------
// pour les expansions de triplets
//---------------------------------------------------------------------



var url_TRIPLETS  = url_path + "/manta/parseRepeatVcf.pl";


var h_gridTRIPLETS = {};
var dataStore_TRIPLETS;

function GetTriplets_Get() 
{
	document.getElementById("span_exp_triplet").innerHTML = "<img src='../../images/polyicons/wait18trans.gif'> TRIPLETS EXPANSION";

		var parameters = location.search.split("&");
		var arg1 = parameters[0].split("=");
   		var arg2 = parameters[1].split("=");
   		
		
    	var project = arg1[1];
    	var patient = arg2[1];
		var dataStore_TRIPLETS = new dojo.data.ItemFileReadStore({ url: url_TRIPLETS+ "?project=" + project + "&patient=" + patient});
		
		
		h_gridTRIPLETS[patient] = dijit.byId("TRIPLETS_gridDetails");
		h_gridTRIPLETS[patient].setStore(dataStore_TRIPLETS);
		h_gridTRIPLETS[patient].store.fetch({
			onComplete: function(items){
				document.getElementById("span_exp_triplet").innerHTML = "TRIPLETS EXPANSION";
			}
		});
	
		return;
}


function GetTriplets_Post() 
{
	document.getElementById("span_exp_triplet").innerHTML = "<img src='../../images/polyicons/wait18trans.gif'> TRIPLETS EXPANSION";
	
		var parameters = location.search.split("&");
		var arg1 = parameters[0].split("=");
   		var arg2 = parameters[1].split("=");
		
    	var project = arg1[1];
    	var patient = arg2[1];
    	
    
 
		var this_grid = dijit.byId("TRIPLETS_gridDetails");
		var emptyCells = { items: "" };
		var emptyStore = new dojo.data.ItemFileWriteStore({data: emptyCells});
		this_grid.setStore(emptyStore);
		this_grid.showMessage("Chargement...");
	
		
		var url = url_TRIPLETS;
		var args = "project=" + project + "&patient=" + patient;
		
		
		var xhr = new XMLHttpRequest();
		xhr.responseType = 'json';
		xhr.open('POST', url, true);
		xhr.setRequestHeader('Content-type', 'application/x-www-form-urlencoded');
		
		xhr.onreadystatechange = function () {
			var json_reponse = xhr.response;
			dataStore_TRIPLETS = new dojo.data.ItemFileReadStore({ data: json_reponse });
			this_grid.setStore(dataStore_TRIPLETS);
			this_grid.store.fetch({ onComplete: function(items){
				document.getElementById("span_exp_triplet").innerHTML = 'TRIPLETS EXPANSION';
				h_gridTRIPLETS[patient] = this_grid;
			} });
		}
		xhr.send(args);
		
	return;
}


function formaterTriplet(value)
{
	alert(value);
	return value;

}

function formaterGeneTriplet(value)
{
	alert('gene = '+value);
	
	return value;
}

//------------------------
// Pour IGV
//------------------------


function launch_web_igv_start_and_end(value) 
{
	var v = value.split(";");
	var first_col =  v[0].split("et");
	var titre = first_col[0];
	var locusstart = first_col[1];
	var locusend = v[1];
	var bamfile =v[2];
	var patient_name = v[3];
	var taille = 1;
	if (v.length > 4) 
	{
		taille = v[4];
	} 
	
 	// pour lancer appli IGV si elle est ouverte
	var tab_bamfile = bamfile.split(',');
	displayListBamInIgvApp(tab_bamfile); 
		
	var locus =locusstart+ ' ' +locusend;
    displayInIGVLocus(locus);
	
	// pour lancer web_igv
	var url_viewer = "PolyCyto_IGV_Viewer.html?titre=" + titre + "&locusstart=" + locusstart + "&locusend=" + locusend + "&bamfile=" + bamfile + "&patient_name=" + patient_name + "&taille=" + taille;
	url_viewer = url_viewer.replace('bamfile=https://www.polyweb.fr','bamfile=');
	if (url_gene_bed) { url_viewer += "&genes_bed=" + url_gene_bed; }
	if (url_fasta) { url_viewer += "&fasta=" + url_fasta; }
	var myWindow = window.open(url_viewer,"_blank",""); 
 
	
 }
 
 
 function igv_view_SV()
{
   	var titre = param("titre");
   	var locusstart = param("locusstart");
   	var locusend = param("locusend");
	var bamfile = param("bamfile");
	var patient_name = param("patient_name");
	var taille = param("taille");
	var genes_bed = param("genes_bed");
	var fasta = param("fasta");
	
	if (taille != 1) {
		titre = patient_name + "  :  " + titre + "  (" + taille + ")";
	}
	
	if (genes_bed) { url_gene_bed = genes_bed; }
	if (fasta) { url_fasta = fasta; }
	document.getElementById('titre').innerHTML=titre;
	var locus = locusstart + " " + locusend;
	view_web_igv_bam('','div_igv_start',locus, bamfile, patient_name);
}

//-------------------------
// Pour le dejavu
//------------------------

function my_view_dejavu_polycyto(value) {
	dijit.byId('waiting').show();
	var values = value.split("+");
	var out = '<span id="span_phenotype_name" ></span>';
	var titre = "&nbsp;&nbsp;DejaVu  : "+ values[0]+"&nbsp;&nbsp;&nbsp";
	document.getElementById('content_dejavu_polycyto').innerHTML=out;
	document.getElementById('dialog_dejavu_polycyto_title').innerHTML=titre;
	view_dejavu_polycyto(values[1]);
}


function formaterDejavuCNV(value)
{
	var dejavu_seuil;
	
	if (value == "-")
	{
			return value;
	}
	
	if ( document.getElementById('dejavu').value != "all")
	{
		dejavu_seuil = parseInt(document.getElementById('dejavu').value);
	}
	else
	{
		dejavu_seuil = 1000000;
	}	
	
	var args = value.split(";");
	var event_id=args[0];
	
	var values = args[1];
	
	var liste_djv_itp="";				// liste des echantillons du projet
	
//	var liste = ""	
//	var tab = liste.split(",");
//	var len = tab.length-1;
//	
//	for (i=0;i<len;i++)
//	{
//		liste_djv_itp = liste_djv_itp + tab[i]+' ';
//	}
//	
	var djv_in_nbother_project = args[1];
	var djv_in_nbother_patient = args[2];
	var liste_djv_iop ="";
	var nbwisecondor=args[3];
	var nbcanvas=args[4];
	var nbmanta=args[5];
	var flagwc = 0;
	var djv_itp = args[6];
	var projectname = args[7];
	
	var out = "";;
   	var boutonID = 'boutondjv' + nb;	
   	
 
   	// pour l'affichage
   var couleur ="white";
   var bgcolor = "#E6E6E6";
   
/*      if (djv_in_nbother_patient < 20  && djv_in_nbother_patient >= 10){
	  bgcolor="#faf0e6;color:black;";
   }
     if (djv_in_nbother_patient < 10  && djv_in_nbother_patient >= 5 ){
	  bgcolor="#7de3b4;color:black;";
   }

   if (djv_in_nbother_patient < 5 ){
	  bgcolor="#52c891;color:black;";
   }*/
   
  
   var boutonType= "btn btn-classic btn-xs";
   
   
   
//   if ( (flagwc == 1) && (nbwisecondor < dejavu_seuil) )
//   {
//   		bgcolor="#FFBF00;color:black;";
//   
//  		 if ( nbwisecondor > (dejavu_seuil/2) )
//   		{
//   			bgcolor="#F5DA81;color:black;";
//  		 }
//  	}
//  
var hcolor = "#E8FFF8";
var th_style = 'style="padding:1px;border: none;text-align: center;"';
th_style = '';
	//<button "type="button" class="btn btn-xs btn-primary " style="border-color:#e74c3c;background-color:white;color:#e74c3c;font-size: 8px;font-family:  Verdana">0.999</button>
	var b2 = '<button "type="button" class="padding:0px;btn btn-xs btn-primary " style="border-color:black;background-color:white;color:#e74c3c;font-size: 10px;font-family:  Verdana">'+djv_in_nbother_patient+'</button>';
   
	//out =   '<center> <table  width="100%"><tr style="font-size:9px;">';
	//out = out + '<td><center><button style="width:95%; background-color:white;"  id="' + boutonID + '" onclick="my_view_dejavu_polycyto(\'' + value + '\')">';
	
	out = out + '<table onClick="alert(\"toto\")" class= "table table-sm table-striped table-condensed table-bordered table-primary table_gnomad" style="box-shadow: 1px 1px 6px #337AB7 ;font-size: 7px;font-family:  Verdana;margin-bottom:0px" >';
	out = out + '<tr>';
	out = out + '<th  '+th_style+'>Prj</td>';
	out = out + '<th  '+th_style+'>Sam</td>';
	out = out + '<th '+th_style+'>Cov</td>';
	out = out + '<th  '+th_style+'>Dep</td>';
	out = out + '<th '+th_style+'>sr</td>';
	out = out + '<th '+th_style+'>itp</td>';
	out = out + '</tr>';
	out = out + '<tr>';
	
	var cmd = "onclick=\"window.open('https://www.polyweb.fr/HG38/polyweb/polydejavu/dejavu_cnv_sv.html?project=" + projectname + "&id=" + event_id + "', '_blank')\"";
	
	var b = return_button (djv_in_nbother_project, cmd);
	out +=  '<td style="vertical-align: middle;">' + b + '</td>';
	b = return_button (djv_in_nbother_patient, cmd);
	out +=  '<td >' + b + '</td>';
	b = return_button (nbwisecondor, cmd);
	out +=  '<td >' + b + '</td>';
	b = return_button (nbcanvas, cmd);
	out = out + '<td>' + b +'</td>';
	b = return_button (nbmanta, cmd);
	out = out + '<td>' + b + '</td>';			
	//out = out + '</tr></table></center></td>'; 
	var b_djv_itp=  '<button type="button" class="btn btn-xs btn-primary " style="background-color:#8D71B4;color:white;font-size: 8px;font-family:  Verdana;border-color:white;">'+djv_itp+'</button>';
	out = out + '<td ><a href="#" data-toggle="tooltip" title=" ' + liste_djv_itp + ' ">' +  b_djv_itp + ' </a></td>';
	out = out + '</tr>';
	out = out + '</table>';
	return out;
}

function return_simple_button(value,color) {
	return '<button type="button" class="btn btn-xs btn-primary " style="box-shadow: 1px 1px 2px rgba(0, 0, 0, 0.1);background-color:'+color+';color:white;font-size: 8px;font-family:  Verdana;border-color:white;">'+value+'</button>';
	}
function return_button(value, cmd) {
	var color;
	if (value > 20) {
		 color="tomato";
	}
	else if (value > 10) {
		color="#E2552D";
	}
	else if (value > 1) {
		color="#6F8D6A";
	}
	else {
		color="#A6BE47";
	}
	
	return '<button '+ cmd + ' type="button" class="btn btn-xs btn-primary " style="box-shadow: 1px 1px 2px rgba(0, 0, 0, 0.1);background-color:'+color+';color:white;font-size: 8px;font-family:  Verdana;border-color:white;">'+value+'</button>';
	
}

var myAlertDialog;
function warning(text){
	if (!myAlertDialog) {
			 myAlertDialog = new dijit.Dialog({
        			title: "Warning",
        			content: text,
        			style: "width: 200px;height:150px;font-size: 13px;"
    				});
    			}
  
    	myAlertDialog.show();
    				return;
	}




function show_lines (value){
	let arr = value.split('!');
	dijit.byId("genes_popup").show();
	 document.getElementById('genes_div').innerHTML = "Loading ...";
	var out ="";
	out = out + '<table class= "table table-sm table-striped table-condensed table-bordered table-primary table_gnomad" style="box-shadow: 1px 1px 6px #337AB7 ;font-size: 7px;font-family:  Verdana;margin-bottom:0px" >';
	out = out + '<tr>';
	 for (let i = 0; i < arr.length; i ++) {
    	 var data =arr.slice(i, i + 1); 
    	  out = out + '<tr>';
    	  data.forEach(gene => {
			  let ds = gene.split('$');
			  var color = return_color_by_score_gene(ds[2]);
			  var name = ds[0];
			  var id = ds[1];
			  var pheno = ds[3];
			  var url = 'https:\/\/www.ensembl.org\/Homo_sapiens\/Gene\/Summary\?db=core;g='+ds[0];
			  var b_djv_itp=  '<button type="button" class="btn btn-xs btn-primary " style="background-color:'+color+';color:white;font-size: 10px;font-family:  Verdana;border-color:white;" onclick="window.open(\'' + url + '\', \'_blank\')">'+name+'</button>';
			  
			 out = out +'<td>'+ b_djv_itp+'</td><td><span style="font-size: 10px;font-family:  Verdana;color:black">'+pheno+"</span></td>";
			});
		out = out + '</tr>';
  	}
  		out = out + '</table>';

	document.getElementById('genes_div').innerHTML = out;
	
	return;
	
}
function formaterDejavuSVeq(value)
{
	//console.log(value);
	
	let obj = JSON.parse(value);
	var dejavu_seuil;

	var djv_in_nbother_project = obj.nb_projects;
	var djv_in_nbother_patient = obj.nb_patients;
	var djv_in_this_project = obj.nb_itp+"/"+obj.total_itp;
	
	var event_id = obj.id;
	var projectname = obj.project;

var hcolor = "#E8FFF8";
var th_style = 'style="padding:1px;border: none;text-align: center;"';
th_style = '';
var out ='';
	//<button "type="button" class="btn btn-xs btn-primary " style="border-color:#e74c3c;background-color:white;color:#e74c3c;font-size: 8px;font-family:  Verdana">0.999</button>
	out = out + '<table class= "table table-sm table-striped table-condensed table-bordered table-primary table_gnomad" style="box-shadow: 1px 1px 6px #337AB7 ;font-size: 7px;font-family:  Verdana;margin-bottom:0px" >';
	out = out + '<tr>';
	out = out + '<th  '+th_style+'>Prj</td>';
	out = out + '<th  '+th_style+'>Sam</td>';
	out = out + '<th '+th_style+'>itp</td>';
	out = out + '</tr>';
	out = out + '<tr>';
	
	var cmd = "onclick=\"window.open('https://www.polyweb.fr/HG38/polyweb/polydejavu/dejavu_cnv_sv.html?project=" + projectname + "&id=" + event_id + "', '_blank')\"";
	
	var b = return_button (djv_in_nbother_project, cmd);
	out +=  '<td style="vertical-align: middle;">' + b + '</td>';
	b = return_button (djv_in_nbother_patient, cmd);
	out +=  '<td >' + b + '</td>';

	var b_djv_itp=  '<button type="button" class="btn btn-xs btn-primary " style="background-color:#8D71B4;color:white;font-size: 8px;font-family:  Verdana;border-color:white;">'+djv_in_this_project+'</button>';
	out = out + '<td ><a href="#" data-toggle="tooltip" title=" ' + djv_in_this_project + ' ">' +  b_djv_itp + ' </a></td>';
	out = out + '</tr>';
	out = out + '</table>';
	return out;
}


function formaterGeneDejavu(value)
{
			var dataStoreGeneDJV;
			var url_GeneDJV = url_path + "/manta/PolyCytoGetGeneDJVinCNV.pl";
			
			var parameters = value.split(";");
			var project = parameters[0];
			var chr = parameters[1];
			var start = parameters[2];
			var end = parameters[3];
			var type = parameters[4];
		
  			var out;
  			var nbproject;
  			var nbpatient;
  			var CnvGeneDejavu;
  			
  			var boutonType= "btn btn-primary btn-xs";
  			var boutonID = 'boutonGenedjv' + nb;
  			

			//dataStoreGeneDJV = new dojo.data.ItemFileReadStore({ url: url_GeneDJV + "?project=" + project + "&gene=" + gene + "&chr=" + chr + "&start=" + start  + "&end=" + end + "&type=" + type });
			dataStoreGeneDJV = new dojo.data.ItemFileReadStore({ url: url_GeneDJV +"?project=NGS2018_2300&chr=3&start=1134342&end=1445901&type=DUP" });
	
			dataStoreGeneDJV.fetch({
				onComplete: function(items,request){
				
						var item = items[0];
   						nbproject = dataStoreGeneDJV.getValue(item, "nbproject");
   				
   						nbpatient = dataStoreGeneDJV.getValue(item, "nbpatient"); 
   						CnvGeneDejavu = dataStoreGeneDJV.getValue(item, "cnvs"); 
   						
   						//out =   '<button  class="'+ boutonType + '" id="' + boutonID + '" onclick="my_view_dejavu_polycyto(\'' + CnvGeneDejavu + '\')">' +  nbproject + ' : ' + nbpatient + ' </button>';  
   						out =   '<button  class="'+ boutonType + '" id="' + boutonID + '>' +  nbproject + ' : ' + nbpatient + ' </button>'; 
   				}
  			});	
  
   			return out;
}

//--------------------------------------------------------------------
// pour polycyto  global view
//---------------------------------------------------------------------

var dataStore_RESUM;
var url_RESUM_Project = url_path + "/manta/PolyCytoGlobalView.pl";


function GetResum() 
{
		var parameters = location.search.split("&");
		var arg1 = parameters[0].split("=");
   		var arg2 = parameters[1].split("=");
		
    	var projectname = arg1[1];
    	var filename = arg2[1];
   
		dataStore_RESUM = new dojo.data.ItemFileReadStore({ url: url_RESUM_Project+ "?projectname=" + projectname + "&filename=" + filename});

	var gridRESUM = dijit.byId("RESUM_gridDetails");
	gridRESUM.setStore(dataStore_RESUM);
	gridRESUM.store.fetch({
		onComplete: function(items){		
		}
	});
	gridRESUM.startup();
	return;
}

function launch_resum() 
{
	return;
	var parameters = location.search.split("&");
	var arg1 = parameters[0].split("=");
  	var arg2 = parameters[1].split("=");
		
  	var projectname = arg1[1];
  	var patientname = arg2[1];
  
 	var url_viewer = "PolyCytoGlobalView.html?projectname=" + projectname + "&patientname=" + patientname;
	var myWindow = window.open(url_viewer,"_blank",""); 
 
 }
 
 function formaterEVENT(value)
{
	var eurl;
	 if ( /DEL/.test(value)) {
	 	eurl = '<font style="color:red"> '+ value +' </font>';
	 }
	 
	 if ( /DUP/.test(value)) {
	 	eurl = '<font style="color:blue"> '+ value +' </font>';
	 }
	 
	 if ( /t/.test(value)) {
			eurl = '<font style="color:DarkOrange"> '+ value +' </font>';	 
	}
	 
	 if ( /inv/.test(value)) {
	 	eurl = '<font style="color:indigo"> '+ value +' </font>';
	 }
	 
	return eurl;
}

//-------------------------------
// En plus
//--------------------------------


function sleep(milliseconds) {
  const start = Date.now();
  while (Date.now() - start < milliseconds);
}

function formatNumber(num) 
{
  return num.toString().replace(/(\d)(?=(\d{3})+(?!\d))/g, '$1,');
}

function changeAccessFilters()
{
	if ( document.getElementById('scoreCNVSelect').disabled == true )
	{
			accessCNVFilters("On");
	}
	else
	{
			accessCNVFilters("Off");
	}
	//GetAllCNV();
}
			

function accessCNVFilters(val)
{
	if (val == "Off")
	{
		// valeurs par defaut pour les bestones
		document.getElementById('scoreCNVSelect').style = "font-size:10px;color:DarkOrange";
		document.getElementById('scoreCNVSelect').value = 1;
		document.getElementById('minlengthSelect').style = "font-size:10px;color:DarkOrange";
		document.getElementById('minlengthSelect').value = 10000;
		document.getElementById('maxlengthSelect').value = "nomax";	
		document.getElementById('transmissionSelect').value = "all";
		document.getElementById('maxfreqSelect').value = "nomax";
		document.getElementById('chromSelect').value = "all";
		document.getElementById('dejavu').style = "font-size:10px;color:DarkOrange";
		document.getElementById('dejavu').value = 10;
		document.getElementById('cytobSelect').value = "all";
		document.getElementById('geneSelect').value = "all";
		document.getElementById('OMIMSelect').checked = false;

		// desactivation	
		document.getElementById('scoreCNVSelect').disabled = true;
		document.getElementById('minlengthSelect').disabled = true;
		document.getElementById('maxlengthSelect').disabled = true;	
		document.getElementById('transmissionSelect').disabled = true;
		document.getElementById('maxfreqSelect').disabled = true;
		document.getElementById('chromSelect').disabled = true;
		document.getElementById('dejavu').disabled = true;
		document.getElementById('cytobSelect').disabled = true;
		document.getElementById('geneSelect').disabled = true;
		document.getElementById('OMIMSelect').disabled = true;
	}
	else
	{
		if (val == "On")
		{
			document.getElementById('scoreCNVSelect').style = "font-size:10px;color:black";
			document.getElementById('scoreCNVSelect').value = 0;
			document.getElementById('minlengthSelect').style = "font-size:10px;color:black";
			document.getElementById('scoreCNVSelect').disabled = false;		
			document.getElementById('minlengthSelect').value = 1000;
			document.getElementById('minlengthSelect').disabled = false;			
			document.getElementById('maxlengthSelect').disabled = false;	
			document.getElementById('transmissionSelect').disabled = false;
			document.getElementById('maxfreqSelect').disabled = false;
			document.getElementById('chromSelect').disabled = false;
			document.getElementById('dejavu').style = "font-size:10px;color:DarkOrange";
			document.getElementById('dejavu').disabled = false;
			document.getElementById('cytobSelect').disabled = false;
			document.getElementById('geneSelect').disabled = false;
			document.getElementById('OMIMSelect').disabled = false;
		}
	}
}

function GetChrsPloidy() 
{
	var parameters = location.search.split("&");
	var par0 = parameters[0].split("=");
	var par1 = parameters[1].split("=");
   	
   	var project =par0[1];
   	var patient = par1[1];
   	
   	url_chrs = url_path + "/manta/PolyCytoGetPloidy.pl";
   	var out;
   	
	dataStorePloidyValue = new dojo.data.ItemFileReadStore({ url: url_chrs + "?project=" + project + "&patient=" + patient});
	
	
	dataStorePloidyValue.fetch({
			onComplete: function(items,request){
				
				var out = '<table class="table-responsive" ><center><tr><td style="padding:5px;border-top-color:white;"> <label style="font-size:16px; color:purple;"> Chromosome View </label></td>';
				
				for (i=0;i<items.length;i++)
				{
                	var item = items[i];
   					out = out + dataStorePloidyValue.getValue(item, "chrmean"); 
   				}
   				out = out + '<td ><div><button  type="button" class="btn-sm  btn-warning" onclick="launch_plotAllChr();" >   ALL  </button></div></td>';
   				out = out + '</tr><center></table>';
  				document.getElementById('chrView').innerHTML = out;
   		}
  	});	
}

function GetChrsDisomie() 
{
	var parameters = location.search.split("&");
	var par0 = parameters[0].split("=");
	var par1 = parameters[1].split("=");
   	
   	var project =par0[1];
   	var patient = par1[1];
   	
   	url_chrs = url_path + "/manta/PolyCytoCheckUniDisomie.pl";
 	
	var dataStoreDisomieValue = new dojo.data.ItemFileReadStore({ url: url_chrs + "?project=" + project + "&patient=" + patient});
	
	dataStoreDisomieValue.fetch({
			onComplete: function(items,request){
						
				var out;
				var item;
				var thevalue;
				
				item = items[0];
				thevalue = dataStoreDisomieValue.getValue(item, "value");
				
				out = '<br>' + thevalue;
				
				for (i=1;i<items.length;i++)
				{
                	item = items[i];
   					out = out + dataStoreDisomieValue.getValue(item, "value"); 
   				}
   				
   				out = out + '</tr><center></table>';
  				document.getElementById('chrView_disomie').innerHTML = out;
  		
   			}
  		});	
}


function get_gene_phenotypes2(this_proj_name,gene_id ) {
	var this_url = url_path + "/json_output_nodb/getGenePhenotypes.pl";
	var this_args;
	if (this_proj_name) { this_args = "project=" + this_proj_name + "&gene=" + gene_id; }
	else { this_args = "gene=" + gene_id; }
	dijit.byId('dialog_list_phenotypes2').show();
	
	var xhr = new XMLHttpRequest();
	xhr.open('POST', this_url, true);
	xhr.setRequestHeader('Content-type', 'application/x-www-form-urlencoded');
	xhr.onreadystatechange = function() {
		if (this.readyState === XMLHttpRequest.DONE && this.status === 200) {
			var data = JSON.parse(this.responseText);
			var store_tmp = new dojo.data.ItemFileWriteStore({
				data: {
					identifier: 'name',
					label: 'name',
					items: data.items
				}
			});
			store_tmp.fetch({
				onItem: function(item, result) {
					document.getElementById("content_table_phenotypes2").innerHTML = item.html;
					if (item.txt_all != '-') {
						add_wordcloud2('container_wordcloud2', item.txt_all);
						document.getElementById("container_wordcloud2").style.display = "block";
					}
					else {
						document.getElementById("container_wordcloud2").style.display = "none";
					}
					
					dijit.byId('dialog_list_phenotypes2').show();
				}
			});
		}
	}
	xhr.send(this_args);
}
function add_wordcloud2(container_name, json_data) {
	document.getElementById(container_name).innerHTML = '';
	anychart.onDocumentReady(function() {
		var data = JSON.parse(json_data);
		var chart = anychart.tagCloud(data);// set a chart title
		chart.title('Most described words in HGMD / HPO / OMIM')
		// set an array of angles at which the words will be laid out
		chart.angles([0, -45, 90]);
		// enable a color range
		chart.colorRange(true);
		// set the color range length
		chart.colorRange().length('80%');// display the word cloud chart
		chart.container(container_name);
		chart.draw();
	});
	return;
}



