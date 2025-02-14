
var dataStorePlotValue;
var gridPlot;
var tabPos ="";
var tabRatio = "";

var nb=0;


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

function formaterDUPDEL(value)
{

	var eurl;
	
	if (value == "-")
	{   
		return value;
	}
	
	 if ( /DEL/.test(value)) {
	 	//eurl = '<b><font style="color:red"> '+ value +' </font></b>';
	 	eurl = '<p style="font-size:10px;color:white;background-color:rgb(255,74,42);padding:8px">' + value + ' </p>';
	 }
	 
	 if ( /DUP/.test(value)) {
	 	//eurl = '<b><font style="color:blue"> '+ value +' </font></b>';
	 	eurl = '<p style="font-size:10px;color:white;background-color:rgb(12,148,230);padding:8px">' + value + ' </p>';
	 }
	 
	 if ( /TRI/.test(value)) {
			//eurl = '<b><font style="color:blue"> '+ value +' </font></b>';	 
			eurl = '<p style="font-size:10px;color:white;background-color:rgb(12,98,228);padding:8px">' + value + ' </p>';	
	}
	 
	 if ( /QUA/.test(value)) {
	 	//eurl = '<b><font style="color:indigo"> '+ value +' </font></b>';
	 	eurl = '<p style="font-size:10px;color:white;background-color:rgb(1,26,255);padding:8px">' + value + ' </p>';	
	 }
	 
	return eurl;
}


function formaterRatio(value)
{
		var out;
		var val;
		var eurl="-";
		
		if (value == "-")
		{
			return value;
		}		
		
		if (value == 0)
		{
			return eurl;
		}

		if (value < 0)
		{
				eurl = '<b><font style="color:red">' + value.toFixed(2) + '</font></b>';
		}
		if (value > 0)
		{
				eurl = '<b><font style="color:blue">' + value.toFixed(2) + '</font></b>';
		}
		return eurl;	
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
	
		var eurl = '<b><font size="4" style="color:green"> - </font></b>';

		if (value == "-")
		{
			return value;
		}

		var val = parseInt(value);
		
		if (val == 0)
		{
				eurl = '<b><font size="4" style="color:green"> - </font></b>';
				return eurl;
		}
		
		
		if (val < 40)
		{
				eurl = '<font style="color:green">' + val + ' % </font>';
		}
		if (val >= 40)
		{
				eurl = '<font style="color:red">' + val +  ' % </font>';
		}
		return eurl;	
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
	
	if (value == 23)
	{
		 eurl ="X";
	}
	if (value == 24)
	{
		 eurl="Y";
	}
	
	eurl = '<b>' + eurl + '</b>';
	
	return eurl;
}

function formaterCytoBand(value)
{
	var out=" ";
	var outshort;
	
	var names = value.split(",");
	var len = names.length;

	if (value == "-")
	{
		return value;
	}


	if (value == "--")
	{
		return "-";
	}
	
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
			out = '<font style="color:black">'+ out + '</font>';
	}
	else
	{
			out = '<a href="#" data-toggle="tooltip" title="' + out + '"> ' + outshort + '</a>';
	}
	return out;
}	

function formaterLocus(value)
{
	if (value == "-")
	{
		return value;
	}
	
	var out=" ";
	var pos = value.split("-");
	out = '<font style="color:grey">' + pos[0] + '</font><font style="color:black">-</font><font style="color:black">' + pos[1] + '</font>';
	return out;
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
   	var out   = '<button class="btn btn-classic btn-m" style="border: 1px solid black;padding:3px" id="' + boutonID + '" onclick="launch_plot(\''+ transmission + '\',' + type + ',' + chr + ',' + debcnv + ',' + fincnv + ' )" style="font-size:16px;">&#128200;</button>';
  
	return out;
}	


function formaterLength(value)
{
	var eurl;
	
	if (value == "-")
	{
		return value;
	}
	
	 if (value > 1000000 ) {
	 	value = value / 1000000;
	 	eurl = '<b>' + value.toFixed(2)  + ' Mb </b>';
	 	return eurl;
	 }
	 
	 if (value > 1000 ) {
	 	value = value / 1000;
	 	eurl = value.toFixed(2)  + ' kb';
	 	return eurl;
	 }	
	 
	 if (value < 1000 ) {
	 	eurl = value  + " pb";
	 	return eurl;
	 }	
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
			
	return etiquette;
}		


function formaterInfoCaller(value)
{
	var out=""; 

	if ( (value == "Sorry") || (value == "No") || (value == "Results") )
	{
		return '<i><font size = "2" color="red">' + value + '</font></i>';
	}
	
	if (value == "-")
	{
		return  value;
	}
	
	if (value == 0)
	{
		etiquette = '<font size = "2" color="red"><span class="glyphicon glyphicon-minus-sign"></span></font>';
		return etiquette;
	}
	else
	{
		
		var tvalue = value.split(";");
		var val = tvalue[0];
		var qual = tvalue[1];
		
	
		var tval = val.split(" ");
		
		var len = tval.length-1;
	
		for ( i = 0 ; i < len-1 ; i++)
		{
				out = out +  tval[i] + " / ";
				
		}
		out = out + tval[i];
		
		if (qual == 1) 
		{
			etiquette = '<font size = "2" color="green"> <span class="glyphicon glyphicon-ok-sign"></span></font>';
		}
		
		if (qual == 0.75)
		{
			etiquette = '<font size = "2" color="lightgreen"> <span class="glyphicon glyphicon-ok-sign"></span></font>';
		}
		
		out = '<a href="#" data-toggle="tooltip" title="' + out + '">' +  etiquette  +'  </a>';
		return out;
	}
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

function formaterScoreCNV(value)
{

	if (value <= -2 )
	{
		return "-";
	}


	var eurl;
	
	var score = Number(value);
	
	 	if (score >= 1 ) {
	 			eurl = '<label style="color:red">' + score + ' </label>';
		 }
	 	else {
	 			eurl = '<label style="color:black">'+ score + '</label>';
	 	}

	 
	return eurl;
	
}	

function formaterScoreGene(value)
{

	if (value <= -2 )
	{
		return "-";
	}


	var eurl;
	
	var score = Number(value);
	
	 	if (score >= 4 ) {
	 			//eurl = '<label style="color:red">' + score + ' </label>';
	 			eurl = '<label style="font-size:10px;color:white;background-color:rgb(255,0,0);padding:3px">' + score + ' </label>';
		}
	 	else 
	 	{
	 		if (score >= 3 ) 
	 		{
	 			//eurl = '<label style="color:Tomato">' + score + ' </label>';
	 			eurl = '<label style="font-size:10px;color:white;background-color:rgb(255,127,79);padding:3px">' + score + ' </label>';	 			
			}
			else
			
				 if (score >= 2 ) 
	 			{
	 				//eurl = '<label style="color:orange">' + score + ' </label>';
	 				eurl = '<label style="font-size:10px;color:white;background-color:rgb(255,165,0);padding:3px">' + score + ' </label>';	 				
				}
				else
				{
					if (score >= 1 ) 
	 				{
	 					//eurl = '<label style="color:black;background-color:yellow">' + score + ' </label>';
	 					eurl = '<label style="font-size:10px;color:black;background-color:rgb(255,255,2);padding:3px">' + score + ' </label>';
					}
					else
					{
	 					//eurl = '<label style="color:black">'+ score + '</label>';
	 					eurl = '<label style="font-size:10px;color:black;background-color:rgb(204,204,204);padding:3px">' + score + ' </label>';
	 				}
	 			}
	 	}
	 	
	return eurl;
	
}	


function formaterTransmission(value)
{

	if (value == "-")
	{
		return value;
	}
	
	if ( value == "X" )
	{
			out ='<b><label  style="color:grey"> X </label></b>';
			return out;
	}
	
	var TransM="";
	var TransF="";
	var valM="";
	var valF="";
	var res;
	
	var tval = value.split(" ");
	var len = tval.length;
		
	for ( i = 1 ; i < len-1 ; i++)
	{
		var val = tval[i];
		
		if ( val == "mother")
		{
			TransM =  '<img src="../../images/polyicons/icons8-person-female-24.png" >';
		}
		if (/maybeM/.test(val))
		{
			valM = val.replace("maybeM","mother");
		}
		
		if ( val == "father")
		{
			TransF = '<img src="../../images/polyicons/icons8-person-24.png" >';
		}
		if (/maybeF/.test(val))
		{
			valF = val.replace("maybeF","father");
		}
	}
	
	if (TransM=="" && valM != "" )
	{ 
		TransM = '<b><i><label style="font-size:10px;color:pink">' + valM  + '</label></i></b>';
	}
	
	if (TransF=="" && valF != "" )
	{ 
		TransF = '<b><i><label style="font-size:10px;color:lightblue">' + valF  + '</label></i></b>';
	}
	
	if (TransM=="" && TransF == "" )
	{
		res = '<b><i><label style="font-size:10px;color:red"> strict-denovo  </label></i></b>';;
	}
	else
	{
		res = TransM+TransF;
	}
	
	var eurl = '<div>' + res +  '</div>';
	
	
	 return eurl;
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

function formaterGENESNew(value)
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
	
	var liste0="";
	var liste1="";
	var liste2="";
	var liste3="";
	var liste4="";

	
	var etiq0="";
	var etiq1="";
	var etiq2="";
	var etiq3="";
	var etiq4="";

	
	var nb0=0;
	var nb1=0;
	var nb2=0;
	var nb3=0;
	var nb4=0;


	locus = values[0];

	for ( i = 1 ; i < len; i++)
	{
		if (values[i] != "-")
		{	
			var geneinfos = values[i].split(";");
			name = geneinfos[0];
			score = geneinfos[1];
			
			if (score < 1 )
			{
				nb0++;
				if(nb0 < 4)
				{
					liste0 = liste0 + '<font style="font-size:10px;color:white;background-color:rgb(128,128,128);padding:1px;border: 1px solid white">' + name + '</font><br>';
				}
			}
			else
			{
				if (score < 2 )
				{
					nb1++;
					if(nb1 < 4)
					{
						liste1 = liste1 + '<font style="font-size:10px;color:black;background-color:rgb(255,255,2);padding:1px;border: 1px solid black">' + name +'</font><br>';
					}
				}
				else
				{
					if (score < 3 )
					{
						nb2++;
						if(nb2 < 4)
						{
							liste2 = liste2 + '<font style="font-size:10px;color:white;background-color:rgb(255,165,0);padding:1px;border: 1px solid black">' + name +'</font><br>';
						}
					}
					else
					{
						if (score < 4 )
						{
							nb3++;
							if(nb3 < 4)
							{
								liste3 = liste3 + '<font style="font-size:10px;color:white;background-color:rgb(255,127,79);padding:1px;border: 1px solid black">' + name + '</font><br>';
							}
						}
						else
						{
						
							nb4++;
							if(nb4 < 6)
							{
								liste4 = liste4 + '<font style="font-size:10px;color:white;background-color:rgb(255,0,0);padding:1px;border: 1px solid black">' + name + '</font><br>';
							}						
						}
					}
				}
			}
		}
	}
	

	
	if(liste4 != "")
	{
		outshort = liste4;
		nb4 -= 5;
	}
	else
	{
		if(liste3 != "")
		{
			outshort = liste3;
			nb3 -= 3;

		}
		else
		{
			if(liste2 != "")
			{
				outshort = liste2;
				nb2 -= 3;
			}
			else
			{
				if(liste1 != "")
				{
					outshort = liste1;
					nb1 -= 3;
				}
				else
				{
					outshort = liste0;
					nb0 -= 3;
				}
			}
		}
	}
		
	out = out + '<button class="btn btn-classic btn-xs" style="color:black;border: 1px solid black;padding:2px" onclick="viewGeneByName(\'' + locus + '\', \''+ project + '\',\'' + patient + '\')">';
	out = out + outshort;
	
	if( nb0+nb1+nb2+nb3+nb4 > 0)
	{
	
	out = out + '<center width="80%"> <table><tr>';
	
	if(nb0>0)
	{
		etiq0 =  '<font style="font-size:9px;color:white;background-color:rgb(128,128,128);padding:3px"> +' + nb0 + '</font>';
		out = out + '<td width="20%" style="padding:2px"><center>' +  etiq0  + '</center></td>';
	}
	if(nb1>0)
	{
		etiq1 = '<font style="font-size:9px;color:black;background-color:rgb(255,255,2);padding:3px"> +' + nb1 +'</font>';
		out = out + '<td width="20%" style="padding:2px"><center>' +  etiq1  + '</center></td>';		
	}
	if(nb2>0)
	{
		etiq2 = '<font style="font-size:9px;color:white;background-color:rgb(255,165,0);padding:3px"> +' + nb2 +'</font>';
		out = out + '<td width="20%" style="padding:2px"><center>' +  etiq2  + '</center></td>';
	}
	if(nb3>0)
	{
		etiq3 = '<font style="font-size:9px;color:white;background-color:rgb(255,127,79);padding:3px"> +' + nb3 +'</font>';
		out = out + '<td width="20%"style="padding:2px" ><center>' +  etiq3  + '</center></td>';
	}
	if(nb4>0)
	{	
		etiq4 = '<font style="font-size:9px;color:white;background-color:rgb(255,0,0);padding:3px"> +' + nb4 +'</font>';
		out = out + '<td width="20%" style="padding:2px"><center>' +  etiq4  + '</center></td>';
	}

	out = out + '</tr></table></center>'; 
	}
	
	out = out + '</button>';
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
	
	var infos = value.split("##");
	var eurl;
	
	if (infos[0] == "inv" )
	{
		eurl = '<p style="font-size:12px;color:purple;background-color:rgb(204,204,204);padding:8px"> <b>inv</b> (' + infos[1] + ',' + infos[3] + ') (' + infos[2] + ',' + infos[4] + ') </p>';
		
	}
	
	else
	{
		eurl = '<p style="font-size:12px;color:white;background-color:rgb(255,165,0);padding:8px"> <b>t </b> (' + infos[1] + ',' + infos[3] + ') (' + infos[2] + ',' + infos[4] + ') </p>';
	}
	
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
		
		alert("coucou_post");
		
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
	
	// pour lancer web_igv
	var url_viewer = "PolyCyto_IGV_Viewer.html?titre=" + titre + "&locusstart=" + locusstart + "&locusend=" + locusend + "&bamfile=" + bamfile + "&patient_name=" + patient_name + "&taille=" + taille;
	if (url_gene_bed) { url_viewer += "&genes_bed=" + url_gene_bed; }
	if (url_fasta) { url_viewer += "&fasta=" + url_fasta; }
	var myWindow = window.open(url_viewer,"_blank",""); 
 
 	// pour lancer appli IGV si elle est ouverte
	var tab_bamfile = bamfile.split(',');
	displayListBamInIgvApp(tab_bamfile); 
		
	var locus =locusstart+ ' ' +locusend;
    displayInIGVLocus(locus);
	
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
	view_web_igv_bam_simple('div_igv_start',locus, bamfile, patient_name);
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
	
	var args = value.split("+");
	var event_id=args[0];
	var values = args[1].split(";");
	
	var djv_itp = values[0];			// nb echantillon
	var liste_djv_itp="";				// liste des echantillons du projet
	
	var liste = values[1];	
	var tab = liste.split(",");
	var len = tab.length-1;
	
	for (i=0;i<len;i++)
	{
		liste_djv_itp = liste_djv_itp + tab[i]+' ';
	}
	
	var djv_in_nbother_project = values[2];
	var djv_in_nbother_patient = values[3];
	var liste_djv_iop =values[4];
	var nbwisecondor=parseInt(values[5]);
	var nbcanvas=values[6];
	var nbmanta=values[7];
	var flagwc =values[8];
	
	var out;
   	var boutonID = 'boutondjv' + nb;	
   	
 
   	// pour l'affichage
   var couleur ="white";
   var bgcolor = "#E6E6E6";
   
   var boutonType= "btn btn-classic btn-xs";
   
   if ( (flagwc == 1) && (nbwisecondor < dejavu_seuil) )
   {
   		bgcolor="#FFBF00;color:black;";
   
  		 if ( nbwisecondor > (dejavu_seuil/2) )
   		{
   			bgcolor="#F5DA81;color:black;";
  		 }
  	}
  
   
	out =   '<center> <table  width="100%"><tr style="font-size:9px;">';
	out = out + '<td><center><button style="width:90%; background-color:' + bgcolor + ';"  id="' + boutonID + '" onclick="my_view_dejavu_polycyto(\'' + value + '\')">';
	out = out + '<table width="100%"><tr style="font-size:9px">';
	out = out + '<td  width="20%" style="padding:4px;border: 1px solid ' + couleur + ';"><center>  Pro </center></td>';
	out = out + '<td width="20%" style="padding:4px;border: 1px solid '  + couleur + ';"><center> Pat </center></td>';
	out = out + '<td width="20%"style="padding:4px;border: 1px solid ' + couleur + ';"><center> Wis  </center></td>';
	out = out + '<td width="20%" style="padding:4px;border: 1px solid ' + couleur + ';"><center> Can </center></td>';
	out = out + '<td width="20%"style="padding:4px;border: 1px solid ' + couleur + ';"><center> Man  </center></td></tr>';
	
	out = out + '<tr style="font-size:9px"><td style="border: 1px solid ' + couleur + ';">' + djv_in_nbother_project + '</td>';
	out = out + '<td style="border: 1px solid ' + couleur + ';">' + djv_in_nbother_patient + ' </td>';
	out = out +  '<td style="border: 1px solid ' + couleur + ';">' +  nbwisecondor + '</td>';
	out = out + '<td style="border: 1px solid ' + couleur + ';">' + nbcanvas + '</td>';
	out = out +' <td style="border: 1px solid ' + couleur + ';">' + nbmanta + ' </td>';
	out = out + '</tr></table></button></center></td>'; 
	out = out + '<td style="padding:4px;"><a href="#" data-toggle="tooltip" title=" ' + liste_djv_itp + ' "><center><button class="btn btn-classic btn-xs" style="width=10%;border: 1px solid black;">' +  djv_itp + ' </button></center></a></td></tr></table></center>';
	return out;
}

function formaterDejavuSVeq(value)
{
	//alert(value);
	if (value == "-")
	{
			return value;
	}
	
	var dejavu_seuil = document.getElementById('dejavu').value;
	
	var args = value.split("+");
	var event_id=args[0];
	

	var values = args[1].split(";");
	
	var djv_itp = values[0];			// nb echantillon
	var liste_djv_itp="";				// liste des echantillons du projet
	
	var liste = values[1];	
	var tab = liste.split(",");
	var len = tab.length-1;
	
	for (i=0;i<len;i++)
	{
		liste_djv_itp = liste_djv_itp + tab[i]+' ';
	}
	
	var djv_in_nbother_project = values[2];
	var djv_in_nbother_patient = values[3];
	var liste_djv_iop =values[4];
	
	var out;
   	var boutonID = 'boutondjv' + nb;	
   	
   	
   
   	// pour l'affichage
   	if ( parseInt(djv_in_nbother_project) > parseInt(dejavu_seuil))
   	{
   		djv_in_nbother_project = '+' + dejavu_seuil;
   	}
   	if (parseInt(djv_in_nbother_patient) > parseInt(dejavu_seuil))
   	{
   		djv_in_nbother_patient = '+' + djv_in_nbother_patient;
   	}
   		
 
   	var boutonType= "btn btn-primary btn-xs";
   
	out =   '<button  class="'+ boutonType + '" id="' + boutonID + '" onclick="my_view_dejavu_polycyto(\'' + value + '\')">' +  djv_in_nbother_project + ' : ' + djv_in_nbother_patient + ' </button>';  
	out = out + '<a href="#" data-toggle="tooltip" title=" ' + liste_djv_itp + ' ">    <button class="btn btn-classic btn-xs">' +  djv_itp + ' </button></a>';
	
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







