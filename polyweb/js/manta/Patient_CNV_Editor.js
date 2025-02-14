
var dataStore_SVCompare;
var gridSVCompare;
var dataStorePlotValue;
var gridPlot;
var tabPos ="";
var tabRatio = "";

var nb=0;






function setTitre()
{
	document.getElementById("waiting").hide;
	var parameters = location.search.split("&");
	var par0 = parameters[0].split("=");
	var par1 = parameters[1].split("=");
   	
   	var project =par0[1];
   	var patient = par1[1];
   	
	titre = '<br>CNV (DUP/DEL) <br> ' + project + '   /   ' + patient + '<br><br>';
	document.getElementById('titre').innerHTML=titre;
}


			
function GetAllCNV(value) 
{
		var parameters = location.search.split("&");
		var arg1 = parameters[0].split("=");
   		var arg2 = parameters[1].split("=");
		
    	var projectname = arg1[1].split(",");
    	var filename = arg2[1].split(",");
    	var thevalue=value;	
		var minlength = document.getElementById('minlengthSelect').value;
		var maxlength = document.getElementById('maxlengthSelect').value;	
		var transmission = document.getElementById('transmissionSelect').value;
		var maxfreq = document.getElementById('maxfreqSelect').value;
		
		// var genotype = document.getElementById('genotypeSelect').value;
		var genotype="both";
		
		var chr = document.getElementById('chromSelect').value;
		var dejavu = document.getElementById('dejavu').value;
		
		var cytoband = document.getElementById('cytobSelect').value;
		if (cytoband == "") { cytoband = "all"};
		
		var genes = document.getElementById('geneSelect').value;
		if (genes == "") { genes = "all"};
		
		
		dataStore_SVCompare = new dojo.data.ItemFileReadStore({ url: url_CNV_Patient+ "?projectname=" + projectname + "&filename=" + filename + "&minlength=" + minlength + "&maxlength=" + maxlength + "&transmission=" + transmission + "&maxfreq=" + maxfreq + "&dejavu=" + dejavu + "&genotype=" + genotype  + "&select_best=" + thevalue + "&chrom=" +chr +"&cytoband=" +cytoband + "&genes=" +genes + "&print=" + 0 });
		
		var gridSVCompare = dijit.byId("SVCompare_gridDetails");
		gridSVCompare.setStore(dataStore_SVCompare);
		gridSVCompare.store.fetch({
			onComplete: function(items){
				if (items.length == 0)
				{
						var out = '<center><i><p style="color:grey;font-size:14px;padding:25px;">  No results for BestOnes... Try with "ALL"  </p></i></center>';
						document.getElementById('dojox_grid__View_1').innerHTML = out;
						
				}
			}
		});
	
		gridSVCompare.startup();
	
		return;
}



function collapse(name){
    var d = document.getElementById(name);
    d.classList.add("in");
}

function collapseResult(){
    var d = document.getElementById("collapse1");
    d.classList.add("in");
    
}

function collapseResult2(){
    var d = document.getElementById("collapse1");
 
    if (d.classList.contains("in")) {
        d.classList.remove("in");
    }
    else {
        d.classList.add("in");
    }
}

function collapse2(name){
    var d = document.getElementById(name);
 
    if (d.classList.contains("in")) {
        d.classList.remove("in");
    }
    else {
        d.classList.add("in");
    }
}



function formaterDUPDEL(value)
{
	var eurl;

	 if ( /DEL/.test(value)) {
	 	eurl = '<b><font style="color:red"> '+ value +' </font></b>';
	 }
	 
	 if ( /DUP/.test(value)) {
	 	eurl = '<b><font style="color:blue"> '+ value +' </font></b>';
	 }
	 
	 if ( /TRI/.test(value)) {
			eurl = '<b><font style="color:blue"> '+ value +' </font></b>';	 
	}
	 
	 if ( /QUA/.test(value)) {
	 	eurl = '<b><font style="color:indigo"> '+ value +' </font></b>';
	 }
	 
	return eurl;
}


function formaterRatio(value)
{
		var out;
		var val;
		var eurl="-";
				
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
				
		if (value == -1)
		{
			return eurl;
		}
		
		eurl = '<b><font style="color:orange">' + value.toFixed(2) + '</font></b>';
		return eurl;	
}	
	
function formaterDUPSEG(value)
{
	
		var eurl = '<b><font size="4" style="color:green"> - </font>';
		
		if (value == "-")
		{
			return eurl;
		}

		var val = parseInt(value);
		
		if (val == 0)
		{
				eurl = '<b><font size="4" style="color:green"> - </font>';
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
	var out=" ";
	var pos = value.split("-");
	
	out = '<font style="color:grey">' + pos[0] + '</font><font style="color:black">-</font><font style="color:black">' + pos[1] + '</font>';
	return out;
}	

function formaterPlot(value)
{
	var val = value.split(":");
	var pos = val[1].split("-");
	
	var chr = val[0];
	
	var deb = pos[0];
	var fin = pos[1];

	var out=" - ";
	
	if ( fin-deb >= 25000)
	{  
			var boutonID = 'boutonPlot'+ nb;				
   			var out   = '<button id="' + boutonID + '" onclick="launch_plot(' + chr + ',' + deb + ',' + fin + ' )"><i class="fa fa-chart-bar" style="font-size:16px; color:orange"></i></button>';
	}
	
	return out;
}	


function formaterLength(value)
{
	var eurl;

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
   			var boutonID = 'boutonIGV'+ nb;					
   			var igv_value = locusstart + ';' + locusend + ';' + bamfiles + ";" + bamnames + ";"+ taille;	
			var etiquette   = '<button class="igvIcon2"  id="' + boutonID + '" onclick="launch_web_igv_start_and_end(\'' + idCNV + "et" + igv_value + '\')"></button>';
	
	return etiquette;
}		


function formaterInfoCaller(value)
{
	var out=""; 
	
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
			etiquette = '<font size = "2" color="orange"> <span class="glyphicon glyphicon-ok-sign"></span></font>';
		}
		
		out = '<a href="#" data-toggle="tooltip" title="' + out + '">' +  etiquette  +'  </a>';
		return out;
	}
}


function formaterScore(value)
{
	var eurl;
	
	var score = value;
	
	 if (score > 3) {
	 	eurl = '<label style="color:red;">' + score + ' </label>';
	 }
	 else {
	 	eurl = '<label style="color:black">'+score + '</label>';
	 }
	 
	return eurl;
	
}	


function formaterTransmission(value)
{
	var eurl="-";
	
	if (/denovo/.test(value) )
	 {
	 	eurl = '<b><i><label style="font-size:10px;color:red"> denovo  </label></i></b>';
	 }
	 
	if ( /mother/.test(value) & !/father/.test(value) )
	{
	 	eurl = '<div class=circle2F><img src="https://img.icons8.com/offices/20/000000/guest-female.png"/></div>';
	 }
	 
	 if  ( /father/.test(value) & !/mother/.test(value) )
	 {
	 	eurl = '<div><img src="https://img.icons8.com/offices/20/000000/person-male.png"/></div>';
	 }
	 
	 if  ( /both/.test(value) | ( /father/.test(value) & /mother/.test(value) ) ) {
	 	eurl = '<div><img src="https://img.icons8.com/offices/20/000000/guest-female.png"/><img src="https://img.icons8.com/offices/20/000000/person-male.png"/></div>';
	 }
	 
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
						out = out + tval[i] + " ";
				}
			}
		}
		return out;	
}


function formaterDejavu(value)
{
	var values = value.split(";");
	
	var djv_itp = values[0];			// nb echantillon
	var liste_djv_itp="";				// liste des echantillons du projet
	
	var liste = values[1];	
	var tab = liste.split(",");
	var len = tab.length-1;
	
	for (i=0;i<len;i++)
	{
		liste_djv_itp = liste_djv_itp + tab[i]+'\n';
	}
	
	var djv_in_nbother_project = values[2];
	var djv_iop = values[3];
	var liste_djv_iop =values[4];
	
	var out;
   	var boutonID = 'boutondjv' + nb;	
		
	out = '<button  class="btn btn-primary btn-xs" id="' + boutonID + '" onclick="my_view_dejavu_polycyto(\'' + value + '\')">' +  djv_in_nbother_project + ' : ' + djv_iop + ' </button>';  
	out = out + '<a href="#" data-toggle="tooltip" title=" ' + liste_djv_itp + ' ">    <button class="btn btn-classic btn-xs">' +  djv_itp + ' </button></a>';
	
	return out;
}


function formaterRANK(value)
{
	var eurl;
	 if (value > 5) {
	 	eurl = '<b><font style="color:red">' + value + ' </font></b>';
	 }
	 else {
	 	eurl = '<b><font style="color:black">'+ value + '</font></b>';
	 }
	 
	return eurl;
	
}	


function formaterGENESOLD(value)
{

	var out=" ";
	var outshort;
	
	var names = value.split(" ");
	var len = names.length-1;
	
	
	if (! /[a-zA-Z0-9]/.test(value))
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
	
	outshort = names[0]+" "+ names[1] + " ... ";
	
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


function formaterGENES(value)
{
	var out;
	var liste="";
	var nbgenes = 0;
	
	var names = value.split(" ");
	var len = names.length;
	
	if (! /[a-zA-Z0-9]/.test(value))
	{
		return "-";
	}
	
	for (i=0; i<len ; i++)
	{
		var name = names[i];
		if (/[a-zA-Z0-9]/.test(name))
		{
			 liste = liste + '<option value=' + i + ' ><b>' + name + '</b></option>';
			 nbgenes++;
		}
	}
		
	liste = liste + '</select>';
	
	if (nbgenes >3)
	{
		out = '<center><select class="form-control" style="width:90%; heigth:50%; font-size:10px;padding:5px"><option value= -1 selected="selected" >' +   names[0] + '  ' +  names[1] + '  ' +  names[2] + '  ... </option>' + liste ;
	}
	else
	{
		if (nbgenes >2)
		{
			out = '<center><select class="form-control" style="width:90%; heigth:50%; font-size:10px;padding:5px"><option value= -1 selected="selected" >' +   names[0] + '  ' +  names[1] + '  ' +  names[2] + ' </option>' + liste ;
		}
		else
		{
				if (nbgenes >1)
				{
					out = '<center><select class="form-control" style="width:90%; heigth:50%; font-size:10px;padding:5px"><option value= -1 selected="selected" >' +   names[0] + '  ' +  names[1]  + '  </option>' + liste ;
				} 
				else
				{
					out = '<center><select class="form-control" style="width:90%; heigth:50%; font-size:10px;padding:5px"><option value= -1 selected="selected" >' +   names[0] +  '  </option>' + liste ;
				}
		}
	}
	

	return out;
}

function formaterOmin(value)
{

	var eurl = "-";
	
	 if (/yes/.test(value) )
	 {
	 	eurl = '<b><font style="color:red">'+value+'</font></b>';
	 }
	return eurl;
}	


function formaterDGV(value)
{
	var eurl = '<b><i><a href="' + value  + '" target="_blank"><font style="color:red;font-size:12px; ">  DGV  </font></a></i></b>';
	return eurl;
}


function formaterFreq(value)
{
	var eurl;

	if (value > 0)
	{
		value = value*100;
		eurl = value.toFixed(3);
		
	 		if (value < 0.05)
			{
	 			eurl = '<font style="color:red">' + eurl + ' % </font>';
	 		}
	 		else
	 		{
	 			if (value < 0.1)
				{
	 				eurl = '<font style="color:orange">' + eurl + ' % </font>';
	 			}
	 			else
	 			{
	 				eurl = eurl + ' %';
	 			}
	 		}	 	
	 		return eurl;
	 }
	 return "-";
}	


function formaterdbvarStatus(value)
{
	
	var eurl=value;
	var out=" ";
	
	var status = value.split(" ");
	var len = status.length;
	
	if (/-/.test(value))
	{
		return "-";
	}
	
	for ( i = 0 ; i < len; i++)
	{
				
				if (status[i] != "-")
				{	
					out = out +"  " + status[i] ;
				}
	}
	
	out = '<a href="#" data-toggle="tooltip" title="' + out + '"><font color="red" size=3> *  </a>';
	return out;
}		


// En plus pour IGV

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
	
	if (taille == 1)
	{
		titre = patient_name;
	}
	else
	{
		titre = patient_name + "  :  " + titre + "  (" + taille + ")";
	}
	
	document.getElementById('titre').innerHTML=titre;
	
	if (genes_bed) { url_gene_bed = genes_bed; }
	if (fasta) { url_fasta = fasta; }
	view_web_igv_bam_simple('div_igv_start', locusstart, bamfile, patient_name);
	view_web_igv_bam_simple('div_igv_end', locusend, bamfile, patient_name);
}

function igv_view_CNV()
{
   	var titre = param("titre");
   	var locusstart = param("locusstart");
   	var locusend = param("locusend");
	var bamfile = param("bamfile");
	var patient_name = param("patient_name");
	var taille = param("taille");
	var genes_bed = param("genes_bed");
	var fasta = param("fasta");
	
	if (taille != 1){
		titre = patient_name + "  :  " + titre + "  (" + taille + ")";
	}
	if (bamfile.match(/HG19/)) { document.getElementById('release').innerHTML = "<span style='color:red;'>Release HG19</span>"; }
	else if (bamfile.match(/HG38/)) { document.getElementById('release').innerHTML = "<span style='color:red;'>!!! Release HG38 !!!</span>"; }
	
	document.getElementById('titre').innerHTML=titre;
	var locus = locusstart + " " + locusend;

	if (genes_bed) { url_gene_bed = genes_bed; }
	if (fasta) { url_fasta = fasta; }
	view_web_igv_bam_simple('div_igv_start',locus, bamfile, patient_name);
	
}

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
	
	
	var url_viewer = "Url_CNV_Viewer.html?titre=" + titre + "&locusstart=" + locusstart + "&locusend=" + locusend + "&bamfile=" + bamfile + "&patient_name=" + patient_name + "&taille=" + taille;
	if (url_gene_bed) { url_viewer += "&genes_bed=" + url_gene_bed; }
	if (url_fasta) { url_viewer += "&fasta=" + url_fasta; }
	var myWindow = window.open(url_viewer,"_blank",""); 
 
 }
 
 
function formatNumber(num) {
  return num.toString().replace(/(\d)(?=(\d{3})+(?!\d))/g, '$1,')
}



function launch_plot(chr,deb,fin)
{

	var parameters = location.search.split("&");
	var par0 = parameters[0].split("=");
	var par1 = parameters[1].split("=");
	
	var project =par0[1];
   	var patient = par1[1];
 
	var debplot = deb - 100000;
	var finplot = fin + 100000;
	

	var url_plot= "Url_CNV_Plot.html?project=" + project + "&patient=" + patient +"&chr=" + chr + "&deb=" + debplot + "&fin=" +finplot;
	var myWindow = window.open(url_plot,"_blank",""); 

}

function PlotValues()
{

	var parameters = location.search.split("&");
	var par0 = parameters[0].split("=");
	var par1 = parameters[1].split("=");	 
	 var par2 = parameters[2].split("=");
	var par3 = parameters[3].split("=");	 
	var par4 = parameters[4].split("=");
	
	var project =par0[1];
	var patient =par1[1];
	var chr=par2[1];
	var debplot = par3[1];
	var finplot = par4[1];
	
	var url_plotCNV = url_path + "/manta/plotCNV.pl";
	
	
	dataStorePlotValue = new dojo.data.ItemFileReadStore({ url: url_plotCNV + "?project=" + project + "&patient=" + patient + "&chr=" + chr + "&deb=" + debplot + "&fin=" + finplot });
	
	
	dataStorePlotValue.fetch({
			onComplete: function(items,request){
  					var i;
                	var item = items[0];
   
   					tabPos = dataStorePlotValue.getValue(item, "POS"); 
   					tabRatio = dataStorePlotValue.getValue(item, "RATIO");
  
   					var tabX = tabPos.split(" ");
   					var tabY = tabRatio.split(" ");
   						
					var trace1 = {
  								x:tabX,
  								y:tabY,
  								mode: 'markers',
					};
						
					var data = [trace1];	
					Plotly.newPlot('myDiv', data);
					
					
   			}
   	});
}

function my_view_dejavu_polycyto(value) {
	dijit.byId('waiting').show();
	var lCol = value.split(';');
	var infos_proj_pat = lCol[-1];
	var build = 'HG19';
	var ltmp = String(value).split(';');
	var args = 'login=' + dojo.cookie("username") + '&pwd=' + dojo.cookie("passwd") + '&values=' + ltmp[4] + '&my_phenotype="" ';
	var xhr = new XMLHttpRequest();
	xhr.open('POST', url_proj_pat_details, true);
	xhr.setRequestHeader('Content-type', 'application/x-www-form-urlencoded');
	xhr.onreadystatechange = function() {
	    if (this.readyState === XMLHttpRequest.DONE && this.status === 200) {
	    	var data = JSON.parse(this.responseText);
			document.getElementById("content_dejavu_polycyto").innerHTML = data.html_table;
			dijit.byId('waiting').hide();
			dijit.byId('dialog_dejavu_polycyto').show();
	    }
	}
	xhr.send(args);
}





	


	








