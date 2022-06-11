var dataStore_SV;
var gridSV;

var dataStore_SVCompare="todo";
var gridSVCompare;
var dataStore_SVFile;


var nb=0.;


var url_allSV_Project = url_path + "/manta/Retrieve_allSV_Project.pl";
var url_BND_Project = url_path + "/manta/Retrieve_Final_BND.pl";


function setTitreBND()
{

	var parameters = location.search.split("&");
	var par0 = parameters[0].split("=");
	var par1 = parameters[1].split("=");
   	
   	var project =par0[1];
   	var patient = par1[1];
   	
	titre = '<br>BreakPoint Viewer  <br> ' + project + '   /   ' + patient + '<br>';
	document.getElementById('titre').innerHTML=titre;
}


			

function GetBNDProject() {

		var parameters = location.search.split("&");
		var arg1 = parameters[0].split("=");
   		var arg2 = parameters[1].split("=");
		
    	var projectname = arg1[1];
    	var filename = arg2[1];
    	
	
		dataStore_BND = new dojo.data.ItemFileReadStore({ url: url_BND_Project+ "?projectname=" + projectname + "&filename=" + filename});

	var gridBND = dijit.byId("BND_gridDetails");
	gridBND.setStore(dataStore_BND);
	gridBND.store.fetch({
		onComplete: function(items){		
		}
	});
	gridBND.startup();
	return;
}


function setTitre()
{

	var parameters = location.search.split("&");
	var par0 = parameters[0].split("=");
	var par1 = parameters[1].split("=");
   	
   	var project =par0[1];
   	var patient = par1[1];
   	
	titre = '<br>CNV (DUP/DEL) <br> ' + project + '   /   ' + patient + '<br>';
	document.getElementById('titre').innerHTML=titre;
}


			

function GetAllSVProject(value) {

		var parameters = location.search.split("&");
		var arg1 = parameters[0].split("=");
   		var arg2 = parameters[1].split("=");
		
    	var projectname = arg1[1].split(",");
    	var filename = arg2[1].split(",");
    	var thevalue=value;	
	    var hide_sd = 1;
		var minlength = document.getElementById('minlengthSelect').value;
		var maxlength = document.getElementById('maxlengthSelect').value;
		var maxfreq = document.getElementById('maxfreqSelect').value;
		var genotype = document.getElementById('genotypeSelect').value;
		var chr = document.getElementById('chromSelect').value;
		var dejavu = document.getElementById('dejavu').value;
		
		var cytoband = document.getElementById('cytobSelect').value;
		if (cytoband == "") { cytoband = "all"};
		
		var genes = document.getElementById('geneSelect').value;
		if (genes == "") { genes = "all"};
		
	
		dataStore_SVCompare = new dojo.data.ItemFileReadStore({ url: url_allSV_Project+ "?projectname=" + projectname + "&filename=" + filename + "&minlength=" + minlength + "&maxlength=" + maxlength + "&maxfreq=" + maxfreq + "&dejavu=" + dejavu + "&genotype=" + genotype + "&transmission=" + 0 + "&hide_sd=" + hide_sd + "&select_best=" + thevalue + "&chrom=" +chr +"&cytoband=" +cytoband + "&genes=" +genes});

	var gridSVCompare = dijit.byId("SVCompare_gridDetails");
	gridSVCompare.setStore(dataStore_SVCompare);
	gridSVCompare.store.fetch({
		onComplete: function(items){		
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


 function formaterSelection(value)
{
	var out;
	var eurl = value.split("_");
	
	var selectID = eurl[0];
	var selection = eurl[1];
	
	if (selection == 1) 
	{
		out = 'X';
	}
	else
	{
		out = '-';
	}
		
	 return out;
}

function formaterDGV(value)
{
	var eurl = '<b><u><a href="' + value  + '" target="_blank"> DGV </a></u></b>';
	return eurl;
}
	
	
function formaterDUPDEL(value)
{
	var eurl;
	
	 if ( /DEL/.test(value)) {
	 	eurl = '<font style="color:red"> '+ value +' </font>';
	 }
	 
	 if ( /DUP/.test(value)) {
	 	eurl = '<font style="color:noir"> '+ value +' </font>';
	 }
	 
	 if ( /INV/.test(value)) {
	 	eurl = '<p style="color:green"> '+ value +' </p>';
	 }
	 
	 if ( /INS/.test(value)) {
	 	eurl = '<p style="color:orange"> '+ value +' </p>';
	 }
	 
	 if ( /LOH/.test(value)) {
	 	eurl = '<p style="color:pink"> '+ value +' </p>';
	 }
	 
	return eurl;
}


function formaterLength(value)
{
	var eurl;
	
	 if (value > 1000000 ) {
	 	value = Math.round(value / 1000000);
	 	eurl = value  + " Mb";
	 	return eurl;
	 }
	 
	 if (value > 1000 ) {
	 	value = Math.round(value / 1000);
	 	eurl = value  + " kb";
	 	return eurl;
	 }	
	 
	 if (value < 1000 ) {
	 	eurl = value  + " pb";
	 	return eurl;
	 }	
}	


function formaterPos(value)
{
	var tval = value.split("-");

	var outshort = tval[0];
	var out = '<a href="#" data-toggle="tooltip" title="' + value + '"> ' + outshort + '</a>';
	return out;
}

function formaterFreq(value)
{
	var eurl;
	
	if (value != "-")
	{
		value = value*100;
		eurl = value.toPrecision(3);
		
	 		if (value < 0.05)
			{
	 			eurl = '<b><font style="color:red">' + eurl + ' % </font></b>';
	 		}
	 		else
	 		{
	 			if (value < 0.1)
				{
	 				eurl = '<b><font style="color:orange">' + eurl + ' % </font></b>';
	 			}
	 			else
	 			{
	 				eurl = '<b>' + eurl + ' % </b>';
	 			}
	 		}	 	
	 		return eurl;
	 }
	 return value;
}	


function formaterDUPSEG(value)
{
	var eurl = "-";
	
	 if (value != "-")
	 {
	 			eurl =  value.toPrecision(3) + ' % ';
	 }

	return eurl;
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
	return eurl;
}


function formaterCHROM_BND1(value)
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
	//return '<b><span class="badge badge-pill badge-warning" style="font-size: 10px;">' + eurl + '</span> '
	return '<span style="font-size: 12px;">' + eurl + '</span> '
}

function formaterCHROM_BND2(value)
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
	//return '<b><span class="badge badge-pill badge-success" style="font-size: 10px;">' + eurl + '</span> '
	return '<span style="font-size: 12px;">' + eurl + '</span> '
	
}

function formaterBreakPoint(value)
{
	var bp =  value.split(" ");
	var len = bp.length;
	var out;
	
	if (  /BND/.test(value) )
	{	
		out = bp[0];
		for ( i = 1 ; i < len; i++)
		{
				if ( bp[i] != "-")
				{	
					out = out +" "+ bp[i] ;
				}	
		}
		out = '<a href="#" data-toggle="tooltip" title="' + out + '"><b><u> BND  </u></b></a>';
	}
	else
	{
		out="-";
	}
	
	return out;
}

function formaterLISTE(value)
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

function formaterGENES(value)
{
	var out=" ";
	var outshort;
	
	var names = value.split(",");
	var len = names.length-1;
	
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
	
function formaterDejavu(value)
{
	var values = value.split(";");
	var out;
	
	out = '<a href="#" data-toggle="tooltip" title="' + values[1] + '"> ' + values[0] + '</a>';
	return out;
}

function formaterOmin(value)
{
	var eurl = "-";
	
	 if (value == "yes") {
	 	eurl = '<b><font style="color:red">'+value+'</font></b>';
	 }
	return eurl;
}	
	
function formaterRANK(value)
{
	var eurl;
	 if (value > 5) {
	 	eurl = '<b><font style="color:red">' + value + ' </font></b>';
	 }
	 else {
	 	eurl = '<b><font style="color:grey">'+ value + '</font></b>';
	 }
	 
	return eurl;
	
}	

function formaterScore(value)
{
	var eurl;
	
	var score = value;
	
	 if (score > 3) {
	 	eurl = '<b><label style="color:red;">' + score + ' </label></b>';
	 }
	 else {
	 	eurl = '<label style="color:grey">'+score + '</label>';
	 }
	 
	return eurl;
	
}	

function formaterTransloc(value)
{
	var chrs = value.split("to");
	//var out = '<label>' + chrs[0] + '         </label> <b><label style="font-size: 12px;"> <span class="badge badge-pill badge-success"><span class="glyphicon glyphicon-transfer"></span></span></label></b><label>            ' +  chrs[1]  + '</label>';
	//var out = '<b><span class="badge badge-pill badge-warning" style="font-size: 10px;">' + chrs[0] + '</span> <i class="fa fa-refresh" style="font-size:12px"></i> <span class="badge badge-pill badge-success" style="font-size: 10px;">' + chrs[1] + '</span></b>';
	var out = '<table style="width:100%;"><tr><td style="border: 2px solid darkgreen;width:30px; vertical-align:middle;font-size: 10px;">' + chrs[0] + '</td><td style="width:30px"><i class="fa fa-refresh" style="font-size:14px"></i></td><td style="border: 2px solid red;width:30px; vertical-align:middle;font-size: 10px;">' + chrs[1] + '</td></tr></table>';
	return out;
}

function formaterTransmission(value)
{
	var eurl="-";

	if (/Denovo/.test(value) )
	 {
	 	eurl = '<b><i><label style="font-size:10px;color:red"> denovo  </label></i></b>';
	 }
	 
	if (/mother/.test(value))
	{
	 	eurl = '<b><label style="font-size:1.5em;color:pink"> <span class="glyphicon glyphicon-user"> </label></b>';
	 }
	 
	 if  (/father/.test(value)){
	 	eurl = '<b><label style="font-size:1.5em;color:blue"> <span class="glyphicon glyphicon-user"> </label></b>';
	 }
	 
	 if  (/both/.test(value)){
	 	eurl = '<b><label style="font-size:1.5em;color:pink"> <span class="glyphicon glyphicon-user"> </label><label style="font-size:1.5em;color:blue"> <span class="glyphicon glyphicon-user"> </label></b>';
	 }
	 
	 return eurl;
}	

function formaterdbvarStatus(value)
{
	var eurl=value;
	var out=" ";

	var status = value.split(",");
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

function igv_view_BND()
{

	var parameters = location.search.split("&");
	var par0 = parameters[0].split("=");
	var par1 = parameters[1].split("=");
   	var par2 = parameters[2].split("=");
   	var par3 = parameters[3].split("=");
   	var par4 = parameters[4].split("=");
   	
   	
   	var titre=par0[1];
   	var locusstart=par1[1];
   	var locusend=par2[1];
	var bamfile=par3[1];
	var patient_name=par4[1];
	
	
	titre = patient_name + "  :  " + titre + "  (" + taille + ")";
	document.getElementById('titre').innerHTML=titre;
	
	view_web_igv_bam_simple('div_igv_start', locusstart, bamfile, patient_name);
	view_web_igv_bam_simple('div_igv_end', locusend, bamfile, patient_name);
}

function igv_view_SV()
{

	var parameters = location.search.split("&");
	var par0 = parameters[0].split("=");
	var par1 = parameters[1].split("=");
   	var par2 = parameters[2].split("=");
   	var par3 = parameters[3].split("=");
   	var par4 = parameters[4].split("=");
   	var par5 = parameters[5].split("=");
   	
   	var titre=par0[1];
   	var locusstart=par1[1];
   	var locusend=par2[1];
	var bamfile=par3[1];
	var patient_name=par4[1];
	var taille=par5[1];
	
	if (taille == 1)
	{
		titre = patient_name;
	}
	else
	{
		titre = patient_name + "  :  " + titre + "  (" + taille + ")";
	}
	
	document.getElementById('titre').innerHTML=titre;
	
	view_web_igv_bam_simple('div_igv_start', locusstart, bamfile, patient_name);
	view_web_igv_bam_simple('div_igv_end', locusend, bamfile, patient_name);
}


 
 function launch_web_igv_start_and_end(value) 
{
	var val = value.split("et");
	var titre = val[0];
	
	var v = val[1].split(";");
	var locusstart = v[0];
	var locusend = v[1];
	var bamfile =v[2];
	var patient_name = v[3];
	var taille = 1;
	
	if (v.length > 4) 
	{
		taille = v[4];
	} 
	
	var url_viewer = "Url_SV_Viewer.html?titre=" + titre + "&locusstart=" + locusstart + "&locusend=" + locusend + "&bamfile=" + bamfile + "&patient_name=" + patient_name + "&taille=" + taille;
	var myWindow = window.open(url_viewer,"_blank",""); 
 
 }
 
function formatNumber(num) {
  return num.toString().replace(/(\d)(?=(\d{3})+(?!\d))/g, '$1,')
}



function formaterInfoCaller(value)
{
	var out;
	var val = " ";
	var t = [" "];
	
	if (value == 0)
	{
		etiquette = '<font size = "2" color="red"><span class="glyphicon glyphicon-minus-sign"></span></font>';
		return etiquette;
	}
	else
	{
			var tab = value.split("ou");
			var nbid = tab.length-1;
			
					
			for (i=0; i<nbid; i++)
			{
			
				var tab2 = tab[i].split("et");
				var tab3 = tab2[1].split(" ");
				
				var tabid =  tab3[0].split("_");
				var val1 = formatNumber(tabid[2]);
				var val2 = formatNumber(tabid[3]);
			
				t[i]  = val1 + "-" + val2;
			}
		
			t.sort();

			 for (i=0; i<nbid; i++)
			{
				if ( i == nbid - 1)
				{
						val  = val + t[i];
				}
				else
				{
						val  = val + t[i] + " / ";
				}
			}
			
		etiquette = '<font size = "2" color="green"> <span class="glyphicon glyphicon-ok-sign"></span></font>';
		out = '<a href="#" data-toggle="tooltip" title="' + val + '">' +  etiquette  +'  </a>';
		return out;
	}
}

function formater_infoPATIENTIGV(value)
{
	
		var parameters = location.search.split("&");
		var arg1 = parameters[0].split("=");
  	 	var arg2 = parameters[1].split("=");
		
  		 var projectname = arg1[1].split(",");
  	  	var patients = arg2[1].split(",");
    	
    
		var out = " ";
	
		var tv =  value.split('et');
		var bampath = tv[0];
		var real_value = tv[1];
	
	
			var infos =  real_value.split(" ");
			var len = infos.length;
				
			var w = infos[0].split("_");
			var chr = w[1];
				
			// 2  kb autour du start
			var s1 = parseInt(w[2]);
			var start1 = s1 - 500;
			var start2   = s1 + 500;
			var locusstart = 'chr' + chr  + ':' + start1 + '-' + start2;
				
			// 2  kb autour du end
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
   			var boutonID = 'bouton'+ nb;				
   			var igv_value = locusstart + ';' + locusend + ';' + bampath + infos[1] + '.bam;' + infos[1] + ";" + taille;				
			var etiquette   = '<button  id="' + boutonID + '" onclick="launch_web_igv_start_and_end(\'' + infos[0] + "et" + igv_value + '\')"><font color="blue"><u> View  </u>  </font></button>';
	
	return etiquette;
}	


function formater_BND_IGV(value)
{
		var parameters = location.search.split("&");
		var arg1 = parameters[0].split("=");
  	 	var arg2 = parameters[1].split("=");
		
  		 var projectname = arg1[1];
  	  	var patientname = arg2[1];
    	
    
		var out = " ";
	
		var tv =  value.split('et');
		var bampath = tv[0];
		var real_value = tv[1];
	
	
			var w =  real_value.split("_");
			var chr1 = w[1];
			var chr2 = w[3];
				
			// 2  kb autour du start
			var s1 = parseInt(w[2]);
			var start1 = s1 - 500;
			var start2   = s1 + 500;
			var locusstart = 'chr' + chr1  + ':' + start1 + '-' + start2;
				
			// 2  kb autour du end
			var s2 = parseInt(w[4]);
			var end1 = s2 - 500;
			var end2   = s2 + 500;
			var locusend = 'chr' + chr2  + ':' + end1 + '-' + end2;
   				
   			nb++;
   			var boutonID = 'bouton'+ nb;				
   			var igv_value = locusstart + ';' + locusend + ';' + bampath + ";" + patientname;				
			var etiquette   = '<button  id="' + boutonID + '" onclick="launch_web_igv_start_and_end(\'' + patientname + "et" + igv_value + '\')"><font color="grey"><u> View  </u>  </font></button>';
	
	return etiquette;
}



function formater_infoTransmission(value)
{
	var out = "";
	
	if (value )
	{
		var lignes =  value.split(";");
		var len = lignes.length;
		
	
		for ( i = 0 ; i < len; i++)
		{					
					out = out + lignes[i] ;
		}
	
		if (out.length > 3 ) {
			out = '<a href="#" data-toggle="tooltip" title="' + out + '">  see details ...  </a>';
		}
		
	}
	
	return out;
}	


function formaterdbvarEvent(value)
{
	var eurl = value.replace(/copy_number/ig,"CNV");
	return eurl;
}	

function progressBar() {
  var elem = document.getElementById("myBar");
  var width = 1;
  var id = setInterval(frame, 30);
  function frame() {
    if (width >= 100) {
      clearInterval(id);
    } else {
      width++;
      elem.style.width = width + '%';
    }
  }
}
