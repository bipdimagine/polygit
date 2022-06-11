var dataStore_SV;
var gridSV;

var dataStore_SVCompare;
var gridSVCompare;
var dataStore_SVFile;


var url_vcf_Compare_trio = url_path + "/manta/compare_SV_Vcf_trio_V3.pl";

function setdefault()
{
	var parameters = location.search.split("&");
	var arg1 = parameters[0].split("=");
    var arg2 = parameters[1].split("=");
    var arg3 = parameters[3].split("=");
    var arg4 = parameters[4].split("=");
    var arg5 = parameters[5].split("=");
    var arg6 = parameters[6].split("=");
   
   var projectname = arg1[1].split(",");
   var filename = arg2[1].split(",");
   var caller = arg3[1].split(",");
   var length = arg4[1].split(",");
   var hide_sd = arg5[1].split(",");
   
   document.getElementById('Callers').value = caller;
   document.getElementById('SDSelect').value = hide_sd;
   document.getElementById('lengthSelect').value = length;
  
}


function getVcfCompare_info_trio() {
	var callers= document.getElementById('Callers').value;
	var hide_sd = document.getElementById('SDSelect').value; 
	var length = document.getElementById('lengthSelect').value;
	
	var parameters = location.search.split("&");
	var arg1 = parameters[0].split("=");
    var arg2 = parameters[1].split("=");
   
   var projectname = arg1[1].split(",");
   var filename = arg2[1].split(",");
	
	var Selectransmission = document.getElementById('transmissionSelect');
	var transmission = "";
  	var i;
 	 for (i = 0; i < Selectransmission.length ;i++) {
 	 		if (Selectransmission.elements[i].checked)
 	 		{
    				transmission += Selectransmission.elements[i].value + " ";
    			}
  }
	
	
	dataStore_SVCompare = new dojo.data.ItemFileReadStore({ url: url_vcf_Compare_trio + "?projectname=" + projectname + "&filename=" + filename + "&length=" + length + "&callers=" + callers + "&transmission=" +transmission+ "&hide_sd=" + hide_sd});
	
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

function formaterDGV(value)
{
	var eurl = '<b><u><a href="' + value  + '" target="_blank"> DGV </a></u></b>';
	if (value == "-") 
	{
		eurl = '<b> - </b>';
	}
	return eurl;
}
	
	
function formaterDUPDEL(value)
{
	var eurl;
	
	 if (value == "DEL") {
	 	eurl = '<b><p style="color:magenta"> DEL </p></b>';
	 }
	 
	if (value == "DUP") {
	 	eurl = '<b><p style="color:blue"> DUP </p></b>';
	 }
	 
	 if (value == "INV") {
	 	eurl = '<b><p style="color:black">INV </p></b>';
	 }
	 
	 if (value == "INS") {
	 	eurl = '<b><p style="color:black"> INS </p></b>';
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
	return value;	
}	


function formaterFreq(value)
{
	var eurl;
	
	if(value != "-")
	{
		value = value*100;
		eurl = value.toPrecision(2);
		
	
	 		if (value < 0.05)
			{
	 			eurl = '<b><p style="color:red">' + eurl + ' % </p></b>';
	 		}
	 		else
	 		{
	 			if (value < 0.1)
				{
	 				eurl = '<b><p style="color:orange">' + eurl + ' % </p></b>';
	 			}
	 			else
	 			{
	 				eurl = '<b><p>' + eurl + ' % </p></b>';
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
	 			eurl =  '<p>' + value.toPrecision(3) + ' %  </p>';
	 }

	return eurl;
}

function formaterBreakPoint(value)
{
	var bp =  value.split("next");
	var len = bp.length;
	var out = " ";
	
	for ( i = 0 ; i < len; i++)
	{
			out = out  + bp[i] + '<br>';

	}
	return out;
	
}

function formaterLISTE(value)
{
	var out=" ";
	var outshort;
	
	var names = value.split(",");
	var len = names.length;
	
	if (value == "- -")
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
	
	outshort = names[0] + " " + names[1];
	
	if (len < 5)
	{
			out = '<p style="color:black">'+ out + '</p>';
			return out;
	}
	else
	{
			out = '<a href="#" data-toggle="tooltip" title="' + out + '"> ' + outshort + '  ...  see more  </a>';
			return out;
	}
}		

function formaterCallers(value)
{
	var out=" ";
	var outshort=" ";
	
	if (  /lumpy/.test(value) )
	{
		outshort = outshort + "LUMPY ";
	}
	if (  /canvas/.test(value) )
	{
		outshort = outshort + "CANVAS ";
	}
	if (  /manta/.test(value) )
	{
		outshort = outshort + "MANTA ";
	}
	if (  /wcondor/.test(value) )
	{
		outshort = outshort + "WCONDOR ";
	}
	
	
	var names = value.split(",");
	var len = names.length;
	
	for ( i = 0 ; i < len; i++)
	{
		out = out + "  " + names[i] ;
	}
	
	out = '<a href="#" data-toggle="tooltip" title="' + out + '"> ' + outshort + ' </a>';
	return out;
	
}		

function formaterOmin(value)
{
	var eurl = "-";
	
	 if (value == "yes") {
	 	eurl = '<b><p style="color:red">'+value+'</p></b>';
	 }
	return eurl;
}	
	
function formaterRANK(value)
{
	var eurl=value;
	
	 if (value < 3) {
	 	eurl = '<b><p style="color:green">'+value+' </p></b>';
	 }
	 if (value > 3) {
	 	eurl = '<b><p style="color:red">'+value+' </p></b>';
	 }
	return eurl;
	
}	

function formaterScoreT(value)
{
	var eurl=value;
	
	 if (/denovo/.test(value)) {
	 	eurl = '<b><p style="color:red">'+value+' </p></b>';
	 }
	 if (/recessive/.test(value)) {
	 	eurl = '<b><p style="color:green">'+value+' </p></b>';
	 }
	  if (/mother/.test(value)) {
	 	eurl = '<b><p style="color:pink">'+value+' </p></b>';
	 }
	  if (/father/.test(value)) {
	 	eurl = '<b><p style="color:blue">'+value+' </p></b>';
	 }
	  if (/error/.test(value)) {
	 	eurl = '<b><p style="color:black">'+value+' </p></b>';
	 }
	return eurl;
	
}	

function formaterdbvarStatus(value)
{
	var eurl=value;
	
	 if ( /athogenic/.test(value) ) {
	 	eurl = '<p style="color:red">' + value + ' </p>';
	 }
	return eurl;
	
}	


function formaterdbvarEvent(value)
{
	var eurl = value.replace(/copy_number/ig,"CNV");
	eurl = '<p>' + eurl + ' </p>';
	
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
