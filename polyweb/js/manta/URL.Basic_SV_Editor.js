var dataStore_SV;
var gridSV;

var dataStore_SVCompare;
var gridSVCompare;
var dataStore_SVFile;


var url_vcf_Compare_basic = url_path + "/manta/compare_SV_Vcf_basic_V3.pl";


var url_parserManta_BND_bycouple = url_path + "/manta/parser_VcfManta_BND_bycouple.pl";
var url_parserLumpy_BND_bycouple = url_path + "/manta/parser_VcfLumpy_BND_bycouple.pl";
var url_parserLumpyManta_BND_bycouple = url_path + "/manta/parser_VcfLumpyManta_BND_bycouple.pl";

var url_get_patients = url_path + "/manta/getVcfFiles.pl";


function getVcfCompare_info_basic() {
	
	var callers= document.getElementById('Callers').value;
	var hide_sd = document.getElementById('SDSelect').value; 
	
	var parameters = location.search.split("&");
	var arg1 = parameters[0].split("=");
    var arg2 = parameters[1].split("=");
    var projectname = arg1[1].split(",");
    var filename = arg2[1].split(",");
   
	var length = 0;
	var minPR = 0 ;
	var minSU = 0; 
	var minQS = 0; 
	
	if (callers != "WCONDOR")
	{
		var length = document.getElementById('lengthSelect').value;
	}
	
	if (callers == "MANTA")
	{
		var minPR = document.getElementById('PRSelect').value;
	}
	
	if(callers == "LUMPY")
	{
		var minSU = document.getElementById('SUSelect').value; 
	}
	
	if (callers == "CANVAS")
	{
		var minQS = document.getElementById('QSSelect').value; 
	}
	
	dataStore_SVCompare = new dojo.data.ItemFileReadStore({ url: url_vcf_Compare_basic + "?projectname=" + projectname + "&filename=" + filename + "&length=" + length + "&callers=" + callers +"&minPR=" + minPR + "&minSU=" + minSU + "&minQS=" + minQS  + "&hide_sd=" + hide_sd});
	
	var gridSVCompare = dijit.byId("SVCompare_gridDetails");
	gridSVCompare.setStore(dataStore_SVCompare);
	gridSVCompare.store.fetch({
		onComplete: function(items){		
		}
	});
	gridSVCompare.startup();
	return;
}

function setCallerOptions()
{
	var caller= document.getElementById('Callers').value;
	var out;
	var id = "collapse1";
	
	// pour tous 
	
	out = "<table><tr>";
	
	if (caller == "CANVAS")
	{
		out +=' <td padding: 15px; style="width:120px";><label><h6>  Min Lenght  : </h6></label>' +
						'<select class="form-control" id="lengthSelect" name="lengthSelect" style="width:100px; height:30px; font-size:12px;padding:5px">  ' +
							' <option value=1000>  1kb  </option>' +
							' <option value=5000>  5kb  </option>' +
							' <option value=15000>  15kb  </option>' +
							' <option value=50000>  	50 kb  </option>' +
							' <option value=100000 selected="selected">  100 kb  </option>' +
							' <option value=200000>  200 kb  </option>' +
							 '<option value=300000>  300 kb  </option>' +
							' <option value=400000>  400 kb  </option>' +
							' <option value=500000>  500 kb  </option>' +
							' </select></td>' +
				
		 					' <td padding: 15px; style="width:120px";><label><h6>  Canvas : QS  &nbsp&nbsp</h6></label>' +
							'<select class="form-control" id="QSSelect" name="QSSelect" style="width:60px; height:30px; font-size:12px;padding:5px">' +
								 '<option value="0" >  0  </option>' +
								 '<option value="10" >  10  </option>' +
								 '<option value="20"  selected="selected" > 20  </option>' +
								 '<option value="30" >  30  </option>' +
								 '<option value="40">  40  </option>' +
							'</select></td>';
	}
	
	if(caller == "LUMPY")
	{
		out +=' <td padding: 15px; style="width:120px";><label><h6>  Min Lenght  : </h6></label>' +
						'<select class="form-control" id="lengthSelect" name="lengthSelect" style="width:100px; height:30px; font-size:12px;padding:5px">  ' +
							 '<option value=1000>  1kb  </option>' +
							 '<option value=5000>  5kb  </option>' +
							 ' <option value=15000>  15kb  </option>' +
							' <option value=50000>  	50 kb  </option>' +
							' <option value=100000  selected="selected" >  100 kb  </option>' +
							' <option value=200000>  200 kb  </option>' +
							 '<option value=300000>  300 kb  </option>' +
							 '<option value=400000>  400 kb  </option>' +
							' <option value=500000>  500 kb  </option>' +
						 '</select></td>' +
					 	 		
			 			' <td padding: 15px; style="width:120px";><label ><h6> Lumpy : SU  </h6></label>' +
  						'<select class="form-control" id="SUSelect" name="SUSelect" style="width:60px; height:30px; font-size:12px;padding:5px">' +
   							 '<option value="0" >  0 </option>' +
							'<option value="20"  selected="selected" >  20  </option>' +
							'<option value="30"  >  30  </option>' +
							'<option value="40" >  40  </option>' +
							'<option value="50" >  50  </option>' +
  						'</select></td>';
	}
	
	if(caller == "MANTA")
	{
		out +=' <td padding: 15px; style="width:120px";><label><h6>  Min Lenght  : </h6></label>' +
						'<select class="form-control" id="lengthSelect" name="lengthSelect" style="width:100px; height:30px; font-size:12px;padding:5px"' +
							' <option value=1000>  1kb  </option>' +
							' <option value=5000>  5kb  </option>' +
							' <option value=15000>  15kb  </option>' +
							' <option value=50000>  	50 kb  </option>' +
							' <option value=100000  selected="selected" >  100 kb  </option>' +
							' <option value=200000>  200 kb  </option>' +
							' <option value=300000>  300 kb  </option>' +
							' <option value=400000>  400 kb  </option>' +
							' <option value=500000>  500 kb  </option>' +
							' </select></td>' +
					
			 			 '<td padding: 15px; style="width:120px";><label "><h6>  Manta : PR %  </h6></label>' +
  						'<select class="form-control" id="PRSelect" name="PRSelect" style="width:60px; height:30px; font-size:12px;padding:5px">' +
   							'<option value="0"   selected="selected" >  0 </option>' +
							'<option value="20"  selected="selected" >  20  </option>' +
							'<option value="30">  30  </option>' +
							'<option value="40" >  40 </option>' +
							'<option value="50">  50  </option>' +
  						'</select></td>';
	}
	
	// pour tous
	
	out += '<td padding: 15px; style="width:250px";><label><h6>  Segmental Duplication  : </h6></label>' +
						'<select class="form-control"  id="SDSelect" name="SDSelect" style="width:100px; height:30px; font-size:12px;padding:5px" onchange="collapseResult();getVcfCompare_info_basic();">' +
							 '<option value="0" >  Shown  </option>' +
							 '<option value="1"  selected="selected" >  Hidden  </option>' +
						'</select><td>' +
				
						'<td padding: 15px; style="width:250px";><label><h6>  Search for SV  </h6></label>' +
						'<div><button type="button" class="btn btn-warning" onclick="collapseResult();getVcfCompare_info_basic();" >  RUN  </button></div></td>' + 
	'</tr></table>';
   	
   		document.getElementById('caller_option').innerHTML = out;
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

function getProjects(){
	
	var projects = ["NGS2018_2300","NGS2018_2047","NGS2018_2291","Strasbourg"];
	var out;
	
	var len = projects.length;
	
	for ( i = 0 ; i < len; i++)
	{
		out += "<option value=" + projects[i]  + ">" + projects[i] + "</option>";
	}
	
	document.getElementById('VcfProjectSelect').innerHTML = out;
}

function getVcfNames(){
	
	var project = document.getElementById('VcfProjectSelect').value;
	var NameStore = new dojo.data.ItemFileReadStore({ url: url_get_patients + "?project=" + project + "&trios=" + 0 });
	var out;
	
	var gotList = function(items, request)
	{
 	 		var itemsList = "";
  			dojo.forEach(items, function(i){
    			itemsList += NameStore.getValue(i, "name") + ",";
 			 });
 	 
  			var tabNames = itemsList.split(",");
			var len =tabNames.length;
	
			for ( i = 0 ; i < len; i++)
			{
				out += "<option value=" + tabNames[i]  + ">" + tabNames[i] + "</option>";
			}
			document.getElementById('VcfFileSelect').innerHTML = out;
	}
	
	
	var gotError = function(error, request){
 	 	alert("The request to the store failed. " +  error);
	}
	
	// Invoke the search 
	
	NameStore.fetch({
  		onComplete: gotList,
  		onError: gotError
	});
	
}
	

function formaterDGV(value)
{
	var eurl = '<b><u><a href="' + value  + '" target="_blank"> DGV </a></u></b>';
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
