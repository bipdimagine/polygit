/**
 * @author pnitschk
 */


function jumpToIgvXml(items) {
	var item = items[0];
//	var url_igv2 ="http://www.broadinstitute.org/igv/dynsession/igv.jnlp?sessionURL="+item.url;
	var url_igv2 ="http://www.broadinstitute.org/igv/projects/current/igv.php?sessionURL="+item.url;
	//var url_igv2 ="http://data.broadinstitute.org/igv/projects/2.4/igv24_hm.jnlp?sessionURL="+item.url;
	var z = item.url;
	//console.log(z+" "+url_igv2);
	window.open(url_igv2);
}
  
function runIGV (){
	var item = items[0];
//	var url_igv2 ="http://www.broadinstitute.org/igv/dynsession/igv.jnlp?sessionURL="+item.url;
	var url_igv2 ="http://www.broadinstitute.org/igv/projects/current/igv.php?sessionURL="+item.url;
	var url_igv2 ="http://data.broadinstitute.org/igv/projects/2.4/igv24_hm.jnlp";

	//console.log(url_igv2);
	window.open(url_igv2);
}
function httpGetFocusOn(url) {
	var xmlHttp = null;
	xmlHttp = new XMLHttpRequest();
	if (xmlHttp) {
		xmlHttp.open( "GET",url, true );
		xmlHttp.send(null);
	}
}

function LoadIGV (){
	
	var url_igvbis =url_igv+"?project="+projectName;
	var store =	new dojo.data.ItemFileWriteStore({url:url_igvbis});
	
	store.fetch({
				onComplete: jumpToIgvXml
			
			});	
	
}
function LoadIGVPatients () {
	
	
	var store =	new dojo.data.ItemFileWriteStore({url:url_igvbis});
	var sel = gridPatient.selection.getSelected();
	var patient_url="&patients=";
	var l = sel.length;
	
	var allids='';
	
	//pour chaque ligne selectionnee
	for (var i = 0; i < l; i++) {
		patient_url +=sel[i]["name"]+",";
		
	}
	
	var url_igvbis =url_igv+"?project="+projectName+patient_url;

	var store =	new dojo.data.ItemFileWriteStore({url:url_igvbis});
	
	store.fetch({
				onComplete: jumpToIgvXml
			
			});
}

function LoadIGVPatients_diag () {
	

	

	var url_igvbis =url_igv+"?project="+projectName;
	

	
var tsamples = dijit.byId("gridPatients").selection.getSelected();
        var p="";
        if (tsamples.length == 0) {
            p="all";
        }
        else {
            
            for (i=0;i<tsamples.length;i++){
                p += tsamples[i].label+",";
            }
        }
	if (p != "all") {
		url_igvbis += "&patients=" + p;
	}
	
	var store =	new dojo.data.ItemFileWriteStore({url:url_igvbis});

	store.fetch({
				onComplete: jumpToIgvXml
			
			});
}

function displayInIGVLocus(locus) {
	 var url = "http://localhost:60151/goto?locus="+locus;
	 httpGetFocusOn(url);
	 return;
	 socket = dojox.socket({ 
         url: "localhost:60151", 
         headers: { 
             "Accept": "application/json", 
             "Content-Type": "text/plain" 
         } 
     }); 
     if (!socket) 
     { 
         console.log("Socket error"); 
         return; 
     } 
     console.log("socket ok"); 
}

function displayInIGV(chr, start, end){
	chr = String(chr).replace('chrchr', 'chr');
	var locus = chr+":"+start+"-"+end;
	displayInIGVLocus(locus);
}


// function httpGetLoadBam(gene,bamFile) {
	// var xmlHttp = null;
	// xmlHttp = new XMLHttpRequest();
	// var bamXmlHttp = null;
	// bamXmlHttp = new XMLHttpRequest();
// 
	// if (xmlHttp) {
		// xmlHttp.open( "GET", "http://localhost:10000/show?request="+gene+"&synchronous=true", true );
		// xmlHttp.send(null);
		// xmlHttp.onreadystatechange=function() {
			// if (xmlHttp.readyState==4) {
				// bamXmlHttp.open( "GET", "http://localhost:10000/show?request=BAM%3C"+bamFile, true );
				// bamXmlHttp.send(null);
			// }
		// }
	// }
// }

function displayOneBAMIGV(file){
    var url = "http://localhost:60151/load?file="+file;
    console.dir(url);
    httpGetFocusOn(url);
	return 1;
}

var old_tracks;
function launch_igv_tool(fasta, tracks, locus) {
	var url = "http://localhost:60151/load";
	if (old_tracks == tracks){
		return;
	} 
	old_tracks = tracks;
	if (fasta) {
		url += "?genome=" + fasta;
		url += "&file=" + tracks;
		url += "&locus=" + locus;
		url += "&merge=false";
	}
	else {
		url += "?file=" + tracks;
		url += "&locus=" + locus;
		url += "&merge=false";
	}
    console.dir(url);
    httpGetFocusOn(url);
	return 1;
}


