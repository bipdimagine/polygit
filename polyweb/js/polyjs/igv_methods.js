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
	 return 1;
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

function getUrlLoadingGenome(file) {
	var url = '';
	var re_hg19 = new RegExp('HG19*');
	var re_hg38 = new RegExp('HG38*');
	var re_mm38 = new RegExp('MM38*');
	var re_mm39 = new RegExp('MM39*');
	
//	if (url_fasta) {
//		url += "&genome=" + window.location.origin + '/' + url_fasta;
//		if (url_gene_bed) {
//			url += "&file=" + window.location.origin + '/' + url_gene_bed;
//		}
//	}
//	else if (re_hg19.test(file)) {

	if (re_hg19.test(file)) {
	    url += "&genome=hg19";
	}
	else if (re_hg38.test(file)) {
	    url += "&genome=hg38_1kg";
	}
	else if (re_mm38.test(file)) {
	    url += "&genome=mm38";
	}
	else if (re_mm39.test(file)) {
	    url += "&genome=mm39";
	}
	return url;
}

function displayOneBAMIGV(file){
    var url = "http://localhost:60151/load?file="+file;
    url += getUrlLoadingGenome(file);
	url += "&merge=true";
    httpGetFocusOn(url);
	return 1;
}

function addOneBAMIGV(file){
    var url = "http://localhost:60151/load?file="+file;
	url += "&merge=true";
    httpGetFocusOn(url);
	return 1;
}

var old_tracks;
function init_igv (){
	var url = "http://localhost:60151/load";
	var polyweb_url =window.location.origin; // url_fasta = "/public-data/"+genome+"/genome/fasta/all.fa";
	 var url_fasta2 = polyweb_url+"/public-data/genome/"+igv_genome+"/fasta/all.fa";
		url += "?genome=hg19"+"&file="+polyweb_url+url_gene_bed;
		url += "&merge=false";
	    httpGetFocusOn(url);
	return 1;		
}

function init_igv_hg38 (){
	var url = "http://localhost:60151/load";
	url += "?genome=hg38";
	url += "&merge=true";
    httpGetFocusOn(url);
	return 1;		
}

function igv_reset () {
	previous_bams = new Array();
}

var rna_last_tracks = '';
function launch_igv_tool_rna(fasta, tracks, locus) {
	var url = "http://localhost:60151/load";
	url += "?file=" + tracks;
	url += "&locus=" + locus;
	if (fasta) { url += "&genome=" + fasta; }
	if (tracks == rna_last_tracks) { url += "&merge=true"; }
	else { url += "&merge=false"; }
	rna_last_tracks = tracks;
    httpGetFocusOn(url);
	return 1;
}


