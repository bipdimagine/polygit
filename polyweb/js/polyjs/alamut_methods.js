function httpGetFocusOn(url) {
	var xmlHttp = null;
	xmlHttp = new XMLHttpRequest();
	if (xmlHttp) {
		xmlHttp.open( "GET",url, true );
		xmlHttp.send(null);
	}
}


function displayInAlamut(chr, start, a){
	alert("coucou");
	chr = chr.replace('-', '');
	var locus;
	if (a[0] == '-') { locus = chr+":"+start; }
	else if (a[1] == '-') { locus = chr+":"+start; }
	else { locus = chr+":"+start+a[0]+">"+a[1]; }
	var list_cmds = [];
	list_cmds.push("http://localhost:10000/search?request="+locus);
	launch_alamut_visual_plus_cmds(list_cmds);
	return;
}

function displayLocusInAlamut(chr, start, end){
	chr = chr.replace('-', '');
	var locus = chr + ":" + start + "-" + end;
	var list_cmds = [];
	list_cmds.push("http://localhost:10000/search?request="+locus);
	launch_alamut_visual_plus_cmds(list_cmds);
	return;
}

function displayInAlamutVector(a){
	console.log("http://localhost:10000/show?request="+a);
	var list_cmds = [];
	list_cmds.push("http://localhost:10000/search?request="+a);
	launch_alamut_visual_plus_cmds(list_cmds);
	return;
}

function httpGetLoadBam(){
	var xmlHttp = null; 
	xmlHttp = new XMLHttpRequest();
	var bamXmlHttp = null;
	bamXmlHttp = new XMLHttpRequest();
	
	var item = dijit.byId("gridGenes").getItem(0);
	var gene = dijit.byId("gridGenes").store.getValue(item,"label");
	
	var patient = return_selected_patient(tab_editor);
	var gridPatients = dijit.byId("gridPatients");
	var alamut_api_key = get_alamut_api_key();
	
	for(var i = 0; i < gridPatients.rowCount; i++){
	 	var item2 = gridPatients.getItem(i);
	    var p = gridPatients.store.getValues(item2, "id");
		if(p==patient){
			var url_bam=dijit.byId("gridPatients").store.getValue(item2,"bam");
			if (xmlHttp){
				if (alamut_api_key) {
					xmlHttp.open("GET", "http://localhost:10000/search?apikey="+alamut_api_key+"&request="+gene+"&synchronous=true", true);
				}
				xmlHttp.send(null);
				xmlHttp.onreadystatechange=function(){
					if (xmlHttp.readyState==4){
						if (alamut_api_key) {
							bamXmlHttp.open("GET", "http://localhost:10000/open?apikey="+alamut_api_key+"&filetype=BAM&path=https://polyweb.fr/"+url_bam, true);
						}
						bamXmlHttp.send(null);
					}
				}
			}
		}
	}
}

function httpGetLoadListBam(string_list_bam){
	var list_bam = String(string_list_bam).split(',');
	var item = dijit.byId("gridGenes").getItem(0);
	var gene = dijit.byId("gridGenes").store.getValue(item,"label");
	var list_cmds = [];
	list_cmds.push("http://localhost:10000/search?request="+gene+"&synchronous=true");
	for(var i = 0; i < list_bam.length; i++){
		list_cmds.push("http://localhost:10000/open?filetype=BAM&path=https://www.polyweb.fr/"+list_bam[i]);
	}
	launch_alamut_visual_plus_cmds(list_cmds);
}

function httpGetLoadOnlyListBam(string_list_bam){
	var list_bam = String(string_list_bam).split(',');
	var list_cmds = [];
	list_cmds.push("http://localhost:10000/search?request=BRCA1&synchronous=true");
	for(var i = 0; i < list_bam.length; i++){
		list_cmds.push("http://localhost:10000/open?filetype=BAM&path=https://www.polyweb.fr/"+list_bam[i]);
	}
	launch_alamut_visual_plus_cmds(list_cmds);
}

function httpGetLoadOnlyListBamInGene(gene_name, string_list_bam){
	var list_bam = String(string_list_bam).split(',');
	var list_cmds = [];
	list_cmds.push("http://localhost:10000/search?request="+gene_name+"&synchronous=true");
	for(var i = 0; i < list_bam.length; i++){
		list_cmds.push("http://localhost:10000/open?filetype=BAM&path=https://www.polyweb.fr/"+list_bam[i]);
	}
	launch_alamut_visual_plus_cmds(list_cmds);
}

function launch_alamut_visual_plus_cmds(list_cmds) {
	dijit.byId('waiting').show();
	var user_name = dojo.cookie("username");
	var url = url_path + "json_output_nodb/get_alamut_api_key.pl";
	url += "?user_name=" + user_name;
	$.getJSON( url, function( data ) {
		dijit.byId('waiting').hide();
		if (data.api_key) {
			var licence = data.licence_alamut;
			var api_key = data.api_key;
			var institution_key = data.institution_key;
			console.log("LOGIN: "+user_name+" - ALAMUT_LICENCE: "+licence+" - API_KEY: "+api_key+" - INSTITUTION: "+institution_key);
			var xmlHttp2 = null;
			xmlHttp2 = new XMLHttpRequest();
			xmlHttp2.open("GET", list_cmds[0]+"&apikey="+api_key+"&institution="+institution_key, true);
			xmlHttp2.send(null);
			xmlHttp2.onreadystatechange=function(){
				if (xmlHttp2.readyState==4){
					for(var i = 1; i < list_cmds.length; i++){
						var url_request = list_cmds[i]+"&apikey="+api_key+"&institution="+institution_key;
						var bamXmlHttp2 = null;
						bamXmlHttp2 = new XMLHttpRequest();
						bamXmlHttp2.open("GET", url_request, true);
						bamXmlHttp2.send(null);
					}
				}
			}
		}
    });
	return;
}

