/**
 * @author pnitschk
 */

var sizew = 800;
var sizeh = 1100;
var surface_trans;
var surface_pat;
var Tsurface_pat;
var Tsurface_trans;
var decx = 2;
var sizeRect = 15;
var sizeRectTranscript = 5;
var space = 4;
var echx = null;
var cadre = 2000;
var oldpos = null;
var echxMap = 0;
var zoom = 1;
var echx =0;
var genomic_start;
var first = true;
var lasthighlight;
var groupVariations;
groupVariations = new Array();
var groupPatients;
groupPatients = new Array();
var groupTranscripts;
groupTranscripts = new Array();
var tabcolor = new Array();
tabcolor[0] = [214,238,240,1];
tabcolor[1] = [245,250,250,1];
var tabcolor2 = new Array();
tabcolor2[0] = [231, 247, 231];
tabcolor2[1] = [255,255,255,1];
var drawingOK = false;

//dojo.addOnLoad(init); 
var container = null;
var container_position = null;   
var container_transcripts = null;
var container_transcripts_position = null;
var container_patients = null;
var container_patients_position = null;
var map_container = null;
var heightmap = 20;
var surfaceMap = null;
var groupMap = null;
var surfaceText =null;
var textw = 100;
     
function initViewer() {
	container_transcripts = dojo.byId('drawing_transcripts'); 
	container_transcripts_position = dojo.coords(container_transcripts, true);  
	container_patients = dojo.byId('drawing_patients'); 
	container_patients_position = dojo.coords(container_patients, true);  
	if (window.innerWidth > (screen.width / 2)) { sizew = window.innerWidth - 70; }
	else { sizew = screen.width -300; }
    sizeh = 1100;
	//var myurl = url_annot+"?project="+projectName+"&selectGene="+selectGene+"&chromosome="+chromosome_name+"&ids="+ids;
	var myurl = url_annot+"?project="+projectName+"&chromosome="+chromosome_name+"&ids="+ids;
    readDataForDrawing(myurl);
}
     
var url_current_covering;

function initCoverViewer(gene,ids,start,end,chromosome){
	document.getElementById('map').innerHTML = "";
	document.getElementById('drawing_transcripts').innerHTML = "";
	document.getElementById('drawing_patients').innerHTML = "";
 	prev_zoom =1;
   	container_text =  dojo.byId('text');
	if (container_text) { container_text.innerHTML = ''; }
    //container = dojo.byId('drawing'); 
	//container.innerHTML ='';
	container_transcripts = dojo.byId('drawing_transcripts'); 
	container_transcripts_position = dojo.coords(container_transcripts, true);  
	container_patients = dojo.byId('drawing_patients'); 
	container_patients_position = dojo.coords(container_patients, true);  
	dojo.byId('map').innerHTML ='' ;
    //container_position = dojo.coords(container, true);  
	var test = document.getElementById("right").offsetWidth;
	sizew = test -50;
	
	
	
	var string_patient_and="";
	for (i=0;i<gridPatient.rowCount;i++) {
		var item = gridPatient.getItem(i);
		if (!item) continue;
		//alert ('Item: '+item+' -> name: '+ item.name +' -> include val: '+item.include);
		if (item.include == 0) {
			string_patient_and += item.name+"+";
		}
		/*	if (item.include == -1) {
			string_patient_not += item.name+"+";
		}
		if (item.include == 2) {
			string_patient_attic += item.name+"+";
		}
		*/
	}
	var myurl = url_annot+"?project="+projectName;
	if (string_patient_and) { myurl += "&patients_and="+string_patient_and; }
	if (chromosome) { myurl +="&start="+start+"&end="+end+"&chromosome="+chromosome; }
	else if (gene) { myurl += "&gene="+gene; }
	if (ids) { myurl +="&ids="+ids; }
	url_current_covering = myurl;
    readDataForDrawing(myurl);
	//var url_covert = "http://mendel.necker.fr/cgi-bin/readCover.pl?project="+projectName;
	//readDataForCover(url_covert);
}
	
//je recupere les donnees a dessiner 
function readDataForDrawing(myurl){
    dojo.xhrGet({ // &acirc;&#382;&sbquo;
	    // The following URL must match that used to test the server.
	    url: myurl,
	    handleAs: "json",
	    timeout: 80000, // Time in milliseconds
	    // The LOAD function will be called on a successful response.
	    load: function(responseObject, ioArgs,evt){
	        dataRef = responseObject.items;
			groupPatients = new Array();
			groupTranscripts = new Array();
			var nb_patient =0;
			genomic_start = dataRef.gstart;
	    	nb_patient = dataRef.patient.length;
	    	sizeh = ((sizeRect+space)*nb_patient)+500;
	    	sizeh += 50;
	        map_container = dojo.byId('map');   
	        surfaceMap = dojox.gfx.createSurface(map_container, sizew + 'px', (heightmap +1) + 'px');
	        dojo.connect(map_container,"onclick",testClick);
	        dojo.connect(map_container,"onmousemove",handleMouseMove1);
	    	Tsurface_trans = dojox.gfx.createSurface(container_transcripts, sizew + 'px', (sizeRectTranscript+space)*(dataRef.transcript.length)+500+50+'px');  
	    	Tsurface_pat = dojox.gfx.createSurface(container_patients, sizew + 'px', sizeh + 'px');  
	    	drawTranscripts(dataRef);
	    	drawPatients(dataRef);
			drawingOK = true;
	        return;
	    },
	    
	    // The ERROR function will be called in an error case.
	    error: function(response, ioArgs){ // &acirc;&#382;&fnof;
	    	alert("error loading data for viewer");
	        // dojo.byId('wait1').src="/icons/myicons/png-24/16-circle-red.png";
			// dojo.byId('wait1').show();
			console.dir(response);
	        console.error("HTTP status code:  _--> ", ioArgs.xhr.status); // &acirc;&#382;&dagger;
	        return response; // &acirc;&#382;&hellip;
	    }
 	});
}

//function initviewer2() {
//    sizeh = (sizeRect+space)*(dataRef.transcript.length+dataRef.patient.length);
//    sizeh+=50;
//    Tsurface = dojox.gfx.createSurface(container, sizew + 'px', sizeh + 'px');  
// 	//surfaceText =  dojox.gfx.createSurface(container_text, textw + 'px',sizeh+ 'px');   
//    map_container = dojo.byId('map');   
//    surfaceMap = dojox.gfx.createSurface(map_container, sizew + 'px', (heightmap +1) + 'px');
//    dojo.connect(map_container,"onclick",testClick);
//    dojo.connect(map_container,"onmousemove",handleMouseMove1);
//    drawReference(dataRef);
//	drawingOK = true;
//}

var makeText = function(surface_pat, text, font, fill, stroke){
    var t = surface_pat.createText(text);
    if(font)   t.setFont(font);
    if(fill)   t.setFill(fill);
    if(stroke) t.setStroke(stroke);
    return t;
};
    	
var groupscontig = new Array();
var dataRef;

function initTranscripts(h) {
    var xref= 1;
	var lref = dataRef.length *1.0;
	var ngene=-1;
	var gene_name;
	var cpt;
	groupTranscripts = new Array();
	if (dataRef.transcript.length > 50) return;
    for (var i=0; i< dataRef.transcript.length;i++) {
        var tr = dataRef.transcript[i];
        if ((tr.x+tr.length)<xref) { continue; }
        if (tr.x>(xref+lref)) { continue; }
        var z = xref+lref;
        if (tr.gene_name && gene_name != tr.gene_name){
			ngene++;
			gene_name = tr.gene_name;
			cpt =0;
		}
		else { cpt ++; }
        var tx = (tr.x-xref);
        var tl = tr.length ;
        if (tx<0){tl += tx;tx=1/echx;}  
        if (tl >lref) {tl=lref;}
		var tname = tr.name;
		var textname = tr.gene_name+":"+tr.name;
		//if (tr.gene_name) {
		//	textname = tr.gene_name;
		//}
        groupTranscripts[tname] = new transcript(tx,h,tl,textname,tabcolor2[ngene%2]);
      	//groupTranscripts[tname].drawText();
        for (var j in tr.exon) {
            var exon = dataRef.transcript[i].exon[j];
            if ((exon.x+exon.length)<xref){
               // continue;
            }
            if (exon.x>(xref+lref)) {
                //continue;
            }
            var ex = (exon.x-xref);
            var el = exon.length;
            if (ex<0){el += ex;ex=1/echx;}  
            if (el >lref) {el=lref;}
			groupTranscripts[tname].addExon(ex,el,exon.strand,exon.utrstart,exon.utrend);
        }
        h += (sizeRect+space);  
		h+=10;
    }
    return h;
}

var legend_height=10;

function change_coverage (){
	for (z in groupPatients){
		groupPatients[z].change_coverage();
	}
}
       
function initVariations (h){
    var xref= 1.0;
	var lref = dataRef.length*1.0;
	legend_height = h;
   	h +=10;
	if (!dataRef.patient) return;
	
	groupPatients = new Array();
	groupVariations = new Array();
    for (var j in dataRef.patient) {
		// if (dataRef.patient[j].trace == 0) continue;
        var c = tabcolor[j%2];
		var namep = dataRef.patient[j].name;
		if (dataRef.patient[j].type == "and") c = [245,250,0,0.5];
		groupPatients[namep] = new patient1(h,dataRef.length,namep,c);
		groupPatients[namep].drawText();
        for (var i in dataRef.patient[j].trace) {
            var tr =  dataRef.patient[j].trace[i];
          	groupPatients[dataRef.patient[j].name].addTrace(tr.x,tr.y,tr.strand);  
        }
        for (var i in dataRef.patient[j].variation) {
            var va = dataRef.patient[j].variation[i];
            if ((va.x+va.length)<xref){ continue; }
            if (va.x>(xref*1.0+lref*1.0)) { continue; }
            var x = (va.x-xref);
            if (x<0) { xx = 1/echx; }
			if (!groupVariations[va.name]) {
				if (va.length == 1) {
					groupVariations[va.name] = new variation(x, h, "red", va);
				}
				else {
					groupVariations[va.name] = new cnv(x, h, "red", va);
				}
			}
            if (va["consequence!all"] == "coding") {  
                groupVariations[va.name].addPatient(namep,h);
            }
			 else if (va["consequence!all"] == "silent") {  
                groupVariations[va.name].addPatient(namep,h);
				groupVariations[va.name].setColor("grey");
            }
            else if (va["consequence!all"] == "utr") {
                groupVariations[va.name].addPatient(namep,h);
                groupVariations[va.name].setColor("purple");
            }
            else {
                groupVariations[va.name].addPatient(namep,h);
                groupVariations[va.name].setColor("blue");
            }
        }
        h += (sizeRect+space);  
    }
    for (var i in groupVariations) {
        if (groupVariations[i].data.filterscore == "weak") statistics["nbweak"]++;
        else if (groupVariations[i].data.filterscore == "medium") statistics["nbmedium"]++;
        else if (groupVariations[i].data.filterscore == "strong") statistics["nbstrong"]++;
        if (groupVariations[i].isCoding())  statistics["nbcoding"]++;
        else if (groupVariations[i].isIntronic())  statistics["nbintronic"]++;
        else if (groupVariations[i].isUTR())  statistics["nbutr"]++;
        if (groupVariations[i].isEnsembl()) statistics["nbdbsnp"]++;
        if (groupVariations[i].data.validate == 1 )  statistics["nbvalidate"]++;
    }
    for (var i in statistics) { statistics2[i] = statistics[i]; }
    update_variation_number();
    first = false;
    return h;
}

function refreshTranscripts(){
	for (var j in groupTranscripts) {
	 	groupTranscripts[j].drawGraph();
	 }
}

function refreshVariations() {
	 for (var j in groupVariations) {
	 	groupVariations[j].draw();
	 }
}

function refreshPatients() {
	 for (var j in groupPatients) {
	 	groupPatients[j].drawGraph();
	 }
}

function drawEverything(){
    echx = (sizew-decx) / dataRef.length; 
    echx *= zoom;
	surface_pat = Tsurface_pat.createGroup();
	refreshPatients();
	refreshVariations();
	surface_trans = Tsurface_trans.createGroup();
	refreshTranscripts();
}

var nbutr;
var drawMapDone = false

function drawTranscripts(data) {
	var h = 0;
	h=initTranscripts(h);
    echx = (sizew-decx) / dataRef.length; 
    echx *= zoom;
	surface_trans = Tsurface_trans.createGroup();
	refreshTranscripts();
}

function drawPatients(data) {
	var h = 0;
	if (drawMapDone == false) {
		drawMap(1,data.length); 
	    drawIndicator(1,data.length);
	}
	h=initVariations(h);
    echx = (sizew-decx) / dataRef.length; 
    echx *= zoom;
	surface_pat = Tsurface_pat.createGroup();
	refreshPatients();
	refreshVariations();
}
    
function drawReference(data) {  
	var h = 0;
	drawMap(1,data.length); 
    drawIndicator(1,data.length);
	h=initTranscripts(h);
	h=initVariations(h);   
	drawEverything();
}
    
var heightrectRef = 20;

function drawMap(xref,lref) {
    xref = xref *1.0;
    lref = lref * 1.0;
    var h = 0;
    echxMap = (sizew-decx) / lref;
    drawRect(surfaceMap,0,0,lref*echxMap-2,heightmap,[200,200,0,0.1],decx);         
    for (var i in dataRef.contig) {
        var contig = dataRef.contig[i];
        var x = contig.x;
        var l = contig.length;
        rect = drawRect(surfaceMap,x*echxMap,h,l*echxMap,heightmap,"orange",decx);  
    }
}
    
var sizeIndicator;
var indicator=null;

function drawIndicator(xref,lref) {
    indicator = new Indicator(surfaceMap,echxMap,lref);
    indicator.draw(zoom);
    groupMap = surfaceMap.createGroup();
    return;
    sizeIndicator = (lref*echxMap)/zoom;
    groupMap = surfaceMap.createGroup();
    if (zoom== 1) return ;
    groupMap.createRect({ x: xref*echxMap +decx  , y:2, height: heightmap-2, width: (lref*echxMap)/zoom})
    .setFill([0,0,254,0.1])
    .setStroke({ color: "grey", width: 2});
    groupMap.createRect({ x: xref*echxMap +decx  , y: 2, height: heightmap-4, width: ((lref*echxMap)/zoom)-2})  
    .setStroke({ color: "black", width: 2});
}
    
function drawRectVariation(surface, x, y, w, h, color, decx){
    var group = surface.createGroup();
    var rectangle = group.createRect({ x: x+decx, y: y, height: h, width: w+1})
    .setFill(color);
	return rectangle;
}

function drawRect (surface,x,y,w,h,color,decx,noborder){
    var group = surface.createGroup();
    if (!noborder) {
        var rectangle = group.createRect({
            x: x + decx,
            y: y,
            height: h,
            width: w
        }).setFill(color);
        var group2 = surface.createGroup();
        var line2 = group.createRect({x:x+decx,y:y,width:w,height:1}).setFill("black");
        group.createRect({x:x+decx,y:h+y,width:w,height:1}).setFill("black");
        group.createRect({x:x+decx,y:y,width:1,height:h}).setFill("black");
        group.createRect({x:x+decx+w,y:y,width:1,height:h}).setFill("black");
        return group ;
    }
    else {
        var rectangle = group.createRect({
            x: x + decx,
            y: y,
            height: h,
            width: w
        }).setFill(color);
    	return group;
    }
}
	
var previousx = 0;
var move_ok_map = false;
var not_move = true;

function handleMouseMove1(event) {
    if (zoom == 1) return;
    if(!indicator.click) return;
    indicator.handleMouseMove(event);
}

function translate(dx) {
    indicator.translate(dx);
    var pixelbas = (dx /echxMap) * echx;
    surface_trans.applyTransform({dx: pixelbas}); 
    surface_pat.applyTransform({dx: pixelbas}); 
}

function testClick(event) {
    if (zoom == 1) return;
    if (indicator.click == true) {
        indicator.click = false;
        dojo.stopEvent(event);
        return;
    }
    else {
        var x = event.clientX - container_transcripts_position.x;
        x += indicator.sizeIndicator / 2;
        indicator.move(x);
        dojo.stopEvent(event);
    }
}
    
var prev_zoom = 1;
function zoomIn(){
	if (prev_zoom < 7) prev_zoom ++ ;
	testSlider(prev_zoom);
}

function zoomOut(){
	if (prev_zoom > 1) prev_zoom -- ;
	testSlider(prev_zoom);
}

function testSlider(prev_zoom) {
    if (indicator.lastx != 2) { indicator.move(2); }
	surface_pat = Tsurface_pat.createGroup();
	testSlider_pat(prev_zoom);
	surface_trans = Tsurface_trans.createGroup();
	testSlider_trans(prev_zoom);
}

function testSlider_trans(z){
    z = parseInt(z);
    zoom = Math.pow(2,(z-1));
    var olds = indicator.sizeIndicator;
    indicator.redraw(zoom);
    echx = (sizew-decx) / dataRef.length; 
    echx *= zoom;
    Tsurface_trans.remove(surface_trans);
	surface_trans = Tsurface_trans.createGroup();
	refreshTranscripts();
    if (zoom > 1) {
        if (previousx > 1) {        
            var dx = (1 - (previousx + olds / 2)) + indicator.sizeIndicator / 2;
            translate(dx);
            previousx = (previousx + olds / 2) - indicator.sizeIndicator / 2;
        }
    }
    else { previousx = 0; }
}

function testSlider_pat(z){
    z = parseInt(z);
    zoom = Math.pow(2,(z-1));
    var olds = indicator.sizeIndicator;
    indicator.redraw(zoom);
    echx = (sizew-decx) / dataRef.length; 
    echx *= zoom;
    Tsurface_pat.remove(surface_pat);
	surface_pat = Tsurface_pat.createGroup();
	refreshPatients();
    redrawAllVariations();
    if (zoom > 1){
        if (previousx > 1) {        
            var dx = (1 - (previousx + olds / 2)) + indicator.sizeIndicator / 2;
            translate(dx);
            previousx = (previousx + olds / 2) - indicator.sizeIndicator / 2;
        }
    }
    else { previousx = 0; }
}
 
    
var showEnsembl = true; 
var showCoding = true;
var showIntronic = true;
var showUTR = true;
var showWeak = true;
var showMedium = true;
var showStrong = true;

function clickEnsembl(){
    showEnsembl = (!showEnsembl);
    redrawAllVariations();
}

function clickCoding(){
    
    showCoding = (!showCoding);
    redrawAllVariations();
}

function clickIntronic(){
    showIntronic = (!showIntronic);
    redrawAllVariations();
}

function clickUTR(){
    showUTR = (!showUTR);
    redrawAllVariations();
}

var filtervalue = 1;

function redrawAllVariations() {
    init_statistics2();
    for (var i in groupVariations) {
        groupVariations[i].hideorshow();
    }
    update_variation_number();
}   
   
var statistics = new Array();
statistics["nbweak"] =0;
statistics["nbmedium"] =0;
statistics["nbstrong"] =0;
statistics["nbcoding"] =0;
statistics["nbutr"] =0;
statistics["nbdbsnp"] =0;
statistics["nbvalidate"] =0;
statistics["nbintronic"] =0;

var statistics2 = new Array();
statistics2["nbweak"] =0;
statistics2["nbmedium"] =0;
statistics2["nbstrong"] =0;
statistics2["nbcoding"] =0;
statistics2["nbutr"] =0;
statistics2["nbintronic"] =0;
statistics2["nbdbsnp"] =0;
statistics2["nbvalidate"] =0;

function update_variation_number(){
    for (var i in statistics) { 
        var st = statistics[i];
		if (!document.getElementById(i)) return; 
        document.getElementById(i).innerHTML = statistics2[i]+"/"+statistics[i];
    }
}

function init_statistics2(){
        for (var i in statistics) {
            statistics2[i] =0;
        }
    
}
 