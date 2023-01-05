/**
 * @author pnitschk
 */

//utiliser pour tracer les electrophoregrammes

//chargement de quelques modules dojo necessaires
dojo.require("dojox.charting.Chart2D");
dojo.require("dojox.charting.Chart3D");
dojo.require("dojox.charting.plot3d.Base");
dojo.require("dojox.charting.plot3d.Bars");
dojo.require("dojo.data.ItemFileWriteStore");


//fonction appelee quand aucun electrophoregramme n'existe
function drawBlank(div){
    var surface2 = dojox.gfx.createSurface(div, 250 + 'px', 180 + 'px');
    var aShape = surface2.createImage({
        width: 100,
        height: 100,
        src: "/icons/Polyicons/null.png"
    }).applyTransform({
        dx: 75,
        dy: 40
    });
    
    
}
var listError = new Array();
var restart = 0;
var nb_connection =0;
var max_connection;
//fonction pour recuperer les donnees permettant de tracer les electrophoregrammes (cf cgi chart_json.pl)
var nb_drawing =0;
function readTraceRawDatasArray(listConnections) {
	
	var args = listConnections.pop();
	if (args["error"]>2){
		return;
	}
	// liste de arguments a passer au cgi
	var url_arg = "?";
	for ( var i in args )
	{
		if (i == "div") continue;
		url_arg +=i+"="+args[i]+"&";
    	
	} 
	
	//var url_arg = "?name=" +  args["trace"]+"&projectName=" + args["project"] + "&contigName=" +args["contig"]  + "&variationName=" + args["variation"]+"&variation_id="+ args["variation_id"]+"&trace_id="+args["trace_id"];
	
	 dojo.xhrGet({ //
        // The following URL must match that used to test the server.
        url: url_chart + url_arg,
        handleAs: "json",
        timeout: 100000, // Time in milliseconds
        // The LOAD function will be called on a successful response.
        load: function(responseObject, ioArgs){
            // Now you can just use the object
			//pour les projets de sequencage classique on trace des chromatogrammes
			
		
			if (responseObject[0].type == 1) {
				TraceCharts(responseObject, args["div"], listConnections);
			}
			//pour les projet de puce de resequencage on trace des histogrammes
		
            else {
				TraceHisto(responseObject, args["div"], listConnections);
			}
        },
        
        // The ERROR function will be called in an error case.
        error: function(response, ioArgs){ // 
    	alert("error");
		//en cas d'erreur
        if(listConnections.length>0){
			 args["error"] ++;
			//listConnections.push(args);
			readTraceRawDatasArray(listConnections);
			
		}
          	args["div"].innerHTML='Loading Error ...';
			console.dir(ioArgs);
			console.error(url_chart+url_arg); // 
            console.error("HTTP status code:  _--> ", ioArgs.xhr.status); // 
            return response; // 
        }
    });
}

// trace for cnv

function TraceRawCnvs(div,array,array2) {
	
	  div.innerHTML=' ';
	  div.style.height= 150+"px";
	  div.style.width= 300+"px";
		
		var chart1 = new dojox.charting.Chart2D(div);
		
		chart1.addPlot("default", {type: "Scatter"});
		chart1.addPlot("grid", {type: "Grid"});
		chart1.addAxis("x");
		chart1.addAxis("y", {vertical: true, fixUpper: "major", includeZero: true});
		chart1.addSeries("Series A", array, {
			stroke: {
				color: "red",width: 0.001
			}
		});
		chart1.addSeries("Series B", array2, {
			stroke: {
				color: "green",width: 0.001
			}
		});
		
		
		chart1.render()

}


var traceHeight = 150;

//fonction pour tracer les histogrammes (puce de resequencage)
function TraceHisto(arrayall, div,list) {
	
	  var  array = arrayall[0].data;
	 
	  var st1 = new Array();
	  var st2 = new Array();
	  var st3 =  new Array();
	 
	  var l = array.length;
		
       for (var i = 0; i < 4; i++) {
	   	st1[i] =array[i]*1.00001;
		st2[(3-i)] =array[i+4]*1.00001;
	
	   }
	
	  div.innerHTML=' ';
	  div.style.height= 150+"px";
	  div.style.width= 300+"px";
	
	  chart1 = new dojox.charting.Chart2D(div);
	  
     chart1.addAxis("y", {vertical: true});
     chart1.addAxis("x", { 
	 labels: [{value: 1, text: "A"}, {value: 2, text: "C"}, 
		{value: 3, text: "G"}, {value: 4, text: "T"}
		]
	 });
	
	 
    chart1.addPlot("default", {
      hAxis: "x",
      vAxis: "y",
      type: "ClusteredColumns",
	  ticks: true,
	  precision: 1,
	  labelOffset: 10,
      areas: true,
	  fixed: true,
      gap: 4,
    });
	
      chart1.addSeries("Series 1", st1,{stroke: {color: "black"}, fill: "green"});
	   chart1.addSeries("Series 2", st2,{stroke: {color: "black"}, fill: "red"});
	    chart1.render();
}

function click_surface (vid,tid){
	alert("click "+tid);
	
} 
//fonction pour tracer les electrophoregrammes (sequencage classique);
function TraceCharts(arrayall, div,list){
		
		
    array = arrayall[0].data;
    var start = 10;
	div.innerHTML=' ';
	var sizew = 210;
	//traceHeight = 75;
    var surface2 = dojox.gfx.createSurface(div, sizew + 'px', traceHeight + 'px');

	var myflag = "/icons/Polyicons/0169_star2.png";

	 surface2.createRect({
            x: 0,
            y: 0,
            height: traceHeight,
            width: sizew
        }).setFill([255, 255, 255]);
		
    if (arrayall[0].data.found == 1) {
		myflag = "/icons/Polyicons/0188_pink_star.png";
   
    }
	
    var surface = surface2.createGroup();
 
    var img = "/icons/Polyicons/16-circle-green-add.png";
    var dx = 220;
    var dy = 0;
   
  

    var polyline2 = surface.createPolyline(array.C).setFill(null).setStroke({
        color: [0, 0, 0, 0.1],
        width: 1
    }).applyTransform({
        dy: start + 2,
        dx: 1
    });
	
    var polyline2 = surface.createPolyline(array.C).setFill(null).setStroke({
        color: "blue",
        width: 1
    }).applyTransform({
        dy: start
    });
    
    var polyline3 = surface.createPolyline(array.G).setFill(null).setStroke({
        color: "black",
        width: 1
    }).applyTransform({
        dy: start
    });
    ;
    var polyline4 = surface.createPolyline(array.T).setFill(null).setStroke({
        color: "red",
        width: 1
    }).applyTransform({
        dy: start
    });
    ;
    var polyline5 = surface.createPolyline(array.A).setFill(null).setStroke({
        color: "green",
        width: 1
    }).applyTransform({
        dy: start
    });
    ;
    var polyline6 = surface.createLine(array.toto).setFill(null).setStroke({
        color: "black",
        width: 0.5
    }).applyTransform({
        dy: start
    });
    var color = "orange";
    var trace_pos = array.name;
	var text_strand = "reverse";
	var color_strand = "red";
	if (array.strand == 1){
	 text_strand = "forward";
	 color_strand = "#10f810";
	} 
    //ecriture de la variable text sur le graph
    t1 = makeText(surface, {
        x: 75,
        y: (traceHeight-10) + start,
        text: text_strand,
        align: "start"
    }, {
        family: "Verdana",
        size: "8pt",
        weight: "bold"
    }, "black");
    
    t1 = makeText(surface, {
        x: 75,
        y: (traceHeight-11) + start,
        text: text_strand,
        align: "start"
    }, {
        family: "Verdana",
        size: "8pt",
        weight: "bold"
    }, color_strand);
	
	  t2 = makeText(surface, {
        x: 90,
        y: (traceHeight-20) + start,
        text: arrayall[0].trace_pos,
        align: "start"
    }, {
        family: "Verdana",
        size: "6pt",
        weight: "bold"
    }, "black");
    var bcolor = new Array();
    bcolor["A"] = "green";
    bcolor["T"] = "red";
    bcolor["C"] = "blue";
    bcolor["G"] = "black";
	bcolor["-"] = "black";
	
	//ecriture base sur le graphe

    for (var i in arrayall[0].data.bases) {
        var base = arrayall[0].data.bases[i];
		
		var color = bcolor[base.lettreRef];
 		
		var font = "8pt";
		var dec = 2;
		var decy= 0;
			var posx_ref = base.electro_pos;
		
		if (base.lettreRef.length > 1) {
			color = bcolor[base.lettre];
			 if (array.strand == -1) {
			 	dec = 0;
			
				var n = (i*1.0+1);
			
				posx_ref = arrayall[0].data.bases[n].electro_pos+5;
				
			 }
			 else {
			 //alert("cocuou");
			 //dec = -5;
			 }
			if (base.lettreRef.length > 2) {
				decy = -10;
				font = "6pt";
			}
			
		}
		
        var t3 = makeText(surface, {
            x: base.electro_pos - 2,
            y: start + 26,
            text: base.lettre,
            align: "start"
        }, {
            family: "Verdana",
            size: "8pt",
        },bcolor[base.lettre]);
	
		
        var t4 = makeText(surface, {
            x: posx_ref ,
            y: start + 10 +decy,
            text: base.lettreRef,
            align: "start"
        }, {
            family: "Verdana",
            size: font,
        }, color);
		
        var img2 = "/icons/Polyicons/bullet_green.png"
        if (base.lettre != base.lettreRef) {
            img2 = "/icons/Polyicons/bullet_red.png"
        }
		if (base.lettreRef.length > 1) {
			
			img2 ="/icons/Polyicons/cross.png";
		}
        var aShape = surface.createImage({
            width: 10,
            height: 10,
            src: img2
        }).applyTransform({
            dx: base.electro_pos - 2,
            dy: 19
        });
        
        if (base.variation == 1) {
        
            var aShape = surface.createImage({
                width: 12,
                height: 12,
                src: myflag
            }).applyTransform({
                dx: base.electro_pos - 3
            });
        }
		
		
    }
	
if (progress_drawing != 0) {
				progress_drawing.update({
					progress: nb_drawing++
				});
				if (nb_drawing >= max_connection) {
					waiting_drawing.hide();
					nb_drawing=0;
					max_connection=0;
				}
			}
	if(list && list.length>0){
		
		delete arrayall;
		readTraceRawDatasArray(list);
	}
	delete arrayall; //arrayall = new Array();
}

var old_tr;	
                          
function out(e){
	if (old_tr) {
		old_tr.setAttribute("style","background: beige;");
	}
	old_tr = null;
}


function overTabVariation(e) {
	if (e.rowNode == null) return ;
		if (old_tr){
			old_tr.setAttribute("style","background: beige;");
			old_tr = null;
		}
	   	var i = e.rowNode.gridRowIndex;	
		var name=gridVar.getCell(1).getNode(i).textContent;		
		var tr = document.getElementById("tr_"+name);
						
		if (tr == null) return ;
			tr.setAttribute("style","background: #FF1050;");
			old_tr = tr;
			document.getElementById("chrom").scrollTop = tr.offsetTop;//(i-1)*traceHeight;
	   }


function initChrom2(url){
	 	var z = new Array();
			z['name'] = "*";
			
			var store = readStore(url);
			store.fetch({onComplete: buildTable2,onError: gotError});
			
}

	
//fonction pour ecrire du texte sur le graphe
	var makeText = function(surface, text, font, fill, stroke){
                				var t = surface.createText(text);
                				if(font)   t.setFont(font);
                				if(fill)   t.setFill(fill);
                				if(stroke) t.setStroke(stroke);
                				return t;
                			};
                
var gotError = function(error, request){
    alert("The request to the store failed. " +  error);
}
	
	//creation du tableau contenant les graphes
   		function buildTable2(tabData,test){
					
                			var reg = new RegExp("[;]+", "g");									
                			var table =  document.getElementById("tableElectro");					       		              			
                			var l = tabData.length;
							
								var listConnections = new Array();		
								var nb_tconnection =0;
							
							nb_drawing =0;
							max_connection = 0;	
							
							for (var i = 0; i < l; i++) {
								var st = tabData[i].traces+"";
            					var tableau = st.split(reg);
								max_connection= max_connection+tableau.length;
							}
							//max_connection = l;
								if (progress_drawing) {
									progress_drawing.update({
										maximum: max_connection,
										progress: 0
									});
								}
							var split_connection =  max_connection /10 ;
							//pour chaque patient
                			for (var i = 0; i < l; i++) {
								nb_tconnection ++;
            					var tr = document.createElement("tr");															
								tr.setAttribute("id","tr_"+tabData[i].real_name);
            					table.appendChild(tr);
            					
								var st = tabData[i].traces+"";
            					var tableau = st.split(reg);
								st = tabData[i].traces_id+"";
								var trids = st.split(reg);
								 st = tabData[i].traces_strand+"";	
								var strandArray =  st.split(reg);
								var nbplus = tabData[i].nbstrandplus;
            					var l2 = tableau.length;
								var nbtd = l2; 
								
								var data = new Array();	
								var plus = -1;	
								var td0 = document.createElement("th");
								
								
							
								td0.innerHTML = '<center><img src="/icons/Polyicons/user_edit.png">' +tabData[i].name+"</span>";
								td0.colSpan=2;
								tr.setAttribute("style","background:beige;");
            					tr.appendChild(td0);
							
								var z= 0;
            					for (var j = 0; j < l2; j++) {
									if (z%2 == 0 ) {
									
										tr = document.createElement("tr");
            							table.appendChild(tr);
									}
								
									var args = new Array();
									
            						var td1 = document.createElement("td");
            						tr.appendChild(td1);
								
									if ((plus== -1 && strandArray[j] == -1) ||(plus== 1 && strandArray[j] == 1) ){
										drawBlank(td1);
										z++;
										if (z%2 == 0 ) {
											tr = document.createElement("tr");								
            								table.appendChild(tr);
										}
										td1 = document.createElement("td");
										tr.appendChild(td1);
										
									} 
									
									z++;
									plus = strandArray[j];
									var table2 = document.createElement("table");
									td1.appendChild(table2);
									trz = document.createElement("tr");
            						table2.appendChild(trz);
									
									var th_title = document.createElement("th");
									 var turl = url_query+"?type=download_trace&trace_id="+trids[j]+"&projectName=" + projectName;
									th_title.setAttribute("style","background:rgb(240, 240, 240);");
									th_title.innerHTML = "'<a href="+turl+" target='_blank'>"+tableau[j]+"</a>";
									
									trz.appendChild(th_title);
									trz1 = document.createElement("tr");
            						table2.appendChild(trz1);
									
									var td_graph =  document.createElement("table");
									
									var td_temp = trz1.appendChild(td_graph);
									args["div"] = td_temp;
									args["project"] = projectName;
									
									if (!tabData[i].variation_id) {
											args["contig_position"] = tabData[i].contig_position; 
									}
									else {
										//args["variation"] = tabData[i].variation_name;
										args["variation_id"] =tabData[i].variation_id ;
									}
									
								
									//args["contig"] = tabData[i].group;
									//args["trace"] =tableau[j] ;
									args["trace_id"] =trids[j] ;
								
									args["error"] =0;
									listConnections.push(args);
									td_graph.innerHTML = '<img id="wait1" src="../../images/polyicons/wait18trans.gif" align="absmiddle">';
									if (nb_tconnection >= split_connection ){
										readTraceRawDatasArray(listConnections);
										nb_tconnection =0;
										listConnections = new Array();
									}
            						//readTraceRawDatas(tableau[j], projectName, tabData[i].group, tabData[i].real_name,td1);
									
            					}
								
								
            				
                			}
							if (listConnections.length >0){
								
									readTraceRawDatasArray(listConnections);
							}
							
                		}	
						

//fonction pour telecharger les fichiers de sequences 				
function downLoad (e){
	
	var sel = e.grid.selection.getSelected();
	var i=0;
	var z = sel[i];
	var prj1 = z['chemin'];
	
	var el = document.getElementById("filename");
   	el.value = prj1;
   	var fo = document.getElementById("formdownload");
	
	fo.action=url_auth;
	fo.submit();
}		
 
 
 function  traceCover (data,div){
 
 	
	 div.innerHTML=' ';
	  div.style.height= 150+"px";
	  div.style.width= 300+"px";
	  var a =0;
	  var t =0;
	  var c =0;
	  var g= 0;
	  var base = new Array();
	  var real_base = ["A","C","T","G","N","-"];
	 
	  var seq_pos  = data[0].base_position;
	
	  for (var i = 0; i < real_base.length; i++) {
	  	var letter = real_base[i];
		base[letter] = 0;
	  }

	for (var i = 0; i < data[0].traces_sequence.length; i++) {
		var letter = data[0].traces_sequence[i][seq_pos];
		base[letter]++;
	}
	

	  var st1 = new Array();
	
	
	for (var i = 0; i < real_base.length; i++) {
		var letter = real_base[i];
		var colors = "red";
		if (letter ==  data[0].reference_sequence[seq_pos]) {
			colors = "green";
		}
		 st1[i] = {y: base[letter], color: colors};
	}
	 

	
	
	  chart1 = new dojox.charting.Chart2D(div);
	  
     chart1.addAxis("y", {vertical: true});
     chart1.addAxis("x", { 
	 labels: [{value: 1, text: "A"}, {value: 2, text: "C"}, 
		{value: 3, text: "T"}, {value: 4, text: "G"}, {
			value: 5,
			text: "N"
		},{value: 6, text: "-"}
		]
	 });
	
	 
    chart1.addPlot("default", {
      hAxis: "x",
      vAxis: "y",
      type: "ClusteredColumns",
	  ticks: true,
	  precision: 1,
	  labelOffset: 10,
      areas: true,
	  fixed: true,
      gap: 4,
    });
	
     chart1.addSeries("Series 1", st1,{stroke: {color: "black"}, fill: "green"});
	

	  chart1.render();
	 
	traceCoverOld(data,div);  
		
 }
 
            
function traceCoverOld (data,td){
	
		var table = document.createElement("table");
		table.setAttribute("style","font-size:9px")
        td.appendChild(table);
	  var seq_pos  = data[0].base_position;					
	var text = data[0].reference_sequence+"<br>";
	
		var tr = document.createElement("tr");
		table.appendChild(tr);
		for (var j = 0; j < data[0].reference_sequence.length; j++) {
			var t  = document.createElement("td");	
       	 	tr.appendChild(t);
			t.setAttribute("style","background:orange");
				if (j==seq_pos){
				t.setAttribute("class","td_align");
			}
			t.innerHTML = data[0].reference_sequence[j]; 
			
		}
	
	var mid = (data[0].traces_sequence.length/2)+1;
	for (var i=0; i<data[0].traces_sequence.length;i++){
		 var tr  = document.createElement("tr");	
        table.appendChild(tr);
	
		for (var j=0;j<data[0].traces_sequence[i].length;j++){
			var t  = document.createElement("td");	
       	 	tr.appendChild(t);
			var base = data[0].traces_sequence[i][j];
			 if (data[0].traces_sequence[i][j] == "+"){
				t.setAttribute("style","color:rgb(50, 50, 50)");
			}
			else if (data[0].traces_sequence[i][j] == "-"){
				t.setAttribute("style","background:rgb(200, 50, 100)");
				
				
			}
			else if (data[0].traces_sequence[i][j] == "N"){
				t.setAttribute("style","background:rgb(50, 50, 200)");
				
				
			}
			 else if (data[0].traces_sequence[i][j] ==  data[0].reference_sequence[j]) {
				t.setAttribute("style"," background:rgb(173, 202, 86)");
				
			}
			
			else if (data[0].traces_sequence[i][j] == " "){
				t.setAttribute("style","color:rgb(50, 50, 50)");
			}
			else if (data[0].traces_sequence[i][j].length > 1){
					t.setAttribute("style","background:rgb(200, 50, 100)");
			}
			else {
				t.setAttribute("style","background:red");
			}
				if (j==seq_pos){
				t.setAttribute("class","td_align");
			}
			t.innerHTML = base; 
		}
	
	}
	
}
