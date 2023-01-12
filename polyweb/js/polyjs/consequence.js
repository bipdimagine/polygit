/**
 * @author pnitschk
 */

//fonction pour connaitre le numro de la ligne
 

     
	 //layout pour mettre en page le tableau de donnes      		
      


function jsonError(e,r) {
	alert("error json");
	console.error(e);
}

function jsonOk(e,r) {
	alert("ok json");
	console.error(e);
}
dojo.require("dojo.data.api.Read");
var methods = new Array("calling","unifiedgenotyper","dibayes","mpileup","gatk","casava");
var methods_abrg = new Array("C","uni","dib","mp","gatk","cas");
 function viewAnnex(items, request){

 	document.getElementById("csmall").style.visibility = "visible";	
	 xpane=0;
	 ypane=0; 
	 nbpane=0;	
			var array = new Array();
	
			var ri= 0;

			for (var i=0;i<items[0].patient_name.length;i++){
			//	if (items[0].patient_name[i] != "6ZD01") continue;
			 var tname= items[0].patient_name[i];
			
				var trace_id =i;
				var trace_name= "None";
				if (items[0].traces_id) {
					 trace_id =items[0].traces_id[i];
					 trace_name= items[0].traces_name[i];
					 
				}
				var st = items[0]["text"];
				st = st+"";
					st = st.replace("/","_");
				
		
						var poly_id = items[0]["chromosome"]+"_"+items[0]["start"]+"_"+st;
						
					array[ri] = {
						fam:items[0].pedigree_fam[i],
						name: items[0].patient_name[i],
					//	score: items[0].patient_score[i],
						base: items[0].patient_base[i],
						variation_id: items[0].id,
						traces_id: trace_id,
						traces_name: trace_name,
						patient_id : items[0].patient_id[i],
						valid: items[0]["valid!"+tname], 
						type:items[0].pedigree_type[i],
						
						poly_id:poly_id,
						
						
					};
					for (var z = 0; z < methods.length; z++) {
					   
						var t = "patient_" + methods[z];
					    // console.dir(items[0]);
			
						if (items[0][t]) {
							var t2 = "patient_" + methods[z] + "2";
							var t4 = "patient_" + methods[z] + "4";
							var t3 = "patient_" + methods[z] + "3";
							var theho =  "patient_" + methods[z] + "heho";
						
							if (items[0][t][i] >= 0 &&  items[0][t3][i]) {
							
							 
								array[ri]["ref"] = items[0][t2][i];
								array[ri]["alt"] = items[0][t + "3"][i];
								array[ri]["status"] = items[0][theho][i];
								array[ri]["dp"] = items[0][t + "4"][i];
								array[ri]["met"] = methods_abrg[z];
								
								break;
								
							}
						}
					}
				
				
			/*	else {
					
					array[ri] = {
						name: items[0].patient_name[i],
						score: items[0].patient_score[i],
						base: items[0].patient_base[i],
						variation_id: items[0].id,
						patient_id : items[0].patient_id[i],
						valid: items[0]["valid!"+tname], 
						
					}
				}*/
				ri ++;
			}
		
			var jsonStore = new dojo.data.ItemFileWriteStore({
				data:  {
				"identifier": "name",
				"items": array
					}
			});
		
		var table =  dojo.byId(current_electro);//document.createElement("table");
		table.innerHTML = ' ';
		
		grid_current_annex.setStore(jsonStore);
		jsonStore.close();
		
		if (patient_name_filter) {
			var filter = new Array();
			filter["name"] = "*" + patient_name_filter + "*";
			grid_current_annex.setQuery(filter);
		}
		
			
			if (items[0]["consequence!all"] == "Intronic" || items[0]["consequence!all"] == "Intergenic") {
				document.getElementById("pconseq").style.visibility = "hidden";	
				return;
			}
		
		
			/*tableau des consequences*/
			document.getElementById("pconseq").style.visibility = "visible";;	
			
			
			document.getElementById("variation_select").innerHTML=items[0].name;
			document.getElementById("consequence_select").innerHTML=items[0]["consequence!all"];
			/*var fil = filterTranscript.getValue();
			if (fil.length == 0) fil="*";
			if (fil == "all") fil="*";
		    */
		
			var array2 = new Array();
		//	array2 = items[0].tab_consequences;
		
		
			for (var i = 0; i < items[0].tab_consequences.length; i++) {
				 if  (items[0].tab_consequences[i] == null) continue;
				array2[i] = {
						AA_position: items[0].tab_consequences[i].AA_position+"",
						name:  items[0].tab_consequences[i].name+i,
						gene: items[0].tab_consequences[i].gene,
						AA_protein:items[0].tab_consequences[i].AA_protein,
						AA_variation:items[0].tab_consequences[i].AA_variation,
						cdna_position:items[0].tab_consequences[i].cdna_position,
						consequence:items[0].tab_consequences[i].consequence,
						description:items[0].tab_consequences[i].description,
						exon:items[0].tab_consequences[i].exon,
						polyphen:items[0].tab_consequences[i].polyphen,
						protein:items[0].tab_consequences[i].protein+"",
						transcript:items[0].tab_consequences[i].transcript+"",
						polyphen_html:items[0].tab_consequences[i].polyphen_html,
						polyphen_status:items[0].tab_consequences[i].polyphen_status,
						sift_status:items[0].tab_consequences[i].sift_status,
						protein_position:items[0].tab_consequences[i].protein_position,
						cds_position:items[0].tab_consequences[i].cds_position,
						nomenclature:items[0].tab_consequences[i].nomenclature,
					};
			
			}
		
			var jsonStore2 = new dojo.data.ItemFileWriteStore({
				data:  {
				"identifier": "name",
				"items": array2
					}
			});
			
			gridConsequence.setStore(jsonStore2);
			resizeDiv('gridConsequence','pconseq');
			
 }
 	

//fonction appele quand on passe la souris sur une ligne du tableau (une variation) qui permet d'afficher les consequence de la variation
	consequenceOver = function(e){
		
	grid_current_annex = gridAnnex;
	current_electro= "tableElectro";
		var i = e.rowNode.gridRowIndex;
		//recupre les donnes de la ligne qui sont dans la 1ere cellule (commence 0)
		var name = gridVar.getCell(1).getNode(i).textContent;
		viewAnnex(gridVar.selection.getSelected());
	
		
	}
	

viewElectro = function(item,table){
		var tr_header = document.createElement("tr");	
		table.appendChild(tr_header);
		var tr_electro = document.createElement("tr");	
		table.appendChild(tr_electro);
		var tids = new Array();
		for (var i = 0; i < item.traces_id.length; i++) {
			 tids[item.traces_name[i]] = item.traces_id[i];
		}
		
		var project;
		if (item.project) { project = item.project; }
		else { project = projectName; }
		
		for (var i = 0; i < item.traces_id.length; i++) {
			var url_arg = "?project=" + project  +"&poly_id="+item.poly_id+"&variation_id=" + item.variation_id + "&trace_id=" + item.traces_id[i] + "&patient_name=" + item.name;
			//var url_arg = "?project=" + projectName  +"&poly_id="+item.poly_id+"&variation_id=" + item.variation_id + "&patient_name=" + item.name;
			dojo.xhrGet({ // 
				// The following URL must match that used to test the server.
				url: url_chart + url_arg,
				handleAs: "json",
				timeout: 100000, // Time in milliseconds
				// The LOAD function will be called on a successful response.
				load: function(responseObject, ioArgs){
					// Now you can just use the object
					//pour les projets de squencage classique on trace des chromatogrammes
					var th = document.createElement("th");
					tr_header.appendChild(th);
					var name= responseObject[0].name;
					if (name == "toto" ) name = "Array";
					//var turl = url_query + "?type=download_trace&trace_id=" + tids[responseObject[0].name] + "&projectName=" + projectName;
					//th.innerHTML = "<a href=" + turl + " target='_blank'>" + name + "</a>";
					if (/[A-Za-a]/.test(project)) { th.innerHTML = "<b>Project: <FONT color=red>" + project + "</FONT> | Patient: <FONT color=red>" + item.name + "</FONT></b>"; }
					else { th.innerHTML = "<b>Project/Patient ID: <FONT color=red>" + item.id + "</FONT></b>"; }
					var td1 = document.createElement("td");					
					tr_electro.appendChild(td1);
					//var tid = responseObject[0].trace_id;
				
					if (responseObject[0].type == 1) {
					//	var td1 = document.getElementById("el_"+tid);
										
						TraceCharts(responseObject, td1, item.traces_id[i]);
						
					}
					else if (responseObject[0].type == 4) {
					
						TraceRawCnvs(td1,responseObject[0].traces_sequence,responseObject[0].other_values);
					}
					else 
					
						if (responseObject[0].type == 2) {
							TraceHisto(responseObject, td1);
							
						}
						//pour les projet de puce de resquenage on trace des histogrammes
						else {
							traceCover(responseObject, td1);
						}
					
				},
				// The ERROR function will be called in an error case.
				error: function(response, ioArgs){ // 
					console.error(response);
					//en cas d'erreur
					td1.innerHTML = 'Loading Error ...';
					return response; // 
				}
			});
			
		}
}	

var pane = new Array();

	var xpane=0;
	var ypane=0; 
	var nbpane =0;
selectPatient1 = function(e){
	

 
 
	
	
		for (var i = 0; i < grid_current_annex.selection.getSelected().length; i++) {
			var name = grid_current_annex.selection.getSelected()[i].name;
			if (createFloatPanel(grid_current_annex.selection.getSelected()[i],xpane,ypane)) {
				if (nbpane%4 == 0){
					xpane = 5;
					ypane =5;
				}
				
				var table = dojo.byId(name + "_table_electro");
				viewElectro(grid_current_annex.selection.getSelected()[i], table);
				if (nbpane%2 == 0) {
					//ypane = 0;
					xpane += (210*grid_current_annex.selection.getSelected()[i].traces_id.length);
				}
				else {
					xpane =0;
					ypane += (205);
				}
				nbpane ++;
			}
		}
	
		
		
   
		
	}
	
function selectPatient(e) { selectPatientGlobal(e); }

function selectPatient_yours(e) { selectPatientGlobal(e, "gridCovYourProj"); }

function selectPatient_others(e) { selectPatientGlobal(e, "gridCovOtherProj"); }

selectPatientGlobal = function(e, gridName){
	var grid_used;
	if (gridName) { grid_used = dijit.byId(gridName); }
	else { grid_used = grid_current_annex; }
	var item = grid_used.getItem(e.rowIndex);
	var table =  dojo.byId(current_electro);//document.createElement("table");
	table.innerHTML = ' ';
	
	var tr = document.createElement("tr");	
																
    table.appendChild(tr);

	
	//alert("BEFORE FOR");
	//var len = grid_used.selection.getSelected().length;
	//alert("LENGTH: "+len);
		//for (var i = 0; i < grid_used.selection.getSelected().length; i++) {
			//alert("IN FOR");
			//var name = grid_used.selection.getSelected()[i].name;
			
			var name = grid_used.store.getValue(item, "name"); 
			
			
				//viewElectro(grid_used.selection.getSelected()[i], table);
				viewElectro(item, table);
		//}
	
		
	}	
	
	
function createFloatPanel(item,x,y) {
	var table;
 	var name= item.name;
 if (!dojo.byId(name+"_floater")) {
 	
   var dive = dojo.byId("electro");
  
 	var div1 = document.createElement("div");
	 dive.appendChild(div1);

 

 table = document.createElement("table");
 table.setAttribute("id",name+"_table_electro");
 div1.appendChild(table);


pane[name] = new dojox.layout.FloatingPane({
 		'title': name,
 		'id': name+"_floater",
 		'closeable': true,
 		'resizable': true,
		'background':"#000000",
 		'dockable': true
 	}, div1);
	var l = item.traces_id.length;
	pane[name].domNode.style.left = x + "px";
	pane[name].domNode.style.zIndex ="100";
  	pane[name].domNode.style.top = y + "px";
  pane[name].resize({ 'w': (210*l), 'h': 200 });
  

 
 dojo.connect(pane[name],"close",function(e){
    pane[name] = undefined;
});

   pane[name].startup(); 
   
 }
 else {
 	pane[name].show();
	return false;
 }
	return true;
}			
                         