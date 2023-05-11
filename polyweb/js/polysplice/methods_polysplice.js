var is_first_launch_polysplice;
var used_dejavu_value = 15;
var used_dejavu_percent_value = 96;
var used_score_value = 10;
var url_last_launch_patient;
var value_dejavu_percent_similar = 96;
var value_score = 10;
var value_dejavu = 15;
var var_load_polysplice = 0;
function load_tab_polysplice(patient_name) {
	alert('POLYSPLICE');
	if (var_load_polysplice == 1) {
		refresh_tab_polysplice();
		return;
	}
	for (var i=0;i<tsamples.length;i++) {
    	container_editor_polysplice.addChild(tab_editor_polysplice[tsamples[i].label]);
		tab_editor_polysplice[tsamples[i].label].add = 1;
	}
	var_load_polysplice = 1;
}

function refresh_tab_polysplice(){
	if ( !(dijit.byId("gridPatients"))) {
		return 1;
	}
	var tsamples2 = dijit.byId("gridPatients").selection.getSelected();
	var selected = new Array();
	 for (var i=0;i<tsamples2.length;i++) {
		selected[tsamples2[i].label] = 1;
	}
	var first = 0;
	for (var i=0;i<tsamples.length;i++) {
		var name = tsamples[i].label
		if (selected[name] == 1){
			if (tab_editor_polysplice[tsamples[i].label].add == 0){
				container_editor_polysplice.addChild(tab_editor_polysplice[tsamples[i].label]);
				tab_editor_polysplice[tsamples[i].label].add = 1;
			}
			if (first < 1){
			tab_editor_polysplice[name].onShow();
			first ++;
		}
		}
		else {
			if (tab_editor_polysplice[tsamples[i].label].add == 1) {
				container_editor_polysplice.removeChild(tab_editor_polysplice[tsamples[i].label]);
				tab_editor_polysplice[tsamples[i].label].add = 0;
				tab_editor_polysplice[name].set('selected',false);
			}
		}
	}
}

function update_score_span_value() {
	var value_score_slider = document.getElementById("slider_score").value;
	if (value_score_slider == 0) {
		document.getElementById("nb_score").innerHTML = '0%';
		value_score = '0';
	}
	else {
		value_score = value_score_slider * 10;
		document.getElementById("nb_score").innerHTML = value_score+'%';
	}
}

function update_dejavu_span_value() {
	value_dejavu = document.getElementById("slider_dejavu").value;
	if (value_dejavu == 51) {
		document.getElementById("nb_max_dejavu_patients").innerHTML = 'ALL';
	}
	else {
		document.getElementById("nb_max_dejavu_patients").innerHTML = value_dejavu;
	}
}

function update_dejavu_percent_span_value() {
	value_dejavu_percent_similar = 90 + parseInt(document.getElementById("slider_dejavu_percent").value);
	document.getElementById("nb_percent_dejavu").innerHTML = value_dejavu_percent_similar;
}

function load_polysplice_from_polyviewer() {
	tab_editor_polysplice[selected_patient_name].toto = 1;
	tab_editor_polysplice[selected_patient_name].onShow();
}
			
function launch_dejavu_span_value(use_patient_name) {
	var only_gene_positions = document.getElementById('input_gene_positions').value;
	if (value_dejavu == 101) { show_waiting("Filter DejaVu - ALL patients"); }
	else { show_waiting("Filter DejaVu max " + value_dejavu + " patients (" + value_dejavu_percent_similar + "% coord similar)"); }
	var url = url_last_launch_patient; 
	if (value_dejavu < 101) { url += "&dejavu=" + value_dejavu; }
	url += "&dejavu_percent=" + value_dejavu_percent_similar;
	url += "&min_score=" + value_score;
	if (is_first_launch_polysplice == 1) {
		url += "&only_dejavu_ratio_10=1";
	}
	else if ($('#b_dejavu_min_ratio_10').prop('checked')) {
		url += "&only_dejavu_ratio_10=1";
	}
	if (only_gene_positions == '') {}
	else if (only_gene_positions.match(':')) {
		url += "&only_positions=" + only_gene_positions;
	}
	else {
		url += "&only_gene=" + only_gene_positions;
	}
	used_dejavu_percent_value = value_dejavu_percent_similar;
	used_score_value = value_score;
	
	var tab_id = 'c1polysplice'+use_patient_name;
	load_polysplice_url(use_patient_name, tab_id, url);
	
}

function view_linked_junctions(patient_name, transcript_name, junction_vector_ids, this_junction, min_ratio) {
	show_waiting("Show linked junctions found");
	var url = url_path + "/rnaseq/rna_junctions_linked_view.pl";
	url += "?project=" + project_name;
	url += "&patient=" + patient_name;
	url += "&transcript=" + transcript_name;
	url += "&min_ratio=" + min_ratio;
	url += "&this_junction=" + this_junction;
	url += "&vector_ids=" + junction_vector_ids;
	$.getJSON( url, function( data ) {
		document.getElementById("div_linked_junction").innerHTML = "";
		document.getElementById("div_linked_junction").innerHTML += data.html;
		for (var i = 0; i < data.list_graphs.length; i++){
			y_norm = data.list_graphs[i].y_norm;
			x_common = data.list_graphs[i].x_common;
			y_common = data.list_graphs[i].y_common;
			plot_patient_name = data.list_graphs[i].patient_name;
			var j = i + 1;
			var plotid = 'plot'+j;
			do_graph(plotid);
		}
		dijit.byId('waiting').hide();
		dijit.byId('dialog_linked_junction').show();
    })
}

var y_norm;
var x_common;
var y_common;
var plot_patient_name;
function do_graph(plotname) {
	var trace1 = {
			x:x_common,
			y:y_common,
			yaxis: 'y1',
			mode: 'lines+markers',
			showlegend:false,
			name:plot_patient_name,
			opacity: 0.6,
			marker: {
				color: 'red',
				size:2
			}
	};
	
	var trace2 = {
			x:x_common,
			y:y_norm,
			yaxis: 'y1',
			mode: 'lines+markers',
			showlegend:false,
			name:'Others',
			opacity: 0.3,
			marker: {
				color: 'blue',
				size:2
			}
	};

	var data=[trace1, trace2];

	var layout = {
		yaxis:   {title : 'Nb Found Reads'},
		xaxis : {title : 'Genomic Position'},
		height: 600,
	};

	var config = { modeBarButtonsToRemove: ['pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d','toggleSpikelines', 'hoverClosestCartesian', 'hoverCompareCartesian'], displayModeBar: false};
	Plotly.newPlot(plotname, data, layout,config);
}


var hash_polysplice_loaded = {};
function load_polysplice_url(patient_name, tab_id, url_polysplice) {
	//dijit.byId('waiting').show();
	//document.getElementById(tab_id).innerHTML = "<img src=\"https://www.polyweb.fr/polyweb/images/polyicons/wait18trans.gif\" align=\"absmiddle\"> Loading Data...";
	$.getJSON( url_polysplice, function( data ) {
    	$.each( data, function( key, val ) {
    		if (key == 'html') {
    			document.getElementById(tab_id).innerHTML = "";
				var myNode = document.getElementById(tab_id);
		    	var myPanel = $(val);
				myPanel.appendTo(myNode);
				
				var list_tables = String(data.tables_id).split(',');
				for (var i = 0; i < list_tables.length; i++){
					enable_table_search(list_tables[i]);
				}
				
				//dijit.byId('waiting').hide();
			}
		});
    })
	hash_polysplice_loaded[patient_name] = 1;
}

function enable_table_search(id) {
	//$('#'+id).bootstrapTable();
}
			
function zoom_file(svg_legend, svg_file) {
	var span = document.createElement("span");
	span.innerHTML = "<br>"+svg_legend;
	var elem = document.createElement("img");
	elem.setAttribute("src", svg_file);
	elem.setAttribute("height", "200%");
	elem.setAttribute("width", "auto");
	elem.setAttribute("style", "transform: translate(0%, -30%);");
	elem.setAttribute("alt", "N.A.");
	document.getElementById("dialog_svg").innerHTML = "";
	document.getElementById("dialog_svg").appendChild(span);
	document.getElementById("dialog_svg").appendChild(elem);
	dijit.byId('dialog_svg').show();
}

function view_pdf_list_files(files) {
	update_sashimi_pdf_file(files, 0);
	dijit.byId('dialog_sashimi').show();
}

function show_waiting(text) {
	dijit.byId('waiting').show();
}
			
function update_sashimi_pdf_file (files, nb) {
	document.getElementById("dialog_sashimi").innerHTML = "<br>Loading...<br>";
	var list_files = files.split(';');
	var max = list_files.length;
	var file = list_files[nb];
	var disabled_previous = 'disabled';
	var nb_previous = 0;
	if (parseInt(nb) > 0) {
		nb_previous = parseInt(nb) - 1;
		disabled_previous = '';
	}
	var disabled_next = 'disabled';
	var nb_next = 4;
	if (parseInt(nb) < 4) {
		nb_next = parseInt(nb) + 1;
		disabled_next = '';
	}
	var elem = "<br><br><div><center><object data='"+file+"' width='95%' height='95%'></object><div style='background-color:#607d8b;'><br><br><label for='zoom_pdf_sashimi' class='form-label' style='color:white;'><nobr><button " + disabled_previous + " style='background-color:#607d8b;border: solid 0px white;border-radius:10px;' onclick=\"update_sashimi_pdf_file('" + files + "','" + nb_previous + "')\"><span class='glyphicon glyphicon-plus-sign' style='color:white;font-size:12px;'></span></button>&nbsp;&nbsp; ZOOM Sashimi Plot &nbsp;&nbsp;<button style='background-color:#607d8b;border: solid 0px white;border-radius:10px;' " + disabled_next + " onclick=\"update_sashimi_pdf_file('" + files + "','" + nb_next + "')\"> &nbsp;&nbsp;<span class='glyphicon glyphicon-minus-sign' style='color:white;font-size:12px;'></span></button></label><input disabled type='range' class='form-range' min='0' max='" + (parseInt(max)-1) + "' value='" + nb + "' style='width:300px;' id='zoom_pdf_sashimi'></nobr><br><br><button style='font-size:16px;margin-bottom:15px;' onClick='dijit.byId(\"dialog_sashimi\").hide();'>Close PLOT</button></center><br><br></div>";
	document.getElementById("dialog_sashimi").innerHTML = elem;
}

function view_pdf_list_files(files) {
	update_sashimi_pdf_file(files, 0);
	dijit.byId('dialog_sashimi').show();
}

function view_deja_vu_rna_junction(project_name, patient_name, junction_vector_id) {
	document.getElementById("div_dejavu_junction").innerHTML = "";
	show_waiting("Search DejaVu");
	var url = url_path + "/polydejavu/dejavu_jonction.pl";
	url += "?project=" + project_name;
	url += "&patient=" + patient_name;
	url += "&vector_id=" + junction_vector_id;
	url += "&dejavu_percent=" + value_dejavu_percent_similar;
	$.getJSON( url, function( data ) {
		document.getElementById("div_dejavu_junction").innerHTML = "";
		document.getElementById("div_dejavu_junction").innerHTML += "<center><span style='margin:5px;'><b>In Patient <span style='color:green;'>" + patient_name + "</span></b></span></center>" + data.html_my_variant;
		document.getElementById("div_dejavu_junction").innerHTML += "<br><center><span style='margin:5px;'><b>In DejaVu Junctions</b></span></center>" + data.html;
		enable_table_search(data.table_id_my_variant);
		enable_table_search(data.table_id);
		dijit.byId('waiting').hide();
		dijit.byId('dialog_dejavu_junction').show();
    })
}

function view_dejavu_nb_int_this_run_patients(project_name, patient_name, junction_vector_id, min_ratio) {
	document.getElementById("div_dejavu_junction").innerHTML = "";
	show_waiting("Search DejaVu (in this run)");
	var url = url_path + "/polydejavu/dejavu_jonction.pl";
	url += "?is_dejavu_inthis_run=1&project=" + project_name;
	url += "&patient=" + patient_name;
	url += "&vector_id=" + junction_vector_id;
	url += "&min_ratio=" + used_score_value;
	$.getJSON( url, function( data ) {
		document.getElementById("div_dejavu_junction").innerHTML = "";
		document.getElementById("div_dejavu_junction").innerHTML += "<center><span style='margin:5px;'><b>In Patient <span style='color:green;'>" + patient_name + "</span></b></span></center>" + data.html_my_variant;
		document.getElementById("div_dejavu_junction").innerHTML += "<br><center><span style='margin:5px;'><b>In DejaVu Junctions</b></span></center>" + data.html;
		enable_table_search(data.table_id);
		enable_table_search(data.table_id_my_variant);
		dijit.byId('waiting').hide();
		dijit.byId('dialog_dejavu_junction').show();
    })
}

function hide_view_table_gene_splices(table_id) {
 	document.getElementById("content_res_splices").innerHTML = "";
	const node = document.getElementById(table_id).lastChild;
	const clone = node.cloneNode(true);
	document.getElementById("content_res_splices").appendChild(clone);
	dijit.byId('dialog_splices').show();
}

function select_patient(elem) {
	console.log(elem);
	var id1 = elem.id;
	var id2 = id1.replace('check_', 'b_selected_patient_');
	document.getElementById(id2).checked = elem.checked;
}

function select_all_patients(elem) {
	if (tsamples.length == 0) { tsamples = storeP; }
	for (var i=0;i<tsamples.length;i++) {
		var label = tsamples[i].label;
		document.getElementById('check_' + label).checked = elem.checked;
		document.getElementById('b_selected_patient_' + label).checked = elem.checked;
	}
}
