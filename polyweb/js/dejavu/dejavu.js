var dataStore;
var dataStore_HeComp;
var dataStore_log;
var dataStore_details;
var urlDejaVu = url_path + "polydejavu/dejavu_cgi.pl";
var urlHeComp = url_path + "polydejavu/hetero_comp_cgi.pl";
var url_checkLogin = url_path + "polydejavu/checkLog_cgi.pl";
var url_varInfos = url_path + "polydejavu/dejavu_var_infos.pl";
var nbResults;
var LoginOk = false;
var infobulle_visible = false;
var hashPatient = new Object();
var anc_onglet = 'projects';
var panelFilters_defaultTitle = "Graphics and grid results filters";
var panelHeComp_defaultTitle = "Heterozygotes composites";
var panelDetails_defaultTitle = "Selected variant consequence(s) by each transcript";
var panelResults_defaultTitle = "Results";
var originLaunchInput_dejavu;
var valId_forHe = false;
var nbAuthProj = "?";
var varId_selected;
var current_electro = "graphAlign";

 

	
function check_slider_cat_annotations() {
	var value = parseInt(dijit.byId("slider_impact").get('value'));
	dijit.byId("filter_intergenic").set('checked', false);
	dijit.byId("filter_intronic").set('checked', false);
	dijit.byId("filter_upstream_downstream").set('checked', false);
	dijit.byId("filter_silent").set('checked', false);
	dijit.byId("filter_utr").set('checked', false);
	dijit.byId("filter_pseudogene").set('checked', false);
	dijit.byId("filter_nonframeshift").set('checked', false);
	dijit.byId("filter_nonsynonymous").set('checked', false);
	dijit.byId("filter_nonsynonymous").set('checked', false);
	dijit.byId("filter_splice_site").set('checked', false);
	dijit.byId("filter_ncrna").set('checked', false);
	dijit.byId("filter_phase").set('checked', false);
	dijit.byId("filter_stop").set('checked', false);
	dijit.byId("filter_frameshift").set('checked', false);
	dijit.byId("filter_essential_splicing").set('checked', false);
	dijit.byId("filter_maturemirna").set('checked', false);
	
	if (value == 1) {
		dijit.byId("filter_phase").set('checked', true);
		dijit.byId("filter_stop").set('checked', true);
		dijit.byId("filter_frameshift").set('checked', true);
		dijit.byId("filter_essential_splicing").set('checked', true);
		dijit.byId("filter_maturemirna").set('checked', true);
	
	}
	else if (value == 2) {
		dijit.byId("filter_nonframeshift").set('checked', true);
		dijit.byId("filter_nonsynonymous").set('checked', true);
		dijit.byId("filter_splice_site").set('checked', true);
		dijit.byId("filter_phase").set('checked', true);
		dijit.byId("filter_stop").set('checked', true);
		dijit.byId("filter_frameshift").set('checked', true);
		dijit.byId("filter_essential_splicing").set('checked', true);
		dijit.byId("filter_maturemirna").set('checked', true);
	}
	else if (value == 3) {
		dijit.byId("filter_ncrna").set('checked', true);
		dijit.byId("filter_intronic").set('checked', true);
		dijit.byId("filter_upstream_downstream").set('checked', true);
		dijit.byId("filter_silent").set('checked', true);
		dijit.byId("filter_utr").set('checked', true);
		dijit.byId("filter_pseudogene").set('checked', true);
		dijit.byId("filter_nonframeshift").set('checked', true);
		dijit.byId("filter_nonsynonymous").set('checked', true);
		dijit.byId("filter_splice_site").set('checked', true);
		dijit.byId("filter_ncrna").set('checked', true);
		dijit.byId("filter_phase").set('checked', true);
		dijit.byId("filter_stop").set('checked', true);
		dijit.byId("filter_frameshift").set('checked', true);
		dijit.byId("filter_essential_splicing").set('checked', true);
		dijit.byId("filter_maturemirna").set('checked', true);
	}
	else if (value == 4) {
		dijit.byId("filter_ncrna").set('checked', true);
		dijit.byId("filter_intronic").set('checked', true);
		dijit.byId("filter_intergenic").set('checked', true);
		dijit.byId("filter_upstream_downstream").set('checked', true);
		dijit.byId("filter_silent").set('checked', true);
		dijit.byId("filter_utr").set('checked', true);
		dijit.byId("filter_pseudogene").set('checked', true);
		dijit.byId("filter_nonframeshift").set('checked', true);
		dijit.byId("filter_nonsynonymous").set('checked', true);
		dijit.byId("filter_splice_site").set('checked', true);
		dijit.byId("filter_ncrna").set('checked', true);
		dijit.byId("filter_phase").set('checked', true);
		dijit.byId("filter_stop").set('checked', true);
		dijit.byId("filter_frameshift").set('checked', true);
		dijit.byId("filter_essential_splicing").set('checked', true);
		dijit.byId("filter_maturemirna").set('checked', true);
	}
}


function defineArrayButton() {
    var array_button = {};
//    array_button["notcosmic"] = "";
//    array_button["cosmic"] = "";
//    array_button["cnv"] = "";
//    array_button["deletion"] = "";
//    array_button["insertion"] = "";
//    array_button["substitution"] = "";

	array_button['maturemirna'] = 'filter_maturemirna';
	array_button['essential_splicing'] = 'filter_essential_splicing';
	array_button['frameshift'] = 'filter_frameshift';
	array_button['stop'] = 'filter_stop';
	array_button['phase'] = 'filter_phase';
	array_button['ncrna'] = 'filter_ncrna';
	array_button['splice_site'] = 'filter_splice_site';
	array_button['nonsynonymous'] = 'filter_nonsynonymous';
	array_button['nonframeshift'] = 'filter_nonframeshift';
	array_button['non-frameshift'] = 'filter_nonframeshift';
	array_button['pseudogene'] = 'filter_pseudogene';
	array_button['utr'] = 'filter_utr';
	array_button['silent'] = 'filter_silent';
	array_button['upstream_downstream'] = 'filter_upstream_downstream';
	array_button['downstream'] = 'filter_upstream_downstream';
	array_button['upstream'] = 'filter_upstream_downstream';
	array_button['intronic'] = 'filter_intronic';
	array_button['intergenic'] = 'filter_intergenic';
    return (array_button);
}

function get_polyquery_var_annotation_checked() {
	var ab = defineArrayButton()
	var filters = '&filters_annot=';
	var has_filters = 0;
	for (annot_name in ab) {
        if (dijit.byId(ab[annot_name]).checked == false){
            filters += annot_name + "+";
            has_filters = 1;
        }
    }
    if (has_filters == 1) { return filters; }
    return '';
}
	
function check_slider_frequence_query() {
	return;
	if (dijit.byId("slider_frequence_query").get('value') > max_slider_frequency_query) {
		dijit.byId("slider_frequence_query").set('value', max_slider_frequency_query);
	}
}
	
function click_annot(value) {

}

function change_slider_dejavu_text(value) {
	text_dejavu.innerHTML = '<i><b>Max</b> Projects: </i>';
    if (value == 0)  {
        text_dejavu.innerHTML = '<i><b><font color="#EA3434">Impossible !</font></i></b>';
        nb_dejavu.innerHTML = '';
    }
    else if (value == 100) {
        text_dejavu.innerHTML = '';
        nb_dejavu.innerHTML = '';
    }
    else if (value <= 10) { nb_dejavu.innerHTML = '<font color="#EA3434">' + value + '</font>'; }
    else if (value <= 30) { nb_dejavu.innerHTML = '<font color="#EAA134">' + value + '</font>'; }
    else if (value <= 50) { nb_dejavu.innerHTML = '<font color="#EACC34">' + value + '</font>'; }
    else                  { nb_dejavu.innerHTML = '<font color="#28A24C">' + value + '</font>'; }
    dijit.byId("slider_dejavu").set('value', value);
}

function get_value_slider_dejavu() {
	var value = parseInt(dijit.byId("slider_dejavu").get('value'));
	return value;
}


function get_value_slider_freq() {
	if      (dijit.byId("slider_frequence_query").get("value") == 1)    { return "&filters_freq=freq_001+freq_01+freq_05+freq_1+"; }
	else if (dijit.byId("slider_frequence_query").get("value") == 2)    { return "&filters_freq=freq_01+freq_05+freq_1+"; }
	else if (dijit.byId("slider_frequence_query").get("value") == 3)    { return "&filters_freq=freq_05+freq_1+"; }
	else if (dijit.byId("slider_frequence_query").get("value") == 4)    { return "&filters_freq=freq_1+"; }
	else {  return ''; }
	return filters;
}

function check_slider_dejavu() {
	var value = get_value_slider_dejavu();
	change_slider_dejavu_text(value);
	return;
	if (value > max_slider_dejavu) {
		change_slider_dejavu_text(max_slider_dejavu);
	}
	else { change_slider_dejavu_text(value); }
}

function checkKeyPressed(event) {
    if (event.keyCode == 13) {
    	launchId(dijit.byId('id').getValue());
    }
}

var to_do_stats = 1; 
function init_dejavu() {
	disabledXlsExport();
	document.getElementById("id").disabled = true;
	if (param("input")) { document.getElementById("id").value = param("input"); }
    else if (document.getElementById("id").value == "") { document.getElementById("id").value = "Write your ID here"; }
    if (param("build")) { document.getElementById("aff_build").innerHTML = param("build"); }
    else { document.getElementById("aff_build").innerHTML = "HG38"; }
    var url_input = window.location.search;
    var username = null;
    var passwd = null;
    if (dojo.cookie("username")) { username = dojo.cookie("username"); }
    if (dojo.cookie("passwd")) { passwd = dojo.cookie("passwd"); }
    if ((!username) || (!passwd)) { dijit.byId('login').show(); }
    var currentTime = new Date();
    var month = currentTime.getMonth() + 1;
    if (month < 10) { month = '0'+month; }
    var day = currentTime.getDate();
    var year = currentTime.getFullYear();
    document.getElementById("stat_description").textContent = day + " - " + month + " - " + year;
    project_stats();
    checkPassword_dejavu();
}

function onClickTextSearch() {
	if (document.getElementById("id").value == "Write your ID here") {
		document.getElementById("id").value = "";
	}
}

function checkPassword_dejavu(tryToLog) {
	var user = '';
	var pass = '';
	if (dojo.cookie("username")) { user = dojo.cookie("username"); }
	if (dojo.cookie("passwd")) { pass = dojo.cookie("passwd"); }
	if (tryToLog == true) {
		user = document.getElementById('username').value;
		pass = document.getElementById('passwd').value;
	}
	dataStore_log = new dojo.data.ItemFileReadStore({ url: url_checkLogin + "?login=" + user + "&pwd=" + pass });
	var gridLog = dijit.byId("gridLog");
	gridLog.setStore(dataStore_log);
	gridLog.store.fetch({
		onComplete: function(items){
			if (checkLogged(user, pass, tryToLog)) {
				nbAuthProj = countUserProjects();
				document.getElementById("aff_nbauthproj").innerHTML = nbAuthProj;
				gridLog.filter({projName: "None"});
				gridLog.startup();
				document.getElementById("aff_log").innerHTML = user;
				dijit.byId('id').set('disabled', false);
				dijit.byId('b_run').set('disabled', false);
				if (document.getElementById("id").value != "Write your ID here") { launchId(dijit.byId("id").getValue()); }
			}
		}
	});
	gridLog.startup();
	return;
}

function checkLogged(user, pass, tryToLog){
	var nb = countUserProjects();
	if (nb > 0) {
		dijit.byId('login').hide();
		dojo.cookie("username", user, { expires: 1 });
		dojo.cookie("passwd", pass, { expires: 1 });
		LoginOk = true;
		return true;
	}
	else {
		if (tryToLog == true) { alert("User '" + user + "' not found or bad password or no project associated to this user..."); }
    	dijit.byId('logind').show();
	}
	return;
}

function countUserProjects() {
	var grid = dijit.byId("gridLog");
	var nb = grid.rowCount;
	if (nb > 0) { LoginOk = true; }
	return nb;
}

function affiche_infobulle(text) {
	if(infobulle_visible == false) {
		document.getElementById("curseur").style.visibility="visible";
		document.getElementById("curseur").innerHTML = text;
		infobulle_visible = true;
	}
}

function cache_infobulle() {
	if(infobulle_visible == true) {
		document.getElementById("curseur").style.visibility="hidden";
		infobulle_visible = false;
	}
}

function popup_help() {
	width = 300;
	height = 200;
	if(window.innerWidth) {
		var left = (window.innerWidth-width)/2;
		var top = (window.innerHeight-height)/2;
	}
	else {
		var left = (document.body.clientWidth-width)/2;
		var top = (document.body.clientHeight-height)/2;
	}
	window.open('help_dejavu.html','HELP DejaVu','menubar=no, scrollbars=no, top='+top+', left='+left+', width='+width+', height='+height+'');
}

function onRowClick_gridResults(value) {
	onRowClick_details(value);
	onRowClick_varInfos(value);
}

var alamut_id;
var filter_var_id_hecomp;
function onRowClick_details(value) {
	var gridLog = dijit.byId("gridLog");
	var gridHe = dijit.byId("gridHeComp");
	gridLog.filter({projName: "None"});
	var gridRes = dijit.byId("gridId");
	var item = gridRes.getItem(value.rowIndex);
	var val_proj = gridRes.store.getValue(item, "nbProj");
	var val_varId = gridRes.store.getValue(item, "varId");
	var val_cons = gridRes.store.getValue(item, "consequence");
	alamut_id = gridRes.store.getValue(item, "alamut");
	var val_varId_toPrint = val_varId;
	if (val_varId_toPrint.match(/;/)) {
		var tab_tmp = val_varId.split(";");
		val_varId_toPrint = tab_tmp[1];
	}
	varId_selected = val_varId_toPrint;
	var toFilter = val_proj.replace(/; /g, "|");
	var id = details(value, 'variant');
	var gridDetails = dijit.byId("gridDetails");
	var build = 'HG19';
	if (document.getElementById("aff_build")) { build = document.getElementById("aff_build").innerHTML; }
	var urlLaunch_dejavu = urlDejaVu + "?login=" + dojo.cookie("username") + "&pwd=" + dojo.cookie("passwd") + "&input=" + "id:" + id + "&details=true&build=" + build;
	if (force_annot_version) {
		urlLaunch_dejavu += "&force_annot_version="+force_annot_version;
	}
	else if (param("origin_project")) {
		urlLaunch_dejavu += "&origin_project="+param("origin_project");
	}
	dataStore_details = new dojo.data.ItemFileReadStore({ url: urlLaunch_dejavu });
	gridDetails.setStore(dataStore_details);
	gridDetails.store.fetch({
		onComplete: function(items){
			doSearch_he('clear');
			valId_forHe = val_varId_toPrint;
			filter_var_id_hecomp = valId_forHe;
			getFiltersForVarSelect_HeComp(val_cons);
			gridLog.filter({projName: "*"});
			gridLog.filter({projName: new RegExp(toFilter)}, true);
			gridHe.filter({id: new RegExp(val_varId)}, true);
			document.getElementById("aff_he").innerHTML = "VARIATION " + val_varId_toPrint;
			printNbResults_dejavu_he();
			gridDetails.filter({transcript: "*"});
			dijit.byId("panel_details").setTitle("<b>" + panelDetails_defaultTitle + "<FONT COLOR=RED> -> for " + val_varId_toPrint + "</FONT></b>");
			dijit.byId("panel_he").setTitle("<b>" + panelHeComp_defaultTitle + "<FONT COLOR=RED> -> for " + val_varId_toPrint + "</FONT></b>");
			gridRes.getRowNode(value.rowIndex).style.backgroundColor = "red"; 
		}
	});
}

function onRowClick_varInfos(value) {
	openPopup();
	dijit.byId('waiting').show();
	var gridLog = dijit.byId("gridLog");
	var gridHe = dijit.byId("gridHeComp");
	gridLog.filter({projName: "None"});
	var gridRes = dijit.byId("gridId");
	var item = gridRes.getItem(value.rowIndex);
	var val_genes = gridRes.store.getValue(item, "gene");
	document.getElementById("align_gene").innerHTML = "<b><u>Gene:</b></u> " + val_genes;
	var val_proj = gridRes.store.getValue(item, "nbProj");
	var val_varId = gridRes.store.getValue(item, "varId");
	var val_cons = gridRes.store.getValue(item, "consequence");
	var val_varId_toPrint = val_varId;
	if (val_varId_toPrint.match(/;/)) {
		var tab_tmp = val_varId.split(";");
		val_varId_toPrint = tab_tmp[1];
	}
	varId_selected = val_varId_toPrint;
	var toFilter = val_proj.replace(/; /g, "|");
	var id = details(value, 'variant');
	var gridInfos = dijit.byId("gridInfos");
	var build = 'HG38';
	if (document.getElementById("aff_build")) { build = document.getElementById("aff_build").innerHTML; }
	var urlLaunch_dejavu = url_varInfos + "?login=" + dojo.cookie("username") + "&pwd=" + dojo.cookie("passwd") + "&input=" + id + "&build=" + build;
	dataStore_details = new dojo.data.ItemFileReadStore({ url: urlLaunch_dejavu });
	gridInfos.setStore(dataStore_details);
	gridInfos.store.fetch({
		onComplete: function(items){
			gridInfos.filter({name: "*"});
			gridRes.getRowNode(value.rowIndex).style.backgroundColor = "red";
			checkDepthAndAuthorizedPatientsFromVarInfos(items);
			dijit.byId('waiting').hide();
		}
	});
}

function checkDepthAndAuthorizedPatientsFromVarInfos(items) {
	var gridInfos = dijit.byId("gridInfos");
	var implicated = 0
	var not_implicated = 0
	var nb_infos = 0;
	var sum_depth = 0;
	var min_depth = 999999999;
	var max_depth = 0;
	dojo.forEach(items, function(item){
		dojo.forEach(gridInfos.store.getAttributes(item), function(attribute) {
			if (attribute == "implicated") {
	        	var value = gridInfos.store.getValues(item, attribute);
				if (value == 'X') { implicated += 1; }
				else { not_implicated += 1; }
			}
			if (attribute == "dp") {
	        	var value = parseFloat(gridInfos.store.getValues(item, attribute));
				if (value > 0) {
					nb_infos += 1;
					sum_depth += value;
					if (value < min_depth) { min_depth = value; }
					if (value > max_depth) { max_depth = value; }
				}
			}
		});
	});
	var mean_depth = sum_depth / parseFloat(nb_infos);
	document.getElementById("align_min").innerHTML = "<b><u>Min DP:</b></u> " + min_depth;
	document.getElementById("align_max").innerHTML = "<b><u>Max DP:</b></u> " + max_depth;
	document.getElementById("align_mean").innerHTML = "<b><u>Mean DP:</b></u> " + mean_depth;
	document.getElementById("details_patOk").innerHTML = " <b><u>Authorized Patient(s): </b></u>" + implicated;
	document.getElementById("details_patNok").innerHTML = " <b><u>Not Authorized Patient(s): </b></u>" + not_implicated;
	document.getElementById("contacts").innerHTML = " <b><u>Contact(s): </b></u> first, click to a patient row...";
}

function selectPatient_dejavu(value) {
	var gridInfos = dijit.byId("gridInfos");
	var item = gridInfos.getItem(value.rowIndex)
	var mails = gridInfos.store.getValue(item, "contacts");
	if (mails == 'NA') { document.getElementById("contacts").innerHTML = " <b><u>Contact(s): </b></u> Not authorized... Please ask to BioInfo TEAM (marc.bras\@institutimagine.org  or  cecile.masson\@gmail.com)"; }
	else { document.getElementById("contacts").innerHTML = " <b><u>Contact(s): </b></u> " + mails; }
	selectPatientGlobal(value, "gridInfos");
}

function onRowDblClick_variant(value) {
	document.getElementById("popup_text").innerHTML = "<u><b><FONT COLOR=red>Variant details</FONT><b></u>";
	openPopup();
	dijit.byId("gridLog").render();
}

function onRowDblClick_hecomp(value) {
	document.getElementById("popup_text").innerHTML = "<u><b><FONT COLOR=red>Heterozygote composite details</FONT><b></u>";
	openPopupLog();
	var gridLog = dijit.byId("gridLog");
	var gridHeComp = dijit.byId("gridHeComp");
	var item = gridHeComp.getItem(value.rowIndex);
	var val_proj = gridHeComp.store.getValue(item, "nbProj");
	var toFilter = val_proj.replace(/; /g, "|");
	var urlToLaunch = url_checkLogin + "?login=" + dojo.cookie("username") + "&pwd=" + dojo.cookie("passwd");
	dataStore_log = new dojo.data.ItemFileReadStore({ url: urlToLaunch });
	gridLog.setStore(dataStore_log);
	gridLog.store.fetch({
		onComplete: function(items){
			details(value, 'hecomp');
			gridLog.filter({projName: new RegExp(toFilter)}, true);
			gridLog.render();
		}
	});
}

var filter_var_ids;
function onRowClick_HeComp(value) {
	var gridHe = dijit.byId("gridHeComp");
	var item = gridHe.getItem(value.rowIndex);
	var varId_a = gridHe.store.getValue(item, "var_1");
	var varId_b = gridHe.store.getValue(item, "var_2");
	filter_var_ids = varId_a + ',' +varId_b;
	var id = varId_a + ',' + varId_b
	var toFilter = varId_a + "|" + varId_b;
	var toWrite = varId_a + " and " + varId_b;
	var gridRes = dijit.byId("gridId");
	gridRes.filter({varId: new RegExp(toFilter)}, true);
	printNbResults_dejavu();
	dijit.byId("panel_results").setTitle("<b>" + panelResults_defaultTitle + "     <FONT COLOR=RED>-> only " + toWrite + "</FONT></b>");
}

function details(value, type) {
	var varId = '';
	var rsName = '';
	var details = new Object();
	var patients = new Object();
	var grid;
	if (type == 'variant') { grid = dijit.byId("gridId"); }
	else if (type == 'hecomp') { grid = dijit.byId("gridHeComp"); }
	var items = grid.getItem(value.rowIndex);
	dojo.forEach(grid.store.getAttributes(items), function(attribute) {
        var value = grid.store.getValues(items, attribute);
		if (value == "") { value = "N.A."; }
		if (attribute == "nbProj") {
			var tmp = new String(value);
			var tab = tmp.split("; ");
			var ok = 0;
			var nok = 0;
			dojo.forEach(tab, function(projName){
				if (projName == "X") { nok += 1; }
				else { ok += 1; }
			});
			document.getElementById("details_projOk").innerHTML = "<b><u>Authorized Project:</b></u> " + ok;
			document.getElementById("details_projNok").innerHTML = "<b><u>Not Authorized Project:</b></u> " + nok;
		}
		if (attribute == "nbPat") {
			var tmp = new String(value);
			var tab = tmp.split("|");
			tmp = tab[1];
			tab = tmp.split("; ");
			dojo.forEach(tab, function(patName){
				var tmp2 = patName.split(" [");
				var project = tmp2[0];
				var pat = tmp2[1].replace("]", "");
				if (dijit.byId("filter_strict_ho").checked || dijit.byId("filter_ho").checked) {
					document.getElementById("details_affAlertHo").innerHTML = "<u><b><FONT COLOR=red>Be careful, only Ho patients given !</FONT><b></u>";
					if (/\(Ho\)/.test(pat)) {
						var patHo = [];
						var tmp3 = pat.split(", ");
						dojo.forEach(tmp3, function(thisPatName){
							if (/\(Ho\)/.test(thisPatName)) {
								thisPatName = thisPatName.replace("(Ho)", "");
								patHo.push(thisPatName);
							}
						});
						hashPatient[project] = patHo.join(', ');
					}
					else { hashPatient[project] = ''; }
				}
				else {
					document.getElementById("details_affAlertHo").innerHTML = "";
					pat = pat.replace("(Ho)", "");
					hashPatient[project] = pat;
				}
			});
		}
		details[attribute] = value;
    });
    if (type == 'variant') {
		document.getElementById("details_id").innerHTML = "<b><u>Id:</b></u> " + details["id"];
		document.getElementById("details_chr").innerHTML = "<b><u>Chr:</b></u> " + details["chr"];
		document.getElementById("details_pos").innerHTML = "<b><u>Position:</b></u> " + details["pos"];
		document.getElementById("details_rsname").innerHTML = "<b><u>rsName: </b></u> " + details["rsName"];
		document.getElementById("details_freq").innerHTML = "<b><u>dbSNP freq:</b></u> " + details["freqDbSnp"];
		document.getElementById("align_id").innerHTML = "<b><u>Id:</b></u> " + details["id"];
		document.getElementById("align_chr").innerHTML = "<b><u>Chr:</b></u> " + details["chr"];
		document.getElementById("align_pos").innerHTML = "<b><u>Position:</b></u> " + details["pos"];
		document.getElementById("align_rsname").innerHTML = "<b><u>rsName: </b></u> " + details["rsName"];
		document.getElementById("align_freq").innerHTML = "<b><u>dbSNP freq:</b></u> " + details["freqDbSnp"];
		return details["varId"];
    }
}

var he_comp_launched = 0;
var last_url_he_comp;
var datastore_gridHeComp_origin;
function launchHeComp(value, gridName, projectName) {
	if (he_comp_launched == 1) { return; }
	he_comp_launched = 1;
	if (value) {  }
	else { value = varId_now; }
	if (value) {  }
	else { return; }
	if (document.getElementById("aff_build")) { document.getElementById("aff_nbhe_dejavu").innerHTML = ''; }
	var user = dojo.cookie("username");
	var pass = dojo.cookie("passwd");
	var gridHeComp;
	if (gridName) { gridHeComp = dijit.byId(gridName); }
	else { gridHeComp = dijit.byId("gridHeComp"); }
	var urlLaunch_HeComp = urlHeComp + "?login=" + user + "&pwd=" + pass + "&input=" + "id:" + value;
	if (dijit.byId("filter_ho") && dijit.byId("filter_strict_ho")) {
		if (dijit.byId("filter_ho").checked || dijit.byId("filter_strict_ho").checked) {
			urlLaunch_HeComp = urlLaunch_HeComp + '&onlyho=1';
		}
		if (dijit.byId("show_large_deletion").checked) {
			urlLaunch_HeComp = urlLaunch_HeComp + '&show_large_del=1';
		}
	}
	if (projectName) { urlLaunch_HeComp += "&thisproject=" + projectName; }
	if (force_annot_version) {
		urlLaunch_HeComp += "&force_annot_version="+force_annot_version;
	}
	else if (param("origin_project")) {
		urlLaunch_HeComp += "&origin_project="+param("origin_project");
	}
	var build = 'HG19';
	if (document.getElementById("aff_build")) { build = document.getElementById("aff_build").innerHTML; }
	urlLaunch_HeComp = urlLaunch_HeComp + "&build=" + build;
	last_url_he_comp = urlLaunch_HeComp;
	dataStore_HeComp = new dojo.data.ItemFileWriteStore({ url: urlLaunch_HeComp });
	gridHeComp.setStore(dataStore_HeComp);
	gridHeComp.store.fetch({
		onComplete: function(items){
			document.getElementById("aff_he").innerHTML = originLaunchInput_dejavu;
			dijit.byId("panel_he").setTitle("<b>" + panelHeComp_defaultTitle + "</b>");
			valId_forHe = false;
			gridHeComp.beginUpdate();
			gridHeComp.filter({id: "*"});
			printNbResults_dejavu_he();
			gridHeComp.endUpdate();
			enabledXlsExport();
			datastore_gridHeComp_origin = dataStore_HeComp;
			dijit.byId('waiting').hide();
		}
	});
}

var varId_now;
var last_url;
var datastore_grid_origin;
function launchId(value) {
	if (LoginOk == true) {
		if (to_do_stats) { project_stats() }
		filter_var_ids = '';
		filter_var_id_hecomp = '';
		filters_sift = [];
		filters_polyphen = [];
		filters_consequence = [];
		doSearch('clear');
		doSearch_he('clear');
		disabledXlsExport();
		dijit.byId('waiting').show();
		dijit.byId("panel_graphics").set('open',false);
		dijit.byId("panel_details").set('open',false);
		dijit.byId("panel_he").set('open',false);
		dijit.byId("panel_results").set('open',true);
		dijit.byId("id").set("disabled", true);
		var gridHe = dijit.byId("gridHeComp");
		var grid = dijit.byId("gridId");
		var user = dojo.cookie("username");
		var pass = dojo.cookie("passwd");
		var varId_input = value;
		varId_now = value;
		//launchHeComp(varId_input);
		if (document.title.match('Deja Vu')) { document.title = "Deja Vu [" + varId_input + "]"; }
		var urlLaunch_dejavu = urlDejaVu + "?login=" + user + "&pwd=" + pass + "&input=" + "id:" + varId_input;
//		if (dijit.byId("filter_only_my_projects").checked) {
//			urlLaunch_dejavu = urlLaunch_dejavu + "&onlyMyProjects=1";
//		}
		if (dijit.byId("filter_ho").checked) {
			dijit.byId("filter_strict_ho").set('checked',false);
			urlLaunch_dejavu = urlLaunch_dejavu + "&onlyHo=1";
		}
		if (dijit.byId("filter_strict_ho").checked) {
			dijit.byId("filter_ho").set('checked',false);
			urlLaunch_dejavu = urlLaunch_dejavu + "&onlyStrictHo=1";
		}
		if (dijit.byId("show_large_deletion").checked) {
			urlLaunch_dejavu = urlLaunch_dejavu + '&show_large_del=1';
		}
		var build = 'HG19';
		if (document.getElementById("aff_build")) { build = document.getElementById("aff_build").innerHTML; }
		urlLaunch_dejavu = urlLaunch_dejavu + "&build=" + build;
		urlLaunch_dejavu += get_polyquery_var_annotation_checked();
		urlLaunch_dejavu += get_value_slider_freq();
		if (force_annot_version) {
			urlLaunch_dejavu += "&force_annot_version="+force_annot_version;
		}
		else if (param("origin_project")) {
			urlLaunch_dejavu += "&origin_project="+param("origin_project");
		}
		var dejavu_value = get_value_slider_dejavu();
		if (dejavu_value == 100) {}
		else if (dejavu_value == 0) {
			dijit.byId('waiting').hide();
			alert('NO Result possible with DejaVu = 0');
			return;
		}
		else {
			urlLaunch_dejavu += "&max_dejavu=" + dejavu_value;
		}
		last_url = urlLaunch_dejavu;
		dataStore = new dojo.data.ItemFileReadStore({ url: urlLaunch_dejavu });
		grid.errorMessage = "No result for " + varId_input + "... (Annotation Version used: " + document.getElementById("aff_annot_version").innerHTML + ")<br><br>Please, reload page for a new research."; 
		grid.noDataMessage = "No result for " + varId_input + "..."; 
		grid.setStore(dataStore);
		grid.store.fetch({
			onComplete: function(items){
				var type;
				if ( /ENST/.test(varId_input) ) { type = "TRANSCRIPT "; }
				else if ( /ENSG/.test(varId_input) ) { type = "GENE "; }
				else if ( /,/.test(varId_input) ) { type = "VARIANTS LIST"; }
				else if ( /_/.test(varId_input) ) { type = "VARIANT "; }
				else if ( /rs/.test(varId_input) ) { type = "VARIANT "; }
				else if ( /:/.test(varId_input) ) { type = "REGION "; }
				else { type = "GENE "; }
				var inputText = type;
				if (type == "VARIANTS LIST") { originLaunchInput_dejavu = type; }
				else {
					originLaunchInput_dejavu = type + varId_input;
					inputText += varId_input
				}
				document.getElementById("aff_log").innerHTML = dojo.cookie("username");
				drawChart_consequence(originLaunchInput_dejavu);
				drawChart_polyphen(originLaunchInput_dejavu);
				drawChart_sift(originLaunchInput_dejavu);
				document.getElementById("affichage").innerHTML = originLaunchInput_dejavu;
				grid.filter({varId: "*"});
				filters_consequence = [];
				filters_cons_hecomp_1 = '';
				filters_cons_hecomp_2 = '';
				doSearch_he('clear');
				refresh();
				printNbResults_dejavu();
				datastore_grid_origin = dataStore;
				dijit.byId('waiting').hide();
				dijit.byId("id").set("disabled", false);
    			dijit.byId("b_export_xls_results").set("disabled", false);
				document.getElementById("aff_nbhe_dejavu").innerHTML = 'Click HeComp tab to run analysis';
				he_comp_launched = 0;
			}
		});
	}
	else { alert("LOGIN FAILED... Reload page please to retry !"); }
}

function updateStatusButons(filterName) {
	if (filterName == "filter_ho") {
		dijit.byId("filter_strict_ho").set('checked',false);
	}
	if (filterName == "filter_strict_ho") {
		dijit.byId("filter_ho").set('checked',false);
	}
} 

function printNbResults_dejavu(){
	var grid = dijit.byId("gridId");
	nbResults = grid.rowCount;
	if (nbResults == 2500) { alert ("WARNING: Only printed first 2500 results !"); }
	document.getElementById("aff_nbres_dejavu").innerHTML = nbResults;
}

function printNbResults_dejavu_he(){
	var grid = dijit.byId("gridHeComp");
	if (dijit.byId("filter_ho").checked) {
		document.getElementById("aff_nbhe_dejavu").innerHTML = "0 (Ho analyse...)";
		dijit.byId("b_export_xls_he_comp").set("disabled", true);
	}
	else if (dijit.byId("filter_strict_ho").checked) {
		document.getElementById("aff_nbhe_dejavu").innerHTML = "0 (Strict Ho analyse...)";
		dijit.byId("b_export_xls_he_comp").set("disabled", true);
	}
	else {
		nbResults = grid.rowCount;
		document.getElementById("aff_nbhe_dejavu").innerHTML = nbResults;
	}
}

function checkAndOrButton(value) {
	if (value == "and") {
		if (dijit.byId("b_launch_and").checked) { dijit.byId("b_launch_or").set('checked',false); }
		else { dijit.byId("b_launch_or").set('checked',true); }
	}
	else if (value == "or") {
		if (dijit.byId("b_launch_or").checked) { dijit.byId("b_launch_and").set('checked',false); }
		else { dijit.byId("b_launch_and").set('checked',true); }
	}
	doSearch("run");
}

function doSearch(type){
	var grid = dijit.byId("gridId");		
	grid.setQuery();
	var isCleaning = false; 
	if (type == "clear") {
		grid.filter({varId: "*"});
		dijit.byId("gridHeComp").filter({id: "*"});
		dijit.byId("gridLog").filter({projName: "None"});
		dijit.byId("gridDetails").filter({transcript: "None"});
		printNbResults_dejavu_he();
		document.getElementById("details_id").innerHTML = "<b><u>Id:</b></u> " + details["varId"];
		document.getElementById("details_chr").innerHTML = "<b><u>Chr:</b></u> " + details["chr"];
		document.getElementById("details_pos").innerHTML = "<b><u>Position:</b></u> " + details["pos"];
		document.getElementById("details_rsname").innerHTML = "<b><u>rsName: </b></u> " + details["rsName"];
		document.getElementById("details_freq").innerHTML = "<b><u>dbSNP freq:</b></u> " + details["freqDbSnp"];
		document.getElementById("aff_he").innerHTML = originLaunchInput_dejavu;
		dijit.byId("panel_details").set('open',false);
		dijit.byId("panel_he").set('open',false);
		isCleaning = true;
		type = "run";
	}
	else if (type == "text"){			
		var field = dijit.byId("search_field");
		var text_search = dijit.byId("search_text");
		if (field.getValue() == "varId") { grid.filter({varId: new RegExp(text_search.getValue())}); }
		else if (field.getValue() == "gene") { grid.filter({gene: new RegExp(text_search.getValue())}); }
		else if (field.getValue() == "transcript") { grid.filter({transcript: new RegExp(text_search.getValue())}); }
		else if (field.getValue() == "rsName") { grid.filter({rsName: new RegExp(text_search.getValue())}); }
		else if (field.getValue() == "patient") { grid.filter({nbPat: new RegExp(text_search.getValue())}); }
	}
	if (type == "run") {		
		if (isCleaning) {
			dijit.byId("panel_graphics").setTitle("<b>" + panelFilters_defaultTitle + "</b>");
			dijit.byId("panel_he").setTitle("<b>" + panelHeComp_defaultTitle + "</b>");
			dijit.byId("panel_details").setTitle("<b>" + panelDetails_defaultTitle + "</b>");
			dijit.byId("panel_results").setTitle("<b>" + panelResults_defaultTitle + "</b>");
		}
		else {
			var regexpFilters;
			var i = 0;
			var j = 0;
			var listRegexpFilters = [];
			if (dijit.byId("b_launch_ho").checked) { listRegexpFilters.push('ho'); i++; }
			if (dijit.byId("b_launch_he").checked) { listRegexpFilters.push('he'); i++; }
			if (dijit.byId("b_launch_strictho").checked) { listRegexpFilters.push('strictHo'); i++; }
			if (dijit.byId("b_launch_stricthe").checked) { listRegexpFilters.push('strictHe'); i++; }
			if (dijit.byId("b_launch_rsname").checked) { listRegexpFilters.push('rsname'); i++; }
			if (dijit.byId("b_launch_norsname").checked) { listRegexpFilters.push('notRsname'); i++; }
			if (dijit.byId("b_launch_implicated").checked) { listRegexpFilters.push('myVar'); i++; }
			if (dijit.byId("b_launch_notmyvar").checked) { listRegexpFilters.push('notMyVar'); i++; }
			if (dijit.byId("filter_dbSnp_l1p").checked) { listRegexpFilters.push('dbSnpLess1%'); j++; }
			if (dijit.byId("filter_dbSnp_l25p").checked) { listRegexpFilters.push('dbSnpLess25%'); j++; }
			if (dijit.byId("filter_dbSnp_l50p").checked) { listRegexpFilters.push('dbSnpLess50%'); j++; }
			if (dijit.byId("filter_dbSnp_l75p").checked) { listRegexpFilters.push('dbSnpLess75%'); j++; }
			if (dijit.byId("filter_dbSnp_g1p").checked) { listRegexpFilters.push('dbSnpSup1%'); j++; }
			if (dijit.byId("filter_dbSnp_g25p").checked) { listRegexpFilters.push('dbSnpSup25%'); j++; }
			if (dijit.byId("filter_dbSnp_g50p").checked) { listRegexpFilters.push('dbSnpSup50%'); j++; }
			if (dijit.byId("filter_dbSnp_g75p").checked) { listRegexpFilters.push('dbSnpSup75%'); j++; }
			listRegexpFilters = listRegexpFilters.sort();
			regexpToWrite = listRegexpFilters.join(', ');
			if (dijit.byId("b_launch_or").checked) {
				regexpToWrite += " with 'OR'";
				regexpFilters = listRegexpFilters.join('|');
				if (regexpFilters == "") { regexpFilters = "xyz"; }
			}
			else if (dijit.byId("b_launch_and").checked) {
				regexpToWrite += " with 'AND'";
				regexpFilters = listRegexpFilters.join('.+');
			}
			if (i == 8 && j == 0 && dijit.byId("b_launch_or").checked) {
				dijit.byId("panel_graphics").setTitle("<b>" + panelFilters_defaultTitle + "</b>");
				dijit.byId("panel_results").setTitle("<b>" + panelResults_defaultTitle + "</b>");
			}
			else if (listRegexpFilters.length == 0) {
				dijit.byId("panel_graphics").setTitle("<b>" + panelFilters_defaultTitle + "<FONT COLOR=RED> -> None...</FONT></b></b>");
				dijit.byId("panel_results").setTitle("<b>" + panelResults_defaultTitle + "<FONT COLOR=RED> -> None...</FONT></b>");
			}
			else {
				dijit.byId("panel_graphics").setTitle("<b>" + panelFilters_defaultTitle + "<FONT COLOR=RED> -> Filters: " + regexpToWrite + " !</FONT></b>");
				dijit.byId("panel_results").setTitle("<b>" + panelResults_defaultTitle + "<FONT COLOR=RED> -> Filters: " + regexpToWrite + " !</FONT></b>");
			}
			grid.filter({filters: new RegExp(regexpFilters)});
		}
	}
	printNbResults_dejavu();	
}
        
function viewAlamut() {
    displayInAlamutVector(alamut_id);
}

var filters_cons_hecomp_1;
var filters_cons_hecomp_2;
function doSearch_he(type){
	var grid = dijit.byId("gridHeComp");
	if (type == "run") {
		filters_cons_hecomp_1 = '';
		filters_cons_hecomp_2 = '';
		var listRegexpFilters = [];
		var listRegexpFilters_2 = [];
		if (dijit.byId("filter_2_intergenic").checked) { listRegexpFilters.push('Intergenic'); }
        if (dijit.byId("filter_2_nocoding").checked) { listRegexpFilters.push('Intronic'); }
        if (dijit.byId("filter_2_pseudogene").checked) { listRegexpFilters.push('Pseudogene'); }
        if (dijit.byId("filter_2_ncRNA").checked) { listRegexpFilters.push('ncRNA'); }
        if (dijit.byId("filter_2_mature").checked) { listRegexpFilters.push('Mature miRNA'); }
        if (dijit.byId("filter_2_utr").checked) { listRegexpFilters.push('Utr'); }
        if (dijit.byId("filter_2_splicing").checked) { listRegexpFilters.push('Splice Region'); }
        if (dijit.byId("filter_2_essential_splicing").checked) { listRegexpFilters.push('Splice Acc/Don'); }
        if (dijit.byId("filter_2_coding").checked) { listRegexpFilters.push('Missense'); }
        if (dijit.byId("filter_2_frameshift").checked) { listRegexpFilters.push('Frameshift'); }
        if (dijit.byId("filter_2_stop").checked) { listRegexpFilters.push('Stop-gained'); }
        if (dijit.byId("filter_2_phase").checked) { listRegexpFilters.push('\(Start/Stop\)-lost'); }
        if (dijit.byId("filter_2_silent").checked) { listRegexpFilters.push('Synonymous'); }
        if (dijit.byId("filter_2_non_frameshift").checked) { listRegexpFilters.push('No-frameshift'); }
		var filter = listRegexpFilters.join('|');
		filters_cons_hecomp_1 = listRegexpFilters.join(',');
		if (dijit.byId("filter_2_intergenic_2").checked) { listRegexpFilters_2.push('Intergenic'); }
        if (dijit.byId("filter_2_nocoding_2").checked) { listRegexpFilters_2.push('Intronic'); }
        if (dijit.byId("filter_2_pseudogene_2").checked) { listRegexpFilters_2.push('Pseudogene'); }
        if (dijit.byId("filter_2_ncRNA_2").checked) { listRegexpFilters_2.push('ncRNA'); }
        if (dijit.byId("filter_2_mature_2").checked) { listRegexpFilters_2.push('Mature miRNA'); }
        if (dijit.byId("filter_2_utr_2").checked) { listRegexpFilters_2.push('Utr'); }
        if (dijit.byId("filter_2_splicing_2").checked) { listRegexpFilters_2.push('Splice Region'); }
        if (dijit.byId("filter_2_essential_splicing_2").checked) { listRegexpFilters_2.push('Splice Acc/Don'); }
        if (dijit.byId("filter_2_coding_2").checked) { listRegexpFilters_2.push('Missense'); }
        if (dijit.byId("filter_2_frameshift_2").checked) { listRegexpFilters_2.push('Frameshift'); }
        if (dijit.byId("filter_2_stop_2").checked) { listRegexpFilters_2.push('Stop-gained'); }
        if (dijit.byId("filter_2_phase_2").checked) { listRegexpFilters_2.push('\(Start/Stop\)-lost'); }
        if (dijit.byId("filter_2_silent_2").checked) { listRegexpFilters_2.push('Synonymous'); }
        if (dijit.byId("filter_2_non_frameshift_2").checked) { listRegexpFilters_2.push('No-frameshift'); }
		var filter_2 = listRegexpFilters_2.join('|');
		filters_cons_hecomp_2 = listRegexpFilters_2.join(',');
		var hash_cons1 = new Object();
		var hash_cons2 = new Object();
		dojo.forEach(listRegexpFilters, function(filter_name) {
			hash_cons1[filter_name] = 1;
		});
		dojo.forEach(listRegexpFilters_2, function(filter_name) {
			hash_cons2[filter_name] = 1;
		});
		if (document.getElementById("aff_nbhe_dejavu").innerHTML == "Click HeComp tab to run analysis") { }
		else if (document.getElementById("aff_nbhe_dejavu").innerHTML == "0") { }
		else {
			dijit.byId("gridHeComp").setStore( datastore_gridHeComp_origin );
			var grid = dijit.byId("gridHeComp");
		    grid.store.fetch({
		        onComplete: function(items, request){
					var n = 0;
					var max = grid.rowCount;
					while (n < max) {
				    	var selectedItem = items[n];
				    	grid.store.setValue(selectedItem, 'pass_filter', '0');
				        var value = grid.store.getValue(selectedItem, 'allConsequences');
						var list_values = value.split('|');
						var lFiltersWithRs = list_values[0].split(';');
						var lFilters = lFiltersWithRs[0].split('_');
						lFilters.pop();
						lFilters.pop();
						lFilters.pop();
						lFilters.pop();
						var lVarFilters1 = [];
						var lVarFilters2 = [];
						dojo.forEach(lFilters, function(filter_name) {
							var lTmp = filter_name.split('.');
							if (lTmp[1] == 1) { lVarFilters1.push(lTmp[0]); }
							else { lVarFilters2.push(lTmp[0]); }
						});
						var is_ok_1_1 = 0;
						var is_ok_1_2 = 0;
						var is_ok_2_1 = 0;
						var is_ok_2_2 = 0;
						dojo.forEach(lVarFilters1, function(filter_name) {
							if (hash_cons1.hasOwnProperty(filter_name)) { is_ok_1_1 = 1; }
							if (hash_cons2.hasOwnProperty(filter_name)) { is_ok_2_1 = 1; }
						});
						dojo.forEach(lVarFilters2, function(filter_name) {
							if (hash_cons1.hasOwnProperty(filter_name)) { is_ok_2_2 = 1; }
							if (hash_cons2.hasOwnProperty(filter_name)) { is_ok_1_2 = 1; }
						});
						if ((is_ok_1_1 + is_ok_1_2) == 2 || (is_ok_2_1 + is_ok_2_2) == 2) {
		        			grid.store.setValue(selectedItem, 'pass_filter', '1');
						}
						else {
		        			grid.store.setValue(selectedItem, 'pass_filter', '0');
						}
						n++;
					}
		        }
		    });
			dijit.byId("gridHeComp").filter({pass_filter: '1'});
		}
		refresh();
		if (valId_forHe) { dijit.byId("panel_he").setTitle("<b>" + panelHeComp_defaultTitle + "<FONT COLOR=RED> -> for <u>" + valId_forHe  + "</u> with filter(s) used !</FONT></b>"); }
		else { dijit.byId("panel_he").setTitle("<b>" + panelHeComp_defaultTitle + "<FONT COLOR=RED> -> filter(s) used !</FONT></b>"); }
	}
    else if (type == "clear") {
		dijit.byId("gridHeComp").filter({allConsequences: '*'});
		dijit.byId("panel_he").setTitle("<b>" + panelHeComp_defaultTitle + "</b>");
		document.getElementById("aff_he").innerHTML = originLaunchInput_dejavu;
		valId_forHe = false;
		getFiltersForVarSelect_HeComp('clear');
	}
	printNbResults_dejavu_he();	
}

function getFiltersForVarSelect_HeComp(value) {
	var tab_idButton = ["intergenic", "nocoding", "pseudogene", "ncRNA", "mature", "utr", "splicing", "essential_splicing", "coding", "frameshift", "stop", "phase", "silent", "non_frameshift"];
	if (value == "clear") {
		for (i=0; i < tab_idButton.length; i++) {
			dijit.byId("filter_2_" + tab_idButton[i]).set('checked', true);
			dijit.byId("filter_2_" + tab_idButton[i]).set('disabled', false);
			dijit.byId("filter_2_" + tab_idButton[i] + '_2').set('checked', true);
			dijit.byId("filter_2_" + tab_idButton[i] + '_2').set('disabled', false);
		}
		document.getElementById("nameA").innerHTML = '';
	}
	else {
		for (i=0; i < tab_idButton.length; i++) {
			var valToCompare = tab_idButton[i];
            dijit.byId("filter_2_" + valToCompare).set('checked', false);
            dijit.byId("filter_2_" + valToCompare).set('disabled', true);
			dijit.byId("filter_2_" + tab_idButton[i] + '_2').set('checked', true);
			dijit.byId("filter_2_" + tab_idButton[i] + '_2').set('disabled', false);
		}
		if (value.match(/[Ii]ntronic/)) { dijit.byId("filter_2_nocoding").set('checked', true); }
		if (value.match(/[Ii]ntergenic/)) { dijit.byId("filter_2_intergenic").set('checked', true); }
        if (value.match(/[Pp]seudogene/)) { dijit.byId("filter_2_pseudogene").set('checked', true); }
        if (value.match(/nc[Rr][Nn][Aa]/)) { dijit.byId("filter_2_ncRNA").set('checked', true); }
        if (value.match(/mature/)) { dijit.byId("filter_2_mature").set('checked', true); }
        if (value.match(/[Uu]tr/)) { dijit.byId("filter_2_utr").set('checked', true); }
        if ((value == "splicing") || (value = "Splice Region")) { dijit.byId("filter_2_splicing").set('checked', true); }
        if (value.match(/essential_splicing/) || value.match(/Splice Acc\/Don/)) { dijit.byId("filter_2_essential_splicing").set('checked', true); }
        if ((value.match(/coding/) )|| (value.match(/Missense/))) { dijit.byId("filter_2_coding").set('checked', true); }
        if ((value == "frameshift") || (value = "Frameshift")) { dijit.byId("filter_2_frameshift").set('checked', true); }
        if (value.match(/Stop-gained/)) { dijit.byId("filter_2_stop").set('checked', true); }
        if (value.match(/phase/)) { dijit.byId("filter_2_phase").set('checked', true); }
        if (value.match(/Start/)) { dijit.byId("filter_2_silent").set('checked', true); }
        if ((value.match(/non-frameshift/)) || (value.match(/No-frameshift/))) { dijit.byId("filter_2_non_frameshift").set('checked', true); }
		document.getElementById("nameA").innerHTML =  "Variant A is <FONT COLOR=red><b><u>" + valId_forHe + "</b></u></FONT>";
	}
}

function showAll(value) {
	filter_var_ids = '';
	filter_var_id_hecomp = '';
	filters_sift = [];
	filters_polyphen = [];
	filters_consequence = [];
	dijit.byId('b_launch_clear_he').set('disabled', false);
	dijit.byId('b_launch_run_he').set('disabled', false);
	dijit.byId("panel_graphics").setTitle("<b>" + panelFilters_defaultTitle + "</b>");
	dijit.byId("panel_results").setTitle("<b>" + panelResults_defaultTitle + "</b>");
	dijit.byId("gridId").setStore( datastore_grid_origin );
	dijit.byId("gridHeComp").setStore( datastore_gridHeComp_origin );
    drawChart_consequence();
    drawChart_polyphen();
    drawChart_sift();
	doSearch('clear');
	printNbResults_dejavu();
	printNbResults_dejavu_he();
	var listCons = ['filter_intergenic', 'filter_nocoding', 'filter_pseudogene', 'filter_ncRNA', 'filter_mature', 'filter_utr', 'filter_splicing', 'filter_essential_splicing', 'filter_coding', 'filter_frameshift', 'filter_stop', 'filter_phase', 'filter_silent', 'filter_non_frameshift'];
	dojo.forEach(listCons, function(filter_name) {
		dijit.byId(filter_name).set('checked', true);
		dijit.byId(filter_name+'_2').set('checked', true);
		dijit.byId(filter_name).set('disabled', false);
		dijit.byId(filter_name+'_2').set('disabled', false);
	});
	dijit.byId('b_launch_clear_he').set('disabled', false);
	dijit.byId('b_launch_run_he').set('disabled', false);
	filters_cons_hecomp_1 = '';
	filters_cons_hecomp_2 = '';
	refresh();
}

function refresh(value) {
	dijit.byId("gridId").render();
	dijit.byId("gridDetails").render();
	dijit.byId("gridHeComp").render();
	dijit.byId("gridInfos").render();
}

function search_project(value) {
	var grid = dijit.byId("gridId");
	grid.setQuery();
	grid.filter({nbProj: new RegExp(value)});
	printNbResults_dejavu();
}

var filters_consequence = [];
function filter_consequence(cons_value) {
	filters_consequence.push(cons_value);
	doSearch('clear');
	var listJson = [];
	var grid = dijit.byId("gridId");
    grid.store.fetch({
        onComplete: function(items, request){
        	var max = items.length;
			var n = 0;
			var max = grid.rowCount;
			while (n < max) {
		    	var selectedItem = items[n];
		    	if (selectedItem) {
			        var value = grid.store.getValue(selectedItem, 'consequence');
					var list_values = value.split(', ');
					dojo.forEach(list_values, function(this_value) {
						if (this_value == cons_value) {
					        var json = create_json_from_item(grid, selectedItem);
					        listJson.push(json);
				    	}
			 	 	});
		       }
			   n++;
			}
		    var data_tmp = { "items": listJson };
		    var newStore = new dojo.data.ItemFileReadStore({ data: data_tmp });
			grid.setStore( newStore );
		    drawChart_consequence();
		    drawChart_polyphen();
		    drawChart_sift();
			grid.filter({varId: "*"});
			doSearch_he('clear');
			refresh();
			printNbResults_dejavu();
			
			var listCons = ['filter_intergenic', 'filter_nocoding', 'filter_pseudogene', 'filter_ncRNA', 'filter_mature', 'filter_utr', 'filter_splicing', 'filter_essential_splicing', 'filter_coding', 'filter_frameshift', 'filter_stop', 'filter_phase', 'filter_silent', 'filter_non_frameshift'];
			dojo.forEach(listCons, function(filter_name) {
				dijit.byId(filter_name).set('checked', false);
				dijit.byId(filter_name+'_2').set('checked', true);
				dijit.byId(filter_name).set('disabled', true);
				dijit.byId(filter_name+'_2').set('disabled', false);
			});
			
			dojo.forEach(filters_consequence, function(filter_name) {
				if      (filter_name == "Stop-gained") 			{ filter_name = 'filter_stop'; }
				else if (filter_name == "(Start/Stop)-lost") 	{ filter_name = 'filter_phase'; }
				else if (filter_name == "Missense") 			{ filter_name = 'filter_coding'; }
				else if (filter_name == "Synonymous") 			{ filter_name = 'filter_silent'; }
				else if (filter_name == "Splice Region") 		{ filter_name = 'filter_splicing'; }
				else if (filter_name == "Splice Acc/Don") 		{ filter_name = 'filter_essential_splicing'; }
				else if (filter_name == "Frameshift") 			{ filter_name = 'filter_frameshift'; }
				else if (filter_name == "No-frameshift") 		{ filter_name = 'filter_non_frameshift'; }
				else if (filter_name == "Utr") 					{ filter_name = 'filter_utr'; }
				else if (filter_name == "Pseudogene") 			{ filter_name = 'filter_pseudogene'; }
				else if (filter_name == "ncRNA") 				{ filter_name = 'filter_ncRNA'; }
	            else if (filter_name == "Mature miRNA") 		{ filter_name = 'filter_mature'; }
				else if (filter_name == "Intronic") 			{ filter_name = 'filter_nocoding'; }
				else if (filter_name == "Intergenic") 			{ filter_name = 'filter_intergenic'; }
				dijit.byId(filter_name).set('checked', true);
			});
			doSearch_he('run');
			
			// filter_consequence_heComp(cons_value);
        }
    });
}

function filter_consequence_heComp(cons_value) {
	filters_consequence.push(cons_value);
	doSearch('clear');
	var listJson_heComp = [];
	var gridHeComp = dijit.byId("gridHeComp");
    gridHeComp.store.fetch({
        onComplete: function(items, request){
        	var max = items.length;
			var n = 0;
			var max = gridHeComp.rowCount;
			while (n < max) {
		    	var selectedItem = items[n];
		    	if (selectedItem) {
			        var value1 = gridHeComp.store.getValue(selectedItem, 'consequence_1');
					var list_values1 = value1.split(', ');
			        var value2 = gridHeComp.store.getValue(selectedItem, 'consequence_2');
			        var list_values2 = value2.split(', ');
					var done = 0;
					dojo.forEach(list_values1, function(this_value1) {
						if (this_value1 == cons_value && done == 0) {
					        var json = create_json_hecomp_from_item(gridHeComp, selectedItem);
					        listJson_heComp.push(json);
					        done = 1;
				    	}
			 	 	});
					dojo.forEach(list_values2, function(this_value2) {
						if (this_value2 == cons_value && done == 0) {
					        var json = create_json_hecomp_from_item(gridHeComp, selectedItem);
					        listJson_heComp.push(json);
					        done = 1;
				    	}
			 	 	});
		       }
			   n++;
			}
		    var data_tmp = { "items": listJson_heComp };
		    var newStore = new dojo.data.ItemFileWriteStore({ data: data_tmp });
			gridHeComp.setStore( newStore );
			document.getElementById("aff_he").innerHTML = originLaunchInput_dejavu;
			dijit.byId("panel_he").setTitle("<b>" + panelHeComp_defaultTitle + "</b>");
			valId_forHe = false;
			// gridHeComp.beginUpdate();
			// gridHeComp.filter({id: "*"});
			printNbResults_dejavu_he();
			// gridHeComp.endUpdate();
			// var listCons = ['filter_intergenic', 'filter_nocoding', 'filter_pseudogene', 'filter_ncRNA', 'filter_mature', 'filter_utr', 'filter_splicing', 'filter_essential_splicing', 'filter_coding', 'filter_frameshift', 'filter_stop', 'filter_phase', 'filter_silent', 'filter_non_frameshift'];
			// dojo.forEach(listCons, function(filter_name) {
				// dijit.byId(filter_name).set('checked', true);
				// dijit.byId(filter_name+'_2').set('checked', true);
				// dijit.byId(filter_name).set('disabled', true);
				// dijit.byId(filter_name+'_2').set('disabled', true);
			// });
			// dijit.byId('b_launch_clear_he').set('disabled', true);
			// dijit.byId('b_launch_run_he').set('disabled', true);
			refresh();
        }
    });
}

var filters_polyphen = [];
function filter_polyphen(polyphen_value) {
	filters_polyphen.push(polyphen_value);
	doSearch('clear');
	var listJson = [];
	var grid = dijit.byId("gridId");
    grid.store.fetch({
        onComplete: function(items, request){
        	var max = items.length;
			var n = 0;
			var max = grid.rowCount;
			while (n < max) {
		    	var selectedItem = items[n];
		    	if (selectedItem) {
			        var polyphen = grid.store.getValue(selectedItem, 'polyphen');
			        var lTmp = polyphen.split('+');
			        if (lTmp[0] == polyphen_value) {
				        var json = create_json_from_item(grid, selectedItem);
				        listJson.push(json);
			       }
		       }
			   n++;
			}
		    var data_tmp = { "items": listJson };
		    var newStore = new dojo.data.ItemFileReadStore({ data: data_tmp });
			grid.setStore( newStore );
		    drawChart_consequence();
		    drawChart_polyphen();
		    drawChart_sift();
			grid.filter({varId: "*"});
			doSearch_he('clear');
			refresh();
			printNbResults_dejavu();
        }
    });
}

var filters_sift = [];
function filter_sift(sift_value) {
	filters_sift.push(sift_value);
	doSearch('clear');
	var listJson = [];
	var grid = dijit.byId("gridId");
    grid.store.fetch({
        onComplete: function(items, request){
        	var max = items.length;
			var n = 0;
			var max = grid.rowCount;
			while (n < max) {
		    	var selectedItem = items[n];
		    	if (selectedItem) {
			        var sift = grid.store.getValue(selectedItem, 'sift');
			        var lTmp = sift.split('+');
			        if (lTmp[0] == sift_value) {
				        var json = create_json_from_item(grid, selectedItem);
				        listJson.push(json);
			       }
		       }
			   n++;
			}
		    var data_tmp = { "items": listJson };
		    var newStore = new dojo.data.ItemFileReadStore({ data: data_tmp });
			grid.setStore( newStore );
		    drawChart_consequence();
		    drawChart_polyphen();
		    drawChart_sift();
			grid.filter({varId: "*"});
			doSearch_he('clear');
			refresh();
			printNbResults_dejavu();
        }
    });
}

function create_json_hecomp_from_item(gridHeComp, item) {
    var json = {
    	"id":            gridHeComp.store.getValue(item, 'id'),
    	"implicated":    gridHeComp.store.getValue(item, 'implicated'),
    	"var_1":         gridHeComp.store.getValue(item, 'var_1'),
    	"consequence_1": gridHeComp.store.getValue(item, 'consequence_1'),
    	"var_2":         gridHeComp.store.getValue(item, 'var_2'),
    	"consequence_2": gridHeComp.store.getValue(item, 'consequence_2'),
    	"nbProj":        gridHeComp.store.getValue(item, 'nbProj'),
    	"nbPat":         gridHeComp.store.getValue(item, 'nbPat'),
    }
	return json;
}

function create_json_from_item(grid, item) {
    var json = {
    	"implicated":  grid.store.getValue(item, 'implicated'),
    	"varId":       grid.store.getValue(item, 'varId'),
    	"type":        grid.store.getValue(item, 'type'),
    	"chr":         grid.store.getValue(item, 'chr'),
    	"pos":         grid.store.getValue(item, 'pos'),
    	"gene":        grid.store.getValue(item, 'gene'),
    	"transcript":  grid.store.getValue(item, 'transcript'),
    	"protein":     grid.store.getValue(item, 'protein'),
    	"nbProj":      grid.store.getValue(item, 'nbProj'),
    	"nbPat":       grid.store.getValue(item, 'nbPat'),
    	"nbHo":        grid.store.getValue(item, 'nbHo'),
    	"nbHe":        grid.store.getValue(item, 'nbHe'),
    	"consequence": grid.store.getValue(item, 'consequence'),
    	"polyphen":    grid.store.getValue(item, 'polyphen'),
    	"sift":        grid.store.getValue(item, 'sift'),
    	"freqDbSnp":   grid.store.getValue(item, 'freqDbSnp'),
    	"alamut":      grid.store.getValue(item, 'alamut'),
    }
	return json;
}

function search_sift(value) {
	filter_sift(value);
	var toWrite;
	if (value.match(/0/))      { toWrite = "sift None"; }
	else if (value.match(/1/)) { toWrite = "sift Benign"; }
	else if (value.match(/2/)) { toWrite = "sift Deleterious"; }
	else if (value.match(/4/)) { toWrite = "sift in error"; }
	dijit.byId("panel_graphics").setTitle("<b>" + panelFilters_defaultTitle + "     <FONT COLOR=RED>-> only " + toWrite + " !</FONT></b>");
	dijit.byId("panel_results").setTitle("<b>" + panelResults_defaultTitle + "     <FONT COLOR=RED>-> only " + toWrite + " !</FONT></b>");
}

function search_polyphen(value) {
	filter_polyphen(value);
	var toWrite;
	if (value.match(/0/))      { toWrite = "polyphen None"; }
	else if (value.match(/1/)) { toWrite = "polyphen Benign"; }
	else if (value.match(/2/)) { toWrite = "polyphen Possibly Damaging"; }
	else if (value.match(/3/)) { toWrite = "polyphen Probably Damaging"; }
	else if (value.match(/4/)) { toWrite = "polyphen in error"; }
	dijit.byId("panel_graphics").setTitle("<b>" + panelFilters_defaultTitle + "     <FONT COLOR=RED>-> only " + toWrite + " !</FONT></b>");
	dijit.byId("panel_results").setTitle("<b>" + panelResults_defaultTitle + "     <FONT COLOR=RED>-> only " + toWrite + " !</FONT></b>");
}

function search_consequence(value) {
	var can_do = 1;
	dojo.forEach(filters_consequence, function(filter_name) {
		if (filter_name == value) { can_do = 0; }
	});
	if (can_do == 1) {
		filter_consequence(value);
		var toWrite = value;
		if (value.match(/stop/)) { toWrite = "stop-gained"; }
		else if (value.match(/phase/)) { toWrite = "(start/stop)-lost"; }
		else if (value.match(/coding/)) { toWrite = "non synonymous"; }
		else if (value == "splicing") { toWrite = "splice_site"; }
		else if (value.match(/essential_splicing,coding/)) { toWrite = "essential_splicing, non synonymous"; }
		else if (value.match(/essential_splicing,stop/)) { toWrite = "essential_splicing, stop-gained"; }
		dijit.byId("panel_graphics").setTitle("<b>" + panelFilters_defaultTitle + "     <FONT COLOR=RED>-> only " + toWrite + " !</FONT></b>");
		dijit.byId("panel_results").setTitle("<b>" + panelResults_defaultTitle + "     <FONT COLOR=RED>-> only " + toWrite + " !</FONT></b>");
	}
}
google.load('visualization', '1.0', {'packages':['corechart']});
google.load('visualization', '1', {packages:['table']});

function drawTableResults(value) {
    	var data = new google.visualization.DataTable();
    	data.addColumn('string', 'Variant(s) Found');
    	data.addColumn('number', value);
    	data.addRows([]);
    	var options = {
            	width:200,
            	height:100,
	};
    	var table = new google.visualization.Table(document.getElementById('table_nbResults'));
    	table.draw(data, options);
}

var lastInputText;
function drawChart_sift(inputText) {
	if (inputText) { lastInputText = inputText; }
	else { inputText = lastInputText; }
	var grid = dijit.byId("gridId");
	grid.filter({sift: new RegExp("4")});
	var nbErrors = grid.rowCount;
	grid.filter({sift: new RegExp("0")});
	var nbNone = grid.rowCount;
	grid.filter({sift: new RegExp("1")});
	var nbBenign = grid.rowCount;
	grid.filter({sift: new RegExp("2")});
	var nbPossibly = grid.rowCount;
	var data = new google.visualization.DataTable();
	data.addColumn('string', 'Topping');
	data.addColumn('number', 'Slices');
	data.addRows([
 	  ['Annotation Error', nbErrors],
 	  ['None', nbNone],
 	  ['Benign', nbBenign],
	  ['Deleterious', nbPossibly],
	]);
	var options = {
		title:'Sift status distribution in ' + inputText,
            	width:380,
            	height:330,
		is3D:true,
		slices:{
			0: {color: 'grey'},
			1: {color: '#77B5FE'},
			2: {color: '#16B84E'},
			3: {color: 'red'},
		}
	};
    var chart_pol = new google.visualization.PieChart(document.getElementById('chart_sift_div'));
	google.visualization.events.addListener(chart_pol, 'select', function() {
		var selectedItem = chart_pol.getSelection()[0];
		if (selectedItem) {
			if (selectedItem.row == "0")      { search_sift("4"); }
			else if (selectedItem.row == "1") { search_sift("0"); }
			else if (selectedItem.row == "2") { search_sift("1"); }
			else if (selectedItem.row == "3") { search_sift("2"); }
		}
	});
	chart_pol.draw(data, options);
}

function drawChart_polyphen(inputText) {
	if (inputText) { lastInputText = inputText; }
	else { inputText = lastInputText; }
	var grid = dijit.byId("gridId");
	grid.filter({polyphen: new RegExp("4")});
	var nbErrors = grid.rowCount;
	grid.filter({polyphen: new RegExp("0")});
	var nbNone = grid.rowCount;
	grid.filter({polyphen: new RegExp("1")});
	var nbBenign = grid.rowCount;
	grid.filter({polyphen: new RegExp("2")});
	var nbPossibly = grid.rowCount;
	grid.filter({polyphen: new RegExp("3")});
	var nbProbably = grid.rowCount;
	var data = new google.visualization.DataTable();
	data.addColumn('string', 'Topping');
	data.addColumn('number', 'Slices');
	data.addRows([
 	  ['Annotation Error', nbErrors],
 	  ['None', nbNone],
 	  ['Benign', nbBenign],
	  ['Possibly Damaging', nbPossibly],
	  ['Probably Damaging', nbProbably]
	]);
	var options = {
		title:'Polyphen status distribution in ' + inputText,
            	width:380,
            	height:330,
		is3D:true,
		slices:{
			0: {color: 'grey'},
			1: {color: '#77B5FE'},
			2: {color: '#16B84E'},
			3: {color: '#FCD21C'},
			4: {color: 'red'}
		}
	};
    var chart_pol = new google.visualization.PieChart(document.getElementById('chart_polyphen_div'));
	google.visualization.events.addListener(chart_pol, 'select', function() {
		var selectedItem = chart_pol.getSelection()[0];
		if (selectedItem) {
			if (selectedItem.row == "0") {search_polyphen("4");}
			else if (selectedItem.row == "1") {search_polyphen("0");}
			else if (selectedItem.row == "2") {search_polyphen("1");}
			else if (selectedItem.row == "3") {search_polyphen("2");}
			else if (selectedItem.row == "4") {search_polyphen("3");}
		}
	});
	chart_pol.draw(data, options);
}

function drawChart_consequence(inputText) {
	if (inputText) { lastInputText = inputText; }
	else { inputText = lastInputText; }
	var grid = dijit.byId("gridId");
	grid.filter({consequence: "Stop-gained"});
	var nbstop = grid.rowCount;
	grid.filter({consequence: new RegExp("lost")});
	var nbphase = grid.rowCount;
	grid.filter({consequence: new RegExp("Missense")});
	var nbnonsynonymous = grid.rowCount;
	grid.filter({consequence: new RegExp("Synonymous")});
	var nbsilent = grid.rowCount;
	grid.filter({consequence: new RegExp("Splice Region")});
	var nbsplicing = grid.rowCount;
	grid.filter({consequence: "Splice Acc/Don"});
	var nbessential_splicing = grid.rowCount;
	grid.filter({consequence: "Frameshift"});
	var nbframeshift = grid.rowCount;
	grid.filter({consequence: "No-frameshift"});
	var nbnonframeshift = grid.rowCount;
	grid.filter({consequence: new RegExp("Utr")});
	var nbutr = grid.rowCount;
	grid.filter({consequence: "Pseudogene"});
	var nbpseudogene = grid.rowCount;
	grid.filter({consequence: "ncRNA"});
	var nbncrna = grid.rowCount;
    grid.filter({consequence: "mature-miRNA"});
    var nbmature = grid.rowCount;
	grid.filter({consequence: "Intronic"});
	var nbintronic = grid.rowCount;
	grid.filter({consequence: "Intergenic"});
	var nbintergenic = grid.rowCount;
	grid.filter({consequence: "Error annotation..."});
	var nberrors = grid.rowCount;
	var txtstop = 'Stop-gained';
	var txtphase = '(Start/Stop)-lost';
	var txtnonsynonymous = 'Missense';
	var txtsilent = 'Synonymous';
	var txtsplicing = 'Splice Region';
	var txtessential_splicing = 'Splice Acc/Don';
	var txtframeshift = 'Frameshift';
	var txtnonframeshift = 'No-frameshift';
	var txtutr = 'Utr';
	var txtpseudogene = 'Pseudogene';
	var txtncrna = 'ncRNA';
	var txtmature = 'Mature miRNA';
	var txtintronic = 'Intronic';
	var txtintergenic = 'Intergenic';
	
	var nbfilters = filters_consequence.length;
	if (nbfilters > 0) {
		var nbtotal = nbstop + nbphase + nbnonsynonymous + nbsilent + nbsplicing + nbessential_splicing + nbframeshift + nbnonframeshift + nbutr + nbpseudogene + nbncrna + nbmature + nbintronic + nbintergenic;
		var i = 0;
		dojo.forEach(filters_consequence, function(filter_name) {
			i++;
			var txtFistFilters = filters_consequence.join(' + ');
			
			if      (filter_name == "Stop-gained") 	{ nbstop = getNbAndTextFilterGraph(nbtotal, nbstop, i, 1); }
			else if (nbfilters == 1)				{ txtstop = txtFistFilters + ' + ' + txtstop; }
			else if (nbfilters == 2)				{ txtstop = txtFistFilters; }
			
			if (filter_name == "(Start/Stop)-lost") { nbphase = getNbAndTextFilterGraph(nbtotal, nbphase, i, 1); }
			else if (nbfilters == 1)				{ txtphase = txtFistFilters + ' + ' + txtphase; }
			else if (nbfilters == 2)				{ txtphase = txtFistFilters; }
			
			if (filter_name == "Missense") 			{ nbnonsynonymous = getNbAndTextFilterGraph(nbtotal, nbnonsynonymous, i, 1); }
			else if (nbfilters == 1)				{ txtnonsynonymous = txtFistFilters + ' + ' + txtnonsynonymous; }
			else if (nbfilters == 2)				{ txtnonsynonymous = txtFistFilters; }
			
			if (filter_name == "Synonymous") 		{ nbsilent = getNbAndTextFilterGraph(nbtotal, nbsilent, i, 1); }
			else if (nbfilters == 1)				{ txtsilent = txtFistFilters + ' + ' + txtsilent; }
			else if (nbfilters == 2)				{ txtsilent = txtFistFilters; }
			
			if (filter_name == "Splice Region") 	{ nbsplicing = getNbAndTextFilterGraph(nbtotal, nbsplicing, i, 1); }
			else if (nbfilters == 1)				{ txtsplicing = txtFistFilters + ' + ' + txtsplicing; }
			else if (nbfilters == 2)				{ txtsplicing = txtFistFilters; }
			
			if (filter_name == "Splice Acc/Don") 	{ nbessential_splicing = getNbAndTextFilterGraph(nbtotal, nbessential_splicing, i, 1); }
			else if (nbfilters == 1)				{ txtessential_splicing = txtFistFilters + ' + ' + txtessential_splicing; }
			else if (nbfilters == 2)				{ txtessential_splicing = txtFistFilters; }
			
			if (filter_name == "Frameshift") 		{ nbframeshift = getNbAndTextFilterGraph(nbtotal, nbframeshift, i, 1); }
			else if (nbfilters == 1)				{ txtframeshift = txtFistFilters + ' + ' + txtframeshift; }
			else if (nbfilters == 2)				{ txtframeshift = txtFistFilters; }
			
			if (filter_name == "No-frameshift") 	{ nbnonframeshift = getNbAndTextFilterGraph(nbtotal, nbnonframeshift, i, 1); }
			else if (nbfilters == 1)				{ txtnonframeshift = txtFistFilters + ' + ' + txtnonframeshift; }
			else if (nbfilters == 2)				{ txtnonframeshift = txtFistFilters; }
			
			if (filter_name == "Utr") 				{ nbutr = getNbAndTextFilterGraph(nbtotal, nbutr, i, 1); }
			else if (nbfilters == 1)				{ txtutr = txtFistFilters + ' + ' + txtutr; }
			else if (nbfilters == 2)				{ txtutr = txtFistFilters; }
			
			if (filter_name == "Pseudogene") 		{ nbpseudogene = getNbAndTextFilterGraph(nbtotal, nbpseudogene, i, 1); }
			else if (nbfilters == 1)				{ txtpseudogene = txtFistFilters + ' + ' + txtpseudogene; }
			else if (nbfilters == 2)				{ txtpseudogene = txtFistFilters; }
			
			if (filter_name == "ncRNA") 			{ nbncrna = getNbAndTextFilterGraph(nbtotal, nbncrna, i, 1); }
			else if (nbfilters == 1)				{ txtncrna = txtFistFilters + ' + ' + txtncrna; }
			else if (nbfilters == 2)				{ txtncrna = txtFistFilters; }
			
	        if (filter_name == "Mature miRNA") 		{ nbmature = getNbAndTextFilterGraph(nbtotal, nbmature, i, 1); }
			else if (nbfilters == 1)				{ txtmature = txtFistFilters + ' + ' + txtmature; }
			else if (nbfilters == 2)				{ txtmature = txtFistFilters; }
			
			if (filter_name == "Intronic") 			{ nbintronic = getNbAndTextFilterGraph(nbtotal, nbintronic, i, 1); }
			else if (nbfilters == 1)				{ txtintronic = txtFistFilters + ' + ' + txtintronic; }
			else if (nbfilters == 2)				{ txtintronic = txtFistFilters; }
			
			if (filter_name == "Intergenic") 		{ nbintergenic = getNbAndTextFilterGraph(nbtotal, nbintergenic, i, 1); }
			else if (nbfilters == 1)				{ txtintergenic = txtFistFilters + ' + ' + txtintergenic; }
			else if (nbfilters == 2)				{ txtintergenic = txtFistFilters; }
		});
	}
	
	var data = new google.visualization.DataTable();
	data.addColumn('string', 'Topping');
	data.addColumn('number', 'Slices');
	data.addRows([
 	  ['Annotation Error', nberrors],
 	  [txtstop, nbstop],
 	  [txtphase, nbphase],
 	  [txtnonsynonymous, nbnonsynonymous],
 	  [txtsilent, nbsilent],
 	  [txtsplicing, nbsplicing],
 	  [txtessential_splicing, nbessential_splicing],
 	  [txtframeshift, nbframeshift],
 	  [txtnonframeshift, nbnonframeshift],
 	  [txtutr, nbutr],
 	  [txtpseudogene, nbpseudogene],
 	  [txtncrna, nbncrna],
      [txtmature, nbmature],
 	  [txtintronic, nbintronic],
 	  [txtintergenic, nbintergenic]
	]);
	var options = {
		title:'Consequence distribution in ' + inputText,
            	width:380,
            	height:330,
		is3D:true
	};
    var chart_cons = new google.visualization.PieChart(document.getElementById('chart_consequence_div'));
	google.visualization.events.addListener(chart_cons, 'select', function() {
		var selectedItem = chart_cons.getSelection()[0];
		if (selectedItem) {
			if (selectedItem.row == "0") {search_consequence("Error annotation...");}
			else if (selectedItem.row == "1") {search_consequence("Stop-gained");}
			else if (selectedItem.row == "2") {search_consequence("(Start/Stop)-lost");}
			else if (selectedItem.row == "3") {search_consequence("Missense");}
			else if (selectedItem.row == "4") {search_consequence("Synonymous");}
			else if (selectedItem.row == "5") {search_consequence("Splice Region");}
			else if (selectedItem.row == "6") {search_consequence("Splice Acc/Don");}
			else if (selectedItem.row == "7") {search_consequence("Frameshift");}
			else if (selectedItem.row == "8") {search_consequence("No-frameshift");}
			else if (selectedItem.row == "9") {search_consequence("Utr");}
			else if (selectedItem.row == "10") {search_consequence("Pseudogene");}
			else if (selectedItem.row == "11") {search_consequence("ncRNA");}
            else if (selectedItem.row == "12") {search_consequence("Mature miRNA");}
			else if (selectedItem.row == "13") {search_consequence("Intronic");}
			else if (selectedItem.row == "14") {search_consequence("Intergenic");}
		}
	});
    chart_cons.draw(data, options);
}

function onlyUnique(value, index, self) { 
    return self.indexOf(value) === index;
}

function getNbAndTextFilterGraph(nbtotal, nbfilter, step, is_filter) {
	if (is_filter == 1) {
		if (step == 1) {
			nbfilter = nbtotal - 2 * (nbtotal - nbfilter);
		}
		else if (step == 2) {
			nbfilter = nbtotal - (nbtotal - nbfilter);
		}
	}
	return nbfilter;
}

function openPopup() { dijit.byId("popupDetails").show(); }

function openPopupLog() { dijit.byId("popupLog").show(); }

function add_arg_filter(list) {
	var lRes = [];
	dojo.forEach(list, function(value) {
		var lTmp = value.split('|');
		dojo.forEach(lTmp, function(value2) {
			lRes.push(value2);
		});
	});
	var arg = lRes.join(',');
	return arg;
}



function downloadFile(url) {
    var downloadPdfIframeName = "downloadPdfIframe"; 
    var iframe = dojo.io.iframe.create(downloadPdfIframeName);
    dojo.io.iframe.setSrc(iframe, url, true);
} 

var deferred_xls;
function downloadFileSession(url){
    dijit.byId("waiting").show();
	require(["dojo/request/xhr"], function(xhr){
        xhr(url, {}).then(function(result) {
            var lTmp = [];
            lTmp = result.split('@@@');
            var sessionId = lTmp.pop();
            var url2 = url + "&xls_load_session=" + String(sessionId);
		    dijit.byId("waiting").hide();
		    dijit.byId("waiting_XLS_OK").show();
		    downloadFile(url2);
        }, function(err) {
        	alert('ERROR XLS');
        }, function(evt) {
        
        });
    });
}

function export_xls_results() {
	var user = dojo.cookie("username");
	var pass = dojo.cookie("passwd");
	var build = 'HG19';
	var urlLaunch_dejavu = last_url + "&xls=1";
	if (filter_var_ids) {
		urlLaunch_dejavu += "&f_ids=" + filter_var_ids;
	}
	else {
		if (filters_sift.length > 0) {
			var arg_f = add_arg_filter(filters_sift);
			urlLaunch_dejavu += "&f_sift=" + arg_f;
		}
		if (filters_polyphen.length > 0) {
			var arg_f = add_arg_filter(filters_polyphen);
			urlLaunch_dejavu += "&f_polyphen=" + arg_f;
		}
		if (filters_consequence.length > 0) {
			var arg_f = add_arg_filter(filters_consequence);
			urlLaunch_dejavu += "&f_cons=" + arg_f;
		}
	}
    alert ("Click OK and wait your file patiently please. Thanks.");
	downloadFileSession(urlLaunch_dejavu);
}

function export_xls_he_comp() {
	var user = dojo.cookie("username");
	var pass = dojo.cookie("passwd");
	var build = 'HG19';
	var urlLaunch_dejavu = last_url_he_comp + "&xls=1";
	if (filters_consequence.length > 0) {
		var arg_f = add_arg_filter(filters_consequence);
		var arg_cons = arg_f.replace(/ /g, '_');
		urlLaunch_dejavu += "&f_cons=" + arg_cons;
	}
	else {
		if (filter_var_id_hecomp) {
			urlLaunch_dejavu += "&f_id=" + filter_var_id_hecomp;
			if (filters_cons_hecomp_2) {
				var arg_cons = filters_cons_hecomp_2.replace(/ /g, '_');
				urlLaunch_dejavu += "&f_cons2=" + arg_cons;
			}
		}
		else {
			if (filters_cons_hecomp_1) {
				var arg_cons = filters_cons_hecomp_1.replace(/ /g, '_');
				urlLaunch_dejavu += "&f_cons1=" + arg_cons;
			}
			if (filters_cons_hecomp_2) {
				var arg_cons = filters_cons_hecomp_2.replace(/ /g, '_');
				urlLaunch_dejavu += "&f_cons2=" + arg_cons;
			}
		}
	}
    alert ("Click OK and wait your file patiently please. Thanks.");
	downloadFileSession(urlLaunch_dejavu);
}

function enabledXlsExport() {
    dijit.byId("b_export_xls_results").set("disabled", false);
    dijit.byId("b_export_xls_he_comp").set("disabled", false);
}

function disabledXlsExport() {
    dijit.byId("b_export_xls_results").set("disabled", true);
    dijit.byId("b_export_xls_he_comp").set("disabled", true);
}

function project_stats(){
	var username = dojo.cookie("username");
	var passwd = dojo.cookie("passwd");
	var url_stat_polyweb1 =url_stat_polyweb+"?user="+username+"&pwd="+passwd;
	if (force_annot_version) {
		url_stat_polyweb1 += "&force_annot_version="+force_annot_version;
	}
	else if (param("origin_project")) {
		url_stat_polyweb1 += "&origin_project="+param("origin_project");
	}
	dojo.xhrGet({
		url:url_stat_polyweb1, handleAs:"json",
		load: function(data){
			document.getElementById("span_projects").innerHTML = " "+data.user.projects.all+" / "+data.global.projects.all;
			document.getElementById("span_samples").innerHTML = " "+data.user.patients.all+" / "+data.global.patients.all;
			document.getElementById("span_exomes").innerHTML = " "+data.user.patients.exome+" / "+data.global.patients.exome;
			document.getElementById("span_genomes").innerHTML = " "+data.user.patients.genome+" / "+data.global.patients.genome;
			document.getElementById("span_ciliomes").innerHTML = " "+data.user.patients.ciliome+" / "+data.global.patients.ciliome;
			document.getElementById("span_diagnostics").innerHTML = " "+data.user.patients.target+" / "+data.global.patients.target;
			document.getElementById("aff_annot_version").innerHTML = " "+data.global.projects.annot_version;
			if (data.global.projects.last_annot_version == data.global.projects.annot_version) {
				document.getElementById("aff_annot_version").innerHTML += " <img src='../images/polyicons/bullet_green.png' alt='Warn' style='width:11px;height:11px;'>";
				document.getElementById("span_tooltip_annot_version").textContent = "This is the LAST annotation version.";
			}
			else {
				document.getElementById("aff_annot_version").innerHTML += " <img src='../images/polyicons/error.png' alt='Warn' style='width:11px;height:11px;'>";
				document.getElementById("span_tooltip_annot_version").textContent = "The LAST annotation version available is " + data.global.projects.last_annot_version;
			}
			last_annot_version = data.global.projects.last_annot_version;
			force_annot_version = last_annot_version;
			var menu = new dijit.Menu( { id: "menu_annot_versions" } );
			var list_all_annot = data.global.projects.all_annot_version;
            for ( i = 0; i < list_all_annot.length; i++ ) {
            	var lGenecode = String(list_all_annot[i]['id']).split('.');
            	if (lGenecode[0] == '32') {}
            	else {
					var onclick_method = "change_annot('" + list_all_annot[i]['id'] + "')";
					var text = "<b><u><span style='color:red'>Version " + list_all_annot[i]['id'] + '</span></b></u>   <i>[<b>Genecode</b> ' + lGenecode[0] + ' - ' + list_all_annot[i]['value'] + ']</i>';
					menu.addChild(new dijit.MenuItem ( { label: text, onClick: onclick_method }));
				}
            }
			menu.placeAt("div_dropdownmenu_annotations");
			menu.startup();
		}
	});
	to_do_stats = '';
}

var last_annot_version;
var force_annot_version;
function change_annot(annot_version) {
	force_annot_version = annot_version;
	document.getElementById("aff_annot_version").innerHTML = " "+annot_version;
	if (annot_version == last_annot_version) {
		document.getElementById("aff_annot_version").innerHTML += " <img src='../images/polyicons/bullet_green.png' alt='Warn' style='width:11px;height:11px;'>";
		document.getElementById("span_tooltip_annot_version").textContent = "This is the LAST annotation version.";
	}
	else {
		document.getElementById("aff_annot_version").innerHTML += " <img src='../images/polyicons/error.png' alt='Warn' style='width:11px;height:11px;'>";
		document.getElementById("span_tooltip_annot_version").textContent = "The LAST annotation version available is " + last_annot_version;
	}
}
        
function zoomHgmd(p,hid,vid){
    dijit.byId("popup_hgmd").show();
    dijit.byId("div_popup_hgmd").setHref(url_hgmd_view+"?polyquery_view=1&hid="+hid+"&vid="+vid+"&project="+p);
}