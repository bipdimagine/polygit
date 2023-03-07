/**
 * 
 */

function construct_content_pane_dude()  {
		require(["dojo/ready", "dijit/layout/ContentPane"], function(ready, ContentPane){
  ready(function(){
	
   var t =  new ContentPane({
      style:"font-size:9px;width:100%;overflow-x:auto;",
	title:"PolyDu<sub>P</sub>De<sub>L</sub>",
	id:"panel_polydude",
	jsId:"panel_polydude",
	href:"include_html/diag/editor_container_dude.html",
	preload:1,
    });
	/*t.attr('href', "include_html/diag/editor_container_dude.html");*/
	/*t.connect(t, 'onLoad', load_tab_polydude());*/
	tabPrincipal.addChild(t);
	t.connect(t, 'onLoad',function(){load_tab_polydude()});
  });
});
}	

function view_variation_validation(id, gene_id){
	var argsPost = new Array();
    argsPost["project"] = project_name;
 	argsPost["vid"] = id;
 	if (gene_id) { argsPost["gid"] = gene_id; }
 	//alert(url_view__variations_validations);
	 post_show_popup_dude(argsPost,url_view__variations_validations);
}

function view_variation_validation_for_project(this_project_name, id){
	var argsPost = new Array();
    argsPost["project"] = this_project_name;
 	argsPost["vid"] = id;
	 post_show_popup_dude(argsPost,url_view__variations_validations);
}

function view_dejavu_polycyto(value) {
	dijit.byId('waiting').show();
	var lCol = value.split(';');
	var infos_proj_pat = lCol[-1];
	var build = 'HG19';
	var ltmp = String(value).split(';');
	var args = "login=" + dojo.cookie("username") + "&pwd=" + dojo.cookie("passwd") + "&values=" + ltmp[4] + "&my_phenotype=" + document.getElementById("span_phenotype_name").innerHTML;
	var xhr = new XMLHttpRequest();
	xhr.open('POST', url_proj_pat_details, true);
	xhr.setRequestHeader('Content-type', 'application/x-www-form-urlencoded');
	xhr.onreadystatechange = function() {
	    if (this.readyState === XMLHttpRequest.DONE && this.status === 200) {
	    	var data = JSON.parse(this.responseText);
			document.getElementById("content_dejavu_polycyto").innerHTML = data.html_table;
			dijit.byId('waiting').hide();
			dijit.byId('dialog_dejavu_polycyto').show();
	    }
	}
	xhr.send(args);
}

var gene_name ;
var first_view= 1;

var last_load_polycyto_argsPost = {};
function load_polycyto_args(mode,dm) {
	var argsPost = {
		projectname: project_name,
    	viewer: 'yes',
		iframe : 1,
	};
	last_load_polycyto_argsPost = argsPost;
	return argsPost;
}

var last_load_polydiag_argsPost = {};
var last_load_polyviewer_argsPost = {};

function load_polydude_args() {
	var argsPost = {project:project_name,
					quality:1
					};
	if ( dijit.byId("slider_dude_ac")){
		argsPost.quality = dijit.byId("slider_dude_ac").value
	}
	
	var v_g = 'all';
	var gene_name = document.getElementById('gene_name_dude').value;
	if (gene_name) { argsPost.quality = gene_name; }
	
	return argsPost;			
}

function load_polydiag_args() {
	var v_qual = dijit.byId("ratio_label").get("value");
	if (v_qual) { v_qual = v_qual * -1; }
	else {
		v_qual = dijit.byId("slider_allele_quality").value;
		dijit.byId("ratio_label").set("value", "");
	}
	var v_panel = '';
	var v_never = dijit.byId("check_never").get('value');
	var v_this = dijit.byId("slider_this").value;
	var v_impact = dijit.byId("slider_consequence").value;
	var v_freq = dijit.byId("slider_frequence").value;
	var v_report_mode = 1;
	var argsPost = {
		project:project_name,
		panel: v_panel,
		edit_mode:1,
		never: v_never,
		this: v_this,
		impact: v_impact,
		frequence: v_freq,
		allele_quality: v_qual,
		report_mode: v_report_mode,
		transcripts: 'all',
	};
	last_load_polydiag_argsPost = argsPost;
	return argsPost;
} 

function load_polyviewer_args(mode,dm) {

	var argsPost = {project:project_name,
    		edit_mode: mode,
			only_DM : dm,
    		never: dijit.byId("check_never_pv").get('value'),
    		ach: dijit.byId("slider_gnomad_ac_ho_pv").value,
    		ac: dijit.byId("slider_gnomad_ac_pv").value,
    		dv: dijit.byId("slider_dv_pv").value,
    		dv_ho: dijit.byId("slider_dv_ho_pv").value,
			in_this_run: dijit.byId("slider_polyviewer_this").value,
    		allele_quality: document.getElementById("htoto").value,
    		report_mode:1,
    		user_name: cookie_username,
    		annot: functional_url_filter_data(),
    	
		
		};
	/*filtering argument*/
	if (document.getElementById("keep_pathogenic")) {
		if (document.getElementById("keep_pathogenic").checked == true)  argsPost["keep_pathogenic"] = 1;
	}
	
	if (document.getElementById("strict_denovo").checked == true) argsPost["strict_denovo"] =1;
	if (document.getElementById("denovo").checked == true) argsPost["denovo"] =1;
	if (document.getElementById("recessive").checked == true) argsPost["recessive"] =1;
	if (document.getElementById("xor_mother").checked == true && document.getElementById("xor_father").checked == true ) {
		argsPost["xor"] =1;
	}
	else if  (document.getElementById("xor_mother").checked == true){
		argsPost["xor_mother"] =1;
	}
	else if  (document.getElementById("xor_father").checked == true){
		argsPost["xor_father"] =1;
	}
	if (document.getElementById("both").checked == true) argsPost["both"] =1;
	
	
	/*Panel */ 
	var panel_name;
	if (document.getElementById("span_panel_name")){
		panel_name = document.getElementById("span_panel_name").innerHTML;
		if (String(panel_name).match(/Gene:/g)) { panel_name = 'all'; }
	}
	else if (vectors_hex) {
		document.getElementById("span_panel_name").innerHTML = "all";
	}
	else {
		panel_name ="sysid";
	}
	if (panel_name == 'All Genes') {
		panel_name = '';
	}
	
	/* phenotype */
	
	
	var phenotype_name;
	if (panel_name == '') {
		phenotype_name = document.getElementById("span_phenotype_name").innerHTML;
		if (phenotype_name == 'All Genes') { phenotype_name = ''; }
		if (phenotype_name == 'No Phenotype Associated') { phenotype_name = ''; }
		if (String(phenotype_name).match(/Gene:/g)) { phenotype_name = ''; }
	}
	
	if (gene_name == '') {}
	else {
		argsPost["gene_name"] = gene_name;
	}
	
    if (panel_name && panel_name != 'all') {
        //url+= "&panel=" + panel_name;
    	argsPost["panel"] = panel_name;
    }
    if (phenotype_name) {	
    	argsPost["phenotype"] = phenotype_name;
    }
	last_load_polyviewer_argsPost = argsPost;
	return argsPost;
}



function refresh_dude() {
   if (tsamples.length == 0) { tsamples = storeP; }
	 for (var i=0;i<tsamples.length;i++) {
	   var label = tsamples[i].label;
	   //tab_editor_polydude[label].setArgsPost = argsPost;
		
	   //tab_editor_polydude[label].transcripts = t;
	   tab_editor_polydude[tsamples[i].label].toto = 1;
   }
 	if (tab_editor_polydude){
     tab_selected_patient = used_patient_in_dropdown_menu;
	   tab_editor_polydude[tab_selected_patient].toto = 1;
	   tab_editor_polydude[tab_selected_patient].onShow();
   }
   return;
}

function cgh(){
	if (tab_editor_polydude){
     tab_selected_patient = used_patient_in_dropdown_menu;
	   tab_editor_polydude[tab_selected_patient].toto = 1;
 		tab_editor_polydude[tab_selected_patient].cgh = 1;
     	tab_editor_polydude[tab_selected_patient].onShow();
   }
}


function load_polyviewer_export_xls(mode) {
	dijit.byId('waiting').show();
	var url;
    var argsPost = load_polyviewer_args(mode);
	var tab_selected_patient = used_patient_in_dropdown_menu;
	argsPost["patients"] = tab_selected_patient;
	argsPost["export_xls"] = 1;
	var args = '';
	var i = 0;
	for (var key in argsPost) {
		var value = argsPost[key];
		if (value) {
			i++;
			if (i == 1) { args += "?" + key + "=" + value; }
			else { args += "&" + key + "=" + value; }
		}
	}
	var url_xls_2 = url_polyviewer + args;
	console.log(url_xls_2);
    var downloadPdfIframeName = "downloadPdfIframe"; 
    var iframe = dojo.io.iframe.create(downloadPdfIframeName);
    dojo.io.iframe.setSrc(iframe, url_xls_2, true);
    return;
}

function load_polydiag() {
	var argsPost = load_polydiag_args();
	if (tsamples.length == 0) { tsamples = storeP; }
	for (var i=0;i<tsamples.length;i++) {
		var label = tsamples[i].label;
		tab_editor_polydiag[label].setArgsPost = argsPost;
		tab_editor_polydiag[tsamples[i].label].toto = 1;
	}
	if (tab_editor_polydiag){
		tab_editor_polydiag[used_patient_in_dropdown_menu].onShow();
	}
	return;
}

function load_polyviewer(mode, only_dude,dm){
	var url;
    var argsPost;
	
    if (only_dude) {
    	argsPost = last_load_polyviewer_argsPost;
    	argsPost["level_dude"] = only_dude;
    }
    else {
		argsPost = load_polyviewer_args(mode,dm);
    }
	
   /*send argument to all tab*/ 
   var tab_selected_patient;

   t = "all";
   if (tsamples.length == 0) { tsamples = storeP; }
   for (var i=0;i<tsamples.length;i++) {
	   var label = tsamples[i].label;
	   tab_editor_polyviewer[label].setArgsPost = argsPost;
	   tab_editor_polyviewer[label].transcripts = t;
	   tab_editor_polyviewer[tsamples[i].label].toto = 1;
   }
	if (tab_editor_polyviewer){
		tab_editor_polyviewer[used_patient_in_dropdown_menu].onShow();
   }
   return;
}

var var_load_polyviewer  = 0;
function load_tab_polyviewer(){
	 for (var i=0;i<tsamples.length;i++) {
	    	 container_editor_polyviewer.addChild(tab_editor_polyviewer[tsamples[i].label]);
			tab_editor_polyviewer[tsamples[i].label].add = 1;
			
	    }
	if (var_load_polyviewer == 1) refresh_tab_polyviewer() ;
}

var var_load_polycyto  = 0;
function load_tab_polycyto(){
	if (var_load_polycyto == 1) {
		refresh_tab_polycyto();
		return;
	}
	for (var i=0;i<tsamples.length;i++) {
    	container_editor_polycyto.addChild(tab_editor_polycyto[tsamples[i].label]);
		tab_editor_polycyto[tsamples[i].label].add = 1;
	}
	var_load_polycyto = 1;
}
function refresh_tab_polycyto(){
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
			if (tab_editor_polycyto[tsamples[i].label].add == 0){
				container_editor_polycyto.addChild(tab_editor_polycyto[tsamples[i].label]);
				tab_editor_polycyto[tsamples[i].label].add = 1;
			}
			if (first < 1){
			tab_editor_polycyto[name].onShow();
			first ++;
		}
		}
		else {
			if (tab_editor_polycyto[tsamples[i].label].add == 1) {
				container_editor_polycyto.removeChild(tab_editor_polycyto[tsamples[i].label]);
				tab_editor_polycyto[tsamples[i].label].add = 0;
				tab_editor_polycyto[name].set('selected',false);
			}
		}
	}
}

var var_load_polydude  = 0;
function load_tab_polydude(){
	if (var_load_polydude == 1) {
		refresh_tab_polydude();
		return;
	}
	var  tsamples = storeP._arrayOfAllItems;
	
	for (var i=0;i<tsamples.length;i++) {
    	container_editor_polydude.addChild(tab_editor_polydude[tsamples[i].label]);
		tab_editor_polydude[tsamples[i].label].add = 1;
	}

	var_load_polydude = 1;
}

function refresh_tab_polydude(){
	if ( !(dijit.byId("gridPatients"))) {
		return 1;
	}
	for (var i=0;i<tsamples.length;i++) {
		tab_editor_polydude[tsamples[i].label].toto = 0;
	}
	var tsamples2 = dijit.byId("gridPatients").selection.getSelected();
	var selected = new Array();
	 for (var i=0;i<tsamples2.length;i++) {
		selected[tsamples2[i].label] = 1;
	}
	var first = 0;
	for (var i=0;i<tsamples.length;i++) {
		var name = tsamples[i].label;
		if (selected[name] == 1){
			if (tab_editor_polydude[tsamples[i].label].add == 0){
				container_editor_polydude.addChild(tab_editor_polydude[tsamples[i].label]);
				tab_editor_polydude[tsamples[i].label].add = 1;
			}
			if (first < 1){
				tab_editor_polydude[name].onShow();
				first ++;
			}
			else {
				tab_editor_polydude[tsamples[i].label].toto = 1;
			}
		}
		else {
			if (tab_editor_polydude[tsamples[i].label].add == 1) {
				container_editor_polydude.removeChild(tab_editor_polydude[tsamples[i].label]);
				tab_editor_polydude[tsamples[i].label].add = 0;
				tab_editor_polydude[name].set('selected',false);
			}
		}
	}
}

//var tab_editor_polyviewer_launch
function refresh_tab_polyviewer(){
	if ( !(dijit.byId("gridPatients"))) {
		return 1;
	}
	for (var i=0;i<tsamples.length;i++) {
		tab_editor_polyviewer[tsamples[i].label].toto = 0;
	}
	var tsamples2 = dijit.byId("gridPatients").selection.getSelected();
	var selected = new Array();
	 for (var i=0;i<tsamples2.length;i++) {
		selected[tsamples2[i].label] = 1;
	}
	var first = 0;
	for (var i=0;i<tsamples.length;i++) {
		var name = tsamples[i].label;
		if (selected[name] == 1){
			if (tab_editor_polyviewer[tsamples[i].label].add == 0){
				container_editor_polyviewer.addChild(tab_editor_polyviewer[tsamples[i].label]);
				tab_editor_polyviewer[tsamples[i].label].add = 1;
			}
			if (first < 1){
				first ++;
				tab_editor_polyviewer[name].onShow();
			}
			else {
				tab_editor_polyviewer[tsamples[i].label].toto = 1;
			}
		}
		else {
			if (tab_editor_polyviewer[tsamples[i].label].add == 1) {
				container_editor_polyviewer.removeChild(tab_editor_polyviewer[tsamples[i].label]);
				tab_editor_polyviewer[tsamples[i].label].add = 0;
				tab_editor_polyviewer[name].set('selected',false);
			}
		}
	}
}



function validation_negative (pname,value,icon){
 	var url = url_valid_acmg+"?vid=0"+"&sample="+pname+"&project="+project_name+"&user="+cookie_username+"&value=-99";
 	 var stamp = document.getElementById( "stamp_"+pname);
 	 dojo.xhrGet({ 
         // The following URL must match that used to test the server.
         url: url,
         handleAs: "json",     
         timeout: 80000, // Time in milliseconds
         // The LOAD function will be called on a successful response.
         
         load: function(responseObject, ioArgs,evt){
             if (responseObject["status"] == "OK") {
            	 stamp.innerHTML = "negative";
				 stamp.style.visibility = "visible";
                  return;
             }       
             else  {
            	 alert("error");
                  return;
             }   
         
             
         },
         
         // The ERROR function will be called in an error case.
         error: function(response, ioArgs){ // &acirc;&#382;&fnof;
             alert("Error !!!!");
            
             return;
         }
     });
   	return;
 	 
	
}

function validation_acmg_cnv (pname,vid,select,tdid){
	var username = dojo.cookie("username");
	//alert("validation disabled");
	//return;
  	var url = url_valid_acmg+"?vid="+vid+"&sample="+pname+"&project="+project_name+"&user="+cookie_username+"&value="+select.value+"&cnv="+1;
    var div = document.getElementById("div_"+tdid);
    div.innerHTML = "<i class='fa fa-spinner fa-spin' style=font-size></i>"; 
    var stamp = document.getElementById( "stamp_"+pname);
    var td = document.getElementById("td_"+tdid);
    var n = document.getElementById("negative_"+pname);
    td.classList.remove("info");
    td.classList.remove("success");
     td.classList.remove("warning");
    td.classList.remove("danger");
    dojo.xhrGet({ 
        // The following URL must match that used to test the server.
        url: url,
        handleAs: "json",     
        timeout: 80000, // Time in milliseconds
        // The LOAD function will be called on a successful response.
        
        load: function(responseObject, ioArgs,evt){
            if (responseObject["status"] == "OK") {
                div.innerHTML = "<i class='fa fa-check-square f' style='font-size:13px;color:white'> SAVED </i>";
               td.classList.add("success");
               stamp.innerHTML = "validated";
               stamp.style.visibility = "visible";
               n.style.visibility ="hidden";
                 return;
            }       
            else  {
                div.innerHTML = "<i class='fa fa-exclamation-triangle style='font-size:13px;color:white''> ERROR</i>";
                td.classList.add("danger");
                 return;
            }   
        
            
        },
        
        // The ERROR function will be called in an error case.
        error: function(response, ioArgs){ // &acirc;&#382;&fnof;
            alert("Error !!!!");
               div.innerHTML = '<i class="fa fa-exclamation-triangle"> Error</i>';
           
            return;
        }
    });
  	return;


}

function validation_acmg_for_project (this_project_name,pname,vid,select,tdid,ensg){
	project_name = this_project_name;
	validation_acmg (pname,vid,select,tdid,ensg);
}

function validation_acmg (pname,vid,select,tdid,ensg){
	var username = dojo.cookie("username");
	//return;
  	var url = url_valid_acmg+"?vid="+vid+"&sample="+pname+"&project="+project_name+"&user="+cookie_username+"&value="+select.value+"&gene="+ensg;
    var div = document.getElementById("div_"+tdid);
    div.innerHTML = "<i class='fa fa-spinner fa-spin' style=font-size></i>"; 
    var stamp = document.getElementById( "stamp_"+pname);
    var td = document.getElementById("td_"+tdid);
    var n = document.getElementById("negative_"+pname);
    td.classList.remove("info");
    td.classList.remove("success");
     td.classList.remove("warning");
    td.classList.remove("danger");
    dojo.xhrGet({ 
        // The following URL must match that used to test the server.
        url: url,
        handleAs: "json",     
        timeout: 80000, // Time in milliseconds
        // The LOAD function will be called on a successful response.
        
        load: function(responseObject, ioArgs,evt){
			console.dir(responseObject);
            if (responseObject["status"] == "OK") {
                div.innerHTML = "<i class='fa fa-check-square f' style='font-size:13px;color:white'> SAVED </i>";
               td.classList.add("success");
				if (stamp) {
               		stamp.innerHTML = "validated";
               		stamp.style.visibility = "visible";
				}
				if (n){
               			n.style.visibility ="hidden";
				}
                 return;
            }       
            else  {
                div.innerHTML = "<i class='fa fa-exclamation-triangle style='font-size:13px;color:white''> ERROR</i>";
                td.classList.add("danger");
                 return;
            }   
        
            
        },
        
        // The ERROR function will be called in an error case.
        error: function(response, ioArgs){ // &acirc;&#382;&fnof;
            alert("Error !!!!");
               div.innerHTML = '<i class="fa fa-exclamation-triangle"> Error</i>';
           
            return;
        }
    });
  	return;


}

function enabled_disabled_buttons_view_interfaces() {
	document.getElementById("b_select_polyviewer").disabled = false;
	document.getElementById("b_select_polysplice").disabled = false;
	if (is_genome_project == 1) {
		document.getElementById("b_select_polycyto").disabled = false;
		document.getElementById("b_select_polydupdel").disabled = true;
	}
	else {
		document.getElementById("b_select_polycyto").disabled = true;
		document.getElementById("b_select_polydupdel").disabled = false;
	}
}

function hide_buttons_view_interfaces() {
	//document.getElementById("b_select_interfaces").style.display = "none";
	document.getElementById("b_select_polyviewer").style.display = "none";
	document.getElementById("b_select_polysplice").style.display = "none";
	document.getElementById("b_select_polycyto").style.display = "none";
	document.getElementById("b_select_polydupdel").style.display = "none";
}

function show_buttons_view_interfaces() {
	document.getElementById("span_legend_view_interface").style.display = "";
	document.getElementById("b_select_polyviewer").style.backgroundColor = "rgb(54, 57, 69)";
	document.getElementById("b_select_polyviewer").style.color = "rgb(127, 128, 131)";
	document.getElementById("b_select_polyviewer").style.display = "";
	
	if (view_polysplice == 1) {
		document.getElementById("b_select_polysplice").style.backgroundColor = "rgb(54, 57, 69)";
		document.getElementById("b_select_polysplice").style.color = "rgb(127, 128, 131)";
		document.getElementById("b_select_polysplice").style.display = "";
	}
	if (is_genome_project == 1) {	
		document.getElementById("b_select_polycyto").style.backgroundColor = "rgb(54, 57, 69)";
		document.getElementById("b_select_polycyto").style.color = "rgb(127, 128, 131)";
		document.getElementById("b_select_polycyto").style.display = "";
	}
	else {
		document.getElementById("b_select_polydupdel").style.backgroundColor = "rgb(54, 57, 69)";
		document.getElementById("b_select_polydupdel").style.color = "rgb(127, 128, 131)";
		document.getElementById("b_select_polydupdel").style.display = "";
	}
	if (is_diag_project == 1) {
		document.getElementById("b_select_polydiag").style.backgroundColor = "rgb(54, 57, 69)";
		document.getElementById("b_select_polydiag").style.color = "rgb(127, 128, 131)";
		document.getElementById("b_select_polydiag").style.display = "";
	}
	if (summary_panel_loaded == 2) {
		document.getElementById("td_b_select_polydiag").style.paddingLeft = "0px";
		document.getElementById("td_b_select_polyviewer").style.paddingLeft = "0px";
		document.getElementById("td_b_select_polycyto").style.paddingLeft = "0px";
		document.getElementById("td_b_select_polydupdel").style.paddingLeft = "0px";
		document.getElementById("td_b_select_polysplice").style.paddingLeft = "0px";
		document.getElementById("b_select_polydiag").style.display = "none";
		document.getElementById("b_select_polyviewer").style.display = "none";
		document.getElementById("b_select_polycyto").style.display = "none";
		document.getElementById("b_select_polydupdel").style.display = "none";
		require(["dojo/on", "dojo/dom", "dojo/ready", "dijit/registry","dojo/dom-construct", "dijit/layout/TabContainer", "dijit/layout/ContentPane","dojo/dom-style"], function(On, dom, ready, registry,domConstruct, TabContainer, ContentPane, domStyle) {
			ready(function() {
			    var tab1 = registry.byId("tabPrincipal_tablist_genes_list");
			    domStyle.set(tab1.domNode, "display", "none");
			    
			    var tab2 = registry.byId("tabPrincipal_tablist_coverage_stack");
			    domStyle.set(tab2.domNode, "display", "none");
			});
		});
	}
	dijit.byId('waiting').hide();
}

function cosmetics_click_button_project_resume() {
	document.getElementById("tabFilters").style.display = "none";
	if (summary_panel_loaded == 2) {
		document.getElementById("navMenuPatient").style.display = "none";
	}
	window.onresize();
}

function click_button_project_resume() {
	document.getElementById("tabPrincipal_tablist_tab_project").click();
	cosmetics_click_button_project_resume();
}

function cosmetics_click_button_coverage() {
	document.getElementById("tabFilters").style.display = "none";
	window.onresize();
}

function click_button_coverage() {
	document.getElementById("tabPrincipal_tablist_coverage_stack").click();
	cosmetics_click_button_coverage();
}

function click_button_view_polydupdel() {
	var id_pane_patient = "tabPrincipal_tablist_panel_" + selected_patient_name;
	document.getElementById(id_pane_patient).click();
	show_buttons_view_interfaces();
	var patient_id = hash_selected_patient_ids[selected_patient_name];
	var id_pane = "dijit_layout_TabContainer_" + patient_id + "_tablist_panel_polydupdel_" + selected_patient_name;
	document.getElementById(id_pane).click();
	show_buttons_view_interfaces();
	document.getElementById("b_select_polydupdel").style.backgroundColor = "white";
	document.getElementById("b_select_polydupdel").style.color = "black";
}

function click_button_view_polyviewer() {
	var id_pane_patient = "tabPrincipal_tablist_panel_" + selected_patient_name;
	document.getElementById(id_pane_patient).click();
	show_buttons_view_interfaces();
	document.getElementById('tabFilters').style.display = '';
	window.onresize();
	var patient_id = hash_selected_patient_ids[selected_patient_name];
	var id_pane = "dijit_layout_TabContainer_" + patient_id + "_tablist_panel_polyviewer_" + selected_patient_name;
	document.getElementById(id_pane).click();
	show_buttons_view_interfaces();
	document.getElementById("b_select_polyviewer").style.backgroundColor = "white";
	document.getElementById("b_select_polyviewer").style.color = "black";
}

function click_button_view_polysplice() {
	var id_pane_patient = "tabPrincipal_tablist_panel_" + selected_patient_name;
	document.getElementById(id_pane_patient).click();
	show_buttons_view_interfaces();
	document.getElementById('tabFilters').style.display = '';
	window.onresize();
	var patient_id = hash_selected_patient_ids[selected_patient_name];
	var id_pane = "dijit_layout_TabContainer_" + patient_id + "_tablist_panel_polysplice_" + selected_patient_name;
	document.getElementById(id_pane).click();
	document.getElementById("b_select_polysplice").style.backgroundColor = "white";
	document.getElementById("b_select_polysplice").style.color = "black";
}

function click_button_view_polycyto() {
	var id_pane_patient = "tabPrincipal_tablist_panel_" + selected_patient_name;
	document.getElementById(id_pane_patient).click();
	show_buttons_view_interfaces();
	document.getElementById('tabFilters').style.display = 'none';
	window.onresize();
	var patient_id = hash_selected_patient_ids[selected_patient_name];
	var id_pane = "dijit_layout_TabContainer_" + patient_id + "_tablist_panel_polycyto_" + selected_patient_name;
	document.getElementById(id_pane).click();
	document.getElementById("b_select_polycyto").style.backgroundColor = "white";
	document.getElementById("b_select_polycyto").style.color = "black";
}

function click_button_view_polydiag() {
	var id_pane_patient = "tabPrincipal_tablist_panel_" + selected_patient_name;
	document.getElementById(id_pane_patient).click();
	show_buttons_view_interfaces();
	document.getElementById('tabFilters').style.display = '';
	window.onresize();
	var patient_id = hash_selected_patient_ids[selected_patient_name];
	var id_pane = "dijit_layout_TabContainer_" + patient_id + "_tablist_panel_polydiag_" + selected_patient_name;
	document.getElementById(id_pane).click();
	document.getElementById("b_select_polydiag").style.backgroundColor = "white";
	document.getElementById("b_select_polydiag").style.color = "black";
}

function get_panels_from_db () {
	document.getElementById("span_phenotype_name").innerHTML = "No Phenotype Associated";
	origin_project_phenotype = 'No Phenotype Associated';
	hash_relation_panel_phenotype['all'] = [];
	hash_relation_panel_phenotype['all'].push(origin_project_phenotype);
	var was_phenotype_project = 0;
	store_tmp.fetch({
		onItem: function(item, result) {
			if (this_phenotype_name == String(item.phenotype_name)) {
				if (item.name in hPanelsFound) { }
				else {
					lPanels.push(item.name);
					lPanelsText.push(item.description);
					lGenesTableHtml.push(item.genes_table_html);
					hPanelsFound[item.name] = 'ok';
					hPanelsNbGenes[item.name] = item.nb_genes;
					hPanelsSource[item.name] = item.source;
					hPanelsCreator[item.name] = item.creator;
				}
			}
			else {
				this_phenotype_name = String(item.phenotype_name);
				//lPanels.push("<li class='divider'></li>");
				if (item.phenotype_project == 1) {
					was_phenotype_project = 1;
					document.getElementById("span_phenotype_name").innerHTML = this_phenotype_name;
					origin_project_phenotype = this_phenotype_name;
					hash_relation_panel_phenotype['all'] = [];
					hash_relation_panel_phenotype['all'].push(origin_project_phenotype);
				}
				else {
					was_phenotype_project = 0;
				}
				if (item.name in hPanelsFound) { }
				else {
					lPanels.push(item.name);
					lPanelsText.push(item.description);
					lGenesTableHtml.push(item.genes_table_html);
					hPanelsFound[item.name] = 'ok';
					hPanelsNbGenes[item.name] = item.nb_genes;
					hPanelsSource[item.name] = item.source;
					hPanelsCreator[item.name] = item.creator;
				}
				var li = document.createElement('li');
				var a = document.createElement('a');
				a.setAttribute("href", "#");
				a.setAttribute("onClick", "update_list_panels_name('" + this_phenotype_name + "');");
				var span = document.createElement('span');
				span.innerHTML = this_phenotype_name;
				a.appendChild(span);
				li.appendChild(a);
				document.getElementById('dropdown-menu-phenotypes').appendChild(li);
			}
			if (item.name in hash_relation_panel_phenotype) { }
			else { hash_relation_panel_phenotype[String(item.name)] = []; }
			hash_relation_panel_phenotype[String(item.name)].push(this_phenotype_name);
		}
	});
	for (var i = 0; i < lPanels.length; i++) {
		var panel_name = lPanels[i];
		var li = document.createElement('li');
		var a = document.createElement('a');
		var span_id = "span_panel_in_list_" + panel_name;
		span_id = span_id.replace(/ /g, '_');
		span_id = span_id.replace(/<label><u>Phenotype:_/, '');
		span_id = span_id.replace(/<\/u><\/label>/, '');
		a.setAttribute("href", "#");
		a.setAttribute("id", span_id);
		a.setAttribute("style", "font-size:9px;");
		a.setAttribute("onClick", "empty_gene_name_selection();update_panel_name('" + panel_name + "');");
		var span = document.createElement('span');
		if (panel_name == 'all') {
			span.innerHTML = 'All Genes';
		} //TODO: ici


		else {
			if (panel_name in hPanelsNbGenes) {
				var panel_text = lPanelsText[i];
				var source = hPanelsSource[panel_name];
				var creator = hPanelsCreator[panel_name];
				var nb_genes = hPanelsNbGenes[panel_name];
				if (nb_genes == 'all') { span.innerHTML = '<b>' + panel_name + '</b>'; }
				else {
					var text = panel_text;
					var table_genes_html = lGenesTableHtml[i];
					var b_genes_id = 'b_genes_' + panel_name + '_' + source;

					//var b_genes = "<button type='button' class='btn btn-info btn-xs' onmouseover=\"alert('" + table_genes_html + "')\" style='background-color:#BC243C;font-size:8px;margin:2px'>" + nb_genes + "</button>";						
					var b_genes = "<button id='" + b_genes_id + "'type='button' class='btn btn-info btn-xs' onmouseout=\"timer_mouse_exit_genes('" + b_genes_id + "', '" + nb_genes + "')\"  onmouseover=\"timer_mouse_over_genes('" + table_genes_html + "', '" + b_genes_id + "', '" + panel_name + "')\" style='background-color:#BC243C;font-size:8px;margin:2px'>" + nb_genes + "</button>";

					if (source == '') { text += ' ('; }
					else { text += ' (Source: ' + source + ' - '; }
					
					if (creator == '') { text += ''; }
					else { text += 'Creator: ' + creator + ' - '; }

					if (nb_genes > 1) { text += b_genes + ' genes)'; }
					else { text += b_genes + ' gene)'; }

					span.innerHTML = text;
				}
				var b_popup_id = 'b_popup_' + panel_name;
				//span.innerHTML += " <button id='" + b_popup_id + "' type='button' class='btn btn-xs btn-danger' data-toggle='popover' title='Test' onmouseover=\"$('#" + b_popup_id + "').popover('show')\" onmouseout=\"$('#" + b_popup_id + "').popover('hide')\" data-content='And here's some amazing content. It's very engaging. Right?'>Click to toggle popover</button>";
				//span.innerHTML += " <button id='" + b_popup_id + "' type='button' class='btn btn-xs btn-danger' data-toggle='popover' container='body' title='Test' onmouseover=\"$('#" + b_popup_id + "').popover('show')\" data-content='And here's some amazing content. It's very engaging. Right?'>Click to toggle popover</button>";
			}
			else { span.innerHTML = panel_name; }
		}
		a.appendChild(span);
		li.appendChild(a);
		document.getElementById('dropdown-menu-panels').appendChild(li);
	}
}

function check_first_patient_click() {
	if (summary_panel_loaded == 2) {
		setTimeout(function () {
			click_button_view_polysplice();
			document.getElementById("b_select_polyviewer").style.display = "none";
			document.getElementById("b_select_polydiag").style.display = "none";
			document.getElementById("b_select_polycyto").style.display = "none";
			document.getElementById("b_select_polydupdel").style.display = "none";
		}, 500);
	}
	else {
		setTimeout(function () {
			click_button_view_polyviewer();
		}, 500);
	}
	return;
}

var selected_patient_name = '';
var hash_fam_name_from_pat = {};
var hash_bam_files = {};
var hash_label_patient = {};
var hash_selected_patient_ids = {};
var hash_tabContainer_patients = {};
function load_panels_patients(items) {
	document.getElementById("tabFilters_tablist").style.display = "none";
	//document.getElementById("div_selector_interface").style.display = "";
	//document.getElementById("tabPrincipal_tablist_tab_project").style.fontSize = "15px";
	//document.getElementById("tabPrincipal_tablist_coverage_stack").style.fontSize = "15px";
    var z = 0;
    var argsPost = new Array();
    tsamples = items;
    tsamples.sort(function(a,b) {
        return b.child > a.child;
    });
	for (var i=0; i<tsamples.length; i++) { 
		var n = tsamples[i].label;
		if (tsamples[i].hasJonction == 1) { view_polysplice = 1; }
		if (tsamples[i].is_genome == 1) { is_genome_project = 1; }
		if (tsamples[i].is_diag == 1) { is_diag_project = 1; }
     	hash_selected_patient_ids[n] = i;
		var ic ='<img src="images/polyicons//iconfinder_baby-boy_s.png" style="width:60%" >';
     	ic ='<img src="images/polyicons/iconfinder_baby-girl_16.png" style="width:60%" >';
       
		if (tsamples[i].child == 0) {
			ic = '<img src="images/polyicons/icons8-person-24.png" style="width:60%" >';
			if (tsamples[i].sex == "F")   ic = '<img src="images/polyicons/icons8-person-female-24.png" style="width:60%">';
		}
		else {
			if (tsamples[i].sex == "F" && tsamples[i].status == "1"){
				ic ='<img src="images/polyicons//iconfinder_baby-girl_s.png" style="width:60%" >';
			}
			if (tsamples[i].sex == "F" && tsamples[i].status == "2"){
				ic ='<img src="images/polyicons//iconfinder_baby-girl_d.png" style="width:60%" >';
			}
			if (tsamples[i].sex == "M" && tsamples[i].status == "1"){
				ic ='<img src="images/polyicons//iconfinder_baby-boy_s.png" style="width:60%" >';
			}
			if (tsamples[i].sex == "M" && tsamples[i].status == "2"){
				ic ='<img src="images/polyicons//iconfinder_baby-boy_d.png" style="width:60%" >';
			}
		}
		var star ="";
		if (tsamples[i].validation_status == -99 ){
			star = "<span class=' glyphicon glyphicon-exclamation-sign' aria-hidden='true' style='position:absolute;top:1px;left:1px;color:#BA0BE0'></span>";
		}
		else if (tsamples[i].validation_status >= 3 ) {
			star = "<span class='glyphicon glyphicon-star' aria-hidden='true' style='position:absolute;top:1px;left:1px;color:#E74C3C'></span>";
		}
		else if (tsamples[i].validation_status >= 1 ) {
			star = "<span class='glyphicon glyphicon-star' aria-hidden='true' style='position:absolute;top:1px;left:1px;color:black'></span>";
		}

		if (!tab_editor[n]) {
	     	var patient_name = tsamples[i].label;
	     	var family_name = tsamples[i].fam;
	     	hash_bam_files[patient_name] = tsamples[i].bam;
	     	hash_fam_name_from_pat[patient_name] = family_name;
	     	
	     	if (selected_patient_name == '') { selected_patient_name = patient_name; }
			var title1;
			if (tsamples[i].label == tsamples[i].fam) { title1 = star+"<table><tr><td rowspan=2>"+ic+ "</td><td style='font-size:10px;'>"+tsamples[i].label+"</td></tr></table>"; }
			else { title1 = star+"<table><tr><td rowspan=2>"+ic+ "</td><td style='font-size:09px;font-weight: bold;'>"+tsamples[i].label+"</td></tr><tr>"+"<td ><span style='font-size:09px;'>"+tsamples[i].fam+"</small></td></tr></table>"; }
			
			hash_label_patient[patient_name] = title1;
			
			var title2 = ic;
			hash_label_patient[patient_name + '_short'] = title2;
			
			require(["dojo/_base/array","dijit/layout/TabContainer", "dijit/layout/ContentPane","dojo/domReady!"], function(array, TabContainer, ContentPane) {
				var t_patient =  new ContentPane({
					style:"font-size:12px;width:100%;overflow-x:auto;overflow-y:scroll;",
					title: title2,
					id:"panel_"+patient_name,
					jsId:"panel_"+patient_name,
					name:patient_name,
					preload: 0,
					onShow:  function() {
						show_buttons_view_interfaces();
						selected_patient_name = this.name;
						get_last_fiters_tab_view_for_patient(this.name);
						document.getElementById("tabFilters").style.display = "";
						//document.getElementById('div_selector_interface').style.display = '';
						window.onresize();
						check_first_patient_click();
					}
				});
				
				hash_tabContainer_patients[patient_name] = new TabContainer({
			        doLayout: false
			    }, 'tabContainer');
			    
				tab_editor_polyviewer[patient_name] = create_content_pane_polyviewer(patient_name);
				hash_tabContainer_patients[patient_name].addChild( tab_editor_polyviewer[patient_name] );
				
				tab_editor_polydiag[patient_name] = create_content_pane_polydiag(patient_name);
				hash_tabContainer_patients[patient_name].addChild( tab_editor_polydiag[patient_name] );
			    
				tab_editor_polysplice[patient_name] =  create_pane_polysplice(patient_name);
				hash_tabContainer_patients[patient_name].addChild( tab_editor_polysplice[patient_name] );
				
				if (is_genome_project == 1) {
					var t_polycyto = create_pane_polycyto(patient_name);
					tab_editor_polycyto[patient_name] = t_polycyto;
					hash_tabContainer_patients[patient_name].addChild(t_polycyto);
				}
				else {
					var t_polydupdel = create_pane_polydupdel(patient_name);
					tab_editor_polydude[patient_name] = t_polydupdel;
					hash_tabContainer_patients[patient_name].addChild(t_polydupdel);
				}
				
				t_patient.addChild(hash_tabContainer_patients[patient_name]);
				tabPrincipal.addChild(t_patient);
				
				document.getElementById("dijit_layout_TabContainer_" + i + "_tablist").style.display = "none";
				require(["dojo/on", "dojo/dom", "dojo/ready", "dijit/registry","dojo/dom-construct", "dijit/layout/TabContainer", "dijit/layout/ContentPane","dojo/dom-style"], function(On, dom, ready, registry,domConstruct, TabContainer, ContentPane, domStyle) {
  					ready(function() {
					    var tab1 = registry.byId("tabPrincipal_tablist_panel_" + patient_name);
					    domStyle.set(tab1.domNode, "display", "none");
					});
				});	
				
				ic = ic.replace('width:60%', 'width:16px;height:16px;');
				add_patient_name_in_dropdown_menu(family_name, patient_name, ic, star);
            });
        }
    }
    
    $('#tabPrincipal_tablist').on('mouseenter', function() {
    	show_patients_panel();
	});
    $('#tabPrincipal_tablist').on('mouseleave', function() {
		hide_patients_panel();
	});
	enabled_disabled_buttons_view_interfaces();
    write_all_sorted_patient_name_in_dropdown_menu();
	write_header_patients_infos(tsamples);
	write_header_project_infos();
    document.getElementById('input_search_patient').value = "";
	document.getElementById("div_selector_interface").style.display = "";
	document.getElementById("tabFilters").style.display = "none";
	window.onresize();
	
	if (is_diag_project == 1) {
		require(["dojo/on", "dojo/dom", "dojo/ready", "dijit/registry","dojo/dom-construct", "dijit/layout/TabContainer", "dijit/layout/ContentPane","dojo/dom-style"], function(On, dom, ready, registry,domConstruct, TabContainer, ContentPane, domStyle) {
			ready(function() {
			    var tab1 = registry.byId("tabPrincipal_tablist_coverage_stack_polydiag");
			    domStyle.set(tab1.domNode, "display", "");
			});
		});
		if (is_diag_project == 1) $.loadScript('js/polyjs/include_OLD.js', function() {
			load_genes_polydiag_tab(items);
		});
	}
	else {
		require(["dojo/on", "dojo/dom", "dojo/ready", "dijit/registry","dojo/dom-construct", "dijit/layout/TabContainer", "dijit/layout/ContentPane","dojo/dom-style"], function(On, dom, ready, registry,domConstruct, TabContainer, ContentPane, domStyle) {
			ready(function() {
			    var tab1 = registry.byId("tabPrincipal_tablist_genes_list");
			    domStyle.set(tab1.domNode, "display", "none");
			    
			    var tab2 = registry.byId("tabPrincipal_tablist_coverage_stack");
			    domStyle.set(tab2.domNode, "display", "");
			});
		});
	}
    dijit.byId('waiting').hide();
}

function write_header_patients_infos(tsamples) {
	document.getElementById('span_header_nb_patients').innerHTML = tsamples.length;
	document.getElementById('span_header_project_name').innerHTML = projectName;
}

function write_header_project_infos() {
	document.getElementById('img_logo_exome').style.display = "none";
	document.getElementById('img_logo_diag').style.display = "none";
	document.getElementById('img_logo_genome').style.display = "none";
	document.getElementById('img_logo_rna').style.display = "none";
	document.getElementById('img_logo_rna_exome').style.display = "none";
	document.getElementById('img_logo_rna_genome').style.display = "none";
	if (view_polysplice == 1) {
		if(summary_panel_loaded == 1) {
			if (is_genome_project == 1) {
				document.getElementById('img_logo_rna_genome').style.display = "";
			}
			else {
				document.getElementById('img_logo_rna_exome').style.display = "";
			}
		}
		else {
			document.getElementById('img_logo_rna').style.display = "";
		}
	}
	else {
		if (is_genome_project == 1) {
			document.getElementById('img_logo_genome').style.display = "";
		}
		else if (is_diag_project == 1) {
			document.getElementById('img_logo_diag').style.display = "";
		}
		else {
			document.getElementById('img_logo_exome').style.display = "";
		}
	}
}

var view_polydiag =0;
function show_low(){
	 var url = url_low+"?project="+projectName+"&utr=0&padding=5&covlimit=30";
	 dijit.byId('waiting').show();
	 document.getElementById("cov50_0").innerHTML = "";
	 $.ajax({ url: url, success: function(data) {
	 	document.getElementById("cov50_0").innerHTML = data;
	 	dijit.byId('waiting').hide();
	 }});
}

var tab_coverage = {};
var first =0;
var panel_coverage_stack_changed = 0;
function preload_coverage_polydiag(items,reload_flag) { 
	if (first >0) {
		reload(5);
		return;
	}
	first ++;
	tsamples = dijit.byId("gridPatients").selection.getSelected();
	var url_cnv = construct_url(url_cnv_overview)+"&test_cache=1";
	require(["dojo/promise/all", "dojo/Deferred", "dojo/request", "dojo/_base/array", "dojo/dom-construct", "dojo/dom", "dojo/json", "dojo/domReady!"], function(all, Deferred, request, arrayUtil, domConstruct, dom, JSON) {
		var p="";
		var samples = new Array();
		var requests = new Array();
		var url = construct_url(url_cached_coverage);
		var progress_id =0;
		var max =0;
		var i,j,temparray,chunk = 4;
		var labels = new Array;
		var z =0;
		for (i=0,j=tsamples.length; i<j; i+=chunk) {
			temparray = tsamples.slice(i,i+chunk);
			var temp = new Array();
			for (var x=0;x<temparray.length;x++){
				temp[x] = temparray[x].label;
			}
			labels[z++] = temp.join(",");
		}
		max = labels.length +1;
		for(i=0;i<Math.round(labels.length/2);i++){
			requests[i] = request.get( url+"&patient_cached="+ labels[i], { handleAs: "json" }).then(function(response){
				progress_id ++;
				myProgressBar.set( { value: progress_id } );
				return 1;
			}) ;
		}           
		myProgressBar.set( { value: 0,maximum:max } );
		dijit.byId('waiting_progress').show();   
    
		var request2 = new Array();
		all(requests).then(function(results){
			for(i=Math.round(labels.length/2);i<labels.length;i++){
				request2[i+labels.length] = request.get( url+"&no_cache=1&patient_cached="+ labels[i], { handleAs: "json" }).then(function(response){
					progress_id ++;
					myProgressBar.set( { value: progress_id } );
					return 1;
				}) ;
			}           
			request2[i+labels.length] = request.get(url_cnv, { handleAs: "html" }).then(function(response){
				progress_id ++;
				myProgressBar.set( { value: progress_id } );
				return 1;
			});
			all(request2).then(function(results){
				if (reload_flag ==1) reload(5);
				dijit.byId('waiting_progress').hide();
			});
		});      
	});
}

function add_random_bam_igv(type_project) {
	return get_random_bam_file('igv', type_project);
}

function add_random_bam_alamut(type_project) {
	return get_random_bam_file('alamut', type_project);
}

function get_random_bam_file(tool, type_project) {
	var bam_url;
	var this_url = url_get_random_bam_file;
	var this_args = "type_project=" + type_project;
	var xhr = new XMLHttpRequest();
	xhr.open('POST', this_url, true);
	xhr.setRequestHeader('Content-type', 'application/x-www-form-urlencoded');
	xhr.onreadystatechange = function() {
		if (this.readyState === XMLHttpRequest.DONE && this.status === 200) {
			var data = JSON.parse(this.responseText);
			bam_url = '/' + data.url_bam;
			console.log(bam_url);
			if (tool == 'igv') {
				var array = [];
				array.push(bam_url)
				addListBamInIgvApp(array);
			}
			if (tool == 'alamut') {
				httpGetLoadOnlyListBam(bam_url);
			}
		}
	}
	xhr.send(this_args);
}

function add_random_patient_name_in_dropdown_menu(type_project) {
	var patient_name = 'random_' + type_project;
	var label_text = "&nbsp;Add Random <b>" + String(type_project).toUpperCase() + "</b>";
	
	// button select bam IGV
	var elem_igv = document.createElement("span");
	elem_igv.setAttribute("id", "span_igv_list_patient_"+patient_name);
	elem_igv.setAttribute("style", "color:black;");
	var html_igv_button = "<a class='dropdown-item text-center' style='color:black;font-size:13px;'><button style='margin:1px;border:solid 1px white;background-color:white;color:black !important;' onClick='add_random_bam_igv(\"" + type_project + "\")'>" + label_text + "</button></a>";
	elem_igv.innerHTML = html_igv_button + "<br>";
	array_list_label_patients_igv[patient_name] = elem_igv;
	
	// button select bam ALAMUT
	var elem_alamut = document.createElement("span");
	elem_alamut.setAttribute("id", "span_alamut_list_patient_"+patient_name);
	elem_alamut.setAttribute("style", "color:black;");
	var html_alamut_button = "<a class='dropdown-item text-center' style='color:black;font-size:13px;'><button style='margin:1px;border:solid 1px white;background-color:white;color:black !important;' onClick='add_random_bam_alamut(\"" + type_project + "\")'>" + label_text + "</button></a>";
	elem_alamut.innerHTML = html_alamut_button + "<br>";
	array_list_label_patients_alamut[patient_name] = elem_alamut;
}

var hash_patient_selected_html = {};
var array_list_label_patients = [];
var array_list_label_patients_div_id = [];
var array_list_label_patients_igv = [];
var array_list_label_patients_alamut = [];
function add_patient_name_in_dropdown_menu(family_name, patient_name, ic, star) {
	if (star) { star = star.replace('position:absolute;top:1px;left:1px;', ''); }
	var patient_name_text = patient_name + "&nbsp;&nbsp;";
	if (String(patient_name).length > 10) {
		patient_name_text = String(patient_name).substring(0, 8) + "...";
	}
	var label_text = "&nbsp;&nbsp;" + ic + " <b>" + family_name + "</b> - " + patient_name_text;
	hash_patient_selected_html[patient_name] = label_text;
	var bam_file = hash_bam_files[patient_name];
	
	// button select patient
	var elem = document.createElement("span");
	elem.setAttribute("id", "span_list_patient_"+patient_name);
	elem.setAttribute("style", "color:black;");
	var display_b_select = "display:none;";
	display_b_select = "";
//	if (is_diag_project == 1) { display_b_select = ""; }
	var html_select_button = "<span class='form-check' style='" + display_b_select + "'><input style='' class='form-check-input' onclick='event.stopPropagation()' type='checkbox' value='' id='b_selected_patient_" + patient_name + "' checked></span>";
	var b_id = "b_choice_patient_" + patient_name;
	var html_text;
	if (family_name && family_name != '-') {
		var length_pat_fam_name = String(family_name).length + String(patient_name).length;
		if (length_pat_fam_name > 25) {
			html_text = "<a class='dropdown-item text-center' style='color:black;font-size:10px !important;text-align:left;'><table onclick=\"select_patient_in_dropdown_menu('" + patient_name + "');\"><td>&nbsp;" + html_select_button + "</td><td>" + ic + " </td><td><b> " + family_name + "</b><br><button id='" + b_id + "' class='btn btn-xs' style='margin:1px;border:solid 1px white;background-color:white;color:black !important;font-size:10px !important;text-align:left;'><span> " + patient_name + " " + star + "</span></td></table></button></a>";
		}
		else {
			html_text = "<a class='dropdown-item text-center' style='color:black;font-size:10px !important;text-align:left;'><table onclick=\"select_patient_in_dropdown_menu('" + patient_name + "');\"><td>&nbsp;" + html_select_button + "</td><td>" + ic + " </td><td><b> " + family_name + "</b> <span style='color:#363945;' class='glyphicon glyphicon-minus'></span> <button id='" + b_id + "' class='btn btn-xs' style='margin:1px;border:solid 1px white;background-color:white;color:black !important;font-size:10px !important;text-align:left;'><span> " + patient_name + " " + star + "</span></td></table></button></a>";
		}
	}
	else {
		html_text = "<a class='dropdown-item text-center' style='color:black;font-size:10px !important;text-align:left;'><table onclick=\"select_patient_in_dropdown_menu('" + patient_name + "');\"><td>&nbsp;" + html_select_button + "</td><td>" + ic + " </td><td> <button id='" + b_id + "' class='btn btn-xs' style='margin:1px;border:solid 1px white;background-color:white;color:black !important;font-size:10px !important;text-align:left;'>" + patient_name + " " + star + "</button></a>";
	}
	elem.innerHTML = html_text;
	array_list_label_patients[family_name + '_' + patient_name] = elem;
	array_list_label_patients_div_id[family_name + '_' + patient_name] = "div_list_patient_"+patient_name;
	array_list_label_patients_div_id[patient_name] = "div_list_patient_"+patient_name;
	
	// button select bam IGV
	var elem_igv = document.createElement("span");
	elem_igv.setAttribute("id", "span_igv_list_patient_"+patient_name);
	elem_igv.setAttribute("style", "color:black;");
	var html_igv_button;
	if (family_name && family_name != '-') {
		html_igv_button = "<a class='dropdown-item text-center' style='color:black;font-size:13px;'><button style='margin:1px;border:solid 1px white;background-color:white;color:black !important;' onClick ='addListBamInIgvApp([\"" + bam_file + "\"])';>" + ic + " <b>" + family_name + "</b> - " + patient_name + " " + star + "</button></a>";
	}
	else {
		html_igv_button = "<a class='dropdown-item text-center' style='color:black;font-size:13px;'><button style='margin:1px;border:solid 1px white;background-color:white;color:black !important;' onClick ='addListBamInIgvApp([\"" + bam_file + "\"])';>" + ic + "&nbsp;&nbsp;" + patient_name + " " + star + "</button></a>";
	}	
	elem_igv.innerHTML = html_igv_button + "<br>";
	array_list_label_patients_igv[family_name + '_' + patient_name] = elem_igv;
	
	// button select bam ALAMUT
	var elem_alamut = document.createElement("span");
	elem_alamut.setAttribute("id", "span_alamut_list_patient_"+patient_name);
	elem_alamut.setAttribute("style", "color:black;");
	var html_alamut_button;
	if (family_name && family_name != '-') {
		html_alamut_button = "<a class='dropdown-item text-center' style='color:black;font-size:13px;'><button style='margin:1px;border:solid 1px white;background-color:white;color:black !important;' onClick ='httpGetLoadOnlyListBam(\"" + bam_file + "\")';>" + ic + " <b>" + family_name + "</b> - " + patient_name + " " + star + "</button></a>";
	}
	else {
		html_alamut_button = "<a class='dropdown-item text-center' style='color:black;font-size:13px;'><button style='margin:1px;border:solid 1px white;background-color:white;color:black !important;' onClick ='httpGetLoadOnlyListBam(\"" + bam_file + "\")';>" + ic + "&nbsp;&nbsp;" + patient_name + " " + star + "</button></a>";
	}
	elem_alamut.innerHTML = html_alamut_button + "<br>";
	array_list_label_patients_alamut[family_name + '_' + patient_name] = elem_alamut;
}

function select_all_patients() {
	for (name in hash_bam_files) {
		var b_id = "b_selected_patient_" + name;
		document.getElementById(b_id).checked = true;
	}
}

function unselect_all_patients() {
	for (name in hash_bam_files) {
		var b_id = "b_selected_patient_" + name;
		document.getElementById(b_id).checked = false;
	}
}

function hide_select_all_patients() {
	for (name in hash_bam_files) {
		var b_id = "b_selected_patient_" + name;
		document.getElementById(b_id).checked = true;
		document.getElementById(b_id).style.display = 'none';
	}
}

function filter_list_dropdown_patients() {
	var text_to_search = document.getElementById('input_search_patient').value;
	var nb_found = 0;
	var nb_total = 0;
	const reguarExp = new RegExp(text_to_search, 'i');
	for (name in hash_bam_files) {
		var fam = hash_fam_name_from_pat[name];
		var fam_pat_name = fam + '_' + name;
		nb_total++;
		var display_value = 'none';
		if (text_to_search.length < 1) {
			display_value = '';
			nb_found++;
		}
		else if (reguarExp.test(fam_pat_name)) {
			display_value = '';
			nb_found++;
		}
		document.getElementById(array_list_label_patients_div_id[name]).style.display = display_value;
		document.getElementById('span_igv_list_patient_' + name).style.display = display_value;
		document.getElementById('span_alamut_list_patient_' + name).style.display = display_value;
	}
	var text_nb_found;
	if (nb_found == nb_total) { document.getElementById('span_header_nb_patients').innerHTML = nb_total; }
	else { document.getElementById('span_header_nb_patients').innerHTML = nb_found + "/" + nb_total; }
}

function write_all_sorted_patient_name_in_dropdown_menu() {
    add_random_patient_name_in_dropdown_menu('genome');
	add_random_patient_name_in_dropdown_menu('exome');
	document.getElementById("b_menu_dropdown_list_patients_bam_igv").appendChild( array_list_label_patients_igv["random_exome"] );
	document.getElementById("b_menu_ddropdown_list_patients_bam_alamut").appendChild( array_list_label_patients_alamut["random_exome"] );
	document.getElementById("b_menu_dropdown_list_patients_bam_igv").appendChild( array_list_label_patients_igv["random_genome"] );
	document.getElementById("b_menu_ddropdown_list_patients_bam_alamut").appendChild( array_list_label_patients_alamut["random_genome"] );
	var max_length_name = 0;
	var array = [];
	for (name in array_list_label_patients) {
		array.push(name);
		var name_filter = String(name).replace('span_list_patient_', '');
		var this_length = String(name_filter).length;
		if (max_length_name < this_length) { max_length_name = this_length; }
	}
	array.sort();
	
	var length_div_row = 2;
	if (max_length_name > 8) { length_div_row = 3; }
	if (max_length_name > 18) { length_div_row = 4; }
	
	var div = document.createElement("div");
	div.setAttribute("style", "width:725px;overflow-x:hidden;");
	var div_row = document.createElement("div");
	div_row.setAttribute("class", "row");
	div.appendChild( div_row );
	for (var i = 0; i < array.length; i++) {
		var name = array[i];
		var div_col = document.createElement("div");
		var div_pat_id = array_list_label_patients_div_id[name];
		div_col.setAttribute("id", div_pat_id);
		div_col.setAttribute("class", "col-sm-" + length_div_row);
		div_col.appendChild( array_list_label_patients[name] );
		div_row.appendChild( div_col );
		document.getElementById("b_menu_dropdown_list_patients_bam_igv").appendChild( array_list_label_patients_igv[name] );
		document.getElementById("b_menu_ddropdown_list_patients_bam_alamut").appendChild( array_list_label_patients_alamut[name] );
	}
	document.getElementById("b_menu_dropdown_list_patients").appendChild( div );
	setTimeout(function () {
		if (is_diag_project == 0) { hide_select_all_patients(); }
		else {
			document.getElementById("b_select_all_pat").style.display = '';
			document.getElementById("b_unselect_all_pat").style.display = '';
		}
		document.getElementById("span_patient_selected_text").innerHTML = "&nbsp;Select Your Sample <b>HERE <span class='glyphicon glyphicon-option-vertical'></span></b>&nbsp;</b>";
		dijit.byId('waiting').hide();
	}, 1000);
	return;
}

var used_patient_in_dropdown_menu;
function select_patient_in_dropdown_menu(patient_name) {
	used_patient_in_dropdown_menu = patient_name;
	for (name in hash_tabContainer_patients) {
		require(["dojo/on", "dojo/dom", "dojo/ready", "dijit/registry","dojo/dom-construct", "dijit/layout/TabContainer", "dijit/layout/ContentPane","dojo/dom-style"], function(On, dom, ready, registry,domConstruct, TabContainer, ContentPane, domStyle) {
			ready(function() {
			    var tab1 = registry.byId("tabPrincipal_tablist_panel_" + name);
//			    domStyle.set(tab1.domNode, "display", "none");
			});
		});	
	}
	require(["dojo/on", "dojo/dom", "dojo/ready", "dijit/registry","dojo/dom-construct", "dijit/layout/TabContainer", "dijit/layout/ContentPane","dojo/dom-style"], function(On, dom, ready, registry,domConstruct, TabContainer, ContentPane, domStyle) {
		ready(function() {
		    var tab1 = registry.byId("tabPrincipal_tablist_panel_" + patient_name);
		    domStyle.set(tab1.domNode, "display", "");
		});
	});	
	
	var html_text = "&nbsp;" + hash_patient_selected_html[patient_name] + "&nbsp;<b><span class='glyphicon glyphicon-option-vertical'></span></b>&nbsp;</b>";
	document.getElementById("span_patient_selected_text").innerHTML = html_text;
	
	var panel_name = "tabPrincipal_tablist_panel_" + patient_name;
	document.getElementById(panel_name).click();
}

jQuery.loadScript = function (url, callback) {
    jQuery.ajax({
        url: url,
        dataType: 'script',
        success: callback,
        async: true
    });
}

function fillTabPatientsDiag(items) {
	nb_patients = items.length;
	preload_coverage_polydiag(items,0);
}  

var layoutGridGenesPolydiag = [
	{ width: '10%', field: 'label', name: 'Name', styles: 'text-align:center;'}, 
	{ width: '15%', field: 'transcript', name: 'Transcript', styles: 'text-align:center;'}, 
	{ width: '5%', field: 'chr', name: 'Chr', styles: 'text-align:center;'}, 
	{ width: '7%', field: 'start', name: 'Start', styles: 'text-align:center;'}, 
	{ width: '7%', field: 'end', name: 'End', styles: 'text-align:center;'}, 
	{ width: '10%', field: 'refseq', name: 'Refseq', styles: 'text-align:center;'}, 
	{ width: '10%', field: 'ccds', name: 'CCDS', styles: 'text-align:center;'}, 
	{ width: '5%', field: 'nb_exons', name: 'Exon', styles: 'text-align:center;'}, 
	{ width: '10%', field: 'bundle', name: 'Group', styles: 'text-align:left;'}, 
	{ width: '20%', field: 'description', name: 'Description', styles: 'text-align:left;'}, 
];

function load_genes_polydiag_tab(items) {
	storeP = items;
	fillTabPatientsDiag(items);
	tabPatients(items);
	var grid2 = new dojox.grid.EnhancedGrid({
    	id: 'gridGenes',
        structure: layoutGridGenesPolydiag,
        selection:{mode:'single'}, 
        plugins: {indirectSelection: {headerSelector:true,  field: 'active',name: 'Select', width:"40px", styles:"text-align: center;",selected:true}}
	}, document.createElement('div'));
	dojo.byId("gridDivGenes").appendChild(grid2.domNode);
	grid2.startup();
	var url2 = construct_url(url_validation_genes);
	storeG = new dojo.data.ItemFileWriteStore({ url: url2 });
	storeG.fetch({onComplete: fillTabGenes});	
	storeG.fetch({onComplete: construct_polydiag_bundle_tree()});
	return;
}

var tree;
function construct_polydiag_bundle_tree() {
    if (!tree) {
        require([
            "dojo/ready",
            "dojo/data/ItemFileWriteStore",
            "cbtree/Tree",                    // Checkbox tree
            "cbtree/models/ForestStoreModel"  // ForestStoreModel
        ], function( ready, ItemFileWriteStore, Tree, ForestStoreModel) {
            var url_bundles = url_cache_gene_vector + "?get_bundles=1&project=" + projectName;
            ready( function() {
                var store_bundles = new ItemFileWriteStore( { url: url_diag_tree+"?project="+project_name+"&panel="+panel });
                var model_bundles = new ForestStoreModel( {   store: store_bundles,
                                                              query: {type: "group"},
                                                              rootLabel: "Groups",
                                                              checkedAll: true,
                                                          });
                tree = new Tree( { model: model_bundles, id: "bundleTree", openOnClick:true, showRoot:true }, "bundle_div" );
                tree.on( "checkBoxClick", polydiag_bundle_tree_change_selected_genes);
                tree.startup();
            });
        });
    }
}
    
function polydiag_bundle_tree_change_selected_genes(item, nodeWidget, event ) {
	var t = new Array();
	if (item.children) {
		for (var l in item.children) {
			var group = item.children[l];
			if (group.children){
				for (var z in group.children){
					var leaf = group.children[z];
					t[leaf.name] = leaf.checked;
				}
			}
			else {
				t[group.name] = group.checked;
			}
		}  
	}	
	else { t[item.name] = item.checked; }
	for (var d in tree.model.root.children){
		var ds=  tree.model.root.children[d];
		for (var l in ds.children){
			var leaf = ds.children[l];
			if (t[leaf.name]){
				if (t[leaf.name]=="false"){
					tree.model.setChecked(leaf,false);
				}
				else {
					tree.model.setChecked(leaf,true);
				}
			}
		}
	}	
	var grid =  dijit.byId("gridGenes");
	for (i = 0; i < grid.rowCount; i++){
		var it = grid.getItem(i);
	    if (t[it.label]){
			var value = t[it.label];
			if (value=="false"){
				grid.rowSelectCell.toggleRow(i, false) ;
			}
			else {
				grid.rowSelectCell.toggleRow(i, true) ;
			}
		}
	}
}	

function show_patients_panel() {
	document.getElementById('tabPrincipal_tablist_tab_project').innerHTML = "<span class='glyphicon glyphicon-home'></span> Project";
	document.getElementById('tabPrincipal_tablist_coverage_stack').innerHTML = "<span class='glyphicon glyphicon-signal'></span> Coverage";
	document.getElementById('tabPrincipal_tablist_coverage_stack_polydiag').innerHTML = "<span class='glyphicon glyphicon-signal'></span> Coverage Diag";
	document.getElementById('tabPrincipal_tablist_genes_list').innerHTML = "<span class='glyphicon glyphicon-list-alt'></span> Genes List";
	for (name in hash_selected_patient_ids) {
		var div_name = "tabPrincipal_tablist_panel_" + name;
		document.getElementById(div_name).innerHTML = hash_label_patient[name];
	}
	window.onresize();
}

function hide_patients_panel() {
	document.getElementById('tabPrincipal_tablist_tab_project').innerHTML = "<span class='glyphicon glyphicon-home'></span>";
	document.getElementById('tabPrincipal_tablist_coverage_stack').innerHTML = "<span class='glyphicon glyphicon-signal'></span>";
	document.getElementById('tabPrincipal_tablist_coverage_stack_polydiag').innerHTML = "<span class='glyphicon glyphicon-signal'></span>";
	document.getElementById('tabPrincipal_tablist_genes_list').innerHTML = "<span class='glyphicon glyphicon-list-alt'></span>";
	for (name in hash_selected_patient_ids) {
		var div_name = "tabPrincipal_tablist_panel_" + name;
		document.getElementById(div_name).innerHTML = hash_label_patient[name + "_short"];
	}
	window.onresize();
}

function select_panel_filters_named(name) {
	document.getElementById("div_header_logos").disabled = true;
	document.getElementById("div_filters_polyviewer").disabled = true;
	document.getElementById("div_filters_polysplice").disabled = true;
	document.getElementById("div_filters_polydupdel").disabled = true;
	document.getElementById("div_filters_polydiag").disabled = true;
	document.getElementById(name).disabled = false; 
	dijit.byId('tabFilters').selectChild(name, true);
}

function select_panel_patients_named(name) {
	dijit.byId('tabPrincipal').selectChild(name, true);
}
var hash_last_interface_selected_patients = {};
function get_last_interface_view_for_patient(patient_name) {
	for (var key in hash_last_interface_selected_patients) {
		var value = hash_last_interface_selected_patients[key];
		if (key == patient_name) { return value; }
	}
	return;
}

function update_last_interface_view_for_patient(patient_name, interface_id) {
	hash_last_interface_selected_patients[patient_name] = interface_id;
}

function get_last_fiters_tab_view_for_patient(patient_name) {
	var interface_id = get_last_interface_view_for_patient(patient_name);
	if (interface_id) {
		if (interface_id == 'panel_polyviewer_' + patient_name) { select_panel_filters_named('div_filters_polyviewer'); }
		if (interface_id == 'panel_polysplice_' + patient_name) { select_panel_filters_named('div_filters_polysplice'); }
		if (interface_id == 'panel_polydupdel_' + patient_name) { select_panel_filters_named('div_filters_polydupdel'); }
		if (interface_id == 'panel_polycyto_' + patient_name)   {  }
		if (interface_id == 'panel_polydiag_' + patient_name)   {  }
	}
}

function httpGet(theUrl) {
    var xmlHttp = new XMLHttpRequest();
    xmlHttp.open( "GET", theUrl, false ); // false for synchronous request
    xmlHttp.send( null );
    var html = xmlHttp.responseText;
    return html;
}

function create_content_pane_polyviewer(patient_name) {
	var t_polyviewer;
	require(["dojo/_base/array","dijit/layout/TabContainer", "dijit/layout/ContentPane","dojo/domReady!"], function(array, TabContainer, ContentPane) {
		t_polyviewer =  new ContentPane({
			style:"font-size:9px;max-height:800px;width:100%;overflow:auto scroll;",
			title: 'PolyViewer',
			id:"panel_polyviewer_"+patient_name,
			jsId:"panel_polyviewer_"+patient_name,
			preload:0,
			ioMethod:dojo.xhrPost,
			toto:1,
			patient:patient_name,
			onShow:  function() {
				update_last_interface_view_for_patient(patient_name, this.id);
				select_panel_filters_named('div_filters_polyviewer');
				document.getElementById("tabFilters").style.height = '145px';
				window.onresize();
				this.toto ++;
				var obj = this;
				if(this.toto == 2) {
					dijit.byId('waiting').show();
					obj.setArgsPost = null;
				}
				setTimeout(function () {
					if (launch_panels_from_db_done == 0) { launch_panels_from_db(); }
					if (!(dijit.byId("check_never_pv"))){ return; }
					if (obj.toto == 2 ){
						if (obj.setArgsPost == null ){ obj.setArgsPost = load_polyviewer_args(); }
						obj.setArgsPost.patients = patient_name;
						obj.ioArgs.content = obj.setArgsPost; 
						obj.setHref(url_polyviewer);
						dijit.byId('waiting').hide();
					}
				}, 2000);
			}
		});
	});
	
	return t_polyviewer;
}

var launch_panels_from_db_done = 0;
function launch_panels_from_db() {
	if (launch_panels_from_db_done == 1) { return; }
	var lButtonGeneSelect = [];
	lButtonGeneSelect.push("<form class='px-4 py-3'><div class='form-group'><label for='dropdown_gene_name'><u>Gene Name</u></label><input id='gene_name_editor' type='text' class='form-control'></div><center><button class='btn btn-xs btn-success' style='width:100%;' onclick='change_gene_name_selection()' type='button'>View Only Gene</button><br></center><center><button class='btn btn-xs btn-warning' style='width:100%;' onclick='view_panels_gene()' type='button'>View Phenotypes / Panels</button><br></center></form>");
	lButtonGeneSelect.push("<li class='divider'></li>");
	lButtonGeneSelect.push("all");
	lButtonGeneSelect.push("<li class='divider'></li>");
	for (var i = 0; i < lButtonGeneSelect.length; i++) {
		var li = document.createElement('li');
		var a = document.createElement('a');
		a.setAttribute("href", "#");
		if (lButtonGeneSelect[i] == 'all') {
			document.getElementById("span_panel_name").innerHTML = "All Genes";
			a.setAttribute("onClick", "empty_gene_name_selection();update_list_panels_name('all');update_panel_name('all');");
		}
		var span = document.createElement('span');
		if (lButtonGeneSelect[i] == 'all') { span.innerHTML = 'All Genes'; }
		else { span.innerHTML = lButtonGeneSelect[i]; }
		a.appendChild(span);
		li.appendChild(a);
		document.getElementById('dropdown-menu-phenotypes').appendChild(li);
	}

	var hPanelsFound = {};
	var hPanelsNbGenes = {};
	var hPanelsSource = {};
	var hPanelsCreator = {};
	var lPanels = [];
	var lPanelsText = [];
	var lGenesTableHtml = [];
	var this_phenotype_name = 'all';
	var this_url = url_path + "/json_output_nodb/getPanelsFromDb.pl";
	var this_args = "project=" + project_name;
	var xhr = new XMLHttpRequest();
	xhr.open('POST', this_url, true);
	xhr.setRequestHeader('Content-type', 'application/x-www-form-urlencoded');
	xhr.onreadystatechange = function() {
		if (this.readyState === XMLHttpRequest.DONE && this.status === 200) {
			var data = JSON.parse(this.responseText);
			var store_tmp = new dojo.data.ItemFileWriteStore({
				data: {
					identifier: 'id',
					label: 'id',
					items: data.items
				}
			});
			document.getElementById("span_phenotype_name").innerHTML = "No Phenotype Associated";
			origin_project_phenotype = 'No Phenotype Associated';
			hash_relation_panel_phenotype['all'] = [];
			hash_relation_panel_phenotype['all'].push(origin_project_phenotype);
			var was_phenotype_project = 0;
			store_tmp.fetch({
				onItem: function(item, result) {
					if (this_phenotype_name == String(item.phenotype_name)) {
						if (item.name in hPanelsFound) { }
						else {
							lPanels.push(item.name);
							lPanelsText.push(item.description);
							lGenesTableHtml.push(item.genes_table_html);
							hPanelsFound[item.name] = 'ok';
							hPanelsNbGenes[item.name] = item.nb_genes;
							hPanelsSource[item.name] = item.source;
							hPanelsCreator[item.name] = item.creator;
						}
					}
					else {
						this_phenotype_name = String(item.phenotype_name);
						//lPanels.push("<li class='divider'></li>");
						if (item.phenotype_project == 1) {
							was_phenotype_project = 1;
							document.getElementById("span_phenotype_name").innerHTML = this_phenotype_name;
							origin_project_phenotype = this_phenotype_name;
							hash_relation_panel_phenotype['all'] = [];
							hash_relation_panel_phenotype['all'].push(origin_project_phenotype);
						}
						else {
							was_phenotype_project = 0;
						}
						if (item.name in hPanelsFound) { }
						else {
							lPanels.push(item.name);
							lPanelsText.push(item.description);
							lGenesTableHtml.push(item.genes_table_html);
							hPanelsFound[item.name] = 'ok';
							hPanelsNbGenes[item.name] = item.nb_genes;
							hPanelsSource[item.name] = item.source;
							hPanelsCreator[item.name] = item.creator;
						}
						var li = document.createElement('li');
						var a = document.createElement('a');
						a.setAttribute("href", "#");
						a.setAttribute("onClick", "update_list_panels_name('" + this_phenotype_name + "');");
						var span = document.createElement('span');
						span.innerHTML = this_phenotype_name;
						a.appendChild(span);
						li.appendChild(a);
						document.getElementById('dropdown-menu-phenotypes').appendChild(li);
					}
					if (item.name in hash_relation_panel_phenotype) { }
					else { hash_relation_panel_phenotype[String(item.name)] = []; }
					hash_relation_panel_phenotype[String(item.name)].push(this_phenotype_name);
				}
			});
			for (var i = 0; i < lPanels.length; i++) {
				var panel_name = lPanels[i];
				var li = document.createElement('li');
				var a = document.createElement('a');
				var span_id = "span_panel_in_list_" + panel_name;
				span_id = span_id.replace(/ /g, '_');
				span_id = span_id.replace(/<label><u>Phenotype:_/, '');
				span_id = span_id.replace(/<\/u><\/label>/, '');
				a.setAttribute("href", "#");
				a.setAttribute("id", span_id);
				a.setAttribute("style", "font-size:9px;");
				a.setAttribute("onClick", "empty_gene_name_selection();update_panel_name('" + panel_name + "');");
				var span = document.createElement('span');
				if (panel_name == 'all') {
					span.innerHTML = 'All Genes';
				} //TODO: ici


				else {
					if (panel_name in hPanelsNbGenes) {
						var panel_text = lPanelsText[i];
						var source = hPanelsSource[panel_name];
						var creator = hPanelsCreator[panel_name];
						var nb_genes = hPanelsNbGenes[panel_name];
						if (nb_genes == 'all') { span.innerHTML = '<b>' + panel_name + '</b>'; }
						else {
							var text = panel_text;
							var table_genes_html = lGenesTableHtml[i];
							var b_genes_id = 'b_genes_' + panel_name + '_' + source;

							//var b_genes = "<button type='button' class='btn btn-info btn-xs' onmouseover=\"alert('" + table_genes_html + "')\" style='background-color:#BC243C;font-size:8px;margin:2px'>" + nb_genes + "</button>";						
							var b_genes = "<button id='" + b_genes_id + "'type='button' class='btn btn-info btn-xs' onmouseout=\"timer_mouse_exit_genes('" + b_genes_id + "', '" + nb_genes + "')\"  onmouseover=\"timer_mouse_over_genes('" + table_genes_html + "', '" + b_genes_id + "', '" + panel_name + "')\" style='background-color:#BC243C;font-size:8px;margin:2px'>" + nb_genes + "</button>";

							if (source == '') { text += ' ('; }
							else { text += ' (Source: ' + source + ' - '; }
							
							if (creator == '') { text += ''; }
							else { text += 'Creator: ' + creator + ' - '; }

							if (nb_genes > 1) { text += b_genes + ' genes)'; }
							else { text += b_genes + ' gene)'; }

							span.innerHTML = text;
						}
						var b_popup_id = 'b_popup_' + panel_name;
						//span.innerHTML += " <button id='" + b_popup_id + "' type='button' class='btn btn-xs btn-danger' data-toggle='popover' title='Test' onmouseover=\"$('#" + b_popup_id + "').popover('show')\" onmouseout=\"$('#" + b_popup_id + "').popover('hide')\" data-content='And here's some amazing content. It's very engaging. Right?'>Click to toggle popover</button>";
						//span.innerHTML += " <button id='" + b_popup_id + "' type='button' class='btn btn-xs btn-danger' data-toggle='popover' container='body' title='Test' onmouseover=\"$('#" + b_popup_id + "').popover('show')\" data-content='And here's some amazing content. It's very engaging. Right?'>Click to toggle popover</button>";
					}
					else { span.innerHTML = panel_name; }
				}
				a.appendChild(span);
				li.appendChild(a);
				document.getElementById('dropdown-menu-panels').appendChild(li);
			}
			if (origin_project_phenotype == 'No Phenotype Associated') { }
			update_list_panels_name(origin_project_phenotype);
			if (dijit.byId('waiting')) {
				dijit.byId('waiting').hide();
			}
			document.getElementById("gene_name_editor").addEventListener('click', function(event) {
				event.stopPropagation();
			});
			launch_panels_from_db_done = 1;
		}
	}
	xhr.send(this_args);
	//https://www.polyweb.fr/cgi-bin/polymorphism-cgi//json_output_nodb/getPanelsFromDb.pl

}


function create_content_pane_polydiag(patient_name) {
	var t_polydiag;
	require(["dojo/_base/array","dijit/layout/TabContainer", "dijit/layout/ContentPane","dojo/domReady!"], function(array, TabContainer, ContentPane) {
		t_polydiag =  new ContentPane({
			style:"font-size:9px;max-height:800px;width:100%;overflow:auto scroll;",
			title: 'PolyDiag',
			id:"panel_polydiag_"+patient_name,
			jsId:"panel_polydiag_"+patient_name,
			preload:0,
			ioMethod:dojo.xhrPost,
			toto:1,
			patient:patient_name,
			onShow:  function() {
				update_last_interface_view_for_patient(patient_name, this.id);
				select_panel_filters_named('div_filters_polydiag');
				document.getElementById("tabFilters").style.height = '80px';
				window.onresize();
				this.toto ++;
				if(this.toto == 2) { dijit.byId('waiting').show(); }
				var obj = this;
				setTimeout(function () {
					if (obj.toto == 2 ){	
						obj.setArgsPost = load_polydiag_args();
						obj.setArgsPost.patients = patient_name;
						obj.ioArgs.content = obj.setArgsPost; 
						obj.setHref(url_report_patient_polydiag);
						dijit.byId('waiting').hide();
					}
				}, 1000);
			}
		});
	});
	return t_polydiag;

}

function create_pane_polysplice(patient_name) {
	var t_polysplice;
	require(["dojo/_base/array","dijit/layout/TabContainer", "dijit/layout/ContentPane","dojo/domReady!"], function(array, TabContainer, ContentPane) {
		t_polysplice =  new ContentPane({
			style:"font-size:9px;max-height:800px;width:100%;overflow:auto scroll;",
			title: 'PolySplice',
			id:"panel_polysplice_"+patient_name,
			jsId:"panel_polysplice_"+patient_name,
			preload:0,
			ioMethod:dojo.xhrPost,
			toto:1,
			patient:patient_name,
			onShow:  function() {
				update_last_interface_view_for_patient(patient_name, this.id);
				select_panel_filters_named('div_filters_polysplice');
				document.getElementById("tabFilters").style.height = '140px';
				window.onresize();
				this.toto ++;
				if(this.toto == 2) { dijit.byId('waiting').show(); }
				var obj = this;
				setTimeout(function () {
					if (obj.toto == 2 ){
		           		var lTmp = String(window.location.href.split('?')[0]).split('/');
		           		lTmp.pop();
		           		var url_polysplice = url_path + '/rnaseq/polysplice.pl';
		           		obj.setArgsPost = {}
						obj.setArgsPost.view_polyviewer=1;
						obj.setArgsPost.project=projectName;
						obj.setArgsPost.patient=patient_name;
						if (value_dejavu < 101) { obj.setArgsPost.dejavu = value_dejavu; }
						obj.setArgsPost.dejavu_percent = value_dejavu_percent_similar;
						obj.setArgsPost.min_score = value_score;
						if ($('#b_dejavu_min_ratio_10').prop('checked')) { obj.setArgsPost.only_dejavu_ratio_10 = 1; }
						var only_gene_positions = document.getElementById('input_gene_positions').value;
						if (only_gene_positions == '') {}
						else if (only_gene_positions.match(':')) { obj.setArgsPost.dejavu_percent = only_gene_positions; }
						else { obj.setArgsPost.dejavu_percent = only_gene_positions; }
						obj.ioArgs.content = obj.setArgsPost; 
						obj.setHref(url_polysplice);
						dijit.byId('waiting').hide();
					}
				}, 1000);
           }
		});
	});
	return t_polysplice;
}

function create_pane_polycyto(patient_name) {
	var t_polycyto;
	require(["dojo/_base/array","dijit/layout/TabContainer", "dijit/layout/ContentPane","dojo/domReady!"], function(array, TabContainer, ContentPane) {
		t_polycyto =  new ContentPane({
			style:"font-size:9px;max-height:1000px;width:100%;overflow:auto scroll;",
			title: 'PolyCyto',
			id:"panel_polycyto_"+patient_name,
			jsId:"panel_polycyto_"+patient_name,
			preload:0,
			ioMethod:dojo.xhrPost,
			toto:1,
			patient:patient_name,
			innerHTML:"<div id='div_polycyto_" + patient_name + "'></div>",
			onShow:  function() {
				update_last_interface_view_for_patient(patient_name, this.id);
				//document.getElementById("div_selector_interface").style.display = "";
				this.toto ++;
				if(this.toto == 2) { dijit.byId('waiting').show(); }
				var obj = this;
				setTimeout(function () {
					if (obj.toto == 2 ){
		           		var lTmp = String(window.location.href.split('?')[0]).split('/');
		           		lTmp.pop();
		           		var url_polycyto = lTmp.join('/') + '/html/manta/PolyCyto.html';
						url_polycyto += "?projectname=" + projectName + "&filename=" + patient_name + "&transloc=yes&iframe=1";
						load_polycyto_url(this.patient, "div_polycyto_"+patient_name, url_polycyto);
						dijit.byId('waiting').hide();
					}
				}, 1500);
           }
		});
	});
	return t_polycyto;
}

function create_pane_polydupdel(patient_name) {
	var t_polydupdel;
	require(["dojo/_base/array","dijit/layout/TabContainer", "dijit/layout/ContentPane","dojo/domReady!"], function(array, TabContainer, ContentPane) {
		t_polydupdel =  new ContentPane({
			style:"font-size:9px;max-height:1000px;width:100%;overflow:auto scroll;",
			title: 'PolyCyto',
			id:"panel_polydupdel_"+patient_name,
			jsId:"panel_polydupdel_"+patient_name,
			preload:0,
			ioMethod:dojo.xhrPost,
			toto:1,
			patient:patient_name,
			innerHTML:"<div id='div_polydupdel_" + patient_name + "'></div>",
			onShow:  function() {
				update_last_interface_view_for_patient(patient_name, this.id);
				select_panel_filters_named('div_filters_polydupdel');
				document.getElementById("tabFilters").style.height = '100px';
				window.onresize();
				this.toto ++;
				if(this.toto == 2) { dijit.byId('waiting').show(); }
				var obj = this;
				setTimeout(function () {
					if (obj.toto == 2 ){
						if (dijit.byId("sw_dude").value == "on"){
							dude_cgh(obj, project_name, patient_name);
						}
						else {
							obj.setArgsPost = load_polydude_args();
							obj.setArgsPost.patients = patient_name;
							obj.ioArgs.content = obj.setArgsPost; 
							obj.setHref(url_dude);
						}
						dijit.byId('waiting').hide();
					}
				}, 1500);
           }
		});
	});
	return t_polydupdel;
}



function launch_polyviewer(patient_name) {
	var url = url_polyviewer;
	var post_args = load_polyviewer_args();
    post_args["patients"] = patient_name;
	alert(url);
	alert(post_args);
}


//var var_load_polyviewer  = 0;
//function load_tab_polyviewer(){
//	 for (var i=0;i<tsamples.length;i++) {
//	    	 container_editor_polyviewer.addChild(tab_editor_polyviewer[tsamples[i].label]);
//			tab_editor_polyviewer[tsamples[i].label].add = 1;
//			
//	    }
//	if (var_load_polyviewer == 1) refresh_tab_polyviewer() ;
//}








function tabPatients(items , request ) {
    var tsamples = dijit.byId("gridPatients").selection.getSelected();
    var z=0;
    var argsPost = new Array();;

    tsamples.sort(function(a,b) {
        return b.child > a.child;
    });
    
    for (var i=0;i<tsamples.length;i++) { 
        var n = tsamples[i].label;
    	
     //   if  (tsamples[i].child == 0) continue;
       var ic ='<img src="images/polyicons//iconfinder_baby-boy_s.png" style="width:75%" >  &nbsp;';
     
       ic ='<img src="images/polyicons/iconfinder_baby-girl_16.png" style="width:75%" >  &nbsp;';
       
   		  if (tsamples[i].child == 0) {
   			 ic = '<img src="images/polyicons/icons8-person-24.png" style="width:75%" >&nbsp;';
   		  	 if (tsamples[i].sex == "F")   ic = '<img src="images/polyicons/icons8-person-female-24.png">&nbsp;';
   		  }
   		  else {
   			  
   			 if (tsamples[i].sex == "F" && tsamples[i].status == "1"){
   				ic ='<img src="images/polyicons//iconfinder_baby-girl_s.png" style="width:75%" >  &nbsp;';
   			 }
   			 if (tsamples[i].sex == "F" && tsamples[i].status == "2"){
    				ic ='<img src="images/polyicons//iconfinder_baby-girl_d.png" style="width:75%" >  &nbsp;';
   			 }
   			 if (tsamples[i].sex == "M" && tsamples[i].status == "1"){
    				ic ='<img src="images/polyicons//iconfinder_baby-boy_s.png" style="width:75%" >  &nbsp;';
    		 }
    		 if (tsamples[i].sex == "M" && tsamples[i].status == "2"){
     				ic ='<img src="images/polyicons//iconfinder_baby-boy_d.png" style="width:75%" >  &nbsp;';
    		}
   			  
   		  }
   		 //#FDAC53
   		  var star ="";
   		  //;
   		  if (tsamples[i].validation_status == -99 ){
   			  
   			  star = "<span class=' glyphicon glyphicon-exclamation-sign' aria-hidden='true' style='position:absolute;top:1px;left:1px;color:#BA0BE0'></span>";
   		  }
   		  else if (tsamples[i].validation_status >= 3 ) {
   			  star = "<span class='glyphicon glyphicon-star' aria-hidden='true' style='position:absolute;top:1px;left:1px;color:#E74C3C'></span>";
   		 }
   		 else if (tsamples[i].validation_status >= 1 ) {
  			  star = "<span class='glyphicon glyphicon-star' aria-hidden='true' style='position:absolute;top:1px;left:1px;color:black'></span>";
  		 }
   		 // else  if (tsamples[i].validation_status >= 3 ) {
   			//  star = "<span class='glyphicon glyphicon-star' aria-hidden='true' style='position:absolute;top:2px;right:3px;font-size:14px;color:black'></span><span class='glyphicon glyphicon-star' aria-hidden='true' style='position:absolute;top:1px;right:5px;font-size:14px;color:#FDAC53'></span>";
   		 // }
        if (!tab_editor_polyviewer[n]){
            z++;
     //       <td><span class='glyphicon glyphicon-star' aria-hidden='true' style='font-size:08px;></span></td>
            //<td><span class='glyphicon glyphicon-star' aria-hidden='true' ></span></td>
	var title1;
		if (tsamples[i].label == tsamples[i].fam){
			 title1 = star+"<table><tr><td rowspan=2>"+ic+ "</td><td style='font-size:10px;'>"+tsamples[i].label+"</td></tr></table>";
			
		}
		else {
			  title1 = star+"<table><tr><td rowspan=2>"+ic+ "</td><td style='font-size:09px;font-weight: bold;'>"+tsamples[i].label+"</td></tr><tr>"+"<td ><span style='font-size:09px;'>"+tsamples[i].fam+"</small></td></tr></table>";
		}
		
         tab_editor_polyviewer[tsamples[i].label] = new dojox.layout.ContentPane({ // new dijit.layout.ContentPane({
            title:title1,
            content: " ",
             closable: true,
             id:"c1pv"+n,
             ioMethod:dojo.xhrPost,
			 toto:1,
			patient:tsamples[i].label,
			
           onShow:  function() {
        	   		
				
				if (!(dijit.byId("check_never_pv"))){
					return;
				}
				this.toto ++;
				
                  if (this.toto == 2 ){	
					if (this.setArgsPost == null ){
						this.setArgsPost = load_polyviewer_args();
					}
                   	   this.setArgsPost.patients=this.patient;
                	   this.ioArgs.content = this.setArgsPost; 
                	   this.setHref(url_polyviewer);
                  }
                  
           }
        });
         tab_editor_polyviewer[tsamples[i].label].patient = tsamples[i].label;
         tab_editor_polyviewer[tsamples[i].label].toto = 0;
         
         tab_editor_polysplice[tsamples[i].label] = new dijit.layout.ContentPane({
            title:title1,
            content: " ",
             closable: true,
             id:"c1polysplice"+n,
             ioMethod:dojo.xhrPost,
			patient:tsamples[i].label,
			
           onShow:  function() {
           		var lTmp = String(window.location.href.split('?')[0]).split('/');
           		lTmp.pop();
           		var url_polysplice = url_path + '/rnaseq/polysplice.pl';
           		url_polysplice += "?view_polyviewer=1&project=" + projectName + "&patient=" + this.patient
           		url_last_launch_patient = url_polysplice;
				url_polysplice += "&dejavu=15&dejavu_percent=96&min_score=10&only_dejavu_ratio_10=1";
				url_polysplice += "&view_polyviewer=1";
				load_polysplice_url(this.patient, this.id, url_polysplice);
           }
        });

         tab_editor_polysplice[tsamples[i].label].patient = tsamples[i].label;
         tab_editor_polysplice[tsamples[i].label].toto = 0;
         
         
         tab_editor_polycyto[tsamples[i].label] = new dijit.layout.ContentPane({
            title:title1,
            content: " ",
             closable: true,
             id:"c1polycyto"+n,
             ioMethod:dojo.xhrPost,
			patient:tsamples[i].label,
			
           onShow:  function() {
           		var lTmp = String(window.location.href.split('?')[0]).split('/');
           		lTmp.pop();
           		var url_polycyto = lTmp.join('/') + '/html/manta/PolyCyto.html';
				url_polycyto += "?projectname=" + projectName + "&filename=" + this.patient + "&transloc=yes&iframe=1";
				load_polycyto_url(this.patient, this.id, url_polycyto);
           }
        });

         tab_editor_polycyto[tsamples[i].label].patient = tsamples[i].label;
         tab_editor_polycyto[tsamples[i].label].toto = 0;



         tab_editor_polydude[tsamples[i].label] = new dojox.layout.ContentPane({
            title:title1,
            content: " ",
             closable: true,
             id:"c5polydude"+n,
             ioMethod:dojo.xhrPost,
			patient:tsamples[i].label,
			 toto:1,
          			
           onShow:  function() {
				this.toto ++;
                  if (this.toto == 2 ){
					if (dijit.byId("sw_dude").value == "on"){
						dude_cgh(this,project_name,this.patient);
					}
					else {
						this.setArgsPost = load_polydude_args();
                   	   this.setArgsPost.patients=this.patient;
                	   this.ioArgs.content = this.setArgsPost; 
                	   this.setHref(url_dude);
					}
                  }
           }
        });

         tab_editor_polydude[tsamples[i].label].patient = tsamples[i].label;
         tab_editor_polydude[tsamples[i].label].toto = 0;

    }
	     if (!tab_editor[n]) {
                z++;
             tab_editor[tsamples[i].label] = new dojox.layout.ContentPane({ // new dijit.layout.ContentPane({
                title:star+"<table><tr><td rowspan=2>"+ic+ "</td><td style='font-size:08px;'>"+tsamples[i].label+"</td></tr><tr>"+"<td ><span style='font-size:08px;'>"+tsamples[i].fam+"</small></td></tr></table>",
                content: " ",
                 closable: true,
                 id:"c1pp"+n,
                 //ioMethod:dojo.xhrPost,
                 onClose: function() {
                 var t = this.get('title');
                 delete tab_patients[t];
                 return true;
              },
               onShow:  function() {
                       this.toto ++;
                       if (this.toto== 2){
						var u = polydiag_url(this.patient);
                         this.setHref(u);
                      }
               }
            });
            
		tab_editor[tsamples[i].label].add = 1;
		tab_editor[tsamples[i].label].patient = tsamples[i].label;
        tab_editor[tsamples[i].label].toto = 0;
        }


    } 
    dijit.byId('waiting').hide();
	

    }

var hash_polycyto_loaded = {};
function load_polycyto_url(patient_name, tab_id, url_polycyto) {
	console.log(tab_id);
	if (!(patient_name in hash_polycyto_loaded)) {
		var div = document.createElement('div');
		var html = "<iframe width='100%' height='1000px' style='' src='" + url_polycyto + "' title='PolyCyto Interface'></iframe>";
		div.innerHTML = html;
		document.getElementById(tab_id).appendChild(div);
		hash_polycyto_loaded[patient_name] = 1;
	}
}

var hash_polydude_loaded = {};

function load_polydude_url(patient_name, tab_id, url_polycyto) {
	if (!(patient_name in hash_polydude_loaded)) {
		var div = document.createElement('div');
		var html = "<iframe width='100%' height='1000px' style='' src='google.fr"+ "' title='PolyCyto Interface'></iframe>";
		div.innerHTML = html;
		document.getElementById(tab_id).appendChild(div);
		hash_polydude_loaded[patient_name] = 1;
	}
}

var timeout;	 
function timer_mouse_over_genes (html_values, b_genes_id, panel_name, only_found) {
	document.getElementById(b_genes_id).innerHTML = '<i>Don\'t click ! Please wait to export XLS</i>';
	$('#' + b_genes_id).animate({backgroundColor:'blue'}, 1500);
    timeout = setTimeout(function() {
    	dijit.byId('waiting').show();
    	export_panel_genes (panel_name);
    	return;
    	
    	
		dijit.byId("dialog_genes_in_panel").hide();
		$('#' + b_genes_id).css('background', '#BC243C');
		document.getElementById("content_genes_in_panel").innerHTML = '';
		var div = document.createElement('div');
		var html = "<div class='row no-gutters' style='height:800px;width:99%;padding-left:15px;'>";
		html += "<br>";
		html += "<center>";
		html += "<div class='btn-group btn-sm'>";
		html += "<button type='button' onclick=\"timer_mouse_over_genes('"+html_values+"', '"+b_genes_id+"', '"+panel_name+"');\" class='btn btn-warning btn-xs'><span style='font-size:8px;'>View ALL Genes</span></button>";
		html += "<button type='button' onclick=\"timer_mouse_over_genes('"+html_values+"', '"+b_genes_id+"', '"+panel_name+"', '1');\" class='btn btn-info btn-xs'><span style='font-size:8px;'>View FOUND Genes</span></button>";
		html += "<button type='button' onclick=\"export_panel_genes('"+panel_name+"');\" class='btn btn-success btn-xs'><span style='font-size:8px;'>Export XLS</span></button>";
		html += "</div>";
		html += "</center>";
		html += "<br>";
		html += "<br>";
		html += "<center>";
		var lGenes = String(html_values).split(',');
		for (var i = 0; i < lGenes.length; i++) {
			var gene_name = lGenes[i];
			var div_id = "div_" + gene_name;
			div_id = div_id.replace(' ', '');
			div_id = div_id.replace('.', '_');
			var myEle = document.getElementById(div_id);
			var bcolor_style = '';
			var bcolor_cnv_style = 'border:solid #607D8B 0.5px;background-color:white;color:black;';
			if(myEle) {
				bcolor_cnv_style = 'background-color:#607D8B;color:white;';
				var classes = myEle.className;
				var list_classes = classes.split(' ');
				var bcolor = list_classes[list_classes.length - 2];
				if (bcolor == 'grey') { bvolor = 'silver'; }
				bcolor_style = "background-color:" + bcolor + ";";
				var cnv_status = list_classes[list_classes.length - 1];
				if (cnv_status == 'cnv_del_medium') { bcolor_cnv_style = "background-color:#FF6961;color:white;"; }
				if (cnv_status == 'cnv_del_high') { bcolor_cnv_style = "background-color:#DD4124;color:white;"; }
				if (cnv_status == 'cnv_del_low') { bcolor_cnv_style = "background-color:#FFB2AE;color:white;"; }
				if (cnv_status == 'cnv_dup_medium') { bcolor_cnv_style = "background-color:#0080FF;color:white;"; }
				if (cnv_status == 'cnv_dup_high') { bcolor_cnv_style = "background-color:#0066FF;color:white;"; }
				if (cnv_status == 'cnv_dup_low') { bcolor_cnv_style = "background-color:#9DFFFF;color:white;"; }
		    }
		    else {
		    	if (only_found == '1') { continue; }
		    }
	        html += "<div class='col-sm-1' style='min-height:60px;border: solid black 1px;" + bcolor_style + "'>";
			html += "<br>";
			html += "<span>" + gene_name + "</span>";
			html += "<br>";
			html += "<div class='btn-group btn-sm'>";
			html += "<button type='button' onclick=\"zoomDude('" + gene_name + "','50');\" class='btn btn-secondary btn-xs' style='" + bcolor_cnv_style + "'><span style='font-size:6px;'><b>CNV</b></span></button>";
			var patientName = used_patient_in_dropdown_menu;
			var geneName = gene_name.replace(' ','');
			var cmd_var = "view_var_from_proj_gene_pat('" + projectName + "', '" + geneName + "', '" + patientName + "', '', 'all', 'nocnv')";
			html += "<button type='button' onclick=\"" + cmd_var + "\" class='btn btn-secondary btn-xs' style='background-color:#3AB795;color:white;'><span style='font-size:6px;'><b>VAR</b></span></button>";
			html += "</div>";
			html += "</div>";
		}
		html += "</div>";
		html += "</center>";
		div.innerHTML = html;
		document.getElementById('content_genes_in_panel').appendChild(div);
		dijit.byId("dialog_genes_in_panel").set('title', "Genes in Panel " + panel_name);
		dijit.byId("dialog_genes_in_panel").show();
	}, 1500);
}

function timer_mouse_exit_genes (b_genes_id, nb_genes) {
	document.getElementById(b_genes_id).innerHTML = nb_genes;
	$('#' + b_genes_id).stop();
	$('#' + b_genes_id).css('background', '#BC243C');
	clearTimeout(timeout);
}

function change_slider_ratio(val){
	var t = Math.floor(val);
	document.getElementById("htoto").value= t;
	//document.getElementById("htoto").value=dojo.number.format(val/100,{places:1,pattern:" #%"});
	}
var w =1400;
function hide_and_show (val) {
	
	var x = document.getElementById(val);
	var x2 =  document.getElementById("up_"+val);
	var x3 =  document.getElementById("right_"+val);
	var l = document.getElementById("image_"+val);
	/*
	 
*/
	if (x.style.display === "none") {
    x.style.display = "block";
/*	x2.style.display = "none";
	x3.style.display = "block";
	*/
	 l.classList.remove("glyphicon-triangle-bottom");
     l.classList.add("glyphicon-triangle-right");
	w +=270;
	table_float.style.width = w+"px";
	
  } else {
  x.style.display = "none";
	w -=270;
	table_float.style.width = w+"px";
	 l.classList.remove("glyphicon-triangle-right");
     l.classList.add("glyphicon-triangle-bottom");
/*  x2.style.display = "block";
 x3.style.display = "none";
*/
  }
	
}



function gene_dude() {
	var gene_name = document.getElementById('gene_name_dude').value;
	 var tab_selected_patient = used_patient_in_dropdown_menu;
	zoomDude(project_name, gene_name,tab_selected_patient, 'force_visualisation');
}
var parent = new Array();
var plot = new Array();
var cgh_div= new Array();
function dude_cgh (containerParent,project,patient) {
	container = containerParent.containerNode;
	container.innerHTML="";
	
		 cgh_div[patient] = dojo.create("div", { style:"display: flex; flex-direction: row;flex-flow: row wrap;"}  ,container);
		for(var i= 1 ;i<=24;i++){
			 var parent  =   dojo.create("div", { style:"height:280px;width:600px; border: 2px solid grey;margin:15px;padding:10px"} ,cgh_div[patient]);
			 dojo.create("div", { style:"margin:1px;padding:1px",innerHTML: '<img src="https://img.icons8.com/external-tal-revivo-green-tal-revivo/32/000000/external-planning-and-new-match-making-of-the-dna-on-a-clipboard-labs-green-tal-revivo.png"/>'+" Chromosome : "+i+"  Patient: "+patient} ,parent);
			 var plot = dojo.create("div", { style:"height:80%;width:99%; "} , parent);
			var url_plotChr = url_path + "/validation_variation/dude_json.pl";
			
			dojo.xhrGet({
   					 url:  url_plotChr + "?project=" + project + "&patient=" + patient + "&chr=" + i,
    				handleAs: "json",
					parent:parent,
					plot: plot,
    			load: function(obj) {
				PlotValuesChr_json(obj.items,this.parent,this.plot);
       					 /* here, obj will already be a JS object deserialized from the JSON response */
    			},
    			error: function(err) {
        /* this will execute if the response couldn't be converted to a JS object,
           or if the request was unsuccessful altogether. */
   				 }
				});
				//var  storeP = new dojo.data.ItemFileWriteStore({
					//	url:  url_plotChr + "?project=" + project + "&patient=" + patient + "&chr=" + i,
				//	});
				   //storeP.fetch({onComplete: PlotValuesChr_json});
		}
}

function PlotValuesChr_json(items,parent,plot) {
		
	var trio;

	var tabPosR;
	var tabRatio;
	var tabRX;
	var tabRY;


	var tabplotRX = [0, 0];
	var tabdelRX = [];
	var tabdupRX = [];
	var tabplotRY = [-1, 1];
	var tabdelRY = [];
	var tabdupRY = [];

		if (items.FY.NO_DATA ==1){
			return;
		}

			tabplotRY = items.FY.all;
			tabdupRY = items.FY.dup;
			tabdelRY = items.FY.del;
				
			var trace1 = {
				x: items.FX.all,
				y: tabplotRY,
				yaxis: 'y1',
				mode: 'markers',
				showlegend: false,
				name: ' ',
				marker: {
					color: "#E0B589",//"#88B04B",//"#F3E0BE",
					// E0B589,//'#AAAAAA',
					size: 2
				}
			};

			var trace2 = {
				x: items.FX.dup,
				y: tabdupRY,
				yaxis: 'y1',
				mode: 'markers',
				showlegend: false,
				name: ' ',
				marker: {
					color: "#0050EF",
					size: 4
				}
			};

			var trace3 = {
				x: items.FX.del,
				y: tabdelRY,
				yaxis: 'y1',
				mode: 'markers',
				showlegend: false,
				name: ' ',
				marker: {
					color: '#FF0097',
					size: 4
				}
			};

			var trace1F = {
				x: items.FX.all,
				y: items.FY.all_father,
				yaxis: 'y2',
				mode: 'markers',
				showlegend: false,
				name: ' ',
				marker: {
					color: '#DEEAEE',
					size: 2
				}
			};
			

			var trace2F = {
				x:  items.FX.dup,
				y: items.FY.dup_father,
				yaxis: 'y2',
				mode: 'markers',
				showlegend: false,
				name: ' ',
				marker: {
					color: "#1BA1E2",
					size: 3
				}
			};

			var trace3F = {
				x: items.FX.del,
				y: items.FY.del_father,
				yaxis: 'y2',
				mode: 'markers',
				showlegend: false,
				name: ' ',
				marker: {
					color: '#1BA1E2',
					size: 3
				}
			};
			var trace1M = {
				x: items.FX.all,
				y: items.FY.all_mother,
				yaxis: 'y3',
				mode: 'markers',
				showlegend: false,
				name: ' ',
				marker: {
					color: '#F5D6C6',
					size: 2
				}
			};
			

			var trace2M = {
				x:  items.FX.dup,
				y: items.FY.dup_mother,
				yaxis: 'y3',
				mode: 'markers',
				showlegend: false,
				name: ' ',
				marker: {
					color: "#D80073",
					size: 3
				}
			};

			var trace3M = {
				x: items.FX.del,
				y: items.FY.del_mother,
				yaxis: 'y3',
				mode: 'markers',
				showlegend: false,
				name: ' ',
				marker: {
					color: '#D80073',
					size: 3
				}
			};
			if (items.FY.all_mother.length == 0){
				//dojo.style(parent[items.CHR],"height","200px");
				dojo.style(parent,"height","200px");
				
				data = [trace1, trace2, trace3];
			layout = {
				yaxis1: { automargin: true, range: [-1.5, 1.5],domain:[0,1],title:{ text:patient ,font: {family: 'Verdana',size: 8,color: 'black'} } 
				 },
				margin: {l: 1,r: 1,b: 0,t: 0,pad: 0 },
				xaxis: {
					title: {
						text: '',font: {family: 'Courier New, monospace',size: 10,color: '#7f7f7f'}
					},
				},
			};
				
			}
			else {
			data = [trace1, trace2, trace3,trace1F,trace2F,trace3F,trace1M,trace2M,trace3M];
			layout = {
				yaxis1: { automargin: true, range: [-1.5, 1.5],domain:[0,0.33],title:{ text:patient ,font: {family: 'Verdana',size: 8,color: 'black'} } 
				 },
				margin: {l: 1,r: 1,b: 0,t: 0,pad: 0 },
				yaxis2: { domain:[0.33,0.66],range: [-1.5,1.5], title:{ text:'Father' ,font: {family: 'Courier New, monospace',size: 10,color: 'blue'}}
				},
				yaxis3: { domain:[0.66,1],range: [-1.5,1.5], title: {
																	text:'Mother' ,font: {family: 'Courier New, monospace',size: 10,color: 'pink'}}
						 },
				xaxis: {
					title: {
						text: '',font: {family: 'Courier New, monospace',size: 10,color: '#7f7f7f'}
					},
				},
			};
			}
			dojo.style(parent,"border-color","grey");
			var config = { modeBarButtonsToRemove: ['pan2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'toggleSpikelines', 'hoverClosestCartesian', 'hoverCompareCartesian'], displayModeBar: true };
			Plotly.newPlot(plot, data, layout, config);
			//Plotly.newPlot(plot_father[item.CHR], [trace1F,trace2F,trace3F], layout, config);
			//Plotly.newPlot(plot_mother[item.CHR], [trace1M,trace2M,trace3M], layout, config);
			
			if (items.FY.red > items.FY.blue){
				//dojo.style(parent[item.PATIENT][item.CHR],"border-color","red");
				dojo.style(parent,"border-color","red");
			}
			else if (items.FY.red < items.FY.blue ){
				//dojo.style(parent[item.PATIENT][items.CHR],"border-color","blue");
				dojo.style(parent,"border-color","blue");
			}
			
			plot.on('plotly_click', function(data) {
				var pts = '';

				for (var i = 0; i < data.points.length; i++) {

					pts = 'x = ' + data.points[i].x + '\ny = ' +

						data.points[i].y.toPrecision(4) + '\n\n';

				}

				alert('Closest point clicked:\n\n' + pts);

			});
			
			plot.on('plotly_selected', function(eventData) {
			if(eventData == null ){
				return;
			}
			dijit.byId("popup_dude").show();
			var selected_patient = used_patient_in_dropdown_menu;
		var argsPost = new Array();
    	argsPost["project"] = project_name;
 		argsPost["patients"] = selected_patient;
		argsPost["chr"] = items.CHR;
		argsPost["start"] = eventData.points[0].x;
		argsPost["end"] = eventData.points[eventData.points.length-1].x;
	
		post_show_popup_dude(argsPost,url_cgh_dude);
			
		//	var url="http://darwin.bipd.fr/cgi-bin/pnitschk/polymorphism-cgi//validation_variation/dude_editor_by_position.pl?project="+project_name+"&patients="+selected_patient+"&chr="+items.CHR+"&start="+eventData.points[0].x+"&end="+eventData.points[eventData.points.length-1].x;			
		//	window.open(url, '_blank');


		});
			
			
		}
		

