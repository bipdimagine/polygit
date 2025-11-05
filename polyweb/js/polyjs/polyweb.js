function view_new_hgmd() {

	if (first_project_new_pathogenic == '') {
		first_project_new_pathogenic = gridProject.getItem(0).name;
	}
	$( "#collapse_" + first_project_new_pathogenic ).addClass('in');

	dijit.byId('dialog_hgmd').show();
}

var is_new_hgmd_running = 0;
var is_copy_html_hgmd_done = 0;
function copy_new_hgmd_table() {
	if (is_new_hgmd_running == 1) {
        setTimeout(function() {
        	return open_new_hgmd_var();
        }, 1000);
	}
	else if (is_copy_html_hgmd_done == 1) {
		dijit.byId("dialog_hgmd").show();
		return;
	}
	else {
		is_copy_html_hgmd_done = 1;
		print_check_new_hgmd = 1;
		is_new_hgmd_running = 1;
		check_new_hgmd_clinvar();
		return open_new_hgmd_var();
	}
	return;
}

function view_new_hgmd_only_project(project_name, release) {
	open_new_hgmd_var();
	var tr_id = 'tr_' + project_name;
	var collapse_id = 'collapse_' + project_name;
	if (release == '') {}
	else {
		release = release.replace('.', '_');
		tr_id += '_' + release;
		collapse_id += '_' + release;
	}
	document.getElementById('content_res_only_project').innerHTML = "";
	$( "#" + tr_id ).clone().appendTo( "#content_res_only_project" );
	document.getElementById('content_res_only_project').innerHTML += "<br>";
	
	$( "#" + collapse_id ).addClass('in');
	$( "#" + collapse_id ).clone().appendTo( "#content_res_only_project" );
	dijit.byId('dialog_hgmd_only_project').show();
}
		
function select_tr(table_id, type_tr, value) {
//	alert('table_id: '+table_id);
	$("#"+table_id).find("tbody>tr").each(function() {
		var classList = this.className.split(/\s+/);
		if (classList && classList[0] === String(type_tr)){
//			alert(this.id+' - '+String(classList[1])+' - '+String(value));
			if(String(classList[1]) == String(value)) {
				this.style.display = "block";
			}
			else {
				this.style.display = "none";
			}
		}
	});
}

var view_others_projects_new_hgmd_clinvar = 0;
function check_new_hgmd_clinvar_others_projects() {
	view_others_projects_new_hgmd_clinvar = 1;
	check_new_hgmd_clinvar();
}

var date_last_days = '';
var date_now_hgmd = '';
var date_release_hgmd = '';
var resume_new_hgmd = '';
var print_check_new_hgmd = '';
var html_table_check_new_hgmd = '';	
var first_project_new_pathogenic;
var job_new_hgmd_clinvar = '';
var launch_new_hgmd = 0;
function check_new_hgmd_clinvar(first_elem, release) {
	if (launch_new_hgmd == 1) { return; }
	launch_new_hgmd = 1;
	var this_user_name = dojo.cookie("username");
	var url = url_path + "json_output_nodb/user_polybtf.pl?login=" + this_user_name + "&pwd=" + dojo.cookie("passwd") + "&count=1";
	dijit.byId('waiting').show();
    job_new_hgmd_clinvar = $.getJSON( url, function( data ) {
    	dijit.byId('waiting').hide();
    	$.each( data, function( key, val ) {
    		if (key == 'count') {
    			resume_new_hgmd = val;
    		}
    		if (key == 'news') {
    			document.getElementById("span_polybtf_wait").innerHTML = "";
    			document.getElementById("span_polybtf_news").innerHTML = "<font style='color:white;font-size:18px;'><i>PolyBTF <b>(" + val + ")</b></i></font>";
    		}
    		if (key == 'date_now') {
    			date_now_hgmd = val;
    		}
    		if (key =='date_last_days') {
    			date_last_days = val;
    		}
    		if (key == 'date_release') {
    			date_release_hgmd = val;
    		}
    		
    		<!--
    		if (key == 'nb_new_in_projects') {
    			if (resume_new_hgmd == '') {
	    			//document.getElementById("span_nb_new_hgmd").innerHTML = "<b>" + val + " New Pathogenic Var !</b>";
	    			document.getElementById("span_nb_new_hgmd").innerHTML = "<b>" + val + " New Var !</b>";
	    			document.getElementById("span_dialog_new_dm").innerHTML = "<i>New Pathogenic Variants found in <u>your project(s)</u>: <font color='red'><b>" + val + "</b></i></font>";
    				resume_new_hgmd = val;
	    		}
    		}
    		if (key == 'max_dejavu') {
    			document.getElementById("span_dialog_max_dejavu").innerHTML = "<span class='glyphicon glyphicon-share-alt' aria-hidden='true'></span> <u><i>variants with DejaVu > <b>" + val + "</b> projects are filtered...</u></i>";
    		}
    		if (key == 'max_gnomadac') {
    			document.getElementById("span_dialog_max_gnomad").innerHTML = "<span class='glyphicon glyphicon-share-alt' aria-hidden='true'></span> <u><i>variants with gnomAD AC > <b>" + val + "</b> are filtered...</u></i>";
    		}
    		if (key == 'html_table') {
    			document.getElementById("content_res").innerHTML = val;
    			update_date_user_connection_index(this_user_name);
    		}
    		if (key == 'current_version') {
    			document.getElementById("hgmd_current_version").innerHTML = "<b><i>" + val + "</b></i>";
				if(document.getElementById("hgmd_current_version_2")){
	    			document.getElementById("hgmd_current_version_2").innerHTML = "<b><i>" + val + "</b></i>";
	    		}	
    		}
    		if (key == 'first_project') {
    			first_project_new_pathogenic = val;
    		}
    		if (key == 'last_hgmd_release') {
    			document.getElementById("span_last_hgmd").innerHTML = val;
    		}
    		if (key == 'last_clinvar_release') {
    			document.getElementById("span_last_clinvar").innerHTML = val;
    		}
    		if (key == 'releases_available') {
				var list_release_dropdown_menu = [];
    			var div_dropdown_menu = document.getElementById('dropdown_releases_new_pathogenics');
    			div_dropdown_menu.innerHTML = "";
    			var list_releases = String(val).split(';');
    			for(i=0; i<list_releases.length; i++) {
    				var lTmp = list_releases[i].split('|');
    				var release_name = lTmp[0];
    				var release_info = lTmp[1];
    				list_release_dropdown_menu.push(release_name);
    				var b_check_id = 'dropdownCheck_' + release_name;
					var text = "<div class='form-check'>";
					text += "<input style='margin-right:8px;' type='checkbox' class='form-check-input' id='" + b_check_id + "'>";
					text += "<label class='form-check-label' for='dropdownCheck'>";
					text += release_info;
					text += "</label>";
					text += "</div>";
    				div_dropdown_menu.innerHTML += text;
				}
				div_dropdown_menu.innerHTML += "<div><button style='color:green;' onClick=\"check_new_hgmd_clinvar('1', '" + list_release_dropdown_menu.join(',') +"')\">Launch</button</div>";
			}
			if (key == 'releases_used') {
				document.getElementById("dropdown_releases_selected").innerHTML = val;
			}
			-->
		});
		if (value_dont_show_dialog_again_index == 1) {
			document.getElementById("img_dont_show_me_again").setAttribute("src", "images/polyicons/12-em-check.png");
		}
		else {
			document.getElementById("img_dont_show_me_again").setAttribute("src", "images/polyicons/cross.png");
		}
		is_new_hgmd_running = 0;
		launch_new_hgmd = 0;
		is_new_polybtf_for_username(this_user_name);
		if(document.getElementById("nb_new_proj_pathogenic_1")){
			count_nb_projects_new_pathogenics();
		}
		if (value_dont_show_dialog_again_index == 0) {
			open_new_hgmd_var();
		}
    });
}

function view_new_var_proj_gene_pat() {
	dijit.byId('dialog_var_from_proj_gene_pat').show();
}

var launch_new_var_proj_gene_pat = 0;
function view_var_from_proj_gene_pat(project_name, gene_tr_name, patient_name, keep_var_id, only_model, only_type) {
	if (launch_new_var_proj_gene_pat == 1) { return; }
	dijit.byId('waiting_variants').show();
	launch_new_var_proj_gene_pat = 1;
	var url = url_path + "json_output_nodb/view_all_var_from_project_patient_gene.pl";
	url += "?project=" + project_name;
	if (String(gene_tr_name).match(/ENST/)) {
		url += "&gene_transcript=" + gene_tr_name;
	}
	else {
		url += "&gene=" + gene_tr_name;
	}
	url += "&patient=" + patient_name;
	url += "&only_model=" + only_model;
	if (only_type) {  url += "&only_type=" + only_type;}
	url += "&keep=" + keep_var_id;
    $.getJSON( url, function( data ) {
    	dijit.byId('waiting_variants').hide();
    	$.each( data, function( key, val ) {
    		if (key == 'html_table') {
				document.getElementById("content_res").innerHTML = val;
				dijit.byId('dialog_hgmd').show();
    		}
		});
		//view_new_var_proj_gene_pat();
        setTimeout(function() {
			if (only_model == 'mother') { select_tr('table_variants', 'tr_variant', 'tr_mother') }
			if (only_model == 'father') { select_tr('table_variants', 'tr_variant', 'tr_father') }
        }, 1000);
		launch_new_var_proj_gene_pat = 0;
    });
}

var open_new_hgmd_var_done = 0;
function open_new_hgmd_var(btf_total_var) {
	if (open_new_hgmd_var_done == 1) {
		dijit.byId("dialog_hgmd").show();
		return;
	}
	var url = url_path + "json_output_nodb/user_polybtf.pl?login=" + dojo.cookie("username") + "&pwd=" + dojo.cookie("passwd");
	dijit.byId('waiting').show();
    job_new_hgmd_clinvar = $.getJSON( url, function( data ) {
    	dijit.byId('waiting').hide();
    	$.each( data, function( key, val ) {
    		if (key == 'html') {
    			document.getElementById("content_res").innerHTML = val;
    			enable_table_search();
    			document.getElementById("span_dialog_new_dm_global").innerHTML = "<i>New Pathogenic Variants found in DataBases and your projects: <b>" + btf_total_var + "</b></i>";
    			open_new_hgmd_var_done = 1;
				dijit.byId("dialog_hgmd").show();
    		}
    		
    		if (key == 'max_dejavu') {
    			document.getElementById("span_dialog_max_dejavu").innerHTML = "<span class='glyphicon glyphicon-share-alt' aria-hidden='true'></span> <u><i>variants with DejaVu > <b>" + val + "</b> samples are filtered...</u></i>";
    		}
    		if (key == 'max_gnomadac') {
    			document.getElementById("span_dialog_max_gnomad").innerHTML = "<span class='glyphicon glyphicon-share-alt' aria-hidden='true'></span> <u><i>variants with gnomAD AC > <b>" + val + "</b> are filtered...</u></i>";
    		}
    		if (key == 'current_version') {
    			document.getElementById("hgmd_current_version").innerHTML = "<b><i>" + val + "</b></i>";
				if(document.getElementById("hgmd_current_version_2")){
	    			document.getElementById("hgmd_current_version_2").innerHTML = "<b><i>" + val + "</b></i>";
	    		}	
    		}
    		if (key == 'last_hgmd_release') {
    			document.getElementById("span_last_hgmd").innerHTML = val;
    		}
    		if (key == 'last_clinvar_release') {
    			document.getElementById("span_last_clinvar").innerHTML = val;
    		}
			if (key == 'releases_used') {
				document.getElementById("dropdown_releases_selected").innerHTML = val;
			}
		});
    });
}

function enable_table_search() {
	$('#table_variants').bootstrapTable();
}

function zoomHgmdWithoutCss(p,hid,vid){
    dijit.byId("popup_hgmd").show();
    dijit.byId("div_popup_hgmd").setHref(url_hgmd_view+"?hid="+hid+"&vid="+vid+"&project="+p+'&polyquery_view=1');
}

function zoomHgmd(p,hid,vid){
    dijit.byId("popup_hgmd").show();
    dijit.byId("div_popup_hgmd").setHref(url_hgmd_view+"?hid="+hid+"&vid="+vid+"&project="+p);
}

function open_dejavu_infos(value,pname,is_in_this_run) {
	var this_url = url_path + "polydejavu/dejavu_var_infos.pl";
	this_url += "?login=" + dojo.cookie("username") + "&pwd=" + dojo.cookie("passwd") + "&input=" + value + "&export_html_bootstrap=1"+"&project="+project_name;
	if (project_name && pname) { this_url += "&patient="+pname; }
	if (is_in_this_run) { this_url += "&in_this_run=1"; }
	//alert (this_url);
	
	//var this_url = url_path + "json_output_nodb/captures_diag_from_transcripts.pl?tr_ids=" + values;
    dijit.byId('waiting').show();
	var xhrArgs = {
        handleAs: "json",
        url: this_url,
        load: function(data) {
			 dijit.byId('dialog_dejavu_var_infos').show();
			document.getElementById("content_dejavu_var_infos").innerHTML = data.html_table;
			dijit.byId('waiting').hide();
        },
        error: function(error) {
            alert("error");
        }
	};
	var deferred = dojo.xhrGet(xhrArgs);
}

function goto_dejavu_in_this_run(vid,pname){
	return open_dejavu_infos(vid,pname,project_name,1);
}

function goto_dejavu(vid,pname){
	return open_dejavu_infos(vid,pname);
	
	var url = url_dejavu+"?input="+vid;
	var win = window.open(url, '_blank');
	win.focus();
}

var old_byPatientDude_patient_name = '';
var old_byPatientDude_level = '';
function byPatientDude(patient,level){
	if (patient == old_byPatientDude_patient_name && level == old_byPatientDude_level) {
		dijit.byId("popup_dude").show();
		return;
	}
	old_byPatientDude_patient_name = patient;
	old_byPatientDude_level = level;
	var argsPost = new Array();
    argsPost["project"] = project_name;
 	argsPost["patients"] = patient;
	argsPost["quality"] = level;
	argsPost["genes"] = "";
	argsPost["force_visualisation"] = "";
	post_show_popup_dude(argsPost,url_cnv_overview_dude);
}
	
function zoomDude(project_name, gene,patient,force_visualisation){
	var argsPost = new Array();
    argsPost["project"] = project_name;
	argsPost["genes"] = gene;
	argsPost["patients"] = patient;
	argsPost["quality"] = "";
	argsPost["force_visualisation"] = "";
	if (force_visualisation) {
		argsPost["force_visualisation"] = 1;
	}
	post_show_popup_dude(argsPost,url_cnv_overview_dude);
}

function check_slider_new_pathogenic() {
	if (dijit.byId("slider_new_pathogenic").get("value") == 1) { onSearchAll('New+++','gridProject'); }
	if (dijit.byId("slider_new_pathogenic").get("value") == 2) { onSearchAll('New++','gridProject'); }
	if (dijit.byId("slider_new_pathogenic").get("value") == 3) { onSearchAll('New+','gridProject'); }
	if (dijit.byId("slider_new_pathogenic").get("value") == 4) { onSearchAll('New','gridProject'); }
	if (dijit.byId("slider_new_pathogenic").get("value") == 5) { onSearchAll('*','gridProject'); }
}

function count_nb_projects_new_pathogenics() {
	onSearchAll('New+++','gridProject');
	document.getElementById("nb_new_proj_pathogenic_1").innerHTML = gridProject.rowCount;
	onSearchAll('New++','gridProject');
	document.getElementById("nb_new_proj_pathogenic_2").innerHTML = gridProject.rowCount;
	onSearchAll('New+','gridProject');
	document.getElementById("nb_new_proj_pathogenic_3").innerHTML = gridProject.rowCount;
	onSearchAll('New','gridProject');
	document.getElementById("nb_new_proj_pathogenic_4").innerHTML = gridProject.rowCount;
	onSearchAll('*','gridProject');
	document.getElementById("nb_new_proj_pathogenic_all").innerHTML = gridProject.rowCount;
}

function clone_div_to_div(div1, div2) {
		setTimeout(function() {
        	//alert ('clone '+div1+' to '+div2);
			document.getElementById('content_detail_hgmg_ng_project').innerHTML = "";
		    $( "#" + div1 ).clone().appendTo( "#content_detail_hgmg_ng_project" );
		    //alert('cloned!');
		    dijit.byId('dialog_detail_hgmg_ng_project').show();
		    document.getElementById(div1).click;
        }, 1000);
	
}

function getRandomInt(max) {
  return Math.floor(Math.random() * Math.floor(max));
}

function is_new_polybtf_for_username(this_user_name) {
	var is_new_polybtf = 0;
	var url = url_path + "json_output_nodb/get_last_date_connection_for_user.pl?user_name=" + this_user_name + "&is_new_polybtf=1";
    job_new_hgmd_clinvar = $.getJSON( url, function( data ) {
    	$.each( data, function( key, val ) {
    		if (key == this_user_name) {
    			if (val == 1) {  is_new_polybtf = 1; }
    			else { is_new_polybtf = 0; }
    		}
    		else { is_new_polybtf = 0; }
		});
		
		
		var color_hgmd = "#222222";
		var l_colors = ['red', '#D35400', '#27AE60', '#669AE1', '#FFCA3A'];
		
		if (is_new_polybtf == 1) {
			color_hgmd = l_colors[getRandomInt(5)];
		}
		
//		if (is_new_polybtf == 1 && parseInt(date_last_days) <= 15) { color_hgmd = 'red'; }
//		else if (is_new_polybtf == 1 && parseInt(date_last_days) <= 20) { color_hgmd = 'orange'; }
//		else if (is_new_polybtf == 1 && parseInt(date_last_days) <= 30) { color_hgmd = 'green'; }
		document.getElementById("span_polybtf_wait").innerHTML = "";
		document.getElementById("span_polybtf_date_release").innerHTML = "<font style='color:white;font-size:18px;'>&nbsp;&nbsp;&nbsp;<span class='glyphicon glyphicon-arrow-right' font style='color:white;font-size:17px;' aria-hidden='true'></span>&nbsp;&nbsp;&nbsp;<b>Latest update " + date_release_hgmd + ", " + date_last_days + " days ago</b></font>";
		document.getElementById("span_polybtf_button").innerHTML = "&nbsp;&nbsp;&nbsp;<span class='glyphicon glyphicon-arrow-right' font style='color:white;font-size:17px;' aria-hidden='true'></span>&nbsp;&nbsp;&nbsp;<button style='font-size:12px;vertical-align:text-bottom;color:" + color_hgmd + ";background-color:white;border-radius:8px;' onclick='open_new_hgmd_var(resume_new_hgmd);'><u><i>View " + resume_new_hgmd + " (potential) New Variants</u></i></button>";
		
		
		
		
		if (color_hgmd) {
			document.getElementById("div_toolbar_polybtf").style.backgroundColor = color_hgmd;
			document.getElementById("span_polybtf_news").style.color = 'white';
		}
		
    });
}

function update_date_user_connection_index(this_user_name) {
	var url = url_path + "json_output_nodb/add_last_date_connection_for_user.pl?user_name=" + this_user_name;
    job_new_hgmd_clinvar = $.getJSON( url, function( data ) {
    	$.each( data, function( key, val ) {
    		if (key == 'ok') {
    			//check_if_user_already_see_popup_new_pathogenic(dojo.cookie("username"));
				//document.getElementById("img_dont_show_me_again").setAttribute("src", "images/polyicons/12-em-check.png");
				dijit.byId('waiting_general').hide();
    			//alert('OK, understood ! Next time, I will be invisible !');
    		}
		});
    });
}

var date_last_connection = '';
var value_dont_show_dialog_again_index = 0;
function dont_show_dialog_again_index() {
    dijit.byId('waiting_general').show();
	document.getElementById("img_dont_show_me_again").setAttribute("src", "images/polyicons/wait18trans.gif");
	if (value_dont_show_dialog_again_index == 0) {
		if (date_last_connection == '') {
			var url = url_path + "json_output_nodb/add_last_date_connection_for_user.pl?user_name=" + dojo.cookie("username");
		    job_new_hgmd_clinvar = $.getJSON( url, function( data ) {
		    	$.each( data, function( key, val ) {
		    		if (key == 'ok') {
		    			check_if_user_already_see_popup_new_pathogenic(dojo.cookie("username"));
						document.getElementById("img_dont_show_me_again").setAttribute("src", "images/polyicons/12-em-check.png");
    					dijit.byId('waiting_general').hide();
		    			alert('OK, understood ! Next time, I will be invisible !');
		    		}
				});
		    });
		}
		value_dont_show_dialog_again_index = 1;
	}
	else {
		if (date_last_connection == '') {}
		else {
			var url = url_path + "json_output_nodb/delete_last_date_connection_for_user.pl?user_name=" + dojo.cookie("username");
		    job_new_hgmd_clinvar = $.getJSON( url, function( data ) {
		    	$.each( data, function( key, val ) {
		    		if (key == 'ok') {
		    			check_if_user_already_see_popup_new_pathogenic(dojo.cookie("username"));
						document.getElementById("img_dont_show_me_again").setAttribute("src", "images/polyicons/cross.png");
    					dijit.byId('waiting_general').hide();
		    			alert('YES! With Pleasure !! Next time, I will appear ! ');
		    		}
				});
		    });
		}
		value_dont_show_dialog_again_index = 0;
	}
}

function check_if_user_already_see_popup_new_pathogenic(user_name) {
	var url = url_path + "json_output_nodb/get_last_date_connection_for_user.pl?user_name=" + user_name;
    job_new_hgmd_clinvar = $.getJSON( url, function( data ) {
    	$.each( data, function( key, val ) {
    		if (key == user_name) {
    			value_dont_show_dialog_again_index = 1;
    			date_last_connection = val;
    			return;
    		}
    		else {
    			value_dont_show_dialog_again_index = 0;
    			date_last_connection = '';
    			return;
    		}
		});
    });
}

function collapse(name){
	var d = document.getElementById(name);
	d.classList.add("in");
}

function preload_coverage(storeP,reload_flag,use_panel_name) {
	if (reload_flag ==1) { 
		reload(5, use_panel_name) ;
		return;
	}

	if (first >0) {
		         reload(5);
		         return;
		     }
		     
		     first ++;
		       var tsamples = storeP._arrayOfAllItems;//dijit.byId("gridPatients").selection.getSelected();
		    	
		    //var url_cnv =    construct_url(url_cnv_overview)+"&test_cache=1";
       // document.getElementById("waiting_text").textContent = "Loading Data... [0/25]";
	 require(["dojo/promise/all", "dojo/Deferred", "dojo/request", "dojo/_base/array", "dojo/dom-construct", "dojo/dom", "dojo/json", "dojo/domReady!"],
	function(all, Deferred, request, arrayUtil, domConstruct, dom, JSON) {
		//function preload_coverage(value){
            var p="";
            var samples = new Array();
            var requests = new Array();
             var url = url_flush+'?project='+project_name;
             var progress_id =0;
        	 var array_chr = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y'];
             var max =0;
              
                   var i,j,temparray,chunk = 1;
                   var labels = new Array;
                   var z =0;
                    for (i=0,j=array_chr.length; i<j; i+=chunk) {
                        temparray = array_chr.slice(i,i+chunk);
                        var temp = new Array();
                        for (var x=0;x<temparray.length;x++){
                            temp[x] = temparray[x].label;
                        }
                        
                        labels[z++] = array_chr[i];
                    }
                     
                          max =    labels.length +1;
                          
                          
                      for(i=0;i<Math.round(labels.length/2);i++){
                         requests[i] = request.get( url+"&patient_cached="+ labels[i], {
                                                        handleAs: "json"
                                                        }).then(function(response){
                                                            
                                                            progress_id ++;
                                                            myProgressBar.set( { value: progress_id } );
                                                          
                                                            return 1;
                                                            }) ;
                           // do whatever
                        }           
                    
                      
                    myProgressBar.set( { value: 0,maximum:max } );
                   dijit.byId('waiting_progress').show();   
                   

    
		var request2 = new Array();
		 all(requests).then(function(results){
		     
		      for(i=Math.round(labels.length/2);i<labels.length;i++){
                         request2[i+labels.length] = request.get( url+"&chromosome="+ labels[i], {
                                                        handleAs: "json"
                                                        }).then(function(response){
                                                            progress_id ++;
                                                            myProgressBar.set( { value: progress_id } );
                                                          
                                                            return 1;
                                                            }) ;
                           // do whatever
                        }           
                 request2[i+labels.length] = request.get(url_cnv, {
                                                        handleAs: "html"
                                                        }).then(function(response){
                                                            progress_id ++;
                                                            myProgressBar.set( { value: progress_id } );
                                                          
                                                            return 1;
                                                            }) ;
                                                             
                all(request2).then(function(results){
                    if (reload_flag ==1) reload(5, use_panel_name) ;
                     dijit.byId('waiting_progress').hide();
                });
                    
                 //  
                   //
                    
                
            });      

		});
		
		}
		

function post_show_popup_dude (setArgsPost,url) {
	dijit.byId("popup_dude").show();
		 var div_popup = 	dijit.byId("div_popup_dude");
		div_popup.ioMethod = dojo.xhrPost;
		div_popup.ioArgs.content = setArgsPost; 
		div_popup.setHref(url);
	
}
		


function reload_cnv(value,use_panel_name){
	//var start_url = url_coverage;
	var url = url_cnv_overview_dude ;
	//var url = construct_url(start_url);
	var widget = dijit.byId("slider_cnv");
	var value = widget.get("value");
	url += "?quality="+value+"&project="+project_name;//+"&patients=BOM_Mad";
	if (use_panel_name) {
	 	url += "&panel_name="+use_panel_name;
	}
	cnv_panel.setHref(url);
}
	

function reload(value) {
	      	var this_order;
			if (dijit.byId("ENST_order").get("value") == "on") { this_order = "name";}
    		else {this_order = "position";}
	      	
		 	var start_url = url_coverage;
		 	if (value == 5 ){
		 		start_url = url_coverage_overview;
		 		
		 		
		 		
		 	}
		 	if (value == 3){
                start_url = url_primer_coverage_overview;
                
            }
            if (value == 4){
            	start_url =url_cnv_overview ;
        		                 
            }
            
		 	var url = construct_url(start_url);
		 	
		 	var nbTranscripts = 1;
		 	
		 	if (value ==1){
		 		dijit.byId("ENST_order").set("value", "on");
				url+= "&only_red=1";
			}
			else if (value ==2){
		 		dijit.byId("ENST_order").set("value", "on");
				url+= "&only_orange=1";
			}
			else if (value ==3){
		 		dijit.byId("ENST_order").set("value", "on");
				url+= "&multiplex=1";
			}
			else if (value ==4){
				url+= "&cnv=1";
				url+='&order='+this_order;
			}
			else if (value ==5){
				url+='&order='+this_order;
		//		if (dijit.byId("only_low")){
				url+="&only_low="+dijit.byId("only_low").get('value');
			
				//url+="&recessive="+dijit.byId("recessive").get('value');
				//url+="&recessive="+dijit.byId("recessive").get('value');
				
			//	}
			}
		var step = 100;
		var first =-1;
		var previous = -1;
		var ii =1;
		for (var n=0;n<nbTranscripts;n+=step){ 
			if (first == -1){ first =ii;} 
		  if (!tab_coverage[ii]){
			 tab_coverage[ii] = new dijit.layout.ContentPane({
         		title: ii,
         		content: " ",
         		 closable: true,
         		 id:"cov50_"+n,
         		 onClose: function() {
         		 var t = this.get('title');
         		 delete tab_patients[t];
                 return true;
              },
               onShow:  function() {
               		var z =this.nextTab +0;
               		if (this.nextTab && tab_coverage[z]){
               			tab_coverage[z].refresh();
               		}
               		
               }
    		});
			
    	table_coverage.addChild(tab_coverage[ii]);
    	}
		
        tab_coverage[ii].setHref(url+"&start="+n+"&step="+step);
		 tab_coverage[ii].nextTab = ii+1; 
			ii++; 	
		}
		table_coverage.selectChild( tab_coverage[first]);
			select_exons = new Array();
			sanger_variations_valid = new Array();
			sanger_exons_valid = new Array();
			type = "todo";
			
		   }

function launch_web_igv_add_bam_and_locus(project, name, file, locus){
	LoadIGVPatient_editor(name, file);
	launch_web_igv_js(project, name, file, locus);
}
var previous_bams = new Array();
var igv_genome;

function launch_web_igv_js(project, patients_names, bams_files, locus){
	const regex = /HG38/g;
	var temp = bams_files.split('/');
	igv_genome = temp[3];
	var genome = 'hg19';
	if (igv_genome.match(regex)) {
		genome = 'hg38';
	}
	
	var locus2 = locus.replace(':', ';');
	locus2 = locus2.replace('-', ';');
	var lTmp = locus2.split(';');
	var chr = 'chr' + lTmp[0];
	var start = lTmp[1];
	var end = start;
	if (lTmp.length == 3) { end = lTmp[2]; }
	var list_bams_files = [];
	if (bams_files.match(/,/)) { list_bams_files = bams_files.split(','); }
	else if (bams_files.match(/;/)) { list_bams_files = bams_files.split(';'); }
	else { list_bams_files.push(bams_files); }
	view_web_igv_bam("dialog_igv", "div_igv", locus, bams_files, patients_names, genome);
	
	for (i in list_bams_files ) {
		var file = list_bams_files[i];
		if (!(file in previous_bams)){
			previous_bams[file] = 1;
			displayOneBAMIGV(file);
		}
		else {
			previous_bams[file] +=1;
		}
		//if (previous_bams[file] )
	}
	displayInIGVLocus(chr+":"+start+"-"+end);
}

function LoadIGVPatient_editor (patients, bams) {
}

var url_gene_bed, url_fasta, url_cyto;
function getGencodeTrackForProject(project_name) {
	var url = url_path + "igv/get_url_files_igv.pl?project=" + project_name;
	$.getJSON( url, function( data ) {
		url_gene_bed = data.url_gene_bed;
		url_fasta = data.url_genome_fasta;
    })
}

function displayListBamInIgvApp(list_bams) {
	for(var i = 0; i < list_bams.length; i++) {
		var bam = list_bams[i];
		displayOneBAMIGV(bam);
    }
}

function addListBamInIgvApp(list_bams) {
	for(var i = 0; i < list_bams.length; i++) {
		var bam = list_bams[i];
		addOneBAMIGV(bam);
    }
}

function click_annot(value) {

}
    		
function change_dialog_dimensions(dialog, padding) {
	if (padding) {}
	else { padding = 35; } 
	var obj_sizes = getViewPortSize();
	var div_width = obj_sizes['width'] - padding;
	var div_heigth = obj_sizes['heigth'] - padding;
	dialog.set("dimensions", [div_width, div_heigth]);
}

function getViewPortSize() {
	var viewPortWidth;
	var viewPortHeight;
	// mozilla/netscape/opera/IE7
	if (typeof window.innerWidth != 'undefined') {
		viewPortWidth = window.innerWidth;
		viewPortHeight = window.innerHeight;
	}
	// IE6 in standards compliant mode
	else if (typeof document.documentElement !== 'undefined' && typeof document.documentElement.clientWidth !== 'undefined' && document.documentElement.clientWidth !== 0) {
		viewPortWidth = document.documentElement.clientWidth;
		viewPortHeight = document.documentElement.clientHeight;
	}
	// older versions of IE fallback
	else {
		viewPortWidth = document.getElementsByTagName('body')[0].clientWidth;
		viewPortHeight = document.getElementsByTagName('body')[0].clientHeight;
	}
	return {
		width: viewPortWidth,
		heigth: viewPortHeight
	};
}
var igv_js = true;
function checkWindowIGV (check) {
	dojo.cookie("igv_js", check.get('value'), { expires:300 ,path: "/" });
	 igv_js =  check.checked;
//alert(check.get('value'));
}
