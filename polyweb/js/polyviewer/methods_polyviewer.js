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

var last_load_polyviewer_argsPost = {};

function load_polydude_args() {
	var argsPost = {project:project_name,
					quality:1
					};
	if ( dijit.byId("slider_dude_ac")){
		argsPost.quality = dijit.byId("slider_dude_ac").value
	}
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
 var  tsamples = storeP._arrayOfAllItems;
	 for (var i=0;i<tsamples.length;i++) {
	   var label = tsamples[i].label;
	   //tab_editor_polydude[label].setArgsPost = argsPost;
		
	   //tab_editor_polydude[label].transcripts = t;
	   tab_editor_polydude[tsamples[i].label].toto = 1;
   }
 	if (tab_editor_polydude){
     tab_selected_patient = return_selected_patient(tab_editor_polydude);
	   tab_editor_polydude[tab_selected_patient].toto = 1;

     tab_editor_polydude[tab_selected_patient].onShow();
   }
   return;
}

function cgh(){
	if (tab_editor_polydude){
     tab_selected_patient = return_selected_patient(tab_editor_polydude);
	   tab_editor_polydude[tab_selected_patient].toto = 1;
 		tab_editor_polydude[tab_selected_patient].cgh = 1;
     	tab_editor_polydude[tab_selected_patient].onShow();
   }
}


function load_polyviewer_export_xls(mode) {
	dijit.byId('waiting').show();
	var url;
    var argsPost = load_polyviewer_args(mode);
	var tab_selected_patient = return_selected_patient(tab_editor_polyviewer);
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
   for (var i=0;i<tsamples.length;i++) {
	   var label = tsamples[i].label;
	   tab_editor_polyviewer[label].setArgsPost = argsPost;
		
	   tab_editor_polyviewer[label].transcripts = t;
	   tab_editor_polyviewer[tsamples[i].label].toto = 1;
   }
 	if (tab_editor_polyviewer){
     tab_selected_patient = return_selected_patient(tab_editor_polyviewer);
     tab_editor_polyviewer[tab_selected_patient].onShow();
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

function force_refresh_tab_polycto_selected(name) {
	if (name in hash_polycyto_loaded) {
	
		//console.log( $('id_cnv_dup_del') );
		//console.log( $('SVCompare_gridDetails') );
		//document.getElementById('SVCompare_gridDetails').render();
		//$('SVCompare_gridDetails').render();
	
//		var div_dupdel = document.getElementById('id_cnv_dup_del');
//		console.log(div_dupdel);
//		if (div_dupdel) {}
//		else {
			//document.getElementById('iframe_polycyto_'+name).src += ''
//		}
	
//		console.log($('iframe_polycyto_'+name));
		
//		console.log($('iframe_polycyto_'+name).context);
		
//		$('iframe_polycyto_'+name).context.URL = $('iframe_polycyto_'+name).context.URL;
		
//		console.log($('c1polycyto'+name));
//		if (hash_polycyto_loaded[name] == 1) {
//			$('c1polycyto'+name).prevObject[0].links[1].click();
//			$('c1polycyto'+name).prevObject[0].links[0].click();
//		}


	}
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
				tab_editor_polyviewer[name].onShow();
				first ++;
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

function tabPatients(items , request ) {
    //var tsamples = storeP._arrayOfAllItems;//dijit.byId("gridPatients").selection.getSelected();
    var z=0;
    var argsPost = new Array();;

   tsamples = items;
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
		var title2;
		if (tsamples[i].label == tsamples[i].fam){
			title1 = star+"<table><tr><td rowspan=2>"+ic+ "</td><td style='font-size:10px;'>"+tsamples[i].label+"</td></tr></table>";
			title2 = star+"<table onclick=\"force_refresh_tab_polycto_selected('"+tsamples[i].label+"');\"><tr><td rowspan=2>"+ic+ "</td><td style='font-size:10px;'>"+tsamples[i].label+"</td></tr></table>";
		}
		else {
			title1 = star+"<table><tr><td rowspan=2>"+ic+ "</td><td style='font-size:09px;font-weight: bold;' >"+tsamples[i].label+"</td></tr><tr>"+"<td ><span style='font-size:09px;'>"+tsamples[i].fam+"</small></td></tr></table>";
			title2 = star+"<table onclick=\"force_refresh_tab_polycto_selected('"+tsamples[i].label+"');\"><tr><td rowspan=2>"+ic+ "</td><td style='font-size:09px;font-weight: bold;'>"+tsamples[i].label+"</td></tr><tr>"+"<td ><span style='font-size:09px;'>"+tsamples[i].fam+"</small></td></tr></table>";
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
         
         tab_editor_polycyto[tsamples[i].label] = new dijit.layout.ContentPane({
            title:title2,
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
           },
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
	if (!(patient_name in hash_polycyto_loaded)) {
		var div = document.createElement('div');
		var id_iframe = "iframe_polycyto_" + patient_name;
		var html = "<iframe id='" + id_iframe + "' width='100%' height='1000px' style='' src='" + url_polycyto + "' title='PolyCyto Interface'></iframe>";
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
			var patientName = return_selected_patient(tab_editor_polyviewer);
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
	 var tab_selected_patient = return_selected_patient(tab_editor_polydude);
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
			var selected_patient = return_selected_patient(tab_editor_polydude);
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
		

