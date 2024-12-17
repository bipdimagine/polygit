var layoutProject = [
	{ field: "name", name: "Name", width: "10%", styles: 'text-align:center;' },
	{ field: "description", name: "Description", width: "35%", styles: 'text-align:center;' },
	{ field: "genome", name: "Genome", width: "5%", styles: 'text-align:center;' },
	{ field: "capture_name", name: "Capture", width: "15%", styles: 'text-align:center;' },
	{ field: "patient", name: "Patient", width: "10%", styles: 'text-align:center;' },
	{ field: "users", name: "Users", width: "20%", styles: 'text-align:center;' },
	//{ field: "date_last_analyse", name: "Analyse Date", width: "10%", styles: 'text-align:center;' },
	{ field: "new_pathogenic", formatter: formaterBadges, name: "PolyBTF", width: "120px", styles: 'text-align:center;' },
];


var layoutProjectRNA = [
	{
		cells: [
			{ field: "name", name: "Name", width: "10%", styles: 'text-align:center;'},
			{ field: "description", name: "Description", width: "30%", styles: 'text-align:center;' },
			{ field: "samples", name: "Samples",formatter: formaterSample, width: "7%", styles: 'text-align:center;' },
			{ field: "capture_name", name: "Capture", width: "15%", styles: 'text-align:center;' },
			//{ field: "phenotype", name: "Phenotype", width: "15%", styles: 'text-align:center;' },
			{ field: "genome", name: "Genome", width: "10%", styles: 'text-align:center;' },
			//{ field: "version", name: "annotation",formatter: formaterVersion, width: "10%", styles: 'text-align:center;'},
			{ field: "creation_date", name: "Creation date", width: "10%", styles: 'text-align:center;'},
			//{ field: "annotation_date", name: "annotation date", width: "10%", styles: 'text-align:center;' },
			//{ field: "latest_view", name: "latest view ", width: "10%", styles: 'text-align:center;' },
			//{ field: "latest_view", name: "latest view", width: "100px",  formatter: formaterGenome,styles: 'text-align:center;'},
			//{ field: "new_pathogenic", formatter: formaterBadgesV, name: "PolyBTF", width: "60px", styles: 'text-align:center;'},
			//{ field: "acmg", name: "Latest Validation", width: "20%", styles: 'text-align:center;', formatter:formaterValidation},
			{ field: "button_splices", name: "PolySplice ", width: "8%", styles: 'text-align:center;',formatter:formaterButtonRNA },
			{ field: "button_rna", name: "PolyRnaView ", width: "8%", styles: 'text-align:center;',formatter:formaterButtonRNA },
			{ field: "link_analyse", name: "Analyse Link ", width: "10%", styles: 'text-align:center;',formatter:formaterButtonRNA },
		],
	}
	];

var layoutProjectPolyviewer = [
	{
		cells: [
			{ field: "name", name: "Name", width: "15%", styles: 'text-align:center;'},
			{ field: "description", name: "Description", width: "35%", styles: 'text-align:center;' },
			{ field: "samples", name: "samples",formatter: formaterSample, width: "7%", styles: 'text-align:center;' },
			{ field: "capture_name", name: "Capture", width: "15%", styles: 'text-align:center;' },
			{ field: "phenotype", name: "phenotype", width: "15%", styles: 'text-align:center;' },
			{ field: "genome", name: "genome", width: "10%", styles: 'text-align:center;' },
			{ field: "version", name: "annotation",formatter: formaterVersion, width: "10%", styles: 'text-align:center;'},
			{ field: "creation_date", name: "creation date", width: "10%", styles: 'text-align:center;'},
			{ field: "annotation_date", name: "annotation date", width: "10%", styles: 'text-align:center;' },
			{ field: "latest_view", name: "latest view ", width: "10%", styles: 'text-align:center;' },
			//{ field: "latest_view", name: "latest view", width: "100px",  formatter: formaterGenome,styles: 'text-align:center;'},
			{ field: "new_pathogenic", formatter: formaterBadgesV, name: "PolyBTF", width: "60px", styles: 'text-align:center;'},
			{ field: "acmg", name: "Latest Validation", width: "20%", styles: 'text-align:center;', formatter:formaterValidation},
			{ field: "button", name: "Viewer ", width: "10%", styles: 'text-align:center;',formatter:formaterButton },
		],
	}
	];


var layoutProjectDefidiag = [
	{
		cells: [
			{ field: "name", name: "Name", width: "15%", styles: 'text-align:center;'},
			{ field: "description", name: "Description", width: "35%", styles: 'text-align:center;' },
			{ field: "samples", name: "samples",formatter: formaterSample, width: "7%", styles: 'text-align:center;' },
			{ field: "capture_name", name: "Capture", width: "15%", styles: 'text-align:center;' },
			{ field: "phenotype", name: "phenotype", width: "15%", styles: 'text-align:center;' },
			{ field: "genome", name: "genome", width: "10%", styles: 'text-align:center;' },
			{ field: "version", name: "annotation",formatter: formaterVersion, width: "10%", styles: 'text-align:center;'},
			{ field: "creation_date", name: "creation date", width: "10%", styles: 'text-align:center;'},
			{ field: "annotation_date", name: "annotation date", width: "10%", styles: 'text-align:center;' },
			{ field: "latest_view", name: "latest view ", width: "10%", styles: 'text-align:center;' },
			//{ field: "latest_view", name: "latest view", width: "100px",  formatter: formaterGenome,styles: 'text-align:center;'},
			{ field: "new_pathogenic", formatter: formaterBadgesV, name: "PolyBTF", width: "60px", styles: 'text-align:center;'},
			{ field: "acmg", name: "Latest Validation", width: "20%", styles: 'text-align:center;', formatter:formaterValidation},
			{ field: "button", name: "Viewer", width: "10%", styles: 'text-align:center;',formatter:formaterButton },
			{ field: "sites", name: "sites", width: "10%", styles: 'text-align:center;'},
		],
	}
	];


/*
	{
		cells:
	[
	//{ field: "users", name: "Users", width: "20%", styles: 'text-align:center;' },
	{ field: "acmg", name: "Latest Validation", width: "20%", styles: 'text-align:center;', formatter:formaterValidation,colSpan:2},
	{ field: "acmg", name: "Latest Validation", width: "20%", styles: 'text-align:center;', formatter:formaterValidation},
	{ field: "acmg", name: "Latest Validation", width: "20%", styles: 'text-align:center;', formatter:formaterValidation},
	{ field: "acmg", name: "Latest Validation", width: "20%", styles: 'text-align:center;', formatter:formaterValidation},
	{ field: "acmg", name: "Latest Validation", width: "20%", styles: 'text-align:center;', formatter:formaterValidation},
	],[
		{ field: "new_pathogenic", formatter: formaterBadges, name: "PolyBTF", width: "120px", styles: 'text-align:center;'},
		{ field: "new_pathogenic", formatter: formaterBadges, name: "PolyBTF", width: "120px", styles: 'text-align:center;'},
		]
	
	},
	/*{
		cells:
	[
	{ field: "new_pathogenic", formatter: formaterBadges, name: "PolyBTF", width: "120px", styles: 'text-align:center;'},
	{ field: "new_pathogenic", formatter: formaterBadges, name: "PolyBTF", width: "120px", styles: 'text-align:center;'},
	]
	}
	
];
*/

function view_rna_junctions(prj) {
	window.open("rna/rna_junctions.html?project="+prj, '_blank');
}

function view_rna_shiny(prj) {
	dijit.byId('waiting').show();
	var url = url_path + "/polyrnaseq/launch_shiny_polyrna.pl";
	$.getJSON( url, function( data ) {
		console.log(data);
		dijit.byId('waiting').hide();
		var url_launch = data.url_docker + '/app/' + prj;
		console.log(url_launch);
	    window.open(url_launch, '_blank');
    })
    return;
}

function view_polyviewer(this_value) {
	var list = this_value.split('::');
	var prj = list[1];
	var genome = list[2];
   var url;
   var regexp = /HG38/g;
   if (String(genome).match(regexp)) {
   		url = 'https://'+window.location.hostname + "/HG38/polyweb/polyviewer.html?project="+prj;
   }
   else {
		url = "polyviewer.html?project="+prj;
	}
	window.open(url,'_blank');
 }

function view_polyquery(this_value) {
	var list = this_value.split('::');
	var prj = list[1];
	var genome = list[2];
   var url;
   var regexp = /HG38/g;
   if (String(genome).match(regexp)) {
   		url = 'https://'+window.location.hostname + "/HG38/polyweb/vector/gene.html?project="+prj;
   }
   else {
		url = "vector/gene.html?project="+prj;
	}
	window.open(url,'_blank');
 }
 
function formaterButtonRNA(this_value) {
	var array_button = [];
	var super_list = this_value.split(';');
	for (i in super_list) {
		var list = super_list[i].split('::');
		var tt = " ";
		var dd = " disabled ";
		var style =" style='margin-top: 3px ;opacity: 0.25;' ";
		// 1 for junctions analyse
		if (list[0] == 1) {
			tt = " onClick=view_rna_junctions('" + list[1]+"') ";
			style =" style='margin-top: 3px ' ";
			dd = ""
			if (list[1] == 'disabled') {
				tt = "disabled";
			}
			style =" style='margin-top: 3px ' ";
			dd = ""
			var button = "<center><table><tr><td><div class=' btn btn-info btn-xs "+dd+"   ' aria-label='Left Align' "+tt+" >PolySplice</div></td></tr></table></center>";
			array_button.push(button);
		} 
		if (list[0] == 2) {
			tt = " onClick=view_rna_shiny('" + list[1]+"') ";
			style =" style='margin-top: 3px ' ";
			dd = ""
			var button = "<center><table><tr><td><div style='background-color:#ca4ff8;color:white;' class=' btn btn-info btn-xs "+dd+"   ' aria-label='Left Align' "+tt+" >PolyRnaSeq</div></td></tr></table></center>";
			array_button.push(button);
		} 
		if (list[0] == 3) {
			var button = "<center><button class='btn btn-success btn-xs' style='padding:3px;' onclick=\" window.open('" + list[2] + "','_blank')\"><span style='font-size:12px;'>" + list[1] + "</span></button><center>";
			array_button.push(button);
		} 
	}
	var html = array_button.join('  ');
	return html;
}

function formaterButton(this_value) {
	var list = this_value.split('::');
	var tt = " ";
	var dd = " disabled ";
	var style =" style='margin-top: 3px ;opacity: 0.25;' ";
	var tt1 = "onClick=view_polyquery('" + this_value + "') ";
	if (list[0] == 1) {
		tt = " onClick=view_polyviewer('" + this_value + "') ";
		style =" style='margin-top: 3px ' ";
		dd = ""
	} 
	
	
	
	return "<table><tr><td><div class=' btn btn-success btn-xs  ' aria-label='Left Align' "+tt1+" >PolyQuery&nbsp;</div></td></tr><tr><td><div class=' btn btn-primary btn-xs "+dd+"  '  aria-label='Left Align' "+tt+style+" >PolyViewer</div></td></tr></table>";
}

var layoutProjectDiag = [
	{ name: "Row", et: getRow, width: "5%", styles: 'text-align:center;' },
	{ field: "name", name: "Name", width: "10%", styles: 'text-align:center;' },
	{ field: "description", name: "Description", width: "35%", styles: 'text-align:center;' },
 	{ field: "genome", name: "Genome", width: "10%", styles: 'text-align:center;' },
	{ field: "capture_name", name: "Capture", width: "15%", styles: 'text-align:center;' },
 	{ field: "capture_type", name: "Capture type", width: "10%", styles: 'text-align:center;' },
	{ field: "nb_runs", name: "Nb runs", width: "10%", styles: 'text-align:center;' },
 	{ field: "nb_transcripts", name: "Genes", width: "10%", styles: 'text-align:center;' },
	{ field: "patient", name: "Patient", width: "10%", styles: 'text-align:center;' },

	{ field: "users", name: "Users", width: "20%", styles: 'text-align:center;' },
	{ field: "machines", name: "Sequencers", width: "10%", styles: 'text-align:center;' },
	{ field: "new_pathogenic", formatter: formaterBadges, name: "PolyBTF", width: "90px", styles: 'text-align:center;' },
];

	function formaterPolyviewer(this_value) {
	if (this_value == 1) {
		return '<img src="https://img.icons8.com/fluent/24/000000/ok.png"/>' ;
	}
	return '<img src="https://img.icons8.com/color/24/000000/box-important--v1.png"/>';
}
function formaterValidation(this_value) {
	if (this_value == "*") {
		return "-" ;
	}
	var words = this_value.split('_');
	var style = "style='background-color:grey'";
	if (words[3] == 5) {
		style = "";
	}
	
	if (words[3] == 4) {
		style = "style='background-color:rgb(255, 202, 58);color:black'";
	}
	
	if (words[3] == 5) {
		style = "style='background-color:beige;color:black'";
	}
	return "<div class=' btn btn-danger btn-xs  ' aria-label='Left Align' "+style+"><small>"+words[0]+" <hr style='margin-bottom: 0px;margin-top: 1px;'>"+words[1]+" <hr style='margin-bottom: 0px;margin-top: 1px;'>"+words[2]+"</div> ";
			
	return "-";
	
}
function formaterSince(this_value) {
	return this_value+" days ago";
	if(this_value > 0){
		return "<img src='images/polyicons/icons8-family-24.png'></img>" 
	}
	
	return "<img src='images/polyicons/iconfinder_baby-boy_d.png'></img>";
}

function formaterSample(this_value) {
	if(this_value == 1){
		return "<div class=' btn btn-primary  ' aria-label='Left Align' style='background-color:rgb(255, 202, 58);border-color:black;color:black;font-size:9px'; >SOLO</div> ";

	} 
	if(this_value == 3){
		return "<div class=' btn btn-primary  ' aria-label='Left Align' style='background-color:rgb(255, 202, 58);border-color:black;color:black;font-size:9px'; >TRIO</div> ";
	} 
	return this_value;

}

function getColorFromVersion(this_value) {
	var list = this_value.split('::');
	 var c1;
	 if(list[0] ==0){
		c1 = "rgb(39, 174, 96)";
	 }
	 else if(list[0] < 2){
		 c1 = "rgb(22, 160, 133)";
	 }
	 else if(list[0] < 4){
		 c1 = "rgb(230, 126, 34)";
	 }
	 else {
		 c1 = "rgb(211, 84, 0)";
	 }
	return c1;
}

function formaterVersion(this_value) {
	var list = this_value.split('::');
	var c1 = getColorFromVersion(this_value);
	var cmd = 'onClick=\'showVersion("' + list[1] + '", "' + c1 + '");\'';
	return '<div class="btn btn-primary" aria-label="Left Align" style="background-color:' + c1 + ';border-color:black;color:white;font-size:9px" ' + cmd + '>' + list[1] + '</div>';
}

function formaterGenome(this_value) {
	var list = this_value.split('::');
	return list[1];
	return "<span class='badge badge-pill badge-primary'  style='background-color: rgb(107, 91, 149);color:white;border-color:black;font-size:11px'>"+list[1]+"</span>";
	return "<div class=' btn btn-danger  ' aria-label='Left Align'  style='background-color: rgb(16, 121, 178);color:white;border-color:black;font-size:9px'>"+list[1]+" </div> ";

	if(list[0] > 0){
		return "<div class=' btn btn-primary  ' aria-label='Left Align' style='background-color: rgb(16, 121, 178);border-color:black;color:black;font-size:9px'"+cmd+"; >"+list[1]+"</div> ";
		}
		return "<div class=' btn btn-danger  ' aria-label='Left Align'  style='background-color: rgb(16, 121, 178);color:white;border-color:black;font-size:9px'>"+list[1]+" </div> ";
	/*
	if(list[0] > 0){
	return "<div class=' btn btn-primary  ' aria-label='Left Align' style='background-color:rgb(255, 202, 58);border-color:black;color:black;font-size:9px'; >"+list[1]+"</div> ";
	}
	return "<div class=' btn btn-danger  ' aria-label='Left Align'  style='background-color:rgb(107, 91, 149);color:white;border-color:black;font-size:9px'>"+list[1]+" </div> ";
	*/
}
function formaterGenome1(this_value) {
	var list = this_value.split('::');
	console.log(this_value);
	if (list[0] ==0){
	return "<span class='badge badge-light' style='background-color:rgb(255, 202, 58);color:white;padding:5px;font-size:9px;'>"+list[1]+"</span>"; 
	}
	else {
		return "<span class='badge badge-light' style='background-color:rgb(107, 91, 149);color:white;padding:5px;font-size:9px;'>"+list[1]+"</span>"; 

	}
}

function formaterBadgesV(this_value) {
	var html = ' ';
	var list = this_value.split(';');
	var nb_versions = list.length;
	//console.log (this_value);
	//html = "<i class='fa fa-flag' aria-hidden='true' style='color:grey;padding:2px;font-size:5px;text-shadow: 2px 4px 3px #333'></i>";

	for(i=0; i<nb_versions; i++) {
		var value = list[i];
		// onClick=\"view_new_hgmd_only_project('" + proj_name + "', '" + version + "');\"
		
		if (String(value) == ' ') { html += "<i class='fa fa-flag' aria-hidden='true' style='color:grey;padding:2px;font-size:5px;text-shadow: 2px 4px 3px #333'></i>"; }
		else if (String(value) == 'F') { html += "<i class='fa fa-flag' aria-hidden='true' style='color:grey;padding:2px;font-size:5px;text-shadow: 2px 4px 3px #333'></i>";  }
		else {
			var lTmp = String(value).split('::');
			var value = lTmp[0];
			var proj_name = lTmp[1];
			var version = lTmp[3];
			var tt = "onClick=view_new_hgmd_only_project('" + proj_name + "','"+ version + "') ";
			var text1 = "<div class=' bcircle ' aria-label='Left Align' "+tt;
			text1 ="";// "<span class='badge badge-light' ";
			
			var text2 =  "<i class='fa fa-flag' aria-hidden='true' style='padding:2px;font-size:20px;text-shadow: 2px 4px 2px #333'></i>";
			
			
			var style = "style='background-color:#EEE;color:rgb(200, 202, 202);'";
			if (value == 'A') {
				style ="style='background-color:#CE0000;border-color:black;color:white;'";
				html += "<i class='fa fa-flag' aria-hidden='true' style='color:red;padding:2px;font-size:25px;text-shadow: 2px 4px 3px #333'"+tt+"></i>";
				 //html += text1+style+text2;
			};
			if (value == 'B') {
				html += "<i class='fa fa-flag' aria-hidden='true' style='color:coral;padding:2px;font-size:25px;text-shadow: 2px 4px 3px #333'"+tt+"></i>";

				//style ="style='background-color:coral;border-color:black;color:white;'";
				 //html += text1+style+text2;
			};
			if (value == 'C') {
				html += "<i class='fa fa-flag' aria-hidden='true' style='color:#EFB73E;padding:2px;font-size:25px;text-shadow: 2px 4px 3px #333'"+tt+"></i>";

				//style = "style='background-color:#EFB73E;border-color:black;color:white;'";
				 //html += text1+style+text2;
			};
			
			if (value == 'D') {
				html += "<i class='fa fa-flag' aria-hidden='true' style='color:#DFCFBE;padding:2px;font-size:25px;text-shadow: 2px 4px 3px #333'"+tt+"></i>";

				//style = "style='background-color:#beige;border-color:black;color:black;'";
				 //html += text1+style+text2;
			};
			if (value == 'E') {
				html += "<i class='fa fa-flag' aria-hidden='true' style='color:#F5F5DC;padding:2px;font-size:25px;text-shadow: 2px 4px 3px #333'"+tt+"></i>";

				//style = "style='background-color:#beige;border-color:black;color:black;'";
				 //html += text1+style+text2;
			};
			
			/*if (value == 'A') { html += "<div class=' btn btn-primary  ' aria-label='Left Align' onClick=view_new_hgmd_only_project("+tt+") style='background-color:rgb(255, 202, 58);border-color:black;color:black;'; >Trio </div>"};
			//if (value == 'A') { html += " <button onClick=\"view_new_hgmd_only_project('" + proj_name + "', '" + version + "');\" style='border-radius:50px 20px;background-color:#CE0000;'><span style='color:white'><i>" + 'New+++' + "</i></span></button>"; }
			if (value == 'B') { html += " <button onClick=\"view_new_hgmd_only_project('" + proj_name + "', '" + version + "');\" style='border-radius:50px 20px;background-color:coral;'><span style='color:white'><i>" + 'New++' + "</i></span></button>"; }
			if (value == 'C') { html += " <button onClick=\"view_new_hgmd_only_project('" + proj_name + "', '" + version + "');\" style='border-radius:50px 20px;background-color:#EFB73E;'><span style='color:white'><i>" + 'New+' + "</i></span></button>"; }
			if (value == 'D') { html += " <button onClick=\"view_new_hgmd_only_project('" + proj_name + "', '" + version + "');\" style='border-radius:50px 20px;background-color:#FFFBAF;'><span style='color:black'><i>" + 'New+' + "</i></span></button>"; }
			if (value == 'E') { html += " <button onClick=\"view_new_hgmd_only_project('" + proj_name + "', '" + version + "');\" style='border-radius:50px 20px;background-color:white;'><span style='color:black'><i>" + 'New' + "</i></span></button>"; }
			*/
		}
	}
	if (html == ' '){
		html = "<i class='fa fa-minus' aria-hidden='true' style='color:#34495E;padding:2px;font-size:8px;text-shadow: 2px 3px 3px black'></i>";
	} 
	return html;
}

function formaterBadges(this_value) {
	var html = ' ';
	var list = this_value.split(';');
	var nb_versions = list.length;
	for(i=0; i<nb_versions; i++) {
		var value = list[i];
		if (String(value) == ' ') { html += ' '; }
		else if (String(value) == 'F') { html += ' '; }
		else {
			var lTmp = String(value).split('::');
			var value = lTmp[0];
			var proj_name = lTmp[1];
			var version = lTmp[3];
			if (value == 'A') { html += " <button onClick=\"view_new_hgmd_only_project('" + proj_name + "', '" + version + "');\" style='border-radius:50px 20px;background-color:#CE0000;'><span style='color:white'><i>" + 'New+++' + "</i></span></button>"; }
			if (value == 'B') { html += " <button onClick=\"view_new_hgmd_only_project('" + proj_name + "', '" + version + "');\" style='border-radius:50px 20px;background-color:coral;'><span style='color:white'><i>" + 'New++' + "</i></span></button>"; }
			if (value == 'C') { html += " <button onClick=\"view_new_hgmd_only_project('" + proj_name + "', '" + version + "');\" style='border-radius:50px 20px;background-color:#EFB73E;'><span style='color:white'><i>" + 'New+' + "</i></span></button>"; }
			if (value == 'D') { html += " <button onClick=\"view_new_hgmd_only_project('" + proj_name + "', '" + version + "');\" style='border-radius:50px 20px;background-color:#FFFBAF;'><span style='color:black'><i>" + 'New+' + "</i></span></button>"; }
			if (value == 'E') { html += " <button onClick=\"view_new_hgmd_only_project('" + proj_name + "', '" + version + "');\" style='border-radius:50px 20px;background-color:white;'><span style='color:black'><i>" + 'New' + "</i></span></button>"; }
		}
	}
	return html;
}
	


 function init_gridProject(url1){
        document.getElementById("title_polyweb").innerHTML = title_page;
              var username = dojo.cookie("username");
              var passwd = dojo.cookie("passwd");
                if (!username) {
                    dijit.byId('login').show();
                }
                else {
            	
				check_if_user_already_see_popup_new_pathogenic(dojo.cookie("username"));
                var url1 = url_listProject+"?type=login&pwd="+passwd+"&login="+username+"&action=list"+option_url;
                var jsonStore = new dojox.data.AndOrReadStore({ url: url1 });
                gridProject.setStore(jsonStore);
                var url_stat_polyweb1 =url_stat_polyweb+"?user="+username+"&pwd="+passwd;
                dojo.xhrGet({
                    url:url_stat_polyweb1, handleAs:"json",
                    load: function(data){
                     document.getElementById("span_projects").innerHTML = data.user.projects.all+"/"+data.global.projects.all;
                     document.getElementById("span_samples").innerHTML = data.user.patients.all+"/"+data.global.patients.all;
                     document.getElementById("span_exomes").innerHTML = data.user.patients.exome+"/"+data.global.patients.exome;
                     document.getElementById("span_genomes").innerHTML = data.user.patients.genome+"/"+data.global.patients.genome;
                     //document.getElementById("span_ciliomes").innerHTML = (data.user.patients.ciliome)+"/"+data.global.patients.ciliome;
                     document.getElementById("span_diagnostics").innerHTML = data.user.patients.target+"/"+data.global.patients.target;
        				check_new_hgmd_clinvar(1);

                    // data is a JavaScript object. The content of foo.php
                    // was passed through dojo.fromJson
             }
            });

            
            
                LoginOk = 1;
               
                }
                
            setTitle(username);
            }
            
           
            
            //executer fonction dblclickGrid si double click sur tableau des projets
            dojo.addOnLoad(function(){
                dojo.connect(gridProject, 'onRowDblClick', dblclickGrid);
                
            });
            
            
        
            
            //fonction permettant de se loger avec un autre login et mdp
             function relogin (e){
                dijit.byId('login').show();
            }
            
     
            
   
            
            //fonction à executer si double click sur une ligne
            dblclickGrid = function(e){
            	 var items = gridProject.selection.getSelected();
            	 
            	 //alert(items[0].name);
            	 
               // var prj = gridProject.getCell(1).getNode(e.rowNode.gridRowIndex).textContent;
                var project_type = "ng";//gridProject.getCell(2).getNode(e.rowNode.gridRowIndex).textContent;
                viewProject(items[0].name, project_type, items[0].genome);
            };
            
            //fonction appelée par dblclickGrid suite au dbl click sur une ligne du tableau
           
            
            //fonction pour donner le nom de l'utilisateur loger dur la page
            function setTitle (uname){
                
                lock_button.set("label", uname);
            }
        
       
            function checkPassword(dialogFields){
                var passwd = document.getElementById('passwd').value;
                var username = document.getElementById('username').value;
                var midnight = new Date();
                midnight.setHours(23,59,59,0);
				dojo.cookie("username", username, { expires:midnight ,path: "/" });
				dojo.cookie("passwd", passwd, { expires: midnight,path: "/" });
                init_gridProject();
                dijit.byId('login').hide();
            }
       
            
        