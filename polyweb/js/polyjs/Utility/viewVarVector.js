function getXlsUrl () {
    var xls_hurl = get_args_vector('other');
    xls_hurl += "&view_all_genes=1";
    xls_hurl += "&xls_by_variants=1";
    var used_only_genes = 0;
    if (args_pat_inattic_genes_filters) {
        if (args_pat_inattic_genes_filters == 'all') {}
        else {
            xls_hurl += "&filter_attic_genes="+args_pat_inattic_genes_filters;
            xls_hurl += "&only_genes="+args_genes_inter+"&filter_genes_intersect="+args_genes_inter;
            used_only_genes = 1;
        }
    }
    else if (lOnlyThesesGenes.length > 0) {
    	if (used_only_genes == 0) {
			xls_hurl += "&only_genes=" + lOnlyThesesGenes.join(',');
		}
	}
    xls_hurl = xls_hurl.replace(/&/g, ';');
    return xls_hurl;
}

function viewVarByGene (value) {
	var item = gridListGenes.getItem(value.rowIndex);
	var locus = item.chromosome+':'+item.start+'-'+item.end;
    var xls_hurl = getXlsUrl() + ";only_genes=" + item.xref;
    var lTmp = xls_hurl.split('?');
    var title = item.xref;
    if (filtering_mode == 'fam') { title += '_-_' + filtering_mode; }
    if (genetic_model == 'None') {}
    else { title = title +"_-_"+genetic_model; }
	var murl = "detailProject.html?project="+projectName+"&chromosome="+item.chromosome+"&gene="+item.name+"&genename="+item.xref+"&type="+project_type+"&vector_chr_size="+item.vector_chr_size+"&vector_ids="+item.vector_ids+"&type_cache="+item.type_cache+"&model="+genetic_model+"&famind="+famind+"&locus="+locus+"&xls_url="+lTmp[1]+"&title="+title;
	if (get_value_slider_regionho() > 1) { murl += "&region_ho="+get_value_slider_regionho(); }
	window.open(murl);
}
	
function goViewVar(items, request){
	var locus = items[0].chromosome+':'+items[0].start+'-'+items[0].end;
	var murl = "detailProject.html?project="+projectName+"&chromosome="+items[0].chromosome+"&gene="+items[0].name+"&genename="+items[0].xref+"&type="+project_type+"&vector_chr_size="+items[0].vector_chr_size+"&vector_ids="+items[0].vector_ids+"&type_cache="+items[0].type_cache+"&model="+genetic_model+"&famind="+famind+"&locus="+locus;
	window.open(murl);
	return;
}

function goViewAllVar() {
	var list = [];
    gridListGenes.store.fetch({
        onComplete: function(items){
            var max = items.length;
            var i = 1;
            var n = 0;
            while (n < max) {
                var chr = gridListGenes.store.getValue(items[n], 'chromosome');
                var vector_ids = gridListGenes.store.getValue(items[n], 'vector_ids');
                var to_push = chr + ':' + vector_ids;
                list.push(to_push);
                n++;
            }
        }
    });
    var xls_hurl = getXlsUrl();
    var lTmp = xls_hurl.split('?');
    var title = "ALL_Genes";
    if (filtering_mode == 'fam') { title += '_-_' + filtering_mode; }
    if (genetic_model == 'None') {}
    else { title = title +"_-_"+genetic_model; }
	var murl = "detailProject.html?project="+projectName+"&chromosome=all&gene=all&type_cache=all&genename=ALL&type="+project_type+"&vector_resume="+list.join('::')+"&model="+genetic_model+"&famind="+famind+"&xls_url="+lTmp[1]+"&title="+title;
	if (get_value_slider_regionho() > 1) { murl += "&region_ho="+get_value_slider_regionho(); }
	window.open(murl);
}

function goViewSelectedlVar() {
	if (document.getElementById("span_nb_sel_genes").innerHTML == 'View <b>ALL</b> Genes') { return goViewAllVar(); }
	var list = [];
	var list_only_genes = [];
    gridListGenes.store.fetch({
        onComplete: function (items, result) {
            var items = gridListGenes.selection.getSelected();
            if (items.length) {
                dojo.forEach(items, function(selectedItem) {
	                if (selectedItem !== null) {
		                var chr = gridListGenes.store.getValue(selectedItem, 'chromosome');
		                var vector_chr_size = gridListGenes.store.getValue(selectedItem, 'vector_chr_size');
		                var vector_ids = gridListGenes.store.getValue(selectedItem, 'vector_ids');
		                var to_push = chr + ':' + vector_chr_size + ':' + vector_ids;
	                	list.push(to_push);
		                var xref = gridListGenes.store.getValue(selectedItem, 'xref');
		                list_only_genes.push(xref);
	               	}
              	});
            }
    	}
	});
    var xls_hurl = getXlsUrl() + ";only_genes=" + list_only_genes.join(',');
    var lTmp = xls_hurl.split('?');
    var title = list_only_genes.length+"_Genes";
    if (filtering_mode == 'fam') { title += '_-_' + filtering_mode; }
    if (genetic_model == 'None') {}
    else { title = title +"_-_"+genetic_model; }
	var murl = "detailProject.html?project="+projectName+"&chromosome=all&gene=all&type_cache=all&genename=ALL&type="+project_type+"&vector_resume="+list.join('::')+"&model="+genetic_model+"&famind="+famind+"&xls_url="+lTmp[1]+"&title="+title;
	if (get_value_slider_regionho() > 1) { murl += "&region_ho="+get_value_slider_regionho(); }
	window.open(murl); 
}

var last_no_gene_selected = 0;
var union_or_intersect = 'first_launch';
viewPatientsImplicated = function(type) {
    isIntersectGenes = 0;
	var need_update_columns;
	if (last_no_gene_selected == 1) {
		need_update_columns = 1;
		last_no_gene_selected = 0;
	}
	if (type == 'no_action') {
		if (union_or_intersect == 'intersection' || union_or_intersect == 'union') { need_update_columns = 1; }
		union_or_intersect = 'no_action';
	}
	else if (type == 'intersection') {
		if (union_or_intersect == 'no_action' || union_or_intersect == 'first_launch') { need_update_columns = 1; }
		union_or_intersect = 'intersection';
	}
	else if (type == 'union') {
		union_or_intersect = 'union';
	}
	else if (union_or_intersect == 'first_launch') {
		return;
	}
    gridPatient.filter({name:'*'});
    gridFamilies.filter({name:'*'});
    var nb_genes = 0;
    var h_genes = {};
    var h_patients = {};
    var h_pat_details = {};
	var lSelectedGenes = gridListGenes.selection.getSelected();
    dojo.forEach(gridListGenes.selection.getSelected(), function (item,i) {
    	if (item == null) {}
    	else {
	    	var id = gridListGenes.store.getValue(item, "id");
	    	h_genes[id] = 1;
	    	var listPat = gridListGenes.store.getValue(item, "patients");
		    var lTmp = listPat.split(',');
		    dojo.forEach(lTmp, function(patient_id) {
		    	if (patient_id in h_patients) {
		        	h_patients[patient_id] += 1;
		        }
		        else {
		        	h_patients[patient_id] = 1;
		        }
		    });
		    nb_genes++;
	   }
    })
    if (nb_genes == 0) { return; }
    var nbInclude = 0;
    gridListGenes.store.fetch({
        onComplete: function(items, request) {
            var n = 0;
            var max = gridListGenes.rowCount;
            while (n < max) {
                var selectedItem = items[n];
                if (selectedItem) {
                    var id = gridListGenes.store.getValue(selectedItem, "id");
		            if (id in h_genes) {
		            	if (union_or_intersect == "intersection") {
		            		gridListGenes.store.setValue(selectedItem, 'include', 0);
		            	}
		            	else if (union_or_intersect == "union") {
		            		gridListGenes.store.setValue(selectedItem, 'include', 1);
		            	}
		        	}
                    var include = gridListGenes.store.getValue(selectedItem, "include");
                    if (include == 1) {}
                    else { nbInclude++; }
               }
               n++;
            }
        }
    });
    if (nbInclude == 0) {
    	isIntersectGenes = 0;
		gridListGenes.layout.setColumnVisibility(1, false);
		gridPatient.layout.setColumnVisibility(0, false);
		gridFamilies.layout.setColumnVisibility(0, false);
		document.getElementById("span_info_pat_commons").innerHTML = "";
		document.getElementById("span_info_fam_commons").innerHTML = "";
    }
    else {
    	isIntersectGenes = 1;
		gridListGenes.layout.setColumnVisibility(1, true);
		gridFamilies.layout.setColumnVisibility(0, true);
		var lCounts = getNbCommonPatientsAfterGeneIntersect();
		var nb_common_pat = lCounts[0];
		var nb_common_fam = lCounts[1];
		if (famind == 'Ind') {
			gridPatient.layout.setColumnVisibility(0, true);
        	dijit.byId("gridPatient").setSortIndex(0, false);
			document.getElementById("span_info_fam_commons").innerHTML = "";
			if (nb_common_pat == 0) { document.getElementById("span_info_pat_commons").innerHTML = "<b><i><font color='red'>"+nb_common_pat+" Patient</font></i></b>"; }
			else if (nb_common_pat == 1) { document.getElementById("span_info_pat_commons").innerHTML = "<b><i><font color='green'>"+nb_common_pat+" Patient</font></i></b>"; }
			else { document.getElementById("span_info_pat_commons").innerHTML = "<b><i><font color='green'>"+nb_common_pat+" Patients</font></i></b>"; }
		}
		else {
			gridPatient.layout.setColumnVisibility(0, false);
        	dijit.byId("gridFamilies").setSortIndex(0, false);
			document.getElementById("span_info_pat_commons").innerHTML = "";
			if (nb_common_fam == 0) { document.getElementById("span_info_fam_commons").innerHTML = "<b><i><font color='red'>"+nb_common_fam+" Family</font></i></b>"; }
			else if (nb_common_fam == 1) { document.getElementById("span_info_fam_commons").innerHTML = "<b><i><font color='green'>"+nb_common_fam+" Family</font></i></b>"; }
			else { document.getElementById("span_info_fam_commons").innerHTML = "<b><i><font color='green'>"+nb_common_fam+" Families</font></i></b>"; }
		}
    }
}

var isIntersectGenes;

function getNbCommonPatientsAfterGeneIntersect() {
	var nb_pat = 0;
	var nb_fam = 0;
    var lFam = [];
    var h_fam_found = {};
	if (gridListGenes.store) {
		var nbGenes_inter = 0;
	    var h_patfam_found = {};
	    gridListGenes.store.fetch({
	        onComplete: function(items, request) {
	            var n = 0;
	            var max = gridListGenes.rowCount;
	            while (n < max) {
	                var selectedItem = items[n];
	                if (selectedItem) {
	                    var id = gridListGenes.store.getValue(selectedItem, "id");
	                    var include = gridListGenes.store.getValue(selectedItem, "include");
	                    // Include == 0 pour Intersection
	                    if (include == 0) {
	                    	nbGenes_inter++;
	                    	var hToCheck = [];
	                    	var to_get;
							if (famind == 'Ind') { to_get = 'patients'; }
							else { to_get = 'families'; }
					    	var listPatFam = gridListGenes.store.getValue(selectedItem, to_get);
						    var lTmp = listPatFam.split(',');
						    dojo.forEach(lTmp, function(patient_text) {
						    	var lTmp2 = patient_text.split(';');
						    	var patient_id = lTmp2[0];
						    	if (nbGenes_inter == 1) { h_patfam_found[patient_id] = 1; }
						    	hToCheck[patient_id] = 1;
						    });
						    for (var patient_id in h_patfam_found) {
						    	if (patient_id in hToCheck) {}
						    	else { delete  h_patfam_found[patient_id]};
						    }
	                    }
	               }
	               n++;
	            }
	        }
	    });
	    
	    var lTOTO = [];
	    for (k in h_patfam_found) {
	    	lTOTO.push(k);
	    }
	    if (nbGenes_inter == 0) { return; }
	    gridPatient.store.fetch({
	        onComplete: function(items, request) {
	            var n = 0;
	            var max = gridPatient.rowCount;
	            while (n < max) {
	                var selectedItem = items[n];
                    var id = gridPatient.store.getValue(selectedItem, "id");
                    var name = gridPatient.store.getValue(selectedItem, "name");
                    var fam = gridPatient.store.getValue(selectedItem, "fam");
                	var to_get;
					if (famind == 'Ind') { to_get = id; }
					else if (famind == 'Fam') { to_get = fam; }
		            if (to_get in h_patfam_found) {
		            	nb_pat++;
		            	lFam.push(fam);
		            	gridPatient.store.setValue(selectedItem, "found", 'yes');
	            	}
	            	else { gridPatient.store.setValue(selectedItem, "found", ''); }
					n++;
	            }
	        }
	    });
	}
	dojo.forEach(lFam, function(fam) {
		if (fam in h_fam_found) {}
		else {
			nb_fam++;
			h_fam_found[fam] = 1;
		}
	});
    gridFamilies.store.fetch({
        onComplete: function(items, request) {
            var n = 0;
            var max = gridFamilies.rowCount;
            while (n < max) {
                var selectedItem = items[n];
                var name = gridFamilies.store.getValue(selectedItem, "name");
	            if (name in h_fam_found) { gridFamilies.store.setValue(selectedItem, "found", 'yes'); }
            	else { gridFamilies.store.setValue(selectedItem, "found", ''); }
				n++;
            }
        }
    });
	return [nb_pat, nb_fam];
}

function get_list_genes_intersected() {
	if (gridListGenes.store) {
		var lGenes = [];
	    gridListGenes.store.fetch({
	        onComplete: function(items, request) {
	            var n = 0;
	            var max = gridListGenes.rowCount;
	            while (n < max) {
	                var selectedItem = items[n];
	                if (selectedItem) {
	                    var xref = gridListGenes.store.getValue(selectedItem, "xref");
	                    var include = gridListGenes.store.getValue(selectedItem, "include");
	                    // Include == 0 pour Intersection
	                    if (include == 0) {
	                    	lGenes.push(xref);
	                    }
	               }
	               n++;
	            }
	        }
	    });
	    return lGenes.join(',');
   }
}

function check_pat_fam_to_intersect_from_genes_filters() {
	if (isIntersectGenes == null || isIntersectGenes == 0) { return; }
	if (gridListGenes.store) {
	    var nbPat_inter = 0;
	    var lPat_inTheAttic = [];
	    gridPatient.store.fetch({
	        onComplete: function(items, request) {
	            var n = 0;
	            var max = gridPatient.rowCount;
	            while (n < max) {
	                var selectedItem = items[n];
	                var id = gridPatient.store.getValue(selectedItem, "id");
	                var name = gridPatient.store.getValue(selectedItem, "name");
	            	var found = gridPatient.store.getValue(selectedItem, "found");
		            if (found == 'yes') { nbPat_inter++; }
	            	else { lPat_inTheAttic.push(name); }
					n++;
	            }
	        }
	    });
	    if (nbPat_inter == 0) { return 'all'; }
	    return lPat_inTheAttic.join('+');
	}
	return;
}

function display_b_viewSelectedGenesVar() {
	var nb = 0;
    dojo.forEach(gridListGenes._by_idx, function(item, index){
        var selectedItem = gridListGenes.getItem(index);
	    if (gridListGenes.selection.isSelected(index)) {
	    	nb++;
	    }
    });
    if (nb >= 1) {
    	document.getElementById("span_nb_sel_genes").innerHTML = "View <font color='green'><b>" + nb + "</b></font> Genes";
    	dijit.byId("b_viewSelectedGenesVar").set("disabled", false);
    }
    else {
    	document.getElementById("span_nb_sel_genes").innerHTML = 'View <b>ALL</b> Genes';
    	if (global_uniq <= 300) {
        	dijit.byId("b_viewSelectedGenesVar").set("disabled", false);
    	}
    	else {
        	dijit.byId("b_viewSelectedGenesVar").set("disabled", true);
    	}
    }
}