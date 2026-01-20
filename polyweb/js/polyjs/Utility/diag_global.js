function check_slider_cat_annotations() {
	var value = parseInt(dijit.byId("slider_impact_pv").get('value'));
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
}




function check_slider_gnomad_ac() {
	return;
	if (dijit.byId("slider_gnomad_ac").get('value') > max_slider_gnomad_ac) {
		dijit.byId("slider_gnomad_ac").set('value', max_slider_gnomad_ac);
	}
}

var max_slider_allele_quality = 5;
var max_slider_gnomad_ac = 6;
var max_slider_dejavu = 100;
var max_slider_frequency_query = 5;

function update_buttons_polyquery() {

	if (document.URL.match('defidiag') !== null) {
		dijit.byId("slider_gnomad_ac_pv").set('value', 5);
		dijit.byId("slider_gnomad_ac_ho_pv").set('value', 2);
		//slider_gnomad_ac.value
	}
	if (is_genome_project ==0 ){
		dijit.byId("slider_dv_pv").set('value', 4);
		
	}
	if (is_genome_project == -1 ){
		dijit.byId("slider_dv_pv").set('value', 4);
		dijit.byId("slider_dv_ho_pv").set('value', 3);
		dijit.byId("slider_gnomad_ac_ho_pv").set('value', 3);
		
		
	}
	if (polyquery_filters) {
		dijit.byId('waiting').show();


		document.getElementById("both").checked = true;
		var lFreqFilters = [];
		text_dejavu.innerHTML = '';

		var lTmp1 = polyquery_filters.split('?');
		var lTmp2 = lTmp1[1].split('&');
		for (var i = 0; i < lTmp2.length; i++) {
			var lTmp3 = lTmp2[i].split('=');
			if (lTmp3[0] == 'filter_type_variation') {
				var lTmp4 = lTmp3[1].split('+');
				for (j = 0; j < lTmp4.length; j++) {
					if (lTmp4[j] == 'intergenic') {
						dijit.byId("filter_intergenic").set('checked', false);
						//dijit.byId("filter_intergenic").set('disabled', true);
					}
					if (lTmp4[j] == 'intronic') {
						dijit.byId("filter_intronic").set('checked', false);
						//dijit.byId("filter_intronic").set('disabled', true);
					}
					if (lTmp4[j] == 'upstream_downstream') {
						dijit.byId("filter_upstream_downstream").set('checked', false);
						//dijit.byId("filter_upstream_downstream").set('disabled', true);
					}
					if (lTmp4[j] == 'silent') {
						dijit.byId("filter_silent").set('checked', false);
						//dijit.byId("filter_silent").set('disabled', true);
					}
					if (lTmp4[j] == 'utr') {
						dijit.byId("filter_utr").set('checked', false);
						//dijit.byId("filter_utr").set('disabled', true);
					}
					if (lTmp4[j] == 'pseudogene') {
						dijit.byId("filter_pseudogene").set('checked', false);
						//dijit.byId("filter_pseudogene").set('disabled', true);
					}
					if (lTmp4[j] == 'nonframeshift') {
						dijit.byId("filter_nonframeshift").set('checked', false);
						//dijit.byId("filter_nonframeshift").set('disabled', true);
					}
					if (lTmp4[j] == 'non-frameshift') {
						dijit.byId("filter_nonframeshift").set('checked', false);
						//dijit.byId("filter_nonframeshift").set('disabled', true);
					}
					if (lTmp4[j] == 'nonsynonymous') {
						dijit.byId("filter_nonsynonymous").set('checked', false);
						//dijit.byId("filter_nonsynonymous").set('disabled', true);
					}
					if (lTmp4[j] == 'nonsynonymous') {
						dijit.byId("filter_nonsynonymous").set('checked', false);
						//dijit.byId("filter_nonsynonymous").set('disabled', true);
					}
					if (lTmp4[j] == 'splice_site') {
						dijit.byId("filter_splice_site").set('checked', false);
						//dijit.byId("filter_splice_site").set('disabled', true);
					}
					if (lTmp4[j] == 'splicing') {
						dijit.byId("filter_splice_site").set('checked', false);
						//dijit.byId("filter_splice_site").set('disabled', true);
					}
					if (lTmp4[j] == 'ncrna') {
						dijit.byId("filter_ncrna").set('checked', false);
						//dijit.byId("filter_ncrna").set('disabled', true);
					}
					if (lTmp4[j] == 'phase') {
						dijit.byId("filter_phase").set('checked', false);
						//dijit.byId("filter_phase").set('disabled', true);
					}
					if (lTmp4[j] == 'stop') {
						dijit.byId("filter_stop").set('checked', false);
						//dijit.byId("filter_stop").set('disabled', true);
					}
					if (lTmp4[j] == 'frameshift') {
						dijit.byId("filter_frameshift").set('checked', false);
						//dijit.byId("filter_frameshift").set('disabled', true);
					}
					if (lTmp4[j] == 'essential_splicing') {
						dijit.byId("filter_essential_splicing").checked = false;
						//dijit.byId("filter_essential_splicing").set('disabled', true);
					}
					if (lTmp4[j] == 'maturemirna') {
						dijit.byId("filter_maturemirna").set('checked', false);
						//dijit.byId("filter_maturemirna").set('disabled', true);
					}
					if (lTmp4[j] == 'freq_1') {
						lFreqFilters.push(lTmp4[j]);
					}
					if (lTmp4[j] == 'freq_05') {
						lFreqFilters.push(lTmp4[j]);
					}
					if (lTmp4[j] == 'freq_01') {
						lFreqFilters.push(lTmp4[j]);
					}
					if (lTmp4[j] == 'freq_001') {
						lFreqFilters.push(lTmp4[j]);
					}
				}
			}

			if (lTmp3[0] == 'dejavu') {
				change_slider_dejavu_text(lTmp3[1]);
				max_slider_dejavu = lTmp3[1];
			}

			if (lTmp3[0] == 'gnomad') {
				dijit.byId("slider_gnomad_ac").set('value', max_slider_gnomad_ac);
				if (lTmp3[1] == 'gnomad_ac_5') { dijit.byId("slider_gnomad_ac").set('value', 1); }
				if (lTmp3[1] == 'gnomad_ac_10') { dijit.byId("slider_gnomad_ac").set('value', 2); }
				if (lTmp3[1] == 'gnomad_ac_20') { dijit.byId("slider_gnomad_ac").set('value', 3); }
				if (lTmp3[1] == 'gnomad_ac_50') { dijit.byId("slider_gnomad_ac").set('value', 4); }
				if (lTmp3[1] == 'gnomad_ac_100') { dijit.byId("slider_gnomad_ac").set('value', 5); }
			}

			//    		if (lTmp3[0] == 'model') {
			//				document.getElementById("strict-denovo").checked = false;
			//				document.getElementById("denovo").checked = false;
			//				document.getElementById("dominant").checked = false;
			//				document.getElementById("recessive").checked = false;
			//				document.getElementById("xor").checked = false;
			//				document.getElementById("mosaic").checked = false;
			//				document.getElementById("both").checked = false;
			//				
			//				document.getElementById("strict-denovo").disabled = true;
			//				document.getElementById("denovo").disabled = true;
			//				document.getElementById("dominant").disabled = true;
			//				document.getElementById("recessive").disabled = true;
			//				document.getElementById("xor").disabled = true;
			//				document.getElementById("mosaic").disabled = true;
			//				document.getElementById("both").disabled = true;
			//				
			//				if (lTmp3[1] == 'strict-denovo') {
			//					document.getElementById("strict-denovo").checked = true;
			//					document.getElementById("strict-denovo").disabled = false;
			//				}
			//				if (lTmp3[1] == 'dominant') {
			//					document.getElementById("dominant").checked = true;
			//					document.getElementById("dominant").disabled = false;
			//				}
			//				if (lTmp3[1] == 'recessif') {
			//					document.getElementById("recessive").checked = true;
			//					document.getElementById("recessive").disabled = false;
			//				}
			//				if (lTmp3[1] == 'compound') {
			//					document.getElementById("xor").checked = true;
			//					document.getElementById("xor").disabled = false;
			//				}
			//				if (lTmp3[1] == 'compound_multi') {
			//					document.getElementById("xor").checked = true;
			//					document.getElementById("xor").disabled = false;
			//				}
			//				if (lTmp3[1] == 'mosaic') {
			//					document.getElementById("mosaic").checked = true;
			//					document.getElementById("mosaic").disabled = false;
			//				}
			//				if (lTmp3[1] == 'denovo') {
			//					document.getElementById("strict-denovo").checked = true;
			//					document.getElementById("denovo").checked = true;
			//					document.getElementById("strict-denovo").disabled = false;
			//					document.getElementById("denovo").disabled = false;
			//				}
			//    			if (lTmp3[1] == 'recessif_compound') {
			//    				document.getElementById("recessive").checked = true;
			//    				document.getElementById("xor").checked = true;
			//					document.getElementById("recessive").disabled = false;
			//					document.getElementById("xor").disabled = false;
			//				}
			//				if (lTmp3[1] == 'recessif_compound_multi') {
			//   				document.getElementById("recessive").checked = true;
			//    				document.getElementById("xor").checked = true;
			//					document.getElementById("recessive").disabled = false;
			//					document.getElementById("xor").disabled = false;
			//				}
			//    		}
		}
		if (lFreqFilters.length == 0) {
			dijit.byId("slider_frequence_query").set('value', 5);
			max_slider_frequency_query = 5;
		}
		if (lFreqFilters.length == 1) {
			dijit.byId("slider_frequence_query").set('value', 4);
			max_slider_frequency_query = 4;
		}
		if (lFreqFilters.length == 2) {
			dijit.byId("slider_frequence_query").set('value', 3);
			max_slider_frequency_query = 3;
		}
		if (lFreqFilters.length == 3) {
			dijit.byId("slider_frequence_query").set('value', 2);
			max_slider_frequency_query = 2;
		}
		if (lFreqFilters.length == 4) {
			dijit.byId("slider_frequence_query").set('value', 1);
			max_slider_frequency_query = 1;
		}
	}
	else {
		load_without_vectors();
	}

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
		}
	}
	xhr.send(this_args);
}

function view_genes_in_panel(html_values) {
	//	dijit.byId('dialog_genes_in_panel').hide();
	//	document.getElementById("content_dejavu_polycyto").innerHTML = html_values;
	//	dijit.byId('dialog_genes_in_panel').show();
}

function update_list_panels_name(phenotype_name) {
	for (panel_name in hash_relation_panel_phenotype) {
		var list_phenotype = hash_relation_panel_phenotype[panel_name];
		var span_panel_id = "span_panel_in_list_" + panel_name;
		span_panel_id = span_panel_id.replace(/ /g, '_');
		if (span_panel_id == "span_panel_in_list_all") {
			span_panel_id = "span_panel_in_list_All_Genes";
		}
		document.getElementById(span_panel_id).style.display = "none";
		for (var i = 0; i < list_phenotype.length; i++) {
			var this_phenotype_name = list_phenotype[i];
			if (this_phenotype_name == phenotype_name) {
				document.getElementById(span_panel_id).style.display = "block";
			}
		}
	}
	document.getElementById("span_phenotype_name").innerHTML = phenotype_name;
	document.getElementById("span_panel_name").innerHTML = 'All Genes';
	document.getElementById("span_panel_in_list_All_Genes").style.display = "block";
	if (phenotype_name == 'No Phenotype Associated') { }
	else if (phenotype_name == 'all') { }
	else { document.getElementById("span_panel_in_list_All_Panels_Genes").style.display = "block"; }
}

var origin_project_phenotype;
var hash_relation_panel_phenotype = {};

function change_gene_name_selection() {
	var this_gene_name = document.getElementById("gene_name_editor").value;
	if (this_gene_name == '') {
		gene_name = '';
		return;
	}
	gene_name = this_gene_name.toUpperCase();
	update_panel_name('Gene: ' + gene_name);
}

function view_panels_gene() {
	var this_gene_name = document.getElementById("gene_name_editor").value;
	if (this_gene_name == '') {
		gene_name = '';
		return;
	}
	dijit.byId('waiting').show();
	gene_name = this_gene_name.toUpperCase();
	var this_url = url_path + "/json_output_nodb/getPanelsFromDb.pl";
	var this_args = "project=" + project_name;
	this_args += "&panels_with_gene=" + this_gene_name;
	var xhr = new XMLHttpRequest();
	xhr.open('POST', this_url, true);
	xhr.setRequestHeader('Content-type', 'application/x-www-form-urlencoded');
	dijit.byId('waiting').show();
	xhr.onreadystatechange = function() {
		dijit.byId('waiting').show();
		if (this.readyState === XMLHttpRequest.DONE && this.status === 200) {
			var data = JSON.parse(this.response);
			dijit.byId('waiting').hide();
			document.getElementById('span_list_panels').innerHTML = data.items[0].genes_table_html;
			dijit.byId('dialog_list_panels').show();
			empty_gene_name_selection();
		}
	}
	xhr.send(this_args);
}

function empty_gene_name_selection() {
	document.getElementById("gene_name_editor").value = '';
}

function load_without_vectors() {
	setTimeout(function() {
		require(["dijit/registry"], function(registry) {
			if (registry.byId("slider_gnomad_ac_pv")) {
				setTimeout(function() {
					load_polyviewer(1);
				}, 200);
			}
			else {
				setTimeout(function() {
					if (registry.byId("slider_gnomad_ac_pv")) {
						setTimeout(function() {
							dijit.byId("slider_gnomad_ac_pv").set('value', 2);
							document.getElementById("td_title_slider_frequence_query").style.display = "none";
							document.getElementById("td_slider_frequence_query").style.display = "none";
							document.getElementById("td_title_slider_dejavu").style.display = "none";
							document.getElementById("td_slider_dejavu").style.display = "none";
							editor(1, 0, 'all');
						}, 1500);
					}
					else {
						setTimeout(function() {
							if (registry.byId("slider_gnomad_ac_pv")) {
								setTimeout(function() {
									dijit.byId("slider_gnomad_ac_pv").set('value', 2);
									document.getElementById("td_title_slider_frequence_query").style.display = "none";
									document.getElementById("td_slider_frequence_query").style.display = "none";
									document.getElementById("td_title_slider_dejavu").style.display = "none";
									document.getElementById("td_slider_dejavu").style.display = "none";
									editor(1, 0, 'all');
								}, 1500);
							}
						}, 1000);
					}
				}, 1000);
			}
		});
	}, 100);
}

function defineArrayButton() {
	var array_button = {};
	//array_button["notcosmic"] = "filter_not_cosmic";
	//array_button["cosmic"] = "filter_cosmic";
	//array_button["cnv"] = "filter_cnv";
	//array_button["deletion"] = "filter_deletion";
	//array_button["insertion"] = "filter_insertion";
	//array_button["substitution"] = "filter_substitution";
	array_button["silent"] = "filter_silent";
	array_button["utr"] = "filter_utr";
	array_button["splicing"] = "filter_splice_site";
	array_button["essential_splicing"] = "filter_essential_splicing";
	array_button["nonsynonymous"] = "filter_nonsynonymous";
	array_button["stop"] = "filter_stop";
	array_button["phase"] = "filter_phase";
	array_button["intergenic"] = "filter_intergenic";
	array_button["upstream_downstream"] = "filter_upstream_downstream";
	array_button["intronic"] = "filter_intronic";
	array_button["ncrna"] = "filter_ncrna";
	array_button["maturemirna"] = "filter_maturemirna";
	array_button["pseudogene"] = "filter_pseudogene";
	array_button["frameshift"] = "filter_frameshift";
	array_button["non-frameshift"] = "filter_nonframeshift";
	array_button["predicted_splice_site"] = "filter_predicted_splice_site";
	array_button["predicted_promoter_ai"] = "filter_predicted_splice_site";
	return (array_button);
}

function functional_url_filter_data() {
	var ab = defineArrayButton();
	var lFilters = [];
	for (z in ab) {
		if (dijit.byId(ab[z]).checked == true) {
			lFilters.push(z);
			//filters += z + "+";
		}
	}
	var filters = lFilters.join('+');
	//	if      (dijit.byId("slider_frequence_query").get("value") == 1)    { filters += "freq_001+freq_01+freq_05+freq_1+"; }
	//	else if (dijit.byId("slider_frequence_query").get("value") == 2)    { filters += "freq_01+freq_05+freq_1+"; }
	//	else if (dijit.byId("slider_frequence_query").get("value") == 3)    { filters += "freq_05+freq_1+"; }
	//	else if (dijit.byId("slider_frequence_query").get("value") == 4)    { filters += "freq_1+"; }
	return filters;
}

var panel_selected = 'All Genes';
function update_panel_name(panel_name) {
	if (String(panel_name).match(/Phenotype: /)) {
		return;
	}
	panel_name = String(panel_name).replace("<label><u>", "");
	panel_name = String(panel_name).replace("</u></label>", "");
	panel_name = String(panel_name).replace("<span style='color:green;'>", "");
	panel_name = String(panel_name).replace("</span>", "");
	if (panel_name.match(/Gene:/g)) {
		panel_selected = 'all';
		document.getElementById("span_phenotype_name").innerHTML = panel_name;
		document.getElementById("span_panel_name").innerHTML = panel_name;
		load_polyviewer(1);
		return;
	}
	gene_name = '';
	panel_selected = panel_name;
	if (panel_name == 'all') {
		document.getElementById("span_phenotype_name").innerHTML = 'All Genes';
		document.getElementById("span_panel_name").innerHTML = 'All Genes';
	}
	else {
		document.getElementById("span_panel_name").innerHTML = panel_name;
	}
	load_polyviewer(1);
	return;
}

function update_grid_gene_phenotypes(gene_id, this_proj_name) {
	document.getElementById("content_table_phenotypes").innerHTML = '';
	dijit.byId('waiting').show();
	get_gene_phenotypes(gene_id, this_proj_name);
}

// Need: <script src="https://cdn.anychart.com/releases/v8/js/anychart-base.min.js"></script>
// Need: <script src="https://cdn.anychart.com/releases/v8/js/anychart-tag-cloud.min.js"></script>	
function add_wordcloud(container_name, json_data) {
	document.getElementById(container_name).innerHTML = '';
	anychart.onDocumentReady(function() {
		var data = JSON.parse(json_data);
		var chart = anychart.tagCloud(data);// set a chart title
		chart.title('Most described words in HGMD / HPO / OMIM')
		// set an array of angles at which the words will be laid out
		chart.angles([0, -45, 90]);
		// enable a color range
		chart.colorRange(true);
		// set the color range length
		chart.colorRange().length('80%');// display the word cloud chart
		chart.container(container_name);
		chart.draw();
	});
	return;
}

function get_gene_phenotypes(gene_id, this_proj_name) {
	var this_url = url_path + "/json_output_nodb/getGenePhenotypes.pl";
	var this_args;
	if (this_proj_name) { this_args = "project=" + this_proj_name + "&gene=" + gene_id; }
	else { this_args = "gene=" + gene_id; }
	var xhr = new XMLHttpRequest();
	xhr.open('POST', this_url, true);
	xhr.setRequestHeader('Content-type', 'application/x-www-form-urlencoded');
	xhr.onreadystatechange = function() {
		if (this.readyState === XMLHttpRequest.DONE && this.status === 200) {
			var data = JSON.parse(this.responseText);
			var store_tmp = new dojo.data.ItemFileWriteStore({
				data: {
					identifier: 'name',
					label: 'name',
					items: data.items
				}
			});
			store_tmp.fetch({
				onItem: function(item, result) {
					document.getElementById("content_table_phenotypes").innerHTML = item.html;
					if (item.txt_all != '-') {
						add_wordcloud('container_wordcloud', item.txt_all);
						document.getElementById("container_wordcloud").style.display = "block";
					}
					else {
						document.getElementById("container_wordcloud").style.display = "none";
					}
					dijit.byId('waiting').hide();
					dijit.byId('dialog_list_phenotypes').show();
				}
			});
		}
	}
	xhr.send(this_args);
}

function construct_url_filter_data() {
	var filter_args = "";
	var args = new Array();
	filter_args = "&filter_type_variation=" + functional_url_filter_data();
	var dejavu_value = get_value_slider_dejavu();
	if (dejavu_value == 0) { filter_args += "&dejavu=uniq"; }
	else if (dejavu_value < 100) { filter_args += "&dejavu=" + dejavu_value; }
	if (dijit.byId("dejavu_ho").get("value") == "off") { filter_args += "&dejavu_ho=1"; }

	if (dijit.byId("slider_gnomad_ac").get("value") == 1) { filter_args += "&gnomad=gnomad_ac_5"; }
	if (dijit.byId("slider_gnomad_ac").get("value") == 2) { filter_args += "&gnomad=gnomad_ac_10"; }
	if (dijit.byId("slider_gnomad_ac").get("value") == 3) { filter_args += "&gnomad=gnomad_ac_20"; }
	if (dijit.byId("slider_gnomad_ac").get("value") == 4) { filter_args += "&gnomad=gnomad_ac_50"; }
	if (dijit.byId("slider_gnomad_ac").get("value") == 5) { filter_args += "&gnomad=gnomad_ac_100"; }

	if (dijit.byId("slider_allele_quality").get("value") == 1) { filter_args += "&ratio=ratio_80"; }
	if (dijit.byId("slider_allele_quality").get("value") == 2) { filter_args += "&ratio=ratio_40"; }
	if (dijit.byId("slider_allele_quality").get("value") == 3) { filter_args += "&ratio=ratio_20"; }
	if (dijit.byId("slider_allele_quality").get("value") == 4) { filter_args += "&ratio=ratio_10"; }

	var list_delete_models = [];

	//	if(document.getElementById('strict-denovo').checked) {}
	//	else { list_delete_models.push("fam_strict_denovo"); }

	if (document.getElementById('denovo').checked) { }
	else { list_delete_models.push("fam_denovo"); }

	//	if(document.getElementById('dominant').checked) {}
	//	else { list_delete_models.push("fam_dominant"); }

	//	if(document.getElementById('mosaic').checked) {}
	//	else { list_delete_models.push("fam_mosaic"); }

	if (document.getElementById('recessive').checked) { }
	else { list_delete_models.push("fam_recessif"); }

	if (document.getElementById('both').checked) { }
	else { list_delete_models.push("fam_both_parents"); }

	//compound

	//PANELS
	if (panel_selected == 'All Genes') { }
	else {
		filter_args += "&panel=" + panel_selected;
	}

	if (list_delete_models.length > 0) {
		filter_args += "&delete_models=" + list_delete_models.join(' ');
	}

	return filter_args;
}

function  reload_tab_polydiag() {
	for (var ii = 0; ii < tsamples.length; ii++) {
		var name = tsamples[ii].label
		var turl = polydiag_url(name);
		tab_editor[name].setUrl = turl;
		tab_editor[name].toto = 1;
		if (tab_editor[name].get('selected')) { tab_editor[name].onShow() }
		}
}


function refresh_tab_polydiag() {
	if (!(dijit.byId("gridPatients"))) {
		return 1;
	}
	
	if (view_polydiag == 0){
		return 1;
	}
	
	var tsamples2 = dijit.byId("gridPatients").selection.getSelected();
	var selected = new Array();
	for (var i = 0; i < tsamples2.length; i++) {
		selected[tsamples2[i].label] = 1;
	}
	var first = 0;
	for (var i = 0; i < tsamples.length; i++) {
		var name = tsamples[i].label
		tab_editor[name].toto = 0;
		}
	for (var i = 0; i < tsamples.length; i++) {
		var name = tsamples[i].label
		tab_editor[name].toto = 1;
		if (selected[name] == 1) {
			
			if (tab_editor[tsamples[i].label].add == 0) {
				//console.log(tsamples[i].label);
				container_editor.addChild(tab_editor[tsamples[i].label]);
				tab_editor[tsamples[i].label].add = 1;
			}
			//tab_editor[tsamples[i].label].setUrl = turl;
			if (first < 1){
				tab_editor[tsamples[i].label].onShow();
				first ++;
			}
			
		}
		else {
			
			if (tab_editor[tsamples[i].label].add == 1) {
				
				tab_editor[name].set('selected',false);
				container_editor.removeChild(tab_editor[tsamples[i].label]);
				
				tab_editor[tsamples[i].label].add = 0;
			}
		}

	}



}


function load_tab_polydiag() {
	for (var i = 0; i < tsamples.length; i++) {
		var name = tsamples[i].label;
		container_editor.addChild(tab_editor[name]);
		tab_editor[name].add = 1;
		tab_editor[name].toto = 1;
	}
	tab_editor[tsamples[0].label].onShow();
	if (var_load_polyviewer == 1){ refresh_tab_polydiag();}
	//else {}
	
}

function polydiag_url (patient) {
	var qual = dijit.byId("ratio_label").get("value");
	if (qual) {
		qual = qual * -1;
	}
	else {
		qual = dijit.byId("slider_allele_quality").value;
		dijit.byId("ratio_label").set("value", "");
	}
	var url = url_report_patient + "?project=" + project_name + "&patients=" + patient+"&panel=" + panel + "&edit_mode=1"+ "&never=" + dijit.byId("check_never").get('value') + "&this=" + dijit.byId("slider_this").value + "&impact=" + dijit.byId("slider_consequence").value + "&frequence=" + dijit.byId("slider_frequence").value + "&allele_quality=" + qual + "&report_mode=1";
    var transcripts = dijit.byId("gridGenes").selection.getSelected();
	var t = "";
	if (transcripts.length == nb_genes) {
		t = "all";

	}
	else {
		for (var i = 0; i < transcripts.length; i++) {
			t += transcripts[i].transcript + ",";
		}
	}
	url += "&transcripts=" + t;
	return url;
}

	function return_selected_patient_polydiag(tab){
		for (var p in tab ){
			//alert(tab[p].get('selected')+" "+p+" polydiag");
			if (tab[p].get('selected')) { return p; }
		}
		return "";
	}


function downloadFile_inprogress() {
//    document.getElementById("wait2").src="/icons/Polyicons/12-em-check.png";
//    document.getElementById("waiting_xls_progress_text").innerHTML = "<b> Downloading XLS file</b><br>Sometimes, your browser can download your XLS files in your <i><b><u>Download</i></b></u> or <i><b><u>Telechargements</i></b></u> directory automatically.";
}

function downloadXlsFiles(base_url, tsamples) {
	var nb_patient = 0;
	downloadFile(base_url, tsamples, nb_patient);
}

function checkBlankOpen(base_url, tsamples, next_nb_patient, w) {
	setTimeout(function(){
		console.log(w);
		if (w == null) {
			downloadFile(base_url, tsamples, next_nb_patient);
			return 1;
		}
		if (w.closed) { 
			downloadFile(base_url, tsamples, next_nb_patient);
			return 1;
		}
		else {
			checkBlankOpen(base_url, tsamples, next_nb_patient, w);
		}
	}, 1000);
}

//dojo.require("dojo.io.iframe");
function downloadFile(base_url, tsamples, nb_patient) {
	if (nb_patient == tsamples.length) {
		alert('Export XLS done'); 
		return;
	}
	var patient= tsamples[nb_patient].label;
	var next_nb_patient = nb_patient + 1;
	
	console.log("\n\nLAUNCH "+patient);
	var promise = new Promise(
		function(resolve, reject) {
			setTimeout(function(){
				var url = base_url + '&patients=' + patient;
				var w = window.open(url);
				resolve (checkBlankOpen(base_url, tsamples, next_nb_patient, w));
			}, 500);
		}
	);
} 

function editor(mode, type) {
	if (is_project_update == 'yes') {
		if (panel_editor_stack_changed == 0) {
			dijit.byId("editor_stack").set("href", "include_html/diag/editor_container_polyviewer.html");
			dijit.byId("editor_stack").refresh();
			update_buttons_polyquery();
			panel_editor_stack_changed = 1;
		}
		return editor_update(mode, type);
	}

	var box1 = dijit.byId("ratio_label");

	var qual = box1.get("value");
	if (qual) {
		qual = qual * -1;


	}
	else {
		qual = dijit.byId("slider_allele_quality").value;
		box1.set("value", "");
	}
	var url = url_report_patient + "?project=" + project_name + "&panel=" + panel + "&edit_mode=" + mode + "&never=" + dijit.byId("check_never").get('value') + "&this=" + dijit.byId("slider_this").value + "&impact=" + dijit.byId("slider_consequence").value + "&frequence=" + dijit.byId("slider_frequence").value + "&allele_quality=" + qual + "&report_mode=1";
	var argsPost = {
		project: project_name,
		"panel": panel,
		edit_mode: mode,
		never: dijit.byId("check_never").get('value'),
		impact: dijit.byId("slider_consequence").value,
		frequence: dijit.byId("slider_frequence").value,
		"this": dijit.byId("slider_this").value,
		allele_quality: qual,
		report_mode: 1,
	};


	//var url = "http://darwin.bipd.fr/cgi-bin/pnitschk/polymorphism-cgi/validation_variation/patient_cache.pl?project="+project_name+"&edit_mode="+mode+"&never="+dijit.byId("check_never").get('value')+"&this="+dijit.byId("slider_this").value+"&impact="+dijit.byId("slider_consequence").value+"&frequence="+dijit.byId("slider_frequence").value+"&allele_quality="+dijit.byId("slider_allele_quality").value+"&report_mode=1";

	if (type == 1) {
		/*report alacarte*/
		url += get_alacarte();
	}
	if (type == 2) {
		/*report alacarte*/
		url += get_alacarte() + "&xls=1";
	}
	if (type == 3) {
		/*report alacarte*/
		url += get_alacarte() + "&list=1";
	}
	else {
	
		url += "&project_summary=1&all_coverage=1&table_variations=1&sanger_variations=1&validated_variations=1&all_variations=1";
		argsPost["project_summary"] = 1;
		argsPost["all_coverage"] = 1;
		argsPost["table_variations"] = 1;
		argsPost["sanger_variations"] = 1;
		argsPost["validated_variations"] = 1;
		argsPost["all_variations"] = 1;
		if (dijit.byId("denovo").checked == true) url += "&denovo=1";
		if (dijit.byId("denovo").checked == true) argsPost["denovo"] = 1;;
		if (dijit.byId("recessive").checked == true) url += "&recessive=1";
		if (dijit.byId("recessive").checked == true) argsPost["recessive"] = 1;
		if (dijit.byId("xor").checked == true) url += "&xor=1";
		if (dijit.byId("xor").checked == true) argsPost["xor"] = 1;
		if (dijit.byId("both").checked == true) url += "&both=1";
		if (dijit.byId("both").checked == true) argsPost["both"] = 1;


	}
	url += return_options();

	var tab_selected_patient;
	if (tab_editor) {
		tab_selected_patient = return_selected_patient_polydiag(tab_editor);
	}
	/* sample */


	var transcripts = dijit.byId("gridGenes").selection.getSelected();
	var t = "";
	var tz = "";
	var total_show = 0;
	if (transcripts.length == nb_genes) {
		t = "all";
		tz = "all";

	}
	else {
		for (var i = 0; i < transcripts.length; i++) {
			t += transcripts[i].transcript + ",";
		}
	}

	//document.getElementById("label_patients").innerHTML = p;
	url += "&transcripts=" + t;
	argsPost["transcripts"] = t;
	var username = dojo.cookie("username");

	url += "&user_name=" + username;
	argsPost["user_name"] = username;
	if (type == 2) {
		window.location = url + "&patients=" + tab_selected_patient;
		return;
	}
	if (type == 3) {
		var tsamples = dijit.byId("gridPatients").selection.getSelected();
		console.log('INIT XLS');
		var base_url = url + "&xls=1";
		console.log('base URL '+base_url);
		downloadXlsFiles(base_url, tsamples);
		return;
		
//		for (var i = 0; i < tsamples.length; i++) {
//			var n = tsamples[i].label;
//			console.log('downloadind '+n);
//			downloadFiles(patients=" + n);
//		}
//		return;

	}

	var tsamples = dijit.byId("gridPatients").selection.getSelected();

	//  return;
	var p = "";
	if (tsamples.length == 0) {
		p = "all";
	}
	else {
		var z = 0;
		var num_samples = 0;
		tsamples.sort(function(a, b) {
			return b.child > a.child;
		});
		for (var i = 0; i < tsamples.length; i++) {
			if (tsamples[i].child == 0) continue;
			num_samples++;
		}

		var n_child = 0;
		for (var i = 0; i < tsamples.length; i++) {
			var n = tsamples[i].label;

			var n = tsamples[i].label;

			//   if  (tsamples[i].child == 0) continue;
			var ic = '<img src="images/polyicons//iconfinder_baby-boy_s.png">  &nbsp;';

			ic = '<img src="images/polyicons/iconfinder_baby-girl_16.png">  &nbsp;';

			if (tsamples[i].child == 0) {
				ic = '<img src="images/polyicons/icons8-person-24.png">&nbsp;';
				if (tsamples[i].sex == "F") ic = '<img src="images/polyicons/icons8-person-female-24.png">&nbsp;';
			}
			else {

				if (tsamples[i].sex == "F" && tsamples[i].status == "1") {
					ic = '<img src="images/polyicons//iconfinder_baby-girl_s.png">  &nbsp;';
				}
				if (tsamples[i].sex == "F" && tsamples[i].status == "2") {
					ic = '<img src="images/polyicons//iconfinder_baby-girl_d.png">  &nbsp;';
				}
				if (tsamples[i].sex == "M" && tsamples[i].status == "1") {
					ic = '<img src="images/polyicons//iconfinder_baby-boy_s.png">  &nbsp;';
				}
				if (tsamples[i].sex == "M" && tsamples[i].status == "2") {
					ic = '<img src="images/polyicons//iconfinder_baby-boy_d.png">  &nbsp;';
				}

			}
			//#FDAC53
			var star = "";
			//;
			if (tsamples[i].validation_status == -99) {

				star = "<span class=' glyphicon glyphicon-exclamation-sign' aria-hidden='true' style='position:absolute;top:1px;left:1px;color:#BA0BE0'></span>";
			}
			else if (tsamples[i].validation_status >= 3) {
				star = "<span class='glyphicon glyphicon-star' aria-hidden='true' style='position:absolute;top:1px;left:1px;color:#E74C3C'></span>";
			}
			else if (tsamples[i].validation_status >= 1) {
				star = "<span class='glyphicon glyphicon-star' aria-hidden='true' style='position:absolute;top:1px;left:1px;color:black'></span>";
			}






			if (tsamples[i].child == 0) continue;
			n_child++;
			var ic = '<i class="fa fa-mars " aria-hidden="true" style="color:blue;vertical-align: top"></i> &nbsp';
			if (tsamples[i].sex == "F") ic = '<i class="fa fa-venus " aria-hidden="true" style="color:#F660AB;vertical-align: top"></i> &nbsp';
			/*     if (!tab_editor[n]) {
					  z++;
				   tab_editor[tsamples[i].label] = new dijit.layout.ContentPane({
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
							 
							 
							// this.setHref(url_report_patient);
							 if (this.toto== 2){
								// this.ioArgs.content = this.setArgsPost; 
								 //console.log(url_report_patient );
								 //console.log(this.setArgsPost );
								 //console.log(this.setUrl);
								 //this.setHref(url_report_patient);
								this.setHref(this.setUrl);
							}
					 }
				  });
				  
			  tab_editor[tsamples[i].label].add = 1;
			  container_editor.addChild(tab_editor[tsamples[i].label]);
			  }
	  */

			tab_editor[tsamples[i].label].setUrl = url + "&patients=" + tsamples[i].label;
			tab_editor[tsamples[i].label].setArgsPost = argsPost;
			tab_editor[tsamples[i].label].setArgsPost.patients = tsamples[i].label;
			tab_editor[tsamples[i].label].transcripts = t;
			tab_editor[tsamples[i].label].toto = 1;


		}

	}

	url += "&patients=" + p;
	//container_editor.selectChild( tab_editor[tsamples[0].label]);

	type = "patient";
	//alert(tab_selected_patient);
	if (tab_selected_patient) {
	//	tab_editor[tab_selected_patient].onShow();

	}
	else {
		//tab_editor[tsamples[0].label].onShow();
	}

}








function patient(mode, type) {
	var url = url_report_patient + "?project=" + project_name + "&panel=" + panel + "&edit_mode=" + mode + "&report_mode=1" + "&this=" + dijit.byId("select_this").get("value") + "&impact=" + dijit.byId("select_impact").get("value") + "&frequence=" + dijit.byId("select_frequence").get("value");
	if (type == 1) {
		/*report alacarte*/
		url += get_alacarte();
	}
	else {
		url += "report_mode=1&project_summary=1&all_coverage=1&table_variations=1&sanger_variations=1&validated_variations=1&all_variations=1";
	}
	url += return_options();

	/* sample */


	var transcripts = dijit.byId("gridGenes").selection.getSelected();
	var t = "";
	var tz = "";
	var total_show = 0;
	if (transcripts.length == nb_genes) {
		t = "all";
		tz = "all";

	}
	else {
		for (var i = 0; i < transcripts.length; i++) {
			t += transcripts[i].transcript + ",";

		}
	}
	//document.getElementById("label_patients").innerHTML = p;
	url += "&transcripts=" + t;

	var username = dojo.cookie("username");

	url += "&user_name=" + username;

	var tsamples = dijit.byId("gridPatients").selection.getSelected();
	var tab_selected_patient;
	if (tab_editor) {
		tab_selected_patient = return_selected_patient(tab_patients);
	}

	//	return;
	var p = "";
	if (tsamples.length == 0) {
		p = "all";
	}
	else {
		var z = 0;
		for (var i = 0; i < tsamples.length; i++) {
			var n = tsamples[i].label;
			if (!tab_patients[n]) {
				z++;
				tab_patients[tsamples[i].label] = new dijit.layout.ContentPane({
					title: tsamples[i].label,

					content: " ",
					closable: true,
					id: "c1p" + n,
					onClose: function() {
						var t = this.get('title');
						delete tab_patients[t];
						return true;
					},
					onShow: function() {

						this.toto++;
						if (this.toto == 2) {


							this.setHref(this.setUrl);
						}
					}
				});


				container_patient.addChild(tab_patients[tsamples[i].label]);
			}
			tab_patients[tsamples[i].label].setUrl = url + "&patients=" + tsamples[i].label;
			tab_patients[tsamples[i].label].setArgsPost = argsPost;//
			tab_patients[tsamples[i].label].setArgsPost.patients = tsamples[i].label;

			if (tsamples.length == 1) {
				tab_patients[tsamples[i].label].ioArgs.content = { test: 'it should work' };
				tab_patients[tsamples[i].label].setHref(url + "&patients=" + tsamples[i].label);
				tab_patients[tsamples[i].label].toto = 2;
			}
			else if ((tsamples.length - i) < 3) {
				tab_patients[tsamples[i].label].toto = 1;

				container_patient.selectChild(tab_patients[tsamples[i].label]);
				//tab_patients[tsamples[i].label].show();
			}
			else {
				console.log("la");
				tab_patients[tsamples[i].label].toto = 1;
			}

		}
	}
	url += "&patients=" + p;

	if (tab_selected_patient) {
		container_patient.selectChild(tab_patients[tab_selected_patient]);
	}

	type = "patient";

}

function editor_update(mode, type) {
	alert("editor_update");
	var box1 = dijit.byId("ratio_label");
	var qual = dijit.byId("slider_allele_quality").value;


	var panel_name;
	if (document.getElementById("span_panel_name")) {
		panel_name = document.getElementById("span_panel_name").innerHTML;
	}
	else if (vectors_hex) {
		document.getElementById("span_panel_name").innerHTML = "all";
	}
	else {
		panel_name = "sysid";
	}
	if (panel_name == 'All Genes') {
		panel_name = 'all';
	}

	var url;
	var argsPost;

	if (vectors_hex) {
		url += "&vectors_hex='" + vectors_hex + "'";
		argsPost = {
			project: project_name,
			edit_mode: 1,
			never: 1,
			this: 6,
			impact: 3,
			frequence: 4,
			allele_quality: 5,
			report_mode: 1,
			project_summary: 1,
			all_coverage: 1,
			table_variations: 1,
			sanger_variations: 1,
			validated_variations: 1,
			all_variations: 1,
			strict_denovo: 1,
			denovo: 1,
			recessive: 1,
			xor: 1,
			both: 1,
			utr: 1,
			intronic: 1,
			all: 1,
			span: 20,
			limit: 30,
		};

		argsPost["vectors_hex"] = vectors_hex;
		argsPost["query_filters"] = "'" + polyquery_filters + "'";


		//dijit.byId("slider_gnomad_ac").set('value', 6);
		//dijit.byId("slider_consequence").set('value', 1);
		//dijit.byId("slider_frequence").set('value', 5);
		//dijit.byId("slider_consequence").set('disabled', 'disabled');
		//dijit.byId("slider_frequence").set('disabled', 'disabled');
		//dijit.byId("slider_this").set('disabled', 'disabled');
		//dijit.byId("slider_gnomad_ac").set('disabled', 'disabled');
		//dijit.byId("slider_allele_quality").set('disabled', 'disabled');

		url = "project=" + project_name + "&edit_mode=1&never=1&this=6&impact=3&frequence=4&allele_quality=5&report_mode=1&project_summary=1&all_coverage=1&table_variations=1&sanger_variations=1&validated_variations=1&all_variations=1&denovo=1&recessive=1&xor=1&both=1&utr=1&intronic=1&all=1";
		url += '&vectors_hex=' + vectors_hex;
	}
	else {
		//url = url_report_patient+"?project="+project_name+"&edit_mode="+mode+"&never="+dijit.byId("check_never").get('value')+"&this="+dijit.byId("slider_this").value+"&impact="+dijit.byId("slider_consequence").value+"&frequence="+dijit.byId("slider_frequence").value+"&allele_quality="+qual+"&report_mode=1";

		url = url_report_patient + "?project=" + project_name;
		url += "&edit_mode=" + mode;
		url += "&never=" + dijit.byId("check_never").get('value');
		//url += "&this="+dijit.byId("slider_this").value;
		//url += "&frequence="+dijit.byId("slider_frequence_query").value;
		url += "&allele_quality=" + qual;
		url += "&report_mode=1";
		url += "&annot=" + functional_url_filter_data();

		argsPost = {
			project: project_name,
			edit_mode: mode,
			never: dijit.byId("check_never").get('value'),
			//impact: dijit.byId("slider_consequence").value,
			//frequence: dijit.byId("slider_frequence_query").value,
			ac: dijit.byId("slider_gnomad_ac").value,
			allele_quality: qual,
			report_mode: 1,
			annot: functional_url_filter_data(),
		};

		if (panel_name) {
			url += "&panel=" + panel_name;
			argsPost["panel"] = panel_name;
		}

		if (type == 1) {
			/*report alacarte*/
			url += get_alacarte();
		}
		if (type == 2) {
			/*report alacarte*/
			url += get_alacarte() + "&xls=1";
		}
		if (type == 3) {
			/*report alacarte*/
			url += get_alacarte() + "&list=1";
		}
		else {
			url += "&project_summary=1&all_coverage=1&table_variations=1&sanger_variations=1&validated_variations=1&all_variations=1";
			argsPost["project_summary"] = 1;
			argsPost["all_coverage"] = 1;
			argsPost["table_variations"] = 1;
			argsPost["sanger_variations"] = 1;
			argsPost["validated_variations"] = 1;
			argsPost["all_variations"] = 1;
			if (document.getElementById("strict_denovo").checked == true) url += "&strict_denovo=1";
			if (document.getElementById("strict_denovo").checked == true) argsPost["strict_denovo"] = 1;
			if (document.getElementById("denovo").checked == true) url += "&denovo=1";
			if (document.getElementById("denovo").checked == true) argsPost["denovo"] = 1;
			if (document.getElementById("recessive").checked == true) url += "&recessive=1";
			if (document.getElementById("recessive").checked == true) argsPost["recessive"] = 1;
			if (document.getElementById("xor").checked == true) url += "&xor=1";
			if (document.getElementById("xor").checked == true) argsPost["xor"] = 1;
			if (document.getElementById("both").checked == true) url += "&both=1";
			if (document.getElementById("both").checked == true) argsPost["both"] = 1;
		}
	}

	url += return_options();

	var tab_selected_patient;
	if (tab_editor) {
		tab_selected_patient = return_selected_patient(tab_editor);
	}
	/* sample */

	if (vectors_hex) {
		url += "&transcripts=all";
		argsPost["transcripts"] = 'all';
	}
	else {
		//var transcripts = dijit.byId("gridGenes").selection.getSelected();
		var t = "";
		var tz = "";
		var total_show = 0;
		//    if (transcripts.length == nb_genes){
		t = "all";
		tz = "all";

		//    }
		/* else {
		 for (i=0;i<transcripts.length;i++){
		 t += transcripts[i].transcript+",";
		 
		 }
		 }
	 */
		//document.getElementById("label_patients").innerHTML = p;
		url += "&transcripts=" + t;
		argsPost["transcripts"] = t;
	}

	url += "&user_name=" + username;
	argsPost["user_name"] = username;
	if (type == 2) {
		window.location = url + "&patients=" + tab_selected_patient;
		return;
	}
	if (type == 3) {
		var tsamples = dijit.byId("gridPatients").selection.getSelected();
		for (var i = 0; i < tsamples.length; i++) {
			var n = tsamples[i].label;
			window.open(url + "&xls=1&patients=" + n);
		}
		return;
	}

	// var tsamples = dijit.byId("gridPatients").selection.getSelected();
	var p = "all";
	var p = "";
	var tsamples = storeP._arrayOfAllItems;//dijit.byId("gridPatients").selection.getSelected();
	if (tsamples.length == 0) {
		p = "all";
	}
	else {
		var z = 0;
		var num_samples = 0;
		for (var i = 0; i < tsamples.length; i++) {
			if (tsamples[i].child == 0) continue;
			num_samples++;
		}

		for (var i = 0; i < tsamples.length; i++) {
			var n = tsamples[i].label;
			if (tsamples[i].child == 0) continue;
			var ic = '<i class="fa fa-mars " aria-hidden="true" style="color:blue;vertical-align: top"></i> &nbsp';
			if (tsamples[i].sex == "F") ic = '<i class="fa fa-venus " aria-hidden="true" style="color:#F660AB;vertical-align: top"></i> &nbsp';
			if (!tab_editor[n]) {
				z++;
				tab_editor[tsamples[i].label] = new dijit.layout.ContentPane({
					title: ic + tsamples[i].label,
					content: " ",
					closable: true,
					id: "c1pp" + n,
					ioMethod: dojo.xhrPost,
					onClose: function() {
						var t = this.get('title');
						delete tab_patients[t];
						return true;
					},
					onShow: function() {
						this.toto++;
						//     this.setHref(url_report_patient_cache);
						if (this.toto == 2) {
							this.setArgsPost.patients = this.patient;
							this.ioArgs.content = this.setArgsPost;
							console.log(url_report_patient_cache);
							console.log(this.setArgsPost);
							console.log(this.setUrl);

							this.setHref(url_report_patient_cache);
							//this.setHref(this.setUrl);
						}
					}
				});


				container_editor.addChild(tab_editor[tsamples[i].label]);
			}
			tab_editor[tsamples[i].label].debug = url + "&patients=" + tsamples[i].label;
			tab_editor[tsamples[i].label].setArgsPost = argsPost;
			tab_editor[tsamples[i].label].patient = tsamples[i].label;
			tab_editor[tsamples[i].label].transcripts = t;
			if (num_samples == 1) {
				tab_editor[tsamples[i].label].setArgsPost = argsPost;

				tab_editor[tsamples[i].label].patient = tsamples[i].label;
				tab_editor[tsamples[i].label].setArgsPost.patients = tsamples[i].label;
				//  tab_editor[tsamples[i].label].setHref(url+"&patients="+tsamples[i].label);

				tab_editor[tsamples[i].label].ioArgs.content = tab_editor[tsamples[i].label].setArgsPost;
				tab_editor[tsamples[i].label].setHref(url_report_patient_cache);
				tab_editor[tsamples[i].label].toto = 2;
			}
			else if ((num_samples - i) < 2) {


				tab_editor[tsamples[i].label].toto = 1;
				container_editor.selectChild(tab_editor[tsamples[i].label]);
			}
			else {
				tab_editor[tsamples[i].label].toto = 1;
			}

		}
	}
	url += "&patients=" + p;



	type = "patient";
	if (tab_selected_patient) {
		container_editor.selectChild(tab_editor[tab_selected_patient]);
	}

}


var selected_panel = 'DI_DI44';

function selectPanel(panel_name) {
	selected_panel = panel_name;
}

var first_panel_name = '';
function load_dropdown_menu_panels_coverage() {
	//var lPanels = ['all'];
	var lPanels = [];
	var this_url = url_path + "/json_output_nodb/getPanelsFromDb.pl";
	var this_args = "project=" + project_name + "&for_coverage_panel=1";
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
			var lPanels = [];
			var lPhenotypes = [];
			var lNbGenes = [];
			store_tmp.fetch({
				onItem: function(item, result) {
					var name = item.name[0];
					if (name == 'All Genes') { }
					else {
						lPanels.push(name);
						lPhenotypes.push(item.phenotype_name[0]);
						lNbGenes.push(item.nb_genes[0]);
					}
				}
			});
			for (var i = 0; i < lPanels.length; i++) {
				var panel_name = lPanels[i];
				var span = document.createElement('span');
				if (panel_name == 'all') {
					continue;
				}
				else if (panel_name == 'Omim-morbid' || panel_name == 'ACMG-Actionable') {
					span.innerHTML = "<span class='glyphicon glyphicon-blackboard' style='color:#cca919'></span> " + panel_name;
				}
				else {
					continue;
				}
				var div = document.createElement('div');
				div.appendChild(span);
				div.setAttribute("onClick", "selectPanel('" + panel_name + "');load_coverage('" + panel_name + "', '" + panel_name + "');");
				//document.getElementById('dropdown-menu-panels_cov').appendChild(div);
			}
			var old_phenotype_name = '';
			for (var i = 0; i < lPanels.length; i++) {
				var phenotype_name = lPhenotypes[i];
				if (phenotype_name == old_phenotype_name) { }
				else {
					var span_p = document.createElement('span');
					span_p.innerHTML = "<br><b><u>Phenotype " + phenotype_name + "</b></u>";
					var div_p = document.createElement('div');
					div_p.appendChild(span_p);
					//					document.getElementById('dropdown-menu-panels_cov').appendChild(div_p);
					old_phenotype_name = phenotype_name;
				}
				var panel_name = lPanels[i];
				var nb_genes = lNbGenes[i];
				var span = document.createElement('span');
				if (panel_name == '') { continue; }
				if (panel_name == 'all') {
					continue;
				}
				else if (panel_name == 'Omim-morbid' || panel_name == 'ACMG-Actionable') {
					continue;
				}
				else {
					if (first_panel_name == '') {
						first_panel_name = panel_name;
						//document.getElementById("span_panel_coverage").innerHTML = first_panel_name;
					}
					span.innerHTML = "<span class='glyphicon glyphicon-blackboard' style='color:#28a7e0'></span> " + panel_name;
					if (nb_genes) { span.innerHTML += " (" + nb_genes + " genes) "; }
				}
				var div = document.createElement('div');
				div.appendChild(span);
				div.setAttribute("onClick", "load_coverage('" + panel_name + "', '" + panel_name + "');");
				if (document.getElementById('dropdown-menu-panels_cov')) {
					document.getElementById('dropdown-menu-panels_cov').appendChild(div);
				}
			}
			if (document.getElementById("span_panel_coverage")) {
				document.getElementById("span_panel_coverage").innerHTML = 'DI_DI44';
			}
			dijit.byId('waiting').hide();
			//			if (first_panel_name != "none") {
			//				preload_coverage();
			//				//selectPanel(first_panel_name);
			//				load_coverage(first_panel_name, first_panel_name);
			//			}
		}
	}
	xhr.send(this_args);
}

//function myFocusFunction(div_id) {
//	document.getElementById(div_id).style.backgroundColor = "yellow";
//}

//function myBlurFunction(div_id) {
//	document.getElementById(div_id).style.backgroundColor = "";
//}

function load_coverage(type, panel) {

	document.getElementById("span_panel_coverage").innerHTML = type;
	var url = construct_url(url_coverage_overview_defidiag);
	var this_order;
	if (dijit.byId("ENST_order").get("value") == "on") { this_order = "name"; }
	else { this_order = "position"; }
	url += '&order=' + this_order;
	url += "&only_low=" + dijit.byId("only_low").get('value');
	url += "&panel_name=" + panel;
	tab_coverage[1].setHref(url);
}

function load_coverage_gene() {
	var gene_name = document.getElementById('gene_name_coverage_editor').value;
	var url = construct_url(url_coverage_overview_defidiag);
	var this_order;
	if (dijit.byId("ENST_order").get("value") == "on") { this_order = "name"; }
	else { this_order = "position"; }
	url += '&order=' + this_order;
	url += "&only_low=" + dijit.byId("only_low").get('value');
	url += "&gene=" + gene_name;
	tab_coverage[1].setHref(url);
	dijit.byId('dialog_cov_view_gene').hide();
}

function collapse_panel(p, list, runid) {
	p = p + "";
	var ps = list.split(",");


	for (var i = 0; i < ps.length; i++) {
		var l;

		if (ps[i] === p) { continue; };

		var d = document.getElementById(ps[i] + "_" + runid);
		if (d == null) {
			continue;
		}
		if (d.classList.contains("in")) {
			d.classList.remove("in");
		}
	}
	var d = document.getElementById(p + "_" + runid);

	if (d.classList.contains("in")) {

		d.classList.remove("in");
	}
	else {
		d.classList.add("in");
	}
}

function project_printer(project, user) {
	var url = url_summary + "?project=" + project + "&user=" + user + "&print=1";
	child = window.open(url); //Open the child in a tiny window.
	window.focus(); //Hide the child as soon as it is opened.

}

function showTranscripts(array_cell, row) {
	dojo.forEach(array_cell, function(cell, i) {
		dojo.setStyle(document.getElementById(cell), "display", "table-row");
	});
	if (row) {
		dojo.setStyle(document.getElementById(row), "display", "none");
	}
}


