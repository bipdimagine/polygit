/**
 * @author pnitschk
 */


		//fonction pour connaitre le numéro de la ligne du tableau  			
        function getRow(inRowIndex){
            return ' ' + inRowIndex;
        }
 


 
 function tutu(value){
 	return value;
 }
 
 function setLayoutVar(type,grid) {
 	var url =url_layout_vector+"?type="+type+"&project="+projectName;
	
 	 dojo.xhrGet({ // ➂
        // The following URL must match that used to test the server.
        url: url,
        handleAs: "json",
        timeout: 100000, // Time in milliseconds
        // The LOAD function will be called on a successful response.
        load: function(responseObject, ioArgs){
            // Now you can just use the object
			//pour les projets de séquencage classique on trace des chromatogrammes
			
			for (var i=0 ; i<responseObject["items"].length;i++){
					if (responseObject["items"][i]["get"]) {
						responseObject["items"][i]["get"] = getRow;
					}
				if(responseObject["items"][i]["formatter"]) {
					var val = responseObject["items"][i]["formatter"] ;
					
					responseObject["items"][i]["formatter"] = formatters[val];
				}
			}
			
		
			grid.setStructure(responseObject["items"]);
			return;		
		},
		
        // The ERROR function will be called in an error case.
        error: function(response, ioArgs){ // ➃
    
			alert("error in layout "+type);
			console.dir(ioArgs);
			console.error(url_chart+url_arg); // ➆
            console.error("HTTP status code:  _--> ", ioArgs.xhr.status); // ➆
            return response; // ➅
        }
    });
	
 }
 
 function setLayoutVar2 (type) {

 	if (type == "ngs") {
		
		 layoutVariation = [
			// view 0
		
		
				{ name: "Row", get: getRow, width:2, styles: 'text-align:center;'},	
				{ name: "Name",field: "name",  formatter: formatSnp,width: 12, styles: 'text-align:center;'},
				{ name: "hp",field: "hapmap",  formatter: formatHapmap,width: 1, styles: 'text-align:center;'},
				{ name: "Nomenclature",field: "nomenclature",  width: 10, styles: 'text-align:center;'},
				{ name: "Allele",field: "text", width: 3 , styles: 'text-align:center;'},
				{ name: "Position",field: "start", width: 8 , styles: 'text-align:center;'},	
				{ name: "Patients",field: "lectures", width: 3 , styles: 'text-align:center;'},
				{ name: "maq",field: "maq",formatter: formatMethods,width: 6 , styles: 'text-align:center;'},
				{ name: "bowtie",field: "bowtie",formatter: formatMethods,width: 6 , styles: 'text-align:center;'},
				{ name: "Consequence",field: "consequence!all",styles: 'text-align:center;'},
				{ name: "Gene",field: "genes", formatter: formatReferenceGene,width: 10 , styles: 'text-align:center;'},
				{ name: "Filter",field: "filter", formatter: formatFilter,width: 5 , styles: 'text-align:center;'},
				{ name: "Valid",field: "valid", formatter: formatValid, width: 3 , styles: 'text-align:center;'},
				{name: "polyphen<br>status",field: "polyphen_status", formatter: formatPolyphen, width: 5,styles:'text-align:center;'} 
		
		];
		
		
	}
	else if (type == "array"){
		 layoutVariation = [
			// view 0
		
				{ name: "Row", get: getRow, width:2, styles: 'text-align:center;'},	
				{ name: "Name",field: "name",  width: 12, styles: 'text-align:center;'},
				{ name: "Nomenclature",field: "nomenclature",  width: 10, styles: 'text-align:center;'},
				{ name: "Allele",field: "text", width: 3 , styles: 'text-align:center;'},
				{ name: "Position",field: "start", width: 8 , styles: 'text-align:center;'},	
				{ name: "Nb traces",field: "lectures", width: 3 , styles: 'text-align:center;'},
				{ name: "Gseq",field: "gseq", width: 6 ,formatter: formatMethods, styles: 'text-align:center;'},
				{ name: "Consequence",field: "consequence",styles: 'text-align:center;'},
				{ name: "Gene",field: "genes", formatter: formatReferenceGene,width: 10 , styles: 'text-align:center;'},
				{ name: "Filter",field: "filter", formatter: formatFilter,width: 5 , styles: 'text-align:center;'},
				{ name: "Valid",field: "valid", formatter: formatValid, width: 3 , styles: 'text-align:center;'},
				{name: "polyphen<br>status",field: "polyphen_status", formatter: formatPolyphen, width: 5,styles:'text-align:center;'} 
			
		];
		
	} 
	else {

		 layoutVariation = [
			// view 0	
				{ name: "Row", get: getRow, width:2, styles: 'text-align:center;'},	
				{ name: "Name",field: "name",  formatter: formatSnp,width: 12, styles: 'text-align:center;'},
				{ name: "Nomenclature",field: "nomenclature",  width: 10, styles: 'text-align:center;'},
				{ name: "Allele",field: "text", width: 3 , styles: 'text-align:center;'},
				{ name: "Position",field: "start", width: 8 , styles: 'text-align:center;'},
				{ name: "contig",field: "group", width: 5 , styles: 'text-align:center;'},	
				{ name: "Nb traces",field: "lectures", width: 3 , styles: 'text-align:center;'},
				{ name: "Bipd",field: "bipd",formatter: formatMethods,width: 6 , styles: 'text-align:center;'},
				{ name: "Polyphred",field: "polyphred", formatter: formatMethods,width: 6 , styles: 'text-align:center;'},
				{ name: "Consequence",field: "consequence!all",styles: 'text-align:center;'},
				{ name: "Gene",field: "genes", formatter: formatReferenceGene,width: 10 , styles: 'text-align:center;'},
				{ name: "Filter",field: "filter", formatter: formatFilter,width: 5 , styles: 'text-align:center;'},
				{ name: "Valid",field: "valid", formatter: formatValid, width: 3 , styles: 'text-align:center;'},
				{name: "polyphen<br>status",field: "polyphen_status", formatter: formatPolyphen, width: 5,styles:'text-align:center;'} 
			
		];
		
	}
 	
	gridVar.setStructure(layoutVariation);
	
 }
 
         //layout pour les variations       		
	
		
		var vardesign = [
				{ name: "Row", get: getRow, width:2, styles: 'text-align:center;'},	
				{ name: "Name",field: "name",  width: 12, styles: 'text-align:center;'},
				{ name: "Allele",field: "text", width: 3 , styles: 'text-align:center;'},
				{ name: "Position",field: "start", width: 8 , styles: 'text-align:center;'},	
				{ name: "Nb traces",field: "lectures", width: 3 , styles: 'text-align:center;'},
				{ name: "Bipd2",field: "bipd",formatter: formatMethods,width: 6 , styles: 'text-align:center;'},
				{ name: "Polyphred",field: "polyphred", formatter: formatMethods,width: 6 , styles: 'text-align:center;'},
				{ name: "Gseq",field: "gseq", width: 6 ,formatter: formatMethods, styles: 'text-align:center;'},
				{ name: "Consequence",field: "consequence",styles: 'text-align:center;'},
				{ name: "Group",field: "group", width: 7 , styles: 'text-align:center;'},
				{ name: "Reference",field: "reference", formatter: formatReferenceGene,width: 10 , styles: 'text-align:center;'},
				{ name: "Filter",field: "filterStar", formatter: formatFilter,width: 5 , styles: 'text-align:center;'},
				{ name: "Valid",field: "valid", formatter: formatValid, width: 3 , styles: 'text-align:center;'},
				{name: "polyphen<br>status",field: "polyphen_status", width: 5,styles:'text-align:center;'} 
		];

			var layoutListGenes = [[
				{
					name: " ",
					field: "include",
					formatter: formatInclude,
					width: "15px",
					styles: 'text-align:center;',
					rowSpan: 2,
					hidden: true,
					noresize:true,
				}, {
					name: "ENSG",
					field: "name",
					formatter: formatReferenceGene,
					width: "7%",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "Gene Name",
					field: "xref",
					width: "7%",
					styles: 'text-align:center;font-size : 10px;bold;font-weight: bold;',
					rowSpan: 2,
					noresize:true,
				},{
					name: "Chr",
					field: "chromosome" ,
					width: "20px",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "Start",
					field: "start",
					width: "5%",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "End",
					field: "end",
					width: "5%",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "Phenotype",
					field: "phenotype",
					width: "24%",
					formatter: format_phenotype,
					styles: 'text-align:center;font-size : 08px;',
					rowSpan: 1,
					noresize:true,
				}, 
				{
					name: "PolyDiag Capture",
					field: "dejavu_capture_diag",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: dejavu_capture_diag,
					rowSpan: 2,
					noresize:true,
					
				}, 
				{
					name: "Omim",
					field: "is_omim",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: format_is_omim,
					rowSpan: 1,
					noresize:true,
					
				}, 
				{
					name: "test",
					field: "gene_atlas_diseases",
					width: "4%",
					formatter: formatFilter,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					hidden: true,
					noresize:true,
				}, 
				
				{
					name: "Nb Pat",
					field: "p_all",
					width: "4%",
					formatter: write_bold_lightgrey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
					},
				{
					name: "All",
					field: "v_all",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, 
				{
					name: "Subs",
					field: "v_substitution",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, 
				{
					name: "Ins/Del",
					field: "v_indel",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "Cnv",
					field: "v_cnv",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "Low",
					field: "v_low",
					width: "4%",
					formatter: write_bold_green,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "Medium",
					field: "v_medium",
					width: "4%",
					formatter: write_bold_orange,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "High",
					field: "v_high",
					width: "4%",
					formatter: write_bold_red,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "Poly Score",
					field: "max_scaled_score",
					width: "4%",
					formatter: formatGeneMaxScoreVariant,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					hidden: true,
					noresize:true,
				},
				], [{
					name: "Description",
					field: "description",
					width: "24%",
					styles: 'text-align:center;font-size : 08px;',
					rowSpan: 1,
					noresize:true,
				}, 
				{
					name: "Is Morbid",
					field: "is_omim_morbid",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: format_is_omim_morbid,
					rowSpan: 1,
					noresize:true,
					
				},
				]
			];

			var layoutListGenesHgmd = [[
				{
					name: " ",
					field: "include",
					formatter: formatInclude,
					width: "15px",
					styles: 'text-align:center;',
					rowSpan: 2,
					hidden: true,
					noresize:true,
				}, {
					name: "ENSG",
					field: "name",
					formatter: formatReferenceGene,
					width: "7%",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "Gene Name",
					field: "xref",
					width: "7%",
					styles: 'text-align:center;font-size : 10px;bold;font-weight: bold;',
					rowSpan: 2,
					noresize:true,
				},{
					name: "Chr",
					field: "chromosome" ,
					width: "20px",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2
				}, {
					name: "Start",
					field: "start",
					width: "5%",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "End",
					field: "end",
					width: "5%",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "Phenotype",
					field: "phenotype",
					width: "24%",
					formatter: format_phenotype,
					styles: 'text-align:center;font-size : 08px;',
					rowSpan: 1,
					noresize:true,
				}, 
				{
					name: "PolyDiag Capture",
					field: "dejavu_capture_diag",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: dejavu_capture_diag,
					rowSpan: 2,
					noresize:true,
					
				}, 
				{
					name: "Omim",
					field: "is_omim",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: format_is_omim,
					rowSpan: 1,
					noresize:true,
					
				}, 
				{
					name: "Gene HGMD DM",
					field: "has_hgmd",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: has_hgmd,
					rowSpan: 2,
					noresize:true,
				}, 
				{
					name: "test",
					field: "gene_atlas_diseases",
					width: "4%",
					formatter: formatFilter,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					hidden: true,
					noresize:true,
				}, 
				
				{
					name: "Nb Pat",
					field: "p_all",
					width: "4%",
					formatter: write_bold_lightgrey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
					}, 
				{
					name: "All",
					field: "v_all",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, 
				{
					name: "Subs",
					field: "v_substitution",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, 
				{
					name: "Ins/Del",
					field: "v_indel",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "Cnv",
					field: "v_cnv",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "Low",
					field: "v_low",
					width: "4%",
					formatter: write_bold_green,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "Medium",
					field: "v_medium",
					width: "4%",
					formatter: write_bold_orange,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "High",
					field: "v_high",
					width: "4%",
					formatter: write_bold_red,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "HGMD",
					field: "v_hgmd_dm",
					width: "4%",
					formatter: write_bold_purple,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "Poly Score",
					field: "max_scaled_score",
					width: "4%",
					formatter: formatGeneMaxScoreVariant,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					hidden: true,
					noresize:true,
				},
				], [  {
					name: "Description",
					field: "description",
					width: "24%",
					styles: 'text-align:center;font-size : 08px;',
					rowSpan: 1,
					noresize:true,
				}, 
				{
					name: "Is Morbid",
					field: "is_omim_morbid",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: format_is_omim_morbid,
					rowSpan: 1,
					noresize:true,
					
				},
				]
			];
			
			
			
			var layoutListGenesWithRegionHo = [[
				{
					name: " ",
					field: "include",
					formatter: formatInclude,
					width: "15px",
					styles: 'text-align:center;',
					rowSpan: 2,
					hidden: true,
					noresize:true,
				}, {
					name: "ENSG",
					field: "name",
					formatter: formatReferenceGene,
					width: "7%",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "Gene Name",
					field: "xref",
					width: "7%",
					styles: 'text-align:center;font-size : 10px;bold;font-weight: bold;',
					rowSpan: 2,
					noresize:true,
				},{
					name: "Chr",
					field: "chromosome" ,
					width: "20px",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "Start",
					field: "start",
					width: "5%",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "End",
					field: "end",
					width: "5%",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "Phenotype",
					field: "phenotype",
					width: "24%",
					formatter: format_phenotype,
					styles: 'text-align:center;font-size : 08px;',
					rowSpan: 1,
					noresize:true,
				}, 
				{
					name: "PolyDiag Capture",
					field: "dejavu_capture_diag",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: dejavu_capture_diag,
					rowSpan: 2,
					noresize:true,
					
				}, 
				{
					name: "Omim",
					field: "is_omim",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: format_is_omim,
					rowSpan: 1,
					noresize:true,
					
				}, 
				{
					name: "test",
					field: "gene_atlas_diseases",
					width: "4%",
					formatter: formatFilter,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					hidden: true,
					noresize:true,
				},  
				{
					name: "Nb Pat HoReg",
					field: "region_ho_all",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: formatRegionHo,
					rowSpan: 2,
					noresize:true,
				},  
				
				{
					name: "Nb Pat",
					field: "p_all",
					width: "4%",
					formatter: write_bold_lightgrey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
					}, 
				{
					name: "All",
					field: "v_all",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, 
				{
					name: "Subs",
					field: "v_substitution",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2
				}, 
				{
					name: "Ins/Del",
					field: "v_indel",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "Cnv",
					field: "v_cnv",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "Low",
					field: "v_low",
					width: "4%",
					formatter: write_bold_green,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "Medium",
					field: "v_medium",
					width: "4%",
					formatter: write_bold_orange,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "High",
					field: "v_high",
					width: "4%",
					formatter: write_bold_red,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "Poly Score",
					field: "max_scaled_score",
					width: "4%",
					formatter: formatGeneMaxScoreVariant,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					hidden: true,
					noresize:true,
				},
				], [{
					name: "Description",
					field: "description",
					width: "24%",
					styles: 'text-align:center;font-size : 08px;',
					rowSpan: 1,
					noresize:true,
				}, 
				{
					name: "Is Morbid",
					field: "is_omim_morbid",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: format_is_omim_morbid,
					rowSpan: 1,
					
				},
				]
			];
			
			var layoutListGenesWithRegionRec = [[
				{
					name: " ",
					field: "include",
					formatter: formatInclude,
					width: "15px",
					styles: 'text-align:center;',
					rowSpan: 2,
					hidden: true,
					noresize:true,
				}, {
					name: "ENSG",
					field: "name",
					formatter: formatReferenceGene,
					width: "7%",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "Gene Name",
					field: "xref",
					width: "7%",
					styles: 'text-align:center;font-size : 10px;bold;font-weight: bold;',
					rowSpan: 2,
					noresize:true,
				},{
					name: "Chr",
					field: "chromosome" ,
					width: "20px",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "Start",
					field: "start",
					width: "5%",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "End",
					field: "end",
					width: "5%",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "Phenotype",
					field: "phenotype",
					width: "24%",
					formatter: format_phenotype,
					styles: 'text-align:center;font-size : 08px;',
					rowSpan: 1,
					noresize:true,
				}, 
				{
					name: "PolyDiag Capture",
					field: "dejavu_capture_diag",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: dejavu_capture_diag,
					rowSpan: 2,
					noresize:true,
					
				}, 
				{
					name: "Omim",
					field: "is_omim",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: format_is_omim,
					rowSpan: 1,
					noresize:true,
					
				}, 
				{
					name: "test",
					field: "gene_atlas_diseases",
					width: "4%",
					formatter: formatFilter,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					hidden: true,
					noresize:true,
				},   
				{
					name: "Fam RecReg",
					field: "region_rec_all",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: formatRegionRec,
					rowSpan: 2,
					noresize:true,
					
				}, 
				
				{
					name: "Nb Pat",
					field: "p_all",
					width: "4%",
					formatter: write_bold_lightgrey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
					},
				{
					name: "All",
					field: "v_all",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, 
				{
					name: "Subs",
					field: "v_substitution",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, 
				{
					name: "Ins/Del",
					field: "v_indel",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "Cnv",
					field: "v_cnv",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "Low",
					field: "v_low",
					width: "4%",
					formatter: write_bold_green,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "Medium",
					field: "v_medium",
					width: "4%",
					formatter: write_bold_orange,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
				},
				{
					name: "High",
					field: "v_high",
					width: "4%",
					formatter: write_bold_red,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				},
				{
					name: "Poly Score",
					field: "max_scaled_score",
					width: "4%",
					formatter: formatGeneMaxScoreVariant,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					hidden: true,
					noresize:true,
				},
				], [{
					name: "Description",
					field: "description",
					width: "24%",
					styles: 'text-align:center;font-size : 08px;',
					rowSpan: 1,
					noresize:true,
				}, 
				{
					name: "Is Morbid",
					field: "is_omim_morbid",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: format_is_omim_morbid,
					rowSpan: 1,
					noresize:true,
					
				},
				]
			
			];
			
			
			
			var layoutListGenesWithRegionHoHgmd = [[
				{
					name: " ",
					field: "include",
					formatter: formatInclude,
					width: "15px",
					styles: 'text-align:center;',
					rowSpan: 2,
					hidden: true,
					noresize:true,
				}, {
					name: "ENSG",
					field: "name",
					formatter: formatReferenceGene,
					width: "7%",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "Gene Name",
					field: "xref",
					width: "7%",
					styles: 'text-align:center;font-size : 10px;bold;font-weight: bold;',
					rowSpan: 2,
					noresize:true,
				},{
					name: "Chr",
					field: "chromosome" ,
					width: "20px",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "Start",
					field: "start",
					width: "5%",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "End",
					field: "end",
					width: "5%",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "Phenotype",
					field: "phenotype",
					width: "24%",
					formatter: format_phenotype,
					styles: 'text-align:center;font-size : 08px;',
					rowSpan: 1,
					noresize:true,
				}, 
				{
					name: "PolyDiag Capture",
					field: "dejavu_capture_diag",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: dejavu_capture_diag,
					rowSpan: 2,
					noresize:true,
					
				}, 
				{
					name: "Omim",
					field: "is_omim",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: format_is_omim,
					rowSpan: 1,
					noresize:true,
					
				}, 
				{
					name: "Gene HGMD DM",
					field: "has_hgmd",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: has_hgmd,
					rowSpan: 2,
					noresize:true,
				}, 
				{
					name: "test",
					field: "gene_atlas_diseases",
					width: "4%",
					formatter: formatFilter,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					hidden: true,
					noresize:true,
				},  
				{
					name: "Nb Pat HoReg",
					field: "region_ho_all",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: formatRegionHo,
					rowSpan: 2,
					noresize:true,
				},  
				
			/*	{
					name: "cover",
					field: "coding_cover",
					width: 5,
					styles: 'text-align:center;',
					rowSpan: 2,
					noresize:true,
				},*/ 
				{
					name: "All",
					field: "v_all",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
				}, 
				{
					name: "Subs",
					field: "v_substitution",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
				}, 
				{
					name: "Ins/Del",
					field: "v_indel",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
				},
				{
					name: "Cnv",
					field: "v_cnv",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
				},
				{
					name: "Low",
					field: "v_low",
					width: "4%",
					formatter: write_bold_green,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
				},
				{
					name: "Medium",
					field: "v_medium",
					width: "4%",
					formatter: write_bold_orange,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
				},
				{
					name: "High",
					field: "v_high",
					width: "4%",
					formatter: write_bold_red,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
				},
				{
					name: "HGMD",
					field: "v_hgmd_dm",
					width: "4%",
					formatter: write_bold_purple,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
				},
				{
					name: "Poly Score",
					field: "max_scaled_score",
					width: "4%",
					formatter: formatGeneMaxScoreVariant,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					hidden: true,
					noresize:true,
				},
				], [{
					name: "Description",
					field: "description",
					width: "24%",
					styles: 'text-align:center;font-size : 08px;',
					rowSpan: 1,
					noresize:true,
				}, 
				{
					name: "Is Morbid",
					field: "is_omim_morbid",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: format_is_omim_morbid,
					rowSpan: 1,
					noresize:true,
					
				},{
					name: "Nb Pat",
					field: "p_all",
					width: "4%",
					formatter: write_bold_lightgrey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
					},
					{
					name: "Nb Pat",
					field: "p_substitution",
					width: "4%",
					formatter: write_bold_lightgrey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
					},
					{
					name: "Nb Pat",
					field: "p_indel",
					width: "4%",
					formatter: write_bold_lightgrey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
					},
					{
					name: "Nb Pat",
					field: "p_cnv",
					width: "4%",
					formatter: write_bold_lightgrey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
					},
				{
					name: "Nb Pat",
					field: "p_low",
					width: "4%",
					formatter: write_bold_lightgreen,
					styles: 'text-align:center;font-size : 09px;',
					colSpan: 1,
					noresize:true,
				},
				{
					name: "Nb Pat",
					field: "p_medium",
					width: "4%",
					formatter: write_bold_lightorange,
					styles: 'text-align:center;font-size : 09px;',
					colSpan: 1,
					noresize:true,
				},
					{
					name: "Nb Pat",
					field: "p_high",
					width: "4%",
					formatter: write_bold_lightred,
					styles: 'text-align:center;font-size : 09px;',
					colSpan: 1,
					noresize:true,
				},
				{
					name: "Nb Pat",
					field: "p_hgmd_dm",
					width: "4%",
					formatter: write_bold_lightpurple,
					styles: 'text-align:center;font-size : 09px;',
					colSpan: 1,
					noresize:true,
				},
				]
			];
			
			var layoutListGenesWithRegionRecHgmd = [[
				{
					name: " ",
					field: "include",
					formatter: formatInclude,
					width: "15px",
					styles: 'text-align:center;',
					rowSpan: 2,
					hidden: true,
					noresize:true,
				}, {
					name: "ENSG",
					field: "name",
					formatter: formatReferenceGene,
					width: "7%",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "Gene Name",
					field: "xref",
					width: "7%",
					styles: 'text-align:center;font-size : 10px;bold;font-weight: bold;',
					rowSpan: 2,
					noresize:true,
				},{
					name: "Chr",
					field: "chromosome" ,
					width: "20px",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "Start",
					field: "start",
					width: "5%",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "End",
					field: "end",
					width: "5%",
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					noresize:true,
				}, {
					name: "Phenotype",
					field: "phenotype",
					width: "24%",
					formatter: format_phenotype,
					styles: 'text-align:center;font-size : 08px;',
					rowSpan: 1,
					noresize:true,
				}, 
				{
					name: "PolyDiag Capture",
					field: "dejavu_capture_diag",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: dejavu_capture_diag,
					rowSpan: 2,
					noresize:true,
					
				}, 
				{
					name: "Omim",
					field: "is_omim",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: format_is_omim,
					rowSpan: 1,
					noresize:true,
					
				}, 
				{
					name: "Gene HGMD DM",
					field: "has_hgmd",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: has_hgmd,
					rowSpan: 2,
					noresize:true,
				}, 
				{
					name: "test",
					field: "gene_atlas_diseases",
					width: "4%",
					formatter: formatFilter,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					hidden: true,
					noresize:true,
				},   
				{
					name: "Fam RecReg",
					field: "region_rec_all",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: formatRegionRec,
					rowSpan: 2,
					noresize:true,
					
				}, 
				
			/*	{
					name: "cover",
					field: "coding_cover",
					width: 5,
					styles: 'text-align:center;',
					rowSpan: 2,
					noresize:true,
				},*/ 
				{
					name: "All",
					field: "v_all",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
				}, 
				{
					name: "Subs",
					field: "v_substitution",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
				}, 
				{
					name: "Ins/Del",
					field: "v_indel",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
				},
				{
					name: "Cnv",
					field: "v_cnv",
					width: "4%",
					formatter: write_bold_grey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
				},
				{
					name: "Low",
					field: "v_low",
					width: "4%",
					formatter: write_bold_green,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
				},
				{
					name: "Medium",
					field: "v_medium",
					width: "4%",
					formatter: write_bold_orange,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
				},
				{
					name: "High",
					field: "v_high",
					width: "4%",
					formatter: write_bold_red,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
				},
				{
					name: "HGMD",
					field: "v_hgmd_dm",
					width: "4%",
					formatter: write_bold_purple,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
				},
				{
					name: "Poly Score",
					field: "max_scaled_score",
					width: "4%",
					formatter: formatGeneMaxScoreVariant,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 2,
					hidden: true,
					noresize:true,
				},
				], [{
					name: "Description",
					field: "description",
					width: "24%",
					styles: 'text-align:center;font-size : 08px;',
					rowSpan: 1,
					noresize:true,
				}, 
				{
					name: "Is Morbid",
					field: "is_omim_morbid",
					width: "4%",
					styles: 'text-align:center;font-size : 09px;',
					formatter: format_is_omim_morbid,
					rowSpan: 1,
					noresize:true,
					
				},{
					name: "Nb Pat",
					field: "p_all",
					width: "4%",
					formatter: write_bold_lightgrey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
					},
					{
					name: "Nb Pat",
					field: "p_substitution",
					width: "4%",
					formatter: write_bold_lightgrey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
					},
					{
					name: "Nb Pat",
					field: "p_indel",
					width: "4%",
					formatter: write_bold_lightgrey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
					},
					{
					name: "Nb Pat",
					field: "p_cnv",
					width: "4%",
					formatter: write_bold_lightgrey,
					styles: 'text-align:center;font-size : 09px;',
					rowSpan: 1,
					noresize:true,
					},
				{
					name: "Nb Pat",
					field: "p_low",
					width: "4%",
					formatter: write_bold_lightgreen,
					styles: 'text-align:center;font-size : 09px;',
					colSpan: 1,
					noresize:true,
				},
				{
					name: "Nb Pat",
					field: "p_medium",
					width: "4%",
					formatter: write_bold_lightorange,
					styles: 'text-align:center;font-size : 09px;',
					colSpan: 1,
					noresize:true,
				},
					{
					name: "Nb Pat",
					field: "p_high",
					width: "4%",
					formatter: write_bold_lightred,
					styles: 'text-align:center;font-size : 09px;',
					colSpan: 1,
					noresize:true,
				},
				{
					name: "Nb Pat",
					field: "p_hgmd_dm",
					width: "4%",
					formatter: write_bold_lightpurple,
					styles: 'text-align:center;font-size : 09px;',
					colSpan: 1,
					noresize:true,
				},
				]
			
			];