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
 	var url =url_layout+"?type="+type+"&project="+projectName;
	
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
					name: "Row",
					get: getRow,
					width: 4,
					styles: 'text-align:center;',
					rowSpan: 2
				}, {
					name: "Name",
					field: "name",
					formatter: formatReferenceGene,
					width: 13,
					styles: 'text-align:center;',
					rowSpan: 2
				}, {
					name: "xref",
					field: "xref",
					width: 9,
					styles: 'text-align:center;',
					rowSpan: 2
				},{
					name: "chr",
					field: "chromosome" ,
					width: 2,
					styles: 'text-align:center;',
					rowSpan: 2
				}, {
					name: "Start",
					field: "start",
					width: 7,
					styles: 'text-align:center;',
					rowSpan: 2
				}, {
					name: "End",
					field: "end",
					width: 7,
					styles: 'text-align:center;',
					rowSpan: 2
				}, {
					name: "Description",
					field: "description",
					width:30,
					styles: 'text-align:center;font-size : 08px;',
					rowSpan: 2
				}, 
				{
					name: "Diseases",
					field: "gene_atlas_diseases",
					width:5,
					formatter: formatFilter,
					styles: 'text-align:center;',
					rowSpan: 2,
					
				}, 
				{
					name: "test",
					field: "gene_atlas_diseases",
					width:5,
					formatter: formatFilter,
					styles: 'text-align:center;',
					rowSpan: 2,
					hidden: true
				}, 
				
			/*	{
					name: "cover",
					field: "coding_cover",
					width: 5,
					styles: 'text-align:center;',
					rowSpan: 2
				},*/ 
				{
					name: "All",
					field: "v_all",
					width: 3,
					styles: 'text-align:center;',
					rowSpan: 1
				}, 
				{
					name: "Subs",
					field: "v_substitution",
					width: 3,
					styles: 'text-align:center;',
					rowSpan: 1
				}, 
				{
					name: "Ins",
					field: "v_insertion",
					width: 3,
					styles: 'text-align:center;',
					rowSpan: 1
				},
				{
					name: "Dels",
					field: "v_deletion",
					width: 3,
					styles: 'text-align:center;',
					rowSpan: 1
				},
				{
					name: "With cons.",
					field: "coding_consequence",
					width: 5,
					styles: 'text-align:center;',
					rowSpan: 1
				}, {
					name: "Syno",
					field: "v_silent",
					width: 3,
					styles: 'text-align:center;',
					rowSpan: 1
				},
				{
					name: "UTR",
					field: "v_utr",
					width: 3,
					styles: 'text-align:center;',
					rowSpan: 1
				},
				{
					name: "Splicing",
					field: "v_splicing",
					width: 5,
					styles: 'text-align:center;',
					rowSpan: 1
				}
				], [{
					name: "Pat.",
					field: "p_all",
					width: 3,
					styles: 'text-align:center;',
					rowSpan: 1
					},
					{
					name: "Pat.",
					field: "p_substitution",
					width: 3,
					styles: 'text-align:center;',
					rowSpan: 1
					},
					{
					name: "Pat.",
					field: "p_insertion",
					width: 3,
					styles: 'text-align:center;',
					rowSpan: 1
					},
					{
					name: "Pat.",
					field: "p_deletion",
					width: 3,
					styles: 'text-align:center;',
					rowSpan: 1
					},
					{
					name: "Pat.",
					field: "p_coding_consequence",
					width: 3,
					styles: 'text-align:center;',
					colSpan: 1
				},
				{
					name: "Pat.",
					field: "p_silent",
					width: 3,
					styles: 'text-align:center;',
					colSpan: 1
				},
				{
					name: "Pat.",
					field: "p_utr",
					width: 3,
					styles: 'text-align:center;',
					colSpan: 1
				},
				{
					name: "Pat.",
					field: "p_splicing",
					width: 3,
					styles: 'text-align:center;',
					colSpan: 1
				},
				]
					
				
			
			];