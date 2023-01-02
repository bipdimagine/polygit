var layout_analysis_step = [{
                    field: "patient",
                    name: "Patient",
                    width: '8',
                    styles: 'text-align:left;'
                },
				{
                    field: "analysis_id",
                    name: "Analysis id",
                    width: '6',
                    styles: 'text-align:center;'
                },
                {
                    field: "step_name",
                    name: "Step name",
                    width: '8',
                    styles: 'text-align:left;'
                },
                {
                    field: "step_order",
                    name: "Step order",
                    width: '6',
                    styles: 'text-align:center;'
                },
                
                {
                    field: "status",
                    name: "Status",
                    width: '6',
                    styles: 'text-align:left;'
                    
                },
                {
                    field: "start",
                    name: "Start",
                    width: '10',
                    styles: 'text-align:left;'
                },
                {
                    field: "end",
                    name: "End",
                    width: '12',
                    styles: 'text-align:left;'
                },
                {
                    field: "command",
                    name: "Command",
                    width: 'auto',
                    styles: 'text-align:left;'
                }];
                
var layout_analysis_version = [{
                    field: "analysis_id",
                    name: "Analysis id",
                    width: '8',
                    styles: 'text-align:center;'
                },
                {
                    field: "version_name",
                    name: "Version name",
                    width: '10',
                    styles: 'text-align:left;'
                },
                {
                    field: "version_path",
                    name: "Version info",
                    width: 'auto',
                    styles: 'text-align:left;'
                }];
                
var layout_analysis_program = [{
                    field: "analysis_id",
                    name: "Analysis id",
                    width: '8',
                    styles: 'text-align:center;'
                },
                {
                    field: "program_name",
                    name: "Program name",
                    width: '10',
                    styles: 'text-align:left;'
                },
                {
                    field: "program_id",
                    name: "Program id",
                    width: '8',
                    styles: 'text-align:center;'
                },
                {
                    field: "program_path",
                    name: "Program path",
                    width: 'auto',
                    styles: 'text-align:left;'
                },
                {
                    field: "program_md5",
                    name: "Program md5",
                    width: '20',
                    styles: 'text-align:left;'
                }];
                
var layout_analysis_reference = [{
                    field: "analysis_id",
                    name: "Analysis id",
                    width: '8',
                    styles: 'text-align:center;'
                },
                {
                    field: "ref_name",
                    name: "Reference name",
                    width: '10',
                    styles: 'text-align:left;'
                },
                {
                    field: "ref_id",
                    name: "Reference id",
                    width: '8',
                    styles: 'text-align:center;'
                },
                {
                    field: "ref_path",
                    name: "Reference path",
                    width: 'auto',
                    styles: 'text-align:left;'
                },
                {
                    field: "ref_md5",
                    name: "Reference md5",
                    width: '20',
                    styles: 'text-align:left;'
                }];
   
                

dojo.addOnLoad(function(){
	
	
    jsonStoreS = new dojo.data.ItemFileReadStore({
        url: url_polypipeline_analysis_step + "?analysis_id="+analysis_id,
        handleAs:"json",
    });
    
	
     var GridAnalysisStep = new dojox.grid.DataGrid({
      	query:{'analysis_id':analysis_id},
        store: jsonStoreS,
        structure: layout_analysis_step ,
    },document.createElement('div'));
    

    dojo.byId("gridAnalysisStep_div").appendChild(GridAnalysisStep.domNode);
    GridAnalysisStep.startup();
    
    
    jsonStoreV = new dojo.data.ItemFileReadStore({
        url: url_polypipeline_analysis_version + "?analysis_id="+analysis_id,
        handleAs:"json",
    });
    
	
     var GridAnalysisVersion = new dojox.grid.DataGrid({
      	query:{'analysis_id':analysis_id},
        store: jsonStoreV,
        structure: layout_analysis_version ,
    },document.createElement('div'));
    

    dojo.byId("gridAnalysisVersion_div").appendChild(GridAnalysisVersion.domNode);
    GridAnalysisVersion.startup();
    
    
    jsonStoreP = new dojo.data.ItemFileReadStore({
        url: url_polypipeline_analysis_program + "?analysis_id="+analysis_id,
        handleAs:"json",
    });
    
	
     var GridAnalysisProgram = new dojox.grid.DataGrid({
      	query:{'analysis_id':analysis_id},
        store: jsonStoreP,
        structure: layout_analysis_program ,
    },document.createElement('div'));
    

    dojo.byId("gridAnalysisProgram_div").appendChild(GridAnalysisProgram.domNode);
    GridAnalysisProgram.startup();
    
    
    jsonStoreR = new dojo.data.ItemFileReadStore({
        url: url_polypipeline_analysis_reference + "?analysis_id="+analysis_id,
        handleAs:"json",
    });
    
	
     var GridAnalysisRef = new dojox.grid.DataGrid({
      	query:{'analysis_id':analysis_id},
        store: jsonStoreR,
        structure: layout_analysis_reference ,
    },document.createElement('div'));
    

    dojo.byId("gridAnalysisRef_div").appendChild(GridAnalysisRef.domNode);
    GridAnalysisRef.startup();
    
 });  
 
 
 function colorCov(value,idx,level,cell,row) {
    if (cell.field == "status" & value=="complete") {
        cell.customStyles.push("background: green");
    } else if (cell.field == "status" & value=="waiting") {
        cell.customStyles.push("background: yellow");   
    } else if (cell.field == "status" & value=="started") {
        cell.customStyles.push("background: blue");     
    } else if (cell.field == "status" & value=="error") {
        cell.customStyles.push("background: red");      
    } 
    return value;
} 