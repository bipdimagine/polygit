var layout_versions = [{
                    field: "version_name",
                    name: "Version name",
                    width: '10',
                    styles: 'text-align:center;'
                },
                {
                    field: "version_id",
                    name: "Id",
                    width: '4',
                    styles: 'text-align:center;'
                },
                {
                    field: "version_path",
                    name: "Path",
                    width: 'auto',
                    styles: 'text-align:left;'
                }];
           
           
var layout_refs = [{
                    field: "ref_name",
                    name: "Reference name",
                    width: '20',
                    styles: 'text-align:center;'
                },
                {
                    field: "ref_id",
                    name: "Id",
                    width: '4',
                    styles: 'text-align:center;'
                },
                {
                    field: "ref_path",
                    name: "Path",
                    width: '45',
                    styles: 'text-align:left;'
                },
                {
                    field: "ref_md5",
                    name: "md5 key",
                    width: 'auto',
                    styles: 'text-align:left;'
                }];
var layout_programs = [{
                    field: "program_name",
                    name: "Program name",
                    width: '20',
                    styles: 'text-align:center;'
                },
                {
                    field: "program_id",
                    name: "Id",
                    width: '4',
                    styles: 'text-align:center;'
                },
                {
                    field: "program_path",
                    name: "Path",
                    width: '45',
                    styles: 'text-align:left;'
                },
                {
                    field: "program_md5",
                    name: "md5 key",
                    width: 'auto',
                    styles: 'text-align:left;'
                }];
                
var layout_projects = [{
                    field: "project_name",
                    name: "Project name",
                    width: '20',
                    styles: 'text-align:center;'
                },
                {
                    field: "status",
                    name: "Status",
                    width: '10',
                    styles: 'text-align:left;'
                },
                {
                    field: "analysis_id",
                    name: "Analysis Ids",
                    width: 'auto',
                    styles: 'text-align:center;'
                },
                ];

var layout_analysis = [{
                    field: "analysis_id",
                    name: "Analysis id",
                    width: '4',
                    styles: 'text-align:center;'
                },
                {
                    field: "project_name",
                    name: "Project name",
                    width: '10',
                    styles: 'text-align:left;'
                },
                {
                    field: "project_id",
                    name: "Project id",
                    width: '4',
                    styles: 'text-align:center;'
                },
                {
                    field: "patients_names",
                    name: "Patients names",
                    width: 'auto',
                    styles: 'text-align:left;'
                },
                {
                    field: "patients_ids",
                    name: "Patients ids",
                    width: '10',
                    styles: 'text-align:left;'
                },
                {
                    field: "user",
                    name: "User",
                    width: '10',
                    styles: 'text-align:left;'
                },
                {
                    field: "date",
                    name: "Date",
                    width: '12',
                    styles: 'text-align:left;'
                },
                {
                    field: "status",
                    name: "Status",
                    width: '6',
                    styles: 'text-align:left;'
                }
                ,
                {
                    field: "type",
                    name: "Type",
                    width: '5',
                    styles: 'text-align:left;'
                },
                {
                    field: "pipeline_steps",
                    name: "Pipeline steps",
                    width: 'auto',
                    styles: 'text-align:left;'
                }
                ];


dojo.addOnLoad(function(){
	
	
    jsonStore = new dojo.data.ItemFileReadStore({
        url: url_polypipeline_version,
        handleAs:"json",
    });

      var GridVersions = new dojox.grid.DataGrid({
      	query:{'version_name':'*'},
        store: jsonStore,
        structure: layout_versions ,
    },document.createElement('div'));
    

    dojo.byId("gridVersions_div").appendChild(GridVersions.domNode);
    GridVersions.startup();
    
    jsonStore2 = new dojo.data.ItemFileReadStore({
        url: url_polypipeline_ref,
        handleAs:"json",
    });
    
    var GridRefs = new dojox.grid.DataGrid({
      	query:{'ref_name':'*'},
        store: jsonStore2,
        structure: layout_refs ,
    },document.createElement('div'));
    

    dojo.byId("gridRefs_div").appendChild(GridRefs.domNode);
    GridRefs.startup();
    
    jsonStore3 = new dojo.data.ItemFileReadStore({
        url: url_polypipeline_prog,
        handleAs:"json",
    });
    
    var GridProgs = new dojox.grid.DataGrid({
      	query:{'program_name':'*'},
        store: jsonStore3,
        structure: layout_programs ,
    },document.createElement('div'));
    

    dojo.byId("gridProgs_div").appendChild(GridProgs.domNode);
    GridProgs.startup();
    
    jsonStore4 = new dojo.data.ItemFileReadStore({
        url: url_polypipeline_analysis,
        handleAs:"json",
    });
    
      jsonStore5 = new dojo.data.ItemFileReadStore({
        url: url_polypipeline_project,
        handleAs:"json",
    });
    var GridProjects = new dojox.grid.DataGrid({
      	query:{'project_name':'*'},
        store: jsonStore5,
        structure: layout_projects ,
    },document.createElement('div'));
    

    dojo.byId("gridProjects_div").appendChild(GridProjects.domNode);
    GridProjects.startup();
    
  
    
});
    
    
       
//loading et paramètres du tableau d'analyse
dojo.addOnLoad(function(){            
    
    var GridAnalysis = new dojox.grid.DataGrid({
      	query:{'analysis_id':'*'},
        store: jsonStore4,
        structure: layout_analysis ,
    },document.createElement('div'));
    
	
    dojo.byId("gridAnalysis_div").appendChild(GridAnalysis.domNode);
    
    
    //executer fonction dblclickGrid si double click sur tableau des projets
    
         
    GridAnalysis.startup();
    
    
	//fonction à executer si double click sur une ligne
	dblclickGrid = function(e){
      		var analysis_id = GridAnalysis.getCell(0).getNode(e.rowNode.gridRowIndex).textContent;
      		viewAnalysis(analysis_id);
            };
            
            function viewAnalysis(analysis_id){                
               this.location.href = "analysis.html?analysis_id="+analysis_id;
           };
    dojo.connect(GridAnalysis, 'onRowDblClick', dblclickGrid);
 });   




            
            

