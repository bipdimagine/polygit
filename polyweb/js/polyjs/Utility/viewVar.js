/**
 * @author pnitschk
 */

	
viewVarByGene = function(e){
		var i = e.rowNode.gridRowIndex;	
		var projectName = param("project");//parent.haut.document.getElementById('project_name').value;
		var name=gridListGenes.getCell(1).getNode(i).textContent;	
			var store =gridListGenes.store;
			store.fetch({
				query: {
					name: name
				},
				onComplete: goViewVar
			});
		//principal(projectName,refname,genename);
	}
	
function goViewVar(items, request){
	var murl = "detailProject.html?project="+projectName+"&chromosome="+items[0].chromosome+"&gene="+items[0].name+"&genename="+items[0].xref+"&type="+project_type+"&ids="+items[0].ids+"&model="+genetic_model+"&famind="+famind;
	window.open(murl);
	return;
}
 
 function goViewAllVariations(){
	dijit.byId('waiting').show();
 	document.location.href = "detailProject.html?project="+projectName+"&gene=all"+"&type="+project_type+"&model="+genetic_model+"&famind="+famind;
}
 