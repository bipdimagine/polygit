/**
 * @author pnitschk
 */
 
dojo.require("dojo.cookie");
  var url_passwd;
	 function checkPassword(project,okf,viewer){
             var  username = dojo.cookie("username");
			  var passwd = dojo.cookie("passwd");
			
			  urlpasswd = url_listProject+"?type=login&pwd="+passwd+"&login="+username+"&action=list";
		

			  if (project){
			  	urlpasswd = urlpasswd +"&project="+project+"&viewer="+viewer;
			  }
			  var jsonStore = new dojo.data.ItemFileWriteStore({
						url: urlpasswd
				});
			  	jsonStore.fetch({
				query: {
					name: project
				},
				onComplete: okf,
				onError: error,		
			});	
				
			
 }
 
 function reloadGene() {
	var prj = filterProject.getValue();
	if (prj == projectName) return;
 	var url = "gene.html?project="+prj;
	this.location.href = url;	
 }
 function reloadVariations() {
 	var gene = filterSelectGenes.getValue();
	if (gene == xref) return;
	dijit.byId('dialogGenes').show();
	var url = "gene.html?project="+projectName+"&gene="+gene+"&reference="+reference;
	var geneurl = url_query+"?type=gene&projectName="+projectName;
	gridListGenes.setStore(readStore(geneurl));
	//this.location.href = url;	
 }
 

   function setPassword(dialogFields){
                var passwd = document.getElementById('passwd').value;
                var username = document.getElementById('username').value;	
                var midnight = new Date();
                midnight.setHours(23,59,59,0);

				dojo.cookie("username", username, { expires:midnight });
				dojo.cookie("passwd", passwd, { expires: midnight });
				dijit.byId('login').hide();
				checkPassword(projectName,okfunction);
				
  }	
  
  function error (items, request){
  	console.dir(items);
	console.dir(request);
	console.error("error in vp.js");
}

function relogin(){
    dijit.byId('login').show();
}
