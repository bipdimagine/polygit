/**
 * @author pnitschk
 */
 
dojo.require("dojo.cookie");
var url_passwd;
var error_vcf_md5 = 0;

	 function checkPassword(project,okf){
	 		   
               var username = dojo.cookie("username");
			  var passwd = dojo.cookie("passwd");
			
			  urlpasswd = url_listProject_vector+"?type=login&pwd="+passwd+"&login="+username+"&action=list";
		

			  if (project){
			  	urlpasswd = urlpasswd +"&project="+project;
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

function checkMD5VcfFiles(username, passwd, projectName) {
	var url_md5 = url_listProject_vector+"?login="+username+"&pwd="+passwd+"&project="+projectName+"&action=checkMD5";
	var jsonStore = new dojo.data.ItemFileWriteStore({ url: url_md5	});
  	jsonStore.fetch({
		onComplete: okfunction_md5,
		onError: error,		
	});	
}

function okfunction_md5 (items, request){
	if (items.length > 0) {
		error_vcf_md5 = 1;
		alert ("ERROR: please contact bioinformatics platform !! ERROR with VCF MD5 check for "+projectName+" project !");
	}
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
				dojo.cookie("username", username, { expires: 1 });
				dojo.cookie("passwd", passwd, { expires: 1 });
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
