/**
 * @author masson***
 */
require(["dojo/parser"]);



//fonction pour changer le statut de validation des objets
function changeStatus(type,grid,url,objet) {
	//recupere la/les lignes selectionnees
	
	var z= grid.selection.getSelected()[0];
	var sel = grid.selection.getSelected();
	
	
	//nombre de lignes selectionnees
	var l = sel.length;
	
	var allids='';
	
	//pour chaque lignes sellectionnees
	for (var i = 0; i < l; i++) {
		var z = sel[i];
		var gid;
	
		//si le parametre objet passe a la fonction est relation alors le type d'objet a valider est une relation entre une variation est un patient
		if(objet == "relation"){
			
			// on donne a gid la valeur du contenu de la colonne nommee relid pour la ligne selectionnae
			gid = sel[i]['relid'];
			
		}
		//sinon on donne a gid la valeur contenue dans la colonne id
		else{
			gid = sel[i]['id'];	
		}
		//on donne a mother_id la valeur contenue dans la colonne mother_id
		var mother_id=sel[i]['mother_id'];
		//on concatene tous les gid dans une liste en les separants par des points-virgules
		allids += ";"+gid;
	}
	
	//on recupere le formulaire appele statusForm dans la page html
	var fo1 = document.getElementById("statusForm");
	//on donne a l'elemnt type du formulaire la valeur type passee en argument a la fonction
	var el = document.getElementById("type");
	el.value = type;
	// on renseigne l'element objet du formulaire de la meme facon
	var el1 = document.getElementById("objet");
	el1.value = objet;
	// on renseigne l'element id du formulaire avec les gid concatenes precedemment
	var el2 = document.getElementById("id");
	el2.value = allids;
	// on renseigne l'element mother_id du formulaire avec la variable mother_id
	var el3 = document.getElementById("mother_id");
	el3.value = mother_id;
	//on envoit les donnees au formulaire
	sendData(url_insert, fo1);
	//on recharge le tableau
	forceRefreshTable(url,grid);
	
}



function updateGlobaleValidation(item){

	return;
	var nbvalid = 0;
	var nbunvalid = 0;
	var nbuncertain = 0;
	for (var i = 0; i < item.patient_name.length; i++) {
		 var tname= item.patient_name[i];
		 if (item["valid!"+tname] >=1) nbvalid ++; 
		  if (item["valid!"+tname] == -1) nbunvalid ++; 
		  if (item["valid!"+tname] == -2) nbuncertain ++; 
	}
	var rep = -8;
	if (nbvalid == item.patient_name.length) rep=1;
	else if (nbunvalid == item.patient_name.length) rep=-1;
	else if (nbuncertain == item.patient_name.length) rep=-2;
	else if ((nbuncertain+nbunvalid+nbvalid) == 0) rep=0;
	gridVar.store.setValue(item, "valid", rep);
	
} 


function updateValidation(value,username,gridsmall,gridbig,bigbutton) {
	//recupere la/les lignes selectionnees
	var varsel= gridbig.selection.getSelected()[0];
	//forceRefreshTableVariation(varurl,gridVar,filterPatientName);
	var pid =  varsel["patient_id"].join(",");
	var pname =  varsel["allpatients"].join(",");
	var allele = varsel["text"][0];
	

	if (bigbutton) {
		alert("disable for now");
		return;
		sendValidations(varsel["id"],varsel["start"],varsel["chromosome"], pid, value,username,allele,pname, gridbig, gridsmall);
		document.getElementById("csmall").style.visibility = "hidden";
	}
	else{
		var l = gridsmall.selection.getSelected().length;
		var patArray = new Array;
		var nameArray = new Array;
		for (var i = 0; i < l; i++) {
			patArray.push(gridsmall.selection.getSelected()[i]["patient_id"]);
			nameArray.push(gridsmall.selection.getSelected()[i]["name"]);
			
		}

		var st = varsel["text"];
		
		
		st = st+"";
		var tt = st.split("/");
		if (varsel["type"] == "insertion"){
			tt[0] = "-";
		}
		st = tt[0]+"_"+tt[1];
		
		var poly_id = varsel["chromosome"]+"_"+varsel["start"]+"_"+st;
		var pid2 =  patArray.join(",");
		var pname2 = nameArray.join(",");
		
		sendValidations(varsel["id"],varsel["start"],varsel["chromosome"], pid2 , value, username,allele,pname2,gridbig,gridsmall,poly_id);
		document.getElementById("csmall").style.visibility = "hidden";
		
	}
	viewAnnex(gridbig.selection.getSelected());
	return;
		
	
}



function sendValidations(var_id,position,chromosome,patient_ids,value,username,allele,pname,gridbig,gridsmall,poly_id){
	
 var url = url_valid+"?variation_id="+var_id+"&samples="+pname+"&project="+projectName+"&user="+username+"&value="+value+"&poly_id="+poly_id;
   dojo.xhrGet({ 
        // The following URL must match that used to test the server.
        url: url,
        handleAs: "json",     
        timeout: 80000, // Time in milliseconds
        // The LOAD function will be called on a successful response.
        
        load: function(responseObject, ioArgs,evt){
        	if (responseObject["status"] == "OK"){
				forceRefreshTableVariation(gridbig,filterPatientName);
				gridbig.store.fetch({ 
				query: { name:"*"},
				onComplete: function(items, result){			
						 gridbig.setQuery(filter);
                }
				});
										 
				 return;
			}		
			alert(responseObject["message"]);
        
            
        },
        
        // The ERROR function will be called in an error case.
        error: function(response, ioArgs){ // &acirc;&#382;&fnof;
        	alert(responseObject["status"]);        
            return;
        }
        
        
    });
 

}



//fonction pour envoyer les donnees d'un formulaire via une url
function sendData(url,fo){
	alert(url);
	dojo.xhrGet({
        url: url,
        method: "post",
        handleAs: "text",
        form: fo ,
        handle: function(data, ioArgs){
            var foo = dojo.fromJson(data);
            if (foo.status == "OK") {		
			
				return 1;
            }
            else {
			
				textError.setContent(foo.message);
				myError.show();
				 return -1;
            }          
        }     
    }); 
}
		


