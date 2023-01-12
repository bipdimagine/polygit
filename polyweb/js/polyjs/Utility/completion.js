
function genecompletion(value){
      var leStore = new dojo.store.Cache(new dojo.store.JsonRest({ target :
            url_find_name+"?project="+projectName, idProperty: "id"
            }), dojo.store.Memory());
 
        var testStore = new dojo.store.Memory({data: []});
        
         if(dojo.byId("listcompletion").value.length > 2 )//&& dijit.byId("listcompletion").get("store") == testStore)
                        dijit.byId("listcompletion").set("store", leStore);
                    if(dojo.byId("listcompletion").value.length <= 2  && dijit.byId("listcompletion").get("store") == leStore)
                        dijit.byId("listcompletion").set("store", testStore);
    
}