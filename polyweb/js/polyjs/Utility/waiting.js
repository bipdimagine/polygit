/**
 * @author pnitschk
 */
	
var old_cache =0;
var nb_cache = 0;
var sleep = 2000;
var first = 0;
function isCache(url) {
    
	run_after_waiting();
	return;
	 dojo.xhrGet({ // ➂
        // The following URL must match that used to test the server.
        url: url,
        handleAs: "json",
        timeout: 100000, // Time in milliseconds
        // The LOAD function will be called on a successful response.
        load: function(responseObject, ioArgs){
            // Now you can just use the object
			//pour les projets de séquencage classique on trace des chromatogrammes
				
		
			if (responseObject["message"] == "OK") {
					//dijit.byId('waiting_cache').hide();
					
					run_after_waiting();
					
			}
		
	
            
        },
        
        // The ERROR function will be called in an error case.
        error: function(response, ioArgs){ // ➃
    
			alert("error in creating cache ");
			console.dir(ioArgs);
			console.error(url_chart); // ➆
            console.error("HTTP status code:  _--> ", ioArgs.xhr.status); // ➆
            return response; // ➅
        }
    });
	
}