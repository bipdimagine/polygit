
      
                var layoutGenes = [
     { field: 'label', name: 'name',width: 8, styles: 'text-align:center;'}, 
      { field: 'transcript', name: 'transcript',width: 14, styles: 'text-align:center;'}, 
       { field: 'chr', name: 'chr',width: 3, styles: 'text-align:center;'}, 
        { field: 'start', name: 'start',width: 9, styles: 'text-align:center;'}, 
         { field: 'end', name: 'end',width: 9, styles: 'text-align:center;'}, 
          { field: 'refseq', name: 'refseq',width: 12, styles: 'text-align:center;'}, 
           { field: 'ccds', name: 'ccds',width: 7, styles: 'text-align:center;'}, 
        { field: 'nb_exons', name: 'ex',width: 3, styles: 'text-align:center;'}, 
       { field: 'bundle', name: 'group',width: 30, styles: 'text-align:left;'}, 
      { field: 'description', name: 'description',width: 25, styles: 'text-align:left;'}, 
    
    ];
 
        var layoutBundles = [
     { field: 'label', name: 'name',width: 8, styles: 'text-align:center;'}, 
      { field: 'transcript', name: 'nb',width: 5, styles: 'text-align:center;'}, 
      { field: 'description', name: 'description',width: 50, styles: 'text-align:left;'}, 
    
    ];
 
   var layoutPatients = [
        { field: 'label', name: 'name',width: 20, styles: 'text-align:center;font-size:8px;'}, 
        { field: 'fam', name: 'familly',width: 15, styles: 'text-align:center;'},
		 { field: 'icon', name: 'status',width: 4, styles: 'text-align:center;',formatter:icon_sex},
        { field: 'validation_term', name: 'validation',width: 14, styles: 'text-align:center;'},
		{ field: 'validation_user', name: 'user',width: 8, styles: 'text-align:center;'},
		{ field: 'validation_date', name: 'date',width: 8, styles: 'text-align:center;'},
    
    ];  
 var nb_patients;
 function icon_sex(value) {
	return "<span class='dijitReset dijitInline dijitIcon  " +value+"'></span>";
	}
 function fillTabGroups(items , request ) {
 
        dijit.byId("gridBundles").setStore(storeGroup);
       
    //  dijit.byId("gridBundles").rowSelectCell.toggleAllSelection(true); 
    }   
    
    function fillTabPatients(items , request ) {
           dijit.byId("gridPatients").setStore(storeP);
            change_loading_button("pwaiting","#3498DB");
           nb_patients = items.length;
           dijit.byId("gridPatients").rowSelectCell.toggleAllSelection(true);
           document.getElementById("label_patients").innerHTML = 'all';
           change_label_patients();
           preload_coverage(storeP,0);
    }        
    var nb_genes;
    function fillTabGenes(items , request ) {
    	
         if (dijit.byId('waiting')) {
            dijit.byId('waiting').hide();
        }
        var tz ="";
         change_loading_button("pwaiting","#3498DB","panel-success");
        change_label_genes(items);
        nb_genes = items.length;

         dijit.byId("gridGenes").setStore(storeG);
         dijit.byId("gridGenes").rowSelectCell.toggleAllSelection(true);
        
          
               dijit.byId("gridGenes").render();
               dijit.byId("gridGenes").startup();
           var tsamples = dijit.byId("gridGenes").selection.getSelected();
    
         

    }
    
    
    
    function change_label_genes(transcripts){
    
        var tt = new Array();
        var tz ="";

        
        if (!transcripts){
            
            transcripts = dijit.byId("gridGenes").selection.getSelected();
            if (transcripts.length == nb_genes){
                transcripts = AllTranscripts ;
            }
        }
        else {
            AllTranscripts = transcripts;
        }
    
        var tz="";
    
        for (i=0;i<transcripts.length;i++){
            
            //console.log(dojo.getById(transcripts[i].label));
            //dojo.getById(transcripts[i].label).removeClass();// = "bulletGreen";
            tz += ""+transcripts[i].label+"&nbsp;";
        } 
    
    
    
        document.getElementById("label_genes").innerHTML = tz;
}


function return_options() {
    var url_option="";
  
    
    if (typeof limit_spinner == 'undefined') {
    	url_option += "&span=20";
        url_option += "&limit=30";
        url_option += "&intronic=0";
        url_option += "&utr=0";
        
        return url_option;
    }
        limit =  limit_spinner.value;
        padding = padding_spinner.value;
        if (utr){
            if (utr.getValue() == 1 ){
                url_option += "&utr=1";
            }
        }
        if (intronic){
            if (intronic.getValue() == 1 ){
                url_option += "&intronic=1";
            }
        }
      
        url_option += "&span="+padding;
        url_option += "&limit="+limit;
        return url_option;
}
 
function construct_url(url,tr) {
		if (typeof(panel) !== 'undefined') {
	        url = url+"?project="+project_name+"&panel="+panel;
	    }
	    else {
	    	url = url+"?project="+project_name;
	    }    
        url += return_options();
         
        /* sample */
        var tsamples = dijit.byId("gridPatients").selection.getSelected();
        var p="";
        if (tsamples.length == nb_patients){
            p="all";
        }   
    
        else if (tsamples.length == 0) {
            p="all";
        }
        else {
            
            for (i=0;i<tsamples.length;i++){
                p += tsamples[i].label+",";
            }
        }
        url += "&patients="+p;
        if (tr){
        	 url += "&transcripts="+tr;
        	 return url;
        }
        
         var transcripts = dijit.byId("gridGenes").selection.getSelected();
         var t="";
         if(transcripts[0] == null ){
        	 
        	 t="all";
         }
         
         
         else if (transcripts.length != nb_genes && transcripts.length > 0){
        	
        	 for (i=0;i<transcripts.length;i++){
        		 t += transcripts[i].transcript+",";
        
        	 } 
        }
        if (t== ""){
            t="all";
        }
        
        document.getElementById("label_patients").innerHTML = p;
        url += "&transcripts="+t;
       
         var username = dojo.cookie("username");
        url += "&user_name="+username;
        return url;
        
     }
     
function change_label_patients (){
    var tsamples = dijit.byId("gridPatients").selection.getSelected();
        var p="";
            
            for (i=0;i<tsamples.length;i++){
               if (!tsamples[i]) continue ;
            
                p += tsamples[i].label+",";
            }
        
        document.getElementById("label_patients").innerHTML = p;
        
    
}