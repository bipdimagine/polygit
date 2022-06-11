
function formatCoverage(value, rowId, cellObj){
    var obj = dojo.fromJson(value);
    var xr = 255;
    var xg = 10;
    var xb = 0;

    var yr = 0;
    var yg = 255;   
    var yb = 0;
    var max = 50;
    var mean = obj.mean;
    var min = obj.min;
    var limit =  limit_spinner.value;
    var flag =  obj.flag*1.0;
    var red;
    var green;
    var blue =50;
    if (flag == -1){
        red = 255;
        green =20;
    }
    else {
        value = value *1.0;
        var red = Math.round((xr + ((min * (yr-xr)) / max) ))-1;
        var green = Math.round((xg + ((min * (yg-xg)) / max) ))-1;
        var blue = 50;
        if (red <0) red =0;
    }
        cellObj.customStyles.push("background:rgb("+red+","+green+","+blue+");color:black;"); 
    
    return  obj.mean+" ("+obj.min+")";
}
 
var table_transcript = new Array();
var title_transcript = new Array();
function coverage_json(arg,transcripts) {
        var tsamples = dijit.byId("gridPatients").selection.getSelected();
        var p="";
   
        var layout = new Array();
        layout[0] = { name: "exon",field: "id",  width: "5em", styles: 'text-align:center;'};
        for (i=0;i<tsamples.length;i++){
                var name =  tsamples[i].label;
                layout[i+1] = { name: name+'',field: name+'', formatter:formatCoverage, width: "auto", styles: 'font-size:09px;text-align:center;'};
        }
        if (!transcripts){
           
         transcripts = dijit.byId("gridGenes").selection.getSelected();
        }
        
         //console.dir(dijit.byId("gridGenes").getRows);//.toggleRow(0, true);
        var t="";
        
        var limit =  limit_spinner.value;
        //transcripts.length
        for (i=0;i<transcripts.length;i++){
            var tname =  transcripts[i].transcript;
            var nb_ex =  transcripts[i].nb_exons;
            var gene = transcripts[i].label;
             var desc = transcripts[i].description;
            var h = 50 + nb_ex*25;
            layout[0] = { name: "exon",field: "id",  width: "5em", styles: 'text-align:center;'};
            
            var url = url_coverage_json+"?project="+project_name+"&patient=all&transcript="+tname+"&limit="+limit+""+arg;
            var store = new dojo.data.ItemFileReadStore({
            url: url,
            autoheight : true,
            autowidth : true,
            });
            store.fetch({
                onComplete: function(a){
                    var tname = a[0].transcript;
                    var exon_one =  a[0].id;
                  
                    if (a[0].flag == 0){
                          title_transcript[tname].titleNode.style.color = "green"; 
                               title_transcript[tname].set("open",0);
                    }
                    else {
                            title_transcript[tname].titleNode.style.color = "red"; 
                    }
                    if (exon_one == "none"){
                         title_transcript[tname].set("open",0);
                          title_transcript[tname].titleNode.style.color = "green"; 
                    }
                      var h = 50 + a.length*25;
                      if (a.length==1){
                          h = 70 + a.length*25;
                      }
                     title_transcript[tname].style.height = h+"px";
                    
                      dojo.byId(tname+"_div").style.height= h+"px";
                       title_transcript[tname].resize();
                }
            
            });
            if (dojo.byId(transcripts[i].transcript+'_grid')){
                 table_transcript[tname].setStructure(layout);
                table_transcript[tname].setStore(store);
                 title_transcript[tname].set("open",true);
            }else {
                var z = "<div id='"+tname+"_div' style='width: 100%;height:"+h+"px;' ></div>";
                 var tid = "title_"+tname;
                title_transcript[tname] = new dijit.TitlePane({id:tid,title:gene+"-"+tname+"-"+desc,content:z,style:'width: 100%;height:"+h+"px;padding-bottom:5px'});
                dojo.byId("table_coverage2").appendChild(title_transcript[tname].domNode);
                title_transcript[tname].startup();
               
                
                table_transcript[tname] = new dojox.grid.EnhancedGrid({
                        name:transcripts[i].transcript,
                         id: transcripts[i].transcript+'_grid',
                         store: store,
                         structure: layout,
                         plugins: {menus:{headerMenu:"headerMenu", cellMenu:"cellMenu"}},
                         rowSelector: '20px'},
                     document.createElement('div'));

           dojo.connect(table_transcript[tname],'onCellDblClick',test_click);
            dojo.connect(table_transcript[tname],'onRowContextMenu',test_click);
          dojo.byId(tname+"_div").appendChild(table_transcript[tname].domNode);

        table_transcript[tname].startup();
            }
        
        } 
    
}
var select_exon;
var select_patient;
var select_transcript;

function click_exon_patient(e){
     load_graph_exon (select_patient,select_transcript,select_exon);
}
function click_exon(e){
     load_graph_one_exon_all_patients(select_transcript,select_exon);
}
function click_patient(e){
     load_graph_transcript (select_patient,select_transcript);
}
function click_transcript(e){
     load_graph_gene (select_transcript);
}

function test_click(e){
    select_patient = e.cell.field; // field name
    var rowIndex = e.rowIndex; // row index
    select_transcript = e.grid.name;
    select_exon = e.grid.getCell(0).getNode(rowIndex).textContent;
}
